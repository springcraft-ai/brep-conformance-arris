//! The arris adapter: answers the B-rep protocol's requests on stdin and stdout, so the
//! conformance suite can drive this kernel.

use std::f64::consts::TAU;
use std::process::ExitCode;
use std::sync::Mutex;

use arris::Control;
use arris::check::{Classification, Classifier, Level, check};
use arris::geom::{
    FrameLaw, NurbsCurve, Path, PathSegment, Profile, ProfileLoop, ProfileSegment,
    SurfaceIntersection, intersect_surfaces,
};
use arris::math::nalgebra::{Matrix3, Rotation3, UnitQuaternion};
use arris::math::{Axis, Frame, Isometry, Meter, Point2, Point3, Precision, UnitVec3, Vec3};
use arris::ops::measure::{mass_properties, surface_properties};
use arris::ops::query::{
    FootParam, Interference, Operand, Ray, RayHits, distance, frame_at, interference, ray_fire,
};
use arris::ops::{BooleanReason, InputReason, LoftKind, OpError, Reason, ShellSide};
use arris::ops::{
    chamfer, common, cut, cut_many, extrude, fillet, fuse, fuse_many, loft, offset_faces,
    primitive_box, primitive_cylinder, revolve, section, shell, sweep, transform,
};
use arris::topo::provenance::Provenance;
use arris::topo::{Body, Edge, Face, Model, Shape};
use brep_protocol::{Meeting, Refusal, model, model::PathSeg, model::Primitive, profile};

struct Arris(Mutex<Model>);

const NONE: &Control<'_> = &Control::NONE;

/// arris's error as the protocol's refusal. An operation's error carries its kind: an input or a
/// result that is no solid is the request's fault, a pair of surfaces arris has no formula for is
/// unsupported, and anything else is a failure.
fn err<E: std::fmt::Display + 'static>(e: E) -> Refusal {
    let reason = e.to_string();
    match (&e as &dyn std::any::Any).downcast_ref::<OpError>() {
        Some(OpError::Unsupported { .. }) => Refusal::unsupported(reason),
        Some(
            OpError::InvalidInput { .. }
            | OpError::Degenerate { .. }
            | OpError::Profile(_)
            | OpError::Path(_),
        ) => Refusal::invalid(reason),
        _ => Refusal::failed(reason),
    }
}

fn p3(a: model::Vec3) -> Point3 {
    Point3::new(a[0], a[1], a[2])
}

fn v3(a: model::Vec3) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

/// The rotation whose columns are `x`, `y` and `z`.
fn rotation(x: Vec3, y: Vec3, z: Vec3) -> UnitQuaternion<f64> {
    let r = Matrix3::from_columns(&[x, y, z]);
    UnitQuaternion::from_rotation_matrix(&Rotation3::from_matrix_unchecked(r))
}

/// The motion taking a solid built about the world origin onto `frame`.
fn placement(f: &model::Frame) -> Result<Isometry, Refusal> {
    let (x, y, z) = (v3(f.x), v3(f.y), v3(f.z));
    if x.cross(&y).dot(&z) < 0.0 {
        return Err("a left-handed frame".into());
    }
    Ok(Isometry::new(rotation(x, y, z), v3(f.origin)))
}

/// The suite's profile as one of arris's, drawn in `plane`: its first loop the outer, the rest
/// holes. arris draws an arc through a point on it, so each arc goes as its end through its middle,
/// and a whole circle as two halves.
fn sketch(p: &profile::Profile, plane: Frame) -> Result<Profile, Refusal> {
    let at = |c: [f64; 2], r: f64, a: f64| Point2::new(c[0] + r * a.cos(), c[1] + r * a.sin());
    let path = |lp: &[profile::Seg]| -> Result<ProfileLoop, Refusal> {
        let mut segments = Vec::new();
        let mut start = None;
        for s in lp {
            match *s {
                profile::Seg::Line { a, b } => {
                    start.get_or_insert(Point2::new(a[0], a[1]));
                    segments.push(ProfileSegment::LineTo(Point2::new(b[0], b[1])));
                }
                profile::Seg::Arc {
                    center,
                    radius,
                    start: a0,
                    sweep,
                } => {
                    start.get_or_insert(at(center, radius, a0));
                    let halves = if sweep.abs() >= TAU - 1e-12 { 2 } else { 1 };
                    let step = sweep / f64::from(halves);
                    for k in 0..halves {
                        let a = a0 + step * f64::from(k);
                        segments.push(ProfileSegment::ArcTo {
                            to: at(center, radius, a + step),
                            via: at(center, radius, a + step / 2.0),
                        });
                    }
                }
                profile::Seg::Spline(_) => {
                    return Err(Refusal::unsupported(
                        "arris profiles have no spline segment",
                    ));
                }
                profile::Seg::Conic { .. } => {
                    return Err(Refusal::unsupported("arris profiles have no conic segment"));
                }
                profile::Seg::Gn { .. } => {
                    return Err(Refusal::unsupported(
                        "arris profiles have no blend-section segment",
                    ));
                }
                _ => {
                    return Err(Refusal::unsupported(
                        "the adapter does not know this kind of profile segment",
                    ));
                }
            }
        }
        Ok(ProfileLoop::Path {
            start: start.ok_or("an empty loop")?,
            segments,
        })
    };
    let mut loops = p.loops.iter();
    let outer = path(loops.next().ok_or("a profile with no loop")?)?;
    let holes = loops.map(|l| path(l)).collect::<Result<Vec<_>, _>>()?;
    Ok(Profile {
        plane,
        outer,
        holes,
    })
}

/// The frame's xy plane, as arris's profiles take it.
fn plane_xy(f: &model::Frame, lift: f64) -> Result<Frame, Refusal> {
    let o = v3(f.origin) + v3(f.z) * lift;
    Frame::new(Point3::from(o), v3(f.z), v3(f.x)).map_err(err)
}

/// The plane a revolution's profile is drawn in: x out from the frame's z axis, y along it.
fn plane_meridian(f: &model::Frame) -> Result<Frame, Refusal> {
    Frame::new(p3(f.origin), -v3(f.y), v3(f.x)).map_err(err)
}

/// The suite's sweep path as arris's, from the frame's origin starting along its z: each arc bends
/// toward the current x about an axis along the current y.
fn trace(f: &model::Frame, segs: &[PathSeg]) -> Path {
    let (mut at, mut x, y, mut z) = (p3(f.origin), v3(f.x), v3(f.y), v3(f.z));
    let mut segments = Vec::new();
    for s in segs {
        match *s {
            PathSeg::Line { length } => {
                at += z * length;
                segments.push(PathSegment::LineTo(at));
            }
            PathSeg::Arc { radius, angle } => {
                let c = at + x * radius;
                let turn = |a: f64| {
                    let (s, co) = a.sin_cos();
                    (
                        c - x * (radius * co) + z * (radius * s),
                        x * co - z * s,
                        z * co + x * s,
                    )
                };
                let (via, _, _) = turn(angle / 2.0);
                let (to, nx, nz) = turn(angle);
                segments.push(PathSegment::ArcTo { via, to });
                (at, x, z) = (to, nx, nz);
            }
        }
    }
    let _ = y;
    Path {
        start: p3(f.origin),
        segments,
    }
}

fn make(m: &mut Model, p: &Primitive) -> Result<Body, Refusal> {
    let o = Point3::origin();
    let (body, frame) = match p {
        Primitive::Box { frame, size } => {
            let max = Point3::new(size[0], size[1], size[2]);
            (primitive_box(m, o, max, NONE).map_err(err)?.0, frame)
        }
        Primitive::Cylinder {
            frame,
            radius,
            height,
        } => (
            primitive_cylinder(m, Axis::z_at(o), *radius, *height, NONE)
                .map_err(err)?
                .0,
            frame,
        ),
        // The profile constructions are built in place, already on their frames.
        Primitive::Extrude {
            frame,
            profile,
            height,
        } => {
            let prof = sketch(profile, plane_xy(frame, 0.0)?)?;
            return Ok(extrude(m, &prof, v3(frame.z), *height, NONE)
                .map_err(err)?
                .0);
        }
        Primitive::Revolve {
            frame,
            profile,
            angle,
        } => {
            let prof = sketch(profile, plane_meridian(frame)?)?;
            let axis = Axis::new(p3(frame.origin), v3(frame.z)).map_err(err)?;
            return Ok(revolve(m, &prof, axis, *angle, NONE).map_err(err)?.0);
        }
        Primitive::Loft {
            frame,
            profile,
            height,
            top,
            ..
        } => {
            let sections = [
                sketch(profile, plane_xy(frame, 0.0)?)?,
                sketch(top, plane_xy(frame, *height)?)?,
            ];
            return Ok(loft(m, &sections, LoftKind::Ruled, false, NONE)
                .map_err(err)?
                .0);
        }
        Primitive::Sweep {
            frame,
            profile,
            path,
        } => {
            let prof = sketch(profile, plane_xy(frame, 0.0)?)?;
            let path = trace(frame, path);
            return Ok(sweep(m, &prof, &path, FrameLaw::RotationMinimizing, NONE)
                .map_err(err)?
                .0);
        }
        // arris's screw: its sweep along a helix, the profile's binormal held along the axis.
        Primitive::Screw {
            frame,
            profile,
            rise,
            angle,
        } => {
            let prof = sketch(profile, plane_meridian(frame)?)?;
            let axis = Axis::new(p3(frame.origin), v3(frame.z)).map_err(err)?;
            // The helix runs from a point of the profile's plane off the axis.
            let reach = profile.bounds().1[0];
            let path = Path {
                start: p3(frame.origin) + v3(frame.x) * reach,
                segments: vec![PathSegment::Helix {
                    axis,
                    pitch: rise * TAU,
                    angle: *angle,
                }],
            };
            let law = FrameLaw::FixedBinormal(UnitVec3::new_normalize(v3(frame.z)));
            return Ok(sweep(m, &prof, &path, law, NONE).map_err(err)?.0);
        }
        Primitive::SweepPath {
            frame,
            profile,
            path,
            twist,
        } => {
            if *twist != 0.0 {
                return Err(Refusal::unsupported(
                    "a twisted sweep: arris's sweep has no twist law",
                ));
            }
            let prof = sketch(profile, plane_xy(frame, 0.0)?)?;
            let (o, x, y, z) = (p3(frame.origin), v3(frame.x), v3(frame.y), v3(frame.z));
            let ctrl: Vec<Point3> = path
                .ctrl
                .iter()
                .map(|c| o + x * c[0] + y * c[1] + z * c[2])
                .collect();
            let weights = vec![1.0; ctrl.len()];
            let start = *ctrl.first().ok_or("a sweep with no path")?;
            let curve =
                NurbsCurve::new(path.degree, path.knots.clone(), ctrl, weights).map_err(err)?;
            let path = Path {
                start,
                segments: vec![PathSegment::Nurbs(curve)],
            };
            return Ok(sweep(m, &prof, &path, FrameLaw::RotationMinimizing, NONE)
                .map_err(err)?
                .0);
        }
        _ => return Err(Refusal::unsupported("this constructor is not supported")),
    };
    Ok(transform(m, body, &placement(frame)?, NONE).map_err(err)?.0)
}

/// A solid: a body, or `None` for the empty set, which arris reports as a refusal.
type Solid = Option<Body>;

#[derive(Clone, Copy)]
enum Bool {
    Fuse,
    Common,
    Cut,
}

/// The boolean of `a` and `b`. arris has no empty body: it answers a result with no material or
/// no thickness by a refusal that says so, which is the empty set.
fn boolean(m: &mut Model, op: Bool, a: Solid, b: Solid) -> Result<Solid, Refusal> {
    let (a, b) = match (op, a, b) {
        (_, Some(a), Some(b)) => (a, b),
        (Bool::Fuse, x, None) | (Bool::Fuse, None, x) => return Ok(x),
        (Bool::Cut, x, None) => return Ok(x),
        _ => return Ok(None),
    };
    let r = match op {
        Bool::Fuse => fuse(m, a, b, NONE),
        Bool::Common => common(m, a, b, NONE),
        Bool::Cut => cut(m, a, b, NONE),
    };
    empty_or(r)
}

/// A boolean's body, or the empty set where arris says the result has no material or no
/// thickness: a regularized boolean of solids that only touch is empty. Any other refusal is
/// passed on.
fn empty_or(r: Result<(Body, Provenance), OpError>) -> Result<Solid, Refusal> {
    match r {
        Ok((body, _)) => Ok(Some(body)),
        Err(OpError::Degenerate {
            reason:
                Reason::Boolean(BooleanReason::Empty) | Reason::Input(InputReason::ZeroThickness),
            ..
        }) => Ok(None),
        Err(e) => Err(err(e)),
    }
}

/// How far `p` is from an edge or face, by arris's own distance query, trims included.
fn gap(m: &Model, p: Point3, s: Shape) -> Option<f64> {
    let d = distance(m, Operand::Point(p), Operand::Shape(s), NONE).ok()?;
    Some(d.distance)
}

/// The edge of `body` each point lies on: the nearest, which must pass within 10·ε.
fn edges_at(m: &Model, body: Body, points: &[model::Vec3]) -> Result<Vec<Edge>, Refusal> {
    let all = m.edges(body).map_err(err)?;
    points
        .iter()
        .map(|&p| {
            all.iter()
                .filter_map(|&e| gap(m, p3(p), Shape::from(e)).map(|d| (d, e)))
                .filter(|(d, _)| *d <= 1e-6)
                .min_by(|a, b| a.0.total_cmp(&b.0))
                .map(|(_, e)| e)
                .ok_or_else(|| Refusal::invalid(format!("no edge passes through {p:?}")))
        })
        .collect()
}

/// The face of `body` each point lies on.
fn faces_at(m: &Model, body: Body, points: &[model::Vec3]) -> Result<Vec<Face>, Refusal> {
    let all = m.faces(body).map_err(err)?;
    points
        .iter()
        .map(|&p| {
            all.iter()
                .filter_map(|&f| gap(m, p3(p), Shape::from(f)).map(|d| (d, f)))
                .filter(|(d, _)| *d <= 1e-6)
                .min_by(|a, b| a.0.total_cmp(&b.0))
                .map(|(_, f)| f)
                .ok_or_else(|| Refusal::invalid(format!("no face passes through {p:?}")))
        })
        .collect()
}

/// `body` moved by the rows `[R | t]`. arris moves a body rigidly in one call: a scale or a mirror
/// is refused.
fn moved(m: &mut Model, body: Body, mat: [[f64; 4]; 3]) -> Result<Body, Refusal> {
    let r = Matrix3::from_fn(|i, j| mat[i][j]);
    if (r.determinant() - 1.0).abs() > 1e-9
        || (r * r.transpose() - Matrix3::identity()).norm() > 1e-9
    {
        return Err(Refusal::unsupported(
            "arris moves bodies rigidly: a scale or a mirror is not supported",
        ));
    }
    let t = Vec3::new(mat[0][3], mat[1][3], mat[2][3]);
    let q = UnitQuaternion::from_rotation_matrix(&Rotation3::from_matrix_unchecked(r));
    Ok(transform(m, body, &Isometry::new(q, t), NONE)
        .map_err(err)?
        .0)
}

/// The plane through `origin` square to `normal`, as a frame with some x.
fn plane_frame(origin: model::Vec3, normal: model::Vec3) -> Result<Frame, Refusal> {
    let z = v3(normal);
    let x = if z.x.abs() < 0.9 * z.norm() {
        Vec3::x()
    } else {
        Vec3::y()
    };
    let x = x - z * (x.dot(&z) / z.dot(&z));
    Frame::new(p3(origin), z, x).map_err(err)
}

impl brep_protocol::Kernel for Arris {
    type Solid = Solid;

    fn name(&self) -> String {
        "arris".into()
    }

    fn version(&self) -> String {
        "0.10.0".into()
    }

    fn tolerance(&self) -> f64 {
        self.0.lock().unwrap().precision().default_tolerance
    }

    fn make(&mut self, p: &Primitive) -> Result<Solid, Refusal> {
        make(&mut self.0.lock().unwrap(), p).map(Some)
    }

    fn union(&mut self, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
        boolean(&mut self.0.lock().unwrap(), Bool::Fuse, *a, *b)
    }

    fn intersection(&mut self, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
        boolean(&mut self.0.lock().unwrap(), Bool::Common, *a, *b)
    }

    fn difference(&mut self, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
        boolean(&mut self.0.lock().unwrap(), Bool::Cut, *a, *b)
    }

    fn union_all(&mut self, solids: &[&Solid]) -> Result<Solid, Refusal> {
        let bodies: Vec<Body> = solids.iter().filter_map(|s| **s).collect();
        let m = &mut *self.0.lock().unwrap();
        match bodies.len() {
            0 => Ok(None),
            1 => Ok(Some(bodies[0])),
            _ => empty_or(fuse_many(m, &bodies, NONE)),
        }
    }

    fn difference_all(&mut self, target: &Solid, tools: &[&Solid]) -> Result<Solid, Refusal> {
        let Some(target) = *target else {
            return Ok(None);
        };
        let tools: Vec<Body> = tools.iter().filter_map(|s| **s).collect();
        if tools.is_empty() {
            return Ok(Some(target));
        }
        let m = &mut *self.0.lock().unwrap();
        empty_or(cut_many(m, target, &tools, NONE))
    }

    fn volume(&mut self, s: &Solid) -> Result<f64, Refusal> {
        let Some(s) = *s else { return Ok(0.0) };
        let m = self.0.lock().unwrap();
        mass_properties(&m, s, NONE).map(|p| p.volume).map_err(err)
    }

    fn area(&mut self, s: &Solid) -> Result<f64, Refusal> {
        let Some(s) = *s else { return Ok(0.0) };
        let m = self.0.lock().unwrap();
        surface_properties(&m, s, NONE).map(|p| p.area).map_err(err)
    }

    fn check(&mut self, s: &Solid) -> Result<Vec<String>, Refusal> {
        let Some(s) = *s else {
            return Ok(Vec::new());
        };
        let m = self.0.lock().unwrap();
        let report = check(&m, s, Level::Full);
        // What the checker could not decide is not shown valid.
        let undecided = report
            .unchecked()
            .iter()
            .map(|u| format!("unchecked: {u:?}"));
        Ok(report
            .violations()
            .iter()
            .map(|v| format!("{v:?}"))
            .chain(undecided)
            .collect())
    }

    fn classify(&mut self, s: &Solid, points: &[model::Vec3]) -> Result<Vec<model::Side>, Refusal> {
        let Some(s) = *s else {
            return Ok(vec![model::Side::Out; points.len()]);
        };
        let m = self.0.lock().unwrap();
        let c = Classifier::of_body(&m, s).map_err(err)?;
        points
            .iter()
            .map(|&p| {
                Ok(match c.classify(p3(p)).map_err(err)? {
                    Classification::Inside => model::Side::In,
                    Classification::On(_) => model::Side::On,
                    Classification::Outside => model::Side::Out,
                })
            })
            .collect()
    }

    fn closest(&mut self, s: &Solid, points: &[model::Vec3]) -> Result<Vec<model::Vec3>, Refusal> {
        let Some(s) = *s else {
            return Err("an empty solid has no nearest point".into());
        };
        let m = self.0.lock().unwrap();
        points
            .iter()
            .map(|&p| {
                let body = Operand::Shape(Shape::from(s));
                let d = distance(&m, Operand::Point(p3(p)), body, NONE).map_err(err)?;
                Ok([d.b.point.x, d.b.point.y, d.b.point.z])
            })
            .collect()
    }

    fn centroid(&mut self, s: &Solid) -> Result<Option<model::Vec3>, Refusal> {
        let Some(s) = *s else { return Ok(None) };
        let m = self.0.lock().unwrap();
        let p = mass_properties(&m, s, NONE).map_err(err)?;
        Ok(Some([p.centroid.x, p.centroid.y, p.centroid.z]))
    }

    fn inertia(&mut self, s: &Solid) -> Result<[[f64; 3]; 3], Refusal> {
        let Some(s) = *s else {
            return Ok([[0.0; 3]; 3]);
        };
        let m = self.0.lock().unwrap();
        let p = mass_properties(&m, s, NONE).map_err(err)?;
        Ok([0, 1, 2].map(|i| [0, 1, 2].map(|j| p.inertia[(i, j)])))
    }

    fn distance(&mut self, a: &Solid, b: &Solid) -> Result<f64, Refusal> {
        let (Some(a), Some(b)) = (*a, *b) else {
            return Err("the distance to an empty solid is undefined".into());
        };
        let m = self.0.lock().unwrap();
        let d = distance(
            &m,
            Operand::Shape(Shape::from(a)),
            Operand::Shape(Shape::from(b)),
            NONE,
        )
        .map_err(err)?;
        Ok(d.distance)
    }

    fn raycast(&mut self, s: &Solid, o: model::Vec3, d: model::Vec3) -> Result<Vec<f64>, Refusal> {
        let Some(s) = *s else { return Ok(Vec::new()) };
        let m = self.0.lock().unwrap();
        let ray = Ray {
            origin: p3(o),
            direction: v3(d),
        };
        let hits = ray_fire(&m, s, &ray, RayHits::All).map_err(err)?;
        Ok(hits.iter().map(|h| h.t).collect())
    }

    fn transform(&mut self, s: &Solid, matrix: [[f64; 4]; 3]) -> Result<Solid, Refusal> {
        let Some(s) = *s else { return Ok(None) };
        moved(&mut self.0.lock().unwrap(), s, matrix).map(Some)
    }

    fn fillet(&mut self, s: &Solid, edges: &[model::Vec3], radius: f64) -> Result<Solid, Refusal> {
        let Some(s) = *s else {
            return Err("nothing to fillet".into());
        };
        let m = &mut *self.0.lock().unwrap();
        let found = edges_at(m, s, edges)?;
        Ok(Some(fillet(m, s, &found, radius, NONE).map_err(err)?.0))
    }

    fn chamfer(&mut self, s: &Solid, edges: &[model::Vec3], d: f64) -> Result<Solid, Refusal> {
        let Some(s) = *s else {
            return Err("nothing to chamfer".into());
        };
        let m = &mut *self.0.lock().unwrap();
        let found = edges_at(m, s, edges)?;
        Ok(Some(chamfer(m, s, &found, d, NONE).map_err(err)?.0))
    }

    fn offset_faces(&mut self, s: &Solid, faces: &[model::Vec3], d: f64) -> Result<Solid, Refusal> {
        let Some(s) = *s else {
            return Err("nothing to offset".into());
        };
        let m = &mut *self.0.lock().unwrap();
        let found = faces_at(m, s, faces)?;
        Ok(Some(offset_faces(m, s, &found, d, NONE).map_err(err)?.0))
    }

    fn shell(&mut self, s: &Solid, open: &[model::Vec3], t: f64) -> Result<Solid, Refusal> {
        let Some(s) = *s else {
            return Err("nothing to shell".into());
        };
        let m = &mut *self.0.lock().unwrap();
        let found = faces_at(m, s, open)?;
        Ok(Some(
            shell(m, s, &found, t, ShellSide::Inward, NONE)
                .map_err(err)?
                .0,
        ))
    }

    fn normal(
        &mut self,
        s: &Solid,
        points: &[model::Vec3],
    ) -> Result<Vec<Option<model::Vec3>>, Refusal> {
        let Some(s) = *s else {
            return Ok(vec![None; points.len()]);
        };
        let m = self.0.lock().unwrap();
        let faces = m.faces(s).map_err(err)?;
        // arris gives a normal on one face at its parameters. The point names the face and the
        // parameters by arris's own distance query; a point several faces pass through, on an
        // edge or a vertex, is one it has no call for.
        points
            .iter()
            .map(|&p| {
                let on: Vec<(Face, FootParam)> = faces
                    .iter()
                    .filter_map(|&f| {
                        let d = distance(&m, Operand::Point(p3(p)), Operand::Shape(f.into()), NONE)
                            .ok()?;
                        (d.distance <= 1e-6).then_some((f, d.b.param))
                    })
                    .collect();
                match on[..] {
                    [] => Ok(None),
                    [(f, FootParam::Face(uv))] => {
                        let frame = frame_at(&m, f, uv).map_err(err)?;
                        let n = frame.z().into_inner();
                        Ok(Some([n.x, n.y, n.z]))
                    }
                    _ => Err(Refusal::unsupported(
                        "a point on an edge or a vertex, where several faces meet",
                    )),
                }
            })
            .collect()
    }

    fn mesh(&mut self, s: &Solid, tolerance: f64) -> Result<brep_protocol::Mesh, Refusal> {
        let Some(s) = *s else {
            return Ok(brep_protocol::Mesh::default());
        };
        let m = self.0.lock().unwrap();
        let mesh = arris::mesh::tessellate(&m, s, tolerance, NONE).map_err(err)?;
        Ok(brep_protocol::Mesh::new(
            mesh.positions().to_vec(),
            mesh.triangles().to_vec(),
        ))
    }

    /// arris's mesher with its render block: each face's own vertices with the surface's outward
    /// normal, and the polyline of each edge with the faces that use it. An edge's polyline runs
    /// through vertices both of its faces have, so the turn along it is read from the normals the
    /// block gives each face there.
    fn display(&mut self, s: &Solid, tolerance: f64) -> Result<brep_protocol::Display, Refusal> {
        let Some(s) = *s else {
            return Ok(brep_protocol::Display::default());
        };
        let m = self.0.lock().unwrap();
        let request = arris::mesh::MeshRequest::new(tolerance).with_corners();
        let mesh = arris::mesh::tessellate_with(&m, s, &request, NONE).map_err(err)?;
        let corners = mesh
            .corners()
            .ok_or_else(|| Refusal::failed("the mesher gave no render block"))?;
        let faces = mesh
            .faces()
            .iter()
            .zip(corners.faces())
            .map(|(face, local)| {
                let first = local.vertices.start as u32;
                brep_protocol::DisplayFace::new(
                    local
                        .vertices
                        .clone()
                        .map(|v| mesh.positions()[corners.positions()[v] as usize])
                        .collect(),
                    corners.normals()[local.vertices.clone()].to_vec(),
                    corners.triangles()[face.triangles.clone()]
                        .iter()
                        .map(|t| t.map(|v| v - first))
                        .collect(),
                )
            })
            .collect();
        let index = |face| mesh.faces().iter().position(|f| f.face == face);
        // The normal face `k` has at the mesh vertex `v`, where it has one: at an apex each of
        // its triangles has a normal of its own.
        let normal = |k: usize, v: u32| -> Option<model::Vec3> {
            let mut at = corners.faces()[k]
                .vertices
                .clone()
                .filter(|&c| corners.positions()[c] == v)
                .map(|c| corners.normals()[c]);
            let first = at.next()?;
            at.all(|n| (0..3).all(|i| (n[i] - first[i]).abs() <= 1e-9))
                .then_some(first)
        };
        let cross = |a: model::Vec3, b: model::Vec3| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let dot = |a: model::Vec3, b: model::Vec3| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        // The angle between the normals of faces `a` and `b` over the polyline, least and
        // greatest. A face's triangles run counter-clockwise seen from outside, so the way one of
        // `a`'s runs along the polyline is the way `a` lies to its left, and a convex edge's
        // normals turn about that direction by the right-hand rule.
        let turn = |a: usize, b: usize, line: &[u32]| -> Option<[f64; 2]> {
            let runs = |p: u32, q: u32| {
                mesh.triangles()[mesh.faces()[a].triangles.clone()]
                    .iter()
                    .find_map(|t| {
                        (0..3).find_map(|i| {
                            let (x, y) = (t[i], t[(i + 1) % 3]);
                            ((x, y) == (p, q))
                                .then_some(1.0)
                                .or(((x, y) == (q, p)).then_some(-1.0))
                        })
                    })
            };
            let way = line.windows(2).find_map(|w| runs(w[0], w[1]))?;
            let mut range: Option<[f64; 2]> = None;
            for (i, &v) in line.iter().enumerate() {
                let (Some(n), Some(o)) = (normal(a, v), normal(b, v)) else {
                    continue;
                };
                let (p, q) = (line[i.saturating_sub(1)], line[(i + 1).min(line.len() - 1)]);
                let along = [0, 1, 2]
                    .map(|k| mesh.positions()[q as usize][k] - mesh.positions()[p as usize][k]);
                let c = cross(n, o);
                let angle = dot(c, c).sqrt().atan2(dot(n, o));
                let t = if way * dot(c, along) < 0.0 {
                    -angle
                } else {
                    angle
                };
                range = Some(range.map_or([t, t], |[lo, hi]: [f64; 2]| [lo.min(t), hi.max(t)]));
            }
            range
        };
        let mut edges = Vec::new();
        let mut seams = Vec::new();
        for e in mesh.edges() {
            let line = &mesh.edge_indices()[e.indices.clone()];
            let points: Vec<model::Vec3> =
                line.iter().map(|&v| mesh.positions()[v as usize]).collect();
            // The model lists an edge's uses in every body that shares it: these are this one's.
            let uses: Vec<&arris::topo::CoedgeRef> = m
                .edge_uses(e.edge)
                .unwrap_or(&[])
                .iter()
                .filter(|u| index(u.face).is_some())
                .collect();
            // An edge one face uses once is a pole or an apex of a solid, with no length.
            if uses.len() < 2 {
                continue;
            }
            // A seam is one face using the edge twice.
            if let [a, b] = uses[..]
                && a.face == b.face
            {
                seams.push(points);
                continue;
            }
            let named: Vec<usize> = uses.iter().filter_map(|u| index(u.face)).collect();
            let pair = match named[..] {
                [a, b] => Some((a, b)),
                _ => None,
            };
            let turned = pair.and_then(|(a, b)| turn(a, b, line));
            // Two faces whose normals agree all along the edge lie on one surface where arris's
            // intersection of their surfaces says the two coincide; a round meeting a flat face
            // smoothly only touches it.
            let smooth = turned.is_some_and(|t| t.iter().all(|a| a.abs() <= 1e-9));
            let split = smooth && {
                let surface = |u: &arris::topo::CoedgeRef| {
                    m.face(u.face)
                        .ok()
                        .and_then(|f| m.surface(f.surface()).ok())
                };
                match (surface(uses[0]), surface(uses[1]), mesh.aabb()) {
                    (Some(a), Some(b), Some(within)) => matches!(
                        intersect_surfaces(
                            a,
                            b,
                            &within,
                            Precision::DEFAULT.tolerance(),
                            &mut Meter::default()
                        ),
                        Ok(SurfaceIntersection::Coincident)
                    ),
                    _ => false,
                }
            };
            edges.push(brep_protocol::DisplayEdge::new(
                points,
                named.iter().map(|&k| k as u32).collect(),
                turned,
                split,
            ));
        }
        Ok(brep_protocol::Display::new(faces, edges, seams))
    }

    fn write_step(&mut self, s: &Solid) -> Result<String, Refusal> {
        let Some(s) = *s else {
            return Err("an empty solid has no STEP".into());
        };
        let m = self.0.lock().unwrap();
        arris::io::step::write(&m, &[s]).map_err(err)
    }

    fn read_step(&mut self, text: &str) -> Result<Solid, Refusal> {
        let m = &mut *self.0.lock().unwrap();
        let read = arris::io::step::read(m, text, &Default::default(), NONE).map_err(err)?;
        let bodies: Vec<Body> = read
            .solids
            .iter()
            .map(|s| {
                s.result
                    .as_ref()
                    .map(|r| r.body)
                    .map_err(|e| Refusal::failed(format!("{e:?}")))
            })
            .collect::<Result<_, _>>()?;
        match bodies[..] {
            [] => Err("the file holds no solid".into()),
            [body] => Ok(Some(body)),
            _ => Err(Refusal::unsupported(format!(
                "the file holds {} solids, and arris reads each as a body of its own",
                bodies.len()
            ))),
        }
    }

    fn loft_through(
        &mut self,
        frame: &model::Frame,
        profiles: &[profile::Profile],
        heights: &[f64],
    ) -> Result<Solid, Refusal> {
        if profiles.len() != heights.len() {
            return Err(Refusal::invalid(format!(
                "{} profiles at {} heights",
                profiles.len(),
                heights.len()
            )));
        }
        let m = &mut *self.0.lock().unwrap();
        let sections = profiles
            .iter()
            .zip(heights)
            .map(|(p, &h)| sketch(p, plane_xy(frame, h)?))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Some(
            loft(m, &sections, LoftKind::Smooth, false, NONE)
                .map_err(err)?
                .0,
        ))
    }

    fn interference(&mut self, a: &Solid, b: &Solid) -> Result<Meeting, Refusal> {
        let (Some(a), Some(b)) = (*a, *b) else {
            return Err("an empty solid meets nothing".into());
        };
        let m = &mut *self.0.lock().unwrap();
        let (volume, gap) = match interference(m, a, b, NONE).map_err(err)? {
            Interference::Overlapping { volume, .. } => (volume, 0.0),
            Interference::Touching { .. } => (0.0, 0.0),
            Interference::Apart(d) => (0.0, d.distance),
        };
        Ok(Meeting { volume, gap })
    }

    fn reset(&mut self) {
        *self.0.lock().unwrap() = Model::default();
    }

    fn split(
        &mut self,
        s: &Solid,
        origin: model::Vec3,
        normal: model::Vec3,
    ) -> Result<(Solid, Solid), Refusal> {
        let Some(s) = *s else { return Ok((None, None)) };
        let plane = plane_frame(origin, normal)?;
        let m = &mut *self.0.lock().unwrap();
        let cut = arris::ops::split(m, s, &plane, NONE).map_err(err)?;
        Ok((Some(cut.positive), Some(cut.negative)))
    }

    fn section(
        &mut self,
        s: &Solid,
        origin: model::Vec3,
        normal: model::Vec3,
        tolerance: f64,
    ) -> Result<Vec<Vec<model::Vec3>>, Refusal> {
        let Some(s) = *s else { return Ok(Vec::new()) };
        let plane = plane_frame(origin, normal)?;
        let m = &mut *self.0.lock().unwrap();
        let Some(cut) = section(m, s, &plane, NONE).map_err(err)? else {
            return Ok(Vec::new());
        };
        // The section's curves are the sheet's edges, and arris's mesh samples each edge once,
        // within the chord: those samples are the pieces.
        let mesh = arris::mesh::tessellate(m, cut.sheet, tolerance, NONE).map_err(err)?;
        Ok(mesh
            .edges()
            .iter()
            .filter_map(|e| mesh.edge_polyline(e.edge))
            .map(|line| line.iter().map(|&i| mesh.positions()[i as usize]).collect())
            .collect())
    }

    fn section_area(
        &mut self,
        s: &Solid,
        origin: model::Vec3,
        normal: model::Vec3,
    ) -> Result<f64, Refusal> {
        let Some(s) = *s else { return Ok(0.0) };
        let plane = plane_frame(origin, normal)?;
        let m = &mut *self.0.lock().unwrap();
        Ok(section(m, s, &plane, NONE)
            .map_err(err)?
            .map_or(0.0, |c| c.area))
    }

    /// What arris 0.10 has one call of the protocol's shape for. Two operations arris reaches in
    /// its own way are left out with those it does not have, since the adapter does not assemble
    /// them: the bounding box, where arris gives boxes round control points and mesh vertices and
    /// the protocol asks for the least box, and the Betti numbers, where arris gives lumps, voids
    /// and a genus and the protocol counts solids that touch as one piece.
    fn supports(&self) -> Vec<String> {
        [
            "box",
            "cylinder",
            "extrude",
            "revolve",
            "screw",
            "loft",
            "loft_through",
            "sweep",
            "sweep_path",
            "union",
            "intersection",
            "difference",
            "union_all",
            "difference_all",
            "fillet",
            "chamfer",
            "offset_faces",
            "shell",
            "write_step",
            "read_step",
            "volume",
            "area",
            "classify",
            "closest",
            "normal",
            "check",
            "centroid",
            "inertia",
            "distance",
            "interference",
            "section_area",
            "section",
            "split",
            "mesh",
            "display",
            "raycast",
            "transform",
        ]
        .map(String::from)
        .to_vec()
    }
}

fn main() -> ExitCode {
    match brep_protocol::serve(&mut Arris(Mutex::new(Model::default()))) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
