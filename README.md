# brep-conformance-arris

The [brep-conformance](https://github.com/springcraft-ai/brep-conformance) adapter for
[arris](https://github.com/arris-labs/arris): it puts arris behind the suite's protocol, so that the
suite can score it.

## Who maintains this

brep-conformance scores a B-rep kernel through an adapter: a small program that answers the suite's
requests with the kernel's own calls. Springcraft, which maintains the suite, wrote this adapter and
maintains it.

- **It is not arris's.** Arris's maintainers did not write it or choose the tests it runs.
- **A failure can be the adapter's.** A case arris fails in the results may be this adapter driving
  it wrongly. Corrections are welcome: open an issue or a pull request.
- **It adds nothing to the kernel.** The adapter is pass-through, as the suite's
  [contract](https://github.com/springcraft-ai/brep-conformance/blob/main/doc/contract.md#what-an-adapter-may-do)
  asks: one call of the kernel for each function, and anything the kernel has no one call for is
  listed as unsupported. What is scored is arris, not this code.
- **Arris's maintainers are welcome here.** They know best how to drive arris, so their review and
  corrections carry the most weight. Ask in an issue to be added as a reviewer.
- **The suite decides what is scored.** An official run uses only a commit of this adapter that the
  suite's maintainers have reviewed as pass-through, whoever wrote it.

## Overview

[`src/main.rs`](src/main.rs) answers each function with arris's one call for it, as
[the contract](https://github.com/springcraft-ai/brep-conformance/blob/main/doc/contract.md#what-an-adapter-may-do)
asks, and lists the rest as unsupported. It builds the suite's box and cylinder with arris's
primitives, each at the origin and moved onto its frame; extrusions, revolutions, lofts and sweeps
with arris's own, sweeps along its paths of lines and arcs. It answers the booleans, the booleans of
many, mass properties, the full validity check, classification, the nearest point, distance, ray
fire, meshes, STEP, fillets and chamfers, face pushes, shells, splits, plane sections, section areas
and the interference of a pair with arris's own, naming edges and faces by arris's distance from the
point. What the checker left unchecked is reported with its violations.

A screw is arris's sweep along a helix with its binormal held along the axis, and a sweep along a
spline is its sweep along a NURBS path, in its rotation-minimizing frame.

A display is one call of arris's mesher with its render block: each face's own vertices with the
surface's outward normals, and the polyline of each edge it meshes between two faces. A seam, which
arris's topology lists as used twice by one face, and an edge between two faces on one surface are
left out, as the protocol's display asks; the choice is arris's own edge uses, with no number of the
adapter's.

Listed as unsupported, since arris 0.10 has no one call of the protocol's shape for it: the sphere,
torus and cone constructors (the suite sends each as a revolution where a case only starts from
one), slabs, drafts, the offset of a whole solid, stretches and shears, blends of unequal or varying
size in one operation, a twisted sweep, face deletes, moves and thickens, sewing, curvatures at a
point, the least bounding box (arris gives boxes round control points and mesh vertices), and Betti
numbers (arris gives lumps, voids and a genus).

arris moves a body rigidly in one call, so a scale and a mirror are refused. It reads each solid of
a STEP file as a body of its own, so a file of several is refused.

arris has no empty body. A boolean whose result has no material, or no thickness, is refused with a
reason that says so, and the adapter reads those two reasons as the empty solid: a regularized
boolean of solids that only touch is empty. Every other refusal is passed on, among them one body as
both operands of a boolean.

## Setup

The Rust toolchain [`rust-toolchain.toml`](rust-toolchain.toml) pins.

## Usage

The suite runs the adapter as a separate process. Build the adapter here, then run the suite from a
checkout of [brep-conformance](https://github.com/springcraft-ai/brep-conformance) beside this one:

```bash
cargo build --release
cd ../brep-conformance
cargo build --release -p brep-conformance
target/release/brep-conformance --adapter ../brep-conformance-arris/target/release/arris-adapter \
    --out arris.json
```

## Licences

The adapter in this repository is under the MIT licence: see [LICENSE](LICENSE).

That licence covers the adapter and nothing else. It does not cover the kernel, and nothing here
relicenses the kernel.

arris is a separate work under its own licence. Its crate declares `MIT OR Apache-2.0`; its
repository is the authority on that, not this page.

Building the adapter downloads arris from crates.io and links it in, so a built adapter contains
arris and is subject to arris's terms as well as this repository's.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).
