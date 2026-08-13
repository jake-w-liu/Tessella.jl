# Tessella.jl architecture and scope

Tessella is a Julia-native mesh generator for the geometry-to-mesh part of the
ASCENT electromagnetics workflow. It is informed by gmsh's architecture but is an
independent implementation. The package roadmap is complete; the live verification
record is [`STATUS.md`](STATUS.md).

## Scope boundary

Tessella owns:

- exact geometric decisions and simplex topology;
- 1-D, surface, and tetrahedral meshing;
- conforming boundary and multi-region recovery;
- uniform and callback-driven size control;
- mesh quality improvement and sliver reduction;
- native primitives, analytical surfaces/imprints, and mesh Boolean CSG;
- linear/P2 mesh validation and solver-consumable MSH/STL I/O.

Tessella does not attempt to reproduce gmsh as a whole. In particular, a pure-Julia
OpenCASCADE/NURBS kernel, GUI, post-processing, FEM solve, and the long tail of mesh
formats are outside package scope. ASCENT remains the solver. Geometry outside the
native analytical/CSG set enters as a triangulated boundary mesh.

## Implemented architecture

```text
Tessella
├── Predicates    adaptive exact orient/incircle/insphere, exact rationals, SoS
├── MeshTypes     compact simplex storage, topology, quality, CRC, validation
├── ExactMesh3D   Rational{BigInt} Delaunay kernel
├── IO            strict/atomic MSH v2.2/v4.1, STL, limited .geo metadata scan
├── Mesh2D        Delaunay, CDT, interior classification, quality refinement
├── SizeField     constant, callable, and pointwise-minimum size fields
├── Mesh1D        metric-length curve and segment discretization
├── MeshSurface   planar, cylindrical, and parametric surface meshing
├── Mesh3D        Delaunay, fills, partitions, sizing, flips, mesh Boolean CSG
├── RecoverCDT    exact conforming-Delaunay boundary/partition recovery
├── Optimize      quality reports, Laplacian/ODT/targeted sliver smoothing
├── Heal          surface defect and meshability diagnostics
├── Geometry      native closed primitive surfaces
├── CAD           analytical surfaces, projection, and imprint curves
└── HighOrder     globally certified quadratic tetrahedra and type-11 I/O
```

`FunctionSize` is the extension point for distance-, curvature-, solution-, and
background-mesh-driven sizing; `MinSize` combines such fields conservatively. This
keeps field policy outside the meshing kernels while preserving a checked `h > 0`
contract.

## Design contracts

1. Exactness is used for topology-changing geometric decisions. Floating filters are
   permitted only with exact fallback or an exact post-certificate.
2. A public meshing operation returns a validated result or a precise blocker. It
   never reports a silent empty region, partial complex, or nonconforming fallback.
3. Boundary recovery is independently certified against the input PLC, including
   manifold vertex links and every required facet/interface.
4. Size-control operations verify their output bound and preserve lower-dimensional
   cells, physical tags, conformity, and region volume.
5. Resource counts are checked before conversion/allocation; connectivity uses Int32,
   compact flat records, reusable scratch, and bounded retry loops.
6. File output is validated before an atomic replacement of the destination.

## Completed stages

| Stage | Exit condition | State |
|---|---|---|
| 0 | exact predicates, mesh validation/CRC, MSH/STL round trip | DONE |
| 1 | arbitrary PSLG Delaunay/CDT and bounded quality refinement | DONE |
| 2 | graded curves and planar/cylinder/parametric surfaces | DONE |
| 3 | 3-D Delaunay, volume fill, exact conforming recovery, enclosure case | DONE |
| 4 | uniform/callback sizing, quality reports, flips and sliver reduction | DONE |
| 5 | native primitives, analytical CAD/imprints, healing gates, mesh CSG | DONE |
| 6 | globally certified P2 elements, solver I/O, 22 geometry regressions | DONE |

The external full-wave rerun of the remaining HFSS guide studies is deliberately
tracked in [`ASCENT.md`](ASCENT.md), because it requires the solver and proprietary
reference artifacts and does not add a missing meshing capability.

## Verification discipline

Every change follows:

`spec → independent oracle → implementation → CRC test → bounds-checked recheck`

Independent evidence includes exact-rational predicates, analytic area/volume,
Euler/manifold/link invariants, empty-circle/sphere checks, PLC facet conservation,
quality monotonicity, MSH round trips, and gmsh cross-validation. The mandatory gate
is:

```sh
julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

See [`DEVELOPMENT.md`](DEVELOPMENT.md) for the non-negotiable anti-false-positive
rules and [`STATUS.md`](STATUS.md) for the last measured gate.

## Standing acceptance case

The enclosure/coax feed-through remains the reason-to-exist regression. gmsh 4.13.1
and 4.15.2-git leave its solid regions empty; Tessella reconstructs a conforming,
tagged, all-filled native mesh which is loadable and solvable in ASCENT. The package
regression is in `test/` and the external proof chain is in `validation/` and
[`ASCENT.md`](ASCENT.md).
