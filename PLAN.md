# Tessella.jl architecture and scope

Tessella is a Julia-native mesh generator for the ASCENT electromagnetics workflow.
It is an independent implementation informed by Gmsh's architecture. The original
simplex-mesher roadmap is complete; the active goal is now full Gmsh 4.15.2 feature
and behavioral parity, prioritized by ASCENT meshing value. The parity goal is not
complete. The live verification record is [`STATUS.md`](STATUS.md).

## Target scope

The parity end state includes:

- Gmsh's geometry/entity model, built-in geometry kernel, OpenCASCADE-facing CAD
  capabilities, transformations, Boolean operations, imports, and `.geo` language;
- all supported element families, orders, structured/unstructured algorithms,
  recombination, boundary layers, adaptation, partitioning, periodicity, embedding,
  sizing fields, and optimization operations;
- compatible model/mesh APIs, options, physical groups, file formats, views, plugins,
  command-line behavior, GUI/post-processing, and parallel workflows;
- differential conformance against the pinned Gmsh release, in addition to Tessella's
  stronger exact-predicate and validated-or-explicit-blocker contracts.

This supersedes the former non-goal boundary around OpenCASCADE/NURBS, mixed elements,
GUI/post-processing, and long-tail formats. Those capabilities are now pending work,
not exclusions. Reimplementation remains independent; calling Gmsh as the production
meshing backend would not resolve the ASCENT failure that motivated Tessella.

## Implemented architecture

```text
Tessella
├── Predicates    adaptive exact orient/incircle/insphere, exact rationals, SoS
├── MeshTypes     compact simplex storage, topology, quality, CRC, validation
├── Transform     validated affine transforms for finalized simplex meshes
├── ExactMesh3D   Rational{BigInt} Delaunay kernel
├── IO            strict/atomic simplex MSH v2.2/v4.1, STL, limited .geo scan
├── Elements      fixed-node Gmsh catalog, mixed metadata, ASCII/binary MSH I/O
├── Mesh2D        Delaunay, CDT, interior classification, quality refinement
├── SizeField     scalar/anisotropic catalog, .geo field graph, context resolvers
├── Mesh1D        metric-length curve and segment discretization
├── MeshSurface   planar, cylindrical, and parametric surface meshing
├── Mesh3D        Delaunay, fills, partitions, sizing, flips, mesh Boolean CSG
├── RecoverCDT    exact conforming-Delaunay boundary/partition recovery
├── Optimize      quality reports, Laplacian/ODT/targeted sliver smoothing
├── Heal          surface defect and meshability diagnostics
├── Geometry      native box/cylinder/cone/geodesic-sphere surfaces
├── CAD           analytical surfaces, projection, and imprint curves
└── HighOrder     globally certified quadratic tetrahedra and type-11 I/O
```

The native field graph includes geometric/composite, analytic and derivative,
coordinate-map, structured-grid, entity-aware, sampled-view, anisotropic,
boundary-layer scalar-law, discrete automatic-sizing analogue, and external-process
fields. Discrete point/segment/triangle distance queries use a deterministic AABB
hierarchy. `FunctionSize` remains the checked callback extension point. Generic fields
may return zero (for example, distance on a target); only `AbstractSizeField` reaches a
meshing kernel, where `size_at` enforces a finite `h > 0` contract.

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

## Completed baseline stages

| Stage | Exit condition | State |
|---|---|---|
| 0 | exact predicates, mesh validation/CRC, MSH/STL round trip | DONE |
| 1 | arbitrary PSLG Delaunay/CDT and bounded quality refinement | DONE |
| 2 | graded curves and planar/cylinder/parametric surfaces | DONE |
| 3 | 3-D Delaunay, volume fill, exact conforming recovery, enclosure case | DONE |
| 4 | uniform/callback sizing, quality reports, flips and sliver reduction | DONE |
| 5 | native primitives, analytical CAD/imprints, healing gates, mesh CSG | DONE |
| 6 | globally certified P2 elements, solver I/O, 22 geometry regressions | DONE |

## Gmsh parity tracks

| Track | Exit condition | State |
|---|---|---|
| P1 | full scalar/isotropic/anisotropic field catalog and field-driven 1-D/2-D/3-D sizing | IN PROGRESS — native catalog, strict field graph, and entity-aware mesher integration shipped |
| P2 | general entity model and every Gmsh element family/order in memory and MSH I/O | IN PROGRESS — 125 fixed-node types, mixed metadata, validation/CRC, and ASCII/binary MSH v2.2/v4.1 shipped |
| P3 | built-in/OCC-equivalent CAD, BREP/NURBS, imports, Booleans, transforms, `.geo` execution | IN PROGRESS — native analytical surfaces/imprints, closed primitives, mesh Booleans, and finalized-mesh affine transforms shipped |
| P4 | structured/unstructured algorithms, recombination, layers, adaptation, periodic/embedded constraints | PENDING |
| P5 | complete API/options/formats, partitioning/parallel paths, views/plugins, CLI/GUI/post-processing | PENDING |
| P6 | tutorial/API corpus and requirement-by-requirement differential conformance to Gmsh 4.15.2 | PENDING |

P1 does not yet claim boundary-layer element topology, Gmsh's global
`AutomaticMeshSizeField` pipeline, non-simplex/higher-order/vector/tensor `PostView`,
direct tensor or metric-meshing parity, full `.geo`/CAD-model execution, or exact CAD
distance queries. P2 does not yet claim element generation/recombination,
integration of mixed blocks into the simplex meshing kernels, variable-connectivity or
internal element types, curved high-order Jacobian certification, preservation of
ancillary/unknown MSH sections (binary readers reject unsupported sections
explicitly), repeated `$Nodes` sections, non-8-byte binary data, internal indexing
beyond `Int32`, or lossless multi-physical-group projection through MSH v2.2.
Some registered fixed tags also require explicit Tessella-only output because Gmsh
4.15.2 cannot re-import them. P3 currently provides the native analytical and
polyhedral path used by ASCENT, including boxes, cylinders, cones, geodesic spheres,
projection/imprint curves, cavities, mesh Booleans, and validated translation,
rotation, dilation, reflection, and general affine transforms of finalized simplex
meshes. It does not yet claim a general entity kernel, OpenCASCADE/BREP/NURBS, CAD
import/export, transformations of analytical/CAD entities, or complete `.geo`
execution. P4–P6 remain pending.

The external HFSS solve campaign remains tracked in [`ASCENT.md`](ASCENT.md); it is a
consumer-side validation track, not a substitute for the parity work above.

## Verification discipline

Every change follows:

`spec → independent oracle → implementation → CRC test → bounds-checked recheck`

Independent evidence includes exact-rational predicates, analytic area/volume,
Euler/manifold/link invariants, empty-circle/sphere checks, PLC facet conservation,
quality monotonicity, MSH round trips, and Gmsh cross-validation. The mandatory gates
are:

```sh
julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'
julia --project --check-bounds=yes validation/run_all.jl
```

The aggregate validation launches the Gmsh 4.15.2 size-field differential as a
required bounds-checked child. A missing or wrong-version Gmsh runtime, a failed probe,
or a parity mismatch makes the aggregate command fail.

See [`DEVELOPMENT.md`](DEVELOPMENT.md) for the non-negotiable anti-false-positive
rules and [`STATUS.md`](STATUS.md) for the last measured gate.

## Standing acceptance case

The enclosure/coax feed-through remains the reason-to-exist regression. gmsh 4.13.1
and 4.15.2-git leave its solid regions empty; Tessella reconstructs a conforming,
tagged, all-filled native mesh which is loadable and solvable in ASCENT. The package
regression is in `test/` and the external proof chain is in `validation/` and
[`ASCENT.md`](ASCENT.md).
