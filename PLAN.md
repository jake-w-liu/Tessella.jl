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
├── IO            strict/atomic MSH v2.2/v4.1, STL, bounded .geo constant scan
├── Elements      fixed/special Gmsh catalog, mixed metadata, ASCII/binary MSH I/O
├── Recombine     deterministic physical-tag-preserving triangle-to-quad pairing
├── Refine        deterministic one-level uniform linear-simplex refinement
├── Transfinite   validated four-sided planar structured triangle patches
├── TransfiniteCurve normalized straight-curve Progression/Bump/Beta laws
├── TransfiniteTriangle validated three-sided structured triangle patches
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
| P2 | general entity model and every Gmsh element family/order in memory and MSH I/O | IN PROGRESS — 125 fixed-node types plus ten serializable cut/border/child/sub-element records, mixed metadata, validation/CRC, and ASCII/binary MSH v2.2/v4.1 shipped |
| P3 | built-in/OCC-equivalent CAD, BREP/NURBS, imports, Booleans, transforms, `.geo` execution | IN PROGRESS — native analytical surfaces/imprints, closed primitives, mesh Booleans, finalized-mesh affine transforms, and bounded `.geo` constant expressions shipped |
| P4 | structured/unstructured algorithms, recombination, layers, adaptation, periodic/embedded constraints | IN PROGRESS — recombination, uniform refinement, straight-curve laws, and three-/four-sided planar transfinite patches shipped |
| P5 | complete API/options/formats, partitioning/parallel paths, views/plugins, CLI/GUI/post-processing | PENDING |
| P6 | tutorial/API corpus and requirement-by-requirement differential conformance to Gmsh 4.15.2 | PENDING |

P1 does not yet claim boundary-layer element topology, Gmsh's global
`AutomaticMeshSizeField` pipeline, high-order/custom-interpolation,
multiple-time-step, or mixed-component `PostView` data, materially warped
quadrangles, `PostView` tensor-to-metric evaluation,
direct tensor or metric-meshing parity, full `.geo`/CAD-model execution, or exact CAD
distance queries. P2 does not yet claim general mixed-element generation or
recombination beyond P4's first-order surface pairing, integration of mixed blocks
into the simplex meshing kernels, basis-selector tags 138/139 as mesh records, curved
high-order Jacobian certification, preservation of
ancillary/unknown MSH sections (binary readers reject unsupported sections
explicitly), repeated `$Nodes` sections, non-8-byte binary data, internal indexing
beyond `Int32`, or lossless multi-physical-group projection through MSH v2.2.
Variable-connectivity types 34/35/69 and parent/domain links are lossless in MSH2
ASCII. Binary MSH2 has fixed widths and supports only fixed special records and parent
links; MSH4 supports fixed unlinked special records. Type 69 and some registered fixed
tags require explicit Tessella-only output because Gmsh 4.15.2 cannot consume them
safely. Pinned Gmsh 4.15.2 also corrupts distinct parent links when it rewrites MSH2
binary, while Tessella's binary round trip and Gmsh's ASCII rewrite preserve them.
Nonzero-physical special MSH4 output requires compatible node/entity classification
metadata for Gmsh-safe rewrites. P3 currently provides the native analytical and
polyhedral path used by ASCENT, including boxes, cylinders, cones, geodesic spheres,
projection/imprint curves, cavities, mesh Booleans, and validated translation,
rotation, dilation, reflection, and general affine transforms of finalized simplex
meshes. It does not yet claim a general entity kernel, OpenCASCADE/BREP/NURBS, CAD
import/export, transformations of analytical/CAD entities, or complete `.geo`
execution. The `.geo` scanner evaluates finite arithmetic constants, pure numeric
functions, prior scalar bindings, and explicit field/physical tags with resource
bounds. It deliberately rejects loops, macros, dynamic tag allocators, option reads,
stateful functions, ranges, logical/ternary syntax, CSG statements, and physical-group
right-hand-side evaluation instead of pretending to be a complete interpreter.

P4 currently covers validated, deterministic pairing of adjacent same-physical-tag
surface triangles into first-order quadrangles, with unpaired triangles and boundary
segments retained. It also covers one-level uniform refinement of linear segments,
triangles, and tetrahedra with Gmsh 4.15.2 child templates, shared lexicographic edge
midpoints, compacted unused nodes, and parent-tag preservation. Normalized affine-line
transfinite parameters cover Gmsh's Progression/Power, Bump, and Beta laws with signed
orientation and representability gates. Three-sided and four-sided planar transfinite
patches implement Gmsh's specific triangular and average-chord Coons interpolation
for already-discretized, count-matched boundary chains. P4 does not yet claim Gmsh's
Blossom/full-quad algorithms, non-affine CAD curve integration, FlexibleTransfinite
or HWall/size-map curve laws,
quasi-transfinite or holed patches, general CAD parameterizations,
transfinite quadrangles or volumes, volume/hybrid recombination, selective or
high-order refinement, coarsening, boundary-layer element topology, or
periodic/embedded model constraints. P5–P6 remain pending.

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

The aggregate validation launches the Gmsh 4.15.2 size-field, uniform-refinement,
four-sided transfinite, straight transfinite curve-law, and three-sided transfinite
differentials as required bounds-checked children. A missing or wrong-version Gmsh
runtime, a failed probe, or a parity mismatch makes the aggregate command fail.

See [`DEVELOPMENT.md`](DEVELOPMENT.md) for the non-negotiable anti-false-positive
rules and [`STATUS.md`](STATUS.md) for the last measured gate.

## Standing acceptance case

The enclosure/coax feed-through remains the reason-to-exist regression. gmsh 4.13.1
and 4.15.2-git leave its solid regions empty; Tessella reconstructs a conforming,
tagged, all-filled native mesh which is loadable and solvable in ASCENT. The package
regression is in `test/` and the external proof chain is in `validation/` and
[`ASCENT.md`](ASCENT.md).
