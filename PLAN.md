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
├── IO            strict/atomic MSH v2.2/v4.1, STL, bounded .geo constant/range scan
├── Elements      fixed/special Gmsh catalog, mixed metadata, ASCII/binary MSH I/O
├── Recombine     greedy or Edmonds-blossom triangle-to-quad pairing, optional full-quad
├── Refine        deterministic one-level uniform linear-simplex refinement
├── Transfinite   validated four-sided planar structured triangle patches
├── TransfiniteCurve normalized straight-curve Progression/Bump/Beta laws
├── TransfiniteTriangle validated three-sided structured triangle patches
├── TransfiniteQuad recombined four-sided structured quadrangle patches
├── TransfiniteVolume affine six-face structured tetrahedral volumes
├── TransfinitePrism affine five-face structured tetrahedral prisms
├── TransfiniteHex affine six-face recombined hexahedral volumes
├── Mesh2D        Delaunay, CDT, interior classification, quality refinement
├── SizeField     scalar/anisotropic catalog, .geo field graph, context resolvers
├── Mesh1D        metric-length curve and segment discretization
├── MeshSurface   planar, cylindrical, and parametric surface meshing
├── Mesh3D        Delaunay, fills, partitions, sizing, flips, mesh Boolean CSG
├── RecoverCDT    exact conforming-Delaunay boundary/partition recovery
├── Optimize      quality reports, Laplacian/ODT/targeted sliver smoothing
├── Heal          surface defect and meshability diagnostics
├── Geometry      native box/cylinder/cone/geodesic-sphere surfaces
├── CAD           analytical surfaces, projection, imprints
├── BRep          ISO-10303-21 STEP / IGES classified-solid import (box/sphere/cylinder/cone) plus NURBS 126/128
├── NURBS         native B-spline/NURBS curve and surface evaluation
├── Model         tagged point/curve/loop/surface/volume entity kernel plus native solids/Booleans/embeds
├── GeoExec       bounded Point/Line/Loop/Surface/Box/Cylinder/Sphere/Cone/Boolean/Translate/Dilate/Rotate `.geo` execution
├── BoundaryLayer prismatic 3-D extrusion, certified remaining-core tet fill,
│                 plus 2-D quad/fan topology
├── Periodic      translation periodic identification
├── Post          list-based views and plugins
├── API           model/mesh/option façade
├── CLI           `tessella file.geo -2|-3` entry
├── GUI           headless command/state machine
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
| P2 | general entity model and every Gmsh element family/order in memory and MSH I/O | IN PROGRESS — 125 fixed-node types plus special records, mixed MSH I/O, and a tagged point/curve/surface/volume kernel |
| P3 | built-in/OCC-equivalent CAD, BREP/NURBS, imports, Booleans, transforms, `.geo` execution | IN PROGRESS — NURBS evaluation and STEP/IGES NURBS import (B_SPLINE / IGES 126/128) with IGES export, classified STEP/IGES box/sphere/cylinder/cone solids, Box/Cylinder/Sphere/Cone/Boolean/Translate/Dilate/90°-Rotate `.geo` execution, mesh Booleans/transforms; unrecognized CAD topology remains an explicit blocker |
| P4 | structured/unstructured algorithms, recombination, layers, adaptation, periodic/embedded constraints | IN PROGRESS — plus blossom/full-quad surface pairing, Point/Line-In-Surface embeddings, Point/Line/Surface-In-Volume recovery, holed plane surfaces, recombined hexahedra, prismatic 3-D layers with certified remaining-core tet fill and cavity walls, 2-D quad/fan layers, and translation periodic identification |
| P5 | complete API/options/formats, partitioning/parallel paths, views/plugins, CLI/GUI/post-processing | IN PROGRESS — model/mesh API, CLI, headless GUI state, views/plugins |
| P6 | tutorial/API corpus and requirement-by-requirement differential conformance to Gmsh 4.15.2 | IN PROGRESS — size-field/transfinite/range differentials plus t1 square, t4 hole, Point/Line-In-Surface, Surface-In-Volume sheet, 2-D boundary-layer quads, API box, OCC cylinder/cone, IGES-128 bilinear patch, and BooleanDifference box corpus |

P1 does not yet claim 3-D multi-wall boundary-layer fans, Gmsh's global
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
meshes. Classified ISO-10303-21 STEP and IGES solids that are axis-aligned blocks,
spheres, right circular cylinders, or right circular cones are imported and filled;
STEP `B_SPLINE_CURVE_WITH_KNOTS` / `B_SPLINE_SURFACE_WITH_KNOTS` (including
complex rational instances) and IGES 126/128 import as native NURBS objects, with
IGES NURBS export. Any other topology is
an explicit blocker listing the seen entity types. The entity kernel records native
boxes, cylinders, spheres, cones, and mesh-Boolean volumes. Bounded `.geo` execution covers
Point/Line/Loop/Surface, Box/Cylinder/Sphere/Cone, BooleanDifference/Union/Intersection
of those solids (with operand `Delete`), Translate of remaining native solids,
Dilate about a center, and coordinate-axis rotations by integer multiples of π/2
that preserve axis-aligned boxes. Point-In-Surface embeddings force the classified
point to appear as a mesh node; Line-In-Surface embeddings recover the curve as an
interior constrained edge chain without creating a hole; Point-In-Volume embeddings
insert a Steiner vertex by a tet split; Line-In-Volume and Surface-In-Volume recover
interior edge chains and open triangular sheets as tet faces without removing volume.
It does not yet claim a general OpenCASCADE BREP kernel, NURBS CAD of
unclassified topology, transformations of arbitrary CAD entities, or complete `.geo`
execution. The `.geo` scanner evaluates finite arithmetic constants, pure numeric
functions, prior scalar bindings, and explicit field/physical tags with resource
bounds. Finite constant `start:end[:increment]` lists are expanded in recognized
numeric field options and field selectors; entirely numeric Physical memberships are
range-checked but remain geometry data. The scanner deliberately rejects loops,
macros, dynamic tag allocators, option reads, stateful functions, dynamic/general
ranges, logical/ternary evaluation, extrusions/fillets/symmetry, and mixed
geometry-derived physical right-hand-side evaluation instead of pretending to be a
complete interpreter.

P4 currently covers validated, deterministic pairing of adjacent same-physical-tag
surface triangles into first-order quadrangles, with unpaired triangles and boundary
segments retained. `:greedy` accepts candidates in shape-score order; `:blossom` runs
Edmonds maximum-cardinality matching on the dual of eligible pairs, and `full_quad`
requires a perfect matching. It also covers one-level uniform refinement of linear segments,
triangles, and tetrahedra with Gmsh 4.15.2 child templates, shared lexicographic edge
midpoints, compacted unused nodes, and parent-tag preservation. Normalized affine-line
transfinite parameters cover Gmsh's Progression/Power, Bump, and Beta laws with signed
orientation and representability gates. Three-sided and four-sided planar transfinite
patches implement Gmsh's specific triangular and average-chord Coons interpolation
for already-discretized, count-matched boundary chains. Four-sided grids can also be
emitted as Gmsh-compatible first-order quadrangles with exact projected
corner-Jacobian certification. Positively ordered affine eight-corner blocks
implement Gmsh's unrecombined six-tetrahedron transfinite volume subdivision;
canonical affine triangular prisms implement Gmsh's legacy collapsed-grid five-face
tetrahedral path; positively ordered affine eight-corner blocks can also be emitted
as first-order recombined hexahedra with type-3 boundary quadrangles. Planar
constant-`z` polylines extrude to type-3 quadrangles along left-normals, with
optional convex-corner fans of first-layer triangles and subsequent ring
quadrangles. P4 does not yet claim
non-affine CAD curve integration, FlexibleTransfinite
or HWall/size-map curve laws,
quasi-transfinite patches, general CAD parameterizations,
curved/warped or compact-TransfiniteTri volumes, recombined three-sided patches,
volume/hybrid recombination, selective or
high-order refinement, coarsening, or
3-D multi-wall boundary-layer fans. The filled extrusion
(`mesh_boundary_layer_filled`) certifies the remaining core with per-wall shell
and global fill volume identities and an interface tiling gate; its core engine
covers Delaunay-friendly caps (planar/primitive walls) directly and smaller
smooth caps through bounded exact-rational recovery — larger smooth caps are an
explicit blocker naming both stages rather than a defective mesh.

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

The aggregate validation launches the Gmsh 4.15.2 size-field, constant-range,
uniform-refinement,
four-sided transfinite, straight transfinite curve-law, three-sided transfinite,
recombined-quadrangle, affine transfinite-volume, five-face-prism, and
recombined-hexahedron differentials
as required bounds-checked children, together with the P6 t1-square, OCC-cylinder,
OCC-cone, IGES-128 bilinear, BooleanDifference box, and 2-D boundary-layer differentials. A
missing or wrong-version Gmsh runtime, a failed probe, or a parity mismatch makes the
aggregate command fail.

See [`DEVELOPMENT.md`](DEVELOPMENT.md) for the non-negotiable anti-false-positive
rules and [`STATUS.md`](STATUS.md) for the last measured gate.

## Standing acceptance case

The enclosure/coax feed-through remains the reason-to-exist regression. gmsh 4.13.1
and 4.15.2-git leave its solid regions empty; Tessella reconstructs a conforming,
tagged, all-filled native mesh which is loadable and solvable in ASCENT. The package
regression is in `test/` and the external proof chain is in `validation/` and
[`ASCENT.md`](ASCENT.md).
