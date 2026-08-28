# Tessella.jl

Tessella is a Julia-native tetrahedral mesh generator built around exact geometric
predicates, certified topology, conforming boundary recovery, and explicit failure
diagnostics. It was created for the ASCENT electromagnetics workflow after gmsh
4.13.1 and 4.15.2-git left the enclosure/coax acceptance geometry with zero volume
elements.

The original simplex-mesher roadmap is complete through Stage 6. Development has now
expanded toward independent Gmsh 4.15.2 feature and behavioral parity, with
ASCENT-relevant meshing capabilities implemented first. That parity target is **not
complete**. Bounds-checked package and focused gates for each implementation
increment are recorded in [`STATUS.md`](STATUS.md).
The separate external ASCENT solve campaign is recorded in [`ASCENT.md`](ASCENT.md).

```julia
using Tessella
using Tessella.IO: write_msh

m = mesh_volume(surface)                 # closed triangle surface → validated tets
ms = mesh_sized(surface; hmax=0.5)       # certified maximum-edge size bound
near = DistanceField(surface)
field = ThresholdField(near; dist_min=0.0, dist_max=0.2,
                       size_min=0.02, size_max=0.5)
mf = mesh_sized(surface; field=field)    # spatially graded volume mesh
write_msh("mesh.msh", ms; version=4.1)   # solver-consumable gmsh MSH
```

## Implemented scope

- adaptive exact 2-D/3-D predicates with one consistent Simulation-of-Simplicity
  convention and exact-rational test oracles;
- compact simplex meshes with finite-input, cell, tag, manifold-link, and quality
  validation;
- 2-D Delaunay/CDT and quality refinement, graded curves, and planar, cylindrical,
  and parametric surface meshing;
- Float64 and exact-coordinate 3-D Delaunay kernels, constrained boundary recovery,
  multi-region partitions, and explicit recovery blockers;
- uniform and field-driven 1-D, 2-D, surface, and 3-D sizing, including geometric and
  composite (`Distance`, `Threshold`, `Box`, `Ball`, finite `Cylinder`, `Frustum`,
  `Min`, `Max`, bounds), analytic/derivative (`MathEval`, `Gradient`, `Laplacian`,
  `Mean`, `Curvature`, `MaxEigenHessian`), coordinate/storage (`LonLat`, `Parametric`,
  `Structured`), entity-aware (`Restrict`, `Constant`), sampled/context (`Extend`,
  `Octree`, first-order scalar/vector `PostView` on every standard point/line/surface/
  volume list element), scalar `BoundaryLayer`,
  discrete `AutomaticMeshSizeField` analogue, and `ExternalProcess` fields;
- anisotropic metric primitives and `MathEvalAniso`, `MinAniso`, `IntersectAniso`, and
  `AttractorAnisoCurve`; discrete distance queries use a deterministic AABB hierarchy,
  and lower-dimensional cells and tags are preserved through volume refinement;
- quality reporting, flips, Laplacian/ODT/targeted sliver smoothing, healing
  diagnostics, native box/cylinder/cone/geodesic-sphere primitives, analytical
  surfaces, imprints, mesh Boolean CSG, and validated affine transformations of
  finalized simplex meshes with orientation and physical-tag preservation;
- globally certified quadratic tetrahedra, plus strict and atomic simplex MSH
  v2.2/v4.1 and STL I/O;
- a resource-bounded `.geo` scanner for finite arithmetic constants, pure numeric
  functions, prior scalar bindings, sizing options, explicit field/physical tags, and
  finite constant `start:end[:increment]` lists in recognized numeric contexts, plus
  read-only tag allocators over tracked Point, shared-region, and Field namespaces and
  factory-aware, checked `SetMaxTag` control of geometric counters; the
  executor applies the same checked semantics to numeric parameters, entity tags,
  and entity lists in every supported geometry statement, and stores positive
  `MeshSize`/`Characteristic Length` constraints on existing explicit Points,
  including recursive explicit-topology `PointsOf` selection;
- a 125-type fixed-node Gmsh catalog plus ten serializable cut/border/child/
  sub-element records, compact variable connectivity and parent/domain references,
  mixed entity/classification metadata, structural validation/CRC, and ASCII/binary
  MSH v2.2/v4.1 mixed-element I/O with opposite-endian decoding and cumulative
  repeated pre-element node sections, plus owned periodic metadata in ASCII and
  binary MSH2/MSH4 files, including entity transforms, node correspondences, and
  declared MSH2 elementary-entity tags; positive, physically consistent MSH2
  entities can be converted to discrete MSH4 entities without changing their
  entity tags;
- deterministic, physical-tag-preserving surface triangle-to-quadrangle
  recombination with convexity, resource-growth, CRC, and Gmsh-load gates;
- oriented-normal extrusion of coherently directed planar polylines into
  first-order quadrangle strips with optional convex-corner fans, exact projected
  corner-Jacobian checks, resource bounds, and rigid-motion certification;
- translation and general finite nonsingular affine periodic node-pair
  certification/snapping using Gmsh's row-major 4×4 convention, with node numbers,
  connectivity, and physical tags preserved; the native model/API owns straight-
  curve relations with one master per slave, reusable masters, and acyclic
  master/slave chains; it synchronizes planar subdivisions for boundary or
  embedded curves and returns detached node maps; it also owns affine-equivalent
  planar boundary-surface pairs on explicit volume shells, synchronizes each slave
  triangulation from its master, and certifies the tetrahedron-face map;
  `model_to_mixed` emits
  classified point/line/triangle blocks, signed hole boundaries, embedded point
  and curve cells, physical memberships, and periodic curve/endpoint
  links for planar surfaces; mixed meshes retain periodic maps and Gmsh's MSH4
  curve-in-surface relation, and the CLI uses this path for periodic or embedded
  `-2` output; the dimension-explicit form also emits classified boundary and
  Point/Line/Surface-In-Volume cells after certifying the selected solid fill,
  including nested Point/Line-In-Surface constraints, holed planar sheets, and
  explicit planar surface-loop volumes with cavity signs; periodic explicit volumes
  retain their derived boundary point/curve forest and surface-node maps through
  MSH2/MSH4 output, and the CLI uses this path for classified `-3` output;
- one-level uniform linear-simplex refinement using Gmsh 4.15.2's ordered 2/4/8
  child templates, shared deterministic edge midpoints, tag preservation, resource
  bounds, and output validation;
- validated four-sided planar transfinite triangle patches using average-chord Coons
  interpolation, all four Gmsh diagonal arrangements, physical tags, bounded
  intersection auditing, and exact orientation postconditions;
- normalized straight-curve transfinite parameters for Gmsh's Progression/Power,
  Bump, and Beta laws plus all three HWall variants, with signed orientation and
  Float64 representability gates;
- validated three-sided planar structured patches using Gmsh's specific triangular
  interpolation, compact unrecombined topology, all four recombined
  triangle/quadrangle arrangements, and exact geometry/boundary postconditions;
- recombined four-sided planar transfinite patches with Gmsh type-3 quadrangles,
  physical tags, and exact projected corner-Jacobian certification;
- affine six-face transfinite volumes using Gmsh's unrecombined six-tetrahedron
  subdivision, with boundary tags, exact topology/orientation postconditions,
  exact-dyadic remote-grid interpolation, and represented-volume certification;
- affine five-face triangular prisms using Gmsh's legacy collapsed-grid tetrahedron
  templates, with exact orientation, boundary, and represented-volume certification.

P1 through P4 remain **in progress**. Current non-claims include 3-D multi-wall
boundary-layer fan topology beyond the certified closed-wall extrusion, the full
Gmsh automatic-sizing pipeline, broader `PostView` data including high-order/custom
interpolation, materially warped quadrangles, mixed component counts, and
tensor-to-metric evaluation, general OpenCASCADE/unclassified NURBS CAD, and full
`.geo` execution (including control-flow loops, macros, option reads,
dynamic/general ranges, allocator reads after topology-changing or untracked
declarations, and geometry-derived physical-group RHSs beyond the documented inline
topology queries),
mixed-element generation beyond the listed first-order surface recombination paths,
non-affine CAD curve integration, FlexibleTransfinite, and size-map laws,
quasi-transfinite or holed transfinite patches, curved/warped or
compact-TransfiniteTri volumes,
selective/high-order refinement, simplex-kernel integration,
MINI basis-selector tags 138/139 as mesh records, curved-cell Jacobian certification,
non-8-byte binary data, and ancillary-section preservation. MSH2 ASCII is the lossless
format for variable connectivity and parent/domain links; binary MSH2 and MSH4 have
explicitly narrower special-record contracts. Type 69 and some registered fixed tags
require Tessella-only output because Gmsh 4.15.2 cannot consume them safely, and
nonzero-physical special MSH4 requires compatible classification metadata for a
Gmsh-safe rewrite. MSH2 retains cell-level elementary ownership but has no entity
topology record for signed boundaries or embedding relations. Native model periodicity
covers straight-curve pairs on one planar triangle-meshed surface and disjoint
affine-equivalent planar boundary surfaces of one explicit surface-loop volume. Each
slave has one master; curve masters may be reused, a curve slave may become a master
in an acyclic chain, and independent relations may share corner points. The bounded
`.geo` executor applies prior scalar bindings, finite arithmetic, pure numeric
functions, bounded numeric list assignment/indexing/selection/mutation, and constant
entity ranges to all supported geometry statements. Entity lists, including periodic
slave/master sets, can reuse whole or selected list variables.
Point `MeshSize` selectors can likewise use `:`, expressions, constant ranges, and
whole or selected numeric-list variables. Inline `PointsOf` blocks additionally
select the recursive boundary Points of explicit Point, Curve/Line, Surface, and
Volume entities; signed entity tags are normalized, hole boundaries participate,
and embeddings do not. Tessella requires existing explicit Points and finite positive
sizes; Gmsh's permissive zero/negative and missing-Point behavior is outside this
bounded contract. The planar surface path extends Point constraints
piecewise-linearly over its deterministic initial constrained triangulation. Generated
straight-curve subdivision nodes interpolate their endpoint sizes, while equal Point
sizes retain the constant-size path and coincident PSLG inputs use the smaller
constraint. Exact Gmsh mesh topology, implicit primitive or Boolean subentities, and
mesh-size selectors other than inline `PointsOf` remain outside this bounded contract.
Physical declarations accept an explicit positive tag, with an optional name, or a
nonempty name with an automatic tag. Automatic Physical tags share one namespace
across entity dimensions and advance later shared-region allocator reads. Physical
Point accepts inline `PointsOf`; Physical Point/Curve/Surface accept inline `Boundary` and
`CombinedBoundary` over Curve/Line, Surface, and explicit Volume entities,
respectively. `Boundary` collects immediate boundaries before physical membership is
deduplicated; `CombinedBoundary` keeps tags with odd multiplicity. Hole and cavity
boundaries participate, while embeddings do not. Empty combined boundaries,
unsupported dimensions, and implicit primitive or Boolean topology are explicit
blockers.
The model and session APIs return detached, sorted group, membership, reverse-membership,
and name-query results. Names are unique within an entity dimension and can be removed
independently of their groups. Selective or all-group removal leaves geometry intact
and does not rewind the global automatic Physical-tag counter.
They also enumerate explicit entities, report model dimension, and return direct or
recursive boundaries and direct adjacencies with deterministic Gmsh-compatible
ordering, orientation, and combined-incidence cancellation. Topology queries preserve
the session mesh cache and exclude embeddings. Primitive and Boolean volumes are
enumerated, but their implicit boundary topology remains an explicit blocker.
It accepts `Periodic Line`, `Periodic Curve`, and `Periodic Surface` with
`Translate`, `Rotate`, and 12- or 16-entry `Affine` transforms. Cyclic curve
dependencies, curved or non-boundary periodic surfaces, periodic volume entities,
and allocator reads after topology-changing or untracked declarations remain
pending. Gmsh 4.15.2 has no
serialized Point-In-Surface or Point/Line/Surface-In-Volume relation:
the classified nodes and elements remain available instead. Its ASCII
MSH4 reader reconstructs curve embeddings, while its binary reader drops that
relation with a warning; Tessella reads both variants losslessly. Classified
volume-embedding and explicit-shell entities and cells, including nested
point/curve/surface classification, survive MSH2/MSH4 output. MSH4 retains signed
surface boundaries for explicit planar volumes; MSH2 retains cell ownership but not
the entity topology.
Complete CAD/BREP, broad file/API compatibility,
GUI, and post-processing remain pending. See
[`PLAN.md`](PLAN.md) rather than treating the completed Stage 0–6 baseline as Gmsh
completeness. ASCENT remains Tessella's primary solver consumer.

## Verification

```sh
julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'
julia --project --check-bounds=yes validation/run_all.jl
```

The first command runs the full package/CRC suite. The second compares analytic
volumes and quality with Gmsh, reproduces the enclosure/coax failure, and runs the
Gmsh 4.15.2 size-field, constant-range, geometry-expression, numeric-list,
spatial Point-size, dynamic-tag/`SetMaxTag`, explicit model-topology,
uniform-refinement, transfinite-patch,
straight transfinite curve-law/HWall, unrecombined/recombined three-sided
transfinite, recombined-quadrangle, affine
transfinite-volume, five-face-prism, native `.geo` and projected single-/two-direction
periodic surface differentials, embedded, reusable-master/chained, and
expression/list-backed periodic curves, low-level translation/rotation-periodic curves,
planar periodic boundaries of an explicit volume, and the MSH2/MSH4 lifecycle, classified
Point/Line-In-Surface MSH4 round trips,
and the classified Surface-In-Volume lifecycle with nested point/curve constraints
and optional holes as mandatory bounds-checked children.
Missing or wrong-version Gmsh and differential failures make the aggregate command
fail. See
[`DEVELOPMENT.md`](DEVELOPMENT.md) for the correctness, robustness, and completeness
rules.

## Repository map

- [`src/`](src/) — package entry point plus domain folders: `core`, `fields`,
  `geometry`, `meshing`, `structured`, and `interfaces`
- [`test/`](test/) — package runner plus matching `core`, `fields`, `geometry`,
  `meshing`, `structured`, and `interfaces` suites, with end-to-end cases in
  `integration`
- [`validation/`](validation/) — workflow-specific external-tool and solver-facing
  fixtures, with shared harness code in `support`
- [`PLAN.md`](PLAN.md) — architecture and scope boundary
- [`STATUS.md`](STATUS.md) — current package verification record
- [`ASCENT.md`](ASCENT.md) — external solver integration and HFSS campaign

## License

No license has been declared by the repository owner.
