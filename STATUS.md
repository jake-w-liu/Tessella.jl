# Tessella.jl status

ASCENT integration and the external HFSS full-wave campaign are tracked separately in
[`ASCENT.md`](ASCENT.md).

## Status

The original Stage 0–6 simplex-mesher roadmap is complete and remains the regression
baseline. The active roadmap now targets independent Gmsh 4.15.2 parity and is **not
complete**. Work is ordered by ASCENT meshing value before UI and post-processing.

| Stage | Capability | State |
|---|---|---|
| 0 | exact predicates, mesh types, validation, MSH/STL I/O | DONE |
| 1 | 2-D Delaunay, constrained Delaunay, quality refinement | DONE |
| 2 | size fields, graded curves, planar/cylinder/parametric surfaces | DONE |
| 3 | 3-D Delaunay, exact-coordinate kernel, volume fill, conforming recovery | DONE |
| 4 | uniform sizing, quality metrics, flips, smoothing, sliver reduction | DONE |
| 5 | healing diagnostics, native primitives, analytical CAD, imprints, mesh CSG | DONE |
| 6 | globally certified P2 tetrahedra and solver-consumable I/O | DONE |

### Active parity increment

| Track | State | Verified implementation increment |
|---|---|---|
| P1 | **IN PROGRESS** | Native scalar/anisotropic catalog, strict `.geo` field graph with injected model/view context, Gmsh-style 1-D policy, and field/entity-aware 2-D, surface, and 3-D refinement |
| P2 | **IN PROGRESS** | 125 fixed-node Gmsh types plus ten serializable cut/border/child/sub-element records, mixed blocks/entities/classification/periodic and embedded-curve metadata, structural validation/CRC, ASCII/binary MSH v2.2/v4.1 read/write with cumulative repeated-node/periodic sections and persistent MSH2 elementary ownership, and classified surface/explicit-shell/embedded-volume model-to-mixed projection |
| P3 | **IN PROGRESS** | Native analytical surfaces/imprints, classified ISO-10303-21 STEP/IGES box/sphere/cylinder/cone import, STEP/IGES NURBS curve and surface import with IGES export, expression-, numeric-list-, and tracked-tag-allocator-backed Point/Line/Surface/Surface Loop/Volume, Box/Cylinder/Sphere/Cone/Boolean/Translate/Dilate/90°-Rotate and straight-curve or planar-surface periodic `.geo` execution, mesh Boolean CSG, and finalized-mesh affine transforms |
| P4 | **IN PROGRESS** | Greedy and Edmonds-blossom surface recombination with optional full-quad, Point/Line-In-Surface embeddings, Point/Line/Surface-In-Volume recovery with nested constraints and holed planar sheets, explicit planar shell/cavity volumes, holed plane surfaces, uniform refinement, Progression/Bump/Beta curve laws and HWall variants, planar triangle/quad transfinite patches including recombined three-sided layouts, affine five-/six-face transfinite volumes, recombined hexahedra, prismatic 3-D layers with certified remaining-core fill/cavity walls, 2-D quad/fan layers, general-affine periodic node-pair certification/snapping, persistent native straight-curve relations for boundary or embedded curves with reusable masters and acyclic chains, synchronized planar periodic boundary surfaces on explicit volumes, expression/list-backed `.geo` periodic entities and transforms, and classified surface/volume projection with MSH2 cell ownership and supported MSH4 periodic/embedding metadata |
| P5–P6 | **IN PROGRESS** | Synchronized model/mesh API with detached cache and periodic-map ownership, non-destructive bounded CLI with periodic/embedded surfaces, embedded volumes, and periodic explicit-shell metadata output, validated headless GUI, owned scalar nodal views, synchronized in-process plugins, plus expression-, numeric-list-, and tracked-tag-allocator-backed geometry/entity lists, t1-square, t4-hole, classified Point/Line-In-Surface, nested and holed Surface-In-Volume, and explicit Surface Loop/Volume MSH lifecycles, native/projected single-/two-direction, embedded, reusable-master/chained, and expression/list-backed periodic checks, planar periodic explicit-volume boundaries, low-level translation/rotation-periodic checks, 2-D boundary-layer quad, API-box, OCC-cylinder/cone, IGES-128 bilinear, and BooleanDifference box Gmsh 4.15.2 differentials |

P1 does not claim 3-D multi-wall boundary-layer fans, the full Gmsh automatic-sizing
pipeline, high-order/custom-interpolation, or mixed-component
`PostView`, materially warped quadrangles, `PostView` tensor-to-metric evaluation, direct tensor or
metric-meshing parity, full `.geo`/CAD-model execution, or exact CAD distance. P2 does
not claim general mixed-element generation or recombination beyond P4's first-order
surface pairing, MINI basis-selector tags 138/139 as mesh records, integration of
mixed blocks into the simplex meshing kernels, curved high-order
Jacobian certification beyond P2 tetrahedra, ancillary/unknown-section preservation (binary readers reject
unsupported sections explicitly), non-8-byte binary data, internal
indices beyond `Int32`, or lossless multi-physical-group MSH v2.2 projection. Some
registered fixed tags and polygon-border type 69 require explicit Tessella-only output
because Gmsh 4.15.2 cannot consume them safely. MSH2 ASCII preserves variable records
and parent/domain links; binary MSH2 and MSH4 have explicitly narrower special-record
contracts. Pinned Gmsh 4.15.2 corrupts distinct parent links in its own binary MSH2
rewrite, and nonzero-physical special MSH4 requires compatible node/entity
classification metadata for a safe rewrite. P3 does not yet claim a general OpenCASCADE BREP kernel, NURBS CAD of
unclassified topology, transformations of arbitrary CAD entities, or full `.geo` execution.
Classified STEP/IGES solids that are axis-aligned blocks, spheres, right circular
cylinders, or right circular cones are imported and filled; STEP B-spline and
IGES 126/128 NURBS import as native curves/surfaces with IGES NURBS export;
other topology is an explicit blocker. Bounded
`.geo` execution covers Point/Line/Loop/Surface/Surface Loop/Volume,
Box/Cylinder/Sphere/Cone,
BooleanDifference/Union/Intersection of those solids, Translate of remaining
native solids, Dilate, and coordinate-axis π/2 rotations of native primitives,
Point/Line-In-Surface embeddings, and Point/Line/Surface-In-Volume recovery. Its
scanner and executor handle finite arithmetic constants, pure numeric functions,
prior scalar bindings, bounded numeric list assignment/indexing/selection/mutation,
explicit field/physical tags, and finite constant ranges in recognized numeric field
and entity lists. Executed geometry parameters, tags, and numeric entity memberships
use those same bounded semantics. It rejects control-flow
loops, macros, option reads, stateful functions, dynamic/general ranges,
logical/ternary evaluation, extrusions/fillets/symmetry, allocator reads after
topology-changing or untracked declarations, and mixed geometry-derived
physical-group right-hand-side evaluation.

P4's uniform-refinement slice applies the exact Gmsh 4.15.2 linear segment, triangle,
and tetrahedron child templates while sharing edge midpoints, compacting unused nodes,
and preserving parent tags. Its straight-curve slice covers normalized affine-line
Progression/Power, Bump, and Beta parameters plus all three HWall variants. Its surface
transfinite slice covers
already-discretized, count-matched, three- and four-sided planar chains using Gmsh's
specific triangular and average-chord Coons interpolation. Four-sided grids can also
be emitted as first-order Gmsh type-3 quadrangles with exact projected
corner-Jacobian certification. Affine eight-corner blocks use Gmsh's unrecombined
six-tetrahedron transfinite volume subdivision; canonical affine triangular prisms
use its legacy collapsed-grid five-face tetrahedral path. Surface recombination now
includes Edmonds blossom matching and a `full_quad` perfect-matching gate.
Planar polylines with an explicit oriented plane normal extrude to type-3
quadrangles along left-normals, with optional convex-corner fans and exact
projected corner-Jacobian checks. Closed manifold walls can use
`mesh_boundary_layer_filled` for certified prism shells, cavity walls, and a
conforming remaining-core tetrahedral fill. Explicit one-to-one translated or
general finite nonsingular affine node pairs can be certified and snapped exactly
without changing node numbering, connectivity, or tags. `MixedPeriodicLink`
retains pair maps and transforms in mixed meshes and through MSH2/MSH4 I/O.
`GeoModel` and the session API own straight-curve relations with one master per
slave, reusable masters, and acyclic master/slave chains. They synchronize boundary
or embedded subdivisions and expose their planar surface-node maps. They also own
affine-equivalent planar boundary-surface pairs of explicit volumes, synchronize
the slave facets, and expose certified tetrahedron-boundary node maps. `model_to_mixed`
projects planar triangle surfaces, including holes, Point/Line-In-Surface entities,
and shared periodic corners, into classified point/line/triangle blocks with
physical ownership, MSH2 elementary ownership, and supported MSH4 embedding and
periodic metadata. Its dimension-explicit volume form certifies the selected native
solid fill and emits classified point, line, surface, and tetrahedron blocks with
MSH2 elementary ownership and MSH4 entity classification. Explicit planar
surface-loop volumes require connected closed shells, support cavity loops, classify
every tetrahedron boundary face exactly once, and retain signed volume boundaries.
Periodic explicit shells retain their surface maps and induced boundary point/curve
forest through MSH2/MSH4.
Gmsh 4.15.2 does not serialize the Point/Line/Surface-In-Volume relation. Nested
Point/Line-In-Surface constraints are certified against each sheet's face complex;
embedded sheets may contain interior loops. MSH4 retains the nested curve relation.
The CLI uses these projections for periodic
or embedded `-2` output and classified `-3` output. P4 does not yet claim
non-affine CAD curve integration, FlexibleTransfinite, or size-map curve laws,
quasi-transfinite or holed transfinite patches,
general CAD parameterizations, curved/warped or
compact-TransfiniteTri volumes, volume/hybrid
recombination, selective or high-order refinement, coarsening,
3-D multi-wall boundary-layer fans, cyclic periodic-curve dependencies, curved
or non-boundary periodic surfaces, periodic volume entities, or allocator reads after
topology-changing or untracked declarations.
Expression/list-backed `Periodic Line`, `Periodic Curve`, and
`Periodic Surface` `Translate`, `Rotate`, and 12- or 16-entry `Affine` statements,
including bounded constant ranges and numeric list variables in entity sets, are in
scope for the bounded `.geo`
executor.

General OpenCASCADE/unclassified NURBS CAD, remaining algorithms/fields, broad
formats and API, GUI, and post-processing are unfinished parity tracks, not
project non-goals.

## Verification history (newest first)

Re-measured on 2026-08-28 with Julia 1.12.7 after adding bounded dynamic `.geo`
tag allocators:

- The scanner and executor now evaluate read-only `newp`, `newf`, and Gmsh's
  shared curve/loop/surface/volume allocator aliases while every preceding
  tag-producing statement remains in the tracked subset. Physical-group tags and
  the hidden Point/Curve/Surface topology of full Box, Cylinder, Sphere, and Cone
  primitives advance the same namespaces Gmsh 4.15.2 advances. Allocator reads
  after Boolean, deletion, or untracked topology remain explicit blockers.
- Bounds-checked focused sets passed under Julia 1.12.7 and 1.11.9: dynamic tags
  50/50, numeric lists 44/44, geometry expressions 49/49, model sets 62/62,
  118/118, and 84/84, periodic `.geo` sets 18/18 and 132/132, IO 422/422, CLI
  117/117, and size fields 6,976/6,976. Public-documentation and recursive
  ambiguity scans returned zero under both versions.
- The required Gmsh 4.15.2 differential confirmed 21 parser variables, repeated
  allocator reads, all aliases, Physical and Field namespaces, full and cone-tip
  OCC primitive allocation, entity tags, two periodic surface pairs, and analytic
  unit volume. Tessella produced 11 nodes and 16 tetrahedra; Gmsh produced 83 nodes
  and 188 tetrahedra. The maximum measured volume error was
  `8.881784197001252e-16`. The native mesh CRC was
  `2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808`,
  and the classified projection CRC was
  `99aeefc2e269090b518f5896f2388bd419c8fb44d06e4743579d50d61aaedf81`.
- The bounds-checked package gate passed 166,642/166,642 assertions in 27m36.6s.
  The bounds-checked aggregate validation gate, including the new allocator
  differential, exited successfully against Gmsh 4.15.2-git.
- The organization ratchet found 138 managed `.jl` files and no repository-root
  `.jl` files. The only top-level Julia entrypoints are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; all other Julia files are in
  domain or workflow subfolders.

Re-measured on 2026-08-27 with Julia 1.12.7 after adding bounded numeric `.geo`
list variables:

- The scanner and executor now cover zero-based indexing, cardinality, copies,
  concatenation, selection, whole-list append/removal, and indexed or selected
  mutation. Known lists expand in field options and entity-list positions,
  including embeddings, Physical groups, and periodic slave/master sets. Scalar
  writes preserve Gmsh's retained list payload and list/scalar mutation mode;
  geometry-derived unknown tails remain unavailable instead of collapsing to a
  one-item list. Bare list right-hand sides are accepted where Gmsh accepts them,
  while bare comma lists are rejected.
- The IO scanner passed 407/407 assertions under Julia 1.12.7 and 1.11.9. The new
  geometry set passed 44/44 under both versions; the neighboring expression sets
  passed 49/49, periodic sets passed 18/18 and 132/132, and the CLI passed 111/111
  under both versions. The related model sets passed 62/62, 118/118, and 84/84,
  and the size-field set passed 6,976/6,976 under Julia 1.12.7. Public-documentation
  and recursive ambiguity scans for Tessella, IO, and GeoExec returned zero under
  both Julia versions.
- The Gmsh 4.15.2 differential confirmed 17 parser lists, geometry and Physical
  entity reuse, field options, two periodic surface pairs, and analytic unit
  volume under both Julia versions. Tessella produced 11 nodes and 16 tetrahedra;
  Gmsh produced 83 nodes and 188 tetrahedra. The maximum measured volume error was
  `6.661338147750939e-16`. Tessella's native mesh CRC was
  `2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808`,
  and the classified projection CRC was
  `27417f652cf93e0d6aad41c2f1b6c65af3751dfb3cb3166432d2e798f25a6493`.
- The bounds-checked package gate passed 166,571/166,571 assertions in 13m37.2s
  (820.97s wall time). The bounds-checked aggregate validation gate, including the
  new numeric-list differential, exited successfully in 1,255.59s against Gmsh
  4.15.2-git.
- The organization ratchet found 136 managed `.jl` files and no repository-root
  `.jl` files. The only top-level Julia entrypoints are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; all other Julia files are in
  domain or workflow subfolders.

Re-measured on 2026-08-27 with Julia 1.12.7 after applying bounded constant
expressions to executable geometry statements:

- Numeric parameters, entity tags, and numeric entity lists now share the checked
  scalar/function/range evaluator across Point/Line/Loop/Surface/Surface Loop/
  Volume, native primitives, Booleans, transformations, embeddings, physical
  groups, and periodic statements. The new geometry-expression set passed 49/49
  assertions under Julia 1.12.7 and 1.11.9. The existing model sets passed 62/62,
  118/118, and 84/84; periodic `.geo` sets passed 18/18 and 132/132; and the CLI
  passed 106/106 under both versions. Public-documentation and recursive ambiguity
  scans for Tessella, GeoExec, and CLI returned zero under both versions.
- The Gmsh 4.15.2 differential confirmed expression-evaluated and
  truncation-toward-zero entity tags, oriented/entity range expansion, physical
  groups, an embedded volume point, and analytic unit volume. Tessella produced
  9 nodes and 12 tetrahedra; the measured Gmsh run produced 81 nodes and 184
  tetrahedra. The maximum volume error was `2.220446049250313e-16`. Tessella's
  native mesh CRC was
  `db4a080cdd8b4cdbd080d3ba42b798475d50a4590e67962c32edb8ac69205f24`,
  and the classified projection CRC was
  `89ee7d39873b202e264917e98fce2756b038d3e6f2f06d1bdbce9c46f7e628cd`.
- The bounds-checked package gate passed 166,498/166,498 assertions in 12m49.8s.
  The bounds-checked aggregate validation gate, including the new differential,
  exited successfully in 1,159.70s against Gmsh 4.15.2-git.
- The organization ratchet found 134 managed `.jl` files and no repository-root
  `.jl` files; the three designated top-level entrypoints remain under `src/`,
  `test/`, and `validation/`.

Re-measured on 2026-08-27 with Julia 1.12.7 after adding planar periodic
boundary surfaces for explicit volumes:

- The native model, bounded `.geo` executor, session API, CLI, and classified
  volume projection passed their focused periodic-surface checks. The new test
  sets passed 60/60 model, 18/18 `.geo`, 12/12 API, and 101/101 CLI assertions
  under Julia 1.12.7; the model, `.geo`, API, and CLI files also passed under
  Julia 1.11.9. Recursive public-documentation and ambiguity scans returned zero
  for the package, Model, GeoExec, API, and CLI under both Julia versions.
- The Gmsh 4.15.2 differential passed translated periodic boundary surfaces in
  two directions, including an embedded face point and ASCII/binary MSH2/MSH4
  round trips. Tessella produced 11 nodes and 16 tetrahedra; Gmsh produced 83
  nodes and 188 tetrahedra. The projected relation forest contained 15 relations,
  and the maximum affine-coordinate error was
  `1.5171197631502764e-13`. Tessella's native mesh CRC was
  `2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808`;
  its projected MSH4 CRC was
  `27417f652cf93e0d6aad41c2f1b6c65af3751dfb3cb3166432d2e798f25a6493`,
  and its MSH2 CRC was
  `9cc65eb95bbcca5508016ff7cc1340a6d1a7311d0482c2759444c16ce4120502`.
- The bounds-checked package gate passed 166,444/166,444 assertions in 18m00.8s.
  The bounds-checked aggregate validation gate, including the new differential,
  exited successfully in 1,126.45s against Gmsh 4.15.2-git.
- The organization ratchet found 132 managed `.jl` files and no repository-root
  `.jl` files; the three designated top-level entrypoints remain under `src/`,
  `test/`, and `validation/`.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding expression-aware
periodic `.geo` execution:

- `Periodic Line`/`Periodic Curve` entity lists and Translate/Rotate/Affine
  entries now evaluate prior scalar bindings, finite arithmetic, and pure
  numeric functions. Entity lists also expand bounded constant ranges and
  truncate positive tags toward zero as Gmsh does. The executor accepts the
  documented 12-entry Affine form by adding the homogeneous row and the full
  16-entry form required by the pinned Gmsh 4.15.2 runtime.
- The focused GeoExec, API, and CLI gates passed under Julia 1.12.7 and 1.11.9:
  132/132 GeoExec assertions, five API sets of 63/63, 26/26, 35/35, 27/27,
  and 49/49, and 94/94 CLI assertions. Recursive ambiguity and public
  documentation scans for the package, Model, API, GeoExec, and CLI returned zero
  under both versions.
- The Gmsh 4.15.2 differential passed expression-backed chained translation,
  rotation, and affine cases. Tessella produced 77 graph nodes and nine pairs;
  Gmsh produced 44 triangles and three pairs, with maximum graph error
  `8.184564212836642e-13`, maximum transform error
  `2.0590196214698153e-12`, and identical affine matrices. The expression graph
  retained native CRC
  `dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338`
  and projected CRC
  `3f98267cc70f9326ebe490c854cb59a9987c638e6aaabcba326d086bfb887ab1`;
  the affine and rotation fixtures produced native CRCs
  `3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5`
  and `f6ad616e56d52d7e10a598a4079db2de9b3d5f2a777f492f5a2366946d8ea990`.
- The bounds-checked package gate passed 166,347/166,347 assertions in 13m14.8s.
  Aggregate bounds-checked validation exited 0 in 18m30.5s against Gmsh
  4.15.2-git. All 131 source-managed Julia files remained categorized with
  only the three designated entry points at their top levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding acyclic
periodic-curve dependency graphs:

- Native boundary or embedded straight curves can reuse a master or form a
  master/slave chain. Graph validation keeps one master per slave and rejects a
  cycle before mutating the model. Subdivision parameters propagate across each
  connected graph, and affine snapping follows dependency order. Modeled points
  that coincide with reconstructed embedded-curve subdivisions retain one mesh
  vertex.
- The native model, model-to-MSH projection, bounded `.geo`, session API, and CLI
  periodic sets passed 118/118, 93/93, 93/93, 49/49, and 94/94 assertions under
  Julia 1.12.7 and 1.11.9. Recursive ambiguity and public-documentation scans for
  the package, Model, API, GeoExec, and CLI returned zero under both versions.
- The Gmsh 4.15.2 branch and chain differential passed for MSH2/MSH4 ASCII and
  binary. Tessella produced 77 nodes and nine pairs per curve; Gmsh produced 44
  triangles and three pairs per curve, with maximum affine error
  `8.184564212836642e-13`. Gmsh returned empty maps for the measured cyclic
  graph. The shared native mesh CRC is
  `dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338`;
  the branch and chain projected MSH4 CRCs are
  `6eb5020b186e9abdd8472312791bf375e1e9512066e7cc67b1fd2c1990a6d82f`
  and `3f98267cc70f9326ebe490c854cb59a9987c638e6aaabcba326d086bfb887ab1`.
- The bounds-checked package gate passed 166,308/166,308 assertions in 13m19.5s.
  Aggregate bounds-checked validation exited 0 in 17m46.6s against Gmsh
  4.15.2-git. `git diff --check` passed, and all 131 source-managed Julia files
  remained categorized with only the three designated entry points at their top
  levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after extending native periodicity
to straight curves embedded in planar surfaces:

- `mesh_model_surface` synchronizes boundary or embedded curve subdivisions,
  certifies the affine node map, and snaps slave nodes exactly. `model_to_mixed`
  retains the embedded curve entities, endpoint relations, and complete curve
  relation through MSH2/MSH4 ASCII and binary round trips.
- The embedded model, bounded `.geo`, API, and CLI test sets passed 51/51,
  38/38, 27/27, and 85/85 assertions under Julia 1.12.7 and 1.11.9. Recursive
  ambiguity and top-level/Model/API public-documentation scans returned zero
  under both versions.
- The Gmsh 4.15.2 differential measured zero affine coordinate error. Tessella
  produced 15 nodes, 20 triangles, and three compact curve-node pairs; Gmsh
  produced 16 triangles and two pairs. Gmsh reopened Tessella's MSH2/MSH4 ASCII
  and binary projections with the periodic curve and endpoint relations. The
  native mesh, projected MSH4, and mode-independent MSH2 CRCs are
  `9794a65ea5402683d0d50612522c2f71f7c98ec2a9f6b9e6b49a61e62cd85cf2`,
  `e32e8317842c099bc4a91cdd94d02d0f816884f0e091d7194bac56e95bbfeade`,
  and `d6da1835be0a570f81b99ccec03acd47bd46722ed69316007d3fc4fa020b2445`.
- The bounds-checked package gate passed 166,138/166,138 assertions in 13m24.9s.
  Aggregate bounds-checked validation exited 0 in 17m16.0s against Gmsh
  4.15.2-git. `git diff --check` passed, and all 130 source-managed Julia files
  remained categorized with only the three designated entry points at their top
  levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after recovering holed planar
Surface-In-Volume sheets:

- `mesh_model_volume` triangulates an embedded sheet's outer and hole loops with
  its nested point/curve constraints, recovers those triangles as tetrahedron
  faces, and rechecks every recovered curve and triangle against the final edge
  and face complexes. Classified projection retains the signed outer/hole
  boundaries and nested Curve-In-Surface relation.
- The classified native-volume, holed-sheet, and explicit-shell test sets passed
  72/72, 13/13, and 78/78 assertions under Julia 1.12.7 and 1.11.9. Recursive
  ambiguity and top-level/Model public-documentation scans returned zero under
  both versions.
- Tessella and Gmsh 4.15.2 measured the sheet area as 0.45 with zero triangle
  centroids inside the hole. Tessella produced 56 nodes and 203 tetrahedra; Gmsh
  produced 1,022 tetrahedra and 26 sheet triangles. Gmsh reopened Tessella's
  MSH2/MSH4 ASCII and binary projections with every classified point, curve,
  surface, and volume element present. The MSH4 and mode-independent MSH2 CRCs are
  `0e92af2702054065d564461a691f4035ab5bace358bd5f37821c9dcb5f54730d`
  and `250f6627ef3712e881a363b0e6d8999a6e77ea263503a0d813c6a0d42169a400`.
- The bounds-checked package gate passed 166,081/166,081 assertions in 13m25.5s.
  Aggregate bounds-checked validation exited 0 in 15m51.0s against Gmsh
  4.15.2-git. `git diff --check` passed, and all 129 source-managed Julia files
  remained categorized with only the three designated entry points at their top
  levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding explicit planar
surface-loop volumes:

- `add_surface_loop!` validates one connected closed shell and keeps Surface Loop
  tags in their own namespace. `add_volume!` records one exterior loop followed by
  surface-disjoint cavity loops; meshing and projection reject cavities outside the
  exterior or overlapping one another. The bounded `.geo` executor, session API,
  and CLI expose the same topology.
- Native meshing supports planar boundary surfaces, cavity shells, and nested
  Point/Line-In-Surface constraints. Classified projection certifies the selected
  fill and assigns every tetrahedron boundary face to exactly one modeled surface.
  MSH4 retains signed volume boundaries; MSH2 retains elementary cell ownership.
- The classified native-volume and explicit-shell test sets passed 72/72 and 78/78
  assertions under Julia 1.12.7 and 1.11.9. The CLI passed 75/75, and the API sets
  passed 63/63, 26/26, and 35/35 under both versions.
  Recursive ambiguity and public-documentation scans returned zero under both.
- Gmsh 4.15.2 and Tessella both measured the explicit unit volume as 1 within
  floating-point tolerance. Gmsh reopened Tessella's MSH2/MSH4 ASCII and binary
  outputs with all classified elements present. The projected cube and hollow-shell
  CRCs are `9bce88e319c67236317df64b876739a62f80982ed86eca028bd1e7bda022bcb6`
  and `83721952195b78f4d18b9e5ff862ef7629b9fbe6b642f33d6649d85e83d0c8b2`.
- The bounds-checked package gate passed 166,068/166,068 assertions in 13m12.0s.
  Aggregate bounds-checked validation exited 0 in 14m49.5s against Gmsh
  4.15.2-git. `git diff --check` passed, and all 128 source-managed Julia files
  remained categorized with only the three designated entry points at their top
  levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after composing nested sheet
constraints inside volume embeddings:

- `mesh_model_volume` now inserts nested Point-In-Surface nodes and recovers nested
  Line-In-Surface chains before recovering the enclosing Surface-In-Volume sheet.
  Its final certificate requires every nested point to be a sheet-face node and
  every nested curve edge to belong to the sheet face complex. A curve cannot be
  both a boundary and an embedded entity of the same sheet.
- The dimension-explicit `model_to_mixed` projection emits the nested point and
  curve cells with point-over-curve-over-surface-over-volume node ownership. MSH4
  retains the nested Curve-In-Surface relation; MSH2 retains elementary ownership.
  Off-sheet nested points and curves are rejected before mixed output is built.
- The classified-volume file passed 72/72 assertions and the CLI file passed 67/67
  under Julia 1.12.7 and Julia 1.11.9. Their nested projection CRCs are
  `e6a1a6de65b65987c543553d6456e4607b43fd3f3127294d926237888c9b5453`
  and `745bc23ab2aa7c0824006a94ef279514a1c6fa97d3cd95960943797b85c6336c`.
  The existing surface projection sets passed 93/93, 38/38, and 9/9 assertions;
  the Model sets passed 62/62, 61/61, and 84/84.
- The Gmsh 4.15.2 differential verified both nested relations in the live source
  model and all nested point/curve/surface/volume elements after Gmsh reopened
  Tessella's four MSH modes. MSH4 ASCII reconstructed Curve-In-Surface; MSH4 binary
  and both MSH2 modes had no relation, as required by the pinned format behavior.
  The nested MSH4 and MSH2 CRCs are
  `785bdb610878978e19cbcf3cfb2402417b9646d3ffc8f23042f922dc9c0b5930`
  and `7a0b70ac205dd985adfa6d2b0a789b791f7bdaab2ce4061c3b08e7eef1df99e4`.
- The bounds-checked package gate passed 165,956/165,956 assertions in 27m31.3s.
  Aggregate bounds-checked validation exited 0 in 20m17.1s against Gmsh
  4.15.2-git. Recursive ambiguity and public-documentation scans returned zero
  under Julia 1.12.7 and 1.11.9. `git diff --check` passed, and all 127
  source-managed Julia files remained categorized with only the three designated
  entry points at their top levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding classified native
volume-to-MSH projection:

- The dimension-explicit `model_to_mixed` path certifies that a validated
  tetrahedron mesh fills the selected native solid, then emits deterministic
  point, curve, embedded-surface, and volume blocks. Point ownership takes
  precedence over curves, curves over surfaces, and surfaces over the volume.
  Unrelated fills, overlapping curve edges or surface faces, nested constraints
  on an embedded sheet, periodic volume relations, and explicit modeled volume
  shells are blocked with diagnostics.
- Physical memberships and names survive the projection. MSH2 retains aligned
  elementary ownership for every emitted cell; MSH4 retains entity and node
  classification. Gmsh's `src/geo/GModelIO_MSH4.cpp` at commit
  `657c8e915f60405e6cad0c8ec7faf812bfff1a60` and the four-mode reopen
  differential verified that Gmsh 4.15.2 has no serialized
  Point/Line/Surface-In-Volume relation, while the classified entities and
  elements remain usable. The bounded CLI selects this path for embedded `-3`
  input.
- The classified-volume fixture passed 66/66 assertions and the CLI file passed
  66/66 under both Julia 1.12.7 and Julia 1.11.9. Their projection CRCs are
  `d12446ff4f9c24254d02ae3938fc51b1010063399ab1a4d42b8947bd6637c1f8`
  and `2fc633ca160054b1c8f86b9981febc79c83ff248afe409aad7fbb8e1e3f27ef4`.
  The Mesh3D file passed 534/534 assertions after exposing oriented covering
  faces for classification.
- The Surface-In-Volume differential passed against Gmsh 4.15.2 in MSH2/MSH4
  ASCII and binary modes. Its MSH4 projection CRC is
  `2ccb9e42be0322ab810c601a99527b241c7e02d7269aa8c51b306f61347b6027`;
  the mode-independent MSH2 CRC is
  `8141bdc49a658a1770944e32527e2971b598f56ea2b0a4bae27e128bee0da640`.
- The bounds-checked package gate passed 165,949/165,949 assertions in 22m33.2s.
  Aggregate bounds-checked validation exited 0 in 25m43.2s against Gmsh
  4.15.2-git. Recursive ambiguity and public-documentation scans returned zero
  under Julia 1.12.7 and 1.11.9. `git diff --check` passed, and the organization
  ratchet found all 127 source-managed Julia files categorized, with only
  `src/Tessella.jl`, `test/runtests.jl`, and `validation/run_all.jl` at their
  respective top levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding classified
Point/Line-In-Surface output:

- `MixedEntity` owns validated surface-embedded curve tags, and `mixed_crc`
  includes them without changing historical digests for meshes that have no
  embedding relation. ASCII and binary MSH4 readers and writers preserve Gmsh's
  encoded Curve-In-Surface record. MSH2 readers now retain declared elementary
  ownership independently of periodic links; conversion preserves positive entity
  tags when every entity has one legacy physical membership and blocks ambiguous
  or lossy layouts atomically.
- `model_to_mixed` emits classified point, boundary/embedded-line, and triangle
  cells for planar surfaces with Point/Line-In-Surface constraints. Embedded points
  take node-classification precedence when they split an embedded curve. Physical
  memberships, names, signed boundaries, curve ownership, and supported periodic
  links survive the projection. The bounded CLI selects this path for embedded or
  periodic `-2` input. Periodic embedded curves are an explicit blocker; classified
  volume-embedding projection remains outside this increment.
- The focused Elements embedding/conversion test sets passed 52/52 assertions,
  the model projection file passed 140/140, and the CLI file passed 56/56 under
  both Julia 1.12.7 and Julia 1.11.9. The combined embedded fixture and CLI CRCs
  are `e762c7c566f1e5768ad1e2849302815dbfd9d19a14c1b3840abcefa4aedcaf43`
  and `6025846e0f58581418401081f092630d2e49999a26e608bf377d1bae4c51dc4b`.
- The Gmsh 4.15.2 differentials passed for classified embedded points and curves.
  Their projection CRCs are
  `222619f8e92298ab72ece09cae6dd9f8300781c9d8a7dc587fc5b3de50219fc3`
  and `0655fe3edb4344be584d2fe12b8d57637f65090524542d2bc1b515e148d55ea5`.
  Gmsh reopened the ASCII curve relation and both point-element variants. Its
  binary reopen retained the curve elements but reported no embedding relation;
  Tessella's ASCII and binary round trips retained the relation and identical CRC.
- The bounds-checked package gate passed 165,873/165,873 assertions in 13m11.8s.
  Aggregate bounds-checked validation exited 0 in 12m57.0s against Gmsh
  4.15.2-git. Recursive ambiguity and public-documentation scans returned zero
  under Julia 1.12.7 and 1.11.9. `git diff --check` passed, and the organization
  ratchet found all 126 Julia files categorized, with only `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl` at their respective top levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding classified native
surface-to-MSH projection:

- `model_to_mixed` accepts validated, unembedded planar triangle meshes whose
  boundary chains represent one native surface. It emits point, boundary-line,
  and triangle blocks; MSH2 elementary tags; MSH4 point/curve/surface ownership;
  physical memberships and names; and every stored periodic curve link. Hole
  loops receive the signed reverse traversal required by Gmsh. Multiple physical
  memberships remain lossless in MSH4; MSH2 uses the lowest tag as its single
  legacy membership.
- Periodic curve nodes must already be exactly snapped. Endpoint metadata uses a
  deterministic one-master-per-slave spanning forest when independent relations
  share corners. Invalid, tagged, volumetric, nonplanar, unsynchronized, damaged,
  or embedded inputs raise explicit diagnostics. The bounded CLI writes this
  classified MSH4 path for periodic `-2` input and blocks periodic `-3` output or
  relations outside the selected surface without changing an existing destination.
- The focused projection suite passed 102/102 assertions and the CLI suite passed
  47/47 under Julia 1.12.7 and Julia 1.11.9. The physical-membership fixture CRCs
  are `12a1eb50575a3af08273b1a0fdefca49d7e4b01b4898573e6346b6b61b4978c3`
  for MSH2 and
  `d9aa0af0ed218f321adea7b7276312583ee31f4747e3771ca410b87be3b628b7`
  for MSH4. The holed projection CRC is
  `654d305c58cfe9db2f87ac1424a31912863a21f75196beabc3f05c37d8f6e73f`.
- The Gmsh 4.15.2 differential reopened Tessella's ASCII and binary MSH2/MSH4
  projections and recovered both curve and endpoint relations. The single-direction
  projected coordinate error was `0.0`; its MSH2/MSH4 CRCs are
  `506ae0fac8562df49231df71f3b12d7259ba44b3fb5618a064a15f97698951a0`
  and `cf03be1a36427f1ef0fbc4e852996bd65d2630b5ac384fa0267dd14e46ea6280`.
  A two-direction fixture recovered three point links and two five-pair curve links
  in all four file modes, with MSH2/MSH4 CRCs
  `bac00f74b86af8d1a6b70de445cdb17a16a9513f0fc4a542bd995d9120923a58`
  and `d5fd8bd6ef46c78772792f0cee0c7b19cdd747f1c2932c2a8992760f19e69b20`.
- The bounds-checked package gate passed 165,771/165,771 assertions in 13m11.4s.
  Aggregate bounds-checked validation exited 0 in 12m19.8s against Gmsh
  4.15.2-git. Recursive ambiguity and public-documentation scans returned zero,
  `git diff --check` passed, and the repository organization ratchet found all 126
  Julia files categorized, with only the three designated entry points at their
  top levels.

Re-measured on 2026-08-26 with Julia 1.12.7 after enabling bounded native `.geo`
periodic-curve execution:

- `execute_geo` accepts legacy `Periodic Line` and current `Periodic Curve`
  statements between straight curves with finite literal `Translate` triples,
  `Rotate` axis/center/angle data, or Gmsh's 12-entry `.geo` `Affine` form.
  Variable entity tags, numeric expressions, curved entities, and periodic
  surfaces or volumes remain explicit blockers.
- Parsed constraints persist in `GeoModel`, mesh through the native planar
  surface path, and survive `API.open_geo!` with detached cached-mesh ownership.
  The translated and non-origin rotated `.geo` mesh SHA-256 values are
  `3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5`
  and `f6ad616e56d52d7e10a598a4079db2de9b3d5f2a777f492f5a2366946d8ea990`.
- The focused `.geo` periodic suite passed 31/31 assertions and the API file
  passed 98/98 under both Julia 1.12.7 and Julia 1.11.9. Invalid dimensions,
  transform arities, nonfinite or nonliteral values, zero rotation axes, and
  malformed entity lists all block.
- The Gmsh 4.15.2 differential now builds Tessella's five native pairs by
  executing a checked-in `.geo` fixture. Their maximum coordinate difference
  from Gmsh's five pairs is `2.0594637106796654e-12`; the existing translated,
  rotated, MSH2, and MSH4 CRCs remain unchanged.
- The bounds-checked package gate passed 165,650/165,650 assertions in 13m12.2s.
  Aggregate bounds-checked validation exited 0 in 12m39.7s against Gmsh
  4.15.2-git. Recursive ambiguity and public-documentation scans returned zero,
  `git diff --check` passed, and the repository organization ratchet kept every
  non-entry-point `.jl` file in a categorized subfolder.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding persistent native
straight-curve periodic constraints:

- `GeoModel` owns atomic, immutable row-major affine curve relations and their
  endpoint orientation. Planar surface meshing uses bounded remeshing to unify
  master/slave boundary parameters, certifies each boundary edge chain, snaps
  every slave node exactly, and rejects incompatible constraints at shared
  corners. Nonperiodic model CRCs remain unchanged.
- The direct and session APIs return detached `Int32` node maps. Successful
  constraint updates invalidate the API mesh cache; rejected updates leave the
  cache and model unchanged. The translated, rotated, and two-direction model
  SHA-256 values are
  `6ea713b4493eeb5b31e7c70ea5312290ef698424fd787bb04f1df4ef55f894cf`,
  `f6ad616e56d52d7e10a598a4079db2de9b3d5f2a777f492f5a2366946d8ea990`,
  and `b82c9f0f4e235e90a754f2ec50b3a373ef0a2d514a79194d9b922873e35f8dd1`.
- The focused Model file passed 207/207 assertions and the API file passed
  93/93 under both Julia 1.12.7 and Julia 1.11.9. Their periodic test sets
  passed 61/61 and 35/35, including invalid meshes, unsynchronized maps,
  nonconvergence, partial-surface relations, multiple directions, and
  inconsistent shared-corner transforms.
- The Gmsh 4.15.2 differential matched five native model pairs to the same five
  Gmsh curve pairs with maximum coordinate difference
  `2.0594637106796654e-12`; the native mesh SHA-256 is
  `3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5`.
- The bounds-checked package gate passed 165,614/165,614 assertions in
  13m56.2s. Aggregate bounds-checked validation exited 0 in 13m26.6s against
  Gmsh 4.15.2-git. Recursive ambiguity and public-documentation scans returned
  zero, `git diff --check` passed, and the repository organization ratchet kept
  every non-entry-point `.jl` file in a categorized subfolder.

Re-measured on 2026-08-26 with Julia 1.12.7 after adding persistent standard
MSH2 periodic sections:

- `MixedMesh.elementary_entities` owns aligned per-cell MSH2 elementary tags.
  ASCII and binary MSH2 readers retain those tags when a periodic section uses
  them; both MSH2 writers preserve the exact classification and emit the
  standard ASCII `$Periodic` payload, including optional 16-entry `Affine`
  transforms.
- The MSH2 reader accepts the format's whitespace-separated grammar, merges
  disjoint repeated sections under cumulative link/pair limits, and rejects
  missing entities, unknown nodes, duplicate slave relations, malformed
  transforms, and metadata-loss conversions. MSH4 behavior and historical CRCs
  remain unchanged; an equal MSH2/MSH4 classification copy is CRC-redundant.
- The fixed MSH2 periodic fixture SHA-256 is
  `a914cf9dd0fb8f7f5c01cea90bbb9457009d579f979ef66fdc1bf924420d830f`.
  The focused Elements suite passed 3,078/3,078 assertions under Julia 1.12.7
  and Julia 1.11.9; its periodic test set passed 124/124.
- The Gmsh 4.15.2 lifecycle differential generated and reopened MSH2 and MSH4
  files in both modes, then recovered the same three entity links and five
  curve-node pairs after Tessella and Gmsh rewrites. The mixed-mesh SHA-256
  values from Gmsh-generated files are
  `cf3f029b790af950d3d8c1e307c99970665bd12de83259afc50448bdf5f4cc6f`
  for MSH2 and
  `461ca91e6359638ebf2be97537660ff0ce760cebf4a06dffd770debe43b62a16`
  for MSH4. **VERIFIED (`getPeriodicNodes` and `getAttributeNames`):** pinned
  Gmsh exposes an MSH2 relation as both live metadata and a raw `Periodic` model
  attribute; the probe removes that duplicate attribute before asking Gmsh to
  serialize the live relation.
- The bounds-checked package gate passed 165,518/165,518 assertions in
  19m43.2s. Aggregate bounds-checked validation exited 0 in 15m59.1s against
  Gmsh 4.15.2-git. Recursive ambiguity and public-documentation scans returned
  zero, `git diff --check` passed, and the repository organization ratchet kept
  every non-entry-point `.jl` file in a categorized subfolder.

Re-measured on 2026-08-25 with Julia 1.12.7 after adding persistent standard
MSH4 periodic sections:

- `MixedPeriodicLink` owns 0-D/1-D/2-D slave/master entity relations, compact
  node pairs, and optional finite nonsingular row-major affine transforms.
  `MixedMesh` validates and CRC-hashes the metadata independently of link or
  pair order. The fixed fixture SHA-256 is
  `9b6f017f0bc019b6d96a12496e67d05046387a1e9d1707939d220a669788f348`.
- ASCII and native-endian binary MSH4 readers and writers preserve the
  relations and sparse external tags. Opposite-endian binary input, repeated
  disjoint sections, cumulative link/pair limits, atomic output blockers, and
  malformed records are covered. Geometry-model persistence and lossless MSH2
  periodic metadata remain outside this increment.
- The focused Elements suite passed 3,022/3,022 assertions under Julia 1.12.7
  and Julia 1.11.9; its periodic test set passed 68/68. The Gmsh lifecycle
  differential generated ASCII and binary periodic files, loaded three links
  in Tessella, rewrote both modes, and recovered the same five curve-node pairs
  in Gmsh. Its MSH4 SHA-256 is
  `461ca91e6359638ebf2be97537660ff0ce760cebf4a06dffd770debe43b62a16`.
- The bounds-checked package gate passed 165,462/165,462 assertions in
  13m53.7s. Aggregate bounds-checked validation exited 0 in 16m05.5s against
  Gmsh 4.15.2-git. Recursive ambiguity and public-documentation scans returned
  zero, `git diff --check` passed, and the repository organization ratchet
  retained only the three designated Julia entry files at their top levels.

Re-measured on 2026-08-25 with Julia 1.12.7 after extending periodic node-pair
certification to general affine transformations:

- `periodic_identify_affine` accepts either a finite 4×4 matrix or Gmsh's
  16-entry row-major representation, requires an exact affine homogeneous row
  and an exactly nonsingular 3×3 linear part, and verifies all disjoint
  one-to-one pairs before snapping a copied mesh. It uses the shared
  exact-dyadic cancellation fallback, blocks unrepresentable outputs, and
  independently validates the completed mesh. Translation behavior and its
  historical CRC remain unchanged.
- The focused periodic suite passed 559/559 assertions under Julia 1.12.7 and
  Julia 1.11.9. It includes 128 fixed-seed exact-rational affine oracles spanning
  rotations, reflections, shear, scaling, and translation. The fixed rotational
  chain SHA-256 is
  `ed4b81783f68a3bb092ac6fa2156196efe7f91fec6c7ceeef1aee154fb85261a`.
- The direct Gmsh 4.15.2 differential certified five translated and five +90°
  rotated curve-node pairs. Maximum pre-snap errors were
  `2.0594637106796654e-12` and `0.0`; both canonical outputs have SHA-256
  `baa96c7ebc0265667209f1940c77d5bdeed5ecb8a12f765d02df9d1945373648`.
- The bounds-checked package gate passed 165,385/165,385 assertions in 12m42.5s.
  Aggregate bounds-checked validation exited 0 in 11m0.8s against Gmsh
  4.15.2-git. Recursive ambiguity and top-level public-documentation scans both
  returned zero.
- The repository organization ratchet passed: no tracked root-level Julia file
  exists; only `src/Tessella.jl`, `test/runtests.jl`, and
  `validation/run_all.jl` occupy their respective top levels, while every other
  tracked `.jl` file remains in a categorized subfolder.

Re-measured on 2026-08-25 with Julia 1.12.7 after adding cumulative repeated
pre-element `$Nodes` sections:

- `read_mixed_msh` now merges disjoint node sections in ASCII and binary MSH
  v2.2/v4.1 while enforcing cumulative node/block limits and global external-tag
  uniqueness. V4 entity classification, sparse external tags, element
  connectivity, and canonical one-section rewrites are preserved.
- Fixed three-node/two-line CRCs are
  `d93f18ff2f3415913e2abd4f31eeb119896dda57755bfd230e616ead6a57c84e`
  for MSH2 and
  `24f2fefdad2a699abf29e9007ebc5c78ff7f80cfa59bbb5c15bdc8a815d3406e`
  for MSH4, identical between ASCII and binary inputs.
- **VERIFIED (`gmsh <file> -check -parse_and_exit -v 5`):** Gmsh 4.15.2
  reports tag 10 from the earlier section as unknown for each of the four raw
  fixtures. Tessella's merged canonical rewrites contain one `$Nodes` section
  and all four are accepted by the same command without an error. Empty repeated
  sections are accepted; duplicates across sections and node sections after
  elements are rejected explicitly.
- The bounds-checked focused Elements suite passed 2,945/2,945 assertions under
  Julia 1.12.7 and Julia 1.11.9. The full package gate passed
  164,844/164,844 assertions in 12m39.6s, and aggregate bounds-checked validation
  exited 0 against Gmsh 4.15.2-git. The tracked Julia layout remains fully
  categorized, with only the three entry-point files at their respective top
  levels.

Re-measured on 2026-08-25 with Julia 1.12.7 after extending 2-D boundary-layer
strips to arbitrary oriented planes:

- `mesh_boundary_layer_2d` now accepts a finite nonzero `plane_normal`, with
  overflow-safe normalization, a deterministic orthonormal frame, scale-aware
  input/output planarity checks, and preserved positive-`z` default behavior.
  Reversing the normal reverses the selected side. The historical flat-strip CRC
  remains `bb9e1fb9a0f0e56de42287ff3f85dd93ea3c7115cc9e51d05a6845d97ce8122b`;
  the exact vertical-plane CRC is
  `49e70bbffb8fd2121a526c35a5ca51ab87ea19790bd715d52a147a21535b3e19`.
- Every emitted triangle and all four corners of every quadrangle are certified
  with exact predicates in the projected plane. The suite rejects nonplanar
  inputs, unrepresentable offsets, and a sharp-turn fixture whose positive total
  shoelace area hides a reversed corner Jacobian.
- A tilted corner-fan oracle matched every rigidly transformed coordinate within
  `2.220446049250313e-16`, and a seeded audit repeated the topology/coordinate
  comparison for 128 random plane orientations and translations. The aggregate
  Gmsh 4.15.2 child retained `tessella_area=0.15960000000000005`,
  `gmsh_quads=10`, and `gmsh_tris=24`, while its added tilted straight-strip
  oracle had zero measured coordinate error.
- The bounds-checked focused suite passed 92/92 assertions under Julia 1.12.7
  and Julia 1.11.9. The full package gate passed 164,769/164,769 assertions in
  12m42.8s, and aggregate bounds-checked validation exited 0 against Gmsh
  4.15.2-git. All tracked Julia files remain organized in categorized subfolders,
  with only the package/test/validation entry points at their respective top
  levels.

Re-measured on 2026-08-25 with Julia 1.12.7 after adding Gmsh-compatible
recombined three-sided transfinite patches:

- `mesh_transfinite_triangle_patch` emits one first-order triangle per logical
  row and `d(d-1)/2` first-order quadrangles for `d` divisions while preserving
  all boundary segments and physical tags. `:left`, `:right`, and the shared
  alternate layout reproduce Gmsh 4.15.2's ordered connectivity. Fixed
  four-division mixed CRCs are
  `09d6619152fe4f42d604c9a95e3805843825a08615d92778ae4e551d85fa2ce3`,
  `b401b6bc71cac6dc3f7f44b20a4128ddf4bb439acec941253b7598af19800205`,
  and `5d200b76825ed699e99b125c28cef49158d30bc56a9e090d8469854cc84dabfa`.
- **VERIFIED (Gmsh 4.15.2 API oracle):** coverage of all four arrangement names,
  five resolutions, and planar/tilted geometries matched Gmsh node placement and
  triangle/quadrangle connectivity across 205 nodes per unrecombined/recombined
  track. Maximum node error was `2.808666774861361e-15`. ASCII/binary MSH
  v2.2/v4.1 round trips preserved coordinates and canonical mixed cell/tag sets.
- Every output is certified against the unrecombined atomic lattice, exact
  projected triangle/quadrangle orientations, and mixed edge incidence. A
  base-valid fixture whose `Left` pairing creates a concave final quadrangle is
  rejected explicitly. A separate seeded audit validated 2,000 affine patches
  across all arrangements with division counts sampled from `1:20`.
- The bounds-checked focused suite passed 241/241 assertions under Julia 1.12.7
  and Julia 1.11.9. The final full package gate passed 164,744/164,744
  assertions in 13m20.3s, and aggregate validation exited 0 against Gmsh
  4.15.2-git. The organized Julia layout remains enforced: only
  `src/Tessella.jl`, `test/runtests.jl`, and `validation/run_all.jl` occupy their
  respective top levels; all other tracked `.jl` files remain in categorized
  subfolders.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening simplex/mixed MSH
I/O, STL ingestion, and the straight-curve/recombined-quad public boundaries:

- ASCII MSH v2.2/v4.1 input is resource-bounded and validated before return.
  MSH4 accepts Gmsh-compatible implicit entities, empty blocks, and repeated
  entity/element sections; repeated physical-name records must agree. A fixed
  repeated-section mesh has SHA-256
  `4c4f930b531f67093077abbde261fe1e48b4e163ce4885b89cd2bca83581bf26`.
  Writers validate mutable connectivity before indexing and replace targets
  atomically. Standard names preserve literal backslashes and tabs, while
  Gmsh-unsafe quotes/line breaks and its measured 128-byte name overflow are
  explicit blockers.
- **VERIFIED (Gmsh 4.15.2 rewrite oracle):** simplex and mixed MSH4 writers now
  classify nodes on an element-owning entity. Gmsh rewrites ASCII and binary
  sources without the former doubled node count, and Tessella rereads the
  rewritten files with identical connectivity/physical-tag CRCs. The same
  oracle preserved literal backslashes and tabs. The aggregate enclosure probe
  now records Gmsh's partial file as structurally invalid because it contains
  duplicate surface cells, while independently confirming zero volume cells.
- Exact STL welding now handles opposite finite Float64 extrema at zero or
  positive tolerance, checks facet/node/file ceilings before corresponding
  growth, and converts malformed text failures to controlled diagnostics. A
  seeded 5,000-case byte-mutation audit of each MSH and STL reader produced only
  valid results or `ArgumentError` (MSH: 65/4,935; STL: 48/4,952); the `.geo`
  scanner likewise returned 850 valid parses and 4,150 controlled blockers.
  The small end-oriented HWall regression has parameter SHA-256
  `757122bf5807435f676f31ffe748f46bc1577a80d13735f9d8f862180bb341b8`;
  its subtraction tolerance is scaled by the represented endpoint.
- Bounds-checked focused suites passed under Julia 1.12.7 and Julia 1.11.9:
  IO 383/383, Elements 2,870/2,870, straight curves 514/514, and recombined
  quads 130/130. The final full package gate passed 164,607/164,607 assertions
  in 13m15.6s, and aggregate validation exited 0 against Gmsh 4.15.2-git.
  Recursive ambiguity detection and the public documentation scan returned
  zero. The organization gate still finds only `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl` at their respective top
  levels; all other Julia files remain in categorized subfolders.

Re-measured on 2026-08-25 with Julia 1.12.7 after bounding exact 3-D
tetrahedralization, PLC recovery, and conforming size refinement:

- **VERIFIED (independent topology oracle):** 200/200 seeded,
  non-cospherical random clouds produced the same coordinate-valued tetrahedron
  sets in the exact-rational kernel and the separately implemented ghost-vertex
  Float64 kernel. The 12 fixed cases retained in
  `test/meshing/mesh3d_test.jl` pass under Julia 1.12.7 and Julia 1.11.9. A
  `3×3×3` exact grid has 48 tetrahedra and connectivity SHA-256
  `998c29691aaadcbeed8f47d61c3729f8836e68f33bbc1437448586103003e623`.
- `delaunay3d_exact` and `is_delaunay_exact` now losslessly accept supported
  generic rational/integer vectors and integer connectivity, reject malformed
  or Boolean inputs explicitly, and bound points, returned/work tetrahedra, and
  cumulative exact-predicate evaluations. An unreadable-vector fixture verifies
  the point ceiling before point access. A seeded 1,000-call malformed-input
  audit returned 1,000 `ArgumentError`s and no unexpected result or exception.
- Exact boundary/partition recovery now bounds input facets, returned/work
  tetrahedra, predicate work, and cumulative edge/region certification work.
  `refine_to_size` separately bounds nodes, live and accumulated tetrahedra,
  segments, and triangles before each growth operation. Successful and failing
  recovery/refinement probes left their input meshes unchanged. Fixed CRCs are
  `f2451e6cb9e424bc520d8b9723fe1d307537f99fd9d8404e98490d2fae4b8bab`
  for the recovered box,
  `0d14f7a477b4d7222dce20b61cb1f7663c1c131dfbf1e945b48e5999e93569c1`
  for the one-region partition,
  `83fbfd93eeff0b9dde4e2f661fc589a87c0294a103ce62f1f75702a3efd4f7e1`
  for its sized CDT, and
  `c3d7c10942ce6348de44d5bb6396a7328c5d035221f1f40fbd34d7064278d8a0`
  for tagged single-tetrahedron refinement.
- The focused Mesh3D suite passed 534/534 assertions under both Julia 1.12.7
  and Julia 1.11.9. The bounds-checked package gate passed
  164,502/164,502 assertions in 18m59.6s; aggregate bounds-checked validation
  exited 0 against Gmsh 4.15.2. Recursive ambiguity detection and the public
  documentation scan both returned zero.
- The enforced Julia-file organization remains satisfied: the only top-level
  Julia files under `src`, `test`, and `validation` are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; all implementation and
  supporting `.jl` files are organized in categorized subfolders.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening analytical CAD,
constructive primitive surfaces, and healing diagnostics:

- **VERIFIED (independent projection oracles):** projecting
  `(floatmax(Float64),floatmax(Float64),0)` onto the unit sphere, z-cylinder,
  and z-disk now returns the diagonal point
  `(0.7071067811865475,0.7071067811865475,0)`, rather than a collapsed center
  or arbitrary radial. A z-cylinder query `(1,0,1e308)` at radius `0.7`
  formerly returned `(0,0.7,1e308)`; exact-dyadic ambiguity recovery now
  returns the nearest `(0.7,0,1e308)`. Opposite `±floatmax` endpoints whose
  direct difference overflows now retain finite plane, disk, sphere, cylinder,
  and circle-imprint projections. Fast finite projections remain
  allocation-free; four separately specialized 100,000-call measurements
  allocated zero bytes.
- CAD vector normalization is exponent-scaled. Projection uses compensated
  arithmetic with an exact `Rational{BigInt}` fallback and a 2,304-bit final
  rounding path, then certifies the represented target surface. Imprints use
  exact parallel/orthogonality decisions, certify both surfaces, and enforce a
  caller-visible point ceiling before allocation. A 3,000-case seeded,
  exponent-varied audit returned 1,149 certified projections and 1,851 precise
  representability blockers, with no unexpected exception and worst
  independent 512-bit residual ratio `1.4916992607748145e-13`. A separate
  5,000-case 1,024-bit nearest-point oracle returned 3,136 projections and
  1,864 blockers with no mismatch; its worst relative coordinate error was
  `2.387918906429682e-14`.
- Primitive builders now reject inappropriate `Bool` and nonnumeric values,
  normalize finite axes whose ordinary norm overflows, preflight cylinder,
  sphere, and cone resource counts before reading point storage, and reject
  radii/levels that collapse at a remote origin or underflow during cone
  interpolation. Existing sphere and cone connectivity CRCs remain
  `2c7bf12222ab5796df858b3ef015349be3fc8acecff5d445963f41215377bd54`
  and `8afc4d9f3bf9740313d9ce099302acb32335a9de392d8940d60cdb42b9465115`.
  A seeded 1,000-case exponent audit returned 825 valid meshable surfaces and
  175 explicit blockers with no unexpected exception; a separate 5,000-case
  malformed audit produced 5,000 `ArgumentError`s.
- **VERIFIED (healing complexity regression):** 5,000 isolated points sharing
  one x-coordinate at `tol=1e-320` took 3.554217458 seconds in the former
  x-only fallback. The exact BigInt spatial grid returns the same zero-pair
  result in 0.243028958 seconds and exactly preserves strict subnormal
  distance comparisons. Mutable connectivity and tag structure is checked
  before indexing, non-finite coordinates remain diagnosable, and a mesh that
  already contains tetrahedra cannot pass the surface meshability gate.
- The focused Heal/Geometry/CAD suites passed 45/45, 88/88, and 8,205/8,205
  assertions under both Julia 1.12.7 and Julia 1.11.9. Related MeshTypes,
  NURBS, BRep, model, and HighOrder suites passed 2,843 assertions. The
  bounds-checked package gate passed 164,449/164,449 assertions in 14m29.0s;
  aggregate bounds-checked validation exited 0 against Gmsh 4.15.2. Recursive
  ambiguity detection and the public documentation scan both returned zero.
- The Julia-file organization policy remains satisfied: the only top-level
  Julia files under `src`, `test`, and `validation` are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; every implementation and
  supporting Julia file remains in a categorized subfolder.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening and organizing the
top-level planar, extrusion, sizing, and volume-meshing pipelines:

- **VERIFIED (endpoint oracle):** for
  `z0=4.097470032826895e-162`, `z1=2.5149445970871698e185`, and
  `hmax=4.710819546778298e182`, the former expression for the last of 755
  layers evaluated to `2.51494459708717e185`, not the requested upper endpoint.
  Exact-dyadic layer counting plus endpoint-pinned convex interpolation now
  produces 756 distinct represented levels, and every node on the final level
  equals `z1` exactly. The layer and edge certificates are isolated in
  `src/meshing/PipelineSupport.jl`.
- Extrusion now chooses a conservative transverse step whose represented
  diagonal is no greater than `hmax`, audits the refined planar edges and every
  output tetrahedron edge, rejects non-finite or zero represented volumes, and
  checks exact node/tetrahedron counts against caller resource ceilings before
  dense 3-D allocation. An unreadable-vector fixture confirmed that the minimum
  `max_nodes`/`max_tets` preflight runs before point access.
- Integer planar coordinates now produce the same fixed topology as Float64
  coordinates. The planar CRC is
  `850fe31fb8b9c7946d716633cfabdfaf13850456a1b53474d21edfcfa9f194f4`;
  the fixed 10-node/12-tetrahedron extrusion CRC is
  `c7783021725d2dfd0b60b83536b5489f556b564af35fe66ef487e0bce15d9e3e`.
  Boolean coordinates, bounds, counts, seeds, tags, and pipeline controls now
  receive explicit `ArgumentError` diagnostics instead of dispatch failures or
  numeric coercion. Input-driven non-finite derived simplex measures are also
  `ArgumentError`s, and input `ArgumentError`s retain that category through the
  top-level fill wrapper.
- A seeded 5,000-case malformed-input audit returned 5,000 bounded
  `ArgumentError`s and no other result or exception type. A separate 100-case
  seeded rectangle/extrusion audit returned 100 valid meshes, pinned every top
  endpoint, and found no edge-bound failure; its worst measured
  `maxedge/hmax` was `0.9994607115469736`.
- The focused pipeline suite passed 82/82 assertions under Julia 1.12.7 and
  Julia 1.11.9. The focused Mesh3D suite passed 481/481 assertions; the affine
  volume and prism suites retained 85/85 and 133/133 assertions and their
  allocation ratchets.
- The bounds-checked package gate passed 164,387/164,387 assertions in
  13m29.3s. Aggregate bounds-checked validation exited 0 in approximately
  11m57s against Gmsh 4.15.2. Recursive package ambiguity detection and the
  public documentation scan both returned zero.
- The Julia-file organization is executable policy: the only top-level Julia
  files under `src`, `test`, and `validation` are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`. Pipeline support now lives in
  the `src/meshing` subfolder, and every other implementation/supporting Julia
  file remains categorized in an enforced domain subfolder.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening affine
transfinite-volume interpolation and structured-input diagnostics:

- **VERIFIED (exact dyadic oracle):** nested Float64 interpolation of a
  `(4,4,4)` block translated to `1e100` produced represented tetrahedron volume
  `4.251693490531567e256` for a corner determinant of
  `4.3388168547720914e256` (ratio `0.979920018024107`). The guarded exact
  affine path now produces `4.3388168547720954e256` (relative error
  `9.22096689867788e-16`) and a valid mesh. Its fixed connectivity SHA-256 is
  `cfdebd9e1af30eb255ed966e95bc3999f89d8062872d5cdb370647aa1737dfa8`.
- The exact path activates only when all eight represented corners satisfy the
  affine identities exactly; a separately accepted one-ULP residual therefore
  remains present at its output corner. Every unrecombined affine block now
  receives the same compensated exponent-scaled determinant audit, with exact
  fallback, previously used by the prism path. The shared implementation lives
  in `src/structured/StructuredNumerics.jl`. A separate noncollapsed remote
  lattice whose normalized determinant sum was `4.829750061035156` instead of
  `4.83782958984375` is now an explicit material-conservation blocker.
- Derived area/volume overflow is now an input `ArgumentError`, not an internal
  validation exception. Across 10,000 seeded remote-lattice candidates per
  generator, all 5,001 canonically oriented block cases and 5,012 prism cases
  either returned a valid mesh (1,359 blocks and 51 prisms) or a documented
  `ArgumentError`; no other exception type escaped.
- Structured patch, triangle, prism, volume, and hexahedron coordinates now
  reject inappropriate `Bool` values explicitly. Patch, volume, and prism
  allocation limits diagnose non-integers with `ArgumentError` instead of a
  keyword-dispatch `TypeError`.
- The seven focused structured suites passed 1,196/1,196 assertions under both
  Julia 1.12.7 and Julia 1.11.9. Under Julia 1.12.7, the volume allocation
  fixtures used 386,016 and 774,160 bytes; the prism fixtures retained 1,997,088
  and 4,008,288 bytes.
- The bounds-checked package gate passed 164,339/164,339 assertions in
  19m08.2s. Aggregate bounds-checked validation exited 0 in 17m20.0s against
  Gmsh 4.15.2, including the unchanged 72-node/144-tetrahedron affine-volume
  differential. Recursive package ambiguity detection and the public
  documentation scan both returned zero.
- The organized source layout remains enforced: the only top-level Julia files
  under `src`, `test`, and `validation` are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; all implementation and
  supporting Julia files remain categorized in subfolders.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening scalar and
anisotropic metric-length evaluation:

- **VERIFIED (exact dyadic oracle):** subtracting the endpoints
  `-floatmax(Float64)` and `floatmax(Float64)` previously overflowed before an
  isotropic size of `1e308` could normalize the edge. The certified result is
  now `3.5953862697246315`. A rotated metric whose two Cholesky products
  overflow and cancel now returns `1.6968532169535264e301`, matching direct
  exact evaluation of the stored `Metric3` quadratic form.
- Floating Cholesky rows retain an allocation-free fast path. Cancellation,
  subnormal products, non-finite intermediates, and ill-conditioned Schur
  complements fall back to exact IEEE-dyadic arithmetic and one outward
  `Float64` rounding. A 20,000-case seeded, exponent-varied audit produced
  13,778 representable valid metric/direction pairs, no mismatch above
  `2e-12` relative error, and maximum relative error
  `3.550982574471733e-13` against an independent 512-bit exact-rational oracle.
- Direction normalization now handles finite vectors whose unscaled Euclidean
  norm overflows. Coordinate midpoints matched a 256-bit oracle in all
  1,000,000 seeded finite pairs; nonzero metric lengths that underflow nearest
  rounding are conservatively represented by the minimum positive subnormal,
  while genuinely unrepresentable upper overflow is diagnosed.
- `DistanceField(mesh)` and `AutomaticMeshSizeField(mesh)` revalidate mutable
  source storage before connectivity indexing. Central field values, points,
  metric entries, and numeric constructor inputs diagnose inappropriate `Bool`
  values explicitly.
- The focused size-field suite passed 6,976/6,976 assertions under Julia 1.12.7
  and Julia 1.11.9. The bounds-checked package gate passed
  164,322/164,322 assertions in 12m15.2s, and aggregate bounds-checked
  validation exited 0 in 10m32.4s against Gmsh 4.15.2. Recursive package
  ambiguity detection and the public documentation scan both returned zero.
- The organized source layout remains enforced: the only top-level Julia files
  under `src`, `test`, and `validation` are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; implementation files remain
  categorized in their subfolders.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening the mixed-element
catalog, containers, and MSH contracts:

- **VERIFIED (ownership regression):** `MixedMesh` construction and
  `add_block!` retained the caller's `ElementBlock` arrays. `Base.mightalias`
  returned true in both paths, and changing the original connectivity made the
  stored mesh invalid. Ordinary and special blocks are now copied on every
  public insertion; coordinates, connectivity, CSR offsets, physical tags,
  parent/domain references, names, and entity metadata are detached from both
  caller storage and independently repeated results. The MSH reader retains a
  token-gated transfer path for its freshly allocated storage, avoiding a
  redundant full-file copy.
- The exported 125-entry `MSH_CATALOG` is now immutable, while an internal
  hash table retains the existing lookup path. The positional raw-storage
  constructor for `MixedEntityData` is sealed so documented construction cannot
  be bypassed accidentally.
- Element types, connectivity, tags, entity metadata, local-node orders, read
  limits, MSH versions, and Boolean controls now diagnose inappropriate `Bool`
  values explicitly. A 20,000-case seeded byte-mutation audit of ASCII and
  binary MSH seeds produced 283 still-readable files and 19,717 bounded
  `ArgumentError` rejections, with no unexpected exception type.
- The installed Gmsh 4.15.2 element differential retained all 950 assertions,
  and the fixed mixed-mesh CRC remains
  `b219f5afde8b589ce8c31c0fb174ebd2373811ab0f4e37a564f18934499e00c0`.
  Binary/ASCII, opposite-endian, sparse-tag, special-record, and all-125-type
  round trips, malformed-input gates, atomic output, and allocation ratchets
  all passed.
- The focused `Elements` suite passed 2,848/2,848 assertions under Julia 1.12.7
  and Julia 1.11.9. The bounds-checked package gate passed
  164,308/164,308 assertions in 12m11.5s, and aggregate bounds-checked
  validation exited 0 in 10m31.6s against Gmsh 4.15.2. The public `Elements`
  documentation scan returned no missing names; recursive package ambiguity
  detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening finalized-mesh
affine transformations:

- **VERIFIED (exact dyadic oracle):** finite cancellation around a remote pivot
  made the identity transform map `(1,0,0)` to `(0,0,0)` when the pivot was
  `(1e16,0,0)`. A conservative accumulated-error filter now sends cancellation,
  subnormal, overflowed-bound, and non-finite cases through exact rational
  evaluation before the result is rounded once to `Float64`; the regression now
  returns `(1,0,0)` exactly.
- Across 100,000 seeded, scale-varied affine coordinate expressions, comparison
  with an independent `Rational{BigInt}` oracle found maximum absolute error
  `3.084631262387123e-16` relative to the conservative expression scale, with
  50,186 results bit-exact. Identity transformation also preserves the minimum
  positive subnormal exactly.
- Transform controls now reject non-Boolean `check` values and non-real angles
  explicitly, every entry point revalidates mutable input storage, and returned
  coordinate, connectivity, and tag arrays are detached from the input and from
  independently repeated results.
- Affine results were invariant after normalization at scales `1e-100`, `1`, and
  `1e100`; an exact 90-degree rotation about a pivot of magnitude `1e100`
  preserved the expected topology and coordinates for a tetrahedron only 16
  coordinate ulps wide.
- The focused `Transform` suite passed 103/103 assertions under Julia 1.12.7 and
  Julia 1.11.9. Its fixed translated-mesh SHA-256 is
  `cfd2502be91e189981fa6a298a188e500c9180ee05866529897fb0a785b59737`.
- The bounds-checked package gate passed 164,244/164,244 assertions in 12m05.3s,
  and the aggregate bounds-checked validation exited 0 in 10m33.0s against Gmsh
  4.15.2. The public `Transform` documentation scan returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening quadratic
tetrahedra and type-11 solver output:

- **VERIFIED (Gmsh 4.15.2 API differential):** type-11 slots 9 and 10 were
  reversed. Tessella emitted edge `(2,4)` before `(3,4)`, while Gmsh requires
  `(3,4)` before `(2,4)`. Generation, shape gradients, edge ownership, curving,
  and MSH output now use Gmsh's ten-slot order, and Gmsh reads every written
  local node coordinate back identically.
- `P2Mesh` now owns and validates tetrahedron tags, `p2_tetmesh` preserves them,
  and `write_msh_p2` uses them by default. The high-order API is available from
  the top-level module, and every public consumer safely revalidates mutable
  coordinate, connectivity, and tag storage before bounds-elided access.
- **VERIFIED (subnormal regression):** the previous half-plus-half midpoint
  turned equal minimum-subnormal coordinates into zero and left a valid linear
  tetrahedron without a positive P2 Jacobian certificate. Correctly rounded
  midpoint construction now yields Jacobian `2e-323` and representable volume
  `5e-324`; node/tet allocation limits are checked before dense output arrays.
- Curving is transactional: a callback failure after one accepted projection
  previously left node 5 changed, while the same regression now restores every
  coordinate. Displacement scaling also handles a finite tetrahedron spanning
  `-floatmax(Float64):floatmax(Float64)` without first overflowing its bounding
  box diagonal.
- An independent exact-rational shape-function oracle matched all 19,600
  evaluations reconstructed from the cubic Bernstein coefficients. Cylinder
  curving from scales `1e-100` through `1e100` had maximum normalized coordinate
  difference `1.1102230246251565e-16`, and a 16-ulp-wide tetrahedron translated
  to `1e100` retained a positive exact certificate.
- The focused `HighOrder` suite passed 561/561 assertions under Julia 1.12.7 and
  Julia 1.11.9. Its fixed Gmsh-readable type-11 file SHA-256 is
  `5a83ebe0386bda71c6761148ed3fe2f964f16c2da2f0b66b6951ef558f4927ab`.
- The bounds-checked package gate passed 164,197/164,197 assertions in 12m09.3s,
  and the aggregate bounds-checked validation, including the new mandatory
  high-order Gmsh child, exited 0 in 10m30.2s against Gmsh 4.15.2. Public module
  and top-level documentation scans returned no missing names; recursive
  ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening one-level uniform
simplex refinement:

- Edge coordinates now use an overflow-safe, correctly rounded midpoint
  calculation. The previous subtract-then-add expression disagreed with a
  256-bit oracle in 1,274 of 200,000 seeded finite cases; the replacement had
  zero disagreements across 2,419,911 random and threshold-focused cases.
- Resource controls now diagnose every non-integer and Boolean value explicitly.
  Refined coordinate, connectivity, and tag arrays are detached from both the
  source mesh and independently repeated results.
- All 14,632 nondegenerate tetrahedra selected from the `3×3×3` integer lattice
  produced eight positive one-eighth-volume children, 16 boundary faces, and
  maximum face incidence two. The maximum relative parent/child volume-sum
  difference was `1.3322676295501878e-16`.
- The child topology and tags were invariant at scales from `1e-300` through
  `1e100` and for a tetrahedron translated to magnitude `1e100` with edges only
  16 coordinate ulps wide.
- The focused `Refine` suite passed 118/118 assertions under Julia 1.12.7 and
  Julia 1.11.9. The Gmsh 4.15.2 differential retained its fixed SHA-256
  `db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031`.
- The bounds-checked package gate passed 164,098/164,098 assertions in 12m02.3s,
  and the aggregate bounds-checked validation exited 0 in 10m26.9s against Gmsh
  4.15.2. The public `Refine` documentation scan returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening surface
triangle-to-quadrangle recombination:

- The Edmonds alternating-tree search now distinguishes an even-tree blossom
  edge from an undiscovered odd vertex. A five-vertex regression that previously
  indexed parent vertex zero now returns a consistent maximum matching.
- An independent subset-search oracle exhaustively checked all 33,868 simple
  undirected graphs through six vertices and a further 24,000 seeded graphs
  through eleven vertices, with zero cardinality or mate-consistency mismatches.
- Recombination now diagnoses non-Symbol algorithms and non-Boolean full-quad
  controls explicitly and rejects the incompatible greedy/full-quad combination
  before candidate construction. Returned coordinates, blocks, tags, and physical
  names are detached from caller storage.
- Strict square pairing retained identical ordered quadrangle connectivity at
  scales `1e-300`, `1e-150`, `1`, and `1e150`, and under three large-translation
  cases whose widths were only 16 coordinate ulps.
- The focused `Recombine` suite passed 92/92 assertions under Julia 1.12.7 and
  Julia 1.11.9. Its fixed 12-by-12 quadrangulation SHA-256 is
  `dbb1bf17965d4e011e7f51a452c6a03e4018a628ffc7e7d9b33d9fc6b922439f`.
- The bounds-checked package gate passed 164,058/164,058 assertions in 12m05.0s,
  and the aggregate bounds-checked validation exited 0 in 10m27.8s against Gmsh
  4.15.2. The public `Recombine` documentation scan returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening Gmsh-style curve
integration and grading:

- `curve_length`, `metric_length`, `mesh_curve`, and `mesh_segment` now reject
  Boolean coordinates and numeric controls, normalize entity contexts, bound
  integration and edge allocations, and validate every static control before
  invoking a caller-supplied curve.
- Uniform/adaptive parameter interpolation, primitive inversion, close-point
  sampling, and straight-segment coordinates use overflow- and cancellation-safe
  convex combinations. Straight segments preserve both supplied endpoints even
  when their magnitudes differ by hundreds of orders, and closed-curve tolerance
  is relative to sampled curve extent instead of a unit-scale floor.
- A 99,906-case finite interpolation audit had maximum error
  `3.769410006981428e-16` relative to the larger weighted term against a 256-bit
  oracle. Segment grading from scales `1e-300` through `1e300` differed after
  normalization by at most `8.881784197001252e-16`; closed-circle parameters at
  scales `1e-200`, `1`, and `1e200` differed by at most
  `2.220446049250313e-16`.
- The focused `Mesh1D` suite passed 104/104 assertions under Julia 1.12.7 and
  Julia 1.11.9. Its fixed graded-chain SHA-256 is
  `c88590e849684b244044860b05a509139b97e422d6a2d65074b19fd73b3f9048`.
- The bounds-checked package gate passed 164,031/164,031 assertions in 12m06.5s,
  and the aggregate bounds-checked validation exited 0 in 10m28.2s against Gmsh
  4.15.2, including all five mesh-observed size-field cases. The public `Mesh1D`
  documentation scan returned no missing names; recursive ambiguity detection
  returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening tetrahedral quality
reporting and mesh optimization:

- `TetQuality` now has one documented, validating construction path. All quality
  and smoother controls reject Boolean, nonfinite, negative, and platform-
  unrepresentable inputs as applicable; floating means are clamped to their
  measured extrema before report construction. `remove_slivers`, including a
  zero-round call, returns detached mesh storage and preserves every cell tag.
- ODT smoothing now solves each circumcenter in a dimensionless local frame,
  omits only shape-negligible tetrahedra, and accumulates physical-volume weights
  in the log domain. Cancellation-safe coordinate reconstruction and convex
  averaging cover finite extreme coordinates, and sorted neighbour traversal
  makes Laplacian updates deterministic.
- A 5,000-tetrahedron scale differential had maximum relative circumcenter error
  `6.258726052278117e-13` and maximum log-weight shift error
  `6.821210263296962e-13`. The same ODT mesh update at scales `1e-100`, `1`, and
  `1e100` differed after normalization by at most
  `5.862947357926637e-14`.
- The focused `Optimize` suite passed 96/96 assertions under Julia 1.12.7 and
  Julia 1.11.9. Its fixed one-step ODT SHA-256 is
  `31a280b6a063b428a11772b752ca1d9a64f70a8d4728029c23064a78e977ee00`.
- The bounds-checked package gate passed 164,001/164,001 assertions in 12m02.4s,
  and the aggregate bounds-checked validation exited 0 in 10m26.5s against Gmsh
  4.15.2. The public `Optimize` documentation scan returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening planar, cylindrical,
and parametric surface meshing:

- `PlaneFrame` now has one validating construction path. `plane_frame` normalizes
  coordinate scales, certifies the Newell filter, and uses an exact-rational
  orientation-preserving fallback. Projection and lifting reject malformed,
  Boolean, nonfinite, or unrepresentable values and use high-precision fallback
  under cancellation. A 200-frame audit from scales `1e-300` through `1e300`
  had maximum normal-direction loss `3.3306690738754696e-16` and maximum relative
  lift/project round-trip residual `1.7826267587088538e-16`; ordinary projection
  and lifting each allocated zero bytes.
- Planar inputs are copied into strict three-coordinate loops, resource-counted,
  and checked with a scale-relative coplanarity tolerance instead of a unit-scale
  floor. Isotropic as well as anisotropic final edges now reach the physical
  metric certificate.
- Cylinder construction validates all controls before meshing, evaluates
  overflow/cancellation-safe coordinates, rejects unrepresentable axial
  subdivisions, and post-certifies physical area, angle, and field-metric bounds.
  A 5,000-case scale-varied quality differential had maximum minimum-angle error
  `8.368306048112117e-13` degrees against a 256-bit oracle.
- General parametric patches now adapt and certify boundary and interior chord
  edges in physical space, preserve curve-vs-face entity context, enforce
  `max_area` as a physical triangle-area contract, check the sampled surface
  Jacobian, and reject mapped inversions. Five affine scale cases from `1e-150`
  through `1e150` had maximum requested-area ratio
  `0.898589065255732` and maximum field-edge metric `0.7544670215115625`.
- The focused `MeshSurface` suite passed 60/60 assertions under Julia 1.12.7 and
  Julia 1.11.9. Its fixed planar-square SHA-256 is
  `a0cfb73fe65d6814802e2df6d534a985c0bbb8d7e69a35eb860029f9d14a48ee`.
- The bounds-checked package gate passed 163,974/163,974 assertions in 13m15.2s,
  and the aggregate bounds-checked validation exited 0 against Gmsh 4.15.2.
  The public `MeshSurface` documentation scan returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening exact predicates and
the planar meshing workspace:

- `orient2` now uses exact dyadic evaluation when a nonzero determinant product
  underflows to zero. All four filtered predicates accept a floating result only
  with a normal finite error bound, and predicate coordinates and SoS indices
  reject Boolean inputs. An independent exact-rational differential covered 480
  cases from coordinate scales `1e-300` through `1e300` with zero sign
  mismatches; a 200-case Delaunay scale audit had zero topology mismatches.
- `Triangulation` now owns its input coordinates and has one validated public
  construction path. Public topology operations diagnose corrupt mutable state
  without bounds errors, validate indices and controls, preserve existing
  constraints during point insertion, and reject crossing or positive-overlap
  constraints before mutation. A nine-mutation corruption audit was rejected
  safely in every case.
- Segment recovery now uses exact collinearity and coordinate-wise betweenness,
  including a finite `floatmax`-scale through-vertex case. Circumcenters and
  radius-edge ratios use a scale-normalized frame with exact-rational fallback;
  2,500 scale-varied triangles had maximum circle residual
  `1.409222853971538e-15` and maximum relative radius-edge error
  `7.899298712677908e-14`.
- The focused predicates suite passed 140,255/140,255 assertions under Julia
  1.12.7 and Julia 1.11.9. The focused Mesh2D suites passed 198/198 assertions
  under both versions. The fixed square SHA-256 is
  `850fe31fb8b9c7946d716633cfabdfaf13850456a1b53474d21edfcfa9f194f4`.
- The bounds-checked package gate passed 163,945/163,945 assertions in 17m32.1s,
  and the aggregate bounds-checked validation exited 0 against Gmsh 4.15.2.
  Public `Predicates` and `Mesh2D` documentation scans returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-25 with Julia 1.12.7 after hardening finalized mesh
storage, raw topology operations, and tetrahedron quality metrics:

- `Mesh` construction now rejects Boolean coordinates, connectivity, and tags;
  checks tag representability with cell-local diagnostics; and retains the
  allocation-free `Int`/`Int32` node-access paths. Public raw-topology
  operations validate one-based shape, Boolean/range errors, and topology-size
  limits before entering bounds-elided loops.
- Finalized mesh arrays remain intentionally mutable, so every public consumer
  now revalidates structural invariants. `validate` returns a diagnostic instead
  of indexing corrupt storage; topology/CRC operations reject it explicitly;
  manifold queries return `false`; and `MeshDiagnostic` owns its messages while
  enforcing a consistent success/failure state.
- Tetrahedron dihedral, circumradius, and radius-edge calculations normalize
  overflowed and subnormal coordinate scales. A finite radius-edge ratio is
  retained when the corresponding physical circumradius legitimately
  overflows. Random scale differentials covered 983 finite huge-coordinate
  cases with maximum angle error `3.0487765090292385e-14`, maximum relative
  radius-edge error `7.172040739078511e-14`, and maximum relative finite-radius
  error `8.659739592076221e-15`.
- The focused `MeshTypes` suite passed 1,993/1,993 assertions under Julia 1.12.7
  and Julia 1.11.9. Its fixed cube SHA-256 is
  `7ea403054f05392f18b404a1f5f78b12d70d45d40c7b04ba8f8dc3e030d8f3f9`.
- The bounds-checked package gate passed 163,864/163,864 assertions in 47m40.5s,
  and the aggregate bounds-checked validation exited 0 against Gmsh 4.15.2.
  The public `MeshTypes` documentation scan returned no missing names and
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening and separating the
API, CLI, and headless GUI interfaces:

- API sessions and model/cache operations are now serialized. Initialization
  resets options, model, and cache; option updates validate atomically; model
  mutations invalidate the cached mesh only after success; generated, cached,
  and `.geo`-returned meshes/models do not share caller-mutable storage.
- The CLI now bounds arguments, rejects duplicate/conflicting or ignored flags,
  derives a safe `.msh` destination for no-extension inputs, and blocks lexical
  and hard-link aliases of the input before executing or writing. The
  no-extension regression leaves its source byte-for-byte unchanged.
- `GuiState` now owns and validates constructor inputs. Commands have byte,
  token, selection, and log bounds; exact arities; finite numeric parsing;
  positive unique selections; atomic state updates; and explicit pre-session
  blockers.
- The focused API, CLI, and GUI suites passed respectively 58/58, 28/28, and
  51/51 assertions under Julia 1.12.7 and Julia 1.11.9. A four-thread API stress
  probe allocated 200 unique automatic point tags without a race. The API cube
  and CLI square deterministic SHA-256 values are
  `e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83`
  and `92e578bac6d8feb3f0f845f100665dcc145edf924965ad72be716da933f34461`.
- The bounds-checked package gate passed 163,815/163,815 assertions in 10m46.2s,
  and the aggregate bounds-checked validation exited 0 against Gmsh 4.15.2.
  Public API, CLI, and GUI documentation scans returned no missing names;
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening and separating the
scalar post-view interface:

- `View` now has one validating construction path for mesh-backed or direct
  `3 × n` coordinates, checks Float64 conversion and finiteness, and owns copies
  of both coordinates and samples. This closes the former exact-field-type
  constructor bypass and prevents caller-array mutation from changing a view.
- `view_value` now rejects Boolean, platform-unrepresentable, and out-of-range
  indices. Plugin registration is synchronized, returns a canonical `String`,
  and documents replacement semantics; plugin lookup releases the registry lock
  before executing user code.
- Post tests now live independently in `test/interfaces/post_test.jl`. The
  focused suite passed 26/26 assertions under Julia 1.12.7 and Julia 1.11.9; its
  deterministic scalar-view SHA-256 is
  `56e766682618b76029640b47caa692205eb967a97438ee383eb057fc47cd96cd`.
- The bounds-checked package gate passed 163,693/163,693 assertions in 10m39.7s,
  and the aggregate bounds-checked validation exited 0 against Gmsh 4.15.2.
  Public `Post` and top-level documentation scans returned no missing names, and
  recursive ambiguity detection returned zero.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening public
tetrahedralization and Steiner insertion:

- `tetrahedralize` now has its public documentation attached to the actual
  method, validates basic surface structure even under the expert `check=false`
  path, rejects tetrahedra and unreferenced surface nodes, bounds and normalizes
  interior-point iterables before meshing, and rejects Boolean random seeds and
  resource limits.
- `insert_steiner3` now validates its input and three-coordinate point contract,
  distinguishes exact duplicates, checks Int32 growth before allocation, and
  certifies both the returned topology and total-volume conservation. Its
  containing-tet tolerance is relative to the tet volume instead of carrying a
  unit-scale absolute floor, so a `1e-15`-edge tet distinguishes an interior
  point from a vertex and rejects a point outside the volume.
- The focused Mesh3D suite passed 473/473 assertions under Julia 1.12.7 and
  Julia 1.11.9. The tiny-tet insertion, certified cube fill, and cube fill with
  one interior point have deterministic SHA-256 values
  `71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b`,
  `e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83`,
  and `4f0d7f17865d02bc785bb2a22b30b7e7de826b771f91ff0b18490657e44fb472`.
- The bounds-checked package gate passed 163,671/163,671 assertions in 10m56.9s,
  and the aggregate bounds-checked validation exited 0, including the point,
  line, and sheet embedding differentials that exercise Steiner recovery.
- The top-level public documentation scan now returns no missing names;
  `tetrahedralize` and `insert_steiner3` are documented in `Mesh3D`, and recursive
  ambiguity detection returned zero.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening the boundary-layer
entry points:

- All three public boundary-layer operations now reject invalid input meshes,
  Boolean or platform-unrepresentable integer controls, nonfinite real controls,
  and overflowed node/cell counts before allocation. Geometric layer offsets are
  accumulated without the cancellation in the former closed-form expression.
- The 2-D path now honors the supplied segment direction when choosing the left
  side, requires one coherently directed chain with no unreferenced nodes, bounds
  fan indices/counts, and checks the predicted node count. Filled layers also
  reject nonfinite wall volumes and coincident offset-cap nodes before core
  recovery.
- The focused boundary-layer suite passed 67/67 assertions under Julia 1.12.7
  and Julia 1.11.9. The deterministic prismatic and 2-D strip SHA-256 values are
  `f37a2141b37471b366cdbcaa5b1ede69c9088833b0f54e142aaa945e4ea23651`
  and `bb9e1fb9a0f0e56de42287ff3f85dd93ea3c7115cc9e51d05a6845d97ce8122b`.
- The bounds-checked package gate passed 163,637/163,637 assertions in 11m32.4s,
  and the aggregate bounds-checked validation exited 0. The direct Gmsh 4.15.2
  boundary-layer differential reported
  `tessella_area=0.15960000000000005`, `gmsh_quads=10`, `gmsh_tris=24`, and
  `tessella_quads=3`.
- The `BoundaryLayer` public documentation scan returned no missing names,
  leaving only `tetrahedralize` undocumented at top level; recursive ambiguity
  detection returned zero.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening translated
periodic node-pair certification:

- `periodic_identify` now validates a finite nonzero translation, a finite
  nonnegative tolerance, representable indices and translated coordinates, and
  a one-to-one mapping whose master and slave sets are unique and disjoint. All
  pairs are certified before the copied mesh is mutated; node numbering,
  connectivity, and tags remain unchanged while slave coordinates are snapped
  exactly.
- The periodic unit suite passed 18/18 assertions under Julia 1.12.7 and Julia
  1.11.9. Its deterministic snapped-mesh SHA-256 is
  `2d4c3e493639ced1a3a2e21a948e73c36b068c18cd37781760faf32eadd8f6f0`.
- The bounds-checked package gate passed 163,604/163,604 assertions in 11m48.0s.
  Recursive ambiguity detection returned zero; the `Periodic` module's public
  documentation scan returned no missing names, leaving only
  `mesh_boundary_layer` and `tetrahedralize` undocumented at top level.
- The aggregate bounds-checked validation exited 0. Its new Gmsh 4.15.2
  translation-periodic curve differential certified five node pairs with maximum
  pre-snap error `2.0594637106796654e-12` and deterministic SHA-256
  `baa96c7ebc0265667209f1940c77d5bdeed5ecb8a12f765d02df9d1945373648`;
  every existing child passed unchanged.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening native `.geo`
execution and primitive translation:

- `translate_volume!` now provides the model-level translation path for boxes,
  cylinders, spheres, and cones. It checks offset shape/finite conversion and
  transformed-coordinate representability before mutation; unsupported Boolean
  encodings remain explicit blockers.
- Boolean operand `Delete` is now applied independently. The selective case
  `{ Volume{1}; Delete; }{ Volume{2}; }` leaves volumes `[2,3]`, matching a direct
  Gmsh 4.15.2 API probe (`[(3,2),(3,3)]`); malformed operand suffixes block.
- The brace-aware executor now handles `/* ... */` comments, rejects unmatched
  braces and unterminated comments, and bounds individual statements and the
  statement count. `mesh_dim` rejects `Bool`, out-of-range integers, and values
  outside `{0,2,3}` before parsing the file.
- Two deleted-temp-file false positives in the `.geo` blocker tests were repaired;
  the intended multi-volume and `Extrude` paths are now exercised while the files
  still exist.
- The bounds-checked package gate passed 163,590/163,590 assertions in 12m25.2s.
  The focused model/`.geo` suite passed 146/146 under Julia 1.12.7 and 1.11.9;
  public Model, GeoExec, and top-level documentation scans now leave only the
  three unrelated meshing exports, and recursive ambiguity detection remains 0.
- The aggregate bounds-checked validation exited 0 against Gmsh 4.15.2-git,
  including every size-field, transfinite, API, CAD, NURBS/IGES, embedding,
  Boolean, boundary-layer, and analytic-volume child.

Re-measured on 2026-08-24 with Julia 1.12.7 after completing the classified
STEP/IGES interoperability increment:

- STEP parsing now rejects duplicate/out-of-range identifiers and nonfinite data;
  point-cloud block classification cannot silently replace mixed topology, and
  multi-solid primitive imports block explicitly. Complex rational STEP surfaces
  now import with their full weight matrix, closing the prior curve-only complex
  rational path.
- IGES input now observes the standard 64-column parameter field instead of
  consuming the directory pointer, accepts `D` exponents, rejects malformed or
  unterminated numeric records, checks integer/count arithmetic before allocation,
  and blocks multiple recognized solids.
- IGES 126/128 export now writes atomic, exact 80-column S/G/D/P/T sections with
  directory entries and an untrimmed type-144 wrapper for type-128 surfaces.
  Mutated/invalid NURBS objects are rejected before replacement. The deterministic
  curve-plus-surface export SHA-256 is
  `ae515df934189f3d0b3cf5614bd39427cd6d015ceb07a9858c44fc32ea7986f5`.
- The bounds-checked package gate passed 163,569/163,569 assertions in 19m18.4s.
  A final checksum-only assertion was then added without changing production code;
  the STEP/IGES suite passed 83/83 under Julia 1.12.7 and Julia 1.11.9.
- The aggregate bounds-checked validation exited 0 against Gmsh 4.15.2-git.
  Tessella recovered the centre of Gmsh's IGES-128 patch, while Gmsh imported
  Tessella's standalone IGES-126 curve and IGES-128/144 surface and meshed the
  latter to area 1 (`export_nodes=50`, `export_tris=66`). Every other required
  differential and parity child passed unchanged.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening the public NURBS
evaluation path:

- Curve/surface degrees now reject Boolean, non-integer, negative, and
  platform-unrepresentable values; knot vectors require a finite, sorted,
  positive-width active interval. Control points have a strict iterable
  three-coordinate contract.
- `bspline_basis` now validates its complete public input contract, handles
  degree zero, clamped endpoints, non-clamped vectors, and out-of-domain
  parameters with an iterative Cox–de Boor evaluation instead of recursive
  failure paths. Homogeneous evaluation normalizes globally scale-invariant
  weights before multiplication and rejects only unrepresentable relative weights
  or results.
- The bounds-checked package gate passed 163,536/163,536 assertions in 20m41.4s.
  Seven degree-zero/non-clamped test-only regressions were then added; the final
  focused NURBS suite passed 60/60 assertions under Julia 1.12.7 and Julia 1.11.9.
  The STEP/IGES suite passed 56/56.
- The aggregate bounds-checked validation exited 0 against Gmsh 4.15.2-git. Its
  IGES-128 differential reported `tessella_centre=0.5`, `gmsh_area=1`, and
  `gmsh_tris=164`; every other required child also passed unchanged.
- The recursive ambiguity scan returned zero, the `NURBS` module's public
  documentation scan returned no missing names, and `git diff --check` passed.

Re-measured on 2026-08-24 with Julia 1.12.7 after hardening and documenting the
native geometry/entity model:

- Signed curve-loop orientation now accepts negative curve references, verifies
  ordered endpoint continuity and closure before mutation, and has a deterministic
  reversed-square mesh CRC. Embeddings validate atomically; physical groups require
  existing entities and allocate independently of entity tags; returned physical
  memberships no longer expose mutable model storage.
- Tags, dimensions, primitive radii/axes, transform vectors, transform results, and
  automatic-tag exhaustion have explicit finite/range checks. Surface refinement now
  includes hole and embedded-point characteristic lengths; the embedded-size
  regression has deterministic SHA-256
  `13917dad18b19e8376640a379a2f1cd338aacca5b01f66e100b5bc372fc91371`.
- The bounds-checked package gate passed 163,499/163,499 assertions in 21m30.3s.
  After the gate, four additional deterministic/overflow regression assertions were
  added; the focused model suite then passed 126/126 assertions under Julia 1.12.7
  and Julia 1.11.9. No executable production code changed after the full gate.
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  exited 0 against Gmsh 4.15.2-git with every required differential and parity child
  passing. The recursive ambiguity scan returned zero, the `Model` module's public
  documentation scan returned no missing names, and `git diff --check` passed.

Re-measured on 2026-08-24 with Julia 1.12.7 after the repository-wide Julia-file
organization:

- The 37 implementation files below `src/Tessella.jl` are divided among `core`,
  `fields`, `geometry`, `interfaces`, `meshing`, and `structured`. The 33 package
  test files below `test/runtests.jl` mirror those domains and add `integration`;
  the shared validation harness is in `validation/support`. The only top-level
  Julia files in those three trees are their entry points.
- A nine-assertion repository-layout ratchet now rejects stray top-level Julia
  files, missing source/test domains, and deeper unclassified Julia paths.
- The final bounds-checked package run passed 163,454/163,454 assertions in
  15m46.4s without method-overwrite warnings. The same reorganized package loaded
  under Julia 1.11.9, whose focused curve suite passed 511/511 assertions.
- Path-sensitive focused gates passed: predicates 140,230/140,230; STEP/IGES 56/56;
  I/O 305/305; Mesh3D 439/439; size fields 6,962/6,962; HFSS 95/95; combined
  structured tests 1,179/1,179; and combined geometry tests 125/125.
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  exited 0 against Gmsh 4.15.2-git after every validation include and fallback
  path was updated. All differential checksums and parity measurements remained
  unchanged.

Re-measured on 2026-08-24 with Julia 1.12.7 after adding the normalized
`Bump_HWall` and `Beta_HWall` curve laws:

- `julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  passed 163,445/163,445 assertions in 16m08.3s.
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  exited 0 against Gmsh 4.15.2-git. The aggregate included 39 straight-curve
  cases, 18 HWall cases, 393 coordinates, and all existing external validation
  children. The HWall differential SHA-256 was
  `7b02b56c7bfc66becafff0793e66df0e24434f15628b24ecb4b555303405cea7`,
  with maximum absolute coordinate error `7.155441150707986e-8`.
- The focused bounds-checked curve suite passed 511/511 assertions under both
  Julia 1.12.7 and Julia 1.11.9. Independent 256-bit primitive-inversion oracles
  covered Bump/Beta HWall laws at 6, 17, and 65 nodes, both Beta orientations,
  near-uniform limits, extreme physical scales, invalid and
  Float64-unrepresentable inputs, and pre-allocation resource rejection. The
  deterministic Bump/Beta HWall test SHA-256 was
  `04169f75cdcfabf540e88477ce54b55d0a6eabcfd5ca48f19b332afaac0fb59a`.
- Allocation ratchets measured 163,904 and 327,744 bytes for 20,000 and 40,000
  nodes for each of the Progression, Bump, and Beta HWall paths.
- The recursive method-ambiguity scan returned zero, the public HWall API was
  documented, and `git diff --check` passed.

Re-measured on 2026-08-24 with Julia 1.12.7 after the normalized
`Progression_HWall` curve-law increment:

- `julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  passed 163,358/163,358 assertions in 15m12.9s. An earlier independent run of
  the same package gate passed 163,356 assertions before the final two robustness
  assertions were added.
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  exited 0 twice against Gmsh 4.15.2-git. The final run preserved the exact flat
  model volumes (box 2, tunnel 24, hollow box 35), the cylinder-prism volume
  62.652572, and the enclosure fixture's measured zero Gmsh volume tetrahedra.
- The focused CRC passed 424/424 bounds-checked assertions under Julia 1.12.7
  and Julia 1.11.9. Its independent 256-bit geometric-sum oracle covered both
  orientations, uniform/large-ratio/near-uniform cases, extreme physical scales,
  invalid and Float64-unrepresentable inputs, pre-allocation resource rejection,
  and linear allocation growth (163,904 and 327,744 bytes for 20,000 and 40,000
  nodes). The deterministic HWall parameter SHA-256 is
  `802ae6dd95259c50b087d03e7b7567b555f6040f8e62b1a2afc3aed6bca22379`.
- The required Gmsh differential passed six `Progression_HWall` cases in both
  orientations as part of 27 total straight-curve cases and 277 coordinates;
  maximum absolute error remained `5.487045007246394e-8`, and the HWall SHA-256
  was `8d79e5323a0c4d7c635b4101bad3f1325f33e0badcded389f5b3bd36ea509213`.
- The recursive method-ambiguity scan returned zero, the new public HWall API was
  documented, and `git diff --check` passed.

Re-measured on 2026-08-21 with Julia 1.12.7 after 2-D boundary-layer
quad/fan topology and the P6 BL-quads corpus. Both bounds-checked package
runs matched:

- `julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 163,235/163,235 assertions passed twice (10m16.5s, then 10m16.1s).
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  — exited 0 against Gmsh 4.15.2-git. Exact flat-model volumes box=2, tunnel=24,
  hollow box=35; cylinder prism 62.652572; enclosure gmsh empty solids reproduced.
- Size-field child: `SIZE_FIELD_DIFFERENTIAL_OK gmsh=4.15.2 plugin_calls=23
  direct_cases=23 direct_samples=63 mesh_cases=5 context_skips=5`.
- Geo-range child: `GEO_RANGE_DIFFERENTIAL_OK gmsh=4.15.2-git float_cases=13
  integer_cases=4 wrapped_cases=2 samples=58 bit_exact=1`.
- Transfinite hex child: `TRANSFINITE_HEX_DIFFERENTIAL_OK gmsh=4.15.2-git cases=8
  nodes=288 hexahedra=96 boundary_quadrangles=256 max_node_error=9.50e-12`.
- P6 box API child: `GMSH_PARITY_BOX_OK gmsh=4.15.2 tessella_volume=1
  gmsh_tets=1158 tessella_tets=12`.
- P6 t1 child: `GMSH_PARITY_T1_OK gmsh=4.15.2 tessella_area=1 gmsh_tris=14
  tessella_tris=16`.
- P6 t4 hole child: `GMSH_PARITY_T4_HOLE_OK gmsh=4.15.2 tessella_area=0.75
  gmsh_tris=12 tessella_tris=12`.
- P6 embed child: `GMSH_PARITY_EMBED_OK gmsh=4.15.2 tessella_area=1 gmsh_tris=16
  tessella_tris=16 tessella_nodes=13`.
- P6 embed-line child: `GMSH_PARITY_EMBED_LINE_OK gmsh=4.15.2 tessella_area=1
  gmsh_tris=22 tessella_tris=22 tessella_nodes=18`.
- P6 embed-sheet child: `GMSH_PARITY_EMBED_SHEET_OK gmsh=4.15.2 tessella_volume=1
  gmsh_tets=904 tessella_tets=44 tessella_nodes=17`.
- P6 2-D boundary-layer child: `GMSH_PARITY_BL2D_OK gmsh=4.15.2
  tessella_area=0.15960000000000005 gmsh_quads=10 gmsh_tris=24 tessella_quads=3`.
- P6 cylinder child: `GMSH_PARITY_CYLINDER_OK gmsh=4.15.2
  tessella_prism=6.211657082460498 gmsh_tets=60 tessella_tets=96`.
- P6 boolean-boxes child: `GMSH_PARITY_BOOLEAN_OK gmsh=4.15.2 tessella_volume=1
  gmsh_tets=100 tessella_tets=12`.
- Focused CRC: model/entity/`.geo` 77/77.
- `git diff --check` passed. `.grok/` is gitignored and absent from the index.

Previous aggregate on the same day with Julia 1.12.7, kept as historical:

- `julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 161,183/161,183 assertions passed in 9m05.3s.
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  — exited 0 against Gmsh 4.15.2-git. Tessella preserved the exact flat-model
  volumes (box 2, tunnel 24, hollow box 35), completed the curved-model comparisons,
  and reproduced the literal enclosure's non-zero Gmsh exit with zero volume
  tetrahedra.
- The required size-field child reported
  `SIZE_FIELD_DIFFERENTIAL_OK gmsh=4.15.2 plugin_calls=23 direct_cases=23
  direct_samples=63 mesh_cases=5 context_skips=5`; the context skips are explicit
  non-claims listed in `validation/size_fields/STATUS.md`.
- Focused bounds-checked gates passed 6,915/6,915 size-field assertions,
  2,398/2,398 fixed-node/mixed-element assertions, and 74/74 Mesh1D assertions.
  The Mesh1D gate also passed under Julia 1.11.
- `detect_ambiguities(Tessella; recursive=true)` and
  `Base.Docs.undocumented_names(Tessella; private=false)` both returned zero.
- `git diff --check` passed.

After that stable aggregate gate, the isolated native-primitive increment passed
69/69 bounds-checked geometry assertions under both Julia 1.12.7 and Julia 1.11,
including deterministic CRCs, analytical/polyhedral volume checks, direct volume
meshing, resource/error paths, and linear allocation-growth ratchets. It will be
included in the next aggregate package gate. The subsequent finalized-mesh transform
increment passed 56/56 bounds-checked assertions under Julia 1.12.7 and Julia 1.11.9.
Its independent Gmsh 4.15.2 `model.mesh.affineTransform` oracle matched all four
fixture nodes exactly; 10,000- and 20,000-node translations allocated 504,984 and
996,504 bytes respectively, and the method-ambiguity scan remained empty. The binary
mixed-MSH increment then passed 2,579/2,579 bounds-checked assertions under both Julia
1.12.7 and Julia 1.11.9, including native and opposite-endian MSH 2.2/4.1, full-width
v4 tags, parametric nodes, atomic/resource failures, and Gmsh 4.15.2 acceptance.
Reading 2,000 and 4,000 nodes allocated 1,670,704 and 4,239,472 bytes respectively.
The following recombination increment passed 44/44 bounds-checked assertions under
Julia 1.12.7 and Julia 1.11.9. Its 12×12 grid produced 144 quadrangles with CRC
`dbb1bf17965d4e011e7f51a452c6a03e4018a628ffc7e7d9b33d9fc6b922439f`; 20×20 and
40×40 grids allocated 4,073,984 and 17,348,496 bytes. Gmsh 4.15.2 accepted both the
ASCII and binary recombined MSH 4.1 fixtures with `-check -parse_and_exit`.

The subsequent `.geo` constant-expression increment passed 231/231 bounds-checked IO
assertions under Julia 1.12.7 and Julia 1.11.9. An independent Gmsh 4.15.2 expression
oracle matched 34 accepted expressions exactly and matched seven error cases; an API
oracle also matched explicit expression-derived physical and field tags. Unsupported
control-flow contexts are rejected when they can affect relevant statements, and
their scalar bindings are invalidated before later use.

The next PostView increment passed 6,953/6,953 bounds-checked size-field assertions
under Julia 1.12.7 and Julia 1.11.9. It covers first-order scalar/vector point, line,
triangle, quadrangle, tetrahedron, hexahedron, prism, and pyramid list data; tensor
views preserve Gmsh's scalar `MAX_LC` result. The Gmsh 4.15.2 pointwise oracle covered
160 samples with maximum absolute error `6.66e-16`, plus eight exact closest-node
fallbacks. Warm 10,000-query loops for quadrangle, hexahedron, prism, pyramid, and
vector-quadrangle fields allocated at most 64 bytes in total.

The subsequent uniform-refinement increment passed 78/78 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 API differential matched
the ordered 2/4/8 segment/triangle/tetrahedron child templates and physical tags; the
combined fixture produced 10 nodes and CRC
`db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031`.
Refining 2,000 and 4,000 connected segments allocated 423,504 and 845,648 bytes,
respectively (1.99679×). The focused module ambiguity and public-doc scans both
returned zero.

The following special-element increment passed 2,822/2,822 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. It covers types 34/35/67/68/69/70/133–136,
compact variable connectivity, parent/domain references, validation and CRC, MSH2
ASCII and fixed-width binary records, and fixed unlinked MSH4 records. Pinned source
hashes and installed Gmsh 4.15.2 probes cover the accepted formats and their explicit
blockers. Rejected 30,000- and 300,000-entry variable records allocated 5,328 and
5,120 bytes after warm-up with `max_connectivity=0`; accepted 2,000- and 4,000-record
fixtures allocated 13,000,112 and 26,462,288 bytes. The tests pin Gmsh's type-69
normal-lifecycle crash and binary distinct-parent rewrite corruption as external
limitations instead of claiming unsafe compatibility.

The subsequent transfinite increment passed 91/91 bounds-checked assertions under
Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 API differential matched all
four triangle arrangements and 80 boundary/interior coordinate samples with maximum
absolute node error `4.44e-16`; it also pins Gmsh's mismatched-opposite-side fallback
at 18 nodes/23 triangles and its holed-surface error. The 64×64 and 128×64 allocation
fixtures used 1,771,392 and 3,342,496 bytes under Julia 1.12.7. Focused ambiguity and
public-documentation scans returned zero.

The subsequent straight-curve increment passed 310/310 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential covered
Progression/Power, Bump, and Beta in 21 cases and 219 coordinates with maximum
absolute error `5.49e-8`; its deterministic parameter SHA-256 is
`b525628945d55f49bc5151ad313d697f9c32f7b3325015e5d0080a69260ffd0e`.
The 20,000- and 40,000-node fixtures allocated 163,904 and 327,744 bytes for each
law; focused ambiguity and public-documentation scans returned zero.

The subsequent three-sided transfinite increment passed 103/103 bounds-checked
assertions under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential
matched exact boundary/triangle topology and 150 nodes across four arrangement names,
four resolutions, and planar/tilted geometries, with maximum coordinate error
`2.81e-15`. Its deterministic mesh SHA-256 is
`5f231f0c22f6812c247514103e21653e15281194911bf456b4c00a9d190a40df`.
The 64- and 128-division fixtures allocated 2,383,144 and 6,154,008 bytes under Julia
1.12.7; focused ambiguity and public-documentation scans returned zero.

The subsequent affine-volume increment passed 74/74 bounds-checked assertions under
Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential matched all 144
tetrahedra, 128 boundary triangles, and 72 mapped nodes across two affine blocks, with
maximum coordinate residual `9.51e-12`. The 8×8×4 and 16×8×4 fixtures allocated
385,984 and 774,128 bytes under Julia 1.12.7. Focused ambiguity and public-documentation
scans returned zero.

The subsequent recombined-quadrangle increment passed 128/128 bounds-checked
assertions under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential
matched ordered type-1/type-3 connectivity and 151 nodes across four arrangement
names, four resolutions, and planar/tilted geometries, with maximum coordinate error
`5.662137425588298e-15`. Its deterministic mixed-mesh SHA-256 is
`d05e9cbc57de975f1f88d8fb1da62e0636ea4bbe63a038d3ca060b16201b0ea2`.
The 64×64 and 128×64 fixtures allocated 3,430,528 and 6,274,688 bytes under Julia
1.12.7. Exact-rational fallback regressions cover valid thin affine patches whose
Float64 plane-normal or normalized-projection calculations cancel; focused ambiguity
and public-documentation scans returned zero.

The subsequent five-face-prism increment passed 130/130 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential matched
the exact canonical tetrahedron vertex order and boundary topology for 117
tetrahedra, 106 boundary triangles, and 63 nodes across three affine prisms; maximum
coordinate error was `6.840577468974138e-12`. The deterministic mesh SHA-256 is
`a16de779890f62f8a09d928cbef67a6f13b09c6765a7d91ce8e86de78c14db6e`.
The 24×12×6 and 48×12×6 fixtures allocated 1,997,056 and 4,008,256 bytes under
Julia 1.12.7. Exact-predicate and exponent-scaled volume gates reject represented
folds or material loss while accepting verified subnormal affine prisms; focused
ambiguity and public-documentation scans returned zero.

Earlier measurements on 2026-08-14 with Julia 1.12.6, kept because they were
not re-run in this gate:

- `julia --project=. --check-bounds=yes validation/run_all.jl`
  — passed in 182.30 s against Gmsh 4.15.2-git, including exact box/tunnel/hollow-box
  volumes, curved-reference comparisons, and the literal enclosure's reproduced
  non-zero Gmsh exit with zero volume tetrahedra.
- Focused gates passed: 983/983 field/refinement assertions, 66/66 existing 1-D and
  surface assertions, and 439/439 Mesh3D assertions.
- Final-code enclosure field acceptance with the literal graph clamped to
  `size_min=0.002`, `size_max=0.012` returned a valid mesh with 706,962 nodes,
  3,838,729 tetrahedra, and 113,350 tagged triangles in 680.61 s, with a
  1,172,013,056-byte maximum resident set. The fixture's literal 0.2667 mm minimum
  remains a separate full-resolution scalability gate and is not claimed here.

## Previous complete baseline gate

Verified on 2026-08-14 with Julia 1.12.6:

- `julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 151,787/151,787 assertions passed in 6m54.1s.
- `julia --project --check-bounds=yes validation/run_all.jl`
  — completed against gmsh 4.15.2-git. Tessella matched the exact model volumes for
  the box (2), box tunnel (24), hollow box (35), and 48-gon cylinder prism
  (62.652572). The external enclosure fixture reproduced gmsh's non-zero exit and
  zero tagged volume tetrahedra.
- A mixed segment/triangle/tetrahedron MSH v4.1 written by Tessella passed
  `gmsh <file> -check -v 3` with exit 0.
- `using Test, Tessella; detect_ambiguities(Tessella; recursive=true)` reported 0
  method ambiguities. `Base.Docs.undocumented_names(Tessella; private=false)`
  reported 0 undocumented public names.
- `git diff --check` passed, and the source/test scan found no `TODO`, `FIXME`,
  `WIP`, placeholder, or unimplemented markers.

The last runtime-bearing baseline head for this gate is
`27010e312455fe3dbb67ce0efd2f51800487c3ad`, and `origin/main` was checked to the
same SHA after its push. Later commits change documentation/docstrings only; their
package-load/focused checks and remote synchronization were also verified.

## Deep-debug closure

The cumulative audit repaired and regression-pinned the following classes:

- exact Delaunay completeness and fallback certification;
- top-level PLC/interface conformity and full triangle/tetrahedron vertex-link
  manifold checks, including pinches and duplicate cells;
- finite/range/resource contracts for public geometry, meshing, refinement,
  optimization, P2, and I/O APIs;
- exact coplanar 3-D circle decisions, consistent bounded SoS evaluation, robust
  extreme-magnitude area/volume signs, and Float64-resolution blockers;
- 2-D CDT/refinement midpoint, encroachment, duplicate-constraint, callback, and
  termination behavior;
- lower-dimensional topology/tag preservation during volume refinement;
- strict MSH section/entity/count/tag validation, atomic writers, and structural
  STL parsing/welding;
- exact global cubic-Bernstein P2 Jacobian certification and exact coefficient
  volume integration;
- analytical projection/imprint edge cases and generated-surface postconditions;
- interrupt propagation through conversion/overflow error paths.

The package's recover-or-block contract intentionally reports an explicit diagnostic
when a requested mesh exceeds representational or configured resource limits; that is
a safety contract, not a silent fallback. The closure statement applies only to the
historical Stage 0–6 scope, not to the active Gmsh parity objective.

## Memory evidence

The reproducible stress measurement uses
`refine_to_size(mesh_box(0,3,0,3,0,3; hmax=3), 0.25)` after one warm-up and `GC.gc()`.
Each measured snapshot produced 22,065 nodes, 98,304 tetrahedra, and a valid mesh:

| Version | Allocated bytes | Timed run |
|---|---:|---:|
| pre-fix `5d3e466` source snapshot | 75,038,688 | 0.3517 s |
| previous hardened implementation | 70,337,184 | 0.2972 s |
| previous field/refinement snapshot | 100,434,688 | 1.1318 s (concurrent enclosure load) |
| current P1/P2 worktree | 72,870,624 | 0.3027–0.3169 s |

The historical hardened version reduced allocation by 6.27% from `5d3e466`. The
previous field/refinement snapshot was 30,097,504 bytes (42.79%) above that version on
the same mesh. The current worktree is 27,564,064 bytes (27.44%) below that snapshot
and 2,533,440 bytes (3.60%) above the historical hardened implementation. All three
current runs produced 22,065 nodes, 98,304 tetrahedra, a valid mesh, and CRC
`7378bf6e460c596aa04f9f3b4a8cc9ce176a69f2386f253f7bd25057b551ceb9`.
The refinement queue keeps one live heap record per long edge via a companion set,
preventing duplicate queued records; the benchmark verifies the resulting mesh rather
than relying on allocation alone.

## Acceptance cases

- Exact predicates remain cross-checked against independent exact-rational oracles,
  including degenerate and randomized configurations.
- General boundary recovery covers the non-star/reflex twisted prism through the
  exact conforming-Delaunay path, while every returned mesh is independently checked
  for input-facet conformity and manifold topology.
- The literal enclosure/coax fixture is reconstructed natively with four filled
  material volumes and tagged boundary groups. The gmsh-impossible case is also
  solver-loadable and solved; those proofs live in [`ASCENT.md`](ASCENT.md).
- All 22 HFSS guide geometry classes have native, valid, watertight, conforming mesh
  regressions in `test/integration/hfss_cases_test.jl`. Re-running the remaining literal cases as
  full-wave ASCENT studies is external solver work, not a missing Tessella feature.

## Provenance

The final hardening commits on `main` were pushed individually, including:

- `5ec9b4e` exact Delaunay completeness;
- `abbfbe1` top-level PLC conformity;
- `9776a51` safe surface diagnostics;
- `bc9f294` manifold topology and PLC fills;
- `a630ab4` 2-D kernel/refinement contracts;
- `9b3b0ac` curve/surface bounds;
- `0b1ec2e` geometry/P2/optimization/I/O contracts;
- `5d3e466` volume meshing, partition recovery, bounded SoS;
- `d4b0dc0` industrial closure hardening;
- `27010e3` bounded flip-optimizer passes and final source marker cleanup.

Historical closure claims above are tied to their dated commands and measurements.
Current P1/P2 claims are limited to the checked source, tests, differential, aggregate
validation, and memory evidence above; the explicit P1/P2 non-claims remain open.
