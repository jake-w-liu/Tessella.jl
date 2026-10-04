# Tessella.jl status

ASCENT integration and the external HFSS full-wave campaign are tracked separately in
[`ASCENT.md`](ASCENT.md).

## Status

The original Stage 0–6 simplex-mesher roadmap is complete and remains the regression
baseline. The active roadmap now targets independent Gmsh 4.15.2 parity and is **not
complete**. Work is ordered by ASCENT meshing value before UI and post-processing.
The active package and verification contract supports Julia 1.12.x and 1.13.x,
with byte-identical mesh output verified across both lines. Older
runtime measurements retained in the dated history below are provenance, not current
support or test requirements.

| Stage | Capability | State |
|---|---|---|
| 0 | exact predicates, mesh types, validation, MSH/STL I/O | DONE |
| 1 | 2-D Delaunay, constrained Delaunay, quality refinement | DONE |
| 2 | size fields, graded curves, planar/cylinder/parametric surfaces | DONE |
| 3 | 3-D Delaunay, exact-coordinate kernel, volume fill, conforming recovery | DONE |
| 4 | uniform sizing, quality metrics, flips, smoothing, sliver reduction | DONE |
| 5 | healing diagnostics, native primitives, analytical CAD, imprints, mesh CSG | DONE |
| 6 | globally certified P2 segments, triangles, and tetrahedra plus solver-consumable I/O | DONE |

### Active parity increment

Arbitrary-length recombined B4 strips and the native GEO Line source-sampling
fix are implemented and verified in `C:/tmp/tessella_nonew_b4_strip`, from verified
and pushed `9806d5b`. Source/plan/emission and the shared boundary/projection
route are integrated. Initial normal112 checks pass24 actual source-to-volume
projection assertions and all twelve saved B4 primary/P2 differentials,
including84 terminal centers and384 retained empty stored Surface UVs.
The preexisting M7 Progression4/count8 source-coordinate error is corrected
at the unchanged2e-11 tolerance; measured maximum source error is5.55e-17.
Final focused geometry checks pass 19,516 assertions on each supported runtime,
and the full strict replay passes all 206 cases. Certified CRC runs pass
24,121 checks each for 24 matching new rows; API lifecycle tests pass 93,523
each. The combined pins now contain 256 NoNew and nine API01 rows, with
40 AddVerts rows tracked separately. Exactly two API01 source
sampling fingerprints change after independently measured Progression/Bump
correction; the other seven API01 rows and prior NoNew rows remain pinned.
The V1 471-input release candidate is rejected after two allocation growth
failures in projection/API and an invalid fixed-ulp source-area harness
assertion. Independent exact arithmetic proves the represented M4000 source
area is exactly one; the harness now checks exact dyadic cancellation.
Full resource checks otherwise pass 3,062,935 assertions; the strict
V1 replay is green with unchanged inputs/index. The V1 package run was
intentionally stopped after the resource failure, before Julia1.13 launched;
it is interrupted evidence, not a package pass. Detached projection capacity
proof passes61 checks with every actual payload/metadata field unchanged;
N500/1000/2000 allocation falls to22.89/46.00/91.27MB and both unchanged
growth bounds pass. API capacity proof passes27 actual-path checks and32
overlapping/coincident-part checks. Promoted fixes pass86 actual API assertions
plus3 provenance guards each on112/113, with exact semantics and unchanged
growth bounds. Actual generic projection AST audits pass23 each over11 bodies.
The471-input V2 strict replay passes206 cases, and both B4 resource gates pass
3,062,938 assertions with matching actual geometry/P2 records and unchanged
growth bounds. V2 is rejected after the package run exposes five stale
Quad-patch atomic-failure assertions: their recombined five-Quad input is now
supported. Its exact471 raw inputs are archived; the package112 run was
deliberately stopped after those failures, with final guard passing and113
not launched. The disproven negative fixtures are corrected without deleting
any atomic assertions. The corrected five GEO suites pass10,000 checks each on
normal bounds-checked112/113; all83 quad-patch assertion tokens and free M5
blockers remain intact. Final V3 freezes471 inputs with new scoped index and
passes both package gates,206-case strict replay and all14 resource gates.
Bounds-checked package suites pass1,311,533 checks each on Julia1.12.7 and
1.13.1, with wrapper times1884.313713s and1503.420155s. All471 inputs and
the scoped index remain unchanged before and after every release gate.
The14 allocation gates pass101,383,766 checks with exact cross-runtime resource
records and unchanged2.15 growth/+65536-byte bounds. Final release evidence
`test/tmp/b4_final_release_evidence_v3.json` has SHA256
B6F204A611C0FBFDFB62ADDFBC296C2FD733C3F53E6A3A206930A1A3850364A4.
All18 old artifact files and immutable existing-corner tables remain pinned,
apart from the two independently justified API01 source fingerprints above.
Free B4 propagation is the next implementation increment. Its physical phase
and final emission drafts pass whole-product checks for nonterminal retained
centers; public Scope/Projection/P2/API integration remains unfinished.
All other parity requirements remain part of the full active goal.

### Latest verified increment (`9806d5b`)

The dynamic rectangular Quad-grid NoNew path is implemented in
`C:/tmp/tessella_nonew_rectangular_grid` on `codex/nonew-rectangular-grid`,
from verified and pushed `78af21d`. Actual regular source incidence, convex
cells and a Jordan boundary certify four native Line chains before strictly
ordered axis-normal columns are accepted. B3 corners, adjacent-B2 edge cells
and B0 interior cells use one coupled physical face plan, indexed emission
and classified projection without an added centroid. The same actual-input
certificate admits skew, trapezoid, rounded and graded straight-edge disks.
Focused normal bounds checks pass 37,773 geometry assertions on each runtime
and 162,172 API assertions on each runtime. The permanent independent planner
suite passes 25,471 on each runtime. Independent
actual P2 Bernstein bounds and exact integrals cover Tet10, Hex27, Prism18 and
rational Pyramid14 maps. Refinement, support ownership, sparse tags, MSH
round trips, coincident independent regions and atomic precision rejection
are covered. Both runtimes produce identical 48 new CRC rows, bringing the
total to 241: 232 NoNew and nine API01. The old 193 rows and both tables remain
unchanged. The V2 full strict replay passes all 194 cases with 459 frozen
inputs unchanged. Both rectangular resource gates pass 8,483,845 assertions,
36 matching geometry/support records and all unchanged allocation-growth
bounds. Actual source-entry F1000/2000/F4000 allocation is approximately
8.6/17.2/34.4 MB after fixing boxed curve-parameter matching, down from the
rejected candidate's 27/87/302 MB. Both runtimes also pass 6,340 independent
actual source replay checks. All twelve frozen resource gates pass 95,257,884
assertions with matching cross-runtime records. Full normal bounds-checked
package tests pass 1,089,298 assertions each on Julia1.12.7 and1.13.1,
in30m05.1s and25m08.6s; all459 frozen inputs and index guards pass.
Upstream pointer-sensitive products are checked for consistent global
rank and admissibility; native source-column ordinals are not claimed to match
raw upstream pointer order. The full mesher/parity goal remains active.

### Previous increment (`78af21d`)

The bounded four-Quad strip is implemented in
`C:/tmp/tessella_nonew_four_quad_strip` on `codex/nonew-four-quad-strip`, from
verified and pushed `0f9601f`. Its ten actual boundary nodes, thirteen source
edges, 81-state coupled relation and four actual recombined terminal centers
are wired into direct indexed output and classified projection. Geometry
passes 5,607 assertions on each supported runtime, including 192 finished
products and all 256 cyclic source-catalog frames. The one-line
`Val(12)` ranking initializer removes a confirmed fixed 13 MB allocation;
detached exact-selection checks pass 629,957 assertions. The independent full
production comparison passes 8,398,080 transitions plus 6,400 scalar rank checks
over separate 256-frame and 1,024-preference domains, not their Cartesian product;
permanent short paths pass
2,429 checks under all 24 source-cell orders. Full strict native replay passes
170 scoped cases and all 24 new CRC rows match across runtimes. The new
resource gate passes 14,958,057 checks on each supported line with identical
corresponding geometry and six P2 audit rows per runtime, unchanged allocation bounds and no
true `Core.Box` in all 41 actual methods. The 446-input freeze and scoped index
are stable. All ten resource gates pass 78,290,194 assertions and 240 measured
P1 rows with matching corresponding geometry and P2 audits. The permanent
new API suite has 133,284 checks; focused runs include nine retained-artifact
checks and pass 133,293 on each line. Normal bounds package verification
passes 862,550/862,550 on each supported runtime, in 30m53.1s/22m32.1s
(wrappers 1,858.0566338s/1,357.1686771s). All 446 frozen inputs and the scoped
index remain unchanged. The 24 new CRC rows bring the total to 193: 184 NoNew
plus nine API01. Both template tables and all prior 169 rows remain byte-identical.
The next bounded implementation is an axis-normal dynamic rectangular Quad-grid
path for B3, adjacent-B2 and B0 cells. Sixteen primary P1/P2 controls and 1,469
source-only checks per supported runtime retain actual coordinates, Curve
parameters and incidence at unchanged `2e-11`; no native general-grid volume
parity is claimed. Pointer minima are not public-tag minima. B1/B4 categories,
shared neighbors/regions and transformed products remain separate pending phases.
The broad parity goal remains active.

### Previous increment (`0f9601f`)

The three-Quad `QuadTriNoNewVerts` boundary strip is implemented and verified
on `codex/nonew-three-quad-strip` from released `e3f18c9`. Its eight-node actual
source catalog, coupled 27-state cap solver, strict source/product certificates,
direct indexed emission and classified boundary planning are wired together.
Normal bounds checks independently verify all 233,280 transition
choices over 64 cyclic frames and 256 exterior preferences against the saved
immutable factory export. Exhaustive short paths pass 762 assertions on each
supported runtime, the new geometry suite passes 1,580, and the new API suite
passes 75,709. Sixteen saved native P1/P2 products and 24 new CRC records
are retained; all 136 prior NoNew and nine API01 records remain byte-identical.
The common centroid overflow fix passes 1,996 permanent regression assertions
on each runtime with healthy output bits and zero warmed allocation preserved.
Final package gates pass 721,230 assertions on each supported runtime; all
eight resource gates, strict 146-case native replay and the shared 32-case
AddVerts differential pass on unchanged 438 frozen inputs and scoped index.
Its verified and pushed release is `0f9601f`. The four-Quad release is verified
above; general source-grid volumes remain unfinished.

| Track | State | Verified implementation increment |
|---|---|---|
| P1 | **IN PROGRESS** | Native scalar/anisotropic catalog, strict `.geo` field graph with injected model/view context, Gmsh-style 1-D policy, model-level `.geo` `Mesh 0`/`Mesh 1` grading with stored `curve_params` discretizations feeding 2-D/3-D boundaries, multi-entity `Mesh 2`/`Mesh 3`/`mesh_dim` generation merged on bitwise coordinates with per-node entity ownership, and field/entity-aware 2-D, surface, and 3-D refinement, plus `.pos`-backed `PostView` scalar/vector/tensor size fields with dominant-component selection, tensor-to-metric `PostViewAnisoField` evaluation, and a documented multi-surface discrete `AutomaticMeshSizeField` analogue — sphere-fit curvature, facing-triangle `nPointsPerGap` local feature size, `hBulk` fallback, and edge-`gradation`/`smoothing` resolved from model surfaces in `.geo` and API-session background-field contexts |
| P2 | **IN PROGRESS** | 125 fixed-node Gmsh types with canonical family/order lookup and detached property metadata plus ten serializable cut/border/child/sub-element records, mixed blocks/entities/classification/dimension-0:3 periodic and embedded-curve metadata, structural validation/CRC, ASCII/binary MSH v2.2/v4.1 read/write with cumulative repeated-node/periodic sections and persistent MSH2 elementary ownership, verbatim ancillary/unknown-section and tag-remapped view-data preservation, structural `$PartitionedEntities`/`$GhostElements` metadata on `MixedMesh`, 4- and 8-byte binary `size_t` decoding plus `size_t_bytes=4` Tessella-only MSH4 binary output, classified surface/explicit-shell/embedded-volume model-to-mixed projection, owned entity names, visibility/color state, attributes, finite Point-coordinate updates, atomic live-reference retagging, dependency-safe recursive removal, explicit topology, spatial, type, plane-property, and nonpartition metadata queries, and native Point/Line/Circle/Ellipse-arc/Spline/BSpline/Bezier/Nurbs/Plane evaluation and surface reparametrization, plus exact tensor-Bernstein minimum-Jacobian certification for P2 quadrangle, hexahedron and prism, and rational collapsed-map Pyramid14 certificates for new API products |
| P3 | **IN PROGRESS** | Native analytical surfaces/imprints, classified ISO-10303-21 STEP/IGES box/sphere/cylinder/cone/ring-torus and closed planar polyhedral shell import, STEP/IGES NURBS curve and surface import with IGES export, expression-, numeric-list-, and tracked-tag-allocator-backed Point/Line/Circle/Ellipse/Spline/BSpline/Bezier/Nurbs/Loop/Plane Surface/Surface/Ruled Surface/Surface Loop/Volume with checked `SetMaxTag`, positive Point `MeshSize`, explicit-topology `PointsOf`, topology-derived Physical groups, global automatic Physical tags, owned operation-time Boolean operands with complete Delete cleanup, N-way multi-operand BooleanDifference/Union/Intersection/Fragments with OCC membership-cell decomposition and preserve-numbering tag rebinding, Box/Cylinder/Sphere/Cone/Torus/Boolean `.geo` solids materializing their Gmsh 4.15.2 OCC boundary layouts (Cylinder/Sphere/Cone behind retained compact encodings for the analytic mesher), Translate/Dilate/90°-Rotate and straight/curved curve or planar-surface periodic `.geo` execution, mesh Boolean CSG, and finalized-mesh affine transforms |
| P4 | **IN PROGRESS** | Greedy and Edmonds-blossom surface recombination with optional full-quad, Point/Line-In-Surface embeddings, Point/Line/Surface-In-Volume recovery with nested constraints and holed planar sheets, explicit planar shell/cavity volumes, holed plane surfaces, piecewise-linear planar Point-size propagation, uniform refinement, Progression/Bump/Beta curve laws and HWall variants on straight and curved (arc/spline/NURBS/OCC) edges via Gmsh's `F_Transfinite` cell-size density integrated over the native parameter, `Mesh.FlexibleTransfinite` count scaling by `Mesh.CharacteristicLengthFactor`/`Mesh.MeshSizeFactor` with the recombined-boundary odd-count rule, planar triangle/quad transfinite patches including recombined three-sided layouts, five-face transfinite prisms — affine or warped/curved boundary face grids mapped onto Gmsh's degenerate-hexahedron slot layout (`s3≡s0`, `s7≡s4`) with collapsed-triangle face fills under `Mesh.TransfiniteTri = 0`, or the compact equal-side triangular lattice under `Mesh.TransfiniteTri = 1` (`transfinite3`: diagonal-expanded slots welded bitwise onto diagonal vertices, distinct diagonal-plane interior evaluations kept unmerged, `SIM_7`–`SIM_12` templates, and GEdgeLoop-style unsigned chaining that canonicalizes surface and volume boundary inputs identically) — with bitwise shared-edge/corner certification and `transfiniteHex` interior interpolation — plus recombined five-face emission via `mesh_transfinite_prism`'s `recombine=` mask (collapsed: wedge prisms + hexahedra or `PRISM_1`/`PRISM_2` pairs; compact: `PRISM_4`/`PRISM_3` pairs; every cell certified against the unrecombined partition; invalid masks fail with Gmsh's wrong-recombination diagnostic) — and affine or face-interpolated warped six-face transfinite volumes, recombined hexahedra, model-level `Recombine`/`Mesh.RecombineAll` plumbing end to end (unstructured `recombine_triangles` post-pass with protected embedded-curve edges, recombined patch-kernel dispatch for transfinite triangle/quad surfaces, boundary-flag-derived `recombine=` masks for five- and six-face transfinite volumes, mixed `Mesh`/`MixedMesh` entity merge bucketing blocks per (MSH type, entity dim, tag) with bitwise coordinate dedup, `execute_geo`/`Save`/`geo_entity_mesh` union threading with `write_mixed_msh` output, `model_to_mixed` classified projection of recombined surface AND volume parts (tet/hex/prism blocks), PLC quadrangle→triangle boundary folding for unstructured volumes, periodic slave copies and attribute passes on mixed parts, and compact-prism interior behind-diagonal slots kept as evaluated orphan nodes matching Gmsh's written node set bitwise), prismatic 3-D layers with certified remaining-core fill/cavity walls, 2-D quad/fan layers, general-affine periodic node-pair certification/snapping, persistent native straight/curved curve relations for boundary or embedded curves with reusable masters and chains or cycles, synchronized planar periodic boundary surfaces on explicit volumes under translation and general affine (rotation) transforms with surface-derived curve masters matching upstream `setMeshMaster`, stored mesh-inert periodic volume relations as a Tessella extension (Gmsh's `setPeriodic` silently ignores dimension 3), expression/list-backed `.geo` periodic entities and transforms with orientation-only curve forms and edge-counterpart surface maps, and classified surface/volume projection with MSH2 cell ownership and supported MSH4 periodic/embedding metadata, plus joined multi-region 3-D boundary-layer fans stitched along shared edges and multi-arc vertices with certified core fill, and warped four-sided transfinite patches on ruled surfaces via 3-D Coons interpolation with exact boundary-simplicity, orientation, and fold audits, and curved-boundary planar surface meshing — native-frame forced/periodic parameter frames, sampled circle/ellipse/spline/NURBS/OCC boundary and embedded-curve chains (including single-vertex full-circle loops), ownership-masked pinch audits, and bitwise periodic affine copies — plus curved `Curve In Volume` and nested `Curve In Surface`-in-Volume embedding: seeded stored 1-D discretizations, per-link straight-chord tetrahedral edge-chain recovery, closed/self-overlapping loops, and chain-based mixed projection, closed native-curve three-segment grading floor, and `Min*`/`Minimum*` mesh-option alias synchronization, and `.geo` `Extrude … Layers` structured meshing — translational, rotational, and twist sweeps of source surface meshes through per-group level parameters with pure `Extrude(u,·)` transform-evaluation corners (upstream `pos.find` semantics), triangle generatrices emitting prisms, recombined quadrilateral generatrices emitting hexahedra, non-recombined triangles subdividing prisms to tetrahedra through the global phase-1/2/3 shared-diagonal selection with lateral-surface remeshing against the shared edge set, lateral quad/tri strips and verbatim top-copy surface meshes welded bitwise at merge, certified collapsed-vertex cells, `setAllVolumesPositive` orientation normalization on tet/hex/prism/pyramid cells, native `QuadTriAddVerts` transitions with optional `RecombLaterals`, and native `QuadTriNoNewVerts` isolated triangle/prism and quadrangle/cap-chain transitions and bounded conforming two-Tri, two-/three-/four-Quad strips and 2-by-2 four-Quad source grids with normalized positive layers, complete boundary/projection planning, actual axis-normal product certificates, and bounded exact cell-hull separation |
| P5–P6 | **IN PROGRESS** | Synchronized model/mesh API with detached cache, session-independent fixed element type/property, bounded fixed-family quadrature and actual- and explicit-order nodal reference functions, atomic whole-cache uniform refinement, affine transformation, and clearing, detached bulk/connectivity-derived data and automatic/manual global edge/triangular/quadrangular-face catalogs, hierarchical H1 bases at orders 1:15 over Point, Line, Triangle, Tetrahedron, Quadrangle, Hexahedron, and Prism families and hierarchical H(curl) bases at orders 0:11 on Line/Triangle/Tetrahedron and 0:10 on Quadrangle/Hexahedron/Prism, lexicographic orientations, and vertex/edge/face/bubble keys, plus robust cached simplex point-location, local-coordinate, forward-map/Jacobian, and element-quality queries, deterministic topology/spatial/type/plane-property/nonpartition queries, Point/Line/Circle/Ellipse-arc/Spline/BSpline/Bezier/Nurbs/Plane evaluation and surface reparametrization, owned visibility/color/attribute state, finite Point-coordinate updates, entity-name/tag/removal lifecycle, Physical-group queries, Point `set_size`, owned Boolean deletion, deterministic contiguous-block task partitioning for detached bulk/connectivity-derived, Jacobian, orientation, and element-quality queries, and periodic-map ownership, non-destructive bounded CLI with periodic/embedded surfaces, embedded volumes, and periodic explicit-shell metadata output, validated headless GUI, owned scalar nodal views, synchronized in-process plugins, plus expression- and numeric-list-backed geometry/entity lists, explicit model-topology, entity-identity/removal, spatial-query, native-metadata, native-evaluation, presentation-state, cached-refinement/affine-transform lifecycle, fixed element type/property, fixed-family quadrature, actual- and explicit-order nodal lookup, bulk/derived mesh-data, automatic/manual global edge/face topology, all-family hierarchical basis/orientation/key queries, point-location, Jacobian/reference-map, and element-quality checks, spatial and explicit-topology Point mesh sizes, topology-derived Physical groups, global automatic Physical tags, tracked tag allocators and `SetMaxTag`, t1-square, t4-hole, classified Point/Line-In-Surface, nested and holed Surface-In-Volume, and explicit Surface Loop/Volume MSH lifecycles, native/projected single-/two-direction, embedded, reusable-master/chained, and expression/list-backed periodic checks, planar periodic explicit-volume boundaries, low-level translation/rotation-periodic checks, 2-D boundary-layer quad, API-box, OCC-cylinder/cone, IGES-128 bilinear, Boolean snapshot/Delete, whole-mesh affine, element-catalog, mesh-query, mesh-entity-topology, mesh-point-location, mesh-Jacobian, mesh-quadrature, mesh-function-space, mesh-element-quality, MSH-section-preservation, MSH-partition-metadata, and MSH-size_t-width Gmsh 4.15.2 differentials, plus discrete-entity storage with `add_discrete_entity`, tag-addressed `add_nodes`/`add_elements`/`add_elements_by_type` records, model-file `import_stl`, mesh-record `create_topology`, dihedral-angle `classify_surfaces`, chord/PCA `create_geometry` with discrete evaluation, GF(2) simplicial `compute_homology` chain generators in new physical groups, element-wise `compute_cross_field` on a session view store, the `mesh.field` submodule (MathEval/Distance/Threshold composition, background and boundary-layer fields, removal semantics) consumed by 2-D/3-D generators, generator consumption of per-entity meshing attributes (transfinite curves/surfaces/volumes, recombine flags, algorithms, smoothing steps, order-2 high-order overlay, reverse and outward orientation, compound entity merging, parametric-point sizes, boundary point-size propagation, and size callbacks), record-based `recombine`/`split_quadrangles` quadrangle round trips, `get_periodic_keys` function-space pairing, entity-scoped `optimize` with Gmsh method names, multi-entity `generate`, Gmsh-parity empty-cache query semantics, `.pos` list-format view read/write including order-2 `X2` records and two-/four-matrix `INTERPOLATION_SCHEME` bindings with exact `PostViewField` evaluation (scalar/vector/tensor components, multiple time steps, curved geometry maps by intrinsic-coordinate Newton inversion, and closest-node fallback), and entity-pair `model_distance`/`get_distance` checks pinned to Gmsh 4.15.2 `occ.getDistance` oracle values |

P1 does not claim octree-identical parity with Gmsh's HXT/p4est
`AutomaticMeshSizeField` internals (the native field remains a documented
closest-vertex discrete analogue, though curvature, `nPointsPerGap`
local-feature-size, `hBulk`, and `gradation`/`smoothing` stages are
implemented and wired into `.geo` and API-session model contexts),
materially warped
quadrangles beyond transfinite ruled-surface patches, direct tensor or
metric-meshing parity, or full
`.geo`/CAD-model execution. P2 does
not claim general mixed-element generation or recombination beyond P4's first-order
surface pairing, integration of
mixed blocks into the simplex meshing kernels, curved high-order
Jacobian certification beyond second-order segments, triangles, tetrahedra,
quadrangles, hexahedra, prisms, and pyramids, internal
indices beyond `Int32`, or lossless multi-physical-group MSH v2.2 projection.
MINI basis-selector tags 138/139 serialize as fixed-width link-free records for
Tessella-only output; pinned Gmsh 4.15.2 has no mesh-record case for them.
Ancillary/unknown MSH sections round trip verbatim and parsed view sections
(`$NodeData`/`$ElementData`/`$ElementNodeData`, including Gmsh's
connectivity-implied dialect) remap tag columns into the output tag space;
explicit-tag `$ElementNodeData` that the connectivity record cannot express is
an explicit blocker under `gmsh_compatible=true`, and non-text binary payloads
are an explicit blocker for ASCII output. Some
registered fixed tags and polygon-border type 69 require explicit Tessella-only output
because Gmsh 4.15.2 cannot consume them safely. MSH2 ASCII preserves variable records
and parent/domain links; binary MSH2 and MSH4 have explicitly narrower special-record
contracts. Pinned Gmsh 4.15.2 corrupts distinct parent links in its own binary MSH2
rewrite, and nonzero-physical special MSH4 requires compatible node/entity
classification metadata for a safe rewrite. P3 does not yet claim a general OpenCASCADE BREP kernel, NURBS CAD of
unclassified topology, transformations of arbitrary CAD entities, or full `.geo` execution.
Classified STEP/IGES solids that are axis-aligned blocks, spheres, right circular
cylinders, right circular cones, or ring tori are imported and filled, as are
closed planar polyhedral shells
(`MANIFOLD_SOLID_BREP`/`FACETED_BREP`/`BREP_WITH_VOIDS` bounded by planar
`FACE`/`FACE_SURFACE`/`ADVANCED_FACE` records, with void shells subtracted);
STEP B-spline and
IGES 126/128 NURBS import as native curves/surfaces with IGES NURBS export;
other topology is an explicit blocker. Bounded
`.geo` execution covers Point/Line/Circle/Ellipse/Spline/BSpline/Bezier/Nurbs/Loop/Plane Surface/Surface/Ruled Surface/Surface Loop/Volume,
Box/Cylinder/Sphere/Cone/Torus,
BooleanDifference/Union/Intersection/Fragments of those solids — operand
groups take multi-entity `Volume{a,b,…}` lists with per-group `Delete;`, N-way
operations decompose into OCC-style membership cells, result pieces preserve or
rebind operand tags under Gmsh 4.15.2 preserve-numbering, and `v[]`/`v()`
capture the result list — Translate of remaining
native solids, Dilate, and coordinate-axis π/2 rotations of native primitives,
Point/Line-In-Surface embeddings, and Point/Line/Surface-In-Volume recovery. Its
scanner and executor handle finite arithmetic constants, pure numeric functions,
comparison/logical/ternary operators, prior scalar bindings, bounded numeric list
assignment/indexing/selection/mutation,
explicit field/physical tags, and finite `start:end[:increment]` ranges — literal or
variable-backed — in recognized numeric field
and entity lists. `If`/`ElseIf`/`Else`/`EndIf` and `For name In
{start:end[:increment]}`/`EndFor` match the built-in kernel (bit-exact against the
pinned Gmsh's entity state over 4 differential cases), with a bounded
`While`/`EndWhile` extension Gmsh lacks. Executed geometry parameters, tags, and
numeric entity memberships
use those same bounded semantics. Translational `Extrude {dx,dy,dz} {..}`
and rotational `Extrude {{axis}, {point}, angle} {..}`
run for points, curves, and planar surfaces as both a statement and a
side-effecting value term, reproducing Gmsh's tag allocation, output lists,
signed-generatrix topology, and post-extrusion merge behavior bit-for-bit
(26 differential cases), with `Layers`/`Recombine`/`ScaleLastLayer`/`QuadTri*`/
`Using` parameters stored per created entity and
`Geometry.ExtrudeReturnLateralEntities` support. `Circle`/`Ellipse`
records reproduce Gmsh's `EndCurve` control-point layout, `Plane{..}`
fallback normal, and `sys2x2` solve bit-for-bit, with `Surface`/`Ruled
Surface` filling typed triangular (3 generatrices) or ruled (4) and
`In Sphere`/`Using Point` sphere-center metadata (8 differential cases
covering entity tags, point coordinates, type strings, bit-exact arc
evaluations, and sorted loop boundaries). `Cylinder`/`Sphere`/`Cone`
materialize their Gmsh 4.15.2 OCC boundary representations — rim/pole
Points, closed-circle, seam, and degenerate edges, `Cylinder`/`Sphere`/`Cone`
and `Plane` faces, signed surface loops — behind the retained compact
encodings that still drive native volume tessellation (9 OCC differential
cases covering entity tags, point coordinates, type strings, signed
boundaries, parametrization bounds, curve evaluations, and bounding boxes).
`Spline`/`BSpline`/`Bezier`/`Nurbs` records reproduce Gmsh's built-in-kernel
control-point lists bit-for-bit — Catmull–Rom with extrapolated or cyclic
ghost endpoints, the piecewise `mat` matrices with right-extremity
reflection, full-list De Casteljau, and float32-rounded knot vectors with
inferred degree — including empty `Knots {}` falling back to the BSpline
record, evaluated-but-ignored `Order`, and `[0,1]` bounds versus the raw
Nurbs knot interval (19 differential cases covering entity tags, point
coordinates, type strings, endpoint wiring, parametrization bounds,
evaluations, derivatives, second derivatives, curvature, bounding boxes,
`parFromPoint`, `getClosestPoint`, `isInside`, and oriented surface
boundaries, all bit-for-bit). Curve-loop tags use their own
allocator namespace and chain-sort by endpoint connectivity. Arc and
spline-family metadata follows transforms, Duplicata, coherence merges,
retagging, and removal;
meshing paths that cannot honor the geometry fail explicitly instead of
treating curved records as chords. `Function`/`Macro name ... Return`/`Call name;`
run Gmsh's zero-argument function semantics with token-level body capture, shared
variable scope, and bounded recursion; named number/string options read and write
through the vendored option table. The meshing-constraint statements reproduce the `GEO_Internals`
setters — `Transfinite Curve` (Progression/Power/Bump/Beta, the `_HWall`
wall-height variants, and the grammar-only `Beta_Symmetrical` forms, with
signed-tag distribution reversal and tag-0 wildcards), `Transfinite Surface`
with `Left`/`Right`/`Alternate*` arrangements and live-target corner
validation, `Transfinite Volume` with 6/8-corner acceptance, `TransfQuadTri`
boundary-diagonal centroid transitions for six-face and collapsed five-face
transfinite volumes, `Recombine`/`Smoother`/`MeshAlgorithm`/`MeshSizeFromBoundary`,
`Reverse`/`ReverseMesh`, `Degenerated`, and `Compound` with the
`MeshAlgorithm` suffix marker — while `RelocateMesh`, `ReorientMesh`, and
`RecombineMesh` validate but store nothing, matching Gmsh's no-mesh state.
`OptimizeMesh "method"` mirrors `GModel::optimizeMesh(how)`: the name is
validated first (unknown optimizers error even without a mesh), `""`/
`"Gmsh"`/`"Optimize"`/`"Relocate3D"` run the boundary-preserving tetrahedral
optimizer when the cache has tets, `"Laplace2D"`/`"Relocate2D"` run Laplacian
smoothing on the triangle cache, and meshes lacking the relevant cells are
silent no-ops.
`Split Curve{c} Point{...}` mirrors `GEO_Internals::splitCurve`: a
`Line`/`Spline`/`BSpline` record breaks at the listed control points into
same-type segments tagged through `NEWCURVE()` — `NEWREG()` under the
default `Geometry.OldNewReg`, so physical-group, loop, and higher-dimension
maxima all feed the counter (verified bit-for-bit against the pinned
binary's 10/11, 21/22, and 6/7 tag choices) — while curve loops rewire
`±c` to the new tags in forward or reversed-negated order, physical
memberships move in the raw group records (a later `-= {c}` touches
nothing), a closed spline reseeds into a single closed curve, and the
original record is deleted. Break vertices off the control list silently
yield one renumbered identical curve; an unknown curve or a non-splittable
type reports `Msg::Error` plus `Could not split curve` and continues. The
deprecated `Split Curve(c) {...}` form parses with its warning, the
`Duplicata { Split ... }` nested-transform form is legal inside shape
lists, and inside an `FExpr_Multi` (`x =`/`x() = Split ...`) the
production's own `tEND` consumption makes the assignment a syntax error
after the split runs — the error lands on the next statement's token and
`error tEND` recovery drops that statement, matching the pinned binary's
statement-poisoning behavior (a non-`;` tail like `, 3` prevents the
reduction, so the split never runs).
The discrete-model statements mirror `Gmsh.y`'s forms:
`Homology`/`Cohomology`/`Betti` queue `addHomologyRequest` requests (bare
`0:3` dimensions, `{dom}`/`{{dom},{sub}}` `ListOfDouble` lists, and the
`(dims){dom,sub}` form that requires both lists), executed at the end of
`Mesh n`/`mesh_dim` generation like `GModel::computeHomology` — an empty
domain covers the model's top-dimensional entities, explicit dims filter to
`0:getDim()`, `Betti` reports ranks without storing, and surviving
generators become discrete entities under `H_k{dom[,sub]}i`/`H^k{…}i`
physical groups matching the pinned binary's names and tag allocation
(requests repeated for a computed space store nothing, requests queued
after meshing never run, and upstream's `.msh` chain-only output quirk is
deliberately not mirrored — the generated mesh is preserved).
`CreateTopology;`/`CreateTopology{a,b}` (simply-connected repair plus
GEO-internals export defaults — the export rebuilds the pending physical
registry from entity memberships like `exportDiscreteGEOInternals`),
`ClassifySurfaces{a,b,c[,d]}`, and
`CreateGeometry;`/`CreateGeometry{shapes}` dispatch to the discrete-entity
implementations. `.geo` `Merge "file.msh"` now imports the file like
`GModel::readMSH`: elementary (MSH2, missing tags → entity 0) or entity
(MSH4) cell classification becomes discrete entities carrying their cells,
MSH4 node/element tags, per-node entity classification, parametric
coordinates, and declared boundaries are preserved, physical memberships and
names join the model, and every element family is accepted while only the
simplex cells fold into the mid-file mesh (meshing discrete entities
themselves remains an explicit non-claim). Discrete tags share the
elementary namespace like upstream's `getMaxElementaryNumber` — `newp`/
`newl`/`news`/`newv`-family reads after `Merge`, `CreateTopology`, or
homology output skip them (bit-exact counters against the pinned binary).
`Delete{…}`/`Recursive Delete{…}` follow `GEO_Internals::remove` exactly:
per-entity list-order attempts, absolute matching for points and curves,
signed matching for surfaces and volumes, boundary-ownership refusal,
pre-collected recursive boundaries, and counter decrement only at the
dimension maximum. Physical groups keep stale member integers that
resurrect on tag re-creation, and queries filter them out like
`GModel::getPhysicalGroups`. The named `Delete` forms cover `Embedded`,
`All`, `Model`, `Physicals`, `Variables`, `Options`, `Meshes`, `Struct`,
`Field[i]`, variables, and `name~{expr}` namespaces. `SetTag` fails
explicitly like Gmsh's mid-parse model retag (28 differential cases over
entity inventory, physical state, mesh node sets, and error paths). It rejects
boundary-layer and pipe (`Using Wire`) extrusion forms,
`Fillet`/`Chamfer`, allocator reads after
untracked topology-changing declarations (tracked Boolean operand `Delete` and
`SetMaxTag` counters stay live, while entity-list `Delete`/`Recursive Delete`
joins the untracked set), and geometry-derived Physical
right-hand sides beyond the documented inline topology queries.
`SetMaxTag Point|Curve|Surface|Volume` follows the active factory: Built-in can set
or lower a geometric counter, while OpenCASCADE only raises it. Allocator reads use
the greatest counter among activated factories; later primitive allocation still
accounts for occupied hidden topology. `Box` and Cylinder/Sphere/Cone/Torus
boundary entities are materialized in the native model like Gmsh's OCC
layouts, and Boolean results materialize their result boundary into explicit
Points, Curves, Curve Loops, Surfaces, Surface Loops, and shell-grouped
Volumes for every supported section pair — disjoint solids materialize as
component volumes, cavities attach as inner shells, a geometrically empty
result binds nothing, and unsupported sections (oblique plane×cylinder,
non-axial plane×cone, off-axis cylinder×sphere, cylinder×cylinder) raise
explicit errors.
`MeshSize` and `Characteristic Length` store finite positive constraints on existing
explicit Points selected by `:`, bounded expressions/ranges, or whole and selected
numeric-list variables. Inline `PointsOf` blocks select recursive boundary Points of
explicit Point, Curve/Line, Surface, and Volume entities; signed tags are normalized,
hole boundaries are included, and embeddings are excluded. Direct updates are atomic,
and the session API invalidates its cached mesh only after success. Planar surface
refinement extends Point constraints piecewise-linearly over the deterministic initial
constrained triangulation, including linearly sized generated straight-curve
subdivision nodes and an exact constant-size path for uniform constraints. Boolean result volumes resolve `PointsOf` through their materialized
boundary like any explicit volume. Exact Gmsh
mesh topology, mesh-size
selectors other than inline `PointsOf`, nonpositive values, and Gmsh's silent
missing-Point behavior are explicit non-claims. Physical declarations accept an
explicit literal tag (zero and negatives included, matching Gmsh's raw `(int)`
cast), with an optional name, or a nonempty name with an automatic tag from the
global Physical namespace. An already-bound name resolves to its existing tag.
Compound `+=`/`-=` modify memberships (`+=` appends without dedup; `-=` deletes
an emptied group and is a silent no-op on a missing group), `*=`/`/=` and `+=` on
a missing group record errors, and duplicate `=` declarations stay recoverable
diagnostics.
Physical Point accepts inline `PointsOf`; Physical Point/Curve/Surface accept inline
`Boundary` and `CombinedBoundary` over Curve/Line, Surface, and explicit Volume
entities, respectively. `Boundary` collects immediate boundaries before group
membership is deduplicated; `CombinedBoundary` keeps tags with odd multiplicity. Hole
and cavity boundaries participate, while embeddings do not; Boolean volumes
answer boundary queries through their materialized result shells. Empty
combined boundaries and unsupported dimensions
are explicit blockers.
The model and session APIs return detached, sorted group, membership,
reverse-membership, and name-query results. Names are dimension-scoped; removing a
name or selected/all groups leaves geometry intact and keeps automatic Physical tags
monotonic.
They also enumerate explicit entities, report model dimension, and return direct or
recursive boundaries and direct adjacencies with deterministic Gmsh-compatible
ordering, orientation, and combined-incidence cancellation. `is_entity_orphan`
reports downward-closure connectivity to the highest-dimension entities,
excluding embeddings, matching Gmsh 4.15.2. Queries preserve the
session mesh cache and exclude embeddings. Box and Cylinder/Sphere/Cone/Torus volumes expose their materialized OCC boundary
topology, and Boolean volumes expose their materialized result boundary the
same way.
Exact bounding boxes cover explicit straight-edge topology, analytical native
primitives, Boolean operation-time result snapshots, and the union over a nonempty
model. Finite containment queries select complete entity boxes and ignore embeddings
when bounding their target; both direct and session queries are read-only. Tessella
does not add OpenCASCADE shape-tolerance padding (`1e-7` in the pinned fixtures),
rejects nonfinite coordinates and invalid filter dimensions.
Native metadata classifies visible entities as `Point`, `Line`, `Plane`,
or `Volume` plus the curved `Circle`/`Ellipse`/`Surface`/`Cylinder`/`Sphere`/
`Cone`/`Unknown` OCC types and `Nurb` for the whole built-in spline family;
Boolean result boundaries classify like any
materialized topology. `get_type` is a compatibility synonym
for `get_entity_type`. Plane-property queries return detached unit-normal
coefficients `[a,b,c,d]` for `a*x+b*y+c*z=d`, oriented by the exterior loop; other
visible native types have empty property vectors. Because `GeoModel` does not own
partition entities, existing entities report parent `(-1,-1)`, empty partition
membership, and a model partition count of zero. Direct and session metadata queries
are read-only and preserve the mesh cache; imported `MixedMesh` files carry
structural `partition_data` (partitioned-entity records, ghost entities, and
ghost-element ownership) round-tripped through MSH 4.1 ASCII and binary output.
Native Point, straight-Line, and explicit-Plane evaluation provides values,
derivatives, zero curvature, Plane normals and principal directions, parametrization
and bounds, containment counts, and closest-point projection. Lines use `[0,1]`;
Planes use the deterministic Gmsh-compatible orthonormal frame derived from the
exterior loop. Built-in Spline/BSpline/Bezier/Nurbs curves evaluate through
bit-exact ports of Gmsh's `InterpolateCurve` family — Catmull–Rom with
extrapolated or cyclic ghost endpoints, the piecewise UBS matrices with
right-extremity reflection, full-list De Casteljau, and `findSpan`/
`basisFuns` over the float32-rounded knot vector — and differentiate by
Gmsh's 1e-8 finite differences for derivatives, second derivatives, and
curvature. Their parameter bounds are `[0,1]` except Nurbs, which reports
its raw first/last knot interval; bounding boxes use Gmsh's 10-sample scan,
parametrization runs the multi-seed damped-Newton `XYZToU`, closest point
runs the 100-sample scan plus golden-section search, and physical
containment chains `XYZToU(relax=1)` with the `parBounds` check — all
matching the pinned binary bit-for-bit.
Physical Plane containment uses exact coplanarity and the strict
trimmed interior, while parametric containment uses inclusive rectangular bounds.
Line projection clamps to the segment and Plane projection is untrimmed. The direct
and session APIs are read-only, preserve the mesh cache, and reject unsupported
implicit, degenerate, nonfinite, malformed, or unrepresentable cases.
Point and straight-Line parameters can be reparametrized on any explicit Plane,
including sources outside that Plane's topology or geometric support. The operation
uses orthogonal Plane parametrization and preserves the session cache. Its `which`
selector is validated but does not alter a native nonperiodic Plane result.
Visibility defaults to `1` and RGBA color to `(0,0,255,0)`. Recursive presentation
setters traverse explicit boundary topology down to Points; state follows retagging,
is cleaned with removed entities, and does not invalidate the session mesh. Global
model attributes own detached NUL-free string vectors and return names in lexical
order. Finite Point-coordinate updates preserve tag-owned metadata and invalidate a
session mesh only after success. Dependent native Line and Plane queries immediately
use the new coordinates; measured Gmsh 4.15.2 Plane parameter bounds remain stale
after the equivalent update. Per-window visibility is stored display state, Box
and Cylinder/Sphere/Cone/Torus volumes present their materialized boundary
entities, and Boolean volumes present their materialized result boundary the
same way.
`API.mesh.refine` atomically refines linear-simplex caches and supported native
mixed P1/P2 cells through the pinned family templates, including four pyramids
and eight tetrahedra per pyramid. Mixed refinement returns linear children;
supported legacy straight-simplex refinement retains its re-elevated query
overlay. Returned storage is detached. Rejected resource bounds or unsupported
curved native CAD placement
leave the prior cache unchanged. `API.mesh.clear` discards
the complete cache and preserves model geometry. Entity-selective clearing removes
the cells classified on the listed entities and drops their owned nodes;
entities owning no cache cells are no-ops. Tagged Point15 cells can be cleared;
boundary entities absent from a legacy tet-only cache remain no-ops.
`API.mesh.affine_transform` accepts a native 4×4 matrix or exactly 12/16 Gmsh
row-major entries, rejects singular and nonfinite maps, and commits only a detached,
validated whole-cache result. Reflections rewind simplex connectivity instead of
leaving inverted cells. The operation does not rewrite model geometry or periodic
relations; the index-aligned classification snapshot is retained through the
transform. Entity-selective transforms move only the nodes classified on the
listed entities and rewind only cells whose every node moved; Gmsh 4.15.2's
per-entity semantics are mirrored on the flat shared-node cache, and outputs
failing validation leave the cache unchanged.
Bulk session queries expose detached node tags and coordinates, element types, tags,
and connectivity, type and whole-dimension filters, and maximum tags. Node and
element tags are dense identifiers for legacy caches and stored public labels
for tagged caches, including sparse tags. Query families use numeric MSH type
order independently of internal block positions. A classification snapshot
stored with the cache supports entity-filtered node, element, type-funnel,
Jacobian, orientation,
key, and `get_element` queries; `include_boundary` emits transitive boundary-entity
nodes after the entity's own, `dim=-1` ignores `tag`, and unknown entities fail
explicitly. Native geometry queries compute available curve or surface
parameters, including ruled and three-sided fills. Tagged raw records preserve
their stored parameters; Gmsh can leave new face-center parameters unstored.
Points, Volumes, and
all-dimension queries emit none, while `get_nodes_by_element_type`
packs each repeated node's parameters on its owning entity in entry order,
matching Gmsh 4.15.2's variable-width emission.
Connectivity-derived queries expose repeated per-element node coordinates,
barycenters, and edge/face nodes in Gmsh's local ordering. Native mixed P1/P2
queries use actual cells. Full quadratic edge/face outputs include interpolation
nodes; `primary=true` selects corners.
Nonfinite fast coordinate sums fail explicitly. Nondefault
`task`/`num_tasks` returns the contiguous Gmsh block slice for `get_elements_by_type`,
barycenters, and edge/face nodes; `task>=num_tasks` is the silently-empty range,
and negative, zero-count, or non-integer task arguments fail explicitly.
Reference quadrature covers all fixed-node Point, Line, Triangle, Quadrangle,
Tetrahedron, Hexahedron, Prism, and Pyramid types without requiring session state.
Detached coordinate triples and weights preserve Gmsh 4.15.2's available economical
rules, including Triangle tables through order 20 and Tetrahedron tables through
order 21. Higher `GaussN` orders use Gmsh's tensor transitions, Prism composes the
matching Triangle and Line rules, and `CompositeGaussN` selects tensor rules
directly. Native Gauss--Legendre, Duffy, and Gauss--Jacobi generation is limited
to 128 points per axis and one million output points. Pinned Gmsh defines no
Trihedron integration rule.
Actual-order Lagrange/isoparametric queries cover every fixed-node Point, Line,
Triangle, Quadrangle, Tetrahedron, Hexahedron, Prism, and Pyramid type. Explicit
numeric names select complete family bases across the catalogued orders.
Hierarchical H1 functions and gradients cover orders 1 through 15 on Line,
Triangle, Tetrahedron, Quadrangle, Hexahedron, and Prism families plus the
order-independent Point basis; hierarchical H(curl) functions and curls cover
orders 0 through 11 on Line, Triangle, and Tetrahedron and orders 0 through 10
on Quadrangle, Hexahedron, and Prism. They preserve Gmsh's orientation-major
layout and the lexicographic orientation rank of primary public node tags.
Node-based keys use public node tags;
bubble keys use public element tags, and edge/face keys use topology identifiers.
Explicit-order queries that require nodes absent from the linear cache fail
instead of inventing keys. Hierarchical key catalogs lay out vertex, edge, face, and bubble keys like
Gmsh 4.15.2's `getKeys`; edge- and face-based keys reuse the global topology
catalogs and lazily add only
requested edges or faces. Key coordinates and all returned arrays are detached. Higher-order
Prism and incomplete
Pyramid paths that Gmsh 4.15.2 cannot construct are native, invariant-certified
extensions. Pyramid/Trihedron hierarchical spaces remain explicit blockers. `get_basis_functions_orientation` accepts nondefault
`task`/`num_tasks` and returns the contiguous Gmsh block slice; `task>=num_tasks`
is empty, where the pinned release segfaults in its unguarded per-entity loop.
Whole-cache edge and face creation assigns positive global identifiers to missing
primary cell topology in first-encounter cache order; legacy simplex caches visit
segments, triangles, then tetrahedra.
`add_edges` and `add_faces` atomically attach explicit positive identifiers to node
pairs, triangles, or quadrangles; face identifiers share one namespace across both
face types. Exact association repeats are idempotent. Later creation preserves
manual entries, begins automatic candidates at the relevant catalog size plus one,
skips identifiers already in use, and fills missing primary edges and faces.
Edge lookup returns the same tag in either direction and orientation `+1` for
ascending node tags or `-1` for descending tags. Both face types preserve Gmsh
4.15.2's zero orientation result. All-entity results are tag-sorted instead of
exposing Gmsh's hash iteration; their tag-to-node maps match. Query results are
detached, and every cache replacement invalidates both catalogs. Tessella rejects
zero or conflicting identifiers, repeated or unknown nodes, and malformed or partly
invalid batches atomically; Gmsh 4.15.2 accepts or partially applies those cases.
Entity-selective creation adds only the cells classified on the listed
entities, matching Gmsh 4.15.2's dimTags selection.
Cached point-location queries use deterministic candidate lookup, decreasing dimension
then increasing public tags, scaled affine inversion, and exact-rational
fallbacks. Native mixed P1/P2 cells use their actual isoparametric maps, with
conservative quadratic candidate bounds; native mixed reference maps and Jacobians
also use the actual P1/P2 cells. The strict contract uses Gmsh 4.15.2's default
`1e-6` reference tolerance; relaxed search widens it by decades through `1.0`.
Segment and triangle off-span
coordinates use orthogonal projection with unused coordinates fixed at zero, avoiding
the pinned implementation's measured inversion artifacts. Every cache replacement or
invalidation discards the hierarchy. `get_element` resolves type, public
connectivity, and entity ownership through the stored classification snapshot.
Cached quality queries preserve public-tag request order and implement the 13
documented Gmsh 4.15.2 measures for linear triangles and tetrahedra. Compatible
segment measures follow the pinned API; its undefined or unreliable 1-D Jacobian,
inverse-gradient-error, and isotropy paths are explicit blockers. Other native
mixed P1/P2 families support `minEdge` and `maxEdge`; remaining quality measures
have precise blockers. Scaled arithmetic
keeps the warmed per-element kernel allocation-free, while exact-rational and fixed
256-bit BigFloat fallbacks distinguish degeneracy and retain finite ratios across
overflowing spans, subnormal scales, and ill-conditioned simplices. Batch queries
allocate detached request-normalization and result arrays and do not mutate the mesh
cache. `get_element_qualities` accepts nondefault `task`/`num_tasks` over the
requested-tag positions and returns the contiguous Gmsh block slice with
slice-scoped validation; `task>=num_tasks` is empty, where the pinned release
reads its request vector out of bounds. Cached `get_jacobians` partitions the same
way over cached type-block positions.
The immutable element catalog owns family/order-to-type lookup and detached property
metadata for all 125 fixed types. Its session-independent API accepts canonical family
names case-insensitively and follows Gmsh's complete-type fallback for unavailable
serendipity types. Verified high-order-prism and trihedron layouts remain available
where Gmsh 4.15.2's property call fails. Compatible special types expose the same
zero-node metadata as Gmsh; special records without nodal properties fail explicitly.
Entity names belong only to existing positive-tag entities and need not be unique.
Atomic retagging moves topology, Point sizes, Physical memberships, embedding sources
and targets, periodic relations, primitive encodings, Boolean-result snapshots, the
entity name, visibility, and color while preserving monotonic allocation. Boolean
operand tags remain
historical provenance. The bounded contract accepts positive `Int32` entity tags in
dimensions 0 through 3; unlike Gmsh 4.15.2, it rejects invalid dimensions and
nonpositive retag targets, does not preload names for missing entities, and migrates
names with their entities.
Ordered entity removal skips surviving boundary and embedding dependencies. Recursive
removal follows explicit boundaries down to Points but leaves embedded entities. It
validates every input before committing, cleans names, visibility, colors, Physical
memberships and empty groups, target embeddings, affected periodic relations,
primitive encodings, Boolean-result snapshots, and newly dangling construction loops,
and keeps allocation monotonic. Box and Cylinder/Sphere/Cone/Torus boundaries
recurse through the materialized shell, and Boolean result boundaries recurse
through their materialized shells the same way — the operation-time operand
snapshots and component links are cleaned with the volume. Tessella owns the native mutation rather than exposing Gmsh's separate
model/CAD synchronization layers, and it does not retain Gmsh's independent names or
stale periodic-master state after deletion.
Boolean volumes own operation-time operand geometry, and API or `.geo` operand
deletion makes volume tags reusable without changing an existing result.

P4's uniform-refinement slice applies the exact Gmsh 4.15.2 linear segment, triangle,
and tetrahedron child templates while sharing edge midpoints, compacting unused nodes,
preserving parent tags, and supporting atomic session-cache refinement. Its
straight-curve slice covers normalized affine-line
Progression/Power, Bump, and Beta parameters plus all three HWall variants. Its surface
transfinite slice covers
already-discretized, count-matched, three- and four-sided planar chains using Gmsh's
specific triangular and average-chord Coons interpolation. Four-sided ruled
surfaces (`Surface`/`Ruled Surface` fillings) with non-coplanar boundaries mesh
as warped transfinite patches — the same Coons interpolation evaluated in 3-D
with an exact orient3 boundary-simplicity audit, per-triangle nonzero-area
certification, an area-weighted orientation check against the ring's Newell
normal, and a fold audit on every shared grid edge. A `Plane Surface` whose
transfinite boundary is not coplanar follows upstream `computeMeanPlane`
semantics instead of rejecting: the declared plane comes from the first
non-collinear on-curve boundary samples (curve control points never enter),
the side chains project onto it for (u,v) bookkeeping and interior
interpolation — so the interior stays exactly planar — while emitted
boundary nodes keep their true positions and the warped-patch audit runs on
the emitted band. `Surface … In Sphere` (and four concentric arc
generatrices, auto-detected like `ruledSurface::checkSphere`) evaluate the
ruled surface's own `S(u,v)` — the `TransfiniteQua` blend projected onto the
sphere at radius |S0−O| — for interior nodes on flat and non-coplanar
boundaries alike, including inside transfinite-volume face grids.
Three-sided ruled surfaces (`Surface` fills, non-coplanar or spherical
boundaries included) interpolate their interior in real space per
`meshGFaceTransfinite`: each `TRAN_TRI` Cartesian point is inverted through
`GFace::XYZtoUV`'s loose off-surface Newton (`Precision = 1e-3`,
`MaxIter = 10`, fixed 9×9 restart grid, silent last-iterate fallback) and
re-evaluated as `point(Up,Vp)` on the ruled parametrization —
`TransfiniteTriB` plus `TransfiniteSph` — in both the collapsed and
compact `TransfiniteTri` kernels. Four-sided grids can also
be emitted as first-order Gmsh type-3 quadrangles with exact projected
corner-Jacobian certification. Eight-corner transfinite volume blocks require
all six boundary surfaces transfinite (Gmsh's incompatible-surface gate) and
mesh them with the patch kernel, reindex each grid into its canonical slot
through the eight dihedral permutations, and interpolate the interior with
Gmsh's `transfiniteHex` Coons volume — six face interpolants minus twelve
edge interpolants plus the trilinear corner term — parameterized by
chord-length ratios along the s0s1/s1s2/s1s5 edge chains. Boundary nodes
reuse the face grids bitwise and the boundary is the canonical conforming
split the six-tet cell subdivision induces, audited strictly outward; warped
and otherwise non-affine blocks are admitted whenever their face grids are
meshable. The kernel's direct `faces=nothing` path still certifies an affine
eight-corner parallelepiped; canonical triangular prisms
use its legacy collapsed-grid five-face tetrahedral path under
`Mesh.TransfiniteTri = 0`, lifted to warped
and curved boundaries through the same meshed face grids on the degenerate
hexahedral slot map (`s3≡s0`, `s7≡s4`), or the compact `transfinite3`
subdivision under `Mesh.TransfiniteTri = 1` — compact equal-side triangular
face lattices expanded into the square slots with upper-triangle slots
welded bitwise onto diagonal vertices, the `SIM_7`–`SIM_12` cell templates,
and unsigned GEdgeLoop-style boundary chaining shared between the standalone
surface meshes and the volume face reads. `mesh_transfinite_prism`'s
`recombine=` mask emits Gmsh's five-face recombined cells on both layouts
as a `MixedMesh` — collapsed: wedge prisms plus interior hexahedra
(`CREATE_HEX`) for all-recombined masks or `CREATE_PRISM_1`/`CREATE_PRISM_2`
pairs for axial-only recombination; compact: `CREATE_PRISM_4` diagonal cells
plus `CREATE_PRISM_3`/`CREATE_PRISM_4` strict-lower pairs with Gmsh 4.15.2's
observed strict ordering — each cell's tetrahedron decomposition certified
against the unrecombined partition, boundary sheets emitted as type-3
quadrangles or outward triangles, and unsupported partial masks rejected
with Gmsh's "Wrong surface recombination in transfinite volume" diagnostic.
`mesh_transfinite_volume`'s `recombine=` mask covers Gmsh's full six-face
decision tree as a `MixedMesh`: all-recombined emits `CREATE_HEX`
hexahedra with quadrangle boundary sheets, a single unrecombined opposite
face pair emits the corresponding prism pair per cell (v/u-spanning pairs
in the orientation-fixed `MPrism` ordering, or the w-spanning
`CREATE_PRISM_1`/`CREATE_PRISM_2` pair), and every other partial mask
rejects with the same diagnostic — each cell carries a certified
shadow-tetrahedron decomposition and the emitted boundary audits exact
shadow-exterior coverage, on both the affine and `transfiniteHex`
warped-face paths.
Surface recombination now
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
slave, reusable masters, and master/slave chains or cycles. They synchronize boundary
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
forest through MSH2/MSH4. Standalone periodic surface pairs — not bounding a volume —
project through the multi-surface `model_to_mixed` overload, which merges the
per-surface meshes on shared nodes, emits each shared entity's cells once, and
serializes the surface map, its induced boundary links, and explicit curve
relations spanning the projected surfaces.
Gmsh 4.15.2 does not serialize the Point/Line/Surface-In-Volume relation. Nested
Point/Line-In-Surface constraints are certified against each sheet's face complex;
embedded sheets may contain interior loops. MSH4 retains the nested curve relation.
The CLI uses these projections for periodic
or embedded `-2` output and classified `-3` output. Volume boundaries may
carry curved edges whose arcs leave an adjacent `Plane Surface`'s declared
plane; the boundary-fold audit certifies each emitted boundary triangle by
the incident tet's apex and rejects only a true edge-through-triangle
pierce — Gmsh 4.15.2 silently emits a self-intersecting mesh on a strongly
inward-bulging boundary, which Tessella refuses to emit. P4 does not yet claim
non-affine CAD curve integration or size-map curve laws,
quasi-transfinite or holed transfinite patches,
general CAD parameterizations,
or corner-reordered fills, volume/hybrid
recombination beyond the five-face transfinite
prism (whose `recombine=` mask covers both collapsed
and compact layouts), selective or high-order refinement, coarsening,
3-D multi-wall boundary-layer fans, curved
periodic surfaces, or allocator reads after
untracked topology-changing declarations. Cyclic periodic relations store and
serialize like upstream (whose deferred copy starves them); Tessella instead
synchronizes cyclic curves through the converged parameter fixpoint and rejects
a cyclic surface dependency explicitly at `mesh_model_surface` rather than
recursing or returning an empty mesh.
Expression/list-backed `Periodic Line`, `Periodic Curve`, and
`Periodic Surface` `Translate`, `Rotate`, and `Affine` statements —
`Periodic Volume` is a `syntax error (Volume)` like upstream — plus the
transform-free orientation-only curve form, signed entity tags, and the
`Periodic Surface j {curves} = k {curves}` edge-counterpart form that derives
its transform from mapped boundary vertices,
including bounded ranges and numeric list variables in entity sets, are in
scope for the bounded `.geo`
executor.

General OpenCASCADE/unclassified NURBS CAD, remaining algorithms/fields, broad
formats and API, GUI, and post-processing are unfinished parity tracks, not
project non-goals.

## Verification history (newest first)

2026-10-04 — Four-quadrangle NoNew strip, verified:

- The bounded native TF5-by-TF2 source retains ten actual boundary nodes,
  four Quad4 cells and three shared faces. Its 81-state joined cap relation,
  strict source/product certificates, direct indexed output and classified
  projection are integrated. Recombined laterals retain four actual terminal
  centers and emit eight Tet4, twenty Pyr5 and `4(N-1)` Hex8 cells. Source
  node/cell order and actual rounded five-/two-node Curve chains are retained.
- The unchanged factory export independently verifies 8,398,080 production
  transitions and 6,400 scalar rank checks. All 256 cyclic-frame and 1,024
  canonical exterior-preference domains are separate, not a Cartesian product.
  Permanent short paths pass 2,429 checks under all 24 source-cell orders.
  Geometry passes 5,607 checks per supported runtime, covering 192 actual
  finished products and 256 catalog-only cyclic frames.
- Twenty-four saved Gmsh 4.15.2 P1/P2 products certify actual interpolation
  maps, typed boundaries, ownership and support identity. The permanent API
  suite passes 133,284 checks per runtime; focused wrappers add nine unchanged
  prior-artifact checks, yielding 133,293. Their times are 294.7669043s and
  269.4085843s, with all 119 input paths stable. The new 24-row CRC artifact is
  byte-identical across runtimes, SHA256
  `2BAEAD51AF5BF3CDE76DF3A26F4E5529DA79BF7B1D0A19505D8A40B30E022D7F`.
  The total is 193 rows: 184 NoNew plus nine API01. All prior 169 rows and both
  immutable 315-/13-record table sources remain byte-identical.
- The twelve-field ranking initializer uses `Val(12)`, removing a confirmed
  fixed 13 MB allocation without changing selected records, costs or tie
  behavior. Detached exact-selection checks pass 629,957 assertions.
  Full strict native replay passes all 170 scoped cases in 312.551s, log SHA256
  `161B29164056DE826F28E31AB86282213902D0CCD26E1825DD6004821A32CB00`.
- All ten resource gates pass 78,290,194 assertions and 240 P1 measured rows
  across the supported runtimes, preserving corresponding geometry and P2
  audits. The new four-Quad gate passes 14,958,057 assertions per runtime,
  including six P2 rows each and 41 actual lowered method bodies with no true
  `Core.Box`. Keyword implementation bodies and boxed positive controls are
  included. The unchanged allocation-growth bound remains mandatory.
- Both final normal bounds package gates pass 862,550/862,550 assertions.
  Julia 1.12.7 takes 30m53.1s (wrapper 1,858.0566338s), log
  `test/tmp/pkg_four_strip_frozen_v1_julia112.log`, SHA256
  `AB63845A963E291C99A2518C9B65EA2BE07C3D7CE5BAB7995059B0B962B50740`.
  Julia 1.13.1 takes 22m32.1s (wrapper 1,357.1686771s), log
  `test/tmp/pkg_four_strip_frozen_v1_julia113.log`, SHA256
  `287F4DCA4B1815E4CE36015D28BDABFBD1268CF8AE574DF1C60395FBE3BCED98`.
  The package count is exactly `721230 + 5607 + 2429 + 133284 = 862550`;
  the nine focused artifact checks are not additional permanent tests.
- Every final gate preserves all 446 frozen paths: 98 production/Project,
  194 test and 154 validation files, plus the scoped index and base HEAD.
  Freeze manifest SHA256:
  `0FF77953A027C73D24F0591AACF9EFAF0212544088DAAFE3348DCFE5BBD81624`.
  Scoped index SHA256:
  `c2614390e877718e1e0744ca397c7dc6967f2b1752140abacf3d478758141678`.
- Next-grid preparation retains sixteen primary 2-by-3/3-by-3 P1/P2 controls
  for N1/N3, both normal directions and both lateral policies. Actual source
  categories include B3, adjacent-B2 and B0. Current source-only comparisons
  pass 1,469 checks per runtime in 152.9688209s/123.7865852s, with sixteen
  cyclic Quad and complete cell-array matches and maximum coordinate/raw
  parameter discrepancy `2.0594637106796654e-12` at unchanged `2e-11`.
  All 117 comparison inputs and the 446 release guards remain unchanged.
  There is no native general-grid volume-parity claim. Pointer minima are not
  public-tag minima; B1/B4 categories, neighbors and transformed grids remain
  separate pending phases.

2026-10-04 — Three-quadrangle NoNew strip and finite centroid overflow fix, verified:

- The native recombined TF4-by-TF2 source retains eight actual boundary nodes,
  three convex Quad4 cells and two opposite shared edges. Original source-cell
  order and rounded native four-/two-node Line chains remain authoritative.
  Ten-edge incidence, exact winding, nonadjacent edge contacts and end-cell
  containment certify the complete source partition. Finite actual stored
  columns certify strictly ordered axis-normal translation in either direction.
- Free laterals use a constant-width 27-state cap solve over the unchanged
  315-record corner relation. Exterior preference changes dominate tentative
  cap alignment; original cell/factory/state order resolves ties. Both shared
  faces are chosen together. Direct indexed output retains `8(N+1)` primary
  nodes and at most `18N` cells. Recombined laterals retain `3(N-1)` Hex8,
  introduce exactly three strict actual terminal centroids and emit six Tet4
  plus fifteen Pyr5. P2 support graphs contain `42N+21` or `42N+45` nodes.
- An independent saved factory export verifies every production transition:
  all 64 cyclic local frames and 256 exterior preference assignments, totaling
  233,280 joined choices. Selected factory records, costs, alignment and both
  shared-face reversals agree. The complete output payload SHA256 is
  `6ab4333a337c2cf349c77363f9d99a7a87426dcffc39c46f9cc3f3353fa9f72e`.
  The permanent short-path suite passes 762/762 on each supported runtime for
  one through four intervals under all six original source-cell permutations.
- Sixteen saved Gmsh 4.15.2 P1/P2 products comprise eight unit fixtures and
  eight independent geometric variants.
  Each order has 1,584 independent actual element-map certificates. Typed
  macro boundaries, opposite internal faces, support identity and classified
  ownership are verified on raw captures before promotion. Actual warped Pri6
  maps use their interpolation maps. Eight native provenance gaps containing
  168 empty upstream UV records are explicit, separately counted non-claims.
  The saved artifact SHA256 is
  `c8c53c5192475adf5006bc812441726ae0fdd5fd32c3d4b899d956b4bbe941e8`.
- The new geometry suite passes 1,580/1,580 on Julia 1.12.7 and 1.13.1,
  including 552 actual local-frame/source-order emissions, graded and signed
  products, corruption and atomic unsupported-category checks. The existing
  two-Quad geometry/control focus passes 911/911. Previous TC4 unsupported
  fixtures advance to TC5 because the three-Quad category is now supported.
- Final new API suites pass 75,709/75,709 on each line. Their bodies take
  3m48.8s and 2m58.1s, with guarded wrappers taking 277.723s and 222.148s.
  All 114 direct input paths and hashes agree before/after and across runtimes.
  Coverage includes lower-dimensional ownership and parameters, P2 supports,
  lifecycle, edits, refinement, file round trips and atomic failures.
  Both suites assert all 24 new CRC rows. Independent artifact generation
  passes 16,831 checks per runtime and produces byte-identical complete text.
  The new CRC artifact SHA256 is
  `1fcf6abebb6de28c9d323f98cf7fb77a52a243f89d35f9594bde4548a00b2966`.
  All 136 previous NoNew rows, nine API01 rows and both template sources remain
  byte-identical; no older artifact is repinned in this increment.
- A valid public large-coordinate witness confirmed a shared centroid bug:
  source `X=1e308` and cap `X=nextfloat(1e308,1000)` have finite positive delta
  `1.9958403095347198e295` and strict finite mean
  `1.0000000000000998e308`, but raw accumulation previously produced `Inf`.
  Build ordinary CAD topology first, move points through the public setter and
  bind the matching retained translation. Direct huge CAD declarations collapse
  at construction and are retained failed fixtures, not this witness. Free
  products succeeded while recombined products rejected a valid actual center.
  Exact full-map and fan certificates independently confirm the corrected mean.
- The common helper retains its original Float64 addition order and fixed-column
  omission. Only nonfinite sums enter a noinline `Rational{BigInt}` fallback;
  finite center components retain their original bits, and nonfinite inputs
  retain downstream rejection. No global precision or production dependency
  changes. The permanent suite passes 1,996/1,996 on both lines, including all
  16 public single-/three-Quad large-coordinate products, exact maps/volume,
  healthy and collapsed six-/eight-corner parity, cancellation, zero warmed
  1k/2k/4k allocation and real `Core.Box` checks with positive controls.
  The focused guarded gates pass 1,997 including their input guard, in
  166.535s/144.216s. No successful huge-coordinate Gmsh parity claim is made.
- Retained development failures were resolved without weakening certificates:
  a probe addressed `.sweep.mesh` instead of `.sweep.volume`; the nonbinary
  test helper used literal `.6` instead of the actual represented layer;
  early API tests used a wrong option function and assumed unavailable bare
  P2 classification. A coincident-region test counted other owners' curves as
  its vertical curves; final tests select the actual four lateral generators.
  An OnlyEmpty P2 bit change is independently confirmed upstream behavior:
  generation first strips high order globally and reconstructs native CAD
  supports afterwards. Exact P1 carriers, P2 parameters and CAD reevaluation
  remain asserted. The prior failed logs are retained under `test/tmp`.
- All eight final resource gates pass, totaling 48,374,080 assertions. Each
  runtime has 96 measured P1 rows; corresponding rows agree in geometry and nodes
  across runtimes and the four paths. Six new-three-Quad P2 audit rows per
  runtime additionally agree in actual supports, owners and digests.
  Each job preserves its 100 direct inputs and all 438 frozen paths/index/HEAD.
  The allocation bound remains `next <= 2.15*previous + 65,536 bytes`.

  | Gate | Assertions per runtime | Body 1.12 / 1.13 | Max doubling 1.12 / 1.13 | Actual methods, no `Core.Box` |
  |---|---:|---|---|---:|
  | Two-Tri | 1,596,393 | 3m13.0s / 2m35.9s | 2.133982296 / 2.133828030 | 17 |
  | 2-by-2 patch | 6,385,557 | 3m49.8s / 3m17.7s | 2.110915215 / 2.111182551 | 29 |
  | Two-Quad strip | 4,846,359 | 3m48.9s / 3m09.0s | 2.105788917 / 2.106187223 | 36 |
  | Three-Quad strip | 11,358,731 | 4m53.3s / 4m04.3s | 2.081875534 / 2.081750514 | 41 |

  Warm ordinary centroid calls still allocate zero bytes. The three-strip scan
  includes both mean helpers among its 41 actual methods and detects both
  deliberately boxed controls. The resource summary retains every warmed
  allocation/time row and all TOML/log/wrapper hashes under
  `test/tmp/three_strip_final_resource_summary_two_tri_quad_patch_quad_strip_three_quad_strip.json`.

- The final implementation, tests and validation are staged/frozen at 438 paths:
  97 production/Project, 188 test and 153 validation. Scoped index SHA256 is
  `78342624ec7a0046099ff47942eac2ad4081b090b2381684ebe286c9d57dd293`.
  Every final job checks all paths, raw hashes, scoped index and base HEAD before
  and after execution. The normal bounds package gate passes
  721,230/721,230 assertions on Julia 1.12.7 in 28m38.4s; its wrapper takes
  1723.0557202s and both whole-freeze guards pass. The sequential Julia 1.13.1
  package gate also passes 721,230/721,230 in 23m01.7s, with wrapper
  1385.9378071s and identical whole-freeze guards. Final package log SHA256s are
  `9c4326cccc349b9a86b31faf93dfc3f22c414e0f7f869963f9db8202ed968d35`
  and `3feea9063640191884980d5ed6b138d6753ce37904d0dcec1312b9e2a5b1a0c8`.
  All eight final resource gates are complete. The common AddVerts primary
  differential passes all 32 existing cases in 166.814s with whole-freeze guards;
  translation, grading, rotation, twist, fixed columns and toroidal sweeps retain
  their existing typed-cell and coordinate contracts. Largest coordinate
  difference is `2.0594637106796654e-12` under the unchanged `2e-11` gate.
  Its log SHA256 is
  `8424d63a25c51593d4ac0814e5f460e37fcf5aa26db3f48409b9705fe36d68cd`.
  Full strict native replay passes all 146 scoped cases in
  286.631s, including the original 86 and all saved 12 two-Tri, 16 pivot-patch,
  16 two-Quad and 16 three-Quad fixtures. Both whole-freeze guards pass; the
  strict log SHA256 is
  `84f2cecdc543fc5a9f9286716c855a7194d67ddc0b351fa223d51684293d456d`. The earlier
  pre-centroid resource results are development evidence, not final gates.


Bounded two-quadrangle NoNew strip and native straight grading correction,
2026-10-04:

- The final V6 freeze contains 96 production/Project, 181 test and 152
  validation paths (429 total). Its exact delta from V4 is five production
  files and two permanent tests; Project and every validation file remain
  unchanged. V5 was invalidated before any release gate because the last
  callback fixture legitimately triggered filtering; its corrected size jump
  still exercises smoothing and full callback traces. All old logs and
  manifests remain retained.
- Final normal bounds-checked `Pkg.test` passes 641,183/641,183 assertions
  on Julia 1.12.7 in 27m03.5s (wrapper 1,628.14s) and Julia 1.13.1 in
  21m38.1s (wrapper 1,302.72s), with all 429
  byte/path and runtime/index guards intact before and after. TOML resolves
  in the isolated test target. The deliberate external-process interruption
  probe emits a child-side broken pipe on both runtimes; each parent suite exits
  successfully. The two full package jobs ran sequentially. Their log SHA-256
  values are `c30a1dd3250130a730f0ae32d66811f81d854ac7040d176a5bba636a03088608`
  and `08dd7b3069b59c4778fbf5d2e52455253f5b8da45a4aa1d11992944491afa321`.
- All six final V6 1k/2k/4k resource gates pass with unchanged bounds:
  strip 4,846,359 assertions per runtime in 214.59s/168.00s; four-quad patch
  6,385,557 in 227.10s/195.48s; two-triangle 1,596,393 in 200.72s/155.32s.
  Each report has 24 measurements and 99 current frozen direct inputs.
  Geometry/node-count digests match across all four paths and both runtimes.
  Maximum doubling ratios are 2.106173, 2.111195 and 2.134001 respectively,
  within `next <= 2.15*previous + 65,536`. Actual keyword bodies and wrappers
  scan 36/29/17 methods, all with zero true boxes and positive controls.
- The final V6 strict NoNew replay passes in 281.37s: original 86 samples,
  twelve two-triangle, sixteen patch and sixteen strip saved P1/P2 products.
  Eight independent strip variants and 126 empty stored UV records retain
  their separate provenance counters. Final public GEO constraint replay
  passes 38 cases with zero documented gaps on both runtimes in
  172.66s/133.32s. Its six nonfatal Gmsh-side temporary merge-file cleanup
  warnings per runtime are retained; no native production defect is established.
- Final API generation 0/1 focused gates pass 3,523/3,523 on both normal
  bounds-checked runtimes in 272.58s/207.84s. Public Gmsh differentials pass
  9,388/9,388 in 256.38s/194.64s, each with 108 fixtures, 196 stages and the
  same one declared source-history blocker. Nine artifact records match the
  staged pins and each other, with data SHA-256
  `77242ec7d81a9c438381b24febe4fc63680481108739edee0adc0a2fe79c8b8a`.
  Every job retains all 429 byte/path and index guards before and after.
- Final affected curve/mesh1d normal focuses pass 565/565 on both runtimes
  in 187.72s/149.94s with four exact input hashes matching V6. The new
  tiny compile-min targeted policy proof passes 203 checks and its eleven actual
  wrapper/generated bodies have zero boxes. Independent ordinary replay
  verifies the eight open-Line products, eight actual closed-Circle products
  and three near-0.75 controls; a separate 21-check proof confirms cached
  smoothing, exact trapezoid order and identical entire callback traces on
  both sides of the generic primitive's rounding threshold.
- Independent staged review confirms both immutable tables and all 112 old
  NoNew CRC rows/complete old artifact blobs unchanged. Twenty-four unique
  strip rows are added; API grading has exactly three justified hash-only
  row changes and six byte-identical data rows. All six final resource reports
  independently match the V6 source hashes, path digests and allocation bounds.
  No confirmed discrepancy remains in this increment's reviewed scope.
- The native recombined source contains six actual boundary nodes and two
  strictly convex Quad4 cells with one shared edge. Four native Line chains
  have opposite two-/three-node widths. Signed source incidence, actual
  sampled chains and strictly ordered stored planes certify the axis-normal
  product in either direction, with positive normalized graded layers.
- A packed nine-state joined cap relation uses the unchanged 315 corner
  templates and joins actual opposite shared-face diagonals. Free laterals
  introduce no nodes and emit at most `12N` cells from `6(N+1)` primary nodes.
  Recombined laterals retain `2(N-1)` Hex8 and emit four Tet4 plus ten Pyr5
  around exactly two terminal centroids. Each actual eight-corner mean must
  be representable strictly inside its macro; one-ULP center failures reject
  atomically. The seven-cell fan keeps its own immutable descriptor instead
  of changing the original corner-template records.
- Direct indexed volume, cap and lateral output and actual two-/three-node
  source/top chains preserve complete boundary and CAD carriers. Actual P2
  graphs have `30N+15` nodes for free laterals and `30N+31` for recombined
  laterals. Dimension 0/1/2 owners are `8`, `8N+12`, `16N-2`; dimension 3 is
  `6N-3` or `6N+13`. Both caps have three owned nodes and fifteen closure nodes.
  Long-chain laterals have `6N-3` owned and `10N+5` closure nodes; short-chain
  laterals have `2N-1` owned and `6N+3` closure nodes.
- Focused normal bounds-checked geometry passes 911 assertions. The completed
  exhaustive joined relation audit passes 10,369 checks. Its 1k/2k/4k packed
  DP allocations are 16,238/29,182/55,182 bytes. All 25 new kernel methods have
  zero true `Core.Box`, with a deliberate boxed-local positive control.
- The permanent API suite passes 48,105/48,105 on normal bounds-checked Julia
  1.12.7 in 3m23.3s and 1.13.1 in 2m29.9s. All 103 input-byte and source-path
  guards pass. The ignored wrappers then exit 1 only while printing a final
  Windows path-key status; both failed-reporting logs are preserved and the
  next wrapper normalizes separators. Twenty-four new P1/P2 CRC data rows are
  byte-identical across both runtimes; their data SHA-256 is
  `7097832e215df543207ae84223f8c0ec1f52f263287ec6d3cdc9e4e5f6fc2b38`.
  All 112 older NoNew rows and both immutable template tables remain unchanged.
  P2/P1/P2, remapping, refinement, selective clearing, atomic precision failures
  and ASCII/binary MSH4 preserve actual carriers and node parameters. Physical
  groups retain their public model contract; the generated volume-only cache
  keeps its existing bare-cache contract, confirmed against clean `a941938`.
- The independent targeted replay passes all sixteen saved Gmsh 4.15.2 P1/P2
  products in 238.68s on normal bounds-checked Julia 1.12.7, with all 104
  discovered inputs unchanged. Whole Prism6 maps use their quadratic normal
  determinants and exact integral; Pyramid5 uses collapsed bilinear base
  determinants and sum/12 integration, including warped internal faces.
  Positive actual maps, convex containment and the oriented typed degree-one
  shell certify each macro partition. The eight recombined fixtures contain
  126 empty stored oracle UV records; computed inverse queries are checked
  separately. Original differential cases and serializers remain intact.
- Six newly saved Progression 4 sources exposed the old native Line analytic
  shortcut, with coordinate errors up to `8.60e-9`. Native nonuniform Line
  density laws now use the existing bounded adaptive primitive and numerical
  inversion from `F_Transfinite`; all eight new source variants pass the
  unchanged `2e-11` gate with maximum error `8.22e-12`. The independent
  thirteen-law capture also covers Bump, Beta, reversals and HWall. Ordinary
  laws pass `2e-11`; HWall Progression passes its established `3e-7` bound
  with maximum error `9.50e-10`. Standalone analytic and nonflexible native
  Circle contracts are unchanged. The decreasing `1e-15` endpoint now has defined
  zero terminal density when its logarithm argument rounds nonpositive;
  unresolved interior partitions still reject explicitly. The independently
  disproved old native fraction expectations were corrected while retaining
  coordinate, ownership, ordering and seeded tiny endpoint-snap checks.
- V2 1k/2k/4k resource gates pass on normal bounds-checked Julia 1.12.7 and
  1.13.1. The strip passes 4,846,359 checks per runtime in 225.04s/185.85s;
  the four-quad regression passes 6,385,557 in 242.83s/201.30s; the two-triangle
  regression passes 1,596,393 in 211.53s/173.83s. Each report contains 24
  measurements and 99 unchanged direct input hashes. Geometry digests match
  across all four paths and both runtimes. The largest allocation doubling
  is 2.106197 for the strip, 2.111191 for the patch and 2.133964 for two
  triangles; every measurement retains `next <= 2.15*previous + 65,536`.
  The V2 strict differential passes in 296.91s: all 86 original samples,
  12 two-triangle products, 16 four-quad products and 16 strip products retain
  their P1/P2 certificates. The new strip accounts for 126 empty stored
  oracle UV records, checked separately from computed inverse parameters.
  All 429 frozen inputs and runtime/index blobs match before and after.
- The corrected native Line route also disproves four old API expectations:
  Progression, Bump and Beta generation fixtures and a five-node periodic
  master. Independent Gmsh 4.15.2 captures differ from the old analytic
  fractions by up to `2.05e-7`, while the corrected output differs by at most
  `2.87e-12`. The periodic sample check keeps its original `3e-12` tolerance;
  the other samples use the established native-coordinate `2e-11` bound.
  Both normal runtimes generate identical nine-record API CRC files: exactly
  the three graded records change, with the other six byte-identical.
  Complete API generation 0/1 unit gates pass 3,523/3,523 on both runtimes
  in 247.32s/196.67s. The entire public Gmsh differential passes 9,388/9,388
  each in 241.47s/191.19s, covering 108 fixtures and 196 generation stages,
  with the same one explicit source-history blocker. All 429 V3 input hashes
  and staged/runtime blobs remain unchanged before and after each run.
  Standalone analytic helper, nonflexible native Circle, Spline and uniform
  Line contracts remain unchanged. The V3 freeze changes exactly the two API test
  files from V2; all 96 production and 152 validation paths, every actual
  resource input and the strict driver's certificate/catalog inputs are
  unchanged. This initially preserved resource and strict provenance; the
  subsequent production fixes below require fresh final runs.
  A full normal 1.12 package gate then ends after 416,427 successful checks
  with one setup error: the new saved-oracle API test imports TOML, but the
  isolated test target did not declare it. TOML is now a standard-library
  test extra; production dependencies remain unchanged. The failed run takes
  827.91s and retains all 429 V3 input/index guards. The next V4 freeze changes
  only that test-target metadata and two API grading documentation hunks;
  the validation driver's noncomment lines are identical.
  Further independent grading controls expose preexisting flexible-count
  defects: a recombined Beta 0.5 Line at twelve declared nodes and factor two
  emits seven nodes instead of Gmsh's six, and the direct API emits twelve
  because dimension 0/1 generation omitted current option mirroring. These
  fixes and their regression checks are being completed. Final production gates
  will be rerun after the source stabilizes; no failed or interrupted package
  run is a pass.
- Further pinned Gmsh 4.15.2 controls confirm and correct the related policy
  defects before the next release freeze. Flexible transfinite odd-N forcing
  depends on the integrated mass being greater than 0.75; it cannot compare
  the raw coefficient. Native GEO Line count evaluation uses Gmsh's bounded
  1e-5 derivative stencil. Seventy-eight saved controls pass 391 checks:
  seventy-six meet the established 2e-11 coordinate bound, while two
  Progression_HWall controls retain the existing 3e-7 bound (maximum error
  5.728e-10). Positive HWall types transform signed or zero wall heights before
  count policy; their native Line uniform fallback uses the original declared
  count and does not solve the law again using the emitted count.
- API generation 0/1 previously returned before mirroring its seven current
  meshing options. The options now reach the planner through an immutable tuple
  and apply after its existing deep copy. Default private planner calls retain
  their contract; source-cache classification still uses retained source
  attributes. The permanent API suite and allocation/AST controls pass
  619/619 on both normal bounds-checked runtimes. Instrumented staging makes
  exactly one deep copy and the measured 1k/2k/4k allocation difference is at
  most 272 bytes. Precision failures preserve model/cache identity, allocators,
  session options and callback. All inspected keyword bodies have zero boxes.
- Ordinary recombination now applies the odd-N rule to the existing integrated
  primitive before filtering, without a second integration or callback pass.
  Gmsh 4.15.2's `increaseN` is an identity, so algorithms 2 and 4 do not add
  two further nodes. Eight saved open-Line controls distinguish this defect
  from the previously masked API-option bug. Closed curves use N=edges+1
  before omitting the repeated endpoint: supported Cylinder Circle controls
  retain seven unique nodes under algorithm 0 and eight under algorithms
  1/2/4, matching Gmsh rather than the former seven/nine outcomes.
- Flexible Circle Progression and Progression_HWall placement formerly used
  the emitted count, giving actual coordinate errors 0.475 and 1.300. The
  scoped correction uses the original HWall count and scaled density-law count;
  both saved public products now satisfy the existing 3e-7 Circle bound.
  Five independent Beta Circle controls pass 21/21 count/parameter checks;
  maximum parameter error is 6.22e-15. The near-0.75 Circle masses differ
  from unit-Line masses and retain their actual Gmsh count decisions.
- A final three-case ordinary Line threshold probe uses fixed mesh sizes
  prevfloat(4/3), 4/3 and nextfloat(4/3). Gmsh emits two nodes and one Line in
  all three; the generic native primitive rounds above 0.75 and emits three
  nodes. The pinned Line derivative primitive rounds below 0.75. This
  correction reuses the existing sample times and smoothed mesh sizes, with
  the exact native trapezoid order and no second integration or field callback.
  All three corrected public products match the saved two-node primary outputs.
  Its permanent controls also verify varying-field single-pass/smoothing traces.
  Failed harness drafts and all earlier gate logs remain retained; final focus
  and the next complete frozen
  package, resource and strict gates are pending.
- The first frozen strict differential passes in 258.85s with all 429 inputs
  and staged/runtime blobs unchanged. The first large strip resource pass has
  4,846,351 successful checks and one recombined GEO allocation-growth failure:
  14,959,403 to 32,631,966 bytes exceeds the unchanged `2.15*previous+65,536`
  bound. Independent profiling isolates the jump to the mixed merge's
  fourfold dictionary expansion. Reserving coordinate lookup, ownership and
  coordinate-vector capacity from the largest retained part reduces complete
  GEO allocations to 14,418,694/28,274,536/56,025,976 bytes. Both growth checks
  pass with the same bounds. A 1,605-check detached comparison covers older
  grids, ordinary triangle/box, duplicate/coincident parts and split buckets;
  a normal 1,140-check whole-GEO comparison preserves all ordered fields,
  retained parts and actual curve parameters. The applied executable body
  matches that proven candidate. The concurrent first package run was stopped
  after 339s before making this production fix; its direct post-run guard
  confirms all 429 inputs and staged blobs unchanged. These logs are retained.
  An additional boxed keyword control exposed skipped implementation bodies
  in the resource scanners. All three scanners now retain wrapper checks and
  inspect generated keyword bodies, including optional positional wrappers.
  Their controls pass on both runtimes. The strip scans 36 actual methods,
  quad patch 29 and two-triangle 17, all with zero boxes. The broader independent
  source audit additionally scans 39 wrappers and nine keyword bodies with
  zero boxes. No production boxing defect was found.
  The broader independent parity objective remains active. Larger bounded
  grids, general transformed products, mixed roots, shared regions,
  copied-source chains, collapsed columns and cyclic sweeps remain unfinished.

Bounded 2-by-2 NoNew quadrangle source and projection allocation fixes,
2026-10-04:

- Final normal, bounds-checked package gates pass 591,635/591,635 assertions
  on Julia 1.12.7 in 26m13.3s and Julia 1.13.1 in 21m22.8s, including all
  38 optional pinned Gmsh source checks. Both wrappers exit successfully and
  verify unchanged 95 production/Project, 176 test and 151 validation inputs
  before and after. All 422 staged blobs match the tested runtime files.
- A native recombined TF3 source retains nine distinct actual nodes and four
  strictly convex Quad4 cells around one existing interior pivot. Its four
  straight native Line curves retain their three-node sampled chains, including
  graded chains and rounded samples that need not be exactly collinear with
  the stored CAD endpoints. The actual source complex and exact axis-normal
  translation product are certified before output, in either direction and
  with positive normalized graded layers.
- The physical pivot fixes shared and cap diagonals; constant phase selections
  preserve the native corner-fan priority. Free laterals emit `(24N-8)` Tet4
  plus four Pyr5 cells; recombined laterals emit `4(N-1)` Hex8 plus twelve Pyr5.
  Both use exactly `9(N+1)` primary nodes without a new body centroid. Direct
  indexed volume, cap and lateral builders retain ordinary constructor copies,
  actual reference-map certification and complete typed-boundary audits with
  opposite outward cycles on every internal face.
- Projection retains actual three-node source/top chains and introduces no CAD
  entity for an internal radial edge. Actual P2 supports give `50N+25` nodes,
  with dimension 0/1/2/3 owner counts `8`, `8N+20`, `24N+6`, `18N-9`.
  Normal bounds-checked focused geometry and API suites pass 1,359 and 69,842
  assertions. Source grading tests distinguish actual source Quad9 centers
  from copied-cap pivot-diagonal midpoint supports. Selective-clear checks
  preserve surviving geometry, carriers and parameters across legacy dense
  tag compaction. All 88 previous NoNew CRC rows and both immutable template
  tables remain unchanged; 24 new API P1/P2 CRC rows retain their baseline pins.
- Both normal bounds-checked resource gates pass 6,385,550 assertions and all
  24 standalone/GEO/projection/API records at 1k/2k/4k layers. Sorted typed
  3-D cell-coordinate signatures agree across paths. The unchanged growth
  bound is `next <= 2.15*previous + 65,536 bytes`; the largest observed doubling
  is 2.141130. All 98 directly used inputs and the full 422-file staged freeze
  remain unchanged. A deliberate boxed closure verifies the lowered-code
  scanner's `GlobalRef(Core, :Box)` detection, and 26 checked helpers have zero
  boxes. The earlier two-triangle resource scanner now has the same positive
  control instead of its former false-negative test. Its final regression gate
  passes 1,596,389 assertions on both runtimes, all 24 records and the unchanged
  growth bound, with zero boxes in 17 helpers and unchanged frozen inputs.
- GEO classification reserves from the actual child-map closure, avoiding
  over-reservation on meshes with many interior nodes. Exact owner and payload
  comparisons cover ordinary surfaces and a transfinite box as well as this
  patch. Projection reserves actual cell/coordinate table sizes and reuses its
  audited input edge and typed external-boundary sets instead of rebuilding
  the volume topology. All actual-cell, chain, duplicate-claim and complete
  shell checks remain. A normal bounds-checked 166-assertion comparison covers
  all four cell families, both policies, node/cell/block permutations, scoped
  projection, exact output metadata and malformed-input/atomicity parity; all
  95 production inputs stay unchanged and four affected methods have zero
  boxes. Julia 1.12.7 at 4k layers allocates 97.48/35.16 MB for
  free/recombined standalone volumes and 675.34/255.75 MB for the complete API
  entry. Recombined projection uses 181.12 MB. Allocation minima include
  construction, validation and downstream copies; extraction, independent
  auditing and digest computation run outside the measured calls.
- The full strict Gmsh 4.15.2 differential passes on normal bounds-checked
  Julia 1.12.7 in 253.59s, retaining all 86 original samples and 12 saved
  two-triangle cases and adding 16 lossless four-quad P1/P2 fixtures, including
  four graded rounded trapezoids. Eight new recombined fixtures contain 160
  empty stored oracle UV records; computed native surface inverses are tested
  separately. All 422 frozen inputs and runtime/index blobs match before and
  after. Pointer-dependent admissible native cell choices remain independent
  geometric certificates rather than production allocation-order pins.
- Independent signed-source checks confirm source/cap orientation, node and
  typed-boundary correspondence, volume families and returned GEO entity
  lists against Gmsh, including its extra source entry for a negative tag.
  All twelve direct checks pass; the mesh sign convention is unchanged and
  all 422 frozen inputs remain stable.
- The two-quadrangle strip is next. Larger/transformed grids, mixed roots,
  shared regions, copied-source chains, collapsed columns and cyclic sweeps
  remain explicit planning blockers. The broad parity and optimization goal
  is not complete.

Bounded two-triangle NoNew grid and precision fixes, 2026-10-04:

- Final normal, bounds-checked package gates pass 520,434/520,434 assertions
  on Julia 1.12.7 in 25m08.6s and Julia 1.13.1 in 20m37.5s, including the
  optional pinned Gmsh source checks. Both wrappers exit successfully and
  confirm unchanged 94 production/Project, 171 test and 150 validation inputs
  before and after; all 415 staged blobs match the tested files.
- One isolated strictly convex four-corner planar transfinite source with
  four native straight endpoint-only curves and two actual conforming Tri3
  supports axis-normal translation in either direction and positive normalized
  graded layers. A joined 13-pattern prism relation fixes the internal swept
  face once for the whole region. Free laterals emit `6N` Tet4; recombined
  laterals emit `2N` Pri6; both use exactly `4(N+1)` primary nodes. No centroid,
  CAD curve or CAD surface is invented for the internal source diagonal.
- Actual source topology, unchanged in-plane columns and strictly ordered
  actual normal planes certify the full product. Actual template partition
  and P1 map checks remain mandatory. Every internal typed face has two
  opposite outward cycles, and the finalized surfaces exactly cover the
  external typed boundary. Volume and lateral matrices use direct source-major
  indices; ordinary constructors, validation and downstream copies remain.
- Actual P2 carriers give exactly `18N+9` nodes, with dimension 0/1/2/3 owner
  counts `8`, `8N+4`, `8N-2` and `2N-1`. Each lateral owns `2N-1` nodes and
  has a `6N+3`-node closure. Both runtime suites retain all 64 older NoNew CRC
  products and verify 24 new native API P1/P2 coordinate, cell and support pins.
- The final normal, bounds-checked resource gate passes 1,596,388 assertions
  on each runtime: 24 complete standalone/GEO/projection/API records at
  1k/2k/4k layers, both lateral policies, equal sorted typed 3-D cell-coordinate
  signatures across all four paths, the unchanged bound
  `next <= 2.15*previous + 65,536 bytes` and zero `Core.Box`
  in 17 helpers. All 97 production/certificate inputs and the full 415-file
  freeze remain unchanged before and after both gates.
  The largest observed doubling is 2.101181. Compared with the initial frozen
  grid implementation, all 48 allocation observations decrease by 6.56–36.36%.
  At 4k layers, Julia 1.12.7 allocates 35.11/23.52 MB for free/recombined
  standalone volumes and 220.49/150.07 MB for the complete API entry.
  These are minima of three warmed calls for this bounded grid; they include
  constructors, validation, classification and retained downstream copies.
  The gate's additional payload extraction, auditing and digest computation
  run outside those measured calls.
- The final strict Gmsh 4.15.2 differential passes on normal bounds-checked
  Julia 1.12.7 in 240.46s. It retains all 86 original samples and adds 12
  saved two-triangle P1/P2 cases. Six recombined cases have 56 empty stored
  oracle UV records; native computed surface inverses are tested separately
  and are not represented as stored-parameter parity. All 415 frozen inputs
  and staged/runtime blobs match before and after.
- Curved classification uses local extent and converts world tolerances to
  the actual native parameter range. Straight and curved chains bound
  endpoint snapping below actual adjacent sample gaps. This fixes large-radius
  curved samples and strongly graded straight samples collapsing into an
  endpoint. Normal bounds-checked focused tests pass 438 assertions, including
  native parameters, curve/point ownership, projection, repeatability, tiny
  scales, closed curves and overlapping admission neighborhoods; both curve
  methods have zero `Core.Box` and stable input hashes.
- Retained NoNew templates use an exact logical centroid witness when the
  stored layer planes have no representable interior Float64 point. Actual
  emitted centroids still must form valid cells. Mixed signed volumes use
  anchored corner differences, reference face winding and compensated sums.
  Independent volume regressions cover Tet4/Hex8/Pri6/Pyr5, both windings,
  large offsets, sub-ULP logical means and small/large scales.
- Coherence uses the unpadded local GEO bounding-box diagonal, including
  native closed-circle bounds, with operation-relative spatial bins and exact
  zero-tolerance duplicate handling. Wide bins and bounded quotients preserve
  nearby-point candidates under tiny explicit tolerances. Factor scaling
  before three-axis norm evaluation avoids overflow. Point extrusion captures
  the source tolerance before moving its copy; padded synchronized factory
  coherence remains a separate unimplemented API context.
- Single-entity simplex refinement inherits actual parent classification and
  supports through successive refinements, matching the existing multi-entity
  route. API classification and merge tables reserve checked actual payload
  bounds. Projection reserves its known face count during the repeated
  boundary audit; the complete audit and its incidence behavior remain.
- The bounded 2-by-2 quadrangle source and two-quadrangle strip are next.
  Larger/transformed grids, mixed roots, shared regions, copied-source chains,
  collapsed columns and cyclic sweeps remain explicit planning blockers.
  The broad parity and optimization goal is not complete.

Triangular NoNew region, transfinite frames and actual P2 queries, 2026-10-04:

- Final normal-compile, bounds-checked package gates pass 490,276/490,276
  assertions on Julia 1.12.7 in 28m54.9s and Julia 1.13.1 in 23m24.1s.
  Both include 38 optional pinned-source provenance checks. All 93 production/
  Project, 165 test and 149 validation inputs remain unchanged throughout;
  every staged input matches the tested file. A final independent review
  confirms no remaining defect within this increment's declared contracts.
- `QuadTriNoNewVerts` now supports one isolated source triangle with three
  distinct boundary vertices. Free laterals emit three Tet4 cells per interval;
  recombined laterals retain one Pri6. Both use exactly `3(N+1)` nodes and
  introduce no centroid or cap diagonal. One operation-owned plan supplies
  the volume, five boundary surfaces, classification and projection through
  standalone, GEO and API entry points.
- The true six-corner prism relation contains 13 admissible patterns among
  27 face states. Its immutable packed table occupies 429 bytes; warmed table
  access and fixed template emission into pre-reserved storage allocate zero
  bytes. Independent exact
  full-map and typed-face tests pass 980 assertions on each Julia version.
  Six/eight-corner global separation shares the existing bounded exact hull
  certificate. New prism global regressions pass 80 assertions, including
  actual thin-helix overlap witnesses and separated sweeps.
- Three- and four-sided public transfinite surfaces reconcile their actual
  winding with the signed CAD boundary before applying `Reverse`. Unpinned
  four-sided frames follow unsigned edge chaining, fixing their deterministic
  diagonal convention. Spherical fills preserve the original CAD frame and
  radius under corner reordering. Internal volume grids keep their existing
  frames. The triangle differential covers 18 orientation cases and 120 cells;
  the extended quadrangle differential covers 244 public cases and 888 cells.
  All 128 straight two-cell diagonal/winding fixtures match exact oriented
  node identities; the earlier eight low-level quadrangle cases are unchanged.
- Ruled-surface inverse evaluation now serves model and node parameter queries,
  including quadratic prism face centers. Geometry queries honor
  `Geometry.OldRuledSurface`. Native caches keep their computed-parameter
  contract; Gmsh can leave newly created face-center parameters unstored.
  That stored-parameter provenance difference is counted separately.
- P2 ownership follows actual classified primary edges and triangle/quadrangle
  faces. Boundary queries include quadratic nodes on the requested entity's
  boundary closure. Supports follow original node identities through selection,
  renumbering, cache adapters, merging and duplicate-node compaction. Native
  mixed construction retains curve and surface carriers after safe point
  compaction, including Hex27 and Prism18 face centers; ambiguous support merges
  reject before cache or attached-record changes.
- Retained legacy Tri6/Tet10 overlays now supply actual geometry and published
  types to Jacobian, inverse-coordinate, locator, basis, key, barycenter and
  edge/face queries. Replaced linear blocks are absent from exact-type queries;
  node-by-element-type queries select actual interpolation nodes by family
  across orders. Lower linear cells and entity/task selection remain supported.
  Quadratic simplex quality queries use actual sampled derivatives and adaptive
  Bernstein bounds, with curved triangle area and corner tetrahedron volume
  following the pinned upstream conventions. Independent edited-geometry
  quality regressions include edited planar and warped surface maps. Normal
  query-resource checks pass 178 assertions: seven quality measures grow
  linearly at 1k/2k/4k cells, single
  queries allocate a constant 128–1,632 bytes, and helper IR has no `Core.Box`.
- Nodal key queries return actual stored interpolation nodes even when the
  requested function-space order differs. Basis counts still describe the
  requested order. Nodal key information emits complete requested-size groups;
  hierarchical information preserves the submitted length and pads an
  incomplete tail with `(0, 0)`, following the two pinned upstream branches.
  Adaptive determinant and isotropy regressions pass 76 checks; independent
  comparison of 64 edited Gmsh fields passes 129 checks. Bounds resource checks
  pass 52 assertions with exact 2x allocation growth at 1k/2k/4k cells and no
  `Core.Box` in the 19 inspected methods.
- Legacy P2 inverse and point-location queries use immutable simplex basis
  tuples instead of allocating basis arrays for every candidate. The measured
  Tet10 inverse query drops from 147,416 to 720 bytes; the 4k-cell location
  fixture drops from 589.8 MB to 2.449 MB. The resource gate passes 177 checks,
  including constant single-query allocation and linear bulk/candidate growth.
- The final normal Julia 1.13.1 Core and eight-file API focus passes 5,425
  assertions in 259.9s, including 432 new hierarchical-tail regressions and
  the actual P2 Jacobian/query/quality/bounds tests. All 93 production hashes
  and 407 staged/runtime inputs remain unchanged. The final v5 metadata
  comparison matches all 192 saved pinned P1/P2 query cases.
  Its 67 resource/IR checks pass with typed returns, no `Core.Box` and
  allocation doubling below 2x at 1,001/2,001/4,001 submitted keys, including
  partial tails. The full normal Julia 1.12.7 function-space differential
  passes in 120.4s: 124 fixed nodal types, 384 actual-order cases and the
  existing orientation/key catalog, plus six explicit-order stored-P1 key
  cases, 36 nodal-tail cases and two hierarchical-padding cases. Its original
  serialized SHA-256 remains
  `b28f429e11b56c08f8b39999b892a7132cdd9d7eed79a5cf2e63835fdf525ac4`.
  All 407 final source/test/validation inputs remain unchanged throughout.
  The final normal Julia 1.12.7 mixed-query differential also passes all
  2,612 checks across 60 cases in 76.3s on those same inputs.
- Simplex and mixed refinement inherit actual child edge and face supports.
  Tagged caches transfer them through the subdivision plan's sparse node-tag
  map. New native/tagged Prism/Hex and raw shared-cell regressions pass 364
  assertions; the final refinement run passes 553 including existing coverage.
  A pinned Prism replay agrees on all four target curve/surface carriers.
- Duplicate-node merging checks retained quadratic geometry before remapping
  carriers. This preserves the existing midpoint-conflict diagnostic while
  keeping strict ambiguous-owner rejection atomic. The unchanged complete
  quadratic certification suite passes 3,061 assertions and the boundary
  ownership suite passes 738 on both normal, bounds-checked Julia versions;
  all 91 production/Project hashes remain unchanged during each run.
- All six affected API differential drivers pass on the v4 frozen source
  with normal Julia 1.12.7, bounds checks and pinned Gmsh 4.15.2: generation
  0D/1D (9,388 assertions, 108 cases, 196 stages), mixed cache (152 checks),
  mixed queries (2,612 checks, 60 cases), mixed refinement (31 cases),
  mesh-data queries and affine transforms. Each wrapper checks all 93 production
  hashes and all 407 source/test/validation inputs before and after each driver.
  Existing explicit upstream locator exceptions and the legacy attached-source
  allocation-history blocker remain separately counted.
- Both Julia versions regenerate the 40 triangle and 24 earlier quadrangle
  artifact records on unchanged source and certificate helpers. The triangle
  files byte-match the tracked pins; all earlier quadrangle rows are preserved.
  Only the two free-lateral P2 digests changed from the pre-ownership draft;
  the other 38 triangle records are unchanged.
- The v4 normal public triangle and boundary suites pass 1,628 assertions on
  Julia 1.12.7 in 238.3s and Julia 1.13.1 in 200.5s (890 triangle and 738
  boundary checks). Artifact regeneration takes 204.8s and 166.4s respectively.
  All 407 source/test/validation inputs remain unchanged before and after each
  run, with independent source and certificate-helper checks. The final
  metadata correction changes no geometry or artifact inputs.
- The v4 strict NoNew differential passes all 86 oracle samples in 234.6s
  under normal Julia 1.12.7 with bounds checks. Six helical oracle gaps and
  two stored-parameter provenance gaps remain separately counted; malformed
  height fixtures still reject. Actual P2 Jacobian queries now satisfy the
  unchanged strict oracle assertions.
- Reserved sweep storage reduces all six measured isolated-triangle allocations
  by 10.2–32.3%. Actual standalone/GEO/projection/API3 1k/2k/4k-layer measurements
  pass 41 checks with maximum allocation doubling 2.323. Actual mixed-part
  support merging drops from quadratic to linear growth in the measured
  100/200/400-part fixture. Support-map and full mixed P2 measurements pass 70
  checks; mixed refinement resource checks pass 47 with maximum doubling 2.214.
  The geometry-first helper comparison passes 33 checks with identical actual
  coordinates, cells and owners, and no new `Core.Box`; staged geometry plus
  owners allocates no more than the previous unified helper in those batches.
  These are bounded measurements, not a universal optimization claim.
- Source grids, mixed roots, collapsed columns, neighboring regions, copied
  chains and cyclic sweeps remain separate NoNew phases. Native curved-CAD P2
  placement, complete higher-dimensional public tag/allocation lifecycle and
  broader mesher/API/format parity remain unfinished. The broad goal stays active.

Isolated NoNew region and actual-map certification increment, 2026-10-04:

- Final normal-compile, bounds-checked package gates pass 484,513/484,513
  assertions on Julia 1.12.7 in 25m34.6s and Julia 1.13.1 in 20m42.0s.
  Both include 38 optional pinned-source provenance checks. All 89 frozen
  production/Project, 151 test and 149 validation input hashes stayed unchanged
  throughout the gates. A trailing blank line in the CRC fixture was removed
  afterward; the native parser verifies all 24 records remain identical.
- `QuadTriNoNewVerts` supports one isolated nondegenerate source quadrangle,
  normalized positive layers, free/recombined laterals, and the documented
  translation/rotation/twist cases. Its operation plan supplies actual volume,
  boundary surfaces and classification. Whole-domain P1 Jacobian proofs and
  bounded global hull separation reject folds and overlapping intervals before
  publication. Grids, triangular/mixed roots, collapsed columns, shared regions,
  copied-source chains and cyclic sweeps remain separate implementations.
- Newly constructed API full P2 products receive speed, Gram or determinant
  certificates for all seven standard families. Pyramid14 uses 144 normalized
  interval/exact tensor-Bernstein coefficients. Independent helper tests pass
  3,408 assertions on each Julia version; the seven-family peer reconstruction
  passes 7,046 checks. Existing imported or user-edited maps remain inspectable.
- Legacy simplex P2 mutations preserve actual midpoints by explicit original
  primary-node identity. Finite primary edits change only the requested node;
  affine maps transform actual stored nodes; repeated order two retains the
  overlay. Coincident components, real compaction, singular maps and conflicting
  midpoint merges have regressions. Curved legacy refine/optimize reject before
  mutation until their placement kernel exists; zero iterations preserve state.
  Exact mixed refinement resets to P1 like Gmsh. Legacy straight-simplex
  refinement retains its existing re-elevated query overlay; exact post-refine
  order parity belongs to the pending higher-dimensional API lifecycle work.
- Final normal public coverage passes 3,604/3,604 assertions in 4m23.5s
  (3,061 P2 and 543 NoNew/certification checks), with source/test hashes stable.
  Independent final P2 peer coverage passes 111/111 with all 88 source paths and
  hashes unchanged. Five affected API files plus allocation/boxing and manifest
  checks pass 585/585 (508 + 76 + 1).
- AddVerts certifies every actual retained Hex/Prism/Pyramid/Tet. Standalone
  grading uses an isolated curve-parameter dictionary and publishes after all
  emission checks; GEO/API also stage failures. Its final regression/resource
  audit passes 436/436, and its 32-case pinned differential passes.
- Malformed GEO scanners and ASCII STL readers close their streams immediately:
  14 immediate-unlink regressions and 416/416 normal IO checks pass without GC.
  The API oracle callback is rooted through unregistration, with forced GC in
  its fixtures. Final API 0D/1D differential passes 9,388/9,388 assertions across
  108 fixtures and 196 stages; mixed-query coverage passes 2,612 checks in 60 cases.
  All six affected API differential drivers exit successfully on the final
  frozen source (generation 0D/1D, mixed cache, mixed queries, mixed refinement,
  mesh-data queries and affine transforms). Strict NoNew coverage passes 20 P1
  and four P2 supported fixtures, two isolated-region cases, four certified
  native helices, six separately counted helical oracle gaps and four expected
  upstream non-unit-height errors, with 44 oracle samples in total.
- Bounded 1k/2k/4k actual-P2 transfer measurements pass 42 checks. Reserved edge
  dictionaries reduce all six measured allocations: Tri6 to
  750,855/1,499,607/2,997,335 bytes and Tet10 to
  1,342,775/2,683,479/5,364,750 bytes, with approximately twofold doubling.
  Cap-chain allocation falls by 87.4–87.8% with unchanged choices. Full NoNew
  volume/GEO/projection/API3 measurements have maximum doubling 2.274;
  regular local certificates allocate zero warmed fast-path bytes. The exact
  correlated Pyramid14 fallback drops from 2,248,112 to 482,872 bytes per cell.
  These are bounded measurements, not a universal optimization claim.
- The 24 native NoNew CRC products include coordinates, oriented typed cells,
  ownership, names and metadata; they agree across Julia 1.12.7 and 1.13.1.
  Pointer-dependent Gmsh splits are certified independently. Six large-angle
  helical oracle gaps are counted separately; excluded oracle products can have
  212 or 213 used nodes, including an actual corner-mean centroid. Those gaps do
  not relax supported fixture checks. Complete API 2D/3D lower-cell publication
  and historical tag allocation, native curved-CAD P2 placement, and broader
  mesher/API/format parity remain unfinished. The broad goal stays active.

API generation-zero/one increment, 2026-10-03:

- The frozen production tree passes the full normal-compile, bounds-checked
  Julia 1.12.7 package gate: 456,607/456,607 assertions in 25m38.0s, including
  38 optional provenance checks against the independently verified Gmsh source.
  Julia 1.13.1 passes the same 456,607/456,607 assertions in 19m33.7s, including
  the same 38 source checks. All 80 frozen production/manifest hashes are unchanged.
- Focused 0D/1D coverage passes 3,522/3,522 assertions in 51.4s. The strict
  Gmsh 4.15.2 differential passes 9,388/9,388 in 36.0s across 108 fixtures
  and 196 stages. One legacy attached-source allocation-history blocker is
  counted separately; the four unsafe raw native P2 DLL fixtures are excluded.
- All six affected pinned differential drivers exit zero under normal Julia
  1.12.7 with bounds checks: mixed cache (152 checks), mixed queries (60 cases,
  2,612 checks), mixed refinement (31 cases), transfinite curves (39 cases plus
  18 HWall cases), mesh-data queries, and mesh affine transforms. The query
  driver's existing upstream exception counts remain explicit.
- Nine complete public API CRC records and allocator maxima match bitwise on
  normal bounds-checked Julia 1.12.7 and 1.13.1. Recipes and pins are in
  `test/interfaces/api_generate01_artifacts.jl` and
  `test/artifacts/api_generate01_crc.txt`. The release allocation audit passes
  76/76 assertions on both Julia versions. Two closure captures discovered
  during the full gates were removed and checked for `Core.Box` regressions.
- Dimensions 0/1 now use detached preparation and atomic model/cache
  publication. Native Point/Line cells, full P2 lines, sparse raw records,
  physical-group priority, `Mesh.MeshOnlyEmpty`, visibility, callbacks,
  fields, graded laws, and immediate order conversion have focused coverage.
- Tagged queries and mutations use public IDs. Mirrored records no longer
  duplicate homology topology. Historical tag maxima persist across clear
  and renumber; attached lower cells and native higher-dimensional source
  boundary provenance survive the verified transitions.
- Immediate order conversion also processes raw meshes before generation.
  Selective clear validates every requested entity before changing records;
  failed mixed valid/unknown selections preserve the entire model/cache.
  Legacy higher-dimensional caches with nonempty native attachments lack
  historical allocation provenance; the dimension-one transition is precisely
  blocked before mutation, with sparse and dense-looking source-tag regressions.
- Periodic nodes support primary/full-P2 correspondence; periodic keys
  follow the pinned seven-output protocol. Node compaction remaps stored
  periodic links and drops pairs involving removed nodes. Identity-preserving
  edits keep global/window visibility and remap flags on element renumbering.
- Six P1 performance fixtures retain identical before/after output hashes.
  Warmed normal-compile Julia 1.12.7 with bounds checks: 333/666/1333 TF3
  curves allocate 12,734,060/21,922,915/49,565,039 bytes, reduced from
  approximately 1.08/2.15/4.31 GB; median times are 7.15/12.08/30.70 ms.
  Single 1000/2000/4000-node curves allocate 2.37/4.46/10.13 MB.
  Full P2 single-curve 1999/3999/7999-node output allocates
  13.27/26.78/53.89 MB; many-curve 1665/3330/6665-node output allocates
  22.05/40.99/86.68 MB. These are measured cases, not a universal bound.

Previous QuadTri increment, 2026-10-03, Julia 1.12.7 with bounds checks and pinned
Gmsh 4.15.2:

- The final frozen-production package gate passes 453,047/453,047 assertions
  in 25m00.2s. Julia 1.13.1 passes 453,085/453,085 in 19m29.6s, including the
  38 optional provenance assertions from the verified Gmsh source checkout.
- Original aggregate coverage is accounted for as prefix 10 + A 11 + B 23 +
  C 23 = 67 original child drivers, all passed through exact resumptions.
  This was not one uninterrupted aggregate process. C's 23 children are
  prefix 9 + sheet 2 + shell 1 + periodic 2D 3 + periodic-volume 1 + final 7;
  the final seven exit zero in 423.053s. The five original analytic benchmark
  rows and final coax forensic probe also completed: flat volumes are exactly
  2/24/35, cylinder volume is 62.652572 and sphere volume is 20.102034. The
  literal coax probe handles its expected upstream nonzero exit, empty volume
  cells and invalid duplicate-surface partial mesh.
- Transfinite QuadTri transitions pass 22,872 assertions over all six-face
  masks, arrangements, prism masks, positive analytic volumes, classified
  projection, physical ownership, deterministic output, and allocation
  limits. A separate 22-assertion certificate suite covers warped, folded,
  flat, duplicate-boundary, and missing-boundary cases.
- QuadTriAddVerts fan certificates and surface/volume boundary coverage
  pass 334 focused assertions. The pinned differential passes 32 fixtures
  with identical node counts, coordinates, per-entity type counts, and every
  typed cell vertex set, including toroidal and neighboring sweeps.
  Transfinite transitions pass 117 independent differential cases.
  Twist keeps Gmsh's distinct orphan
  control-point nodes: 90 nodes with free laterals, 74 with recombined
  laterals. Fixed-column revolutions preserve source/copy curve identity.
- Parser/geometry extrusion regressions pass 455 assertions, including
  supported nested numeric selectors, lexical variants, syntax precedence,
  side-effect rejection, and fixed-edge coherence.
- Deleted-curve discretization and per-model field allocation fixes have
  focused regressions. Rotation checks preserve 720 ordinary-case bitwise
  comparisons and allocate zero bytes for 1,000 warmed scalar rotations;
  extreme-axis and near-parallel cases use an analytic scale oracle.
- Fixed-size transition-cell extraction reduced warmed kernel allocations
  from 174,450 to 130,034 bytes at n=2 and from 5,947,070 to 3,872,894 bytes
  at n=8. AddVerts n=4/n=8 allocations are 945,398/2,259,490 bytes for three
  layers. No Core.Box remains in the audited rotation/transition helpers.
- Flat or zero-volume sweeps and unrecombined quad-source sweeps now fail
  explicitly instead of returning invalid output, even where Gmsh emits
  zero-volume elements. `QuadTriNoNewVerts` remains an explicit blocker.
- The API retains native mixed blocks/classification, including linear and
  full quadratic standard families, through order changes, transforms,
  duplicate removal, deterministic partitioning, and uniform refinement.
  Reference/Jacobian/location, basis-key/orientation, and full edge/face node
  queries operate on the native element maps. Mixed query tests pass 319
  assertions; the pinned differential passes 60 fixtures and 2,612 checks.
  Global primary topology passes 12 mixed-family oracle fixtures.
- Cache lifecycle tests pass 160 assertions and 152 pinned oracle checks,
  including independent coincident triangular/quadrangular surfaces at orders
  one/two and corrected entity-filtered simplex P2 connectivity.
- Mixed refinement passes 189 assertions and 31 pinned oracle cases, comparing
  all typed cell vertex sets and node coordinates, repeated refinement, warped
  cells, shared hex/prism and pyramid/tet faces, orphan pruning, and positive
  volume normalization. Node identity passes eighteen pinned geometry fixtures;
  focused tests also cover imported homology tags and stale mesh records.
- The query oracle records Gmsh's standalone-Point octree omission, unused-w
  Newton artifact on warped quadrangles, and inconsistent strict P2 pyramid
  octree results separately. Tessella's analytic membership stays enabled;
  no tolerance was weakened. P2 location uses conservative candidate bounds.
  Unsupported orders and native quality names fail explicitly.
- Curved native CAD elevation/refinement fails before cache mutation, pending
  curve/surface placement stencils. Straight/discrete P2 conversion remains
  available. The native quarter-annulus Q9 chord center differs from Gmsh's
  curved stencil center, independently confirming this required guard.
- P2 basis-derivative reuse reduced 100-hexahedron, one-point bulk Jacobian
  allocations from 734,174 to 89,262 bytes (87.8%). Native disconnected-cell
  partition allocations were 392,732/784,796 bytes at 2,000/4,000 cells.
  These are bounded measurements; curved P2 point inversion remains costly.
- Six full QuadTri CRC records are pinned in `test/artifacts/quadtri_crc.txt`
  and tests. They match bitwise under normal bounds-checked Julia 1.12.7 and
  1.13.1. Raw record nodes retain their owning entity, even at model endpoints;
  coincident raw node tags within one record require a future tag-identity
  merge and are rejected precisely rather than silently welded.
- Windows external size-field launch now uses the absolute command shell and
  its native quote rules. Twenty-seven protocol checks and a three-coordinate
  pinned oracle pass, including quoted helper/batch paths, missing System32
  from PATH, abnormal peer exits, and cancellation. Validation helpers also
  capture Gmsh's version portably and forward the requested algorithm; real
  CLI algorithm-one/ten option guards pass.
- An independently reproduced overflow reflection now uses the exact
  determinant sign for native mixed caches and mesh records. The matrix
  `1e150*[-3 2 1;4 1 2;1 3 5]` has negative determinant while its floating
  cofactor sum becomes NaN; five regressions verify finite transformed
  coordinates and positive normalized tetrahedron volumes.
- Native Point attachments replace the synthesized geometry-position node.
  The pinned oracle retains one mesh node at `(2,0,0)` for a Point at the
  origin with that stored mesh vertex; native output and H0 rank now match.
  Forty-seven identity-helper and 87 homology/identity assertions pass.
  Homology preserves explicit stored Point elements and otherwise creates a
  single Point cell on the last stored mesh vertex, matching upstream even
  when other vertices remain classified on that Point.
  Recovered curve samples share explicit co-embedded Point constraints with
  the actual carrier cells. Canonical Point incidence uses the first stored
  vertex; extra stored vertices remain separate. Raw curves and unrelated
  coincident entities preserve their original identities.
  Displaced native mesh vertices used by curves or embeddings reject before
  meshing until their boundary placement semantics are implemented.
- The initial full bounds-checked Julia 1.12 gate passed 452,978 assertions,
  with 20 stale checksum assertions and no test errors. Independent clean
  `43183bd` comparisons reproduce all eleven distinct expected-hash changes
  with bitwise identical coordinates, connectivity, tags, and classification.
  The three API box/refinement records also match across Julia 1.12 and 1.13;
  the updated API/lifecycle/transform files pass 675 assertions on Julia 1.13.
- Two later embedded-sheet validator assertions were separately corrected:
  the classified-projection and MSH2 checksum pins. Clean `43183bd` and current
  Julia 1.12 snapshots match every field and coordinate bit across 14 stages
  for the sheet and holed-sheet fixtures, including merged/pure/projected
  meshes and ASCII/binary MSH2/MSH4 rereads. The holed-sheet hashes were already
  accepted and remain unchanged. Both full pinned drivers pass (176.700s and
  181.128s). These two validator pins are separate from the initial twenty
  stale unit assertions.
- Periodic validation drivers now pass the entity mesh to the existing
  surface CRC/projection checks and validate the merged global mesh separately.
  Six fixtures at three stages match clean `43183bd` exactly; their existing
  surface-mesh checksum pins remain unchanged. This corrects the validation input
  contract rather than weakening an identity or connectivity assertion.
- A third later validator pin, the periodic-volume MSH2 checksum, was corrected
  after a separate clean `43183bd` comparison. Global, pure-volume, classified
  projection and four ASCII/binary MSH2/MSH4 rereads match every metadata field
  and coordinate bit. Both MSH2 modes produce
  `b8251bd55dc17ab832e64c1cb976936fb22bdb4777587095334c6d4b3716546f`.
  This pin is separate from the two sheet-validator pins and the initial
  twenty unit assertions. The full periodic-volume driver passes in 167.977s,
  retaining its geometric/periodic assertions and all four file modes.
- The aggregate function-space driver exhausted memory while retaining its
  entire multi-gigabyte checksum stream. Incremental SHA encoding preserves
  the original checksum and every full orientation comparison with a fixed
  65,536-byte buffer (65,800 bytes retained). Twenty-six independent old/new
  byte-protocol checks pass; warmed 250,000-float encoding allocates zero bytes.
- The available Gmsh 4.15.2 source checkout passes all 38 optional provenance
  assertions, plus an independent version assertion. The full Julia 1.13
  gate uses that verified checkout explicitly and passes all 453,085 assertions.

The following record describes the prior `43183bd` increment:

Re-verified on 2026-10-03 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing the `.geo` `Extrude … Layers` sweep:

- `meshGRegionExtruded` parity: a `Layers`-marked extruded volume now sweeps
  its source surface mesh through every layer-group level instead of
  falling back to an unstructured fill. Level parameters come from the
  parsed `Layers` records (uniform counts, per-group counts, and explicit
  height fractions normalized to `u`), and every swept corner is a pure
  `_extrude_at` transform evaluation — upstream's
  `pos.find(Extrude(u, pos))` semantics — so translation, rotation, and
  twist (`{{T},{axis},{P},a}`) cases weld bitwise at merge time. Twist
  connectors keep their off-path spline endpoints as curve endpoints only
  (identical to Gmsh's own connectivity); connector curve parts emit the
  transform evaluations like `extrudeMesh(GVertex)`/`copyMesh`.
- Element emission matches the oracle type-for-type: triangle generatrices
  produce prisms under `Recombine` (21/21 connectivity sets bitwise on the
  verification fixture) and subdivide to three tetrahedra per prism without
  it through the global `SubdivideExtrudedMesh` phase-1/2/3 shared-diagonal
  pass, recombined quadrilateral generatrices produce hexahedra, lateral
  surfaces emit quad strips (recombined) or the shared-diagonal triangle
  remesh, and top surfaces copy the source mesh verbatim. Layer groups and
  nonuniform heights reproduce the oracle's 5-level prism stack exactly.
- `setAllVolumesPositive` normalizes cell orientation across tets, hexes,
  prisms, and the newly face-tabled pyramids (`_VOLUME_CELL_FACES[7]`).
- Verified against the oracle: element-type counts and node counts
  identical on translation (tri/quad sources, `Recombine` on/off),
  rotation, twist, multi-layer, and no-`Layers` fallback cases; merged
  output is conforming (`validate` clean, zero nonconforming boundary
  faces) on every nondegenerate case. The one degenerate case — a flat
  in-plane rotation about an axis through the source — emits Gmsh's
  identical 160 zero-volume tets, which `validate` correctly flags.
- Attainable-parity bound recorded: non-recombined tet splits can pick a
  different valid member of the diagonal-assignment class than Gmsh because
  the source triangles' cyclic vertex order is mesher-internal, and
  recombined-quad hex connectivity cannot match until the surface
  recombiner reproduces Gmsh's interior quad layout; counts, node sets,
  and total volume are all identical in both cases.
- `QuadTriAddVerts`/`QuadTriNoNewVerts` extrusions raise an explicit
  blocker naming the QuadToTri kernel instead of silently sweeping wrong
  output.
- Gates green: `geo_constraints_test` extrusion set (308 assertions) plus
  the full `Pkg.test()` run — the 20 pre-existing stale CRC pins it
  surfaced were verified bitwise-identical between `d53490f` and this tree
  and repinned (`test/geometry/*` and `test/interfaces/*` files only).

Re-verified on 2026-10-01 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing three-sided ruled transfinite fills:

- Three-generatrix `Surface` loops (gmsh `MSH_SURF_TRIC`, geomType
  RuledSurface) — including non-coplanar boundaries and concentric-arc
  sphere patches — now mesh instead of rejecting on the planar gate. The
  kernels keep their real-space `TRAN_TRI` interior; a new `project`
  callback inverts each interior point through `_ruled_xyz_to_uv`'s new
  `on_surface=false` variant (the exact `XYZtoUV(..., 1.0, false)`
  convention: `Precision = 1e-3`, `MaxIter = 10`, 9×9 restart grid,
  last-iterate return on non-convergence) and re-evaluates the ruled
  `point(Up,Vp)`. On the sphere-octant fixture (three unit great-circle
  arcs, `Transfinite Curve = 5`) both kernels reproduce Gmsh 4.15.2's
  node sets within ~2.3e-9: collapsed 22 nodes / 28 triangles, compact
  `Mesh.TransfiniteTri=1` 16 / 16 — every interior node on the unit
  sphere, boundary nodes bitwise on the true arcs.
- The planar gate is untouched when `project === nothing` (flat patches
  byte-identical), and existing on-surface `XYZtoUV` callers keep the
  strict 1e-8/25-iteration recursive-relaxation behavior.

Re-verified on 2026-10-01 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing `In Sphere` transfinite fills:

- Four-sided ruled surfaces carrying sphere geometry — `Surface …
  In Sphere{p}` or four concentric arc generatrices (the same
  `checkSphere` detection upstream applies) — now mesh through the
  surface's own `S(u,v)` evaluation: `mesh_transfinite_patch`'s new
  `interpolate` callback runs the `TransfiniteQua` generatrix blend plus
  `TransfiniteSph` projection (radius |S0−O|) on the averaged-chord grid,
  replacing the Coons interior while boundary nodes keep their welded
  side-chain positions. Coplanar and non-coplanar boundaries both mesh
  and bitwise/ulp-level agree with Gmsh 4.15.2: the flat-corner fixture
  reproduces gmsh's 26-node/32-triangle patch (interior nodes at
  2.1213203436 from `(0.5,0.5,-2)`, e.g. (0.2427816723, 0.2427816723,
  0.0898989123)) — differences ≤ ~2e-12, gmsh's own transfinite-curve
  Newton-spacing noise — and the four great-circle-arc unit-sphere patch
  matches within ~1e-9.
- The volume kernel consumes the same evaluation through
  `_transfinite_volume_face_grid`: an `In Sphere` top face of a
  transfinite cube produces gmsh's bitwise-identical projected node
  (0.5, 0.5, 1.1213203435596424) and a volume interior node matching to
  ~1.5e-12.
- Unsafe frame mismatches reject precisely: a degenerated-curve skip or
  a corner reorder would desynchronize the generatrix evaluation from
  the kernel grid (identity pinning stays supported).
- Gates green: `transfinite_test` 226/226 (new sphere assertions),
  `transfinite_volume_test`, `transfinite_hex_test`,
  `geo_curved_test`, `model_test` 104/104, `geo_mesh_dim_test`,
  `geo_geometry_expression_test`, `geo_periodic_test`, `api_test`,
  `geo_constraints_test`.

Re-verified on 2026-10-01 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing the all-dimension `execute_geo`
emission increment:

- `execute_geo(path; mesh_dim=n)` now routes every model through
  `_geo_mesh_model` like the `Mesh n` statement: vertex parts (orphan
  `Point` entities included — e.g. an arc's control-point entity),
  graded curve parts, surface parts, and volume parts merge on bitwise
  coordinates. Verified kwarg-vs-statement bitwise-identical on the
  periodic-surface volume, its rotated variant, geometry-expression,
  list-variable, SetMaxTag, PointsOf, and all periodic dim-2 fixtures;
  Gmsh 4.15.2 emits the same element mix (126 nodes / 633 elements on
  the curved-box oracle fixture).
- The per-entity decomposition stays on `GeoExecution.mesh_parts`;
  `geo_entity_mesh(execution, dim, tag)` recovers each entity's own
  mesh for classified `model_to_mixed` projection, which keeps its
  strict pure-entity contract and reproduces the pre-unification output
  bitwise — every `mixed_crc` pin across the geometry and CLI suites is
  unchanged. The CLI's `.msh` writer now serializes vertex/edge/face/
  region cells like `gmsh -n`.
- OCC `gp_Lin` parameters carry native (possibly nonzero-origin) bounds:
  `_periodic_curve_point` no longer evaluates them as chord fractions,
  and `_curve_parameter_nodes`/`_embedded_line_curve_nodes` emit native
  frames with endpoint-exact snapping — a latent bug that broke OCC
  Boolean-difference face meshing (constraint overlap / degenerate
  triangles / broken endpoint chains). The BooleanDifference fixture now
  meshes all six faces plus the volume (299 nodes / 72 segs / 592 tris /
  12 tets, `validate` green, unit volume).
- OCC-curved boundary faces (cylinder/sphere/cone/torus side walls that
  cannot standalone-mesh) skip their surface part when they bound a
  volume — the facets ride volume boundary recovery — while the OCC
  Cylinder kwarg path now matches the `Mesh 3` statement path.
- The `Include`/`Merge` relative-path helper's shadowed `in` binding (a
  pre-existing `MethodError` on absolute Windows child paths) is fixed.
- Caveat recorded: `mesh_model_volume`/`model_to_mixed` output is
  sensitive to `--check-bounds` codegen (borderline FP decisions in the
  Delaunay/refinement path); all pins are therefore verified under the
  suite's `--check-bounds=yes` mode, and CRC probes must run with that
  flag.
- Gates green: `model_test` 104/104, `geo_periodic_test`,
  `geo_mesh_dim_test` 159/159, `geo_curved_test`, `geo_constraints_test`,
  `geo_mesh_size_test` 242/242, `geo_geometry_expression_test` 58/58,
  `geo_list_variable_test` 66/66, `geo_set_max_tag_test` 95/95,
  `geo_dynamic_tag_test` 143/143, `geo_discrete_statements_test`,
  `transfinite_test` 193/193, `model_volume_io_test`,
  `model_periodic_io_test`, `occ_primitives_test`, boolean multi/snapshot
  and lifecycle files, `sizefield_test`, `cli_test`, `io_test`,
  `pipeline_test`, and the remaining geo files.

Re-verified on 2026-10-01 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing the non-coplanar `Plane Surface`
transfinite increment:

- `_transfinite_declared_plane` reproduces upstream `computeMeanPlane`:
  two on-curve samples per boundary edge at ~1/3 and ~2/3 of each curve's
  native parameter bounds, then the plane through the first non-collinear
  triple — off-plane arc control points never veto the fit.
  `mesh_transfinite_patch`'s new `project_plane` mode projects side chains
  for (u,v) bookkeeping and Coons interpolation (interior lands exactly on
  the declared plane) while emitted boundary nodes keep true positions;
  the warped-patch audit runs on the emitted band. Verified bitwise: the
  lifted-corner `Plane Surface` interior lands on `z = 0.5y` with all 9
  interior nodes bitwise-identical to Gmsh 4.15.2.
- The transfinite-volume boundary-fold audit now certifies each emitted
  boundary triangle by the incident tet's apex and escalates a sign
  straddle to an exact edge-through-triangle pierce test over the cell's
  19 subdivision edges (`_segment_crosses_triangle`, exact `orient3`).
  Outward and shallow-inward curved-face bulges mesh with min tet volume
  matching Gmsh to ~9 digits; the strong inward bulge stays rejected —
  Gmsh 4.15.2's own output there is self-intersecting (verified: a tet
  edge pierces a boundary triangle), which Tessella refuses to emit.
- New coverage: `transfinite_test` asserts the lifted-corner patch meshes
  planar-interior with true boundary positions (25 nodes / 32 tris like
  Gmsh) plus an off-plane-arc `Plane Surface`; `geo_constraints_test`
  adds outward + shallow-inward curved-front-face volume fixtures
  (384 tets, all-positive, `validate` green, ≥3 off-plane arc nodes) and
  pins the strong-inward rejection.
- Gates green: `transfinite_test` 192/192, `transfinite_volume_test`
  274/274, `geo_constraints_test`, `geo_curved_test`, `geo_periodic_test`,
  `model_test`, `geo_mesh_dim_test`, `mesh3d_test`, structured
  quad/prism/hex/triangle/curve files, and `allocation_audit_test` all
  pass; `transfinite` (4 arrangements, max err 1.6e-15),
  `transfinite_volume` (2 cases, 9.5e-12), `geo_constraints` (32),
  `geo_curved` (8), and `geo_splines` (20) differentials green against
  Gmsh 4.15.2; full suite 427,308/427,308 under `--check-bounds=yes`.

Re-verified on 2026-10-01 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing the closed-curve grading floor and
`Min*`/`Minimum*` mesh-option alias increment:

- `_model_minimum_curve_segments` now applies the upstream closed native
  curve floor: `max(np, 3)` segments for `a==b` non-`Line` native curves
  (Gmsh's `N = minimumMeshSegments + 1` node target becomes `N` edges on a
  closed loop — a two-segment digon is unrealizable). Verified against
  Gmsh 4.15.2: a coarse closed spline emits 3 elements (nodes at arc-length
  thirds, matched to ~1e-9), an open spline keeps its 2-element floor, and
  `Mesh.MinimumCurvePoints = 5` raises the closed spline to 4 segments —
  the raised regular floor still wins over the closed floor.
- `_geo_store_option_number!` synchronizes the Gmsh-canonical spellings:
  `MinCircleNodes`/`MinimumCircleNodes`/`MinimumCirclePoints`,
  `MinCurveNodes`/`MinimumCurveNodes`/`MinimumCurvePoints`, and
  `MinLineNodes`/`MinimumLineNodes` each address one `CTX.mesh` slot
  upstream (`opt_mesh_min_*_nodes`), so a write to any spelling now
  mirrors to the siblings the graders read.
- New `geo_spline_test` coverage (7 assertions): default closed-spline
  floor, open-spline floor, and raised floor via the canonical
  `Mesh.MinimumCurvePoints` spelling.
- Gates green: `geo_spline_test`, `model_test`, `geo_mesh_dim`,
  `model_volume_io`, `geo_periodic`, `geo_curved`, `cli_test`, and
  `allocation_audit` files all pass; `geo_splines` (20 cases),
  `geo_curved` (8), and `geo_constraints` (32) differentials green against
  Gmsh 4.15.2; full suite 427,276/427,276 under `--check-bounds=yes`.

Re-verified on 2026-10-01 on Windows with Julia 1.13.1 against the pinned
Gmsh 4.15.2 binary after completing the curved `Curve In Volume` increment:

- `Curve{t} In Volume{v}` and `Curve{t} In Surface{s}`-in-volume no longer
  require `Line` curves. `mesh_model_volume` seeds the full stored 1-D
  discretization of an embedded non-`Line` curve as interior points
  (`restoreEmbeddedEdges` semantics), lazily graded through
  `_volume_embedded_curve_params!` on direct calls. Embedded sheets recover
  before embedded curves so a nested curve's links arrive as sheet edges —
  the per-link `recover_segment3` loop degenerates to a coverage audit and
  registers realized edges into the protected-cell registry.
- Closed curves (`a==b`) work end to end: the chain wraps to the shared
  vertex, the duplicated `t1` endpoint is dropped from seeding, and all
  three `model_to_mixed` projection paths guard the far-point map and emit
  the wraparound edge once. Tet-mesh curve classification runs through
  `_segment_chain3` (BFS corridor path over tet edges on each chord), so
  recovery subdivision Steiner nodes ride the chain and foreign corridor
  vertices are skipped; sibling-curve-owned vertices are masked via
  `_curve_owned_coordinates`. A sheet-nested curve prefers its host sheet's
  face-edge complex with a tet-edge fallback that registers realized edges
  into the sheet complex for the nested validator. Self-overlapping closed
  curves emit each covered edge once (`edge_owners` in all three line-cell
  loops, atomic two-phase `classify_curve!`).
- `Mesh3D` repairs: `_snap_to_plane3` leaves vertices within interpolation
  noise of the plane unmoved (per-triangle Float64 planes differ ~1e-15 and
  re-snapping walked shared vertices off exact model coordinates, defeating
  bitwise dedupe); `recover_triangle3` snaps crossing candidates along the
  dominant normal axis so boundary-face coordinates stay bit-exact; both
  segment hit paths skip crossings within `1e-6` of an existing on-segment
  station (Zeno cascade that minted sub-resolution twin vertices); a pocket
  seed within `1e-9` of an existing vertex no longer duplicates it.
- New `model_test.jl` coverage (19 assertions): open arc in volume, closed
  spline in volume, self-overlapping closed spline, nested open arc and
  nested closed spline on an embedded triangular sheet, and the `.geo`
  `Curve{…} In Volume{…}` statement — all with `validate` and
  `model_to_mixed` checks plus duplicate-line-cell audits.
- New `geo_constraints` differential case `curve_in_volume_arc`: Gmsh 4.15.2
  and Tessella both carry off-chord nodes on an embedded circular arc in a
  built-in-geometry volume (`GEO_CONSTRAINTS_DIFFERENTIAL_OK cases=32`).
- Known limitation (pre-existing on the pre-change worktree): cavity refill
  can starve on coplanar vertex soup — a closed spline inside a square
  embedded sheet fails in `recover_segment3`/`recover_triangle3`; plain
  square sheets, open curves on them, and closed curves on triangular
  sheets all pass.
- Focused gates green under Julia 1.13.1: `model_test` 345/345 (104 in the
  embedded-constraints testset), `geo_mesh_dim` 159/159, and the
  constraints/curved/spline/periodic files clean.

Re-verified on 2026-09-29 on Windows with Julia 1.13.1 and 1.12.7 against the
pinned Gmsh 4.15.2 binary after completing the holed embedded-sheet recovery
and projection-audit repair increment:

- `recover_segment3` dead-end escape: a Zeno-like cluster of near-duplicate
  Steiner vertices (each insertion minting a new crossing ~1e-8 past the
  `_node_at3` gate) previously starved every legal cavity fill because the
  fill must retain all interior vertices. `_absorb_segment_pair3` now refills
  the pair's union star with the whole absorbable cluster removed from the
  fill pool (`_refill_segment_cavity3` `absorb_verts`), then
  `_compact_nodes3` drops the orphan coordinate columns and remaps cells and
  the protected registry. `_absorbable3` pins any vertex referenced by
  classified segments, triangles, or protected cells; absorption is rejected
  when a vertex is referenced outside the cavity or is required by an
  explicit edge/face, and `absorb_p`/`absorb_q` flags prevent a call from
  absorbing its own constraint endpoints. The chain loop rescans from the
  cursor (`_segment_chain_points3`) after each sub-recovery because
  compaction renumbers node ids, and whole-segment coverage is re-verified
  after every sub-segment.
- `_model_projection_volume_surface_faces!` now audits the GENERATED
  surface triangulation (what recovery enforces — upstream
  `allEmbeddedFaces` semantics) rather than the coarse boundary+embedded-
  points CDT whose internal edges need not exist in the tet face complex.
  `model_to_mixed` rebuilds targets from `mesh_model_surface`;
  `_mesh_model_volume` passes the actually recovered sheet triangles via
  the new `targets` kwarg.
- `embed_sheet_hole`: full differential green —
  `tessella_volume=0.9999999999999992`, `sheet_area=0.45`, `hole_centroid_hits=0`,
  26/26 sheet triangles recovered, MSH2/MSH4 ASCII+binary round trips. The
  Windows recovery path legitimately produces a different (structurally
  valid) mesh than the recorded pin, so the projected/MSH2 CRCs are explicit
  whitelists.
- `Pkg.test()` Julia 1.12.7: **425,608/425,608 in 17m17.5s**.
- `Pkg.test()` Julia 1.13.1: **425,608/425,608** (two `mixed_crc` sha pins
  repinned after every structural assertion passed independently; the
  embedded-sheet face sets follow the generated triangulation).
- Full `validation/run_all.jl` green under 1.13.1 + Gmsh 4.15.2 on Windows —
  all 58 driver gates including the embedded-surface, periodic, and
  differential matrix. Note for reproducers: the driver runs children with
  `--check-bounds=yes`, which changes mesh FP output versus a bare run —
  always reproduce failures with that flag. Windows oracle/Julia
  subprocesses use `Sys.BINDIR` discovery (`Sys.which` and POSIX quoting are
  unreliable on Windows).

Re-verified on 2026-09-27 with Julia 1.13.0 and 1.12.7 after completing the
embedded-recovery determinism and API-review repair increment:

- `recover_triangle3`/`recover_segment3` now converge on holed and foreign-
  vertex embedded sheets: covering faces are registered in the protected-cell
  registry every iteration, pocket dispatch processes all distinct pockets per
  pass, graft candidates require strict 2-D enclosure with no cavity vertex on
  their edges or interior (vertex-on-face legality, matching `inside_cavity`'s
  rule), an unseeded pool fill covers pockets with no clean seed face, and a
  monotone-coverage acceptance gate rejects any fill or insertion that
  regresses target coverage — eliminating the oscillation the registry alone
  could not prevent. Sheet fills pass `soft_keepfaces=true` so coverage, not
  constraint count, is the progress measure; curve recovery keeps hard
  constraints.
- Output is byte-identical across Julia 1.12.7 and 1.13.x on every exercised
  fixture after sorting every `Dict`/`Set` iteration that influences output
  order: the refinement edge queue and deferred requeue, CDT region and
  missing-crease order, piercing-candidate edges, exact-Delaunay cavity faces,
  segment-cavity boundary faces with a `claimseq` DFS ledger, cavity-growth
  incidence, and the sheet gap edge set.
- `mesh_model_volume` now scopes the task-local protected-cell registry:
  emptied on entry and in `finally`. A registry left populated by an earlier
  call in the same task previously leaked foreign keep-constraints into
  `refine_to_size`'s cavity splits — same canonical `mesh_crc`, different raw
  cell order, different `mixed_crc` per process. Verified by poisoning the
  registry before `execute_geo`: the projected mesh still produces the
  expected `ffd2559d…` CRC.
- `Mesh.MeshSizeExtendFromBoundary` stores integer semantics; the documented
  `0` disable setting is accepted and distinguishable from an absent value,
  and float inputs truncate toward integers. Reversed-volume refinement
  preserves parent orientation and `tet_tag`, and the classification skeleton
  carries segments, triangles, and tags through.
- Focused gates green: `mesh3d` 146,746/146,746 on both versions, closure-
  boxing audit 76/76, refinement-classification 18/18, `.geo` constraints
  suite clean. A full `Pkg.test()` under Julia 1.13 on this final tree is
  still pending — earlier suite runs surfaced only stale-precompile failures
  that each reproduce expected values on current code; see HANDOFF.md.

Re-verified on 2026-09-26 with Julia 1.13.0 against the pinned Gmsh 4.15.2
binary and vendored 4.15.2 source after completing the periodic explicit-volume
conformity increment and the 2-D/3-D interior-density work beneath it:

- `Periodic Surface` declarations now install derived curve masters per
  resolved boundary pair, matching `GFace::setMeshMaster`'s per-edge
  `GEdge::setMeshMaster` calls — the same induced pairs the pinned binary logs
  as "Setting curve master". Slave curves mirror master `curve_params`
  positionally instead of grading independently, and `ModelPeriodicConstraint`
  carries a `derived` flag so `model_periodic_constraints` still reports only
  declared relations while the internal graph sees the full upstream set.
- Slave curve evaluation routes through the stored master-parameter pairing
  (`affine(master(u_paired))`), so adjacent native faces and periodic copies
  emit bitwise-identical boundary nodes; cyclic dim-1 relations detect the
  closing chain and evaluate natively instead of composing affines around the
  cycle. Periodic-surface boundary nodes are emitted through each slave
  curve's own relation — upstream `copyMesh` semantics — while slave endpoints
  keep their model `GVertex` coordinates (upstream behavior; the affine snap
  applies to interior nodes only).
- `mesh_curve` pins open-curve endpoint parameters to `t0`/`t1` — the
  inversion at `total*nedge/nedge` rounded below `total` and produced `1−ε`
  params that wrote back as 1-ulp-duplicate corner nodes. Writeback and PSLG
  emission normalize near-endpoint stored parameters so stale `1−ε` entries
  cannot demote a true endpoint.
- The surface PSLG builders record canonical `(u,v)→3-D` coordinates for
  entity-derived vertices and snap lifted output nodes to them, closing the
  per-face `_plane_dropped_coordinate` bitwise divergence that cracked shared
  curves on explicit-volume shells.
- `_model_projection_boundary_surface_faces!` accepts near-coplanar faces
  under a scale-relative plane-offset tolerance alongside exact `orient3==0`,
  admitting affine-copied slave nodes sitting ~1 ulp off the rotated model
  plane.
- The rotated explicit periodic cube (`Periodic Surface 4 ← 6 Rotate π`
  about the volume center) meshes natively and validates: 138 nodes/459 tets,
  `validate` clean; the translational fixture meshes 142 nodes/484 tets.
  Slave periodic surfaces carry 21 nodes each — matching the pinned binary's
  `.msh` output (4 corners + 8 curve + 8 interior + 1 embedded) — and each
  rotation variant pins its own connectivity CRC since Gmsh's own rotated run
  produces different connectivity than the translated one.
- Supporting repairs in the same increment: `mesh_curve`/`_invert_primitive`
  endpoint pinning, `norm`→`hypot` in `_model_mesh_bbox`, a restored
  `_model_surface_plane` gate before lazy curve grading so corrupt curved
  geometry still rejects with `ArgumentError`, boxed-closure audit fixes via
  a callable LL-resolver struct and single-assignment captures, and the
  guarded all-zero-subvolume branch in `_insert_steiner3`.

Re-verified on 2026-09-25 with Julia 1.13.0 against the pinned Gmsh 4.15.2
binary and vendored 4.15.2 source after completing the native model-level
`.geo` `Mesh 0`/`Mesh 1` increment (the `meshGEdge` port):

- `Mesh 0` synchronizes the physical view and returns like upstream's
  no-op `GModel::mesh(0)` arm; `Mesh 1` grades every curve through the
  `BGM_MeshSize` composition (`l1` vertex min, `l2` boundary edge-length,
  `l3` background field, `l4` entity size, `l5` parametric sizes, then the
  `lc`/`lcMin`/`lcMax` clamp and `lcFactor`), the `minimumMeshSegments`
  floors (`Mesh.MinLineNodes`/`MinCircleNodes`/`MinCurveNodes`), and the
  recombination odd-node/`increaseN` adjustments. `filterPoints` matches
  upstream's `(lc/d)`-sorted removal including the `recombineAll`-only
  even-count truncation. Curves store native-range parameters in
  `GeoModel.curve_params` (the `GEdge::mesh_vertices` analogue) with
  `deMeshGEdge` lifecycle coverage across entity removal, retagging,
  OCC rollback, `Delete`/`NewModel` resets, and `mesh.clear()`.
- `Degenerated Curve`, coincident-endpoint `Line{p,p}` self-loops
  (warning text and single-element emission match the binary), discrete
  curves (parametrized re-grade; non-parametrized nodes kept), extruded
  curves (`ExtrudeParams::u` interpolation for point generatrices and
  signed-source `copyMesh` for top copies through the new
  `extrude_sources` map), and periodic slaves (`copyMesh` forward-shift
  and `u_max - u + u_min` mirrored forms) all reproduce upstream.
- `Mesh 2`/`Mesh 3` reuse the stored discretization — boundary
  `forced` PSLG parameters seed from `curve_params` and the face
  refinement splits write back onto the curves like `meshGFace`, so the
  emitted line elements are exactly the triangulation's boundary edges.
  Periodic slave/master surfaces and periodic-involved curves are
  excluded from writeback so slave meshes stay bitwise-identical mapped
  copies, and nodes on periodic chains are protected from unrelated
  endpoint snapping.
- `.geo` `Save` emits MSH 2.2 type-15 point elements (vertex-tag order,
  physical from the model view) matching upstream's per-vertex `MPoint`
  serialization; the MSH 4.1 entity-block layout rejects them
  explicitly. `write_msh` tolerates degenerate self-loop segments for
  serialization while `validate` keeps them invalid by default.
- Homology cell chains consume the stored `curve_params` evaluated
  through the part evaluator instead of reconstructing attributes.
- `Mesh.MeshSizeFactor`/`CharacteristicLengthFactor` alias mirroring
  extended to `MeshSizeMin`/`CharacteristicLengthMin` and
  `MeshSizeMax`/`CharacteristicLengthMax`; `Background Field`/
  `BoundaryLayer Field` selections persist on the numeric context and
  feed the raw `l3` field into curve grading.
- Differential results against the pinned binary: `Mesh 1` node/segment
  counts and parameter sets match exactly on plain, transfinite,
  point-sized, minimum-node, periodic, re-graded, and zero-length
  fixtures (6/5, 8/6, 13/12, 17/16, 18/16, 91/90, 2/1 node/segment
  pairs); `Mesh 2`/`Mesh 3` boundary element counts match exactly
  (32/64 line elements) with interior simplex density remaining the
  documented 2-D/3-D mesher gap. The full repository gate (423,971
  assertions) surfaced exactly four regressions — allocation-audit
  boxed closures, two stale CRC pins, and a `mesh.clear()`
  `curve_params` lifecycle miss — each fixed and re-verified green in
  its own testset (audit 76/76, mesh-transform 92/92, periodic-io
  153/153, CLI 131/131) with no other suite affected.

Re-verified on 2026-09-21 with Julia 1.13.0 against the pinned Gmsh 4.15.2
binary and vendored 4.15.2 source after completing the `.geo` dynamic-allocator
and option-expression parity increment:

- `NumberOption` expression reads and writes now cover the whole `x[i].y`
  grammar in both the scanner and the executor, including `Field[i].x`,
  `View[i].x`, and indexed members. Storage-kind transforms were generated
  mechanically from the vendored option callbacks (252 entries): `(int)val`
  truncation, `(unsigned)` clamps, `val ? 1 : 0` bool stores, and
  `max(trunc(val),1)` floors all apply at write time, so
  `Geometry.OldNewReg = 0.7` reads back as per-dim mode exactly like upstream.
- `Geometry.OldNewReg` (default 1) selects the shared `NEWREG()` counter for
  `news`/`newl`/`newll`/`newsl`/`newv`; `= 0` reads each symbol's own
  dimension counter, matching the lexer bindings at `Gmsh.l`. `newreg`
  always reads the shared allocator. Per-dimension reads take the
  cross-kernel maximum when OCC internals exist, incrementing each kernel
  separately before the max so the Int32 wrap order matches
  `NEWCURVE()`/`NEWPOINT()`/`NEWVOLUME()` in `Gmsh.y`.
- `SetMaxTag` applies to the active factory kernel for dims -2..3
  (`GeoEntity{-1}`/`{-2}` are the curve-loop and surface-loop counters),
  emits the non-aborting "out of range" `yymsg` diagnostic for those loop
  dims while still applying — verified against the binary — and genuinely
  out-of-range dims error without applying. `SetFactory` synchronizes all
  counters across dims -2..3 in both directions, matching `Gmsh.y:302-315`.
- Model-lifecycle counters: `NewModel` (`new GModel`) and
  `Delete Model`/`Delete All` (`destroy` → `_freeAll`+`_allocateAll`)
  re-initialize every entity counter from `Geometry.FirstEntityTag - 1` and
  the physical counter from `Geometry.FirstPhysicalTag - 1`, reading the
  options at reset time — `Geometry.FirstEntityTag = 5; NewModel` yields
  `newp == 5` on both sides. `Delete Model` keeps the OCC internals object
  alive (`Box` still works, verified against the binary); `NewModel` and
  `Delete All` create a fresh `GModel` with null OCC internals, so
  OCC-gated statements fall back to the built-in diagnostic.
  `Delete Physicals` matches `resetPhysicalGroups`: memberships and raw
  groups clear while the physical counter and name bindings survive.
- Factory gating parity: `Box`/`Cylinder`/`Cone`/`Torus` and the 4-7
  parameter `Sphere` forms are OCC-only upstream; under the built-in
  factory they emit "… only available with OpenCASCADE geometry kernel",
  create nothing, and consume no allocator counters (the scan-side
  observer mirrors the same gate). The two-point `Sphere`
  (`newGeometrySphere`) and `PolarSphere` stay built-in under either
  factory. Booleans reproduce the grammar's asymmetry: tagged
  `BooleanUnion(3) = …` statements silently no-op under built-in while
  standalone `BooleanUnion{…}` terms report the diagnostic and yield an
  empty list, skipping operand resolution entirely.
- `Normal Surface {tag} Parametric {u,v}` and
  `Parametric Point {point} In Surface {surface}` evaluation cover ruled
  and triangular surface kinds with upstream's forward/backward finite
  differences; remaining digits-level residuals are fd noise and the
  LAPACK-vendor `uv error` bound already documented.
- Test fixtures previously encoding pre-parity permissive behavior were
  corrected to Gmsh-faithful sources: `.geo` `Box`/`Sphere`/Boolean cases
  now open with `SetFactory("OpenCASCADE")`, a five-component `Point`
  `VExpr` and three-component `Symmetry` plane are accepted as legal, and
  diagnostics match upstream's exact text ("Unknown action", "syntax
  error"). Eighteen CRC pins across `api_test.jl`/`cli_test.jl`/
  `geo_dynamic_tag_test.jl`/`geo_geometry_expression_test.jl`/
  `geo_mesh_size_test.jl`/`geo_set_max_tag_test.jl` were re-measured and
  proven stale by re-running every fixture on the pre-change commit
  (`2f1b16a`), which produced bit-identical SHAs — the drift predates this
  increment.

Re-verified on 2026-09-19 with Julia 1.12.7 against the pinned Gmsh 4.15.2
binary and vendored 4.15.2 source after a physical-group/allocator parity audit
of the `.geo` scanner and executor:

- Raw physical declarations now live in a parser registry while observable
  groups derive from signed entity-level memberships at gated synchronization
  points (`GEO_Internals::_changed` is mirrored by an executor `geo_changed`
  flag): `SyncModel`, `BoundingBox`, `Mesh`, `Save`/`Print`/`Merge`,
  `Delete`, `Show`/`Hide`/`Color`, `Boundary`/`PointsOf` actions,
  `Physical{}` selectors, `SetFactory`, `RelocateMesh`/`ReorientMesh`/
  `ClassifySurfaces`/`AdaptMesh`/`RefineMesh`, and expression-level reads
  all match the grammar's sync sites; `Include` does not sync.
- Compound `+=`/`-=`/`*=`/`/=` semantics verified statement-for-statement
  against `Gmsh.y`/`GModelIO_GEO.cpp`: `+=` appends without dedup on existing
  groups and is a recoverable error on missing ones, `-=` deletes emptied raw
  groups but leaves entity memberships of a negative-tag group intact through
  Gmsh's `abs(p) != tag` filter (the `-4` quirk reproduced bit-for-bit), and
  `*=`/`/=` are recoverable errors.
- Signed and zero tags: literal `Physical Point(0)`/`(-4)` tags, `(int)`
  truncation for fractional tags and members (`{1.9}`→1, `{-2.9}`→-2 with the
  negative sign preserved), member `0` resolving entity `0` to a `sign(0)`
  membership, `("name",0)` resolving `getMaxPhysicalNumber` from entity-level
  physicals (empty pre-sync state yields tag 1, verified), name-only
  declarations bumping `setMaxPhysicalTag(maxTag+1)` unconditionally with
  Int32 wrap to `-2147483648`, and the signed-subtraction `ComparePhysicalGroup`
  wrap ordering all reproduce Gmsh exactly.
- Lifecycle parity: `NewModel` clears names/internals but keeps parser
  symbols, `Delete Model` preserves `_physicalNames`/`_elementaryNames`,
  `Delete All` additionally clears function/string symbol tables, and
  `Delete Physicals` clears memberships while keeping name bindings.
- Recoverable diagnostics: duplicate field declarations, unknown/undeclared
  field options, multi-value `Background Field`, and unresolvable
  `BoundaryLayer Field` ids record `scan_errors`/`scan_warnings`/
  `scan_msg_error_count` on `GeoParams` instead of aborting; range/resource
  and syntax errors still throw.
- Two real defects fixed: the ten-field positional `GeoParams` constructor was
  missing after the diagnostic fields were added (restored), and
  `remove_physical_groups!` erased groupless-tag name bindings without
  counting them, so the API mesh cache stayed stale — a selected tag now
  counts when it carried a group or a name binding, matching
  `GModel::removePhysicalGroup`'s unconditional `_physicalNames` erase.
- `model_set_tag!` keeps its deliberate positive-tag contract at the public
  layer even though `GModel::changeEntityTag` performs no sign validation
  (probed `setTag(0,1,0)` retags to `(0,0)`): `.geo` reserves nonpositive
  tags for auto-allocation, so `model.set_tag` rejecting tag 0 is a bounded
  divergence the `model_entity_identity` differential asserts explicitly.
  Tag-0 entities themselves remain reachable through `.geo` literals
  (`Point(0)`). A `.geo`-deleted entity tag resurrected by a later
  definition regains its physical memberships through tag resolution —
  verified in Gmsh's own MSH output (`Sphere(1)` inherits physical 10 after
  `BooleanDifference` deletes box 1).
- Closing gates on the final tree (Julia 1.12.7): first full-suite run found
  423268/423281 with 12 failures + 1 error, every one root-caused to a stale
  pre-parity expectation or the two defects above. The `run_all.jl`
  `model_entity_identity` differential then caught that the tag-0
  `set_tag`/`entity_name` relaxations had broken a deliberate divergence —
  restored to `>0` at the public layer (tag-0 entities still reachable via
  `.geo` literals, which bypass the public validator). Re-verified
  file-by-file —
  `entity model validation` 130/130, `synchronized API entity identity
  lifecycle` 41/41, `owned Physical-group lifecycle through API` 66/66,
  `bounded .geo Boolean snapshots` 13/13, `bounded .geo geometry expressions`
  49/49, `Gmsh-compatible size fields` 6976/6976, `.geo` dynamic-tag 123/123,
  IO 433/433, removal 105/105, constraints 131/131, OCC primitives 110/110 —
  and `validation/gmsh_parity/geo_dynamic_tags.jl` passed end-to-end with
  identical CRCs against the pinned binary. The closing full-suite re-run on
  the final tree (after the `>0` revert) passed 423358/423358 in 153m30.6s
  with zero failures and zero errors, and `validation/run_all.jl` completed
  end-to-end: all 58 differential scripts green against the pinned Gmsh
  4.15.2 binary (including `GMSH_PARITY_MODEL_IDENTITY_OK` after the
  revert), all five primitive case volumes exact, the ASCENT enclosure/coax
  acceptance case reproducing Gmsh's documented empty-volume failure, and
  `validation/REPORT.md` regenerated.

Re-measured on 2026-09-17 with Julia 1.12.7 after completing the built-in-kernel
spline family (`Spline`/`BSpline`/`Bezier`/`Nurbs`):

- `src/geometry/ModelSplines.jl` adds all four curve kinds, every one
  reporting `getType` `"Nurb"` like Gmsh. `Spline` ports `InterpolateCurve`'s
  Catmull–Rom branch (ghost-endpoint extrapolation, cyclic first==last
  handling); `BSpline` ports the uniform `InterpolateUBS` matrix branches;
  `Bezier` ports De Casteljau; `Nurbs` ports `findSpan`/`basisFuns`/
  `InterpolateNurbs` with Float32-rounded knots, degree inferred as
  `nknots − npts − 1`, and the raw first/last knot values as the parameter
  interval. `Order` is parsed and ignored, and empty `Knots {}` records
  convert to `BSpline` — all matching `GEO_Internals::addBSpline`.
- The evaluators are bit-exact against the pinned Gmsh 4.15.2 binary:
  the compiled-C++ FMA contraction patterns were recovered by disassembling
  `libgmsh.dylib` (Catmull–Rom coefficients and accumulation, UBS `vec[i]`
  fma chains with plain mul+add `T·vec` reduction, De Casteljau
  `fma(w,c1,u*c2)`, `basisFuns` `fmadd`, NURBS fused `fmla` accumulation)
  and reproduced in Julia.
- Derivatives, curvature, bounds, bounding boxes, projection, and
  containment port the Gmsh code paths: `getDerivative`/`getSecondDerivative`
  use the same 1e-8 forward/central finite differences (`point(u±eps)`),
  `getCurvature` the `‖d1×d2‖·‖d1‖^-3` formula, `parFromPoint` the
  coarse-scan + `XYZToU` Newton refine (with `std::min(uMax,std::max(uMin,·))`
  NaN clamp semantics), `getClosestPoint` the golden-section search,
  `getBoundingBox` the 10-sample scan with `SBoundingBox3d`'s
  skip-non-finite comparison semantics, and `isInside` physical/parametric
  modes (parametric uses the raw knot interval for Nurbs).
- `.geo` execution handles `Spline`/`BSpline`/`Bezier`/`Nurbs` including
  the full `ListOfDouble` grammar at both list positions: literal `{...}`,
  list variables `pts[]`/`kk[]`, indexed lists, `-{...}`, `expr*{...}`,
  and scalar expressions, via a new top-level keyword splitter
  (`_geo_split_at_keyword`). `Order <expr>` is evaluated for validity and
  ignored semantically.
- Deliberate divergences where Gmsh itself is undefined: inferred NURBS
  degrees above `npts−1` are rejected at construction (Gmsh's `findSpan`
  then indexes control points out of bounds through `List_Read`), and a
  reversed Nurbs in a surface loop is record-only — Gmsh 4.15.2 segfaults
  inside `findSpan` on that path, so there is no upstream parity target.
- Duplication semantics: GEO `Duplicata` on a negative loop member
  materializes a positive-tag copy of the `CreateReversedCurve` record
  (knots index-mirrored, bounds `1−uend`/`1−ubeg`), while OCC
  `BRepBuilderAPI_Copy` preserves signed members and shares memoized
  copied edges — `_duplicate_surface!` now branches on the caller.
- Focused differential `validation/geo_splines/differential.jl` compares
  Tessella against the pinned binary on 20 cases (open/closed
  Spline/BSpline/Bezier/Nurbs, variable-list grammar, negated/multiplied
  entity lists, `Duplicata` reversal) across `getType`, `getValue`,
  `getDerivative`, `getSecondDerivative`, `getCurvature`, bounds, bounding
  boxes, `getParametrization`, `getClosestPoint`, `isInside`, and oriented
  surface boundaries — `GEO_SPLINES_DIFFERENTIAL_OK gmsh=4.15.2-git
  cases=20 samples=1130`, bit-for-bit (NaN positions compared
  positionally).
- Closing gates on the final tree (Julia 1.12.7): `Pkg.test()`
  423161 passed / 0 failed / 0 errored (65m45s); `validation/run_all.jl`
  re-verified end-to-end with every section reporting OK against the
  pinned Gmsh 4.15.2 and `validation/REPORT.md` byte-identical. An earlier
  full-suite run on the same milestone tree caught a
  `_duplicate_surface!` regression where the GEO reversal convention had
  been applied to the OCC copy path (423156/4/1, all five non-passes in
  OCC primitive lifecycle); the caller-split fix re-passed
  `occ_primitives_test.jl` 219/219 standalone before the closing run.
  `git diff --check` clean.

Re-measured on 2026-09-16 with Julia 1.12.7 after completing the N-way
Boolean milestone and adding `.geo` `Function`/`Call` execution:

- `boolean_volumes_multi!` and the `.geo` `BooleanDifference/Union/
  Intersection/Fragments` statements accept multi-entity `Volume{a,b,…}`
  object and tool operand lists with optional per-group `Delete` suffixes,
  empty tool groups, and explicit result tags. The arrangement decomposes
  into membership-labeled cells (UInt64 masks, ≤62 operands), materializing
  one result volume per disconnected kept cell with partition faces shared
  between adjacent cells. Gmsh 4.15.2 preserve-numbering is reproduced:
  untouched operands stay pseudo-preserved in place, unique sole-image
  modifications rebind the freed operand tag, shared-image and multi-source
  pieces allocate fresh tags after the surviving maximum, explicit tags
  reject multi-piece results, and empty results bind nothing. `v[]`/`v()`
  captures and arbitrary list variables hold the output list for later
  statements. Differentials against the pinned binary on the box battery —
  cut3/frag3/int3/frag2 and containment variants — match output tags, piece
  order, and entity counts exactly (frag3: 24 points, 44 curves, 26
  surfaces, 5 volumes), including shared cavity-wall reuse when a result
  piece borders a preserved operand.
- Boolean operand records store `(op, operands)` plus `(meshes, cell)`
  operation-time snapshots; bounds, transforms, identity/retag, removal,
  topology, and spatial-query paths consume both the legacy binary tuple
  shape and the N-way record. The `.geo` tag allocator resynchronizes from
  the materialized model after Boolean statements, and its predictive
  observer only tracks simple single-operand forms.
- `.geo` `Function name ... Return`/`Call name;` run Gmsh's zero-argument
  function semantics: token-level body capture, shared variable scope,
  `Call` inside `If`/`For`/`While`, quoted string names, Gmsh's
  `Unknown function`/`Redefinition of function` errors, and a bounded
  128-deep call stack. Verified against pinned Gmsh 4.15.2 probes and
  `geo_control_flow_test.jl` (25 focused tests).
- Full suite on the Boolean-milestone tree: 423023 passed / 0 failed /
  0 errored (Julia 1.12.7, ~66 min). The Function/Call feature landed
  after that run started.
- Closing run on the combined tree (2026-09-17, Julia 1.12.7):
  `julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'` —
  423051 passed / 0 failed / 0 errored (~82 min, concurrent with the
  validation driver). The aggregate `validation/run_all.jl` gate
  re-verified end-to-end: every section reported OK against the pinned
  Gmsh 4.15.2, including `GMSH_PARITY_BOOLEAN_OK` (tessella_volume=1,
  retained_snapshot=1, deleted_reuse=1, tessella_tets=12) and the
  documented `06_enclosure_coax` gmsh-failure/Tessella-pass row; the
  rewritten `validation/REPORT.md` is byte-identical. The aggregate run
  exposed one regression the unit suite missed: Boolean-materialized
  result vertices carried `_add_occ_point!`'s placeholder `point_size`,
  which fed `_volume_boundary_size_field` and silently refined the
  result mesh (24 tets, V=0.9999999999999996 instead of 12 tets, V=1.0).
  Both materializers now scrub the placeholder like every other OCC
  materializer; `model_boolean_multi_test.jl` gained a regression test
  (empty `point_size`, 12-tet unit difference). Focused re-verification:
  `model_boolean_multi_test.jl` 67/67, `geo_control_flow_test.jl` 85/85
  (incl. 25 Function/Call), `model_topology_query_test.jl` 828/828 (incl.
  742 materialized-Boolean pcurve/incidence assertions), plus
  dynamic-tag/expression/transform/identity/removal/model_test files
  green. README.md updated to the materialized-Boundary contract.
  `git diff --check` clean.

Re-measured on 2026-09-15 with Julia 1.12.7 after implementing rotational
`Extrude` (Gmsh's built-in-kernel `revolve`) for `.geo` and the model API:

- `revolve_entities!` and the `.geo` motion group `Extrude {{axis},
  {point}, angle} {..}` rotate-extrude points, curves, and planar surfaces
  with `ExtrudeShapes(ROTATE, ...)` semantics: swept vertices produce
  `Circle` arcs wired `[start, axis-center, end]` (the center is the
  source's orthogonal axis projection, allocated after the arc tag like
  `DuplicateVertex` inside `ExtrudePoint`, with the `Circle.n=(0,0,1)`
  `EndCurve` fallback seed), curve generatrices produce ruled/triangular
  laterals, surface generatrices produce volumes, and signed generatrices
  follow the reversed-record rules. On-axis, zero-angle, and full-turn
  collapses return the source tag with the unmerged copies left behind —
  `[src]` in `out` — matching Gmsh's skip-coherence path; `out` otherwise
  appends the top/body pair plus laterals under
  `Geometry.ExtrudeReturnLateralEntities`. OCC generatrix records transform
  with their copies through the shared rigid-transform helpers.
- The motion group is decoded by evaluated element shape (`[3,3,1]` revolve,
  `[3,3,3,1]` twist, `[1,1,1]`/single length-3 list translate) like the
  grammar's `VExpr`/`FExpr` alternatives; twist binds translation, axis
  direction, and point-on-axis in `GEO_Internals::twist`'s forwarded order,
  sweeping a `Spline` generatrix through `Geometry.ExtrudeSplinePoints`
  (default 5) rotate-translate steps, while mixed revolve+translate groups
  raise explicit `ArgumentError`s.
- `twist_entities!` and the `.geo` motion group `Extrude {{delta},
  {axis}, {point}, angle} {..}` run `ExtrudeShapes(TRANSLATE_ROTATE, ...)`
  semantics on points, curves, and planar surfaces: swept vertices produce
  `Spline` helices through `Geometry.ExtrudeSplinePoints` generated
  vertices — each `DuplicateVertex`-copied from the previous and stepped
  `angle/d` about the axis plus `delta/d` — while chapeau curves/surfaces
  take the full rotate-translate copy and endpoint generatrices extrude
  recursively. The 15-case `geo_twist` differential covers point, curve,
  and surface twists, `ExtrudeSplinePoints` = 1/3/5, mixed entity lists,
  lateral suppression, extrude parameters, and spline-knot evaluations:
  13/15 bit-exact, two arbitrary-axis cases inside the documented 64-ulp
  FMA-contraction band.
- Bit parity against the shipped 4.15.2 binary required per-site contraction
  matching: `SetRotationMatrix`'s Gram-Schmidt/`norme`/`prodve`/`prosca` and
  matrix products fuse (`fma`), while `vecmat4x4` application rounds each
  product-add separately — verified by extracting Gmsh's effective affine
  maps on probe points. `EndCurve`-equivalent arc helpers and the OCC record
  evaluators were brought to the same contract, including `ElSLib`
  `A1·X+A2·Y+A3·Z+P` associations, `semiAngle = atan((r2-r1)/h)` cone
  evaluation, and the OCC620 torus epsilon clamp.
- Volume meshing through revolved non-planar laterals now fails explicitly
  via the straight-curve gate in `_model_planar_surface_mesh` instead of
  reaching CDT internals.
- Verified in `geo_extrude/differential.jl` (26 cases/350 samples) against
  Gmsh 4.15.2: entity tags per dimension, point coordinates, curve/surface/
  volume boundary wiring, and `out[]` result lists match bit-for-bit across
  point/curve/surface revolves on and off the axis origin, negative
  generatrices and angles, endpoint-on-axis collapses, arc generatrices,
  signed surfaces, on-axis/zero-angle/full-turn collapses, >π arcs,
  `ExtrudeReturnLateralEntities=0`, non-axis-aligned axes, and mixed
  translate+revolve allocation interleaving.
  `geo_extrude_test.jl` covers the API surface, OCC-generatrix records,
  closed-curve generatrices, degenerate paths, and allocator plumbing
  (142 tests); `geo_transform_test.jl`, `geo_curved_test.jl`, and
  `occ_primitives_test.jl` all pass.

Re-measured on 2026-09-14 with Julia 1.12.7 after materializing the Torus
primitive and its analytic OCC records:

- `add_torus!` and the `.geo` `Torus(id) = {x,y,z,r1,r2[,angle]};` statement
  build `BRepPrimAPI_MakeTorus`'s exact entity layout: a full torus is one
  rim Point, a closed outer-equator Circle, a closed meridian Circle, and a
  `Torus` face wired `[-equator,+meridian,+equator,-meridian]` under shell
  `[face]`; a partial torus adds the second rim vertex, trims the equator to
  `[0,angle]`, closes both end meridians, and caps the ends with `Plane`
  faces under shell `[torus,+start_cap,-end_cap]`. Every face and edge owns
  an analytic record, so entity, boundary, type, evaluation, derivative,
  parametrization-bound, and bounding-box queries answer through the
  materialized topology; the (possibly partial) face box is an exact
  sweep maximization including spindle-torus stationary roots.
- Whole-solid similarities rewrite the stored frames atomically; a
  reflection flips the torus axis and each trimmed arc's circle normal so
  `p'(t) = T·p(t)` keeps `t0` on the start vertex, closed circles keep
  `+T·n` and negate their range, and anisotropic dilations reject
  explicitly. Independently moved rim vertices leave records that fail
  queries via the new endpoint-satisfaction check instead of answering
  with stale geometry.
- Verified against Gmsh 4.15.2 OCC in `geo_primitives/differential.jl`
  extended to 15 cases/307 samples: entity tags per dimension, point
  coordinates, curve/surface type strings, signed boundary wiring,
  parametrization bounds, and curve/surface evaluations match bit-for-bit
  or under 1e-12 for full, partial, explicit-2π, off-center, and spindle
  toruses; OCC's coarse polyhedral torus box is checked for containment.
  `occ_primitives_test.jl` covers materialization, transforms, reflection,
  stale-record failure, shared-boundary duplication, retagging, and
  removal; `geo_dynamic_tag_test.jl` covers the full-vs-partial allocator
  counts.

Re-measured on 2026-09-14 with Julia 1.12.7 after implementing the legacy
`Mesh.TransfiniteTri=0` three-sided algorithm:

- `set_transfinite_surface` on a 3-curve loop now defaults to
  `mesh_transfinite_triangle_collapsed` — Gmsh's collapsed-quadrilateral
  `TransfiniteTri=0` path — with auto-rotated or explicitly pinned collapsed
  corners; `Mesh.TransfiniteTri=1` (via `set_transfinite_tri!`,
  `option("Mesh.TransfiniteTri", 1)`, or `.geo` `Mesh.TransfiniteTri = 1;`)
  selects the compact triangular lattice.
- Verified against Gmsh 4.15.2 (`Mesh.TransfiniteTri=0`): identical collapsed
  grids including the corner-rotation rule, chord-averaged interior placement
  (max node error ≈2.5e-15), fan-plus-cell element order, and all four
  diagonal arrangements; unequal-side boundaries like (5,5,8) rotate to the
  matching corner exactly as Gmsh's `findTransfiniteCorners` does.
- Four-sided transfinite surfaces now route through `mesh_transfinite_patch`,
  fixing a `Left`/`Right` diagonal inversion in the previous inline grid and
  adding the `AlternateLeft`/`AlternateRight` parity Gmsh emits.
- The `.geo` executor accepts `Mesh.TransfiniteTri = 0|1;` and reports it on
  `GeoExecution.transfinite_tri`; `open_geo!` propagates it into the session
  option like Gmsh's global option store.
- `api_test.jl`, `cli_test.jl`, and `transfinite_triangle_test.jl` updated and
  passing; `validation/transfinite_triangle/differential.jl` extended with a
  collapsed-algorithm section.
- The collapsed kernel's recombined arm is
  `mesh_transfinite_triangle_collapsed_patch` — the `Recombine`/`RecombineAll`
  arm of `Mesh.TransfiniteTri=0`, returning a `MixedMesh` with the apex fan
  kept triangular and one quadrangle per remaining grid cell. Emission is
  arrangement-independent (upstream's recombine branch precedes the diagonal
  dispatch), certified against the unrecombined kernel's geometry, audited
  for atomic-triangle coverage and boundary conservation, and differential-
  verified cell-for-cell against Gmsh 4.15.2 including corner rotation and
  pinned-corner boundaries.

Re-measured on 2026-09-14 with Julia 1.12.7 after materializing `add_box!`
boundary topology and generalizing planar surface meshing off z=0:

- `add_box!` now owns Gmsh `addBox`'s exact boundary representation: 8 corner
  Points, 12 edge Curves, 6 planar Surfaces, and one Surface Loop with Gmsh's
  oriented shell signs `[-1,2,-3,4,-5,6]` and entity numbering. Box volumes
  answer boundary, adjacency, spatial-query, and MSH-entity-classification
  queries like explicit shells; recursive removal descends through them;
  transforms resynchronize the stored corner Points; and the retained
  `box_extents` encoding is validated geometrically against the shell. Corner
  Points carry no explicit size so box meshes stay bit-identical to the former
  primitive path, matching Gmsh OCC corner sizing.
- Planar surface meshing, embedded-entity checks, projection classification,
  and `model_to_mixed` now derive each surface's coordinate-axis plane instead
  of assuming z=0: a box's six faces mesh standalone, and full-3D
  point-to-segment tests replace projected 2-D checks so out-of-plane nodes
  are rejected rather than silently flattened.
- Verified against Gmsh 4.15.2: identical `addBox` entity counts, shell signs,
  tag-allocation side effects (`newp`/`newl`/`news`/`newv` advance through the
  shared counter exactly as the `.geo` simulator reserves), `Point(1)`
  collision errors, `CombinedBoundary`/`PointsOf`/`MeshSize` resolution on box
  volumes, and non-recursive `removeEntities` retaining the 26 subentities.

Re-measured on 2026-09-14 with Julia 1.12.7 after wiring 3-sided transfinite
surfaces into the model attribute path:

- `set_transfinite_surface` on a 3-curve loop now routes through
  `mesh_transfinite_triangle` — the dedicated `Mesh.TransfiniteTri=1` patch —
  instead of failing; the legacy collapsed-quadrilateral `TransfiniteTri=0`
  algorithm remains a documented non-claim. The entity cache keeps the
  untagged simplex contract used by the four-sided Coons path.
- Verified against Gmsh 4.15.2 (`Mesh.TransfiniteTri=1`, three 5-point
  transfinite edges): identical 15-node triangular lattice including interior
  node placement; mismatched side counts and non-planar loops still fail
  explicitly, and reversed-sign loops mesh.
- Full suite `421,474/421,474` passed under bounds checking;
  `validation/run_all.jl` re-run after the change.

Re-measured on 2026-09-14 with Julia 1.12.7 after the hot-kernel allocation
audit:

- `orient2`/`orient3` now run Shewchuk expansion stages B–D
  (`PredicatesAdaptive.jl`): exact Float64 error-free transformations in
  per-thread scratch decide every case inside a conservative magnitude band
  with zero allocation; the exact dyadic BigInt path remains only for inputs
  outside the band. A new exact `diametral_sign` predicate replaces the
  `Rational{BigInt}` encroachment test in 2-D refinement, and the 3-D
  conformity gate's region side/pierce decisions use adaptive `orient3`/
  `orient2` instead of exact rationals, with the rational reference retained
  as the test oracle.
- `insert_point3!` reuses triangulation-owned cavity, epoch-stamped mark,
  boundary, stack, and spoke scratch instead of fresh containers per point.
  Measured on a 20,000-point `delaunay3d`: 24,687,168 bytes and 0.243 s versus
  111,840,800 bytes and 0.389 s before the change (−77.9% allocation, −37.5%
  time), bit-identical output.
- The 2-D Ruppert loop resumes its quality scan from the lowest slot that may
  have changed (`_newtri!` slot reuse and interior reclassification lower the
  bound) and caches the sorted constraint list, reproducing the full-scan
  insertion order exactly; a dedicated testset checks slot-for-slot parity
  with the uncached scan.
- Closure-boxing sweep across `MeshSurface`, `Optimize`, `Recombine`,
  `Mesh3D`, `Mesh2D`, `SizeField`, `SizeFieldCatalog`, `MeshPointLocation`,
  `HierarchicalBases`, `HigherOrderNodal`, `CAD`, `Heal`, `IO`,
  `BoundaryLayer`, `ExactMesh3D`, `TransfiniteHex`, `API`, and
  `ModelEntityEvaluation`: reassigned captures were rewritten as single
  assignments or explicit tuples, and the recursive local BVH builders were
  lifted to top-level builder records shared by the distance- and
  view-hierarchy paths.
- `test/core/allocation_audit_test.jl` gates the result: a lowered-code scan
  rejects any new `Core.Box` outside a documented cold-path allowlist, and
  per-call `@allocated` checks confirm the predicates, quality kernels, and
  conformity-gate decisions allocate nothing after warm-up at ordinary,
  exactly degenerate, tiny, and huge scales.

Re-measured on 2026-09-11 with Julia 1.12.7 after implementing entity-filtered
session mesh queries through a `model_to_mixed` classification snapshot stored
with the cache:

- `_MeshClassification` records per-node owning entities, per-cell entity tags
  for each simplex block, and the entity boundary map, built through canonical
  `model_to_mixed` at generation, re-derived on `refine`, rebound across
  `affine_transform` (connectivity is index-invariant), and dropped by every
  other cache replacement including model-invalidating mutations.
- `get_nodes`, `get_elements`, `get_element_types`, `get_elements_by_type`,
  `get_nodes_by_element_type`, `get_barycenters`, `get_element_edge_nodes`,
  `get_element_face_nodes`, `get_jacobians`, `get_basis_functions_orientation`,
  and `get_keys` now accept nonnegative entity tags: `dim=-1` ignores `tag`,
  `include_boundary` appends transitive boundary-entity nodes breadth-first
  after the entity's own (duplicates across the per-entity outer loop match
  Gmsh), filtered task slices partition the entity subset, and unknown or
  dimension-mismatched entities fail explicitly. New position-based kernels
  evaluate only the selected columns and keep hierarchical key element tags
  global.
- `get_element` resolves a dense element tag to `(element_type, node_tags,
  entity_dimension, entity_tag)` in Gmsh's result order; out-of-range tags and
  caches without classification fail explicitly. `get_node` likewise resolves a
  dense node tag to `(coordinates, parametric_coordinates, entity_dimension,
  entity_tag)`, reparametrizing the node on its owning entity.
  `get_nodes_for_physical_group` emits the sorted unique node set over each
  member's own, boundary, and transitively embedded entities;
  `get_embedded` reports stored entity embeddings; `get_sizes` returns Point
  mesh sizes with Gmsh-compatible zeros for other and unknown entities.
- Entity-selective `clear`, `affine_transform`, `create_edges`, and
  `create_faces` operate on classified nodes and cells; entities owning no
  cache cells no-op. `remove_elements` drops listed or all cells on an entity
  while retaining nodes; `reverse`/`reverse_elements` flip first-order simplex
  orientation with Gmsh's vertex conventions; `reorder_elements` permutes an
  entity's element block with Gmsh's zero-based source-position ordering;
  `set_node`, `renumber_nodes`, and `renumber_elements` apply validated
  coordinate and dense-tag updates; `remove_embedded` drops embedding records;
  `get_duplicate_nodes` scans owned nodes for exact-coordinate duplicates and
  `remove_duplicate_nodes`/`remove_duplicate_elements` merge or drop them.
  `get_periodic` reports each entity's periodic master (or itself),
  `compute_renumbering` returns a reverse Cuthill-McKee node renumbering over
  the shared-element adjacency graph, `optimize` runs the validated
  boundary-preserving tetrahedral optimizer on the cache, `remove_constraints`
  is a validated no-op matching Gmsh 4.15.2's per-entity-attribute scope, and
  `set_visibility`/`get_visibility` track raw per-element display state.
  Primitive `add_box`
  models expose no boundary entities, so their boundary queries correctly
  reject. The session owns a Gmsh-parity multi-model list: `model.add` appends
  a fresh named model and selects it (duplicate names allowed), `set_current`
  selects the first matching slot and restores its geometry, mesh, and
  visibility state, `remove` deletes the current slot and selects the last
  remaining one, and `get_file_name`/`set_file_name` plus
  `set_visibility_per_window` track per-model state.
- The mesh-data differential now generates the same square in both engines and
  checks own-node-first ordering, transitive boundary closure against a
  recursive `getBoundary` walk, filtered type funnels, filtered task unions,
  `get_element` parity, `get_node` ownership/parametrization parity,
  Physical-group closure, embedding, size, `remove_elements`, `reverse`,
  `reorder_elements`, `set_node`, renumbering, `remove_embedded`,
  `remove_constraints`, RCMK `compute_renumbering` contract, and
  duplicate-node/element parity, and unknown/phantom-entity rejection:
  `entity_filtered=tris(4, 26)`.
- Focused bounds-checked suites: mesh-data 345/345 (including the new
  66-assertion entity-filtered testset), Jacobian 64/64, and function-space
  126/126. The complete bounds-checked package gate passed 212,533/212,533
  assertions and the aggregate validation exited zero with every mandatory
  child, analytic case, and the enclosure/coax acceptance probe completed.

Re-measured on 2026-09-11 with Julia 1.12.7 after extending hierarchical
function spaces from order-one H1 plus linear-simplex lowest-order H(curl) to
the full Gmsh 4.15.2 hierarchical contract: `H1LegendreN`/`GradH1LegendreN`
orders 1:15 on Line, Triangle, Tetrahedron, Quadrangle, Hexahedron, and Prism
reference families (Point stays order-independent), and
`HcurlLegendreN`/`CurlHcurlLegendreN` orders 0:11 on Line/Triangle/Tetrahedron
and 0:10 on Quadrangle/Hexahedron/Prism:

- The new `src/core/HierarchicalBases.jl` ports the pinned release's
  `HierarchicalBasisH1*`/`HierarchicalBasisHcurl*` generators verbatim,
  including Horner-form `OrthogonalPoly` polynomials, Solin edge/face
  conventions, all orientation sign tables, and Gmsh's Float32 recurrence
  coefficients in the H(curl) non-simplex paths. `MeshFunctionSpaces` routes
  hierarchical names through them for counts, orientation-major evaluation,
  orientation indices, key dimensions, key information, and vertex/edge/face/
  bubble key catalogs laid out like `getKeys`; `get_keys` and
  `get_keys_for_element` lazily populate the global edge and face catalogs only
  for spaces whose basis owns those entities.
- The mesh-function-space differential now checks every supported
  family/order/space combination through the public API against Gmsh 4.15.2:
  basis values, gradients, and curls at asymmetric and seeded random reference
  points, all orientations for orientations ≤ 720 and a spread beyond,
  orientation counts and indices, actual- and explicit-order key counts, key
  information on fabricated keys, and the type-4 dense-key catalog layout.
  Every hierarchical comparison is strict (8e-14); the reported
  max_abs_difference=1.19e-6 comes from the pre-existing high-order nodal
  comparator covering Gmsh's measured 8.5e-7 P9-hexahedron Vandermonde
  artifact, not the hierarchical paths. Stream checksum
  `45af63358fafcbc3b914260d66b0b1ba61d198560f5fa7a4ddc707ae3ef19076`
  reproduced identically across two runs and inside the aggregate gate.
- Fixes applied during verification: prism quadrilateral-face H1 orientation
  used Lobatto indexing on kernel polynomials (kernel functions index 0:pf-2,
  not degrees 2:pf); the public key path referenced an undefined mesh variable
  for element-level queries; and hierarchical key-information queries now
  ignore submitted type-key values like the pinned release does. Test
  expectations asserting rejection of the now-supported spaces were updated
  and a new contract testset pins Gmsh-verified key counts, orientation
  counts, slice consistency, and order boundaries.
- Focused bounds-checked suites: core function-space 3,434/3,434 and
  session-API function-space 126/126. The complete bounds-checked package gate
  passed 211,822/211,822 assertions in 48m14.7s and the aggregate
  bounds-checked validation exited zero with every mandatory child, analytic
  case, and the enclosure/coax acceptance probe completed, regenerating the
  untracked `validation/REPORT.md`. `git diff --check` passed.
- Residual: Pyramid/Trihedron hierarchical spaces (undefined in the pinned
  release) and entity-filtered orientation/key results remain pending per
  PLAN.md.

Re-measured on 2026-09-11 with Julia 1.12.7 after extending order-one
`H1Legendre1`/`GradH1Legendre1` reference evaluation, counts, and key metadata
to every fixed Quadrangle, Hexahedron, and Prism type (at any Lagrange order)
and Point type 15:

- `MeshFunctionSpaces` attaches H1 order 1 to the reference family, not the
  input type's Lagrange order: values repeat the family's vertex functions once
  per orientation (24/40320/720 for Quadrangle/Hexahedron/Prism, 1 for Point),
  keys are the family vertex counts (4/8/6/1), and Point H1 key metadata reports
  order 0 since the Point vertex function is constant. Lowest-order H(curl)
  stays linear-simplex-only; Pyramid/Trihedron hierarchical spaces stay
  rejected. Simplex hierarchical paths are bit-identical: every legacy
  differential comparison reproduces its prior value.
- Three pinned-release divergences are documented, not copied. Hexahedron
  `getNumberOfOrientations` reads uninitialized memory (observed 1833382193
  across three probe processes, 1918128693 inside the validation process) while
  `getBasisFunctions` counts 8! every time; Tessella uses 8!. Point
  `GradH1Legendre1` answers `[1,0,0]` per point where the pinned release's own
  nodal gradient correctly answers zeros; Tessella returns the verified zero
  gradient. Pyramid hierarchical queries fail inside Gmsh 4.15.2 itself
  (`Unknown familyType 6`), so they remain explicit blockers on both sides.
- Allocation scales with the selected slice on two-point quadrangle blocks
  (seven-run `@allocated` medians, bit-identical repeats): full 24-orientation
  block 2,096 bytes against a 1,536-byte payload, single-orientation selection
  720 bytes. New ratchets require the selection to cost less than the full
  payload and the full block to stay within payload + 4,096 bytes, rejecting
  full-compute-then-slice patches.
- Focused bounds-checked suites: core function-space 3,406/3,406 (101 new
  non-simplex H1 assertions) and session-API function-space 133/133.
  Neighboring element-catalog, reference-geometry, Jacobian, and data suites
  re-run with the full gate below. The Gmsh 4.15.2 mesh-function-space
  differential passes with 14 new non-simplex H1 cases (full/selected bases,
  components, orientation counts, key counts, key information, pinned
  divergences, Gmsh-side pyramid rejection) and SHA-256
  `86d2308f75da9d99ea76f8479d264f85ceaade9fcd46e36b0d34211e15fad0c5`
  (stable across two runs; legacy sections byte-identical).
- The complete bounds-checked package gate passed 211,801/211,801 assertions in
  29m19.6s and the aggregate bounds-checked validation exited zero with every
  mandatory child, analytic case, and the enclosure/coax acceptance probe
  completed, regenerating the untracked `validation/REPORT.md`.
  `git diff --check`, layout, recursive ambiguity, public-documentation, and
  Markdown-link gates passed.
- Residual: non-simplex H(curl) spaces, higher-order hierarchical spaces, and
  entity-filtered orientation/key results remain pending per PLAN.md.

Re-measured on 2026-09-11 with Julia 1.12.7 after implementing deterministic
contiguous-block `task`/`num_tasks` partitioning for the seven detached
read-only mesh-data queries (`get_elements_by_type`, `get_barycenters`,
`get_element_edge_nodes`, `get_element_face_nodes`, `get_jacobians`,
`get_basis_functions_orientation`, `get_element_qualities`):

- The contract follows Gmsh 4.15.2's `begin=(task*count)/numTasks`,
  `end=((task+1)*count)/numTasks` blocks (`src/common/gmsh.cpp:2582-2583`),
  probed through the pinned Julia binding. Tessella returns exactly the computed
  slice (no zero padding, which would fabricate invalid zero tags); qualities
  partition over requested-tag positions with slice-scoped validation while the
  other six partition over cached type-block positions. `task>=num_tasks` is
  the silently-empty range; negative, zero-count, or non-integer task arguments
  still fail explicitly. Two pinned-release divergences are documented, not
  copied: `getElementQualities` indexes its request vector unguarded, so
  `task>=num_tasks` reads out of bounds (observed `Unknown element 0`), and
  hierarchical `getBasisFunctionsOrientation` indexes per-entity elements
  unguarded, so `task>=num_tasks` segfaults the pinned release (reproduced in
  isolation before the differential was restricted to `task<num_tasks` for
  those two queries). Int128 block products survive adversarial counts
  (`typemax(Int)-1` of `typemax(Int)` selects the last of five segments).
- Allocation scales with the slice on a 10,000-segment chain (five-run
  `@allocated` medians, all runs bit-identical): elements 245,888→131,200
  bytes, barycenters 245,824→131,136, edge nodes 163,968→82,048, qualities
  566,016→307,712; wall time 54μs→36.7μs (elements) and 1.00ms→0.43ms
  (qualities). New `4part<=3full` ratchets reject full-compute-then-slice
  patches; hardcoded N=5/n=3 splits reject strided `i%num_tasks==task`
  patches; exact lengths/strictly-positive tags reject zero-padded patches.
- Focused bounds-checked suites: mesh-data 204/204 plus 76/76 new partition
  assertions, element-quality 36/36 plus 19/19, Jacobian 64/64,
  function-space 111/111. The four touched Gmsh 4.15.2 differentials
  (mesh-data five-segment chain with zero-padding placement checks,
  Jacobian two-element blocks, quality five-tag requests, orientation 2/6/24
  exhaustive models) pass with all pre-existing SHA-256 checksums unchanged.
- The complete bounds-checked package gate passed 211,667/211,667 assertions
  in 16m04.4s. Aggregate bounds-checked validation exited 0 in 1,767 s with
  39 `*_OK` child sections, zero failure markers, and the standing
  enclosure/coax probe (pinned gmsh leaves all three solid volumes empty;
  Tessella natively meshes them). `git diff --check`, layout, recursive
  ambiguity, public-documentation, and Markdown-link gates passed.
- The organized tree contains 205 tracked Julia files and no repository-root
  `.jl`; all 75 organized `*_test.jl` files match 75 runner includes.
  Compatibility is exactly `1.12 - 1.12`; no Julia 1.11 test was run for this
  increment.

Re-measured on 2026-09-10 with Julia 1.12.7 after completing bounded `GaussN`
coverage for the Triangle, Tetrahedron, and Prism families:

- `MeshQuadrature` preserves every Gmsh 4.15.2 economical table: Triangle
  through order 20 and Tetrahedron through order 21. Higher `GaussN` orders use
  Gmsh's tensor transitions, and Prism composes the matching Triangle and Line
  rules. The simplex tables moved to `src/core/MeshQuadratureSimplex.jl` with
  symmetry-orbit expansion; `CompositeGaussN` still selects bounded
  tensor/Duffy/Gauss--Jacobi rules directly. Limits remain 128 points per axis
  and one million output points. Trihedra still fail explicitly.
- The quadrature-focused bounds-checked gate passed 36,914/36,914 assertions:
  analytic moments through Triangle/Prism order 20 and Tetrahedron order 21
  (2e-10 tolerance only at Triangle order 20 and Tetrahedron orders 10 and
  above, where the published decimals accumulate), exact point counts per order,
  even-order lattice aliasing, Gauss21/Gauss22-to-composite tensor-transition
  equality, detached lattice results, first-rejected-order bounds
  (`Gauss255`/`Gauss198`/`Gauss199`), and recursive ambiguity checks.
  Neighboring element-catalog, reference-geometry, function-space, and API
  quadrature sets passed 7,724/7,724.
- The complete bounds-checked package gate passed 211,536/211,536 assertions in
  15m14.8s. The mandatory Gmsh 4.15.2 differential matched 457
  coordinate/weight cases (Triangle/Prism `Gauss6` through `Gauss30`,
  Tetrahedron `Gauss6` through `Gauss29`, plus every fixed type and
  representative composite orders) with SHA-256
  `eed6c09d0cc9af974b030cb12ff9eba5892fa148dcc447f4ca8439fe04cbfeb1`.
- The organized tree contains 205 tracked Julia files and no repository-root
  `.jl`; all 75 organized `*_test.jl` files match 75 runner includes, and all 55
  source-include edges resolve and reach all 56 source files. Compatibility is
  exactly `1.12 - 1.12`; active configuration has no Julia 1.11 support or test
  target, and no Julia 1.11 test was run for this increment. `git diff --check`
  is clean.
- Residual: the aggregate `validation/run_all.jl` driver finished rewriting
  `validation/REPORT.md` on 2026-09-09 during the pre-shutdown session, but its
  exit status could not be re-verified after the PC shutdown; re-running the
  aggregate gate is the remaining step before the next increment.

Re-measured on 2026-09-09 with Julia 1.12.7 after implementing actual- and
explicit-order nodal reference functions:

- `MeshFunctionSpaces` and `API.mesh` now evaluate the actual complete or
  serendipity basis of every supported fixed type through unqualified
  Lagrange/isoparametric names. Canonical `Lagrange0` through `Lagrange10` and
  gradient names select the complete basis of the input family, subject to the
  fixed Gmsh 4.15.2 catalog. Numeric key queries reject interpolation nodes that
  the linear-simplex mesh cache does not own. Trihedron remains an explicit
  blocker.
- The focused core and session-API sets passed 3,294/3,294 and 97/97
  bounds-checked assertions. Across all 124 non-Trihedron fixed types they cover
  exact-node interpolation, partition and gradient sums, coordinate and Jacobian
  reproduction, explicit-order family selection, nodal metadata, malformed and
  unavailable orders, empty and detached results, concurrent first cache
  construction, and exact or near-node cancellation. Independent product,
  polynomial-mode, edge-trace, finite-difference, and apex-limit oracles cover
  the high-order Prism and incomplete Pyramid paths unavailable from Gmsh.
- A fixed-seed 372-case catalog audit measured worst errors of
  `1.970645868709653e-12` for partition, `1.4415135751733033e-11` for gradient
  partition, `7.595035711460696e-13` for coordinate reproduction,
  `5.970890448736554e-12` for coordinate-gradient reproduction, and
  `9.594913095156699e-7` against centered finite differences.
- The Gmsh 4.15.2 differential matched 384 supported actual-order and 16
  explicit-order cases in addition to the 30 basis, 30 orientation, 30 key, 248
  fixed-order-one, and 20 alias cases. Its maximum absolute difference was
  `1.190163588954641e-6`; SHA-256 is
  `32a6443fb95b938d0c32da3840847f2f3fac08f1ecf92a2d2b033f6397ed9fa7`.
- The complete bounds-checked package gate passed 186,642/186,642 assertions in
  50m52.4s (51m22.4s including package setup). Aggregate bounds-checked
  validation returned zero in 3h44m34.7s after every mandatory child, analytic
  case, and enclosure/coax acceptance probe completed; it rewrote the verified
  `validation/REPORT.md` content without changing the tracked file.
- Warmed 256-point order-nine gradient queries had seven-run median
  time/allocation pairs of 0.020284s/6,493,936 bytes for complete Hexahedron,
  0.061764833s/657,152 bytes for incomplete Hexahedron,
  0.03502775s/3,741,072 bytes for complete Prism,
  0.032282875s/508,896 bytes for incomplete Prism,
  0.8799345s/6,703,392 bytes for complete Pyramid, and
  0.008698416s/2,930,688 bytes for incomplete Pyramid. The complete order-nine
  Pyramid gradient at `prevfloat(1.0)` differed from its apex limit by at most
  `2.7000623958883807e-13`; its five-run one-point median was 0.263376375s and
  87,864,800 bytes.
- The organized tree contains 204 tracked Julia files and no repository-root
  `.jl`; all 75 organized `*_test.jl` files match 75 runner includes, and all 54
  source-include edges resolve and reach all 55 source files. Compatibility is
  exactly `1.12 - 1.12`; active configuration has no Julia 1.11 support or test
  target, and no Julia 1.11 test was run for this increment.

Re-measured on 2026-09-08 with Julia 1.12.7 after extending order-one nodal
reference functions to every supported fixed family:

- `MeshFunctionSpaces` and `API.mesh` now evaluate `Lagrange1` and
  `GradLagrange1` for all 124 fixed Point, Line, Triangle, Quadrangle,
  Tetrahedron, Hexahedron, Prism, and Pyramid catalog types. Unqualified
  Lagrange/isoparametric aliases cover Point and first-order types 1--7.
  Existing simplex H1/H(curl) orientation and node/edge-key behavior is
  preserved; unavailable interpolation orders, non-simplex hierarchical spaces,
  and the Gmsh-unsupported Trihedron basis fail explicitly.
- Independent checks cover vertex Kronecker interpolation, partition of unity,
  zero gradient sums, central finite-difference gradients, exact and near-apex
  pyramid behavior, exact-rational extreme-coordinate pyramid values and
  gradients, all fixed catalog dispatches, empty/detached results, metadata,
  bounds, and error paths. The focused core/API command passed 1,355/1,355
  assertions, and the neighboring catalog/reference-geometry/quadrature/API
  command passed 18,119/18,119.
- The complete bounds-checked package gate passed 184,606/184,606 assertions in
  59m37.9s. Direct and aggregate runs of the final Gmsh 4.15.2 differential
  matched 30 simplex basis, 30 orientation, 30 key, 248 fixed nodal, and 20
  linear-family alias cases. Each fixed nodal case used 20 points; the maximum
  absolute difference was `8.881784197001252e-16`, and SHA-256 remained
  `26d28400c434fa81836b6c0e83cf8afeabc8f812c86f6176d464ac7520816612`.
  Aggregate bounds-checked validation returned zero and wrote
  `validation/REPORT.md` after every mandatory child and acceptance probe ran.
- Written mutant checks reject wrong interpolation-order dispatch, vertex
  permutations, polynomial substitution for the rational pyramid, omitted prism
  or tensor factors, wrong tetrahedron edge order, parity-based orientation,
  point-major flattening, and accidental whole-cache edge creation.
- Warmed 100,000-point gradient evaluations had seven-run medians of
  0.036166291s for Hexahedron and 0.041223208s for Pyramid: 2.765 and 2.426
  million points/s, allocating 21,610,752 and 16,440,128 bytes (216.10752 and
  164.40128 bytes per point), respectively. The pyramid sample includes exact
  boundary factors that exercise its rational cancellation fallback.
- The organized tree contains 203 tracked Julia files and no repository-root
  `.jl`; all 75 organized tests exactly match 75 runner includes, and all 53
  source-include edges resolve and reach all 54 source files. Compatibility is
  exactly `1.12 - 1.12`; active configuration has no Julia 1.11 support or test
  target, and no Julia 1.11 test was run. Documentation, recursive ambiguity,
  new-placeholder, layout, and `git diff --check` gates passed.

Re-measured on 2026-09-08 with Julia 1.12.7 after extending bounded reference
quadrature to every fixed-node family with a Gmsh 4.15.2 rule:

- `MeshQuadrature` and `API.mesh.get_integration_points` now cover Point, Line,
  Triangle, Quadrangle, Tetrahedron, Hexahedron, Prism, and Pyramid across all
  124 supported fixed catalog types. Returned coordinate triples and weights are
  detached. Economical Gmsh tables, Cartesian tensors, simplex Duffy maps, prism
  composition, and the pyramid Gauss--Jacobi/Duffy construction all retain checked
  per-axis and total-point bounds. Trihedra and unavailable higher economical
  simplex/prism rules fail explicitly.
- Analytic moments passed for all eight families through economical orders zero to
  five and composite orders 0–6, 12, and 20, with the pinned anisotropic
  Quadrangle `Gauss2` exception checked against its exact Gmsh array. Boundary
  checks covered the 128-point order-255 line, 16,384-point order-254 triangle and
  order-255 quadrangle, and one-million-point order-197 tetrahedron, order-199
  hexahedron/pyramid, and order-198 prism rules, plus each first rejected order and
  allocation-before-preflight guard.
- The focused core/API gate passed 12,020/12,020 bounds-checked assertions;
  neighboring element-catalog, reference-geometry/Jacobian, and function-space
  sets passed 988/988. The complete bounds-checked package gate passed
  183,424/183,424 in 35m42.7s.
- The mandatory Gmsh 4.15.2 differential matched 377 coordinate/weight cases
  across all 124 supported fixed types and retained SHA-256
  `0ff395a225821835fad550c4388545faf45637f9485d5511e1ebc5c80e837562`.
  Aggregate bounds-checked validation returned zero after every mandatory child,
  analytic case, and the enclosure/coax acceptance probe completed and wrote
  `validation/REPORT.md`.
- Written mutant checks reject reversed tensor or prism nesting, omitted pyramid
  measure/scale factors, the wrong odd-order prism count, interpolation-order
  dispatch, and output allocation before family-specific preflight.
- The organized tree still contains 203 managed Julia files and no repository-root
  `.jl`; all 75 organized `*_test.jl` files match 75 runner includes, and all 53
  source-include edges resolve and reach all 54 source files. Julia compatibility
  is exactly `1.12 - 1.12`; the active tracked tree has no Julia 1.11 runtime,
  support, or test target, and no Julia 1.11 test was run for this increment.

Re-measured on 2026-09-08 with Julia 1.12.7 after adding bounded reference
quadrature for finite elements:

- `MeshQuadrature` and `API.mesh.get_integration_points` now return detached
  Gmsh-shaped coordinates and weights for all 50 fixed-node Point, Line, Triangle,
  and Tetrahedron catalog types. Economical `Gauss0` through `Gauss5` simplex
  tables are pinned to Gmsh 4.15.2; native Gauss--Legendre and Duffy construction
  supplies checked `CompositeGaussN` rules through the documented resource limits.
  Malformed names, excessive point counts, higher economical simplex orders, and
  non-simplex families fail explicitly before output allocation.
- Independent analytic moments passed through every requested degree for economical
  orders 0–5 and composite orders 0–6, 12, and 20. Boundary checks covered the
  128-point order-255 line, 16,384-point order-254 triangle, and 1,000,000-point
  order-197 tetrahedron rules, plus the first rejected order for each bound. The
  focused core/API gate passed 3,470/3,470 bounds-checked assertions; neighboring
  element-catalog, Jacobian, and function-space sets passed 988/988. The complete
  bounds-checked package gate passed 174,874/174,874 in 109m10.8s.
- Two direct runs and the aggregate run of the mandatory Gmsh 4.15.2 differential
  matched 160 coordinate/weight cases across all 50 supported fixed types and
  retained SHA-256
  `a7e167c24bde2a871b1c9e4e4e5ae6c6bbbec0675ca2c56eda0602760e2888c2`.
  Written mutant checks reject Cartesian simplex points, wrong point-count laws,
  transposed loop ordering, interpolation-order dispatch, shared output storage,
  omitted Duffy Jacobians, and permissive/overflowing rule parsing. Aggregate
  bounds-checked validation returned zero in 87m46.3s after every mandatory child,
  analytic case, and the enclosure/coax acceptance probe completed.
- Warmed order-20, order-40, and order-100 composite tetrahedron queries produced
  1,728, 10,648, and 140,608 points, allocated 63,936, 361,056, and 4,506,688
  bytes, and had five-run medians of 0.000014834, 0.000079083, and 0.001589708
  seconds, respectively.
- The organized tree contains 203 managed Julia files and no repository-root `.jl`;
  all 75 organized `*_test.jl` files match 75 runner includes, and all 53
  source-include edges resolve and reach all 54 source files. Julia compatibility
  remains exactly `1.12 - 1.12`; active configuration and code contain no Julia
  1.11 support or test target. Recursive ambiguity, public-documentation, tracked
  local Markdown-link, new-placeholder, layout, and `git diff --check` gates passed.

Re-measured on 2026-09-07 with Julia 1.12.7 after adding first-order
finite-element function spaces and degree-of-freedom keys:

- `MeshFunctionSpaces` and the synchronized `API.mesh` façade now evaluate the
  supported type-1, type-2, and type-4 Lagrange, H1, and lowest-order H(curl)
  bases and curls; report lexicographic element orientations; and return detached
  node or global-edge keys with stable coordinates and metadata. Edge keys lazily
  extend only the requested topology, reuse explicit identifiers, and avoid
  copying an already complete global catalog.
- The focused core and session-API sets passed 105/105 and 68/68 bounds-checked
  assertions. They cover aliases, multi-point and selected-orientation layout,
  all 2/6/24 simplex orientations, exact-cancellation fallbacks, sequential and
  manual edge identifiers, absent element types, detached results, cache rollback,
  and no-copy repeated element queries. The complete bounds-checked package gate
  passed 171,404/171,404 assertions in 46m22.6s.
- The mandatory Gmsh 4.15.2 differential matched 30 basis, 30 orientation, and
  30 key cases; exhaustively checked 32 hierarchical orientation vectors; and
  retained SHA-256
  `faa6db6b9191fd1722bfa46768d9959ba4d1a237b1e2670d547a08e6abad7c91`.
  Its written mutant checks reject wrong tetrahedron edge order, parity-based
  orientation, point-major flattening, and accidental whole-cache edge creation.
  Aggregate bounds-checked validation returned zero after this child, every other
  mandatory Gmsh child, all analytic cases, and the enclosure/coax acceptance
  probe completed and wrote `validation/REPORT.md`.
- Selected-orientation H(curl) evaluation for 10,000, 20,000, and 100,000 points
  allocated 1,688,352, 3,375,904, and 16,810,784 bytes; the five-run medians were
  0.003184875, 0.007372167, and 0.048284958 seconds. A repeated single-element key
  query allocated 1,152 bytes against catalogs containing 6,000, 60,000, and
  300,000 edges, with respective five-run medians of 1.459e-6, 1.583e-6, and
  1.375e-6 seconds.
- The organized tree contains 199 managed Julia files and no repository-root
  `.jl`; all 73 organized `*_test.jl` files match 73 runner includes, and all 52
  source-include edges resolve and reach all 53 source files. Julia compatibility
  remains exactly `1.12 - 1.12`; the active tracked tree and this increment contain
  no Julia 1.11 support or test target. Public-documentation, recursive ambiguity,
  tracked local Markdown-link, new-placeholder, layout, and `git diff --check`
  gates passed.

Re-measured on 2026-09-01 with Julia 1.12.7 after adding atomic explicit global
edge and face insertion:

- `MeshEntityTopology`, `API.mesh.add_edges`, and `API.mesh.add_faces` now attach
  positive caller-owned identifiers to node pairs, triangles, and quadrangles.
  Edge identifiers are unique within the edge catalog; face identifiers share one
  namespace across triangles and quadrangles. Exact association repeats are
  idempotent. Validation and copy-on-success replacement make every batch atomic.
  Existing `create_edges` and `create_faces` now preserve manual entries and fill
  missing simplex topology with unused deterministic identifiers.
- The focused core testsets passed 172/172 automatic and 88/88 manual assertions;
  the session API passed 196/196. The 2,843 neighboring mesh-type, mesh-data,
  lifecycle, affine-transform, point-location, Jacobian, element-quality, and public
  API assertions also passed. Coverage includes all 24 quadrangle permutations,
  tag collisions across automatic/manual and triangle/quadrangle namespaces,
  before/after-create insertion, detached output, fixed-seed randomized maps,
  cache invalidation, and rollback after malformed or conflicting batches.
- The Gmsh 4.15.2 differential matched automatic and manual edge, triangle, and
  quadrangle maps plus orientations before and after automatic completion. It also
  pins six deliberate safety divergences: Tessella rejects zero identifiers,
  repeated nodes, conflicting geometry/tag associations, reused tags, and partially
  invalid batches without mutation. The unchanged automatic SHA-256 is
  `7cb14f79831c1fb7f3420e34759c16513f64e4d8a5cad82e89bd4e71507f57d0`;
  the manual SHA-256 is
  `00110d2815d99ced60d72a8344958d0f0797fc5b85c938142dc4061f0abc8b06`.
- Batched insertion of 10,000, 20,000, and 100,000 entries allocated, respectively,
  754,272, 1,507,936, and 10,027,616 bytes for edges and 1,049,328, 2,097,904,
  and 13,730,544 bytes for quadrangles. Five-run 100,000-entry medians were
  0.012734708 seconds for edges and 0.016819292 seconds for quadrangles.
- The full bounds-checked package gate passed 171,231/171,231 assertions in
  23m57.6s. Aggregate bounds-checked validation returned zero after the updated
  topology differential, every other mandatory Gmsh child, all analytic cases, and
  the enclosure/coax acceptance probe completed and wrote `validation/REPORT.md`.
- The tree still contains 195 managed Julia files and no repository-root `.jl`;
  all 71 organized `*_test.jl` files match 71 runner includes, and all 51 source
  include edges resolve and reach all 52 source files. Julia compatibility remains
  exactly `1.12 - 1.12`, with no Julia 1.11 test or active target. Public
  documentation, recursive ambiguity, tracked local Markdown-link, placeholder,
  duplicate-pattern, layout, and `git diff --check` gates passed.

Re-measured on 2026-09-01 with Julia 1.12.7 after adding cached global simplex
edge and face topology:

- `MeshEntityTopology` and `API.mesh.create_edges`, `create_faces`, `get_edges`,
  `get_faces`, `get_all_edges`, and `get_all_faces` now own deterministic
  whole-cache catalogs and return detached query arrays. Tags follow first encounter
  across segments, triangles, then tetrahedra; edge orientation is `+1` for ascending
  node tags and `-1` for descending tags. Triangular-face orientation remains zero,
  matching Gmsh 4.15.2. Refinement, affine transformation, replacement, and clearing
  all invalidate both catalogs.
- The focused core and session-API sets passed 172/172 and 101/101 bounds-checked
  assertions. Neighboring mesh-type, mesh-data, lifecycle, transform,
  point-location, Jacobian, quality, and public-API sets passed all 3,116
  assertions. The tests cover exact insertion order, all triangular permutations,
  reversed edges, cache lifecycle and detachment, degenerate cells, invalid input,
  fixed-seed randomized topology, and allocation growth.
- The required Gmsh 4.15.2 differential matched nine edges, seven triangular
  faces, eight edge queries, and eight face queries after normalizing Gmsh's
  unordered all-entity output by tag. Its fixed SHA-256 is
  `7cb14f79831c1fb7f3420e34759c16513f64e4d8a5cad82e89bd4e71507f57d0`.
  Tessella deliberately rejects malformed partial node groups that Gmsh silently
  truncates. Entity-selective catalog creation adds only the cells classified
  on the listed entities since the cache owns classification metadata.
- Combined edge-and-face creation for 10,000, 20,000, and 100,000 disjoint
  tetrahedra allocated 4,588,000, 9,142,752, and 49,447,392 bytes. The
  100,000-tetrahedron five-run median was 0.041943292 seconds, confirming linear
  work and output-size-proportional storage.
- The full bounds-checked package gate passed 171,048/171,048 assertions in
  38m24.8s. Aggregate bounds-checked validation completed through its final report
  after the new topology child, every existing mandatory Gmsh child, and the
  enclosure/coax acceptance probe passed. The tree contains 195 managed Julia
  files and no repository-root `.jl`; all 71 organized `*_test.jl` files are
  matched by 71 runner includes, and all 51 source-include edges resolve and reach
  all 52 source files. Julia compatibility is exactly `1.12 - 1.12`, with no
  Julia 1.11 test or active target. Public-documentation, recursive ambiguity,
  source-include, tracked local Markdown-link, placeholder, duplicate-pattern,
  layout, and `git diff --check` gates passed.

Re-measured on 2026-08-31 with Julia 1.12.7 after adding cached linear-simplex
forward-map and Jacobian queries:

- `MeshReferenceGeometry`, `API.mesh.get_jacobians`, and
  `API.mesh.get_jacobian` now return detached, element-then-point Jacobians,
  determinants, and physical coordinates for cached type-1 segments, type-2
  triangles, and type-4 tetrahedra. Matrices use Gmsh's column layout;
  low-dimensional frames use compatible regularization, and tetrahedron
  determinants retain orientation. Exact-rational/`BigFloat` fallbacks reject
  degenerate, nonfinite, overflowed, and under-resolved results explicitly.
- The fixed-seed affine, focused robustness, and synchronized API sets passed
  416/416, 68/68, and 50/50 bounds-checked assertions. The neighboring point,
  quality, mesh-data, lifecycle, transform, and public-API sets passed all 1,180
  assertions. Malformed partial triples, cache mutation, output aliasing,
  out-of-order results, inverted orientation, extreme scales, and allocation
  growth are covered.
- The required Gmsh 4.15.2 differential matched six elements across all three
  simplex types, three bulk evaluation points, and four single-element requests.
  Its fixed SHA-256 is
  `85803b9ce40eed41e887957e73a58a1a3a1c1251a72ef3c1e5315b165772d631`.
  Tessella deliberately rejects incomplete/nonfinite coordinate triples and
  degenerate maps that Gmsh silently truncates or returns as zero maps.
- Warmed 10,000-, 20,000-, and 100,000-tetrahedron single-point batch queries
  allocated 1,048,848, 2,097,424, and 10,420,496 bytes. The 100,000-element
  five-run median was 0.019878541 seconds, confirming output-size-proportional
  allocation and linear work.
- The full bounds-checked package gate passed 170,775/170,775 assertions in
  15m28.6s. Aggregate bounds-checked validation completed through its final report
  after the new Jacobian child, every existing mandatory Gmsh child, and the
  enclosure/coax acceptance probe passed. The tree contains 191 managed Julia
  files and no repository-root `.jl`; the only subtree entrypoints outside domain
  folders are `src/Tessella.jl`, `test/runtests.jl`, and `validation/run_all.jl`.
  All 69 organized `*_test.jl` files are matched by 69 runner includes; all 50
  source-include edges resolve and reach all 51 source files. Julia compatibility is exactly `1.12 - 1.12`, with
  no active Julia 1.11 target. Public-documentation, recursive ambiguity,
  source-include, tracked local Markdown-link, placeholder, layout, and
  `git diff --check` gates passed.

Re-measured on 2026-08-31 with Julia 1.12.7 after adding cached linear-simplex
element-quality queries:

- `MeshElementQuality` and `API.mesh.get_element_qualities` now return the 13
  documented Gmsh 4.15.2 measures for dense linear triangles and tetrahedra, plus
  the compatible segment subset. Results preserve request order, duplicates,
  signed inverted-tetrahedron behavior, and cache ownership. Invalid requests are
  rejected before computation; degenerate shape measures return zero and undefined
  circumradii return `Inf`.
- The focused core and session-API sets passed 172/172 and 34/34 bounds-checked
  assertions. The adjacent cache-lifecycle, affine-transform, mesh-data,
  point-location, and public API sets passed 69/69, 92/92, 198/198, 76/76, and
  331/331 assertions. An additional fixed-seed affine probe passed 2,560
  shape-invariance and length-covariance checks across triangles, tetrahedra, and
  scales `2^-200` and `2^200`.
- The required Gmsh 4.15.2 differential matched all 13 measures on mixed,
  duplicate, positive, and inverted triangle/tetrahedron requests and all nine
  reliable segment measures. The four unreliable segment combinations are explicit
  blockers. Its fixed SHA-256 is
  `4b32e86fee56ef55e7ad571c640155d88f40ceaff857186768554acb4250f3e0`.
- Warmed 10,000-, 20,000-, and 100,000-tetrahedron batch queries allocated
  164,096, 327,936, and 1,605,888 bytes; the 100,000-element median was 0.02850
  seconds. The warmed per-element quality kernel allocated zero bytes.
- The full bounds-checked package gate passed 170,241/170,241 assertions in
  1,929.80 seconds. Aggregate bounds-checked validation, including the new quality
  child and enclosure/coax probe, exited 0 in 2,546.83 seconds. The tree contains
  187 managed Julia files and no repository-root `.jl`; the only subtree entrypoints
  outside domain folders are `src/Tessella.jl`, `test/runtests.jl`, and
  `validation/run_all.jl`. All 67 organized `*_test.jl` files are matched by 67
  runner includes. Julia compatibility is exactly `1.12 - 1.12`, with no active
  Julia 1.11 target. Public-documentation, recursive ambiguity, source-include,
  local Markdown-link, and `git diff --check` gates passed.

Re-measured on 2026-08-31 with Julia 1.12.7 after adding cached simplex
point-location queries:

- `API.mesh` now locates all or the deterministic first cached segment, triangle,
  or tetrahedron at a finite point and returns representable local coordinates for
  any nondegenerate dense element tag. A lazily reused AABB hierarchy is invalidated
  with the mesh cache. Exact-rational fallbacks diagnose degenerate maps and solve
  ill-conditioned ones. Float64-unrepresentable local coordinates fail explicitly.
  Strict and decade-relaxed reference tolerances follow the pinned Gmsh contract,
  while lower-dimensional off-span coordinates use stable orthogonal projection.
- The focused core sets passed 75/75 and 384/384 bounds-checked assertions; the
  fixed-seed affine oracle covers 64 segments, 64 triangles, and 64 tetrahedra.
  The session API set passed 76/76, including cache identity, invalidation,
  detachment, invalid-input preservation, degeneracy, and allocation scaling. The
  adjacent mesh-data, cache-lifecycle, affine, and public API sets passed 198/198,
  69/69, 92/92, and 331/331 assertions.
- The required Gmsh 4.15.2 differential matched six location and four local-coordinate
  cases, strict rejection, relaxed acceptance, and no-match errors. Its fixed SHA-256
  is `afdef11844fca534f32a62914a8c21432c282788c6c0504c74cfb2a2205dbdb7`.
  Orthogonal off-span segment and triangle results are the documented bounded
  difference from Gmsh's inversion artifacts.
- A seven-build, 1,000-query measurement over 100,000 segments recorded a 0.188 s
  median hierarchy build, 52.38256 retained index bytes per element, a 0.458 μs
  median strict query, and 256 bytes allocated per warmed query.
- The full bounds-checked package gate passed 170,035/170,035 assertions in
  1,087.73 seconds. Aggregate bounds-checked validation, including the new
  point-location child and enclosure/coax probe, exited 0 in 2,524.51 seconds. The
  tree contains 183 managed Julia files and no repository-root `.jl`; the only
  subtree entrypoints outside domain folders are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`. All 65 organized `*_test.jl`
  files are matched by 65 runner includes. Julia compatibility is exactly
  `1.12 - 1.12`, with no active Julia 1.11 target. Public-documentation, recursive
  ambiguity, source-include, local Markdown-link, and `git diff --check` gates passed.

Re-measured on 2026-08-31 with Julia 1.12.7 after adding fixed element catalog
queries:

- `msh_type` and `API.mesh.get_element_type` now resolve every fixed catalog
  family/order pair, including Gmsh's complete-type fallback for unavailable
  serendipity types. `msh_properties` and its API wrapper return detached names,
  dimensions, orders, node counts, local coordinates, and primary-node counts for
  all 125 fixed types. The compatible zero-node special metadata is exposed; other
  special records fail explicitly.
- The fixed-catalog, special/bounds, installed-oracle, and focused API sets passed
  1,643/1,643, 54/54, 1,213/1,213, and 281/281 bounds-checked assertions. The
  neighboring public API suite passed 331/331. The fixed API result checksum is
  `4c8e4086febda2bc3f482087d589eb157dbe57394f30781449a17a64e1b2a65c`.
- The required Gmsh 4.15.2 differential matched all 125 type lookups, all 110
  fixed types for which Gmsh returns properties, and seven compatible special
  property records. It confirmed 15 upstream property-call gaps; Tessella retains
  the catalog's separately checked high-order-prism and trihedron layouts for those
  types.
  Unknown family names are a deliberate explicit error instead of Gmsh's zero
  sentinel, and returned coordinate arrays do not alias catalog state.
- The full bounds-checked package gate passed 169,500/169,500 assertions in
  1,782.16 seconds. Aggregate bounds-checked validation, including the new element
  catalog child and enclosure/coax probe, exited 0 in 3,400.07 seconds. The tree
  contains 179 Julia files, no repository-root `.jl`, and 63 organized
  `*_test.jl` files matched by 63 runner includes. Julia compatibility is exactly
  `1.12 - 1.12`; the active configuration/runtime search found no Julia 1.11
  target. Public-documentation and recursive API ambiguity scans, local Markdown
  links, and `git diff --check` all passed.

Re-measured on 2026-08-31 with Julia 1.12.7 after adding connectivity-derived
mesh queries:

- `API.mesh` now returns detached per-element node coordinates, barycenters, and
  edge/face node arrays for cached linear segments, triangles, and tetrahedra.
  Shared nodes repeat in element order, local edge/face ordering matches Gmsh,
  and the primary-node selector is exact because all cached elements are linear.
- The focused mesh-data suite passed 198/198 bounds-checked assertions, including
  mutation isolation, invalid-input/cache preservation, refinement lifecycle,
  fixed derived-output SHA
  `04e09b72ebf17bdc7ab2f9f96da2927c5a6e892e5c2d98c9ddbb8313dc4cab13`,
  and linear allocation scaling. The public API, cached-lifecycle, and affine
  session suites passed 331/331, 69/69, and 92/92 assertions.
- The Gmsh 4.15.2 differential matched every type-1/2/4 node, coordinate,
  normalized/fast barycenter, edge, and triangular-face array for both primary
  settings. It also fixed the bounded differences: entity/task metadata remains
  unavailable, special MSH records remain distinct from fixed-node elements,
  invalid face-node counts are rejected instead of returning empty, and Tessella
  returns the finite true barycenter at `floatmax(Float64)` while rejecting an
  unrepresentable fast sum; Gmsh returned infinity for both calculations.
- The full bounds-checked package gate passed 168,423/168,423 assertions in
  1,436.13 seconds. Aggregate bounds-checked validation, including the extended
  mesh-data child and enclosure/coax probe, exited 0 in 2,492.67 seconds. The
  tree contains 177 Julia files, no root-level `.jl`, and 62 organized
  `*_test.jl` files matched by 62 runner includes. Julia compatibility is exactly
  `1.12 - 1.12`; the active-tree 1.11 support sweep, public mesh-doc scan,
  recursive API ambiguity scan, local Markdown links, and `git diff --check`
  all passed.

Re-measured on 2026-08-30 with Julia 1.12.7 after adding detached bulk mesh
queries:

- `API.mesh` now returns detached Gmsh-shaped node and linear-simplex element
  arrays, whole-dimension and type filters, and maximum dense tags. The session
  rebuilds node and element tags for each current cache and rejects entity
  classification, parametric coordinates, and nondefault task partitioning that
  the cache cannot represent.
- The focused mesh-data suite passed 110/110 bounds-checked assertions, including
  mutation-isolation, invalid-input/cache-preservation, affine/refinement/clear
  lifecycle, fixed refinement CRC
  `db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031`,
  and allocation-scaling ratchets. The full API, cached-lifecycle, and affine
  session suites passed 325/325, 69/69, and 92/92 assertions.
- The Gmsh 4.15.2 differential matched types 1/2/4, flattened coordinates and
  normalized connectivity, whole-dimension filters, absent and invalid type
  behavior, and detached results. It recorded Gmsh's explicit maxima 40/300
  against Tessella's dense maxima 4/3 and exercised the documented no-cache,
  classification, and task-partition blockers.
- The full bounds-checked package gate passed 168,329/168,329 assertions in
  857.18 seconds. Aggregate bounds-checked validation, including the new
  mesh-data child and enclosure/coax probe, exited 0 in 1,858.91 seconds. The
  tree contains 177 Julia files, no root-level `.jl`, and 62 organized
  `*_test.jl` files matched by 62 runner includes. Julia compatibility is exactly
  `1.12 - 1.12`; the active-tree 1.11 support sweep, public mesh-doc scan,
  recursive API ambiguity scan, local Markdown links, and `git diff --check`
  all passed.

Re-measured on 2026-08-30 with Julia 1.12.7 after adding whole-cache affine
transformation:

- `API.mesh.affine_transform` now normalizes strict 12-/16-entry Gmsh row-major
  input or a native 4×4 matrix through the canonical transform kernel. It commits
  only after validation, returns storage detached from the cache, preserves model
  geometry, and rejects entity-selective requests until classification metadata is
  available. Reflections rewind triangle/tetrahedron connectivity; singular,
  nonfinite, malformed, overflow-producing, and unsupported selections preserve the
  prior cache.
- The focused affine session suites passed 92/92 bounds-checked assertions. The
  neighboring direct transform, cached-refinement, and full API suites passed
  103/103, 69/69, and 315/315 assertions. The box connectivity CRC remained
  `e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83`;
  analytic transformed bounds were `((1,-2,0.5),(3,1,4.5))`, and clear/regenerate
  restored the original mesh.
- The required Gmsh 4.15.2 differential reported zero coordinate error for both a
  positive affine map and a reflection. It also re-measured the bounded differences:
  Gmsh retains reflected connectivity and accepts singular/surplus-entry input,
  while Tessella preserves positive simplex orientation and rejects those inputs.
- The full bounds-checked package gate passed 168,209/168,209 assertions in
  965.16 seconds. The aggregate differential/enclosure gate passed in 1,717.47
  seconds and included the new affine child. The tree contains 175 Julia files,
  no root-level `.jl`, and 61 organized `*_test.jl` files matched by 61 runner
  includes. Julia compatibility is exactly `1.12 - 1.12`; the active-tree 1.11
  support sweep, public mesh-doc scan, recursive API ambiguity scan, local Markdown
  links, and `git diff --check` all passed.

Re-measured on 2026-08-30 with Julia 1.12.7 after adding cached-mesh refinement
and clearing:

- `API.mesh.refine` now routes the complete cached linear-simplex mesh through
  `Refine.refine_uniform`, commits the new cache only after successful validation
  and resource preflight, and returns independently owned storage. Missing meshes
  fail explicitly, and rejected resource limits preserve the prior cache.
- `API.mesh.clear` idempotently removes the complete cache without changing model
  geometry. Nonempty entity selections are explicit blockers because the simplex
  cache does not own entity-classification metadata.
- The focused lifecycle suites passed 69/69 bounds-checked assertions, including
  fixed first- and second-refinement CRCs
  `6fb8a362968e08263e38a6c59444f7b5503ffdb74b7db54cb90dfd59b683cc69`
  and `09fd5ced56aba7a5b1b0380f8f9189dc95d3676793430f61fd01269baca1445c`.
  Refined curve and surface periodic maps retained 9, 13, and 9 node pairs. The
  neighboring periodic-volume API set passed 17/17 assertions, the direct refinement
  suite retained 118/118, and public-documentation and recursive ambiguity scans
  returned zero.
- The Gmsh 4.15.2 differential retained the exact segment/triangle/tetrahedron
  child-template CRC
  `db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031`,
  verified whole-mesh clearing, and matched the session result to the canonical
  kernel with the first-refinement CRC above.
- The complete bounds-checked package gate passed 168,117/168,117 assertions in
  13m43.3s. Aggregate bounds-checked validation exited 0 in 27m19.45s against Gmsh
  4.15.2, including the extended refinement differential and enclosure/coax probe.
  The organized layout contains 173 tracked Julia files with zero repository-root
  `.jl` files; all 60 `*_test.jl` files are included by the domain-organized test
  entrypoint.

Re-measured on 2026-08-30 with Julia 1.12.7 after adding owned model-entity state:

- `Tessella.Model` and `API.model` now expose entity visibility and RGBA color,
  finite Point-coordinate updates, and global string-attribute lifecycle operations.
  Visibility defaults to `1`, color to `(0,0,255,0)`, recursive setters traverse
  explicit boundaries, attribute values are detached, and names are lexically ordered.
- Visibility and colors follow atomic retagging and are cleaned on removal; attributes
  remain model-global. Presentation and attribute operations preserve the session
  mesh cache. A successful coordinate change invalidates it, while a rejected change
  leaves both geometry and cache untouched. Dependent native geometry immediately
  reflects updated Point coordinates.
- The direct state suites passed 46/46 bounds-checked assertions, the session API
  suite passed 40/40, and the neighboring identity/removal suites passed 294/294.
  Public-documentation and recursive ambiguity scans returned zero. The complete
  bounds-checked package gate passed 168,043/168,043 assertions in 13m56.1s.
- The Gmsh 4.15.2 differential matched 27 visibility queries, 18 color queries,
  three coordinate-update query families, three attributes, and one retagged entity.
  The bounded differences are strict dimensions, finite coordinates,
  documented RGBA ranges, NUL-free strings, unavailable implicit boundaries, and
  immediately updated Plane bounds instead of Gmsh's measured stale result.
- Aggregate bounds-checked validation exited 0 in 27m00.89s against Gmsh 4.15.2,
  including the new model-state differential and enclosure/coax probe. The
  organization ratchet covers 172 Julia source, test, and validation files with zero
  repository-root `.jl` files; all 59 `*_test.jl` files are included by the
  subfolder-organized test entrypoint.

Re-measured on 2026-08-30 with Julia 1.12.7 after completing native surface
reparametrization:

- `Tessella.Model.model_reparametrize_on_surface` and synchronized
  `API.model.reparametrize_on_surface` map Point or straight-Line parameters into an
  explicit Plane's deterministic parameters. Sources need not belong to the Plane
  and off-plane coordinates are orthogonally projected. Native Planes are
  nonperiodic, so the validated `which` selector does not change the result.
- Retagged target Planes retain the operation, malformed or nonfinite inputs and
  unsupported entities fail explicitly, and session success or failure preserves
  the existing mesh cache. The direct evaluation suite passed 63/63 bounds-checked
  assertions and the session API suite passed 44/44; documentation and recursive
  ambiguity scans returned zero.
- The Gmsh 4.15.2 differential passed 6 Point queries, 12 Line query families, 12
  axis-aligned Plane query families, and 5 tilted-Plane query families, including
  boundary-independent and off-plane reparametrization. The complete bounds-checked
  package gate passed 167,957/167,957 assertions in 19m22.6s.
- Aggregate bounds-checked validation exited 0 in 29m36.24s against Gmsh 4.15.2,
  including the surface-reparametrization differential and enclosure/coax probe. The
  organization ratchet remains 168 Julia source, test, and validation files with
  zero repository-root `.jl` files and all 57 `*_test.jl` files included by the
  subfolder-organized test entrypoint.

Re-measured on 2026-08-30 with Julia 1.12.7 after adding native Point/Line/Plane
geometry evaluation:

- `Tessella.Model` and `API.model` now expose value, first- and second-derivative,
  curvature, principal-curvature, normal, parametrization, parameter-bound,
  containment, and closest-point queries for the supported native entity families.
  Retagging preserves evaluation, returned arrays are detached, and session queries
  preserve the mesh cache on success and failure.
- The implementation uses exact-rational straight-Line arithmetic and exact
  coplanarity and loop-membership predicates for physical Plane containment. Plane
  frames follow Gmsh 4.15.2's deterministic axis selection; Line closest points are
  clamped, Plane closest points are untrimmed projections, and invalid or
  unrepresentable inputs fail explicitly.
- The direct evaluation suite passed 53/53 bounds-checked assertions, and the session
  API suite passed 39/39. The complete bounds-checked package gate passed
  167,942/167,942 assertions in 24m35.9s. The neighboring metadata suites remained
  109/109, and public-documentation and recursive ambiguity scans returned zero.
- The Gmsh 4.15.2 differential passed 4 Point queries, 10 Line query families, 12
  axis-aligned Plane query families, and 5 tilted-Plane query families. The bounded
  differences are exact membership instead of a shape tolerance and exact zero
  straight-Line second derivatives instead of Gmsh's numerical noise.
- Aggregate bounds-checked validation exited 0 in 52m14.32s against Gmsh 4.15.2,
  including the new evaluation differential and the enclosure/coax probe. The
  organization ratchet covers 168 Julia source, test, and validation files with zero
  repository-root `.jl` files; all 57 `*_test.jl` files are included by the
  subfolder-organized test entrypoint.

Re-measured on 2026-08-30 with Julia 1.12.7 after adding native entity-property
queries:

- `Tessella.Model.model_entity_properties` and synchronized
  `API.model.get_entity_properties` return detached property vectors for every
  visible native entity. Explicit planes return unit-normal coefficients
  `[a,b,c,d]` for `a*x+b*y+c*z=d`, oriented by their exterior loop; Points, Lines,
  Volumes, and visible primitive or Boolean Volumes return empty vectors.
- Plane equations use exact-predicate coplanarity and exact-rational normal
  construction before the final `Float64` normalization. Queries reject missing
  topology, collinear or
  noncoplanar boundaries, and unrepresentable results; retagging preserves the
  equation and session queries preserve the mesh cache.
- The direct metadata suite passed 109/109 bounds-checked assertions and the session
  API suite passed 36/36. The complete bounds-checked package gate passed
  167,850/167,850 assertions in 46m04.0s. Public-documentation and recursive
  ambiguity scans returned zero.
- The Gmsh 4.15.2 differential matched 20 explicit, retagged, and primitive-volume
  cases through 20 property queries, alongside the existing type, parent,
  partition-membership, and partition-count checks.
- Aggregate bounds-checked validation exited 0 in 102m12.94s against Gmsh
  4.15.2-git, including the property differential and enclosure/coax probe. The
  unusually long wall time includes simultaneous user-owned Julia workloads; no
  such processes were stopped. The organization ratchet covers 164 Julia source,
  test, and validation files with zero repository-root `.jl` files; all 55
  `*_test.jl` files are included by the subfolder-organized test entrypoint.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding native entity metadata
queries:

- `Tessella.Model.model_entity_type` and synchronized `API.model.get_entity_type`
  classify every existing native entity as `Point`, `Line`, `Plane`, or `Volume`;
  `API.model.get_type` is the compatibility synonym. Retagging preserves the type,
  removal makes the old entity unavailable, and primitive solids expose only their
  visible `Volume` entity.
- Native `GeoModel` entities are explicitly nonpartitioned: parent queries return
  `(-1,-1)`, partition-membership queries return detached empty vectors, and the
  model partition count is zero. All metadata queries validate the entity and
  preserve the synchronized mesh cache.
- The direct metadata suite passed 81/81 bounds-checked assertions and the session
  API suite passed 29/29. The complete bounds-checked package gate passed
  167,815/167,815 assertions in 14m27.2s. Public-documentation and recursive
  ambiguity scans returned zero.
- The Gmsh 4.15.2 differential matched 20 explicit, retagged, and primitive-volume
  cases through 40 type/type-alias queries, 20 parent queries, 20 partition-membership
  queries, and two model partition counts.
- Aggregate bounds-checked validation exited 0 in 30m54.96s against Gmsh
  4.15.2-git, including the metadata differential and enclosure/coax probe. The
  organization ratchet covers 164 Julia source, test, and validation files with zero
  repository-root `.jl` files; all 55 `*_test.jl` files are included by the
  subfolder-organized test entrypoint.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding analytical model spatial
queries:

- `Tessella.Model.model_bounding_box` and the synchronized
  `API.model.get_bounding_box` return exact bounds for explicit straight-edge
  topology and native boxes, spheres, oriented cylinders, and cones. Whole-model
  queries union every live entity. Boolean bounds come from the owned operation-time
  result snapshot, so later operand transforms or removal do not change them.
- `model_entities_in_bounding_box` and its session API return detached, sorted
  entities whose complete bounds lie inside a finite query box. Reversed boxes return
  an empty result, embeddings do not enlarge their targets, queries preserve the
  session mesh cache, and corrupt or multiply encoded live geometry is rejected.
- The direct model suite passed 61/61 bounds-checked assertions and the synchronized
  API suite passed 26/26. The complete bounds-checked package gate passed
  167,705/167,705 assertions in 21m51.7s. Public-documentation and recursive
  ambiguity scans returned zero.
- The Gmsh 4.15.2 differential matched 16 exact explicit-model boxes and 506 seeded
  containment queries, spanning 11 nonempty queries and ten result signatures. Four
  analytical primitive boxes remained inside Gmsh's OCC-padded bounds; the largest
  measured padding was `1.0000000116860974e-7`.
- Aggregate bounds-checked validation exited 0 in 28m33.85s against Gmsh
  4.15.2-git, including the spatial-query differential and enclosure/coax probe.
- The organization ratchet covers 160 Julia source, test, and validation files with
  zero repository-root `.jl` files; all 53 `*_test.jl` files are included by the
  subfolder-organized test entrypoint.

Re-measured on 2026-08-28 after narrowing the runtime contract to Julia 1.12.x:

- `Project.toml` now accepts only the Julia 1.12 minor series. A Julia 1.12.7 Pkg
  resolve/instantiate completed, and the parsed compatibility range contains 1.12.7
  while excluding 1.11.9 and 1.13.0. Active user, developer, plan, and validation
  instructions now name 1.12.x as the sole supported runtime; dated older-runtime
  measurements remain labeled as historical provenance.
- The Julia 1.11 allocation allowance was removed from the PostView nearest-point
  ratchet. Its warmed 10,000-query loop allocated zero bytes on Julia 1.12.7, and the
  bounds-checked size-field suite passed 6,976/6,976 assertions in 5m16.9s.
- The complete bounds-checked package gate passed 167,618/167,618 assertions in
  14m21.0s on Julia 1.12.7. The Gmsh 4.15.2 finite-range differential matched all 58
  samples bit-exactly, and the transfinite-hexahedron differential matched eight
  cases, 288 nodes, 96 hexahedra, and 256 boundary quadrangles.
- Recursive public-documentation and ambiguity scans returned zero. The active
  repository sweep found no Julia 1.11 runtime command, compatibility entry, or
  support statement outside the explicitly historical verification record.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding ordered,
dependency-safe model-entity removal:

- `Tessella.Model.remove_entities!` and synchronized `API.model.remove_entities`
  validate a complete request before mutation, preserve request order, block surviving
  boundary or embedding dependencies, and recursively remove explicit boundaries
  while leaving embedded source entities intact. Successful session mutations
  invalidate the detached mesh cache; no-ops and failures preserve it.
- Removal atomically cleans entity names, Physical memberships and emptied groups,
  embedding targets, affected periodic relations, native primitive encodings,
  Boolean-result snapshots, and construction loops made invalid by a deleted
  boundary entity. Entity and Physical allocators remain monotonic. Recursive removal
  stops at primitive and Boolean Volumes because their boundary topology is implicit.
- The Gmsh 4.15.2 differential matched five ordered shared-topology cases, two
  embedding cases, recursive deletion of all 15 explicit tetrahedron entities, and
  864 seeded single/pair/triple request comparisons across recursive and
  non-recursive modes. It also verifies Tessella's bounded choices to reject invalid
  dimensions/nonpositive tags and safely clean native metadata and periodic state.
- The 142 new assertions passed, and the expanded bounds-checked focused gate passed
  495/495 assertions on Julia 1.12.7. The complete bounds-checked package gate passed
  167,618/167,618 assertions in 14m42.0s. Public-documentation and recursive
  ambiguity scans returned zero.
- Aggregate bounds-checked validation exited 0 in 29m44.4s against Gmsh
  4.15.2-git, including the entity-removal differential and the enclosure/coax probe.
  The organization ratchet covers 156 versioned Julia source, test, and validation
  files with zero repository-root `.jl` files; all 51 `*_test.jl` files are included
  by the subfolder-organized test entrypoint.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding owned model-entity
identity:

- `Tessella.Model` and the synchronized `API.model` now get, set, clear, and remove
  dimension-scoped entity names, and atomically move an existing entity to an unused
  positive tag. A successful mutation invalidates the session mesh cache; getters,
  no-ops, and failed atomic changes preserve it.
- Retagging migrates Point/Curve/Surface/Volume topology, Point sizes, Physical
  memberships, embedding sources and targets, periodic slave/master relations,
  primitive encodings, Boolean-result snapshots, and the entity name. Automatic tags
  remain monotonic. Boolean operand tags retain their operation-time provenance.
- The direct Gmsh 4.15.2 differential matched entity sets, boundaries, Physical
  memberships, embeddings, periodic pairs `(101,102)` and `(102,103)`, retagged box
  bounds, and nonempty volume meshes. The Tessella box mesh was unchanged at CRC
  `e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83`.
  The same probe records the bounded differences: Tessella requires dimensions 0–3
  and positive `Int32` tags, stores names only for existing entities, and moves names
  on retag; Gmsh accepts nonpositive targets, ignores an invalid dimension, stores a
  missing-entity name, and leaves names on old tags.
- The 152 new assertions passed, and an expanded bounds-checked focused gate passed
  1,362/1,362 assertions under both Julia 1.12.7 and Julia 1.11.9. Public-documentation
  and recursive ambiguity scans for `Tessella.Model`, `Tessella.API`, and
  `Tessella.GeoExec`, plus direct checks for the four new `API.model` bindings,
  returned zero under both versions.
- The bounds-checked package gate passed 167,476/167,476 assertions in 13m35.0s.
  Aggregate bounds-checked validation exited 0 in 28m48.1s against Gmsh
  4.15.2-git, including the new model-identity differential and the enclosure/coax
  probe.
- The organization ratchet covers 152 versioned source, test, and validation Julia
  files with zero repository-root `.jl` files. All 49 `*_test.jl` files are included
  by the subfolder-organized test entrypoint.

Re-measured on 2026-08-28 with Julia 1.12.7 after making Boolean volumes own
operation-time operand geometry:

- Native Boolean creation now snapshots both operand surfaces before model mutation.
  Retained-operand transforms, deleting either operand, nested Boolean deletion, and
  reusing deleted tags for a different primitive leave earlier Boolean geometry
  unchanged. Provenance remains `(op, a, b)` and no longer acts as a live dependency.
- API and `.geo` Boolean deletion share one model-owned cleanup path. It removes the
  visible volume, every native solid encoding for its tag, target volume embeddings,
  and dimension-3 Physical memberships. Nonempty groups retain their name; empty
  groups and names are removed; automatic entity and Physical counters stay
  monotonic. Failed Boolean creation preserves model state and the API mesh cache.
- The exact Gmsh 4.15.2 differential matched retained-operand result stability,
  deleted-tag reuse, surviving Physical group `(3,10)` with member `[4]` and name
  `mixed`, and final Volume entities `[1,2,3,4]`. Both kernels retained result volume
  `1`; Gmsh generated 100 tetrahedra and Tessella 12.
- The base box-difference regression checksum is `(nodes=9, segments=0, triangles=0,
  tetrahedra=12, bbox=((1.0,0.0,0.0),(2.0,1.0,1.0)), dihedral=(0.7853981633974484,
  0.7853981633974484), radius_edge=(0.8660254037844387,0.8660254037844387),
  boundary_faces=12, sha=8c402b5d221eb617f78f508d089f61c8471e20f3b9f24fcc3ed6bac90ea0ea30)`.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9:
  713/713 assertions per version, including 82/82 new model, `.geo`, nested-Boolean,
  deletion, tag-reuse, Physical/embedding, atomicity, API-cache, and diagnostic
  assertions. Documentation and recursive ambiguity scans for `Tessella.Model`,
  `Tessella.GeoExec`, and `Tessella.API` each returned zero.
- The bounds-checked package gate passed 167,324/167,324 assertions in 13m35.9s.
  Aggregate bounds-checked validation exited 0 against Gmsh 4.15.2-git, including
  the expanded Boolean snapshot/Delete differential and the enclosure/coax probe.
- The organization ratchet covers 148 source, test, and validation Julia files,
  with zero repository-root `.jl` files. All 47 `*_test.jl` files are included by
  the subfolder-organized test entrypoint.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding explicit model-topology
queries:

- `Tessella.Model` and the synchronized `API.model` now enumerate explicit
  entities, report model dimension, return direct/combined/oriented/recursive
  boundaries, and return upward/downward adjacencies. Results are detached and
  deterministic; queries do not invalidate the session mesh cache. Hole and cavity
  orientation participates, while embeddings do not become topology. Primitive and
  Boolean volumes are enumerated but block boundary and downward-adjacency queries
  because their subentities remain implicit.
- The exact Gmsh 4.15.2 differential matched 15 nonmonotonic-tag entities, 20 boundary
  cases, and 10 adjacency cases. It covers signed direct incidence, combined parity
  cancellation, recursive Point closure, mixed dimensions, holes, cavities, and
  embedded entities. The query fixture retained 5 nodes, 4 tetrahedra, and mesh CRC
  `71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b`.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9: model
  topology 80/80, session topology/cache ownership 39/39, model lifecycle/entity
  behavior 119/119, native periodic relations 118/118, model/`.geo` execution 84/84,
  topology/Point sizing 240/240, API session/option behavior 64/64, API Physical
  lifecycle 62/62, API Point sizing 35/35, explicit volumes 26/26, periodic volume
  boundaries 12/12, periodic ownership 35/35, embedded periodic curves 27/27,
  periodic dependency graphs 49/49, periodic projection 93/93, 51/51, and 9/9, and
  volume projection 72/72, 13/13, 80/80, and 60/60. Public-documentation and
  recursive ambiguity scans for `Tessella.Model` and `Tessella.API`, plus direct
  documentation checks for all 13 query/lifecycle bindings in `API.model`, returned
  zero under both versions.
- The bounds-checked package gate passed 167,242/167,242 assertions in 13m50.6s.
  The aggregate bounds-checked validation exited 0 in 26m31.0s against Gmsh
  4.15.2-git, including the model-topology differential and unchanged CRC above.
- The organization ratchet covers 146 Julia files with no repository-root `.jl`
  files. Source, test, and validation code remain divided among domain subfolders;
  the only top-level Julia entrypoints are `src/Tessella.jl`, `test/runtests.jl`, and
  `validation/run_all.jl`.

Re-measured on 2026-08-28 with Julia 1.12.7 after completing the native
Physical-group lifecycle:

- The model and synchronized session API now query sorted, detached Physical groups,
  memberships, reverse memberships, names, and name-selected entities. Names are
  unique within a dimension and reusable across dimensions. Gmsh 4.15.2-compatible
  name assignment, cross-dimensional name removal, selective/all-group removal, and
  monotonic automatic tags are covered. Group removal leaves geometry intact; failed
  or no-op mutations are atomic and preserve the session mesh cache.
- The exact Tessella/Gmsh API differential matched automatic tags `1,2,3`,
  cross-dimensional names, same-dimension duplicate-name omission, all query results,
  name mutation/removal, selective removal, three final groups, and replacement tag
  `4`. It retained dynamic mesh CRC
  `2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808`,
  projected CRC
  `99aeefc2e269090b518f5896f2388bd419c8fb44d06e4743579d50d61aaedf81`,
  and `SetMaxTag` CRC
  `aa3127e5c1a302ebdb98890a4d7a62d5cb385e3cf73aa024372f809a7542a45d`.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9: model
  lifecycle/entity behavior 119/119, native periodic relations 118/118, model/`.geo`
  execution 84/84, `.geo` dynamic tags 67/67, topology/Point sizing 240/240,
  API session/option behavior 64/64, API lifecycle 62/62, API Point sizing 35/35,
  explicit volumes 26/26, periodic volume boundaries 12/12, periodic ownership
  35/35, embedded periodic curves 27/27, periodic dependency graphs 49/49, IO
  426/426, CLI 131/131, periodic projection 93/93, 51/51, and 9/9, and volume
  projection 72/72, 13/13, 80/80, and 60/60.
  Public-documentation and recursive ambiguity scans for `Tessella.Model` and
  `Tessella.API`, plus direct documentation checks for the nine `API.model` lifecycle
  functions, returned zero under both versions.
- The bounds-checked package gate passed 167,123/167,123 assertions in 14m02.6s.
  The aggregate bounds-checked validation exited 0 in 34m54.0s against Gmsh
  4.15.2-git, including the lifecycle differential and unchanged CRCs above.
- The organization ratchet covers 143 tracked `.jl` files with no repository-root
  `.jl` files. Source, test, and validation code remain divided among their existing
  domain subfolders; the only top-level Julia entrypoints are `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding global automatic
Physical-group tags:

- Name-only `.geo` declarations such as `Physical Surface("skin") = {...};`
  allocate from one Physical counter across Point, Curve, Surface, and Volume.
  Explicit named and unnamed positive tags remain supported. The native model and
  session API use the same global counter, while explicit tags remain keyed by
  dimension. Empty automatic names, duplicate names within one dimension, duplicate
  dimension/tag pairs, malformed headers, and signed-32-bit exhaustion are atomic
  blockers.
- The Gmsh 4.15.2 parser and Tessella both assigned tags 1, 2, 20, 21, and 22 to
  the mixed-dimension parser corpus and made subsequent `newreg` lines 4 and 23.
  Tessella and Gmsh API calls both returned global automatic tags `1,2,3`. The
  differential reported five parser groups and retained dynamic mesh CRC
  `2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808`,
  projected CRC
  `99aeefc2e269090b518f5896f2388bd419c8fb44d06e4743579d50d61aaedf81`,
  and `SetMaxTag` CRC
  `aa3127e5c1a302ebdb98890a4d7a62d5cb385e3cf73aa024372f809a7542a45d`.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9:
  geometry/Physical allocators 67/67, IO 426/426, Point sizing 240/240, CLI
  131/131, and the API automatic-tag set 8/8. Neighboring `SetMaxTag`,
  geometry-expression, numeric-list, and periodic-executor sets passed 49/49,
  49/49, 44/44, 18/18, and 132/132 under both versions. Public-documentation and
  recursive ambiguity scans for `Tessella.Model`, `Tessella.GeoExec`, `Tessella.IO`,
  and `Tessella.API` returned zero under both versions.
- The bounds-checked package gate passed 167,023/167,023 assertions in 16m30.9s.
  The aggregate bounds-checked validation exited successfully in 26m45.3s against
  Gmsh 4.15.2-git, including the automatic Physical parser/API differential and the
  unchanged topology-query projected CRC
  `608dcd81b4ecab3138fd8da610ec109e1972c799ccb2c9d755917c9901250905`.
- The organization ratchet still covers 143 tracked `.jl` files with no
  repository-root `.jl` files. The only top-level Julia entrypoints remain
  `src/Tessella.jl`, `test/runtests.jl`, and `validation/run_all.jl`.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding topology-derived
Physical groups:

- The bounded `.geo` executor accepts inline `Boundary` and `CombinedBoundary`
  queries for Physical Point/Curve/Surface groups and `PointsOf` for Physical Point.
  `Boundary` joins immediate explicit boundaries before group membership is
  deduplicated; `CombinedBoundary` retains odd-multiplicity tags. Signed selectors,
  multiple entity blocks, `:`, expressions, ranges, and numeric-list variables are
  supported. Each query accepts at most 65,536 input entities and traversed boundary
  occurrences. Hole and cavity boundaries participate; embeddings do not. Invalid
  dimensions, unknown or implicit primitive topology, and empty combined results are
  atomic, explicit blockers.
- Named and unnamed explicit-tag Physical declarations execute. The metadata scanner
  applies the same bounded checks to both forms and records only nonempty names. The
  explicit tetrahedron matched Gmsh 4.15.2 memberships `[1,3]`, `[1,2,3]`,
  `[1,2,3,4]`, `[2,3,4,5]`,
  `[1,2,3,4,5]`, and recursive Volume Points `[1,2,3,4]`. Its native mesh remained
  6 nodes and 6 tetrahedra with CRC
  `cf091ac13ba325f5f68650b192598a5159679b40e19dbea744d66c58962357b5`;
  the expanded classified projection CRC is
  `608dcd81b4ecab3138fd8da610ec109e1972c799ccb2c9d755917c9901250905`.
  Gmsh exposed all six hidden OCC Box boundary surfaces; Tessella returned the
  documented explicit-topology blocker.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9: topology
  queries and Point sizing 240/240, IO 425/425, CLI 131/131, and volume projection
  72/72, 13/13, 80/80, and 60/60. Neighboring `SetMaxTag`, geometry-expression,
  numeric-list, periodic-executor, and periodic-projection sets passed 49/49, 49/49,
  44/44, 18/18, 132/132, 93/93, 51/51, and 9/9 under both versions.
  Public-documentation and recursive ambiguity scans for `Tessella.Model`,
  `Tessella.GeoExec`, and `Tessella.IO` returned zero under both versions.
- The bounds-checked package gate passed 166,986/166,986 assertions in 13m08.9s.
  The aggregate bounds-checked validation exited successfully in 24m13.9s against
  Gmsh 4.15.2-git, including the updated topology-query differential.
- The organization ratchet still covers 143 managed `.jl` files with no
  repository-root `.jl` files. Topology queries remain in
  `src/geometry/ModelTopologyQueries.jl`; the only top-level Julia entrypoints are
  `src/Tessella.jl`, `test/runtests.jl`, and `validation/run_all.jl`.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding explicit-topology
`PointsOf` mesh-size selectors:

- The bounded `.geo` executor accepts inline `PointsOf` blocks for Point,
  Curve/Line, Surface, and explicit Volume entities in `MeshSize` and
  `Characteristic Length` statements. Multiple blocks, signed tags, `:`, bounded
  expressions/ranges, and numeric-list variables resolve to sorted, unique recursive
  boundary Points. Hole boundaries participate; embeddings do not. Implicit primitive
  and Boolean volume topology remains a precise blocker.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9: Point
  mesh-size 197/197 and CLI 131/131. The explicit-volume mesh has 6 nodes and 6
  tetrahedra with native CRC
  `cf091ac13ba325f5f68650b192598a5159679b40e19dbea744d66c58962357b5`;
  its CLI MSH projection CRC is
  `a03e62cdb5b5049f0c3ac79f829707cab1ce7575ade8a2b246808a7e222e50ca`.
- The Gmsh 4.15.2 differential matched stored Point sizes
  `[0.4, 0.5, 0.2, 0.4, 1.0]` and recursive Volume boundary Points
  `[1, 2, 3, 4]`. Its OCC Box exposed eight hidden Points with the requested size;
  Tessella returned the documented explicit-topology blocker for that volume.
- The bounds-checked package gate passed 166,938/166,938 assertions in 15m05.9s.
  The bounds-checked aggregate validation gate exited successfully against Gmsh
  4.15.2-git, including the updated direct/topology Point-size child.
  Public-documentation and recursive ambiguity scans for `Tessella.Model` and
  `Tessella.GeoExec` returned zero under both Julia versions.
- The organization ratchet covers 143 managed `.jl` files with no repository-root
  `.jl` files. `ModelTopologyQueries.jl` is in `src/geometry`; the only top-level
  Julia entrypoints remain `src/Tessella.jl`, `test/runtests.jl`, and
  `validation/run_all.jl`.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding spatial planar Point-size
propagation:

- `mesh_model_surface` now carries participating positive Point constraints into
  surface refinement.
  Generated straight-curve subdivision nodes receive linear endpoint values, and
  nonuniform constraints are interpolated over the deterministic initial constrained
  triangulation. Uniform constraints retain the previous constant-size path exactly;
  coincident PSLG inputs use the smaller constraint.
- The bounds-checked Point mesh-size set passed 166/166 assertions under Julia
  1.12.7 and Julia 1.11.9. Public-documentation and recursive ambiguity scans for
  `Tessella.Model` returned zero under both versions. Sixty seeded nonuniform square
  cases and a holed surface with embedded Point and Line constraints all produced
  valid meshes; the measured longest-edge/centroid-size ratio never exceeded 1.
- The Gmsh 4.15.2 differential retained the four stored constraints and added a
  spatial square with sizes `[0.1, 0.8, 0.8, 0.8]`. Tessella produced 72 triangles
  with quadrant counts `[45, 8, 11, 8]`; Gmsh produced 68 with
  `[41, 11, 8, 8]`. Both localized the fine constraint, and both areas were 4 within
  `128eps(Float64)`. The updated native CRC is
  `b3f1bf410e917d050eacceab998b0fdf7b4cd61d1d9f263805b5120c06f1f4df`;
  the classified projection CRC is
  `b7202dfa1cfb7469e7541c34e2b1bfae404c66f2462abc1953fa0b9374e5a010`.
  Exact Gmsh topology is not claimed.
- The bounds-checked package gate passed 166,902/166,902 assertions in 17m48.8s.
  The bounds-checked aggregate validation gate exited successfully against Gmsh
  4.15.2-git, including the updated spatial Point-size child.
- The organization ratchet covers 142 managed `.jl` files with no repository-root
  `.jl` files. `SurfacePointSizing.jl` is in `src/geometry`; the only top-level Julia
  entrypoints remain `src/Tessella.jl`, `test/runtests.jl`, and
  `validation/run_all.jl`.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding bounded Point mesh-size
constraints:

- `set_point_mesh_size!` atomically validates existing Point tags and one finite,
  positive size before mutation. `API.mesh.set_size` supports dimension-tag pairs
  for Points and invalidates its detached cache only after success. The bounded
  `.geo` executor supports `MeshSize` and `Characteristic Length` with `:`, numeric
  expressions/ranges, and whole or selected numeric-list variables.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9: Point
  mesh-size 67/67, API 248/248, and CLI 126/126. Public-documentation and recursive
  ambiguity scans returned zero under both versions.
- The required Gmsh 4.15.2 differential matched the four stored constraints
  `[0.8, 0.4, 0.6000000000000001, 0.4]` and the API update `[0.35, 0.35]`.
  Tessella produced 23 nodes and 28 triangles; Gmsh produced 17 nodes and 22
  triangles. Their maximum unit-square area error was
  `1.1102230246251565e-16`. The native CRC is
  `f1bfa8a1cc61158cc6293540ad6ce6d7ce48054a616a3df57d93597857f1b089`;
  the classified projection CRC is
  `1489fd244841d079350f33439a97b6f33205bcff32e4b99ac672e79a4387eaca`.
  The differential also pins the narrower native contract: Gmsh accepts zero or a
  negative size and silently ignores a pre-creation or missing-Point assignment, while
  Tessella requires positive sizes on Points that already exist.
- The bounds-checked package gate passed 166,803/166,803 assertions in 17m24.7s.
  The bounds-checked aggregate validation gate exited successfully against Gmsh
  4.15.2-git, including the new Point mesh-size child.
- The organization ratchet covers 141 managed `.jl` files with no repository-root
  `.jl` files. The only top-level Julia entrypoints remain `src/Tessella.jl`,
  `test/runtests.jl`, and `validation/run_all.jl`; every other managed Julia file is
  in a domain or workflow subfolder.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding bounded,
factory-aware `.geo` `SetMaxTag` execution:

- `SetMaxTag Point|Curve|Surface|Volume` now uses independent Built-in and
  OpenCASCADE counters. Built-in can lower its counter, OpenCASCADE only raises
  its counter, and allocator reads use the greatest counter among activated
  factories. Fractional conversion, allocator-backed expressions, negative
  counters, factory activation, signed-32-bit exhaustion, explicit topology,
  Physical groups, Fields, and full primitive allocation are checked. Primitive
  boundary tags remain implicit, preserving existing modeled-subentity tag reuse.
- Bounds-checked focused sets passed under Julia 1.12.7 and Julia 1.11.9:
  `SetMaxTag` 49/49, dynamic tags 50/50, model sets 62/62, 118/118, and 84/84,
  IO 422/422, and CLI 122/122. Public-documentation and recursive ambiguity scans
  returned zero under both versions.
- The required Gmsh 4.15.2 differential covered all four `SetMaxTag` categories,
  Built-in lowering, OpenCASCADE raise-only behavior, cross-factory maxima,
  fractional and allocator-expression values, negative factory activation, and
  occupied primitive allocation. The new tetrahedron fixture produced 5 nodes and
  4 tetrahedra in Tessella and 17 nodes and 23 tetrahedra in Gmsh. Its native CRC is
  `71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b`; its
  classified projection CRC is
  `aa3127e5c1a302ebdb98890a4d7a62d5cb385e3cf73aa024372f809a7542a45d`.
  The combined dynamic-tag maximum volume error was
  `8.881784197001252e-16`.
- The bounds-checked package gate passed 166,696/166,696 assertions in 22m36.3s.
  The bounds-checked aggregate validation gate passed in 2,171.35s against Gmsh
  4.15.2-git, including the new `SetMaxTag` differential.
- The organization ratchet passed with 139 managed `.jl` files and no
  repository-root `.jl` files. The only top-level Julia entrypoints remain
  `src/Tessella.jl`, `test/runtests.jl`, and `validation/run_all.jl`; every other
  Julia file is in a domain or workflow subfolder.

Re-measured on 2026-08-28 with Julia 1.12.7 after adding bounded dynamic `.geo`
tag allocators:

- The scanner and executor now evaluate read-only `newp`, `newf`, and Gmsh's
  shared curve/loop/surface/volume allocator aliases while every preceding
  tag-producing statement remains in the tracked subset. Physical-group tags and
  the hidden Point/Curve/Surface topology of full Box, Cylinder, Sphere, Cone,
  and Torus primitives advance the same namespaces Gmsh 4.15.2 advances. Allocator reads
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
  numeric functions. Entity lists also expand bounded ranges — literal or
  variable-backed — and
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
