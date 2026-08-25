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
| P2 | **IN PROGRESS** | 125 fixed-node Gmsh types plus ten serializable cut/border/child/sub-element records, mixed blocks/entities/classification metadata, structural validation/CRC, and ASCII/binary MSH v2.2/v4.1 read/write |
| P3 | **IN PROGRESS** | Native analytical surfaces/imprints, classified ISO-10303-21 STEP/IGES box/sphere/cylinder/cone import, STEP/IGES NURBS curve and surface import with IGES NURBS export, Box/Cylinder/Sphere/Cone/Boolean/Translate/Dilate/90°-Rotate `.geo` execution, mesh Boolean CSG, and finalized-mesh affine transforms |
| P4 | **IN PROGRESS** | Greedy and Edmonds-blossom surface recombination with optional full-quad, Point/Line-In-Surface embeddings, Point/Line/Surface-In-Volume recovery, holed plane surfaces, uniform refinement, Progression/Bump/Beta curve laws and HWall variants, planar triangle/quad transfinite patches, affine five-/six-face transfinite volumes, recombined hexahedra, prismatic 3-D layers with certified remaining-core fill/cavity walls, 2-D quad/fan layers, and translation-periodic node-pair certification/snapping |
| P5–P6 | **IN PROGRESS** | Synchronized model/mesh API with detached cache ownership, non-destructive bounded CLI, validated headless GUI, owned scalar nodal views, synchronized in-process plugins, plus t1-square, t4-hole, Point/Line-In-Surface, Surface-In-Volume sheet, translation-periodic curve, 2-D boundary-layer quad, API-box, OCC-cylinder/cone, IGES-128 bilinear, and BooleanDifference box Gmsh 4.15.2 differentials |

P1 does not claim 3-D multi-wall boundary-layer fans, the full Gmsh automatic-sizing
pipeline, high-order/custom-interpolation, or mixed-component
`PostView`, materially warped quadrangles, `PostView` tensor-to-metric evaluation, direct tensor or
metric-meshing parity, full `.geo`/CAD-model execution, or exact CAD distance. P2 does
not claim general mixed-element generation or recombination beyond P4's first-order
surface pairing, MINI basis-selector tags 138/139 as mesh records, integration of
mixed blocks into the simplex meshing kernels, curved high-order
Jacobian certification beyond P2 tetrahedra, ancillary/unknown-section preservation (binary readers reject
unsupported sections explicitly), repeated `$Nodes`, non-8-byte binary data, internal
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
`.geo` execution covers Point/Line/Loop/Surface, Box/Cylinder/Sphere/Cone,
BooleanDifference/Union/Intersection of those solids, Translate of remaining
native solids, Dilate, coordinate-axis π/2 rotations of AABB boxes,
Point/Line-In-Surface embeddings, and Point/Line/Surface-In-Volume recovery. Its
scanner handles finite arithmetic constants,
pure numeric functions, prior scalar bindings, explicit field/physical tags, and
finite constant ranges in recognized numeric field lists/selectors. Entirely numeric
Physical memberships are range-checked but remain geometry data. It rejects loops,
macros, dynamic tags, option reads, stateful functions, dynamic/general ranges,
logical/ternary evaluation, extrusions/fillets/symmetry, and mixed geometry-derived
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
Planar constant-`z` polylines extrude to type-3 quadrangles along left-normals,
with optional convex-corner fans. Closed manifold walls can use
`mesh_boundary_layer_filled` for certified prism shells, cavity walls, and a
conforming remaining-core tetrahedral fill. Explicit one-to-one translated node
pairs can be certified and snapped exactly without changing connectivity or tags;
the pair map remains caller-owned rather than persistent model metadata. P4 does
not yet claim
non-affine CAD curve integration, FlexibleTransfinite, or size-map curve laws,
quasi-transfinite or holed patches,
general CAD parameterizations, recombined three-sided patches, curved/warped or
compact-TransfiniteTri volumes, volume/hybrid
recombination, selective or high-order refinement, coarsening,
3-D multi-wall boundary-layer fans, or persistent model/I/O periodic entity metadata.

General OpenCASCADE/unclassified NURBS CAD, remaining algorithms/fields, broad
formats and API, GUI, and post-processing are unfinished parity tracks, not
project non-goals.

## Current worktree verification

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
