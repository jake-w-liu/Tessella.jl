# Implementation audit handoff to Devin — 2026-10-11

The user requested a complete unfinished-work handoff. New implementation, Julia runs and SDK probes have stopped. This file records the implementation_audit agent's concrete continuation scope. The final verified batch and this handoff are published together; the next owner resumes implementation from current Main. No owned runtime remains active.

This design was prepared against parent Main `da92208`. The successor224-check
batch integrates shared periodic, STEP SI units and view representation replacement.
Its final ordinary tests pass on both runtimes, actual exit0,45.21/38.38 s.
Resume from current `origin/main` and [the Devin checklist](../../DEVIN_HANDOFF.md);
older candidate paths below preserve provenance. All P1–P6 remain open.

## Last completed packet: shared translated periodic hybrid Regions

Frozen manifest: `C:/tmp/tessella_hybrid_shared_periodic_20261011/ready.json`. Patch: `shared_periodic_source.patch` in that directory. Root has integrated these three selected paths:

- `src/geometry/ModelNativeHybridBatch.jl`: saved base SHA256 `62b3b75401e4eaef977611ccb0d1fb8838bcb765543dda591dbc62489e67c40e`; final `a9f1518a2c29c915f62e206c8ad4fcdffd0216845d3834c7f9e8c3e5b837a53c`.
- `test/interfaces/api_hybrid_shared_periodic_functional_test.jl`: SHA256 `c8d5a4bef0240d40cf847b7878a41967efbc5800f32fb09cb048057c50feb60f`.
- `docs/shared_periodic_hybrids.md`: historical candidate SHA256 `1aca94734bc2f5ee7ea62f19b584e12d89c780bdc0f5c682482ef0637530e7c2`; the integrated documentation adds final qualification details.

Actual final qualification: Julia 1.12.7, `--startup-file=no --compile=min -O0 --check-bounds=yes --threads=1`, **63/63, actual exit 0, 36.9840977 s**. Receipt and stdout/stderr: `receipts/attempt2.*` in the candidate. PID 25356/session 94899 drained. This is functional evidence, not an optimized JIT performance claim.

Two adjacent unit Regions share two actual Tri3 cells and each retain an exterior translated copied Surface pair. Both emit 21 Tet4 + 5 Pyramid5, positive material volume 1 each; shared Tri facets have one carrier per Region and opposite sides. Public maps retain four primary and nine P2 pairs. Owned MSH4 roundtrip preserves actual UID, owner and ordered connectivity. Explicit `Mesh.Renumber=0` preserves source node labels; the ordinary SetOrder1 factory still renews linear cell labels by design.

The new admission checks uniquely owned exterior periodic Faces, existing complete per-Face constructor/copy-event witnesses, and consistent actual masterUID/affine bindings for repeated shared slaveUIDs before RNG. It synthesizes no copy pairs. Root and concise_tests independently reviewed that Source and its called complete certificate.

Scope remains Pyramid mode 0, translations, one pair per Region. The final two-check invalid fixture exercises an earlier generate2 lower periodic-Curve guard; it does **not** qualify the new joint-copy inconsistency gate. Attempt1 (actual exit1/70.1031068 s, positives59 + Pyramid guard2 green) is preserved. No private contradictory map was fabricated to make a test pass.

## Implemented design: proper rigid periodic hybrid pair

Status: **implemented and verified (2026-10-11)** per the design below. `_surface_periodic_constraints` now routes actual Surface-derived cross-side Curve relations (`ModelPeriodicConstraint.derived`) to the owned lower producer instead of rejecting them; directly declared cross-surface relations keep the explicit guard. The hybrid admission gate certifies the applied 3×4 map as proper rigid — finite coefficients, orthonormal linear part, determinant +1 within a scale-aware `1e-12` bound — while the unchanged remainder still requires the actual copy UID events, complete constructor rows, Point/Curve/Surface ownership, affine coordinates and ordered cell incidence. Focused coverage: `test/interfaces/api_rigid_periodic_hybrid_functional_test.jl` (61 assertions: quarter- and half-turn 24Tet+6Pyr cubes, P1/P2 maps, MSH4 round-trip, atomic shear/reflection and direct-relation rejections) plus the `rigid_periodic_smoke.jl` default check. Cross-Region and mixed-direction copies moved their atomic rejection to the volume certificate layer. Original probe drivers and observations remain under `C:/tmp/tessella_hybrid_next_scope_20261011`; portable JSON observations are linked in the repository evidence section below:

- `native_rigid_periodic_probe.py` / `.json`: pinned public Gmsh 4.15.2 quarter-turn and half-turn unit-cube controls.
- `native_shared_periodic_probe.py` / `.json`: valid shared translated two-Region control for the integrated feature.
- `native_shared_periodic_probe.v1.py` / `.json`: earlier transfinite slave override produced no copy maps; unqualified as a periodic-copy witness.
- `native_quad_sheet_probe.py` / `.json`: retained Quad4 sheet control described below; Native returned but did not produce a valid constrained volume.

The rigid controls copy top Surface2 from bottom Surface1, leaving the slave to the copy producer. Rz90 affine is `[0,-1,0,1, 1,0,0,0, 0,0,1,1, 0,0,0,1]`; Rz180 is `[-1,0,0,1, 0,-1,0,1, 0,0,1,1, 0,0,0,1]`. Both actual SDK outputs have 24 Tet4 + 6 Pyramid5, every typed cell positive, volume1, exact six-Quad exterior, facet multiplicities1/2, coordinate bounds, proper orthogonal determinant-positive affine, four primary and nine P2 actual pairs. Both JSON `qualified_control` fields are true. SDK bodies took 0.01346159994 and 0.00243679993 seconds; final Python wrapper exited0. No Native implementation Source was read.

There are two distinct production blockers and a bounded implementation plan:

1. `src/geometry/Model.jl`, `_surface_periodic_constraints` (around3352), rejects a dim1 relation when only one Curve belongs the current Surface. A Surface-derived rigid relation can couple Curves on different side Surfaces. **Existing provenance is `ModelPeriodicConstraint.derived::Bool` (around116)**; `set_periodic!` marks actual Surface-induced boundary correspondences true. Do not invent a flag or infer provenance from XYZ. Inspect the actual generation/lower-authority call around8499 before periodic Surface copy dispatch around8518. A repair should let the owned global lower producer enforce a certified Surface-derived cross-side relation while preserving direct independent Curve constraints and unrelated public paths. The proposed rigid route reaching this guard is a Source-based inference; rotation has not been run in Tessella. The mixed-direction negative did actually reach it.
2. `src/geometry/ModelNativeHybridVolume.jl`, `_model_native_hybrid_periodic_boundary_certificate` (around90–145), explicitly accepts translations only at line98. The remaining certificate already checks complete actual constructor rows, bijective actual UID events, Point/Curve/Surface node ownership, affine coordinates and ordered actual cell incidence. Replace only that admission with a finite **proper rigid** certificate: orthogonality, determinant+1 and homogeneous final row, with a scale-aware rounding bound separate from geometric `atol`. Do not admit shear/reflection by widening to all affine maps, alter stored coordinates/order, or infer pairs.
3. Reuse `ModelGmshSurfaceProducer.jl::_model_periodic_surface_production`, `ModelPeriodicUIDHistory.jl`, Higher retained Surface history and P2 lifting unchanged unless an actual failure proves a transport bug. The existing Surface producer already transforms owned master cells/rows using the stored affine and actual lower UID binding.

First meaningful target: one quarter-turn cube, actual lower retention, 24Tet+6Pyr positivity/volume/exterior, actual4→9 P2 maps, MSH4 transport, and precise atomic admission errors for unsupported affine/direct foreign Curve context. Review the derived-Curve scope and its callers before edits. No default Volume workload or generic rotation/multiple-direction parity claim. Multiple periodic directions and cross-Region periodic Face ownership remain a later design; no valid native control was established for them here.

## Remaining hybrid production gaps

Implemented: ordinary constructor transport, separate and connected multi-Region production with literal shared Tri interfaces; negative normal Pyramid mode, nonnegative relocation modes, Winslow1; separated/shared-Tri nonzero batches; retained Point/Curve collections; one endpoint-only planar disk with3–6 vertices and `n-2` actual Tri3; default-mode translated periodic pair per Region including shared-Tri components.

- **Broader retained sheet topology:** `ModelNativeHybridEmbeddedSurface.jl` still requires one Region and one embedded planar Surface, one boundary loop,3–6 straight Lines, one actual Line2 per boundary Curve, endpoint-only actual Point rows, and `n-2` Tri3. Larger disks, extra boundary subdivisions, interior Surface nodes, holes, multiple sheets, curved sheets and Quad sheets are genuine additional producer scopes. Next substantial increment should use an actual disk triangulation with interior Surface vertices, extending topology/Euler/dual and exact nonoverlap proofs without reconstructing caller triangulation. Preserve sizing witness in `ModelNativeHybridSizing.jl`, source UID/owner binding, protected facet/edge unions, full node spherical stars, cyclic edge links and exactly-two opposite Tet facet incidence through recovery/carve/insertion/optimizer.
- **Quad sheet is not a demonstrated successful Native context:** the opaque Quad4 probe generated only6 boundary Pyramids, zero Tets, no Quad typed carrier and lost its original UID/order. A Native return alone does not justify admission. Keep the Quad guard until a genuine positive constrained volume or an independently specified owned producer exists.
- **Optimization/algorithm contexts:** `ModelNativeHybridVolume.jl` and Batch explicitly reject Netgen optimization; the independent equivalent kernel is missing. Retained constructors are actual linear Tri/Quad/Line cells. Curved/nonlinear retained source contexts need genuine producer and placement certificates, rather than simply widening catalog admission. `ModelNativeHybridComponents.jl` retains selected planner/connected scheduling guards. Native Region reverse/smoothing were observed no-ops and are not missing native Region features. Shared-Quad nonzero two-Region composition is unsupported by the pinned native control; retain its explicit diagnostic.

The last polygon packet is integrated: `C:/tmp/tessella_hybrid_polygon_sheet_20261011/ready.json`, final48/48 Julia112 exit0/26.823576 s (`receipts/attempt5.*`). Selected Surface helper SHA256 `49a4648c3190c86e60f66aa3720819c828b697775a4a4d02bf5409b6d2e9f2f5`; `Mesh3DGmshOptimize.jl` `62e97a1073eef0103bc985760935dcb7a1367be2920aa04c1f7c259984d1fceb`. The optimizer fix canonicalizes only validated internal facet/edge lookup keys. It fixes a real post-optimizer lost embedded facet while preserving original owned cells/UID/XYZ/order and final certificates. Original failed receipts and pre/post-carve/optimizer evidence remain; do not restore the old keys.

## Remaining ordinary layer/extrusion production gaps

Ordinary Layers and numerous QuadTriNoNew single-source catalogs already exist; do not restart the stale ordinary transport problem as though no producer had been integrated. Source domains are explicit, and operation-local plans prevent stale cached geometry/layer products.

- `APIGenerateHigher.jl` around765/775 requires actual linear Curve/Surface cells for ordinary Layers; native order elevation after generation exists, but direct nonlinear retained lower-source sweeps remain additional implementation.
- `ModelExtrudeNoNewScope.jl::_extrude_nonew_scope_contract` guards transfinite cap/lateral overrides (constrained-boundary planner), reversed Surface/Curve/Volume (oriented-boundary planner), periodic/compound boundary entities (neighboring-region planner). The actual scoped builder around463 and Surface lookup around497 reject a Surface owned by multiple selected sweeps: a **shared-region planner** is still missing for this NoNew route.
- `ModelExtrudeQuadTriNoNew.jl` and the catalog helpers encode real one-cell/two-Tri, Quad patches/strips, rectangular grids and B4 scopes. Do not generalize a catalog by choosing source UID or column chronology from sorted XYZ. An extension needs actual source-column/side Surface provenance, one joint shared Face product, diagonal compatibility and complete exterior/volume certificates.
- Boundary-layer Extrude and pipe frontends beyond ordinary Layers remain distinct missing frontends/production; start from `ModelExtrude.jl`, `GeoExec.jl` and existing `BoundaryLayer` kernels after a tiny actual native-success context is selected. Existing native boundary-layer kernels do not alone implement those frontends.

Suggested bounded layer successor: two actual adjacent NoNew sweeps sharing one literal lateral Surface, initially linear and translation only. Bind the shared actual Surface constructor once and transport its corner/column UID incidences into both Region plans. Certify opposite child incidences, matching side diagonals, positive typed volumes and external boundary. Keep other planner guards until their actual effects are implemented.

## Remaining tensor/automatic-sizing production gaps

Implemented and integrated packets:

- `C:/tmp/tessella_takeover_constant_metric_box_20261011/ready.json`:43green. `src/meshing/Mesh3DMetricBox.jl` gives explicit constant `Metric3` boxes. The original physical `hmax` method remains untouched.
- `C:/tmp/tessella_takeover_constant_metric_regions_20261011/ready.json`:40green. Shared literal Region planes, tags/last-priority/void semantics and full bounding-grid resource preflight; actual rounded gaps and all7Kuhn displacement classes certified.
- `C:/tmp/tessella_takeover_metric_surface_20261011/ready.json`:33green. `ModelMetricSurfaceProducer.jl` plus Model include/endpoint. Variable tensor, actual linear Plane lower authority, conforming interior bisection, preserved boundary cells and actual endpoint/midpoint edge-norm certificate.
- `C:/tmp/tessella_takeover_api_metric_surface_20261011/ready.json`:47green. `APIMetricSurface.jl` + narrow API/Higher composition exposes `API.mesh.generate_metric_surface` using the selected tensor background, field caps/factor once and atomic classified publication. The first-order retention and regeneration birth-row fixes are intentional; ordinary native branches remain unchanged.

Genuine gaps: general variable-metric Volume generation, boundary adaptation for the explicit tensor Surface route, curved/periodic/recombined/embedded tensor Surface contexts, and native automatic-sizing/octree equivalence. Existing Surface kernels and Volume `refine_to_size` already consume directional metrics. The ordinary background adapters scalarize in some generation routes; **algorithm6 is not thereby a confirmed native parity defect** (tiny SDK direction swaps retained the same mesh), while native7/BAMG changes interiors and its actual kernel is missing. Keep algorithm7 guarded; do not label the explicit bisection endpoint BAMG.

Paths: `src/fields/SizeFieldCatalog.jl` (`AutomaticMeshSizeField` around4377 is the documented discrete sphere-fit/gap/gradation analogue), `ModelMetricSurfaceProducer.jl`, `APIMetricSurface.jl`, `APIGenerateHigher.jl` field adapters, `ModelNativeHybridSizing.jl`, `Mesh3D.jl` and refinement consumers. First choose either boundary splitting with actual lower ownership transport or an explicit bounded tensor Volume endpoint. For the latter, reuse existing owned exact/refinement kernels only where they actually consume the tensor; certify actual directional edge samples, material conformity, positive volume, enclosure and resource limits. Do not guess a field is constant or reinterpret scalar Gmsh algorithms as anisotropic ones.

## Remaining partition production gaps

Implemented: bare Tet4 interface-only topology/ghosts; discrete Tri3/Tet4 lower BRep with exterior Surface children and complete signed Region facets; closed seams and nonbranching open seams preserving actual original Surface borders; genuine pre-existing Line2 Curve and Point15 Point BRep topology with one partition child per original entity, signed Surface→Curve→Point child closures, shared `partitions` from literal Tri3 carriers, unchanged original cell/node UIDs, physical memberships, declared Curve edges on interface borders and lowest-dimension node reclassification. Ghost copies are query aliases of one authoritative actual Tet, not duplicate material cells. Parents/partitions/physicals and actual original UID/connectivity are retained. Actual MSH4 `$PartitionedEntities`/`$GhostElements`, public import and unpartition are qualified. MSH2 partition export rejects before opening output. Default legacy assignment behavior for other mixed families remains compatible; labels alone are not topology production.

Frozen packets: `C:/tmp/tessella_takeover_partition_20261011/ready.json`, `...partition_brep_20261011/ready.json` (56green), `...partition_open_seams_20261011/ready.json` (28green). Current Source: `APIPartitionProducer.jl`, `APIPartitionNative.jl`, `APIPartitionBRep.jl`, `APINativeViewPartition.jl` and MSH writer/metadata transport. Docs: `docs/partitioned_tets.md`, `partitioned_brep.md`, `partitioned_open_seams.md`, `partitioned_declared_brep.md`.

Actual guards in `APIPartitionBRep.jl`: only linear Point15/Line2/Tri3/Tet4; declared Curve topology must be one nonbranching directed chain whose signed Point bounds own the chain ends' Point15 cells; declared border Curves must lie on their Surface's literal border with matching orientation and complete coverage; declared Surfaces/Regions must be bound; referenced vertices need consistent lowest-dimension classification; no periodic, parametric or already partitioned data; no seam branches/triple junctions; no geometric support that simultaneously represents an original Tri3 and a new interface Tri3. The last limitation is a real distinct-role representation problem. Preserve normal duplicate validation until a certified representation exists.

Landed (2026-10-12): the bounded increment above — actual pre-existing Line2 Curve and Point15 Point BRep topology is now retained. `generate` always materializes a vertex node plus Point15 cell for each dim-0 entity, so a missing Point carrier is only reachable through non-generated models; `A.open` normalizes original non-partition boundary signs to `abs`, which blocks partition certificate reuse on imported declared Curve bounds — both pre-existing behaviors.

Best next bounded increment: mixed/high-order partition cells or periodic partition production. Triple junctions, coincident-support roles and actual parallel partition algorithms remain separate later increments.

## Publication and continuation notes

- The current HANDOFF backlog records the connected shared-Tri translation scope. SI-unit transport52 and view representation replacement26 are integrated and qualified in the successor224 batch; exact current scope is in the main Devin checklist.
- Native high-order Plane/transfinite Hex production and catalogued tagged conversion1–10 already exist. Do not repeat the old "set_order only1/2" limitation for the native tagged path; legacy compact overlay remains1/2. Remaining nonlinear source/curved placement contexts are real scopes to extend.
- The independent prismatic fallback packet is integrated: `C:/tmp/tessella_takeover_prismatic_fill_20261011/ready.json`. `Mesh3DPrismaticFill.jl` recognizes paired actual cap supports/side diagonal precedence and accepts only through the unchanged full PLC fill certificate. It fixes the retained rotated IGES hole case; unsupported/cyclic/different-cap controls use the old fallback. Do not bypass the full certificate or restore the expensive path for recognized inputs.
- All previous failed attempts remain in candidate receipts. Their failure is not a release gate once an actual corrected final source/fixture is qualified. Historical dated campaigns and source-pin receipts are archival evidence, not instructions to rerun multi-hour matrices.
- Continue with one concrete producer increment and proportional public/analytic tests, preserving exact ownership and full relevant geometry certificates. No production Gmsh backend, restricted matching implementation Source, chronology guessed from sorted XYZ, silent empty-success, or metadata-only replacement for missing topology.

This file documents unimplemented designs and remaining gaps. No production Source was changed for those designs, and no new Julia or native process was launched after the user's handoff request. Previously running final ordinary tests completed before publication.

## Repository evidence

- [Rigid periodic controls](../../test/artifacts/handoff_20261011/native_rigid_periodic_probe.json).
- [Shared translated control](../../test/artifacts/handoff_20261011/native_shared_periodic_probe.json).
- [Invalid Quad-sheet observation](../../test/artifacts/handoff_20261011/native_quad_sheet_probe.json).

Only the shared translated feature has a Tessella implementation here. Rigid
rotation observations remain design evidence; native return status alone does
not establish a valid producer. No new code or runtime remains active.
