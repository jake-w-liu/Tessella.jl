# Devin handoff — complete Tessella's remaining implementation

This is the current continuation checklist requested by the user on 2026-10-11.
**Tessella is not complete. All P1–P6 parity tracks remain open.** The original
Stage0–6 simplex baseline is complete; it is not the broader completion target.
Finish every open requirement below, update the ledger as work lands, and stop
when the complete authorized implementation is finished. Do not stop after only
the first two prepared designs or call a bounded test pass full parity.

## Start from the published repository

Use current `origin/main` in `https://github.com/jake-w-liu/Tessella.jl`.
The final batch follows parent `da92208`; its commit is recorded by Git rather
than a self-referential hash in this file. Read this checklist, [HANDOFF](HANDOFF.md),
[STATUS](STATUS.md), [PLAN](PLAN.md), [AGENTS](AGENTS.md) and [test instructions](test/README.md).
Resume from checked-out Main, not an old C:/tmp candidate. No implementation
agent or Julia reader remains active at handoff. Existing unrelated untracked
`nul` and `Ctmppel_repro1.jl` were preserved; do not remove them.

The final ordinary suite passes **224/224 on both supported runtimes**, actual
exit 0: Julia 1.12.7 **45.21 s**, body 11.5 s; Julia 1.13.1 **38.38 s**, body 10.6 s.
Both read the same 252 Source files, aggregate raw SHA256
`6174da34f24ab27c6c3c0afe214d249b78e0a2b564eadb5f7101e9d34d3b2bf7`.
Portable [Source map](test/artifacts/handoff_20261011/verified_source.json)
includes normalized LF hashes for fresh checkouts. Actual receipts:
[Julia1.12](test/artifacts/handoff_20261011/ordinary_default_224_112_20261011.receipt.json)
and [Julia1.13](test/artifacts/handoff_20261011/ordinary_default_224_113_20261011.receipt.json).
Their separate stdout/stderr logs are in that same committed artifact directory.
The final documentation audit covers all 70 repository Markdown files:
335 relative file links resolve, with no missing files or heading anchors.
Current instructions supersede the explicitly dated historical evidence.

## Work already finished — preserve it

The latest batch implements and qualifies:

- Direct STEP product-usage graphs and [SI source/placement length transport](docs/step_relationship_si_units.md):52 focused checks. Each source and both frame origins use their actual context factors; ordinary/mapped-root units remain separate work.
- [Connected hybrid Regions with shared Tri interfaces and translated exterior periodic pairs](docs/shared_periodic_hybrids.md):63 focused checks,21Tet4+5Pyramid5 and positive volume1 per Region, actual primary/P2 maps and MSH4 identity transport.
- [Atomic list/model view replacement](docs/view_representations.md):26 checks, obsolete representation rows removed, display/tag/name/held models preserved, serialization valid. Two cheap checks enter the default.

Earlier landed features include atomic GEO/MSH/plain-linear POS open/merge,
current-model MSH export, actual view-time metadata, list/model CSV/VTK,
complete catalogued high-order VTK layouts, unit-scale Cartesian STEP targets,
recognized independent CAD roots and planar holes/IGES143 boundaries,
explicit IGES190 UV charts/rigid124 chains, rigid-prism fill, Point/Curve
collections and endpoint-only3–6-vertex Tri3 sheet recovery, Pyramid relocation/
Winslow and Region batches, tensor box/shared-material grids, explicit variable-
tensor Surface/API construction, and Tet/Tri3 partition topology with ghosts and
closed/open seams. Consult their docs before reclassifying them as absent.
The known lost-sheet optimizer key bug and model-view kind-reset bounds bug
are fixed. Linear cell factory renewal and optional renumbering are intentional.

## Next two prepared implementation designs

1. **Nested STEP product usages.** [Design with exact blockers](docs/handoff/nested_step_assemblies.md). Two portable STEP inputs and opaque public observations are committed under `test/artifacts/handoff_20261011/nested_relationships`. They demonstrate noncommuting rigid transforms and root-m/intermediate-cm/leaf-mm origins. No nested Tessella implementation exists. Extend the complete reachable acyclic product graph, validate each edge's actual parent/context/frame, compose physical transforms, scale leaf geometry once, and emit separate root-to-leaf occurrence identities. Preserve full graph, shell, separation, volume and resource certificates.
2. **Proper rigid periodic hybrid pair.** [Meshing design and remaining guards](docs/handoff/meshing_remaining.md). Portable [Rz90/Rz180 observations](test/artifacts/handoff_20261011/native_rigid_periodic_probe.json) show valid24Tet4+6Pyramid5, volume1 and4/9 primary/P2 pairs. No rotated Tessella implementation exists. Repair actual Surface-derived cross-side Curve planning using existing `ModelPeriodicConstraint.derived` provenance, then replace only the hybrid translation gate with a certified proper-rigid admission. Do not bypass direct Curve, copy-history, ownership or final topology checks.

These are design/evidence packets, not already qualified Source. Review and
implement isolated deltas before promoting them. Local C:/tmp originals are
provenance; the repository copies allow a different machine/session to resume.

## Complete remaining checklist

Each row is a required work area, not a single trivial edit. Maintain subitems
for the actual pinned-release requirements as they are identified. Expand the
ledger rather than silently narrowing the target. Prioritize meshing/CAD value.

| ID | Open implementation and completion requirement | Starting Source / evidence |
|---|---|---|
| CAD01 | Nested product-usage DAGs with repeated paths, dimensional edge transforms, complete reachability/cycle/orphan checks and separate occurrence ownership. | `BRepStepRelationships.jl`; prepared nested design above. |
| CAD02 | Mixed mapped/product workflows, multiple roots/solids per leaf, geometry plus subassembly usages, and general occurrence transforms with actual context/item provenance. | `BRepStepRelationships.jl`, `BRepStepOccurrences.jl`, `BRepAssembly.jl`. |
| CAD03 | Ordinary/mapped-root unit conversion, conversion-based length/angular units, compatible uncertainty aliases and target-unit contexts. Direct relationship SI scaling already works. | `BRep.jl`, `BRepStepOccurrences.jl`, `BRepStepRelationships.jl`. |
| CAD04 | Curved trims, periodic/general NURBS topology, nonlinear UV charts and complete imported analytic Face/Curve parameterization/registry. Mesh import alone is not a CAD kernel. | `BRep.jl`, `BRepPlanarFaces.jl`, `BRepPlanePcurves.jl`, Model CAD/evaluation helpers. |
| CAD05 | Contact/touching/overlapping assemblies with genuine interfaces; arbitrary CAD transformations, Booleans/imprints/intersections beyond materialized supported primitives. Current disjoint-box admission is a sufficient guard, not a contact producer. | `BRepAssembly.jl`, `Model.jl`; PLAN's oblique/non-axial/curved intersection guards. |
| GEO01 | Remaining `.geo` language and genuine `model.geo`/`model.occ` frontend coverage, allocator/synchronization/lifecycle and full CAD execution. Reuse implemented Model operations before adding wrappers. | `GeoExec.jl`, `Model.jl`, `API.jl`. |
| PER01 | Proper-rigid rotated hybrid copies with certified derived cross-side Curve transport. | `Model.jl::_surface_periodic_constraints`, `ModelNativeHybridVolume.jl`; prepared rigid design. |
| PER02 | Multiple periodic directions/pairs, general affine and cross-Region/shared periodic Face contexts, consistent corner/edge relations and embedding combinations. Shared-Tri components with uniquely owned exterior translated pairs already work. | `ModelNativeHybridBatch.jl`, `ModelPeriodicUIDHistory.jl`, Surface producers, API periodic history/order helpers. |
| EMB01 | Larger/subdivided/interior-vertex/holed/multiple/curved sheets and supported mixed sheet families; nonlinear Curve and broader multi-Region embedded composition. | `ModelNativeHybridEmbeddedSurface/Curve/Constraints.jl`, Recovery/Carve/Insertion/Sizing; retain exact facet/star/link certificates. |
| ALG01 | Remaining successful native ordinary surface/volume algorithms and mixed-block unstructured integration/recombination. Enumerate algorithm dispatch requirements rather than removing unsupported selectors. | `MeshSurface.jl`, `Mesh3D.jl`, `ModelSurfaceAlgorithms.jl`, `APIGenerateHigher.jl`. |
| MET01 | Native anisotropic algorithm7/BAMG producer. The explicit conforming tensor bisection API is implemented and must not be renamed BAMG. | `ModelMetricSurfaceProducer.jl`, `APIMetricSurface.jl`, SizeField modules. |
| MET02 | Variable-metric Volume generation; tensor boundary adaptation and curved/periodic/recombined/embedded Surface contexts with actual lower ownership transport. | Metric producers, `Mesh3DMetricBox.jl`, refinement and hybrid sizing consumers. |
| MET03 | Full native automatic-sizing/octree behavior. Curvature/gap/bulk/gradation stages and a discrete analogue exist; HXT/p4est equivalence does not. | `SizeFieldCatalog.jl::AutomaticMeshSizeField`, field/background context adapters. |
| LAY01 | General NoNew source grids/categories and neighboring/shared-Region propagation, periodic/compound/reversed/constrained boundaries, nonlinear retained sweep sources. | `ModelExtrudeNoNewScope.jl`, NoNew/QuadTri catalog helpers, `APIGenerateHigher.jl`; meshing handoff has concrete shared-sweep plan. |
| LAY02 | Boundary-layer Extrude and pipe frontends beyond existing ordinary Layers; remaining curved layer/recombination contexts. | `ModelExtrude.jl`, `GeoExec.jl`, `BoundaryLayer.jl`. |
| HO01 | Remaining curved/mixed/hybrid high-order generation, projection/refinement placement and global Jacobian certification beyond the current globally certified domains. Native tagged high-order generation/order conversion already exists. | `HigherOrderNodal.jl`, Model high-order helpers, API order/refine producers; preserve every UID/node. |
| OPT01 | Mixed Surface/nonlinear/primitive relocation, general high-order/hybrid optimization and independent Netgen-equivalent behavior. Tet4/Tet10 Relocate3D and bounded optimizers are existing scopes. | `APINativeOptimize/Relocate3D.jl`, hybrid optimizers, `Mesh3DGmshOptimize.jl`. |
| PART01 | Partition children/interfaces/ghosts with actual original Point/Curve classification, mixed/high-order cells, periodic/parametric topology, triple junctions and distinct coincident support roles. | `APIPartitionBRep/Producer/Native.jl`, `APINativeViewPartition.jl`; preserve parents/UIDs and import/unpartition. |
| PAR01 | Actual parallel producer and partition workflows. Deterministic query task slicing does not implement parallel meshing. | Model/mesh orchestration and partition producers; P5 scope in PLAN. |
| IO01 | Additional true model filename imports/merges/exports and long-tail formats; bridge real standalone STL/CAD readers into session lifecycles without fabricating models. | `APIModelExport.jl`, `APIModelImport/ModelMerge.jl`, format readers/writers. |
| IO02 | Lossless higher/incomplete pyramid and remaining high-order target mappings, including Pyramid14 VTK only if a truthful representation exists; pyramid model-view interpolation beyond P2. | `APINativeModelViewFormats.jl`, `APINativeViewWriter.jl`; keep no-node-dropping guard. |
| VIEW01 | POS text/interpolation/high-order session storage/query/writer transport. The standalone parser already supports richer structures; plain24 session import and view transitions are done. | `APIListViewImport.jl`, `_SessionView`, `PostViewIO.jl`, view writers. |
| VIEW02 | Remaining view endpoints: homogeneous model data, list strings, interpolation matrices, aliases, combine/probe and per-window visibility; preserve real owned data and time semantics. | `API.jl::view`; official pinned SDK wrapper/API catalog. |
| CLI01 | General conversion, partition/adaptation/post action frontends. GEO generation, sizing/order, help/version, binary MSH and read-only checking are done. | `CLI.jl`, true model/mesh/post producers. |
| UI01 | Gmsh-style plugin options/run catalog and remaining parser/onelab/graphics/fltk session services; actual interactive GUI/post controls beyond headless command state. | `Post.jl`, `GUI.jl`, absent API namespaces; P5 is required, not a non-goal. |
| MATCH01 | Remaining native matching/optimizer chronology and tie contexts through independent algorithms and opaque public evidence. Preserve the existing seeded odd-cycle fixes. | `OrderedMatching.jl`, `OrderedFractionalMatching.jl`, `NativeRecombineMatching.jl`, Model matching context; never inspect restricted native matching code. |
| P601 | Requirement-by-requirement tutorial/API/option/format corpus and differential conformance across P1–P6. Record every missing endpoint/guard as a classified requirement; include remaining quality/basis/query/lifecycle contexts. | PLAN P1–P6, official pinned public API/manual, existing bounded regressions. |

## Distinguish missing features from invalid witnesses

- Shared-Quad two-Region Pyramid finalization is a pinned native unsupported case.
- The native Quad-sheet probe returned only6 boundary Pyramids, zero Tets and no interior typed carrier. It is not evidence of successful constrained meshing. [Observation](test/artifacts/handoff_20261011/native_quad_sheet_probe.json).
- Region reverse/smoothing and certain Cartesian-scale observations were native no-ops. Do not invent a missing successful feature from their spelling.
- Native formats cannot encode every embedding/physical/variable-connectivity relation. Preserve supported actual records and precise format limits; never drop nodes or fabricate metadata.
- Atomic import/conversion errors intentionally improve on native partial publication. Preserve that contract.
- The tiny PostView display-selection control leaves its first-step field mesh unchanged. It does not justify changing field evaluation or proving broader field parity.
- A historical failure or missing receipt is not proof of an absent implementation. Reproduce the actual current public context before changing Source.

## Verification, publication and finish criteria

1. Use independent Julia production. Gmsh4.15.2 is an opaque public comparator, never the production mesher. Do not read/port restricted Gmsh/Concorde matching implementation.
2. Stage real producer data and publication. Preserve ownership, literal UIDs/ordered cells/coordinates, actual construction/copy histories, allocator clocks and relevant geometry/topology/volume certificates. Implement unsupported transport rather than weakening certificates.
3. Keep default tests small public functional checks. **No multi-hour native sweeps, source/AST campaigns, memory stress, benchmarks, HFSS or solver runs in `Pkg.test()`.** Feature-specific tiny native/volume readers stay explicit opt-in. Do not rerun unchanged campaigns simply for a handoff or MD edit.
4. For changed final Source, run ordinary `Pkg.test()` on Julia1.12 and1.13 against the same Source, record actual child exits and scoped meaningful checks. See test/README; normal invocation automatically chooses minimal compilation. Do not add elapsed-time assertions that become flaky on another host.
5. Review and normally commit/push verified implementations to **main**, never force-push. Verify exact remote SHA. Preserve unrelated work. Update all current Markdown and keep dated failures/history labeled historical.
6. Maintain a requirement ledger under each checklist ID. A gate is either a verified implemented domain, demonstrated native/format limitation, intentional documented stricter contract, or a real open implementation. Do not move required features into non-goals to make the list look finished.
7. Complete all actual P1–P6 requirements/subitems and relevant conformance before marking Tessella complete. Publish final Source/docs, ensure no required work remains, then stop. Large ASCENT/HFSS external campaigns remain separately opt-in; their old exhaustive gates do not replace the user's concise-test instruction.

Suggested instruction to Devin: **Read this file, continue from origin/main,
complete every open requirement and newly confirmed gap, preserve concise default
tests, push verified implementations to main, keep Markdown current, and stop
only after the full implementation is complete.**
