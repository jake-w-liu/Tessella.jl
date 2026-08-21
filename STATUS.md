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
| P2 | **IN PROGRESS** | 125 fixed-node Gmsh types and local ordering, mixed blocks/entities/classification metadata, structural validation/CRC, and ASCII/binary MSH v2.2/v4.1 read/write |
| P3 | **IN PROGRESS** | Native analytical surfaces/imprints, closed box/cylinder/cone/geodesic-sphere primitives, cavities, mesh Boolean CSG, and finalized-mesh affine transforms |
| P4–P6 | **PENDING** | No state change |

P1 does not claim boundary-layer element topology, the full Gmsh automatic-sizing
pipeline, non-simplex/higher-order/vector/tensor `PostView`, direct tensor or
metric-meshing parity, full `.geo`/CAD-model execution, or exact CAD distance. P2 does
not claim element generation/recombination, variable-connectivity/internal types,
integration of mixed blocks into the simplex meshing kernels, curved high-order
Jacobian certification, ancillary/unknown-section preservation (binary readers reject
unsupported sections explicitly), repeated `$Nodes`, non-8-byte binary data, internal
indices beyond `Int32`, or lossless multi-physical-group MSH v2.2 projection. Some
registered fixed tags require explicit Tessella-only output because
Gmsh 4.15.2 cannot re-import them. P3 does not yet claim a general entity kernel,
OpenCASCADE/BREP/NURBS, CAD import/export, transformations of analytical/CAD
entities, or full `.geo` execution.

OpenCASCADE/NURBS, remaining algorithms/fields, broad formats and API, GUI, and
post-processing are unfinished parity tracks, not project non-goals.

## Current worktree verification

Verified on 2026-08-21 with Julia 1.12.7:

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
  regressions in `test/hfss_cases_test.jl`. Re-running the remaining literal cases as
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
