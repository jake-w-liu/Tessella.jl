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

- Implemented native scalar/size field separation and Gmsh-compatible `Distance`,
  `Threshold` (linear, sigmoid, and `StopAtDistMax`), `Box`, `Ball`, finite
  `Cylinder`, `Frustum`, `Min`, `Max`, and final clamp/factor semantics.
- `DistanceField` evaluates discrete points, segments, and triangles; generic distance
  values can be zero without violating the positive final-size contract. Its
  deterministic AABB hierarchy was differentially checked against exhaustive
  primitive scans and allocates zero bytes per warmed query in the enclosure graph.
- `refine_to_size` and `mesh_sized` now accept spatial size fields and verify every
  output tet edge against the minimum field value at its endpoints and midpoint.
- Focused field tests, existing 1-D/surface tests, the full Mesh3D suite, the
  installed-Gmsh 4.15.2 differential, and the full bounds-checked package gate are
  green.

OpenCASCADE/NURBS, mixed elements, remaining algorithms/fields, broad formats and API,
GUI, and post-processing are now unfinished parity tracks, not project non-goals.

## Current verified parity gate

Re-measured on 2026-08-19 with Julia 1.12.7 before the increment was pushed to
`origin/main`. Both bounds-checked package runs matched:

- `julia --project=. --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 152,789/152,789 assertions passed twice (3m26.3s, then 3m28.1s).
- `julia --project=. validation/size_fields/differential.jl`
  — `SIZE_FIELD_DIFFERENTIAL_OK` against installed Gmsh 4.15.2 for
  `Distance`→`Threshold`, `Box`, `Ball`, finite `Cylinder`, and `Frustum`
  line-mesh grading.
- Fresh-process `using Tessella` plus public `mesh_volume` / `size_at` /
  `mesh_sized` succeeded twice on `box_surface(0,1,0,1,0,1)`: volume fill
  9 nodes / 12 tets, `size_at(BoxField(...),0.25,0.5,0.5)==0.5`, field-sized
  mesh 60 nodes / 167 tets, both `validate` ok.

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

The 2026-08-14 stress measurement used
`refine_to_size(mesh_box(0,3,0,3,0,3; hmax=3), 0.25)` after one warm-up and `GC.gc()`.
Both versions produced 22,065 nodes, 98,304 tetrahedra, and a valid mesh:

| Version | Allocated bytes | Timed run |
|---|---:|---:|
| pre-fix `5d3e466` source snapshot | 75,038,688 | 0.3517 s |
| previous hardened implementation | 70,337,184 | 0.2972 s |
| current field/refinement increment | 100,434,688 | 1.1318 s (concurrent enclosure load) |

The historical hardened version reduced allocation by 6.27% from `5d3e466`. The
current increment is 30,097,504 bytes (42.79%) above that version on the same mesh;
this is an open performance regression, not a resolved result. The refinement queue
keeps one live heap record per long edge via a companion set, preventing duplicate
queued records; the benchmark verifies the resulting mesh rather than relying on
allocation alone.

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

All code claims above are tied to the commands and measurements in “Current verified
gate”; historical exploratory findings were removed from this live status file to
avoid presenting superseded diagnoses as current state.
