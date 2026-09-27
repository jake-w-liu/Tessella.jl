# Handoff — Tessella.jl work state

Continuation instructions for resuming this work on another machine.
Branch: `main` (this state is pushed). Goal: independent Gmsh 4.15.2 parity —
never use Gmsh as the production mesher; it is only a differential oracle.

## Environment

- Julia compat: `1.12 - 1.13` (Project.toml). Verified on 1.12.7, 1.13.0, 1.13.1.
- Pinned oracle: Gmsh 4.15.2 (`/opt/homebrew/bin/gmsh` on this machine;
  install 4.15.2 on the new machine for validation differentials).
- Run tests: `julia --project=. --check-bounds=yes -e 'using Pkg; Pkg.test()'`
- Validation driver: `julia --project=. validation/run_all.jl`

## What this push contains (increment just landed)

1. **`Mesh.MeshSizeExtendFromBoundary` integer semantics** — storage moved
   Bool→Int; documented `0` (disable) now accepted and distinguished from
   absent; float inputs truncate toward integer; Bool/NaN/Inf/out-of-range
   reject without mutating prior state. Files: `ModelMeshingAttributes.jl`,
   `GeoExec.jl`, `API.jl`, `GeoOptionTables.jl`, `ModelIdentity.jl`, tests in
   `geo_constraints_test.jl`, `api_refinement_classification_test.jl`.

2. **Reversed-volume refinement classification fix** — `Refine.jl`
   `_write_positive_tet!` → `_write_oriented_tet!`: children preserve parent
   orientation (checked via `orient3`) and `tet_tag`; the classification
   skeleton in `API.jl` carries `segs`/`tris`/tags through instead of
   discarding them. `refine_uniform` gained `require_positive_tets=false`.
   Covered by `api_refinement_classification_test.jl` (18 tests).

3. **Embedded curve/surface recovery overhaul** (`Mesh3D.jl`, ~1800 lines) —
   `recover_segment3`/`recover_triangle3` now converge on the hard sheet
   fixture (triangle at z=0.5, verts (0.2,0.2),(0.8,0.2),(0.5,0.8)): protected
   face/edge registries per iteration, multi-pocket batch dispatch, geometric
   dedupe of detector output, strict-enclosure grafting, foreign-vertex
   seed-face rejection (vertex-on-edge AND vertex-on-face interior tests),
   unseeded pool fallback, cavity-growth retry (`_refill_with_growth3`),
   monotone-coverage acceptance gate (`_accept_sheet3`), `soft_keepfaces`
   mode for sheet fills, and `_snap_to_plane3` for 1-ulp off-plane Steiner
   vertices. Exact coverage certificate + `validate(mesh)` enforced.

4. **Cross-version/per-process determinism** — sorted every hash-order-
   sensitive iteration that influences output: `refine_to_size` initial edge
   queue + deferred requeue, `RecoverCDT.build_regions` component order,
   missing-crease edges, piercing-candidate edge set, `delaunay3d_exact`
   cavity faces, `_refill_segment_cavity3` boundary faces + `claimseq` DFS
   ledger, `_grow_cavity3` incidence, `_sheet_chain_edge_gaps3` edge set.
   **Result: byte-identical coords/tets across Julia 1.12.7 and 1.13.x** on
   every fixture (direct 5474 tets, classified 5934, holed, .geo 1067).

5. **Task-local registry scoping fix** — `_protected_faces3`/`_protected_edges3`
   live in `task_local_storage`. `mesh_model_volume` is now a thin wrapper
   (`_mesh_model_volume` is the body) that empties the registry on entry and
   in `finally`. Previously a registry left populated by an earlier call in
   the same task leaked foreign keep-constraints into `refine_to_size`'s
   cavity fills (via `_repair_steiner_split3`→`_refill_segment_cavity3`),
   changing raw cell order → same canonical `mesh_crc` but different
   `mixed_crc` per process. Verified: a deliberately poisoned registry now
   produces the expected `ffd2559d…` CRC.

## Verified results (this machine, this code)

- `test/meshing/mesh3d_test.jl`: **146,746/146,746 on both 1.12.7 and 1.13**
- `test/core/allocation_audit_test.jl`: **76/76** (closure-boxing clean after
  extracting `tryaccept`→`_accept_sheet3`, `onplane_host`→`_onplane_host3`,
  and the `sortperm` comprehension → explicit loop)
- `test/interfaces/api_refinement_classification_test.jl`: 18/18
- `test/geometry/geo_constraints_test.jl`: all green
- Fixture CRCs stable across versions: fine `7d82449b…`, sized `e314a5a0…`,
  geom-expr mesh `7527493c…`, mixed `ffd2559d…`, list-var `14168011…`,
  dyn-tag `75365136…`, set-max `d5d07bb1…`
- `.geo` embedded sheet: `validate=true`, `covers=true`, 1067 tets

## Pending gate before calling this increment done

**A fresh full `Pkg.test()` on the final code under Julia 1.13 was never
completed** — the run was killed because the machine was under memory
pressure (many parallel suite processes). Prior suite runs showed only
stale-build artifacts (failures from processes that precompiled before the
last edits): every such failure was verified to produce the expected value
on the final code. Run it fresh first thing:

```bash
cd Tessella.jl
julia +1.13 --project=. --check-bounds=yes -e 'using Pkg; Pkg.test()' 2>&1 | tee /tmp/pkgtest.log
```

If a `mixed_crc`/`mesh_crc` expectation fails, FIRST verify the value is not
a stale-build artifact: recompute the CRC in a fresh process on current code
before touching the expectation. The geo `mixed_crc` pins (`ffd2559d`,
`14168011`, `75365136`, `d5d07bb1`) were all confirmed correct on both
versions — do not change them.

Also rerun under 1.12.7 when convenient (`julia +1.12`).

Then write the dated entry in `STATUS.md` "Verification history" (the support
statement already says 1.12.x and 1.13.x — keep it accurate: only claim
full-suite green once it actually is).

## Remaining parity work (PLAN.md — all IN PROGRESS tracks)

- **P1**: boundary-layer 3-D multi-wall fans, full Gmsh automatic-sizing
  pipeline, high-order/custom-interpolation size fields, `PostView`
  metric/tensor fields, exact CAD distance.
- **P3**: general OpenCASCADE BREP kernel, NURBS CAD of unclassified
  topology, transforms of arbitrary CAD entities, full `.geo` execution,
  unrecognized CAD topology.
- **P4/P5/P6**: mixed-element generation/recombination beyond first-order
  surface pairing, non-simplex hierarchical bases (Pyramid/Trihedron),
  partitioning/parallel paths, views/plugins depth, CLI/GUI/post-processing,
  long-tail formats, and the standing requirement-by-requirement differential
  corpus vs Gmsh 4.15.2.

Workflow per user instruction: implement each parity item, **deep-debug it**
(verify against pinned Gmsh 4.15.2 where applicable, keep output
deterministic across Julia versions — audit every Dict/Set iteration that
can reach output), run focused tests + allocation audit, then commit/push.

## Gotchas learned the hard way

- **Stale precompile artifacts**: a test process that started before an edit
  reports failures that vanish on a fresh process. Always re-probe in a new
  process before classifying a failure.
- **Hash-order leaks**: `for (k,v) in dict` / `for x in set` anywhere near
  mesh output = cross-version AND per-process (hash seed) nondeterminism.
  `mesh_crc` canonicalizes so it can hide raw-order diffs; `mixed_crc`
  hashes raw order. Sort only iterations that influence decisions/output.
- **Task-local state**: the protected-cell registry must stay driver-scoped
  (see `mesh_model_volume` wrapper). Any new task-local/global mutable state
  needs the same scoping discipline.
- **Boxed closures**: `allocation_audit_test.jl` flags closures capturing
  reassigned locals — extract to module functions or precompute.
- **Memory**: do NOT run multiple full suites in parallel — each Julia test
  process peaks at several GB; that is what OOM-killed the machine.
- Debug instrumentation is fully removed; `DBG_*`/`REFILL_NOCLEAN` env vars
  no longer exist.
