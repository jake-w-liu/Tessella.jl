# Development

Read [AGENTS.md](AGENTS.md) and the [requirement checklist](CHECKLIST.md) before
implementation. Resume from the current `main` checkout and preserve unrelated
work. Verified changes are normally committed and pushed to `main` without
force-pushing.

The active objective covers correctness, resource efficiency, complete native
meshing functionality, repository cleanup, and registration of version 0.1.0.
Full Gmsh 4.15.2 compatibility remains open. Each required domain must be
implemented and verified, or supported by evidence of an upstream/format
limitation or an intentional stricter contract. An unexplained rejection is
unfinished work.

Use small functional fixtures with independent analytic or topology checks.
Follow [test/README.md](../../test/README.md) for the default suite. Record actual
process exits and distinguish minimal compilation from normal optimized checks.
Full native sweeps, source/AST campaigns, allocation stress, benchmarks, HFSS,
and external solver campaigns are opt-in.

## Current evidence

### 2026-11-02 fresh correctness and release pass

Commits `1a722d5`, `91bf429`, and `4454182` are pushed. Verified implementation checkpoint:
`4454182b2f9a76391976a8c2f21e29535cc16590`. Clean detached checkouts of that
commit passed ordinary `Pkg.test()` with bounds checks and automatic minimal
compilation on both supported lines:

| Runtime | Assertions | Actual exit | Test body |
| --- | ---: | ---: | ---: |
| Julia 1.12.7 | 749/749 | 0 | 13m07.1s |
| Julia 1.13.1 | 749/749 | 0 | 12m48.9s |

Verification directories: `C:/Users/User/Tessella-verification-4454182`
and `C:/Users/User/Tessella-verification-4454182-113`; logs are
`C:/Users/User/Tessella-clean-4454182-112.log` and
`C:/Users/User/Tessella-clean-4454182-113.log`. Their tracked sources stayed
clean. This qualifies that exact published snapshot, not later main-worktree
edits or the whole completion objective. Registration remains open.

The following additional fixes have fresh scoped evidence:

- Model-based PostView backgrounds now transport owned full-order geometry and
  NodeData/ElementData/ElementNodeData from the actual step-0 source model.
  Curved search bounds use conservative polynomial/rational intervals with an
  FP/subnormal envelope; linear and high-order cells share shape precedence.
  Actual target meshing, aliases, sparse UIDs, tensor data, warped quads, cubic
  pyramids and large-coordinate bulges are covered by 31 public checks.
- Model plugin sampling now uses actual family/full-order bases and per-step
  owned source slots. It respects dimensions and UID authority, including
  inactive models. Per-run snapshots are reused and cleaned on nested/error
  exits. Julia 1.12.7 minimal: 49/49, exit 0; Julia 1.13.1 normal `-O2`:
  new and existing plugin checks 104/104, exit 0.
- Replaced a model-view export-plan comprehension with a typed loop. The same
  cached three-node check fell from 166.9 to 1.15 seconds, compiler allocation
  from 57.77 GB to 610.9 MB. Targeted writer/lifecycle checks: 202/202 on each
  supported Julia line, minimal compilation, exit 0. A normal 1.13 `-O2` cold
  check passed 4/4 but still took 335.6 seconds; optimized startup is not solved.
- Safeguarded high-order field inversion now backtracks Newton steps and checks
  full parameter correction plus actual map residual. A globally regular
  monotone Line6 formerly cycled and returned 1 instead of 0.5. Independent
  known-root, outside, normal-offset and pyramid checks pass: 1.12.7 minimal
  model/inverse 48/48; 1.13.1 minimal inverse 17/17; 1.13.1 normal `-O2`
  inverse/list/scalar/tensor checks 204/204. All exits 0. Three warmed inverse
  query paths allocate 0 bytes; existing warmed high-order checks also pass.

Older full checks reached a torus/box Boolean integration case and were
deliberately superseded (exit 1, not product-failure/pass receipts). That exact
case is now opt-in alongside full GEO adaptation/matching and NURBS volume
recovery; ordinary checks retain small box Booleans and CAD boundaries.
Initial main-worktree 445 checks also overlapped an uncommitted probe-helper
parse error and exited 1; they are not evidence about committed 445. The clean
checkout checks above replace those invalid/partial runs. Continue checking
immutable commits separately from parallel main-worktree implementation.

Canonical probing, thin geometry, and embedded metric transport have the
additional scoped evidence below. Confirmed raw-basis arithmetic and general
non-partition geometry-map cases are the next active fixes. The remaining
whole P1–P6 checklist and fresh independent audit remain open.

### Probe and embedded metric implementation checkpoint

The public canonical `view.probe(tag,x,y,z,...)` now shares full-order geometry
and value bases with the legacy model-name overload. It supports component
selection, gradients, dimensions, nearest-node distance limits, ordered element
coordinates, aliases and sparse/empty time steps. All-step model queries and
the public Probe plugin retain the first frame's material geometry and actual
node/element UIDs; explicit steps retain their own geometry. Thin regular
elements use scaled rank/solve checks instead of an absolute determinant cutoff.

Final-source scoped checks: Julia 1.12.7 minimal compilation, 205/205; Julia
1.13.1 normal `-O2`, 478/478; both actual exit 0. These include the new 106-check
probe fixture, seven thin-geometry assertions, and nearby model, plugin, inverse
and high-order fixtures. They do not qualify the subsequent arithmetic fixes.

Planar algorithm-7 meshing now transports actual embedded Point and Curve
products through grading, cavity refinement, recombination and full-quad
subdivision. Lower Curve edges remain protected with their required incidence;
supplied Points split the detached Curve chain before coarsening, preserving
their UIDs. New Curve middles reserve UIDs before interior births. Strict metric
and ownership certificates remain enforced. The new fixture passed 191/191
on Julia 1.13.1 minimal compilation before expansion; final-source normal
`-O2` checks passed 328/328 (251 new and 77 neighboring assertions), exit 0.
Julia 1.12 and whole-package verification of this increment remain pending.

### Bounded parallel dispatch increment

Inactive `@threads` expansions forced compilation of whole serial partition
producers. The shared bounded dispatcher isolates that machinery and preserves
real parallel callbacks. A cold two-Tet minimal-compilation check fell from
202.57 seconds / 58.44 GB total allocation to 0.17–0.22 seconds. These are
compiler/startup observations, not a steady-state mesh-memory or whole optimized
producer claim; the separate whole-fixture `-O2` diagnostic was deliberately
terminated after 498 CPU seconds without assertions and is not a pass/failure.

The same investigation confirmed `General.NumThreads=2` allowed four concurrent
callbacks in partitioning, structure validation, and quality scans. All now use
joined chunks capped by request, default-pool size and item count; zero selects
the system/default pool. Source extraction, ownership certificates and ordering
remain intact. New 53-check regression includes caps, empty/tiny ranges, joined
failures, and byte-identical surface/volume MSH4 products with ghosts.

Scoped checks on final dispatcher/Core source: Julia 1.12.7 and 1.13.1,
`-t4 --compile=min -O0 --check-bounds=yes`, 198 assertions each, actual exit 0;
Julia 1.13.1 `-t4 -O2`, 39 dispatcher/scan checks, exit 0. Actual callbacks ran
on multiple default-pool threads. The separate declared-surface cold compiler
hotspot still requires investigation; global optimization is not claimed.

#### Earlier checks in this pass

- Fixed rotation/reflection axes with subnormal or overflowing norms by scaling
  before normalization; ordinary-axis arithmetic is unchanged. Analytic
  coordinate, length, orientation, ownership, and involution checks:
  `test/core/transform_axis_scale_test.jl`, 35/35, exit 0 on Julia 1.13.1
  both `--compile=min -O0` and normal `-O2`, bounds checks enabled.
- Fixed triangle distance at large translations: the plane-interior distance
  uses a local signed residual instead of rounded world-coordinate projection.
  Uncertain cancellation and overflowed region tests use the existing exact
  fallback. `test/fields/distance_triangle_translation_test.jl`: 18/18,
  exit 0 on Julia 1.12.7 minimal compilation and Julia 1.13.1 normal `-O2`.
  Independent near-plane checks also passed (3/3). Warmed 100-query observations
  of ordinary 12-primitive point, segment, and triangle BVHs each allocated
  0 bytes on Julia 1.13.1 `-O2`; exact fallback allocations remain intentional.
- Direct-model Plane faces bounded by one closed Curve now share the native
  constructor's existing QR/Jacobi frame fallback in meshing and evaluation.
  Focused exact-final-source minimal check: 20/20, exit 0 on Julia 1.12.7
  (14.7 seconds). Normal `-O2` also passed 20/20, exit 0, in 7m32.8s including
  native dispatcher compilation; this run preceded removal of one unreachable
  check, with the same implementation. This closes the reproduced producer gap
  without weakening boundary or ownership certificates.
- The historical minimal-compiler GmshLibm diagnosis is stale: a fresh 16-check
  scalar/CLI task-local probe passed with exit 0 on Julia 1.12.7
  `--compile=min -O0 --check-bounds=yes`.
- Early ordinary package checks on Julia 1.12.7 and 1.13.1 were superseded. The
  initial attempt stopped before tests because the ignored local manifest
  lacked Sockets; it was backed up outside the repository and resolved. A
  superseded baseline and first 1.12 final-source run were stopped after
  replacement checks started; their exit 1 is deliberate termination, not a
  passing receipt. Per-file diagnosis confirmed accidental default volume
  recovery work on 130-vertex and 440-vertex NURBS boundaries. Those two fills
  remain in the opt-in geometry regressions; the default keeps boundary imports.
  Clean immutable-commit checks of that routine scope passed on both runtimes
  as recorded above.
- Prepared version 0.1.0, bounded stdlib compatibility, the owner-selected MIT
  license, and Windows/Linux CI for both supported Julia lines. Root tracked
  Markdown is now only README; feature scopes, full historical ledgers and
  prepared designs are retained below .github. Third-party notices remain
  unchanged under LICENSES. Registration is not complete; General's human
  maintainer-review prerequisite is also unverified.

The formerly missing model-record PostView adapter is implemented and verified
above. AutomaticMeshSizeField's translated closest-point helper remains a
separate unverified candidate. Full P1–P6 implementation and the fresh
whole-codebase/resource audit remain open.

The last repository state before documentation consolidation was
[`a3c4ac3`](https://github.com/jake-w-liu/Tessella.jl/commit/a3c4ac37af44bcfe30b5bc7976a284b65f53a9a4).
Its latest published changes fixed coincident/tangent Boolean materialization
and full-quad shared-endpoint transport. The
[historical handoff](history/HANDOFF.md) contains their dated test receipts.
Those receipts qualify the recorded revisions and scopes, not subsequent edits
or the whole completion goal. New fixes must carry fresh evidence.

## Preserved records

Detailed dated evidence remains available for checking specific requirements:

- [Complete previous checklist and prepared designs](history/DEVIN_HANDOFF.md)
- [Architecture and P1–P6 scope](history/PLAN.md)
- [Dated verification ledger](history/STATUS.md)
- [Implementation handoffs](history/HANDOFF.md)
- [CRC development discipline](history/DEVELOPMENT.md)
- [ASCENT integration and external campaigns](history/ASCENT.md)

These files are explicitly historical. Read their dates and current code before
reusing a diagnosis, test count, support claim, or blocker.
