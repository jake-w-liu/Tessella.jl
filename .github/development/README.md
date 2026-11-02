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

### 2026-11-02 fresh correctness and release pass (in progress)

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
- Full ordinary package checks on Julia 1.12.7 and 1.13.1 are running. The
  initial attempt stopped before tests because the ignored local manifest
  lacked Sockets; it was backed up outside the repository and resolved. A
  superseded baseline and first 1.12 final-source run were stopped after
  replacement checks started; their exit 1 is deliberate termination, not a
  passing receipt. Per-file diagnosis confirmed accidental default volume
  recovery work on 130-vertex and 440-vertex NURBS boundaries. Those two fills
  remain in the opt-in geometry regressions; the default keeps boundary imports.
  Fresh package checks of that routine scope are now running on both runtimes.
- Prepared version 0.1.0, bounded stdlib compatibility, the owner-selected MIT
  license, and Windows/Linux CI for both supported Julia lines. Root tracked
  Markdown is now only README; feature scopes, full historical ledgers and
  prepared designs are retained below .github. Third-party notices remain
  unchanged under LICENSES. Registration is not complete; General's human
  maintainer-review prerequisite is also unverified.

Next confirmed implementation gap: a separate model's three-node constant
NodeData view used as a public PostView sizing field rejects with `view has no
element records` instead of returning 0.5 (Julia 1.13.1 reproduction, exit 0).
An owned model-record adapter is being implemented. AutomaticMeshSizeField's
translated closest-point helper remains a separate unverified candidate. Full
P1–P6 implementation and the fresh whole-codebase/resource audit remain open.

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
