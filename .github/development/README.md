# Development

The active goal is a complete native Julia mesher, fresh correctness and
resource audits, a clean repository, and registration of version 0.1.0.
The [requirement checklist](CHECKLIST.md) preserves all required domains and
links to their detailed specifications. Full Gmsh compatibility and the final
independent audit remain open. The user revoked the previous AGENTS.md
instructions on 2026-11-06.

Use [the test guide](../../test/README.md) for ordinary and extended checks.
Verify a concrete failure before fixing it, retain geometry/ownership
certificates, and qualify results by source, runtime, compiler mode and actual
process exit. Production meshing must remain independent of external meshers.
Push verified changes normally to main and verify the remote commit.

## Verified checkpoints

These receipts establish their stated scopes. A combined final package run is
still required after the current implementation work.

| Source | Verification |
| --- | --- |
| `e18b989` Opencode checkpoint | Default `Pkg.test()`: 1638/1638 on Julia 1.13.3; an independent detached checkout also passed 1638/1638 on Julia 1.12.7 (22m53.9s). Both use minimal compilation and bounds checks, actual exit 0. [Windows/Linux CI on both Julia lines](https://github.com/jake-w-liu/Tessella.jl/actions/runs/42519129163) also passed. |
| `11aee42` finite-normal inverse | Small intrinsic parameters beside ordinary coordinates: Julia 1.12.7 minimal, 162/162; Julia 1.13.3 normal `-O2`, 352/352 including existing warmed zero-allocation checks. Bounds enabled, exits 0. |
| `fc0f599` GEO/view publication | Julia 1.12.7 and 1.13.3 minimal each passed the 256-check view scope and 162-check import/merge/CLI scope. Normal Julia 1.13.3 passed the 52 new regressions. Bounds enabled, exits 0. The first normal plugin block took 11m13s; no compiler improvement is claimed. |
| `8a1f177` extreme nonlinear inverse | Subnormal derivatives and mixed-scale residuals: Julia 1.12.7 and 1.13.3 minimal, 168/168 each; Julia 1.13.3 normal `-O2`, 358/358 including existing ordinary warmed zero-allocation checks. Bounds enabled, exits 0. Exact recovery can allocate. |
| `fbf3451` factory Box sizing | Five production assertions passed on both Julia lines with minimal compilation. The unit box retained its boundary facets and volume, produced 1350 valid tets, and had zero oversized edges among 1818 (largest normalized length 0.98943). The same fixture on clean `e18b989` failed two assertions, including normalized length 1.63031. |
| Simplex-coordinate increment | Julia 1.12.7 and 1.13.3 normal `-O2` each passed 1026 coordinate/reference/locator assertions; minimal functional checks passed 80/80 each. A separate validation, quality, topology and transform audit passed 3124/3124 on each runtime. Bounds enabled, actual exits 0. |
| Triangle distance and automatic sizing | Julia 1.12.7 minimal passed 71/71 regressions. Julia 1.13.3 normal `-O2` passed 39,643/39,643 regressions, neighboring distance/automatic-constructor checks and resource checks. Bounds enabled, actual exits 0. |

Triangle region predicates now recover uncertain cancellation and underflow
through the existing exact path. Automatic sizing retains local closest-point
displacements and certifies scaled feature cones, including nearly parallel
tangents and the 60-degree boundary. Sampled ordinary point/segment/triangle
BVHs, automatic closest-point and gap queries each allocate zero warmed bytes;
the exact long-edge fallback intentionally allocates 30,448 bytes. The normal
neighbor scope stops before the unrelated GEO frontend campaign; unsuccessful
broad compiler/interpreter runs are excluded from these receipts.

Linear-simplex inversion now filters cancelled projection/Cramer numerators and
subnormal denominators before using the existing exact solve. Tiny nonzero
coordinates survive underflow in intermediate products. Explicitly binding the
scaled-frame vararg length also removed pre-existing triangle/tetrahedron
query allocations: sampled ordinary segment, triangle and tetrahedron loops
each allocate zero bytes under normal optimization. Exact fallback can allocate.
Unchanged extended allocation/digest checks assume normal optimization and
were qualified in that mode; the minimal receipts cover functional checks.

The nonlinear inverse now certifies the remaining intrinsic Newton correction,
damps exact oversized directions before Float64 conversion, and resolves
rounded merit-norm ties with exact squared norms. The regressions exercise
Line, Quadrangle and Tetrahedron maps and public view probing.

GEO plugin publication preserves aliases, explicit replacement tags, actual
model owners and imported node UIDs. CLI combination transports existing
PostView tag/index selections. The import rollback fixture uses explicit fatal
`Exit 2` and retains all state assertions; its old recoverable-error fixture
also failed 66/74 on clean `e18b989`.

Factory Boxes again extend boundary sizing without an explicitly selected
field. The private boundary view caches at most 4096 exact query results under
a lock; user fields and callbacks are evaluated independently. Cache checks
passed 13/13 on Julia 1.12.7 minimal and 14/14 on Julia 1.13.3 normal `-O2`,
both with four threads and bounds checks. The normal warmed query loop allocates
zero bytes. A 640-query minimal-compilation observation fell from 5.63s/298MB
total allocation to 0.65s/31MB with identical values. The production regression
took approximately 149s/147s on Julia 1.12.7/1.13.3. These observations do not
establish whole-mesher memory or performance parity.

## Local verification receipts

The following logs under `C:/Users/User/` have corresponding persisted `.exit`
files. They are local artifacts, not package dependencies:

- `Tessella-default-1133-final6-20261109.log` and `Tessella-verify-e18b989-112-20261109.log`.
- `Tessella-field-audit-gate-optimized-20261109.log` and `Tessella-field-focused-gate-112-min-20261109.log`.
- `Tessella-view-scoped-gates-{1127,1133}-min-20261109.log`, `Tessella-view-import-cli-{1127,1133}-min-20261109.log`, and `Tessella-view-final-gates-1133-normal-20261109.log`.
- `Tessella-field-audit-gate-optimized-overflow-20261109.log`, `Tessella-field-focused-gate-112-min-overflow-final-20261109.log`, and `Tessella-field-focused-gate-113-min-overflow-20261109.log`.
- `Tessella-box-production-sizing-112-final-20261109.log`, `Tessella-box-production-sizing-1133-20261109.log`, `Tessella-boundary-cache-112-min-20261109.log`, and `Tessella-boundary-cache-1133-optimized2-20261109.log`.
- `Tessella-core-coordinate-112-optimized-20261109.log`, `Tessella-core-coordinate-1133-optimized2-20261109.log`, `Tessella-core-coordinate-{112,1133}-min-functional-20261109.log`, and `Tessella-core-independent-{112,1133}-optimized-20261109.log`.
- `Tessella-distance-closest-regression-112-min-final-20261109.log` and `Tessella-closest-geometry-gate-113-optimized-finalsource-20261109.log`.

## Release and earlier evidence

Version 0.1.0 declares Julia 1.12–1.13 compatibility and the owner-selected MIT
license. Third-party notices remain under `LICENSES/`. AI assistance is disclosed
in the root README. Registration is unfinished: complete the implementation and
final checks, establish the maintainer-review prerequisite, then trigger
Registrator and verify the actual General merge and package installation.

Superseded audit logs, handoffs, repository instructions and design snapshots
were removed from the checkout. Git retains the
[complete earlier record](https://github.com/jake-w-liu/Tessella.jl/tree/fbf34516da30408248e21a592ea599fdb9736b5c/.github/development/history),
including detailed P1–P6 requirements and prior receipts. Their dates and source
revisions matter; they do not establish current behavior or override the user's
active instructions. Current API documentation remains in
[the reference directory](../reference/).
