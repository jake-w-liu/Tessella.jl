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
| `09bf184` combined checkpoint | Default `Pkg.test()` on detached, unchanged checkouts passed 2582/2582 on Julia 1.12.7 (29m27.8s) and 2582/2582 on Julia 1.13.3 (31m48.1s). Both use minimal compilation and bounds checks, actual exits 0. Later increments require their own final combined package run. |
| `acd32c0` combined checkpoint | Default `Pkg.test()` on detached, unchanged checkouts passed 3237/3237 on Julia 1.12.7 (26m58.2s) and 3237/3237 on Julia 1.13.3 (29m17.8s). Both use minimal compilation and bounds checks, actual exits 0. Later incidence, pyramid-binding and ongoing increments retain separate checks until the next combined run. |
| `11aee42` finite-normal inverse | Small intrinsic parameters beside ordinary coordinates: Julia 1.12.7 minimal, 162/162; Julia 1.13.3 normal `-O2`, 352/352 including existing warmed zero-allocation checks. Bounds enabled, exits 0. |
| `fc0f599` GEO/view publication | Julia 1.12.7 and 1.13.3 minimal each passed the 256-check view scope and 162-check import/merge/CLI scope. Normal Julia 1.13.3 passed the 52 new regressions. Bounds enabled, exits 0. The first normal plugin block took 11m13s; no compiler improvement is claimed. |
| `8a1f177` extreme nonlinear inverse | Subnormal derivatives and mixed-scale residuals: Julia 1.12.7 and 1.13.3 minimal, 168/168 each; Julia 1.13.3 normal `-O2`, 358/358 including existing ordinary warmed zero-allocation checks. Bounds enabled, exits 0. Exact recovery can allocate. |
| `fbf3451` factory Box sizing | Five production assertions passed on both Julia lines with minimal compilation. The unit box retained its boundary facets and volume, produced 1350 valid tets, and had zero oversized edges among 1818 (largest normalized length 0.98943). The same fixture on clean `e18b989` failed two assertions, including normalized length 1.63031. |
| Simplex-coordinate increment | Julia 1.12.7 and 1.13.3 normal `-O2` each passed 1026 coordinate/reference/locator assertions; minimal functional checks passed 80/80 each. A separate validation, quality, topology and transform audit passed 3124/3124 on each runtime. Bounds enabled, actual exits 0. |
| Triangle distance and automatic sizing | Julia 1.12.7 minimal passed 71/71 regressions. Julia 1.13.3 normal `-O2` passed 39,643/39,643 regressions, neighboring distance/automatic-constructor checks and resource checks. Bounds enabled, actual exits 0. |
| Native binary STL | Julia 1.13.3 minimal passed 58/58 format/frontend checks; Julia 1.12.7 minimal passed 113/113 including neighboring frontend checks. Julia 1.13.3 normal `-O2` passed 63/63, including five zero-allocation normal checks. Bounds enabled, actual exits 0. |
| STEP/IGES file frontends | 74/74 CAD lifecycle checks passed on Julia 1.12.7 minimal and Julia 1.13.3 normal `-O2`. The initial 70-check scope also passed on 1.13.3 minimal; a 1.12.7 minimal runner passed those 70 plus 143 existing import/merge/native-CAD checks. Bounds enabled, actual exits 0. |
| OBJ relative indices | Julia 1.12.7 and 1.13.3 minimal each passed 120/120 OBJ and neighboring surface-format checks, bounds enabled, actual exits 0. Negative references now select the correct vertices at the face's file position; invalid zero/out-of-range references reject before publication. |
| Finite scalar transitions | Julia 1.12.7 minimal passed 72/72 and Julia 1.13.3 normal `-O2` passed 76/76, including four zero-allocation ordinary queries. Bounds enabled, actual exits 0. Exact recovery handles overflowed interpolation spans and Ball/Box distances, including subnormal transitions; sigmoid semantics are preserved. |
| OFF/PLY2 surface orders | Julia 1.12.7 and 1.13.3 minimal each passed 510/510 surface-order and neighboring format checks, bounds enabled, actual exits 0. Catalog-based corner projection covers complete/incomplete orders 1–10, excludes volume cells, and rejects corrupt sources before replacing the destination. |
| Closed MeshAdapt and recovery | Both Julia lines passed 129 closed/Box/quotient checks and 32 final edited-source/P2 checks. Nearby checks passed 74 MeshAdapt assertions on 1.12.7 and 87 algorithm-7 assertions on 1.13.3, with further current first-order/P2 rechecks. All minimal compilation with bounds, actual exits 0. Direct default Model factories passed 9/9 on 1.13.3. |
| Finite blends and anisotropic Restrict | Julia 1.12.7 minimal passed 705/705 and Julia 1.13.3 normal `-O2` passed 711/711, including 576 exact-dyadic oracle cases and six normal zero-allocation checks. Bounds enabled, actual exits 0. Finite cancellation is recovered, and selected entity restrictions retain the tensor metric and original scalar operator. |
| Surface reader records/resources | Julia 1.12.7 and 1.13.3 minimal each passed 704/704 reader/neighbor checks; focused normal `-O2` checks passed 60/60 each. Bounds enabled, actual exits 0. OFF accepts separated headers/comments and validates face counts/indices; tiny truncated OFF/PLY2/PLY/VTK files reject before reserving their declared coordinate buffers. |
| File dispatch inference boundary | Julia 1.12.7 minimal passed 259/259 CAD, GEO-view, import and merge checks, bounds enabled, actual exit 0. The same two dispatch boundaries completed the public STEP-open/mesh/query path on Julia 1.13.3 normal `-O2`, actual exit 0. |
| BDF large fields/connectivity | Julia 1.12.7 and 1.13.3 minimal each passed 220/220 BDF/SU2 checks; the new 43 checks also passed under normal `-O2` on both runtimes. Bounds enabled, actual exits 0. Standard sixteen-column fields retain coordinates and large identifiers, and invalid quadratic connectivity cannot silently become a linear element. |
| Extreme metric spectra | Julia 1.12.7 minimal passed 745/745 and Julia 1.13.3 normal `-O2` passed 753/753 scalar, Restrict and spectral checks, including existing anisotropy-cap neighbors and eight sampled ordinary zero-allocation loops. Bounds enabled, actual exits 0. Physical sizes and bounded tensors survive eigenvalues outside Float64's range; uniform scaling retains exact finite tensor entries. |
| Other native reader resources | Julia 1.12.7 and 1.13.3 minimal each passed 538/538 reader/neighbor checks; focused normal `-O2` checks passed 79/79 each. Bounds enabled, actual exits 0. P3D, MAIL, IR3, UNV and DIFF bound impossible declared payloads before reserving storage and check header arithmetic before overflow. The eight tiny-file probes fell from 0.8–6.6MB to 1.98–4.24KB normal allocation, with permanent 100KB caps. |
| PLY count validation | Julia 1.12.7 and 1.13.3 minimal each passed 199/199 reader/neighbor checks; focused normal `-O2` checks passed 22/22 each. Bounds enabled, actual exits 0. Explicit negative counts reject in ASCII and both binary byte orders, while missing/zero face counts remain valid and failed open/merge preserve the current model. |
| OCC explicit format dispatch | Julia 1.12.7 and 1.13.3 minimal each passed 121/121 import/CAD checks, bounds enabled, actual exits 0. Canonical `step`/`iges` selectors force the corresponding native reader despite the filename suffix; invalid selectors and incompatible payloads preserve held geometry. The genuine native BREP implementation remains active work. |
| Native incomplete pyramid views | Julia 1.12.7 minimal passed 854/854 pyramid/writer checks; Julia 1.13.3 normal `-O2` passed 858/858 with four threads. A separate eight-check normal gate verified sparse construction, zero-allocation cache hits, concurrent cache creation, copied public schemes and sampled near-apex gradients. Bounds enabled, actual exits 0. The subsequent P9 query failure is covered by the coordinate-binding increment below. |
| Cylinder/Frustum geometry | Julia 1.12.7 minimal passed 812/812 and Julia 1.13.3 normal `-O2` passed 822/822 scalar, metric and axial checks, including ten sampled ordinary zero-allocation loops. Bounds enabled, actual exits 0. Predicates retain original finite axes/endpoints, strict Cylinder and inclusive Frustum boundaries, subnormal radial/axial information and representable geometry whose intermediate length overflows. |
| X3D surface/export frontends | Julia 1.12.7 and 1.13.3 minimal each passed 750/750 format/frontend/neighbor checks; focused normal `-O2` checks passed 143/143 each. Bounds enabled, actual exits 0. Every supported triangle/quadrangle order retains its full primary surface; source and both target paths are checked before replacement. Canonical Print.X3d selection works through API, GEO Save and CLI. |
| OCC lower incidence results | Julia 1.13.3 minimal passed 104/104 incidence, selector and existing IGES checks; Julia 1.12.7 minimal passed 39/39 focused incidence/selector checks. Bounds enabled, actual exits 0. `highestDimOnly=false` includes each imported root's actual lower topology with repeated shared incidences, ordered by dimension and excluding unrelated existing entities. |
| Native pyramid coordinate binding | Julia 1.13.3 normal `-O2` passed 203/203 native/reopened query, gradient and field checks plus a 547/547 neighbor/resource scope. Julia 1.12.7 minimal passed 597/597 public/neighbor/scaled checks. Bounds enabled, actual exits 0. Twelve sampled stored-field queries allocate zero warmed bytes. Exact dyadic accumulation binds native pyramid coordinate offsets once, fixing P9 point-location cancellation while preserving the strict shared inverse certificates. |
| Native pyramid collapse denominator | Julia 1.12.7 minimal passed 691/691 public/neighbor checks; Julia 1.13.3 normal `-O2` passed 836/836 component/neighbor/resource checks, with twelve sampled stored queries allocating zero warmed bytes. Bounds enabled, actual exits 0. Native geometry and values retain their actual collapse denominator through the nearest Float64 reference to the apex; independent flags preserve the existing explicit POS matrix convention. The new 86-check default regression covers scalar/vector/tensor values, gradients, both model-data kinds, exact apex behavior and explicit two/four-matrix controls. |

Normal allocation on the bounded truncated-file probes fell from about 2.4MB
to 2.5–4.7KB. The resource checks remain below 100KB for declared counts through
`typemax(Int)` and binary PLY; no giant allocation probe was used.

Incomplete Pyramid P3–P9 model queries and PostView fields now use an exact
symbolic expansion of the existing native edge-trace basis. The seven reference
schemes are cached under a lock; sparse construction reduced the sampled P9
build from 169.4MB to 25.9MB total allocation and warmed lookups allocate zero
bytes. The new checks retain curved geometry and gradients, all-order nodal
values, and ASCII/binary data sections. This is Tessella's existing native
incomplete-pyramid extension; it does not establish compatibility with Gmsh
readers that reject these catalog types. Complete-pyramid and P2 behavior is
covered by the neighboring checks.

The subsequent P9 failure was in native point location before export as well
as after binary reopen, with identical coordinates. Binding its nodal-offset
geometry polynomial with exact dyadic accumulation fixes the cancellation.
The default regression now includes native P9 plus value/gradient round trips
through `API.open` in both data kinds and byte formats. The shared Newton
solver and acceptance certificates are unchanged. Native probes now also retain
the actual collapse denominator near the apex. Separate geometry/value producer
flags preserve explicit POS matrix behavior. Further strict inverse false-miss
cases remain active investigation.

Algorithm 1 now adapts actual sphere, cone and torus charts, and Box volumes
consume their actual adapted faces. Auto selects the native producer on closed
curved charts. Periodic/pole preimages have explicit quotient topology guards;
inherited boundary sizes come from retained cells rather than artificial seed
spokes. Recovery transports actual node parameters and complete UID row
receipts, including retained P2 and legally edited periodic sources, while the
existing ownership audits remain enforced. Default direct Box1 produced 4197
valid positive tets; Sphere1 and AutoSphere produced 204 each. That direct
minimal run took 19m29.8s, so it establishes correctness without a meshing
performance claim.

`API.open` and `API.merge` now dispatch `.step`, `.stp`, `.igs` and `.iges`
directly to the existing native CAD importer and synchronize actual geometry.
Staging preserves held models/views, existing mesh UIDs, fields and forced
bounding context, with rollback on malformed input. The checks also cover
subsequent meshing/transforms, suffix aliases and Unicode filenames. The first
normal block took 9m23.8s cold; no CAD compiler improvement is claimed.

A later measured inference boundary keeps the GEO execution graph out of
unrelated file dispatch. In fresh Julia 1.13.3 normal `-O2` processes, the
sampled public STEP-open call fell from 423.25s/96.87GB total allocation to
253.58s/68.77GB. The variant also completed native curve meshing and node
queries. Those allocation totals include compilation and are not resident
memory or steady-state mesh costs. Cold compilation remains substantial.

STL output now honors `Mesh.Binary` through API, GEO Save and CLI paths.
Native binary output uses explicit little-endian records, Float32 geometry,
and checked facet counts. Invalid input, overflow and geometry collapsed by
Float32 rounding preserve the destination. ASCII remains the default, with
scaled/exact normal recovery for extreme finite coordinates. The normal
frontend block took 9m12.8s cold; this establishes format correctness and the
sampled helper allocation bounds, without a frontend compiler improvement claim.

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
- `Tessella-default-09bf184-{112,1133}-20261109.log`.
- `Tessella-default-acd32c0-{112,1133}-20261109.log`.
- `Tessella-field-audit-gate-optimized-20261109.log` and `Tessella-field-focused-gate-112-min-20261109.log`.
- `Tessella-view-scoped-gates-{1127,1133}-min-20261109.log`, `Tessella-view-import-cli-{1127,1133}-min-20261109.log`, and `Tessella-view-final-gates-1133-normal-20261109.log`.
- `Tessella-field-audit-gate-optimized-overflow-20261109.log`, `Tessella-field-focused-gate-112-min-overflow-final-20261109.log`, and `Tessella-field-focused-gate-113-min-overflow-20261109.log`.
- `Tessella-box-production-sizing-112-final-20261109.log`, `Tessella-box-production-sizing-1133-20261109.log`, `Tessella-boundary-cache-112-min-20261109.log`, and `Tessella-boundary-cache-1133-optimized2-20261109.log`.
- `Tessella-core-coordinate-112-optimized-20261109.log`, `Tessella-core-coordinate-1133-optimized2-20261109.log`, `Tessella-core-coordinate-{112,1133}-min-functional-20261109.log`, and `Tessella-core-independent-{112,1133}-optimized-20261109.log`.
- `Tessella-distance-closest-regression-112-min-final-20261109.log` and `Tessella-closest-geometry-gate-113-optimized-finalsource-20261109.log`.
- `Tessella-io-stl-binary-1133-min-20261109.log`, `Tessella-io-stl-gate-112-min-20261109.log`, and `Tessella-io-stl-binary-1133-optimized-20261109.log`.
- `Tessella-cad-file-frontends-112-min-final-20261109.log`, `Tessella-cad-file-frontends-1133-optimized-20261109.log`, `Tessella-cad-file-frontends-1133-min-final-20261109.log`, and `Tessella-cad-file-frontends-112-min-neighbors-20261109.log`.
- `Tessella-io-obj-gate-{112,1133}-min-20261109.log`.
- `Tessella-scalar-interpolation-gate-112-min-final-20261109.log` and `Tessella-scalar-interpolation-gate-113-optimized-final-20261109.log`.
- `Tessella-io-surface-order-gate-112-min-final-20261109.log` and `Tessella-io-surface-order-gate-1133-min-public-20261109.log`.
- `Tessella-meshadapt-source-gate-1127-20261109.log`, `Tessella-meshadapt-metric-source-gate-1133-20261109.log`, `Tessella-meshadapt-receipt-final-gate-{1127,1133}-20261109.log`, and `Tessella-meshadapt-direct-factory-gate-1133-20261109.log`.
- `Tessella-scalar-metric-gate-{112-min,113-optimized}-20261109.log`.
- `Tessella-io-surface-records-gate-{112,1133}-{min,normal}-20261109.log`.
- `Tessella-file-dispatch-gate-112-min-20261109.log` and `Tessella-cad-dispatch-barrier-20261109.log`.
- `Tessella-io-bdf-large-gate-{112,1133}-{min,normal}-20261109.log`.
- `Tessella-metric-spectral-gate-{112-min,113-optimized}-final-20261109.log`.
- `Tessella-io-foreign-resource-gate-{112,1133}-{min-final,normal}-20261109.log`.
- `Tessella-io-ply-header-gate-{112,1133}-{min-final,normal}-20261109.log`.
- `Tessella-occ-format-regression-{1127,1133}-20261109.log`.
- `Tessella-pyramid-final-gate-{112-min,1133-normal}-20261109.log` and `Tessella-pyramid-cache-gate-1133-normal-20261109.log`.
- `Tessella-axial-field-gate-{112-min,113-optimized}-20261109.log`.
- `Tessella-io-x3d-order-gate-{112,1133}-{min-final,normal}-20261109.log`.
- `Tessella-occ-incidence-gate-{1127,1133}-20261109.log`.
- `Tessella-pyramid-public-gate-1133-O2-20261109.log` and `Tessella-pyramid-neighbor-gate-{112-min,1133-O2}-20261109.log`.
- `Tessella-native-pyramid-apex-final-gate-112-min-20261109.log` and `Tessella-native-pyramid-apex-component-gate-1133-O2-20261109.log`.

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
