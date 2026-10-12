# Testing

The default suite checks functionality and correctness with small real public
calls: geometric predicates, mesh validation and ownership, transforms, entity
topology, NURBS and fields, `.geo` parameter reading, bounded 1D/2D/3D meshing,
structured patches, refinement and recombination, quadratic tetrahedra,
MSH/STL I/O, API state, query ownership and failure atomicity, plus small matching
and compact surface coverage regressions, list/model-view format checks and
independent multi-root CAD boundaries, nested and mixed mapped/product STEP
product-usage forest paths
and planar Face holes, IGES143
boundaries and explicit plane UV/rigid-transform chains, complete quadratic
Hex, Pyramid13 and orders3–10 Curve VTK export, actual two-Tet partition
interfaces, closed Tri3-BRep seam transport and ghost queries, declared
Curve/Point lower-BRep partition topology, a mixed linear family
(Tet4+Prism6 over Tri3+Quad4) BRep partition, a complete second-order
(two-Tet10 over Tri6) BRep partition, a three-Tet4 junction BRep partition
with junction seam and points, a coincident-Surface BRep partition, a
complete third-order (two-Tet20 over Tri10) BRep partition, a periodic
link and parametric-node BRep partition transport, complete
catalogued high-order Surface/Volume VTK, directional constant-tensor
box and shared-material grids, floating-point Pyramid options, CLI help/version
and read-only MSH checking, plus atomic current-model opening and merging,
including one POS point import and an empty-file merge, plus derived two-step
view metadata and display selection, two list/model replacement checks, and a
quarter-turn periodic Surface copy. It needs
no external mesher or solver.

Run the default suite from the repository root. The runner automatically uses
minimal compilation for small functional checks. Ordinary `Pkg.test()` and
standalone invocation relaunch once with the same active project when needed,
before loading Tessella. Child failures propagate. Fixed custom coverage
filenames retain the original invocation; normal coverage flags are inherited.
This keeps large production dispatcher optimization out of routine tests:

```powershell
julia --project=. --startup-file=no --check-bounds=yes --compile=min -O0 test/runtests.jl
```

`Pkg.test()` runs the same default suite in Julia's package test environment:

```powershell
julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

For checks of the normal optimized compiler path, explicitly opt in to the
same assertions with normal compilation:

```powershell
julia --project=. --startup-file=no -O2 -e 'ENV["TESSELLA_TEST_COMPILER"]="optimized"; using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "-O2"])'
```

Local receipts on 2026-10-11, with bounds checks and one Julia thread:

| Runtime | Mode and entry point | Passing assertions | Elapsed |
| --- | --- | ---: | ---: |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 42.3 s | 279 | ~45 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 46.5 s | 279 | ~49 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 44.4 s | 279 | ~47 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 47.2 s | 279 | ~50 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 41.3 s | 272 | ~42 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 45.2 s | 272 | ~46 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 41.0 s | 272 | ~42 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 43.5 s | 272 | ~44 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 40.5 s | 269 | ~41 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 44.0 s | 269 | ~45 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 40.5 s | 269 | ~41 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 43.3 s | 269 | ~44 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 40.7 s | 266 | ~41 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 43.9 s | 266 | ~44 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 40.7 s | 263 | ~41 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 41.4 s | 263 | ~42 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 38.1 s | 259 | ~39 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 41.5 s | 259 | ~42 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 37.5 s | 256 | ~70 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 40.8 s | 256 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 36.3 s | 249 | ~70 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 41.7 s | 249 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 34.2 s | 246 | ~70 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 38.7 s | 246 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 31.9 s | 243 | ~65 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 35.6 s | 243 | ~70 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 24.4 s | 240 | ~60 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 27.2 s | 240 | ~60 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 17.2 s | 237 | ~50 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 16.9 s | 237 | ~40 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 15.7 s | 234 | ~50 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 15.2 s | 234 | ~40 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 15.4 s | 231 | ~50 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 14.8 s | 231 | ~40 s |
| Julia 1.12.7 | Ordinary `Pkg.test()`, automatic minimal compilation | 224 | 45.21 s |
| Julia 1.13.1 | Ordinary `Pkg.test()`, automatic minimal compilation | 224 | 38.38 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 222 | 45.18 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 222 | 38.74 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 220 | 47.05 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 220 | 40.34 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 218 | 45.41 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 218 | 38.66 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 214 | 49.28 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 214 | 38.93 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 209 | 46.49 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 209 | 39.23 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 207 | 45.23 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 207 | 37.67 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 203 | 46.00 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 203 | 39.12 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 201 | 44.88 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 201 | 36.53 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 192 | 47.14 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 192 | 39.75 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 190 | 47.66 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 190 | 37.80 s |
| Julia 1.12.7 | `--compile=min -O0`, preceding `Pkg.test()` | 180 | 38.97 s |
| Julia 1.13.1 | `--compile=min -O0`, preceding `Pkg.test()` | 180 | 42.20 s |
| Julia 1.13.1 | `--compile=min -O0`, earlier CLI `test/runtests.jl`, warm | 173 | 5.71 s |
| Julia 1.13.1 | `--compile=min -O0`, preceding CAD/model/periodic `Pkg.test()` | 171 | 40.61 s |
| Julia 1.13.1 | `--compile=min -O0`, earlier CLI/list increment `Pkg.test()` | 163 | 36.41 s |
| Julia 1.12.7 | `--compile=min -O0`, earlier `Pkg.test()` | 157 | 37.03 s |
| Julia 1.13.1 | Normal `-O2`, earlier `Pkg.test()` | 157 | 217.47 s |
| Julia 1.12.7 | Normal `-O2`, earlier standalone diagnostic | 157 | 301.75 s |

The current 279-assertion bodies took44.4/47.2 s onJulia1.12.7/1.13.1.
Cold elapsed times include about30/25 s of precompilation.
These ordinary invocations verify the automatic runner in each package test
environment on the same final Source. Logs and actual exit-zero receipts are
`test/tmp/ordinary_default_224_{112,113}_20261011.*`. The same252 Source files
remained frozen throughout both calls. The latest increment adds a
periodic/parametric BRep partition testset (3 assertions: an imported
`$Periodic` link and UV-parametric vertices survive partition on their
original entity tags, `validate`, unpartition restore) to the
declared-BRep default smoke plus a 28-assertion periodic set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
364 assertions on both runtimes (a discrete periodic prism fixture — an
imported dim-2 `$Periodic` link with an affine transform, six
UV-parametric boundary vertices, Tet4+Tri3 cells — verifies verbatim link
transport through partition, ghosts, MSH4 round-trips and `unpartition`,
same-dimension same-owner parametric reclassification, and the bare-Tet4
producer path). The preceding increment adds a
third-order BRep partition testset (4 assertions: two Tet20 sharing a Tri10
facet — Tri10 interface cell, three generated Line4 seam cells, `validate`,
unpartition restore) to the declared-BRep default smoke plus a 40-assertion
complete third-order set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
336 assertions on both runtimes (Line4 curves, Tri10/Quad16 surfaces and
Tet20/Hex64/Prism40/Pyr30 volumes partition through the same declared-BRep
certificate — a Tet20 pair rimmed by a declared closed Line4 Curve, a
Hex64 pair with a closed four-cell Line4 seam, a Prism40 pair mixing
Tri10/Quad16 exteriors and a Pyr30 pair sharing the Quad16 base, plus
ghosts, MSH4 round-trips, unpartition and atomic rejection of Tri15,
Quad12 and Tet35 blocks). The preceding increment adds a
coincident-support BRep partition testset (3 assertions: two Regions sharing
a facet carried by one declared Surface — the partition child merges [1,2]
membership, binds both Region children with opposite signs, no duplicate
interface entity or cell, unpartition restore) to the declared-BRep default
smoke plus a 22-assertion coincident set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
296 assertions on both runtimes (the merged child keeps cell and UID,
declared Curve/Point rim roles persist, MSH4 round-trips and the atomic
one-sided rejection are covered). The preceding increment adds a junction
BRep partition testset (3 assertions: three Tet4 meeting at a junction rim
edge, two dim-0 junction entities, one junction seam carrying all three
partitions, unpartition restore) to the declared-BRep default smoke plus a
40-assertion junction set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
274 assertions on both runtimes (multi-volume junction fans partition
through the manifold rim certificate; junction seams share all incident
interfaces under Surface-or-Region parents; Point15 junction entities reuse
declared Point children; a carrier-free four-volume fan, a declared
Curve/Point junction and a Tet10 junction with a Line3 seam are covered
alongside MSH4 round-trips and atomicity). The preceding increment adds a second-order
lower-BRep partition testset (3 assertions: two Tet10 sharing a Tri6 facet,
Tri6 interface cell, three generated Line3 seam cells, unpartition restore) to
the declared-BRep default smoke plus a 46-assertion complete second-order set
in `test/interfaces/api_partition_lower_brep_functional_test.jl`, now
234 assertions on both runtimes (Line3 curves, Tri6/Quad9 surfaces and
Tet10/Hex27/Prism18/Pyr14 volumes partition through the same declared-BRep
certificate, including declared Line3 border edge-node checks, ghosts, MSH4
round-trips, unpartition and atomic violations for foreign mids, serendipity
families and reversed volumes). The preceding increment adds a mixed
linear family BRep partition testset (4 assertions: Tet4+Prism6 sharing a
Tri3 face, Tri3 interface, generated seams, unpartition restore) to the
declared-BRep default smoke plus a 61-assertion mixed-family set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
188 assertions on both runtimes (Quad4 surfaces and Hex8/Prism6/Pyr5
volumes partition through the same declared-BRep certificate, including
ghosts, MSH4 round-trips and atomic violations). The preceding increment
adds three
cheap MSH4 signed-boundary import checks to the declared-BRep default
smoke (10 assertions) plus a 27-assertion import testset in
`test/interfaces/api_partition_lower_brep_functional_test.jl`. The preceding increment adds a declared
Curve/Point lower-BRep partition testset (7 assertions) on that fixture's
original 100 assertions. The preceding increment added three cheap
meridian-certified `TOROIDAL_SURFACE` segment checks;
`test/geometry/brep_step_torus_segment_functional_test.jl` is the opt-in
segment fixture (boundary-only — toroidal fills remain a hard tetrahedralizer
input). The preceding increment adds three cheap
rim/apex-certified `CONICAL_SURFACE` checks; `test/geometry/brep_step_cone_functional_test.jl`
is the opt-in frustum/cone fixture set. The preceding increment adds three cheap
rim-certified `CYLINDRICAL_SURFACE` checks; `test/geometry/brep_step_cylinder_functional_test.jl`
is the opt-in cylinder fixture. The preceding increment adds three cheap
closed-face `SPHERICAL_SURFACE` checks; `test/geometry/brep_step_sphere_functional_test.jl`
is the opt-in sphere fixture. The preceding increment adds three cheap
STEP unit-conversion checks; `test/geometry/brep_step_units_functional_test.jl`
is the opt-in conversion fixture set. The preceding increment added three cheap
mixed STEP forest checks; `test/geometry/brep_step_relationship_mixed_functional_test.jl`
is the opt-in mixed fixture. The preceding increment added four cheap
quarter-turn periodic Surface-copy checks; `test/interfaces/api_rigid_periodic_hybrid_functional_test.jl`
is the opt-in rigid/hybrid fixture. The preceding224-check increment added two cheap
list/model representation checks. The preceding218-check increment added two
file-merge and two read-only CLI assertions; periodic Volume, tensor production, native and solver
campaigns remain opt-in. Portable final224 receipts and Source hashes are also
committed under `test/artifacts/handoff_20261011`; see the
[Devin continuation checklist](../DEVIN_HANDOFF.md). The preceding214-check bodies took11.1/10.9 s and
qualify their own Source. The preceding 209-check
bodies took 11.1/10.7 s and qualify their own Source. The preceding 207-check
bodies took 10.1 s on both runtimes and qualify their own Source. The preceding203-check
bodies took10.4 s on both; their receipts remain separate. The preceding 201-check
bodies took 9.7/9.1 s; their receipts remain separate. The earlier 192-check
bodies took 9.3/9.7 s and remain recorded separately.
The earlier190-check
bodies took10.3/8.8 s. The earlier 180-assertion
bodies took 7.2/6.2 s; those explicit minimal-compilation receipts are historical.
The earlier 173-assertion standalone run took 5.71 s overall, with a 5.3 s
body. The preceding 171-assertion package run also had a 5.3 s body. The earlier 163-assertion check took 1.4 s.
The earlier 157-assertion routine check took
1.6 s. The earlier optimized
package receipt includes 24.3 s of precompilation and 188.4 s of tests. The
earlier standalone diagnostic used the same assertions with extra API progress
messages and included 36.14 s of precompilation. Those optimized receipts
precede the final narrow native API corrections, which were checked separately
with focused public regressions. API generation dominated normal compilation
time. Both commands exercise the same public operations and assertions; these
timings include their compilation costs.

Keep default fixtures small. Check observable behavior against analytic values,
valid topology, ownership and explicit errors. Exhaustive combinations, large
native meshes, resource/benchmark campaigns, pinned parity and source audits,
and application/solver campaigns belong in opt-in validation.

The complete historical regression files remain in their existing directories.
An extended runner combines the original and newer regression lists:

```powershell
julia --project=. --startup-file=no --check-bounds=yes test/extended/runtests.jl
```

This extended suite contains exhaustive predicate/configuration sweeps,
allocation and source audits, native meshing stress cases and all 22 HFSS
geometries. It can take hours. Run a specific existing regression file when
investigating its feature instead of automatically running every campaign.
Full `.geo` execution is opt-in too: even the two-point For/If fixture in
`extended/geo_functionality.jl` compiles the full native dispatcher. Default
geometry tests exercise entity operations and the parameter reader directly.
`extended/legacy_runtests.jl.txt` and `extended/devin_runtests.jl.txt` preserve
the previous runners verbatim for historical reference. The extended runner
omits repository-layout enforcement and the static development-stage assertion.

External Gmsh comparison campaigns are separately opt-in:

```powershell
julia --project=. --check-bounds=yes validation/run_all.jl
```

Consult the relevant validation directory for its dependencies and resource
requirements. Archived artifacts and validation evidence are retained.
