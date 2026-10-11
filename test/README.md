# Testing

The default suite checks functionality and correctness with small real public
calls: geometric predicates, mesh validation and ownership, transforms, entity
topology, NURBS and fields, `.geo` parameter reading, bounded 1D/2D/3D meshing,
structured patches, refinement and recombination, quadratic tetrahedra,
MSH/STL I/O, API state, query ownership and failure atomicity, plus small matching
and compact surface coverage regressions, list/model-view format checks and
independent multi-root CAD boundaries and planar Face holes, complete quadratic
Hex VTK export, plus CLI help/version. It needs
no external mesher or solver.

Run the default suite from the repository root. For routine functional checks,
minimal compilation avoids compiling large native planning paths that these
small fixtures do not use:

```powershell
julia --project=. --startup-file=no --check-bounds=yes --compile=min -O0 test/runtests.jl
```

`Pkg.test()` runs the same default suite in Julia's package test environment:

```powershell
julia --project=. --startup-file=no --compile=min -O0 -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "--compile=min", "-O0"])'
```

For checks of the normal optimized compiler path, run the same assertions with
normal compilation:

```powershell
julia --project=. --startup-file=no -O2 -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "-O2"])'
```

Local receipts on 2026-10-11, with bounds checks and one Julia thread:

| Runtime | Mode and entry point | Passing assertions | Elapsed |
| --- | --- | ---: | ---: |
| Julia 1.12.7 | `--compile=min -O0`, current `Pkg.test()` | 180 | 38.97 s |
| Julia 1.13.1 | `--compile=min -O0`, current `Pkg.test()` | 180 | 42.20 s |
| Julia 1.13.1 | `--compile=min -O0`, earlier CLI `test/runtests.jl`, warm | 173 | 5.71 s |
| Julia 1.13.1 | `--compile=min -O0`, preceding CAD/model/periodic `Pkg.test()` | 171 | 40.61 s |
| Julia 1.13.1 | `--compile=min -O0`, earlier CLI/list increment `Pkg.test()` | 163 | 36.41 s |
| Julia 1.12.7 | `--compile=min -O0`, final `Pkg.test()` | 157 | 37.03 s |
| Julia 1.13.1 | Normal `-O2`, earlier `Pkg.test()` | 157 | 217.47 s |
| Julia 1.12.7 | Normal `-O2`, earlier standalone diagnostic | 157 | 301.75 s |

The current 180-assertion package bodies took 7.2 s on Julia 1.12.7 and 6.2 s
on Julia 1.13.1. Cold elapsed times include about 30/34 s of precompilation.
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
