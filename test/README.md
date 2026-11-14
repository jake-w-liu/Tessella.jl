# Testing

CI tests Julia 1.12 and 1.13 on Linux, Windows, and macOS. The default suite checks small real fixtures
with analytic values, topology, ownership, and error-path checks. It needs no
external mesher or solver.

Run ordinary package tests from the repository root:

```powershell
julia +1.12 --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
julia +1.13.3 --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

The runner selects minimal compilation automatically before loading Tessella.
Child failures propagate. Each functional file reports start and completion,
without timing assertions. Standalone invocation runs the same checks:

```powershell
julia --project=. --startup-file=no --check-bounds=yes --compile=min -O0 test/runtests.jl
```

For the full suite with normal optimized compilation, explicitly opt in:

```powershell
julia --project=. --startup-file=no -O2 -e 'ENV["TESSELLA_TEST_COMPILER"]="optimized"; using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "-O2"])'
```

CI collects coverage in a separate Julia 1.13 Linux job with normal compilation.
The compiled suite exercises predicates, topology, coordinate queries, quality,
transforms, 1-D/2-D/3-D meshing, refinement, periodicity, structured meshes, NURBS,
healing, fields and MSH/STL I/O. It verifies executed function-body counters,
attaches an LCOV report and records total and covered source-line counts. The
report records lines exercised by this suite, rather than the separately run
session API suite.
Minimal compilation does not reliably emit function coverage counters.

Run the same coverage suite locally:

```sh
julia --project=. --startup-file=no --compiled-modules=no --check-bounds=yes --code-coverage=user -O2 test/coverage/runtests.jl
python3 test/coverage/verify.py
```

The default covers predicates, fields, small 1D/2D/3D kernels, structured meshes,
refinement, high order, I/O, API state, CAD boundaries, partitions, views and
plugins. It also includes native Boolean seam/fuse-image controls, finite
natural BREP faces and solids, centroidal CAD measures, wedge orientation, and
the native ONELAB process server (including finite deadline/backpressure
cleanup). Torus/box curved Boolean mass integration, full NURBS pillow and
spherical-cap volume recovery, full `.geo`
adaptation and geometry-mesh matching, native-oracle sweeps, allocation stress,
benchmarks, and HFSS/external solver campaigns are opt-in:

```powershell
julia --project=. --startup-file=no --check-bounds=yes test/extended/runtests.jl
```

The extended runner includes those cases and larger regression tests.
It can take hours. Run the individual relevant file when checking a changed
feature; documentation-only changes need a documentation build. External
differential drivers under local `validation/` have their own tool requirements
and are never production dependencies.

See the [contribution guide](../docs/src/contributing.md) for documentation builds
and change-verification guidance.
