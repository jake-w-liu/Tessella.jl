# Testing

Requires Julia 1.12.x or 1.13.x. The default suite checks small real fixtures
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

Fixed custom coverage filenames keep their original invocation; normal coverage
flags are inherited. For normal optimized compiler coverage, explicitly opt in:

```powershell
julia --project=. --startup-file=no -O2 -e 'ENV["TESSELLA_TEST_COMPILER"]="optimized"; using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "-O2"])'
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

The extended runner retains those cases and larger historical regressions.
It can take hours. Run the individual relevant file when checking a changed
feature; don't run unrelated campaigns for documentation updates. External
differential drivers under local `validation/` have their own tool requirements
and are never production dependencies.

[Current development evidence](../.github/development/README.md) records actual
scopes and exits. [Historical testing receipts](https://github.com/jake-w-liu/Tessella.jl/blob/fbf34516da30408248e21a592ea599fdb9736b5c/.github/development/history/TESTING.md)
qualify their recorded revisions, not current source. A bounded suite is not
proof of complete mesher compatibility or global optimization.
