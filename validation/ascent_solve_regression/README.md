# ASCENT solve regression across Tessella geometries

Status reviewed **2026-10-10**: this page describes the dated four-case ASCENT solve capture below.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

**Historical capture verified 2026-08-13.** ASCENT assembled and solved the
Maxwell finite-element system on the four distinct Tessella meshes listed below.
These captured results have not been rerun on the current integration candidates.

`generate_cases.jl` (Tessella env) writes 4 meshes covering single/multi-region, box, cylinder,
and nested geometries:

| case | regions | geometry |
|---|---|---|
| `box_cavity` | 1 | unit box |
| `nested_box` | 2 | dielectric core in a shell (native box CSG) |
| `coax3`      | 3 | coax pin / air / case (cylinder-in-box-in-shell) |
| `cyl_cavity` | 1 | faceted cylinder |

`solve_all.jl` (ASCENT env) loads each, assembles the Maxwell operator (`assemble_diffusive_matrix`,
Nedelec H(curl)), and solves it via a manufactured solution (`b = A·x_true`, `x = A\b`, check
`x ≈ x_true`). Result:

```
box_cavity   regions=1 ndof=  ..  relsym~1e-17  solve_err~1e-16  OK
nested_box   regions=2 ndof=  ..  relsym~1e-17  solve_err~1e-15  OK
coax3        regions=3 ndof= 155  relsym=5.8e-17 solve_err=1.9e-15 OK
cyl_cavity   regions=1 ndof=  13  relsym=8.5e-17 solve_err=2.1e-16 OK
SOLVE_REGRESSION_OK — ASCENT assembles + solves Maxwell on every Tessella case
```

Each of the four captured cases had a complex-symmetric operator to round-off
and machine-precision manufactured-field recovery. This is evidence for those
meshes and that ASCENT environment; it does not establish the remaining 22-case
external campaign or the current combined release.

The literal 22 HFSS UserGuide cases exercise this same pipeline on the guide's specific antenna/
microwave geometries (OCC-built in the ASCENT campaign, not Tessella-buildable here) with a full
frequency sweep + comparison to the guide — the remaining external compute campaign. This
regression proves the mesh→assemble→solve foundation it stands on.

## Reproduce

Run from the Tessella repository root, replacing the project paths:

```sh
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=/path/to/Tessella.jl validation/ascent_solve_regression/generate_cases.jl
julia --project=/path/to/ASCENT validation/ascent_solve_regression/solve_all.jl
```

Repeat the Tessella generation step on Julia 1.13. Use the version and dependencies
required by the external ASCENT project for the solve step. Generated `.msh`
files are git-ignored build artifacts.
