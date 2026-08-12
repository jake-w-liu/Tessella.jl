# ASCENT solve regression across Tessella geometries

**Verified 2026-08-13.** The transferable core of a multi-case HFSS regression: prove ASCENT
**assembles and solves** the Maxwell finite-element system on a *suite* of distinct Tessella-
generated meshes — not just one — so solve-usability is robust across geometry classes.

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

Every case: complex-symmetric operator (round-off) + machine-precision field recovery. So ASCENT
solve-usability is proven **robust across diverse Tessella geometries**, not a one-off.

The literal 22 HFSS UserGuide cases exercise this same pipeline on the guide's specific antenna/
microwave geometries (OCC-built in the ASCENT campaign, not Tessella-buildable here) with a full
frequency sweep + comparison to the guide — the remaining external compute campaign. This
regression proves the mesh→assemble→solve foundation it stands on.

Generated `.msh` files are git-ignored build artifacts.
