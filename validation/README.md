# Validation — Tessella vs external mesh tools

Independent cross-validation of Tessella against established open-source meshers
(currently **gmsh**). The convention: for every geometry, mesh the *same domain*
with both tools and compare against a shared oracle — the **analytic volume** — plus
element quality and wall-clock time. This checks that Tessella is not merely
self-consistent but agrees with an independent implementation, and quantifies where
it is more accurate or faster.

## Layout

```
validation/
  common.jl              # helpers: gmsh runner, tet metrics, comparison
  run_all.jl             # driver: runs every case, writes REPORT.md
  REPORT.md              # generated results table (git-ignored until you run it)
  size_fields/
    differential.jl      # required Gmsh 4.15.2 field differential
    STATUS.md            # exact coverage and explicit non-claims
  uniform_refine/
    differential.jl      # required Gmsh 4.15.2 linear-simplex template differential
  transfinite/
    differential.jl      # required Gmsh 4.15.2 four-sided patch differential
  transfinite_curve/
    differential.jl      # required Gmsh 4.15.2 straight-curve-law differential
  transfinite_triangle/
    differential.jl      # required Gmsh 4.15.2 three-sided patch differential
  transfinite_quad/
    differential.jl      # required Gmsh 4.15.2 recombined-quad differential
  transfinite_volume/
    differential.jl      # required Gmsh 4.15.2 affine-volume differential
  transfinite_prism/
    differential.jl      # required Gmsh 4.15.2 five-face-prism differential
  cases/
    01_box/box.geo               # reference gmsh script (retained)
    02_cylinder/cylinder.geo
    03_box_tunnel/box_tunnel.geo
    04_hollow_box/hollow_box.geo
    05_sphere/sphere.geo
    06_enclosure_coax/enclosure_coax_junction.geo   # ASCENT acceptance fixture
```

Each case folder keeps its reference `.geo` script. Generated `gmsh_out.msh` files
are git-ignored.

## Run

```sh
julia --project=. --check-bounds=yes validation/run_all.jl
```

The aggregate gate requires the Gmsh 4.15.2 CLI and matching Julia API. It launches
the size-field, uniform-refinement, four-sided transfinite, straight transfinite
curve-law, three-sided transfinite, recombined-quadrangle, affine
transfinite-volume, and five-face-prism differentials as mandatory bounds-checked
children; missing or wrong-version Gmsh, failed probes, and parity mismatches make
the aggregate command fail. Mesh-case results print to the terminal and to
`validation/REPORT.md`.

## What each case checks

| Case | Geometry | Oracle | Point |
|------|----------|--------|-------|
| 01_box | axis-aligned box | V = 2 exact | both meshers must conform to a flat solid |
| 02_cylinder | solid cylinder | V = πR²H | curved-surface fidelity trade-off (Tessella inscribed N-gon vs gmsh true circle) |
| 03_box_tunnel | box with a through-tunnel | V = 24 exact | genus-1 flat solid |
| 04_hollow_box | box minus interior cavity | V = 35 exact | Boolean-difference (CSG) solid |
| 05_sphere | ball | V = 4/3·πR³ | curved-surface fidelity |
| 06_enclosure_coax | ASCENT coax feed-through | volumes non-empty | the acceptance case: gmsh 4.13/4.15 leave the air/case/pin volumes **empty**; Tessella's native pipeline targets filling them |

Flat solids give an *exact* analytic volume, so a passing row is a hard correctness
cross-check. Curved solids expose the geometric-fidelity trade-off honestly (each
tool is exact for its own surface model).

## Adding a case

Drop a `.geo` in a new `cases/NN_name/` folder and add a matching Tessella surface
builder + analytic volume to the `cases` list in `run_all.jl`.

## Extending to other tools

The convention is tool-agnostic: add a `run_<tool>` helper in `common.jl` (mirroring
`run_gmsh`) that meshes to a Tessella-readable format, and a column in the report.
Retain each tool's input script alongside the case.
