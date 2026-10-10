# Transfinite QuadTri differential

Status reviewed **2026-10-10**: this page describes the bounded transfinite QuadTri driver and recorded artifacts.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/quadtri_transfinite/differential.jl
```

The pinned Gmsh 4.15.2 gate compares 117 fixtures: all 64 six-face and 32
collapsed-prism recombination masks, arrangements, graded spacing, warped grids,
and compact-prism behavior. It checks node identities, per-entity type counts,
coordinates, and every typed cell vertex set. Compact `Mesh.TransfiniteTri=1`
prisms follow upstream's ordinary branch, which ignores `TransfQuadTri`.

The focused tests additionally check analytic positive volumes, boundary coverage,
face incidence, classification, physical ownership, malformed/folded inputs,
resource limits, allocation growth, and deterministic CRC records. Six records
shared with the AddVerts gate match Julia 1.12.7 and 1.13.1; see
`test/artifacts/quadtri_crc.txt`. Production meshing never calls Gmsh.
