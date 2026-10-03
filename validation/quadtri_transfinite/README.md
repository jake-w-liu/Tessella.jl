# Transfinite QuadTri differential

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia --project=. --check-bounds=yes validation/quadtri_transfinite/differential.jl
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
