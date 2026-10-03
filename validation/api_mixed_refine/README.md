# Native mixed API refinement differential

Run with Julia 1.12.7 and the pinned Gmsh 4.15.2 Julia API:

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia --project=. --check-bounds=yes validation/api_mixed_refine/differential.jl
```

The gate compares every node and typed cell vertex set for linear and quadratic
inputs, repeat refinement, warped cells, a shared hex/prism quadrangle, a shared
pyramid/tetrahedron triangle, reversed volume winding, and a coincident unused
primary node that is pruned during order-two support creation. Ordinary
refinement yields linear children: two segments, four triangles/quadrangles,
eight tetrahedra/hexahedra/prisms, or four pyramids plus eight tetrahedra.
Off-midpoint quadratic quad/tet records confirm Gmsh resets existing support
nodes to first order before recreating them on discrete entities.
Every volume child is normalized to positive winding, matching `RefineMesh`.
Production uses the native templates and never loads Gmsh.
