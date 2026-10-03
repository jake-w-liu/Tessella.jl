# Native mixed cache lifecycle differential

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia --project=. --check-bounds=yes validation/api_mixed_cache/differential.jl
```

The pinned Gmsh 4.15.2 gate checks standard Point/Line/Triangle/Quadrangle/
Tetrahedron/Hexahedron/Prism/Pyramid blocks, complete order-two node ordering,
linear round trips, affine transforms, and independent coincident triangular and
quadrangular surfaces at both orders, including per-entity node ownership.
It currently performs 152 checks. Production caches preserve native blocks and
classification; the oracle is a development-only dependency.

Curved native CAD order elevation/refinement rejects before cache mutation until
curve/surface placement stencils are available. Straight planar extruded faces and
unrelated curved entities have explicit regression coverage.
