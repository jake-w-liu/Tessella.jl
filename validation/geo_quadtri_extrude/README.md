# QuadTriAddVerts extrusion differential

Run the native extrusion cases against the pinned Gmsh 4.15.2 Julia API:

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia --project=. --check-bounds=yes validation/geo_quadtri_extrude/differential.jl
```

The gate covers quadrilateral and triangular sources, untouched interior
hexahedra and prisms, transition tetrahedra and pyramids, free and recombined
laterals, uniform and grouped layers, translation, rotation, and twist,
fixed rotation columns, and closed revolutions with mixed ordinary/QuadTri sweeps.
It compares node counts, coordinates, every surface/volume element-type count,
and every cell's vertex set with its element type. It preserves unused final-layer
centroids in toroidal regions, matching the upstream vertex preallocation pass.
The native output must pass validation. The production mesher never calls Gmsh.

The geometry unit tests additionally certify exact surface/volume face
conformity, shared laterals against neighboring QuadTri, hexahedral, and
tetrahedral sweeps, standalone classified projection, folded-fan rejection,
and output-sized allocation growth. `QuadTriNoNewVerts` remains an explicit
blocker; its distinct diagonal-generation algorithm is not replaced by a fan.
