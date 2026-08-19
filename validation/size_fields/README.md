# Gmsh size-field differential

`differential.jl` meshes the same one-dimensional domain with Tessella and the
installed Gmsh Julia API. It exercises five independent field graphs:

- `Distance → Threshold` with a linear 0.1-to-0.5 transition;
- `Box` with `VIn=0.1`, `VOut=0.5`, and a 0.5-thick exterior transition;
- `Ball` with the equivalent spherical exterior transition;
- finite `Cylinder` with strict radial and cap boundaries;
- annular `Frustum` with radial interpolation and an unconstrained exterior.

The gate compares node counts and fine/coarse regional spacing. Exact node
coordinates are not expected: Tessella uses equal metric-length placement while
Gmsh uses its own curve discretizer. The tolerances require both implementations
to reproduce the same field-driven density and transition behavior.

Run:

```sh
julia --project=. validation/size_fields/differential.jl
```

If `gmsh.jl` is not installed next to the `gmsh` executable, set
`GMSH_JULIA_API` to its absolute path.

The field definitions follow the official Gmsh 4.15.2 field reference:
<https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-mesh-size-fields>.
