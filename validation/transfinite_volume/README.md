# Affine transfinite-volume differential

This gate checks the bounded `mesh_transfinite_volume` implementation against
the installed Gmsh 4.15.2 API, entirely in memory:

```sh
julia --project=. --startup-file=no --check-bounds=yes \
  validation/transfinite_volume/differential.jl
```

The implementation was traced against Gmsh 4.15.2
`src/mesh/meshGRegionTransfinite.cpp` (SHA-256
`3781c0d782eaa1345713fe835d7f430dd2896db532f0bfd90c59f691d962e7cd`):

- lines 30–41 define the canonical `s0`…`s7` corners and `f0`…`f5` faces;
- lines 43–55 state the 5/6-face and transfinite-boundary requirements;
- lines 128–146 define the volume transfinite interpolation;
- lines 394–420 obtain logical counts and chord-length parameters from the
  oriented boundary faces;
- lines 599–610 select `CREATE_SIM_1`…`CREATE_SIM_6` for six unrecombined
  faces, whose connectivity is defined at lines 80–102.

The differential builds two affine six-face volumes, including a sheared and
rotated case, then compares every mapped node, every tetrahedron, and every
boundary triangle arrangement. Gmsh's straight-curve parametric inversion
leaves roughly `1e-12` coordinate residuals on this build, so node matching uses
`65536eps(Float64) * max(coordinate_scale, 1)`; connectivity comparisons are
exact after node mapping.

This gate does not cover five-face degeneracies, curved or independently
discretized faces, nonuniform curve laws, recombined hexahedra/prisms, QuadTri,
holes, multiple blocks, periodic seams, or high-order elements. Those remain
explicitly unsupported by this bounded increment.
