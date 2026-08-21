# Five-face transfinite-prism differential

`differential.jl` checks `mesh_transfinite_prism` against the installed Gmsh
4.15.2 Julia API without writing geometry or mesh files. It builds affine
five-face prisms in canonical corner order, forces the legacy collapsed
three-sided surface path (`Mesh.TransfiniteTri = 0`), leaves all faces
unrecombined with the `Left` arrangement, and compares:

- every node after coordinate-based bijective matching;
- canonical type-4 tetrahedron connectivity;
- canonical type-2 connectivity over all five boundary faces; and
- the analytic node, tetrahedron, and boundary-triangle counts.

The production path additionally preserves Gmsh's single canonical
tetrahedron vertex order, rejects zero or reversed exact-predicate signs,
requires exact combinatorial agreement between emitted and extracted boundary
triangles, and compares a compensated, exponent-scaled determinant sum with
the analytic affine-prism volume. Independently scaled edge columns avoid
Float64 exponent underflow/overflow inside this relative audit;
ill-conditioned determinants use an exact dyadic fallback. Final `Mesh`
validation still requires represented triangle areas and tetrahedron volumes
to be finite Float64 values. These gates reject inputs whose intermediate
Float64 layers fold or lose material volume even when the six input corners
are exactly affine, without introducing a cell-count cliff below the Float64
volume range.

The implementation was traced against Gmsh tag `gmsh_4_15_2`, commit
`657c8e915f60405e6cad0c8ec7faf812bfff1a60`. The checked source SHA-256 values
are `3781c0d782eaa1345713fe835d7f430dd2896db532f0bfd90c59f691d962e7cd`
for `src/mesh/meshGRegionTransfinite.cpp` and
`59e045f19b8118c4522f2056d5357f24319560005fd65e809104e79e70a12ee2`
for `src/mesh/meshGFaceTransfinite.cpp`.

Run:

```sh
julia --project=. --startup-file=no --check-bounds=yes \
  validation/transfinite_prism/differential.jl
```

Success prints `TRANSFINITE_PRISM_DIFFERENTIAL_OK` with the runtime version,
aggregate counts, and maximum matched-node error.

This validation does not claim the compact `Mesh.TransfiniteTri = 1` volume
path, curved or independently discretized faces, nonuniform curve laws,
recombination, QuadTri, holes, multiple blocks, periodic seams, or high-order
elements. It also does not claim coordinate scales whose derived cell measures
overflow Float64, even when the coordinates themselves remain finite.
