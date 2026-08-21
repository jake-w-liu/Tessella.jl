# Recombined affine transfinite-hexahedron differential

`differential.jl` checks the bounded `mesh_transfinite_hex` implementation
against the installed Gmsh 4.15.2 Julia API without writing geometry or mesh
files. It builds axis-aligned and sheared affine six-face volumes for every
four-sided transfinite arrangement and compares:

- every node through a tolerance-bounded coordinate bijection;
- exact ordered type-5 connectivity after that node mapping;
- exact ordered type-3 connectivity over the six canonical face entities;
- face and volume physical-tag correspondence; and
- analytic node, hexahedron, and boundary-quadrangle counts.

The implementation was traced against Gmsh tag `gmsh_4_15_2`, commit
`657c8e915f60405e6cad0c8ec7faf812bfff1a60`. Relevant pinned sources are:

- `src/mesh/meshGRegionTransfinite.cpp`, SHA-256
  `3781c0d782eaa1345713fe835d7f430dd2896db532f0bfd90c59f691d962e7cd`:
  canonical corners/faces and `CREATE_HEX` are at lines 30–61; logical node
  interpolation is at lines 394–505; fully recombined selection and `i,j,k`
  element iteration are at lines 520–560.
- `src/mesh/meshGFaceTransfinite.cpp`, SHA-256
  `59e045f19b8118c4522f2056d5357f24319560005fd65e809104e79e70a12ee2`:
  recombined four-sided faces emit `(v1,v2,v3,v4)` at lines 898–907 before
  any triangle-arrangement branch.
- `src/geo/MHexahedron.h`, SHA-256
  `4899375e38ced6992db6c1b1a670073d14def3c84bb60bf5b4dca1a4d38478a2`:
  linear local nodes and the six outward face orders are defined at lines
  28–63 and 191–202.
- `src/geo/MHexahedron.cpp`, SHA-256
  `b35d172b3a382e4a22d17d2646e445cba788efa9c5259596e3c9bfabca948392`:
  Gmsh's local volume-sign convention is implemented at lines 48–67.

Gmsh obtains straight-curve nodes through a parametric inversion path; the
pinned Homebrew build leaves coordinate residuals around `1e-12`. Node matching
therefore uses `65536eps(Float64) * max(coordinate_scale, 1)`. Connectivity
comparisons remain exact after mapping.

Run the persistent gate with either supported Julia runtime:

```sh
julia +1.12 --project=. --startup-file=no --check-bounds=yes \
  validation/transfinite_hex/differential.jl
julia +1.11 --project=. --startup-file=no --check-bounds=yes \
  validation/transfinite_hex/differential.jl
```

The production path additionally certifies every represented trilinear
hexahedron with outward-rounded degree-(2,2,2) Bernstein Jacobian intervals and
an exact-rational fallback over the represented dyadic nodes. It integrates the
certified Jacobian polynomial, encloses every physical cell measure with an
exponent-scaled interval, and uses a constant-memory balanced sum to require
agreement with the exact affine volume to a `65_536`-ulp normalized error bound,
without Float64 underflow or cell-count cliffs. Exact corner orientations,
finite positive six-tetrahedron subcell measures and boundary areas,
constant-memory internal-face pairing, exact boundary-face equality,
deterministic CRCs, pre-allocation count limits, and linear allocation growth
are separate gates.

This bounded increment does not support five-face degeneracies, curved/warped
or independently discretized faces, nonuniform curve laws, partial
recombination into prisms or pyramids, QuadTri, holes, multiple blocks,
periodic seams, embedded entities, or high-order elements.
