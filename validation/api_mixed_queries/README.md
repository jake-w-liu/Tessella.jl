# Native mixed mesh query differential

Status reviewed **2026-10-10**: this page describes the bounded mixed-query differential and recorded allocation measurement.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

Run with the pinned Gmsh 4.15.2 Julia API:

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/api_mixed_queries/differential.jl
```

The recorded gate scope covers 60 fixtures and 2,612 checks over standard linear/full quadratic
families: affine and warped maps, genuinely curved P2 volume records, permuted
node tags, cache and record-only queries, Jacobians, inversion/location,
orientations and basis keys, full/primary edge and face nodes, and both barycenter
modes. The loaded library must report version 4.15.2.

Upstream exceptions remain visible: standalone Point octree omissions, a warped
quadrangle Newton update that accumulates an unused third coordinate, and
inconsistent strict quadratic-pyramid octree results. Tessella verifies analytic
membership and keeps strict lookup enabled with the existing reference tolerance.
Curved quadratic location uses conservative candidate bounds. Geometry orders
beyond two and unimplemented native quality names fail explicitly.

Production queries use native element maps and never load Gmsh. Bulk P2
Jacobians reuse basis derivatives; a recorded 100-hexahedron/one-point query
reduced normalized allocations from 734,174 to 89,262 bytes. Single-point curved
inversion remains a separate optimization target. This historical allocation
comparison has not been rerun on Root V24 or Release V6.
