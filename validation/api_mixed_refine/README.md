# Native mixed API refinement differential

Status reviewed **2026-10-10**: this page describes the bounded mixed-refinement differential.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

Run on both supported runtimes with the pinned Gmsh 4.15.2 Julia API:

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/api_mixed_refine/differential.jl
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

Curved raw-refinement owner/parameter preservation has a separate corrected
Source candidate and native causal controls. Its complete public qualification
is pending; the template comparisons above do not establish that result.
