# Native mixed cache lifecycle differential

Status reviewed **2026-10-10**: this page describes the bounded mixed-cache differential.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/api_mixed_cache/differential.jl
```

The pinned Gmsh 4.15.2 gate checks standard Point/Line/Triangle/Quadrangle/
Tetrahedron/Hexahedron/Prism/Pyramid blocks, complete order-two node ordering,
linear round trips, affine transforms, and independent coincident triangular and
quadrangular surfaces at both orders, including per-entity node ownership.
Its recorded bounded scope performs 152 checks. Production caches preserve native blocks and
classification; the oracle is a development-only dependency.

The published-baseline scope rejects unsupported curved native CAD elevation
or refinement before cache mutation; straight planar extruded faces and unrelated
curved entities have explicit regression coverage. New curved P2 placement and
raw-refinement context repairs have separate Source/native evidence and pending
whole public gates in the handoff; this 152-check driver does not qualify them.
