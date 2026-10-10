# QuadTriAddVerts extrusion differential

Status reviewed **2026-10-10**: this page describes the bounded AddVerts driver and recorded geometry certificates.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

Run the native extrusion cases against the pinned Gmsh 4.15.2 Julia API:

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/geo_quadtri_extrude/differential.jl
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
and output-sized allocation growth. The separate `QuadTriNoNewVerts` planner
is covered by `validation/quadtri_nonew/`; its diagonal-selection algorithm
uses a certified cap chain for the isolated source-quadrangle slice.

These bounded geometry checks do not establish full public operator-history or
source-provenance parity for all sweep categories. The handoff retains separate
AddVerts and actual-cap NoNew qualifications, plus general hybrid generation,
shared/global identity and sizing work. Their pending gates remain required.
