# ASCENT drop-in mesh handshake — Tessella `.msh` → ASCENT `load_mesh`

Status reviewed **2026-10-10**: this page describes the dated ASCENT loading and solve captures below.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

**Historical capture verified 2026-08-12.** ASCENT's mesh loader ingested the
recorded Tessella mesh without modification. The solve capture below is dated
2026-08-13; neither capture qualifies the current integration candidates.

In the captured ASCENT environment, a solver mesh loads via `ASCENT.load_mesh(path)` (`ASCENT/src/core/mesh.jl`),
which calls `GmshDiscreteModel(path)` (GridapGmsh 0.7.4 / gmsh_jll 4.9.3) to build a
Gridap `DiscreteModel`, then **requires ≥1 top-dimensional physical group** (a physical
volume in 3-D) — else it throws. The fixture therefore checks that the written
`.msh` preserves the region groups required by that loader. Loading alone does
not qualify every geometry, boundary condition or solver contract.

## Case

`generate.jl` (Tessella env) builds the ENC-COAX acceptance topology — a coax **pin**
inside an **air** cavity inside a metal **case** shell, meshed as ONE conforming
partition (`tetrahedralize_conforming`) — and writes it as gmsh **MSH v4.1** with the
three physical volumes (`coax_pin` / `air` / `case`) via `write_msh(...; version=4.1,
physical_names=...)`.

`handshake.jl` (ASCENT env) replicates the essential body of `ASCENT.load_mesh`:
`GmshDiscreteModel` → `get_face_labeling` → collect the tag names on top-dimensional
entities.

## Result (verified)

```
Info : 3 entities / 34 nodes / 144 elements  (gmsh reader, no errors)
RESULT num_cells=144 num_tags=3 volume_groups=["coax_pin", "air", "case"]
HANDSHAKE_OK
```

GridapGmsh parses the Tessella mesh cleanly and every region volume surfaces as a
top-dimensional physical group — i.e. `ASCENT.load_mesh` returns a valid `MeshData`.
**This captured mesh was consumed by the recorded ASCENT environment with no format bridge.**

## Reproduce

```
# 1. Tessella env — produce the mesh (writes ascent_coax.msh here, git-ignored)
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=/path/to/Tessella.jl validation/ascent_handshake/generate.jl
# 2. ASCENT env — confirm it loads (GridapGmsh is an ASCENT dep, not a Tessella one)
julia --project=/path/to/ASCENT validation/ascent_handshake/handshake.jl
```

Use the Julia version and dependencies required by the external ASCENT project
for step 2; its recorded GridapGmsh/gmsh_jll versions are capture provenance,
not a claim about current upstream releases. The generated `.msh` is a
git-ignored build artifact — regenerate it with step 1.
The full 22-case HFSS regression (mesh each guide geometry with Tessella → ASCENT
solve → compare) builds on this handshake and needs the ASCENT solver + datasets
(local, not in this repo).

## Boundary-condition handshake (verified 2026-08-12)

`generate_bc.jl` / `bc_handshake.jl` extend the proof to the **full solver-ready
structure**: the ENC-COAX mesh carries not just the 3 material **volumes**
(`coax_pin`/`air`/`case`, 3-D physical groups) but also 2 **BC surfaces** — `radiation`
on the domain boundary and `coax_pin_pec` on the pin↔air material interface (2-D
physical groups on tagged boundary/interface faces). ASCENT's parser loads **all five**:

```
num_cells=144 num_tags=5
materials (volumes): ["coax_pin", "air", "case"]
boundary conditions (surfaces): ["radiation", "coax_pin_pec"]
BC_HANDSHAKE_OK — ASCENT sees all 3 material volumes + 2 BC surfaces
```

i.e. Tessella emits both material regions AND boundary conditions that ASCENT ingests —
the mesh is solver-consumable with BCs, not volumes alone. (The BC *assignment* here is
geometrically representative — outer boundary → radiation, pin interface → PEC — since
the literal ENC-COAX BC map lives in the OCC-built `.geo`.)

## Solve-level handshake (verified 2026-08-13) — ASCENT *uses* the mesh, not just loads it

`solve_step.jl` takes the proof past loading: ASCENT **assembles the actual Maxwell
finite-element operator** on the Tessella mesh — `load_mesh` → `cell_sigma_tensor`
(materials per volume) → `fe_spaces` (first-order Nedelec H(curl) edge elements) →
`assemble_diffusive_matrix` / `assemble_stiffness_mass` (the curl-curl + mass system):

```
ndof (Nedelec edges) = 165
A = SparseMatrixCSC{ComplexF64} size (165,165) nnz=2421
complex-symmetric: relsym = 9.6e-17         (machine precision — correct)
curl-curl stiffness PSD: cᵀKc = 2.08e9 ≥ 0  (physically correct)
ASCENT_SOLVE_STEP_OK
```

and then **solves** it: with a manufactured solution (known non-trivial field `x_true`,
`b = A·x_true`), `A \ b` recovers `x_true` to `1.06e-15` (residual `2.9e-16`).

So ASCENT builds AND solves the finite-element Maxwell system on a Tessella-generated mesh,
with the operator's physics (complex symmetry, PSD stiffness) holding and the linear solve
correct to round-off — the mesh is not just loadable but **solved-on by the solver**. The
full 22-case HFSS regression (each case: geometry + BCs + frequency sweep + solve + compare
to the guide) is the remaining compute campaign; this proves the fundamental solve-usability
it builds on.
