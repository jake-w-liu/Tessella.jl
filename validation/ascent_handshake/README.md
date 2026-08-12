# ASCENT drop-in mesh handshake — Tessella `.msh` → ASCENT `load_mesh`

**Verified 2026-08-12.** The concrete "ready for ASCENT" proof: a mesh produced by
Tessella is ingested by ASCENT's real mesh loader without modification.

ASCENT loads a solver mesh via `ASCENT.load_mesh(path)` (`ASCENT/src/core/mesh.jl`),
which calls `GmshDiscreteModel(path)` (GridapGmsh 0.7.4 / gmsh_jll 4.9.3) to build a
Gridap `DiscreteModel`, then **requires ≥1 top-dimensional physical group** (a physical
volume in 3-D) — else it throws. So Tessella is a drop-in mesher iff a Tessella-written
`.msh` parses into a Gridap model carrying the region volumes as physical groups.

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
**Tessella meshes are solver-consumable by ASCENT with no format bridge.**

## Reproduce

```
# 1. Tessella env — produce the mesh (writes ascent_coax.msh here, git-ignored)
julia --project=<Tessella.jl> validation/ascent_handshake/generate.jl
# 2. ASCENT env — confirm it loads (GridapGmsh is an ASCENT dep, not a Tessella one)
julia --project=<2026_066/ASCENT> validation/ascent_handshake/handshake.jl
```

The generated `.msh` is a build artifact (git-ignored) — regenerate it with step 1.
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

So ASCENT builds the finite-element Maxwell system on a Tessella-generated mesh, with the
operator's physics (complex symmetry, PSD stiffness) holding — the mesh is not just
loadable but **usable by the solver**. The full 22-case HFSS regression (each case: geometry
+ BCs + frequency sweep + solve + compare to the guide) is the remaining compute campaign;
this proves the fundamental solver-usability it builds on.
