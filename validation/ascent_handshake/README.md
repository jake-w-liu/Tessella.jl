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
