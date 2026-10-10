# Literal ENC-COAX geometry — meshed natively from the `.geo` primitives

Status reviewed **2026-10-10**: this page describes the dated literal-enclosure reconstruction and external solver capture.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

**Historical capture verified 2026-08-13.** Tested gmsh 4.13/4.15 runs on the
enclosure coax fixture produced **0 volume tets**. This is a captured oracle
outcome, not a claim that every current Gmsh configuration fails. Its `.geo`
(`test/fixtures/enclosure_coax_junction.geo`)
is fully parametric — every solid is a gmsh `Box` or `Cylinder` primitive.

`reconstruct.jl` **parses those literal primitives directly** (no OpenCASCADE evaluation) and
reconstructs four physical volumes — **air cavity, metal case shell, coax pin, and slot** — with
Tessella's native primitives at the **exact literal dimensions** from the fixture, then meshes them
as ONE conforming partition (`tetrahedralize_conforming_exact`, exact-kernel — robust to the
thin-pin cosphericity). Result:

```
air_inside = (0,0,0, 0.22,0.14,0.3)   case_outer = (-5e-4,-5e-4,-5e-4, 0.221,0.141,0.301)   pin r=8e-4 len=0.1589
mesh: tets=134 valid=true regions = pin/slot/air/case (all four filled)
surface BC groups: radiation=68 pin_pec=16 case_pec=26 resistor=1 p1=1
```

and the written MSH v4.1 carries the **complete literal physical-group structure** — 4 volumes
(coax_pin/slot/air/case) + 5 BC surfaces (radiation, coax_pin_pecskin, case_pecskin, resistor,
p1_surface) — all of which load into ASCENT's parser (`LITERAL_HANDSHAKE_OK`, 9/9 groups) —
solver-consumable in that captured environment, where the tested gmsh runs left volumes empty.

## Scope

This meshes the **four literal physical volumes** (pin/slot/air/case — every gmsh Box/Cylinder
*volume* primitive, at the exact fixture dimensions) and tags the **five literal BC surfaces**
(radiation, the pin/case PEC skins, and the resistor + p1 port patches), so the mesh carries the
`.geo`'s complete physical-group structure and loaded whole in the captured ASCENT
environment. The tested gmsh 4.13/4.15 runs produced 0 tets. What is *not* reproduced is the **exact OpenCASCADE `BooleanFragments` imprint
interfaces** (the precise curved shield/bore boundaries and the exact resistor/p1 surface geometry).
Those remain pending CAD-parity work and are represented here as topologically tagged BC surface
patches rather than OCC-imprinted geometry.

The script also resolves the literal `.geo` entity names and constructs its
`Distance`/`Threshold`/`Box`/`Min` background-field graph. Pass `--apply-fields` to refine with the
literal graph before writing. That full-resolution path is a scalability acceptance, not part of
the fast topology reconstruction reported above.

## Reproduce

```
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=/path/to/Tessella.jl validation/enclosure_literal/reconstruct.jl # parse + mesh + write .msh
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=/path/to/Tessella.jl validation/enclosure_literal/reconstruct.jl --apply-fields
julia --project=/path/to/ASCENT validation/enclosure_literal/handshake.jl      # confirm it loads in ASCENT
```

Use the Julia version and dependencies required by the external ASCENT project
for its handshake command. The generated `.msh` is a git-ignored build artifact.
The dated topology reconstruction and group-loading results do not qualify
BooleanFragments CAD parity, the field-refinement scalability path, or the
current combined release.
