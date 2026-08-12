# Literal ENC-COAX geometry — meshed natively from the `.geo` primitives

**Verified 2026-08-13.** The enclosure coax feed-through is the project's reason to exist
(gmsh 4.13/4.15 produce **0 volume tets**). Its `.geo` (`test/fixtures/enclosure_coax_junction.geo`)
is fully parametric — every solid is a gmsh `Box` or `Cylinder` primitive.

`reconstruct.jl` **parses those literal primitives directly** (no OpenCASCADE evaluation) and
reconstructs the three main physical volumes — **air cavity, metal case shell, coax pin** — with
Tessella's native primitives at the **exact literal dimensions** from the fixture, then meshes them
as ONE conforming partition (`tetrahedralize_conforming_exact`, exact-kernel — robust to the
thin-pin cosphericity). Result:

```
air_inside = (0,0,0, 0.22,0.14,0.3)   case_outer = (-5e-4,-5e-4,-5e-4, 0.221,0.141,0.301)   pin r=8e-4 len=0.1589
mesh: tets=100 valid=true regions=Dict(pin=>24, air=>52, case=>24)   all three volumes filled = true
```

and the written MSH v4.1 loads into ASCENT's parser with all three volume physical groups
(`LITERAL_HANDSHAKE_OK`) — solver-consumable, the geometry gmsh leaves empty.

## Scope

This meshes the **three main physical volumes at literal scale** (the physics-carrying regions,
and the exact place gmsh fails). The full `.geo` additionally carries, via OpenCASCADE
`BooleanFragments`, tiny sub-features — the coax shield/bore cylinders, the 1 mm coupling slot, the
resistor rectangle, and the p1 port disk — plus their exact curved imprint interfaces. Producing
those *exact OCC-fragment* interfaces is the pure-Julia OpenCASCADE-library path that PLAN §1/§6
lists as an explicit non-goal; they are represented here (and in `validation/ascent_handshake/`) as
tagged BC surfaces rather than OCC-imprinted volumes.

## Reproduce

```
julia --project=<Tessella.jl> validation/enclosure_literal/reconstruct.jl        # parse + mesh + write .msh
julia --project=<2026_066/ASCENT> validation/enclosure_literal/handshake.jl      # confirm it loads in ASCENT
```

The generated `.msh` is a git-ignored build artifact.
