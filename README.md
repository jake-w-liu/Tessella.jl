# Tessella.jl

A Julia-native, robust, memory-efficient **mesh generator** — built so that
*design → mesh always works*, the way HFSS meshing "just works". It replaces the
external `gmsh` dependency in the ASCENT / ascent-studio electromagnetics
toolchain, where gmsh's OpenCASCADE boundary recovery fails on hard geometries
(thin features, multi-way Boolean junctions) that a solver must nevertheless
mesh.

> **Status: scaffolding only.** No meshing is implemented yet. This repository
> currently holds the architecture (`PLAN.md`), the development/CRC discipline
> (`DEVELOPMENT.md`), the status tracker (`STATUS.md`), and a package skeleton.
> Development begins from **`startup.md`**.

## Why

During the ASCENT HFSS-UserGuide validation campaign, gmsh 4.13.1 **and**
4.15.2-git both fail to mesh a shielded-enclosure coax feed-through
(`Invalid boundary mesh (overlapping facets)` → every volume empty). It is a
geometry/boundary-recovery defect, **not** a memory limit (peak 3.3 GB), and it
resists every standard remedy (algorithm switches, OCC healing, tolerances). A
solver cannot depend on a mesher that fails on valid geometry. Tessella owns the
meshing pipeline in Julia so robustness is ours to guarantee.

## Scope (honest)

Tessella is **not** a literal port of gmsh (~340k LOC core; >1.5M with its
OpenCASCADE/TetGen/Netgen contrib), and it does **not** reimplement OpenCASCADE.
It owns the part that matters for solver robustness — the **meshing algorithms
and boundary recovery** — and consumes geometry via OCC interop initially. See
`PLAN.md` §1 for the full scope statement and staged roadmap.

## Layout

```
PLAN.md         architecture + rigorous staged port plan (read first)
DEVELOPMENT.md  CRC (Correctness–Robustness–Completeness) discipline — mandatory
STATUS.md       stage board + acceptance cases + log (update every session)
startup.md      entry point for a fresh development session
src/Tessella.jl  package skeleton
test/           test harness (CRC suites arrive per stage)
```

## License

TBD by the author.
