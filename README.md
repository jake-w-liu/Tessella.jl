# Tessella.jl

A Julia-native, robust, memory-efficient **mesh generator** — built so that
*design → mesh always works*, the way HFSS meshing "just works". It replaces the
external `gmsh` dependency in the ASCENT / ascent-studio electromagnetics
toolchain, where gmsh's OpenCASCADE boundary recovery fails on hard geometries
(thin features, multi-way Boolean junctions) that a solver must nevertheless
mesh.

> **Status: a working pipeline (142k+ CRC-verified assertions), built from
> `startup.md` under the mandatory CRC discipline.** Implemented and gated:
> exact adaptive predicates (Shewchuk + SoS, vs an independent exact-rational
> oracle); compact mesh types + topology + deterministic checksums; gmsh `.msh`
> v2/v4 + STL I/O; **2-D** Delaunay + constrained Delaunay + Ruppert refinement;
> size fields + graded 1-D + planar/cylinder/parametric **surface** meshing;
> a robust **3-D** Delaunay kernel + **volume filling** (convex, non-convex,
> genus-1, thin, and multi-region — including a representative coax junction with
> **all volumes filled** where gmsh leaves them empty); tet quality reporting +
> Laplacian smoothing; surface healing (defect detection); and a top-level
> `mesh_volume` with a *validated-or-explicit-blocker* contract.
>
> **Remaining** (tracked honestly in `STATUS.md`): conforming interface recovery
> and sliver exudation for the *literal* enclosure `.geo`; the OpenCASCADE / CSG
> geometry kernel; high-order elements; ASCENT integration. See `STATUS.md` for
> the live stage board.

```julia
using Tessella
m = mesh_volume(surface)          # closed triangle surface → validated tet mesh
                                  #   (throws with a precise report if the surface is defective)
```

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
