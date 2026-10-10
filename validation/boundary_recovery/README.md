# Boundary recovery — `recover_boundary` validation

Status reviewed **2026-10-10**: this page describes the covered boundary-recovery fixtures and historical prototype findings.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

`recover_boundary(surface)` attempts to recover a closed PLC surface (possibly
**non-convex**) as a **geometrically conforming** interior tetrahedral mesh, or
throws an **explicit blocker** naming the first unrecovered facet. The geometric
gate permits a coplanar region to be re-triangulated with the same vertices;
that contract does not require each original triangle to survive as a tet face. It never returns a silently non-conforming mesh
(PLAN principle #4). Conformity/validity are decided by an exact `Rational{BigInt}`
geometric gate; the domain is filled by insertion-order retry (`perturb=false`)
plus a gated flat-tet drop, which sidesteps the cospherical zero-volume-tet
degeneracy that stops the base kernel on `box_tunnel`.

| script | checks |
|---|---|
| `conformity_check.jl` | box, non-convex **genus-1 through-tunnel**, **hollow shell** — tet-mesh boundary area == input surface area, exact domain volume, `validate.ok`, closed-manifold, ~1–2 s |
| `schonhardt_blocker_check.jl` | the **Schönhardt** polyhedron (non-tetrahedralizable without a Steiner point) raises the explicit blocker by default; the convex triangulation of the same 6 points meshes — the blocker is exactly discriminating. With **`steiner=true`** the (star-shaped) Schönhardt is meshed by fan-tetrahedralization from an interior kernel point |

The regression form of both lives in `test/meshing/mesh3d_test.jl` ("recover_boundary:
conforming tetrahedralization of arbitrary PLCs"), which also covers a star-shaped
L-prism and a faceted (octagonal) cylinder.

**Coverage boundary (measured):** `recover_boundary` conforms convex, non-convex
Delaunay-recoverable (tunnel, shell, L-prism, faceted cylinder, **and non-star
U-channel / comb / star**), and star-shaped non-tetrahedralizable (Schönhardt via
`steiner=true`). A **measured unhandled class** is *non-star AND reflex* (a **twisted
non-convex prism**, a concrete constructed case), which raises the **explicit
blocker** under both `steiner` modes (never a silent bad mesh; regression-pinned in
`mesh3d_test.jl`). A historical recovery prototype investigated boundary-Steiner insertion using
Gabriel encroachment, segments before faces, and rejection of circumcenters that
encroach on subsegments. The cited Murphy–Mount–Gable /
Cohen-Steiner–Colin-de-Verdière–Yvinec packing bounds were design references,
not a qualification of a completed production algorithm. The prototype
**measured a representation blocker**: Steiner points *on slanted creases* can
require exact positions that are not representable in Float64; their coordinates
need not be irrational. The **Float64 kernel rounds them off the feature** (5/11
crease midpoints of the twisted U-prism become non-collinear), so the exact
`_rb_gate` correctly rejects them and the loop never converges (68+ Steiner points,
injects non-manifold invalidity). The historical proposed remedy was an
exact-coordinate (`Rational{BigInt}`) 3-D Delaunay kernel, a substantial new
subsystem. This experiment does not prove that it is the only possible remedy.
A narrow extruded-prism column decomposition using strictly interior Steiner
points meshed twisted prisms up to θ≈25° in the prototype, but did not generalize
and was not integrated. The published-baseline explicit blocker is regression-pinned.

Run either:

```
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/boundary_recovery/conformity_check.jl
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/boundary_recovery/schonhardt_blocker_check.jl
```

The table's ~1–2 s figure is a historical observation, not a current performance
qualification. The covered fixtures and geometric conformity gate do not imply
exact native UID/cell parity or support for every closed PLC. Follow the handoff
for the remaining general recovery and hybrid-generation work.
