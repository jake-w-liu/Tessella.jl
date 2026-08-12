# Boundary recovery — `recover_boundary` validation

`recover_boundary(surface)` recovers an arbitrary closed PLC surface (possibly
**non-convex**) as a **conforming** interior tetrahedral mesh — every input facet
appears as a tet face (by vertex identity; a coplanar region re-triangulated with
the same vertices still conforms) — or throws an **explicit blocker** naming the
first unrecovered facet. It never returns a silently non-conforming mesh
(PLAN principle #4). Conformity/validity are decided by an exact `Rational{BigInt}`
geometric gate; the domain is filled by insertion-order retry (`perturb=false`)
plus a gated flat-tet drop, which sidesteps the cospherical zero-volume-tet
degeneracy that stops the base kernel on `box_tunnel`.

| script | checks |
|---|---|
| `conformity_check.jl` | box, non-convex **genus-1 through-tunnel**, **hollow shell** — tet-mesh boundary area == input surface area, exact domain volume, `validate.ok`, closed-manifold, ~1–2 s |
| `schonhardt_blocker_check.jl` | the **Schönhardt** polyhedron (non-tetrahedralizable without a Steiner point) raises the explicit blocker by default; the convex triangulation of the same 6 points meshes — the blocker is exactly discriminating. With **`steiner=true`** the (star-shaped) Schönhardt is meshed by fan-tetrahedralization from an interior kernel point |

The regression form of both lives in `test/mesh3d_test.jl` ("recover_boundary:
conforming tetrahedralization of arbitrary PLCs"), which also covers a star-shaped
L-prism and a faceted (octagonal) cylinder.

**Coverage boundary (measured):** `recover_boundary` conforms convex, non-convex
Delaunay-recoverable (tunnel, shell, L-prism, faceted cylinder, **and non-star
U-channel / comb / star**), and star-shaped non-tetrahedralizable (Schönhardt via
`steiner=true`). The **one unhandled class** — *non-star AND reflex* (a **twisted
non-convex prism**, a concrete constructed case) — correctly raises the **explicit
blocker** under both `steiner` modes (never a silent bad mesh; regression-pinned in
`mesh3d_test.jl`). Closing it needs TetGen-style boundary-Steiner recovery. A design→implement→verify
workflow **designed the complete algorithm** (conforming-Delaunay by Gabriel-
encroachment Steiner insertion — segments before faces, reject-circumcenter-that-
encroaches-a-subsegment; provably terminating by the Murphy–Mount–Gable /
Cohen-Steiner–Colin-de-Verdière–Yvinec packing bound) and **measured the precise
blocker**: the algorithm needs Steiner points *on slanted creases*, whose exact
positions are irrational; the **Float64 kernel rounds them off the feature** (5/11
crease midpoints of the twisted U-prism become non-collinear), so the exact
`_rb_gate` correctly rejects them and the loop never converges (68+ Steiner points,
injects non-manifold invalidity). **Closing this class therefore requires an
exact-coordinate (`Rational{BigInt}`) 3-D Delaunay kernel** — a substantial new
subsystem — not a localized addition. (A narrow extruded-prism column-decomposition
using strictly-*interior* Steiner points does mesh twisted prisms up to θ≈25°, but
does not generalize; not integrated.) The current safe behaviour — an **explicit
blocker, never a silent bad mesh** — is regression-pinned.

Run either:

```
julia --project=. validation/boundary_recovery/<script>.jl
```
