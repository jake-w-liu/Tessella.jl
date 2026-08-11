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
| `schonhardt_blocker_check.jl` | the **Schönhardt** polyhedron (non-tetrahedralizable without a Steiner point) raises the explicit blocker; the convex triangulation of the same 6 points meshes — the blocker is exactly discriminating |

The regression form of both lives in `test/mesh3d_test.jl` ("recover_boundary:
conforming tetrahedralization of arbitrary PLCs"), which also covers a star-shaped
L-prism and a faceted (octagonal) cylinder.

**Open:** Steiner-point recovery for genuinely non-tetrahedralizable
(Schönhardt-type) inputs — those currently block rather than mesh.

Run either:

```
julia --project=. validation/boundary_recovery/<script>.jl
```
