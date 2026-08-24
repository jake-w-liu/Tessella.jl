# Native mesh-CSG — `mesh_boolean` validation

`mesh_boolean(A, B, op)` computes the `:union` / `:intersection` / `:difference`
of two closed, watertight, manifold triangulated solids **natively — no
OpenCASCADE** — returning the closed, watertight, outward-oriented result surface
(fillable by `recover_boundary`), or throwing an explicit blocker on an
unsupported degeneracy (never a silently wrong surface).

Two exact engines, auto-dispatched:
- **axis-aligned** inputs (boxes / box-unions, incl. every coplanar shared-face
  case) → exact plane-arrangement path (grid-cell classification);
- **general position** (e.g. box × cylinder) → Cork/libigl-style **exact `orient3`
  tri-tri intersection** with `Rational{BigInt}` seam points (bit-identical shared
  seam ⇒ watertight) → `Mesh2D` constrained-Delaunay re-triangulation of cut
  regions → ray-cast inside/outside classification → per-op assembly.

`boolean_volumes.jl` checks (independent divergence-theorem volume oracle):

| case | union | inter | diff |
|---|---|---|---|
| box∩box (x-offset, coplanar faces) | 96 | 32 | 32 |
| box∩box (xy-offset) | 112 | 16 | 48 |
| box∩box (general position) | 1875 | 125 | 875 |
| box − cylinder (tri-tri+CDT path) | — | — | exact faceted |

Every result: **exact volume, watertight, `recover_boundary`-fillable**. The
regression form lives in `test/meshing/mesh3d_test.jl` ("mesh_boolean: native mesh-CSG").

Run:

```
julia --project=. validation/csg/boolean_volumes.jl
```
