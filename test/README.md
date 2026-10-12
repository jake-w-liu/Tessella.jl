# Testing

The default suite checks functionality and correctness with small real public
calls: geometric predicates, mesh validation and ownership, transforms, entity
topology, NURBS and fields, `.geo` parameter reading, bounded 1D/2D/3D meshing,
structured patches, refinement and recombination, quadratic tetrahedra,
MSH/STL I/O, API state, query ownership and failure atomicity, plus small matching
and compact surface coverage regressions, list/model-view format checks and
independent multi-root CAD boundaries, nested and mixed mapped/product STEP
product-usage forest paths
and planar Face holes, IGES143
boundaries and explicit plane UV/rigid-transform chains, complete quadratic
Hex, Pyramid13 and orders3–10 Curve VTK export, actual two-Tet partition
interfaces, closed Tri3-BRep seam transport and ghost queries, declared
Curve/Point lower-BRep partition topology, a mixed linear family
(Tet4+Prism6 over Tri3+Quad4) BRep partition, a complete second-order
(two-Tet10 over Tri6) BRep partition, a three-Tet4 junction BRep partition
with junction seam and points, a coincident-Surface BRep partition, a
complete third-order (two-Tet20 over Tri10) BRep partition, a periodic
link and parametric-node BRep partition transport, a parallel
`General.NumThreads` BRep partition, a complete fourth-order
(two-Tet35 over Tri15) BRep partition, complete
catalogued high-order Surface/Volume VTK, directional constant-tensor
box and shared-material grids, floating-point Pyramid options, CLI help/version
and read-only MSH checking, plus atomic current-model opening and merging,
including one POS point import and an empty-file merge, plus derived two-step
view metadata and display selection, two list/model replacement checks, a
quarter-turn periodic Surface copy, OCC circle
`add_circle_arc`/`add_circle` creation, synchronization and evaluation,
and OCC wire/planar-face `add_wire`/`add_plane_surface`/`add_rectangle`
creation, synchronization, bounds and boundary checks. It needs
no external mesher or solver.

Run the default suite from the repository root. The runner automatically uses
minimal compilation for small functional checks. Ordinary `Pkg.test()` and
standalone invocation relaunch once with the same active project when needed,
before loading Tessella. Child failures propagate. Fixed custom coverage
filenames retain the original invocation; normal coverage flags are inherited.
This keeps large production dispatcher optimization out of routine tests:

```powershell
julia --project=. --startup-file=no --check-bounds=yes --compile=min -O0 test/runtests.jl
```

`Pkg.test()` runs the same default suite in Julia's package test environment:

```powershell
julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

For checks of the normal optimized compiler path, explicitly opt in to the
same assertions with normal compilation:

```powershell
julia --project=. --startup-file=no -O2 -e 'ENV["TESSELLA_TEST_COMPILER"]="optimized"; using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "-O2"])'
```

Local receipts on 2026-10-13, with bounds checks; only the marked
`test/runtests.jl` rows ran under `julia -t4` to exercise
`General.NumThreads`, the `Pkg.test()` rows on one Julia thread:

| Runtime | Mode and entry point | Passing assertions | Elapsed |
| --- | --- | ---: | ---: |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 85.5 s | 369 | ~88 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 85.0 s | 369 | ~120 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 88.7 s | 369 | ~92 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 90.8 s | 369 | ~124 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 83.1 s | 364 | ~86 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 83.8 s | 364 | ~119 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 87.2 s | 364 | ~90 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 89.1 s | 364 | ~122 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 80.8 s | 357 | ~84 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 82.6 s | 357 | ~117 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 82.8 s | 357 | ~86 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 84.6 s | 357 | ~119 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 82.0 s | 350 | ~86 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 83.3 s | 350 | ~118 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 86.6 s | 350 | ~90 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 90.3 s | 350 | ~124 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 81.1 s | 342 | ~85 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 81.7 s | 342 | ~118 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 89.8 s | 342 | ~93 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 87.8 s | 342 | ~122 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 95.6 s | 329 | ~128 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 105.3 s | 329 | ~140 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 95.5 s | 329 | ~128 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 94.9 s | 329 | ~127 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 78.6 s | 321 | ~112 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 79.8 s | 321 | ~115 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 88.4 s | 321 | ~92 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 87.4 s | 321 | ~90 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 75.2 s | 313 | ~110 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 72.5 s | 313 | ~110 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 74.6 s | 313 | ~80 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 76.5 s | 313 | ~82 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 77.1 s | 306 | ~110 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 78.6 s | 306 | ~110 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 86.2 s | 306 | ~88 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 88.2 s | 306 | ~90 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 78.9 s | 302 | ~110 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 76.4 s | 302 | ~110 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 79.7 s | 302 | ~82 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 79.3 s | 302 | ~82 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 82.2 s | 301 | ~85 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 83.8 s | 301 | ~87 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 76.8 s | 301 | ~80 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 79.6 s | 301 | ~82 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 91.8 s | 300 | ~94 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 83.5 s | 300 | ~86 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 87.7 s | 300 | ~90 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 85.2 s | 300 | ~88 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 74.0 s | 298 | ~77 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 76.7 s | 298 | ~79 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 78.7 s | 298 | ~81 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 77.4 s | 298 | ~79 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 78.6 s | 296 | ~81 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 78.0 s | 296 | ~80 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 78.9 s | 296 | ~81 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 77.2 s | 296 | ~79 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 73.8 s | 294 | ~76 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 66.2 s | 294 | ~69 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 73.2 s | 294 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 74.0 s | 294 | ~76 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 52.3 s | 292 | ~55 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 51.3 s | 292 | ~53 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 51.1 s | 292 | ~53 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 50.0 s | 292 | ~52 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 55.8 s | 291 | ~58 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 53.8 s | 291 | ~56 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 51.6 s | 291 | ~53 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 50.6 s | 291 | ~52 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 58.1 s | 288 | ~60 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 55.3 s | 288 | ~57 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 54.8 s | 288 | ~56 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 51.8 s | 288 | ~53 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 45.9 s | 286 | ~48 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 44.4 s | 286 | ~46 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation `-t4`, body 48.8 s | 286 | ~51 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation `-t4`, body 48.1 s | 286 | ~50 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 53.3 s | 282 | ~56 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 52.8 s | 282 | ~55 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 52.2 s | 282 | ~54 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 50.8 s | 282 | ~53 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 42.3 s | 279 | ~45 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 46.5 s | 279 | ~49 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 44.4 s | 279 | ~47 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 47.2 s | 279 | ~50 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 41.3 s | 272 | ~42 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 45.2 s | 272 | ~46 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 41.0 s | 272 | ~42 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 43.5 s | 272 | ~44 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 40.5 s | 269 | ~41 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 44.0 s | 269 | ~45 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 40.5 s | 269 | ~41 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 43.3 s | 269 | ~44 s |
| Julia 1.12.7 | `Pkg.test()`, automatic minimal compilation, body 40.7 s | 266 | ~41 s |
| Julia 1.13.1 | `Pkg.test()`, automatic minimal compilation, body 43.9 s | 266 | ~44 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 40.7 s | 263 | ~41 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 41.4 s | 263 | ~42 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 38.1 s | 259 | ~39 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 41.5 s | 259 | ~42 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 37.5 s | 256 | ~70 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 40.8 s | 256 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 36.3 s | 249 | ~70 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 41.7 s | 249 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 34.2 s | 246 | ~70 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 38.7 s | 246 | ~75 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 31.9 s | 243 | ~65 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 35.6 s | 243 | ~70 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 24.4 s | 240 | ~60 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 27.2 s | 240 | ~60 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 17.2 s | 237 | ~50 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 16.9 s | 237 | ~40 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 15.7 s | 234 | ~50 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 15.2 s | 234 | ~40 s |
| Julia 1.12.7 | `test/runtests.jl`, automatic minimal compilation, body 15.4 s | 231 | ~50 s |
| Julia 1.13.1 | `test/runtests.jl`, automatic minimal compilation, body 14.8 s | 231 | ~40 s |
| Julia 1.12.7 | Ordinary `Pkg.test()`, automatic minimal compilation | 224 | 45.21 s |
| Julia 1.13.1 | Ordinary `Pkg.test()`, automatic minimal compilation | 224 | 38.38 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 222 | 45.18 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 222 | 38.74 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 220 | 47.05 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 220 | 40.34 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 218 | 45.41 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 218 | 38.66 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 214 | 49.28 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 214 | 38.93 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 209 | 46.49 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 209 | 39.23 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 207 | 45.23 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 207 | 37.67 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 203 | 46.00 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 203 | 39.12 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 201 | 44.88 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 201 | 36.53 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 192 | 47.14 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 192 | 39.75 s |
| Julia 1.12.7 | Preceding ordinary package invocation | 190 | 47.66 s |
| Julia 1.13.1 | Preceding ordinary package invocation | 190 | 37.80 s |
| Julia 1.12.7 | `--compile=min -O0`, preceding `Pkg.test()` | 180 | 38.97 s |
| Julia 1.13.1 | `--compile=min -O0`, preceding `Pkg.test()` | 180 | 42.20 s |
| Julia 1.13.1 | `--compile=min -O0`, earlier CLI `test/runtests.jl`, warm | 173 | 5.71 s |
| Julia 1.13.1 | `--compile=min -O0`, preceding CAD/model/periodic `Pkg.test()` | 171 | 40.61 s |
| Julia 1.13.1 | `--compile=min -O0`, earlier CLI/list increment `Pkg.test()` | 163 | 36.41 s |
| Julia 1.12.7 | `--compile=min -O0`, earlier `Pkg.test()` | 157 | 37.03 s |
| Julia 1.13.1 | Normal `-O2`, earlier `Pkg.test()` | 157 | 217.47 s |
| Julia 1.12.7 | Normal `-O2`, earlier standalone diagnostic | 157 | 301.75 s |

The current assertion bodies ran onJulia1.12.7/1.13.1.
Cold elapsed times include about30/25 s of precompilation.
These ordinary invocations verify the automatic runner in each package test
environment on the same final Source. The latest increment adds a
compact OCC wedge check (`occ_wedge_smoke`: a standard `add_wedge`
synchronized into its trapezoid-prism volume and centroid, plus a
truncated `ltx>0` wedge) — the
`test/interfaces/api_occ_wedge_test.jl` file covers the
`BRepPrim_Wedge` profile: the `ltx==0` triangular-ridge 5-face
layout, the `ltx≠0` truncated 6-face layout, the `ltx==dx` box limit,
`ltx>dx` overhangs, signed extents with outward-consistent winding,
`zAxis` frames, analytic `get_mass`/`get_center_of_mass` values
(`dy·dz·(dx+ltx)/2`), outward normal sets, manifold edge sharing,
"Degenerate wedge"/finite/`zAxis`/tag validation with failure
atomicity, recursive-removal surface-loop cleanup, kernel-side
transform visibility before `synchronize` and end-to-end `generate(3)`
tetrahedralization. The preceding increment adds a
compact OCC measure check (`occ_measures_smoke`: an `add_box` and
`add_cylinder` with analytic volumes and centroids plus a line's
arclength midpoint) — the
`test/interfaces/api_occ_measures_test.jl` file covers kernel-side
BRepGProp measures matching upstream `getMass`/`getCenterOfMass`
semantics: line/circle closed forms, ellipse/spline/iso curve panel
quadrature, plane shoelace with conic segment corrections and the
`|outer|−|holes|` convention, quadric face patch integrals, NURBS
knot-span tensor quadrature, primitive-record and torus-sector solid
formulas, the all-plane divergence path (volume plus second-moment
centroid), partial torus, degenerate/empty/non-1:3 entities,
transformed-entity visibility before `synchronize`, the conic-faced
get_mass/get_center_of_mass split and upstream error messages. The
preceding increment adds a
compact OCC query check (`occ_queries_smoke`: an `add_box` with
per-dimension bounding boxes, entity containment filtering,
surface/curve loop traversal order and a `mesh.set_size` point
attribute transported on sync) — the
`test/interfaces/api_occ_queries_test.jl` file covers kernel-side
registry reads matching upstream `OCC_Internals` semantics:
per-dimension and curved-record bounding boxes (sphere ball box,
torus extrema, circle/ellipse boxes), transform visibility before
`synchronize`, containment and dimension filters, multi-wire loop
order, standalone loop/volume reporting, silent/overwriting
`set_size` behavior and upstream error messages. The preceding
increment adds a
compact OCC transform check (`occ_transform_smoke`: an `add_box`
translated with its sub-entities re-created under fresh tags while a
holder curve keeps its tag, `copy` rebinding the solid under tag 2,
`dilate` on the copy, and removal pruning its surface loop) — the
`test/interfaces/api_occ_transform_test.jl` file covers
`translate`/`rotate`/`dilate`/`mirror`/`symmetrize`/`affine_transform`/
`copy` upstream `_transform`/`_multiBind` semantics: tag churn with a
holder entity, ownership-protected shared-vertex detach, rigid motion
and rotated compact records, isotropy gates and atomic rejection,
reflected sphere poles through degenerate edges, shear absorption by
plane records, copy explosion (wire→edges, point→fresh tag, duplicate
input dedup, highest-dim-only output, compact-encoding duplication),
guard paths, surface-loop lifecycle pruning, idempotent
re-synchronization and end-to-end `generate(3)` containment on a
translated cylinder. The preceding increment adds a
compact OCC primitive-solid check (`occ_primitives_smoke`: a cylinder
plus sphere synchronized into typed faces with retained compact
encodings, and a standalone torus volume) — the
`test/interfaces/api_occ_primitives_test.jl` file covers the
`BRepPrim_OneAxis` topology (rim/pole vertices, rim circles, seam lines,
degenerate apex/pole edges, signed wires and shells), analytic record
contents, pcurve→3-D round-trips, outward normals, parametrization
bounds, frustum and both apex variants, `zAxis`, tag/collision/finite
guards, failure atomicity, ownership-aware removal with encoding
cleanup, idempotent re-synchronization and end-to-end `generate(3)`
containment on cylinder and sphere volumes. The preceding increment adds a
compact OCC shell/solid check (`occ_volume_smoke`: an `add_box`
synchronized into one volume, one surface loop, six surfaces, twelve
curves and eight points, plus a standalone `add_surface_loop`/`add_volume`
pair) — the `test/interfaces/api_occ_volume_test.jl` file covers
box topology with shared edges and outward normals, opposite-corner
signed extents, `Degenerate box` and finite/tag guards, cavity shell
ordering, idempotent loop synchronization, ownership-aware removal and
end-to-end `generate(3)` tetrahedral meshing. The preceding increment adds a
compact OCC ellipse check (`occ_ellipse_smoke`: an `add_ellipse_arc`
quarter-arc, a full `add_ellipse` and an `add_disk` face synchronized
into `:ellipse`/`:plane` records with parametrization bounds and
point evaluation through the generic arithmetic paths) — the
`test/interfaces/api_occ_ellipse_test.jl` file covers the
`EllipseD1`/`EllipseD2` evaluators, axis forms, the two-point arc
solver's null/symmetric/axis-swap behavior, disk topology with
`ellips2d` pcurve round-trips, dilation and reversed-hole loops,
coherence merging, trimmed-extrema projection and end-to-end disk
meshing, plus the `add_plane_surface` `abs`-wire-tag and empty-input
parity fixes. The preceding increment adds a
compact OCC wire/planar-face check (`occ_plane_smoke`: a lines-only
square wire synchronized into a `:plane` face with UV bounds,
evaluation and boundary queries, plus `add_rectangle`'s published
sub-entities) — the 153-assertion
`test/interfaces/api_occ_plane_test.jl` covers signed chaining,
first-tag reversal, `checkClosed`, hole loops, annulus `circ2d`
pcurves, rounded corners and failure atomicity. The
preceding increment adds a
same-pair mixed-fragment check (compact: a crossing plus chord imprint
on one solid pair decomposes into {tet}, {box} and {tet∩box} Regions)
to `test/functional/brep_assembly_smoke.jl` and extends the assembly
functional file to 397 assertions on both runtimes — the fragment's
region classifier now resolves every collected contact class spatially
(welds, coplanar overlays, imprints, enclosures, interpenetrating
volumes and pinched region vertices self-touching at isolated points),
composing mixed arrangements under both fill modes. The
preceding increment adds an
off-plane T-junction check (compact: a tet vertex inside the box's top
face retessellates the facet around a Steiner vertex, unwelded) to
`test/functional/brep_assembly_smoke.jl` and extends the assembly
functional file to 342 assertions on both runtimes — isolated
non-coplanar contacts imprint points as Steiner vertices and chords as
constrained sub-edges into every touched facet lacking them, staying
unwelded with each solid's own coincident node, while same-pair
crossing+imprint and enclosed-contact mixtures still reject. The
preceding increment adds a
coplanar-imprint check (compact: a partial-face overlap welds the exact
patch into a shared Surface bound by both Regions) to
`test/functional/brep_assembly_smoke.jl` and extends the assembly
functional file to 315 assertions on both runtimes — non-conforming
contacts on one coplanar support plane (partial faces, T-junctions,
crossed rim edges, mismatched shared-quad diagonals) retessellate over a
shared exact planar arrangement into welded coincident pieces, while
aligned-winding overlap, crossing/enclosure mixtures and near-coplanar
inputs still reject. The preceding increment adds a
scaled-Cartesian-operator check (compact: scale-2 mapped target produces
the 2-edge box at the operator origin, scale-0 rejects) to
`test/functional/brep_mapped_occurrences_smoke.jl`, extends the assembly
functional file to 241 assertions on both runtimes (IGES solid-root
type-124 affines — scale, mirror, shear — admit exactly, singular
rejects), and extends the Cartesian file to 54 assertions — nonzero
Cartesian-operator scale applies `scale × R_axes` (mirrored scale flips
winding), nested compositions certify scaled-rigid. The preceding
increment adds an
assembly transversal-crossing check (compact: three Regions welded at two
shared Surfaces, `valid`,
`test/fixtures/brep_crossing_boxes.step`) to
`test/functional/brep_assembly_smoke.jl` and extends the assembly
functional file to 236 assertions on both runtimes — transversally
crossing boundaries fragment into containment-signature Regions
(A only, B only, A∩B) whose coincident cut interfaces weld under both
fill modes, while coplanar overlaps, T-junctions and
crossing-plus-weld/enclosure mixtures still reject. The preceding
increment adds an
assembly nested-enclosed check (compact: filled product is `valid`,
enclosed boundary welds into one shared Surface bound by both Regions,
`test/fixtures/brep_nested.step`) to
`test/functional/brep_assembly_smoke.jl` and extends the assembly
functional file to 202 assertions on both runtimes — strictly enclosed
solids admit as nested Regions (verbatim under `fill=false`; one joint
conforming fill under `fill=true` tiling each Region exactly), partial
interpenetration and enclosed-contact still reject. The preceding
increment adds an
assembly conforming vertex/edge-contact check (compact: admitted
unwelded, each solid keeps its coincident boundary node,
`test/fixtures/brep_contact_vertex.step`) to
`test/functional/brep_assembly_smoke.jl` and extends the assembly
functional file to 178 assertions on both runtimes — point/segment
contacts are conforming exactly at coincident discrete topology while
T-junctions and partial-edge overlap still reject. The preceding
increment adds an
assembly conforming-weld testset (compact: one shared Surface entity bound
by both Regions, one emitted interface cell, deduplicated coordinates,
`valid`) to `test/functional/brep_assembly_smoke.jl` and extended the
assembly functional file to 164 assertions on both runtimes — vertex-
identical coincident facets with opposed winding weld into a shared
Surface (`test/fixtures/brep_contact_twin.step` twin tetrahedra plus a
same-diagonal quad-face weld over two coincident facet pairs),
mismatched-diagonal/partial-quad fills, three-solid facets and
rim-crossing cells still reject atomically, and an MSH4 round-trip
preserves the shared entity bindings. The preceding increment adds an
assembly admit-plus-reject testset (2 assertions: genuinely disjoint
solids with overlapping bounding boxes import as classified independent
mixed meshes, a coincident-shell twin rejects) to
`test/functional/brep_assembly_smoke.jl` plus a 30-assertion geometric
independence set in `test/geometry/brep_assembly_functional_test.jl`
(exact triangle-intersection and
ray-parity enclosure checks in `Rational{BigInt}` — vertex/edge/partial
face/coincident-facet contact, interpenetration and enclosure all reject
atomically, verified by shifted-box mutants of
`test/fixtures/brep_disjoint_overlap.step`). The preceding increment adds an
order ≥4 BRep partition testset (4 assertions: two Tet35 sharing a Tri15
facet — Tri15 interface cell, three generated Line5 seam cells,
`validate`, unpartition restore) to the declared-BRep default smoke
plus a 39-assertion complete order ≥4 set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
412 assertions on both runtimes (admitted types derive from
`MSH_CATALOG` so every catalogued complete family partitions at any
order — Tet35/Hex125/Prism75/Pyr55 pairs with Tri15/Quad25 interfaces
and Line5 seams plus an order-5 Tet56 pair with Tri21/Line6 output,
canonical facet orientation via monomial-lattice permutations, ghosts,
MSH4 round-trips, unpartition, byte-identical threaded output and
atomic rejection of Quad12, Tri12 and Tet22 blocks). The preceding
increment adds a
parallel BRep partition testset (3 assertions: `General.NumThreads=3`
lower-BRep partition validates, the Region-parented interface carries
both partitions, unpartition restore) to the declared-BRep default smoke
plus a 9-assertion parallel-determinism set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
373 assertions on both runtimes (`NumThreads=1` vs 3/4 produce identical
`partition_records` and byte-identical MSH4 on both producers under
`julia -t4`). The preceding increment adds a
periodic/parametric BRep partition testset (3 assertions: an imported
`$Periodic` link and UV-parametric vertices survive partition on their
original entity tags, `validate`, unpartition restore) to the
declared-BRep default smoke plus a 28-assertion periodic set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
364 assertions on both runtimes (a discrete periodic prism fixture — an
imported dim-2 `$Periodic` link with an affine transform, six
UV-parametric boundary vertices, Tet4+Tri3 cells — verifies verbatim link
transport through partition, ghosts, MSH4 round-trips and `unpartition`,
same-dimension same-owner parametric reclassification, and the bare-Tet4
producer path). The preceding increment adds a
third-order BRep partition testset (4 assertions: two Tet20 sharing a Tri10
facet — Tri10 interface cell, three generated Line4 seam cells, `validate`,
unpartition restore) to the declared-BRep default smoke plus a 40-assertion
complete third-order set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
336 assertions on both runtimes (Line4 curves, Tri10/Quad16 surfaces and
Tet20/Hex64/Prism40/Pyr30 volumes partition through the same declared-BRep
certificate — a Tet20 pair rimmed by a declared closed Line4 Curve, a
Hex64 pair with a closed four-cell Line4 seam, a Prism40 pair mixing
Tri10/Quad16 exteriors and a Pyr30 pair sharing the Quad16 base, plus
ghosts, MSH4 round-trips, unpartition and atomic rejection of Tri15,
Quad12 and Tet35 blocks). The preceding increment adds a
coincident-support BRep partition testset (3 assertions: two Regions sharing
a facet carried by one declared Surface — the partition child merges [1,2]
membership, binds both Region children with opposite signs, no duplicate
interface entity or cell, unpartition restore) to the declared-BRep default
smoke plus a 22-assertion coincident set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
296 assertions on both runtimes (the merged child keeps cell and UID,
declared Curve/Point rim roles persist, MSH4 round-trips and the atomic
one-sided rejection are covered). The preceding increment adds a junction
BRep partition testset (3 assertions: three Tet4 meeting at a junction rim
edge, two dim-0 junction entities, one junction seam carrying all three
partitions, unpartition restore) to the declared-BRep default smoke plus a
40-assertion junction set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
274 assertions on both runtimes (multi-volume junction fans partition
through the manifold rim certificate; junction seams share all incident
interfaces under Surface-or-Region parents; Point15 junction entities reuse
declared Point children; a carrier-free four-volume fan, a declared
Curve/Point junction and a Tet10 junction with a Line3 seam are covered
alongside MSH4 round-trips and atomicity). The preceding increment adds a second-order
lower-BRep partition testset (3 assertions: two Tet10 sharing a Tri6 facet,
Tri6 interface cell, three generated Line3 seam cells, unpartition restore) to
the declared-BRep default smoke plus a 46-assertion complete second-order set
in `test/interfaces/api_partition_lower_brep_functional_test.jl`, now
234 assertions on both runtimes (Line3 curves, Tri6/Quad9 surfaces and
Tet10/Hex27/Prism18/Pyr14 volumes partition through the same declared-BRep
certificate, including declared Line3 border edge-node checks, ghosts, MSH4
round-trips, unpartition and atomic violations for foreign mids, serendipity
families and reversed volumes). The preceding increment adds a mixed
linear family BRep partition testset (4 assertions: Tet4+Prism6 sharing a
Tri3 face, Tri3 interface, generated seams, unpartition restore) to the
declared-BRep default smoke plus a 61-assertion mixed-family set in
`test/interfaces/api_partition_lower_brep_functional_test.jl`, now
188 assertions on both runtimes (Quad4 surfaces and Hex8/Prism6/Pyr5
volumes partition through the same declared-BRep certificate, including
ghosts, MSH4 round-trips and atomic violations). The preceding increment
adds three
cheap MSH4 signed-boundary import checks to the declared-BRep default
smoke (10 assertions) plus a 27-assertion import testset in
`test/interfaces/api_partition_lower_brep_functional_test.jl`. The preceding increment adds a declared
Curve/Point lower-BRep partition testset (7 assertions) on that fixture's
original 100 assertions. The preceding increment added three cheap
meridian-certified `TOROIDAL_SURFACE` segment checks;
`test/geometry/brep_step_torus_segment_functional_test.jl` is the opt-in
segment fixture (boundary-only — toroidal fills remain a hard tetrahedralizer
input). The preceding increment adds three cheap
rim/apex-certified `CONICAL_SURFACE` checks; `test/geometry/brep_step_cone_functional_test.jl`
is the opt-in frustum/cone fixture set. The preceding increment adds three cheap
rim-certified `CYLINDRICAL_SURFACE` checks; `test/geometry/brep_step_cylinder_functional_test.jl`
is the opt-in cylinder fixture. The preceding increment adds three cheap
closed-face `SPHERICAL_SURFACE` checks; `test/geometry/brep_step_sphere_functional_test.jl`
is the opt-in sphere fixture. The preceding increment adds three cheap
STEP unit-conversion checks; `test/geometry/brep_step_units_functional_test.jl`
is the opt-in conversion fixture set. The preceding increment added three cheap
mixed STEP forest checks; `test/geometry/brep_step_relationship_mixed_functional_test.jl`
is the opt-in mixed fixture. The preceding increment added four cheap
quarter-turn periodic Surface-copy checks; `test/interfaces/api_rigid_periodic_hybrid_functional_test.jl`
is the opt-in rigid/hybrid fixture. The preceding224-check increment added two cheap
list/model representation checks. The preceding218-check increment added two
file-merge and two read-only CLI assertions; periodic Volume, tensor production, native and solver
campaigns remain opt-in. Portable final224 receipts and Source hashes are also
committed under `test/artifacts/handoff_20261011`; see the
[Devin continuation checklist](../DEVIN_HANDOFF.md). The preceding214-check bodies took11.1/10.9 s and
qualify their own Source. The preceding 209-check
bodies took 11.1/10.7 s and qualify their own Source. The preceding 207-check
bodies took 10.1 s on both runtimes and qualify their own Source. The preceding203-check
bodies took10.4 s on both; their receipts remain separate. The preceding 201-check
bodies took 9.7/9.1 s; their receipts remain separate. The earlier 192-check
bodies took 9.3/9.7 s and remain recorded separately.
The earlier190-check
bodies took10.3/8.8 s. The earlier 180-assertion
bodies took 7.2/6.2 s; those explicit minimal-compilation receipts are historical.
The earlier 173-assertion standalone run took 5.71 s overall, with a 5.3 s
body. The preceding 171-assertion package run also had a 5.3 s body. The earlier 163-assertion check took 1.4 s.
The earlier 157-assertion routine check took
1.6 s. The earlier optimized
package receipt includes 24.3 s of precompilation and 188.4 s of tests. The
earlier standalone diagnostic used the same assertions with extra API progress
messages and included 36.14 s of precompilation. Those optimized receipts
precede the final narrow native API corrections, which were checked separately
with focused public regressions. API generation dominated normal compilation
time. Both commands exercise the same public operations and assertions; these
timings include their compilation costs.

Keep default fixtures small. Check observable behavior against analytic values,
valid topology, ownership and explicit errors. Exhaustive combinations, large
native meshes, resource/benchmark campaigns, pinned parity and source audits,
and application/solver campaigns belong in opt-in validation.

The complete historical regression files remain in their existing directories.
An extended runner combines the original and newer regression lists:

```powershell
julia --project=. --startup-file=no --check-bounds=yes test/extended/runtests.jl
```

This extended suite contains exhaustive predicate/configuration sweeps,
allocation and source audits, native meshing stress cases and all 22 HFSS
geometries. It can take hours. Run a specific existing regression file when
investigating its feature instead of automatically running every campaign.
Full `.geo` execution is opt-in too: even the two-point For/If fixture in
`extended/geo_functionality.jl` compiles the full native dispatcher. Default
geometry tests exercise entity operations and the parameter reader directly.
`extended/legacy_runtests.jl.txt` and `extended/devin_runtests.jl.txt` preserve
the previous runners verbatim for historical reference. The extended runner
omits repository-layout enforcement and the static development-stage assertion.

External Gmsh comparison campaigns are separately opt-in:

```powershell
julia --project=. --check-bounds=yes validation/run_all.jl
```

Consult the relevant validation directory for its dependencies and resource
requirements. Archived artifacts and validation evidence are retained.
