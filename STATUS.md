# Tessella.jl — Development Status Tracker

**Single source of truth for build status. Update every session. Nothing dropped.**

Package: `Tessella` · Julia ≥ 1.11 · goal: robust Julia-native mesh generator
replacing the gmsh dependency (see `PLAN.md`). CRC discipline mandatory
(`DEVELOPMENT.md`).

## Stage board

| Stage | Scope | State | CRC gate |
|---|---|---|---|
| 0 | Foundations: repo, CI, mesh types, `.msh` I/O, exact predicates | **DONE — gate green** | predicates vs exact-rational oracle ✓; `.msh` round-trip (v2↔v4) CRC-preserving ✓ |
| 1 | 2-D Delaunay + CDT + quality refinement | **DONE — gate green** | exact empty-circumcircle oracle ✓; 2n−2−h count + hull ✓; CDT constraints present + locally Delaunay ✓; Ruppert angle/area bound achieved, domain area preserved ✓ |
| 2 | 1-D edge meshing + parametric surface meshing | **DONE — gate green** | size-graded edge mesh vs analytic arc length ✓; planar face exact area + quality ✓; cylinder watertight + area convergence ✓; parametric area convergence ✓ |
| 3 | 3-D Delaunay + **robust boundary recovery** + slivers | **kernel + filling DONE; interface recovery + slivers WIP** | 3-D empty-circumsphere oracle ✓; convex/non-convex/genus-1/thin/multi-region fills validated ✓; representative coax junction all volumes filled ✓; conforming interface recovery + exact enclosure geometry OPEN |
| 4 | size fields + optimization | **partial** | tet quality report + Laplacian smoothing ✓ (volume/validity preserved, mean dihedral up, slivers down); 3-D sliver exudation + domain-bounded refinement WIP |
| 5 | geometry kernel (OCC interop / native CSG) + heal | **partial** | `Heal` surface-defect detection ✓; native primitives (box/cylinder/box-tunnel) ✓ + fill to exact volume; Boolean CSG + OCC interop WIP |
| 6 | high-order + ASCENT integration + 22-case regression | **partial** | quadratic (P2) tet generation ✓ (shared mid-nodes, exact-midpoint volume) + gmsh type-11 I/O ✓ + curved P2 onto a cylinder primitive ✓; general-geometry curving + ASCENT drop-in + 22-case regression WIP |

## Standing acceptance cases (regression, CRC-stamped when they pass)

| id | geometry | why | status |
|---|---|---|---|
| ENC-COAX | enclosure coax feed-through (HFSS 9.2), surface 86 shield/case/air junction | gmsh 4.13.1 **and** 4.15.2 fail (`overlapping facets`, all volumes empty); the project's reason to exist | **PARTIAL** — a *representative* 3-region coax junction (pin/case/air-gap meeting at the bore walls) is filled with all three volumes non-empty & validated (`mesh3d_test.jl`), the exact failure mode gmsh cannot. Meshing the *literal* `.geo` still needs the geometry kernel (Stage 5) to produce its boundary surfaces + conforming interface recovery. |
| THIN-SLOT | enclosure 1 mm coupling slot | thin-feature Boolean sliver | open |
| SPIRAL | 10.1 silicon spiral (thin swept traces, layered stack) | high-aspect thin conductors | open |
| ARRAY-PML | 5.7 endfire unit cell (periodic + PML) | conformal periodic faces | open |

Reference artifact for ENC-COAX: the failing gmsh `.geo` and the gmsh-API
diagnosis (surface 86 bbox ≈ (168.4,140,148.4)–(171.6,160.5,151.6) mm; empty
volumes = pin/case/air) are recorded in the session that created this repo; copy
the `.geo` into `test/fixtures/` at Stage 0.

## Verified facts carried in from the ASCENT campaign (2026-08-11)

- The enclosure meshing failure is **NOT memory** — a no-time-limit run peaked at
  3.3 GB. It is a geometry/boundary-recovery defect at the coax junction.
- gmsh `occ.healShapes()` is a **false positive** here: it dissolves the OCC
  volumes to 0 tets, so "no empty volumes" is vacuous. Always check the *actual*
  per-volume tetrahedron count. (This is a CRC-discipline lesson, encoded in
  `DEVELOPMENT.md`.)
- 5 gmsh 2-D/3-D algorithms, OCC fix options, `Geometry.ToleranceBoolean`, and
  geometry protrusion all fail to mesh the enclosure. Standard tools do not fix
  it → the robustness must come from Tessella's own boundary recovery + heal.

## CRC regression artifacts (stamped when green)

| artifact | checksum | oracle |
|---|---|---|
| Predicates | 142,141 assertions green; every `orient2/3`, `incircle`, `insphere` sign matches the independent exact-rational homogeneous-determinant oracle over exhaustive integer-grid degeneracies (collinear/coplanar/cocircular/cospherical) + fixed-seed random floats | `test/oracles.jl` (generic Laplace determinant, `Rational{BigInt}`) |
| Unit cube (6-tet Kuhn) | `nodes=8 tets=6 bbox=[(0,0,0)→(1,1,1)] dihedral(min,mean)=(0.7854,0.7854) radedge(min,mean)=(0.86603,0.86603) bfaces=12` · SHA-256 `7ea403054f05392f18b404a1f5f78b12d70d45d40c7b04ba8f8dc3e030d8f3f9` | analytic volume 1/6·6=1, Euler χ=1 (ball), boundary χ=2 (sphere), 12 boundary tris |
| 2-D Delaunay (random/grid/cocircular) | exact empty-circumcircle holds (0 violations, `incircle_sos`); triangle count = 2n−2−h vs independent monotone-chain hull; seed-independent (SoS) | `is_delaunay`, hull oracle, Euler χ=1 |
| 10×10 square refined to 20° (max_area 1.0, seed 1) | `nodes=91 tris=148 bbox=[(0,0,0)→(10,10,0)]` · SHA-256 `583c615df1862c8518bbda409347f109dc25f7f8f5362562badf160fe6af30c1` | min angle ≥ 20°, max tri area ≤ 1.0, domain area = 100 preserved, CDT locally Delaunay |

## Log

- **2026-08-11** — repo scaffolded (PLAN, STATUS, DEVELOPMENT, README, skeleton,
  startup). gmsh 4.13.1 source studied and archived as reference. Package name
  `Tessella.jl` chosen. Local git initialized; GitHub remote pending `gh auth
  login`. **No code implemented yet** — Stage 0 is the first development task
  (start from `startup.md`).
- **2026-08-11** — **Stage 0 complete, gate green.** Implemented:
  - `src/Predicates.jl` — adaptive exact `orient2/orient3/incircle/insphere`
    (Shewchuk A-stage float filter → `Rational{BigInt}` exact fallback) plus
    Simulation-of-Simplicity `*_sos` variants (never return 0). Fast path is
    allocation-free (verified with `@allocated == 0`). Verified against an
    independent Laplace-determinant oracle on exhaustive degeneracies.
  - `src/MeshTypes.jl` — compact column-major `Mesh` (3×N coords, Int32 cells),
    allocation-free `node()`, area/volume/dihedral/circumradius/radius-edge
    quality, boundary-facet extraction, Euler characteristic, manifold check,
    `validate` (positive-tet / degeneracy / non-manifold diagnostics), and the
    deterministic order-invariant SHA-256 `mesh_crc` checksum.
  - `src/IO.jl` — gmsh `.msh` v2.2 **and** v4.1 ASCII read/write (round-trip
    preserves connectivity CRC across versions), ASCII+binary STL ingest with
    vertex welding, `.geo` parameter/physical-group scanner (no OCC eval — that
    is Stage 5). `read_geo_params` on the enclosure fixture recovers all three
    volume groups (air/case/coax_pin) + sizing + seed.
  - Full `Pkg.test()` green: **142,141 assertions**. `stage()` bumped to 1.
- **2026-08-11** — **Stage 1 complete, gate green.** `src/Mesh2D.jl`:
  - **Delaunay** (incremental Bowyer–Watson) closed by a single GHOST vertex at
    infinity — *not* a finite super-triangle, which is not robust (a point outside
    the evolving hull can be wrongly excluded from a super-adjacent circumcircle,
    cracking the boundary; verified failure on ~half of random seeds). Ghost
    in-circle = orientation test on the real hull edge; collinear-aware so
    collinear hull points never spawn flat triangles. Remembering stochastic walk.
  - **CDT**: flip-based (Sloan) constrained segment insertion with `_flip!`
    (verified: double-flip = identity), through-vertex splitting, Lawson
    legalization; parity flood-fill `classify_interior` (nested holes correct).
  - **Ruppert refinement**: encroached-subsegment splitting + skinny-triangle
    circumcenter insertion (constrained), achieving the min-angle/max-area bound
    while preserving the domain and constraints.
  - Bug caught by `--check-bounds=yes` (Pkg.test): a Steiner-insertion OOB read
    (`_pt(T, nreal+1)`) masked by `@inbounds` in normal mode; fixed. Lesson:
    verify under bounds-checking. Determinism bug (Set iteration order in
    encroachment scans) fixed by sorted iteration.
  - Oracles: exact empty-circumcircle, 2n−2−h + hull cross-check, constrained-
    Delaunay locally-Delaunay check, min-angle/area quality, domain-area
    conservation. Full `Pkg.test()` green: **142,240 assertions**. `stage()` → 2.
- **2026-08-11** — **Stage 2 complete, gate green.**
  - `src/SizeField.jl` — `ConstantSize`/`FunctionSize`/`MinSize`, `size_at`.
  - `src/Mesh1D.jl` — graded 1-D edge meshing: nodes at equal increments of the
    metric length ∫|γ'|/h. Verified vs analytic arc length (segment/circle/helix),
    uniform + graded spacing, closed loops (no duplicate node).
  - `src/MeshSurface.jl` — planar faces (Newell frame → project → CDT + size-
    driven Ruppert → lift; exact area, quality preserved), cylinder lateral
    surface (**graded structured**, watertight by construction — replaced a
    Ruppert-unrolled variant whose seam gapped under a varying size field, a bug
    caught by a seam-gap-edge count), and parametric patches (isotropic-metric
    approximation; area converges to analytic).
  - `refine!` extended with a `size` callback (refine where an edge exceeds the
    local target). Full `Pkg.test()` green: **142,282 assertions**. `stage()` → 3.
- **2026-08-11** — **Stage 3 core (3-D kernel + volume filling).** `src/Mesh3D.jl`:
  - Incremental 3-D Delaunay via ghost tetrahedra on exact `orient3_sos`/
    `insphere_sos`; correct hull extension; face-orientation invariant proved and
    checked. Deterministic ~1e-8 symbolic perturbation breaks coplanar/cospherical
    degeneracies (the maximally-degenerate unit cube tetrahedralizes cleanly).
    Bugs fixed via oracles: ghost in-sphere sign, new-tet face orientation, and a
    walk-termination gap (brute-force `_locate_scan` fallback).
  - `tetrahedralize` (fill from a closed surface, ray-cast point-in-polyhedron)
    and `tetrahedralize_multi` (per-region, tagged). Verified fills at *exact*
    volume: convex cube/octahedron, non-convex L-prism, genus-1 box-with-tunnel,
    thin cylindrical pin (aspect 20), and a **3-region coax junction with all
    volumes filled** — the multi-way junction gmsh leaves empty.
  - Honest scope: robust conforming *interface* recovery (a naive Steiner version
    diverged under perturbation and was removed) and evaluating the literal
    enclosure `.geo` (needs the Stage-5 geometry kernel) remain. Full `Pkg.test()`
    green: **142,335 assertions**.
- **2026-08-11** — **Stages 4–5 partials + integration.**
  - `src/Optimize.jl` — `mesh_quality` (dihedral/radius-edge/sliver report) and
    boundary-preserving `smooth_laplacian` (positive-volume guard; volume/validity/
    tags preserved; mean dihedral up, slivers down). `refine3d!` was prototyped and
    **removed** (unconstrained circumcentre insertion raised mean circumradius by
    creating slivers — CRC: don't ship a tool that degrades its output).
  - `src/Heal.jl` — `surface_diagnostics`/`is_meshable`: detect open, non-manifold,
    degenerate, duplicate, mis-oriented, and near-coincident defects (the detection
    half of "heal, don't fail").
  - `src/Geometry.jl` — native `box_surface`/`cylinder_surface`/`box_tunnel_surface`
    primitives (closed, manifold, fill to exact analytic volume).
  - Top-level `Tessella.mesh_volume` — surface→volume pipeline with the
    validated-or-explicit-blocker contract (Heal gate → fill → smooth → validate).
  - Deep-debug pass on the 3-D+ modules: default path (`perturb=true`) robust; the
    `perturb=false` coplanar case is the known SoS-artifact (perturbation is the
    fix), documented. Full `Pkg.test()` green: **142,388 assertions**.
- **2026-08-11** — **Stage 6 partial (high-order) + integration.** `src/HighOrder.jl`:
  quadratic (P2) tet generation (`p2_tetmesh`, one shared node per edge, straight
  P2 volume == linear), gmsh **type-11** writer + reader (round-trips to the same
  linear connectivity CRC), and `curve_to_cylinder!` (projects lateral-wall edge
  nodes onto the exact cylinder — genuine curved high-order for a primitive,
  verified to land at radius R within 1e-12). End-to-end integration test:
  primitive → `mesh_volume` → `.msh` (v2/v4) → read reproduces the CRC; multi-
  region tags + physical names survive. All 7 stages now have working, CRC-gated
  implementations; the remaining work is the research-hard / out-of-near-term-scope
  tail (conforming interface recovery, sliver exudation, Boolean CSG + OCC interop,
  general-geometry curving, ASCENT drop-in + 22-case regression). Full `Pkg.test()`
  green: **142,414 assertions**.
