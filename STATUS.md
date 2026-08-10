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
| 2 | 1-D edge meshing + parametric surface meshing | in progress | HFSS flat/patch/coax faces |
| 3 | 3-D Delaunay + **robust boundary recovery** + slivers | not started | **mesh enclosure coax junction (9.2)** |
| 4 | size fields + optimization | not started | quality ≥ gmsh on 22 cases |
| 5 | geometry kernel (OCC interop / native CSG) + heal | not started | ingest ASCENT `solid_model` |
| 6 | high-order + ASCENT integration + 22-case regression | not started | all HFSS cases solve via Tessella |

## Standing acceptance cases (regression, CRC-stamped when they pass)

| id | geometry | why | status |
|---|---|---|---|
| ENC-COAX | enclosure coax feed-through (HFSS 9.2), surface 86 shield/case/air junction | gmsh 4.13.1 **and** 4.15.2 fail (`overlapping facets`, all volumes empty); the project's reason to exist | **OPEN — reference input captured** |
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
