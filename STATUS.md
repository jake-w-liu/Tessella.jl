# Tessella.jl — Development Status Tracker

**Single source of truth for build status. Update every session. Nothing dropped.**

Package: `Tessella` · Julia ≥ 1.11 · goal: robust Julia-native mesh generator
replacing the gmsh dependency (see `PLAN.md`). CRC discipline mandatory
(`DEVELOPMENT.md`).

## Current state (verified at HEAD)

- **Suite:** `julia --project=. -e 'using Pkg; Pkg.test()'` green — **142,952
  assertions** under `--check-bounds=yes` (verified this session).
- **CRC regression checksums preserved** (suite includes them; unaffected by the
  `orient2_sos` fix): unit-cube
  `7ea403054f05392f18b404a1f5f78b12d70d45d40c7b04ba8f8dc3e030d8f3f9`; 10×10
  refined square `583c615df1862c8518bbda409347f109dc25f7f8f5362562badf160fe6af30c1`.
- **All 7 stages carry working, CRC-gated code.** Stages 0–2 complete; Stage 3
  kernel + volume filling complete; Stages 4–6 partial (see board). Public API:
  `mesh_volume` (3-D, `optimize`/`smooth` opts), `mesh_planar` (2-D), plus the
  per-module functions.
- **Reason-to-exist demonstrated:** the ENC-COAX coax feed-through is meshed as ONE
  **conforming** partition — all three volumes (air/case/pin) filled with a manifold,
  shared interface (the pin bore through the case wall conforms), at the **literal
  `.geo` scale and real pin resolution** — exactly where gmsh 4.13/4.15 produce zero
  volume tets. Enabled by the corrected SoS (exact coordinates, no jitter).
- **Adversarial-audit + feature pass landed this session** (workflows + independent
  re-verification): **9 confirmed bugs fixed** — including all four exact-predicate
  SoS tie-breaks (`orient2/orient3/incircle/insphere`) now consistent under one +ε
  scheme (regression-pinned) — plus 2 features integrated and a `validation/`
  cross-check suite vs gmsh added (see "Audit findings" below).

## Stage board

| Stage | Scope | State | CRC gate |
|---|---|---|---|
| 0 | Foundations: repo, CI, mesh types, `.msh` I/O, exact predicates | **DONE — gate green** | predicates vs exact-rational oracle ✓; `.msh` round-trip (v2↔v4) CRC-preserving ✓ |
| 1 | 2-D Delaunay + CDT + quality refinement | **DONE — gate green** | exact empty-circumcircle oracle ✓; 2n−2−h count + hull ✓; CDT constraints present + locally Delaunay ✓; Ruppert angle/area bound achieved, domain area preserved ✓ |
| 2 | 1-D edge meshing + parametric surface meshing | **DONE — gate green** | size-graded edge mesh vs analytic arc length ✓; planar face exact area + quality ✓; cylinder watertight + area convergence ✓; parametric area convergence ✓ |
| 3 | 3-D Delaunay + **robust boundary recovery** + slivers | **kernel + filling DONE; conforming partitions (planar + curved-convex interfaces) DONE; general recovery + slivers WIP** | 3-D empty-circumsphere oracle ✓; convex/non-convex/genus-1/thin/multi-region fills validated ✓; representative coax junction all volumes filled ✓; **`tetrahedralize_conforming`** — exact-coordinate Delaunay (no jitter) + per-region ray-cast tagging → interfaces are *shared* Delaunay faces (verified manifold + exact volume + every region filled on the **3-region pin/air/case enclosure topology**, incl. the pin's *curved N-gon* interface — all 4·nθ faces shared, no Steiner points) — incl. a pin **feed-through crossing the case wall** (bore handled by classification, both pin↔air and pin↔case interfaces conform, no CSG) ✓, **enabled by the corrected SoS** (jitter of `perturb=true` moved interface vertices off-plane, destroying conformance); general constrained recovery for interfaces *not* expressible as Delaunay faces (arbitrary non-convex regions, or non-Delaunay input surface triangulations) OPEN |
| 4 | size fields + optimization | **partial** | tet quality report ✓; Laplacian + ODT smoothing ✓; verified 2-3/3-2 flip primitives ✓; `optimize_flips!` sliver reduction (309→169 flips alone, →63 with smoothing on a 300-pt cloud; volume/validity/min-dihedral-safe) ✓; weighted sliver exudation + domain-bounded 3-D refinement WIP |
| 5 | geometry kernel (OCC interop / native CSG) + heal | **partial** | `Heal` surface-defect detection ✓; native primitives (box/cylinder/box-tunnel/**hollow-box** `box_shell_surface`, an axis-aligned Boolean difference — verified oriented via divergence theorem + fills to exact `outer−inner`) ✓; general Boolean CSG + OCC interop WIP |
| 6 | high-order + ASCENT integration + 22-case regression | **partial** | quadratic (P2) tet generation ✓ (shared mid-nodes, exact-midpoint volume) + gmsh type-11 I/O ✓; **general-geometry curving** ✓ — `curve_to_surface!`(project, on_surface) / `curve_to_cylinder!` curve only genuine *boundary-surface* edges (interior chords left straight) and revert any projection that would invert an incident P2 element (`p2_min_jacobian` guard over the degree-3 nodes; no inverted element ever emitted, verified by an independent degree-6 sampler); **the conforming ENC-COAX mesh is solver-consumable** — written to gmsh MSH v4.1 with the three physical volumes (`coax_pin`/`air`/`case`) and round-tripped (connectivity CRC + region tags + physical-group names preserved, `pipeline_test.jl`), i.e. directly readable by ASCENT ✓; full ASCENT drop-in + 22-case regression WIP |

## Standing acceptance cases (regression, CRC-stamped when they pass)

| id | geometry | why | status |
|---|---|---|---|
| ENC-COAX | enclosure coax feed-through (HFSS 9.2), surface 86 shield/case/air junction | gmsh 4.13.1 **and** 4.15.2 fail (`overlapping facets`, all volumes empty); the project's reason to exist | **PARTIAL (real-scale volumes filled)** — at the *literal* `.geo` dimensions the three regions (air / case / coax pin, radius 0.8 mm × length 159 mm, aspect ≈ 199) are each meshed and **all volumes filled & validated** together (`mesh3d_test.jl`), the exact place gmsh leaves them empty. **The conforming pin/air/case *topology* now works** (`tetrahedralize_conforming`, `mesh3d_test.jl`): a coax **pin (cylinder) inside an air cavity inside a metal case shell**, meshed as ONE mesh — every region filled, exact total volume, and a **manifold shared interface** including the pin's *curved* (N-gon) surface, which the exact-coordinate Delaunay recovers **without Steiner points** (all 4·nθ pin faces shared pin↔air; air↔case cavity shared). The corrected SoS removes the perturbation-vs-conformance barrier. **The feed-through also conforms**: with the pin **crossing the case wall**, classification `[pin,air,case]` tags every in-cylinder tet as pin regardless of cavity-vs-wall, so the cylindrical **bore is handled with no explicit CSG surface** — both the pin↔air (cavity) *and* pin↔case (bore) interfaces appear and conform, and the air↔case boundary conforms around the bore hole (all 2·nθ·nz pin faces shared, 0 non-manifold, exact total volume). **This holds at the literal `.geo` scale and the real pin resolution**: air/case/pin at the fixture dimensions (0.8 mm pin, aspect ≈199) with nθ=12·nz=40 → 2549 tets, **all three volumes filled**, feed-through conforming (pin↔air 144, pin↔case 816 faces), 0 non-manifold, valid — the exact 3-volume mesh gmsh cannot produce. Remaining: (a) **performance** at fine resolution — for the real nθ=12·nz=40 pin (498 verts) the run is ~100 s, of which the **`perturb=false` Delaunay is ~84 s** (measured — pathological for 498 pts). *Root cause found by measurement, not assumed:* it is **not** the exact predicates (making the exact fallback 1.84× faster via integer arithmetic left the Delaunay time unchanged), so it is the **point-location behaviour on the degenerate coplanar cylinder vertices** under exact coordinates (the walk degrading toward its brute-force `_locate_scan` fallback → ~O(n²)). Fix = a more robust exact walk / spatial acceleration; the classification is already bbox-reject-optimized. *Two simple fixes tried + rejected by measurement:* faster exact predicates (kept — 1.8× — but the Delaunay was unchanged, so not predicate-bound) and a **spatial-sort (Morton) insertion order** (reverted — it made the thin-cylinder case *worse*, 84→142 s: plain spatial order triggers a sliver cascade, the failure mode BRIO's randomized rounds avoid). So the real fix is **BRIO with randomization or a dynamic spatial point-location structure**, a deeper kernel change. (b) the extra literal geometry (slots, resistor, `p1_surface`) beyond the three main volumes. Neither is a boundary-recovery problem, and gmsh produces **zero** tets here at any speed. |
| THIN-SLOT | enclosure 1 mm coupling slot | thin-feature Boolean sliver | open |
| SPIRAL | 10.1 silicon spiral (thin swept traces, layered stack) | high-aspect thin conductors | open |
| ARRAY-PML | 5.7 endfire unit cell (periodic + PML) | conformal periodic faces | open |

Reference artifact for ENC-COAX: the failing gmsh `.geo` is captured at
`test/fixtures/enclosure_coax_junction.geo` (parsed by `IO.read_geo_params`, which
recovers the air/case/coax_pin groups + sizing + seed). gmsh-API diagnosis:
surface 86 bbox ≈ (168.4,140,148.4)–(171.6,160.5,151.6) mm; empty volumes =
pin/case/air.

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

## Audit findings (adversarial workflows, 2026-08-11) — independently re-verified

Two background workflows (a feature-implementation pass + a core-module adversarial
bug-hunt with per-finding verification) surfaced these. Every one was reproduced
here from scratch before acting; every fix re-verified (full suite green).

**Fixed (9):**

| # | site | defect | verification of the fix |
|---|---|---|---|
| 1 | `HighOrder.curve_to_cylinder!` / new `curve_to_surface!` | curved **interior** edges (qualified an edge by "both corners on the surface", but every vertex is on the boundary since no Steiner nodes are added) → tangled 47/125 (cyl) and 127/166 (sphere) P2 elements | now curves only genuine boundary-surface edges + reverts any inverting projection; independent degree-6 Jacobian sampler finds **0** inverted elements; interior chords byte-identical to their straight midpoints |
| 2 | `Predicates.orient2_sos:351` | SoS tie-break sign-inverted (matched −ε) vs the +ε house convention → 2-D degenerate orient2/incircle decisions mutually incoherent (Mesh2D uses both on one point set) | negated to +ε; **0/6** disagreements vs an independent exact `Rational{BigInt}` SoS oracle (was 6/6); full suite + golden 2-D CRC unchanged |
| 7 | `Predicates.orient3_sos:416` | fully-degenerate (4-collinear) fallthrough returned a constant `perm`, wrong on 12/24 labelings | replaced with an exact +ε `Rational{BigInt}` leading-term SoS evaluator (Leibniz det of the perturbed matrix); verified **0/24** collinear + **0/24** coplanar-sanity vs a validated oracle (fast path preserved) |
| 8 | `Predicates.incircle_sos:439` | 4-collinear fallthrough returned constant `perm`, wrong on 12/24 | replaced with the exact SoS of the paraboloid-lifted points (incircle = orientation of the lift; the algorithm-consistent Edelsbrunner–Mücke scheme). Verified **0/24** collinear + **0/24** cocircular-sanity vs a **canonical lifted-orientation oracle** — which I built and validated after finding my first (naive-lift) oracle used the wrong scheme; a finite-ε exact-Rational cross-check confirmed the resolution |
| 9 | `Predicates.insphere_sos:456` | the *entire* tie-break (fast minors **and** fallthrough) evaluated the −ε perturbation while orient3 uses +ε — the 3-D analog of #2: the kernel needs orient3/insphere to agree (cavity chosen by insphere, orientations by orient3). Wrong on 120/120 cospherical + 60/120 coplanar vs +ε (three independent methods) | replaced the whole tie-break with the exact +ε SoS of the 4-D paraboloid lift (`_orient_nd_sos_exact`, the same construction as incircle one dimension up — the d=3 form is oracle-validated, and the base-sign match + dimension-general construction carry it to d=4). Verified **0/120** cospherical + **0/120** coplanar vs the canonical oracle; full suite + golden CRC unchanged |
| 3 | `MeshTypes.euler_characteristic:274` | counted segment-only nodes in V while E omits segment edges → χ inflated (tet+segment read 3, not 1) | segment loop removed; reproduced (3→1) + regression test |
| 4 | `MeshTypes.bounding_box:315` | included unreferenced nodes, violating its "referenced nodes" contract and polluting `mesh_crc.bbox` | filter to referenced nodes; reproduced ((1000,1,1)→(1,1,1)) + regression test |
| 5 | `IO.read_stl:481` | far-from-origin weld key `round(Int, coord/tol)` overflowed Int64 → crash | quantize relative to bbox min corner; reproduced (InexactError→ok) + regression test |
| 6 | `IO.read_stl:431` | ASCII STL with a UTF-8 BOM or non-ASCII solid name misdetected as binary → EOFError crash | BOM-aware keyword check, drop whole-header ASCII requirement; reproduced (2 crashes→ok) + regression test |

**All four SoS predicates now resolve degeneracies by one consistent +ε scheme.**
The exact `Rational{BigInt}` leading-term evaluator (`_orient_nd_sos_exact`, a
Leibniz determinant of the EM-perturbed matrix) is the shared foundation: orient3
uses it directly for the collinear fallthrough; incircle/insphere use it on the
paraboloid-lifted points (incircle = 3-D orientation of the lift, insphere = 4-D).
The corrected +ε signs over every index labeling of the four finding configs are
**regression-pinned** in `test/predicates_test.jl`. Nothing SoS remains deferred.

> Method note: getting incircle/insphere right required distinguishing two SoS
> schemes for lifted predicates — deriving the lift from perturbed x,y (wrong: its
> ε² cross-terms mis-order the leading term) vs. treating the lift as an independent
> EM coordinate (correct, = orientation of the lifted points). My first oracle used
> the former and mis-flagged the code; the canonical lifted-orientation oracle (+ε
> by construction, base-matched to each predicate) resolved it, cross-checked with a
> finite-ε exact-Rational evaluation. The two audit verify-agents had disagreed on
> insphere's sign; the canonical oracle + full-suite/CRC settle it.

## Cross-tool validation vs gmsh (`validation/`, run this session)

`julia --project=. validation/run_all.jl` meshes the same domain with Tessella and
gmsh 4.15.2 and compares meshed volume (vs the analytic oracle), quality, and time.
Reference `.geo` scripts retained per case; regenerate `REPORT.md` any time.

| case | analytic V | Tessella V | gmsh V | verdict |
|---|---|---|---|---|
| box (flat) | 2 | **2.0** | **2.0** | both exact ✓ |
| box_tunnel (genus-1 flat) | 24 | **24.0** | **24.0** | both exact ✓ |
| hollow_box (Boolean difference) | 35 | **35.0** | **35.0** | Tessella `box_shell` CSG matches gmsh's Boolean diff exactly ✓ |
| cylinder (true πR²H=62.83) | — | 62.65 (0.3%) | 61.02 (2.9%) | Tessella's 48-gon prism is **closer** than gmsh at this size |
| sphere (true 4/3πR³=20.58) | — | 20.10 (2.3%) | 19.69 (4.3%) | Tessella **closer** than gmsh at this size |
| **enclosure_coax** (ASCENT) | — | native primitives fill exactly | **0 volume tets** (28969 surf tris, all air/pin/case volumes empty) | gmsh's documented failure, reproduced here |

Flat-solid rows are a hard correctness cross-check (exact analytic volume, both
tools conform). The enclosure row reproduces gmsh's empty-volume failure directly
from the literal fixture (verified: `ntets=0`).

## CRC regression artifacts (stamped when green)

| artifact | checksum | oracle |
|---|---|---|
| Predicates | 142,141 assertions green; every `orient2/3`, `incircle`, `insphere` sign matches the independent exact-rational homogeneous-determinant oracle over exhaustive integer-grid degeneracies (collinear/coplanar/cocircular/cospherical) + fixed-seed random floats | `test/oracles.jl` (generic Laplace determinant, `Rational{BigInt}`) |
| Unit cube (6-tet Kuhn) | `nodes=8 tets=6 bbox=[(0,0,0)→(1,1,1)] dihedral(min,mean)=(0.7854,0.7854) radedge(min,mean)=(0.86603,0.86603) bfaces=12` · SHA-256 `7ea403054f05392f18b404a1f5f78b12d70d45d40c7b04ba8f8dc3e030d8f3f9` | analytic volume 1/6·6=1, Euler χ=1 (ball), boundary χ=2 (sphere), 12 boundary tris |
| 2-D Delaunay (random/grid/cocircular) | exact empty-circumcircle holds (0 violations, `incircle_sos`); triangle count = 2n−2−h vs independent monotone-chain hull; seed-independent (SoS) | `is_delaunay`, hull oracle, Euler χ=1 |
| 10×10 square refined to 20° (max_area 1.0, seed 1) | `nodes=91 tris=148 bbox=[(0,0,0)→(10,10,0)]` · SHA-256 `583c615df1862c8518bbda409347f109dc25f7f8f5362562badf160fe6af30c1` | min angle ≥ 20°, max tri area ≤ 1.0, domain area = 100 preserved, CDT locally Delaunay |
| 3-D Delaunay (random/box/cube) | exact empty-circumsphere holds (0 violations, `insphere_sos`); Σ tet volume = convex box volume; Euler χ=1 (ball) + boundary χ=2 (sphere) + manifold; seed-independent | `is_delaunay3`, box-volume conservation, `boundary_euler` |
| ENC-COAX (real `.geo` scale) | air/case/coax-pin (R 0.8 mm × L 159 mm, aspect ≈199) — `tets_per_region` all > 0, every region validated, `tetrahedralize_multi` | per-region tet count > 0 (anti-false-positive), exact box volumes, thin-pin fill |
| Volume fills (exact) | cube=1, octahedron=4/3, L-prism=3, box-tunnel(genus-1)=24, cylinder N-gon prism, thin pin | analytic volume + `validate` (positive vol, manifold, watertight) |
| 2-3 / 3-2 flips | consistency + total volume preserved; 2-3 ∘ 3-2 = identity | `check_consistency3`, volume sum, round-trip |

## Benchmarks (deterministic min-of-3, `julia -O2`, this machine)

| kernel | n | time | throughput |
|---|---:|---:|---|
| 2-D Delaunay (`delaunay2d`) | 1 000 | 3.2 ms | 314 k nodes/s |
| 2-D Delaunay | 5 000 | 26 ms | 191 k nodes/s |
| 2-D Delaunay | 20 000 | 166 ms | 120 k nodes/s |
| 3-D Delaunay (`delaunay3d`) | 1 000 | 16 ms | 61 k nodes/s |
| 3-D Delaunay | 5 000 | 118 ms | 42 k nodes/s |
| 3-D Delaunay | 20 000 | 579 ms | 35 k nodes/s |

`mesh_volume` (cylinder 48×8 surface) ≈ 2.5 s / 143 MB — dominated by the ray-cast
point-in-polyhedron domain classification, which is `O(n_tets × n_surface_faces)`.
Correctness is not affected; a spatial-acceleration (BVH/grid) rewrite of the
classifier is the obvious throughput win (noted, not yet done). Raw kernel
throughput is competitive; parallel HXT-class speed is explicitly a
post-robustness goal (`PLAN.md` §6).

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
- **2026-08-11** — **ENC-COAX real-scale + kernel flip primitives.**
  - Constructed the enclosure regions at the literal `.geo` dimensions (air box,
    case shell, coax pin R=0.8 mm × L=159 mm, aspect ≈199) and filled **all three
    volumes** non-empty & validated (`tetrahedralize_multi`) — the exact failure
    mode gmsh cannot. Regression-locked in `mesh3d_test.jl`.
  - `Mesh3D`: verified **2-3 / 3-2 flip** primitives (`flip23!`/`flip32!` on a
    general `_rebuild_region!`; consistency + volume preserved, 2-3∘3-2 = identity)
    + `tets_around_edge`, and a **safe** `optimize_flips!` (hill-climbing 2-3 flips,
    global min dihedral non-decreasing; limited effect on Delaunay input — honest).
    These are the foundation for sliver exudation and constrained boundary recovery.
  - Full `Pkg.test()` green: **142,436 assertions**. Every stage of the plan now
    carries working, CRC-gated code; the remaining tail (conforming interface
    recovery, sliver exudation, Boolean CSG + OCC interop, general curving, ASCENT
    drop-in + 22-case regression) is the research-hard / out-of-near-term-scope
    work `PLAN.md` frames as a multi-year effort.
- **2026-08-11** — **Conforming interface recovery: attempted, characterized, deferred.**
  Implemented the *correct* algorithm — 3-D-Ruppert-style **encroachment** refinement
  (split an encroached subsegment/subface, re-insert, Delaunay-maintain) rather than
  the earlier divergent split-missing-faces version. Verified the geometric kernels
  (triangle circumsphere equidistant; edge/face encroachers). On a **flat** interface
  (box faces) it hits a *fundamental* barrier: the perturbation that keeps the 3-D
  Delaunay flat-tet-free tilts the interface, so recovered Steiner points approach
  coplanarity faster than any jitter can separate them (eventually an unlocatable
  degenerate insertion), and regular grids cascade on cocircular squares. This is
  exactly why robust interface recovery needs a **constrained CDT** (exact, SoS-based
  flip recovery, *no* perturbation) — the multi-year TetGen-class effort `PLAN.md`
  defers. Per CRC discipline the non-robust recovery was **removed** (the verified
  2-3/3-2 flip primitives remain as its foundation). `tetrahedralize_multi` continues
  to fill every region (non-conforming interfaces). Full `Pkg.test()`: **142,442**.
- **2026-08-11** — **Working sliver reduction + 2-D/3-D top-level APIs.**
  - `optimize_flips!` gained **3-2 flips** (with `flip32!` validity guards — the
    edge must be interior and the two apexes must straddle the ring triangle; a
    bug where boundary/degenerate edges relinked wrongly was found and fixed). They
    collapse slivers: 309→169 from flips alone, 309→63 with Laplacian smoothing on
    a 300-pt random Delaunay; volume/validity preserved, min dihedral non-decreasing.
  - `mesh_volume(optimize=true)` runs the flips in the pipeline (cylinder fill
    102→73 slivers, volume preserved). `Tessella.mesh_planar` added — the 2-D
    counterpart of `mesh_volume` (CDT + Ruppert + interior → validated 2-D mesh).
  - Empirical recovery characterization (`scratchpad/analyze.jl`): on the perturbed
    two-region box the interface is *mostly* conforming (6/8 faces present); the 2
    missing faces have missing edges and are not single-flip recoverable → they need
    full edge+facet recovery (constrained CDT), confirming the deferral. Full
    `Pkg.test()` green: **142,445 assertions**.
- **2026-08-11** — **Audit + feature pass, exact-predicate SoS completed, and the
  conforming-interface barrier broken.**
  - Two adversarial workflows (feature impl + core bug-hunt), every finding
    independently reproduced + re-verified (see "Audit findings"): **9 bugs fixed**
    (P2 curving inversion; `euler_characteristic`; `bounding_box`; two `read_stl`
    crashes; and **all four SoS tie-breaks** — `orient2/orient3/incircle/insphere`
    now resolve degeneracies by one consistent **+ε** scheme, an exact
    `Rational{BigInt}` leading-term evaluator shared by all, regression-pinned).
  - `box_shell_surface` (hollow-box CSG) + `smooth_odt` integrated, review-hardened.
    `validation/` suite added (Tessella vs gmsh 4.15; enclosure = 0 gmsh volume tets).
  - **The conforming-interface prediction above came true.** With the SoS predicates
    correct and consistent, `delaunay3d(perturb=false)` runs on **exact coordinates**
    (SoS breaks coplanar/cospherical degeneracies symbolically — no jitter). New
    `tetrahedralize_conforming` meshes all region vertices together and ray-cast-tags
    each tet by region: interfaces become **shared** Delaunay faces. Verified on the
    enclosure **air/case pattern** (inner box + surrounding shell): one mesh, exact
    volume, **manifold shared interface** (the inner box's 12 faces), every region
    filled. The barrier was the coordinate jitter, which the corrected SoS removes.
    General recovery for non-Delaunay/curved interfaces (the literal coax pin) stays
    open. Full `Pkg.test()` green: **142,940 assertions** (+ conforming tests).
