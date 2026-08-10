# Tessella.jl — Development Status Tracker

**Single source of truth for build status. Update every session. Nothing dropped.**

Package: `Tessella` · Julia ≥ 1.11 · goal: robust Julia-native mesh generator
replacing the gmsh dependency (see `PLAN.md`). CRC discipline mandatory
(`DEVELOPMENT.md`).

## Stage board

| Stage | Scope | State | CRC gate |
|---|---|---|---|
| 0 | Foundations: repo, CI, mesh types, `.msh` I/O, exact predicates | **NOT STARTED** | predicates vs exact-rational oracle; `.msh` round-trip vs gmsh |
| 1 | 2-D Delaunay + CDT + quality refinement | not started | angle/area guarantees; CRC vs analytic |
| 2 | 1-D edge meshing + parametric surface meshing | not started | HFSS flat/patch/coax faces |
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

## Log

- **2026-08-11** — repo scaffolded (PLAN, STATUS, DEVELOPMENT, README, skeleton,
  startup). gmsh 4.13.1 source studied and archived as reference. Package name
  `Tessella.jl` chosen. Local git initialized; GitHub remote pending `gh auth
  login`. **No code implemented yet** — Stage 0 is the first development task
  (start from `startup.md`).
