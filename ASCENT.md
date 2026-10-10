# ASCENT integration & HFSS solve campaign — status tracker

## Current continuation status (2026-10-10)

The solver measurements and campaign audit below are **historical evidence from
2026-08-11 through 2026-08-13**. They were not rerun during this documentation
review and do not qualify Root V24, release V6, or a new full-wave campaign.
Tessella currently supports Julia 1.12.x and 1.13.x; the older external solver
results retain their original environment and provenance.

Use [`HANDOFF.md`, continuation 2026-10-10](HANDOFF.md#continuation-handoff-2026-10-10)
for current implementation/verification, live reader ownership, and the required
normal push to `main`. The full parity goal remains unfinished; Root V24 and
release V6 remain unqualified at that handoff. No new ASCENT/HFSS result is claimed.

**Scope of this file:** everything about Tessella *feeding* the external ASCENT FEM
solver and the HFSS validation campaign. Tessella's own development (the mesher:
kernels, recovery, CSG, sizing, I/O) is tracked in `STATUS.md`; this file tracks only
the ASCENT-facing work.

ASCENT is an **external** solver. Its historical source and proprietary HFSS data
were recorded at `/Users/jake/EMPIRE/projects/ongoing/2026_066`; those artifacts
are never pushed to this repo. The
2026-08-12 instruction treated the wider ASCENT campaign as a future external
verification step after mesher completion. Continuing it requires the real
ASCENT binary and proprietary HFSS datasets in that project's environment;
their availability on the current machine has not been reverified here.
Tessella's job is to produce solver-consumable meshes; ASCENT's job is to solve on them.

Nothing here is faked — a fabricated solve result would violate `DEVELOPMENT.md`'s CRC
bar exactly as a silent bad mesh would.

## Historical campaign snapshot (2026-08-12/13)

- **Mesh drop-in: VERIFIED.** A Tessella MSH v4.1 loads whole into ASCENT's real parser
  (`GmshDiscreteModel`, GridapGmsh 0.7.4) with material volumes *and* boundary-condition
  surfaces as physical groups.
- **Solve-on: VERIFIED at the physics level.** ASCENT assembles and solves the Maxwell
  FEM system on Tessella meshes, recovers manufactured solutions to ~1e-15, and computes
  a PEC cavity's resonant frequency to within **0.046 %** of the closed-form analytic
  value — a full geometry→mesh→solve→compare-to-reference regression with an independent
  oracle.
- **Flagship case 9.2 (the enclosure gmsh cannot mesh): meshed natively AND solved.**
- **Remaining:** the other 21 HFSS guide cases as a full-wave solve campaign (external
  ASCENT compute; the geometries are already meshed natively by Tessella — see
  `STATUS.md` / `validation/hfss_cases/`).

## Recorded handshake + solve proofs (`validation/`)

| proof | tag | what it shows |
|---|---|---|
| Mesh handshake | `HANDSHAKE_OK` (`validation/ascent_handshake/`) | a Tessella MSH v4.1 (ENC-COAX pin/air/case conforming partition) loads into `GmshDiscreteModel` with all three region volumes as top-dimensional physical groups → `ASCENT.load_mesh` returns a valid `MeshData` |
| BC handshake | `BC_HANDSHAKE_OK` | a Tessella mesh carrying 3 material **volumes** AND boundary-condition **surfaces** (`radiation` on the domain boundary, `coax_pin_pec` on the pin↔air interface, as 2-D physical groups on tagged faces) loads with all groups visible — solver-consumable *with* BCs, not volumes alone |
| Literal ENC-COAX handshake | `LITERAL_HANDSHAKE_OK` 9/9 (`validation/enclosure_literal/`) | the literal enclosure mesh — 4 volumes + 5 tagged BC surfaces — loads whole in ASCENT (all 9 physical groups) |
| Solve step | `ASCENT_SOLVE_STEP_OK` (`validation/ascent_handshake/solve_step.jl`) | ASCENT assembles AND solves the Maxwell FEM system on a Tessella mesh: `load_mesh`→materials→Nedelec H(curl)→`assemble_diffusive_matrix`→`A\b`; a 165-DOF complex operator (complex-symmetric to 9.6e-17, curl-curl stiffness PSD) whose linear solve recovers a known manufactured field to 1.06e-15 |
| Solve regression | `SOLVE_REGRESSION_OK` (`validation/ascent_solve_regression/`) | solve-usability holds robustly across a 4-geometry-class suite |
| Cavity eigenmode | `CAVITY_EIGENMODE_OK` (`validation/ascent_cavity_eigenmode/`) | ASCENT's eigenmode solver computes a PEC cavity's resonant frequency on a Tessella mesh to within **0.046 %** of the closed-form analytic value (249.711 vs 249.827 MHz) — a complete geometry→mesh→solve→compare-to-reference regression with an independent analytic oracle, the exact shape of an HFSS cavity example |
| Case 9.2 solve | `CASE_9_2_OK` (`validation/enclosure_literal/solve_case_9_2.jl`) | the enclosure — the one case gmsh cannot mesh (0 tets) — is meshed natively by Tessella and ASCENT assembles + solves the Maxwell FEM system on that mesh (107 DOF, complex-symmetric, manufactured-solution field recovered to 1.5e-12) |

## Historical 22-case HFSS audit

The HFSS v10 User Guide worked examples (ch. 5–10). Independent audit at
`2026_066/HFSS_MASTER_TRACKER.md` / `HFSS_22CASE_COMPARISON.html` recorded by
2026-08-13: **all 22 cases were solved + audited in the ASCENT project** (all PASS vs. the guide — 9 exact, 8
qualitative, 4 verified-residual, 1 unscored; every ASCENT result traces to a real
solve-output file, none fabricated). Those audited solves used ASCENT's standard
(gmsh-era) mesh pipeline.

**Tessella's role in #12:**
- **Only case 9.2 (the enclosure) is a gmsh-impossible case** — the one that motivates
  Tessella's existence ("the geo-emitter + gmsh **cannot mesh** the coax"). Tessella
  meshes it natively (gmsh: 0 tets) and ASCENT solves it (`CASE_9_2_OK`).
- **The other 21 are NOT gmsh failures** — ASCENT's standard pipeline already meshed
  them. That historical observation does not complete an end-to-end campaign
  using Tessella meshes.
- Tessella nonetheless meshes **all 22 case geometry classes natively** (representative
  geometries; `validation/hfss_cases/`, tracked in `STATUS.md`), so the mesher covers
  the whole guide.

**Remaining external work:** running the literal 22 guide cases end-to-end on Tessella
meshes — build each antenna geometry + ports/sources/radiation-BC/frequency-sweep in
ASCENT + solve + post-process (S-params / gain / far-field) + compare to the guide
figure. Examples: 5.1 UHF probe, 5.2 conical horn, 5.3 probe-fed patch, 10.1
silicon-spiral inductor. The historical campaign recorded guide references in
`HFSS_22CASE_COMPARISON.html` / `hfss/ug.txt` / the guide PDF. Check those artifacts
and the actual solver environment before continuing; no fresh reference-data
availability or full-wave result is claimed here. This external campaign remains
separate from, and does not substitute for, Tessella's current package parity gates.

## Why the enclosure needs Tessella (gmsh-failure diagnosis, carried from the ASCENT campaign 2026-08-11)

- The enclosure meshing failure is **NOT memory** — a no-time-limit run peaked at 3.3 GB.
  It is a geometry/boundary-recovery defect at the coax junction.
- gmsh `occ.healShapes()` is a **false positive** here: it dissolves the OCC volumes to 0
  tets, so "no empty volumes" is vacuous. Always check the *actual* per-volume tet count.
  (A CRC-discipline lesson, encoded in `DEVELOPMENT.md`.)
- 5 gmsh 2-D/3-D algorithms, OCC fix options, `Geometry.ToleranceBoolean`, and geometry
  protrusion all fail to mesh the enclosure → the robustness must come from Tessella's
  own boundary recovery + heal.

## History (ASCENT-integration milestones)

- **2026-08-12** — ASCENT mesh drop-in verified (`HANDSHAKE_OK`): a Tessella MSH v4.1
  (ENC-COAX pin/air/case conforming partition) loads into `GmshDiscreteModel` (GridapGmsh
  0.7.4) with all three region volumes as physical groups. Extended to full BC structure
  (`BC_HANDSHAKE_OK`): volumes AND boundary-condition surfaces load together.
- **2026-08-13** — Literal ENC-COAX loads whole (`LITERAL_HANDSHAKE_OK` 9/9): 4 volumes +
  5 tagged BC surfaces. ASCENT assembles + solves on Tessella meshes (`ASCENT_SOLVE_STEP_OK`,
  `SOLVE_REGRESSION_OK`); PEC cavity eigenmode to 0.046 % of analytic (`CAVITY_EIGENMODE_OK`);
  flagship case 9.2 meshed natively + solved (`CASE_9_2_OK`).
- **2026-08-13** — Confirmed the 22-case audit is already complete in the ASCENT project
  (all PASS vs. guide); only 9.2 is a gmsh failure. Tessella meshes all 22 geometry classes
  natively (`validation/hfss_cases/`, tracked in `STATUS.md`). The remaining 21-case
  full-wave solve on Tessella meshes is the external ASCENT campaign.
