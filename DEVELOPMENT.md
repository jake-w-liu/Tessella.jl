# Tessella.jl — Development & CRC Discipline

This project meshes geometry that a FEM solver trusts. **A mesh is evidence, not
scaffolding.** These rules are mandatory and mirror the ASCENT research-code
standard.

Development and verification use Julia 1.12.x and 1.13.x. `Project.toml` owns the
machine-readable runtime requirement.

At the 2026-10-10 handoff, full Gmsh 4.15.2 parity is unfinished. Root V24 and
release V6 are unqualified integration candidates. Read the
[current continuation handoff](HANDOFF.md#continuation-handoff-2026-10-10)
for exact source qualifications, live reader ownership, remaining gates, and the
required normal push of verified implementations to `main`. Historical checks
do not qualify a changed candidate.

## Core loop (every change)

`spec → independent oracle → implement smallest slice → CRC test → reverify →
fix at source → reverify`. No task is done after implementation alone.

## CRC (Correctness–Robustness–Completeness) gate

Every nontrivial function ships with **all three**:

1. **Correctness** — an *independent oracle*, never the code checking itself:
   - exact-rational recomputation (predicates),
   - analytic mesh (unit cube, sphere, known Delaunay triangulation),
   - invariants: Euler characteristic, Delaunay empty-sphere, positive volumes,
     manifold/closed boundary, boundary-facet conservation,
   - cross-check against a reference gmsh mesh (counts, bbox, quality, boundary
     hash) where a reference is legitimate.
2. **Robustness** — realistic and degenerate inputs: cospherical/coplanar points,
   near-coincident faces, slivers, thin features, multi-way junctions. Fixed RNG
   seeds. No happy-path-only tests.
3. **Completeness** — real error handling and a **validated or explicit-blocker**
   contract: return a mesh that passes validation, or a precise diagnostic. Never
   a silently empty region.

## Mesh-CRC checksum

Each accepted mesh emits a deterministic checksum recorded in `STATUS.md` /
`test/artifacts/`:
`(n_nodes, n_edges, n_tris, n_tets, bbox, min/mean dihedral, min/mean radius-edge,
 boundary-facet count, SHA-256 of sorted connectivity)`.
Regression = re-derive the checksum and diff. A change requires a justified note.

## Anti-false-positive rules (hard-won)

- **Count the actual elements.** "No empty volumes" is meaningless if the volume
  list is empty. Report tets-per-region; a valid volume mesh has tets in *every*
  region. (This exact trap sank `occ.healShapes()` in the ASCENT campaign.)
- Suspicious success, suspicious failure, flat/degenerate output, zero-count
  regions, NaN/Inf, negative volumes, or too-good quality = a bug until an
  independent oracle disproves it.
- Never weaken a tolerance, delete a check, or change an expected value without
  first proving the prior expectation wrong.
- No `@test true`, `x == x`, or self-`isapprox`.

## Predicate policy (foundational)

3-D meshing robustness *is* predicate robustness. `orient2/orient3/incircle/
insphere` are **adaptive exact** (Shewchuk-style staged precision) with a
**Simulation of Simplicity** tie-break. They get an exhaustive degenerate test
against exact rationals before any mesher uses them. Non-negotiable.

## Test & harness gate (before "done")

- Before claiming closure, the complete package gate must be green on
  **both Julia 1.12 and Julia 1.13**,
  using the same final source bytes. The handoff records exact tested versions
  (1.12.7 and 1.13.1); `Project.toml` owns the supported range.
- Record the effective child-runtime flags as well as the parent command.
  A closure command, run with each supported runtime, is:

  ```sh
  julia --project --startup-file=no --history-file=no --check-bounds=yes -O2 --threads=1 --gcthreads=1 -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no","--history-file=no","--check-bounds=yes","-O2","--threads=1","--gcthreads=1"])'
  ```

  Preserve the original explicit four-thread resource tests. Scoped Source,
  extracted-AST, parser, or native-oracle checks retain their stated scope and
  do not replace the complete package/public gates.
- Stage regression meshes re-checksum-match.
- The standing acceptance cases (`STATUS.md`) for the current stage pass.
- Benchmarks recorded (nodes/s, memory/node) — deterministic cost, not one-shot
  wall-clock noise.

For the current parity continuation, first adopt and drain the existing reader
listed in `HANDOFF.md`. Use the bound runners, fresh physical and available
commit memory above 30 GiB, normal `-O2`, default inlining, bounds checks, and
full pre/post source, fixture, runtime, and archive hashes. One package/facade
reader owns the slot at a time. A prepared or scoped green result does not
authorize a second reader or establish combined release qualification.

Keep literal native rows, bit patterns, diagnostics, counters, state assertions,
and original allocation/growth checks. Diagnose a failed source/fixture hash
against its actual byte lineage; preserve exact EOL transport where pinned.
Changed method bodies need new qualification before their provenance pins can
change. A confirmed harness correction needs bound before/after evidence.

## Reproducibility

Explicit RNG objects, pinned deps, rerunnable fixtures, documented tolerances
with justification comments. Non-deterministic results are unverified until the
variability is quantified.

## Git & CRC provenance

Small, reviewed commits. Each commit that changes a meshing kernel notes the
oracle it was verified against and the regression checksum delta. Record the
exact qualified source and runtime receipts; preserve unrelated changes.
After the required combined gates pass, commit and **normally push the verified
implementations to `main`**. Documentation commits and handoffs do not publish
an uncommitted implementation or complete the parity goal. Do not force-push
or describe Root V24/release V6 as shipped before integration and qualification.
