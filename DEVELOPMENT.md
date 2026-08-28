# Tessella.jl — Development & CRC Discipline

This project meshes geometry that a FEM solver trusts. **A mesh is evidence, not
scaffolding.** These rules are mandatory and mirror the ASCENT research-code
standard.

Development and verification use Julia 1.12.x only. `Project.toml` owns the
machine-readable runtime requirement.

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

- `julia --project -e 'using Pkg; Pkg.test()'` green.
- Stage regression meshes re-checksum-match.
- The standing acceptance cases (`STATUS.md`) for the current stage pass.
- Benchmarks recorded (nodes/s, memory/node) — deterministic cost, not one-shot
  wall-clock noise.

## Reproducibility

Explicit RNG objects, pinned deps, rerunnable fixtures, documented tolerances
with justification comments. Non-deterministic results are unverified until the
variability is quantified.

## Git & CRC provenance

Small, reviewed commits. Each commit that changes a meshing kernel notes the
oracle it was verified against and the regression checksum delta. `main` stays
green.
