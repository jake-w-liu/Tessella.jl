# Tessella.jl — Development & CRC Discipline

This project meshes geometry that a FEM solver trusts. **A mesh is evidence, not
scaffolding.** These rules are mandatory and mirror the ASCENT research-code
standard.

Development and verification use Julia 1.12.x and 1.13.x. `Project.toml` owns the
machine-readable runtime requirement.

The 2026-10-11 Main integration includes Devin's implementation and reviewed
curved/periodic/matching/view packets. Full Gmsh parity remains unfinished.
Read the [current handoff](HANDOFF.md#takeover-2026-10-11) for actual integration
and verification status. Use the concise default suite described below;
historical exhaustive campaign requirements do not override the user's current
testing instruction. Normally push verified implementations to `main`.

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
**Simulation of Simplicity** tie-break. Their historical exhaustive
exact-rational controls remain available in the opt-in suite. Routine tests use
small exact/degenerate representatives; predicate changes need focused independent
oracles appropriate to the changed arithmetic.

## Test & harness gate (before "done")

The default package suite must remain concise: functionality and correctness
through small real public calls, analytic expectations, topology, ownership,
state transitions and meaningful invalid inputs. Do not add production-size
runs, exhaustive native sweeps, source-layout/hash/AST checks, allocation
campaigns, benchmarks or solver workloads to `test/runtests.jl`.

Run the default suite on **both Julia 1.12 and Julia 1.13** using the same final
source. Record the runtime, effective flags and measured duration:

```sh
julia --project=. --startup-file=no --compile=min -O0 -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes", "--compile=min", "-O0"])'
```

Use this compilation mode for routine functional checks. For compiler-focused
checks, run the same assertions with normal `-O2`, omitting `--compile=min` from
both Julia and `julia_args`. Keep functional and compiler timing receipts clear
about which mode was measured.

Run additional focused regressions relevant to a change when needed. The
historical comprehensive runner is explicitly opt-in:

```sh
julia --project=. --startup-file=no --check-bounds=yes test/extended/runtests.jl
```

External mesher differentials and application, stress, performance and solver
campaigns remain under `validation/` or the extended regression files. They are
separate evidence and may take hours. Keep their original assertions, fixtures,
allocation checks and archived receipts. Run these campaigns when validating
their claimed scope; an ordinary package check does not require every campaign.

Historical release/parity instructions in `HANDOFF.md` describe those campaigns
and their evidence. They do not override the concise default suite. Any claim
of parity or production qualification still needs evidence appropriate to that
claim. See [`test/README.md`](test/README.md) for suite boundaries and commands.

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
