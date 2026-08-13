# Tessella.jl session handoff

Tessella's package roadmap is complete. Start a maintenance or development session by
reading, in order:

1. [`PLAN.md`](PLAN.md) — current architecture and scope boundary;
2. [`DEVELOPMENT.md`](DEVELOPMENT.md) — mandatory CRC discipline;
3. [`STATUS.md`](STATUS.md) — latest verified package gate;
4. [`ASCENT.md`](ASCENT.md) only for external solver/HFSS work.

## Project contract

Tessella is an independent Julia mesh generator, not a literal gmsh port and not an
OpenCASCADE/NURBS reimplementation. It owns exact predicates, simplex meshing,
conforming recovery, sizing, quality improvement, native analytical/CSG geometry,
P2 elements, and solver-consumable mesh I/O.

Every public operation must return a validated result or an explicit diagnostic.
Count actual cells per region; an empty list of regions is not proof that no region is
empty. Never weaken a geometric or topology check to make a fixture pass.

## Before changing code

- Reproduce the issue against the current `main`; historical notes are not proof.
- Identify an independent oracle or invariant that will fail before the fix.
- Preserve unrelated working-tree changes and the explicit scope in `PLAN.md`.
- Treat the enclosure/coax fixture as a standing regression, not as a one-off special
  case.

## Required closure gate

```sh
julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

For I/O or cross-tool changes, also run:

```sh
julia --project --check-bounds=yes validation/run_all.jl
```

Update `STATUS.md` only with commands and measurements actually run. Keep `main`
green, commit fixes in oracle-justified units, and push each accepted fix to `main`.

The external ASCENT full-wave campaign has its own authority, data, and verification
requirements. Do not fabricate solver/reference results or treat that campaign as an
unimplemented Tessella package feature.
