# Tessella.jl session handoff

The published Stage 0–6 simplex baseline has recorded verification; full Gmsh
4.15.2 parity remains incomplete. Start the next session by reading, in order:

1. [`HANDOFF.md`, continuation 2026-10-10](HANDOFF.md#continuation-handoff-2026-10-10)
   — preserved candidates, exact receipts, live reader adoption, and unfinished work;
2. [`STATUS.md`](STATUS.md) — which source snapshot each verified gate covers;
3. [`PLAN.md`](PLAN.md) — architecture and the complete parity scope;
4. [`DEVELOPMENT.md`](DEVELOPMENT.md) — mandatory CRC and publication discipline;
5. [`ASCENT.md`](ASCENT.md) for the separate, historical external solver/HFSS work.

Root V24 and release V6 are **unqualified and unpublished candidates** at this
handoff. Adopt the existing Auto reader by its recorded process identity before
starting another package/facade run. Preserve all frozen sources and old red
receipts; follow the handoff's ownership and fresh-memory guards. This
documentation update did not run tests or publish an implementation.

## Project contract

Tessella is an independent Julia mesh generator pursuing Gmsh 4.15.2 feature and
behavioral parity. It must not delegate production meshing back to Gmsh: the ASCENT
acceptance geometry is precisely a Gmsh failure. Implement meshing capabilities first,
then the remaining CAD/API/UI/post-processing tracks in `PLAN.md`. Supported
development and verification runtimes are Julia **1.12.x and 1.13.x**, as specified
by `Project.toml`; qualify the same final source on both.

The verified baseline owns exact predicates, simplex meshing, conforming recovery,
sizing, quality improvement, native analytical/CSG geometry, P2 elements, and
solver-consumable mesh I/O. Do not confuse that baseline with parity completion.

Every public operation must return a validated result or an explicit diagnostic.
Count actual cells per region; an empty list of regions is not proof that no region is
empty. Never weaken a geometric or topology check to make a fixture pass.

## Before changing code

- Identify the exact published or restored source snapshot first. Reproduce on
  that snapshot; a historical gate or a prepared candidate is not current proof.
- Identify an independent oracle or invariant that will fail before the fix.
- Preserve unrelated working-tree changes and the full parity scope in `PLAN.md`.
- Treat the enclosure/coax fixture as a standing regression, not as a one-off special
  case.

## Required closure gate

```sh
julia --project --startup-file=no --history-file=no --check-bounds=yes -O2 --threads=1 --gcthreads=1 -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no","--history-file=no","--check-bounds=yes","-O2","--threads=1","--gcthreads=1"])'
```

Run the package gate with each supported runtime, after the handoff's sole-reader
and memory guards permit it. Preserve explicit four-thread resource controls;
scoped Source/AST/parser results do not replace whole public/package checks.

For I/O or cross-tool changes, also run with both runtimes and the pinned oracle:

```sh
julia --project --startup-file=no --history-file=no --check-bounds=yes -O2 --threads=1 --gcthreads=1 validation/run_all.jl
```

Update `STATUS.md` only with commands and measurements actually run, exact source
hashes, actual process completion/drain, and the limits of each result. Finish all
required combined gates before committing and normally pushing accepted
implementations to `main`. A documentation-only push does not satisfy the pending
implementation push. Preserve all original native/assertion/allocation controls;
do not waive historical source pins without qualifying the changed methods.

The external ASCENT full-wave campaign has its own authority, data, and verification
requirements. Do not fabricate solver/reference results or treat that campaign as an
unimplemented Tessella package feature.
