# Tessella.jl continuation

Read [HANDOFF.md](HANDOFF.md#takeover-2026-10-11), [STATUS.md](STATUS.md),
[PLAN.md](PLAN.md), and [test/README.md](test/README.md) before continuing.
The 2026-10-11 user request replaces older exhaustive test requirements:
default tests must be concise functional checks on small real fixtures.
Native parity sweeps, source audits, stress, benchmarks, HFSS and external
ASCENT solves are explicit opt-in work.

Tessella remains an independent Julia mesher. Do not use Gmsh as its production
backend or read restricted matching implementation code. Support Julia 1.12.x
and 1.13.x as specified by Project.toml. Public operations must produce valid
results or precise errors, preserve ownership and caller state, and never hide
an empty/invalid region by skipping a geometric certificate.

Continue actual implementation tasks, diagnose confirmed bugs, and verify each
change with proportionate public calls and independent invariants. Preserve
historical fixtures/receipts and unrelated changes. Use subagents for independent
work and coordinate file ownership. Normally push verified implementations to
`main`; do not force-push. A documentation-only commit does not publish code.

Full P1–P6 parity remains unfinished. The current handoff records integrated and
pending packets and actual runtime outcomes. Do not claim completion from
prepared code, scoped checks, historical receipts, or a fast default suite.
When all authorized work is complete, stop.
