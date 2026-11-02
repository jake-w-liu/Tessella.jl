# Repository continuation instructions

Read [the current development index](README.md), [CHECKLIST.md](CHECKLIST.md),
and [test/README.md](../../test/README.md) before implementation. Consult the
[preserved complete ledger](history/DEVIN_HANDOFF.md), latest dated
[handoff](history/HANDOFF.md), [status](history/STATUS.md), and
[plan](history/PLAN.md) for requirement details and prior evidence.

Resume from the current `main` checkout. Older temporary candidates and dated
receipts are historical references. Preserve unrelated user changes and
untracked files. Normally push verified implementations to `main`; never
force-push. Continue genuine unfinished features after each verified increment.
Do not declare full parity complete from a bounded suite.

Complete the full P1–P6 scope and newly confirmed requirements before declaring
Tessella complete. Keep the current checklist and verification evidence accurate.
The prior Devin transfer is historical context, not a reason to leave work undone.

Routine tests must be concise functional checks on small real fixtures. The
2026-10-11 user instruction supersedes older exhaustive publication gates.
Native-oracle sweeps, source/AST audits, allocation stress, benchmarks, HFSS,
and external solver campaigns are opt-in. Use independent analytic or topology
checks where they establish correctness. Record actual exits and distinguish
minimal compilation from optimized compiler checks. Run ordinary `Pkg.test()`
on both supported Julia lines for changed final source.

Production meshing stays independent of Gmsh. Do not read or port restricted
Gmsh/Concorde matching implementations. Public SDK observations and retained
opaque fixtures may be used for comparison. Preserve ownership and certificates;
fix missing producer transport rather than bypassing validation.

Coordinate subagent file and runtime ownership. Identify a process before
stopping it and leave unrelated jobs alone. Stop when all authorized work is
complete.
