# Weighted matching source provenance

The Julia weighted Blossom solver adapts the Galil/Edmonds primal-dual algorithm
and bookkeeping described by NetworkX 3.6.1 `max_weight_matching`.

- Source: https://github.com/networkx/networkx/blob/networkx-3.6.1/networkx/algorithms/matching.py
- Raw source SHA256: `5e08becf2f42da936f1b04e2c53155eeeb7b5be97cb3bcd480985af6ff03d2cf`.
- The official tagged source bytes exactly match the installed NetworkX 3.6.1 wheel.
- License: BSD-3-Clause; complete upstream text retained in `networkx-BSD-3-Clause.txt`.
- The older unpromoted Tessella Julia draft was audited against that source before
  adapting indexed state, bounded reusable blossom slots and workspaces.

Minimum-cost graph interpretation for Gmsh 4.15.2 is a separate layer, based on
`meshGFaceOptimize.cpp`; no Gmsh or NetworkX backend executes in production.
Minimum-weight correctness and the pinned primary solver own tie choices are
separate requirements and must be tested independently.

`OrderedMatching.jl` is an independent singleton-dual and alternating-tree
implementation for bipartite graphs. Its initialization and traversal policy are
checked against observed Gmsh 4.15.2 directed pair identities. Integer dual
feasibility, matched-edge tightness, and equality of primal and dual objectives
certify each perfect result. `OrderedFractionalMatching.jl` independently
extends singleton searches to degree-one half-integral solutions: matched
edges and disjoint odd cycles of half-edges. Feasible singleton duals and
equality with the fractional primal certify the relaxation. Fully integral
completions also certify the perfect optimum and retain the input-order
policy. Remaining half-cycles use the BSD weighted Blossom kernel; native
ties in that category are not established by minimum-cost correctness.

No Concorde or Blossom IV implementation is adapted into production. Gmsh's
`contrib/blossom/README.txt` grants its special license only within the Gmsh
system; old unpromoted literal translations remain excluded from this increment.

## Seeded continuation (2026-10-11)

`SeededWeightedMatching.jl` continues the independently constructed fractional
singleton state through the BSD weighted solver. It preserves directed odd-ring
normalization, validates the seed and retains optimum certificates. A bounded
independent pairing oracle and four stored native seed cases check this
increment. Those checks do not establish all native chronological tie choices;
full matching parity remains open.
