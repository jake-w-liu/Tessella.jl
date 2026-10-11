# Multiple retained Point and Curve constraints

One ordinary Delaunay hybrid Region can consume a complete inventory of
multiple embedded Points and connected linear Line2 Curves. Every declared
entity has its own actual lower-cell witness, bound to the same current model,
Region, mesh authority, and complete target inventory. The producer consumes
the union of their referenced actual node UIDs. Shared UIDs must retain the
same lower owner and coordinates across every witness.

Each Point and Curve remains strictly inside a certified convex shell. The
literal Curve edges enter the existing recovery algorithm and remain protected
through insertion and final Tet optimization. The final certificates check
every declared entity: each Point and Curve vertex must have a closed interior
Tet star, and every actual Curve segment must have a closed cyclic Tet link.
An unreferenced lower row, omitted entity, or partial Curve cannot count as
successful embedding. Existing positive typed-cell, exterior identity,
Pyramid, ownership, and analytic-volume certificates remain required.

The original single Point, single Curve, and single retained Tri3 Surface
routes remain separate. Inventories containing a Surface, multi-Region
embedding, and recovery needing unavailable full-search or Steiner branches
remain outside this collection admission.

The opt-in fixture combines two direct Points with two retained Curve chains
containing four actual Line2 segments in one unit cube. It checks complete
lower UID/owner/XYZ and ordered-cell retention, positive typed volume, exterior
identity, protected edge incidence, and MSH4 ownership. Outside Point,
boundary Curve, and removed actual Point/Line-cell constraints reject
atomically. No volume generation is added to the default test suite.
