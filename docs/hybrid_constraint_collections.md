# Multiple retained Point, Curve and Surface constraints

One ordinary Delaunay hybrid Region can consume a complete inventory of
multiple embedded Points, connected linear Line2 Curves, and retained planar
Tri3 sheet Surfaces. Every declared entity has its own actual lower-cell
witness, bound to the same current model, Region, mesh authority, and complete
target inventory. The producer consumes the union of their referenced actual
node UIDs. Shared UIDs must retain the same lower owner and coordinates across
every witness.

Each Point, Curve and Surface remains strictly inside a certified convex
shell — inside the tetrahedral region, not only the pyramid bumps. The literal
Curve edges and each sheet's literal facets enter the existing recovery
algorithm and remain protected through insertion and final Tet optimization.
Sheet facets are recovered as interior constraints with incidence two; every
literal sheet edge stays a protected segment. The final certificates check
every declared entity: each Point, Curve and Surface vertex must have a closed
interior Tet star, every actual Curve segment and every literal sheet edge
must have a closed cyclic Tet link, and every sheet facet must keep two-sided
positive Tet incidence. An unreferenced lower row, omitted entity, or partial
Curve cannot count as successful embedding. Existing positive typed-cell,
exterior identity, Pyramid, ownership, and analytic-volume certificates remain
required.

The original single Point, single Curve, and single retained Tri3 Surface
routes remain separate. Multi-Region embedding, nonlinear embedded curves,
Quad sheets, and recovery needing unavailable full-search or Steiner branches
remain outside this collection admission.

The opt-in fixture combines two direct Points with two retained Curve chains
containing four actual Line2 segments in one unit cube, plus a second fixture
adding a retained Tri3 sheet to the same inventory. They check complete lower
UID/owner/XYZ and ordered-cell retention, positive typed volume, exterior
identity, protected edge incidence, two-sided sheet facets, and MSH4
ownership. Outside Point, boundary Curve, and removed actual Point/Line-cell
constraints reject atomically. No volume generation is added to the default
test suite.
