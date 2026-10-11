# Retained embedded Curves in a hybrid Region

An ordinary Delaunay hybrid Region can consume one embedded Curve whose current
lower mesh consists of connected linear Line2 cells. Its declared endpoint
Point cells and every referenced Curve node retain their actual UIDs, owners,
and coordinates. The complete Curve must lie strictly inside a certified convex
shell. Boundary or outside constraints and incomplete lower discretizations
are rejected before mesh publication.

The producer inserts the actual Curve vertices into the initial constructor and
recovers each literal lower segment with the existing segment-flip algorithm.
Recovered edges remain protected during carving, quality improvement, Tet
insertion, and final Tet optimization. Every Curve vertex must have a closed
interior Tet star, and every segment must have one closed cyclic Tet edge link.
The existing positive typed-cell, shell-incidence, Pyramid, ownership, and
analytic-volume certificates remain required.

This increment covers one Curve in one convex hybrid Region. Embedded Surface
recovery, multiple embedded entities, multi-Region embedding, and recovery that
needs the unavailable full-search or Steiner branches remain outside this
admission. A Curve cannot succeed by leaving its nodes or segments unreferenced.

The opt-in functional fixture uses one unit cube with two retained Line2
segments, verifies positive Tet4/Pyramid5 volume and exterior identity, checks
the exact lower cell/node transport and MSH4 ownership, and exercises outside,
boundary, and removed-cell atomic rejection. It is excluded from the default
suite because it performs actual volume generation.
