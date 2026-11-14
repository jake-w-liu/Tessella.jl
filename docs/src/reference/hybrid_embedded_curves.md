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

This document covers one Curve in one convex hybrid Region. Separate producers
support [retained Tri3 Surfaces](hybrid_embedded_surfaces.md) and
[multiple Point, Curve and Surface constraints](hybrid_constraint_collections.md).
Multi-Region embedding and recovery requiring unavailable full-search or
Steiner branches are unsupported. Every embedded Curve node and segment must
be referenced by the resulting mesh.
