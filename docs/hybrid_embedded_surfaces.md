# One retained interior Tri3 Surface

An ordinary Delaunay hybrid Region can consume one retained planar Surface
containing one actual linear Tri3 cell, bounded by three actual Line2 cells and
their three endpoint Point cells. The sheet must lie strictly inside a
certified convex shell. The producer consumes its literal lower node and cell
UIDs, ownership, and coordinates; it does not reconstruct a sheet from nearby
mesh coordinates.

The actual sheet edges pass through segment recovery, and its literal triangle
passes through facet recovery. Carving uses the exterior shell facets while
the interior sheet remains a protected constraint. Tet insertion and final
optimization protect its face and edges. The final sheet must have exactly two
incident Tet cells on opposite geometric sides. Its vertices require closed
interior Tet stars, and its boundary edges require closed cyclic Tet links.
Existing positive typed-cell, exterior identity, Pyramid, ownership, and analytic
volume certificates remain required.

This first increment covers one Tri3 sheet in one convex hybrid Region. Larger
sheet triangulations, Quad sheets, nested embedded entities, multi-Region
embedding, and the unavailable full-search or Steiner recovery branches remain
outside this admission. An unreferenced sheet or lower row cannot count as
successful embedding.

The opt-in functional fixture performs one real unit-cube volume generation,
checks lower UID/owner/XYZ retention, two-sided Tet incidence, positive typed
volume and exterior identity, and verifies an MSH4 ownership round-trip.
Boundary, outside, and missing actual Tri3-cell constraints reject atomically.
No volume fixture is added to the default suite.
