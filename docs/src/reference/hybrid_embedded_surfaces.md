# Retained planar Tri3 sheet Surfaces

An ordinary Delaunay hybrid Region can consume retained planar Surfaces made
of actual linear Tri3 cells with declared straight Line boundary loops. Each
boundary Curve keeps its complete actual Line2 cell chain — subdivided chains
with interior chain vertices are admitted — and every declared boundary or
nested Point vertex retains its actual Point cell. Sheets may carry more than
one boundary loop (holes appear as reversed inner loops), a subdivided actual
triangulation with interior vertices, and declared Point/Curve entities nested
inside the Surface. The sheet must lie strictly inside a certified convex
shell. The producer consumes literal lower node and cell UIDs, ownership,
ordered corners and coordinates.

Each sheet must form one coherent planar sheet per declared topology: the
actual boundary is exactly `L` connected cycles for `L` declared loops, the
triangle dual graph is connected, the Euler value is `2-L`, and shared edges
have opposite directed incidences with triangles on opposite sides. Its
complete exterior edge set equals the actual boundary Line cells; nested Curve
chains must be actual sheet edges. Exact projected checks reject crossing,
overlap, unmatched edge vertices, folded or disconnected triangles; nonplanar
sheets also reject. See the [polygonal scope](hybrid_polygon_sheets.md).

Every actual sheet edge, including internal diagonals and interior-vertex
spokes, passes through segment recovery, and each literal triangle passes
through facet recovery. Carving uses the exterior shell facets while the
interior sheets remain protected constraints. Tet insertion and final
optimization protect all sheet faces and edges. Every sheet facet must have
exactly two incident Tet cells on opposite geometric sides. Its vertices
require closed interior Tet stars, and every sheet edge requires a closed
cyclic Tet link. Existing positive typed-cell, exterior identity, Pyramid,
ownership and analytic volume certificates remain required.

Multiple coexisting sheets travel the [constraint collection](hybrid_constraint_collections.md)
together with embedded Points and Curves;
each Surface keeps its own witness. The single-Surface route also remains
available for a Region inventory holding exactly one Surface.

The current admission covers planar Tri3 sheets with any number of declared
straight-Line boundary loops in one convex hybrid Region. Quad sheets,
nonlinear sheet boundary curves, multi-Region embedding, constraints that land
inside the Pyramid bump region rather than the tetrahedral interior, and
unavailable full-search or Steiner recovery branches remain outside this
admission. An unreferenced sheet or lower row cannot count as successful
embedding.

Planarity uses the adaptive determinant magnitude and a finite dimensional
error bound. Exact orientation signs remain topology predicates. A changed
nonplanar corner rejects before volume production.

With `Mesh.Renumber=0`, lower node UIDs and coordinates remain unchanged.
Ordinary first-order generation can renew linear cell labels before binding
Volume inputs; owner-bound cyclic incidence is preserved.
