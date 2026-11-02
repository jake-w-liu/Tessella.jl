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

Multiple coexisting sheets travel the [constraint collection](
hybrid_constraint_collections.md) together with embedded Points and Curves;
each Surface keeps its own witness. The single-Surface route also remains
available for a Region inventory holding exactly one Surface.

The current admission covers planar Tri3 sheets with any number of declared
straight-Line boundary loops in one convex hybrid Region. Quad sheets,
nonlinear sheet boundary curves, multi-Region embedding, constraints that land
inside the Pyramid bump region rather than the tetrahedral interior, and
unavailable full-search or Steiner recovery branches remain outside this
admission. An unreferenced sheet or lower row cannot count as successful
embedding.

The opt-in fixtures check lower UID/owner/XYZ retention, each two-sided Tet
facet and closed edge link, positive typed volume and exterior identity, and
an MSH4 ownership round-trip. Boundary, outside, and missing actual Tri3 or
Line2 constraints reject atomically. The functional file covers the original
single-Tri3 sheet, a transfinite-subdivided quad sheet with an interior
vertex, a two-loop annulus, and two coexisting sheets. No volume fixture is
added to the default suite.

Planarity uses the adaptive determinant magnitude and a finite dimensional
error bound. Exact orientation signs remain topology predicates. The retained
rotated lower-Surface regression checks rounded planar coordinates through the
actual public generation/constructor binding; a changed nonplanar corner
rejects before volume production. Its current gate passes14 checks onJulia1.13,
actual zero/41.41 s including compilation, body8.17 s. The preceding axis-aligned
49 and original single-Tri40 receipts qualify their preceding gate Source;
the packet retains that distinction. No volume campaign was repeated solely
for this numerical guard repair.

The subsequent polygonal regression passes48 checks, Julia1.12 actual zero/
26.82 s, including two actual5/6-gon Volume fills. It exposed and fixed a real
optimizer protection defect: literal cyclic facets/directed edges are now
canonicalized only for internal lookup keys. Full final facet/star/link checks
remain required. With `Mesh.Renumber=0`, literal lower node UIDs and XYZ remain;
ordinary first-order cell factories can renew linear cell labels before binding
Volume inputs. The test checks actual resulting cell identities and unchanged
owner-bound cyclic incidence. No Volume fixture enters the default suite.
