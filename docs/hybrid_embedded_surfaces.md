# One retained planar Tri3 Surface

An ordinary Delaunay hybrid Region can consume one retained planar Surface
containing either one actual linear Tri3 cell with three boundary Lines, or
two connected Tri3 cells with four boundary Lines. Each Line retains one
actual Line2 cell, and every sheet vertex retains its declared endpoint Point
cell. The sheet must lie strictly inside a certified convex shell. The
producer consumes literal lower node and cell UIDs, ownership, ordered corners
and coordinates.

The two-triangle sheet must form one coherent planar disk: its triangles meet
at one oppositely directed shared edge, lie on opposite sides of that diagonal,
and have a complete exterior edge set equal to the four actual boundary Line
cells. Degenerate, folded, disconnected or nonplanar sheets reject.

Every actual sheet edge, including the shared diagonal, passes through segment
recovery, and each literal triangle passes through facet recovery. Carving
uses the exterior shell facets while the interior sheet remains a protected
constraint. Tet insertion and final optimization protect all sheet faces and
edges. Every sheet facet must have exactly two incident Tet cells on opposite
geometric sides. Its vertices require closed interior Tet stars, and every
sheet edge requires a closed cyclic Tet link. Existing positive typed-cell,
exterior identity, Pyramid, ownership and analytic volume certificates remain
required.

The current admission covers one- or two-cell Tri3 sheets with only their
three or four boundary Point vertices, in one convex hybrid Region. Larger
triangulations, additional Surface or Curve vertices, Quad sheets, nested
embedded entities, multi-Region embedding, and unavailable full-search or
Steiner recovery branches remain outside this admission. An unreferenced
sheet or lower row cannot count as successful embedding.

The opt-in fixtures check lower UID/owner/XYZ retention, each two-sided Tet
facet and closed edge link, positive typed volume and exterior identity, and
an MSH4 ownership round-trip. Boundary, outside, and missing actual Tri3 or
Line2 constraints reject atomically. The two-triangle fixture uses one small
unit-cube volume generation. No volume fixture is added to the default suite.

Planarity uses the adaptive determinant magnitude and a finite dimensional
error bound. Exact orientation signs remain topology predicates. The retained
rotated lower-Surface regression checks rounded planar coordinates through the
actual public generation/constructor binding; a changed nonplanar corner
rejects before volume production. Its current gate passes14 checks onJulia1.13,
actual zero/41.41 s including compilation, body8.17 s. The preceding axis-aligned
49 and original single-Tri40 receipts qualify their preceding gate Source;
the packet retains that distinction. No volume campaign was repeated solely
for this numerical guard repair.
