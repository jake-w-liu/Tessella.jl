# One retained planar Tri3 disk Surface

An ordinary Delaunay hybrid Region can consume one retained planar Surface
containing `n-2` actual linear Tri3 cells with three to six boundary Lines.
Each Line retains one
actual Line2 cell, and every sheet vertex retains its declared endpoint Point
cell. The sheet must lie strictly inside a certified convex shell. The
producer consumes literal lower node and cell UIDs, ownership, ordered corners
and coordinates.

The sheet must form one coherent planar disk: its actual boundary is one
connected cycle, the triangle dual graph is connected, and shared edges have
opposite directed incidences with triangles on opposite sides. Its complete
exterior edge set equals the actual boundary Line cells. Exact projected
checks reject crossing, overlap, unmatched edge vertices, folded or disconnected
triangles; nonplanar sheets also reject. See the [polygonal scope](hybrid_polygon_sheets.md).

Every actual sheet edge, including the shared diagonal, passes through segment
recovery, and each literal triangle passes through facet recovery. Carving
uses the exterior shell facets while the interior sheet remains a protected
constraint. Tet insertion and final optimization protect all sheet faces and
edges. Every sheet facet must have exactly two incident Tet cells on opposite
geometric sides. Its vertices require closed interior Tet stars, and every
sheet edge requires a closed cyclic Tet link. Existing positive typed-cell,
exterior identity, Pyramid, ownership and analytic volume certificates remain
required.

The current admission covers one endpoint-only Tri3 disk with three to six
boundary Point vertices, in one convex hybrid Region. Larger polygons, holes,
additional Surface or Curve vertices, Quad sheets, nested
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

The subsequent polygonal regression passes48 checks, Julia1.12 actual zero/
26.82 s, including two actual5/6-gon Volume fills. It exposed and fixed a real
optimizer protection defect: literal cyclic facets/directed edges are now
canonicalized only for internal lookup keys. Full final facet/star/link checks
remain required. With `Mesh.Renumber=0`, literal lower node UIDs and XYZ remain;
ordinary first-order cell factories can renew linear cell labels before binding
Volume inputs. The test checks actual resulting cell identities and unchanged
owner-bound cyclic incidence. No Volume fixture enters the default suite.
