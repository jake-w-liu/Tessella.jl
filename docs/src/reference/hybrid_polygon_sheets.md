# Retained polygonal hybrid sheets

The retained planar hybrid Surface producer accepts one endpoint-only disk with
three to six declared straight boundary Curves and exactly `n-2` actual Tri3
cells. This adds pentagon and hexagon sheets to the existing triangle and
quadrilateral scope. Each endpoint and boundary segment must have its actual
Point and Line cell; the actual Surface constructor, UIDs, classifications,
coordinates and cell order remain bound to those retained products.

Admission certifies one connected boundary cycle and a connected triangle dual
graph. Every shared edge has two opposite directed incidences and projected
opposite-side triangles. Exact orientation predicates reject proper edge
crossings, collinear overlaps, unmatched vertices on edges, and triangle
interior overlap. The existing finite determinant planarity bound remains;
coordinates are not flattened or snapped.

The existing hybrid pipeline recovers every actual facet and protects every
sheet edge. Final checks retain the source UIDs/coordinates, spherical vertex
stars, cyclic edge links, and exactly two oppositely sided Tet incidences for
each sheet Tri. Optimizer protection uses canonical internal facet/edge support
keys while retaining the actual cyclic constructor and node/cell presentation.
Positive cell orientation, closed exterior shell and analytic
Region volume remain mandatory. This scope is one retained planar Surface in
one certified convex hybrid Region; it does not add holes, Surface interior
vertices, larger polygons, multiple sheets, or general native ordering parity.

`Mesh.Renumber=0` preserves lower node UIDs. Ordinary generation can renew
linear cell labels before binding Volume inputs while preserving owner-bound
cyclic incidence. The public renumber pass can change labels according to
physical-group priority while preserving cells, owners and coordinates.

Reversing one retained Tri through the public mutation API invalidates the
owned Surface constructor. The subsequent Volume request rejects at that
ownership guard and leaves the retained public mesh unchanged.
