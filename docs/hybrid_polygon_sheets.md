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

The opt-in [public fixture](../test/interfaces/api_hybrid_polygon_sheet_functional_test.jl)
generates actual pentagon/hexagon lower meshes through the public API, verifies
actual lower UID/Physical incidence and independent facet/edge/volume checks,
and rejects one publicly reversed actual Tri without publishing a Volume. No
default Volume generation workload is added. The fixture selects
`Mesh.Renumber=0` to check literal lower node UID retention. Ordinary generation
renews linear cell factory identities before binding the Volume input; the
fixture checks the resulting actual nonzero cell UIDs and unchanged cyclic
owner incidence. The ordinary public renumber pass remains available; its
physical-priority numbering may also change labels while preserving the actual
owner-bound cells and coordinates.

Reversing one retained Tri through the public mutation API invalidates the
owned Surface constructor. The subsequent Volume request rejects at that
ownership guard and leaves the retained public mesh unchanged.

The final bounded Julia1.12.7 reader passed48/48, actual exit0 in26.823576 s
with minimal compilation and bounds checks. Its45 positive and3 rejection
checks preserve full final certificates. Earlier failures and the actual
pre/post-optimizer missing-facet trace remain in the qualification packet.
