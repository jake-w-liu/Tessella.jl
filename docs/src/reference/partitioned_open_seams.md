# Partition interface seams across original Surfaces

The Tri3/Tet4 partition producer also accepts a new interface whose border
crosses several original discrete Surfaces. Its actual lower Tri3 incidences
select the seam's original Surface parent and signed child Surface incidences.
Each nonbranching connected open path or closed cycle emits actual Line2 Curve
cells. Original Tri3/Tet4 UIDs, connectivity, parents and physical memberships
remain retained on their actual partition children.

The input can have a partial lower BRep: original Surface triangle borders
without declared Curve entities. Those borders are recorded from the literal
oriented source edges before partitioning. A child Surface border must then be
covered by an actual signed new seam or exactly retain the corresponding
original border edge and orientation. New interface Surfaces have no inherited
border allowance: every actual interface border edge must be covered by its
Curve cells. Complete child Region outward facet coverage stays required.
For an originally closed Surface there are no inherited border edges, so the
previous complete child Surface coverage certificate remains required.

Open Curve paths have no Point boundary when the original input has none. This
preserves the partial lower topology seen through the public Gmsh4.15.2 API;
it does not reconstruct missing Point or original Curve topology. Curve paths
must be connected, have at most one incoming and outgoing edge per node, and
use only actual seam incidences. Disconnected components emit separate Curves.
Branches, triple junctions, missing lower carriers, coincident original/new
Tri3 roles, pre-existing Curve/Point cell or node roles, periodic data and
high-order/parametric data still reject before publication. Normal mesh
validation and the duplicate-cell guard stay enabled.

Entity and cell ordering follows interface edge discovery and chain emission;
it is not inferred from sorted coordinates. General Gmsh partition ordering
and complete lower-BRep parity are not guaranteed.
