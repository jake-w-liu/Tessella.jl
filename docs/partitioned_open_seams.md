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

Entity/cell births follow the producer's actual interface edge discovery and
chain emission. They do not infer native chronology from sorted coordinates.
The tiny opaque native observation has two Tets and six exterior Tri3 cells,
each on its own original Surface. Native creates three open one-Line2 Curve
entities with original Surface parents, opposite Surface signs and no Point
boundary. The owned producer matches this bounded incidence/retention behavior
while certifying Region winding against its actual Tet cells. General native
partition ordering or complete lower BRep parity is not claimed.

The functional target checks those six Surface children, the three real Curve
cells, original UIDs/physical memberships, full interface edge coverage,
closed child Region facets, positive material volumes, detached outputs,
ghost aliases, MSH4 serialization, true public Merge and unpartition. The
unchanged two-check closed-BRep default fixture is a continuity control.
The frozen Julia 1.13.1 run passed 28/28 checks with actual exit 0 in
35.9762 seconds including loading, with bounds checking, one thread and
`--compile=min -O0`. The new target passed 26 checks in 2.8 seconds; the
unchanged closed default passed 2 in 0.3 seconds. This is no optimized JIT
performance claim. The exact receipt is
`C:/tmp/tessella_takeover_partition_open_seams_20261011/receipts/open_seams_attempt1.receipt.json`.
The public-only SDK observation/script are in
`C:/tmp/tessella_partition_brep_next_design_20261011/native_multisurface_probe.json`
and the adjacent Python script. No native implementation Source was read.
