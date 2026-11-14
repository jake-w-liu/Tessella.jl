# Lower Tri3 BRep partition production

The first-order Tet4 partition producer supports an owned, complete discrete
Tri3 Surface boundary. The same public
`API.mesh.partition` call consumes a complete explicit UID/partition assignment
for every original Tri3 and Tet4. Assignment-only defaults remain unchanged.

With `Mesh.PartitionCreateTopology=1`, accepted inputs retain the original
exterior Tri3 cells on actual Surface children with their original Surface
parent, UID, connectivity and physical memberships. Tet4 cells move to actual
Region children with the same guarantees. New inter-partition Tri3 facets have
opposite signed Region incidences. A supported interface boundary creates
actual closed Line2 Curve cells; the Curve's parent is its original Surface,
and child Surfaces retain the corresponding signed Curve incidences. Node
classification follows actual cell incidence, including the new lowest
dimension Curve ownership. No coordinates or sorted XYZ are used to infer
entity ownership or cell identity.

Every declared original Region boundary is checked against the actual outward
Tet facets. Topology publication then checks complete child Region facet
coverage and complete oriented Surface border coverage by Curve cells. The
normal duplicate-cell validation stays enabled.

The supported input has positive manifold Tet4 cells, complete discrete Tri3
Region boundaries, no pre-existing Curve boundaries/cells or Point/Curve
vertex classifications, and no periodic, high-order or parametric storage.
Each original Tri3's assigned partition must match an actual Tet facet carrier,
and every partition must own a Tet4. Closed seams on one original Surface are supported. The [open-seam producer](partitioned_open_seams.md) also admits
nonbranching paths across original Surfaces, while recording and preserving
every literal original border edge. New interfaces retain complete Curve-edge
coverage and complete child Region coverage. Branches, triple junctions,
missing lower carriers and pre-existing Curve/Point roles remain guarded.
Orphan node classifications remain retained.

When an interface coincides with an original Tri3 support, topology production
rejects before publication: distinct overlapping geometric roles are not yet
represented by a certified owned mesh. This is an explicit scope limit;
Gmsh's richer partition BRep is not claimed. Ghost-only production can retain
the original BRep without producing a new duplicate interface or Curve.

Ghost Regions expose aliases to the actual neighboring Tet UIDs. The owned
mesh still carries each material Tet once. MSH4 ASCII/binary export and public
`.geo` Merge import retain parents, physical memberships, signed partition
boundaries, the seam cells and ghost metadata. Canonical homogeneous block
copies remap flattened ghost and retained data references with the literal
source cell permutation. Imported partition children retain signed boundary
incidences. Ordinary unpartition removes the generated Curve/interface/ghost
entities and restores the original Tri3/Tet4 UIDs to their original parents.
Discrete entity creation validates the positive registry entity for a signed
boundary tag while retaining its incidence sign, as independently observed
through the public Gmsh4.15.2 API.

See the linked partitioning pages for supported lower-dimensional topology.
