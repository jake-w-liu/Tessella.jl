# Tet partition interfaces and ghost queries

`API.mesh.partition` supports a bounded producer for bare, positive first-order Tet4
discrete Regions. The existing assignment-only default is unchanged.
`Mesh.PartitionCreateTopology=1` creates child Regions and one shared Tri3
interface entity per connected component between two actual child Regions.
Both Regions reference that interface with opposite signs. Exterior facets
and Curve/Point boundary representations are not reconstructed. The source
Tet order, actual node/cell UIDs and material volumes are retained.

`Mesh.PartitionCreateGhostCells=1` adds actual Ghost Region query aliases for
cells incident to nodes shared with another partition. Physical cells occur
once in the owned mesh and in the MSH4 Elements section. Ghost metadata and
entity queries refer to those same actual cell UIDs. Original physical groups
are inherited by child Regions; no new physical groups are invented.

Both options default to zero to preserve Tessella's existing partition
assignment behavior. Gmsh's topology option defaults to one. This producer does
not claim automatic METIS assignments or general-family partition parity.
See the [partition tutorial](https://gmsh.info/doc/texinfo/gmsh.html#t21_003a-Mesh-partitioning)
and [partition API](https://gmsh.info/doc/texinfo/gmsh.html#index-gmsh_002fmodel_002fmesh_002fpartition)
for the corresponding Gmsh interface.

This document records the bare-Tet producer. The later
[complete Tri3 BRep producer](partitioned_brep.md) also admits existing Surface
children and a closed interface seam on one original Surface. Other lower
cell/BRep contexts, referenced unsupported lower-classified Tet vertices, P2 or
higher cells, periodic constraints, existing partition metadata, uncovered
caller records and empty declared partitions reject before publication.
Orphan nodes remain owned and retained. Standard MSH4 partition and ghost
sections survive both the owned writer and API model-view writer, with ghost
indices remapped by actual UID after serialization reorders cell blocks.
These products require MSH4; MSH2 export rejects before opening the target file.
