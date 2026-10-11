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
The actual opaque public Gmsh4.15.2 two-Tet observations confirm child/parent
bindings, an interface with no lower boundary, opposite Region signs, node
ownership, inherited physical groups and reciprocal ghost query aliases.
Primary public references are [the partition tutorial](https://gmsh.info/doc/texinfo/gmsh.html#t21_003a-Mesh-partitioning)
and [the partition API](https://gmsh.info/doc/texinfo/gmsh.html#index-gmsh_002fmodel_002fmesh_002fpartition).
Opaque observations are in C:/tmp/tessella_partition_design_20261011; no
implementation Source from Gmsh or any matching library was read or ported.

Lower cell families/BRep, referenced lower-classified Tet vertices, P2 or
higher cells, periodic constraints, existing partition metadata, uncovered
caller records and empty declared partitions reject before publication.
Orphan nodes remain owned and retained. Standard MSH4 partition and ghost
sections survive both the owned writer and API model-view writer, with ghost
indices remapped by actual UID after serialization reorders cell blocks.
These products require MSH4; MSH2 export rejects before opening the target file.

The functional fixture is two Tets sharing one triangular face, total volume
1/3. It checks actual public topology, owner/parent/classification queries,
cached lookup lifetime, ghost aliases, positive signed volumes, detached
outputs, ASCII/binary MSH4 transport, unpartitioning, legacy assignments and
atomic unsupported requests. The final bounded Julia 1.13.1 run passed 66/66
checks with actual exit 0 in 37.0038 seconds including loading. The four groups
reported 4.7, 0.2, 0.3 and 0.2 seconds. The last group is the two-check default
fixture. This is a minimal-compilation/O0 result, not an optimized JIT timing.
The exact receipt is
`C:/tmp/tessella_takeover_partition_20261011/receipts/partition_113_attempt4.receipt.json`.

Prior receipts remain preserved. Attempt 1 caught a wrong module qualification
in the producer, and attempt 3 caught the missing binary GhostElements footer.
The final Source fixes both and also permits the registered binary flags to
toggle 0/1/0 through checked native integer coercion. The fixture was corrected
to inspect NodeData's node references and to import through the actual public
`open_geo!`/`Merge` path. Its earlier `API.open` call resolved to Base.open;
that attempt did not prove mesh import and leaked a test-owned stream. The
final public import passed without a file cleanup error.
