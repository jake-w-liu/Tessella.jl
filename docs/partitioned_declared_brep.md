# Declared Curve/Point lower-BRep partition production

This increment extends the lower BRep partition producer to inputs that
carry genuine pre-existing lower-dimensional BRep: actual `Line2` cells on
discrete Curve entities and `Point15` cells on discrete Point entities,
declared through the ordinary signed `add_discrete_entity` boundary lists.
Earlier increments ([closed seams](partitioned_brep.md),
[open seams](partitioned_open_seams.md)) admitted only Tri3/Tet4 meshes whose
Surfaces had no declared Curve boundaries; this slice retains the declared
topology instead of rejecting it. A follow-on slice generalizes the cell
catalog to every linear family — Quad4 surfaces beside Tri3 and
Hex8/Prism6/Pyr5 volumes beside Tet4 — and a later slice extends it to
every complete second-order family — Line3 curves, Tri6/Quad9 surfaces and
Tet10/Hex27/Prism18/Pyr14 volumes — and later slices complete the
catalogue through complete third order — Line4 curves, Tri10/Quad16
surfaces and Tet20/Hex64/Prism40/Pyr30 volumes — and then through every
catalogued complete order, so fourth-order families (Line5, Tri15/Quad25,
Tet35/Hex125/Prism75/Pyr55) and higher partition identically.

## Input certificate

The snapshot admits every complete MSH cell family the catalogue carries —
Point15, Line(p+1), Tri and Quad surfaces and Tet/Hex/Prism/Pyr volumes at
any catalogued order — and requires every cell's classified parent to exist
as a discrete entity. Surface and volume facet/edge cycles come from the
canonical `MeshEntityTopology` outward tables (`_simplex_edge_patterns`,
`_simplex_face_patterns`); an order-p facet derives its full node tuple by
appending `p-1` edge nodes per consecutive corner pair — taken along the
facet's own traversal, so a facet edge running against the canonical edge
table picks the volume's edge nodes in reverse — then the facet's interior
slots, whose offsets accumulate triangle interiors before quadrangle
interiors in `_FACES` order exactly as the MSH pointsGenerators layout
appends them. Every volume cell must carry positive signed volume under
`_volume_cell_signed_volume` on its corner subset. Cycle equality,
reversal and support keys are layout-generic: each admits a source-index
permutation derived from the monomial lattice — a corner map σ induces the
unique affine map sending old vertex k's lattice point to new vertex σ(k)'s,
which carries edge and interior nodes to their exact slots at any order —
so no per-order tables exist. Serendipity/incomplete families (Quad8/Quad12,
Tri12, Tet22 and peers) reject. Certification is atomic — the model
and cache are untouched on failure:

- Each declared Curve must be one nonbranching directed chain of its actual
  Line cells at a single order. An open chain requires exactly two signed
  Point boundaries; the declared bound Point must own the Point15 cell on
  the corresponding literal chain end node. A closed chain must declare no
  Point boundaries.
- Every declared Point entity needs exactly one actual Point15 cell, and no
  two Points may share a vertex.
- Each declared Surface boundary lists distinct Curves whose cells cover the
  Surface's literal oriented border completely; orientation must agree with
  the declared sign. A declared higher-order Line Curve must additionally
  reuse the surface border's actual edge nodes — adjacent surface cells
  must agree on the edge's interior-node sequence up to traversal reversal,
  and a foreign sequence rejects.
- Every Curve must be declared on some Surface boundary and every Point on
  some Curve boundary; unbound lower entities reject.
- Every node referenced by a cell must carry consistent lowest-dimension
  classification: a Point-classified vertex must be that Point's cell vertex,
  a Curve-classified vertex must lie on that Curve's cells, and Surface/Region
  classifications must name an actual parent.

## Production

One partition child is created per original Curve and Point entity —
`parent=(dim, original_tag)`, `partitions` derived from literal Tri3 carrier
incidence (a Curve shared by both partition sides yields `[1,2]`; a Point
inherits the union of its declaring Curves' memberships), physical groups
copied and child membership appended, bounding box from the actual cell
rows. Original cells keep their UIDs and are emitted once under the child
owner. Child Curve boundaries bind the corresponding Point children with the
declared signs; child Surface boundaries bind the corresponding Curve
children signed by actual border orientation; declared Curve edges on an
interface border require actual lower Surface-cell carriers on both partition
sides and stay bound to the interface child. Shared nodes reclassify to the
lowest
actual incident dimension, so curve-end vertices owned by Points classify to
the Point children. No coordinates or sorted XYZ are used to infer entity
ownership or cell identity, and no Point roles are manufactured without an
actual Point15 carrier.

Ghosts, `$PartitionedEntities`/`$GhostElements` MSH4 transport, public
import and `mesh.unpartition()` follow the established producer semantics.
MSH4 import preserves the signed `boundary` tags the `$Entities` section
stores (sign encodes orientation, matching `add_discrete_entity!`'s signed
storage) for partition children and original non-partition entities alike,
so an imported declared-BRep mesh partitions through the unchanged
certificate and oriented `get_boundary` parity holds after `open`/`merge`;
unpartition restores the original entities, cells and declared signed
boundaries.

Each connected interface component produces one dim-2 entity carrying a
Tri and/or Quad block at the volume facets' own order, wound to the
primary partition side, so Hex↔Pyramid interfaces emit quads and
Tet↔Prism interfaces emit triangles at whatever order the volume facets
carry; the two sides must agree on shared edge-node sequences up to
traversal reversal. Generated seams group by parent entity and the
complete incident partition membership, then chain nonbranchingly per the
open-seam rules and emit Line(p+1) cells by the border arity;
a declared Curve carrier on an interface border binds the shared Curve
child instead.

## Junction seams and junction points

Interface-border edges are indexed globally across every generated
component, so a rim edge is no longer required to sit between exactly two
lower Surface carriers: the certificate only requires every volume incident
on the edge to contribute exactly two of its own facets (interfaces,
same-part interior facets, or exterior Surface carriers). Three or more
partitions meeting along the same edge therefore share one junction seam —
grouped by its parent entity and the complete incident partition set — whose
parent is the carrier Surface when carriers exist and otherwise the incident
Region, matching a public Gmsh 4.15.2 observation of a four-volume fan
parented by the Region. A declared Curve on a junction edge is shared by
every incident interface and its child's partition list is extended to the
full incident set.

Nodes where three or more seam edges meet — including declared Curve edge
ends in the valence count — materialize as dim-0 entities carrying Point15
cells partitioned by the union of incident seam memberships and parented by
the lowest-dimension incident seam parent. A vertex already carried by a
declared Point child reuses that child and merges the partition union into
it. Open junction seam chains bind their junction endpoints with signed
bounds (+start, −stop); closed chains keep empty bounds. Seam edge
directions are normalized at chain time so a junction arc mixing carrier and
interface-derived orientations cannot falsely trip the branching guard.

## Coincident interface support

A declared Surface cell can share its node support with a cross-partition
volume facet. Public Gmsh 4.15.2 resolves that configuration with two
Surface entities carrying duplicate cells on the same support — the
original child and a new interface entity — which collides with this
model's structural duplicate-cell certificate. The producer instead lets
the original Surface child carry the interface role itself: the cell keeps
its actual ownership and UID, the child's partition membership extends to
the complete incident set and both Region children bind it through the
declared boundary with opposite signs. No interface entity or second cell
is emitted for the coincident support. Both Regions must declare the
shared Surface — each side's coverage certificate requires the facet —
and the declared-sign winding checks force the two declarations to carry
opposite signs automatically. Declared Curve and Point children on the
coincident rim keep their roles, and a junction seam whose edge is carried
by the coincident facet binds the Surface child through the ordinary
border certificate.

## Periodic links and parametric nodes

Periodic `MixedPeriodicLink` records transport through partition unchanged —
bound to the original slave/master entity tags and the same node rows —
matching the Gmsh 4.15.2 convention: a partitioned `$Periodic` section keeps
its links on the original entities, which persist in the partition model
alongside their children. `model.periodic` constraints and
`meshing.periodic_uid_history` persist on the staged model for the same
reason, and `unpartition()` rebuilds each surviving link against the
restored node rows exactly as it already did for imported partitions.

Parametric `node_parametric` rows travel with their node rows. Classification
still retargets each vertex to its lowest-dimension incident partition
entity, but a parametric row may only land on a same-dimension child *of its
own original owner* — its parameters describe that entity's parametrization
domain — so a Surface-classified UV vertex keeps a Surface child (never a
seam Curve or junction Point) and falls back to the original entity when no
eligible child exists. The `MixedEntityData` length invariant
(`length(parameters) == entity dim`) therefore holds end to end.

## Parallel production

`General.NumThreads` requests parallel workers for the producers' per-cell
extraction stages — volume-facet derivation, signed-volume positivity,
surface-border and lower-edge extraction, interface incidence, per-block
incident-entity reduction and ghost destinations. `Threads.@threads` caps
actual concurrency at the runtime thread pool, so the option names a
requested worker count rather than a guarantee. Every parallel stage
writes disjoint per-index slots whose partial results merge serially in
the original traversal order (or by order-independent `isless` reduction
for incident entities), so the produced mesh — entities, bounds, signed
boundaries, block and UID layout — is identical under any worker count.
The default `General.NumThreads=1` keeps production serial, matching the
observed Gmsh 4.15.2 default. Like Gmsh's int option, `option()` stores a
truncated signed value (`2.7` stores `2`); non-positive or fractional
requests below `2` select serial production, and values outside the signed
32-bit range reject rather than wrapping.

## Guards and scope

Incomplete/serendipity high-order cells (Quad8/Quad12, Tri12, Tet22 and
peers), already partitioned data and seam branches within one
partition-membership group remain guarded or open. Bare meshes without lower
cells still route to the Tet4-only producer and bare non-Tet4 volume meshes
remain unsupported. `generate` always
materializes a vertex node plus Point15 cell for every dim-0 entity, so a
declared Point without a carrier is reachable only through non-generated
models (imports or direct model edits); the certificate rejects it
atomically either way.

The qualification uses two Tets, six exterior Tri3 on six Surfaces, nine
declared open Curves and five declared Points (shared and local), plus a
closed declared Curve bounding a partition interface; a Hex8+Pyr5 fixture
over five Quad4 and four Tri3 surfaces with a declared Curve/Point border
on one Quad4 side (Quad4 interface); a Tet4+Prism6 fixture sharing a Tri3
face across the two volume families (Tri3 interface); a Tet10 pair under a
declared closed Line3 Curve (Tri6 interface plus a generated Line3 seam);
a Hex27 pair (Quad9 interface and a closed four-cell Line3 seam); a
Tet10+Prism18 pair sharing a Tri6 across the two volume families (a Surface
child keeping mixed Tri6+Quad9 blocks); a three-Tet4 junction whose rim
edge closes two pairwise interfaces; the same junction with the shared
triangle split onto two declared Surfaces and its edges declared as
Curves/Points; a four-Tet4 fan whose carrier-free junction edge takes the
Region parent; a three-Tet10 junction emitting a Line3 junction seam; a
Tet20 pair whose Tri10 interface border is a declared closed Line4 Curve;
a Hex64 pair (Quad16 interface and a closed four-cell Line4 seam); a
Prism40 pair mixing Tri10 and Quad16 exterior blocks around a Quad16
interface; a Pyr30 pair sharing the Quad16 base; a Tet35 pair whose
Tri15 interface border is a declared closed Line5 Curve; a Hex125 pair
(Quad25 interface and a closed four-cell Line5 seam); a Prism75 pair
mixing Tri15 and Quad25 exterior blocks around a Quad25 interface; a
Pyr55 pair sharing the Quad25 base; a Tet56 order-5 pair (Tri21 interface
and Line6 seam); and a periodic prism
whose imported dim-2 `MixedPeriodicLink` and six UV-parametric boundary
vertices survive partition, ghosts, MSH4 round-trip and `unpartition`.
Checks cover child
parents/partitions, signed Curve→Point and Surface→Curve closures, interface
rim curves, junction seam sharing, junction point emission and reuse,
original cell/node UIDs, physical groups, ghost records,
ASCII/binary MSH4 and public `write`/`open` round-trips, `unpartition`, and
atomic topology violations across three rejection testsets including a
missing Point15 carrier, a Line3 Curve carrying a foreign edge node, a
reversed Tet10, a serendipity Quad8 block and serendipity/incomplete blocks
(Quad12, Tri12, Tet22) riding on a valid Tet20 mesh. Actual test results
are recorded in `test/README.md`.
