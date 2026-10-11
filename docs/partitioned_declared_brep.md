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
Hex8/Prism6/Pyr5 volumes beside Tet4.

## Input certificate

The snapshot admits linear MSH blocks only — Point15, Line2, Tri3, Quad4,
Tet4, Hex8, Prism6 and Pyr5 — and requires every cell's classified parent
to exist as a discrete entity. Surface and volume facet/edge cycles come
from the canonical `MeshEntityTopology` outward tables
(`_simplex_edge_patterns`, `_simplex_face_patterns`), and every volume cell
must carry positive signed volume under `_volume_cell_signed_volume`.
Certification is atomic — the model and cache are untouched on failure:

- Each declared Curve must be one nonbranching directed chain of its actual
  Line2 cells. An open chain requires exactly two signed Point boundaries;
  the declared bound Point must own the Point15 cell on the corresponding
  literal chain end node. A closed chain must declare no Point boundaries.
- Every declared Point entity needs exactly one actual Point15 cell, and no
  two Points may share a vertex.
- Each declared Surface boundary lists distinct Curves whose cells cover the
  Surface's literal oriented border completely; orientation must agree with
  the declared sign.
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
Tri3 and/or Quad4 block per its facet types, wound to the primary partition
side, so Hex↔Pyramid interfaces emit Quad4 cells and Tet↔Prism interfaces
emit Tri3 cells. Generated seams group by the adjacent original Surface
parent and partition membership, then chain nonbranchingly per the open-seam
rules; a declared Curve carrier on an interface border binds the shared
Curve child instead.

## Guards and scope

High-order cells, periodic or parametric topology, already partitioned
data, seam branches/triple junctions, coincident original/interface support
roles and actual parallel partition algorithms remain guarded or open.
Bare meshes without lower cells still route to the Tet4-only producer and
bare non-Tet4 volume meshes remain unsupported. `generate` always
materializes a vertex node plus Point15 cell for every dim-0 entity, so a
declared Point without a carrier is reachable only through non-generated
models (imports or direct model edits); the certificate rejects it
atomically either way.

The qualification uses two Tets, six exterior Tri3 on six Surfaces, nine
declared open Curves and five declared Points (shared and local), plus a
closed declared Curve bounding a partition interface; a Hex8+Pyr5 fixture
over five Quad4 and four Tri3 surfaces with a declared Curve/Point border
on one Quad4 side (Quad4 interface); and a Tet4+Prism6 fixture sharing a
Tri3 face across the two volume families (Tri3 interface). Checks cover child
parents/partitions, signed Curve→Point and Surface→Curve closures, interface
rim curves, original cell/node UIDs, physical groups, ghost records,
ASCII/binary MSH4 and public `write`/`open` round-trips, `unpartition`, and
eleven atomic topology violations including a missing Point15 carrier. Actual
test results are recorded in `test/README.md`.
