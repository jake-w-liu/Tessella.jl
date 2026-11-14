# STEP mapped solid occurrences

`import_step` follows one explicit `SHAPE_DEFINITION_REPRESENTATION` root that
contains `MAPPED_ITEM` occurrences through actual `REPRESENTATION_MAP` records.
The mapped source can contain supported planar BRep solids or recognized
primitives. Nested maps expand in literal representation-item order; source
definitions are not emitted as extra solids at their original positions.

The mapping origin remains an actual `AXIS2_PLACEMENT_3D`. An AXIS2 target
uses both actual Cartesian origins and normalized rigid frames; optional
directions use the STEP [axis-placement schema](https://www.steptools.com/docs/stp_aim/html/t_axis2_placement_3d.html)
and [first projected axis](https://www.steptools.com/docs/stp_aim/html/t_first_proj_axis.html).
For axis-column matrices A/B and origins a/b, this importer uses the pinned
native coordinate-frame law `R=AᵀB, t=Aᵀ(b-a)`. Each parent acts after its
child. This follows the observed Gmsh4.15.2/OCCT import behavior, including a
translated mapping origin followed by a rotated target; it does not infer a
different placement law from the bounding boxes. The public nonmatching
[OCCT mapped-item translator](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_8_0/src/StepToTopoDS/StepToTopoDS_MakeTransformed.cxx)
specifies the target-to-origin coordinate change.
When the shared root representation context declares units, every source row
and transform origin converts by its resolved SI or conversion-based length
factor to millimetres; contexts without unit declarations keep raw rows.

One occurrence returns `Mesh`; multiple independent occurrences return a
classified `MixedMesh` with separate node/cell identities, Surface/Region
owners and physical names containing the actual mapped-item path and source
solid ID. The existing complete shell, positive-fill and exact geometric
independence certificates remain required. `fill=false` returns the retained
closed boundaries. No source nodes are welded across occurrences.

This first producer requires one explicit shape-definition root and identical
representation-context references. Proper eight-field Cartesian transformation
operator targets (rigid, uniformly scaled, or mirrored) have a separate
[target admission](step_cartesian_rigid_targets.md). The
[direct product-usage producer](step_relationship_assemblies.md) owns
NAUO/CDSR/IDT product graphs, including multiple roots, per-definition solids
and `MAPPED_ITEM` occurrences inside product representations with per-context
SI transport. Unreachable mapped
items, cycles, missing references, invalid directions, empty representations,
touching or intersecting occurrences and unsupported items reject explicitly. Limits are128 nesting
levels,4096 solid occurrences and the existing total Surface storage caps.
The operation returns an owned mesh; it does not populate analytic Model/OCC
Face maps or implement general CAD assemblies.
