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

One occurrence returns `Mesh`; multiple independent occurrences return a
classified `MixedMesh` with separate node/cell identities, Surface/Region
owners and physical names containing the actual mapped-item path and source
solid ID. The existing complete shell, positive-fill and strictly disjoint
bounding-box certificates remain required. `fill=false` returns the retained
closed boundaries. No source nodes are welded across occurrences.

This first producer requires one explicit shape-definition root and identical
representation-context references. Proper eight-field Cartesian transformation
operator targets with unit/default scale have a separate
[rigid target admission](step_cartesian_rigid_targets.md). Multiple product
roots, representation-relationship transformations and
context-dependent shape representations remain guarded. Unreachable mapped
items, cycles, missing references, invalid directions, empty representations,
touching boxes and unsupported items reject explicitly. Limits are128 nesting
levels,4096 solid occurrences and the existing total Surface storage caps.
The operation returns an owned mesh; it does not populate analytic Model/OCC
Face maps or implement general CAD assemblies.

The small retained fixture starts from the public OCC box export and adds two
explicit mapped occurrences: a translated unit box at x=3:4 and a rotated one
at x=5:6. The opaque SDK imports both with unit masses; shifted-origin and
nested-map controls also succeed. Exact observations and earlier failed
receipts are under `C:/tmp/tessella_step_mapped_occurrences_20261011`.
The functional target passed31 checks plus2 default smoke checks on Julia1.12.7,
actual exit0 in37.27 s including30.65 s precompilation. Its small test bodies
total5.8 s, with no SDK calls. The initial optional-direction transport
failure is preserved. Tests check positive analytic volumes, literal
placement, ownership, unique identities, MSH4 ASCII/binary round trips and
precise unsupported-input errors.
