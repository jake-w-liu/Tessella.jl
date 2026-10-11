# Direct STEP product-usage assemblies

The importer stages literal direct placements through
`CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`, its placement
`PRODUCT_DEFINITION_SHAPE`, and `NEXT_ASSEMBLY_USAGE_OCCURRENCE`. Each usage
binds the actual parent and child product definitions to their own
`SHAPE_DEFINITION_REPRESENTATION` and supported shape representation.

The transformed shape relationship must name child representation first and
parent representation second. Its `ITEM_DEFINED_TRANSFORMATION` binds an
actual `AXIS2_PLACEMENT_3D` item in each corresponding representation. The
proper rigid transform uses the retained coordinate-frame law
`R = AᵀB`. Both representation contexts must be distinct three-dimensional
instances with verified actual SI declarations. Source geometry and each frame
origin now use their own declared length factor; output coordinates use
millimetres. [SI length transport](step_relationship_si_units.md) documents
`t = Aᵀ(s_parent*b − s_child*a)` and the source-row conversion. Angular
declarations remain unprefixed radians/steradians. An optional uncertainty
measure must retain the literal selected context length-unit reference.

Every root product may have direct or nested leaf usages, including repeated
uses of the same solid or intermediate definition. The usage graph must be a
reachable acyclic forest over literal `PRODUCT_DEFINITION` nodes; independent
roots — assembly trees and standalone products alike — expand side by side in
identifier order after the same reachability and cycle checks. Any definition,
whether root, intermediate or leaf, may also place literal closed solids and
`MAPPED_ITEM` occurrences next to its placement frames; a leaf must carry at
least one solid-bearing item. Each edge binds the actual child and parent
representations, context witnesses and member frames, and physical rigid
transforms compose in parent-to-child order; the shared intermediate length
factor converts only its literal frame origins, while the leaf's raw source
rows scale once by its own declared factor. A `MAPPED_ITEM` inside a product
representation binds an actual `REPRESENTATION_MAP` whose source representation
carries its own declared SI context, expands transitively through nested maps,
and composes under the same `R = AᵀB, t = Aᵀ(s_target*b − s_map*a)` law —
including literal `CARTESIAN_TRANSFORMATION_OPERATOR_3D` targets. Every
root-to-leaf usage path emits a separate occurrence, so a DAG can emit more
occurrences than its usage count. The importer consumes the complete product,
shape, usage, relationship, transform, map and solid inventories before
emission. Missing, ambiguous or orphan bindings, cyclic usages or mapped
representations, unreached definitions, unselected maps/items/representations,
unplaced solids, geometry-free leaves, reversed relationships, nonmember
transform frames, unsupported contexts/transforms and touching or overlapping
occurrence boxes reject explicitly. Conversion-based units, prefixed angular
declarations, uncertainty aliases and ordinary/mapped root unit conversion are
covered by the [SI unit doc](step_relationship_si_units.md).

Every occurrence receives separate mesh rows, UIDs and discrete Surface/Region
owners. Names retain the root product and the literal usage/product/solid path
in path order; deeper paths compact product labels and identifier segments to
stay inside the bounded MSH physical-name transport. Existing closed-shell,
independent bounding-box, positive-fill and storage certificates apply.
Definitions are not emitted at their original placement as extra solids.
This mesh import does not reconstruct analytic CAD parameterization.

## Evidence and qualification

Two opaque Gmsh 4.15.2 public OCC imports use the same literal unit-solid
definition through two product-usage placement paths. Both report exactly two
unit solids, with 16 Points, 24 Curves and 12 Faces. Base bounds are
`[3,4]×[0,1]×[0,1]` and `[5,6]×[0,1]×[0,1]`. Moving the actual child frame
origin to `(1,0,0)` gives x bounds `[2,3]` and `[4,5]`, independently checking
the coordinate-frame law. Together these imports took 0.0265 seconds without
meshing. Artifacts are retained in
`C:/tmp/tessella_step_relationship_assembly_20261011/native_relationship_assembly_probe.json`.

A separate 0.01695-second opaque centimetre control imports the same stored
geometry as side-10 solids, with volumes 1000 and x bounds `[30,40]` and
`[50,60]`. This confirmed that identical context-unit IDs do not justify raw
coordinate admission outside millimetres. Its receipt is
`native_relationship_unit_control.json`. The subsequent SI producer implements
the real source/frame conversions; it no longer rejects valid centimetres.

The opt-in functional module is
`test/geometry/brep_step_relationship_functional_test.jl`. It checks both
literal imports, analytic corners/volumes, distinct ownership/UIDs, one actual
fill, MSH4 ASCII/binary transport and malformed graph rejection. The preceding
mm-only Source passed all 54 assertions on Julia 1.13.1 with actual exit 0 in 44.03 seconds including
loading; the test body took 7.82 seconds. The command used `--compile=min -O0`,
bounds checks and one thread. The actual receipt is
`C:/tmp/tessella_step_relationship_assembly_20261011/relationship_113_v1.receipt.json`.
The obsolete centimetre-rejection control now checks an undeclared STEP prefix,
after the separate 52-assertion SI conversion qualification passed. This bounded
history does not establish optimized compiler or general assembly/unit-conversion
parity. It adds no default volume run.

The nested-DAG module `test/geometry/brep_step_relationship_nested_functional_test.jl`
imports the retained two-level millimetre and metre/centimetre/millimetre
fixtures, checks the composed noncommuting bounds against the opaque public
observations, a shared-intermediate diamond expansion, one actual fill with
MSH4 ASCII/binary identity transport, and eight malformed graph/cycle/orphan
rejections. The direct and mapped producers remain unchanged controls.

The mixed module `test/geometry/brep_step_relationship_mixed_functional_test.jl`
imports the retained `brep_relationship_mixed_roots.step` fixture — a
geometry-bearing root assembly with two product usages (a two-solid leaf and a
mapped leaf), a transitive `REPRESENTATION_MAP` source rep and a standalone
second root — and checks the five literal occurrence bounds and volumes, the
path-qualified `MAPPED_ITEM` names, a shared-map second occurrence, a
`CARTESIAN_TRANSFORMATION_OPERATOR_3D` mapped target, one actual fill with
MSH4 ASCII/binary identity transport, and eleven orphan/malformed rejections
including a mapped representation cycle and an unplaced solid.

The [STEP relationship schema](https://downloads.steptools.com/docs/stp_aim/html/t_representation_relationship_with_transformation.html)
declares the distinct-context and transform-item membership requirements.
The [CDSR schema](https://downloads.steptools.com/docs/stp_aim/html/t_context_dependent_shape_representation.html)
binds the represented product relation; the
[NAUO schema](https://downloads.steptools.com/docs/stp_aim/html/t_next_assembly_usage_occurrence.html)
provides parent and child definitions. Public nonmatching
[OCCT assembly construction](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_8_0/src/STEPConstruct/STEPConstruct_Assembly.cxx)
builds this literal child-to-parent representation and usage chain.
