# Rigid Cartesian targets in STEP mapped occurrences

`import_step` accepts a proper eight-field `CARTESIAN_TRANSFORMATION_OPERATOR_3D`
as a mapped-item target with unit or omitted scale. The explicit inherited fields
are two names and an optional description, followed by axis1, axis2, local origin, scale,
and axis3, as specified by the [STEP entity schema](https://downloads.steptools.com/docs/stp_aim/html/t_cartesian_transformation_operator_3d.html).
The [inherited description is optional](https://downloads.steptools.com/docs/stp_aim/html/t_functionally_defined_transformation.html);
it may be omitted. The actual
mapping origin remains a validated `AXIS2_PLACEMENT_3D`.

A Cartesian target applies its own rigid basis and origin directly to the
source coordinates. Its transform does not subtract or rotate through the
mapping-origin placement. This is the distinct Cartesian-target branch in the
public [OCCT mapped-item translator](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_8_0/src/StepToTopoDS/StepToTopoDS_MakeTransformed.cxx).
It differs from the retained AXIS2-to-AXIS2 coordinate-frame rule.

Missing axis1/axis2 use X/Y defaults. If axis3 is absent, the actual first and
second directions determine it by their cross product. X is projected
perpendicular to the resulting Z, and Y is Z crossed with X. The
[public OCCT transform conversion](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_8_0/src/StepToGeom/StepToGeom.cxx)
establishes this rigid frame. Zero, parallel or nonfinite directions, missing
references and malformed records reject. Every nested composition still passes
the scaled-rigid certificate — orthogonal columns of equal norm — so composed
products of unit/rigid and scaled operators admit while shear and non-uniform
scale reject; complete shell, disjoint occurrence,
positive fill, ownership and unique identity checks remain required.

Nonzero scale applies `scale × R_axes` about the local origin. A negative
scale mirrors the occurrence and reverses emitted winding; zero scale is
singular and rejects. Gmsh 4.15.2 can ignore this scale field, so scaled imports
need not have identical bounds. Multiple product roots and general assembly
relationships are outside this mapped-target interface; a separate
[direct product-usage producer](step_relationship_assemblies.md) now admits
complete rigid NAUO/CDSR/IDT paths with actual
[SI source/placement length conversion](step_relationship_si_units.md).
