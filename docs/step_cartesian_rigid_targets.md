# Rigid Cartesian targets in STEP mapped occurrences

`import_step` accepts a proper eight-field `CARTESIAN_TRANSFORMATION_OPERATOR_3D`
as a mapped-item target with unit or omitted scale. The explicit inherited fields
are two names and an optional description, followed by axis1, axis2, local origin, scale,
and axis3, as specified by the [STEP entity schema](https://downloads.steptools.com/docs/stp_aim/html/t_cartesian_transformation_operator_3d.html).
The [inherited description is optional](https://downloads.steptools.com/docs/stp_aim/html/t_functionally_defined_transformation.html);
its omission is covered by an additional public native control. The actual
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
the existing proper-rotation certificate; complete shell, disjoint occurrence,
positive fill, ownership and unique identity checks remain required.

This admission covers unit/default scale. A retained proper eight-field native
scale2 observation shows that the pinned importer leaves the box unit-sized;
it supplies no general scaling result. Nonunit scales remain explicitly guarded
in this bounded implementation. The earlier malformed six-field observation is
not a successful transform witness. Multiple product roots, general assembly
relationships remain outside this mapped-target admission; a separate
[direct product-usage producer](step_relationship_assemblies.md) now admits
complete rigid NAUO/CDSR/IDT paths with actual millimetre contexts.

Five opaque Gmsh4.15.2 import-only controls cover default axes, a quarter-turn,
inferred axis3, a shifted mapping origin, and an arbitrary rigid frame. Each
imports two separate unit-mass boxes. The shifted-origin control distinguishes
the direct Cartesian rule from the AXIS2 rule. Observations are retained in
`C:/tmp/tessella_step_cartesian_rigid_20261011/native_cartesian_rigid_probe.json`.

The focused [functional regression](../test/geometry/brep_cartesian_rigid_functional_test.jl)
checks those five import-only cases, one actual quarter-turn fill, classified
MSH4 ASCII/binary transport and error controls. Julia1.12.7 passed43/43,
actual exit0, warm wall5.4328898s and body4.9577108s, with checked bounds,
minimal compilation and one thread. Exact UID-to-owner/family/ordered-corner
records survive MSH4 family regrouping. Earlier test-only failures and their
receipts remain archived; the selected Source stayed unchanged throughout
qualification. This fixture is opt-in; default tests gain no volume generation
or SDK calls. These timings describe the small fixture, with no production
performance claim.
