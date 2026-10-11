# Read-only nested STEP assembly design and Devin handoff

The user requested an unfinished-work handoff on 2026-10-11. The next owner should implement the graph producer described below. This packet contains design and two opaque import-only
controls. No production Source change, Main mutation or Julia run was made for
nested assemblies. Preserve the qualified direct54 and SI52 packets. Current
Main dependencies are relationship354b6a93, Cartesian56cb893a and BRepbc3ba588.

## Actual successful graph and dimensional observations

The two retained STEP files place intermediatePD2002 under rootPD1002.
IntermediatePD2002 then uses leafPD5 twice. Root usage2020 binds
intermediate-to-root through CDSR2024, placementPDS2021 and IDT2022. Inner
usages1020/1030 bind leaf-to-intermediate through their actual CDSR/PDS/IDT
objects. Every definition retains an actual SDR/PDS binding and literal frame
membership. The inventory and paths are in nested_reference_trace.json.

The outer rotation has rows (0,0,1), (1,0,0), (0,1,0). The right inner rotation
is a positive quarter-turn about Z. They do not commute, so this controls
composition order as well as translation.

| Actual contexts and origins | Native physical unit-box bounds |
| --- | --- |
| All mm; intermediate and leaf origins zero | [10,11]x[5,6]x[0,1] and [10,11]x[7,8]x[0,1] |
| Root m; intermediate cm with origin1cm; leaf mm with origin1mm | [9990,9991]x[2029,2030]x[0,1] and [9990,9991]x[2058,2059]x[0,1] |

Both Gmsh4.15.2 public OCC imports returned exactly two solids, 16 Points,
24 Curves and 12 Faces, with unit masses. Times were 0.0081773s and 0.0042535s,
0.0124308s total, without meshing. Actual receipts, hashes and labels are in
native_nested_import_observations.json. Additional subassembly logger paths
are not evidence of extra emitted solids; actual entities and physical bounds
establish the two placed solids.

The intermediate factor affects its literal frame origins, not the normalized
leaf geometry. It must not become another scale matrix on descendants. Retain
the current per-edge physical law R=A^T B and
t=A^T(s_parent*b-s_child*a). Compose parent and child physical transforms using
the existing proper-rigid certificate. Scale a leaf's raw source rows once.

Primary provenance uses the [NAUO schema](https://downloads.steptools.com/docs/stp_aim/html/t_next_assembly_usage_occurrence.html),
[CDSR](https://downloads.steptools.com/docs/stp_aim/html/t_context_dependent_shape_representation.html),
[RRWT](https://downloads.steptools.com/docs/stp_aim/html/t_representation_relationship_with_transformation.html),
[IDT](https://downloads.steptools.com/docs/stp_aim/html/t_item_defined_transformation.html)
and [public nonmatching OCCT construction](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_8_0/src/STEPConstruct/STEPConstruct_Assembly.cxx).
No restricted matching implementation was read.

## Source blockers and next increment

Current BRepStepRelationships.jl354b6a93 blockers:

- Line161 permits one parent and no parent/child intersection. A legitimate
  intermediate is both parent and child.
- Line207 selects only(parents), without the required unique actual graph
  root, reachability and cycle checks.
- Lines212/216 require every child to bind a solid, rejecting the actual
  frame-only intermediate representation.
- Line244 binds rep2 to the global root instead of the usage's actual parent.
  Lines247/260 likewise select the global root units for every parent frame.
- Lines274/275 emit a solid for every usage. Dropping guards would try to emit
  the frame-only intermediate and lose the actual root-to-leaf path.

The next implementation should own only this helper plus proportional opt-in
tests, fixtures and docs, based on current Main after Root's final publication.
Do not edit historical candidates or replace the Cartesian helper.

1. Build adjacency from literal NAUO parent/child references. Select one actual
   root from parent-minus-child IDs. Certify every definition/usage reachable
   and detect cycles before expansion, fill or factory use.
2. Preserve complete SDR/PDS/PD, product/formation, representation, CDSR/PDS/
   usage, relationship and IDT inventories. Each edge must bind its own actual
   parent/child representation, frame membership and dimensional witnesses.
3. For this first slice, nodes with outgoing usages contain placement frames
   only. Leaves bind one supported closed solid. Mixed geometry/children and
   multi-solid leaves retain explicit guards.
4. Expand literal usage paths deterministically. Compose physical edge rigid
   transforms in parent-to-child order; scale final leaf raw rows once. Retain
   root/intermediate/usage/leaf/solid labels and separate owned rows, UIDs and
   Surface/Region identities for every occurrence.
5. Repeated definitions may expand through multiple paths. Cache raw source
   boundaries, not placed products. Track active product IDs and consumed
   usages; cap depth128, emitted solid occurrences4096 and checked aggregate
   storage. A DAG can emit more paths than its usage-record count.
6. Preserve proper-composition, finite-coordinate, closed-shell, strict
   separation, positive-fill, ownership and storage certificates. Reject
   orphan definitions/placements/solids, cycles, multiple roots, incomplete
   references, reversed relationships, nonmember frames and unsupported units.

Qualification should use the two retained real imports, analytic corners and
unit volumes, distinct UIDs/owners/path labels, one tiny nested-mm fill and
MSH4 owner/connectivity records, plus cheap graph failures. One direct-mm and
one SI mixed-unit Surface are appropriate continuity controls. All volume and
native checks remain opt-in. Coordinate runtime ownership before any Julia.

## Further CAD work

These are remaining contexts, not newly qualified behavior:

- Mixed mapped/product workflows reject at relationship143. The mapped
  helper also assumes one SDR root and equal contexts (occurrence139). Start
  from one proper Native-successful control, then expand literal map items
  inside the graph producer with actual origin/context transport.
- Multi-solid leaves and representations with geometry plus usages reject
  at relationship208/216. They need complete source-item and occurrence
  inventories. Establish one real compound/part control before implementing.
- Conversion-based length units reject at relationship59; prefixed angular
  conversion at72; uncertainty-unit aliases at92. SI52 covers only direct
  relationship SI length transport. Ordinary/mapped root unit conversion is
  still absent. Trace actual unit refs/dimensions/factors for each new context.
- Touching/overlapping boxes reject in BRepAssembly's initial separation
  check. A genuine contact workflow needs geometric/interface ownership
  certificates; do not remove the sufficient separation guard alone.
- BRep.jl's planar straight-edge classification around1060, full-torus-only
  guard1075 and fallback1079 reject other nonlinear/trimmed topology.
  Standalone NURBS import does not construct a certified trimmed closed BRep.
  Start with one bounded curved Face/shell and actual UV/model/edge ownership.
- Retained proper8-field Cartesian scale2 is ignored by the pinned importer;
  it is not evidence of general scaling. Malformed6-field observations remain
  unqualified. Native no-ops and shared-Quad rejection must not be recast as
  fabricated missing successful features.

The next owner resumes from published Main and owns isolated implementation,
proportionate checks and normal Git publication. No runtime is held here. Full CAD/P1-P6 parity
remains unfinished.

## Repository evidence

Portable inputs and opaque observations are committed with this handoff:

- [Reference graph](../../test/artifacts/handoff_20261011/nested_relationships/nested_reference_trace.json).
- [Native observations](../../test/artifacts/handoff_20261011/nested_relationships/native_nested_import_observations.json).
- [Millimetre input](../../test/artifacts/handoff_20261011/nested_relationships/nested_mm_boxes.step).
- [Mixed-unit input](../../test/artifacts/handoff_20261011/nested_relationships/nested_m_cm_mm_shifted_boxes.step).

The STEP inputs have only trailing whitespace/newline normalization; non-whitespace
STEP tokens are unchanged. The native JSON retains original input provenance.
