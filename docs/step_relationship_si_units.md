# SI length units in direct STEP product usages

Direct product-usage assemblies carry the actual selected representation
context into each occurrence. The importer verifies the declared SI unit
categories, literal references and supported STEP prefixes. Output coordinates
use millimetres. A METRE declaration with decimal prefix exponent e has length
factor `10^(e+3)`; an omitted prefix means metres. RADIAN and STERADIAN
declarations may carry supported SI prefixes with factor `10^e`. A
`CONVERSION_BASED_UNIT` declaration resolves recursively: its
`*_MEASURE_WITH_UNIT` typed measure value times its base unit's resolved
factor, where the measure's dimension, the enclosing category part and every
chained base must agree. Optional uncertainty declarations alias any unit
entity whose resolved length factor equals the selected context's — entity
identity is no longer required.

Each leaf's literal source geometry uses its own length factor. The actual
child and parent placement-frame origins use their own respective factors.
With child origin/basis `(a,A)` and parent `(b,B)`, physical transport is
`R=AᵀB`, `t=Aᵀ(s_parent*b−s_child*a)`, and
`x_output=R*(s_child*x)+t`. No factor or ownership is inferred from output XYZ.
Literal graph completeness, frame membership, proper rotation, finite
coordinates, closed-shell, geometric independence, storage, positive-fill and
UID/owner certificates remain mandatory after conversion. Ordinary roots,
multi-root assembly records and the separate mapped-item/Cartesian producer
convert the same declared contexts: every unit-bearing shape representation
contributes its length factor to the member solids it wraps (a solid shared by
incompatible factors rejects as ambiguous), mapped occurrences scale source
rows and transform origins by the shared root context, and representations
without unit declarations keep raw coordinates.

Opaque Gmsh 4.15.2 import-only controls establish distinct source/frame
behavior: both centimetre contexts produce side-10 solids at x30:40/x50:60;
metres produce side-1000 solids at x3000:4000/x5000:6000; child-mm/parent-cm
keeps unit sides at x30:31/x59:60. Moving the actual child origin to1mm gives
x29:30/x58:59. These controls used the actual empty `Geometry.OCCTargetUnit`
default. The three latter imports took0.04054s, without meshing. Their retained
receipt is `C:/tmp/tessella_step_relationship_assembly_20261011/si_conversion_design/native_si_conversion_trace.json`.

The proportional opt-in module
`test/geometry/brep_step_relationship_si_functional_test.jl` uses those actual
files, analytic per-owner corners/volumes, one centimetre fill, MSH4 records,
the unchanged millimetre Surface control and malformed/overflow errors. It
adds no default volume generation. All52 checks passed onJulia1.13.1 with actual
exit0 in40.37s including loading; the body took5.54s with `--compile=min -O0`,
bounds checks and one thread. The actual receipt is
`C:/tmp/tessella_step_relationship_si_20261011/relationship_si_113_v1.receipt.json`.
The preceding mm-only54 receipt is historical; its old CENTI rejection is now a
malformed STEP-prefix control. No optimized compiler/general unit parity claim
follows from this bounded qualification.

[STEP SI_UNIT](https://downloads.steptools.com/docs/stp_aim/html/t_si_unit.html)
declares the actual optional prefix, unit name and derived dimensions;
[STEP SI_PREFIX](https://downloads.steptools.com/docs/stp_aim/html/t_si_prefix.html)
defines the schema's admitted prefix enumeration. [BIPM decimal factors](https://www.bipm.org/en/measurement-units/si-prefixes/)
establish the conversion values. The CAD03 increment qualifies this transport
across ordinary, mapped and relationship producers — including nested
conversion chains (YARD→FOOT→INCH→mm), factor-preserving conversion/SI
equivalence, metre fills and MSH4 records — on Julia1.12.7 and1.13.1. It still
does not claim general CAD unit-conversion parity (derived-unit measure
contexts beyond the three categories remain explicit blockers).
