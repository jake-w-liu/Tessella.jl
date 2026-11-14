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

[STEP SI_UNIT](https://downloads.steptools.com/docs/stp_aim/html/t_si_unit.html)
declares the actual optional prefix, unit name and derived dimensions;
[STEP SI_PREFIX](https://downloads.steptools.com/docs/stp_aim/html/t_si_prefix.html)
defines the schema's admitted prefix enumeration. [BIPM decimal factors](https://www.bipm.org/en/measurement-units/si-prefixes/)
establish the conversion values. Nested conversion chains, including
YARD→FOOT→INCH→mm, use the same resolution rules. Derived-unit measure
contexts outside length, angle and solid angle remain unsupported.
