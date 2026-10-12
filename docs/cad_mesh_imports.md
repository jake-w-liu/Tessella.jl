# Independent solid mesh imports

`import_step(path)` and `import_iges(path)` import the supported recognized
solid layouts into meshes. Existing single-solid imports return `Mesh`.
Multiple explicit independent source roots return a classified `MixedMesh`.
Use `fill=false` to obtain their closed triangle boundaries without filling.

Each imported solid retains one discrete boundary Surface and Region, separate
node/element UIDs, and source identity in physical names. STEP entity IDs or
IGES directory-entry IDs determine solid ordering. Filled boundary triangles
come from actual tetrahedron incidence; boundary/interior node owners follow
that incidence. IGES primitive and solid
root rigid directory-entry transforms are applied.

Strictly separated bounding boxes certify independent solids directly. When
boxes overlap an exact geometric certificate decides: facet pairs sharing
exactly the same vertex positions with opposed outward windings are
conforming contact and weld — each coincident pair becomes one shared
Surface entity bound by both Regions, emitted once, with rim nodes
deduplicated by exact position. Vertex/edge-only coincident contacts are
conforming too but weld nothing: each solid keeps its own coincident
boundary node, matching OCC import semantics. Filled solids weld only
when both
independent fills emit the identical interface triangulation; a facet
shared by three solids, a fill crossing the interface rim, and
non-coincident touching or intersecting boundary facets (partial overlap,
T-junction or partial-edge contact, same-region duplication) all reject,
as does one
solid enclosing or dipping into the other. Genuinely disjoint solids are
admitted. Boxes are never used to discover, split or
replace source roots. Point-cloud box imports retain their original single-solid behavior.
Multi-root IGES identities must be positive, odd, unique and consistent with
the directory-entry type. Unsupported geometric entity types are rejected;
supported type-406 property records can accompany the geometry. A single-face
`MANIFOLD_SOLID_BREP` on an untrimmed `TOROIDAL_SURFACE` imports as a full
torus; a single-face shell on a `SPHERICAL_SURFACE` whose boundary loops
collapse onto the two placement poles imports as a full sphere (boundary
vertices off either pole mean a genuine trim and reject). A three-face shell
of one `CYLINDRICAL_SURFACE` lateral face plus two `PLANE` caps imports as a
full cylinder when the edge-curve inventory certifies it: lateral bounds hold
only closed rim `CIRCLE`s centred on the axis with matching radius and axial
normal plus optional axial seam segments, and each cap bounds closed circles
sharing one rim inside its plane; partial arcs, off-axis or wrong-radius
circles, unclosed rims and non-axial seams are genuine trims and reject. A
`CONICAL_SURFACE` lateral face plus `PLANE` caps imports as a full cone or
frustum when the same inventory certifies it: rim `CIRCLE`s centred on the
axis with axial normals and the radius the semi-angle predicts at that axial
coordinate, a closed degenerate edge at the apex for a pointed cone, and
seam segments whose carrier line passes through the apex; off-law radii,
non-generator seams, missing caps and unclosed rim/apex sets reject as
malformed or trimmed. A `TOROIDAL_SURFACE` lateral face plus two `PLANE`
caps imports as a partial-angle torus segment when the inventory certifies
it: rim `CIRCLE`s of the minor radius centred on the tube circle inside
meridian planes containing the axis, optional open seam arcs on surface
parallels joining the rims, and caps each binding closed circles of one
rim in its meridian plane. The cap face orientations choose between the
two complementary sweeps; inconsistent orientations, off-surface parallels
and non-meridian rims reject. Toroidal boundaries (full torus or segment)
are hard inputs for the current tetrahedralizer: `fill=true` can spend a
long recovery before its explicit blocker, so segment checks exercise
`fill=false`. Every STEP
solid or primitive wrapped by a unit-bearing shape representation converts by
its declared SI or conversion-based length factor into millimetres; a solid
shared by incompatible factors rejects as ambiguous and representations
without unit declarations keep raw coordinates.

This operation returns a mesh. It does not populate the Model/OCC registry,
reconstruct analytic CAD Face/Curve parameter maps or provide general
assembly/Boolean topology. A [mapped STEP occurrence producer](step_mapped_occurrences.md)
supports one explicit representation root, rigid axis-placement targets,
[unit-scale Cartesian targets](step_cartesian_rigid_targets.md) and nested maps
of supported source solids. A separate [direct product-usage producer](step_relationship_assemblies.md)
admits complete single-root NAUO/CDSR/IDT assembly paths with actual rigid
frames and [SI source/placement length conversion](step_relationship_si_units.md)
into millimetre output coordinates. Nested product usages and other assembly transformations remain
guarded. Trimmed and periodic
NURBS topology remains a separate implementation track.

Supported planar STEP/IGES Faces can have one outer wire and disjoint inner
wires. STEP uses FACE_OUTER_BOUND when present; otherwise one wire must uniquely
contain the others. IGES144 and MSBO510 retain their literal outer/inner wire
roles; IGES143 preserves multiple declared boundary groups and certifies one
unique containing outer wire. Tessella CDT triangulates the domain without adding unbound vertices.
Crossing, touching, overlap, wrong winding, nested holes, ambiguous outer roles
and nonplanarity reject. The existing closed-shell and positive-fill checks
certify the result; single-wire Faces retain their previous path.

The default hole check imports a static public-OCC-exported STEP square tube
with fill=false. Full functional checks also cover filled analytic volume 0.75
and IGES144/510 layouts. The own POLY_LOOP STEP and MSBO fixtures test supported
mesh-parser records; those two fixtures do not claim native CAD import parity.
IGES141 admits representation types 0/1. IGES142 uses its 3-D model curve
pointer, rather than interpreting its UV pointer as XYZ. UV-only bounded
line/composite wires on explicit PlaneSurface190 Form1 use the declared
parameter-to-model frame, including finite proper rigid124 Form0 surface transform chains;
see [IGES plane charts](iges_plane_charts.md). Implicit108, nonlinear UV and
general curved Face trimming remain unfinished. These field meanings follow the primary
[Boundary schema](https://occt3d.com/dev/doc/refman/html/class_i_g_e_s_geom___boundary.html)
and [CurveOnSurface schema](https://occt3d.com/dev/doc/refman/html/class_i_g_e_s_geom___curve_on_surface.html).
