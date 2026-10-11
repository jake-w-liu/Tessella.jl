# Independent solid mesh imports

`import_step(path)` and `import_iges(path)` import the supported recognized
solid layouts into meshes. Existing single-solid imports return `Mesh`.
Multiple explicit independent source roots return a classified `MixedMesh`.
Use `fill=false` to obtain their closed triangle boundaries without filling.

Each imported solid retains one discrete boundary Surface and Region, separate
node/element UIDs, and source identity in physical names. STEP entity IDs or
IGES directory-entry IDs determine solid ordering. Filled boundary triangles
come from actual tetrahedron incidence; boundary/interior node owners follow
that incidence. Nodes are not welded across solids. IGES primitive and solid
root rigid directory-entry transforms are applied.

Strictly disjoint bounding boxes certify independent solids. Touching or
overlapping boxes are rejected, including physically disjoint solids whose
boxes overlap. Boxes are never used to discover, split or replace source
roots. Point-cloud box imports retain their original single-solid behavior.
Multi-root IGES identities must be positive, odd, unique and consistent with
the directory-entry type. Unsupported geometric entity types are rejected;
supported type-406 property records can accompany the geometry.

This operation returns a mesh. It does not populate the Model/OCC registry,
reconstruct analytic CAD Face/Curve parameter maps or STEP occurrence
transforms, or provide general assembly/Boolean topology. Trimmed and periodic
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
pointer, rather than interpreting its UV pointer as XYZ; UV-only curves reject
until a parameter-to-model mapping is implemented. General curved Face trimming
remains unfinished. These field meanings follow the primary
[Boundary schema](https://occt3d.com/dev/doc/refman/html/class_i_g_e_s_geom___boundary.html)
and [CurveOnSurface schema](https://occt3d.com/dev/doc/refman/html/class_i_g_e_s_geom___curve_on_surface.html).
