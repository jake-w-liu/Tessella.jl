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
