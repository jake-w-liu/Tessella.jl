# Atomic merge into the current model

`Tessella.API.merge(filename)` merges a `.geo`, `.msh` or `.pos` file into the current
model and returns `nothing`. Existing models and the current model's slot keep
their identities, so held model views continue to follow that model. MSH merge
retains its name; GEO `SetName` changes it when explicitly requested. The
associated model filename remains unchanged.

[POS list merging](session_pos_import.md) adds supported named list views and
physical times while retaining current geometry, mesh and model identity.
It prepares the complete file and tag allocation before publishing any view.

The current actual cache and caller-owned records form one detached classified
mesh, including classified unused nodes. Supported MSH2/4 ASCII/binary imports
retain their actual UIDs, owners, physical metadata, periodic pairs, partition
parents/memberships, ghost records, and data sections. Only incoming data
sections create new session views; existing mesh data and held views are not
reimported as duplicates. Metadata indices follow actual UID maps. Conflicting
partition namespaces, duplicate UID occurrences, unsupported element layouts,
and existing format/import admission failures reject before publication.

GEO executes against the actual current geometry, mesh, fields and options. It
can reference existing entities and uses their tags for `newp` and related
allocators. Literal unchanged production inputs and Surface constructor vectors
retain their receipts after revalidation. Products changed by imported nodes or
cells lose those receipts rather than receiving invented emission history.
Prior GEO numeric/string variables are not stored in the API session, so their
persistence across independent file calls is outside this contract.

Parsing, options, diagnostics, mesh/lookup construction, and imported view
transport are staged before the current slot changes. Failed preparation keeps
models, current contents, caches, fields, options, view bindings and diagnostics
unchanged. Explicit external GEO actions such as `Save` or `System` cannot be
rolled back by the session transaction. Internal NamedTuple/Dict calls still
delegate to `Base.merge`; no filename method is installed globally in Base.
Geometry-only periodic declarations do not require realized mesh pairs until
their mesh is produced. Supported GEO `Delete Meshes` and `Delete Model` publish
their actual empty state. A script that changes retained data rows while adding
other data sections currently rejects without an explicit section-occurrence
trace; it does not invent which sections are newly imported.

Native read-time lookup caches can omit caller-owned nodes added after an
earlier import, even though bulk node queries and held model data retain them.
A tiny Gmsh 4.15.2 control reproduces `Unknown node 501 in data` on post-merge
view writing and succeeds after explicit `rebuildNodeCache(False)`. This route
retains that lookup behavior; the writer continuity control calls the matching
public `API.mesh.rebuild_node_cache(false)` before writing the held view.

The focused Julia 1.13.1 qualification passed 48/48 assertions: 46 opt-in
functional checks and two no-generation default smoke checks. Actual child exit
was zero in 44.516 seconds including cold precompilation; the functional groups
took 6.5, 0.8, 0.3, 1.2, and 1.0 seconds. Its small public fixtures exercise all
four MSH encodings, a cache plus caller-owned records, held older/current views,
new NodeData/ElementData/ElementNodeData, literal periodic/ghost transport,
seeded GEO entities and fields, one five-node line generation, and atomic
failure controls. Four additional controls exercise periodic geometry and
clearing; two retain a genuine four-node transfinite Quad product. The two
default smoke assertions generate no mesh. Public
Gmsh 4.15.2 controls use caller-owned Line2 records and a GEO existing-Point
reference without generation; they confirm model/view/name/filename behavior.
The [Gmsh API manual](https://www.gmsh.info/doc/texinfo/gmsh.html) defines merging
file data into the current model separately from opening a model.

The preserved first reader failed before any assertions because a ternary
collection in a new loop needed parentheses. The second reader passed its first
five assertions, then hit the SDK-matching retained-cache writer behavior above.
The final fixture uses the public explicit cache rebuild. Review also identified
and corrected two real merge preparation failures: geometry-only periodic
declarations were incorrectly subjected to the export-only realized-link guard,
and valid GEO clearing statements had been assumed to leave a nonempty mesh.
The final controls verify both. The candidate includes the separately reviewed
shared atomic-import homology callback correction; other open defaults remain.
