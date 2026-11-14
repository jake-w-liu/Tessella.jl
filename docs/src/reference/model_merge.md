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
allocators. Unchanged production inputs and Surface constructors retain their ownership
provenance after revalidation. Imported nodes or cells invalidate affected
provenance.
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

Read-time lookup caches can omit caller-owned nodes added after an earlier
import, even though bulk node queries and held model data retain them. If a
held view cannot resolve such a node when writing after a merge, call
`API.mesh.rebuild_node_cache(false)` before writing the view.
