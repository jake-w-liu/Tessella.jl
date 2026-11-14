# Atomic model file import

`Tessella.API.open(filename)` imports GEO, supported mesh/post-processing
formats, and native STEP (`.step`/`.stp`) or IGES (`.igs`/`.iges`) CAD files.
It returns `nothing`. A successful open makes the imported model current while
retaining nonempty prior models and their views. It reuses an empty current
model, including the initial unnamed model. Repeated opens of a nonempty file
can create models with the same name, following the existing duplicate-name
session contract.

```julia
using Tessella
API.initialize()
API.open("model.msh")
API.option("Mesh.SaveAll", 1)
API.write("copy.msh")
API.finalize()
```

The model name defaults to the filename stem. GEO `SetName` overrides it,
including an explicitly empty name. The associated model filename is the
filename passed to `open`. Nonempty older models keep their filenames, public
mesh records, and held model-view bindings when the current model changes.

[POS list imports](session_pos_import.md) publish supported list views
into a real empty current model. They retain names, physical times, text and
bound interpolation schemes; malformed
or unsupported files leave the session unchanged.

STEP and IGES files use the same native importer as
`API.model.occ.import_shapes`, then synchronize their actual Point, Curve,
Surface and Volume entities. They retain the importer's geometry coverage and
diagnostics. No mesh or discrete CAD substitute is invented during open;
ordinary mesh generation remains available afterwards. CAD merge adds these
entities while retaining the current mesh, fields, filename and held view
owners. Parsing, synchronization and bounding-context updates occur on
detached state before publication. File suffix aliases are case insensitive,
and filenames are passed directly to the reader as data.

MSH 2.2 and 4.1 imports accept ASCII and binary files through the existing
literal mesh importer. They retain supported actual public UIDs, entity owners,
physical names, periodic pairs, partition parents/memberships, ghost records,
and mesh data sections. NodeData, ElementData, and ElementNodeData sections
create session views bound directly to the imported model slot. Data values,
components, and time steps follow the literal indexed rows. Implicit nodal data
uses its cell's ordered connectivity; explicit nodal data is admitted when its
declared node order matches that connectivity.

Existing import guards remain: unsupported special cells or repeated occurrence
identities are rejected, as are MSH2 data rows referring to unclassified orphan
nodes. MSH4 can retain classified orphan-node data. The existing binary data
section tag width is signed 32-bit; large ASCII data UIDs and large mesh UIDs
have separate supported paths. These are mesh imports, without a CAD backend
or invented geometric classification. GEO files follow the bounded syntax and
production contexts of `GeoExec.execute_geo`.

Parsing, option events, diagnostics, mesh/classification construction, lookup
validation, and view transport are prepared before publication. Failed
preparation leaves the model list/current model, caches, options, native
bounding-box context, lookup state, and session views unchanged. API diagnostic
events are also staged. This intentionally gives a stricter failure contract
than the observed Gmsh 4.15.2 missing-file open, which creates an empty current
model before returning an error. This session contract does not roll back a GEO
script's explicit external actions such as `Save` or `System`.

`API.open_geo!(path; mesh_dim=0)` replaces the current model's contents and
returns a detached `GeoExecution`. It uses the same atomic publication boundary.

The `open` binding is local to API. Internal view-writer stream calls explicitly
use `Base.open`, and no filename method is added to `Base.open`. `Elements`
readers remain available for detached mesh reads. The separate
[session `merge(filename)` frontend](model_merge.md) stages changes against the
current model while retaining its slot; GEO `Merge` statements remain supported.
