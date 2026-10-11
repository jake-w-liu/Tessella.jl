# Atomic model file import

`Tessella.API.open(filename)` is a session importer for `.geo` and `.msh` files.
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

`API.open_geo!(path; mesh_dim=0)` retains its success behavior: it replaces the
current model's contents and returns a detached `GeoExecution`. It now shares
the same preparation/publication boundary. Earlier native failure witnesses
that expected option writes or selector4 Face prefixes after a failed
`open_geo!` describe the previous partial-publication contract. They are not
evidence for the current atomic behavior.

The `open` binding is local to API. Internal view-writer stream calls explicitly
use `Base.open`, and no filename method is added to `Base.open`. `Elements`
readers remain available for detached mesh reads. A session `merge(filename)`
frontend remains separate work; the existing GEO `Merge` statement is supported.

The opt-in regression `test/interfaces/api_model_import_functional_test.jl`
passed 74 checks. A second opt-in regression,
`test/interfaces/api_model_import_generation_functional_test.jl`, passed 11
checks on a real GEO file containing two Points, one transfinite Line, and
`Mesh 1`. It verifies three actual node UIDs, four Point/Line cell UIDs, literal
Line ownership, and detached legacy results. The default smoke adds three
small checks. The final Julia 1.13.1 run passed 88/88 with bounds checking and
`--compile=min -O0` in 43.448 seconds including cold precompilation; the three
test groups took 5.9, 0.6, and 0.1 seconds. The legacy native midpoint and its
detached analytic carrier differ by about 1.3e-12, so their geometry comparison
uses an absolute tolerance of 1e-11 while UID/count/classification checks remain
exact. Earlier midpoint assertion failures were fixture precision errors.
Preparation also performs the ordinary cache publisher's record-mirror refresh
on the detached model before publication, preserving that existing invariant.
These checks use small actual files and caller-owned records with only this
three-node mesh generation. Public Gmsh 4.15.2 controls cover model/filename/view lifetimes,
empty-current reuse, duplicate names, and declared GEO names. The
[Gmsh API manual](https://www.gmsh.info/doc/texinfo/gmsh.html) distinguishes opening
a new model from merging model data into the current model.
