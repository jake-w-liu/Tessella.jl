# Current-model MSH export

`Tessella.API.write("model.msh")` exports the current model without creating a
post-processing view. It returns `nothing` and leaves the live model, mesh
cache, public tags, views, counters, and model filename unchanged.

```julia
using Tessella
API.initialize()
# Populate or generate the current model using the public model/mesh API.
API.option("Mesh.MshFileVersion", 4.1)
API.option("Mesh.Binary", 1)
API.option("Mesh.SaveAll", 1)
API.write("model.msh")
API.finalize()
```

MSH 2.2 and 4.1 support ASCII and binary output. Physical groups select cells
unless `Mesh.SaveAll` is enabled or the model has no physical groups. MSH4 keeps
actual public node and element UIDs and includes unused nodes on selected
physical entities. MSH2 retains referenced nodes and uses dense file numbering;
`Mesh.SaveAll=1` clears element physical associations while keeping physical
names, matching the observed Gmsh 4.15.2 public writer behavior.

The detached export combines the cached mesh with actual caller-owned records.
Retained periodic pairs, partition/ghost records, and mesh data sections are
transported through literal UID maps. MSH4 supports partition and ghost
metadata. Physical selection that would discard retained periodic or data
references is rejected; use `Mesh.SaveAll=1` to keep those records. A declared
periodic Curve or Surface must have an actual retained mesh link. Other file
extensions and unsupported metadata/version combinations are rejected before
the existing atomic writer replaces the destination.

`API.view.write` remains the endpoint for exporting a selected view. Its binary
MSH writer continues to use `Base.write` through an API-local stream method;
the new filename endpoint does not extend `Base.write(::AbstractString)`.

The opt-in regression in `test/interfaces/api_model_export_functional_test.jl`
checks eight format/encoding/selection combinations, actual large UIDs, retained
periodic/partition/ghost/data records, destination atomicity, detached state, and
binary view continuity. With the two small default smoke assertions, the
isolated Julia 1.13.1 run passed 70/70 checks with bounds checking and
`--compile=min -O0` in 5.223 seconds. This uses small caller-owned/imported
meshes; it does not generate a mesh. Public Gmsh 4.15.2 observations cover the
physical-selection and numbering behaviors on a two-Line2 fixture. The public
endpoint is documented in the [Gmsh API manual](https://www.gmsh.info/doc/texinfo/gmsh.html).
