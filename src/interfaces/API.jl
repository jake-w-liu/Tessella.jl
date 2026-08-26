"""
    API

Gmsh-style model/mesh/option façade over Tessella's native kernels, including
explicit planar surface-loop volumes and persistent affine relations between
straight periodic boundary or embedded curves.
Production meshing is never delegated to Gmsh.
"""
module API

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_surface_loop!, add_volume!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
using ..Model: embed!
using ..Model: add_physical_group!, set_periodic!, model_periodic_nodes
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..GeoExec: execute_geo

export initialize, finalize, option, model, mesh, open_geo!

const CURRENT = Ref{Union{Nothing,GeoModel}}(nothing)
const DEFAULT_OPTIONS = Dict{String,Float64}(
    "Mesh.MeshSizeMin"=>0.0,
    "Mesh.MeshSizeMax"=>1.0e22,
    "Mesh.MeshSizeFactor"=>1.0)
const OPTIONS = copy(DEFAULT_OPTIONS)
const LAST_MESH = Ref{Union{Nothing,Mesh}}(nothing)
const STATE_LOCK = ReentrantLock()

"""
    initialize()

Start a fresh process-global API session, discarding any prior model or cached
mesh and restoring supported options to their defaults.
"""
function initialize()
    lock(STATE_LOCK) do
        CURRENT[]=GeoModel()
        LAST_MESH[]=nothing
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
    end
    return nothing
end

"""End the API session, discard its model and cached mesh, and restore option defaults."""
function finalize()
    lock(STATE_LOCK) do
        CURRENT[]=nothing
        LAST_MESH[]=nothing
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
    end
    return nothing
end

function _model_locked()
    CURRENT[]===nothing && throw(ArgumentError("API: call initialize() first"))
    return CURRENT[]
end

function _with_model(f::Function;invalidate::Bool=false)
    return lock(STATE_LOCK) do
        current=_model_locked()
        result=f(current)
        invalidate && (LAST_MESH[]=nothing)
        result
    end
end

function _copy_mesh(mesh::Mesh)
    return Mesh(mesh.coords;segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
end

"""
    option(name) -> Float64
    option(name, value) -> Float64

Get or set one supported process-global mesh option in an initialized session.
`MeshSizeMin` is nonnegative; `MeshSizeMax` and `MeshSizeFactor` are positive;
the minimum may not exceed the maximum. Boolean and nonfinite values are
rejected, and a failed update leaves all options unchanged.
"""
function option(name::AbstractString)
    key=String(name)
    return lock(STATE_LOCK) do
        _model_locked()
        haskey(OPTIONS,key) || throw(ArgumentError("API.option: unknown option $name"))
        OPTIONS[key]
    end
end
function option(name::AbstractString, value::Real)
    key=String(name)
    return lock(STATE_LOCK) do
        _model_locked()
        haskey(OPTIONS,key) || throw(ArgumentError("API.option: unknown option $name"))
        value isa Bool && throw(ArgumentError("API.option: value must not be Bool"))
        v=try
            Float64(value)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("API.option: value must be Float64-representable"))
        end
        isfinite(v) || throw(ArgumentError("API.option: value must be finite"))
        if key=="Mesh.MeshSizeMin"
            v>=0 || throw(ArgumentError("API.option: MeshSizeMin must be nonnegative"))
            v<=OPTIONS["Mesh.MeshSizeMax"] || throw(ArgumentError(
                "API.option: MeshSizeMin must not exceed MeshSizeMax"))
        elseif key=="Mesh.MeshSizeMax"
            v>0 || throw(ArgumentError("API.option: MeshSizeMax must be positive"))
            v>=OPTIONS["Mesh.MeshSizeMin"] || throw(ArgumentError(
                "API.option: MeshSizeMax must not be below MeshSizeMin"))
        else
            v>0 || throw(ArgumentError("API.option: MeshSizeFactor must be positive"))
        end
        OPTIONS[key]=v
        v
    end
end

"""Gmsh-style geometry-model operations for the active [`API`](@ref) session."""
module model
using ..API: _with_model, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..API: add_surface_loop!, add_volume!
using ..API: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!, embed!, add_physical_group!
add_point(x,y,z;tag=0,meshSize=1.0)=_with_model(invalidate=true) do m
    add_point!(m,x,y,z;tag=tag,mesh_size=meshSize)
end
add_line(a,b;tag=0)=_with_model(invalidate=true) do m
    add_line!(m,a,b;tag=tag)
end
add_curve_loop(curves;tag=0)=_with_model(invalidate=true) do m
    add_curve_loop!(m,curves;tag=tag)
end
add_plane_surface(loops;tag=0)=_with_model(invalidate=true) do m
    add_plane_surface!(m,loops;tag=tag)
end
"""Add a connected closed planar-surface shell to the active model."""
add_surface_loop(surfaces;tag=0)=_with_model(invalidate=true) do m
    add_surface_loop!(m,surfaces;tag=tag)
end
"""Add a volume from an exterior surface loop and optional cavity loops."""
add_volume(surface_loops;tag=0)=_with_model(invalidate=true) do m
    add_volume!(m,surface_loops;tag=tag)
end
add_box(x,y,z,dx,dy,dz;tag=0)=_with_model(invalidate=true) do m
    add_box!(m,x,y,z,dx,dy,dz;tag=tag)
end
add_cylinder(x,y,z,dx,dy,dz,r;tag=0)=_with_model(invalidate=true) do m
    add_cylinder!(m,x,y,z,dx,dy,dz,r;tag=tag)
end
add_sphere(x,y,z,r;tag=0)=_with_model(invalidate=true) do m
    add_sphere!(m,x,y,z,r;tag=tag)
end
add_cone(x,y,z,dx,dy,dz,r1,r2;tag=0)=_with_model(invalidate=true) do m
    add_cone!(m,x,y,z,dx,dy,dz,r1,r2;tag=tag)
end
embed(dim,tags,target_dim,target_tag)=_with_model(invalidate=true) do m
    embed!(m,dim,tags,target_dim,target_tag)
end
function _boolean(op,a,b; tag=0)
    return _with_model(invalidate=true) do m
        t=boolean_volumes!(m,op,a,b;tag=tag)
        # OCC cut/fuse/common remove object and tool by default.
        delete!(m.volumes,Int(a));delete!(m.volumes,Int(b))
        t
    end
end
boolean_difference(a,b; tag=0)=_boolean(:difference,a,b; tag=tag)
boolean_union(a,b; tag=0)=_boolean(:union,a,b; tag=tag)
boolean_intersection(a,b; tag=0)=_boolean(:intersection,a,b; tag=tag)
add_physical_group(dim,tags;tag=0,name="")=_with_model(invalidate=true) do m
    add_physical_group!(m,dim,tags;tag=tag,name=name)
end
end

function _generate(dim::Integer)
    dim isa Bool && throw(ArgumentError("API.mesh.generate: dim must not be Bool"))
    dimension=try
        Int(dim)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("API.mesh.generate: dim exceeds the platform Int range"))
    end
    return lock(STATE_LOCK) do
        m=_model_locked()
        generated=if dimension==2
            isempty(m.surfaces) && throw(ArgumentError("API.mesh.generate: no surfaces"))
            length(m.surfaces)==1 || throw(ArgumentError(
                "API.mesh.generate: Mesh 2 with multiple remaining surfaces is a blocker"))
            mesh_model_surface(m,only(keys(m.surfaces)))
        elseif dimension==3
            isempty(m.volumes) && throw(ArgumentError("API.mesh.generate: no volumes"))
            length(m.volumes)==1 || throw(ArgumentError(
                "API.mesh.generate: Mesh 3 with multiple remaining volumes is a blocker — Boolean Delete the operands or mesh a single volume"))
            mesh_model_volume(m,only(keys(m.volumes)))
        else
            throw(ArgumentError("API.mesh.generate: dim must be 2 or 3"))
        end
        LAST_MESH[]=_copy_mesh(generated)
        generated
    end
end

function _get_mesh()
    return lock(STATE_LOCK) do
        _model_locked()
        LAST_MESH[]===nothing && throw(ArgumentError("API.mesh.get: no mesh"))
        _copy_mesh(LAST_MESH[])
    end
end

function _set_periodic(dim,slave_entities,master_entities,affine;atol=1e-12)
    return _with_model(invalidate=true) do current
        set_periodic!(
            current,dim,slave_entities,master_entities,affine;atol=atol)
    end
end

function _get_periodic_nodes(dim,slave_entity)
    return lock(STATE_LOCK) do
        current=_model_locked()
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "API.mesh.get_periodic_nodes: no mesh"))
        model_periodic_nodes(current,cached,dim,slave_entity)
    end
end

"""Gmsh-style mesh generation, retrieval, and periodic-curve operations."""
module mesh
using ..API: _generate,_get_mesh,_set_periodic,_get_periodic_nodes
generate(dim::Integer)=_generate(dim)
get()=_get_mesh()

"""
    set_periodic(dim, slave_entities, master_entities, affine; atol=1e-12)

Store validated straight-curve relations in the active model and invalidate any
cached mesh. `affine` maps each master curve to its corresponding slave curve in
Gmsh row-major 4×4 order. Their endpoints must be disjoint, and both curves in
each pair must belong to the same planar surface, as boundary or embedded curves,
when meshed.
"""
set_periodic(dim,slave_entities,master_entities,affine;atol=1e-12)=
    _set_periodic(dim,slave_entities,master_entities,affine;atol=atol)

"""
    get_periodic_nodes(dim, slave_entity)

Return the master entity, detached slave/master node arrays, and affine transform
for one boundary or embedded curve relation in the cached surface mesh.
"""
get_periodic_nodes(dim,slave_entity)=
    _get_periodic_nodes(dim,slave_entity)
end

"""
    open_geo!(path; mesh_dim=0)

Execute a bounded `.geo` file into the initialized API session. Stored model and
mesh state are detached from the returned execution result.
"""
function open_geo!(path::AbstractString; mesh_dim::Integer=0)
    return lock(STATE_LOCK) do
        _model_locked()
        result=execute_geo(path;mesh_dim=mesh_dim)
        stored_model=deepcopy(result.model)
        stored_mesh=result.mesh===nothing ? nothing : _copy_mesh(result.mesh)
        CURRENT[]=stored_model
        LAST_MESH[]=stored_mesh
        result
    end
end

end # module
