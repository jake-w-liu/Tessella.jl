"""
    API

Gmsh-style model/mesh/option façade over Tessella's native kernels. Production
meshing is never delegated to Gmsh.
"""
module API

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_box!, add_physical_group!, mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..GeoExec: execute_geo

export initialize, finalize, option, model, mesh, open_geo!

const CURRENT = Ref{Union{Nothing,GeoModel}}(nothing)
const OPTIONS = Dict{String,Float64}("Mesh.MeshSizeMin"=>0.0,"Mesh.MeshSizeMax"=>1.0e22,
                                     "Mesh.MeshSizeFactor"=>1.0)
const LAST_MESH = Ref{Union{Nothing,Mesh}}(nothing)

function initialize()
    CURRENT[]=GeoModel()
    LAST_MESH[]=nothing
    return nothing
end

function finalize()
    CURRENT[]=nothing
    LAST_MESH[]=nothing
    return nothing
end

function _model()
    CURRENT[]===nothing && throw(ArgumentError("API: call initialize() first"))
    return CURRENT[]
end

function option(name::AbstractString)
    haskey(OPTIONS,name) || throw(ArgumentError("API.option: unknown option $name"))
    return OPTIONS[name]
end
function option(name::AbstractString, value::Real)
    haskey(OPTIONS,name) || throw(ArgumentError("API.option: unknown option $name"))
    v=Float64(value); isfinite(v) || throw(ArgumentError("API.option: value must be finite"))
    OPTIONS[name]=v
    return v
end

module model
using ..API: _model, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..API: add_box!, add_physical_group!
add_point(x,y,z; tag=0, meshSize=1.0)=add_point!(_model(),x,y,z; tag=tag, mesh_size=meshSize)
add_line(a,b; tag=0)=add_line!(_model(),a,b; tag=tag)
add_curve_loop(curves; tag=0)=add_curve_loop!(_model(),curves; tag=tag)
add_plane_surface(loops; tag=0)=add_plane_surface!(_model(),loops; tag=tag)
add_box(x,y,z,dx,dy,dz; tag=0)=add_box!(_model(),x,y,z,dx,dy,dz; tag=tag)
add_physical_group(dim,tags; tag=0, name="")=add_physical_group!(_model(),dim,tags; tag=tag, name=name)
end

module mesh
using ..API: _model, LAST_MESH, mesh_model_surface, mesh_model_volume
function generate(dim::Integer)
    m=_model()
    if dim==2
        isempty(m.surfaces) && throw(ArgumentError("API.mesh.generate: no surfaces"))
        LAST_MESH[]=mesh_model_surface(m, minimum(keys(m.surfaces)))
    elseif dim==3
        isempty(m.volumes) && throw(ArgumentError("API.mesh.generate: no volumes"))
        LAST_MESH[]=mesh_model_volume(m, minimum(keys(m.volumes)))
    else
        throw(ArgumentError("API.mesh.generate: dim must be 2 or 3"))
    end
    return LAST_MESH[]
end
get() = (LAST_MESH[]===nothing ? throw(ArgumentError("API.mesh.get: no mesh")) : LAST_MESH[])
end

function open_geo!(path::AbstractString; mesh_dim::Integer=0)
    result=execute_geo(path; mesh_dim=mesh_dim)
    CURRENT[]=result.model
    LAST_MESH[]=result.mesh
    return result
end

end # module
