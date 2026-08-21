"""
    GeoExec

Execute a bounded subset of Gmsh `.geo`: Point/Line/Line Loop/Plane Surface,
Box, Physical groups, and Mesh 2/3 via the native [`Model`](@ref) kernel.
Loops, macros, OCC factories, and CSG remain explicit blockers.
"""
module GeoExec

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_box!, add_physical_group!, set_physical_name!
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..IO: read_geo_params

export execute_geo, GeoExecution

struct GeoExecution
    model::GeoModel
    mesh::Union{Nothing,Mesh}
    params
end

function execute_geo(path::AbstractString; mesh_dim::Integer=0)
    isfile(path) || throw(ArgumentError("execute_geo: missing file $path"))
    params=read_geo_params(path)
    model=GeoModel()
    dim=Int(mesh_dim)
    dim in (0,2,3) || throw(ArgumentError("execute_geo: mesh_dim must be 0, 2, or 3"))
    for raw in eachline(path)
        line=strip(first(split(raw,"//";limit=2)))
        isempty(line) && continue
        occursin(r"\b(For|While|Macro|Function|If|Boolean|Extrude|Disk|Cylinder|Sphere|Torus)\b",
                 line) && throw(ArgumentError(
            "execute_geo: unsupported statement $(line) — loops, macros, OCC, and CSG are blockers"))
        _exec_line!(model,line)
    end
    mesh=nothing
    if dim==2
        isempty(model.surfaces) && throw(ArgumentError("execute_geo: Mesh 2 requested but no surfaces exist"))
        tag=minimum(keys(model.surfaces))
        mesh=mesh_model_surface(model,tag)
    elseif dim==3
        isempty(model.volumes) && throw(ArgumentError("execute_geo: Mesh 3 requested but no volumes exist"))
        tag=minimum(keys(model.volumes))
        mesh=mesh_model_volume(model,tag)
    end
    return GeoExecution(model,mesh,params)
end

function _exec_line!(m::GeoModel, line::AbstractString)
    if (mm=match(r"^Point\s*\(\s*([0-9]+)\s*\)\s*=\s*\{\s*([^,]+),\s*([^,]+),\s*([^,}]+)(?:,\s*([^}]+))?\s*\}\s*;$", line)) !== nothing
        tag=parse(Int,mm.captures[1])
        x=parse(Float64,mm.captures[2]); y=parse(Float64,mm.captures[3]); z=parse(Float64,mm.captures[4])
        lc=mm.captures[5]===nothing ? 1.0 : parse(Float64,mm.captures[5])
        add_point!(m,x,y,z; tag=tag, mesh_size=lc)
        return
    elseif (mm=match(r"^Line\s*\(\s*([0-9]+)\s*\)\s*=\s*\{\s*([0-9]+)\s*,\s*([0-9]+)\s*\}\s*;$", line)) !== nothing
        add_line!(m, parse(Int,mm.captures[2]), parse(Int,mm.captures[3]); tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^(?:Line\s+Loop|Curve\s+Loop)\s*\(\s*([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        ids=parse.(Int, split(mm.captures[2],','))
        add_curve_loop!(m, ids; tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^Plane\s+Surface\s*\(\s*([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        ids=parse.(Int, split(mm.captures[2],','))
        add_plane_surface!(m, ids; tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^Box\s*\(\s*([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        nums=parse.(Float64, split(mm.captures[2],','))
        length(nums)==6 || throw(ArgumentError("execute_geo: Box needs six numbers"))
        add_box!(m, nums[1],nums[2],nums[3],nums[4],nums[5],nums[6]; tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^Physical\s+(Point|Curve|Line|Surface|Volume)\s*\(\s*(?:\"([^\"]*)\"\s*,\s*)?([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        dim=Dict("Point"=>0,"Curve"=>1,"Line"=>1,"Surface"=>2,"Volume"=>3)[mm.captures[1]]
        name=mm.captures[2]===nothing ? "" : mm.captures[2]
        tag=parse(Int,mm.captures[3])
        ids=parse.(Int, split(mm.captures[4],','))
        add_physical_group!(m, dim, ids; tag=tag, name=name)
        return
    elseif startswith(line,"Mesh.") || startswith(line,"SetFactory") ||
           startswith(line,"Field") || startswith(line,"Background") ||
           occursin(r"^Mesh\s+[0-9]\s*;", line)
        return
    end
    throw(ArgumentError("execute_geo: unrecognized statement: $line"))
end

end # module
