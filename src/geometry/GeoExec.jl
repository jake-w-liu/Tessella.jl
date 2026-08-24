"""
    GeoExec

Execute a bounded subset of Gmsh `.geo`: Point/Line/Line Loop/Plane Surface,
Box/Cylinder/Sphere/Cone, Boolean union/difference/intersection, Translate/Dilate/
90°-Rotate of those solids, Point-In-Surface embeddings, Physical groups, and
Mesh 2/3 via the native [`Model`](@ref) kernel. Loops, macros, extrusions,
fillets, and general OCC BREP remain explicit blockers.
"""
module GeoExec

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
using ..Model: embed!, translate_volume!, dilate_volume!, rotate_volume!
using ..Model: add_physical_group!
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..IO: read_geo_params

export execute_geo, GeoExecution

"""
    GeoExecution

Result of [`execute_geo`](@ref): the populated native `model`, an optional
validated `mesh`, and the sizing/field `params` recovered by
[`Tessella.IO.read_geo_params`](@ref).
"""
struct GeoExecution
    model::GeoModel
    mesh::Union{Nothing,Mesh}
    params
end

const _MAX_GEO_EXEC_STATEMENT_BYTES=1_000_000
const _MAX_GEO_EXEC_STATEMENTS=1_000_000

# Brace-aware statement split so BooleanDifference `{ Volume{1}; Delete; }{...};`
# is one statement. Quoted strings and line/block comments are respected.
function _geo_exec_statements(path::AbstractString)
    statements=String[]
    buf=IOBuffer()
    depth=0
    quote_char='\0'
    block_comment=false
    for raw in eachline(path)
        i=firstindex(raw); last=lastindex(raw)
        while i<=last
            c=raw[i]
            nxt=nextind(raw,i)
            nextc=nxt<=last ? raw[nxt] : '\0'
            if block_comment
                if c=='*' && nextc=='/'
                    block_comment=false
                    i=nextind(raw,nxt)
                    continue
                end
            elseif quote_char=='\0' && c=='/' && nextc=='/'
                break
            elseif quote_char=='\0' && c=='/' && nextc=='*'
                write(buf,' ')
                block_comment=true
                i=nextind(raw,nxt)
                continue
            elseif quote_char!='\0'
                write(buf,c)
                c==quote_char && (quote_char='\0')
            elseif c=='"' || c=='\''
                quote_char=c; write(buf,c)
            elseif c=='{'
                depth+=1; write(buf,c)
            elseif c=='}'
                depth>0 || throw(ArgumentError(
                    "execute_geo: unmatched closing brace"))
                depth-=1; write(buf,c)
            elseif c==';' && depth==0
                write(buf,c)
                statement=strip(String(take!(buf)))
                if !isempty(statement)
                    length(statements)<_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                        "execute_geo: input exceeds $_MAX_GEO_EXEC_STATEMENTS statements"))
                    push!(statements,statement)
                end
            else
                write(buf,c)
            end
            position(buf)<=_MAX_GEO_EXEC_STATEMENT_BYTES || throw(ArgumentError(
                "execute_geo: statement exceeds $_MAX_GEO_EXEC_STATEMENT_BYTES bytes"))
            i=nxt
        end
        if quote_char!='\0'
            write(buf,'\n')
        elseif position(buf)>0
            write(buf,' ')
        end
    end
    quote_char=='\0' || throw(ArgumentError("execute_geo: unterminated quoted string"))
    block_comment && throw(ArgumentError("execute_geo: unterminated block comment"))
    depth==0 || throw(ArgumentError("execute_geo: unmatched opening brace"))
    tail=strip(String(take!(buf)))
    isempty(tail) || throw(ArgumentError("execute_geo: unterminated statement: $tail"))
    return statements
end

"""
    execute_geo(path; mesh_dim=0) -> GeoExecution

Execute Tessella's documented bounded `.geo` subset into a new [`GeoModel`](@ref).
Use `mesh_dim=2` or `3` to mesh the single remaining surface or volume;
`mesh_dim=0` only builds the model. Geometry statements outside the bounded
subset and malformed input raise `ArgumentError` instead of being partially
accepted.
"""
function execute_geo(path::AbstractString; mesh_dim::Integer=0)
    isfile(path) || throw(ArgumentError("execute_geo: missing file $path"))
    mesh_dim isa Bool && throw(ArgumentError("execute_geo: mesh_dim must not be Bool"))
    dim=try
        Int(mesh_dim)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("execute_geo: mesh_dim exceeds the platform Int range"))
    end
    dim in (0,2,3) || throw(ArgumentError("execute_geo: mesh_dim must be 0, 2, or 3"))
    params=read_geo_params(path)
    model=GeoModel()
    for line in _geo_exec_statements(path)
        occursin(r"\b(For|While|Macro|Function|If|Extrude|Torus|Fillet|Chamfer|Symmetry)\b",
                 line) && throw(ArgumentError(
            "execute_geo: unsupported statement $(line) — loops, macros, and advanced OCC features are blockers"))
        _exec_line!(model,line)
    end
    mesh=nothing
    if dim==2
        isempty(model.surfaces) && throw(ArgumentError("execute_geo: Mesh 2 requested but no surfaces exist"))
        length(model.surfaces)==1 || throw(ArgumentError(
            "execute_geo: Mesh 2 with multiple remaining surfaces $(sort(collect(keys(model.surfaces)))) is a blocker"))
        mesh=mesh_model_surface(model, only(keys(model.surfaces)))
    elseif dim==3
        isempty(model.volumes) && throw(ArgumentError("execute_geo: Mesh 3 requested but no volumes exist"))
        length(model.volumes)==1 || throw(ArgumentError(
            "execute_geo: Mesh 3 with multiple remaining volumes $(sort(collect(keys(model.volumes)))) is a blocker — Boolean Delete the operands or mesh a single volume"))
        mesh=mesh_model_volume(model, only(keys(model.volumes)))
    end
    return GeoExecution(model,mesh,params)
end

function _boolean_delete_operand(raw::AbstractString)
    suffix=String(strip(raw))
    match(r"^;?\s*(?:Delete\s*;?)?$",suffix)===nothing && throw(ArgumentError(
        "execute_geo: Boolean operand suffix must contain only optional Delete; got $(repr(suffix))"))
    return occursin(r"\bDelete\b",suffix)
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
    elseif (mm=match(r"^Cylinder\s*\(\s*([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        nums=parse.(Float64, split(mm.captures[2],','))
        length(nums)==7 || throw(ArgumentError("execute_geo: Cylinder needs x,y,z,dx,dy,dz,r"))
        add_cylinder!(m, nums[1],nums[2],nums[3],nums[4],nums[5],nums[6],nums[7];
                      tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^Sphere\s*\(\s*([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        nums=parse.(Float64, split(mm.captures[2],','))
        length(nums)==4 || throw(ArgumentError("execute_geo: Sphere needs x,y,z,r"))
        add_sphere!(m, nums[1],nums[2],nums[3],nums[4]; tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^Cone\s*\(\s*([0-9]+)\s*\)\s*=\s*\{([^}]+)\}\s*;$", line)) !== nothing
        nums=parse.(Float64, split(mm.captures[2],','))
        length(nums)==8 || throw(ArgumentError("execute_geo: Cone needs x,y,z,dx,dy,dz,r1,r2"))
        add_cone!(m, nums[1],nums[2],nums[3],nums[4],nums[5],nums[6],nums[7],nums[8];
                  tag=parse(Int,mm.captures[1]))
        return
    elseif (mm=match(r"^Boolean(Difference|Union|Intersection)\s*\(\s*([0-9]+)\s*\)\s*=\s*\{\s*Volume\{([0-9]+)\}([^}]*)\}\s*\{\s*Volume\{([0-9]+)\}([^}]*)\}\s*;$", line)) !== nothing
        op=Dict("Difference"=>:difference,"Union"=>:union,"Intersection"=>:intersection)[mm.captures[1]]
        a=parse(Int,mm.captures[3]); b=parse(Int,mm.captures[5])
        delete_a=_boolean_delete_operand(mm.captures[4])
        delete_b=_boolean_delete_operand(mm.captures[6])
        boolean_volumes!(m, op, a, b; tag=parse(Int,mm.captures[2]))
        delete_a && delete!(m.volumes,a)
        delete_b && delete!(m.volumes,b)
        return
    elseif (mm=match(r"^Translate\s*\{\s*([^}]+)\s*\}\s*\{\s*Volume\{([0-9]+)\}\s*;?\s*\}\s*;$", line)) !== nothing
        nums=parse.(Float64, split(mm.captures[1],','))
        length(nums)==3 || throw(ArgumentError("execute_geo: Translate needs three offsets"))
        tag=parse(Int,mm.captures[2])
        translate_volume!(m,tag,(nums[1],nums[2],nums[3]))
        return
    elseif (mm=match(r"^Dilate\s*\{\s*\{\s*([^}]+)\s*\}\s*,\s*([^}]+)\s*\}\s*\{\s*Volume\{([0-9]+)\}\s*;?\s*\}\s*;$", line)) !== nothing
        center=parse.(Float64, split(mm.captures[1],','))
        length(center)==3 || throw(ArgumentError("execute_geo: Dilate center needs three coordinates"))
        scale=parse(Float64, mm.captures[2])
        dilate_volume!(m, parse(Int,mm.captures[3]), (center[1],center[2],center[3]), scale)
        return
    elseif (mm=match(r"^Rotate\s*\{\s*\{\s*([^}]+)\s*\}\s*,\s*\{\s*([^}]+)\s*\}\s*,\s*([^}]+)\s*\}\s*\{\s*Volume\{([0-9]+)\}\s*;?\s*\}\s*;$", line)) !== nothing
        axis=parse.(Float64, split(mm.captures[1],','))
        origin=parse.(Float64, split(mm.captures[2],','))
        (length(axis)==3 && length(origin)==3) || throw(ArgumentError(
            "execute_geo: Rotate needs axis and origin triples"))
        angle=parse(Float64, mm.captures[3])
        rotate_volume!(m, parse(Int,mm.captures[4]),
                       (axis[1],axis[2],axis[3]), (origin[1],origin[2],origin[3]), angle)
        return
    elseif (mm=match(r"^(Point|Line|Curve)\s*\{\s*([^}]+)\s*\}\s+In\s+Surface\s*\{\s*([0-9]+)\s*\}\s*;$", line)) !== nothing
        dim=mm.captures[1]=="Point" ? 0 : 1
        ids=parse.(Int, split(mm.captures[2],','))
        embed!(m, dim, ids, 2, parse(Int,mm.captures[3]))
        return
    elseif (mm=match(r"^(Point|Line|Curve|Surface)\s*\{\s*([^}]+)\s*\}\s+In\s+Volume\s*\{\s*([0-9]+)\s*\}\s*;$", line)) !== nothing
        dim=Dict("Point"=>0,"Line"=>1,"Curve"=>1,"Surface"=>2)[mm.captures[1]]
        ids=parse.(Int, split(mm.captures[2],','))
        embed!(m, dim, ids, 3, parse(Int,mm.captures[3]))
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
           startswith(line,"Coherence") || occursin(r"^Mesh\s+[0-9]\s*;", line)
        return
    end
    throw(ArgumentError("execute_geo: unrecognized statement: $line"))
end

end # module
