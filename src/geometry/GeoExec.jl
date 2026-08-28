"""
    GeoExec

Execute a bounded subset of Gmsh `.geo`: Point/Line/Line Loop/Plane Surface/
Surface Loop/Volume, Box/Cylinder/Sphere/Cone, Boolean union/difference/intersection,
Translate/Dilate/90°-Rotate of those solids, Point/Line-In-Surface and
Point/Line/Surface-In-Volume
embeddings with nested point/curve sheet constraints, Physical groups, and
Translate/Rotate/Affine periodic straight curves or explicit-volume planar boundary
surfaces with reusable masters and acyclic dependency chains. Numeric parameters,
entity tags, and entity lists use bounded constant-expression evaluation. Numeric
list variables provide zero-based indexing, cardinality, copying, concatenation,
selection, and checked mutation; entity lists can expand whole variables. Mesh 2/3
runs through the native [`Model`](@ref) kernel. Control-flow loops, macros,
extrusions, fillets, and general OCC BREP remain explicit blockers.
"""
module GeoExec

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_surface_loop!, add_volume!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
using ..Model: embed!, translate_volume!, dilate_volume!, rotate_volume!
using ..Model: add_physical_group!, set_periodic!
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..IO: read_geo_params, _GeoNumericContext, _geo_eval_numeric
using ..IO: _geo_split_list, _geo_numeric_list_terms, _geo_numeric_list_values
using ..IO: _geo_positive_gmsh_tag
using ..IO: _geo_signed_gmsh_int_value
using ..IO: _GEO_SIDE_EFFECT_SYMBOLS
using ..IO: _geo_context_set_scalar!, _geo_apply_list_assignment!
using ..IO: _GeoTagAllocatorState, _geo_context_refresh_allocators!
using ..IO: _geo_allocator_observe_statement!
using ..Transform: _affine_coordinate

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
accepted. Every numeric parameter, entity tag, and numeric entity-list entry in a
supported statement accepts finite arithmetic, prior scalar bindings, and pure
numeric functions. Numeric lists support zero-based scalar indexing, `#name[]`
cardinality, bounded ranges, copies, concatenation, whole-list append/removal, and
indexed or selected mutation. Entity-list positions expand whole or selected list
variables as well as constant ranges. Tags
follow Gmsh's truncation toward zero into positive 32-bit values; oriented Curve and
Surface Loop entries may instead be nonzero signed 32-bit values. `Periodic Line`,
`Periodic Curve`, and `Periodic Surface` accept `Translate`, `Rotate`, and
12- or 16-entry `.geo` `Affine` transforms. Curves must be straight and surfaces
must be planar boundaries of one explicit volume. Multiple periodic statements may
reuse a master or form an acyclic master/slave chain. Read-only `newp`, the shared
curve/loop/surface/volume/Physical-group allocator aliases, and `newf` follow the
tracked explicit topology and supported full Box/Cylinder/Sphere/Cone primitives.
`SetMaxTag Point|Curve|Surface|Volume` follows the active factory: Built-in sets the
checked counter, while OpenCASCADE only raises it. Reads use the greatest counter
among activated factories. Later primitive allocation still accounts for occupied
hidden topology. Primitive boundary entities remain implicit in `GeoModel`, so an
explicit modeled subentity may reuse one of their numeric tags.
An allocator read after a topology-changing or untracked declaration is rejected.
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
    context=_GeoNumericContext()
    allocator_state=_GeoTagAllocatorState()
    for line in _geo_exec_statements(path)
        occursin(r"\b(For|While|Macro|Function|If|Extrude|Torus|Fillet|Chamfer|Symmetry)\b",
                 line) && throw(ArgumentError(
            "execute_geo: unsupported statement $(line) — control-flow loops, " *
            "macros, and advanced OCC features are blockers"))
        _geo_context_refresh_allocators!(context,allocator_state)
        _exec_line!(model,line,context)
        _geo_allocator_observe_statement!(
            allocator_state,line,context,"execute_geo")
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

function _geo_periodic_expressions(raw::AbstractString,count::Int,
                                   context::_GeoNumericContext,
                                   caller::AbstractString)
    return _geo_exec_numeric_values(raw,count,context,caller)
end

function _geo_periodic_affine(raw::AbstractString,context::_GeoNumericContext,
                              caller::AbstractString)
    entries=_geo_exec_numeric_values(raw,context,caller)
    length(entries) in (12,16) || throw(ArgumentError(
        "$caller: expected 12 or 16 numeric values after list expansion"))
    return length(entries)==12 ?
        (entries...,0.0,0.0,0.0,1.0) : Tuple(entries)
end

function _geo_exec_numeric_values(raw::AbstractString,
                                  context::_GeoNumericContext,
                                  caller::AbstractString)
    pieces=_geo_split_list("{"*String(raw)*"}",caller)
    terms,total=_geo_numeric_list_terms(pieces,context,caller)
    values=Float64[];sizehint!(values,total)
    for term in terms
        numeric=term.first
        for _ in 1:term.count
            push!(values,numeric)
            numeric+=term.step
        end
    end
    return values
end

function _geo_exec_numeric_values(raw::AbstractString,count::Int,
                                  context::_GeoNumericContext,
                                  caller::AbstractString)
    values=_geo_exec_numeric_values(raw,context,caller)
    length(values)==count || throw(ArgumentError(
        "$caller: expected $count numeric values after range expansion; " *
        "got $(length(values))"))
    return values
end

function _geo_exec_entity_tags(raw::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString;
                               signed::Bool=false,
                               wrap::Bool=true)
    source=wrap ? "{"*String(raw)*"}" : String(raw)
    values=_geo_numeric_list_values(
        source,context,caller;allow_multiplier=true)
    isempty(values) && throw(ArgumentError(
        "$caller: entity list is empty"))
    tags=Int[];sizehint!(tags,length(values))
    for numeric in values
        tag=if signed
            value=_geo_signed_gmsh_int_value(numeric,"$caller entry")
            iszero(value) && throw(ArgumentError(
                "$caller entry must evaluate to a nonzero signed tag"))
            value
        else
            _geo_positive_gmsh_tag(numeric,"$caller entry")
        end
        push!(tags,tag)
    end
    return tags
end

function _geo_exec_entity_tag(raw::AbstractString,
                              context::_GeoNumericContext,
                              caller::AbstractString)
    return _geo_positive_gmsh_tag(
        _geo_eval_numeric(raw,context,caller),caller)
end

function _geo_exec_single_entity(raw::AbstractString,
                                 context::_GeoNumericContext,
                                 caller::AbstractString)
    tags=_geo_exec_entity_tags(raw,context,caller)
    length(tags)==1 || throw(ArgumentError(
        "$caller: expected exactly one entity; got $(length(tags))"))
    return only(tags)
end

function _geo_exec_entity_rhs_tags(raw::AbstractString,
                                   context::_GeoNumericContext,
                                   caller::AbstractString;
                                   signed::Bool=false)
    source=String(strip(raw))
    isempty(source) && throw(ArgumentError("$caller: entity list must not be empty"))
    has_open=startswith(source,"{")
    has_close=endswith(source,"}")
    has_open==has_close || throw(ArgumentError(
        "$caller: malformed brace-delimited entity list"))
    return _geo_exec_entity_tags(
        source,context,caller;signed=signed,wrap=false)
end

_geo_periodic_tags(raw::AbstractString,context::_GeoNumericContext,
                   caller::AbstractString)=
    _geo_exec_entity_tags(raw,context,caller)

function _geo_exec_scalar!(context::_GeoNumericContext,name::AbstractString,
                           raw::AbstractString)
    variable=String(name)
    variable=="Pi" && throw(ArgumentError(
        "execute_geo: Pi is a reserved numeric constant; use a different scalar name"))
    variable in _GEO_SIDE_EFFECT_SYMBOLS && throw(ArgumentError(
        "execute_geo: dynamic tag allocator $variable is read-only"))
    caller="execute_geo: scalar variable $variable"
    value=_geo_eval_numeric(raw,context,caller)
    _geo_context_set_scalar!(context,variable,value,caller)
    return nothing
end

function _geo_periodic_rotation(axis,origin,angle::Float64,
                                caller::AbstractString)
    axis_scale=max(abs(axis[1]),abs(axis[2]),abs(axis[3]))
    axis_scale>0 || throw(ArgumentError(
        "$caller: rotation axis must have positive length"))
    scaled=(axis[1]/axis_scale,axis[2]/axis_scale,axis[3]/axis_scale)
    scaled_length=hypot(scaled...)
    (isfinite(scaled_length) && scaled_length>0) || throw(ArgumentError(
        "$caller: rotation axis is not normalizable"))
    x=scaled[1]/scaled_length;y=scaled[2]/scaled_length
    z=scaled[3]/scaled_length
    sine,cosine=sincos(angle);one_minus=1-cosine
    r11=muladd(x*x,one_minus,cosine)
    r12=muladd(x*y,one_minus,-z*sine)
    r13=muladd(x*z,one_minus,y*sine)
    r21=muladd(y*x,one_minus,z*sine)
    r22=muladd(y*y,one_minus,cosine)
    r23=muladd(y*z,one_minus,-x*sine)
    r31=muladd(z*x,one_minus,-y*sine)
    r32=muladd(z*y,one_minus,x*sine)
    r33=muladd(z*z,one_minus,cosine)
    rows=((r11,r12,r13),(r21,r22,r23),(r31,r32,r33))
    translation=ntuple(3) do row
        coefficients=rows[row]
        _affine_coordinate(
            origin[row],0.0,coefficients[1],coefficients[2],coefficients[3],
            0.0,0.0,0.0,origin[1],origin[2],origin[3],caller,row)
    end
    return (r11,r12,r13,translation[1],
            r21,r22,r23,translation[2],
            r31,r32,r33,translation[3],
            0.0,0.0,0.0,1.0)
end

function _exec_periodic!(m::GeoModel,line::AbstractString,
                         context::_GeoNumericContext)
    entity=match(r"^Periodic\s+([A-Za-z]+)",line)
    entity===nothing && throw(ArgumentError(
        "execute_geo: malformed periodic statement $line"))
    entity_name=entity.captures[1]
    dim=entity_name in ("Line","Curve") ? 1 : entity_name=="Surface" ? 2 :
        throw(ArgumentError(
            "execute_geo: only straight Line/Curve and planar Surface " *
            "periodicity is implemented"))
    caller="execute_geo: Periodic $(dim==1 ? "Curve" : "Surface")"
    statement=match(
        r"^Periodic\s+(?:Line|Curve|Surface)\s*\{\s*([^}]*)\s*\}\s*=\s*\{\s*([^}]*)\s*\}\s*(.*?)\s*;$",
        line)
    statement===nothing && throw(ArgumentError(
        "$caller: malformed periodic statement $line"))
    slaves=_geo_periodic_tags(statement.captures[1],context,caller)
    masters=_geo_periodic_tags(statement.captures[2],context,caller)
    transform=strip(statement.captures[3])
    affine=if (matched=match(r"^Translate\s*\{\s*([^}]*)\s*\}$",transform)) !== nothing
        delta=_geo_periodic_expressions(
            matched.captures[1],3,context,"$caller Translate")
        (1.0,0.0,0.0,delta[1],
         0.0,1.0,0.0,delta[2],
         0.0,0.0,1.0,delta[3],
         0.0,0.0,0.0,1.0)
    elseif (matched=match(r"^Affine\s*\{\s*([^}]*)\s*\}$",transform)) !== nothing
        _geo_periodic_affine(matched.captures[1],context,"$caller Affine")
    elseif (matched=match(
            r"^Rotate\s*\{\s*\{\s*([^}]*)\s*\}\s*,\s*\{\s*([^}]*)\s*\}\s*,\s*([^}]*)\s*\}$",
            transform)) !== nothing
        axis=_geo_periodic_expressions(
            matched.captures[1],3,context,"$caller Rotate axis")
        origin=_geo_periodic_expressions(
            matched.captures[2],3,context,"$caller Rotate center")
        angle=_geo_eval_numeric(
            matched.captures[3],context,"$caller Rotate angle")
        _geo_periodic_rotation(axis,origin,angle,"$caller Rotate")
    else
        throw(ArgumentError(
            "$caller: expected a Translate, Rotate, or Affine transform"))
    end
    set_periodic!(m,dim,slaves,masters,affine)
    return nothing
end

function _exec_line!(m::GeoModel,line::AbstractString,
                     context::_GeoNumericContext)
    if occursin(r"^Periodic(?:\s|$)",line)
        _exec_periodic!(m,line,context)
        return
    elseif (mm=match(
            r"^Point\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Point"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller coordinates")
        length(values) in (3,4) || throw(ArgumentError(
            "$caller coordinates: expected three coordinates and optional " *
            "mesh size; got $(length(values)) values after range expansion"))
        mesh_size=length(values)==4 ? values[4] : 1.0
        add_point!(m,values[1],values[2],values[3];
                   tag=tag,mesh_size=mesh_size)
        return
    elseif (mm=match(
            r"^Line\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Line"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller endpoints")
        length(points)==2 || throw(ArgumentError(
            "$caller: expected two endpoint tags; got $(length(points))"))
        add_line!(m,points[1],points[2];tag=tag)
        return
    elseif (mm=match(
            r"^(?:Line\s+Loop|Curve\s+Loop)\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Curve Loop"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        curves=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller curves";signed=true)
        add_curve_loop!(m,curves;tag=tag)
        return
    elseif (mm=match(
            r"^Plane\s+Surface\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Plane Surface"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        loops=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller loops")
        add_plane_surface!(m,loops;tag=tag)
        return
    elseif (mm=match(
            r"^Surface\s+Loop\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Surface Loop"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        surfaces=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller surfaces";signed=true)
        add_surface_loop!(m,surfaces;tag=tag)
        return
    elseif (mm=match(
            r"^Volume\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Volume"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        shells=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller surface loops")
        add_volume!(m,shells;tag=tag)
        return
    elseif (mm=match(
            r"^Box\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Box"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],6,context,"$caller parameters")
        add_box!(m,values...;tag=tag)
        return
    elseif (mm=match(
            r"^Cylinder\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Cylinder"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],7,context,"$caller parameters")
        add_cylinder!(m,values...;tag=tag)
        return
    elseif (mm=match(
            r"^Sphere\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Sphere"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],4,context,"$caller parameters")
        add_sphere!(m,values...;tag=tag)
        return
    elseif (mm=match(
            r"^Cone\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Cone"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],8,context,"$caller parameters")
        add_cone!(m,values...;tag=tag)
        return
    elseif (mm=match(
            r"^Boolean(Difference|Union|Intersection)\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*Volume\s*\{\s*(.*?)\s*\}([^}]*)\}\s*\{\s*Volume\s*\{\s*(.*?)\s*\}([^}]*)\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Boolean$(mm.captures[1])"
        op=Dict("Difference"=>:difference,"Union"=>:union,"Intersection"=>:intersection)[mm.captures[1]]
        tag=_geo_exec_entity_tag(mm.captures[2],context,"$caller result tag")
        a=_geo_exec_single_entity(
            mm.captures[3],context,"$caller first operand")
        b=_geo_exec_single_entity(
            mm.captures[5],context,"$caller second operand")
        delete_a=_boolean_delete_operand(mm.captures[4])
        delete_b=_boolean_delete_operand(mm.captures[6])
        boolean_volumes!(m,op,a,b;tag=tag)
        delete_a && delete!(m.volumes,a)
        delete_b && delete!(m.volumes,b)
        return
    elseif (mm=match(
            r"^Translate\s*\{\s*(.*?)\s*\}\s*\{\s*Volume\s*\{\s*(.*?)\s*\}\s*;?\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Translate"
        offset=_geo_exec_numeric_values(
            mm.captures[1],3,context,"$caller offset")
        tag=_geo_exec_single_entity(
            mm.captures[2],context,"$caller volume")
        translate_volume!(m,tag,Tuple(offset))
        return
    elseif (mm=match(
            r"^Dilate\s*\{\s*\{\s*(.*?)\s*\}\s*,\s*(.*?)\s*\}\s*\{\s*Volume\s*\{\s*(.*?)\s*\}\s*;?\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Dilate"
        center=_geo_exec_numeric_values(
            mm.captures[1],3,context,"$caller center")
        scale=_geo_eval_numeric(mm.captures[2],context,"$caller scale")
        tag=_geo_exec_single_entity(
            mm.captures[3],context,"$caller volume")
        dilate_volume!(m,tag,Tuple(center),scale)
        return
    elseif (mm=match(
            r"^Rotate\s*\{\s*\{\s*(.*?)\s*\}\s*,\s*\{\s*(.*?)\s*\}\s*,\s*(.*?)\s*\}\s*\{\s*Volume\s*\{\s*(.*?)\s*\}\s*;?\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Rotate"
        axis=_geo_exec_numeric_values(
            mm.captures[1],3,context,"$caller axis")
        origin=_geo_exec_numeric_values(
            mm.captures[2],3,context,"$caller origin")
        angle=_geo_eval_numeric(mm.captures[3],context,"$caller angle")
        tag=_geo_exec_single_entity(
            mm.captures[4],context,"$caller volume")
        rotate_volume!(m,tag,Tuple(axis),Tuple(origin),angle)
        return
    elseif (mm=match(
            r"^(Point|Line|Curve)\s*\{\s*(.*?)\s*\}\s+In\s+Surface\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: $(mm.captures[1]) In Surface"
        dim=mm.captures[1]=="Point" ? 0 : 1
        ids=_geo_exec_entity_tags(mm.captures[2],context,"$caller entities")
        target=_geo_exec_single_entity(
            mm.captures[3],context,"$caller target")
        embed!(m,dim,ids,2,target)
        return
    elseif (mm=match(
            r"^(Point|Line|Curve|Surface)\s*\{\s*(.*?)\s*\}\s+In\s+Volume\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: $(mm.captures[1]) In Volume"
        dim=Dict("Point"=>0,"Line"=>1,"Curve"=>1,"Surface"=>2)[mm.captures[1]]
        ids=_geo_exec_entity_tags(mm.captures[2],context,"$caller entities")
        target=_geo_exec_single_entity(
            mm.captures[3],context,"$caller target")
        embed!(m,dim,ids,3,target)
        return
    elseif (mm=match(
            r"^Physical\s+(Point|Curve|Line|Surface|Volume)\s*\(\s*(?:\"([^\"]*)\"\s*,\s*)?(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Physical $(mm.captures[1])"
        dim=Dict("Point"=>0,"Curve"=>1,"Line"=>1,"Surface"=>2,"Volume"=>3)[mm.captures[1]]
        name=mm.captures[2]===nothing ? "" : mm.captures[2]
        tag=_geo_exec_entity_tag(mm.captures[3],context,"$caller tag")
        ids=_geo_exec_entity_rhs_tags(
            mm.captures[4],context,"$caller entities")
        add_physical_group!(m,dim,ids;tag=tag,name=name)
        return
    elseif match(
            r"^SetMaxTag\s+(?:Point|Curve|Surface|Volume)\s*\(\s*.+\s*\)\s*;$",
            line) !== nothing
        return
    elseif startswith(line,"Mesh.") || startswith(line,"SetFactory") ||
           startswith(line,"Field") || startswith(line,"Background") ||
           startswith(line,"BoundaryLayer") || startswith(line,"Coherence") ||
           occursin(r"^Mesh\s+[0-9]\s*;", line)
        return
    elseif occursin(r"^[A-Za-z_][A-Za-z0-9_]*\s*\[",line)
        body=String(strip(line[firstindex(line):prevind(line,lastindex(line))]))
        _geo_apply_list_assignment!(
            context,body,"execute_geo: numeric list assignment") || throw(
                ArgumentError("execute_geo: malformed numeric list assignment: $line"))
        return
    elseif (mm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*;$",line)) !== nothing
        _geo_exec_scalar!(context,mm.captures[1],mm.captures[2])
        return
    end
    throw(ArgumentError("execute_geo: unrecognized statement: $line"))
end

end # module
