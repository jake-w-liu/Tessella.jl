"""
    GeoExec

Execute a bounded subset of Gmsh `.geo`: Point/Line/Line Loop/Plane Surface/
Surface Loop/Volume, Box/Cylinder/Sphere/Cone, Boolean union/difference/intersection,
Translate/Dilate/90°-Rotate of those solids, Point/Line-In-Surface and
Point/Line/Surface-In-Volume
embeddings with nested point/curve sheet constraints, Physical groups, and
Translate/Rotate/Affine periodic straight curves, explicit-volume planar boundary
surfaces, or stored mesh-inert volume relations with reusable masters and acyclic
dependency chains. `MeshSize` and
`Characteristic Length` update existing explicit Point constraints directly or
through recursive `PointsOf` boundaries of explicit Point/Curve/Surface/Volume
entities. Physical declarations accept an explicit tag, with an optional name, or
a nonempty name with a tag allocated from the global Physical-group namespace. Their
entity lists accept numeric selectors, `PointsOf` for Physical Point, and immediate
`Boundary` or `CombinedBoundary` queries over explicit entities of the next higher
dimension.
Numeric parameters, entity tags, and entity lists use bounded
constant-expression evaluation with Gmsh's comparison, logical, and ternary
operators. `If`/`ElseIf`/`Else`/`EndIf` and `For name In {a:b[:c]}` follow the
built-in kernel (the For range is evaluated once at loop entry and the loop
variable is left at its first out-of-range value); a bounded
`While (expr) ... EndWhile` is a Tessella extension — pinned Gmsh has no While
keyword. Control headers must complete on one line. Numeric
list variables provide zero-based indexing, cardinality, copying, concatenation,
selection, and checked mutation; entity lists can expand whole variables. Mesh 2/3
runs through the native [`Model`](@ref) kernel. Boolean results own operation-time
operand geometry. `Delete` also removes the operand's native encoding, target
embeddings, and Physical memberships; a group and its name are removed when no
members remain. Macros, extrusions, fillets, and general OCC
BREP remain explicit blockers.
"""
module GeoExec

using ..Model: GeoModel, add_point!, set_point_mesh_size!
using ..Model: add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_surface_loop!, add_volume!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
using ..Model: _remove_volume_entity!
using ..Model: embed!, translate_volume!, dilate_volume!, rotate_volume!
using ..Model: add_physical_group!, set_periodic!, set_transfinite_tri!
using ..Model: _model_boundary, _model_points_of
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..IO: read_geo_params, _GeoNumericContext, _geo_eval_numeric
using ..IO: _geo_split_list, _geo_split_range, _geo_range_count
using ..IO: _geo_numeric_list_terms, _geo_numeric_list_values
using ..IO: _geo_positive_gmsh_tag
using ..IO: _geo_signed_gmsh_int_value
using ..IO: _GEO_SIDE_EFFECT_SYMBOLS
using ..IO: _MAX_GEO_LIST_ITEMS
using ..IO: _geo_context_set_scalar!, _geo_apply_list_assignment!
using ..IO: _GeoTagAllocatorState, _geo_context_refresh_allocators!
using ..IO: _geo_allocator_observe_statement!
using ..IO: _geo_physical_declaration
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
    transfinite_tri::Union{Nothing,Int}
end

const _MAX_GEO_EXEC_STATEMENT_BYTES=1_000_000
const _MAX_GEO_EXEC_STATEMENTS=1_000_000
const _MAX_GEO_LOOP_ITERATIONS=1_000_000

const _GEO_CONTROL_BARE=Dict(
    "Else"=>:else,"EndIf"=>:endif,"EndWhile"=>:endwhile,"EndFor"=>:endfor)

# Scan a balanced open/close group starting at raw[i] (the opener). Returns the
# index of the matching close on the same line; control headers that spill to a
# second line are a bounded-subset blocker.
function _geo_scan_balanced(raw::AbstractString,i::Int,last::Int,
                            open::Char,close::Char)
    depth=0;qc='\0';k=i
    while k<=last
        c=raw[k]
        if qc!='\0'
            c==qc && (qc='\0')
        elseif c=='"' || c=='\''
            qc=c
        elseif c==open
            depth+=1
        elseif c==close
            depth-=1
            depth==0 && return k
        end
        k=nextind(raw,k)
    end
    throw(ArgumentError(
        "execute_geo: control statement must complete on a single line"))
end

# If raw[i:last] begins a control construct (If/ElseIf/While headers with a
# parenthesized expression, `For name In {a:b[:c]}`, or the bare Else/EndIf/
# EndWhile/EndFor markers), return (statement_text, last_consumed_index).
# Gmsh 4.15.2's built-in kernel accepts only `For name In {a:b[:c]}` — comma
# lists, single values, and C-style `For (init; cond; incr)` are rejected —
# and has no While keyword; `While (expr) ... EndWhile` is a bounded
# Tessella extension.
function _geo_control_statement(raw::AbstractString,i::Int,last::Int)
    rest=SubString(raw,i,last)
    matched=match(
        r"^(If|ElseIf|While|Else|EndIf|EndWhile|For|EndFor)\b",rest)
    matched===nothing && return nothing
    word=matched.captures[1]
    haskey(_GEO_CONTROL_BARE,word) && return (word,i+sizeof(word)-1)
    j=firstindex(rest)+sizeof(word)
    jlast=lastindex(rest)
    while j<=jlast && isspace(rest[j])
        j=nextind(rest,j)
    end
    if word=="For"
        name_match=match(r"^[A-Za-z_][A-Za-z0-9_]*",SubString(rest,j))
        name_match===nothing && throw(ArgumentError(
            "execute_geo: malformed For header; use `For name In {a:b[:c]}`"))
        j+=sizeof(name_match.match)
        while j<=jlast && isspace(rest[j])
            j=nextind(rest,j)
        end
        match(r"^In\b",SubString(rest,j))===nothing && throw(ArgumentError(
            "execute_geo: malformed For header; use `For name In {a:b[:c]}`"))
        j+=2
        while j<=jlast && isspace(rest[j])
            j=nextind(rest,j)
        end
        (j<=jlast && rest[j]=='{') || throw(ArgumentError(
            "execute_geo: malformed For header; use `For name In {a:b[:c]}`"))
        close=_geo_scan_balanced(rest,j,jlast,'{','}')
        return (String(rest[1:close]),i+close-1)
    end
    (j<=jlast && rest[j]=='(') || throw(ArgumentError(
        "execute_geo: malformed $word header; use `$word (expression)`"))
    close=_geo_scan_balanced(rest,j,jlast,'(',')')
    return (String(rest[1:close]),i+close-1)
end

# Classify a standalone control statement; returns nothing for ordinary
# `;`-terminated statements.
function _geo_control_parse(line::AbstractString)
    if (matched=match(r"^(If|ElseIf|While)\s*\((.*)\)$",line))!==nothing
        kind=Dict("If"=>:if,"ElseIf"=>:elseif,"While"=>:while)[matched.captures[1]]
        return (kind=kind,cond=String(matched.captures[2]))
    elseif (matched=match(
            r"^For\s+([A-Za-z_][A-Za-z0-9_]*)\s+In\s*\{(.*)\}$",line))!==nothing
        return (kind=:for,var=String(matched.captures[1]),
                range=String(matched.captures[2]))
    elseif haskey(_GEO_CONTROL_BARE,line)
        return (kind=_GEO_CONTROL_BARE[line],)
    end
    return nothing
end

# Brace-aware statement split so BooleanDifference `{ Volume{1}; Delete; }{...};`
# is one statement. Quoted strings and line/block comments are respected.
# Control constructs (If/ElseIf/Else/EndIf, For/EndFor, While/EndWhile) carry no
# `;`; each is emitted as its own statement whenever it starts at a statement
# boundary.
function _geo_exec_statements(path::AbstractString)
    statements=String[]
    buf=IOBuffer()
    depth=0
    quote_char='\0'
    block_comment=false
    buf_has_content=false
    for raw in eachline(path)
        i=firstindex(raw); last=lastindex(raw)
        while i<=last
            if !block_comment && quote_char=='\0' && !buf_has_content
                control=_geo_control_statement(raw,i,last)
                if control!==nothing
                    take!(buf)
                    length(statements)<_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                        "execute_geo: input exceeds $_MAX_GEO_EXEC_STATEMENTS statements"))
                    push!(statements,control[1])
                    i=nextind(raw,control[2])
                    continue
                end
            end
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
                quote_char=c; write(buf,c); buf_has_content=true
            elseif c=='{'
                depth+=1; write(buf,c); buf_has_content=true
            elseif c=='}'
                depth>0 || throw(ArgumentError(
                    "execute_geo: unmatched closing brace"))
                depth-=1; write(buf,c); buf_has_content=true
            elseif c==';' && depth==0
                write(buf,c)
                statement=strip(String(take!(buf)))
                if !isempty(statement)
                    length(statements)<_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                        "execute_geo: input exceeds $_MAX_GEO_EXEC_STATEMENTS statements"))
                    push!(statements,statement)
                end
                buf_has_content=false
            else
                write(buf,c)
                isspace(c) || (buf_has_content=true)
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
supported statement accepts finite arithmetic, prior scalar bindings, pure
numeric functions, and Gmsh's comparison, logical, and ternary operators.
`If (expr) ... ElseIf (expr) ... Else ... EndIf` and
`For name In {start:end[:increment]} ... EndFor` follow the built-in kernel: the
For range is a single colon term evaluated once at loop entry (the implicit
two-term increment is +1, so `{5:0}` iterates zero times), the loop variable is
assigned per iteration and remains at its first out-of-range value, and control
headers must complete on one line. `While (expr) ... EndWhile` is a bounded
Tessella extension — pinned Gmsh has no While keyword — capped at
$_MAX_GEO_LOOP_ITERATIONS iterations; total executed statements stay within
$_MAX_GEO_EXEC_STATEMENTS. Numeric lists support zero-based scalar indexing, `#name[]`
cardinality, bounded ranges, copies, concatenation, whole-list append/removal, and
indexed or selected mutation. Entity-list positions expand whole or selected list
variables as well as constant ranges. Tags
follow Gmsh's truncation toward zero into positive 32-bit values; oriented Curve and
Surface Loop entries may instead be nonzero signed 32-bit values. `Periodic Line`,
`Periodic Curve`, `Periodic Surface`, and `Periodic Volume` accept `Translate`,
`Rotate`, and
12- or 16-entry `.geo` `Affine` transforms. Curves must be straight and surfaces
must be planar boundaries of one explicit volume; volume relations are stored
but mesh-inert, as in Gmsh 4.15.2. Multiple periodic statements may
reuse a master or form an acyclic master/slave chain. Read-only `newp`, the shared
curve/loop/surface/volume/Physical-group allocator aliases, and `newf` follow the
tracked explicit topology and supported full Box/Cylinder/Sphere/Cone primitives.
`SetMaxTag Point|Curve|Surface|Volume` follows the active factory: Built-in sets the
checked counter, while OpenCASCADE only raises it. Reads use the greatest counter
among activated factories. Later primitive allocation still accounts for occupied
hidden topology. `Box` materializes Gmsh's exact boundary entities (eight
points, twelve curves, six surfaces) on the model; Cylinder/Sphere/Cone
boundary entities remain implicit, so an explicit modeled subentity may reuse
one of their numeric tags. `Mesh.TransfiniteTri = 0|1` sets the model's
three-sided transfinite surface algorithm like Gmsh's option of the same name.
`MeshSize {points} = value` and its `Characteristic Length` alias update existing
explicit Points through `:`, numeric expressions and ranges, numeric-list variables,
or inline `PointsOf` blocks for Point, Curve/Line, Surface, and explicit Volume
entities. `PointsOf` recursively follows entity boundaries, including hole loops;
embedded entities are not boundaries. Their values must be finite and positive.
Physical declarations accept an explicit positive tag, with an optional name, or a
nonempty name with an automatic tag. Automatic Physical tags share one namespace
across entity dimensions and contribute to later shared-region allocator reads.
Physical Point accepts inline `PointsOf`; Physical Point/Curve/Surface accept inline
`Boundary` and `CombinedBoundary` over Curve/Line, Surface, and explicit Volume
entities, respectively. `Boundary` joins the selected immediate boundaries and
physical-group insertion removes duplicate tags. `CombinedBoundary` retains tags
with odd multiplicity. Hole and cavity boundaries participate; embeddings do not.
Unknown entities, empty combined boundaries, implicit primitive or Boolean volume
topology, unsupported geometry-derived selectors, and invalid query dimensions are
explicit blockers.
An allocator read after an untracked topology-changing declaration is rejected;
tracked `Boolean` operand `Delete` and `SetMaxTag` counters stay live.
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
    statements=_geo_exec_statements(path)
    executed=Ref(0)
    transfinite_tri=_exec_geo_statements!(
        model,statements,firstindex(statements),lastindex(statements),
        context,allocator_state,executed)
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
    return GeoExecution(model,mesh,params,transfinite_tri)
end

const _GEO_CONTROL_CLOSE=Dict(:if=>:endif,:for=>:endfor,:while=>:endwhile)
const _GEO_CONTROL_CLOSE_NAME=Dict(
    :endif=>"EndIf",:endfor=>"EndFor",:endwhile=>"EndWhile")
const _GEO_CONTROL_OPEN_NAME=Dict(
    :if=>"If",:for=>"For",:while=>"While")

# Locate the closer matching the opener statements[i] within
# statements[i+1:hi], honoring nested blocks of any family. For an If block the
# collected depth-1 ElseIf/Else midpoints are returned as well.
function _geo_exec_find_block_end(statements::Vector{String},i::Int,hi::Int)
    opener=_geo_control_parse(statements[i])
    closer=_GEO_CONTROL_CLOSE[opener.kind]
    depth=1;mids=Int[];mid_kinds=Symbol[]
    j=i+1
    while j<=hi
        control=_geo_control_parse(statements[j])
        if control!==nothing
            if haskey(_GEO_CONTROL_CLOSE,control.kind)
                depth+=1
            elseif haskey(_GEO_CONTROL_CLOSE_NAME,control.kind)
                depth-=1
                if depth==0
                    control.kind==closer || throw(ArgumentError(
                        "execute_geo: $(statements[j]) cannot close a " *
                        "$(_GEO_CONTROL_OPEN_NAME[opener.kind]) block"))
                    return j,mids,mid_kinds
                end
            elseif depth==1 && control.kind in (:elseif,:else)
                opener.kind===:if || throw(ArgumentError(
                    "execute_geo: $(statements[j]) cannot appear inside a " *
                    "$(_GEO_CONTROL_OPEN_NAME[opener.kind]) block"))
                push!(mids,j);push!(mid_kinds,control.kind)
            end
        end
        j+=1
    end
    throw(ArgumentError(
        "execute_geo: $(statements[i]) has no matching " *
        "$(_GEO_CONTROL_CLOSE_NAME[closer])"))
end

function _geo_exec_if!(m::GeoModel,statements::Vector{String},i::Int,hi::Int,
                       context::_GeoNumericContext,
                       allocator_state::_GeoTagAllocatorState,
                       executed::Base.RefValue{Int})
    opener=_geo_control_parse(statements[i])
    close,mids,mid_kinds=_geo_exec_find_block_end(statements,i,hi)
    caller="execute_geo: If"
    # Validate the whole branch structure up front, as Gmsh's parser does: a
    # malformed trailing branch is an error even when an earlier branch ran.
    seen_else=false
    for kind in mid_kinds
        if kind===:else
            seen_else && throw(ArgumentError(
                "execute_geo: If block has more than one Else"))
            seen_else=true
        else
            seen_else && throw(ArgumentError(
                "execute_geo: ElseIf cannot follow Else in an If block"))
        end
    end
    branch_lo=i+1
    pending_cond=opener.cond
    for k in eachindex(mids)
        if pending_cond!==nothing &&
           _geo_eval_numeric(pending_cond,context,caller)!=0
            return _exec_geo_statements!(
                m,statements,branch_lo,mids[k]-1,context,allocator_state,executed),
                close+1
        end
        pending_cond=mid_kinds[k]===:else ? nothing :
            _geo_control_parse(statements[mids[k]]).cond
        branch_lo=mids[k]+1
    end
    if seen_else || (pending_cond!==nothing &&
                     _geo_eval_numeric(pending_cond,context,caller)!=0)
        return _exec_geo_statements!(
            m,statements,branch_lo,close-1,context,allocator_state,executed),
            close+1
    end
    return nothing,close+1
end

function _geo_exec_for!(m::GeoModel,statements::Vector{String},i::Int,hi::Int,
                        context::_GeoNumericContext,
                        allocator_state::_GeoTagAllocatorState,
                        executed::Base.RefValue{Int})
    opener=_geo_control_parse(statements[i])
    close,=_geo_exec_find_block_end(statements,i,hi)
    caller="execute_geo: For"
    variable=opener.var
    (variable=="Pi" || variable in _GEO_SIDE_EFFECT_SYMBOLS) && throw(ArgumentError(
        "$caller: loop variable $variable is reserved"))
    pieces=_geo_split_list("{$(opener.range)}",caller)
    range_pieces=length(pieces)==1 ?
        _geo_split_range(pieces[1],caller) : nothing
    range_pieces===nothing && throw(ArgumentError(
        "$caller: Gmsh For ranges require a single `start:end[:increment]` " *
        "term — comma lists and single values are rejected"))
    # The For range is evaluated once at loop entry. Unlike a `a:b` list
    # expression, the implicit two-term increment is always +1, so `{5:0}` is
    # an empty loop (verified against pinned Gmsh 4.15.2's unrolled output).
    first=_geo_eval_numeric(range_pieces[1],context,"$caller range start")
    last=_geo_eval_numeric(range_pieces[2],context,"$caller range end")
    step=length(range_pieces)==2 ? 1.0 :
         _geo_eval_numeric(range_pieces[3],context,"$caller range increment")
    count=_geo_range_count(first,last,step,_MAX_GEO_LIST_ITEMS,caller)
    transfinite_tri=nothing
    value=first
    for _ in 1:count
        _geo_context_set_scalar!(context,variable,value,caller)
        assigned=_exec_geo_statements!(
            m,statements,i+1,close-1,context,allocator_state,executed)
        assigned===nothing || (transfinite_tri=assigned)
        value+=step
    end
    # Gmsh's unroller leaves the loop variable at the first out-of-range value
    # (i = start + count*step), and at the start value for an empty range.
    _geo_context_set_scalar!(context,variable,value,caller)
    return transfinite_tri,close+1
end

function _geo_exec_while!(m::GeoModel,statements::Vector{String},i::Int,hi::Int,
                          context::_GeoNumericContext,
                          allocator_state::_GeoTagAllocatorState,
                          executed::Base.RefValue{Int})
    opener=_geo_control_parse(statements[i])
    close,=_geo_exec_find_block_end(statements,i,hi)
    caller="execute_geo: While"
    transfinite_tri=nothing
    iterations=0
    while _geo_eval_numeric(opener.cond,context,caller)!=0
        iterations+=1
        iterations<=_MAX_GEO_LOOP_ITERATIONS || throw(ArgumentError(
            "$caller: loop exceeds $_MAX_GEO_LOOP_ITERATIONS iterations"))
        assigned=_exec_geo_statements!(
            m,statements,i+1,close-1,context,allocator_state,executed)
        assigned===nothing || (transfinite_tri=assigned)
    end
    return transfinite_tri,close+1
end

# Execute statements[lo:hi]; returns the last Mesh.TransfiniteTri assignment, if
# any. Control constructs recurse into their block ranges; If evaluates
# constant-expression conditions, For iterates `name In {a:b[:c]}` evaluated
# once at loop entry, and While (a bounded Tessella extension) re-evaluates its
# condition each iteration.
function _exec_geo_statements!(m::GeoModel,statements::Vector{String},
                               lo::Int,hi::Int,
                               context::_GeoNumericContext,
                               allocator_state::_GeoTagAllocatorState,
                               executed::Base.RefValue{Int})
    transfinite_tri=nothing
    i=lo
    while i<=hi
        line=statements[i]
        control=_geo_control_parse(line)
        if control===nothing
            executed[]+=1
            executed[]<=_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                "execute_geo: control flow exceeds $_MAX_GEO_EXEC_STATEMENTS " *
                "executed statements"))
            occursin(
                r"\b(Macro|Function|Extrude|Torus|Fillet|Chamfer|Symmetry)\b",
                line) && throw(ArgumentError(
                "execute_geo: unsupported statement $(line) — macros, " *
                "extrusions, and advanced OCC features are blockers"))
            _geo_context_refresh_allocators!(context,allocator_state)
            assigned=_exec_line!(m,line,context)
            assigned===nothing || (transfinite_tri=assigned)
            _geo_allocator_observe_statement!(
                allocator_state,line,context,"execute_geo")
            i+=1
        elseif control.kind===:if
            assigned,i=_geo_exec_if!(
                m,statements,i,hi,context,allocator_state,executed)
            assigned===nothing || (transfinite_tri=assigned)
        elseif control.kind===:for
            assigned,i=_geo_exec_for!(
                m,statements,i,hi,context,allocator_state,executed)
            assigned===nothing || (transfinite_tri=assigned)
        elseif control.kind===:while
            assigned,i=_geo_exec_while!(
                m,statements,i,hi,context,allocator_state,executed)
            assigned===nothing || (transfinite_tri=assigned)
        else
            throw(ArgumentError(
                "execute_geo: $line without a matching opener"))
        end
    end
    return transfinite_tri
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

function _geo_exec_topology_query_payload(selector::AbstractString,
                                          query::AbstractString,
                                          caller::AbstractString)
    source=String(strip(selector))
    matched=match(r"^([A-Za-z][A-Za-z0-9]*)\s*\{",source)
    (matched!==nothing && matched.captures[1]==query) || throw(ArgumentError(
        "$caller: malformed $query selector"))
    opening=findfirst(==('{'),source)
    opening===nothing && throw(ArgumentError(
        "$caller: malformed $query selector"))
    depth=0;closing=nothing;i=opening
    while i<=lastindex(source)
        character=source[i]
        if character=='{'
            depth+=1
        elseif character=='}'
            depth-=1
            depth>=0 || throw(ArgumentError(
                "$caller: $query selector has an unmatched closing brace"))
            if depth==0
                closing=i
                break
            end
        end
        i=nextind(source,i)
    end
    closing===nothing && throw(ArgumentError(
        "$caller: $query selector has an unmatched opening brace"))
    after=nextind(source,closing)
    trailing=after>lastindex(source) ? "" : String(strip(source[after:end]))
    isempty(trailing) || throw(ArgumentError(
        "$caller: unexpected text after the $query selector"))
    first_payload=nextind(source,opening)
    payload=first_payload==closing ? "" :
        String(source[first_payload:prevind(source,closing)])
    return payload
end

function _geo_exec_topology_query_blocks(payload::AbstractString,
                                         query::AbstractString,
                                         caller::AbstractString)
    blocks=String[];buffer=IOBuffer();depth=0
    for character in payload
        if character=='{'
            depth+=1
            write(buffer,character)
        elseif character=='}'
            depth-=1
            depth>=0 || throw(ArgumentError(
                "$caller: $query entity blocks have an unmatched closing brace"))
            write(buffer,character)
        elseif character==';' && depth==0
            block=String(strip(String(take!(buffer))))
            isempty(block) && throw(ArgumentError(
                "$caller: $query contains an empty entity block"))
            length(blocks)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: $query exceeds $_MAX_GEO_LIST_ITEMS entity blocks"))
            push!(blocks,block)
        else
            write(buffer,character)
        end
    end
    depth==0 || throw(ArgumentError(
        "$caller: $query entity blocks have an unmatched opening brace"))
    tail=String(strip(String(take!(buffer))))
    isempty(tail) || throw(ArgumentError(
        "$caller: every $query entity block must end with a semicolon"))
    isempty(blocks) && throw(ArgumentError(
        "$caller: $query must contain at least one entity block"))
    return blocks
end

function _geo_exec_all_entity_tags(m::GeoModel,dimension::Int,
                                   caller::AbstractString)
    entities=dimension==0 ? m.points : dimension==1 ? m.curves :
             dimension==2 ? m.surfaces : m.volumes
    length(entities)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
        "$caller: entity wildcard expands beyond $_MAX_GEO_LIST_ITEMS entities"))
    return sort!(collect(keys(entities)))
end

function _geo_exec_points_of(m::GeoModel,selector::AbstractString,
                             context::_GeoNumericContext,
                             caller::AbstractString)
    query="PointsOf"
    payload=_geo_exec_topology_query_payload(selector,query,caller)
    blocks=_geo_exec_topology_query_blocks(payload,query,caller)
    entities=NTuple{2,Int}[]
    for block in blocks
        matched=match(
            r"^(Point|Curve|Line|Surface|Volume)\s*\{\s*(.*)\s*\}$",block)
        matched===nothing && throw(ArgumentError(
            "$caller: PointsOf supports Point, Curve/Line, Surface, and Volume blocks"))
        kind=String(matched.captures[1])
        dimension=kind=="Point" ? 0 : kind in ("Curve","Line") ? 1 :
                  kind=="Surface" ? 2 : 3
        raw_tags=String(strip(matched.captures[2]))
        tags=if raw_tags==":"
            _geo_exec_all_entity_tags(m,dimension,"$caller PointsOf $kind selector")
        else
            signed_tags=_geo_exec_entity_tags(
                raw_tags,context,"$caller PointsOf $kind selector";signed=true)
            normalized_tags=Int[];sizehint!(normalized_tags,length(signed_tags))
            for tag in signed_tags
                tag==typemin(Int32) && throw(ArgumentError(
                    "$caller PointsOf $kind selector entry magnitude exceeds Int32"))
                push!(normalized_tags,abs(tag))
            end
            normalized_tags
        end
        length(entities)<=_MAX_GEO_LIST_ITEMS-length(tags) || throw(ArgumentError(
            "$caller: PointsOf expands beyond $_MAX_GEO_LIST_ITEMS entities"))
        append!(entities,((dimension,tag) for tag in tags))
    end
    return _model_points_of(m,entities,caller)
end

function _geo_exec_physical_topology(m::GeoModel,selector::AbstractString,
                                     context::_GeoNumericContext,
                                     caller::AbstractString,
                                     physical_dimension::Int)
    source=String(strip(selector))
    if match(r"^PointsOf(?:\s|\{)",source)!==nothing
        physical_dimension==0 || throw(ArgumentError(
            "$caller: PointsOf produces Point entities and can only populate " *
            "a Physical Point group"))
        return _geo_exec_points_of(m,source,context,caller)
    end

    matched=match(r"^(Boundary|CombinedBoundary)(?:\s|\{)",source)
    matched===nothing && throw(ArgumentError(
        "$caller: unsupported topology query; use inline Boundary, " *
        "CombinedBoundary, or PointsOf for Physical Point"))
    query=String(matched.captures[1])
    physical_dimension<3 || throw(ArgumentError(
        "$caller: $query cannot populate a Physical Volume group because " *
        "there is no dimension-4 entity boundary"))
    source_dimension=physical_dimension+1
    source_kind=source_dimension==1 ? "Curve/Line" :
                source_dimension==2 ? "Surface" : "Volume"
    payload=_geo_exec_topology_query_payload(source,query,caller)
    blocks=_geo_exec_topology_query_blocks(payload,query,caller)
    entities=NTuple{2,Int}[]
    for block in blocks
        entity_match=match(
            r"^(Curve|Line|Surface|Volume)\s*\{\s*(.*)\s*\}$",block)
        entity_match===nothing && throw(ArgumentError(
            "$caller: $query for this Physical group requires $source_kind blocks"))
        kind=String(entity_match.captures[1])
        dimension=kind in ("Curve","Line") ? 1 : kind=="Surface" ? 2 : 3
        dimension==source_dimension || throw(ArgumentError(
            "$caller: $query for this Physical group requires $source_kind blocks; " *
            "got $kind"))
        raw_tags=String(strip(entity_match.captures[2]))
        tags=raw_tags==":" ? _geo_exec_all_entity_tags(
            m,dimension,"$caller $query $kind selector") :
            _geo_exec_entity_tags(
                raw_tags,context,"$caller $query $kind selector";signed=true)
        for tag in tags
            tag==typemin(Int32) && throw(ArgumentError(
                "$caller $query $kind selector entry magnitude exceeds Int32"))
        end
        length(entities)<=_MAX_GEO_LIST_ITEMS-length(tags) || throw(ArgumentError(
            "$caller: $query expands beyond $_MAX_GEO_LIST_ITEMS entities"))
        append!(entities,((dimension,tag) for tag in tags))
    end
    tags=_model_boundary(
        m,entities,caller;combined=query=="CombinedBoundary",
        max_entities=_MAX_GEO_LIST_ITEMS)
    isempty(tags) && throw(ArgumentError(
        "$caller: $query matched no entities after boundary cancellation; " *
        "select topology with a nonempty boundary"))
    return sort!(unique!(tags))
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
        entity_name=="Volume" ? 3 :
        throw(ArgumentError(
            "execute_geo: only Line/Curve, Surface, and Volume " *
            "periodicity is implemented"))
    caller="execute_geo: Periodic $(dim==1 ? "Curve" : dim==2 ? "Surface" : "Volume")"
    statement=match(
        r"^Periodic\s+(?:Line|Curve|Surface|Volume)\s*\{\s*([^}]*)\s*\}\s*=\s*\{\s*([^}]*)\s*\}\s*(.*?)\s*;$",
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
        delete_a && _remove_volume_entity!(m,a)
        delete_b && _remove_volume_entity!(m,b)
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
    elseif (physical=_geo_physical_declaration(line)) !== nothing
        caller="execute_geo: Physical $(physical.kind)"
        dim=Dict("Point"=>0,"Curve"=>1,"Line"=>1,
                 "Surface"=>2,"Volume"=>3)[physical.kind]
        name=physical.name===nothing ? "" : physical.name
        physical.tag_source===nothing && isempty(name) && throw(ArgumentError(
            "$caller: an automatic group requires a nonempty name"))
        if !isempty(name)
            for ((existing_dim,existing_tag),existing_name) in m.physical_names
                if existing_dim==dim && existing_name==name
                    throw(ArgumentError(
                        "$caller: name $(repr(name)) already identifies " *
                        "Physical $(physical.kind)[$existing_tag]"))
                end
            end
        end
        tag=physical.tag_source===nothing ? 0 : _geo_exec_entity_tag(
            physical.tag_source,context,"$caller tag")
        membership=physical.membership
        ids=if match(
                r"^(?:Boundary|CombinedBoundary|PointsOf)(?:\s|\{)",
                membership)!==nothing
            _geo_exec_physical_topology(
                m,membership,context,"$caller entities",dim)
        else
            occursin(r"[A-Za-z_][A-Za-z0-9_]*\s*\{",membership) &&
                throw(ArgumentError(
                    "$caller entities: unsupported topology query; use inline " *
                    "Boundary, CombinedBoundary, or PointsOf for Physical Point"))
            _geo_exec_entity_rhs_tags(
                membership,context,"$caller entities")
        end
        add_physical_group!(m,dim,ids;tag=tag,name=name)
        return
    elseif occursin(r"^Physical(?:\s|$)",line)
        throw(ArgumentError(
            "execute_geo: malformed Physical declaration; use " *
            "Physical Kind(tag), Physical Kind(\"name\", tag), or " *
            "Physical Kind(\"name\")"))
    elseif (mm=match(
            r"^(MeshSize|Characteristic\s+Length)\s*\{\s*(.*?)\s*\}\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: $(mm.captures[1])"
        selector=String(strip(mm.captures[2]))
        point_tags=if selector==":"
            sort!(collect(keys(m.points)))
        elseif match(r"^PointsOf(?:\s|\{)",selector)!==nothing
            _geo_exec_points_of(m,selector,context,caller)
        else
            occursin(r"[A-Za-z_][A-Za-z0-9_]*\s*\{",selector) &&
                throw(ArgumentError(
                    "$caller: unsupported topology query; use inline PointsOf " *
                    "with Point, Curve/Line, Surface, or explicit Volume blocks"))
            _geo_exec_entity_tags(selector,context,"$caller Point selector")
        end
        isempty(point_tags) && throw(ArgumentError(
            "$caller: Point selector matched no explicit modeled points"))
        mesh_size=_geo_eval_numeric(
            mm.captures[3],context,"$caller value")
        set_point_mesh_size!(m,point_tags,mesh_size)
        return
    elseif match(
            r"^SetMaxTag\s+(?:Point|Curve|Surface|Volume)\s*\(\s*.+\s*\)\s*;$",
            line) !== nothing
        return
    elseif (mm=match(
            r"^Mesh\.TransfiniteTri\s*=\s*(.*?)\s*;$",line)) !== nothing
        value=_geo_eval_numeric(mm.captures[1],context,
                                "execute_geo: Mesh.TransfiniteTri")
        isinteger(value) || throw(ArgumentError(
            "execute_geo: Mesh.TransfiniteTri must be 0 or 1 (got $value)"))
        set_transfinite_tri!(m,Int(value))
        return Int(value)
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
