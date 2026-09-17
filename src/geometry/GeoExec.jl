"""
    GeoExec

Execute a bounded subset of Gmsh `.geo`: Point/Line/Circle/Ellipse/
Spline/BSpline/Bezier/Nurbs (with `Knots`/`Order`)/Line Loop/
Plane Surface/Surface/Ruled Surface/Surface Loop/Volume,
Box/Cylinder/Sphere/Cone/Torus,
Boolean union/difference/intersection,
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
using ..Model: add_circle_arc!, add_ellipse_arc!, add_ruled_surface!
using ..Model: add_spline!, add_bspline!, add_bezier!, add_nurbs!
using ..Model: add_surface_loop!, add_volume!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, add_torus!, boolean_volumes!, boolean_volumes_multi!
using ..Model: _remove_volume_entity!
using ..Model: embed!, translate_volume!, dilate_volume!, rotate_volume!
using ..Model: transform_entities!, duplicate_entities!, coherence!
using ..Model: merge_vertices!, extrude_entities!, revolve_entities!
using ..Model: _GeoExtrudeParams
using ..Model: _affine_translation, _affine_dilation
using ..Model: _affine_rotation, _affine_symmetry, _entity_label
using ..Model: add_physical_group!, set_periodic!, set_transfinite_tri!
using ..Model: _model_boundary, _model_points_of, _model_direct_boundary
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..IO: read_geo_params, _GeoNumericContext, _geo_eval_numeric
using ..IO: _geo_split_list, _geo_split_range, _geo_range_count
using ..IO: _geo_numeric_list_terms, _geo_numeric_list_values
using ..IO: _geo_positive_gmsh_tag
using ..IO: _geo_signed_gmsh_int_value
using ..IO: _GEO_SIDE_EFFECT_SYMBOLS
using ..IO: _MAX_GEO_LIST_ITEMS
using ..IO: _geo_context_set_scalar!, _geo_context_set_list!
using ..IO: _geo_apply_list_assignment!, _geo_list_term_values
using ..IO: _geo_brace_terminated_statement, _geo_extrude_shape_group
using ..IO: _GeoTagAllocatorState, _geo_context_refresh_allocators!
using ..IO: _geo_allocator_observe_statement!,_geo_allocator_resync_model!
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
    # Final scalar/list variables from the `.geo` program (`out[]`,
    # `x[]`, ...), for inspection and differential validation.
    lists::Dict{String,Vector{Float64}}
end

const _MAX_GEO_EXEC_STATEMENT_BYTES=1_000_000
const _MAX_GEO_EXEC_STATEMENTS=1_000_000
const _MAX_GEO_LOOP_ITERATIONS=1_000_000
const _MAX_GEO_CALL_DEPTH=128

const _GEO_CONTROL_BARE=Dict(
    "Else"=>:else,"EndIf"=>:endif,"EndWhile"=>:endwhile,"EndFor"=>:endfor,
    "Return"=>:return)

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
        r"^(If|ElseIf|While|Else|EndIf|EndWhile|For|EndFor|Function|Return)\b",
        rest)
    matched===nothing && return nothing
    word=matched.captures[1]
    haskey(_GEO_CONTROL_BARE,word) && return (word,i+sizeof(word)-1)
    j=firstindex(rest)+sizeof(word)
    jlast=lastindex(rest)
    while j<=jlast && isspace(rest[j])
        j=nextind(rest,j)
    end
    if word=="Function"
        # Gmsh's `Function name` header carries an identifier or quoted
        # string-expression name on the same line and no `;`.
        name_match=match(
            r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\"|'[^']*')",SubString(rest,j))
        name_match===nothing && throw(ArgumentError(
            "execute_geo: malformed Function header; use `Function name`"))
        consumed=j-1+sizeof(name_match.match)
        return (String(rest[firstindex(rest):consumed]),i+consumed-1)
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
    elseif (matched=match(r"^Function\s+(.+)$",line))!==nothing
        return (kind=:function,name=String(matched.captures[1]))
    elseif haskey(_GEO_CONTROL_BARE,line)
        return (kind=_GEO_CONTROL_BARE[line],)
    end
    return nothing
end

# Brace-aware statement split so BooleanDifference `{ Volume{1}; Delete; }{...};`
# is one statement. Quoted strings and line/block comments are respected.
# Control constructs (If/ElseIf/Else/EndIf, For/EndFor, While/EndWhile) carry no
# `;`; each is emitted as its own statement whenever it starts at a statement
# boundary. Transform and query statements may end at their final `}` (the Gmsh
# grammar needs no `;` there) or at a `;`.
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
                if depth==0 && _geo_brace_terminated_statement(
                        String(buf.data[1:position(buf)]))
                    statement=strip(String(take!(buf)))
                    length(statements)<_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                        "execute_geo: input exceeds $_MAX_GEO_EXEC_STATEMENTS statements"))
                    push!(statements,statement)
                    buf_has_content=false
                end
            elseif c==';' && depth==0
                write(buf,c)
                statement=strip(String(take!(buf)))
                # A `;` left after a `}`-terminated transform, or a bare `;`
                # between statements, is an empty statement.
                if !all(==(';'),statement)
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
$_MAX_GEO_EXEC_STATEMENTS. `Function name ... Return` registers a zero-argument
body — the statements up to the first `Return` marker, matching Gmsh's
token-level capture — and `Call name;` re-executes it in the shared variable
scope with recursion bounded by $_MAX_GEO_CALL_DEPTH; names are identifiers or
quoted string literals, redefinition and call-before-definition fail like Gmsh.
Numeric lists support zero-based scalar indexing, `#name[]`
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
tracked explicit topology and supported full Box/Cylinder/Sphere/Cone/Torus
primitives.
`SetMaxTag Point|Curve|Surface|Volume` follows the active factory: Built-in sets the
checked counter, while OpenCASCADE only raises it. Reads use the greatest counter
among activated factories. Later primitive allocation still accounts for occupied
hidden topology. `Box` materializes Gmsh's exact boundary entities (eight
points, twelve curves, six surfaces) on the model; Cylinder/Sphere/Cone/Torus
materialize the matching OCC layouts — analytic records on every face and
edge, Plane caps, degenerate apex and seam edges, and signed shells.
`Mesh.TransfiniteTri = 0|1` sets the model's
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
Unknown entities, empty combined boundaries, unmaterialized volume topology,
unsupported geometry-derived selectors, and invalid query dimensions are
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
    # `Extrude{...}{...}` and `BooleanX{...}{...}` are side-effecting value
    # terms in the `.geo` grammar; the numeric evaluator calls back through
    # this hook.
    context.exec_hook=src->_geo_exec_value_term(model,src,context)
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
    return GeoExecution(model,mesh,params,transfinite_tri,context.lists)
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
            if match(r"^Call\b",line)!==nothing
                assigned=_geo_exec_call!(
                    m,statements,line,context,allocator_state,executed)
                assigned===nothing || (transfinite_tri=assigned)
                i+=1
                continue
            end
            occursin(
                r"\b(Macro|Fillet|Chamfer)\b",
                line) && throw(ArgumentError(
                "execute_geo: unsupported statement $(line) — macros and " *
                "advanced OCC features are blockers"))
            # `Extrude {shapes} Using Wire {n}` splits at the shape-list
            # brace; reject the OCC pipe continuation before the extrusion
            # itself executes.
            if match(r"^Extrude\b",line)!==nothing && i+1<=hi &&
               match(r"^Using\b",statements[i+1])!==nothing
                throw(ArgumentError(
                    "execute_geo: pipe extrusion `Extrude {..} Using Wire " *
                    "{..}` is OpenCASCADE-only and not implemented"))
            end
            _geo_context_refresh_allocators!(context,allocator_state)
            assigned=_exec_line!(m,line,context)
            assigned===nothing || (transfinite_tri=assigned)
            _geo_allocator_observe_statement!(
                allocator_state,line,context,"execute_geo")
            occursin(r"\bBoolean(?:Difference|Union|Intersection|Fragments)?\b",
                     line) &&
                _geo_allocator_resync_model!(allocator_state,m)
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
        elseif control.kind===:function
            i=_geo_exec_function_def!(statements,i,hi,context)
        elseif control.kind===:return
            throw(ArgumentError(
                "execute_geo: Return without an enclosing Function"))
        else
            throw(ArgumentError(
                "execute_geo: $line without a matching opener"))
        end
    end
    return transfinite_tri
end

# A `.geo` Function/Call name is a bare identifier or a quoted string literal
# in the bounded subset (Gmsh accepts general string expressions; Tessella
# has no string variables, so anything else is an explicit blocker).
function _geo_function_name(raw::AbstractString)
    name=strip(raw)
    match(r"^[A-Za-z_][A-Za-z0-9_]*$",name)!==nothing && return String(name)
    quoted=match(r"^\"(.*)\"$",name)
    quoted===nothing && (quoted=match(r"^'(.*)'$",name))
    quoted!==nothing && return String(quoted.captures[1])
    throw(ArgumentError(
        "execute_geo: Function/Call name must be an identifier or a " *
        "quoted string literal"))
end

# `Function name` registers its body — the statements up to the first `Return`
# marker, matching Gmsh's token-level capture — without executing it. The body
# may itself contain a `Function` statement, which registers when the outer
# body runs (Gmsh's `Call` semantics resolve names at call time).
function _geo_exec_function_def!(statements::Vector{String},i::Int,hi::Int,
                                 context::_GeoNumericContext)
    control=_geo_control_parse(statements[i])
    name=_geo_function_name(control.name)
    haskey(context.functions,name) && throw(ArgumentError(
        "execute_geo: Redefinition of function $name"))
    j=i+1
    while j<=hi
        inner=_geo_control_parse(statements[j])
        (inner!==nothing && inner.kind===:return) && break
        j+=1
    end
    j<=hi || throw(ArgumentError(
        "execute_geo: Function $name has no matching Return"))
    context.functions[name]=(i+1):(j-1)
    return j+1
end

# `Call name;` re-executes the registered body in the caller's scope with a
# bounded recursion depth (Gmsh recurses by re-parsing the body text; an
# unbounded call stack is a resource-bound violation here).
function _geo_exec_call!(m::GeoModel,statements::Vector{String},
                         line::AbstractString,
                         context::_GeoNumericContext,
                         allocator_state::_GeoTagAllocatorState,
                         executed::Base.RefValue{Int})
    matched=match(r"^Call\s+(.+?)\s*;?\s*$",line)
    matched===nothing && throw(ArgumentError(
        "execute_geo: malformed Call statement; use `Call name;`"))
    name=_geo_function_name(matched.captures[1])
    body=get(context.functions,name,nothing)
    body===nothing && throw(ArgumentError(
        "execute_geo: Unknown function '$name'"))
    context.call_depth+=1
    try
        context.call_depth<=_MAX_GEO_CALL_DEPTH || throw(ArgumentError(
            "execute_geo: Call depth exceeds $_MAX_GEO_CALL_DEPTH"))
        return _exec_geo_statements!(
            m,statements,first(body),last(body),
            context,allocator_state,executed)
    finally
        context.call_depth-=1
    end
end

function _boolean_delete_operand(raw::AbstractString)
    suffix=String(strip(raw))
    match(r"^;?\s*(?:Delete\s*;?)?$",suffix)===nothing && throw(ArgumentError(
        "execute_geo: Boolean operand suffix must contain only optional Delete; got $(repr(suffix))"))
    return occursin(r"\bDelete\b",suffix)
end

const _GEO_BOOLEAN_OPS=Dict("Difference"=>:difference,"Union"=>:union,
                            "Intersection"=>:intersection,
                            "Fragments"=>:fragments)

# Consecutive `{ ... }` operand groups of a `BooleanX ...` term: each group's
# content is returned without its braces. The statement splitter already
# keeps the whole form together; anything that is not a group is malformed.
function _geo_boolean_groups(raw::AbstractString,caller::AbstractString)
    s=String(strip(raw))
    endswith(s,";") &&
        (s=String(strip(s[firstindex(s):prevind(s,lastindex(s))])))
    groups=String[]
    while !isempty(s)
        startswith(s,"{") || throw(ArgumentError(
            "$caller: expected a `{ ... }` operand group; got $(repr(s))"))
        (content,s)=_geo_balanced_group(s,caller)
        push!(groups,content)
    end
    return groups
end

# One operand group `{ Volume{tags}; Delete; }` — `;`-separated clauses of
# `Kind{...}` entity lists and `Delete` markers. Only Volume operands are
# supported; other dimension clauses fail explicitly.
function _geo_boolean_operand_group(raw::AbstractString,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    tags=Int[];del=false
    for clause in split(raw,';')
        text=strip(String(clause))
        isempty(text) && continue
        text=="Delete" && (del=true;continue)
        mm=match(r"^([A-Za-z]+)\s*\{(.*)\}$",text)
        mm===nothing && throw(ArgumentError(
            "$caller: malformed operand clause $(repr(text))"))
        mm.captures[1]=="Volume" || throw(ArgumentError(
            "$caller: only Volume Boolean operands are supported; " *
            "got $(mm.captures[1]){$(mm.captures[2])}"))
        append!(tags,_geo_exec_entity_tags(
            mm.captures[2],context,"$caller Volume operands"))
    end
    return (tags=tags,delete=del)
end

function _geo_boolean_operands(groups::Vector{String},
                               context::_GeoNumericContext,
                               caller::AbstractString)
    length(groups)==2 || throw(ArgumentError(
        "$caller: expected `{ objects }{ tools }` operand groups; " *
        "got $(length(groups)) group$(length(groups)==1 ? "" : "s")"))
    objects=_geo_boolean_operand_group(groups[1],context,caller)
    tools=_geo_boolean_operand_group(groups[2],context,caller)
    isempty(objects.tags) && throw(ArgumentError(
        "$caller: the object operand group is empty"))
    return (objects=objects,tools=tools)
end

# A `BooleanX{...}{...}` value term — the form allowed inside `v() =` /
# `name[] =` captures and entity lists. Returns nothing when `raw` is not a
# Boolean term (the `(t) =` tagged form is a statement, not a term).
function _geo_exec_boolean_term(m::GeoModel,raw::AbstractString,
                                context::_GeoNumericContext)
    source=String(strip(raw))
    mm=match(r"^Boolean(Difference|Union|Intersection|Fragments)\b(.*)$",source)
    mm===nothing && return nothing
    caller="execute_geo: Boolean$(mm.captures[1])"
    rest=String(strip(mm.captures[2]))
    startswith(rest,"(") && return nothing
    groups=_geo_boolean_groups(rest,caller)
    operands=_geo_boolean_operands(groups,context,caller)
    out=boolean_volumes_multi!(m,_GEO_BOOLEAN_OPS[mm.captures[1]],
        operands.objects.tags,operands.tools.tags;
        remove_object=operands.objects.delete,
        remove_tool=operands.tools.delete,caller=caller)
    return Float64.(out)
end

function _geo_exec_value_term(m::GeoModel,raw::AbstractString,
                              context::_GeoNumericContext)
    extruded=_geo_exec_extrude_term(m,raw,context)
    extruded!==nothing && return extruded
    return _geo_exec_boolean_term(m,raw,context)
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
        append!(values,_geo_list_term_values(term))
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

# Extract one balanced `{...}` group from the front of `raw`; returns
# `(content, rest)`.
function _geo_balanced_group(raw::AbstractString, caller::AbstractString)
    s=String(strip(raw))
    startswith(s,"{") || throw(ArgumentError(
        "$caller: expected a `{...}` group; got $(repr(s))"))
    depth=0;closing=0;i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if c=='{'
            depth+=1
        elseif c=='}'
            depth-=1
            depth==0 && (closing=i; break)
            depth<0 && throw(ArgumentError(
                "$caller: unmatched closing brace"))
        end
        i=nextind(s,i)
    end
    closing==0 && throw(ArgumentError("$caller: unmatched opening brace"))
    content=closing>2 ? String(s[2:prevind(s,closing)]) : ""
    rest=closing<last ? String(strip(s[nextind(s,closing):end])) : ""
    return (content,rest)
end

_geo_word_char(c::AbstractChar) = isletter(c) || isdigit(c) || c=='_'

# Split `raw` at the first top-level occurrence of `keyword` — a word-boundary
# identifier outside every `{}`, `()`, and `[]` pair — returning the stripped
# text before and after it, or `nothing` when the keyword never appears at top
# level (an empty `after` then still distinguishes a trailing keyword).
function _geo_split_at_keyword(raw::AbstractString, keyword::AbstractString,
                               caller::AbstractString)
    s=String(raw);depth=0;parens=0;brackets=0
    i=firstindex(s);last=lastindex(s);klen=ncodeunits(keyword)
    while i<=last
        c=s[i]
        if c=='{'
            depth+=1
        elseif c=='}'
            depth-=1
            depth<0 && throw(ArgumentError(
                "$caller: unmatched closing brace"))
        elseif c=='('
            parens+=1
        elseif c==')'
            parens-=1
            parens<0 && throw(ArgumentError(
                "$caller: unmatched closing parenthesis"))
        elseif c=='['
            brackets+=1
        elseif c==']'
            brackets-=1
            brackets<0 && throw(ArgumentError(
                "$caller: unmatched closing bracket"))
        elseif depth==0 && parens==0 && brackets==0
            j=i+klen
            if j-1<=last && isvalid(s,j-1) && s[i:j-1]==keyword &&
               (i==firstindex(s) || !_geo_word_char(s[prevind(s,i)])) &&
               (j>last || (isvalid(s,j) && !_geo_word_char(s[j])))
                return (String(strip(s[firstindex(s):prevind(s,i)])),
                        String(strip(s[j:last])))
            end
        end
        i=nextind(s,i)
    end
    return nothing
end

# Split `raw` on commas outside every `{}`, `()`, and `[]` pair.
function _geo_split_top_commas(raw::AbstractString, caller::AbstractString)
    s=String(raw);parts=String[];depth=0;parens=0;brackets=0
    start=firstindex(s);i=start;last=lastindex(s)
    while i<=last
        c=s[i]
        if c=='{';depth+=1
        elseif c=='}';depth-=1;depth>=0 || throw(ArgumentError(
            "$caller: unmatched closing brace"))
        elseif c=='(';parens+=1
        elseif c==')';parens-=1;parens>=0 || throw(ArgumentError(
            "$caller: unmatched closing parenthesis"))
        elseif c=='[';brackets+=1
        elseif c==']';brackets-=1;brackets>=0 || throw(ArgumentError(
            "$caller: unmatched closing bracket"))
        elseif c==',' && depth==0 && parens==0 && brackets==0
            push!(parts,String(strip(s[start:prevind(s,i)])))
            start=nextind(s,i)
        end
        i=nextind(s,i)
    end
    (depth==0 && parens==0 && brackets==0) || throw(ArgumentError(
        "$caller: unmatched opening delimiter"))
    push!(parts,String(strip(s[start:last])))
    return parts
end

_geo_shape_kind_dim(name::AbstractString) =
    name=="Point" ? 0 : name in ("Curve","Line") ? 1 :
    name=="Surface" ? 2 : 3

# Entity-block tag entries: a `:` wildcard or a tag list (signed when the
# enclosing construct allows it — `Extrude` curves carry orientation).
function _geo_shape_tags(m::GeoModel, dim::Int, raw::AbstractString,
                         context::_GeoNumericContext, caller::AbstractString;
                         signed_tags::Bool=false)
    r=String(strip(raw))
    r==":" && return _geo_exec_all_entity_tags(m,dim,"$caller wildcard")
    return _geo_exec_entity_tags(r,context,caller;signed=signed_tags)
end

function _geo_shape_physical_tags(m::GeoModel, dim::Int, raw::AbstractString,
                                  context::_GeoNumericContext,
                                  caller::AbstractString)
    r=String(strip(raw))
    r==":" && return sort!([gt for (d,gt) in keys(m.physical) if d==dim])
    return _geo_exec_entity_tags(r,context,caller)
end

function _geo_require_list_semicolon(s::AbstractString, caller, name)
    rest=String(strip(s))
    startswith(rest,";") || throw(ArgumentError(
        "$caller: $name{...} entries in a transform list must end with `;`"))
    return String(strip(rest[nextind(rest,1):end]))
end

# Evaluate a `.geo` MultipleShape, mirroring the Gmsh grammar
# (`MultipleShape : ListOfShapes | Transform`). A ListOfShapes is a sequence
# of `;`-terminated entries: entity blocks (`Point{...}`, `Curve{...}`,
# `Surface{...}`, `Volume{...}`), `Kind{:}` wildcards, `Physical Kind{...}`
# and `Parent Kind{...}` selectors, and inline Shape definitions
# (`Point(7) = {..};`). A Transform — a named transform or an action string
# like `Duplicata`/`Boundary`/`PointsOf` applied to a nested MultipleShape —
# is valid only as the whole list, so nested transforms and actions cannot be
# mixed with other entries. Copies and nested transforms execute as side
# effects; the resulting entity list is returned for the enclosing transform.
function _geo_shape_list_entities!(m::GeoModel, raw::AbstractString,
                                   context::_GeoNumericContext,
                                   caller::AbstractString)
    entities=NTuple{2,Int}[]
    s=String(strip(raw))
    while !isempty(s)
        (entries,rest,transform)=
            _geo_multiple_shape_element!(m,s,context,caller)
        if transform && (!isempty(entities) || !isempty(rest))
            throw(ArgumentError(
                "$caller: a transform or shape action cannot be combined " *
                "with other transform-list entries"))
        end
        append!(entities,entries)
        s=rest
        length(entities)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
            "$caller: shape list expands beyond $_MAX_GEO_LIST_ITEMS entities"))
    end
    return entities
end

# One MultipleShape element. Returns (entities, remaining_text, is_transform).
# `signed_tags` allows negative entity tags (`Curve{-1}`) — Gmsh's
# `RecursiveListOfDouble` permits them and `Extrude` reads the sign as the
# reversed-record orientation.
function _geo_multiple_shape_element!(m::GeoModel, s0::AbstractString,
                                      context::_GeoNumericContext,
                                      caller::AbstractString;
                                      signed_tags::Bool=false)
    mm=match(r"^([A-Za-z_][A-Za-z0-9_]*)",s0)
    mm===nothing && throw(ArgumentError(
        "$caller: malformed shape list element near $(repr(s0))"))
    name=mm.captures[1]
    s=String(strip(s0[nextind(s0,firstindex(s0),ncodeunits(mm.match)):end]))
    if name in ("Translate","Rotate","Dilate","Symmetry","Affine","Closest")
        (params,s)=_geo_balanced_group(s,caller)
        (inner_list,s)=_geo_balanced_group(s,caller)
        (name=="Affine" || name=="Closest") && throw(ArgumentError(
            "$caller: $name transforms require the OpenCASCADE geometry " *
            "kernel, which Tessella does not implement"))
        nested="$caller: $name"
        t=_geo_transform_params(kind=name,params=params,
                                context=context,caller=nested)
        inner=_geo_shape_list_entities!(m,inner_list,context,caller)
        transform_entities!(m,t,inner;caller=nested)
        # Gmsh returns the input shape list (`$$ = $MultipleShape`).
        return (inner,s,true)
    end
    if name=="Physical" || name=="Parent"
        km=match(r"^(Point|Curve|Line|Surface|Volume)\b",s)
        if km===nothing && (em=match(r"^GeoEntity\b",s))!==nothing
            s=String(strip(s[nextind(s,firstindex(s),ncodeunits(em.match)):end]))
            (dgroup,s)=_geo_balanced_group(s,caller)
            dim_raw=_geo_eval_numeric(dgroup,context,
                                      "$caller $name GeoEntity dimension")
            isinteger(dim_raw) || throw(ArgumentError(
                "$caller: $name GeoEntity dimension must be an integer"))
            dim=Int(dim_raw)
        else
            km===nothing && throw(ArgumentError(
                "$caller: $name requires Point, Curve/Line, Surface, " *
                "Volume, or GeoEntity{d} inside a shape list"))
            dim=_geo_shape_kind_dim(km.captures[1])
            s=String(strip(s[nextind(s,firstindex(s),ncodeunits(km.match)):end]))
        end
        (group,s)=_geo_balanced_group(s,caller)
        s=_geo_require_list_semicolon(s,caller,name)
        # Tessella entities have no parent entities; the selector is empty.
        name=="Parent" && return (NTuple{2,Int}[],s,false)
        entries=NTuple{2,Int}[]
        for gtag in _geo_shape_physical_tags(m,dim,group,context,caller)
            haskey(m.physical,(dim,gtag)) || throw(ArgumentError(
                "$caller: unknown Physical $(_entity_label(dim))[$gtag]"))
            for tag in m.physical[(dim,gtag)]
                push!(entries,(dim,tag))
            end
        end
        return (entries,s,false)
    end
    if name=="GeoEntity"
        # `GeoEntity{dim}{tags};` — the generic `tGeoEntity` selector.
        (dgroup,s)=_geo_balanced_group(s,caller)
        dim_raw=_geo_eval_numeric(dgroup,context,"$caller GeoEntity dimension")
        isinteger(dim_raw) || throw(ArgumentError(
            "$caller: GeoEntity dimension must be an integer; got $dim_raw"))
        dim=Int(dim_raw)
        0<=dim<=3 || throw(ArgumentError(
            "$caller: GeoEntity dim out of range [0,3]"))
        (group,s)=_geo_balanced_group(s,caller)
        s=_geo_require_list_semicolon(s,caller,name)
        return (NTuple{2,Int}[(dim,t) for t in _geo_shape_tags(
                    m,dim,group,context,caller;signed_tags=signed_tags)],
                s,false)
    end
    if name in ("Point","Curve","Line","Surface","Volume")
        dim=_geo_shape_kind_dim(name)
        if startswith(s,"{")
            (group,s)=_geo_balanced_group(s,caller)
            s=_geo_require_list_semicolon(s,caller,name)
            return (NTuple{2,Int}[(dim,t) for t in _geo_shape_tags(
                        m,dim,group,context,caller;signed_tags=signed_tags)],
                    s,false)
        end
        # `Kind(tag) = rhs;` — an inline Shape definition.
        return _geo_shape_definition_element!(m,s0,context,caller)
    end
    if startswith(s,"{")
        # tSTRING '{' MultipleShape '}' — Duplicata, a boundary query, or an
        # unknown action.
        (group,s)=_geo_balanced_group(s,caller)
        inner=_geo_shape_list_entities!(m,group,context,caller)
        return (_geo_shape_action!(m,name,inner,caller),s,true)
    end
    (name=="Split" || name=="Intersect") && throw(ArgumentError(
        "$caller: $name requires built-in kernel curve splitting, which " *
        "Tessella does not implement"))
    (startswith(s,"(") || match(r"^[A-Za-z_]",s)!==nothing) &&
        return _geo_shape_definition_element!(m,s0,context,caller)
    throw(ArgumentError("$caller: unsupported shape list element '$name'"))
end

# A `Name(tag) = rhs;` definition (possibly multi-word, e.g. `Plane Surface`)
# inside a shape list executes through the normal statement executor and adds
# the created entity to the list, as in the Gmsh `Shape` production. Loop
# records and non-entity definitions are not transformable.
const _GEO_SHAPE_DEFINITION_DIMS=Dict(
    "Point"=>0,
    "Line"=>1,"Curve"=>1,"Spline"=>1,"BSpline"=>1,"Bezier"=>1,"Nurbs"=>1,
    "Circle"=>1,"Ellipse"=>1,"Wire"=>1,"Compound Spline"=>1,
    "Compound BSpline"=>1,"Compound Curve"=>1,
    "Surface"=>2,"Plane Surface"=>2,"Ruled Surface"=>2,"BSpline Surface"=>2,
    "Bezier Surface"=>2,"Parametric Surface"=>2,"Rectangle"=>2,"Disk"=>2,
    "Compound Surface"=>2,
    "Volume"=>3,"Sphere"=>3,"PolarSphere"=>3,"Box"=>3,"Torus"=>3,
    "Cylinder"=>3,"Cone"=>3,"Wedge"=>3,"ThickSolid"=>3,"Compound Volume"=>3)

function _geo_shape_definition_element!(m::GeoModel, s0::AbstractString,
                                        context::_GeoNumericContext,
                                        caller::AbstractString)
    depth=0;parens=0;i=firstindex(s0);last=lastindex(s0);cut=0
    while i<=last
        c=s0[i]
        c=='{' && (depth+=1)
        c=='}' && (depth-=1)
        c=='(' && (parens+=1)
        c==')' && (parens-=1)
        (c==';' && depth==0 && parens==0) && (cut=i; break)
        i=nextind(s0,i)
    end
    cut==0 && throw(ArgumentError(
        "$caller: shape list definition entries must end with `;`"))
    stmt=String(s0[firstindex(s0):cut])
    rest=String(strip(s0[nextind(s0,cut):end]))
    head=match(r"^((?:[A-Za-z_][A-Za-z0-9_]*\s+)*[A-Za-z_][A-Za-z0-9_]*)\s*\(",
               stmt)
    head===nothing && throw(ArgumentError(
        "$caller: malformed shape list element near $(repr(stmt))"))
    kind=String(strip(replace(head.captures[1],r"\s+"=>" ")))
    dim=get(_GEO_SHAPE_DEFINITION_DIMS,kind,-1)
    dim<0 && throw(ArgumentError(
        "$caller: '$kind' entries in a transform list do not produce " *
        "transformable entities"))
    tm=match(r"\(([^()]*)\)",stmt)
    tag=_geo_exec_entity_tag(tm.captures[1],context,"$caller $kind tag")
    _exec_line!(m,stmt,context)
    return (NTuple{2,Int}[(dim,tag)],rest,false)
end

# A `Name{ MultipleShape }` action element: Duplicata copies, boundary
# queries, and PointsOf. Returns the action's output entities.
function _geo_shape_action!(m::GeoModel, name::AbstractString,
                            inner::Vector{NTuple{2,Int}},
                            caller::AbstractString)
    if name=="Duplicata"
        return duplicate_entities!(m,inner;caller="$caller: Duplicata")
    elseif name in ("Boundary","CombinedBoundary","OrientedBoundary",
                    "OrientedCombinedBoundary")
        combined=occursin("Combined",name)
        oriented=occursin("Oriented",name)
        entries=NTuple{2,Int}[]
        for d in sort!(unique!(first.(inner));rev=true)
            d==0 && continue   # Points have no boundary (Gmsh returns none)
            group=[e for e in inner if e[1]==d]
            if oriented
                for (_,signed_tag) in Iterators.flatten(
                        _model_direct_boundary(
                            m,d,t,caller;canonical_orientation=false)
                        for (_,t) in group)
                    push!(entries,(d-1,abs(signed_tag)))
                end
            else
                for tag in _model_boundary(
                        m,group,caller;
                        combined=combined,max_entities=_MAX_GEO_LIST_ITEMS)
                    push!(entries,(d-1,tag))
                end
            end
        end
        return entries
    elseif name=="PointsOf"
        return NTuple{2,Int}[(0,t) for t in _model_points_of(
            m,inner,"$caller: PointsOf")]
    end
    throw(ArgumentError(
        "$caller: unknown action on multiple shapes '$name'"))
end

# ── Extrude ──────────────────────────────────────────────────────────────
#
# `.geo` translational `Extrude {dx,dy,dz} { ...; }`, mirroring the Gmsh
# grammar (`tExtrude VExpr '{' ListOfShapes ExtrudeParameters '}'`) and the
# `GEO_Internals::extrude` → `ExtrudeShapes` semantics: the expression value
# is the flat `[top, body, laterals...]` tag list (laterals only under
# `Geometry.ExtrudeReturnLateralEntities`, which defaults on). Rotational
# (`{{axis}, {point}, angle}`), twist (`{{axis}, {point}, {delta}, angle}`),
# boundary-layer (`Extrude {shapes; params}`), and `Using Wire` pipe forms are
# recognized and rejected with explicit errors.

# Scan a balanced `(...)` group — the VExpr grammar also allows parenthesized
# vectors — returning `(content, rest)`.
function _geo_balanced_paren(raw::AbstractString, caller::AbstractString)
    s=String(strip(raw))
    startswith(s,"(") || throw(ArgumentError(
        "$caller: expected a `(...)` group; got $(repr(s))"))
    depth=0;closing=0;i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if c=='('
            depth+=1
        elseif c==')'
            depth-=1
            depth==0 && (closing=i; break)
            depth<0 && throw(ArgumentError(
                "$caller: unmatched closing parenthesis"))
        end
        i=nextind(s,i)
    end
    closing==0 && throw(ArgumentError("$caller: unmatched opening parenthesis"))
    content=closing>2 ? String(s[2:prevind(s,closing)]) : ""
    rest=closing<last ? String(strip(s[nextind(s,closing):end])) : ""
    return (content,rest)
end

# One `VExpr` member of a revolve/twist motion group evaluated to a flat
# numeric vector; a braced member loses its braces, a bare term evaluates as
# a list term (list variables expand).
function _geo_extrude_vector_values(part::AbstractString,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    s=String(strip(part))
    if startswith(s,"{") || startswith(s,"(")
        content,rest=s[1]=='{' ? _geo_balanced_group(s,caller) :
                                 _geo_balanced_paren(s,caller)
        isempty(strip(rest)) || throw(ArgumentError(
            "$caller: malformed displacement expression $(repr(s))"))
        return _geo_exec_numeric_values(content,context,caller)
    end
    return _geo_exec_numeric_values(s,context,caller)
end

# The revolve/twist `{{A}, {X}, alpha}` / `{{A}, {X}, {T}, alpha}` motion
# group, decoded from the group's evaluated element lengths: the grammar's
# VExpr members are length-3 vectors and the trailing FExpr a scalar.
function _geo_extrude_revolve(values::Vector{Vector{Float64}},
                              caller::AbstractString)
    lengths=length.(values)
    if length(lengths)==4 && lengths[1:3]==[3,3,3] && lengths[4]==1
        throw(ArgumentError(
            "$caller: twist extrusion `Extrude {{axis}, {point}, {delta}, " *
            "angle}` is not implemented"))
    end
    if !(length(lengths)==3 && lengths==[3,3,1])
        throw(ArgumentError(
            "$caller: rotational extrusion requires `{{axis}, {point}, " *
            "angle}`"))
    end
    a,x,α=values
    return (axis=(a[1],a[2],a[3]),origin=(x[1],x[2],x[3]),angle=α[1])
end

# Accumulate one VExpr group into `delta` (`sign` folds `VExpr '+' VExpr` and
# leading `tMINUS VExpr`), or decode the revolve `{{axis}, {point}, angle}` /
# twist `{{axis}, {point}, {delta}, angle}` motion group. The group's
# evaluated element lengths decide: `[1,1,1]` (or a single list variable
# evaluating to three numbers) is the translational displacement, `[3,3,1]`
# is revolve and `[3,3,3,1]` is twist — the same shapes the grammar's
# `VExpr, VExpr, FExpr` / `VExpr, VExpr, VExpr, FExpr` alternatives produce.
# Returns `(delta, revolve)`: `revolve` is `nothing` for translation.
function _geo_extrude_vector!(delta::NTuple{3,Float64}, sign::Float64,
                              body::AbstractString,
                              context::_GeoNumericContext,
                              caller::AbstractString)
    parts=_geo_split_top_commas(body,caller)
    if length(parts)<=4
        values=[_geo_extrude_vector_values(p,context,caller) for p in parts]
        lengths=length.(values)
        revolve_shape=length(lengths)==3 ? lengths==[3,3,1] :
            length(lengths)==4 && lengths==[3,3,3,1]
        if revolve_shape || (length(parts) in (3,4) &&
                             any(v->length(v)!=1,values))
            (sign==1.0 && delta==(0.0,0.0,0.0)) || throw(ArgumentError(
                "$caller: malformed displacement expression $(repr(body))"))
            return (delta,_geo_extrude_revolve(values,caller))
        elseif length(lengths)==3 && all(v->length(v)==1,values)
            return (delta .+ sign .* ntuple(i->values[i][1],3),nothing)
        elseif length(lengths)==1 && length(values[1])==3
            return (delta .+ sign .* Tuple(values[1]),nothing)
        end
    end
    length(parts)==3 || throw(ArgumentError(
        "$caller: displacement requires exactly three components; got " *
        "$(length(parts)) in $(repr(body))"))
    return (delta .+ sign .* ntuple(
        i->_geo_eval_numeric(parts[i],context,"$caller displacement"),3),
        nothing)
end

const _GEO_EXTRUDE_PARAMS=(layers=Int[],heights=Float64[],scale_last=false,
                           recombine=false,quad_to_tri=:none,
                           recomb_laterals=false)

# Parse one `;`-stripped element of an `Extrude` shape list. Returns
# `entities` for shape elements (a `Vector{NTuple{2,Int}}`) or `nothing` for
# `ExtrudeParameter` elements, which update `params` in place.
function _geo_extrude_element!(m::GeoModel, element::AbstractString,
                               context::_GeoNumericContext,
                               params::_GeoExtrudeParams,
                               caller::AbstractString)
    mm=match(r"^([A-Za-z_][A-Za-z0-9_]*)",element)
    mm===nothing && throw(ArgumentError(
        "$caller: malformed shape list element near $(repr(element))"))
    name=mm.captures[1]
    s=String(strip(element[nextind(element,firstindex(element),
                             ncodeunits(mm.match)):end]))
    if name=="Layers"
        (group,s)=_geo_balanced_group(s,caller)
        isempty(s) || throw(ArgumentError(
            "$caller: unexpected text after Layers parameter"))
        # `Layers{n}` | `Layers{counts, heights}` | `Layers{{c..},{h..}}`:
        # the grammar is `tLayers '{' (FExpr | ListOfDouble ',' ListOfDouble)
        # '}'`, so the group's top-level items split after the first item —
        # `Layers{3,5}` is one layer of 3 elements at height 5, while
        # `Layers{{2,4},{0.2,0.8}}` is the documented two-layer form.
        items=_geo_split_top_commas(group,caller)
        expand(item)=_geo_exec_numeric_values(
            startswith(strip(item),"{") ?
            String(strip(item)[2:prevind(strip(item),end)]) : item,
            context,"$caller Layers")
        if length(items)==1
            item=String(strip(only(items)))
            startswith(item,"{") && throw(ArgumentError(
                "$caller: Layers requires `{counts, heights}` or a single " *
                "element count"))
            # `Layers{n}`: a single layer of unit height; `n==0` leaves the
            # parameters untouched (Gmsh accepts it to make disabling easy).
            n=abs(_geo_signed_gmsh_int_value(
                _geo_eval_numeric(item,context,"$caller Layers"),
                "$caller Layers"))
            n==0 && return (nothing,params)
            layers=Int[n];heights=Float64[1.0]
        else
            counts=expand(items[1])
            heights=Float64[]
            for item in items[2:end]
                append!(heights,expand(item))
            end
            length(heights)==length(counts) || throw(ArgumentError(
                "$caller: wrong layer definition {$(length(counts)), " *
                "$(length(heights))}"))
            layers=Int[]
            for value in counts
                layer=_geo_signed_gmsh_int_value(value,"$caller Layers entry")
                push!(layers,layer>0 ? layer : 1)
            end
        end
        return (nothing,(layers=layers,heights=heights,
            scale_last=params.scale_last,recombine=params.recombine,
            quad_to_tri=params.quad_to_tri,
            recomb_laterals=params.recomb_laterals))
    elseif name=="ScaleLast"
        isempty(s) || throw(ArgumentError(
            "$caller: ScaleLast takes no value"))
        return (nothing,(layers=params.layers,heights=params.heights,
            scale_last=true,recombine=params.recombine,
            quad_to_tri=params.quad_to_tri,
            recomb_laterals=params.recomb_laterals))
    elseif name=="Recombine"
        value=isempty(s) ? true :
            _geo_eval_numeric(s,context,"$caller Recombine")!=0
        return (nothing,(layers=params.layers,heights=params.heights,
            scale_last=params.scale_last,recombine=value,
            quad_to_tri=params.quad_to_tri,
            recomb_laterals=params.recomb_laterals))
    elseif name in ("QuadTriAddVerts","QuadTriNoNewVerts")
        recomb=startswith(s,"RecombLaterals")
        (recomb || isempty(s)) || throw(ArgumentError(
            "$caller: unexpected text after $name parameter"))
        recomb && !isempty(strip(s[nextind(s,firstindex(s),14):end])) &&
            throw(ArgumentError(
                "$caller: unexpected text after $name RecombLaterals"))
        kind=name=="QuadTriAddVerts" ? :add_verts : :no_new_verts
        return (nothing,(layers=params.layers,heights=params.heights,
            scale_last=params.scale_last,recombine=params.recombine,
            quad_to_tri=kind,recomb_laterals=recomb))
    elseif name=="Using"
        # `Using name[i]` is a legal ExtrudeParameter; Gmsh only acts on
        # `Index`/`View` (boundary-layer mesh metadata) and silently drops
        # every other name — inert either way for translational extrusion.
        um=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\[",s)
        um===nothing && throw(ArgumentError(
            "$caller: expected `Using name[i]`"))
        rest=String(strip(s[nextind(s,firstindex(s),ncodeunits(um.match)):end]))
        endswith(rest,"]") || throw(ArgumentError(
            "$caller: malformed Using $(um.captures[1]) parameter"))
        _geo_eval_numeric(rest[firstindex(rest):prevind(rest,end)],
                          context,"$caller Using $(um.captures[1])")
        return (nothing,params)
    elseif name=="Hole"
        throw(ArgumentError(
            "$caller: Hole extrusion parameters apply to boundary-layer " *
            "extrusion, which is not implemented"))
    elseif name=="Extrude"
        throw(ArgumentError(
            "$caller: Extrude terms cannot be nested inside an Extrude " *
            "shape list"))
    elseif name in ("Translate","Rotate","Dilate","Symmetry","Affine",
                    "Duplicata","Boundary","CombinedBoundary",
                    "OrientedBoundary","OrientedCombinedBoundary",
                    "PointsOf","Split","Intersect","Closest")
        throw(ArgumentError(
            "$caller: $name cannot appear inside an Extrude shape list"))
    end
    # A `Shape` element: `Kind{...}`, `Entity{d}{...}`, `Physical`/`Parent`
    # selectors, or an inline `Kind(tag) = rhs` definition. Feed the element
    # its stripped semicolon so the shared parsers see the full form.
    entries,rest,transform=_geo_multiple_shape_element!(
        m,element*";",context,caller;signed_tags=true)
    transform && throw(ArgumentError(
        "$caller: transforms cannot appear inside an Extrude shape list"))
    isempty(strip(rest)) || throw(ArgumentError(
        "$caller: unexpected text in Extrude shape list near $(repr(rest))"))
    return (entries,params)
end

function _geo_extrude_shape_list!(m::GeoModel, body::AbstractString,
                                  context::_GeoNumericContext,
                                  caller::AbstractString)
    params=_GEO_EXTRUDE_PARAMS
    entities=NTuple{2,Int}[]
    isempty(strip(body)) && return (entities,params)
    for element in _geo_exec_topology_query_blocks(
            body,"Extrude shape list",caller)
        entries,params=_geo_extrude_element!(
            m,element,context,params,caller)
        entries===nothing || append!(entities,entries)
        length(entities)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
            "$caller: shape list expands beyond $_MAX_GEO_LIST_ITEMS entities"))
    end
    return (entities,params)
end

# The `Extrude` value term: `nothing` when `raw` is not an extrusion (the
# numeric-list evaluator then treats it as an ordinary term), else the flat
# `[top, body, laterals...]` list as Float64. Statement-level `Extrude`
# statements route here through `_exec_line!` and discard the value.
function _geo_exec_extrude_term(m::GeoModel, raw::AbstractString,
                                context::_GeoNumericContext)
    source=String(strip(raw))
    match(r"^Extrude\b",source)===nothing && return nothing
    caller="execute_geo: Extrude"
    endswith(source,";") &&
        (source=String(strip(source[firstindex(source):prevind(source,end)])))
    rest=String(strip(source[nextind(source,firstindex(source),7):end]))
    delta=(0.0,0.0,0.0);sign=1.0;shapes=nothing;nvec=0;revolve=nothing
    while !isempty(rest)
        c=rest[firstindex(rest)]
        if c=='+' || c=='-'
            c=='-' && (sign=-sign)
            rest=String(strip(rest[nextind(rest,firstindex(rest)):end]))
            continue
        elseif c=='{'
            (group,rest)=_geo_balanced_group(rest,caller)
            if _geo_extrude_shape_group(group)
                shapes=group
                break
            end
            delta,motion=_geo_extrude_vector!(delta,sign,group,context,caller)
            if motion!==nothing
                revolve===nothing || throw(ArgumentError(
                    "$caller: multiple extrusion motion groups"))
                revolve=motion
            end
            sign=1.0;nvec+=1
        elseif c=='('
            (group,rest)=_geo_balanced_paren(rest,caller)
            delta,motion=_geo_extrude_vector!(delta,sign,group,context,caller)
            if motion!==nothing
                revolve===nothing || throw(ArgumentError(
                    "$caller: multiple extrusion motion groups"))
                revolve=motion
            end
            sign=1.0;nvec+=1
        else
            break
        end
    end
    match(r"^Using\b",rest)!==nothing && throw(ArgumentError(
        "$caller: pipe extrusion `Extrude {..} Using Wire {..}` is " *
        "OpenCASCADE-only and not implemented"))
    isempty(rest) || throw(ArgumentError(
        "$caller: unexpected text after the shape list"))
    shapes===nothing && throw(ArgumentError(
        "$caller: expected a `{shape list}` group"))
    nvec==0 && throw(ArgumentError(
        "$caller: boundary-layer extrusion `Extrude {shapes; params}` is " *
        "not implemented"))
    (revolve===nothing || delta==(0.0,0.0,0.0)) || throw(ArgumentError(
        "$caller: malformed displacement expression"))
    entities,params=_geo_extrude_shape_list!(m,shapes,context,caller)
    tags=revolve===nothing ?
        extrude_entities!(m,entities,delta;params=params,
            return_lateral=context.extrude_return_lateral,caller=caller) :
        revolve_entities!(m,entities,revolve.axis,revolve.origin,
            revolve.angle;params=params,
            return_lateral=context.extrude_return_lateral,caller=caller)
    return Float64.(tags)
end

function _geo_transform_params(;kind::AbstractString,params::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString)
    if kind=="Translate"
        return _affine_translation(Tuple(_geo_exec_numeric_values(
            params,3,context,"$caller delta")),caller)
    elseif kind=="Symmetry"
        return _affine_symmetry(_geo_exec_numeric_values(
            params,4,context,"$caller plane coefficients")...,caller)
    elseif kind=="Dilate"
        parts=_geo_split_top_commas(params,caller)
        length(parts)==2 || throw(ArgumentError(
            "$caller: Dilate requires `{center, scale}` parameters"))
        (cinner,crest)=_geo_balanced_group(parts[1],caller)
        isempty(crest) || throw(ArgumentError(
            "$caller: Dilate center must be a single `{...}` group"))
        center=Tuple(_geo_exec_numeric_values(
            cinner,3,context,"$caller center"))
        scale_raw=String(strip(parts[2]))
        scale=if startswith(scale_raw,"{")
            (sinner,srest)=_geo_balanced_group(scale_raw,caller)
            isempty(srest) || throw(ArgumentError(
                "$caller: Dilate scales must be a single `{...}` group"))
            Tuple(_geo_exec_numeric_values(
                sinner,3,context,"$caller scales"))
        else
            _geo_eval_numeric(scale_raw,context,"$caller scale")
        end
        return _affine_dilation(center,scale,caller)
    else
        parts=_geo_split_top_commas(params,caller)
        length(parts)==3 || throw(ArgumentError(
            "$caller: Rotate requires `{{axis}, {origin}, angle}` parameters"))
        (ainner,arest)=_geo_balanced_group(parts[1],caller)
        (oinner,orest)=_geo_balanced_group(parts[2],caller)
        (isempty(arest) && isempty(orest)) || throw(ArgumentError(
            "$caller: Rotate axis and origin must be single `{...}` groups"))
        axis=Tuple(_geo_exec_numeric_values(ainner,3,context,"$caller axis"))
        origin=Tuple(_geo_exec_numeric_values(oinner,3,context,"$caller origin"))
        angle=_geo_eval_numeric(parts[3],context,"$caller angle")
        return _affine_rotation(axis,origin,angle,caller)
    end
end

function _geo_exec_transform_statement!(m::GeoModel,line::AbstractString,
                                        context::_GeoNumericContext)
    source=String(strip(line))
    endswith(source,";") && (source=String(strip(source[1:prevind(source,end)])))
    mm=match(r"^(Translate|Rotate|Dilate|Symmetry)\s*",source)
    kind=String(mm.captures[1])
    caller="execute_geo: $kind"
    rest=String(strip(source[nextind(source,firstindex(source),
                                   ncodeunits(kind)):end]))
    (params,rest)=_geo_balanced_group(rest,caller)
    (shape_list,rest)=_geo_balanced_group(rest,caller)
    isempty(rest) || throw(ArgumentError(
        "$caller: unexpected text after the shape list"))
    t=_geo_transform_params(kind=kind,params=params,context=context,caller=caller)
    entities=_geo_shape_list_entities!(m,shape_list,context,caller)
    transform_entities!(m,t,entities;caller=caller)
    return nothing
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
    # Gmsh parses every entity RHS as `ListOfDouble`: `{...}` groups, bare
    # `name[]`/`name[{..}]` references, `-{...}` negation, `expr * {...}`
    # multipliers and plain scalars — `_geo_numeric_list_values` covers the
    # whole grammar and reports malformed brace forms itself.
    source=String(strip(raw))
    isempty(source) && throw(ArgumentError("$caller: entity list must not be empty"))
    return _geo_exec_entity_tags(
        source,context,caller;signed=signed,wrap=false)
end

_geo_periodic_tags(raw::AbstractString,context::_GeoNumericContext,
                   caller::AbstractString)=
    _geo_exec_entity_tags(raw,context,caller)

# Evaluate a Gmsh `VExpr` — `{a,b,c[,d[,e]]}` or `(a,b,c)` groups composed by
# unary/binary `+`/`-` — returning the first three components (all consumers
# here read only 0..2, matching `addCircleArc`'s use of `CircleOptions`).
function _geo_exec_vexpr3(raw::AbstractString,context::_GeoNumericContext,
                          caller::AbstractString)
    s=String(strip(raw))
    isempty(s) && throw(ArgumentError(
        "$caller: expected a vector expression"))
    acc=(0.0,0.0,0.0);sign=1.0;expect_term=true
    while !isempty(s)
        if expect_term
            while startswith(s,"+") || startswith(s,"-")
                s[1]=='-' && (sign=-sign)
                s=String(strip(s[nextind(s,firstindex(s)):end]))
            end
            isempty(s) && throw(ArgumentError(
                "$caller: malformed vector expression $(repr(raw))"))
            group,rest=if s[1]=='{'
                _geo_balanced_group(s,caller)
            elseif s[1]=='('
                _geo_balanced_paren(s,caller)
            else
                throw(ArgumentError(
                    "$caller: expected a `{...}` or `(...)` vector group; " *
                    "got $(repr(s))"))
            end
            parts=_geo_split_top_commas(group,caller)
            length(parts) in 3:5 || throw(ArgumentError(
                "$caller: a vector expression needs 3 to 5 components; got " *
                "$(length(parts))"))
            acc=acc .+ sign .* ntuple(
                i->_geo_eval_numeric(parts[i],context,caller),3)
            s=String(strip(rest));sign=1.0;expect_term=false
        else
            (s[1]=='+' || s[1]=='-') || throw(ArgumentError(
                "$caller: expected `+` or `-` in vector expression; got " *
                "$(repr(s))"))
            sign=s[1]=='+' ? 1.0 : -1.0
            s=String(strip(s[nextind(s,firstindex(s)):end]))
            expect_term=true
        end
    end
    expect_term && throw(ArgumentError(
        "$caller: vector expression ends with an operator"))
    return acc
end

# Parse the shared `Circle`/`Ellipse` RHS: `{point tags}` plus the optional
# `Plane VExpr` normal override (`CircleOptions` in the Gmsh grammar).
function _geo_circle_rhs(raw::AbstractString,context::_GeoNumericContext,
                         caller::AbstractString)
    (group,rest)=_geo_balanced_group(String(strip(raw)),caller)
    points=_geo_exec_entity_tags(group,context,"$caller points")
    isempty(rest) && return (points,nothing)
    match(r"^Plane\b",rest)===nothing && throw(ArgumentError(
        "$caller: expected a `Plane` option after the point list; got " *
        "$(repr(rest))"))
    normal=_geo_exec_vexpr3(
        String(strip(rest[nextind(rest,firstindex(rest),5):end])),
        context,"$caller Plane")
    return (points,normal)
end

# Parse the optional `Surface`/`Ruled Surface` trailing constraint
# (`SurfaceConstraints`): `In Sphere{p}` or a single-element `Using Point{p}`
# supply the sphere-center point tag; `Using GeoEntity{..}` and multi-point
# lists are ignored by Gmsh's built-in kernel. Returns the center tag or
# `nothing` (a nonpositive value is ignored exactly like Gmsh's
# `sphereCenterTag >= 0` guard).
function _geo_surface_constraint(rest::AbstractString,
                                 context::_GeoNumericContext,
                                 caller::AbstractString)
    s=String(strip(rest))
    isempty(s) && return nothing
    if (mm=match(r"^In\s+Sphere\b",s)) !== nothing
        (group,r2)=_geo_balanced_group(
            String(strip(s[nextind(s,firstindex(s),ncodeunits(mm.match)):end])),
            caller)
        isempty(r2) || throw(ArgumentError(
            "$caller: unexpected text after In Sphere constraint"))
        parts=_geo_split_top_commas(group,caller)
        length(parts)==1 || throw(ArgumentError(
            "$caller: In Sphere takes a single expression"))
        raw=_geo_signed_gmsh_int_value(
            _geo_eval_numeric(parts[1],context,"$caller In Sphere"),
            "$caller In Sphere")
        return raw>=0 ? raw : nothing
    elseif (mm=match(r"^Using\s+(Point|GeoEntity)\b",s)) !== nothing
        (group,r2)=_geo_balanced_group(
            String(strip(s[nextind(s,firstindex(s),ncodeunits(mm.match)):end])),
            caller)
        isempty(r2) || throw(ArgumentError(
            "$caller: unexpected text after Using $(mm.captures[1]) constraint"))
        isempty(strip(group)) && return nothing
        values=_geo_numeric_list_values(
            "{"*group*"}",context,"$caller Using $(mm.captures[1])")
        mm.captures[1]=="Point" && length(values)==1 || return nothing
        raw=_geo_signed_gmsh_int_value(values[1],"$caller Using Point")
        return raw>=0 ? raw : nothing
    end
    throw(ArgumentError(
        "$caller: unsupported surface constraint $(repr(s)); expected " *
        "`In Sphere{..}` or `Using Point{..}`"))
end

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
            r"^Circle\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Circle"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        (points,normal)=_geo_circle_rhs(mm.captures[2],context,caller)
        length(points)==3 || throw(ArgumentError(
            "$caller: Circle requires 3 points"))
        add_circle_arc!(m,points[1],points[2],points[3];
                        tag=tag,plane_normal=normal)
        return
    elseif (mm=match(
            r"^Ellipse\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Ellipse"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        (points,normal)=_geo_circle_rhs(mm.captures[2],context,caller)
        length(points) in (3,4) || throw(ArgumentError(
            "$caller: Ellipse requires 4 points"))
        # The built-in 3-tag form mirrors the OCC backward-compatibility
        # record: the start point doubles as the major-axis point
        # (`addEllipseArc(num, tags[0], tags[1], tags[0], tags[2], ..)`).
        length(points)==3 &&
            (points=[points[1],points[2],points[1],points[3]])
        add_ellipse_arc!(m,points[1],points[2],points[3],points[4];
                         tag=tag,plane_normal=normal)
        return
    elseif (mm=match(
            r"^Spline\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Spline"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller control points")
        add_spline!(m,points;tag=tag)
        return
    elseif (mm=match(
            r"^BSpline\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: BSpline"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller control points")
        add_bspline!(m,points;tag=tag)
        return
    elseif (mm=match(
            r"^Bezier\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Bezier"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller control points")
        add_bezier!(m,points;tag=tag)
        return
    elseif (mm=match(
            r"^Nurbs\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        # `tNurbs (FExpr) = ListOfDouble tNurbsKnots ListOfDouble tNurbsOrder
        # FExpr` — the built-in kernel routes through `addBSpline(num, tags,
        # seqknots)` and never reads the parsed `Order` expression.
        caller="execute_geo: Nurbs"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        split_knots=_geo_split_at_keyword(mm.captures[2],"Knots",caller)
        split_knots===nothing && throw(ArgumentError(
            "$caller: expected `Knots <list>` after the control-point list; " *
            "got $(repr(strip(mm.captures[2])))"))
        (points_raw,knots_tail)=split_knots
        points=_geo_exec_entity_rhs_tags(
            points_raw,context,"$caller control points")
        split_order=_geo_split_at_keyword(knots_tail,"Order",caller)
        split_order===nothing && throw(ArgumentError(
            "$caller: expected `Order <expr>` after the knot list; got " *
            "$(repr(knots_tail))"))
        (knots_raw,order_raw)=split_order
        isempty(order_raw) && throw(ArgumentError(
            "$caller: Nurbs requires an `Order` expression"))
        knots=_geo_numeric_list_values(knots_raw,context,"$caller knots")
        _geo_eval_numeric(order_raw,context,"$caller Order")
        add_nurbs!(m,points,knots;tag=tag)
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
            r"^(Ruled\s+)?Surface\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        # `Surface` and the deprecated `Ruled Surface` alias both execute
        # Gmsh's `addSurfaceFilling`: the first wire's curve count selects the
        # patch kind (4 -> ruled, 3 -> triangular).
        caller="execute_geo: $(mm.captures[1]===nothing ? "" : "Ruled ")Surface"
        tag=_geo_exec_entity_tag(mm.captures[2],context,"$caller tag")
        (group,rest)=_geo_balanced_group(
            String(strip(mm.captures[3])),caller)
        wires=_geo_exec_entity_tags(group,context,"$caller loops")
        sphere_center=_geo_surface_constraint(rest,context,caller)
        add_ruled_surface!(m,wires;tag=tag,sphere_center=sphere_center)
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
            r"^Torus\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Torus"
        tag=_geo_exec_entity_tag(mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        (length(values)==5 || length(values)==6) || throw(ArgumentError(
            "$caller: Torus requires 5 or 6 numeric parameters; " *
            "got $(length(values)) after range expansion"))
        add_torus!(m,values[1:5]...;tag=tag,
                   angle=length(values)==6 ? values[6] : 2π)
        return
    elseif (mm=match(
            r"^Boolean(Difference|Union|Intersection|Fragments)\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Boolean$(mm.captures[1])"
        tag=_geo_exec_entity_tag(mm.captures[2],context,"$caller result tag")
        groups=_geo_boolean_groups(mm.captures[3],caller)
        operands=_geo_boolean_operands(groups,context,caller)
        boolean_volumes_multi!(m,_GEO_BOOLEAN_OPS[mm.captures[1]],
            operands.objects.tags,operands.tools.tags;tag=tag,
            remove_object=operands.objects.delete,
            remove_tool=operands.tools.delete,caller=caller)
        return
    elseif match(r"^(Translate|Rotate|Dilate|Symmetry)\s*\{",line)!==nothing
        _geo_exec_transform_statement!(m,line,context)
        return
    elseif match(r"^Affine\s*\{",line)!==nothing
        throw(ArgumentError(
            "execute_geo: Affine transforms require the OpenCASCADE geometry " *
            "kernel, which Tessella does not implement"))
    elseif match(r"^Extrude\b",line)!==nothing
        _geo_exec_extrude_term(m,line,context)
        return
    elseif match(
            r"^(Duplicata|Boundary|CombinedBoundary|OrientedBoundary|OrientedCombinedBoundary|PointsOf)\s*\{",
            line)!==nothing
        # Standalone MultipleShape statements are legal in the Gmsh grammar;
        # Duplicata copies are the observable side effect.
        s=String(strip(line))
        endswith(s,";") &&
            (s=String(strip(s[firstindex(s):prevind(s,end)])))
        _geo_shape_list_entities!(m,s,context,"execute_geo")
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
            r"^Coherence(?:\s+(Geometry|Mesh)|\s+Point\s*\{\s*(.*?)\s*\})?\s*;$",
            line)) !== nothing
        if mm.captures[2] !== nothing
            merge_vertices!(m,_geo_exec_entity_tags(
                mm.captures[2],context,"execute_geo: Coherence Point");
                caller="execute_geo: Coherence Point")
        elseif mm.captures[1] != "Mesh"
            coherence!(m)
        end
        # `Coherence Mesh` deduplicates mesh vertices; no mesh exists during
        # .geo execution, so there is nothing to merge.
        return
    elseif (mm=match(
            r"^Mesh\.TransfiniteTri\s*=\s*(.*?)\s*;$",line)) !== nothing
        value=_geo_eval_numeric(mm.captures[1],context,
                                "execute_geo: Mesh.TransfiniteTri")
        isinteger(value) || throw(ArgumentError(
            "execute_geo: Mesh.TransfiniteTri must be 0 or 1 (got $value)"))
        set_transfinite_tri!(m,Int(value))
        return Int(value)
    elseif startswith(line,"Coherence")
        throw(ArgumentError(
            "execute_geo: unknown coherence command: $line"))
    elseif startswith(line,"Geometry.ExtrudeReturnLateralEntities")
        m2=match(r"^Geometry\.ExtrudeReturnLateralEntities\s*=\s*(.*?)\s*;?\s*$",
                 line)
        m2===nothing && throw(ArgumentError(
            "unrecognized statement: $line"))
        vals=_geo_numeric_list_values(strip(m2[1]),context,
                                      "Geometry.ExtrudeReturnLateralEntities")
        length(vals)==1 || throw(ArgumentError(
            "Geometry.ExtrudeReturnLateralEntities expects a scalar value"))
        context.extrude_return_lateral=!iszero(vals[1])
        return
    elseif startswith(line,"Mesh.") || startswith(line,"SetFactory") ||
           startswith(line,"Field") || startswith(line,"Background") ||
           startswith(line,"BoundaryLayer") ||
           occursin(r"^Mesh\s+[0-9]\s*;", line)
        return
    elseif (mm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*=\s*(.*?)\s*;$",line)) !== nothing
        # `v() =` is Gmsh's ListOfDouble affectation — equivalent to `v[] =`;
        # it is the capture form for side-effecting terms such as
        # `BooleanX{...}{...}` and `Extrude`.
        name=String(mm.captures[1])
        (name=="Pi" || name in _GEO_SIDE_EFFECT_SYMBOLS) && throw(ArgumentError(
            "execute_geo: $name is read-only and cannot be assigned as a list"))
        values=_geo_numeric_list_values(mm.captures[2],context,
            "execute_geo: list variable $name";allow_multiplier=false)
        _geo_context_set_list!(context,name,values,
            "execute_geo: list variable $name")
        return
    elseif occursin(r"^[A-Za-z_][A-Za-z0-9_]*\s*\[",line)
        body=String(strip(line[firstindex(line):prevind(line,lastindex(line))]))
        _geo_apply_list_assignment!(
            context,body,"execute_geo: numeric list assignment") || throw(
                ArgumentError("execute_geo: malformed numeric list assignment: $line"))
        return
    elseif (mm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*;$",line)) !== nothing
        if match(r"^Extrude\b",strip(mm.captures[2]))!==nothing
            # `x = Extrude{..}{..}` stores the flat result list (Gmsh keeps a
            # ListOfDouble), so the name becomes a list variable.
            values=_geo_numeric_list_values(
                mm.captures[2],context,
                "execute_geo: variable $(mm.captures[1])";
                allow_multiplier=false)
            _geo_context_set_list!(context,String(mm.captures[1]),values,
                "execute_geo: variable $(mm.captures[1])")
            return
        end
        _geo_exec_scalar!(context,mm.captures[1],mm.captures[2])
        return
    end
    throw(ArgumentError("execute_geo: unrecognized statement: $line"))
end

end # module
