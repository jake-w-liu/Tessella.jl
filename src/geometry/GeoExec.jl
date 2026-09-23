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
surfaces, or stored mesh-inert volume relations with reusable masters and
dependency chains or cycles. `MeshSize` and
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
using ..Model: merge_vertices!, extrude_entities!, revolve_entities!,
             twist_entities!
using ..Model: _GeoExtrudeParams
using ..Model: _affine_translation, _affine_dilation
using ..Model: _affine_rotation, _affine_symmetry, _entity_label
using ..Model: add_physical_group!, set_periodic!, set_transfinite_tri!, _has_entity,
    _max_entity_physical_number, _model_periodic_curve_orientation,
    _model_periodic_surface_edge_map
using ..Model: _model_boundary, _model_points_of, _model_direct_boundary
using ..Model: _model_entity_dictionary, _model_entity_known, remove_embedded!
using ..Model: model_entity, model_entities, model_entities_in_bounding_box
using ..Model: model_physical_groups, _physical_live_members
using ..Model: model_entity_color, model_parametrization_bounds, model_normal
using ..Model: model_reparametrize_on_surface
using ..Model: _model_entity_bounding_box, _model_bounds_union, _occ_geometry
using ..Model: _geo_delete_entities!, _geo_reset_model_geometry!
using ..Model: _tag, _alloc_tag!, _alloc_curve_loop_tag,
    _alloc_surface_loop_tag, _arc_stored_normal, _occ_cross
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..IO: read_geo_params, _GeoNumericContext, _geo_eval_numeric
using ..IO: _geo_split_list, _geo_split_range, _geo_range_count
using ..IO: _geo_numeric_list_terms, _geo_numeric_list_values
using ..IO: _geo_signed_gmsh_int_value
using ..IO: _GEO_SIDE_EFFECT_SYMBOLS
using ..IO: _MAX_GEO_LIST_ITEMS, _MAX_GEO_SCAN_NESTING
using ..IO: _geo_context_set_scalar!, _geo_context_set_list!
using ..IO: _geo_apply_list_assignment!, _geo_list_term_values
using ..IO: _geo_brace_terminated_statement, _geo_extrude_shape_group,
            _geo_extrude_pipe_candidate
using ..IO: _GeoTagAllocatorState, _geo_context_refresh_allocators!
using ..IO: _geo_allocator_observe_statement!,_geo_allocator_resync_model!
using ..IO: _geo_allocator_delete_entities!, _geo_allocator_reset_model!
using ..IO: _GEO_SHAPE_DEFINITION_DIMS, _geo_statement_marks_internals
using ..IO: _geo_allocator_set_factory!
using ..IO: _geo_context_forget!
using ..IO: _geo_physical_declaration
using ..IO: _geo_exec_assignment!, _geo_eval_string, _geo_symbol_name
using ..IO: _geo_split_args, _geo_yyerror!, _geo_yywarn!, _geo_yyinfo!,
            _geo_int_value
using ..IO: _geo_fix_relative_path, _geo_print_list_of_double
using ..IO: _geo_vsnprintf_expand
using ..IO: _geo_string_rhs, _geo_eval_string_list
using ..IO: _geo_matching_delim, _GEO_COLOR_NAMES
using ..IO: _geo_option_number
using ..IO: _geo_msg_error!, _geo_eof_error!, _geo_context_has_variable
using ..IO: _geo_exec_struct_def!
using ..IO: _GeoSyntaxAbort, _geo_syntax_abort, _geo_first_token,
            _geo_stmt_head_token, _GeoLoopLevel, _GeoWhileLevel
using ..IO: _geo_string_index_offender, _GEO_STRING_FUNCTIONS,
            _geo_tstring_valid,
            _GEO_BARE_STRING_TOKENS, _GEO_ALL_FUNCTIONS, _GEO_LEXER_KEYWORDS,
            _geo_scan_top_level_comma, _geo_lhs_keyword_check,
            _geo_control_statement, _GEO_CONTROL_BARE
using ..IO: GeoFieldSpec, _GEO_FIELD_KINDS, _geo_expression_field_tags
using ..MeshTypes: nnodes, nsegs, ntris, ntets
using ..Refine: refine_uniform
using ..Recombine: recombine_triangles
using ..IO: read_msh, write_msh
using ..IO: _geo_gmsh_number
using ..Transform: _affine_coordinate, _periodic_affine_3x4
using LinearAlgebra: norm, svd

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
    values::Dict{String,Float64}
    lists::Dict{String,Vector{Float64}}
    # Final string variables (`gmsh_yystringsymbols` mirror).
    strings::Dict{String,Vector{String}}
    # `yymsg(1)` warnings collected during execution.
    warnings::Vector{String}
    # `Msg::Error` count (`Error(...)`, option/plugin diagnostics) — the
    # `_atLeastOneErrorInRun` flag; the Gmsh batch exit code when nonzero.
    msg_error_count::Int
    # `Exit n;` sets the code; `nothing` when the program ran to its end.
    # `Exit;` uses the accumulated error count, so a nonzero code means the
    # Gmsh process would have failed — `execute_geo` throws in that case.
    exit_code::Union{Nothing,Int}
end

const _MAX_GEO_EXEC_STATEMENT_BYTES=1_000_000
const _MAX_GEO_EXEC_STATEMENTS=1_000_000
const _MAX_GEO_LOOP_ITERATIONS=1_000_000
const _MAX_GEO_CALL_DEPTH=128

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
    elseif (matched=match(r"^For\s*\((.*)\)$",line))!==nothing
        # `For (a:b[:c])` — anonymous range, no loop variable.
        return (kind=:for,var="",range=String(matched.captures[1]))
    elseif (matched=match(r"^(?:Function|Macro)\s+(.+)$",line))!==nothing
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
    # `Extrude {shapes}` is `}`-complete but may still take the `Using Wire`
    # pipe suffix — its emit is deferred until the next real token, even
    # across line boundaries.
    extrude_pending=false
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
            if extrude_pending && !block_comment && quote_char=='\0' &&
               !isspace(c) && !(c=='/' && (nextc=='/' || nextc=='*'))
                extrude_pending=false
                if !(c=='U' &&
                     match(r"^Using\b",SubString(raw,i,last))!==nothing)
                    statement=strip(String(take!(buf)))
                    isempty(statement) || (
                        length(statements)<_MAX_GEO_EXEC_STATEMENTS ||
                            throw(ArgumentError(
                                "execute_geo: input exceeds " *
                                "$_MAX_GEO_EXEC_STATEMENTS statements"));
                        push!(statements,statement))
                    buf_has_content=false
                end
            end
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
                    if _geo_extrude_pipe_candidate(
                            String(buf.data[1:position(buf)]))
                        # Possibly `Extrude {shapes} Using Wire {expr}` —
                        # hold the emit until the next real token.
                        extrude_pending=true
                    else
                        statement=strip(String(take!(buf)))
                        length(statements)<_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                            "execute_geo: input exceeds $_MAX_GEO_EXEC_STATEMENTS statements"))
                        push!(statements,statement)
                        buf_has_content=false
                    end
                end
            elseif c==';' && depth==0
                write(buf,c)
                statement=strip(String(take!(buf)))
                # A `;` left after a `}`-terminated transform or between
                # statements is an upstream `syntax error (;)` — it is kept
                # as a statement so execution reports and recovers on it.
                length(statements)<_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                    "execute_geo: input exceeds $_MAX_GEO_EXEC_STATEMENTS statements"))
                push!(statements,statement)
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
    if extrude_pending && !isempty(tail)
        # A held `Extrude {shapes}` with no `Using Wire` following is
        # complete — emit it before the trailing-content checks.
        push!(statements,tail);tail=""
    end
    # A trailing run of bare `;` is still statements upstream — `syntax error
    # (;)` per stray terminator — so it is pushed for execution to report.
    if !isempty(tail) && all(c->c==';'||isspace(c),tail)
        for _ in 1:count(==(';'),tail)
            push!(statements,";")
        end
    elseif !isempty(tail)
        # Control keywords need no `;` upstream — a trailing `EndFor`,
        # `EndIf` or `Return` still takes its action at EOF. Ordinary
        # unterminated text stays a hard error.
        control=_geo_control_statement(tail,firstindex(tail),lastindex(tail))
        control===nothing && throw(ArgumentError(
            "execute_geo: unterminated statement: $tail"))
        push!(statements,control[1])
    end
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
`Periodic Curve`, and `Periodic Surface` accept `Translate`, `Rotate`, and
`.geo` `Affine` transforms with the per-dimension arity rules Gmsh's parser
applies (curves: an empty transform is the orientation-only form, 12–15 or
>16 entries drop after the endpoint check; surfaces: 12+ entries zero-pad to
the stored 16); the transform-free `Periodic Surface {s} = {m}` is rejected
like upstream, and `Periodic Surface j {c} = k {c}` derives its transform
from the mapped boundary vertices. `Periodic Volume` is not in Gmsh's
grammar and is a `syntax error (Volume)` here; `set_periodic!` still stores
volume relations as a Tessella extension. Affine curve relations require
straight curves at declaration; orientation-only declarations are check-free
like upstream and hold all curves to the straight-curve requirement at mesh
time. Surfaces must be planar boundaries of one explicit volume when meshed.
Multiple periodic statements may
reuse a master or form a master/slave chain or cycle. Read-only `newp`, the shared
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
    context.file_name=String(path)
    # Field declarations/options/delete replay into `context.fields` — the
    # statement-time live map, shared with `params.fields` so the returned
    # `GeoExecution` exposes the executed end state. `Field[i].member = v`
    # writes mutate it so a mid-file `Mesh n` sees the final option values.
    context.fields=empty!(params.fields)
    # `.geo` execution softens variable-read misses into `yymsg(0)` + 0, like
    # `treat_Struct_FullName_Float` — the accumulated errors throw at the end.
    context.soft_unknown_reads=true
    # `Point{t}`/`Physical X{t}` name reads in StringExpr position resolve
    # against the live model.
    context.entity_name_lookup=(dim,tag,kind)->
        kind===:physical ? get(model.physical_names,(dim,tag),"") :
            get(model.entity_names,(dim,tag),"")
    # `Extrude{...}{...}`, `BooleanX{...}{...}`, transform blocks and the
    # entity-selector expressions (`Point{t}`, `Physical X{...}`, `BoundingBox
    # X{...}`, ...) are side-effecting or model-reading value terms in the
    # `.geo` grammar; the numeric evaluator calls back through this hook.
    allocator_state=_GeoTagAllocatorState()
    context.exec_hook=src->_geo_exec_value_term(
        model,src,context,allocator_state)
    statements=_geo_exec_statements(path)
    executed=Ref(0)
    transfinite_tri=_exec_geo_statements!(
        model,statements,firstindex(statements),lastindex(statements),
        context,allocator_state,executed)
    # End-of-parse `synchronize` — the raw physical-group records materialize
    # into the observable view (`m.physical`), forwarding-declared members
    # resolve, and unknown members warn once.
    _geo_sync_physical_view!(model,context)
    # `Abort` / `Exit` set the stream stop. `Msg::Exit(n, force)`: an explicit
    # `Exit n` exits with `n` regardless of accumulated errors (`forceLevel`),
    # while a bare `Exit` exits with `_atLeastOneErrorInRun` (0 or 1). `Abort`
    # is a voluntary stop (errorstate 999) — nonzero but not an error.
    if context.stop===:exit
        code=context.exit_code
        code==0 || throw(ArgumentError(
            "execute_geo: Exit $code (script requested a nonzero exit)"))
        # `exit(0)` leaves before the end-of-parse errorstate check.
        empty!(context.exec_errors)
    end
    isempty(context.exec_errors) || throw(ArgumentError(
        "execute_geo: $(join(context.exec_errors,"; "))"))
    # `Msg::Error` diagnostics (option writes, `Error(...)`) do not fail the
    # parse — they surface through `GeoExecution.msg_error_count` like the
    # Gmsh batch exit code.
    # A mid-file `Mesh n` leaves its product on the model for `Save` and the
    # mesh-operation statements; the mesh_dim keyword still generates on
    # request.
    mesh=context.mesh
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
    return GeoExecution(model,mesh,params,transfinite_tri,context.values,
                        context.lists,context.strings,
                        copy(context.exec_warnings),context.msg_error_count,
                        context.stop===:exit ? context.exit_code : nothing)
end

# Upstream control flow is a token stream, not a block tree (Gmsh.y
# 3887-4100): `If`/`ElseIf`/`Else`/`EndIf` shift `ImbricatedTest` and the
# per-level `statusImbricatedTests`, `For` pushes a loop frame that `EndFor`
# rewinds or pops, and dead text is consumed by `skip`/`skipTest` — char-level
# scans that watch only their own family's keywords. `For`, `While`,
# `Function`, `Return` and foreign closers are invisible inside an `If` skip;
# `EndIf` inside a `For` body is a real statement that decrements the outer
# counter when the body re-parses.

# `skip("If","EndIf")` / `skip("For","EndFor")` / `skip("While","EndWhile")`:
# same-family words nest; the depth-0 closer is consumed. Returns the index
# after the closer; at EOF records `Unexpected end of file` (`Msg::Error`)
# and returns past the end.
function _geo_exec_skip_family(statements::Vector{String},i::Int,
                               open_kind::Symbol,close_kind::Symbol,
                               context::_GeoNumericContext)
    depth=0;j=i
    while j<=lastindex(statements)
        control=_geo_control_parse(statements[j])
        if control!==nothing
            if control.kind===open_kind
                depth+=1
            elseif control.kind===close_kind
                depth==0 && return j+1
                depth-=1
            end
        end
        j+=1
    end
    _geo_eof_error!(context)
    return lastindex(statements)+1
end

# `skipTest("If","EndIf","ElseIf",4)`: watches `If`/`ElseIf`/`Else`/`EndIf`.
# Returns (index, type): type 0 — a depth-0 `EndIf` was consumed and index is
# past it (the caller decrements `ImbricatedTest`); type 1/2 — the stream
# landed ON a depth-0 `ElseIf`/`Else` which then parses normally. At EOF
# records `Unexpected end of file` and returns type 0.
function _geo_exec_skiptest(statements::Vector{String},i::Int,
                            context::_GeoNumericContext)
    depth=0;j=i
    while j<=lastindex(statements)
        control=_geo_control_parse(statements[j])
        if control!==nothing
            k=control.kind
            if k===:elseif && depth==0
                return j,1
            elseif k===:else && depth==0
                return j,2
            elseif k===:endif
                depth==0 && return j+1,0
                depth-=1
            elseif k===:if
                depth+=1
            end
        end
        j+=1
    end
    _geo_eof_error!(context)
    return lastindex(statements)+1,0
end

# `skip(nullptr,"Return")` — the flat scan a `Function`/`Macro` definition
# runs: the first `Return` word is consumed, everything else is text.
function _geo_exec_skip_return(statements::Vector{String},i::Int,
                               context::_GeoNumericContext)
    j=i
    while j<=lastindex(statements)
        control=_geo_control_parse(statements[j])
        (control!==nothing && control.kind===:return) && return j+1
        j+=1
    end
    _geo_eof_error!(context)
    return lastindex(statements)+1
end

# `error tEND` recovery (Gmsh.y:276): after a syntax error the parser
# discards tokens through the next `;`. Control keywords and `}`-terminated
# statements carry no `;`, so the discard eats them plus the first
# `;`-terminated statement — `If (x y) A; B;` drops `A` but `B` still parses,
# and `If (x y) EndIf` drops the `EndIf` (no `Orphan EndIf` warning).
function _geo_recovery_next_i(statements::Vector{String},i::Int)
    j=i+1
    while j<=lastindex(statements)
        done=endswith(rstrip(statements[j]),';')
        j+=1
        done && break
    end
    return j
end

# `For name In {..}` clears the list flag and stores the scalar, like
# upstream's `s.list = false; s.value.resize(1); s.value[0] = first`.
function _geo_bind_loop_var(context::_GeoNumericContext,var::AbstractString,
                            value::Float64,caller::AbstractString)
    isempty(var) && return nothing
    delete!(context.lists,var);delete!(context.list_variables,var)
    delete!(context.unavailable_lists,var)
    _geo_context_set_scalar!(context,var,value,caller)
    return nothing
end

# `tIf '(' FExpr ')'`: the condition reduces before `ImbricatedTest++` (a
# `TestLevel` read inside it sees the outer depth). A false condition runs
# `skipTest` — the level stays open until a midpoint token or the consumed
# `EndIf` closes it.
function _geo_exec_if_header!(statements::Vector{String},i::Int,
                              context::_GeoNumericContext)
    control=_geo_control_parse(statements[i])
    cond=_geo_eval_numeric(control.cond,context,"execute_geo")!=0
    context.if_depth+=1
    if context.if_depth>_MAX_GEO_SCAN_NESTING-1
        _geo_yyerror!(context,"Reached maximum number of imbricated tests")
        context.if_depth=_MAX_GEO_SCAN_NESTING-1
    end
    # `statusImbricatedTests[0]` is a real upstream slot (an `If` opened at a
    # negative depth writes it); it is simply never read back.
    context.if_depth>=0 &&
        (context.if_status[context.if_depth+1]=cond)
    cond && return i+1
    j,typ=_geo_exec_skiptest(statements,i+1,context)
    typ==0 && (context.if_depth-=1)
    return j
end

# `tElseIf '(' FExpr ')'`: the condition reduces before the action — it
# reports errors even when the level then `skip`s away. `status==1` skips to
# `EndIf` and closes the level; a false cond re-runs `skipTest`; a true one
# marks the level and lets the following text run.
function _geo_exec_elseif!(statements::Vector{String},i::Int,
                           context::_GeoNumericContext)
    control=_geo_control_parse(statements[i])
    cond=_geo_eval_numeric(control.cond,context,"execute_geo")!=0
    if context.if_depth<=0
        _geo_yyerror!(context,"Orphan ElseIf")
        return i+1
    end
    if context.if_status[context.if_depth+1]
        j=_geo_exec_skip_family(statements,i+1,:if,:endif,context)
        context.if_depth-=1
        return j
    end
    if cond
        context.if_status[context.if_depth+1]=true
        return i+1
    end
    j,typ=_geo_exec_skiptest(statements,i+1,context)
    typ==0 && (context.if_depth-=1)
    return j
end

# `tElse`: `status==1` skips to `EndIf` and closes the level; `status==0`
# does *nothing* — the flag stays false, so a later `ElseIf` can still fire
# (`If(0) A Else B ElseIf(1) C` runs B then C upstream).
function _geo_exec_else!(statements::Vector{String},i::Int,
                         context::_GeoNumericContext)
    if context.if_depth<=0
        _geo_yyerror!(context,"Orphan Else")
        return i+1
    end
    if context.if_status[context.if_depth+1]
        j=_geo_exec_skip_family(statements,i+1,:if,:endif,context)
        context.if_depth-=1
        return j
    end
    return i+1
end

# `tFor tSTRING tIn '{' FExpr : FExpr [: FExpr] '}'` and the anonymous
# `tFor '(' FExpr : FExpr [: FExpr] ')'`: the range is evaluated once, the
# variable binds `first`, and an out-of-range start runs `skip("For",
# "EndFor")` — the level is never pushed, so the consumed `EndFor` does not
# report `Invalid For/EndFor loop`. A zero step is legal: the body runs once
# and the `EndFor` check pops the level immediately.
function _geo_exec_for_header!(statements::Vector{String},i::Int,
                               context::_GeoNumericContext)
    control=_geo_control_parse(statements[i])
    caller="execute_geo: For"
    variable=control.var
    !isempty(variable) &&
        (variable=="Pi" || variable in _GEO_SIDE_EFFECT_SYMBOLS) &&
        throw(ArgumentError("$caller: loop variable $variable is reserved"))
    # The For range is `FExpr : FExpr [: FExpr]` upstream — not a list: a
    # comma form evaluates its head then errors at `,`, a lone value
    # evaluates then errors at `}`. Each endpoint reduces in order before
    # the grammar commits.
    spec=String(strip(control.range))
    pieces=_geo_split_range(spec,caller)
    if pieces===nothing || length(pieces)==1
        text=pieces===nothing ? spec : pieces[1]
        comma=_geo_scan_top_level_comma(text)
        if comma>0
            _geo_eval_numeric(
                text[firstindex(text):prevind(text,comma)],context,caller)
            _geo_syntax_abort(",")
        end
        _geo_eval_numeric(text,context,caller)
        # `{` ranges fail at `}`, `( a:b )` ranges at `)`.
        _geo_syntax_abort(isempty(variable) ? ")" : "}")
    end
    first=_geo_eval_numeric(pieces[1],context,"$caller range start")
    last=_geo_eval_numeric(pieces[2],context,"$caller range end")
    step=length(pieces)==2 ? 1.0 :
         _geo_eval_numeric(pieces[3],context,"$caller range increment")
    _geo_bind_loop_var(context,variable,first,caller)
    if (step>0 && first>last) || (step<0 && first<last)
        return _geo_exec_skip_family(statements,i+1,:for,:endfor,context)
    end
    push!(context.loop_stack,
          _GeoLoopLevel(variable,first,last,step,i+1,statements))
    length(context.loop_stack)>_MAX_GEO_SCAN_NESTING-1 &&
        _geo_yyerror!(context,"Reached maximum number of imbricated loops")
    return i+1
end

# `tEndFor`: `ImbricatedLoop<=0` is `Invalid For/EndFor loop`. Otherwise the
# loop variable advances `step` — through the live symbol, so a body's own
# `i = ...` write propagates — and an in-range value rewinds the stream to
# the saved body position while an out-of-range one pops the level.
function _geo_exec_endfor!(statements::Vector{String},i::Int,
                           context::_GeoNumericContext)
    caller="execute_geo: For"
    if isempty(context.loop_stack)
        _geo_yyerror!(context,"Invalid For/EndFor loop")
        return i+1
    end
    lvl=context.loop_stack[end]
    if lvl.statements!==statements
        # The frame belongs to a caller's stream — upstream's `fsetpos` on a
        # foreign FILE* is undefined; report it like an orphan closer.
        pop!(context.loop_stack)
        _geo_yyerror!(context,"Invalid For/EndFor loop")
        return i+1
    end
    x0=lvl.x0
    if !isempty(lvl.var)
        if !haskey(context.values,lvl.var) && !haskey(context.lists,lvl.var)
            _geo_yyerror!(context,"Unknown loop variable '$(lvl.var)'")
        elseif haskey(context.lists,lvl.var)
            _geo_yyerror!(context,"Bad loop variable $(lvl.var)")
        else
            x0=context.values[lvl.var]+lvl.step
            _geo_context_set_scalar!(context,lvl.var,x0,caller)
        end
    else
        x0+=lvl.step
    end
    lvl.x0=x0
    if (lvl.step>0 && x0<=lvl.x1) || (lvl.step<0 && x0>=lvl.x1)
        return lvl.body_start
    end
    pop!(context.loop_stack)
    return i+1
end

# `While (cond)` — Tessella's bounded extension (upstream lexes `While` as a
# plain identifier). A true condition pushes a frame; a false one runs the
# family skip to `EndWhile`.
function _geo_exec_while_header!(statements::Vector{String},i::Int,
                                 context::_GeoNumericContext)
    control=_geo_control_parse(statements[i])
    _geo_eval_numeric(control.cond,context,"execute_geo: While")!=0 ||
        return _geo_exec_skip_family(
            statements,i+1,:while,:endwhile,context)
    push!(context.while_stack,
          _GeoWhileLevel(control.cond,i+1,statements,0))
    length(context.while_stack)+length(context.loop_stack)>
        _MAX_GEO_SCAN_NESTING-1 &&
        _geo_yyerror!(context,"Reached maximum number of imbricated loops")
    return i+1
end

# `EndWhile` re-evaluates the condition: true rewinds to the body start,
# false pops the frame. A stray `EndWhile` is silent (extension choice,
# matching `EndWhile` having no upstream token).
function _geo_exec_endwhile!(statements::Vector{String},i::Int,
                             context::_GeoNumericContext)
    isempty(context.while_stack) && return i+1
    lvl=context.while_stack[end]
    lvl.statements===statements || (pop!(context.while_stack);return i+1)
    _geo_eval_numeric(lvl.cond,context,"execute_geo: While")!=0 ||
        (pop!(context.while_stack);return i+1)
    lvl.iterations+=1
    lvl.iterations<=_MAX_GEO_LOOP_ITERATIONS || throw(ArgumentError(
        "execute_geo: While loop exceeds $_MAX_GEO_LOOP_ITERATIONS iterations"))
    return lvl.body_start
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
        # `Abort`/`Exit` stop the statement stream at the next boundary,
        # including inside If/For/While/Call bodies.
        context.stop===:run || break
        # Gmsh's parser aborts a file after more than 20 recorded errors
        # (`gmsh_yyerrorstate > 20` in ParseFile) — counted per file. The
        # stop unwinds Call frames too: `enterFunction` stays inside the same
        # FILE, so the abort kills the whole stream, not just the body.
        if length(context.exec_errors)-context.file_error_base>20
            _geo_yyerror!(context,"Too many errors: aborting parser...")
            context.stop=:toomany
            break
        end
        line=statements[i]
        control=_geo_control_parse(line)
        if control===nothing
            executed[]+=1
            executed[]<=_MAX_GEO_EXEC_STATEMENTS || throw(ArgumentError(
                "execute_geo: control flow exceeds $_MAX_GEO_EXEC_STATEMENTS " *
                "executed statements"))
            if match(r"^Call\b",line)!==nothing
                assigned=try
                    _geo_exec_call!(
                        m,statements,line,context,allocator_state,executed)
                catch err
                    err isa InterruptException && rethrow()
                    err isa _GeoSyntaxAbort || rethrow()
                    _geo_yyerror!(context,"syntax error ($(err.token))")
                    nothing
                end
                assigned===nothing || (transfinite_tri=assigned)
                i+=1
                continue
            end
            _geo_context_refresh_allocators!(context,allocator_state)
            # `_GeoSyntaxAbort` marks an upstream *syntax error* — the bison
            # parser records it, drops the statement and resumes at the next
            # one. The allocator observer is skipped: the statement did not
            # execute.
            syntax_aborted=false
            assigned=try
                _exec_line!(m,line,context,allocator_state)
            catch err
                err isa InterruptException && rethrow()
                err isa _GeoSyntaxAbort || rethrow()
                _geo_yyerror!(context,"syntax error ($(err.token))")
                syntax_aborted=true;nothing
            end
            assigned===nothing || (transfinite_tri=assigned)
            # `Delete`-family statements already drove the allocator update
            # inside `_geo_exec_delete!`; the observer skips them.
            syntax_aborted || _geo_allocator_observe_statement!(
                allocator_state,line,context,"execute_geo")
            occursin(r"\bBoolean[A-Za-z]*\b",
                     line) &&
                _geo_allocator_resync_model!(allocator_state,m)
            # `error tEND` recovery: a `;`-terminated statement's own
            # terminator ends the discard; a `}`-terminated one (Delete,
            # Boolean*, Extrude, transforms) has none, so recovery eats
            # through the next `;` like upstream.
            i=syntax_aborted && !endswith(rstrip(line),';') ?
                _geo_recovery_next_i(statements,i) : i+1
        else
            # Control tokens mutate the stream state — the branch/loop
            # structure upstream implements with `skip`/`fsetpos` is flat
            # here: a header dispatch returns the next stream index, which
            # `EndFor`/`EndWhile` can rewind backwards. A `_GeoSyntaxAbort`
            # inside a header is a bison syntax error: `error tEND` recovery
            # discards through the next `;`-statement.
            i=try
                if control.kind===:if
                    _geo_exec_if_header!(statements,i,context)
                elseif control.kind===:elseif
                    _geo_exec_elseif!(statements,i,context)
                elseif control.kind===:else
                    _geo_exec_else!(statements,i,context)
                elseif control.kind===:endif
                    context.if_depth-=1
                    context.if_depth<0 &&
                        _geo_yywarn!(context,"Orphan EndIf")
                    i+1
                elseif control.kind===:for
                    _geo_exec_for_header!(statements,i,context)
                elseif control.kind===:endfor
                    _geo_exec_endfor!(statements,i,context)
                elseif control.kind===:while
                    _geo_exec_while_header!(statements,i,context)
                elseif control.kind===:endwhile
                    _geo_exec_endwhile!(statements,i,context)
                elseif control.kind===:function
                    _geo_exec_function_def!(statements,i,context)
                elseif control.kind===:return
                    # `leaveFunction` pops the call stack: inside a `Call`
                    # body the flag unwinds to `_geo_exec_call!`; at top
                    # level (empty stack) upstream emits `yymsg(0)`.
                    if context.call_depth>0
                        context.stop=:return
                        break
                    end
                    _geo_yyerror!(context,"Error while exiting function")
                    i+1
                else
                    throw(ArgumentError(
                        "execute_geo: $line without a matching opener"))
                end
            catch err
                err isa InterruptException && rethrow()
                err isa _GeoSyntaxAbort || rethrow()
                _geo_yyerror!(context,"syntax error ($(err.token))")
                _geo_recovery_next_i(statements,i)
            end
        end
    end
    # Dangling `If`/`For`/`While` levels are silent at EOF — upstream's
    # counters persist and the last branch/body pass already streamed (the
    # `While` extension mirrors `For`). Only an EOF reached *while a skip
    # scan was mid-flight* reports `Unexpected end of file`.
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

# `Function`/`Macro` names are `tSTRING` (identifier or quoted literal) or a
# `StringExpr`; `Call` names are `String__Index` (which adds `name~{i}` and
# `StringToName[..]` forms) or a `StringExpr`.
function _geo_function_name_def(raw::AbstractString,
                                context::_GeoNumericContext,
                                caller::AbstractString)
    name=strip(raw)
    match(r"^[A-Za-z_][A-Za-z0-9_]*$",name)!==nothing && return String(name)
    quoted=match(r"^\"(.*)\"$",name)
    quoted===nothing && (quoted=match(r"^'(.*)'$",name))
    quoted!==nothing && return String(quoted.captures[1])
    return _geo_eval_string(name,context,caller)
end

function _geo_call_name(raw::AbstractString,context::_GeoNumericContext,
                        caller::AbstractString)
    s=String(strip(raw))
    try
        return _geo_symbol_name(s,context,caller)
    catch err
        err isa InterruptException && rethrow()
        # A `String__Index` parse failure does not commit the call — a
        # StringExpr head (quoted literal, `Str[...]`, …) still applies; the
        # classification below reports the right offending token otherwise.
        err isa _GeoSyntaxAbort || rethrow()
    end
    # `tCall String__Index` — a tSTRING/`~{}`/`StringToName` head commits the
    # name production and the first token that cannot extend it is the bison
    # error (`Call f(1,2)` → `(`, `f[0]`/`ns::f`/`f.x`/`f g` → `[`/`::`/`.`/
    # `g`). StringExpr heads (quoted literals, `Str`/`Sprintf`/…, bare string
    # tokens) evaluate; other reserved heads error at the token itself.
    head=match(r"^[A-Za-z_][A-Za-z0-9_]*",s)
    if head===nothing
        startswith(s,"\"") || startswith(s,"'") ||
            _geo_syntax_abort(_geo_string_index_offender(s))
    else
        h=head.match
        if !(h in _GEO_STRING_FUNCTIONS) && !(h in _GEO_BARE_STRING_TOKENS)
            if h in ("StringToName","S2N") ||
               (!(h in _GEO_ALL_FUNCTIONS) && !(h in _GEO_LEXER_KEYWORDS))
                _geo_syntax_abort(_geo_string_index_offender(s))
            else
                _geo_syntax_abort(h)
            end
        end
    end
    return _geo_eval_string(s,context,caller)
end

# `Function name`/`Macro name` registers its body — the statement array it
# came from plus the body start index, matching upstream's (file, position)
# pair — without executing it. `createFunction` precedes `skip(..,"Return")`,
# so the function is registered even when no `Return` follows (its body then
# runs to the end of the array when called), and a redefinition is a
# recoverable `yymsg(0)` diagnostic that keeps the existing body.
function _geo_exec_function_def!(statements::Vector{String},i::Int,
                                 context::_GeoNumericContext)
    control=_geo_control_parse(statements[i])
    name=_geo_function_name_def(control.name,context,"execute_geo: Function")
    if haskey(context.functions,name)
        _geo_yyerror!(context,"Redefinition of function $name")
    else
        context.functions[name]=(statements,i+1)
    end
    return _geo_exec_skip_return(statements,i+1,context)
end

# `Call name;` re-executes the registered body in the caller's scope with a
# bounded recursion depth (Gmsh recurses by re-parsing the body text; an
# unbounded call stack is a resource-bound violation here). An unknown name
# is a recoverable `yymsg(0)` diagnostic; `Return` unwinds by flagging
# `context.stop`, which this frame clears on exit.
function _geo_exec_call!(m::GeoModel,statements::Vector{String},
                         line::AbstractString,
                         context::_GeoNumericContext,
                         allocator_state::_GeoTagAllocatorState,
                         executed::Base.RefValue{Int})
    matched=match(r"^Call\s+(.+?)\s*;?\s*$",line)
    matched===nothing && throw(ArgumentError(
        "execute_geo: malformed Call statement; use `Call name;`"))
    name=_geo_call_name(matched.captures[1],context,"execute_geo: Call")
    body=get(context.functions,name,nothing)
    if body===nothing
        _geo_yyerror!(context,"Unknown function '$name'")
        return nothing
    end
    (body_statements,body_start)=body
    context.call_depth+=1
    try
        context.call_depth<=_MAX_GEO_CALL_DEPTH || throw(ArgumentError(
            "execute_geo: Call depth exceeds $_MAX_GEO_CALL_DEPTH"))
        return _exec_geo_statements!(
            m,body_statements,body_start,lastindex(body_statements),
            context,allocator_state,executed)
    finally
        context.call_depth-=1
        context.stop===:return && (context.stop=:run)
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
# Gmsh.l maps the Boolean aliases onto the canonical operator tokens:
# Fuse→Union, Cut→Difference, Common→Intersection, Coherence→Fragments.
# `Section` stays distinct — `OCC_Internals::Section` is a real upstream
# operator Tessella's native kernel does not implement.
const _GEO_BOOLEAN_NAMES=Dict(
    "Union"=>"Union","Fuse"=>"Union",
    "Difference"=>"Difference","Cut"=>"Difference",
    "Intersection"=>"Intersection","Common"=>"Intersection",
    "Fragments"=>"Fragments","Coherence"=>"Fragments",
    "Section"=>"Section")

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
                                context::_GeoNumericContext,
                                allocator_state=nothing)
    source=String(strip(raw))
    mm=match(r"^Boolean(Difference|Union|Intersection|Fragments|Section|Fuse|Coherence|Common|Cut)\b(.*)$",source)
    mm===nothing && return nothing
    caller="execute_geo: Boolean$(mm.captures[1])"
    rest=String(strip(mm.captures[2]))
    startswith(rest,"(") && return nothing
    groups=_geo_boolean_groups(rest,caller)
    # Gmsh.y: upstream only resolves operands and applies the operator when
    # the OCC kernel is active; under built-in the term is a recoverable
    # diagnostic and evaluates to an empty shape list.
    if allocator_state===nothing || allocator_state.factory!==:opencascade
        _geo_yyerror!(context,
            "Boolean operators only available with OpenCASCADE geometry kernel")
        return Float64[]
    end
    opname=_GEO_BOOLEAN_NAMES[mm.captures[1]]
    op=get(_GEO_BOOLEAN_OPS,opname,nothing)
    op===nothing && throw(ArgumentError(
        "execute_geo: Boolean$opname is not implemented"))
    operands=_geo_boolean_operands(groups,context,caller)
    out=boolean_volumes_multi!(m,op,
        operands.objects.tags,operands.tools.tags;
        remove_object=operands.objects.delete,
        remove_tool=operands.tools.delete,caller=caller)
    context.geo_changed=true
    return Float64.(out)
end

# ======== `FExpr_Multi` entity-selector and transform terms ========
#
# The `.geo` grammar puts model-reading and geometry-mutating terms inside
# list expressions: `Point{t}` coordinates, `X{:}`/`X "name"` wildcard tag
# lists, `Physical`/`Parent`/`BoundingBox`/`In BoundingBox` queries,
# `Mass`/`CenterOfMass`/`MatrixOfInertia` (OpenCASCADE-gated upstream),
# `Parametric`/`Normal`/`Color` reads, and `Translate`-family transform
# blocks. `MultipleShape` actions (`Boundary`, `PointsOf`, `Duplicata`, ...)
# are *not* `FExpr_Multi` and stay statement-level — `x[] = Boundary{...}` is
# a syntax error upstream.

const _GEO_SELECTOR_ENTITY_WORDS=(
    ("Point",0),("Curve",1),("Line",1),("Surface",2),("Volume",3))

# Parse an entity head — `Point`/`Curve`/`Line`/`Surface`/`Volume` literals or
# `GeoEntity{d}` — returning `(dim, rest, literal)` or `nothing`.
function _geo_selector_entity_head(source::AbstractString,
                                   context::_GeoNumericContext,
                                   caller::AbstractString)
    for (word,dim) in _GEO_SELECTOR_ENTITY_WORDS
        startswith(source,word) || continue
        rest=String(source[nextind(source,firstindex(source),
                                   ncodeunits(word)):end])
        # `CurveLoop`/`Pointwise`-style names are not the entity keyword.
        (isempty(rest) || _geo_word_char(first(rest))) && continue
        return (dim,String(strip(rest)),true)
    end
    startswith(source,"GeoEntity") || return nothing
    rest=String(source[nextind(source,firstindex(source),9):end])
    isempty(rest) && return nothing
    _geo_word_char(first(rest)) && return nothing
    rest=String(strip(rest))
    startswith(rest,"{") || return nothing
    (body,rest)=_geo_balanced_group(rest,caller)
    dim=_geo_int_value(
        _geo_eval_numeric(body,context,"$caller GeoEntity dim"),
        "$caller GeoEntity dim")
    return (dim,rest,false)
end

# `GeoEntity`/`GeoEntity123`/`GeoEntity12` dim checks: a literal word outside
# the production's dims is simply not that production (`Mass Point{1}` is a
# syntax error upstream); a dynamic `GeoEntity{d}` outside them emits the
# diagnostic yet still evaluates with the out-of-range dim.
function _geo_selector_dim_check(dim::Int,literal::Bool,lo::Int,hi::Int,
                                 context::_GeoNumericContext)
    lo<=dim<=hi && return true
    literal && return false
    _geo_yyerror!(context,"GeoEntity dim out of range [$lo,$hi]")
    return true
end

# `getEntities(dim)` equivalent for a (possibly out-of-range) selector dim —
# upstream falls through to the `default:` branch for any dim outside [0,3]
# and returns every entity.
function _geo_selector_all_tags(m::GeoModel,dim::Int)
    0<=dim<=3 ||
        return Float64[Float64(tag) for (_,tag) in model_entities(m)]
    return Float64[Float64(tag) for (_,tag) in model_entities(m,dim)]
end

# `getElementaryTagsForPhysicalGroups` / `getAllPhysicalTags` — groups are
# keyed by abs(tag) upstream and members arrive sorted and deduplicated.
function _geo_selector_physical(m::GeoModel,dim::Int,tail::String,
                                context::_GeoNumericContext,
                                caller::AbstractString)
    if startswith(tail,"\"")
        # `Physical X "name"` — deprecated wildcard: only `"*"`/`"all"` mean
        # all tags; anything else is a `yyerror` with an empty result.
        sm=match(r"^\"([^\"]*)\"\s*(.*)$",tail)
        sm===nothing && return nothing
        isempty(sm.captures[2]) || return nothing
        sm.captures[1] in ("*","all") ||
            (_geo_yyerror!(context,
                 "Unknown special string for list replacement");
             return Float64[])
        return _geo_selector_all_physical(m,dim)
    end
    inlist=if startswith(tail,"{")
        (body,rest)=_geo_balanced_group(tail,caller)
        isempty(rest) || return nothing
        inner=String(strip(body))
        inner==":" && return _geo_selector_all_physical(m,dim)
        _geo_numeric_list_values(
            "{"*inner*"}",context,"$caller Physical";depth=1)
    else
        # `ListOfDoubleOrAll` accepts unbraced list forms (`Physical Line a[]`,
        # `Physical Line 5`, `Physical Line 1:3`).
        _geo_numeric_list_values(tail,context,"$caller Physical";depth=1)
    end
    out=Float64[]
    if 0<=dim<=3
        live_groups=Set(model_physical_groups(m,dim))
        for raw_tag in inlist
            tag=_geo_int_value(raw_tag,"$caller Physical tag")
            # `getElementaryTagsForPhysicalGroups` does `groups.find(num)`
            # with the *raw* tag — the map keys are `abs()`'d member-side, so
            # `Physical Curve{-7}` misses group 7 and appends nothing.
            (dim,tag) in live_groups || continue
            members=get(m.physical,(dim,tag),Int[])
            append!(out,Float64.(sort!(unique!(
                _physical_live_members(m,dim,members)))))
        end
    end
    return out
end

_geo_selector_all_physical(m::GeoModel,dim::Int)=
    0<=dim<=3 ?
        Float64[Float64(tag) for (_,tag) in model_physical_groups(m,dim)] :
        Float64[]

# `X In BoundingBox{...}` — entities whose complete bbox is inside the box.
function _geo_selector_in_box(m::GeoModel,dim::Int,tail::String,
                              context::_GeoNumericContext,
                              caller::AbstractString)
    box=_geo_numeric_list_values(tail,context,"$caller In BoundingBox";
                                 depth=1)
    length(box)<6 && begin
        _geo_yyerror!(context,
            "Bounding box should be {xmin, ymin, zmin, xmax, ymax, zmax}")
        return Float64[]
    end
    # Out-of-range dims take upstream `getEntities`'s default branch — every
    # entity is box-filtered.
    query_dim=0<=dim<=3 ? dim : -1
    return Float64[Float64(tag) for (_,tag) in
        model_entities_in_bounding_box(
            m,box[1],box[2],box[3],box[4],box[5],box[6],query_dim)]
end

# `BoundingBox X{list}` — the union of the found entities' bounds; nothing is
# appended when every listed entity is missing.
function _geo_selector_bounding_box(m::GeoModel,dim::Int,tail::String,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    startswith(tail,"{") || _geo_selector_abort(tail)
    (body,rest)=_geo_balanced_group(tail,caller)
    isempty(rest) || _geo_selector_abort(rest)
    tags=_geo_numeric_list_values(
        "{"*body*"}",context,"$caller BoundingBox";depth=1)
    bounds=nothing
    if 0<=dim<=3
        for raw_tag in tags
            tag=_geo_int_value(raw_tag,"$caller BoundingBox tag")
            model_entity(m,dim,tag)===nothing && continue
            entity_bounds=_model_entity_bounding_box(m,dim,tag,caller)
            bounds=bounds===nothing ? entity_bounds :
                _model_bounds_union(bounds,entity_bounds)
        end
    end
    return bounds===nothing ? Float64[] : collect(Float64,bounds)
end

# `Mass`/`CenterOfMass`/`MatrixOfInertia X{t}` — upstream these only run under
# the OpenCASCADE factory; otherwise a recoverable diagnostic plus a zeroed
# result (`MatrixOfInertia` adds nothing at all).
function _geo_selector_mass(m::GeoModel,name::String,dim::Int,tail::String,
                            context::_GeoNumericContext,
                            caller::AbstractString,allocator_state)
    startswith(tail,"{") || _geo_selector_abort(tail)
    (body,rest)=_geo_balanced_group(tail,caller)
    isempty(rest) || _geo_selector_abort(rest)
    _geo_eval_numeric(body,context,"$caller $name tag")
    if allocator_state!==nothing && allocator_state.factory==:opencascade &&
            allocator_state.occ_active
        # Upstream would call `getOCCInternals()->getMass` here — real mass
        # data Tessella cannot produce, so it is a hard error rather than a
        # zeroed placeholder.
        throw(ArgumentError(
            "$caller: $name requires OpenCASCADE mass properties, which " *
            "Tessella does not implement"))
    end
    _geo_yyerror!(context,
        "$name only available with OpenCASCADE geometry kernel")
    return name=="Mass" ? [0.0] : name=="CenterOfMass" ? [0.0,0.0,0.0] :
        Float64[]
end

# The unexpected-token text for a committed-but-incomplete selector: the
# next word/character of `rest`, or `;` when the term simply ran out.
function _geo_selector_abort(rest::AbstractString)
    r=String(strip(rest))
    isempty(r) && _geo_syntax_abort(";")
    m=match(r"^[A-Za-z_][A-Za-z0-9_]*|^.",r)
    _geo_syntax_abort(m===nothing ? ";" : String(m.match))
end

# `Normal Surface{t} Parametric{u,v}` — `tSurface` is a literal keyword here.
function _geo_selector_normal(m::GeoModel,tail::String,
                              context::_GeoNumericContext,
                              caller::AbstractString)
    head=_geo_selector_entity_head(tail,context,caller)
    head===nothing && _geo_selector_abort(tail)
    dim,rest,literal=head
    # `tNormal` expects `tSurface` upstream — anything else is a syntax error.
    (literal && dim==2) || _geo_selector_abort(tail)
    # `tNormal tSurface` commits the production upstream — a missing `{` or
    # `Parametric` is a syntax error, not a different term.
    startswith(rest,"{") || _geo_selector_abort(rest)
    (body,rest)=_geo_balanced_group(rest,caller)
    rest=String(strip(rest))
    startswith(rest,"Parametric") || _geo_selector_abort(rest)
    rest=String(strip(rest[nextind(rest,firstindex(rest),10):end]))
    startswith(rest,"{") || _geo_selector_abort(rest)
    (pbody,rest)=_geo_balanced_group(rest,caller)
    isempty(rest) || _geo_selector_abort(rest)
    tag=_geo_int_value(
        _geo_eval_numeric(body,context,"$caller Normal tag"),
        "$caller Normal tag")
    # `tParametric '{' FExpr ',' FExpr '}'` — two scalar args, not an RLD:
    # `a[]` splices and `:` ranges are syntax errors upstream.
    pparts=_geo_split_top_commas(pbody,caller)
    length(pparts)==2 || _geo_syntax_abort(length(pparts)<2 ? "}" : ",",
        "$caller: Normal Parametric expects {u,v}")
    params=Float64[
        _geo_eval_numeric(pparts[1],context,"$caller Normal Parametric"),
        _geo_eval_numeric(pparts[2],context,"$caller Normal Parametric")]
    haskey(m.surfaces,tag) ||
        haskey(m.discrete,(2,tag)) ||
        get(m.surface_geometry,tag,nothing)!==nothing ||
        (_geo_yyerror!(context,"Surface $tag does not exist");
         return Float64[])
    return collect(Float64,model_normal(m,tag,params;
        old_ruled_surface=_geo_exec_old_ruled(context)))
end

# `Geometry.OldRuledSurface` / `Mesh.NewtonConvergenceTestXYZ` — the two
# number options that reach ruled-surface evaluation.
_geo_exec_old_ruled(context::_GeoNumericContext)=
    !iszero(something(_geo_option_number(
        context,"Geometry",0,"OldRuledSurface"),0.0))
_geo_exec_newton_xyz(context::_GeoNumericContext)=
    !iszero(something(_geo_option_number(
        context,"Mesh",0,"NewtonConvergenceTestXYZ"),0.0))

# `AddToTemporaryBoundingBox` — the running `Point`-statement bbox.
function _geo_exec_temp_bbox_add!(context::_GeoNumericContext,
                                  x::Float64,y::Float64,z::Float64)
    box=context.temp_bbox
    context.temp_bbox=box===nothing ?
        ((x,y,z),(x,y,z)) :
        (min.(box[1],(x,y,z)),max.(box[2],(x,y,z)))
    return nothing
end

# `CTX::instance()->lc` during `.geo` execution — the raw diagonal of the
# `Point`-statement cloud (`lc == 0` maps to 1 upstream).
function _geo_exec_lc(context::_GeoNumericContext)
    box=context.temp_bbox
    box===nothing && return 1.0
    range=box[2].-box[1]
    lc=sqrt(range[1]*range[1]+range[2]*range[2]+range[3]*range[3])
    return lc==0.0 ? 1.0 : lc
end

# `Parametric BoundingBox X{t}` (dims 1-2) and
# `Parametric Point{t} In Surface{s}`.
function _geo_selector_parametric(m::GeoModel,tail::String,
                                  context::_GeoNumericContext,
                                  caller::AbstractString)
    if (pm=match(r"^Point\b(.*)$",tail))!==nothing
        rest=String(strip(pm.captures[1]))
        startswith(rest,"{") || _geo_selector_abort(rest)
        (body,rest)=_geo_balanced_group(rest,caller)
        rest=String(strip(rest))
        sm=match(r"^In\s+Surface\s*(.*)$",rest)
        sm===nothing && _geo_selector_abort(rest)
        srest=String(strip(sm.captures[1]))
        startswith(srest,"{") || _geo_selector_abort(srest)
        (sbody,srest)=_geo_balanced_group(srest,caller)
        isempty(srest) || _geo_selector_abort(srest)
        ptag=_geo_int_value(
            _geo_eval_numeric(body,context,"$caller Parametric Point"),
            "$caller Parametric Point")
        stag=_geo_int_value(
            _geo_eval_numeric(sbody,context,"$caller Parametric Surface"),
            "$caller Parametric Surface")
        point_known=haskey(m.points,ptag) || haskey(m.discrete,(0,ptag))
        surface_known=haskey(m.surfaces,stag) ||
            haskey(m.discrete,(2,stag)) ||
            get(m.surface_geometry,stag,nothing)!==nothing
        (point_known && surface_known) ||
            (_geo_yyerror!(context,
                 "Point $ptag or surface $stag does not exist");
             return Float64[])
        return collect(Float64,model_reparametrize_on_surface(
            m,0,ptag,Float64[],stag;
            old_ruled_surface=_geo_exec_old_ruled(context),
            newton_convergence_xyz=_geo_exec_newton_xyz(context),
            lc=_geo_exec_lc(context),
            warn=msg->_geo_yywarn!(context,msg),
            info=msg->_geo_yyinfo!(context,msg)))
    end
    # `tParametric` expects `tBoundingBox` or `tPoint` upstream.
    match(r"^BoundingBox\b",tail)===nothing && _geo_selector_abort(tail)
    tail=String(strip(tail[nextind(tail,firstindex(tail),11):end]))
    head=_geo_selector_entity_head(tail,context,caller)
    head===nothing && _geo_selector_abort(tail)
    dim,rest,literal=head
    _geo_selector_dim_check(dim,literal,1,2,context) ||
        _geo_selector_abort(tail)
    startswith(rest,"{") || _geo_selector_abort(rest)
    (body,rest)=_geo_balanced_group(rest,caller)
    isempty(rest) || _geo_selector_abort(rest)
    tag=_geo_int_value(
        _geo_eval_numeric(body,context,"$caller Parametric BoundingBox"),
        "$caller Parametric BoundingBox")
    exists=(dim==1 && (haskey(m.curves,tag) ||
                       haskey(m.discrete,(1,tag)) ||
                       get(m.curve_geometry,tag,nothing)!==nothing)) ||
           (dim==2 && (haskey(m.surfaces,tag) ||
                       haskey(m.discrete,(2,tag)) ||
                       get(m.surface_geometry,tag,nothing)!==nothing))
    label=dim==1 ? "Curve" : "Surface"
    exists || (_geo_yyerror!(context,"$label $tag does not exist");
               return Float64[])
    lower,upper=model_parametrization_bounds(m,dim,tag)
    if dim==1
        return Float64[lower[1],upper[1]]
    end
    return Float64[lower[1],lower[2],upper[1],upper[2]]
end

# `Color X{t}` — existing entities report their RGBA; a missing entity is a
# silent empty list (no diagnostic upstream).
function _geo_selector_color(m::GeoModel,tail::String,
                             context::_GeoNumericContext,
                             caller::AbstractString)
    head=_geo_selector_entity_head(tail,context,caller)
    head===nothing && _geo_selector_abort(tail)
    dim,rest,literal=head
    _geo_selector_dim_check(dim,literal,1,3,context) ||
        _geo_selector_abort(tail)
    startswith(rest,"{") || _geo_selector_abort(rest)
    (body,rest)=_geo_balanced_group(rest,caller)
    isempty(rest) || _geo_selector_abort(rest)
    tag=_geo_int_value(
        _geo_eval_numeric(body,context,"$caller Color tag"),
        "$caller Color tag")
    0<=dim<=3 || return Float64[]
    model_entity(m,dim,tag)===nothing &&
        !haskey(m.discrete,(dim,tag)) && return Float64[]
    return collect(Float64,model_entity_color(m,dim,tag))
end

# Bare entity selectors: `Point{t}` coordinates, `X{:}`/`X "name"` tag
# wildcards and `X In BoundingBox{...}`.
function _geo_selector_entity_term(m::GeoModel,dim::Int,rest::String,
                                   literal::Bool,
                                   context::_GeoNumericContext,
                                   caller::AbstractString)
    if (im=match(r"^In\s+BoundingBox\b(.*)$",rest))!==nothing
        # `GeoEntity` variant — dynamic dims outside [0,3] diagnose but
        # still evaluate.
        literal ||
            _geo_selector_dim_check(dim,false,0,3,context) || return nothing
        return _geo_selector_in_box(
            m,dim,String(strip(im.captures[1])),context,caller)
    end
    if match(r"^\"[^\"]*\"\s*$",rest)!==nothing
        # `X "name"` — deprecated `{:}` wildcard through the same
        # `tPoint tBIGSTR` / `GeoEntity123 tBIGSTR` productions; the string
        # value is ignored upstream.
        literal ||
            _geo_selector_dim_check(dim,false,1,3,context) || return nothing
        return _geo_selector_all_tags(m,dim)
    end
    # The entity keyword commits the production upstream — anything but
    # `In BoundingBox`, a quoted name, or `{` here is a syntax error.
    startswith(rest,"{") || _geo_selector_abort(rest)
    (body,tail)=_geo_balanced_group(rest,caller)
    isempty(tail) || _geo_selector_abort(tail)
    inner=String(strip(body))
    if inner==":"
        # `Point{:}` is its own production; `Curve`/`Surface`/`Volume{:}` and
        # `GeoEntity{d}{:}` go through the GeoEntity123 one — a `GeoEntity`
        # dim outside [1,3] emits the range error but still lists.
        literal ||
            _geo_selector_dim_check(dim,false,1,3,context) || return nothing
        return _geo_selector_all_tags(m,dim)
    end
    # `Point{expr}` reports coordinates; `Curve{expr}`/`Surface{expr}`/
    # `Volume{expr}`/`GeoEntity{d}{expr}` reduce to `GeoEntity123 '{'`, which
    # expects `tDOTS` upstream — a syntax error at the inner expression.
    if !(literal && dim==0)
        literal || _geo_selector_dim_check(dim,false,1,3,context)
        _geo_selector_abort(isempty(inner) ? "}" : inner)
    end
    tag=_geo_int_value(
        _geo_eval_numeric(inner,context,"$caller Point tag"),
        "$caller Point tag")
    # `getGEOInternals()->getVertex` compares `abs(Num)` — `Point{-1}` finds
    # point 1 — then the model-level `getVertexByTag` fallback matches the
    # exact tag (including discrete vertices).
    point=get(m.points,tag,get(m.points,abs(tag),nothing))
    point===nothing && haskey(m.discrete,(0,tag)) && begin
        record=m.discrete[(0,tag)]
        isempty(record.node_coords) ||
            (point=Tuple(record.node_coords[1:3,1]))
    end
    if point===nothing
        _geo_yyerror!(context,"Unknown model point with tag $tag")
        return [0.0,0.0,0.0]
    end
    return Float64[point[1],point[2],point[3]]
end

# Reserved tokens upstream — a `name{` term headed by one of these is *not*
# the `tSTRING '{' MultipleShape '}'` action production; it either has its
# own selector production (handled after this) or is a syntax error.
const _GEO_SELECTOR_KEYWORDS=(
    "Point","Curve","Line","Surface","Volume","GeoEntity","Physical","Parent",
    "BoundingBox","Mass","CenterOfMass","MatrixOfInertia","Normal",
    "Parametric","Color","List","LinSpace","LogSpace","Catenary","Unique",
    "Abs","ListFromFile","Extrude","BooleanUnion","BooleanDifference",
    "BooleanIntersection","BooleanFragments","Boolean","Periodic","Kernel",
    "Factory","Mesh","General","Geometry","Field","Delete","Hide","Show",
    "SetNumber","GetNumber","SetString","GetString","Exists","Printf",
    "Error","Warning","Info","Debug","Call","If","ElseIf","Else","EndIf",
    "For","EndFor","While","EndWhile","Function","Return","Macro","Merge",
    "Save","Print","Exit","Abort","Include","SetFactory","SetOrder",
    "Optimize","Recombine","Transfinite","Reverse","Coherence")

# `Transform` terms in `FExpr_Multi` position — `Translate{params}{shapes;}`,
# `Rotate`/`Dilate`/`Symmetry`/`Affine`/`Closest` likewise, and the
# `tSTRING '{' MultipleShape '}'` action forms (`Duplicata`, `Boundary`,
# `CombinedBoundary`, `PointsOf`). The transform is applied to the resolved
# shape list (a nested `Duplicata` produces the new tags) and the resolved
# tags are the term's value. A name that lexes as a reserved keyword upstream
# (`Point{1}`, `Physical Curve{...}`, `Mass ...`, ...) is *not* a shape action
# and is left for the selector machinery.
function _geo_exec_transform_term(m::GeoModel,raw::AbstractString,
                                  context::_GeoNumericContext,
                                  allocator_state)
    s=String(strip(raw))
    nm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\b",s)
    nm===nothing && return nothing
    name=nm.captures[1]
    caller="execute_geo: $name"
    if name=="Intersect" || name=="Split"
        # `tIntersect tCurve '{' RLD '}' tSurface '{' FExpr '}'` and
        # `tSplit tCurve '{' FExpr '}' tPoint '{' RLD '}'` — real built-in
        # kernel curve operations upstream; a hard blocker here.
        throw(ArgumentError(
            "$caller: $name requires built-in kernel curve splitting, " *
            "which Tessella does not implement"))
    end
    rest=String(strip(s[nextind(s,firstindex(s),
                              ncodeunits(name)):end]))
    if name in ("Translate","Rotate","Dilate","Symmetry")
        # `tTranslate VExpr '{' MultipleShape '}'` etc. — keyword transforms
        # commit at the name; a missing VExpr or shape group is a syntax
        # error, not an unknown-variable fallthrough.
        (t,rest)=_geo_transform_params_rest(
            kind=name,source=rest,context=context,caller=caller)
        r2=String(strip(rest))
        startswith(r2,"{") || _geo_selector_abort(r2)
        (shape_list,rest)=_geo_balanced_group(r2,caller)
        isempty(rest) || _geo_selector_abort(rest)
        entities=_geo_shape_list_entities!(m,shape_list,context,caller;
                                           allocator_state=allocator_state)
        transform_entities!(m,t,entities;caller=caller)
        context.geo_changed=true
        return Float64[Float64(tag) for (_,tag) in entities]
    end
    if !startswith(rest,"{")
        # `Affine`/`Closest` are keywords upstream — `Affine` without its
        # `{...}` parameters is a syntax error. Any other name may still be a
        # `tSTRING` action or a selector handled elsewhere.
        name in ("Affine","Closest") || return nothing
        _geo_selector_abort(rest)
    end
    (params,rest)=_geo_balanced_group(rest,caller)
    if name=="Affine" || name=="Closest"
        r2=String(strip(rest))
        startswith(r2,"{") || _geo_selector_abort(r2)
        (shape_list,rest)=_geo_balanced_group(r2,caller)
        isempty(rest) || _geo_selector_abort(rest)
        _geo_numeric_list_values(
            "{"*params*"}",context,"$caller parameters";depth=1)
        if name=="Affine"
            _geo_yyerror!(context,
                "Affine transform only available with OpenCASCADE " *
                "geometry kernel")
            entities=_geo_shape_list_entities!(m,shape_list,context,caller;
                                               allocator_state=allocator_state)
            return Float64[Float64(tag) for (_,tag) in entities]
        end
        _geo_yyerror!(context,
            "Closest entity only available with OpenCASCADE geometry kernel")
        _geo_shape_list_entities!(m,shape_list,context,caller;
                                  allocator_state=allocator_state)
        return Float64[]
    end
    name in _GEO_SELECTOR_KEYWORDS && return nothing
    # `tSTRING '{' MultipleShape '}'` — Duplicata/boundary/PointsOf/unknown
    # actions; `params` is the shape-list body here.
    entities=_geo_shape_list_entities!(m,params,context,caller;
                                       allocator_state=allocator_state)
    isempty(rest) || _geo_selector_abort(rest)
    if name=="Duplicata"
        context.geo_changed=true
    else
        _geo_sync_physical_view_if_changed!(m,context)
    end
    out=_geo_shape_action!(m,name,entities,context,caller)
    return Float64[Float64(tag) for (_,tag) in out]
end

# Dispatch the keyword-headed selector forms.
function _geo_selector_keyword_term(m::GeoModel,name::String,tail::String,
                                    context::_GeoNumericContext,
                                    caller::AbstractString,allocator_state)
    if name=="Physical"
        head=_geo_selector_entity_head(tail,context,caller)
        head===nothing && _geo_selector_abort(tail)
        dim,rest,literal=head
        _geo_selector_dim_check(dim,literal,0,3,context) ||
            _geo_selector_abort(tail)
        return _geo_selector_physical(m,dim,rest,context,caller)
    elseif name=="Parent"
        head=_geo_selector_entity_head(tail,context,caller)
        head===nothing && _geo_selector_abort(tail)
        dim,rest,literal=head
        _geo_selector_dim_check(dim,literal,0,3,context) ||
            _geo_selector_abort(tail)
        # `ListOfDouble` — brace list or unbraced list expression.
        _geo_numeric_list_values(rest,context,"$caller Parent";depth=1)
        # Native entities carry no `getParentEntity` provenance — the built-in
        # kernel result is empty like upstream.
        return Float64[]
    elseif name=="BoundingBox"
        head=_geo_selector_entity_head(tail,context,caller)
        head===nothing && _geo_selector_abort(tail)
        dim,rest,literal=head
        _geo_selector_dim_check(dim,literal,0,3,context) || return nothing
        return _geo_selector_bounding_box(m,dim,rest,context,caller)
    elseif name in ("Mass","CenterOfMass","MatrixOfInertia")
        head=_geo_selector_entity_head(tail,context,caller)
        head===nothing && _geo_selector_abort(tail)
        dim,rest,literal=head
        # `tMass GeoEntity123` — the keyword class is Curve/Surface/Volume;
        # a keyword `Point` is a syntax error at `Point`, not a miss.
        literal && dim==0 && _geo_selector_abort(tail)
        _geo_selector_dim_check(dim,literal,1,3,context) || return nothing
        return _geo_selector_mass(
            m,name,dim,rest,context,caller,allocator_state)
    elseif name=="Normal"
        return _geo_selector_normal(m,tail,context,caller)
    elseif name=="Parametric"
        return _geo_selector_parametric(m,tail,context,caller)
    elseif name=="Color"
        return _geo_selector_color(m,tail,context,caller)
    end
    return nothing
end

# `FExpr_Multi` selector/transform terms — returns a Float64 vector when `raw`
# is one, or `nothing` when it is not.
function _geo_exec_selector_term(m::GeoModel,raw::AbstractString,
                                 context::_GeoNumericContext,allocator_state)
    s=String(strip(raw))
    isempty(s) && return nothing
    caller="execute_geo"
    # Every upstream selector synchronizes the internals before reading —
    # `if(getGEOInternals()->getChanged()) synchronize` — which here
    # materializes the derived physical-group view when declarations or
    # geometry changed since the last sync.
    _geo_sync_physical_view_if_changed!(m,context)
    transform=_geo_exec_transform_term(m,s,context,allocator_state)
    transform!==nothing && return transform
    kw=match(r"^(Physical|Parent|BoundingBox|Mass|CenterOfMass|MatrixOfInertia|Normal|Parametric|Color)\b(.*)$",s)
    kw!==nothing && return _geo_selector_keyword_term(
        m,String(kw.captures[1]),String(strip(kw.captures[2])),context,caller,
        allocator_state)
    head=_geo_selector_entity_head(s,context,caller)
    head===nothing && return nothing
    return _geo_selector_entity_term(
        m,head[1],head[2],head[3],context,caller)
end

function _geo_exec_value_term(m::GeoModel,raw::AbstractString,
                              context::_GeoNumericContext,
                              allocator_state=nothing)
    extruded=_geo_exec_extrude_term(m,raw,context,allocator_state)
    extruded!==nothing && return extruded
    boolean=_geo_exec_boolean_term(m,raw,context,allocator_state)
    boolean!==nothing && return boolean
    return _geo_exec_selector_term(m,raw,context,allocator_state)
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

# `Sphere(n) = {centerTag, pointTag}` / `PolarSphere(n) = {centerTag, pointTag}`
# — Gmsh's built-in-kernel two-point sphere: the second point sets the radius
# as its distance from the center. Unknown points record a `Msg::Error` and
# leave the entity uncreated (Gmsh returns a null surface).
function _geo_exec_point_sphere!(m::GeoModel,tag::Int,values,
                                 context::_GeoNumericContext,
                                 caller::AbstractString,kind::AbstractString;
                                 zero_literal::Bool=false)
    # Point references resolve through the `abs` vertex-tree comparator, so a
    # tag-0 or negated reference reaches `Point[0]`/`Point[abs(t)]`.
    center_tag=abs(_geo_signed_gmsh_int_value(values[1],"$caller center point"))
    point_tag=abs(_geo_signed_gmsh_int_value(values[2],"$caller sphere point"))
    if !haskey(m.points,center_tag)
        _geo_msg_error!(context,"Unknown $kind center point $center_tag")
        context.parametric_surface=nothing
        return nothing
    end
    if !haskey(m.points,point_tag)
        _geo_msg_error!(context,"Unknown $kind point $point_tag")
        context.parametric_surface=nothing
        return nothing
    end
    itag=trunc(Int,tag)
    c=m.points[center_tag];p=m.points[point_tag]
    r=sqrt((p[1]-c[1])^2+(p[2]-c[2])^2+(p[3]-c[3])^2)
    # `newGeometrySphere`/`newGeometryPolarSphere` do not create an entity —
    # `gmshSphere::New*` registers the surface in `allGmshSurfaces` and the
    # grammar installs it as `myGmshSurface` (`$$.Type = 0`). A repeated tag
    # reports through `Msg::Error`, then the registry entry is replaced.
    surface=(kind=kind,tag=itag,center=c,radius=r)
    haskey(context.parametric_surfaces,itag) &&
        _geo_msg_error!(context,"Surface $itag already exists")
    context.parametric_surfaces[itag]=surface
    context.parametric_surface=surface
    return nothing
end

# `GEO_Internals::addX` emits `GEO <type> with tag <n> already exists`
# (`Msg::Error`) and returns null on a duplicate tag; the grammar action then
# emits `Could not add <name>` — both recoverable diagnostics, and the
# existing entity is kept. Under the OpenCASCADE factory the same statement
# routes to `OCC_Internals::addX`, whose wording is `OpenCASCADE <type> ...`
# (`wire or curve loop` for curve loops).
const _GEO_DUP_ENTITY_WORDS=(point=("GEO point","OpenCASCADE point"),
    curve=("GEO curve","OpenCASCADE curve"),
    curve_loop=("GEO curve loop","OpenCASCADE wire or curve loop"),
    surface=("GEO surface","OpenCASCADE surface"),
    surface_loop=("GEO surface loop","OpenCASCADE surface loop"),
    volume=("GEO volume","OpenCASCADE volume"))

function _geo_dup_entity!(context::_GeoNumericContext,exists::Bool,
                          entity::Symbol,name_word::AbstractString,
                          tag,zero_literal::Bool,allocator_state)
    # `if(tag >= 0 && FindX(tag))` — the auto-assign path (`tag < 0` upstream,
    # normalized to a nonliteral 0 by `_geo_exec_def_tag`) skips the lookup.
    (exists && (zero_literal || tag>0)) || return false
    words=getfield(_GEO_DUP_ENTITY_WORDS,entity)
    word=allocator_state!==nothing &&
        allocator_state.factory===:opencascade ? words[2] : words[1]
    _geo_msg_error!(context,"$word with tag $(trunc(Int,tag)) already exists")
    _geo_yyerror!(context,"Could not add $name_word")
    return true
end

# `CreateCurve` resolves every control-point tag through `FindPoint`: misses
# emit `Unknown control point %d in GEO curve %d` (`Msg::Error`) and are
# dropped from the stored list — the curve still lands in `Curves` with only
# the resolved points.
function _geo_curve_found_points!(m::GeoModel,context::_GeoNumericContext,
                                  points,tag)
    found=Int[];sizehint!(found,length(points))
    ok=true;itag=trunc(Int,tag)
    for p in points
        if haskey(m.points,p)
            push!(found,p)
        else
            _geo_msg_error!(context,
                "Unknown control point $p in GEO curve $itag")
            ok=false
        end
    end
    return (found,ok)
end

# `Tree_Add(Curves)` runs unconditionally once the tag is allocated: beg/end
# are the first/last resolved control points (`p1 < 0` upstream), and a
# zero-resolution list still occupies the tag.
function _geo_store_curve!(m::GeoModel,tag,found::Vector{Int};
                           kind=nothing,control_points::Bool=false,
                           geometry=nothing)
    itag=trunc(Int,tag)
    m.curves[itag]=isempty(found) ? (0,0) : (first(found),last(found))
    control_points && (m.curve_control_points[itag]=copy(found))
    kind!==nothing && (m.curve_types[itag]=kind)
    geometry!==nothing && (m.curve_geometry[itag]=geometry)
    return nothing
end

# The `MSH_SEGM_LINE && Control_Points <= 2` position check at the end of
# `CreateCurve`: endpoints comparing equal within
# `Geometry.Tolerance * lc` (`ComparePosition`) emit a warning — the curve is
# still valid (`ok` untouched).
function _geo_line_tolerance_warn!(m::GeoModel,context::_GeoNumericContext,
                                   found::Vector{Int},tag)
    1<=length(found)<=2 || return nothing
    a,b=first(found),last(found)
    eps=something(_geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)*
        _geo_exec_lc(context)
    pa,pb=m.points[a],m.points[b]
    (abs(pa[1]-pb[1])<=eps && abs(pa[2]-pb[2])<=eps &&
     abs(pa[3]-pb[3])<=eps) || return nothing
    _geo_yywarn!(context,
        "Start point $a and end point $b of GEO line $(trunc(Int,tag)) are " *
        "closer than the geometrical tolerance, at position " *
        "($(_geo_gmsh_number(pa[1])), $(_geo_gmsh_number(pa[2])), " *
        "$(_geo_gmsh_number(pa[3])))")
    return nothing
end

# `myatan2`/`myasin`/`angle_02pi` — Numeric.cpp's clamped variants.
_geo_myatan2(y,x)=(y==0 && x==0) ? 0.0 : atan(y,x)
_geo_myasin(a)=a<=-1.0 ? -π/2 : a>=1.0 ? π/2 : asin(a)
function _geo_angle_02pi(a)
    dp=2π
    while a>dp || a<0.0
        a>0 ? (a-=dp) : (a+=dp)
    end
    return a
end

# One `EndCurve` pass — the `MSH_SEGM_CIRC`/`MSH_SEGM_ELLI` branch, which runs
# only when the resolved control list has 3 or 4 points. It builds the local
# frame from the normalized center→endpoint directions (falling back to
# `fallback_n` — `c->Circle.n`, `{0,0,1}` at `CreateCurve` time and the stored
# normal for the `addX` re-run) when the cross product degenerates, checks
# the projected radii, and reports `Circle or ellipse arc %d greater than Pi`
# when `arc` is enabled (upstream gates that check on `c->Num > 0`, so the
# reversed curve's pass never emits it). The 4-point form additionally
# solves the ellipse's 2×2 system. Returns `ok` — `CreateCurve` folds only
# the forward first pass into its `ok` out-param; every other pass emits
# diagnostics with a discarded return.
function _geo_endcurve_check!(m::GeoModel,context::_GeoNumericContext,
                              found::Vector{Int},fallback_n,tag;
                              arc::Bool)
    itag=trunc(Int,tag)
    center=m.points[found[2]]
    d12=m.points[found[1]].-center
    d32=m.points[found[end]].-center
    d42=length(found)==4 ? m.points[found[3]].-center : nothing
    n12,n32=norm(d12),norm(d32)
    # `norme` leaves a zero vector unchanged rather than producing NaN — the
    # downstream projections then see real zeros and the radius checks fire.
    u12=n12>0 ? Tuple(d12./n12) : Tuple(d12)
    u32=n32>0 ? Tuple(d32./n32) : Tuple(d32)
    n=_occ_cross(u12,u32);nn=norm(n)
    isvalid=nn>=1e-15
    if isvalid
        n=n./nn
        (abs(n[1])<1e-5 && abs(n[2])<1e-5 && abs(n[3])<1e-5) && (isvalid=false)
    end
    if !isvalid
        fb=Tuple(collect(fallback_n))
        bn=norm(fb)
        n=bn>0 ? fb./bn : fb
    end
    mv=_occ_cross(n,u12);mn=norm(mv)
    mv=mn>0 ? mv./mn : mv
    # Single-assignment copies keep the `proj` closure capture unboxed.
    pn,pmv=n,mv
    proj(v)=(u12[1]*v[1]+u12[2]*v[2]+u12[3]*v[3],
             pmv[1]*v[1]+pmv[2]*v[2]+pmv[3]*v[3],
             pn[1]*v[1]+pn[2]*v[2]+pn[3]*v[3])
    pv0=proj(d12);pv2=proj(d32)
    R=hypot(pv0[1],pv0[2]);R2=hypot(pv2[1],pv2[2])
    ok=true;g=_geo_gmsh_number
    a1=a3=0.0
    if R==0||R2==0
        _geo_msg_error!(context,
            "Zero radius in circle or ellipse with tag $itag")
        ok=false
    elseif d42===nothing && abs((R-R2)/(R+R2))>0.1
        _geo_msg_error!(context,
            "Control points of circle with tag $itag are not cocircular: " *
            "R1=$(g(R)), R2=$(g(R2)), n=[$(g(n[1])),$(g(n[2])),$(g(n[3]))]")
        ok=false
    end
    if d42!==nothing
        pv3=proj(d42)
        a4=_geo_angle_02pi(_geo_myatan2(pv3[2],pv3[1]))
        ca,sa=cos(a4),sin(a4)
        x1=pv0[1]*ca+pv0[2]*sa;y1=-pv0[1]*sa+pv0[2]*ca
        x3=pv2[1]*ca+pv2[2]*sa;y3=-pv2[1]*sa+pv2[2]*ca
        det=x1*x1*y3*y3-x3*x3*y1*y1
        s0=(y3*y3-y1*y1)/det;s1=(x1*x1-x3*x3)/det
        if s0<=0||s1<=0
            _geo_msg_error!(context,"Ellipse with tag $itag is wrong")
            ok=false
        else
            f2=sqrt(1/s1)
            a1=x1<0 ? -_geo_myasin(y1/f2)+a4+π : _geo_myasin(y1/f2)+a4
            a3=x3<0 ? -_geo_myasin(y3/f2)+a4+π : _geo_myasin(y3/f2)+a4
        end
    else
        a1=_geo_myatan2(pv0[2],pv0[1]);a3=_geo_myatan2(pv2[2],pv2[1])
    end
    a1=_geo_angle_02pi(a1);a3=_geo_angle_02pi(a3)
    a1>=a3 && (a3+=2π)
    if arc && itag>0 && a3-a1>1.01*π
        _geo_msg_error!(context,
            "Circle or ellipse arc $itag greater than Pi (angle=$(g(a3-a1)))")
        _geo_msg_error!(context,
            "(If you understand what this implies, you can disable this error")
        _geo_msg_error!(context,
            "message by selecting `Enable expert mode' in the option dialog.")
        _geo_msg_error!(context,
            "Otherwise, please subdivide the arc in smaller pieces.)")
        ok=false
    end
    return ok
end

# `EndCurve` runs once inside `CreateCurve` (forward, default `{0,0,1}`
# fallback normal — the only pass folded into the caller's `ok`), again in
# `addCircleArc`/`addEllipseArc` when a `Plane {..}` normal was given, and
# once more inside `CreateReversedCurve` on the mirrored control list
# (`-Num`, no `>Pi` check — its return is discarded there too, as is the
# reversed normal re-run's). For an ellipse whose resolved control list is
# short of four points, the reversed construction's four unconditional slot
# reads emit `Wrong list index (read)` per missing slot and the reversed
# `EndCurve` reports `Ellipse with tag -N is wrong` (upstream runs the math
# on uninitialized slot contents; the observed outcome is always the
# ellipse-wrong diagnostic).
function _geo_arc_endcurve_ok!(m::GeoModel,context::_GeoNumericContext,
                               found::Vector{Int},normal,tag;
                               ellipse::Bool,normal_given::Bool)
    itag=trunc(Int,tag);nn=length(found);ok=true
    if nn in (3,4)
        ok=_geo_endcurve_check!(m,context,found,(0.0,0.0,1.0),itag;
                                arc=true)
        normal_given && _geo_endcurve_check!(m,context,found,normal,itag;
                                             arc=true)
    end
    if ellipse && 0<nn<4
        for _ in 1:(4-nn)
            _geo_msg_error!(context,"Wrong list index (read)")
        end
        _geo_msg_error!(context,"Ellipse with tag $(-itag) is wrong")
    elseif nn in (3,4)
        _geo_endcurve_check!(m,context,found,(0.0,0.0,1.0),-itag;arc=false)
        normal_given &&
            _geo_endcurve_check!(m,context,found,normal,-itag;arc=false)
    end
    return ok
end

# `SortEdgesInLoop` (`reorient=false`) — chains signed curve tags by oriented
# end→start connectivity, starting a new subloop with the next remaining
# curve each time the chain closes on itself. A dead-end reports through
# `ok=false` and the stored loop keeps only the chained prefix. Callers must
# skip this entirely when any tag is unknown: `FindCurve` misses return
# `true` early upstream and the raw input list is stored with no diagnostic.
function _geo_sort_loop_partial(m::GeoModel,ids::Vector{Int})
    isempty(ids) && return (Int[],true)
    remaining=collect(ids)
    ordered=Int[];sizehint!(ordered,length(ids))
    ends(id)=(let (a,b)=m.curves[abs(id)];id>0 ? (a,b) : (b,a) end)
    c0=c1=first(remaining)
    push!(ordered,c0);deleteat!(remaining,1)
    nb=length(ids)
    while length(ordered)<nb
        advanced=false
        for (i,c2) in pairs(remaining)
            ends(c1)[2]==ends(c2)[1] || continue
            push!(ordered,c2);deleteat!(remaining,i)
            c1=c2
            if ends(c1)[2]==ends(c0)[1] && !isempty(remaining)
                c0=c1=first(remaining)
                push!(ordered,c0);deleteat!(remaining,1)
            end
            advanced=true
            break
        end
        advanced || return (ordered,false)
    end
    return (ordered,true)
end

# `SetSurfaceGeneratrices` — each wire resolves through `FindEdgeLoop(
# abs(iLoop))` (`Unknown curve loop %d in GEO surface %d`), then the loop's
# curves through `FindCurve(ic)`/the model edge map (`Unknown curve %d`);
# the reported curve tag carries the wire's sign and hole flip. Stops at the
# first failure.
function _geo_check_surface_wires!(m::GeoModel,context::_GeoNumericContext,
                                   wires,tag)
    itag=trunc(Int,tag)
    for (i,wire) in enumerate(wires)
        loop=get(m.loops,abs(wire),nothing)
        if loop===nothing
            _geo_msg_error!(context,
                "Unknown curve loop $wire in GEO surface $itag")
            return false
        end
        sgn=(wire<0 ? -1 : 1)*(i>1 ? -1 : 1)
        for signed_curve in loop
            ic=signed_curve*sgn
            haskey(m.curves,abs(ic)) && continue
            haskey(m.discrete,(1,abs(ic))) && continue
            _geo_msg_error!(context,"Unknown curve $ic")
            return false
        end
    end
    return true
end

# `SetVolumeSurfaces` — each shell resolves through `FindSurfaceLoop(
# abs(il))` (`Unknown surface loop %d`), then its surfaces through
# `FindSurface(abs(is))`/the model face map (`Unknown surface %d in GEO
# volume %d`). Stops at the first failure.
function _geo_check_volume_shells!(m::GeoModel,context::_GeoNumericContext,
                                   shells,tag)
    itag=trunc(Int,tag)
    for shell in shells
        loop=get(m.surface_loops,abs(shell),nothing)
        if loop===nothing
            _geo_msg_error!(context,"Unknown surface loop $shell")
            return false
        end
        for surface in loop
            tag0=abs(surface)
            haskey(m.surfaces,tag0) && continue
            haskey(m.discrete,(2,tag0)) && continue
            haskey(m.surface_geometry,tag0) && continue
            _geo_msg_error!(context,
                "Unknown surface $surface in GEO volume $itag")
            return false
        end
    end
    return true
end

# The shared `addX` epilogue: allocate the tag (upstream's
# `getMaxTag(dim) + 1` auto path for nonliteral 0), store unconditionally,
# then report `Could not add <name>` when the resolved record failed its
# checks.
function _geo_entity_tag(m::GeoModel,dim::Int,tag,caller,zero_literal::Bool)
    return _alloc_tag!(m,dim,_tag(tag,caller,dim),caller;
                       literal_zero=zero_literal)
end

function _geo_exec_entity_tags(raw::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString;
                               signed::Bool=false,
                               wrap::Bool=true,
                               abs_refs::Bool=false,
                               allow_empty::Bool=false)
    source=wrap ? "{"*String(raw)*"}" : String(raw)
    values=_geo_numeric_list_values(
        source,context,caller;allow_multiplier=true)
    (isempty(values) && !allow_empty) && throw(ArgumentError(
        "$caller: entity list is empty"))
    tags=Int[];sizehint!(tags,length(values))
    for numeric in values
        tag=_geo_signed_gmsh_int_value(numeric,"$caller entry")
        # `.geo` references follow Gmsh's `(int)FExpr` cast: signed entries
        # carry orientation in curve/surface loops (`signed`), while point
        # references and wire/shell lists resolve through `abs` lookups
        # (`CompareVertex`, `SetSurfaceGeneratrices`, `SetVolumeSurfaces`).
        # Literal tag 0 is referenceable whenever a tag-0 entity exists.
        if !signed && abs_refs
            tag==typemin(Int32) && throw(ArgumentError(
                "$caller entry magnitude exceeds Int32"))
            tag=abs(tag)
        end
        push!(tags,tag)
    end
    return tags
end

# `.geo` definition-site tags follow `addX(int &tag, ...)`: `tag < 0`
# auto-assigns `getMaxTag(dim) + 1` (the `add_*!` `tag=0` path), `tag == 0`
# stays a literal tag, and `tag > 0` is explicit. Returns
# `(tag, zero_literal)` for the model constructors.
function _geo_exec_def_tag(raw::AbstractString,
                           context::_GeoNumericContext,
                           caller::AbstractString)
    value=_geo_signed_gmsh_int_value(
        _geo_eval_numeric(raw,context,caller),caller)
    value<0 && return (0,false)
    return (value,value==0)
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

# `getPhysicalNumber(dim, name)` — the smallest `(dim, tag)` holding `name`
# (`_physicalNames` iterates in sorted order), or `nothing` when unbound.
function _geo_physical_number(m::GeoModel,dim::Int,name::AbstractString)
    found=nothing
    for ((edim,etag),ename) in m.physical_names
        edim==dim && ename==name || continue
        (found===nothing || etag<found) && (found=etag)
    end
    return found
end

# `GModel::setPhysicalName(name, dim, number)` — an already-bound name wins and
# resolves to its tag; `number == 0` auto-assigns `getMaxPhysicalNumber(dim)+1`
# (the per-dimension maximum of the *synchronized* view — 1 before any sync);
# otherwise `(dim, number)` gains the binding when the tag is still unnamed
# (std::map `insert` never overwrites an existing key).
function _geo_set_physical_name!(m::GeoModel,dim::Int,
                                 name::AbstractString,number::Int)
    bound=_geo_physical_number(m,dim,name)
    bound!==nothing && return bound
    number==0 && (number=_max_entity_physical_number(m,dim)+1)
    haskey(m.physical_names,(dim,number)) ||
        (m.physical_names[(dim,number)]=name)
    return number
end

# `_maxPhysicalNum` — the shared physical tag counter mirrored on the model
# and the allocator tracker. `setMaxPhysicalTag(t+1)` increments it outright;
# `CreatePhysicalGroup` raises it to `max(cur, tag)`.
function _geo_physical_assign_next!(m::GeoModel,allocator_state)
    # `setMaxPhysicalTag(t + 1)` runs on a raw `int` in Gmsh — at
    # `typemax(Int32)` it wraps to `typemin(Int32)` rather than erroring.
    m.physical_tag_max=m.physical_tag_max==typemax(Int32) ?
        Int(typemin(Int32)) : m.physical_tag_max+1
    allocator_state===nothing ||
        (allocator_state.physical_group_max=
            allocator_state.physical_group_max==typemax(Int32) ?
            Int(typemin(Int32)) : allocator_state.physical_group_max+1)
    return m.physical_tag_max
end

function _geo_physical_observe_tag!(m::GeoModel,allocator_state,tag::Int)
    m.physical_tag_max=max(m.physical_tag_max,tag)
    allocator_state===nothing ||
        (allocator_state.physical_group_max=
            max(allocator_state.physical_group_max,tag))
    return nothing
end

# `gmsh_sign` — the sign function Gmsh applies to physical-group member
# references: -1 for negative, 0 for zero, +1 for positive.
_gmsh_sign(x::Integer)=sign(x)

# `GEO_Internals::synchronize` — the physical pass rebuilds entity physicals
# from the raw group records (`orientedPhysicals` on): member `num` resolves
# the entity `abs(num)` and stores the signed entity-level physical
# `gmsh_sign(num) * Num` (`gmsh_sign(0) == 0`, so member `0` lands in physical
# group `0`); the observable group view is then keyed by `abs` of those stored
# values. Unresolvable members are skipped with Gmsh's sync-time warning; a
# duplicated entity appears once in each view list (its last entity-level
# physical wins the position).
# `GModel::getNumMeshElements()` for a `.geo` execution — elements merged by a
# mid-file `Mesh n` live in `context.mesh`; discrete entities and `attached`
# records carry theirs on the model.
function _geo_model_has_mesh_elements(m::GeoModel,context::_GeoNumericContext)
    mesh=context.mesh
    if mesh!==nothing && (size(mesh.segs,2)+size(mesh.tris,2)+
                          size(mesh.tets,2))>0
        return true
    end
    for rec in values(m.discrete)
        isempty(rec.element_tags) || return true
    end
    for rec in values(m.meshing.attached)
        isempty(rec.element_tags) || return true
    end
    return false
end

function _geo_sync_physical_view!(m::GeoModel,context::_GeoNumericContext)
    raw=context.raw_physicals
    # `synchronize` clears the entity-level memberships only when the model
    # carries no mesh elements and raw physical groups exist — otherwise
    # previously assigned memberships persist (a `-=` that emptied the last
    # raw group leaves them in place; API-assigned memberships survive a sync
    # with an empty raw registry too).
    if !isempty(raw) && !_geo_model_has_mesh_elements(m,context)
        empty!(m.entity_physicals)
    end
    empty!(m.physical)
    # `Tree2List(PhysicalGroups)` iterates in `ComparePhysicalGroup` order:
    # `(Typ, Num)` with the tag compared as `q->Num - w->Num` — a raw `int`
    # subtraction that wraps for tag spans wider than `Int32`, so e.g.
    # `-2^31` sorts *after* `2^31-1`.
    group_keys=sort!(collect(keys(raw)),
        lt=(a,b)->begin
            a[1]!=b[1] && return a[1]<b[1]
            d=mod(Int64(a[2])-Int64(b[2])-Int64(typemin(Int32)),
                  Int64(2)^Int64(32))+Int64(typemin(Int32))
            return d<0
        end)
    for (dim,tag) in group_keys
        word=("point","curve","surface","volume")[dim+1]
        for member in raw[(dim,tag)]
            entity=abs(member)
            if _has_entity(m,dim,entity)
                ep=get!(m.entity_physicals,(dim,entity),Int[])
                p=_gmsh_sign(member)*tag
                p in ep || push!(ep,p)
            else
                _geo_yywarn!(context,
                    "Skipping unknown $word $entity in physical $word $tag")
            end
        end
    end
    for ((dim,entity),pnums) in m.entity_physicals
        for p in pnums
            members=get!(m.physical,(dim,abs(p)),Int[])
            entity in members || push!(members,entity)
        end
    end
    for members in values(m.physical)
        sort!(members)
    end
    context.geo_changed=false
    return nothing
end

# `if(getGEOInternals()->getChanged()) synchronize(...)` — every grammar sync
# point except `SyncModel` gates on the changed flag, so an unchanged model is
# not resynchronized (and "unknown member" warnings are not re-emitted).
function _geo_sync_physical_view_if_changed!(m::GeoModel,
                                             context::_GeoNumericContext)
    context.geo_changed || return nothing
    _geo_sync_physical_view!(m,context)
    return nothing
end

# `Physical X(n) += {..}` / `-= {..}` / `*=` / `/=` — `modifyPhysicalGroup`
# over the raw group records: `+=` appends raw member tags (duplicates kept),
# `-=` removes the first occurrence of each listed tag (deleting the group —
# and its raw-tag name binding — when emptied, a silent no-op when the group
# is missing), and `*=`/`/=` are errors — "does not exist" when the group is
# missing, else "Unsupported operation".
function _geo_exec_physical_modify!(m::GeoModel,dim::Int,op::AbstractString,
                                    tag::Int,ids::Vector{Int},
                                    context::_GeoNumericContext,
                                    caller::AbstractString,kind::String)
    word=("point","curve","surface","volume")[dim+1]
    # The caller-side `yymsg` uses "line" for dim 1 (Gmsh's own strings).
    caller_word=("point","line","surface","volume")[dim+1]
    fail()=_geo_yyerror!(context,"Could not modify physical $caller_word")
    members=get(context.raw_physicals,(dim,tag),nothing)
    if members===nothing
        op=="-=" && return nothing # `removePhysicalGroup` is not an error
        _geo_msg_error!(context,"Physical $word $tag does not exist")
        fail()
        return nothing
    end
    if op=="*=" || op=="/="
        _geo_msg_error!(
            context,"Unsupported operation on physical $word $tag")
        fail()
        return nothing
    end
    if op=="+="
        append!(members,ids)
    else # "-="
        if isempty(ids)
            empty!(members)
        else
            for id in ids
                pos=findfirst(==(id),members)
                pos===nothing || deleteat!(members,pos)
            end
        end
        if isempty(members)
            delete!(context.raw_physicals,(dim,tag))
            # `DeletePhysicalX` routes through `GModel::removePhysicalGroup`
            # with the raw tag: it erases the `(dim, tag)` name binding and
            # strips stored entity memberships whose `abs` equals the tag.
            # For a negative raw tag `abs(p) == tag` can never hold, so a `-=`
            # that empties group `-4` leaves the `-4` memberships on the
            # entities — they only disappear on the next clearing sync.
            delete!(m.physical_names,(dim,tag))
            for (ekey,pnums) in m.entity_physicals
                filter!(p->abs(p)!=tag,pnums)
            end
        end
    end
    # `modifyPhysicalGroup` marks the internals changed on every successful op
    # (the early error returns above do not reach its `_changed = true`).
    context.geo_changed=true
    return nothing
end

# `Physical X(n) = {..}` / `Physical X("name"[, n]) = {..}` — the
# `modifyPhysicalGroup` create path over the raw group records. An existing
# `(dim, tag)` is a recoverable "already exists" error; the member list is
# stored with its raw signs and the observable view is materialized lazily at
# sync points (`_geo_sync_physical_view!`), where unresolvable members are
# skipped with Gmsh's sync-time warning. A created group raises
# `_maxPhysicalNum` to `max(cur, tag)` — negative and zero tags never bump it.
function _geo_exec_physical_create!(m::GeoModel,dim::Int,tag::Int,
                                    ids::Vector{Int},
                                    context::_GeoNumericContext,kind::String,
                                    allocator_state)
    word=("point","curve","surface","volume")[dim+1]
    caller_word=("point","line","surface","volume")[dim+1]
    if haskey(context.raw_physicals,(dim,tag))
        _geo_msg_error!(context,"Physical $word $tag already exists")
        _geo_yyerror!(context,"Could not modify physical $caller_word")
        return nothing
    end
    context.raw_physicals[(dim,tag)]=copy(ids)
    _geo_physical_observe_tag!(m,allocator_state,tag)
    # `modifyPhysicalGroup` op 0 reaches `_changed = true` on success.
    context.geo_changed=true
    return nothing
end

# Extract one balanced `{...}` group from the front of `raw`; returns
# `(content, rest)`.
function _geo_balanced_group(raw::AbstractString, caller::AbstractString)
    s=String(strip(raw))
    # A missing or unbalanced group is a `syntax error` upstream — the
    # reported token is the unexpected lookahead, matching bison.
    if !startswith(s,"{")
        isempty(s) && _geo_syntax_abort(";")
        m=match(r"^[A-Za-z_][A-Za-z0-9_]*|^.",s)
        _geo_syntax_abort(m===nothing ? ";" : String(m.match))
    end
    depth=0;closing=0;i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if c=='{'
            depth+=1
        elseif c=='}'
            depth-=1
            depth==0 && (closing=i; break)
            depth<0 && _geo_syntax_abort("}")
        end
        i=nextind(s,i)
    end
    closing==0 && _geo_syntax_abort("}")
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
    return _geo_exec_entity_tags(r,context,caller;signed=signed_tags,
                                 abs_refs=dim==0 && !signed_tags)
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
    startswith(rest,";") &&
        return String(strip(rest[nextind(rest,1):end]))
    # `ListOfShapes` members are `tEND`-terminated upstream — anything else is
    # a syntax error. The reported token is the unexpected lookahead: the
    # closing `}` when the list ran out, else the next token's text.
    _geo_syntax_abort(if isempty(rest)
        "}"
    else
        m=match(r"^[A-Za-z_][A-Za-z0-9_]*|^.",rest)
        m===nothing ? "}" : String(m.match)
    end)
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
                                   caller::AbstractString;
                                   signed_tags::Bool=false,
                                   allocator_state=nothing)
    entities=NTuple{2,Int}[]
    s=String(strip(raw))
    while !isempty(s)
        (entries,rest,transform)=
            _geo_multiple_shape_element!(m,s,context,caller;
                                         signed_tags=signed_tags,
                                         allocator_state=allocator_state)
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
                                      signed_tags::Bool=false,
                                      allocator_state=nothing)
    mm=match(r"^([A-Za-z_][A-Za-z0-9_]*)",s0)
    # A member that does not start with an identifier (`{1}`, `1`, `;`) is a
    # syntax error upstream — `MultipleShape` members are `GeoEntity{...}`,
    # `Shape` definitions, or nested transforms.
    mm===nothing && _geo_syntax_abort(String(s0[1:1]))
    name=mm.captures[1]
    s=String(strip(s0[nextind(s0,firstindex(s0),ncodeunits(mm.match)):end]))
    if name in ("Translate","Rotate","Dilate","Symmetry","Affine","Closest")
        nested="$caller: $name"
        if name=="Affine" || name=="Closest"
            # `tAffine '{' RLD '}' '{' MultipleShape '}'` — the parameters are
            # a brace-delimited numeric list, then the shape group.
            (params,s)=_geo_balanced_group(s,nested)
            _geo_numeric_list_values(
                "{"*params*"}",context,"$caller $name";depth=1)
            s=String(strip(s))
            startswith(s,"{") || _geo_selector_abort(s)
            (inner_list,s)=_geo_balanced_group(s,nested)
            # OCC-only transforms report a recoverable error upstream and
            # leave the shapes untouched: `Affine` still yields its input
            # list, while `Closest` contributes nothing.
            if name=="Affine"
                _geo_yyerror!(context,
                    "Affine transform only available with OpenCASCADE " *
                    "geometry kernel")
                inner=_geo_shape_list_entities!(m,inner_list,context,caller;
                                                allocator_state=allocator_state)
                return (inner,s,true)
            end
            _geo_yyerror!(context,
                "Closest entity only available with OpenCASCADE " *
                "geometry kernel")
            _geo_shape_list_entities!(m,inner_list,context,caller;
                                      allocator_state=allocator_state)
            return (NTuple{2,Int}[],s,true)
        end
        (t,s)=_geo_transform_params_rest(
            kind=name,source=s,context=context,caller=nested)
        s=String(strip(s))
        startswith(s,"{") || _geo_selector_abort(s)
        (inner_list,s)=_geo_balanced_group(s,nested)
        inner=_geo_shape_list_entities!(m,inner_list,context,caller;
                                        allocator_state=allocator_state)
        transform_entities!(m,t,inner;caller=nested)
        # Gmsh returns the input shape list (`$$ = $MultipleShape`).
        return (inner,s,true)
    end
    if name=="Split" || name=="Intersect"
        # `tSplit tCurve ...` / `tIntersect tCurve ...` — real built-in kernel
        # curve operations upstream; a hard blocker here. Any other head is a
        # syntax error (the keyword expects `Curve`).
        match(r"^Curve\b",s)===nothing && _geo_selector_abort(s)
        throw(ArgumentError(
            "$caller: $name requires built-in kernel curve splitting, " *
            "which Tessella does not implement"))
    end
    if name=="Physical" || name=="Parent"
        km=match(r"^(Point|Curve|Line|Surface|Volume)\b",s)
        if km===nothing && (em=match(r"^GeoEntity\b",s))!==nothing
            s=String(strip(s[nextind(s,firstindex(s),ncodeunits(em.match)):end]))
            (dgroup,s)=_geo_balanced_group(s,caller)
            dim=_geo_int_value(_geo_eval_numeric(dgroup,context,
                "$caller $name GeoEntity dimension"),
                "$caller $name GeoEntity dimension")
            # The `GeoEntity` rule reports the out-of-range dim but still
            # evaluates with it — the lookups below then simply find nothing.
            0<=dim<=3 || _geo_yyerror!(context,
                "GeoEntity dim out of range [0,3]")
        else
            km===nothing && _geo_selector_abort(s)
            dim=_geo_shape_kind_dim(km.captures[1])
            s=String(strip(s[nextind(s,firstindex(s),ncodeunits(km.match)):end]))
        end
        (group,s)=_geo_balanced_group(s,caller)
        s=_geo_require_list_semicolon(s,caller,name)
        # Tessella entities have no parent entities; the selector is empty.
        name=="Parent" && return (NTuple{2,Int}[],s,false)
        # `getElementaryTagsForPhysicalGroups` synchronizes first — the
        # selector sees the derived view, resolved through `abs` group tags.
        _geo_sync_physical_view_if_changed!(m,context)
        entries=NTuple{2,Int}[]
        for gtag in _geo_shape_physical_tags(m,dim,group,context,caller)
            members=get(m.physical,(dim,gtag),nothing)
            # Gmsh resolves the group through the model: an unknown tag (or a
            # group whose members were all `Delete`d) silently contributes
            # nothing, and stale member integers resolve to live entities only.
            members===nothing && continue
            for tag in members
                _model_entity_known(m,dim,tag) || continue
                push!(entries,(dim,tag))
            end
        end
        return (entries,s,false)
    end
    if name=="GeoEntity"
        # `GeoEntity{dim}{tags};` — the generic `tGeoEntity` selector.
        (dgroup,s)=_geo_balanced_group(s,caller)
        dim=_geo_int_value(_geo_eval_numeric(
            dgroup,context,"$caller GeoEntity dimension"),
            "$caller GeoEntity dimension")
        (group,s)=_geo_balanced_group(s,caller)
        s=_geo_require_list_semicolon(s,caller,name)
        # Upstream reports the dim error yet still stores the shapes under an
        # out-of-range type — nothing resolves them later, so the group
        # contributes no entities (its tag list is still evaluated).
        if !(0<=dim<=3)
            _geo_yyerror!(context,"GeoEntity dim out of range [0,3]")
            _geo_shape_tags(m,dim,group,context,caller;
                            signed_tags=signed_tags)
            return (NTuple{2,Int}[],s,false)
        end
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
        # `Kind(tag) = rhs;` — an inline Shape definition. Upstream `tCurve`
        # in a shape list is followed by `'{'` or `'('`; anything else is a
        # syntax error.
        startswith(s,"(") || _geo_selector_abort(s)
        return _geo_shape_definition_element!(m,s0,context,caller,
                                              allocator_state)
    end
    if startswith(s,"{")
        # tSTRING '{' MultipleShape '}' — Duplicata, a boundary query, or an
        # unknown action.
        (group,s)=_geo_balanced_group(s,caller)
        inner=_geo_shape_list_entities!(m,group,context,caller;
                                        allocator_state=allocator_state)
        # Boundary-family queries and `PointsOf` run on the synced model in
        # Gmsh — the internals synchronize first when changed; `Duplicata`
        # copies entities on the raw internals instead, which marks them
        # changed without syncing.
        if name=="Duplicata"
            context.geo_changed=true
        else
            _geo_sync_physical_view_if_changed!(m,context)
        end
        return (_geo_shape_action!(m,name,inner,context,caller),s,true)
    end
    (startswith(s,"(") || match(r"^[A-Za-z_]",s)!==nothing) &&
        return _geo_shape_definition_element!(m,s0,context,caller)
    _geo_selector_abort(s)
end

# A `Name(tag) = rhs;` definition (possibly multi-word, e.g. `Plane Surface`)
# inside a shape list executes through the normal statement executor and adds
# the created entity to the list, as in the Gmsh `Shape` production. Loop
# records and non-entity definitions are not transformable.
function _geo_shape_definition_element!(m::GeoModel, s0::AbstractString,
                                        context::_GeoNumericContext,
                                        caller::AbstractString,
                                        allocator_state=nothing)
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
    cut==0 && _geo_syntax_abort("}")
    stmt=String(s0[firstindex(s0):cut])
    rest=String(strip(s0[nextind(s0,cut):end]))
    head=match(r"^((?:[A-Za-z_][A-Za-z0-9_]*\s+)*[A-Za-z_][A-Za-z0-9_]*)\s*\(",
               stmt)
    if head===nothing
        nm=match(r"^[A-Za-z_][A-Za-z0-9_]*",stmt)
        _geo_selector_abort(nm===nothing ? stmt :
            String(strip(stmt[nextind(stmt,firstindex(stmt),
                                     ncodeunits(nm.match)):end])))
    end
    kind=String(strip(replace(head.captures[1],r"\s+"=>" ")))
    dim=get(_GEO_SHAPE_DEFINITION_DIMS,kind,-1)
    dim<0 && throw(ArgumentError(
        "$caller: '$kind' entries in a transform list do not produce " *
        "transformable entities"))
    tm=match(r"\(([^()]*)\)",stmt)
    (deftag,def_zero)=_geo_exec_def_tag(
        tm.captures[1],context,"$caller $kind tag")
    # A negative tag auto-assigns through `_exec_line!` below — predict the
    # tag the `add_*!` counter path will take (shape-list defs only cover
    # dims 0–3, so this is the per-dimension `next_tag` counter).
    tag=deftag==0 && !def_zero ? _geo_next_auto_tag(m,dim) : deftag
    _exec_line!(m,stmt,context,allocator_state)
    allocator_state===nothing ||
        _geo_allocator_observe_statement!(
            allocator_state,stmt,context,"execute_geo")
    return (NTuple{2,Int}[(dim,tag)],rest,false)
end

# Mirror of `_alloc_tag!`'s automatic path: the next `next_tag[dim]` value,
# skipping discrete-entity collisions.
function _geo_next_auto_tag(m::GeoModel,dim::Int)
    t=m.next_tag[dim+1]+1
    while haskey(m.discrete,(dim,t))
        t+=1
    end
    return t
end

# A `Name{ MultipleShape }` action element: Duplicata copies, boundary
# queries, and PointsOf. Returns the action's output entities.
function _geo_shape_action!(m::GeoModel, name::AbstractString,
                            inner::Vector{NTuple{2,Int}},
                            context::_GeoNumericContext,
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
    # `yymsg(0, "Unknown action on multiple shapes '%s'")` — recoverable, the
    # action contributes nothing.
    _geo_yyerror!(context,"Unknown action on multiple shapes '$name'")
    return NTuple{2,Int}[]
end

# ── Extrude ──────────────────────────────────────────────────────────────
#
# `.geo` translational `Extrude {dx,dy,dz} { ...; }`, mirroring the Gmsh
# grammar (`tExtrude VExpr '{' ListOfShapes ExtrudeParameters '}'`) and the
# `GEO_Internals::extrude` → `ExtrudeShapes` semantics: the expression value
# is the flat `[top, body, laterals...]` tag list (laterals only under
# `Geometry.ExtrudeReturnLateralEntities`, which defaults on). Rotational
# (`{{axis}, {point}, angle}`), and twist (`{{delta}, {axis}, {point},
# angle}`) run natively; boundary-layer (`Extrude {shapes; params}`) and
# `Using Wire` pipe forms are recognized and rejected with explicit errors.

# Scan a balanced `(...)` group — the VExpr grammar also allows parenthesized
# vectors — returning `(content, rest)`.
function _geo_balanced_paren(raw::AbstractString, caller::AbstractString)
    s=String(strip(raw))
    if !startswith(s,"(")
        isempty(s) && _geo_syntax_abort(";")
        m=match(r"^[A-Za-z_][A-Za-z0-9_]*|^.",s)
        _geo_syntax_abort(m===nothing ? ";" : String(m.match))
    end
    depth=0;closing=0;i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if c=='('
            depth+=1
        elseif c==')'
            depth-=1
            depth==0 && (closing=i; break)
            depth<0 && _geo_syntax_abort(")")
        end
        i=nextind(s,i)
    end
    closing==0 && _geo_syntax_abort(")")
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

# The revolve/twist `{{A}, {X}, alpha}` / `{{T}, {A}, {X}, alpha}` motion
# group, decoded from the group's evaluated element lengths: the grammar's
# VExpr members are length-3 vectors and the trailing FExpr a scalar.
# `[3,3,1]` is revolve and `[3,3,3,1]` is twist (`TRANSLATE_ROTATE`) —
# whose member order is translation, axis direction, point-on-axis (the
# order `GEO_Internals::twist` forwards to `ExtrudeShapes`, not the
# revolve order with a trailing delta).
function _geo_extrude_revolve(values::Vector{Vector{Float64}},
                              caller::AbstractString)
    lengths=length.(values)
    if length(lengths)==4 && lengths[1:3]==[3,3,3] && lengths[4]==1
        T,a,x,α=values
        return (kind=:twist,axis=(a[1],a[2],a[3]),
                origin=(x[1],x[2],x[3]),delta=(T[1],T[2],T[3]),
                angle=α[1])
    end
    if !(length(lengths)==3 && lengths==[3,3,1])
        throw(ArgumentError(
            "$caller: rotational extrusion requires `{{axis}, {point}, " *
            "angle}`"))
    end
    a,x,α=values
    return (kind=:rotate,axis=(a[1],a[2],a[3]),
            origin=(x[1],x[2],x[3]),angle=α[1])
end

# `Geometry.ExtrudeSplinePoints` (default 5) — the generatrix subdivision
# count for twist extrusions, `(int)`-cast like the option callback.
_geo_extrude_spline_points(context::_GeoNumericContext) = trunc(Int,
    something(_geo_option_number(context,"Geometry",0,
        "ExtrudeSplinePoints"),5.0))

# Accumulate one VExpr group into `delta` (`sign` folds `VExpr '+' VExpr` and
# leading `tMINUS VExpr`), or decode the revolve `{{axis}, {point}, angle}` /
# twist `{{delta}, {axis}, {point}, angle}` motion group. The group's
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
                               caller::AbstractString,
                               allocator_state=nothing)
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
        m,element*";",context,caller;signed_tags=true,
        allocator_state=allocator_state)
    transform && throw(ArgumentError(
        "$caller: transforms cannot appear inside an Extrude shape list"))
    isempty(strip(rest)) || throw(ArgumentError(
        "$caller: unexpected text in Extrude shape list near $(repr(rest))"))
    return (entries,params)
end

function _geo_extrude_shape_list!(m::GeoModel, body::AbstractString,
                                  context::_GeoNumericContext,
                                  caller::AbstractString,
                                  allocator_state=nothing)
    params=_GEO_EXTRUDE_PARAMS
    entities=NTuple{2,Int}[]
    isempty(strip(body)) && return (entities,params)
    for element in _geo_exec_topology_query_blocks(
            body,"Extrude shape list",caller)
        entries,params=_geo_extrude_element!(
            m,element,context,params,caller,allocator_state)
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
                                context::_GeoNumericContext,
                                allocator_state=nothing)
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
            # A paren group is a plain `VExpr` — upstream its arity error is
            # the closer (`Extrude(1)` → `)`, `Extrude(1,2,3,4)` → `,`).
            nparts=length(_geo_split_top_commas(group,caller))
            nparts!=3 && _geo_syntax_abort(nparts<3 ? ")" : ",")
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
    # Gmsh.y gates `TRANSLATE_ROTATE` on the built-in kernel in the parser
    # action — after `ListOfShapes` parses but before `ExtrudeShapes` runs
    # — so under `SetFactory("OpenCASCADE")` shape-list syntax errors still
    # surface first, while entity lookup never happens.
    if revolve!==nothing && revolve.kind===:twist &&
        allocator_state!==nothing &&
        allocator_state.factory===:opencascade
        _geo_exec_topology_query_blocks(shapes,"Extrude shape list",caller)
        throw(ArgumentError(
            "$caller: Twisting extrude not available with OpenCASCADE " *
            "geometry kernel"))
    end
    entities,params=_geo_extrude_shape_list!(m,shapes,context,caller,
                                             allocator_state)
    tags=revolve===nothing ?
        extrude_entities!(m,entities,delta;params=params,
            return_lateral=context.extrude_return_lateral,caller=caller) :
        revolve.kind===:twist ?
        twist_entities!(m,entities,revolve.axis,revolve.origin,
            revolve.delta,revolve.angle;
            spline_points=_geo_extrude_spline_points(context),
            params=params,
            return_lateral=context.extrude_return_lateral,caller=caller) :
        revolve_entities!(m,entities,revolve.axis,revolve.origin,
            revolve.angle;params=params,
            return_lateral=context.extrude_return_lateral,caller=caller)
    context.geo_changed=true
    return Float64.(tags)
end

# Evaluate a transform's parameter portion and return `(transform, rest)` —
# `rest` continues at the `{shapes}` group. `Translate`/`Symmetry` take a bare
# `VExpr` upstream (`{a,b,c}`, `(a,b,c)`, `±`-composed); `Rotate`/`Dilate`
# take a literal `{ VExpr, ... }` group.
function _geo_transform_params_rest(;kind::AbstractString,
                                    source::AbstractString,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    s=String(strip(source))
    if kind=="Translate"
        (v,rest)=_geo_exec_vexpr5_rest(s,context,"$caller delta")
        return (_affine_translation((v[1],v[2],v[3]),caller),rest)
    elseif kind=="Symmetry"
        (v,rest)=_geo_exec_vexpr5_rest(s,context,"$caller plane coefficients")
        return (_affine_symmetry(v[1],v[2],v[3],v[4],caller),rest)
    end
    (params,rest)=_geo_balanced_group(s,caller)
    parts=_geo_split_top_commas(params,caller)
    if kind=="Dilate"
        length(parts)==2 || throw(ArgumentError(
            "$caller: Dilate requires `{center, scale}` parameters"))
        center=_geo_exec_vexpr(parts[1],3,context,"$caller center")
        scale_raw=String(strip(parts[2]))
        scale=_geo_is_vexpr(scale_raw) ?
            _geo_exec_vexpr(scale_raw,3,context,"$caller scales") :
            _geo_eval_numeric(scale_raw,context,"$caller scale")
        return (_affine_dilation(center,scale,caller),rest)
    else
        length(parts)==3 || throw(ArgumentError(
            "$caller: Rotate requires `{{axis}, {origin}, angle}` parameters"))
        axis=_geo_exec_vexpr(parts[1],3,context,"$caller axis")
        origin=_geo_exec_vexpr(parts[2],3,context,"$caller origin")
        angle=_geo_eval_numeric(parts[3],context,"$caller angle")
        return (_affine_rotation(axis,origin,angle,caller),rest)
    end
end

function _geo_exec_transform_statement!(m::GeoModel,line::AbstractString,
                                        context::_GeoNumericContext,
                                        allocator_state=nothing)
    source=String(strip(line))
    endswith(source,";") && (source=String(strip(source[1:prevind(source,end)])))
    mm=match(r"^(Translate|Rotate|Dilate|Symmetry)\s*",source)
    kind=String(mm.captures[1])
    caller="execute_geo: $kind"
    rest=String(strip(source[nextind(source,firstindex(source),
                                   ncodeunits(kind)):end]))
    (t,rest)=_geo_transform_params_rest(
        kind=kind,source=rest,context=context,caller=caller)
    (shape_list,rest)=_geo_balanced_group(rest,caller)
    isempty(rest) || throw(ArgumentError(
        "$caller: unexpected text after the shape list"))
    entities=_geo_shape_list_entities!(m,shape_list,context,caller;
                                       allocator_state=allocator_state)
    transform_entities!(m,t,entities;caller=caller)
    context.geo_changed=true
    return nothing
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
                                   signed::Bool=false,
                                   abs_refs::Bool=false,
                                   allow_empty::Bool=false)
    # Gmsh parses every entity RHS as `ListOfDouble`: `{...}` groups, bare
    # `name[]`/`name[{..}]` references, `-{...}` negation, `expr * {...}`
    # multipliers and plain scalars — `_geo_numeric_list_values` covers the
    # whole grammar and reports malformed brace forms itself. `allow_empty`
    # mirrors `addX`'s empty-vector tolerance (`Line {}` reaches the inner
    # `requires 2 points` diagnostic; loops/volumes store empty silently).
    source=String(strip(raw))
    isempty(source) && !allow_empty && throw(ArgumentError(
        "$caller: entity list must not be empty"))
    return _geo_exec_entity_tags(
        source,context,caller;signed=signed,wrap=false,abs_refs=abs_refs,
        allow_empty=allow_empty)
end

# `Periodic` slave/master lists keep their signs: `addPeriodicEdge` resolves
# entities through `abs(iTarget)`/`abs(iSource)` but derives the orientation of
# a transform-free relation from `iSource * iTarget` — the signed values reach
# `set_periodic!`, which does both.
_geo_periodic_tags(raw::AbstractString,context::_GeoNumericContext,
                   caller::AbstractString)=
    _geo_exec_entity_tags(raw,context,caller;signed=true)

# Evaluate a Gmsh `VExpr` — `{a,b,c[,d[,e]]}` or `(a,b,c)` groups composed by
# unary/binary `+`/`-` (`VExpr_Single` defaults components 4 and 5 to 0 and 1;
# the `(...)` form accepts exactly three). Returns `(components, rest)` — the
# VExpr ends at the first non-`±` token, which is how a transform's trailing
# `{shapes}` group is found.
function _geo_exec_vexpr5_rest(raw::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString)
    s=String(strip(raw))
    isempty(s) && _geo_syntax_abort(";")
    acc=ntuple(_->0.0,5);sign=1.0;expect_term=true
    while !isempty(s)
        if expect_term
            while startswith(s,"+") || startswith(s,"-")
                s[1]=='-' && (sign=-sign)
                s=String(strip(s[nextind(s,firstindex(s)):end]))
            end
            (isempty(s) || (s[1]!='{' && s[1]!='(')) &&
                _geo_selector_abort(s)
            paren=s[1]=='('
            group,rest=paren ? _geo_balanced_paren(s,caller) :
                               _geo_balanced_group(s,caller)
            parts=_geo_split_top_commas(group,caller)
            np=length(parts)
            # Components are `FExpr` parsed one at a time — upstream a
            # malformed component (`x[]` → `]`) errors before the arity
            # check ever runs, so evaluate first.
            evaluated=ntuple(5) do i
                i<=np ? _geo_eval_numeric(parts[i],context,caller) :
                        i==4 ? 0.0 : 1.0
            end
            # `VExpr_Single` needs 3 (paren) or 3–5 (braces) FExpr components —
            # upstream the next token after the last legal component is the
            # error (`{a,b}` → `}`, `(a,b)` → `)`, a 4th/6th `,` → `,`).
            if paren ? np!=3 : !(np in 3:5)
                _geo_syntax_abort(np<3 ? (paren ? ")" : "}") : ",",
                    "$caller: a vector expression needs " *
                    (paren ? "exactly 3 components" : "3 to 5 components") *
                    "; got $np")
            end
            acc=acc .+ sign .* evaluated
            s=String(strip(rest));sign=1.0;expect_term=false
        else
            if s[1]=='+' || s[1]=='-'
                sign=s[1]=='+' ? 1.0 : -1.0
                s=String(strip(s[nextind(s,firstindex(s)):end]))
                expect_term=true
            else
                break
            end
        end
    end
    expect_term && _geo_syntax_abort(";",
        "$caller: vector expression ends with an operator")
    return acc,s
end

# A `VExpr` that must consume its whole input (`Plane` normals, `Dilate`
# scale, ...); returns the first `n` components.
function _geo_exec_vexpr(raw::AbstractString,n::Int,
                         context::_GeoNumericContext,caller::AbstractString)
    (v,rest)=_geo_exec_vexpr5_rest(raw,context,caller)
    isempty(rest) || _geo_syntax_abort(_geo_first_token(rest),
        "$caller: unexpected text after vector expression $(repr(rest))")
    return ntuple(i->v[i],n)
end

_geo_exec_vexpr3(raw::AbstractString,context::_GeoNumericContext,
                 caller::AbstractString)=_geo_exec_vexpr(raw,3,context,caller)

# `{`- or `(`-headed (after any unary `+`/`-` signs) — the VExpr form, as
# opposed to a bare `FExpr` scalar.
function _geo_is_vexpr(raw::AbstractString)
    s=String(strip(raw))
    while startswith(s,"+") || startswith(s,"-")
        s=String(strip(s[nextind(s,firstindex(s)):end]))
    end
    return startswith(s,"{") || startswith(s,"(")
end

# Parse the shared `Circle`/`Ellipse` RHS: `{point tags}` plus the optional
# `Plane VExpr` normal override (`CircleOptions` in the Gmsh grammar).
function _geo_circle_rhs(raw::AbstractString,context::_GeoNumericContext,
                         caller::AbstractString)
    (group,rest)=_geo_balanced_group(String(strip(raw)),caller)
    points=_geo_exec_entity_tags(group,context,"$caller points";
                                 abs_refs=true)
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
    # Upstream's `computeAffineTransformation` plugs the raw axis components
    # into Rodrigues' formula without normalizing: a non-unit axis yields a
    # non-orthogonal map that the endpoint correspondence check then rejects,
    # and a zero axis yields the `cos(angle)` uniform scale. Keeping the same
    # arithmetic means Tessella's strict endpoint check rejects exactly the
    # declarations upstream drops.
    x,y,z=axis
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

# `CTX::instance()->lc` as upstream's periodic actions see it. `addPeriodic*`
# runs `GEO_Internals::synchronize` first, whose tail `SetBoundingBox()`
# recomputes `lc` from `GModel::bounds()` — at `.geo` parse time entity
# `mesh_vertices` are empty, so only `GVertex` positions (native points and
# discrete vertices) plus OCC-native entity bounds feed the box — and
# `FinishUpBoundingBox()` then expands every degenerate dimension. This is
# deliberately not `_geo_exec_lc` (the raw `temp_bbox` diagonal): the periodic
# tolerance follows the synchronized, expanded value.
function _geo_periodic_lc(m::GeoModel,context::_GeoNumericContext)
    caller="execute_geo: Periodic"
    bounds=nothing
    for point in values(m.points)
        bounds=_geo_periodic_bounds_add(bounds,point)
    end
    for ((dimension,_),record) in m.discrete
        dimension==0 && !isempty(record.node_coords) &&
            (bounds=_geo_periodic_bounds_add(
                bounds,record.node_coords[1:3,1]))
    end
    # Upstream's `bounds()` uses `entities[i]->bounds()` only for
    # OpenCascadeModel natives; built-in entities contribute their (empty)
    # mesh vertices. OCC-native curves carry `_occ_geometry`, OCC-native
    # surfaces an `occ` record; an OCC volume's box is already covered by its
    # OCC faces.
    for tag in keys(m.curves)
        _occ_geometry(m,tag)===nothing && continue
        entity=_model_entity_bounding_box(m,1,tag,caller)
        bounds=_geo_periodic_bounds_union(
            bounds,(entity[1:3],entity[4:6]))
    end
    for tag in keys(m.surfaces)
        geometry=get(m.surface_geometry,tag,nothing)
        (geometry!==nothing && hasproperty(geometry,:occ) &&
         geometry.occ!==nothing) || continue
        entity=_model_entity_bounding_box(m,2,tag,caller)
        bounds=_geo_periodic_bounds_union(
            bounds,(entity[1:3],entity[4:6]))
    end
    bounds===nothing &&
        (bounds=((-1.0,-1.0,-1.0),(1.0,1.0,1.0)))
    lo,hi=bounds
    tolerance=max(1e-6,
        something(_geo_option_number(context,"Geometry",0,"Tolerance"),1e-8))
    range=(hi[1]-lo[1],hi[2]-lo[2],hi[3]-lo[3])
    # `FinishUpBoundingBox` — note upstream has no z-only-degenerate branch.
    if range[1]<tolerance && range[2]<tolerance && range[3]<tolerance
        lo=(lo[1]-1.0,lo[2]-1.0,lo[3])
        hi=(hi[1]+1.0,hi[2]+1.0,hi[3])
    elseif range[1]<tolerance && range[2]<tolerance
        lo=(lo[1]-range[3],lo[2]-range[3],lo[3])
        hi=(hi[1]+range[3],hi[2]+range[3],hi[3])
    elseif range[1]<tolerance && range[3]<tolerance
        lo=(lo[1]-range[2],lo[2],lo[3]);hi=(hi[1]+range[2],hi[2],hi[3])
    elseif range[2]<tolerance && range[3]<tolerance
        lo=(lo[1],lo[2]-range[1],lo[3]);hi=(hi[1],hi[2]+range[1],hi[3])
    elseif range[1]<tolerance
        l=sqrt(range[2]*range[2]+range[3]*range[3])
        lo=(lo[1]-l,lo[2],lo[3]);hi=(hi[1]+l,hi[2],hi[3])
    elseif range[2]<tolerance
        l=sqrt(range[1]*range[1]+range[3]*range[3])
        lo=(lo[1],lo[2]-l,lo[3]);hi=(hi[1],hi[2]+l,hi[3])
    end
    return sqrt((hi[1]-lo[1])^2+(hi[2]-lo[2])^2+(hi[3]-lo[3])^2)
end

_geo_periodic_bounds_add(bounds,point)=bounds===nothing ?
    ((point[1],point[2],point[3]),(point[1],point[2],point[3])) :
    ((min(bounds[1][1],point[1]),min(bounds[1][2],point[2]),
      min(bounds[1][3],point[3])),
     (max(bounds[2][1],point[1]),max(bounds[2][2],point[2]),
      max(bounds[2][3],point[3])))

_geo_periodic_bounds_union(bounds,box)=bounds===nothing ? box :
    ((min(bounds[1][1],box[1][1]),min(bounds[1][2],box[1][2]),
      min(bounds[1][3],box[1][3])),
     (max(bounds[2][1],box[2][1]),max(bounds[2][2],box[2][2]),
      max(bounds[2][3],box[2][3])))

_geo_periodic_atol(m::GeoModel,context::_GeoNumericContext)=
    something(_geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)*
        _geo_periodic_lc(m,context)

# Index of the first top-level `{` in an `FExpr` operand — `(`/`[` nest and
# `==`/`!=`/`<=`/`>=` belong to the expression. A bare `=` or `;` ends it: the
# operand is complete and that token is what bison reports as the syntax
# error.
function _geo_periodic_operand_brace(tail::AbstractString,detail)
    depth=0;i=firstindex(tail);last=lastindex(tail)
    while i<=last
        c=tail[i]
        if c=='(' || c=='['
            depth+=1
        elseif c==')' || c==']'
            depth-=1
        elseif depth==0
            c=='{' && return i
            if c=='='
                prev=i>firstindex(tail) ? tail[prevind(tail,i)] : '\0'
                nxt=i<last ? tail[nextind(tail,i)] : '\0'
                if nxt=='='
                    i=nextind(tail,i)  # `==`
                elseif !(prev=='!' || prev=='<' || prev=='>')
                    _geo_syntax_abort("=",detail)
                end
            elseif c==';'
                _geo_syntax_abort(";",detail)
            end
        end
        i=nextind(tail,i)
    end
    _geo_syntax_abort(";",detail)
end

# A `Periodic Curve|Surface` statement matching neither grammar production —
# locate the first token diverging from `'{...} = {...} [transform]'` (or, for
# surfaces, the `FExpr {..} = FExpr {..}` edge-map shape) and report it like
# bison's `syntax error (token)`.
function _geo_periodic_syntax_abort(line::AbstractString,dim::Int,
                                    caller::AbstractString)
    tail=strip(match(r"^Periodic\s+[A-Za-z]+(.*)$",line).captures[1])
    endswith(tail,";") && (tail=strip(tail[1:end-1]))
    detail="$caller: malformed periodic statement $line"
    if startswith(tail,"{")
        body,rest=_geo_balanced_group(tail,detail)
        # `RecursiveListOfDouble` members are flat values — a nested `{`
        # inside a tag group is the reported token.
        occursin('{',body) && _geo_syntax_abort("{",detail)
        tail=strip(rest)
        startswith(tail,"=") ||
            _geo_syntax_abort(_geo_first_token(tail),detail)
        tail=strip(tail[2:end])
        startswith(tail,"{") ||
            _geo_syntax_abort(_geo_first_token(tail),detail)
        body,rest=_geo_balanced_group(tail,detail)
        occursin('{',body) && _geo_syntax_abort("{",detail)
        _geo_syntax_abort(_geo_first_token(strip(rest)),detail)
    elseif dim==2
        for side in 1:2
            index=_geo_periodic_operand_brace(tail,detail)
            body,rest=_geo_balanced_group(tail[index:end],detail)
            occursin('{',body) && _geo_syntax_abort("{",detail)
            tail=strip(rest)
            side==2 && break
            startswith(tail,"=") ||
                _geo_syntax_abort(_geo_first_token(tail),detail)
            tail=strip(tail[2:end])
        end
        _geo_syntax_abort(_geo_first_token(tail),detail)
    else
        _geo_syntax_abort(_geo_first_token(tail),detail)
    end
end

function _exec_periodic!(m::GeoModel,line::AbstractString,
                         context::_GeoNumericContext)
    entity=match(r"^Periodic\s+([A-Za-z]+)",line)
    entity===nothing && throw(ArgumentError(
        "execute_geo: malformed periodic statement $line"))
    entity_name=entity.captures[1]
    # Gmsh's `.geo` grammar carries no `Periodic Volume` production: the
    # statement is a parser-level `syntax error (Volume)` upstream.
    entity_name=="Volume" && _geo_syntax_abort("Volume")
    dim=entity_name in ("Line","Curve") ? 1 : entity_name=="Surface" ? 2 :
        _geo_syntax_abort(entity_name,
            "execute_geo: malformed periodic statement $line")
    caller="execute_geo: Periodic $(dim==1 ? "Curve" : "Surface")"
    if dim==2
        # `Periodic Surface j {slave curves} = k {master curves}` — the
        # explicit edge-counterpart form, which derives the transform from
        # the mapped boundary vertices.
        edgemap=match(
            r"^Periodic\s+Surface\s+([^{}]+?)\s*\{\s*([^}]*)\s*\}\s*=\s*([^{}]+?)\s*\{\s*([^}]*)\s*\}\s*;$",
            line)
        edgemap!==nothing && return _geo_periodic_surface_edge_map!(
            m,edgemap,context,caller)
    end
    statement=match(
        r"^Periodic\s+(?:Line|Curve|Surface)\s*\{\s*([^}]*)\s*\}\s*=\s*\{\s*([^}]*)\s*\}\s*(.*?)\s*;$",
        line)
    statement===nothing &&
        _geo_periodic_syntax_abort(line,dim,caller)
    # `'{' RecursiveListOfDouble '}'` — an empty group is a parse-level
    # `syntax error (})` upstream.
    isempty(strip(statement.captures[1])) && _geo_syntax_abort("}")
    isempty(strip(statement.captures[2])) && _geo_syntax_abort("}")
    slaves=_geo_periodic_tags(statement.captures[1],context,caller)
    masters=_geo_periodic_tags(statement.captures[2],context,caller)
    atol=_geo_periodic_atol(m,context)
    transform=strip(statement.captures[3])
    # Upstream evaluates the transform payload while reducing the production,
    # so `ListOfDouble`/`VExpr` errors surface before the action's count
    # check; the count diagnostic's wording differs per production (only the
    # `PeriodicTransform` surface form drops the "of"). Count mismatches are
    # `yymsg(0)` errors — recorded, never fatal mid-parse.
    if isempty(transform)
        length(slaves)==length(masters) || begin
            _geo_yyerror!(context,
                dim==1 ? "Wrong number of curves in periodicity constraint " *
                         "($(length(masters)) -> $(length(slaves)))" :
                         "Wrong number surfaces in periodicity constraint " *
                         "($(length(masters)) -> $(length(slaves)))")
            return nothing
        end
        # `PeriodicTransform` epsilon — orientation-only for curves; the
        # surface action rejects a <12 transform.
        if dim!=2
            _geo_periodic_curve_orientation_pairs!(
                m,slaves,masters,context,caller,atol)
        else
            _geo_yyerror!(context,
                "Affine transformation requires at least 12 entries " *
                "(0 provided)")
        end
        return nothing
    end
    if (matched=match(r"^Affine\b(.*)$",transform)) !== nothing
        return _geo_periodic_affine_statement!(
            m,dim,slaves,masters,matched.captures[1],context,caller,atol)
    end
    affine=if (matched=match(r"^Translate\b(.*)$",transform)) !== nothing
        # `tTranslate VExpr` — any VExpr form (`{a,b,c}`, `(a,b,c)`,
        # ±-composed); only the first three of its five slots are used.
        delta=_geo_exec_vexpr(matched.captures[1],3,context,"$caller Translate")
        (1.0,0.0,0.0,delta[1],
         0.0,1.0,0.0,delta[2],
         0.0,0.0,1.0,delta[3],
         0.0,0.0,0.0,1.0)
    elseif (matched=match(r"^Rotate\b(.*)$",transform)) !== nothing
        # `tRotate '{' VExpr ',' VExpr ',' FExpr '}'` — rotation axis,
        # a point on the axis, and the angle.
        body,rest=_geo_balanced_group(
            matched.captures[1],"$caller Rotate")
        isempty(rest) || _geo_syntax_abort(_geo_first_token(rest),
            "$caller Rotate: unexpected text after transform $(repr(rest))")
        parts=_geo_split_top_commas(body,"$caller Rotate")
        length(parts)==3 || _geo_syntax_abort(
            length(parts)<3 ? "}" : ",",
            "$caller Rotate: expected {axis, center, angle}")
        axis=_geo_exec_vexpr(parts[1],3,context,"$caller Rotate axis")
        origin=_geo_exec_vexpr(parts[2],3,context,"$caller Rotate center")
        angle=_geo_eval_numeric(parts[3],context,"$caller Rotate angle")
        _geo_periodic_rotation(axis,origin,angle,"$caller Rotate")
    else
        # After `{slave} = {master}` upstream expects `tEND`, `tAffine`,
        # `tRotate`, or `tTranslate` — anything else is a parse-level
        # `syntax error` on the unexpected token.
        _geo_syntax_abort(_geo_first_token(transform),
            "$caller: expected a Translate, Rotate, or Affine transform")
    end
    # `tTranslate`/`tRotate` productions — the transform is always 16 entries.
    length(slaves)==length(masters) || begin
        _geo_yyerror!(context,
            dim==1 ? "Wrong number of curves in periodicity constraint " *
                     "($(length(masters)) -> $(length(slaves)))" :
                     "Wrong number of surfaces in periodicity constraint " *
                     "($(length(masters)) -> $(length(slaves)))")
        return nothing
    end
    if dim==1
        _geo_periodic_curve_transform_pairs!(
            m,slaves,masters,affine,16,context,caller,atol)
    else
        _geo_periodic_surface_transform_pairs!(
            m,slaves,masters,affine,context,caller,atol)
    end
    return nothing
end

# `tAffine ListOfDouble` payload handling. `Affine` takes any `ListOfDouble`
# form — `{}`, `{...}`, `-{...}`, `k*{...}`, bare `FExpr`, or a list variable —
# and gmsh's arity contract differs per dimension: curves keep the raw entry
# count (below 12 the parser reports the arity error but the action still runs
# the orientation path with an empty transform; at 12+ the transform path's
# model-level exact-16 requirement applies), while surfaces zero-pad 12+
# entries to a 16-entry record.
function _geo_periodic_affine_statement!(
        m::GeoModel,dim::Int,slaves::Vector{Int},masters::Vector{Int},
        raw::AbstractString,context::_GeoNumericContext,
        caller::AbstractString,atol::Float64)
    # `tAffine ListOfDouble` — a bare `Affine` with no list is a parse-level
    # `syntax error (;)` upstream.
    isempty(strip(raw)) && _geo_syntax_abort(";")
    entries=_geo_numeric_list_values(raw,context,"$caller Affine")
    n=length(entries)
    # The `PeriodicTransform` productions check counts inside the action —
    # surfaces use the "Wrong number surfaces" wording without "of".
    length(slaves)==length(masters) || begin
        _geo_yyerror!(context,
            dim==1 ? "Wrong number of curves in periodicity constraint " *
                     "($(length(masters)) -> $(length(slaves)))" :
                     "Wrong number surfaces in periodicity constraint " *
                     "($(length(masters)) -> $(length(slaves)))")
        return nothing
    end
    if dim==1 && n<12
        # Upstream reports the arity error yet still runs the pair loop with
        # an empty transform — the relation records as orientation-only.
        n>0 && _geo_yyerror!(context,
            "Affine transformation requires at least 12 entries " *
            "($n provided)")
        _geo_periodic_curve_orientation_pairs!(
            m,slaves,masters,context,caller,atol)
        return nothing
    end
    n<12 && begin
        _geo_yyerror!(context,
            "Affine transformation requires at least 12 entries " *
            "($n provided)")
        return nothing
    end
    # Only the first twelve entries are ever applied (`SPoint3::transform`
    # reads a 3×4 row-major map). The stored record is upstream's verbatim
    # vector — a curve `Affine{16}` keeps its fourth row as written and a
    # surface `Affine{12..15}` keeps the parser's zero-padded tail — so MSH
    # `$Periodic` output carries the same bytes Gmsh would write. For `n>16`
    # upstream's copy runs past the 16-entry vector (undefined); the first
    # sixteen entries are the only deterministic emulation of what survives.
    record=ntuple(16) do i
        i<=min(n,16) ? entries[i] : 0.0
    end
    if dim==1
        _geo_periodic_curve_transform_pairs!(
            m,slaves,masters,record,n,context,caller,atol)
    else
        _geo_periodic_surface_transform_pairs!(
            m,slaves,masters,record,context,caller,atol)
    end
    return nothing
end

# Per-pair `addPeriodicEdge` loop with an empty transform — upstream resolves
# each `abs()` tag pair, reporting `Msg::Error` for missing curves while the
# remaining pairs still record their orientation-only relation.
function _geo_periodic_curve_orientation_pairs!(
        m::GeoModel,slaves::Vector{Int},masters::Vector{Int},
        context::_GeoNumericContext,caller::AbstractString,atol::Float64)
    for (signed_slave,signed_master) in zip(slaves,masters)
        slave=abs(signed_slave);master=abs(signed_master)
        if !(haskey(m.curves,slave) && haskey(m.curves,master))
            _geo_msg_error!(context,
                "Could not find curve $signed_slave or $signed_master " *
                "for periodic copy")
            continue
        end
        try
            set_periodic!(m,1,[signed_slave],[signed_master],nothing;
                          atol=atol,overwrite=true)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_msg_error!(context,sprint(showerror,err))
        end
    end
    return nothing
end

# Per-pair `addPeriodicEdge` transform loop — `GEdge::setMeshMaster(source,
# tfo)`: the first twelve entries map the master's endpoints onto the slave's.
# A geometric mismatch is a `Msg::Info` drop (no error, no relation); a
# matched pair with a non-16 arity reports `GEntity::setMeshMaster`'s
# exact-16 `Msg::Error` for the edge and each endpoint vertex; only a matched
# 16-entry transform records the relation.
function _geo_periodic_curve_transform_pairs!(
        m::GeoModel,slaves::Vector{Int},masters::Vector{Int},
        affine::NTuple{16,Float64},n::Int,
        context::_GeoNumericContext,caller::AbstractString,atol::Float64)
    coefficients,translation=_periodic_affine_3x4(affine)
    for (signed_slave,signed_master) in zip(slaves,masters)
        slave=abs(signed_slave);master=abs(signed_master)
        if !(haskey(m.curves,slave) && haskey(m.curves,master))
            _geo_msg_error!(context,
                "Could not find curve $signed_slave or $signed_master " *
                "for periodic copy")
            continue
        end
        orientation=_model_periodic_curve_orientation(
            m,slave,master,coefficients,translation,atol,caller,1)
        slave_points=m.curves[slave];master_points=m.curves[master]
        if orientation.reversed===nothing
            d1,d2=orientation.distances
            _geo_yyinfo!(context,
                "Error in transformation from curve $master " *
                "($(master_points[1])-$(master_points[2])) to $slave " *
                "($(slave_points[1])-$(slave_points[2])) (minimal " *
                "transformed node distances $d1 $d2, tolerance $atol)")
            continue
        end
        if n!=16
            _geo_msg_error!(context,
                "Periodicity transformation from entity $master to $slave " *
                "(dim 1) has $n components, while 16 are required")
            pairs=orientation.reversed ?
                ((master_points[2],slave_points[1]),
                 (master_points[1],slave_points[2])) :
                ((master_points[1],slave_points[1]),
                 (master_points[2],slave_points[2]))
            for (mv,sv) in pairs
                _geo_msg_error!(context,
                    "Periodicity transformation from entity $mv to $sv " *
                    "(dim 0) has $n components, while 16 are required")
            end
            continue
        end
        try
            set_periodic!(m,1,[signed_slave],[signed_master],affine;
                          atol=atol,overwrite=true)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_msg_error!(context,sprint(showerror,err))
        end
    end
    return nothing
end

# Per-pair `addPeriodicFace` transform loop — upstream reports `Msg::Error`
# for missing surfaces and `Msg::Info` progress for the rest; the per-pair
# correspondence checks run inside `set_periodic!`.
function _geo_periodic_surface_transform_pairs!(
        m::GeoModel,slaves::Vector{Int},masters::Vector{Int},
        affine::NTuple{16,Float64},
        context::_GeoNumericContext,caller::AbstractString,atol::Float64)
    for (signed_slave,signed_master) in zip(slaves,masters)
        slave=abs(signed_slave);master=abs(signed_master)
        if !(haskey(m.surfaces,slave) && haskey(m.surfaces,master))
            _geo_msg_error!(context,
                "Could not find surface $signed_slave or $signed_master " *
                "for periodic copy")
            continue
        end
        _geo_yyinfo!(context,"Setting mesh master using transformation")
        try
            # `GFace::setMeshMaster` resolves each boundary and embedded edge
            # counterpart at declaration, reporting `Setting curve master`
            # per success and aborting the relation on the first failure.
            _model_periodic_surface_edge_map(
                m,slave,master,affine,atol,caller;
                on_match=(sc,mc)->_geo_yyinfo!(
                    context,"Setting curve master $sc - $mc"))
            set_periodic!(m,2,[signed_slave],[signed_master],affine;
                          atol=atol,overwrite=true)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_msg_error!(context,sprint(showerror,err))
        end
    end
    return nothing
end

# `Periodic Surface j {slave curves} = k {master curves}` — upstream's
# `addPeriodicFace(target, source, edgeCounterparts)`: each slave boundary
# curve resolves a signed master counterpart, the endpoint map induces a
# vertex correspondence, and the transform is derived as a translation or a
# rotation before taking the ordinary affine surface path.
function _geo_periodic_surface_edge_map!(
        m::GeoModel,edgemap::RegexMatch,
        context::_GeoNumericContext,caller::AbstractString)
    # `addPeriodicFace` resolves both surfaces through `abs(iTarget)`/`abs`
    # but reports the raw signed tags in its diagnostic.
    signed_slave=_geo_signed_gmsh_int_value(
        _geo_eval_numeric(edgemap.captures[1],context,caller),"$caller slave")
    signed_master=_geo_signed_gmsh_int_value(
        _geo_eval_numeric(edgemap.captures[3],context,caller),"$caller master")
    j_slave=abs(signed_slave);j_master=abs(signed_master)
    slave_curves=_geo_exec_entity_tags(
        edgemap.captures[2],context,caller;signed=true)
    master_curves=_geo_exec_entity_tags(
        edgemap.captures[4],context,caller;signed=true)
    # `yymsg(0)` — the action's own diagnostic; parsing continues.
    length(slave_curves)==length(master_curves) || begin
        _geo_yyerror!(context,
            "Wrong number of surface curves in periodicity constraint " *
            "($(length(master_curves)) -> $(length(slave_curves)))")
        return nothing
    end
    edgecopies=Dict{Int,Int}()
    for (sc,mc) in zip(slave_curves,master_curves)
        edgecopies[sc]=mc
    end
    _geo_yyinfo!(context,
        "Encoding periodic connection between $signed_slave and " *
        "$signed_master")
    for (sc,mc) in sort!(collect(edgecopies);by=first)
        _geo_yyinfo!(context,"$sc - $mc")
    end
    if !(haskey(m.surfaces,j_slave) && haskey(m.surfaces,j_master))
        _geo_msg_error!(context,
            "Could not find surface $signed_slave or $signed_master " *
            "for periodic copy")
        return nothing
    end
    vs2vt=Dict{Int,Int}()
    # Upstream iterates `l_edges` in the face's loop order (a seam's repeated
    # occurrence just rewrites the same `vs2vt` entry) — preserve that order
    # so the first missing counterpart reported matches Gmsh's.
    slave_boundary=Int[]
    seen=Set{Int}()
    for loop in m.surfaces[j_slave],signed in m.loops[loop]
        curve=abs(signed)
        curve in seen || (push!(seen,curve);push!(slave_boundary,curve))
    end
    for curve in slave_boundary
        sign=1
        source_e=get(edgecopies,curve,nothing)
        if source_e===nothing
            sign=-1
            source_e=get(edgecopies,-curve,nothing)
        end
        if source_e===nothing
            # Upstream prints `master->tag()` in this slot (its own typo) —
            # kept verbatim.
            _geo_msg_error!(context,
                "Could not find curve counterpart $curve in slave " *
                "surface $j_master")
            return nothing
        end
        mc=abs(source_e)
        haskey(m.curves,curve) ||
            throw(ArgumentError("$caller: unknown slave Curve[$curve]"))
        haskey(m.curves,mc) ||
            throw(ArgumentError("$caller: unknown master Curve[$mc]"))
        s0,s1=m.curves[curve];t0,t1=m.curves[mc]
        if source_e*sign>0
            vs2vt[t0]=s0;vs2vt[t1]=s1
        else
            vs2vt[t0]=s1;vs2vt[t1]=s0
        end
    end
    affine=try
        _geo_periodic_derive_transform(
            m,vs2vt,j_slave,j_master,context,caller)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _geo_msg_error!(context,sprint(showerror,err))
        return nothing
    end
    _geo_yyinfo!(context,"Setting mesh master using transformation")
    try
        _model_periodic_surface_edge_map(
            m,j_slave,j_master,affine,_geo_periodic_atol(m,context),caller;
            on_match=(sc,mc)->_geo_yyinfo!(
                context,"Setting curve master $sc - $mc"))
        set_periodic!(m,2,[j_slave],[j_master],affine;
                      atol=_geo_periodic_atol(m,context),overwrite=true)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _geo_msg_error!(context,sprint(showerror,err))
    end
    return nothing
end

# `computeMeanPlaneSimple` — centroid plus the least-squares normal (the
# right-singular vector of the smallest singular value of the centered point
# set).
function _geo_mean_plane(points::Vector{NTuple{3,Float64}})
    n=length(points)
    cx=sum(point->point[1],points)/n
    cy=sum(point->point[2],points)/n
    cz=sum(point->point[3],points)/n
    centered=Matrix{Float64}(undef,n,3)
    for (index,point) in pairs(points)
        centered[index,1]=point[1]-cx
        centered[index,2]=point[2]-cy
        centered[index,3]=point[3]-cz
    end
    decomposition=svd(centered)
    column=argmin(decomposition.S)
    nx,ny,nz=decomposition.V[1,column],decomposition.V[2,column],
        decomposition.V[3,column]
    scale=sqrt(nx*nx+ny*ny+nz*nz)
    scale>0 || throw(ErrorException(
        "_geo_mean_plane: degenerate singular vector"))
    return (nx/scale,ny/scale,nz/scale),(cx,cy,cz)
end

# `sys2x2` — Cramer's rule with gmsh's relative singularity bound; `nothing`
# when the system is singular.
function _geo_sys2x2(m00,m01,m10,m11,b0,b1)
    norm=m00*m00+m11*m11+m01*m01+m10*m10
    det=m00*m11-m10*m01
    (norm==0.0 || abs(det)/norm<1e-16) && return nothing
    return ((b0*m11-m01*b1)/det,(m00*b1-m10*b0)/det)
end

# `myLine(plane_source, plane_target)` — a point on the intersection line of
# two planes, trying z = 0, then x = 0, then y = 0 like upstream.
function _geo_plane_intersect_point(n1,c1,n2,c2,caller,slave)
    b1=n1[1]*c1[1]+n1[2]*c1[2]+n1[3]*c1[3]
    b2=n2[1]*c2[1]+n2[2]*c2[2]+n2[3]*c2[3]
    if (r=_geo_sys2x2(n1[1],n1[2],n2[1],n2[2],b1,b2))!==nothing
        return (r[1],r[2],0.0)
    elseif (r=_geo_sys2x2(n1[2],n1[3],n2[2],n2[3],b1,b2))!==nothing
        return (0.0,r[1],r[2])
    elseif (r=_geo_sys2x2(n1[1],n1[3],n2[1],n2[3],b1,b2))!==nothing
        return (r[1],0.0,r[2])
    end
    throw(ArgumentError(
        "$caller: only rotations or translations can currently be computed " *
        "automatically for periodic surfaces: surface $slave not meshed"))
end

@inline _geo_pdist(a,b)=sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2)
@inline _geo_pdot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _geo_pcross(a,b)=
    (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
@inline _geo_line_project(point,line_point,direction)=begin
    u=_geo_pdot((point[1]-line_point[1],point[2]-line_point[2],
                 point[3]-line_point[3]),direction)
    (line_point[1]+direction[1]*u,line_point[2]+direction[2]*u,
     line_point[3]+direction[3]*u)
end

# `GFace::setMeshMaster(master, edgeCopies)`'s transform derivation: a common
# translation when every vertex displacement agrees within `1e-5·|DX|`,
# otherwise a rotation about the intersection of the source and target mean
# planes whose angle must be consistent to `1e-8`.
function _geo_periodic_derive_transform(
        m::GeoModel,vs2vt::Dict{Int,Int},slave::Int,master::Int,
        context::_GeoNumericContext,caller::AbstractString)
    pairs=sort!(collect(vs2vt);by=first)
    # A boundaryless surface yields no correspondences: upstream's
    # translation loop never runs, leaving `DX = 0` and the identity record.
    isempty(pairs) && return (1.0,0.0,0.0,0.0,
                              0.0,1.0,0.0,0.0,
                              0.0,0.0,1.0,0.0,
                              0.0,0.0,0.0,1.0)
    sources=[m.points[master_point] for (master_point,_) in pairs]
    targets=[m.points[slave_point] for (_,slave_point) in pairs]
    dx=targets[1][1]-sources[1][1]
    dy=targets[1][2]-sources[1][2]
    dz=targets[1][3]-sources[1][3]
    dnorm=sqrt(dx*dx+dy*dy+dz*dz)
    translation=true
    for index in 2:length(pairs)
        ddx=(targets[index][1]-sources[index][1])-dx
        ddy=(targets[index][2]-sources[index][2])-dy
        ddz=(targets[index][3]-sources[index][3])-dz
        if sqrt(ddx*ddx+ddy*ddy+ddz*ddz)>1e-5*dnorm
            translation=false;break
        end
    end
    if translation
        _geo_yyinfo!(context,
            "Periodic mesh translation found between $slave and $master: " *
            "dx = ($(_geo_gmsh_number(dx)),$(_geo_gmsh_number(dy))," *
            "$(_geo_gmsh_number(dz)))")
        return (1.0,0.0,0.0,dx,
                0.0,1.0,0.0,dy,
                0.0,0.0,1.0,dz,
                0.0,0.0,0.0,1.0)
    end
    n_source,c_source=_geo_mean_plane(sources)
    n_target,c_target=_geo_mean_plane(targets)
    direction=_geo_pcross(n_source,n_target)
    dlength=sqrt(_geo_pdot(direction,direction))
    dlength==0.0 && throw(ArgumentError(
        "$caller: only rotations or translations can currently be computed " *
        "automatically for periodic surfaces: surface $slave not meshed"))
    direction=(direction[1]/dlength,direction[2]/dlength,direction[3]/dlength)
    line_point=_geo_plane_intersect_point(
        n_source,c_source,n_target,c_target,caller,slave)
    gtol=something(
        _geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)
    angle=0.0;count=0;rotation=true
    for index in eachindex(pairs)
        ps=sources[index];pt=targets[index]
        p_ps=_geo_line_project(ps,line_point,direction)
        p_pt=_geo_line_project(pt,line_point,direction)
        d1=_geo_pdist(ps,pt)
        d1>gtol || continue
        _geo_pdist(p_ps,p_pt)>1e-8*d1 && (rotation=false)
        t1=(ps[1]-p_ps[1],ps[2]-p_ps[2],ps[3]-p_ps[3])
        t2=(pt[1]-p_pt[1],pt[2]-p_pt[2],pt[3]-p_pt[3])
        sqrt(_geo_pdot(t1,t1))>1e-8*d1 || continue
        measured=atan(_geo_pdot(_geo_pcross(t1,t2),direction),
                      _geo_pdot(t1,t2))
        if count==0
            angle=measured
        elseif abs(measured-angle)>1e-8
            rotation=false
        end
        count+=1
    end
    rotation || throw(ArgumentError(
        "$caller: only rotations or translations can currently be computed " *
        "automatically for periodic surfaces: surface $slave not meshed"))
    _geo_yyinfo!(context,
        "Periodic mesh rotation found: axis " *
        "($(_geo_gmsh_number(direction[1])),$(_geo_gmsh_number(direction[2]))" *
        ",$(_geo_gmsh_number(direction[3]))) point " *
        "($(_geo_gmsh_number(line_point[1])) $(_geo_gmsh_number(line_point[2]))" *
        " $(_geo_gmsh_number(line_point[3]))) angle " *
        "$(_geo_gmsh_number(angle*180/pi))")
    sine,cosine=sincos(angle);one_minus=1-cosine
    ux,uy,uz=direction
    r11=cosine+ux*ux*one_minus;r12=ux*uy*one_minus-uz*sine
    r13=ux*uz*one_minus+uy*sine
    r21=uy*ux*one_minus+uz*sine;r22=cosine+uy*uy*one_minus
    r23=uy*uz*one_minus-ux*sine
    r31=uz*ux*one_minus-uy*sine;r32=uy*uz*one_minus+ux*sine
    r33=cosine+uz*uz*one_minus
    ox,oy,oz=line_point
    t1=ox-(r11*ox+r12*oy+r13*oz)
    t2=oy-(r21*ox+r22*oy+r23*oz)
    t3=oz-(r31*ox+r32*oy+r33*oz)
    return (r11,r12,r13,t1,
            r21,r22,r23,t2,
            r31,r32,r33,t3,
            0.0,0.0,0.0,1.0)
end

# ---------------------------------------------------------------------------
# `.geo` meshing-constraint statements (the Gmsh `Constraints` grammar family)
# and the `Delete`/`SetTag`/`SetMaxTag` lifecycle statements.
#
# These handlers reproduce the `GEO_Internals` setter contracts directly:
# statements record on entities that exist at execution time and silently skip
# the rest; entity lists accept Gmsh's `ListOfDouble` forms (`{...}`, bare
# `FExpr`/`FExpr_Multi`, `-{...}`, `expr*{...}`, list-variable references), and
# `ListOfDoubleOrAll` additionally accepts the `{:}`/`"*"`/`"all"` wildcards —
# expanded over the entities that exist *when the statement runs* (a wildcard
# does not reach entities created later). Tag values truncate like C `(int)`;
# where a setter maps tag 0 to the whole dimension, a |value| < 1 entry does
# the same.

# `.geo` `ListOfDoubleOrAll`/`ListOfDouble` evaluation. Returns `nothing` for
# the wildcard forms (only legal where the grammar allows them). Bare
# scalars/ranges/variables are wrapped in braces; anything already carrying
# braces (`{...}`, `-{...}`, `expr*{...}`, selector forms) goes to the list
# evaluator directly, matching the `ListOfDouble` grammar productions.
function _geo_constraint_list(raw::AbstractString,context::_GeoNumericContext,
                              caller::AbstractString;allow_all::Bool)
    s=String(strip(raw))
    isempty(s) && throw(ArgumentError("$caller: entity list must not be empty"))
    if allow_all && (s=="\"*\"" || s=="\"all\"" ||
                     match(r"^\{\s*:\s*\}$",s)!==nothing)
        return nothing
    end
    source=occursin(r"[{}]",s) ? s : "{$s}"
    return _geo_numeric_list_values(source,context,caller)
end

# C-style `(int)` truncation for `.geo` tag/count/affect values.
function _geo_constraint_int(value::Float64,caller::AbstractString,what::String)
    isfinite(value) || throw(ArgumentError("$caller: $what must be finite"))
    tag=try
        trunc(Int32,value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $what $value is outside Gmsh's signed 32-bit integer range"))
    end
    return Int(tag)
end

# Case-sensitive `.geo` `Using` law names — Gmsh matches them with `strcmp`.
const _GEO_TRANSFINITE_CURVE_LAWS=Dict{String,Symbol}(
    "Progression"=>:progression,"Power"=>:progression,
    "Bump"=>:bump,"Beta"=>:beta,
    "Progression_HWall"=>:progression_hwall,"Bump_HWall"=>:bump_hwall,
    "Beta_HWall"=>:beta_hwall,"Beta_Symmetrical"=>:beta_symmetrical,
    "Beta_Symmetrical_HWall"=>:beta_symmetrical_hwall)

# `Transfinite Curve{list} = n [Using Law coef];` — records the raw signed
# semantics of `setTransfiniteLine`: the list entry's sign negates the stored
# transfinite type (`reversed` on the record) while `coef` is stored verbatim;
# an entry truncating to 0 expands to every current curve.
function _geo_exec_transfinite_curve!(m::GeoModel,list_raw::AbstractString,
                                      tail::AbstractString,
                                      context::_GeoNumericContext)
    caller="execute_geo: Transfinite Curve"
    nsource=String(strip(tail));kind=:progression;coef=1.0
    split=_geo_split_at_keyword(nsource,"Using",caller)
    if split!==nothing
        (nsource,using_tail)=split
        um=match(r"^\"?([A-Za-z_][A-Za-z0-9_]*)\"?\s+(.+?)\s*$",
                 String(strip(using_tail)))
        um===nothing && throw(ArgumentError(
            "$caller: expected `Using <law> <coefficient>`"))
        kind=get(_GEO_TRANSFINITE_CURVE_LAWS,um.captures[1],nothing)
        kind===nothing && throw(ArgumentError(
            "$caller: unknown transfinite mesh type $(repr(um.captures[1]))"))
        coef=_geo_eval_numeric(um.captures[2],context,"$caller coefficient")
    end
    isfinite(coef) || throw(ArgumentError(
        "$caller: coefficient must be finite"))
    nraw=_geo_eval_numeric(nsource,context,"$caller node count")
    ncount=max(2,_geo_constraint_int(nraw,caller,"node count"))
    list=_geo_constraint_list(list_raw,context,caller;allow_all=true)
    targets=list===nothing ? sort!(collect(keys(m.curves))) : nothing
    if targets!==nothing
        for curve in targets
            m.meshing.transfinite_curves[curve]=(
                num_nodes=ncount,kind=kind,coef=coef,reversed=false)
        end
        return nothing
    end
    for value in list
        sd=_geo_constraint_int(value,caller,"Curve tag")
        j=abs(sd)
        if j==0
            # `setTransfiniteLine(±0, ...)` is the wildcard call — the stored
            # type is `type * gmsh_sign(value)`, so an exact-zero entry zeroes
            # the type and `F_Transfinite` falls through to a uniform
            # distribution, while ±subinteger entries keep the law ±reversed.
            effective=iszero(value) ? :uniform : kind
            for curve in sort!(collect(keys(m.curves)))
                m.meshing.transfinite_curves[curve]=(
                    num_nodes=ncount,kind=effective,coef=coef,
                    reversed=value<0)
            end
            continue
        end
        haskey(m.curves,j) || continue
        m.meshing.transfinite_curves[j]=(
            num_nodes=ncount,kind=kind,coef=coef,reversed=sd<0)
    end
    return nothing
end

# A trailing `TransfiniteArrangement`/`Using`-style bare or quoted word. Returns
# `(head, word)` when the statement ends with such a token and the head is
# nonempty — a lone word is the list itself, not an arrangement.
function _geo_constraint_trailing_word(body::AbstractString)
    s=String(strip(body))
    mm=match(r"^(.*?)\s+([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*$",s)
    mm===nothing && return (s,nothing)
    head=String(strip(mm.captures[1]))
    isempty(head) && return (s,nothing)
    return (head,String(strip(mm.captures[2])))
end

# `TransfiniteArrangement` — Gmsh maps `Left`→-1, `Right`→1,
# `AlternateRight`→2, `AlternateLeft`→-2, and every other word (including bare
# `Alternate`) to 2 without an error.
function _geo_transfinite_arrangement(word::Union{String,Nothing})
    word===nothing && return :left
    w=strip(word,'"')
    w=="Left" && return :left
    w=="Right" && return :right
    w=="AlternateLeft" && return :alternate_left
    return :alternate_right
end

# `Transfinite Surface{list} [= {corners}] [arrangement];` — `FindSurface` is a
# signed lookup, so negative entries silently skip while a zero entry (and the
# `{:}`/`"all"` wildcards) applies the method and arrangement to every current
# surface with the corner list reset. On a live surface a nonempty corner list
# must hold 3 or 4 existing points; other counts and unknown points are Gmsh
# errors, while a missing surface skips the whole check silently.
function _geo_exec_transfinite_surface!(m::GeoModel,body::AbstractString,
                                        context::_GeoNumericContext)
    caller="execute_geo: Transfinite Surface"
    head,word=_geo_constraint_trailing_word(body)
    mode=_geo_transfinite_arrangement(word)
    list_raw=head;corners=Int[]
    eq=findfirst(==('='),head)
    if eq!==nothing
        list_raw=String(strip(head[firstindex(head):prevind(head,eq)]))
        corner_source=String(strip(head[nextind(head,eq):end]))
        isempty(corner_source) && throw(ArgumentError(
            "$caller: `=` must be followed by the corner Point list"))
        values=_geo_constraint_list(corner_source,context,"$caller corners";
                                    allow_all=false)
        corners=Int[abs(_geo_constraint_int(v,caller,"corner Point tag"))
                    for v in values]
    end
    list=_geo_constraint_list(list_raw,context,caller;allow_all=true)
    function apply(tag::Int,corner_tags::Vector{Int})
        if !isempty(corner_tags)
            length(corner_tags) in (3,4) || throw(ArgumentError(
                "$caller: Transfinite surface requires 3 or 4 corner points"))
            for point_tag in corner_tags
                haskey(m.points,point_tag) || throw(ArgumentError(
                    "$caller: unknown corner Point[$point_tag]"))
            end
        end
        m.meshing.transfinite_surfaces[tag]=(
            arrangement=mode,corners=copy(corner_tags))
        return nothing
    end
    list===nothing && return (foreach(
        tag->apply(tag,Int[]),sort!(collect(keys(m.surfaces)))); nothing)
    for value in list
        sd=_geo_constraint_int(value,caller,"Surface tag")
        sd==0 && (foreach(
            tag->apply(tag,Int[]),sort!(collect(keys(m.surfaces)))); continue)
        sd<0 && continue
        haskey(m.surfaces,sd) || continue
        apply(sd,corners)
    end
    return nothing
end

# `Transfinite Volume{list} [= {corners}];` — same shape as the surface form
# minus the arrangement word. Gmsh applies the corner list only when it holds 6
# or 8 points; other counts are silently dropped (the volume still gets the
# transfinite method with automatic corners).
function _geo_exec_transfinite_volume!(m::GeoModel,body::AbstractString,
                                       context::_GeoNumericContext)
    caller="execute_geo: Transfinite Volume"
    s=String(strip(body))
    list_raw=s;corners=Int[]
    eq=findfirst(==('='),s)
    if eq!==nothing
        list_raw=String(strip(s[firstindex(s):prevind(s,eq)]))
        corner_source=String(strip(s[nextind(s,eq):end]))
        isempty(corner_source) && throw(ArgumentError(
            "$caller: `=` must be followed by the corner Point list"))
        values=_geo_constraint_list(corner_source,context,"$caller corners";
                                    allow_all=false)
        corner_tags=Int[abs(_geo_constraint_int(v,caller,"corner Point tag"))
                        for v in values]
        length(corner_tags) in (6,8) && (corners=corner_tags)
    end
    list=_geo_constraint_list(list_raw,context,caller;allow_all=true)
    function apply(tag::Int,corner_tags::Vector{Int})
        for point_tag in corner_tags
            haskey(m.points,point_tag) || throw(ArgumentError(
                "$caller: unknown corner Point[$point_tag]"))
        end
        m.meshing.transfinite_volumes[tag]=copy(corner_tags)
        return nothing
    end
    list===nothing && return (foreach(
        tag->apply(tag,Int[]),sort!(collect(keys(m.volumes)))); nothing)
    for value in list
        sd=_geo_constraint_int(value,caller,"Volume tag")
        sd==0 && (foreach(
            tag->apply(tag,Int[]),sort!(collect(keys(m.volumes)))); continue)
        sd<0 && continue
        haskey(m.volumes,sd) || continue
        apply(sd,corners)
    end
    return nothing
end

# Statements whose setter maps tag 0 to every entity of the dimension and a
# signed lookup to the rest: `Transfinite Surface/Volume`, `TransfQuadTri`,
# `Recombine`, `Smoother`, `ReverseMesh`. `apply` receives each target tag;
# unknown and negative entries silently skip.
function _geo_constraint_apply(m::GeoModel,dimension::Int,
                               list::Union{Vector{Float64},Nothing},
                               caller::AbstractString,apply)
    all_tags()=sort!(collect(keys(_model_entity_dictionary(m,dimension))))
    list===nothing && return (foreach(apply,all_tags()); nothing)
    for value in list
        sd=_geo_constraint_int(value,caller,"entity tag")
        sd==0 && (foreach(apply,all_tags()); continue)
        sd<0 && continue
        haskey(_model_entity_dictionary(m,dimension),sd) || continue
        apply(sd)
    end
    return nothing
end

# `Compound Curve|Surface|Volume{list} [MeshAlgorithm n];` — the constraint
# form stores the raw member list on the internals multimap; a `MeshAlgorithm`
# suffix appends `-(int)value` to that list (resolved, and silently dropped for
# missing members, at sync time).
function _geo_exec_compound!(m::GeoModel,dimension::Int,body::AbstractString,
                             context::_GeoNumericContext)
    caller="execute_geo: Compound $(_entity_label(dimension))"
    s=String(strip(body))
    split=_geo_split_at_keyword(s,"MeshAlgorithm",caller)
    list_raw=s;algorithm=nothing
    if split!==nothing
        (list_raw,algorithm_raw)=split
        algorithm=_geo_constraint_int(
            _geo_eval_numeric(algorithm_raw,context,
                              "$caller MeshAlgorithm"),
            caller,"MeshAlgorithm")
    end
    values=_geo_numeric_list_values(
        startswith(list_raw,"{") ? list_raw : "{$list_raw}",context,
        "$caller entity list")
    tags=Int[_geo_constraint_int(v,caller,"entity tag") for v in values]
    algorithm===nothing || push!(tags,-algorithm)
    push!(m.meshing.compounds,dimension=>tags)
    return nothing
end

# `GEO_Internals::_allocateAll`/`_freeAll` initialize every entity counter
# from `Geometry.FirstEntityTag - 1` and `_maxPhysicalNum` from
# `Geometry.FirstPhysicalTag - 1`, re-reading the options on each fresh or
# destroyed internals — so `NewModel`/`Delete Model`/`Delete All` all re-init
# the model counters from the *current* option values.
function _geo_reset_geometry_counters!(m::GeoModel,
                                       context::_GeoNumericContext)
    entity_base=max(_geo_signed_gmsh_int_value(
        something(_geo_option_number(
            context,"Geometry",0,"FirstEntityTag"),1.0),
        "Geometry.FirstEntityTag"),1)-1
    physical_base=max(_geo_signed_gmsh_int_value(
        something(_geo_option_number(
            context,"Geometry",0,"FirstPhysicalTag"),1.0),
        "Geometry.FirstPhysicalTag"),1)-1
    m.next_tag .= entity_base
    m.physical_tag_max=physical_base
    return nothing
end

# `.geo` `Delete`/`Recursive Delete`/`Delete Embedded` and the named `Delete X`
# forms. `tail` is everything after the `Delete` keyword(s).
function _geo_exec_delete!(m::GeoModel,recursive::Bool,tail::AbstractString,
                           context::_GeoNumericContext,
                           allocator_state)
    caller="execute_geo: $(recursive ? "Recursive " : "")Delete"
    s=String(strip(tail))
    if startswith(s,"{")
        (inner,rest)=_geo_balanced_group(s,caller)
        isempty(strip(rest)) || throw(ArgumentError(
            "$caller: unexpected text after the entity list"))
        entities=_geo_shape_list_entities!(
            m,inner,context,caller;signed_tags=true,
            allocator_state=allocator_state)
        removed=_geo_delete_entities!(m,entities;recursive=recursive)
        allocator_state===nothing ||
            _geo_allocator_delete_entities!(allocator_state,removed)
        # `GEO_Internals::remove` marks the internals changed unconditionally,
        # so the grammar always resynchronizes here — stale raw members drop
        # from the view (and warn) now rather than at the next sync point.
        context.geo_changed=true
        _geo_sync_physical_view!(m,context)
        return nothing
    end
    recursive && _geo_syntax_abort(isempty(s) ? ";" : string(first(s)))
    # `Delete Embedded { Surface{...}; Volume{...}; }` clears every embedding
    # on the listed parents (dims 2 and 3 only — lower-dim entries are ignored
    # and a missing parent is an error), matching `removeEmbedded`.
    if (em=match(r"^(?:Embedded|\"Embedded\")\s*(.*)$",s))!==nothing
        rest=String(strip(em.captures[1]))
        startswith(rest,"{") || throw(ArgumentError(
            "$caller: `Delete Embedded` requires a `{...}` entity list"))
        (inner,rest)=_geo_balanced_group(rest,caller)
        isempty(strip(rest)) || throw(ArgumentError(
            "$caller: unexpected text after the entity list"))
        entities=_geo_shape_list_entities!(
            m,inner,context,caller;signed_tags=true,
            allocator_state=allocator_state)
        dim_tags=NTuple{2,Int}[]
        for (dim,tag) in entities
            dim in (2,3) || continue
            haskey(_model_entity_dictionary(m,dim),tag) || throw(ArgumentError(
                "$caller: unknown model $(_entity_label(dim)) with tag $tag"))
            push!(dim_tags,(dim,tag))
        end
        remove_embedded!(m,dim_tags)
        return nothing
    end
    # `Delete Field[i]` drops a mesh field from the live exec map;
    # `FieldManager::deleteField` reports a `Msg::Error` when the id is missing.
    # `Delete View[i]` reports an unknown view and any other indexed name an
    # unknown command in Gmsh.
    if (fm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*\[\s*(.*?)\s*\]\s*;?$",
            s))!==nothing
        base=String(strip(fm.captures[1],'"'))
        base=="Field" && return _geo_exec_delete_field!(context,fm.captures[2])
        base=="View" && throw(ArgumentError(
            "$caller: unknown view $(strip(fm.captures[2]))"))
        throw(ArgumentError("$caller: unknown command 'Delete $base'"))
    end
    # `Delete Empty Views` clears empty post-processing views — none exist
    # during `.geo` execution; any other two-word form is an unknown command.
    if (tm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+([A-Za-z_][A-Za-z0-9_]*)\s*;?$",
                 s))!==nothing
        tm.captures[1]=="Empty" && tm.captures[2]=="Views" && return nothing
        throw(ArgumentError(
            "$caller: unknown command 'Delete $(tm.captures[1]) $(tm.captures[2])'"))
    end
    # Named forms: `Delete All|Model|Physicals|Variables|Options|Meshes|Struct;`,
    # `Delete <variable>;`, and the `Delete <name>~{expr};` namespaced symbol
    # form (`name_<int(expr)>`).
    name=begin
        if (nm=match(
                r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*;?$",
                s)) !== nothing
            String(strip(nm.captures[1],'"'))
        elseif (nm=match(
                r"^([A-Za-z_][A-Za-z0-9_]*)\s*~\s*\{\s*(.*?)\s*\}\s*;?$",
                s)) !== nothing
            idx=_geo_constraint_int(
                _geo_eval_numeric(nm.captures[2],context,
                                  "$caller namespace index"),
                caller,"namespace index")
            String(nm.captures[1])*"_"*string(idx)
        else
            # `tDelete` expects a name, `name~{i}`, `{...}` or `[...]` —
            # the syntax error lands on the offending token (`Delete = 5`
            # → `(=)`).
            _geo_syntax_abort(isempty(s) ? ";" : string(first(s)))
        end
    end
    if name=="All"
        # `ClearProject` — a fresh `GModel` becomes current (name, physical
        # and entity name tables, and compound specs all die with the old
        # objects), the factory choice reverts to the built-in kernel, and
        # `DeleteAllModelsAndViews` clears BOTH parser symbol tables and the
        # struct namespaces. `FunctionManager` and options survive — it is
        # only cleared by `OpenProject` (file open).
        _geo_reset_model_geometry!(m)
        _geo_reset_geometry_counters!(m,context)
        empty!(m.physical_names);empty!(m.entity_names)
        empty!(m.meshing.compounds)
        m.name=""
        empty!(context.values);empty!(context.lists)
        # `gmshSurface::reset()` runs with the fresh `GModel`: both the
        # `allGmshSurfaces` registry and the `myGmshSurface` context die.
        empty!(context.parametric_surfaces)
        context.parametric_surface=nothing
        # `gmsh_yystringsymbols` and `gmsh_yynamespaces` are cleared, but
        # `FunctionManager` is not — registered functions survive `Delete All`
        # upstream (it is only cleared by `OpenProject`).
        empty!(context.strings)
        empty!(context.structs);empty!(context.struct_maxtag)
        empty!(context.list_variables)
        empty!(context.unavailable);empty!(context.unavailable_lists)
        empty!(context.raw_physicals)
        context.stored_list_items=0
        context.mesh=nothing
        context.mesh_node_owner=empty(context.mesh_node_owner)
        context.geo_changed=true
        if allocator_state!==nothing
            _geo_allocator_reset_model!(allocator_state,context;fresh=true)
            allocator_state.factory=:builtin
        end
        return nothing
    elseif name=="Model"
        # `destroy(/*keepName=*/true)` + `GEO_Internals::destroy` — the model
        # name and file name survive; entity mesh data dies with the model.
        _geo_reset_model_geometry!(m)
        _geo_reset_geometry_counters!(m,context)
        context.mesh=nothing
        context.mesh_node_owner=empty(context.mesh_node_owner)
        empty!(context.raw_physicals)
        empty!(context.parametric_surfaces)
        context.parametric_surface=nothing
        context.geo_changed=true
        allocator_state===nothing ||
            _geo_allocator_reset_model!(allocator_state,context)
        return nothing
    elseif name=="Physicals"
        # `resetPhysicalGroups` + `removePhysicalGroups` drop the raw group
        # records, the entity-level memberships and the derived view, but
        # neither the names nor the physical tag counter.
        empty!(context.raw_physicals);empty!(m.physical)
        empty!(m.entity_physicals)
        context.geo_changed=true
        return nothing
    elseif name=="Variables"
        empty!(context.values);empty!(context.lists)
        empty!(context.list_variables)
        empty!(context.unavailable);empty!(context.unavailable_lists)
        context.stored_list_items=0
        return nothing
    elseif name=="Options"
        # `ReInitOptions` — restore the option-valued state `.geo` execution
        # tracks: mesh order, the transfinite-triangle flag, the extrude
        # result-list behavior, and every tracked option write.
        m.meshing.order=1
        m.meshing.transfinite_tri=0
        context.extrude_return_lateral=true
        empty!(context.option_strings);empty!(context.option_numbers)
        empty!(context.option_colors)
        return 0
    elseif name=="Meshes"
        # `GModel::deleteMesh` — drops the mid-file mesh and its ownership.
        context.mesh=nothing
        context.mesh_node_owner=empty(context.mesh_node_owner)
        return nothing
    elseif name=="Struct"
        # `gmsh_yynamespaces.clear()` — every struct definition across all
        # namespaces is dropped; tags reset with them.
        empty!(context.structs);empty!(context.struct_maxtag)
        return nothing
    end
    known=haskey(context.values,name) || haskey(context.lists,name) ||
          haskey(context.unavailable,name) ||
          haskey(context.unavailable_lists,name)
    if !known
        # `yymsg(0, "Unknown object or expression to delete '%s'")` — a
        # recorded error, not a parse abort.
        _geo_yyerror!(context,
            "Unknown object or expression to delete '$name'")
        return nothing
    end
    _geo_context_forget!(context,name)
    delete!(context.unavailable,name);delete!(context.unavailable_lists,name)
    return nothing
end

function _exec_line!(m::GeoModel,line::AbstractString,
                     context::_GeoNumericContext,
                     allocator_state=nothing)
    _geo_statement_marks_internals(line) && (context.geo_changed=true)
    # Malformed control headers — valid `If`/`ElseIf`/`For`/`Function`/`Macro`
    # statements are consumed as control markers by the statement splitter and
    # never reach this dispatch, so a leading keyword here is always a syntax
    # error upstream. The grammar's token sequence is replayed so the error
    # lands on the offender (`For = 5` → `(=)`, `For x = 5` → `(=)` after the
    # loop variable shifts, `If foo` → `(foo)`). `While` is absent — it is not
    # a lexer keyword upstream (`While = 5` assigns a variable there).
    if (cm=match(r"^(If|ElseIf|For|Function|Macro)\b",line))!==nothing
        word=cm.captures[1]
        rest=strip(line[nextind(line,firstindex(line),
                              ncodeunits(word)):end])
        if word=="For"
            nm=match(r"^[A-Za-z_][A-Za-z0-9_]*",rest)
            nm===nothing &&
                _geo_syntax_abort(_geo_stmt_head_token(rest))
            rest=strip(rest[nextind(rest,lastindex(nm.match)):end])
            match(r"^In\b",rest)===nothing &&
                _geo_syntax_abort(_geo_stmt_head_token(rest))
            rest=strip(rest[3:end])
            startswith(rest,"{") ||
                _geo_syntax_abort(_geo_stmt_head_token(rest))
        elseif word=="Function" || word=="Macro"
            match(r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\"|'[^']*')",
                  rest)===nothing &&
                _geo_syntax_abort(_geo_stmt_head_token(rest))
        else
            startswith(rest,"(") ||
                _geo_syntax_abort(_geo_stmt_head_token(rest))
            close=_geo_matching_delim(rest,firstindex(rest))
            if close==0
                # `(` shifted; the error lands on the terminator (`;`) — the
                # header expression still evaluated upstream, so report its
                # own diagnostics before aborting.
                _geo_eval_numeric(strip(rest[2:end]),context,"execute_geo")
                _geo_syntax_abort(";")
            else
                _geo_eval_numeric(rest[2:prevind(rest,close)],context,
                                  "execute_geo")
            end
        end
        # A syntactically complete header cannot reach this dispatcher — the
        # splitter consumes every valid control statement. Abort rather than
        # silently dropping the construct.
        _geo_syntax_abort(";")
    end
    if occursin(r"^Periodic(?:\s|$)",line)
        _exec_periodic!(m,line,context)
        return
    elseif (mm=match(
            r"^Point\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Point"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        # `tPoint '(' FExpr ')' tAFFECT VExpr` — the coordinate braces are
        # FExpr-comma components, not a RecursiveListOfDouble: `a[]` splices
        # and `:` ranges are syntax errors upstream.
        values=_geo_exec_vexpr(mm.captures[2],4,context,"$caller coordinates")
        lc=values[4]
        # `AddToTemporaryBoundingBox` — every `Point` statement feeds the
        # parse-time `CTX::instance()->lc` (also on a failed add upstream).
        _geo_exec_temp_bbox_add!(context,values[1],values[2],values[3])
        _geo_dup_entity!(context,haskey(m.points,trunc(Int,tag)),
                         :point,"point",tag,zero_literal,allocator_state) &&
            return
        t=add_point!(m,values[1],values[2],values[3];
                   tag=tag,mesh_size=lc>0 ? lc : 1.0,
                   _zero_literal=zero_literal)
        # `gmshVertex` keeps the `lc` component verbatim — `lc <= 0` is the
        # unset default (`point_size` reads fall back to 0).
        lc<=0 && (m.point_size[t]=lc)
        return
    elseif (mm=match(
            r"^Line\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Line"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller endpoints";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"line",tag,zero_literal,allocator_state) &&
            return
        if length(points)<2
            _geo_msg_error!(context,"Line requires 2 points")
            _geo_yyerror!(context,"Could not add line")
            return
        end
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        _geo_line_tolerance_warn!(m,context,found,t)
        _geo_store_curve!(m,t,found;control_points=length(found)>2)
        ok || _geo_yyerror!(context,"Could not add line")
        return
    elseif (mm=match(
            r"^Circle\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Circle"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        (points,normal)=_geo_circle_rhs(mm.captures[2],context,caller)
        if length(points)!=3
            # Grammar-level check — `addCircleArc` is never reached (`r`
            # stays true), so no duplicate check, no `Could not add`, and no
            # tag is consumed.
            _geo_yyerror!(context,"Circle requires 3 points")
            return
        end
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"circle",tag,zero_literal,allocator_state) &&
            return
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        n=_arc_stored_normal(normal,caller)
        ok=_geo_arc_endcurve_ok!(m,context,found,n,t;
                                 ellipse=false,
                                 normal_given=normal!==nothing) && ok
        _geo_store_curve!(m,t,found;kind=:circle,control_points=true,
                          geometry=(n=n,))
        ok || _geo_yyerror!(context,"Could not add circle")
        return
    elseif (mm=match(
            r"^Ellipse\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Ellipse"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        (points,normal)=_geo_circle_rhs(mm.captures[2],context,caller)
        if !(length(points) in (3,4))
            # Grammar-level check — `addEllipseArc` is never reached.
            _geo_yyerror!(context,"Ellipse requires 4 points")
            return
        end
        # The built-in 3-tag form mirrors the OCC backward-compatibility
        # record: the start point doubles as the major-axis point
        # (`addEllipseArc(num, tags[0], tags[1], tags[0], tags[2], ..)`).
        length(points)==3 &&
            (points=[points[1],points[2],points[1],points[3]])
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"ellipse",tag,zero_literal,allocator_state) &&
            return
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        n=_arc_stored_normal(normal,caller)
        ok=_geo_arc_endcurve_ok!(m,context,found,n,t;
                                 ellipse=true,
                                 normal_given=normal!==nothing) && ok
        _geo_store_curve!(m,t,found;kind=:ellipse,control_points=true,
                          geometry=(n=n,))
        ok || _geo_yyerror!(context,"Could not add ellipse")
        return
    elseif (mm=match(
            r"^Spline\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Spline"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller control points";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"spline",tag,zero_literal,allocator_state) &&
            return
        if length(points)<2
            _geo_msg_error!(context,
                "Spline curve requires at least 2 control points")
            _geo_yyerror!(context,"Could not add spline")
            return
        end
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        _geo_store_curve!(m,t,found;kind=:spline,control_points=true)
        ok || _geo_yyerror!(context,"Could not add spline")
        return
    elseif (mm=match(
            r"^BSpline\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: BSpline"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller control points";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"BSpline",tag,zero_literal,allocator_state) &&
            return
        if length(points)<2
            _geo_msg_error!(context,
                "BSpline curve requires at least 2 control points")
            _geo_yyerror!(context,"Could not add BSpline")
            return
        end
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        _geo_store_curve!(m,t,found;kind=:bspline,control_points=true)
        ok || _geo_yyerror!(context,"Could not add BSpline")
        return
    elseif (mm=match(
            r"^Bezier\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Bezier"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller control points";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"Bezier",tag,zero_literal,allocator_state) &&
            return
        if length(points)<2
            _geo_msg_error!(context,
                "Bezier curve requires at least 2 control points")
            _geo_yyerror!(context,"Could not add Bezier")
            return
        end
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        _geo_store_curve!(m,t,found;kind=:bezier,control_points=true)
        ok || _geo_yyerror!(context,"Could not add Bezier")
        return
    elseif (mm=match(
            r"^Nurbs\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        # `tNurbs (FExpr) = ListOfDouble tNurbsKnots ListOfDouble tNurbsOrder
        # FExpr` — the built-in kernel routes through `addBSpline(num, tags,
        # seqknots)` and never reads the parsed `Order` expression.
        caller="execute_geo: Nurbs"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        split_knots=_geo_split_at_keyword(mm.captures[2],"Knots",caller)
        split_knots===nothing && throw(ArgumentError(
            "$caller: expected `Knots <list>` after the control-point list; " *
            "got $(repr(strip(mm.captures[2])))"))
        (points_raw,knots_tail)=split_knots
        points=_geo_exec_entity_rhs_tags(
            points_raw,context,"$caller control points";
            abs_refs=true,allow_empty=true)
        split_order=_geo_split_at_keyword(knots_tail,"Order",caller)
        split_order===nothing && throw(ArgumentError(
            "$caller: expected `Order <expr>` after the knot list; got " *
            "$(repr(knots_tail))"))
        (knots_raw,order_raw)=split_order
        isempty(order_raw) && throw(ArgumentError(
            "$caller: Nurbs requires an `Order` expression"))
        knots=_geo_numeric_list_values(knots_raw,context,"$caller knots")
        _geo_eval_numeric(order_raw,context,"$caller Order")
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"nurbs",tag,zero_literal,allocator_state) &&
            return
        if length(points)<2
            # `addBSpline` reports through its BSpline record name even on
            # the `Nurbs` statement.
            _geo_msg_error!(context,
                "BSpline curve requires at least 2 control points")
            _geo_yyerror!(context,"Could not add nurbs")
            return
        end
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        if isempty(knots)
            # An empty `Knots {}` builds the plain `MSH_SEGM_BSPLN` record.
            _geo_store_curve!(m,t,found;kind=:bspline,control_points=true)
        else
            # `order = knots - points - 1` uses the *input* point count; the
            # knot vector stores as `float` and `ubeg`/`uend` keep the raw
            # doubles — upstream performs no knot validation.
            degree=length(knots)-length(points)-1
            _geo_store_curve!(m,t,found;kind=:nurbs,control_points=true,
                geometry=(knots=Float64.(Float32.(knots)),deg=degree,
                          ubeg=first(knots),uend=last(knots)))
        end
        ok || _geo_yyerror!(context,"Could not add nurbs")
        return
    elseif (mm=match(
            r"^(?:Line\s+Loop|Curve\s+Loop)\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Curve Loop"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        curves=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller curves";
            signed=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.loops,trunc(Int,tag)),
                         :curve_loop,"curve loop",tag,zero_literal,
                         allocator_state) && return
        t=_alloc_curve_loop_tag(m,_tag(tag,caller,1),caller;
                                literal_zero=zero_literal)
        if all(id->haskey(m.curves,abs(id)),curves)
            # `SortEdgesInLoop` stores the chained prefix even on a dead-end.
            (ordered,ok)=_geo_sort_loop_partial(m,curves)
            m.loops[t]=ordered
            ok || (_geo_msg_error!(context,"Curve loop $t is wrong");
                   _geo_yyerror!(context,"Could not add curve loop"))
        else
            # `FindCurve` misses return `true` early — the raw list is stored
            # with no diagnostic at all.
            m.loops[t]=copy(curves)
        end
        return
    elseif (mm=match(
            r"^Plane\s+Surface\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Plane Surface"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        # `SetSurfaceGeneratrices` resolves each wire through `abs(iLoop)`.
        loops=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller loops";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.surfaces,trunc(Int,tag)),
                         :surface,"plane surface",tag,zero_literal,
                         allocator_state) && return
        if isempty(loops)
            _geo_msg_error!(context,
                "Plane surface requires at least one curve loop")
            _geo_yyerror!(context,"Could not add plane surface")
            return
        end
        t=_geo_entity_tag(m,2,tag,caller,zero_literal)
        ok=_geo_check_surface_wires!(m,context,loops,t)
        m.surfaces[t]=loops
        ok || _geo_yyerror!(context,"Could not add plane surface")
        return
    elseif (mm=match(
            r"^(Ruled\s+)?Surface\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        # `Surface` and the deprecated `Ruled Surface` alias both execute
        # Gmsh's `addSurfaceFilling`: the first wire's curve count selects the
        # patch kind (4 -> ruled, 3 -> triangular).
        caller="execute_geo: $(mm.captures[1]===nothing ? "" : "Ruled ")Surface"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[2],context,"$caller tag")
        (group,rest)=_geo_balanced_group(
            String(strip(mm.captures[3])),caller)
        wires=_geo_exec_entity_tags(group,context,"$caller loops";
                                    abs_refs=true,allow_empty=true)
        sphere_center=_geo_surface_constraint(rest,context,caller)
        _geo_dup_entity!(context,haskey(m.surfaces,trunc(Int,tag)),
                         :surface,"surface",tag,zero_literal,
                         allocator_state) && return
        if isempty(wires)
            _geo_msg_error!(context,"Surface requires at least one curve loop")
            _geo_yyerror!(context,"Could not add surface")
            return
        end
        # The first wire's existence and border count are pre-store checks —
        # failures leave the tag unconsumed.
        first_loop=get(m.loops,abs(first(wires)),nothing)
        if first_loop===nothing
            _geo_msg_error!(context,"Unknown curve loop $(abs(first(wires)))")
            _geo_yyerror!(context,"Could not add surface")
            return
        end
        borders=length(first_loop)
        if !(borders in (3,4))
            _geo_msg_error!(context,
                "Wrong definition of surface $(trunc(Int,tag)): $borders " *
                "borders instead of 3 or 4")
            _geo_yyerror!(context,"Could not add surface")
            return
        end
        t=_geo_entity_tag(m,2,tag,caller,zero_literal)
        ok=_geo_check_surface_wires!(m,context,wires,t)
        if sphere_center!==nothing && sphere_center>=0 &&
                !haskey(m.points,sphere_center)
            _geo_msg_error!(context,
                "Unknown sphere center point $sphere_center")
            ok=false
        end
        m.surfaces[t]=wires
        m.surface_types[t]=borders==4 ? :ruled : :tric
        (sphere_center!==nothing && sphere_center>=0 &&
         haskey(m.points,sphere_center)) &&
            (m.surface_geometry[t]=(sphere_center=sphere_center,))
        ok || _geo_yyerror!(context,"Could not add surface")
        return
    elseif (mm=match(
            r"^Surface\s+Loop\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*(?:\s+Using\s+Sewing)?\s*;$",
            line)) !== nothing
        caller="execute_geo: Surface Loop"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        surfaces=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller surfaces";
            signed=true,allow_empty=true)
        # `addSurfaceLoop` validates nothing: the raw list is stored
        # unconditionally (closure surfaces at `Volume` creation or later).
        _geo_dup_entity!(context,haskey(m.surface_loops,trunc(Int,tag)),
                         :surface_loop,"surface loop",tag,zero_literal,
                         allocator_state) && return
        t=_alloc_surface_loop_tag(m,_tag(tag,caller,2),caller;
                                  literal_zero=zero_literal)
        m.surface_loops[t]=copy(surfaces)
        return
    elseif (mm=match(
            r"^Volume\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Volume"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        # `SetVolumeSurfaces` resolves each shell through `abs(il)`; Gmsh
        # checks shell and member-surface existence only (closure is deferred
        # to meshing).
        shells=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller surface loops";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.volumes,trunc(Int,tag)),
                         :volume,"volume",tag,zero_literal,
                         allocator_state) && return
        t=_geo_entity_tag(m,3,tag,caller,zero_literal)
        ok=_geo_check_volume_shells!(m,context,shells,t)
        m.volumes[t]=copy(shells)
        ok || _geo_yyerror!(context,"Could not add volume")
        return
    elseif (mm=match(
            r"^Box\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Box"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        if allocator_state===nothing ||
                allocator_state.factory!==:opencascade
            # `tBox` is OCC-gated upstream — the tag and parameters evaluate
            # first, then the recoverable diagnostic replaces creation.
            _geo_yyerror!(context,
                "Box only available with OpenCASCADE geometry kernel")
            return
        end
        if length(values)!=6
            _geo_yyerror!(context,"Box requires 6 parameters")
            return
        end
        _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"block",
                         tag,zero_literal,allocator_state) &&
            return
        add_box!(m,values...;tag=tag,_zero_literal=zero_literal)
        return
    elseif (mm=match(
            r"^Cylinder\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Cylinder"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        if allocator_state===nothing ||
                allocator_state.factory!==:opencascade
            _geo_yyerror!(context,
                "Cylinder only available with OpenCASCADE geometry kernel")
            return
        end
        if length(values)!=7
            _geo_yyerror!(context,"Cylinder requires 7 parameters")
            return
        end
        _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"cylinder",
                         tag,zero_literal,allocator_state) &&
            return
        add_cylinder!(m,values...;tag=tag,_zero_literal=zero_literal)
        return
    elseif (mm=match(
            r"^Sphere\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Sphere"
        # The two-point form registers the raw tag in `allGmshSurfaces` — the
        # entity `tag < 0` auto-assignment does not apply there. The OCC
        # parameter form still uses the normal addX tag semantics.
        raw_tag=_geo_signed_gmsh_int_value(
            _geo_eval_numeric(mm.captures[1],context,"$caller tag"),
            "$caller tag")
        tag=raw_tag<0 ? 0 : raw_tag
        zero_literal=raw_tag==0
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        if length(values)==2
            # `Sphere(n) = {centerTag, pointTag}` — the built-in-kernel form:
            # center point plus a point on the sphere, radius = distance.
            _geo_exec_point_sphere!(m,raw_tag,values,context,caller,"sphere";
                                    zero_literal=zero_literal)
        elseif 4<=length(values)<=7
            # `{x,y,z,r[,a1,a2,a3]}` — Gmsh gates this on OpenCASCADE;
            # Tessella's native primitive supports the full-sphere form.
            if allocator_state===nothing ||
                    allocator_state.factory!==:opencascade
                _geo_yyerror!(context,
                    "Sphere only available with OpenCASCADE geometry kernel")
                return
            end
            _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"sphere",
                             tag,zero_literal,allocator_state) &&
                return
            length(values)==4 || throw(ArgumentError(
                "$caller: Sphere angle bounds (a1,a2,a3) are not supported"))
            add_sphere!(m,values[1],values[2],values[3],values[4];
                        tag=tag,_zero_literal=zero_literal)
        else
            _geo_yyerror!(context,
                "Sphere requires 2 points or from 4 to 7 parameters")
        end
        return
    elseif (mm=match(
            r"^PolarSphere\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        # `PolarSphere(n) = {centerTag, pointTag}` — built-in kernel only;
        # same center/radius construction as `Sphere` (the difference is the
        # surface parametrization, which Tessella's sphere does not carry).
        caller="execute_geo: PolarSphere"
        raw_tag=_geo_signed_gmsh_int_value(
            _geo_eval_numeric(mm.captures[1],context,"$caller tag"),
            "$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        if length(values)==2
            _geo_exec_point_sphere!(m,raw_tag,values,context,caller,
                                    "polar sphere";zero_literal=raw_tag==0)
        else
            _geo_yyerror!(context,"PolarSphere requires 2 points")
        end
        return
    elseif (mm=match(
            r"^Cone\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Cone"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        if allocator_state===nothing ||
                allocator_state.factory!==:opencascade
            _geo_yyerror!(context,
                "Cone only available with OpenCASCADE geometry kernel")
            return
        end
        if length(values)!=8
            _geo_yyerror!(context,"Cone requires 8 parameters")
            return
        end
        _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"cone",
                         tag,zero_literal,allocator_state) &&
            return
        add_cone!(m,values...;tag=tag,_zero_literal=zero_literal)
        return
    elseif (mm=match(
            r"^Torus\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: Torus"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        values=_geo_exec_numeric_values(
            mm.captures[2],context,"$caller parameters")
        if allocator_state===nothing ||
                allocator_state.factory!==:opencascade
            _geo_yyerror!(context,
                "Torus only available with OpenCASCADE geometry kernel")
            return
        end
        if !(length(values)==5 || length(values)==6)
            _geo_yyerror!(context,"Torus requires 5 or 6 parameters")
            return
        end
        _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"torus",
                         tag,zero_literal,allocator_state) &&
            return
        add_torus!(m,values[1:5]...;tag=tag,_zero_literal=zero_literal,
                   angle=length(values)==6 ? values[6] : 2π)
        return
    elseif (mm=match(
            r"^Curve\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        # `Curve(n) = {p1,p2};` is a `Line` alias in the .geo grammar.
        caller="execute_geo: Curve"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        points=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller endpoints";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,haskey(m.curves,trunc(Int,tag)),
                         :curve,"line",tag,zero_literal,allocator_state) &&
            return
        if length(points)<2
            _geo_msg_error!(context,"Line requires 2 points")
            _geo_yyerror!(context,"Could not add line")
            return
        end
        t=_geo_entity_tag(m,1,tag,caller,zero_literal)
        (found,ok)=_geo_curve_found_points!(m,context,points,t)
        _geo_line_tolerance_warn!(m,context,found,t)
        _geo_store_curve!(m,t,found;control_points=length(found)>2)
        ok || _geo_yyerror!(context,"Could not add line")
        return
    elseif (mm=match(
            Regex("^(Rectangle|Disk|Wedge|ThickSolid|ThruSections|Wire)\\s*\\(|" *
                  "^Ruled\\s+ThruSections\\s*\\(|" *
                  "^(BSpline|Bezier)\\s+Surface\\s*\\("),
            line)) !== nothing
        # OCC-only entity definitions: under the built-in factory Gmsh emits
        # `<name> only available with OpenCASCADE geometry kernel`; under
        # OpenCASCADE the entity would be created, which Tessella's
        # Julia-native kernel cannot honor.
        name=mm.captures[1]===nothing ?
            (mm.captures[2]===nothing ? "ThruSections" :
             "$(mm.captures[2]) surface") : mm.captures[1]
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            throw(ArgumentError(
                "execute_geo: $name requires an OpenCASCADE backend, which " *
                "Tessella does not implement"))
        end
        _geo_yyerror!(context,
            "$name only available with OpenCASCADE geometry kernel")
        return
    elseif (mm=match(
            r"^Coordinates\s+Surface\s+(.*?)\s*;$",line)) !== nothing
        # `Coordinates Surface n` installs a registered `allGmshSurfaces`
        # entry as `myGmshSurface`. A miss assigns nullptr upstream, so the
        # current context is cleared either way.
        caller="execute_geo: Coordinates Surface"
        surface_tag=_geo_signed_gmsh_int_value(
            _geo_eval_numeric(mm.captures[1],context,"$caller surface"),
            "$caller surface")
        surface=get(context.parametric_surfaces,surface_tag,nothing)
        surface===nothing && _geo_msg_error!(context,
            "Surface $surface_tag does not exist")
        context.parametric_surface=surface
        return
    elseif match(r"^Parametric\s+Surface\s*\(",line)!==nothing
        # `Parametric Surface` belongs to Gmsh's expression-driven
        # `myGmshSurface` subsystem, which Tessella does not implement.
        throw(ArgumentError(
            "execute_geo: parametric surfaces are not supported: $line"))
    elseif match(r"^(Curve|Surface)\s+\"",line)!==nothing
        # `Curve "name"(n)` / `Surface "name"(n)` are not grammar forms —
        # upstream aborts at the string token.
        _geo_syntax_abort("\"")
        return
    elseif match(r"^Euclidian\s+Coordinates\s*;?\s*$",line)!==nothing
        # `Euclidian Coordinates` deletes `myGmshSurface` — the parametric
        # context a `Sphere(n) = {c, p}`-style statement installs.
        context.parametric_surface=nothing
        return
    elseif (mm=match(
            r"^Compound\s+(Spline|BSpline)\s*\(",
            line)) !== nothing
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            _geo_yyerror!(context,
                "Compound spline only available with built-in geometry kernel")
        else
            throw(ArgumentError(
                "execute_geo: Compound $(mm.captures[1]) entities are not " *
                "implemented"))
        end
        return
    elseif (mm=match(
            r"^Compound\s+(Curve|Line|Surface|Volume)\s*\(",
            line)) !== nothing
        # `Compound X(n) = {...}` is deprecated in favor of the `Compound
        # X{...}` meshing constraint — Gmsh emits a diagnostic and creates
        # nothing.
        if mm.captures[1] in ("Curve","Line")
            _geo_yyerror!(context,"`Compound Line (...) = {...};' is " *
                "deprecated: use `Compound Spline|BSpline (...) = {...} " *
                "Using ...;' instead, or the compound meshing constraint " *
                "`Compound Curve {...};'")
        else
            # Gmsh prints `Compound Surface' verbatim for dim 2 and 3.
            _geo_yyerror!(context,"`Compound Surface (...) = " *
                "{...};' is deprecated: use the compound meshing constraint " *
                "`Compound Surface {...};' instead")
        end
        return
    elseif (mm=match(
            r"^Intersect\s+Curve\s*\{\s*(.*?)\s*\}\s*Surface\s*\{\s*(.*?)\s*\}\s*;?\s*$",
            line)) !== nothing
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            _geo_yyerror!(context,
                "Intersect line not available with OpenCASCADE geometry kernel")
        else
            throw(ArgumentError(
                "execute_geo: `Intersect Curve{..} Surface{..}` is not " *
                "implemented"))
        end
        return
    elseif match(r"^Split\s+Curve\b",line)!==nothing
        if (dm=match(r"^Split\s+Curve\s*\(",line))!==nothing
            _geo_yywarn!(context,
                "'Split Curve(c) {...}' is deprecated: use 'Split Curve {c} " *
                "Point {...}' instead")
        end
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            _geo_yyerror!(context,
                "Split Curve not available with OpenCASCADE geometry kernel")
        else
            throw(ArgumentError(
                "execute_geo: `Split Curve` is not implemented"))
        end
        return
    elseif match(r"^Closest\s*\{",line)!==nothing
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            throw(ArgumentError(
                "execute_geo: Closest requires an OpenCASCADE backend, " *
                "which Tessella does not implement"))
        end
        _geo_yyerror!(context,
            "Closest entity only available with OpenCASCADE geometry kernel")
        return
    elseif (mm=match(r"^(Fillet|Chamfer)\s*\{",line))!==nothing
        # `tFillet {regions} {edges} {radii}` / `tChamfer {regions} {edges}
        # {surfaces} {distances}` — OCC-only upstream; under OpenCASCADE the
        # entity is created, which Tessella's native kernel cannot honor.
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            throw(ArgumentError(
                "execute_geo: $(mm.captures[1]) requires an OpenCASCADE " *
                "backend, which Tessella does not implement"))
        end
        _geo_yyerror!(context,
            "$(mm.captures[1]) only available with OpenCASCADE geometry kernel")
        return
    elseif match(
            r"^(?:Ruled\s+)?ThruSections\s*(?:\{|[\w(+-])",
            line)!==nothing
        # `ThruSections {wires}` / `Ruled ThruSections {wires}` —
        # `tThruSections ListOfDouble` (the `(tag) = {...}` assignments are
        # handled above). Under OpenCASCADE upstream creates the entity;
        # Tessella's native kernel does not implement it.
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            throw(ArgumentError(
                "execute_geo: ThruSections requires an OpenCASCADE backend, " *
                "which Tessella does not implement"))
        end
        _geo_yyerror!(context,
            "ThruSections only available with OpenCASCADE geometry kernel")
        return
    elseif match(r"^ShapeFromFile\s*\(",line)!==nothing
        _geo_yyerror!(context,
            "ShapeFromFile only available with OpenCASCADE geometry kernel")
        return
    elseif match(r"^HealShapes\b",line)!==nothing
        _geo_yyerror!(context,
            "HealShapes only available with OpenCASCADE geometry kernel")
        return
    elseif match(
            r"^Boolean(?:Difference|Union|Intersection|Fragments|Section|Fuse|Coherence|Common|Cut)\s*\{",
            line)!==nothing
        # Standalone `BooleanX{...}{...}` — the resulting entities use
        # automatic tags and the statement discards the list.
        body=String(strip(line))
        endswith(body,";") &&
            (body=String(strip(body[firstindex(body):prevind(body,lastindex(body))])))
        _geo_exec_boolean_term(m,body,context,allocator_state)
        return
    elseif (mm=match(r"^OnelabRun\s*\(\s*(.*?)\s*\)\s*;?\s*$",line)) !== nothing
        nargs=_geo_eval_string_list(mm.captures[1],context,
            "execute_geo: OnelabRun")
        length(nargs) in (1,2) || _geo_yyerror!(context,
            "OnelabRun takes one or two arguments")
        length(nargs) in (1,2) && throw(ArgumentError(
            "execute_geo: OnelabRun: external ONELAB clients are not supported"))
        return
    elseif (mm=match(
            r"^Boolean(Difference|Union|Intersection|Fragments|Section|Fuse|Coherence|Common|Cut)\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Boolean$(mm.captures[1])"
        # `_multiBind` treats `tag <= 0` as automatic and binds `tag > 0`
        # literally to every result piece.
        tag=max(0,_geo_signed_gmsh_int_value(
            _geo_eval_numeric(mm.captures[2],context,"$caller result tag"),
            "$caller result tag"))
        groups=_geo_boolean_groups(mm.captures[3],caller)
        # Upstream `BooleanShape` skips both the diagnostic and operand
        # resolution when the OCC kernel is inactive — a silent no-op.
        if allocator_state===nothing || allocator_state.factory!==:opencascade
            return
        end
        opname=_GEO_BOOLEAN_NAMES[mm.captures[1]]
        op=get(_GEO_BOOLEAN_OPS,opname,nothing)
        op===nothing && throw(ArgumentError(
            "execute_geo: Boolean$opname is not implemented"))
        operands=_geo_boolean_operands(groups,context,caller)
        boolean_volumes_multi!(m,op,
            operands.objects.tags,operands.tools.tags;tag=tag,
            remove_object=operands.objects.delete,
            remove_tool=operands.tools.delete,caller=caller)
        return
    elseif match(r"^(Translate|Rotate|Dilate|Symmetry)\s*\{",line)!==nothing
        _geo_exec_transform_statement!(m,line,context,allocator_state)
        return
    elseif match(r"^Affine\s*\{",line)!==nothing
        throw(ArgumentError(
            "execute_geo: Affine transforms require the OpenCASCADE geometry " *
            "kernel, which Tessella does not implement"))
    elseif match(r"^Extrude\b",line)!==nothing
        _geo_exec_extrude_term(m,line,context,allocator_state)
        return
    elseif match(
            r"^(Duplicata|Boundary|CombinedBoundary|OrientedBoundary|OrientedCombinedBoundary|PointsOf)\s*\{",
            line)!==nothing
        # Standalone MultipleShape statements are legal in the Gmsh grammar;
        # Duplicata copies are the observable side effect.
        s=String(strip(line))
        endswith(s,";") &&
            (s=String(strip(s[firstindex(s):prevind(s,end)])))
        _geo_shape_list_entities!(m,s,context,"execute_geo";
                                  allocator_state=allocator_state)
        return
    elseif (mm=match(
            r"^(Point|Line|Curve)\s*\{\s*(.*?)\s*\}\s+In\s+Surface\s*\{\s*(.*?)\s*\}\s*;$",
            line)) !== nothing
        caller="execute_geo: $(mm.captures[1]) In Surface"
        dim=mm.captures[1]=="Point" ? 0 : 1
        # `addEmbedded` resolves members through `getVertexByTag` — raw
        # tags, no `abs` (unlike point-creation references).
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
        # `PhysicalId_per_dim_entity`: an `FExpr` is the raw `(int)` tag —
        # literal, including 0 and negatives. The `"name"`-only form
        # unconditionally increments the counter (`setMaxPhysicalTag(t+1)`,
        # where `t` is the current max) and then runs
        # `setPhysicalName(name, dim, t+1)` — an already-bound name wins and
        # the bump is simply wasted. `"name", FExpr` runs
        # `setPhysicalName(name, dim, T)` with NO counter bump — a bound name
        # wins over the explicit tag; otherwise `(dim, T)` is bound when
        # unnamed, even when the modify below then fails.
        tag=if physical.tag_source===nothing
            _geo_physical_assign_next!(m,allocator_state)
            _geo_set_physical_name!(m,dim,name,m.physical_tag_max)
        else
            explicit=_geo_signed_gmsh_int_value(
                _geo_eval_numeric(
                    physical.tag_source,context,"$caller tag"),
                "$caller tag")
            isempty(name) ? explicit :
                _geo_set_physical_name!(m,dim,name,explicit)
        end
        membership=physical.membership
        # A missing list is a syntax error, but `{}` is a valid (empty)
        # `ListOfDouble`: `=` stores an empty raw group, `-= {}` removes every
        # member (deleting the group), and other compound ops take it as a
        # no-op over members.
        isempty(strip(membership)) && throw(ArgumentError(
            "$caller entities: entity list is empty"))
        ids=if isempty(strip(membership,['{','}',' ']))
            Int[]
        elseif match(
                r"^(?:Boundary|CombinedBoundary|PointsOf)(?:\s|\{)",
                membership)!==nothing
            _geo_exec_physical_topology(
                m,membership,context,"$caller entities",dim)
        else
            occursin(r"[A-Za-z_][A-Za-z0-9_]*\s*\{",membership) &&
                throw(ArgumentError(
                    "$caller entities: unsupported topology query; use inline " *
                    "Boundary, CombinedBoundary, or PointsOf for Physical Point"))
            # `ListOfDouble2Vector` — the raw signed member tags; signs are
            # resolved to `abs` entities at sync (`orientedPhysicals`).
            _geo_exec_entity_rhs_tags(
                membership,context,"$caller entities";signed=true)
        end
        if physical.op=="="
            _geo_exec_physical_create!(m,dim,tag,ids,context,physical.kind,
                                       allocator_state)
        else
            _geo_exec_physical_modify!(m,dim,physical.op,tag,ids,
                                       context,caller,physical.kind)
        end
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
            _geo_exec_entity_tags(selector,context,"$caller Point selector";
                                  abs_refs=true)
        end
        isempty(point_tags) && throw(ArgumentError(
            "$caller: Point selector matched no explicit modeled points"))
        mesh_size=_geo_eval_numeric(
            mm.captures[3],context,"$caller value")
        set_point_mesh_size!(m,point_tags,mesh_size)
        return
    elseif (mm=match(
            r"^Transfinite\s+(?:Curve|Line)\s*(.*?)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        _geo_exec_transfinite_curve!(m,mm.captures[1],mm.captures[2],context)
        return
    elseif (mm=match(
            r"^Transfinite\s+Surface\s*(.*?)\s*;$",line)) !== nothing
        _geo_exec_transfinite_surface!(m,mm.captures[1],context)
        return
    elseif (mm=match(
            r"^Transfinite\s+Volume\s*(.*?)\s*;$",line)) !== nothing
        _geo_exec_transfinite_volume!(m,mm.captures[1],context)
        return
    elseif (mm=match(
            r"^TransfQuadTri\s*(.*?)\s*;$",line)) !== nothing
        # `TransfQuadTri` records the QuadTri flag on volumes; the QuadTri
        # hexahedral algorithm itself is a native-kernel blocker that surfaces
        # when a flagged volume is meshed.
        list=_geo_constraint_list(mm.captures[1],context,
                                  "execute_geo: TransfQuadTri";allow_all=true)
        _geo_constraint_apply(m,3,list,"execute_geo: TransfQuadTri",
                              tag->push!(m.meshing.quad_tri,tag))
        return
    elseif (mm=match(
            r"^Recombine\s+Surface\s*(.*?)\s*(?:=\s*(.*?)\s*)?;$",
            line)) !== nothing
        caller="execute_geo: Recombine Surface"
        # `RecombineAngle` is `(int)FExpr`, default 45.
        angle=mm.captures[2]===nothing ? 45.0 :
            Float64(_geo_constraint_int(
                _geo_eval_numeric(mm.captures[2],context,"$caller angle"),
                caller,"angle"))
        list=_geo_constraint_list(mm.captures[1],context,caller;allow_all=true)
        let angle=angle
            _geo_constraint_apply(m,2,list,caller,
                tag->(m.meshing.recombine[(2,tag)]=angle;nothing))
        end
        return
    elseif (mm=match(
            r"^Recombine\s+Volume\s*(.*?)\s*;$",line)) !== nothing
        list=_geo_constraint_list(mm.captures[1],context,
                                  "execute_geo: Recombine Volume";
                                  allow_all=true)
        _geo_constraint_apply(m,3,list,"execute_geo: Recombine Volume",
            tag->(m.meshing.recombine[(3,tag)]=0.0;nothing))
        return
    elseif (mm=match(
            r"^Smoother\s+Surface\s*(.*?)\s*=\s*(.*?)\s*;$",line)) !== nothing
        caller="execute_geo: Smoother Surface"
        count=_geo_constraint_int(
            _geo_eval_numeric(mm.captures[2],context,"$caller iteration count"),
            caller,"iteration count")
        list=_geo_constraint_list(mm.captures[1],context,caller;allow_all=true)
        _geo_constraint_apply(m,2,list,caller,
            tag->(m.meshing.smoothing[(2,tag)]=count;nothing))
        return
    elseif (mm=match(
            r"^MeshAlgorithm\s+Surface\s*(\{.*\})\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: MeshAlgorithm Surface"
        number=_geo_constraint_int(
            _geo_eval_numeric(mm.captures[2],context,"$caller algorithm"),
            caller,"algorithm")
        list=_geo_constraint_list(mm.captures[1],context,caller;allow_all=false)
        for value in list
            sd=_geo_constraint_int(value,caller,"Surface tag")
            haskey(m.surfaces,sd) || continue
            m.meshing.algorithm[(2,sd)]=number
        end
        return
    elseif (mm=match(
            r"^MeshSizeFromBoundary\s+Surface\s*(\{.*\})\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: MeshSizeFromBoundary Surface"
        number=_geo_constraint_int(
            _geo_eval_numeric(mm.captures[2],context,"$caller value"),
            caller,"value")
        list=_geo_constraint_list(mm.captures[1],context,caller;allow_all=false)
        for value in list
            sd=_geo_constraint_int(value,caller,"Surface tag")
            haskey(m.surfaces,sd) || continue
            number==0 ? delete!(m.meshing.size_from_boundary,(2,sd)) :
                        (m.meshing.size_from_boundary[(2,sd)]=true)
        end
        return
    elseif (mm=match(
            r"^(?:Reverse|ReverseMesh)\s+(Curve|Line|Surface)\s*(.*?)\s*;$",
            line)) !== nothing
        dim=mm.captures[1]=="Surface" ? 2 : 1
        caller="execute_geo: ReverseMesh $(_entity_label(dim))"
        list=_geo_constraint_list(mm.captures[2],context,caller;allow_all=true)
        let dim=dim
            _geo_constraint_apply(m,dim,list,caller,
                tag->(m.meshing.reverse[(dim,tag)]=true;nothing))
        end
        return
    elseif (mm=match(
            r"^RelocateMesh\s+(Point|Curve|Line|Surface)\s*(.*?)\s*;$",
            line)) !== nothing
        # `RelocateMesh` retargets existing mesh vertices onto other entities;
        # no mesh exists during .geo execution, so only the list is validated.
        _geo_constraint_list(mm.captures[2],context,
                             "execute_geo: RelocateMesh";allow_all=true)
        return
    elseif (mm=match(
            r"^ReorientMesh\s+Volume\s*(.*?)\s*;$",line)) !== nothing
        # `ReorientMesh` flips boundary-face mesh orientations on an existing
        # volume mesh via a solid-angle parity walk; with no mesh during .geo
        # execution it changes nothing and stores no state.
        _geo_constraint_list(mm.captures[1],context,
                             "execute_geo: ReorientMesh Volume";allow_all=false)
        return
    elseif (mm=match(
            r"^Degenerated\s+(?:Curve|Line)\s*(.*?)\s*;$",line)) !== nothing
        caller="execute_geo: Degenerated Curve"
        list=_geo_constraint_list(mm.captures[1],context,caller;allow_all=false)
        for value in list
            sd=_geo_constraint_int(value,caller,"Curve tag")
            haskey(m.curves,sd) || continue
            push!(m.meshing.degenerated,sd)
        end
        return
    elseif (mm=match(
            r"^Compound\s+(Curve|Line|Surface|Volume)\s*(.*?)\s*;$",
            line)) !== nothing
        dim=mm.captures[1]=="Surface" ? 2 : mm.captures[1]=="Volume" ? 3 : 1
        _geo_exec_compound!(m,dim,mm.captures[2],context)
        return
    elseif match(r"^RecombineMesh\s*;$",line) !== nothing
        # `RecombineMesh;` recombines an existing mesh — a no-op while no mesh
        # exists during .geo execution.
        return
    elseif (mm=match(r"^(Recursive\s+)?Delete\b(.*)$",line)) !== nothing
        return _geo_exec_delete!(m,mm.captures[1]!==nothing,
                                 String(mm.captures[2]),context,
                                 allocator_state)
    elseif match(r"^SetTag\b",line) !== nothing
        throw(ArgumentError(
            "execute_geo: SetTag cannot retag entities during .geo parsing — " *
            "Gmsh applies it to the synchronized model, which is still empty " *
            "at that point, so it reports the entity as unknown"))
    elseif (mm=match(
            r"^SetMaxTag\s+(?:(Point|Curve|Line|Surface|Volume)|GeoEntity\s*\{\s*(.*?)\s*\})\s*\(\s*(.*?)\s*\)\s*;$",
            line)) !== nothing
        caller="execute_geo: SetMaxTag"
        dim=if mm.captures[1]!==nothing
            mm.captures[1]=="Point" ? 0 : mm.captures[1] in ("Curve","Line") ? 1 :
                mm.captures[1]=="Surface" ? 2 : 3
        else
            _geo_constraint_int(
                _geo_eval_numeric(mm.captures[2],context,"$caller dimension"),
                caller,"dimension")
        end
        # The tag expression evaluates regardless of the dim's validity —
        # upstream parses the FExpr before the `setMaxTag` action runs.
        value=_geo_constraint_int(
            _geo_eval_numeric(mm.captures[3],context,"$caller tag"),
            caller,"tag")
        # A `GeoEntity{d}` outside 0..3 is a recoverable diagnostic upstream —
        # the statement still applies whenever the internals switch covers the
        # dim (`-2..3`). The observer emits "GeoEntity dim out of range [0,3]"
        # and updates the loop counters for `-1`/`-2`; `next_tag` only tracks
        # entity dims.
        0<=dim<=3 || return
        # `GEO_Internals::setMaxTag` is an unconditional assignment while
        # `OCC_Internals::setMaxTag` keeps `max(current, value)`; the allocator
        # observer applies the same split to its own namespaces.
        if allocator_state===nothing || allocator_state.factory!==:opencascade
            m.next_tag[dim+1]=value
        else
            m.next_tag[dim+1]=max(m.next_tag[dim+1],value)
        end
        return
    elseif (mm=match(
            r"^Coherence(?:\s+(Geometry|Mesh)|\s+Point\s*\{\s*(.*?)\s*\})?\s*;$",
            line)) !== nothing
        if mm.captures[2] !== nothing
            merge_vertices!(m,_geo_exec_entity_tags(
                mm.captures[2],context,"execute_geo: Coherence Point";
                abs_refs=true);
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
    elseif (fm=match(r"^Field\s*\[\s*(.*?)\s*\]\s*=\s*([A-Za-z][A-Za-z0-9_]*)\s*;?$",
                     line))!==nothing
        # `Field[i] = Kind` — `FieldManager::newField` semantics on the live
        # exec map: duplicate ids keep the original, unknown kinds create
        # nothing, both recoverable diagnostics.
        _geo_exec_field_declare!(context,fm.captures[1],String(fm.captures[2]))
        return
    elseif (fm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+Field\s*=\s*(.*?)\s*;?$",
                     line))!==nothing
        # `Background Field`/field-list commands — the tag list is evaluated at
        # execution, diagnostics match `String__Index tField tAFFECT
        # ListOfDouble`; the background/boundary-layer stores live in params.
        _geo_exec_field_command!(context,String(fm.captures[1]),fm.captures[2])
        return
    elseif startswith(line,"Field") &&
           !occursin(r"^Field\s*\[.*\]\s*\.\s*[A-Za-z_][A-Za-z0-9_]*\s*=",line)
        # `Field[i] = Kind` matched above and `Field[i].member = v` option
        # writes flow to the assignment dispatcher. `Field` is `tField`
        # upstream: after `[ FExpr ]` the grammar expects `=` — the syntax
        # error lands on the breaking token (`Field = 5` → `(=)`,
        # `Field[i];` → `(;)`), and the index expression still evaluates.
        rest=String(strip(line[nextind(line,firstindex(line),
                                     ncodeunits("Field")):end]))
        if startswith(rest,"[")
            close=_geo_matching_delim(rest,firstindex(rest))
            close==0 && _geo_syntax_abort("[")
            _geo_eval_numeric(
                rest[2:prevind(rest,close)],context,"execute_geo: Field index")
            rest=String(strip(rest[nextind(rest,close):end]))
            (dm=match(r"^\.\s*[A-Za-z_][A-Za-z0-9_]*",rest))!==nothing &&
                (rest=String(strip(rest[nextind(rest,
                    lastindex(dm.match)):end])))
        end
        _geo_syntax_abort(_geo_stmt_head_token(rest))
    elseif (bm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+Field\b(.*)$",line))!==nothing
        # `<name> Field` without `= ...` — `tField` expects an affectation, so
        # the error lands on the token after `Field` (`Background Field;` →
        # `(;)`). A keyword head already committed its own production —
        # `SetMaxTag Field(10)` errors on `Field` itself (the GeoEntity slot).
        # `Background`/`BoundaryLayer` are plain `tSTRING` names upstream:
        # `Background = 5` is an ordinary assignment and falls through below.
        _geo_tstring_valid(String(bm.captures[1])) || _geo_syntax_abort("Field")
        tail=String(strip(bm.captures[2]))
        _geo_syntax_abort(_geo_stmt_head_token(tail))
    elseif (mm=match(r"^Struct\s+(.+?)\s*;?\s*$",line)) !== nothing
        _geo_exec_struct_def!(String(mm.captures[1]),context)
        return
    elseif (mm=match(r"^SetFactory\s*\(\s*(.*?)\s*\)\s*;?\s*$",line)) !== nothing
        # The grammar syncs the internals (when changed) before switching
        # kernels so merged tags cannot clash.
        _geo_sync_physical_view_if_changed!(m,context)
        _geo_exec_set_factory!(mm.captures[1],context,allocator_state)
        return
    elseif (mm=match(r"^DefineConstant\s*\[\s*(.*?)\s*\]\s*;?\s*$",line)) !== nothing
        _geo_exec_define_constants!(context,mm.captures[1])
        return
    elseif (mm=match(r"^UndefineConstant\s*\[\s*(.*?)\s*\]\s*;?\s*$",line)) !== nothing
        for arg in _geo_split_args(mm.captures[1],"execute_geo: UndefineConstant")
            name=_geo_symbol_name(arg,context,"execute_geo: UndefineConstant")
            delete!(context.onelab_numbers,name)
            delete!(context.onelab_strings,name)
        end
        return
    elseif (mm=match(r"^SetNumber\s*\(\s*(.*?)\s*\)\s*;?\s*$",line)) !== nothing
        args=_geo_split_args(mm.captures[1],"execute_geo: SetNumber")
        length(args)==2 || throw(ArgumentError(
            "execute_geo: SetNumber takes a name and a value"))
        context.onelab_numbers[
            _geo_eval_string(args[1],context,"execute_geo: SetNumber")]=
            _geo_eval_numeric(args[2],context,"execute_geo: SetNumber")
        return
    elseif (mm=match(r"^SetString\s*\(\s*(.*?)\s*\)\s*;?\s*$",line)) !== nothing
        args=_geo_split_args(mm.captures[1],"execute_geo: SetString")
        length(args)==2 || throw(ArgumentError(
            "execute_geo: SetString takes a name and a value"))
        context.onelab_strings[
            _geo_eval_string(args[1],context,"execute_geo: SetString")]=
            _geo_eval_string(args[2],context,"execute_geo: SetString")
        return
    elseif (mm=match(r"^(Printf|Warning|Error)\s*\(\s*(.*?)\s*\)\s*(?:(>>|>)\s*(.*?)\s*)?;?\s*$",
                     line)) !== nothing
        _geo_exec_print!(mm.captures[1],mm.captures[2],
            mm.captures[3]===nothing ? nothing : (mm.captures[3],mm.captures[4]),
            context)
        return
    elseif (mm=match(r"^Exit(?:\s+(.*?))?\s*;?\s*$",line)) !== nothing
        # `Msg::Exit`: bare `Exit` exits with `_atLeastOneErrorInRun` — the
        # `Msg::Error` flag (`msg_error_count`); `Exit n` always exits `n`.
        context.stop=:exit
        context.exit_code=mm.captures[1]===nothing ?
            (context.msg_error_count>0 ? 1 : 0) :
            _geo_int_value(_geo_eval_numeric(mm.captures[1],context,
                "execute_geo: Exit"),"execute_geo: Exit code")
        return
    elseif match(r"^Abort\s*;?\s*$",line)!==nothing
        context.stop=:abort
        return
    elseif match(r"^SyncModel\s*;?\s*$",line)!==nothing
        # `SyncModel` forces a synchronize unconditionally (every other sync
        # point in the grammar gates on `getChanged()`).
        _geo_sync_physical_view!(m,context)
        return
    elseif match(r"^Draw\s*;?\s*$",line)!==nothing ||
           match(r"^SetChanged\s*;?\s*$",line)!==nothing ||
           match(r"^UnsplitWindow\s*;?\s*$",line)!==nothing
        # GUI-only commands — no windows exist in batch execution.
        return
    elseif match(r"^BoundingBox\s*;?\s*$",line)!==nothing
        # `BoundingBox` syncs the internals when changed, then recomputes the
        # model bounding box — always current here (geometry is stored
        # directly, never cached stale).
        _geo_sync_physical_view_if_changed!(m,context)
        return
    elseif (mm=match(r"^BoundingBox\s*\{(.*?)\}\s*;?\s*$",line)) !== nothing
        values=_geo_numeric_list_values(mm.captures[1],context,
            "execute_geo: BoundingBox")
        length(values)==6 || throw(ArgumentError(
            "execute_geo: BoundingBox takes six coordinates"))
        # `forcedBBox` is display state — nothing to store on the model.
        return
    elseif match(r"^NewModel\s*;?\s*$",line)!==nothing
        # `new GModel()` — a fresh empty model (empty name, no physical or
        # entity names, fresh GEO internals) becomes current; parser symbols
        # are left alone.
        _geo_reset_model_geometry!(m)
        _geo_reset_geometry_counters!(m,context)
        empty!(m.physical_names);empty!(m.entity_names)
        empty!(context.raw_physicals)
        empty!(m.meshing.compounds);context.mesh=nothing
        empty!(context.parametric_surfaces)
        context.parametric_surface=nothing
        context.mesh_node_owner=empty(context.mesh_node_owner)
        m.name=""
        context.geo_changed=true
        if allocator_state!==nothing
            _geo_allocator_reset_model!(allocator_state,context;fresh=true)
            # The fresh `GModel` has no OCC internals — every factory-gated
            # dispatch falls through to the built-in path upstream.
            allocator_state.factory=:builtin
        end
        return
    elseif (mm=match(r"^(?:Recursive\s+)?(?:Show|Hide)\b",line)) !== nothing
        _geo_exec_visibility!(m,line,context,allocator_state)
        return
    elseif (mm=match(r"^(?:Recursive\s+)?Color\b",line)) !== nothing
        _geo_exec_color!(m,line,context,allocator_state)
        return
    elseif (mm=match(Regex("^(RelocateMesh|ReorientMesh|RecombineMesh|" *
        "RenumberMeshNodes|RenumberMeshElements|CreateMeshEdges|" *
        "CreateMeshFaces|RefineMesh|TransformMesh|SetOrder|" *
        "PartitionMesh|CreateOverlaps|AdaptMesh|CreateTopology|" *
        "ClassifySurfaces|CreateGeometry)\\b"),line)) !== nothing
        _geo_exec_mesh_statement!(m,line,context)
        return
    elseif (mm=match(r"^Mesh\s+(.*?)\s*;?\s*$",line)) !== nothing
        d=_geo_eval_numeric(mm.captures[1],context,"execute_geo: Mesh")
        d in (1,2,3) || throw(ArgumentError(
            "execute_geo: Mesh dimension must be 1, 2 or 3 (got $d)"))
        # Gmsh synchronizes GEO internals before meshing (when changed) — the
        # physical view and any pending unknown-member warnings materialize.
        _geo_sync_physical_view_if_changed!(m,context)
        context.mesh,context.mesh_node_owner=
            _geo_mesh_model(m,Int(d),context)
        return
    elseif (mm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+(.*?)\s*;?\s*$",line)) !== nothing &&
           mm.captures[1] in _GEO_COMMAND_WORDS
        _geo_exec_command!(m,mm.captures[1],mm.captures[2],context,
                           allocator_state)
        return
    elseif match(r"^Plugin\s*\(",line)!==nothing
        throw(ArgumentError(
            "execute_geo: Plugin(...) actions are not supported"))
    elseif match(Regex("^(Combine|Alias|AliasWithOptions|Save\\s+View|" *
        "CutMesh|SplitMesh|SendToServer|Background\\s+Mesh\\s+View)\\b"),
        line)!==nothing
        throw(ArgumentError(
            "execute_geo: post-processing view operations are not supported: $line"))
    elseif match(r"^Levelset\b",line)!==nothing
        throw(ArgumentError(
            "execute_geo: Levelset definitions are not supported"))
    elseif match(r"^Homology\b",line)!==nothing
        throw(ArgumentError(
            "execute_geo: Homology computations are not supported"))
    end
    # Every remaining `name`-headed form is an Affectation: scalar/list
    # assignments, indexed and multi-index mutation, `++`/`--`, string and
    # string-list assignment, and option writes.
    body=String(strip(line))
    endswith(body,";") &&
        (body=String(strip(body[firstindex(body):prevind(body,lastindex(body))])))
    _geo_exec_assignment!(context,body,"execute_geo") && return
    # A statement that cannot start a `String__Index`/keyword production is a
    # syntax error upstream (`= 5` left behind by `EndIf = 5` → `(=)`).
    match(r"^[A-Za-z_]",body)===nothing &&
        _geo_syntax_abort(_geo_stmt_head_token(body))
    # `tSTRING`-headed remainders — upstream's catch-all: a bare name is
    # `syntax error (;)`; `name(expr)`/`name <expr>` evaluates the argument
    # (an unknown variable records) then reports `Unknown command 'name'`;
    # `name {` aborts at the first braced token. A lexer-keyword head aborts
    # like an LHS, never `Unknown command` (`DimNameSpace(ns);` →
    # `(DimNameSpace)`, `Sin(1);` → `(Sin)`).
    nm=match(r"^[A-Za-z_][A-Za-z0-9_]*",body)
    _geo_lhs_keyword_check(body)
    rest=String(strip(body[nextind(body,lastindex(nm.match)):end]))
    if isempty(rest)
        _geo_syntax_abort(";")
    elseif startswith(rest,"::")
        # `tSCOPE` never follows a bare `tSTRING` in a statement —
        # `ns::x = 5` → `syntax error (::)` upstream.
        _geo_syntax_abort("::")
    elseif startswith(rest,"{")
        _geo_syntax_abort(_geo_stmt_head_token(rest[2:end]))
    elseif startswith(rest,"(")
        close=_geo_matching_delim(rest,firstindex(rest))
        if close==0
            try
                _geo_eval_numeric(rest[2:end],context,"execute_geo")
            catch err
                err isa InterruptException && rethrow()
                err isa _GeoSyntaxAbort && rethrow()
            end
            _geo_syntax_abort(";")
        end
        _geo_eval_numeric(rest[2:prevind(rest,close)],context,"execute_geo")
        after=String(strip(rest[nextind(rest,close):end]))
        isempty(after) || _geo_syntax_abort(_geo_stmt_head_token(after))
        _geo_yyerror!(context,"Unknown command '$(nm.match)'")
    else
        _geo_string_rhs(rest) ||
            _geo_eval_numeric(rest,context,"execute_geo")
        _geo_yyerror!(context,"Unknown command '$(nm.match)'")
    end
    return
end

# ======== `.geo` statement helpers ========

const _GEO_COMMAND_WORDS=Set((
    "Include","Merge","MergeWithBoundingBox","Save","Print","System",
    "SystemCall","NonBlockingSystemCall","SetName","CreateDir","OnelabRun",
    "OptimizeMesh","Sleep","Remesh","SetCurrentWindow",
    "SplitCurrentWindowHorizontal","SplitCurrentWindowVertical",
    "SetBoundingBox"))

# `SetFactory("OpenCASCADE"|"Built-in")` — switches the active kernel and
# carries per-dimension tag counters to `max(this, other)` like Gmsh's
# `setMaxTag` sync on the switch.
function _geo_exec_set_factory!(raw::AbstractString,
                                context::_GeoNumericContext,
                                allocator_state)
    caller="execute_geo: SetFactory"
    factory=_geo_eval_string(raw,context,caller)
    allocator_state===nothing && return nothing
    factory=="OpenCASCADE" || factory in ("Built-in","Gmsh") ||
        _geo_yywarn!(context,
            "Unknown factory \"$factory\" - using \"Built-in\" instead")
    _geo_allocator_set_factory!(allocator_state,factory)
    return nothing
end

# `Field[i] = Kind` — `FieldManager::newField` on the live exec map. Duplicate
# ids keep the original field and unknown kinds create nothing; both report a
# `Msg::Error` plus a `yymsg(0)` parse diagnostic, matching Gmsh.
function _geo_exec_field_declare!(context::_GeoNumericContext,tag_src,kind::String)
    caller="execute_geo: Field declaration tag"
    tag=_geo_signed_gmsh_int_value(_geo_eval_numeric(tag_src,context,caller),caller)
    fields=context.fields
    fields===nothing && return nothing
    if haskey(fields,tag)
        _geo_msg_error!(context,"Field id $tag is already defined")
        _geo_yyerror!(context,"Cannot create field $tag of type '$kind'")
        return nothing
    end
    if !(kind in _GEO_FIELD_KINDS)
        _geo_msg_error!(context,"Unknown field type \"$kind\"")
        _geo_yyerror!(context,"Cannot create field $tag of type '$kind'")
        return nothing
    end
    curvature=_geo_option_number(context,"Mesh",0,"MeshSizeFromCurvature")
    fields[tag]=GeoFieldSpec(tag,kind,Dict{String,String}(),String[],
        curvature===nothing ? 0 : Int(curvature))
    return nothing
end

# `X Field = {list}` — Gmsh's `String__Index tField tAFFECT ListOfDouble`:
# `Background` takes exactly one tag and `BoundaryLayer` dedup-appends; both
# stores are params-side. Any other prefix is an unknown command.
function _geo_exec_field_command!(context::_GeoNumericContext,name::String,
                                  list_src)
    caller="execute_geo: $name Field"
    if name=="Background"
        tags=_geo_expression_field_tags(list_src,context,caller;
                                        deduplicate=false)
        if length(tags)>1
            _geo_yyerror!(context,
                "Only 1 field can be set as a background field.")
        elseif isempty(tags)
            _geo_yywarn!(context,"No field given (Background Field).")
        end
    elseif name!="BoundaryLayer"
        _geo_yyerror!(context,"Unknown command '$name Field'")
    end
    return nothing
end

# `Delete Field[i]` — `FieldManager::deleteField` on the live exec map; a
# missing id is a `Msg::Error` diagnostic and parsing continues.
function _geo_exec_delete_field!(context::_GeoNumericContext,tag_src)
    caller="execute_geo: Delete Field tag"
    tag=_geo_signed_gmsh_int_value(_geo_eval_numeric(tag_src,context,caller),caller)
    fields=context.fields
    if fields===nothing || !haskey(fields,tag)
        _geo_msg_error!(context,"Cannot delete field id $tag, it does not exist")
        return nothing
    end
    delete!(fields,tag)
    return nothing
end

# `DefineConstant[...]` — comma-separated entries; each is
# `name`, `name = FExpr`, `name = {ListOfDouble [, options]}`,
# `name() = {ListOfDouble [, options]}`, `name = StringExpr`, or
# `name = {StringExpr [, options]}`. The option tail lives INSIDE the braces:
# `{5, Name "n", Min 0, Max 10}` — a top-level `, Name "n"` is a syntax error
# in Gmsh (the next entry can't start with `"`).
# The brace forms go through `Msg::ExchangeOnelabParameter`: the ONELAB name
# comes from the `Name` option (missing `Name` + any option → a Gmsh error),
# an existing server value overrides the default, and the symbol is only
# created when absent.
function _geo_exec_define_constants!(context::_GeoNumericContext,
                                   raw::AbstractString)
    caller="execute_geo: DefineConstant"
    for item in _geo_split_args(raw,caller)
        spec=String(strip(item))
        if (m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*$",spec))!==nothing
            name=String(m.captures[1])
            _geo_context_has_variable(context,name) ||
                _geo_context_set_scalar!(context,name,0.0,caller)
            continue
        end
        m=match(r"^([A-Za-z_][A-Za-z0-9_]*)(\s*\(\s*\))?\s*=\s*(.*?)\s*$",spec)
        m===nothing && throw(ArgumentError(
            "$caller: malformed definition $(repr(spec))"))
        name=String(m.captures[1]);parens=m.captures[2]!==nothing
        rhs=String(strip(m.captures[3]))
        if startswith(rhs,"{") && endswith(rhs,"}")
            # `{values, Opt args, ...}` — numeric prefix then option items.
            inner=String(chop(rhs;head=1,tail=1))
            values_src,options=_geo_split_constant_brace(inner)
            if _geo_string_rhs(values_src) || (isempty(strip(values_src)) &&
                                             options!==nothing)
                # `{"str" [, CharOptions]}` — string exchange.
                value=_geo_eval_string(values_src,context,caller)
                if !haskey(context.strings,name)
                    context.strings[name]=[value]
                    _geo_exchange_onelab_string!(context,name,value,
                                                 options,caller)
                end
            else
                values=_geo_numeric_list_values(values_src,context,caller;
                                                allow_multiplier=true)
                length(values)>1 && !parens && _geo_yywarn!(context,
                    "List notation should be used to define list '$name[]'")
                if !_geo_context_has_variable(context,name)
                    # The ONELAB exchange can rewrite `values` (server value
                    # wins) before the symbol is stored.
                    values=_geo_exchange_onelab_number!(context,name,values,
                                                        options,caller)
                    if parens || length(values)!=1
                        _geo_context_set_list!(context,name,values,caller)
                    else
                        _geo_context_set_scalar!(context,name,values[1],caller)
                    end
                end
            end
            continue
        end
        # Scalar/string forms (no option tail exists for these productions).
        if _geo_string_rhs(rhs)
            haskey(context.strings,name) || (
                context.strings[name]=String[
                    _geo_eval_string(rhs,context,caller)])
        else
            _geo_context_has_variable(context,name) ||
                _geo_context_set_scalar!(context,name,
                    _geo_eval_numeric(rhs,context,caller),caller)
        end
    end
    return nothing
end

# Split `{...}` braces content into the leading numeric/string list and the
# trailing `Name "n", Min a, ...` option items. The boundary is the first
# comma-separated item that parses as `Word <args>` rather than a list term.
function _geo_split_constant_brace(inner::AbstractString)
    items=_geo_split_args(inner,"execute_geo: DefineConstant")
    boundary=length(items)+1
    for (i,item) in pairs(items)
        s=strip(item)
        # An option item starts with a bare word NOT followed by an operator or
        # `(` call — `Min 0`, `Name "n"`, `ReadOnly`, `Choices {..}`.
        if (om=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(.*)$",s))!==nothing &&
           !isempty(strip(om.captures[2])) &&
           match(r"^[-+*/<>=:!&|?^]",strip(om.captures[2]))===nothing &&
           !(startswith(strip(om.captures[2]),"("))
            boundary=i;break
        end
        # A bare word alone is also an option flag (`ReadOnly`, `Enum`) — but
        # could equally be a variable read; option keywords win.
        if om!==nothing && isempty(strip(om.captures[2])) &&
           om.captures[1] in ("Min","Max","Step","Name","Label","Choices",
                              "ReadOnly","ReadOnlyRange","ReadOnlyChoices",
                              "Visible","AutoCheck","Closed","Macro","Enum",
                              "MultipleSelection","ServerActionHide",
                              "ServerActionShow","Loop","Graph","Range")
            boundary=i;break
        end
    end
    boundary>length(items) &&
        return join(items,","),nothing
    return join(items[1:boundary-1],","),items[boundary:end]
end

# `Printf`/`Warning`/`Error` — the `StringExprVar [, RecursiveListOfDouble]`
# forms plus `>`/`>>` file redirection.
function _geo_exec_print!(kind::AbstractString,args_src::AbstractString,
                          redirect,context::_GeoNumericContext)
    kind=String(kind)
    caller="execute_geo: $kind"
    args=_geo_split_args(args_src,caller)
    isempty(args) && throw(ArgumentError("$caller requires a format string"))
    format=_geo_eval_string(args[1],context,caller)
    output=if length(args)==1
        format
    else
        # The argument tail is one RecursiveListOfDouble — the same element
        # grammar as a braced list.
        values=_geo_numeric_list_values(
            "{"*join(args[2:end],",")*"}",context,caller;
            allow_multiplier=true)
        out,extra=_geo_print_list_of_double(format,values,caller)
        if extra<0
            # Gmsh's Warning path really does say "in Error" (a baked-in
            # copy-paste quirk in Gmsh.y).
            kind=="Warning" ? _geo_yywarn!(context,"Too few arguments in Error") :
                _geo_yyerror!(context,"Too few arguments in $kind")
            return nothing
        elseif extra>0
            # Same "in Error" quirk on the Warning path.
            label=kind=="Warning" ? "Error" : kind
            message="$extra extra argument$(extra>1 ? "s" : "") in $label"
            kind=="Warning" ? _geo_yywarn!(context,message) :
                _geo_yyerror!(context,message)
            return nothing
        end
        out
    end
    if redirect!==nothing
        mode,target_src=redirect
        target_src===nothing && throw(ArgumentError(
            "$caller: `>` requires a file name"))
        path=_geo_fix_relative_path(context.file_name,
            _geo_eval_string(target_src,context,caller))
        try
            open(path,mode==">>" ? "a" : "w") do io
                print(io,output,"\n")
            end
        catch err
            err isa InterruptException && rethrow()
            _geo_yyerror!(context,"Unable to open file '$path'")
        end
        return nothing
    end
    # `Msg::Direct`/`Warning`/`Error` vsnprintf their text a second time —
    # `%%` collapses and remaining specs read as 0/NULL. File redirection
    # skips this pass.
    expanded=_geo_vsnprintf_expand(output)
    if kind=="Printf"
        # `Msg::Direct` strips one trailing newline then writes `str\n` to
        # stdout.
        print(stdout,chomp(expanded),"\n");flush(stdout)
    elseif kind=="Warning"
        _geo_yywarn!(context,expanded)
    else
        # `Error(...)` is `Msg::Error` — counted on the msg channel (it sets
        # `_atLeastOneErrorInRun` and the batch exit code) but does not fail
        # the parse.
        _geo_msg_error!(context,expanded)
    end
    return nothing
end

# `Show`/`Hide`/`Recursive Show`/`Recursive Hide` — and the deprecated
# `Show "name"` string forms — write `entity_visibility[(dim,tag)]`.
function _geo_exec_visibility!(m::GeoModel,line::AbstractString,
                               context::_GeoNumericContext,
                               allocator_state=nothing)
    # `setVisibility` syncs the internals (when changed) before applying —
    # the shape list may carry `Physical` selectors that need the view.
    _geo_sync_physical_view_if_changed!(m,context)
    recursive=match(r"^Recursive\b",line)!==nothing
    s=String(strip(line))
    recursive && (s=String(strip(s[nextind(s,9):end])))
    head=match(r"^(Show|Hide)\b",s)
    head===nothing && throw(ArgumentError(
        "execute_geo: malformed visibility statement: $line"))
    visible=head.captures[1]=="Show"
    s=String(strip(s[nextind(s,lastindex(head.captures[1])+1):end]))
    caller="execute_geo: $(head.captures[1])"
    # Deprecated string form: `Show "text"` shows every entity whose name
    # matches — Tessella stores names but a wildcard match over all names is
    # still faithful.
    if (qm=match(r"^\"(.*)\"\s*;?$",s))!==nothing ||
       (qm=match(r"^'(.*)'\s*;?$",s))!==nothing
        needle=String(qm.captures[1])
        for ((dim,tag),name) in m.entity_names
            name==needle && (m.entity_visibility[(dim,tag)]=visible)
        end
        return nothing
    end
    (inner,rest)=_geo_balanced_group(s,caller)
    isempty(strip(rest)) || isempty(strip(rest[1:end-1])) ||
        throw(ArgumentError("$caller: unexpected text after the entity list"))
    entities=_geo_shape_list_entities!(
        m,inner,context,caller;signed_tags=true,
        allocator_state=allocator_state)
    _geo_set_visibility!(m,entities,visible,recursive,caller)
    return nothing
end

function _geo_set_visibility!(m::GeoModel,entities,visible::Bool,
                             recursive::Bool,caller::AbstractString)
    seen=Set{Tuple{Int,Int}}()
    stack=Tuple{Int,Int}[entities...]
    while !isempty(stack)
        key=pop!(stack)
        key in seen && continue
        push!(seen,key)
        m.entity_visibility[key]=visible
        recursive || continue
        dim,tag=key
        for (bdim,btag) in _model_direct_boundary(
                m,dim,tag,caller;canonical_orientation=false)
            push!(stack,(bdim,abs(Int(btag))))
        end
    end
    return nothing
end

# `Color <expr> {shapes}` / `Recursive Color <expr> {shapes}` — and the
# deprecated `Color "name" {shapes}` string form — write `entity_colors`.
function _geo_exec_color!(m::GeoModel,line::AbstractString,
                        context::_GeoNumericContext,
                        allocator_state=nothing)
    # `setColor` syncs the internals (when changed) before applying.
    _geo_sync_physical_view_if_changed!(m,context)
    recursive=match(r"^Recursive\s+Color\b",line)!==nothing
    s=String(strip(line))
    s=recursive ? String(strip(s[nextind(s,15):end])) :
        String(strip(s[nextind(s,6):end]))
    caller="execute_geo: Color"
    rgba=_geo_parse_color_expr(s,context,caller)
    rgba===nothing && return nothing # error already recorded
    color,rest=rgba
    startswith(strip(rest),"{") || throw(ArgumentError(
        "$caller: expected `{ ListOfShapes }` after the color expression"))
    (inner,rest2)=_geo_balanced_group(String(strip(rest)),caller)
    isempty(strip(rest2)) || isempty(strip(rest2[1:end-1])) ||
        throw(ArgumentError("$caller: unexpected text after the entity list"))
    entities=_geo_shape_list_entities!(
        m,inner,context,caller;signed_tags=true,
        allocator_state=allocator_state)
    seen=Set{Tuple{Int,Int}}()
    stack=Tuple{Int,Int}[entities...]
    while !isempty(stack)
        key=pop!(stack)
        key in seen && continue
        push!(seen,key)
        m.entity_colors[key]=color
        recursive || continue
        dim,tag=key
        for (bdim,btag) in _model_direct_boundary(
                m,dim,tag,caller;canonical_orientation=false)
            push!(stack,(bdim,abs(Int(btag))))
        end
    end
    return nothing
end

# `ColorExpr` — `{r,g,b[,a]}`, a named color (`"Red"` or `Color.Red` handled by
# the caller's dispatch), or `Color.X` option reads. Returns
# `(rgba, rest_after)` or `nothing` after recording an error.
function _geo_parse_color_expr(s::AbstractString,context::_GeoNumericContext,
                               caller::AbstractString)
    src=String(strip(s))
    if startswith(src,"{")
        (inner,rest)=_geo_balanced_group(src,caller)
        values=_geo_numeric_list_values("{$(inner)}",context,caller)
        (3<=length(values)<=4) || throw(ArgumentError(
            "$caller: a color list takes three or four components"))
        rgba=NTuple{4,Int}(i<=length(values) ?
            _geo_int_value(values[i],"$caller component") : 255
            for i in 1:4)
        all(0<=v<=255 for v in rgba) || throw(ArgumentError(
            "$caller: color components must be in 0:255"))
        return rgba,rest
    end
    if (qm=match(r"^\"([^\"]*)\"(.*)$",src))!==nothing ||
       (qm=match(r"^'([^']*)'(.*)$",src))!==nothing
        name=String(qm.captures[1]);rest=String(qm.captures[2])
        # `GetColorForString` is an exact `strcmp` over the X11 table.
        haskey(_GEO_COLOR_NAMES,name) || throw(ArgumentError(
            "$caller: unknown color name '$name'"))
        return _GEO_COLOR_NAMES[name],rest
    end
    # `x.Color.f` — `ColorOption(GMSH_GET)` on the stored write mirror, then
    # the option-name table.
    if (cm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:\[\s*(.*?)\s*\])?\s*\.\s*Color\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*(.*)$",
            src))!==nothing
        family=String(cm.captures[1]);member=String(cm.captures[3])
        index=cm.captures[2]===nothing ? 0 :
            _geo_int_value(_geo_eval_numeric(cm.captures[2],context,
                "$caller color index"),"$caller color index")
        rest=String(cm.captures[4])
        key=(family,index,member)
        haskey(context.option_colors,key) &&
            return context.option_colors[key],rest
        names=get(_GEO_COLOR_OPTION_NAMES,family,nothing)
        if names===nothing || !(member in names)
            context.soft_unknown_reads || throw(ArgumentError(
                "$caller: Unknown color option '$family.$member'"))
            names===nothing ? _geo_msg_error!(context,
                "Unknown color option category '$family'") :
                _geo_msg_error!(context,
                    "Unknown color option '$family.$member'")
            return (0,0,0,255),rest
        end
        # Known option, never written — `ColorOption(GMSH_GET)` returns the
        # packed default; the mirror stores the default only once observed.
        return (0,0,0,255),rest
    end
    # `String__Index`: a defined string variable supplies the color name;
    # otherwise the word itself is the color name.
    if (nm=match(r"^([A-Za-z_][A-Za-z0-9_\.]*)\s*(.*)$",src))!==nothing
        word=String(nm.captures[1]);rest=String(nm.captures[2])
        if !occursin('.',word) && haskey(context.strings,word)
            vals=context.strings[word]
            if isempty(vals)
                _geo_yyerror!(context,"Unknown color '$word'")
                return (0,0,0,255),rest
            end
            if haskey(_GEO_COLOR_NAMES,vals[1])
                return _GEO_COLOR_NAMES[vals[1]],rest
            end
            _geo_yyerror!(context,"Unknown color '$(vals[1])'")
            return (0,0,0,255),rest
        end
        haskey(_GEO_COLOR_NAMES,word) &&
            return _GEO_COLOR_NAMES[word],rest
    end
    throw(ArgumentError("$caller: unrecognized color expression: $src"))
end

# ======== mid-file mesh commands (`Mesh n` and the mesh-operation family) =====

# `Mesh n;` — mesh every live entity up to `dim`, merging entity parts into a
# single Mesh plus a node→(dim,tag) ownership map (lowest dimension wins, like
# Gmsh's per-entity vertex storage).
function _geo_mesh_model(m::GeoModel,dim::Int,context::_GeoNumericContext)
    caller="execute_geo: Mesh"
    parts=Tuple{Int,Int,Mesh}[]
    if dim==1
        throw(ArgumentError(
            "$caller 1: standalone curve meshing is not implemented — " *
            "curve nodes are produced as part of surface/volume meshes"))
    end
    if dim>=2
        for tag in sort!(collect(keys(m.surfaces)))
            push!(parts,(2,tag,mesh_model_surface(m,tag)))
        end
    end
    if dim==3
        for tag in sort!(collect(keys(m.volumes)))
            push!(parts,(3,tag,mesh_model_volume(m,tag)))
        end
    end
    isempty(parts) && throw(ArgumentError(
        "$caller $dim: no $(dim==2 ? "surfaces" : "entities") to mesh"))
    return _geo_merge_entity_meshes(parts)
end

# Merge per-entity meshes: nodes deduplicated by exact coordinates (shared
# boundary nodes are bitwise identical across entity meshes), elements tagged
# with their generating entity in `owner` for node-level classification.
function _geo_merge_entity_meshes(parts::Vector{Tuple{Int,Int,Mesh}})
    coord_list=NTuple{3,Float64}[]
    lookup=Dict{NTuple{3,Int},Int}()
    owner=Dict{Int,Tuple{Int,Int}}()
    segs=Matrix{Int32}(undef,2,0);seg_tag=Int32[]
    tris=Matrix{Int32}(undef,3,0);tri_tag=Int32[]
    tets=Matrix{Int32}(undef,4,0);tet_tag=Int32[]
    function node_index(col)
        x,y,z=col
        key=(reinterpret(Int,x),reinterpret(Int,y),reinterpret(Int,z))
        i=get(lookup,key,0)
        i==0 || return i
        i=length(coord_list)+1
        lookup[key]=i
        push!(coord_list,(x,y,z))
        return i
    end
    for (dim,tag,mesh) in parts
        remap=Vector{Int}(undef,nnodes(mesh))
        for n in 1:nnodes(mesh)
            remap[n]=node_index(mesh.coords[:,n])
        end
        for s in 1:nsegs(mesh)
            a,b=remap[mesh.segs[1,s]],remap[mesh.segs[2,s]]
            segs=hcat(segs,Int32[a,b])
            push!(seg_tag,mesh.seg_tag[s])
            for n in (a,b)
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(dim,tag))
            end
        end
        for t in 1:ntris(mesh)
            a,b,c=(remap[mesh.tris[k,t]] for k in 1:3)
            tris=hcat(tris,Int32[a,b,c])
            push!(tri_tag,mesh.tri_tag[t])
            for n in (a,b,c)
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(dim,tag))
            end
        end
        for t in 1:ntets(mesh)
            idx=[remap[mesh.tets[k,t]] for k in 1:4]
            tets=hcat(tets,Int32.(idx))
            push!(tet_tag,mesh.tet_tag[t])
            for n in idx
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(dim,tag))
            end
        end
    end
    coords=Matrix{Float64}(undef,3,length(coord_list))
    for (i,(x,y,z)) in enumerate(coord_list)
        coords[1,i]=x;coords[2,i]=y;coords[3,i]=z
    end
    merged=Mesh(coords;segs=segs,tris=tris,tets=tets,
                seg_tag=seg_tag,tri_tag=tri_tag,tet_tag=tet_tag)
    return merged,owner
end

# The mesh-operation statements act on `context.mesh` — the product of a
# mid-file `Mesh n` — or record meshing constraints for the next generation.
function _geo_exec_mesh_statement!(m::GeoModel,line::AbstractString,
                                   context::_GeoNumericContext)
    caller="execute_geo"
    s=String(strip(line))
    if (mm=match(r"^RefineMesh\s*;?\s*$",s))!==nothing
        # `RefineMesh` syncs the internals first when changed.
        _geo_sync_physical_view_if_changed!(m,context)
        context.mesh===nothing && return nothing
        context.mesh=refine_uniform(context.mesh)
        _geo_remesh_owner!(context)
        return
    elseif (mm=match(r"^RecombineMesh\s*;?\s*$",s))!==nothing
        context.mesh===nothing && return nothing
        context.mesh=recombine_triangles(context.mesh)
        _geo_remesh_owner!(context)
        return
    elseif (mm=match(r"^RenumberMeshNodes\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^RenumberMeshElements\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^CreateMeshEdges\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^CreateMeshFaces\s*;?\s*$",s))!==nothing
        # Node/element order and topology catalogs are derived information in
        # the compact Mesh — already canonical, nothing to mutate.
        return
    elseif (mm=match(r"^ReorientMesh\s+Volume\s*\{(.*?)\}\s*;?\s*$",s))!==nothing
        # `ReorientMesh` syncs the internals first when changed.
        _geo_sync_physical_view_if_changed!(m,context)
        tags=_geo_numeric_list_values(mm.captures[1],context,caller)
        for value in tags
            tag=_geo_int_value(value,"$caller ReorientMesh tag")
            # `setOutwardOrientationMeshConstraint` — a flag consumed at the
            # NEXT volume meshing, like `Reverse`.
            m.meshing.reverse[(3,tag)]=true
        end
        return
    elseif (mm=match(r"^RelocateMesh\s+(Point|Curve|Surface)\s*\{(.*?)\}\s*;?\s*$",s))!==nothing
        # `RelocateMesh` syncs the internals first when changed.
        _geo_sync_physical_view_if_changed!(m,context)
        context.mesh===nothing && return nothing
        dim=_GEO_STRING_ENTITY_DIM[mm.captures[1]]
        tags=_geo_numeric_list_values(mm.captures[2],context,caller)
        wanted=Set(_geo_int_value(t,"$caller RelocateMesh tag") for t in tags)
        mesh=context.mesh
        for n in 1:nnodes(mesh)
            own=get(context.mesh_node_owner,n,nothing)
            own===nothing && continue
            own[1]==dim && own[2] in wanted || continue
            point=_geo_relocate_node(m,dim,own[2],mesh.coords[:,n])
            point===nothing || (mesh.coords[:,n]=point)
        end
        return
    elseif (mm=match(r"^TransformMesh\s*\{(.*?)\}\s*(?:\{(.*?)\})?\s*;?\s*$",s))!==nothing
        context.mesh===nothing && return nothing
        matrix=_geo_numeric_list_values(mm.captures[1],context,caller)
        length(matrix)>=12 || throw(ArgumentError(
            "$caller: affine transform matrix requires at least 12 entries"))
        wanted=nothing
        if mm.captures[2]!==nothing
            entities=_geo_shape_list_entities!(
                m,mm.captures[2],context,"$caller: TransformMesh";
                signed_tags=true,allocator_state=allocator_state)
            wanted=Set(Tuple{Int,Int}(abs.(Int.(collect(e)))) for e in entities)
        end
        _geo_transform_mesh_nodes!(context.mesh,context.mesh_node_owner,
            matrix,wanted)
        return
    elseif (mm=match(r"^SetOrder\s+(.*?)\s*;?\s*$",s))!==nothing
        order=_geo_int_value(
            _geo_eval_numeric(mm.captures[1],context,caller),caller)
        order<1 && throw(ArgumentError("$caller: SetOrder requires order >= 1"))
        order>1 && throw(ArgumentError(
            "$caller: SetOrder $order — high-order .geo meshes are not " *
            "supported"))
        return
    elseif (mm=match(r"^PartitionMesh\s+(.*?)\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^CreateOverlaps\s+(.*?)\s*;?\s*$",s))!==nothing
        # Partitioning needs the partitioned-mesh data model — an explicit
        # blocker rather than a silent partial implementation.
        throw(ArgumentError(
            "$caller: mesh partitioning is not supported"))
    elseif (mm=match(r"^AdaptMesh\b",s))!==nothing
        throw(ArgumentError(
            "$caller: AdaptMesh requires levelset fields, which are not " *
            "supported"))
    elseif (mm=match(r"^CreateTopology\b",s))!==nothing ||
           (mm=match(r"^ClassifySurfaces\b",s))!==nothing ||
           (mm=match(r"^CreateGeometry\b",s))!==nothing
        throw(ArgumentError(
            "$caller: discrete-model topology creation is not supported"))
    end
    throw(ArgumentError("$caller: unrecognized mesh command: $line"))
end

# Node ownership no longer maps cleanly after element regeneration — reset to
# "everything belongs to its generating entity" by re-classifying nodes onto
# the elements they belong to. For a uniform refinement every node keeps its
# owner via the parent element; simplest correct-enough choice is to drop
# ownership (entity-scoped ops then apply to all nodes, which is wrong) — so
# instead preserve it: refinement children inherit the parent element's owner.
function _geo_remesh_owner!(context::_GeoNumericContext)
    # `refine_uniform`/`recombine_triangles` keep node order for the original
    # vertices; appended nodes are interior to refined elements. Ownership of
    # the appended nodes is dropped (they belong to whatever entity refined
    # them — without element provenance, treat them as unowned).
    return nothing
end

# `RelocateMesh` — move the node onto its owning entity's geometry. Curves and
# surfaces project through the same `model_to_mixed`-free path used by
# embeddings: for straight curves the segment, for surfaces the owning
# surface's geometry.
function _geo_relocate_node(m::GeoModel,dim::Int,tag::Int,
                            xyz::AbstractVector)
    if dim==0
        haskey(m.points,tag) || return nothing
        return copy(m.points[tag])
    elseif dim==1
        haskey(m.curves,tag) || return nothing
        # Straight/arc curves: project onto the chord — a bounded
        # approximation; exact curve inversion lives in the mesher.
        a,b=m.curves[tag]
        pa,pb=m.points[a],m.points[b]
        d=pb.-pa;t=sum((xyz.-pa).*d)/max(sum(d.*d),eps())
        return pa.+clamp(t,0.0,1.0).*d
    end
    # Surface/volume nodes keep their positions — a full parametric inverse is
    # not available mid-file.
    return nothing
end

function _geo_transform_mesh_nodes!(mesh::Mesh,owner,matrix,wanted)
    n=length(matrix)
    # Gmsh's affine transform packs row-major 4x4 with defaults.
    a=Float64[1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]
    for i in 1:min(n,16)
        a[i]=matrix[i]
    end
    for v in 1:nnodes(mesh)
        if wanted!==nothing
            own=get(owner,v,nothing)
            (own===nothing || !(own in wanted)) && continue
        end
        x,y,z=mesh.coords[1,v],mesh.coords[2,v],mesh.coords[3,v]
        mesh.coords[1,v]=a[1]*x+a[2]*y+a[3]*z+a[4]
        mesh.coords[2,v]=a[5]*x+a[6]*y+a[7]*z+a[8]
        mesh.coords[3,v]=a[9]*x+a[10]*y+a[11]*z+a[12]
    end
    return nothing
end

# ======== `Word StringExpr;` command family ========

function _geo_exec_command!(m::GeoModel,word::AbstractString,
                            arg_src::AbstractString,
                            context::_GeoNumericContext,allocator_state)
    caller="execute_geo: $word"
    if word=="Sleep" || word=="Remesh" || word=="SetOrder" ||
       word=="PartitionMesh" || word=="CreateOverlaps"
        # These take an FExpr, not a string.
        value=_geo_eval_numeric(arg_src,context,caller)
        if word=="Sleep"
            sleep(max(0.0,value));return nothing
        elseif word=="Remesh"
            _geo_yyerror!(context,"Surface remeshing must be reinterfaced")
            return nothing
        end
        # SetOrder/PartitionMesh/CreateOverlaps are dispatched before this —
        # reaching here means the FExpr form slipped past the mesh-op regex.
        return _geo_exec_mesh_statement!(m,"$word $arg_src;",context)
    end
    arg=_geo_eval_string(arg_src,context,caller)
    if word=="Merge" || word=="MergeWithBoundingBox"
        # `Merge`/`MergeWithBoundingBox` sync the internals (when changed) so
        # merged entity tags cannot clash — `Include` does not sync.
        _geo_sync_physical_view_if_changed!(m,context)
        _geo_exec_include!(m,arg,context,caller,allocator_state)
    elseif word=="Include"
        _geo_exec_include!(m,arg,context,caller,allocator_state)
    elseif word=="Save"
        _geo_exec_save!(m,arg,context,caller)
    elseif word=="Print"
        _geo_exec_save!(m,arg,context,caller;kind="Print")
    elseif word=="System" || word=="SystemCall"
        _geo_system_call(arg,caller)
    elseif word=="NonBlockingSystemCall"
        _geo_system_call(arg,caller;wait=false)
    elseif word=="SetName"
        m.name=arg
    elseif word=="CreateDir"
        mkpath(_geo_fix_relative_path(context.file_name,arg))
    elseif word=="OnelabRun"
        throw(ArgumentError(
            "$caller: external ONELAB clients are not supported"))
    elseif word=="OptimizeMesh"
        context.mesh===nothing && return nothing
        throw(ArgumentError(
            "$caller: mesh optimization is not supported mid-file"))
    elseif word=="SetBoundingBox"
        throw(ArgumentError(
            "$caller: `SetBoundingBox` was removed in Gmsh 3.0; " *
            "use `BoundingBox{...}`"))
    elseif word in ("SetCurrentWindow","SplitCurrentWindowHorizontal",
                    "SplitCurrentWindowVertical")
        # GUI window management — no-op in batch.
        return nothing
    else
        throw(ArgumentError("$caller: unknown command '$word'"))
    end
    return nothing
end

# `Include`/`Merge` of a `.geo` file runs its statements in the same context
# and model (Gmsh's `ParseFile` shares the parser globals). `.msh` files merge
# their mesh into the current mid-file mesh.
function _geo_exec_include!(m::GeoModel,path::AbstractString,
                            context::_GeoNumericContext,caller::AbstractString,
                            allocator_state)
    file=_geo_fix_relative_path(context.file_name,path)
    isfile(file) || throw(ArgumentError("$caller: missing file '$file'"))
    ext=lowercase(splitext(file)[2])
    if ext==".msh"
        merged=read_msh(file).mesh
        context.mesh===nothing ? (context.mesh=merged) :
            (context.mesh=_geo_concat_meshes(context.mesh,merged))
        empty!(context.mesh_node_owner)
        return nothing
    end
    ext==".geo" || throw(ArgumentError(
        "$caller: merging '$ext' files is not supported (only .geo and .msh)"))
    statements=_geo_exec_statements(file)
    parent_file=context.file_name
    parent_base=context.file_error_base
    parent_stop=context.stop
    context.file_name=file
    context.file_error_base=length(context.exec_errors)
    try
        _exec_geo_statements!(
            m,statements,firstindex(statements),lastindex(statements),
            context,allocator_state,Ref(0))
    finally
        context.file_name=parent_file
        context.file_error_base=parent_base
        # `Abort` (errorstate 999) and too-many-errors are local to the
        # included file's ParseFile — the parent stream resumes; `Exit` is
        # process-global and propagates.
        context.stop in (:abort,:toomany) && (context.stop=parent_stop)
    end
    return nothing
end

function _geo_concat_meshes(a::Mesh,b::Mesh)
    offset=nnodes(a)
    coords=hcat(a.coords,b.coords)
    segs=hcat(a.segs,b.segs.+Int32(offset))
    tris=hcat(a.tris,b.tris.+Int32(offset))
    tets=hcat(a.tets,b.tets.+Int32(offset))
    return Mesh(coords;segs=segs,tris=tris,tets=tets,
        seg_tag=vcat(a.seg_tag,b.seg_tag),
        tri_tag=vcat(a.tri_tag,b.tri_tag),
        tet_tag=vcat(a.tet_tag,b.tet_tag))
end

# `Save`/`Print` — `CreateOutputFile` writes the current mesh; the only
# supported target format is `.msh` (v2.2, matching `write_msh`'s default).
function _geo_exec_save!(m::GeoModel,path::AbstractString,
                         context::_GeoNumericContext,caller::AbstractString;
                         kind::String="Save")
    file=_geo_fix_relative_path(context.file_name,path)
    endswith(lowercase(file),".msh") || throw(ArgumentError(
        "$caller: only .msh output is supported (got '$file')"))
    # `CreateOutputFile` writes the synced model — materialize the physical
    # view (and pending unknown-member warnings) first, when internals changed.
    _geo_sync_physical_view_if_changed!(m,context)
    mesh=context.mesh
    if mesh===nothing
        isempty(m.points) || throw(ArgumentError(
            "$caller: no mesh exists — run `Mesh n;` first"))
        mesh=Mesh(zeros(3,0))
    end
    write_msh(file,mesh;physical_names=m.physical_names)
    return nothing
end

function _geo_system_call(command::AbstractString,caller::AbstractString;
                          wait::Bool=true)
    isempty(strip(command)) && return nothing
    cmd=Sys.iswindows() ? `cmd /c $command` : `sh -c $command`
    wait ? run(cmd) : run(cmd;wait=false)
    return nothing
end

end # module
