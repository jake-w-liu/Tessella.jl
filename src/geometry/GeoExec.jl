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
Translate/Rotate/Affine periodic straight or curved curves, explicit-volume
planar boundary
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

import ..Model
import ..Elements

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
using ..Model: _model_boundary, _model_points_of, _model_direct_boundary,
               _automatic_context_field
using ..Model: _model_entity_dictionary, _model_entity_known, remove_embedded!
using ..Model: model_entity, model_entities, model_entities_in_bounding_box
using ..Model: model_physical_groups, _physical_live_members
using ..Model: model_entity_color, model_parametrization_bounds, model_normal
using ..Model: model_reparametrize_on_surface
using ..Model: model_parametrization, model_value
using ..Model: _model_entity_bounding_box, _model_bounds_union, _occ_geometry
using ..Model: _surface_type, _model_volume_boundary_surfaces
using ..Model: _geo_delete_entities!, _geo_reset_model_geometry!
using ..Model: _tag, _alloc_tag!, _alloc_curve_loop_tag,
    _alloc_surface_loop_tag, _arc_stored_normal, _occ_cross,
    _geo_newreg_alloc!
using ..Model: mesh_model_surface, mesh_model_volume
using ..Model: create_topology!, classify_surfaces!, create_geometry!,
    add_homology_request!, compute_homology!, _discrete_mesh_records_model,
    _curve_type, _periodic_curve_point, _transfinite_parameters,
    model_dimension
using ..MeshTypes: Mesh
using ..IO: read_geo_params, _GeoNumericContext, _GeoGeometryParts, _geo_eval_numeric
using ..IO: _geo_exchange_onelab_number!, _geo_exchange_onelab_string!
using ..IO: _geo_split_list, _geo_split_range, _geo_range_count
using ..IO: _geo_numeric_list_terms, _geo_numeric_list_values
using ..IO: _geo_signed_gmsh_int_value
using ..IO: _geo_store_option_number!
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
using ..IO: _geo_matching_delim, _GEO_COLOR_NAMES, _geo_braced_end
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
using ..IO: GeoFieldSpec, GeoParams, _GEO_FIELD_KINDS,
            _geo_expression_field_tags
using ..SizeField: _build_geo_field, AbstractField, AbstractSizeField
import ..SizeField: field_value
using ..Model: _ModelMesh1DOptions, _model_mesh_lc, _model_mesh_bbox,
               _model_mesh_dim01!, _model_mesh_dim01_grade!,
               _model_mesh_dim01_parts!, _model_point_mesh_parts,
               _model_curve_mesh_parts,
               _model_curve_part_point, _model_planar_surface_mesh,
               _model_curve_param_bounds,
               _model_projection_physical_tags, _model_projection_legacy_tag,
               _model_projection_surface_curves, _model_mesh_part_node_entities,
               _extrude_volume_pass, _extrude_gate, _extrude_is_nonew,
               _extrude_is_quadtri,
               _extrude_nonew_scope
using ..Model: _extrude_level_us, _extrude_curve_nodes,
               _extrude_curve_params_stage, _extrude_bl_plane_normal,
               _extrude_bl_normal_field, _extrude_bl_kernel,
               _extrude_bl_refresh!,
               _duplicate_curve!, _fresh_point_copy!, _coherence_merge!,
               _geo_derived_loop_tag
using ..MeshTypes: nnodes, nsegs, ntris, ntets
using ..Refine: refine_uniform
import ..Refine
using ..Recombine: recombine_triangles
using ..IO: write_msh
using ..Elements: read_mixed_msh, write_mixed_msh, MixedMesh, ElementBlock, msh_spec,
                  SpecialElementBlock, MixedElementBlock,
                  msh_dimension, msh_family, msh_num_nodes
using ..Model: add_discrete_entity!, DiscreteEntity, _record_append_node!,
               _record_append_element!, _model_fresh_element_tag
using ..MeshTypes: validate
using ..Optimize: smooth_optimize, _laplacian_smooth_tri_cache
using ..IO: _geo_gmsh_number
using ..Transform: _affine_coordinate, _periodic_affine_3x4
using LinearAlgebra: norm, svd

# Separate emitted geometry entity parts from actual native products. Dense
# carrier rows and native identities are explicit; coordinates never pick IDs.
function _geo_geometry_blocks(mesh::Mesh)
    return ((1,mesh.segs,mesh.seg_tag),(2,mesh.tris,mesh.tri_tag),(4,mesh.tets,mesh.tet_tag))
end
function _geo_geometry_blocks(mesh::MixedMesh)
    return [(Int(block.msh),block.nodes,block.tags) for block in mesh.blocks]
end

function _geo_geometry_part_coordinates(part::Mesh,coordinates)
    return Mesh(coordinates;segs=part.segs,tris=part.tris,tets=part.tets,
        seg_tag=part.seg_tag,tri_tag=part.tri_tag,tet_tag=part.tet_tag)
end
function _geo_geometry_part_coordinates(part::MixedMesh,coordinates)
    part.entity_data===nothing || throw(ArgumentError(
        "execute_geo: geometry entity parts must not claim actual native metadata"))
    return MixedMesh(coordinates,part.blocks;physical_names=part.physical_names,
        elementary_entities=part.elementary_entities,periodic_links=part.periodic_links,
        ancillary_sections=part.ancillary_sections,data_sections=part.data_sections,
        partition_data=part.partition_data)
end

function _geo_geometry_receipt(carrier,parts,rows,node_tags,caller)
    length(parts)==length(rows) || throw(ArgumentError("$caller: geometry part rows do not align"))
    length(node_tags)==nnodes(carrier) && allunique(node_tags) || throw(ArgumentError(
        "$caller: geometry carrier identities do not align"))
    keys=Set{Tuple{Int,Int}}()
    emitted=Dict{Tuple,Int}()
    for ((dim,tag,part),mapping) in zip(parts,rows)
        owner=(dim,tag)
        owner in keys && throw(ArgumentError("$caller: repeated geometry entity part $owner"))
        push!(keys,owner)
        length(mapping)==nnodes(part) && allunique(mapping) &&
            all(row->0<row<=nnodes(carrier),mapping) || throw(ArgumentError(
                "$caller: geometry entity $owner has invalid carrier rows"))
        isequal(part.coords,carrier.coords[:,mapping]) || throw(ArgumentError(
            "$caller: geometry entity $owner disagrees with its carrier coordinates"))
        part isa MixedMesh && part.entity_data!==nothing && throw(ArgumentError(
            "$caller: geometry entity $owner claims actual native metadata"))
        for (msh,cells,_) in _geo_geometry_blocks(part),cell in eachcol(cells)
            msh_dimension(msh)==dim || throw(ArgumentError(
                "$caller: geometry entity $owner contains a foreign dimension"))
            key=(msh,Tuple(mapping[cell]))
            emitted[key]=get(emitted,key,0)+1
        end
    end
    merged=Dict{Tuple,Int}()
    for (msh,cells,_) in _geo_geometry_blocks(carrier),cell in eachcol(cells)
        key=(msh,Tuple(cell));merged[key]=get(merged,key,0)+1
    end
    emitted==merged || throw(ArgumentError(
        "$caller: geometry entity cells do not partition the returned carrier"))
    return _GeoGeometryParts(carrier,parts,rows,copy(node_tags))
end

function _geo_geometry_carrier_tags(authority,carrier,caller)
    carrier===authority.carrier && length(authority.carrier_rows)==nnodes(carrier) ||
        throw(ArgumentError("$caller: geometry carrier lacks its actual identity rows"))
    return authority.mesh.entity_data.external_node_tags[authority.carrier_rows]
end

function _geo_geometry_coordinate_receipt(receipt,carrier,caller)
    receipt===nothing && return nothing
    nnodes(carrier)==nnodes(receipt.carrier) || throw(ArgumentError(
        "$caller: coordinate edit changes geometry carrier row count"))
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
        (dim,tag,_geo_geometry_part_coordinates(part,carrier.coords[:,rows]))
        for ((dim,tag,part),rows) in zip(receipt.parts,receipt.rows)]
    return _geo_geometry_receipt(carrier,parts,receipt.rows,receipt.node_tags,caller)
end

# Reclassification, renumbering and changed cell orientations use real native
# cell tags. An unlabelled proxy cell can retain its old explicitly bound owner;
# ambiguity is rejected before any model/context publication.
function _geo_geometry_rebind_receipt(authority,carrier,previous,plan,caller)
    carrier===authority.mesh && return nothing
    tags=_geo_geometry_carrier_tags(authority,carrier,caller)
    native=authority.mesh.entity_data
    cell_owners=Dict{UInt64,Tuple{Int,Int}}()
    for (block,owners,ids) in zip(authority.mesh.blocks,native.block_entities,native.external_element_tags)
        for (owner,id) in zip(owners,ids)
            cell_owners[id]=(msh_dimension(block.msh),Int(owner))
        end
    end
    node_map=hasproperty(plan,:node_tag_map) ? plan.node_tag_map : nothing
    mapped(tag)=node_map===nothing ? tag : get(node_map,tag,nothing)
    proxy=Dict{Tuple,Union{Nothing,Tuple{Int,Int}}}()
    preferred=Dict{Tuple{Int,Int},Vector{UInt64}}()
    previous_kind=Dict{Tuple{Int,Int},Bool}()
    if previous!==nothing
        for ((dim,tag,part),rows) in zip(previous.parts,previous.rows)
            owner=(dim,tag);previous_kind[owner]=part isa Mesh
            old=mapped.(previous.node_tags[rows])
            preferred[owner]=UInt64[value for value in old if value!==nothing]
            for (msh,cells,_) in _geo_geometry_blocks(part),cell in eachcol(cells)
                values=old[cell];any(isnothing,values) && continue
                key=(msh,Tuple(sort!(UInt64[value for value in values])))
                existing=get(proxy,key,owner)
                proxy[key]=existing==owner ? owner : nothing
            end
        end
    end
    geometry_blocks=collect(_geo_geometry_blocks(carrier))
    length(authority.carrier_cells)==length(geometry_blocks) || throw(ArgumentError(
        "$caller: geometry carrier lacks its explicit cell identity blocks"))
    selected=Dict{Tuple{Int,Int},Vector{Vector{Int}}}()
    needed=Dict{Tuple{Int,Int},Set{Int32}}()
    for (bi,(msh,cells,_)) in enumerate(geometry_blocks)
        ids=authority.carrier_cells[bi]
        length(ids)==size(cells,2) || throw(ArgumentError(
            "$caller: geometry carrier cell identities do not align"))
        for column in axes(cells,2)
            id=ids[column]
            owner=if id===nothing
                key=(msh,Tuple(sort!(copy(tags[cells[:,column]]))))
                get(proxy,key,nothing)
            else
                get(cell_owners,id,nothing)
            end
            owner!==nothing && owner[1]==msh_dimension(msh) || throw(ArgumentError(
                "$caller: geometry cell lacks its explicitly bound entity owner"))
            blocks=get!(()->[Int[] for _ in geometry_blocks],selected,owner)
            push!(blocks[bi],column)
            union!(get!(()->Set{Int32}(),needed,owner),cells[:,column])
        end
    end
    for (row,actual) in enumerate(authority.carrier_rows)
        dim,tag=native.node_entities[actual];owner=(dim,Int(tag))
        push!(get!(()->Set{Int32}(),needed,owner),Int32(row))
    end
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(tags))
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[];remaps=Vector{Int32}[]
    for owner in sort!(collect(keys(needed)))
        keep=needed[owner];rows=Int32[];seen=Set{Int32}()
        for tag in get(preferred,owner,UInt64[])
            row=get(positions,tag,nothing)
            row===nothing || !(row in keep) || row in seen || (push!(rows,row);push!(seen,row))
        end
        for row in sort!(collect(keep))
            row in seen || (push!(rows,row);push!(seen,row))
        end
        inverse=Dict(row=>Int32(local_row) for (local_row,row) in enumerate(rows))
        indices=get(selected,owner,[Int[] for _ in geometry_blocks])
        coords=carrier.coords[:,rows]
        part=if carrier isa Mesh && get(previous_kind,owner,true)
            cells=ntuple(bi->begin
                nodes=geometry_blocks[bi][2][:,indices[bi]]
                map(row->inverse[row],nodes)
            end,3)
            labels=ntuple(bi->geometry_blocks[bi][3][indices[bi]],3)
            Mesh(coords;segs=cells[1],tris=cells[2],tets=cells[3],
                seg_tag=labels[1],tri_tag=labels[2],tet_tag=labels[3])
        else
            blocks=ElementBlock[]
            for (bi,(msh,cells,labels)) in enumerate(geometry_blocks)
                isempty(indices[bi]) && continue
                local_nodes=map(row->inverse[row],cells[:,indices[bi]])
                push!(blocks,ElementBlock(msh,local_nodes,labels[indices[bi]]))
            end
            MixedMesh(coords,blocks)
        end
        push!(parts,(owner[1],owner[2],part));push!(remaps,rows)
    end
    return _geo_geometry_receipt(carrier,parts,remaps,tags,caller)
end

function _geo_geometry_plan_receipt(authority,carrier,previous,plan,caller)
    carrier===authority.mesh && return nothing
    if hasproperty(plan,:geometry_parts) && plan.geometry_parts!==nothing
        return _geo_geometry_receipt(carrier,plan.geometry_parts,plan.geometry_part_rows,
            _geo_geometry_carrier_tags(authority,carrier,caller),caller)
    end
    return _geo_geometry_rebind_receipt(authority,carrier,previous,plan,caller)
end

function _geo_execution_geometry_parts(m,context)
    receipt=context.mesh_geometry_parts
    if receipt!==nothing
        receipt.carrier===context.mesh || throw(ArgumentError(
            "execute_geo: geometry entity parts belong to a stale carrier"))
        return receipt.parts
    end
    authority=m.meshing.mesh_authority
    authority===nothing || context.mesh===nothing || context.mesh===authority.mesh ||
        throw(ArgumentError("execute_geo: independent geometry carrier lacks its entity parts"))
    return context.mesh_parts
end

# The global and local refinement traces name primary rows and each actual
# inserted midpoint. A part never guesses a support vertex from its position.
function _geo_geometry_refined_parts(previous,trace,carrier,caller)
    previous===nothing && return nothing,nothing
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[];remaps=Vector{Int32}[]
    for ((dim,tag,part),old_rows) in zip(previous.parts,previous.rows)
        part isa Mesh || throw(ArgumentError("$caller: simplex refinement has a non-simplex geometry part"))
        if dim==0
            nsegs(part)==ntris(part)==ntets(part)==0 || throw(ArgumentError(
                "$caller: Point geometry part contains non-Point cells"))
            # Point parts have no local incident cells. Their liveness comes
            # from the global carrier trace, including incident Curve cells.
            rows=Int32[trace.primary_map[old] for old in old_rows
                if trace.primary_map[old]>0]
            allunique(rows) || throw(ArgumentError("$caller: refined Point rows repeat"))
            push!(parts,(dim,tag,Mesh(carrier.coords[:,rows])));push!(remaps,rows)
            continue
        end
        local_trace=Refine._refine_uniform_traced(part)
        rows=zeros(Int32,nnodes(local_trace.mesh))
        compact_old=zeros(Int32,count(!iszero,local_trace.primary_map))
        for (old,local_row) in enumerate(local_trace.primary_map)
            local_row==0 && continue
            compact_old[local_row]=Int32(old)
            global_row=trace.primary_map[old_rows[old]]
            global_row>0 || throw(ArgumentError("$caller: refined geometry part loses a primary row"))
            rows[local_row]=global_row
        end
        for ((a,b),local_row) in local_trace.midpoint_ids
            ca=trace.primary_map[old_rows[compact_old[a]]]
            cb=trace.primary_map[old_rows[compact_old[b]]]
            global_row=get(trace.midpoint_ids,minmax(ca,cb),nothing)
            global_row===nothing && throw(ArgumentError("$caller: refined geometry part loses a midpoint event"))
            rows[local_row]=global_row
        end
        all(>(0),rows) && allunique(rows) || throw(ArgumentError(
            "$caller: refined geometry part rows do not form an injection"))
        refined=_geo_geometry_part_coordinates(local_trace.mesh,carrier.coords[:,rows])
        push!(parts,(dim,tag,refined));push!(remaps,rows)
    end
    return parts,remaps
end

export execute_geo, GeoExecution, geo_entity_mesh

"""
    GeoExecution

Result of [`execute_geo`](@ref): the populated native `model`, an optional
validated `mesh`, and the sizing/field `params` recovered by
[`Tessella.IO.read_geo_params`](@ref).
"""
struct GeoExecution
    model::GeoModel
    mesh::Union{Nothing,Mesh,MixedMesh}
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
    # The `(dim, tag, mesh)` entity parts `mesh` was merged from — each part
    # carries the entity's own mesh in its own node space, so entity-level
    # consumers (`model_to_mixed` classified projection) recover the pure
    # per-entity view instead of the all-dimension product. Empty when no
    # mesh was generated through `_geo_mesh_model`.
    mesh_parts::Vector{Tuple{Int,Int,Union{Mesh,MixedMesh}}}
    # Actual native products retain their own UID/owner/UV metadata.
    actual_mesh_parts::Vector{Tuple{Int,Int,Union{Mesh,MixedMesh}}}
    # Factory maxima can exceed every live tag after renumbering or deletion.
    mesh_node_max::UInt64
    mesh_element_max::UInt64
    native_context_lc::Float64
    native_context_forced_bbox::Bool
    native_element_lookup_receipt::Any
end

GeoExecution(model,mesh,params,transfinite_tri,values,lists,strings,warnings,
             msg_error_count,exit_code,mesh_parts)=GeoExecution(model,mesh,params,
    transfinite_tri,values,lists,strings,warnings,msg_error_count,exit_code,
    mesh_parts,mesh_parts,UInt64(0),UInt64(0),1.0,false,nothing)

GeoExecution(model,mesh,params,transfinite_tri,values,lists,strings,warnings,
    msg_error_count,exit_code,mesh_parts,node_max,element_max,lc,forced_bbox)=
    GeoExecution(model,mesh,params,transfinite_tri,values,lists,strings,warnings,
        msg_error_count,exit_code,mesh_parts,mesh_parts,node_max,element_max,lc,forced_bbox,nothing)

GeoExecution(model,mesh,params,transfinite_tri,values,lists,strings,warnings,
    msg_error_count,exit_code,mesh_parts,actual_mesh_parts,node_max,element_max,lc,forced_bbox)=
    GeoExecution(model,mesh,params,transfinite_tri,values,lists,strings,warnings,
        msg_error_count,exit_code,mesh_parts,actual_mesh_parts,node_max,element_max,lc,forced_bbox,nothing)

"""
    geo_entity_mesh(execution, dim, tag; geometry=false) -> Union{Mesh,MixedMesh}

Return the `(dim, tag)` entity's actual native mesh from
`execution.actual_mesh_parts`, retaining its public identities, owners and
parameters. Set `geometry=true` to obtain the entity's emitted geometry view
from `execution.mesh_parts`; those parts partition `execution.mesh` and can
have separate coordinate arithmetic without claiming actual native metadata.
Entity-level consumers select the view required by their projection or
certificate instead of mixing the two row spaces. Throws
`ArgumentError` when no matching part exists (for example when the file
produced no `Mesh n`/`mesh_dim` product).
"""
function geo_entity_mesh(execution::GeoExecution,dim::Integer,tag::Integer;geometry::Bool=false)
    parts=geometry ? execution.mesh_parts : execution.actual_mesh_parts
    for (d,t,mesh) in parts
        (d==dim && t==tag) && return mesh
    end
    throw(ArgumentError(
        "geo_entity_mesh: no ($dim, $tag) part in the execution mesh"))
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
    open(path,"r") do stream
        _geo_exec_statements(stream)
    end
end

function _geo_exec_statements(stream::Base.IO)
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
    for raw in eachline(stream)
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
Use `mesh_dim=1`, `2` or `3` to mesh curves, surfaces or volumes — multiple
entities merge on bitwise coordinates like a `Mesh n` statement; `mesh_dim=0`
only builds the model. Geometry statements outside the bounded
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
volume relations as a Tessella extension. Affine curve relations accept
straight and curved curves — the declaration check is the endpoint
correspondence like upstream's `GEdge::setMeshMaster`, and the slave copies
the master's stored parameters at mesh time; orientation-only declarations
are check-free like upstream. Surfaces must be planar boundaries of one
explicit volume when meshed.
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
function execute_geo(path::AbstractString; mesh_dim::Integer=0,
                     _initial_mesh_algorithms=(algorithm_2d=6,algorithm_3d=1),
                     _initial_context_lc::Float64=1.0,
                     _initial_context_forced_bbox::Bool=false,
                     _initial_abort_on_error=nothing,
                     _initial_mesh_error_callback=nothing,
                     _initial_mesh_begin_callback=nothing,
                     _initial_abort_option_callback=nothing,
                     _initial_numeric_options=nothing,
                     _forced_numeric_options=nothing,
                     _initial_numeric_option_callback=nothing,
                     _initial_recombine_partial_callback=nothing,
                     _initial_recombine_diagnostics=nothing,
                     _initial_model_name=nothing,
                     _initial_model=nothing,
                     _initial_mesh_context=nothing,
                     _initial_fields=nothing,
                     _initial_onelab_numbers::Union{Nothing,Dict{String,Float64}}=nothing,
                     _initial_onelab_strings::Union{Nothing,Dict{String,String}}=nothing)
    isfile(path) || throw(ArgumentError("execute_geo: missing file $path"))
    mesh_dim isa Bool && throw(ArgumentError("execute_geo: mesh_dim must not be Bool"))
    dim=try
        Int(mesh_dim)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("execute_geo: mesh_dim exceeds the platform Int range"))
    end
    dim in 0:3 || throw(ArgumentError("execute_geo: mesh_dim must be in 0:3"))
    params=read_geo_params(path)
    # A model-file merge supplies an operation-owned detached model and context.
    # Ordinary standalone/open execution still starts with an empty model.
    model=_initial_model===nothing ? GeoModel() : _initial_model
    model isa GeoModel || throw(ArgumentError("execute_geo: initial model must be a GeoModel"))
    _initial_model_name===nothing || (model.name=String(_initial_model_name))
    context=_GeoNumericContext()
    _initial_mesh_context===nothing || _geo_seed_mesh_context!(context,_initial_mesh_context)
    # Session-bound runs (API `open`/`merge`/`open_geo!`) pass the session's
    # detached ONELAB scalar mirrors so `SetNumber`/`GetNumber`/`DefineNumber`
    # and their string forms exchange with the session parameter space.
    _initial_onelab_numbers===nothing ||
        (context.onelab_numbers=_initial_onelab_numbers)
    _initial_onelab_strings===nothing ||
        (context.onelab_strings=_initial_onelab_strings)
    context.native_context_lc=_initial_context_lc
    context.native_context_forced_bbox=_initial_context_forced_bbox
    context.native_mesh_error_callback=_initial_mesh_error_callback
    context.native_mesh_begin_callback=_initial_mesh_begin_callback
    context.native_abort_option_callback=_initial_abort_option_callback
    context.native_numeric_option_callback=_initial_numeric_option_callback
    context.native_recombine_partial_callback=_initial_recombine_partial_callback
    context.native_recombine_diagnostics=_initial_recombine_diagnostics
    _initial_abort_on_error===nothing ||
        (context.option_numbers[("General",0,"AbortOnError")]=Float64(_initial_abort_on_error))
    context.option_numbers[("Mesh",0,"Algorithm")]=
        Float64(_initial_mesh_algorithms.algorithm_2d)
    context.option_numbers[("Mesh",0,"Algorithm3D")]=
        Float64(_initial_mesh_algorithms.algorithm_3d)
    for (name,field,default) in (("RecombineOptimizeTopology",:recombine_topology_passes,5.0),
            ("RecombineNodeRepositioning",:recombine_node_repositioning,1.0),
            ("RecombineMinimumQuality",:recombine_minimum_quality,0.01),
            ("Smoothing",:native_smoothing,1.0),
            ("SecondOrderLinear",:second_order_linear,0.0))
        context.option_numbers[("Mesh",0,name)]=hasproperty(_initial_mesh_algorithms,field) ?
            Float64(getproperty(_initial_mesh_algorithms,field)) : default
    end
    if _initial_numeric_options!==nothing
        for (name,value) in _initial_numeric_options
            family,member=split(String(name),'.';limit=2)
            context.option_numbers[(family,0,member)]=Float64(value)
        end
    end
    # CLI options retain priority over subsequent script writes, including
    # aliases. Seed through the existing option setter to acquire all aliases;
    # only explicitly forced values are restored after each option assignment.
    if _forced_numeric_options!==nothing && !isempty(_forced_numeric_options)
        forced=Dict{Tuple{String,Int,String},Float64}()
        for (name,value) in _forced_numeric_options
            family,member=split(String(name),'.';limit=2)
            before=copy(context.option_numbers)
            _geo_store_option_number!(context,String(family),0,String(member),
                Float64(value),"execute_geo: command-line option")
            for (key,stored) in context.option_numbers
                if key==(family,0,member) || get(before,key,nothing)!=stored
                    forced[key]=stored
                end
            end
        end
        callback=context.native_numeric_option_callback
        context.native_numeric_option_callback=(family,member,value)->begin
            merge!(context.option_numbers,forced)
            callback===nothing || callback(family,member,
                get(forced,(family,0,member),value))
            nothing
        end
    end
    context.file_name=String(path)
    # Field declarations/options/delete replay into `context.fields` — the
    # statement-time live map, shared with `params.fields` so the returned
    # `GeoExecution` exposes the executed end state. `Field[i].member = v`
    # writes mutate it so a mid-file `Mesh n` sees the final option values.
    context.fields=empty!(params.fields)
    if _initial_fields!==nothing
        merge!(context.fields,deepcopy(_initial_fields.fields))
        context.background_field=_initial_fields.background_field
        context.boundary_layer_fields=copy(_initial_fields.boundary_layer_fields)
    end
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
    if _initial_model!==nothing
        _geo_allocator_resync_model!(allocator_state,model)
        allocator_state.physical_group_max=model.physical_tag_max
        allocator_state.live_fields=Set(Int(tag) for tag in keys(context.fields))
        for ((dimension,entity),physicals) in model.entity_physicals,physical in physicals
            push!(get!(() -> Int[],context.raw_physicals,(dimension,physical)),entity)
        end
    end
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
    if dim in 1:3
        _geo_sync_meshing_options!(model,context)
    end
    if dim==1
        mesh,context.mesh_node_owner=_geo_mesh_model(model,1,context)
        context.mesh=mesh
    elseif dim==2
        (isempty(model.surfaces) && !any(key->key[1]==2,
            Iterators.flatten((keys(model.discrete),keys(model.meshing.attached))))) &&
            throw(ArgumentError("execute_geo: Mesh 2 requested but no surfaces exist"))
        # Upstream `Mesh n` meshes every entity of dimension <= n regardless
        # of entity count: vertex parts (including orphan `Point` entities),
        # graded curve parts, and the surface parts merge on bitwise
        # coordinates with per-node ownership for the mesh operations.
        mesh,context.mesh_node_owner=_geo_mesh_model(model,2,context)
        context.mesh=mesh
    elseif dim==3
        (isempty(model.volumes) && !any(key->key[1]==3,
            Iterators.flatten((keys(model.discrete),keys(model.meshing.attached))))) &&
            throw(ArgumentError("execute_geo: Mesh 3 requested but no volumes exist"))
        mesh,context.mesh_node_owner=_geo_mesh_model(model,3,context)
        context.mesh=mesh
    end
    _geo_sync_mesh_algorithms!(model,context)
    _geo_sync_native_surface_controls!(model,context)
    if _initial_fields!==nothing
        # The execution map is shared with params.fields; selection commands
        # live on the context and need their literal final values transported.
        values=ntuple(length(fieldnames(GeoParams))) do index
            name=fieldnames(GeoParams)[index]
            name===:background_field ? context.background_field :
                name===:boundary_layer_fields ? copy(context.boundary_layer_fields) :
                getfield(params,index)
        end
        params=GeoParams(values...)
    end
    return GeoExecution(model,mesh,params,transfinite_tri,context.values,
                        context.lists,context.strings,
                        copy(context.exec_warnings),context.msg_error_count,
                        context.stop===:exit ? context.exit_code : nothing,
                        _geo_execution_geometry_parts(model,context),context.mesh_parts,
                        context.mesh_node_max,context.mesh_element_max,
                        context.native_context_lc,context.native_context_forced_bbox,
                        context.native_element_lookup_receipt)
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
            abort_ate_end=false
            # OCC-routed statements bind and publish whole shape graphs —
            # `Rectangle`'s boundary entities, loft faces, transform rebuilds
            # — that the text observer cannot count. A kernel signature
            # change under the OpenCASCADE factory marks the statement for
            # the model-driven resync below.
            occ_sig=allocator_state!==nothing &&
                allocator_state.factory===:opencascade ?
                (m.occ_kernel.next_shape,length(m.occ_kernel.bindings),
                 copy(m.occ_kernel.maxima)) : nothing
            assigned=try
                _exec_line!(m,line,context,allocator_state)
            catch err
                err isa InterruptException && rethrow()
                err isa _GeoSyntaxAbort || rethrow()
                # A `header` abort surfaces at the *next* token upstream:
                # the nested production already consumed this statement's
                # `;` (the `Split` Transform inside an `x =`/`x() =` list),
                # so the error lands on the following statement's head —
                # `$end` at the stream tail.
                token=err.header ?
                    (i+1<=lastindex(statements) ?
                        _geo_stmt_head_token(statements[i+1]) : "\$end") :
                    err.token
                _geo_yyerror!(context,"syntax error ($token)")
                syntax_aborted=true;abort_ate_end=err.header;nothing
            end
            assigned===nothing || (transfinite_tri=assigned)
            # `Delete`-family statements already drove the allocator update
            # inside `_geo_exec_delete!`; the observer skips them.
            syntax_aborted || _geo_allocator_observe_statement!(
                allocator_state,line,context,"execute_geo")
            occ_sig!==nothing &&
                (m.occ_kernel.next_shape,length(m.occ_kernel.bindings),
                 m.occ_kernel.maxima)!=occ_sig &&
                _geo_occ_allocator_resync!(allocator_state,m)
            # `occ_sig` snapshots the kernel before the statement; a pure
            # lazy-bind (e.g. `Closest` operand binding) also changes it —
            # the resync rebuild is idempotent on those statements.
            occursin(r"\bBoolean[A-Za-z]*\b",
                     line) &&
                _geo_allocator_resync_model!(allocator_state,m)
            # Statements that create entities outside the tracked allocator
            # path resync so `newp`-family reads skip their tags upstream-
            # style: `.msh` `Merge`/`Include` materialize discrete entities,
            # the discrete-model statements derive boundary entities, and
            # `Mesh n` runs queued homology requests whose chains share the
            # elementary namespace.
            match(Regex("^(?:Merge|Include|MergeWithBoundingBox|" *
                "CreateTopology|ClassifySurfaces|CreateGeometry|" *
                "Homology|Cohomology|Betti|Mesh)\\b"),line)!==nothing &&
                _geo_allocator_resync_model!(allocator_state,m)
            # `Split Curve{..} Point{..}` reallocates through `NEWCURVE()`
            # (`NEWREG` under `OldNewReg`) and rewires loops/physicals — the
            # observer invalidates on it, so resync to keep `newX` reads
            # live. The unanchored match also covers `s() = Split Curve…`.
            occursin(r"\bSplit\s+Curve\b",line) &&
                _geo_allocator_resync_model!(allocator_state,m)
            # `error tEND` recovery: a `;`-terminated statement's own
            # terminator ends the discard; a `}`-terminated one (Delete,
            # Boolean*, Extrude, transforms) has none, so recovery eats
            # through the next `;` like upstream. A `header` abort already
            # consumed the `;` (its production owned the terminator), so
            # the discard starts with the next statement.
            i=syntax_aborted &&
                (abort_ate_end || !endswith(rstrip(line),';')) ?
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

# `v() = ShapeFromFile("x.step")` — `OCC_Internals::shapeFromFile`: the
# same native STEP import `Merge` takes, returning the highest-dimensional
# imported tags so `Volume{v()}` can select them. Gated on the OCC factory
# like the boolean terms; the bare-statement form keeps its diagnostic.
function _geo_exec_shape_from_file(m::GeoModel,raw::AbstractString,
                                   context::_GeoNumericContext,
                                   allocator_state=nothing)
    source=String(strip(raw))
    mm=match(r"^ShapeFromFile\s*\(\s*(.*?)\s*\)\s*$",source)
    mm===nothing && return nothing
    caller="execute_geo: ShapeFromFile"
    if allocator_state===nothing || allocator_state.factory!==:opencascade
        _geo_yyerror!(context,
            "ShapeFromFile only available with OpenCASCADE geometry kernel")
        return Float64[]
    end
    file=_geo_fix_relative_path(context.file_name,
        _geo_eval_string(mm.captures[1],context,caller))
    out=try
        Model._model_occ_import_shapes(m,file;caller="OCC.importShapes")
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _geo_occ_fail!(context,err,"Could not import shapes")
        return Float64[]
    end
    _geo_occ_finish!(m,context)
    return Float64[tag for (_,tag) in out]
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
    context.native_context_lc=_geo_temp_bbox_lc(context)
    return nothing
end

# `CTX::instance()->lc` during `.geo` execution — the raw diagonal of the
# `Point`-statement cloud (`lc == 0` maps to 1 upstream).
function _geo_temp_bbox_lc(context::_GeoNumericContext)
    box=context.temp_bbox
    box===nothing && return 1.0
    range=box[2].-box[1]
    lc=sqrt(range[1]*range[1]+range[2]*range[2]+range[3]*range[3])
    return lc==0.0 ? 1.0 : lc
end

_geo_exec_lc(context::_GeoNumericContext)=context.native_context_lc

function _geo_sync_native_context!(m::GeoModel,context::_GeoNumericContext)
    context.native_context_forced_bbox && return nothing
    tolerance=something(_geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)
    context.native_context_lc=_model_mesh_lc(m,tolerance,"execute_geo: synchronize")
    return nothing
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
    if name=="Split"
        # `tSplit tCurve '{' FExpr '}' tPoint '{' RLD '}' tEND` — reachable
        # through `FExpr_Multi`, but the Transform production's own `tEND`
        # is the statement `;`. The enclosing `x =`/`x() =` production is
        # then left without a terminator: upstream reduces the transform
        # (the split side effect commits), reports the *next* statement's
        # token as a syntax error, and `error tEND` recovery discards that
        # statement. Without the `;` the Transform cannot reduce at all —
        # `x() = Split Curve{1} Point{2}, 3` errors at the `,` and the
        # split does not run — so the arguments are validated first.
        tail=String(strip(s[nextind(s,firstindex(s),5):end]))
        (curve_tag,point_tags,rest)=_geo_split_curve_args(
            tail,context,caller)
        isempty(rest) || _geo_selector_abort(rest)
        _geo_split_curve_run!(m,curve_tag,point_tags,context,caller,
                              allocator_state)
        throw(_GeoSyntaxAbort("Split","",true))
    end
    if name=="Intersect"
        # `tIntersect tCurve '{' RLD '}' tSurface '{' FExpr '}'` — a real
        # built-in kernel curve operation upstream; a hard blocker here.
        throw(ArgumentError(
            "$caller: Intersect requires built-in kernel curve splitting, " *
            "which Tessella does not implement"))
    end
    rest=String(strip(s[nextind(s,firstindex(s),
                              ncodeunits(name)):end]))
    ruled=false
    if name=="Ruled" && match(r"^ThruSections\b",rest)!==nothing
        # `tRuled tThruSections ListOfDouble` — the ruled-loft keyword pair.
        ruled=true
        name="ThruSections"
        rest=String(strip(rest[nextind(rest,firstindex(rest),12):end]))
    end
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
        if allocator_state!==nothing && allocator_state.factory===:opencascade
            entities=_geo_occ_apply_transform!(
                m,entities,t,context,caller,_GEO_OCC_TRANSFORM_WHAT[name])
        else
            transform_entities!(m,t,entities;caller=caller,
                                auto_coherence=_geo_auto_coherence(context))
            context.geo_changed=true
        end
        return Float64[Float64(tag) for (_,tag) in entities]
    end
    if !startswith(rest,"{")
        # `Affine`/`Closest`/`Fillet`/`Chamfer`/`ThruSections` are keywords
        # upstream — `Affine` without its `{...}` parameters is a syntax
        # error. Any other name may still be a `tSTRING` action or a
        # selector handled elsewhere.
        ruled && _geo_selector_abort(rest)
        name in ("Affine","Closest","Fillet","Chamfer","ThruSections") ||
            return nothing
        _geo_selector_abort(rest)
    end
    (params,rest)=_geo_balanced_group(rest,caller)
    if name=="Fillet" || name=="Chamfer" || name=="ThruSections"
        # `tFillet RLD RLD RLD` / `tChamfer RLD RLD RLD RLD` /
        # `tThruSections RLD` — every operand is a `ListOfDouble` group.
        groups=[params]
        ngroups=name=="Fillet" ? 3 : name=="Chamfer" ? 4 : 1
        while length(groups)<ngroups
            r2=String(strip(rest))
            startswith(r2,"{") || _geo_selector_abort(r2)
            (params,rest)=_geo_balanced_group(r2,caller)
            push!(groups,params)
        end
        isempty(rest) || _geo_selector_abort(rest)
        if allocator_state===nothing || allocator_state.factory!==:opencascade
            for g in groups
                _geo_numeric_list_values(
                    "{"*g*"}",context,"$caller parameters";depth=1)
            end
            _geo_yyerror!(context,
                "$name only available with OpenCASCADE geometry kernel")
            return Float64[]
        end
        kernel=_geo_occ_prepare_kernel!(m)
        if name=="ThruSections"
            wires=_geo_exec_entity_rhs_tags(
                "{"*groups[1]*"}",context,"$caller wires";abs_refs=true)
            out=_geo_occ_thru_sections_run(m,wires,context,caller,ruled)
            return Float64[Float64(tag) for (_,tag) in out]
        end
        out=_geo_occ_corner_run(m,name,groups,context,caller)
        out===nothing && return Float64[]
        return Float64[Float64(tag) for (_,tag) in out]
    end
    if name=="Affine" || name=="Closest"
        r2=String(strip(rest))
        startswith(r2,"{") || _geo_selector_abort(r2)
        (shape_list,rest)=_geo_balanced_group(r2,caller)
        isempty(rest) || _geo_selector_abort(rest)
        vals=_geo_numeric_list_values(
            "{"*params*"}",context,"$caller parameters";depth=1)
        entities=_geo_shape_list_entities!(m,shape_list,context,caller;
                                           allocator_state=allocator_state)
        if allocator_state===nothing || allocator_state.factory!==:opencascade
            # OCC-only transforms report a recoverable error upstream and
            # leave the shapes untouched: `Affine` still yields its input
            # list, while `Closest` contributes nothing.
            if name=="Affine"
                _geo_yyerror!(context,
                    "Affine transform only available with OpenCASCADE " *
                    "geometry kernel")
                return Float64[Float64(tag) for (_,tag) in entities]
            end
            _geo_yyerror!(context,
                "Closest entity only available with OpenCASCADE " *
                "geometry kernel")
            return Float64[]
        end
        kernel=_geo_occ_prepare_kernel!(m)
        if name=="Affine"
            try
                _geo_occ_bind_entities!(m,entities,kernel,caller)
                Model._model_occ_affine_transform!(m,entities,vals)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                _geo_occ_fail!(context,err,"Could not apply affine transform")
                return Float64[Float64(tag) for (_,tag) in entities]
            end
            _geo_occ_finish!(m,context)
            return Float64[Float64(tag) for (_,tag) in entities]
        end
        length(vals)==3 || throw(ArgumentError(
            "$caller: Closest requires `{x, y, z}` parameters"))
        out=try
            _geo_occ_bind_entities!(m,entities,kernel,caller)
            first(Model._model_occ_get_closest_entities(
                m,vals[1],vals[2],vals[3],entities,1))
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not get closest entity")
            Tuple{Int,Int}[]
        end
        return Float64[Float64(tag) for (_,tag) in out]
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
    out=_geo_shape_action!(m,name,entities,context,caller;
                           allocator_state=allocator_state)
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
    transform=_geo_exec_transform_term(m,s,context,allocator_state)
    transform!==nothing && return transform
    kw=match(r"^(Physical|Parent|BoundingBox|Mass|CenterOfMass|MatrixOfInertia|Normal|Parametric|Color)\b(.*)$",s)
    if kw!==nothing
        _geo_sync_physical_view_if_changed!(m,context)
        return _geo_selector_keyword_term(
            m,String(kw.captures[1]),String(strip(kw.captures[2])),context,caller,
            allocator_state)
    end
    head=_geo_selector_entity_head(s,context,caller)
    head===nothing && return nothing
    # The hook also sees ordinary numeric terms and option reads. Only an
    # admitted model selector triggers the native synchronization event.
    _geo_sync_physical_view_if_changed!(m,context)
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
    shapefile=_geo_exec_shape_from_file(m,raw,context,allocator_state)
    shapefile!==nothing && return shapefile
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
                           geometry=nothing,native_geo::Bool=true)
    itag=trunc(Int,tag)
    m.curves[itag]=isempty(found) ? (0,0) : (first(found),last(found))
    control_points && (m.curve_control_points[itag]=copy(found))
    kind!==nothing && (m.curve_types[itag]=kind)
    geometry!==nothing && (m.curve_geometry[itag]=geometry)
    native_geo && Model._native_geo_kernel_admit_curve!(m,itag)
    return nothing
end

function _geo_occ_line_geometry(m,found)
    length(found)==2 || return nothing
    origin=m.points[first(found)];endpoint=m.points[last(found)]
    delta=ntuple(axis->endpoint[axis]-origin[axis],3)
    edge_length=sqrt(delta[1]*delta[1]+delta[2]*delta[2]+delta[3]*delta[3])
    isfinite(edge_length) && edge_length>0 || return nothing
    direction=ntuple(axis->delta[axis]/edge_length,3)
    return (occ=:line,origin=origin,dir=direction,t0=0.0,t1=edge_length)
end

# BRepLib_FindSurface samples both endpoints of each straight edge, weighted
# by that edge's length. An axis-aligned wire has an exact covariance null
# axis; preserve that constructor instead of fitting a built-in mean Plane.
function _geo_occ_axis_plane_geometry(m,loops)
    curves=get(m.loops,first(loops),Int[])
    isempty(curves) && return nothing
    samples=NTuple{3,Float64}[];weights=Float64[]
    for signed in curves
        curve=abs(signed);geometry=get(m.curve_geometry,curve,nothing)
        geometry!==nothing && hasproperty(geometry,:occ) && geometry.occ===:line ||
            return nothing
        first_point=Model._occ_line_point(geometry,geometry.t0)
        last_point=Model._occ_line_point(geometry,geometry.t1)
        delta=ntuple(axis->last_point[axis]-first_point[axis],3)
        weight=sqrt(delta[1]*delta[1]+delta[2]*delta[2]+delta[3]*delta[3])
        isfinite(weight) && weight>0 || return nothing
        push!(samples,first_point);push!(samples,last_point)
        push!(weights,weight);push!(weights,weight)
    end
    axis=findfirst(axis->all(point->point[axis]==samples[1][axis],samples),1:3)
    axis===nothing && return nothing
    weighted=zeros(3);total=0.0
    for (point,weight) in zip(samples,weights)
        for component in 1:3
            weighted[component]+=weight*point[component]
        end
        total+=weight
    end
    center=(weighted[1]/total,weighted[2]/total,weighted[3]/total)
    first_axis=mod1(axis+1,3);second_axis=mod1(axis+2,3)
    area=0.0
    for signed in curves
        begin_tag,end_tag=m.curves[abs(signed)]
        signed<0 && ((begin_tag,end_tag)=(end_tag,begin_tag))
        begin_point=m.points[begin_tag];end_point=m.points[end_tag]
        area+=(begin_point[first_axis]-center[first_axis])*(end_point[second_axis]-center[second_axis])-
              (begin_point[second_axis]-center[second_axis])*(end_point[first_axis]-center[first_axis])
    end
    isfinite(area) && !iszero(area) || return nothing
    normal=ntuple(component->component==axis ? 1.0 : 0.0,3)
    area<0 && (normal=(-normal[1],-normal[2],-normal[3]))
    first_direction,second_direction=Model._occ_ax2_frame(normal)
    pcurves=Dict{Int,NamedTuple}()
    base_geometry=(occ=:plane,center=center,axis=normal,X=first_direction,Y=second_direction,
        reversed=false,pcurves=pcurves)
    parameters=[Model._occ_plane_uv(base_geometry,point) for point in samples]
    lower=Float64[minimum(point[component] for point in parameters) for component in 1:2]
    upper=Float64[maximum(point[component] for point in parameters) for component in 1:2]
    geometry=merge(base_geometry,(uvb=(lower,upper),))
    for loop in loops,signed in m.loops[loop]
        curve=abs(signed);edge,low,high=Model._brep_edge_curve(m,curve)
        pcurves[curve]=(fwd=Model._brep_synth_pcurve(geometry,edge,low,high,
            "Plane Surface"),rev=nothing)
    end
    return geometry
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
    if mesh isa Mesh
        (size(mesh.segs,2)+size(mesh.tris,2)+size(mesh.tets,2))>0 &&
            return true
    elseif mesh isa MixedMesh
        for block in mesh.blocks
            block isa ElementBlock && size(block.nodes,2)>0 &&
                return true
            block isa SpecialElementBlock &&
                length(block.tags)>0 && return true
        end
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
    _geo_restore_orientation_constraints!(m)
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
    _geo_sync_native_context!(m,context)
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
            vals=_geo_numeric_list_values(
                "{"*params*"}",context,"$caller $name";depth=1)
            s=String(strip(s))
            startswith(s,"{") || _geo_selector_abort(s)
            (inner_list,s)=_geo_balanced_group(s,nested)
            inner=_geo_shape_list_entities!(m,inner_list,context,caller;
                                            allocator_state=allocator_state)
            if allocator_state===nothing ||
                    allocator_state.factory!==:opencascade
                # OCC-only transforms report a recoverable error upstream and
                # leave the shapes untouched: `Affine` still yields its input
                # list, while `Closest` contributes nothing.
                if name=="Affine"
                    _geo_yyerror!(context,
                        "Affine transform only available with OpenCASCADE " *
                        "geometry kernel")
                    return (inner,s,true)
                end
                _geo_yyerror!(context,
                    "Closest entity only available with OpenCASCADE " *
                    "geometry kernel")
                return (NTuple{2,Int}[],s,true)
            end
            kernel=_geo_occ_prepare_kernel!(m)
            if name=="Affine"
                try
                    _geo_occ_bind_entities!(m,inner,kernel,nested)
                    Model._model_occ_affine_transform!(m,inner,vals)
                catch err
                    err isa InterruptException && rethrow()
                    err isa ArgumentError || rethrow()
                    _geo_occ_fail!(context,err,
                                   "Could not apply affine transform")
                    return (inner,s,true)
                end
                _geo_occ_finish!(m,context)
                return (inner,s,true)
            end
            length(vals)==3 || throw(ArgumentError(
                "$nested: Closest requires `{x, y, z}` parameters"))
            out=try
                _geo_occ_bind_entities!(m,inner,kernel,nested)
                first(Model._model_occ_get_closest_entities(
                    m,vals[1],vals[2],vals[3],inner,1))
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                _geo_occ_fail!(context,err,"Could not get closest entity")
                Tuple{Int,Int}[]
            end
            return (NTuple{2,Int}[e for e in out],s,true)
        end
        (t,s)=_geo_transform_params_rest(
            kind=name,source=s,context=context,caller=nested)
        s=String(strip(s))
        startswith(s,"{") || _geo_selector_abort(s)
        (inner_list,s)=_geo_balanced_group(s,nested)
        inner=_geo_shape_list_entities!(m,inner_list,context,caller;
                                        allocator_state=allocator_state)
        if allocator_state!==nothing && allocator_state.factory===:opencascade
            _geo_occ_apply_transform!(
                m,inner,t,context,nested,_GEO_OCC_TRANSFORM_WHAT[name])
        else
            transform_entities!(m,t,inner;caller=nested,
                                auto_coherence=_geo_auto_coherence(context))
        end
        # Gmsh returns the input shape list (`$$ = $MultipleShape`).
        return (inner,s,true)
    end
    if name=="Fillet" || name=="Chamfer" || name=="ThruSections" ||
            (name=="Ruled" && match(r"^ThruSections\b",s)!==nothing)
        # `tFillet`/`tChamfer`/`(tRuled) tThruSections` Transform productions
        # — the output entity list, not the input shapes, is yielded.
        ruled=false
        if name=="Ruled"
            ruled=true;name="ThruSections"
            s=String(strip(s[nextind(s,firstindex(s),12):end]))
        end
        nested="$caller: $name"
        ngroups=name=="Fillet" ? 3 : name=="Chamfer" ? 4 : 1
        groups=String[]
        for _ in 1:ngroups
            startswith(s,"{") || _geo_selector_abort(s)
            (params,s)=_geo_balanced_group(s,nested)
            push!(groups,params)
            s=String(strip(s))
        end
        s=_geo_require_list_semicolon(s,caller,name)
        if allocator_state===nothing || allocator_state.factory!==:opencascade
            for g in groups
                _geo_numeric_list_values(
                    "{"*g*"}",context,"$caller $name";depth=1)
            end
            _geo_yyerror!(context,
                "$name only available with OpenCASCADE geometry kernel")
            return (NTuple{2,Int}[],s,true)
        end
        out=if name=="ThruSections"
            wires=_geo_exec_entity_rhs_tags(
                "{"*groups[1]*"}",context,"$caller wires";abs_refs=true)
            _geo_occ_thru_sections_run(m,wires,context,nested,ruled)
        else
            _geo_occ_corner_run(m,name,groups,context,nested)
        end
        out===nothing && return (NTuple{2,Int}[],s,true)
        return (NTuple{2,Int}[e for e in out],s,true)
    end
    if name=="Split" || name=="Intersect"
        # `tSplit tCurve '{' FExpr '}' tPoint '{' RLD '}'` /
        # `tIntersect tCurve '{' RLD '}' tSurface '{' FExpr '}'` — any other
        # head is a syntax error (the keyword expects `Curve`). `Split`
        # yields the new curve records; `Intersect` is a built-in kernel
        # blocker.
        match(r"^Curve\b",s)===nothing && _geo_selector_abort(s)
        name=="Intersect" && throw(ArgumentError(
            "$caller: Intersect requires built-in kernel curve splitting, " *
            "which Tessella does not implement"))
        (new_tags,s)=_geo_exec_split_curve!(
            m,s,context,"$caller: $name",allocator_state)
        # The Transform production terminates on `tEND` — the member `;` is
        # part of the `Split` production, like every other shape-list entry.
        s=_geo_require_list_semicolon(s,caller,name)
        return ([(1,tag) for tag in new_tags],s,true)
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
        return (_geo_shape_action!(m,name,inner,context,caller;
                                   allocator_state=allocator_state),s,true)
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
                            caller::AbstractString;
                            allocator_state=nothing)
    if name=="Duplicata"
        if allocator_state!==nothing && allocator_state.factory===:opencascade
            # `OCC_Internals::copy` — `_multiBind` binds the inputs at their
            # tags, duplicates, and reports the new highest-dimension tags.
            kernel=_geo_occ_prepare_kernel!(m)
            out=try
                _geo_occ_bind_entities!(m,inner,kernel,"$caller: Duplicata")
                Model._model_occ_copy_entities!(m,inner)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                _geo_occ_fail!(context,err,"Could not copy shapes")
                nothing
            end
            out===nothing && return NTuple{2,Int}[]
            _geo_occ_finish!(m,context)
            return NTuple{2,Int}[(d,t) for (d,t) in out]
        end
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

# Automatic built-in geometry coherence follows the live option; explicit
# Coherence commands call the unconditional model merge.
_geo_auto_coherence(context::_GeoNumericContext) =
    trunc(Int,something(_geo_option_number(context,"Geometry",0,"AutoCoherence"),1.0))

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
                               allocator_state=nothing;
                               allow_shapes::Bool=true)
    mm=match(r"^([A-Za-z_][A-Za-z0-9_]*)",element)
    mm===nothing && throw(ArgumentError(
        "$caller: malformed shape list element near $(repr(element))"))
    name=mm.captures[1]
    # Gmsh.l intentionally permits these selected letter-case variants,
    # without making keywords generally case-insensitive.
    if startswith(name,"Quad")
        occursin(r"^Quad[tT]ri[aA]dd[vV]erts$",name) &&
            (name="QuadTriAddVerts")
        occursin(r"^Quad[tT]ri[nN]o[nN]ew[vV]erts$",name) &&
            (name="QuadTriNoNewVerts")
    end
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
    elseif name=="ScaleLastLayer"
        isempty(s) || throw(ArgumentError(
            "$caller: ScaleLastLayer takes no value"))
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
        recomb=match(r"^Recomb[lL]aterals\b",s)!==nothing
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
            "$caller: Hole extrusion parameters on boundary-layer " *
            "extrusion are not implemented"))
    elseif name in ("Extrude","Translate","Rotate","Dilate","Symmetry","Affine",
                    "Duplicata","Boundary","CombinedBoundary",
                    "OrientedBoundary","OrientedCombinedBoundary",
                    "PointsOf","Split","Intersect","Closest")
        # `ListOfShapes` contains entity selectors and Shape definitions,
        # never Transform or Extrude productions. Extrude expressions are
        # still valid inside an entity's numeric tag list (`Point{Extrude
        # {...}{Point{1};}};`), where the value-term hook evaluates them.
        _geo_syntax_abort(name)
    end
    # Gmsh's Extrude production is `ListOfShapes ExtrudeParameters`: after
    # the first parameter the parser cannot return to shape entries. Check
    # the head before evaluating its tag expressions or inline definitions,
    # since a rejected entry must not commit those side effects.
    allow_shapes || _geo_syntax_abort(name)
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
    allow_shapes=true
    isempty(strip(body)) && return (entities,params)
    for element in _geo_exec_topology_query_blocks(
            body,"Extrude shape list",caller)
        entries,params=_geo_extrude_element!(
            m,element,context,params,caller,allocator_state;
            allow_shapes=allow_shapes)
        if entries===nothing
            allow_shapes=false
        else
            append!(entities,entries)
        end
        length(entities)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
            "$caller: shape list expands beyond $_MAX_GEO_LIST_ITEMS entities"))
    end
    return (entities,params)
end

# ── `Extrude {shapes; params}` — `BOUNDARY_LAYER` ────────────────────────────
# Gmsh's no-displacement extrusion mode: each curve generatrix grows a
# boundary-layer surface strip along the smoothed source-edge normal field,
# with `ExtrudePoint` connector curves and a `DuplicateCurve` top copy wiring
# the strip's boundary like an ordinary curve extrusion. The
# `mesh_boundary_layer_2d` kernel emits the certified strip; helpers shared
# with the mesh pass live in `ModelExtrude.jl` (`_extrude_bl_*`).

# `ExtrudePoint` under `BOUNDARY_LAYER` — `DuplicateVertex` + connecting
# `CreateCurve`. The copy's position and the connector's spec record are set
# when the wall's kernel run writes them back (`_extrude_bl_refresh!`).
function _geo_bl_point_copy!(m::GeoModel,src_pt::Int,params,caller)
    dest=_fresh_point_copy!(m,src_pt,caller)
    curve=_geo_newreg_alloc!(m,1,caller)
    m.curves[curve]=(src_pt,dest)
    m.curve_control_points[curve]=Int[src_pt]
    m.meshing.extrude[(1,curve)]=params
    m.meshing.extrude_sources[(1,curve)]=(0,src_pt)
    return (curve,dest)
end

# `ExtrudeCurve` under `BOUNDARY_LAYER`: `DuplicateCurve` chapeau (the signed
# generatrix's reversed record duplicates end-to-begin), `ExtrudePoint`
# connectors off each signed endpoint, and the lateral surface wired
# `pc, CurveEnd, -chapeau, -CurveBeg`. A closed wall still extrudes both
# coincident connectors like upstream — the dest vertices merge in the
# operation's coherence pass.
function _geo_bl_curve_lateral!(m::GeoModel,c::Int,params,closed,kernel,
                                working,caller)
    src=abs(c)
    a,b=m.curves[src]
    sbeg,send=c>0 ? (a,b) : (b,a)
    chapeau=_duplicate_curve!(m,src,caller;reversed=c<0)
    m.meshing.extrude[(1,chapeau)]=params
    # `geo.Source = ic` — the signed generatrix, like the ordinary path.
    m.meshing.extrude_sources[(1,chapeau)]=(1,c)
    (cbeg,_)=_geo_bl_point_copy!(m,sbeg,params,caller)
    (cend,_)=_geo_bl_point_copy!(m,send,params,caller)
    surf=_geo_newreg_alloc!(m,2,caller)
    generatrices=Int[c,cend,-chapeau,-cbeg]
    lt=_geo_derived_loop_tag(m,surf)
    m.loops[lt]=generatrices
    m.surfaces[surf]=Int[lt]
    m.surface_types[surf]=:ruled
    m.meshing.extrude[(2,surf)]=params
    m.meshing.extrude_sources[(2,surf)]=(1,c)
    spec=(type=:boundary_layer,connectors=(cbeg,cend),chapeau=chapeau)
    m.meshing.extrude_specs[(2,surf)]=spec
    _extrude_bl_refresh!(m,spec,src,c,closed,kernel,caller;
                         _params_model=working)
    return (chapeau,surf)
end

# The merge that welds sibling dest vertices also collapses coincident
# connector/chapeau records; references inside `extrude_sources` values and
# the boundary-layer `extrude_specs` still name the dropped tags. Rewrite
# them through the coherence maps like `ReplaceAllDuplicates` implies
# upstream, so every stored reference names a live entity.
function _geo_bl_rebind!(m::GeoModel,maps)
    for (key,src) in collect(m.meshing.extrude_sources)
        sdim,stag=src
        if sdim==0
            keep=get(maps.points,stag,0)
            keep!=0 && (m.meshing.extrude_sources[key]=(0,keep))
        elseif sdim==1
            keep=get(maps.curves,abs(stag),stag)
            nk=sign(stag)*keep
            nk!=stag && (m.meshing.extrude_sources[key]=(1,nk))
        end
    end
    for (key,spec) in collect(m.meshing.extrude_specs)
        # Only the dim-2 lateral specs carry `connectors`/`chapeau` — the
        # dim-1 chapeau record shares `type=:boundary_layer` but stores
        # `tops`, which has no tags to rebind.
        (key[1]==2 && spec isa NamedTuple &&
         get(spec,:type,nothing)===:boundary_layer) || continue
        cb,cn=spec.connectors
        cb2=abs(get(maps.curves,abs(cb),cb))
        cn2=abs(get(maps.curves,abs(cn),cn))
        ch2=abs(get(maps.curves,abs(spec.chapeau),spec.chapeau))
        (cb2,cn2,ch2)==(cb,cn,spec.chapeau) && continue
        m.meshing.extrude_specs[key]=(type=:boundary_layer,
            connectors=(cb2,cn2),chapeau=ch2)
    end
    return nothing
end

function _geo_boundary_layer_extrude!(m::GeoModel,
                                      entities::AbstractVector{<:NTuple{2,Int}},
                                      params,context,caller)
    for (dim,tag) in entities
        dim==1 || throw(ArgumentError(
            "$caller: boundary-layer extrusion of " *
            "$(_entity_label(dim))[$tag] is not implemented"))
    end
    isempty(entities) && return Int[]
    isempty(params.layers) && throw(ArgumentError(
        "$caller: boundary-layer extrusion requires a Layers{...} parameter"))
    params.scale_last && throw(ArgumentError(
        "$caller: ScaleLastLayer on boundary-layer extrusion is not " *
        "implemented"))
    params.quad_to_tri!==:none && throw(ArgumentError(
        "$caller: QuadTri parameters on boundary-layer extrusion are not " *
        "implemented"))
    # Snapshot every wall's discretization on a stage copy so `Extrude` keeps
    # no mesh side effects; the mesh pass re-grades the live chains and
    # refreshes the sibling specs the same way. The kernel runs eagerly here
    # so malformed walls fail at the statement like upstream's `Extrude`.
    working=_extrude_curve_params_stage(m)
    walls=Tuple{Vector{NTuple{3,Float64}},Bool,NTuple{3,Float64}}[]
    chains=Dict{Int,Vector{NTuple{3,Float64}}}()
    for (_,c) in entities
        src=abs(c)
        haskey(m.curves,src) || throw(ArgumentError(
            "$caller: unknown Curve[$src]"))
        a,b=m.curves[src]
        (haskey(m.points,a) && haskey(m.points,b)) || throw(ArgumentError(
            "$caller: cannot extrude Curve[$src] with no begin or end point"))
        closed=a==b
        chain=_extrude_curve_nodes(working,src,caller)
        chains[c]=chain
        push!(walls,(chain,closed,
            _extrude_bl_plane_normal(chain,closed,c>0,caller)))
    end
    field=_extrude_bl_normal_field(walls)
    offsets=_extrude_level_us(params)[2:end]
    made=NTuple{3,Int}[]
    for (i,(_,c)) in enumerate(entities)
        chain=chains[c]
        closed=walls[i][2]
        kernel=_extrude_bl_kernel(chain,closed,walls[i][3],field,offsets,caller)
        (chapeau,surf)=_geo_bl_curve_lateral!(m,c,params,closed,kernel,
            working,caller)
        push!(made,(c,chapeau,surf))
    end
    maps=_coherence_merge!(m;enabled=!iszero(_geo_auto_coherence(context)),
                           remove_degenerated=_geo_auto_coherence(context)>1)
    _geo_bl_rebind!(m,maps)
    out=Int[]
    for (c,chapeau,surf) in made
        top=get(maps.curves,chapeau,chapeau)
        haskey(m.curves,abs(top)) && push!(out,top)
        haskey(m.surfaces,surf) || continue
        push!(out,surf)
        # `Geometry.ExtrudeReturnLateralEntities` — the connector
        # generatrices. Upstream compares `abs(record tag)` against the
        # *signed* input, so a negative generatrix keeps its own `-tag`.
        iszero(context.extrude_return_lateral) && continue
        for g in m.loops[only(m.surfaces[surf])]
            (abs(g)==c || abs(g)==abs(top)) || push!(out,g)
        end
    end
    return out
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
    if match(r"^Using\b",rest)!==nothing
        return _geo_exec_extrude_pipe!(m,rest,shapes,nvec,context,caller,
                                       allocator_state)
    end
    isempty(rest) || throw(ArgumentError(
        "$caller: unexpected text after the shape list"))
    shapes===nothing && throw(ArgumentError(
        "$caller: expected a `{shape list}` group"))
    (nvec==0 || revolve===nothing || delta==(0.0,0.0,0.0)) ||
        throw(ArgumentError("$caller: malformed displacement expression"))
    # Gmsh.y's `ListOfShapes` parses before the parser action runs — surface
    # its errors ahead of the boundary-layer/`TRANSLATE_ROTATE` gates like
    # upstream, whose OCC gate also sits between the two.
    entities,params=_geo_extrude_shape_list!(m,shapes,context,caller,
                                             allocator_state)
    if nvec==0
        # `BOUNDARY_LAYER` (`Extrude {shapes; params}`): upstream rejects it
        # in the parser action under OCC before `ExtrudeShapes` runs.
        allocator_state!==nothing &&
            allocator_state.factory===:opencascade && throw(ArgumentError(
                "$caller: Boundary layer extrusion not available with " *
                "OpenCASCADE geometry kernel"))
        tags=_geo_boundary_layer_extrude!(m,entities,params,context,caller)
        context.geo_changed=true
        return Float64.(tags)
    end
    if revolve!==nothing && revolve.kind===:twist &&
        allocator_state!==nothing &&
        allocator_state.factory===:opencascade
        throw(ArgumentError(
            "$caller: Twisting extrude not available with OpenCASCADE " *
            "geometry kernel"))
    end
    tags=revolve===nothing ?
        extrude_entities!(m,entities,delta;params=params,
            return_lateral=context.extrude_return_lateral,caller=caller,
            auto_coherence=_geo_auto_coherence(context)) :
        revolve.kind===:twist ?
        twist_entities!(m,entities,revolve.axis,revolve.origin,
            revolve.delta,revolve.angle;
            spline_points=_geo_extrude_spline_points(context),
            params=params,
            return_lateral=context.extrude_return_lateral,caller=caller,
            auto_coherence=_geo_auto_coherence(context)) :
        revolve_entities!(m,entities,revolve.axis,revolve.origin,
            revolve.angle;params=params,
            return_lateral=context.extrude_return_lateral,caller=caller,
            auto_coherence=_geo_auto_coherence(context))
    context.geo_changed=true
    return Float64.(tags)
end

# ── `Extrude {shapes} Using Wire {expr}` — OCC pipe extrusion ────────────────
# `.geo` OCC-factory statements keep entities as GEO records carrying OCC
# geometry (the same world `ThruSections`/`Pipe` bind through their own
# paths). The pipe itself runs in the OCC kernel, so the referenced profile
# and spine entities bind into `m.occ_kernel` at their existing tags first —
# vertices/edges/faces from `m.points`/`m.curves`/`m.surfaces`, wires from
# `m.loops` — then `Model._model_occ_add_pipe!` executes the parallel-
# transport loft and `_model_occ_synchronize!` publishes the generated shapes
# back as GEO entities. Under the built-in factory the suffix still rejects
# (the statement parses first, matching upstream's action-order gates).

# Bind a `.geo` point into the OCC kernel as a vertex at its own tag.
function _geo_occ_pipe_vertex!(m::GeoModel,ptag::Int,kernel,caller)
    haskey(kernel.bindings,(0,ptag)) && return kernel.bindings[(0,ptag)]
    haskey(m.points,ptag) || throw(ArgumentError(
        "$caller: unknown Point[$ptag]"))
    vid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,0,ptag,vid,
        Model._ModelOCCShape(0,
            (point=m.points[ptag],mesh_size=get(m.point_size,ptag,0.0)),
            UInt64[]))
    kernel.next_shape=vid+UInt64(1)
    return vid
end

# Bind a `.geo` curve into the OCC kernel as an edge. Curves that already
# carry an OCC geometry record (`Line`, `Circle` under the OCC factory)
# bind verbatim; `Spline`/`BSpline`/`Bezier` materialize through the OCC
# add functions (which interpolate poles the way `occ.addSpline` does);
# a bare two-point curve falls back to the `:line` chord record. Kinds
# with no OCC representation (`Ellipse`, `Nurbs`, built-in fillings)
# reject rather than silently degrading the profile.
function _geo_occ_pipe_edge!(m::GeoModel,tag::Int,kernel,caller)
    haskey(kernel.bindings,(1,tag)) && return kernel.bindings[(1,tag)]
    haskey(m.curves,tag) || throw(ArgumentError(
        "$caller: unknown Curve[$tag]"))
    kind=get(m.curve_types,tag,:line)
    g=get(m.curve_geometry,tag,nothing)
    if g===nothing || !hasproperty(g,:occ)
        if kind===:spline || kind===:bspline || kind===:bezier
            control=get(m.curve_control_points,tag,Int[])
            length(control)>=2 || throw(ArgumentError(
                "$caller: Curve[$tag] has too few control points to bind"))
            for ptag in control
                _geo_occ_pipe_vertex!(m,ptag,kernel,caller)
            end
            if kind===:spline
                Model._model_occ_add_spline!(m,copy(control);tag=tag)
            elseif kind===:bspline
                Model._model_occ_add_bspline!(m,copy(control);tag=tag)
            else
                Model._model_occ_add_curve!(m,copy(control),:bezier;tag=tag)
            end
            return kernel.bindings[(1,tag)]
        end
        kind===:line || throw(ArgumentError(
            "$caller: Curve[$tag] kind $(kind) has no OCC geometry " *
            "record for pipe extrusion"))
        g=_geo_occ_line_geometry(m,collect(m.curves[tag]))
        g===nothing && throw(ArgumentError(
            "$caller: Curve[$tag] endpoints do not yield a valid line"))
    end
    endpoints=m.curves[tag]
    vertex_ids=UInt64[]
    for ptag in endpoints
        push!(vertex_ids,_geo_occ_pipe_vertex!(m,ptag,kernel,caller))
    end
    eid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,1,tag,eid,
        Model._ModelOCCShape(1,(kind=kind,geometry=g),vertex_ids))
    kernel.next_shape=eid+UInt64(1)
    return eid
end

# Bind a `.geo` `Curve Loop` into the kernel as a wire. The stored loop is the
# source of truth for traversal order and signs (the `Curve Loop` sorter
# already ordered the chain); a negative curve tag walks the edge in reverse.
# Outer wires bind at the loop's own tag so `Using Wire {n}` references them;
# a hole use binds a fresh-tagged wire carrying the reversed traversal, the
# same shape convention `_model_occ_synchronize!` round-trips.
function _geo_occ_pipe_wire!(m::GeoModel,tag::Int,kernel,caller;
                             reversed::Bool=false)
    !reversed && haskey(kernel.bindings,(-1,tag)) &&
        return kernel.bindings[(-1,tag)]
    haskey(m.loops,tag) || throw(ArgumentError(
        "$caller: unknown Wire/Curve Loop[$tag]"))
    signed=Int.(m.loops[tag])
    seq=reversed ? Int[-s for s in Iterators.reverse(signed)] : signed
    uses=Tuple{UInt64,Bool}[]
    for s in seq
        eid=_geo_occ_pipe_edge!(m,abs(s),kernel,caller)
        push!(uses,(eid,s<0))
    end
    wtag=reversed ? Model._model_occ_next_tag(kernel,-1,-1,caller) : tag
    wid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,-1,wtag,wid,
        Model._ModelOCCShape(-1,(uses=uses,),[u[1] for u in uses]))
    kernel.next_shape=wid+UInt64(1)
    return wid
end

# A single `Curve{t}` profile becomes a fresh (open) kernel wire — upstream
# `OCCGeoInternals` wraps the lone edge the same way.
function _geo_occ_pipe_profile_wire!(m::GeoModel,tag::Int,kernel,caller)
    eid=_geo_occ_pipe_edge!(m,tag,kernel,caller)
    wtag=Model._model_occ_next_tag(kernel,-1,-1,caller)
    wid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,-1,wtag,wid,
        Model._ModelOCCShape(-1,(uses=((eid,false),),),[eid]))
    kernel.next_shape=wid+UInt64(1)
    return wtag
end

# Bind a `.geo` `Surface` (plane surface with boundary loops) into the kernel
# as a face. `Plane Surface` already records an OCC plane geometry under the
# OCC factory; otherwise an axis-aligned loop can still produce one through
# `_geo_occ_axis_plane_geometry`, and anything else rejects.
function _geo_occ_pipe_face!(m::GeoModel,tag::Int,kernel,caller)
    haskey(kernel.bindings,(2,tag)) && return kernel.bindings[(2,tag)]
    haskey(m.surfaces,tag) || throw(ArgumentError(
        "$caller: unknown Surface[$tag]"))
    loops=m.surfaces[tag]
    isempty(loops) && throw(ArgumentError(
        "$caller: Surface[$tag] has no boundary loops"))
    wire_ids=UInt64[]
    for (i,lt) in enumerate(loops)
        push!(wire_ids,_geo_occ_pipe_wire!(m,Int(lt),kernel,caller;
            reversed=i>1))
    end
    kind=get(m.surface_types,tag,:plane)
    g=get(m.surface_geometry,tag,nothing)
    if g===nothing && kind===:plane
        g=_geo_occ_axis_plane_geometry(m,loops)
    end
    (g===nothing || !hasproperty(g,:occ)) && throw(ArgumentError(
        "$caller: Surface[$tag] has no OCC geometry record for pipe " *
        "extrusion"))
    fid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,2,tag,fid,
        Model._ModelOCCShape(2,(kind=kind,geometry=g),wire_ids))
    kernel.next_shape=fid+UInt64(1)
    return fid
end

function _geo_exec_extrude_pipe!(m::GeoModel,rest::AbstractString,shapes,
        nvec::Int,context::_GeoNumericContext,caller::AbstractString,
        allocator_state)
    nvec>0 && throw(ArgumentError(
        "$caller: `Extrude .. Using Wire` does not take a motion or " *
        "displacement group"))
    urest=String(strip(rest[nextind(rest,firstindex(rest),5):end]))
    match(r"^Wire\b",urest)===nothing && throw(ArgumentError(
        "$caller: unsupported `Extrude .. Using` suffix (only " *
        "`Using Wire {expr}` is supported)"))
    wrest=String(strip(urest[nextind(urest,firstindex(urest),4):end]))
    startswith(wrest,"{") || throw(ArgumentError(
        "$caller: expected `Using Wire {expr}`"))
    (wire_group,tail)=_geo_balanced_group(wrest,caller)
    isempty(strip(tail)) || throw(ArgumentError(
        "$caller: unexpected text after `Using Wire {..}`"))
    shapes===nothing && throw(ArgumentError(
        "$caller: expected a `{shape list}` group"))
    # `ListOfShapes` resolves before the factory gate, matching the upstream
    # parser action order (shape-list syntax and entity lookups surface first).
    entities,params=_geo_extrude_shape_list!(m,shapes,context,caller,
                                             allocator_state)
    (allocator_state!==nothing &&
        allocator_state.factory===:opencascade) || throw(ArgumentError(
        "$caller: pipe extrusion `Extrude {..} Using Wire {..}` is " *
        "OpenCASCADE-only"))
    (isempty(params.layers) && isempty(params.heights) && !params.scale_last &&
        !params.recombine && params.quad_to_tri===:none &&
        !params.recomb_laterals) || throw(ArgumentError(
        "$caller: `Extrude .. Using Wire` with mesh parameters " *
        "(Layers/Recombine/QuadTri*/ScaleLast) is not implemented"))
    wire_values=_geo_exec_numeric_values(wire_group,context,
                                         "$caller Using Wire")
    length(wire_values)==1 || throw(ArgumentError(
        "$caller: `Using Wire {expr}` expects one wire tag, got " *
        "$(length(wire_values))"))
    wire_tag=_geo_signed_gmsh_int_value(wire_values[1],"$caller Using Wire")
    wire_tag>0 || throw(ArgumentError(
        "$caller: `Using Wire` tag must be positive (got $wire_tag)"))
    kernel=m.occ_kernel
    # Keep kernel tag allocation ahead of every `.geo` namespace it may emit:
    # entity maxima per dimension, loop tags for wires, and surface-loop tags
    # for the shells a solid pipe creates.
    for d in 0:3
        kernel.maxima[d+3]=max(kernel.maxima[d+3],m.next_tag[d+1])
    end
    isempty(m.loops) ||
        (kernel.maxima[2]=max(kernel.maxima[2],maximum(keys(m.loops))))
    isempty(m.surface_loops) ||
        (kernel.maxima[1]=max(kernel.maxima[1],maximum(keys(m.surface_loops))))
    _geo_occ_pipe_wire!(m,wire_tag,kernel,caller)
    out=NTuple{2,Int}[]
    for (dim,tag) in entities
        if dim==2
            _geo_occ_pipe_face!(m,tag,kernel,caller)
            append!(out,Model._model_occ_add_pipe!(m,-1,(2,tag),wire_tag))
        elseif dim==1
            wtag=_geo_occ_pipe_profile_wire!(m,tag,kernel,caller)
            append!(out,Model._model_occ_add_pipe!(m,-1,wtag,wire_tag))
        else
            throw(ArgumentError(
                "$caller: cannot pipe-extrude $(_entity_label(dim))[$tag]"))
        end
    end
    Model._model_occ_synchronize!(m)
    context.geo_changed=true
    return Float64[Float64(t) for (_,t) in out]
end

# ── OCC kernel routing for `.geo` OCC-only statements ────────────────────────
# Gmsh executes `Rectangle`, `Disk`, `Wedge`, `Wire`, `ThickSolid`,
# `ThruSections`/`Ruled ThruSections`, `BSpline|Bezier Surface`, `Closest`,
# `Fillet`, `Chamfer` and the `Affine`/`Closest` transforms through the
# OpenCASCADE internals when `SetFactory("OpenCASCADE")` is active; under the
# built-in factory each reports the recoverable `<name> only available with
# OpenCASCADE geometry kernel`. `.geo` entities stay GEO records carrying OCC
# geometry records — the kernel endpoints below bind referenced entities into
# `m.occ_kernel` at their existing tags first (`_geo_occ_pipe_*!` for
# vertices/edges/wires/faces, `_geo_occ_bind_*!` for shells and volumes), run
# the shared `Model._model_occ_*` endpoint, then `_model_occ_synchronize!`
# publishes created shapes back as GEO entities.

# Keep kernel tag allocation ahead of every `.geo` namespace it may emit:
# entity maxima per dimension, loop tags for wires, and surface-loop tags
# for the shells a solid binds or creates.
function _geo_occ_prepare_kernel!(m::GeoModel,allocator_state=nothing)
    kernel=m.occ_kernel
    for d in 0:3
        kernel.maxima[d+3]=max(kernel.maxima[d+3],m.next_tag[d+1])
    end
    isempty(m.loops) ||
        (kernel.maxima[2]=max(kernel.maxima[2],maximum(keys(m.loops))))
    isempty(m.surface_loops) ||
        (kernel.maxima[1]=max(kernel.maxima[1],maximum(keys(m.surface_loops))))
    return kernel
end

# Bind a `.geo` `Surface Loop` into the kernel as a shell at its own tag —
# the `(-2, tag)` binding records the GModel loop mapping `synchronize`
# round-trips through `kernel.shell_loops`.
function _geo_occ_bind_shell!(m::GeoModel,tag::Int,kernel,caller)
    haskey(kernel.bindings,(-2,tag)) && return kernel.bindings[(-2,tag)]
    haskey(m.surface_loops,tag) || throw(ArgumentError(
        "$caller: unknown Surface Loop[$tag]"))
    uses=Tuple{UInt64,Bool}[]; face_ids=UInt64[]
    for s in m.surface_loops[tag]
        st=Int(s); st!=0 || throw(ArgumentError(
            "$caller: Surface Loop[$tag] contains an invalid Surface[0]"))
        fid=_geo_occ_pipe_face!(m,abs(st),kernel,caller)
        push!(uses,(fid,st<0)); push!(face_ids,fid)
    end
    isempty(uses) && throw(ArgumentError(
        "$caller: Surface Loop[$tag] is empty"))
    sid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,-2,tag,sid,
        Model._ModelOCCShape(-2,(uses=uses,),face_ids))
    kernel.next_shape=sid+UInt64(1)
    kernel.shell_loops[sid]=tag
    return sid
end

# Bind a `.geo` `Volume` into the kernel as a solid at its own tag — each
# member `Surface Loop` becomes a bound shell.
function _geo_occ_bind_volume!(m::GeoModel,tag::Int,kernel,caller)
    haskey(kernel.bindings,(3,tag)) && return kernel.bindings[(3,tag)]
    haskey(m.volumes,tag) || throw(ArgumentError(
        "$caller: unknown Volume[$tag]"))
    shell_ids=UInt64[]
    for st in m.volumes[tag]
        s=Int(st); s!=0 || throw(ArgumentError(
            "$caller: Volume[$tag] contains an invalid Surface Loop[0]"))
        push!(shell_ids,_geo_occ_bind_shell!(m,abs(s),kernel,caller))
    end
    isempty(shell_ids) && throw(ArgumentError(
        "$caller: Volume[$tag] is empty"))
    vid=only(Model._model_occ_reserve_shapes(kernel,1,caller))
    Model._model_occ_bind_shape!(kernel,3,tag,vid,
        Model._ModelOCCShape(3,(kind=:solid,),shell_ids))
    kernel.next_shape=vid+UInt64(1)
    return vid
end

# Lazily bind one `(dim, tag)` entity into the kernel, dim -1 resolving a
# `Wire`/`Curve Loop` record and dim -2 a `Surface Loop`.
function _geo_occ_bind_dim_tag!(m::GeoModel,dim::Int,tag::Int,kernel,caller)
    d=Int(dim); t=Int(tag)
    d==0 && return _geo_occ_pipe_vertex!(m,t,kernel,caller)
    d==1 && return _geo_occ_pipe_edge!(m,t,kernel,caller)
    d==2 && return _geo_occ_pipe_face!(m,t,kernel,caller)
    d==-1 && return _geo_occ_pipe_wire!(m,t,kernel,caller)
    d==-2 && return _geo_occ_bind_shell!(m,t,kernel,caller)
    d==3 && return _geo_occ_bind_volume!(m,t,kernel,caller)
    throw(ArgumentError("$caller: cannot bind dimension $d into OpenCASCADE"))
end

function _geo_occ_bind_entities!(m::GeoModel,entities,kernel,caller)
    for (dim,tag) in entities
        _geo_occ_bind_dim_tag!(m,dim,tag,kernel,caller)
    end
    return entities
end

# A kernel wire referenced where `.geo` expects a `Curve Loop` gets a GEO
# loop record carrying the wire's signed traversal, like the `Wire`
# statement's own materialization.
function _geo_occ_wire_geo_loop!(m::GeoModel,wire_tag::Int,kernel,caller)
    haskey(m.loops,wire_tag) && return m.loops[wire_tag]
    wid=get(kernel.bindings,(-1,wire_tag),nothing)
    wid===nothing && throw(ArgumentError(
        "$caller: unknown Wire/Curve Loop[$wire_tag]"))
    m.loops[wire_tag]=Int[
        (rev ? -1 : 1)*Model._model_occ_shape_tag(kernel,eid,caller)
        for (eid,rev) in kernel.shapes[wid].geometry.uses]
    return m.loops[wire_tag]
end

# Resolve a `.geo` wire operand (a `Curve Loop` record or an existing OCC
# wire binding) to a kernel wire tag, binding GEO loops lazily.
function _geo_occ_wire_tag!(m::GeoModel,tag,kernel,caller)
    t=abs(Int(tag))
    haskey(kernel.bindings,(-1,t)) && return t
    _geo_occ_pipe_wire!(m,t,kernel,caller)
    return t
end

# Statement tail after an OCC mutation: publish the kernel graph into the
# GEO records and mark the internals changed (the enclosing statement loop
# resynchronizes the allocator from the synced model).
function _geo_occ_finish!(m::GeoModel,context::_GeoNumericContext)
    Model._model_occ_synchronize!(m)
    context.geo_changed=true
    return nothing
end

# Kernel endpoint failures are recoverable `.geo` diagnostics upstream —
# `OCC_Internals` returns `false` after its own error and the grammar emits
# `Could not add <what>` (or the operation-specific wording). Surface the
# endpoint's own message first so both diagnostics appear.
function _geo_occ_fail!(context::_GeoNumericContext,err,what::AbstractString)
    _geo_msg_error!(context,sprint(showerror,err))
    _geo_yyerror!(context,what)
    return nothing
end

# `.geo` definition tags feed `addX(int tag, ...)` — `(0)` and negatives are
# automatic like every other `.geo` definition (the `addX` `tag=-1` path);
# a positive tag binds verbatim.
_geo_occ_requested_tag(tag::Int,zero_literal::Bool)=tag>0 ? tag : -1

# After statements the text observer cannot model (OCC-only creates and the
# topology rewrites below), rebuild live sets and counters from the synced
# model and carry the kernel loop maxima across — upstream `newcl`/`newsl`
# read the OCC loop namespaces.
function _geo_occ_allocator_resync!(allocator_state,m::GeoModel)
    _geo_allocator_resync_model!(allocator_state,m)
    kernel=m.occ_kernel
    allocator_state.occ_curveloop_max=max(
        allocator_state.occ_curveloop_max,kernel.maxima[2])
    allocator_state.occ_surfaceloop_max=max(
        allocator_state.occ_surfaceloop_max,kernel.maxima[1])
    return nothing
end

# Run a kernel transform endpoint (`translate`/`rotate`/`dilate`/`symmetry`)
# on resolved `.geo` entities: bind lazily, transform, sync. Endpoint
# failures report upstream's recoverable "Could not ..." wording.
function _geo_occ_apply_transform!(m::GeoModel,entities,t::Model._AffineTransform,
                                   context::_GeoNumericContext,
                                   caller::AbstractString,what::AbstractString)
    kernel=_geo_occ_prepare_kernel!(m)
    try
        _geo_occ_bind_entities!(m,entities,kernel,caller)
        Model._model_occ_transform!(m,entities,t,caller)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _geo_occ_fail!(context,err,"Could not $what shapes")
        return entities
    end
    _geo_occ_finish!(m,context)
    return entities
end

const _GEO_OCC_TRANSFORM_WHAT=Dict("Translate"=>"translate",
    "Rotate"=>"rotate","Dilate"=>"dilate","Symmetry"=>"apply symmetry to")

# ── `.geo` OCC-only entity definitions ───────────────────────────────────────
# `Name(tag) = {params}` dispatch while the OpenCASCADE factory is active —
# the `(FExpr) = ListOfDouble` actions upstream runs as `OCC_Internals::addX`
# calls. Each definition evaluates its tag and parameters first (the grammar
# actions resolve them before the factory gate fires upstream), reports the
# OCC-namespace duplicate diagnostic, runs the `Model._model_occ_*` endpoint,
# and publishes the new shapes through `_geo_occ_finish!`. Endpoint
# `ArgumentError`s surface the recoverable `Could not add <name>` pair
# upstream emits; any other exception propagates.
function _geo_exec_occ_shape_def!(m::GeoModel,s::AbstractString,
                                  context::_GeoNumericContext,
                                  allocator_state)
    kernel=_geo_occ_prepare_kernel!(m,allocator_state)
    if (mm=match(
            r"^Rectangle\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",s))!==nothing
        caller="execute_geo: Rectangle"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        v=_geo_exec_numeric_values(mm.captures[2],context,"$caller parameters")
        if !(length(v) in (5,6))
            _geo_yyerror!(context,"Rectangle requires 5 or 6 parameters")
            return
        end
        _geo_dup_entity!(context,
            haskey(m.surfaces,tag)||haskey(kernel.bindings,(2,tag)),
            :surface,"rectangle",tag,zero_literal,allocator_state) && return
        try
            Model._model_occ_add_rectangle!(m,v[1],v[2],v[3],v[4],v[5];
                tag=_geo_occ_requested_tag(tag,zero_literal),
                rounded_radius=length(v)==6 ? v[6] : 0.0)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add rectangle")
            return
        end
        _geo_occ_finish!(m,context)
        return
    end
    if (mm=match(
            r"^Disk\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",s))!==nothing
        caller="execute_geo: Disk"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        v=_geo_exec_numeric_values(mm.captures[2],context,"$caller parameters")
        if !(length(v) in (4,5))
            _geo_yyerror!(context,"Disk requires 4 or 5 parameters")
            return
        end
        _geo_dup_entity!(context,
            haskey(m.surfaces,tag)||haskey(kernel.bindings,(2,tag)),
            :surface,"disk",tag,zero_literal,allocator_state) && return
        try
            Model._model_occ_add_disk!(m,v[1],v[2],v[3],v[4],
                length(v)==5 ? v[5] : v[4];
                tag=_geo_occ_requested_tag(tag,zero_literal))
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add disk")
            return
        end
        _geo_occ_finish!(m,context)
        return
    end
    if (mm=match(
            r"^Wedge\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",s))!==nothing
        caller="execute_geo: Wedge"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        v=_geo_exec_numeric_values(mm.captures[2],context,"$caller parameters")
        if !(length(v) in (6,7))
            _geo_yyerror!(context,"Wedge requires 6 or 7 parameters")
            return
        end
        _geo_dup_entity!(context,
            haskey(m.volumes,tag)||haskey(kernel.bindings,(3,tag)),
            :volume,"wedge",tag,zero_literal,allocator_state) && return
        try
            Model._model_occ_add_wedge!(m,v[1],v[2],v[3],v[4],v[5],v[6];
                tag=_geo_occ_requested_tag(tag,zero_literal),
                ltx=length(v)==7 ? v[7] : 0.0)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add wedge")
            return
        end
        _geo_occ_finish!(m,context)
        return
    end
    if (mm=match(
            r"^Wire\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",s))!==nothing
        caller="execute_geo: Wire"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        curves=_geo_exec_entity_rhs_tags(
            mm.captures[2],context,"$caller curves";
            signed=true,allow_empty=true)
        _geo_dup_entity!(context,
            haskey(m.loops,tag)||haskey(kernel.bindings,(-1,tag)),
            :curve_loop,"wire",tag,zero_literal,allocator_state) && return
        allocated=try
            for c in curves
                _geo_occ_pipe_edge!(m,abs(c),kernel,caller)
            end
            Model._model_occ_add_wire!(m,curves;
                tag=_geo_occ_requested_tag(tag,zero_literal))
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add wire")
            return
        end
        _geo_occ_finish!(m,context)
        # A kernel wire doubles as a `.geo` curve-loop record — the same
        # signed-traversal list `_geo_occ_wire_geo_loop!` materializes — so
        # `Plane Surface`/`Surface`/`Volume` references see it.
        _geo_occ_wire_geo_loop!(m,allocated,kernel,caller)
        return
    end
    if (mm=match(
            r"^ThickSolid\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}\s*;$",s))!==nothing
        caller="execute_geo: ThickSolid"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        v=_geo_exec_numeric_values(mm.captures[2],context,"$caller parameters")
        if length(v)<2
            _geo_yyerror!(context,
                "ThickSolid requires a volume tag and an offset")
            return
        end
        _geo_dup_entity!(context,
            haskey(m.volumes,tag)||haskey(kernel.bindings,(3,tag)),
            :volume,"thick solid",tag,zero_literal,allocator_state) && return
        volume_tag=_geo_signed_gmsh_int_value(v[1],"$caller volume")
        excludes=Int[_geo_signed_gmsh_int_value(x,"$caller exclude surface")
                     for x in v[3:end]]
        try
            _geo_occ_bind_volume!(m,abs(volume_tag),kernel,caller)
            Model._model_occ_thick_solid(m,
                _geo_occ_requested_tag(tag,zero_literal),abs(volume_tag),
                excludes,v[2])
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add thick solid")
            return
        end
        _geo_occ_finish!(m,context)
        return
    end
    if (mm=match(
            r"^(Ruled\s+)?ThruSections\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",s))!==nothing
        ruled=mm.captures[1]!==nothing
        caller="execute_geo: $(ruled ? "Ruled ThruSections" : "ThruSections")"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[2],context,"$caller tag")
        wires=_geo_exec_entity_rhs_tags(
            mm.captures[3],context,"$caller wires";
            abs_refs=true,allow_empty=true)
        _geo_dup_entity!(context,
            haskey(m.volumes,tag)||haskey(m.surfaces,tag)||
            haskey(kernel.bindings,(3,tag))||haskey(kernel.bindings,(2,tag)),
            :volume,"thru sections",tag,zero_literal,allocator_state) && return
        try
            wire_tags=Int[_geo_occ_wire_tag!(m,wt,kernel,caller) for wt in wires]
            Model._model_occ_thru_sections!(m,
                _geo_occ_requested_tag(tag,zero_literal),wire_tags;
                make_ruled=ruled)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add thru sections")
            return
        end
        _geo_occ_finish!(m,context)
        return
    end
    if (mm=match(
            r"^(BSpline|Bezier)\s+Surface\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",s))!==nothing
        family=mm.captures[1]=="BSpline" ? :bspline : :bezier
        caller="execute_geo: $(mm.captures[1]) Surface"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[2],context,"$caller tag")
        # `{{row},{row},...}` — upstream accepts a `ListOfLists` of control
        # point rows (each row a `{tags}` group) as well as a flat list.
        raw=String(strip(mm.captures[3]))
        points=Int[]
        if startswith(raw,"{") && startswith(strip(raw[2:end]),"{")
            inner=String(strip(raw[2:prevind(raw,end)]))
            while !isempty(inner)
                startswith(inner,"{") || _geo_selector_abort(inner)
                (row,inner)=_geo_balanced_group(inner,caller)
                append!(points,_geo_exec_entity_rhs_tags(
                    "{"*row*"}",context,"$caller control points";
                    abs_refs=true,allow_empty=true))
                inner=String(strip(inner))
                startswith(inner,",") &&
                    (inner=String(strip(inner[2:end])))
            end
        else
            append!(points,_geo_exec_entity_rhs_tags(
                raw,context,"$caller control points";
                abs_refs=true,allow_empty=true))
        end
        _geo_dup_entity!(context,
            haskey(m.surfaces,tag)||haskey(kernel.bindings,(2,tag)),
            :surface,"$(family) surface",tag,zero_literal,allocator_state) && return
        # `addBSplineSurface(tag, pointTags)` derives the u count as
        # `sqrt(|points|)` — a non-square grid is an upstream addX failure.
        nu=isqrt(length(points))
        if nu*nu!=length(points)
            _geo_msg_error!(context,
                "Number of control points $(length(points)) is not a square")
            _geo_yyerror!(context,"Could not add $(family) surface")
            return
        end
        try
            for p in points
                _geo_occ_pipe_vertex!(m,p,kernel,caller)
            end
            family===:bspline ?
                Model._model_occ_add_bspline_surface!(m,points,nu;
                    tag=_geo_occ_requested_tag(tag,zero_literal)) :
                Model._model_occ_add_bezier_surface!(m,points,nu;
                    tag=_geo_occ_requested_tag(tag,zero_literal))
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add $(family) surface")
            return
        end
        _geo_occ_finish!(m,context)
        return
    end
    throw(ArgumentError(
        "execute_geo: unsupported OpenCASCADE entity definition: $s"))
end

# `ThruSections {wires}` / `Ruled ThruSections {wires}` (term and statement
# forms) — `tThruSections ListOfDouble` runs `addThruSections(-1, wires)`
# upstream, returning the created `(dim, tag)` entities.
function _geo_occ_thru_sections_run(m::GeoModel,wire_tags,
                                    context::_GeoNumericContext,
                                    caller::AbstractString,ruled::Bool)
    kernel=_geo_occ_prepare_kernel!(m)
    out=try
        tags=Int[_geo_occ_wire_tag!(m,wt,kernel,caller) for wt in wire_tags]
        Model._model_occ_thru_sections!(m,-1,tags;make_ruled=ruled)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _geo_occ_fail!(context,err,"Could not add thru sections")
        return Tuple{Int,Int}[]
    end
    _geo_occ_finish!(m,context)
    return Tuple{Int,Int}[e for e in out]
end

# `Fillet{regions}{edges}{radii}` / `Chamfer{regions}{edges}{surfaces}
# {distances}` shared by the term and shape-list-element paths — returns the
# output `(dim, tag)` entities, or `nothing` after a recoverable failure.
function _geo_occ_corner_run(m::GeoModel,name::AbstractString,groups,
                             context::_GeoNumericContext,
                             caller::AbstractString)
    kernel=_geo_occ_prepare_kernel!(m)
    vols=_geo_exec_entity_rhs_tags(
        "{"*groups[1]*"}",context,"$caller volumes";abs_refs=true)
    cs=_geo_exec_entity_rhs_tags(
        "{"*groups[2]*"}",context,"$caller curves";abs_refs=true)
    out=if name=="Fillet"
        radii=_geo_numeric_list_values(
            "{"*groups[3]*"}",context,"$caller radii";depth=1)
        length(radii)==1 && (radii=fill(radii[1],length(cs)))
        try
            for v in vols
                _geo_occ_bind_volume!(m,v,kernel,caller)
            end
            for c in cs
                _geo_occ_pipe_edge!(m,c,kernel,caller)
            end
            Model._model_occ_fillet(m,vols,cs,radii)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not fillet shapes")
            return nothing
        end
    else
        # `Chamfer{regions}{edges}{surfaces}{distances}` — the surface list
        # is a hint upstream resolves against the solid; Tessella derives
        # the adjacent faces from the edge map and only needs the group
        # evaluated.
        _geo_exec_entity_rhs_tags(
            "{"*groups[3]*"}",context,"$caller surfaces";abs_refs=true)
        dist=_geo_numeric_list_values(
            "{"*groups[4]*"}",context,"$caller distances";depth=1)
        if !(length(dist) in (1,2))
            _geo_msg_error!(context,
                "Chamfer distances must have 1 or 2 elements")
            _geo_yyerror!(context,"Could not chamfer shapes")
            return nothing
        end
        try
            for v in vols
                _geo_occ_bind_volume!(m,v,kernel,caller)
            end
            for c in cs
                _geo_occ_pipe_edge!(m,c,kernel,caller)
            end
            Model._model_occ_chamfer(
                m,vols,cs,dist[1],length(dist)==2 ? dist[2] : dist[1])
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not chamfer shapes")
            return nothing
        end
    end
    _geo_occ_finish!(m,context)
    return Tuple{Int,Int}[e for e in out]
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
    if allocator_state!==nothing && allocator_state.factory===:opencascade
        _geo_occ_apply_transform!(
            m,entities,t,context,caller,_GEO_OCC_TRANSFORM_WHAT[kind])
        return nothing
    end
    transform_entities!(m,t,entities;caller=caller,auto_coherence=_geo_auto_coherence(context))
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
    acc=nothing;sign=1.0;expect_term=true
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
            # A lone `VExpr_Single` is assigned, not accumulated — its signed
            # zeros survive (`{-0,...}` keeps -0.0); later terms still add.
            acc=acc===nothing ? sign .* evaluated :
                                acc .+ sign .* evaluated
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

# OCC `Circle`/`Ellipse` accept a center+parameters form whose `{...}` arity
# differs from the point-tag grammar (`Circle`: 4, 6, 9 or 12 numbers for
# `{x,y,z,r[,a1,a2[,zAxis...,xAxis...]]}`; `Ellipse`: 5, 7, 10 or 13 with the
# `rx,ry` pair). Upstream disambiguates purely by parameter count, so scan
# the raw RHS before the tag parser runs: returns `(values, rest)` with
# `values===nothing` whenever `raw` is not a plain `{...}` list of one of the
# `allowed` arities — the caller then falls back to `_geo_circle_rhs`.
function _geo_occ_center_params(raw::AbstractString,
                                context::_GeoNumericContext,
                                caller::AbstractString,
                                allowed)
    s=String(strip(raw))
    startswith(s,"{") || return (nothing,"")
    (group,rest)=_geo_balanced_group(s,caller)
    isempty(rest) || return (nothing,rest)
    parts=_geo_split_top_commas(group,caller)
    length(parts) in allowed || return (nothing,rest)
    values=Float64[_geo_eval_numeric(p,context,caller) for p in parts]
    return (values,rest)
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
    cx=0.0; cy=0.0; cz=0.0
    @inbounds for i in 1:n
        point=points[i]
        cx+=point[1]; cy+=point[2]; cz+=point[3]
    end
    cx/=n; cy/=n; cz/=n
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
            Model._native_geo_kernel_declare_transfinite!(
                m,curve,m.meshing.transfinite_curves[curve])
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
                Model._native_geo_kernel_declare_transfinite!(
                    m,curve,m.meshing.transfinite_curves[curve])
            end
            continue
        end
        haskey(m.curves,j) || continue
        m.meshing.transfinite_curves[j]=(
            num_nodes=ncount,kind=kind,coef=coef,reversed=sd<0)
        Model._native_geo_kernel_declare_transfinite!(
            m,j,m.meshing.transfinite_curves[j])
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
# `tSplit tCurve '{' FExpr '}' tPoint '{' RecursiveListOfDouble '}'` — parse
# the two argument groups, returning `(curve_tag, point_tags, rest)`. The
# deprecated `Split Curve(c) {...}` form parses identically after its
# warning, like upstream.
function _geo_split_curve_args(s::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString)
    s=String(strip(s))
    match(r"^Curve\b",s)===nothing && _geo_selector_abort(s)
    s=String(strip(s[nextind(s,firstindex(s),5):end]))
    deprecated=startswith(s,"(")
    if deprecated
        _geo_yywarn!(context,
            "'Split Curve(c) {...}' is deprecated: use 'Split Curve {c} " *
            "Point {...}' instead")
        (inner,s)=_geo_balanced_paren(s,caller)
    else
        startswith(s,"{") || _geo_selector_abort(s)
        (inner,s)=_geo_balanced_group(s,caller)
    end
    curve_tag=_geo_constraint_int(_geo_eval_numeric(
        inner,context,"$caller curve"),caller,"curve tag")
    s=String(strip(s))
    if !deprecated
        match(r"^Point\b",s)===nothing && _geo_selector_abort(s)
        s=String(strip(s[nextind(s,firstindex(s),5):end]))
    end
    startswith(s,"{") || _geo_selector_abort(s)
    (plist,s)=_geo_balanced_group(s,caller)
    point_tags=Int[_geo_constraint_int(value,caller,"point tag") for
        value in _geo_numeric_list_values(
            "{"*plist*"}",context,"$caller points";depth=1)]
    return (curve_tag,point_tags,String(strip(s)))
end

# `GEO_Internals::splitCurve` — break a LINE/SPLN/BSPLN curve at the listed
# control points: each run of consecutive control points becomes a new
# same-type curve, curve loops and physical-group memberships are rewired
# in place (`-c` expands to the new tags in reversed, negated order), and
# the original is deleted. Break vertices not on the control-point list
# are ignored; a periodic curve (beg == end) wraps its first and last
# runs into a single curve. Segment tags come from `NEWCURVE()` — `NEWREG()`
# under `Geometry.OldNewReg` (the default), the per-dimension curve counter
# at 0 — and physical memberships update the raw group records, which the
# lazy view resync then re-derives. Returns the new curve tags in creation
# order.
# The two upstream `SplitCurve` failure modes are recoverable `Msg::Error`s
# (followed by the grammar's `yymsg(0, "Could not split curve")`), not aborts:
# a missing curve or a non-LINE/SPLN/BSPLN type reports and yields nothing.
# `nothing` marks that path; `caller` formats only the cause messages.
function _geo_split_curve_impl!(m::GeoModel,curve_tag::Int,
                                point_tags::Vector{Int},
                                context::_GeoNumericContext,
                                caller::AbstractString)
    if !haskey(m.curves,curve_tag)
        _geo_msg_error!(context,"Unknown curve $curve_tag to split")
        return nothing
    end
    kind=get(m.curve_types,curve_tag,:line)
    if !(kind in (:line,:spline,:bspline))
        # `c->Typ` prints the `MSH_SEGM_*` code: CIRC 202, ELLI 204,
        # NURBS 208, BEZIER 209 (LINE/SPLN/BSPLN never reach this branch).
        type_id=kind==:circle ? 202 :
            kind==:ellipse ? 204 :
            kind==:nurbs ? 208 :
            kind==:bezier ? 209 : 0
        _geo_msg_error!(context,
            "Cannot split curve $curve_tag with type $type_id")
        return nothing
    end
    controls=kind==:line ? collect(m.curves[curve_tag]) :
             copy(m.curve_control_points[curve_tag])
    breaks=Set(point_tags)
    a,b=m.curves[curve_tag]
    periodic=a==b
    first_periodic=true
    last_periodic=false
    runs=Vector{Int}[]
    run=Int[]
    i=1;n=length(controls)
    while i<=n
        pv=controls[i]
        push!(run,pv)
        if pv in breaks && length(run)>1
            last_periodic && break
            (periodic && first_periodic) || push!(runs,copy(run))
            first_periodic=false
            run=[pv]
        end
        if i==n && periodic && !first_periodic
            # Upstream's `i = 0` followed by the loop `i++` resumes at the
            # second control point (the first repeats the last on a
            # periodic list).
            i=1
            last_periodic=true
        end
        i+=1
    end
    length(run)>1 && push!(runs,copy(run))
    old_new_reg=!iszero(something(_geo_option_number(
        context,"Geometry",0,"OldNewReg"),1.0))
    new_tags=Int[]
    for nodes in runs
        tag=old_new_reg ? _geo_newreg_alloc!(m,1,caller) :
            _alloc_tag!(m,1,0,caller)
        _geo_store_curve!(m,tag,nodes;
            kind=kind==:line ? nothing : kind,
            control_points=kind!=:line)
        push!(new_tags,tag)
    end
    for (loop,signed_curves) in m.loops
        any(sc->abs(sc)==curve_tag,signed_curves) || continue
        rebuilt=Int[]
        for sc in signed_curves
            if sc==curve_tag
                append!(rebuilt,new_tags)
            elseif sc==-curve_tag
                for tag in Iterators.reverse(new_tags)
                    push!(rebuilt,-tag)
                end
            else
                push!(rebuilt,sc)
            end
        end
        m.loops[loop]=rebuilt
    end
    # Upstream rewrites the physical group's `Entities` list in place: the
    # raw records are the canonical store here, and the derived view follows
    # at the next sync (memberships stay immediately visible either way).
    carried=get(m.entity_physicals,(1,curve_tag),Int[])
    for ((dim,ptag),members) in context.raw_physicals
        dim==1 || continue
        any(==(curve_tag),members) || continue
        rebuilt=Int[]
        for member in members
            member==curve_tag ? append!(rebuilt,new_tags) :
                push!(rebuilt,member)
        end
        context.raw_physicals[(dim,ptag)]=rebuilt
    end
    for ((dim,ptag),members) in m.physical
        dim==1 || continue
        any(==(curve_tag),members) || continue
        rebuilt=Int[]
        for member in members
            member==curve_tag ? append!(rebuilt,new_tags) :
                push!(rebuilt,member)
        end
        m.physical[(dim,ptag)]=rebuilt
    end
    for tag in new_tags
        m.entity_physicals[(1,tag)]=copy(carried)
    end
    _geo_delete_entities!(m,[(1,curve_tag)])
    return new_tags
end

# The `SplitCurve` call itself — separated from argument parsing so the
# `FExpr_Multi` Transform path can check the terminator before committing
# the side effect (upstream only runs `splitCurve` when the production
# reduces).
function _geo_split_curve_run!(m::GeoModel,curve_tag::Int,
                               point_tags::Vector{Int},
                               context::_GeoNumericContext,
                               caller::AbstractString,allocator_state)
    if allocator_state!==nothing && allocator_state.factory==:opencascade
        _geo_yyerror!(context,
            "Split Curve not available with OpenCASCADE geometry kernel")
        return Int[]
    end
    new_tags=_geo_split_curve_impl!(m,curve_tag,point_tags,context,caller)
    # `GEO_Internals::splitCurve` marks the internals changed even when the
    # operation fails — the `_changed = true` runs before `return ok`.
    context.geo_changed=true
    if new_tags===nothing
        _geo_yyerror!(context,"Could not split curve")
        return Int[]
    end
    return new_tags
end

function _geo_exec_split_curve!(m::GeoModel,tail::AbstractString,
                                context::_GeoNumericContext,
                                caller::AbstractString,allocator_state)
    (curve_tag,point_tags,rest)=_geo_split_curve_args(tail,context,caller)
    return (_geo_split_curve_run!(
        m,curve_tag,point_tags,context,caller,allocator_state),rest)
end

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
        # `GModel::remove` unbinds every OCC-resident entity it deletes —
        # leaving the binding would resurrect the tag at the next
        # `_model_occ_synchronize!`. `_geo_delete_entities!` already
        # cascaded boundary children, so each binding drops directly
        # (the `_model_occ_shape_is_owned` guard does not apply to an
        # explicit `.geo` delete); unreachable shapes are collected.
        kernel=m.occ_kernel
        bound=sort!(Tuple{Int,Int}[e for e in removed
                                  if haskey(kernel.bindings,e)];
                    by=e->-e[1])
        for entry in bound
            id=kernel.bindings[entry]
            Model._model_occ_entity_map_unbind!(
                kernel.entity_maps[first(entry)+3],last(entry))
            delete!(kernel.bindings,entry)
            delete!(kernel.shape_tags,id)
            push!(kernel.removed,entry)
            kernel.maxima[first(entry)+3]=maximum(
                (t for (d,t) in keys(kernel.bindings) if d==first(entry));
                init=0)
            kernel.changed=true
        end
        isempty(bound) || Model._model_occ_collect_shapes!(kernel)
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
        context.mesh_parts=empty(context.mesh_parts)
        context.mesh_geometry_parts=nothing
        _geo_forget_mesh_products!(m,context;reset_tags=true)
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
        context.mesh_parts=empty(context.mesh_parts)
        context.mesh_geometry_parts=nothing
        _geo_forget_mesh_products!(m,context;reset_tags=true)
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
        # result-list behavior, the generation-scoped option mirrors, and
        # every tracked option write.
        m.meshing.order=1
        m.meshing.transfinite_tri=0
        m.meshing.flexible_transfinite=false
        m.meshing.lc_factor=1.0
        m.meshing.aniso_max=1e33
        m.meshing.recombine_all=false
        m.meshing.recombine_algo=1
        m.meshing.lc_extend_from_boundary=1
        context.extrude_return_lateral=true
        empty!(context.option_strings);empty!(context.option_numbers)
        empty!(context.option_colors)
        return 0
    elseif name=="Meshes"
        cleared_authority=Model._model_mesh_cleared_authority(m.meshing.mesh_authority)
        # `GModel::deleteMesh` — drops the mid-file mesh and its ownership.
        context.mesh=nothing
        context.mesh_node_owner=empty(context.mesh_node_owner)
        context.mesh_parts=empty(context.mesh_parts)
        context.mesh_geometry_parts=nothing
        _geo_forget_mesh_products!(m,context)
        m.meshing.mesh_authority=cleared_authority
        context.native_element_lookup_receipt=nothing
        empty!(m.curve_params);empty!(m.meshing.attached)
        for record in values(m.discrete)
            empty!(record.node_tags);record.node_coords=zeros(3,0)
            record.node_params=zeros(0,0)
            empty!(record.element_types);empty!(record.element_tags);empty!(record.element_nodes)
            empty!(record.aux_params)
        end
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
        if allocator_state===nothing || allocator_state.factory!==:opencascade
            Model._native_geo_kernel_admit_point!(m,t)
        end
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
        occ=allocator_state!==nothing && allocator_state.factory===:opencascade
        _geo_store_curve!(m,t,found;control_points=length(found)>2,
            geometry=occ ? _geo_occ_line_geometry(m,found) : nothing,native_geo=!occ)
        ok || _geo_yyerror!(context,"Could not add line")
        return
    elseif (mm=match(
            r"^Circle\s*\(\s*(.*?)\s*\)\s*=\s*(.*?)\s*;$",
            line)) !== nothing
        caller="execute_geo: Circle"
        (tag,zero_literal)=_geo_exec_def_tag(
            mm.captures[1],context,"$caller tag")
        occ=allocator_state!==nothing && allocator_state.factory===:opencascade
        if occ
            # `Circle(tag) = {x,y,z,r[,a1,a2[,zAxis...,xAxis...]]}` — the OCC
            # `addCircle` center+radius form (4, 6, 9 or 12 parameters);
            # upstream selects it purely by parameter count.
            (vals,)=_geo_occ_center_params(
                mm.captures[2],context,caller,(4,6,9,12))
            if vals!==nothing
                kernel=_geo_occ_prepare_kernel!(m)
                _geo_dup_entity!(context,
                    haskey(m.curves,trunc(Int,tag)) ||
                    haskey(kernel.bindings,(1,tag)),
                    :curve,"circle",tag,zero_literal,allocator_state) && return
                try
                    Model._model_occ_add_circle!(m,
                        vals[1],vals[2],vals[3],vals[4];
                        tag=_geo_occ_requested_tag(tag,zero_literal),
                        angle1=length(vals)>=6 ? vals[5] : 0.0,
                        angle2=length(vals)>=6 ? vals[6] : 2pi,
                        zaxis=length(vals)>=9 ? vals[7:9] : nothing,
                        xaxis=length(vals)>=12 ? vals[10:12] : nothing)
                catch err
                    err isa InterruptException && rethrow()
                    err isa ArgumentError || rethrow()
                    _geo_occ_fail!(context,err,"Could not add circle")
                    return
                end
                _geo_occ_finish!(m,context)
                return
            end
        end
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
        if occ
            for point in points
                if !haskey(m.points,point)
                    _geo_msg_error!(context,"Unknown OpenCASCADE point with tag $point")
                    _geo_yyerror!(context,"Could not add circle")
                    return
                end
            end
            geometry,error=Model._model_occ_circle_arc_geometry(
                m.points[points[1]],m.points[points[2]],m.points[points[3]])
            if error!==nothing
                _geo_msg_error!(context,error)
                _geo_yyerror!(context,"Could not add circle")
                return
            end
            _geo_store_curve!(m,t,points;kind=:circle,control_points=true,geometry)
            return
        end
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
        occ=allocator_state!==nothing && allocator_state.factory===:opencascade
        if occ
            # `Ellipse(tag) = {x,y,z,rx,ry[,a1,a2[,zAxis...,xAxis...]]}` —
            # the OCC `addEllipse` center form (5, 7, 10 or 13 parameters),
            # again selected purely by parameter count upstream.
            (vals,)=_geo_occ_center_params(
                mm.captures[2],context,caller,(5,7,10,13))
            if vals!==nothing
                kernel=_geo_occ_prepare_kernel!(m)
                _geo_dup_entity!(context,
                    haskey(m.curves,trunc(Int,tag)) ||
                    haskey(kernel.bindings,(1,tag)),
                    :curve,"ellipse",tag,zero_literal,allocator_state) &&
                    return
                try
                    Model._model_occ_add_ellipse!(m,
                        vals[1],vals[2],vals[3],vals[4],vals[5];
                        tag=_geo_occ_requested_tag(tag,zero_literal),
                        angle1=length(vals)>=7 ? vals[6] : 0.0,
                        angle2=length(vals)>=7 ? vals[7] : 2pi,
                        zaxis=length(vals)>=10 ? vals[8:10] : nothing,
                        xaxis=length(vals)>=13 ? vals[11:13] : nothing)
                catch err
                    err isa InterruptException && rethrow()
                    err isa ArgumentError || rethrow()
                    _geo_occ_fail!(context,err,"Could not add ellipse")
                    return
                end
                _geo_occ_finish!(m,context)
                return
            end
        end
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
        if ok && allocator_state!==nothing && allocator_state.factory===:opencascade
            geometry=_geo_occ_axis_plane_geometry(m,loops)
            geometry===nothing || (m.surface_geometry[t]=geometry)
        end
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
        if !(length(values) in (7,8))
            _geo_yyerror!(context,"Cylinder requires 7 or 8 parameters")
            return
        end
        _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"cylinder",
                         tag,zero_literal,allocator_state) &&
            return
        if length(values)==8
            # `{x,y,z,dx,dy,dz,r,angle}` — a partial sweep the native
            # encoding cannot store; the OCC endpoint builds it.
            kernel=_geo_occ_prepare_kernel!(m)
            try
                Model._model_occ_add_cylinder!(m,values[1:7]...;
                    tag=_geo_occ_requested_tag(tag,zero_literal),
                    angle=values[8])
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                _geo_occ_fail!(context,err,"Could not add cylinder")
                return
            end
            _geo_occ_finish!(m,context)
            return
        end
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
            if length(values)==4
                add_sphere!(m,values[1],values[2],values[3],values[4];
                            tag=tag,_zero_literal=zero_literal)
            else
                # `{x,y,z,r,a1[,a2[,a3]]}` — the angle-bound forms bind a
                # partial sphere through the OCC endpoint (`add_sphere!`
                # only encodes the full sweep).
                kernel=_geo_occ_prepare_kernel!(m)
                try
                    Model._model_occ_add_sphere!(m,values[1],values[2],
                        values[3],values[4];
                        tag=_geo_occ_requested_tag(tag,zero_literal),
                        angle1=values[5],
                        angle2=length(values)>=6 ? values[6] : pi/2,
                        angle3=length(values)==7 ? values[7] : 2pi)
                catch err
                    err isa InterruptException && rethrow()
                    err isa ArgumentError || rethrow()
                    _geo_occ_fail!(context,err,"Could not add sphere")
                    return
                end
                _geo_occ_finish!(m,context)
            end
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
        if !(length(values) in (8,9))
            _geo_yyerror!(context,"Cone requires 8 or 9 parameters")
            return
        end
        _geo_dup_entity!(context,haskey(m.volumes,tag),:volume,"cone",
                         tag,zero_literal,allocator_state) &&
            return
        if length(values)==9
            # `{x,y,z,dx,dy,dz,r1,r2,angle}` — a partial cone sweep, built
            # through the OCC endpoint like `Cylinder`'s angle form.
            kernel=_geo_occ_prepare_kernel!(m)
            try
                Model._model_occ_add_cone!(m,values[1:8]...;
                    tag=_geo_occ_requested_tag(tag,zero_literal),
                    angle=values[9])
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                _geo_occ_fail!(context,err,"Could not add cone")
                return
            end
            _geo_occ_finish!(m,context)
            return
        end
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
        kernel=_geo_occ_prepare_kernel!(m)
        _geo_dup_entity!(context,
            haskey(m.volumes,tag) || haskey(kernel.bindings,(3,tag)),
            :volume,"torus",tag,zero_literal,allocator_state) &&
            return
        # The `.geo` Torus is an OCC-only entity; the kernel endpoint owns the
        # sweep topology so OCC ops (Boundary, Duplicata, Fillet, …) consume
        # it. `angle` stays optional like upstream `addTorus`.
        try
            Model._model_occ_add_torus!(m,values[1],values[2],values[3],
                values[4],values[5];
                tag=_geo_occ_requested_tag(tag,zero_literal),
                angle=length(values)==6 ? values[6] : 2pi)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            _geo_occ_fail!(context,err,"Could not add torus")
            return
        end
        _geo_occ_finish!(m,context)
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
        _geo_store_curve!(m,t,found;control_points=length(found)>2,
            native_geo=allocator_state===nothing || allocator_state.factory!==:opencascade)
        ok || _geo_yyerror!(context,"Could not add line")
        return
    elseif (mm=match(
            Regex("^(Rectangle|Disk|Wedge|ThickSolid|ThruSections|Wire)\\s*\\(|" *
                  "^(Ruled\\s+)ThruSections\\s*\\(|" *
                  "^(BSpline|Bezier)\\s+Surface\\s*\\("),
            line)) !== nothing
        # OCC-only entity definitions: under the built-in factory Gmsh emits
        # `<name> only available with OpenCASCADE geometry kernel`; under
        # OpenCASCADE the statement runs its `OCC_Internals::addX` analog.
        name=mm.captures[1]!==nothing ? mm.captures[1] :
            (mm.captures[2]!==nothing ? "ThruSections" :
             "$(mm.captures[3]) surface")
        if allocator_state!==nothing && allocator_state.factory==:opencascade
            _geo_exec_occ_shape_def!(m,line,context,allocator_state)
            return
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
        splitcaller="execute_geo: Split"
        (new_tags,rest)=_geo_exec_split_curve!(
            m,String(strip(line[nextind(line,firstindex(line),5):end])),
            context,splitcaller,allocator_state)
        _geo_require_statement_end(rest,splitcaller)
        return
    elseif match(
            r"^(Closest|Fillet|Chamfer|Ruled\s+ThruSections|ThruSections)\s*(?:\{|[\w(+-])",
            line)!==nothing
        # `Closest{x,y,z}{shapes}`, `Fillet{regions}{edges}{radii}`,
        # `Chamfer{regions}{edges}{surfaces}{distances}` and
        # `(Ruled) ThruSections{wires}` are Transform productions upstream:
        # the statement form discards the produced list. The term executor
        # owns both the built-in OCC-gated diagnostic and the OpenCASCADE
        # kernel routing. An unbraced `ThruSections` operand list (a bare
        # `FExpr` upstream) is wrapped into the `{wires}` group the term
        # handler parses.
        s=String(strip(line))
        endswith(s,";") &&
            (s=String(strip(s[firstindex(s):prevind(s,end)])))
        if (tm=match(r"^(Ruled\s+)?ThruSections\b",s))!==nothing &&
                !startswith(strip(s[nextind(s,firstindex(s),
                                          ncodeunits(tm.match)):end]),"{")
            tail=strip(s[nextind(s,firstindex(s),ncodeunits(tm.match)):end])
            s="$(tm.match){$(tail)}"
        end
        _geo_exec_transform_term(m,s,context,allocator_state)
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
        # `Affine{matrix}{shapes};` — a Transform production like the
        # keyword transforms above: the statement discards the list.
        s=String(strip(line))
        endswith(s,";") &&
            (s=String(strip(s[firstindex(s):prevind(s,end)])))
        _geo_exec_transform_term(m,s,context,allocator_state)
        return
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
            # Upstream stores `(int)val` verbatim — 0 disables and negatives
            # behave as unset, deferring to lcExtendFromBoundary at mesh time.
            m.meshing.size_from_boundary[(2,sd)]=number
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
                tag->begin
                    m.meshing.reverse[(dim,tag)]=true
                    dim==2 && haskey(m.meshing.orientation_geo_reverse,tag) &&
                        (m.meshing.orientation_geo_reverse[tag]=true)
                    nothing
                end)
        end
        return
    elseif (mm=match(
            r"^RelocateMesh\s+(Point|Curve|Line|Surface)\s*(.*?)\s*;$",
            line)) !== nothing
        _geo_exec_relocate_mesh!(m,mm.captures[1],mm.captures[2],context)
        return
    elseif (mm=match(
            r"^ReorientMesh\s+Volume\s*(.*?)\s*;$",line)) !== nothing
        _geo_exec_reorient_mesh!(m,mm.captures[1],context)
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
        return _geo_exec_mesh_statement!(m,line,context,allocator_state)
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
                caller="execute_geo: Coherence Point",
                remove_degenerated=_geo_auto_coherence(context)>1)
        elseif mm.captures[1] != "Mesh"
            coherence!(m;remove_degenerated=_geo_auto_coherence(context)>1)
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
        context.native_context_forced_bbox=false
        _geo_sync_physical_view_if_changed!(m,context)
        _geo_sync_native_context!(m,context)
        return
    elseif (mm=match(r"^BoundingBox\s*\{(.*?)\}\s*;?\s*$",line)) !== nothing
        values=_geo_numeric_list_values("{"*mm.captures[1]*"}",context,
            "execute_geo: BoundingBox")
        length(values)==6 || throw(ArgumentError(
            "execute_geo: BoundingBox takes six coordinates"))
        context.native_context_forced_bbox=true
        tolerance=something(_geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)
        context.native_context_lc=Model._gmsh_bbox_characteristic_length(
            Tuple(values),tolerance)
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
        context.mesh_parts=empty(context.mesh_parts)
        context.mesh_geometry_parts=nothing
        _geo_forget_mesh_products!(m,context;reset_tags=true)
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
        _geo_exec_mesh_statement!(m,line,context,allocator_state)
        return
    elseif (mm=match(r"^Mesh\s+(.*?)\s*;?\s*$",line)) !== nothing
        d=_geo_eval_numeric(mm.captures[1],context,"execute_geo: Mesh")
        d in (0,1,2,3) || throw(ArgumentError(
            "execute_geo: Mesh dimension must be 0, 1, 2 or 3 (got $d)"))
        # Gmsh synchronizes GEO internals before meshing (when changed) — the
        # physical view and any pending unknown-member warnings materialize.
        # Dimension zero skips element births, but GModel::mesh still runs
        # generation setup, renumbering and queued homology on living cells.
        _geo_sync_physical_view_if_changed!(m,context)
        if d==0
            context.native_element_lookup_receipt=nothing
            context.native_mesh_begin_callback===nothing || context.native_mesh_begin_callback()
            context.mesh,context.mesh_node_owner=
                _geo_mesh_actual_higher(m,0,context,"execute_geo: Mesh")
            return
        end
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
    elseif (mm=match(
            r"^Background\s+Mesh\s+View\s*\[\s*([^\]]*)\s*\]\s*;$",
            line))!==nothing
        _geo_exec_background_mesh_view!(
            mm.captures[1],context,allocator_state)
        return
    elseif match(Regex("^(Combine|Alias|AliasWithOptions|Save\\s+View|" *
        "CutMesh|SplitMesh|SendToServer)\\b"),line)!==nothing
        throw(ArgumentError(
            "execute_geo: post-processing view operations are not supported: $line"))
    elseif match(r"^Levelset\b",line)!==nothing
        throw(ArgumentError(
            "execute_geo: Levelset definitions are not supported"))
    elseif (mm=match(r"^(Homology|Cohomology|Betti)\b(.*)$",line))!==nothing
        # `Homology`/`Cohomology`/`Betti` queue a request like
        # `GModel::addHomologyRequest`; requests run inside the next
        # `Mesh n` (and the implicit final generation), after the mesh
        # exists — the upstream `GModel::mesh` placement.
        _geo_exec_homology_statement!(m,mm.captures[1],mm.captures[2],
                                      context)
        return
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

# `Background Mesh View[i]` — `GModel::setBackgroundMesh(i)` creates a new
# `PostView` field bound to view index `i` and selects it as the background
# field. An out-of-range index reports `Msg::Error("Unknown view %d")` and
# leaves the previous selection untouched.
function _geo_exec_background_mesh_view!(index_src,
                                         context::_GeoNumericContext,
                                         allocator_state)
    caller="execute_geo: Background Mesh View"
    index=_geo_signed_gmsh_int_value(
        _geo_eval_numeric(index_src,context,caller),caller)
    store=get(_GEO_EXEC_VIEWS,context,nothing)
    if store===nothing || !(0<=index<length(store))
        _geo_msg_error!(context,"Unknown view $index")
        return nothing
    end
    fields=context.fields
    fields===nothing && return nothing
    tag=isempty(fields) ? 1 : maximum(keys(fields))+1
    curvature=_geo_option_number(context,"Mesh",0,"MeshSizeFromCurvature")
    fields[tag]=GeoFieldSpec(tag,"PostView",
        Dict{String,String}("ViewIndex"=>string(index)),String[],
        curvature===nothing ? 0 : Int(curvature))
    context.background_field=tag
    if allocator_state!==nothing
        push!(allocator_state.live_fields,tag)
    end
    return nothing
end

# `X Field = {list}` — Gmsh's `String__Index tField tAFFECT ListOfDouble`:
# `Background` takes exactly one tag (`setBackgroundFieldId` stores the raw
# id unvalidated — a dangling id only fails when the field is used) and
# `BoundaryLayer` dedup-appends. Any other prefix is an unknown command.
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
        else
            context.background_field=tags[1]
        end
    elseif name=="BoundaryLayer"
        for tag in _geo_expression_field_tags(list_src,context,caller)
            tag in context.boundary_layer_fields ||
                push!(context.boundary_layer_fields,tag)
        end
    else
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

# Generation-scoped `Mesh.*` options are upstream global CTX state read at
# mesh time; mirror them into the model's meshing attributes before every
# generation so `Mesh n` statements see the current option values. The
# `MeshSizeFactor`/`CharacteristicLengthFactor` alias mirror in
# `_geo_store_option_number!` keeps a single `lcFactor` visible under both
# names, so reading either here picks up writes to the other.
function _geo_sync_mesh_algorithms!(m::GeoModel,context::_GeoNumericContext)
    m.meshing.algorithm_2d=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"Algorithm"),6.0),
        "execute_geo: Mesh.Algorithm")
    algorithm_3d=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"Algorithm3D"),1.0),
        "execute_geo: Mesh.Algorithm3D")
    m.meshing.algorithm_3d=algorithm_3d==2 ? 1 : algorithm_3d
    return nothing
end

function _geo_sync_meshing_options!(m::GeoModel,context::_GeoNumericContext)
    _geo_sync_mesh_algorithms!(m,context)
    m.meshing.flexible_transfinite=!iszero(something(
        _geo_option_number(context,"Mesh",0,"FlexibleTransfinite"),0.0))
    m.meshing.lc_factor=something(
        _geo_option_number(context,"Mesh",0,"CharacteristicLengthFactor"),1.0)
    m.meshing.aniso_max=something(
        _geo_option_number(context,"Mesh",0,"AnisoMax"),1e33)
    m.meshing.recombine_all=!iszero(something(
        _geo_option_number(context,"Mesh",0,"RecombineAll"),0.0))
    m.meshing.recombine_algo=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"RecombinationAlgorithm"),1.0),
        "execute_geo: Mesh.RecombinationAlgorithm")
    _geo_sync_native_surface_controls!(m,context)
    m.meshing.lc_extend_from_boundary=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"MeshSizeExtendFromBoundary"),1.0),
        "execute_geo: Mesh.MeshSizeExtendFromBoundary")
    return nothing
end

function _geo_sync_native_surface_controls!(m::GeoModel,context::_GeoNumericContext)
    m.meshing.recombine_topology_passes=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"RecombineOptimizeTopology"),5.0),
        "execute_geo: Mesh.RecombineOptimizeTopology")
    m.meshing.recombine_node_repositioning=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"RecombineNodeRepositioning"),1.0),
        "execute_geo: Mesh.RecombineNodeRepositioning")
    m.meshing.recombine_minimum_quality=something(
        _geo_option_number(context,"Mesh",0,"RecombineMinimumQuality"),0.01)
    m.meshing.native_smoothing=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"Smoothing"),1.0),
        "execute_geo: Mesh.Smoothing")
    m.meshing.second_order_linear=!iszero(something(
        _geo_option_number(context,"Mesh",0,"SecondOrderLinear"),0.0))
    return nothing
end

# Entity references inside `Distance`/`Restrict`/`Constant` field options —
# the same `(dim, tag) -> Mesh` contract `build_geo_size_field`'s resolver
# documents: a model point becomes a single-node mesh, a curve its stored
# `curve_params` nodes when meshed (upstream samples the edge's mesh
# vertices) else its endpoints, a surface its boundary PSLG mesh.
function _geo_field_entity_mesh(m::GeoModel,dim::Int,name::AbstractString,
                                caller::AbstractString)
    tag=tryparse(Int,String(name))
    if tag===nothing
        numeric=tryparse(Float64,String(name))
        (numeric===nothing || !isinteger(numeric) ||
         numeric<typemin(Int) || numeric>typemax(Int)) && return nothing
        tag=Int(numeric)
    end
    if dim==0
        haskey(m.points,tag) || return nothing
        return Mesh(reshape(collect(Float64.(m.points[tag])),3,1))
    elseif dim==1
        params=get(m.curve_params,tag,nothing)
        if params!==nothing
            coords=Matrix{Float64}(undef,3,length(params))
            for (i,u) in enumerate(params)
                p=_model_curve_part_point(m,tag,Float64(u),caller)
                coords[1,i],coords[2,i],coords[3,i]=p
            end
            return Mesh(coords)
        end
        if haskey(m.curves,tag)
            a,b=m.curves[tag]
            (haskey(m.points,a) && haskey(m.points,b)) || return nothing
            coords=hcat(collect(Float64.(m.points[a])),
                        collect(Float64.(m.points[b])))
            return Mesh(coords;segs=reshape(Int32[1;2],2,1))
        end
        return nothing
    elseif dim==2
        haskey(m.surfaces,tag) || return nothing
        return _model_planar_surface_mesh(
            m,tag,caller;include_embeddings=false)
    end
    return nothing
end

# `context_fields` resolver for the five model/view-backed field kinds —
# `AutomaticMeshSizeField` resolves against the model surfaces through the
# provided entity callback; `PostView` resolves against the `.pos` views the
# execution merged (`Merge`/`Include`), mirroring `GModel`'s view list; the
# other kinds keep the explicit diagnostic.
function _geo_mesh_context_fields(m::GeoModel,
                                  context::_GeoNumericContext,
                                  spec,config,entities,params)
    if lowercase(spec.kind)=="postview"
        fields_mod=getfield(parentmodule(@__MODULE__),:SizeField)
        store=get(_GEO_EXEC_VIEWS,context,nothing)
        (store===nothing || isempty(store)) && throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] references view " *
            "index $(config.view_index) but no views are loaded"))
        result=fields_mod.geo_postview_context(store)(
            spec,config,entities,params)
        # Scalar/vector views resolve to the `AbstractField` sampler —
        # rewrap as the size field upstream's PostView operator is.
        result isa fields_mod.PostViewField &&
            (result=_GeoPostViewSizeField(result))
        return result
    elseif lowercase(spec.kind)=="automaticmeshsizefield"
        return _automatic_context_field(
            m,spec,config,entities,"execute_geo")
    elseif lowercase(spec.kind)=="boundarylayer"
        # Scalar BoundaryLayer over the discretized wall targets: upstream
        # feeds the boundary edges' mesh vertices (and endpoint chords) to
        # `FieldManager::BoundaryLayer`, then `hwall`/`ratio`/`thickness`/
        # `hfar` grow the prescribed size away from the wall. `_geo_field_
        # entity_mesh` supplies that same discretization per (dim, tag).
        fields_mod=getfield(parentmodule(@__MODULE__),:SizeField)
        points=NTuple{3,Float64}[]
        segments=Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
        triangles=Tuple{NTuple{3,Float64},NTuple{3,Float64},NTuple{3,Float64}}[]
        pt(em,i)=(em.coords[1,i],em.coords[2,i],em.coords[3,i])
        for (dim,tags) in ((0,config.points),(1,config.curves))
            for tag in tags
                mesh=entities(dim,string(tag))
                mesh===nothing && throw(ArgumentError(
                    "build_geo_size_field: Field[$(spec.tag)] references " *
                    "unknown entity ($dim,$tag)"))
                meshes=mesh isa AbstractVector ? mesh : (mesh,)
                for em in meshes
                    for i in axes(em.coords,2)
                        push!(points,pt(em,i))
                    end
                    for i in axes(em.segs,2)
                        a,b=Int(em.segs[1,i]),Int(em.segs[2,i])
                        push!(segments,(pt(em,a),pt(em,b)))
                    end
                    for i in axes(em.tris,2)
                        a,b,c=Int.(em.tris[:,i])
                        push!(triangles,(pt(em,a),pt(em,b),pt(em,c)))
                    end
                end
            end
        end
        isempty(points) && isempty(segments) && isempty(triangles) &&
            throw(ArgumentError(
                "build_geo_size_field: Field[$(spec.tag)] resolves no " *
                "boundary-layer distance targets"))
        distance=fields_mod.DistanceField(;points,segments,triangles)
        return fields_mod.BoundaryLayerField(distance;
            hwall=config.size,ratio=config.ratio,
            thickness=config.thickness,hfar=config.size_far)
    end
    throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)] kind $(spec.kind) " *
        "requires model context that is not available to .geo meshing"))
end

# The exec-time background field — `BGM_MeshSize` consumes the raw field
# (`l3`) and applies the global lcMin/lcMax/lcFactor itself, so the builder
# must not wrap it in `_gmsh_background`. `FieldManager::getBackgroundField`
# returning an undeclared id only fails when the field is queried — mirror
# that laziness by erroring at build time (the first use is the query).
function _geo_mesh_background_field(m::GeoModel,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    tag=context.background_field
    tag<=0 && return nothing
    fields=context.fields
    (fields===nothing || !haskey(fields,tag)) && throw(ArgumentError(
        "$caller: background field $tag is not defined"))
    opt(name)=_geo_option_number(context,"Mesh",0,name)
    geometry_tolerance=something(
        _geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)
    params=GeoParams(
        something(opt("MeshSizeMin"),0.0),
        something(opt("MeshSizeMax"),1e22),
        m.meshing.lc_factor,
        1,Dict{Tuple{Int,Int},String}(),fields,tag,
        copy(context.boundary_layer_fields),geometry_tolerance)
    bbox=try
        _model_mesh_bbox(m,caller)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    return _build_geo_field(
        params,(d,name)->_geo_field_entity_mesh(m,d,name,caller),tag;
        model_bbox=bbox,geometry_tolerance=geometry_tolerance,
        context_fields=(spec,config,entities,params)->
            _geo_mesh_context_fields(m,context,spec,config,entities,params))
end

# Mesh-time `_ModelMesh1DOptions` — upstream reads global `CTX` state inside
# `meshGEdge`/`BGM_MeshSize`; this bundle captures the equivalent `.geo`
# option values, the padded model-bounds `lc` (`SetBoundingBox` runs at
# `synchronize` before every `Mesh` statement), the selected background
# field, the size callback, and the `Msg::Warning`/`Msg::Error` sinks.
function _geo_mesh_1d_options(m::GeoModel,context::_GeoNumericContext,
                              caller::AbstractString)
    opt(name)=_geo_option_number(context,"Mesh",0,name)
    iopt(name,default)=_geo_signed_gmsh_int_value(
        something(opt(name),default),"$caller: Mesh.$name")
    geometry_tolerance=something(
        _geo_option_number(context,"Geometry",0,"Tolerance"),1e-8)
    return _ModelMesh1DOptions(
        context.native_context_lc,
        something(opt("MeshSizeMin"),0.0),
        something(opt("MeshSizeMax"),1e22),
        m.meshing.lc_factor,
        !iszero(something(opt("MeshSizeFromPoints"),1.0)),
        something(opt("MeshSizeFromCurvature"),0.0),
        something(opt("LcIntegrationPrecision"),1e-9),
        iopt("MinLineNodes",2.0),iopt("MinCircleNodes",7.0),
        iopt("MinCurveNodes",3.0),
        something(opt("ToleranceEdgeLength"),0.0),
        geometry_tolerance,
        !iszero(something(opt("MeshOnlyEmpty"),0.0)),
        !iszero(something(opt("MeshOnlyVisible"),0.0)),
        iopt("MaxRetries",10.0),
        _geo_mesh_background_field(m,context,caller),
        m.meshing.size_callback,
        msg->_geo_yywarn!(context,msg),
        msg->_geo_msg_error!(context,msg),
        Ref(Model._model_native_random_seed(Int(something(opt("RandomSeed"),1.0)))),
        something(opt("RandomFactor3D"),1e-12),
        Int(something(opt("Optimize"),1.0)),
        something(opt("OptimizeThreshold"),0.3),
        Int(something(opt("OptimizeNetgen"),0.0)),
        m.meshing.aniso_max)
end

# `Mesh n;` — mesh every live entity up to `dim`, merging entity parts into a
# single Mesh plus a node→(dim,tag) ownership map (lowest dimension wins, like
# Gmsh's per-entity vertex storage). `Mesh 1` re-grades every curve
# (`deMeshGEdge` upstream); `Mesh 2`/`Mesh 3` reuse the stored
# `curve_params` discretizations, re-running the 1-D pass only when a curve
# is missing parameters (upstream's `old < 1` arm).
include("GeoDim01Completion.jl")

function _geo_mesh_model(m::GeoModel,dim::Int,context::_GeoNumericContext)
    context.native_element_lookup_receipt=nothing
    context.native_mesh_begin_callback===nothing || context.native_mesh_begin_callback()
    caller="execute_geo: Mesh"
    if dim>=2 && !isempty(context.boundary_layer_fields)
        # `BoundaryLayer Field =` selections feed the wall-extrusion pass
        # upstream runs inside `meshGFace`/`meshGRegion`; the Tessella
        # surface/volume meshers have no boundary-layer consumer, so the
        # selection must reject loudly instead of silently meshing on.
        fields=context.fields
        for tag in context.boundary_layer_fields
            spec=fields===nothing ? nothing : get(fields,tag,nothing)
            spec===nothing && throw(ArgumentError(
                "$caller $dim: BoundaryLayer Field[$tag] is not declared"))
            lowercase(spec.kind)=="boundarylayer" || throw(ArgumentError(
                "$caller $dim: BoundaryLayer Field[$tag] has kind " *
                "$(spec.kind), not BoundaryLayer"))
        end
        throw(ArgumentError(
            "$caller $dim: `BoundaryLayer Field =` selects boundary-layer " *
            "mesh generation, which is not implemented"))
    end
    only_empty=!iszero(something(_geo_option_number(context,"Mesh",0,"MeshOnlyEmpty"),0.0))
    actual=(context.mesh isa MixedMesh && context.mesh.entity_data!==nothing) ||
        (m.meshing.mesh_authority!==nothing && m.meshing.mesh_authority.mesh!==nothing &&
         m.meshing.mesh_authority.carrier===context.mesh)
    high_order=dim>=1 && _geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"ElementOrder"),1.0),caller)>1
    records=any(!isempty(record.node_tags) || !isempty(record.element_tags) for
        (_,record) in Model._discrete_mesh_records_model(m))
    if dim>=2 || actual || high_order || records || (dim>=1 && only_empty) || Model._model_mesh_has_native_sources(m) || (dim>=2 &&
            context.mesh!==nothing && m.meshing.mesh_authority!==nothing)
        return _geo_mesh_actual_higher(m,dim,context,caller)
    end
    original=m
    # First-order Curve emission and its identity completion are one detached
    # operation. Geometry rows retain their existing emitter/merge order.
    dim==1 && (m=Model._model_mesh_transaction_copy(m))
    # For an ordinary native transfinite first-order producer, record its
    # actual Point/Curve/Surface scheduling before the merged carrier exists.
    # Extruded and unstructured factory traces have their separate planner.
    bind_completed=dim==1 || (dim==2 && context.mesh===nothing && !isempty(m.surfaces) &&
        something(_geo_option_number(context,"Mesh",0,"ElementOrder"),1.0)==1)
    if bind_completed && dim==2
        for tag in keys(m.surfaces)
            if !haskey(m.meshing.transfinite_surfaces,tag) || _extrude_gate(m,2,tag)!==nothing
                bind_completed=false
                break
            end
        end
    end
    producer_order=Tuple{Int,Int}[]
    bind_completed && append!(producer_order,((0,tag) for tag in sort!(collect(keys(m.points)))))
    # QuadTri certificates can reject after source grading. Stage the complete
    # operation, including 1-D grading and homology, so a failed request does
    # not publish model discretizations or mesh parts.
    staged=false
    if dim>=2
        for tag in keys(m.volumes)
            params=_extrude_gate(m,3,tag)
            if _extrude_is_nonew(params) ||
                    (dim==3 && _extrude_is_quadtri(params))
                staged=true
                break
            end
        end
    end
    staged && (m=deepcopy(m))
    _geo_sync_meshing_options!(m,context)
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[]
    if dim>=1
        options=_geo_mesh_1d_options(m,context,caller)
        if dim==1
            append!(parts,_model_mesh_dim01!(m,options,caller))
        else
            # Grade first (`meshDone1D` upstream); the parts are emitted only
            # after surface meshing so boundary splits written back to
            # `curve_params` by `meshGFace`-equivalent writeback appear in the
            # line elements.
            _model_mesh_dim01_grade!(m,options,caller)
        end
        bind_completed && append!(producer_order,((1,tag) for tag in sort!(collect(keys(m.curve_params)))))
    end
    # NoNew finalizes its region and boundary together. Keep this completed
    # plan within the current Mesh operation, including Mesh 2 queries of its
    # cap/laterals, instead of retaining a geometry-sensitive session cache.
    extrude_scope=dim>=2 ?
        _extrude_nonew_scope(m,caller;_working_model=staged) : nothing
    if dim>=2
        # OCC-curved primitive faces (Cylinder/Sphere/Cone/Torus side walls)
        # cannot standalone-mesh; when the face bounds a volume its facets
        # ride the volume's boundary recovery instead of emitting tri
        # elements (upstream emits the face's own tris — a documented gap).
        # A standalone unmeshable surface still fails through
        # `mesh_model_surface`'s own diagnostic.
        boundary_surfaces=Set{Int}()
        if dim==3
            for tag in keys(m.volumes)
                # Only the declared boundary list is needed for the OCC-curved
                # side-wall skip; the strict shell certification belongs to the
                # volume-mesh path, which admits valid unmerged sweep topology.
                for signed in _model_volume_boundary_surfaces(m,tag,caller;
                        certify_shells=false)
                    push!(boundary_surfaces,abs(signed))
                end
            end
        end
        for tag in sort!(collect(keys(m.surfaces)))
            if dim==3 && tag in boundary_surfaces &&
                    _surface_type(m,tag) in (:cylinder,:sphere,:cone,:torus)
                continue
            end
            push!(parts,(2,tag,mesh_model_surface(
                m,tag;_extrude_scope=extrude_scope,_mesh1d_options=options)))
            bind_completed && push!(producer_order,(2,tag))
        end
    end
    if dim==3
        # `meshGRegionExtruded` sweeps every `Layers`-marked volume ahead of
        # the unstructured pass, then `SubdivideExtrudedMesh` picks the
        # shared prism-face diagonals globally and remeshes the affected
        # lateral surfaces against that edge set.
        pass=_extrude_volume_pass(m,caller;scope=extrude_scope)
        for tag in sort!(collect(keys(m.volumes)))
            if pass!==nothing && haskey(pass.parts,tag)
                push!(parts,(3,tag,pass.parts[tag]))
            else
                push!(parts,(3,tag,mesh_model_volume(
                    m,tag;_extrude_scope=extrude_scope,_mesh1d_options=options)))
            end
        end
        if pass!==nothing && !isempty(pass.remesh)
            for (i,(pdim,ptag,_)) in enumerate(parts)
                pdim==2 && ptag in pass.remesh || continue
                parts[i]=(2,ptag,mesh_model_surface(
                    m,ptag;_extrude_edges=pass.edges,_extrude_scope=extrude_scope,
                    _mesh1d_options=options))
            end
        end
    end
    if dim>=2
        point_parts=_model_point_mesh_parts(m,caller)
        # NoNew retains CAD control points as well as the completed sweep.
        # Give those point parts a stable order within this scoped product.
        extrude_scope!==nothing && sort!(point_parts;by=part->part[2])
        append!(parts,point_parts)
        append!(parts,_model_curve_mesh_parts(m,caller))
    end
    isempty(parts) && throw(ArgumentError(
        "$caller $dim: no entities to mesh"))
    node_entities=_geo_mesh_part_node_entities(m,parts,caller)
    merged,owner,rows=_geo_merge_entity_meshes(parts;node_entities,return_remaps=true)
    if dim==1 && bind_completed
        produced=_geo_dim01_producer_trace(m,merged,parts,rows,caller)
        births=_geo_dim01_native_curve_births!(m,merged,parts,rows,produced,options,caller)
        plan=Model._model_mesh_generate_higher(m,merged,parts,rows,produced.owners,
            dim,nothing,caller;operation=:bind_completed,
            producer_order=produced.order,producer_parameters=births.parameters,
            producer_coordinates=births.coordinates,
            producer_geometry_parameters=births.geometry_parameters,
            renumber=!iszero(something(_geo_option_number(context,"Mesh",0,"Renumber"),1.0)),
            save_all=!iszero(something(_geo_option_number(context,"Mesh",0,"SaveAll"),0.0)),
            initial_max_node_tag=context.mesh_node_max,
            initial_max_element_tag=context.mesh_element_max)
        plan=merge(plan,(geometry_parts=parts,geometry_part_rows=rows))
        return _geo_publish_actual_plan!(original,context,plan,caller;
            homology=true,carrier=merged)
    end
    parameters=_geo_mesh_stored_parameters(m,merged,owner,context)
    # Retain the per-entity decomposition so classified projections
    # (`model_to_mixed`) can recover each entity's own mesh from the merged
    # product — upstream serializes vertex/edge/face/region cells to their
    # owning entity, and the part mesh is that view.
    _geo_run_homology!(m,merged,parts,caller;node_owner=owner,
        part_node_entities=node_entities,error_callback=_geo_homology_error_callback(context))
    completed_plan=bind_completed ? Model._model_mesh_generate_higher(m,merged,parts,rows,
        owner,dim,nothing,caller;operation=:bind_completed,producer_order,
        initial_max_node_tag=context.mesh_node_max,
        initial_max_element_tag=context.mesh_element_max) : nothing
    if completed_plan!==nothing
        parts,rows=_geo_actual_product_parts(completed_plan.model.meshing.mesh_authority,caller)
        owner=_geo_actual_carrier_owners(completed_plan.model.meshing.mesh_authority,merged)
        for name in fieldnames(GeoModel)
            setfield!(m,name,getfield(completed_plan.model,name))
        end
    end
    if staged
        for name in fieldnames(GeoModel)
            setfield!(original,name,getfield(m,name))
        end
    end
    context.mesh_parts=parts
    context.mesh_part_rows=rows
    context.mesh_geometry_parts=nothing
    if completed_plan===nothing
        context.mesh_node_max=max(context.mesh_node_max,UInt64(nnodes(merged)))
        context.mesh_element_max=max(context.mesh_element_max,UInt64(sum(
            _geo_mesh_dim_cells(merged,dimension) for dimension in 1:3)))
        original.meshing.mesh_authority=Model._model_mesh_completed(original,dim)
    else
        context.mesh_node_max=completed_plan.max_node_tag
        context.mesh_element_max=completed_plan.max_element_tag
    end
    context.mesh_node_parametric=parameters
    return merged,owner
end

function _geo_mesh_working_model(m,context)
    working=Model._model_mesh_transaction_copy(m)
    previous=m.meshing.mesh_authority
    if previous!==nothing && previous.mesh!==nothing && previous.carrier===context.mesh
        detached=working.meshing.mesh_authority
        # The geometry carrier is an input owned by this context. Retain its
        # identity while copying the actual products that the plan can mutate.
        # The carrier swap keeps `detached.mesh` and its tag space, so the
        # retained Surface constructor inputs stay bound to the same products.
        working.meshing.mesh_authority=Model._ModelMeshAuthority(detached.mesh,
            detached.products,detached.done,detached.private_rank,detached.native_inputs,
            context.mesh,detached.carrier_rows,detached.producer_rows,detached.carrier_cells,
            detached.surface_constructor_inputs)
    end
    _geo_sync_meshing_options!(working,context)
    return working
end

function _geo_mesh_actual_higher(m,dim,context,caller)
    working=_geo_mesh_working_model(m,context)
    order=_geo_signed_gmsh_int_value(something(
        _geo_option_number(context,"Mesh",0,"ElementOrder"),1.0),caller)
    settings=(;algorithm_2d=working.meshing.algorithm_2d,
        algorithm_3d=working.meshing.algorithm_3d,
        flexible_transfinite=working.meshing.flexible_transfinite,
        lc_factor=working.meshing.lc_factor,recombine_all=working.meshing.recombine_all,
        recombine_algo=working.meshing.recombine_algo,
        recombine_topology_passes=working.meshing.recombine_topology_passes,
        recombine_node_repositioning=working.meshing.recombine_node_repositioning,
        recombine_minimum_quality=working.meshing.recombine_minimum_quality,
        native_smoothing=working.meshing.native_smoothing,
        optimize_pyramids=something(_geo_option_number(context,"Mesh",0,"OptimizePyramids"),0.0),
        second_order_linear=working.meshing.second_order_linear,
        transfinite_tri=working.meshing.transfinite_tri,
        lc_extend_from_boundary=working.meshing.lc_extend_from_boundary,
        aniso_max=working.meshing.aniso_max,order)
    only_empty=!iszero(something(_geo_option_number(context,"Mesh",0,"MeshOnlyEmpty"),0.0))
    previous=m.meshing.mesh_authority
    preserve_geometry=dim>=2 && !only_empty &&
        (context.mesh===nothing || context.mesh isa Mesh ||
         (previous!==nothing && previous.carrier isa Mesh) ||
         # A retained non-simplex geometry carrier is still the geometry view.
         (previous!==nothing && previous.mesh!==nothing &&
          previous.carrier!==nothing && previous.carrier!==previous.mesh &&
          context.mesh===previous.carrier))
    plan=Model._model_mesh_generate_higher(working,context.mesh,context.mesh_parts,
        context.mesh_part_rows,context.mesh_node_owner,dim,
        _geo_mesh_1d_options(working,context,caller),caller;
        generation_meshing=settings,element_order=order,
        linear=!iszero(something(_geo_option_number(context,"Mesh",0,"SecondOrderLinear"),0.0)),
        incomplete=!iszero(something(_geo_option_number(context,"Mesh",0,"SecondOrderIncomplete"),0.0)),
        renumber=!iszero(something(_geo_option_number(context,"Mesh",0,"Renumber"),1.0)),
        save_all=!iszero(something(_geo_option_number(context,"Mesh",0,"SaveAll"),0.0)),
        size_field=_geo_mesh_background_field(working,context,caller),
        initial_max_node_tag=context.mesh_node_max,
        initial_max_element_tag=context.mesh_element_max,preserve_geometry,geometry_layout=:all,
        geometry_merge=_geo_geometry_carrier_merge)
    if dim==1 && order==1
        plan=_geo_dim01_generation_carrier!(plan,previous,caller)
    end
    return _geo_publish_actual_plan!(m,context,plan,caller;
        homology=true,carrier=plan.model.meshing.mesh_authority.carrier)
end

function _geo_actual_product_parts(authority,caller)
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[];rows=Vector{Int32}[]
    for ((dimension,tag),product) in sort!(collect(authority.products);by=first)
        push!(parts,(dimension,Int(tag),Model._model_mesh_product(authority,dimension,Int(tag);caller)))
        push!(rows,copy(product.rows))
    end
    return parts,rows
end

function _geo_actual_carrier_owners(authority,carrier)
    rows=carrier===authority.mesh ? collect(Int32,1:nnodes(carrier)) : authority.carrier_rows
    length(rows)==nnodes(carrier) || throw(ArgumentError(
        "execute_geo: actual carrier rows do not align with its nodes"))
    data=authority.mesh.entity_data
    return Dict(row=>(data.node_entities[actual][1],Int(data.node_entities[actual][2]))
        for (row,actual) in enumerate(rows))
end

function _geo_publish_actual_plan!(m,context,plan,caller;homology=false,carrier=plan.mesh)
    if homology && !isempty(plan.model.meshing.homology_requests)
        plan=_geo_actual_homology_plan(plan,caller;
            error_callback=_geo_homology_error_callback(context))
        carrier=plan.model.meshing.mesh_authority.carrier
    end
    authority=plan.model.meshing.mesh_authority
    parts,rows=_geo_actual_product_parts(authority,caller)
    owner=_geo_actual_carrier_owners(authority,carrier)
    actual_rows=carrier===authority.mesh ? collect(Int32,1:nnodes(carrier)) : authority.carrier_rows
    parameters=Model._model_mesh_runtime_parameters(plan.model,authority.mesh)
    carrier_parameters=hasproperty(plan,:geometry_node_parametric) ?
        plan.geometry_node_parametric[actual_rows] : parameters[actual_rows]
    geometry_receipt=_geo_geometry_plan_receipt(authority,carrier,
        context.mesh_geometry_parts,plan,caller)
    for name in fieldnames(GeoModel)
        setfield!(m,name,getfield(plan.model,name))
    end
    context.mesh_parts=parts;context.mesh_part_rows=rows
    context.mesh_geometry_parts=geometry_receipt
    context.mesh_node_max=plan.max_node_tag;context.mesh_element_max=plan.max_element_tag
    context.mesh=carrier;context.mesh_node_owner=owner
    context.mesh_node_parametric=carrier_parameters
    return carrier,owner
end

# `PartitionMesh n;` — partition the current mesh into `n` parts, the
# session-API `mesh.partition(n)` equivalent (`GModel::partitionMesh`
# upstream). The session API is reached lazily because `GeoExec` is loaded
# before `API`: deterministic element assignment reuses
# `API._partition_elements` plus the upstream lower-cell inheritance rule
# (each lower-dimension cell follows an adjacent carrier at the nearest higher
# dimension). `Mesh.PartitionCreateTopology`/`Mesh.PartitionCreateGhostCells`
# route through the Tet4/lower-BRep production plan; when both are zero the
# assignment-only producer still materializes the per-(entity, partition)
# child entities, so the actual published mesh carries complete
# `MixedPartitionData` — `Save` then emits `$PartitionedEntities` and the
# partition element columns.
function _geo_exec_partition_mesh!(m::GeoModel,args::AbstractString,
                                   context::_GeoNumericContext,caller::AbstractString)
    caller="$caller: PartitionMesh"
    count=_geo_int_value(_geo_eval_numeric(String(strip(args)),context,caller),caller)
    count>=1 || throw(ArgumentError(
        "$caller: PartitionMesh requires a positive partition count"))
    _geo_sync_physical_view_if_changed!(m,context)
    context.mesh===nothing && throw(ArgumentError(
        "$caller: PartitionMesh requires a current mesh — run `Mesh n;` first"))
    _geo_check_part_rows(context)
    topology=!iszero(something(
        _geo_option_number(context,"Mesh",0,"PartitionCreateTopology"),1.0))
    ghosts=!iszero(something(
        _geo_option_number(context,"Mesh",0,"PartitionCreateGhostCells"),0.0))
    api=getfield(parentmodule(@__MODULE__),:API)
    source,actual,class=api._higher_frontend_input(m,context.mesh,
        context.mesh_parts,context.mesh_part_rows,context.mesh_node_owner,caller;
        initial_max_node_tag=context.mesh_node_max,
        initial_max_element_tag=context.mesh_element_max)
    actual===nothing && throw(ArgumentError(
        "$caller: PartitionMesh requires a current mesh — run `Mesh n;` first"))
    _,_,total=api._mesh_element_offsets(actual)
    total>0 || throw(ArgumentError(
        "$caller: PartitionMesh requires a mesh with elements"))
    assignment=api._partition_lower_inherit(api._partition_elements(actual,count),actual)
    plan=topology||ghosts ?
        api._partition_tet_production_plan(source,actual,class,count,assignment,
            caller;topology,ghosts) :
        api._partition_assign_production_plan(source,actual,class,count,assignment,caller)
    _geo_publish_actual_plan!(m,context,plan,caller)
    return nothing
end

function _geo_refined_simplex_carrier!(plan,previous,caller;source_authority=nothing,geometry_receipt=nothing)
    trace=Refine._refine_uniform_traced(previous)
    tags=plan.refine_carrier_tags
    length(tags)==nnodes(previous) || throw(ArgumentError(
        "$caller: original carrier identities do not align with refinement input"))
    rows=Vector{Int32}(undef,nnodes(trace.mesh))
    old_rows=zeros(Int32,count(!iszero,trace.primary_map))
    for (old,compact) in enumerate(trace.primary_map)
        compact==0 && continue
        old_rows[compact]=Int32(old)
        rows[compact]=plan.refine_primary_rows[tags[old]]
    end
    for ((a,b),midpoint) in trace.midpoint_ids
        key=minmax(tags[old_rows[a]],tags[old_rows[b]])
        rows[midpoint]=plan.refine_edge_rows[key]
    end
    allunique(rows) && all(row->0<row<=nnodes(plan.mesh),rows) || throw(ArgumentError(
        "$caller: refined carrier identities do not form a valid injection"))
    refined=trace.mesh
    carrier=Mesh(plan.mesh.coords[:,rows];segs=refined.segs,tris=refined.tris,tets=refined.tets,
        seg_tag=refined.seg_tag,tri_tag=refined.tri_tag,tet_tag=refined.tet_tag)
    previous_authority=plan.model.meshing.mesh_authority
    cells=source_authority===nothing ? Vector{UInt64}[] :
        Model._model_mesh_refined_cell_tags(plan,source_authority;caller)
    authority=Model._model_mesh_authority(plan.mesh,previous_authority.done;
        caller,private_rank=previous_authority.private_rank,carrier,carrier_rows=rows,carrier_cells=cells)
    plan.model.meshing.mesh_authority=Model._ModelMeshAuthority(authority.mesh,
        authority.products,authority.done,authority.private_rank,previous_authority.native_inputs,
        carrier,rows,previous_authority.producer_rows,authority.carrier_cells)
    parts,part_rows=_geo_geometry_refined_parts(geometry_receipt,trace,carrier,caller)
    return carrier,parts,part_rows
end

function _geo_operate_actual_mesh!(m,context,operation,caller;element_order=1)
    working=_geo_mesh_working_model(m,context)
    plan=Model._model_mesh_generate_higher(working,context.mesh,context.mesh_parts,
        context.mesh_part_rows,context.mesh_node_owner,0,nothing,caller;
        operation,element_order,
        linear=!iszero(something(_geo_option_number(context,"Mesh",0,"SecondOrderLinear"),0.0)),
        incomplete=!iszero(something(_geo_option_number(context,"Mesh",0,"SecondOrderIncomplete"),0.0)),
        renumber=!iszero(something(_geo_option_number(context,"Mesh",0,"Renumber"),1.0)),
        save_all=!iszero(something(_geo_option_number(context,"Mesh",0,"SaveAll"),0.0)),
        initial_max_node_tag=context.mesh_node_max,
        initial_max_element_tag=context.mesh_element_max)
    carrier=if operation===:refine && context.mesh isa Mesh
        refined,parts,part_rows=_geo_refined_simplex_carrier!(plan,context.mesh,caller;
            source_authority=working.meshing.mesh_authority,
            geometry_receipt=context.mesh_geometry_parts)
        parts===nothing || (plan=merge(plan,(geometry_parts=parts,geometry_part_rows=part_rows)))
        refined
    elseif operation in (:geo_renumber_nodes,:geo_renumber_elements)
        plan.model.meshing.mesh_authority.carrier
    else
        plan.mesh
    end
    _geo_publish_actual_plan!(m,context,plan,caller;carrier)
    # RefineMesh and SetOrder destroy both native tag caches before their
    # conversion. A read-time import catalogue cannot survive the operation.
    operation in (:refine,:order) && (context.native_element_lookup_receipt=nothing)
    return nothing
end

include("GeoRecombineBipartite.jl")

function _geo_recombine_actual_plan!(m,context,caller)
    previous=m.meshing.mesh_authority
    original=context.mesh
    working=_geo_mesh_working_model(m,context)
    algorithm=working.meshing.recombine_algo in (1,3) ? :blossom : :greedy
    plan=Model._model_mesh_generate_higher(working,context.mesh,context.mesh_parts,
        context.mesh_part_rows,context.mesh_node_owner,0,nothing,caller;
        operation=:recombine,algorithm,
        initial_max_node_tag=context.mesh_node_max,
        initial_max_element_tag=context.mesh_element_max)
    if original isa MixedMesh && original.entity_data===nothing &&
            previous!==nothing && previous.carrier===original
        carrier=_geo_recombined_mixed_carrier!(plan,previous,original,caller)
        _geo_publish_actual_plan!(m,context,plan,caller;carrier)
    elseif plan.recombined || context.mesh isa MixedMesh && context.mesh.entity_data!==nothing
        _geo_publish_actual_plan!(m,context,plan,caller)
    else
        # A strict angle rejection keeps its original geometry carrier. The
        # first-pair temporary factories still advance the historical maximum.
        _geo_publish_actual_plan!(m,context,plan,caller;carrier=context.mesh)
    end
    return nothing
end

function _geo_recombined_mixed_carrier!(plan,previous,original,caller)
    actual=plan.mesh;data=actual.entity_data
    old_tags=previous.mesh.entity_data.external_node_tags[previous.carrier_rows]
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(data.external_node_tags))
    rows=Int32[get(positions,tag,Int32(0)) for tag in old_tags]
    length(rows)==nnodes(original) && all(>(0),rows) && allunique(rows) ||
        throw(ArgumentError("$caller: recombined geometry loses its bound vertex identities"))
    actual.coords[:,rows]==original.coords || throw(ArgumentError(
        "$caller: fixed-node recombination changes bound geometry coordinates"))
    carrier=original
    carrier_cells=Model._model_mesh_rebound_cell_tags(previous,plan.element_tag_map)
    if plan.recombined
        original.elementary_entities===nothing && throw(ArgumentError(
            "$caller: native Mixed geometry lacks its producer cell owners"))
        order=Tuple{Int,Int32}[]
        physical=Dict{Tuple{Int,Int32},Int32}()
        for (block,owners) in zip(original.blocks,original.elementary_entities)
            block isa ElementBlock || throw(ArgumentError(
                "$caller: native Mixed geometry contains an unsupported special cell"))
            for column in axes(block.nodes,2)
                owner=(msh_dimension(block.msh),owners[column])
                if !haskey(physical,owner)
                    push!(order,owner);physical[owner]=block.tags[column]
                elseif physical[owner]!=block.tags[column]
                    throw(ArgumentError(
                        "$caller: native geometry needs explicit parent binding for varying cell tags"))
                end
            end
        end
        inverse=zeros(Int32,nnodes(actual));inverse[rows]=Int32.(1:length(rows))
        blocks=ElementBlock[];elementary=Vector{Int32}[];carrier_cells=Vector{UInt64}[]
        for owner in order,(bi,(block,owners)) in enumerate(zip(actual.blocks,data.block_entities))
            msh_dimension(block.msh)==owner[1] || continue
            columns=Int[column for column in axes(block.nodes,2) if owners[column]==owner[2]]
            isempty(columns) && continue
            nodes=inverse[block.nodes[:,columns]]
            all(>(0),nodes) || throw(ArgumentError(
                "$caller: recombined cell references a vertex outside its geometry carrier"))
            push!(blocks,ElementBlock(block.msh,nodes,fill(physical[owner],length(columns))))
            push!(elementary,fill(owner[2],length(columns)))
            push!(carrier_cells,copy(data.external_element_tags[bi][columns]))
        end
        carrier=MixedMesh(copy(original.coords),blocks;
            physical_names=original.physical_names,elementary_entities=elementary,
            periodic_links=original.periodic_links,ancillary_sections=original.ancillary_sections,
            data_sections=original.data_sections,partition_data=original.partition_data)
    end
    authority=plan.model.meshing.mesh_authority
    # The geometry carrier swap keeps `authority.mesh` and its products, so the
    # retained Surface constructor inputs stay bound to the same products.
    plan.model.meshing.mesh_authority=Model._ModelMeshAuthority(authority.mesh,
        authority.products,authority.done,authority.private_rank,authority.native_inputs,
        carrier,rows,authority.producer_rows,carrier_cells,
        authority.surface_constructor_inputs)
    return carrier
end

function _geo_rebind_product_coordinates!(m,context)
    mesh=context.mesh
    previous=m.meshing.mesh_authority
    if (mesh isa Mesh || mesh isa MixedMesh && mesh.entity_data===nothing) &&
            previous!==nothing && previous.mesh!==nothing && previous.carrier===mesh
        source=previous.mesh;coords=copy(source.coords)
        coords[:,previous.carrier_rows]=mesh.coords
        mesh=MixedMesh(coords,source.blocks;physical_names=source.physical_names,
            entity_data=source.entity_data,elementary_entities=source.elementary_entities,
            periodic_links=source.periodic_links,ancillary_sections=source.ancillary_sections,
            data_sections=source.data_sections,partition_data=source.partition_data)
    end
    for ((_,_,part),rows) in zip(context.mesh_parts,context.mesh_part_rows)
        length(rows)==nnodes(part) || throw(ArgumentError(
            "execute_geo: stored entity rows do not align with their actual product"))
        part.coords[:,:]=mesh.coords[:,rows]
    end
    if mesh isa MixedMesh && mesh.entity_data!==nothing
        data=mesh.entity_data
        positions=Dict(tag=>row for (row,tag) in enumerate(data.external_node_tags))
        for (_,record) in Model._discrete_mesh_records_model(m),column in eachindex(record.node_tags)
            row=get(positions,UInt64(record.node_tags[column]),0)
            row==0 || (record.node_coords[:,column]=mesh.coords[:,row])
        end
        ranks=previous!==nothing && length(previous.private_rank)==nnodes(mesh) ?
              previous.private_rank : collect(UInt64,1:nnodes(mesh))
        done=previous===nothing ? Set{Tuple{Int,Int32}}() : previous.done
        authority=Model._model_mesh_authority(mesh,done;private_rank=ranks,
            carrier=context.mesh,carrier_rows=previous!==nothing && previous.carrier===context.mesh ?
                previous.carrier_rows : collect(Int32,1:nnodes(mesh)),
            carrier_cells=previous!==nothing && previous.carrier===context.mesh ?
                previous.carrier_cells : nothing)
        if previous!==nothing && previous.mesh!==nothing
            authority=Model._model_mesh_with_producer_rows(authority,previous.producer_rows)
        end
        m.meshing.mesh_authority=authority
    end
    return nothing
end

function _geo_forget_mesh_products!(m,context;reset_tags=false)
    empty!(context.mesh_parts);empty!(context.mesh_part_rows)
    context.mesh_geometry_parts=nothing
    empty!(context.mesh_node_parametric)
    m.meshing.mesh_authority=nothing
    if reset_tags
        context.mesh_node_max=UInt64(0);context.mesh_element_max=UInt64(0)
    end
    return nothing
end

# Run the queued `Homology`/`Cohomology`/`Betti` requests — upstream calls
# `GModel::computeHomology` at the end of `GModel::mesh`, once every entity's
# mesh vertices exist. The generated mesh is preserved and each surviving
# generator is stored as a discrete entity under a named physical group;
# pinned Gmsh 4.15.2's `.msh` output drops the mesh and serializes only the
# chains, which Tessella deliberately does not mirror.
function _geo_run_homology!(m::GeoModel,mesh,
                          parts::Vector{Tuple{Int,Int,Union{Mesh,MixedMesh}}},
                          caller::AbstractString;node_owner=nothing,part_node_entities=nothing,
                          error_callback=nothing)
    isempty(m.meshing.homology_requests) && return nothing
    coordinates=Dict{Int32,NTuple{3,Float64}}()
    cells=_geo_homology_cells(m,mesh,parts,caller;node_owner,part_node_entities,
                              homology_coordinates=coordinates)
    compute_homology!(m,cells;requests=m.meshing.homology_requests,coordinates,error_callback)
    return nothing
end

# Entity→element map for `compute_homology!`, in one node-index space keyed by
# entity ownership and exact coordinates: the generated mesh's nodes first, then any
# extra nodes the curve/discrete records need (homology indices only identify
# cells — they need not exist in `mesh`).
function _geo_homology_cells(m::GeoModel,mesh,
                             parts::Vector{Tuple{Int,Int,Union{Mesh,MixedMesh}}},
                             caller::AbstractString;node_owner=nothing,part_node_entities=nothing,
                             homology_coordinates=nothing)
    lookup=Dict{Tuple{Int,Int,NTuple{3,Int}},Int32}()
    legacy_lookup=Dict{NTuple{3,Int},Int32}()
    mcoords=mesh.coords
    for n in 1:nnodes(mesh)
        entity=node_owner===nothing ? (-1,0) : get(node_owner,n,(3,0))
        x,y,z=mcoords[1,n],mcoords[2,n],mcoords[3,n]
        homology_coordinates===nothing || (homology_coordinates[Int32(n)]=(x,y,z))
        coordinate=_geo_mesh_coordinate_key(x,y,z)
        lookup[(entity...,coordinate)]=Int32(n)
        legacy_lookup[coordinate]=Int32(n)
    end
    # One-element counter — mutating `next[1]` inside `index_for` leaves the
    # captured binding itself assigned once, keeping it unboxed.
    next=Int32[nnodes(mesh)+1]
    function index_for(x,y,z,entity)
        coordinate=_geo_mesh_coordinate_key(x,y,z)
        key=(entity...,coordinate)
        i=entity==(-1,0) ? get(legacy_lookup,coordinate,Int32(0)) : get(lookup,key,Int32(0))
        i!=0 && return i
        i=next[1];lookup[key]=i;legacy_lookup[coordinate]=i;next[1]=i+Int32(1)
        homology_coordinates===nothing || (homology_coordinates[i]=(Float64(x),Float64(y),Float64(z)))
        return i
    end
    cells=Dict{Tuple{Int,Int},Vector{Tuple{Int32,Vector{Int32}}}}()
    part_nodes=Dict{Tuple{Int,Int},Dict{NTuple{3,Int},Int32}}()
    emitted=Set{Tuple{Int,Int}}()
    point_sources=Set{Tuple{Int,Int}}()
    function entry(key)
        get!(cells,key,Tuple{Int32,Vector{Int32}}[])
    end
    # Meshed entities own their top-dimension elements; boundary elements
    # belong to the boundary entities (curves below, discrete records after).
    for (part_index,(dim,tag,part)) in enumerate(parts)
        pcoords=part.coords
        remap=Int32[index_for(pcoords[1,n],pcoords[2,n],pcoords[3,n],
                    part_node_entities===nothing ? (-1,0) : part_node_entities[part_index][n])
                  for n in 1:nnodes(part)]
        elements=entry((dim,tag))
        part_nodes[(dim,tag)]=Dict(_geo_mesh_coordinate_key(
            pcoords[1,n],pcoords[2,n],pcoords[3,n])=>remap[n]
            for n in 1:nnodes(part))
        if part isa Mesh
            if dim==0
                isempty(remap) || push!(point_sources,(0,tag))
                record=get(m.discrete,(0,tag),get(m.meshing.attached,(0,tag),nothing))
                # Mesh0D keeps explicitly stored MPoint elements. When none
                # exist it creates one on the last stored vertex; the other
                # classified vertices remain mesh nodes without point cells.
                if !isempty(remap) &&
                   (record===nothing || !(Int32(15) in record.element_types))
                    push!(elements,(Int32(15),Int32[remap[end]]))
                end
            elseif dim==1
                for t in 1:nsegs(part)
                    push!(elements,(Int32(1),
                        Int32[remap[part.segs[k,t]] for k in 1:2]))
                end
            elseif dim==2
                for t in 1:ntris(part)
                    push!(elements,(Int32(2),
                        Int32[remap[part.tris[k,t]] for k in 1:3]))
                end
            elseif dim==3
                for t in 1:ntets(part)
                    push!(elements,(Int32(4),
                        Int32[remap[part.tets[k,t]] for k in 1:4]))
                end
            end
        else
            for block in part.blocks
                block isa ElementBlock || continue
                msh_dimension(block.msh)==dim || continue
                width=size(block.nodes,1)
                for cell in axes(block.nodes,2)
                    push!(elements,(Int32(block.msh),
                        Int32[remap[block.nodes[k,cell]] for k in 1:width]))
                end
            end
        end
        isempty(elements) || push!(emitted,(dim,tag))
    end
    # Native curve chains above use the emitted segment parts, preserving
    # topological identities at coincident points and the actual 1-D mesh.
    for (tag,point) in m.points
        ((0,tag) in emitted || (0,tag) in point_sources) && continue
        record=get(m.discrete,(0,tag),get(m.meshing.attached,(0,tag),nothing))
        record!==nothing && Int32(15) in record.element_types && continue
        index=index_for(point...,node_owner===nothing ? (-1,0) : (0,tag))
        push!(entry((0,tag)),(Int32(15),Int32[index]))
    end
    # Discrete entities (merged meshes, classified topology, meshing
    # attachments) already store their elements against record node tags.
    # Connectivity addresses the records' shared tag space — a node may live
    # on a different entity's record (shared boundary vertices, and the
    # child-entity migration `create_topology!` performs) — so tags resolve
    # against a model-wide map, not the owning record.
    global_nodes=Dict{Int32,Int32}()
    records=sort!(_discrete_mesh_records_model(m);by=first)
    for (key,record) in records
        rcoords=record.node_coords
        actual=get(part_nodes,key,nothing)
        for (i,t) in enumerate(record.node_tags)
            get!(global_nodes,t) do
                x,y,z=rcoords[1,i],rcoords[2,i],rcoords[3,i]
                coordinate=_geo_mesh_coordinate_key(x,y,z)
                # Resolve through the current part, including its boundary
                # owners. Distinct coincident record tags keep their entity
                # identity; shared tags use this model-wide map.
                mapped=actual===nothing ? Int32(0) :
                    get(actual,coordinate,Int32(0))
                mapped==0 ? index_for(x,y,z,key) : mapped
            end
        end
    end
    for (key,record) in records
        # Generated parts replace stored mesh connectivity, which may have
        # changed during remeshing. Re-adding records would duplicate cells
        # or attach a stale chain to another coincident entity.
        key in emitted && continue
        elements=entry(key)
        for (msh_type,nodes) in zip(record.element_types,
                                    record.element_nodes)
            mapped=Vector{Int32}(undef,length(nodes))
            for (k,t) in enumerate(nodes)
                i=get(global_nodes,t,Int32(0))
                i==0 && throw(ArgumentError(
                    "$caller: discrete entity $key references unknown " *
                    "node tag $t"))
                mapped[k]=i
            end
            push!(elements,(Int32(msh_type),mapped))
        end
    end
    return cells
end

@inline _geo_mesh_coordinate_key(x,y,z)=(
    reinterpret(Int,x==0 ? 0.0 : x),reinterpret(Int,y==0 ? 0.0 : y),
    reinterpret(Int,z==0 ? 0.0 : z))

function _geo_mesh_children(m,dim,tag,caller)
    children=Tuple{Int,Int}[]
    if dim==1 && haskey(m.curves,tag)
        a,b=m.curves[tag];push!(children,(0,a),(0,b))
    elseif dim==2 && haskey(m.surfaces,tag)
        append!(children,((1,c) for c in _model_projection_surface_curves(m,tag)))
    elseif dim==3
        append!(children,((2,abs(s)) for s in _model_volume_boundary_surfaces(m,tag,caller)))
    end
    append!(children,get(m.embeds,(dim,tag),Tuple{Int,Int}[]))
    return children
end

# Classify from actual lower-dimensional mesh parts. This avoids geometric
# projection, preserves independent coincident entities, and shares boundary
# identities through the model's topology. Distinct interior/tab evaluations
# remain volume nodes unless they exactly share a lower-dimensional node.
function _geo_mesh_point_incidence_coordinate(m,point)
    haskey(m.points,point) && return m.points[point]
    record=get(m.discrete,(0,point),get(m.meshing.attached,(0,point),nothing))
    (record===nothing || isempty(record.node_tags)) && return nothing
    return (record.node_coords[1,1],record.node_coords[2,1],record.node_coords[3,1])
end

function _geo_mesh_constraint_point_owners(m,parts,caller)
    isempty(m.embeds) && return nothing
    aliases=Dict{Int,Dict{NTuple{3,Int},Tuple{Int,Int}}}()
    entities=Dict((dim,tag)=>mesh for (dim,tag,mesh) in parts)
    for ((dim,tag),constraints) in m.embeds
        dim in (2,3) || continue
        mesh=get(entities,(dim,tag),nothing)
        mesh===nothing && continue
        points=Int[p for (d,p) in constraints if d==0]
        curves=Int[c for (d,c) in constraints if d==1 && haskey(m.curve_params,c)]
        (isempty(points) || isempty(curves)) && continue
        # Only coordinates referenced by actual carrier cells can transmit
        # recovered constraint incidence. Unused interpolation/tab vertices
        # and unrelated coincident model points keep their own identities.
        candidates=Dict{NTuple{3,Int},Vector{Int}}()
        for point in points
            coordinate=_geo_mesh_point_incidence_coordinate(m,point)
            coordinate===nothing && continue
            key=_geo_mesh_coordinate_key(coordinate...)
            matches=get!(candidates,key,Int[])
            point in matches || push!(matches,point)
        end
        used=Set{NTuple{3,Int}}()
        cells=mesh isa Mesh ? (dim==2 ? (mesh.tris,) : (mesh.tets,)) :
            (b.nodes for b in mesh.blocks if b isa ElementBlock && msh_dimension(b.msh)==dim)
        for block in cells,node in block
            key=_geo_mesh_coordinate_key(mesh.coords[1,node],
                mesh.coords[2,node],mesh.coords[3,node])
            haskey(candidates,key) && push!(used,key)
        end
        for curve in curves
            curve_mesh=get(entities,(1,curve),nothing)
            curve_mesh===nothing && continue
            source=_model_mesh_part_node_entities(m,1,curve,curve_mesh,caller)
            for node in eachindex(source)
                source[node]==(1,curve) || continue
                key=_geo_mesh_coordinate_key(curve_mesh.coords[1,node],
                    curve_mesh.coords[2,node],curve_mesh.coords[3,node])
                key in used || continue
                matches=get(candidates,key,nothing)
                matches===nothing && continue
                length(matches)==1 || throw(ArgumentError(
                    "$caller: embedded Points[$(matches[1])] and [$(matches[2])] " *
                    "coincide on Curve[$curve]; preserving their separate " *
                    "constraint identities is not implemented"))
                point=only(matches)
                owners=get!(aliases,curve,Dict{NTuple{3,Int},Tuple{Int,Int}}())
                previous=get(owners,key,(0,point))
                previous==(0,point) || throw(ArgumentError(
                    "$caller: embedded Points[$(previous[2])] and [$point] " *
                    "coincide on Curve[$curve]; preserving their separate " *
                    "constraint identities is not implemented"))
                owners[key]=(0,point)
            end
        end
    end
    return aliases
end

function _geo_mesh_part_node_entities(m,parts,caller)
    result=Vector{Vector{Tuple{Int,Int}}}(undef,length(parts))
    maps=Dict{Tuple{Int,Int},Dict{NTuple{3,Int},Tuple{Int,Int}}}()
    sizehint!(maps,length(parts))
    constraint_points=_geo_mesh_constraint_point_owners(m,parts,caller)
    for dim in 0:3, (index,(pdim,tag,mesh)) in enumerate(parts)
        pdim==dim || continue
        owners=fill((dim,tag),nnodes(mesh))
        if dim<=1
            source=_model_mesh_part_node_entities(m,dim,tag,mesh,caller)
            for n in eachindex(owners);owners[n]=(source[n][1],Int(source[n][2]));end
            aliases=dim==1 && constraint_points!==nothing ? get(constraint_points,tag,nothing) : nothing
            if aliases!==nothing
                for n in eachindex(owners)
                    owners[n]==(1,tag) || continue
                    key=_geo_mesh_coordinate_key(mesh.coords[1,n],mesh.coords[2,n],mesh.coords[3,n])
                    owners[n]=get(aliases,key,owners[n])
                end
            end
        else
            pending=_geo_mesh_children(m,dim,tag,caller)
            seen=Set{Tuple{Int,Int}}()
            entries=Dict{NTuple{3,Int},Tuple{Int,Int}}[]
            boundary_capacity=0
            while !isempty(pending)
                child=pop!(pending)
                child in seen && continue
                push!(seen,child)
                child[1]<dim || continue
                entry=get(maps,child,nothing)
                if entry!==nothing
                    push!(entries,entry)
                    boundary_capacity=Base.checked_add(boundary_capacity,length(entry))
                end
                child[1]>0 && append!(pending,_geo_mesh_children(m,child...,caller))
            end
            # Interior-heavy parts need capacity for their actual child maps,
            # rather than for every interior vertex in the carrier mesh.
            boundary=Dict{NTuple{3,Int},Tuple{Int,Int}}()
            sizehint!(boundary,min(nnodes(mesh),boundary_capacity))
            for entry in entries
                for (key,owner) in entry
                    boundary[key]=min(get(boundary,key,(4,0)),owner)
                end
            end
            for node in eachindex(owners)
                key=_geo_mesh_coordinate_key(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
                owners[node]=get(boundary,key,owners[node])
            end
        end
        result[index]=owners
        entry=Dict{NTuple{3,Int},Tuple{Int,Int}}()
        # Each actual part contributes at most one key per stored vertex.
        sizehint!(entry,nnodes(mesh))
        incidence=dim==0 ? _geo_mesh_point_incidence_coordinate(m,tag) : nothing
        incidence_key=incidence===nothing ? nothing : _geo_mesh_coordinate_key(incidence...)
        for node in eachindex(owners)
            key=_geo_mesh_coordinate_key(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
            dim==0 && incidence_key!==nothing && key!=incidence_key && continue
            entry[key]=min(get(entry,key,(4,0)),owners[node])
        end
        maps[(dim,tag)]=entry
    end
    return result
end

# Merge per-entity meshes: nodes deduplicated by entity ownership and exact coordinates (shared
# boundary nodes are bitwise identical across entity meshes), elements tagged
# with their generating entity in `owner` for node-level classification. A
# `MixedMesh` part (recombined surfaces/volumes) switches the product to the
# mixed merge — quadrangles and polyhedra have no simplex container.
# Geometry keeps its established simplex carrier even when an actual producer
# stores those constructors in MixedMesh blocks. Public native products remain
# in the authority; only this detached geometry merge changes representation.
function _geo_simplex_geometry_part(part::MixedMesh,bindings)
    all(block->block isa ElementBlock && (isempty(block.nodes) || block.msh in (1,2,4)),part.blocks) ||
        throw(ArgumentError("execute_geo: geometry part contains non-simplex constructors"))
    bindings===nothing || length(bindings)==length(part.blocks) || throw(ArgumentError(
        "execute_geo: geometry constructor blocks do not align"))
    matrices=Matrix{Int32}[];tags=Vector{Int32}[]
    groups=bindings===nothing ? nothing : [Model._model_mesh_cell_binding_type([bindings])[] for _ in 1:3]
    for (kind,width) in ((1,2),(2,3),(4,4))
        count=sum(size(block.nodes,2) for block in part.blocks if block.msh==kind;init=0)
        cells=Matrix{Int32}(undef,width,count);physical=Vector{Int32}(undef,count)
        cursor=0
        for (index,block) in enumerate(part.blocks)
            block.msh==kind || continue
            columns=size(block.nodes,2)
            cells[:,cursor+1:cursor+columns]=block.nodes
            physical[cursor+1:cursor+columns]=block.tags
            if groups!==nothing
                length(bindings[index])==columns || throw(ArgumentError(
                    "execute_geo: simplex geometry constructor identities do not align"))
                append!(groups[length(matrices)+1],bindings[index])
            end
            cursor+=columns
        end
        push!(matrices,cells);push!(tags,physical)
    end
    mesh=Mesh(part.coords;segs=matrices[1],tris=matrices[2],tets=matrices[3],
              seg_tag=tags[1],tri_tag=tags[2],tet_tag=tags[3])
    return mesh,groups
end


function _geo_geometry_carrier_merge(parts,node_entities,node_identities,cell_identities,caller)
    simplex=all(part->part[3] isa Mesh ||
        all(block->block isa ElementBlock && (isempty(block.nodes) || block.msh in (1,2,4)),
            part[3].blocks),parts)
    if simplex && any(part->part[3] isa MixedMesh,parts)
        converted=Tuple{Int,Int,Mesh}[]
        binding_type=Model._model_mesh_cell_binding_type(cell_identities)
        converted_bindings=cell_identities===nothing ? nothing : Vector{Vector{binding_type}}[]
        for (index,(dimension,tag,part)) in enumerate(parts)
            bindings=cell_identities===nothing ? nothing : cell_identities[index]
            geometry,groups=part isa Mesh ? (part,bindings) : _geo_simplex_geometry_part(part,bindings)
            push!(converted,(dimension,tag,geometry))
            converted_bindings===nothing ||
                push!(converted_bindings,[Vector{binding_type}(group) for group in groups])
        end
        parts=converted;cell_identities=converted_bindings
    end
    merged,_,rows,cells=_geo_merge_entity_meshes(parts;
        node_entities,node_identities,cell_identities,return_remaps=true)
    return merged,rows,cells
end

function _geo_merge_entity_meshes(parts;node_entities=nothing,node_identities=nothing,
                                  cell_identities=nothing,return_remaps=false)
    caller="execute_geo: Mesh"
    all(part->part[3] isa Mesh,parts) ||
        return _geo_merge_entity_meshes_mixed(parts,caller;
            node_entities,node_identities,cell_identities,return_remaps)
    remaps=Vector{Int32}[]
    coord_list=NTuple{3,Float64}[]
    lookup=Dict{Tuple{Int,Int,NTuple{3,Int}},Int}()
    owner=Dict{Int,Tuple{Int,Int}}()
    nseg=sum(part->nsegs(part[3]),parts)
    ntri=sum(part->ntris(part[3]),parts)
    ntet=sum(part->ntets(part[3]),parts)
    segs=Matrix{Int32}(undef,2,nseg);seg_tag=Vector{Int32}(undef,nseg)
    tris=Matrix{Int32}(undef,3,ntri);tri_tag=Vector{Int32}(undef,ntri)
    tets=Matrix{Int32}(undef,4,ntet);tet_tag=Vector{Int32}(undef,ntet)
    seg_at=tri_at=tet_at=0
    binding_type=Model._model_mesh_cell_binding_type(cell_identities)
    bindings=cell_identities===nothing ? nothing : [binding_type[],binding_type[],binding_type[]]
    function node_index(x,y,z,entity)
        key=(entity...,_geo_mesh_coordinate_key(x,y,z))
        i=get(lookup,key,0)
        i==0 || return i
        i=length(coord_list)+1
        lookup[key]=i
        push!(coord_list,(x,y,z))
        return i
    end
    for (index,(dim,tag,mesh)) in enumerate(parts)
        if cell_identities!==nothing
            length(cell_identities[index])==3 && all(length(tags)==size(cells,2) for
                (tags,cells) in zip(cell_identities[index],(mesh.segs,mesh.tris,mesh.tets))) ||
                throw(ArgumentError("$caller: produced cell identities do not align with geometry"))
            for bi in 1:3
                append!(bindings[bi],cell_identities[index][bi])
            end
        end
        node_identities===nothing || length(node_identities[index])==nnodes(mesh) ||
            throw(ArgumentError("$caller: produced vertex identities do not align with geometry rows"))
        remap=Vector{Int}(undef,nnodes(mesh))
        pcoords=mesh.coords
        for n in 1:nnodes(mesh)
            entity=node_entities===nothing ? (-1,0) : node_entities[index][n]
            identity=node_identities===nothing ? nothing : node_identities[index][n]
            identity===nothing || identity>0 || throw(ArgumentError(
                "$caller: produced vertex identity must be positive"))
            key=identity===nothing ? entity : (typemin(Int),Int(identity))
            remap[n]=node_index(pcoords[1,n],pcoords[2,n],pcoords[3,n],key)
            if node_entities!==nothing
                own=node_entities[index][n]
                owner[remap[n]]=min(get(owner,remap[n],(4,0)),own)
            end
        end
        return_remaps && push!(remaps,Int32.(remap))
        if dim==0
            # Vertex nodes classify on the vertex even when the part carries
            # no elements — upstream `Mesh0D` stores each vertex's mesh
            # vertex on its GVertex regardless of MPoint emission.
            for n in remap
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(0,tag))
            end
        end
        @inbounds for s in 1:nsegs(mesh)
            a,b=remap[mesh.segs[1,s]],remap[mesh.segs[2,s]]
            segs[1,seg_at+=1]=a;segs[2,seg_at]=b
            seg_tag[seg_at]=mesh.seg_tag[s]
            for n in (a,b)
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(dim,tag))
            end
        end
        @inbounds for t in 1:ntris(mesh)
            a,b,c=(remap[mesh.tris[k,t]] for k in 1:3)
            tris[1,tri_at+=1]=a;tris[2,tri_at]=b;tris[3,tri_at]=c
            tri_tag[tri_at]=mesh.tri_tag[t]
            for n in (a,b,c)
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(dim,tag))
            end
        end
        @inbounds for t in 1:ntets(mesh)
            tet_at+=1
            @inbounds for k in 1:4
                n=remap[mesh.tets[k,t]]
                tets[k,tet_at]=n
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(dim,tag))
            end
            tet_tag[tet_at]=mesh.tet_tag[t]
        end
    end
    coords=Matrix{Float64}(undef,3,length(coord_list))
    for (i,(x,y,z)) in enumerate(coord_list)
        coords[1,i]=x;coords[2,i]=y;coords[3,i]=z
    end
    merged=Mesh(coords;segs=segs,tris=tris,tets=tets,
                seg_tag=seg_tag,tri_tag=tri_tag,tet_tag=tet_tag)
    return bindings===nothing ? (return_remaps ? (merged,owner,remaps) : (merged,owner)) :
        return_remaps ? (merged,owner,remaps,bindings) : (merged,owner,bindings)
end

# Mixed counterpart: parts carrying quadrangle/hexahedron/prism blocks merge
# into a `MixedMesh` whose blocks bucket by (MSH type, entity dim, entity tag)
# in first-seen order, so each entity's cells serialize to their own block
# with the entity tag as the elementary field — upstream's per-entity store.
# Simplex parts contribute their segs/tris/tets as line/triangle/tetrahedron
# cells; `tags` keep each cell's physical tag and `owner` keeps the same
# lowest-dimension node classification the simplex merge computes.
function _geo_merge_entity_meshes_mixed(parts,caller::AbstractString;
                                       node_entities=nothing,node_identities=nothing,
                                       cell_identities=nothing,return_remaps=false)
    remaps=Vector{Int32}[]
    coord_list=NTuple{3,Float64}[]
    lookup=Dict{Tuple{Int,Int,NTuple{3,Int}},Int}()
    owner=Dict{Int,Tuple{Int,Int}}()
    # Lower-dimensional parts often repeat the largest retained part's nodes.
    # Use that part to reserve the initial merge tables and coordinate vector.
    capacity=maximum(part->nnodes(part[3]),parts;init=0)
    sizehint!(lookup,capacity)
    sizehint!(owner,capacity)
    sizehint!(coord_list,capacity)
    order=Tuple{Int,Int,Int}[]
    buckets=Dict{Tuple{Int,Int,Int},Tuple{Vector{Int32},Vector{Int32},Int}}()
    binding_type=Model._model_mesh_cell_binding_type(cell_identities)
    bindings=cell_identities===nothing ? nothing : Dict{Tuple{Int,Int,Int},Vector{binding_type}}()
    function node_index(x,y,z,entity)
        key=(entity...,_geo_mesh_coordinate_key(x,y,z))
        i=get(lookup,key,0)
        i==0 || return i
        i=length(coord_list)+1
        lookup[key]=i
        push!(coord_list,(x,y,z))
        return i
    end
    function bucket(msh,dim,tag)
        key=(Int(msh),dim,tag)
        entry=get(buckets,key,nothing)
        entry===nothing || return entry
        width=msh_num_nodes(msh)
        width===nothing && throw(ArgumentError(
            "$caller: cannot merge variable-width MSH type $msh"))
        entry=(Int32[],Int32[],Int(width))
        buckets[key]=entry
        push!(order,key)
        return entry
    end
    for (index,(dim,tag,mesh)) in enumerate(parts)
        if cell_identities!==nothing
            input=Model._model_mesh_carrier_blocks(mesh)
            length(cell_identities[index])==length(input) && all(length(tags)==size(cells,2)
                for (tags,(_,cells)) in zip(cell_identities[index],input)) ||
                throw(ArgumentError("$caller: produced cell identities do not align with geometry"))
        end
        node_identities===nothing || length(node_identities[index])==nnodes(mesh) ||
            throw(ArgumentError("$caller: produced vertex identities do not align with geometry rows"))
        remap=Vector{Int}(undef,nnodes(mesh))
        pcoords=mesh.coords
        for n in 1:nnodes(mesh)
            entity=node_entities===nothing ? (-1,0) : node_entities[index][n]
            identity=node_identities===nothing ? nothing : node_identities[index][n]
            identity===nothing || identity>0 || throw(ArgumentError(
                "$caller: produced vertex identity must be positive"))
            key=identity===nothing ? entity : (typemin(Int),Int(identity))
            remap[n]=node_index(pcoords[1,n],pcoords[2,n],pcoords[3,n],key)
            if node_entities!==nothing
                own=node_entities[index][n]
                owner[remap[n]]=min(get(owner,remap[n],(4,0)),own)
            end
        end
        return_remaps && push!(remaps,Int32.(remap))
        if dim==0
            for n in remap
                old=get(owner,n,(3,0))
                dim<=old[1] && (owner[n]=(0,tag))
            end
        end
        if mesh isa Mesh
            for (bi,(msh,nodes,tags)) in enumerate(((1,mesh.segs,mesh.seg_tag),
                                     (2,mesh.tris,mesh.tri_tag),
                                     (4,mesh.tets,mesh.tet_tag)))
                size(nodes,2)==0 && continue
                flat,cell_tags,width=bucket(msh,dim,tag)
                size(nodes,1)==width || throw(ErrorException(
                    "$caller: MSH type $msh node width $width does not " *
                    "match $(size(nodes,1))"))
                base=length(flat)
                ncells=size(nodes,2)
                resize!(flat,base+ncells*width)
                @inbounds for cell in 1:ncells,row in 1:width
                    flat[base+(cell-1)*width+row]=Int32(remap[nodes[row,cell]])
                end
                append!(cell_tags,tags)
                bindings===nothing || append!(get!(()->binding_type[],bindings,(msh,dim,tag)),
                    cell_identities[index][bi])
                @inbounds for cell in 1:ncells,row in 1:width
                    n=remap[nodes[row,cell]]
                    old=get(owner,n,(3,0))
                    dim<=old[1] && (owner[n]=(dim,tag))
                end
            end
        else
            for (bi,block) in enumerate(mesh.blocks)
                block isa ElementBlock || throw(ArgumentError(
                    "$caller: cannot merge a special-element block " *
                    "(MSH type $(block.msh)) into the model mesh"))
                cells=block.nodes
                size(cells,2)==0 && continue
                flat,cell_tags,width=bucket(block.msh,dim,tag)
                size(cells,1)==width || throw(ErrorException(
                    "$caller: MSH type $(block.msh) node width $width " *
                    "does not match $(size(cells,1))"))
                base=length(flat)
                ncells=size(cells,2)
                resize!(flat,base+ncells*width)
                @inbounds for cell in 1:ncells,row in 1:width
                    flat[base+(cell-1)*width+row]=Int32(remap[cells[row,cell]])
                end
                append!(cell_tags,block.tags)
                bindings===nothing || append!(get!(()->binding_type[],bindings,(block.msh,dim,tag)),
                    cell_identities[index][bi])
                @inbounds for cell in 1:ncells,row in 1:width
                    n=remap[cells[row,cell]]
                    old=get(owner,n,(3,0))
                    dim<=old[1] && (owner[n]=(dim,tag))
                end
            end
        end
    end
    blocks=ElementBlock[]
    elementary=Vector{Vector{Int32}}()
    for (msh,dim,tag) in order
        flat,tags,width=buckets[(msh,dim,tag)]
        nodes=reshape(flat,width,:)
        push!(blocks,ElementBlock(msh,nodes,tags))
        push!(elementary,fill(Int32(tag),size(nodes,2)))
    end
    coords=Matrix{Float64}(undef,3,length(coord_list))
    for (i,(x,y,z)) in enumerate(coord_list)
        coords[1,i]=x;coords[2,i]=y;coords[3,i]=z
    end
    merged=MixedMesh(coords,blocks;elementary_entities=elementary)
    output=bindings===nothing ? nothing : [bindings[key] for key in order]
    return bindings===nothing ? (return_remaps ? (merged,owner,remaps) : (merged,owner)) :
        return_remaps ? (merged,owner,remaps,output) : (merged,owner,output)
end

# Cell count of a given dimension on either mesh representation — the
# `ntris`/`ntets` presence checks the mesh-operation guards need.
_geo_mesh_dim_cells(mesh::Mesh,dim::Int)=
    dim==2 ? ntris(mesh) : dim==3 ? ntets(mesh) :
    dim==1 ? nsegs(mesh) : 0
function _geo_mesh_dim_cells(mesh::MixedMesh,dim::Int)
    total=0
    for block in mesh.blocks
        block isa ElementBlock || continue
        msh_dimension(block.msh)==dim || continue
        total+=size(block.nodes,2)
    end
    return total
end

# Group triangle rows by actual Surface identity before pairing. Actual
# boundary and embedded Curve cells protect only their owning Surface.
function _geo_recombine_mixed_blocks(mesh::MixedMesh;algorithm::Symbol,
                                     recombine_angles=nothing,model=nothing)
    out=MixedElementBlock[]
    elementary=mesh.elementary_entities
    elementary_out=elementary===nothing ? nothing : Vector{Int32}[]
    did_recombine=false
    for (bi,block) in enumerate(mesh.blocks)
        if block isa ElementBlock && Int(block.msh)==2
            # Generated mixed products carry actual elementary entity rows;
            # physical block tags do not identify the generating Surface.
            ncells=size(block.nodes,2)
            ncells==0 && begin
                push!(out,block)
                elementary_out!==nothing && push!(elementary_out,elementary[bi])
                continue
            end
            owners=elementary!==nothing ? elementary[bi] : mesh.entity_data!==nothing ?
                mesh.entity_data.block_entities[bi] : fill(Int32(0),ncells)
            order=Int32[];groups=Dict{Int32,Vector{Int}}()
            for cell in 1:ncells
                owner=owners[cell]
                if !haskey(groups,owner)
                    push!(order,owner);groups[owner]=Int[]
                end
                push!(groups[owner],cell)
            end
            for owner in order
                positions=groups[owner]
                angle=_geo_recombine_entity_angle(recombine_angles,Int(owner))
                sub=Mesh(mesh.coords;tris=block.nodes[:,positions],tri_tag=block.tags[positions])
                protected=_geo_recombine_mixed_edges(mesh,model,Int(owner))
                rec=recombine_triangles(sub;algorithm=algorithm,recombine_angle=angle,
                                        protected_edges=protected)
                for rblock in rec.blocks
                    rblock isa ElementBlock && rblock.msh==3 &&
                        size(rblock.nodes,2)>0 && (did_recombine=true)
                    rblock isa ElementBlock && rblock.msh in (2,3) || continue
                    push!(out,rblock)
                    elementary_out!==nothing && push!(elementary_out,
                        fill(owner,size(rblock.nodes,2)))
                end
            end
        else
            push!(out,block)
            elementary_out!==nothing &&
                push!(elementary_out,elementary[bi])
        end
    end
    did_recombine || return mesh
    changed=length(out)!=length(mesh.blocks)
    if !changed
        for index in eachindex(out)
            old=mesh.blocks[index];new=out[index]
            old===new && continue
            if typeof(old)!==typeof(new) || old.msh!=new.msh ||
                    old.nodes!=new.nodes || old.tags!=new.tags
                changed=true;break
            end
        end
    end
    if changed && (mesh.entity_data!==nothing ||
            any(section->!isempty(section.elements),mesh.data_sections) ||
            (mesh.partition_data!==nothing && !isempty(mesh.partition_data.ghost_elements)))
        throw(ArgumentError("execute_geo: RecombineMesh cannot preserve " *
            "cell-attached MSH metadata while changing triangle blocks"))
    end
    return MixedMesh(mesh.coords,out;
                     physical_names=mesh.physical_names,
                     entity_data=mesh.entity_data,
                     elementary_entities=elementary_out,
                     periodic_links=mesh.periodic_links,
                     ancillary_sections=mesh.ancillary_sections,
                     data_sections=mesh.data_sections,
                     partition_data=mesh.partition_data)
end

function _geo_recombine_curve_tags(m,owner::Int)
    owner>0 || throw(ArgumentError(
        "execute_geo: RecombineMesh lacks elementary Surface identity " *
        "for actual Curve exclusions"))
    curves=Set{Int32}()
    for (dim,tag) in _geo_mesh_children(m,2,owner,"execute_geo: RecombineMesh")
        dim==1 && push!(curves,Int32(abs(tag)))
    end
    record=get(m.discrete,(2,owner),nothing)
    if record!==nothing
        for (dim,tag) in record.boundary
            dim==1 && push!(curves,Int32(abs(tag)))
        end
    end
    return curves
end

function _geo_recombine_mixed_edges(mesh::MixedMesh,m,owner::Int)
    m===nothing && return nothing
    curves=_geo_recombine_curve_tags(m,owner)
    protected=Set{NTuple{2,Int32}}()
    isempty(curves) && return protected
    elementary=mesh.elementary_entities;data=mesh.entity_data
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock && msh_family(block.msh)===:lin || continue
        owners=elementary!==nothing ? elementary[bi] : data!==nothing ?
            data.block_entities[bi] : nothing
        owners===nothing && !isempty(block.nodes) && throw(ArgumentError(
            "execute_geo: RecombineMesh lacks actual Curve element identities"))
        for cell in axes(block.nodes,2)
            owners[cell] in curves || continue
            push!(protected,minmax(block.nodes[1,cell],block.nodes[2,cell]))
        end
    end
    return protected
end

function _geo_recombine_angles(m::GeoModel)
    angles=Dict{Tuple{Int,Int},Float64}()
    for ((dim,tag),angle) in m.meshing.recombine
        dim==2 && (angles[(2,tag)]=angle)
    end
    for tag in keys(m.surfaces);get!(angles,(2,tag),45.0);end
    for ((dim,tag),_) in m.discrete
        dim==2 && get!(angles,(2,tag),45.0)
    end
    return angles
end

# Unknown imported triangle owners can use a common angle; an actual
# per-Surface difference requires explicit elementary entity provenance.
function _geo_recombine_entity_angle(angles,owner::Int)
    angles===nothing && return nothing
    owner>0 && return get(angles,(2,owner),45.0)
    angle_values=Float64[value for ((dim,_),value) in angles if dim==2]
    isempty(angle_values) && return 45.0
    all(value->value==first(angle_values),angle_values) || throw(ArgumentError(
        "execute_geo: RecombineMesh requires elementary Surface identities " *
        "to apply different recombination angles"))
    return first(angle_values)
end

function _geo_recombine_simplex_blocks(mesh::Mesh,m::GeoModel,context;
                                       algorithm::Symbol)
    angles=_geo_recombine_angles(m)
    # Preserve the existing native rejection of volumetric simplex inputs.
    ntets(mesh)>0 && return recombine_triangles(mesh;algorithm=algorithm)
    parts=context.mesh_parts
    isempty(parts) && throw(ArgumentError(
        "execute_geo: RecombineMesh lacks generating Surface identities"))
    all(part->part[3] isa Mesh,parts) || throw(ArgumentError(
        "execute_geo: inconsistent simplex generating products"))
    node_entities=_geo_mesh_part_node_entities(m,parts,"execute_geo: RecombineMesh")
    source,_=_geo_merge_entity_meshes(parts;node_entities=node_entities)
    owners=Int32[];segment_owners=Int32[]
    for (dim,tag,part) in parts
        append!(owners,fill(Int32(dim==2 ? tag : 0),ntris(part)))
        append!(segment_owners,fill(Int32(dim==1 ? tag : 0),nsegs(part)))
    end
    # RefineMesh has deterministic parent-contiguous children. Reproduce only
    # its topology and verify every dense slot before accepting the inheritance.
    while ntris(source)>0 && ntris(source)<ntris(mesh) &&
          ntris(source)<=div(ntris(mesh),4)
        source=refine_uniform(source)
        owners=repeat(owners;inner=4)
        segment_owners=repeat(segment_owners;inner=2)
    end
    source.tris==mesh.tris && source.segs==mesh.segs &&
        length(owners)==ntris(mesh) && all(>(0),owners) || throw(ArgumentError(
            "execute_geo: RecombineMesh cannot resolve current triangle " *
            "Surface identities from generating products"))
    blocks=MixedElementBlock[];elementary=Vector{Int32}[]
    did_recombine=false
    if nsegs(mesh)>0
        push!(blocks,ElementBlock(1,mesh.segs,mesh.seg_tag))
        push!(elementary,segment_owners)
    end
    order=Int32[];groups=Dict{Int32,Vector{Int}}()
    for cell in eachindex(owners)
        owner=owners[cell]
        if !haskey(groups,owner)
            push!(order,owner);groups[owner]=Int[]
        end
        push!(groups[owner],cell)
    end
    for owner in order
        positions=groups[owner]
        curves=_geo_recombine_curve_tags(m,Int(owner))
        protected=Set{NTuple{2,Int32}}()
        for cell in axes(mesh.segs,2)
            segment_owners[cell] in curves || continue
            push!(protected,minmax(mesh.segs[1,cell],mesh.segs[2,cell]))
        end
        sub=Mesh(mesh.coords;segs=mesh.segs,seg_tag=mesh.seg_tag,
                 tris=mesh.tris[:,positions],tri_tag=mesh.tri_tag[positions])
        rec=recombine_triangles(sub;algorithm=algorithm,
            recombine_angle=get(angles,(2,Int(owner)),45.0),
            protected_edges=protected)
        for block in rec.blocks
            block isa ElementBlock && block.msh==3 &&
                size(block.nodes,2)>0 && (did_recombine=true)
            block isa ElementBlock && block.msh in (2,3) || continue
            push!(blocks,block);push!(elementary,fill(owner,size(block.nodes,2)))
        end
    end
    did_recombine || return mesh
    return MixedMesh(mesh.coords,blocks;elementary_entities=elementary)
end

# The mesh-operation statements act on `context.mesh` — the product of a
# mid-file `Mesh n` — or record meshing constraints for the next generation.
function _geo_exec_mesh_statement!(m::GeoModel,line::AbstractString,
                                   context::_GeoNumericContext,
                                   allocator_state=nothing)
    caller="execute_geo"
    s=String(strip(line))
    if (mm=match(r"^RefineMesh\s*;?\s*$",s))!==nothing
        # `RefineMesh` syncs the internals first when changed.
        _geo_sync_physical_view_if_changed!(m,context)
        context.mesh===nothing && return nothing
        if !isempty(context.mesh_parts) ||
                (context.mesh isa MixedMesh && context.mesh.entity_data!==nothing)
            _geo_operate_actual_mesh!(m,context,:refine,"$caller: RefineMesh")
            return nothing
        end
        context.mesh isa Mesh || throw(ArgumentError(
            "$caller: RefineMesh requires a simplex mesh — the current " *
            "mesh carries non-simplex elements"))
        context.mesh=refine_uniform(context.mesh)
        _geo_remesh_owner!(context)
        return
    elseif (mm=match(r"^RecombineMesh\s*;?\s*$",s))!==nothing
        _geo_sync_meshing_options!(m,context)
        if m.meshing.recombine_algo==4
            _geo_recombine_bipartite_plan!(m,context,"$caller: RecombineMesh")
            return nothing
        end
        context.mesh===nothing && return nothing
        if !isempty(context.mesh_parts) ||
                (context.mesh isa MixedMesh && context.mesh.entity_data!==nothing)
            _geo_recombine_actual_plan!(m,context,"$caller: RecombineMesh")
            return nothing
        end
        _geo_sync_meshing_options!(m,context)
        # Upstream `RecombineMesh` (Generator.cpp) reads `algoRecombine`
        # live: blossom matching for 1/3, simple pairing otherwise — the
        # statement has no full-quad arm.
        algorithm=m.meshing.recombine_algo in (1,3) ? :blossom : :greedy
        context.mesh=if context.mesh isa Mesh
            _geo_recombine_simplex_blocks(context.mesh,m,context;algorithm=algorithm)
        else
            angles=_geo_recombine_angles(m)
            _geo_recombine_mixed_blocks(context.mesh;algorithm=algorithm,
                                        recombine_angles=angles,model=m)
        end
        _geo_remesh_owner!(context)
        return
    elseif (mm=match(r"^RenumberMeshNodes\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^RenumberMeshElements\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^CreateMeshEdges\s*;?\s*$",s))!==nothing ||
           (mm=match(r"^CreateMeshFaces\s*;?\s*$",s))!==nothing
        authority=m.meshing.mesh_authority
        if context.mesh isa MixedMesh && context.mesh.entity_data!==nothing ||
                authority!==nothing && authority.mesh!==nothing && authority.carrier===context.mesh
            if startswith(s,"RenumberMeshNodes")
                _geo_operate_actual_mesh!(m,context,:geo_renumber_nodes,"$caller: RenumberMeshNodes")
                context.native_element_lookup_receipt=nothing
            elseif startswith(s,"RenumberMeshElements")
                _geo_operate_actual_mesh!(m,context,:geo_renumber_elements,"$caller: RenumberMeshElements")
                context.native_element_lookup_receipt=nothing
            end
        end
        # Node/element order and topology catalogs are derived information in
        # the compact Mesh — already canonical, nothing to mutate.
        return
    elseif (mm=match(r"^ReorientMesh\s+Volume\s*\{(.*?)\}\s*;?\s*$",s))!==nothing
        _geo_exec_reorient_mesh!(m,"{"*mm.captures[1]*"}",context)
        return
    elseif (mm=match(r"^RelocateMesh\s+(Point|Curve|Line|Surface)\s*(.*?)\s*;?\s*$",s))!==nothing
        _geo_exec_relocate_mesh!(m,mm.captures[1],mm.captures[2],context)
        return
    elseif (mm=match(r"^TransformMesh\s*\{(.*?)\}\s*(?:\{(.*?)\})?\s*;?\s*$",s))!==nothing
        context.mesh===nothing && return nothing
        matrix=_geo_numeric_list_values("{"*mm.captures[1]*"}",context,caller)
        length(matrix)>=12 || throw(ArgumentError(
            "$caller: affine transform matrix requires at least 12 entries"))
        wanted=nothing
        if mm.captures[2]!==nothing
            entities=_geo_shape_list_entities!(
                m,mm.captures[2],context,"$caller: TransformMesh";
                signed_tags=true,allocator_state=allocator_state)
            wanted=Set(Tuple{Int,Int}(abs.(Int.(collect(e)))) for e in entities)
        end
        source,owner,_=_geo_mesh_coordinate_source(m,context)
        _geo_check_part_rows(context,source)
        edited=deepcopy(source)
        _geo_transform_mesh_nodes!(edited,owner,matrix,wanted)
        _geo_publish_mesh_coordinates!(m,context,edited.coords,"$caller: TransformMesh")
        return
    elseif (mm=match(r"^SetOrder\s+(.*?)\s*;?\s*$",s))!==nothing
        order=_geo_int_value(
            _geo_eval_numeric(mm.captures[1],context,caller),caller)
        order<1 && throw(ArgumentError("$caller: SetOrder requires order >= 1"))
        if !isempty(context.mesh_parts) ||
                context.mesh isa MixedMesh && context.mesh.entity_data!==nothing
            _geo_operate_actual_mesh!(m,context,:order,"$caller: SetOrder";element_order=order)
            return nothing
        end
        order>1 && throw(ArgumentError(
            "$caller: SetOrder $order — high-order .geo meshes are not " *
            "supported"))
        return
    elseif (mm=match(r"^PartitionMesh\b\s*(.+?)\s*;?\s*$",s))!==nothing
        _geo_exec_partition_mesh!(m,mm.captures[1],context,caller)
        return
    elseif match(r"^CreateOverlaps\b",s)!==nothing
        # `CreateOverlaps` builds inter-partition overlap element layers — a
        # separate data flow upstream (`GModel::createOverlap*` builds them
        # through the parallel-mesh machinery, not the partition producers),
        # so there is no existing Tessella producer to reuse. Keep the
        # explicit blocker rather than a silent partial implementation.
        throw(ArgumentError(
            "$caller: CreateOverlaps is not supported"))
    elseif (mm=match(r"^AdaptMesh\b",s))!==nothing
        throw(ArgumentError(
            "$caller: AdaptMesh requires levelset fields, which are not " *
            "supported"))
    elseif (mm=match(r"^CreateTopology\b(.*)$",s))!==nothing
        _geo_exec_create_topology!(m,mm.captures[1],context,caller)
        return
    elseif (mm=match(r"^ClassifySurfaces\b(.*)$",s))!==nothing
        _geo_exec_classify_surfaces!(m,mm.captures[1],context,caller)
        return
    elseif (mm=match(r"^CreateGeometry\b(.*)$",s))!==nothing
        _geo_exec_create_geometry!(m,mm.captures[1],context,caller,
                                   allocator_state)
        return
    end
    throw(ArgumentError("$caller: unrecognized mesh command: $line"))
end

# Statements taking a balanced `{...}` argument must end right after it —
# trailing text is the token the upstream parser would flag next.
function _geo_require_statement_end(rest::AbstractString,
                                    caller::AbstractString)
    r=strip(rest)
    (isempty(r) || r==";") && return nothing
    _geo_syntax_abort(_geo_first_token(r))
end

# Fixed-arity `{...}` argument lists fail upstream while parsing — `}` when
# the list is short, `,` when it is long.
function _geo_arity_abort(count::Int,allowed,caller::AbstractString)
    count in allowed && return nothing
    _geo_syntax_abort(count<minimum(allowed) ? "}" : ",")
end

# `CreateTopology [ { makeSimplyConnected, exportDiscrete } ]` — derive a
# boundary representation from the discrete entities like
# `GModel::createTopologyFromMesh` (plus the simply-connected repair and the
# GEO-internals export in the braced form).
function _geo_exec_create_topology!(m::GeoModel,rest::AbstractString,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    caller="$caller: CreateTopology"
    s=String(strip(rest))
    if isempty(s) || s==";"
        create_topology!(m)
        _geo_export_discrete_internals!(m,context)
        return nothing
    end
    (content,rest)=_geo_balanced_group(s,caller)
    values=_geo_numeric_list_values("{"*content*"}",context,caller)
    _geo_arity_abort(length(values),(2,),caller)
    _geo_require_statement_end(rest,caller)
    create_topology!(m;make_simply_connected=!iszero(values[1]),
                     export_discrete=!iszero(values[2]))
    iszero(values[2]) || _geo_export_discrete_internals!(m,context)
    return nothing
end

# `GModel::exportDiscreteGEOInternals` — "Warning: this clears
# GEO_Internals!" — rebuilds the internals from the model after
# `createTopologyFromMesh`. Tessella's geometric entities live in the model
# eagerly, so the only pending parser state is `raw_physicals`: pending
# declarations that were never synchronized are dropped upstream, and the
# registry is rebuilt from the entities' current physical memberships.
function _geo_export_discrete_internals!(m::GeoModel,
                                       context::_GeoNumericContext)
    empty!(context.raw_physicals)
    for ((dim,entity),pnums) in m.entity_physicals
        for p in pnums
            push!(get!(() -> Int[],context.raw_physicals,(dim,p)),entity)
        end
    end
    # The derived `m.physical` view is stale until the next sync point.
    context.geo_changed=true
    return nothing
end

# `ClassifySurfaces { angle, boundary, forReparametrization [, curveAngle] }`
# — `GModel::classifySurfaces`, which also runs `createTopologyFromMesh` when
# the discrete model has no boundary representation yet.
function _geo_exec_classify_surfaces!(m::GeoModel,rest::AbstractString,
                                      context::_GeoNumericContext,
                                      caller::AbstractString)
    caller="$caller: ClassifySurfaces"
    s=String(strip(rest))
    startswith(s,"{") || _geo_syntax_abort(_geo_first_token(s))
    (content,rest)=_geo_balanced_group(s,caller)
    values=_geo_numeric_list_values("{"*content*"}",context,caller)
    _geo_arity_abort(length(values),(3,4),caller)
    _geo_require_statement_end(rest,caller)
    classify_surfaces!(m;angle=values[1],boundary=!iszero(values[2]),
                       for_reparametrization=!iszero(values[3]),
                       curve_angle=length(values)==4 ? values[4] : pi)
    return nothing
end

# `CreateGeometry [ { entities } ]` — `GModel::createGeometryOfDiscreteEntities`,
# on every discrete curve/surface when the entity list is omitted.
function _geo_exec_create_geometry!(m::GeoModel,rest::AbstractString,
                                    context::_GeoNumericContext,
                                    caller::AbstractString,allocator_state)
    caller="$caller: CreateGeometry"
    s=String(strip(rest))
    if isempty(s) || s==";"
        create_geometry!(m)
        return nothing
    end
    (content,rest)=_geo_balanced_group(s,caller)
    entities=_geo_shape_list_entities!(m,content,context,caller;
                                       allocator_state=allocator_state)
    _geo_require_statement_end(rest,caller)
    create_geometry!(m,entities)
    return nothing
end

# `Homology`/`Cohomology`/`Betti` queue a homology request — the upstream
# `addHomologyRequest` calls in parser/Gmsh.y. Three forms:
#   `Cmd;`                       — empty domain/subdomain, dims 0:3
#   `Cmd { dom [, sub] };`       — one or two `ListOfDouble`s
#   `Cmd (dims) { dom, sub };`   — dims then both lists (required upstream)
function _geo_exec_homology_statement!(m::GeoModel,word::AbstractString,
                                       rest::AbstractString,
                                       context::_GeoNumericContext)
    caller="execute_geo: $word"
    s=String(strip(rest))
    dims=Int[0,1,2,3]
    if startswith(s,"(")
        (inner,s)=_geo_balanced_paren(s,caller)
        dims=_geo_homology_int_list(inner,context,caller)
        (content,s)=_geo_balanced_group(s,caller)
        (domain,subdomain)=_geo_homology_domains(
            content,context,caller;require_two=true)
        _geo_require_statement_end(s,caller)
    elseif isempty(s) || s==";"
        domain=Int[];subdomain=Int[]
    elseif startswith(s,"{")
        (content,s)=_geo_balanced_group(s,caller)
        (domain,subdomain)=_geo_homology_domains(content,context,caller)
        _geo_require_statement_end(s,caller)
    else
        _geo_syntax_abort(_geo_first_token(s))
    end
    add_homology_request!(m;kind=String(word),domain_tags=domain,
                          subdomain_tags=subdomain,dims=dims)
    return nothing
end

# `'{' ListOfDouble '}'` / `'{' ListOfDouble ',' ListOfDouble '}'` — a second
# list only follows a complete braced/negated/scaled first list; a bare
# `a,b` always parses as one flat `RecursiveListOfDouble` upstream.
function _geo_homology_domains(content::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString;require_two::Bool=false)
    s=String(strip(content))
    isempty(s) &&
        (require_two && _geo_syntax_abort("}"))
    isempty(s) && return (Int[],Int[])
    (first_src,rest)=_geo_homology_list_arg(s,caller)
    domain=_geo_homology_int_list(first_src,context,caller)
    if isempty(rest)
        require_two && _geo_syntax_abort("}")
        return (domain,Int[])
    end
    startswith(rest,",") || _geo_syntax_abort(_geo_first_token(rest))
    (second_src,rest)=_geo_homology_list_arg(
        String(strip(rest[nextind(rest,firstindex(rest)):end])),caller)
    _geo_require_statement_end(rest,caller)
    return (domain,_geo_homology_int_list(second_src,context,caller))
end

# One `ListOfDouble` at the head of `s` → `(source, rest)`. A `{...}` group
# (optionally `-`-negated or `FExpr*`-scaled) ends at its balancing `}`;
# anything else consumes the whole string as one flat list and the shared
# evaluator reports embedded `{` misuse at the same token upstream would.
function _geo_homology_list_arg(s::AbstractString,caller::AbstractString)
    s=String(strip(s))
    isempty(s) && _geo_syntax_abort("}")
    brace=findfirst(==('{'),s)
    brace===nothing && return (s,"")
    prefix=String(strip(s[firstindex(s):prevind(s,brace)]))
    (isempty(prefix) || prefix=="-" || endswith(prefix,"*")) ||
        return (s,"")
    group=String(s[brace:lastindex(s)])
    closing=_geo_braced_end(group)
    closing===nothing && _geo_syntax_abort("{")
    src_end=nextind(s,brace,closing-1)  # index of the balancing `}` in `s`
    src=String(s[firstindex(s):src_end])
    rest=src_end<lastindex(s) ?
        String(strip(s[nextind(s,src_end):end])) : ""
    return (src,rest)
end

# Physical-group tags are C ints upstream — `(int)` truncation of the
# evaluated scalar, with the grammar's int-range check.
function _geo_homology_int_list(raw::AbstractString,
                                context::_GeoNumericContext,
                                caller::AbstractString)
    s=String(strip(raw))
    # A braced/negated/scaled group evaluates verbatim; flat content wraps in
    # braces for the shared `RecursiveListOfDouble` evaluator.
    brace=findfirst(==('{'),s)
    src=if brace===nothing
        "{"*s*"}"
    else
        prefix=String(strip(s[firstindex(s):prevind(s,brace)]))
        (isempty(prefix) || prefix=="-" || endswith(prefix,"*")) ?
            s : "{"*s*"}"
    end
    return Int[_geo_signed_gmsh_int_value(v,caller)
               for v in _geo_numeric_list_values(src,context,caller)]
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

function _geo_mesh_stored_parameters(m,mesh,owners,context)
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    if mesh isa MixedMesh && mesh.entity_data!==nothing
        return Model._model_mesh_runtime_parameters(m,mesh)
    end
    old_ruled_surface=_geo_exec_old_ruled(context)
    for (node,(dimension,tag)) in owners
        dimension in (1,2) || continue
        haskey(_model_entity_dictionary(m,dimension),tag) || continue
        parameters[node]=model_parametrization(m,dimension,tag,
            collect(@view mesh.coords[:,node]);old_ruled_surface)
    end
    return parameters
end

function _geo_check_part_rows(context,mesh=context.mesh)
    length(context.mesh_parts)==length(context.mesh_part_rows) || throw(ArgumentError(
        "execute_geo: entity products do not align with their stored node rows"))
    for ((_,_,part),rows) in zip(context.mesh_parts,context.mesh_part_rows)
        length(rows)==nnodes(part) && all(row->1<=row<=nnodes(mesh),rows) ||
            throw(ArgumentError("execute_geo: invalid stored entity node rows"))
    end
    return nothing
end

function _geo_publish_part_coordinates!(context,mesh=context.mesh)
    _geo_check_part_rows(context,mesh)
    for ((_,_,part),rows) in zip(context.mesh_parts,context.mesh_part_rows)
        for (column,row) in enumerate(rows),axis in 1:3
            part.coords[axis,column]=mesh.coords[axis,row]
        end
    end
    return nothing
end

function _geo_mesh_coordinate_source(m,context)
    previous=m.meshing.mesh_authority
    if previous!==nothing && previous.mesh!==nothing && previous.carrier===context.mesh
        source=previous.mesh
        owner=Dict(row=>(dimension,Int(tag)) for (row,(dimension,tag)) in
            enumerate(source.entity_data.node_entities))
        return source,owner,Model._model_mesh_runtime_parameters(m,source)
    end
    return context.mesh,context.mesh_node_owner,context.mesh_node_parametric
end

function _geo_coordinate_mesh(mesh,coordinates)
    if mesh isa Mesh
        return Mesh(coordinates;segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
            seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
    end
    return MixedMesh(coordinates,mesh.blocks;physical_names=mesh.physical_names,
        entity_data=mesh.entity_data,elementary_entities=mesh.elementary_entities,
        periodic_links=mesh.periodic_links,ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
end

function _geo_publish_mesh_coordinates!(m,context,coordinates,caller)
    source,_,parameters=_geo_mesh_coordinate_source(m,context)
    size(coordinates)==size(source.coords) && all(isfinite,coordinates) ||
        throw(ArgumentError("$caller: coordinates must be finite and align with actual nodes"))
    _geo_check_part_rows(context,source)
    mesh=_geo_coordinate_mesh(source,coordinates)
    previous=m.meshing.mesh_authority
    bound=previous!==nothing && previous.mesh!==nothing && previous.carrier===context.mesh
    carrier=bound && context.mesh!==source ?
        _geo_coordinate_mesh(context.mesh,coordinates[:,previous.carrier_rows]) : mesh
    parts=Matrix{Float64}[coordinates[:,rows] for rows in context.mesh_part_rows]
    record_updates=Tuple{DiscreteEntity,Matrix{Float64}}[]
    authority=nothing
    if bound
        authority=Model._model_mesh_authority(mesh,previous.done;caller,
            private_rank=previous.private_rank,carrier,carrier_rows=previous.carrier_rows)
        inputs=Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}()
        for (owner,input) in previous.native_inputs
            product=Model._model_mesh_product(authority,owner[1],Int(owner[2]);caller)
            inputs[owner]=product===nothing ? input : product
        end
        # The coordinate rebuild shares the previous entity data and tag space,
        # so retained Surface constructor inputs stay bound to the same products.
        # Their coordinates ride the same transform — a receipt carrying stale
        # coordinates would fail its next product binding and get dropped.
        positions=Dict(tag=>row for (row,tag) in enumerate(mesh.entity_data.external_node_tags))
        constructors=Dict{Tuple{Int,Int32},MixedMesh}()
        for (owner,ctor) in previous.surface_constructor_inputs
            updated=copy(ctor.coords);ok=true
            for (column,tag) in enumerate(ctor.entity_data.external_node_tags)
                row=get(positions,tag,0)
                row==0 && (ok=false;break)
                updated[:,column]=coordinates[:,row]
            end
            constructors[owner]= ok ? _geo_coordinate_mesh(ctor,updated) : ctor
        end
        authority=Model._ModelMeshAuthority(authority.mesh,authority.products,
            authority.done,authority.private_rank,inputs,carrier,authority.carrier_rows,
            previous.producer_rows,previous.carrier_cells,constructors)
        for (_,record) in Model._discrete_mesh_records_model(m)
            updated=copy(record.node_coords);changed=false
            for (column,tag) in enumerate(record.node_tags)
                row=get(positions,UInt64(tag),0)
                row==0 && continue
                updated[:,column]=coordinates[:,row];changed=true
            end
            changed && push!(record_updates,(record,updated))
        end
    end
    carrier_parameters=bound ? deepcopy(parameters[authority.carrier_rows]) : deepcopy(parameters)
    geometry_receipt=_geo_geometry_coordinate_receipt(context.mesh_geometry_parts,carrier,caller)
    # Evaluation, metadata validation and all coordinate allocations precede
    # publication to the context and the model's actual entity products.
    for ((_,_,part),updated) in zip(context.mesh_parts,parts)
        part.coords[:,:]=updated
    end
    for (record,updated) in record_updates
        record.node_coords=updated
    end
    bound && (m.meshing.mesh_authority=authority)
    context.mesh=carrier
    context.mesh_geometry_parts=geometry_receipt
    context.mesh_node_parametric=carrier_parameters
    return nothing
end

function _geo_exec_relocate_mesh!(m,kind,list_source,context)
    caller="execute_geo: RelocateMesh $kind"
    list=_geo_constraint_list(list_source,context,caller;allow_all=true)
    _geo_sync_physical_view_if_changed!(m,context)
    dimension=kind=="Point" ? 0 : kind in ("Curve","Line") ? 1 : 2
    wanted=list===nothing ? nothing :
        Set(_geo_constraint_int(value,caller,"entity tag") for value in list)
    mesh,owners,stored_parameters=_geo_mesh_coordinate_source(m,context)
    mesh===nothing && return nothing
    length(stored_parameters)==nnodes(mesh) || throw(ArgumentError(
        "$caller: stored node parameters do not align with the current mesh"))
    _geo_check_part_rows(context,mesh)
    coordinates=copy(mesh.coords)
    for (node,(owner_dimension,tag)) in owners
        owner_dimension==dimension && (wanted===nothing || tag in wanted) || continue
        haskey(_model_entity_dictionary(m,dimension),tag) || continue
        parameters=stored_parameters[node]
        dimension>0 && parameters===nothing && continue
        point=model_value(m,dimension,tag,dimension==0 ? Float64[] : parameters;
            old_ruled_surface=_geo_exec_old_ruled(context))
        length(point)==3 && all(isfinite,point) || throw(ArgumentError(
            "$caller: entity evaluation must produce three finite coordinates"))
        coordinates[:,node]=point
    end
    _geo_publish_mesh_coordinates!(m,context,coordinates,caller)
    return nothing
end

function _geo_transform_mesh_nodes!(mesh,owner,matrix,wanted)
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
        return _geo_exec_mesh_statement!(m,"$word $arg_src;",context,
                                         allocator_state)
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
        # `GModel::optimizeMesh(how)` — the method is validated before any
        # entity iteration, so an unknown name errors even without a mesh;
        # a mesh without the relevant cells is a silent no-op.
        arg in ("","Gmsh","Optimize","Relocate3D","Laplace2D","Relocate2D") ||
            throw(ArgumentError(
                "$caller: unknown or unsupported mesh optimization " *
                "method \"$arg\""))
        context.mesh===nothing && return nothing
        mesh=context.mesh
        if arg in ("","Gmsh","Optimize","Relocate3D")
            _geo_mesh_dim_cells(mesh,3)>0 || return nothing
            mesh isa Mesh || throw(ArgumentError(
                "$caller: OptimizeMesh \"$arg\" requires a simplex mesh — " *
                "the current mesh carries non-simplex elements"))
            replacement=smooth_optimize(mesh;iters=1,require_positive_tets=false)
            receipt=_geo_geometry_coordinate_receipt(context.mesh_geometry_parts,replacement,caller)
            context.mesh=replacement;context.mesh_geometry_parts=receipt
        else
            _geo_mesh_dim_cells(mesh,2)>0 || return nothing
            mesh isa Mesh || throw(ArgumentError(
                "$caller: OptimizeMesh \"$arg\" requires a simplex mesh — " *
                "the current mesh carries non-simplex elements"))
            replacement=_laplacian_smooth_tri_cache(mesh,1,nothing)
            receipt=_geo_geometry_coordinate_receipt(context.mesh_geometry_parts,replacement,caller)
            context.mesh=replacement;context.mesh_geometry_parts=receipt
        end
        return nothing
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

include("GeoMeshImport.jl")

# Lowest free Surface tag for an imported STL discrete surface — the same
# `API.mesh.import_stl` allocation rule, extended over imported mesh entities
# and the live mesh authority so a merged STL can never shadow an entity a
# previous `Include`/`Merge` created.
function _geo_stl_surface_tag(m::GeoModel)
    occupied=Set{Int}()
    for (dim,tag) in Model.model_entities(m)
        dim==2 && push!(occupied,Int(tag))
    end
    for (dim,tag) in keys(m.discrete)
        dim==2 && push!(occupied,Int(tag))
    end
    for (dim,tag) in keys(m.meshing.imported_entities)
        dim==2 && push!(occupied,Int(tag))
    end
    authority=m.meshing.mesh_authority
    if authority!==nothing
        for mesh in (authority.mesh,authority.carrier)
            if mesh isa MixedMesh && mesh.entity_data!==nothing
                for (dim,tag) in keys(mesh.entity_data.entities)
                    dim==2 && push!(occupied,Int(tag))
                end
            end
        end
    end
    tag=1
    while tag in occupied
        tag+=1
    end
    return tag
end

# Parsed `.pos` views for the running execution, keyed on the numeric context.
# `Merge`/`Include "x.pos"` appends in file order (upstream `GModel._views`);
# `Background Mesh View[i]` and `Field[i] = "PostView"` resolve into it.
# Weak keys drop the store with the context, so no per-run cleanup is needed.
const _GEO_EXEC_VIEWS=WeakKeyDict()

# Scalar `PostView` background adapter — upstream's `PostViewField` *is* the
# size field (its scalar operator interpolates the view values). Tessella's
# `PostViewField` is an `AbstractField` sampler, so the `.geo` background
# slot — typed `AbstractSizeField` — needs this façade, the scalar
# counterpart of `PostViewAnisoField`.
struct _GeoPostViewSizeField{F<:AbstractField} <: AbstractSizeField
    field::F
end
field_value(field::_GeoPostViewSizeField,x,y,z)=
    field_value(field.field,x,y,z)
field_value(field::_GeoPostViewSizeField,x,y,z,::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=field_value(field.field,x,y,z)

# `Include`/`Merge` of a `.geo` file runs its statements in the same context
# and model (Gmsh's `ParseFile` shares the parser globals). `.msh` files merge
# their mesh into the current mid-file mesh.
function _geo_exec_include!(m::GeoModel,path::AbstractString,
                            context::_GeoNumericContext,caller::AbstractString,
                            allocator_state)
    file=_geo_fix_relative_path(context.file_name,path)
    isfile(file) || throw(ArgumentError("$caller: missing file '$file'"))
    ext=lowercase(splitext(file)[2])
    if ext==".pos"
        # `Merge`/`Include` of a `.pos` file appends its views to the model's
        # view list (upstream `PView` load order); `Background Mesh View[i]`
        # and `Field[i] = "PostView"` index into the same ordering. The store
        # keys on the execution context — a `WeakKeyDict` drops entries with
        # the run.
        post=getfield(parentmodule(@__MODULE__),:Post)
        store=get!(_GEO_EXEC_VIEWS,context,Any[])
        append!(store,post.read_pos(file))
        return nothing
    end
    if ext==".msh" || ext==".stl"
        # `GModel::readMSH` imports the file's elementary entities as discrete
        # model entities carrying their cells, and the cells join the mid-file
        # mesh as well. `read_mixed_msh` keeps every element family plus the
        # MSH2 elementary / MSH4 entity classification the simplex `read_msh`
        # folds away. `.stl` takes the same channel: `read_stl` welds the
        # facet soup into a surface `Mesh`, a temporary MSH 2.2 bridge encodes
        # it under the lowest free Surface tag, and the importer publishes it
        # as a discrete Surface — matching `GModel::readSTL`'s discreteFace.
        if ext==".stl"
            io_mod=getfield(parentmodule(@__MODULE__),:IO)
            temporary=io_mod._stl_bridge_temp(file,_geo_stl_surface_tag(m))
            try
                _geo_import_mesh_file!(m,context,temporary,caller)
            finally
                rm(temporary;force=true)
            end
            return nothing
        end
        _geo_import_mesh_file!(m,context,file,caller)
        return nothing
    end
    if ext in (".step",".stp")
        # `GModel::readOCCSTEP` — STEP always imports through the OCC
        # kernel upstream, whatever the active factory. The native
        # importer binds into the OCC kernel; synchronize publishes the
        # merged solids immediately (same `geo_changed` receipt a GEO
        # create carries).
        Model._model_occ_import_shapes(m,file;caller="OCC.importShapes")
        _geo_occ_finish!(m,context)
        return nothing
    end
    ext==".geo" || throw(ArgumentError(
        "$caller: merging '$ext' files is not supported " *
        "(only .geo, .msh, .stl, .pos, .step and .stp)"))
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

# `Merge "x.msh"` against a recombined/mixed mid-file mesh: both sides join
# as `MixedMesh` blocks with the simplex part expanded into line/triangle/
# tetrahedron blocks. Per-entity classification is already dropped on this
# path (`mesh_parts`/`mesh_node_owner` reset), so the concat keeps cell tags
# and physical names only.
function _geo_concat_meshes(a,b)
    offset=Int32(nnodes(a))
    coords=hcat(a.coords,b.coords)
    function shifted_blocks(mesh)
        shift=mesh===a ? Int32(0) : offset
        out=ElementBlock[]
        if mesh isa Mesh
            for (msh,nodes,tags) in ((1,mesh.segs,mesh.seg_tag),
                                     (2,mesh.tris,mesh.tri_tag),
                                     (4,mesh.tets,mesh.tet_tag))
                size(nodes,2)>0 &&
                    push!(out,ElementBlock(msh,nodes.+shift,tags))
            end
        else
            for block in mesh.blocks
                block isa ElementBlock || throw(ArgumentError(
                    "execute_geo: cannot merge a special-element block " *
                    "(MSH type $(block.msh)) into the current mesh"))
                push!(out,ElementBlock(block.msh,block.nodes.+shift,
                                       block.tags))
            end
        end
        return out
    end
    blocks=append!(shifted_blocks(a),shifted_blocks(b))
    names=merge(a isa MixedMesh ? a.physical_names :
                    Dict{Tuple{Int,Int},String}(),
                b isa MixedMesh ? b.physical_names :
                    Dict{Tuple{Int,Int},String}())
    return MixedMesh(coords,blocks;physical_names=names)
end

# `Merge "x.msh"` — fold the line/triangle/tetrahedron families into the
# mid-file simplex mesh. Higher-order cells collapse to their corner vertices
# (matching the `read_msh` fold); other families and point cells are accepted
# by `GModel::readMSH` upstream but have no place in the simplex buffer — they
# still land on the discrete entity records in `_geo_merge_discrete!`.
function _geo_mixed_simplex_mesh(mixed::MixedMesh,caller::AbstractString)
    ns=0;nf=0;nt=0
    for block in mixed.blocks
        block isa SpecialElementBlock && continue
        family=msh_family(block.msh)
        # Constant/order-zero cells have one node and no simplex corners.
        msh_spec(block.msh).order==0 && continue
        if family===:lin
            ns=Base.checked_add(ns,size(block.nodes,2))
        elseif family===:tri
            nf=Base.checked_add(nf,size(block.nodes,2))
        elseif family===:tet
            nt=Base.checked_add(nt,size(block.nodes,2))
        end
    end
    segs=Matrix{Int32}(undef,2,ns);st=Vector{Int32}(undef,ns)
    tris=Matrix{Int32}(undef,3,nf);tt=Vector{Int32}(undef,nf)
    tets=Matrix{Int32}(undef,4,nt);qt=Vector{Int32}(undef,nt)
    js=0;jf=0;jt=0
    @inbounds for block in mixed.blocks
        block isa SpecialElementBlock && continue
        family=msh_family(block.msh)
        # Constant/order-zero cells have one node and no simplex corners.
        msh_spec(block.msh).order==0 && continue;n=size(block.nodes,2);n==0 && continue
        if family===:lin
            segs[:,js+1:js+n].=view(block.nodes,1:2,:)
            copyto!(st,js+1,block.tags,1,n);js+=n
        elseif family===:tri
            tris[:,jf+1:jf+n].=view(block.nodes,1:3,:)
            copyto!(tt,jf+1,block.tags,1,n);jf+=n
        elseif family===:tet
            tets[:,jt+1:jt+n].=view(block.nodes,1:4,:)
            copyto!(qt,jt+1,block.tags,1,n);jt+=n
        end
    end
    mesh=Mesh(mixed.coords;segs=segs,tris=tris,tets=tets,
        seg_tag=st,tri_tag=tt,tet_tag=qt)
    # `read_msh` validates the folded mesh; keep the same gate on this path.
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: merged mesh is invalid — " * join(diagnostic.messages,"; ")))
    return mesh
end

# `GModel::readMSH` — each elementary (MSH2) or entity (MSH4) tag becomes a
# discrete model entity holding its cells; declared entity metadata adds
# boundary lists, physical memberships, and physical names. MSH2 supplies no
# node/element tag space or node classification, so nodes take the compact
# 1:N tags and each vertex lands on the first entity whose cells reference it
# (upstream's first-claim vertex map); MSH4 keeps the file's node/element tags,
# per-node entity classification, and parametric coordinates.
# Batch the import's node append so large UID maps do not copy the entire
# coordinate/parameter prefix once per vertex. Actual row order, IEEE bits and
# the existing duplicate-label append policy are retained.
function _geo_append_discrete_import_nodes!(record,mixed,positions,data)
    old_count=length(record.node_tags)
    count=length(positions)
    width=size(record.node_params,1)
    for position in positions
        params=data===nothing ? nothing : data.node_parametric[position]
        params===nothing || (width=max(width,length(params)))
    end
    coordinates=Matrix{Float64}(undef,3,old_count+count)
    coordinates[:,1:old_count].=record.node_coords
    parameters=width==0 ? record.node_params : zeros(Float64,width,old_count+count)
    width==0 || (parameters[1:size(record.node_params,1),1:old_count].=record.node_params)
    sizehint!(record.node_tags,old_count+count)
    for (offset,position) in enumerate(positions)
        tag=data===nothing ? UInt64(position) : data.external_node_tags[position]
        params=data===nothing ? nothing : data.node_parametric[position]
        push!(record.node_tags,tag)
        coordinates[:,old_count+offset].=view(mixed.coords,:,position)
        values=params===nothing ? Float64[] : copy(params)
        record.node_native_params[tag]=values
        width==0 || isempty(values) ||
            (parameters[1:length(values),old_count+offset].=values)
    end
    record.node_coords=coordinates
    record.node_params=parameters
    return nothing
end

function _geo_merge_discrete!(m::GeoModel,mixed::MixedMesh,
                              caller::AbstractString)
    blocks=mixed.blocks
    data=mixed.entity_data
    elementary=mixed.elementary_entities
    nnode=size(mixed.coords,2)
    # Per-cell owner entity; entity order follows first appearance upstream.
    order=Tuple{Int,Int}[]
    cells=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (bi,block) in enumerate(blocks)
        block isa SpecialElementBlock && continue
        dim=msh_dimension(block.msh)
        n=size(block.nodes,2)
        for j in 1:n
            tag=Int(data!==nothing ? data.block_entities[bi][j] :
                    elementary!==nothing ? elementary[bi][j] : Int32(0))
            key=(dim,tag)
            if !haskey(cells,key)
                cells[key]=NTuple{2,Int}[]
                push!(order,key)
            end
            push!(cells[key],(bi,j))
        end
    end
    entity_records=data===nothing ? Dict{Tuple{Int,Int},Elements.MixedEntity}() :
        mixed.partition_data===nothing ? data.entities : merge(data.entities,mixed.partition_data.entities)
    data===nothing || for key in keys(entity_records)
        haskey(cells,key) && continue
        cells[key]=NTuple{2,Int}[]
        push!(order,key)
    end
    # Node ownership: the MSH4 node-block classification is authoritative;
    # remaining vertices are claimed by the first entity whose cells use them.
    owner=fill((-1,Int32(0)),nnode)
    data===nothing ||
        for i in 1:nnode
            owner[i]=data.node_entities[i]
            key=(owner[i][1],Int(owner[i][2]))
            if !haskey(cells,key)
                cells[key]=NTuple{2,Int}[];push!(order,key)
            end
        end
    for (bi,block) in enumerate(blocks)
        block isa SpecialElementBlock && continue
        dim=msh_dimension(block.msh)
        for j in 1:size(block.nodes,2)
            tag=Int32(data!==nothing ? data.block_entities[bi][j] :
                      elementary!==nothing ? elementary[bi][j] : Int32(0))
            for ni in view(block.nodes,:,j)
                owner[ni][1]<0 && (owner[ni]=(dim,tag))
            end
        end
    end
    isempty(order) && return nothing
    # Records: existing discrete entities append; a collision with a real
    # model entity leaves its cells mesh-level only (the entity itself is not
    # reclassified); a zero tag is stored literally like upstream.
    records=Dict{Tuple{Int,Int},Union{DiscreteEntity,Nothing}}()
    for key in order
        dim,tag=key
        records[key]=if haskey(m.discrete,key)
            m.discrete[key]
        elseif _has_entity(m,dim,tag)
            nothing
        elseif tag==0
            m.discrete[key]=DiscreteEntity()
        else
            m.discrete[(dim,add_discrete_entity!(m,dim,tag))]
        end
    end
    node_tag=i->data!==nothing ? data.external_node_tags[i] : UInt64(i)
    materialized=Dict{Int,Tuple{Int,Int}}()
    planned_nodes=Dict{Tuple{Int,Int},Vector{Int32}}()
    function materialize(i::Int,key::Tuple{Int,Int})
        haskey(materialized,i) && return nothing
        record=get(records,key,nothing)
        record===nothing && return nothing
        push!(get!(Vector{Int32},planned_nodes,key),Int32(i))
        materialized[i]=key
        return nothing
    end
    # Explicitly classified nodes materialize on their owner first.
    data===nothing || for i in 1:nnode
        owner[i][1]<0 && continue
        materialize(i,(owner[i][1],Int(owner[i][2])))
    end
    next_elem=nothing
    used=Set{UInt64}()
    for (_,record) in _discrete_mesh_records_model(m)
        union!(used,record.element_tags)
    end
    for key in order
        record=records[key]
        record===nothing && continue
        for (bi,j) in cells[key]
            block=blocks[bi]
            nodes=block.nodes[:,j]
            for ni in nodes
                materialize(Int(ni),key)
            end
            # MSH4 element tags are file tags; MSH2 and colliding tags fall
            # back to fresh allocations.
            element_tag=data!==nothing ? data.external_element_tags[bi][j] : nothing
            if element_tag===nothing || element_tag in used
                next_elem===nothing && (next_elem=_model_fresh_element_tag(m))
                while next_elem in used
                    next_elem+=UInt64(1)
                end
                element_tag=next_elem;next_elem+=UInt64(1)
            end
            push!(used,element_tag)
            _record_append_element!(record,block.msh,element_tag,
                UInt64[node_tag(Int(ni)) for ni in nodes])
        end
    end
    for key in order
        positions=get(planned_nodes,key,nothing)
        positions===nothing && continue
        _geo_append_discrete_import_nodes!(records[key],mixed,positions,data)
    end
    # Physical memberships live on the entity (upstream `GEntity::physicals`):
    # MSH4 declares them per entity, MSH2 derives them from the cell tags.
    for key in order
        dim,tag=key
        # `GEntity::physicals` keeps declaration order — MSH4's per-entity
        # `physicalTags` list order or MSH2's first-seen cell-tag order —
        # so collect into a deduplicating vector, not a Set whose hash
        # order scrambles e.g. [5,6] into [6,5].
        physicals=Int32[]
        seen_physicals=Set{Int32}()
        if data!==nothing
            entity=get(entity_records,key,nothing)
            entity===nothing || append!(physicals,entity.physical_tags)
        else
            for (bi,j) in cells[key]
                cell_tag=blocks[bi].tags[j]
                cell_tag>0 && !(cell_tag in seen_physicals) &&
                    (push!(physicals,cell_tag);push!(seen_physicals,cell_tag))
            end
        end
        for ptag in physicals
            ep=get!(() -> Int[],m.entity_physicals,key)
            member=tag>0 ? Int(ptag) : -Int(ptag)
            member in ep || push!(ep,member)
            members=get!(() -> Int[],m.physical,(dim,Int(ptag)))
            if !(tag in members)
                push!(members,tag);sort!(members)
            end
            m.physical_tag_max=max(m.physical_tag_max,Int(ptag))
        end
    end
    # Declared boundary lists attach to the new records only.
    data===nothing || for key in order
        record=records[key]
        record===nothing && continue
        entity=get(entity_records,key,nothing)
        entity===nothing && continue
        dim=key[1];dim==0 && continue
        record.boundary=NTuple{2,Int}[(dim-1,Int(b)) for b in entity.boundaries]
    end
    for (key,name) in mixed.physical_names
        get!(m.physical_names,key,name)
    end
    return nothing
end

# `Save`/`Print` — `CreateOutputFile` writes the current mesh. `.msh` keeps
# its historical MSH 2.2 output; the versioned aliases `.msh2`/`.msh22` and
# `.msh4`/`.msh41` select that dialect explicitly, while `.msh1`/`.msh3`/
# `.msh40` reject (MSH 1.x/3.x never shipped writable support and 4.0 is not
# a stable public dialect). `.stl` writes ASCII surface facets through
# `IO.write_stl`.
function _geo_exec_save!(m::GeoModel,path::AbstractString,
                         context::_GeoNumericContext,caller::AbstractString;
                         kind::String="Save")
    file=_geo_fix_relative_path(context.file_name,path)
    extension=lowercase(splitext(file)[2])
    extension in (".msh1",".msh3",".msh40") && throw(ArgumentError(
        "$caller: unsupported MSH extension `$extension` (supported aliases: " *
        ".msh, .msh2/.msh22 for MSH 2.2, .msh4/.msh41 for MSH 4.1)"))
    version=if extension in (".msh",".msh2",".msh22")
        2.2
    elseif extension in (".msh4",".msh41")
        4.1
    else
        -1.0
    end
    (version>0 || extension==".stl") || throw(ArgumentError(
        "$caller: supported output extensions are .msh, .msh2/.msh22, " *
        ".msh4/.msh41 and .stl (got '$file')"))
    # `CreateOutputFile` writes the synced model — materialize the physical
    # view (and pending unknown-member warnings) first, when internals changed.
    _geo_sync_physical_view_if_changed!(m,context)
    mesh=context.mesh
    authority=m.meshing.mesh_authority
    if authority!==nothing && authority.mesh!==nothing && authority.carrier===mesh
        mesh=authority.mesh
    end
    if mesh===nothing
        isempty(m.points) || throw(ArgumentError(
            "$caller: no mesh exists — run `Mesh n;` first"))
        mesh=Mesh(zeros(3,0))
    end
    if extension==".stl"
        io_mod=getfield(parentmodule(@__MODULE__),:IO)
        io_mod.write_stl(file,mesh)
        return nothing
    end
    # `.msh` keeps MSH 2.2 for backward compatibility, but an explicit
    # `Mesh.MshFileVersion` wins; partition and ghost metadata cannot be
    # encoded before MSH 4.1, so a partitioned mesh always upgrades rather
    # than silently dropping or rejecting its own metadata.
    version==2.2 && haskey(context.option_numbers,("Mesh",0,"MshFileVersion")) &&
        (version=Float64(_geo_option_number(context,"Mesh",0,"MshFileVersion")))
    version<4.1 && mesh isa MixedMesh && mesh.partition_data!==nothing &&
        (version=4.1)
    if mesh isa MixedMesh && mesh.entity_data!==nothing
        # The actual product already owns its Point cells and public labels.
        # Serialize it directly instead of appending a second Point block or
        # replacing its node/element identity tables with dense proxy labels.
        write_mixed_msh(file,mesh;version)
        return nothing
    end
    # Upstream serializes an `MPoint` element per model vertex — the dim-0
    # owners the merge recorded carry them, emitted in vertex-tag order.
    point_elements=Tuple{Int,Int,Int}[]
    owner=context.mesh_node_owner
    dim0=sort!([(tag,node) for (node,(dim,tag)) in owner if dim==0])
    for (tag,node) in dim0
        physical=_model_projection_legacy_tag(
            _model_projection_physical_tags(m,0,tag))
        push!(point_elements,(node,Int(physical),tag))
    end
    if mesh isa MixedMesh
        # `write_msh` is simplex-only — the recombined product serializes
        # through `write_mixed_msh`, keeping each entity's blocks and its
        # elementary tags. MPoint cells join as a trailing point block.
        blocks=Vector{MixedElementBlock}(mesh.blocks)
        elementary=mesh.elementary_entities===nothing ? nothing :
            Vector{Int32}[copy(entities) for entities in
                          mesh.elementary_entities]
        if !isempty(point_elements)
            npoints=length(point_elements)
            pnodes=Matrix{Int32}(undef,1,npoints)
            ptags=Vector{Int32}(undef,npoints)
            pent=Vector{Int32}(undef,npoints)
            for (i,(node,physical,etag)) in enumerate(point_elements)
                pnodes[1,i]=Int32(node)
                ptags[i]=Int32(physical)
                pent[i]=Int32(etag)
            end
            push!(blocks,ElementBlock(15,pnodes,ptags))
            elementary!==nothing && push!(elementary,pent)
        end
        out=MixedMesh(mesh.coords,blocks;
                      physical_names=m.physical_names,
                      elementary_entities=elementary)
        write_mixed_msh(file,out;version)
    else
        write_msh(file,mesh;version,physical_names=m.physical_names,
                  point_elements=point_elements)
    end
    return nothing
end

function _geo_system_call(command::AbstractString,caller::AbstractString;
                          wait::Bool=true)
    isempty(strip(command)) && return nothing
    cmd=Sys.iswindows() ? `cmd /c $command` : `sh -c $command`
    wait ? run(cmd) : run(cmd;wait=false)
    return nothing
end

include("GeoMeshOrientation.jl")
include("GeoHomologyNative.jl")

end # module
