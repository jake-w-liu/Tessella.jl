struct _NativeQuadSplitTrace
    parent_nodes::NTuple{4,UInt64}
    child_tags::NTuple{2,UInt64}
    child_nodes::NTuple{2,NTuple{3,UInt64}}
end

# Bind carrier columns only through actual public IDs or an unambiguous set of
# primary-node public IDs. Geometry never establishes node or cell identity.
function _native_quad_split_carrier_bindings(previous,caller)
    carrier=previous.carrier
    actual=previous.mesh
    actual!==nothing && actual.entity_data!==nothing || throw(ArgumentError(
        "$caller: retained carrier samples lack actual public node identities"))
    length(previous.carrier_rows)==nnodes(carrier) &&
        all(row->1<=row<=nnodes(actual),previous.carrier_rows) &&
        allunique(previous.carrier_rows) || throw(ArgumentError(
        "$caller: retained carrier node bindings are incomplete or ambiguous"))
    node_tags=actual.entity_data.external_node_tags[previous.carrier_rows]
    catalog=Dict{UInt64,Tuple{Symbol,Tuple}}()
    primary_cells=Dict{Tuple{Symbol,Tuple},Vector{UInt64}}()
    for (bi,block) in enumerate(actual.blocks),column in axes(block.nodes,2)
        family=msh_family(block.msh)
        primary=Elements.msh_num_nodes(msh_type(family,1))
        key=Tuple(sort(actual.entity_data.external_node_tags[block.nodes[1:primary,column]]))
        tag=actual.entity_data.external_element_tags[bi][column]
        catalog[tag]=(family,key)
        push!(get!(()->UInt64[],primary_cells,(family,key)),tag)
    end
    blocks=Model._model_mesh_carrier_blocks(carrier)
    supplied=previous.carrier_cells
    isempty(supplied) || length(supplied)==length(blocks) || throw(ArgumentError(
        "$caller: retained carrier cell bindings do not align with its blocks"))
    data=carrier isa MixedMesh ? carrier.entity_data : nothing
    resolved=Vector{Vector{UInt64}}()
    for (bi,(msh,cells)) in enumerate(blocks)
        family=msh_family(msh)
        primary=Elements.msh_num_nodes(msh_type(family,1))
        isempty(supplied) || length(supplied[bi])==size(cells,2) || throw(ArgumentError(
            "$caller: retained carrier cell bindings do not align with its columns"))
        tags=UInt64[]
        for column in axes(cells,2)
            key=(family,Tuple(sort(node_tags[cells[1:primary,column]])))
            tag=isempty(supplied) ? nothing : supplied[bi][column]
            if tag!==nothing
                get(catalog,tag,nothing)==key || throw(ArgumentError(
                    "$caller: retained carrier cell $tag disagrees with its primary-node public IDs"))
            elseif data!==nothing
                proposed=data.external_element_tags[bi][column]
                get(catalog,proposed,nothing)==key && (tag=proposed)
            end
            if tag===nothing
                candidates=get(primary_cells,key,UInt64[])
                length(candidates)==1 || throw(ArgumentError(
                    "$caller: retained carrier cell lacks an unambiguous public-ID binding"))
                tag=only(candidates)
            end
            push!(tags,tag)
        end
        push!(resolved,tags)
    end
    return node_tags,resolved
end

function _native_quad_split_carrier(plan,previous,splits,caller)
    if previous===nothing || previous.carrier===nothing
        return plan.mesh,collect(Int32,1:nnodes(plan.mesh)),
            plan.mesh.entity_data.external_element_tags
    end
    carrier=previous.carrier
    node_tags,bindings=_native_quad_split_carrier_bindings(previous,caller)
    # `rebound_previous` keeps the previous mesh and products, so retained
    # Surface constructor inputs stay bound to the same products.
    rebound_previous=Model._ModelMeshAuthority(previous.mesh,previous.products,
        previous.done,previous.private_rank,previous.native_inputs,carrier,
        previous.carrier_rows,previous.producer_rows,bindings,
        previous.surface_constructor_inputs)
    affected=any(tag->haskey(splits,tag),Iterators.flatten(bindings))
    if !affected
        rebound,rows=Model._model_mesh_rebound_carrier(rebound_previous,plan.mesh;
            node_tag_map=plan.node_tag_map,element_tag_map=plan.element_tag_map)
        rebound!==nothing || throw(ArgumentError(
            "$caller: unchanged carrier samples cannot be rebound by public IDs"))
        return rebound,rows,Model._model_mesh_rebound_cell_tags(rebound_previous,plan.element_tag_map)
    end
    carrier isa MixedMesh || throw(ArgumentError(
        "$caller: a changed Quad carrier requires explicit Mixed Quad connectivity"))
    positions=Dict(tag=>Int32(row) for (row,tag) in
        enumerate(plan.mesh.entity_data.external_node_tags))
    rows=Int32[get(positions,get(plan.node_tag_map,tag,nothing),Int32(0)) for tag in node_tags]
    all(>(0),rows) && allunique(rows) || throw(ArgumentError(
        "$caller: Quad split loses retained carrier node identities"))
    carrier_positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(node_tags))
    blocks=ElementBlock[];output_bindings=Vector{Union{Nothing,UInt64}}[]
    owners=Vector{Int32}[]
    old_data=carrier.entity_data
    elementary=carrier.elementary_entities
    have_owners=old_data!==nothing || elementary!==nothing
    for (bi,block) in enumerate(carrier.blocks)
        kept=Int[];child_nodes=Vector{Int32}[];child_tags=UInt64[]
        child_physical=Int32[];child_owners=Int32[]
        source_owners=old_data!==nothing ? old_data.block_entities[bi] :
            elementary===nothing ? nothing : elementary[bi]
        for column in axes(block.nodes,2)
            parent=bindings[bi][column]
            trace=get(splits,parent,nothing)
            if trace===nothing
                push!(kept,column)
                continue
            end
            msh_family(block.msh)===:qua || throw(ArgumentError(
                "$caller: changed parent $parent has no explicit carrier Quad column"))
            Tuple(sort(node_tags[block.nodes[1:4,column]]))==Tuple(sort(collect(trace.parent_nodes))) ||
                throw(ArgumentError("$caller: changed carrier Quad $parent loses its primary-node public IDs"))
            for (tag,nodes) in zip(trace.child_tags,trace.child_nodes)
                mapped=Int32[get(carrier_positions,node,Int32(0)) for node in nodes]
                all(>(0),mapped) || throw(ArgumentError(
                    "$caller: child Triangle $tag lacks a retained carrier corner sample"))
                push!(child_nodes,mapped);push!(child_tags,tag)
                push!(child_physical,block.tags[column])
                source_owners===nothing || push!(child_owners,source_owners[column])
            end
        end
        if !isempty(kept)
            push!(blocks,ElementBlock(block.msh,block.nodes[:,kept],block.tags[kept]))
            retained=Union{Nothing,UInt64}[get(plan.element_tag_map,bindings[bi][column],nothing) for column in kept]
            all(!isnothing,retained) || throw(ArgumentError(
                "$caller: unchanged carrier cells lose their actual identities"))
            push!(output_bindings,retained)
            have_owners && push!(owners,source_owners[kept])
        end
        if !isempty(child_nodes)
            push!(blocks,ElementBlock(2,hcat(child_nodes...),child_physical))
            push!(output_bindings,Union{Nothing,UInt64}[tag for tag in child_tags])
            have_owners && push!(owners,child_owners)
        end
    end
    data=old_data===nothing ? nothing : Elements._native_mixed_entity_data(old_data.entities;
        node_entities=old_data.node_entities,node_parametric=old_data.node_parametric,
        external_node_tags=plan.mesh.entity_data.external_node_tags[rows],
        block_entities=owners,external_element_tags=[UInt64[tag for tag in tags] for tags in output_bindings])
    replacement=MixedMesh(carrier.coords,blocks;physical_names=carrier.physical_names,
        entity_data=data,elementary_entities=have_owners ? owners : nothing,
        periodic_links=carrier.periodic_links,ancillary_sections=carrier.ancillary_sections)
    if !isempty(carrier.data_sections) || carrier.partition_data!==nothing
        old_data!==nothing || throw(ArgumentError(
            "$caller: retained carrier file rows require explicit public-ID metadata"))
        bound_data=Elements._native_mixed_entity_data(old_data.entities;
            node_entities=old_data.node_entities,node_parametric=old_data.node_parametric,
            external_node_tags=node_tags,block_entities=old_data.block_entities,
            external_element_tags=bindings)
        bound_source=MixedMesh(carrier.coords,carrier.blocks;
            physical_names=carrier.physical_names,entity_data=bound_data,
            elementary_entities=carrier.elementary_entities,periodic_links=carrier.periodic_links,
            ancillary_sections=carrier.ancillary_sections,data_sections=carrier.data_sections,
            partition_data=carrier.partition_data)
        sections=_reclassify_data_sections(bound_source,replacement,plan,caller)
        partition=_reclassify_partition_data(bound_source,replacement,plan,caller)
        replacement=MixedMesh(replacement.coords,replacement.blocks;
            physical_names=replacement.physical_names,entity_data=replacement.entity_data,
            elementary_entities=replacement.elementary_entities,periodic_links=replacement.periodic_links,
            ancillary_sections=replacement.ancillary_sections,data_sections=sections,partition_data=partition)
    end
    return replacement,rows,output_bindings
end

# Native quadsToTriangles consumes all four temporary triangle identities.
# Only corner vertices enter the replacements; interpolation nodes stay owned.
function _native_quad_split_diagonal(points)
    a,b,c,d=points
    first=min(_native_cell_triangle_gamma(a,b,c),_native_cell_triangle_gamma(c,d,a))
    second=min(_native_cell_triangle_gamma(b,c,d),_native_cell_triangle_gamma(d,a,b))
    normal_first=_native_cell_unit(_native_cell_cross(
        _native_cell_difference(c,b),_native_cell_difference(d,b)))
    normal_second=_native_cell_unit(_native_cell_cross(
        _native_cell_difference(a,d),_native_cell_difference(b,d)))
    return first>second || _native_cell_dot(normal_first,normal_second)<0.0
end

function _native_quad_split_tags(maximum_element,points,caller)
    # Default-argument MTriangle constructors consume the full size_t clock.
    # All four candidates are born, including the two discarded triangles.
    diagonal=_native_quad_split_diagonal(points)
    first=maximum_element+UInt64(diagonal ? 1 : 3)
    return first,first+UInt64(1),maximum_element+UInt64(4),diagonal
end

function _native_quad_split_node_point(index,cached,public,overlay,tag,caller)
    if cached!==nothing
        count=overlay===nothing ? nnodes(cached) : nnodes(overlay)
        cache_row=public===nothing ? (1<=tag<=count ? Int(tag) : 0) :
            get(public.node_indices,UInt64(tag),0)
        if cache_row!=0
            geometry=cache_row<=nnodes(cached) ? cached : overlay
            geometry===nothing || return ntuple(axis->geometry.coords[axis,cache_row],3)
        end
    end
    binding=get(index,tag,nothing)
    binding===nothing && throw(ArgumentError("$caller: node $tag has no coordinates"))
    record,record_row=binding
    return ntuple(axis->record.node_coords[axis,record_row],3)
end

function _native_quad_split_record_plan!(model,cached,keys,limit,maximum_element,caller)
    changed=false
    # The node index is operation-owned and resolved before any publication.
    index=_record_node_index(model)
    public=cached===nothing ? nothing : _mesh_public_tags(cached)
    overlay=cached===nothing ? nothing : _high_order_overlay(cached)
    deleted=UInt64[]
    counter_wrapped=ELEMENT_COUNTER_WRAPPED[]
    for key in keys
        record=_entity_record(model,key[1],key[2])
        record===nothing && continue
        removed=falses(length(record.element_types))
        additions=Tuple{UInt64,Vector{UInt64}}[]
        for position in eachindex(record.element_types)
            msh_family(record.element_types[position])===:qua || continue
            nodes=record.element_nodes[position]
            length(nodes)==msh_spec(record.element_types[position]).nnodes || throw(ArgumentError(
                "$caller: malformed quadrangle connectivity"))
            points=ntuple(4) do corner
                _native_quad_split_node_point(index,cached,public,overlay,nodes[corner],caller)
            end
            MeshElementQuality._quadrangle_native_eta(points)<limit+1e-12 || continue
            previous=maximum_element
            first,second,maximum_element,diagonal=_native_quad_split_tags(previous,points,caller)
            counter_wrapped |= maximum_element<previous
            if diagonal
                push!(additions,(first,UInt64[nodes[1],nodes[2],nodes[3]]))
                push!(additions,(second,UInt64[nodes[3],nodes[4],nodes[1]]))
            else
                push!(additions,(first,UInt64[nodes[2],nodes[3],nodes[4]]))
                push!(additions,(second,UInt64[nodes[4],nodes[1],nodes[2]]))
            end
            removed[position]=true
            push!(deleted,UInt64(record.element_tags[position]))
        end
        any(removed) || continue
        changed=true
        kept=findall(!,removed)
        record.element_types=record.element_types[kept]
        record.element_tags=record.element_tags[kept]
        record.element_nodes=record.element_nodes[kept]
        for (tag,nodes) in additions
            _record_append_element!(record,2,tag,nodes)
        end
    end
    return changed,maximum_element,deleted,counter_wrapped
end

# Each closure receives one operation-owned cell through an immutable
# function argument, rather than capturing locals rebound by the plan loop.
function _native_quad_split_trace(assembly,cell,first_nodes,second_nodes,first,second)
    return _NativeQuadSplitTrace(
        ntuple(corner->assembly.nodes[cell.nodes[corner]].tag,4),(first,second),
        (ntuple(corner->assembly.nodes[first_nodes[corner]].tag,3),
         ntuple(corner->assembly.nodes[second_nodes[corner]].tag,3)))
end

function _native_quad_split_mixed_plan(model,cached,class,keys,limit,maximum_element,caller)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],maximum_element,typemax(Int32),typemax(Int32),String(caller),
        NODE_COUNTER_WRAPPED[],ELEMENT_COUNTER_WRAPPED[])
    source_positions=_dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    selected=Set((dim,Int32(tag)) for (dim,tag) in keys)
    output=_Dim01Cell[]
    sizehint!(output,length(assembly.cells))
    changed=false
    splits=Dict{UInt64,_NativeQuadSplitTrace}()
    # The entity order and stable source order match GFace's quadrangle vector.
    order=sortperm(assembly.cells;by=cell->cell.owner,alg=Base.Sort.MergeSort)
    for position in order
        cell=assembly.cells[position]
        if !(cell.owner in selected) || msh_family(cell.msh)!==:qua
            push!(output,cell)
            continue
        end
        points=ntuple(corner->assembly.nodes[cell.nodes[corner]].point,4)
        if !(MeshElementQuality._quadrangle_native_eta(points)<limit+1e-12)
            push!(output,cell)
            continue
        end
        first,second,next,diagonal=_native_quad_split_tags(assembly.cell_max,points,caller)
        assembly.cell_counter_wrapped |= next<assembly.cell_max
        assembly.cell_max=next
        nodes=cell.nodes
        if diagonal
            first_nodes=Int32[nodes[1],nodes[2],nodes[3]]
            second_nodes=Int32[nodes[3],nodes[4],nodes[1]]
        else
            first_nodes=Int32[nodes[2],nodes[3],nodes[4]]
            second_nodes=Int32[nodes[4],nodes[1],nodes[2]]
        end
        splits[cell.tag]=_native_quad_split_trace(
            assembly,cell,first_nodes,second_nodes,first,second)
        push!(output,_Dim01Cell(2,cell.owner,first_nodes,first,nothing,cell.physical))
        push!(output,_Dim01Cell(2,cell.owner,second_nodes,second,nothing,cell.physical))
        changed=true
    end
    changed || return nothing
    assembly.cells=output
    plan=_dim01_finish!(assembly,model,model_entities(model),Set{Tuple{Int,Int32}}(),
        falses(length(assembly.nodes)),false,false;source_mesh=cached,
        source_node_positions=source_positions)
    previous=model.meshing.mesh_authority
    # Reconstruct retained file rows through public identities. Deleted Quad
    # rows leave these mesh-attached sections; session views retain their UIDs.
    if !isempty(cached.ancillary_sections) || !isempty(cached.data_sections) || cached.partition_data!==nothing
        cached.entity_data!==nothing ||
            (isempty(cached.data_sections) && cached.partition_data===nothing) ||
            throw(ArgumentError("$caller: retained file rows require explicit public-ID metadata"))
        sections=_reclassify_data_sections(cached,plan.mesh,plan,caller)
        partition=_reclassify_partition_data(cached,plan.mesh,plan,caller)
        replacement=MixedMesh(plan.mesh.coords,plan.mesh.blocks;
            physical_names=plan.mesh.physical_names,entity_data=plan.mesh.entity_data,
            elementary_entities=plan.mesh.elementary_entities,periodic_links=plan.mesh.periodic_links,
            ancillary_sections=cached.ancillary_sections,data_sections=sections,partition_data=partition)
        replacement_class=_higher_classification_from_data(replacement,caller)
        plan=merge(plan,(mesh=replacement,class=replacement_class))
    end
    carrier,carrier_rows,carrier_cells=_native_quad_split_carrier(plan,previous,splits,caller)
    plan=merge(plan,(quad_split_trace=splits,))
    return _higher_rebind_plan_authority!(plan,previous;carrier,carrier_rows,carrier_cells)
end

function _native_split_quadrangles(limit::Float64,tag_value,caller)
    original=_model_locked()
    cached=LAST_MESH[]
    class=cached===nothing ? nothing : _cached_classification_locked(cached)
    if cached isa MixedMesh && class===nothing && cached.entity_data!==nothing
        class=_higher_classification_from_data(cached,caller)
    end
    keys=Tuple{Int,Int}[]
    if tag_value<0
        append!(keys,[(2,Int(tag)) for (dimension,tag) in model_entities(original) if dimension==2])
    else
        tag_value<=typemax(Int32) || throw(ArgumentError("$caller: tag is out of range"))
        _model_entity_known(original,2,Int(tag_value)) || throw(ArgumentError(
            "$caller: unknown entity (2, $tag_value)"))
        push!(keys,(2,Int(tag_value)))
    end
    sort!(keys)
    isempty(keys) && return nothing
    maximum_element=_get_max_element_tag()
    staged=Model._model_mesh_transaction_copy(original)
    if cached isa MixedMesh
        class===nothing && throw(ArgumentError("$caller: existing mesh needs entity classification"))
        plan=_native_quad_split_mixed_plan(staged,cached,class,keys,limit,maximum_element,caller)
        plan===nothing && return nothing
        _commit_tagged_plan_locked!(plan;preserve_visibility=true,element_tag_map=plan.visibility_tag_map)
    else
        changed,maximum_element,deleted,counter_wrapped=_native_quad_split_record_plan!(
            staged,cached,keys,limit,maximum_element,caller)
        changed || return nothing
        # Existing simplex coordinates and any P2 overlay retain their identity.
        _replace_mesh_cache_with_overlay_locked!(cached,class,caller;
            overlay=cached===nothing ? nothing : _high_order_overlay(cached),certify=false,
            preserve_visibility=true,staged_model=staged)
        if counter_wrapped
            ELEMENT_TAG_MAX[]=maximum_element
            ELEMENT_COUNTER_WRAPPED[]=true
        else
            ELEMENT_TAG_MAX[]=max(ELEMENT_TAG_MAX[],maximum_element)
        end
        for tag in deleted
            delete!(ELEMENT_VISIBILITY[],tag)
        end
        index=_current_slot_index_locked()
        if index!=0
            for flags in values(MODEL_SLOTS[index].window_element_visibility),tag in deleted
                delete!(flags,tag)
            end
        end
    end
    return nothing
end
