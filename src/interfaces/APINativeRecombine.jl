# Detached public recombination planning. This adapter uses actual public
# identities and raw mesh-vertex parameters; strict MSH metadata keeps its
# separate owner-dimensional parameter projection. It never publishes state.
using ..Recombine: _primary_recombine_mixed, _PrimarySurfaceContext,
    _PrimaryPlaneEvaluation, _PrimaryBlockedSurfaceEvaluation,
    _primary_plane_parameters, _primary_bipartite_raw, _PRIMARY_BIPARTITE_EMPTY_ERROR

function _api_native_recombine_option_integer(options,name,default,caller)
    value=get(options,name,default)
    value isa Real && !(value isa Bool) && isfinite(value) || throw(ArgumentError(
        "$caller: $name must be finite for native recombination"))
    return try
        trunc(Int,value)
    catch error
        error isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name exceeds the supported integer range"))
    end
end

function _api_native_recombine_settings(options,caller)
    algorithm=_api_native_recombine_option_integer(options,"Mesh.RecombinationAlgorithm",1,caller)
    algorithm in (0,1,2,3,4) || throw(ArgumentError(
        "$caller: this public recombination selector is not implemented by this adapter"))
    algorithm==4 && return (algorithm=:bipartite,reposition=false,
        topology_passes=0,smoothing=0,minimum_quality=0.0)
    reposition=get(options,"Mesh.RecombineNodeRepositioning",1.0)
    reposition isa Real && !(reposition isa Bool) && isfinite(reposition) ||
        throw(ArgumentError("$caller: Mesh.RecombineNodeRepositioning must be finite"))
    quality=get(options,"Mesh.RecombineMinimumQuality",.01)
    quality isa Real && !(quality isa Bool) && isfinite(quality) ||
        throw(ArgumentError("$caller: Mesh.RecombineMinimumQuality must be finite"))
    return (algorithm=algorithm in (1,3) ? :blossom : :greedy,
        reposition=!iszero(reposition),
        topology_passes=_api_native_recombine_option_integer(options,"Mesh.RecombineOptimizeTopology",5,caller),
        smoothing=_api_native_recombine_option_integer(options,"Mesh.Smoothing",1,caller),
        minimum_quality=Float64(quality))
end

function _api_native_recombine_source_payload(cached,class,node_mapping,cell_mapping,caller)
    cached isa MixedMesh || return (periodic_links=Elements.MixedPeriodicLink[],
        ancillary_sections=Elements.MshAncillarySection[],data_sections=Elements.MshDataSection[],
        partition_data=nothing)
    data=cached.entity_data
    data!==nothing || throw(ArgumentError("$caller: a mixed cache needs entity metadata"))
    sections=Elements.MshDataSection[]
    for section in cached.data_sections
        push!(sections,Elements.MshDataSection(section.name,section.strings,section.reals,
            section.header,cell_mapping[section.elements],node_mapping[section.nodes],
            section.row_nodes,section.implicit_nodes,section.values,section.anchor))
    end
    links=Elements.MixedPeriodicLink[]
    for link in cached.periodic_links
        push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
            node_mapping[link.slave_nodes],node_mapping[link.master_nodes];affine=link.affine))
    end
    partition=cached.partition_data
    if partition!==nothing
        ghosts=Elements.MixedGhostElement[Elements.MixedGhostElement(cell_mapping[ghost.element],
            ghost.partition,ghost.ghost_partitions) for ghost in partition.ghost_elements]
        partition=Elements.MixedPartitionData(partition.num_partitions;entities=partition.entities,
            ghost_entities=partition.ghost_entities,ghost_elements=ghosts)
    end
    return (periodic_links=links,ancillary_sections=cached.ancillary_sections,
        data_sections=sections,partition_data=partition)
end

# Reuse the classified source-ingestion rules. They join exact public UIDs,
# reject conflicting coordinates/owners, and retain raw record parameters.
# Existing cache rows and cells precede only the additional record objects.
function _api_native_recombine_source(m,cached,class,caller;
        initial_max_node_tag,initial_max_element_tag,max_nodes,max_cells)
    cached===nothing || class!==nothing && class.public_tags!==nothing ||
        throw(ArgumentError("$caller: an existing cache needs classified public identities"))
    cached isa Mesh && _high_order_overlay(cached)!==nothing && throw(ArgumentError(
        "$caller: native recombination with a high-order overlay is not implemented"))
    if cached isa MixedMesh
        all(block->block isa ElementBlock,cached.blocks) || throw(ArgumentError(
            "$caller: native recombination requires fixed-connectivity cache blocks"))
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),Int(max_nodes),Int(max_cells),String(caller))
    cache_rows=_dim01_ingest!(assembly,m,cached,class;preserve_source_order=true)
    allunique(cache_rows) || throw(ArgumentError("$caller: source cache repeats a public node identity"))
    rows_seen=Set(cache_rows)
    row_order=vcat(Int.(cache_rows),Int[index for index in eachindex(assembly.nodes) if !(index in rows_seen)])
    row_inverse=zeros(Int32,length(row_order))
    for (new,old) in enumerate(row_order);row_inverse[old]=Int32(new);end
    cell_positions=Dict(cell.tag=>index for (index,cell) in enumerate(assembly.cells))
    length(cell_positions)==length(assembly.cells) || throw(ArgumentError(
        "$caller: source union repeats a public element identity"))
    cell_order=Int[]
    if cached!==nothing
        for (bi,(_,_,offset,cells,_)) in enumerate(_cache_catalog(cached,class))
            physical=_cache_native_blocks(cached)[bi][3]
            for column in axes(cells,2)
                tag=_cache_element_tag(class.public_tags,offset+column)
                position=get(cell_positions,tag,0)
                position>0 || throw(ArgumentError("$caller: cached element $tag lost its source identity"))
                # An authoritative cache can have a richer legacy scalar-tag
                # payload than its mirrored entity record. Preserve that field.
                assembly.cells[position].physical=physical[column]
                push!(cell_order,position)
            end
        end
    end
    cells_seen=Set(cell_order)
    append!(cell_order,Int[index for index in eachindex(assembly.cells) if !(index in cells_seen)])
    allunique(cell_order) || throw(ArgumentError("$caller: source cache repeats a public cell identity"))
    node_tags=UInt64[assembly.nodes[index].tag for index in row_order]
    all(!iszero,node_tags) || throw(ArgumentError(
        "$caller: native zero-node UID contexts are not implemented by this adapter"))
    all(cell->!iszero(cell.tag),assembly.cells) || throw(ArgumentError(
        "$caller: native zero-element UID contexts are not implemented by this adapter"))
    coordinates=Matrix{Float64}(undef,3,length(row_order))
    owners=Tuple{Int,Int32}[assembly.nodes[index].owner for index in row_order]
    runtime=Union{Nothing,Vector{Float64}}[deepcopy(assembly.nodes[index].parameter) for index in row_order]
    for (new,old) in enumerate(row_order),axis in 1:3
        coordinates[axis,new]=assembly.nodes[old].point[axis]
    end
    projected=Union{Nothing,Vector{Float64}}[
        Model._model_mesh_projected_parameter(owner,value) for (owner,value) in zip(owners,runtime)]
    entities,boundaries=_dim01_entities(m,model_entities(m),assembly,collect(eachindex(assembly.nodes)))
    names=copy(m.physical_names)
    if cached isa MixedMesh
        data=cached.entity_data
        data!==nothing || throw(ArgumentError("$caller: a mixed cache needs entity metadata"))
        merge!(entities,data.entities)
        merge!(names,cached.physical_names)
        for (old,new) in enumerate(row_inverse[cache_rows])
            projected[new]=deepcopy(data.node_parametric[old])
        end
    end
    groups=Dict{Int,Vector{Int}}();families=Int[]
    for index in cell_order
        family=assembly.cells[index].msh
        haskey(groups,family) || push!(families,family)
        push!(get!(()->Int[],groups,family),index)
    end
    blocks=ElementBlock[];block_owners=Vector{Int32}[];external=Vector{UInt64}[]
    output_cell_positions=Dict{UInt64,Int32}();flat=0
    for family in families
        selected=groups[family]
        nodes=Matrix{Int32}(undef,Elements.msh_num_nodes(family),length(selected))
        physical=Int32[];classified=Int32[];tags=UInt64[]
        for (column,index) in enumerate(selected)
            cell=assembly.cells[index]
            nodes[:,column]=row_inverse[cell.nodes]
            push!(physical,cell.physical);push!(classified,cell.owner[2]);push!(tags,cell.tag)
            flat=Base.checked_add(flat,1)
            flat<=typemax(Int32) || throw(ArgumentError("$caller: source cells exceed Int32 dense indexing"))
            output_cell_positions[cell.tag]=Int32(flat)
        end
        push!(blocks,ElementBlock(family,nodes,physical));push!(block_owners,classified);push!(external,tags)
    end
    old_cells=Int32[]
    if cached!==nothing
        for (_,_,offset,cells,_) in _cache_catalog(cached,class),column in axes(cells,2)
            push!(old_cells,output_cell_positions[_cache_element_tag(class.public_tags,offset+column)])
        end
    end
    payload=_api_native_recombine_source_payload(cached,class,row_inverse[cache_rows],old_cells,caller)
    metadata=Elements._native_mixed_entity_data(entities;node_entities=owners,node_parametric=projected,
        external_node_tags=node_tags,block_entities=block_owners,external_element_tags=external)
    mesh=MixedMesh(coordinates,blocks;entity_data=metadata,elementary_entities=block_owners,
        physical_names=names,periodic_links=payload.periodic_links,
        ancillary_sections=payload.ancillary_sections,data_sections=payload.data_sections,
        partition_data=payload.partition_data)
    return (mesh,runtime_parameters=runtime,boundaries,
        max_node_tag=assembly.node_max,max_element_tag=assembly.cell_max)
end

function _api_native_recombine_topology_admitted(mesh)
    data=mesh.entity_data
    for (bi,block) in enumerate(mesh.blocks)
        msh_dimension(block.msh)==2 || continue
        for column in axes(block.nodes,2),row in @view(block.nodes[:,column])
            owner=data.node_entities[row]
            owner[1]==2 && owner[2]!=data.block_entities[bi][column] && return false
        end
    end
    return true
end

function _api_native_recombine_primitive_surface(m,surface)
    for collection in (m.box_extents,m.cylinders,m.spheres,m.cones)
        for volume in keys(collection),shell in get(m.volumes,volume,Int[]),
            face in get(m.surface_loops,abs(shell),Int[])
            abs(face)==surface && return true
        end
    end
    return false
end

function _api_native_recombine_surface_context(m,source,owner,rows,settings,
        topology_admitted,boundary_layer_context,caller)
    data=source.mesh.entity_data
    owners=data.node_entities[rows];tags=data.external_node_tags[rows]
    parameters=zeros(2,length(rows))
    parameter_present=trues(length(rows));project=nothing;plane=nothing
    discrete=haskey(m.discrete,(2,Int(owner)))
    evaluate=if discrete
        nothing
    elseif model_entity_type(m,2,Int(owner))!="Plane" ||
           haskey(m.surface_geometry,Int(owner)) || _api_native_recombine_primitive_surface(m,Int(owner))
        _PrimaryBlockedSurfaceEvaluation("$caller: native relocation for this nonlinear or primitive Surface context is not implemented")
    else
        plane=Model._model_plane_frame(m,Int(owner),caller)
        _PrimaryPlaneEvaluation(plane.origin,plane.first_direction,plane.second_direction)
    end
    if settings.reposition && settings.smoothing>0 && !discrete && any(value->value[1]==2,owners)
        evaluate isa _PrimaryPlaneEvaluation || throw(ArgumentError(evaluate.message))
        boundary_layer_context && throw(ArgumentError(
            "$caller: native relocation with active boundary-layer context is not implemented"))
        missing_face_parameters=any(enumerate(rows)) do (local_row,row)
            owners[local_row][1]==2 || return false
            parameter=source.runtime_parameters[row]
            return parameter===nothing || isempty(parameter)
        end
        if missing_face_parameters
            lower,upper=Model._model_plane_parameter_bounds(plane,caller)
            frame=plane
            project=point->first(Model._model_closest_plane_native(frame,point,lower,upper))
        end
        for (local_row,row) in enumerate(rows)
            dimension,entity=owners[local_row]
            parameter=source.runtime_parameters[row]
            uv=if dimension==2
                entity==owner || throw(ArgumentError(
                    "$caller: native relocation across Surface parameter owners is not implemented"))
                if parameter===nothing || isempty(parameter)
                    parameter_present[local_row]=false
                    _primary_plane_parameters(evaluate,(source.mesh.coords[1,row],
                        source.mesh.coords[2,row],source.mesh.coords[3,row]))
                else
                    length(parameter)==2 && all(isfinite,parameter) || throw(ArgumentError(
                        "$caller: native relocation requires two finite stored Surface parameters"))
                    (parameter[1],parameter[2])
                end
            elseif dimension==0
                point=get(m.points,Int(entity),nothing)
                point!==nothing || throw(ArgumentError(
                    "$caller: native relocation with a non-CAD Point boundary is not implemented"))
                _primary_plane_parameters(evaluate,(source.mesh.coords[1,row],
                    source.mesh.coords[2,row],source.mesh.coords[3,row]))
            elseif dimension==1
                (parameter===nothing || isempty(parameter) ||
                 (length(parameter)==1 && all(isfinite,parameter))) ||
                    throw(ArgumentError("$caller: native relocation requires one finite stored Curve parameter"))
                haskey(m.curves,Int(entity)) && Model._curve_type(m,Int(entity))===:line &&
                    Model._occ_geometry(m,Int(entity))===nothing &&
                    isempty(get(m.curve_control_points,Int(entity),Int[])) || throw(ArgumentError(
                    "$caller: native relocation for this Curve evaluation context is not implemented"))
                first,last=m.curves[Int(entity)];a=m.points[first];b=m.points[last]
                # Base native MVertex::getParameter writes zero before returning
                # false. The relocation caller still uses Curve reparametrization
                # at that value; the missing raw parameter storage stays empty.
                u=parameter===nothing || isempty(parameter) ? 0.0 : parameter[1]
                point=ntuple(axis->a[axis]+u*(b[axis]-a[axis]),3)
                _primary_plane_parameters(evaluate,point)
            else
                throw(ArgumentError("$caller: native relocation with Volume-owned Surface support is not implemented"))
            end
            parameters[1,local_row]=uv[1];parameters[2,local_row]=uv[2]
        end
    end
    return _PrimarySurfaceContext(owner,owners,tags,parameters,evaluate;
        reposition=settings.reposition,topology_passes=settings.topology_passes,
        smoothing=settings.smoothing,topology_admitted,parameter_present,project)
end

function _api_native_recombine_surface_signature(mesh,owner)
    data=mesh.entity_data
    cells=Tuple{Int,UInt64,Vector{UInt64},Int32}[]
    rows=Set{Int32}()
    for (bi,block) in enumerate(mesh.blocks)
        msh_dimension(block.msh)==2 || continue
        for column in axes(block.nodes,2)
            data.block_entities[bi][column]==owner || continue
            nodes=@view block.nodes[:,column]
            push!(cells,(block.msh,data.external_element_tags[bi][column],
                copy(data.external_node_tags[nodes]),block.tags[column]))
            union!(rows,nodes)
        end
    end
    sort!(cells;by=cell->(cell[1],cell[2]))
    ordered=sort!(collect(rows);by=row->data.external_node_tags[row])
    return (cells,nodes=data.external_node_tags[ordered],coordinates=mesh.coords[:,ordered],
        parameters=data.node_parametric[ordered])
end

function _api_native_recombine_has_surfaces(m,cached)
    # Recombination pairs first-order Surface elements; a Surface whose stored
    # cells are all high-order has nothing to pair and keeps its mesh.
    cached isa Mesh && ntris(cached)>0 && return true
    cached isa MixedMesh && any(block->msh_dimension(block.msh)==2 &&
                                    msh_spec(block.msh).order==1 &&
                                    !isempty(block.nodes),
                               cached.blocks) && return true
    return any(pair->pair[1][1]==2 &&
                   any(type->msh_dimension(type)==2 && msh_spec(type).order==1,
                       pair[2].element_types),
               Model._discrete_mesh_records_model(m))
end

# GModel's public recombine entry numbers vertices and cells after the raw
# stage, including when a prior Msg error skipped that stage. Keep its born
# maxima and detached runtime parameters while composing the two label maps.
function _api_native_recombine_renumber_plan(plan,caller;save_all)
    nodes=_higher_geo_renumber_plan(plan.model,plan.mesh,plan.class,:geo_renumber_nodes,caller;
        initial_max_node_tag=plan.max_node_tag,initial_max_element_tag=plan.max_element_tag,save_all)
    cells=_higher_geo_renumber_plan(nodes.model,nodes.mesh,nodes.class,:geo_renumber_elements,caller;
        initial_max_node_tag=nodes.max_node_tag,initial_max_element_tag=nodes.max_element_tag,save_all)
    node_mapping=Dict(old=>nodes.node_tag_map[current] for (old,current) in plan.node_tag_map)
    element_mapping=Dict(old=>cells.element_tag_map[current] for (old,current) in plan.element_tag_map)
    visibility_mapping=Dict(old=>cells.element_tag_map[current] for (old,current) in plan.visibility_tag_map)
    previous=plan.model.meshing.mesh_authority
    if previous!==nothing && !isempty(previous.surface_constructor_inputs)
        constructors=Dict{Tuple{Int,Int32},MixedMesh}()
        for (owner,source) in previous.surface_constructor_inputs
            constructor=deepcopy(source)
            data=constructor.entity_data
            data.external_node_tags[:]=UInt64[nodes.node_tag_map[tag] for tag in data.external_node_tags]
            for labels in data.external_element_tags
                labels[:]=UInt64[cells.element_tag_map[tag] for tag in labels]
            end
            constructors[owner]=constructor
        end
        cells.model.meshing.mesh_authority=Model._model_mesh_with_surface_constructor_inputs(
            cells.model.meshing.mesh_authority,constructors;caller)
    end
    # Both planners relabel identities without changing any dense row. Raw
    # Point/Discrete emptiness and Surface runtime UV remain separate from
    # the strict projected MSH metadata carried by the mesh itself.
    return merge(plan,cells,(node_tag_map=node_mapping,element_tag_map=element_mapping,
        visibility_tag_map=visibility_mapping,
        native_view_node_tag_map=nodes.node_tag_map,
        native_view_element_tag_map=cells.element_tag_map))
end

function _api_native_bipartite_splice(input,raw,caller)
    data=input.entity_data
    blocks=ElementBlock[];owners=Vector{Int32}[];labels=Vector{UInt64}[]
    old_to_new=Int32[];position=0
    for (bi,block) in enumerate(input.blocks)
        keep=raw.keep_columns[bi];kept=Set(keep)
        for column in axes(block.nodes,2)
            if column in kept
                position+=1;push!(old_to_new,Int32(position))
            else
                push!(old_to_new,0)
            end
        end
        isempty(keep) && continue
        push!(blocks,ElementBlock(block.msh,block.nodes[:,keep],block.tags[keep]))
        push!(owners,data.block_entities[bi][keep]);push!(labels,data.external_element_tags[bi][keep])
    end
    updated=Elements._native_mixed_entity_data(data.entities;
        node_entities=data.node_entities,node_parametric=data.node_parametric,
        external_node_tags=data.external_node_tags,block_entities=owners,external_element_tags=labels)
    sections=Elements.MshDataSection[]
    for section in input.data_sections
        mapped=old_to_new[section.elements]
        all(!=(0),mapped) || throw(ArgumentError(
            "$caller: experimental transfer consumes an element with unsupported view-data transfer"))
        push!(sections,Elements.MshDataSection(section.name,section.strings,section.reals,
            section.header,mapped,section.nodes,section.row_nodes,section.implicit_nodes,
            section.values,section.anchor))
    end
    partition=input.partition_data
    if partition!==nothing
        ghosts=Elements.MixedGhostElement[]
        for ghost in partition.ghost_elements
            mapped=old_to_new[ghost.element]
            mapped!=0 || throw(ArgumentError(
                "$caller: experimental transfer consumes an element with unsupported ghost transfer"))
            push!(ghosts,Elements.MixedGhostElement(mapped,ghost.partition,ghost.ghost_partitions))
        end
        partition=Elements.MixedPartitionData(partition.num_partitions;entities=partition.entities,
            ghost_entities=partition.ghost_entities,ghost_elements=ghosts)
    end
    return MixedMesh(input.coords,blocks;entity_data=updated,elementary_entities=owners,
        physical_names=input.physical_names,periodic_links=input.periodic_links,
        ancillary_sections=input.ancillary_sections,data_sections=sections,partition_data=partition)
end

function _api_native_recombine_plan(m,cached,class,caller;
        options,initial_max_node_tag,initial_max_element_tag,
        node_counter_wrapped=false,element_counter_wrapped=false,
        boundary_layer_context=false,max_nodes=typemax(Int32),max_cells=typemax(Int32),
        skip_recombination::Bool=false)
    selector=_api_native_recombine_option_integer(options,"Mesh.RecombinationAlgorithm",1,caller)
    has_surfaces=_api_native_recombine_has_surfaces(m,cached) ||
        (selector==4 && !isempty(model_entities(m,2)))
    renumber=!iszero(get(options,"Mesh.Renumber",1.0))
    run_recombination=has_surfaces && !skip_recombination
    (run_recombination || renumber) || return nothing
    if cached===nothing && !any(pair->!isempty(pair[2].node_tags) || !isempty(pair[2].element_types),
            Model._discrete_mesh_records_model(m))
        if run_recombination && _api_native_recombine_option_integer(
                options,"Mesh.RecombinationAlgorithm",1,caller)==4
            owner=first(filter(entity->entity[1]==2,model_entities(m)))[2]
            return (native_recombine_events=[(:trace,"$owner  Quadrangulation of face"),
                (:fault,_PRIMARY_BIPARTITE_EMPTY_ERROR)],
                native_recombine_finished=false,publication_required=false)
        end
        return nothing
    end
    (node_counter_wrapped || element_counter_wrapped) && throw(ArgumentError(
        "$caller: wrapped UInt64 factory contexts are not implemented by this adapter"))
    settings=run_recombination ? _api_native_recombine_settings(options,caller) : nothing
    staged=Model._model_mesh_transaction_copy(m)
    source=_api_native_recombine_source(staged,cached,class,caller;
        initial_max_node_tag,initial_max_element_tag,max_nodes,max_cells)
    input=source.mesh;data=input.entity_data
    surface_owners=Set{Int32}()
    for (bi,block) in enumerate(input.blocks)
        msh_dimension(block.msh)==2 || continue
        union!(surface_owners,data.block_entities[bi])
    end
    angles=Dict(Int(owner)=>Float64(get(m.meshing.recombine,(2,Int(owner)),45.0)) for owner in surface_owners)
    boundaries=Dict{Int,Vector{Int32}}();embedded=Dict{Int,Vector{Int32}}()
    for owner in surface_owners
        boundaries[Int(owner)]=Int32[abs(child) for (dimension,child) in
            model_boundary(m,[(2,Int(owner))],false,false,false) if dimension==1]
        embedded[Int(owner)]=Int32[abs(child) for (dimension,child) in
            get(m.embeds,(2,Int(owner)),Tuple{Int,Int}[]) if dimension==1]
    end
    topology_admitted=_api_native_recombine_topology_admitted(input)
    contexts=Dict{Int32,_PrimarySurfaceContext}()
    function surface_context(owner,rows)
        context=_api_native_recombine_surface_context(m,source,owner,rows,settings,
            topology_admitted,boundary_layer_context,caller)
        contexts[Int32(owner)]=context
        return context
    end
    native_events=Tuple{Symbol,String}[];native_finished=true
    result=if run_recombination && settings.algorithm===:bipartite
        faces=Int32[tag for (dimension,tag) in model_entities(m,2)]
        raw=_primary_bipartite_raw(input,faces,caller;abort=get(options,"General.AbortOnError",2.0))
        native_events=raw.events;native_finished=raw.finished
        (mesh=_api_native_bipartite_splice(input,raw,caller),maximum_element_tag=source.max_element_tag)
    elseif run_recombination
        _primary_recombine_mixed(input;algorithm=settings.algorithm,angles,
        minimum_quality=settings.minimum_quality,
        initial_max_element_tag=source.max_element_tag,maximum_element_tag=UInt64(typemax(Int32)),
        boundary_curves=boundaries,embedded_curves=embedded,surface_contexts=surface_context)
    else
        (mesh=input,maximum_element_tag=source.max_element_tag)
    end
    mesh=result.mesh
    nnodes(mesh)<=max_nodes && sum(size(block.nodes,2) for block in mesh.blocks;init=0)<=max_cells ||
        throw(ArgumentError("$caller: recombined result exceeds the configured resource limits"))
    output_data=mesh.entity_data
    source_positions=Dict(tag=>row for (row,tag) in enumerate(data.external_node_tags))
    runtime=Union{Nothing,Vector{Float64}}[]
    for (row,tag) in enumerate(output_data.external_node_tags)
        owner=output_data.node_entities[row]
        context=owner[1]==2 ? get(contexts,owner[2],nothing) : nothing
        moved=context!==nothing && context.reposition && context.smoothing>0 &&
            context.evaluate isa _PrimaryPlaneEvaluation
        push!(runtime,moved ? deepcopy(output_data.node_parametric[row]) :
            deepcopy(source.runtime_parameters[source_positions[tag]]))
    end
    classification=_higher_classification_from_data(mesh,caller)
    classification=_mixed_rebind_class(classification,mesh;node_parametric=runtime)
    identity=Dict(tag=>tag for tag in output_data.external_node_tags)
    retained=Set(Iterators.flatten(output_data.external_element_tags))
    elements=Dict(tag=>tag for tag in Iterators.flatten(data.external_element_tags) if tag in retained)
    covered=Tuple{Int,Int}[(dimension,Int(tag)) for (dimension,tag) in classification.public_tags.authority]
    # CAD-owned corner vertices remain on Point/Curve records after the raw
    # Face transfer consumes every cell. Those changed Face records still
    # need reconciliation even when the result has no Face-owned support.
    append!(covered,Tuple{Int,Int}[(2,Int(owner)) for owner in surface_owners])
    unique!(covered)
    _dim01_reconcile_records!(staged,covered,mesh,identity;runtime_parameters=runtime)
    # Empty native parameter registries are valid even with owned Point or
    # Discrete vertices. Reconciliation must not fabricate a populated matrix.
    for (entity,previous_record) in Model._discrete_mesh_records_model(m)
        size(previous_record.node_params,1)==0 || continue
        record=_dim01_record(staged,entity)
        record===nothing && continue
        record.node_params=zeros(0,size(previous_record.node_params,2)==0 ? 0 : length(record.node_tags))
    end
    changed_surfaces=Set{Tuple{Int,Int32}}()
    for owner in surface_owners
        isequal(_api_native_recombine_surface_signature(input,owner),
            _api_native_recombine_surface_signature(mesh,owner)) || push!(changed_surfaces,(2,owner))
    end
    initial_plan=(model=staged,mesh,class=classification,authority=classification.public_tags.authority,
        node_tag_map=identity,element_tag_map=elements,visibility_tag_map=elements,
        max_node_tag=source.max_node_tag,max_element_tag=result.maximum_element_tag,
        node_counter_wrapped=false,element_counter_wrapped=false,
        recombined=!isempty(changed_surfaces),native_runtime_parameters=runtime)
    previous=staged.meshing.mesh_authority
    carrier,carrier_rows=previous===nothing ? (nothing,Int32[]) :
        Model._model_mesh_rebound_carrier(previous,mesh;node_tag_map=identity,element_tag_map=elements)
    plan=carrier===nothing ? _higher_rebind_plan_authority!(initial_plan,previous) :
        _higher_rebind_plan_authority!(initial_plan,previous;carrier,carrier_rows,
            carrier_cells=Model._model_mesh_rebound_cell_tags(previous,elements))
    if previous!==nothing
        unchanged=Dict(owner=>constructor for (owner,constructor) in previous.surface_constructor_inputs
            if !(owner in changed_surfaces) && haskey(plan.model.meshing.mesh_authority.products,owner))
        plan.model.meshing.mesh_authority=Model._model_mesh_with_surface_constructor_inputs(
            plan.model.meshing.mesh_authority,unchanged;caller)
    end
    final_plan=renumber && native_finished ? _api_native_recombine_renumber_plan(plan,caller;
        save_all=!iszero(get(options,"Mesh.SaveAll",0.0))) : plan
    return merge(final_plan,(native_recombine_events=native_events,
        native_recombine_finished=native_finished,
        publication_required=isempty(native_events) || !isempty(changed_surfaces) || (renumber && native_finished)))
end
