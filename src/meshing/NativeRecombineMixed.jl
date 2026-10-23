# Detached recombination of actual classified mixed Surface cells.
# Untouched node and cell metadata keep their original public identities.
function _primary_recombine_mixed(mesh::Elements.MixedMesh;algorithm=:blossom,angles=Dict{Int,Float64}(),
        minimum_quality=.01,initial_max_element_tag=UInt64(0),
        maximum_element_tag=UInt64(typemax(Int32)),boundary_curves=nothing,
        embedded_curves=nothing,recombine_coordinates=mesh.coords,
        recombine_node_dimensions=nothing,protected_edges_by_surface=nothing,
        embedded_edges_by_surface=nothing,surface_contexts=nothing)
    data=mesh.entity_data
    data===nothing && throw(ArgumentError("actual classified entity data is required"))
    size(recombine_coordinates)==size(mesh.coords) || throw(ArgumentError("actual recombination coordinates must align with mixed node rows"))
    all(isfinite,recombine_coordinates) || throw(ArgumentError("actual recombination coordinates must be finite"))
    dimensions=recombine_node_dimensions===nothing ? first.(data.node_entities) : recombine_node_dimensions
    length(dimensions)==size(mesh.coords,2) || throw(ArgumentError("actual recombination ownership must align with mixed node rows"))
    # Distinct native Point-cell identities can reference one Point vertex.
    # Preserve that supported lower state while retaining duplicate checks for
    # every non-Point family.
    Elements.validate(mesh;reject_duplicate_point_cells=false).ok || throw(ArgumentError("invalid actual mixed mesh"))
    initial_max_element_tag isa Integer && !(initial_max_element_tag isa Bool) &&
        0<=initial_max_element_tag<=typemax(UInt64) || throw(ArgumentError("invalid initial factory tag"))
    maximum_element_tag isa Integer && !(maximum_element_tag isa Bool) &&
        0<=maximum_element_tag<=typemax(UInt64) || throw(ArgumentError("invalid factory tag limit"))
    maximum_tag=max(UInt64(initial_max_element_tag),
        maximum((maximum(tags;init=UInt64(0)) for tags in data.external_element_tags);init=UInt64(0)))
    limit=UInt64(maximum_element_tag)
    groups=Dict{Int32,Vector{Tuple{Int,Int}}}()
    lines=Dict{Int32,Set{Tuple{Int32,Int32}}}()
    for (bi,block) in enumerate(mesh.blocks)
        block isa Elements.ElementBlock || throw(ArgumentError("fixed connectivity is required"))
        spec=Elements.msh_spec(block.msh)
        if spec.dim==1
            size(block.nodes,1)>=2 || throw(ArgumentError("actual Curve elements require two endpoints"))
            for column in axes(block.nodes,2)
                owner=data.block_entities[bi][column]
                push!(get!(()->Set{Tuple{Int32,Int32}}(),lines,owner),minmax(block.nodes[1,column],block.nodes[2,column]))
            end
        elseif spec.dim==2 && spec.family in (:tri,:qua)
            # Native recombination pairs first-order Surface elements only;
            # cells already above order one stay stored untouched, so a
            # Surface holding exclusively high-order cells returns unchanged.
            spec.order==1 || continue
            for column in axes(block.nodes,2)
                owner=data.block_entities[bi][column]
                haskey(data.entities,(2,Int(owner))) || throw(ArgumentError("missing actual Surface entity"))
                push!(get!(()->Tuple{Int,Int}[],groups,owner),(bi,column))
            end
        end
    end
    # An entirely high-order Surface has no pairable cells: return the stored
    # mesh itself so callers keep the cache object identity.
    isempty(groups) && return (mesh=mesh,maximum_element_tag=maximum_tag)
    blocks=Elements.ElementBlock[];block_owners=Vector{Int32}[];element_tags=Vector{UInt64}[]
    coordinates_out=copy(mesh.coords)
    parameters_out=deepcopy(data.node_parametric)
    removed_nodes=falses(size(mesh.coords,2))
    # Untouched lower, Volume and high-order Surface cells remain in their
    # stored order.
    for (bi,block) in enumerate(mesh.blocks)
        spec=Elements.msh_spec(block.msh)
        spec.dim==2 && spec.family in (:tri,:qua) && spec.order==1 && continue
        push!(blocks,block);push!(block_owners,data.block_entities[bi]);push!(element_tags,data.external_element_tags[bi])
    end
    for owner in sort!(collect(keys(groups)))
        selected=groups[owner]
        rows=Int32[]
        for (bi,column) in selected;append!(rows,@view mesh.blocks[bi].nodes[:,column]);end
        unique!(rows);sort!(rows;by=row->data.external_node_tags[row])
        local_index=Dict(row=>Int32(index) for (index,row) in enumerate(rows))
        triangles=Tuple{Int,Int}[entry for entry in selected if mesh.blocks[entry[1]].msh==2]
        quadrangles=Tuple{Int,Int}[entry for entry in selected if mesh.blocks[entry[1]].msh==3]
        tri_cells=Matrix{Int32}(undef,3,length(triangles));tri_physical=Int32[]
        quad_cells=Matrix{Int32}(undef,4,length(quadrangles));quad_physical=Int32[]
        previous_tri=Dict{NTuple{3,Int32},Vector{UInt64}}()
        previous_quad=Dict{NTuple{4,Int32},Vector{UInt64}}()
        for (column,(bi,source_column)) in enumerate(triangles)
            source=mesh.blocks[bi];nodes=source.nodes[:,source_column]
            for row in 1:3;tri_cells[row,column]=local_index[nodes[row]];end
            push!(tri_physical,source.tags[source_column])
            push!(get!(()->UInt64[],previous_tri,Tuple(sort(nodes))),data.external_element_tags[bi][source_column])
        end
        for (column,(bi,source_column)) in enumerate(quadrangles)
            source=mesh.blocks[bi];nodes=source.nodes[:,source_column]
            for row in 1:4;quad_cells[row,column]=local_index[nodes[row]];end
            push!(quad_physical,source.tags[source_column])
            push!(get!(()->UInt64[],previous_quad,Tuple(sort(nodes))),data.external_element_tags[bi][source_column])
        end
        entity=data.entities[(2,Int(owner))]
        boundary=boundary_curves===nothing ? entity.boundaries : get(boundary_curves,Int(owner),Int32[])
        embedded=embedded_curves===nothing ? entity.embedded_curves : get(embedded_curves,Int(owner),Int32[])
        protected=Set{Tuple{Int32,Int32}}();embedded_edges=Set{Tuple{Int32,Int32}}()
        for curve in boundary
            for (a,b) in get(lines,Int32(abs(curve)),Set{Tuple{Int32,Int32}}())
                haskey(local_index,a) && haskey(local_index,b) || continue
                push!(protected,minmax(local_index[a],local_index[b]))
            end
        end
        for curve in embedded
            for (a,b) in get(lines,Int32(abs(curve)),Set{Tuple{Int32,Int32}}())
                haskey(local_index,a) && haskey(local_index,b) || continue
                edge=minmax(local_index[a],local_index[b]);push!(protected,edge);push!(embedded_edges,edge)
            end
        end
        # Public Curve records can supplement the stored MLine blocks.
        # Convert only actual endpoint labels; an empty Curve contributes no edge.
        public_positions=Dict(data.external_node_tags[row]=>local_index[row] for row in rows)
        if protected_edges_by_surface!==nothing
            for (a,b) in get(protected_edges_by_surface,Int(owner),())
                haskey(public_positions,a) && haskey(public_positions,b) || continue
                push!(protected,minmax(public_positions[a],public_positions[b]))
            end
        end
        if embedded_edges_by_surface!==nothing
            for (a,b) in get(embedded_edges_by_surface,Int(owner),())
                haskey(public_positions,a) && haskey(public_positions,b) || continue
                edge=minmax(public_positions[a],public_positions[b])
                push!(protected,edge);push!(embedded_edges,edge)
            end
        end
        coordinates=recombine_coordinates[:,rows]
        context=surface_contexts===nothing ? nothing : surface_contexts(Int(owner),rows)
        result=_recombine_triangles(MeshTypes.Mesh(coordinates;tris=tri_cells,tri_tag=tri_physical);
            primary_native=true,algorithm,recombine_angle=get(angles,Int(owner),45.0),
            node_owner_dims=dimensions[rows],primary_minimum_quality=minimum_quality,
            primary_quadrangles=quad_cells,primary_quadrangle_tags=quad_physical,
            protected_edges=protected,embedded_edges,primary_surface_context=context)
        if context!==nothing
            coordinates_out[:,rows]=result.mesh.coords
            removed_nodes[rows].|=result.removed_nodes
            for (local_row,row) in enumerate(rows)
                context.node_entities[local_row]==(2,owner) || continue
                (!context.reposition || context.smoothing<=0 || context.evaluate===nothing) && continue
                context.parameter_present[local_row] || continue
                parameters_out[row]=copy(result.node_parameters[:,local_row])
            end
        end
        # Public element tags stay Int32-bounded while the label history is in
        # that domain; an existing UInt64 history keeps native UInt64 counter
        # arithmetic and rejects only an actual UInt64 reservation overflow.
        ceiling=maximum_tag>typemax(Int32) ? typemax(UInt64) : limit
        result.primary_factory_births<=ceiling-maximum_tag ||
            throw(ArgumentError(ceiling>typemax(Int32) ?
                "recombined element tags exceed UInt64" :
                "recombined element tags exceed Int32"))
        origin=maximum_tag;maximum_tag+=result.primary_factory_births
        for (bi,block) in enumerate(result.mesh.blocks)
            nodes=rows[block.nodes];external=UInt64[]
            previous=block.msh==2 ? previous_tri : previous_quad
            for column in axes(nodes,2)
                birth=result.element_births[bi][column]
                if birth==0
                    key=Tuple(sort(nodes[:,column]));labels=get(previous,key,nothing)
                    labels!==nothing && !isempty(labels) || throw(ArgumentError("retained cell lost its actual identity"))
                    push!(external,popfirst!(labels))
                else
                    push!(external,Base.checked_add(origin,birth))
                end
            end
            push!(blocks,Elements.ElementBlock(block.msh,nodes,block.tags))
            push!(block_owners,fill(owner,size(nodes,2)));push!(element_tags,external)
        end
    end
    kept_nodes=findall(!,removed_nodes)
    node_positions=zeros(Int32,size(mesh.coords,2))
    for (new,old) in enumerate(kept_nodes);node_positions[old]=Int32(new);end
    if any(removed_nodes)
        isempty(mesh.ancillary_sections) || throw(ArgumentError(
            "native topology node removal with opaque ancillary sections is not implemented"))
        for index in eachindex(blocks)
            block=blocks[index]
            all(node->node_positions[node]!=0,block.nodes) || throw(ArgumentError(
                "native topology removes a node still used by another mesh consumer"))
            blocks[index]=Elements.ElementBlock(block.msh,node_positions[block.nodes],block.tags)
        end
    end
    updated=Elements.MixedEntityData(data.entities;node_entities=data.node_entities[kept_nodes],
        node_parametric=parameters_out[kept_nodes],external_node_tags=data.external_node_tags[kept_nodes],
        block_entities=block_owners,external_element_tags=element_tags)
    # Element-associated payloads use flat PRIVATE indices, not public labels.
    # Rebind surviving identities after block splicing. A consumed identity
    # needs a separately proven native data-transfer policy; do not fabricate
    # values for its replacement cells.
    new_positions=Dict{UInt64,Int32}()
    flat=0
    for tags in element_tags,tag in tags
        flat+=1;flat<=typemax(Int32) || throw(ArgumentError("mixed result exceeds Int32 cell indexing"))
        new_positions[tag]=Int32(flat)
    end
    old_to_new=Int32[]
    for tags in data.external_element_tags,tag in tags
        push!(old_to_new,get(new_positions,tag,Int32(0)))
    end
    sections=Elements.MshDataSection[]
    for section in mesh.data_sections
        mapped=Int32[]
        for index in section.elements
            target=old_to_new[index]
            target!=0 || throw(ArgumentError("recombination consumes an element with unsupported view-data transfer"))
            push!(mapped,target)
        end
        mapped_nodes=node_positions[section.nodes]
        all(!=(0),mapped_nodes) || throw(ArgumentError(
            "native topology consumes a node with unsupported view-data transfer"))
        push!(sections,Elements.MshDataSection(section.name,section.strings,section.reals,section.header,
            mapped,mapped_nodes,section.row_nodes,section.implicit_nodes,section.values,section.anchor))
    end
    partition=mesh.partition_data
    if partition!==nothing
        ghosts=Elements.MixedGhostElement[]
        for ghost in partition.ghost_elements
            target=old_to_new[ghost.element]
            target!=0 || throw(ArgumentError("recombination consumes an element with unsupported ghost transfer"))
            push!(ghosts,Elements.MixedGhostElement(target,ghost.partition,ghost.ghost_partitions))
        end
        partition=Elements.MixedPartitionData(partition.num_partitions;entities=partition.entities,
            ghost_entities=partition.ghost_entities,ghost_elements=ghosts)
    end
    links=Elements.MixedPeriodicLink[]
    for link in mesh.periodic_links
        slaves=node_positions[link.slave_nodes];masters=node_positions[link.master_nodes]
        all(!=(0),slaves) && all(!=(0),masters) || throw(ArgumentError(
            "native topology consumes a node with unsupported periodic transfer"))
        push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
            slaves,masters;affine=link.affine))
    end
    output=Elements.MixedMesh(coordinates_out[:,kept_nodes],blocks;entity_data=updated,elementary_entities=block_owners,
        physical_names=mesh.physical_names,periodic_links=links,
        ancillary_sections=mesh.ancillary_sections,data_sections=sections,
        partition_data=partition)
    Elements.validate(output;reject_duplicate_point_cells=false).ok || throw(ArgumentError("invalid detached mixed output"))
    return (mesh=output,maximum_element_tag=maximum_tag)
end
