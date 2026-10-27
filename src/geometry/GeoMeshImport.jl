# Public import identities belong to the actual mesh. The simplex carrier keeps
# dense rows for geometry operations and binds each retained corner cell to its
# file identity; it never becomes a second public mesh.
function _geo_import_entity_receipts!(m,actual;ghost_dimension=nothing)
    data=actual.entity_data
    partition=actual.partition_data
    entities=partition===nothing ? data.entities : merge(data.entities,partition.entities)
    top=ghost_dimension===nothing ? maximum(first,keys(data.entities);init=-1) : ghost_dimension
    ghosts=partition===nothing ? Dict{Int,Int}() :
        Dict(Int(tag)=>Int(part) for (tag,part) in partition.ghost_entities)
    for (key,entity) in entities
        dimension,tag=key
        ghost=top>0 && dimension==top && haskey(ghosts,tag)
        separated=partition!==nothing && haskey(partition.entities,key)
        kind=ghost ? ("Ghost point","Ghost curve","Ghost surface","Ghost volume")[dimension+1] :
            separated ? ("Partition point","Partition curve","Partition surface","Partition volume")[dimension+1] :
            ("Point","Discrete curve","Discrete surface","Discrete volume")[dimension+1]
        if !ghost && !separated && !haskey(m.discrete,key)
            kind=Model.model_entity_type(m,dimension,tag)
        end
        parent=ghost ? (-1,-1) : (Int(entity.parent[1]),Int(entity.parent[2]))
        parts=ghost ? Int[ghosts[tag]] : Int.(entity.partitions)
        m.meshing.imported_entities[key]=(kind=kind,parent=parent,partitions=parts)
    end
    if partition!==nothing && !isempty(partition.ghost_elements)
        # GModel traverses entities by dimension and tag when it selects the
        # Ghost entity for each partition. Existing Ghost cells stay on their
        # original entity when a later Merge selects a different entity.
        selected=Dict{Int32,Tuple{Int,Int}}()
        for key in sort!(collect(keys(m.meshing.imported_entities)))
            receipt=m.meshing.imported_entities[key]
            startswith(receipt.kind,"Ghost ") || continue
            for part in receipt.partitions
                selected[Int32(part)]=key
            end
        end
        element_tags=reduce(vcat,data.external_element_tags;init=UInt64[])
        element_positions=Dict(uid=>position for (position,uid) in enumerate(element_tags))
        additions=Dict{Tuple{Int,Int},Vector{Tuple{UInt64,Int32}}}()
        for ghost in partition.ghost_elements,part in ghost.ghost_partitions
            key=get(selected,part,nothing)
            key===nothing && continue
            members=get!(()->Tuple{UInt64,Int32}[],additions,key)
            uid=element_tags[Int(ghost.element)]
            any(member->member[1]==uid,members) || push!(members,(uid,ghost.partition))
        end
        for (key,added) in additions
            # Native pointer-ordered Ghost collections can vary between runs.
            # Use stable source-element birth order consistently for both
            # Ghost UID queries and bulk collections, throughout their life.
            sort!(added;by=member->element_positions[member[1]],alg=Base.Sort.MergeSort)
            members=get!(()->Tuple{UInt64,Int32}[],m.meshing.imported_ghost_elements,key)
            for member in added
                any(existing->existing[1]==member[1],members) || push!(members,member)
            end
        end
    end
    partition===nothing || (m.meshing.imported_num_partitions=partition.num_partitions)
    return nothing
end

function _geo_import_actual_mesh(mixed,identities,caller;ghost_dimension=nothing)
    for block in mixed.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH type $(block.msh) requires an actual special-cell authority"))
    end
    if mixed.entity_data!==nothing
        # MSH4 declares Ghost entities before its partition-entity records.
        # An empty Ghost can have only its (tag,partition) declaration: no
        # entity record, node block, or element block is required.
        data=mixed.entity_data
        partition=mixed.partition_data
        top=ghost_dimension===nothing ? maximum(first,keys(data.entities);init=-1) : ghost_dimension
        if partition!==nothing && top>0
            for (tag,_) in partition.ghost_entities
                key=(top,Int(tag))
                (haskey(data.entities,key) || haskey(partition.entities,key)) && continue
                data.entities[key]=Elements._native_mixed_entity(top,tag,ntuple(_->0.0,6))
            end
        end
        return mixed
    end
    n=nnodes(mixed)
    owners=fill((-1,Int32(0)),n)
    cell_owners=mixed.elementary_entities===nothing ?
        [zeros(Int32,size(block.nodes,2)) for block in mixed.blocks] : mixed.elementary_entities
    referenced=Dict{Tuple{Int,Int},Set{Int32}}()
    physical=Dict{Tuple{Int,Int},Vector{Int32}}()
    for (block,tags) in zip(mixed.blocks,cell_owners),column in axes(block.nodes,2)
        key=(msh_dimension(block.msh),Int(tags[column]))
        rows=get!(()->Set{Int32}(),referenced,key)
        values=get!(()->Int32[],physical,key)
        block.tags[column]>0 && !(block.tags[column] in values) && push!(values,block.tags[column])
        for row in view(block.nodes,:,column)
            push!(rows,row)
            # Native MSH2 classification follows the dimension/entity traversal,
            # independently of source element order and the MSH type buckets.
            (owners[row][1]<0 || isless(key,(owners[row][1],Int(owners[row][2])))) &&
                (owners[row]=(key[1],Int32(key[2])))
        end
    end
    retained=Int32[row for row in 1:n if owners[row][1]>=0]
    mapping=zeros(Int32,n);mapping[retained]=Int32.(1:length(retained))
    entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    for (key,selected) in referenced
        rows=sort!(collect(selected));dim,tag=key
        first_row=first(rows)
        box=dim==0 ? ntuple(i->mixed.coords[mod1(i,3),first_row],6) :
            ntuple(i->i<=3 ? minimum(view(mixed.coords,i,rows)) :
                maximum(view(mixed.coords,i-3,rows)),6)
        entities[key]=Elements._native_mixed_entity(dim,tag,box;physical_tags=physical[key])
    end
    blocks=ElementBlock[]
    for (block,tags) in zip(mixed.blocks,cell_owners)
        projected=Int32[isempty(physical[(msh_dimension(block.msh),Int(tag))]) ? 0 :
            first(physical[(msh_dimension(block.msh),Int(tag))]) for tag in tags]
        push!(blocks,ElementBlock(block.msh,mapping[block.nodes],projected))
    end
    links=Elements.MixedPeriodicLink[]
    for link in mixed.periodic_links
        slave=mapping[link.slave_nodes];master=mapping[link.master_nodes]
        all(>(0),slave) && all(>(0),master) || throw(ArgumentError(
            "$caller: an MSH2 periodic relation references an unclassified orphan node"))
        push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
            slave,master;affine=link.affine))
    end
    sections=Elements.MshDataSection[]
    for section in mixed.data_sections
        nodes=mapping[section.nodes]
        all(>(0),nodes) || throw(ArgumentError(
            "$caller: an MSH2 data view references an unclassified orphan node"))
        push!(sections,Elements.MshDataSection(section.name,section.strings,section.reals,
            section.header,section.elements,nodes,section.row_nodes,section.implicit_nodes,
            section.values,section.anchor))
    end
    data=Elements._native_mixed_entity_data(entities;node_entities=owners[retained],
        node_parametric=fill(nothing,length(retained)),
        external_node_tags=identities.node_tags[retained],block_entities=cell_owners,
        external_element_tags=identities.element_tags)
    return MixedMesh(mixed.coords[:,retained],blocks;physical_names=mixed.physical_names,
        entity_data=data,elementary_entities=cell_owners,periodic_links=links,
        ancillary_sections=mixed.ancillary_sections,data_sections=sections,
        partition_data=mixed.partition_data)
end

function _geo_import_entity_union(a,b,caller)
    output=copy(a)
    for (key,incoming) in b
        previous=get(output,key,nothing)
        previous===nothing && (output[key]=incoming;continue)
        previous.partitions==incoming.partitions && previous.parent==incoming.parent ||
            throw(ArgumentError("$caller: overlapping partition entities require a partition identity map"))
        box=previous.dim==0 ? previous.bbox : ntuple(i->i<=3 ?
            min(previous.bbox[i],incoming.bbox[i]) : max(previous.bbox[i],incoming.bbox[i]),6)
        output[key]=Elements._native_mixed_entity(previous.dim,previous.tag,box;
            physical_tags=unique(vcat(previous.physical_tags,incoming.physical_tags)),
            boundaries=unique(vcat(previous.boundaries,incoming.boundaries)),
            embedded_curves=unique(vcat(previous.embedded_curves,incoming.embedded_curves)),
            partitions=previous.partitions,parent=previous.parent)
    end
    return output
end

function _geo_join_import_actual(a,b,caller)
    a===nothing && return b
    ad=a.entity_data;bd=b.entity_data
    isempty(intersect(Set(ad.external_node_tags),Set(bd.external_node_tags))) ||
        throw(ArgumentError("$caller: duplicate imported node UID objects require occurrence identities"))
    old_tags=Set(tag for tags in ad.external_element_tags for tag in tags)
    any(tag->tag in old_tags,(tag for tags in bd.external_element_tags for tag in tags)) &&
        throw(ArgumentError("$caller: duplicate imported element UID objects require occurrence identities"))
    total=Base.checked_add(nnodes(a),nnodes(b))
    total<=typemax(Int32) || throw(ArgumentError("$caller: merged node rows exceed Int32 indexing"))
    node_offset=Int32(nnodes(a))
    element_offset=sum(size(block.nodes,2) for block in a.blocks;init=0)
    added_elements=sum(size(block.nodes,2) for block in b.blocks;init=0)
    Base.checked_add(element_offset,added_elements)<=typemax(Int32) ||
        throw(ArgumentError("$caller: merged element rows exceed Int32 indexing"))
    entities=_geo_import_entity_union(ad.entities,bd.entities,caller)
    partition=nothing;ap=a.partition_data;bp=b.partition_data
    if ap!==nothing || bp!==nothing
        count=ap===nothing ? bp.num_partitions : ap.num_partitions
        bp===nothing || bp.num_partitions==count || throw(ArgumentError(
            "$caller: merging different partition namespaces requires a partition identity map"))
        empty_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
        parts=_geo_import_entity_union(ap===nothing ? empty_entities : ap.entities,
            bp===nothing ? empty_entities : bp.entities,caller)
        ghosts=ap===nothing ? Elements.MixedGhostElement[] : copy(ap.ghost_elements)
        bp===nothing || append!(ghosts,[Elements.MixedGhostElement(
            Base.checked_add(element_offset,Int(record.element)),record.partition,
            record.ghost_partitions) for record in bp.ghost_elements])
        partition=Elements.MixedPartitionData(count;entities=parts,
            ghost_entities=unique(vcat(ap===nothing ? Tuple{Int32,Int32}[] : ap.ghost_entities,
                bp===nothing ? Tuple{Int32,Int32}[] : bp.ghost_entities)),ghost_elements=ghosts)
    end
    effective=partition===nothing ? entities : merge(entities,partition.entities)
    blocks=ElementBlock[]
    for (mesh,data,offset) in ((a,ad,Int32(0)),(b,bd,node_offset)),
            (block,owners) in zip(mesh.blocks,data.block_entities)
        projected=Vector{Int32}(undef,length(owners))
        for i in eachindex(owners)
            entity=get(effective,(msh_dimension(block.msh),Int(owners[i])),nothing)
            projected[i]=(entity===nothing || isempty(entity.physical_tags)) ? 0 : first(entity.physical_tags)
        end
        push!(blocks,ElementBlock(block.msh,block.nodes.+offset,projected))
    end
    links=vcat(a.periodic_links,[Elements.MixedPeriodicLink(link.dim,
        link.slave_entity,link.master_entity,link.slave_nodes.+node_offset,
        link.master_nodes.+node_offset;affine=link.affine) for link in b.periodic_links])
    sections=vcat(a.data_sections,[Elements.MshDataSection(section.name,section.strings,
        section.reals,section.header,section.elements.+Int32(element_offset),
        section.nodes.+node_offset,section.row_nodes,section.implicit_nodes,
        section.values,section.anchor) for section in b.data_sections])
    data=Elements._native_mixed_entity_data(entities;
        node_entities=vcat(ad.node_entities,bd.node_entities),
        node_parametric=vcat(ad.node_parametric,bd.node_parametric),
        external_node_tags=vcat(ad.external_node_tags,bd.external_node_tags),
        block_entities=vcat(ad.block_entities,bd.block_entities),
        external_element_tags=vcat(ad.external_element_tags,bd.external_element_tags))
    names=merge(copy(b.physical_names),a.physical_names)
    return MixedMesh(hcat(a.coords,b.coords),blocks;physical_names=names,entity_data=data,
        elementary_entities=data.block_entities,periodic_links=links,
        ancillary_sections=vcat(a.ancillary_sections,b.ancillary_sections),
        data_sections=sections,partition_data=partition)
end

function _geo_import_simplex_cell_bindings(mesh)
    output=[Union{Nothing,UInt64}[],Union{Nothing,UInt64}[],Union{Nothing,UInt64}[]]
    for (block,tags) in zip(mesh.blocks,mesh.entity_data.external_element_tags)
        family=msh_family(block.msh)
        msh_spec(block.msh).order==0 && continue
        slot=family===:lin ? 1 : family===:tri ? 2 : family===:tet ? 3 : 0
        slot==0 || append!(output[slot],tags)
    end
    return output
end

function _geo_import_retained_cell_bindings(previous)
    blocks=Model._model_mesh_carrier_blocks(previous.carrier)
    length(previous.carrier_cells)==length(blocks) && return previous.carrier_cells
    # Old completed buffers can lack provenance; preserve that explicit absence
    # instead of deriving an identity from equal coordinates or dense indices.
    return [fill!(Vector{Union{Nothing,UInt64}}(undef,size(cells,2)),nothing)
        for (_,cells) in blocks]
end

# Native aggregate node emission traverses entities in (dimension, tag)
# order and retains source order within each entity. Reader storage stays in
# source order; only the public import authority adopts this emission order.
function _geo_import_order_nodes(mesh)
    data=mesh.entity_data
    order=sortperm(1:nnodes(mesh);by=row->data.node_entities[row],alg=MergeSort)
    mapping=zeros(Int32,nnodes(mesh));mapping[order]=Int32.(1:length(order))
    order==collect(1:nnodes(mesh)) && return mesh,mapping
    blocks=ElementBlock[ElementBlock(block.msh,mapping[block.nodes],block.tags) for block in mesh.blocks]
    metadata=Elements._native_mixed_entity_data(data.entities;
        node_entities=data.node_entities[order],node_parametric=data.node_parametric[order],
        external_node_tags=data.external_node_tags[order],block_entities=data.block_entities,
        external_element_tags=data.external_element_tags)
    links=[Elements.MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
        mapping[link.slave_nodes],mapping[link.master_nodes];affine=link.affine)
        for link in mesh.periodic_links]
    sections=[Elements.MshDataSection(section.name,section.strings,section.reals,
        section.header,section.elements,mapping[section.nodes],section.row_nodes,
        section.implicit_nodes,section.values,section.anchor) for section in mesh.data_sections]
    return MixedMesh(mesh.coords[:,order],blocks;physical_names=mesh.physical_names,
        entity_data=metadata,elementary_entities=mesh.elementary_entities,
        periodic_links=links,ancillary_sections=mesh.ancillary_sections,
        data_sections=sections,partition_data=mesh.partition_data),mapping
end

# MSH4 installs its read-time element cache directly, including raw elements
# whose Ghost owner has no public bulk collection. Preserve the owner at this
# event: a later Unpartition changes the live mesh but does not rebuild cache.
function _geo_import_element_lookup!(context,incoming,identities)
    (identities.element_lookup_primed || identities.node_lookup_primed) || return nothing
    receipt=context.native_element_lookup_receipt
    if receipt===nothing
        receipt=(vector=Dict{UInt64,Tuple{Int32,Int,Int}}(),
            vector_max=nothing,map=Dict{UInt64,Tuple{Int32,Int,Int}}(),
            node_vector=Dict{UInt64,Bool}(),node_vector_max=nothing,node_map=Dict{UInt64,Bool}())
    end
    data=incoming.entity_data
    count=sum(length,identities.element_tags;init=0)
    maximum_tag=maximum((maximum(tags;init=UInt64(0)) for tags in identities.element_tags);
        init=UInt64(0))
    minimum_tag=minimum((minimum(tags;init=typemax(UInt64)) for tags in identities.element_tags);
        init=typemax(UInt64))
    dense=(minimum_tag==1 && maximum_tag==count) || maximum_tag<10*UInt64(count)
    vector=copy(receipt.vector);map=copy(receipt.map)
    if dense
        # Native vector resize can also shrink an existing imported cache.
        filter!(pair->first(pair)<=maximum_tag,vector)
    end
    bindings=dense ? vector : map
    for (block,tags,owners) in zip(incoming.blocks,identities.element_tags,data.block_entities),
            column in eachindex(tags)
        uid=tags[column]
        get!(bindings,uid) do
            (Int32(block.msh),msh_dimension(block.msh),Int(owners[column]))
        end
    end
    node_vector=copy(receipt.node_vector);node_map=copy(receipt.node_map)
    node_vector_max=receipt.node_vector_max
    if identities.node_lookup_primed
        node_count=length(identities.node_tags)
        node_max=maximum(identities.node_tags;init=UInt64(0))
        node_min=minimum(identities.node_tags;init=typemax(UInt64))
        node_dense=(node_min==1 && node_max==node_count) || node_max<10*UInt64(node_count)
        if node_dense
            filter!(pair->first(pair)<=node_max,node_vector)
            node_vector_max=node_max
        end
        node_bindings=node_dense ? node_vector : node_map
        for uid in identities.node_tags
            node_bindings[uid]=true
        end
    end
    context.native_element_lookup_receipt=(vector=vector,
        vector_max=dense ? maximum_tag : receipt.vector_max,map=map,
        node_vector=node_vector,node_vector_max=node_vector_max,node_map=node_map)
    return nothing
end

# Bounded native CAD import admission. It owns its MSH constructor event stream,
# detached actual input and final UID maps; generic/discrete imports retain their
# existing parser and strict Volume provenance guards.
function _geo_native_lower_factory_scope(m,context,file)
    context.mesh===nothing || return false
    previous=m.meshing.mesh_authority
    previous===nothing || previous.mesh===nothing || return false
    isempty(m.points) && return false
    available=Set((dim,Int32(tag)) for (dim,tag) in Model.model_entities(m)
        if !haskey(m.discrete,(dim,tag)))
    open(file,"r") do io
        eof(io) && return false
        strip(readline(io))=="\$MeshFormat" || return false
        split(readline(io))==["4.1","0","8"] || return false
        strip(readline(io))=="\$EndMeshFormat" || return false
        have_nodes=false;have_elements=false;have_entities=false
        while !eof(io)
            header=strip(readline(io));isempty(header) && continue
            if header=="\$PhysicalNames"
                count=parse(Int,strip(readline(io)))
                count>=0 || return false
                for _ in 1:count;readline(io);end
                strip(readline(io))=="\$EndPhysicalNames" || return false
            elseif header=="\$Entities"
                counts=parse.(Int,split(readline(io)))
                length(counts)==4 && all(>=(0),counts) || return false
                for dim in 0:3,_ in 1:counts[dim+1]
                    record=split(readline(io));isempty(record) && return false
                    owner=(dim,parse(Int32,record[1]));owner in available || return false
                end
                strip(readline(io))=="\$EndEntities" || return false
                have_entities=true
            elseif header=="\$Nodes"
                # Header probing never consumes or maps a raw identity. The actual
                # bounded reader owns parsing and allocation after scope admission.
                while !eof(io) && strip(readline(io))!="\$EndNodes";end
                have_nodes=true
            elseif header=="\$Elements"
                declaration=parse.(UInt64,split(readline(io)))
                length(declaration)==4 || return false
                declaration[1]<=typemax(Int32) || return false
                for _ in 1:Int(declaration[1])
                    block=parse.(Int64,split(readline(io)))
                    length(block)==4 && block[1] in 0:2 && block[3] in (15,1,2,3) &&
                        block[4]>=0 && block[4]<=typemax(Int32) || return false
                    owner=(Int(block[1]),Int32(block[2]))
                    owner in available && msh_dimension(Int(block[3]))==owner[1] || return false
                    for _ in 1:Int(block[4]);readline(io);end
                end
                strip(readline(io))=="\$EndElements" || return false
                have_elements=true
            else
                return false
            end
        end
        return have_entities && have_nodes && have_elements
    end
end

function _geo_import_cleared_native_factory!(m,context,file,caller)
    _geo_native_lower_factory_scope(m,context,file) || return false
    previous=m.meshing.mesh_authority
    switch_tags=_geo_mesh_switch_tags(context)
    mixed,identities=Elements._read_mixed_msh_factory_with_identities(file,
        Elements._mixed_msh_factory_context(context.mesh_node_max,context.mesh_element_max);
        switch_tags)
    plan=Model._model_mesh_factory_authority(m,mixed,identities,previous;caller)
    actual=plan.authority.mesh
    carrier=_geo_mixed_simplex_mesh(actual,caller)
    rows=collect(Int32,1:nnodes(actual))
    authority=Model._model_mesh_authority(actual,plan.authority.done;caller,carrier,
        carrier_rows=rows,carrier_cells=_geo_import_simplex_cell_bindings(actual),
        private_rank=plan.authority.private_rank)
    authority=Model._model_mesh_with_surface_constructor_inputs(authority,
        plan.authority.surface_constructor_inputs;caller)
    parts,part_rows=_geo_actual_product_parts(authority,caller)
    owners=_geo_actual_carrier_owners(authority,carrier)
    geometry=_geo_geometry_rebind_receipt(authority,carrier,nothing,nothing,caller)
    parameters=Model._model_mesh_runtime_parameters(m,actual)[rows]
    lookup=Model._model_mesh_factory_lookup_receipt(m,mixed,identities,previous;caller)
    # Everything above is detached and certified before model/context publication.
    _geo_import_entity_receipts!(m,actual)
    m.meshing.mesh_authority=authority
    context.mesh=carrier;context.mesh_parts=parts;context.mesh_part_rows=part_rows
    context.mesh_geometry_parts=geometry;context.mesh_node_owner=owners
    context.mesh_node_parametric=parameters
    context.mesh_node_max=plan.node_high_water
    context.mesh_element_max=plan.element_high_water
    context.native_element_lookup_receipt=lookup
    return true
end

# A real file frontend can stage this same importer without generating a GEO
# wrapper or installing a temporary model in the live API session.
function _geo_stage_mesh_file(file;context_lc=1.0,forced_bbox=false,
                              switch_tags=false)
    model=GeoModel();context=_GeoNumericContext()
    model.name=splitext(basename(file))[1]
    context.file_name=String(file)
    context.native_context_lc=context_lc
    context.native_context_forced_bbox=forced_bbox
    # `Mesh.SwitchElementTags` is a stored option like any other — the API
    # session's option table seeds it into the staged context.
    switch_tags &&
        (context.option_numbers[("Mesh",0,"SwitchElementTags")]=1.0)
    _geo_import_mesh_file!(model,context,file,"API.open")
    return (model=model,mesh=context.mesh,params=nothing,transfinite_tri=nothing,
        mesh_node_max=context.mesh_node_max,mesh_element_max=context.mesh_element_max,
        native_context_lc=context.native_context_lc,
        native_context_forced_bbox=context.native_context_forced_bbox,
        native_element_lookup_receipt=context.native_element_lookup_receipt)
end

# Seed only actual operation-owned mesh state, retaining carrier/event bindings.
function _geo_seed_mesh_context!(context,seed)
    for name in (:mesh,:mesh_parts,:mesh_part_rows,:mesh_geometry_parts,
            :mesh_node_owner,:mesh_node_parametric,:mesh_node_max,:mesh_element_max,
            :native_element_lookup_receipt)
        setfield!(context,name,getproperty(seed,name))
    end
    return nothing
end

function _geo_stage_mesh_merge(file,model,seed;context_lc=1.0,forced_bbox=false,
                               switch_tags=false)
    context=_GeoNumericContext();_geo_seed_mesh_context!(context,seed)
    context.file_name=String(file)
    context.native_context_lc=context_lc
    context.native_context_forced_bbox=forced_bbox
    switch_tags &&
        (context.option_numbers[("Mesh",0,"SwitchElementTags")]=1.0)
    _geo_import_mesh_file!(model,context,file,"API.merge")
    return (model=model,mesh=context.mesh,params=nothing,transfinite_tri=nothing,
        mesh_node_max=context.mesh_node_max,mesh_element_max=context.mesh_element_max,
        native_context_lc=context.native_context_lc,
        native_context_forced_bbox=context.native_context_forced_bbox,
        native_element_lookup_receipt=context.native_element_lookup_receipt)
end

# `Mesh.SwitchElementTags` — stored like any other `.geo` numeric option,
# consumed when an MSH 2.2 file is merged/imported through this context.
_geo_mesh_switch_tags(context)=
    !iszero(_geo_option_number(context,"Mesh",0,"SwitchElementTags"))

function _geo_import_mesh_file!(m,context,file,caller)
    _geo_import_cleared_native_factory!(m,context,file,caller) && return nothing
    switch_tags=_geo_mesh_switch_tags(context)
    mixed,identities=Elements._read_mixed_msh_with_identities(file;switch_tags)
    # Native MSH4 creates ghosts at the live model dimension after ordinary
    # file entities are read, before partition entities are constructed.
    file_dimension=mixed.entity_data===nothing ? -1 :
        maximum(first,keys(mixed.entity_data.entities);init=-1)
    ghost_dimension=max(file_dimension,maximum(first,Model.model_entities(m);init=-1))
    incoming=_geo_import_actual_mesh(mixed,identities,caller;ghost_dimension)
    previous=m.meshing.mesh_authority
    if context.mesh!==nothing && (previous===nothing || previous.mesh===nothing ||
            previous.carrier!==context.mesh)
        working=_geo_mesh_working_model(m,context)
        plan=Model._model_mesh_generate_higher(working,context.mesh,context.mesh_parts,
            context.mesh_part_rows,context.mesh_node_owner,0,nothing,caller;
            operation=:bind_completed,initial_max_node_tag=context.mesh_node_max,
            initial_max_element_tag=context.mesh_element_max)
        _geo_publish_actual_plan!(m,context,plan,caller;carrier=context.mesh)
        previous=m.meshing.mesh_authority
    end
    actual=_geo_join_import_actual(previous===nothing ? nothing : previous.mesh,incoming,caller)
    actual,node_mapping=_geo_import_order_nodes(actual)
    imported_carrier=_geo_mixed_simplex_mesh(incoming,caller)
    carrier=context.mesh===nothing ? imported_carrier : _geo_concat_meshes(context.mesh,imported_carrier)
    old_rows=previous===nothing ? Int32[] : previous.carrier_rows
    actual_offset=previous===nothing ? 0 : nnodes(previous.mesh)
    rows=node_mapping[vcat(old_rows,Int32.(actual_offset.+(1:nnodes(incoming))))]
    incoming_cells=_geo_import_simplex_cell_bindings(incoming)
    cells=if previous===nothing
        incoming_cells
    elseif previous.carrier isa Mesh
        [vcat(old,new) for (old,new) in zip(_geo_import_retained_cell_bindings(previous),incoming_cells)]
    else
        vcat(_geo_import_retained_cell_bindings(previous),
            [tags for (tags,(_,nodes)) in zip(incoming_cells,
                Model._model_mesh_carrier_blocks(imported_carrier)) if size(nodes,2)>0])
    end
    done=previous===nothing ? Set{Tuple{Int,Int32}}() : copy(previous.done)
    union!(done,incoming.entity_data.node_entities)
    for (block,owners) in zip(incoming.blocks,incoming.entity_data.block_entities)
        union!(done,((msh_dimension(block.msh),owner) for owner in owners))
    end
    ranks=previous===nothing ? UInt64[] : copy(previous.private_rank)
    next_rank=maximum(ranks;init=UInt64(0))
    UInt64(nnodes(incoming))<=typemax(UInt64)-next_rank ||
        throw(ArgumentError("$caller: imported private vertex rank overflows UInt64"))
    append!(ranks,[next_rank+UInt64(row) for row in 1:nnodes(incoming)])
    ordered_ranks=similar(ranks);ordered_ranks[node_mapping]=ranks
    authority=Model._model_mesh_authority(actual,done;caller,carrier,carrier_rows=rows,
        carrier_cells=cells,private_rank=ordered_ranks)
    if previous!==nothing
        native_inputs=copy(previous.native_inputs)
        producer_rows=Dict(owner=>node_mapping[positions] for (owner,positions) in previous.producer_rows)
        touched=union(Set(incoming.entity_data.node_entities),Set(
            (msh_dimension(block.msh),owner) for (block,owners) in
                zip(incoming.blocks,incoming.entity_data.block_entities) for owner in owners))
        for owner in touched
            delete!(native_inputs,owner);delete!(producer_rows,owner)
        end
        authority=Model._ModelMeshAuthority(authority.mesh,authority.products,authority.done,
            authority.private_rank,native_inputs,carrier,rows,producer_rows,authority.carrier_cells)
        # Untouched Surface products keep their actual producer constructor
        # vectors. Their literal UID/owner/cyclic connectivity is revalidated;
        # imports that add nodes or cells on that owner invalidate its receipt.
        constructors=Dict(owner=>constructor for (owner,constructor) in
            previous.surface_constructor_inputs if !(owner in touched))
        authority=Model._model_mesh_with_surface_constructor_inputs(authority,constructors;caller)
    end
    # Construct the complete candidate and all product slices before records
    # are changed; unsupported import guards cannot leave half-imported cells.
    parts,part_rows=_geo_actual_product_parts(authority,caller)
    carrier_owners=_geo_actual_carrier_owners(authority,carrier)
    geometry_receipt=_geo_geometry_rebind_receipt(authority,carrier,
        context.mesh_geometry_parts,nothing,caller)
    if previous!==nothing && previous.emitted_geometry!==nothing && geometry_receipt!==nothing
        authority=Model._model_mesh_with_emitted_receipt(authority,
            geometry_receipt.parts,geometry_receipt.rows;caller)
    end
    _geo_merge_discrete!(m,incoming,caller)
    # Later imports can raise the model dimension; previous Ghost entities
    # retain the kind and partition chosen at their own read-time event.
    _geo_import_entity_receipts!(m,incoming;ghost_dimension)
    m.meshing.mesh_authority=authority
    context.mesh=carrier
    context.mesh_parts=parts;context.mesh_part_rows=part_rows
    context.mesh_geometry_parts=geometry_receipt
    context.mesh_node_owner=carrier_owners
    context.mesh_node_parametric=Model._model_mesh_runtime_parameters(m,actual)[rows]
    context.mesh_node_max=max(context.mesh_node_max,maximum(identities.node_tags;init=UInt64(0)))
    context.mesh_element_max=max(context.mesh_element_max,
        maximum((maximum(tags;init=UInt64(0)) for tags in identities.element_tags);init=UInt64(0)))
    _geo_import_element_lookup!(context,incoming,identities)
    return nothing
end
