# Native SDK clear calls GModel::deleteMesh, which leaves entity status intact.
# Every owned mesh/provenance object is dropped; only completed CAD status lives.
function _model_mesh_cleared_authority(previous)
    previous===nothing && return nothing
    return _ModelMeshAuthority(nothing,Dict{Tuple{Int,Int32},_ModelMeshProduct}(),
        copy(previous.done),UInt64[])
end

function _model_mesh_factory_lookup_receipt(m,mesh,identities,previous;caller="MSH native lookup")
    _model_mesh_factory_authority(m,mesh,identities,previous;caller)
    return _model_mesh_factory_lookup_fields(mesh,identities,caller)
end

function _model_mesh_factory_lookup_fields(mesh,identities,caller)
    factory=_model_mesh_factory_events(mesh,identities,caller)
    node_tags=UInt64[event[3] for event in factory.node_events]
    cell_tags=UInt64[event[3] for event in factory.cell_events]
    allunique(node_tags) && allunique(cell_tags) || throw(ArgumentError(
        "$caller: repeated raw lookup keys require an explicit native occurrence receipt"))
    function dense(tags)
        maximum_tag=maximum(tags;init=UInt64(0))
        minimum_tag=minimum(tags;init=typemax(UInt64))
        count=UInt64(length(tags))
        return (minimum_tag==1 && maximum_tag==count) || maximum_tag<10*count
    end
    node_dense=dense(node_tags);cell_dense=dense(cell_tags)
    nodes=Dict{UInt64,Bool}(raw=>true for (_,_,raw,_) in factory.node_events)
    cells=Dict{UInt64,Tuple{Int32,Int,Int}}(raw=>(Int32(msh),msh_dimension(msh),Int(owner))
        for (msh,owner,raw,_) in factory.cell_events)
    return (vector=cell_dense ? cells : Dict{UInt64,Tuple{Int32,Int,Int}}(),
        vector_max=cell_dense ? maximum(cell_tags;init=UInt64(0)) : nothing,
        map=cell_dense ? Dict{UInt64,Tuple{Int32,Int,Int}}() : cells,
        node_vector=node_dense ? nodes : Dict{UInt64,Bool}(),
        node_vector_max=node_dense ? maximum(node_tags;init=UInt64(0)) : nothing,
        node_map=node_dense ? Dict{UInt64,Bool}() : nodes,
        node_targets=Dict{UInt64,UInt64}(raw=>final for (_,_,raw,final) in factory.node_events),
        element_targets=Dict{UInt64,UInt64}(raw=>final for (_,_,raw,final) in factory.cell_events),
        factory_owner=mesh,factory_identities=identities)
end

function _model_mesh_factory_lookup_binding(receipt,actual::MixedMesh,caller)
    source=receipt.factory_owner
    factory=_model_mesh_factory_events(source,receipt.factory_identities,caller)
    expected=_model_mesh_factory_lookup_fields(source,receipt.factory_identities,caller)
    all(field->isequal(getproperty(receipt,field),getproperty(expected,field)),
        (:vector,:vector_max,:map,:node_vector,:node_vector_max,:node_map,
         :node_targets,:element_targets)) || throw(ArgumentError(
            "$caller: imported raw lookup keys change constructor identities"))
    data=actual.entity_data
    data!==nothing || throw(ArgumentError("$caller: imported lookup has no current actual metadata"))
    positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    length(positions)==nnodes(actual) || throw(ArgumentError(
        "$caller: imported lookup has duplicate current node identities"))
    for (owner,uid,point,parameter) in factory.binding.nodes
        row=get(positions,uid,0)
        row>0 && data.node_entities[row]==owner && isequal(Tuple(actual.coords[:,row]),point) &&
            isequal(data.node_parametric[row]===nothing ? nothing : Tuple(data.node_parametric[row]),parameter) ||
            throw(ArgumentError("$caller: imported lookup loses current node $uid"))
    end
    cells=Dict{UInt64,Tuple}()
    for (block,owners,labels) in zip(actual.blocks,data.block_entities,data.external_element_tags),column in axes(block.nodes,2)
        uid=labels[column]
        haskey(cells,uid) && throw(ArgumentError("$caller: imported lookup has duplicate current cell $uid"))
        cells[uid]=(Int(block.msh),owners[column],block.tags[column],Tuple(data.external_node_tags[block.nodes[:,column]]))
    end
    for (msh,owner,uid,physical,nodes) in factory.binding.cells
        get(cells,uid,nothing)==(msh,owner,physical,nodes) || throw(ArgumentError(
            "$caller: imported lookup loses current ordered cell $uid"))
    end
    return nothing
end

function _model_mesh_factory_events(mesh,identities,caller)
    factory=identities.factory_events
    Elements._mixed_msh_factory_binding(mesh,factory) || throw(ArgumentError(
        "$caller: imported mesh lacks its owned exact constructor receipt"))
    data=mesh.entity_data
    state=Elements._mixed_msh_factory_state(Elements._mixed_msh_factory_context(
        factory.initial_node_max,factory.initial_element_max))
    length(factory.node_events)==nnodes(mesh) &&
        length(factory.cell_events)==sum(size(block.nodes,2) for block in mesh.blocks;init=0) ||
        throw(ArgumentError("$caller: native constructor events do not cover actual input"))
    for (row,(dim,owner,raw,final)) in enumerate(factory.node_events)
        data.node_entities[row]==(dim,owner) && data.external_node_tags[row]==final &&
            Elements._mixed_msh_factory_tag!(state,:node,(dim,owner),raw)==final ||
            throw(ArgumentError("$caller: native node constructor event loses actual identity"))
    end
    events=Dict{Int,Vector{Tuple{Int32,UInt64}}}()
    for (msh,owner,raw,final) in factory.cell_events
        Elements._mixed_msh_factory_tag!(state,:cell,(msh,owner),raw)==final ||
            throw(ArgumentError("$caller: native cell constructor event loses actual identity"))
        push!(get!(()->Tuple{Int32,UInt64}[],events,msh),(owner,final))
    end
    for (block,owners,labels) in zip(mesh.blocks,data.block_entities,data.external_element_tags)
        get(events,Int(block.msh),nothing)==collect(zip(owners,labels)) || throw(ArgumentError(
            "$caller: native cell constructor order differs from actual family vector"))
        delete!(events,Int(block.msh))
    end
    isempty(events) && state.node_max==factory.final_node_max &&
        state.element_max==factory.final_element_max &&
        identities.cell_events==Tuple((msh,owner,raw) for (msh,owner,raw,_) in factory.cell_events) ||
        throw(ArgumentError("$caller: native constructor events or final counters disagree"))
    return factory
end

# An ingestion-owned factory receipt certifies constructor order independently
# of presence/completion. CAD import never makes a pending entity completed.
function _model_mesh_factory_authority(m,mesh::MixedMesh,identities,previous;
                                        caller="MSH native factory")
    previous===nothing || previous.mesh===nothing || throw(ArgumentError(
        "$caller: constructor replay requires a cleared actual mesh"))
    factory=_model_mesh_factory_events(mesh,identities,caller)
    data=mesh.entity_data
    available=Set((dim,Int32(tag)) for (dim,tag) in model_entities(m))
    owners=union(Set(data.node_entities),Set((msh_dimension(block.msh),owner)
        for (block,tags) in zip(mesh.blocks,data.block_entities) for owner in tags))
    all(owner->owner[1] in 0:2 &&
        owner in available &&
        !haskey(m.discrete,(owner[1],Int(owner[2]))),owners) || throw(ArgumentError(
        "$caller: constructor import requires existing native lower CAD entities"))
    all(block->block.msh in (15,1,2,3),mesh.blocks) || throw(ArgumentError(
        "$caller: constructor import requires linear lower Point/Line/Tri/Quad cells"))
    done=previous===nothing ? Set{Tuple{Int,Int32}}() : copy(previous.done)
    authority=_model_mesh_authority(mesh,done;caller)
    constructors=Dict{Tuple{Int,Int32},MixedMesh}()
    for owner in keys(authority.products)
        owner[1]==2 || continue
        constructors[owner]=_model_mesh_product(authority,2,Int(owner[2]);caller)
    end
    authority=_model_mesh_with_surface_constructor_inputs(authority,constructors;caller)
    return (authority=authority,factory=factory,
        node_high_water=factory.final_node_max,
        element_high_water=factory.final_element_max)
end
