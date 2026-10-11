# Tag lookup memo is separate from the authoritative mesh and spatial locator.
# Native element lookup stores the owner integer at rebuild; vertex lookup
# follows the current living vertex. Model slots survive detached plan commits.
struct _NativeElementLookupBinding
    family::Int32
    dimension::Int
    owner::Int
end

mutable struct _NativeLookupCatalogue
    nodes::Dict{UInt64,Bool}
    elements::Dict{UInt64,Union{Nothing,_NativeElementLookupBinding}}
    node_dense_max::Union{Nothing,UInt64}
    element_dense_max::Union{Nothing,UInt64}
    nodes_primed::Bool
    elements_primed::Bool
    node_targets::Dict{UInt64,UInt64}
    element_targets::Dict{UInt64,UInt64}
end
_NativeLookupCatalogue(nodes,elements,node_dense_max,element_dense_max,nodes_primed,elements_primed)=
    _NativeLookupCatalogue(nodes,elements,node_dense_max,element_dense_max,nodes_primed,elements_primed,
        Dict{UInt64,UInt64}(),Dict{UInt64,UInt64}())
_NativeLookupCatalogue()=_NativeLookupCatalogue(Dict{UInt64,Bool}(),
    Dict{UInt64,Union{Nothing,_NativeElementLookupBinding}}(),nothing,nothing,false,false)
const NATIVE_LOOKUP_CATALOGUES=IdDict{_ModelSlot,_NativeLookupCatalogue}()
# Msg's error counter is process-global, including across initialize/finalize.
# This counter records errors from the native tag-query paths ported here.
const NATIVE_LOOKUP_QUERY_ERRORS=Ref(0)
const NATIVE_MSG_LAST_ERROR=Ref("")

function _native_lookup_catalogue_locked()
    _model_locked()
    index=_current_slot_index_locked()
    index!=0 || throw(ArgumentError("API: current model has no session slot"))
    return get!(_NativeLookupCatalogue,NATIVE_LOOKUP_CATALOGUES,MODEL_SLOTS[index])
end

function _native_lookup_reset_locked!()
    empty!(NATIVE_LOOKUP_CATALOGUES)
    return nothing
end

function _native_lookup_remove_slot_locked!(slot::_ModelSlot)
    delete!(NATIVE_LOOKUP_CATALOGUES,slot)
    return nothing
end

function _native_lookup_invalidate_locked!()
    index=_current_slot_index_locked()
    index==0 && return nothing
    catalogue=get(NATIVE_LOOKUP_CATALOGUES,MODEL_SLOTS[index],nothing)
    catalogue===nothing && return nothing
    empty!(catalogue.nodes);empty!(catalogue.elements)
    empty!(catalogue.node_targets);empty!(catalogue.element_targets)
    catalogue.node_dense_max=nothing;catalogue.element_dense_max=nothing
    catalogue.nodes_primed=false;catalogue.elements_primed=false
    return nothing
end

function _native_lookup_import_receipt_locked!(receipt)
    if receipt!==nothing && hasproperty(receipt,:factory_owner)
        authority=_model_locked().meshing.mesh_authority
        authority!==nothing && authority.mesh!==nothing || throw(ArgumentError(
            "API: imported lookup lacks its current actual authority"))
        Model._model_mesh_factory_lookup_binding(receipt,authority.mesh,"API imported lookup")
    end
    _native_lookup_invalidate_locked!()
    receipt===nothing && return nothing
    catalogue=_native_lookup_catalogue_locked()
    # Native getMeshElementByTag consults its vector whenever it is present,
    # even if later sparse imports also populated a separate map.
    bindings=receipt.vector_max===nothing ? receipt.map : receipt.vector
    for (uid,(family,dimension,owner)) in bindings
        catalogue.elements[uid]=_NativeElementLookupBinding(family,dimension,owner)
    end
    catalogue.element_dense_max=receipt.vector_max
    catalogue.elements_primed=receipt.vector_max!==nothing || !isempty(bindings)
    nodes=receipt.node_vector_max===nothing ? receipt.node_map : receipt.node_vector
    merge!(catalogue.nodes,nodes)
    catalogue.node_dense_max=receipt.node_vector_max
    catalogue.nodes_primed=receipt.node_vector_max!==nothing || !isempty(nodes)
    hasproperty(receipt,:node_targets) && merge!(catalogue.node_targets,receipt.node_targets)
    hasproperty(receipt,:element_targets) && merge!(catalogue.element_targets,receipt.element_targets)
    return nothing
end

function _native_lookup_node_tags_locked()
    m=_model_locked();tags=Set{UInt64}()
    cached=LAST_MESH[]
    if cached!==nothing
        public=_mesh_public_tags(cached)
        overlay=_high_order_overlay(cached)
        count=overlay===nothing ? nnodes(cached) : nnodes(overlay)
        for row in 1:count
            push!(tags,_cache_node_tag(public,row))
        end
    end
    for (_,_,record) in _discrete_mesh_records(m),tag in record.node_tags
        push!(tags,UInt64(tag))
    end
    return tags
end

function _native_lookup_element_bindings_locked()
    m=_model_locked()
    bindings=Dict{UInt64,Union{Nothing,_NativeElementLookupBinding}}()
    occurrences=0
    ghost_aliases=_imported_ghost_uid_owners(m)
    cached=LAST_MESH[]
    if cached!==nothing
        class=_cached_classification_locked(cached)
        public=_mesh_public_tags(cached)
        overlay=_high_order_overlay(cached)
        for (msh,dim,offset,cells,owners) in _cache_catalog(cached,class)
            family=overlay!==nothing && msh==_p2_skeleton_etype(overlay) ?
                _p2_etype(overlay) : msh
            for column in axes(cells,2)
                tag=_cache_element_tag(public,offset+column)
                owner=isempty(owners) ? -1 : Int(owners[column])
                receipt=get(m.meshing.imported_entities,(dim,owner),nothing)
                if receipt!==nothing && startswith(receipt.kind,"Ghost ")
                    if !any(row->row[1]==tag,get(m.meshing.imported_ghost_elements,
                            (dim,owner),Tuple{UInt64,Int32}[]))
                        isempty(get(ghost_aliases,tag,Tuple{Int,Int32}[])) && continue
                        # The raw object supplies its payload, while only its
                        # assigned Ghost collection contributes a cache entry.
                        owner=-1
                    end
                end
                bindings[tag]=_NativeElementLookupBinding(Int32(family),dim,owner)
                owner<0 || (occurrences+=1)
            end
        end
    end
    for (dimension,owner,record) in _discrete_mesh_records(m)
        for row in eachindex(record.element_tags)
            tag=UInt64(record.element_tags[row])
            haskey(bindings,tag) && throw(ArgumentError(
                "API: element lookup has ambiguous public tag $tag"))
            bindings[tag]=_NativeElementLookupBinding(record.element_types[row],dimension,owner)
            occurrences+=1
        end
    end
    # GModel::rebuildMeshElementCache visits every entity collection, including
    # Ghost aliases of an existing element. The last entity tag supplies the
    # cached owner; it does not change the stored element or its living nodes.
    for (tag,aliases) in ghost_aliases
        binding=get(bindings,tag,nothing)
        binding===nothing && continue
        for (dimension,owner) in aliases
            dimension==binding.dimension || throw(ArgumentError(
                "API: Ghost element $tag changes its stored dimension"))
            owner==binding.owner && continue
            occurrences+=1
            if owner>binding.owner
                binding=_NativeElementLookupBinding(binding.family,dimension,Int(owner))
                bindings[tag]=binding
            end
        end
    end
    return bindings,occurrences
end

function _native_lookup_rebuild_nodes_locked!(only_if_necessary::Bool)
    catalogue=_native_lookup_catalogue_locked()
    only_if_necessary && catalogue.nodes_primed && return catalogue
    tags=_native_lookup_node_tags_locked()
    maximum_tag=max(NODE_TAG_MAX[],maximum(tags;init=UInt64(0)))
    dense=maximum_tag==length(tags) || maximum_tag<10*UInt64(length(tags))
    empty!(catalogue.nodes);empty!(catalogue.node_targets)
    for tag in tags;catalogue.nodes[tag]=true;end
    catalogue.node_dense_max=dense ? maximum_tag : nothing
    catalogue.nodes_primed=dense || !isempty(catalogue.nodes)
    return catalogue
end

function _native_lookup_rebuild_elements_locked!(only_if_necessary::Bool)
    catalogue=_native_lookup_catalogue_locked()
    only_if_necessary && catalogue.elements_primed && return catalogue
    bindings,occurrences=_native_lookup_element_bindings_locked()
    maximum_tag=max(ELEMENT_TAG_MAX[],maximum(keys(bindings);init=UInt64(0)))
    dense=maximum_tag==occurrences || maximum_tag<10*UInt64(occurrences)
    empty!(catalogue.element_targets)
    catalogue.elements=bindings
    catalogue.element_dense_max=dense ? maximum_tag : nothing
    catalogue.elements_primed=dense || !isempty(bindings)
    return catalogue
end

function _native_lookup_has_node_locked(tag::Integer)
    catalogue=_native_lookup_rebuild_nodes_locked!(true)
    tag>=0 || return false
    uid=UInt64(tag)
    present=get(catalogue.nodes,uid,false)
    if !present && (catalogue.node_dense_max===nothing || uid>catalogue.node_dense_max)
        # The native sparse map records a null entry for an unknown tag.
        catalogue.nodes[uid]=false;catalogue.nodes_primed=true
    end
    return present
end

function _native_lookup_element_locked(tag::Integer)
    catalogue=_native_lookup_rebuild_elements_locked!(true)
    tag>=0 || return nothing
    uid=UInt64(tag)
    binding=get(catalogue.elements,uid,nothing)
    if binding===nothing && (catalogue.element_dense_max===nothing || uid>catalogue.element_dense_max)
        catalogue.elements[uid]=nothing;catalogue.elements_primed=true
    end
    # A removed native cached pointer is undefined memory. Retain catalogue
    # misses, but never reconstruct or dereference a deleted element payload.
    binding===nothing && return nothing
    current=get(catalogue.element_targets,uid,uid)
    _mesh_has_public_element(_model_locked(),LAST_MESH[],current) || return nothing
    return binding
end

function _native_lookup_query_error_locked(caller,kind,tag)
    NATIVE_LOOKUP_QUERY_ERRORS[]+=1
    # getNode's native diagnostic narrows size_t to int; getElement retains
    # the complete size_t word. Msg mode 1 records the diagnostic and returns;
    # mode 2 throws that same diagnostic, including getNode's int narrowing.
    diagnostic_tag=kind=="node" ? reinterpret(Int32,UInt32(UInt64(tag)&UInt64(0xffffffff))) : tag
    NATIVE_MSG_LAST_ERROR[]="Unknown $kind $diagnostic_tag"
    _native_diagnostic_emit_locked(:error,NATIVE_MSG_LAST_ERROR[])
    abort=OPTIONS["General.AbortOnError"]
    (iszero(abort) || abort==1) || throw(ErrorException(NATIVE_MSG_LAST_ERROR[]))
    return nothing
end

function _native_lookup_require_element_locked(tag::Integer,caller)
    binding=_native_lookup_element_locked(tag)
    binding===nothing && _native_lookup_query_error_locked(caller,"element",tag)
    return binding
end

function _native_lookup_require_node_locked(tag::Integer,caller)
    _native_lookup_has_node_locked(tag) && return true
    _native_lookup_query_error_locked(caller,"node",tag)
    return false
end

function _native_lookup_mutation(f,args...;kwargs...)
    return lock(STATE_LOCK) do
        result=f(args...;kwargs...)
        _native_lookup_invalidate_locked!()
        return result
    end
end

function _native_lookup_set_order(order)
    return lock(STATE_LOCK) do
        _model_locked()
        value=Model._mesh_attr_positive_int(order,"set_order!","order")
        value>=1 || throw(ArgumentError("set_order!: order must be at least 1 (got $order)"))
        if value>1 && !iszero(OPTIONS["General.AbortOnError"]) &&
                NATIVE_LOOKUP_QUERY_ERRORS[]>0
            # GModel::setOrderN still runs its requested renumbering after
            # HighOrder's error guard, so that branch also clears lookup memo.
            if !iszero(OPTIONS["Mesh.Renumber"])
                _renumber_nodes();_renumber_elements()
            end
            return nothing
        end
        result=_set_order_native_lookup_live(order)
        _native_lookup_invalidate_locked!()
        return result
    end
end

function _get_node(node_tag)
    caller="API.mesh.get_node"
    return lock(STATE_LOCK) do
        _model_locked()
        tag=_mesh_query_uid(node_tag,caller,"node_tag")
        _native_lookup_require_node_locked(tag,caller) ||
            return Float64[],Float64[],-1,-1
        current=get(_native_lookup_catalogue_locked().node_targets,tag,tag)
        return _get_node_native_lookup_live(current)
    end
end

function _get_element(element_tag)
    caller="API.mesh.get_element"
    return lock(STATE_LOCK) do
        _model_locked()
        tag=_mesh_query_uid(element_tag,caller,"element_tag")
        binding=_native_lookup_require_element_locked(tag,caller)
        binding===nothing && return Int32(-1),UInt64[],-1,-1
        current=get(_native_lookup_catalogue_locked().element_targets,tag,tag)
        family,nodes,dimension,owner=_get_element_native_lookup_live(current)
        family==binding.family && dimension==binding.dimension || throw(ArgumentError(
            "$caller: cached public tag $tag no longer binds the original living element"))
        return family,nodes,dimension,binding.owner<0 ? owner : binding.owner
    end
end

_set_order(order)=_native_lookup_set_order(order)
_add_nodes(args...;kwargs...)=_native_lookup_mutation(_add_nodes_native_lookup_live,args...;kwargs...)
_add_elements(args...;kwargs...)=_native_lookup_mutation(_add_elements_native_lookup_live,args...;kwargs...)
_add_elements_by_type(args...;kwargs...)=_native_lookup_mutation(_add_elements_by_type_native_lookup_live,args...;kwargs...)
_reverse_mesh(args...;kwargs...)=_native_lookup_mutation(_reverse_mesh_native_lookup_live,args...;kwargs...)
_reverse_elements(args...;kwargs...)=_native_lookup_mutation(_reverse_elements_native_lookup_live,args...;kwargs...)
_renumber_nodes(args...;kwargs...)=_native_lookup_mutation(_renumber_nodes_native_lookup_live,args...;kwargs...)
_renumber_elements(args...;kwargs...)=_native_lookup_mutation(_renumber_elements_native_lookup_live,args...;kwargs...)
_remove_elements(args...;kwargs...)=_native_lookup_mutation(_remove_elements_native_lookup_live,args...;kwargs...)
_clear_mesh(args...;kwargs...)=_native_lookup_mutation(_clear_mesh_native_lookup_live,args...;kwargs...)
_refine(args...;kwargs...)=_native_lookup_mutation(_refine_native_lookup_live,args...;kwargs...)

function _generate(dim::Integer)
    return lock(STATE_LOCK) do
        _model_locked()
        dim isa Bool && throw(ArgumentError("API.mesh.generate: dim must not be Bool"))
        0<=dim<=3 || throw(ArgumentError("API.mesh.generate: dim must be in 0:3"))
        # Generator.cpp resets Msg errors before its meshing work begins.
        _native_mesh_generation_begin_locked()
        result=_generate_native_lookup_live(dim)
        _native_lookup_invalidate_locked!()
        return result
    end
end

# getPeriodicKeys consults the vertex catalog only for actual emitted exact
# nodal keys. Bulk Point adapters above use live readers and never prime it.
function _get_periodic_keys(args...;kwargs...)
    return lock(STATE_LOCK) do
        result=_get_periodic_keys_native_lookup_live(args...;kwargs...)
        master,types,master_types,keys,master_keys,master_coords,slave_coords=result
        isempty(keys) && return result
        missing=Bool[!_native_lookup_has_node_locked(uid) for uid in keys]
        any(missing) || return result
        mapped=copy(master_keys);coordinates=copy(master_coords)
        for row in eachindex(keys)
            missing[row] || continue
            mapped[row]=keys[row]
            isempty(coordinates) || (coordinates[3row-2:3row]=slave_coords[3row-2:3row])
        end
        return master,types,master_types,keys,mapped,coordinates,slave_coords
    end
end
