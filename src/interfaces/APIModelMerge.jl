# File merge consumes the literal cache/record union on an operation-owned
# model. Publication retains the current slot object held by existing views.
import ..IO

function _model_merge_complete(model,source,class,caller)
    source,class=_dim01_actual_cache(source,class;physical_names=model.physical_names)
    if source isa MixedMesh && !isempty(model.meshing.periodic_uid_history)
        source,class=_periodic_uid_history_source(model,source,class,caller)
    end
    complete=_native_view_writer_mesh((;model,mesh=source,class),caller)
    original=_model_export_source_tags(source,class,caller)
    node_map=Int32[];cell_map=Int32[]
    if original!==nothing
        data=complete.entity_data
        nodes=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
        cells=Dict(uid=>Int32(row) for (row,uid) in
            enumerate(reduce(vcat,data.external_element_tags;init=UInt64[])))
        node_map=Int32[get(nodes,uid,0) for uid in original.nodes]
        cell_map=Int32[get(cells,uid,0) for uid in original.cells]
    end
    # Merge seeds geometry declarations even before their first mesh exists.
    # Only retained actual mesh metadata needs row transport; the export-only
    # requirement for a realized link on every periodic declaration stays in
    # API.write and is not imposed on this model execution input.
    links=_model_export_periodic(source,node_map,caller)
    sections=_model_export_data_sections(source,node_map,cell_map,caller)
    ancillary=source isa MixedMesh ? source.ancillary_sections : Elements.MshAncillarySection[]
    return MixedMesh(complete.coords,complete.blocks;entity_data=complete.entity_data,
        elementary_entities=complete.elementary_entities,physical_names=complete.physical_names,
        partition_data=complete.partition_data,periodic_links=links,
        data_sections=sections,ancillary_sections=ancillary)
end

function _model_merge_section_signature(section,mesh)
    cells=reduce(vcat,mesh.entity_data.external_element_tags;init=UInt64[])
    return (section.name,section.strings,section.reals,section.header,
        cells[section.elements],mesh.entity_data.external_node_tags[section.nodes],
        section.row_nodes,section.implicit_nodes,section.values,section.anchor)
end

function _model_merge_incoming_sections(previous,current,caller)
    current isa MixedMesh || return Elements.MshDataSection[]
    isempty(current.data_sections) && return Elements.MshDataSection[]
    count=length(previous.data_sections)
    length(current.data_sections)>=count || throw(ArgumentError(
        "$caller: changed retained data requires an explicit section occurrence trace"))
    for index in 1:count
        isequal(_model_merge_section_signature(previous.data_sections[index],previous),
            _model_merge_section_signature(current.data_sections[index],current)) || throw(ArgumentError(
            "$caller: changed retained data requires an explicit section occurrence trace"))
    end
    return current.data_sections[count+1:end]
end

function _model_merge_same_product(before,after,owner,caller)
    left=Model._model_mesh_product(before,owner[1],Int(owner[2]);caller)
    right=Model._model_mesh_product(after,owner[1],Int(owner[2]);caller)
    left===nothing || right===nothing || begin
        ld=left.entity_data;rd=right.entity_data
        Set(ld.external_node_tags)==Set(rd.external_node_tags) || return false
        positions=Dict(uid=>row for (row,uid) in enumerate(rd.external_node_tags))
        for (row,uid) in enumerate(ld.external_node_tags)
            other=positions[uid]
            isequal(left.coords[:,row],right.coords[:,other]) &&
                ld.node_entities[row]==rd.node_entities[other] &&
                isequal(ld.node_parametric[row],rd.node_parametric[other]) || return false
        end
        length(left.blocks)==length(right.blocks) || return false
        for (lb,rb,lt,rt) in zip(left.blocks,right.blocks,
                ld.external_element_tags,rd.external_element_tags)
            lb.msh==rb.msh && lt==rt &&
                ld.external_node_tags[lb.nodes]==rd.external_node_tags[rb.nodes] || return false
        end
        return true
    end
    return false
end

function _model_merge_authority(model,complete,caller)
    previous=model.meshing.mesh_authority
    actual=complete.entity_data
    done=previous===nothing ? Set{Tuple{Int,Int32}}() : copy(previous.done)
    union!(done,actual.node_entities)
    for (block,owners) in zip(complete.blocks,actual.block_entities)
        union!(done,((msh_dimension(block.msh),owner) for owner in owners))
    end
    if previous===nothing || previous.mesh===nothing
        return Model._model_mesh_authority(complete,done;caller)
    end
    rebound=Model._model_mesh_rebind_authority(previous,complete)
    products=copy(rebound.products);produced=copy(rebound.producer_rows)
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(actual.external_node_tags))
    old_tags=previous.mesh.entity_data.external_node_tags
    for (owner,product) in previous.products
        haskey(products,owner) && continue
        isempty(product.blocks) || throw(ArgumentError(
            "$caller: retained actual product $owner lost its literal cells"))
        rows=Int32[get(positions,old_tags[row],0) for row in product.rows]
        all(>(0),rows) || throw(ArgumentError(
            "$caller: retained empty product $owner lost its actual support UIDs"))
        products[owner]=Model._ModelMeshProduct(owner,rows,ElementBlock[],Vector{UInt64}[])
        haskey(previous.producer_rows,owner) &&
            (produced[owner]=Int32[positions[old_tags[row]] for row in previous.producer_rows[owner]])
    end
    rebound=Model._ModelMeshAuthority(complete,products,done,rebound.private_rank,
        rebound.native_inputs,rebound.carrier,rebound.carrier_rows,produced,rebound.carrier_cells)
    inputs=Dict(owner=>input for (owner,input) in previous.native_inputs
        if _model_merge_same_product(previous,rebound,owner,caller))
    rebound=Model._ModelMeshAuthority(complete,products,done,rebound.private_rank,
        inputs,rebound.carrier,rebound.carrier_rows,produced,rebound.carrier_cells)
    constructors=Dict(owner=>input for (owner,input) in previous.surface_constructor_inputs
        if _model_merge_same_product(previous,rebound,owner,caller))
    rebound=Model._model_mesh_with_surface_constructor_inputs(rebound,constructors;caller)
    if previous.emitted_geometry!==nothing && rebound.carrier!==nothing
        parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[];rows=Vector{Int32}[]
        for (owner,(part,binding)) in sort!(collect(previous.emitted_geometry.parts);by=first)
            push!(parts,(owner[1],Int(owner[2]),part));push!(rows,binding)
        end
        rebound=Model._model_mesh_with_emitted_receipt(rebound,parts,rows;caller)
    end
    return rebound
end

function _model_merge_lookup_seed(slot)
    catalogue=get(NATIVE_LOOKUP_CATALOGUES,slot,nothing)
    catalogue===nothing && return nothing
    bindings=Dict(uid=>(binding.family,binding.dimension,binding.owner)
        for (uid,binding) in catalogue.elements if binding!==nothing)
    dense=catalogue.element_dense_max
    node_dense=catalogue.node_dense_max
    return (vector=dense===nothing ? empty(bindings) : copy(bindings),
        vector_max=dense,map=dense===nothing ? copy(bindings) : empty(bindings),
        node_vector=node_dense===nothing ? empty(catalogue.nodes) : copy(catalogue.nodes),
        node_vector_max=node_dense,node_map=node_dense===nothing ? copy(catalogue.nodes) : empty(catalogue.nodes),
        node_targets=copy(catalogue.node_targets),element_targets=copy(catalogue.element_targets))
end

function _model_merge_seed(model,source,caller)
    authority=model.meshing.mesh_authority
    carrier=authority.carrier===nothing ? authority.mesh : authority.carrier
    parts,rows=GeoExec._geo_actual_product_parts(authority,caller)
    geometry=GeoExec._geo_geometry_rebind_receipt(authority,carrier,nothing,nothing,caller)
    if authority.emitted_geometry!==nothing && carrier!==authority.mesh
        emitted=authority.emitted_geometry.parts
        entries=sort!(collect(emitted);by=first)
        geometry=GeoExec._geo_geometry_receipt(carrier,
            Tuple{Int,Int,Union{Mesh,MixedMesh}}[(owner[1],Int(owner[2]),part)
                for (owner,(part,_)) in entries],
            Vector{Int32}[binding for (_,(part,binding)) in entries],
            authority.mesh.entity_data.external_node_tags[authority.carrier_rows],caller)
    end
    binding=carrier===authority.mesh ? collect(Int32,1:nnodes(carrier)) : authority.carrier_rows
    parameters=Model._model_mesh_runtime_parameters(model,authority.mesh)
    return (mesh=carrier,mesh_parts=parts,mesh_part_rows=rows,mesh_geometry_parts=geometry,
        mesh_node_owner=GeoExec._geo_actual_carrier_owners(authority,carrier),
        mesh_node_parametric=parameters[binding],mesh_node_max=NODE_TAG_MAX[],
        mesh_element_max=ELEMENT_TAG_MAX[],native_element_lookup_receipt=_model_merge_lookup_seed(source))
end

function _model_merge_fields()
    fields=Dict(tag=>GeoExec.GeoFieldSpec(tag,spec.kind,copy(spec.options),
        copy(spec.order))
        for (tag,spec) in FIELD_SPECS[])
    return (;fields,background_field=BACKGROUND_FIELD[],
        boundary_layer_fields=copy(BOUNDARY_LAYER_FIELDS[]))
end

function _model_merge_prepare(path)
    caller="API.merge";current=_model_locked();index=_current_slot_index_locked()
    index!=0 || throw(ArgumentError("$caller: current model has no session slot"))
    source=MODEL_SLOTS[index]
    model=Model._model_mesh_transaction_copy(current)
    model.name=MODEL_NAME[]
    complete=_model_merge_complete(current,LAST_MESH[],LAST_MESH_CLASS[],caller)
    model.meshing.mesh_authority=_model_merge_authority(model,complete,caller)
    seed=_model_merge_seed(model,source,caller)
    options=copy(OPTIONS);diagnostics=_model_import_diagnostics(options)
    fields=_model_merge_fields()
    extension=lowercase(splitext(path)[2])
    result=if extension==".msh"
        GeoExec._geo_stage_mesh_merge(path,model,seed;
            context_lc=NATIVE_CONTEXT_LC[],forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[])
    elseif extension in (".unv",".mesh",".su2",".bdf",".nas",".inp",".neu",".p3d",".mail",".ir3",".diff",".vtk")
        IO._stage_foreign_mesh_merge(path,model,seed;
            context_lc=NATIVE_CONTEXT_LC[],forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[])
    elseif extension==".stl"
        IO._stage_stl_mesh_merge(path,model,seed;
            context_lc=NATIVE_CONTEXT_LC[],forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[])
    else
        _model_import_execute_geo(path,0,options,diagnostics;legacy=true,
            initial_model=model,initial_mesh_context=seed,initial_fields=fields)
    end
    final_authority=result.model.meshing.mesh_authority
    final_mesh=final_authority===nothing ? nothing : final_authority.mesh
    incoming=_model_merge_incoming_sections(complete,final_mesh,caller)
    plan=_model_import_prepare_result(result,path,options,diagnostics;
        caller,detached_result=true,view_sections=incoming)
    slot=plan.slot;slot.name=result.model.name;slot.file_name=MODEL_FILE_NAME[]
    if result.params===nothing
        slot.fields=deepcopy(FIELD_SPECS[]);slot.background_field=BACKGROUND_FIELD[]
        slot.boundary_layer_fields=copy(BOUNDARY_LAYER_FIELDS[])
        slot.field_tag_max=FIELD_TAG_MAX[]
    end
    slot.node_tag_max=result.mesh_node_max;slot.element_tag_max=result.mesh_element_max
    slot.node_counter_wrapped=final_mesh===nothing ? false : NODE_COUNTER_WRAPPED[]
    slot.element_counter_wrapped=final_mesh===nothing ? false : ELEMENT_COUNTER_WRAPPED[]
    live=slot.class===nothing || slot.class.public_tags===nothing ? Set{UInt64}() :
        Set(slot.class.public_tags.element_tags)
    slot.visibility=Dict(uid=>value for (uid,value) in ELEMENT_VISIBILITY[] if uid in live)
    slot.window_visibility=deepcopy(source.window_visibility)
    slot.window_element_visibility=Dict(window=>Dict(uid=>value for (uid,value) in values if uid in live)
        for (window,values) in source.window_element_visibility)
    return plan
end

"""
    merge(filename::AbstractString)

Merge a `.geo`, `.msh`, `.unv`, `.mesh`, `.stl`, `.su2`, `.bdf`, `.nas`,
`.inp`, `.neu`, `.p3d`, `.mail`, `.ir3`, `.diff`, `.vtk` or supported linear-list
`.pos` file atomically. Existing models,
the current model slot, associated filename and held views survive. GEO executes
against current entities, fields and options; an explicit `SetName` is honored.
Supported literal MSH identities/periodic/partition/data metadata are retained;
duplicate UID occurrences or incompatible partition namespaces reject before
publication. STL welds through `read_stl` and lands on the lowest free discrete
Surface tag.
Explicit external GEO actions are outside the rollback contract.
POS adds list views only, leaving the current model and associated filename unchanged.
"""
function merge(filename::AbstractString)
    path=_model_import_filename(filename,"API.merge")
    return lock(STATE_LOCK) do
        if lowercase(splitext(path)[2])==".pos"
            _model_locked()
            return _list_model_merge!(path,"API.merge")
        end
        plan=_model_merge_prepare(path)
        _model_import_publish!(plan;legacy=true)
        return nothing
    end
end
