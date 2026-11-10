# Model-based geometry plugins update the source model's actual mesh vertices.
# Scalar/vector/tensor rows remain owned by their original view and time step.

struct _PluginModelGeometryUpdate
    updates::Dict{UInt64,NTuple{3,Float64}}
    seen::Set{UInt64}
    reversed::Set{UInt64}
end

_plugin_model_affine_point(point,matrix,translation)=
    ntuple(d->sum(matrix[d][a]*point[a] for a in 1:3)+translation[d],3)

function _plugin_model_radial_point(point,center,amount)
    radial=ntuple(d->point[d]-center[d],3)
    radius=hypot(radial...)
    radius>0 || return point
    return ntuple(d->point[d]+amount*(radial[d]/radius),3)
end

function _plugin_model_geometry_plan(slot,updates,reversed,caller)
    model=slot.model
    cached,class=_dim01_actual_cache(slot.mesh,slot.class;
        physical_names=model.physical_names)
    source=Model._model_mesh_transaction_copy(model)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        slot.node_tag_max,slot.element_tag_max,typemax(Int32),typemax(Int32),String(caller))
    mapping=_dim01_ingest!(assembly,source,cached,class;preserve_source_order=true)
    for i in eachindex(assembly.nodes)
        node=assembly.nodes[i]
        point=get(updates,node.tag,nothing)
        point===nothing && continue
        all(isfinite,point) || throw(ArgumentError("$caller: nonfinite transformed node"))
        assembly.nodes[i]=_Dim01Node(node.tag,point,node.owner,node.parameter)
    end
    for cell in assembly.cells
        cell.tag in reversed || continue
        cell.nodes=cell.nodes[_mixed_reverse_permutation(cell.msh)]
    end
    plan=_dim01_finish!(assembly,source,model_entities(source),Set{Tuple{Int,Int32}}(),
        falses(length(assembly.nodes)),false,false;save_all=true,source_mesh=cached,
        source_node_positions=mapping)
    if cached isa MixedMesh && (!isempty(cached.ancillary_sections) ||
        !isempty(cached.data_sections) || cached.partition_data!==nothing)
        sections=_reclassify_data_sections(cached,plan.mesh,plan,caller)
        partition=_reclassify_partition_data(cached,plan.mesh,plan,caller)
        preserved=MixedMesh(plan.mesh.coords,plan.mesh.blocks;
            physical_names=plan.mesh.physical_names,entity_data=plan.mesh.entity_data,
            elementary_entities=plan.mesh.elementary_entities,periodic_links=plan.mesh.periodic_links,
            ancillary_sections=cached.ancillary_sections,data_sections=sections,
            partition_data=partition)
        plan=merge(plan,(mesh=preserved,class=_mixed_rebind_class(plan.class,preserved),))
    end
    return plan
end

function _plugin_model_publish_geometry!(work,caller)
    # Build every candidate mesh before publishing source models.
    plans=Tuple{_ModelSlot,Any}[]
    for (slot,entry) in work
        isempty(entry.updates) && isempty(entry.reversed) && continue
        push!(plans,(slot,_plugin_model_geometry_plan(slot,entry.updates,entry.reversed,caller)))
    end
    current=_current_slot_index_locked()
    previous=current==0 ? nothing : MODEL_SLOTS[current]
    try
        for (slot,plan) in plans
            _slot_load_locked!(slot)
            _commit_tagged_plan_locked!(plan;preserve_visibility=true,
                element_tag_map=plan.visibility_tag_map)
            _slot_capture_locked!()
        end
    finally
        previous===nothing || _slot_load_locked!(previous)
    end
    return nothing
end

function _plugin_model_geometry!(view,caller;matrix=nothing,translation=nothing,
        swap=false,center=nothing,raise=0.0,offset=0.0,timestep=0)
    _slot_capture_locked!()
    keys=_plugin_model_keys(view)
    if center!==nothing
        selected=findfirst(key->key[1]==timestep,keys)
        if selected===nothing
            _native_mesh_input_error_locked("Invalid TimeStep ($timestep) in view")
            return nothing
        end
        # Model data validates TimeStep but visits each stored mesh step.
        # NodeData de-duplicates shared vertices; element data visits each
        # occurrence, including distinct source models.
    end
    work=IdDict{_ModelSlot,_PluginModelGeometryUpdate}()
    contexts=_plugin_model_contexts()
    for key in keys
        slot=view.model_contexts[key].model_slot
        context=_plugin_model_mesh_context(slot,contexts,caller)
        mesh=context.mesh
        entry=get!(work,slot) do
            _PluginModelGeometryUpdate(Dict{UInt64,NTuple{3,Float64}}(),
                Set{UInt64}(),Set{UInt64}())
        end
        elements,_=_plugin_model_elements(view,key,caller)
        for element in elements
            index=_plugin_element_index(mesh,context.public,element.etag)
            index==0 && continue
            record=mesh_element_record(mesh,index)
            for node in axes(element.vcoords,2)
                uid=_cache_node_tag(context.public,Int(record.node_tags[node]))
                # Nodal model data visits each vertex once across all elements
                # and steps. Element/element-node data visits each occurrence.
                if view.model_kind=="NodeData"
                    uid in entry.seen && continue
                    push!(entry.seen,uid)
                end
                point=get(entry.updates,uid,
                    ntuple(d->mesh.coords[d,Int(record.node_tags[node])],3))
                if center===nothing
                    point=_plugin_model_affine_point(point,matrix,translation)
                else
                    value=_plugin_scalar_rep(element.vals,(node-1)*element.ncomp,
                        element.ncomp)
                    point=_plugin_model_radial_point(point,center,offset+raise*value)
                end
                all(isfinite,point) || throw(ArgumentError("$caller: nonfinite transformed node"))
                entry.updates[uid]=point
            end
            # Public controls show the model-data reverse operation is supplied
            # for nodal views, once per element, including shared time steps.
            swap && view.model_kind=="NodeData" && push!(entry.reversed,element.etag)
        end
    end
    return _plugin_model_publish_geometry!(work,caller)
end
