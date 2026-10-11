using ..SizeField: AbstractAnisoField
using ..PipelineSupport: pipeline_float, pipeline_limit, PIPELINE_DEFAULT_MAX_NODES

function _generate_metric_surface(tag;max_metric_edge=1.0,
        max_nodes=PIPELINE_DEFAULT_MAX_NODES,max_triangles=20_000_000)
    caller="API.mesh.generate_metric_surface"
    surface=_mesh_query_integer(tag,caller,"tag")
    target=pipeline_float(max_metric_edge,caller,"max_metric_edge")
    target>0 || throw(ArgumentError("$caller: max_metric_edge must be positive"))
    node_limit=pipeline_limit(max_nodes,caller,"max_nodes",typemax(Int32))
    triangle_limit=pipeline_limit(max_triangles,caller,"max_triangles",typemax(Int32)÷2)
    return lock(STATE_LOCK) do
        source_model=_model_locked()
        entities=model_entities(source_model)
        count(entity->entity[1]==2,entities)==1 && haskey(source_model.surfaces,surface) &&
            all(entity->entity[1]!=3,entities) || throw(ArgumentError(
                "$caller: requires the model's single Plane Surface and no Regions"))
        Int(OPTIONS["Mesh.ElementOrder"])==1 || throw(ArgumentError(
            "$caller: Mesh.ElementOrder must be one"))
        iszero(OPTIONS["Mesh.MeshOnlyEmpty"]) && iszero(OPTIONS["Mesh.MeshOnlyVisible"]) ||
            throw(ArgumentError("$caller: MeshOnlyEmpty and MeshOnlyVisible are not supported"))
        isempty(source_model.meshing.homology_requests) || throw(ArgumentError(
            "$caller: queued homology requests are not supported"))
        Model._metric_surface_admit(source_model,surface,caller)
        field=_session_size_field_locked(source_model)
        field isa AbstractAnisoField || throw(ArgumentError(
            "$caller: select a tensor-valued background field first"))
        previous=LAST_MESH[]
        class=previous===nothing ? nothing : _cached_classification_locked(previous)
        actual,actual_class=_dim01_actual_cache(previous,class;
            physical_names=source_model.physical_names)
        request=(tag=surface,max_metric_edge=target,max_nodes=node_limit,max_triangles=triangle_limit)
        plan=Base.inferencebarrier(_generate_higher_plan)(source_model,actual,actual_class,2,
            _session_mesh1d_options(source_model,caller),caller;
            generation_meshing=_generation_meshing_options(1),element_order=1,
            renumber=!iszero(OPTIONS["Mesh.Renumber"]),
            save_all=!iszero(OPTIONS["Mesh.SaveAll"]),size_field=field,metric_surface=request,
            initial_max_node_tag=NODE_TAG_MAX[],initial_max_element_tag=ELEMENT_TAG_MAX[],
            preserve_geometry=true)
        _commit_tagged_plan_locked!(plan;preserve_visibility=true,
            element_tag_map=plan.visibility_tag_map)
        _native_lookup_invalidate_locked!()
        return _copy_mesh(plan.model.meshing.mesh_authority.carrier)
    end
end
