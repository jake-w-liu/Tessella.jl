# Chain publication is another detached operation on the completed actual
# products. Its owner and row maps are rebuilt from the resulting authority
# before the GEO model or context receives any state.
function _geo_homology_error_callback(context)
    key=("General",0,"AbortOnError")
    haskey(context.option_numbers,key) ||
        context.native_mesh_error_callback!==nothing || return nothing
    return function (message)
        abort=get(context.option_numbers,key,2.0)
        if context.native_mesh_error_callback!==nothing
            # A staged API importer owns diagnostic publication. Preserve
            # this execution's count without printing before it can commit.
            context.suppress_diagnostics || (context.msg_error_count+=1)
            context.native_mesh_error_callback(String(message),abort)
        else
            _geo_msg_error!(context,message)
            (abort==2 || abort==3) && throw(ErrorException(message))
            abort==4 && exit(1)
        end
        return nothing
    end
end

function _geo_actual_homology_plan(plan,caller;error_callback=nothing)
    return Model._model_mesh_generate_higher(plan.model,plan.mesh,
        Tuple{Int,Int,Union{Mesh,MixedMesh}}[],Vector{Int32}[],nothing,
        0,nothing,caller;operation=:homology,error_callback,
        initial_max_node_tag=plan.max_node_tag,
        initial_max_element_tag=plan.max_element_tag)
end
