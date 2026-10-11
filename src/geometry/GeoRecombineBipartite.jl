# RecombineMesh executes the same detached GFace transfer and optional label
# maps as the public entry, using the current GEO numeric options.
function _geo_recombine_bipartite_plan!(m,context,caller)
    options=Dict{String,Float64}("Mesh.RecombinationAlgorithm"=>4.0)
    for (family,member,default) in (("General","AbortOnError",2.0),
            ("General","Terminal",1.0),("General","Verbosity",5.0),
            ("Mesh","SaveAll",0.0),("Mesh","Renumber",1.0))
        options[string(family,'.',member)]=something(
            _geo_option_number(context,family,0,member),default)
    end
    diagnostics=context.native_recombine_diagnostics===nothing ?
        Model._MODEL_NATIVE_RECOMBINE_DIAGNOSTICS[] : context.native_recombine_diagnostics
    diagnostics!==nothing || throw(ArgumentError(
        "$caller: native experimental diagnostic publication is unavailable"))
    skipped=!iszero(options["General.AbortOnError"]) && diagnostics.errors()>0
    wall=time_ns();cpu=diagnostics.cpu()
    working=_geo_mesh_working_model(m,context)
    plan=Model._model_mesh_generate_higher(working,context.mesh,context.mesh_parts,
        context.mesh_part_rows,context.mesh_node_owner,0,nothing,caller;
        operation=:native_bipartite_recombine,native_recombine_options=options,
        native_recombine_skip=skipped,initial_max_node_tag=context.mesh_node_max,
        initial_max_element_tag=context.mesh_element_max)
    skipped || diagnostics.emit(:info,"Recombining 2D mesh...",options)
    if plan!==nothing
        get(plan,:publication_required,true) && _geo_publish_actual_plan!(m,context,plan,caller)
        for (kind,message) in get(plan,:native_recombine_events,())
            if kind===:error
                context.msg_error_count+=1
                if get(options,"General.AbortOnError",2.0) in (2,3,4)
                    callback=context.native_recombine_partial_callback
                    callback===nothing || callback(m,context)
                end
                diagnostics.error(message,options)
            elseif kind===:fault
                callback=context.native_recombine_partial_callback
                callback===nothing || callback(m,context)
                throw(ErrorException(message))
            else
                diagnostics.emit(kind,message,options)
            end
        end
    end
    skipped || diagnostics.done(wall,cpu,options;empty=isempty(model_entities(m,2)))
    # Automatic numbering clears native node/element lookup state. Do not
    # reinstall the earlier Merge receipt over newly numbered live objects.
    if !iszero(options["Mesh.Renumber"]) &&
            (plan===nothing || get(plan,:native_recombine_finished,true))
        context.native_element_lookup_receipt=nothing
    end
    return nothing
end
