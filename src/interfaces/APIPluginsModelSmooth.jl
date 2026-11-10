# Model smoothing projects incident element values to their actual mesh UIDs.
# NodeData is already continuous on those identities and needs no projection.

function _plugin_model_smooth_mean(previous,value,count)
    isfinite(value) || throw(ArgumentError("Smooth: nonfinite model value"))
    count==1 && return value
    difference=value-previous
    return isfinite(difference) ? muladd(difference,1.0/count,previous) :
        muladd(previous,(count-1)/count,value/count)
end

function _plugin_model_smooth!(view,caller)
    view.model_kind=="NodeData" && return nothing
    view.model_kind in ("ElementData","ElementNodeData") || throw(ArgumentError(
        "$caller: unsupported model data type '$(view.model_kind)'"))
    _slot_capture_locked!()
    contexts=_plugin_model_contexts()
    staged_data=empty(view.model_data)
    staged_contexts=empty(view.model_contexts)
    for key in _plugin_model_keys(view)
        context=view.model_contexts[key]
        source=_plugin_model_mesh_context(context.model_slot,contexts,caller)
        mesh=source.mesh
        tags,values,nc=view.model_data[key]
        nc>0 || throw(ArgumentError("$caller: malformed model component count"))
        widths=get(context,:counts,Int32[])
        means=Dict{UInt64,Vector{Float64}}()
        counts=Dict{UInt64,Int}()
        offset=1
        for (row,tag) in enumerate(tags)
            width=row<=length(widths) ? Int(widths[row]) : length(values)÷length(tags)
            next=Base.checked_add(offset,width)
            next<=length(values)+1 || throw(ArgumentError(
                "$caller: malformed model data row extents"))
            view.model_kind=="ElementNodeData" && width%nc!=0 &&
                throw(ArgumentError(
                    "$caller: malformed model data row extents"))
            index=_plugin_element_index(mesh,source.public,tag)
            if index!=0 && (view.model_kind=="ElementNodeData" || width>=nc)
                record=mesh_element_record(mesh,index)
                # Short ElementNodeData rows leave trailing node positions
                # undefined; skip those positions rather than averaging a
                # fabricated zero into a real node mean.
                full_positions=view.model_kind=="ElementNodeData" ? width÷nc : 0
                for (position,node) in enumerate(record.node_tags)
                    view.model_kind=="ElementNodeData" &&
                        position>full_positions && continue
                    uid=_cache_node_tag(source.public,Int(node))
                    mean=get!(means,uid) do
                        zeros(nc)
                    end
                    count=Base.checked_add(get(counts,uid,0),1)
                    counts[uid]=count
                    base=view.model_kind=="ElementData" ? offset : offset+(position-1)*nc
                    for component in 1:nc
                        mean[component]=_plugin_model_smooth_mean(
                            mean[component],values[base+component-1],count)
                    end
                end
            end
            offset=next
        end
        offset==length(values)+1 || throw(ArgumentError(
            "$caller: malformed model data row extents"))
        output_tags=sort!(collect(keys(means)))
        output=Vector{Float64}(undef,Base.checked_mul(length(output_tags),nc))
        for (row,uid) in enumerate(output_tags),component in 1:nc
            output[(row-1)*nc+component]=means[uid][component]
        end
        newkey=(key[1],"NodeData")
        # Projection must not shrink the view's logical sparse storage.
        extent=max(get(context,:extent,UInt64(0)),
            _native_model_view_extent(output_tags,nnodes(mesh)))
        staged_data[newkey]=(output_tags,output,nc)
        staged_contexts[newkey]=merge(context,(
            counts=fill(Int32(nc),length(output_tags)),extent=extent))
    end
    view.model_data=staged_data
    view.model_contexts=staged_contexts
    view.model_kind="NodeData"
    return nothing
end
