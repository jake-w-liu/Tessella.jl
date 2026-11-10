# GEO declarations update complete parameter records on their operation-owned
# store. Values remain shared until a writer replaces them; metadata has its
# own containers so failed file preparation cannot mutate session records.
function _onelab_metadata_copy(parameter)
    return _OnelabParameter(parameter.name,parameter.kind,parameter.values_f,parameter.values_s,
        parameter.min,parameter.max,parameter.step,copy(parameter.choices_f),copy(parameter.choices_s),
        copy(parameter.value_labels),copy(parameter.attributes),parameter.label,parameter.short_help,
        parameter.long_help,parameter.read_only,parameter.visible,parameter.closed,parameter.kind_name)
end

function _onelab_geo_declare!(runner,kind,defaults,fopts,copts,labels,caller)
    _onelab_geo_sync!(runner)
    store=runner.store;name=first(copts["Name"]);key=(kind,name)
    runner.remote===nothing || _onelab_remote_query!(runner,kind,name)
    return _onelab_store_locked(store,runner.epoch) do
        all(values->all(isfinite,values),Base.values(fopts)) &&
            (kind!==:number || all(isfinite,defaults)) || throw(ArgumentError(
                "$caller: ONELAB numeric values and options must be finite"))
        old=get(store.parameters,key,nothing)
        parameter=old===nothing ? _OnelabParameter(name,kind) : _onelab_metadata_copy(old)
        if kind===:number
            parameter.min===nothing && (parameter.min=-1e200)
            parameter.max===nothing && (parameter.max=1e200)
            parameter.step===nothing && (parameter.step=0.)
        elseif old===nothing
            parameter.kind_name="generic"
        end
        for (option,values) in fopts
            isempty(values) && continue
            value=first(values)
            if option=="ReadOnly"
                parameter.read_only=!iszero(value)
            elseif option=="Visible"
                parameter.visible=!iszero(value)
            elseif option=="Closed"
                parameter.closed=!iszero(value)
                parameter.attributes[option]=_onelab_wire_f64(value)
            elseif kind===:number && option in ("Min","Max","Step")
                setfield!(parameter,option=="Min" ? :min : option=="Max" ? :max : :step,value)
            elseif kind===:number && option=="Choices"
                parameter.choices_f=copy(values)
            elseif option!="Range"
                parameter.attributes[option]=_onelab_wire_f64(value)
            end
        end
        for (option,values) in copts
            option=="Name" && continue
            isempty(values) && continue
            value=first(values)
            if option=="Label"
                parameter.label=value
            elseif option in ("Help","ShortHelp")
                parameter.short_help=value
            elseif option=="LongHelp"
                parameter.long_help=value
            else
                parameter.attributes[option]=value
                if kind===:string && option=="Choices"
                    parameter.choices_s=copy(values)
                elseif kind===:string && option=="Kind"
                    parameter.kind_name=value
                end
            end
        end
        isempty(labels) || (parameter.value_labels=Dict(_onelab_wire_f64(value)=>label for (value,label) in labels))
        uses_default=old===nothing || parameter.read_only
        if uses_default
            if kind===:number
                parameter.values_f=copy(defaults)
            else
                parameter.values_s=copy(defaults)
            end
        end
        changed=old===nothing || old.values_f!=parameter.values_f || old.values_s!=parameter.values_s
        store.parameters[key]=parameter
        changed && _onelab_store_mark_changed!(store,key)
        previous=get(store.wire,key,nothing)
        new_client=previous===nothing || !haskey(previous.clients,"Gmsh")
        wire=_onelab_store_register_client!(store,parameter,"Gmsh")
        if changed || new_client
            wire.clients["Gmsh"]=wire.changed_value
            store.client_changed["Gmsh"]=get(store.client_changed,"Gmsh",false) || wire.changed_value!=0
        end
        runner.remote===nothing || _onelab_remote_publish!(store,runner.remote,key)
        if kind===:number
            projection=isempty(parameter.values_f) ? 0. : first(parameter.values_f)
            runner.numbers[name]=projection;runner.seed_numbers[name]=projection
            return uses_default || isempty(parameter.values_f) ? copy(defaults) : Float64[projection]
        else
            projection=isempty(parameter.values_s) ? "" : first(parameter.values_s)
            runner.strings[name]=projection;runner.seed_strings[name]=projection
            return projection
        end
    end
end
