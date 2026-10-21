# Plain list records have a literal API storage layout. The full POS decoder
# remains richer; unsupported text/schemes/high-order records never disappear.
const _LIST_IMPORT_TYPE_ORDER=String[
    head*shape for shape in ("P","L","T","Q","S","H","I","Y")
    for head in ("S","V","T")]

function _list_model_views_prepare(path,caller)
    # Session options feed `Family.Option` reads inside record expressions
    # (e.g. `T2(General.GraphicsWidth/4,…)`), then the DefaultOptions table.
    parsed=Post.read_pos(path;options=OPTIONS)
    views=_SessionView[]
    for source in parsed
        isempty(source.interpolations) &&
            all(element->isempty(element.suffix) && element.scheme===nothing,source.elements) ||
            throw(ArgumentError("$caller: POS session import supports plain linear list records without interpolation schemes"))
        view=_SessionView(source.name)
        view.list_time=copy(source.time)
        # Native empty views have no data steps, even when TIME is declared.
        view.numbers["NbTimeStep"]=isempty(source.elements) ? 0.0 :
            Float64(size(first(source.elements).values,3))
        # T2/T3 text records keep upstream's PViewDataList layout: one list
        # entry per record, the record's per-step strings '\0'-joined (T2C/
        # T3C), and 4 coords per record — x, y, z (0 for T2), style (T2D/T3D
        # minus the char-buffer index the string side retains implicitly).
        for text in source.texts
            key=text.dim==2 ? "T2" : "T3"
            if !(key in view.list_string_types)
                push!(view.list_string_types,key)
                view.list_strings[key]=String[]
                view.list_string_coords[key]=Float64[]
                view.list_string_counts[key]=0
            end
            push!(view.list_strings[key],join(text.strings,'\0'))
            append!(view.list_string_coords[key],
                    (text.coords[1],text.coords[2],text.coords[3],text.style))
            view.list_string_counts[key]+=1
        end
        for element in source.elements
            components=size(element.values,1)
            type_key=Post._POS_RECORD_NAME[(element.kind,components)]
            kind,nodes,expected_components=Post._POS_RECORD_TYPES[type_key]
            kind==element.kind && size(element.coords)==(3,nodes) &&
                components==expected_components || throw(ArgumentError(
                    "$caller: POS record does not match its actual linear list layout"))
            count=Base.checked_add(get(view.list_counts,type_key,0),1)
            count<=typemax(Int32) || throw(ArgumentError("$caller: list element count is out of range"))
            view.list_counts[type_key]=count
            values=get!(()->Float64[],view.list_data,type_key)
            for axis in 1:3;append!(values,@view(element.coords[axis,:]));end
            append!(values,vec(element.values))
        end
        view.list_types=String[type_key for type_key in _LIST_IMPORT_TYPE_ORDER
            if haskey(view.list_counts,type_key)]
        push!(views,view)
    end
    tags=Int[];next_tag=VIEW_TAG_MAX[]
    for _ in views
        next_tag=Base.checked_add(next_tag,1)
        while haskey(SESSION_VIEWS,next_tag);next_tag=Base.checked_add(next_tag,1);end
        0<=next_tag<=typemax(Int32) || throw(ArgumentError("$caller: view tag is out of range"))
        push!(tags,next_tag)
    end
    return (;views,tags,next_tag)
end

function _list_model_open_prepare(path,caller)
    prepared=_list_model_views_prepare(path,caller)
    model=GeoModel();model.name=splitext(basename(path))[1]
    # This is an actual empty model import result, without executing GEO.
    result=(model=model,mesh=nothing,params=nothing,transfinite_tri=nothing,
        mesh_node_max=UInt64(0),mesh_element_max=UInt64(0),
        native_context_lc=NATIVE_CONTEXT_LC[],
        native_context_forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[],
        native_element_lookup_receipt=nothing)
    options=copy(OPTIONS);diagnostics=_model_import_diagnostics(options)
    plan=_model_import_prepare_result(result,path,options,diagnostics;caller,detached_result=true)
    return Base.merge(plan,(views=prepared.views,view_tags=prepared.tags,next_tag=prepared.next_tag))
end

function _list_model_merge!(path,caller)
    prepared=_list_model_views_prepare(path,caller)
    for (tag,view) in zip(prepared.tags,prepared.views);SESSION_VIEWS[tag]=view;end
    VIEW_TAG_MAX[]=prepared.next_tag
    return nothing
end
