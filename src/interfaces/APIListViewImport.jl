# List records keep upstream's flat API storage layout — per type name, all x
# coordinates, all y, all z, then component/node values per time step. Suffixed
# record names keep their suffix so distinct bound layouts stay separate; the
# bound family scheme lives in the view's interpolation-matrices storage.
const _LIST_IMPORT_TYPE_ORDER=String[
    head*shape for shape in ("P","L","T","Q","S","H","I","Y")
    for head in ("S","V","T")]

# A suffix-bound family stores the resolved scheme's full four-matrix form so
# its geometry width survives the flat storage; exponent tables keep the
# family's parametric dimension like upstream's built-in schemes.
function _list_import_suffix_matrices(kind,scheme::Post.PosScheme)
    dim=Post._POS_FAMILY_DIMS[kind]
    matrices=Matrix{Float64}[copy(scheme.coefval),
        copy(scheme.expval[:,1:dim])]
    isempty(scheme.coefgeo) || append!(matrices,Matrix{Float64}[
        copy(scheme.coefgeo),copy(scheme.expgeo[:,1:dim])])
    return matrices
end

# One parsed `.pos` view (`Post.PosView`) as the session `_SessionView`
# layout — shared by the `.pos` import path and the `.geo` `Merge`/`Plugin`
# view bridge.
function _session_view_from_pos(source,caller)
    view=_SessionView(source.name)
    view.list_time=copy(source.time)
    # Native empty views have no data steps, even when TIME is declared.
    view.numbers["NbTimeStep"]=isempty(source.elements) ? 0.0 :
        Float64(size(first(source.elements).values,3))
    # INTERPOLATION_SCHEME records and record-name suffixes bind a family
    # scheme in file order — first binding wins. `scheme_source` keeps the
    # interpolation record index per family (0 for suffix-bound schemes)
    # so session storage mirrors what upstream's binding stores.
    bound,scheme_source=Post._pos_resolve_schemes(
        Symbol[element.kind for element in source.elements],
        String[element.suffix for element in source.elements],
        Union{Nothing,Post.PosScheme}[nothing for _ in source.elements],
        source.interpolations,source.scheme_positions,caller)
    for (kind,scheme) in bound
        index=scheme_source[kind]
        view.interp[Post._POS_FAMILY_CODES[kind]]=index==0 ?
            _list_import_suffix_matrices(kind,scheme) :
            Matrix{Float64}[copy(matrix) for matrix in
                source.interpolations[index]]
    end
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
        base_key=Post._POS_RECORD_NAME[(element.kind,components)]
        kind,arity,expected_components=Post._POS_RECORD_TYPES[base_key]
        scheme=element.scheme
        ngeo=scheme===nothing ? arity : Post._pos_scheme_geom_nodes(scheme,kind)
        nval=scheme===nothing ? arity : Post._pos_scheme_value_nodes(scheme)
        kind==element.kind && components==expected_components &&
            size(element.coords)==(3,ngeo) &&
            size(element.values,2)==nval || throw(ArgumentError(
                "$caller: POS record does not match its bound interpolation scheme"))
        # Verbatim record names keep suffixed (high-order or otherwise
        # rebound) payload layouts distinguishable inside the flat lists.
        type_key=base_key*element.suffix
        count=Base.checked_add(get(view.list_counts,type_key,0),1)
        count<=typemax(Int32) || throw(ArgumentError("$caller: list element count is out of range"))
        view.list_counts[type_key]=count
        values=get!(()->Float64[],view.list_data,type_key)
        for axis in 1:3;append!(values,@view(element.coords[axis,:]));end
        append!(values,vec(element.values))
    end
    # Canonical type ordering first (upstream's fixed list table), then
    # suffixed variants ordered after their base name.
    base_order=Dict(name=>index for (index,name) in
        enumerate(_LIST_IMPORT_TYPE_ORDER))
    view.list_types=sort!(collect(keys(view.list_counts));
        by=key->(get(base_order,ncodeunits(key)>=2 ? key[1:2] : key,
                     typemax(Int)),
                 ncodeunits(key)>2 ? String(key[3:end]) : ""))
    return view
end

function _list_model_views_prepare(path,caller)
    # Session options feed `Family.Option` reads inside record expressions
    # (e.g. `T2(General.GraphicsWidth/4,…)`), then the DefaultOptions table.
    parsed=Post.read_pos(path;options=OPTIONS)
    views=_SessionView[_session_view_from_pos(source,caller)
        for source in parsed]
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
