# File execution/import prepares detached state before touching model slots,
# options, lookup catalogues, or session views.
import ..IO

function _model_import_option!(options,family,member,value)
    name=string(family,'.',member)
    key=name=="Mesh.MinimumCirclePoints" ? "Mesh.MinCircleNodes" :
        name=="Mesh.MinimumCurvePoints" ? "Mesh.MinCurveNodes" :
        get(_OPTION_ALIASES,name,name)
    haskey(options,key) || return nothing
    options[key]=value
    if key in ("Mesh.MeshSizeFactor","Mesh.CharacteristicLengthFactor")
        options[key=="Mesh.MeshSizeFactor" ? "Mesh.CharacteristicLengthFactor" :
            "Mesh.MeshSizeFactor"]=value
    end
    return nothing
end

function _model_import_diagnostics(options)
    errors=Ref(NATIVE_LOOKUP_QUERY_ERRORS[]);last_error=Ref(NATIVE_MSG_LAST_ERROR[])
    events=Tuple{Symbol,String,Dict{String,Float64}}[]
    emit=(kind,message,values)->(push!(events,(kind,String(message),copy(values)));nothing)
    reject=(message,values)->begin
        errors[]+=1;last_error[]=String(message);emit(:error,message,values)
        abort=get(values,"General.AbortOnError",2.0)
        (abort==2 || abort==3) && throw(ErrorException(message))
        abort==4 && exit(1)
        nothing
    end
    begin_mesh=()->(errors[]=0;last_error[]="";nothing)
    done=(wall,cpu,values;empty=false)->begin
        elapsed=empty ? 0.0 : max(0.0,Float64(time_ns()-wall)/1.e9)
        consumed=empty ? 0.0 : max(0.0,_native_diagnostic_cpu_time()-cpu)
        emit(:info,@sprintf("Done recombining 2D mesh (Wall %gs, CPU %gs)",elapsed,consumed),values)
    end
    hooks=(emit=emit,error=reject,errors=()->errors[],cpu=_native_diagnostic_cpu_time,done=done)
    input_error=function (message,abort=get(options,"General.AbortOnError",2.0))
        values=copy(options)
        values["General.AbortOnError"]=Float64(abort)
        return reject(message,values)
    end
    return (;errors,last_error,events,hooks,begin_mesh,input_error)
end

function _model_import_execute_geo(path,dimension,options,diagnostics;legacy=false,
        initial_model=nothing,initial_mesh_context=nothing,initial_fields=nothing)
    # Session ONELAB bridge: the script's ONELAB exchanges operate on detached
    # scalar mirrors seeded from the parameter space (existing values win the
    # script's DefineNumber defaults, like the upstream localGmsh client);
    # only a successful run adopts the mirrors back into the session records,
    # preserving the file frontends' failure atomicity.
    onelab=_onelab_geo_seed()
    result=execute_geo(path;mesh_dim=dimension,
        _initial_mesh_algorithms=(algorithm_2d=Int(options["Mesh.Algorithm"]),
            algorithm_3d=Int(options["Mesh.Algorithm3D"]),
            recombine_topology_passes=Int(options["Mesh.RecombineOptimizeTopology"]),
            recombine_node_repositioning=Int(options["Mesh.RecombineNodeRepositioning"]),
            recombine_minimum_quality=options["Mesh.RecombineMinimumQuality"],
            native_smoothing=Int(options["Mesh.Smoothing"]),
            second_order_linear=!iszero(options["Mesh.SecondOrderLinear"])),
        _initial_context_lc=NATIVE_CONTEXT_LC[],
        _initial_context_forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[],
        _initial_abort_on_error=options["General.AbortOnError"],
        _initial_mesh_error_callback=diagnostics.input_error,
        _initial_mesh_begin_callback=diagnostics.begin_mesh,
        _initial_abort_option_callback=value->(options["General.AbortOnError"]=value;nothing),
        _initial_numeric_options=options,
        _initial_numeric_option_callback=(family,member,value)->
            _model_import_option!(options,family,member,value),
        _initial_recombine_diagnostics=diagnostics.hooks,
        _initial_model_name=legacy ? nothing : splitext(basename(path))[1],
        _initial_model=initial_model,_initial_mesh_context=initial_mesh_context,
        _initial_fields=initial_fields,
        _initial_onelab_numbers=onelab[1],_initial_onelab_strings=onelab[2])
    _onelab_geo_adopt!(onelab[1],onelab[2],onelab[3],onelab[4])
    return result
end

function _model_import_lookup(receipt,model,caller)
    catalogue=_NativeLookupCatalogue()
    receipt===nothing && return catalogue
    if hasproperty(receipt,:factory_owner)
        authority=model.meshing.mesh_authority
        authority!==nothing && authority.mesh!==nothing || throw(ArgumentError(
            "$caller: imported lookup lacks its actual mesh authority"))
        Model._model_mesh_factory_lookup_binding(receipt,authority.mesh,caller)
    end
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
    return catalogue
end

function _model_import_view_name(section,caller)
    isempty(section.strings) && return section.name
    literal=strip(first(section.strings))
    length(literal)>=2 && first(literal)=='"' && last(literal)=='"' ||
        throw(ArgumentError("$caller: data-section name must be quoted"))
    encoded=length(literal)==2 ? "" :
        SubString(literal,nextind(literal,firstindex(literal)),prevind(literal,lastindex(literal)))
    return Elements._unescape_mixed_name(encoded)
end

function _model_import_views(slot,caller;sections=nothing)
    mesh=slot.mesh
    mesh isa MixedMesh || return _SessionView[]
    views=_SessionView[];positions=Dict{Tuple{String,String},Int}()
    element_tags=reduce(vcat,mesh.entity_data.external_element_tags;init=UInt64[])
    cell_nodes=[copy(block.nodes[:,column]) for block in mesh.blocks for column in axes(block.nodes,2)]
    for section in (sections===nothing ? mesh.data_sections : sections)
        name=_model_import_view_name(section,caller);kind=section.name
        key=(name,kind)
        position=get(positions,key,0)
        if position==0
            push!(views,_SessionView(name));position=length(views);positions[key]=position
        end
        view=views[position]
        step=Int(section.header[1]);components=Int(section.header[2])
        kind=="NodeData" ? (tags=mesh.entity_data.external_node_tags[section.nodes]) :
            (tags=element_tags[section.elements])
        counts=kind=="ElementNodeData" ? Int32[Base.checked_mul(components,Int(count)) for count in section.row_nodes] :
            fill(Int32(components),length(tags))
        if kind=="ElementNodeData" && !section.implicit_nodes
            offset=0
            for (element,count) in zip(section.elements,section.row_nodes)
                nodes=section.nodes[offset+1:offset+count];offset+=count
                nodes==cell_nodes[element] || throw(ArgumentError(
                    "$caller: explicit ElementNodeData nodes must follow literal cell connectivity"))
            end
        end
        storage=_native_view_storage(view,kind,slot,step,components)
        context_key=(step,kind);previous=get(storage.contexts,context_key,nothing)
        old_tags,old_values=previous===nothing ? (UInt64[],Float64[]) : storage.data[context_key][1:2]
        old_counts=previous===nothing ? Int32[] : previous.counts
        rows=_native_view_merge_rows(old_tags,old_counts,old_values,tags,counts,section.values)
        storage.data[context_key]=(rows.tags,rows.values,components)
        partition=length(section.header)>=4 ? Int(section.header[4]) : 0
        partitions=_native_view_partition_plan(previous===nothing ? Set{Int}() : previous.partitions,
            partition;reset=false)
        extent=_native_model_view_extent(tags,_native_view_model_count(slot,kind);
            previous=previous===nothing ? UInt64(0) : previous.extent)
        storage.contexts[context_key]=(model_slot=slot,
            time=isempty(section.reals) ? 0.0 : first(section.reals),counts=rows.counts,
            extent=extent,partitions=partitions,interpolation_msh=_native_view_interpolation_types(slot))
        view.model_data=storage.data;view.model_contexts=storage.contexts
        view.model_step_ranges=storage.ranges;view.model_kind=kind
    end
    return views
end

function _model_import_prepare(path,dimension;legacy=false)
    caller=legacy ? "API.open_geo!" : "API.open"
    options=copy(OPTIONS);diagnostics=_model_import_diagnostics(options)
    extension=lowercase(splitext(path)[2])
    extension==".pos" && return _list_model_open_prepare(path,caller)
    result=if extension==".msh"
        GeoExec._geo_stage_mesh_file(path;
            context_lc=NATIVE_CONTEXT_LC[],forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[])
    elseif extension in (".unv",".mesh",".meshb",".su2",".bdf",".nas",".inp",".neu",".p3d",".mail",".ir3",".diff",".vtk")
        IO._stage_foreign_mesh_file(path;
            context_lc=NATIVE_CONTEXT_LC[],forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[])
    elseif extension==".stl"
        IO._stage_stl_mesh_file(path;
            context_lc=NATIVE_CONTEXT_LC[],forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[])
    else
        _model_import_execute_geo(path,dimension,options,diagnostics;legacy)
    end
    return _model_import_prepare_result(result,path,options,diagnostics;legacy,caller)
end

function _model_import_prepare_result(result,path,options,diagnostics;legacy=false,
        caller="API.open",detached_result=false,view_sections=nothing)
    extension=lowercase(splitext(path)[2])
    model=detached_result ? result.model : deepcopy(result.model)
    authority=model.meshing.mesh_authority
    mesh=authority!==nothing && authority.mesh!==nothing ? authority.mesh :
        result.mesh===nothing ? nothing : _copy_mesh(result.mesh)
    class=if authority!==nothing && authority.mesh!==nothing
        _mixed_rebind_class(_higher_classification_from_data(mesh,caller),mesh;
            node_parametric=Model._model_mesh_runtime_parameters(model,mesh))
    else
        nothing
    end
    # Apply the ordinary cache publisher's record-mirror refresh while the
    # model is detached, so covered records follow the final literal cache.
    if class!==nothing && class.public_tags!==nothing && !isempty(class.public_tags.authority)
        identity=Dict(uid=>uid for uid in class.public_tags.node_tags)
        _dim01_reconcile_records!(model,collect(class.public_tags.authority),mesh,identity;
            runtime_parameters=class.node_parametric)
    end
    if result.transfinite_tri!==nothing
        options["Mesh.TransfiniteTri"]=Float64(result.transfinite_tri)
    end
    model.meshing.transfinite_tri=Int(options["Mesh.TransfiniteTri"])
    if extension==".geo"
        for (name,field) in (("Mesh.Algorithm",:algorithm_2d),("Mesh.Algorithm3D",:algorithm_3d),
                ("Mesh.RecombineOptimizeTopology",:recombine_topology_passes),
                ("Mesh.RecombineNodeRepositioning",:recombine_node_repositioning),
                ("Mesh.RecombineMinimumQuality",:recombine_minimum_quality),
                ("Mesh.Smoothing",:native_smoothing),("Mesh.SecondOrderLinear",:second_order_linear))
            options[name]=Float64(getproperty(model.meshing,field))
        end
    end
    name=legacy ? MODEL_NAME[] : model.name
    slot=_new_slot(name);slot.file_name=String(path);slot.model=model;slot.mesh=mesh;slot.class=class
    slot.node_tag_max=max(result.mesh_node_max,mesh===nothing ? UInt64(0) :
        class===nothing || class.public_tags===nothing ? UInt64(nnodes(mesh)) :
        maximum(class.public_tags.node_tags;init=UInt64(0)))
    slot.element_tag_max=max(result.mesh_element_max,mesh===nothing ? UInt64(0) :
        class===nothing || class.public_tags===nothing ? UInt64(last(_mesh_element_offsets(mesh))) :
        maximum(class.public_tags.element_tags;init=UInt64(0)))
    if legacy
        slot.fields=FIELD_SPECS[];slot.field_tag_max=FIELD_TAG_MAX[]
        slot.background_field=BACKGROUND_FIELD[];slot.boundary_layer_fields=BOUNDARY_LAYER_FIELDS[]
        slot.node_tag_max=max(slot.node_tag_max,NODE_TAG_MAX[])
        slot.element_tag_max=max(slot.element_tag_max,ELEMENT_TAG_MAX[])
        slot.node_counter_wrapped=NODE_COUNTER_WRAPPED[];slot.element_counter_wrapped=ELEMENT_COUNTER_WRAPPED[]
        index=_current_slot_index_locked()
        index==0 || (slot.window_visibility=MODEL_SLOTS[index].window_visibility)
    elseif result.params!==nothing
        slot.fields=Dict(tag=>_ApiField(spec.kind,copy(spec.options),copy(spec.option_order))
            for (tag,spec) in result.params.fields)
        slot.field_tag_max=maximum(keys(slot.fields);init=0)
        slot.background_field=result.params.background_field
        slot.boundary_layer_fields=copy(result.params.boundary_layer_fields)
    end
    lookup=_model_import_lookup(result.native_element_lookup_receipt,model,caller)
    views=legacy ? _SessionView[] : _model_import_views(slot,caller;sections=view_sections)
    view_tags=Int[];next_tag=VIEW_TAG_MAX[]
    for _ in views
        next_tag=Base.checked_add(next_tag,1)
        while haskey(SESSION_VIEWS,next_tag);next_tag=Base.checked_add(next_tag,1);end
        push!(view_tags,next_tag)
    end
    return (;result,slot,options,diagnostics,lookup,views,view_tags,next_tag)
end

function _model_import_rebind_views!(views,source,target)
    source===target && return nothing
    for view in views
        view.model_step_ranges=[_NativeViewStepRange(range.first,range.last,
            range.model_slot===source ? target : range.model_slot,range.components)
            for range in view.model_step_ranges]
        for (key,context) in view.model_contexts
            context.model_slot===source && (view.model_contexts[key]=merge(context,(model_slot=target,)))
        end
    end
    return nothing
end

function _model_import_publish!(plan;legacy=false)
    current=_model_locked();index=_current_slot_index_locked()
    reuse=legacy || isempty(Model.model_entities(current)) &&
        isempty(_native_lookup_node_tags_locked()) && isempty(first(_native_lookup_element_bindings_locked()))
    slot=reuse && index!=0 ? MODEL_SLOTS[index] : plan.slot
    _model_import_rebind_views!(plan.views,plan.slot,slot)
    # All parsing, row transport, lookup validation and allocation plans are
    # complete. Existing held view slots retain their identities.
    _slot_capture_locked!()
    if slot===plan.slot
        push!(MODEL_SLOTS,slot)
    else
        for field in fieldnames(_ModelSlot);setfield!(slot,field,getfield(plan.slot,field));end
    end
    _slot_load_locked!(slot)
    NATIVE_LOOKUP_CATALOGUES[slot]=plan.lookup
    merge!(OPTIONS,plan.options)
    NATIVE_CONTEXT_LC[]=plan.result.native_context_lc
    NATIVE_CONTEXT_FORCED_BBOX[]=plan.result.native_context_forced_bbox
    for (tag,view) in zip(plan.view_tags,plan.views);SESSION_VIEWS[tag]=view;end
    VIEW_TAG_MAX[]=plan.next_tag
    NATIVE_LOOKUP_QUERY_ERRORS[]=plan.diagnostics.errors[]
    NATIVE_MSG_LAST_ERROR[]=plan.diagnostics.last_error[]
    for (kind,message,options) in plan.diagnostics.events
        _native_diagnostic_emit_locked(kind,message;options)
    end
    return nothing
end

function _model_import_filename(filename,caller;geo_only=false)
    filename isa AbstractString && !isempty(filename) && !occursin('\0',filename) ||
        throw(ArgumentError("$caller: filename must be a nonempty NUL-free string"))
    extension=lowercase(splitext(filename)[2])
    extension in (geo_only ? (".geo",) : (".geo",".msh",".pos",".unv",".mesh",".meshb",".stl",".su2",".bdf",".nas",".inp",".neu",".p3d",".mail",".ir3",".diff",".vtk")) || throw(ArgumentError(
        "$caller: supported file extensions are $(geo_only ? ".geo" : ".geo, .msh, .pos, .unv, .mesh, .meshb, .stl, .su2, .bdf, .nas, .inp, .neu, .p3d, .mail, .ir3, .diff and .vtk")"))
    isfile(filename) || throw(ArgumentError("$caller: missing file $filename"))
    return String(filename)
end

"""
    open(filename)

Import a `.geo`, `.msh`, `.unv`, `.mesh`, `.meshb`, `.stl`, `.su2`, `.bdf`,
`.nas`, `.inp`, `.neu`, `.p3d`, `.mail`, `.ir3`, `.diff`, `.vtk` or supported
linear-list `.pos` file into a new current model, reusing an empty
current model. Nonempty prior models and session views remain alive. MSH
imports keep actual UIDs, entity ownership and retained metadata; UNV and
Medit/MESH imports (ASCII `.mesh` and binary `.meshb` alike) publish through
the same channel with their elementary and
physical tags mapped to entity ownership, as do SU2 markers (mapped to
physical groups and names), Nastran property ids (mapped to elementary
entity and physical tags), Abaqus element sets (mapped to elementary and
physical tags plus physical names), Gambit element groups/boundary sets
(mapped the same way) and PLOT3D structured blocks (mapped to
per-dimension entities carrying quads or hexes); MAIL imports publish their
triangles on one implicit surface and IR3 imports publish their
self-describing element records on the tag-column entities; DIFF imports
publish their subdomain-column records on per-subdomain entities like the
pinned reader; VTK imports
publish unstructured-grid/polydata cells on per-dimension entities like
upstream;
STL imports weld through
`read_stl` and publish one discrete Surface; ordinary
model data sections create views bound to the imported model. The model name
is the filename stem unless GEO `SetName` supplies one. Failed preparation
leaves model slots, caches, options and views unchanged; this stricter failure
contract differs from Gmsh's missing-file partial publication.
POS imports preserve named list views and literal time values; text,
interpolation schemes and higher-order list records reject before publication.
"""
function open(filename)
    path=_model_import_filename(filename,"API.open")
    return lock(STATE_LOCK) do
        _model_locked();plan=_model_import_prepare(path,0)
        _model_import_publish!(plan)
        nothing
    end
end
