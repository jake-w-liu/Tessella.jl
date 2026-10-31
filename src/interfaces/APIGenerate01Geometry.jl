# The public node producer and Tessella's established geometry view have
# separate arithmetic contracts. Bind each geometry sample at its actual birth
# identity; never infer these bindings from coordinate coincidence.
mutable struct _Dim01GeometryBirth
    nodes::Dict{UInt64,_Dim01Node}
    parameters::Dict{Int,Vector{Float64}}
end

function _dim01_geometry_birth(previous)
    birth=_Dim01GeometryBirth(Dict{UInt64,_Dim01Node}(),Dict{Int,Vector{Float64}}())
    previous===nothing && return birth
    previous.mesh!==nothing && previous.carrier isa Union{Mesh,MixedMesh} || return birth
    data=previous.carrier isa MixedMesh ? previous.carrier.entity_data : nothing
    actual_data=previous.mesh.entity_data
    actual_data===nothing && return birth
    for (slot,row) in enumerate(previous.carrier_rows)
        tag=previous.mesh.entity_data.external_node_tags[row]
        owner=data===nothing ? actual_data.node_entities[row] : data.node_entities[slot]
        parameter=data===nothing ? actual_data.node_parametric[row] : data.node_parametric[slot]
        birth.nodes[tag]=_Dim01Node(tag,ntuple(axis->previous.carrier.coords[axis,slot],3),
            owner,parameter===nothing ? nothing : copy(parameter))
    end
    return birth
end

function _dim01_geometry_curve_parameters!(birth,m,curve,actual,options,caller)
    # A stateful callback or opaque background evaluator is called solely by
    # the native producer. Its samples are also the geometry view's samples.
    result=actual
    if options.callback===nothing && options.inner_field===nothing &&
            haskey(m.curves,curve) && !haskey(m.discrete,(1,curve)) &&
            !(Model._curve_type(m,curve) in (:tensor_bezier,:tensor_bspline)) &&
            !haskey(m.meshing.extrude,(1,curve)) &&
            Model._model_curve_degenerate_kind(m,curve) == -1
        constraint=get(m.periodic,(1,curve),nothing)
        if constraint!==nothing
            master=abs(constraint.master_entity)
            master_parameters=get(birth.parameters,master,nothing)
            master_parameters===nothing || (result=Model._model_curve_periodic_params(
                m,curve,master,master_parameters,constraint.reversed,caller))
        end
    end
    # A differing birth count has no ordinal identity correspondence. Such a
    # producer retains its actual vector rather than inventing a geometry map.
    length(result)==length(actual) || (result=actual)
    birth.parameters[curve]=copy(result)
    return result
end

# copyMesh transforms existing master mesh vertices, projects each onto the
# slave, and stores the slave parameter. Public IDs bind these samples even
# when two distinct vertices have identical coordinates.
function _dim01_native_periodic_parameters(assembly,m,curve,params,curve_births)
    constraint=get(m.periodic,(1,curve),nothing)
    constraint!==nothing && constraint.affine!==nothing || return params
    master=Int(abs(constraint.master_entity))
    adapter=Model._model_high_order_curve_geometry(m,curve,assembly.caller)
    adapter===nothing && throw(ArgumentError(
        "$(assembly.caller): periodic Curve[$curve] lacks native projection geometry"))
    master_tags=curve_births===nothing ? nothing : get(curve_births,(1,Int32(master)),nothing)
    rows=Int32[]
    if master_tags!==nothing
        for tag in master_tags
            row=get(assembly.node_positions,tag,Int32(0))
            row>0 || throw(ArgumentError("$(assembly.caller): periodic master loses vertex $tag"))
            push!(rows,row)
        end
    else
        # A skipped master retains its actual samples. Its stored parameters
        # identify interior vertices; endpoint Point identity is explicit.
        endpoints=get(m.curves,master,nothing)
        endpoints===nothing && throw(ArgumentError(
            "$(assembly.caller): retained periodic master needs explicit endpoint identities"))
        master_params=get(m.curve_params,master,nothing)
        master_params===nothing && throw(ArgumentError(
            "$(assembly.caller): retained periodic master lacks stored parameters"))
        point_rows=Dict{Int32,Vector{Int32}}()
        parameter_rows=Dict{Float64,Vector{Int32}}()
        for (row,node) in enumerate(assembly.nodes)
            if node.owner[1]==0
                push!(get!(()->Int32[],point_rows,node.owner[2]),Int32(row))
            elseif node.owner==(1,Int32(master)) && node.parameter!==nothing && length(node.parameter)==1
                push!(get!(()->Int32[],parameter_rows,node.parameter[1]),Int32(row))
            end
        end
        for (ordinal,parameter) in enumerate(master_params)
            candidates=ordinal==1 ? get(point_rows,Int32(endpoints[1]),Int32[]) :
                ordinal==length(master_params) ? get(point_rows,Int32(endpoints[2]),Int32[]) :
                get(parameter_rows,parameter,Int32[])
            length(candidates)==1 || throw(ArgumentError(
                "$(assembly.caller): retained periodic master has ambiguous vertex identity"))
            push!(rows,only(candidates))
        end
    end
    length(rows)==length(params) || throw(ArgumentError(
        "$(assembly.caller): periodic master/slave vertex counts differ"))
    constraint.reversed && reverse!(rows)
    result=copy(params)
    for ordinal in 2:length(rows)-1
        source=assembly.nodes[rows[ordinal]].point
        target=Model._model_native_curve_transform(source,constraint.affine)
        result[ordinal]=Model._model_native_curve_copy_parameter(m,curve,adapter,target,assembly.caller)
    end
    return result
end

function _dim01_geometry_curve_node!(birth,assembly,m,curve,node,parameter)
    source=assembly.nodes[node]
    point=Model._curve_type(m,curve) in (:tensor_bezier,:tensor_bspline) ?
        source.point : Model._model_curve_part_point(m,curve,parameter,assembly.caller;
        _native_producer=false,_parameter_source=birth.parameters)
    birth.nodes[source.tag]=_Dim01Node(source.tag,point,source.owner,Float64[parameter])
    return nothing
end

function _dim01_geometry_support!(birth,assembly,m,node,key;linear=false)
    source=assembly.nodes[node]
    owner=source.owner
    if (owner[1]==1 && Model._curve_type(m,Int(owner[2])) in (:tensor_bezier,:tensor_bspline)) ||
            (owner[1]==2 && Model._model_nurbs_surface(get(m.surface_geometry,Int(owner[2]),nothing)))
        # Newly constructed tensor geometry has no inherited alternate numeric
        # producer. Its existing public birth is the geometry carrier's birth.
        birth.nodes[source.tag]=_Dim01Node(source.tag,source.point,owner,source.parameter)
        return nothing
    end
    native_curve=source.owner[1]==1 && haskey(m.curves,Int(source.owner[2]))
    native_curve || any(term->haskey(birth.nodes,assembly.nodes[term[1]].tag),key) || return nothing
    point=ntuple(axis->sum(weight*get(birth.nodes,assembly.nodes[index].tag,
        assembly.nodes[index]).point[axis] for (index,weight) in key),3)
    parameter=source.parameter
    if !linear && native_curve
        curve=Int(source.owner[2])
        values=Model.model_parametrization(m,1,curve,collect(point))
        parameter=Float64[values[1]]
        point=Model._model_curve_part_point(m,curve,values[1],assembly.caller;
            _native_producer=false,_parameter_source=birth.parameters)
    end
    birth.nodes[source.tag]=_Dim01Node(source.tag,point,source.owner,parameter)
    return nothing
end

# MSH4 writes Curve begin/end with +1/-1 roles. The established geometry
# descriptor instead follows getBoundaryTags' unsigned ordered vertex list.
# Only the current known Curve's explicit topology can bind these two views;
# raw/discrete descriptors and the signed actual mesh remain untouched.
function _dim01_geometry_entities(mesh,m)
    entities=mesh.entity_data.entities
    replacement=nothing
    for (key,record) in entities
        record.dim==1 || continue
        endpoints=get(m.curves,Int(record.tag),nothing)
        endpoints===nothing && continue
        first_point,last_point=endpoints
        signed=Int32[first_point,-last_point]
        record.boundaries==signed || continue
        public=Int32[first_point,last_point]
        replacement===nothing && (replacement=copy(entities))
        replacement[key]=Elements._native_mixed_entity(record.dim,record.tag,record.bbox;
            physical_tags=record.physical_tags,boundaries=public,
            embedded_curves=record.embedded_curves,partitions=record.partitions,parent=record.parent)
    end
    return replacement===nothing ? entities : replacement
end

function _dim01_geometry_carrier(plan,birth)
    mesh=plan.mesh
    entities=_dim01_geometry_entities(mesh,plan.model)
    isempty(birth.nodes) && entities===mesh.entity_data.entities && return mesh
    data=mesh.entity_data
    coords=copy(mesh.coords)
    parameters=deepcopy(data.node_parametric)
    rows=Dict(tag=>row for (row,tag) in enumerate(data.external_node_tags))
    for (tag,node) in birth.nodes
        row=get(rows,get(plan.node_tag_map,tag,UInt64(0)),0)
        row>0 && data.node_entities[row]==node.owner || continue
        coords[:,row].=node.point
        parameters[row]=Model._model_mesh_projected_parameter(node.owner,node.parameter)
    end
    geometry_data=Elements._native_mixed_entity_data(entities;node_entities=data.node_entities,
        node_parametric=parameters,external_node_tags=data.external_node_tags,
        block_entities=data.block_entities,external_element_tags=data.external_element_tags)
    return MixedMesh(coords,mesh.blocks;physical_names=mesh.physical_names,
        entity_data=geometry_data,elementary_entities=mesh.elementary_entities,
        periodic_links=mesh.periodic_links,ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
end
