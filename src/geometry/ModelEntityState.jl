const _MODEL_DEFAULT_ENTITY_COLOR=(0,0,255,0)

function _model_entity_state_integer(
    value,caller::AbstractString,what::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$caller: $what must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: $what must not be Bool"))
    typemin(Int32)<=value<=typemax(Int32) || throw(ArgumentError(
        "$caller: $what exceeds the Int32 range"))
    return Int(value)
end

function _model_entity_state_targets(
    m::GeoModel,dim_tags,recursive,caller::AbstractString)
    recursive isa Bool || throw(ArgumentError(
        "$caller: recursive must be Bool"))
    requested=_model_removal_dim_tags(dim_tags,caller)
    targets=Set{Tuple{Int,Int}}()
    function add_entity!(dimension::Int,tag::Int)
        entity=(dimension,tag)
        entity in targets && return nothing
        _model_entity_known(m,dimension,tag) || return nothing
        push!(targets,entity)
        if recursive && dimension>0
            for (boundary_dimension,boundary_tag) in
                    _model_removal_boundary(m,dimension,tag,caller)
                add_entity!(boundary_dimension,boundary_tag)
            end
        end
        return nothing
    end
    for (dimension,tag) in requested
        add_entity!(dimension,tag)
    end
    return sort!(collect(targets))
end

"""
    set_entity_visibility!(model, dim_tags, value, recursive=false)

Set an `Int32` visibility value on existing native entities. Missing entities are
unchanged. With `recursive=true`, explicit boundary entities are updated down to
Points; implicit primitive boundaries remain unavailable.
"""
function set_entity_visibility!(m::GeoModel,dim_tags,value,recursive=false)
    caller="set_entity_visibility!"
    visibility=_model_entity_state_integer(value,caller,"value")
    targets=_model_entity_state_targets(m,dim_tags,recursive,caller)
    for entity in targets
        if visibility==1
            delete!(m.entity_visibility,entity)
        else
            m.entity_visibility[entity]=visibility
        end
    end
    return nothing
end

"""Return an existing native entity's visibility, which defaults to `1`."""
function model_entity_visibility(m::GeoModel,dim,tag)
    dimension,entity_tag=_model_metadata_entity(
        m,dim,tag,"model_entity_visibility")
    return get(m.entity_visibility,(dimension,entity_tag),1)
end

function _model_entity_color(r,g,b,a,caller::AbstractString)
    channels=ntuple(4) do index
        name=("red","green","blue","alpha")[index]
        value=(r,g,b,a)[index]
        channel=_model_entity_state_integer(value,caller,"$name channel")
        0<=channel<=255 || throw(ArgumentError(
            "$caller: $name channel must be between 0 and 255"))
        channel
    end
    return channels
end

"""
    set_entity_color!(model, dim_tags, r, g, b, a=255, recursive=false)

Set an RGBA color on existing native entities. Channels are integers from 0 to
255. With `recursive=true`, explicit boundary entities are updated down to Points.
"""
function set_entity_color!(m::GeoModel,dim_tags,r,g,b,a=255,recursive=false)
    caller="set_entity_color!"
    color=_model_entity_color(r,g,b,a,caller)
    targets=_model_entity_state_targets(m,dim_tags,recursive,caller)
    for entity in targets
        if color==_MODEL_DEFAULT_ENTITY_COLOR
            delete!(m.entity_colors,entity)
        else
            m.entity_colors[entity]=color
        end
    end
    return nothing
end

"""Return an existing entity's RGBA color, defaulting to `(0,0,255,0)`."""
function model_entity_color(m::GeoModel,dim,tag)
    dimension,entity_tag=_model_metadata_entity(
        m,dim,tag,"model_entity_color")
    return get(m.entity_colors,(dimension,entity_tag),
               _MODEL_DEFAULT_ENTITY_COLOR)
end

function _model_point_position_update!(updates,coords::Matrix{Float64},row::Int,caller)
    size(coords,1)==3 && 1<=row<=size(coords,2) || throw(ArgumentError(
        "$caller: attached Point vertex does not align with its coordinates"))
    any(binding->binding[1]===coords && binding[2]==row,updates) ||
        push!(updates,(coords,row))
    return nothing
end

# Point mesh vertices retain their entity vector order. Element slots can name
# another vertex, and SetOrderN can move an old endpoint onto its Curve owner.
function _model_point_position_plan(m::GeoModel,tag,x,y,z;first_node=nothing)
    caller="set_point_coordinates!"
    _,point_tag=_model_metadata_entity(m,0,tag,caller)
    coordinate=_finite3(x,y,z,caller)
    record=get(m.discrete,(0,point_tag),nothing)
    record!==nothing && isempty(record.node_tags) && throw(ArgumentError(
        "$caller: discrete point has no node"))
    authority=m.meshing.mesh_authority
    actual=authority===nothing ? nothing : authority.mesh
    owner=(0,Int32(point_tag))
    node=first_node===nothing ? nothing : UInt64(first_node)
    actual_row=0
    if actual!==nothing
        data=actual.entity_data
        if first_node===nothing
            product=get(authority.products,owner,nothing)
            if product!==nothing
                position=findfirst(row->data.node_entities[row]==owner,product.rows)
                if position!==nothing
                    actual_row=Int(product.rows[position])
                    node=data.external_node_tags[actual_row]
                end
            end
        else
            position=findfirst(==(node),data.external_node_tags)
            if position!==nothing
                data.node_entities[position]==owner || throw(ArgumentError(
                    "$caller: attached Point vertex has a different current owner"))
                actual_row=position
            end
        end
    end
    native_record=record===nothing ? get(m.meshing.attached,(0,point_tag),nothing) : record
    if node===nothing && native_record!==nothing && !isempty(native_record.node_tags)
        node=UInt64(first(native_record.node_tags))
    end
    updates=Tuple{Matrix{Float64},Int}[]
    if node!==nothing && native_record!==nothing
        position=findfirst(==(UInt64(node)),native_record.node_tags)
        position===nothing || _model_point_position_update!(updates,native_record.node_coords,
            position,caller)
    end
    if actual_row>0
        _model_point_position_update!(updates,actual.coords,actual_row,caller)
        if authority.carrier!==nothing
            position=findfirst(==(Int32(actual_row)),authority.carrier_rows)
            position===nothing || _model_point_position_update!(updates,
                authority.carrier.coords,position,caller)
        end
        for input in values(authority.native_inputs)
            input isa MixedMesh && input.entity_data!==nothing || continue
            position=findfirst(==(node),input.entity_data.external_node_tags)
            position===nothing || _model_point_position_update!(updates,input.coords,position,caller)
        end
        for input in values(authority.surface_constructor_inputs)
            input isa MixedMesh && input.entity_data!==nothing || continue
            position=findfirst(==(node),input.entity_data.external_node_tags)
            position===nothing || _model_point_position_update!(updates,input.coords,position,caller)
        end
    end
    return (;point_tag,coordinate,node,updates,discrete=record!==nothing)
end

function _model_commit_point_position!(m,plan)
    for (coordinates,row) in plan.updates,axis in 1:3
        coordinates[axis,row]=plan.coordinate[axis]
    end
    if !plan.discrete
        m.points[plan.point_tag]=plan.coordinate
        _native_geo_kernel_update_point!(m,plan.point_tag)
    end
    # A producer relocation receipt describes a particular previous CAD state.
    # A legitimate Point setter supersedes that fact; regeneration can create
    # a new receipt using the current Curve and attached vertex identities.
    hasproperty(m.meshing,:native_curve_endpoints) &&
        delete!(m.meshing.native_curve_endpoints,plan.point_tag)
    return nothing
end

"""
    set_point_coordinates!(model, tag, x, y, z)

Replace one existing Point's finite coordinates and the position of its first
current mesh vertex, when present. A discrete Point must own a mesh node.
Topology, mesh size, names, and other entity metadata retain their identities.
"""
function set_point_coordinates!(m::GeoModel,tag,x,y,z)
    return _model_commit_point_position!(m,_model_point_position_plan(m,tag,x,y,z))
end

function _model_attribute_name(name,caller::AbstractString)
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    result=String(name)
    occursin('\0',result) && throw(ArgumentError(
        "$caller: name must not contain NUL"))
    return result
end

function _model_attribute_values(values,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: values must be a vector or tuple of strings"))
    result=String[]
    sizehint!(result,length(values))
    for value in values
        value isa AbstractString || throw(ArgumentError(
            "$caller: every value must be a string"))
        entry=String(value)
        occursin('\0',entry) && throw(ArgumentError(
            "$caller: values must not contain NUL"))
        push!(result,entry)
    end
    return result
end

"""Set detached string `values` for a model attribute, replacing prior values."""
function set_model_attribute!(m::GeoModel,name,values)
    caller="set_model_attribute!"
    attribute_name=_model_attribute_name(name,caller)
    attribute_values=_model_attribute_values(values,caller)
    m.attributes[attribute_name]=attribute_values
    return nothing
end

"""Return a detached model-attribute value list, or an empty list if absent."""
function model_attribute(m::GeoModel,name)
    attribute_name=_model_attribute_name(name,"model_attribute")
    return copy(get(m.attributes,attribute_name,String[]))
end

"""Return detached model-attribute names in deterministic lexical order."""
model_attribute_names(m::GeoModel)=sort!(collect(keys(m.attributes)))

"""Remove one model attribute and report whether it existed."""
function remove_model_attribute!(m::GeoModel,name)
    attribute_name=_model_attribute_name(name,"remove_model_attribute!")
    return pop!(m.attributes,attribute_name,nothing)!==nothing
end
