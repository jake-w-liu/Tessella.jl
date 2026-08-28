function _model_points_of(m::GeoModel,entities,caller::AbstractString)
    (entities isa AbstractVector || entities isa Tuple) || throw(ArgumentError(
        "$caller: PointsOf entities must be a vector or tuple of (dimension, tag) pairs"))
    normalized=NTuple{2,Int}[]
    sizehint!(normalized,length(entities))
    for entry in entities
        entry isa Tuple && length(entry)==2 || throw(ArgumentError(
            "$caller: each PointsOf entity must be a (dimension, tag) pair"))
        dimension=_dimension(entry[1],caller)
        tag=_tag(entry[2],caller,dimension)
        tag>0 || throw(ArgumentError(
            "$caller: PointsOf entity tags must be positive"))
        push!(normalized,(dimension,tag))
    end

    # PointsOf follows recursive entity boundaries; embeddings are not boundary
    # topology and therefore do not participate.
    points=Set{Int}()
    function add_point!(point::Int,owner::AbstractString)
        haskey(m.points,point) || throw(ArgumentError(
            "$caller: $owner references unknown Point[$point]"))
        push!(points,point)
        return nothing
    end
    function add_curve!(curve::Int)
        haskey(m.curves,curve) || throw(ArgumentError(
            "$caller: unknown Curve[$curve]"))
        for point in m.curves[curve]
            add_point!(point,"Curve[$curve]")
        end
        return nothing
    end
    function add_surface!(surface::Int)
        haskey(m.surfaces,surface) || throw(ArgumentError(
            "$caller: unknown Surface[$surface]"))
        for loop in m.surfaces[surface]
            haskey(m.loops,loop) || throw(ArgumentError(
                "$caller: Surface[$surface] references unknown Loop[$loop]"))
            for signed_curve in m.loops[loop]
                add_curve!(abs(signed_curve))
            end
        end
        return nothing
    end
    function add_volume!(volume::Int)
        haskey(m.volumes,volume) || throw(ArgumentError(
            "$caller: unknown Volume[$volume]"))
        shells=m.volumes[volume]
        isempty(shells) && throw(ArgumentError(
            "$caller: Volume[$volume] has no explicit surface-loop topology; " *
            "define it from Surface Loop entities or select explicit Points"))
        for shell in shells
            haskey(m.surface_loops,shell) || throw(ArgumentError(
                "$caller: Volume[$volume] references unknown Surface Loop[$shell]"))
            for signed_surface in m.surface_loops[shell]
                add_surface!(abs(signed_surface))
            end
        end
        return nothing
    end

    for (dimension,tag) in normalized
        if dimension==0
            add_point!(tag,"Point[$tag]")
        elseif dimension==1
            add_curve!(tag)
        elseif dimension==2
            add_surface!(tag)
        else
            add_volume!(tag)
        end
    end
    return sort!(collect(points))
end

@inline function _model_entity_dictionary(m::GeoModel,dimension::Int)
    dimension==0 && return m.points
    dimension==1 && return m.curves
    dimension==2 && return m.surfaces
    return m.volumes
end

function _model_topology_dim_tags(dim_tags,caller::AbstractString)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    normalized=Tuple{Int,Int}[]
    sizehint!(normalized,length(dim_tags))
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        dimension=_dimension(pair[1],caller)
        value=pair[2]
        value isa Integer || throw(ArgumentError(
            "$caller: entity tag must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: entity tag must not be Bool"))
        value!=0 || throw(ArgumentError(
            "$caller: entity tags must be nonzero"))
        (-typemax(Int32)<=value<=typemax(Int32)) || throw(ArgumentError(
            "$caller: entity tag magnitude exceeds Int32"))
        push!(normalized,(dimension,Int(value)))
    end
    return normalized
end

function _model_surface_boundary_curves(
    m::GeoModel,surface::Int,caller::AbstractString;
    orient_holes::Bool)
    haskey(m.surfaces,surface) || throw(ArgumentError(
        "$caller: unknown Surface[$surface]"))
    boundaries=Int[]
    for (loop_index,loop) in pairs(m.surfaces[surface])
        haskey(m.loops,loop) || throw(ArgumentError(
            "$caller: Surface[$surface] references unknown Loop[$loop]"))
        curves=m.loops[loop]
        for signed_curve in curves
            curve=abs(signed_curve)
            haskey(m.curves,curve) || throw(ArgumentError(
                "$caller: Loop[$loop] references unknown Curve[$curve]"))
        end
        if orient_holes && loop_index>1
            for signed_curve in Iterators.reverse(curves)
                push!(boundaries,-signed_curve)
            end
        else
            append!(boundaries,curves)
        end
    end
    return boundaries
end

function _model_direct_boundary(
    m::GeoModel,dimension::Int,signed_tag::Int,caller::AbstractString;
    canonical_orientation::Bool)
    tag=abs(signed_tag)
    if dimension==0
        haskey(m.points,tag) || throw(ArgumentError(
            "$caller: unknown Point[$tag]"))
        return Tuple{Int,Int}[]
    elseif dimension==1
        haskey(m.curves,tag) || throw(ArgumentError(
            "$caller: unknown Curve[$tag]"))
        endpoints=m.curves[tag]
        for point in endpoints
            haskey(m.points,point) || throw(ArgumentError(
                "$caller: Curve[$tag] references unknown Point[$point]"))
        end
        first_point,last_point=signed_tag>0 ? endpoints : reverse(endpoints)
        return Tuple{Int,Int}[(0,first_point),(0,last_point)]
    elseif dimension==2
        curves=_model_surface_boundary_curves(
            m,tag,caller;orient_holes=canonical_orientation)
        return Tuple{Int,Int}[(1,curve) for curve in curves]
    end
    haskey(m.volumes,tag) || throw(ArgumentError(
        "$caller: unknown Volume[$tag]"))
    isempty(m.volumes[tag]) && throw(ArgumentError(
        "$caller: Volume[$tag] has no explicit surface-loop topology; " *
        "define it from Surface Loop entities before querying its boundary"))
    surfaces=_model_volume_boundary_surfaces(
        m,tag,caller;orient_cavities=canonical_orientation)
    return Tuple{Int,Int}[(2,surface) for surface in surfaces]
end

function _model_combined_boundary(boundaries::Vector{Tuple{Int,Int}})
    odd=Dict{Tuple{Int,Int},Int}()
    for (dimension,signed_tag) in boundaries
        key=(dimension,abs(signed_tag))
        if haskey(odd,key)
            delete!(odd,key)
        else
            odd[key]=signed_tag
        end
    end
    combined=Tuple{Int,Int}[
        (dimension,signed_tag)
        for ((dimension,_),signed_tag) in odd]
    return sort!(combined;by=entry->(entry[1],abs(entry[2])))
end

function _model_recursive_boundary(
    m::GeoModel,dimension::Int,signed_tag::Int,caller::AbstractString)
    tag=abs(signed_tag)
    if dimension==0
        haskey(m.points,tag) || throw(ArgumentError(
            "$caller: unknown Point[$tag]"))
        return Tuple{Int,Int}[(0,tag)]
    elseif dimension==1
        return _model_direct_boundary(
            m,dimension,signed_tag,caller;canonical_orientation=true)
    end
    return Tuple{Int,Int}[
        (0,point) for point in _model_points_of(
            m,[(dimension,tag)],caller)]
end

"""
    model_entities(model, dim=-1) -> Vector{Tuple{Int,Int}}

Return detached, sorted `(dimension, tag)` pairs for every native model entity.
`dim=-1` selects all dimensions; dimensions 0 through 3 filter the result. Curve
loops and surface loops are construction records, not model entities.
"""
function model_entities(m::GeoModel,dim=-1)
    selected=_query_dimension(dim,"model_entities")
    result=Tuple{Int,Int}[]
    for dimension in 0:3
        selected==-1 || selected==dimension || continue
        append!(result,((dimension,tag)
                        for tag in keys(_model_entity_dictionary(m,dimension))))
    end
    return sort!(result)
end

"""
    model_dimension(model) -> Int

Return the greatest dimension containing a native model entity, or `-1` when the
model is empty.
"""
function model_dimension(m::GeoModel)
    for dimension in 3:-1:0
        isempty(_model_entity_dictionary(m,dimension)) || return dimension
    end
    return -1
end

"""
    model_boundary(model, dim_tags, combined=true, oriented=false,
                   recursive=false) -> Vector{Tuple{Int,Int}}

Return the explicit topological boundary of the selected signed entities. With
`combined=false`, individual boundaries retain input and incidence order. With
`combined=true`, even incidences cancel and the remaining entries are sorted.
`oriented=true` retains signed Curve and Surface tags. `recursive=true` returns
the dimension-0 closure, including input Points. Queries that require implicit
primitive or Boolean subentities raise an error.
"""
function model_boundary(
    m::GeoModel,dim_tags,combined=true,oriented=false,recursive=false)
    caller="model_boundary"
    combined isa Bool || throw(ArgumentError(
        "$caller: combined must be Bool"))
    oriented isa Bool || throw(ArgumentError(
        "$caller: oriented must be Bool"))
    recursive isa Bool || throw(ArgumentError(
        "$caller: recursive must be Bool"))
    normalized=_model_topology_dim_tags(dim_tags,caller)
    boundaries=Tuple{Int,Int}[]
    for (dimension,signed_tag) in normalized
        current=recursive ?
            _model_recursive_boundary(m,dimension,signed_tag,caller) :
            _model_direct_boundary(
                m,dimension,signed_tag,caller;canonical_orientation=true)
        append!(boundaries,current)
    end
    combined && (boundaries=_model_combined_boundary(boundaries))
    oriented || (boundaries=Tuple{Int,Int}[
        (dimension,abs(tag)) for (dimension,tag) in boundaries])
    return boundaries
end

"""
    model_adjacencies(model, dim, tag) -> (upward, downward)

Return detached topology-only adjacency tags for an existing positive entity.
`upward` is sorted and contains entities of dimension `dim + 1`; `downward`
retains direct boundary order and contains dimension `dim - 1`. Mesh embeddings
do not create topological adjacencies. Downward queries for primitive and Boolean
Volumes raise an error because they have no explicit Surface Loop topology.
"""
function model_adjacencies(m::GeoModel,dim,tag)
    caller="model_adjacencies"
    dimension=_dimension(dim,caller)
    entity_tag=_tag(tag,caller,dimension)
    entity_tag>0 || throw(ArgumentError(
        "$caller: entity tag must be positive"))
    downward=Int[abs(boundary_tag) for (_,boundary_tag) in
        _model_direct_boundary(
            m,dimension,entity_tag,caller;canonical_orientation=true)]
    unique!(downward)
    upward=Int[]
    if dimension==0
        for (curve,endpoints) in m.curves
            entity_tag in endpoints && push!(upward,curve)
        end
    elseif dimension==1
        for surface in keys(m.surfaces)
            curves=_model_surface_boundary_curves(
                m,surface,caller;orient_holes=false)
            any(curve->abs(curve)==entity_tag,curves) && push!(upward,surface)
        end
    elseif dimension==2
        for (volume,shells) in m.volumes
            isempty(shells) && continue
            surfaces=_model_volume_boundary_surfaces(
                m,volume,caller;orient_cavities=false)
            any(surface->abs(surface)==entity_tag,surfaces) && push!(upward,volume)
        end
    end
    return sort!(unique!(upward)),downward
end

function _model_boundary(m::GeoModel,entities,caller::AbstractString;
                         combined::Bool=false,
                         max_entities::Int=typemax(Int))
    query=combined ? "CombinedBoundary" : "Boundary"
    (entities isa AbstractVector || entities isa Tuple) || throw(ArgumentError(
        "$caller: $query entities must be a vector or tuple of " *
        "(dimension, signed tag) pairs"))
    isempty(entities) && throw(ArgumentError(
        "$caller: $query must contain at least one entity"))

    normalized=NTuple{2,Int}[]
    sizehint!(normalized,length(entities))
    source_dimension=nothing
    for entry in entities
        entry isa Tuple && length(entry)==2 || throw(ArgumentError(
            "$caller: each $query entity must be a (dimension, signed tag) pair"))
        dimension=_dimension(entry[1],caller)
        dimension>0 || throw(ArgumentError(
            "$caller: $query input entities must have dimension 1, 2, or 3"))
        value=entry[2]
        value isa Integer || throw(ArgumentError(
            "$caller: $query entity tag must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: $query entity tag must not be Bool"))
        value!=0 || throw(ArgumentError(
            "$caller: $query entity tags must be nonzero"))
        (-typemax(Int32)<=value<=typemax(Int32)) || throw(ArgumentError(
            "$caller: $query entity tag magnitude exceeds Int32"))
        source_dimension===nothing || source_dimension==dimension ||
            throw(ArgumentError(
                "$caller: $query input entities must share one dimension"))
        source_dimension=dimension
        push!(normalized,(dimension,Int(value)))
    end

    boundary=Int[]
    function add_boundary!(tag::Int)
        length(boundary)<max_entities || throw(ArgumentError(
            "$caller: $query expands beyond $max_entities entities"))
        push!(boundary,tag)
        return nothing
    end
    for (dimension,signed_tag) in normalized
        for (_,boundary_tag) in _model_direct_boundary(
                m,dimension,signed_tag,caller;canonical_orientation=false)
            add_boundary!(abs(boundary_tag))
        end
    end
    combined || return boundary

    odd=Set{Int}()
    for tag in boundary
        tag in odd ? delete!(odd,tag) : push!(odd,tag)
    end
    return sort!(collect(odd))
end
