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
        tag=abs(signed_tag)
        if dimension==1
            haskey(m.curves,tag) || throw(ArgumentError(
                "$caller: unknown Curve[$tag]"))
            endpoints=m.curves[tag]
            for point in endpoints
                haskey(m.points,point) || throw(ArgumentError(
                    "$caller: Curve[$tag] references unknown Point[$point]"))
            end
            first_point,last_point=endpoints
            if signed_tag>0
                add_boundary!(first_point)
                add_boundary!(last_point)
            else
                add_boundary!(last_point)
                add_boundary!(first_point)
            end
        elseif dimension==2
            haskey(m.surfaces,tag) || throw(ArgumentError(
                "$caller: unknown Surface[$tag]"))
            for loop in m.surfaces[tag]
                haskey(m.loops,loop) || throw(ArgumentError(
                    "$caller: Surface[$tag] references unknown Loop[$loop]"))
                for signed_curve in m.loops[loop]
                    curve=abs(signed_curve)
                    haskey(m.curves,curve) || throw(ArgumentError(
                        "$caller: Loop[$loop] references unknown Curve[$curve]"))
                    add_boundary!(curve)
                end
            end
        else
            haskey(m.volumes,tag) || throw(ArgumentError(
                "$caller: unknown Volume[$tag]"))
            shells=m.volumes[tag]
            isempty(shells) && throw(ArgumentError(
                "$caller: Volume[$tag] has no explicit surface-loop topology; " *
                "define it from Surface Loop entities before querying its boundary"))
            for shell in shells
                haskey(m.surface_loops,shell) || throw(ArgumentError(
                    "$caller: Volume[$tag] references unknown Surface Loop[$shell]"))
                for signed_surface in m.surface_loops[shell]
                    surface=abs(signed_surface)
                    haskey(m.surfaces,surface) || throw(ArgumentError(
                        "$caller: Surface Loop[$shell] references unknown Surface[$surface]"))
                    add_boundary!(surface)
                end
            end
        end
    end
    combined || return boundary

    odd=Set{Int}()
    for tag in boundary
        tag in odd ? delete!(odd,tag) : push!(odd,tag)
    end
    return sort!(collect(odd))
end
