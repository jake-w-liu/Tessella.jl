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
