struct _ConstantSurfacePointSize
    value::Float64
end

@inline (field::_ConstantSurfacePointSize)(x,y)=field.value

struct _InterpolatedSurfacePointSize
    normalized::PostViewField
    minimum::Float64
    maximum::Float64
    span::Float64
    surface::Int
    caller::String
end

@inline function (field::_InterpolatedSurfacePointSize)(x,y)
    normalized=field_value(field.normalized,x,y,0.0)
    tolerance=2*field.normalized.reference_tolerance+128eps(Float64)
    (isfinite(normalized) && -tolerance<=normalized<=1+tolerance) ||
        throw(ErrorException(
            "$(field.caller): Point-size query ($x,$y) is outside the " *
            "initial triangulation of Surface[$(field.surface)]"))
    fraction=clamp(normalized,0.0,1.0)
    mesh_size=fraction<=0.5 ? muladd(fraction,field.span,field.minimum) :
                             muladd(-(1-fraction),field.span,field.maximum)
    (isfinite(mesh_size) && mesh_size>0) || throw(ErrorException(
        "$(field.caller): Surface[$(field.surface)] Point-size interpolation " *
        "is not finite and positive at ($x,$y)"))
    return mesh_size
end

@inline _surface_size_key(x::Float64,y::Float64)=
    (iszero(x) ? 0.0 : x,iszero(y) ? 0.0 : y)

function _surface_point_size_field(T,xs,ys,mesh_sizes,t::Int,
                                   caller::AbstractString)
    length(xs)==length(ys)==length(mesh_sizes) || throw(ErrorException(
        "$caller: Surface[$t] Point-size data is inconsistent"))
    isempty(mesh_sizes) && throw(ErrorException(
        "$caller: Surface[$t] has no Point-size data"))
    all(size->isfinite(size) && size>0,mesh_sizes) || throw(ErrorException(
        "$caller: Surface[$t] has a Point size that is not finite and positive"))
    minimum_size,maximum_size=extrema(mesh_sizes)
    minimum_size==maximum_size && return _ConstantSurfacePointSize(minimum_size)

    initial=to_mesh(T;interior=classify_interior(T))
    ntris(initial)>0 || throw(ErrorException(
        "$caller: Surface[$t] has no initial interior triangles"))
    source_sizes=Dict{NTuple{2,Float64},Float64}()
    for vertex in eachindex(xs)
        key=_surface_size_key(xs[vertex],ys[vertex])
        source_sizes[key]=min(get(source_sizes,key,Inf),mesh_sizes[vertex])
    end
    nodal_sizes=Vector{Float64}(undef,nnodes(initial))
    @inbounds for node in 1:nnodes(initial)
        key=_surface_size_key(initial.coords[1,node],initial.coords[2,node])
        haskey(source_sizes,key) || throw(ErrorException(
            "$caller: Surface[$t] initial node $node has no Point-size value"))
        nodal_sizes[node]=source_sizes[key]
    end

    span=maximum_size-minimum_size
    normalized=(nodal_sizes .- minimum_size)./span
    interpolation=PostViewField(
        initial.coords,normalized;triangles=initial.tris,crop_negative=false,
        use_closest=false,max_nodes=nnodes(initial),max_elements=ntris(initial))
    return _InterpolatedSurfacePointSize(
        interpolation,minimum_size,maximum_size,span,t,String(caller))
end
