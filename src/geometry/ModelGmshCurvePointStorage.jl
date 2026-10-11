# Native RecursiveIntegration has depth <= 26 and emits two points per leaf.
# Keep small operations in memory and spill larger tables to a temporary mapped
# file. Sampling, insertion order and every Float64 word remain unchanged.
const _GMSH_CURVE_INTEGRATION_CAPACITY=(1<<26)+1

mutable struct _GmshCurveIntegrationPoints <: AbstractVector{_IntegrationPoint}
    memory::Vector{_IntegrationPoint}
    mapping::Vector{_IntegrationPoint}
    stream::Union{Nothing,IOStream}
    path::String
    count::Int
    memory_limit::Int
    closed::Bool
end

function _gmsh_curve_point_store(memory_limit::Int=_DEFAULT_MAX_INTEGRATION_POINTS)
    2<=memory_limit<=_GMSH_CURVE_INTEGRATION_CAPACITY || throw(ArgumentError(
        "native Curve integration memory capacity is out of range"))
    memory=_IntegrationPoint[]
    sizehint!(memory,min(memory_limit,257))
    return _GmshCurveIntegrationPoints(memory,_IntegrationPoint[],nothing,"",0,memory_limit,false)
end

Base.IndexStyle(::Type{_GmshCurveIntegrationPoints})=IndexLinear()
Base.size(points::_GmshCurveIntegrationPoints)=(points.count,)
Base.length(points::_GmshCurveIntegrationPoints)=points.count

@inline function Base.getindex(points::_GmshCurveIntegrationPoints,index::Int)
    @boundscheck checkbounds(points,index)
    points.closed && throw(ArgumentError("native Curve integration storage is closed"))
    return isempty(points.mapping) ? points.memory[index] : points.mapping[index]
end

@inline function Base.setindex!(points::_GmshCurveIntegrationPoints,value::_IntegrationPoint,index::Int)
    @boundscheck checkbounds(points,index)
    points.closed && throw(ArgumentError("native Curve integration storage is closed"))
    if isempty(points.mapping)
        points.memory[index]=value
    else
        points.mapping[index]=value
    end
    return value
end

function _gmsh_curve_spill_points!(points::_GmshCurveIntegrationPoints)
    points.closed && throw(ArgumentError("native Curve integration storage is closed"))
    isempty(points.mapping) || return nothing
    path,stream=mktemp(;cleanup=true)
    points.path=path
    points.stream=stream
    try
        truncate(stream,_GMSH_CURVE_INTEGRATION_CAPACITY*sizeof(_IntegrationPoint))
        seekstart(stream)
        write(stream,reinterpret(UInt8,points.memory))
        flush(stream)
        mapping=Mmap.mmap(stream,Vector{_IntegrationPoint},_GMSH_CURVE_INTEGRATION_CAPACITY,0;grow=false)
        points.mapping=mapping
        points.memory=_IntegrationPoint[]
    catch
        _gmsh_curve_close_points!(points)
        rethrow()
    end
    return nothing
end

function Base.push!(points::_GmshCurveIntegrationPoints,value::_IntegrationPoint)
    points.closed && throw(ArgumentError("native Curve integration storage is closed"))
    points.count<_GMSH_CURVE_INTEGRATION_CAPACITY || throw(ArgumentError(
        "native Curve integration exceeded its recursion bound"))
    if isempty(points.mapping) && points.count==points.memory_limit
        _gmsh_curve_spill_points!(points)
    end
    next=points.count+1
    if isempty(points.mapping)
        push!(points.memory,value)
    else
        points.mapping[next]=value
    end
    points.count=next
    return points
end

# Explicit, idempotent close. No finalizer is registered: GC finalizers must not
# run I/O, locks, or explicit Memory finalization, and the creator is always the
# owner that closes in a finally. The spill file is mktemp(cleanup=true), so an
# abnormal exit still deletes it.
function _gmsh_curve_close_points!(points::_GmshCurveIntegrationPoints)
    points.closed && return nothing
    points.closed=true
    if !isempty(points.mapping)
        # Mmap attaches the unmap/handle finalizer to Array's Memory on both
        # supported Julia versions. Run it before closing/deleting on Windows.
        Base.finalize(points.mapping.ref.mem)
        points.mapping=_IntegrationPoint[]
    end
    if points.stream!==nothing
        close(points.stream)
        points.stream=nothing
    end
    if !isempty(points.path)
        rm(points.path;force=true)
        points.path=""
    end
    points.memory=_IntegrationPoint[]
    points.count=0
    return nothing
end

_gmsh_curve_point_capacity(points::_GmshCurveIntegrationPoints)=_GMSH_CURVE_INTEGRATION_CAPACITY
_gmsh_curve_point_capacity(points::AbstractVector)=_DEFAULT_MAX_INTEGRATION_POINTS
