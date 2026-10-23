# Native Delaunay/Frontal-Delaunay ordering for Gmsh algorithm ids 5 and 6.
# Coordinates are physical orthonormal plane coordinates supplied by Model.
# The existing Ruppert/legacy insertion path remains a separate refine! mode.
const _GMSH_SURFACE_RADIUS_LIMIT=sqrt(2.0)/2.0

include("Mesh2DGmshNumeric.jl")
include("Mesh2DGmshSurfaceTrace.jl")

# The MTriangle overload of inCircumCircleAniso used by insertAPoint and
# recurFindCavityAniso compares UV metric squared distances strictly. A
# cocircular point is outside; symbolic tie breaking belongs to the generic
# CDT constructor and cannot enlarge this native insertion cavity.
@inline function _gmsh_in_circumcircle(T,t,vid,trace=nothing)
    _is_ghost_tri(T,t) && return false
    first,second,third=_gmsh_surface_triangle_nodes(T,t,trace)
    a=_pt(T,first);b=_pt(T,second);c=_pt(T,third)
    metric=_gmsh_surface_metric(trace,a,b,c)
    center,radius=_gmsh_circumcenter_metric(a,b,c,metric)
    point=_pt(T,vid)
    distance=_gmsh_metric_distance2(center,point,metric)
    return distance<radius
end

# Exact depth-first incoming-edge traversal of recurFindCavityAniso.
# Each explicit frame stores the next incoming edge; visiting a neighbour
# finishes its traversal before the current triangle continues. The generic
# CDT constructor retains its own cavity traversal.
function _gmsh_delaunay_cavity(T::Triangulation,t0,vid,trace=nothing)
    root=Int32(t0)
    cavity=Int32[root]
    visited=Set{Int32}((root,))
    boundary=Tuple{Int32,Int32,Int32}[]
    stack=Tuple{Int32,Int8}[(root,Int8(1))]
    while !isempty(stack)
        triangle,next=last(stack)
        if next>3
            pop!(stack)
            continue
        end
        stack[end]=(triangle,next+Int8(1))
        slot=_gmsh_surface_incoming_slots(T,triangle,trace)[next]
        neighbour=_nbr(T,triangle,slot)
        a,b=_edge(T,triangle,slot)
        if neighbour==0 || _constrained_edge(T,_skey(a,b))
            push!(boundary,(a,b,neighbour))
        elseif !(neighbour in visited)
            if _gmsh_in_circumcircle(T,neighbour,vid,trace)
                push!(visited,neighbour)
                push!(cavity,neighbour)
                push!(stack,(neighbour,Int8(1)))
            else
                push!(boundary,(a,b,neighbour))
            end
        end
    end
    return cavity,boundary
end

@inline function _gmsh_surface_area2(a,b,c)
    return 0.5*abs((b[1]-a[1])*(c[2]-a[2])-
                   (c[1]-a[1])*(b[2]-a[2]))
end

struct _GmshTriangleEntry
    radius::Float64
    key::NTuple{3,Int32}
    priority::NTuple{3,UInt64}
    triangle::Int32
end

mutable struct _GmshTriangleQueue
    heap::Vector{_GmshTriangleEntry}
    pending::Set{NTuple{3,Int32}}
    vertex_numbers::Vector{UInt64}
    vertex_high_water::UInt64
end

struct _GmshSurfaceInsertionProduct
    interior::Vector{Bool}
    vertex_numbers::Vector{UInt64}
    vertex_high_water::UInt64
end

@inline function _gmsh_surface_insertion_result(interior,queue,return_numbers)
    return_numbers || return interior
    return _GmshSurfaceInsertionProduct(
        interior,queue.vertex_numbers,queue.vertex_high_water)
end

@inline function _gmsh_mean_size3(a,b,c)
    total=a+b+c
    return isfinite(total) || isinf(a) || isinf(b) || isinf(c) ?
        total*(1.0/3.0) : a*(1.0/3.0)+b*(1.0/3.0)+c*(1.0/3.0)
end

struct _GmshShellLL{S}
    sizes::S
    propagated::Float64
    background::Float64
end

@inline function (shell::_GmshShellLL)(u::Integer,v::Integer)
    return min(_gmsh_mean_size3(shell.sizes.vsize[u],shell.sizes.vsize[v],shell.propagated),
               _gmsh_mean_size3(shell.sizes.vbgm[u],shell.sizes.vbgm[v],shell.background))
end

function _gmsh_vertex_numbers(T,numbers,high_water)
    values=if numbers===nothing
        UInt64.(1:T.nreal)
    else
        numbers isa AbstractVector{<:Integer} &&
            length(numbers)==T.nreal || throw(ArgumentError(
            "refine!: gmsh_vertex_numbers must contain one integer label per initial vertex"))
        Base.require_one_based_indexing(numbers)
        all(number->!(number isa Bool) && 0<number<=typemax(UInt64),numbers) ||
            throw(ArgumentError("refine!: initial vertex numbers must be positive UInt64 labels"))
        UInt64.(numbers)
    end
    allunique(values) || throw(ArgumentError(
        "refine!: initial vertex numbers must be unique"))
    largest=maximum(values;init=UInt64(0))
    limit=if high_water===nothing
        largest
    else
        high_water isa Integer && !(high_water isa Bool) &&
            largest<=high_water<=typemax(UInt64) || throw(ArgumentError(
            "refine!: gmsh_vertex_high_water must include all initial labels and fit UInt64"))
        UInt64(high_water)
    end
    return values,limit
end

function _gmsh_next_vertex_number!(queue::_GmshTriangleQueue)
    queue.vertex_high_water<typemax(UInt64) || throw(ArgumentError(
        "refine!: native vertex labels overflow UInt64"))
    queue.vertex_high_water+=UInt64(1)
    return queue.vertex_high_water
end

@inline _gmsh_triangle_priority(queue,T,t)=_gmsh_triangle_priority_numbers(queue.vertex_numbers,T,t)
@inline function _gmsh_triangle_priority_numbers(numbers,T,t)
    a=numbers[_vert(T,t,1)]
    b=numbers[_vert(T,t,2)]
    c=numbers[_vert(T,t,3)]
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end

@inline function _gmsh_triangle_key(T,t)
    a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3)
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end

@inline _gmsh_entry_before(a::_GmshTriangleEntry,b::_GmshTriangleEntry)=
    a.radius>b.radius || (a.radius==b.radius && isless(a.priority,b.priority))

function _gmsh_heap_push!(queue::_GmshTriangleQueue,entry::_GmshTriangleEntry)
    heap=queue.heap
    push!(heap,entry)
    slot=length(heap)
    @inbounds while slot>1
        parent=slot>>1
        _gmsh_entry_before(entry,heap[parent]) || break
        heap[slot]=heap[parent]
        slot=parent
    end
    @inbounds heap[slot]=entry
    return nothing
end

function _gmsh_heap_pop!(queue::_GmshTriangleQueue)
    heap=queue.heap
    isempty(heap) && return nothing
    first=heap[1]
    last=pop!(heap)
    if !isempty(heap)
        slot=1
        @inbounds while 2slot<=length(heap)
            child=2slot
            child<length(heap) && _gmsh_entry_before(heap[child+1],heap[child]) &&
                (child+=1)
            _gmsh_entry_before(heap[child],last) || break
            heap[slot]=heap[child]
            slot=child
        end
        @inbounds heap[slot]=last
    end
    delete!(queue.pending,first.key)
    return first
end

@inline function _gmsh_triangle_radius(T,t,gsizes,extend_boundary)
    # Only a full trace carries persistent birth radii; a metric-only carrier
    # keeps the untraced scalar radius computation.
    gsizes.trace isa _GmshSurfaceTrace &&
        return gsizes.trace.radii[_gmsh_triangle_key(T,t)]
    a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3)
    background=_gmsh_mean_size3(gsizes.vbgm[a],gsizes.vbgm[b],gsizes.vbgm[c])
    lc=extend_boundary ? min(background,
        _gmsh_mean_size3(gsizes.vsize[a],gsizes.vsize[b],gsizes.vsize[c])) : background
    first=_gmsh_flat_xyz(_pt(T,a))
    center=_gmsh_circumcenter_xyz(first,_gmsh_flat_xyz(_pt(T,b)),
                                _gmsh_flat_xyz(_pt(T,c)))
    return _gmsh_norm3(_gmsh_sub3(first,center))/lc
end

# Gmsh scans the incoming edges (c,a), (a,b), (b,c). Tessella neighbours
# are stored opposite local vertices, so this is the slot order 2,3,1.
function _gmsh_active_edge(T,t,interior,skipped,gsizes,extend_boundary)
    @inbounds for slot in _gmsh_surface_incoming_slots(T,t,gsizes.trace)
        neighbour=_nbr(T,t,slot)
        if neighbour==0 || !T.alive[neighbour] ||
           neighbour>length(interior) || !interior[neighbour] ||
           _is_ghost_tri(T,neighbour)
            return slot
        end
        _gmsh_triangle_key(T,neighbour) in skipped && continue
        radius=_gmsh_triangle_radius(T,neighbour,gsizes,extend_boundary)
        0<radius<_GMSH_SURFACE_RADIUS_LIMIT && return slot
    end
    return 0
end

function _gmsh_queue_triangle!(queue,T,t,interior,skipped,gsizes,
                               algorithm,extend_boundary)
    (t!=0 && t<=length(T.alive) && T.alive[t] &&
     t<=length(interior) && interior[t] && !_is_ghost_tri(T,t)) || return nothing
    key=_gmsh_triangle_key(T,t)
    (key in skipped || key in queue.pending) && return nothing
    radius=_gmsh_triangle_radius(T,t,gsizes,extend_boundary)
    isnan(radius) && throw(ErrorException(
        "refine!: native normalized triangle radius is NaN"))
    if algorithm===:frontal
        radius>_GMSH_SURFACE_RADIUS_LIMIT || return nothing
        _gmsh_active_edge(T,t,interior,skipped,gsizes,extend_boundary)!=0 ||
            return nothing
    else
        radius>=_GMSH_SURFACE_RADIUS_LIMIT || return nothing
    end
    push!(queue.pending,key)
    _gmsh_heap_push!(queue,_GmshTriangleEntry(
        radius,key,_gmsh_triangle_priority(queue,T,t),Int32(t)))
    return nothing
end

function _gmsh_queue_init(T,interior,skipped,gsizes,algorithm,extend_boundary,
                          numbers,high_water)
    queue=_GmshTriangleQueue(_GmshTriangleEntry[],Set{NTuple{3,Int32}}(),
                             numbers,high_water)
    sizehint!(queue.heap,length(T.alive))
    sizehint!(queue.pending,length(T.alive))
    for t in eachindex(T.alive)
        _gmsh_queue_triangle!(queue,T,t,interior,skipped,gsizes,
                              algorithm,extend_boundary)
    end
    return queue
end

function _gmsh_queue_select!(queue,T,interior,skipped,gsizes,
                             algorithm,extend_boundary)
    while true
        entry=_gmsh_heap_pop!(queue)
        entry===nothing && return (Int32(0),0)
        t=entry.triangle
        (T.alive[t] && t<=length(interior) && interior[t] &&
         !_is_ghost_tri(T,t) && _gmsh_triangle_key(T,t)==entry.key &&
         !(entry.key in skipped)) || continue
        if algorithm===:frontal
            slot=_gmsh_active_edge(T,t,interior,skipped,gsizes,extend_boundary)
            slot==0 && continue
            return (t,slot)
        end
        return (t,0)
    end
end

# optimalPointFrontal supplies the distance along the Voronoi edge. The B
# variant then intersects the normal circle with the CAD plane from angle
# zero, placing the result toward the triangle's opposite vertex. These
# directions differ for an obtuse triangle whose circumcenter is outside.
function _gmsh_frontal_point(T,t,slot,gsizes,extend_boundary)
    a,b,c=_gmsh_surface_triangle_nodes(T,t,gsizes.trace)
    metric=_gmsh_surface_metric(gsizes.trace,_pt(T,a),_pt(T,b),_pt(T,c))
    center,_=_gmsh_circumcenter_metric(_pt(T,a),_pt(T,b),_pt(T,c),
                                     metric)
    u,v=_edge(T,t,slot)
    pu=_pt(T,u);pv=_pt(T,v)
    midpoint=(0.5*(pu[1]+pv[1]),0.5*(pu[2]+pv[2]))
    dx=center[1]-midpoint[1];dy=center[2]-midpoint[2]
    direction_norm=sqrt(dx*dx+dy*dy)
    q=sqrt(dx*metric[1]*dx+2.0*dx*metric[2]*dy+dy*metric[3]*dy)
    (isfinite(q) && q>0) || return nothing
    dirx=dx/direction_norm;diry=dy/direction_norm
    ratio=sqrt(dirx*dirx*metric[1]+2.0*diry*dirx*metric[2]+diry*diry*metric[3])
    background=0.5*(gsizes.vbgm[u]+gsizes.vbgm[v])
    lc=extend_boundary ? min(background,
        0.5*(gsizes.vsize[u]+gsizes.vsize[v])) : background
    distance=min(lc*sqrt(3.0)*0.5,q)
    fast=(midpoint[1]+distance*dirx/ratio,
          midpoint[2]+distance*diry/ratio)
    opposite=_pt(T,_vert(T,t,slot))
    ex=pv[1]-pu[1];ey=pv[2]-pu[2]
    cross=ex*(opposite[2]-midpoint[2])-ey*(opposite[1]-midpoint[1])
    # optimalPointFrontalB retains its fast point when the cross-product
    # normal is below the pinned physical 1e-12 threshold.
    if abs(cross)>=1e-12
        tangent=_gmsh_normalize3((-cross*ey,cross*ex,0.0))
        return _gmsh_plane_circle_intersection(fast,midpoint,tangent,distance)
    end
    return fast
end

@inline function _gmsh_contains_surface_point(T,t,point)
    a=_pt(T,_vert(T,t,1));b=_pt(T,_vert(T,t,2));c=_pt(T,_vert(T,t,3))
    w=_gmsh_barycentric2(a,b,c,point)
    tolerance=1e-8
    return w[2]>=-tolerance && w[3]>=-tolerance &&
           w[2]<=1.0+tolerance && w[3]<=1.0+tolerance && w[1]>-tolerance
end

@inline function _gmsh_intersect_segments(a,b,c,d)
    orient2(a,b,c)*orient2(a,b,d)>0 && return false
    orient2(c,d,a)*orient2(c,d,b)>0 && return false
    return true
end

# search4Triangle uses invMapUV's parametric tolerance before walking from
# the selected triangle. Near-boundary attempts can therefore allocate a
# native birth label before insertVertexB rejects them. The generic exact
# CDT locator cannot decide this operation's allocation order.
function _gmsh_locate_surface(T,point,interior,start)
    triangle=Int32(start)
    iterations=0
    maximum=length(T.alive)
    while triangle!=0 && triangle<=length(interior) && interior[triangle] &&
          T.alive[triangle] && !_is_ghost_tri(T,triangle)
        _gmsh_contains_surface_point(T,triangle,point) && return triangle
        a=_pt(T,_vert(T,triangle,1))
        b=_pt(T,_vert(T,triangle,2))
        c=_pt(T,_vert(T,triangle,3))
        center=((a[1]+b[1]+c[1])*(1.0/3.0),
                (a[2]+b[2]+c[2])*(1.0/3.0))
        next=Int32(0)
        for slot in (2,3,1)
            u,v=_edge(T,triangle,slot)
            if _gmsh_intersect_segments(_pt(T,u),_pt(T,v),point,center)
                next=_nbr(T,triangle,slot)
                break
            end
        end
        triangle=next
        iterations>maximum && return Int32(0)
        iterations+=1
    end
    return Int32(0)
end

# Updating only the replacement cavity and its outside neighbours gives
# the same active-front admission as the reference set insertion. Heap
# entries carry immutable vertex triples so reused triangle slots are safe.
function _gmsh_queue_update!(queue,T,newtris,interior,skipped,gsizes,
                             algorithm,extend_boundary)
    for t in newtris
        _gmsh_queue_triangle!(queue,T,t,interior,skipped,gsizes,
                              algorithm,extend_boundary)
        algorithm===:frontal || continue
        for slot in 1:3
            _gmsh_queue_triangle!(queue,T,_nbr(T,t,slot),interior,skipped,
                                  gsizes,algorithm,extend_boundary)
        end
    end
    return nothing
end
