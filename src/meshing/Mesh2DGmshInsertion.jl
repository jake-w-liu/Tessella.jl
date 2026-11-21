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
    # The cavity is its own visited list: every visited cell is a member, and
    # a cavity has a handful of cells, so a hash set per insertion cost more
    # than the search it saved.
    cavity=Int32[root]
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
        elseif !(neighbour in cavity)
            if _gmsh_in_circumcircle(T,neighbour,vid,trace)
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

# Scalar local bound folded into a 2×2 SPD metric: `min(lc_bgm, lc_vsize)`
# composed in metric space — the local `h` caps every directional size, so
# each eigendirection's eigenvalue floors at `1/h²` (the canonical metric
# intersection under simultaneous diagonalization).
@inline function _gmsh_metric_floor(metric,h)
    a,b,d=metric
    (h>0 && isfinite(h)) || return metric
    tr=0.5*(a+d)
    det_half=sqrt(0.25*(a-d)*(a-d)+b*b)
    lmin=tr-det_half;lmax=tr+det_half
    floor_=1.0/(h*h)
    lmin>=floor_ && return metric
    if det_half<=1e-14*max(tr,1.0)
        return (max(a,floor_),0.0,max(d,floor_))
    end
    # Eigenvector of lmin: both null-space candidates (b, lmin-a) and
    # (lmin-d, -b) satisfy the eigen equation; keep the better-conditioned
    # one (the first degenerates when lmin==a, the second when lmin==d).
    vx=b;vy=lmin-a
    wx=lmin-d;wy=-b
    if wx*wx+wy*wy>vx*vx+vy*vy
        vx=wx;vy=wy
    end
    nv=sqrt(vx*vx+vy*vy)
    vx/=nv;vy/=nv
    ux=-vy;uy=vx
    return (floor_*vx*vx+lmax*ux*ux,
            floor_*vx*vy+lmax*ux*uy,
            floor_*vy*vy+lmax*uy*uy)
end

# The same scalar-bound fold for a parameter domain whose UV coordinates are
# not world units: `h` is a world length, so the bound pulls back to
# `N = G/h²` (first fundamental form over the bound), and the fold is the
# simultaneous-diagonalization intersection `max(M, N)` — each generalized
# eigendirection of `(M, N)` floors at 1. `M' = M + (1−λmin)·(Nv)(Nv)ᵀ` for
# `v` N-normalized realizes exactly that: `Mv = λmin·Nv` becomes `M'v = Nv`
# while the N-orthogonal eigendirection keeps `λmax` — exact when `M ∥ N`
# (isotropic background), the canonical intersection approximation otherwise.
# Unlike `_gmsh_metric_floor`, which retains a strong `lmax` below the bound
# (the scalar-radius approximation upstream applies beside the metric), this
# fold returns `N` itself when `λmax ≤ 1`: `G/h²` is a genuine metric cap, so
# every generalized direction floors at the world-space bound.
function _gmsh_metric_ground_floor(metric,ground,h)
    a,b,d=metric
    (h>0 && isfinite(h)) || return metric
    inv2=1.0/(h*h)
    n11=ground[1]*inv2;n12=ground[2]*inv2;n22=ground[3]*inv2
    detn=n11*n22-n12*n12
    detn>0 || return _gmsh_metric_floor(metric,h)
    # det(M-λN)=0 → λ²·detN − λ·(a·n22 + d·n11 − 2b·n12) + detM = 0
    cB=a*n22+d*n11-2.0*b*n12
    detm=a*d-b*b
    disc=cB*cB-4.0*detn*detm
    # det(M-λN) has nonnegative roots for this SPD pair; a tiny negative
    # discriminant is rounding on a repeated eigenvalue, which is the common
    # case when `metric` and `ground/h²` are both proportional to the first
    # fundamental form. Clamping keeps the pulled-back bound; only a grossly
    # negative discriminant indicates a degenerate input and falls back.
    disc < -64.0*eps(Float64)*max(cB*cB,4.0*detn*abs(detm),1.0) &&
        return _gmsh_metric_floor(metric,h)
    r=sqrt(max(disc,0.0))
    lmin=(cB-r)/(2.0*detn);lmax=(cB+r)/(2.0*detn)
    lmax<=1.0 && return (n11,n12,n22)
    lmin>=1.0 && return metric
    # Generalized eigenvector of the weak direction: (M − λmin·N)x = 0,
    # null-space candidates (b − λmin·n12, λmin·n11 − a) and
    # (λmin·n22 − d, b − λmin·n12) — keep the better-conditioned one.
    vx=b-lmin*n12;vy=lmin*n11-a
    wx=lmin*n22-d;wy=b-lmin*n12
    wx*wx+wy*wy>vx*vx+vy*vy && ((vx,vy)=(wx,wy))
    nv=sqrt(vx*vx+vy*vy)
    nv>0 || return _gmsh_metric_floor(metric,h)
    vx/=nv;vy/=nv
    nx=n11*vx+n12*vy;ny=n12*vx+n22*vy
    nvN=sqrt(vx*nx+vy*ny)
    nvN>0 || return _gmsh_metric_floor(metric,h)
    nx/=nvN;ny/=nvN
    s=1.0-lmin
    return (a+s*nx*nx,b+s*nx*ny,d+s*ny*ny)
end

@inline function _gmsh_triangle_radius(T,t,gsizes,extend_boundary)
    # Only a full trace carries persistent birth radii; a metric-only carrier
    # keeps the untraced scalar radius computation.
    gsizes.trace isa _GmshSurfaceTrace &&
        return gsizes.trace.radii[_gmsh_triangle_key(T,t)]
    a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3)
    if gsizes.trace isa _GmshSurfaceMetricCarrier
        # Directional queue radius: the metric-normalized circumradius
        # `sqrt(R'M R)` measures the cell in units of `h`, so weak-direction
        # elongation no longer fires the scalar min-direction bound.
        pa=_pt(T,a);pb=_pt(T,b);pc=_pt(T,c)
        if gsizes.trace.cartesian_close
            # Isotropic background on a warped chart: upstream's default
            # `MTri3` radius (`radiusNorm == 2`) is the Cartesian
            # circumradius over the scalar size, exactly as on a plane. A
            # centroid metric cannot stand in for it on a cell that touches
            # a pole, where one parameter direction has no world length.
            wa=_gmsh_surface_world(gsizes.trace,T,a)
            wb=_gmsh_surface_world(gsizes.trace,T,b)
            wc=_gmsh_surface_world(gsizes.trace,T,c)
            # A cell resting on a collapsed boundary row (a pole/apex edge:
            # two parameter representatives of one world point) has no world
            # extent and is dropped by the weld, so it is neither refined
            # nor an active front.
            ((wa==wb && _constrained_edge(T,_skey(a,b))) ||
             (wb==wc && _constrained_edge(T,_skey(b,c))) ||
             (wc==wa && _constrained_edge(T,_skey(c,a)))) && return 0.0
            background=_gmsh_mean_size3(gsizes.vbgm[a],gsizes.vbgm[b],gsizes.vbgm[c])
            lc=extend_boundary ? min(background,
                _gmsh_mean_size3(gsizes.vsize[a],gsizes.vsize[b],gsizes.vsize[c])) : background
            # Two seam representatives of one world point across the sheet:
            # `circumCenterXY` reports collinear input with its -99999
            # sentinel, which queues the cell ahead of every real one.
            (wa==wb || wb==wc || wc==wa) && return 99999.0*sqrt(2.0)/lc
            radius=_gmsh_cartesian_radius(wa,wb,wc,lc)
            # Second pass only (`_gmsh_surface_second_pass!`): a cell the
            # Cartesian radius calls resolved is asked a second question. A
            # parameter edge that runs more than half a turn around the
            # surface has a small flat world image however large the region
            # it covers (the far side of a sphere seen from one seam
            # column), so that cell is queued like the collinear ones
            # whatever the size. A cell the radius already queues keeps it.
            (gsizes.trace.wrap[] && radius<_GMSH_SURFACE_RADIUS_LIMIT &&
             _gmsh_cell_wrapped(gsizes.trace,pa,pb,pc,wa,wb,wc)) &&
                return 99999.0*sqrt(2.0)/lc
            return radius
        end
        metric=_gmsh_surface_metric(gsizes.trace,pa,pb,pc)
        if extend_boundary
            vs=_gmsh_mean_size3(gsizes.vsize[a],gsizes.vsize[b],gsizes.vsize[c])
            metric=_gmsh_surface_metric_floor(gsizes.trace,metric,pa,pb,pc,vs)
        end
        _,r2=_gmsh_circumcenter_metric(pa,pb,pc,metric)
        return sqrt(max(r2,0.0))
    end
    background=_gmsh_mean_size3(gsizes.vbgm[a],gsizes.vbgm[b],gsizes.vbgm[c])
    lc=extend_boundary ? min(background,
        _gmsh_mean_size3(gsizes.vsize[a],gsizes.vsize[b],gsizes.vsize[c])) : background
    first=_gmsh_flat_xyz(_pt(T,a))
    center=_gmsh_circumcenter_xyz(first,_gmsh_flat_xyz(_pt(T,b)),
                                _gmsh_flat_xyz(_pt(T,c)))
    return _gmsh_norm3(_gmsh_sub3(first,center))/lc
end

# Cartesian circumradius over the scalar size (`MTri3`, `radiusNorm == 2`).
@inline function _gmsh_cartesian_radius(wa,wb,wc,lc)
    center=_gmsh_circumcenter_xyz(wa,wb,wc)
    radius=_gmsh_norm3(_gmsh_sub3(wa,center))/lc
    return isfinite(radius) ? radius : 99999.0*sqrt(2.0)/lc
end

# Wrapped-edge test of the Cartesian mode. The parameter segment `pa`–`pb`
# lifts to a surface path; its chord `wa`–`wb` represents that path only while
# the lifted parameter midpoint stays within half a chord of the chord's own
# midpoint — a geodesic arc reaches that bound at a half turn. Only edges
# spanning more than a quarter of the chart are examined, so a resolved mesh
# pays nothing. At a collapsed endpoint (pole/apex: one parameter direction
# has no world length there) that parameter is free and follows the other
# endpoint, i.e. the path is the meridian actually drawn on the surface.
# Endpoint parameters of the surface path an edge draws. At a collapsed
# endpoint (pole/apex) the parameter direction without world length is free:
# it follows the other endpoint, so the path is the meridian, not the
# parameter-plane segment to whichever representative the cell happens to use.
@inline function _gmsh_edge_parameters(carrier,pa,pb)
    ua=pa[1];va=pa[2];ub=pb[1];vb=pb[2]
    ga=carrier.ground(pa[1],pa[2]);gb=carrier.ground(pb[1],pb[2])
    ga[1]<=1e-10*ga[3] && (ua=pb[1])
    ga[3]<=1e-10*ga[1] && (va=pb[2])
    gb[1]<=1e-10*gb[3] && (ub=pa[1])
    gb[3]<=1e-10*gb[1] && (vb=pa[2])
    return (ua,va),(ub,vb)
end

function _gmsh_edge_wrapped(carrier,pa,pb,wa,wb)
    bounds=carrier.bounds
    (abs(pa[1]-pb[1])>0.25*(bounds[2]-bounds[1]) ||
     abs(pa[2]-pb[2])>0.25*(bounds[4]-bounds[3])) || return false
    qa,qb=_gmsh_edge_parameters(carrier,pa,pb)
    middle=carrier.lift(0.5*(qa[1]+qb[1]),0.5*(qa[2]+qb[2]))
    dx=middle[1]-0.5*(wa[1]+wb[1]);dy=middle[2]-0.5*(wa[2]+wb[2])
    dz=middle[3]-0.5*(wa[3]+wb[3])
    chord=_gmsh_sub3(wb,wa)
    return dx*dx+dy*dy+dz*dz>0.25*_gmsh_dot3(chord,chord)
end

@inline function _gmsh_cell_wrapped(carrier,pa,pb,pc,wa,wb,wc)
    (isfinite(carrier.bounds[2]-carrier.bounds[1]) ||
     isfinite(carrier.bounds[4]-carrier.bounds[3])) || return false
    return _gmsh_edge_wrapped(carrier,pa,pb,wa,wb) ||
           _gmsh_edge_wrapped(carrier,pb,pc,wb,wc) ||
           _gmsh_edge_wrapped(carrier,pc,pa,wc,wa)
end

# End of the first insertion pass on a Cartesian carrier. If interior cells
# that wrap around the surface are left — a size larger than the model leaves
# the far side of a closed face covered by a few of them — the wrapped rule is
# armed, those cells are released from the skip list and `true` asks the
# caller to rebuild its queue and continue. Runs once per refinement; any
# other trace, and a carrier without a bounded chart, has no second pass.
_gmsh_surface_second_pass!(trace,T,interior,skipped)=false
function _gmsh_surface_second_pass!(carrier::_GmshSurfaceMetricCarrier,T,
                                    interior,skipped)
    (carrier.cartesian_close && !carrier.wrap[]) || return false
    (isfinite(carrier.bounds[2]-carrier.bounds[1]) ||
     isfinite(carrier.bounds[4]-carrier.bounds[3])) || return false
    found=false
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && t<=length(interior) && interior[t] &&
         !_is_ghost_tri(T,t)) || continue
        a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3)
        wa=_gmsh_surface_world(carrier,T,a)
        wb=_gmsh_surface_world(carrier,T,b)
        wc=_gmsh_surface_world(carrier,T,c)
        # Collapsed-row cells weld away; seam-degenerate cells are queued by
        # their own sentinel in every pass.
        (wa==wb || wb==wc || wc==wa) && continue
        _gmsh_cell_wrapped(carrier,_pt(T,a),_pt(T,b),_pt(T,c),wa,wb,wc) ||
            continue
        found=true
        delete!(skipped,_gmsh_triangle_key(T,t))
    end
    found && (carrier.wrap[]=true)
    return found
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
    u,v=_edge(T,t,slot)
    if gsizes.trace isa _GmshSurfaceMetricCarrier && gsizes.trace.cartesian_close
        return _gmsh_frontal_point_cartesian(T,t,slot,gsizes,extend_boundary,
                                             a,b,c,u,v)
    end
    # The carrier's `q` is metric-normalized: `extend_boundary` composes
    # the edge's local scalar bound as the same eigenvalue floor used by
    # the queue radius.
    extend_boundary && gsizes.trace isa _GmshSurfaceMetricCarrier &&
        (metric=_gmsh_surface_metric_floor(gsizes.trace,metric,
            _pt(T,a),_pt(T,b),_pt(T,c),
            0.5*(gsizes.vsize[u]+gsizes.vsize[v])))
    center,_=_gmsh_circumcenter_metric(_pt(T,a),_pt(T,b),_pt(T,c),
                                     metric)
    pu=_pt(T,u);pv=_pt(T,v)
    midpoint=(0.5*(pu[1]+pv[1]),0.5*(pu[2]+pv[2]))
    dx=center[1]-midpoint[1];dy=center[2]-midpoint[2]
    direction_norm=sqrt(dx*dx+dy*dy)
    q=sqrt(dx*metric[1]*dx+2.0*dx*metric[2]*dy+dy*metric[3]*dy)
    (isfinite(q) && q>0) || return nothing
    dirx=dx/direction_norm;diry=dy/direction_norm
    ratio=sqrt(dirx*dirx*metric[1]+2.0*diry*dirx*metric[2]+diry*diry*metric[3])
    distance=if gsizes.trace isa _GmshSurfaceMetricCarrier
        # In metric units a perfectly-sized triangle puts its circumcenter
        # `sqrt(3)/2` normalized units off the edge — the lc-free cap.
        min(sqrt(3.0)*0.5,q)
    else
        background=0.5*(gsizes.vbgm[u]+gsizes.vbgm[v])
        lc=extend_boundary ? min(background,
            0.5*(gsizes.vsize[u]+gsizes.vsize[v])) : background
        min(lc*sqrt(3.0)*0.5,q)
    end
    fast=(midpoint[1]+distance*dirx/ratio,
          midpoint[2]+distance*diry/ratio)
    opposite=_pt(T,_vert(T,t,slot))
    if gsizes.trace isa _GmshSurfaceMetricCarrier
        # optimalPointFrontalB's real surface circle: the circle is built in
        # world space from the lifted vertices, not in the parametric plane.
        wa=_gmsh_surface_world(gsizes.trace,T,u)
        wb=_gmsh_surface_world(gsizes.trace,T,v)
        wc=_gmsh_surface_world(gsizes.trace,T,_vert(T,t,slot))
        return _gmsh_surface_circle_intersection(fast,gsizes.trace.lift,
            wa,wb,wc)
    end
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

# `optimalPointFrontal`/`optimalPointFrontalB` for an isotropic background
# on a warped chart. The step along the Voronoi direction is the world length
# `L = min(ρ·√3/2, q)`: `ρ` the edge's scalar size, `q` the first-fundamental
# length from the edge midpoint to the metric circumcentre. The point is then
# the intersection of the surface with the world circle of radius `L` about
# the edge's world midpoint, Newton-solved from the parametric guess.
function _gmsh_frontal_point_cartesian(T,t,slot,gsizes,extend_boundary,a,b,c,u,v)
    carrier=gsizes.trace
    pa=_pt(T,a);pb=_pt(T,b);pc=_pt(T,c)
    pu=_pt(T,u);pv=_pt(T,v);po=_pt(T,_vert(T,t,slot))
    midpoint=(0.5*(pu[1]+pv[1]),0.5*(pu[2]+pv[2]))
    centroid=((pa[1]+pb[1]+pc[1])/3.0,(pa[2]+pb[2]+pc[2])/3.0)
    background=0.5*(gsizes.vbgm[u]+gsizes.vbgm[v])
    lc=extend_boundary ? min(background,
        0.5*(gsizes.vsize[u]+gsizes.vsize[v])) : background
    wa=_gmsh_surface_world(carrier,T,u)
    wb=_gmsh_surface_world(carrier,T,v)
    wc=_gmsh_surface_world(carrier,T,_vert(T,t,slot))
    # The wrapped construction below belongs to the second pass and to a
    # cell queued by the wrapped rule alone: one the Cartesian radius calls
    # resolved. Every other cell — in particular the open sector of a pole
    # fan, whose flat triangle is obtuse and therefore queued by its radius —
    # takes the ordinary construction, which builds the fan from equilateral
    # world triangles. A cell joining two seam representatives of one world
    # point does too (its frontal circle is oriented by the guess).
    second_pass=carrier.wrap[]
    wrapped=false
    if second_pass && !(wa==wb || wb==wc || wc==wa)
        cell_background=_gmsh_mean_size3(gsizes.vbgm[a],gsizes.vbgm[b],
                                         gsizes.vbgm[c])
        cell_lc=extend_boundary ? min(cell_background,_gmsh_mean_size3(
            gsizes.vsize[a],gsizes.vsize[b],gsizes.vsize[c])) : cell_background
        wrapped=_gmsh_cartesian_radius(wa,wb,wc,cell_lc)<
                    _GMSH_SURFACE_RADIUS_LIMIT &&
                _gmsh_cell_wrapped(carrier,pu,pv,po,wa,wb,wc)
    end
    if !wrapped
        ground=carrier.ground(centroid[1],centroid[2])
        center,_=_gmsh_circumcenter_metric(pa,pb,pc,ground)
        dx=center[1]-midpoint[1];dy=center[2]-midpoint[2]
        direction_norm=sqrt(dx*dx+dy*dy)
        q=sqrt(dx*ground[1]*dx+2.0*dx*ground[2]*dy+dy*ground[3]*dy)
        (isfinite(q) && q>0 && direction_norm>0) || return nothing
        dirx=dx/direction_norm;diry=dy/direction_norm
        ratio=sqrt(dirx*dirx*ground[1]+2.0*diry*dirx*ground[2]+
                   diry*diry*ground[3])
        distance=min(lc*sqrt(3.0)*0.5,q)
        fast=(midpoint[1]+distance*dirx/ratio,midpoint[2]+distance*diry/ratio)
        point=_gmsh_surface_circle_intersection(fast,carrier.lift,wa,wb,wc;
                                                radius=distance)
        # The solve can settle on a parameter alias outside the chart. That
        # candidate is simply rejected downstream; in the second pass, on a
        # cell that also wraps, the rejection would leave the far side
        # unfilled, so that one continues with the wrapped construction.
        (!second_pass || _gmsh_carrier_inside(carrier,point) ||
         wa==wb || wb==wc || wc==wa ||
         !_gmsh_cell_wrapped(carrier,pu,pv,po,wa,wb,wc)) &&
            return (point[1],point[2])
    end
    # A cell that wraps around the surface has no usable centroid metric or
    # circumcentre. The front still advances from its active edge: step the
    # world length `ρ·√3/2` from the edge midpoint along the parameter
    # direction that the first fundamental form at that midpoint makes
    # orthogonal to the edge, toward the opposite vertex, and solve the same
    # world circle from there. The edge is taken along the path it draws on
    # the surface (`_gmsh_edge_parameters`), so an edge into a pole is its
    # meridian and the step leaves it sideways.
    qu,qv=_gmsh_edge_parameters(carrier,pu,pv)
    middle=(0.5*(qu[1]+qv[1]),0.5*(qu[2]+qv[2]))
    ground=carrier.ground(middle[1],middle[2])
    ex=qv[1]-qu[1];ey=qv[2]-qu[2]
    dx=-(ground[2]*ex+ground[3]*ey);dy=ground[1]*ex+ground[2]*ey
    dx*(po[1]-middle[1])+dy*(po[2]-middle[2])<0 && (dx=-dx;dy=-dy)
    extent=sqrt(dx*ground[1]*dx+2.0*dx*ground[2]*dy+dy*ground[3]*dy)
    distance=lc*sqrt(3.0)*0.5
    if isfinite(extent) && extent>0
        guess=(middle[1]+distance*dx/extent,middle[2]+distance*dy/extent)
        if _gmsh_carrier_inside(carrier,guess)
            point=_gmsh_surface_circle_intersection(guess,carrier.lift,
                wa,wb,wc;radius=distance)
            # Only a point that is on the frontal circle is a frontal point:
            # at the world distance `L` from the edge midpoint and as far
            # from one endpoint as from the other. A failed solve hands back
            # its guess, which is neither.
            if _gmsh_carrier_inside(carrier,point)
                w=carrier.lift(point[1],point[2])
                wm=(0.5*(wa[1]+wb[1]),0.5*(wa[2]+wb[2]),0.5*(wa[3]+wb[3]))
                da=_gmsh_norm3(_gmsh_sub3(w,wa));db=_gmsh_norm3(_gmsh_sub3(w,wb))
                abs(_gmsh_norm3(_gmsh_sub3(w,wm))-distance)<=0.05*distance &&
                    abs(da-db)<=0.05*distance && return (point[1],point[2])
            end
        end
    end
    # No frontal point. Splitting the cell at its parameter centroid still
    # makes progress, but only where that lands clear of the cell's own
    # vertices: beside a pole it would plant a vertex half a size from it.
    w=carrier.lift(centroid[1],centroid[2])
    min(_gmsh_norm3(_gmsh_sub3(w,wa)),_gmsh_norm3(_gmsh_sub3(w,wb)),
        _gmsh_norm3(_gmsh_sub3(w,wc)))>=0.7*lc && return centroid
    return nothing
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
