# The public native producer retains the literal MTriangle slots, Cartesian
# radius arithmetic and every MElement construction, including rejected cells.
# Generic/raw refinement keeps its existing storage and arithmetic contract.
mutable struct _GmshSurfaceTrace{F,G}
    world::Vector{NTuple{3,Float64}}
    lift::F
    metric::G
    # Per triangulation slot: the literal MTriangle node order, the element
    # label and the queue radius of the cell stored there. (These were three
    # dictionaries keyed by the sorted node triple, hashed on every read.)
    triangles::Vector{NTuple{3,Int32}}
    elements::Vector{UInt64}
    radii::Vector{Float64}
    element_high_water::UInt64
    extend_boundary::Bool
    # The cells of the candidate under evaluation, one per shell edge in
    # shell order. They are constructed (and numbered) before the insertion
    # is accepted; `_gmsh_surface_trace_bind!` moves them to their slots once
    # the kernel has created the cells.
    child_triangles::Vector{NTuple{3,Int32}}
    child_elements::Vector{UInt64}
    child_radii::Vector{Float64}
end

@inline _gmsh_surface_triangle_nodes(T,t,::Nothing)=
    (_vert(T,t,1),_vert(T,t,2),_vert(T,t,3))
@inline _gmsh_surface_triangle_nodes(T,t,trace::_GmshSurfaceTrace)=
    trace.triangles[t]
@inline _gmsh_surface_incoming_slots(T,t,::Nothing)=(2,3,1)
@inline function _gmsh_surface_incoming_slots(T,t,trace::_GmshSurfaceTrace)
    logical=_gmsh_surface_triangle_nodes(T,t,trace)
    storage=(_vert(T,t,1),_vert(T,t,2),_vert(T,t,3))
    return ntuple(i->something(findfirst(==(logical[(2,3,1)[i]]),storage)),3)
end
@inline function _gmsh_surface_world(trace,T,row)
    row<=length(trace.world) && return trace.world[row]
    return trace.lift(T.x[row],T.y[row])
end
function _gmsh_surface_trace_next_element!(trace)
    trace.element_high_water<typemax(UInt64) || throw(ArgumentError(
        "native Surface element construction exhausts UInt64 public labels"))
    trace.element_high_water+=UInt64(1)
    return trace.element_high_water
end
function _gmsh_surface_trace_radius(trace,T,nodes,lc)
    a,b,c=nodes
    first=_gmsh_surface_world(trace,T,a)
    center=_gmsh_circumcenter_xyz(first,_gmsh_surface_world(trace,T,b),
                                 _gmsh_surface_world(trace,T,c))
    return _gmsh_norm3(_gmsh_sub3(first,center))/lc
end
function _gmsh_surface_trace(T,inside,world,lift,high_water,sizes,background,extend;
                             temporary_cells=0,metric=(1.0,0.0,1.0),initial=nothing)
    0<=high_water<=typemax(UInt64) || throw(ArgumentError(
        "native Surface element high water must be UInt64 representable"))
    length(world)==T.nreal || throw(ArgumentError(
        "native Surface Cartesian birth rows must align with initial vertices"))
    slots=length(T.alive)
    trace=_GmshSurfaceTrace(copy(world),lift,metric,
        fill((Int32(0),Int32(0),Int32(0)),slots),zeros(UInt64,slots),zeros(slots),
        UInt64(high_water),extend,NTuple{3,Int32}[],UInt64[],Float64[])
    temporary_cells isa Integer && !(temporary_cells isa Bool) && temporary_cells>=0 ||
        throw(ArgumentError("native temporary Surface cell count must be a nonnegative integer"))
    # meshGenerator constructs the recovered initial faces, delaunayizes BDS,
    # deletes the temporary GFace mesh, then BDS2GMSH constructs the queue input.
    for _ in 1:temporary_cells
        _gmsh_surface_trace_next_element!(trace)
    end
    # `initial` lists the cells as BDS2GMSH constructs them after
    # `delaunayizeBDS` — (slot, MTriangle node order) — when that pass may have
    # swapped edges; otherwise the slots are the recovered PolyMesh faces.
    cells=initial===nothing ?
        ((triangle,(_vert(T,triangle,1),_vert(T,triangle,3),_vert(T,triangle,2)))
         for triangle in eachindex(T.alive) if inside[triangle]) : initial
    for (triangle,nodes) in cells
        trace.triangles[triangle]=nodes
        trace.elements[triangle]=_gmsh_surface_trace_next_element!(trace)
        a,b,c=nodes
        # buildMeshGenerationDataStructures uses this decimal literal for the
        # initial queue; inserted MTri3 cells use ONE_THIRD instead.
        lc=0.3333333333*(sizes[a]+sizes[b]+sizes[c])
        bg=0.3333333333*(background[a]+background[b]+background[c])
        trace.radii[triangle]=_gmsh_surface_trace_radius(trace,T,nodes,extend ? min(lc,bg) : bg)
    end
    return trace
end
@inline _gmsh_surface_metric(::Nothing)=(1.0,0.0,1.0)
@inline _gmsh_surface_metric(trace::_GmshSurfaceTrace)=trace.metric
# Per-triangle tensor arm: a constant (a,b,d) pullback stays constant, while a
# callable metric source evaluates at the cell's parameter centroid so the
# cavity in-circle test and circumcenter use the local background tensor.
@inline _gmsh_surface_metric(::Nothing,pu,pv,pw)=(1.0,0.0,1.0)
@inline _gmsh_surface_metric(trace::_GmshSurfaceTrace,pu,pv,pw)=
    _gmsh_trace_metric(trace.metric,pu,pv,pw)
@inline _gmsh_trace_metric(metric::NTuple{3,Float64},pu,pv,pw)=metric
@inline _gmsh_trace_metric(metric,pu,pv,pw)=metric(
    (pu[1]+pv[1]+pw[1])/3.0,(pu[2]+pv[2]+pw[2])/3.0)

# `buildMetric(gf, uv)` of a curved chart: the first fundamental form
# `(E,F,G)` of `form(u,v)`. Its cells read it at the centroid exactly as
# `optimalPointFrontal` and `inCircumCircleAniso` form it, `* ONE_THIRD`; the
# frontal point then lies on the surface (`optimalPointFrontalB`) and every
# test of one insertion uses the seed cell's tensor (`metricb`).
struct _GmshFirstFundamental{F}
    form::F
end
(metric::_GmshFirstFundamental)(u,v)=metric.form(u,v)
const _GMSH_ONE_THIRD=1.0/3.0
@inline _gmsh_trace_metric(metric::_GmshFirstFundamental,pu,pv,pw)=metric.form(
    (pu[1]+pv[1]+pw[1])*_GMSH_ONE_THIRD,(pu[2]+pv[2]+pw[2])*_GMSH_ONE_THIRD)
@inline _gmsh_trace_curved(trace)=false
@inline _gmsh_trace_curved(trace::_GmshSurfaceTrace)=trace.metric isa _GmshFirstFundamental
# The tensor every cavity test of one insertion uses on a curved chart, or
# `nothing` where each cell reads its own (a constant tensor makes both equal).
@inline _gmsh_seed_metric(trace,T,a,b,c)=_gmsh_trace_curved(trace) ?
    _gmsh_surface_metric(trace,_pt(T,a),_pt(T,b),_pt(T,c)) : nothing
@inline _gmsh_surface_delaunay_point(T,a,b,c,::Nothing)=
    _gmsh_circumcenter_xy(_pt(T,a),_pt(T,b),_pt(T,c))
@inline function _gmsh_surface_delaunay_point(T,a,b,c,trace::_GmshSurfaceTrace)
    # bowyerWatson recomputes circumCenterMetric after circUV. The earlier
    # Euclidean XY construction is overwritten and must not decide a birth.
    pa=_pt(T,a);pb=_pt(T,b);pc=_pt(T,c)
    metric=_gmsh_surface_metric(trace,pa,pb,pc)
    center,_=_gmsh_circumcenter_metric(pa,pb,pc,metric)
    return center
end
# The three node rows in ascending order: the key of `_gmsh_triangle_key`
# for a cell that is not stored yet, without sorting a temporary vector.
@inline function _gmsh_sorted_key(a::Int32,b::Int32,c::Int32)
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end
@inline _gmsh_sorted_key(a,b,c)=_gmsh_sorted_key(Int32(a),Int32(b),Int32(c))
function _gmsh_surface_trace_children!(trace,T,cavity,boundary,vid,gsizes,lc1,bg)
    trace===nothing && return nothing
    empty!(trace.child_triangles);empty!(trace.child_elements)
    empty!(trace.child_radii)
    length(cavity)>1 && length(boundary)==length(cavity)+2 || return nothing
    for (a,b,_) in boundary
        nodes=(b,a,vid)
        push!(trace.child_triangles,nodes)
        push!(trace.child_elements,_gmsh_surface_trace_next_element!(trace))
        lc=(gsizes.vsize[b]+gsizes.vsize[a]+lc1)*(1.0/3.0)
        background=(gsizes.vbgm[b]+gsizes.vbgm[a]+bg)*(1.0/3.0)
        push!(trace.child_radii,_gmsh_surface_trace_radius(trace,T,nodes,
            trace.extend_boundary ? min(lc,background) : background))
    end
    return nothing
end

# The kernel created one cell per shell edge of the accepted candidate, in
# shell order: `slots[i]` holds the i-th child registered above.
_gmsh_surface_trace_bind!(trace,slots)=nothing
function _gmsh_surface_trace_bind!(trace::_GmshSurfaceTrace,slots::Vector{Int32})
    length(slots)==length(trace.child_triangles) || throw(ErrorException(
        "native Surface trace: new cells do not match the accepted shell"))
    highest=maximum(slots;init=Int32(0))
    if highest>length(trace.triangles)
        resize!(trace.triangles,highest);resize!(trace.elements,highest)
        resize!(trace.radii,highest)
    end
    @inbounds for (index,slot) in enumerate(slots)
        trace.triangles[slot]=trace.child_triangles[index]
        trace.elements[slot]=trace.child_elements[index]
        trace.radii[slot]=trace.child_radii[index]
    end
    return nothing
end

# `GPoint::succeeded` of `GFace::point(u,v)`; a chart that can fail to
# evaluate (a discreteFace outside its parameter triangles) overrides it.
_gmsh_surface_point_succeeded(lift,u,v)=true

# insertVertexB evaluates the retained MTriangle slots in UV, but all close
# point distances use the actual Cartesian MVertex coordinates. The raw
# insertion path retains its separate existing admission contract.
function _gmsh_native_trace_cavity_accept(T,cavity,boundary,parameter,onmodel,
                                         length_scale,trace,verify_star)
    # insertAPoint calls insertVertexB only when `gf->point` succeeded.
    _gmsh_surface_point_succeeded(trace.lift,parameter[1],parameter[2]) || return false
    length(cavity)==1 && return false
    length(boundary)==length(cavity)+2 || return false
    old_volume=0.0
    for triangle in cavity
        _is_ghost_tri(T,triangle) && return false
        a,b,c=_gmsh_surface_triangle_nodes(T,triangle,trace)
        old_volume+=_gmsh_surface_area2(_pt(T,a),_pt(T,b),_pt(T,c))
    end
    point=trace.lift(parameter...)
    new_volume=0.0;too_close=false
    for (a,b,_) in boundary
        first=_gmsh_surface_world(trace,T,b)
        second=_gmsh_surface_world(trace,T,a)
        d1=_gmsh_norm3(_gmsh_sub3(first,point))
        d2=_gmsh_norm3(_gmsh_sub3(second,point))
        d3=_gmsh_norm3(_gmsh_sub3(first,second))
        limit=length_scale(b,a)
        cosine=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2)
        perpendicular=1e22
        if b in onmodel && a in onmodel
            cross=_gmsh_cross3(_gmsh_sub3(second,first),_gmsh_sub3(point,first))
            perpendicular=_gmsh_norm3(cross)/d3
        end
        if trace isa _GmshSurfaceMetricCarrier && !trace.cartesian_close
            # Metric-only carrier: the local `vsize` bound folds into the
            # tensor as the same isotropic size floor the queue radius uses
            # (`min(lc_bgm, lc_vsize)` in metric space — the tensor already
            # carries `lc_bgm`); distances are read in normalized units so
            # the `.5`/`.4` bounds apply to the directional local length
            # scale, not the min-direction scalar.
            pb=(T.x[b],T.y[b]);pa=(T.x[a],T.y[a])
            vs=length_scale isa _GmshShellLL ?
                _gmsh_mean_size3(length_scale.sizes.vsize[b],
                                length_scale.sizes.vsize[a],
                                length_scale.propagated) : limit
            mn=_gmsh_surface_metric_floor(trace,
                _gmsh_surface_metric(trace,pb,pa,parameter),pb,pa,parameter,vs)
            d1n=sqrt(max(_gmsh_metric_distance2(pb,parameter,mn),0.0))
            d2n=sqrt(max(_gmsh_metric_distance2(pa,parameter,mn),0.0))
            close=d1n<.5 || d2n<.5
            if b in onmodel && a in onmodel
                ex=pa[1]-pb[1];ey=pa[2]-pb[2];e2=ex*ex+ey*ey
                if e2>0
                    proj=((parameter[1]-pb[1])*ex+(parameter[2]-pb[2])*ey)/e2
                    qx=parameter[1]-pb[1]-proj*ex
                    qy=parameter[2]-pb[2]-proj*ey
                    close |= sqrt(max(qx*qx*mn[1]+qy*qy*mn[3]+
                                      2.0*qx*qy*mn[2],0.0))<.4
                end
            end
            (close || cosine<-.9999) && (too_close=true)
        else
            (d1<limit*.5 || d2<limit*.5 || perpendicular<limit*.4 ||
             cosine<-.9999) && (too_close=true)
        end
        area=_gmsh_surface_area2(_pt(T,b),_pt(T,a),parameter)
        area<1e-25 && (area=1e22)
        new_volume+=area
    end
    epsilon=verify_star ? 1e-12 : 1e12
    return abs(old_volume-new_volume)<epsilon*old_volume && !too_close
end
# The `metricb` overload of inCircumCircleAniso evaluates every cavity
# member against the candidate's metric — the tensor `buildMetric` returns
# for the seed triangle's centroid in `insertAPoint`'s caller — not a fresh
# per-member tensor. `metric[]` is armed per candidate by `refine!`; a
# `nothing` payload keeps the per-member evaluation for other callers.
mutable struct _GmshSurfaceCircleTest{S}
    trace::S
    metric::Base.RefValue{Any}
end
_GmshSurfaceCircleTest(trace)=_GmshSurfaceCircleTest{typeof(trace)}(trace,Ref{Any}(nothing))
function (test::_GmshSurfaceCircleTest)(T,t,vid)
    test.metric[]===nothing && return _gmsh_in_circumcircle(T,t,vid,test.trace)
    _is_ghost_tri(T,t) && return false
    first,second,third=_gmsh_surface_triangle_nodes(T,t,test.trace)
    a=_pt(T,first);b=_pt(T,second);c=_pt(T,third)
    center,radius=_gmsh_circumcenter_metric(a,b,c,test.metric[])
    point=_pt(T,vid)
    return _gmsh_metric_distance2(center,point,test.metric[])<radius
end

# Metric-only carrier for the ordinary (non-trace) Surface algorithms: a
# tensor background rides the same dispatch points as a full trace —
# cavity in-circle, circumcenter, frontal direction, the Cartesian
# admission distances, queue radii and the close-point `LL` bound —
# while triangle slot order, incoming-edge order and element bookkeeping
# keep the untraced scalar semantics. Scalar local bounds (`vsize`, `LL`)
# fold into the tensor through `_gmsh_metric_floor` — the metric-space
# analogue of `min(lc_bgm, lc_vsize)`.
# `world` covers the birth rows exactly like `_gmsh_surface_trace`, so the
# admission path's `trace.world`/`trace.lift` reads stay consistent as
# inserted vertices append world rows.
struct _GmshSurfaceMetricCarrier{L,G,G2}
    world::Vector{NTuple{3,Float64}}
    lift::L
    metric::G
    # First-fundamental tensor `G(u,v)=JᵀJ` of the lifted domain — `nothing`
    # for isometric (planar) coordinates where UV units are world units.
    # Armed on genuinely warped parameter domains so a world scalar bound
    # `vs` folds into the UV metric as `G/vs²` (the `I/vs²` pullback), not
    # the dimensionless `I/vs²` eigenvalue floor `_gmsh_metric_floor` uses.
    ground::G2
    # `insertVertexB` measures its close-point bounds on the actual Cartesian
    # vertices. An isotropic background on a warped chart keeps that world
    # test (`true`): a metric read at the cell centroid overestimates the
    # distance to a collapsed pole/apex row, whose parameter representatives
    # lie far apart while their world points coincide. Tensor backgrounds
    # keep the metric-normalized bounds (`false`).
    cartesian_close::Bool
    # Parameter box `(umin, umax, vmin, vmax)` of the chart for the Cartesian
    # mode: the wrapped-cell test (`_gmsh_cell_wrapped`) and the admission of
    # a Newton-solved frontal point. An unbounded box disables both.
    bounds::NTuple{4,Float64}
    # Second-pass switch of the Cartesian mode. The insertion first runs on
    # the Cartesian radius alone; only if cells that wrap around the surface
    # are left when its queue empties (`_gmsh_surface_second_pass!`) is the
    # wrapped rule armed and the queue rebuilt. A resolved mesh never sees it.
    wrap::Base.RefValue{Bool}
end
const _GMSH_CARRIER_UNBOUNDED=(-Inf,Inf,-Inf,Inf)
@inline _GmshSurfaceMetricCarrier(world,lift,metric)=
    _GmshSurfaceMetricCarrier(world,lift,metric,nothing,false,
                              _GMSH_CARRIER_UNBOUNDED,Ref(false))
@inline _GmshSurfaceMetricCarrier(world,lift,metric,ground)=
    _GmshSurfaceMetricCarrier(world,lift,metric,ground,false,
                              _GMSH_CARRIER_UNBOUNDED,Ref(false))
@inline _GmshSurfaceMetricCarrier(world,lift,metric,ground,cartesian_close::Bool)=
    _GmshSurfaceMetricCarrier(world,lift,metric,ground,cartesian_close,
                              _GMSH_CARRIER_UNBOUNDED,Ref(false))
@inline _GmshSurfaceMetricCarrier(world,lift,metric,ground,cartesian_close::Bool,
                                  bounds::NTuple{4,Float64})=
    _GmshSurfaceMetricCarrier(world,lift,metric,ground,cartesian_close,
                              bounds,Ref(false))
@inline _gmsh_carrier_inside(carrier::_GmshSurfaceMetricCarrier,point)=
    isfinite(point[1]) && isfinite(point[2]) &&
    carrier.bounds[1]<=point[1]<=carrier.bounds[2] &&
    carrier.bounds[3]<=point[2]<=carrier.bounds[4]
@inline _gmsh_surface_triangle_nodes(T,t,::_GmshSurfaceMetricCarrier)=
    (_vert(T,t,1),_vert(T,t,2),_vert(T,t,3))
@inline _gmsh_surface_incoming_slots(T,t,::_GmshSurfaceMetricCarrier)=(2,3,1)
@inline _gmsh_surface_metric(carrier::_GmshSurfaceMetricCarrier)=carrier.metric
@inline _gmsh_surface_metric(carrier::_GmshSurfaceMetricCarrier,pu,pv,pw)=
    _gmsh_trace_metric(carrier.metric,pu,pv,pw)
@inline function _gmsh_surface_world(carrier::_GmshSurfaceMetricCarrier,T,row)
    row<=length(carrier.world) && return carrier.world[row]
    return carrier.lift(T.x[row],T.y[row])
end
@inline function _gmsh_surface_delaunay_point(T,a,b,c,carrier::_GmshSurfaceMetricCarrier)
    pa=_pt(T,a);pb=_pt(T,b);pc=_pt(T,c)
    metric=_gmsh_surface_metric(carrier,pa,pb,pc)
    center,_=_gmsh_circumcenter_metric(pa,pb,pc,metric)
    return center
end
# Scalar local bound folded into the UV metric at the evaluation point
# `(pu,pv,pw)` — the trace centroid like `_gmsh_surface_metric`. World-aligned
# traces keep `_gmsh_metric_floor`; a carrier armed with `ground` folds the
# pulled-back bound `G/h²` through `_gmsh_metric_ground_floor`.
@inline _gmsh_surface_metric_floor(trace,metric,pu,pv,pw,h)=
    _gmsh_metric_floor(metric,h)
@inline function _gmsh_surface_metric_floor(
        carrier::_GmshSurfaceMetricCarrier,metric,pu,pv,pw,h)
    carrier.ground===nothing && return _gmsh_metric_floor(metric,h)
    u=(pu[1]+pv[1]+pw[1])/3.0;v=(pu[2]+pv[2]+pw[2])/3.0
    return _gmsh_metric_ground_floor(metric,carrier.ground(u,v),h)
end
function _gmsh_surface_trace_children!(::_GmshSurfaceMetricCarrier,T,cavity,
                                       boundary,vid,gsizes,lc1,bg)
    return nothing
end
