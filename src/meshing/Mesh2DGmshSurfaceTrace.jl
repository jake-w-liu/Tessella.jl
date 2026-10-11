# The public native producer retains the literal MTriangle slots, Cartesian
# radius arithmetic and every MElement construction, including rejected cells.
# Generic/raw refinement keeps its existing storage and arithmetic contract.
mutable struct _GmshSurfaceTrace{F}
    world::Vector{NTuple{3,Float64}}
    lift::F
    metric::NTuple{3,Float64}
    triangles::Dict{NTuple{3,Int32},NTuple{3,Int32}}
    elements::Dict{NTuple{3,Int32},UInt64}
    radii::Dict{NTuple{3,Int32},Float64}
    element_high_water::UInt64
    extend_boundary::Bool
end

@inline _gmsh_surface_triangle_nodes(T,t,::Nothing)=
    (_vert(T,t,1),_vert(T,t,2),_vert(T,t,3))
@inline _gmsh_surface_triangle_nodes(T,t,trace::_GmshSurfaceTrace)=
    trace.triangles[_gmsh_triangle_key(T,t)]
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
                             temporary_cells=0,metric=(1.0,0.0,1.0))
    0<=high_water<=typemax(UInt64) || throw(ArgumentError(
        "native Surface element high water must be UInt64 representable"))
    length(world)==T.nreal || throw(ArgumentError(
        "native Surface Cartesian birth rows must align with initial vertices"))
    trace=_GmshSurfaceTrace(copy(world),lift,metric,Dict{NTuple{3,Int32},NTuple{3,Int32}}(),
        Dict{NTuple{3,Int32},UInt64}(),Dict{NTuple{3,Int32},Float64}(),UInt64(high_water),extend)
    temporary_cells isa Integer && !(temporary_cells isa Bool) && temporary_cells>=0 ||
        throw(ArgumentError("native temporary Surface cell count must be a nonnegative integer"))
    # meshGenerator constructs the recovered initial faces, delaunayizes BDS,
    # deletes the temporary GFace mesh, then BDS2GMSH constructs the queue input.
    for _ in 1:temporary_cells
        _gmsh_surface_trace_next_element!(trace)
    end
    for triangle in eachindex(T.alive)
        inside[triangle] || continue
        nodes=(_vert(T,triangle,1),_vert(T,triangle,3),_vert(T,triangle,2))
        key=_gmsh_triangle_key(T,triangle);trace.triangles[key]=nodes
        trace.elements[key]=_gmsh_surface_trace_next_element!(trace)
        a,b,c=nodes
        # buildMeshGenerationDataStructures uses this decimal literal for the
        # initial queue; inserted MTri3 cells use ONE_THIRD instead.
        lc=0.3333333333*(sizes[a]+sizes[b]+sizes[c])
        bg=0.3333333333*(background[a]+background[b]+background[c])
        trace.radii[key]=_gmsh_surface_trace_radius(trace,T,nodes,extend ? min(lc,bg) : bg)
    end
    return trace
end
@inline _gmsh_surface_metric(::Nothing)=(1.0,0.0,1.0)
@inline _gmsh_surface_metric(trace::_GmshSurfaceTrace)=trace.metric
@inline _gmsh_surface_delaunay_point(T,a,b,c,::Nothing)=
    _gmsh_circumcenter_xy(_pt(T,a),_pt(T,b),_pt(T,c))
@inline function _gmsh_surface_delaunay_point(T,a,b,c,trace::_GmshSurfaceTrace)
    # bowyerWatson recomputes circumCenterMetric after circUV. The earlier
    # Euclidean XY construction is overwritten and must not decide a birth.
    center,_=_gmsh_circumcenter_metric(_pt(T,a),_pt(T,b),_pt(T,c),trace.metric)
    return center
end
function _gmsh_surface_trace_children!(trace,T,cavity,boundary,vid,gsizes,lc1,bg)
    trace===nothing && return nothing
    length(cavity)>1 && length(boundary)==length(cavity)+2 || return nothing
    for (a,b,_) in boundary
        nodes=(b,a,vid);key=Tuple(sort!(collect(nodes)))
        trace.triangles[key]=nodes
        trace.elements[key]=_gmsh_surface_trace_next_element!(trace)
        lc=(gsizes.vsize[b]+gsizes.vsize[a]+lc1)*(1.0/3.0)
        background=(gsizes.vbgm[b]+gsizes.vbgm[a]+bg)*(1.0/3.0)
        trace.radii[key]=_gmsh_surface_trace_radius(trace,T,nodes,
            trace.extend_boundary ? min(lc,background) : background)
    end
    return nothing
end

# insertVertexB evaluates the retained MTriangle slots in UV, but all close
# point distances use the actual Cartesian MVertex coordinates. The raw
# insertion path retains its separate existing admission contract.
function _gmsh_native_trace_cavity_accept(T,cavity,boundary,parameter,onmodel,
                                         length_scale,trace,verify_star)
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
        (d1<limit*.5 || d2<limit*.5 || perpendicular<limit*.4 || cosine<-.9999) &&
            (too_close=true)
        area=_gmsh_surface_area2(_pt(T,b),_pt(T,a),parameter)
        area<1e-25 && (area=1e22)
        new_volume+=area
    end
    epsilon=verify_star ? 1e-12 : 1e12
    return abs(old_volume-new_volume)<epsilon*old_volume && !too_close
end
struct _GmshSurfaceCircleTest{S}
    trace::S
end
(test::_GmshSurfaceCircleTest)(T,t,vid)=_gmsh_in_circumcircle(T,t,vid,test.trace)
