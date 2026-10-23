# Physical plane coordinates for the native algorithms. The normal's dropped
# component is the best-conditioned projection selected by Model. Cholesky
# factors of the induced metric convert kept world axes to Euclidean lengths.
struct _SurfaceMetricFrame{P}
    plane::P
    identity::Bool
    alpha::Float64
    beta::Float64
    gamma::Float64
end

function _surface_metric_frame(plane)
    axes=plane.axes
    du=-plane.normal[axes[1]]/plane.normal[plane.k]
    dv=-plane.normal[axes[2]]/plane.normal[plane.k]
    alpha=hypot(1.0,du)
    beta=du*(dv/alpha)
    gamma=hypot(1.0,dv/alpha)
    return _SurfaceMetricFrame(plane,iszero(du) && iszero(dv),alpha,beta,gamma)
end

@inline function _surface_metric_point(frame::_SurfaceMetricFrame,u,v)
    frame.identity && return (u,v)
    plane=frame.plane
    x=u-plane.anchor[plane.axes[1]]
    y=v-plane.anchor[plane.axes[2]]
    return (frame.alpha*x+frame.beta*y,frame.gamma*y)
end

@inline function _surface_metric_unproject(frame::_SurfaceMetricFrame,x,y)
    frame.identity && return (x,y)
    plane=frame.plane
    v=y/frame.gamma
    u=(x-frame.beta*v)/frame.alpha
    return (u+plane.anchor[plane.axes[1]],v+plane.anchor[plane.axes[2]])
end

@inline function _surface_metric_lift(frame::_SurfaceMetricFrame,x,y)
    plane=frame.plane
    u,v=_surface_metric_unproject(frame,x,y)
    dropped=_plane_dropped_coordinate(plane,u,v)
    return ntuple(axis->axis==plane.k ? dropped :
        axis==plane.axes[1] ? u : v,3)
end

function _surface_metric_pslg(frame,xs,ys,canonical,caller,t)
    frame.identity && return (xs,ys,canonical)
    px=Vector{Float64}(undef,length(xs))
    py=Vector{Float64}(undef,length(ys))
    original=Dict{NTuple{2,Float64},NTuple{2,Float64}}()
    for node in eachindex(xs)
        x,y=_surface_metric_point(frame,xs[node],ys[node])
        (isfinite(x) && isfinite(y)) || throw(ArgumentError(
            "$caller: physical plane coordinates of Surface[$t] overflow"))
        key=_surface_size_key(x,y)
        source=_surface_size_key(xs[node],ys[node])
        haskey(original,key) && original[key]!=source && throw(ArgumentError(
            "$caller: distinct Surface[$t] vertices coincide at physical " *
            "plane Float64 resolution"))
        original[key]=source
        px[node]=x;py[node]=y
    end
    lifted=Dict{NTuple{2,Float64},NTuple{3,Float64}}()
    for (uv,point) in canonical
        x,y=_surface_metric_point(frame,uv...)
        key=_surface_size_key(x,y)
        haskey(lifted,key) && lifted[key]!=point && throw(ArgumentError(
            "$caller: distinct canonical vertices of Surface[$t] coincide " *
            "at physical plane Float64 resolution"))
        lifted[key]=point
    end
    return px,py,lifted
end

# This admission runs only for the ordinary surface mesher. TF, extrusion
# and periodic slave copies have already selected their own construction.
# `metric=true` admits upstream's algorithm 7 selector (BAMG), which the
# ordinary dispatcher resolves to the native conforming tensor-bisection
# construction in ModelMetricSurfaceProducer; the OCC-curved and primitive
# factory kernels keep the unextended set because 7 names machinery they do
# not have.
function _normal_surface_mesh_algorithm(m::GeoModel,t::Int,caller;metric::Bool=false)
    algorithm=_surface_mesh_algorithm(m,t)
    algorithm==2 && (algorithm=_surface_type(m,t)==:plane ? 5 : 1)
    algorithm in (metric ? (3,5,6,7) : (3,5,6)) || throw(ArgumentError(
        "$caller: mesh algorithm $algorithm on Surface[$t] is not " *
        "implemented by the ordinary native surface mesher; supported: " *
        "2 (Auto on Plane), 3 (initial only), 5 (Delaunay), " *
        "6 (Frontal-Delaunay)" *
        (metric ? ", 7 (metric-driven bisection)" : "")))
    return algorithm
end

# Tensor background arm for the ordinary Delaunay/Frontal refinement: when
# the resolved size selection is an `AbstractAnisoField` the insertion kernel
# receives a metric-only carrier, so every cavity in-circle test,
# circumcenter and frontal direction step uses the tensor pulled back at the
# triangle's lifted centroid. The directional metric gets the same global
# size policy as `_metric_surface_algorithm_field` — the CTX::lc cap, the
# [lcMin,lcMax] clamp and lcFactor applied per direction — while queue radii
# and size propagation keep the untraced scalar semantics (`size_at` on a
# tensor scalarizes to the strongest direction). A scalar or absent
# selection returns `nothing`, leaving the Euclidean arm unchanged.
function _surface_metric_trace_carrier(m,options,size_field,frame,T,t,caller)
    raw=options.inner_field===nothing ? size_field : options.inner_field
    raw isa AbstractAnisoField || return nothing
    lc_min=isnan(options.lc_min) ? 0.0 : options.lc_min
    lc_max=isnan(options.lc_max) ? GMSH_MAX_SIZE : options.lc_max
    lo=min(lc_min,lc_max)
    hi=max(lo,min(lc_max,options.ctx_lc))
    wrapped=BoundedSize(raw;size_min=lo,size_max=hi,factor=options.lc_factor)
    # The PSLG frame's lift is affine and isometric, so the world-metric
    # pullback onto its coordinate columns e1,e2 gives the exact uv tensor.
    origin=_surface_metric_lift(frame,0.0,0.0)
    e1=ntuple(k->_surface_metric_lift(frame,1.0,0.0)[k]-origin[k],3)
    e2=ntuple(k->_surface_metric_lift(frame,0.0,1.0)[k]-origin[k],3)
    entity=(2,t)
    metric=function(u,v)
        w=_surface_metric_lift(frame,u,v)
        mm=metric_at(wrapped,w[1],w[2],w[3],entity)
        return (mm.m11*e1[1]*e1[1]+mm.m22*e1[2]*e1[2]+mm.m33*e1[3]*e1[3]+
            2*(mm.m12*e1[1]*e1[2]+mm.m13*e1[1]*e1[3]+mm.m23*e1[2]*e1[3]),
            mm.m11*e1[1]*e2[1]+mm.m22*e1[2]*e2[2]+mm.m33*e1[3]*e2[3]+
            mm.m12*(e1[1]*e2[2]+e1[2]*e2[1])+mm.m13*(e1[1]*e2[3]+e1[3]*e2[1])+
            mm.m23*(e1[2]*e2[3]+e1[3]*e2[2]),
            mm.m11*e2[1]*e2[1]+mm.m22*e2[2]*e2[2]+mm.m33*e2[3]*e2[3]+
            2*(mm.m12*e2[1]*e2[2]+mm.m13*e2[1]*e2[3]+mm.m23*e2[2]*e2[3]))
    end
    lift=(u,v)->_surface_metric_lift(frame,u,v)
    world=[lift(T.x[row],T.y[row]) for row in 1:Int(T.nreal)]
    return Mesh2D._GmshSurfaceMetricCarrier(world,lift,metric)
end

function _surface_algorithm_refine!(m,T,xs,ys,segs,markers,embedded,index,
                                     bgm_size,algorithm,t,caller;
                                     _gmsh_surface_trace=nothing)
    interior=Vector{Bool}(collect(classify_interior(T)))
    algorithm==3 && return interior
    embedded_ids=_surface_embedded_vertex_ids(embedded,index)
    sizes=_surface_vertex_sizes(T,xs,ys,segs,markers,embedded_ids;interior)
    vbgm=copy(sizes)
    embedded_vids=Int32[]
    if !isempty(embedded_ids)
        map=_pslg_vertex_map(T,xs,ys)
        loose=_pslg_default_size(xs,ys,caller,t)
        for pid in embedded_ids
            vid=map[_surface_size_key(xs[pid],ys[pid])]
            vbgm[vid]=bgm_size===nothing ? loose : bgm_size(xs[pid],ys[pid])
            push!(embedded_vids,vid)
        end
    end
    extend=_resolved_extend_from_boundary(m,2,t)
    return refine!(T;min_angle_deg=0,gmsh_insertion=true,
        gmsh_algorithm=algorithm==5 ? :delaunay : :frontal,
        gmsh_extend_boundary=extend>0 || extend==-2,
        vertex_sizes=sizes,vertex_bgm=vbgm,bgm_size=bgm_size,
        model_vertices=embedded_vids,_gmsh_surface_trace)
end

# Box factories carry their explicit encoding provenance in each typed face.
# Their established geometry representation retains the prior refinement
# protocol; this is separate from native Algorithm6 on user-assembled PLCs.
function _model_factory_surface_refine!(m,T,xs,ys,segs,markers,embedded,index,
                                        plane,size_field,t,caller)
    callback=m.meshing.size_callback
    axes=plane.axes;dropped=plane.k
    bgm_size=if size_field===nothing && callback===nothing
        nothing
    else
        loose=_pslg_default_size(xs,ys,caller,t)
        function background(x,y)
            size=loose
            point=ntuple(3) do axis
                axis==dropped ? _plane_dropped_coordinate(plane,x,y) :
                    axis==axes[1] ? x : y
            end
            size_field===nothing || (size=min(size,size_at(size_field,point...,(2,t))))
            callback===nothing || (size=_apply_size_callback(callback,2,t,point...,size,caller))
            return size
        end
    end
    extend=_resolved_extend_from_boundary(m,2,t)
    if extend>0 || extend==-2
        embedded_ids=_surface_embedded_vertex_ids(embedded,index)
        sizes=_surface_vertex_sizes(T,xs,ys,segs,markers,embedded_ids)
        background_sizes=nothing;embedded_vertices=nothing
        if !isempty(embedded_ids)
            mapping=_pslg_vertex_map(T,xs,ys)
            loose_vertex=_pslg_default_size(xs,ys,caller,t)
            bg=copy(sizes)
            for pid in embedded_ids
                vertex=mapping[_surface_size_key(xs[pid],ys[pid])]
                bg[vertex]=bgm_size===nothing ? loose_vertex : bgm_size(xs[pid],ys[pid])
            end
            background_sizes=bg
            embedded_vertices=[mapping[_surface_size_key(xs[pid],ys[pid])] for pid in embedded_ids]
        end
        return refine!(T;min_angle_deg=0,gmsh_insertion=true,
            vertex_sizes=sizes,vertex_bgm=background_sizes,bgm_size,
            model_vertices=embedded_vertices)
    end
    radius=bgm_size===nothing ? (x,y)->_pslg_default_size(xs,ys,caller,t) : bgm_size
    return refine!(T;min_angle_deg=0,size_radius=radius,gmsh_insertion=true)
end
