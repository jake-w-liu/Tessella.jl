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
# Upstream `Mesh.Algorithm=2` (AUTO) resolves per face: Delaunay on planes
# and MeshAdapt on every non-planar class — verified on the Gmsh 4.15.2
# binary, where Cylinder/Sphere/ruled Surface/BSpline faces all report
# "Meshing surface N (<kind>, MeshAdapt)". No native MeshAdapt producer
# exists, so the implemented Frontal-Delaunay construction serves the
# off-plane admission instead of rejecting: it is the closest native
# producer for upstream's curved-face outcome. On the planar PSLG path
# `:ruled`/`:tric` faces mesh through the `:frontal` arm; on OCC-curved
# faces the parametric kernel treats 5/6 identically.
function _normal_surface_mesh_algorithm(m::GeoModel,t::Int,caller;metric::Bool=false)
    algorithm=_surface_mesh_algorithm(m,t)
    algorithm==2 && (algorithm=_surface_type(m,t)==:plane ? 5 : 6)
    algorithm in (metric ? (3,5,6,7) : (3,5,6)) || throw(ArgumentError(
        "$caller: mesh algorithm $algorithm on Surface[$t] is not " *
        "implemented by the ordinary native surface mesher; supported: " *
        "2 (Auto), 3 (initial only), 5 (Delaunay), " *
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
# measure the metric-normalized circumradius and scalar local bounds fold in
# through `_gmsh_metric_floor`; `vbgm` itself stays scalarized to the
# strongest direction for the vertex-size records. A scalar or absent
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

# ── Non-planar built-in ruled/tric Surface fill ────────────────────────────
#
# A declared `Surface`/`RuledSurface` record whose boundary does not fit one
# plane cannot take the projected-PSLG path: upstream meshes these through the
# ruled patch's own `(u,v)` parameter domain (`parBounds` is the unit square
# for four generatrices, the `v<=u` triangle for three) and evaluates
# `GFace::point(u,v)` for every interior node, measuring distances through the
# first fundamental form. The construction below keeps the same contract with
# native pieces:
#
#   * boundary curve parameters map onto the patch edges through
#     `gmshEdge::reparamOnFace` (`_ruled_curve_parameters` — the direct
#     edge-coordinate map for generatrices, `XYZtoUV`/`parFromPoint` for
#     off-boundary and embedded curves);
#   * corner and embedded vertices map through `gmshVertex::reparamOnFace`
#     (`_ruled_vertex_parameters`);
#   * the triangulation runs in `(u,v)` with the metric-only carrier armed on
#     every triangle: `metric(u,v) = Jᵀ M_world(p(u,v)) J`, where `J` is the
#     ruled patch's `d1` Jacobian and `M_world` is the resolved background
#     tensor — an anisotropic field's `metric_at` pullback, or the isotropic
#     `I/h²` of the scalar background. The isotropic arm reduces to `G/h²`
#     (first fundamental form over local size squared), so a metric-normalized
#     distance of 1 is exactly one `h`-sized element measured along the
#     surface;
#   * vertex target sizes (`setLcs` analogues) measure incident edge lengths
#     in world coordinates, not parameter space;
#   * recombination angle/priority tests and the final lift both evaluate
#     `_ruled_surface_point`, while entity-derived vertices keep their
#     canonical bitwise curve-mesh coordinates exactly as the planar path
#     publishes them.
#
# Coplanar declared ruled surfaces keep the established planar projection
# path — their existing oracle outputs are pinned; only genuinely warped
# patches (the configuration the plane fit cannot represent) route here.

function _add_ruled_surface_point!(xs,ys,mesh_sizes,index,canonical,m::GeoModel,
                                   t::Int,pid::Int,lc::Float64,caller)
    haskey(index,pid) && return index[pid]
    haskey(m.points,pid) || throw(ArgumentError("$caller: unknown Point[$pid]"))
    silent=Returns(nothing)
    u,v=_ruled_vertex_parameters(m,t,pid,lc,false,silent,silent,caller,false)
    (isfinite(u) && isfinite(v)) || throw(ArgumentError(
        "$caller: Point[$pid] has no finite parameter on Surface[$t]"))
    push!(xs,u);push!(ys,v)
    push!(mesh_sizes,get(m.point_size,pid,0.0))
    canonical[_surface_size_key(u,v)]=m.points[pid]
    index[pid]=length(xs)
    return index[pid]
end

function _add_ruled_surface_curve_point!(xs,ys,mesh_sizes,index,canonical,
                                       m::GeoModel,t::Int,point,
                                       mesh_size::Float64,curve::Int,
                                       parameter::Float64,lc::Float64,caller)
    # A subdivision node that resolves to a model vertex keeps that vertex —
    # the bitwise-identical rule the planar `_add_surface_curve_point!`
    # applies before emitting a new row.
    scale=max(1.0,hypot(point...))
    tolerance=128eps(Float64)*scale
    matched_vertex=0
    matched_point=0
    for point_tag in sort!(collect(keys(index)))
        candidate=m.points[point_tag]
        distance=hypot(candidate[1]-point[1],candidate[2]-point[2],
                       candidate[3]-point[3])
        distance<=tolerance || continue
        vertex=index[point_tag]
        if matched_vertex!=0 && matched_vertex!=vertex
            throw(ArgumentError(
                "$caller: Curve[$curve] subdivision coincides with " *
                "Points[$matched_point] and [$point_tag]"))
        end
        matched_vertex=vertex
        matched_point=point_tag
    end
    if matched_vertex!=0
        existing=mesh_sizes[matched_vertex]
        mesh_sizes[matched_vertex]=
            existing>0 ? (mesh_size>0 ? min(existing,mesh_size) : existing) :
                         mesh_size
        return matched_vertex
    end
    silent=Returns(nothing)
    u,v=_ruled_curve_parameters(m,t,curve,parameter,lc,false,silent,silent,
                                caller,false)
    (isfinite(u) && isfinite(v)) || throw(ArgumentError(
        "$caller: Curve[$curve] subdivision has no finite parameter on " *
        "Surface[$t]"))
    push!(xs,u);push!(ys,v)
    push!(mesh_sizes,mesh_size)
    canonical[_surface_size_key(u,v)]=point
    return length(xs)
end

# `_surface_pslg` over the ruled parameter domain. The loop walk, forced-
# parameter enumeration, marker bookkeeping and segment/embedded records are
# identical to the planar builder; only the vertex map differs.
function _surface_ruled_pslg(m::GeoModel,t::Int,forced,caller::AbstractString;
                             param_sizes::Dict{Tuple{Int,Float64},Float64}=
                                 Dict{Tuple{Int,Float64},Float64}(),
                             _periodic_positions=nothing)
    lc=_model_lc(m)
    xs=Float64[];ys=Float64[];mesh_sizes=Float64[]
    segs=Tuple{Int,Int}[]
    index=Dict{Int,Int}()
    canonical=Dict{Tuple{Float64,Float64},NTuple{3,Float64}}()
    for loop_id in m.surfaces[t]
        _verify_loop_closed(m,loop_id,caller,"Surface[$t]")
        loop_idx=Int[]
        for signed in m.loops[loop_id]
            curve=abs(signed)
            a,b=m.curves[curve]
            t0,t1=_model_curve_param_bounds(m,curve,caller)
            for parameter in _surface_curve_parameters(m,forced,curve,signed,
                                                       caller)
                vertex=if parameter==t0
                    _add_ruled_surface_point!(xs,ys,mesh_sizes,index,canonical,
                                              m,t,a,lc,caller)
                elseif parameter==t1
                    _add_ruled_surface_point!(xs,ys,mesh_sizes,index,canonical,
                                              m,t,b,lc,caller)
                else
                    point=_periodic_forced_curve_point(
                        m,curve,parameter,_periodic_positions,caller)
                    _add_ruled_surface_curve_point!(
                        xs,ys,mesh_sizes,index,canonical,m,t,point,
                        get(param_sizes,(curve,parameter),
                            _surface_curve_mesh_size(
                                m,curve,parameter,caller)),
                        curve,parameter,lc,caller)
                end
                push!(loop_idx,vertex)
            end
        end
        nloop=length(loop_idx)
        nloop>=3 || throw(ArgumentError(
            "$caller: Loop[$loop_id] needs at least three points"))
        for k in 1:nloop
            i,j=loop_idx[k],loop_idx[mod1(k+1,nloop)]
            i==j || (xs[i]==xs[j] && ys[i]==ys[j]) || push!(segs,(i,j))
        end
    end
    embedded=get(m.embeds,(2,t),NTuple{2,Int}[])
    internal=Tuple{Int,Int}[]
    for (edim,etag) in embedded
        if edim==0
            _add_ruled_surface_point!(xs,ys,mesh_sizes,index,canonical,m,t,
                                      etag,lc,caller)
        elseif edim!=1
            throw(ArgumentError(
                "$caller: unsupported embedding dimension $edim"))
        end
    end
    for (edim,etag) in embedded
        edim==1 || continue
        haskey(m.curves,etag) || throw(ArgumentError(
            "$caller: unknown embedded Curve[$etag]"))
        a,b=m.curves[etag]
        t0,t1=_model_curve_param_bounds(m,etag,caller)
        parameters=etag in m.meshing.degenerated ? nothing :
                   get(forced,etag,nothing)
        curve_nodes=Int[]
        if parameters===nothing
            push!(curve_nodes,_add_ruled_surface_point!(
                xs,ys,mesh_sizes,index,canonical,m,t,a,lc,caller))
            push!(curve_nodes,_add_ruled_surface_point!(
                xs,ys,mesh_sizes,index,canonical,m,t,b,lc,caller))
        else
            for parameter in parameters
                vertex=if parameter==t0
                    _add_ruled_surface_point!(xs,ys,mesh_sizes,index,canonical,
                                              m,t,a,lc,caller)
                elseif parameter==t1
                    _add_ruled_surface_point!(xs,ys,mesh_sizes,index,canonical,
                                              m,t,b,lc,caller)
                else
                    point=_periodic_forced_curve_point(
                        m,etag,parameter,_periodic_positions,caller)
                    _add_ruled_surface_curve_point!(
                        xs,ys,mesh_sizes,index,canonical,m,t,point,
                        get(param_sizes,(etag,parameter),
                            _surface_curve_mesh_size(
                                m,etag,parameter,caller)),
                        etag,parameter,lc,caller)
                end
                push!(curve_nodes,vertex)
            end
        end
        length(curve_nodes)>=2 || throw(ArgumentError(
            "$caller: embedded Curve[$etag] needs two distinct endpoints"))
        for segment_index in 1:(length(curve_nodes)-1)
            first_node=curve_nodes[segment_index]
            second_node=curve_nodes[segment_index+1]
            (first_node==second_node ||
             (xs[first_node]==xs[second_node] &&
              ys[first_node]==ys[second_node])) && throw(ArgumentError(
                "$caller: embedded Curve[$etag] has coincident subdivision nodes"))
            push!(internal,(first_node,second_node))
        end
    end
    return xs,ys,mesh_sizes,segs,embedded,internal,index,canonical
end

# `bgm_size(u,v)` on the ruled domain — the resolved model background
# evaluated at the lifted world point, the same contract the planar
# `_ModelSurfaceBackgroundSize` provides through its affine lift. The
# background member stays parametric: `_ModelAlgorithmBackgroundSize` is
# defined in a later include.
struct _ModelSurfaceParamBackgroundSize{F,B}
    lift::F
    background::B
end

function (background::_ModelSurfaceParamBackgroundSize)(u,v)
    point=background.lift(u,v)
    return background.background(point...)
end

# The patch Jacobian: `d1` finite-difference columns of `_ruled_surface_point`.
@inline function _ruled_jacobian(m::GeoModel,t::Int,u::Float64,v::Float64,
                                 caller)
    return (_ruled_surface_d1(m,t,u,v,1,caller),
            _ruled_surface_d1(m,t,u,v,2,caller))
end

@inline function _ruled_metric_pullback(mm,du,dv)
    return (mm.m11*du[1]*du[1]+mm.m22*du[2]*du[2]+mm.m33*du[3]*du[3]+
            2*(mm.m12*du[1]*du[2]+mm.m13*du[1]*du[3]+mm.m23*du[2]*du[3]),
            mm.m11*du[1]*dv[1]+mm.m22*du[2]*dv[2]+mm.m33*du[3]*dv[3]+
            mm.m12*(du[1]*dv[2]+du[2]*dv[1])+mm.m13*(du[1]*dv[3]+du[3]*dv[1])+
            mm.m23*(du[2]*dv[3]+du[3]*dv[2]),
            mm.m11*dv[1]*dv[1]+mm.m22*dv[2]*dv[2]+mm.m33*dv[3]*dv[3]+
            2*(mm.m12*dv[1]*dv[2]+mm.m13*dv[1]*dv[3]+mm.m23*dv[2]*dv[3]))
end

# Metric-only carrier for the ruled fill, armed on every scalar or tensor
# background — the carrier is what makes cavity, circumcenter, queue radius
# and frontal-step tests read surface distances instead of raw `(u,v)`
# lengths. `world` covers the birth rows like `_gmsh_surface_trace`: boundary
# rows resolve through `canonical` so entity-derived coordinates stay bitwise
# identical to the owning curve meshes.
function _surface_ruled_trace_carrier(m,options,size_field,lift,T,t,caller,
                                      canonical)
    raw=options.inner_field===nothing ? size_field : options.inner_field
    entity=(2,t)
    metric=if raw isa AbstractAnisoField
        lc_min=isnan(options.lc_min) ? 0.0 : options.lc_min
        lc_max=isnan(options.lc_max) ? GMSH_MAX_SIZE : options.lc_max
        lo=min(lc_min,lc_max)
        hi=max(lo,min(lc_max,options.ctx_lc))
        wrapped=BoundedSize(raw;size_min=lo,size_max=hi,
                            factor=options.lc_factor)
        function(u,v)
            w=lift(u,v)
            mm=metric_at(wrapped,w[1],w[2],w[3],entity)
            du,dv=_ruled_jacobian(m,t,u,v,caller)
            return _ruled_metric_pullback(mm,du,dv)
        end
    else
        background=_ModelAlgorithmBackgroundSize(m,2,t,options,size_field,
                                                 caller)
        function(u,v)
            w=lift(u,v)
            h=background(w...)
            du,dv=_ruled_jacobian(m,t,u,v,caller)
            inv2=1.0/(h*h)
            return (inv2*Mesh2D._gmsh_dot3(du,du),
                    inv2*Mesh2D._gmsh_dot3(du,dv),
                    inv2*Mesh2D._gmsh_dot3(dv,dv))
        end
    end
    world=Vector{NTuple{3,Float64}}(undef,Int(T.nreal))
    @inbounds for row in 1:Int(T.nreal)
        key=_surface_size_key(T.x[row],T.y[row])
        stored=get(canonical,key,nothing)
        world[row]=stored===nothing ? lift(T.x[row],T.y[row]) : stored
    end
    # UV units are parameter units here, not world lengths: the carrier's
    # `ground` supplies `G=JᵀJ` so a scalar `vs` bound folds through
    # `G/vs²` — the actual `I/vs²` pullback — instead of the world-aligned
    # eigenvalue floor `_gmsh_metric_floor` applies on planar domains.
    ground=function(u,v)
        du,dv=_ruled_jacobian(m,t,u,v,caller)
        return (Mesh2D._gmsh_dot3(du,du),
                Mesh2D._gmsh_dot3(du,dv),
                Mesh2D._gmsh_dot3(dv,dv))
    end
    return Mesh2D._GmshSurfaceMetricCarrier(world,lift,metric,ground)
end

# `setLcs` analogue on the warped patch: incident edge lengths are the world
# chords of the initial triangulation (`world` rows carry the canonical
# entity coordinates, matching the vertices the final mesh publishes).
function _ruled_vertex_sizes(T,xs,ys,floor_segs,markers,embedded_ids,world;
                             interior=nothing)
    vid_of=_pslg_vertex_map(T,xs,ys)
    vid(i::Int)=vid_of[_surface_size_key(xs[i],ys[i])]
    n=T.nreal
    sizes=fill(Inf,n)
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tri(T,t)) || continue
        interior===nothing || interior[t] || continue
        for (i,j) in ((1,2),(2,3),(3,1))
            a=_vert(T,t,i); b=_vert(T,t,j)
            (_is_ghost_v(a) || _is_ghost_v(b)) && continue
            wa=world[Int(a)];wb=world[Int(b)]
            len=hypot(wb[1]-wa[1],wb[2]-wa[2],wb[3]-wa[3])
            len>0 || continue
            len<sizes[a] && (sizes[a]=len)
            len<sizes[b] && (sizes[b]=len)
        end
    end
    @inbounds for (a,b) in floor_segs
        va=vid(a);vb=vid(b)
        wa=world[Int(va)];wb=world[Int(vb)]
        half=0.5*hypot(wb[1]-wa[1],wb[2]-wa[2],wb[3]-wa[3])
        sizes[va]<half && (sizes[va]=half)
        sizes[vb]<half && (sizes[vb]=half)
    end
    @inbounds for vertex in embedded_ids
        v=vid(vertex)
        prescribed=markers[vertex]
        (isfinite(prescribed) && prescribed>0 && prescribed<sizes[v]) &&
            (sizes[v]=prescribed)
    end
    if any(isinf,sizes)
        xmin=ymin=zmin=Inf
        xmax=ymax=zmax=-Inf
        @inbounds for w in world
            w[1]<xmin && (xmin=w[1]);w[1]>xmax && (xmax=w[1])
            w[2]<ymin && (ymin=w[2]);w[2]>ymax && (ymax=w[2])
            w[3]<zmin && (zmin=w[3]);w[3]>zmax && (zmax=w[3])
        end
        fallback=hypot(xmax-xmin,ymax-ymin,zmax-zmin)
        (!isfinite(fallback) || fallback<=0) && (fallback=1.0)
        @inbounds for vertex in eachindex(sizes)
            isinf(sizes[vertex]) && (sizes[vertex]=fallback)
        end
    end
    return sizes
end

# `_surface_algorithm_refine!` over the ruled parameter domain — identical
# structure, with `setLcs`-world sizes and the metric carrier supplied by the
# caller.
function _surface_algorithm_ruled_refine!(m,T,xs,ys,segs,markers,embedded,
                                         index,bgm_size,algorithm,t,caller,
                                         world;
                                         _gmsh_surface_trace=nothing)
    interior=Vector{Bool}(collect(classify_interior(T)))
    algorithm==3 && return interior
    embedded_ids=_surface_embedded_vertex_ids(embedded,index)
    sizes=_ruled_vertex_sizes(T,xs,ys,segs,markers,embedded_ids,world;
                              interior)
    vbgm=copy(sizes)
    embedded_vids=Int32[]
    if !isempty(embedded_ids)
        map=_pslg_vertex_map(T,xs,ys)
        loose=_ruled_pslg_default_size(world,caller,t)
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

# Characteristic length of the lifted boundary, the `_pslg_default_size`
# analogue measured on world rows.
function _ruled_pslg_default_size(world,caller::AbstractString,t::Int)
    isempty(world) && throw(ArgumentError(
        "$caller: Surface[$t] has no boundary points"))
    xmin=ymin=zmin=Inf
    xmax=ymax=zmax=-Inf
    @inbounds for w in world
        w[1]<xmin && (xmin=w[1]);w[1]>xmax && (xmax=w[1])
        w[2]<ymin && (ymin=w[2]);w[2]>ymax && (ymax=w[2])
        w[3]<zmin && (zmin=w[3]);w[3]>zmax && (zmax=w[3])
    end
    lc=hypot(xmax-xmin,ymax-ymin,zmax-zmin)
    (isfinite(lc) && lc>0) || throw(ArgumentError(
        "$caller: Surface[$t] has a degenerate bounding box"))
    return lc
end

function _ruled_param_physical_coords(mesh,lift,canonical)
    coords=Matrix{Float64}(undef,3,size(mesh.coords,2))
    @inbounds for node in axes(mesh.coords,2)
        u=mesh.coords[1,node];v=mesh.coords[2,node]
        lifted=get(canonical,_surface_size_key(u,v),nothing)
        point=lifted===nothing ? lift(u,v) : lifted
        coords[1,node]=point[1]
        coords[2,node]=point[2]
        coords[3,node]=point[3]
    end
    return coords
end

# The ordinary fill for a non-planar built-in `Surface`: constrained
# Delaunay of the ruled `(u,v)` domain, Gmsh-family refinement under the
# patch metric, upstream's recombines-before-attributes ordering, then the
# `_ruled_surface_point` lift — the direct counterpart of the planar fallback
# in `_mesh_model_surface_once`.
function _model_ruled_surface_mesh(m::GeoModel,t::Int,forced,
                                   caller::AbstractString;
                                   size_field::Union{Nothing,AbstractSizeField}=
                                       nothing,
                                   _mesh1d_options=nothing,
                                   param_sizes::Dict{Tuple{Int,Float64},
                                                     Float64}=
                                       Dict{Tuple{Int,Float64},Float64}(),
                                   _periodic_positions=nothing,
                                   algorithm::Int=6)
    xs,ys,markers,segs,embedded,internal,index,canonical=
        _surface_ruled_pslg(m,t,forced,caller;param_sizes=param_sizes,
            _periodic_positions)
    T=constrained_delaunay(xs,ys,segs;internal_segments=internal)
    options=_model_algorithm_size_options(m,_mesh1d_options,caller)
    lift=(u,v)->_ruled_surface_point(m,t,u,v,caller)
    bgm_size=_ModelSurfaceParamBackgroundSize(lift,
        _ModelAlgorithmBackgroundSize(m,2,t,options,size_field,caller))
    trace=_surface_ruled_trace_carrier(m,options,size_field,lift,T,t,caller,
                                       canonical)
    world=trace.world
    interior=_surface_algorithm_ruled_refine!(
        m,T,xs,ys,segs,markers,embedded,index,bgm_size,algorithm,t,caller,
        world;_gmsh_surface_trace=trace)
    mesh=to_mesh(T; interior=interior)
    if _model_surface_recombined(m,t)
        algo=m.meshing.recombine_algo
        algo>=2 && throw(ArgumentError(
            "$caller: Mesh.RecombinationAlgorithm=$algo requests the " *
            "full-quad subdivision pipeline on Surface[$t], which is not " *
            "implemented — supported values are 0 (simple) and 1 " *
            "(blossom)"))
        angle=get(m.meshing.recombine,(2,t),45.0)
        physical_coords=_ruled_param_physical_coords(mesh,lift,canonical)
        mesh=recombine_triangles(
            mesh;
            algorithm=algo==0 ? :greedy : :blossom,
            recombine_angle=angle,
            recombine_angle_coords=physical_coords,
            recombine_priority_coords=physical_coords,
            protected_edges=_embedded_protected_edges(mesh,xs,ys,internal))
    end
    mesh=_consume_surface_attributes(m,t,mesh,caller)
    @inbounds for node in axes(mesh.coords,2)
        u=mesh.coords[1,node];v=mesh.coords[2,node]
        lifted=get(canonical,_surface_size_key(u,v),nothing)
        point=lifted===nothing ? lift(u,v) : lifted
        mesh.coords[1,node]=point[1]
        mesh.coords[2,node]=point[2]
        mesh.coords[3,node]=point[3]
    end
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    _surface_cell_count(mesh)>0 || throw(ErrorException(
        "$caller: Surface[$t] produced no cells"))
    return mesh,embedded
end
