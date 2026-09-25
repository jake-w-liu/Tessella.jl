# ============================================================================
# Model-level 1-D mesh generation — Gmsh `Mesh0D` + `Mesh1D` over `GeoModel`
# curves. `_model_mesh_curves!` fills `m.curve_params` (each curve's stored
# parameter discretization, mirroring `GEdge::mesh_vertices`+`lines`);
# `_model_point_mesh_parts`/`_model_curve_mesh_parts` turn the stored
# parameters into per-entity `(dim, tag, Mesh)` triples for
# `_geo_merge_entity_meshes`. The size field reproduces `BGM_MeshSize` for
# dim 0/1 entities: `min(CTX::lc, l1 points, l2 curvature, l3 background
# field, l4 entity size, l5 parametric sizes, lc-callback)` clamped to
# `[MeshSizeMin, MeshSizeMax]` and scaled by `lcFactor` — with the `F_Lc`
# endpoint rule `min(BGM(vertex), BGM(curve))` delivered through the kernel's
# `endpoint_entities` channel.
# ============================================================================

import ..SizeField: field_value

# Mesh-time option bundle — values mirrored from the `.geo`/API option store
# each time a model mesh runs, mirroring upstream's global `CTX::mesh` reads.
struct _ModelMesh1DOptions
    ctx_lc::Float64            # padded model-bounds diagonal (CTX::lc)
    lc_min::Float64            # MeshSizeMin / CharacteristicLengthMin
    lc_max::Float64            # MeshSizeMax / CharacteristicLengthMax
    lc_factor::Float64         # MeshSizeFactor / CharacteristicLengthFactor
    from_points::Bool          # MeshSizeFromPoints
    from_curvature::Float64    # MeshSizeFromCurvature (elements per 2π; 0 = off)
    integration_precision::Float64
    min_line_nodes::Int
    min_circle_nodes::Int
    min_curve_nodes::Int
    tolerance_edge_length::Float64
    geometry_tolerance::Float64
    mesh_only_empty::Bool
    mesh_only_visible::Bool
    max_retries::Int
    inner_field::Any           # background AbstractSizeField or nothing
    callback::Any              # m.meshing.size_callback or nothing
    warn::Any                  # Msg::Warning sink or nothing
    report_error::Any          # Msg::Error sink or nothing
end

function _model_mesh_warn(options::_ModelMesh1DOptions, message::AbstractString)
    options.warn===nothing || options.warn(message)
    return nothing
end

function _model_mesh_error(options::_ModelMesh1DOptions,
                           message::AbstractString)
    options.report_error===nothing || options.report_error(message)
    return nothing
end

# ── CTX::lc at mesh time ─────────────────────────────────────────────────────
# `SetBoundingBox` = `GModel::bounds()` then `FinishUpBoundingBox`. The vertex
# cloud covers model points, generated curve nodes (stored `curve_params` —
# upstream bounds() includes generated mesh vertices), discrete record nodes
# and the stored OCC primitive extents; empty models use upstream's [-1,1]^3
# placeholder box.
function _model_mesh_bbox(m::GeoModel, caller::AbstractString)
    lo=fill(Inf,3); hi=fill(-Inf,3)
    function grow!(p)
        for axis in 1:3
            v=p[axis]
            isfinite(v) || throw(ArgumentError(
                "$caller: model bounds hit a non-finite coordinate"))
            v<lo[axis] && (lo[axis]=v)
            v>hi[axis] && (hi[axis]=v)
        end
        return nothing
    end
    for point in values(m.points)
        grow!(point)
    end
    for (curve,params) in m.curve_params
        haskey(m.curves,curve) || continue
        for u in params
            grow!(_model_curve_point(m,curve,u,caller))
        end
    end
    for (_,record) in _discrete_mesh_records_model(m)
        for i in axes(record.node_coords,2)
            grow!((record.node_coords[1,i],record.node_coords[2,i],
                   record.node_coords[3,i]))
        end
    end
    for (xmin,xmax,ymin,ymax,zmin,zmax) in values(m.box_extents)
        grow!((xmin,ymin,zmin)); grow!((xmax,ymax,zmax))
    end
    # Finite-primitive extents — the same vertex-cloud role OCC geometric
    # bounds play upstream (a free cylinder/sphere/cone must not shrink lc).
    for cyl in values(m.cylinders)
        c=cyl.center; ax=cyl.axis; r=cyl.radius
        n=norm((ax[1],ax[2],ax[3]))
        n>0 || continue
        u=(ax[1]/n,ax[2]/n,ax[3]/n); h=cyl.height/2
        for s in (-1,1)
            grow!(ntuple(i->c[i]+s*(h*abs(u[i])+r*sqrt(max(0.0,1-u[i]^2))),3))
        end
    end
    for sph in values(m.spheres)
        c=sph.center; r=sph.radius
        grow!((c[1]-r,c[2]-r,c[3]-r)); grow!((c[1]+r,c[2]+r,c[3]+r))
    end
    for cone in values(m.cones)
        c=cone.center; ax=cone.axis; n=norm((ax[1],ax[2],ax[3]))
        n>0 || continue
        u=(ax[1]/n,ax[2]/n,ax[3]/n); h=cone.height/2
        for s in (-1,1)
            r=s<0 ? cone.r1 : cone.r2
            grow!(ntuple(i->c[i]+s*(h*abs(u[i])+r*sqrt(max(0.0,1-u[i]^2))),3))
        end
    end
    isfinite(lo[1]) || return (-1.0,1.0,-1.0,1.0,-1.0,1.0)
    return (lo[1],hi[1],lo[2],hi[2],lo[3],hi[3])
end

function _model_mesh_lc(m::GeoModel, geometry_tolerance::Float64,
                      caller::AbstractString)
    bbox=_model_mesh_bbox(m,caller)
    # `FinishUpBoundingBox` padding + diagonal — shared with the field builder.
    return _gmsh_bbox_characteristic_length(bbox,geometry_tolerance)
end

# ── BGM_MeshSize for dim 0/1 entities ────────────────────────────────────────
# `prescribedMeshSizeAtVertex`: stored 0 means unsized → MAX_LC upstream
# (`addVertex` maps `!lc` to MAX_LC); negative values propagate like upstream.
function _model_vertex_lc(m::GeoModel,point::Int)
    lc=get(m.point_size,point,0.0)
    return lc==0.0 ? GMSH_MAX_SIZE : lc
end

# `LC_MVertex_PNTS` for a vertex: its prescribed size, or `CTX::lc/10` when
# unsized.
function _model_vertex_l1(m::GeoModel,point::Int,ctx_lc::Float64)
    lc=_model_vertex_lc(m,point)
    return lc>=GMSH_MAX_SIZE ? ctx_lc/10 : lc
end

# `LC_MVertex_PNTS` for an edge: `CTX::lc/10` when both endpoints are unsized,
# else the linear interpolation `(1-a)lc1 + a*lc2` over the parameter range.
function _model_curve_l1(m::GeoModel,curve::Integer,u::Float64,t0::Float64,
                         t1::Float64,ctx_lc::Float64)
    a,b=m.curves[curve]
    lc1=_model_vertex_lc(m,a); lc2=_model_vertex_lc(m,b)
    (lc1>=GMSH_MAX_SIZE && lc2>=GMSH_MAX_SIZE) && return ctx_lc/10
    alpha=t1==t0 ? 0.0 : (u-t0)/(t1-t0)
    return (1-alpha)*lc1 + alpha*lc2
end

# `LC_MVertex_CURV` for an edge: `2π/(κ·ne)` from the curve's own curvature.
# The upstream `max_surf_curvature` term needs adjacent-face curvature, which
# Tessella's built-in surface evaluation does not provide — plane faces
# contribute zero, matching this truncation.
function _model_curve_l2(m::GeoModel,curve::Integer,u::Float64,ne::Float64)
    ne>0 || return GMSH_MAX_SIZE
    ne<1 && (ne=1.0)
    curvature=only(model_curvature(m,1,curve,[u]))
    return curvature>0 ? 2π/curvature/ne : GMSH_MAX_SIZE
end

# `LC_MVertex_CURV` for a vertex: max incident-edge curvature at the
# endpoint touching the vertex (`max_edge_curvature`).
function _model_vertex_l2(m::GeoModel,point::Int,ne::Float64)
    ne>0 || return GMSH_MAX_SIZE
    ne<1 && (ne=1.0)
    curvature=0.0
    for (curve,(a,b)) in m.curves
        a==point || b==point || continue
        # Lines have no curvature — and a coincident-endpoint `Line` cannot
        # be evaluated at all.
        get(m.curve_types,curve,:line)===:line &&
            !haskey(m.curve_geometry,curve) && continue
        lo,hi=_model_curve_param_bounds(m,curve,
                                        "_model_vertex_l2")
        # `max_edge_curvature` — a closed curve touches the vertex at both
        # ends, so both parameter extrema contribute.
        a==point &&
            (c=only(model_curvature(m,1,curve,[lo]));
             c>curvature && (curvature=c))
        b==point &&
            (c=only(model_curvature(m,1,curve,[hi]));
             c>curvature && (curvature=c))
    end
    return curvature>0 ? 2π/curvature/ne : GMSH_MAX_SIZE
end

# Native parameter bounds that tolerate degenerate curves: a `Line` is always
# [0,1]-parametrized, and `model_parametrization_bounds` throws on coincident
# endpoints (which are legal, meshable zero-length curves).
function _model_curve_param_bounds(m::GeoModel,curve::Integer,
                                   caller::AbstractString)
    occ=get(m.curve_geometry,curve,nothing)
    (occ!==nothing && hasproperty(occ,:t0)) &&
        return Float64(occ.t0),Float64(occ.t1)
    get(m.curve_types,curve,:line)===:line && return 0.0,1.0
    lo,hi=model_parametrization_bounds(m,1,curve)
    return Float64(lo[1]),Float64(hi[1])
end

# `prescribedMeshSizeAtParam`: linear interpolation over the sorted stored
# (u, lc) pairs; clamps to the first entry below the range and linearly
# extrapolates past the last like upstream.
function _model_curve_l5(m::GeoModel,curve::Integer,u::Float64)
    entries=get(m.meshing.size_at_params,(1,curve),nothing)
    (entries===nothing || isempty(entries)) && return GMSH_MAX_SIZE
    pairs=sort!(Tuple{Float64,Float64}[(u,lc) for (plist,lc) in entries
                                       for u in plist]; by=first)
    us=Float64[p[1] for p in pairs]
    n=length(us)
    it=searchsortedfirst(us,u)          # first index with us[it] >= u (n+1 past end)
    i1=min(it-1,n-1)                    # 0-based lower_bound position, capped
    i0=max(1,i1)-1
    u0,l0=pairs[i0+1]; u1,l1=pairs[i1+1]
    (i1==i0 || u0==u1) && return l0
    alpha=(u-u0)/(u1-u0)
    return l0*(1-alpha)+l1*alpha
end

# `BGM_MeshSize` for a dim-0 entity (the `F_Lc` endpoint term): l1 vertex
# size, l2 incident-edge curvature, l3 background field, callback, then the
# CTX::lc cap, lcMin/lcMax clamp and lcFactor scaling.
function _model_vertex_bgm(m::GeoModel,point::Int,xyz,
                           options::_ModelMesh1DOptions,caller::AbstractString)
    l1=options.from_points ? _model_vertex_l1(m,point,options.ctx_lc) :
        GMSH_MAX_SIZE
    l2=_model_vertex_l2(m,point,options.from_curvature)
    l3=options.inner_field===nothing ? GMSH_MAX_SIZE :
        field_value(options.inner_field,xyz[1],xyz[2],xyz[3],(0,point))
    lc=min(l1,l2,l3,GMSH_MAX_SIZE)
    if options.callback!==nothing
        lc=_apply_size_callback(options.callback,0,point,xyz[1],xyz[2],
                                xyz[3],lc,caller)
    end
    lc=min(options.ctx_lc,lc)
    lc=max(lc,options.lc_min)
    lc=min(lc,options.lc_max)
    if lc<=0
        _model_mesh_error(options,
            "Wrong mesh element size lc = $lc on point $point")
        lc=options.ctx_lc
    end
    return lc*options.lc_factor
end

# `BGM_MeshSize` for a dim-1 entity at parameter `u` (the `F_Lc`/`F_LcB`
# curve term). `l4` (entity mesh size) is MAX_LC — `.geo` cannot set a
# per-curve size — and the per-entity `meshSizeFactor` is likewise unset.
function _model_curve_bgm(m::GeoModel,curve::Integer,u::Float64,xyz,t0::Float64,
                          t1::Float64,options::_ModelMesh1DOptions,
                          caller::AbstractString)
    l1=options.from_points ?
        _model_curve_l1(m,curve,u,t0,t1,options.ctx_lc) : GMSH_MAX_SIZE
    l2=_model_curve_l2(m,curve,u,options.from_curvature)
    l3=options.inner_field===nothing ? GMSH_MAX_SIZE :
        field_value(options.inner_field,xyz[1],xyz[2],xyz[3],(1,curve))
    l5=_model_curve_l5(m,curve,u)
    lc=min(l1,l2,l3,GMSH_MAX_SIZE,l5)
    if options.callback!==nothing
        lc=_apply_size_callback(options.callback,1,curve,xyz[1],xyz[2],
                                xyz[3],lc,caller)
    end
    lc=min(options.ctx_lc,lc)
    lc=max(lc,options.lc_min)
    lc=min(lc,options.lc_max)
    if lc<=0
        _model_mesh_error(options,
            "Wrong mesh element size lc = $lc on curve $curve")
        lc=options.ctx_lc
    end
    return lc*options.lc_factor
end

# Entity-aware size-field adapter passed to `mesh_curve` as the `field`
# argument: dim-0 entities evaluate `BGM(vertex)` for the `F_Lc` endpoint
# rule; curve/unknown contexts fall back to the background field (the
# driver's `size_function` supplies the real curve evaluation).
struct _ModelCurveEntityField <: AbstractSizeField
    m::GeoModel
    options::_ModelMesh1DOptions
    caller::String
end

function field_value(f::_ModelCurveEntityField,x,y,z)
    f.options.inner_field===nothing && return GMSH_MAX_SIZE
    return field_value(f.options.inner_field,x,y,z)
end
field_value(f::_ModelCurveEntityField,x,y,z,::Nothing)=
    field_value(f,x,y,z)
function field_value(f::_ModelCurveEntityField,x,y,z,
                     entity::Tuple{T,U}) where {T<:Integer,U<:Integer}
    dim=Int(entity[1]); tag=Int(entity[2])
    dim==0 && return _model_vertex_bgm(f.m,tag,(x,y,z),f.options,f.caller)
    return field_value(f,x,y,z)
end

# ── Per-curve meshing ────────────────────────────────────────────────────────
# `.geo Degenerated` sets the record flag `degenerate(1)`: meshGEdge marks the
# curve DONE at entry and keeps any existing elements. OCC `:degenerate`
# records are `degenerate(0)`, which still emits the single endpoint edge
# (`N=1` inside `meshGEdgeProcessing`).
function _model_curve_degenerate_kind(m::GeoModel,curve::Integer)
    curve in m.meshing.degenerated && return 1
    get(m.curve_types,curve,:line)===:degenerate && return 0
    return -1
end

# `gmshEdge`/`OCCEdge::minimumMeshSegments` — the curve-specific segment
# floor: `MinLineNodes-1` for lines, arc-fraction of `MinCircleNodes` for
# circles/ellipses, `MinCurveNodes-1` otherwise; OCC closed edges force ≥4.
function _model_minimum_curve_segments(m::GeoModel,curve::Integer,
                                       options::_ModelMesh1DOptions,
                                       caller::AbstractString)
    kind=get(m.curve_types,curve,:line)
    occ=get(m.curve_geometry,curve,nothing)
    is_occ=occ!==nothing && hasproperty(occ,:occ)
    if kind===:line || (is_occ && occ.occ===:line)
        np=options.min_line_nodes-1
    elseif kind in (:circle,:ellipse) ||
           (is_occ && occ.occ in (:circle,:ellipse))
        a=if is_occ
            abs(occ.t1-occ.t0)
        else
            g=_arc_geometry(m,curve,caller)
            abs(g.t1-g.t2)
        end
        n=Float64(options.min_circle_nodes)
        np=a>6.28 ? trunc(Int,n) : trunc(Int,0.99+(n-1)*a/(2π))
    else
        np=options.min_curve_nodes-1
    end
    # `OCCEdge`: a closed edge generates at least 4 segments (one degenerate).
    if is_occ && haskey(m.curves,curve)
        a,b=m.curves[curve]
        a==b && (np=max(4,np))
    end
    # A curve bounding a two-generatrix surface needs at least 2 segments —
    # `gmshEdge::minimumMeshSegments`'s `Generatrices==2` rule (the
    # generatrices are the total boundary-curve count across the face's
    # loops).
    if np<2
        for face in _model_curve_faces(m,curve)
            sum(l->length(m.loops[l]),m.surfaces[face];init=0)==2 &&
                (np=2;break)
        end
    end
    return np
end

# A parametrized discrete curve re-meshes like `discreteEdge::mesh` (frames
# carry the stored parameters); records without parameters keep their
# elements, like `_discretization.empty()` upstream.
function _model_discrete_curve_parametrized(record::DiscreteEntity)
    return size(record.node_params,1)>=1 && size(record.node_params,2)>0
end

function _model_discrete_curve_params(m::GeoModel,curve::Integer,
                                      record::DiscreteEntity,
                                      caller::AbstractString)
    _model_discrete_curve_parametrized(record) ||
        return nothing
    us=sort!(unique(Float64.(vec(record.node_params[1,:]))))
    lo,hi=_model_curve_param_bounds(m,curve,caller)
    isempty(us) && return nothing
    return Float64[lo;filter(u->lo<u<hi,us);hi]
end

# The `degenerate(0)`/zero-length single-edge result (upstream `N=1` → one
# MLine between the two endpoints).
function _model_curve_single_edge_params(m::GeoModel,curve::Integer,
                                         caller::AbstractString)
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    return Float64[t0,t1]
end

# Extruded-curve params (`ExtrudeMesh` is only set when `Layers` ran):
# point-extruded generatrices take the cumulative layer heights
# (`ep->u(j,k+1)` scaled into the parameter range); face/volume top copies
# take their source curve's interior parameters like `copyMesh`, mirrored
# into the copy's own parameter frame when the source record is signed.
function _model_curve_extrude_params(m::GeoModel,curve::Integer,
                                     params::_GeoExtrudeParams,
                                     options::_ModelMesh1DOptions,
                                     caller::AbstractString)
    src=get(m.meshing.extrude_sources,(1,curve),nothing)
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    if src===nothing
        # An attached extrusion record with no recorded source (hand-built
        # model): the point-generatrix convention — upstream `ExtrudePoint`
        # wires the source vertex as the curve's begin vertex.
        haskey(m.curves,curve) || throw(ArgumentError(
            "$caller: Curve[$curve] extrusion has no recorded source"))
        src=(0,first(m.curves[curve]))
    end
    if src[1]==0
        # `extrudeMesh`: interior nodes at `ep->u(j,k+1)`, skipping the last
        # element's last node (the end vertex owns it).
        heights=params.heights;layers=params.layers
        length(heights)==length(layers) || throw(ArgumentError(
            "$caller: Curve[$curve] extrusion layers/heights length mismatch"))
        values=Float64[t0]
        for j in eachindex(layers)
            h0=j==1 ? 0.0 : heights[j-1]
            h1=heights[j]
            for k in 1:layers[j]
                (j==length(layers) && k==layers[j]) && continue
                u=h0+k/layers[j]*(h1-h0)
                push!(values,t0+u*(t1-t0))
            end
        end
        push!(values,t1)
        return values
    end
    # `copyMesh` — `direction = geo.Source > 0 ? 1 : -1`; a pending source
    # retries the curve like upstream's PENDING status.
    source=abs(src[2])
    haskey(m.curves,source) || haskey(m.discrete,(1,source)) ||
        throw(ArgumentError(
            "$caller: Curve[$curve] extrusion source $source is missing"))
    src_params=get(m.curve_params,source,nothing)
    if src_params===nothing && haskey(m.discrete,(1,source))
        src_params=_model_discrete_curve_params(
            m,source,m.discrete[(1,source)],caller)
    end
    src_params===nothing && return :pending
    interior=length(src_params)>2 ? src_params[2:end-1] : Float64[]
    if src[2]>0
        return Float64[t0;interior;t1]
    end
    # Reversed source: upstream stores `newu = u_max - u + u_min`. Tessella's
    # reversed copies mirror the parameter interval itself, so the equivalent
    # stored value on the copy is `v = t1 - (u - u_min)` — identical for
    # [0,1] ranges, and consistent for OCC sign-flipped and NURBS
    # knot-mirrored bounds. Upstream pushes them ascending.
    u_min=_model_curve_param_bounds(m,source,caller)[1]
    mapped=reverse!(t1 .- interior .+ u_min)
    return Float64[t0;mapped;t1]
end

# Periodic slave copy — `copyMesh(gef, ge, masterOrientation)`: forward
# relations shift the master's parameters by the bound difference; reversed
# relations mirror them inside the master's bounds (upstream stores
# `newu = u_max - u + to_u_min`, pushed ascending).
function _model_curve_periodic_params(m::GeoModel,curve::Integer,
                                      master::Integer,master_params,
                                      reversed::Bool,
                                      caller::AbstractString)
    isempty(master_params) && return Float64[]
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    m0,m1=_model_curve_param_bounds(m,master,caller)
    if reversed
        # `newu = u_max - u + to_u_min` — descending for an ascending master
        # list, so push them back in upstream's ascending order.
        return reverse!(m1 .- master_params .+ t0)
    end
    return (master_params .- m0) .+ t0
end

# Ordinary (non-transfinite) grading: `F_Lc` integration through
# `mesh_curve`, the `minimumMeshSegments` floor and the recombination
# odd-count/`increaseN` adjustments; `filterPoints` runs inside the kernel
# with the BGM field and the `forceOdd` truncation rule.
function _model_curve_grade_params(m::GeoModel,curve::Integer,t0::Float64,
                                   t1::Float64,
                                   options::_ModelMesh1DOptions,
                                   caller::AbstractString)
    γ=u->_model_curve_point(m,curve,u,caller)
    field=_ModelCurveEntityField(m,options,String(caller))
    a,b=m.curves[curve]
    closed=a==b
    size_function=u->_model_curve_bgm(m,curve,u,γ(u),t0,t1,options,caller)
    minimum=_model_minimum_curve_segments(m,curve,options,caller)
    endpoints=((0,a),(0,b))
    # Recombination odd-N forcing (upstream counts N = elements + 1):
    # `(method != TRANSFINITE || flexible) && algoRecombine != 0`, then
    # `N%2==0 → N++` under RecombineAll or a recombined adjacent face, and
    # `increaseN` for blossom algorithms 2/4. filterPoints truncates its
    # removal count to even under the same condition.
    odd_nodes=m.meshing.recombine_algo!=0 &&
        (m.meshing.recombine_all ||
         any(face->haskey(m.meshing.recombine,(2,face)),
             _model_curve_faces(m,curve)))
    # `filterPoints` truncates its removal count to even only under
    # `recombineAll` — the face-recombine arm does not force it.
    force_odd=m.meshing.recombine_algo!=0 && m.meshing.recombine_all
    # `Integration` precision upstream is `lcIntegrationPrecision * lc`.
    precision=options.integration_precision*options.ctx_lc
    function grade(exact)
        _,graded=mesh_curve(
            γ,field;t0=t0,t1=t1,closed=closed,
            entity=(1,curve),endpoint_entities=endpoints,
            size_function=size_function,
            integration_precision=precision,
            minimum_segments=minimum,exact_edges=exact,
            force_odd=force_odd)
        return graded
    end
    parameters=grade(nothing)
    nedge=closed ? length(parameters) : length(parameters)-1
    # `N` counts nodes upstream (edges+1 for open curves; a closed curve's
    # shared vertex makes its node count equal its edge count). The odd-N
    # and `increaseN` adjustments only apply when recombination is actually
    # requested — `RecombineAll` or a recombined adjacent face (upstream
    # `meshGEdgeProcessing` lines ~696-715).
    if odd_nodes
        nnode=closed ? nedge : nedge+1
        iseven(nnode) && (nnode+=1)
        if m.meshing.recombine_algo in (2,4) && isodd(div(nnode+1,2)-1)
            nnode+=2
        end
        target=closed ? nnode : nnode-1
        target!=nedge && return grade(target)
    end
    return parameters
end

# Transfinite path — `F_Transfinite` positions on the stored law; flexible
# transfinite and the recombination odd-count rules are folded into
# `_transfinite_parameters` (which returns node positions on [0,1]).
function _model_curve_transfinite_params(m::GeoModel,curve::Integer,t0::Float64,
                                         t1::Float64,spec,
                                         caller::AbstractString)
    params=_transfinite_parameters(m,spec.num_nodes,spec.kind,spec.coef,
                                   caller,curve;reversed=spec.reversed)
    return t0 .+ params .* (t1-t0)
end

# `BGM_MeshSize` for a dim-1 discrete entity — the same `meshGEdgeProcessing`
# composition as `_model_curve_bgm`. `l1` interpolates the boundary-vertex
# prescribed sizes (discrete vertices carry none → `CTX::lc/10`, like
# upstream's unsized `discreteVertex`), `l2` is MAX_LC (discrete curvature is
# not evaluated upstream either — `discreteEdge::curvature` returns 0), `l4`
# is unset; `l5` still applies through `setSizeAtParametricPoints`.
function _model_discrete_curve_bgm(m::GeoModel,curve::Integer,u::Float64,xyz,
                                   t0::Float64,t1::Float64,boundary,
                                   options::_ModelMesh1DOptions,
                                   caller::AbstractString)
    l1=GMSH_MAX_SIZE
    if options.from_points && length(boundary)>=2
        function lc_of(key)
            key[1]==0 || return GMSH_MAX_SIZE
            return _model_vertex_lc(m,key[2])
        end
        lc1=lc_of(boundary[1]); lc2=lc_of(boundary[2])
        l1=(lc1>=GMSH_MAX_SIZE && lc2>=GMSH_MAX_SIZE) ? options.ctx_lc/10 :
            begin
                alpha=t1==t0 ? 0.0 : (u-t0)/(t1-t0)
                (1-alpha)*lc1+alpha*lc2
            end
    end
    l3=options.inner_field===nothing ? GMSH_MAX_SIZE :
        field_value(options.inner_field,xyz[1],xyz[2],xyz[3],(1,curve))
    l5=_model_curve_l5(m,curve,u)
    lc=min(l1,l3,GMSH_MAX_SIZE,l5)
    if options.callback!==nothing
        lc=_apply_size_callback(options.callback,1,curve,xyz[1],xyz[2],
                                xyz[3],lc,caller)
    end
    lc=min(options.ctx_lc,lc)
    lc=max(lc,options.lc_min)
    lc=min(lc,options.lc_max)
    if lc<=0
        _model_mesh_error(options,
            "Wrong mesh element size lc = $lc on curve $curve")
        lc=options.ctx_lc
    end
    return lc*options.lc_factor
end

# Discrete-curve re-meshing (`discreteEdge::mesh` when `_discretization` is
# non-empty): grade through the stored frame parametrization. Closed
# discrete loops (a single boundary vertex) mesh `closed` like upstream.
function _model_discrete_curve_grade_params(m::GeoModel,curve::Integer,
                                            record::DiscreteEntity,
                                            options::_ModelMesh1DOptions,
                                            caller::AbstractString)
    _,frames=_discrete_eval_setup(m,1,curve,caller)
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    γ=u->_discrete_curve_point(frames,u,caller)
    field=_ModelCurveEntityField(m,options,String(caller))
    boundary=record.boundary
    endpoints=length(boundary)>=2 ?
        (boundary[1],boundary[2]) : (nothing,nothing)
    closed=length(boundary)==1
    size_function=u->_model_discrete_curve_bgm(
        m,curve,u,γ(u),t0,t1,boundary,options,caller)
    _,parameters=mesh_curve(
        γ,field;t0=t0,t1=t1,closed=closed,
        entity=(1,curve),endpoint_entities=endpoints,
        size_function=size_function,
        integration_precision=options.integration_precision*options.ctx_lc,
        minimum_segments=_model_minimum_curve_segments(m,curve,options,caller))
    return parameters
end

# One curve's `meshGEdge`, in `meshGEdge::operator()` order: the
# `discreteEdge::mesh` early return for non-parametrized records, the skip
# gates (`degenerate(1)`, `Mesh.OnlyVisible`, `Mesh.OnlyEmpty` — all keep the
# stored mesh, like upstream returning before `deMeshGEdge`), then
# `MeshExtrudedCurve`, the periodic master copy, and
# `meshGEdgeProcessing` (`degenerate(0)`/zero-length `N=1`, transfinite,
# ordinary grading). Returns the stored parameter vector, `:keep` when the
# existing discretization is preserved, or `:pending` when a source/master
# has not been meshed yet (upstream leaves the curve PENDING and retries).
function _model_mesh_curve!(m::GeoModel,curve::Integer,
                            options::_ModelMesh1DOptions,
                            caller::AbstractString)
    record=get(m.discrete,(1,curve),nothing)
    if record!==nothing && !_model_discrete_curve_parametrized(record)
        return :keep   # `_discretization.empty()` upstream — elements kept
    end
    degenerate=_model_curve_degenerate_kind(m,curve)
    degenerate==1 && return :keep
    if options.mesh_only_visible
        model_entity_visibility(m,1,curve)==0 && return :keep
    end
    if options.mesh_only_empty
        params=get(m.curve_params,curve,nothing)
        params!==nothing && !isempty(params) && return :keep
        # Discrete/attached element tables are the entity's stored `lines`
        # upstream — a record that still carries elements stays untouched.
        record!==nothing && !isempty(record.element_types) && return :keep
    end
    # `deMeshGEdge` runs before meshing — dropping stale params here keeps a
    # failed grade from leaving an old discretization behind.
    delete!(m.curve_params,curve)
    # `MeshExtrudedCurve` runs before the periodic master copy upstream.
    extrude=get(m.meshing.extrude,(1,curve),nothing)
    if extrude!==nothing && !isempty(extrude.layers)
        return _model_curve_extrude_params(m,curve,extrude,options,caller)
    end
    # Periodic slave: copy the master's stored parameters once it is DONE.
    constraint=get(m.periodic,(1,curve),nothing)
    if constraint!==nothing
        master=abs(constraint.master_entity)
        haskey(m.curves,master) || haskey(m.discrete,(1,master)) ||
            throw(ArgumentError(
                "$caller: Curve[$curve] periodic master $master is missing"))
        master_params=get(m.curve_params,master,nothing)
        if master_params===nothing && haskey(m.discrete,(1,master))
            master_params=_model_discrete_curve_params(
                m,master,m.discrete[(1,master)],caller)
        end
        master_params===nothing && return :pending
        return _model_curve_periodic_params(
            m,curve,master,master_params,constraint.reversed,caller)
    end
    record!==nothing &&
        return _model_discrete_curve_grade_params(
            m,curve,record,options,caller)
    a,b=m.curves[curve]
    # Upstream warns "Skipping curve with no begin or end point" at the
    # vertex-creation stage and leaves the curve PENDING (retried, never
    # DONE) when a vertex is missing.
    if !(haskey(m.points,a) && haskey(m.points,b))
        _model_mesh_warn(options,
            "Skipping curve $curve with no begin or end point")
        return :pending
    end
    kind=get(m.curve_types,curve,:line)
    if kind===:line && a==b
        # A Line with coincident endpoints has no parametrization — upstream's
        # closed-`Line` shortcut evaluates `position(0.5)` to the vertex, so
        # `length=0` unconditionally and bounds stay the [0,1] line frame.
        t0,t1=0.0,1.0
        degenerate==0 && return Float64[t0,t1]
        if options.tolerance_edge_length==0.0
            return Float64[t0,t1]
        end
        minimum=_model_minimum_curve_segments(m,curve,options,caller)
        return collect(range(t0,t1,length=minimum+1))
    end
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    degenerate==0 && return _model_curve_single_edge_params(m,curve,caller)
    γ=u->_model_curve_point(m,curve,u,caller)
    length=curve_length(γ;t0=t0,t1=t1,
                        integration_precision=options.integration_precision *
                            options.ctx_lc)
    if length==0.0
        if options.tolerance_edge_length==0.0
            # `!length && !toleranceEdgeLength` → N=1 upstream.
            return _model_curve_single_edge_params(m,curve,caller)
        end
        # With a positive edge-length tolerance upstream still integrates the
        # zero-length metric: `a=0` → `N = minimumMeshSegments + 1` nodes at
        # the collapsed position.
        minimum=_model_minimum_curve_segments(m,curve,options,caller)
        return collect(range(t0,t1,length=minimum+1))
    end
    spec=get(m.meshing.transfinite_curves,curve,nothing)
    spec!==nothing &&
        return _model_curve_transfinite_params(m,curve,t0,t1,spec,caller)
    return _model_curve_grade_params(m,curve,t0,t1,options,caller)
end

# `Mesh1D` — every curve starts PENDING; periodic slaves and curve-extruded
# tops wait for their source, mirroring upstream's retry loop bounded by
# `Mesh.MaxRetries`. Curves still pending after the cap (periodic cycles)
# keep no stored mesh, matching upstream's silently-starved status.
function _model_mesh_curves!(m::GeoModel,options::_ModelMesh1DOptions,
                           caller::AbstractString)
    pending=Set{Int}(keys(m.curves))
    for (dim,tag) in keys(m.discrete)
        dim==1 && push!(pending,tag)
    end
    iter=0
    while !isempty(pending)
        for curve in sort!(collect(pending))
            result=_model_mesh_curve!(m,curve,options,caller)
            result===:pending && continue
            result===:keep || (m.curve_params[curve]=result)
            delete!(pending,curve)
        end
        isempty(pending) && break
        iter+=1
        iter>options.max_retries && break
    end
    return nothing
end

# The `old < 1` arm of `GModel::mesh(ask>1)`: re-run `Mesh1D` when the stored
# 1-D state is incomplete — any native curve without stored parameters that
# is not a `.geo Degenerated` keep (upstream's DONE-at-entry), or any
# non-parametrized discrete curve (upstream leaves it PENDING forever, so
# `meshDone1D` never holds).
function _model_mesh_needs_dim1(m::GeoModel)
    for curve in keys(m.curves)
        haskey(m.curve_params,curve) && continue
        curve in m.meshing.degenerated && continue
        return true
    end
    for ((dim,tag),record) in _discrete_mesh_records_model(m)
        dim==1 || continue
        haskey(m.curve_params,tag) && continue
        _model_discrete_curve_parametrized(record) || return true
    end
    return false
end

# ── Parts for `_geo_merge_entity_meshes` ─────────────────────────────────────
# One `(0, point, Mesh)` part per model point — upstream `Mesh0D` classifies
# a mesh vertex on every model vertex, positioned at the vertex's own
# coordinates (curve parts snap their endpoints to it, mirroring upstream's
# shared `MVertex`).
function _model_point_mesh_parts(m::GeoModel,caller::AbstractString)
    parts=Tuple{Int,Int,Mesh}[]
    for (tag,point) in m.points
        coords=Matrix{Float64}(undef,3,1)
        coords[1,1],coords[2,1],coords[3,1]=point
        push!(parts,(0,tag,Mesh(coords)))
    end
    # Discrete point entities keep their stored nodes (type-15 elements
    # upstream); the point-element MSH representation lands with the writer.
    for ((dim,tag),record) in _discrete_mesh_records_model(m)
        dim==0 || continue
        n=size(record.node_coords,2)
        n==0 && continue
        push!(parts,(0,tag,Mesh(Matrix{Float64}(record.node_coords[:,1:n]))))
    end
    return parts
end

# Evaluate a stored curve parameter for element emission. `:line` curves go
# through `_periodic_curve_point` — the same affine evaluator the surface
# PSLG uses — so curve-part nodes deduplicate bitwise against the boundary
# nodes surface meshing generates (the merge keys on `reinterpret`ed
# coordinates).
function _model_curve_part_point(m::GeoModel,curve::Integer,u::Float64,
                                 caller::AbstractString)
    _curve_type(m,curve)===:line &&
        return _periodic_curve_point(m,curve,u,caller)
    return _model_curve_point(m,curve,u,caller)
end

# Segment parts from the stored `curve_params` plus any discrete (1,*)
# records' line elements. Open-curve endpoints snap to the model points —
# upstream shares the vertex `MVertex`, so endpoint nodes must deduplicate
# bitwise against the `(0, point)` parts. Closed curves (`a==b`) connect
# their last node back to the first, like upstream's shared end vertex.
function _model_curve_mesh_parts(m::GeoModel,caller::AbstractString)
    parts=Tuple{Int,Int,Mesh}[]
    for (curve,params) in sort!(collect(m.curve_params))
        haskey(m.curves,curve) || continue
        a,b=m.curves[curve]
        closed=a==b
        us=Float64.(params)
        # Stored parameters can carry the duplicated closed endpoint
        # (transfinite and periodic copies span [t0,t1] inclusive); the
        # shared vertex is only stored once upstream.
        if closed && length(us)>1
            _,hi=_model_curve_param_bounds(m,curve,caller)
            us[end]>=hi && (us=us[1:end-1])
        end
        length(us)<(closed ? 1 : 2) && continue
        γ=u->_model_curve_part_point(m,curve,u,caller)
        coords=Matrix{Float64}(undef,3,length(us))
        for (i,u) in enumerate(us)
            p=γ(u)
            coords[1,i],coords[2,i],coords[3,i]=p
        end
        if closed
            # The shared vertex is the model vertex upstream — snap the
            # first node so it deduplicates bitwise with the (0, tag) part.
            haskey(m.points,a) &&
                (coords[:,1]=collect(m.points[a]))
        else
            haskey(m.points,a) &&
                (coords[:,1]=collect(m.points[a]))
            haskey(m.points,b) &&
                (coords[:,end]=collect(m.points[b]))
        end
        nseg=closed ? length(us) : length(us)-1
        nseg<1 && continue
        segs=Matrix{Int32}(undef,2,nseg)
        reversed=get(m.meshing.reverse,(1,curve),false)
        for i in 1:nseg
            j=closed ? mod1(i+1,length(us)) : i+1
            if reversed
                segs[1,i]=Int32(j); segs[2,i]=Int32(i)
            else
                segs[1,i]=Int32(i); segs[2,i]=Int32(j)
            end
        end
        tag=_model_projection_legacy_tag(
            _model_projection_physical_tags(m,1,curve))
        push!(parts,(1,curve,Mesh(coords;segs=segs,
                                 seg_tag=fill(Int32(tag),size(segs,2)))))
    end
    # Discrete (1,*) records: a parametrized record re-meshed this pass emits
    # its freshly graded discretization through the stored frames
    # (`discreteEdge::mesh` replaces the elements); a record without stored
    # parameters keeps its element table verbatim, like upstream's
    # `_discretization.empty()` early return. Native curves whose parts the
    # loop above already emitted are skipped here — their attached records
    # only surface elements when the curve kept no stored mesh
    # (`degenerate(1)`/`MeshOnlyEmpty` keeps upstream).
    for ((dim,tag),record) in _discrete_mesh_records_model(m)
        dim==1 || continue
        params=get(m.curve_params,tag,nothing)
        physical=_model_projection_legacy_tag(
            _model_projection_physical_tags(m,1,tag))
        if params!==nothing && !haskey(m.curves,tag)
            _,frames=_discrete_eval_setup(m,1,tag,caller)
            γ=u->_discrete_curve_point(frames,u,caller)
            us=Float64.(params)
            closed=length(record.boundary)==1
            if closed && length(us)>1
                _,hi=_model_curve_param_bounds(m,tag,caller)
                us[end]>=hi && (us=us[1:end-1])
            end
            length(us)<(closed ? 1 : 2) && continue
            coords=Matrix{Float64}(undef,3,length(us))
            for (i,u) in enumerate(us)
                p=γ(u)
                coords[1,i],coords[2,i],coords[3,i]=p
            end
            nseg=closed ? length(us) : length(us)-1
            nseg<1 && continue
            segs=Matrix{Int32}(undef,2,nseg)
            for i in 1:nseg
                j=closed ? mod1(i+1,length(us)) : i+1
                segs[1,i]=Int32(i); segs[2,i]=Int32(j)
            end
            push!(parts,(1,tag,Mesh(coords;segs=segs,
                                   seg_tag=fill(physical,nseg))))
            continue
        end
        haskey(m.curves,tag) && haskey(m.curve_params,tag) && continue
        isempty(record.element_types) && continue
        local_index=Dict{Int32,Int}(nt=>i for (i,nt) in
                                    enumerate(record.node_tags))
        segs=Matrix{Int32}(undef,2,0);seg_tag=Int32[]
        for (i,msh_type) in enumerate(record.element_types)
            msh_dimension(Int(msh_type))==1 || continue
            conn=record.element_nodes[i]
            length(conn)==2 || continue
            all(n->haskey(local_index,n),conn) || continue
            segs=hcat(segs,Int32[local_index[conn[1]],
                                 local_index[conn[2]]])
            push!(seg_tag,physical)
        end
        size(segs,2)>0 || continue
        push!(parts,(1,tag,Mesh(record.node_coords;segs=segs,
                               seg_tag=seg_tag)))
    end
    return parts
end

# ── `GModel::mesh` dim-1 entry ───────────────────────────────────────────────
# `Mesh0D` classifies a mesh vertex on every model vertex and `Mesh1D` fills
# `m.curve_params` through `meshGEdge`. Returns the `(dim, tag, Mesh)` parts
# the merge step consumes — upstream `Mesh 1` also deletes face/volume
# meshes, which the caller mirrors by replacing the stored merged mesh.
# `Mesh 0` is an upstream no-op and never reaches here.
function _model_mesh_dim01!(m::GeoModel,options::_ModelMesh1DOptions,
                          caller::AbstractString)
    _model_mesh_curves!(m,options,caller)
    parts=_model_point_mesh_parts(m,caller)
    append!(parts,_model_curve_mesh_parts(m,caller))
    return parts
end

# Dim-0/dim-1 parts without re-grading — `Mesh 2`/`Mesh 3` reuse the stored
# discretizations (`meshDone1D` upstream), re-running `Mesh1D` only when a
# curve's parameters are missing (the `old < 1` arm). Grading is split from
# part emission so `_geo_mesh_model` can interleave surface meshing between
# them: upstream `meshGFace` writes refined boundary subdivisions back onto
# each `GEdge`'s mesh, so the emitted line elements must come from the
# post-surface `curve_params`, not the pre-refinement ones.
function _model_mesh_dim01_grade!(m::GeoModel,
                                  options::_ModelMesh1DOptions,
                                  caller::AbstractString)
    _model_mesh_needs_dim1(m) && _model_mesh_curves!(m,options,caller)
    return nothing
end

function _model_mesh_dim01_parts!(m::GeoModel,
                                options::_ModelMesh1DOptions,
                                caller::AbstractString)
    _model_mesh_dim01_grade!(m,options,caller)
    parts=_model_point_mesh_parts(m,caller)
    append!(parts,_model_curve_mesh_parts(m,caller))
    return parts
end

# Edges incident to exactly one triangle are the domain boundary; their
# endpoints are the only nodes a boundary curve's chain can visit. Restricting
# `_curve_parameter_nodes` to them keeps an interior node that happens to lie
# collinear with a curve from corrupting the classification.
function _mesh_boundary_classification(mesh::Mesh)
    counts=Dict{NTuple{2,Int32},Int}()
    @inbounds for cell in axes(mesh.tris,2),(u,v) in ((1,2),(2,3),(3,1))
        a=mesh.tris[u,cell]; b=mesh.tris[v,cell]
        key=a<b ? (a,b) : (b,a)
        counts[key]=get(counts,key,0)+1
    end
    eligible=falses(nnodes(mesh))
    edges=Set{NTuple{2,Int32}}()
    for (key,n) in counts
        n==1 || continue
        push!(edges,key)
        eligible[key[1]]=true; eligible[key[2]]=true
    end
    return eligible,edges
end

# Upstream `meshGFace` stores the refined boundary subdivision back on each
# `GEdge`'s mesh — the emitted line elements then share the face's vertices
# exactly. Mirror that: recover the final (parameter,node) chain on every
# boundary curve, record it in `m.curve_params`, and snap the boundary nodes
# onto the analytic curve point — every sibling part evaluates the same
# `_periodic_curve_point` expression, so the merge deduplicates them
# bitwise. `Degenerated` (tooSmall) curves emit no elements upstream and are
# skipped, as are zero-length `a==b` lines, which have no boundary chain.
# Embedded curves classify over all triangle edges — their recovered
# constraint chain is interior — and get the same writeback.
function _model_surface_boundary_writeback!(m::GeoModel,t::Int,mesh::Mesh,
                                            caller::AbstractString)
    # A periodic slave surface is a verbatim mapped copy of its master — its
    # boundary nodes are owned by the mapping — and a periodic master's
    # boundary discretization is shared with the slave's copy. Upstream
    # keeps the edge meshes on both sides in lockstep, so the whole surface
    # skips the per-surface writeback.
    haskey(m.periodic,(2,t)) && return mesh
    for constraint in values(m.periodic)
        constraint.dim==2 && Int(constraint.master_entity)==t && return mesh
    end
    boundary_curves=Set{Int}()
    for loop in m.surfaces[t],signed in m.loops[loop]
        push!(boundary_curves,abs(signed))
    end
    embedded_curves=Int[]
    for (edim,etag) in get(m.embeds,(2,t),NTuple{2,Int}[])
        edim==1 && push!(embedded_curves,etag)
    end
    isempty(boundary_curves) && isempty(embedded_curves) && return mesh
    eligible_nodes,eligible_edges=_mesh_boundary_classification(mesh)
    all_edges=_model_projection_triangle_edges(mesh)
    # Nodes lying on a periodic-constrained curve — including vertices shared
    # with an unconstrained neighbour — keep the positions the periodic copy/
    # snap produced (`affine(master)`); snapping them to the analytic curve
    # point would break the bitwise slave↔master correspondence.
    protected_nodes=falses(nnodes(mesh))
    for curve in sort!(collect(boundary_curves))
        _model_curve_periodic_involved(m,curve) || continue
        haskey(m.curves,curve) || continue
        a,b=m.curves[curve]
        a==b && continue
        entries=_model_curve_chain_entries(
            m,mesh,curve,eligible_nodes,eligible_edges,caller)
        entries===nothing && continue
        for (_,node) in entries
            protected_nodes[node]=true
        end
    end
    for curve in sort!(collect(boundary_curves))
        _model_surface_curve_writeback!(
            m,curve,mesh,eligible_nodes,eligible_edges,
            protected_nodes,caller)
    end
    for curve in sort!(embedded_curves)
        _model_surface_curve_writeback!(
            m,curve,mesh,trues(nnodes(mesh)),all_edges,
            protected_nodes,caller)
    end
    return mesh
end

# A curve participates in a dim-1 periodic relation as slave or master —
# upstream keeps both sides' edge meshes in lockstep, so the stored
# parameters and node positions must not diverge per-surface.
function _model_curve_periodic_involved(m::GeoModel,curve::Int)
    haskey(m.periodic,(1,curve)) && return true
    for constraint in values(m.periodic)
        constraint.dim==1 && Int(constraint.master_entity)==curve &&
            return true
    end
    return false
end

# Soft variant of `_curve_parameter_nodes`: returns `nothing` instead of
# throwing when the curve has no usable mesh-edge chain.
function _model_curve_chain_entries(m::GeoModel,mesh::Mesh,curve::Int,
                                    eligible_nodes,eligible_edges,
                                    caller::AbstractString)
    try
        entries,_=_curve_parameter_nodes(
            m,mesh,curve,eligible_nodes,eligible_edges,0.0,caller)
        return entries
    catch
        return nothing
    end
end

function _model_surface_curve_writeback!(m::GeoModel,curve::Int,mesh::Mesh,
                                         eligible_nodes,eligible_edges,
                                         protected_nodes,
                                         caller::AbstractString)
    curve in m.meshing.degenerated && return nothing
    haskey(m.curves,curve) || return nothing
    a,b=m.curves[curve]
    a==b && return nothing
    _model_curve_periodic_involved(m,curve) && return nothing
    entries,_=_curve_parameter_nodes(
        m,mesh,curve,eligible_nodes,eligible_edges,0.0,caller)
    m.curve_params[curve]=Float64[entry[1] for entry in entries]
    @inbounds for (parameter,node) in entries
        protected_nodes[node] && continue
        point=_periodic_curve_point(m,curve,parameter,caller)
        mesh.coords[1,node]=point[1]
        mesh.coords[2,node]=point[2]
        mesh.coords[3,node]=point[3]
    end
    return nothing
end
