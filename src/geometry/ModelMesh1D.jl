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
        n=hypot(ax[1],ax[2],ax[3])
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
        c=cone.center; ax=cone.axis; n=hypot(ax[1],ax[2],ax[3])
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

# `_ModelMesh1DOptions` for the model-API path (`mesh_model_surface` and
# friends): upstream `CTX` defaults where the `.geo` context would supply
# option values — `MeshSizeMin` 0, `MeshSizeMax` 1e22, `MeshSizeFromPoints` 1,
# `MeshSizeFromCurvature` 0, `LcIntegrationPrecision` 1e-9, `MinLineNodes` 2,
# `MinCircleNodes` 7, `MinCurveNodes` 3, `ToleranceEdgeLength` 0,
# `Geometry.Tolerance` 1e-8, `MaxRetries` 10. There is no model-level field
# store (`inner_field`) and no Msg channel — warning/error sinks stay
# `nothing`, so upstream's print-and-continue paths continue silently.
function _model_default_mesh1d_options(m::GeoModel,caller::AbstractString)
    geometry_tolerance=1e-8
    return _ModelMesh1DOptions(
        _model_mesh_lc(m,geometry_tolerance,caller),
        0.0,1e22,m.meshing.lc_factor,
        true,0.0,1e-9,
        2,7,3,0.0,geometry_tolerance,
        false,false,10,
        nothing,m.meshing.size_callback,
        nothing,nothing)
end

# Curves a `_model_mesh_curve!` `:pending` result may be waiting on — the
# periodic master and the dim-1 extrusion source (a point-generatrix source
# has no curve dependency).
function _model_curve_mesh_dependencies(m::GeoModel,curve::Integer)
    deps=Int[]
    constraint=get(m.periodic,(1,Int(curve)),nothing)
    constraint!==nothing && push!(deps,Int(abs(constraint.master_entity)))
    src=get(m.meshing.extrude_sources,(1,Int(curve)),nothing)
    src!==nothing && src[1]==1 && push!(deps,abs(src[2]))
    return deps
end

# Lazily run `meshGEdge` for boundary curves a `mesh_model_surface` call
# still lacks stored discretizations for — upstream `generate(2)` runs the
# whole `Mesh1D` pass first, so `meshGFace` always finds
# `GEdge::mesh_vertices` populated. Periodic masters and extrusion sources
# join the pending set (upstream's global pass reaches them too); curves
# that stay PENDING past `Mesh.MaxRetries` keep no stored parameters, like
# upstream's silently-starved status.
function _model_surface_mesh_curves!(m::GeoModel,curves,caller::AbstractString)
    pending=Set{Int}()
    for curve in curves
        haskey(m.curve_params,Int(curve)) || push!(pending,Int(curve))
    end
    isempty(pending) && return nothing
    options=_model_default_mesh1d_options(m,caller)
    for _ in 1:options.max_retries
        isempty(pending) && break
        progressed=false
        for curve in sort!(collect(pending))
            result=_model_mesh_curve!(m,curve,options,caller)
            if result===:pending
                for dep in _model_curve_mesh_dependencies(m,curve)
                    haskey(m.curve_params,dep) && continue
                    dep in pending && continue
                    push!(pending,dep); progressed=true
                end
                continue
            end
            result===:keep || (m.curve_params[curve]=result)
            delete!(pending,curve); progressed=true
        end
        progressed || break
    end
    return nothing
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
# extrapolates past the last like upstream. The API stores normalized [0,1]
# positions, so the native `u` is mapped into the same frame first.
function _model_curve_l5(m::GeoModel,curve::Integer,u::Float64,
                        caller::AbstractString)
    entries=get(m.meshing.size_at_params,(1,curve),nothing)
    (entries===nothing || isempty(entries)) && return GMSH_MAX_SIZE
    pairs=sort!(Tuple{Float64,Float64}[(u,lc) for (plist,lc) in entries
                                       for u in plist]; by=first)
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    v=(u-t0)/(t1-t0)
    us=Float64[p[1] for p in pairs]
    n=length(us)
    it=searchsortedfirst(us,v)          # first index with us[it] >= v (n+1 past end)
    i1=min(it-1,n-1)                    # 0-based lower_bound position, capped
    i0=max(1,i1)-1
    u0,l0=pairs[i0+1]; u1,l1=pairs[i1+1]
    (i1==i0 || u0==u1) && return l0
    alpha=(v-u0)/(u1-u0)
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
    l5=_model_curve_l5(m,curve,u,caller)
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
    # A closed native curve needs at least three segments upstream: its
    # repeated endpoint is omitted from the `N = minimumMeshSegments + 1`
    # target, leaving N−1 unique nodes and edges. Gmsh emits a three-edge loop on a coarse
    # closed spline (still four under `MinCurvePoints 5`, matching the
    # regular floor) — a two-segment digon only arises on degenerate input.
    if !is_occ && kind!==:line && haskey(m.curves,curve)
        a,b=m.curves[curve]
        a==b && (np=max(np,3))
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
# mass-gated odd-count adjustment; `filterPoints` runs inside the kernel
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
    # `N%2==0 && a>0.75 → N++` under RecombineAll or a recombined adjacent
    # face. Gmsh 4.15.2's `increaseN` is an identity for all algorithms.
    odd_nodes=m.meshing.recombine_algo!=0 &&
        (m.meshing.recombine_all ||
         any(face->haskey(m.meshing.recombine,(2,face)),
             _model_curve_faces(m,curve)))
    # `filterPoints` truncates its removal count to even only under
    # `recombineAll` — the face-recombine arm does not force it.
    force_odd=m.meshing.recombine_algo!=0 && m.meshing.recombine_all
    native_line=_curve_type(m,curve)===:line &&
        get(m.curve_geometry,curve,nothing)===nothing
    # Reuse the existing size samples/smoothing for the count-only native
    # Line stencil; a different differentiation residual can straddle 0.75.
    # Do not invoke field callbacks or rebuild the placement primitive.
    count_primitive=odd_nodes && native_line ?
        points->_model_native_line_recombination_mass(
            m.points[a],m.points[b],points) : nothing
    # `Integration` precision upstream is `lcIntegrationPrecision * lc`.
    precision=options.integration_precision*options.ctx_lc
    _,parameters=mesh_curve(
        γ,field;t0=t0,t1=t1,closed=closed,
        entity=(1,curve),endpoint_entities=endpoints,
        size_function=size_function,
        integration_precision=precision,
        minimum_segments=minimum,force_odd=force_odd,
        _recombine_odd_nodes=odd_nodes,
        _recombination_primitive=count_primitive)
    return parameters
end

# ── Arc-length inversion ─────────────────────────────────────────────────────
#
# Upstream reparametrizes every built-in curve by geometric length: the
# transfinite law positions are length fractions, and `GEdge::point(par)`
# inverts them into the native parameter through the edge's discretization.
# `Line` and `Circle` parametrizations already advance at constant speed, so
# their fractions map linearly; the other kinds (`Ellipse`, `Spline`,
# `BSpline`, `Bezier`, `Nurbs`, OCC records) invert through an adaptively
# refined cumulative-length table — the same `_adaptive_points`/trapezoid
# primitive `mesh_curve` uses for its node placement.
function _model_curve_length_table(m::GeoModel,curve::Int,
                                   caller::AbstractString)
    kind=_curve_type(m,curve)
    # Uniform-speed parametrizations need no table.
    (kind===:line || kind===:circle) && return nothing
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    γ=function(u)
        point=_model_curve_point(m,curve,u,caller)
        all(isfinite,point) || throw(ArgumentError(
            "$caller: Curve[$curve] evaluation is not Float64-representable"))
        return point
    end
    evaluate(t,_,_)=_length_point(γ,nothing,t,t0,t1,caller)
    points=_adaptive_points(evaluate,t0,t1,_GMSH_INTEGRATION_PRECISION,
                            _DEFAULT_MAX_INTEGRATION_POINTS,
                            _GMSH_MIN_INTEGRATION_DEPTH,
                            _GMSH_MAX_INTEGRATION_DEPTH,caller)
    points[end].p>0 || throw(ArgumentError(
        "$caller: Curve[$curve] has zero geometric length"))
    return points
end

# Map the transfinite length fraction `t∈[0,1]` to the stored/native
# parameter. `table` comes from `_model_curve_length_table`; uniform-speed
# kinds take the linear map (bitwise identical to the legacy behavior).
function _model_curve_length_param(m::GeoModel,curve::Int,t::Float64,
                                   caller::AbstractString,
                                   table=_model_curve_length_table(
                                       m,curve,caller))
    if table===nothing
        t0,t1=_model_curve_param_bounds(m,curve,caller)
        return t0+t*(t1-t0)
    end
    u=_invert_primitive(table,t*table[end].p)
    isfinite(u) || throw(ArgumentError(
        "$caller: Curve[$curve] arc-length inversion is not Float64-representable"))
    return u
end

# ── Transfinite density (`F_Transfinite` + `Integration`) ────────────────────
#
# Gmsh does NOT reparametrize curved transfinite curves by arc length: the
# stored law defines a cell-size DENSITY in the normalized parameter —
# `F_Transfinite` evaluates `val(u) = ‖C′(u)‖/cellsize((u−t0)/(t1−t0))`, and
# `Integration`/`RecursiveIntegration` builds its primitive with the same
# adaptive trapezoid `F_Lc` uses. Nodes then sit at equal primitive marks,
# linearly inverted on (t, primitive) (`meshGEdge` NUMP walk). For
# native Lines the sampled density and primitive inversion can differ from
# the closed-form grading law, even with constant geometric speed. Circles
# without FlexibleTransfinite retain their established analytic positions;
# flexible Circles, nonuniform native Lines and
# non-uniform-speed curves dispatch through the same three `val` arms
# upstream takes:
#
#   * default arm — `coef <= 0`, `coef == 1`, or a beta coefficient < 1:
#     `val ∝ ‖C′‖` → the primitive is arc length → uniform LENGTH fractions,
#     recovered by the `_model_curve_length_param` inversion;
#   * unknown-type arm — types 8/9 (`Beta_Symmetrical*`), reversed HWall
#     records (the `type == 5..7` transforms match positive types only), and
#     the ±0 wildcard (`:uniform` → type 0): `val = 1` → uniform PARAMETER
#     fractions;
#   * law arm — progression/bump/beta, plus HWall coefficients solved through
#     the upstream bounded searches below: integrate `‖C′‖/cellsize` and
#     invert equal mass marks.
#
# `coeffTransfinite <= 0` can only reach the record through a verbatim `.geo
# Using` coefficient; upstream folds it into the default arm where the
# resulting mass is nonpositive — the existing closed-form paths already fall
# back to uniform there, and the curved path does the same through the
# `:length` arm.

# Stored kind → upstream `typeTransfinite`. `Power` aliases `Progression` at
# parse time; types 8/9 and the ±0 wildcard have no `val` case.
function _transfinite_law_type(kind::Symbol)
    kind===:progression && return 1
    kind===:bump && return 2
    kind===:beta && return 3
    kind===:progression_hwall && return 5
    kind===:bump_hwall && return 6
    kind===:beta_hwall && return 7
    kind===:beta_symmetrical && return 8
    kind===:beta_symmetrical_hwall && return 9
    return 0
end

# `dfbeta` (meshGEdge.cpp) — the beta-law density.
function _transfinite_dfbeta(t::Float64,beta::Float64)
    zlog=log((1.0+beta)/(beta-1.0))
    return 2.0beta/((1.0+beta-t)*(-1.0+beta+t)*zlog)
end

# `f_prog` + `newton_get_r` (a bounded bissection despite the name): solve
# `hw·(r^s − 1)/(r − 1) == length` for the progression ratio on [1,4] with
# `s = n − 1` segments; out-of-bracket inputs return 1.0 like upstream.
function _transfinite_hwall_ratio(hw::Float64,length::Float64,n::Int)
    segments=n-1
    function f(r)
        r==1.0 && return segments*hw/length
        return hw*(r^segments-1.0)/(r-1.0)/length
    end
    (f(1.0)<1.0 && f(4.0)>1.0) || return 1.0
    r1,r2=1.0,4.0
    while true
        r3=0.5*(r1+r2)
        f3=f(r3)
        abs(f3-1.0)<1e-12 && return r3
        f3<1.0 ? (r1=r3) : (r2=r3)
    end
end

# `f_bump` + `bissection_get_a`: solve the bump coefficient on [1e-8,100] so
# the first integrated unit lands at `hw/length`; out-of-bracket → 1.0.
function _transfinite_hwall_bump(hw::Float64,length::Float64,n::Int)
    t=hw/length
    function f(coef)
        if coef>1.0
            a=atan(1.0,sqrt(coef-1.0))/sqrt(coef-1.0)/n
            A=coef-1.0
            return (atan(sqrt(A))-atan(sqrt(A)*(1.0-2.0t)))/2.0/sqrt(A)/a
        end
        a=atanh(sqrt(1.0-coef))/sqrt(1.0-coef)/n
        A=1.0-coef
        return (atanh(sqrt(A))-atanh(sqrt(A)*(1.0-2.0t)))/2.0/sqrt(A)/a
    end
    (f(1e-8)>1.0 && f(100.0)<1.0) || return 1.0
    alpha1,alpha2=1e-8,100.0
    while true
        alpha3=0.5*(alpha1+alpha2)
        f3=f(alpha3)
        abs(f3-1.0)<1e-8 && return alpha3
        f3>1.0 ? (alpha1=alpha3) : (alpha2=alpha3)
    end
end

# `f_beta` + `bissection_get_beta`: solve beta on [1+1e-8,5] so the first
# integrated unit lands at `hw/length`; out-of-bracket → 100.0.
function _transfinite_hwall_beta(hw::Float64,length::Float64,n::Int)
    t=hw/length
    f(beta)=n*(1.0+atanh((t-1.0)/beta)/atanh(1.0/beta))
    (f(1.0+1e-8)>1.0 && f(5.0)<1.0) || return 100.0
    beta1,beta2=1.0+1e-8,5.0
    while true
        beta3=0.5*(beta1+beta2)
        f3=f(beta3)
        abs(f3-1.0)<1e-8 && return beta3
        f3>1.0 ? (beta1=beta3) : (beta2=beta3)
    end
end

# The `nbpt` `F_Transfinite` builds its cell partition from — the stored node
# count with the `FlexibleTransfinite` lcFactor truncation applied, WITHOUT
# the recombination odd bump (upstream applies that to N only, after the
# primitive is integrated).
function _transfinite_law_nodes(m::GeoModel,num_nodes::Int,
                                caller::AbstractString)
    m.meshing.flexible_transfinite || return num_nodes
    factor=m.meshing.lc_factor
    factor==0.0 && return num_nodes
    q=num_nodes/factor
    adjusted=isfinite(q) ?
        trunc(Int,clamp(q,-2.2e9,2.2e9)) :
        (q>0 ? Int64(2.2e9) : 0)
    return max(2,adjusted)
end

# The `(type, coef)` the `F_Transfinite` transforms leave behind: positive
# HWall types solve their wall height into an ordinary coefficient (type 5
# reciprocates a negative wall through the ratio; type 7 mirrors a negative
# wall through the signed type). Reversed HWall records (type −5..−7) do not
# transform and fall to the unknown-type arm — the same place grammar-only
# types 8/9 and the type-0 wildcard land.
function _transfinite_effective_law(spec,length::Float64,nbpt::Int)
    type=_transfinite_law_type(spec.kind)
    spec.reversed && (type=-type)
    coef=spec.coef
    if type==5
        coef=_transfinite_hwall_ratio(abs(coef),length,nbpt)
        spec.coef>0.0 || (coef=1.0/coef)
        type=1
    elseif type==6
        coef=_transfinite_hwall_bump(abs(coef),length,nbpt)
        type=2
    elseif type==7
        coef=_transfinite_hwall_beta(abs(coef),length,nbpt)
        type=3*(spec.coef>0.0 ? 1 : -1)
    end
    return type,coef
end

# Which `val` arm the effective law lands in.
function _transfinite_val_arm(type::Int,coef::Float64)
    (coef<=0.0 || coef==1.0 || (abs(type)==3 && coef<1.0)) &&
        return :length
    abs(type) in (1,2,3) && return :density
    return :param
end

# `F_Transfinite::operator()` — `val(t_normalized, speed)`. The progression
# cell index `i` partitions the normalized parameter scaled by the geometric
# length (upstream's `t * length / a`), and a reversed type reciprocates the
# ratio / mirrors the beta argument.
function _transfinite_val(type::Int,coef::Float64,length::Float64,nbpt::Int)
    atype=abs(type)
    if coef<=0.0 || coef==1.0 || (atype==3 && coef<1.0)
        return (t,d)->d*coef/length
    elseif atype==1
        r=type>=0 ? coef : 1.0/coef
        a=length*(r-1.0)/(r^(nbpt-1.0)-1.0)
        lgr=log(r)
        return function (t,d)
            argument=t*length/a*(r-1.0)+1.0
            # An extreme decreasing ratio can round the terminal logarithm
            # argument to zero. Native F_Transfinite contributes zero density
            # at that endpoint; make this explicit without converting Inf to
            # an integer. An unresolved interior partition is an error.
            if argument<=0.0 && t==1.0 && r<1.0
                return 0.0
            end
            (isfinite(argument) && argument>0.0) || throw(ArgumentError(
                "transfinite progression has an unrepresentable interior density partition"))
            index=log(argument)/lgr
            (isfinite(index) && typemin(Int)<=index<typemax(Int)) ||
                throw(ArgumentError(
                    "transfinite progression has an unrepresentable density partition index"))
            i=floor(Int,index)
            return d/(a*r^i)
        end
    elseif atype==2
        aa=if coef>1.0
            -4.0*sqrt(coef-1.0)*atan(1.0,sqrt(coef-1.0))/(nbpt*length)
        else
            2.0*sqrt(1.0-coef)*
                log(abs((1.0+1.0/sqrt(1.0-coef))/
                        (1.0-1.0/sqrt(1.0-coef))))/(nbpt*length)
        end
        b=-aa*length*length/(4.0*(coef-1.0))
        return (t,d)->d/(-aa*(t*length-length*0.5)^2+b)
    elseif atype==3
        return type<0 ? ((t,d)->_transfinite_dfbeta(1.0-t,coef)) :
                        ((t,d)->_transfinite_dfbeta(t,coef))
    else
        return (t,d)->1.0
    end
end

# Native GEO Line derivatives in `InterpolateCurve` use a bounded 1e-5
# difference, including the reciprocal multiplication shown here. The
# recombination threshold reads the numerically integrated density, so a
# coefficient near 0.75 cannot be compared directly to that threshold.
# This count-only evaluator leaves the established node-placement stencil
# and exact public Line value/derivative contracts unchanged.
function _transfinite_native_line_speed(first::NTuple{3,Float64},
                                         last::NTuple{3,Float64},t::Float64)
    left_step=t<1e-5 ? 0.0 : 1e-5
    right_step=t>1.0-1e-5 ? 0.0 : 1e-5
    reciprocal=1.0/(left_step+right_step)
    squared=0.0
    @inbounds for axis in 1:3
        delta=last[axis]-first[axis]
        left=first[axis]+(t-left_step)*delta
        right=first[axis]+(t+right_step)*delta
        derivative=(right-left)*reciprocal
        squared+=derivative*derivative
    end
    return sqrt(squared)
end

function _model_native_line_recombination_mass(first::NTuple{3,Float64},
                                               last::NTuple{3,Float64},points)
    previous=points[1]
    previous_density=_transfinite_native_line_speed(first,last,previous.t)/previous.h
    mass=0.0
    @inbounds for index in 2:length(points)
        current=points[index]
        density=_transfinite_native_line_speed(first,last,current.t)/current.h
        mass+=0.5*(previous_density+density)*(current.t-previous.t)
        previous=current
        previous_density=density
    end
    return mass
end

function _model_curve_transfinite_mass(m::GeoModel,curve::Int,spec,
                                       caller::AbstractString)
    # A caller with no stored law has the ordinary uniform transfinite
    # density. Nonpositive coefficients never satisfy the positive mass
    # threshold; the uniform coefficient's mass is safely above it.
    spec===nothing && return 1.0
    raw_type=_transfinite_law_type(spec.kind)*(spec.reversed ? -1 : 1)
    hwall=raw_type in (5,6,7)
    # Positive HWall types interpret a negative coefficient as the wall
    # side and transform it into a positive ordinary law first.
    (!hwall && spec.coef<=0.0) && return 0.0
    (!hwall && spec.coef==1.0) && return 1.0
    t0,t1=_model_curve_param_bounds(m,curve,caller)
    _transfinite_val_arm(raw_type,spec.coef)===:param &&
        !hwall && return t1-t0

    native_line=_curve_type(m,curve)===:line &&
        get(m.curve_geometry,curve,nothing)===nothing
    endpoints=native_line ?
        (m.points[first(m.curves[curve])],m.points[last(m.curves[curve])]) :
        ((0.0,0.0,0.0),(0.0,0.0,0.0))
    gamma=u->_model_curve_point(m,curve,u,caller)
    speed=t->native_line ?
        _transfinite_native_line_speed(endpoints[1],endpoints[2],t) :
        _length_point(gamma,nothing,t,t0,t1,caller).xp
    length_evaluate=(t,_,_)->begin
        d=speed(t)
        _IntegrationPoint(t,d,0.0,d,1.0)
    end
    length_points=_adaptive_points(length_evaluate,t0,t1,
        _GMSH_INTEGRATION_PRECISION,_DEFAULT_MAX_INTEGRATION_POINTS,
        _GMSH_MIN_INTEGRATION_DEPTH,_GMSH_MAX_INTEGRATION_DEPTH,caller)
    geometric_length=length_points[end].p
    geometric_length==0.0 && return 0.0
    type,coef=_transfinite_effective_law(spec,geometric_length,spec.num_nodes)
    nbpt=_transfinite_law_nodes(m,spec.num_nodes,caller)
    val=_transfinite_val(type,coef,geometric_length,nbpt)
    density_evaluate=(t,_,_)->begin
        d=speed(t)
        density=val((t-t0)/(t1-t0),d)
        _IntegrationPoint(t,density,0.0,d,1.0)
    end
    density_points=_adaptive_points(density_evaluate,t0,t1,
        _GMSH_INTEGRATION_PRECISION,_DEFAULT_MAX_INTEGRATION_POINTS,
        _GMSH_MIN_INTEGRATION_DEPTH,_GMSH_MAX_INTEGRATION_DEPTH,caller)
    return density_points[end].p
end

# `Integration` + the `meshGEdge` node-mark walk: integrate the transfinite
# density over [t0,t1] with the shared adaptive trapezoid, then place `count`
# nodes at equal primitive marks (linear inversion on (t,p)). The speed `d`
# is the same bounded finite difference the length table uses — upstream's
# `ge->firstDer` norm.
function _model_curve_transfinite_density_params(m::GeoModel,curve::Int,
                                                 t0::Float64,t1::Float64,
                                                 val,count::Int,
                                                 caller::AbstractString)
    γ=function(u)
        point=_model_curve_point(m,curve,u,caller)
        all(isfinite,point) || throw(ArgumentError(
            "$caller: Curve[$curve] evaluation is not Float64-representable"))
        return point
    end
    evaluate=function (t,_,_)
        lp=_length_point(γ,nothing,t,t0,t1,caller)
        tn=(t-t0)/(t1-t0)
        return _IntegrationPoint(t,val(tn,lp.xp),0.0,lp.xp,1.0)
    end
    points=_adaptive_points(evaluate,t0,t1,_GMSH_INTEGRATION_PRECISION,
                            _DEFAULT_MAX_INTEGRATION_POINTS,
                            _GMSH_MIN_INTEGRATION_DEPTH,
                            _GMSH_MAX_INTEGRATION_DEPTH,caller)
    a=points[end].p
    (isfinite(a) && a>0.0) || throw(ArgumentError(
        "$caller: Curve[$curve] transfinite law integrates to a nonpositive " *
        "mass"))
    params=Vector{Float64}(undef,count)
    params[1]=t0
    params[end]=t1
    b=a/(count-1)
    @inbounds for k in 1:count-2
        u=_invert_primitive(points,k*b)
        isfinite(u) || throw(ArgumentError(
            "$caller: Curve[$curve] transfinite inversion is not " *
            "Float64-representable"))
        params[k+1]=u
    end
    return params
end

# The native-parameter node list for a transfinite curve — the single source
# both `curve_params` and the transfinite surface/volume side chains consume,
# so every part emits bitwise-identical boundary nodes. Native straight
# nonuniform laws use F_Transfinite's sampled density primitive as upstream
# does; a constant geometric speed does not make its stepwise progression
# density or numerical inversion equal to the analytic grading formula.
# The standalone analytic grading helper remains a separate contract.
function _model_curve_transfinite_native_params(m::GeoModel,curve::Int,
                                                t0::Float64,t1::Float64,spec,
                                                caller::AbstractString)
    kind=_curve_type(m,curve)
    if kind===:line || (kind===:circle && m.meshing.flexible_transfinite)
        count=_flexible_transfinite_nodes(m,spec.num_nodes,curve,caller)
        count==2 && return Float64[t0,t1]
        raw_type=_transfinite_law_type(spec.kind)*(spec.reversed ? -1 : 1)
        if _transfinite_val_arm(raw_type,spec.coef)===:density || raw_type in (5,6,7)
            a,b=m.curves[curve]
            length=kind===:line ?
                _model_point_distance(m.points[a],m.points[b]) :
                _model_curve_length(m,curve,caller)
            (isfinite(length) && length>0.0) || throw(ArgumentError(
                "$caller: Curve[$curve] has nonpositive or nonfinite geometric length"))
            type,coef=_transfinite_effective_law(spec,length,spec.num_nodes)
            if _transfinite_val_arm(type,coef)===:density
                n_law=_transfinite_law_nodes(m,spec.num_nodes,caller)
                val=_transfinite_val(type,coef,length,n_law)
                return _model_curve_transfinite_density_params(
                    m,curve,t0,t1,val,count,caller)
            elseif _transfinite_val_arm(type,coef)===:length
                # HWall's original declared count can put its transform
                # outside the solver bracket (ordinary coefficient 1).
                # Do not solve it again with the recombination output count.
                return Float64[_convex_coordinate(t0,t1,p)
                    for p in _transfinite_uniform_parameters(count)]
            end
        end
    end
    if kind===:line || kind===:circle
        params=_transfinite_parameters(m,spec.num_nodes,spec.kind,
                                       spec.coef,caller,curve;
                                       reversed=spec.reversed)
        return Float64[t0+p*(t1-t0) for p in params]
    end
    # One shared integration: the converged primitive IS `ge->length()`, and
    # the :length arm inverts the same table — a second `curve_length` call
    # would re-integrate the identical evaluator bitwise.
    table=_model_curve_length_table(m,curve,caller)
    length=table[end].p
    # The HWall transforms run on `nbPointsTransfinite` BEFORE `nbpt` absorbs
    # the `lcFactor` division — the divided count only enters the `val`
    # partition (`F_Transfinite` line order).
    type,coef=_transfinite_effective_law(spec,length,spec.num_nodes)
    n_law=_transfinite_law_nodes(m,spec.num_nodes,caller)
    arm=_transfinite_val_arm(type,coef)
    count=_flexible_transfinite_nodes(m,spec.num_nodes,curve,caller)
    if arm===:length
        # `val ∝ speed` → uniform geometric-length fractions.
        return Float64[_model_curve_length_param(m,curve,f,caller,table)
                       for f in _transfinite_uniform_parameters(count)]
    elseif arm===:param
        # `val = 1` → the primitive is the parameter itself.
        return Float64[_convex_coordinate(t0,t1,k/(count-1))
                       for k in 0:count-1]
    end
    val=_transfinite_val(type,coef,length,n_law)
    return _model_curve_transfinite_density_params(m,curve,t0,t1,val,
                                                   count,caller)
end

# Transfinite path — `F_Transfinite` positions on the stored law; flexible
# transfinite and the recombination odd-count rules are folded into the count
# (`_flexible_transfinite_nodes`). The node list is the shared
# `_model_curve_transfinite_native_params` so stored parameters and boundary
# side chains agree bitwise.
function _model_curve_transfinite_params(m::GeoModel,curve::Integer,t0::Float64,
                                         t1::Float64,spec,
                                         caller::AbstractString)
    return _model_curve_transfinite_native_params(m,curve,t0,t1,spec,caller)
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
    l5=_model_curve_l5(m,curve,u,caller)
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
    spec=get(m.meshing.transfinite_curves,curve,nothing)
    length=if spec!==nothing && kind===:line
        # The F_One prepass integrates a constant derivative on a Line. Only
        # its zero-length gate is consumed by this transfinite branch; the
        # prescribed law computes its own line length where needed. Avoid
        # rebuilding exact line evaluators for every integration sample.
        chord=_model_point_distance(m.points[a],m.points[b])
        isfinite(chord) || throw(ArgumentError(
            "$caller: Curve[$curve] has non-finite geometric length"))
        chord
    else
        γ=u->_model_curve_point(m,curve,u,caller)
        curve_length(γ;t0=t0,t1=t1,
                     integration_precision=options.integration_precision *
                         options.ctx_lc)
    end
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
    _model_require_point_mesh_geometry(m,caller)
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
function _model_require_point_mesh_geometry(m::GeoModel,caller::AbstractString)
    displaced=Set{Int}()
    for ((dim,tag),record) in m.meshing.attached
        dim==0 && haskey(m.points,tag) && !isempty(record.node_tags) || continue
        point=ntuple(axis->record.node_coords[axis,1],3)
        point==m.points[tag] || push!(displaced,tag)
    end
    isempty(displaced) && return nothing
    for curve in sort!(collect(keys(m.curves)))
        for point in m.curves[curve]
            point in displaced || continue
            throw(ArgumentError("$caller: native Point[$point] has an attached " *
                "mesh vertex off its geometry and is incident to Curve[$curve]; " *
                "boundary meshing from displaced native vertices is not implemented"))
        end
    end
    for ((dim,tag),record) in _discrete_mesh_records_model(m)
        dim==1 || continue
        for (boundary_dim,point) in record.boundary
            boundary_dim==0 && point in displaced || continue
            throw(ArgumentError("$caller: native Point[$point] has an attached " *
                "mesh vertex off its geometry and is incident to Curve[$tag]; " *
                "boundary meshing from displaced native vertices is not implemented"))
        end
    end
    for (entity,embedded) in m.embeds, (dim,point) in embedded
        dim==0 && point in displaced || continue
        throw(ArgumentError("$caller: native Point[$point] has an attached " *
            "mesh vertex off its geometry and is embedded in entity $entity; " *
            "embedding displaced native vertices is not implemented"))
    end
    return nothing
end

# The current entity-part merger identifies vertices by entity and exact
# coordinates. Raw records can contain distinct tags at the same position,
# even for disconnected components of one curve. Reject that unresolved
# identity case before emitting a part rather than silently welding its nodes.
function _model_assert_raw_node_identity(record::DiscreteEntity,dim::Int,
                                          tag::Int,caller::AbstractString)
    seen=Dict{NTuple{3,UInt64},Int32}()
    for (node,node_tag) in enumerate(record.node_tags)
        key=ntuple(3) do axis
            coordinate=record.node_coords[axis,node]
            reinterpret(UInt64,coordinate==0.0 ? 0.0 : coordinate)
        end
        previous=get(seen,key,Int32(0))
        if previous!=0 && previous!=node_tag
            entity=dim==0 ? "Point" : "Curve"
            throw(ArgumentError("$caller: raw $entity[$tag] has distinct node " *
                "tags $previous and $node_tag at identical coordinates; " *
                "preserving coincident nodes within one entity is not implemented"))
        end
        seen[key]=node_tag
    end
    return nothing
end

# One `(0, point, Mesh)` part per model point — upstream `Mesh0D` classifies
# a mesh vertex on every model vertex, positioned at the vertex's own
# coordinates (curve parts snap their endpoints to it, mirroring upstream's
# shared `MVertex`).
function _model_mesh_part_node_entities(m::GeoModel,dim::Int,tag::Int,mesh,
                                         caller::AbstractString)
    owners=fill((dim,Int32(tag)),nnodes(mesh))
    dim==0 && return owners
    if dim==1
        record=get(m.discrete,(1,tag),get(m.meshing.attached,(1,tag),nothing))
        if haskey(m.curves,tag) &&
           (record===nothing || haskey(m.curve_params,tag))
            a,b=m.curves[tag]
            if !isempty(owners)
                owners[1]=(0,Int32(a))
                a!=b && (owners[end]=(0,Int32(b)))
            end
            return owners
        end
        record===nothing && throw(ArgumentError("$caller: unknown discrete curve $tag"))
        # MeshOnlyEmpty keeps raw node classification verbatim. Coincident
        # record nodes are not the model vertex's own mesh vertex upstream.
        if !haskey(m.curve_params,tag)
            _model_assert_raw_node_identity(record,1,tag,caller)
            return owners
        end
        endpoints=if haskey(m.curves,tag)
            unique!(collect(m.curves[tag]))
        else
            Int[point for (boundary_dim,point) in record.boundary if boundary_dim==0]
        end
        # Remeshed discrete curves use their boundary vertices, while the
        # stored record need not enumerate those vertices in curve order.
        degree=zeros(Int,nnodes(mesh))
        for cell in axes(mesh.segs,2),node in @view mesh.segs[:,cell]
            degree[node]+=1
        end
        for point in endpoints
            coordinates=if haskey(m.points,point)
                (m.points[point],)
            else
                point_record=get(m.discrete,(0,point),get(m.meshing.attached,(0,point),nothing))
                point_record===nothing ? () : Tuple(
                    ntuple(axis->point_record.node_coords[axis,node],3)
                    for node in axes(point_record.node_coords,2))
            end
            for node in eachindex(owners)
                (degree[node]==1 || length(endpoints)==1 && degree[node]>0) || continue
                p=ntuple(axis->mesh.coords[axis,node],3)
                any(==(p),coordinates) || continue
                owners[node]=(0,Int32(point))
            end
        end
        return owners
    end
    dim in (2,3) || throw(ArgumentError("$caller: invalid mesh part dimension $dim"))
    projected=model_to_mixed(m,mesh,dim,tag)
    data=projected.entity_data
    data!==nothing && length(data.node_entities)==nnodes(mesh) ||
        throw(ErrorException("$caller: mesh part projection changed the node set"))
    return copy(data.node_entities)
end

function _model_point_mesh_parts(m::GeoModel,caller::AbstractString)
    _model_require_point_mesh_geometry(m,caller)
    parts=Tuple{Int,Int,Mesh}[]
    for (tag,point) in m.points
        record=get(m.meshing.attached,(0,tag),nothing)
        record!==nothing && !isempty(record.node_tags) && continue
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
        _model_assert_raw_node_identity(record,dim,tag,caller)
        push!(parts,(0,tag,Mesh(Matrix{Float64}(record.node_coords[:,1:n]))))
    end
    return parts
end

# Evaluate a stored curve parameter for element emission. Every kind goes
# through `_periodic_curve_point` — the same affine evaluator the surface
# PSLG uses — so curve-part nodes deduplicate bitwise against the boundary
# nodes surface meshing generates (the merge keys on `reinterpret`ed
# coordinates). A periodic slave's curved part must emit `affine(master)`
# nodes, not the slave's native evaluation: the two differ by ulps on arcs
# and splines, and emitting both would leave unpaired near-duplicate nodes.
#
# Extruded curves are the exception: upstream places their mesh nodes at
# `ep->Extrude` transform evaluations, not on the curve's own geometry —
# a point-generatrix connector emits `extrudeMesh(GVertex)` (the source
# vertex swept to each level's `u`) and a top-copy chapeau emits
# `copyMesh(GEdge)` (each generatrix node swept to `u_top`). Evaluating the
# recorded spec keeps the curve's nodes bitwise equal to the positions the
# lateral surface and volume sweeps compute; a connector whose spline end
# vertex sits off the eval path (a twist `{{T},{axis},{P},a}` connector)
# keeps that vertex only as the curve's own endpoint, like upstream.
function _model_curve_part_point(m::GeoModel,curve::Integer,u::Float64,
                                 caller::AbstractString)
    src=get(m.meshing.extrude_sources,(1,Int(curve)),nothing)
    if src!==nothing
        params=get(m.meshing.extrude,(1,Int(curve)),nothing)
        spec=get(m.meshing.extrude_specs,(1,Int(curve)),nothing)
        if params!==nothing && !isempty(params.layers) && spec!==nothing
            if src[1]==0 && haskey(m.points,src[2])
                t0,t1=_model_curve_param_bounds(m,curve,caller)
                level=t1==t0 ? 0.0 : (u-t0)/(t1-t0)
                return _extrude_at(spec,m.points[src[2]],level,caller)
            elseif src[1]==1 && (haskey(m.curves,abs(src[2])) ||
                                 haskey(m.discrete,(1,abs(src[2]))))
                gen=abs(src[2])
                usrc=u
                if src[2]<0
                    _,t1=_model_curve_param_bounds(m,curve,caller)
                    u_min=_model_curve_param_bounds(m,gen,caller)[1]
                    usrc=t1+u_min-u
                end
                u_top=_extrude_level_us(params)[end]
                return _extrude_at(spec,
                    _periodic_curve_point(m,gen,usrc,caller),u_top,caller)
            end
        end
    end
    return _periodic_curve_point(m,curve,u,caller)
end

# Segment parts from the stored `curve_params` plus any discrete (1,*)
# records' line elements. Open-curve endpoints snap to the model points —
# upstream shares the vertex `MVertex`, so endpoint nodes must deduplicate
# bitwise against the `(0, point)` parts. Closed curves (`a==b`) connect
# their last node back to the first, like upstream's shared end vertex.
function _model_curve_mesh_parts(m::GeoModel,caller::AbstractString)
    _model_require_point_mesh_geometry(m,caller)
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
        _model_assert_raw_node_identity(record,dim,tag,caller)
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
    return _mesh_boundary_from_counts(mesh,counts)
end

# Mixed variant: quadrangle perimeter edges pair identically under the
# degree-1 boundary rule — recombination never touches boundary edges.
function _mesh_boundary_classification(mesh::MixedMesh)
    counts=Dict{NTuple{2,Int32},Int}()
    for block in mesh.blocks
        block isa ElementBlock || continue
        msh_dimension(block.msh)==2 || continue
        slots=_surface_cell_edge_slots(msh_num_nodes(block.msh))
        slots===nothing && continue
        @inbounds for cell in axes(block.nodes,2),(u,v) in slots
            a=block.nodes[u,cell]; b=block.nodes[v,cell]
            key=a<b ? (a,b) : (b,a)
            counts[key]=get(counts,key,0)+1
        end
    end
    return _mesh_boundary_from_counts(mesh,counts)
end

function _mesh_boundary_from_counts(mesh,counts)
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
function _model_surface_boundary_writeback!(m::GeoModel,t::Int,mesh,
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
    # Vertices a sibling curve owns outright are excluded per curve — at a
    # boundary pinch the two sides' subdivision vertices sit inside the
    # audit's projection tolerance without sharing an edge.
    owned=_curve_owned_coordinates(
        m,Iterators.flatten((boundary_curves,embedded_curves)),caller)
    boundary_eligible(curve)=_curve_chain_eligible(
        mesh,eligible_nodes,owned,curve)
    all_eligible=trues(nnodes(mesh))
    for curve in sort!(collect(boundary_curves))
        _model_curve_periodic_involved(m,curve) || continue
        haskey(m.curves,curve) || continue
        a,b=m.curves[curve]
        a==b && continue
        entries=_model_curve_chain_entries(
            m,mesh,curve,boundary_eligible(curve),eligible_edges,caller)
        entries===nothing && continue
        for (_,node) in entries
            protected_nodes[node]=true
        end
    end
    for curve in sort!(collect(boundary_curves))
        _model_surface_curve_writeback!(
            m,curve,mesh,boundary_eligible(curve),eligible_edges,
            protected_nodes,caller)
    end
    for curve in sort!(embedded_curves)
        _model_surface_curve_writeback!(
            m,curve,mesh,
            _curve_chain_eligible(mesh,all_eligible,owned,curve),all_edges,
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
function _model_curve_chain_entries(m::GeoModel,mesh,curve::Int,
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

function _model_surface_curve_writeback!(m::GeoModel,curve::Int,mesh,
                                         eligible_nodes,eligible_edges,
                                         protected_nodes,
                                         caller::AbstractString)
    curve in m.meshing.degenerated && return nothing
    haskey(m.curves,curve) || return nothing
    a,b=m.curves[curve]
    a==b && return nothing
    _model_curve_periodic_involved(m,curve) && return nothing
    entries,parameter_tolerance=_curve_parameter_nodes(
        m,mesh,curve,eligible_nodes,eligible_edges,0.0,caller)
    # Inverse projection drifts a few ulps from the parameter the node was
    # created at; each face's writeback would ratchet the stored value and the
    # next face would emit a bitwise-different node — cracking the shared
    # boundary in a volume PLC. Snap recomputed parameters onto the stored
    # values so shared nodes keep bitwise-stable parameters; only genuinely
    # new subdivisions contribute new entries.
    bounds=_model_curve_param_bounds(m,curve,caller)
    existing=get(m.curve_params,curve,nothing)
    if existing!==nothing
        entries=[begin
            parameter=entry[1]
            match=findfirst(value->abs(value-parameter)<=parameter_tolerance,
                            existing)
            match===nothing || (parameter=existing[match])
            # A stored value within endpoint tolerance still denotes the
            # endpoint — mapping an exact bound entry back onto bound−eps
            # would evaluate one ulp off the shared corner and crack the
            # boundary.
            parameter-bounds[1]<=parameter_tolerance &&
                (parameter=bounds[1])
            bounds[2]-parameter<=parameter_tolerance &&
                (parameter=bounds[2])
            (parameter,entry[2])
        end for entry in entries]
    end
    m.curve_params[curve]=Float64[entry[1] for entry in entries]
    @inbounds for (parameter,node) in entries
        protected_nodes[node] && continue
        # Endpoints carry the model points bitwise — evaluating one ulp off
        # would crack the boundary shared with the adjacent face, which
        # emits the corner's stored coordinates.
        point=parameter==bounds[1] ? m.points[a] :
              parameter==bounds[2] ? m.points[b] :
              _periodic_curve_point(m,curve,parameter,caller)
        mesh.coords[1,node]=point[1]
        mesh.coords[2,node]=point[2]
        mesh.coords[3,node]=point[3]
    end
    return nothing
end
