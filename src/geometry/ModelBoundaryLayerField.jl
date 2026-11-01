# `BoundaryLayer Field =` topological surface pass — upstream's
# `Mesh2DWithBoundaryLayers` (`Generator.cpp`) drives `buildAdditionalPoints2D`
# (`boundaryLayersData.cpp`), which grows vertex columns off each face's
# wall-listed boundary edges before `modifyInitialMeshForBoundaryLayers`
# (`meshGFace.cpp`) emits the strip and re-fills the shrunken face. This file
# transports the resolved field parameters (thickness/ratio/hwall/fan targets)
# into the same bounded construction: `_bl_2d_offsets` emits the certified
# strip and the ordinary constrained-Delaunay fill covers the remainder.
#
# Bounded scope, matching the upstream option set the field record carries:
# - Plane surfaces (and `allow_ruled` projections) with straight *or* curved
#   wall curves — the strip is built in the face's planar coordinates, so a
#   curved wall's inward-normal columns work exactly like upstream's
#   `parFromPoint` construction (verified on circle walls: oracle ring
#   radii match the geometric schedule). Genuinely non-planar parametric
#   (u,v) column construction stays a loud rejection.
# - `Size`/`Ratio`/`Thickness`/`Quads`/`FanPointsList`/`ExcludedSurfacesList`
#   and the `BetaLaw`/`NbLayers`/`Beta` schedule apply. `IntersectMetrics`
#   and `AnisoMax` are metric-operator options that never fire on `.geo`
#   BoundaryLayer fields (no `SurfacesList`; the scalar size field ignores
#   `iIntersect` upstream too) — accepted as inert.
# - `SizesList`/`hwall_n_nodes` per-node wall heights: `getLocalInfoAtNode`
#   resolves each vertex's `hWall` (point-tag lookup, quadratic endpoint
#   interpolation along the owning curve), the column grows that vertex's
#   geometric/beta schedule, and the strip tiles each segment with `min`
#   quadrangle rows plus upstream's closing triangle. Open runs emit the
#   same `min`-paired depth: their end columns ride the adjacent edge's
#   own mesh chain (`addColumnAtTheEndOfTheBL`), so an uninjected
#   adjacent edge contributes its ordinary nodes only.
# - `PointsList`/`NodesList` end nodes drive `addBoundaryLayerPoints`
#   grading in `ModelMesh1D` (geometric schedule capped at
#   `min(thickness, 0.4·chord)` on straight edges, reduced-interval
#   regrade, `filterPoints` skipped like upstream); the injected chains
#   are what the run-end columns consume.
# - `FanPointsSizesList` per-point fan budgets apply positionally to
#   `FanPointsList` (default `Mesh.BoundaryLayerFanElements` on a length
#   mismatch); signed sweeps, `fan_wedges` and `fan_skips` cover the
#   ordinary fan and the upstream skip types (−4/auto). The cross/mix/
#   tail fan types (−1/−2/−3) stay loud rejections.
# - Transfinite, extruded, periodic-slave and parametric (OCC) surfaces
#   keep their own producers; a `BoundaryLayer` field never alters them,
#   matching upstream's dispatch order (transfinite/extruded paths run
#   before the boundary-layer pass and return first). A swept lateral
#   (`EXTRUDED_ENTITY`) counts as "not another face" when a wall's
#   `setupFor2d` adjacency is audited, so a BL wall shared by a source
#   plane and its extruded sibling still grows the source's layer —
#   and a transfinite *curve* composes with `PointsList` injection
#   (`addBoundaryLayerPoints` runs before `meshGEdgeProcessing`, which
#   keeps `N = nbPointsTransfinite` over the reduced window).

# Per-selected-field normalized record (`_ModelBoundaryLayerField` is
# declared in `Model.jl` next to `ModelMeshingAttributes`, which stores the
# resolved map). Entity tags are literal model tags; the per-face curve
# intersection upstream computes in `setupFor2d` happens at mesh time in
# `_model_surface_bl_specs` below.
function _ModelBoundaryLayerField(tag::Integer,config,caller::AbstractString)
    # Upstream's `fanSize(iV)` reads `FanPointsSizesList[iV]` only when the
    # list matches `FanPointsList` positionally; a length mismatch silently
    # falls back to `Mesh.BoundaryLayerFanElements`, so the record keeps the
    # raw list and the consumer applies the same match rule.
    return _ModelBoundaryLayerField(
        Int(tag),Int.(collect(config.curves)),Int.(collect(config.points)),
        Int.(collect(config.fan_points)),Int.(collect(config.fan_sizes)),
        Int.(collect(config.excluded_surfaces)),
        Float64(config.size),Float64.(collect(config.sizes)),
        Float64(config.ratio),Float64(config.size_far),
        Float64(config.thickness),!iszero(config.quads),
        Int(config.fan_elements),!iszero(config.beta_law),
        Float64(config.beta),Int(config.nb_layers),
        !iszero(config.intersect_metrics),Float64(config.aniso_max),
        Set{Int}(Int.(config.curves)),Set{Int}(Int.(config.fan_points)),
        Set{Int}(Int.(config.excluded_surfaces)))
end

# Entity retag for `_model_identity_*` moves — entity-carrying option lists
# remap like every other attribute; the numeric schedule crosses unchanged.
function _model_identity_bl_field(bl::_ModelBoundaryLayerField,
                                  dimension::Int,old_tag::Int,new_tag::Int)
    remap(list,dim)=Int[v==old_tag && dim==dimension ? new_tag : v for v in list]
    curves=remap(bl.curves,1);fans=remap(bl.fan_points,0)
    excluded=remap(bl.excluded_surfaces,2);points=remap(bl.points,0)
    return _ModelBoundaryLayerField(
        bl.tag,curves,points,fans,copy(bl.fan_sizes),excluded,
        bl.hwall,copy(bl.sizes),bl.ratio,bl.hfar,bl.thickness,bl.quads,
        bl.fan_elements,bl.beta_law,bl.beta,bl.nb_layers,
        bl.intersect_metrics,bl.aniso_max,
        Set{Int}(curves),Set{Int}(fans),Set{Int}(excluded))
end

# Geometric row schedule (`L = hWall; while L <= thickness { push L;
# L += hWall·ratio^ith }` — boundaryLayersData.cpp's non-beta branch) and its
# Beta-law alternative (`t[i] = 1+β(1−p)/(1+p)`, `L = hWall·t[i]/t[0]`).
# `hwall` is the local wall height upstream's `getLocalInfoAtNode` resolves
# per column vertex — the field's `Size` unless `SizesList` grading applies.
function _bl_field_offsets(bl::_ModelBoundaryLayerField,hwall::Float64,
                           caller::AbstractString)
    if bl.beta_law
        beta=bl.beta
        (isfinite(beta) && beta>1) || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(bl.tag)].Beta must exceed 1"))
        nl=bl.nb_layers
        nl>=1 || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(bl.tag)].NbLayers must be positive"))
        zlog=_gm_log((1+beta)/(beta-1))
        t=Vector{Float64}(undef,nl)
        for i in 1:nl
            p=_gm_exp(zlog*(1-i/nl));t[i]=1+beta*(1-p)/(1+p)
        end
        return Float64[hwall*t[i]/t[1] for i in 1:nl]
    end
    hw=hwall
    (isfinite(hw) && hw>0) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)].Size must be positive"))
    (isfinite(bl.ratio) && bl.ratio>1) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)].Ratio must exceed 1"))
    (isfinite(bl.thickness) && bl.thickness>0) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)].Thickness must be positive"))
    offsets=Float64[]
    L=hw;k=0
    while L<=bl.thickness
        push!(offsets,L)
        k+=1
        # `pow(ratio, ith)` — C `pow` semantics, not Julia `^`.
        L+=hw*_gm_pow(bl.ratio,Float64(k))
        length(offsets)>100_000 && throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(bl.tag)] schedule does not " *
            "reach Thickness=$(bl.thickness)"))
    end
    isempty(offsets) && throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)] Thickness=" *
        "$(bl.thickness) admits no layer row"))
    return offsets
end

# `SizesList`/`hwall_n_nodes` — upstream's flat `{tag,h,...}` pair list feeding
# `BoundaryLayerField::hWall(iV)`; a point tag not listed falls back to
# `hWallN` (the field's `Size`). Like upstream's list scan, a repeated tag
# keeps its first value.
function _bl_hwall_map(spec::_ModelBoundaryLayerField,caller::AbstractString)
    sizes=spec.sizes
    isempty(sizes) && return Dict{Int,Float64}()
    iseven(length(sizes)) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(spec.tag)].SizesList " *
        "(hwall_n_nodes) must list point-tag/height pairs"))
    hmap=Dict{Int,Float64}()
    for i in 1:2:length(sizes)
        tagf=sizes[i];h=sizes[i+1]
        (isfinite(tagf) && abs(tagf)<9.0e15) || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].SizesList " *
            "point tag $tagf is not a finite integer"))
        ptag=round(Int,tagf)
        ptag==tagf && ptag>=1 || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].SizesList " *
            "point tag $tagf is not a positive integer"))
        (isfinite(h) && h>0) || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].SizesList " *
            "height $h must be positive"))
        haskey(hmap,ptag) || (hmap[ptag]=h)
    end
    return hmap
end

# `getLocalInfoAtNode` (boundaryLayersData.cpp): a dim-0 vertex reads
# `hWall(point_tag)` directly; a dim-1 subdivision vertex interpolates its
# curve's endpoint heights *quadratically* through the minimum midpoint —
# `hbeg·(1−x)² + 2·min(hbeg,hend)·x(1−x) + hend·x²` over the normalized curve
# parameter. Anything else falls back to the field's `Size`.
function _bl_vertex_hwall(v::Int,point_of::Dict{Int,Int},
                          vert_param::Dict{Int,Float64},
                          vert_curve::Dict{Int,Int},
                          hmap::Dict{Int,Float64},default::Float64,
                          m::GeoModel,caller::AbstractString)
    pid=get(point_of,v,0)
    pid!=0 && return get(hmap,pid,default)
    parameter=get(vert_param,v,nothing)
    if parameter!==nothing
        curve=vert_curve[v]
        a,b=m.curves[curve]
        t0,t1=_model_curve_param_bounds(m,curve,caller)
        t1>t0 || throw(ArgumentError(
            "$caller: Curve[$curve] has a degenerate parameter range"))
        x=(parameter-t0)/(t1-t0)
        hb=get(hmap,a,default);he=get(hmap,b,default)
        return hb*(1-x)*(1-x)+min(hb,he)*2*x*(1-x)+he*x*x
    end
    return default
end

# `qmTriangle::gamma` (qualityMeasures.cpp) — `2·r_in/R` — evaluated with
# upstream's own arithmetic so the comparison bits match: `norm3` is a plain
# `sqrt(x²+y²+z²)` (not a scaled `hypot`), `norme` multiplies by the
# reciprocal `1/n`, and `prodve` keeps its component order
# `c=(a1·b2-a2·b1, a2·b0-a0·b2, a0·b1-a1·b0)`. `quadsToTriangles` evaluates
# each candidate split on the quad's *stored* vertex tuple, so this routine
# is deliberately order-sensitive: on the upstream-faithful tuples emitted
# below it reproduces upstream's comparisons — including near-ties decided
# by rounding — instead of manufacturing exact ties upstream never saw.
function _bl_tri_gamma(p,q,r)
    ax,ay,az=r[1]-q[1],r[2]-q[2],r[3]-q[3]
    bx,by,bz=p[1]-r[1],p[2]-r[2],p[3]-r[3]
    cx,cy,cz=q[1]-p[1],q[2]-p[2],q[3]-p[3]
    na=sqrt(ax*ax+ay*ay+az*az);nb=sqrt(bx*bx+by*by+bz*bz)
    nc=sqrt(cx*cx+cy*cy+cz*cz)
    (na>0 && nb>0 && nc>0) || return 0.0
    oa=1.0/na;ax*=oa;ay*=oa;az*=oa
    ob=1.0/nb;bx*=ob;by*=ob;bz*=ob
    oc=1.0/nc;cx*=oc;cy*=oc;cz*=oc
    p1x,p1y,p1z=by*cz-bz*cy,bz*cx-bx*cz,bx*cy-by*cx
    p2x,p2y,p2z=cy*az-cz*ay,cz*ax-cx*az,cx*ay-cy*ax
    p3x,p3y,p3z=ay*bz-az*by,az*bx-ax*bz,ax*by-ay*bx
    sina=sqrt(p1x*p1x+p1y*p1y+p1z*p1z)
    sinb=sqrt(p2x*p2x+p2y*p2y+p2z*p2z)
    sinc=sqrt(p3x*p3x+p3y*p3y+p3z*p3z)
    (sina==0.0 && sinb==0.0 && sinc==0.0) && return 0.0
    return 2*(2*sina*sinb*sinc/(sina+sinb+sinc))
end

# `MFace::normal` (MFace.cpp via `normal3points`): `norme(prodve(v1-v0,
# v2-v0))` — normalized, so `quadsToTriangles`' `dot(normal,normal) < 0`
# test sees upstream's bits rather than an unnormalized cross product.
function _bl_tri_normal(p,q,r)
    ux,uy,uz=q[1]-p[1],q[2]-p[2],q[3]-p[3]
    vx,vy,vz=r[1]-p[1],r[2]-p[2],r[3]-p[3]
    nx,ny,nz=uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx
    n=sqrt(nx*nx+ny*ny+nz*nz)
    n==0.0 && return (0.0,0.0,0.0)
    on=1.0/n
    return (nx*on,ny*on,nz*on)
end

# ── 1-D end-node injection (`addBoundaryLayerPoints`/`createPoints`) ────────
#
# Upstream runs `addBoundaryLayerPoints` at the top of `meshGEdge` for every
# curve: a curve listed in *any* field's `CurvesList` skips injection (its
# endpoints belong to the wall), then each `PointsList` end node grows a
# geometric/beta point schedule along the endpoint chord capped by the
# layer `Thickness` and `0.4·LEdge`, and the ordinary discretization fills
# only the reduced middle interval (`filterPoints` is skipped whenever BL
# points exist).
function _model_curve_bl_end_nodes(m::GeoModel,curve::Integer,a::Int,b::Int,
                                   caller::AbstractString)
    fields=m.meshing.boundary_layer_fields
    isempty(fields) && return nothing
    spec_a=nothing;spec_b=nothing
    for spec in values(fields)
        # `setupFor1d(iE)`: a field whose `CurvesList` already walls this
        # edge leaves `_pointTags` empty — a wall never receives the same
        # field's `PointsList` injection. Another field's point end can
        # still land on the edge, matching upstream's per-field gate.
        curve in spec.curve_set && continue
        if a in spec.points
            spec_a===nothing || throw(ArgumentError(
                "$caller: Different boundary layers cannot touch each " *
                "other (curve $curve begin point $a)"))
            spec_a=spec
        end
        if b in spec.points
            spec_b===nothing || throw(ArgumentError(
                "$caller: Different boundary layers cannot touch each " *
                "other (curve $curve end point $b)"))
            spec_b=spec
        end
    end
    (spec_a===nothing && spec_b===nothing) && return nothing
    return (spec_a,spec_b)
end

# `createPoints`: points at distance `L` along the chord through the end
# vertex — the geometric law emits `L` while `L ≤ min(Thickness, 0.4·chord)`
# with `L += hWall·ratio^ith`; the beta law emits all `NbLayers` rows
# uncapped. For a `Line` the `parFromPoint` projection is linear in the
# chord fraction.
function _model_curve_bl_point_params(spec,hwall::Float64,chord::Float64,
                                      from_begin::Bool,t0::Float64,
                                      t1::Float64,caller::AbstractString)
    dists=if spec.beta_law
        _bl_field_offsets(spec,hwall,caller)
    else
        out=Float64[]
        L=hwall;k=0
        while L<=spec.thickness && L<=0.4*chord
            push!(out,L)
            k+=1
            L+=hwall*_gm_pow(spec.ratio,Float64(k))
            length(out)>100_000 && throw(ArgumentError(
                "$caller: BoundaryLayer Field[$(spec.tag)] end-point " *
                "schedule does not converge"))
        end
        out
    end
    return Float64[from_begin ? t0+d/chord*(t1-t0) :
                              t1-d/chord*(t1-t0) for d in dists]
end

function _model_curve_bl_grade_params(m::GeoModel,curve::Integer,
                                      spec_begin,spec_end,
                                      t0::Float64,t1::Float64,chord::Float64,
                                      options::_ModelMesh1DOptions,
                                      caller::AbstractString;
                                      transfinite=nothing,
                                      _native_producer::Bool=false)
    a,b=m.curves[curve]
    begin_params=Float64[];end_params=Float64[]
    if spec_begin!==nothing
        hmap=_bl_hwall_map(spec_begin,caller)
        begin_params=_model_curve_bl_point_params(
            spec_begin,get(hmap,a,spec_begin.hwall),chord,true,t0,t1,caller)
    end
    if spec_end!==nothing
        hmap=_bl_hwall_map(spec_end,caller)
        end_params=_model_curve_bl_point_params(
            spec_end,get(hmap,b,spec_end.hwall),chord,false,t0,t1,caller)
    end
    tb=isempty(begin_params) ? t0 : begin_params[end]
    te=isempty(end_params) ? t1 : end_params[end]
    tb<te || throw(ArgumentError(
        "$caller: boundary layers overlap on Curve[$curve] — upstream " *
        "requires separated layers"))
    # `meshGEdgeProcessing` on the reduced interval — the vertex-size
    # interpolation keeps the *original* param bounds (upstream's `F_Lc`
    # reads `ge->parBounds`), so only the integration window shrinks, and
    # `filterPoints` is skipped outright once BL points exist on the edge.
    # A transfinite edge keeps `N = nbPointsTransfinite` over the reduced
    # window (the BL chains then prepend/append like upstream's
    # `mesh_vertices` splice).
    middle=if transfinite===nothing
        _model_curve_grade_params(m,curve,tb,te,options,caller;
            _native_producer=_native_producer,_l1_bounds=(t0,t1),
            _skip_filter=true)
    else
        _model_curve_transfinite_params(m,curve,tb,te,transfinite,caller,
                                        options;_native_producer)
    end
    inner=length(middle)>2 ? middle[2:end-1] : Float64[]
    return Float64[t0;begin_params;inner;reverse(end_params);t1]
end

# Fields whose `CurvesList` intersects this face's boundary and whose
# `ExcludedSurfacesList` does not suppress it — `setupFor2d(face)`'s
# `_curveTags ∩ face.edges` test. A field listing only `PointsList` end
# nodes contributes nothing topologically here. Like upstream
# `setupFor2d`, a wall curve shared by more than one un-excluded,
# non-extruded surface is dropped from the face's copy of that field's
# `_curveTags` — `Msg::Error` ("Only 2D Boundary Layers are supported")
# reports the skip while the mesh continues on the surviving curves. The
# drop is per field: each spec counts adjacencies against its own
# `excluded_set`. Returns `(specs, drops)` — `drops[i]` holds the curve
# tags `specs[i]` lost on this face — or `nothing` when no field keeps a
# surviving wall curve here.
function _model_surface_bl_specs(m::GeoModel,t::Int;report_error=nothing)
    fields=m.meshing.boundary_layer_fields
    isempty(fields) && return nothing
    # A swept lateral (`EXTRUDED_ENTITY` — a `(1, ·)` extrude-sources link)
    # carries BL curves in its boundary but upstream's column pass produces
    # no effective layer there: its end-of-BL tangents land on `GVertex`s
    # (the generatrix endpoints), so no end columns exist and the per-edge
    # `min`-depth pairing degenerates to filtered self-paired cells. Skip
    # the layer on the lateral itself — the plane source face still grows
    # its strip thanks to the adjacency exemption below.
    src=get(m.meshing.extrude_sources,(2,t),nothing)
    src!==nothing && src[1]==1 && return nothing
    # Upstream iterates the face's own edge list (`gf->edges()`), so the
    # per-curve diagnostics follow boundary order, not field order.
    boundary=Int[]
    for loop in m.surfaces[t],signed in m.loops[loop]
        curve=abs(signed)
        curve in boundary || push!(boundary,curve)
    end
    specs=_ModelBoundaryLayerField[]
    drops=Set{Int}[]
    for tag in sort!(collect(keys(fields)))
        spec=fields[tag]
        t in spec.excluded_set && continue
        shared=[curve for curve in boundary if curve in spec.curve_set]
        isempty(shared) && continue
        # `setupFor2d` counts the wall curve's adjacent faces that are
        # neither extruded nor excluded; more than one means the layer
        # would grow inside a volume — upstream logs the error and leaves
        # the curve out of that face's `_curveTags` rather than failing
        # the mesh. Faces carrying an `EXTRUDED_ENTITY` mode (a `(1, ·)`
        # `extrude_sources` link — a swept lateral, not a `COPIED_ENTITY`
        # top) get a free pass: they are meshed by the extrusion path,
        # which never builds columns, so the layer grows on the single
        # genuine face only.
        drop=Set{Int}()
        for curve in shared
            numf=0;nadj=0
            for (stag,loops) in m.surfaces
                any(loop->any(signed->abs(signed)==curve,
                              m.loops[loop]),loops) || continue
                nadj+=1  # upstream prints `fc.size()` — every adjacency
                stag in spec.excluded_set && continue
                src=get(m.meshing.extrude_sources,(2,stag),nothing)
                src!==nothing && src[1]==1 && continue
                numf+=1
            end
            numf>1 || continue
            report_error===nothing || report_error(
                "Only 2D Boundary Layers are supported (curve $curve " *
                "is adjacent to $nadj surfaces)")
            push!(drop,curve)
        end
        # A field whose wall curves all fell out of this face's `_curveTags`
        # contributes no layer here — upstream's emptied field leaves the
        # face's strip pass a no-op, so take the ordinary fill path.
        any(curve->!(curve in drop),shared) || continue
        push!(specs,spec);push!(drops,drop)
    end
    return isempty(specs) ? nothing : (specs,drops)
end

# `_surface_pslg` twin that additionally records the ordered wire vertex
# lists and each boundary segment's source curve — the wall-run extraction
# below needs the curve ownership of every boundary segment.
function _bl_surface_pslg(m::GeoModel,t::Int,forced,caller::AbstractString;
                          param_sizes::Dict{Tuple{Int,Float64},Float64}=
                              Dict{Tuple{Int,Float64},Float64}(),
                          plane=_model_surface_plane(m,t,caller;
                                                     allow_ruled=true),
                          _periodic_positions=nothing,
                          uv::Bool=false)
    xs=Float64[];ys=Float64[];mesh_sizes=Float64[]
    index=Dict{Int,Int}()
    canonical=Dict{Tuple{Float64,Float64},NTuple{3,Float64}}()
    wires=Vector{Int}[]
    wire_curves=Vector{Int}[]
    # Per-segment sign of the owning curve's loop traversal — upstream's
    # `MQuadrangle(v11,v21,v22,v12)` eval tuple follows each edge's `lines`
    # in parameter order, so whether the kernel's walk agrees with that
    # direction (per segment, through the chain reversal) fixes which
    # stored order the `quadsToTriangles` split must evaluate.
    wire_signs=Vector{Int}[]
    vert_curve=Dict{Int,Int}()
    # Subdivision vertex → its raw curve parameter — upstream's
    # `getLocalInfoAtNode` interpolates `hWall` along the curve parameter for
    # dim-1 mesh vertices.
    vert_param=Dict{Int,Float64}()
    # Upstream always builds `Mesh2DWithBoundaryLayers` on the face's own
    # parametrization — `reparamMeshEdgeOnFace`/`parFromPoint` give the
    # working coordinates, and the lifted strip lands back through the
    # surface evaluator. Plane surfaces take the identical working frame
    # through the projection; ruled/tric patches take their `(u,v)` chart
    # (`uv=true`), where every generatrix edge is straight — the case the
    # plane projection cannot represent.
    lc=uv ? _model_lc(m) : 0.0
    add_point! = function (pid::Int)
        uv && return _add_ruled_surface_point!(
            xs,ys,mesh_sizes,index,canonical,m,t,pid,lc,caller)
        return _add_surface_point!(xs,ys,mesh_sizes,index,canonical,m,pid,
                                   caller,plane)
    end
    add_curve_point! = function (point,mesh_size::Float64,curve::Int,
                                 parameter::Float64)
        uv && return _add_ruled_surface_curve_point!(
            xs,ys,mesh_sizes,index,canonical,m,t,point,mesh_size,curve,
            parameter,lc,caller)
        return _add_surface_curve_point!(xs,ys,mesh_sizes,index,canonical,m,
                                         point,mesh_size,curve,caller,plane)
    end
    for loop_id in m.surfaces[t]
        _verify_loop_closed(m,loop_id,caller,"Surface[$t]")
        loop_idx=Int[]
        loop_curve=Int[]
        loop_sign=Int[]
        for signed in m.loops[loop_id]
            curve=abs(signed)
            a,b=m.curves[curve]
            t0,t1=_model_curve_param_bounds(m,curve,caller)
            for parameter in _surface_curve_parameters(m,forced,curve,signed,
                                                       caller)
                vertex=if parameter==t0
                    add_point!(a)
                elseif parameter==t1
                    add_point!(b)
                else
                    point=_periodic_forced_curve_point(
                        m,curve,parameter,_periodic_positions,caller)
                    uv || begin
                        scale=max(1.0,hypot(point...))
                        abs(_plane_offset(plane,point))<=1e-12*scale ||
                            throw(ArgumentError(
                                "$caller: Curve[$curve] subdivision is not " *
                                "coplanar with Surface[$t]"))
                    end
                    add_curve_point!(
                        point,
                        get(param_sizes,(curve,parameter),
                            _surface_curve_mesh_size(
                                m,curve,parameter,caller)),
                        curve,parameter)
                end
                if !isempty(loop_idx) && loop_idx[end]==vertex
                    # Corner shared by two boundary curves: the slot
                    # owns the outgoing segment, so the continuing
                    # curve's tag wins.
                    loop_curve[end]=curve
                    loop_sign[end]=signed>0 ? 1 : -1
                else
                    push!(loop_idx,vertex)
                    push!(loop_curve,curve)
                    push!(loop_sign,signed>0 ? 1 : -1)
                end
                parameter==t0 || parameter==t1 || begin
                    haskey(vert_curve,vertex) || (vert_curve[vertex]=curve)
                    haskey(vert_param,vertex) || (vert_param[vertex]=parameter)
                end
            end
        end
        while length(loop_idx)>1 && loop_idx[end]==loop_idx[1]
            pop!(loop_idx);pop!(loop_curve);pop!(loop_sign)
        end
        nloop=length(loop_idx)
        nloop>=3 || throw(ArgumentError(
            "$caller: Loop[$loop_id] needs at least three points"))
        push!(wires,loop_idx)
        push!(wire_curves,loop_curve)
        push!(wire_signs,loop_sign)
    end
    embedded=get(m.embeds,(2,t),NTuple{2,Int}[])
    internal=Tuple{Int,Int}[]
    for (edim,etag) in embedded
        if edim==0
            add_point!(etag)
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
            push!(curve_nodes,add_point!(a))
            push!(curve_nodes,add_point!(b))
        else
            for parameter in parameters
                vertex=if parameter==t0
                    add_point!(a)
                elseif parameter==t1
                    add_point!(b)
                else
                    point=_periodic_forced_curve_point(
                        m,etag,parameter,_periodic_positions,caller)
                    uv || begin
                        scale=max(1.0,hypot(point...))
                        abs(_plane_offset(plane,point))<=1e-12*scale ||
                            throw(ArgumentError(
                                "$caller: embedded Curve[$etag] subdivision " *
                                "is not coplanar with Surface[$t]"))
                    end
                    add_curve_point!(
                        point,
                        get(param_sizes,(etag,parameter),
                            _surface_curve_mesh_size(
                                m,etag,parameter,caller)),
                        etag,parameter)
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
        for node in curve_nodes
            haskey(vert_curve,node) || (vert_curve[node]=etag)
        end
    end
    return xs,ys,mesh_sizes,wires,wire_curves,wire_signs,embedded,internal,
           index,canonical,vert_curve,vert_param
end

# Even-odd containment over the face's wires — winding-insensitive, so hole
# wires and reversed outer wires classify identically (`inside` = contained
# by an odd number of wires).
function _bl_ray_cast(px::Float64,py::Float64,
                      xs::Vector{Float64},ys::Vector{Float64},
                      wire::Vector{Int})
    n=length(wire);inside=false
    for k in 1:n
        i=wire[k];j=wire[mod1(k+1,n)]
        y1=ys[i];y2=ys[j]
        (y1>py)!=(y2>py) || continue
        xint=xs[i]+(py-y1)*(xs[j]-xs[i])/(y2-y1)
        xint>px && (inside=!inside)
    end
    return inside
end

function _bl_point_in_domain(px::Float64,py::Float64,
                             xs::Vector{Float64},ys::Vector{Float64},
                             wires::Vector{Vector{Int}})
    count=0
    for wire in wires
        _bl_ray_cast(px,py,xs,ys,wire) && (count+=1)
    end
    return isodd(count)
end

# Collect a run end's column chain along the wire — upstream's
# `addColumnAtTheEndOfTheBL` makes the end vertex's column the adjacent
# edge's own node chain verbatim (the boundary write-back BL-graded the
# adjacent curve near the end through `addBoundaryLayerPoints`, so the
# chain already holds the injected layer points). The walk consumes every
# wire position interior to `adjacent_curve` — a vertex stays in the wire
# when the segment after it in walk direction leaves the curve (the edge's
# far endpoint or a curve junction). Returns `(positions, dists, dir3)`:
# consumed rotated positions nearest→farthest, their cumulative 2-D
# distances from the end vertex (the column's `vertex_offsets` entry — on
# a straight edge the 2-D and 3-D distances agree), and the unit 3-D
# direction of the first consumed segment — the `vertex_dirs` override.
# Bounded: the column direction must stay straight within tolerance —
# bending end columns are a documented gap.
function _bl_run_end_chain(dsign::Int,p_end::Int,adjacent_curve::Int,
                           spec,caller::AbstractString,t::Int,
                           seg_spec::Vector{Int},loop::Vector{Int},
                           loop_curve::Vector{Int},xs::Vector{Float64},
                           ys::Vector{Float64},
                           coordsw::Vector{NTuple{3,Float64}},
                           start0::Int,nseg::Int)
    pmod(p)=mod1(p,nseg)
    rvert(p)=loop[mod1(start0+p-1,nseg)]
    positions=Int[];dists=Float64[]
    dir3=nothing
    traveled=0.0;d2x=0.0;d2y=0.0
    pos=pmod(p_end)
    while true
        npos=pmod(pos+dsign)
        sidx=dsign>0 ? mod1(start0+pos-1,nseg) :
                      mod1(start0+npos-1,nseg)
        (loop_curve[sidx]==adjacent_curve && seg_spec[sidx]==0) || break
        a=rvert(pos);b=rvert(npos)
        s2x=xs[b]-xs[a];s2y=ys[b]-ys[a]
        sl2=hypot(s2x,s2y)
        sl2>0 || throw(ArgumentError(
            "$caller: Surface[$t] boundary has a coincident vertex"))
        if dir3===nothing
            s3x=coordsw[b][1]-coordsw[a][1]
            s3y=coordsw[b][2]-coordsw[a][2]
            s3z=coordsw[b][3]-coordsw[a][3]
            sl3=hypot(s3x,s3y,s3z)
            dir3=(s3x/sl3,s3y/sl3,s3z/sl3)
            d2x=s2x/sl2;d2y=s2y/sl2
        else
            cross=s2x*d2y-s2y*d2x
            abs(cross)<=1e-9*sl2 || throw(ArgumentError(
                "$caller: BoundaryLayer Field[$(spec.tag)] end column on " *
                "Surface[$t] bends along Curve[$adjacent_curve] — " *
                "bending end columns are not implemented"))
        end
        # `npos` joins the column only when the segment after it (in walk
        # direction) stays on the same non-wall curve — otherwise it is
        # the adjacent edge's far endpoint and remains boundary.
        farpos=pmod(npos+dsign)
        nsidx=dsign>0 ? mod1(start0+npos-1,nseg) :
                        mod1(start0+farpos-1,nseg)
        (loop_curve[nsidx]==adjacent_curve && seg_spec[nsidx]==0) || break
        traveled+=sl2
        push!(positions,npos);push!(dists,traveled)
        pos=npos
    end
    return positions,dists,dir3
end

# Mark an end column's *emitted* nodes as covered — upstream removes the
# used chain prefix from the fill wire (it stays domain boundary since the
# strip cells ride it). `emit` is the strip's used depth
# `min(K, neighbour_depth+1)`; chain positions beyond it keep their place
# in the reduced wire. Two run ends meeting on one adjacent edge may
# share exactly their deepest node (each column a single level — the
# `meet_nodes` case upstream writes when the covered spans partition the
# edge); any deeper sharing is an overlapping-layers rejection.
function _bl_mark_end_covered!(positions::Vector{Int},emit::Int,
                               covered_positions::Set{Int},
                               claimed::Set{Int},meet_nodes::Set{Int},
                               caller::AbstractString,t::Int)
    for (j,pos) in enumerate(positions)
        j>emit && break
        if pos in claimed
            (emit==1 && pos in meet_nodes) && continue
            throw(ArgumentError(
                "$caller: Surface[$t] boundary-layer run ends consume " *
                "overlapping spans — upstream requires separated layers"))
        end
        push!(claimed,pos);push!(covered_positions,pos)
    end
    emit>0 && push!(meet_nodes,positions[emit])
    return nothing
end

# Strip + fill construction for one face. Returns Mesh (all triangles) or
# MixedMesh (triangles + kept quadrangles when the field's `Quads` or the
# face's recombine flag retains them, like upstream's
# `splitElementsInBoundaryLayerIfNeeded` gate).
function _model_bl_surface_mesh(m::GeoModel,t::Int,
                                specs::Vector{_ModelBoundaryLayerField},
                                forced,
                                param_sizes::Dict{Tuple{Int,Float64},Float64},
                                plane,caller::AbstractString;
                                _periodic_positions=nothing,
                                uv::Bool=false,
                                dropped_curves::Vector{Set{Int}}=Set{Int}[])
    # `IntersectMetrics`/`AnisoMax` are metric-operator options — upstream
    # reads `iIntersect` only inside the anisotropic `operator()(SMetric3&)`
    # (attractor-metric intersection instead of closest-only) and
    # `tgtAnisoRatio` only inside its dim-2 *surface*-attractor branch.
    # `.geo` BoundaryLayer fields admit no `SurfacesList`, so `AnisoMax`
    # can never fire, and the strip pass itself never evaluates either
    # operator — both options are inert for `BoundaryLayer Field =`
    # meshes exactly as upstream. The Tessella `BoundaryLayerField` size
    # field is scalar like upstream's scalar `operator()`, which likewise
    # ignores `iIntersect`, so `Background Field`/`Min`-style scalar
    # composition matches too.
    xs,ys,mesh_sizes,wires,wire_curves,wire_signs,embedded,internal,index,
        canonical,vert_curve,vert_param=
        _bl_surface_pslg(m,t,forced,caller;param_sizes=param_sizes,
                         plane=plane,_periodic_positions=_periodic_positions,
                         uv=uv)
    ax=plane.axes
    coords3d=Vector{NTuple{3,Float64}}(undef,length(xs))
    point_of=Dict{Int,Int}()
    for (pid,v) in index
        point_of[v]=pid
    end
    for i in eachindex(xs)
        key=(xs[i]==0 ? 0.0 : xs[i],ys[i]==0 ? 0.0 : ys[i])
        coordinate=get(canonical,key,nothing)
        coordinate===nothing && throw(ErrorException(
            "$caller: Surface[$t] boundary vertex has no canonical " *
            "coordinate"))
        coords3d[i]=coordinate
    end
    # Working-frame coordinates the strip kernel sees: the plane projection
    # keeps the physical 3-D positions (its `(x,y)` is just a dropped-axis
    # read), while the ruled `(u,v)` frame gives the kernel `(u,v,0)` nodes
    # and lifts them back through the surface evaluator — upstream's
    # `Mesh2DWithBoundaryLayers` runs the identical pass on the face's own
    # parameter coordinates (`parFromPoint`/`reparamMeshEdgeOnFace`).
    coordsw=uv ?
        NTuple{3,Float64}[(xs[i],ys[i],0.0) for i in eachindex(xs)] :
        coords3d
    lift=uv ?
        c->_ruled_surface_point(m,t,c[1],c[2],caller) : identity
    # Component pair that recovers a working node coordinate's `(x,y)` —
    # dropped-plane axes in the plane frame, the literal `(u,v)` in the
    # ruled frame.
    axw=uv ? (1,2) : ax
    # The surface's true unit normal, signed so that a positive rotation
    # around it matches the left-hand side in the (x,y) projection frame.
    # The `(u,v)` frame is literally two-dimensional — its normal is the
    # parametric +z, and the domain-side probe below decides which side of
    # a wall the strip grows on regardless of the parametrization's
    # handedness.
    plane_normal=if uv
        (0.0,0.0,1.0)
    else
        nvec=plane.normal
        nlen=plane.norm
        nsurf=(nvec[1]/nlen,nvec[2]/nlen,nvec[3]/nlen)
        nuv=ax==(1,2) ? (0.0,0.0,1.0) :
            ax==(1,3) ? (0.0,-1.0,0.0) : (1.0,0.0,0.0)
        flip=nsurf[1]*nuv[1]+nsurf[2]*nuv[2]+nsurf[3]*nuv[3]
        flip>0 ? nsurf : (-nsurf[1],-nsurf[2],-nsurf[3])
    end
    # Upstream's `parFromPoint`/`point` parametrization for builtin planar
    # faces — `GFace::computeMeanPlane`'s edge-sample fast path plus
    # `fillMeanPlane` — so the strip kernel's working coordinates, emitted
    # nodes and stored surface params all live in upstream's frame rather
    # than a locally orthonormalized one. `nothing` when the fast path
    # cannot run (upstream falls back to its SVD mean plane, which the
    # kernel's own orthonormal basis approximates).
    pframe=uv ? nothing : _model_surface_meanplane_frame(m,t,caller)
    # Mean-plane `(u,v)` of a physical point — upstream `parFromPoint`
    # (`prosca` component order) — for `(2,t)`-owned vertices' params.
    _mp_param=pframe===nothing ? nothing : function (c)
        O,u,v=pframe
        d1=c[1]-O[1]; d2=c[2]-O[2]; d3=c[3]-O[3]
        return (d1*u[1]+d2*u[2]+d3*u[3],
                d1*v[1]+d2*v[2]+d3*v[3])
    end
    blmod=getfield(parentmodule(@__MODULE__),:BoundaryLayer)
    new_xs=Float64[];new_ys=Float64[]
    new_coords=NTuple{3,Float64}[]
    new_owners=Tuple{Int,Int32}[]
    new_params=Vector{Union{Nothing,Vector{Float64}}}()
    vertex_index=Dict{NTuple{2,Float64},Int}()
    # Nodes a run end column emits on the adjacent non-wall segment —
    # upstream's `addColumnAtTheEndOfTheBL` makes that column the adjacent
    # edge's own mesh chain, so those coordinates are owned by the
    # adjacent curve, not the surface.
    covered_owner=Dict{NTuple{3,Float64},Int32}()
    # Exact wire-vertex (or interpolated top) coordinates that strip
    # column nodes must weld onto bitwise — column nodes computed from
    # `vertex_dirs` may differ by a rounding ulp from the covered wire
    # vertices the boundary write-back produced.
    snapped_coords=NTuple{3,Float64}[]
    function vertex_id!(i::Int)
        key=(xs[i]==0 ? 0.0 : xs[i],ys[i]==0 ? 0.0 : ys[i])
        return get!(vertex_index,key) do
            push!(new_xs,xs[i]);push!(new_ys,ys[i])
            push!(new_coords,coords3d[i])
            pid=get(point_of,i,0)
            if pid!=0
                push!(new_owners,(0,Int32(pid)))
                push!(new_params,nothing)
            else
                curve=get(vert_curve,i,0)
                if curve!=0
                    push!(new_owners,(1,Int32(curve)))
                    push!(new_params,nothing)
                else
                    push!(new_owners,(2,Int32(t)))
                    push!(new_params,_mp_param===nothing ?
                        [xs[i],ys[i]] : collect(_mp_param(coords3d[i])))
                end
            end
            length(new_coords)
        end
    end
    function new_vertex!(x::Float64,y::Float64,coordinate)
        key=(x==0 ? 0.0 : x,y==0 ? 0.0 : y)
        return get!(vertex_index,key) do
            push!(new_xs,x);push!(new_ys,y);push!(new_coords,coordinate)
            ctag=get(covered_owner,coordinate,Int32(0))
            if ctag!=0
                push!(new_owners,(1,ctag))
                push!(new_params,nothing)
            else
                push!(new_owners,(2,Int32(t)))
                push!(new_params,_mp_param===nothing ?
                    [x,y] : collect(_mp_param(coordinate)))
            end
            length(new_coords)
        end
    end
    strip_coords=NTuple{3,Float64}[]
    # Working-frame `(x,y)` of each strip node — the surface parameters the
    # ownership records publish for `(2,t)`-owned vertices (`(u,v)` on the
    # ruled frame; dropped-plane components on the plane frame).
    strip_params=Vector{Float64}[]
    strip_tris=NTuple{3,Int}[]
    strip_quads=NTuple{4,Int}[]
    # Parallel upstream-order tuples for `quadsToTriangles`' split
    # evaluation — the same strip-node ids as `strip_quads`, permuted to
    # each quad's `MQuadrangle`-stored order.
    strip_quads_eval=NTuple{4,Int}[]
    reduced_wires=Vector{Int}[]
    # Map each wire's segment index to its field (0 = ordinary boundary
    # segment). Segments carry the tag of the curve that emitted their first
    # vertex — `_bl_surface_pslg`'s per-run bookkeeping. `hwall_maps` caches
    # each spec's `SizesList` pair-map across its runs.
    hwall_maps=Dict{Int,Dict{Int,Float64}}()
    for (w,loop) in enumerate(wires)
        loop_curve=wire_curves[w]
        loop_sign=wire_signs[w]
        nseg=length(loop)
        seg_spec=Vector{Int}(undef,nseg)
        for k in 1:nseg
            curve=loop_curve[k]
            owner=0
            for (si,spec) in enumerate(specs)
                # Curves this face's `setupFor2d` dropped from spec `si`
                # (shared across un-excluded, non-extruded surfaces) stay
                # ordinary boundary — upstream removes them from that
                # field's `_curveTags`, so they neither grow columns nor
                # mark their segments as wall.
                if curve in spec.curve_set &&
                        !(si<=length(dropped_curves) &&
                          curve in dropped_curves[si])
                    owner=si
                    break
                end
            end
            seg_spec[k]=owner
        end
        all(==(0),seg_spec) && begin
            push!(reduced_wires,[vertex_id!(v) for v in loop])
            continue
        end
        # Seed every wire vertex so wall-chain and run-end nodes keep
        # their point/curve ownership through the strip pool — upstream
        # reuses the same MVertex objects along the covered boundary.
        for v in loop
            vertex_id!(v)
        end
        if all(!=(0),seg_spec) && !all(==(seg_spec[1]),seg_spec)
            throw(ArgumentError(
                "$caller: Surface[$t] wall segments split across " *
                "boundary-layer fields around a closed wire is not " *
                "implemented"))
        end
        # Rotate to a non-wall segment (a wholly-walled wire stays put and
        # produces the closed-ring run).
        start0=1
        if !all(!=(0),seg_spec)
            start0=findfirst(==(0),seg_spec)::Int+1
            start0>nseg && (start0=1)
        end
        # Maximal cyclic runs of same-field wall segments: (spec_index,
        # first_seg, last_seg, closed) in rotated walk order.
        runs=NTuple{4,Int}[]
        k=1
        while k<=nseg
            sk=mod1(start0+k-1,nseg)
            owner=seg_spec[sk]
            if owner==0
                k+=1
                continue
            end
            first=k
            last=k
            while k+1<=nseg && seg_spec[mod1(start0+k,nseg)]==owner
                k+=1;last=k
            end
            push!(runs,(owner,first,last,
                        first==1 && last==nseg ? 1 : 0))
            k+=1
        end
        # A closed run covers the whole wire; otherwise run segs first..last
        # span wall vertices at rotated positions first..last+1.
        top_of=Dict{Int,Int}()
        # Per-wall-vertex exposed boundary chains — the strip's free top
        # edges plus the uncancelled column rungs between two neighboring
        # segments of unequal depth (`SizesList` grading). Each entry lists
        # the fill-boundary nodes a wire vertex contributes in walk order.
        wire_nodes=Dict{Int,Vector{Int}}()
        # Run-end covered spans: wire positions consumed by the strip side
        # (skipped in the reduced wire), positions already claimed by an
        # end walk (two ends covering one vertex = touching layers), exact
        # snap targets for strip column nodes, and per-end-vertex top
        # vertex overrides.
        covered_positions=Set{Int}()
        claimed=Set{Int}()
        meet_nodes=Set{Int}()
        top_override=Dict{Int,Int}()
        # Runs whose strips emit no cells at all — every wall segment pairs
        # two empty columns and no fan exists — keep their boundary
        # verbatim (upstream's `modifyInitialMeshForBoundaryLayers` simply
        # writes no BL cells there either).
        dead_runs=Set{Int}()
        rvert(p)=loop[mod1(start0+p-1,nseg)]
        for (ri,(si,first,last,closed_flag)) in enumerate(runs)
            spec=specs[si]
            closed_run=!iszero(closed_flag)
            chain_positions=closed_run ? collect(1:nseg) :
                collect(first:last+1)
            chain=Int[loop[mod1(start0+p-1,nseg)] for p in chain_positions]
            nchain=length(chain)
            # Domain-side orientation: probe a point off the first wall
            # segment's left normal; flip the chain when the wall's left is
            # not the face interior.
            a=chain[1];b=chain[2]
            dx=xs[b]-xs[a];dy=ys[b]-ys[a]
            seglen=hypot(dx,dy)
            seglen>0 || throw(ArgumentError(
                "$caller: Surface[$t] wall curve has a coincident vertex"))
            offsets=_bl_field_offsets(spec,spec.hwall,caller)
            probe=min(0.5*offsets[1],0.05*seglen)
            lx=-dy/seglen;ly=dx/seglen
            mx=0.5*(xs[a]+xs[b])+probe*lx
            my=0.5*(ys[a]+ys[b])+probe*ly
            inside=_bl_point_in_domain(mx,my,xs,ys,wires)
            chain_rev=false
            if !inside
                reverse!(chain);chain_rev=true
                a=chain[1];b=chain[2]
            end
            # `SizesList` grading — upstream's `getLocalInfoAtNode`
            # resolves a per-vertex `hWall` (point-tag lookup, quadratic
            # interpolation along the owning curve), so each wall vertex
            # grows its own column schedule. On an *open* run the two end
            # vertices take no schedule of their own — upstream's
            # `addColumnAtTheEndOfTheBL` substitutes the adjacent edge's
            # own node chain for their column, collected below.
            hmap=get!(hwall_maps,si) do
                _bl_hwall_map(spec,caller)
            end
            # Upstream's `getLocalInfoAtNode` resolves a per-vertex `hWall`
            # for EVERY column vertex, not only graded ones: a dim-1
            # subdivision vertex always runs the quadratic endpoint
            # interpolation `hb·(1-x)² + 2·min(hb,he)·x(1-x) + he·x²` at its
            # stored curve parameter — with equal endpoint heights that
            # Bernstein still drifts a rounding ulp off the constant, and
            # those ulp shifts are observable in the emitted mesh (they
            # break `quadsToTriangles` near-ties identically upstream).
            hws=Float64[_bl_vertex_hwall(v,point_of,vert_param,vert_curve,
                                         hmap,spec.hwall,m,caller)
                        for v in chain]
            # Interpolating equal endpoint heights can drift by an ulp —
            # upstream's `getLocalInfoAtNode` emits the same values, so the
            # schedules below carry the drift per vertex rather than
            # collapsing it; only genuinely `SizesList`-scale variation
            # marks the run `variable` (it is far larger than the ulp band,
            # which would otherwise manufacture degenerate closing
            # triangles in the topology checks downstream).
            variable=any(
                h->abs(h-hws[1])>1e-12*max(abs(hws[1]),eps(Float64)),hws)
            if hws!==nothing && !variable
                # A uniformly-resolved map still overrides `Size` — every
                # vertex reports the same `getLocalInfoAtNode` hWall, so the
                # run keeps the shared schedule built from that value
                # (end columns and probes read `offsets`).
                offsets=_bl_field_offsets(spec,hws[1],caller)
                isempty(offsets) && throw(ArgumentError(
                    "$caller: BoundaryLayer Field[$(spec.tag)] hWall=" *
                    "$(hws[1]) at a wall vertex exceeds Thickness=" *
                    "$(spec.thickness) — the column admits no layer row"))
            end
            # Every wall vertex carries its own schedule — `voff[i]` is the
            # cumulative offset vector `getLocalInfoAtNode`'s hWall resolves
            # to. On an *open* run the two end slots hold the adjacent
            # edge's chain distances (overwritten below), not the shared
            # `offsets`; on a closed run the per-vertex schedules differ
            # only where `hws` does — ulp noise included, matching
            # upstream's per-column `bl_depth`/`L` accumulation.
            voff=Vector{Vector{Float64}}(undef,nchain)
            for i in 1:nchain
                hw=hws[i]
                voff[i]=_bl_field_offsets(spec,hw,caller)
                isempty(voff[i]) && throw(ArgumentError(
                    "$caller: BoundaryLayer Field[$(spec.tag)] hWall=" *
                    "$hw at a wall vertex exceeds Thickness=" *
                    "$(spec.thickness) — the column admits no layer row"))
            end
            wall_coords=Matrix{Float64}(undef,3,nchain)
            for (i,v) in enumerate(chain)
                wall_coords[1,i],wall_coords[2,i],wall_coords[3,i]=coordsw[v]
            end
            wall_segs=Matrix{Int32}(undef,2,closed_run ? nchain : nchain-1)
            for i in 1:(nchain-1)
                wall_segs[1,i]=i;wall_segs[2,i]=i+1
            end
            if closed_run
                wall_segs[1,nchain]=nchain;wall_segs[2,nchain]=1
            end
            wall_mesh=Mesh(wall_coords;segs=wall_segs)
            # The strip's base nodes lift the wall vertices back through
            # the surface evaluator — on the ruled `(u,v)` frame that
            # round-trip lands an ulp off the canonical curve-mesh
            # coordinate and would weld a *new* `(2,t)`-owned node next to
            # the seeded point/curve entry (upstream keeps the same
            # MVertex). Every wall position is therefore a bitwise snap
            # target like the covered run-end spans.
            for v in chain
                push!(snapped_coords,coords3d[v])
            end
            # Run-end columns ride the adjacent edge — upstream's
            # `addColumnAtTheEndOfTheBL` makes the end vertex's column the
            # adjacent edge's own node chain verbatim (the boundary
            # write-back BL-graded that edge near the end through
            # `addBoundaryLayerPoints`, so the injected layer points are
            # already interior wire vertices). The strip then uses the
            # chain prefix up to `min(K, neighbour+1)` levels; chain nodes
            # beyond the used depth keep their place in the reduced wire.
            # Bounded: the column direction must stay straight and two
            # ends may share only their deepest node when the covered
            # spans partition one edge.
            vdirs=nothing
            end_cols=nothing
            if !closed_run
                p_start=first
                p_end=last+1
                v_start=loop[mod1(start0+p_start-1,nseg)]
                v_end=loop[mod1(start0+p_end-1,nseg)]
                prev_seg=mod1(start0+first-2,nseg)
                next_seg=mod1(start0+last,nseg)
                seg_spec[prev_seg]==0 || throw(ArgumentError(
                    "$caller: Surface[$t] has adjacent boundary-layer " *
                    "runs — upstream requires separated layers"))
                seg_spec[next_seg]==0 || throw(ArgumentError(
                    "$caller: Surface[$t] has adjacent boundary-layer " *
                    "runs — upstream requires separated layers"))
                ch_s,dist_s,dir_s=_bl_run_end_chain(
                    -1,p_start,loop_curve[prev_seg],spec,caller,t,
                    seg_spec,loop,loop_curve,xs,ys,coordsw,start0,nseg)
                ch_e,dist_e,dir_e=_bl_run_end_chain(
                    +1,p_end,loop_curve[next_seg],spec,caller,t,
                    seg_spec,loop,loop_curve,xs,ys,coordsw,start0,nseg)
                # `chain` may have been reversed for orientation; attach
                # each chain to the slot holding its end vertex and let
                # the slot's wall-side neighbour bound the used depth
                # (`min(K, n+1)` — upstream pairs `min(na,nb)` quadrangle
                # rows plus the closing triangle's next level).
                slot_s=0;slot_e=0
                for (i,v) in enumerate(chain)
                    v==v_start && (slot_s=i)
                    v==v_end && (slot_e=i)
                end
                (slot_s==0 || slot_e==0) && throw(ErrorException(
                    "$caller: Surface[$t] boundary-layer run end is not " *
                    "on its wall chain"))
                voff[slot_s]=dist_s
                voff[slot_e]=dist_e
                nb_s=length(voff[slot_s==1 ? min(2,nchain) : slot_s-1])
                nb_e=length(voff[slot_e==nchain ? max(nchain-1,1) : slot_e+1])
                end_cols=((ch_s,min(length(dist_s),nb_s+1),v_start,
                           loop_curve[prev_seg],slot_s,dir_s),
                          (ch_e,min(length(dist_e),nb_e+1),v_end,
                           loop_curve[next_seg],slot_e,dir_e))
            end
            fans=Int[]
            for (i,v) in enumerate(chain)
                pid=get(point_of,v,0)
                pid in spec.fan_set || continue
                interior=closed_run || (i>1 && i<nchain)
                interior || continue
                push!(fans,i)
            end
            # `FanPointsSizesList` — upstream `fanSize(iV)` reads the list
            # positionally per `FanPointsList` entry, with a length mismatch
            # falling back to `Mesh.BoundaryLayerFanElements`
            # (`spec.fan_elements`). Nonnegative sizes select the standard
            # fan (`fanType -1`, `round(frac·size)+1` wedges); -1/-2/-3 select
            # upstream's cross/mix/tail topologies (`fanType` 1/2/3 — the
            # bilinear / rounded-fold / `l·n±ll·nc` column machinery through
            # `fan_types`); -4 and every other negative fall to `fanType 4`
            # (skip: boundary rays only, no wedge cells).
            vfans=nothing;fskips=nothing;ftypes=nothing
            if !isempty(fans)
                matched=length(spec.fan_sizes)==length(spec.fan_points)
                vfans=Dict{Int,Int}();fskips=Set{Int}();ftypes=Dict{Int,Int}()
                for i in fans
                    pid=get(point_of,chain[i],0)
                    f=spec.fan_elements
                    if matched
                        fpos=findfirst(==(pid),spec.fan_points)
                        fpos===nothing || (f=spec.fan_sizes[fpos])
                    end
                    if f>=0
                        vfans[i]=f
                    elseif f==-1 || f==-2 || f==-3
                        ftypes[i]=-f
                    else
                        push!(fskips,i)
                    end
                end
                isempty(vfans) && (vfans=nothing)
                isempty(fskips) && (fskips=nothing)
                isempty(ftypes) && (ftypes=nothing)
            end
            # Whether this run emits any strip cell: one wall segment with
            # two non-empty side columns, or a wedge-emitting fan (a skip
            # fan contributes no cells of its own). A dead run leaves its
            # wire untouched.
            emits=fskips===nothing ? !isempty(fans) :
                any(i->!(i in fskips),fans)
            nl_off=length(offsets)
            if !emits
                for i in 1:(closed_run ? nchain : nchain-1)
                    na=voff===nothing ? nl_off : length(voff[i])
                    nb=voff===nothing ? nl_off : length(voff[i==nchain ? 1 : i+1])
                    min(na,nb)>0 && (emits=true;break)
                end
            end
            if !emits
                push!(dead_runs,ri)
                continue
            end
            if end_cols!==nothing
                for (ch,emit,v_end_v,adj,slot,dir) in end_cols
                    _bl_mark_end_covered!(ch,emit,covered_positions,
                                          claimed,meet_nodes,caller,t)
                    for pos in ch[1:emit]
                        wv=rvert(pos)
                        push!(snapped_coords,coords3d[wv])
                        covered_owner[coords3d[wv]]=Int32(adj)
                    end
                    top_override[v_end_v]=emit>0 ?
                        vertex_id!(rvert(ch[emit])) : vertex_id!(v_end_v)
                    if dir!==nothing
                        vdirs===nothing &&
                            (vdirs=Dict{Int,NTuple{3,Float64}}())
                        vdirs[slot]=dir
                    end
                end
            end
            # Upstream `treat2Connections`: at a wall vertex connected to
            # two wall segments the growth direction is the *biased*
            # averaged inward normal — `normalize(1.01·N1 + N2)` — scaled
            # by 1/cos(θ/2) where θ is the two normals' angle, so a corner
            # column grows proportionally deeper into the wedge. Straight
            # continuations reduce to the plain averaged normal (factor
            # 1); run-end vertices keep the covered-span direction
            # upstream's `addColumnAtTheEndOfTheBL` assigns; fans sweep
            # their own rays. `N1` is the normal of `_connections[0]` —
            # `getEdgesData` iterates the face's edges in tag order and
            # equal multimap keys keep insertion order, so the biased side
            # is the *lower-tagged* adjacent curve; interior to one curve
            # it is the earlier line index (the smaller-parameter side).
            fan_pos=Set(fans)
            nseg3=closed_run ? nchain : nchain-1
            segleft=Vector{NTuple{3,Float64}}(undef,max(nseg3,1))
            for i in 1:nseg3
                a3=coordsw[chain[i]]
                b3=coordsw[chain[i==nchain ? 1 : i+1]]
                s=(b3[1]-a3[1],b3[2]-a3[2],b3[3]-a3[3])
                sl=hypot(s[1],s[2],s[3])
                sl>0 || throw(ArgumentError(
                    "$caller: Surface[$t] wall curve has a coincident vertex"))
                ux,uy,uz=s[1]/sl,s[2]/sl,s[3]/sl
                segleft[i]=(plane_normal[2]*uz-plane_normal[3]*uy,
                            plane_normal[3]*ux-plane_normal[1]*uz,
                            plane_normal[1]*uy-plane_normal[2]*ux)
            end
            outcurve_raw=Int[loop_curve[mod1(start0+chain_positions[j]-1,
                                             nseg)] for j in 1:nseg3]
            outsign_raw=Int[loop_sign[mod1(start0+chain_positions[j]-1,
                                           nseg)] for j in 1:nseg3]
            seg_curve=Vector{Int}(undef,max(nseg3,1))
            seg_param_rev=falses(max(nseg3,1))
            for i in 1:nseg3
                sidx=chain_rev ?
                    (i<nchain ? nchain-i : nchain) : i
                seg_curve[i]=outcurve_raw[sidx]
                # Upstream iterates each edge's `lines` in parameter order
                # when it stores `MQuadrangle(v11,v21,v22,v12)`; a kernel
                # segment that walks param-decreasing (a loop-reversed
                # curve the orientation probe kept, or a loop-forward
                # curve the probe flipped) must evaluate its quads with
                # the column order swapped.
                kernel_param_positive=(outsign_raw[sidx]>0) != chain_rev
                seg_param_rev[i]=!kernel_param_positive
            end
            # Curve parameter of a chain vertex on a specific curve —
            # subdivision vertices carry it directly; model endpoints
            # resolve through the curve's point tags and param bounds.
            _p_on_curve=function (v::Int,c::Int)
                if get(vert_curve,v,0)==c && haskey(vert_param,v)
                    return vert_param[v]
                end
                pid=get(point_of,v,0)
                if pid!=0 && haskey(m.curves,c)
                    pa,pb=m.curves[c]
                    (pid==pa || pid==pb) || return NaN
                    t0,t1=_model_curve_param_bounds(m,c,caller)
                    return pid==pa ? t0 : t1
                end
                return NaN
            end
            vscale=nothing
            for i in 1:nchain
                interior=closed_run || (i>1 && i<nchain)
                interior || continue
                i in fan_pos && continue
                iin=closed_run ? mod1(i-1,nchain) : i-1
                nin=segleft[iin];nout=segleft[i]
                c_in=seg_curve[iin];c_out=seg_curve[i]
                if c_in==c_out
                    p_in=_p_on_curve(chain[iin],c_in)
                    p_out=_p_on_curve(chain[mod1(i+1,nchain)],c_out)
                    bias_in=(isnan(p_in)||isnan(p_out)) || p_in<p_out
                else
                    bias_in=c_in<c_out
                end
                n1=bias_in ? nin : nout
                n2=bias_in ? nout : nin
                x=(1.01*n1[1]+n2[1],1.01*n1[2]+n2[2],
                   1.01*n1[3]+n2[3])
                xl=hypot(x[1],x[2],x[3])
                xl>0 || continue
                d=abs(nin[1]*nout[1]+nin[2]*nout[2]+nin[3]*nout[3])
                d=min(d,1.0)
                ct2=cos(acos(d)/2)
                vdirs===nothing &&
                    (vdirs=Dict{Int,NTuple{3,Float64}}())
                vdirs[i]=(x[1]/xl,x[2]/xl,x[3]/xl)
                vscale===nothing && (vscale=Dict{Int,Float64}())
                vscale[i]=ct2>0 ? 1.0/ct2 : 1.0
            end
            # Tail fans measure `step` off the physical incident wall-edge
            # lengths — upstream reads `MEdge::length()`, an
            # `SVector3::norm` (`sqrt(Σ squares)`, not `hypot`); in the
            # ruled and mean-plane frames alike the kernel's projected
            # length is a working-frame distance, so pass the physical
            # mean through explicitly.
            tsteps=nothing
            if ftypes!==nothing && any(==(3),values(ftypes))
                tsteps=Dict{Int,Float64}()
                for (i,tp) in ftypes
                    tp==3 || continue
                    pa=coords3d[chain[i]]
                    tot=0.0;cnt=0
                    if closed_run || i>1
                        pb=coords3d[chain[i==1 ? nchain : i-1]]
                        dx=pa[1]-pb[1];dy=pa[2]-pb[2];dz=pa[3]-pb[3]
                        tot+=sqrt(dx*dx+dy*dy+dz*dz)
                        cnt+=1
                    end
                    if closed_run || i<nchain
                        pb=coords3d[chain[i==nchain ? 1 : i+1]]
                        dx=pa[1]-pb[1];dy=pa[2]-pb[2];dz=pa[3]-pb[3]
                        tot+=sqrt(dx*dx+dy*dy+dz*dz)
                        cnt+=1
                    end
                    tsteps[i]=cnt==0 ? 0.0 : tot/cnt
                end
            end
            strip=blmod._bl_2d_offsets(
                wall_mesh,offsets;fans=fans,
                fan_elements=spec.fan_elements,
                vertex_dirs=vdirs,
                vertex_offsets=voff,
                vertex_scale_last=vscale,
                vertex_fans=vfans,
                fan_skips=fskips,
                fan_types=ftypes,
                tail_steps=tsteps,
                seg_param_rev=seg_param_rev,
                plane_frame=pframe,
                plane_normal=plane_normal,caller=caller)
            mesh=strip.mesh
            # tops keys on the wall node's exact 3-D coordinate. Run-end
            # tops are the covered-span top vertices the walk resolved —
            # snap the strip's computed top onto that exact coordinate.
            for (i,v) in enumerate(chain)
                vid=get(top_override,v,0)
                if vid!=0
                    push!(snapped_coords,new_coords[vid])
                    top_of[v]=vid
                else
                    top=strip.tops[coordsw[v]]
                    top_of[v]=new_vertex!(top[axw[1]],top[axw[2]],lift(top))
                end
            end
            if variable || !closed_run
                # The strip's exposed boundary through each wall vertex:
                # the top free edge enters at the deepest level the
                # preceding segment used; the uncancelled column rungs walk
                # to the level the next segment exits at (upstream's
                # `bedges` symmetric difference). A zero-depth neighbour
                # emits no cells on the segment — the free side descends
                # the whole column to the wall vertex. A fan vertex
                # instead runs up its entry ray, across the level-n ring,
                # and down the exit ray. Open-run end vertices take the
                # covered-chain top directly (`top_override`).
                fan_set_chain=Set(fans)
                for (i,v) in enumerate(chain)
                    !closed_run && (i==1 || i==nchain) && continue
                    iprev=i==1 ? nchain : i-1
                    inext=i==nchain ? 1 : i+1
                    ni=strip.nlevels[i]
                    nb_p=strip.nlevels[iprev];nb_n=strip.nlevels[inext]
                    e_i=(ni==0 || nb_p==0) ? 0 : (ni>nb_p ? nb_p+1 : ni)
                    x_i=(ni==0 || nb_n==0) ? 0 : (ni>nb_n ? nb_n+1 : ni)
                    entry=chain_rev ? x_i : e_i
                    last_level=chain_rev ? e_i : x_i
                    ids=Int[]
                    vn=strip.verts[i]
                    if i in fan_set_chain
                        fanw=strip.fan_wedges[i]
                        rs=chain_rev ? fanw : 0
                        re=chain_rev ? 0 : fanw
                        if fskips!==nothing && i in fskips
                            # Skip fan: no wedge cells cover the sector —
                            # upstream's `bedges` keep every one-sided rung
                            # free, so the boundary walks down the entry
                            # ray to the wall vertex and back up the exit
                            # ray (the fill then touches the vertex itself).
                            for k in entry:-1:1
                                push!(ids,strip.id_of[(vn,k,rs)])
                            end
                            push!(ids,strip.id_of[(vn,0,0)])
                            for k in 1:last_level
                                push!(ids,strip.id_of[(vn,k,re)])
                            end
                        elseif ftypes!==nothing && haskey(ftypes,i)
                            # Cross/mix/tail: the kernel's `fan_free` chain
                            # threads the interior column tops from ray 0 to
                            # ray `fanw` (kernel order — reverse it when the
                            # loop walks the chain backwards); the entry
                            # boundary column climbs to its top and the exit
                            # side descends to the pair depth like upstream's
                            # `bedges` keep the uncovered rungs free.
                            fc=strip.fan_free[vn]
                            ordered=chain_rev ? reverse(fc) : fc
                            for k in entry:ni-1
                                push!(ids,strip.id_of[(vn,k,k==0 ? 0 : rs)])
                            end
                            append!(ids,ordered)
                            for k in (ni-1):-1:last_level
                                push!(ids,strip.id_of[(vn,k,k==0 ? 0 : re)])
                            end
                        else
                            for k in entry:ni
                                push!(ids,strip.id_of[(vn,k,k==0 ? 0 : rs)])
                            end
                            dr=re>rs ? 1 : -1
                            for r in (rs+dr):dr:re
                                push!(ids,strip.id_of[(vn,ni,r)])
                            end
                            for k in (ni-1):-1:last_level
                                push!(ids,strip.id_of[(vn,k,k==0 ? 0 : re)])
                            end
                        end
                    else
                        dk=last_level>=entry ? 1 : -1
                        for k in entry:dk:last_level
                            push!(ids,strip.id_of[(vn,k,0)])
                        end
                    end
                    nodes=Int[]
                    for nid in ids
                        c=(mesh.coords[1,nid],mesh.coords[2,nid],
                           mesh.coords[3,nid])
                        push!(nodes,new_vertex!(c[axw[1]],c[axw[2]],lift(c)))
                    end
                    wire_nodes[v]=nodes
                end
            end
            base=length(strip_coords)
            for i in axes(mesh.coords,2)
                kc=(mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])
                push!(strip_coords,lift(kc))
                # The `(2,t)` params are upstream's `parFromPoint` working
                # coordinates — the kernel's mean-plane pair in `pframe`
                # mode; the dropped-plane/`(u,v)` working pair otherwise.
                push!(strip_params,pframe===nothing ?
                    [kc[axw[1]],kc[axw[2]]] :
                    collect(strip.work_points[i]))
            end
            for block in mesh.blocks
                if block.msh==2
                    for c in axes(block.nodes,2)
                        push!(strip_tris,(base+Int(block.nodes[1,c]),
                            base+Int(block.nodes[2,c]),
                            base+Int(block.nodes[3,c])))
                    end
                elseif block.msh==3
                    for c in axes(block.nodes,2)
                        push!(strip_quads,(base+Int(block.nodes[1,c]),
                            base+Int(block.nodes[2,c]),
                            base+Int(block.nodes[3,c]),
                            base+Int(block.nodes[4,c])))
                        # `quads_eval` columns line up 1:1 with the emitted
                        # quad block — upstream's stored vertex order for
                        # the `quadsToTriangles` split comparisons.
                        push!(strip_quads_eval,(
                            base+Int(strip.quads_eval[1,c]),
                            base+Int(strip.quads_eval[2,c]),
                            base+Int(strip.quads_eval[3,c]),
                            base+Int(strip.quads_eval[4,c])))
                    end
                else
                    throw(ErrorException(
                        "$caller: boundary-layer strip emitted unexpected " *
                        "cell type $(block.msh)"))
                end
            end
        end
        # Rebuild the wire: walk the rotated order; at each run emit only
        # the top polyline. Run-end vertices are covered by the strip
        # side (their columns ride the adjacent segments), so they drop
        # out of the fill boundary — the top of each end column lands
        # on the adjacent segment and the wire resumes past it.
        new_loop=Int[]
        k=1
        while k<=nseg
            sk=mod1(start0+k-1,nseg)
            run_index=0
            for (ri,run) in enumerate(runs)
                if run[2]<=k<=run[3]
                    run_index=ri
                    break
                end
            end
            if run_index==0
                k in covered_positions ||
                    push!(new_loop,vertex_id!(loop[sk]))
                k+=1
                continue
            end
            si,first,last,closed_flag=runs[run_index]
            if !iszero(closed_flag)
                for p in 1:nseg
                    v=loop[mod1(start0+p-1,nseg)]
                    wn=get(wire_nodes,v,nothing)
                    if wn===nothing
                        push!(new_loop,top_of[v])
                    else
                        append!(new_loop,wn)
                    end
                end
                break
            end
            for p in first:last+1
                v=loop[mod1(start0+p-1,nseg)]
                if run_index in dead_runs
                    push!(new_loop,vertex_id!(v))
                    continue
                end
                wn=get(wire_nodes,v,nothing)
                if wn===nothing
                    push!(new_loop,top_of[v])
                else
                    append!(new_loop,wn)
                end
            end
            # The run's end vertex sits at rotated position `last+1`,
            # consumed by the strip side — resume on the next segment.
            k=last+2
        end
        # Two run ends meeting on one adjacent edge both resolve their top
        # to the shared chain node — drop the consecutive duplicate.
        if !isempty(new_loop)
            dedup=Int[]
            for vid in new_loop
                (!isempty(dedup) && dedup[end]==vid) || push!(dedup,vid)
            end
            length(dedup)>1 && dedup[1]==dedup[end] && pop!(dedup)
            new_loop=dedup
        end
        push!(reduced_wires,new_loop)
    end
    segs2=Tuple{Int,Int}[]
    for wire in reduced_wires
        nloop=length(wire)
        nloop>=3 || throw(ArgumentError(
            "$caller: Surface[$t] boundary layer leaves a degenerate wire"))
        for k in 1:nloop
            i,j=wire[k],wire[mod1(k+1,nloop)]
            i==j || (new_xs[i]==new_xs[j] && new_ys[i]==new_ys[j]) ||
                push!(segs2,(i,j))
        end
    end
    # Embedded entities ride the same constraints the ordinary fill carries —
    # remap their PSLG ids through `vertex_id!`.
    embedded_vertex_remap=Dict{Int,Int}()
    internal2=Tuple{Int,Int}[]
    for (a,b) in internal
        na=get!(embedded_vertex_remap,a) do
            vertex_id!(a)
        end
        nb=get!(embedded_vertex_remap,b) do
            vertex_id!(b)
        end
        na==nb && continue
        push!(internal2,(na,nb))
    end
    for (edim,etag) in embedded
        edim==0 || continue
        haskey(m.points,etag) || throw(ArgumentError(
            "$caller: unknown embedded Point[$etag]"))
        # `_bl_surface_pslg` already added the point vertex; find it through
        # `index` so it lands in the fill's vertex set.
        vertex_id!(index[etag])
    end
    triangulation=try
        constrained_delaunay(
            new_xs,new_ys,segs2;internal_segments=internal2)
    catch err
        # A strip whose free boundary self-intersects — a convex-corner fan
        # reaching onto the neighbour wall, or an over-thick layer whose
        # opposing strips meet — cannot be triangulated; upstream's
        # `modifyInitialMeshForBoundaryLayers` fails the same configuration
        # with "edge not recovered", so surface the cause plainly.
        if err isa ArgumentError &&
           occursin("crosses existing constraint",sprint(showerror,err))
            throw(ArgumentError(
                "$caller: Surface[$t] boundary-layer strip free boundary " *
                "self-intersects — upstream fails this configuration " *
                "identically (edge not recovered)"))
        end
        rethrow()
    end
    fill=to_mesh(
        triangulation;interior=classify_interior(triangulation))
    # Merge every vertex class into one coordinate pool: boundary and top
    # vertices first (stable order), then strip interiors.
    pool=Dict{NTuple{3,Float64},Int32}()
    out_coords=NTuple{3,Float64}[]
    out_owners=Tuple{Int,Int32}[]
    out_params=Union{Nothing,Vector{Float64}}[]
    function pool_id!(coordinate,owner,param)
        return get!(pool,coordinate) do
            push!(out_coords,coordinate)
            push!(out_owners,owner);push!(out_params,param)
            Int32(length(out_coords))
        end
    end
    new_to_pool=Vector{Int32}(undef,length(new_coords))
    for i in eachindex(new_coords)
        new_to_pool[i]=pool_id!(new_coords[i],new_owners[i],new_params[i])
    end
    fill_to_pool=Vector{Int32}(undef,nnodes(fill))
    for node in 1:nnodes(fill)
        key=(fill.coords[1,node]==0 ? 0.0 : fill.coords[1,node],
             fill.coords[2,node]==0 ? 0.0 : fill.coords[2,node])
        vid=get(vertex_index,key,nothing)
        vid===nothing && throw(ErrorException(
            "$caller: Surface[$t] fill introduced an unmapped coordinate"))
        fill_to_pool[node]=new_to_pool[vid]
    end
    snap_tol=1e-9*max(1.0,maximum(xs)-minimum(xs),
                      maximum(ys)-minimum(ys))
    strip_to_pool=Vector{Int32}(undef,length(strip_coords))
    for i in eachindex(strip_coords)
        coordinate=strip_coords[i]
        for s in snapped_coords
            (abs(coordinate[1]-s[1])<=snap_tol &&
             abs(coordinate[2]-s[2])<=snap_tol &&
             abs(coordinate[3]-s[3])<=snap_tol) &&
                (coordinate=s;break)
        end
        ctag=get(covered_owner,coordinate,Int32(0))
        strip_to_pool[i]=ctag==0 ?
            pool_id!(coordinate,(2,Int32(t)),strip_params[i]) :
            pool_id!(coordinate,(1,ctag),nothing)
    end
    out_tris=NTuple{3,Int32}[]
    out_quads=NTuple{4,Int32}[]
    for tri in strip_tris
        push!(out_tris,(strip_to_pool[tri[1]],strip_to_pool[tri[2]],
                        strip_to_pool[tri[3]]))
    end
    keep_quads=_model_surface_recombined(m,t) ||
        any(spec->spec.quads,specs)
    for (qi,quad) in enumerate(strip_quads)
        a,b,c,d=quad
        if keep_quads
            push!(out_quads,(strip_to_pool[a],strip_to_pool[b],
                             strip_to_pool[c],strip_to_pool[d]))
        else
            # Upstream `quadsToTriangles` (meshGFaceOptimize.cpp) splits
            # each emitted quad on whichever diagonal gives the better
            # worst-triangle `gammaShapeMeasure`, evaluated on the quad's
            # *stored* vertex order — `t11=(v0,v1,v2)`, `t12=(v2,v3,v0)`
            # (option 1, the `(v0,v2)` diagonal) against `t21=(v1,v2,v3)`,
            # `t22=(v3,v0,v1)` (option 2, the `(v1,v3)` diagonal). Option 1
            # wins on `qual1>qual2` and is forced when option 2's
            # triangles face opposite ways (a non-convex quad, `ori2<0`);
            # an exact tie takes option 2. `strip_quads_eval` holds
            # upstream's stored tuple so the fp-level comparisons —
            # including near-ties — reproduce upstream's, while `quad`
            # keeps the positive ring order for emission.
            e1,e2,e3,e4=strip_quads_eval[qi]
            p1,p2,p3,p4=strip_coords[e1],strip_coords[e2],
                        strip_coords[e3],strip_coords[e4]
            qual1=min(_bl_tri_gamma(p1,p2,p3),_bl_tri_gamma(p3,p4,p1))
            qual2=min(_bl_tri_gamma(p2,p3,p4),_bl_tri_gamma(p4,p1,p2))
            n1x,n1y,n1z=_bl_tri_normal(p2,p3,p4)
            n2x,n2y,n2z=_bl_tri_normal(p4,p1,p2)
            ori2=n1x*n2x+n1y*n2y+n1z*n2z
            d1,d2=if qual1>qual2 || ori2<0
                (e1,e3)
            else
                (e2,e4)
            end
            if (d1==a && d2==c) || (d1==c && d2==a)
                # (a,c) diagonal — positive halves (a,b,c) and (a,c,d).
                push!(out_tris,(strip_to_pool[a],strip_to_pool[b],
                                strip_to_pool[c]))
                push!(out_tris,(strip_to_pool[a],strip_to_pool[c],
                                strip_to_pool[d]))
            else
                ((d1==b && d2==d) || (d1==d && d2==b)) || throw(
                    ErrorException(
                        "$caller: quad eval tuple is not a diagonal of "*
                        "the stored quadrangle"))
                # (b,d) diagonal — positive halves (a,b,d) and (b,c,d).
                push!(out_tris,(strip_to_pool[a],strip_to_pool[b],
                                strip_to_pool[d]))
                push!(out_tris,(strip_to_pool[b],strip_to_pool[c],
                                strip_to_pool[d]))
            end
        end
    end
    for c in 1:ntris(fill)
        push!(out_tris,(fill_to_pool[Int(fill.tris[1,c])],
                        fill_to_pool[Int(fill.tris[2,c])],
                        fill_to_pool[Int(fill.tris[3,c])]))
    end
    coords_out=Matrix{Float64}(undef,3,length(out_coords))
    for i in eachindex(out_coords)
        coords_out[1,i],coords_out[2,i],coords_out[3,i]=out_coords[i]
    end
    isempty(out_tris) && isempty(out_quads) && throw(ErrorException(
        "$caller: Surface[$t] boundary-layer pass produced no cells"))
    result=if isempty(out_quads)
        tris=Matrix{Int32}(undef,3,length(out_tris))
        for (i,tri) in enumerate(out_tris)
            tris[1,i],tris[2,i],tris[3,i]=tri
        end
        Mesh(coords_out;tris=tris)
    else
        tris=Matrix{Int32}(undef,3,length(out_tris))
        for (i,tri) in enumerate(out_tris)
            tris[1,i],tris[2,i],tris[3,i]=tri
        end
        quads=Matrix{Int32}(undef,4,length(out_quads))
        for (i,quad) in enumerate(out_quads)
            quads[1,i],quads[2,i],quads[3,i],quads[4,i]=quad
        end
        blocks=ElementBlock[]
        isempty(out_tris) || push!(blocks,ElementBlock(2,tris))
        push!(blocks,ElementBlock(3,quads))
        MixedMesh(coords_out,blocks)
    end
    diagnostic=validate(result)
    diagnostic.ok || throw(ErrorException(
        "$caller: Surface[$t] boundary-layer mesh is invalid — " *
        join(diagnostic.messages,"; ")))
    return result,out_owners,out_params
end

# `_occ_curved_surface_production` for `Mesh|MixedMesh` strip output —
# boundary vertices reuse the consumed lower product's actual UIDs
# (matched by owning entity and bitwise coordinate, like upstream's
# `meshGFace` vertex reassociation) while strip/fill interiors mint fresh
# labels from the authority's high water.
function _bl_surface_production(m::GeoModel,t::Int,mesh,
                                owners::Vector{Tuple{Int,Int32}},
                                params,input,limit,caller::AbstractString)
    tags=Vector{UInt64}(undef,nnodes(mesh))
    retained=falses(nnodes(mesh))
    lookup=Dict{Tuple{Int,Int32,NTuple{3,Float64}},UInt64}()
    source=input===nothing ? nothing : input.mesh
    if source!==nothing && source.entity_data!==nothing
        data=source.entity_data
        for row in 1:nnodes(source)
            entity=data.node_entities[row]
            key=(entity[1],entity[2],
                 (source.coords[1,row],source.coords[2,row],
                  source.coords[3,row]))
            haskey(lookup,key) || (lookup[key]=data.external_node_tags[row])
        end
    end
    next=UInt64(limit)
    for k in 1:nnodes(mesh)
        owner=owners[k]
        point=(mesh.coords[1,k],mesh.coords[2,k],mesh.coords[3,k])
        tag=get(lookup,(owner[1],owner[2],point),UInt64(0))
        if tag!=0
            tags[k]=tag;retained[k]=true
        else
            next+=UInt64(1)
            iszero(next) && throw(ArgumentError(
                "$caller: Surface[$t] vertex births exhaust UInt64 labels"))
            tags[k]=next
        end
    end
    return _ModelSurfaceProduction(mesh,tags,owners,params,retained,next)
end
