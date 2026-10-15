# `OCC_Internals::convertToNURBS` — `BRepBuilderAPI_NurbsConvert` applied per
# bound entity: geometry records are replaced by exact NURBS forms under the
# same tags and topology. Lines become degree-1 B-splines on the identical
# (arc-length) parametrization; circle and ellipse edges become quadratic-
# rational B-splines carrying `_ModelOCCConicReparam`, so the public
# parameter map — and every stored pcurve — is preserved bitwise (upstream
# refits 2-D curves approximately instead). Spline edges are already NURBS;
# vertices, wires, degenerate and `surface_iso` edges pass through. Surface
# and solid/shell patch conversion is a follow-up half of this slice — dim
# -2, 2 and 3 inputs reject loudly for now rather than half-converting.

# Rational-quadratic image of a circle/ellipse arc: equal spans of at most
# π/2, interior poles on the tangent intersections `C + sec(Δ/2)·(p(m)−C)`,
# weights cos(Δ/2). `frame` resolves the eval-consistent record so
# native-backed circles convert through the frame evaluators actually use.
function _model_occ_nurbs_conic(geometry,caller)
    frame=_occ_2df_frame(geometry)
    C=frame.center
    eval=geometry.occ===:circle ? _occ_circle_point : _occ_ellipse_point
    t0,t1=geometry.t0,geometry.t1
    span=t1-t0
    span>0.0 || throw(ArgumentError("$caller: degenerate conic range"))
    n=max(1,ceil(Int,span/(π/2)))
    δ=span/n
    quarter=_gm_tan(δ/4)
    weight=_gm_cos(δ/2)
    secant=1/weight
    poles=NTuple{3,Float64}[eval(geometry,t0)]
    weights=Float64[1.0]
    tanq=Float64[]
    for i in 1:n
        a=t0+(i-1)*δ
        b=i==n ? t1 : t0+i*δ
        m=(a+b)/2
        pm=eval(geometry,m)
        push!(tanq,quarter)
        push!(poles,ntuple(k->fma(pm[k]-C[k],secant,C[k]),3))
        push!(weights,weight)
        push!(poles,eval(geometry,b))
        push!(weights,1.0)
    end
    knots=Float64[t0,t0,t0]
    for i in 1:n-1
        knot=t0+i*δ
        push!(knots,knot,knot)
    end
    append!(knots,(t1,t1,t1))
    curve=NURBS.NURBSCurve(2,knots,poles,weights)
    inner=_ModelOCCNURBSCurve(:bspline,curve)
    reparam=_ModelOCCConicReparam(Float64[t0+i*δ for i in 0:n],tanq)
    return _ModelOCCReparamCurve(:bspline,curve,inner.nets,reparam)
end

# Degree-1 B-spline on the identical arc-length parametrization — no
# reparametrization wrapper required.
function _model_occ_nurbs_line(geometry,caller)
    t0,t1=geometry.t0,geometry.t1
    t1>t0 || throw(ArgumentError("$caller: degenerate line range"))
    return _model_occ_curve_geometry(:bspline,
        NTuple{3,Float64}[_occ_line_point(geometry,t0),
                          _occ_line_point(geometry,t1)];
        degree=1,knots=[t0,t1],multiplicities=[2,2],caller=caller)
end

function _model_occ_nurbs_curve_record(geometry,caller)
    occ=get(geometry,:occ,nothing)
    if occ===:line
        # Analytic fields stay as inert metadata: eval dispatches on
        # `occ===:copied_curve` → `.curve`, while materialized-primitive
        # consistency checks can still read center/radius parameters.
        return merge(geometry,(occ=:copied_curve,converted_from=:line,
                curve=_model_occ_nurbs_line(geometry,caller)))
    elseif occ in (:circle,:ellipse)
        return merge(geometry,(occ=:copied_curve,converted_from=occ,
                curve=_model_occ_nurbs_conic(geometry,caller)))
    elseif occ===:copied_curve
        curve=geometry.curve
        curve isa _ModelOCCReparamCurve && return geometry
        # `GeomConvert` turns a Bezier into a single-span B-spline.
        curve.family===:bezier ||
            return geometry
        return merge(geometry,(curve=_ModelOCCNURBSCurve(:bspline,curve.curve),))
    elseif occ in (:degenerate,:surface_iso)
        # No standalone 3-D curve — the degenerate edge evaluates through its
        # owner face's pcurve and the iso edge through the surface record.
        return geometry
    end
    throw(ArgumentError(
        "$caller: OCC curve kind $(occ) cannot be converted to NURBS"))
end

function _model_occ_convert_to_nurbs!(m::GeoModel,dim_tags,caller)
    kernel=m.occ_kernel
    # Validate and stage everything before the first kernel mutation —
    # failures leave the registry byte-identical (a superset of upstream,
    # which converts greedily per entity).
    roots=UInt64[]
    for dt in dim_tags
        length(dt)==2 || throw(ArgumentError(
            "$caller: dimTag pairs must have two entries"))
        dim=_model_occ_dimension(dt[1],caller)
        tag=_model_occ_requested_tag(dt[2],caller)
        id=get(kernel.bindings,(dim,tag),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $dim with tag $tag"))
        push!(roots,id)
    end
    # Convert each edge/face record in the transitive shape closure once;
    # the public parametrization is preserved so shared sub-shapes stay
    # consistent for every referencing face.
    staged=Dict{UInt64,NamedTuple}()
    stack=copy(roots)
    seen=Set{UInt64}()
    order=UInt64[]
    while !isempty(stack)
        id=pop!(stack)
        id in seen && continue
        push!(seen,id)
        push!(order,id)
        shape=kernel.shapes[id]
        append!(stack,shape.children)
        if shape.dimension==1
            record=shape.geometry
            converted=_model_occ_nurbs_curve_record(record.geometry,caller)
            converted===record.geometry && continue
            staged[id]=merge(record,(kind=:tensor_bspline,geometry=converted))
        elseif shape.dimension==2
            record=shape.geometry
            converted=_model_occ_nurbs_surface_record(record.geometry,caller)
            converted===record.geometry && continue
            staged[id]=merge(record,(geometry=converted,))
        end
    end
    # Member degenerate/surface-iso edges embed their owner face record —
    # rebind each to the converted record so the pcurve lookup and the
    # surface evaluation stay on the same object.
    for id in order
        shape=kernel.shapes[id]
        shape.dimension==2 && haskey(staged,id) || continue
        newface=staged[id].geometry
        for wire_id in shape.children
            wire=kernel.shapes[wire_id]
            haskey(wire.geometry,:uses) || continue
            for (eid,_) in wire.geometry.uses
                eshape=kernel.shapes[eid]
                rec=eshape.geometry.geometry
                hasproperty(rec,:face) || continue
                entry=get(staged,eid,eshape.geometry)
                staged[eid]=merge(entry,(geometry=merge(
                    entry.geometry,(face=newface,)),))
            end
        end
    end
    for (id,record) in staged
        shape=kernel.shapes[id]
        kernel.shapes[id]=_ModelOCCShape(shape.dimension,record,shape.children)
        kernel.changed=true
        tag=get(kernel.shape_tags,id,nothing)
        tag===nothing && continue
        t=tag[2]
        if shape.dimension==1
            # The rebound edge is re-imported on the next synchronize — drop
            # the previously materialized model entries like the corner ops.
            delete!(m.curves,t);delete!(m.curve_types,t)
            delete!(m.curve_geometry,t)
        elseif shape.dimension==2
            delete!(m.surfaces,t);delete!(m.surface_types,t)
            delete!(m.surface_geometry,t)
        end
    end
    return nothing
end

# ── surface conversion ─────────────────────────────────────────────────────
# A converted face keeps its `occ` kind as analytic metadata and gains
# `patch`+`reparam`; every evaluator checks `_model_nurbs_surface` first, so
# the NURBS is authoritative while the preserved external parametrization
# keeps pcurves, u-periods, measures and materializers bitwise consistent.

# One conic axis' pole layout over [t0,t1]: per-pole parameter value, radial
# factor (1 or sec(Δ/2)) and weight — shared by every quadratic-rational
# surface builder.
function _occ_nurbs_conic_axis(t0,t1)
    span=t1-t0
    n=max(1,ceil(Int,span/(π/2)))
    δ=span/n
    secant=1/_gm_cos(δ/2); weight=_gm_cos(δ/2)
    vals=Float64[];qs=Float64[];ws=Float64[]
    for i in 0:2n
        if iseven(i)
            push!(vals,t0+(i>>1)*δ);push!(qs,1.0);push!(ws,1.0)
        else
            push!(vals,t0+(i>>1)*δ+δ/2);push!(qs,secant);push!(ws,weight)
        end
    end
    knots=Float64[t0,t0,t0]
    for i in 1:n-1
        knot=t0+i*δ
        push!(knots,knot,knot)
    end
    append!(knots,(t1,t1,t1))
    reparam=_ModelOCCConicReparam(Float64[t0+i*δ for i in 0:n],
        fill(_gm_tan(δ/4),n))
    return vals,qs,ws,knots,reparam
end

function _model_occ_nurbs_plane(g)
    lo,hi=g.uvb
    U=Float64[lo[1],lo[1],hi[1],hi[1]]
    V=Float64[lo[2],lo[2],hi[2],hi[2]]
    C,X,Y=g.center,g.X,_occ_frame_y(g)
    us=(lo[1],hi[1]);vs=(lo[2],hi[2])
    controls=Matrix{NTuple{3,Float64}}(undef,2,2)
    for i in 1:2,j in 1:2
        controls[i,j]=ntuple(k->fma(us[i],X[k],fma(vs[j],Y[k],C[k])),3)
    end
    patch=_ModelNURBSSurface(:bspline,
        NURBS.NURBSSurface(1,1,U,V,controls,ones(2,2)))
    return patch,nothing
end

function _model_occ_nurbs_cylinder(g)
    angle=get(g,:angle,2π)
    vals,qs,ws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    V=Float64[0.0,0.0,g.height,g.height]
    C,X,Y,Z=g.center,g.X,_occ_frame_y(g),g.axis
    r=g.radius
    nu=length(vals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,2)
    W=Matrix{Float64}(undef,nu,2)
    for i in 1:nu,j in 1:2
        v=j==1 ? 0.0 : g.height
        s,c=_gm_sincos(vals[i])
        controls[i,j]=ntuple(k->fma(qs[i]*r*c,X[k],
            fma(qs[i]*r*s,Y[k],fma(v,Z[k],C[k]))),3)
        W[i,j]=ws[i]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,1,U,V,controls,W))
    return patch,(u=urep,v=nothing)
end

function _model_occ_nurbs_cone(g)
    angle=get(g,:angle,2π)
    vals,qs,ws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    slant=sqrt(g.height*g.height+(g.r1-g.r2)*(g.r1-g.r2))
    V=Float64[0.0,0.0,slant,slant]
    sa=_gm_atan((g.r2-g.r1)/g.height);ssa,csa=_gm_sincos(sa)
    C,X,Y,Z=g.center,g.X,_occ_frame_y(g),g.axis
    nu=length(vals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,2)
    W=Matrix{Float64}(undef,nu,2)
    for i in 1:nu,j in 1:2
        v=j==1 ? 0.0 : slant
        R=fma(v,ssa,g.r1);A3=v*csa
        s,c=_gm_sincos(vals[i])
        controls[i,j]=ntuple(k->fma(qs[i]*R*c,X[k],
            fma(qs[i]*R*s,Y[k],fma(A3,Z[k],C[k]))),3)
        W[i,j]=ws[i]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,1,U,V,controls,W))
    return patch,(u=urep,v=nothing)
end

function _model_occ_nurbs_sphere(g)
    angle=get(g,:angle,2π)
    lat1,lat2=get(g,:lat1,-π/2),get(g,:lat2,π/2)
    uvals,uqs,uws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    vvals,vqs,vws,V,vrep=_occ_nurbs_conic_axis(lat1,lat2)
    C,X,Y,Z=g.center,g.X,g.Y,g.axis
    r=g.radius
    nu,nv=length(uvals),length(vvals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    W=Matrix{Float64}(undef,nu,nv)
    for i in 1:nu,j in 1:nv
        su,cu=_gm_sincos(uvals[i]);sv,cv=_gm_sincos(vvals[j])
        # Row j is a circle of radius r·qv_j·cos(v̂_j) centered at
        # C + r·qv_j·sin(v̂_j)·Z — the mid-v factor applies to both
        # components, the mid-u factor only to the radial direction.
        radial=r*vqs[j]*cv*uqs[i]
        controls[i,j]=ntuple(k->fma(radial*cu,X[k],
            fma(radial*su,Y[k],fma(r*vqs[j]*sv,Z[k],C[k]))),3)
        W[i,j]=uws[i]*vws[j]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,2,U,V,controls,W))
    return patch,(u=urep,v=vrep)
end

function _model_occ_nurbs_torus(g)
    angle=get(g,:angle,2π)
    uvals,uqs,uws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    vvals,vqs,vws,V,vrep=_occ_nurbs_conic_axis(0.0,2π)
    C,X,Y,Z=g.center,g.X,g.Y,g.axis
    r1,r2=g.r1,g.r2
    nu,nv=length(uvals),length(vvals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    W=Matrix{Float64}(undef,nu,nv)
    for i in 1:nu,j in 1:nv
        su,cu=_gm_sincos(uvals[i]);sv,cv=_gm_sincos(vvals[j])
        rho=uqs[i]*(r1+r2*vqs[j]*cv)
        A3=r2*vqs[j]*sv
        controls[i,j]=ntuple(k->fma(rho*cu,X[k],
            fma(rho*su,Y[k],fma(A3,Z[k],C[k]))),3)
        W[i,j]=uws[i]*vws[j]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,2,U,V,controls,W))
    return patch,(u=urep,v=vrep)
end

function _model_occ_nurbs_surface_record(geometry,caller)
    occ=get(geometry,:occ,nothing)
    if hasproperty(geometry,:patch)
        # Already-NURBS tensor record: a bezier patch renormalizes to a
        # B-spline like `GeomConvert`; anything else passes through.
        geometry.patch.family===:bezier || return geometry
        return merge(geometry,(occ=:tensor_bspline,
            patch=_ModelNURBSSurface(:bspline,geometry.patch.surface),))
    end
    patch,reparam=if occ===:plane
        _model_occ_nurbs_plane(geometry)
    elseif occ===:cylinder
        _model_occ_nurbs_cylinder(geometry)
    elseif occ===:cone
        _model_occ_nurbs_cone(geometry)
    elseif occ===:sphere
        _model_occ_nurbs_sphere(geometry)
    elseif occ===:torus
        _model_occ_nurbs_torus(geometry)
    else
        throw(ArgumentError(
            "$caller: OCC surface kind $(occ) cannot be converted to NURBS"))
    end
    return merge(geometry,(patch=patch,reparam=reparam))
end
