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
        return (occ=:copied_curve,
                curve=_model_occ_nurbs_line(geometry,caller),
                t0=geometry.t0,t1=geometry.t1)
    elseif occ in (:circle,:ellipse)
        return (occ=:copied_curve,
                curve=_model_occ_nurbs_conic(geometry,caller),
                t0=geometry.t0,t1=geometry.t1)
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
        dim in (-2,2,3) && throw(ArgumentError(
            "$caller: surface and solid NURBS conversion is not implemented"))
        push!(roots,id)
    end
    # Convert each edge record in the transitive shape closure once; the
    # public parametrization is preserved so shared sub-shapes stay
    # consistent for every referencing face.
    staged=Dict{UInt64,NamedTuple}()
    stack=copy(roots)
    seen=Set{UInt64}()
    while !isempty(stack)
        id=pop!(stack)
        id in seen && continue
        push!(seen,id)
        shape=kernel.shapes[id]
        append!(stack,shape.children)
        shape.dimension==1 || continue
        record=shape.geometry
        converted=_model_occ_nurbs_curve_record(record.geometry,caller)
        converted===record.geometry && continue
        staged[id]=merge(record,(kind=:tensor_bspline,geometry=converted))
    end
    for (id,record) in staged
        shape=kernel.shapes[id]
        kernel.shapes[id]=_ModelOCCShape(1,record,shape.children)
        kernel.changed=true
        tag=kernel.shape_tags[id][2]
        # The rebound edge is re-imported on the next synchronize — drop the
        # previously materialized model entries like the corner ops do.
        delete!(m.curves,tag);delete!(m.curve_types,tag)
        delete!(m.curve_geometry,tag)
    end
    return nothing
end
