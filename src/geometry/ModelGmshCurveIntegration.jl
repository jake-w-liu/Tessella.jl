# Built-in GEO Curve integration follows meshGEdge.cpp. The public standalone
# Mesh1D kernel keeps its robust callable-curve arithmetic; this path preserves
# the native derivative, primitive, smoothing and NUMP expression order.
@inline function _gmsh_curve_trapezoid(first, last, caller)
    value=0.5*(first.lc+last.lc)*(last.t-first.t)
    isfinite(value) || throw(ArgumentError("$caller: native Curve integral overflowed"))
    return value
end

function _gmsh_curve_recursive_integration!(points, first, last, evaluate,
                                            precision, depth, caller)
    t=0.5*(first.t+last.t)
    isfinite(t) || (t=_convex_coordinate(first.t,last.t,0.5))
    first.t<t<last.t || throw(ArgumentError(
        "$caller: native Curve integration cannot refine its parameter interval"))
    middle=evaluate(t,false,false)
    coarse=_gmsh_curve_trapezoid(first,last,caller)
    left=_gmsh_curve_trapezoid(first,middle,caller)
    right=_gmsh_curve_trapezoid(middle,last,caller)
    error=abs(coarse-left-right)
    if (error<precision && depth>6) || depth>25
        length(points)<=_gmsh_curve_point_capacity(points)-2 || throw(ArgumentError(
            "$caller: native Curve integration exceeded its point capacity"))
        middle_p=points[end].p+left
        last_p=middle_p+right
        isfinite(middle_p) && isfinite(last_p) || throw(ArgumentError(
            "$caller: native Curve primitive overflowed"))
        push!(points,_IntegrationPoint(middle.t,middle.lc,middle_p,middle.xp,middle.h))
        push!(points,_IntegrationPoint(last.t,last.lc,last_p,last.xp,last.h))
    else
        _gmsh_curve_recursive_integration!(points,first,middle,evaluate,
                                           precision,depth+1,caller)
        _gmsh_curve_recursive_integration!(points,middle,last,evaluate,
                                           precision,depth+1,caller)
    end
    return nothing
end

function _gmsh_curve_integration(evaluate,t0::Float64,t1::Float64,
                                 precision::Float64,caller)
    isfinite(precision) && precision>0 || throw(ArgumentError(
        "$caller: native Curve integration precision must be positive and finite"))
    first=evaluate(t0,true,false)
    last=evaluate(t1,false,true)
    points=_gmsh_curve_point_store()
    push!(points,first)
    try
        _gmsh_curve_recursive_integration!(points,first,last,evaluate,precision,1,caller)
    catch
        _gmsh_curve_close_points!(points)
        rethrow()
    end
    return points
end

function _gmsh_curve_smooth_primitive!(points,caller)
    alpha=sqrt(1.8)
    @inbounds for _ in 1:2001
        changed=false
        for i in 2:length(points)
            current,previous=points[i],points[i-1]
            dh=current.xp/current.lc-previous.xp/previous.lc
            dt=current.t-previous.t
            if dh/dt/current.xp>(alpha-1.0)*1.01
                hnew=previous.xp/previous.lc+dt*(alpha-1.0)*current.xp
                lc=current.xp/hnew
                isfinite(lc) && lc>0 || throw(ArgumentError(
                    "$caller: native Curve smoothing produced an invalid density"))
                points[i]=_IntegrationPoint(current.t,lc,current.p,current.xp,hnew)
                changed=true
            end
        end
        for i in length(points):-1:2
            current,previous=points[i],points[i-1]
            dh=previous.xp/previous.lc-current.xp/current.lc
            dt=abs(previous.t-current.t)
            if dh/dt/previous.xp>(alpha-1.0)*1.01
                hnew=current.xp/current.lc+dt*(alpha-1.0)*current.xp
                lc=previous.xp/hnew
                isfinite(lc) && lc>0 || throw(ArgumentError(
                    "$caller: native Curve smoothing produced an invalid density"))
                points[i-1]=_IntegrationPoint(previous.t,lc,previous.p,previous.xp,hnew)
                changed=true
            end
        end
        changed || break
    end
    @inbounds for i in 2:length(points)
        current,previous=points[i],points[i-1]
        p=previous.p+(current.t-previous.t)*0.5*(current.lc+previous.lc)
        isfinite(p) || throw(ArgumentError("$caller: native Curve smoothed primitive overflowed"))
        points[i]=_IntegrationPoint(current.t,current.lc,p,current.xp,current.h)
    end
    return points[end].p
end

function _gmsh_curve_parameters(points,count::Int,caller)
    2<=count<=10_000_001 || throw(ArgumentError(
        "$caller: native Curve node count exceeds its allocation capacity"))
    mass=points[end].p
    isfinite(mass) && mass!=0.0 || throw(ArgumentError(
        "$caller: native Curve integrates to a zero or nonfinite mass"))
    params=Vector{Float64}(undef,count)
    params[1],params[end]=points[1].t,points[end].t
    step=mass/(count-1)
    interval=2
    @inbounds for birth in 1:count-2
        target=birth*step
        while interval<=length(points)
            first,last=points[interval-1],points[interval]
            if abs(last.p)>=abs(target) && abs(first.p)<abs(target)
                t=first.t+(last.t-first.t)/(last.p-first.p)*(target-first.p)
                isfinite(t) || throw(ArgumentError(
                    "$caller: native Curve primitive inversion is not representable"))
                params[birth+1]=t
                break
            end
            interval+=1
        end
        interval<=length(points) || throw(ArgumentError(
            "$caller: native Curve primitive does not bracket a node mark"))
    end
    return params
end

@inline function _model_native_line_point(first,last,t::Float64)
    point=ntuple(axis->first[axis]+t*(last[axis]-first[axis]),3)
    # Preserve representable robust geometry at extreme coordinate scales.
    return all(isfinite,point) ? point :
        ntuple(axis->_convex_coordinate(first[axis],last[axis],t),3)
end

# Grading filters retain original producer rows. Native MVertex allocation
# precedes filtering, so a removed row still advances the public UID maximum.
struct _ModelNativeCurveGrade
    parameters::Vector{Float64}
    rows::Vector{Int}
    count::Int
end

function _model_native_line_grade_params(m::GeoModel,curve::Int,
                                         options::_ModelMesh1DOptions,caller;
                                         _native_producer::Bool=false,_birth_rows::Bool=false)
    a,b=m.curves[curve]
    first,last=m.points[a],m.points[b]
    gamma=t->_model_native_line_point(first,last,t)
    tangent=(last[1]-first[1],last[2]-first[2],last[3]-first[3])
    # meshGEdge's `F_Lc_aniso` dispatch (`algo2d == ALGO_2D_BAMG`, resolved
    # through `_model_curve_aniso_grading`): the directional `BGM_MeshMetric`
    # metric replaces `speed/h`, endpoint samples take the vertex metric, and
    # neither `smoothPrimitive` nor `filterPoints` applies.
    aniso=_model_curve_aniso_grading(m,curve)
    evaluate=function(t,at_begin,at_end)
        point=gamma(t)
        speed=_model_native_line_sampling_speed(first,last,t,gamma,caller)
        if aniso
            metric=if at_begin
                _model_vertex_bgm_metric(m,a,point,options,caller;_native_producer)
            elseif at_end
                _model_vertex_bgm_metric(m,b,point,options,caller;_native_producer)
            else
                _model_curve_bgm_metric(m,curve,t,point,0.0,1.0,options,caller;
                    _native_producer)
            end
            lc=_metric_displacement(metric,tangent[1],tangent[2],tangent[3],caller)
            isfinite(lc) && lc>0 || throw(ArgumentError(
                "$caller: native Line metric density is nonpositive or not representable"))
            return _IntegrationPoint(t,lc,0.0,speed,1.0)
        end
        endpoint=at_begin ? _model_vertex_bgm(m,a,point,options,caller;_native_producer) :
                 at_end ? _model_vertex_bgm(m,b,point,options,caller;_native_producer) : GMSH_MAX_SIZE
        h=min(endpoint,_model_curve_bgm(m,curve,t,point,0.0,1.0,options,caller;_native_producer))
        lc=speed/h
        isfinite(lc) && lc>0 || throw(ArgumentError(
            "$caller: native Line density is nonpositive or not representable"))
        return _IntegrationPoint(t,lc,0.0,speed,h)
    end
    points=_gmsh_curve_integration(evaluate,0.0,1.0,options.integration_precision,caller)
    try
    # Upstream applies `smoothPrimitive` only away from BAMG (`algo2d !=
    # ALGO_2D_BAMG`); the aniso primitive is used as integrated.
    mass=aniso ? points[end].p : _gmsh_curve_smooth_primitive!(points,caller)
    minimum=_model_minimum_curve_segments(m,curve,options,caller)
    nedge=_gmsh_edge_count(mass,minimum,10_000_000,caller)
    force_odd=m.meshing.recombine_algo!=0 && m.meshing.recombine_all
    odd_nodes=m.meshing.recombine_algo!=0 && (m.meshing.recombine_all ||
        any(face->haskey(m.meshing.recombine,(2,face)),_model_curve_faces(m,curve)))
    if odd_nodes && iseven(nedge+1) && mass>0.75
        nedge<10_000_000 || throw(ArgumentError("$caller: native Line edge count exceeds capacity"))
        nedge+=1
    end
    params=_gmsh_curve_parameters(points,nedge+1,caller)
    coordinates=NTuple{3,Float64}[gamma(t) for t in params]
    field=_ModelCurveEntityField(m,options,String(caller))
    size_function=t->_model_curve_bgm(m,curve,t,gamma(t),0.0,1.0,options,caller;_native_producer)
    count=length(params)-2
    rows=_birth_rows ? collect(eachindex(params)) : nothing
    if aniso
        # Equal-primitive marks invert a few ulps past the strict
        # `metric_edge_length` certificate — re-bisect interior nodes to the
        # conservative share (the endpoint stays at t1=1.0).
        metric_of=function(u,xyz,entity)
            dim,tag=entity===nothing ? (1,curve) :
                (Int(entity[1]),Int(entity[2]))
            dim==0 ? _model_vertex_bgm_metric(
                         m,tag,xyz,options,caller;_native_producer) :
                     _model_curve_bgm_metric(
                         m,tag,u,xyz,0.0,1.0,options,caller;_native_producer)
        end
        _undershoot_metric_edges!(gamma,metric_of,params,coordinates,1.0,
                                  (1,curve),mass/nedge,false,caller)
    else
        _filter_close_points!(coordinates,params,gamma,field,(1,curve),
                              minimum,false,force_odd,size_function,caller;row_map=rows)
    end
    return _birth_rows ? _ModelNativeCurveGrade(params,rows,count) : params
    finally
        _gmsh_curve_close_points!(points)
    end
end
