# Actual Curve vertices consume the sampled derivative and F_Transfinite
# primitive. Legacy geometry views keep their existing parameter arithmetic.
function _model_native_curve_endpoint_matches(m,tag,record)
    birth=get(m.meshing.native_curve_endpoints,tag,nothing)
    birth===nothing && return false
    isempty(record.node_tags) && return false
    UInt64(record.node_tags[1])==birth.node || return false
    for axis in 1:3
        reinterpret(UInt64,record.node_coords[axis,1])==
            reinterpret(UInt64,birth.point[axis]) || return false
    end
    return true
end

function _model_native_curve_endpoint_remap!(m,mapping)
    bindings=m.meshing.native_curve_endpoints
    for tag in collect(keys(bindings))
        birth=bindings[tag]
        if !haskey(mapping,birth.node)
            delete!(bindings,tag)
        else
            bindings[tag]=_ModelNativeCurveEndpoint(birth.curve,mapping[birth.node],birth.point)
        end
    end
    return nothing
end

function _model_native_curve_mesh_geometry(m::GeoModel,curve::Int,caller)
    kind=_curve_type(m,curve)
    if kind in (:tensor_bezier,:tensor_bspline)
        return _model_high_order_curve_geometry(m,curve,caller)
    end
    if kind===:line
        occ=get(m.curve_geometry,curve,nothing)
        if occ!==nothing
            hasproperty(occ,:occ) && occ.occ===:line || return nothing
        end
    else
        (kind in (:circle,:ellipse) || kind in _SPLINE_CURVE_TYPES) || return nothing
    end
    haskey(m.meshing.extrude_sources,(1,curve)) && return nothing
    return _model_high_order_curve_geometry(m,curve,caller)
end

function _model_native_curve_transfinite_params(m::GeoModel,curve::Int,
                                                spec,options,caller)
    adapter=_model_high_order_curve_geometry(m,curve,caller)
    return _model_native_curve_transfinite_with_geometry(m,curve,spec,options,adapter,caller)
end

function _model_native_curve_transfinite_with_geometry(m::GeoModel,curve::Int,
                                                       spec,options,
                                                       adapter::_ModelHighOrderCurveGeometry,caller)
    lower,upper=adapter.lower,adapter.upper
    count=_flexible_transfinite_nodes(m,spec.num_nodes,curve,caller)
    count==2 && return Float64[lower,upper]
    speed=parameter->begin
        derivative=_model_high_order_curve_derivative(adapter,parameter,caller)
        sqrt(derivative[1]*derivative[1]+derivative[2]*derivative[2]+derivative[3]*derivative[3])
    end
    length_evaluate=(parameter,_,_)->begin
        value=speed(parameter)
        _IntegrationPoint(parameter,value,0.0,value,1.0)
    end
    length_points=_gmsh_curve_integration(length_evaluate,lower,upper,
        options.integration_precision*options.ctx_lc,caller)
    geometric_length=try
        length_points[end].p
    finally
        _gmsh_curve_close_points!(length_points)
    end
    type,coefficient=_transfinite_effective_law(spec,geometric_length,spec.num_nodes)
    law_nodes=_transfinite_law_nodes(m,spec.num_nodes,caller)
    val=_transfinite_val(type,coefficient,geometric_length,law_nodes)
    density_evaluate=(parameter,_,_)->begin
        value=speed(parameter)
        density=val((parameter-lower)/(upper-lower),value)
        _IntegrationPoint(parameter,density,0.0,value,1.0)
    end
    points=_gmsh_curve_integration(density_evaluate,lower,upper,
        options.integration_precision,caller)
    try
        return _gmsh_curve_parameters(points,count,caller)
    finally
        _gmsh_curve_close_points!(points)
    end
end

function _model_native_arc_mesh_geometry(m::GeoModel,curve::Int,caller)
    _curve_type(m,curve) in (:circle,:ellipse) || return nothing
    return _model_native_curve_mesh_geometry(m,curve,caller)
end

_model_native_arc_transfinite_params(m::GeoModel,curve::Int,spec,options,caller)=
    _model_native_curve_transfinite_params(m,curve,spec,options,caller)

# Actual Curve vertices consume native GEdge point/firstDer and F_Lc.
# The standalone callable-curve and legacy geometry paths retain their kernel.
function _model_native_curve_grade_params(m::GeoModel,curve::Int,
                                          options::_ModelMesh1DOptions,caller;
                                          _birth_rows::Bool=false)
    adapter=_model_high_order_curve_geometry(m,curve,caller)
    return _model_native_curve_grade_with_geometry(m,curve,options,adapter,caller;_birth_rows)
end

function _model_native_curve_grade_with_geometry(m::GeoModel,curve::Int,
                                                 options::_ModelMesh1DOptions,
                                                 adapter::_ModelHighOrderCurveGeometry,caller;
                                                 _birth_rows::Bool=false)
    lower,upper=adapter.lower,adapter.upper
    a,b=m.curves[curve]
    gamma=t->_model_high_order_curve_point(adapter,t,caller)
    evaluate=function(t,at_begin,at_end)
        point=gamma(t)
        endpoint=at_begin ? _model_vertex_bgm(m,a,point,options,caller;_native_producer=true) :
                 at_end ? _model_vertex_bgm(m,b,point,options,caller;_native_producer=true) : GMSH_MAX_SIZE
        h=min(endpoint,_model_curve_bgm(m,curve,t,point,lower,upper,options,caller;
            _native_producer=true,_native_geometry=adapter))
        derivative=_model_high_order_curve_derivative(adapter,t,caller)
        speed=sqrt(derivative[1]*derivative[1]+derivative[2]*derivative[2]+
                   derivative[3]*derivative[3])
        lc=speed/h
        isfinite(lc) && lc>=0 || throw(ArgumentError(
            "$caller: native Curve density is negative or not representable"))
        return _IntegrationPoint(t,lc,0.0,speed,h)
    end
    points=_gmsh_curve_integration(evaluate,lower,upper,options.integration_precision,caller)
    try
    mass=_gmsh_curve_smooth_primitive!(points,caller)
    minimum=_model_minimum_curve_segments(m,curve,options,caller)
    nedge=_gmsh_edge_count(mass,minimum,10_000_000,caller)
    force_odd=m.meshing.recombine_algo!=0 && m.meshing.recombine_all
    odd_nodes=m.meshing.recombine_algo!=0 && (m.meshing.recombine_all ||
        any(face->haskey(m.meshing.recombine,(2,face)),_model_curve_faces(m,curve)))
    if odd_nodes && iseven(nedge+1) && mass>0.75
        nedge<10_000_000 || throw(ArgumentError("$caller: native Curve edge count exceeds capacity"))
        nedge+=1
    end
    params=_gmsh_curve_parameters(points,nedge+1,caller)
    coordinates=NTuple{3,Float64}[gamma(t) for t in params]
    field=_ModelCurveEntityField(m,options,String(caller))
    size_function=t->_model_curve_bgm(m,curve,t,gamma(t),lower,upper,options,caller;
        _native_producer=true,_native_geometry=adapter)
    count=length(params)-2
    rows=_birth_rows ? collect(eachindex(params)) : nothing
    # `params` is always open-style `[t0,…,t1]` here —
    # `_gmsh_curve_parameters` pins `params[end]=points[end].t=t1` even when
    # the curve is closed (`a==b`).  `filterPoints` never removes a model
    # vertex endpoint; passing `closed=true` would expose the trailing `t1`
    # station as a removable close-point, silently truncating the stored
    # discretization against the `[t0,…,t1]` contract every consumer holds
    # (`_dim01_emit_curve!`, `_surface_curve_polygon_params`, write-back
    # span certificates).  The honest list-convention flag is `false`.
    _filter_close_points!(coordinates,params,gamma,field,(1,curve),
                          minimum,false,force_odd,size_function,caller;row_map=rows)
    return _birth_rows ? _ModelNativeCurveGrade(params,rows,count) : params
    finally
        _gmsh_curve_close_points!(points)
    end
end
