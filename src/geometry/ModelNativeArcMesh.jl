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
    # meshGEdge's `F_Lc_aniso` dispatch — the same arm
    # `_model_native_line_grade_params` carries: the directional
    # `BGM_MeshMetric` metric replaces `speed/h`, endpoint samples take the
    # vertex metric, and neither `smoothPrimitive` nor `filterPoints` applies.
    aniso=_model_curve_aniso_grading(m,curve)
    inner=options.inner_field
    # A Curve feeding an algorithm-7 face certifies every emitted edge by its
    # world-chord `metric_edge_length` (the BAMG unit bound 1.0) through the
    # surface producer's strict retained-boundary check — see
    # `_model_curve_grade_params`, whose tightening this mirrors. The emitted
    # discretization is measured against the producer's own certificate field
    # and re-placed under a progressively tightened grading metric while a
    # remnant persists, upstream of the Curve products it can no longer edit.
    cert_fields=Tuple{Int,Any}[]
    if aniso && inner isa AbstractAnisoField && !(inner isa _MetricScaledAniso)
        for face in _model_curve_faces(m,curve)
            _surface_mesh_algorithm(m,face)==7 || continue
            push!(cert_fields,(Int(face),
                _metric_surface_algorithm_field(m,face,options,nothing,caller)))
        end
    end
    minimum=_model_minimum_curve_segments(m,curve,options,caller)
    force_odd=m.meshing.recombine_algo!=0 && m.meshing.recombine_all
    odd_nodes=m.meshing.recombine_algo!=0 && (m.meshing.recombine_all ||
        any(face->haskey(m.meshing.recombine,(2,face)),_model_curve_faces(m,curve)))
    runopts=options
    scale=1.0
    local params,count,rows
    done=false
    for _attempt in 1:8
        evaluate=function(t,at_begin,at_end)
            point=gamma(t)
            derivative=_model_high_order_curve_derivative(adapter,t,caller)
            speed=sqrt(derivative[1]*derivative[1]+derivative[2]*derivative[2]+
                       derivative[3]*derivative[3])
            if aniso
                metric=if at_begin
                    _model_vertex_bgm_metric(m,a,point,runopts,caller;_native_producer=true)
                elseif at_end
                    _model_vertex_bgm_metric(m,b,point,runopts,caller;_native_producer=true)
                else
                    _model_curve_bgm_metric(m,curve,t,point,lower,upper,runopts,caller;
                        _native_producer=true,_native_geometry=adapter)
                end
                lc=_metric_displacement(metric,derivative[1],derivative[2],derivative[3],caller)
                isfinite(lc) && lc>0 || throw(ArgumentError(
                    "$caller: native Curve metric density is nonpositive or not representable"))
                return _IntegrationPoint(t,lc,0.0,speed,1.0)
            end
            endpoint=at_begin ? _model_vertex_bgm(m,a,point,runopts,caller;_native_producer=true) :
                     at_end ? _model_vertex_bgm(m,b,point,runopts,caller;_native_producer=true) : GMSH_MAX_SIZE
            h=min(endpoint,_model_curve_bgm(m,curve,t,point,lower,upper,runopts,caller;
                _native_producer=true,_native_geometry=adapter))
            lc=speed/h
            isfinite(lc) && lc>=0 || throw(ArgumentError(
                "$caller: native Curve density is negative or not representable"))
            return _IntegrationPoint(t,lc,0.0,speed,h)
        end
        points=_gmsh_curve_integration(evaluate,lower,upper,options.integration_precision,caller)
        try
        # Upstream applies `smoothPrimitive` only away from BAMG
        # (`algo2d != ALGO_2D_BAMG`); the aniso primitive is used as integrated.
        mass=aniso ? points[end].p : _gmsh_curve_smooth_primitive!(points,caller)
        nedge=_gmsh_edge_count(mass,minimum,10_000_000,caller)
        if odd_nodes && iseven(nedge+1) && mass>0.75
            nedge<10_000_000 || throw(ArgumentError("$caller: native Curve edge count exceeds capacity"))
            nedge+=1
        end
        params=_gmsh_curve_parameters(points,nedge+1,caller)
        coordinates=NTuple{3,Float64}[gamma(t) for t in params]
        field=_ModelCurveEntityField(m,runopts,String(caller))
        size_function=t->_model_curve_bgm(m,curve,t,gamma(t),lower,upper,runopts,caller;
            _native_producer=true,_native_geometry=adapter)
        count=length(params)-2
        rows=_birth_rows ? collect(eachindex(params)) : nothing
        if aniso
            # Equal-primitive marks invert a few ulps past the strict
            # `metric_edge_length` certificate — re-bisect interior nodes to
            # the conservative share. `params` is always open-style
            # `[t0,…,t1]` here — `_gmsh_curve_parameters` pins
            # `params[end]=points[end].t=t1` even when the curve is closed
            # (`a==b`) — so the closing edge is the last listed pair and the
            # list-convention flag is `false` (see the scalar arm below).
            metric_of=function(u,xyz,entity)
                dim,tag=entity===nothing ? (1,curve) :
                    (Int(entity[1]),Int(entity[2]))
                dim==0 ? _model_vertex_bgm_metric(
                             m,tag,xyz,runopts,caller;_native_producer=true) :
                         _model_curve_bgm_metric(
                             m,tag,u,xyz,lower,upper,runopts,caller;
                             _native_producer=true,_native_geometry=adapter)
            end
            _undershoot_metric_edges!(gamma,metric_of,params,coordinates,upper,
                                      (1,curve),mass/nedge,false,caller)
        else
            # `filterPoints` never removes a model vertex endpoint; passing
            # `closed=true` would expose the trailing `t1` station as a
            # removable close-point, silently truncating the stored
            # discretization against the `[t0,…,t1]` contract every consumer
            # holds (`_dim01_emit_curve!`, `_surface_curve_polygon_params`,
            # write-back span certificates).
            _filter_close_points!(coordinates,params,gamma,field,(1,curve),
                                  minimum,false,force_odd,size_function,caller;row_map=rows)
        end
        if isempty(cert_fields)
            done=true
        else
            worst=0.0
            for (face,cert) in cert_fields
                previous=gamma(first(params))
                for u in Iterators.drop(params,1)
                    point=gamma(u)
                    worst=max(worst,metric_edge_length(cert,previous,point;
                                                     entity=(2,face)))
                    previous=point
                end
            end
            if worst<=1.0
                done=true
            else
                scale=max(scale*1.02,worst*worst)
                runopts=_mesh1d_options_with_inner_field(options,
                    _MetricScaledAniso(inner,scale))
            end
        end
        finally
            _gmsh_curve_close_points!(points)
        end
        done && break
    end
    return _birth_rows ? _ModelNativeCurveGrade(params,rows,count) : params
end
