# Pending OCC edges use copied geometry and explicit endpoint shape references.
# Public Point coordinates and control tags are not their evaluator storage.
function _model_occ_add_curve!(m::GeoModel,points,family;tag=-1,degree=3,
    weights=Float64[],knots=Float64[],multiplicities=Int[],tolerance=1e-8)
    caller=family===:bezier ? "OCC.add_bezier" : "OCC.add_bspline"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,1,_model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(1,allocated)) && throw(ArgumentError(
        "$caller: OCC Curve[$allocated] already exists"))
    tags=Int[_model_occ_requested_tag(t,caller) for t in points]
    length(tags)>=2 || throw(ArgumentError("$caller: at least two controls are required"))
    poles=NTuple{3,Float64}[]
    point_ids=UInt64[]
    for point in tags
        id=get(kernel.bindings,(0,point),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Point[$point]"))
        push!(point_ids,id)
        push!(poles,kernel.shapes[id].geometry.point)
    end
    family===:bspline && first(tags)==last(tags) && throw(ArgumentError(
        "$caller: periodic copied OCC Curve construction is not implemented"))
    curve=_model_occ_curve_geometry(family,poles;degree,weights,knots,multiplicities,caller)
    endpoint_tolerance=_finite_scalar(tolerance,caller,"geometric tolerance")
    endpoint_tolerance>=0.0 || throw(ArgumentError("$caller: geometric tolerance must be nonnegative"))
    lower,upper=_model_occ_curve_bounds(curve)
    geometry=(occ=:copied_curve,curve=curve,t0=lower,t1=upper)
    reuse=family===:bezier ||
        (sqrt(sum((first(_model_occ_curve_d2(curve,lower))[axis]-first(poles)[axis])^2 for axis in 1:3))<=endpoint_tolerance &&
         sqrt(sum((first(_model_occ_curve_d2(curve,upper))[axis]-last(poles)[axis])^2 for axis in 1:3))<=endpoint_tolerance &&
         sqrt(sum((first(poles)[axis]-last(poles)[axis])^2 for axis in 1:3))>endpoint_tolerance)
    shape_ids=_model_occ_reserve_shapes(kernel,reuse ? 1 : 3,caller)
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    endpoints=if reuse
        UInt64[first(point_ids),last(point_ids)]
    else
        maximum=kernel.maxima[3]
        maximum<=typemax(Int32)-2 || throw(ArgumentError(
            "$caller: OCC endpoint tag allocation exceeds Int32"))
        size=kernel.shapes[first(point_ids)].geometry.mesh_size
        for (ordinal,parameter) in enumerate((lower,upper))
            point=first(_model_occ_curve_d2(curve,parameter))
            push!(pending,(0,maximum+ordinal,shape_ids[ordinal],
                _ModelOCCShape(0,(point=point,mesh_size=size),UInt64[])))
        end
        UInt64[shape_ids[1],shape_ids[2]]
    end
    kind=family===:bezier ? :tensor_bezier : :tensor_bspline
    push!(pending,(1,allocated,shape_ids[end],
        _ModelOCCShape(1,(kind=kind,geometry=geometry),endpoints)))
    for (dimension,current,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,current,id,shape)
    end
    kernel.next_shape=shape_ids[end]+UInt64(1)
    return allocated
end

_model_occ_add_bezier!(m,points;kwargs...)=_model_occ_add_curve!(m,points,:bezier;kwargs...)
_model_occ_add_bspline!(m,points;kwargs...)=_model_occ_add_curve!(m,points,:bspline;kwargs...)

function _model_occ_add_line!(m::GeoModel,start,stop;tag=-1)
    caller="OCC.add_line"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,1,_model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(1,allocated)) && throw(ArgumentError(
        "$caller: OCC Curve[$allocated] already exists"))
    endpoints=UInt64[]
    for raw in (start,stop)
        current=_model_occ_requested_tag(raw,caller)
        id=get(kernel.bindings,(0,current),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Point[$current]"))
        push!(endpoints,id)
    end
    first_point,last_point=(kernel.shapes[id].geometry.point for id in endpoints)
    difference=ntuple(axis->last_point[axis]-first_point[axis],3)
    length=_occ_modulus(difference)
    length>_OCC_CONFUSION || throw(ArgumentError("$caller: cannot construct a coincident edge"))
    geometry=(occ=:line,origin=first_point,dir=_occ_dir(difference),t0=0.0,t1=length)
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,1,allocated,id,
        _ModelOCCShape(1,(kind=:line,geometry=geometry),endpoints))
    kernel.next_shape=id+UInt64(1)
    return allocated
end
