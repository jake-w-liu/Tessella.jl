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

# `gmsh.model.occ.addCircleArc(start, middle, end, center=true)` — the
# middle point is the circle's center; `center=false` makes it a point
# the arc passes through (`gce_MakeCirc` 3-point form upstream). The
# middle point never becomes an edge vertex either way.
function _model_occ_add_circle_arc!(m::GeoModel,start,middle,stop;
                                    tag=-1,center=true)
    caller="OCC.add_circle_arc"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,1,_model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(1,allocated)) && throw(ArgumentError(
        "$caller: OCC Curve[$allocated] already exists"))
    center isa Bool || throw(ArgumentError("$caller: center must be a boolean"))
    ids=UInt64[]; points=NTuple{3,Float64}[]
    for raw in (start,middle,stop)
        current=_model_occ_requested_tag(raw,caller)
        id=get(kernel.bindings,(0,current),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Point[$current]"))
        push!(ids,id)
        push!(points,kernel.shapes[id].geometry.point)
    end
    (geometry,error)=center ?
        _model_occ_circle_arc_geometry(points[1],points[2],points[3]) :
        _model_occ_circle_arc_through_geometry(points[1],points[2],points[3])
    error===nothing || throw(ArgumentError("$caller: $error"))
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,1,allocated,id,
        _ModelOCCShape(1,(kind=:circle,geometry=geometry),UInt64[ids[1],ids[3]]))
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# `gp_Dir` raises when the input modulus is at/below Precision::Confusion;
# `_occ_dir` would emit NaNs instead, so gate it like OCCT does.
function _model_occ_gp_dir(vector::NTuple{3,Float64},caller::AbstractString)
    _occ_modulus(vector)>_OCC_CONFUSION || throw(ArgumentError(
        "$caller: OpenCASCADE exception gp_Dir() : input vector has null modulus"))
    return _occ_dir(vector)
end

function _model_occ_axis_arg(value,what::AbstractString,caller::AbstractString)
    value===nothing && return nothing
    value isa AbstractVector || throw(ArgumentError(
        "$caller: $what must be a 3-element vector"))
    length(value) in (0,3) || throw(ArgumentError(
        "$caller: $what must be a 3-element vector"))
    isempty(value) && return nothing
    return ntuple(i->_finite_scalar(value[i],caller,what*"[$i]"),3)
end

# `gmsh.model.occ.addCircle(x,y,z,r,angle1=0,angle2=2π,zAxis,xAxis)` —
# center/radius/axis form. Full circles (`angle1==0 && angle2==2π` exactly,
# matching upstream's `MakeEdge(circ)` branch) collapse to one vertex;
# trimmed ranges emit one vertex per endpoint after OCCT's parameter swap.
function _model_occ_add_circle!(m::GeoModel,x,y,z,radius;tag=-1,
                                angle1=0.0,angle2=2π,
                                zaxis=nothing,xaxis=nothing)
    caller="OCC.add_circle"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,1,_model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(1,allocated)) && throw(ArgumentError(
        "$caller: OCC Curve[$allocated] already exists"))
    center=(_finite_scalar(x,caller,"x"),_finite_scalar(y,caller,"y"),
            _finite_scalar(z,caller,"z"))
    r=_finite_scalar(radius,caller,"radius")
    r>0.0 || throw(ArgumentError("$caller: Circle radius should be positive"))
    a1=_finite_scalar(angle1,caller,"angle1")
    a2=_finite_scalar(angle2,caller,"angle2")
    nv=_model_occ_axis_arg(zaxis,"zAxis",caller)
    xv=_model_occ_axis_arg(xaxis,"xAxis",caller)
    n,X,Y=if nv!==nothing
        normal=_model_occ_gp_dir(nv,caller)
        if xv!==nothing
            # gp_Ax2(center, N, Vx): X is Vx projected onto the plane ⊥ N;
            # a parallel pair projects to the null vector and OCCT raises.
            _occ_modulus(_occ_cross(normal,_model_occ_gp_dir(xv,caller)))>
                _OCC_CONFUSION || throw(ArgumentError(
                "$caller: OpenCASCADE exception gp_Ax2() : x direction is parallel to axis"))
            X,Y=_occ_ax2_setx(normal,_occ_dir(xv))
        else
            X,Y=_occ_ax2_frame(normal)
        end
        (normal,X,Y)
    else
        normal=(0.0,0.0,1.0)
        (normal,_occ_ax2_setx(normal,(1.0,0.0,0.0))...)
    end
    full=a1==0.0 && a2==2π
    full || a1!=a2 || throw(ArgumentError(
        "$caller: OpenCASCADE exception Geom_TrimmedCurve::U1 == U2"))
    t0,t1=full ? (0.0,2π) : minmax(a1,a2)
    geometry=(occ=:circle,center=center,n=n,X=X,Y=Y,r=r,t0=t0,t1=t1)
    shape_ids=_model_occ_reserve_shapes(kernel,full ? 2 : 3,caller)
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    kernel.maxima[3]<=typemax(Int32)-(full ? 1 : 2) || throw(ArgumentError(
        "$caller: OCC endpoint tag allocation exceeds Int32"))
    vertex_tags=[kernel.maxima[3]+i for i in 1:(full ? 1 : 2)]
    for (ordinal,parameter) in enumerate(full ? (t0,) : (t0,t1))
        push!(pending,(0,vertex_tags[ordinal],shape_ids[ordinal],
            _ModelOCCShape(0,(point=_occ_circle_point(geometry,parameter),
                              mesh_size=0.0),UInt64[])))
    end
    endpoints=full ? UInt64[shape_ids[1],shape_ids[1]] :
                     UInt64[shape_ids[1],shape_ids[2]]
    push!(pending,(1,allocated,shape_ids[end],
        _ModelOCCShape(1,(kind=:circle,geometry=geometry),endpoints)))
    for (dimension,current,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,current,id,shape)
    end
    kernel.next_shape=shape_ids[end]+UInt64(1)
    return allocated
end

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
