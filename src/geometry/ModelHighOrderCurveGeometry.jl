# A high-order operation owns its geometry adapter. Mesh-vertex parameters and
# coordinates stay independent of the captured Curve parametrization.
struct _ModelHighOrderCurveGeometry{G}
    tag::Int
    kind::Symbol
    geometry::G
    lower::Float64
    upper::Float64
end

function _model_high_order_curve_geometry(m::GeoModel,tag::Int,caller)
    discrete=_model_discrete_curve_geometry(m,tag)
    discrete!==nothing && return _ModelHighOrderCurveGeometry(
        tag,:discrete,discrete,discrete.t0,discrete.t1)
    haskey(m.curves,tag) || return nothing
    occ=_occ_geometry_checked(m,tag,caller)
    if occ!==nothing
        hasproperty(occ,:native_geometry) && (occ=occ.native_geometry)
        lower,upper=_model_curve_param_bounds(m,tag,caller)
        return _ModelHighOrderCurveGeometry(tag,Symbol("occ_",occ.occ),occ,lower,upper)
    end
    kind=_curve_type(m,tag)
    if kind in (:line,:degenerate)
        first,last=m.curves[tag]
        frame=(first=m.points[first],last=m.points[last])
        return _ModelHighOrderCurveGeometry(tag,:line,frame,0.0,1.0)
    elseif kind in _SPLINE_CURVE_TYPES
        frame=_spline_geometry(m,tag,caller)
        return _ModelHighOrderCurveGeometry(tag,:spline,frame,frame.ubeg,frame.uend)
    elseif kind in (:circle,:ellipse)
        return _ModelHighOrderCurveGeometry(tag,:arc,_model_high_order_arc_geometry(m,tag,caller),0.0,1.0)
    end
    return nothing
end

# MinGW's paired sin/cos uses FSINCOS on Windows x86. Keep this native
# Julia-compiled evaluation local to actual high-order geometry; separate CRT
# calls differ by an ulp, which finite differences amplify in the Newton solve.
@inline function _model_high_order_curve_sincos(x::Float64)
    @static if Sys.iswindows() && Sys.ARCH in (:x86_64,:i686)
        if isfinite(x) && abs(x)<0x1p63
            return Base.llvmcall(raw"""
                %input = alloca double, align 8
                %output = alloca [2 x double], align 8
                store double %0, ptr %input, align 8
                %sin = getelementptr [2 x double], ptr %output, i32 0, i32 0
                %cos = getelementptr [2 x double], ptr %output, i32 0, i32 1
                call void asm sideeffect "fldl ($0); fsincos; fstpl ($1); fstpl ($2)", "r,r,r,~{st},~{st(1)},~{memory}"(ptr %input, ptr %cos, ptr %sin)
                %result = load [2 x double], ptr %output, align 8
                ret [2 x double] %result
            """,NTuple{2,Float64},Tuple{Float64},x)
        end
    end
    return _gm_sincos(x)
end

function _model_high_order_arc_point(frame,u)
    @static if Sys.iswindows()
        theta=frame.t1-(frame.t1-frame.t2)*u-frame.incl
        s,c=_model_high_order_curve_sincos(theta)
        si,ci=_model_high_order_curve_sincos(frame.incl)
        x=frame.f1*c*ci-frame.f2*s*si
        y=frame.f1*c*si+frame.f2*s*ci
        return ntuple(axis->(x*frame.u[axis]+y*frame.v[axis])+0.0*frame.n[axis]+frame.center[axis],3)
    else
        theta=fma(-(frame.t1-frame.t2),u,frame.t1)-frame.incl
        s,c=_model_high_order_curve_sincos(theta)
        si,ci=_model_high_order_curve_sincos(frame.incl)
        x=fma(frame.f1*c,ci,-(frame.f2*s*si))
        y=fma(frame.f1*c,si,frame.f2*s*ci)
        return ntuple(axis->fma(0.0,frame.n[axis],fma(x,frame.u[axis],y*frame.v[axis]))+frame.center[axis],3)
    end
end

function _model_high_order_curve_point(adapter::_ModelHighOrderCurveGeometry,u,caller)
    frame=adapter.geometry;kind=adapter.kind
    kind===:discrete && return _discrete_curve_point(frame.frames,u,caller)
    kind===:line && return ntuple(axis->frame.first[axis]+u*(frame.last[axis]-frame.first[axis]),3)
    if kind===:spline
        @static if Sys.iswindows()
            return _model_high_order_spline_point(frame,u)
        else
            return _spline_point(frame,u)
        end
    end
    kind===:arc && return _model_high_order_arc_point(frame,u)
    kind===:occ_circle && return _occ_circle_point(frame,u)
    kind===:occ_surface_iso && return _model_surface_iso_point(frame,u)
    kind===:occ_copied_curve && return first(_model_occ_curve_d2(frame.curve,u))
    kind===:occ_line && return _occ_line_point(frame,u)
    kind===:occ_degenerate && return _occ_degenerate_point(frame,adapter.tag,u)
    throw(ArgumentError("$caller: unsupported high-order Curve geometry $kind"))
end

function _model_high_order_curve_derivative(adapter::_ModelHighOrderCurveGeometry,u,caller)
    frame=adapter.geometry;kind=adapter.kind
    kind===:discrete && return _discrete_curve_derivative(frame.frames,u,caller)
    if kind in (:arc,:spline)
        curved_left_step=u<1e-8 ? 0.0 : 1e-8
        curved_right_step=u>1.0-1e-8 ? 0.0 : 1e-8
        curved_left=_model_high_order_curve_point(adapter,u-curved_left_step,caller)
        curved_right=_model_high_order_curve_point(adapter,u+curved_right_step,caller)
        curved_reciprocal=1.0/(curved_left_step+curved_right_step)
        return ntuple(axis->(curved_right[axis]-curved_left[axis])*curved_reciprocal,3)
    end
    kind===:occ_circle && return _occ_circle_derivative(frame,u)
    kind===:occ_surface_iso && return _model_surface_iso_d1(frame,u)
    kind===:occ_copied_curve && return _model_occ_curve_d2(frame.curve,u)[2]
    kind===:occ_line && return frame.dir
    kind===:occ_degenerate && return _occ_degenerate_d1(frame,adapter.tag,u)
    if kind===:line
        left_step=u<1e-5 ? 0.0 : 1e-5
        right_step=u>1.0-1e-5 ? 0.0 : 1e-5
        left=_model_high_order_curve_point(adapter,u-left_step,caller)
        right=_model_high_order_curve_point(adapter,u+right_step,caller)
        reciprocal=1.0/(left_step+right_step)
        return ntuple(axis->(right[axis]-left[axis])*reciprocal,3)
    end
    throw(ArgumentError("$caller: unsupported high-order Curve derivative $kind"))
end

# Built-in gmshEdge uses a finite difference of firstDer; OCCEdge inherits
# GEdge's larger, bounds-aware step instead.
function _model_native_curve_second_derivative(adapter::_ModelHighOrderCurveGeometry,u,caller)
    kind=adapter.kind
    kind===:discrete && return _model_discrete_curve_second_derivative(adapter.geometry,u,caller)
    if kind in (:line,:arc,:spline)
        step=kind===:line ? 1e-5 : 1e-8
        left_step=u<step ? 0.0 : step
        right_step=u>1.0-step ? 0.0 : step
        first=_model_high_order_curve_derivative(adapter,u-left_step,caller)
        last=_model_high_order_curve_derivative(adapter,u+right_step,caller)
        reciprocal=1.0/(left_step+right_step)
        return ((last[1]-first[1])*reciprocal,
                (last[2]-first[2])*reciprocal,
                (last[3]-first[3])*reciprocal)
    end
    if u-1e-3<=adapter.lower
        first=_model_high_order_curve_derivative(adapter,u,caller)
        last=_model_high_order_curve_derivative(adapter,u+1e-3,caller)
        scale=1000.0
    elseif u+1e-3>=adapter.upper
        first=_model_high_order_curve_derivative(adapter,u-1e-3,caller)
        last=_model_high_order_curve_derivative(adapter,u,caller)
        scale=1000.0
    else
        first=_model_high_order_curve_derivative(adapter,u-1e-3,caller)
        last=_model_high_order_curve_derivative(adapter,u+1e-3,caller)
        scale=500.0
    end
    return (scale*(last[1]-first[1]),
            scale*(last[2]-first[2]),
            scale*(last[3]-first[3]))
end

function _model_native_curve_curvature(m::GeoModel,adapter::_ModelHighOrderCurveGeometry,u,caller)
    kind=adapter.kind
    kind===:line && return 0.0
    if kind===:discrete
        boundary=m.discrete[(1,adapter.tag)].boundary
        closed=!isempty(boundary) && first(boundary)[1]==last(boundary)[1] &&
            abs(first(boundary)[2])==abs(last(boundary)[2])
        return _model_discrete_curve_curvature(adapter.geometry,u,closed,caller)
    end
    if kind in (:occ_circle,:occ_line,:occ_degenerate)
        return _occ_curve_curvature(m,adapter.geometry,adapter.tag,u,caller)
    end
    d1=_model_high_order_curve_derivative(adapter,u,caller)
    d2=_model_native_curve_second_derivative(adapter,u,caller)
    cross=(d1[2]*d2[3]-d2[2]*d1[3],
           -(d1[1]*d2[3]-d2[1]*d1[3]),
           d1[1]*d2[2]-d2[1]*d1[2])
    cross_norm=sqrt(cross[1]*cross[1]+cross[2]*cross[2]+cross[3]*cross[3])
    reciprocal=1.0/sqrt(d1[1]*d1[1]+d1[2]*d1[2]+d1[3]*d1[3])
    return cross_norm*((reciprocal*reciprocal)*reciprocal)
end

# Public built-in Curve queries share the actual producer geometry. Other
# entities retain their existing surface/discrete/OCC evaluation dispatch.
function _model_native_curve_evaluate(m::GeoModel,dim,tag,parameters,operation)
    dim==1 || return nothing
    caller="model_$(operation)"
    _,entity=_model_evaluation_entity(m,dim,tag,(1,),caller)
    adapter=_model_high_order_curve_geometry(m,entity,caller)
    adapter===nothing && return nothing
    (adapter.kind in (:line,:arc,:spline) ||
        (adapter.kind===:occ_circle && hasproperty(adapter.geometry,:native_arithmetic))) || return nothing
    values=_model_evaluation_values(parameters,1,caller,"parametric coordinates")
    output=Float64[]
    sizehint!(output,(operation===:curvature ? 1 : 3)*length(values))
    for parameter in values
        if operation===:curvature
            push!(output,_model_native_curve_curvature(m,adapter,parameter,caller))
        elseif operation===:value
            _model_append_point!(output,_model_high_order_curve_point(adapter,parameter,caller))
        elseif operation===:derivative
            append!(output,_model_high_order_curve_derivative(adapter,parameter,caller))
        else
            append!(output,_model_native_curve_second_derivative(adapter,parameter,caller))
        end
    end
    return output
end

# GEdge::XYZToU uses the first converging seed and its signed du guard.
# Actual periodic copies consume this protocol rather than exact chord
# projection; the slave owns the resulting parameter and geometry sample.
function _model_native_curve_transform(point,affine)
    return ntuple(3) do axis
        result=0.0
        for component in 1:4
            value=component==4 ? 1.0 : point[component]
            result+=affine[4*(axis-1)+component]*value
        end
        result
    end
end

function _model_native_curve_refine_projection(adapter,target,parameter,relaxation,lc,caller)
    point=_model_high_order_curve_point(adapter,parameter,caller)
    difference=(point[1]-target[1],point[2]-target[2],point[3]-target[3])
    error=sqrt(difference[1]*difference[1]+difference[2]*difference[2]+difference[3]*difference[3])
    threshold=1e-8*lc
    for iteration in 1:25
        error>threshold || break
        derivative=_model_high_order_curve_derivative(adapter,parameter,caller)
        correction=(difference[1]*derivative[1]+difference[2]*derivative[2]+difference[3]*derivative[3])/
            (derivative[1]*derivative[1]+derivative[2]*derivative[2]+derivative[3]*derivative[3])
        correction<1e-8 && error>threshold && (correction=1.0)
        parameter=min(adapter.upper,max(adapter.lower,parameter-relaxation*correction))
        point=_model_high_order_curve_point(adapter,parameter,caller)
        difference=(point[1]-target[1],point[2]-target[2],point[3]-target[3])
        error=sqrt(difference[1]*difference[1]+difference[2]*difference[2]+difference[3]*difference[3])
    end
    return error<=threshold,parameter,error
end

function _model_native_curve_xyz_to_u(adapter,target,lc,caller;relaxation=1.0)
    errors=Dict{Float64,Float64}()
    for index in 0:20
        parameter=adapter.lower+(adapter.upper-adapter.lower)/20*index
        converged,result,error=_model_native_curve_refine_projection(
            adapter,target,parameter,relaxation,lc,caller)
        converged && return true,result
        errors[error]=result
    end
    if relaxation>0.1
        converged,result=_model_native_curve_xyz_to_u(
            adapter,target,lc,caller;relaxation=0.75*relaxation)
        converged && return true,result
        point=_model_high_order_curve_point(adapter,result,caller)
        difference=ntuple(axis->point[axis]-target[axis],3)
        error=sqrt(difference[1]*difference[1]+difference[2]*difference[2]+difference[3]*difference[3])
        errors[error]=result
    end
    return false,errors[minimum(keys(errors))]
end

# Built-in Line queries use GEdge::parFromPoint, including its native seed,
# clamping and convergence protocol. Exact Model line inversion stays separate.
function _model_native_curve_parametrization(m::GeoModel,dim,tag,coordinates)
    dim==1 || return nothing
    caller="model_parametrization"
    _,entity=_model_evaluation_entity(m,dim,tag,(1,),caller)
    haskey(m.discrete,(1,entity)) && return nothing
    _curve_type(m,entity)===:line || return nothing
    adapter=_model_high_order_curve_geometry(m,entity,caller)
    adapter===nothing && return nothing
    adapter.kind===:line || return nothing
    values=_model_evaluation_values(coordinates,3,caller,"coordinates")
    output=Vector{Float64}(undef,length(values)÷3)
    lc=_model_lc(m)
    for index in eachindex(output)
        offset=3*(index-1)
        point=(values[offset+1],values[offset+2],values[offset+3])
        output[index]=_model_native_curve_xyz_to_u(adapter,point,lc,caller)[2]
    end
    return output
end

function _model_native_curve_copy_parameter(m,tag,adapter,point,caller)
    if adapter.kind in (:line,:arc,:spline)
        return _model_native_curve_xyz_to_u(adapter,point,_model_lc(m),caller)[2]
    end
    return _model_high_order_curve_inverse(m,tag,adapter,point,caller)
end

function _model_high_order_curve_inverse(m,tag,adapter,point,caller)
    frame=adapter.geometry;kind=adapter.kind
    kind===:discrete && return _model_discrete_curve_xyz_to_u(frame,point,_model_lc(m),caller)[2]
    kind===:arc && return _arc_xyz_to_u(frame,point,_model_lc(m))[2]
    kind===:spline && return _spline_xyz_to_u(frame,point,_model_lc(m))[2]
    if kind===:line
        line=_model_line_geometry(m,tag,caller)
        parameter,_=_model_line_parameter_exact(line,point)
        return _model_rational_float(parameter,caller,"Curve[$tag] support parameter")
    end
    projection=_occ_curve_project(m,frame,tag,point,caller)
    return projection===nothing ? _occ_curve_xyz_to_u(m,frame,tag,point,_model_lc(m),caller)[2] : projection[1]
end

# Literal GEdge::length(...,10), followed by the scalar N=3 Newton protocol in
# HighOrder.cpp. Do not split a piecewise Curve at geometric knots.
function _model_high_order_curve_length10(adapter,u0,u1,caller)
    total=0.0;jacobian=(u1-u0)*0.5
    for index in 1:10
        abscissa=_DISCRETE_CURVE_GL10_POINTS[index]
        u=u0*0.5*(1.0-abscissa)+u1*0.5*(1.0+abscissa)
        derivative=_model_high_order_curve_derivative(adapter,u,caller)
        speed=sqrt(derivative[1]*derivative[1]+derivative[2]*derivative[2]+derivative[3]*derivative[3])
        total+=speed*_DISCRETE_CURVE_GL10_WEIGHTS[index]*jacobian
    end
    return total
end

function _model_high_order_curve_equidistant(adapter,u0,u1,relaxation,caller)
    increment=(u1-u0)/2.0;parameter=u0+increment;upper=parameter+increment
    epsilon=(u1-u0)*1e-5
    for iteration in 2:50
        residual=_model_high_order_curve_length10(adapter,parameter,upper,caller)-
                 _model_high_order_curve_length10(adapter,u0,parameter,caller)
        parameter+=epsilon
        perturbed=_model_high_order_curve_length10(adapter,parameter,upper,caller)-
                  _model_high_order_curve_length10(adapter,u0,parameter,caller)
        jacobian=(perturbed-residual)/epsilon
        parameter-=epsilon
        correction=residual/jacobian;parameter-=relaxation*correction
        (parameter<u0 || parameter>u1) && break
        abs(correction)<1e-6 && return true,parameter
    end
    return false,parameter
end

function _model_high_order_curve_p2_parameter(adapter,first,last,caller)
    lower,upper=minmax(first,last);relaxation=1.0
    while relaxation>1e-2
        converged,parameter=_model_high_order_curve_equidistant(adapter,lower,upper,relaxation,caller)
        converged && return parameter
        relaxation*=0.5
    end
    return lower+(upper-lower)/2.0
end

function _model_high_order_curve_residual!(residual,lengths,adapter,parameters,caller)
    for index in eachindex(lengths)
        lengths[index]=_model_high_order_curve_length10(
            adapter,parameters[index],parameters[index+1],caller)
    end
    for index in eachindex(residual)
        residual[index]=lengths[index+1]-lengths[index]
    end
    return residual
end

# HighOrder.cpp uses a row for each perturbed parameter and a column for each
# residual. Preserve that orientation, the accumulated endpoint, and the
# relaxation schedule when creating more than one interior Curve vertex.
function _model_high_order_curve_equidistant_parameters(adapter,lower,upper,
                                                        order,relaxation,caller)
    parameters=Vector{Float64}(undef,order+1)
    parameters[1]=lower
    increment=(upper-lower)/order
    for index in 2:order+1
        parameters[index]=parameters[index-1]+increment
    end
    count=order-1
    lengths=Vector{Float64}(undef,order)
    residual=Vector{Float64}(undef,count)
    perturbed=similar(residual)
    jacobian=Matrix{Float64}(undef,count,count)
    epsilon=(upper-lower)*1e-5
    for iteration in 2:50
        _model_high_order_curve_residual!(residual,lengths,adapter,parameters,caller)
        for index in 1:count
            parameters[index+1]+=epsilon
            _model_high_order_curve_residual!(perturbed,lengths,adapter,parameters,caller)
            for column in 1:count
                jacobian[index,column]=(perturbed[column]-residual[column])/epsilon
            end
            parameters[index+1]-=epsilon
        end
        # The pinned Windows Gmsh build selects Eigen's column-pivoted QR
        # implementation of fullMatrix::luSolve. Julia owns this factorization.
        correction=_model_high_order_curve_qr_solve(jacobian,residual)
        norm_squared=0.0
        for index in 1:count
            parameters[index+1]-=relaxation*correction[index]
            norm_squared+=correction[index]*correction[index]
        end
        (parameters[2]<lower || parameters[order]>upper) && break
        sqrt(norm_squared)<1e-6 && return true,parameters
    end
    return false,parameters
end

function _model_high_order_curve_parameters(adapter,first,last,order::Int,caller)
    2<=order<=10 || throw(ArgumentError("$caller: Curve order must be between 2 and 10"))
    order==2 && return Float64[_model_high_order_curve_p2_parameter(adapter,first,last,caller)]
    lower,upper=minmax(first,last)
    relaxation=1.0
    parameters=Float64[]
    converged=false
    while relaxation>1e-2
        converged,parameters=_model_high_order_curve_equidistant_parameters(
            adapter,lower,upper,order,relaxation,caller)
        converged && break
        relaxation*=0.5
    end
    if !converged
        parameters[1]=lower
        increment=(upper-lower)/order
        for index in 2:order
            parameters[index]=parameters[index-1]+increment
        end
    end
    interior=parameters[2:order]
    first<last || reverse!(interior)
    return interior
end

# Small column-pivoted Householder QR for the native high-order Newton system.
# Keep two scalar lanes and the same update order as the pinned SSE2 build.
function _model_high_order_curve_qr_dot(first,second)
    count=length(first)
    count==0 && return 0.0
    count==1 && return first[1]*second[1]
    lane0=first[1]*second[1]
    lane1=first[2]*second[2]
    aligned=2*(count÷2)
    aligned4=4*(count÷4)
    if aligned>2
        other0=first[3]*second[3]
        other1=first[4]*second[4]
        for index in 5:4:aligned4
            lane0+=first[index]*second[index]
            lane1+=first[index+1]*second[index+1]
            other0+=first[index+2]*second[index+2]
            other1+=first[index+3]*second[index+3]
        end
        lane0+=other0
        lane1+=other1
        if aligned>aligned4
            lane0+=first[aligned4+1]*second[aligned4+1]
            lane1+=first[aligned4+2]*second[aligned4+2]
        end
    end
    result=lane0+lane1
    for index in aligned+1:count
        result+=first[index]*second[index]
    end
    return result
end

# Eigen's matrix reflector uses a matrix-vector product: each packet lane
# accumulates consecutive pairs. Norms and the right-hand-side reflector use
# Redux's two packet accumulators instead, so retain their separate dot above.
function _model_high_order_curve_qr_matrix_dot(first,second)
    count=length(first)
    count==0 && return 0.0
    count==1 && return first[1]*second[1]
    lane0=0.0
    lane1=0.0
    aligned=2*(count÷2)
    for index in 1:2:aligned
        lane0+=first[index]*second[index]
        lane1+=first[index+1]*second[index+1]
    end
    result=lane0+lane1
    for index in aligned+1:count
        result+=first[index]*second[index]
    end
    return result
end

function _model_high_order_curve_qr_reflect!(matrix,pivot,column,tau)
    count=size(matrix,1)
    pivot==count && return nothing
    tau==0.0 && return nothing
    essential=@view matrix[pivot+1:count,pivot]
    tail=@view matrix[pivot+1:count,column]
    temporary=_model_high_order_curve_qr_matrix_dot(essential,tail)+matrix[pivot,column]
    matrix[pivot,column]-=tau*temporary
    for row in pivot+1:count
        matrix[row,column]-=(tau*matrix[row,pivot])*temporary
    end
    return nothing
end

function _model_high_order_curve_qr_solve(matrix,residual)
    count=length(residual)
    qr=copy(matrix)
    direct=Vector{Float64}(undef,count)
    updated=similar(direct)
    coefficients=similar(direct)
    permutation=collect(1:count)
    for column in 1:count
        values=@view qr[:,column]
        direct[column]=sqrt(_model_high_order_curve_qr_dot(values,values))
        updated[column]=direct[column]
    end
    threshold=(maximum(updated)*eps(Float64))^2/count
    downdate_threshold=sqrt(eps(Float64))
    nonzero=count
    for pivot in 1:count
        largest=pivot
        for column in pivot+1:count
            updated[column]>updated[largest] && (largest=column)
        end
        nonzero==count && updated[largest]^2<threshold*(count-pivot+1) &&
            (nonzero=pivot-1)
        if largest!=pivot
            for row in 1:count
                qr[row,pivot],qr[row,largest]=qr[row,largest],qr[row,pivot]
            end
            updated[pivot],updated[largest]=updated[largest],updated[pivot]
            direct[pivot],direct[largest]=direct[largest],direct[pivot]
            permutation[pivot],permutation[largest]=permutation[largest],permutation[pivot]
        end
        tail=@view qr[pivot+1:count,pivot]
        squared=_model_high_order_curve_qr_dot(tail,tail)
        first=qr[pivot,pivot]
        tau=0.0
        beta=first
        if squared<=floatmin(Float64)
            fill!(tail,0.0)
        else
            beta=sqrt(first*first+squared)
            first>=0.0 && (beta=-beta)
            divisor=first-beta
            for row in pivot+1:count
                qr[row,pivot]/=divisor
            end
            tau=(beta-first)/beta
        end
        coefficients[pivot]=tau
        qr[pivot,pivot]=beta
        for column in pivot+1:count
            _model_high_order_curve_qr_reflect!(qr,pivot,column,tau)
            if updated[column]!=0.0
                temporary=abs(qr[pivot,column])/updated[column]
                temporary=(1.0+temporary)*(1.0-temporary)
                temporary<0.0 && (temporary=0.0)
                accuracy=temporary*(updated[column]/direct[column])^2
                if accuracy<=downdate_threshold
                    values=@view qr[pivot+1:count,column]
                    direct[column]=sqrt(_model_high_order_curve_qr_dot(values,values))
                    updated[column]=direct[column]
                else
                    updated[column]*=sqrt(temporary)
                end
            end
        end
    end
    right=copy(residual)
    for pivot in 1:nonzero
        tau=coefficients[pivot]
        if pivot==count
            right[pivot]*=1.0-tau
        elseif tau!=0.0
            essential=@view qr[pivot+1:count,pivot]
            tail=@view right[pivot+1:count]
            temporary=_model_high_order_curve_qr_dot(essential,tail)+right[pivot]
            right[pivot]-=tau*temporary
            for row in pivot+1:count
                right[row]-=(tau*qr[row,pivot])*temporary
            end
        end
    end
    # Eigen substitutes within eight-column panels, then applies the solved
    # panel to earlier rows with one matrix-vector accumulation. This matters
    # for the nine-unknown order-ten Newton system.
    panel_end=nonzero
    while panel_end>0
        panel_start=max(1,panel_end-7)
        for column in panel_end:-1:panel_start
            if right[column]!=0.0
                right[column]/=qr[column,column]
                for row in panel_start:column-1
                    right[row]-=right[column]*qr[row,column]
                end
            end
        end
        for row in 1:panel_start-1
            product=0.0
            for column in panel_start:panel_end
                product+=qr[row,column]*right[column]
            end
            right[row]+=(-1.0)*product
        end
        panel_end=panel_start-1
    end
    solution=zeros(Float64,count)
    for index in 1:nonzero
        solution[permutation[index]]=right[index]
    end
    return solution
end
