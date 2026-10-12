# OCC.add_spline mirrors GeomAPI_Interpolate's interpolation path: chord-length
# parameters, Lagrange endpoint tangents, and a B-spline collocation solve.
# Periodic input (first tag == last tag) uses the wrapped knot vector OCCT
# solves internally, stored flat so the whole curve stays evaluable.

# Nonzero basis values and first derivatives on a flat knot vector. Returns the
# 1-based span index s with entries aligned to pole indices s-degree..s.
function _model_occ_bspline_basis_all(knots::Vector{Float64},degree::Int,u::Float64)
    count=length(knots)
    poles=count-degree-1
    high=knots[poles+1]
    s=if u>=high
        poles
    elseif u<=knots[degree+1]
        degree+1
    else
        clamp(searchsortedlast(knots,u),degree+1,poles)
    end
    function values_at(order)
        result=Vector{Float64}(undef,order+1)
        left=Vector{Float64}(undef,order+1)
        right=Vector{Float64}(undef,order+1)
        result[1]=1.0
        for j in 1:order
            left[j]=u-knots[s+1-j]
            right[j]=knots[s+j]-u
            saved=0.0
            for r in 0:j-1
                denominator=right[r+1]+left[j-r]
                temp=denominator==0.0 ? 0.0 : result[r+1]/denominator
                result[r+1]=saved+right[r+1]*temp
                saved=left[j-r]*temp
            end
            result[j+1]=saved
        end
        return result
    end
    values=values_at(degree)
    if degree==0
        return s,values,zeros(Float64,degree+1)
    end
    lower=values_at(degree-1)
    derivatives=Vector{Float64}(undef,degree+1)
    for k in 0:degree
        j=s-degree+k
        first_denom=knots[j+degree]-knots[j]
        second_denom=knots[j+degree+1]-knots[j+1]
        first_term=first_denom==0.0 ? 0.0 :
            degree*(k>=1 ? lower[k] : 0.0)/first_denom
        second_term=second_denom==0.0 ? 0.0 :
            degree*(k<degree ? lower[k+1] : 0.0)/second_denom
        derivatives[k+1]=first_term-second_term
    end
    return s,values,derivatives
end

# Derivative at t of the degree-`degree` Lagrange interpolant through the first
# degree+1 (point, parameter) pairs. Mirrors PLib::EvalLagrange order 1.
function _model_occ_lagrange_tangent(points,parameters,degree,t)
    result=(0.0,0.0,0.0)
    for j in 1:degree+1
        denominator=1.0
        for k in 1:degree+1
            k==j || (denominator*=parameters[j]-parameters[k])
        end
        derivative=0.0
        for m in 1:degree+1
            m==j && continue
            term=1.0
            for k in 1:degree+1
                (k==j || k==m) && continue
                term*=t-parameters[k]
            end
            derivative+=term/denominator
        end
        point=points[j]
        result=(result[1]+derivative*point[1],result[2]+derivative*point[2],
                result[3]+derivative*point[3])
    end
    return result
end

# Solve the B-spline collocation system: one row per constraint, where
# orders[i]==0 interpolates a point and orders[i]==1 interpolates a derivative.
# Returns the poles as NTuple{3,Float64}; throws on a singular system.
function _model_occ_spline_solve(flatknots,degree,taus,orders,values,caller)
    poles=length(flatknots)-degree-1
    length(taus)==poles && length(values)==poles ||
        throw(ArgumentError("$caller: OCC spline constraint count mismatch"))
    matrix=zeros(Float64,poles,poles)
    rhs=zeros(Float64,poles,3)
    for i in 1:poles
        s,basis,derivatives=_model_occ_bspline_basis_all(flatknots,degree,taus[i])
        row=orders[i]==0 ? basis : derivatives
        for k in 1:degree+1
            matrix[i,s-degree+k-1]=row[k]
        end
        value=values[i]
        rhs[i,1]=value[1];rhs[i,2]=value[2];rhs[i,3]=value[3]
    end
    controls=try
        matrix\rhs
    catch err
        err isa InterruptException && rethrow()
        err isa LinearAlgebra.SingularException ||
            err isa LinearAlgebra.LAPACK.LAPACKException || rethrow()
        throw(ArgumentError("$caller: could not interpolate spline"))
    end
    all(isfinite,controls) || throw(ArgumentError(
        "$caller: could not interpolate spline"))
    return NTuple{3,Float64}[(controls[i,1],controls[i,2],controls[i,3])
                             for i in 1:poles]
end

function _model_occ_add_spline!(m::GeoModel,points;tag=-1,tolerance=1e-8)
    caller="OCC.add_spline"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,1,_model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(1,allocated)) && throw(ArgumentError(
        "$caller: OCC Curve[$allocated] already exists"))
    tags=Int[_model_occ_requested_tag(t,caller) for t in points]
    length(tags)>=2 || throw(ArgumentError(
        "$caller: OCC spline requires at least two points"))
    ids=UInt64[]
    coordinates=NTuple{3,Float64}[]
    for current in tags
        id=get(kernel.bindings,(0,current),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Point[$current]"))
        push!(ids,id)
        push!(coordinates,kernel.shapes[id].geometry.point)
    end
    geometric_tolerance=_finite_scalar(tolerance,caller,"geometric tolerance")
    geometric_tolerance>=0.0 || throw(ArgumentError(
        "$caller: geometric tolerance must be nonnegative"))
    periodic=first(tags)==last(tags)
    count=periodic ? length(tags)-1 : length(tags)
    count<2 && throw(ArgumentError(
        "$caller: could not interpolate spline"))
    samples=coordinates[1:count]
    for i in 1:count-1
        gap=samples[i+1][1]-samples[i][1]
        squared=gap*gap
        for axis in 2:3
            gap=samples[i+1][axis]-samples[i][axis]
            squared+=gap*gap
        end
        squared>=geometric_tolerance^2 || throw(ArgumentError(
            "$caller: OCC spline points must be at least $geometric_tolerance apart"))
    end
    parameters=Vector{Float64}(undef,periodic ? count+1 : count)
    parameters[1]=0.0
    for i in 1:count-1
        delta=samples[i+1][1]-samples[i][1]
        squared=delta*delta
        for axis in 2:3
            delta=samples[i+1][axis]-samples[i][axis]
            squared+=delta*delta
        end
        parameters[i+1]=parameters[i]+sqrt(squared)
    end
    if periodic
        delta=samples[count][1]-samples[1][1]
        squared=delta*delta
        for axis in 2:3
            delta=samples[count][axis]-samples[1][axis]
            squared+=delta*delta
        end
        parameters[count+1]=parameters[count]+sqrt(squared)
    end
    flatknots=Float64[]
    degree=0
    controls=NTuple{3,Float64}[]
    if !periodic
        if count==2
            degree=1
            append!(flatknots,(parameters[1],parameters[1],
                               parameters[2],parameters[2]))
            append!(controls,samples)
        elseif count==3
            degree=2
            append!(flatknots,(parameters[1],parameters[1],parameters[1],
                               parameters[3],parameters[3],parameters[3]))
            controls=_model_occ_spline_solve(flatknots,degree,parameters,
                                             fill(0,count),samples,caller)
        else
            degree=3
            resize!(flatknots,count+6)
            flatknots[1]=flatknots[2]=flatknots[3]=flatknots[4]=parameters[1]
            for i in 2:count-1
                flatknots[3+i]=parameters[i]
            end
            flatknots[count+3]=flatknots[count+4]=
                flatknots[count+5]=flatknots[count+6]=parameters[count]
            first_tangent=_model_occ_lagrange_tangent(
                samples[1:4],parameters[1:4],3,parameters[1])
            last_tangent=_model_occ_lagrange_tangent(
                samples[count-3:count],parameters[count-3:count],3,parameters[count])
            taus=Vector{Float64}(undef,count+2)
            orders=Vector{Int}(undef,count+2)
            values=Vector{NTuple{3,Float64}}(undef,count+2)
            taus[1]=parameters[1];orders[1]=0;values[1]=samples[1]
            taus[2]=parameters[1];orders[2]=1;values[2]=first_tangent
            for i in 2:count-1
                taus[i+1]=parameters[i];orders[i+1]=0;values[i+1]=samples[i]
            end
            taus[count+1]=parameters[count];orders[count+1]=1
            values[count+1]=last_tangent
            taus[count+2]=parameters[count];orders[count+2]=0
            values[count+2]=samples[count]
            controls=_model_occ_spline_solve(flatknots,degree,taus,orders,
                                             values,caller)
        end
    else
        period=parameters[count+1]
        if count==2
            degree=1
            append!(flatknots,(parameters[1],parameters[1],parameters[2],
                               period,period))
            append!(controls,(samples[1],samples[2],samples[1]))
        else
            degree=3
            shift=parameters[count]-period
            extension=parameters[2]+period
            resize!(flatknots,count+7)
            flatknots[1]=flatknots[2]=shift
            flatknots[3]=flatknots[4]=parameters[1]
            for i in 2:count
                flatknots[3+i]=parameters[i]
            end
            flatknots[count+4]=flatknots[count+5]=period
            flatknots[count+6]=flatknots[count+7]=extension
            lagrange_degree=count==3 ? 2 : 3
            tangent=_model_occ_lagrange_tangent(
                samples[1:lagrange_degree+1],parameters[1:lagrange_degree+1],
                lagrange_degree,parameters[1])
            taus=Vector{Float64}(undef,count+3)
            orders=Vector{Int}(undef,count+3)
            values=Vector{NTuple{3,Float64}}(undef,count+3)
            taus[1]=parameters[1];orders[1]=0;values[1]=samples[1]
            taus[2]=parameters[1];orders[2]=1;values[2]=tangent
            for i in 2:count
                taus[i+1]=parameters[i];orders[i+1]=0;values[i+1]=samples[i]
            end
            taus[count+2]=period;orders[count+2]=1;values[count+2]=tangent
            taus[count+3]=period;orders[count+3]=0;values[count+3]=samples[1]
            controls=_model_occ_spline_solve(flatknots,degree,taus,orders,
                                             values,caller)
        end
    end
    curve=_ModelOCCNURBSCurve(:bspline,
        NURBS.NURBSCurve(degree,flatknots,controls,ones(Float64,length(controls))))
    lower,upper=_model_occ_curve_bounds(curve)
    geometry=(occ=:copied_curve,curve=curve,t0=lower,t1=upper)
    shape_id=only(_model_occ_reserve_shapes(kernel,1,caller))
    endpoints=UInt64[first(ids),last(ids)]
    _model_occ_bind_shape!(kernel,1,allocated,shape_id,
        _ModelOCCShape(1,(kind=:tensor_bspline,geometry=geometry),endpoints))
    kernel.next_shape=shape_id+UInt64(1)
    return allocated
end
