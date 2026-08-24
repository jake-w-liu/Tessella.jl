"""
    NURBS

Native B-spline / NURBS curve and surface evaluation (De Boor / Cox–de Boor).
This is the OCC-equivalent query path for parametric CAD. STEP/IGES NURBS
entities and classified solids are imported by [`BRep`](@ref); unrecognized
topology remains an explicit blocker.
"""
module NURBS

export NURBSCurve, NURBSSurface, nurbs_eval, bspline_basis

function _finite(value, caller, name)
    v=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    isfinite(v) || throw(ArgumentError("$caller: $name must be finite"))
    return v
end

function _point3(raw, caller)
    values=try
        Tuple(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: control point must be an iterable with 3 coordinates"))
    end
    length(values)==3 || throw(ArgumentError("$caller: control point must have 3 coordinates"))
    return (_finite(values[1],caller,"x"),_finite(values[2],caller,"y"),
            _finite(values[3],caller,"z"))
end

function _degree(value, caller, name; allow_zero=false)
    value isa Integer || throw(ArgumentError("$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    minimum=allow_zero ? 0 : 1
    value>=minimum || throw(ArgumentError("$caller: $name must be ≥ $minimum"))
    value<=typemax(Int) || throw(ArgumentError("$caller: $name exceeds the platform Int range"))
    return Int(value)
end

function _validate_knots(U::Vector{Float64}, p::Int, ncontrols::Int, caller, name)
    ncontrols>p || throw(ArgumentError("$caller: need at least degree+1 controls in $name"))
    expected=try
        Base.checked_add(Base.checked_add(ncontrols,p),1)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: knot count overflows the platform Int range"))
    end
    length(U)==expected || throw(ArgumentError(
        "$caller: $name knot vector length must be n+p+1 ($expected), got $(length(U))"))
    issorted(U) || throw(ArgumentError("$caller: $name knots must be nondecreasing"))
    U[p+1]<U[end-p] || throw(ArgumentError(
        "$caller: $name active knot interval must have positive width"))
    return nothing
end

"""Open NURBS curve of degree `p` with `n` weighted controls and `n+p+1` knots."""
struct NURBSCurve
    degree::Int
    knots::Vector{Float64}
    controls::Vector{NTuple{3,Float64}}
    weights::Vector{Float64}
    function NURBSCurve(degree, knots, controls, weights=nothing)
        caller="NURBSCurve"
        p=_degree(degree,caller,"degree")
        U=Float64[_finite(k,caller,"knot") for k in knots]
        P=[_point3(c,caller) for c in controls]
        n=length(P)
        _validate_knots(U,p,n,caller,"curve")
        W=if weights===nothing
            ones(Float64,n)
        else
            Float64[_finite(w,caller,"weight") for w in weights]
        end
        length(W)==n || throw(ArgumentError("$caller: weight count mismatch"))
        all(>(0),W) || throw(ArgumentError("$caller: weights must be positive"))
        return new(p,U,P,W)
    end
end

"""Tensor-product NURBS surface."""
struct NURBSSurface
    degree_u::Int
    degree_v::Int
    knots_u::Vector{Float64}
    knots_v::Vector{Float64}
    controls::Matrix{NTuple{3,Float64}}   # n_u × n_v
    weights::Matrix{Float64}
    function NURBSSurface(degree_u, degree_v, knots_u, knots_v, controls, weights=nothing)
        caller="NURBSSurface"
        pu=_degree(degree_u,caller,"degree_u")
        pv=_degree(degree_v,caller,"degree_v")
        U=Float64[_finite(k,caller,"knot_u") for k in knots_u]
        V=Float64[_finite(k,caller,"knot_v") for k in knots_v]
        C=controls isa AbstractMatrix ? controls :
            throw(ArgumentError("$caller: controls must be an n_u × n_v matrix"))
        nu,nv=size(C)
        _validate_knots(U,pu,nu,caller,"u")
        _validate_knots(V,pv,nv,caller,"v")
        P=Matrix{NTuple{3,Float64}}(undef,nu,nv)
        @inbounds for j in 1:nv, i in 1:nu
            P[i,j]=_point3(C[i,j],caller)
        end
        W=if weights===nothing
            ones(Float64,nu,nv)
        else
            size(weights)==(nu,nv) || throw(ArgumentError("$caller: weight size mismatch"))
            Float64[_finite(weights[i,j],caller,"weight") for i in 1:nu, j in 1:nv]
        end
        all(>(0),W) || throw(ArgumentError("$caller: weights must be positive"))
        return new(pu,pv,U,V,P,W)
    end
end

function _span(p, U, u)
    n=length(U)-p-2  # last control index (0-based n)
    u>=U[end-p] && return length(U)-p-2
    low=p; high=length(U)-p-1
    mid=(low+high)÷2
    while u<U[mid+1] || u>=U[mid+2]
        u<U[mid+1] ? (high=mid) : (low=mid)
        mid=(low+high)÷2
    end
    return mid
end

function _deboor(p, U, Pw, u)
    h=_deboor_h(p,U,Pw,u)
    (isfinite(h[4]) && h[4]>0) || throw(ArgumentError(
        "nurbs_eval: non-positive or unrepresentable weight at parameter $u"))
    result=(h[1]/h[4], h[2]/h[4], h[3]/h[4])
    all(isfinite,result) || throw(ArgumentError(
        "nurbs_eval: evaluated coordinates are not Float64-representable"))
    return result
end

function _deboor_h(p, U, Pw, u)
    i=_span(p,U,u)
    d=Vector{NTuple{4,Float64}}(undef,p+1)
    @inbounds for j in 0:p
        d[j+1]=Pw[i-p+j+1]
    end
    @inbounds for r in 1:p, j in p:-1:r
        ulo=U[i-p+j+1]; uhi=U[i+j-r+2]
        den=uhi-ulo
        (isfinite(den) && den>0) || throw(ArgumentError(
            "nurbs_eval: knot span width is zero or not Float64-representable"))
        α=(u-ulo)/den
        isfinite(α) || throw(ArgumentError(
            "nurbs_eval: knot interpolation factor is not Float64-representable"))
        a=d[j]; b=d[j+1]
        blended=((1-α)*a[1]+α*b[1], (1-α)*a[2]+α*b[2],
                 (1-α)*a[3]+α*b[3], (1-α)*a[4]+α*b[4])
        all(isfinite,blended) || throw(ArgumentError(
            "nurbs_eval: homogeneous coordinates are not Float64-representable"))
        d[j+1]=blended
    end
    return d[p+1]
end

"""
    nurbs_eval(curve, u) -> (x, y, z)
    nurbs_eval(surface, u, v) -> (x, y, z)

Evaluate a [`NURBSCurve`](@ref) or [`NURBSSurface`](@ref) with the De Boor
algorithm. Finite parameters outside the active knot interval are clamped to
the nearest endpoint. An `ArgumentError` is raised when the homogeneous result
cannot be represented by finite `Float64` coordinates.
"""
function nurbs_eval(c::NURBSCurve, u)
    uu=_finite(u,"nurbs_eval","u")
    U=c.knots
    uu<U[c.degree+1] && (uu=U[c.degree+1])
    uu>U[end-c.degree] && (uu=U[end-c.degree])
    scale=maximum(c.weights)
    Pw=map(eachindex(c.controls)) do i
        w=c.weights[i]/scale
        w>0 || throw(ArgumentError(
            "nurbs_eval: relative control weight is not Float64-representable"))
        (w*c.controls[i][1],w*c.controls[i][2],w*c.controls[i][3],w)
    end
    return _deboor(c.degree,c.knots,Pw,uu)
end

function nurbs_eval(s::NURBSSurface, u, v)
    uu=_finite(u,"nurbs_eval","u"); vv=_finite(v,"nurbs_eval","v")
    return _nurbs_surface_eval(s,uu,vv)
end

function _nurbs_surface_eval(s::NURBSSurface, u, v)
    pu,pv=s.degree_u,s.degree_v
    uu=clamp(u,s.knots_u[pu+1],s.knots_u[end-pu])
    vv=clamp(v,s.knots_v[pv+1],s.knots_v[end-pv])
    nu=size(s.controls,1)
    tmp=Vector{NTuple{4,Float64}}(undef,nu)
    scale=maximum(s.weights)
    @inbounds for i in 1:nu
        Pw=map(1:size(s.controls,2)) do j
            w=s.weights[i,j]/scale
            w>0 || throw(ArgumentError(
                "nurbs_eval: relative control weight is not Float64-representable"))
            (w*s.controls[i,j][1],w*s.controls[i,j][2],w*s.controls[i,j][3],w)
        end
        # deboor returns cartesian; keep homogeneous by evaluating as a curve of 4D
        tmp[i]=_deboor_h(pv,s.knots_v,Pw,vv)
    end
    p=_deboor(pu,s.knots_u,tmp,uu)
    return p
end

"""
    bspline_basis(i, degree, knots, u) -> Float64

Evaluate the 1-based Cox–de Boor basis function `N[i,degree](u)`. The knot
vector must be finite, nondecreasing, and define a positive-width active
interval with at least `degree+1` basis functions. Parameters outside that
active interval return zero. Invalid indices, degrees, knots, and nonfinite
parameters raise `ArgumentError`.
"""
function bspline_basis(i, degree, knots, u)
    caller="bspline_basis"
    p=_degree(degree,caller,"degree"; allow_zero=true)
    U=Float64[_finite(k,caller,"knot") for k in knots]
    p<=length(U)-2 || throw(ArgumentError(
        "$caller: knot vector is too short for degree $p"))
    n=length(U)-p-1
    n>p || throw(ArgumentError(
        "$caller: need at least degree+1 basis functions"))
    _validate_knots(U,p,n,caller,"basis")
    i isa Integer || throw(ArgumentError("$caller: control index must be an integer"))
    i isa Bool && throw(ArgumentError("$caller: control index must not be Bool"))
    (1<=i<=n) || throw(ArgumentError("$caller: control index out of range"))
    ii=Int(i)
    uu=_finite(u,caller,"u")
    active_lo,active_hi=U[p+1],U[end-p]
    (active_lo<=uu<=active_hi) || return 0.0
    if uu==active_hi && active_hi==U[end]
        return ii==n ? 1.0 : 0.0
    end

    work=Vector{Float64}(undef,p+1)
    @inbounds for offset in 0:p
        k=ii+offset
        work[offset+1]=(U[k]<=uu<U[k+1] ||
                        (uu==U[end] && k==length(U)-1)) ? 1.0 : 0.0
    end
    @inbounds for d in 1:p, offset in 0:p-d
        k=ii+offset
        left_den=U[k+d]-U[k]
        right_den=U[k+d+1]-U[k+1]
        (isfinite(left_den) && isfinite(right_den)) || throw(ArgumentError(
            "$caller: knot span width is not Float64-representable"))
        left=left_den==0 ? 0.0 : (uu-U[k])/left_den*work[offset+1]
        right=right_den==0 ? 0.0 : (U[k+d+1]-uu)/right_den*work[offset+2]
        value=left+right
        isfinite(value) || throw(ArgumentError(
            "$caller: basis value is not Float64-representable"))
        work[offset+1]=value
    end
    return work[1]
end

end # module
