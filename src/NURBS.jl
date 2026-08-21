"""
    NURBS

Native B-spline / NURBS curve and surface evaluation (De Boor / Cox–de Boor).
This is the OCC-equivalent query path for parametric CAD. Classified STEP/IGES
solids are imported by [`BRep`](@ref); unrecognized topology remains an explicit
blocker.
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
    length(raw)==3 || throw(ArgumentError("$caller: control point must have 3 coordinates"))
    return (_finite(raw[1],caller,"x"),_finite(raw[2],caller,"y"),_finite(raw[3],caller,"z"))
end

"""Open NURBS curve of degree `p` with `n` weighted controls and `n+p+1` knots."""
struct NURBSCurve
    degree::Int
    knots::Vector{Float64}
    controls::Vector{NTuple{3,Float64}}
    weights::Vector{Float64}
    function NURBSCurve(degree, knots, controls, weights=nothing)
        caller="NURBSCurve"
        p=Int(degree)
        p>=1 || throw(ArgumentError("$caller: degree must be ≥ 1"))
        U=Float64[_finite(k,caller,"knot") for k in knots]
        issorted(U) || throw(ArgumentError("$caller: knots must be nondecreasing"))
        P=[_point3(c,caller) for c in controls]
        n=length(P)
        n>=p+1 || throw(ArgumentError("$caller: need at least degree+1 controls"))
        length(U)==n+p+1 || throw(ArgumentError(
            "$caller: knot vector length must be n+p+1 ($(n+p+1)), got $(length(U))"))
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
        pu,pv=Int(degree_u),Int(degree_v)
        (pu>=1 && pv>=1) || throw(ArgumentError("$caller: degrees must be ≥ 1"))
        U=Float64[_finite(k,caller,"knot_u") for k in knots_u]
        V=Float64[_finite(k,caller,"knot_v") for k in knots_v]
        (issorted(U) && issorted(V)) || throw(ArgumentError("$caller: knots must be nondecreasing"))
        C=controls isa AbstractMatrix ? controls :
            throw(ArgumentError("$caller: controls must be an n_u × n_v matrix"))
        nu,nv=size(C)
        nu>=pu+1 && nv>=pv+1 || throw(ArgumentError("$caller: too few controls"))
        length(U)==nu+pu+1 || throw(ArgumentError("$caller: U knot length must be nu+pu+1"))
        length(V)==nv+pv+1 || throw(ArgumentError("$caller: V knot length must be nv+pv+1"))
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
    h[4]>0 || throw(ArgumentError("nurbs_eval: non-positive weight at parameter $u"))
    return (h[1]/h[4], h[2]/h[4], h[3]/h[4])
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
        den==0 && throw(ArgumentError("nurbs_eval: repeated knot span has zero width"))
        α=(u-ulo)/den
        a=d[j]; b=d[j+1]
        d[j+1]=((1-α)*a[1]+α*b[1], (1-α)*a[2]+α*b[2],
                (1-α)*a[3]+α*b[3], (1-α)*a[4]+α*b[4])
    end
    return d[p+1]
end

function nurbs_eval(c::NURBSCurve, u)
    uu=_finite(u,"nurbs_eval","u")
    U=c.knots
    uu<U[c.degree+1] && (uu=U[c.degree+1])
    uu>U[end-c.degree] && (uu=U[end-c.degree])
    Pw=NTuple{4,Float64}[(c.weights[i]*c.controls[i][1],
                          c.weights[i]*c.controls[i][2],
                          c.weights[i]*c.controls[i][3],
                          c.weights[i]) for i in eachindex(c.controls)]
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
    @inbounds for i in 1:nu
        Pw=NTuple{4,Float64}[(s.weights[i,j]*s.controls[i,j][1],
                              s.weights[i,j]*s.controls[i,j][2],
                              s.weights[i,j]*s.controls[i,j][3],
                              s.weights[i,j]) for j in 1:size(s.controls,2)]
        # deboor returns cartesian; keep homogeneous by evaluating as a curve of 4D
        tmp[i]=_deboor_h(pv,s.knots_v,Pw,vv)
    end
    p=_deboor(pu,s.knots_u,tmp,uu)
    return p
end

# Cox–de Boor scalar basis N_{i,p}(u) for independent-oracle tests.
function bspline_basis(i::Int, p::Int, U::Vector{Float64}, u::Float64)
    # i is 1-based control index
    n=length(U)-p-1
    1<=i<=n || throw(ArgumentError("bspline_basis: control index out of range"))
    u==U[end] && i==n && return 1.0
    if p==0
        return (U[i]<=u<U[i+1] || (u==U[end] && i==n)) ? 1.0 : 0.0
    end
    left=U[i+p]==U[i] ? 0.0 : (u-U[i])/(U[i+p]-U[i])*bspline_basis(i,p-1,U,u)
    right=U[i+p+1]==U[i+1] ? 0.0 : (U[i+p+1]-u)/(U[i+p+1]-U[i+1])*bspline_basis(i+1,p-1,U,u)
    return left+right
end

end # module
