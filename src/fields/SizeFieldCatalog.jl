# Extended Gmsh 4.15.2 field catalog. Included from SizeField after the core
# Distance/Threshold/Box/Ball/Cylinder/Frustum/Min/Max/Bounded types. Some
# geometry- or view-backed fields below intentionally expose native discrete
# constructors; the .geo builder rejects them until their required model context
# is available instead of pretending that a coordinate-only substitute is parity.

using ..Predicates: orient2, orient3

# ── Anisotropic metric ────────────────────────────────────────────────────────

"""Symmetric 3×3 metric `xᵀ M x`. Isotropic size `h` corresponds to `M = I / h²`."""
struct Metric3
    m11::Float64; m22::Float64; m33::Float64
    m12::Float64; m13::Float64; m23::Float64
    function Metric3(m11::Real,m22::Real,m33::Real,m12::Real,m13::Real,m23::Real)
        any(value -> value isa Bool,(m11,m22,m33,m12,m13,m23)) &&
            throw(ArgumentError("Metric3: entries must not be Bool"))
        values = try
            (Float64(m11),Float64(m22),Float64(m33),
             Float64(m12),Float64(m13),Float64(m23))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("Metric3: entries must be Float64-representable"))
        end
        all(isfinite,values) || throw(ArgumentError("Metric3: entries must be finite"))
        _metric_positive_definite(values...) || throw(ArgumentError(
            "Metric3: matrix must be symmetric positive definite"))
        return new(values...)
    end
end

"""Return the isotropic metric `I / h²` for a finite, representable size `h`."""
function isotropic_metric(h::Real)
    hh=_positive_value(h,"isotropic_metric","h")
    ih=inv(hh); s=ih*ih
    (isfinite(s) && s>0) || throw(ArgumentError(
        "isotropic_metric: h=$hh is outside the representable metric range"))
    return Metric3(s,s,s,0.0,0.0,0.0)
end

@inline function _metric_quad(m::Metric3,x,y,z)
    return m.m11*x*x + m.m22*y*y + m.m33*z*z + 2*(m.m12*x*y + m.m13*x*z + m.m23*y*z)
end

@inline function _jacobi_rotation(app,aqq,apq)
    threshold=16eps(Float64)*sqrt(abs(app))*sqrt(abs(aqq))
    abs(apq)<=threshold && return (false,1.0,0.0)
    tau=(aqq-app)/(2*apq)
    t=tau>=0 ? 1/(tau+hypot(1,tau)) : -1/(-tau+hypot(1,tau))
    cosine=inv(sqrt(1+t*t))
    return (true,cosine,t*cosine)
end

@inline function _jacobi12(A::NTuple{6,Float64},V::NTuple{9,Float64})
    a11,a22,a33,a12,a13,a23=A
    changed,c,s=_jacobi_rotation(a11,a22,a12)
    changed || return A,V,false
    matrix=(c*c*a11-2s*c*a12+s*s*a22,
            s*s*a11+2s*c*a12+c*c*a22,a33,0.0,
            c*a13-s*a23,s*a13+c*a23)
    vectors=(c*V[1]-s*V[2],s*V[1]+c*V[2],V[3],
             c*V[4]-s*V[5],s*V[4]+c*V[5],V[6],
             c*V[7]-s*V[8],s*V[7]+c*V[8],V[9])
    return matrix,vectors,true
end

@inline function _jacobi13(A::NTuple{6,Float64},V::NTuple{9,Float64})
    a11,a22,a33,a12,a13,a23=A
    changed,c,s=_jacobi_rotation(a11,a33,a13)
    changed || return A,V,false
    matrix=(c*c*a11-2s*c*a13+s*s*a33,a22,
            s*s*a11+2s*c*a13+c*c*a33,
            c*a12-s*a23,0.0,s*a12+c*a23)
    vectors=(c*V[1]-s*V[3],V[2],s*V[1]+c*V[3],
             c*V[4]-s*V[6],V[5],s*V[4]+c*V[6],
             c*V[7]-s*V[9],V[8],s*V[7]+c*V[9])
    return matrix,vectors,true
end

@inline function _jacobi23(A::NTuple{6,Float64},V::NTuple{9,Float64})
    a11,a22,a33,a12,a13,a23=A
    changed,c,s=_jacobi_rotation(a22,a33,a23)
    changed || return A,V,false
    matrix=(a11,c*c*a22-2s*c*a23+s*s*a33,
            s*s*a22+2s*c*a23+c*c*a33,
            c*a12-s*a13,s*a12+c*a13,0.0)
    vectors=(V[1],c*V[2]-s*V[3],s*V[2]+c*V[3],
             V[4],c*V[5]-s*V[6],s*V[5]+c*V[6],
             V[7],c*V[8]-s*V[9],s*V[8]+c*V[9])
    return matrix,vectors,true
end

@inline _jacobi_sweeps(A,V,::Val{0})=(A,V)
@inline function _jacobi_sweeps(A,V,::Val{N}) where {N}
    A,V,c12=_jacobi12(A,V)
    A,V,c13=_jacobi13(A,V)
    A,V,c23=_jacobi23(A,V)
    (c12||c13||c23) || return A,V
    return _jacobi_sweeps(A,V,Val(N-1))
end

@inline function _sym3_eigh(a11,a22,a33,a12,a13,a23)
    if a12==0 && a13==0 && a23==0
        return ((Float64(a11),Float64(a22),Float64(a33)),
                ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)))
    end
    scale=max(abs(a11),abs(a22),abs(a33),abs(a12),abs(a13),abs(a23))
    scale==0 && return ((0.0,0.0,0.0),
                        ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)))
    isfinite(scale) || throw(ArgumentError("metric eigensystem: non-finite matrix"))
    A=(Float64(a11/scale),Float64(a22/scale),Float64(a33/scale),
       Float64(a12/scale),Float64(a13/scale),Float64(a23/scale))
    V=(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)
    A,V=_jacobi_sweeps(A,V,Val(12))
    values=(scale*A[1],scale*A[2],scale*A[3])
    spectral=max(abs(values[1]),abs(values[2]),abs(values[3]))
    if spectral>0 && minimum(abs,values)<=4096eps(Float64)*spectral
        return _sym3_eigh_high_precision(a11,a22,a33,a12,a13,a23)
    end
    vectors=((V[1],V[4],V[7]),(V[2],V[5],V[8]),(V[3],V[6],V[9]))
    return values,vectors
end

function _sym3_eigh_high_precision(a11,a22,a33,a12,a13,a23)
    return setprecision(BigFloat,256) do
        A=BigFloat[BigFloat(a11) BigFloat(a12) BigFloat(a13);
                   BigFloat(a12) BigFloat(a22) BigFloat(a23);
                   BigFloat(a13) BigFloat(a23) BigFloat(a33)]
        V=fill(BigFloat(0),3,3)
        V[1,1]=V[2,2]=V[3,3]=BigFloat(1)
        tolerance=ldexp(BigFloat(1),-220)
        for _ in 1:80
            changed=false
            for (p,q) in ((1,2),(1,3),(2,3))
                app=A[p,p]; aqq=A[q,q]; apq=A[p,q]
                abs(apq)<=tolerance*sqrt(abs(app))*sqrt(abs(aqq)) && continue
                tau=(aqq-app)/(2apq)
                t=tau>=0 ? inv(tau+hypot(BigFloat(1),tau)) :
                           -inv(-tau+hypot(BigFloat(1),tau))
                cosine=inv(sqrt(1+t*t)); sine=t*cosine
                for k in 1:3
                    (k==p || k==q) && continue
                    akp=A[k,p]; akq=A[k,q]
                    A[k,p]=A[p,k]=cosine*akp-sine*akq
                    A[k,q]=A[q,k]=sine*akp+cosine*akq
                end
                A[p,p]=cosine*cosine*app-2sine*cosine*apq+sine*sine*aqq
                A[q,q]=sine*sine*app+2sine*cosine*apq+cosine*cosine*aqq
                A[p,q]=A[q,p]=0
                for k in 1:3
                    vkp=V[k,p]; vkq=V[k,q]
                    V[k,p]=cosine*vkp-sine*vkq
                    V[k,q]=sine*vkp+cosine*vkq
                end
                changed=true
            end
            changed || break
        end
        tofloat(value::BigFloat)=begin
            result=Float64(value)
            result==0 && !iszero(value) ?
                (signbit(value) ? -nextfloat(0.0) : nextfloat(0.0)) : result
        end
        values=(tofloat(A[1,1]),tofloat(A[2,2]),tofloat(A[3,3]))
        vectors=((Float64(V[1,1]),Float64(V[2,1]),Float64(V[3,1])),
                 (Float64(V[1,2]),Float64(V[2,2]),Float64(V[3,2])),
                 (Float64(V[1,3]),Float64(V[2,3]),Float64(V[3,3])))
        return values,vectors
    end
end

function _metric_cholesky_components(a,b,c,d,e,f)::Union{Nothing,NTuple{6,Float64}}
    all(isfinite,(a,b,c,d,e,f)) || return nothing
    (a>0 && b>0 && c>0) || return nothing
    l11=sqrt(a); l21=d/l11; l31=e/l11
    l21sq=l21*l21
    r22=b-l21sq
    # Once a Schur complement loses roughly half the Float64 significand, its
    # square root is not accurate enough for directional length certificates.
    # Reconstruct those uncommon, ill-conditioned factors from exact dyadics.
    guard=sqrt(eps(Float64))
    uncertain22=max(abs(b),abs(l21sq))>0 &&
                abs(r22)<=guard*max(abs(b),abs(l21sq))
    if r22>0 && !uncertain22
        product=l31*l21; numerator=f-product
        uncertain32=max(abs(f),abs(product))>0 &&
                    abs(numerator)<=guard*max(abs(f),abs(product))
        l22=sqrt(r22); l32=numerator/l22
        l31sq=l31*l31; l32sq=l32*l32
        r33=c-l31sq-l32sq
        uncertain33=max(abs(c),abs(l31sq)+abs(l32sq))>0 &&
                    abs(r33)<=guard*max(abs(c),abs(l31sq)+abs(l32sq))
        if r33>0 && !uncertain32 && !uncertain33
            result=(l11,l21,l22,l31,l32,sqrt(r33))
            all(isfinite,result) && return result
        end
    end
    return _metric_positive_definite(a,b,c,d,e,f) ?
           _metric_cholesky_exact(a,b,c,d,e,f) : nothing
end

@inline function _sort3(v::NTuple{3,Float64})
    a,b,c=v
    a>b && ((a,b)=(b,a)); b>c && ((b,c)=(c,b)); a>b && ((a,b)=(b,a))
    return (a,b,c)
end
"""Return the three eigenvalues of `m` in ascending order."""
function metric_eigenvalues(m::Metric3)
    λ,_=_sym3_eigh(m.m11,m.m22,m.m33,m.m12,m.m13,m.m23)
    return _sort3(λ)
end

"""Return the smallest directional size `1 / √λₘₐₓ` represented by `m`."""
function metric_size(m::Metric3)
    λ=metric_eigenvalues(m)
    minimum(λ)>0 || throw(ArgumentError("metric_size: metric is not positive definite (λ=$λ)"))
    h=inv(sqrt(maximum(λ)))
    (isfinite(h) && h>0) || throw(ArgumentError("metric_size: size is not representable (λ=$λ)"))
    return h
end

function metric_from_axes(λ1,λ2,λ3,v1,v2,v3)
    # M = Σ λ_i v_i v_iᵀ with unit axes
    n1=_norm3(v1); n2=_norm3(v2); n3=_norm3(v3)
    (n1>0 && n2>0 && n3>0) || throw(ArgumentError("metric_from_axes: axes must have positive length"))
    e1=ntuple(i->v1[i]/n1,3); e2=ntuple(i->v2[i]/n2,3); e3=ntuple(i->v3[i]/n3,3)
    tol=1024eps(Float64)
    (abs(_dot3(e1,e2))<=tol && abs(_dot3(e1,e3))<=tol && abs(_dot3(e2,e3))<=tol) ||
        throw(ArgumentError("metric_from_axes: axes must be mutually orthogonal"))
    a=_positive_value(λ1,"metric_from_axes","λ1")
    b=_positive_value(λ2,"metric_from_axes","λ2")
    c=_positive_value(λ3,"metric_from_axes","λ3")
    m11=a*e1[1]*e1[1]+b*e2[1]*e2[1]+c*e3[1]*e3[1]
    m22=a*e1[2]*e1[2]+b*e2[2]*e2[2]+c*e3[2]*e3[2]
    m33=a*e1[3]*e1[3]+b*e2[3]*e2[3]+c*e3[3]*e3[3]
    m12=a*e1[1]*e1[2]+b*e2[1]*e2[2]+c*e3[1]*e3[2]
    m13=a*e1[1]*e1[3]+b*e2[1]*e2[3]+c*e3[1]*e3[3]
    m23=a*e1[2]*e1[3]+b*e2[2]*e2[3]+c*e3[2]*e3[3]
    return Metric3(m11,m22,m33,m12,m13,m23)
end

function _metric_from_raw_axes(λ1::Float64,λ2::Float64,λ3::Float64,
                               v1::NTuple{3,Float64},v2::NTuple{3,Float64},
                               v3::NTuple{3,Float64})::Metric3
    m11=λ1*v1[1]^2+λ2*v2[1]^2+λ3*v3[1]^2
    m22=λ1*v1[2]^2+λ2*v2[2]^2+λ3*v3[2]^2
    m33=λ1*v1[3]^2+λ2*v2[3]^2+λ3*v3[3]^2
    m12=λ1*v1[1]*v1[2]+λ2*v2[1]*v2[2]+λ3*v3[1]*v3[2]
    m13=λ1*v1[1]*v1[3]+λ2*v2[1]*v2[3]+λ3*v3[1]*v3[3]
    m23=λ1*v1[2]*v1[3]+λ2*v2[2]*v2[3]+λ3*v3[2]*v3[3]
    return Metric3(m11,m22,m33,m12,m13,m23)
end

@inline _metric_matrix(m::Metric3) =
    (m.m11,m.m12,m.m13, m.m12,m.m22,m.m23, m.m13,m.m23,m.m33)
@inline _transpose3(a::NTuple{9,Float64}) =
    (a[1],a[4],a[7], a[2],a[5],a[8], a[3],a[6],a[9])
@inline function _mul3(a::NTuple{9,Float64},b::NTuple{9,Float64})
    return (a[1]*b[1]+a[2]*b[4]+a[3]*b[7],
            a[1]*b[2]+a[2]*b[5]+a[3]*b[8],
            a[1]*b[3]+a[2]*b[6]+a[3]*b[9],
            a[4]*b[1]+a[5]*b[4]+a[6]*b[7],
            a[4]*b[2]+a[5]*b[5]+a[6]*b[8],
            a[4]*b[3]+a[5]*b[6]+a[6]*b[9],
            a[7]*b[1]+a[8]*b[4]+a[9]*b[7],
            a[7]*b[2]+a[8]*b[5]+a[9]*b[8],
            a[7]*b[3]+a[8]*b[6]+a[9]*b[9])
end
@inline function _symmetric_metric(a::NTuple{9,Float64})
    return Metric3(a[1],a[5],a[9],(a[2]+a[4])/2,(a[3]+a[7])/2,(a[6]+a[8])/2)
end

"""
    intersection_alauzet(m1, m2) -> Metric3

Return the symmetric-positive-definite metric intersection obtained by a
generalized eigendecomposition. The result is independently certifiable in the
Loewner order: it dominates both input metrics in every direction.
"""
@inline function intersection_alauzet(m1::Metric3, m2::Metric3)
    m1==m2 && return m1
    # Avoid rebuilding a metric only after a guarded or exact certificate. A
    # bare floating generalized ratio of one can be a false positive for a
    # nearly singular residual at the top of the Float64 range.
    _metric_dominates(m1,m2) && return m1
    _metric_dominates(m2,m1) && return m2
    diagonal1=(m1.m12==0 && m1.m13==0 && m1.m23==0)
    diagonal2=(m2.m12==0 && m2.m13==0 && m2.m23==0)
    if diagonal1 && diagonal2
        return Metric3(max(m1.m11,m2.m11),max(m1.m22,m2.m22),max(m1.m33,m2.m33),
                       0.0,0.0,0.0)
    end
    # Normalize both metrics by the same exact power of two. This prevents the
    # Cholesky/generalized-eigen reconstruction from overflowing when the
    # representable answer lies near `floatmax`.
    scale=max(abs(m1.m11),abs(m1.m22),abs(m1.m33),abs(m1.m12),abs(m1.m13),abs(m1.m23),
              abs(m2.m11),abs(m2.m22),abs(m2.m33),abs(m2.m12),abs(m2.m13),abs(m2.m23))
    exponent2=exponent(scale)
    a=_metric_ldexp(m1,-exponent2)
    b=_metric_ldexp(m2,-exponent2)
    if a===nothing || b===nothing
        return _metric_boundary_upper(nothing,m1,m2)
    end
    reduced_result=try
        _intersection_alauzet_core(a::Metric3,b::Metric3)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return _metric_boundary_upper(nothing,m1,m2)
    end
    result=_metric_ldexp(reduced_result,exponent2)
    result===nothing && return _metric_boundary_upper(nothing,m1,m2)
    return _inflate_to_dominate(result::Metric3,m1,m2)
end

@inline function _intersection_alauzet_core(m1::Metric3,m2::Metric3)
    l11,l21,l22,l31,l32,l33=_metric_cholesky_components(
        m1.m11,m1.m22,m1.m33,m1.m12,m1.m13,m1.m23)::NTuple{6,Float64}
    invL=(inv(l11),0.0,0.0,
          -l21/(l11*l22),inv(l22),0.0,
          (l21*l32-l31*l22)/(l11*l22*l33),-l32/(l22*l33),inv(l33))
    reduced=_mul3(_mul3(invL,_metric_matrix(m2)),_transpose3(invL))
    c=_symmetric_metric(reduced)
    λ,Q=_sym3_eigh(c.m11,c.m22,c.m33,c.m12,c.m13,c.m23)
    q=(Q[1][1],Q[2][1],Q[3][1],
       Q[1][2],Q[2][2],Q[3][2],
       Q[1][3],Q[2][3],Q[3][3])
    d=(max(1.0,λ[1]),0.0,0.0, 0.0,max(1.0,λ[2]),0.0,
       0.0,0.0,max(1.0,λ[3]))
    reduced_intersection=_mul3(_mul3(q,d),_transpose3(q))
    L=(l11,0.0,0.0, l21,l22,0.0, l31,l32,l33)
    return _symmetric_metric(_mul3(_mul3(L,reduced_intersection),_transpose3(L)))
end

function intersection_conserve_m1(m1::Metric3, m2::Metric3)
    _metric_dominates(m1,m2) && return m1
    _metric_dominates(m2,m1) && return m2
    _,V=_sym3_eigh(m1.m11,m1.m22,m1.m33,m1.m12,m1.m13,m1.m23)
    v1,v2,v3=V
    n1=_norm3(v1); n2=_norm3(v2); n3=_norm3(v3)
    (n1>0 && n2>0 && n3>0) || return intersection_alauzet(m1,m2)
    e1=ntuple(i->v1[i]/n1,3); e2=ntuple(i->v2[i]/n2,3); e3=ntuple(i->v3[i]/n3,3)
    l0=max(_metric_quad(m1,e1...),_metric_quad(m2,e1...))
    l1=max(_metric_quad(m1,e2...),_metric_quad(m2,e2...))
    l2=max(_metric_quad(m1,e3...),_metric_quad(m2,e3...))
    candidate=try
        metric_from_axes(l0,l1,l2,e1,e2,e3)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return _metric_boundary_upper(nothing,m1,m2)
    end
    return _inflate_to_dominate(candidate,m1,m2)
end

function intersection_conserve_mostaniso(m1::Metric3, m2::Metric3)
    λ1=metric_eigenvalues(m1); λ2=metric_eigenvalues(m2)
    r1=minimum(λ1)/maximum(λ1); r2=minimum(λ2)/maximum(λ2)
    return r1<r2 ? intersection_conserve_m1(m1,m2) : intersection_conserve_m1(m2,m1)
end

@inline function _metric_dominates(a::Metric3,b::Metric3)
    a==b && return true
    if a.m12==0 && a.m13==0 && a.m23==0 &&
       b.m12==0 && b.m13==0 && b.m23==0
        return a.m11>=b.m11 && a.m22>=b.m22 && a.m33>=b.m33
    end
    # A floating generalized eigenvalue is not a one-sided certificate at high
    # condition number. The outward interval check is allocation-free and can
    # only accept when every exact residual principal minor is nonnegative.
    return _metric_interval_dominates(a,b)
end

@inline function _metric_residual_sdd_guarded(a::Metric3,b::Metric3)
    d1=a.m11-b.m11; d2=a.m22-b.m22; d3=a.m33-b.m33
    u=a.m12-b.m12; v=a.m13-b.m13; w=a.m23-b.m23
    all(isfinite,(d1,d2,d3,u,v,w)) || return false
    factor=1+64sqrt(eps(Float64))
    return _metric_sdd_row(d1,u,v,factor) && _metric_sdd_row(d2,u,w,factor) &&
           _metric_sdd_row(d3,v,w,factor)
end

@inline function _metric_sdd_row(diagonal,x,y,factor)
    diagonal>=0 || return false
    off=abs(x)+abs(y)
    isfinite(off) || return false
    off==0 && return true
    return diagonal>=factor*off
end

@inline function _metric_generalized_min(a::Metric3,b::Metric3)
    l11,l21,l22,l31,l32,l33=_metric_cholesky_components(
        b.m11,b.m22,b.m33,b.m12,b.m13,b.m23)::NTuple{6,Float64}
    invL=(inv(l11),0.0,0.0,
          -l21/(l11*l22),inv(l22),0.0,
          (l21*l32-l31*l22)/(l11*l22*l33),-l32/(l22*l33),inv(l33))
    reduced=_mul3(_mul3(invL,_metric_matrix(a)),_transpose3(invL))
    all(isfinite,reduced) || return -Inf
    values,_=_sym3_eigh(reduced[1],reduced[5],reduced[9],
                        (reduced[2]+reduced[4])/2,
                        (reduced[3]+reduced[7])/2,
                        (reduced[6]+reduced[8])/2)
    return minimum(values)
end

# Gmsh's `intersection` has two source-visible numerical quirks: it stores the
# lower triangle of `inv(m1) * m2` back into an SMetric3 before diagonalizing,
# and it assumes the resulting Euclidean eigenvectors are the metric axes. Keep
# this behavior private to field scalar overloads that require exact 4.15.2
# parity; the public Alauzet helper retains Tessella's certified Loewner-upper
# contract.
@inline function _gmsh_metric_intersection_core(m1::Metric3,m2::Metric3,
                                                proportional_shortcut::Bool)
    chol=_metric_cholesky_components(
        m1.m11,m1.m22,m1.m33,m1.m12,m1.m13,m1.m23)
    chol===nothing && return nothing
    l11,l21,l22,l31,l32,l33=chol::NTuple{6,Float64}
    invL=(inv(l11),0.0,0.0,
          -l21/(l11*l22),inv(l22),0.0,
          (l21*l32-l31*l22)/(l11*l22*l33),-l32/(l22*l33),inv(l33))
    all(isfinite,invL) || return nothing
    inverse=_mul3(_transpose3(invL),invL)
    product=_mul3(inverse,_metric_matrix(m2))
    all(isfinite,product) || return nothing
    # SMetric3::setMat visits rows in order, so the later lower-triangle values
    # overwrite the upper-triangle entries in packed symmetric storage.
    values,vectors=_sym3_eigh(product[1],product[5],product[9],
                              product[4],product[7],product[8])
    all(isfinite,values) || return nothing
    largest=max(values...);smallest=min(values...)
    if proportional_shortcut && largest!=0 && isfinite(largest) &&
       abs((largest-smallest)/largest)<1e-2
        return largest>=1 ? m2 : m1
    end
    weights=ntuple(i->max(_metric_quad(m1,vectors[i]...),
                          _metric_quad(m2,vectors[i]...)),3)
    all(isfinite,weights) || return nothing
    return try
        _metric_from_raw_axes(weights[1],weights[2],weights[3],
                              vectors[1],vectors[2],vectors[3])
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        nothing
    end
end

@inline function _gmsh_metric_intersection_impl(m1::Metric3,m2::Metric3,
                                                proportional_shortcut::Bool)
    scale=max(abs(m1.m11),abs(m1.m22),abs(m1.m33),abs(m1.m12),abs(m1.m13),abs(m1.m23),
              abs(m2.m11),abs(m2.m22),abs(m2.m33),abs(m2.m12),abs(m2.m13),abs(m2.m23))
    exponent2=exponent(scale)
    a=_metric_ldexp(m1,-exponent2);b=_metric_ldexp(m2,-exponent2)
    (a===nothing || b===nothing) && return intersection_alauzet(m1,m2)
    candidate=_gmsh_metric_intersection_core(a::Metric3,b::Metric3,
                                             proportional_shortcut)
    candidate===nothing && return intersection_alauzet(m1,m2)
    result=_metric_ldexp(candidate::Metric3,exponent2)
    return result===nothing ? intersection_alauzet(m1,m2) : result::Metric3
end

@inline _gmsh_metric_intersection(m1::Metric3,m2::Metric3)=
    _gmsh_metric_intersection_impl(m1,m2,true)
@inline _gmsh_metric_intersection_alauzet(m1::Metric3,m2::Metric3)=
    _gmsh_metric_intersection_impl(m1,m2,false)

@inline function _metric_scale(m::Metric3,factor::Float64)
    return Metric3(factor*m.m11,factor*m.m22,factor*m.m33,
                   factor*m.m12,factor*m.m13,factor*m.m23)
end

@inline function _metric_ldexp(m::Metric3,power::Int)
    values=(ldexp(m.m11,power),ldexp(m.m22,power),ldexp(m.m33,power),
            ldexp(m.m12,power),ldexp(m.m13,power),ldexp(m.m23,power))
    all(isfinite,values) || return nothing
    return try
        Metric3(values...)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        nothing
    end
end

@inline _metric_rational(x::Float64)=Rational{BigInt}(x)
@inline _metric_rational(x::Rational{BigInt})=x

const _MetricInterval=NTuple{2,Float64}

@inline function _metric_iv_scaled(x::Float64,power::Int)::_MetricInterval
    value=ldexp(x,power)
    return ldexp(value,-power)==x ? (value,value) :
           (prevfloat(value),nextfloat(value))
end

@inline function _metric_iv_twosum(x::Float64,y::Float64)::_MetricInterval
    value=x+y
    residual=(x-(value-(value-x)))+(y-(value-x))
    return residual>0 ? (value,nextfloat(value)) :
           (residual<0 ? (prevfloat(value),value) : (value,value))
end

@inline function _metric_iv_sub(x::_MetricInterval,y::_MetricInterval)::_MetricInterval
    x[1]==x[2] && y[1]==y[2] && return _metric_iv_twosum(x[1],-y[1])
    return (prevfloat(x[1]-y[2]),nextfloat(x[2]-y[1]))
end
@inline _metric_iv_add(x::_MetricInterval,y::_MetricInterval) =
    (prevfloat(x[1]+y[1]),nextfloat(x[2]+y[2]))
@inline function _metric_iv_mul(x::_MetricInterval,y::_MetricInterval)::_MetricInterval
    p1=x[1]*y[1]; p2=x[1]*y[2]; p3=x[2]*y[1]; p4=x[2]*y[2]
    return (prevfloat(min(p1,p2,p3,p4)),nextfloat(max(p1,p2,p3,p4)))
end
@inline function _metric_iv_square(x::_MetricInterval)::_MetricInterval
    if x[1]<=0<=x[2]
        return (0.0,nextfloat(max(x[1]*x[1],x[2]*x[2])))
    end
    return (prevfloat(min(x[1]*x[1],x[2]*x[2])),
            nextfloat(max(x[1]*x[1],x[2]*x[2])))
end
@inline _metric_iv_minor(x,y,z)=_metric_iv_sub(_metric_iv_mul(x,y),
                                                _metric_iv_square(z))

# Allocation-free one-sided certificate. Every operation encloses its exact
# IEEE-dyadic value after a common power-of-two normalization; accepting the
# lower bounds of all seven principal minors therefore cannot produce a false
# Loewner-dominance result.
@inline function _metric_interval_dominates(candidate::Metric3,input::Metric3)
    scale=max(abs(candidate.m11),abs(candidate.m22),abs(candidate.m33),
              abs(candidate.m12),abs(candidate.m13),abs(candidate.m23),
              abs(input.m11),abs(input.m22),abs(input.m33),
              abs(input.m12),abs(input.m13),abs(input.m23))
    power=-exponent(scale)
    a=_metric_iv_sub(_metric_iv_scaled(candidate.m11,power),
                     _metric_iv_scaled(input.m11,power))
    b=_metric_iv_sub(_metric_iv_scaled(candidate.m22,power),
                     _metric_iv_scaled(input.m22,power))
    c=_metric_iv_sub(_metric_iv_scaled(candidate.m33,power),
                     _metric_iv_scaled(input.m33,power))
    d=_metric_iv_sub(_metric_iv_scaled(candidate.m12,power),
                     _metric_iv_scaled(input.m12,power))
    e=_metric_iv_sub(_metric_iv_scaled(candidate.m13,power),
                     _metric_iv_scaled(input.m13,power))
    f=_metric_iv_sub(_metric_iv_scaled(candidate.m23,power),
                     _metric_iv_scaled(input.m23,power))
    a[1]>=0 && b[1]>=0 && c[1]>=0 || return false
    _metric_iv_minor(a,b,d)[1]>=0 && _metric_iv_minor(a,c,e)[1]>=0 &&
        _metric_iv_minor(b,c,f)[1]>=0 || return false
    determinant=_metric_iv_add(_metric_iv_mul(_metric_iv_mul(a,b),c),
        _metric_iv_mul((2.0,2.0),_metric_iv_mul(_metric_iv_mul(d,e),f)))
    determinant=_metric_iv_sub(determinant,_metric_iv_mul(a,_metric_iv_square(f)))
    determinant=_metric_iv_sub(determinant,_metric_iv_mul(b,_metric_iv_square(e)))
    determinant=_metric_iv_sub(determinant,_metric_iv_mul(c,_metric_iv_square(d)))
    return determinant[1]>=0
end

function _metric_exact_positive_definite(a::Float64,b::Float64,c::Float64,
                                         d::Float64,e::Float64,f::Float64)
    ra=Rational{BigInt}(a); rb=Rational{BigInt}(b); rc=Rational{BigInt}(c)
    rd=Rational{BigInt}(d); re=Rational{BigInt}(e); rf=Rational{BigInt}(f)
    ra>0 || return false
    ra*rb-rd*rd>0 || return false
    return ra*rb*rc+2rd*re*rf-ra*rf*rf-rb*re*re-rc*rd*rd>0
end

# Sound strict Sylvester certificate. The interval path is allocation-free;
# exact dyadic arithmetic settles only the rounding-indeterminate boundary.
@inline function _metric_positive_definite(a::Float64,b::Float64,c::Float64,
                                           d::Float64,e::Float64,f::Float64)::Bool
    all(isfinite,(a,b,c,d,e,f)) || return false
    (a>0 && b>0 && c>0) || return false
    scale=max(abs(a),abs(b),abs(c),abs(d),abs(e),abs(f))
    power=-exponent(scale)
    ia=_metric_iv_scaled(a,power); ib=_metric_iv_scaled(b,power)
    ic=_metric_iv_scaled(c,power); id=_metric_iv_scaled(d,power)
    ie=_metric_iv_scaled(e,power); iff=_metric_iv_scaled(f,power)
    minor2=_metric_iv_minor(ia,ib,id)
    determinant=_metric_iv_add(_metric_iv_mul(_metric_iv_mul(ia,ib),ic),
        _metric_iv_mul((2.0,2.0),_metric_iv_mul(_metric_iv_mul(id,ie),iff)))
    determinant=_metric_iv_sub(determinant,_metric_iv_mul(ia,_metric_iv_square(iff)))
    determinant=_metric_iv_sub(determinant,_metric_iv_mul(ib,_metric_iv_square(ie)))
    determinant=_metric_iv_sub(determinant,_metric_iv_mul(ic,_metric_iv_square(id)))
    minor2[1]>0 && determinant[1]>0 && return true
    (minor2[2]<=0 || determinant[2]<=0) && return false
    return _metric_exact_positive_definite(a,b,c,d,e,f)
end

function _metric_cholesky_exact(a::Float64,b::Float64,c::Float64,
                                d::Float64,e::Float64,
                                f::Float64)::Union{Nothing,NTuple{6,Float64}}
    ra=Rational{BigInt}(a); rb=Rational{BigInt}(b); rc=Rational{BigInt}(c)
    rd=Rational{BigInt}(d); re=Rational{BigInt}(e); rf=Rational{BigInt}(f)
    l21=rd/ra; l31=re/ra
    d2=rb-l21*rd
    l32=(rf-l31*rd)/d2
    d3=rc-l31*re-l32*l32*d2
    (ra>0 && d2>0 && d3>0) || return nothing
    return setprecision(BigFloat,256) do
        s1=sqrt(BigFloat(ra)); s2=sqrt(BigFloat(d2)); s3=sqrt(BigFloat(d3))
        positive_float(value)=begin
            result=Float64(value)
            result==0 && value>0 ? nextfloat(0.0) : result
        end
        result=(positive_float(s1),Float64(BigFloat(l21)*s1),positive_float(s2),
                Float64(BigFloat(l31)*s1),Float64(BigFloat(l32)*s2),positive_float(s3))
        return all(isfinite,result) ? result : nothing
    end
end

# Exact IEEE-dyadic Loewner certificate for a symmetric 3x3 residual. All seven
# principal minors are nonnegative iff the residual is positive semidefinite.
function _metric_exact_dominates(candidate::Metric3,input::Metric3)
    a=_metric_rational(candidate.m11)-_metric_rational(input.m11)
    b=_metric_rational(candidate.m22)-_metric_rational(input.m22)
    c=_metric_rational(candidate.m33)-_metric_rational(input.m33)
    d=_metric_rational(candidate.m12)-_metric_rational(input.m12)
    e=_metric_rational(candidate.m13)-_metric_rational(input.m13)
    f=_metric_rational(candidate.m23)-_metric_rational(input.m23)
    a>=0 && b>=0 && c>=0 || return false
    a*b-d*d>=0 && a*c-e*e>=0 && b*c-f*f>=0 || return false
    return a*b*c+2*d*e*f-a*f*f-b*e*e-c*d*d>=0
end

@inline function _metric_midpoint(a::Float64,b::Float64)
    return signbit(a)==signbit(b) ? a+(b-a)/2 : a/2+b/2
end

function _metric_round_up(r::Rational{BigInt})
    r>0 || return nothing
    r<=_metric_rational(floatmax(Float64)) || return nothing
    value=Float64(r)
    _metric_rational(value)<r && (value=nextfloat(value))
    return isfinite(value) ? value : nothing
end

function _metric_centered_sdd(a::Metric3,b::Metric3)
    d=_metric_midpoint(a.m12,b.m12)
    e=_metric_midpoint(a.m13,b.m13)
    f=_metric_midpoint(a.m23,b.m23)
    rd=_metric_rational(d); re=_metric_rational(e); rf=_metric_rational(f)
    function diagonal(input::Metric3,which::Int)
        if which==1
            return _metric_rational(input.m11)+abs(rd-_metric_rational(input.m12))+
                   abs(re-_metric_rational(input.m13))
        elseif which==2
            return _metric_rational(input.m22)+abs(rd-_metric_rational(input.m12))+
                   abs(rf-_metric_rational(input.m23))
        end
        return _metric_rational(input.m33)+abs(re-_metric_rational(input.m13))+
               abs(rf-_metric_rational(input.m23))
    end
    r1=max(diagonal(a,1),diagonal(b,1))
    r2=max(diagonal(a,2),diagonal(b,2))
    r3=max(diagonal(a,3),diagonal(b,3))
    c1=_metric_round_up(r1); c1===nothing && return nothing
    c2=_metric_round_up(r2); c2===nothing && return nothing
    c3=_metric_round_up(r3); c3===nothing && return nothing
    candidate=try
        Metric3(c1::Float64,c2::Float64,c3::Float64,d,e,f)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return nothing
    end
    return (_metric_exact_dominates(candidate,a) &&
            _metric_exact_dominates(candidate,b)) ? candidate : nothing
end

const _METRIC_BOUNDARY_PAIRS=((1,2),(1,3),(2,3))

@inline _metric_diagonal(m::Metric3,i::Int)=
    i==1 ? m.m11 : (i==2 ? m.m22 : m.m33)
@inline _metric_offdiagonal(m::Metric3,k::Int)=
    k==1 ? m.m12 : (k==2 ? m.m13 : m.m23)
@inline _metric_boundary_scaled(x::Float64)=ldexp(x,-1023)

function _metric_boundary_interval_disjoint(a::Metric3,b::Metric3,k::Int)
    i,j=_METRIC_BOUNDARY_PAIRS[k]
    F=_metric_rational(floatmax(Float64))
    ra=(F-_metric_rational(_metric_diagonal(a,i)))*
       (F-_metric_rational(_metric_diagonal(a,j)))
    rb=(F-_metric_rational(_metric_diagonal(b,i)))*
       (F-_metric_rational(_metric_diagonal(b,j)))
    delta=_metric_rational(_metric_offdiagonal(a,k))-
          _metric_rational(_metric_offdiagonal(b,k))
    distance2=delta*delta
    excess=distance2-ra-rb
    return distance2>ra+rb && excess*excess>4ra*rb
end

function _metric_boundary_box(a::Metric3,b::Metric3)
    top=_metric_boundary_scaled(floatmax(Float64))
    function bounds(m::Metric3,k::Int)
        i,j=_METRIC_BOUNDARY_PAIRS[k]
        center=_metric_boundary_scaled(_metric_offdiagonal(m,k))
        ri=max(0.0,top-_metric_boundary_scaled(_metric_diagonal(m,i)))
        rj=max(0.0,top-_metric_boundary_scaled(_metric_diagonal(m,j)))
        radius=sqrt(ri*rj)
        return prevfloat(center-radius),nextfloat(center+radius)
    end
    ab=ntuple(k->bounds(a,k),3); bb=ntuple(k->bounds(b,k),3)
    lo=ntuple(k->max(ab[k][1],bb[k][1]),3)
    hi=ntuple(k->min(ab[k][2],bb[k][2]),3)
    return lo,hi
end

@inline function _metric_boundary_forces(m::Metric3,k::Int)
    i,j=_METRIC_BOUNDARY_PAIRS[k]
    F=floatmax(Float64)
    return _metric_diagonal(m,i)==F || _metric_diagonal(m,j)==F
end

function _metric_boundary_inverse_logdet(m::Metric3,y::NTuple{3,Float64},shift::Float64)
    top=_metric_boundary_scaled(floatmax(Float64))
    a=top-_metric_boundary_scaled(m.m11)+shift
    b=top-_metric_boundary_scaled(m.m22)+shift
    c=top-_metric_boundary_scaled(m.m33)+shift
    d=y[1]-_metric_boundary_scaled(m.m12)
    e=y[2]-_metric_boundary_scaled(m.m13)
    f=y[3]-_metric_boundary_scaled(m.m23)
    chol=_metric_cholesky_components(a,b,c,d,e,f)
    chol===nothing && return nothing
    l11,_,l22,_,_,l33=chol::NTuple{6,Float64}
    determinant=a*b*c+2*d*e*f-a*f*f-b*e*e-c*d*d
    (determinant>0 && isfinite(determinant)) || return nothing
    inverse=((b*c-f*f)/determinant,(a*c-e*e)/determinant,
             (a*b-d*d)/determinant,(e*f-c*d)/determinant,
             (d*f-b*e)/determinant,(d*e-a*f)/determinant)
    all(isfinite,inverse) || return nothing
    logdet=2*(log(l11)+log(l22)+log(l33))
    return isfinite(logdet) ? (logdet,inverse) : nothing
end

function _metric_boundary_barrier(a::Metric3,b::Metric3,
                                  y::NTuple{3,Float64},shift::Float64)
    va=_metric_boundary_inverse_logdet(a,y,shift); va===nothing && return nothing
    vb=_metric_boundary_inverse_logdet(b,y,shift); vb===nothing && return nothing
    loga,ga=va; logb,gb=vb
    a11,a22,a33,a12,a13,a23=ga
    b11,b22,b33,b12,b13,b23=gb
    gradient=(2*(a12+b12),2*(a13+b13),2*(a23+b23))
    K=(2*(a11*a22+a12*a12+b11*b22+b12*b12),
       2*(a11*a33+a13*a13+b11*b33+b13*b13),
       2*(a22*a33+a23*a23+b22*b33+b23*b23),
       2*(a11*a23+a13*a12+b11*b23+b13*b12),
       2*(a12*a23+a13*a22+b12*b23+b13*b22),
       2*(a12*a33+a13*a23+b12*b33+b13*b23))
    all(isfinite,gradient) && all(isfinite,K) || return nothing
    return (loga+logb,gradient,K)
end

function _metric_boundary_newton(K::NTuple{6,Float64},gradient::NTuple{3,Float64},
                                 fixed::NTuple{3,Bool})
    a,b,c,d,e,f=K; g1,g2,g3=gradient
    fixed[1] && ((a,d,e,g1)=(1.0,0.0,0.0,0.0))
    fixed[2] && ((b,d,f,g2)=(1.0,0.0,0.0,0.0))
    fixed[3] && ((c,e,f,g3)=(1.0,0.0,0.0,0.0))
    chol=_metric_cholesky_components(a,b,c,d,e,f)
    chol===nothing && return nothing
    l11,l21,l22,l31,l32,l33=chol::NTuple{6,Float64}
    z1=g1/l11; z2=(g2-l21*z1)/l22; z3=(g3-l31*z1-l32*z2)/l33
    x3=z3/l33; x2=(z2-l32*x3)/l22; x1=(z1-l21*x2-l31*x3)/l11
    step=(x1,x2,x3)
    return all(isfinite,step) ? step : nothing
end

@inline function _metric_boundary_clamp(y,lo,hi,fixed,forced)
    return ntuple(k->fixed[k] ? forced[k] : clamp(y[k],lo[k],hi[k]),3)
end

function _metric_boundary_optimize(a,b,initial,shift,lo,hi,fixed,forced)
    y=_metric_boundary_clamp(initial,lo,hi,fixed,forced)
    current=_metric_boundary_barrier(a,b,y,shift)
    current===nothing && return nothing
    for _ in 1:30
        objective,gradient,K=current
        step=_metric_boundary_newton(K,gradient,fixed); step===nothing && return y
        ascent=gradient[1]*step[1]+gradient[2]*step[2]+gradient[3]*step[3]
        (isfinite(ascent) && ascent>16eps(Float64)) || return y
        alpha=1.0; accepted=false
        for _ in 1:40
            trial=ntuple(k->y[k]+alpha*step[k],3)
            trial=_metric_boundary_clamp(trial,lo,hi,fixed,forced)
            displacement=ntuple(k->trial[k]-y[k],3)
            directional=gradient[1]*displacement[1]+gradient[2]*displacement[2]+
                        gradient[3]*displacement[3]
            if directional>0
                value=_metric_boundary_barrier(a,b,trial,shift)
                if value!==nothing && value[1]>=objective+1e-4*directional
                    y=trial; current=value; accepted=true; break
                end
            end
            alpha*=0.5
        end
        accepted || return y
    end
    return y
end

function _metric_boundary_candidate(y,a::Metric3,b::Metric3)
    values=(ldexp(y[1],1023),ldexp(y[2],1023),ldexp(y[3],1023))
    all(isfinite,values) || return nothing
    candidate=try
        Metric3(floatmax(Float64),floatmax(Float64),floatmax(Float64),values...)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return nothing
    end
    return _metric_exact_dominates(candidate,a) && _metric_exact_dominates(candidate,b) ?
           candidate : nothing
end

function _metric_fixed_diagonal_search(a::Metric3,b::Metric3)
    for k in 1:3
        _metric_boundary_interval_disjoint(a,b,k) && return nothing,:impossible
    end
    lo,hi=_metric_boundary_box(a,b)
    any(k->lo[k]>hi[k],1:3) && return nothing,:search
    fixed=ntuple(k->_metric_boundary_forces(a,k)||_metric_boundary_forces(b,k),3)
    forced=ntuple(3) do k
        source=_metric_boundary_forces(a,k) ? a : b
        _metric_boundary_scaled(_metric_offdiagonal(source,k))
    end
    lo=ntuple(k->fixed[k] ? forced[k] : lo[k],3)
    hi=ntuple(k->fixed[k] ? forced[k] : hi[k],3)
    offa=ntuple(k->_metric_boundary_scaled(_metric_offdiagonal(a,k)),3)
    offb=ntuple(k->_metric_boundary_scaled(_metric_offdiagonal(b,k)),3)
    starts=(ntuple(k->_metric_midpoint(lo[k],hi[k]),3),
            ntuple(k->clamp(_metric_midpoint(offa[k],offb[k]),lo[k],hi[k]),3),
            ntuple(k->clamp(offa[k],lo[k],hi[k]),3),
            ntuple(k->clamp(offb[k],lo[k],hi[k]),3))
    for start in starts
        y=_metric_boundary_clamp(start,lo,hi,fixed,forced)
        candidate=_metric_boundary_candidate(y,a,b)
        candidate===nothing || return candidate,:ok
        completed=true
        for k in 0:58
            next=_metric_boundary_optimize(a,b,y,ldexp(8.0,-k),lo,hi,fixed,forced)
            if next===nothing; completed=false; break; end
            y=next
        end
        completed || continue
        next=_metric_boundary_optimize(a,b,y,0.0,lo,hi,fixed,forced)
        next===nothing || (y=next)
        candidate=_metric_boundary_candidate(y,a,b)
        candidate===nothing || return candidate,:ok
    end
    return nothing,:search
end

function _metric_boundary_upper(candidate::Union{Nothing,Metric3},a::Metric3,b::Metric3)
    if candidate!==nothing
        current=candidate::Metric3
        (_metric_exact_dominates(current,a) && _metric_exact_dominates(current,b)) &&
            return current
        for _ in 1:16
            d1=nextfloat(current.m11); d2=nextfloat(current.m22)
            d3=nextfloat(current.m33)
            (isfinite(d1) || isfinite(d2) || isfinite(d3)) || break
            d1=isfinite(d1) ? d1 : current.m11
            d2=isfinite(d2) ? d2 : current.m22
            d3=isfinite(d3) ? d3 : current.m33
            current=Metric3(d1,d2,d3,current.m12,current.m13,current.m23)
            (_metric_exact_dominates(current,a) &&
             _metric_exact_dominates(current,b)) && return current
        end
    end
    sdd=_metric_centered_sdd(a,b)
    sdd===nothing || return sdd
    maximal=Metric3(floatmax(Float64),floatmax(Float64),floatmax(Float64),0.0,0.0,0.0)
    (_metric_exact_dominates(maximal,a) && _metric_exact_dominates(maximal,b)) &&
        return maximal
    searched,status=_metric_fixed_diagonal_search(a,b)
    searched===nothing || return searched
    status===:impossible && throw(ErrorException(
        "intersection_alauzet: no finite Float64 metric can dominate both inputs (exact range-boundary interval proof)"))
    throw(ErrorException(
        "intersection_alauzet: unable to certify a common upper metric at the Float64 range boundary"))
end

@inline function _inflate_to_dominate(candidate::Metric3,a::Metric3,b::Metric3)
    guard=64sqrt(eps(Float64))
    for _ in 1:6
        (_metric_interval_dominates(candidate,a) &&
         _metric_interval_dominates(candidate,b)) && return candidate
        ratio=min(_metric_generalized_min(candidate,a),
                  _metric_generalized_min(candidate,b))
        (isfinite(ratio) && ratio>0) || return _metric_boundary_upper(candidate,a,b)
        # A ratio that claims dominance while the one-sided interval proof does
        # not is a numerical-uncertainty signal, never an acceptance proof.
        ratio>=1+guard/2 && return _metric_boundary_upper(candidate,a,b)
        # Never shrink an uncertified candidate, even if the generalized ratio
        # is slightly greater than one through rounding.
        factor=max(1+guard,(1+guard)/ratio)
        isfinite(factor) || return _metric_boundary_upper(candidate,a,b)
        limit=floatmax(Float64)/factor
        maxentry=max(abs(candidate.m11),abs(candidate.m22),abs(candidate.m33),
                     abs(candidate.m12),abs(candidate.m13),abs(candidate.m23))
        maxentry<=limit || return _metric_boundary_upper(candidate,a,b)
        nextcandidate=try
            _metric_scale(candidate,factor)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            return _metric_boundary_upper(candidate,a,b)
        end
        candidate=nextcandidate
    end
    return _metric_boundary_upper(candidate,a,b)
end

"""
Anisotropic size field. Unless a Gmsh field type defines a distinct scalar
operator, `field_value` returns `1/√λ_max(M)`.
"""

"""
    metric_at(field, x, y, z) -> Metric3

Evaluate the anisotropic metric of `field`; isotropic size fields are converted
to `I / h²`.
"""
function metric_at end

function metric_at(field::AbstractSizeField,x,y,z)
    field isa AbstractAnisoField &&
        throw(ErrorException("metric_at: $(typeof(field)) must implement metric_at"))
    return isotropic_metric(size_at(field,x,y,z))
end

function field_value(field::AbstractAnisoField,x,y,z)
    return metric_size(metric_at(field,x,y,z))
end

"""
    directional_size(field, point, direction) -> Float64

Return the physical length corresponding to unit metric length in `direction`.
For an isotropic field this is exactly `size_at(field, point)`.
"""
function directional_size(field::AbstractSizeField,p,direction;entity=nothing)
    q=_point3(p,"directional_size")
    d=_point3(direction,"directional_size")
    scale=max(abs(d[1]),abs(d[2]),abs(d[3]))
    scale>0 || throw(ArgumentError("directional_size: direction must have positive length"))
    field isa AbstractAnisoField || return size_at(field,q...,entity)
    scaled=(d[1]/scale,d[2]/scale,d[3]/scale)
    n=_norm3(scaled)
    u=(scaled[1]/n,scaled[2]/n,scaled[3]/n)
    m=metric_at(field,q...,entity)
    m isa Metric3 || throw(ArgumentError(
        "directional_size: metric_at returned $(typeof(m)), not Metric3"))
    value=_metric_displacement(m,u...,"directional_size")
    value>0 || throw(ArgumentError(
        "directional_size: metric is not positive in direction $u at $q"))
    h=inv(value)
    (isfinite(h) && h>0) || throw(ArgumentError(
        "directional_size: directional size is not representable at $q"))
    return h
end

@inline function _metric_linear3(a::Float64,x::Float64,b::Float64,y::Float64,
                                 c::Float64,z::Float64)
    p1=a*x; p2=b*y; p3=c*z
    all(isfinite,(p1,p2,p3)) || return (0.0,false)
    # Subnormal products and rounded-to-zero products have no useful relative
    # error bound. They are rare and go through the exact dyadic path below.
    for (coefficient,component,product) in ((a,x,p1),(b,y,p2),(c,z,p3))
        coefficient!=0 && component!=0 &&
            (product==0 || abs(product)<floatmin(Float64)) && return (0.0,false)
    end
    magnitude=abs(p1)+abs(p2)+abs(p3)
    isfinite(magnitude) || return (0.0,false)
    value=p1+p2+p3
    isfinite(value) || return (0.0,false)
    magnitude==0 && return (value,true)
    abs(value)>64eps(Float64)*magnitude || return (value,false)
    return (value,true)
end

@inline function _metric_displacement_fast(m::Metric3,dx::Float64,dy::Float64,
                                            dz::Float64)
    l11,l21,l22,l31,l32,l33=_metric_cholesky_components(
        m.m11,m.m22,m.m33,m.m12,m.m13,m.m23)::NTuple{6,Float64}
    u,safeu=_metric_linear3(l11,dx,l21,dy,l31,dz)
    v,safev=_metric_linear3(0.0,dx,l22,dy,l32,dz)
    w,safew=_metric_linear3(0.0,dx,0.0,dy,l33,dz)
    safeu && safev && safew || return nothing
    result=hypot(u,v,w)
    return isfinite(result) && (result>0 || (dx==0 && dy==0 && dz==0)) ? result : nothing
end

function _metric_positive_sqrt(value::Rational{BigInt},caller::AbstractString)
    value>0 || throw(ErrorException(
        "$caller: a positive-definite metric produced a non-positive squared length"))
    result=setprecision(BigFloat,256) do
        setrounding(BigFloat,RoundUp) do
            Float64(sqrt(BigFloat(value)),RoundUp)
        end
    end
    isfinite(result) || throw(ArgumentError(
        "$caller: metric displacement is outside the representable Float64 range"))
    return result
end

function _metric_displacement_exact(m::Metric3,dx,dy,dz,caller::AbstractString)
    x=_metric_rational(dx); y=_metric_rational(dy); z=_metric_rational(dz)
    value=_metric_rational(m.m11)*x*x+
          _metric_rational(m.m22)*y*y+
          _metric_rational(m.m33)*z*z+
          2*(_metric_rational(m.m12)*x*y+
             _metric_rational(m.m13)*x*z+
             _metric_rational(m.m23)*y*z)
    return _metric_positive_sqrt(value,caller)
end

function _metric_endpoint_displacement_exact(m::Metric3,p,q,caller::AbstractString)
    dx=_metric_rational(q[1])-_metric_rational(p[1])
    dy=_metric_rational(q[2])-_metric_rational(p[2])
    dz=_metric_rational(q[3])-_metric_rational(p[3])
    return _metric_displacement_exact(m,dx,dy,dz,caller)
end

@inline function _metric_displacement(m::Metric3,dx::Float64,dy::Float64,dz::Float64,
                                      caller::AbstractString)
    # For M=L*Lᵀ, evaluate ‖Lᵀd‖ with `hypot`. Expanding dᵀMd can
    # become spuriously negative through cancellation for rotated metrics and
    # can overflow while the final norm is still representable.
    result=_metric_displacement_fast(m,dx,dy,dz)
    result===nothing || return result::Float64
    return _metric_displacement_exact(m,dx,dy,dz,caller)
end

@inline function _metric_endpoint_displacement(m::Metric3,p,q,caller::AbstractString)
    dx=q[1]-p[1]; dy=q[2]-p[2]; dz=q[3]-p[3]
    if all(isfinite,(dx,dy,dz))
        result=_metric_displacement_fast(m,dx,dy,dz)
        result===nothing || return result::Float64
    end
    return _metric_endpoint_displacement_exact(m,p,q,caller)
end

function _isotropic_endpoint_ratio_exact(p,q,h::Float64,caller::AbstractString)
    dx=_metric_rational(q[1])-_metric_rational(p[1])
    dy=_metric_rational(q[2])-_metric_rational(p[2])
    dz=_metric_rational(q[3])-_metric_rational(p[3])
    rh=_metric_rational(h)
    value=(dx*dx+dy*dy+dz*dz)/(rh*rh)
    return _metric_positive_sqrt(value,caller)
end

@inline function _isotropic_endpoint_ratio(p,q,dx,dy,dz,h::Float64,
                                           caller::AbstractString)
    distance=hypot(dx,dy,dz)
    value=distance/h
    isfinite(value) && value>0 && return value
    return _isotropic_endpoint_ratio_exact(p,q,h,caller)
end

@inline function _coordinate_midpoint(a::Float64,b::Float64)
    if signbit(a)!=signbit(b) ||
       (abs(a)<=floatmax(Float64)/2 && abs(b)<=floatmax(Float64)/2)
        return (a+b)/2
    end
    return a/2+b/2
end

"""
    metric_edge_length(field, a, b) -> Float64

Maximum metric length of edge `a`–`b`, sampling the metric at both endpoints
and the midpoint. A value at most one satisfies the field's directional edge
size contract.
"""
function metric_edge_length(field::AbstractSizeField,a,b;entity=nothing)
    p=_point3(a,"metric_edge_length")
    q=_point3(b,"metric_edge_length")
    dx=q[1]-p[1]; dy=q[2]-p[2]; dz=q[3]-p[3]
    p==q && return 0.0
    mx=_coordinate_midpoint(p[1],q[1]); my=_coordinate_midpoint(p[2],q[2])
    mz=_coordinate_midpoint(p[3],q[3])
    if !(field isa AbstractAnisoField)
        h0=size_at(field,p...,entity); h1=size_at(field,q...,entity)
        hm=size_at(field,mx,my,mz,entity)
        return max(_isotropic_endpoint_ratio(p,q,dx,dy,dz,h0,"metric_edge_length"),
                   _isotropic_endpoint_ratio(p,q,dx,dy,dz,h1,"metric_edge_length"),
                   _isotropic_endpoint_ratio(p,q,dx,dy,dz,hm,"metric_edge_length"))
    end
    m0=metric_at(field,p...,entity); m1=metric_at(field,q...,entity)
    mm=metric_at(field,mx,my,mz,entity)
    (m0 isa Metric3 && m1 isa Metric3 && mm isa Metric3) || throw(ArgumentError(
        "metric_edge_length: metric_at must return Metric3"))
    return max(_metric_endpoint_displacement(m0,p,q,"metric_edge_length"),
               _metric_endpoint_displacement(m1,p,q,"metric_edge_length"),
               _metric_endpoint_displacement(mm,p,q,"metric_edge_length"))
end

# Trapezoidal metric increment used by curve integration. Keeping this separate
# from the max-sampled edge certificate preserves Gmsh-style arc-length grading.
function _metric_curve_increment(field::AbstractSizeField,a,b,entity=nothing)
    p=_point3(a,"metric curve increment"); q=_point3(b,"metric curve increment")
    dx=q[1]-p[1]; dy=q[2]-p[2]; dz=q[3]-p[3]
    p==q && return 0.0
    if !(field isa AbstractAnisoField)
        first=_isotropic_endpoint_ratio(p,q,dx,dy,dz,size_at(field,p...,entity),
                                        "metric curve increment")
        second=_isotropic_endpoint_ratio(p,q,dx,dy,dz,size_at(field,q...,entity),
                                         "metric curve increment")
        return first/2+second/2
    end
    m0=metric_at(field,p...,entity); m1=metric_at(field,q...,entity)
    (m0 isa Metric3 && m1 isa Metric3) || throw(ArgumentError(
        "metric curve increment: metric_at must return Metric3"))
    first=_metric_endpoint_displacement(m0,p,q,"metric curve increment")
    second=_metric_endpoint_displacement(m1,p,q,"metric curve increment")
    return first/2+second/2
end

# Entity-aware evaluation. 4-arg methods remain the default; the 5-arg form is
# used by Restrict/Constant. Existing fields ignore the entity.
field_value(field::AbstractField,x,y,z,::Nothing)=field_value(field,x,y,z)
field_value(field::AbstractField,x,y,z,::Tuple{T,U}) where {T<:Integer,U<:Integer}=
    field_value(field,x,y,z)
metric_at(field::AbstractAnisoField,x,y,z,::Nothing)=metric_at(field,x,y,z)
metric_at(field::AbstractAnisoField,x,y,z,::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=metric_at(field,x,y,z)
function metric_at(field::AbstractSizeField,x,y,z,entity)
    field isa AbstractAnisoField && throw(ErrorException(
        "metric_at: $(typeof(field)) must implement entity-aware metric_at"))
    return isotropic_metric(size_at(field,x,y,z,entity))
end
field_value(field::AbstractAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=metric_size(metric_at(field,x,y,z,entity))
function size_at(field::AbstractSizeField,x,y,z,::Nothing)
    return _checked_field_result(field_value(field,x,y,z,nothing),"size_at",x,y,z;positive=true)
end
function size_at(field::AbstractSizeField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}
    dim=Int(entity[1]); tag=Int(entity[2])
    dim in 0:3 || throw(ArgumentError("size_at: entity dimension must be in 0:3"))
    tag>0 || throw(ArgumentError("size_at: entity tag must be positive"))
    return _checked_field_result(field_value(field,x,y,z,entity),"size_at",x,y,z;positive=true)
end
size_at(field::AbstractSizeField,x,y,z,entity)=throw(ArgumentError(
    "size_at: entity context must be nothing or a (dimension, tag) integer tuple (got $entity)"))

# ── MathEval ──────────────────────────────────────────────────────────────────

const _MATH_FUN1=Dict{String,UInt8}(
    "abs"=>1,"fabs"=>1,"acos"=>2,"asin"=>3,"atan"=>4,"cos"=>5,"cosh"=>6,
    "deg"=>7,"exp"=>8,"fac"=>9,"log"=>10,"log10"=>11,"rad"=>12,
    "round"=>13,"sign"=>14,"sin"=>15,"sinh"=>16,"sqrt"=>17,"step"=>18,
    "tan"=>19,"tanh"=>20,"atanh"=>21,"trunc"=>22,"floor"=>23,"ceil"=>24)
const _MATH_FUNN=Dict{String,UInt8}("sum"=>1,"max"=>2,"min"=>3,"med"=>4)

abstract type _MENode end
struct _MELit <: _MENode; v::Float64; end
struct _MEVar <: _MENode; which::UInt8; end          # 0=x,1=y,2=z,3=pi
struct _MEField <: _MENode; tag::Int; end
struct _MEUn <: _MENode; op::UInt8; a::_MENode; end  # 0=-, 1=+
struct _MEBin <: _MENode; op::UInt8; a::_MENode; b::_MENode; end # 0+ 1- 2* 3/ 4^ 5%
struct _MECall1 <: _MENode; f::UInt8; a::_MENode; end
struct _MECallN <: _MENode; f::UInt8; args::Vector{_MENode}; end

struct _MEInstr
    op::UInt8
    a::Int32
    b::Int32
    value::Float64
end
struct _MEProgram
    code::Vector{_MEInstr}
    root::Int32
    field_tags::Vector{Int}
end

mutable struct _MELex
    s::String
    i::Int
    tok::Symbol
    num::Float64
    ident::String
    nodes::Int
    depth::Int
end

@inline _me_ascii_digit(c::Char)=('0'<=c<='9')
@inline _me_ascii_alpha(c::Char)=(('a'<=c<='z')||('A'<=c<='Z'))

@inline function _me_node!(L::_MELex,node::_MENode)
    L.nodes+=1
    L.nodes<=1024 || throw(ArgumentError("MathEval: expression exceeds 1024 syntax nodes"))
    return node
end

function _me_lex!(L::_MELex)
    s=L.s; n=lastindex(s)
    while L.i<=n && isspace(s[L.i]); L.i=nextind(s,L.i); end
    if L.i>n; L.tok=:eof; return; end
    c=s[L.i]
    if c in ('+','-','*','/','%','^','(',')',',')
        L.tok=Symbol(c); L.i=nextind(s,L.i); return
    end
    if c=='.' || _me_ascii_digit(c)
        j=L.i
        while j<=n && (_me_ascii_digit(s[j]) || s[j]=='.' ||
                       (s[j] in ('e','E','+','-') && j>L.i))
            # allow exponent sign only after e/E
            if s[j] in ('+','-')
                previous=prevind(s,j)
                s[previous] in ('e','E') || break
            end
            j=nextind(s,j)
        end
        raw=s[L.i:prevind(s,j)]; v=tryparse(Float64,raw)
        v===nothing && throw(ArgumentError("MathEval: invalid number $raw"))
        L.num=v; L.tok=:num; L.i=j; return
    end
    if _me_ascii_alpha(c) || c=='_'
        j=L.i
        while j<=n && (_me_ascii_alpha(s[j]) || _me_ascii_digit(s[j]) || s[j]=='_')
            j=nextind(s,j)
        end
        L.ident=s[L.i:prevind(s,j)]; L.tok=:ident; L.i=j; return
    end
    throw(ArgumentError("MathEval: unexpected character $(repr(c)) at position $(L.i)"))
end

function parse_matheval(expr::AbstractString)
    raw=String(strip(expr))
    isempty(raw) && throw(ArgumentError("MathEval: empty expression"))
    ncodeunits(raw)<=65536 || throw(ArgumentError(
        "MathEval: expression exceeds the 65536-byte resource limit"))
    L=_MELex(raw,1,:eof,0.0,"",0,0)
    _me_lex!(L)
    node=_me_expr!(L)
    L.tok===:eof || throw(ArgumentError("MathEval: trailing input after expression"))
    return _me_compile(node)
end

function _me_accept!(L,t); L.tok===t || return false; _me_lex!(L); true; end
function _me_expect!(L,t)
    L.tok===t || throw(ArgumentError("MathEval: expected $t, got $(L.tok)"))
    _me_lex!(L)
end

function _me_expr!(L)
    a=_me_term!(L)
    while L.tok===:+ || L.tok===:-
        op=L.tok===:+ ? UInt8(0) : UInt8(1); _me_lex!(L)
        (L.tok===:+ || L.tok===:-) && throw(ArgumentError(
            "MathEval: signed operands after binary operators must be parenthesized"))
        a=_me_node!(L,_MEBin(op,a,_me_term!(L)))
    end
    return a
end
function _me_term!(L)
    a=_me_power!(L)
    while L.tok===:* || L.tok===:/ || L.tok===Symbol("%")
        op=L.tok===:* ? UInt8(2) : L.tok===:/ ? UInt8(3) : UInt8(5); _me_lex!(L)
        (L.tok===:+ || L.tok===:-) && throw(ArgumentError(
            "MathEval: signed operands after binary operators must be parenthesized"))
        a=_me_node!(L,_MEBin(op,a,_me_power!(L)))
    end
    return a
end
function _me_power!(L)
    a=_me_unary!(L)
    if L.tok===:^
        _me_lex!(L)
        (L.tok===:+ || L.tok===:-) && throw(ArgumentError(
            "MathEval: parenthesize a signed exponent"))
        # Gmsh MathEx accepts one power at this level; `2^3^2` is rejected.
        a=_me_node!(L,_MEBin(UInt8(4),a,_me_unary!(L)))
    end
    return a
end
function _me_unary!(L)
    if L.tok===:-
        _me_lex!(L)
        (L.tok===:+ || L.tok===:-) && throw(ArgumentError(
            "MathEval: consecutive unary signs are invalid"))
        return _me_node!(L,_MEUn(UInt8(0),_me_primary!(L)))
    elseif L.tok===:+
        _me_lex!(L)
        (L.tok===:+ || L.tok===:-) && throw(ArgumentError(
            "MathEval: consecutive unary signs are invalid"))
        return _me_node!(L,_MEUn(UInt8(1),_me_primary!(L)))
    end
    return _me_primary!(L)
end
const _ME_LPAREN=Symbol("(")
const _ME_RPAREN=Symbol(")")
const _ME_COMMA=Symbol(",")

function _me_primary!(L)
    if L.tok===:num
        v=L.num; _me_lex!(L); return _me_node!(L,_MELit(v))
    elseif L.tok===:ident
        name=L.ident; _me_lex!(L)
        if L.tok===_ME_LPAREN
            _me_lex!(L)
            L.depth+=1
            L.depth<=256 || throw(ArgumentError(
                "MathEval: expression nesting exceeds the 256-level resource limit"))
            args=_MENode[]
            try
                if L.tok!==_ME_RPAREN
                    push!(args,_me_expr!(L))
                    while L.tok===_ME_COMMA
                        _me_lex!(L); push!(args,_me_expr!(L))
                    end
                end
                _me_expect!(L,_ME_RPAREN)
            finally
                L.depth-=1
            end
            key=lowercase(name)
            gmsh_case=(name==key || name==uppercasefirst(key))
            gmsh_case && key=="rand" && throw(ArgumentError(
                "MathEval: rand() requires nondeterministic shared state and is unsupported; use FunctionSize with an explicit RNG"))
            if gmsh_case && haskey(_MATH_FUN1,key)
                length(args)==1 || throw(ArgumentError(
                    "MathEval: function $name requires exactly one argument"))
                return _me_node!(L,_MECall1(_MATH_FUN1[key],args[1]))
            elseif gmsh_case && haskey(_MATH_FUNN,key)
                isempty(args) && throw(ArgumentError(
                    "MathEval: function $name requires at least one argument"))
                return _me_node!(L,_MECallN(_MATH_FUNN[key],args))
            end
            throw(ArgumentError("MathEval: unknown function $name"))
        end
        name=="x" && return _me_node!(L,_MEVar(UInt8(0)))
        name=="y" && return _me_node!(L,_MEVar(UInt8(1)))
        name=="z" && return _me_node!(L,_MEVar(UInt8(2)))
        name in ("Pi","pi") && return _me_node!(L,_MEVar(UInt8(3)))
        name=="e" && return _me_node!(L,_MELit(Float64(MathConstants.e)))
        if length(name)>=2 && name[1]=='F' && all(_me_ascii_digit,name[2:end])
            tag=tryparse(Int,name[2:end])
            tag===nothing && throw(ArgumentError("MathEval: field tag is outside the Int range"))
            tag>=0 || throw(ArgumentError("MathEval: field tag $tag is invalid"))
            return _me_node!(L,_MEField(tag))
        end
        throw(ArgumentError("MathEval: unknown identifier $name"))
    elseif L.tok===_ME_LPAREN
        _me_lex!(L); L.depth+=1
        L.depth<=256 || throw(ArgumentError(
            "MathEval: expression nesting exceeds the 256-level resource limit"))
        try
            a=_me_expr!(L); _me_expect!(L,_ME_RPAREN); return a
        finally
            L.depth-=1
        end
    end
    throw(ArgumentError("MathEval: unexpected token $(L.tok)"))
end

function _me_emit!(code::Vector{_MEInstr},tags::Vector{Int},node::_MENode)::Int32
    if node isa _MELit
        push!(code,_MEInstr(0,0,0,node.v))
    elseif node isa _MEVar
        push!(code,_MEInstr(1,Int32(node.which),0,0.0))
    elseif node isa _MEField
        node.tag<=typemax(Int32) || throw(ArgumentError(
            "MathEval: field tag $(node.tag) exceeds the Int32 evaluator limit"))
        push!(tags,node.tag); push!(code,_MEInstr(2,Int32(node.tag),0,0.0))
    elseif node isa _MEUn
        a=_me_emit!(code,tags,node.a)
        node.op==1 && return a
        push!(code,_MEInstr(3,a,0,0.0))
    elseif node isa _MEBin
        a=_me_emit!(code,tags,node.a); b=_me_emit!(code,tags,node.b)
        push!(code,_MEInstr(UInt8(4+node.op),a,b,0.0))
    elseif node isa _MECall1
        a=_me_emit!(code,tags,node.a)
        push!(code,_MEInstr(UInt8(32+node.f),a,0,0.0))
    else
        args=(node::_MECallN).args
        root=_me_emit!(code,tags,args[1])
        foldop=node.f==2 ? UInt8(10) : node.f==3 ? UInt8(11) : UInt8(4)
        @inbounds for i in 2:length(args)
            b=_me_emit!(code,tags,args[i]); push!(code,_MEInstr(foldop,root,b,0.0))
            root=Int32(length(code))
        end
        if node.f==4
            push!(code,_MEInstr(0,0,0,Float64(length(args)))); n=Int32(length(code))
            push!(code,_MEInstr(7,root,n,0.0)); root=Int32(length(code))
        end
        return root
    end
    return Int32(length(code))
end

function _me_compile(node::_MENode)
    code=_MEInstr[]; tags=Int[]; root=_me_emit!(code,tags,node)
    sort!(unique!(tags))
    return _MEProgram(code,root,tags)
end
matheval_field_tags(program::_MEProgram)=copy(program.field_tags)

@inline function _me_fac(x::Float64)
    (isfinite(x) && 0<=x<=170) || throw(ArgumentError("MathEval: fac argument is outside [0,170]"))
    x<2 && return 1.0
    n=trunc(Int,x+0.5)
    n==2 && return 2.0
    value=1.0
    @inbounds for k in 2:n; value*=k; end
    return value
end
@inline function _me_call1(f::UInt8,x::Float64)
    f==1 && return abs(x); f==2 && return acos(x); f==3 && return asin(x)
    f==4 && return atan(x); f==5 && return cos(x); f==6 && return cosh(x)
    f==7 && return x*180/π; f==8 && return exp(x); f==9 && return _me_fac(x)
    f==10 && return log(x); f==11 && return log10(x); f==12 && return x*π/180
    f==13 && return trunc(x+0.5); f==14 && return x<0 ? -1.0 : 1.0
    f==15 && return sin(x); f==16 && return sinh(x); f==17 && return sqrt(x)
    f==18 && return x<0 ? 0.0 : 1.0; f==19 && return tan(x)
    f==20 && return tanh(x); f==21 && return atanh(x); f==22 && return trunc(x)
    f==23 && return floor(x); return ceil(x)
end

struct _MEFieldBinding{F<:AbstractField}
    tag::Int
    field::F
end

function _me_field_bindings(fields)
    table=try
        Dict{Int,AbstractField}(fields)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("MathEval: fields must map integer tags to AbstractField values"))
    end
    tags=sort!(collect(keys(table)))
    return Tuple(_MEFieldBinding(tag,table[tag]) for tag in tags)
end

@noinline function _me_unknown_field(tag::Int)
    # Gmsh warns and substitutes MAX_LC for an unresolved F<tag> reference.
    return GMSH_MAX_SIZE
end

@generated function _me_field_value(bindings::T,tag::Int,x::Float64,y::Float64,
                                    z::Float64,entity::E)::Float64 where
                                    {T<:Tuple,E}
    body=:(return _me_unknown_field(tag))
    for i in fieldcount(T):-1:1
        body=quote
            binding=getfield(bindings,$i)
            if binding.tag==tag
                return _checked_field_result(field_value(binding.field,x,y,z,entity),
                                             "MathEval field reference",x,y,z)
            end
            $body
        end
    end
    return body
end

function _me_eval(program::_MEProgram,index::Int32,x::Float64,y::Float64,z::Float64,
                  fields::Tuple,entity)::Float64
    ins=@inbounds program.code[index]; op=ins.op
    op==0 && return ins.value
    if op==1
        ins.a==0 && return x; ins.a==1 && return y; ins.a==2 && return z
        return Float64(π)
    elseif op==2
        return _me_field_value(fields,Int(ins.a),x,y,z,entity)
    elseif op==3
        return -_me_eval(program,ins.a,x,y,z,fields,entity)
    elseif op>=33
        return _me_call1(UInt8(op-32),_me_eval(program,ins.a,x,y,z,fields,entity))
    end
    a=_me_eval(program,ins.a,x,y,z,fields,entity)
    b=_me_eval(program,ins.b,x,y,z,fields,entity)
    op==4 && return a+b; op==5 && return a-b; op==6 && return a*b
    op==7 && (b==0 && throw(ArgumentError("MathEval: division by zero")); return a/b)
    op==8 && return a^b
    op==9 && return rem(a,b)
    op==10 && return max(a,b)
    return min(a,b)
end

eval_matheval(program::_MEProgram,x,y,z,fields::Tuple)=
    _me_eval(program,program.root,Float64(x),Float64(y),Float64(z),fields,nothing)
eval_matheval(program::_MEProgram,x,y,z,fields)=
    eval_matheval(program,x,y,z,_me_field_bindings(fields))

@inline _eval_matheval_entity(program::_MEProgram,x,y,z,fields::Tuple,entity)=
    _me_eval(program,program.root,Float64(x),Float64(y),Float64(z),fields,entity)

"""
    MathEvalField(expr; fields=Dict{Int,AbstractField}())

Gmsh `MathEval`. `expr` may use `x,y,z,Pi` and `F<tag>` references into `fields`.
"""
struct MathEvalField{F<:Tuple} <: AbstractSizeField
    expr::String
    ast::_MEProgram
    fields::F
end
function MathEvalField(expr::AbstractString; fields=Dict{Int,AbstractField}())
    ast=parse_matheval(expr)
    bindings=_me_field_bindings(fields)
    return MathEvalField{typeof(bindings)}(String(expr),ast,bindings)
end
function _matheval_value(field::MathEvalField,x,y,z,entity)
    x=_float_value(x,"MathEvalField","x"); y=_float_value(y,"MathEvalField","y")
    z=_float_value(z,"MathEvalField","z")
    v=try
        _eval_matheval_entity(field.ast,x,y,z,field.fields,entity)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("MathEvalField: cannot evaluate $(field.expr) at ($x,$y,$z): $err"))
    end
    return _checked_field_result(v,"MathEvalField",x,y,z)
end
field_value(field::MathEvalField,x,y,z)=_matheval_value(field,x,y,z,nothing)
field_value(field::MathEvalField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_matheval_value(field,x,y,z,entity)

"""Gmsh `MathEvalAniso`: six expressions for the packed SPD metric."""
struct MathEvalAnisoField{F<:Tuple} <: AbstractAnisoField
    m11::_MEProgram; m22::_MEProgram; m33::_MEProgram
    m12::_MEProgram; m13::_MEProgram; m23::_MEProgram
    fields::F
end
function MathEvalAnisoField(; m11,m22,m33,m12="0",m13="0",m23="0",
                            fields=Dict{Int,AbstractField}())
    bindings=_me_field_bindings(fields)
    return MathEvalAnisoField{typeof(bindings)}(
        parse_matheval(m11),parse_matheval(m22),parse_matheval(m33),
        parse_matheval(m12),parse_matheval(m13),parse_matheval(m23),bindings)
end
function _matheval_aniso_metric(field::MathEvalAnisoField,x,y,z,entity)
    x=_float_value(x,"MathEvalAnisoField","x"); y=_float_value(y,"MathEvalAnisoField","y")
    z=_float_value(z,"MathEvalAnisoField","z")
    ev(n)=_eval_matheval_entity(n,x,y,z,field.fields,entity)
    m=Metric3(ev(field.m11),ev(field.m22),ev(field.m33),ev(field.m12),ev(field.m13),ev(field.m23))
    λ=metric_eigenvalues(m)
    minimum(λ)>0 || throw(ArgumentError("MathEvalAnisoField: metric is not positive definite at ($x,$y,$z)"))
    return m
end
metric_at(field::MathEvalAnisoField,x,y,z)=
    _matheval_aniso_metric(field,x,y,z,nothing)
metric_at(field::MathEvalAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_matheval_aniso_metric(field,x,y,z,entity)

# Gmsh's scalar MathEvalAniso operator is deliberately not a mesh length: it
# evaluates M11 only. MathEval `F<tag>` references and scalar differential fields
# call this operator, while anisotropic mesh consumers call `metric_at` above.
function _matheval_aniso_value(field::MathEvalAnisoField,x,y,z,entity)
    x=_float_value(x,"MathEvalAnisoField","x");y=_float_value(y,"MathEvalAnisoField","y")
    z=_float_value(z,"MathEvalAnisoField","z")
    value=_eval_matheval_entity(field.m11,x,y,z,field.fields,entity)
    return _checked_field_result(value,"MathEvalAnisoField",x,y,z)
end
field_value(field::MathEvalAnisoField,x,y,z)=
    _matheval_aniso_value(field,x,y,z,nothing)
field_value(field::MathEvalAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_matheval_aniso_value(field,x,y,z,entity)

# ── Finite-difference operators ───────────────────────────────────────────────

function _fd_delta(delta,caller)
    d=_float_value(delta,caller,"delta")
    d>0 || throw(ArgumentError("$caller: delta must be positive"))
    return d
end

"""Gmsh `Gradient`: component 0/1/2 = ∂/∂x,∂/∂y,∂/∂z; 3 = Euclidean norm."""
struct GradientField{F<:AbstractField} <: AbstractField
    input::F; kind::Int; delta::Float64
end
function GradientField(input::AbstractField; kind::Integer=3, delta::Real=1e-4)
    k=Int(kind); k in 0:3 || throw(ArgumentError("GradientField: kind must be 0,1,2, or 3"))
    return GradientField{typeof(input)}(input,k,_fd_delta(delta,"GradientField"))
end
function _gradient_value(field::GradientField,x,y,z,entity)
    x=_float_value(x,"GradientField","x"); y=_float_value(y,"GradientField","y")
    z=_float_value(z,"GradientField","z"); d=field.delta; g=field.input
    gx()=(field_value(g,x+d/2,y,z,entity)-field_value(g,x-d/2,y,z,entity))/d
    gy()=(field_value(g,x,y+d/2,z,entity)-field_value(g,x,y-d/2,z,entity))/d
    gz()=(field_value(g,x,y,z+d/2,entity)-field_value(g,x,y,z-d/2,entity))/d
    value = field.kind==0 ? gx() : field.kind==1 ? gy() : field.kind==2 ? gz() :
            hypot(gx(),gy(),gz())
    return _checked_field_result(value,"GradientField",x,y,z)
end
field_value(field::GradientField,x,y,z)=_gradient_value(field,x,y,z,nothing)
field_value(field::GradientField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_gradient_value(field,x,y,z,entity)

"""Gmsh `Laplacian` of the input field."""
struct LaplacianField{F<:AbstractField} <: AbstractField
    input::F; delta::Float64
end
LaplacianField(input::AbstractField; delta::Real=1e-4) =
    LaplacianField{typeof(input)}(input,_fd_delta(delta,"LaplacianField"))
function _laplacian_value(field::LaplacianField,x,y,z,entity)
    x=_float_value(x,"LaplacianField","x"); y=_float_value(y,"LaplacianField","y")
    z=_float_value(z,"LaplacianField","z"); d=field.delta; g=field.input
    value=(field_value(g,x+d,y,z,entity)+field_value(g,x-d,y,z,entity)+
            field_value(g,x,y+d,z,entity)+field_value(g,x,y-d,z,entity)+
            field_value(g,x,y,z+d,entity)+field_value(g,x,y,z-d,entity)-
            6*field_value(g,x,y,z,entity))/(d*d)
    return _checked_field_result(value,"LaplacianField",x,y,z)
end
field_value(field::LaplacianField,x,y,z)=_laplacian_value(field,x,y,z,nothing)
field_value(field::LaplacianField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_laplacian_value(field,x,y,z,entity)

"""Gmsh `Mean` of seven samples of the input field."""
struct MeanField{F<:AbstractField} <: AbstractField
    input::F; delta::Float64
end
MeanField(input::AbstractField; delta::Real=1e-4) =
    MeanField{typeof(input)}(input,_fd_delta(delta,"MeanField"))
function _mean_value(field::MeanField,x,y,z,entity)
    x=_float_value(x,"MeanField","x"); y=_float_value(y,"MeanField","y")
    z=_float_value(z,"MeanField","z"); d=field.delta; g=field.input
    value=(field_value(g,x+d,y,z,entity)+field_value(g,x-d,y,z,entity)+
            field_value(g,x,y+d,z,entity)+field_value(g,x,y-d,z,entity)+
            field_value(g,x,y,z+d,entity)+field_value(g,x,y,z-d,entity)+
            field_value(g,x,y,z,entity))/7
    return _checked_field_result(value,"MeanField",x,y,z)
end
field_value(field::MeanField,x,y,z)=_mean_value(field,x,y,z,nothing)
field_value(field::MeanField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_mean_value(field,x,y,z,entity)

"""Gmsh `Curvature`: `div(grad(F)/|grad(F)|)` by finite differences."""
struct CurvatureOpField{F<:AbstractField} <: AbstractField
    input::F; delta::Float64
end
CurvatureOpField(input::AbstractField; delta::Real=1e-4) =
    CurvatureOpField{typeof(input)}(input,_fd_delta(delta,"CurvatureField"))
"""Alias for [`CurvatureOpField`](@ref), matching Gmsh's `Curvature` name."""
const CurvatureField=CurvatureOpField
function _grad_norm(f,x,y,z,d,entity)
    g=((field_value(f,x+d/2,y,z,entity)-field_value(f,x-d/2,y,z,entity)),
       (field_value(f,x,y+d/2,z,entity)-field_value(f,x,y-d/2,z,entity)),
       (field_value(f,x,y,z+d/2,entity)-field_value(f,x,y,z-d/2,entity)))
    n=_norm3(g)
    n==0 && return nothing
    return ntuple(i->g[i]/n,3)
end
function _curvature_value(field::CurvatureOpField,x,y,z,entity)
    x=_float_value(x,"CurvatureField","x"); y=_float_value(y,"CurvatureField","y")
    z=_float_value(z,"CurvatureField","z"); d=field.delta; g=field.input
    gp=_grad_norm(g,x+d/2,y,z,d,entity); gm=_grad_norm(g,x-d/2,y,z,d,entity)
    hp=_grad_norm(g,x,y+d/2,z,d,entity); hm=_grad_norm(g,x,y-d/2,z,d,entity)
    kp=_grad_norm(g,x,y,z+d/2,d,entity); km=_grad_norm(g,x,y,z-d/2,d,entity)
    (gp===nothing || gm===nothing || hp===nothing || hm===nothing ||
     kp===nothing || km===nothing) && return GMSH_MAX_SIZE
    gp=gp::NTuple{3,Float64}; gm=gm::NTuple{3,Float64}
    hp=hp::NTuple{3,Float64}; hm=hm::NTuple{3,Float64}
    kp=kp::NTuple{3,Float64}; km=km::NTuple{3,Float64}
    value=(gp[1]-gm[1]+hp[2]-hm[2]+kp[3]-km[3])/d
    return _checked_field_result(value,"CurvatureField",x,y,z)
end
field_value(field::CurvatureOpField,x,y,z)=_curvature_value(field,x,y,z,nothing)
field_value(field::CurvatureOpField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_curvature_value(field,x,y,z,entity)

"""Gmsh `MaxEigenHessian`: largest eigenvalue of the Hessian of the input."""
struct MaxEigenHessianField{F<:AbstractField} <: AbstractField
    input::F; delta::Float64
end
MaxEigenHessianField(input::AbstractField; delta::Real=1e-4) =
    MaxEigenHessianField{typeof(input)}(input,_fd_delta(delta,"MaxEigenHessianField"))
function _maxeigenhessian_value(field::MaxEigenHessianField,x,y,z,entity)
    x=_float_value(x,"MaxEigenHessianField","x"); y=_float_value(y,"MaxEigenHessianField","y")
    z=_float_value(z,"MaxEigenHessianField","z"); d=field.delta; g=field.input
    f0=field_value(g,x,y,z,entity)
    hxx=field_value(g,x+d,y,z,entity)+field_value(g,x-d,y,z,entity)-2*f0
    hyy=field_value(g,x,y+d,z,entity)+field_value(g,x,y-d,z,entity)-2*f0
    hzz=field_value(g,x,y,z+d,entity)+field_value(g,x,y,z-d,entity)-2*f0
    hxy=field_value(g,x+d/2,y+d/2,z,entity)+field_value(g,x-d/2,y-d/2,z,entity)-
        field_value(g,x-d/2,y+d/2,z,entity)-field_value(g,x+d/2,y-d/2,z,entity)
    hxz=field_value(g,x+d/2,y,z+d/2,entity)+field_value(g,x-d/2,y,z-d/2,entity)-
        field_value(g,x-d/2,y,z+d/2,entity)-field_value(g,x+d/2,y,z-d/2,entity)
    hyz=field_value(g,x,y+d/2,z+d/2,entity)+field_value(g,x,y-d/2,z-d/2,entity)-
        field_value(g,x,y-d/2,z+d/2,entity)-field_value(g,x,y+d/2,z-d/2,entity)
    λ,_=_sym3_eigh(hxx,hyy,hzz,hxy,hxz,hyz)
    return _checked_field_result(maximum(λ)/(d*d),"MaxEigenHessianField",x,y,z)
end
field_value(field::MaxEigenHessianField,x,y,z)=
    _maxeigenhessian_value(field,x,y,z,nothing)
field_value(field::MaxEigenHessianField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_maxeigenhessian_value(field,x,y,z,entity)

# ── Coordinate maps ───────────────────────────────────────────────────────────

"""Gmsh `LonLat`: evaluate the input in geographic coordinates."""
struct LonLatField{F<:AbstractField} <: AbstractField
    input::F; from_stereo::Bool; radius::Float64
end
function LonLatField(input::AbstractField; from_stereo::Bool=false, radius::Real=6371e3)
    r=_positive_value(radius,"LonLatField","radius")
    return LonLatField{typeof(input)}(input,from_stereo,r)
end
function _lonlat_value(field::LonLatField,x,y,z,entity)
    x=_float_value(x,"LonLatField","x"); y=_float_value(y,"LonLatField","y")
    z=_float_value(z,"LonLatField","z")
    if field.from_stereo
        r2=field.radius*field.radius
        den=4*r2+x*x+y*y
        den==0 && return GMSH_MAX_SIZE
        X=4*r2*x/den; Y=4*r2*y/den
        Z=field.radius*(4*r2-x*x-y*y)/den
        x,y,z=X,Y,Z
    end
    az=z/field.radius
    # Gmsh's asin produces an invalid sample outside the sphere; represent the
    # resulting unconstrained field explicitly instead of leaking NaN.
    abs(az)<=1 || return GMSH_MAX_SIZE
    value=field_value(field.input,atan(y,x),asin(az),0.0,entity)
    return _checked_field_result(value,"LonLatField",x,y,z)
end
field_value(field::LonLatField,x,y,z)=_lonlat_value(field,x,y,z,nothing)
field_value(field::LonLatField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_lonlat_value(field,x,y,z,entity)

"""Gmsh `Param`: evaluate the input at `(FX(x,y,z), FY, FZ)`."""
struct ParametricField{F<:AbstractField,B<:Tuple} <: AbstractField
    input::F; fx::_MEProgram; fy::_MEProgram; fz::_MEProgram
    fields::B
end
@inline _parse_parametric_expression(expr::AbstractString)=
    isempty(expr) ? _me_compile(_MELit(GMSH_MAX_SIZE)) : parse_matheval(expr)
function ParametricField(input::AbstractField; fx::AbstractString="", fy::AbstractString="",
                         fz::AbstractString="", fields=Dict{Int,AbstractField}())
    bindings=_me_field_bindings(fields)
    return ParametricField{typeof(input),typeof(bindings)}(
        input,_parse_parametric_expression(fx),_parse_parametric_expression(fy),
        _parse_parametric_expression(fz),bindings)
end
function _parametric_value(field::ParametricField,x,y,z,entity)
    x=_float_value(x,"ParametricField","x"); y=_float_value(y,"ParametricField","y")
    z=_float_value(z,"ParametricField","z")
    X=_checked_field_result(_eval_matheval_entity(field.fx,x,y,z,field.fields,entity),
                            "ParametricField FX",x,y,z)
    Y=_checked_field_result(_eval_matheval_entity(field.fy,x,y,z,field.fields,entity),
                            "ParametricField FY",x,y,z)
    Z=_checked_field_result(_eval_matheval_entity(field.fz,x,y,z,field.fields,entity),
                            "ParametricField FZ",x,y,z)
    return _checked_field_result(field_value(field.input,X,Y,Z,entity),
                                 "ParametricField input",x,y,z)
end
field_value(field::ParametricField,x,y,z)=_parametric_value(field,x,y,z,nothing)
field_value(field::ParametricField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_parametric_value(field,x,y,z,entity)

# ── Structured grid ───────────────────────────────────────────────────────────

"""Gmsh `Structured`: trilinear interpolation on a regular 3-D grid."""
struct StructuredField <: AbstractField
    origin::NTuple{3,Float64}
    delta::NTuple{3,Float64}
    nx::Int; ny::Int; nz::Int
    data::Array{Float64,3}            # (nx, ny, nz) with v[i,j,k] = v(i-1,j-1,k-1)
    outside::Union{Nothing,Float64}
end
function StructuredField(origin, delta, data::AbstractArray{<:Real,3};
                         outside=nothing)
    o=_point3(origin,"StructuredField")
    draw=_point3(delta,"StructuredField")
    nx,ny,nz=size(data)
    nx>=1 && ny>=1 && nz>=1 || throw(ArgumentError("StructuredField: grid must be non-empty"))
    dims=(nx,ny,nz); names=("x","y","z")
    d=ntuple(3) do i
        draw[i]!=0 && return draw[i]
        dims[i]==1 && return 1.0
        throw(ArgumentError(
            "StructuredField: D$(names[i]) must be non-zero when n$(names[i]) > 1"))
    end
    grid=Array{Float64}(undef,nx,ny,nz)
    @inbounds for i in 1:nx, j in 1:ny, k in 1:nz
        v=Float64(data[i,j,k]); isfinite(v) ||
            throw(ArgumentError("StructuredField: non-finite sample at ($i,$j,$k)"))
        grid[i,j,k]=v
    end
    out=if outside===nothing
        nothing
    else
        vo=_float_value(outside,"StructuredField","outside"); vo
    end
    return StructuredField(o,d,nx,ny,nz,grid,out)
end
function _structured_count(nx::Int,ny::Int,nz::Int,max_samples::Int)
    (nx>0 && ny>0 && nz>0) || throw(ArgumentError("StructuredField: invalid grid size"))
    count=try
        Base.checked_mul(Base.checked_mul(nx,ny),nz)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("StructuredField: grid dimensions overflow the platform Int limit"))
    end
    count<=max_samples || throw(ArgumentError(
        "StructuredField: $count samples exceed max_samples=$max_samples"))
    return count
end

function _structured_dim(value::Float64,name::AbstractString)
    (isfinite(value) && value==trunc(value)) || throw(ArgumentError(
        "StructuredField: $name must be an integer"))
    return try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("StructuredField: $name exceeds the platform Int range"))
    end
end

function StructuredField(path::AbstractString; text::Bool=true, outside=nothing,
                         max_samples::Integer=10_000_000)
    isfile(path) || throw(ArgumentError("StructuredField: cannot open $path"))
    limit=try
        Int(max_samples)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("StructuredField: max_samples exceeds the platform Int limit"))
    end
    limit>0 || throw(ArgumentError("StructuredField: max_samples must be positive"))
    if text
        maxbytes=try
            Base.checked_add(1_048_576,Base.checked_mul(64,limit))
        catch err
            err isa InterruptException && rethrow()
            err isa OverflowError || rethrow()
            typemax(Int)
        end
        filesize(path)<=maxbytes || throw(ArgumentError(
            "StructuredField: ASCII file exceeds the configured resource limit"))
        nums=Float64[]
        open(path,"r") do io
            for raw in eachline(io)
                line=strip(first(split(raw,'#';limit=2)))
                isempty(line) && continue
                for tok in split(line)
                    v=tryparse(Float64,tok)
                    v===nothing && throw(ArgumentError("StructuredField: non-numeric token $tok"))
                    push!(nums,v)
                    length(nums)<=limit+9 || throw(ArgumentError(
                        "StructuredField: ASCII grid exceeds max_samples=$limit"))
                end
            end
        end
        length(nums)>=9 || throw(ArgumentError("StructuredField: truncated ASCII grid"))
        o=(nums[1],nums[2],nums[3]); d=(nums[4],nums[5],nums[6])
        nx=_structured_dim(nums[7],"nx"); ny=_structured_dim(nums[8],"ny")
        nz=_structured_dim(nums[9],"nz")
        count=_structured_count(nx,ny,nz,limit)
        length(nums)==9+count || throw(ArgumentError("StructuredField: sample count mismatch"))
        data=Array{Float64}(undef,nx,ny,nz)
        t=10
        # Gmsh stores k-fastest then j then i.
        @inbounds for i in 1:nx, j in 1:ny, k in 1:nz
            data[i,j,k]=nums[t]; t+=1
        end
        return StructuredField(o,d,data;outside=outside)
    end
    filesize(path)>=60 || throw(ArgumentError("StructuredField: truncated binary grid header"))
    try
        return open(path,"r") do io
            o=ntuple(_->read(io,Float64),3)
            d=ntuple(_->read(io,Float64),3)
            n=ntuple(_->Int(read(io,Int32)),3)
            count=_structured_count(n[1],n[2],n[3],limit)
            expected=try
                Base.checked_add(60,Base.checked_mul(8,count))
            catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError("StructuredField: binary grid byte count overflow"))
            end
            filesize(path)==expected || throw(ArgumentError(
                "StructuredField: binary file size mismatch (expected $expected bytes)"))
            data=Array{Float64}(undef,n)
            @inbounds for i in 1:n[1], j in 1:n[2], k in 1:n[3]
                data[i,j,k]=read(io,Float64)
            end
            return StructuredField(o,d,data;outside=outside)
        end
    catch err
        err isa InterruptException && rethrow()
        err isa EOFError && throw(ArgumentError("StructuredField: truncated binary grid"))
        rethrow()
    end
end
function field_value(field::StructuredField,x,y,z)
    x=_float_value(x,"StructuredField","x"); y=_float_value(y,"StructuredField","y")
    z=_float_value(z,"StructuredField","z")
    xyz=(x,y,z); id0=Vector{Int}(undef,3); id1=Vector{Int}(undef,3); ξ=Vector{Float64}(undef,3)
    @inbounds for i in 1:3
        t=(xyz[i]-field.origin[i])/field.delta[i]
        i0=floor(Int,t); i1=i0+1; n=(field.nx,field.ny,field.nz)[i]
        if field.outside!==nothing && n>1 && (i0<0 || i1>=n)
            return field.outside::Float64
        end
        i0=clamp(i0,0,n-1); i1=clamp(i1,0,n-1)
        id0[i]=i0+1; id1[i]=i1+1
        ξ[i]=clamp((xyz[i]-(field.origin[i]+i0*field.delta[i]))/field.delta[i],0.0,1.0)
    end
    v=0.0
    @inbounds for a in 0:1, b in 0:1, c in 0:1
        ia=a==0 ? id0[1] : id1[1]
        ib=b==0 ? id0[2] : id1[2]
        ic=c==0 ? id0[3] : id1[3]
        wa=a==0 ? 1-ξ[1] : ξ[1]
        wb=b==0 ? 1-ξ[2] : ξ[2]
        wc=c==0 ? 1-ξ[3] : ξ[3]
        v+=field.data[ia,ib,ic]*wa*wb*wc
    end
    return _checked_field_result(v,"StructuredField",x,y,z)
end

# ── Restrict / Constant (entity-aware) ────────────────────────────────────────

function _tagset(values,caller)
    s=Set{Int}()
    for v in values
        t=Int(v); t>0 || throw(ArgumentError("$caller: entity tags must be positive"))
        push!(s,t)
    end
    return s
end

function _topology_links(topology,key,caller)
    topology isa AbstractDict || throw(ArgumentError(
        "$caller: entity topology must be an AbstractDict"))
    raw=get(topology,key,())
    links=Tuple{Int,Int}[]
    for item in raw
        (item isa Tuple && length(item)==2 && item[1] isa Integer &&
         item[2] isa Integer) || throw(ArgumentError(
            "$caller: topology values must contain (dimension, tag) integer tuples"))
        dim=Int(item[1]); tag=Int(item[2])
        dim in 0:3 && tag>0 || throw(ArgumentError(
            "$caller: invalid topology entity ($dim,$tag)"))
        push!(links,(dim,tag))
    end
    return links
end

function _expanded_entity_sets(points,curves,surfaces,volumes,include_boundary,
                               include_embedded,boundaries,embedded,caller)
    sets=(_tagset(points,caller),_tagset(curves,caller),
          _tagset(surfaces,caller),_tagset(volumes,caller))
    seeds=Tuple{Int,Int}[]
    for dim in 0:3,tag in sets[dim+1]; push!(seeds,(dim,tag)); end
    function add_boundary!(root)
        queue=Tuple{Int,Int}[root]; seen=Set{Tuple{Int,Int}}()
        while !isempty(queue)
            key=pop!(queue); key in seen && continue; push!(seen,key)
            for child in _topology_links(boundaries,key,caller)
                push!(sets[child[1]+1],child[2]); push!(queue,child)
            end
        end
    end
    include_boundary && foreach(add_boundary!,seeds)
    if include_embedded
        for seed in seeds
            seed[1] in (2,3) || continue
            for child in _topology_links(embedded,seed,caller)
                child[1]<seed[1] || throw(ArgumentError(
                    "$caller: embedded entity $child must have lower dimension than $seed"))
                push!(sets[child[1]+1],child[2])
                # Gmsh always includes the bounding vertices of embedded curves
                # and the bounding curves/vertices of embedded surfaces,
                # independently of IncludeBoundary on the explicit seeds.
                child[1] in (1,2) && add_boundary!(child)
            end
        end
    end
    return sets
end

function _entity_hit(dim,tag,points,curves,surfaces,volumes,include_boundary)
    dim==0 && tag in points && return true
    dim==1 && tag in curves && return true
    dim==2 && tag in surfaces && return true
    dim==3 && tag in volumes && return true
    include_boundary || return false
    # Without a CAD adjacency graph, boundary inclusion is the caller's job via
    # passing the bounding entity tag. Geometric membership is handled separately.
    return false
end

"""Gmsh `Restrict`: apply `input` only on listed entity tags."""
struct RestrictField{F<:AbstractField} <: AbstractSizeField
    input::F
    points::Set{Int}; curves::Set{Int}; surfaces::Set{Int}; volumes::Set{Int}
    include_boundary::Bool; include_embedded::Bool
end
function RestrictField(input::AbstractField; points=(), curves=(), surfaces=(), volumes=(),
                       include_boundary::Bool=true,include_embedded::Bool=true,
                       entity_boundaries=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}(),
                       entity_embedded=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}())
    sets=_expanded_entity_sets(points,curves,surfaces,volumes,include_boundary,
        include_embedded,entity_boundaries,entity_embedded,"RestrictField")
    return RestrictField{typeof(input)}(input,sets...,include_boundary,include_embedded)
end
function field_value(field::RestrictField,x,y,z)
    return field_value(field.input,x,y,z)
end
field_value(field::RestrictField,x,y,z,::Nothing)=field_value(field.input,x,y,z)
function field_value(field::RestrictField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}
    dim=Int(entity[1]); tag=Int(entity[2])
    _entity_hit(dim,tag,field.points,field.curves,field.surfaces,field.volumes,
                field.include_boundary) && return field_value(field.input,x,y,z,entity)
    return GMSH_MAX_SIZE
end

"""Gmsh `Constant`: `VIn` on listed entities, `VOut` elsewhere."""
struct ConstantField <: AbstractSizeField
    vin::Float64; vout::Float64
    points::Set{Int}; curves::Set{Int}; surfaces::Set{Int}; volumes::Set{Int}
    include_boundary::Bool; include_embedded::Bool
end
function ConstantField(; vin::Real=GMSH_MAX_SIZE, vout::Real=GMSH_MAX_SIZE,
                       points=(), curves=(), surfaces=(), volumes=(),
                       include_boundary::Bool=true,include_embedded::Bool=true,
                       entity_boundaries=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}(),
                       entity_embedded=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}())
    sets=_expanded_entity_sets(points,curves,surfaces,volumes,include_boundary,
        include_embedded,entity_boundaries,entity_embedded,"ConstantField")
    return ConstantField(_float_value(vin,"ConstantField","vin"),
                         _float_value(vout,"ConstantField","vout"),
                         sets...,include_boundary,include_embedded)
end
function field_value(field::ConstantField,x,y,z)
    return GMSH_MAX_SIZE
end
field_value(field::ConstantField,x,y,z,::Nothing)=GMSH_MAX_SIZE
function field_value(field::ConstantField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}
    dim=Int(entity[1]); tag=Int(entity[2])
    _entity_hit(dim,tag,field.points,field.curves,field.surfaces,field.volumes,
                field.include_boundary) && return field.vin
    return field.vout
end

# ── Extend ────────────────────────────────────────────────────────────────────

struct _ExtendSamples
    points::Matrix{Float64}     # 3 × n
    sizes::Vector{Float64}
    tree::DistanceField
end

function _extend_samples(seeds,sizes,caller)
    pts=_point_matrix(seeds,"ExtendField")
    sz=Float64[ _positive_value(s,"ExtendField","size") for s in sizes ]
    size(pts,2)==length(sz) || throw(ArgumentError("ExtendField: seed and size counts differ"))
    size(pts,2)>0 || throw(ArgumentError("ExtendField: need at least one seed"))
    tree=DistanceField(;points=[(pts[1,j],pts[2,j],pts[3,j]) for j in 1:size(pts,2)])
    return _ExtendSamples(pts,sz,tree)
end

"""
    ExtendField(; curve_seeds, curve_sizes, surface_seeds, surface_sizes,
                dist_max=1, size_max=1, power=1, global_factor=1,
                entity_factors=Dict())

Gmsh `Extend`: extend sizes from already-meshed curve elements while meshing a
surface, and from already-meshed surface elements while meshing a volume. Samples
are boundary-element barycentres and their edge length or mean triangle-edge
length. Queries without a geometric entity, in dimensions 0/1, or without the
corresponding seed set return `GMSH_MAX_SIZE`.

Boundary sizes are divided by the target entity's factor and `global_factor`
before interpolation, matching Gmsh's later application of those factors.
The positional `ExtendField(seeds, sizes; ...)` convenience uses the same discrete
samples for both surface and volume queries; it still requires an entity context.
"""
struct ExtendField <: AbstractSizeField
    curves::Union{Nothing,_ExtendSamples}
    surfaces::Union{Nothing,_ExtendSamples}
    dist_max::Float64
    size_max::Float64
    power::Float64
    global_factor::Float64
    entity_factors::Dict{Tuple{Int,Int},Float64}
end

function ExtendField(;curve_seeds=nothing,curve_sizes=(),surface_seeds=nothing,
                     surface_sizes=(),dist_max::Real=1.0,size_max::Real=1.0,
                     power::Real=1.0,global_factor::Real=1.0,
                     entity_factors=Dict{Tuple{Int,Int},Float64}())
    curves=curve_seeds===nothing ? nothing :
           _extend_samples(curve_seeds,curve_sizes,"ExtendField")
    surfaces=surface_seeds===nothing ? nothing :
             _extend_samples(surface_seeds,surface_sizes,"ExtendField")
    dm=_positive_value(dist_max,"ExtendField","dist_max")
    sm=_positive_value(size_max,"ExtendField","size_max")
    pw=_float_value(power,"ExtendField","power")
    gf=_positive_value(global_factor,"ExtendField","global_factor")
    factors=Dict{Tuple{Int,Int},Float64}()
    entity_factors isa AbstractDict || throw(ArgumentError(
        "ExtendField: entity_factors must be a dictionary"))
    for (key,value) in entity_factors
        key isa Tuple && length(key)==2 && key[1] isa Integer && key[2] isa Integer ||
            throw(ArgumentError("ExtendField: entity-factor keys must be (dimension, tag) tuples"))
        dim=Int(key[1]);tag=Int(key[2])
        dim in 0:3 && tag>0 || throw(ArgumentError(
            "ExtendField: invalid entity-factor key $key"))
        factors[(dim,tag)]=_positive_value(value,"ExtendField","entity factor")
    end
    return ExtendField(curves,surfaces,dm,sm,pw,gf,factors)
end

function ExtendField(seeds,sizes;kwargs...)
    data=_extend_samples(seeds,sizes,"ExtendField")
    field=ExtendField(;kwargs...)
    return ExtendField(data,data,field.dist_max,field.size_max,field.power,
                       field.global_factor,field.entity_factors)
end

field_value(::ExtendField,x,y,z)=GMSH_MAX_SIZE
field_value(::ExtendField,x,y,z,::Nothing)=GMSH_MAX_SIZE
function field_value(field::ExtendField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}
    p=(_float_value(x,"ExtendField","x"),_float_value(y,"ExtendField","y"),
       _float_value(z,"ExtendField","z"))
    dim=Int(entity[1]);tag=Int(entity[2])
    samples=dim==2 ? field.curves : dim==3 ? field.surfaces : nothing
    samples===nothing && return GMSH_MAX_SIZE
    best,ibest=_nearest_point(samples.tree,p)
    best>=field.dist_max && return field.size_max
    f=(field.dist_max-best)/field.dist_max
    field.power!=1 && (f=f^field.power)
    factor=get(field.entity_factors,(dim,tag),1.0)*field.global_factor
    boundary_size=samples.sizes[ibest]/factor
    value=f*boundary_size+(1-f)*field.size_max
    return _checked_field_result(value,"ExtendField",p...;positive=true)
end

# ── Octree acceleration ───────────────────────────────────────────────────────

mutable struct _OctCell
    leaf::Bool
    value::Float64
    children::Vector{_OctCell}
end
function _oct_build(field,x0,y0,z0,l,level,budget::Base.RefValue{Int})
    budget[]>0 || throw(ArgumentError(
        "OctreeField: adaptive tree exceeds max_cells resource limit"))
    budget[]-=1
    vc=_checked_field_result(field_value(field,x0+l/2,y0+l/2,z0+l/2),
                             "OctreeField input",x0+l/2,y0+l/2,z0+l/2)
    split=level>0
    if level>-4 && !split
        dl=l/2
        dmax=0.0; vmin=vc
        @inbounds for i in 0:2, j in 0:2, k in 0:2
            sx=x0+i*dl; sy=y0+j*dl; sz=z0+k*dl
            w=_checked_field_result(field_value(field,sx,sy,sz),
                                    "OctreeField input",sx,sy,sz)
            dmax=max(dmax,abs(vc-w)); vmin=min(vmin,w)
            if dmax/vmin>0.2 && vmin<l
                split=true; break
            end
        end
    end
    if split && level> -4
        kids=_OctCell[]
        l2=l/2
        for (dx,dy,dz) in ((0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1))
            push!(kids,_oct_build(field,x0+dx*l2,y0+dy*l2,z0+dz*l2,l2,level-1,budget))
        end
        return _OctCell(false,vc,kids)
    end
    return _OctCell(true,vc,_OctCell[])
end
function _oct_eval(cell::_OctCell,x,y,z)
    cell.leaf && return cell.value
    i=x>0.5 ? 1 : 0; j=y>0.5 ? 1 : 0; k=z>0.5 ? 1 : 0
    return _oct_eval(cell.children[i*4+j*2+k+1],2x-i,2y-j,2z-k)
end

"""Gmsh `Octree`: precompute another field on an adaptive octree over `bbox`."""
struct OctreeField{F<:AbstractField} <: AbstractField
    input::F
    origin::NTuple{3,Float64}
    length::Float64
    root::_OctCell
end
function OctreeField(input::AbstractField, xmin,xmax,ymin,ymax,zmin,zmax;
                     max_level::Integer=4,max_cells::Integer=1_000_000)
    xa=_float_value(xmin,"OctreeField","xmin"); xb=_float_value(xmax,"OctreeField","xmax")
    ya=_float_value(ymin,"OctreeField","ymin"); yb=_float_value(ymax,"OctreeField","ymax")
    za=_float_value(zmin,"OctreeField","zmin"); zb=_float_value(zmax,"OctreeField","zmax")
    (xb>=xa && yb>=ya && zb>=za) || throw(ArgumentError(
        "OctreeField: bounding-box maxima must not be below minima"))
    l=max(xb-xa,yb-ya,zb-za)
    (isfinite(l) && l>0) || throw(ArgumentError("OctreeField: empty bounding box"))
    lv=try
        Int(max_level)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("OctreeField: max_level exceeds the platform Int range"))
    end
    lv>=0 || throw(ArgumentError("OctreeField: max_level must be non-negative"))
    lv<=32 || throw(ArgumentError("OctreeField: max_level must be at most 32"))
    cells=try
        Int(max_cells)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("OctreeField: max_cells exceeds the platform Int range"))
    end
    cells>0 || throw(ArgumentError("OctreeField: max_cells must be positive"))
    budget=Ref(cells)
    root=_oct_build(input,xa,ya,za,l,lv,budget)
    return OctreeField{typeof(input)}(input,(xa,ya,za),l,root)
end
function field_value(field::OctreeField,x,y,z)
    x=_float_value(x,"OctreeField","x"); y=_float_value(y,"OctreeField","y")
    z=_float_value(z,"OctreeField","z")
    xn=(x-field.origin[1])/field.length
    yn=(y-field.origin[2])/field.length
    zn=(z-field.origin[3])/field.length
    return _oct_eval(field.root,xn,yn,zn)
end

# ── PostView ──────────────────────────────────────────────────────────────────

"""
    PostViewField(coords, values; points=nothing, lines=(), triangles=(),
                  quadrangles=(), tetrahedra=(), hexahedra=(), prisms=(),
                  pyramids=(), crop_negative=true, use_closest=true,
                  reference_tolerance=1e-6, max_nodes=1_000_000,
                  max_elements=1_000_000)

Gmsh `PostView` on first-order list data. `coords` is `3×n`. Scalar `values`
can be a length-`n` vector or a `1×n` matrix; vector and tensor values are `3×n`
and `9×n` matrices; multiple time steps pass a `c×n×steps` array and select one
with `time` (all steps are validated, the selected step is stored). Connectivity
is 1-based and can be an
arity-by-element matrix or an iterable of index tuples. When no connectivity is
supplied, every column remains a point element, preserving the original API.
With explicit topology, point elements are absent unless supplied through
`points`.

Queries use Gmsh's first-order reference-element interpolants for points, lines,
triangles, planar quadrangles, tetrahedra, hexahedra, prisms and pyramids. If no
element contains the query, `use_closest=true` retries at the closest active
node. Scalar data returns the interpolant; vector data returns the Euclidean norm
of the interpolated vector, matching the scalar operator of Gmsh's `PostView`
field. The scalar operator of tensor data is [`GMSH_MAX_SIZE`](@ref), as Gmsh
does; tensor-to-metric evaluation is not exposed by this scalar constructor.
Non-positive scalar results become `GMSH_MAX_SIZE` when `crop_negative=true`.

Curved/high-order list interpolation matrices, mixed component
counts, and non-planar quadrangles are rejected or not represented. Gmsh itself
requires an adapted visualization grid before querying non-adapted high-order list
data. Vector norms use overflow-safe arithmetic; results that are not representable
as finite `Float64` values are rejected by the field-result contract.
"""
struct PostViewField <: AbstractField
    coords::Matrix{Float64}
    values::Vector{Float64}
    num_components::UInt8
    points::Matrix{Int32}
    lines::Matrix{Int32}
    triangles::Matrix{Int32}
    quadrangles::Matrix{Int32}
    tetrahedra::Matrix{Int32}
    hexahedra::Matrix{Int32}
    prisms::Matrix{Int32}
    pyramids::Matrix{Int32}
    crop_negative::Bool
    use_closest::Bool
    reference_tolerance::Float64
    tree::DistanceField
    bvh_order::Vector{Int}
    bvh_lo::Vector{NTuple{3,Float64}}
    bvh_hi::Vector{NTuple{3,Float64}}
    bvh_left::Vector{Int}
    bvh_right::Vector{Int}
    bvh_first::Vector{Int}
    bvh_count::Vector{Int}
end

const _POSTVIEW_BVH_LEAF_SIZE=8

function _postview_limit(value::Integer,name::AbstractString)
    value isa Bool && throw(ArgumentError("PostViewField: $name must be an integer"))
    limit=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("PostViewField: $name exceeds the platform Int range"))
    end
    limit>0 || throw(ArgumentError("PostViewField: $name must be positive"))
    return limit
end

@inline function _postview_index(value,n::Int,what::AbstractString)
    (value isa Integer && !(value isa Bool)) || throw(ArgumentError(
        "PostViewField: $what indices must be integers"))
    index=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("PostViewField: $what index $value exceeds the platform Int range"))
    end
    1<=index<=n || throw(ArgumentError(
        "PostViewField: $what index $value is outside 1:$n"))
    return Int32(index)
end

function _postview_connectivity(input,arity::Int,n::Int,what::AbstractString,
                                remaining::Int)
    if input isa AbstractMatrix
        size(input,1)==arity || throw(ArgumentError(
            "PostViewField: $what connectivity must have $arity rows"))
        count=size(input,2)
        count<=remaining || throw(ArgumentError(
            "PostViewField: element count exceeds max_elements"))
        result=Matrix{Int32}(undef,arity,count)
        row_axes=axes(input,1); column_axes=axes(input,2)
        @inbounds for (j,source_j) in enumerate(column_axes),
                      (i,source_i) in enumerate(row_axes)
            result[i,j]=_postview_index(input[source_i,source_j],n,what)
        end
        return result
    end
    count=try
        length(input)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "PostViewField: $what connectivity must be a matrix or a finite iterable with length"))
    end
    count>=0 || throw(ArgumentError(
        "PostViewField: $what connectivity reported a negative length"))
    count<=remaining || throw(ArgumentError(
        "PostViewField: element count exceeds max_elements"))
    result=Matrix{Int32}(undef,arity,count)
    seen=0
    for cell in input
        seen+=1
        seen<=count || throw(ArgumentError(
            "PostViewField: $what connectivity yielded more entries than its length"))
        if arity==1 && cell isa Integer
            result[1,seen]=_postview_index(cell,n,what)
        else
            cell_length=try
                length(cell)
            catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError(
                    "PostViewField: each $what entry must contain $arity indices"))
            end
            cell_length==arity || throw(ArgumentError(
                "PostViewField: each $what entry must contain $arity indices"))
            for i in 1:arity
                result[i,seen]=_postview_index(cell[i],n,what)
            end
        end
    end
    seen==count || throw(ArgumentError(
        "PostViewField: $what connectivity yielded fewer entries than its length"))
    return result
end

@inline _postview_point(coords,i::Integer)=
    (coords[1,i],coords[2,i],coords[3,i])

@inline function _postview_cross(a,b)
    return (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],
            a[1]*b[2]-a[2]*b[1])
end

@inline _postview_det(a,b,c)=_dot3(a,_postview_cross(b,c))

@inline function _postview_shape(kind::UInt8,slot::Int,u,v,w)
    if kind==7 # quadrangle
        slot==1 && return 0.25*(1-u)*(1-v)
        slot==2 && return 0.25*(1+u)*(1-v)
        slot==3 && return 0.25*(1+u)*(1+v)
        return 0.25*(1-u)*(1+v)
    elseif kind==4 # hexahedron
        slot==1 && return 0.125*(1-u)*(1-v)*(1-w)
        slot==2 && return 0.125*(1+u)*(1-v)*(1-w)
        slot==3 && return 0.125*(1+u)*(1+v)*(1-w)
        slot==4 && return 0.125*(1-u)*(1+v)*(1-w)
        slot==5 && return 0.125*(1-u)*(1-v)*(1+w)
        slot==6 && return 0.125*(1+u)*(1-v)*(1+w)
        slot==7 && return 0.125*(1+u)*(1+v)*(1+w)
        return 0.125*(1-u)*(1+v)*(1+w)
    elseif kind==5 # prism
        slot==1 && return 0.5*(1-u-v)*(1-w)
        slot==2 && return 0.5*u*(1-w)
        slot==3 && return 0.5*v*(1-w)
        slot==4 && return 0.5*(1-u-v)*(1+w)
        slot==5 && return 0.5*u*(1+w)
        return 0.5*v*(1+w)
    elseif kind==6 # pyramid
        r=(w!=1 && slot!=5) ? u*v*w/(1-w) : zero(u+v+w)
        slot==1 && return 0.25*((1-u)*(1-v)-w+r)
        slot==2 && return 0.25*((1+u)*(1-v)-w-r)
        slot==3 && return 0.25*((1+u)*(1+v)-w+r)
        slot==4 && return 0.25*((1-u)*(1+v)-w-r)
        return w
    end
    throw(ErrorException("PostViewField: internal unsupported shape kind $kind"))
end

@inline function _postview_shape_gradient(kind::UInt8,slot::Int,u,v,w)
    if kind==7
        slot==1 && return (-0.25*(1-v),-0.25*(1-u),0.0)
        slot==2 && return ( 0.25*(1-v),-0.25*(1+u),0.0)
        slot==3 && return ( 0.25*(1+v), 0.25*(1+u),0.0)
        return (-0.25*(1+v),0.25*(1-u),0.0)
    elseif kind==4
        slot==1 && return (-0.125*(1-v)*(1-w),-0.125*(1-u)*(1-w),-0.125*(1-u)*(1-v))
        slot==2 && return ( 0.125*(1-v)*(1-w),-0.125*(1+u)*(1-w),-0.125*(1+u)*(1-v))
        slot==3 && return ( 0.125*(1+v)*(1-w), 0.125*(1+u)*(1-w),-0.125*(1+u)*(1+v))
        slot==4 && return (-0.125*(1+v)*(1-w), 0.125*(1-u)*(1-w),-0.125*(1-u)*(1+v))
        slot==5 && return (-0.125*(1-v)*(1+w),-0.125*(1-u)*(1+w), 0.125*(1-u)*(1-v))
        slot==6 && return ( 0.125*(1-v)*(1+w),-0.125*(1+u)*(1+w), 0.125*(1+u)*(1-v))
        slot==7 && return ( 0.125*(1+v)*(1+w), 0.125*(1+u)*(1+w), 0.125*(1+u)*(1+v))
        return (-0.125*(1+v)*(1+w),0.125*(1-u)*(1+w),0.125*(1-u)*(1+v))
    elseif kind==5
        slot==1 && return (-0.5*(1-w),-0.5*(1-w),-0.5*(1-u-v))
        slot==2 && return ( 0.5*(1-w),0.0,-0.5*u)
        slot==3 && return (0.0,0.5*(1-w),-0.5*v)
        slot==4 && return (-0.5*(1+w),-0.5*(1+w),0.5*(1-u-v))
        slot==5 && return ( 0.5*(1+w),0.0,0.5*u)
        return (0.0,0.5*(1+w),0.5*v)
    elseif kind==6
        if w==1
            slot==1 && return (-0.25,-0.25,-0.25)
            slot==2 && return ( 0.25,-0.25,-0.25)
            slot==3 && return ( 0.25, 0.25,-0.25)
            slot==4 && return (-0.25, 0.25,-0.25)
            return (0.0,0.0,1.0)
        end
        denominator=1-w
        rr=u*v/denominator+u*v*w/(denominator*denominator)
        slot==1 && return (0.25*(-(1-v)+v*w/denominator),
                            0.25*(-(1-u)+u*w/denominator),0.25*(-1+rr))
        slot==2 && return (0.25*((1-v)-v*w/denominator),
                            0.25*(-(1+u)-u*w/denominator),0.25*(-1-rr))
        slot==3 && return (0.25*((1+v)+v*w/denominator),
                            0.25*((1+u)+u*w/denominator),0.25*(-1+rr))
        slot==4 && return (0.25*(-(1+v)-v*w/denominator),
                            0.25*((1-u)-u*w/denominator),0.25*(-1-rr))
        return (0.0,0.0,1.0)
    end
    throw(ErrorException("PostViewField: internal unsupported gradient kind $kind"))
end

@inline function _postview_relative(coords,index::Integer,base::Integer,d::Int,
                                    scale::Float64)
    delta=coords[d,index]-coords[d,base]
    return isfinite(delta) ? delta/scale :
           coords[d,index]/scale-coords[d,base]/scale
end

function _postview_cell_scale(coords,cells,j,arity::Int)
    base=cells[1,j]; scale=0.0
    @inbounds for slot in 2:arity, d in 1:3
        delta=coords[d,cells[slot,j]]-coords[d,base]
        if isfinite(delta)
            scale=max(scale,abs(delta))
        else
            scale=max(scale,abs(coords[d,cells[slot,j]]),abs(coords[d,base]))
        end
    end
    return scale
end

@inline function _postview_reference_jacobian(coords,cells,j,arity::Int,
                                               kind::UInt8,u,v,w,scale)
    du=(0.0,0.0,0.0);dv=(0.0,0.0,0.0);dw=(0.0,0.0,0.0)
    base=cells[1,j]
    @inbounds for slot in 1:arity
        gradient=_postview_shape_gradient(kind,slot,u,v,w)
        index=cells[slot,j]
        q=(_postview_relative(coords,index,base,1,scale),
           _postview_relative(coords,index,base,2,scale),
           _postview_relative(coords,index,base,3,scale))
        du=(muladd(q[1],gradient[1],du[1]),muladd(q[2],gradient[1],du[2]),
            muladd(q[3],gradient[1],du[3]))
        dv=(muladd(q[1],gradient[2],dv[1]),muladd(q[2],gradient[2],dv[2]),
            muladd(q[3],gradient[2],dv[3]))
        dw=(muladd(q[1],gradient[3],dw[1]),muladd(q[2],gradient[3],dw[2]),
            muladd(q[3],gradient[3],dw[3]))
    end
    if kind==7
        a=( _postview_relative(coords,cells[2,j],base,1,scale),
            _postview_relative(coords,cells[2,j],base,2,scale),
            _postview_relative(coords,cells[2,j],base,3,scale))
        b=( _postview_relative(coords,cells[3,j],base,1,scale),
            _postview_relative(coords,cells[3,j],base,2,scale),
            _postview_relative(coords,cells[3,j],base,3,scale))
        normal=_postview_cross(a,b)
        physical_scale=(normal[1]*scale,normal[2]*scale,normal[3]*scale)
        # Gmsh inserts the physical (length-squared) surface normal as the
        # artificial third Jacobian row. Restore that relative scaling after
        # normalizing coordinates; for extreme data, retain the finite
        # normalized normal instead of manufacturing an infinite Jacobian.
        usable_physical_scale=all(isfinite,physical_scale) &&
            (!iszero(physical_scale[1]) || !iszero(physical_scale[2]) ||
             !iszero(physical_scale[3]))
        dw=usable_physical_scale ? physical_scale : normal
    end
    return du,dv,dw
end

function _postview_triangle_nondegenerate(coords,cell,j)
    a=_postview_point(coords,cell[1,j]);b=_postview_point(coords,cell[2,j])
    c=_postview_point(coords,cell[3,j])
    return _postview_orient2_nonzero(a[1],a[2],b[1],b[2],c[1],c[2]) ||
           _postview_orient2_nonzero(a[1],a[3],b[1],b[3],c[1],c[3]) ||
           _postview_orient2_nonzero(a[2],a[3],b[2],b[3],c[2],c[3])
end

@inline function _postview_orient2_nonzero(ax,ay,bx,by,cx,cy)
    orient2((ax,ay),(bx,by),(cx,cy))!=0 && return true
    # `orient2` deliberately optimizes zero-product cases. At the extreme ends
    # of Float64, however, a finite subtraction can overflow and turn 0*Inf
    # into NaN; confirm an apparent zero exactly before rejecting list data.
    axr=_metric_rational(ax);ayr=_metric_rational(ay)
    bxr=_metric_rational(bx);byr=_metric_rational(by)
    cxr=_metric_rational(cx);cyr=_metric_rational(cy)
    return (bxr-axr)*(cyr-ayr)-(byr-ayr)*(cxr-axr)!=0
end

function _postview_tetrahedron_nondegenerate(coords,cell,j)
    a=_postview_point(coords,cell[1,j]);b=_postview_point(coords,cell[2,j])
    c=_postview_point(coords,cell[3,j]);d=_postview_point(coords,cell[4,j])
    return orient3(a,b,c,d)!=0
end

function _postview_unique_cell(cell,j,arity::Int,what::AbstractString)
    @inbounds for a in 1:arity-1, b in a+1:arity
        cell[a,j]!=cell[b,j] || throw(ArgumentError(
            "PostViewField: $what $j repeats a node index"))
    end
    return nothing
end

function _postview_nonlinear_nondegenerate(coords,cell,j,arity::Int,kind::UInt8)
    scale=_postview_cell_scale(coords,cell,j,arity)
    scale>0 || return false
    du,dv,dw=_postview_reference_jacobian(
        coords,cell,j,arity,kind,0.0,0.0,0.0,scale)
    determinant=_postview_det(du,dv,dw)
    return isfinite(determinant) && !iszero(determinant)
end

function _postview_quadrangle_planar(coords,cell,j)
    scale=_postview_cell_scale(coords,cell,j,4)
    scale>0 || return false
    base=cell[1,j]
    a=ntuple(d->_postview_relative(coords,cell[2,j],base,d,scale),3)
    b=ntuple(d->_postview_relative(coords,cell[3,j],base,d,scale),3)
    c=ntuple(d->_postview_relative(coords,cell[4,j],base,d,scale),3)
    normal=_postview_cross(a,b);normal_length=hypot(normal...)
    (isfinite(normal_length) && normal_length>0) || return false
    distance=abs(_dot3(c,normal))/normal_length
    # Affine construction in Float64 can leave an exact-predicate residue of a
    # few ulps. Accept only that rounding envelope, not geometrically warped
    # quadrangles whose Gmsh reference lookup is not reliable.
    return isfinite(distance) && distance<=128eps(Float64)
end

function _validate_postview_cells(coords,lines,triangles,quadrangles,
                                  tetrahedra,hexahedra,prisms,pyramids)
    @inbounds for j in axes(lines,2)
        a=lines[1,j];b=lines[2,j]
        a!=b && _postview_point(coords,a)!=_postview_point(coords,b) ||
            throw(ArgumentError("PostViewField: line $j is degenerate"))
    end
    @inbounds for j in axes(triangles,2)
        a=triangles[1,j];b=triangles[2,j];c=triangles[3,j]
        (a!=b && a!=c && b!=c) || throw(ArgumentError(
            "PostViewField: triangle $j repeats a node index"))
        _postview_triangle_nondegenerate(coords,triangles,j) || throw(ArgumentError(
            "PostViewField: triangle $j is degenerate"))
    end
    @inbounds for j in axes(quadrangles,2)
        _postview_unique_cell(quadrangles,j,4,"quadrangle")
        _postview_quadrangle_planar(coords,quadrangles,j) || throw(ArgumentError(
            "PostViewField: quadrangle $j is non-planar; warped first-order list quadrangles are unsupported"))
        _postview_nonlinear_nondegenerate(coords,quadrangles,j,4,UInt8(7)) ||
            throw(ArgumentError("PostViewField: quadrangle $j is degenerate"))
    end
    @inbounds for j in axes(tetrahedra,2)
        a=tetrahedra[1,j];b=tetrahedra[2,j]
        c=tetrahedra[3,j];d=tetrahedra[4,j]
        (a!=b && a!=c && a!=d && b!=c && b!=d && c!=d) || throw(ArgumentError(
            "PostViewField: tetrahedron $j repeats a node index"))
        _postview_tetrahedron_nondegenerate(coords,tetrahedra,j) ||
            throw(ArgumentError("PostViewField: tetrahedron $j is degenerate"))
    end
    @inbounds for (cells,arity,kind,what) in
            ((hexahedra,8,UInt8(4),"hexahedron"),
             (prisms,6,UInt8(5),"prism"),(pyramids,5,UInt8(6),"pyramid"))
        for j in axes(cells,2)
            _postview_unique_cell(cells,j,arity,what)
            _postview_nonlinear_nondegenerate(coords,cells,j,arity,kind) ||
                throw(ArgumentError("PostViewField: $what $j is degenerate"))
        end
    end
    return nothing
end

@inline function _postview_cell_kind(field,id::Int)
    nt=size(field.tetrahedra,2);nh=size(field.hexahedra,2)
    ni=size(field.prisms,2);ny=size(field.pyramids,2)
    ntri=size(field.triangles,2);nq=size(field.quadrangles,2);nl=size(field.lines,2)
    id<=nt && return UInt8(3),id
    id-=nt;id<=nh && return UInt8(4),id
    id-=nh;id<=ni && return UInt8(5),id
    id-=ni;id<=ny && return UInt8(6),id
    id-=ny;id<=ntri && return UInt8(2),id
    id-=ntri;id<=nq && return UInt8(7),id
    id-=nq;id<=nl && return UInt8(1),id
    return UInt8(0),id-nl
end

@inline function _postview_cell_indices(tetrahedra,hexahedra,prisms,pyramids,
                                        triangles,quadrangles,lines,points,id::Int)
    nt=size(tetrahedra,2);nh=size(hexahedra,2);ni=size(prisms,2)
    ny=size(pyramids,2);ntri=size(triangles,2);nq=size(quadrangles,2)
    nl=size(lines,2)
    if id<=nt
        return tetrahedra,id,4
    end
    id-=nt;id<=nh && return hexahedra,id,8
    id-=nh;id<=ni && return prisms,id,6
    id-=ni;id<=ny && return pyramids,id,5
    id-=ny;id<=ntri && return triangles,id,3
    id-=ntri;id<=nq && return quadrangles,id,4
    id-=nq;id<=nl && return lines,id,2
    return points,id-nl,1
end

@inline function _postview_saturating_sub(x::Float64,y::Float64)
    value=x-y
    return isfinite(value) ? value : -floatmax(Float64)
end

@inline function _postview_saturating_add(x::Float64,y::Float64)
    value=x+y
    return isfinite(value) ? value : floatmax(Float64)
end

function _build_postview_bvh(coords,tetrahedra,hexahedra,prisms,pyramids,
                             triangles,quadrangles,lines,points)
    ncells=sum(size(cells,2) for cells in
        (tetrahedra,hexahedra,prisms,pyramids,triangles,quadrangles,lines,points))
    primitive_lo=Matrix{Float64}(undef,3,ncells)
    primitive_hi=Matrix{Float64}(undef,3,ncells)
    centroid=Matrix{Float64}(undef,3,ncells)
    @inbounds for id in 1:ncells
        cells,j,arity=_postview_cell_indices(tetrahedra,hexahedra,prisms,pyramids,
                                             triangles,quadrangles,lines,points,id)
        for k in 1:3
            lo=Inf;hi=-Inf
            for slot in 1:arity
                value=coords[k,cells[slot,j]]
                lo=min(lo,value);hi=max(hi,value)
            end
            primitive_lo[k,id]=lo;primitive_hi[k,id]=hi
        end
        dx=primitive_hi[1,id]-primitive_lo[1,id]
        dy=primitive_hi[2,id]-primitive_lo[2,id]
        dz=primitive_hi[3,id]-primitive_lo[3,id]
        diagonal=hypot(dx,dy,dz)
        padding=isfinite(diagonal) ? 0.01diagonal : floatmax(Float64)
        for k in 1:3
            lo=_postview_saturating_sub(primitive_lo[k,id],padding)
            hi=_postview_saturating_add(primitive_hi[k,id],padding)
            primitive_lo[k,id]=lo;primitive_hi[k,id]=hi
            centroid[k,id]=signbit(lo)==signbit(hi) ? lo+(hi-lo)/2 : lo/2+hi/2
        end
    end
    order=collect(1:ncells)
    lows=NTuple{3,Float64}[];highs=NTuple{3,Float64}[]
    left=Int[];right=Int[];firsts=Int[];counts=Int[]
    function range_bounds(first,last)
        lo=(Inf,Inf,Inf);hi=(-Inf,-Inf,-Inf)
        @inbounds for pos in first:last
            id=order[pos]
            lo=ntuple(k->min(lo[k],primitive_lo[k,id]),3)
            hi=ntuple(k->max(hi[k],primitive_hi[k,id]),3)
        end
        return lo,hi
    end
    function build!(first,last)
        node=length(left)+1
        push!(lows,(0.0,0.0,0.0));push!(highs,(0.0,0.0,0.0))
        push!(left,0);push!(right,0);push!(firsts,0);push!(counts,0)
        count=last-first+1
        if count<=_POSTVIEW_BVH_LEAF_SIZE
            lows[node],highs[node]=range_bounds(first,last)
            firsts[node]=first;counts[node]=count
            return node
        end
        clo=(Inf,Inf,Inf);chi=(-Inf,-Inf,-Inf)
        @inbounds for pos in first:last
            id=order[pos]
            clo=ntuple(k->min(clo[k],centroid[k,id]),3)
            chi=ntuple(k->max(chi[k],centroid[k,id]),3)
        end
        ext=ntuple(k->chi[k]-clo[k],3)
        axis=ext[1]>=ext[2] ? (ext[1]>=ext[3] ? 1 : 3) :
             (ext[2]>=ext[3] ? 2 : 3)
        sort!(view(order,first:last);by=id->centroid[axis,id],alg=QuickSort)
        mid=(first+last)>>>1
        l=build!(first,mid);r=build!(mid+1,last)
        left[node]=l;right[node]=r
        lows[node]=ntuple(k->min(lows[l][k],lows[r][k]),3)
        highs[node]=ntuple(k->max(highs[l][k],highs[r][k]),3)
        return node
    end
    build!(1,ncells)
    return order,lows,highs,left,right,firsts,counts
end

function _postview_values(values,n::Int,time::Int)
    if values isa AbstractVector
        length(values)==n || throw(ArgumentError("PostViewField: value count mismatch"))
        time==1 || throw(ArgumentError(
            "PostViewField: time step $time is outside the single-step data"))
        return _postview_finite(Matrix{Float64}(reshape(values,1,n))),1
    elseif values isa AbstractMatrix
        size(values,2)==n || throw(ArgumentError("PostViewField: value count mismatch"))
        size(values,1) in (1,3,9) || throw(ArgumentError(
            "PostViewField: values must have 1, 3, or 9 components per node"))
        time==1 || throw(ArgumentError(
            "PostViewField: time step $time is outside the single-step data"))
        return _postview_finite(Matrix{Float64}(values)),1
    elseif values isa AbstractArray && ndims(values)==3
        size(values,2)==n || throw(ArgumentError("PostViewField: value count mismatch"))
        size(values,1) in (1,3,9) || throw(ArgumentError(
            "PostViewField: values must have 1, 3, or 9 components per node"))
        steps=size(values,3)
        steps<=typemax(Int32) || throw(ArgumentError(
            "PostViewField: time-step count exceeds the Int32 limit"))
        steps>0 || throw(ArgumentError("PostViewField: empty time-step array"))
        1<=time<=steps || throw(ArgumentError(
            "PostViewField: time step $time is outside 1:$steps"))
        # Validate every step once; only the selected step is retained.
        full=try
            Matrix{Float64}(reshape(values,size(values,1),n*steps))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("PostViewField: values must be Float64-representable"))
        end
        all(isfinite,full) || throw(ArgumentError("PostViewField: values must be finite"))
        lo=(time-1)*n+1
        return full[:,lo:lo+n-1],steps
    else
        throw(ArgumentError(
            "PostViewField: values must be a vector, a c×n matrix, or a c×n×steps array"))
    end
end

@inline function _postview_finite(result::Matrix{Float64})
    all(isfinite,result) || throw(ArgumentError("PostViewField: values must be finite"))
    return result
end

function PostViewField(coords::AbstractMatrix{<:Real}, values;
                       points=nothing,lines=(),triangles=(),quadrangles=(),
                       tetrahedra=(),hexahedra=(),prisms=(),pyramids=(),
                       time::Integer=1,
                       crop_negative::Bool=true,use_closest::Bool=true,
                       reference_tolerance::Real=1e-6,
                       max_nodes::Integer=1_000_000,
                       max_elements::Integer=1_000_000)
    size(coords,1)==3 || throw(ArgumentError("PostViewField: coords must be 3×n"))
    size(coords,2)>0 || throw(ArgumentError("PostViewField: empty view"))
    node_limit=_postview_limit(max_nodes,"max_nodes")
    element_limit=_postview_limit(max_elements,"max_elements")
    n=size(coords,2)
    n<=node_limit || throw(ArgumentError(
        "PostViewField: node count $n exceeds max_nodes=$node_limit"))
    n<=typemax(Int32) || throw(ArgumentError(
        "PostViewField: node count $n exceeds the Int32 connectivity limit"))
    C=try
        Matrix{Float64}(coords)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("PostViewField: coordinates must be Float64-representable"))
    end
    time==0 && throw(ArgumentError("PostViewField: time step 0 is outside range"))
    tsel=Int(time); tsel>=1 || throw(ArgumentError(
        "PostViewField: negative time step $time"))
    V,_steps=_postview_values(values,n,tsel)
    @inbounds for i in axes(C,2), d in 1:3
        isfinite(C[d,i]) || throw(ArgumentError("PostViewField: non-finite coordinate"))
    end
    tolerance=_float_value(reference_tolerance,"PostViewField","reference_tolerance")
    tolerance>=0 || throw(ArgumentError(
        "PostViewField: reference_tolerance must be non-negative"))

    remaining=element_limit
    L=_postview_connectivity(lines,2,n,"line",remaining);remaining-=size(L,2)
    T=_postview_connectivity(triangles,3,n,"triangle",remaining);remaining-=size(T,2)
    Q=_postview_connectivity(quadrangles,4,n,"quadrangle",remaining);remaining-=size(Q,2)
    S=_postview_connectivity(tetrahedra,4,n,"tetrahedron",remaining);remaining-=size(S,2)
    H=_postview_connectivity(hexahedra,8,n,"hexahedron",remaining);remaining-=size(H,2)
    I=_postview_connectivity(prisms,6,n,"prism",remaining);remaining-=size(I,2)
    Y=_postview_connectivity(pyramids,5,n,"pyramid",remaining);remaining-=size(Y,2)
    explicit_topology=any(!isempty,(L,T,Q,S,H,I,Y))
    P=if points===nothing
        explicit_topology ? Matrix{Int32}(undef,1,0) :
            reshape(Int32.(1:n),1,n)
    else
        _postview_connectivity(points,1,n,"point",remaining)
    end
    size(P,2)<=remaining || throw(ArgumentError(
        "PostViewField: element count exceeds max_elements=$element_limit"))
    ncells=sum(size(cells,2) for cells in (P,L,T,Q,S,H,I,Y))
    ncells>0 || throw(ArgumentError("PostViewField: view contains no elements"))
    _validate_postview_cells(C,L,T,Q,S,H,I,Y)

    used=falses(n)
    @inbounds for cells in (P,L,T,Q,S,H,I,Y), index in cells
        used[index]=true
    end
    nearest_nodes=Int32[i for i in 1:n if used[i]]
    nearest_coords=Matrix{Float64}(undef,3,length(nearest_nodes))
    @inbounds for (j,index) in enumerate(nearest_nodes),d in 1:3
        nearest_coords[d,j]=C[d,index]
    end
    tree=DistanceField(;points=nearest_coords)
    bvh=_build_postview_bvh(C,S,H,I,Y,T,Q,L,P)
    return PostViewField(C,vec(V),UInt8(size(V,1)),P,L,T,Q,S,H,I,Y,
                         crop_negative,use_closest,tolerance,tree,bvh...)
end

@inline function _postview_scaled_query(p,a,scale)
    delta=_sub3(p,a)
    all(isfinite,delta) && return (delta[1]/scale,delta[2]/scale,delta[3]/scale)
    return (p[1]/scale-a[1]/scale,p[2]/scale-a[2]/scale,
            p[3]/scale-a[3]/scale)
end

@inline _postview_component(field::PostViewField,component::Int,index::Integer)=
    field.values[component+Int(field.num_components)*(Int(index)-1)]

@inline function _postview_weighted_component(field,cells,j,weights,arity::Int,
                                              component::Int)
    value=0.0
    @inbounds for slot in 1:arity
        value=muladd(weights[slot],_postview_component(
            field,component,cells[slot,j]),value)
    end
    return value
end

@inline function _postview_weighted_scalar_operator(field,cells,j,weights,arity::Int)
    if field.num_components==1
        return _postview_weighted_component(field,cells,j,weights,arity,1)
    end
    x=_postview_weighted_component(field,cells,j,weights,arity,1)
    y=_postview_weighted_component(field,cells,j,weights,arity,2)
    z=_postview_weighted_component(field,cells,j,weights,arity,3)
    return hypot(x,y,z)
end

@inline function _postview_node_scalar_operator(field,index::Integer)
    field.num_components==1 && return _postview_component(field,1,index)
    return hypot(_postview_component(field,1,index),
                 _postview_component(field,2,index),
                 _postview_component(field,3,index))
end

function _postview_cell_value_big(field::PostViewField,kind::UInt8,j::Int,p,
                                  component::Int)::Tuple{Bool,Float64}
    cells=kind==3 ? field.tetrahedra : kind==2 ? field.triangles : field.lines
    arity=Int(kind)+1
    points=ntuple(slot->ntuple(k->_metric_rational(field.coords[k,cells[slot,j]]),3),arity)
    query=ntuple(k->_metric_rational(p[k]),3)
    tolerance=_metric_rational(field.reference_tolerance)
    if kind==1
        a,b=points;d=_sub3(b,a)
        midpoint=ntuple(k->(a[k]+b[k])/2,3);q=_sub3(query,midpoint)
        u=2*_dot3(q,d)/_dot3(d,d)
        artificial=if (abs(d[1])>=abs(d[2]) && abs(d[1])>=abs(d[3])) ||
                      (abs(d[2])>=abs(d[1]) && abs(d[2])>=abs(d[3]))
            (d[2],-d[1],zero(tolerance))
        else
            (zero(tolerance),d[3],-d[2])
        end
        normal=_postview_cross(d,artificial)
        v=_dot3(q,artificial)/_dot3(artificial,artificial)
        w=_dot3(q,normal)/_dot3(normal,normal)
        (u<-(1+tolerance)||u>1+tolerance||abs(v)>tolerance||abs(w)>tolerance) &&
            return false,0.0
        weights=((1-u)/2,(1+u)/2)
    elseif kind==2
        a,b,c=points;e1=_sub3(b,a);e2=_sub3(c,a);q=_sub3(query,a)
        jxy=e1[1]*e2[2]-e1[2]*e2[1]
        jxz=e1[1]*e2[3]-e1[3]*e2[1]
        jyz=e1[2]*e2[3]-e1[3]*e2[2]
        if abs(jxy)>abs(jxz) && abs(jxy)>abs(jyz)
            u=(q[1]*e2[2]-q[2]*e2[1])/jxy
            v=(q[2]*e1[1]-q[1]*e1[2])/jxy
        elseif abs(jxz)>abs(jyz)
            u=(q[1]*e2[3]-q[3]*e2[1])/jxz
            v=(q[3]*e1[1]-q[1]*e1[3])/jxz
        else
            u=(q[2]*e2[3]-q[3]*e2[2])/jyz
            v=(q[3]*e1[2]-q[2]*e1[3])/jyz
        end
        (u < -tolerance || v < -tolerance || u > (1+tolerance)-v) &&
            return false,0.0
        weights=(1-u-v,u,v)
    else
        a,b,c,d=points;e1=_sub3(b,a);e2=_sub3(c,a);e3=_sub3(d,a)
        q=_sub3(query,a);det=_postview_det(e1,e2,e3)
        u=_postview_det(q,e2,e3)/det
        v=_postview_det(e1,q,e3)/det
        w=_postview_det(e1,e2,q)/det
        (u < -tolerance || v < -tolerance || w < -tolerance ||
         u > (1+tolerance)-v-w) && return false,0.0
        weights=(1-u-v-w,u,v,w)
    end
    value=zero(tolerance)
    @inbounds for slot in 1:arity
        value+=weights[slot]*_metric_rational(
            _postview_component(field,component,cells[slot,j]))
    end
    return true,Float64(value)
end

@inline function _postview_big_value_or_nan(field,kind,j,p)::Float64
    found,x=_postview_cell_value_big(field,kind,j,p,1)
    found || return NaN
    field.num_components==1 && return x
    _,y=_postview_cell_value_big(field,kind,j,p,2)
    _,z=_postview_cell_value_big(field,kind,j,p,3)
    return hypot(x,y,z)
end

@inline function _postview_line_value(field::PostViewField,j::Int,p)::Float64
    cells=field.lines
    a=_postview_point(field.coords,cells[1,j])
    b=_postview_point(field.coords,cells[2,j])
    d=_sub3(b,a)
    if all(isfinite,d)
        scale=max(abs(d[1]),abs(d[2]),abs(d[3]))
        scale>0 && (d=(d[1]/scale,d[2]/scale,d[3]/scale))
    else
        scale=max(abs(a[1]),abs(a[2]),abs(a[3]),
                  abs(b[1]),abs(b[2]),abs(b[3]))
        if scale>0
            abase=(a[1]/scale,a[2]/scale,a[3]/scale)
            d=(b[1]/scale-abase[1],b[2]/scale-abase[2],b[3]/scale-abase[3])
        end
    end
    scale>0 || return _postview_big_value_or_nan(field,UInt8(1),j,p)
    q=_postview_scaled_query(p,a,scale)
    all(isfinite,q) || return _postview_big_value_or_nan(field,UInt8(1),j,p)
    midpoint_offset=(q[1]-d[1]/2,q[2]-d[2]/2,q[3]-d[3]/2)
    denominator=_dot3(d,d)
    u=2*_dot3(midpoint_offset,d)/denominator
    artificial=if (abs(d[1])>=abs(d[2]) && abs(d[1])>=abs(d[3])) ||
                  (abs(d[2])>=abs(d[1]) && abs(d[2])>=abs(d[3]))
        (d[2],-d[1],0.0)
    else
        (0.0,d[3],-d[2])
    end
    normal=_postview_cross(d,artificial)
    v=_dot3(midpoint_offset,artificial)/_dot3(artificial,artificial)
    # Gmsh forms the second artificial vector from two physical-length
    # vectors. Its coefficient therefore has one extra inverse-length factor.
    w=_dot3(midpoint_offset,normal)/(scale*_dot3(normal,normal))
    all(isfinite,(u,v,w)) || return _postview_big_value_or_nan(field,UInt8(1),j,p)
    tolerance=field.reference_tolerance
    guard=64eps(Float64)*(1+abs(u)+abs(v)+abs(w)+tolerance)
    min(abs(u+1+tolerance),abs(1+tolerance-u),
        abs(abs(v)-tolerance),abs(abs(w)-tolerance))<=guard &&
        return _postview_big_value_or_nan(field,UInt8(1),j,p)
    (u<-(1+tolerance)||u>1+tolerance||abs(v)>tolerance||abs(w)>tolerance) &&
        return NaN
    value=_postview_weighted_scalar_operator(
        field,cells,j,((1-u)/2,(1+u)/2),2)
    isnan(value) && throw(ArgumentError(
        "PostViewField: line $j produced a NaN interpolant"))
    return value
end

@inline function _postview_triangle_value(field::PostViewField,j::Int,p)::Float64
    cells=field.triangles
    a=_postview_point(field.coords,cells[1,j])
    b=_postview_point(field.coords,cells[2,j])
    c=_postview_point(field.coords,cells[3,j])
    e1=_sub3(b,a);e2=_sub3(c,a)
    if all(isfinite,e1) && all(isfinite,e2)
        scale=max(abs(e1[1]),abs(e1[2]),abs(e1[3]),
                  abs(e2[1]),abs(e2[2]),abs(e2[3]))
        if scale>0
            e1=(e1[1]/scale,e1[2]/scale,e1[3]/scale)
            e2=(e2[1]/scale,e2[2]/scale,e2[3]/scale)
        end
    else
        scale=max(abs(a[1]),abs(a[2]),abs(a[3]),
                  abs(b[1]),abs(b[2]),abs(b[3]),
                  abs(c[1]),abs(c[2]),abs(c[3]))
        if scale>0
            abase=(a[1]/scale,a[2]/scale,a[3]/scale)
            e1=(b[1]/scale-abase[1],b[2]/scale-abase[2],b[3]/scale-abase[3])
            e2=(c[1]/scale-abase[1],c[2]/scale-abase[2],c[3]/scale-abase[3])
        end
    end
    scale>0 || return _postview_big_value_or_nan(field,UInt8(2),j,p)
    q=_postview_scaled_query(p,a,scale)
    all(isfinite,q) || return _postview_big_value_or_nan(field,UInt8(2),j,p)
    jxy=e1[1]*e2[2]-e1[2]*e2[1]
    jxz=e1[1]*e2[3]-e1[3]*e2[1]
    jyz=e1[2]*e2[3]-e1[3]*e2[2]
    if abs(jxy)>abs(jxz) && abs(jxy)>abs(jyz)
        abs(jxy)<=16eps(Float64)*(abs(e1[1]*e2[2])+abs(e1[2]*e2[1])) &&
            return _postview_big_value_or_nan(field,UInt8(2),j,p)
        u=(q[1]*e2[2]-q[2]*e2[1])/jxy
        v=(q[2]*e1[1]-q[1]*e1[2])/jxy
    elseif abs(jxz)>abs(jyz)
        abs(jxz)<=16eps(Float64)*(abs(e1[1]*e2[3])+abs(e1[3]*e2[1])) &&
            return _postview_big_value_or_nan(field,UInt8(2),j,p)
        u=(q[1]*e2[3]-q[3]*e2[1])/jxz
        v=(q[3]*e1[1]-q[1]*e1[3])/jxz
    elseif !iszero(jyz)
        abs(jyz)<=16eps(Float64)*(abs(e1[2]*e2[3])+abs(e1[3]*e2[2])) &&
            return _postview_big_value_or_nan(field,UInt8(2),j,p)
        u=(q[2]*e2[3]-q[3]*e2[2])/jyz
        v=(q[3]*e1[2]-q[2]*e1[3])/jyz
    else
        return _postview_big_value_or_nan(field,UInt8(2),j,p)
    end
    all(isfinite,(u,v)) || return _postview_big_value_or_nan(field,UInt8(2),j,p)
    tolerance=field.reference_tolerance
    guard=64eps(Float64)*(1+abs(u)+abs(v)+tolerance)
    min(abs(u+tolerance),abs(v+tolerance),
        abs(1+tolerance-u-v))<=guard &&
        return _postview_big_value_or_nan(field,UInt8(2),j,p)
    (u < -tolerance || v < -tolerance || u > (1+tolerance)-v) && return NaN
    value=_postview_weighted_scalar_operator(field,cells,j,(1-u-v,u,v),3)
    isnan(value) && throw(ArgumentError(
        "PostViewField: triangle $j produced a NaN interpolant"))
    return value
end

@inline function _postview_tetrahedron_value(field::PostViewField,j::Int,p)::Float64
    cells=field.tetrahedra
    a=_postview_point(field.coords,cells[1,j])
    b=_postview_point(field.coords,cells[2,j])
    c=_postview_point(field.coords,cells[3,j])
    d=_postview_point(field.coords,cells[4,j])
    e1=_sub3(b,a);e2=_sub3(c,a);e3=_sub3(d,a)
    if all(isfinite,e1) && all(isfinite,e2) && all(isfinite,e3)
        scale=max(abs(e1[1]),abs(e1[2]),abs(e1[3]),
                  abs(e2[1]),abs(e2[2]),abs(e2[3]),
                  abs(e3[1]),abs(e3[2]),abs(e3[3]))
        if scale>0
            e1=(e1[1]/scale,e1[2]/scale,e1[3]/scale)
            e2=(e2[1]/scale,e2[2]/scale,e2[3]/scale)
            e3=(e3[1]/scale,e3[2]/scale,e3[3]/scale)
        end
    else
        scale=max(abs(a[1]),abs(a[2]),abs(a[3]),
                  abs(b[1]),abs(b[2]),abs(b[3]),
                  abs(c[1]),abs(c[2]),abs(c[3]),
                  abs(d[1]),abs(d[2]),abs(d[3]))
        if scale>0
            abase=(a[1]/scale,a[2]/scale,a[3]/scale)
            e1=(b[1]/scale-abase[1],b[2]/scale-abase[2],b[3]/scale-abase[3])
            e2=(c[1]/scale-abase[1],c[2]/scale-abase[2],c[3]/scale-abase[3])
            e3=(d[1]/scale-abase[1],d[2]/scale-abase[2],d[3]/scale-abase[3])
        end
    end
    scale>0 || return _postview_big_value_or_nan(field,UInt8(3),j,p)
    q=_postview_scaled_query(p,a,scale)
    all(isfinite,q) || return _postview_big_value_or_nan(field,UInt8(3),j,p)
    p1=e1[1]*e2[2]*e3[3];p2=e1[2]*e2[3]*e3[1]
    p3=e1[3]*e2[1]*e3[2];p4=e1[3]*e2[2]*e3[1]
    p5=e1[2]*e2[1]*e3[3];p6=e1[1]*e2[3]*e3[2]
    determinant=_postview_det(e1,e2,e3)
    iszero(determinant) && return _postview_big_value_or_nan(field,UInt8(3),j,p)
    abs(determinant)<=32eps(Float64)*
        (abs(p1)+abs(p2)+abs(p3)+abs(p4)+abs(p5)+abs(p6)) &&
        return _postview_big_value_or_nan(field,UInt8(3),j,p)
    u=_postview_det(q,e2,e3)/determinant
    v=_postview_det(e1,q,e3)/determinant
    w=_postview_det(e1,e2,q)/determinant
    all(isfinite,(u,v,w)) || return _postview_big_value_or_nan(field,UInt8(3),j,p)
    tolerance=field.reference_tolerance
    guard=64eps(Float64)*(1+abs(u)+abs(v)+abs(w)+tolerance)
    min(abs(u+tolerance),abs(v+tolerance),abs(w+tolerance),
        abs(1+tolerance-u-v-w))<=guard &&
        return _postview_big_value_or_nan(field,UInt8(3),j,p)
    (u < -tolerance || v < -tolerance || w < -tolerance ||
     u > (1+tolerance)-v-w) && return NaN
    value=_postview_weighted_scalar_operator(field,cells,j,(1-u-v-w,u,v,w),4)
    isnan(value) && throw(ArgumentError(
        "PostViewField: tetrahedron $j produced a NaN interpolant"))
    return value
end

@inline function _postview_nonlinear_cells(field::PostViewField,kind::UInt8)
    kind==4 && return field.hexahedra,8,"hexahedron"
    kind==5 && return field.prisms,6,"prism"
    kind==6 && return field.pyramids,5,"pyramid"
    kind==7 && return field.quadrangles,4,"quadrangle"
    throw(ErrorException("PostViewField: internal unsupported nonlinear cell kind $kind"))
end

@inline function _postview_nonlinear_scalar_operator(field,cells,j,arity::Int,
                                                     kind::UInt8,u,v,w)
    x=0.0;y=0.0;z=0.0
    if field.num_components==1
        @inbounds for slot in 1:arity
            weight=_postview_shape(kind,slot,u,v,w)
            x=muladd(weight,_postview_component(field,1,cells[slot,j]),x)
        end
        return x
    end
    @inbounds for slot in 1:arity
        weight=_postview_shape(kind,slot,u,v,w);index=cells[slot,j]
        x=muladd(weight,_postview_component(field,1,index),x)
        y=muladd(weight,_postview_component(field,2,index),y)
        z=muladd(weight,_postview_component(field,3,index),z)
    end
    return hypot(x,y,z)
end

@inline function _postview_reference_inside(kind::UInt8,u,v,w,tolerance)
    if kind==7
        return u>=-(1+tolerance) && u<=1+tolerance &&
               v>=-(1+tolerance) && v<=1+tolerance && abs(w)<=tolerance
    elseif kind==4
        return u>=-(1+tolerance) && u<=1+tolerance &&
               v>=-(1+tolerance) && v<=1+tolerance &&
               w>=-(1+tolerance) && w<=1+tolerance
    elseif kind==5
        return w>=-(1+tolerance) && w<=1+tolerance && u>=-tolerance &&
               v>=-tolerance && u<=(1+tolerance)-v
    end
    return u>=w-(1+tolerance) && u<=(1+tolerance)-w &&
           v>=w-(1+tolerance) && v<=(1+tolerance)-w &&
           w>=-tolerance && w<=1+tolerance
end

@inline function _postview_nonlinear_value(field::PostViewField,kind::UInt8,
                                           j::Int,p)::Float64
    cells,arity,what=_postview_nonlinear_cells(field,kind)
    scale=_postview_cell_scale(field.coords,cells,j,arity)
    scale>0 || return NaN
    base=cells[1,j];a=_postview_point(field.coords,base)
    query=_postview_scaled_query(p,a,scale)
    all(isfinite,query) || return NaN
    u=0.0;v=0.0;w=0.0;error=1.0;iteration=1
    while error>1e-6 && iteration<20
        mapped=(0.0,0.0,0.0)
        @inbounds for slot in 1:arity
            weight=_postview_shape(kind,slot,u,v,w);index=cells[slot,j]
            mapped=(muladd(weight,_postview_relative(
                        field.coords,index,base,1,scale),mapped[1]),
                    muladd(weight,_postview_relative(
                        field.coords,index,base,2,scale),mapped[2]),
                    muladd(weight,_postview_relative(
                        field.coords,index,base,3,scale),mapped[3]))
        end
        residual=_sub3(query,mapped)
        du,dv,dw=_postview_reference_jacobian(
            field.coords,cells,j,arity,kind,u,v,w,scale)
        determinant=_postview_det(du,dv,dw)
        (isfinite(determinant) && !iszero(determinant)) || return NaN
        delta_u=_postview_det(residual,dv,dw)/determinant
        delta_v=_postview_det(du,residual,dw)/determinant
        delta_w=_postview_det(du,dv,residual)/determinant
        all(isfinite,(delta_u,delta_v,delta_w)) || return NaN
        un=u+delta_u;vn=v+delta_v;wn=w+delta_w
        all(isfinite,(un,vn,wn)) || return NaN
        error=hypot(un-u,vn-v,wn-w)
        u=un;v=vn;w=wn;iteration+=1
    end
    _postview_reference_inside(kind,u,v,w,field.reference_tolerance) || return NaN
    value=_postview_nonlinear_scalar_operator(field,cells,j,arity,kind,u,v,w)
    isnan(value) && throw(ArgumentError(
        "PostViewField: $what $j produced a NaN interpolant"))
    return value
end

@inline function _postview_cell_value(field::PostViewField,id::Int,p)::Float64
    kind,j=_postview_cell_kind(field,id)
    if kind==0
        index=field.points[1,j]
        found=p[1]==field.coords[1,index] && p[2]==field.coords[2,index] &&
              p[3]==field.coords[3,index]
        return found ? _postview_node_scalar_operator(field,index) : NaN
    end
    kind==1 && return _postview_line_value(field,j,p)
    kind==2 && return _postview_triangle_value(field,j,p)
    kind==3 && return _postview_tetrahedron_value(field,j,p)
    kind in UInt8(4):UInt8(7) && return _postview_nonlinear_value(field,kind,j,p)
    throw(ErrorException("PostViewField: internal unsupported cell kind $kind"))
end

@inline function _postview_in_box(p,lo,hi)
    return lo[1]<=p[1]<=hi[1] && lo[2]<=p[2]<=hi[2] && lo[3]<=p[3]<=hi[3]
end

function _postview_bvh_value(field::PostViewField,node::Int,p,best_id::Int,
                             best_value::Float64,found::Bool)::Tuple{Int,Float64,Bool}
    _postview_in_box(p,field.bvh_lo[node],field.bvh_hi[node]) ||
        return best_id,best_value,found
    count=field.bvh_count[node]
    if count>0
        first=field.bvh_first[node]
        @inbounds for pos in first:first+count-1
            id=field.bvh_order[pos]
            id>=best_id && continue
            value=_postview_cell_value(field,id,p)
            if !isnan(value)
                best_id=id;best_value=value;found=true
            end
        end
        return best_id,best_value,found
    end
    best_id,best_value,found=_postview_bvh_value(
        field,field.bvh_left[node],p,best_id,best_value,found)
    return _postview_bvh_value(
        field,field.bvh_right[node],p,best_id,best_value,found)
end

function field_value(field::PostViewField,x,y,z)
    p=(_float_value(x,"PostViewField","x"),_float_value(y,"PostViewField","y"),
       _float_value(z,"PostViewField","z"))
    # `PostViewField::operator()`, unlike its metric overload, deliberately has
    # no scalar interpretation for tensor views in Gmsh 4.15.2.
    field.num_components==9 && return GMSH_MAX_SIZE
    _,v,found=_postview_bvh_value(field,1,p,typemax(Int),0.0,false)
    if !found
        field.use_closest || return GMSH_MAX_SIZE
        _,nearest=_nearest_point(field.tree,p)
        closest=(field.tree.points[1,nearest],field.tree.points[2,nearest],
                 field.tree.points[3,nearest])
        _,v,found=_postview_bvh_value(field,1,closest,typemax(Int),0.0,false)
        found || throw(ErrorException(
            "PostViewField: closest active node is not contained in any view element"))
    end
    (v<=0 && field.crop_negative) && return GMSH_MAX_SIZE
    return _checked_field_result(v,"PostViewField",p...)
end

# ── Anisotropic composition and attractor ─────────────────────────────────────

# Gmsh's anisotropic list operators use a default-constructed identity metric
# for unresolved/self references in their scalar overloads, while their metric
# overloads skip some entries. Preserve list position explicitly.
struct _GmshAnisoReferencePlaceholder <: AbstractAnisoField end
@inline metric_at(::_GmshAnisoReferencePlaceholder,x,y,z)=
    Metric3(1.0,1.0,1.0,0.0,0.0,0.0)
@inline field_value(::_GmshAnisoReferencePlaceholder,x,y,z)=1.0

"""Anisotropic minimum: intersect every child metric in the Loewner order."""
struct MinAnisoField{F<:Tuple} <: AbstractAnisoField
    fields::F
end
function MinAnisoField(fields::Union{Tuple,AbstractVector})
    values=Tuple(fields)
    all(f->f isa AbstractField,values) ||
        throw(ArgumentError("MinAnisoField: every input must be an AbstractField"))
    return MinAnisoField{typeof(values)}(values)
end
function _minaniso_metric(field::MinAnisoField,x,y,z,entity)
    # Gmsh seeds MinAniso with SMetric3(1 / MAX_LC), corresponding to
    # directional size sqrt(MAX_LC).
    seed=inv(GMSH_MAX_SIZE)
    m=Metric3(seed,seed,seed,0.0,0.0,0.0)
    for child in field.fields
        child isa _GmshAnisoReferencePlaceholder && continue
        m=intersection_conserve_mostaniso(m,_as_metric(child,x,y,z,entity))
    end
    return m
end
metric_at(field::MinAnisoField,x,y,z)=_minaniso_metric(field,x,y,z,nothing)
metric_at(field::MinAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_minaniso_metric(field,x,y,z,entity)

# Gmsh's scalar MinAniso operator uses its ordinary `intersection` fold, which
# is intentionally distinct from the most-anisotropic conservative metric fold
# above. This is the value observed by MathEval F-references.
function _minaniso_value(field::MinAnisoField,x,y,z,entity)
    seed=inv(GMSH_MAX_SIZE)
    metric=Metric3(seed,seed,seed,0.0,0.0,0.0)
    for child in field.fields
        metric=_gmsh_metric_intersection(metric,_as_metric(child,x,y,z,entity))
    end
    return metric_size(metric)
end
field_value(field::MinAnisoField,x,y,z)=_minaniso_value(field,x,y,z,nothing)
field_value(field::MinAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_minaniso_value(field,x,y,z,entity)

"""Pairwise conservative anisotropic metric intersection."""
struct IntersectAnisoField{F<:Tuple} <: AbstractAnisoField
    fields::F
end
function IntersectAnisoField(fields::Union{Tuple,AbstractVector})
    values=Tuple(fields)
    all(f->f isa AbstractField,values) ||
        throw(ArgumentError("IntersectAnisoField: every input must be an AbstractField"))
    return IntersectAnisoField{typeof(values)}(values)
end
function _as_metric(child,x,y,z,entity)
    child isa AbstractAnisoField && return metric_at(child,x,y,z,entity)
    # Gmsh's anisotropic composers square raw scalar values. Intermediate
    # scalar fields may therefore be negative; only zero is unrepresentable.
    h=_checked_field_result(field_value(child,x,y,z,entity),
                            "anisotropic field composition",x,y,z)
    h==0 && throw(ArgumentError(
        "anisotropic field composition: zero scalar value at ($x,$y,$z)"))
    return isotropic_metric(abs(h))
end
function _intersect_aniso_metric(field::IntersectAnisoField,x,y,z,entity)
    isempty(field.fields) && return Metric3(1.0,1.0,1.0,0.0,0.0,0.0)
    value=Metric3(1.0,1.0,1.0,0.0,0.0,0.0)
    @inbounds for i in eachindex(field.fields)
        child=field.fields[i]
        child isa _GmshAnisoReferencePlaceholder && continue
        metric=_as_metric(child,x,y,z,entity)
        value=i==firstindex(field.fields) ? metric :
              _gmsh_metric_intersection_alauzet(value,metric)
    end
    return value
end
metric_at(field::IntersectAnisoField,x,y,z)=
    _intersect_aniso_metric(field,x,y,z,nothing)
metric_at(field::IntersectAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_intersect_aniso_metric(field,x,y,z,entity)

function _intersect_aniso_value(field::IntersectAnisoField,x,y,z,entity)
    metric=Metric3(1.0,1.0,1.0,0.0,0.0,0.0)
    @inbounds for i in eachindex(field.fields)
        child_metric=_as_metric(field.fields[i],x,y,z,entity)
        metric=i==firstindex(field.fields) ? child_metric :
               _gmsh_metric_intersection_alauzet(metric,child_metric)
    end
    return metric_size(metric)
end
field_value(field::IntersectAnisoField,x,y,z)=
    _intersect_aniso_value(field,x,y,z,nothing)
field_value(field::IntersectAnisoField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_intersect_aniso_value(field,x,y,z,entity)

"""
    AttractorAnisoCurveField(points, tangents; dist_min, dist_max,
                             size_min_tangent, size_max_tangent,
                             size_min_normal, size_max_normal)

Gmsh `AttractorAnisoCurve` on an already-sampled curve.
"""
struct AttractorAnisoCurveField <: AbstractAnisoField
    points::Matrix{Float64}
    tangents::Matrix{Float64}
    dist_min::Float64; dist_max::Float64
    lmin_t::Float64; lmax_t::Float64; lmin_n::Float64; lmax_n::Float64
    tree::DistanceField
end
function AttractorAnisoCurveField(points, tangents; dist_min::Real, dist_max::Real,
                                  size_min_tangent::Real, size_max_tangent::Real,
                                  size_min_normal::Real, size_max_normal::Real)
    P=_point_matrix(points,"AttractorAnisoCurveField")
    T=_point_matrix(tangents,"AttractorAnisoCurveField")
    size(P,2)==size(T,2) && size(P,2)>0 ||
        throw(ArgumentError("AttractorAnisoCurveField: need matching non-empty samples"))
    dmin=_float_value(dist_min,"AttractorAnisoCurveField","dist_min")
    dmax=_float_value(dist_max,"AttractorAnisoCurveField","dist_max")
    dmax>dmin || throw(ArgumentError("AttractorAnisoCurveField: require dist_max > dist_min"))
    @inbounds for j in axes(T,2)
        _norm3((T[1,j],T[2,j],T[3,j]))>0 || throw(ArgumentError(
            "AttractorAnisoCurveField: zero tangent at sample $j"))
    end
    tree=DistanceField(;points=[(P[1,j],P[2,j],P[3,j]) for j in axes(P,2)])
    return AttractorAnisoCurveField(P,T,dmin,dmax,
        _positive_value(size_min_tangent,"AttractorAnisoCurveField","size_min_tangent"),
        _positive_value(size_max_tangent,"AttractorAnisoCurveField","size_max_tangent"),
        _positive_value(size_min_normal,"AttractorAnisoCurveField","size_min_normal"),
        _positive_value(size_max_normal,"AttractorAnisoCurveField","size_max_normal"),tree)
end
function _lerp_size(d,dmin,dmax,amin,amax)
    d<=dmin && return amin
    d>=dmax && return amax
    return amin+(amax-amin)*(d-dmin)/(dmax-dmin)
end
function metric_at(field::AttractorAnisoCurveField,x,y,z)
    p=(_float_value(x,"AttractorAnisoCurveField","x"),
       _float_value(y,"AttractorAnisoCurveField","y"),
       _float_value(z,"AttractorAnisoCurveField","z"))
    best,ibest=_nearest_point(field.tree,p)
    lt=_lerp_size(best,field.dist_min,field.dist_max,field.lmin_t,field.lmax_t)
    ln=_lerp_size(best,field.dist_min,field.dist_max,field.lmin_n,field.lmax_n)
    tangent=(field.tangents[1,ibest],field.tangents[2,ibest],field.tangents[3,ibest])
    nt=_norm3(tangent); nt>0 || throw(ArgumentError("AttractorAnisoCurveField: zero tangent"))
    unit_tangent=(tangent[1]/nt,tangent[2]/nt,tangent[3]/nt)
    helper=abs(unit_tangent[1])>abs(unit_tangent[2]) ?
           (0.0,1.0,0.0) : (1.0,0.0,0.0)
    n0=_cross3(unit_tangent,helper)
    if !(_norm3(n0)>0)
        helper=(0.0,0.0,1.0)
        n0=_cross3(unit_tangent,helper)
    end
    n1=_cross3(unit_tangent,n0)
    # Match Gmsh 4.15.2: the tangent is normalized during curve sampling, but
    # the two cross-product normals are passed to SMetric3 without renormalizing.
    return _metric_from_raw_axes(1/(lt*lt),1/(ln*ln),1/(ln*ln),unit_tangent,n0,n1)
end
@inline _cross3(a,b)=(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])

# Gmsh exposes a separate scalar operator for F-references: the closest sampled
# curve distance, cropped from below at 0.05, independent of the metric sizes.
function field_value(field::AttractorAnisoCurveField,x,y,z)
    p=(_float_value(x,"AttractorAnisoCurveField","x"),
       _float_value(y,"AttractorAnisoCurveField","y"),
       _float_value(z,"AttractorAnisoCurveField","z"))
    distance,_=_nearest_point(field.tree,p)
    return max(distance,0.05)
end
field_value(field::AttractorAnisoCurveField,x,y,z,
            ::Tuple{T,U}) where {T<:Integer,U<:Integer}=field_value(field,x,y,z)

# ── Boundary layer size ───────────────────────────────────────────────────────

"""
    BoundaryLayerField(distance; hwall, ratio, thickness, hfar)

Gmsh `BoundaryLayer` scalar law: `h = min(hFar, dist*(ratio-1) + hWall)`
inside `thickness*ratio`, otherwise `GMSH_MAX_SIZE`. Options are retained as
finite raw scalars: non-positive combinations can be no-ops or intermediate
non-positive values, with positivity enforced only when the field is consumed
as a final mesh size.
"""
struct BoundaryLayerField{F<:AbstractField} <: AbstractSizeField
    distance::F
    hwall::Float64
    ratio::Float64
    thickness::Float64
    hfar::Float64
end
function BoundaryLayerField(distance::AbstractField; hwall::Real, ratio::Real,
                            thickness::Real, hfar::Real)
    hw=_float_value(hwall,"BoundaryLayerField","hwall")
    ra=_float_value(ratio,"BoundaryLayerField","ratio")
    th=_float_value(thickness,"BoundaryLayerField","thickness")
    hf=_float_value(hfar,"BoundaryLayerField","hfar")
    return BoundaryLayerField{typeof(distance)}(distance,hw,ra,th,hf)
end
function _boundary_layer_value(field::BoundaryLayerField,x,y,z,entity)
    d=_checked_field_result(field_value(field.distance,x,y,z,entity),
                            "BoundaryLayerField",x,y,z)
    d<0 && throw(ArgumentError("BoundaryLayerField: distance is negative at ($x,$y,$z)"))
    d>field.thickness*field.ratio && return GMSH_MAX_SIZE
    value=min(field.hfar, d*(field.ratio-1)+field.hwall)
    return _checked_field_result(value,"BoundaryLayerField",x,y,z)
end
field_value(field::BoundaryLayerField,x,y,z)=
    _boundary_layer_value(field,x,y,z,nothing)
field_value(field::BoundaryLayerField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_boundary_layer_value(field,x,y,z,entity)

# ── Automatic mesh size from discrete curvature ───────────────────────────────

"""
    AutomaticMeshSizeField(surface; n_nodes_per_circle=20, hmin=nothing, hmax=nothing)

Native analogue of Gmsh `AutomaticMeshSizeField` on a triangulated surface:
per-vertex sphere-fit curvature `κ` yields `h = clamp(2π /(n κ), hmin, hmax)`,
evaluated by closest-vertex interpolation.
"""
struct AutomaticMeshSizeField <: AbstractSizeField
    coords::Matrix{Float64}
    h::Vector{Float64}
    tree::DistanceField
end
function AutomaticMeshSizeField(surface::Mesh; n_nodes_per_circle::Real=20.0,
                                hmin=nothing, hmax=nothing)
    n=_positive_value(n_nodes_per_circle,"AutomaticMeshSizeField","n_nodes_per_circle")
    diagnostic=validate(surface;require_positive_tets=false,require_manifold_tris=true)
    diagnostic.ok || throw(ArgumentError(
        "AutomaticMeshSizeField: invalid surface mesh: $(join(diagnostic.messages, "; "))"))
    size(surface.tris,2)>0 || throw(ArgumentError("AutomaticMeshSizeField: surface has no triangles"))
    xmin,xmax=extrema(@view surface.coords[1,:])
    ymin,ymax=extrema(@view surface.coords[2,:])
    zmin,zmax=extrema(@view surface.coords[3,:])
    # Gmsh's automatic-size backend defines its model scale as the largest
    # bounding-box span (distinct from CTX::lc's diagonal convention).
    model_length=max(xmax-xmin,ymax-ymin,zmax-zmin)
    (isfinite(model_length) && model_length>0) || throw(ArgumentError(
        "AutomaticMeshSizeField: maximum surface bounding-box span must be positive and finite"))
    hn=hmin===nothing ? model_length/1000 :
       _positive_value(hmin,"AutomaticMeshSizeField","hmin")
    hx=hmax===nothing ? model_length/20 :
       _positive_value(hmax,"AutomaticMeshSizeField","hmax")
    (isfinite(hn) && hn>0 && isfinite(hx) && hx>0) || throw(ArgumentError(
        "AutomaticMeshSizeField: default sizes are not representable"))
    hx>=hn || throw(ArgumentError("AutomaticMeshSizeField: require hmax >= hmin"))
    nv=size(surface.coords,2)
    rings=[Int[] for _ in 1:nv]
    @inbounds for t in axes(surface.tris,2)
        a=Int(surface.tris[1,t]); b=Int(surface.tris[2,t]); c=Int(surface.tris[3,t])
        push!(rings[a],b,c); push!(rings[b],a,c); push!(rings[c],a,b)
    end
    h=Vector{Float64}(undef,nv)
    @inbounds for i in 1:nv
        nbr=unique(rings[i])
        κ=_sphere_fit_curvature(surface.coords,i,nbr)
        raw=κ>0 ? 2*π/(n*κ) : hx
        h[i]=clamp(raw,hn,hx)
        isfinite(h[i]) && h[i]>0 || throw(ArgumentError(
            "AutomaticMeshSizeField: non-positive size at vertex $i"))
    end
    coords=copy(surface.coords)
    tree=DistanceField(;points=[(coords[1,j],coords[2,j],coords[3,j]) for j in axes(coords,2)])
    return AutomaticMeshSizeField(coords,h,tree)
end
function _sphere_fit_curvature(coords,i,nbr)
    isempty(nbr) && return 0.0
    # Algebraic fit of xᵀx + d·x + e = 0 through vertex i and its 1-ring.
    pts=Vector{NTuple{3,Float64}}(undef,length(nbr)+1)
    pts[1]=(coords[1,i],coords[2,i],coords[3,i])
    @inbounds for (k,j) in enumerate(nbr)
        pts[k+1]=(coords[1,j],coords[2,j],coords[3,j])
    end
    # Shift to vertex i for conditioning.
    c0=pts[1]
    scale=0.0
    @inbounds for p in pts
        scale=max(scale,abs(p[1]-c0[1]),abs(p[2]-c0[2]),abs(p[3]-c0[3]))
    end
    (isfinite(scale) && scale>0) || return 0.0
    N=zeros(Float64,4,4); b=zeros(Float64,4)
    @inbounds for k in eachindex(pts)
        x=(pts[k][1]-c0[1])/scale
        y=(pts[k][2]-c0[2])/scale
        z=(pts[k][3]-c0[3])/scale
        row=(x,y,z,1.0); rhs=-(x*x+y*y+z*z)
        for i in 1:4, j in 1:4
            N[i,j]+=row[i]*row[j]
        end
        for i in 1:4; b[i]+=row[i]*rhs; end
    end
    d=_solve4(N,b)
    d===nothing && return 0.0
    center=(-0.5d[1],-0.5d[2],-0.5d[3])
    r2=center[1]^2+center[2]^2+center[3]^2-d[4]
    r2>0 || return 0.0
    radius=scale*sqrt(r2)
    radius>0 || return Inf
    return inv(radius)
end
function _solve4(A::Matrix{Float64},b::Vector{Float64})
    M=copy(A); x=copy(b)
    @inbounds for k in 1:4
        piv=k; amax=abs(M[k,k])
        for i in k+1:4
            if abs(M[i,k])>amax; amax=abs(M[i,k]); piv=i; end
        end
        amax<=eps(Float64) && return nothing
        if piv!=k
            for j in k:4; M[k,j],M[piv,j]=M[piv,j],M[k,j]; end
            x[k],x[piv]=x[piv],x[k]
        end
        for i in k+1:4
            f=M[i,k]/M[k,k]; M[i,k]=0
            for j in k+1:4; M[i,j]-=f*M[k,j]; end
            x[i]-=f*x[k]
        end
    end
    @inbounds for i in 4:-1:1
        s=x[i]
        for j in i+1:4; s-=M[i,j]*x[j]; end
        abs(M[i,i])<=eps(Float64) && return nothing
        x[i]=s/M[i,i]
    end
    return x
end
function field_value(field::AutomaticMeshSizeField,x,y,z)
    p=(_float_value(x,"AutomaticMeshSizeField","x"),
       _float_value(y,"AutomaticMeshSizeField","y"),
       _float_value(z,"AutomaticMeshSizeField","z"))
    _,ibest=_nearest_point(field.tree,p)
    return field.h[ibest]
end

# ── External process ──────────────────────────────────────────────────────────

"""
    ExternalProcessField(command)

Gmsh `ExternalProcess` protocol: write three native `Float64` coordinates and
read one `Float64` size. An empty command or a failed handshake is a blocker.
"""
mutable struct ExternalProcessField <: AbstractSizeField
    command::String
    io::Union{Nothing,Base.Process}
    lock::ReentrantLock
    closed::Bool
end
function ExternalProcessField(command::AbstractString)
    cmd=String(strip(command))
    isempty(cmd) && throw(ArgumentError("ExternalProcessField: command must be non-empty"))
    field=ExternalProcessField(cmd,nothing,ReentrantLock(),false)
    finalizer(field) do value
        _ext_close_io!(value,false)
    end
    return field
end
function _ext_start!(field::ExternalProcessField)
    field.closed && throw(ArgumentError("ExternalProcessField: field is closed"))
    field.io===nothing || return
    try
        field.io=open(Cmd(["/bin/sh","-c",field.command]); read=true, write=true)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("ExternalProcessField: failed to launch $(field.command): $err"))
    end
end

function _ext_close_io!(field::ExternalProcessField,throw_errors::Bool)
    io=field.io
    field.io=nothing
    io===nothing && return nothing
    first_error=nothing
    # Some failed or partially initialized process streams throw from `isopen`.
    # Attempting the protocol terminator directly keeps finalization best-effort
    # while explicit `close` still reports a real transport error.
    try
        write(io,NaN); write(io,NaN); write(io,NaN); flush(io)
    catch err
        err isa InterruptException && throw_errors && rethrow()
        # Once the NaN shutdown record reaches the child it may close its read
        # end before Julia finishes flushing/closing the bidirectional process.
        # EPIPE is therefore an expected shutdown race; `close(io)` below still
        # reports an abnormal child exit through its own non-EPIPE exception.
        !(err isa Base.IOError && err.code==Base.UV_EPIPE) && (first_error=err)
    end
    # Close only the child's stdin before waiting. Closing the read side here
    # can race with the child's final response flush and make an otherwise clean
    # protocol peer exit with EPIPE.
    try
        close(Base.pipe_writer(io))
    catch err
        err isa InterruptException && throw_errors && rethrow()
        if !(err isa Base.IOError && err.code==Base.UV_EPIPE)
            first_error===nothing && (first_error=err)
        end
    end
    if throw_errors
        status=timedwait(() -> process_exited(io),5.0;pollint=0.01)
        if status==:timed_out
            try
                kill(io)
                wait(io)
            catch err
                err isa InterruptException && rethrow()
                first_error===nothing && (first_error=err)
            end
            first_error===nothing && (first_error=ErrorException(
                "external process did not exit within 5 seconds after its shutdown record"))
        else
            try
                success(io) || (first_error===nothing &&
                    (first_error=Base.ProcessFailedException(io)))
            catch err
                err isa InterruptException && rethrow()
                first_error===nothing && (first_error=err)
            end
        end
    elseif process_running(io)
        # A finalizer or failed/cancelled transaction cannot safely reuse this
        # transport. Terminate the process owned by this field before dropping
        # the last handle; libuv will reap it asynchronously.
        try
            kill(io)
        catch
        end
    end
    try
        close(Base.pipe_reader(io))
    catch err
        err isa InterruptException && throw_errors && rethrow()
        if !(err isa Base.IOError && err.code==Base.UV_EPIPE)
            first_error===nothing && (first_error=err)
        end
    end
    if throw_errors && first_error!==nothing
        throw(ArgumentError("ExternalProcessField: failed to close process: " *
                            sprint(showerror,first_error)))
    end
    return nothing
end

function Base.close(field::ExternalProcessField)
    lock(field.lock)
    try
        field.closed && return nothing
        field.closed=true
        return _ext_close_io!(field,true)
    finally
        unlock(field.lock)
    end
end
Base.isopen(field::ExternalProcessField)=!field.closed

function field_value(field::ExternalProcessField,x,y,z)::Float64
    x=_float_value(x,"ExternalProcessField","x")
    y=_float_value(y,"ExternalProcessField","y")
    z=_float_value(z,"ExternalProcessField","z")
    lock(field.lock)
    try
        _ext_start!(field)
        io=field.io::Base.Process
        write(io,x); write(io,y); write(io,z); flush(io)
        v=read(io,Float64)
        isfinite(v) || throw(ArgumentError("ExternalProcessField: non-finite value at ($x,$y,$z)"))
        return v
    catch err
        _ext_close_io!(field,false)
        field.closed=true
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("ExternalProcessField: protocol failure at ($x,$y,$z): $err"))
    finally
        unlock(field.lock)
    end
end

# ── .geo builders for the remaining kinds ─────────────────────────────────────

function _geo_math_fields(spec,built)
    # MathEval F<tag> references are resolved from already-built fields; the
    # caller must build dependencies first (the graph walker does that).
    return Dict{Int,AbstractField}(built)
end
