# Whole-domain sign certificate for the actual Gmsh type-14 rational map.
# In collapsed coordinates u=(1-w)r, v=(1-w)s, r,s in [-1,1], w in [0,1],
# its polynomial map is
# G = (1-w)^2 Q9 + w(1-w)(4M4-B4) + w(2w-1)P,
# where Q9 is the quadratic base, B4 its bilinear corner interpolant, M4 the
# bilinear interpolant of the four sloping-edge midpoints, and P the apex.
# A=G_r/(1-w), B=G_s/(1-w), C=G_w have tensor degrees (1,2,1), (2,1,1),
# (2,2,1). det(A,B,C) equals the physical reference Jacobian determinant;
# adding r*A+s*B to C recovers the physical w derivative without changing it.
# Its degree-(5,5,3) Bernstein coefficients therefore bound the entire pyramid,
# including all limiting apex directions. No sample-only admission is used.
#
# A translation-free outward interval derivative-box certificate handles
# regular elements without heap allocation. Otherwise 144 exact dyadic integer
# Bernstein coefficients provide a sufficient certificate. An inconclusive
# bound is rejected without asserting that every rejected map is folded.
# Storage and work per element are bounded; there is no subdivision search.

using ..Predicates: _common_ints

struct _P2JacInterval
    lo::Float64
    hi::Float64
end

const _P2_JAC_ZERO = _P2JacInterval(0.0,0.0)
const _P2_JAC_UNKNOWN = _P2JacInterval(-Inf,Inf)

@inline function _p2_jac_interval_bounds(lo,hi)
    (isnan(lo) || isnan(hi)) && return _P2_JAC_UNKNOWN
    return _P2JacInterval(prevfloat(lo),nextfloat(hi))
end

@inline function _p2_jac_interval_difference(a::Float64,b::Float64)
    a==b && return _P2_JAC_ZERO
    value=a-b
    return _p2_jac_interval_bounds(value,value)
end

# Binary scaling is exact unless a result overflows or loses subnormal bits.
# Normal results retain the original significand; an exactly reversible
# subnormal scaling is also admitted as a point interval. Other results need
# outward endpoints. No floating scale factor is formed, even at exponent1074.
@inline function _p2_jac_interval_scaled_value(value::Float64,power::Int)
    scaled=ldexp(value,power)
    if isfinite(scaled) && (value==0.0 || abs(scaled)>=floatmin(Float64) ||
                           ldexp(scaled,-power)==value)
        return _P2JacInterval(scaled,scaled)
    end
    return _p2_jac_interval_bounds(scaled,scaled)
end

@inline function _p2_jac_interval_ldexp(value,power::Int)
    value.lo==value.hi && return _p2_jac_interval_scaled_value(value.lo,power)
    return _p2_jac_interval_bounds(ldexp(value.lo,power),ldexp(value.hi,power))
end

@inline _p2_jac_interval_iszero(a)=a.lo==0.0 && a.hi==0.0
@inline _p2_jac_interval_neg(a)=_P2JacInterval(-a.hi,-a.lo)

@inline function _p2_jac_interval_add(a,b)
    _p2_jac_interval_iszero(a) && return b
    _p2_jac_interval_iszero(b) && return a
    return _p2_jac_interval_bounds(a.lo+b.lo,a.hi+b.hi)
end

@inline _p2_jac_interval_sub(a,b)=
    _p2_jac_interval_add(a,_p2_jac_interval_neg(b))

@inline function _p2_jac_interval_scale(a,k::Float64)
    k==0.0 && return _P2_JAC_ZERO
    k==1.0 && return a
    k<0.0 && return _p2_jac_interval_neg(_p2_jac_interval_scale(a,-k))
    return _p2_jac_interval_bounds(a.lo*k,a.hi*k)
end

@inline function _p2_jac_interval_mul(a,b)
    (_p2_jac_interval_iszero(a) || _p2_jac_interval_iszero(b)) &&
        return _P2_JAC_ZERO
    aa=a.lo*b.lo;ab=a.lo*b.hi;ba=a.hi*b.lo;bb=a.hi*b.hi
    return _p2_jac_interval_bounds(min(aa,ab,ba,bb),max(aa,ab,ba,bb))
end

@inline function _p2_jac_interval_det(a,b,c)
    first=_p2_jac_interval_mul(a[1],_p2_jac_interval_sub(
        _p2_jac_interval_mul(b[2],c[3]),_p2_jac_interval_mul(b[3],c[2])))
    second=_p2_jac_interval_mul(a[2],_p2_jac_interval_sub(
        _p2_jac_interval_mul(b[1],c[3]),_p2_jac_interval_mul(b[3],c[1])))
    third=_p2_jac_interval_mul(a[3],_p2_jac_interval_sub(
        _p2_jac_interval_mul(b[1],c[2]),_p2_jac_interval_mul(b[2],c[1])))
    return _p2_jac_interval_add(_p2_jac_interval_sub(first,second),third)
end

const _P2PYR_BASE_NODES = ((1,1,1),(2,3,1),(3,3,3),(4,1,3),
    (6,2,1),(9,3,2),(11,2,3),(7,1,2),(14,2,2))
const _P2PYR_CORNER_NODES = ((1,1,1),(2,2,1),(3,2,2),(4,1,2))
const _P2PYR_MIDDLE_NODES = ((8,1,1),(10,2,1),(12,2,2),(13,1,2))
const _P2PYR_BILINEAR_VAL = ((2,1,0),(0,1,2))

# All derivative coefficients are integer weighted node sums divided by four.
# Each row sums to zero, so exact translation removal is valid before either
# interval or common-integer arithmetic. These immutable tables are private.
function _p2pyr_derivative_weights(axis::Int,index::Int)
    weights=zeros(Int16,14)
    nr=axis==1 ? 2 : 3
    ns=axis==2 ? 2 : 3
    ir=(index-1)÷(2ns)
    js=(index-1)÷2%ns
    kw=(index-1)%2
    0<=ir<nr || throw(ErrorException("invalid Pyramid14 derivative index"))
    if kw==0
        for (node,a,b) in _P2PYR_BASE_NODES
            value=axis==1 ? _P2L2_DER[a][ir+1]*_P2L2_VAL[b][js+1] :
                  axis==2 ? _P2L2_VAL[a][ir+1]*_P2L2_DER[b][js+1] :
                            -2*_P2L2_VAL[a][ir+1]*_P2L2_VAL[b][js+1]
            weights[node]+=value
        end
    end
    if axis==3 || kw==1
        for (node,a,b) in _P2PYR_CORNER_NODES
            value=axis==1 ? (2a-3)*_P2PYR_BILINEAR_VAL[b][js+1] :
                  axis==2 ? _P2PYR_BILINEAR_VAL[a][ir+1]*(2b-3) :
                            _P2PYR_BILINEAR_VAL[a][ir+1]*_P2PYR_BILINEAR_VAL[b][js+1]
            weights[node]+=(axis==3 && kw==1 ? 1 : -1)*value
        end
        for (node,a,b) in _P2PYR_MIDDLE_NODES
            value=axis==1 ? (2a-3)*_P2PYR_BILINEAR_VAL[b][js+1] :
                  axis==2 ? _P2PYR_BILINEAR_VAL[a][ir+1]*(2b-3) :
                            _P2PYR_BILINEAR_VAL[a][ir+1]*_P2PYR_BILINEAR_VAL[b][js+1]
            weights[node]+=(axis==3 && kw==1 ? -4 : 4)*value
        end
    end
    axis==3 && (weights[5]+=kw==0 ? -4 : 12)
    sum(weights)==0 || throw(ErrorException("Pyramid14 derivative weights do not conserve constants"))
    return Tuple(weights)
end

const _P2PYR_DR_WEIGHTS=ntuple(i->_p2pyr_derivative_weights(1,i),12)
const _P2PYR_DS_WEIGHTS=ntuple(i->_p2pyr_derivative_weights(2,i),12)
const _P2PYR_DW_WEIGHTS=ntuple(i->_p2pyr_derivative_weights(3,i),18)

@inline function _p2pyr_interval_row(points,weights)
    return ntuple(Val(3)) do d
        value=_P2_JAC_ZERO
        for node in 2:14
            weight=weights[node]
            weight==0 && continue
            value=_p2_jac_interval_add(value,
                _p2_jac_interval_scale(points[node][d],Float64(weight)))
        end
        _p2_jac_interval_scale(value,0.25)
    end
end

function _p2pyr_interval_scale(v,axis)
    magnitude=0.0
    for node in 1:14
        magnitude=max(magnitude,abs(v[node][axis]))
    end
    magnitude==0.0 && return 0
    _,power=frexp(magnitude)
    return 1-power
end

@inline function _p2pyr_interval_difference(a,b,power)
    a==b && return _P2_JAC_ZERO
    first=_p2_jac_interval_scaled_value(a,power)
    second=_p2_jac_interval_scaled_value(b,power)
    if first.lo==first.hi && second.lo==second.hi
        return _p2_jac_interval_difference(first.lo,second.lo)
    end
    return _p2_jac_interval_sub(first,second)
end

function _p2pyr_interval_derivatives(v,normalize::Bool=false)
    powers=normalize ? ntuple(d->_p2pyr_interval_scale(v,d),Val(3)) : (0,0,0)
    points=ntuple(Val(14)) do node
        ntuple(d->_p2pyr_interval_difference(v[node][d],v[1][d],powers[d]),Val(3))
    end
    first=ntuple(i->_p2pyr_interval_row(points,_P2PYR_DR_WEIGHTS[i]),Val(12))
    second=ntuple(i->_p2pyr_interval_row(points,_P2PYR_DS_WEIGHTS[i]),Val(12))
    third=ntuple(i->_p2pyr_interval_row(points,_P2PYR_DW_WEIGHTS[i]),Val(18))
    return first,second,third
end

function _p2pyr_interval_box(values)
    return ntuple(Val(3)) do d
        lo=values[1][d].lo;hi=values[1][d].hi
        for i in 2:length(values)
            lo=min(lo,values[i][d].lo)
            hi=max(hi,values[i][d].hi)
        end
        _P2JacInterval(lo,hi)
    end
end

function _p2pyr_interval_bound(v,normalize::Bool=false)
    first,second,third=_p2pyr_interval_derivatives(v,normalize)
    return _p2_jac_interval_det(_p2pyr_interval_box(first),
        _p2pyr_interval_box(second),_p2pyr_interval_box(third))
end

function _p2pyr_exact_derivatives(points,rows)
    values=Matrix{BigInt}(undef,3,length(rows))
    for index in eachindex(rows),d in 1:3
        value=zero(BigInt)
        for node in 2:14
            weight=rows[index][node]
            weight==0 && continue
            value+=weight*points[3(node-1)+d]
        end
        values[d,index]=value
    end
    return values
end

# The returned S[i,j,k] carries the positive denominator
# 64*C(5,i-1)*C(5,j-1)*C(3,k-1), and a common positive dyadic coordinate
# scale cubed. Keeping the numerator integral preserves even subnormal signs.
@noinline function _p2pyr_exact_coefficients(v)
    integers=_common_ints(ntuple(i->v[(i-1)÷3+1][(i-1)%3+1],Val(42)))
    points=ntuple(i->integers[i]-integers[(i-1)%3+1],Val(42))
    first=_p2pyr_exact_derivatives(points,_P2PYR_DR_WEIGHTS)
    second=_p2pyr_exact_derivatives(points,_P2PYR_DS_WEIGHTS)
    third=_p2pyr_exact_derivatives(points,_P2PYR_DW_WEIGHTS)
    # Form A x B once (degree3,3,2), then its dot product with C. Deferring
    # each product's binomial denominators avoids 2592 separate determinants.
    cross=zeros(BigInt,3,4,4,3)
    for ir in 0:1,js in 0:2,kw in 0:1
        a=(first[1,ir*6+js*2+kw+1],first[2,ir*6+js*2+kw+1],
           first[3,ir*6+js*2+kw+1])
        for ir2 in 0:2,js2 in 0:1,kw2 in 0:1
            b=(second[1,ir2*4+js2*2+kw2+1],second[2,ir2*4+js2*2+kw2+1],
               second[3,ir2*4+js2*2+kw2+1])
            product=_bigint_cross(a,b)
            weight=_P2B2[ir2+1]*_P2B2[js+1]
            for d in 1:3
                cross[d,ir+ir2+1,js+js2+1,kw+kw2+1]+=weight*product[d]
            end
        end
    end
    coefficients=zeros(BigInt,6,6,4)
    for ir in 0:3,js in 0:3,kw in 0:2
        a=(cross[1,ir+1,js+1,kw+1],cross[2,ir+1,js+1,kw+1],cross[3,ir+1,js+1,kw+1])
        for ir2 in 0:2,js2 in 0:2,kw2 in 0:1
            position=ir2*6+js2*2+kw2+1
            value=a[1]*third[1,position]+a[2]*third[2,position]+a[3]*third[3,position]
            weight=_P2B2[ir2+1]*_P2B2[js2+1]
            coefficients[ir+ir2+1,js+js2+1,kw+kw2+1]+=weight*value
        end
    end
    return coefficients
end

function _p2_pyramid_jacobian_certify(v::NTuple{14,NTuple{3,Float64}},
                                      orientation::Int,caller)
    (orientation==1 || orientation==-1) || throw(ArgumentError(
        "$caller: Pyramid14 Jacobian certificate needs an orientation sign"))
    all(p->all(isfinite,p),v) || throw(ArgumentError(
        "$caller: Pyramid14 Jacobian certificate needs finite coordinates"))
    # Independent positive power-of-two coordinate scales preserve the
    # determinant sign and keep ordinary extreme-scale cells on the fast path.
    bound=_p2pyr_interval_bound(v,true)
    (orientation==1 ? bound.lo>0 : bound.hi<0) && return nothing
    coefficients=_p2pyr_exact_coefficients(v)
    all(value->orientation==1 ? value>0 : value<0,coefficients) && return nothing
    throw(ArgumentError("$caller: full P2 Pyramid14 Jacobian positivity is not " *
        "certified by tensor Bernstein bounds; the general isoparametric " *
        "placement planner is pending"))
end
