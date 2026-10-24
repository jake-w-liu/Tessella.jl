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

# ══════════════════════════════════════════════════════════════════════════════
# Mesh.HighOrderOptimize — Pyr14 admission into the certificate-gated interior
# relocation engine (OPT01 bounded slice).
#
# optimize_p2_interior! admits a block family through exactly one gate — the
# _p2opt_family lookup defined earlier in this module — evaluates per-element
# acceptance through _p2opt_certify, and samples the minimum scaled Jacobian
# through _p2opt_grads.  This include is loaded after every certified-family
# table, so the Pyr14 pieces live here beside the sign certificate they reuse:
#
#   * _P2OPT_PYR14 is the family descriptor in the same layout as the
#     _P2OPT_FAMILIES entries.  Node roles follow the gmsh type-14 order —
#     primaries 1:5 (four base corners + apex), edge mids 6:13 in the
#     interleaved order (1,2),(1,4),(1,5),(2,3),(2,5),(3,4),(3,5),(4,5), and
#     the base-quad center 14 as the sole face-center slot; triangles carry
#     no face-center mid, so there are no interior slots.
#   * _p2opt_family is extended below on every builtin concrete Integer
#     subtype — new, strictly-more-specific methods (the existing generic
#     ::Integer method stays in place and still answers for custom Integer
#     subtypes, where type 14 resolves to `nothing` exactly as before).
#     This include is loaded after the generic method was defined, and
#     module precompilation forbids overwriting it; dispatching on the
#     concrete tag types the gmsh catalog stores is the lawful extension
#     point that admits the descriptor into the single admission gate.
# ══════════════════════════════════════════════════════════════════════════════

# Per-node (i,j) index tables of _p2opt_grads(::Val{14}), (0,0) = unused:
# the base quad's quadratic-Lagrange index on {-1,0,1}² (nodes 1:4,6,7,9,11,
# 14 per _P2PYR_BASE_NODES), the bilinear index on {-1,1}² of the base
# corners feeding B4 (_P2PYR_CORNER_NODES), and of the sloping-edge mids
# feeding M4 (_P2PYR_MIDDLE_NODES).
const _P2OPT_PYR14_Q9 = ((1,1),(3,1),(3,3),(1,3),(0,0),(2,1),(1,2),
                        (0,0),(3,2),(0,0),(2,3),(0,0),(0,0),(2,2))
const _P2OPT_PYR14_B4 = ((1,1),(2,1),(2,2),(1,2),(0,0),(0,0),(0,0),
                        (0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0))
const _P2OPT_PYR14_M4 = ((0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),
                        (1,1),(0,0),(2,1),(0,0),(2,2),(1,2),(0,0))

# Collapsed coordinates (r,s,w) of the fourteen lagrange nodes — the sample
# set of the relocation engine.  u=(1-w)r, v=(1-w)s maps the physical node
# positions onto r,s ∈ [-1,1]; the apex maps every (r,s,1) to the apex, and
# the (0,0,1) axis direction is used as its sample point.
const _P2OPT_PYR14_SAMPLES = (
    (-1.,-1.,0.), (1.,-1.,0.), (1.,1.,0.), (-1.,1.,0.), (0.,0.,1.),
    (0.,-1.,0.), (-1.,0.,0.), (-1.,-1.,.5), (1.,0.,0.), (1.,-1.,.5),
    (0.,1.,0.), (1.,1.,.5), (-1.,1.,.5), (0.,0.,0.))

const _P2OPT_PYR14 = (
    kind=Val(14), msh=Int32(14), rows=14, primary=5,
    edge_slots=((6,1,2),(7,1,4),(8,1,5),(9,2,3),(10,2,5),
                (11,3,4),(12,3,5),(13,4,5)),
    face_slots=((14,(1,2,3,4)),), interior_slots=(),
    facets=((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5)),
    samples=_P2OPT_PYR14_SAMPLES, volume=true)

# Physical (u,v,w) derivative weights of the collapsed type-14 map at the
# collapsed point (r,s,w) — the per-node coefficients of the Jacobian
# columns, so Σ_k x_k·g_k yields F_u, F_v, F_w:
#   A = F_u = ∂_rG/(1-w) = (1-w)·∂_rQ9 + w·∂_r(4M4-B4)
#   B = F_v = ∂_sG/(1-w) = (1-w)·∂_sQ9 + w·∂_s(4M4-B4)
#   F_w   = G_w + r·A + s·B,
#   G_w = C = -2(1-w)·Q9 + (1-2w)·(4M4-B4) + (4w-1)·P
# All terms are polynomial in (r,s,w), so the apex limit w→1 stays finite
# for every approach direction (det(A,B,C) equals the physical Jacobian
# determinant — the r·A+s·B column operation leaves det invariant while
# giving F_w its true directional value for the scaled-Jacobian norms).
@inline function _p2opt_grads(::Val{14}, r, s, w)
    q = 1.0 - w
    lr = _p2opt_lagrange1(r); dlr = _p2opt_dlagrange1(r)
    ls = _p2opt_lagrange1(s); dls = _p2opt_dlagrange1(s)
    br = (0.5 - 0.5*r, 0.5 + 0.5*r); dbr = (-0.5, 0.5)
    bs = (0.5 - 0.5*s, 0.5 + 0.5*s); dbs = (-0.5, 0.5)
    w2 = 1.0 - 2.0*w
    return ntuple(Val(14)) do k
        au = 0.0; av = 0.0; cw = 0.0
        qi = _P2OPT_PYR14_Q9[k]
        if qi[1] != 0
            i, j = qi
            au += q*dlr[i]*ls[j]
            av += q*lr[i]*dls[j]
            cw += -2.0*q*lr[i]*ls[j]
        end
        ci = _P2OPT_PYR14_B4[k]
        if ci[1] != 0
            i, j = ci
            au += -w*dbr[i]*bs[j]
            av += -w*br[i]*dbs[j]
            cw += -w2*br[i]*bs[j]
        end
        mi = _P2OPT_PYR14_M4[k]
        if mi[1] != 0
            i, j = mi
            au += 4.0*w*dbr[i]*bs[j]
            av += 4.0*w*br[i]*dbs[j]
            cw += 4.0*w2*br[i]*bs[j]
        end
        k == 5 && (cw += 4.0*w - 1.0)
        (au, av, cw + r*au + s*av)
    end
end

# The certified relocation families by gmsh type.  The single admission
# gate is _p2opt_family(::Integer) defined earlier in this module and it
# must not be overwritten during module precompilation, so Pyr14 is
# admitted through strictly-more-specific methods on every builtin concrete
# Integer subtype — the concrete tag types a caller can actually hand the
# optimizer — while the generic method keeps serving every custom Integer
# subtype with the original _P2OPT_FAMILIES scan (type 14 stays `nothing`
# there, the pre-admission fail-closed behavior).
@inline function _p2opt_pyr14_family(msh::Integer)
    msh == 14 && return _P2OPT_PYR14
    for fam in _P2OPT_FAMILIES
        fam.msh == msh && return fam
    end
    return nothing
end
for T in (Bool, Int8, Int16, Int32, Int64, Int128,
          UInt8, UInt16, UInt32, UInt64, UInt128, BigInt)
    @eval _p2opt_family(msh::$T) = _p2opt_pyr14_family(msh)
end

# The family's exact positivity certificate — the relocation gate.  The
# existing sufficient bound (interval derivative box first, then the 144
# exact dyadic Bernstein coefficients) is evaluated at strict positive
# orientation; its throw contract becomes `false` for the optimizer, so an
# inconclusive bound freezes the element without ever asserting a folded
# element is valid.
@inline function _p2opt_certify(::Val{14}, coords, conn, t)
    v = ntuple(Val(14)) do k
        n = conn[k, t]
        (coords[1, n], coords[2, n], coords[3, n])
    end
    try
        _p2_pyramid_jacobian_certify(v, 1, "optimize_p2_interior!")
        return true
    catch err
        (err isa InterruptException || err isa OutOfMemoryError) && rethrow()
        return false
    end
end
