# Full-domain positivity of an actual first-order Hex8 isoparametric map.
# Its determinant has degree at most two in each reference coordinate. Values
# on the 3x3x3 grid {0,1/2,1} determine that polynomial; conversion to tensor
# Bernstein coefficients bounds its value everywhere on the reference cube.
# Strictly positive coefficients therefore certify the entire cell, rather
# than only sampled Jacobians or a tetrahedral volume decomposition.
#
# A cheap interval derivative-box test handles affine and small-turn cells.
# Otherwise outward-rounded interval Bernstein coefficients provide the fast
# certificate. Undecidable intervals use exact common-denominator BigInts.
# Nonpositive exact coefficients are conservatively rejected: this sufficient
# certificate does not claim every rejected map is folded. Work and storage
# per cell are constant (27 samples/coefficients), with no subdivision tree.

using ..Predicates: _common_ints

struct _ExtrudeNoNewJacInterval
    lo::Float64
    hi::Float64
end

const _EXTRUDE_NONEW_JAC_ZERO = _ExtrudeNoNewJacInterval(0.0,0.0)
const _EXTRUDE_NONEW_JAC_UNKNOWN = _ExtrudeNoNewJacInterval(-Inf,Inf)
const _EXTRUDE_NONEW_JAC_EDGES = ((2,1),(3,4),(6,5),(7,8),
    (4,1),(3,2),(8,5),(7,6),(5,1),(6,2),(7,3),(8,4))

@inline function _extrude_nonew_jac_bounds(lo,hi)
    (isnan(lo) || isnan(hi)) && return _EXTRUDE_NONEW_JAC_UNKNOWN
    return _ExtrudeNoNewJacInterval(prevfloat(lo),nextfloat(hi))
end

@inline function _extrude_nonew_jac_difference(a::Float64,b::Float64)
    a==b && return _EXTRUDE_NONEW_JAC_ZERO
    value=a-b
    return _extrude_nonew_jac_bounds(value,value)
end

@inline _extrude_nonew_jac_iszero(a) = a.lo==0.0 && a.hi==0.0
@inline _extrude_nonew_jac_neg(a) = _ExtrudeNoNewJacInterval(-a.hi,-a.lo)

@inline function _extrude_nonew_jac_add(a,b)
    _extrude_nonew_jac_iszero(a) && return b
    _extrude_nonew_jac_iszero(b) && return a
    return _extrude_nonew_jac_bounds(a.lo+b.lo,a.hi+b.hi)
end

@inline _extrude_nonew_jac_sub(a,b) =
    _extrude_nonew_jac_add(a,_extrude_nonew_jac_neg(b))

@inline function _extrude_nonew_jac_scale(a,k::Float64)
    k==0.0 && return _EXTRUDE_NONEW_JAC_ZERO
    k==1.0 && return a
    k<0.0 && return _extrude_nonew_jac_neg(_extrude_nonew_jac_scale(a,-k))
    return _extrude_nonew_jac_bounds(a.lo*k,a.hi*k)
end

@inline function _extrude_nonew_jac_mul(a,b)
    (_extrude_nonew_jac_iszero(a) || _extrude_nonew_jac_iszero(b)) &&
        return _EXTRUDE_NONEW_JAC_ZERO
    aa=a.lo*b.lo;ab=a.lo*b.hi;ba=a.hi*b.lo;bb=a.hi*b.hi
    return _extrude_nonew_jac_bounds(min(aa,ab,ba,bb),max(aa,ab,ba,bb))
end

@inline function _extrude_nonew_jac_det(a,b,c)
    first=_extrude_nonew_jac_mul(a[1],_extrude_nonew_jac_sub(
        _extrude_nonew_jac_mul(b[2],c[3]),_extrude_nonew_jac_mul(b[3],c[2])))
    second=_extrude_nonew_jac_mul(a[2],_extrude_nonew_jac_sub(
        _extrude_nonew_jac_mul(b[1],c[3]),_extrude_nonew_jac_mul(b[3],c[1])))
    third=_extrude_nonew_jac_mul(a[3],_extrude_nonew_jac_sub(
        _extrude_nonew_jac_mul(b[1],c[2]),_extrude_nonew_jac_mul(b[2],c[1])))
    return _extrude_nonew_jac_add(_extrude_nonew_jac_sub(first,second),third)
end

function _extrude_nonew_jac_interval_edges(v)
    return ntuple(Val(12)) do i
        a,b=_EXTRUDE_NONEW_JAC_EDGES[i]
        ntuple(d->_extrude_nonew_jac_difference(v[a][d],v[b][d]),Val(3))
    end
end

# A reference derivative is a convex combination of these four edge vectors.
# Their component bounds give a safe derivative box over the whole cube.
function _extrude_nonew_jac_derivative_box(edges,offset)
    return ntuple(Val(3)) do d
        _ExtrudeNoNewJacInterval(
            min(edges[offset+1][d].lo,edges[offset+2][d].lo,
                edges[offset+3][d].lo,edges[offset+4][d].lo),
            max(edges[offset+1][d].hi,edges[offset+2][d].hi,
                edges[offset+3][d].hi,edges[offset+4][d].hi))
    end
end

function _extrude_nonew_jac_derivative(edges,offset,weights)
    return ntuple(Val(3)) do d
        first=_extrude_nonew_jac_add(
            _extrude_nonew_jac_scale(edges[offset+1][d],weights[1]),
            _extrude_nonew_jac_scale(edges[offset+2][d],weights[2]))
        second=_extrude_nonew_jac_add(
            _extrude_nonew_jac_scale(edges[offset+3][d],weights[3]),
            _extrude_nonew_jac_scale(edges[offset+4][d],weights[4]))
        _extrude_nonew_jac_add(first,second)
    end
end

function _extrude_nonew_jac_sample(edges,index)
    r=Float64((index-1)%3)*0.5
    s=Float64((index-1)÷3%3)*0.5
    t=Float64((index-1)÷9)*0.5
    dr=_extrude_nonew_jac_derivative(edges,0,
        ((1-s)*(1-t),s*(1-t),(1-s)*t,s*t))
    ds=_extrude_nonew_jac_derivative(edges,4,
        ((1-r)*(1-t),r*(1-t),(1-r)*t,r*t))
    dt=_extrude_nonew_jac_derivative(edges,8,
        ((1-r)*(1-s),r*(1-s),r*s,(1-r)*s))
    return _extrude_nonew_jac_det(dr,ds,dt)
end

function _extrude_nonew_jac_bernstein_axis(values,stride)
    return ntuple(Val(27)) do i
        digit=(i-1)÷stride%3
        digit==1 || return values[i]
        base=i-stride
        # b1 = 2*f(1/2) - (f(0)+f(1))/2; b0/b2 are endpoint values.
        _extrude_nonew_jac_sub(_extrude_nonew_jac_scale(values[i],2.0),
            _extrude_nonew_jac_scale(_extrude_nonew_jac_add(
                values[base],values[base+2stride]),0.5))
    end
end

function _extrude_nonew_jac_exact_derivative(edges,offset,weights)
    return ntuple(d->weights[1]*edges[offset+1][d]+
        weights[2]*edges[offset+2][d]+weights[3]*edges[offset+3][d]+
        weights[4]*edges[offset+4][d],Val(3))
end

function _extrude_nonew_jac_exact_sample(edges,index)
    r=(index-1)%3;s=(index-1)÷3%3;t=(index-1)÷9
    # These integer weights represent four times each derivative. The common
    # positive coordinate scale and determinant factor 4^3 preserve signs.
    dr=_extrude_nonew_jac_exact_derivative(edges,0,
        ((2-s)*(2-t),s*(2-t),(2-s)*t,s*t))
    ds=_extrude_nonew_jac_exact_derivative(edges,4,
        ((2-r)*(2-t),r*(2-t),(2-r)*t,r*t))
    dt=_extrude_nonew_jac_exact_derivative(edges,8,
        ((2-r)*(2-s),r*(2-s),r*s,(2-r)*s))
    return dr[1]*(ds[2]*dt[3]-ds[3]*dt[2])-
           dr[2]*(ds[1]*dt[3]-ds[3]*dt[1])+
           dr[3]*(ds[1]*dt[2]-ds[2]*dt[1])
end

function _extrude_nonew_jac_exact_bernstein_axis(values,stride)
    return ntuple(Val(27)) do i
        digit=(i-1)÷stride%3
        digit==1 || return 2values[i]
        base=i-stride
        # Twice the Bernstein transform keeps every intermediate integral.
        # Applying it on three axes scales all coefficients by positive 2^3.
        4values[i]-values[base]-values[base+2stride]
    end
end

@noinline function _extrude_nonew_jac_exact_positive(v,orientation)
    ints=_common_ints(ntuple(i->v[(i-1)÷3+1][(i-1)%3+1],Val(24)))
    edges=ntuple(Val(12)) do i
        a,b=_EXTRUDE_NONEW_JAC_EDGES[i]
        ntuple(d->ints[3(a-1)+d]-ints[3(b-1)+d],Val(3))
    end
    samples=ntuple(i->_extrude_nonew_jac_exact_sample(edges,i),Val(27))
    first=_extrude_nonew_jac_exact_bernstein_axis(samples,1)
    second=_extrude_nonew_jac_exact_bernstein_axis(first,3)
    coefficients=_extrude_nonew_jac_exact_bernstein_axis(second,9)
    return all(value->orientation==1 ? value>0 : value<0,coefficients)
end

function _extrude_nonew_hex_jacobian_certify(v::NTuple{8,NTuple{3,Float64}},
        orientation::Int,caller)
    (orientation==1 || orientation==-1) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts Hex8 Jacobian certificate needs an orientation sign"))
    all(p->all(isfinite,p),v) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts Hex8 Jacobian certificate needs finite coordinates"))
    edges=_extrude_nonew_jac_interval_edges(v)
    whole=_extrude_nonew_jac_det(
        _extrude_nonew_jac_derivative_box(edges,0),
        _extrude_nonew_jac_derivative_box(edges,4),
        _extrude_nonew_jac_derivative_box(edges,8))
    (orientation==1 ? whole.lo>0 : whole.hi<0) && return nothing
    samples=ntuple(i->_extrude_nonew_jac_sample(edges,i),Val(27))
    first=_extrude_nonew_jac_bernstein_axis(samples,1)
    second=_extrude_nonew_jac_bernstein_axis(first,3)
    coefficients=_extrude_nonew_jac_bernstein_axis(second,9)
    all(value->orientation==1 ? value.lo>0 : value.hi<0,coefficients) &&
        return nothing
    _extrude_nonew_jac_exact_positive(v,orientation) && return nothing
    throw(ArgumentError("$caller: QuadTriNoNewVerts full P1 Hex8 Jacobian " *
        "positivity is not certified by tensor Bernstein bounds; " *
        "the general isoparametric cell planner is pending"))
end
