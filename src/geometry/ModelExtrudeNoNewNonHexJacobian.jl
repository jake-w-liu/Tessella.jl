# Positivity of the complete P1 pyramid and prism reference maps. These
# certificates use actual emitted corners; a positive tetrahedral partition
# alone does not certify a warped isoparametric cell.
#
# In collapsed pyramid coordinates X(r,s,t)=(1-t)B(r,s)+tA, the physical
# reference Jacobian is (B_r x B_s).(A-B). It is independent of t and bilinear
# in r,s. Its four corner values are oriented tetrahedral determinants, so
# exact positive corner signs certify the whole reference domain.
#
# For a prism, X(u,v,t) interpolates its two triangles. Its determinant is
# affine in u,v and quadratic in t. Three triangle-vertex quadratics therefore
# determine positivity everywhere. Their degree-two Bernstein coefficients
# are b0, h/2, b2. Positive endpoints and h>=0 suffice; if h<0, positivity is
# equivalent to 4*b0*b2-h^2>0. Outward intervals give the usual allocation-free
# certificate; ambiguous intervals use exact common-denominator integers.

const _EXTRUDE_NONEW_PYRAMID_JAC_CORNERS =
    ((1,2,4),(2,3,1),(3,4,2),(4,1,3))

@inline function _extrude_nonew_nonhex_finite(v,caller,family)
    for point in v,value in point
        isfinite(value) || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts $family Jacobian certificate " *
            "needs finite coordinates"))
    end
    return nothing
end

function _extrude_nonew_pyramid_jacobian_certify(
        v::NTuple{5,NTuple{3,Float64}},caller)
    _extrude_nonew_nonhex_finite(v,caller,"Pyr5")
    first=-orient3(v[1],v[2],v[3],v[5])
    orientation=first>0 ? 1 : first<0 ? -1 : 0
    if orientation!=0
        for (a,b,c) in _EXTRUDE_NONEW_PYRAMID_JAC_CORNERS
            -orient3(v[a],v[b],v[c],v[5])*orientation>0 ||
                throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts full P1 Pyr5 Jacobian " *
                    "positivity is not certified by exact base-corner signs"))
        end
        return nothing
    end
    throw(ArgumentError(
        "$caller: QuadTriNoNewVerts full P1 Pyr5 Jacobian " *
        "positivity is not certified by exact base-corner signs"))
end

@inline _extrude_nonew_prism_jac_difference(a,b) =
    ntuple(k->_extrude_nonew_jac_difference(a[k],b[k]),Val(3))

@inline function _extrude_nonew_prism_jac_oriented(value,orientation)
    return orientation==1 ? value : _extrude_nonew_jac_neg(value)
end

@inline function _extrude_nonew_prism_jac_interval_positive(
        bottom_first,bottom_second,top_first,top_second,column,orientation)
    b0=_extrude_nonew_prism_jac_oriented(_extrude_nonew_jac_det(
        bottom_first,bottom_second,column),orientation)
    b2=_extrude_nonew_prism_jac_oriented(_extrude_nonew_jac_det(
        top_first,top_second,column),orientation)
    (b0.lo>0 && b2.lo>0) || return false
    h=_extrude_nonew_prism_jac_oriented(_extrude_nonew_jac_add(
        _extrude_nonew_jac_det(bottom_first,top_second,column),
        _extrude_nonew_jac_det(top_first,bottom_second,column)),orientation)
    h.lo>=0 && return true
    h.hi<0 || return false
    discriminant=_extrude_nonew_jac_sub(
        _extrude_nonew_jac_scale(_extrude_nonew_jac_mul(b0,b2),4.0),
        _extrude_nonew_jac_mul(h,h))
    return discriminant.lo>0
end

@inline function _extrude_nonew_prism_jac_exact_det(a,b,c)
    return a[1]*(b[2]*c[3]-b[3]*c[2])-
           a[2]*(b[1]*c[3]-b[3]*c[1])+
           a[3]*(b[1]*c[2]-b[2]*c[1])
end

@noinline function _extrude_nonew_prism_jac_exact_positive(v,orientation)
    ints=_common_ints(ntuple(i->v[(i-1)÷3+1][(i-1)%3+1],Val(18)))
    bottom_first=ntuple(k->ints[3+k]-ints[k],Val(3))
    bottom_second=ntuple(k->ints[6+k]-ints[k],Val(3))
    top_first=ntuple(k->ints[12+k]-ints[9+k],Val(3))
    top_second=ntuple(k->ints[15+k]-ints[9+k],Val(3))
    for vertex in 1:3
        column=ntuple(k->ints[3(vertex+2)+k]-ints[3(vertex-1)+k],Val(3))
        b0=_extrude_nonew_prism_jac_exact_det(
            bottom_first,bottom_second,column)*orientation
        b2=_extrude_nonew_prism_jac_exact_det(
            top_first,top_second,column)*orientation
        (b0>0 && b2>0) || return false
        h=(_extrude_nonew_prism_jac_exact_det(
               bottom_first,top_second,column)+
           _extrude_nonew_prism_jac_exact_det(
               top_first,bottom_second,column))*orientation
        (h>=0 || 4b0*b2>h*h) || return false
    end
    return true
end

function _extrude_nonew_prism_jacobian_certify(
        v::NTuple{6,NTuple{3,Float64}},caller)
    _extrude_nonew_nonhex_finite(v,caller,"Pri6")
    first=-orient3(v[1],v[2],v[3],v[6])
    orientation=first>0 ? 1 : first<0 ? -1 : 0
    if orientation!=0
        bottom_first=_extrude_nonew_prism_jac_difference(v[2],v[1])
        bottom_second=_extrude_nonew_prism_jac_difference(v[3],v[1])
        top_first=_extrude_nonew_prism_jac_difference(v[5],v[4])
        top_second=_extrude_nonew_prism_jac_difference(v[6],v[4])
        certified=true
        for vertex in 1:3
            column=_extrude_nonew_prism_jac_difference(v[vertex+3],v[vertex])
            if !_extrude_nonew_prism_jac_interval_positive(
                    bottom_first,bottom_second,top_first,top_second,column,orientation)
                certified=false
                break
            end
        end
        (certified || _extrude_nonew_prism_jac_exact_positive(v,orientation)) &&
            return nothing
    end
    throw(ArgumentError(
        "$caller: QuadTriNoNewVerts full P1 Pri6 Jacobian " *
        "positivity is not certified by the exact vertex-quadratic criterion"))
end
