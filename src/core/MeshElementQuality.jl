"""
    MeshElementQuality

Scale-robust Gmsh-shaped quality measures for dense linear-simplex elements in a
finalized [`Mesh`](@ref). The implementation follows the Gmsh 4.15.2 definitions
for type-1 segments, type-2 triangles, and type-4 tetrahedra while rejecting the
four 1-D measures that Gmsh itself does not define reliably.
"""
module MeshElementQuality

using ..MeshTypes: Mesh, node, triangle_area, tet_signed_volume
using ..MeshPointLocation: mesh_element_offsets

export mesh_element_quality, mesh_element_qualities

const ELEMENT_QUALITY_NAMES=(
    "minDetJac","maxDetJac","minSJ","minSICN","minSIGE","gamma",
    "innerRadius","outerRadius","minIsotropy","angleShape","minEdge",
    "maxEdge","volume")
const _UNDEFINED_SEGMENT_QUALITIES=("minDetJac","maxDetJac","minSIGE",
                                    "minIsotropy")
const _Exact=Rational{BigInt}
const _FILTER_FACTOR=sqrt(eps(Float64))

struct _FrameScale
    coefficient::Float64
    exponent::Int
end

@inline _sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _scale3(a,s)=(a[1]*s,a[2]*s,a[3]*s)
@inline _dot3(a,b)=muladd(a[1],b[1],muladd(a[2],b[2],a[3]*b[3]))
@inline _norm3(a)=hypot(a[1],a[2],a[3])
@inline _cross3(a,b)=(a[2]*b[3]-a[3]*b[2],
                      a[3]*b[1]-a[1]*b[3],
                      a[1]*b[2]-a[2]*b[1])
@inline _det3(a,b,c)=_dot3(a,_cross3(b,c))
@inline _maxabs3(a)=max(abs(a[1]),abs(a[2]),abs(a[3]))
@inline _finite3(a)=isfinite(a[1])&&isfinite(a[2])&&isfinite(a[3])

function _quality_name(value,caller::AbstractString)
    value isa AbstractString || throw(ArgumentError(
        "$caller: quality_name must be a string"))
    name=String(value)
    occursin('\0',name) && throw(ArgumentError(
        "$caller: quality_name must not contain NUL"))
    name in ELEMENT_QUALITY_NAMES || throw(ArgumentError(
        "$caller: unknown quality_name $(repr(name)); expected one of " *
        join(repr.(ELEMENT_QUALITY_NAMES),", ")))
    return name
end

function _quality_tags(mesh::Mesh,values,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: element_tags must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    _,_,total=mesh_element_offsets(mesh)
    tags=Vector{Int}(undef,length(values))
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: element tag $index must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: element tag $index must not be Bool"))
        tag=try
            Int(value)
        catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: element tag $index exceeds the platform Int range"))
        end
        1<=tag<=total || throw(ArgumentError(
            "$caller: unknown element tag $value; expected a dense tag in 1:$total"))
        tags[index]=tag
    end
    return tags
end

function _exact_frame(points)
    anchor=ntuple(i->_Exact(points[1][i]),3)
    edges=map(Base.tail(points)) do point
        ntuple(i->_Exact(point[i])-anchor[i],3)
    end
    scale=maximum(_maxabs3,edges;init=zero(_Exact))
    scale==0 && return _FrameScale(0.0,0),map(_->(0.0,0.0,0.0),edges)
    normalized=setprecision(BigFloat,256) do
        map(edges) do edge
            ntuple(i->Float64(BigFloat(edge[i]/scale)),3)
        end
    end
    frame_scale=setprecision(BigFloat,256) do
        coefficient,exponent=frexp(BigFloat(scale))
        _FrameScale(Float64(coefficient),exponent)
    end
    # The exponent form can represent a span above `floatmax`: a shorter edge can
    # still have finite physical length even when the longest difference does not.
    return frame_scale,normalized
end

function _normalized_frame(points)
    anchor=points[1]
    edges=map(point->_sub3(point,anchor),Base.tail(points))
    if all(_finite3,edges)
        scale=maximum(_maxabs3,edges;init=0.0)
        scale==0 && return _FrameScale(0.0,0),map(_->(0.0,0.0,0.0),edges)
        reciprocal=inv(scale)
        normalized=map(edge->_scale3(edge,reciprocal),edges)
        all(_finite3,normalized) && return _FrameScale(scale,0),normalized
    end
    return _exact_frame(points)
end

function _exact_frame_scale(points)
    anchor=ntuple(i->_Exact(points[1][i]),3)
    scale=zero(_Exact)
    @inbounds for point in points[2:end],dimension in 1:3
        scale=max(scale,abs(_Exact(point[dimension])-anchor[dimension]))
    end
    return scale
end

function _exact_cross_norm_scaled(a,b,c,frame)
    ae=ntuple(i->_Exact(a[i]),3)
    edge1=ntuple(i->_Exact(b[i])-ae[i],3)
    edge2=ntuple(i->_Exact(c[i])-ae[i],3)
    cross=(edge1[2]*edge2[3]-edge1[3]*edge2[2],
           edge1[3]*edge2[1]-edge1[1]*edge2[3],
           edge1[1]*edge2[2]-edge1[2]*edge2[1])
    squared=_dot3(cross,cross)
    squared==0 && return 0.0
    scale=_exact_frame_scale(frame)
    scale>0 || return 0.0
    return setprecision(BigFloat,256) do
        value=Float64(sqrt(BigFloat(squared))/BigFloat(scale)^2)
        value==0.0 ? nextfloat(0.0) : value
    end
end

function _cross_norm_filtered(edge1,edge2,a,b,c,frame)
    cross=_cross3(edge1,edge2)
    safe=true
    @inbounds for (component,i,j,k,l) in
        ((1,2,3,3,2),(2,3,1,1,3),(3,1,2,2,1))
        permanent=abs(edge1[i]*edge2[j])+abs(edge1[k]*edge2[l])
        safe &= isfinite(permanent) &&
                (permanent==0 ||
                 abs(cross[component])>16eps(Float64)*permanent)
    end
    magnitude=_norm3(cross)
    if safe && isfinite(magnitude) && magnitude>0
        return magnitude
    end
    return _exact_cross_norm_scaled(a,b,c,frame)
end

function _exact_det_scaled(points)
    anchor=ntuple(i->_Exact(points[1][i]),3)
    edges=ntuple(3) do index
        ntuple(i->_Exact(points[index+1][i])-anchor[i],3)
    end
    determinant=_det3(edges...)
    determinant==0 && return 0.0
    scale=maximum(_maxabs3,edges)
    return setprecision(BigFloat,256) do
        value=Float64(BigFloat(determinant)/BigFloat(scale)^3)
        if value==0.0
            return determinant<0 ? -nextfloat(0.0) : nextfloat(0.0)
        end
        value
    end
end

function _big_frame(points)
    anchor=ntuple(i->_Exact(points[1][i]),3)
    exact_edges=ntuple(length(points)-1) do index
        ntuple(i->_Exact(points[index+1][i])-anchor[i],3)
    end
    scale=maximum(_maxabs3,exact_edges;init=zero(_Exact))
    scale==0 && return scale,ntuple(
        _->(BigFloat(0),BigFloat(0),BigFloat(0)),length(exact_edges))
    edges=setprecision(BigFloat,256) do
        denominator=BigFloat(scale)
        ntuple(length(exact_edges)) do index
            ntuple(i->BigFloat(exact_edges[index][i])/denominator,3)
        end
    end
    return scale,edges
end

function _float_preserve(value::BigFloat)
    result=Float64(value)
    if result==0.0 && value!=0
        return signbit(value) ? -nextfloat(0.0) : nextfloat(0.0)
    end
    return result
end

function _big_rescale(scale,value::BigFloat)
    value==0 && return copysign(0.0,Float64(value))
    isinf(value) && return Float64(value)
    return setprecision(BigFloat,256) do
        _float_preserve(BigFloat(scale)*value)
    end
end

@inline function _big_clamp(value,lower,upper)
    return min(max(value,BigFloat(lower)),BigFloat(upper))
end

function _triangle_quality_big(name,a,b,c)
    return setprecision(BigFloat,256) do
    scale,edges=_big_frame((a,b,c))
    edge1,edge2=edges;edge3=_sub3(edge2,edge1)
    lengths=(_norm3(edge1),_norm3(edge2),_norm3(edge3))
    area2=_norm3(_cross3(edge1,edge2))
    area2>0 || return name=="outerRadius" ? Inf : 0.0
    if name in ("minSICN","minIsotropy")
        value=2sqrt(BigFloat(3))*area2/sum(abs2,lengths)
        return _float_preserve(_big_clamp(value,0,1))
    elseif name=="minSIGE"
        l1,l2,l3=lengths
        value=(2/sqrt(BigFloat(3)))*area2*
              (inv(l1*l2)+inv(l1*l3)+inv(l2*l3))/3
        return _float_preserve(_big_clamp(value,0,1))
    elseif name=="gamma"
        l1,l2,l3=lengths
        value=4area2^2/((l1+l2+l3)*l1*l2*l3)
        return _float_preserve(_big_clamp(value,0,1))
    elseif name=="innerRadius"
        return _big_rescale(scale,area2/sum(lengths))
    elseif name=="outerRadius"
        return _big_rescale(scale,prod(lengths)/(2area2))
    elseif name=="angleShape"
        angles=(atan(_norm3(_cross3(edge1,edge2)),_dot3(edge1,edge2)),
                atan(_norm3(_cross3(_scale3(edge1,-1),edge3)),
                     _dot3(_scale3(edge1,-1),edge3)),
                atan(_norm3(_cross3(_scale3(edge2,-1),_scale3(edge3,-1))),
                     _dot3(_scale3(edge2,-1),_scale3(edge3,-1))))
        sharpness=BigFloat(500);pi_big=BigFloat(pi)
        denominator=2atan(sharpness*pi_big/9)
        value=minimum(angles) do angle
            offset=angle-pi_big/3
            (atan(sharpness*(offset+pi_big/9))+
             atan(sharpness*(pi_big/9-offset)))/denominator
        end
        return _float_preserve(value)
    elseif name=="minEdge"
        return _big_rescale(scale,minimum(lengths))
    elseif name=="maxEdge"
        return _big_rescale(scale,maximum(lengths))
    end
    error("internal BigFloat triangle quality dispatch failure for $name")
    end
end

function _tet_quality_big(name,a,b,c,d)
    return setprecision(BigFloat,256) do
    scale,edges=_big_frame((a,b,c,d))
    edge1,edge2,edge3=edges
    edge4=_sub3(edge2,edge1);edge5=_sub3(edge3,edge1)
    edge6=_sub3(edge3,edge2)
    lengths=(_norm3(edge1),_norm3(edge2),_norm3(edge3),
             _norm3(edge4),_norm3(edge5),_norm3(edge6))
    determinant=_det3(edge1,edge2,edge3)
    determinant!=0 || return name=="outerRadius" ? Inf : 0.0
    if name=="minSICN"
        jacobian1=edge1
        jacobian2=ntuple(i->(-edge1[i]+2edge2[i])/sqrt(BigFloat(3)),3)
        jacobian3=ntuple(i->sqrt(BigFloat(3)/2)*edge3[i]-
                            (3/(2sqrt(BigFloat(6))))*jacobian1[i]-
                            (1/(2sqrt(BigFloat(2))))*jacobian2[i],3)
        determinant_ideal=sqrt(BigFloat(2))*determinant
        norm_j=sum(abs2,jacobian1)+sum(abs2,jacobian2)+sum(abs2,jacobian3)
        norm_cofactor=sum(abs2,_cross3(jacobian2,jacobian3))+
                       sum(abs2,_cross3(jacobian3,jacobian1))+
                       sum(abs2,_cross3(jacobian1,jacobian2))
        value=3determinant_ideal/sqrt(norm_j*norm_cofactor)
        return _float_preserve(_big_clamp(value,-1,1))
    elseif name=="minSIGE"
        l0,l1,l2,l3,l4,l5=lengths
        reciprocal_sum=
            inv(l0*l5*l1)+inv(l0*l5*l2)+inv(l0*l5*l3)+inv(l0*l5*l4)+
            inv(l1*l4*l0)+inv(l1*l4*l2)+inv(l1*l4*l3)+inv(l1*l4*l5)+
            inv(l2*l3*l0)+inv(l2*l3*l1)+inv(l2*l3*l4)+inv(l2*l3*l5)
        value=sqrt(BigFloat(2))*determinant*reciprocal_sum/12
        return _float_preserve(_big_clamp(value,-1,1))
    elseif name=="minIsotropy"
        determinant>0 || return 0.0
        value=12*(determinant/2)^(BigFloat(2)/3)/sum(abs2,lengths)
        return _float_preserve(_big_clamp(value,0,1))
    end
    face_area2=(_norm3(_cross3(edge1,edge2)),
                _norm3(_cross3(edge2,edge3)),
                _norm3(_cross3(edge1,edge3)),
                _norm3(_cross3(edge4,edge5)))
    inradius=determinant/sum(face_area2)
    numerator=ntuple(3) do i
        _dot3(edge1,edge1)*_cross3(edge2,edge3)[i]+
        _dot3(edge2,edge2)*_cross3(edge3,edge1)[i]+
        _dot3(edge3,edge3)*_cross3(edge1,edge2)[i]
    end
    circumradius=_norm3(_scale3(numerator,inv(2determinant)))
    if name=="gamma"
        value=3abs(inradius)/circumradius
        return _float_preserve(_big_clamp(value,0,1))
    elseif name=="innerRadius"
        return _big_rescale(scale,inradius)
    elseif name=="outerRadius"
        return _big_rescale(scale,circumradius)
    elseif name=="minEdge"
        return _big_rescale(scale,minimum(lengths))
    elseif name=="maxEdge"
        return _big_rescale(scale,maximum(lengths))
    end
    error("internal BigFloat tetrahedron quality dispatch failure for $name")
    end
end

function _det_filtered(edge1,edge2,edge3,points)
    determinant=_det3(edge1,edge2,edge3)
    permanent=
        abs(edge1[1])*(abs(edge2[2]*edge3[3])+abs(edge2[3]*edge3[2]))+
        abs(edge1[2])*(abs(edge2[1]*edge3[3])+abs(edge2[3]*edge3[1]))+
        abs(edge1[3])*(abs(edge2[1]*edge3[2])+abs(edge2[2]*edge3[1]))
    if isfinite(determinant)&&isfinite(permanent)&&
       abs(determinant)>_FILTER_FACTOR*permanent
        return determinant,false
    end
    return _exact_det_scaled(points),true
end

@inline function _rescale(scale::_FrameScale,value::Float64)
    value==0 && return copysign(0.0,value)
    isinf(value) && return value
    scale.exponent==0 && return scale.coefficient*value
    return setprecision(BigFloat,256) do
        _float_preserve(ldexp(BigFloat(scale.coefficient)*BigFloat(value),
                              scale.exponent))
    end
end

function _segment_geometry(a,b)
    scale,edges=_normalized_frame((a,b))
    normalized_length=_norm3(only(edges))
    length=_rescale(scale,normalized_length)
    return length,normalized_length>0
end

function _triangle_geometry(a,b,c)
    points=(a,b,c)
    scale,edges=_normalized_frame(points)
    edge1,edge2=edges
    edge3=_sub3(edge2,edge1)
    lengths=(_norm3(edge1),_norm3(edge2),_norm3(edge3))
    area2=_cross_norm_filtered(edge1,edge2,a,b,c,points)
    return scale,lengths,area2,area2>0
end

function _tetrahedron_geometry(a,b,c,d)
    points=(a,b,c,d)
    scale,edges=_normalized_frame(points)
    edge1,edge2,edge3=edges
    edge4=_sub3(edge2,edge1)
    edge5=_sub3(edge3,edge1)
    edge6=_sub3(edge3,edge2)
    lengths=(_norm3(edge1),_norm3(edge2),_norm3(edge3),
             _norm3(edge4),_norm3(edge5),_norm3(edge6))
    determinant,exact=_det_filtered(edge1,edge2,edge3,points)
    face_area2=(
        _cross_norm_filtered(edge1,edge2,a,b,c,points),
        _cross_norm_filtered(edge2,edge3,a,c,d,points),
        _cross_norm_filtered(edge1,edge3,a,b,d,points),
        _cross_norm_filtered(edge4,edge5,b,c,d,points))
    return scale,edges,lengths,determinant,face_area2,exact,determinant!=0
end

@inline function _triangle_isotropy(lengths,area2)
    denominator=sum(abs2,lengths)
    denominator>0 || return 0.0
    return clamp(2sqrt(3.0)*area2/denominator,0.0,1.0)
end

function _triangle_sige(lengths,area2)
    l1,l2,l3=lengths
    minimum(lengths)>0 || return 0.0
    value=(2/sqrt(3.0))*area2*
          (inv(l1*l2)+inv(l1*l3)+inv(l2*l3))/3
    return clamp(value,0.0,1.0)
end

function _triangle_gamma(lengths,area2)
    l1,l2,l3=lengths
    denominator=(l1+l2+l3)*l1*l2*l3
    denominator>0 || return 0.0
    return clamp(4area2^2/denominator,0.0,1.0)
end

function _triangle_angle_shape(edge1,edge2,edge3)
    # Interior angles are evaluated with atan2 to retain accuracy near 0 and pi.
    angles=(atan(_norm3(_cross3(edge1,edge2)),_dot3(edge1,edge2)),
            atan(_norm3(_cross3(_scale3(edge1,-1.0),edge3)),
                 _dot3(_scale3(edge1,-1.0),edge3)),
            atan(_norm3(_cross3(_scale3(edge2,-1.0),_scale3(edge3,-1.0))),
                 _dot3(_scale3(edge2,-1.0),_scale3(edge3,-1.0))))
    any(iszero,angles) && return 0.0
    sharpness=500.0
    denominator=2atan(sharpness*(Float64(pi)/9))
    return minimum(angles) do angle
        offset=angle-Float64(pi)/3
        (atan(sharpness*(offset+Float64(pi)/9))+
         atan(sharpness*(Float64(pi)/9-offset)))/denominator
    end
end

function _tet_sicn(edges,determinant)
    determinant==0 && return 0.0
    edge1,edge2,edge3=edges
    jacobian1=edge1
    jacobian2=ntuple(i->(-edge1[i]+2edge2[i])/sqrt(3.0),3)
    jacobian3=ntuple(i->sqrt(1.5)*edge3[i]-
                        (3/(2sqrt(6.0)))*jacobian1[i]-
                        (1/(2sqrt(2.0)))*jacobian2[i],3)
    determinant_ideal=sqrt(2.0)*determinant
    norm_j=sum(abs2,jacobian1)+sum(abs2,jacobian2)+sum(abs2,jacobian3)
    cofactor1=_cross3(jacobian2,jacobian3)
    cofactor2=_cross3(jacobian3,jacobian1)
    cofactor3=_cross3(jacobian1,jacobian2)
    norm_cofactor=sum(abs2,cofactor1)+sum(abs2,cofactor2)+
                  sum(abs2,cofactor3)
    denominator=sqrt(norm_j*norm_cofactor)
    denominator>0 || return 0.0
    return clamp(3determinant_ideal/denominator,-1.0,1.0)
end

function _tet_sige(lengths,determinant)
    determinant==0 && return 0.0
    minimum(lengths)>0 || return 0.0
    l0,l1,l2,l3,l4,l5=lengths
    reciprocal_sum=
        inv(l0*l5*l1)+inv(l0*l5*l2)+inv(l0*l5*l3)+inv(l0*l5*l4)+
        inv(l1*l4*l0)+inv(l1*l4*l2)+inv(l1*l4*l3)+inv(l1*l4*l5)+
        inv(l2*l3*l0)+inv(l2*l3*l1)+inv(l2*l3*l4)+inv(l2*l3*l5)
    return clamp(sqrt(2.0)*determinant*reciprocal_sum/12,-1.0,1.0)
end

function _tet_isotropy(lengths,determinant)
    determinant>0 || return 0.0
    denominator=sum(abs2,lengths)
    denominator>0 || return 0.0
    return clamp(12*(abs(determinant)/2)^(2/3)/denominator,0.0,1.0)
end

function _tet_circumradius_normalized(edges,determinant)
    determinant==0 && return Inf
    edge1,edge2,edge3=edges
    numerator=ntuple(3) do i
        _dot3(edge1,edge1)*_cross3(edge2,edge3)[i]+
        _dot3(edge2,edge2)*_cross3(edge3,edge1)[i]+
        _dot3(edge3,edge3)*_cross3(edge1,edge2)[i]
    end
    radius=_norm3(_scale3(numerator,inv(2determinant)))
    return isfinite(radius) ? radius : Inf
end

function _segment_quality(name,length,valid)
    name in ("minEdge","maxEdge","volume") && return length
    name=="innerRadius" && return length/2
    name in ("gamma","outerRadius","minSICN") && return 0.0
    name=="angleShape" && return 1.0
    name=="minSJ" && return valid ? 1.0 : 0.0
    error("internal segment quality dispatch failure for $name")
end

function _triangle_quality(name,a,b,c)
    scale,lengths,area2,valid=_triangle_geometry(a,b,c)
    name in ("minDetJac","maxDetJac") && return 2triangle_area(a,b,c)
    name=="volume" && return triangle_area(a,b,c)
    name=="minSJ" && return valid ? 1.0 : 0.0
    l1,l2,l3=lengths
    ill_conditioned=valid &&
        (minimum(lengths)==0 || any(length->!isfinite(inv(length)),lengths) ||
         area2<=_FILTER_FACTOR*max(l1*l2,l1*l3,l2*l3))
    ill_conditioned && return _triangle_quality_big(name,a,b,c)
    name in ("minSICN","minIsotropy") &&
        return valid ? _triangle_isotropy(lengths,area2) : 0.0
    name=="minSIGE" && return valid ? _triangle_sige(lengths,area2) : 0.0
    name=="gamma" && return valid ? _triangle_gamma(lengths,area2) : 0.0
    if name=="innerRadius"
        perimeter=sum(lengths)
        return valid && perimeter>0 ? _rescale(scale,area2/perimeter) : 0.0
    end
    if name=="outerRadius"
        return valid ? _rescale(scale,prod(lengths)/(2area2)) : Inf
    end
    if name=="angleShape"
        edge1,edge2=_normalized_frame((a,b,c))[2]
        return valid ? _triangle_angle_shape(edge1,edge2,_sub3(edge2,edge1)) : 0.0
    end
    name=="minEdge" && return _rescale(scale,minimum(lengths))
    name=="maxEdge" && return _rescale(scale,maximum(lengths))
    error("internal triangle quality dispatch failure for $name")
end

function _tetrahedron_quality(name,a,b,c,d)
    scale,edges,lengths,determinant,face_area2,exact_fallback,valid=
        _tetrahedron_geometry(a,b,c,d)
    if name in ("minDetJac","maxDetJac")
        return 6tet_signed_volume(a,b,c,d)
    end
    name=="volume" && return tet_signed_volume(a,b,c,d)
    name=="minSJ" && return valid ? copysign(1.0,determinant) : 0.0
    name=="angleShape" && return 1.0
    ill_conditioned=valid &&
        (exact_fallback || minimum(lengths)==0 ||
         any(length->!isfinite(inv(length)),lengths))
    ill_conditioned && return _tet_quality_big(name,a,b,c,d)
    name=="minSICN" && return valid ? _tet_sicn(edges,determinant) : 0.0
    name=="minSIGE" && return valid ? _tet_sige(lengths,determinant) : 0.0
    name=="minIsotropy" && return _tet_isotropy(lengths,determinant)
    surface_area2=sum(face_area2)
    inradius=valid && surface_area2>0 ? determinant/surface_area2 : 0.0
    circumradius=valid ? _tet_circumradius_normalized(edges,determinant) : Inf
    if name=="gamma"
        isfinite(circumradius) && circumradius>0 || return 0.0
        return clamp(3abs(inradius)/circumradius,0.0,1.0)
    end
    name=="innerRadius" && return _rescale(scale,inradius)
    name=="outerRadius" && return _rescale(scale,circumradius)
    name=="minEdge" && return _rescale(scale,minimum(lengths))
    name=="maxEdge" && return _rescale(scale,maximum(lengths))
    error("internal tetrahedron quality dispatch failure for $name")
end

@inline function _element_kind(tag::Int,triangle_offset::Int,
                               tetrahedron_offset::Int)
    tag<=triangle_offset && return 1
    tag<=tetrahedron_offset && return 2
    return 3
end

function _validate_quality_kinds(tags,name,triangle_offset,tetrahedron_offset,
                                 caller::AbstractString)
    name in _UNDEFINED_SEGMENT_QUALITIES || return nothing
    for tag in tags
        if _element_kind(tag,triangle_offset,tetrahedron_offset)==1
            throw(ArgumentError(
                "$caller: quality $(repr(name)) is not defined for cached " *
                "linear segments; request triangle or tetrahedron tags"))
        end
    end
    return nothing
end

function _element_quality(mesh::Mesh,tag::Int,name::String,
                          triangle_offset::Int,tetrahedron_offset::Int)
    kind=_element_kind(tag,triangle_offset,tetrahedron_offset)
    if kind==1
        a=node(mesh,mesh.segs[1,tag]);b=node(mesh,mesh.segs[2,tag])
        length,valid=_segment_geometry(a,b)
        return _segment_quality(name,length,valid)
    elseif kind==2
        cell=tag-triangle_offset
        a=node(mesh,mesh.tris[1,cell]);b=node(mesh,mesh.tris[2,cell])
        c=node(mesh,mesh.tris[3,cell])
        return _triangle_quality(name,a,b,c)
    end
    cell=tag-tetrahedron_offset
    a=node(mesh,mesh.tets[1,cell]);b=node(mesh,mesh.tets[2,cell])
    c=node(mesh,mesh.tets[3,cell]);d=node(mesh,mesh.tets[4,cell])
    return _tetrahedron_quality(name,a,b,c,d)
end

"""
    mesh_element_qualities(mesh, element_tags, quality_name="minSICN")

Return one `Float64` quality per dense element tag, preserving request order and
duplicates. Supported names match the documented Gmsh 4.15.2 quality list.
`minDetJac`, `maxDetJac`, `minSIGE`, and `minIsotropy` reject segment tags because
Gmsh 4.15.2 has no reliable 1-D definition for those measures. Degenerate simplex
shape measures return zero; an undefined circumradius returns `Inf`.
"""
function mesh_element_qualities(mesh::Mesh,element_tags,
                                quality_name="minSICN")
    caller="mesh_element_qualities"
    name=_quality_name(quality_name,caller)
    tags=_quality_tags(mesh,element_tags,caller)
    triangle_offset,tetrahedron_offset,_=mesh_element_offsets(mesh)
    _validate_quality_kinds(
        tags,name,triangle_offset,tetrahedron_offset,caller)
    qualities=Vector{Float64}(undef,length(tags))
    @inbounds for index in eachindex(tags)
        qualities[index]=_element_quality(
            mesh,tags[index],name,triangle_offset,tetrahedron_offset)
    end
    return qualities
end

"""Return one named quality for a dense linear-simplex element tag."""
function mesh_element_quality(mesh::Mesh,element_tag,
                              quality_name="minSICN")
    return only(mesh_element_qualities(mesh,(element_tag,),quality_name))
end

end # module MeshElementQuality
