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
using ..GmshLibm: _gm_atan, _gm_acos, _gm_atan2
using LinearAlgebra: svd

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

# MQuadrangle::getVolume uses the physical six-distance formula. For a
# warped quadrangle this differs from a bilinear reference-map integral.
@inline function _quadrangle_native_distance(a,b)
    x=a[1]-b[1];y=a[2]-b[2];z=a[3]-b[3]
    return sqrt(x*x+y*y+z*z)
end
function _quadrangle_native_volume(a,b,c,d)
    ab=_quadrangle_native_distance(a,b);bc=_quadrangle_native_distance(b,c)
    cd=_quadrangle_native_distance(c,d);da=_quadrangle_native_distance(d,a)
    ac=_quadrangle_native_distance(a,c);bd=_quadrangle_native_distance(b,d)
    mn=2.0*ac*bd
    abcd=ab*ab-bc*bc+cd*cd-da*da
    squared=mn*mn-abcd*abcd
    return squared<0 ? NaN : sqrt(squared)/4.0
end

@inline _quadrangle_native_dot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _quadrangle_native_norm(a)=sqrt(_quadrangle_native_dot(a,a))
@inline function _quadrangle_native_unit(a)
    magnitude=_quadrangle_native_norm(a)
    magnitude==0.0 && return a
    inverse=1.0/magnitude
    return (a[1]*inverse,a[2]*inverse,a[3]*inverse)
end
@inline function _quadrangle_native_angle(a,b,c)
    first=_sub3(a,b);second=_sub3(c,b)
    return _gm_atan2(_quadrangle_native_norm(_cross3(first,second)),
        _quadrangle_native_dot(first,second))
end
function _quadrangle_native_eta(points)
    edges=ntuple(i->_sub3(points[mod1(i+1,4)],points[i]),4)
    normals=ntuple(i->_cross3(edges[i],edges[mod1(i+1,4)]),4)
    sign=1.0
    for i in 2:4
        _quadrangle_native_dot(normals[1],normals[i])<0 && (sign=-1.0)
    end
    deviation=0.0
    for i in 1:4
        angle=min(180.0,180.0*_quadrangle_native_angle(
            points[i],points[mod1(i+1,4)],points[mod1(i+2,4)])/pi)
        deviation=i==1 ? abs(90.0-angle) : max(abs(90.0-angle),deviation)
    end
    return sign*(1.0-deviation*(1.0/90.0))
end
function _quadrangle_native_gradients(points,u,v)
    du=(-0.25+0.25*v,0.25-0.25*v,0.25+0.25*v,-0.25-0.25*v)
    dv=(-0.25+0.25*u,-0.25-0.25*u,0.25+0.25*u,0.25-0.25*u)
    first=ntuple(axis->begin
        total=0.0
        for i in 1:4
            total+=points[i][axis]*du[i]
        end
        total
    end,3)
    second=ntuple(axis->begin
        total=0.0
        for i in 1:4
            total+=points[i][axis]*dv[i]
        end
        total
    end,3)
    return first,second
end
function _quadrangle_native_angle_shape(points)
    coefficient=100.0
    denominator=_gm_atan(coefficient*(pi/4.0))+
        _gm_atan(coefficient*(2.0*pi/4.0-(pi/4.0)))
    quality=floatmax(Float64)
    for (u,v) in ((-1.,-1.),(-1.,1.),(1.,1.),(1.,-1.),
            (0.,-1.),(0.,1.),(1.,0.),(-1.,0.),(0.,0.))
        first,second=_quadrangle_native_gradients(points,u,v)
        cosine=_quadrangle_native_dot(_quadrangle_native_unit(first),
            _quadrangle_native_unit(second))
        angle=abs(_gm_acos(cosine))-pi/2.0
        value=(_gm_atan(coefficient*(angle+pi/4.0))+
            _gm_atan(coefficient*(2.0*pi/4.0-(angle+pi/4.0))))/denominator
        # std::min retains its first argument when the second is NaN.
        value<quality && (quality=value)
    end
    return quality
end
function _quadrangle_native_tangent_radius(x,y,i)
    a=ntuple(j->y[mod1(i+j-1,4)]-y[mod1(i+j,4)],3)
    b=ntuple(j->x[mod1(i+j,4)]-x[mod1(i+j-1,4)],3)
    c=ntuple(j->y[mod1(i+j,4)]*x[mod1(i+j-1,4)]-
        y[mod1(i+j-1,4)]*x[mod1(i+j,4)],3)
    lengths=ntuple(j->sqrt(a[j]*a[j]+b[j]*b[j]),3)
    a12=a[1]/lengths[1]-a[2]/lengths[2]
    a23=a[2]/lengths[2]-a[3]/lengths[3]
    b12=b[1]/lengths[1]-b[2]/lengths[2]
    b23=b[2]/lengths[2]-b[3]/lengths[3]
    c12=c[1]/lengths[1]-c[2]/lengths[2]
    c23=c[2]/lengths[2]-c[3]/lengths[3]
    center_x=(c12*b23-c23*b12)/(a23*b12-a12*b23)
    center_y=b12!=0.0 ? -a12/b12*center_x-c12/b12 :
        -a23/b23*center_x-c23/b23
    return (a[1]*center_x+b[1]*center_y+c[1])/lengths[1]
end
function _quadrangle_native_inner_radius(points)
    center=ntuple(axis->(points[1][axis]+points[2][axis]+
        points[3][axis]+points[4][axis])/4.0,3)
    matrix=Matrix{Float64}(undef,4,3)
    for i in 1:4,axis in 1:3
        matrix[i,axis]=points[i][axis]-center[axis]
    end
    decomposition=svd(matrix)
    sigma=decomposition.S
    column=abs(sigma[1])<abs(sigma[2]) && abs(sigma[1])<abs(sigma[3]) ? 1 :
        abs(sigma[2])<abs(sigma[1]) && abs(sigma[2])<abs(sigma[3]) ? 2 : 3
    a,b,c=decomposition.V[:,column]
    d=-(center[1]*a+center[2]*b+center[3]*c)
    magnitude=sqrt(a*a+b*b+c*c)
    projected=ntuple(i->begin
        px,py,pz=points[i]
        (((b*b+c*c)*px-a*b*py-a*c*pz-d*a)/magnitude,
         (-a*b*px+(a*a+c*c)*py-b*c*pz-d*b)/magnitude,
         (-a*c*px-b*c*py+(a*a+b*b)*pz-d*c)/magnitude)
    end,4)
    edge_first=_sub3(projected[2],projected[1])
    edge_second=_sub3(projected[3],projected[1])
    normal=_cross3(edge_first,edge_second)
    axis_first=_quadrangle_native_unit(edge_first)
    axis_second=_quadrangle_native_unit(_cross3(normal,edge_first))
    x=ntuple(i->_quadrangle_native_dot(_sub3(projected[i],projected[1]),axis_first),4)
    y=ntuple(i->_quadrangle_native_dot(_sub3(projected[i],projected[1]),axis_second),4)
    radius=1e22
    for i in 1:4
        value=_quadrangle_native_tangent_radius(x,y,i)
        value<radius && (radius=value)
    end
    return radius
end

@inline function _quadrangle_native_determinant(first,second,normal)
    return first[1]*(second[2]*normal[3]-second[3]*normal[2])-
        second[1]*(first[2]*normal[3]-first[3]*normal[2])+
        normal[1]*(first[2]*second[3]-first[3]*second[2])
end
function _quadrangle_native_normal(points)
    first,second=_quadrangle_native_gradients(points,0.,0.)
    normal=_cross3(first,second)
    magnitude=_quadrangle_native_norm(normal)
    inverse=1.0/magnitude
    return ntuple(i->normal[i]*inverse,3),magnitude
end
function _quadrangle_native_sicn(first,second,normal)
    x,y,z=first;u,v,w=second
    xs=x*x;ys=y*y;zs=z*z;us=u*u;vs=v*v;ws=w*w
    dx=xs-us;dy=ys-vs;cx=x*u;cy=y*v
    s1=ws*ws;s2=(zs-dy-dx)*ws;s3=(cy+cx)*z*w
    s4=zs*zs;s5=(dy+dx)*zs;s6=cx*cy;s7=ys*ys
    s8=dy*ys;s9=xs*xs;s10=dx*xs;s11=dy*dy;s12=dx*dy;s13=dx*dx
    discriminant=2.0*(s2+s5+s12)+4.0*(s7-s8+s9-s10)+
        8.0*(s3+s6)+s1+s4+s11+s13
    total=xs+us+ys+vs+zs+ws
    root=discriminant>0 ? sqrt(discriminant) : 0.0
    largest=0.5*(total+root);smallest=0.5*(total-root)
    squared=largest*smallest
    value=2.0*(squared<0 ? NaN : sqrt(squared))/(largest+smallest)
    local_normal=_cross3(first,second)
    return _quadrangle_native_dot(local_normal,normal)>=0 ? value : -value
end
function _quadrangle_native_jacobian_quality(name,points)
    normal,scale=_quadrangle_native_normal(points)
    if name=="minSJ" && scale==0.0
        return 0.0
    end
    # MFace::normal supplies the validity sign for CondNumBasis. It uses the
    # first three slots; JacobianBasis instead uses the barycenter normal.
    sicn_normal=_quadrangle_native_unit(_cross3(_sub3(points[2],points[1]),
        _sub3(points[3],points[1])))
    minimum=floatmax(Float64);maximum=-minimum
    for (u,v) in ((-1.,-1.),(1.,-1.),(1.,1.),(-1.,1.))
        first,second=_quadrangle_native_gradients(points,u,v)
        determinant=_quadrangle_native_determinant(first,second,normal)
        value=name=="minSICN" ? _quadrangle_native_sicn(first,second,sicn_normal) :
            name=="minSIGE" ? determinant/_quadrangle_native_norm(first)/
                _quadrangle_native_norm(second) :
            name=="minSJ" ? determinant/scale : determinant
        # The native range scans use std::min/max, which retain the existing
        # extrema when a subsequent sample is NaN. CondNumBasis uses the first
        # sampling value as its min_element seed, retaining a first NaN.
        if name=="minSICN" && u==-1.0 && v==-1.0
            minimum=value
        else
            value<minimum && (minimum=value)
        end
        value>maximum && (maximum=value)
    end
    return name=="maxDetJac" ? maximum : minimum
end

function _quadrangle_native_icn_domain(points,normal,bounds)
    u0,u1,v0,v1=bounds
    gradients=Array{Float64}(undef,2,2,6)
    minimum=floatmax(Float64)
    for j in 1:2,i in 1:2
        first,second=_quadrangle_native_gradients(points,i==1 ? u0 : u1,j==1 ? v0 : v1)
        determinant=_quadrangle_native_determinant(first,second,normal)
        total=0.0
        for (k,value) in enumerate((first...,second...))
            gradients[i,j,k]=value
            total+=value*value
        end
        value=2.0*determinant/total
        value<minimum && (minimum=value)
    end
    lower=0.0
    if !(minimum-lower<1e-3)
        lower=floatmax(Float64)
        for j in 0:2,i in 0:2
            u=u0+(u1-u0)*(i/2.0);v=v0+(v1-v0)*(j/2.0)
            first,second=_quadrangle_native_gradients(points,u,v)
            numerator=_quadrangle_native_determinant(first,second,normal)
            if numerator<0
                lower=0.0;break
            end
            denominator=0.0
            for k in 1:6,b in 0:1,a in 0:1
                c=i-a;d=j-b
                0<=c<=1 && 0<=d<=1 || continue
                weight=1.0/((i==1 ? 2 : 1)*(j==1 ? 2 : 1))
                denominator+=gradients[a+1,b+1,k]*gradients[c+1,d+1,k]*weight
            end
            denominator<=0.0 && continue
            ratio=2.0*numerator/denominator
            ratio<lower && (lower=ratio)
        end
    end
    return (bounds=bounds,minimum=minimum,lower=lower)
end
function _quadrangle_native_isotropy(points)
    normal,_=_quadrangle_native_normal(points)
    minimum_jac=_quadrangle_native_jacobian_quality("minDetJac",points)
    maximum_jac=_quadrangle_native_jacobian_quality("maxDetJac",points)
    (minimum_jac<=0 && maximum_jac>=0) && return 0.0
    # Unknown validity with an undefined primary normal yields no positive
    # bound in the pinned ICN subdivision routine.
    all(isfinite,normal) || return 0.0
    domains=[_quadrangle_native_icn_domain(points,normal,(-1.,1.,-1.,1.))]
    minimum=domains[1].minimum
    for _ in 1:999
        index=1
        for i in 2:length(domains)
            domains[i].lower<domains[index].lower && (index=i)
        end
        current=domains[index]
        minimum-current.lower<1e-3+(1e-2-1e-3)*max(current.lower,0.0) && break
        deleteat!(domains,index)
        u0,u1,v0,v1=current.bounds
        um=0.5*(u0+u1);vm=0.5*(v0+v1)
        for bounds in ((u0,um,v0,vm),(um,u1,v0,vm),(u0,um,vm,v1),(um,u1,vm,v1))
            next=_quadrangle_native_icn_domain(points,normal,bounds)
            next.minimum<minimum && (minimum=next.minimum)
            push!(domains,next)
        end
    end
    lower=minimum
    for domain in domains
        domain.lower<lower && (lower=domain.lower)
    end
    factor=0.5*(lower+minimum)
    return factor*minimum+(1.0-factor)*lower
end

function _quality_name(value,caller::AbstractString;allow_quad_eta=false)
    value isa AbstractString || throw(ArgumentError(
        "$caller: quality_name must be a string"))
    name=String(value)
    occursin('\0',name) && throw(ArgumentError(
        "$caller: quality_name must not contain NUL"))
    (name in ELEMENT_QUALITY_NAMES || (allow_quad_eta && name=="eta")) || throw(ArgumentError(
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
