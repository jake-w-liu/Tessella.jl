# Native N>2 public measures; the standalone FEM quality evaluator is unchanged.
@inline _native_cell_dot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _native_cell_cross(a,b)=
    (a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
@inline _native_cell_norm(a)=sqrt(_native_cell_dot(a,a))
@inline _native_cell_difference(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
function _native_cell_unit(a)
    norm=_native_cell_norm(a)
    factor=norm==0.0 ? 1.0 : 1.0/norm
    return (a[1]*factor,a[2]*factor,a[3]*factor)
end

function _native_cell_triangle_gamma(a,b,c)
    side_a=_native_cell_unit(_native_cell_difference(c,b))
    side_b=_native_cell_unit(_native_cell_difference(a,c))
    side_c=_native_cell_unit(_native_cell_difference(b,a))
    sina=_native_cell_norm(_native_cell_cross(side_b,side_c))
    sinb=_native_cell_norm(_native_cell_cross(side_c,side_a))
    sinc=_native_cell_norm(_native_cell_cross(side_a,side_b))
    return sina==0.0 && sinb==0.0 && sinc==0.0 ? 0.0 :
        2*(2*sina*sinb*sinc/(sina+sinb+sinc))
end

function _native_cell_quadrangle_gamma(a,b,c,d)
    v01=_native_cell_difference(b,a);v12=_native_cell_difference(c,b)
    v23=_native_cell_difference(d,c);v30=_native_cell_difference(a,d)
    ca=_native_cell_cross(v01,v12);cb=_native_cell_cross(v12,v23)
    cc=_native_cell_cross(v23,v30);cd=_native_cell_cross(v30,v01)
    sign=_native_cell_dot(ca,cb)<0 || _native_cell_dot(ca,cc)<0 || _native_cell_dot(ca,cd)<0 ? -1.0 : 1.0
    corners=(a,b,c,d);worst=0.0
    for index in 1:4
        center=corners[mod1(index+1,4)]
        first=_native_cell_difference(corners[index],center)
        last=_native_cell_difference(corners[mod1(index+2,4)],center)
        angle=180*atan(_native_cell_norm(_native_cell_cross(first,last)),_native_cell_dot(first,last))/Float64(pi)
        angle=min(180.0,angle)
        worst=max(abs(90.0-angle),worst)
    end
    return sign*(1.0-worst*(1.0/90.0))
end

function _native_cell_orient3_value(a,b,c,d)
    ad=_native_cell_difference(a,d);bd=_native_cell_difference(b,d);cd=_native_cell_difference(c,d)
    determinant=ad[3]*(bd[1]*cd[2]-cd[1]*bd[2])+
        bd[3]*(cd[1]*ad[2]-ad[1]*cd[2])+cd[3]*(ad[1]*bd[2]-bd[1]*ad[2])
    permanent=(abs(bd[1]*cd[2])+abs(cd[1]*bd[2]))*abs(ad[3])+
        (abs(cd[1]*ad[2])+abs(ad[1]*cd[2]))*abs(bd[3])+
        (abs(ad[1]*bd[2])+abs(bd[1]*ad[2]))*abs(cd[3])
    abs(determinant)>O3DERRBOUND_A*permanent && return determinant
    # Mesh quality needs the determinant value, not just the predicate sign.
    ar=ntuple(i->Rational{BigInt}(a[i])-Rational{BigInt}(d[i]),3)
    br=ntuple(i->Rational{BigInt}(b[i])-Rational{BigInt}(d[i]),3)
    cr=ntuple(i->Rational{BigInt}(c[i])-Rational{BigInt}(d[i]),3)
    return Float64(ar[3]*(br[1]*cr[2]-cr[1]*br[2])+
        br[3]*(cr[1]*ar[2]-ar[1]*cr[2])+cr[3]*(ar[1]*br[2]-br[1]*ar[2]))
end

@inline _native_cell_triangle_area(a,b,c)=
    0.5*_native_cell_norm(_native_cell_cross(_native_cell_difference(c,b),_native_cell_difference(a,b)))

function _native_cell_tetrahedron_gamma(a,b,c,d)
    volume=_native_cell_orient3_value(a,b,c,d)/6.0
    abs(volume)==0.0 && return 0.0
    ab=_native_cell_difference(b,a);ac=_native_cell_difference(c,a);ad=_native_cell_difference(d,a)
    dc=_native_cell_difference(d,c);db=_native_cell_difference(d,b);cb=_native_cell_difference(c,b)
    la=_native_cell_dot(ab,ab);lb=_native_cell_dot(ac,ac);lc=_native_cell_dot(ad,ad)
    lA=_native_cell_dot(dc,dc);lB=_native_cell_dot(db,db);lC=_native_cell_dot(cb,cb)
    lalA=sqrt(la*lA);lblB=sqrt(lb*lB);lclC=sqrt(lc*lC)
    inside=(lalA+lblB+lclC)*(lalA+lblB-lclC)*(lalA-lblB+lclC)*(-lalA+lblB+lclC)
    inside<=0.0 && return 0.0
    partR=sqrt(inside)*(1.0/24.0)
    s1=abs(_native_cell_triangle_area(a,b,c));s2=abs(_native_cell_triangle_area(a,c,d))
    s3=abs(_native_cell_triangle_area(a,b,d));s4=abs(_native_cell_triangle_area(b,c,d))
    rho=3*3*volume/(s1+s2+s3+s4)
    return rho*volume/partR
end

function _native_cell_tetrahedron_volume(a,b,c,d)
    # MTetrahedron::getVolume remains the corner determinant for every order.
    # Its MElement quadrature implementation is not the virtual override.
    e1=_native_cell_difference(b,a)
    e2=_native_cell_difference(c,a)
    e3=_native_cell_difference(d,a)
    determinant=e1[1]*(e2[2]*e3[3]-e3[2]*e2[3])-
        e2[1]*(e1[2]*e3[3]-e3[2]*e1[3])+
        e3[1]*(e1[2]*e2[3]-e2[2]*e1[3])
    return determinant*0.166666666666666666
end

function _native_cell_prism_ncj(a,b,c,d)
    v1=_native_cell_difference(b,a);v2=_native_cell_difference(c,a);v3=_native_cell_difference(d,a)
    return (2.0/sqrt(3.0))*abs(_native_cell_dot(v1,_native_cell_cross(v2,v3)))/
        (_native_cell_norm(v1)*_native_cell_norm(v2)*_native_cell_norm(v3))
end

function _native_cell_quality_gamma(mesh,cells,cell,family,caller)
    a=_mixed_node(mesh,cells,cell,1);b=_mixed_node(mesh,cells,cell,2);c=_mixed_node(mesh,cells,cell,3)
    family===:tri && return _native_cell_triangle_gamma(a,b,c)
    d=_mixed_node(mesh,cells,cell,4)
    family===:qua && return _native_cell_quadrangle_gamma(a,b,c,d)
    family===:tet && return _native_cell_tetrahedron_gamma(a,b,c,d)
    family in (:hex,:pyr) && return 0.0
    if family===:pri
        e=_mixed_node(mesh,cells,cell,5);f=_mixed_node(mesh,cells,cell,6)
        values=(_native_cell_prism_ncj(a,b,c,d),_native_cell_prism_ncj(b,a,c,e),
            _native_cell_prism_ncj(c,a,b,f),_native_cell_prism_ncj(d,a,e,f),
            _native_cell_prism_ncj(e,b,d,f),_native_cell_prism_ncj(f,c,d,e))
        # qmPrism::minNCJ uses std::min_element: later NaNs do not
        # replace the first candidate, while a first-candidate NaN survives.
        result=first(values)
        for i in 2:length(values)
            values[i]<result && (result=values[i])
        end
        return result
    end
    throw(ArgumentError("$caller: unsupported native gamma family $family"))
end

struct _NativeCellMeasureRule
    weights::Vector{Float64}
    shapes::Vector{Tuple{Vector{Float64},Vector{NTuple{3,Float64}}}}
    dimension::Int
end

function _native_cell_measure_rule(msh,caller)
    spec=msh_spec(msh)
    _native_cell_numeric_contract(spec,caller)
    coordinates,weights=mesh_integration_points(msh,"Gauss$(spec.order+3)";caller=caller)
    data=spec.family in (:tri,:qua,:tet,:hex,:pri) ? _native_cell_polynomial(msh,caller) : nothing
    shapes=Tuple{Vector{Float64},Vector{NTuple{3,Float64}}}[]
    sizehint!(shapes,length(weights))
    for point in eachindex(weights)
        offset=3point-2
        uvw=(coordinates[offset],coordinates[offset+1],coordinates[offset+2])
        push!(shapes,_native_cell_query_shape(msh,uvw,data,caller))
    end
    return _NativeCellMeasureRule(weights,shapes,spec.dim)
end

function _native_cell_quality_volume(mesh,cells,cell,rule,caller,tag)
    result=0.0
    for point in eachindex(rule.weights)
        _,determinant,_=_native_cell_public_frame(mesh,cells,cell,rule.dimension,
            rule.shapes[point],caller,tag)
        result+=determinant*rule.weights[point]
    end
    return result
end

# Order-agnostic MElement measures that upstream defines on the primary
# vertices of hexahedra and prisms.
function _native_cell_hexahedron_angles(mesh,cells,cell)
    # qmHexahedron::angles: smallest corner angle over the six quad faces.
    angle_max=0.0;angle_min=pi
    for face in MeshEntityTopology._HEXAHEDRON_FACES
        corners=ntuple(i->_mixed_node(mesh,cells,cell,face[i]),4)
        for j in 1:4
            a=_native_cell_difference(corners[mod1(j+2,4)],corners[mod1(j+1,4)])
            b=_native_cell_difference(corners[mod1(j+1,4)],corners[j])
            angle=acos(_native_cell_dot(a,b)/(_native_cell_norm(a)*_native_cell_norm(b)))
            angle>angle_max && (angle_max=angle)
            angle<angle_min && (angle_min=angle)
        end
    end
    return 1.0-max((angle_max-0.5*pi)/(0.5*pi),(0.5*pi-angle_min)/(0.5*pi))
end

function _native_cell_hexahedron_inner_radius(mesh,cells,cell)
    # MHexahedron::getInnerRadius: the smallest face tangent radius.
    radius=Inf
    for face in MeshEntityTopology._HEXAHEDRON_FACES
        points=ntuple(i->_mixed_node(mesh,cells,cell,face[i]),4)
        value=MeshElementQuality._quadrangle_native_inner_radius(points)
        value<radius && (radius=value)
    end
    return radius
end

function _native_cell_prism_inner_radius(mesh,cells,cell)
    # MPrism::getInnerRadius: min(inradius of the bottom triangular edges
    # {0,1,3}, half the (0,3) vertical edge length).
    a=_mixed_node(mesh,cells,cell,1);b=_mixed_node(mesh,cells,cell,2)
    c=_mixed_node(mesh,cells,cell,3);d=_mixed_node(mesh,cells,cell,4)
    dist=(_native_cell_norm(_native_cell_difference(b,a)),
          _native_cell_norm(_native_cell_difference(c,a)),
          _native_cell_norm(_native_cell_difference(c,b)))
    k=0.5*(dist[1]+dist[2]+dist[3])
    rad_tri=sqrt(k*(k-dist[1])*(k-dist[2])*(k-dist[3]))/k
    rad_vert=0.5*_native_cell_norm(_native_cell_difference(d,a))
    return min(rad_tri,rad_vert)
end
