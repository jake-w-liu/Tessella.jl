# Gmsh 4.15.2's simplex quality extrema use Bernstein bounds on the actual
# affine derivative field, followed by midpoint subdivision. Corner values
# guide convergence; the reported determinant extrema remain polynomial bounds.
struct _LegacyP2QualityDomain{V}
    vertices::V
    min_l::Float64
    max_l::Float64
    min_b::Float64
    max_b::Float64
end

@inline _legacy_p2_frame_midpoint(a,b)=ntuple(column->ntuple(axis->
    (a[column][axis]+b[column][axis])/2,3),3)

function _legacy_p2_quality_children(v::NTuple{3})
    a,b,c=v
    ab=_legacy_p2_frame_midpoint(a,b)
    ac=_legacy_p2_frame_midpoint(a,c)
    bc=_legacy_p2_frame_midpoint(b,c)
    return ((a,ab,ac),(bc,ab,ac),(bc,c,ac),(bc,ab,b))
end

function _legacy_p2_quality_children(v::NTuple{4})
    a,b,c,d=v
    ab=_legacy_p2_frame_midpoint(a,b);ac=_legacy_p2_frame_midpoint(a,c)
    ad=_legacy_p2_frame_midpoint(a,d);bc=_legacy_p2_frame_midpoint(b,c)
    bd=_legacy_p2_frame_midpoint(b,d);cd=_legacy_p2_frame_midpoint(c,d)
    # The middle octahedron uses the AD--BC diagonal, matching bezierBasis.cpp.
    return ((a,ab,ac,ad),(bc,ab,ac,ad),(bc,cd,ac,ad),(bc,ab,bd,ad),
            (bc,cd,bd,ad),(bc,cd,ac,c),(bc,ab,bd,b),(d,cd,bd,ad))
end

@inline function _legacy_p2_corner_det(frame,normal)
    du,dv,dw=frame
    return normal===nothing ? MeshElementQuality._det3(du,dv,dw) :
        MeshElementQuality._dot3(normal,MeshElementQuality._cross3(du,dv))
end

function _legacy_p2_jac_domain(vertices,normal)
    coefficients=_legacy_p2_det_coefficients(vertices,normal)
    corners=ntuple(vertex->_legacy_p2_corner_det(vertices[vertex],normal),length(vertices))
    return _LegacyP2QualityDomain(vertices,minimum(corners),maximum(corners),
        minimum(coefficients),maximum(coefficients))
end

@inline function _legacy_p2_ideal_frame(frame,dimension)
    du,dv,dw=frame
    ideal_v=ntuple(axis->(2dv[axis]-du[axis])/sqrt(3.0),3)
    ideal_w=dimension==2 ? (0.0,0.0,0.0) : ntuple(axis->
        sqrt(1.5)*dw[axis]-(du[axis]+dv[axis])/sqrt(6.0),3)
    return (du,ideal_v,ideal_w)
end

@inline _legacy_p2_frobenius_dot(a,b)=sum(
    a[column][axis]*b[column][axis] for column in 1:3 for axis in 1:3)

function _legacy_p2_icn_denominator(vertices::NTuple{3})
    coefficients=zeros(Float64,6)
    for first in 1:3,second in 1:3
        group=_LEGACY_P2_TRI_DET_GROUP[first][second]
        coefficients[group]+=_legacy_p2_frobenius_dot(vertices[first],vertices[second])
    end
    for index in eachindex(coefficients)
        coefficients[index]/=HighOrder._B2TRI_WEIGHT[index]
    end
    return coefficients
end

function _legacy_p2_icn_denominator(vertices::NTuple{4})
    # The norm of the affine matrix field is bounded above by the affine
    # combination of corner norms. Cube that positive degree-one polynomial.
    norms=ntuple(vertex->sqrt(_legacy_p2_frobenius_dot(vertices[vertex],vertices[vertex])),4)
    coefficients=zeros(Float64,20)
    for first in 1:4,second in 1:4,third in 1:4
        coefficients[HighOrder._B3_GROUP[first,second,third]]+=
            norms[first]*norms[second]*norms[third]
    end
    for index in eachindex(coefficients)
        coefficients[index]/=HighOrder._B3_NPERM[index]
    end
    return coefficients
end

function _legacy_p2_rational_lower(numerator,denominator)
    lower=-Inf
    upper=Inf
    for index in eachindex(numerator)
        den=denominator[index];num=numerator[index]
        if den==0
            num<0 && return -Inf
        elseif den>0
            upper=min(upper,num/den)
        else
            lower=max(lower,num/den)
        end
    end
    return lower>upper ? -Inf : upper
end

@inline function _legacy_p2_icn_corner(frame,normal,dimension)
    determinant=_legacy_p2_corner_det(frame,normal)
    norm_squared=_legacy_p2_frobenius_dot(frame,frame)
    norm_squared>0 || return 0.0
    return dimension==2 ? 2determinant/norm_squared :
        3cbrt(determinant*determinant)/norm_squared
end

@inline _legacy_p2_icn_bounds_ok(domain,min_l)=
    min_l-domain.min_b < .001+.009*max(domain.min_b,0.0)

function _legacy_p2_icn_domain(vertices,normal)
    dimension=length(vertices)-1
    ideal=ntuple(vertex->_legacy_p2_ideal_frame(vertices[vertex],dimension),length(vertices))
    corners=ntuple(vertex->_legacy_p2_icn_corner(ideal[vertex],normal,dimension),length(vertices))
    min_l=minimum(corners);max_l=maximum(corners)
    initial=_LegacyP2QualityDomain(vertices,min_l,max_l,0.0,0.0)
    _legacy_p2_icn_bounds_ok(initial,min_l) && return initial
    numerator=_legacy_p2_det_coefficients(ideal,normal)
    min_b=if any(<(0.0),numerator)
        0.0
    else
        denominator=_legacy_p2_icn_denominator(ideal)
        bound=_legacy_p2_rational_lower(numerator,denominator)
        dimension==2 ? 2bound : 3cbrt(bound*bound)
    end
    return _LegacyP2QualityDomain(vertices,min_l,max_l,min_b,0.0)
end

@inline function _legacy_p2_jac_bounds_ok(domain,min_l,max_l)
    tolerance=max(abs(min_l),abs(max_l))*.001
    return (min_l<=0 || domain.min_b>0) &&
           (max_l>=0 || domain.max_b<0) &&
           min_l-domain.min_b<tolerance && domain.max_b-max_l<tolerance
end

@inline _legacy_p2_quality_priority(domain,maximum_mode)=
    maximum_mode ? -domain.max_b : domain.min_b

function _legacy_p2_heap_sift!(domains,index,maximum_mode)
    count=length(domains)
    while 2index<=count
        child=2index
        if child<count && _legacy_p2_quality_priority(domains[child+1],maximum_mode)<
                          _legacy_p2_quality_priority(domains[child],maximum_mode)
            child+=1
        end
        _legacy_p2_quality_priority(domains[index],maximum_mode)<=
            _legacy_p2_quality_priority(domains[child],maximum_mode) && break
        domains[index],domains[child]=domains[child],domains[index]
        index=child
    end
    return nothing
end

function _legacy_p2_heap_push!(domains,domain,maximum_mode)
    push!(domains,domain)
    index=length(domains)
    while index>1
        parent=index÷2
        _legacy_p2_quality_priority(domains[parent],maximum_mode)<=
            _legacy_p2_quality_priority(domains[index],maximum_mode) && break
        domains[parent],domains[index]=domains[index],domains[parent]
        index=parent
    end
    return nothing
end

function _legacy_p2_subdivide_quality!(domains,normal,min_l,max_l,
        isotropy::Bool,maximum_mode::Bool)
    for index in (length(domains)÷2):-1:1
        _legacy_p2_heap_sift!(domains,index,maximum_mode)
    end
    for iteration in 1:999
        domain=first(domains)
        ready=isotropy ? _legacy_p2_icn_bounds_ok(domain,min_l) :
            _legacy_p2_jac_bounds_ok(domain,min_l,max_l)
        ready && break
        tail=pop!(domains)
        if !isempty(domains)
            domains[1]=tail
            _legacy_p2_heap_sift!(domains,1,maximum_mode)
        end
        for vertices in _legacy_p2_quality_children(domain.vertices)
            child=isotropy ? _legacy_p2_icn_domain(vertices,normal) :
                _legacy_p2_jac_domain(vertices,normal)
            min_l=min(min_l,child.min_l);max_l=max(max_l,child.max_l)
            _legacy_p2_heap_push!(domains,child,maximum_mode)
        end
    end
    return min_l,max_l
end

function _legacy_p2_jac_extrema(vertices,normal)
    initial=_legacy_p2_jac_domain(vertices,normal)
    # A zero polynomial has exact extrema without subdivision. This also
    # avoids a zero tolerance spending the entire upstream subdivision budget.
    initial.min_b==initial.max_b==0 && return 0.0,0.0
    domains=[initial]
    min_l,max_l=_legacy_p2_subdivide_quality!(domains,normal,
        initial.min_l,initial.max_l,false,false)
    _legacy_p2_subdivide_quality!(domains,normal,min_l,max_l,false,true)
    return minimum(domain.min_b for domain in domains),
           maximum(domain.max_b for domain in domains)
end

function _legacy_p2_bounded_quality(mesh,vertices,normal,scale,name,caller)
    min_det,max_det=_legacy_p2_jac_extrema(vertices,normal)
    if name=="minDetJac" || name=="maxDetJac"
        value=name=="minDetJac" ? min_det : max_det
        return _legacy_p2_scaled_quantity(value,scale,length(vertices)-1)
    end
    ((min_det<=0 && max_det>=0) || max_det<0) && return 0.0
    initial=_legacy_p2_icn_domain(vertices,normal)
    domains=[initial]
    _legacy_p2_subdivide_quality!(domains,normal,initial.min_l,initial.max_l,true,false)
    min_b=minimum(domain.min_b for domain in domains)
    min_l=minimum(domain.min_l for domain in domains)
    factor=.5*(min_b+min_l)
    return factor*min_l+(1-factor)*min_b
end
