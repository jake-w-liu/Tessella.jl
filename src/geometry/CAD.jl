"""
    CAD

Native analytical geometry ("CAD-lite") — implemented in Julia, no OpenCASCADE. Provides
the exact analytical surfaces the ASCENT primitive set uses (planes, cylinders, spheres,
disks), exact **surface membership** and **projection** (place a mesh node exactly on the
true surface), and exact **intersection curves** (the boolean *imprint* where two surfaces
meet — e.g. a cylinder bore piercing a planar wall is an exact circle, not a facet polygon).

This is the analytical geometry model gmsh gets from OpenCASCADE, built natively: the mesh
nodes lie on the exact surfaces (to round-off), boolean interface curves are exact, and
high-order (P2) elements curve onto the true surface ([`HighOrder.curve_to_cylinder!`](@ref)).
The volume mesh, as always, is a piecewise approximation that converges to the exact model.
"""
module CAD

export PlaneS, CylinderS, SphereS, DiskS
export on_surface, project_to, surface_residual, imprint_circle, imprint_ellipse

@inline function _v(a)
    count = try
        length(a)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("CAD: a point/vector needs three coordinates"))
    end
    count >= 3 || throw(ArgumentError("CAD: a point/vector needs three coordinates"))
    raw = try
        (a[1], a[2], a[3])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "CAD: point/vector coordinates must support integer indexing"))
    end
    all(value -> value isa Real, raw) || throw(ArgumentError(
        "CAD: coordinates must be real"))
    any(value -> value isa Bool, raw) && throw(ArgumentError(
        "CAD: coordinates must not be Bool"))
    q = try
        (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("CAD: coordinates must be Float64-representable: $(sprint(showerror, err))"))
    end
    (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) ||
        throw(ArgumentError("CAD: coordinates must be finite (got $q)"))
    return q
end
@inline _dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _sub(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _add(a,b) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
@inline _scale(a,s) = (a[1]*s, a[2]*s, a[3]*s)
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _norm(a) = hypot(a[1],a[2],a[3])
@inline _maxabs(a) = max(abs(a[1]),abs(a[2]),abs(a[3]))
@inline function _two_sum(a::Float64,b::Float64)
    total=a+b
    right=total-a
    error=(a-(total-right))+(b-right)
    return total,error
end
@inline function _compensated_dot(a,b)
    p1=a[1]*b[1];p2=a[2]*b[2];p3=a[3]*b[3]
    all(isfinite,(p1,p2,p3)) || return p1+p2+p3
    s12,e12=_two_sum(p1,p2)
    total,e3=_two_sum(s12,p3)
    product_error=muladd(a[1],b[1],-p1)+
                  muladd(a[2],b[2],-p2)+muladd(a[3],b[3],-p3)
    return total+(e12+e3+product_error)
end
@inline function _unit(a)
    scale=_maxabs(a)
    (isfinite(scale) && scale > 0) || throw(ArgumentError(
        "CAD: vector must have finite positive length"))
    scaled=(a[1]/scale,a[2]/scale,a[3]/scale)
    magnitude=_norm(scaled)
    (isfinite(magnitude) && magnitude > 0) || throw(ArgumentError(
        "CAD: vector must have finite positive length"))
    return (scaled[1]/magnitude,scaled[2]/magnitude,scaled[3]/magnitude)
end
@inline function _radius(r, caller)
    r isa Real || throw(ArgumentError("$caller: radius must be real"))
    r isa Bool && throw(ArgumentError("$caller: radius must not be Bool"))
    rf = try
        Float64(r)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: radius must be representable as Float64"))
    end
    (isfinite(rf) && rf > 0) || throw(ArgumentError("$caller: radius must be finite and positive (got $r)"))
    return rf
end
@inline function _result_point(p, caller)
    (isfinite(p[1]) && isfinite(p[2]) && isfinite(p[3])) ||
        throw(ArgumentError("$caller: computed point is non-finite"))
    return p
end
@inline function _result_scalar(x, caller)
    isfinite(x) || throw(ArgumentError("$caller: computed value is non-finite"))
    return x
end

@inline _radial_tolerance(radius::Float64) =
    max(128eps(Float64)*radius,nextfloat(0.0))

const _FALLBACK_PRECISION=2304

@inline _exact_float(value::Float64)=Rational{BigInt}(value)

function _exact_parallel(a,b)
    exact_a=ntuple(index->_exact_float(a[index]),3)
    exact_b=ntuple(index->_exact_float(b[index]),3)
    return exact_a[2]*exact_b[3]-exact_a[3]*exact_b[2]==0 &&
           exact_a[3]*exact_b[1]-exact_a[1]*exact_b[3]==0 &&
           exact_a[1]*exact_b[2]-exact_a[2]*exact_b[1]==0
end

function _exact_dot_iszero(a,b)
    return sum(_exact_float(a[index])*_exact_float(b[index])
               for index in 1:3)==0
end

function _exact_cylinder_components(surface,point)
    axis=ntuple(index->_exact_float(surface.axis[index]),3)
    displacement=ntuple(index->_exact_float(point[index])-
                               _exact_float(surface.base[index]),3)
    axis_norm_squared=sum(component^2 for component in axis)
    axial=sum(displacement[index]*axis[index] for index in 1:3)/
          axis_norm_squared
    radial=ntuple(
        index->displacement[index]-axial*axis[index],3)
    return axis,axial,radial
end

function _exact_cylinder_distance(surface,point)
    _,_,radial=_exact_cylinder_components(surface,point)
    squared=sum(component^2 for component in radial)
    return setprecision(BigFloat,_FALLBACK_PRECISION) do
        Float64(sqrt(BigFloat(squared)))
    end
end

function _exact_sphere_distance(surface,point)
    displacement=ntuple(index->_exact_float(point[index])-
                               _exact_float(surface.center[index]),3)
    squared=sum(component^2 for component in displacement)
    return setprecision(BigFloat,_FALLBACK_PRECISION) do
        Float64(sqrt(BigFloat(squared)))
    end
end

function _project_sphere_exact(surface,point,caller)::NTuple{3,Float64}
    direction=ntuple(index->_exact_float(point[index])-
                           _exact_float(surface.center[index]),3)
    if all(iszero,direction)
        direction=(Rational{BigInt}(1),Rational{BigInt}(0),
                   Rational{BigInt}(0))
    end
    squared=sum(component^2 for component in direction)
    projected=setprecision(BigFloat,_FALLBACK_PRECISION) do
        magnitude=sqrt(BigFloat(squared))
        ntuple(3) do index
            value=BigFloat(surface.center[index])+BigFloat(surface.r)*
                  BigFloat(direction[index])/magnitude
            Float64(value)
        end
    end
    result=_result_point(projected,caller)
    return _certify_sphere_point(surface,result,caller)
end

function _exact_point_distance(origin,point)
    displacement=ntuple(index->_exact_float(point[index])-
                               _exact_float(origin[index]),3)
    squared=sum(component^2 for component in displacement)
    return setprecision(BigFloat,_FALLBACK_PRECISION) do
        Float64(sqrt(BigFloat(squared)))
    end
end

function _project_disk_radial_exact(surface,planar_point,caller)::NTuple{3,Float64}
    direction=ntuple(index->_exact_float(planar_point[index])-
                           _exact_float(surface.center[index]),3)
    squared=sum(component^2 for component in direction)
    squared<=_exact_float(surface.r)^2 && return planar_point
    projected=setprecision(BigFloat,_FALLBACK_PRECISION) do
        magnitude=sqrt(BigFloat(squared))
        ntuple(3) do index
            value=BigFloat(surface.center[index])+BigFloat(surface.r)*
                  BigFloat(direction[index])/magnitude
            Float64(value)
        end
    end
    return _result_point(projected,caller)
end

function _project_cylinder_exact(surface,point,caller)::NTuple{3,Float64}
    axis,axial,radial=_exact_cylinder_components(surface,point)
    direction=radial
    if all(iszero,direction)
        a1=abs(surface.axis[1]);a2=abs(surface.axis[2]);a3=abs(surface.axis[3])
        reference=a1<=a2 ?
            (a1<=a3 ? (1,0,0) : (0,0,1)) :
            (a2<=a3 ? (0,1,0) : (0,0,1))
        reference_exact=ntuple(index->Rational{BigInt}(reference[index]),3)
        direction=(axis[2]*reference_exact[3]-axis[3]*reference_exact[2],
                   axis[3]*reference_exact[1]-axis[1]*reference_exact[3],
                   axis[1]*reference_exact[2]-axis[2]*reference_exact[1])
    end
    squared=sum(component^2 for component in direction)
    squared>0 || throw(ErrorException(
        "$caller: internal cylinder projection direction invariant failed"))
    projected=setprecision(BigFloat,_FALLBACK_PRECISION) do
        magnitude=sqrt(BigFloat(squared))
        ntuple(3) do index
            foot=_exact_float(surface.base[index])+axis[index]*axial
            value=BigFloat(foot)+BigFloat(surface.r)*
                  BigFloat(direction[index])/magnitude
            Float64(value)
        end
    end
    result=_result_point(projected,caller)
    return _certify_cylinder_point(surface,result,caller)
end

function _exact_plane_residual(origin,normal,point)
    return sum((_exact_float(point[index])-_exact_float(origin[index]))*
               _exact_float(normal[index]) for index in 1:3)
end

function _plane_residual(origin,normal,point,caller)
    displacement=_sub(point,origin)
    if !all(isfinite,displacement)
        return _result_scalar(
            Float64(_exact_plane_residual(origin,normal,point)),caller)
    end
    products=ntuple(index->displacement[index]*normal[index],3)
    residual=_compensated_dot(displacement,normal)
    magnitude_sum=sum(abs,products)
    if !isfinite(residual) || !isfinite(magnitude_sum)
        residual=Float64(_exact_plane_residual(origin,normal,point))
    end
    return _result_scalar(residual,caller)
end

function _project_to_plane_exact(origin,normal,point,caller)::NTuple{3,Float64}
    normal_exact=(_exact_float(normal[1]),_exact_float(normal[2]),
                  _exact_float(normal[3]))
    displacement_exact=(_exact_float(point[1])-_exact_float(origin[1]),
                        _exact_float(point[2])-_exact_float(origin[2]),
                        _exact_float(point[3])-_exact_float(origin[3]))
    denominator_exact=normal_exact[1]^2+normal_exact[2]^2+
                      normal_exact[3]^2
    numerator_exact=displacement_exact[1]*normal_exact[1]+
                    displacement_exact[2]*normal_exact[2]+
                    displacement_exact[3]*normal_exact[3]
    factor_exact=numerator_exact/denominator_exact
    projected=(Float64(_exact_float(point[1])-
                       normal_exact[1]*factor_exact),
               Float64(_exact_float(point[2])-
                       normal_exact[2]*factor_exact),
               Float64(_exact_float(point[3])-
                       normal_exact[3]*factor_exact))
    return _result_point(projected,caller)
end

function _project_to_plane(origin,normal,point,caller)::NTuple{3,Float64}
    displacement=_sub(point,origin)
    all(isfinite,displacement) ||
        return _project_to_plane_exact(origin,normal,point,caller)
    products=ntuple(index->displacement[index]*normal[index],3)
    numerator=_compensated_dot(displacement,normal)
    magnitude_sum=sum(abs,products)
    if !isfinite(numerator) || !isfinite(magnitude_sum)
        return _project_to_plane_exact(origin,normal,point,caller)
    end
    factor=numerator/_dot(normal,normal)
    candidate=_sub(point,_scale(normal,factor))
    all(isfinite,candidate) ||
        return _project_to_plane_exact(origin,normal,point,caller)
    return candidate
end

function _line_plane_point_exact(start,direction,origin,normal,
                                 caller)::NTuple{3,Float64}
    start_exact=ntuple(index->_exact_float(start[index]),3)
    direction_exact=ntuple(index->_exact_float(direction[index]),3)
    origin_exact=ntuple(index->_exact_float(origin[index]),3)
    normal_exact=ntuple(index->_exact_float(normal[index]),3)
    numerator=sum(normal_exact[index]*
                  (start_exact[index]-origin_exact[index]) for index in 1:3)
    denominator=sum(normal_exact[index]*direction_exact[index]
                    for index in 1:3)
    denominator!=0 || throw(ArgumentError(
        "$caller: line is parallel to the plane"))
    parameter=-numerator/denominator
    point=ntuple(index->Float64(
        start_exact[index]+parameter*direction_exact[index]),3)
    return _result_point(point,caller)
end

function _line_plane_point(start,direction,origin,normal,
                           caller)::NTuple{3,Float64}
    displacement=_sub(start,origin)
    if all(isfinite,displacement)
        numerator=_compensated_dot(normal,displacement)
        denominator=_compensated_dot(normal,direction)
        parameter=-numerator/denominator
        if isfinite(parameter)
            candidate=_add(start,_scale(direction,parameter))
            all(isfinite,candidate) && return candidate
        end
    end
    return _line_plane_point_exact(
        start,direction,origin,normal,caller)
end

function _certify_sphere_point(surface,point,caller)::NTuple{3,Float64}
    radial=_sub(point,surface.center)
    distance=all(isfinite,radial) ? _norm(radial) :
             _exact_sphere_distance(surface,point)
    tolerance=_radial_tolerance(surface.r)
    (isfinite(distance) && distance>0 &&
     abs(distance-surface.r)<=tolerance) || throw(ArgumentError(
        "$caller: projected sphere point is not representable within " *
        "the radius tolerance"))
    return point
end

function _certify_cylinder_point(surface,point,caller)::NTuple{3,Float64}
    displacement=_sub(point,surface.base)
    axis_norm_squared=_dot(surface.axis,surface.axis)
    tolerance=_radial_tolerance(surface.r)
    axial=_compensated_dot(displacement,surface.axis)/axis_norm_squared
    certified=false
    if all(isfinite,displacement) && isfinite(axial)
        radial=_sub(displacement,_scale(surface.axis,axial))
        distance=_norm(radial)
        certified=isfinite(distance) && distance>0 &&
                  abs(distance-surface.r)<=tolerance
    end
    if !certified
        exact_distance=_exact_cylinder_distance(surface,point)
        certified=isfinite(exact_distance) && exact_distance>0 &&
                  abs(exact_distance-surface.r)<=tolerance
    end
    certified || throw(ArgumentError(
        "$caller: projected cylinder point is not representable within " *
        "the radius tolerance"))
    return point
end

function _certify_plane_components(origin,normal,point,caller)::NTuple{3,Float64}
    displacement=_sub(point,origin)
    finite_displacement=all(isfinite,displacement)
    residual=finite_displacement ? _compensated_dot(displacement,normal) : Inf
    scale=finite_displacement ? _maxabs(displacement) :
        max(_maxabs(point),_maxabs(origin))
    scale=max(scale,nextfloat(0.0))
    tolerance=256eps(Float64)*scale
    certified=isfinite(residual) && abs(residual)<=tolerance
    if !certified
        certified=abs(_exact_plane_residual(origin,normal,point))<=
                  _exact_float(tolerance)
    end
    certified ||
        throw(ArgumentError(
            "$caller: projected point is not representable on the plane"))
    return point
end


function _certify_plane_point(surface,point,caller)::NTuple{3,Float64}
    return _certify_plane_components(surface.p0,surface.n,point,caller)
end

function _certify_disk_point(surface,point,caller)::NTuple{3,Float64}
    _certify_plane_components(
        surface.center,surface.n,point,caller)
    radial=_sub(point,surface.center)
    distance=all(isfinite,radial) ? _norm(radial) :
             _exact_point_distance(surface.center,point)
    tolerance=_radial_tolerance(surface.r)
    (isfinite(distance) && distance<=surface.r+tolerance) ||
        throw(ArgumentError(
            "$caller: projected point is not representable within the disk"))
    return point
end

# ── analytical surfaces ──────────────────────────────────────────────────────────
# Inner constructors normalize/convert (an OUTER constructor of the same field types
# would be shadowed by the auto-generated one, silently storing a non-unit normal).
"Plane through `p0` with unit normal `n` (n·(x−p0) = 0)."
struct PlaneS
    p0::NTuple{3,Float64}; n::NTuple{3,Float64}
    PlaneS(p0, n) = new(_v(p0), _unit(_v(n)))
end

"Infinite cylinder: `base` on the axis, unit `axis`, radius `r` (points with radial distance == r)."
struct CylinderS
    base::NTuple{3,Float64}; axis::NTuple{3,Float64}; r::Float64
    CylinderS(base, axis, r) = new(_v(base), _unit(_v(axis)), _radius(r,"CylinderS"))
end

"Sphere of `center`, radius `r`."
struct SphereS
    center::NTuple{3,Float64}; r::Float64
    SphereS(center, r) = new(_v(center), _radius(r,"SphereS"))
end

"Disk: a bounded planar circle — center, unit normal, radius (the plane region ≤ r)."
struct DiskS
    center::NTuple{3,Float64}; n::NTuple{3,Float64}; r::Float64
    DiskS(center, n, r) = new(_v(center), _unit(_v(n)), _radius(r,"DiskS"))
end

# ── exact membership + projection ────────────────────────────────────────────────
"Signed/absolute residual of `x` from surface `s` (0 exactly on the surface)."
surface_residual(s::PlaneS,x)=_plane_residual(
    s.p0,s.n,_v(x),"surface_residual(::PlaneS)")
function surface_residual(s::CylinderS, x)
    point=_v(x)
    d = _sub(point,s.base)
    all(isfinite,d) || return _result_scalar(
        _exact_cylinder_distance(s,point)-s.r,
        "surface_residual(::CylinderS)")
    products=ntuple(index->d[index]*s.axis[index],3)
    numerator=_compensated_dot(d,s.axis)
    magnitude_sum=sum(abs,products)
    distance=if !isfinite(numerator) || !isfinite(magnitude_sum)
        _exact_cylinder_distance(s,point)
    else
        axl=numerator/_dot(s.axis,s.axis)
        radial=_sub(d,_scale(s.axis,axl))
        if _maxabs(radial)<=64eps(Float64)*
                            max(_maxabs(d),nextfloat(0.0))
            _exact_cylinder_distance(s,point)
        else
            _norm(radial)
        end
    end
    _result_scalar(distance-s.r,"surface_residual(::CylinderS)")
end
function surface_residual(s::SphereS,x)
    point=_v(x)
    displacement=_sub(point,s.center)
    distance=all(isfinite,displacement) ? _norm(displacement) :
             _exact_sphere_distance(s,point)
    return _result_scalar(distance-s.r,"surface_residual(::SphereS)")
end
surface_residual(s::DiskS,x)=_plane_residual(
    s.center,s.n,_v(x),"surface_residual(::DiskS)")

@inline function _tol(tol, caller)
    tol isa Real || throw(ArgumentError("$caller: tol must be real"))
    tol isa Bool && throw(ArgumentError("$caller: tol must not be Bool"))
    t = try
        Float64(tol)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: tol must be representable as Float64"))
    end
    (isfinite(t) && t >= 0) || throw(ArgumentError("$caller: tol must be finite and non-negative (got $tol)"))
    return t
end
on_surface(s, x; tol=1e-9) = abs(surface_residual(s, x)) <= _tol(tol,"on_surface")
function on_surface(s::DiskS, x; tol=1e-9)
    t = _tol(tol,"on_surface(::DiskS)")
    point=_v(x)
    axial=_plane_residual(s.center,s.n,point,"on_surface(::DiskS)")
    planar=_project_to_plane(s.center,s.n,point,"on_surface(::DiskS)")
    radial=_sub(planar,s.center)
    rl=all(isfinite,radial) ? _norm(radial) :
       _exact_point_distance(s.center,planar)
    radial_ok=rl<=s.r || rl-s.r<=t
    return abs(axial) <= t && radial_ok
end

"Project `x` exactly onto surface `s` (nearest point on the true surface)."
function project_to(s::PlaneS, x)
    point=_v(x)
    projected=_project_to_plane(
        s.p0,s.n,point,"project_to(::PlaneS)")
    return _certify_plane_point(s,projected,"project_to(::PlaneS)")
end
function project_to(s::CylinderS, x)
    point=_v(x)
    d = _sub(point,s.base)
    all(isfinite,d) ||
        return _project_cylinder_exact(s,point,"project_to(::CylinderS)")
    products=ntuple(index->d[index]*s.axis[index],3)
    numerator=_compensated_dot(d,s.axis)
    product_sum=sum(abs,products)
    if !isfinite(numerator) || !isfinite(product_sum)
        return _project_cylinder_exact(s,point,"project_to(::CylinderS)")
    end
    axl = numerator/_dot(s.axis,s.axis)
    footpt = _add(s.base, _scale(s.axis, axl))       # closest axis point
    rad = _sub(d,_scale(s.axis,axl))
    all(isfinite,rad) ||
        return _project_cylinder_exact(s,point,"project_to(::CylinderS)")
    # A radial component inside this round-off envelope is ambiguous in ordinary
    # Float64 arithmetic. Exact dyadic decomposition distinguishes a true on-axis
    # query from a genuine small radial displacement, even at a huge axial offset.
    if _maxabs(rad)<=64eps(Float64)*
                     max(_maxabs(d),s.r,nextfloat(0.0))
        return _project_cylinder_exact(s,point,"project_to(::CylinderS)")
    end
    candidate=_add(footpt,_scale(_unit(rad),s.r))
    all(isfinite,candidate) ||
        return _project_cylinder_exact(s,point,"project_to(::CylinderS)")
    projected=_result_point(candidate,"project_to(::CylinderS)")
    return _certify_cylinder_point(s,projected,"project_to(::CylinderS)")
end
function project_to(s::SphereS, x)
    point=_v(x)
    d = _sub(point,s.center)
    all(isfinite,d) ||
        return _project_sphere_exact(s,point,"project_to(::SphereS)")
    direction=all(iszero,d) ? (1.0,0.0,0.0) : _unit(d)
    projected=_result_point(
        _add(s.center,_scale(direction,s.r)),"project_to(::SphereS)")
    return _certify_sphere_point(s,projected,"project_to(::SphereS)")
end
function project_to(s::DiskS, x)
    point=_v(x)
    planar_point=_project_to_plane(
        s.center,s.n,point,"project_to(::DiskS)")
    planar=_sub(planar_point,s.center)
    all(isfinite,planar) || begin
        projected=_project_disk_radial_exact(
            s,planar_point,"project_to(::DiskS)")
        return _certify_disk_point(s,projected,"project_to(::DiskS)")
    end
    radial_distance = _norm(planar)
    q = radial_distance <= s.r ? planar_point :
        _add(s.center,_scale(_unit(planar),s.r))
    projected=_result_point(q,"project_to(::DiskS)")
    return _certify_disk_point(s,projected,"project_to(::DiskS)")
end

function _imprint_count(nseg,max_points,caller::AbstractString)
    nseg isa Integer || throw(ArgumentError("$caller: nseg must be an integer"))
    nseg isa Bool && throw(ArgumentError("$caller: nseg must not be Bool"))
    max_points isa Integer || throw(ArgumentError(
        "$caller: max_points must be an integer"))
    max_points isa Bool && throw(ArgumentError(
        "$caller: max_points must not be Bool"))
    nseg >= 3 || throw(ArgumentError("$caller: nseg must be at least 3 (got $nseg)"))
    point_limit=try
        Int(max_points)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: max_points is outside the platform Int range"))
    end
    point_limit>=1 || throw(ArgumentError(
        "$caller: max_points must be positive (got $max_points)"))
    n = try
        Int(nseg)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: nseg exceeds the platform Int limit"))
    end
    n<=point_limit || throw(ArgumentError(
        "$caller: nseg=$n exceeds max_points=$point_limit"))
    try
        Base.checked_mul(n, sizeof(NTuple{3,Float64}))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: requested point storage overflows the platform Int limit"))
    end
    return n
end

# ── exact intersection curves (boolean imprints) ─────────────────────────────────
"""
    imprint_circle(cyl, plane; nseg=48, max_points=10_000_000) -> (center, points)

EXACT intersection of a cylinder with a plane **perpendicular** to the cylinder axis: a
circle of radius `cyl.r` centred where the axis pierces the plane. Every returned point
lies exactly (to round-off) on BOTH the cylinder and the plane — the native boolean imprint
of a coax bore/shield through a flat case wall. Throws if the plane is not ⊥ the axis
(oblique ⇒ use [`imprint_ellipse`](@ref)).
"""
function imprint_circle(cyl::CylinderS,plane::PlaneS;
                        nseg=48,max_points=10_000_000)
    nseg = _imprint_count(nseg,max_points,"imprint_circle")
    ax = cyl.axis
    _exact_parallel(ax,plane.n) ||
        throw(ArgumentError("imprint_circle: plane must be perpendicular to the cylinder axis (use imprint_ellipse)"))
    center=_line_plane_point(
        cyl.base,ax,plane.p0,plane.n,"imprint_circle")
    a0 = abs(ax[1]) < 0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
    e1 = _unit(_cross(ax, a0)); e2 = _cross(ax, e1)
    _certify_plane_point(plane,center,"imprint_circle")
    pts = [_certify_plane_point(
        plane,
        _certify_cylinder_point(
            cyl,
            _result_point(
                _add(center,_add(_scale(e1,cyl.r*cos(2π*k/nseg)),
                                 _scale(e2,cyl.r*sin(2π*k/nseg)))),
                "imprint_circle"),
            "imprint_circle"),
        "imprint_circle") for k in 0:nseg-1]
    (center, pts)
end

"""
    imprint_ellipse(cyl, plane; nseg=48, max_points=10_000_000) -> (center, points)

EXACT intersection of a cylinder with a general (non-parallel) plane: an ellipse. Every
returned point lies exactly on both the cylinder and the plane. (Parallel plane ⇒ no
bounded intersection ⇒ throws.)
"""
function imprint_ellipse(cyl::CylinderS,plane::PlaneS;
                         nseg=48,max_points=10_000_000)
    nseg = _imprint_count(nseg,max_points,"imprint_ellipse")
    ax = cyl.axis
    _exact_dot_iszero(plane.n,ax) && throw(ArgumentError(
        "imprint_ellipse: plane is parallel to the cylinder axis"))
    a0 = abs(ax[1]) < 0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
    e1 = _unit(_cross(ax, a0)); e2 = _cross(ax, e1)
    pts = NTuple{3,Float64}[]
    sizehint!(pts,nseg)
    for k in 0:nseg-1
        θ = 2π*k/nseg
        # a point on the cylinder surface at angle θ, axial coord chosen so it lands on the plane
        rvec=_result_point(
            _add(_scale(e1,cyl.r*cos(θ)),_scale(e2,cyl.r*sin(θ))),
            "imprint_ellipse")
        p_ax0=_result_point(_add(cyl.base,rvec),"imprint_ellipse")
        point=_line_plane_point(
            p_ax0,ax,plane.p0,plane.n,"imprint_ellipse")
        _certify_cylinder_point(cyl,point,"imprint_ellipse")
        _certify_plane_point(plane,point,"imprint_ellipse")
        push!(pts,point)
    end
    # center = axis∩plane
    center=_line_plane_point(
        cyl.base,ax,plane.p0,plane.n,"imprint_ellipse")
    _certify_plane_point(plane,center,"imprint_ellipse")
    return center,pts
end

end # module CAD
