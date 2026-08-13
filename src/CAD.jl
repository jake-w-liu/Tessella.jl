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
    length(a) >= 3 || throw(ArgumentError("CAD: a point/vector needs three coordinates"))
    q = try
        (Float64(a[1]), Float64(a[2]), Float64(a[3]))
    catch err
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
@inline _norm(a) = sqrt(_dot(a,a))
@inline function _unit(a)
    l=_norm(a)
    (isfinite(l) && l > 0) || throw(ArgumentError("CAD: vector must have finite positive length"))
    (a[1]/l,a[2]/l,a[3]/l)
end
@inline function _radius(r, caller)
    rf = Float64(r)
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
surface_residual(s::PlaneS, x) =
    _result_scalar(_dot(s.n, _sub(_v(x), s.p0)), "surface_residual(::PlaneS)")
function surface_residual(s::CylinderS, x)
    d = _sub(_v(x), s.base); axl = _dot(d, s.axis)
    rad = _sub(d, _scale(s.axis, axl))
    _result_scalar(_norm(rad) - s.r, "surface_residual(::CylinderS)")
end
surface_residual(s::SphereS, x) =
    _result_scalar(_norm(_sub(_v(x), s.center)) - s.r, "surface_residual(::SphereS)")
surface_residual(s::DiskS, x) =
    _result_scalar(_dot(s.n, _sub(_v(x), s.center)), "surface_residual(::DiskS)")

@inline function _tol(tol, caller)
    t = Float64(tol)
    (isfinite(t) && t >= 0) || throw(ArgumentError("$caller: tol must be finite and non-negative (got $tol)"))
    return t
end
on_surface(s, x; tol=1e-9) = abs(surface_residual(s, x)) <= _tol(tol,"on_surface")
function on_surface(s::DiskS, x; tol=1e-9)
    t = _tol(tol,"on_surface(::DiskS)")
    d = _sub(_v(x), s.center); axial = _dot(d,s.n)
    radial = _sub(d,_scale(s.n,axial))
    return abs(axial) <= t && _norm(radial) <= s.r + t
end

"Project `x` exactly onto surface `s` (nearest point on the true surface)."
project_to(s::PlaneS, x) = _result_point(
    _sub(_v(x), _scale(s.n, _dot(s.n, _sub(_v(x), s.p0)))), "project_to(::PlaneS)")
function project_to(s::CylinderS, x)
    d = _sub(_v(x), s.base); axl = _dot(d, s.axis)
    footpt = _add(s.base, _scale(s.axis, axl))       # closest axis point
    rad = _sub(_v(x), footpt); rl = _norm(rad)
    # On (or numerically on) the axis: `rad` is not exactly zero when `axis` has irrational
    # unit components (e.g. (1,1,1)/√3) — the `axl·axis` foot-point rounding leaves a tiny
    # axis-PARALLEL residual, whose magnitude is ~ε·|axl|. A `== 0` test misses it and the
    # normal branch would then normalise that parallel residual and move along the axis
    # (radial distance 0, residual −r). A scale-relative tolerance detects the on-axis case.
    # The scale is tied to actual floating-point roundoff, not a loose decimal
    # tolerance that would misclassify genuine near-axis points far from `base`.
    if rl <= 4eps(Float64) * (_norm(d) + s.r + 1.0)    # on axis: pick ANY radial direction
        # reference the coordinate axis LEAST aligned with s.axis so the cross is never
        # zero (a fixed (1,0,0) fails when the cylinder axis is itself ∥ x, e.g. an
        # x-directed coax bore).
        a1=abs(s.axis[1]); a2=abs(s.axis[2]); a3=abs(s.axis[3])
        ref = a1<=a2 ? (a1<=a3 ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) : (a2<=a3 ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
        return _result_point(_add(footpt, _scale(_unit(_cross(s.axis, ref)), s.r)),
                             "project_to(::CylinderS)")
    end
    _result_point(_add(footpt, _scale(rad, s.r/rl)), "project_to(::CylinderS)")
end
function project_to(s::SphereS, x)
    d = _sub(_v(x), s.center); dl = _norm(d)
    dl == 0 && return _result_point(_add(s.center, (s.r,0.0,0.0)), "project_to(::SphereS)")
    _result_point(_add(s.center, _scale(d, s.r/dl)), "project_to(::SphereS)")
end
function project_to(s::DiskS, x)
    d = _sub(_v(x), s.center); axial = _dot(d,s.n)
    planar = _sub(d,_scale(s.n,axial)); rl = _norm(planar)
    q = rl <= s.r ? _add(s.center,planar) : _add(s.center,_scale(planar,s.r/rl))
    return _result_point(q, "project_to(::DiskS)")
end

function _imprint_count(nseg::Integer, caller::AbstractString)
    nseg >= 3 || throw(ArgumentError("$caller: nseg must be at least 3 (got $nseg)"))
    nseg <= typemax(Int) || throw(ArgumentError("$caller: nseg exceeds the platform Int limit"))
    n = Int(nseg)
    try
        Base.checked_mul(n, sizeof(NTuple{3,Float64}))
    catch
        throw(ArgumentError("$caller: requested point storage overflows the platform Int limit"))
    end
    return n
end

# ── exact intersection curves (boolean imprints) ─────────────────────────────────
"""
    imprint_circle(cyl, plane; nseg=48) -> (center, points)

EXACT intersection of a cylinder with a plane **perpendicular** to the cylinder axis: a
circle of radius `cyl.r` centred where the axis pierces the plane. Every returned point
lies exactly (to round-off) on BOTH the cylinder and the plane — the native boolean imprint
of a coax bore/shield through a flat case wall. Throws if the plane is not ⊥ the axis
(oblique ⇒ use [`imprint_ellipse`](@ref)).
"""
function imprint_circle(cyl::CylinderS, plane::PlaneS; nseg::Integer=48)
    nseg = _imprint_count(nseg,"imprint_circle")
    ax = cyl.axis
    abs(abs(_dot(ax, plane.n)) - 1.0) < 1e-12 ||
        throw(ArgumentError("imprint_circle: plane must be perpendicular to the cylinder axis (use imprint_ellipse)"))
    denom = _dot(plane.n, ax)
    t = -_dot(plane.n, _sub(cyl.base, plane.p0)) / denom
    center = _add(cyl.base, _scale(ax, t))
    a0 = abs(ax[1]) < 0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
    e1 = _unit(_cross(ax, a0)); e2 = _cross(ax, e1)
    pts = [ _add(center, _add(_scale(e1, cyl.r*cos(2π*k/nseg)), _scale(e2, cyl.r*sin(2π*k/nseg)))) for k in 0:nseg-1 ]
    (center, pts)
end

"""
    imprint_ellipse(cyl, plane; nseg=48) -> (center, points)

EXACT intersection of a cylinder with a general (non-parallel) plane: an ellipse. Every
returned point lies exactly on both the cylinder and the plane. (Parallel plane ⇒ no
bounded intersection ⇒ throws.)
"""
function imprint_ellipse(cyl::CylinderS, plane::PlaneS; nseg::Integer=48)
    nseg = _imprint_count(nseg,"imprint_ellipse")
    ax = cyl.axis; denom = _dot(plane.n, ax)
    abs(denom) > 1e-14 || throw(ArgumentError("imprint_ellipse: plane is parallel to the cylinder axis"))
    a0 = abs(ax[1]) < 0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
    e1 = _unit(_cross(ax, a0)); e2 = _cross(ax, e1)
    pts = NTuple{3,Float64}[]
    for k in 0:nseg-1
        θ = 2π*k/nseg
        # a point on the cylinder surface at angle θ, axial coord chosen so it lands on the plane
        rvec = _add(_scale(e1, cyl.r*cos(θ)), _scale(e2, cyl.r*sin(θ)))
        p_ax0 = _add(cyl.base, rvec)                      # on cylinder, at axial 0
        t = -_dot(plane.n, _sub(p_ax0, plane.p0)) / denom # slide along axis onto the plane
        push!(pts, _add(p_ax0, _scale(ax, t)))
    end
    # center = axis∩plane
    tc = -_dot(plane.n, _sub(cyl.base, plane.p0)) / denom
    (_add(cyl.base, _scale(ax, tc)), pts)
end

end # module CAD
