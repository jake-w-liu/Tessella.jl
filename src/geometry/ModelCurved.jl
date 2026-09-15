# ── Curved entity kernel ─────────────────────────────────────────────────────
#
# Gmsh built-in-kernel curved geometry. Arcs (`:circle`, `:ellipse`) carry an
# ordered control-point list in `curve_control_points` — `[start, center, end]`
# or `[start, center, major, end]` — plus the `Plane{..}` hint normal in
# `curve_geometry`, mirroring Gmsh's `Curve` records (`MSH_SEGM_CIRC`/`ELLI`).
# Surfaces carry a type tag in `surface_types` — `:plane` (default), `:ruled`
# (4-border filling), or `:tric` (3-border filling) — matching Gmsh's
# `MSH_SURF_PLAN`/`REGL`/`TRIC` records, with auxiliary parameters such as the
# `In Sphere` center point in `surface_geometry`.
#
# Derived arc parameters are recomputed on demand exactly like Gmsh's
# `EndCurve`, so point moves, transforms, and coherence passes need no cache
# invalidation.

_curve_type(m::GeoModel, tag::Int) = get(m.curve_types, tag, :line)
_surface_type(m::GeoModel, tag::Int) = get(m.surface_types, tag, :plane)

const _CURVE_TYPE_NAMES = Dict{Symbol,String}(
    :line=>"Line", :circle=>"Circle", :ellipse=>"Ellipse",
    :degenerate=>"Unknown")
const _SURFACE_TYPE_NAMES = Dict{Symbol,String}(
    :plane=>"Plane", :ruled=>"Surface", :tric=>"Surface",
    :cylinder=>"Cylinder", :sphere=>"Sphere", :cone=>"Cone",
    :torus=>"Torus", :disk=>"Disk")

_curve_type_name(t::Symbol) = get(_CURVE_TYPE_NAMES, t, "Unknown")
_surface_type_name(t::Symbol) = get(_SURFACE_TYPE_NAMES, t, "Unknown")

# Gmsh evaluates and meshes only straight lines through these paths; curved
# kinds fail explicitly rather than degrading to the chord.
function _model_require_line_curve(m::GeoModel, curve::Int,
                                   caller::AbstractString,
                                   what::AbstractString="this operation")
    kind=_curve_type(m, curve)
    kind==:line && return nothing
    throw(ArgumentError(
        "$caller: $what requires a straight Line; Curve[$curve] is a " *
        "$(_curve_type_name(kind))"))
end

function _model_require_plane_surface(m::GeoModel, surface::Int,
                                      caller::AbstractString,
                                      what::AbstractString="this operation")
    kind=_surface_type(m, surface)
    kind==:plane && return nothing
    throw(ArgumentError(
        "$caller: $what requires a Plane surface; Surface[$surface] is a " *
        "$(_surface_type_name(kind))"))
end

# Every boundary curve tag of a surface, unsigned.
function _model_surface_curves(m::GeoModel, surface::Int)
    tags=Int[]
    for loop in m.surfaces[surface]
        for signed in m.loops[loop]
            push!(tags, abs(signed))
        end
    end
    return tags
end

# ── Arc geometry (Gmsh `EndCurve`) ───────────────────────────────────────────

# Gmsh's `angle_02pi`: wrap into [0, 2π]; the comparison is strict so 2π is
# preserved rather than folded to 0.
function _angle_02pi(a::Float64)
    twopi=2π
    while a>twopi || a<0.0
        a = a>0.0 ? a-twopi : a+twopi
    end
    return a
end

# Gmsh's `myasin`: clamp outside [-1,1] to ±π/2.
_myasin(a::Float64) = a<=-1.0 ? -π/2 : a>=1.0 ? π/2 : asin(a)

# Gmsh's `myatan2`: (0,0) maps to +0.0 rather than atan2's signed zero.
_myatan2(y::Float64,x::Float64) = (y==0.0 && x==0.0) ? 0.0 : atan(y,x)

# `norme`: normalize by the `norm3` magnitude (sqrt of the summed squares,
# not hypot); a zero vector stays zero.
@inline function _arc_norme(v::NTuple{3,Float64})
    n=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])
    n==0.0 && return v
    inv=1.0/n
    return (v[1]*inv, v[2]*inv, v[3]*inv)
end

@inline _arc_cross(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _arc_dot(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _arc_sub(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    (a[1]-b[1], a[2]-b[2], a[3]-b[3])

# Derived arc parameters, ported from Gmsh's `EndCurve`. `cps` are the ordered
# control-point coordinates — [start, center, end] or
# [start, center, major, end] — `stored_n` is the `Plane{..}` hint normal used
# only when the start/end directions are too degenerate to define the arc
# plane. Returns (center, u, v, n, t1, t2, f1, f2, incl): the orthonormal
# frame (u = center→start, v = n×u, n = arc-plane normal), the angle interval
# [t1, t2] swept counterclockwise in that frame, the semi-axis lengths f1/f2,
# and the in-plane major-axis inclination `incl` (0 for circles).
function _arc_geometry(cps::Vector{NTuple{3,Float64}}, kind::Symbol,
                       stored_n::NTuple{3,Float64}, tag::Int,
                       caller::AbstractString)
    length(cps) in (3,4) || throw(ErrorException(
        "$caller: Curve[$tag] has $(length(cps)) control points; " *
        "an arc needs 3 or 4"))
    has_major=length(cps)==4
    P0,C=cps[1],cps[2]
    P2=has_major ? cps[4] : cps[3]
    raw0=_arc_sub(P0,C); raw2=_arc_sub(P2,C)
    u=_arc_norme(raw0); e=_arc_norme(raw2)
    n=_arc_cross(u,e)
    if sqrt(n[1]*n[1]+n[2]*n[2]+n[3]*n[3])<1e-15
        n=_arc_norme(stored_n)
    else
        n=_arc_norme(n)
        # Gmsh rechecks the normalized vector before falling back to the
        # stored `Plane{..}` normal.
        if abs(n[1])<1e-5 && abs(n[2])<1e-5 && abs(n[3])<1e-5
            n=_arc_norme(stored_n)
        end
    end
    v=_arc_norme(_arc_cross(n,u))
    x0=_arc_dot(raw0,u); y0=_arc_dot(raw0,v)
    x2=_arc_dot(raw2,u); y2=_arc_dot(raw2,v)
    R=sqrt(x0*x0+y0*y0); R2=sqrt(x2*x2+y2*y2)
    (R==0.0 || R2==0.0) && throw(ArgumentError(
        "$caller: zero radius in circle or ellipse with tag $tag"))
    if !has_major && abs((R-R2)/(R+R2))>0.1
        throw(ArgumentError(
            "$caller: control points of circle with tag $tag are not " *
            "cocircular (R1=$R, R2=$R2)"))
    end
    local A1, A3, A4, f1, f2
    if has_major
        raw3=_arc_sub(cps[3],C)
        x3=_arc_dot(raw3,u); y3=_arc_dot(raw3,v)
        A4=_angle_02pi(_myatan2(y3,x3))
        s,c=sin(A4),cos(A4)
        x1=x0*c+y0*s; y1=-x0*s+y0*c
        xe=x2*c+y2*s; ye=-x2*s+y2*c
        # sys2x2 [x1² y1²; xe² ye²]·sol = (1,1) → sol = (1/f1², 1/f2²),
        # evaluated in Gmsh's exact association: pre-rounded squared
        # products and a multiply-by-reciprocal solve.
        det=(x1*x1)*(ye*ye) - (xe*xe)*(y1*y1)
        matnorm=(x1*x1)^2+(ye*ye)^2+(y1*y1)^2+(xe*xe)^2
        s0=s1=0.0
        if !(matnorm==0.0 || abs(det)/matnorm<1e-16)
            ud=1.0/det
            s0=(1.0*(ye*ye)-(y1*y1)*1.0)*ud
            s1=((x1*x1)*1.0-(xe*xe)*1.0)*ud
        end
        (s0<=0.0 || s1<=0.0) && throw(ArgumentError(
            "$caller: ellipse with tag $tag is wrong"))
        f1=sqrt(1.0/s0); f2=sqrt(1.0/s1)
        A1= x1<0.0 ? -_myasin(y1/f2)+A4+π : _myasin(y1/f2)+A4
        A3= xe<0.0 ? -_myasin(ye/f2)+A4+π : _myasin(ye/f2)+A4
    else
        A1=_myatan2(y0,x0); A3=_myatan2(y2,x2); A4=0.0; f1=f2=R
    end
    A1=_angle_02pi(A1); A3=_angle_02pi(A3)
    A1>=A3 && (A3+=2π)
    span=A3-A1
    span>1.01π && throw(ArgumentError(
        "$caller: circle or ellipse arc $tag greater than Pi (angle=$span)"))
    return (center=C, u=u, v=v, n=n, t1=A1, t2=A3, f1=f1, f2=f2, incl=A4)
end

function _arc_geometry(m::GeoModel, tag::Int, caller::AbstractString)
    cps=get(m.curve_control_points, tag, nothing)
    cps===nothing && throw(ErrorException(
        "$caller: Curve[$tag] has no arc control points; rebuild the model"))
    kind=_curve_type(m, tag)
    kind in (:circle, :ellipse) || throw(ArgumentError(
        "$caller: Curve[$tag] is a $(_curve_type_name(kind)), not an arc"))
    coords=NTuple{3,Float64}[]
    for p in cps
        haskey(m.points, p) || throw(ArgumentError(
            "$caller: Curve[$tag] references unknown Point[$p]"))
        push!(coords, m.points[p])
    end
    stored=get(m.curve_geometry, tag, nothing)
    n=stored===nothing ? (0.0,0.0,1.0) : stored.n
    return _arc_geometry(coords, kind, n, tag, caller)
end

# Gmsh's `InterpolateCurve` arc branch: θ sweeps from t1 down to t2 at u=1,
# evaluated in the local frame as an inclined ellipse, then mapped back to
# world coordinates through invmat (whose columns are u, v, n).
function _arc_point(g, u::Float64)
    θ=g.t1-(g.t1-g.t2)*u-g.incl
    lx=g.f1*cos(θ)*cos(g.incl)-g.f2*sin(θ)*sin(g.incl)
    ly=g.f1*cos(θ)*sin(g.incl)+g.f2*sin(θ)*cos(g.incl)
    C,U,V=g.center,g.u,g.v
    # Association order follows `Projette` + the center add so values agree
    # bit-for-bit with `InterpolateCurve`.
    return ((lx*U[1]+ly*V[1])+C[1],
            (lx*U[2]+ly*V[2])+C[2],
            (lx*U[3]+ly*V[3])+C[3])
end

# Analytic first derivative dP/du. Gmsh's getDerivative uses a finite
# difference; the analytic form agrees within that finite-difference error.
function _arc_first_derivative(g, u::Float64)
    θp=g.t1-(g.t1-g.t2)*u
    s,c=sin(θp),cos(θp)
    si,ci=sin(g.incl),cos(g.incl)
    dlx=-g.f1*s*ci-g.f2*c*si
    dly=-g.f1*s*si+g.f2*c*ci
    k=g.t2-g.t1
    U,V=g.u,g.v
    return (k*(dlx*U[1]+dly*V[1]),
            k*(dlx*U[2]+dly*V[2]),
            k*(dlx*U[3]+dly*V[3]))
end

# d²P/du² = -(t2-t1)²·(P(u) - center) — the local-frame position vector.
function _arc_second_derivative(g, u::Float64)
    θp=g.t1-(g.t1-g.t2)*u
    s,c=sin(θp),cos(θp)
    si,ci=sin(g.incl),cos(g.incl)
    lx=g.f1*c*ci-g.f2*s*si
    ly=g.f1*c*si+g.f2*s*ci
    k2=-(g.t2-g.t1)^2
    U,V=g.u,g.v
    return (k2*(lx*U[1]+ly*V[1]),
            k2*(lx*U[2]+ly*V[2]),
            k2*(lx*U[3]+ly*V[3]))
end

# Gmsh's `GEdge::curvature`: |d1 × d2| / |d1|³ (unsigned).
function _arc_curvature(g, u::Float64)
    d1=_arc_first_derivative(g,u)
    d2=_arc_second_derivative(g,u)
    cr=_arc_cross(d1,d2)
    n1=hypot(d1[1],d1[2],d1[3])
    return hypot(cr[1],cr[2],cr[3])/n1^3
end

# Exact arc bounding box: endpoints plus each coordinate-axis extremum, where
# d/dθ'(f1 u_j cosθ' + f2 v_j sinθ') = 0 → tanθ' = f2 v_j / (f1 u_j).
function _arc_bounding_box(g)
    p0=_arc_point(g,0.0); p1=_arc_point(g,1.0)
    lo=[min(p0[i],p1[i]) for i in 1:3]
    hi=[max(p0[i],p1[i]) for i in 1:3]
    for axis in 1:3
        a1=g.f1*g.u[axis]; a2=g.f2*g.v[axis]
        (a1==0.0 && a2==0.0) && continue
        base=atan(a2,a1)
        for k in -2:4
            θp=base+k*π
            (g.t2-1e-13<=θp<=g.t1+1e-13) || continue
            u=(g.t1-θp)/(g.t1-g.t2)
            (0.0-1e-12<=u<=1.0+1e-12) || continue
            p=_arc_point(g,clamp(u,0.0,1.0))
            for i in 1:3
                lo[i]=min(lo[i],p[i]); hi[i]=max(hi[i],p[i])
            end
        end
    end
    return (lo[1],lo[2],lo[3],hi[1],hi[2],hi[3])
end

# ── Arc construction ─────────────────────────────────────────────────────────

# Effective stored normal for a new arc: Gmsh keeps the default (0,0,1) unless
# `Plane{..}` supplies a nonzero vector.
function _arc_stored_normal(plane_normal, caller::AbstractString)
    plane_normal===nothing && return (0.0,0.0,1.0)
    n=_finite_vector3(plane_normal,caller,"Plane normal")
    return any(!iszero,n) ? n : (0.0,0.0,1.0)
end

"""
    add_circle_arc!(model, start, center, stop; tag=0, plane_normal=nothing) -> tag

Add a circular arc from existing Point `start` to existing Point `stop` around
the Point `center`, mirroring Gmsh's `Circle(tag) = {start, center, stop}`.
The arc spans at most π and travels counterclockwise in the plane defined by
the three points; `plane_normal` supplies the plane normal Gmsh uses when the
three points are collinear.
"""
function add_circle_arc!(m::GeoModel, start, center, stop; tag::Integer=0,
                         plane_normal=nothing)
    caller="add_circle_arc!"
    p1=_tag(start,caller,1); pc=_tag(center,caller,1); p2=_tag(stop,caller,1)
    for p in (p1,pc,p2)
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: unknown Point[$p]"))
    end
    n=_arc_stored_normal(plane_normal,caller)
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller)
    (haskey(m.curves,t) || haskey(m.discrete,(1,t))) && throw(ArgumentError(
        "$caller: Curve[$t] already exists"))
    # Gmsh validates inside `EndCurve` at construction; Tessella rejects the
    # whole statement instead of storing a broken arc.
    _arc_geometry([m.points[p1],m.points[pc],m.points[p2]],
                  :circle,n,t,caller)
    m.curves[t]=(p1,p2)
    m.curve_control_points[t]=[p1,pc,p2]
    m.curve_types[t]=:circle
    m.curve_geometry[t]=(n=n,)
    return t
end

"""
    add_ellipse_arc!(model, start, center, major, stop; tag=0,
                     plane_normal=nothing) -> tag

Add an elliptical arc mirroring Gmsh's `Ellipse(tag) = {start, center, major,
stop}`: `major` is a Point on the major axis, and the arc travels
counterclockwise from `start` to `stop` in the defined plane.
"""
function add_ellipse_arc!(m::GeoModel, start, center, major, stop;
                        tag::Integer=0, plane_normal=nothing)
    caller="add_ellipse_arc!"
    p1=_tag(start,caller,1); pc=_tag(center,caller,1)
    pm=_tag(major,caller,1); p2=_tag(stop,caller,1)
    for p in (p1,pc,pm,p2)
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: unknown Point[$p]"))
    end
    n=_arc_stored_normal(plane_normal,caller)
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller)
    (haskey(m.curves,t) || haskey(m.discrete,(1,t))) && throw(ArgumentError(
        "$caller: Curve[$t] already exists"))
    _arc_geometry([m.points[p1],m.points[pc],m.points[pm],m.points[p2]],
                  :ellipse,n,t,caller)
    m.curves[t]=(p1,p2)
    m.curve_control_points[t]=[p1,pc,pm,p2]
    m.curve_types[t]=:ellipse
    m.curve_geometry[t]=(n=n,)
    return t
end

"""
    add_ruled_surface!(model, loops; tag=0, sphere_center=nothing) -> tag

Add a Gmsh surface-filling entity (`Surface`/`Ruled Surface` in `.geo`): the
first boundary loop's curve count selects the kind — four curves give a ruled
`:ruled` patch and three give a `:tric` patch, matching Gmsh's
`MSH_SURF_REGL`/`MSH_SURF_TRIC`. `sphere_center` records the `In Sphere{p}` /
single `Using Point{p}` center constraint.
"""
function add_ruled_surface!(m::GeoModel, loops; tag::Integer=0,
                            sphere_center=nothing)
    caller="add_ruled_surface!"
    ids=Int[_tag(ℓ,caller,2) for ℓ in loops]
    isempty(ids) && throw(ArgumentError(
        "$caller: a surface requires at least one curve loop"))
    for id in ids
        haskey(m.loops,id) || throw(ArgumentError(
            "$caller: unknown Loop[$id]"))
    end
    borders=length(m.loops[first(ids)])
    borders in (3,4) || throw(ArgumentError(
        "$caller: wrong definition of surface: $borders border curves " *
        "instead of 3 or 4"))
    sc=sphere_center===nothing ? nothing : _tag(sphere_center,caller,0)
    sc!==nothing && (haskey(m.points,sc) || throw(ArgumentError(
        "$caller: unknown sphere center Point[$sc]")))
    t=_alloc_tag!(m,2,_tag(tag,caller,2),caller)
    (haskey(m.surfaces,t) || haskey(m.discrete,(2,t))) && throw(ArgumentError(
        "$caller: Surface[$t] already exists"))
    m.surfaces[t]=ids
    m.surface_types[t]=borders==4 ? :ruled : :tric
    sc!==nothing && (m.surface_geometry[t]=(sphere_center=sc,))
    return t
end

# ── Curve loop ordering ──────────────────────────────────────────────────────

# Gmsh's `SortEdgesInLoop` (reorient=false): reorder the signed curve tags into
# connectivity order, starting a new subloop with the next remaining curve
# whenever the chain closes. A chain that dead-ends before consuming every
# curve is an error; closure itself is not required (verified at mesh time).
function _sort_curve_loop(m::GeoModel, ids::Vector{Int},
                          caller::AbstractString)
    isempty(ids) && return Int[]
    remaining=collect(ids)
    ordered=Int[]
    ends(id)=(let (a,b)=m.curves[abs(id)]; id>0 ? (a,b) : (b,a) end)
    c0=c1=first(remaining)
    push!(ordered,c0); deleteat!(remaining,1)
    nb=length(ids)
    while length(ordered)<nb
        advanced=false
        for (i,c2) in pairs(remaining)
            ends(c1)[2]==ends(c2)[1] || continue
            push!(ordered,c2); deleteat!(remaining,i)
            c1=c2
            if ends(c1)[2]==ends(c0)[1] && !isempty(remaining)
                c0=c1=first(remaining)
                push!(ordered,c0); deleteat!(remaining,1)
            end
            advanced=true
            break
        end
        advanced || throw(ArgumentError("$caller: curve loop is wrong"))
    end
    return ordered
end

# A loop is closed for meshing when its signed curves form closed chains:
# the multiset of oriented end points equals the multiset of start points.
function _verify_loop_closed(m::GeoModel, loop::Int,
                             caller::AbstractString,
                             owner::AbstractString)
    starts=Dict{Int,Int}(); ends=Dict{Int,Int}()
    for signed in m.loops[loop]
        a,b=m.curves[abs(signed)]
        s,e=signed>0 ? (a,b) : (b,a)
        starts[s]=get(starts,s,0)+1
        ends[e]=get(ends,e,0)+1
    end
    starts==ends || throw(ArgumentError(
        "$caller: Loop[$loop] is not closed on $owner"))
    return nothing
end

# ── OCC-style primitive entities ─────────────────────────────────────────────
#
# Gmsh's solid primitives live in its OpenCASCADE kernel, which materializes an
# explicit boundary representation: closed-circle edges sharing a single seam
# vertex, arc-length-parametrized seam lines, degenerate point edges at poles
# and apices, and typed Cylinder/Sphere/Cone surfaces. Tessella materializes
# the same entity layout for `add_cylinder!`/`add_sphere!`/`add_cone!` while
# keeping the analytic encodings that drive native meshing.
#
# OCC curves carry their parametrization in `curve_geometry` (distinct from the
# `(n=,)` records of built-in arcs):
#   (occ=:circle,      center, n, X, r, t0, t1) — p(t)=center+r(cos t·X+sin t·Y),
#                                               Y=n×X; t0..t1 = 0..2π for closed
#                                               edges, the OCC trim for meridians
#   (occ=:line,        t0, t1)                  — arc-length parameter range
#   (occ=:degenerate,  t0, t1)                  — collapsed pole/apex edge
# The circle center is stored as coordinates, not a Point entity: OCC exposes
# only the rim vertices, so centers must not appear in the dim-0 entity list.
# Stored geometry is authoritative; transforms rewrite it, and the
# `_reconcile_curved_encodings!` pass drops a volume's compact encoding when
# independently moved rim vertices no longer satisfy the encoded circles.
#
# OCC surfaces carry matching analytic records in `surface_geometry`:
#   (occ=:cylinder, center, axis, X, radius, height)
#   (occ=:sphere,   center, radius, axis, X)
#   (occ=:cone,     center, axis, X, r1, r2, height)
#   (occ=:torus,    center, axis, X, r1, r2, angle)
# `axis`/`X` are the OCC construction frame; the second in-plane direction is
# always axis×X. They parametrize `model_value`/`model_parametrization_bounds`
# for the face: u sweeps around `axis`, v is the profile direction (arc length
# for Cylinder/Cone, latitude for Sphere, tube angle for Torus).

const _OCC_TWO_PI = 2π

# The XDirection OCC assigns to `gp_Ax2(origin, axis)` (same heuristic as
# `gp_Pln(P,V)`): the in-plane direction lying in the coordinate plane of the
# axis's smallest component, signed by the other two components so the frame
# stays continuous across quadrant changes.
function _occ_reference_direction(axis::NTuple{3,Float64})
    a,b,c=axis
    aa,bb,cc=abs(a),abs(b),abs(c)
    x=if bb<=aa && bb<=cc
        aa>cc ? (-c,0.0,a) : (c,0.0,-a)
    elseif aa<=bb && aa<=cc
        bb>cc ? (0.0,-c,b) : (0.0,c,-b)
    else
        aa>bb ? (-b,a,0.0) : (b,-a,0.0)
    end
    l=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
    l>0 || throw(ErrorException(
        "_occ_reference_direction: degenerate axis $axis"))
    return (x[1]/l,x[2]/l,x[3]/l)
end

# Second in-plane direction of an OCC circle record.
@inline _occ_circle_y(g) = _arc_cross(g.n,g.X)

function _occ_circle_point(g,t::Float64)
    Y=_occ_circle_y(g)
    c,s=cos(t),sin(t)
    C,X=g.center,g.X
    return (C[1]+g.r*(c*X[1]+s*Y[1]),
            C[2]+g.r*(c*X[2]+s*Y[2]),
            C[3]+g.r*(c*X[3]+s*Y[3]))
end

function _occ_circle_derivative(g,t::Float64)
    Y=_occ_circle_y(g)
    c,s=cos(t),sin(t)
    X=g.X
    return (g.r*(-s*X[1]+c*Y[1]),
            g.r*(-s*X[2]+c*Y[2]),
            g.r*(-s*X[3]+c*Y[3]))
end

function _occ_circle_second_derivative(g,t::Float64)
    Y=_occ_circle_y(g)
    c,s=cos(t),sin(t)
    X=g.X
    return (-g.r*(c*X[1]+s*Y[1]),
            -g.r*(c*X[2]+s*Y[2]),
            -g.r*(c*X[3]+s*Y[3]))
end

# Exact OCC-circle bounding box. A full circle's extent along axis i is
# r·sqrt(1−nᵢ²); a trimmed range adds its in-range axis extrema to the
# endpoint box.
function _occ_circle_bounding_box(g)
    if g.t1-g.t0>=_OCC_TWO_PI-1e-12
        C=g.center
        return (ntuple(i->C[i]-g.r*sqrt(max(0.0,1.0-g.n[i]*g.n[i])),3)...,
                ntuple(i->C[i]+g.r*sqrt(max(0.0,1.0-g.n[i]*g.n[i])),3)...)
    end
    p0=_occ_circle_point(g,g.t0); p1=_occ_circle_point(g,g.t1)
    lo=[min(p0[i],p1[i]) for i in 1:3]
    hi=[max(p0[i],p1[i]) for i in 1:3]
    Y=_occ_circle_y(g)
    for axis in 1:3
        (g.X[axis]==0.0 && Y[axis]==0.0) && continue
        base=atan(Y[axis],g.X[axis])
        for k in -4:4
            t=base+k*π
            (g.t0-1e-12<=t<=g.t1+1e-12) || continue
            p=_occ_circle_point(g,clamp(t,g.t0,g.t1))
            for i in 1:3
                lo[i]=min(lo[i],p[i]); hi[i]=max(hi[i],p[i])
            end
        end
    end
    return (lo[1],lo[2],lo[3],hi[1],hi[2],hi[3])
end

# Inverse parametrization: the curve parameter of the closest point, mapped
# into the stored [t0,t1] interval.
function _occ_circle_parameter(g,point::NTuple{3,Float64})
    Y=_occ_circle_y(g)
    rel=_arc_sub(point,g.center)
    θ=atan(_arc_dot(rel,Y),_arc_dot(rel,g.X))
    span=_OCC_TWO_PI
    θ=θ-floor((θ-g.t0)/span)*span
    return θ
end

# The OCC record of a curve, or `nothing` for built-in entities.
_occ_geometry(m::GeoModel,tag::Int) =
    _occ_geometry(get(m.curve_geometry,tag,nothing))
_occ_geometry(g::NamedTuple) = hasproperty(g,:occ) ? g : nothing
_occ_geometry(::Nothing) = nothing

# Whether the curve's endpoints still satisfy its OCC record. Materialized
# geometry is authoritative: a vertex moved independently of its edge leaves a
# stale record that must fail queries rather than answer with old geometry.
function _occ_endpoints_satisfied(m::GeoModel, tag::Int, g)
    a,b=m.curves[tag]
    pa,pb=m.points[a],m.points[b]
    scale=max(1.0,maximum(abs,pa),maximum(abs,pb))
    tol=1e-9*scale
    if g.occ===:circle
        a==b && begin
            rel=_arc_sub(pa,g.center)
            return abs(sqrt(_arc_dot(rel,rel))-g.r)<=tol &&
                   abs(_arc_dot(rel,g.n))<=tol
        end
        return _points_close(_occ_circle_point(g,g.t0),pa,tol) &&
               _points_close(_occ_circle_point(g,g.t1),pb,tol)
    elseif g.occ===:line
        chord=_arc_sub(pb,pa)
        return abs(sqrt(_arc_dot(chord,chord))-(g.t1-g.t0))<=tol
    end
    return true
end

# The OCC record of a curve for query paths — `nothing` for built-in entities,
# an explicit failure when the record no longer matches its endpoints.
function _occ_geometry_checked(m::GeoModel, tag::Int, caller::AbstractString)
    g=_occ_geometry(m,tag)
    g===nothing && return nothing
    _occ_endpoints_satisfied(m,tag,g) || throw(ErrorException(
        "$caller: Curve[$tag] geometry is inconsistent with its endpoints; " *
        "rebuild the model"))
    return g
end

# Rim samples standing in for the vertices an OCC circle collapses: a plane
# bounded by a single closed circle still yields three non-collinear points
# for the exact fit. Returns (coordinates, stand-in point tags) pairs.
function _occ_surface_samples(m::GeoModel, surface::Int,
                              caller::AbstractString)
    samples=NTuple{3,Float64}[]
    tags=Int[]
    for loop in m.surfaces[surface], signed in m.loops[loop]
        curve=abs(signed)
        occ=_occ_geometry_checked(m,curve,caller)
        (occ===nothing || occ.occ!==:circle) && continue
        a,_=m.curves[curve]
        for k in 0:2
            push!(samples,_occ_circle_point(
                occ,occ.t0+(occ.t1-occ.t0)*(k/3)))
            push!(tags,a)
        end
    end
    return samples,tags
end

# ── OCC entity construction ──────────────────────────────────────────────────

# A closed circle edge (`[rim,rim]` endpoints) or OCC-trimmed circle. No public
# equivalent exists: callers materialize whole solids and roll back as a unit.
function _add_occ_circle!(m::GeoModel, p1::Int, p2::Int,
                         center::NTuple{3,Float64}, n::NTuple{3,Float64},
                         X::NTuple{3,Float64}, r::Float64,
                         t0::Float64, t1::Float64)
    t=_alloc_tag!(m,1,0,"_add_occ_circle!")
    m.curves[t]=(p1,p2)
    m.curve_types[t]=:circle
    m.curve_geometry[t]=(occ=:circle,center=center,n=n,X=X,r=r,t0=t0,t1=t1)
    return t
end

# A seam line with an OCC arc-length parameter range [0,length].
function _add_occ_line!(m::GeoModel, p1::Int, p2::Int)
    t=_alloc_tag!(m,1,0,"_add_occ_line!")
    a,b=m.points[p1],m.points[p2]
    len=sqrt((b[1]-a[1])^2+(b[2]-a[2])^2+(b[3]-a[3])^2)
    m.curves[t]=(p1,p2)
    m.curve_geometry[t]=(occ=:line,t0=0.0,t1=len)
    return t
end

# A degenerate OCC edge: a zero-length edge collapsed on a pole or apex vertex.
function _add_occ_degenerate!(m::GeoModel, p::Int)
    t=_alloc_tag!(m,1,0,"_add_occ_degenerate!")
    m.curves[t]=(p,p)
    m.curve_types[t]=:degenerate
    m.curve_geometry[t]=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI)
    return t
end

# A typed OCC face: `loops` are existing curve-loop tags, `kind` one of
# :cylinder/:sphere/:cone/:torus, `geometry` the analytic record.
function _add_occ_surface!(m::GeoModel, kind::Symbol, loops::Vector{Int},
                           geometry)
    t=_alloc_tag!(m,2,0,"_add_occ_surface!")
    m.surfaces[t]=loops
    m.surface_types[t]=kind
    m.surface_geometry[t]=geometry
    return t
end

# ── Materialization consistency ──────────────────────────────────────────────
#
# A materialized primitive's compact encoding is only valid while its boundary
# entities still describe the encoded solid. These predicates re-derive that
# satisfaction from the OCC geometry records and rim/pole vertices — shared by
# `_reconcile_curved_encodings!` (drops a stale encoding after independent
# sub-entity moves) and `_model_volume_bounds` (rejects corrupted models).

# The volume's boundary curves classified by OCC role, deduplicated — seam and
# meridian edges occur twice in their periodic face's wire.
function _materialized_boundary_records(m::GeoModel, tag::Int)
    circles=Tuple{Int,NamedTuple}[]
    seams=Int[]
    degenerates=Int[]
    seen=Set{Int}()
    for shell in m.volumes[tag], ssurf in m.surface_loops[shell],
            loop in m.surfaces[abs(ssurf)], signed in m.loops[loop]
        curve=abs(signed)
        curve in seen && continue
        push!(seen,curve)
        occ=_occ_geometry(m,curve)
        occ===nothing && continue
        if occ.occ===:circle
            push!(circles,(curve,occ))
        elseif occ.occ===:line
            push!(seams,curve)
        else
            a,_=m.curves[curve]
            push!(degenerates,a)
        end
    end
    return circles,seams,degenerates
end

# A circle edge satisfies its expected (center, n, r) when the stored record
# matches and its rim vertex still lies on the circle.
function _occ_circle_consistent(m::GeoModel, curve::Int, g, center, n,
                                r::Float64, tol::Float64)
    _points_close(g.center,center,tol) || return false
    _points_close(g.n,n,tol) || return false
    abs(g.r-r)>tol && return false
    a,_=m.curves[curve]
    rim=_arc_sub(m.points[a],g.center)
    radial=sqrt(_arc_dot(rim,rim))
    abs(radial-r)>tol && return false
    abs(_arc_dot(rim,n))>tol && return false
    return true
end

function _materialized_cylinder_consistent(m::GeoModel, tag::Int, rec)
    h=rec.height
    n=(rec.axis[1]/h,rec.axis[2]/h,rec.axis[3]/h)
    top=rec.center .+ rec.axis
    scale=max(1.0,rec.radius,
              maximum(abs,rec.center),maximum(abs,rec.axis))
    tol=1e-9*scale
    circles,seams,degenerates=_materialized_boundary_records(m,tag)
    (length(circles)==2 && length(seams)==1 && isempty(degenerates)) ||
        return false
    ends_found=falses(2)
    for (curve,g) in circles
        if _points_close(g.center,rec.center,tol)
            ends_found[1] && return false
            ends_found[1]=true
            _occ_circle_consistent(m,curve,g,rec.center,n,rec.radius,tol) ||
                return false
        elseif _points_close(g.center,top,tol)
            ends_found[2] && return false
            ends_found[2]=true
            _occ_circle_consistent(m,curve,g,top,n,rec.radius,tol) ||
                return false
        else
            return false
        end
    end
    return all(ends_found)
end

function _materialized_sphere_consistent(m::GeoModel, tag::Int, rec)
    scale=max(1.0,rec.radius,maximum(abs,rec.center))
    tol=1e-9*scale
    circles,seams,degenerates=_materialized_boundary_records(m,tag)
    (length(circles)==1 && isempty(seams) && length(degenerates)==2) ||
        return false
    curve,g=first(circles)
    _points_close(g.center,rec.center,tol) || return false
    abs(g.r-rec.radius)>tol && return false
    north=(rec.center[1],rec.center[2],rec.center[3]+rec.radius)
    south=(rec.center[1],rec.center[2],rec.center[3]-rec.radius)
    found_north=found_south=false
    for p in degenerates
        point=m.points[p]
        if _points_close(point,north,tol)
            found_north && return false
            found_north=true
        elseif _points_close(point,south,tol)
            found_south && return false
            found_south=true
        else
            return false
        end
    end
    return found_north && found_south
end

function _materialized_cone_consistent(m::GeoModel, tag::Int, rec)
    h=rec.height
    n=(rec.axis[1]/h,rec.axis[2]/h,rec.axis[3]/h)
    top=rec.center .+ rec.axis
    scale=max(1.0,rec.r1,rec.r2,
              maximum(abs,rec.center),maximum(abs,rec.axis))
    tol=1e-9*scale
    circles,seams,degenerates=_materialized_boundary_records(m,tag)
    (length(seams)==1 &&
     length(circles)==(rec.r1>0)+(rec.r2>0) &&
     length(degenerates)==(rec.r1==0)+(rec.r2==0)) || return false
    top_found=bottom_found=false
    for (curve,g) in circles
        if rec.r2>0 && _points_close(g.center,top,tol)
            top_found && return false
            top_found=true
            _occ_circle_consistent(m,curve,g,top,n,rec.r2,tol) || return false
        elseif rec.r1>0 && _points_close(g.center,rec.center,tol)
            bottom_found && return false
            bottom_found=true
            _occ_circle_consistent(m,curve,g,rec.center,n,rec.r1,tol) ||
                return false
        else
            return false
        end
    end
    apex_found=falses(2)
    for p in degenerates
        point=m.points[p]
        if rec.r2==0 && _points_close(point,top,tol)
            apex_found[1] && return false
            apex_found[1]=true
        elseif rec.r1==0 && _points_close(point,rec.center,tol)
            apex_found[2] && return false
            apex_found[2]=true
        else
            return false
        end
    end
    return top_found==(rec.r2>0) && bottom_found==(rec.r1>0) &&
           apex_found[1]==(rec.r2==0) && apex_found[2]==(rec.r1==0)
end

# Whether a materialized primitive volume's compact encoding still describes
# its boundary entities.
function _materialized_curved_consistent(m::GeoModel, tag::Int)
    haskey(m.cylinders,tag) &&
        return _materialized_cylinder_consistent(m,tag,m.cylinders[tag])
    haskey(m.spheres,tag) &&
        return _materialized_sphere_consistent(m,tag,m.spheres[tag])
    haskey(m.cones,tag) &&
        return _materialized_cone_consistent(m,tag,m.cones[tag])
    return false
end

# Roll back the partially materialized entity set on a mid-build failure,
# mirroring `add_box!`'s transactional construction.
function _occ_materialize_rollback!(m::GeoModel, points, curves, loops,
                                    surfaces, shell)
    shell!=0 && delete!(m.surface_loops,shell)
    for surface in surfaces
        delete!(m.surfaces,surface)
        delete!(m.surface_types,surface)
        delete!(m.surface_geometry,surface)
    end
    for loop in loops
        delete!(m.loops,loop)
    end
    for curve in curves
        delete!(m.curves,curve)
        delete!(m.curve_types,curve)
        delete!(m.curve_geometry,curve)
        delete!(m.curve_control_points,curve)
    end
    for point in points
        delete!(m.points,point)
        delete!(m.point_size,point)
    end
    return nothing
end

# OCC's cylinder layout: top rim Point, bottom rim Point; top Circle, seam
# Line, bottom Circle; Cylinder face ([-top,-seam,bottom,seam]), Plane caps;
# shell [lateral, +top, -bottom]. Both circles wind CCW about +axis from the
# OCC `gp_Ax2` reference direction X (`_occ_reference_direction`).
function _materialize_cylinder!(m::GeoModel, base::NTuple{3,Float64},
                                axis::NTuple{3,Float64}, r::Float64,
                                h::Float64)
    n=(axis[1]/h,axis[2]/h,axis[3]/h)
    X=_occ_reference_direction(n)
    top=(base[1]+axis[1],base[2]+axis[2],base[3]+axis[3])
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        push!(points,add_point!(m,top[1]+r*X[1],top[2]+r*X[2],top[3]+r*X[3]))
        push!(points,add_point!(m,base[1]+r*X[1],base[2]+r*X[2],base[3]+r*X[3]))
        for p in points
            delete!(m.point_size,p)
        end
        p_top,p_bot=points
        push!(curves,_add_occ_circle!(m,p_top,p_top,top,n,X,r,0.0,_OCC_TWO_PI))
        push!(curves,_add_occ_line!(m,p_bot,p_top))
        push!(curves,_add_occ_circle!(m,p_bot,p_bot,base,n,X,r,0.0,_OCC_TWO_PI))
        c_top,c_seam,c_bot=curves
        push!(loops,add_curve_loop!(m,[-c_top,-c_seam,c_bot,c_seam]))
        push!(surfaces,_add_occ_surface!(m,:cylinder,[last(loops)],
              (occ=:cylinder,center=base,axis=n,X=X,radius=r,height=h)))
        push!(loops,add_curve_loop!(m,[c_top]))
        push!(surfaces,add_plane_surface!(m,[last(loops)]))
        push!(loops,add_curve_loop!(m,[c_bot]))
        push!(surfaces,add_plane_surface!(m,[last(loops)]))
        s_lat,s_top,s_bot=surfaces
        shell=add_surface_loop!(m,[s_lat,s_top,-s_bot])
    catch
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,shell)
        rethrow()
    end
    return shell
end

# OCC's sphere layout: north/south pole Points; a degenerate edge on each pole
# and a meridian Circle trimmed to [3π/2,5π/2] in the xz-plane through +x̂;
# one Sphere face ([-degN,-meridian,degS,meridian]); shell [face].
function _materialize_sphere!(m::GeoModel, center::NTuple{3,Float64},
                              r::Float64)
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        push!(points,add_point!(m,center[1],center[2],center[3]+r))
        push!(points,add_point!(m,center[1],center[2],center[3]-r))
        for p in points
            delete!(m.point_size,p)
        end
        p_n,p_s=points
        push!(curves,_add_occ_degenerate!(m,p_n))
        # The OCC meridian frame is X=+x̂, Y=+ẑ (n=X×Y=-ŷ).
        push!(curves,_add_occ_circle!(m,p_s,p_n,center,(0.0,-1.0,0.0),
              (1.0,0.0,0.0),r,1.5π,2.5π))
        push!(curves,_add_occ_degenerate!(m,p_s))
        c_n,c_mer,c_s=curves
        push!(loops,add_curve_loop!(m,[-c_n,-c_mer,c_s,c_mer]))
        push!(surfaces,_add_occ_surface!(m,:sphere,[last(loops)],
              (occ=:sphere,center=center,radius=r,
               axis=(0.0,0.0,1.0),X=(1.0,0.0,0.0))))
        shell=add_surface_loop!(m,[last(surfaces)])
    catch
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,shell)
        rethrow()
    end
    return shell
end

# OCC's cone layout mirrors the cylinder's, with a degenerate edge replacing
# whichever end has zero radius and only the non-degenerate end receiving a
# Plane cap. Shells: [lateral,+top,-bottom] / [lateral,+top] (r1=0) /
# [lateral,-bottom] (r2=0).
function _materialize_cone!(m::GeoModel, base::NTuple{3,Float64},
                          axis::NTuple{3,Float64}, r1::Float64, r2::Float64,
                          h::Float64)
    n=(axis[1]/h,axis[2]/h,axis[3]/h)
    X=_occ_reference_direction(n)
    top=(base[1]+axis[1],base[2]+axis[2],base[3]+axis[3])
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        push!(points,add_point!(m,top[1]+r2*X[1],top[2]+r2*X[2],top[3]+r2*X[3]))
        push!(points,add_point!(m,base[1]+r1*X[1],base[2]+r1*X[2],base[3]+r1*X[3]))
        for p in points
            delete!(m.point_size,p)
        end
        p_top,p_bot=points
        push!(curves,r2>0 ?
              _add_occ_circle!(m,p_top,p_top,top,n,X,r2,0.0,_OCC_TWO_PI) :
              _add_occ_degenerate!(m,p_top))
        push!(curves,_add_occ_line!(m,p_bot,p_top))
        push!(curves,r1>0 ?
              _add_occ_circle!(m,p_bot,p_bot,base,n,X,r1,0.0,_OCC_TWO_PI) :
              _add_occ_degenerate!(m,p_bot))
        c_top,c_seam,c_bot=curves
        push!(loops,add_curve_loop!(m,[-c_top,-c_seam,c_bot,c_seam]))
        push!(surfaces,_add_occ_surface!(m,:cone,[last(loops)],
              (occ=:cone,center=base,axis=n,X=X,r1=r1,r2=r2,height=h)))
        s_lat=last(surfaces)
        shell_signs=Int[s_lat]
        if r2>0
            push!(loops,add_curve_loop!(m,[c_top]))
            push!(surfaces,add_plane_surface!(m,[last(loops)]))
            push!(shell_signs,last(surfaces))
        end
        if r1>0
            push!(loops,add_curve_loop!(m,[c_bot]))
            push!(surfaces,add_plane_surface!(m,[last(loops)]))
            push!(shell_signs,-last(surfaces))
        end
        shell=add_surface_loop!(m,shell_signs)
    catch
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,shell)
        rethrow()
    end
    return shell
end

# OCC's torus layout (always the z axis — `gmsh.model.occ.addTorus` and the
# `.geo` `Torus` statement share `BRepPrimAPI_MakeTorus`). A full torus
# (angle == 2π) is a single rim Point at `center + (r1+r2)·x̂`, an outer-equator
# Circle of radius r1+r2 and a meridian Circle of radius r2 both closed on it,
# and one Torus face wired `[-equator, +meridian, +equator, -meridian]`;
# shell [face]. A partial torus materializes the sweep's end vertex first,
# then the start vertex, a trimmed equator arc `start→end` over `[0,angle]`,
# the closed end and start meridian circles, the Torus face
# `[-arc, -start_meridian, +arc, +end_meridian]`, and Plane caps on each end;
# shell [torus, +start_cap, -end_cap].
function _materialize_torus!(m::GeoModel, center::NTuple{3,Float64},
                             r1::Float64, r2::Float64, angle::Float64)
    n=(0.0,0.0,1.0); X=(1.0,0.0,0.0); Y=(0.0,1.0,0.0)
    ca,sa=cos(angle),sin(angle)
    Xe=(ca*X[1]+sa*Y[1],ca*X[2]+sa*Y[2],ca*X[3]+sa*Y[3])
    # A meridian's plane contains the axis: X_dir = R(u)·X, Y_dir = axis, so
    # its normal is R(u)·X × axis.
    meridian_normal=(dir)->_arc_cross(dir,n)
    full=angle>=_OCC_TWO_PI
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        if full
            push!(points,add_point!(m,
                center[1]+(r1+r2)*X[1],center[2]+(r1+r2)*X[2],
                center[3]+(r1+r2)*X[3]))
        else
            push!(points,add_point!(m,
                center[1]+(r1+r2)*Xe[1],center[2]+(r1+r2)*Xe[2],
                center[3]+(r1+r2)*Xe[3]))
            push!(points,add_point!(m,
                center[1]+(r1+r2)*X[1],center[2]+(r1+r2)*X[2],
                center[3]+(r1+r2)*X[3]))
        end
        for p in points
            delete!(m.point_size,p)
        end
        if full
            p_rim=points[1]
            push!(curves,_add_occ_circle!(m,p_rim,p_rim,center,n,X,r1+r2,
                  0.0,_OCC_TWO_PI))
            push!(curves,_add_occ_circle!(m,p_rim,p_rim,
                  (center[1]+r1*X[1],center[2]+r1*X[2],center[3]+r1*X[3]),
                  meridian_normal(X),X,r2,0.0,_OCC_TWO_PI))
            c_eq,c_mer=curves
            push!(loops,add_curve_loop!(m,[-c_eq,c_mer,c_eq,-c_mer]))
            push!(surfaces,_add_occ_surface!(m,:torus,[last(loops)],
                  (occ=:torus,center=center,axis=n,X=X,r1=r1,r2=r2,
                   angle=angle)))
            shell=add_surface_loop!(m,[last(surfaces)])
        else
            p_end,p_start=points
            push!(curves,_add_occ_circle!(m,p_start,p_end,center,n,X,r1+r2,
                  0.0,angle))
            push!(curves,_add_occ_circle!(m,p_end,p_end,
                  (center[1]+r1*Xe[1],center[2]+r1*Xe[2],center[3]+r1*Xe[3]),
                  meridian_normal(Xe),Xe,r2,0.0,_OCC_TWO_PI))
            push!(curves,_add_occ_circle!(m,p_start,p_start,
                  (center[1]+r1*X[1],center[2]+r1*X[2],center[3]+r1*X[3]),
                  meridian_normal(X),X,r2,0.0,_OCC_TWO_PI))
            c_arc,c_end,c_start=curves
            push!(loops,add_curve_loop!(m,[-c_arc,-c_start,c_arc,c_end]))
            push!(surfaces,_add_occ_surface!(m,:torus,[last(loops)],
                  (occ=:torus,center=center,axis=n,X=X,r1=r1,r2=r2,
                   angle=angle)))
            s_torus=last(surfaces)
            push!(loops,add_curve_loop!(m,[c_start]))
            push!(surfaces,add_plane_surface!(m,[last(loops)]))
            s_start=last(surfaces)
            push!(loops,add_curve_loop!(m,[c_end]))
            push!(surfaces,add_plane_surface!(m,[last(loops)]))
            s_end=last(surfaces)
            shell=add_surface_loop!(m,[s_torus,s_start,-s_end])
        end
    catch
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,shell)
        rethrow()
    end
    return shell
end

# ── OCC surface evaluation ───────────────────────────────────────────────────
#
# `model_value`/`model_parametrization_bounds` on an OCC face read the stored
# analytic record, matching the parametrizations OCC reports for
# BRepPrimAPI-made primitives: u is the azimuth about `axis` and v the profile
# parameter (arc length along the axis for Cylinder, slant length for Cone,
# latitude for Sphere, tube angle for Torus).

@inline _occ_frame_y(g) = _arc_cross(g.axis,g.X)

function _occ_surface_point(g, u::Float64, v::Float64)
    Y=_occ_frame_y(g)
    c,s=cos(u),sin(u)
    if g.occ===:cylinder
        return (g.center[1]+g.radius*(c*g.X[1]+s*Y[1])+v*g.axis[1],
                g.center[2]+g.radius*(c*g.X[2]+s*Y[2])+v*g.axis[2],
                g.center[3]+g.radius*(c*g.X[3]+s*Y[3])+v*g.axis[3])
    elseif g.occ===:sphere
        cv,sv=cos(v),sin(v)
        return (g.center[1]+g.radius*(cv*(c*g.X[1]+s*Y[1])+sv*g.axis[1]),
                g.center[2]+g.radius*(cv*(c*g.X[2]+s*Y[2])+sv*g.axis[2]),
                g.center[3]+g.radius*(cv*(c*g.X[3]+s*Y[3])+sv*g.axis[3]))
    elseif g.occ===:cone
        slant=sqrt(g.height*g.height+(g.r1-g.r2)*(g.r1-g.r2))
        r=g.r1-v*(g.r1-g.r2)/slant
        z=v*g.height/slant
        return (g.center[1]+r*(c*g.X[1]+s*Y[1])+z*g.axis[1],
                g.center[2]+r*(c*g.X[2]+s*Y[2])+z*g.axis[2],
                g.center[3]+r*(c*g.X[3]+s*Y[3])+z*g.axis[3])
    elseif g.occ===:torus
        cv,sv=cos(v),sin(v)
        w=g.r1+g.r2*cv
        return (g.center[1]+w*(c*g.X[1]+s*Y[1])+g.r2*sv*g.axis[1],
                g.center[2]+w*(c*g.X[2]+s*Y[2])+g.r2*sv*g.axis[2],
                g.center[3]+w*(c*g.X[3]+s*Y[3])+g.r2*sv*g.axis[3])
    end
    throw(ArgumentError(
        "_occ_surface_point: unsupported OCC surface kind $(g.occ)"))
end

# OCC-reported parametrization bounds per surface kind.
function _occ_surface_bounds(g)
    if g.occ===:cylinder
        return [0.0,0.0],[2π,g.height]
    elseif g.occ===:sphere
        return [0.0,-π/2],[2π,π/2]
    elseif g.occ===:cone
        return [0.0,0.0],
               [2π,sqrt(g.height*g.height+(g.r1-g.r2)*(g.r1-g.r2))]
    elseif g.occ===:torus
        return [0.0,0.0],[g.angle,2π]
    end
    throw(ArgumentError(
        "_occ_surface_bounds: unsupported OCC surface kind $(g.occ)"))
end

# Exact axis-aligned box of a (possibly partial) torus face. Along coordinate
# axis i the surface coordinate is
#   C_i + r1·w_i(u) + r2·(w_i(u)·cos v + n_i·sin v),  w_i(u) = ρ_i·cos(u−φ_i)
# with ρ_i = √(1−n_i²) the in-plane projection of axis i and φ_i the azimuth
# that aligns R(u)·X with axis i. Maximizing over v at fixed u leaves the
# one-variable function f(θ) = r1·ρ·cosθ + r2·√(1−ρ²·sin²θ) on
# θ ∈ [−φ_i, angle−φ_i]; the lower box side is the same maximization shifted
# by π. Interior extrema sit at θ = kπ or at the spindle-torus stationary
# points given by the closed-form root below.
function _occ_torus_bounding_box(g)
    Y=_occ_frame_y(g)
    lo=Vector{Float64}(undef,3); hi=Vector{Float64}(undef,3)
    for i in 1:3
        n_i=g.axis[i]
        rho2=max(0.0,1.0-n_i*n_i)
        rho=sqrt(rho2)
        if rho<=1e-15
            lo[i]=g.center[i]-g.r2*abs(n_i)
            hi[i]=g.center[i]+g.r2*abs(n_i)
            continue
        end
        phi=atan(Y[i],g.X[i])
        fmax=_occ_torus_axis_max(g.r1,g.r2,rho,-phi,g.angle-phi)
        fmin=_occ_torus_axis_max(g.r1,g.r2,rho,-phi+π,g.angle-phi+π)
        lo[i]=g.center[i]-fmin
        hi[i]=g.center[i]+fmax
    end
    return (lo[1],lo[2],lo[3],hi[1],hi[2],hi[3])
end

# max over θ∈[a,b] of r1·ρ·cosθ + r2·√(1−ρ²·sin²θ). The candidate set is
# exact: interval endpoints, every kπ in the range (θ=0 is the outer-equator
# maximum, θ=π the opposite-side extremum), and the stationary roots
# cosθ = −r1·√(1−ρ²·s)/(r2·ρ), sin²θ = s = (r1²−r2²ρ²)/(ρ²(r1²−r2²)), which
# only exist for self-intersecting spindles (r1 < r2, ρ ≤ r1/r2).
function _occ_torus_axis_max(r1::Float64, r2::Float64, rho::Float64,
                             a::Float64, b::Float64)
    f(θ)=r1*rho*cos(θ)+r2*sqrt(max(0.0,1.0-rho*rho*sin(θ)*sin(θ)))
    best=max(f(a),f(b))
    lo_k=ceil(Int,a/π-1e-12); hi_k=floor(Int,b/π+1e-12)
    for k in lo_k:hi_k
        best=max(best,f(k*π))
    end
    if r1!=r2
        s=(r1*r1-r2*r2*rho*rho)/(rho*rho*(r1*r1-r2*r2))
        if 0.0<=s<1.0
            root=acos(-sqrt(1.0-s))
            for base in (root,-root), k in -1:1
                θ=base+2k*π
                a-1e-12<=θ<=b+1e-12 || continue
                best=max(best,f(clamp(θ,a,b)))
            end
        end
    end
    return best
end
