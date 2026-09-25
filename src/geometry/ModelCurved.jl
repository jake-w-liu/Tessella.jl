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

# `gmshEdge::geomType()` — all four spline-family records report `Nurb`.
const _CURVE_TYPE_NAMES = Dict{Symbol,String}(
    :line=>"Line", :circle=>"Circle", :ellipse=>"Ellipse",
    :spline=>"Nurb", :bspline=>"Nurb", :bezier=>"Nurb", :nurbs=>"Nurb",
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
_myasin(a::Float64) = a<=-1.0 ? -π/2 : a>=1.0 ? π/2 : _gm_asin(a)

# Gmsh's `myatan2`: (0,0) maps to +0.0 rather than atan2's signed zero.
_myatan2(y::Float64,x::Float64) = (y==0.0 && x==0.0) ? 0.0 : _gm_atan2(y,x)

# `norme`: normalize by the `norm3` magnitude (sqrt of the summed squares,
# not hypot); a zero vector stays zero. `fma` replicates the fused
# multiply-add that the compiled Gmsh `norm3`/`prodve`/`prosca` expressions
# produce (`x²+y²+z²` → `fma(z,z,fma(x,x,y*y))`; the first product of each
# cross/dot component fused, later products rounded).
@inline function _arc_norme(v::NTuple{3,Float64})
    n=sqrt(fma(v[1],v[1],fma(v[2],v[2],v[3]*v[3])))
    n==0.0 && return v
    inv=1.0/n
    return (v[1]*inv, v[2]*inv, v[3]*inv)
end

# `prodve` (Gmsh) writes its y component as `-a[0]*b[2] + a[2]*b[0]` while
# OCCT `gp_XYZ::Crossed` writes `z*b.x - x*b.z`; the fused multiply-add lands
# on a different operand, so the two ports stay distinct.
@inline _arc_cross(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    (fma(a[2],b[3],-(a[3]*b[2])), fma(-a[1],b[3],a[3]*b[1]),
     fma(a[1],b[2],-(a[2]*b[1])))
@inline _arc_dot(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    fma(a[3],b[3],fma(a[1],b[1],a[2]*b[2]))
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
    if sqrt(fma(n[1],n[1],fma(n[2],n[2],n[3]*n[3])))<1e-15
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
    R=sqrt(fma(x0,x0,y0*y0)); R2=sqrt(fma(x2,x2,y2*y2))
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
        s,c=_gm_sincos(A4)
        x1=fma(x0,c,y0*s); y1=fma(-x0,s,y0*c)
        xe=fma(x2,c,y2*s); ye=fma(-x2,s,y2*c)
        # sys2x2 [x1² y1²; xe² ye²]·sol = (1,1) → sol = (1/f1², 1/f2²),
        # evaluated in Gmsh's exact association: pre-rounded squared
        # products and a multiply-by-reciprocal solve. `det`/`res` contract
        # `a*b - c*d` to `fma(a,b,-(c*d))`; `pow(x,2)` in `matnorm` is not
        # contracted (the call lowers to `x*x` after the contraction pass).
        det=fma(x1*x1,ye*ye,-((xe*xe)*(y1*y1)))
        matnorm=(x1*x1)^2+(ye*ye)^2+(y1*y1)^2+(xe*xe)^2
        s0=s1=0.0
        if !(matnorm==0.0 || abs(det)/matnorm<1e-16)
            ud=1.0/det
            s0=fma(1.0,ye*ye,-((y1*y1)*1.0))*ud
            s1=fma(x1*x1,1.0,-((xe*xe)*1.0))*ud
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
    θ=fma(-(g.t1-g.t2),u,g.t1)-g.incl
    # `cos(theta)`/`sin(theta)` (and the incl pair) fuse into `__sincos_stret`
    # in the compiled `InterpolateCurve` — the fused sin can differ one ulp
    # from standalone `sin`, which the derivative FDs amplify visibly.
    sθ,cθ=_gm_sincos(θ); si,ci=_gm_sincos(g.incl)
    lx=fma(g.f1*cθ,ci,-(g.f2*sθ*si))
    ly=fma(g.f1*cθ,si,g.f2*sθ*ci)
    C,U,V,N=g.center,g.u,g.v,g.n
    # `fma` reproduces the compiled contractions: `theta` fuses the sweep
    # product, the local frame fuses the first product of each `f·trig·trig`
    # term, and `Projette` is a three-product sum (V.Pos.Z is exactly 0.0)
    # followed by the separate center add — bit-for-bit `InterpolateCurve`.
    return (fma(0.0,N[1],fma(lx,U[1],ly*V[1]))+C[1],
            fma(0.0,N[2],fma(lx,U[2],ly*V[2]))+C[2],
            fma(0.0,N[3],fma(lx,U[3],ly*V[3]))+C[3])
end

# `InterpolateCurve(der=2)` — a finite difference of the der=1 finite
# difference (`gmshEdge::secondDer`), same 1e-8 step and one-sided bounds.
function _arc_second_derivative_fd(g,u::Float64)
    eps=1e-8
    eps1=u<eps ? 0.0 : eps
    eps2=u>1.0-eps ? 0.0 : eps
    d0=_arc_first_derivative_fd(g,u-eps1)
    d1=_arc_first_derivative_fd(g,u+eps2)
    inv=1.0/(eps1+eps2)
    return ((d1[1]-d0[1])*inv,(d1[2]-d0[2])*inv,(d1[3]-d0[3])*inv)
end

# `GEdge::curvature`: norm(d1 × secondDer)·pow(1/norm(d1),3) on the
# InterpolateCurve finite differences (`norm` = sqrt of the summed squares,
# `crossprod` the SVector3 contraction — same association as `gp_XYZ::Cross`).
function _arc_curvature(g, u::Float64)
    d1=_arc_first_derivative_fd(g,u)
    d2=_arc_second_derivative_fd(g,u)
    cr=_occ_cross(d1,d2)
    return sqrt(_sqlen(cr))*_gm_pow(1.0/sqrt(_sqlen(d1)),3.0)
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
        base=_gm_atan2(a2,a1)
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
                         plane_normal=nothing, _zero_literal::Bool=false)
    caller="add_circle_arc!"
    p1=_tag(start,caller,1); pc=_tag(center,caller,1); p2=_tag(stop,caller,1)
    for p in (p1,pc,p2)
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: unknown Point[$p]"))
    end
    n=_arc_stored_normal(plane_normal,caller)
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller;literal_zero=_zero_literal)
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
                        tag::Integer=0, plane_normal=nothing,
                        _zero_literal::Bool=false)
    caller="add_ellipse_arc!"
    p1=_tag(start,caller,1); pc=_tag(center,caller,1)
    pm=_tag(major,caller,1); p2=_tag(stop,caller,1)
    for p in (p1,pc,pm,p2)
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: unknown Point[$p]"))
    end
    n=_arc_stored_normal(plane_normal,caller)
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller;literal_zero=_zero_literal)
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
                            sphere_center=nothing, _zero_literal::Bool=false)
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
    t=_alloc_tag!(m,2,_tag(tag,caller,2),caller;literal_zero=_zero_literal)
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
#   (occ=:circle,      center, n, X, Y, r, t0, t1) — p(t)=center+r(cos t·X+sin t·Y);
#                                               Y stored explicitly (transformed
#                                               frames need not satisfy Y=n×X)
#   (occ=:line,        origin, dir, t0, t1)     — p(t)=origin+t·dir, the stored
#                                               gp_Lin parametrization
#   (occ=:degenerate,  t0, t1)                  — collapsed pole/apex edge, no 3-D
#                                               curve; evaluates through the
#                                               pcurve stored on its face
# The circle center is stored as coordinates, not a Point entity: OCC exposes
# only the rim vertices, so centers must not appear in the dim-0 entity list.
# Stored geometry is authoritative; transforms rewrite it, and the
# `_reconcile_curved_encodings!` pass drops a volume's compact encoding when
# independently moved rim vertices no longer satisfy the encoded circles.
#
# OCC surfaces carry matching analytic records in `surface_geometry`:
#   (occ=:plane,    center, axis, X, Y, reversed, pcurves)
#   (occ=:cylinder, center, axis, X, Y, radius, height, pcurves)
#   (occ=:sphere,   center, radius, axis, X, Y, pcurves)
#   (occ=:cone,     center, axis, X, Y, r1, r2, height, pcurves)
#   (occ=:torus,    center, axis, X, Y, r1, r2, angle, pcurves)
# `axis`/`X`/`Y` are the OCC construction frame (a `gp_Ax3`; transforms keep
# the stored directions verbatim — under a reflection they go left-handed, so
# Y is stored, not derived). `reversed` mirrors the `TopAbs_REVERSED`
# orientation flag that `OCCFace::normal` flips for. `pcurves` maps each
# boundary edge tag to `(fwd=, rev=)` — the `Geom2d_Curve` equivalents OCC
# stores for the edge's FORWARD and (for periodic seams) REVERSED
# orientations: `(lin2d=(o,d),)` or `(circ2d=(o,x,r),)` records evaluated by
# `_occ_pcurve_value`/`_occ_surface_bounds`. They drive
# `model_reparametrize_on_surface` (OCC's native pcurve path) and the trimmed
# evaluation of degenerate edges.
# The records parametrize `model_value`/`model_parametrization_bounds` for the
# face: u sweeps around `axis`, v is the profile direction (arc length for
# Cylinder/Cone, latitude for Sphere, tube angle for Torus).

const _OCC_TWO_PI = 2π

# ── OCCT arithmetic kernel ─────────────────────────────────────────────────
#
# Bit-exact ports of the `gp_*` operations that build and transform OCC
# primitive geometry (verified against the installed OCCT 7.9.3 headers and
# compiled output). Every helper reproduces the C++ expression's floating-
# point contraction: `gp_XYZ` products fuse the first product of each
# `a·b − c·d` pair, `gp_XYZ::Multiply(gp_Mat)` contracts
# `m3·z + (m1·x + m2·y)`, and `gp_Dir`'s cross/normalize operations are
# componentwise `x/|x|`.

@inline _occ_modulus(v::NTuple{3,Float64}) =
    sqrt(fma(v[3],v[3],fma(v[1],v[1],v[2]*v[2])))
@inline _occ_dir(v::NTuple{3,Float64}) =
    (m=_occ_modulus(v); (v[1]/m,v[2]/m,v[3]/m))

# `gp_XYZ::Cross` — compiled as fma(a2,b3, −(a3·b2)) per component.
@inline _occ_cross(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    (fma(a[2],b[3],-(a[3]*b[2])), fma(a[3],b[1],-(a[1]*b[3])),
     fma(a[1],b[2],-(a[2]*b[1])))
# `gp_Dir::Crossed`/`CrossCrossed` normalize after the product.
@inline _occ_crossed(a,b) = _occ_dir(_occ_cross(a,b))

# `gp_XYZ::CrossCross`: a × (b × c) written as the direct triple-product
# expansion — t12/t23/t31 are the (b×c) components.
@inline function _occ_crosscross(a,b,c)
    t12=fma(b[1],c[2],-(b[2]*c[1]))
    t23=fma(b[2],c[3],-(b[3]*c[2]))
    t31=fma(b[3],c[1],-(b[1]*c[3]))
    (fma(a[2],t12,-(a[3]*t31)), fma(a[3],t23,-(a[1]*t12)),
     fma(a[1],t31,-(a[2]*t23)))
end
@inline _occ_crosscrossed(a,b,c) = _occ_dir(_occ_crosscross(a,b,c))

# `gp_XYZ::Multiply(gp_Mat)` — per-row `m1·x + m2·y + m3·z` contracted as
# fma(m3,z, fma(m1,x, m2·y)).
@inline _occ_matvec(m,v) =
    (fma(m[1][3],v[3],fma(m[1][1],v[1],m[1][2]*v[2])),
     fma(m[2][3],v[3],fma(m[2][1],v[1],m[2][2]*v[2])),
     fma(m[3][3],v[3],fma(m[3][1],v[1],m[3][2]*v[2])))
@inline _occ_mul(v::NTuple{3,Float64},s::Float64) = (v[1]*s,v[2]*s,v[3]*s)
@inline _add3(a::NTuple{3,Float64},b::NTuple{3,Float64}) =
    (a[1]+b[1],a[2]+b[2],a[3]+b[3])
# `gp_Dir2d` normalization (`gp_XY::Modulus`).
@inline _occ_dir2(v) =
    (m=sqrt(fma(v[1],v[1],v[2]*v[2])); (v[1]/m,v[2]/m))

# `gp_Trsf::SetRotation(Ax1(c,axis), ang)` — M = sin·K + I + (1−cos)·(aaᵀ−I)
# on the renormalized axis, translation c − M·c. Returns (M rows, loc).
function _occ_trsf_rot(c,axis,ang)
    aV=_occ_dir(axis); A,B,C=aV
    s,co=_gm_sincos(ang); om=1.0-co
    M=((0.0*s+1.0+(-C*C-B*B)*om, -C*s+(A*B)*om,      B*s+(A*C)*om),
       (    C*s+(A*B)*om,      0.0*s+1.0+(-A*A-C*C)*om, -A*s+(B*C)*om),
       (   -B*s+(A*C)*om,          A*s+(B*C)*om,    0.0*s+1.0+(-A*A-B*B)*om))
    mc=_occ_matvec(M,c)
    loc=(c[1]-mc[1],c[2]-mc[2],c[3]-mc[3])
    return M,loc
end

# `gp_Pnt::Transform`/`gp_Dir::Transform` on a rotation/reflection Trsf —
# `gp_Dir` renormalizes after the multiply (`gp_Lin`/`gp_Ax1` transforms).
_occ_xform_pnt(M,loc,p) = _add3(_occ_matvec(M,p),loc)
_occ_xform_dir(M,d) = _occ_dir(_occ_matvec(M,d))
# `gp_Dir::Rotate` — the matvec only, no renormalization.
_occ_rot_dir(M,d) = _occ_matvec(M,d)

# `gp_Ax2(P,N,Vx)` 3-arg ctor and `SetXDirection`: vx = n ^ (Vx ^ n)
# normalized, vy = n ^ vx normalized. `gp_Ax3` uses the same rule.
function _occ_ax2_setx(n,vx0)
    X=_occ_crosscrossed(n,vx0,n)
    Y=_occ_crossed(n,X)
    return X,Y
end

# `gp_Ax2(P,V)` 2-arg ctor: the reference direction lies in the coordinate
# plane of the axis's smallest component, `gp_Dir::SetCoord`-normalized,
# then `SetXDirection`.
function _occ_ax2_frame(n)
    a,b,c=n; aa,bb,cc=abs(a),abs(b),abs(c)
    d=if bb<=aa && bb<=cc
        aa>cc ? (-c,0.0,a) : (c,0.0,-a)
    elseif aa<=bb && aa<=cc
        bb>cc ? (0.0,-c,b) : (0.0,c,-b)
    else
        aa>bb ? (-b,a,0.0) : (b,-a,0.0)
    end
    D=_occ_dir(d)
    return _occ_ax2_setx(n,D)
end

# `gp_Ax2::Transform` — loc and both directions transform (dirs
# renormalized), then axis = X′ × Y′ recomputed (unlike `gp_Ax3`, which
# keeps its axis). `occ_adapt_ax2` is the BRepAdaptor identity-transform
# view: only the axis recompute is observable.
function _occ_xform_ax2(M,loc,ax2)
    l=_occ_xform_pnt(M,loc,ax2.loc)
    X=_occ_xform_dir(M,ax2.X); Y=_occ_xform_dir(M,ax2.Y)
    n=_occ_crossed(X,Y)
    return (loc=l,axis=n,X=X,Y=Y)
end
_occ_adapt_ax2(ax2) = _occ_xform_ax2(
    ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)),(0.0,0.0,0.0),ax2)

# `gp_Ax2::Rotate(Ax1, ang)` — the location rotates like `gp_Pnt::Rotate`
# (matvec + translation, same arithmetic as `Transform`), but X/Y use
# `gp_Dir::Rotate`: the raw matvec, no renormalization. The axis is then
# recomputed as `dir(X′ × Y′)`. Used by `EndFace` for sweep-end caps.
function _occ_rotate_ax2(M,loc,ax2)
    l=_occ_xform_pnt(M,loc,ax2.loc)
    X=_occ_rot_dir(M,ax2.X); Y=_occ_rot_dir(M,ax2.Y)
    n=_occ_crossed(X,Y)
    return (loc=l,axis=n,X=X,Y=Y)
end

# `BRepPrim_OneAxis` end-vertex construction: w = Loc + n·mp.Y + X·mp.X in
# the build frame, then `Transformed(rotation)` — `Added`/`Multiplied` are
# plain per-component operations.
function _occ_vertex(c,n,X,mp,M,loc)
    w=ntuple(i->(c[i]+n[i]*mp[2])+X[i]*mp[1],3)
    return _occ_xform_pnt(M,loc,w)
end

# `add_point!` minus the finite-coordinate check: OCC materializers store
# whatever the OCCT math produces — e.g. pole vertices of an oversized
# sphere overflow to ±Inf/NaN, which the real OCC vertex records also hold.
function _add_occ_point!(m::GeoModel,p::NTuple{3,Float64})
    t=_alloc_tag!(m,0,0,"_add_occ_point!")
    (haskey(m.points,t) || haskey(m.discrete,(0,t))) && throw(ArgumentError(
        "_add_occ_point!: Point[$t] already exists"))
    m.points[t]=p; m.point_size[t]=1.0
    return t
end

# `ElCLib2d::LineValue`/`CircleValue` — the 2-D meridian pcurve evals.
@inline _occ_lin2d_eval(o,d,t) = (fma(t,d[1],o[1]),fma(t,d[2],o[2]))
@inline function _occ_circ2d_eval(o,x,r,t)
    # gp_Circ2d XDirection x and derived YDirection (−x.y, x.x).
    s,c=_gm_sincos(t)
    A1,A2=r*c,r*s
    return (fma(A1,x[1],-(A2*x[2]))+o[1],fma(A1,x[2],A2*x[1])+o[2])
end

# The reference direction OCC assigns to `gp_Ax2(origin, axis)` — the
# XDirection of `_occ_ax2_frame`.
_occ_reference_direction(axis) = _occ_ax2_frame(axis)[1]

# Second in-plane direction of an OCC circle record — stored explicitly
# (`gp_Ax2` frames need not satisfy Y = n×X after transforms).
@inline _occ_circle_y(g) = g.Y

# `ElCLib::CircleValue`: A1·X + A2·Y + PLoc with A1 = R·cos, A2 = R·sin;
# `CircleD1`/`CircleD2` use the same `SetLinearForm` association — the first
# product of each component fused, the second rounded.
function _occ_circle_point(g,t::Float64)
    Y=_occ_circle_y(g)
    s,c=_gm_sincos(t)
    A1,A2=g.r*c,g.r*s
    C,X=g.center,g.X
    return (fma(A1,X[1],A2*Y[1])+C[1],
            fma(A1,X[2],A2*Y[2])+C[2],
            fma(A1,X[3],A2*Y[3])+C[3])
end

function _occ_circle_derivative(g,t::Float64)
    Y=_occ_circle_y(g)
    s,c=_gm_sincos(t)
    Xc,Yc=g.r*c,g.r*s
    X=g.X
    return (fma(-Yc,X[1],Xc*Y[1]),
            fma(-Yc,X[2],Xc*Y[2]),
            fma(-Yc,X[3],Xc*Y[3]))
end

function _occ_circle_second_derivative(g,t::Float64)
    Y=_occ_circle_y(g)
    s,c=_gm_sincos(t)
    Xc,Yc=g.r*c,g.r*s
    X=g.X
    return (fma(-Xc,X[1],(-Yc)*Y[1]),
            fma(-Xc,X[2],(-Yc)*Y[2]),
            fma(-Xc,X[3],(-Yc)*Y[3]))
end

# `ElCLib::LineValue`/`LineD1` — `U·dir + loc` on the stored `gp_Lin`
# parametrization (not the endpoint chord).
@inline _occ_line_point(g,t::Float64) =
    (fma(t,g.dir[1],g.origin[1]), fma(t,g.dir[2],g.origin[2]),
     fma(t,g.dir[3],g.origin[3]))

# `ElCLib2d` pcurve eval — the `Geom2d_Line`/`Geom2d_Circle` `D0`/`D1`
# bodies used by `Adaptor2d_Curve2d` on a stored face pcurve.
function _occ_pcurve_d0(pc,t::Float64)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return _occ_lin2d_eval(o,d,t)
    else
        o,x,r=pc.circ2d
        return _occ_circ2d_eval(o,x,r,t)
    end
end
function _occ_pcurve_d1(pc,t::Float64)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return _occ_lin2d_eval(o,d,t), d
    else
        o,x,r=pc.circ2d
        # `Geom2d_Circle::D1` — SetLinearForm(−Yc, X, Xc, Y) with the derived
        # YDirection (−x.y, x.x).
        s,c=_gm_sincos(t)
        A1,A2=r*c,r*s
        uv=_occ_circ2d_eval(o,x,r,t)
        duv=(fma(-A2,x[1],-(A1*x[2])), fma(-A2,x[2],A1*x[1]))
        return uv,duv
    end
end

# `ElCLib2d::CircleD2`/`LineD2` on the stored pcurve kinds — returns
# (point, d1, d2). The circle's D2 negates the point's linear form before
# the location add, exactly as `CircleD2` builds `V2` then `P`.
function _occ_pcurve_d2(pc,t::Float64)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return _occ_lin2d_eval(o,d,t), d, (0.0,0.0)
    end
    o,x,r=pc.circ2d
    s,c=_gm_sincos(t)
    Xc,Yc=r*c,r*s
    # YDirection is the derived (−x.y, x.x); `V2 = −(Xc·X + Yc·Y)`.
    p=_occ_circ2d_eval(o,x,r,t)
    d1=(fma(-Yc,x[1],-(Xc*x[2])), fma(-Yc,x[2],Xc*x[1]))
    d2=(-fma(Xc,x[1],-(Yc*x[2])), -fma(Xc,x[2],Yc*x[1]))
    return p,d1,d2
end

# `SearchForExtremum` — the Newton step `t -= (D1·dir)/(D2·dir)`, at most 10
# iterations; converges when |D2·dir| < 1e-10 or |Δt| < Precision::PConfusion
# (1e-9), fails after a third excursion outside [f, l] or a repeated clamp.
# On convergence `res` is the D0 evaluated at the parameter used in the last
# `D2` call (the pre-step parameter on the |Δt| exit).
function _occ_pcurve_extremum(pc,f,l,dir,par)
    nbOut=0; res=(0.0,0.0)
    for _ in 1:10
        prev=par
        res,d1,d2=_occ_pcurve_d2(pc,par)
        det=fma(d2[1],dir[1],d2[2]*dir[2])
        abs(det)<1e-10 && return res
        par=prev-(fma(d1[1],dir[1],d1[2]*dir[2]))/det
        abs(par-prev)<1e-9 && return res
        if par<f
            nbOut+=1
            (nbOut>2 || prev==f) && return nothing
            par=f
        end
        if par>l
            nbOut+=1
            (nbOut>2 || prev==l) && return nothing
            par=l
        end
    end
    return res
end

# `ShapeAnalysis_Curve::FillBndBox(c2d, f, l, 20, Exact=true)` — elementary
# pcurves are a single C2 interval, so nbSamples = 19: twenty sample points
# `f + i·(l−f)/19` are boxed, then a Newton extremum along x̂ and ŷ is seeded
# at each interval midpoint. Returns ((umin,vmin),(umax,vmax)).
function _occ_pcurve_fillbndbox(pc,f::Float64,l::Float64)
    step=(l-f)/19
    umin=vmin=Inf; umax=vmax=-Inf
    for i in 0:19
        p=_occ_pcurve_d0(pc,fma(Float64(i),step,f))
        umin=min(umin,p[1]); vmin=min(vmin,p[2])
        umax=max(umax,p[1]); vmax=max(vmax,p[2])
    end
    for i in 0:18
        a1=fma(Float64(i),step,f); a2=fma(Float64(i+1),step,f)
        mid=(a1+a2)*0.5
        for dir in ((1.0,0.0),(0.0,1.0))
            e=_occ_pcurve_extremum(pc,a1,a2,dir,mid)
            e===nothing && continue
            umin=min(umin,e[1]); vmin=min(vmin,e[2])
            umax=max(umax,e[1]); vmax=max(vmax,e[2])
        end
    end
    return (umin,vmin),(umax,vmax)
end

# `ShapeAnalysis::GetFaceUVBounds` — the union of `FillBndBox` boxes over
# every wire pcurve. Each dict entry carries the forward pcurve and, for a
# seam occurrence on a closed surface, the reverse one (`PCurve2`); both are
# boxed since `TopExp_Explorer` yields the seam twice. `f,l` are the edge
# parameter ranges (`BRep_Tool::CurveOnSurface` out-parameters).
function _occ_face_uv_bounds(m::GeoModel,g)
    umin=vmin=Inf; umax=vmax=-Inf; found=false
    for (ct,e) in g.pcurves
        rec=m.curve_geometry[ct]
        for pc in (e.fwd,e.rev)
            pc===nothing && continue
            found=true
            lo,hi=_occ_pcurve_fillbndbox(pc,rec.t0,rec.t1)
            umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
            umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
        end
    end
    found || throw(ArgumentError(
        "_occ_face_uv_bounds: face carries no pcurves"))
    return [umin,vmin],[umax,vmax]
end

# `BRep_Tool::CurveOnSurface(E, F)` as stored on `face_g`: `which=1` reads the
# FORWARD-orientation pcurve (`PCurve`), `which=0` the REVERSED one
# (`PCurve2` on a closed-surface record, otherwise the same single pcurve).
# `nothing` when the edge carries no pcurve on this face.
function _occ_edge_pcurve_on_face(face_g, edge_tag::Int, which::Int)
    e=get(face_g.pcurves,edge_tag,nothing)
    e===nothing && return nothing
    return which==1 ? e.fwd : (e.rev===nothing ? e.fwd : e.rev)
end

# `OCCEdge::reparamOnFace`: the edge's face pcurve evaluated at the edge
# parameter directly — no 3-D evaluation, no robustness recheck
# (`reparamOnFaceRobust` defaults to 0). Edges without a pcurve on the face
# take the `closestPoint(sp, {0,0})` fallback.
function _occ_edge_reparam_on_face(m::GeoModel, face_g, edge_tag::Int,
                                   epar::Float64, which::Int,
                                   lc::Float64, caller::AbstractString)
    pc=_occ_edge_pcurve_on_face(face_g,edge_tag,which)
    pc!==nothing && return _occ_pcurve_d0(pc,epar)
    pt=_model_curve_point(m,edge_tag,epar,caller,0)
    u,v,_=_occ_surface_closest(face_g,pt,lc)
    return (u,v)
end

# `OCCVertex::reparamOnFace`: the first `l_edges` edge that also bounds the
# face supplies its pcurve range (`s0` at the begin vertex, `s1` at the end);
# the vertex then delegates to the edge's `reparamOnFace`. Incident-edge
# order is curve-creation order, matching the `TopExp_Explorer` order in
# which OCC edges adopt their vertices. Falls back to
# `OCCFace::parFromPoint` when no shared OCC edge exists.
function _occ_vertex_reparam_on_face(m::GeoModel, surface_tag::Int,
                                     face_g, vtag::Int, which::Int,
                                     lc::Float64)
    incident=Int[]
    for (ct,(a,b)) in m.curves
        (a==vtag || b==vtag) && push!(incident,ct)
    end
    sort!(incident)
    for ct in incident
        onface=false
        for loop in m.surfaces[surface_tag], signed in m.loops[loop]
            abs(signed)==ct && (onface=true; break)
        end
        onface || continue
        rec=get(m.curve_geometry,ct,nothing)
        pc=rec!==nothing && hasproperty(rec,:occ) ?
            _occ_edge_pcurve_on_face(face_g,ct,which) : nothing
        pc===nothing && break   # non-OCC shared edge → parFromPoint below
        a,_=m.curves[ct]
        return _occ_pcurve_d0(pc, a==vtag ? rec.t0 : rec.t1)
    end
    return _occ_surface_parameter_on_face(face_g,m.points[vtag],lc)
end

# `BRepAdaptor_Curve`'s pcurve path for a degenerate edge: the edge's stored
# pcurve on `g.face` (`_trimmed`/`_curve2d` bound by `_occ_bind_degenerate!`).
# `OCCEdge::point` returns `_trimmed->point(u,v)`;
# `Adaptor3d_CurveOnSurface::D1` applies `D1 = du·S_u + dv·S_v` (the
# `EvalKPart` `GeomAbs_OtherCurve` branch — pole pcurves are VIsos at the
# singular latitude, which EvalKPart refuses to re-recognize).
function _occ_degenerate_pcurve(g,tag::Int)
    hasproperty(g,:face) || throw(ArgumentError(
        "_occ_degenerate_pcurve: degenerate edge has no trimmed face"))
    pc=get(g.face.pcurves,tag,nothing)
    pc!==nothing && pc.fwd!==nothing ||
        throw(ArgumentError(
            "_occ_degenerate_pcurve: degenerate edge has no pcurve on its face"))
    return pc.fwd
end
function _occ_degenerate_point(g,tag::Int,t::Float64)
    u,v=_occ_pcurve_d0(_occ_degenerate_pcurve(g,tag),t)
    return _occ_surface_point(g.face,u,v)
end
function _occ_degenerate_d1(g,tag::Int,t::Float64)
    (u,v),duv=_occ_pcurve_d1(_occ_degenerate_pcurve(g,tag),t)
    _,su,sv=_occ_surface_d1(g.face,u,v)
    return ntuple(i->fma(duv[1],su[i],duv[2]*sv[i]),3)
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
        base=_gm_atan2(Y[axis],g.X[axis])
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
        return _points_close(
                   _add3(g.origin,_occ_mul(g.dir,g.t0)),pa,tol) &&
               _points_close(
                   _add3(g.origin,_occ_mul(g.dir,g.t1)),pb,tol)
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
                         X::NTuple{3,Float64}, Y::NTuple{3,Float64},
                         r::Float64, t0::Float64, t1::Float64)
    t=_alloc_tag!(m,1,0,"_add_occ_circle!")
    m.curves[t]=(p1,p2)
    m.curve_types[t]=:circle
    m.curve_geometry[t]=(occ=:circle,center=center,n=n,X=X,Y=Y,r=r,
                         t0=t0,t1=t1)
    return t
end

# A seam line: `origin`/`dir` are the stored `gp_Lin` parametrization (with
# the R(2π) transform residues OCC bakes in); the range is [0,length].
function _add_occ_line!(m::GeoModel, p1::Int, p2::Int,
                        origin::NTuple{3,Float64}, dir::NTuple{3,Float64},
                        t1::Float64)
    t=_alloc_tag!(m,1,0,"_add_occ_line!")
    m.curves[t]=(p1,p2)
    m.curve_geometry[t]=(occ=:line,origin=origin,dir=dir,t0=0.0,t1=t1)
    return t
end

# A degenerate OCC edge: a zero-length edge collapsed on a pole or apex
# vertex. `face` binds the surface record OCC's `setTrimmed` attaches (the
# first face adopting the edge); the edge's pcurve is looked up in
# `face.pcurves` at evaluation time.
function _add_occ_degenerate!(m::GeoModel, p::Int)
    t=_alloc_tag!(m,1,0,"_add_occ_degenerate!")
    m.curves[t]=(p,p)
    m.curve_types[t]=:degenerate
    m.curve_geometry[t]=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI)
    return t
end

# `OCCFace` ctor `setTrimmed`: the first face adopting a non-3D edge binds it
# (`_trimmed` + `_curve2d`). Called by the materializers once the surface
# record exists.
function _occ_bind_degenerate!(m::GeoModel, edge::Int, face)
    m.curve_geometry[edge]=merge(m.curve_geometry[edge],(face=face,))
    return edge
end

# A typed OCC face: `loops` are existing curve-loop tags, `kind` one of
# :cylinder/:sphere/:cone/:torus/:plane, `geometry` the analytic record.
# A `:plane` record additionally caches `uvb` — the `GetFaceUVBounds` result
# `OCCFace`'s constructor stores as `_umin.._vmax` (the face bounds are fixed
# once the pcurves exist; transforms leave them invariant).
function _add_occ_surface!(m::GeoModel, kind::Symbol, loops::Vector{Int},
                           geometry)
    if geometry.occ===:plane
        lo,hi=_occ_face_uv_bounds(m,geometry)
        geometry=merge(geometry,(uvb=(lo,hi),))
    end
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
        delete!(m.curve_params,curve)
    end
    for point in points
        delete!(m.points,point)
        delete!(m.point_size,point)
    end
    return nothing
end

# A full-rotation `gp_Trsf` residue: OCCT's `EndEdge`/`TopEndVertex` chain
# transforms the unrotated construction by R(2π) about the build axis, so
# seam origins and rim vertices carry `sin(2π)` residuals.
_occ_end_trsf(base, n) = _occ_trsf_rot(base, n, _OCC_TWO_PI)

# `gp_Pnt::Translated` — plain componentwise add of a scaled direction.
@inline _occ_translated(p,d,s) =
    (p[1]+d[1]*s, p[2]+d[2]*s, p[3]+d[3]*s)

# OCC's cylinder layout, built in `BRepPrim_OneAxis` order: top rim Point
# (TopEndVertex), bottom rim Point (BottomEndVertex); top Circle (ETOP),
# seam Line (ESTART/EEND shared), bottom Circle (EBOTTOM); the Cylinder
# lateral face ([-top,-seam,bottom,seam]), Plane caps, shell
# [lateral,+top,-bottom]. Rim frames come from the `gp_Ax2(P,N,X)` 3-arg
# ctor; the seam line carries the `Geom_Line::Transform(R(2π))` origin/dir.
function _materialize_cylinder!(m::GeoModel, base::NTuple{3,Float64},
                                axis::NTuple{3,Float64}, r::Float64,
                                h::Float64)
    n=_occ_dir(axis)
    X,Y=_occ_ax2_frame(n)
    M,loc=_occ_end_trsf(base,n)
    # Rim/cap frames: `gp_Ax2(Loc + n·mp.Y, n, X)` 3-arg ctor; the stored
    # circle axis is the adaptor view X×Y.
    Xc,Yc=_occ_ax2_setx(n,X)
    nc=_occ_crossed(Xc,Yc)
    # Meridian profile `Lin2d((r,0),(0,1))`: mp(v) = (r, v).
    w_top=_occ_vertex(base,n,X,(r,h),M,loc)
    w_bot=_occ_vertex(base,n,X,(r,0.0),M,loc)
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        push!(points,_add_occ_point!(m,w_top))
        push!(points,_add_occ_point!(m,w_bot))
        for p in points
            delete!(m.point_size,p)
        end
        p_top,p_bot=points
        top_center=_occ_translated(base,n,h)
        push!(curves,_add_occ_circle!(m,p_top,p_top,top_center,nc,Xc,Yc,r,
                                      0.0,_OCC_TWO_PI))
        # `Geom_Line(Ax1(base + X·r, n)).Transform(R(2π))`, range [0,h].
        seam_origin=_occ_xform_pnt(M,loc,_occ_translated(base,X,r))
        seam_dir=_occ_xform_dir(M,n)
        push!(curves,_add_occ_line!(m,p_bot,p_top,seam_origin,seam_dir,h))
        push!(curves,_add_occ_circle!(m,p_bot,p_bot,base,nc,Xc,Yc,r,
                                      0.0,_OCC_TWO_PI))
        c_top,c_seam,c_bot=curves
        push!(loops,add_curve_loop!(m,[-c_top,-c_seam,c_bot,c_seam]))
        push!(surfaces,_add_occ_surface!(m,:cylinder,[last(loops)],
              (occ=:cylinder,center=base,axis=n,X=X,Y=Y,radius=r,height=h,
               pcurves=Dict{Int,NamedTuple}(
                   c_top=>(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing),
                   c_bot=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing),
                   c_seam=>(fwd=(lin2d=((_OCC_TWO_PI,0.0),(0.0,1.0)),),
                            rev=(lin2d=((0.0,0.0),(0.0,1.0)),))))))
        # Cap planes: `gp_Pln(gp_Ax3(Loc + n·mp.Y, n, X))` top FORWARD and
        # bottom REVERSED; each rim edge's pcurve is `Circ2d((0,0),x̂,r)`.
        push!(loops,add_curve_loop!(m,[c_top]))
        push!(surfaces,_add_occ_surface!(m,:plane,[last(loops)],
              (occ=:plane,center=top_center,axis=n,X=Xc,Y=Yc,reversed=false,
               pcurves=Dict{Int,NamedTuple}(
                   c_top=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),r),),
                           rev=nothing)))))
        push!(loops,add_curve_loop!(m,[c_bot]))
        push!(surfaces,_add_occ_surface!(m,:plane,[last(loops)],
              (occ=:plane,center=base,axis=n,X=Xc,Y=Yc,reversed=true,
               pcurves=Dict{Int,NamedTuple}(
                   c_bot=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),r),),
                           rev=nothing)))))
        s_lat,s_top,s_bot=surfaces
        shell=add_surface_loop!(m,[s_lat,s_top,-s_bot])
    catch
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,shell)
        rethrow()
    end
    return shell
end

# OCC's sphere layout: north/south pole Points (`TopEndVertex`/`BottomEndVertex`
# — the `Circ2d` profile eval at ±π/2 rotated by R(2π)); a degenerate edge on
# each pole and a meridian Circle — `gp_Ax2(Loc, -Y, X)` transformed by
# R(2π), trimmed to [3π/2,5π/2] running south→north; one Sphere face
# ([-degN,-meridian,degS,meridian]); shell [face]. The surface frame is the
# `gp_Ax2(Loc, ẑ, x̂)` 3-arg ctor output.
function _materialize_sphere!(m::GeoModel, center::NTuple{3,Float64},
                              r::Float64)
    n=(0.0,0.0,1.0)
    X,Y=_occ_ax2_setx(n,(1.0,0.0,0.0))
    M,loc=_occ_end_trsf(center,n)
    # Profile `Circ2d((0,0),x̂,r)`: mp(v) = (r·cos v, r·sin v); pole vertices
    # evaluate at the raw profile extrema ±π/2.
    mp(v)=_occ_circ2d_eval((0.0,0.0),(1.0,0.0),r,v)
    w_n=_occ_vertex(center,n,X,mp(π/2),M,loc)
    w_s=_occ_vertex(center,n,X,mp(-π/2),M,loc)
    # Meridian ctor `gp_Ax2(Loc, -Y, X)`, then `gp_Circ::Transform(R(2π))`.
    Dm=(-Y[1],-Y[2],-Y[3])
    Xm,Ym=_occ_ax2_setx(Dm,X)
    mer=_occ_xform_ax2(M,loc,(loc=center,axis=Dm,X=Xm,Y=Ym))
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        push!(points,_add_occ_point!(m,w_n))
        push!(points,_add_occ_point!(m,w_s))
        for p in points
            delete!(m.point_size,p)
        end
        p_n,p_s=points
        push!(curves,_add_occ_degenerate!(m,p_n))
        push!(curves,_add_occ_circle!(m,p_s,p_n,mer.loc,mer.axis,mer.X,mer.Y,
                                      r,1.5π,2.5π))
        push!(curves,_add_occ_degenerate!(m,p_s))
        c_n,c_mer,c_s=curves
        push!(loops,add_curve_loop!(m,[-c_n,-c_mer,c_s,c_mer]))
        push!(surfaces,_add_occ_surface!(m,:sphere,[last(loops)],
              (occ=:sphere,center=center,radius=r,
               axis=n,X=X,Y=Y,
               pcurves=Dict{Int,NamedTuple}(
                   c_n=>(fwd=(lin2d=((0.0,π/2),(1.0,0.0)),),rev=nothing),
                   c_s=>(fwd=(lin2d=((0.0,-π/2),(1.0,0.0)),),rev=nothing),
                   c_mer=>(fwd=(lin2d=((_OCC_TWO_PI,-_OCC_TWO_PI),
                                       (0.0,1.0)),),
                           rev=(lin2d=((0.0,-_OCC_TWO_PI),(0.0,1.0)),))))))
        shell=add_surface_loop!(m,[last(surfaces)])
        face=m.surface_geometry[last(surfaces)]
        _occ_bind_degenerate!(m,c_n,face)
        _occ_bind_degenerate!(m,c_s,face)
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
    n=_occ_dir(axis)
    X,Y=_occ_ax2_frame(n)
    M,loc=_occ_end_trsf(base,n)
    dr=r2-r1
    sa=_gm_atan(dr/h)
    len=sqrt(fma(h,h,dr*dr))
    # Meridian profile `Lin2d((r1,0),Dir2d(sin sa,cos sa))`; `gp_Dir2d`
    # renormalizes, so mp.x differs from r2 by ~1ulp — the rim circle stores
    # the evaluated mp.x, which is why OCC's cone curvature can read
    # 1/0.9999999999999998.
    dir2d=_occ_dir2(_gm_sincos(sa))
    mp(v)=_occ_lin2d_eval((r1,0.0),dir2d,v)
    mpT=mp(len); mpB=mp(0.0)
    Xc,Yc=_occ_ax2_setx(n,X)
    nc=_occ_crossed(Xc,Yc)
    # `A.Rotate(Ax1(Loc,Y),sa)` — `gp_Dir::Rotate`, matvec only, no
    # renormalization.
    MY,_=_occ_trsf_rot(base,Y,sa)
    d0=_occ_rot_dir(MY,n)
    w_top=_occ_vertex(base,n,X,mpT,M,loc)
    w_bot=_occ_vertex(base,n,X,mpB,M,loc)
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        push!(points,_add_occ_point!(m,w_top))
        push!(points,_add_occ_point!(m,w_bot))
        for p in points
            delete!(m.point_size,p)
        end
        p_top,p_bot=points
        top_center=_occ_translated(base,n,mpT[2])
        push!(curves,r2>0 ?
              _add_occ_circle!(m,p_top,p_top,top_center,nc,Xc,Yc,mpT[1],
                               0.0,_OCC_TWO_PI) :
              _add_occ_degenerate!(m,p_top))
        # `Geom_Line(Ax1(base + X·r1, d0)).Transform(R(2π))`.
        seam_origin=_occ_xform_pnt(M,loc,_occ_translated(base,X,r1))
        seam_dir=_occ_xform_dir(M,d0)
        push!(curves,_add_occ_line!(m,p_bot,p_top,seam_origin,seam_dir,len))
        push!(curves,r1>0 ?
              _add_occ_circle!(m,p_bot,p_bot,base,nc,Xc,Yc,mpB[1],
                               0.0,_OCC_TWO_PI) :
              _add_occ_degenerate!(m,p_bot))
        c_top,c_seam,c_bot=curves
        push!(loops,add_curve_loop!(m,[-c_top,-c_seam,c_bot,c_seam]))
        push!(surfaces,_add_occ_surface!(m,:cone,[last(loops)],
              (occ=:cone,center=base,axis=n,X=Xc,Y=Yc,r1=r1,r2=r2,height=h,
               pcurves=Dict{Int,NamedTuple}(
                   c_top=>(fwd=(lin2d=((0.0,len),(1.0,0.0)),),rev=nothing),
                   c_bot=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing),
                   c_seam=>(fwd=(lin2d=((_OCC_TWO_PI,-0.0),(0.0,1.0)),),
                            rev=(lin2d=((0.0,-0.0),(0.0,1.0)),))))))
        s_lat=last(surfaces)
        lat_face=m.surface_geometry[s_lat]
        r2==0 && _occ_bind_degenerate!(m,c_top,lat_face)
        r1==0 && _occ_bind_degenerate!(m,c_bot,lat_face)
        shell_signs=Int[s_lat]
        if r2>0
            push!(loops,add_curve_loop!(m,[c_top]))
            push!(surfaces,_add_occ_surface!(m,:plane,[last(loops)],
                  (occ=:plane,center=top_center,axis=n,X=Xc,Y=Yc,
                   reversed=false,
                   pcurves=Dict{Int,NamedTuple}(
                       c_top=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),mpT[1]),),
                               rev=nothing)))))
            push!(shell_signs,last(surfaces))
        end
        if r1>0
            push!(loops,add_curve_loop!(m,[c_bot]))
            push!(surfaces,_add_occ_surface!(m,:plane,[last(loops)],
                  (occ=:plane,center=base,axis=n,X=Xc,Y=Yc,
                   reversed=true,
                   pcurves=Dict{Int,NamedTuple}(
                       c_bot=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),mpB[1]),),
                               rev=nothing)))))
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
    n=(0.0,0.0,1.0)
    X,Y=_occ_ax2_setx(n,(1.0,0.0,0.0))
    M,loc=_occ_end_trsf(center,n)
    # Profile `Circ2d((r1,0),x̂,r2)`: mp(v) = (r1 + r2·cos v, r2·sin v). The
    # rim vertex is `TopEndVertex`/`BottomStartVertex` — `mp(VMax)` carries
    # the `r2·sin(2π)` residue.
    mp(v)=_occ_circ2d_eval((r1,0.0),(1.0,0.0),r2,v)
    mpE=mp(_OCC_TWO_PI)
    w_rim=_occ_vertex(center,n,X,mpE,M,loc)
    # Meridian ctor `gp_Ax2(Loc + X·r1, -Y, X)`; the stored circle is
    # `gp_Circ::Transform`'d by R(2π) (full) or R(angle)/R(0) (partial).
    Dm=(-Y[1],-Y[2],-Y[3])
    Xm,Ym=_occ_ax2_setx(Dm,X)
    locm=_occ_translated(center,X,r1)
    # Equator `gp_Circ(Ax2(Loc + n·mpE.y, n, X), mpE.x)` — the
    # `r2·sin(2π)`-offset center.
    Xe,Ye=_occ_ax2_setx(n,X)
    ne=_occ_crossed(Xe,Ye)
    ce=_occ_translated(center,n,mpE[2])
    full=angle>=_OCC_TWO_PI
    points=Int[]; curves=Int[]; loops=Int[]; surfaces=Int[]; shell=0
    try
        if full
            push!(points,_add_occ_point!(m,w_rim))
        else
            # VTOPEND (`mp(VMax)` rotated by the sweep angle), then
            # VTOPSTART — the `MeridianClosed` deduction makes it the
            # unrotated `mp(VMax)` vertex, carrying the `r2·sin(2π)`
            # residue on n.
            Ma,loca=_occ_trsf_rot(center,n,angle)
            push!(points,_add_occ_point!(m,
                _occ_vertex(center,n,X,mpE,Ma,loca)))
            push!(points,_add_occ_point!(m,
                ntuple(i->(center[i]+n[i]*mpE[2])+X[i]*mpE[1],3)))
        end
        for p in points
            delete!(m.point_size,p)
        end
        if full
            p_rim=points[1]
            mer=_occ_xform_ax2(M,loc,(loc=locm,axis=Dm,X=Xm,Y=Ym))
            push!(curves,_add_occ_circle!(m,p_rim,p_rim,ce,ne,Xe,Ye,mpE[1],
                  0.0,_OCC_TWO_PI))
            push!(curves,_add_occ_circle!(m,p_rim,p_rim,mer.loc,mer.axis,
                  mer.X,mer.Y,r2,0.0,_OCC_TWO_PI))
            c_eq,c_mer=curves
            push!(loops,add_curve_loop!(m,[-c_eq,c_mer,c_eq,-c_mer]))
            push!(surfaces,_add_occ_surface!(m,:torus,[last(loops)],
                  (occ=:torus,center=center,axis=n,X=X,Y=Y,r1=r1,r2=r2,
                   angle=angle,
                   pcurves=Dict{Int,NamedTuple}(
                       c_eq=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),
                              rev=(lin2d=((0.0,_OCC_TWO_PI),(1.0,0.0)),)),
                       c_mer=>(fwd=(lin2d=((_OCC_TWO_PI,-0.0),(0.0,1.0)),),
                               rev=(lin2d=((0.0,-0.0),(0.0,1.0)),))))))
            shell=add_surface_loop!(m,[last(surfaces)])
        else
            p_end,p_start=points
            Ma,loca=_occ_trsf_rot(center,n,angle)
            mer_e=_occ_xform_ax2(Ma,loca,(loc=locm,axis=Dm,X=Xm,Y=Ym))
            # The start meridian goes through `gp_Ax2::Transform` with the
            # angle-0 rotation Trsf (M=I, loc=0): directions normalize to
            # themselves, only the axis is recomputed as X×Y.
            mer_s=_occ_xform_ax2(
                ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)),(0.0,0.0,0.0),
                (loc=locm,axis=Dm,X=Xm,Y=Ym))
            push!(curves,_add_occ_circle!(m,p_start,p_end,ce,ne,Xe,Ye,mpE[1],
                  0.0,angle))
            push!(curves,_add_occ_circle!(m,p_end,p_end,mer_e.loc,mer_e.axis,
                  mer_e.X,mer_e.Y,r2,0.0,_OCC_TWO_PI))
            push!(curves,_add_occ_circle!(m,p_start,p_start,mer_s.loc,
                  mer_s.axis,mer_s.X,mer_s.Y,r2,0.0,_OCC_TWO_PI))
            c_arc,c_end,c_start=curves
            push!(loops,add_curve_loop!(m,[-c_arc,-c_start,c_arc,c_end]))
            push!(surfaces,_add_occ_surface!(m,:torus,[last(loops)],
                  (occ=:torus,center=center,axis=n,X=X,Y=Y,r1=r1,r2=r2,
                   angle=angle,
                   pcurves=Dict{Int,NamedTuple}(
                       c_arc=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),
                               rev=(lin2d=((0.0,_OCC_TWO_PI),(1.0,0.0)),)),
                       c_start=>(fwd=(lin2d=((0.0,-0.0),(0.0,1.0)),),
                                  rev=nothing),
                       c_end=>(fwd=(lin2d=((angle,-0.0),(0.0,1.0)),),
                                rev=nothing)))))
            s_torus=last(surfaces)
            # Caps: `gp_Pln(Ax2(Loc,-Y,X))` forward on the start side, the
            # `gp_Ax2::Rotate`d plane + `ReverseFace` on the end side; each
            # meridian's cap pcurve is the profile `Circ2d((r1,0),x̂,r2)`.
            push!(loops,add_curve_loop!(m,[c_start]))
            push!(surfaces,_add_occ_surface!(m,:plane,[last(loops)],
                  (occ=:plane,center=center,axis=Dm,X=Xm,Y=Ym,reversed=false,
                   pcurves=Dict{Int,NamedTuple}(
                       c_start=>(fwd=(circ2d=((r1,0.0),(1.0,0.0),r2),),
                                 rev=nothing)))))
            s_start=last(surfaces)
            cap_e=_occ_rotate_ax2(Ma,loca,(loc=center,axis=Dm,X=Xm,Y=Ym))
            push!(loops,add_curve_loop!(m,[c_end]))
            push!(surfaces,_add_occ_surface!(m,:plane,[last(loops)],
                  (occ=:plane,center=cap_e.loc,axis=cap_e.axis,X=cap_e.X,
                   Y=cap_e.Y,reversed=true,
                   pcurves=Dict{Int,NamedTuple}(
                       c_end=>(fwd=(circ2d=((r1,0.0),(1.0,0.0),r2),),
                               rev=nothing)))))
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

# Stored surface YDirection — `gp_Ax3` keeps all three directions verbatim
# through transforms (a reflection makes the frame left-handed).
@inline _occ_frame_y(g) = g.Y

# `ElSLib::*Value`: A1·X + A2·Y + A3·Z + PLoc with the coefficient products
# pre-rounded — `(A1·X + A2·Y) + A3·Z` chains the fused multiply-adds in
# written order before the location add. `TorusValue` additionally zeroes
# coefficients below `10·(r_minor+r_major)·eps` (the OCC620 clamp).
function _occ_surface_point(g, u::Float64, v::Float64)
    Y=_occ_frame_y(g)
    C,X,Z=g.center,g.X,g.axis
    if g.occ===:plane
        # `ElSLib::PlaneValue` — U·X + V·Y + Loc.
        return (fma(u,X[1],v*Y[1])+C[1],
                fma(u,X[2],v*Y[2])+C[2],
                fma(u,X[3],v*Y[3])+C[3])
    end
    s,c=_gm_sincos(u)
    if g.occ===:cylinder
        A1,A2,A3=g.radius*c,g.radius*s,v
    elseif g.occ===:sphere
        sv,cv=_gm_sincos(v)
        R,A3=g.radius*cv,g.radius*sv
        A1,A2=R*c,R*s
    elseif g.occ===:cone
        # OCCT stores the half-angle atan((R2−R1)/H) (BRepPrim_Cone::
        # SetParameters) and ConeValue calls sin/cos on it per eval; the
        # slant-ratio forms differ by ~1ulp from the libm trig values.
        sa=_gm_atan((g.r2-g.r1)/g.height)
        ssa,csa=_gm_sincos(sa)
        R,A3=fma(v,ssa,g.r1),v*csa
        A1,A2=R*c,R*s
    elseif g.occ===:torus
        sv,cv=_gm_sincos(v)
        R,A3=fma(g.r2,cv,g.r1),g.r2*sv
        A1,A2=R*c,R*s
        clamp_eps=10.0*(g.r2+g.r1)*eps()
        abs(A1)<=clamp_eps && (A1=0.0)
        abs(A2)<=clamp_eps && (A2=0.0)
        abs(A3)<=clamp_eps && (A3=0.0)
    else
        throw(ArgumentError(
            "_occ_surface_point: unsupported OCC surface kind $(g.occ)"))
    end
    return (fma(A3,Z[1],fma(A1,X[1],A2*Y[1]))+C[1],
            fma(A3,Z[2],fma(A1,X[2],A2*Y[2]))+C[2],
            fma(A3,Z[3],fma(A1,X[3],A2*Y[3]))+C[3])
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
    elseif g.occ===:plane
        return g.uvb
    end
    throw(ArgumentError(
        "_occ_surface_bounds: unsupported OCC surface kind $(g.occ)"))
end

# `ElSLib::*D1` ports — `mySurface->D1` behind `BRepAdaptor_Surface::D1` /
# `OCCFace::firstDer`, with each function's own coefficient contractions and
# the OCC620 torus clamp. Returns (point, ∂P/∂u, ∂P/∂v).
function _occ_surface_d1(g,u::Float64,v::Float64)
    Y=_occ_frame_y(g)
    X,Z=g.X,g.axis
    p=_occ_surface_point(g,u,v)
    if g.occ===:plane
        # `ElSLib::PlaneD1` — Vu = XDirection, Vv = YDirection.
        return p,X,Y
    end
    s,c=_gm_sincos(u)
    # The `du` coefficients and `dv` are bound once as an if-expression so the
    # ntuple closures capture single-assignment locals (no Core.Box).
    A1,A2,dv=if g.occ===:cylinder
        (g.radius*c,g.radius*s,Z)
    elseif g.occ===:sphere
        let (sv,cv)=_gm_sincos(v),R1=g.radius*cv,R2=g.radius*sv,
            A3=R2*c,A4=R2*s
            (R1*c,R1*s,
             ntuple(i->fma(R1,Z[i],fma(-A3,X[i],-(A4*Y[i]))),3))
        end
    elseif g.occ===:cone
        let sa=_gm_atan((g.r2-g.r1)/g.height),(SinA,CosA)=_gm_sincos(sa),
            R=fma(v,SinA,g.r1),R1=SinA*c,R2=SinA*s
            (R*c,R*s,
             ntuple(i->fma(CosA,Z[i],fma(R1,X[i],R2*Y[i])),3))
        end
    elseif g.occ===:torus
        let (sv,cv)=_gm_sincos(v),R1=g.r2*cv,R2=g.r2*sv,
            R=fma(g.r2,cv,g.r1),clamp_eps=10.0*(g.r2+g.r1)*eps(),
            a1=abs(R*c)<=clamp_eps ? 0.0 : R*c,
            a2=abs(R*s)<=clamp_eps ? 0.0 : R*s,
            a3=abs(R2*c)<=clamp_eps ? 0.0 : R2*c,
            a4=abs(R2*s)<=clamp_eps ? 0.0 : R2*s
            (a1,a2,ntuple(i->fma(R1,Z[i],fma(-a3,X[i],-(a4*Y[i]))),3))
        end
    else
        throw(ArgumentError(
            "_occ_surface_d1: unsupported OCC surface kind $(g.occ)"))
    end
    du=ntuple(i->fma(-A2,X[i],A1*Y[i]),3)
    return p,du,dv
end

# `ElSLib::*D2` ports — returns (point, du, dv, duu, dvv, duv). `Som1` is the
# shared `A1·X + A2·Y` term; `Vuu` is its negation.
function _occ_surface_d2(g,u::Float64,v::Float64)
    Y=_occ_frame_y(g)
    X,Z=g.X,g.axis
    p,du,dv=_occ_surface_d1(g,u,v)
    if g.occ===:plane
        # `ElSLib::PlaneDN` — all second-order partials vanish.
        z=(0.0,0.0,0.0)
        return p,du,dv,z,z,z
    end
    s,c=_gm_sincos(u)
    if g.occ===:cylinder
        let A1=g.radius*c,A2=g.radius*s,
            som1=ntuple(i->fma(A1,X[i],A2*Y[i]),3)
            return p,du,dv,.-som1,(0.0,0.0,0.0),(0.0,0.0,0.0)
        end
    elseif g.occ===:sphere
        let (sv,cv)=_gm_sincos(v),R1=g.radius*cv,R2=g.radius*sv,
            A1=R1*c,A2=R1*s,A3=R2*c,A4=R2*s,
            som1=ntuple(i->fma(A1,X[i],A2*Y[i]),3)
            return p,du,dv,.-som1,
                   ntuple(i->-som1[i]-R2*Z[i],3),
                   ntuple(i->fma(A4,X[i],-(A3*Y[i])),3)
        end
    elseif g.occ===:cone
        let sa=_gm_atan((g.r2-g.r1)/g.height),(SinA,CosA)=_gm_sincos(sa),
            R=fma(v,SinA,g.r1),A1=R*c,A2=R*s,R1=SinA*c,R2=SinA*s,
            som1=ntuple(i->fma(A1,X[i],A2*Y[i]),3)
            return p,du,dv,.-som1,(0.0,0.0,0.0),
                   ntuple(i->fma(-R2,X[i],R1*Y[i]),3)
        end
    elseif g.occ===:torus
        let (sv,cv)=_gm_sincos(v),R1=g.r2*cv,R2=g.r2*sv,
            R=fma(g.r2,cv,g.r1),clamp_eps=10.0*(g.r2+g.r1)*eps(),
            A1=abs(R*c)<=clamp_eps ? 0.0 : R*c,
            A2=abs(R*s)<=clamp_eps ? 0.0 : R*s,
            A3=abs(R2*c)<=clamp_eps ? 0.0 : R2*c,
            A4=abs(R2*s)<=clamp_eps ? 0.0 : R2*s,
            A5=abs(R1*c)<=clamp_eps ? 0.0 : R1*c,
            A6=abs(R1*s)<=clamp_eps ? 0.0 : R1*s,
            som1=ntuple(i->fma(A1,X[i],A2*Y[i]),3)
            return p,du,dv,.-som1,
                   ntuple(i->fma(-A5,X[i],-(A6*Y[i]))-R2*Z[i],3),
                   ntuple(i->fma(A4,X[i],-(A3*Y[i])),3)
        end
    end
    throw(ArgumentError(
        "_occ_surface_d2: unsupported OCC surface kind $(g.occ)"))
end

# `gp_XYZ::SquareModulus` — compiles to fma(z,z, fma(x,x, y·y)).
@inline _sqlen(v) = fma(v[3],v[3],fma(v[1],v[1],v[2]*v[2]))
@inline _dot3(a,b) = fma(a[3],b[3],fma(a[1],b[1],a[2]*b[2]))

# `OCCFace::normal`: the `SVector3` cross product of the D1 tangents,
# multiplied by the reciprocal (a zero cross stays zero, like Gmsh's
# `SVector3::normalize`), then negated on a `TopAbs_REVERSED` face.
function _occ_surface_normal(g,u::Float64,v::Float64)
    _,du,dv=_occ_surface_d1(g,u,v)
    # `SVector3::crossprod` — a.y·b.z − a.z·b.y, a.z·b.x − a.x·b.z,
    # a.x·b.y − a.y·b.x (each second product rounded, the first fused).
    n=(fma(du[2],dv[3],-(dv[2]*du[3])),
       fma(du[3],dv[1],-(du[1]*dv[3])),
       fma(du[1],dv[2],-(dv[1]*du[2])))
    len=sqrt(_sqlen(n))
    len==0.0 && return (0.0,0.0,0.0)
    inv=1.0/len
    out=(n[1]*inv,n[2]*inv,n[3]*inv)
    (hasproperty(g,:reversed) && g.reversed) && (out=(.-out))
    return out
end

# `math_DirectPolynomialRoots`'s quadratic `Solve(A,B,C)` with the `Improve`
# Newton polish. Returns the (root1, root2) pair when `NbSolutions() == 2`,
# `nothing` otherwise.
function _occ_quadratic_roots(a::Float64,b::Float64,c::Float64)
    abs(a)<=1e-30 && return nothing   # Solve(b,c): at most one root
    epsd=3.0*eps()*fma(b,b,abs(4.0*a*c))
    disc=fma(b,b,-((4.0*a)*c))
    abs(disc)<=epsd && (disc=0.0)
    disc<0.0 && return nothing
    if disc==0.0
        r=_occ_improve_root(a,b,c,-0.5*b/a)
        return (r,r)
    end
    r0=b>0.0 ? -(b+sqrt(disc))/(2.0*a) : -(b-sqrt(disc))/(2.0*a)
    r0=_occ_improve_root(a,b,c,r0)
    r1=_occ_improve_root(a,b,c,c/(a*r0))
    return (r0,r1)
end

# `Improve(3, {A,B,C}, ini)`: up to nine Newton steps on `A·x² + B·x + C`
# evaluated by OCCT's Horner `Values`, keeping the better of polished/start.
function _occ_improve_root(a::Float64,b::Float64,c::Float64,ini::Float64)
    ini_val=fma(fma(a,ini,b),ini,c)
    sol=ini; val=ini_val
    for _ in 1:9
        val=fma(fma(a,sol,b),sol,c)
        der=fma(a,sol,fma(a,sol,b))
        abs(der)<=1e-30 && break
        delta=-val/der
        abs(delta)<=eps()*abs(sol) && break
        sol+=delta
    end
    return abs(val)<=abs(ini_val) ? sol : ini
end

# `BRepLProp_SLProps` (myCN=4, LinTol=eps) curvature block on the D2
# derivatives: `CSLib::Normal` gate, the ombilic shortcut, and the two
# `math_DirectPolynomialRoots` branches. Returns
# (defined, cmax, cmin, dirmax, dirmin); `defined=false` mirrors
# `IsCurvatureDefined()` failing (singular or parallel tangents).
function _occ_surface_curvatures(g,u::Float64,v::Float64)
    linTol=1e-12
    _,du,dv,duu,dvv,duv=_occ_surface_d2(g,u,v)
    # `CSLib::Normal(D1U, D1V, SinTol, status, Normal)` — `D1UW2` (the cross
    # square magnitude) has its own `gp::Resolution()` gate ahead of the
    # sin² ratio, and the ratio test is `<=`.
    d1u_mag=_sqlen(du); d1v_mag=_sqlen(dv)
    cross=_occ_cross(du,dv)
    (d1u_mag<=floatmin(Float64) || d1v_mag<=floatmin(Float64)) &&
        return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
    cross_sq=_sqlen(cross)
    cross_sq<=floatmin(Float64) &&
        return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
    cross_sq/(d1u_mag*d1v_mag)<=linTol*linTol &&
        return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
    normal=cross./sqrt(cross_sq)          # gp_Dir(gp_Vec): divide by modulus
    # `IsTangentUDefined`/`IsTangentVDefined` at orders 1 then 2
    tol2=linTol*linTol
    (d1u_mag>tol2 || _sqlen(duu)>tol2) ||
        return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
    (d1v_mag>tol2 || _sqlen(dvv)>tol2) ||
        return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
    E,F,G=d1u_mag,_dot3(du,dv),d1v_mag
    L,M,N=_dot3(normal,duu),_dot3(normal,duv),_dot3(normal,dvv)
    A=fma(E,M,-(F*L)); B=fma(E,N,-(G*L)); C=fma(F,N,-(G*M))
    max_abc=max(abs(A),abs(B),abs(C))
    dir(v)=v./sqrt(_sqlen(v))
    if max_abc<eps()                      # ombilic
        c=N/G
        return (true,c,c,dir(_occ_cross(du,normal)),dir(du))
    end
    A/=max_abc; B/=max_abc; C/=max_abc
    if abs(A)>eps()
        roots=_occ_quadratic_roots(A,B,C)
        roots===nothing &&
            return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
        let r1=roots[1],r2=roots[2]
            curv1=fma(fma(L,r1,2.0*M),r1,N)/fma(fma(E,r1,2.0*F),r1,G)
            curv2=fma(fma(L,r2,2.0*M),r2,N)/fma(fma(E,r2,2.0*F),r2,G)
            v1=ntuple(i->fma(r1,du[i],dv[i]),3)
            v2=ntuple(i->fma(r2,du[i],dv[i]),3)
        end
    elseif abs(C)>eps()
        roots=_occ_quadratic_roots(C,B,A)
        roots===nothing &&
            return (false,0.0,0.0,(0.0,0.0,0.0),(0.0,0.0,0.0))
        let r1=roots[1],r2=roots[2]
            curv1=fma(fma(N,r1,2.0*M),r1,L)/fma(fma(G,r1,2.0*F),r1,E)
            curv2=fma(fma(N,r2,2.0*M),r2,L)/fma(fma(G,r2,2.0*F),r2,E)
            v1=ntuple(i->fma(r1,dv[i],du[i]),3)
            v2=ntuple(i->fma(r2,dv[i],du[i]),3)
        end
    else
        curv1,curv2=L/E,N/G
        v1,v2=du,dv
    end
    if curv1<curv2
        return (true,curv2,curv1,dir(v2),dir(v1))
    end
    return (true,curv1,curv2,dir(v1),dir(v2))
end

# ── OCCT extrema/projection ports ──────────────────────────────────────────
#
# `OCCFace::_project`/`OCCEdge::_project` wrap `GeomAPI_ProjectPointOnSurf`/
# `ProjectPointOnCurve`. For the analytic surface/curve kinds these run
# `Extrema_ExtPElS`/`ExtPElC`'s closed-form extremum solvers, then keep only
# candidates whose (period-mapped) parameters fall inside the projector's
# bounds (`Extrema_ExtPS::TreatSolution`, `Extrema_GExtPC`'s post-check).
# The projector bounds are the face/edge parameter ranges padded by
# `max(range·1e-8, 1e-12)` on non-periodic directions; the accepted
# parameter is returned unclamped. Empty candidate lists fall back to the
# generic `GFace::XYZtoUV`/`GEdge::XYZToU` Newton scans, ported below.

const _OCC_CONFUSION=1e-7           # Precision::Confusion
const _OCC_EXT_EPS=eps(2π)          # ExtPElS_MyEps = Epsilon(2·π)
const _OCC_ANGULAR=1e-12            # Precision::Angular

@inline _sqdist(a::NTuple{3,Float64},b::NTuple{3,Float64}) =
    _sqlen(_arc_sub(a,b))

# `gp_Dir::AngleWithRef` — `acos` between ±45°, `asin` outside; the sign
# comes from `(self × other)·vref`. All three vectors normalize first (the
# `gp_Dir` conversions).
function _occ_angle_with_ref(a,b,vref)
    da,db,dv=_occ_dir(a),_occ_dir(b),_occ_dir(vref)
    xyz=_occ_cross(da,db)
    cosinus=_dot3(da,db)
    sinus=_occ_modulus(xyz)
    ang=(cosinus > -0.70710678118655 && cosinus < 0.70710678118655) ?
        _gm_acos(cosinus) : cosinus<0.0 ? π-_gm_asin(sinus) : _gm_asin(sinus)
    return _dot3(xyz,dv)>=0.0 ? ang : -ang
end

# `gp_Dir::Angle` — the unsigned angle in [0,π].
function _occ_dir_angle(a,b)
    da,db=_occ_dir(a),_occ_dir(b)
    cosinus=_dot3(da,db)
    (cosinus > -0.70710678118655 && cosinus < 0.70710678118655) &&
        return _gm_acos(cosinus)
    sinus=_occ_modulus(_occ_cross(da,db))
    return cosinus<0.0 ? π-_gm_asin(sinus) : _gm_asin(sinus)
end

# `ElCLib::InPeriod(u, a, b)` — maps u into [a, a+(b−a)] by `ceil`.
function _occ_in_period(u::Float64,a::Float64,b::Float64)
    period=b-a
    period<eps(b) && return u
    return max(a,u+period*ceil((a-u)/period))
end

# `ElCLib::AdjustPeriodic(ufirst, ulast, preci, u1, u2)` — shifts u1 into
# [ufirst, ufirst+p) (or one period below ulast when within `preci`), then
# u2 into [u1, u1+p] keeping `u2-u1 >= preci`.
function _occ_adjust_periodic(ufirst,ulast,preci,u1,u2)
    period=ulast-ufirst
    period<eps(ulast) && return u1,u2
    u1-=floor((u1-ufirst)/period)*period
    ulast-u1<preci && (u1-=period)
    u2-=floor((u2-u1)/period)*period
    u2-u1<preci && (u2+=period)
    return u1,u2
end

# `Extrema_ExtPElS::Perform(P, gp_Cylinder, Tol)` — the two angular extrema
# at the query's axial coordinate; `nothing` on the axis (continuum).
function _occ_extpe_cylinder(g,P)
    O,X,Z=g.center,g.X,g.axis
    myZ=_occ_cross(X,_occ_frame_y(g))
    V=_dot3(_arc_sub(P,O),Z)
    Pp=_arc_sub(P,_occ_mul(Z,V))
    OPp=_arc_sub(Pp,O)
    _occ_modulus(OPp)<_OCC_CONFUSION && return nothing
    U1=_occ_angle_with_ref(X,OPp,myZ)
    (-_OCC_EXT_EPS < U1 < _OCC_EXT_EPS) && (U1=0.0)
    U2=U1+π
    U1<0.0 && (U1+=2π)
    ps1=_occ_surface_point(g,U1,V); ps2=_occ_surface_point(g,U2,V)
    return [(U1,V,ps1,_sqdist(ps1,P)),(U2,V,ps2,_sqdist(ps2,P))]
end

# `Extrema_ExtPElS::Perform(P, gp_Sphere, Tol)` — (U1,V) minimum and
# (U2,−V) maximum; `nothing` at the center.
function _occ_extpe_sphere(g,P)
    O,X,Z=g.center,g.X,g.axis
    OP=_arc_sub(P,O)
    _sqlen(OP)<_OCC_CONFUSION*_OCC_CONFUSION && return nothing
    Zp=_dot3(OP,Z)
    Pp=_arc_sub(P,_occ_mul(Z,Zp))
    OPp=_arc_sub(Pp,O)
    if _sqlen(OPp)<_OCC_CONFUSION*_OCC_CONFUSION
        U1=U2=0.0
        V=Zp<0.0 ? -π/2 : π/2
    else
        myZ=_occ_cross(X,_occ_frame_y(g))
        U1=_occ_angle_with_ref(X,OPp,myZ)
        (-_OCC_EXT_EPS < U1 < _OCC_EXT_EPS) && (U1=0.0)
        U2=U1+π
        U1<0.0 && (U1+=2π)
        V=_occ_dir_angle(OP,OPp)
        Zp<0.0 && (V=-V)
    end
    ps1=_occ_surface_point(g,U1,V); ps2=_occ_surface_point(g,U2,-V)
    return [(U1,V,ps1,_sqdist(ps1,P)),(U2,-V,ps2,_sqdist(ps2,P))]
end

# `Extrema_ExtPElS::Perform(P, gp_Cone, Tol)` — the apex-shortcut and the
# two meridian extrema; `nothing` on the axis below the apex side.
function _occ_extpe_cone(g,P)
    O,X,Z=g.center,g.X,g.axis
    A=_gm_atan((g.r2-g.r1)/g.height)
    M=_add3(O,_occ_mul(Z,-g.r1/_gm_tan(A)))     # gp_Cone::Apex
    myZ=_occ_cross(X,_occ_frame_y(g))
    MP=_arc_sub(P,M)
    L2=_sqlen(MP)
    Vm=-g.r1/_gm_sin(A)
    L2<_OCC_CONFUSION*_OCC_CONFUSION && return [(0.0,Vm,M,L2)]
    DirZ=_sqdist(M,O)<_OCC_CONFUSION*_OCC_CONFUSION ?
        (A<0.0 ? _occ_mul(Z,-1.0) : Z) : _arc_sub(O,M)
    Zp=_dot3(_arc_sub(P,O),Z)
    Pp=_add3(P,_occ_mul(Z,-Zp))
    OPp=_arc_sub(Pp,O)
    _sqlen(OPp)<_OCC_CONFUSION*_OCC_CONFUSION && return nothing
    Same=_dot3(DirZ,MP)>=0.0
    U1=_occ_angle_with_ref(X,OPp,myZ)
    (-_OCC_EXT_EPS < U1 < _OCC_EXT_EPS) && (U1=0.0)
    Same || (U1+=π)
    U2=U1+π
    U1<0.0 && (U1+=2π)
    U2>2π && (U2-=2π)
    B=_occ_dir_angle(MP,DirZ)
    Aa=abs(A)
    L=sqrt(L2)
    if Same
        V1=L*_gm_cos(B-Aa); V2=L*_gm_cos(B+Aa)
    else
        B=π-B
        V1=-L*_gm_cos(B-Aa); V2=-L*_gm_cos(B+Aa)
    end
    Sense=_dot3(Z,_occ_dir(DirZ))
    V1=V1*Sense+Vm; V2=V2*Sense+Vm
    ps1=_occ_surface_point(g,U1,V1); ps2=_occ_surface_point(g,U2,V2)
    return [(U1,V1,ps1,_sqdist(ps1,P)),(U2,V2,ps2,_sqdist(ps2,P))]
end

# `Extrema_ExtPElS::Perform(P, gp_Torus, Tol)` — the four meridian-circle
# extrema; `nothing` on the axis or on a tube-center circle.
function _occ_extpe_torus(g,P)
    O,X,Z=g.center,g.X,g.axis
    myZ=_occ_cross(X,_occ_frame_y(g))
    tol2=_OCC_CONFUSION*_OCC_CONFUSION
    Pp=_add3(P,_occ_mul(Z,-_dot3(_arc_sub(P,O),Z)))
    OPp=_arc_sub(Pp,O)
    R2=_sqlen(OPp)
    R2<tol2 && return nothing
    U1=_occ_angle_with_ref(X,OPp,myZ)
    (-_OCC_EXT_EPS < U1 < _OCC_EXT_EPS) && (U1=0.0)
    U2=U1+π
    U1<0.0 && (U1+=2π)
    R=sqrt(R2)
    # `OPp.Divided(R)` is componentwise division, not a reciprocal multiply.
    OO1=_occ_mul(ntuple(i->OPp[i]/R,3),g.r1)
    O1=_add3(O,OO1); O2=_add3(O,_occ_mul(OO1,-1.0))
    _sqdist(O1,P)<tol2 && return nothing
    _sqdist(O2,P)<tol2 && return nothing
    V1=_occ_angle_with_ref(OPp,_arc_sub(P,O1),_occ_cross(OPp,Z))
    (-_OCC_EXT_EPS < V1 < _OCC_EXT_EPS) && (V1=0.0)
    NOp=_occ_mul(OPp,-1.0)
    V2=_occ_angle_with_ref(NOp,_arc_sub(O2,P),_occ_cross(NOp,Z))
    (-_OCC_EXT_EPS < V2 < _OCC_EXT_EPS) && (V2=0.0)
    V1<0.0 && (V1+=2π)
    V2<0.0 && (V2+=2π)
    V1p=V1+π; V2p=V2+π
    ps1=_occ_surface_point(g,U1,V1); ps2=_occ_surface_point(g,U1,V1p)
    ps3=_occ_surface_point(g,U2,V2); ps4=_occ_surface_point(g,U2,V2p)
    return [(U1,V1,ps1,_sqdist(ps1,P)),(U1,V1p,ps2,_sqdist(ps2,P)),
            (U2,V2,ps3,_sqdist(ps3,P)),(U2,V2p,ps4,_sqdist(ps4,P))]
end

# `ElSLib::PlaneParameters` — `gp_Trsf::SetTransformation(Ax3)` maps P into
# the face frame: rows of M are the X/Y/Z directions, `loc = −(M·Loc)`, and
# `Transformed` applies `M·P + loc` (`gp_XYZ::Multiply` + `Add`).
function _occ_plane_uv(g,P)
    M=(g.X,_occ_frame_y(g),g.axis)
    loc=_occ_mul(_occ_matvec(M,g.center),-1.0)
    ploc=_add3(_occ_matvec(M,P),loc)
    return ploc[1],ploc[2]
end

# `Extrema_ExtPElS::Perform(P, gp_Pln, Tol)` — always done with a single
# extremum: the orthogonal projection `Pp = P − (OP·Z)·Z` and the
# `ElSLib::Parameters` frame uv.
function _occ_extpe_plane(g,P)
    O,Z=g.center,g.axis
    v0=_dot3(_arc_sub(P,O),Z)
    Pp=_add3(P,_occ_mul(Z,-v0))
    u,v=_occ_plane_uv(g,P)
    return [(u,v,Pp,_sqdist(Pp,P))]
end

# `Extrema_ExtPS::TreatSolution` per candidate — period-map into
# [inf, inf+period), pull back one period when still outside
# [inf−tol, sup+tol], then the bounds acceptance test. `u` is periodic on
# the four revolution kinds (`UPeriod` = 2π); `v` is periodic on the torus;
# the plane is periodic in neither.
function _occ_treat_solutions(g,cands,umin,umax,vmin,vmax,tolu,tolv)
    accepted=Tuple{Float64,Float64,NTuple{3,Float64},Float64}[]
    for (U,V,xyz,sqd) in cands
        if g.occ!==:plane
            U=_occ_in_period(U,umin,umin+2π)
            U>umax+tolu && (U-=2π)
            U<umin-tolu && (U+=2π)
        end
        if g.occ===:torus
            V=_occ_in_period(V,vmin,vmin+2π)
            V>vmax+tolv && (V-=2π)
            V<vmin-tolv && (V+=2π)
        end
        if (umin-U)<=tolu && (U-umax)<=tolu &&
           (vmin-V)<=tolv && (V-vmax)<=tolv
            push!(accepted,(U,V,xyz,sqd))
        end
    end
    return accepted
end

# The `Extrema_ExtPElS` candidates after `TreatSolution` filtering against
# the caller's bounds, or `nothing` when the solver is not done (degenerate
# query points).
function _occ_surface_extrema(
    g,P,umin,umax,vmin,vmax,tolu,tolv)
    cands=g.occ===:cylinder ? _occ_extpe_cylinder(g,P) :
          g.occ===:sphere   ? _occ_extpe_sphere(g,P) :
          g.occ===:cone     ? _occ_extpe_cone(g,P) :
          g.occ===:torus    ? _occ_extpe_torus(g,P) :
          g.occ===:plane    ? _occ_extpe_plane(g,P) :
          throw(ArgumentError(
              "_occ_surface_extrema: unsupported OCC surface kind $(g.occ)"))
    cands===nothing && return nothing
    return _occ_treat_solutions(g,cands,umin,umax,vmin,vmax,tolu,tolv)
end

# `OCCFace::_project` — padded face bounds, `Precision::Confusion`
# tolerances, first-minimum selection. Returns (u, v, xyz) or `nothing`.
function _occ_surface_project(g,P)
    lo,hi=_occ_surface_bounds(g)
    umin,vmin=lo[1],lo[2]; umax,vmax=hi[1],hi[2]
    # `OCCFace`'s constructor pads each non-periodic direction of the
    # projector bounds by max(|range|·1e-8, 1e-12): u is padded only on the
    # plane, v on everything but the torus.
    if g.occ===:plane
        du=umax-umin
        ut=max(abs(du)*1e-8,1e-12)
        umin-=ut; umax+=ut
    end
    if g.occ!==:torus
        dv=vmax-vmin
        vt=max(abs(dv)*1e-8,1e-12)
        vmin-=vt; vmax+=vt
    end
    acc=_occ_surface_extrema(
        g,P,umin,umax,vmin,vmax,_OCC_CONFUSION,_OCC_CONFUSION)
    (acc===nothing || isempty(acc)) && return nothing
    best=acc[1]
    for cand in Iterators.drop(acc,1)
        cand[4]<best[4] && (best=cand)
    end
    return (best[1],best[2],best[3])
end

# ── BRepClass 2-D face classifier ────────────────────────────────────────────
#
# Port of the OCCT 7.9.3 `BRepClass_FaceClassifier` chain restricted to the
# pcurve kinds the OCC materializers emit (`lin2d`, `circ2d`) and to
# `UseBndBox == false` (Gmsh's `OCCFace::containsPoint`/`containsParam`
# invocation). The mirrored pieces:
#   `BRepClass_FClassifier` / `TopClass_FaceClassifier` — the driver;
#   `BRepClass_FaceExplorer` — probing-segment construction;
#   `TopClass_Classifier2d` — closest-intersection state machine;
#   `TopTrans_CurveTransition` — complex transition at edge ends;
#   `BRepClass_Intersector` — `CheckOn` + `Geom2dInt_GInter` + `CheckSkip`;
#   `IntCurve_IntConicConic` — the lin×lin / lin×circ solvers;
#   `Extrema_ExtPC2d`/`ExtPElC2d` — the `CheckOn` projector;
#   `Geom2dLProp_CLProps2d` — `LocalGeometry` tangent/normal/curvature;
#   `BRepTools::UVBounds` — pcurve-derived face bounds.
# TopAbs states are :in/:out/:on/:unknown; IntRes2d positions
# :head/:middle/:end; transition types :in/:out/:touch/:undecided; situations
# :inside/:outside/:unknown; orientations :forward/:reversed/:internal/:external.

const _OCC_PCONFUSION=1e-9          # Precision::PConfusion()
const _OCC_PINTERSECTION=1e-11      # Precision::PIntersection()
const _OCC_SQPARTOL=1e-18           # Precision::PConfusion()²
const _OCC_TOLANG=1e-8              # TOLERANCE_ANGULAIRE
const _OCC_GAPCHECK=0.1             # BRepClass_FaceExplorer::myMaxTolerance
const _OCC_PROBING=(0.123,0.7,0.2111)  # Probing_Start/End/Step
const _OCC_INFINITE=2e100           # Precision::Infinite()

# `IntRes2d_Transition` — (tan, pos, tr, situ, oppos).
_itr()=(tan=true,pos=:middle,tr=:undecided,situ=:unknown,oppos=false)
_itr(pos::Symbol)=(tan=true,pos=pos,tr=:undecided,situ=:unknown,oppos=false)
_itr(tan::Bool,pos::Symbol,tr::Symbol)=
    (tan=tan,pos=pos,tr=tr,situ=:unknown,oppos=false)
_itr(tan::Bool,pos::Symbol,situ::Symbol,oppos::Bool)=
    (tan=tan,pos=pos,tr=:touch,situ=situ,oppos=oppos)

# `IntRes2d_IntersectionPoint` — parameters/transitions already in caller-curve
# order (the `ReversedFlag` swap applied at construction).
_ipt(P,u1,u2,t1,t2,rev)=rev ? (pt=P,p1=u2,p2=u1,t1=t2,t2=t1) :
                              (pt=P,p1=u1,p2=u2,t1=t1,t2=t2)

# `IntRes2d_IntersectionSegment`.
function _iseg(p1,p2,oppos,rev)
    rev && oppos && return (p1=p2,p2=p1,first=true,last=true,oppos=oppos)
    return (p1=p1,p2=p2,first=true,last=true,oppos=oppos)
end
function _iseg_half(P,first,oppos,rev)
    if rev && oppos
        first=!first
    end
    return (p1=first ? P : (pt=(0.0,0.0),p1=Inf,p2=Inf,t1=_itr(),t2=_itr()),
            p2=first ? (pt=(0.0,0.0),p1=Inf,p2=Inf,t1=_itr(),t2=_itr()) : P,
            first=first,last=!first,oppos=oppos)
end
_iseg(oppos::Bool)=(p1=(pt=(0.0,0.0),p1=Inf,p2=Inf,t1=_itr(),t2=_itr()),
                    p2=(pt=(0.0,0.0),p1=Inf,p2=Inf,t1=_itr(),t2=_itr()),
                    first=false,last=false,oppos=oppos)

# `IntRes2d_Domain` — (fpt,fpar,ftol,hasf, lpt,lpar,ltol,hasl, e0,ep,hasper).
# `LimitInfinite` clamps magnitudes beyond `Precision::Infinite()` (2e100).
_liminf(v)=abs(v)>_OCC_INFINITE ? (v>0 ? _OCC_INFINITE : -_OCC_INFINITE) : v
function _idom(fpt,fpar,ftol,lpt,lpar,ltol)
    return (fpt=(_liminf(fpt[1]),_liminf(fpt[2])),fpar=_liminf(fpar),ftol=ftol,
            hasf=true,
            lpt=(_liminf(lpt[1]),_liminf(lpt[2])),lpar=_liminf(lpar),ltol=ltol,
            hasl=true,e0=0.0,ep=0.0,hasper=false)
end
function _idom(pt,par,tol,first)
    return (fpt=first ? (_liminf(pt[1]),_liminf(pt[2])) : (0.0,0.0),
            fpar=first ? _liminf(par) : 0.0,ftol=first ? tol : 0.0,hasf=first,
            lpt=first ? (0.0,0.0) : (_liminf(pt[1]),_liminf(pt[2])),
            lpar=first ? 0.0 : _liminf(par),ltol=first ? 0.0 : tol,hasl=!first,
            e0=0.0,ep=0.0,hasper=false)
end
_idom()=(fpt=(0.0,0.0),fpar=0.0,ftol=0.0,hasf=false,
         lpt=(0.0,0.0),lpar=0.0,ltol=0.0,hasl=false,e0=0.0,ep=0.0,hasper=false)
_idom_per(d,e0,ep)=(fpt=d.fpt,fpar=d.fpar,ftol=d.ftol,hasf=d.hasf,
                   lpt=d.lpt,lpar=d.lpar,ltol=d.ltol,hasl=d.hasl,
                   e0=e0,ep=ep,hasper=true)

# ── gp 2-D helpers ───────────────────────────────────────────────────────────
@inline _dot2(a,b)=fma(a[1],b[1],a[2]*b[2])
@inline _cross2(a,b)=a[1]*b[2]-a[2]*b[1]      # gp_XY::Crossed
@inline _sqmag2(a)=fma(a[1],a[1],a[2]*a[2])
@inline _mag2(a)=sqrt(_sqmag2(a))
@inline _sub2(a,b)=(a[1]-b[1],a[2]-b[2])
@inline _add2(a,b)=(a[1]+b[1],a[2]+b[2])
@inline _mul2(a,s)=(a[1]*s,a[2]*s)
@inline _sqdist2(a,b)=_sqmag2(_sub2(a,b))
@inline _dist2(a,b)=sqrt(_sqdist2(a,b))
@inline _dir2(v)=let m=_mag2(v); (v[1]/m,v[2]/m) end   # gp_Dir2d ctor

# gp_Lin2d — (loc, dir).
@inline _lin2d_dist(L,P)=abs(_cross2(_sub2(P,L[1]),L[2]))
@inline _lin2d_sqdist(L,P)=let d=_cross2(_sub2(P,L[1]),L[2]); d*d end
@inline _lin2d_coefs(L)=(L[2][2],-L[2][1],-(L[2][2]*L[1][1]-L[2][1]*L[1][2]))
@inline _lin2d_param(L,P)=_dot2(_sub2(P,L[1]),L[2])       # ElCLib::Parameter
@inline _lin2d_value(L,t)=(fma(t,L[2][1],L[1][1]),fma(t,L[2][2],L[1][2]))
@inline _lin2d_d1(L,t)=(_lin2d_value(L,t),L[2])           # ElCLib::LineD1

# gp_Circ2d — (loc, xdir, r) with derived ydir = (−x.y, x.x) (always direct).
@inline _circ2d_axis(pc)=let (o,x,r)=pc; ((o,x,(-x[2],x[1])),r) end
@inline _circ2d_value(ax,r,t)=_occ_circ2d_eval(ax[1],ax[2],r,t)
@inline _circ2d_d1(ax,r,t)=let (A2_,A1_)=_gm_sincos(t),A1=r*A1_,A2=r*A2_,x=ax[2],y=ax[3]
    (_occ_circ2d_eval(ax[1],x,r,t),
     (fma(-A2,x[1],A1*y[1]),fma(-A2,x[2],A1*y[2])))
end
@inline _circ2d_d2(ax,r,t)=let (A2_,A1_)=_gm_sincos(t),A1=r*A1_,A2=r*A2_,x=ax[2],y=ax[3]
    (_occ_circ2d_eval(ax[1],x,r,t),
     (fma(-A2,x[1],A1*y[1]),fma(-A2,x[2],A1*y[2])),
     (-fma(A1,x[1],A2*y[1]),-fma(A1,x[2],A2*y[2])))
end
# ElCLib::CircleParameter — `XDirection.Angle(P−loc)` signed by `IsDirect`
# (always true here), with the −1e-16/0 snap.
function _circ2d_param(ax,P)
    Teta=_vec2d_angle(ax[2],_sub2(P,ax[1]))
    Teta=_cross2(ax[2],ax[3])>=0.0 ? Teta : -Teta
    if Teta<-1e-16
        Teta+=_OCC_TWO_PI
    elseif Teta<0
        Teta=0.0
    end
    return Teta
end

# gp_Vec2d::Angle — the acos/asin branch formula; both operands need not be
# unit (the norms are divided out).
function _vec2d_angle(a,b)
    na,nb=_mag2(a),_mag2(b)
    D=na*nb
    Cosinus=_dot2(a,b)/D; Sinus=_cross2(a,b)/D
    if Cosinus>-0.70710678118655 && Cosinus<0.70710678118655
        return Sinus>0.0 ? _gm_acos(Cosinus) : -_gm_acos(Cosinus)
    elseif Cosinus>0.0
        return _gm_asin(Sinus)
    else
        return Sinus>0.0 ? π-_gm_asin(Sinus) : -π-_gm_asin(Sinus)
    end
end

# `Precision::IsInfinite` — |v| ≥ 0.5·Infinite() = 1e100.
@inline _occ_isinf(v)=abs(v)>=0.5*_OCC_INFINITE
@inline _occ_isinfpos(v)=v>=0.5*_OCC_INFINITE
@inline _occ_isinfneg(v)=v<=-0.5*_OCC_INFINITE
# `Epsilon(v)` — the nextafter gap in the direction of sign(v); negative for
# v < 0, as in OCCT (`nextafter` toward −Inf).
@inline _occ_epsilon(v)=v>0.0 ? nextfloat(v)-v :
    v<0.0 ? prevfloat(v)-v : nextfloat(0.0)

# ── pcurve adaptor helpers ───────────────────────────────────────────────────
_occ_pcurve_kind(pc)=hasproperty(pc,:lin2d) ? :lin2d :
    hasproperty(pc,:circ2d) ? :circ2d :
    throw(ArgumentError("_occ_pcurve_kind: unsupported pcurve record"))
@inline _occ_pcurve_isperiodic(pc)=_occ_pcurve_kind(pc)===:circ2d
# `Geom2dAdaptor_Curve::Resolution(Ruv)` — line: Ruv; circle: 2·asin(Ruv/2r)
# or 2π when r ≤ Ruv/2.
function _occ_pcurve_resolution(pc,ruv)
    _occ_pcurve_kind(pc)===:lin2d && return ruv
    r=pc.circ2d[3]
    return r>ruv/2.0 ? 2.0*_gm_asin(ruv/(2.0*r)) : _OCC_TWO_PI
end
@inline _occ_pcurve_value(pc,t)=_occ_pcurve_d0(pc,t)
function _occ_pcurve_d2eval(pc,t)
    p,d1,d2=_occ_pcurve_d2(pc,t)
    return p,d1,d2
end

# `BRep_Tool::CurveOnSurface(E,F)` — `PCurve2` only when the edge is reversed
# AND carries two pcurves on the face (`IsCurveOnClosedSurface`); otherwise
# `PCurve`. Returns (pcurve_record, deb, fin); `nothing` when no pcurve exists.
function _occ_face_pcurve(m::GeoModel,g,etag::Int)
    ent=get(g.pcurves,abs(etag),nothing)
    ent===nothing && return nothing
    pc=(etag<0 && ent.rev!==nothing) ? ent.rev : ent.fwd
    rec=m.curve_geometry[abs(etag)]
    return pc,rec.t0,rec.t1
end

# ── IntRes2d domain helpers ─────────────────────────────────────────────────
# `IntImpParGen::DeterminePosition` — point-distance based Head/End/Middle.
function _occ_determine_position(D,P,param)
    pos=:middle
    if D.hasf && _dist2(P,D.fpt)<=D.ftol
        pos=:head
    end
    if D.hasl && _dist2(P,D.lpt)<=D.ltol
        if pos===:head
            abs(param-D.lpar)<abs(param-D.fpar) && (pos=:end)
        else
            pos=:end
        end
    end
    return pos
end

# `DomainIntersection` — clip [Uinf,Usup] against a domain with end
# tolerances; empty when the interval lies entirely outside a bound.
function _occ_domain_intersection(D,Uinf,Usup)
    if D.hasf
        Usup < D.fpar-D.ftol && return (1.0,-1.0,:middle,:middle)
        if Uinf > D.fpar+D.ftol
            resinf,posinf=Uinf,:middle
        else
            resinf,posinf=D.fpar,:head
        end
    else
        resinf,posinf=Uinf,:middle
    end
    if D.hasl
        Uinf > D.lpar+D.ltol && return (1.0,-1.0,:middle,:middle)
        if Usup < D.lpar-D.ltol
            ressup,possup=Usup,:middle
        else
            ressup,possup=D.lpar,:end
        end
    else
        ressup,possup=Usup,:middle
    end
    if resinf>ressup
        if possup===:middle
            ressup=resinf
        else
            resinf=ressup
        end
    end
    return resinf,ressup,posinf,possup
end

# `FindPositionLL` — snap Param to a bound inside its tolerance; returns
# (pos, adjusted param). Head wins unless the End is strictly closer.
function _occ_findpositionll(param,D)
    dpar=Inf; pos=:middle; res=param
    if D.hasf
        dpar=abs(param-D.fpar)
        if dpar<=D.ftol
            res=D.fpar; pos=:head
        end
    end
    if D.hasl
        d2=abs(param-D.lpar)
        if d2<=D.ltol && (pos===:middle || d2<dpar)
            res=D.lpar; pos=:end
        end
    end
    return pos,res
end

# ── PeriodicInterval / Interval ──────────────────────────────────────────────
const _PIPI=_OCC_TWO_PI
_pint()=(binf=0.0,bsup=0.0,isnull=true)
_pint(a,b)=let bi,bs
    if b-a<_PIPI
        bi,bs=_pint_normalize(a,b)
    else
        bi,bs=a,b
    end
    (binf=bi,bsup=bs,isnull=false)
end
function _pint_normalize(binf,bsup)
    while binf>_PIPI; binf-=_PIPI end
    while binf<0.0; binf+=_PIPI end
    while bsup<binf; bsup+=_PIPI end
    while bsup>=binf+_PIPI; bsup-=_PIPI end
    return binf,bsup
end
_pint_len(p)=p.isnull ? -100.0 : abs(p.bsup-p.binf)
function _pint_complement(p)
    p.isnull && return p
    t=p.binf; binf=p.bsup; bsup=t+_PIPI
    if binf>_PIPI
        binf-=_PIPI; bsup-=_PIPI
    end
    return (binf=binf,bsup=bsup,isnull=false)
end
_pint_set(a,b)=_pint(a,b)
# `PeriodicInterval(Domain)` — first/last params or the (−1, 20) defaults.
function _pint(D)
    binf=D.hasf ? D.fpar : -1.0
    bsup=D.hasl ? D.lpar : 20.0
    return (binf=binf,bsup=bsup,isnull=false)
end

# `Interval` (bounded-domain interval) — (binf,bsup,hf,hl,isnull).
_ivl()=(binf=0.0,bsup=0.0,hf=false,hl=false,isnull=true)
_ivl(a,b)=(binf=min(a,b),bsup=max(a,b),hf=true,hl=true,isnull=false)
_ivl_raw(a,hf,b,hl)=(binf=a,bsup=b,hf=hf,hl=hl,isnull=false)
function _ivl(D)
    hf,hl=D.hasf,D.hasl
    binf=hf ? D.fpar-D.ftol : 0.0
    bsup=hl ? D.lpar+D.ltol : 0.0
    return (binf=binf,bsup=bsup,hf=hf,hl=hl,isnull=false)
end
_ivl_len(i)=i.isnull ? -1.0 : abs(i.bsup-i.binf)
function _ivl_intersect_bounded(self,other)
    (self.isnull || other.isnull) && return _ivl()
    if !(self.hf || self.hl)
        return _ivl(other.binf,other.bsup)
    end
    if self.hf
        other.bsup<self.binf && return _ivl()
        a=other.binf<self.binf ? self.binf : other.binf
    else
        a=other.binf
    end
    if self.hl
        other.binf>self.bsup && return _ivl()
        b=other.bsup>self.bsup ? self.bsup : other.bsup
    else
        b=other.bsup
    end
    return _ivl(a,b)
end

# `PeriodicInterval::FirstIntersection` — mutates PInter (normalizes it
# toward this interval); returns (result, mutated PInter).
function _pint_first_intersection(self,P)
    (P.isnull || self.isnull) && return _pint(),P
    _pint_len(self)>=_PIPI && return _pint(P.binf,P.bsup),P
    _pint_len(P)>=_PIPI && return _pint(self.binf,self.bsup),P
    binf,bsup=P.binf,P.bsup
    if bsup<=self.binf
        while binf<=self.binf && bsup<=self.binf
            binf+=_PIPI; bsup+=_PIPI
        end
    end
    if binf>=self.bsup
        while binf>=self.bsup && bsup>=self.bsup
            binf-=_PIPI; bsup-=_PIPI
        end
    end
    (bsup<self.binf || binf>self.bsup) && return _pint(),P
    a=max(binf,self.binf); b=min(bsup,self.bsup)
    return _pint(a,b),(binf=binf,bsup=bsup,isnull=false)
end
function _pint_second_intersection(self,P)
    (P.isnull || self.isnull || _pint_len(self)>=_PIPI ||
     _pint_len(P)>=_PIPI) && return _pint()
    inf,sup=P.binf+_PIPI,P.bsup+_PIPI
    if inf>self.bsup
        inf,sup=P.binf-_PIPI,P.bsup-_PIPI
    end
    (sup<self.binf || inf>self.bsup) && return _pint()
    return _pint(max(inf,self.binf),min(sup,self.bsup))
end

# `IntRes2d_Intersection::Insert` — sorted insertion by ParamOnFirst with
# `PARAMEQUAL` (|a−b| < 1e-8) deduplication on equal params AND transitions.
_parmeq(a,b)=abs(a-b)<1e-8
function _tr_equal(T1,T2)
    T1.pos!==T2.pos && return false
    T1.tr!==T2.tr && return false
    T1.tr===:touch || return true
    T1.tan!==T2.tan && return false
    T1.situ!==T2.situ && return false
    return T1.oppos===T2.oppos
end
function _ires_insert!(lpnt,P)
    n=length(lpnt)
    n==0 && (push!(lpnt,P); return)
    u=P.p1; b=n+1
    for i in 1:n
        Pnti=lpnt[i]; ui=Pnti.p1
        ui>=u && (b=i)
        if _parmeq(ui,u) && _parmeq(P.p2,Pnti.p2) &&
           _tr_equal(P.t1,Pnti.t1) && _tr_equal(P.t2,Pnti.t2)
            b=0
        end
        b==0 && break
        ui>=u && break
    end
    b>n ? push!(lpnt,P) : b>0 && insert!(lpnt,b,P)
    return
end

# `SegmentToPoint` — merge a degenerate segment into one of its endpoints.
function _occ_segmenttopoint(Pa,T1a,T2a,Pb,T1b,T2b)
    (T1b.pos===:middle && T2b.pos===:middle) && return Pa
    (T1a.pos===:middle && T2a.pos===:middle) && return Pb
    t1=T1a; t2=T2a; u1=Pa.p1; u2=Pa.p2
    if t1.pos===:middle
        t1=(tan=t1.tan,pos=T1b.pos,tr=t1.tr,situ=t1.situ,oppos=t1.oppos)
        u1=Pb.p1
    end
    if t2.pos===:middle
        t2=(tan=t2.tan,pos=T2b.pos,tr=t2.tr,situ=t2.situ,oppos=t2.oppos)
        u2=Pb.p2
    end
    return _ipt(Pa.pt,u1,u2,t1,t2,false)
end

# `Determine_Transition_LC` — line/circle transition from tangents and the
# conic normal; mutates Tan1 (normalizes it). Norm2 is the zero vector for
# the line side.
function _occ_transition_lc(Pos1,Tan1,Norm1,Pos2,Tan2,Norm2)
    sgn=_cross2(Tan1,Tan2)
    norm=_mag2(Tan1)*_mag2(Tan2)
    if abs(sgn)<=_OCC_TOLANG*norm
        opos=_dot2(Tan1,Tan2)<0.0
        Tan1=_dir2(Tan1)
        Norm=(-Tan1[2],Tan1[1])
        Val1=_dot2(Norm,Norm1); Val2=_dot2(Norm,Norm2)
        if abs(Val1-Val2)<=floatmin(Float64)          # gp::Resolution()
            T1=_itr(true,Pos1,:unknown,opos); T2=_itr(true,Pos2,:unknown,opos)
        elseif Val2>Val1
            T2=_itr(true,Pos2,:inside,opos)
            T1=_itr(true,Pos1,opos ? :inside : :outside,opos)
        else
            T2=_itr(true,Pos2,:outside,opos)
            T1=_itr(true,Pos1,opos ? :outside : :inside,opos)
        end
    elseif sgn<0.0
        T1=_itr(false,Pos1,:in); T2=_itr(false,Pos2,:out)
    else
        T1=_itr(false,Pos1,:out); T2=_itr(false,Pos2,:in)
    end
    return T1,T2
end

# `NormalizeOnCircleDomain` — wrap Param into [First, Last] by 2π shifts.
function _occ_normalize_on_circle_domain(param,D)
    while param<D.fpar; param+=_PIPI end
    while param>D.lpar; param-=_PIPI end
    return param
end

# `LineLineGeometricIntersection` — (U1,U2,sindemi,nbsol). nbsol: 1
# crossing, 2 coincident within Tol, 0 parallel disjoint.
function _occ_linlin_geometric(L1,L2,Tol)
    U1x,U1y=L1[2]; U2x,U2y=L2[2]
    Uo21x,Uo21y=L2[1][1]-L1[1][1],L2[1][2]-L1[1][2]
    D=U1y*U2x-U1x*U2y
    if abs(D)<_OCC_TOLANG
        D=U1y*Uo21x-U1x*Uo21y
        return 0.0,0.0,0.0,(abs(D)<=Tol ? 2 : 0)
    end
    U1=(Uo21y*U2x-Uo21x*U2y)/D
    U2=(Uo21y*U1x-Uo21x*U1y)/D
    D<0.0 && (D=-D)
    D>1.0 && (D=1.0)
    return U1,U2,_gm_sin(0.5*_gm_asin(D)),1
end

# `CheckLLCoincidence` — both trimmed curves coincide within tolerance when
# all four endpoints lie on the other line.
function _occ_checkllcoincidence(L1,L2,D1,D2,Tol)
    isf1=D1.hasf && _lin2d_dist(L2,D1.fpt)<Tol
    isl1=D1.hasl && _lin2d_dist(L2,D1.lpt)<Tol
    (isf1 && isl1) && return true
    isf2=D2.hasf && _lin2d_dist(L1,D2.fpt)<Tol
    isl2=D2.hasl && _lin2d_dist(L1,D2.lpt)<Tol
    return isf2 && isl2
end

# `computeIntPoint` — the single-point construction when a domain is
# narrower than the intersection tolerance. Returns the point or `nothing`.
function _occ_compute_int_point(curDom,othDom,curL,othL,cosT1T2,parCur,
                                parOther,resInf,resSup,num,curTrans)
    abs(resSup-parCur)>abs(resInf-parCur) && (resSup=resInf)
    aRes2=parOther+(resSup-parCur)*cosT1T2
    aFirst2=othDom.hasf ? othDom.fpar : -Inf
    aLast2=othDom.hasl ? othDom.lpar : Inf
    aTol21=othDom.hasf ? othDom.ftol : 0.0
    aTol22=othDom.hasl ? othDom.ltol : 0.0
    (aRes2<aFirst2-aTol21 || aRes2>aLast2+aTol22) && return nothing
    pos1,resSup=_occ_findpositionll(resSup,curDom)
    pos2,aRes2=_occ_findpositionll(aRes2,othDom)
    othTrans=curTrans===:out ? :in : curTrans===:in ? :out : :undecided
    if curTrans!==:undecided
        aT1=_itr(false,pos1,curTrans); aT2=_itr(false,pos2,othTrans)
    else
        opp=cosT1T2<0.0
        aT1=_itr(false,pos1,:unknown,opp)
        aT2=_itr(false,pos2,:unknown,opp)
    end
    aResU1=parCur; aResU2=parOther
    aFirst1=curDom.hasf ? curDom.fpar : -Inf
    aLast1=curDom.hasl ? curDom.lpar : Inf
    aTol11=curDom.hasf ? curDom.ftol : 0.0
    aTol12=curDom.hasl ? curDom.ltol : 0.0
    in1=parCur>=aFirst1 && parCur<=aLast1
    in2=parOther>=aFirst2 && parOther<=aLast2
    if !in1 || !in2
        if in1
            Pt1=_lin2d_value(othL,aRes2)
            aResU2=aRes2
            p1=_lin2d_param(curL,Pt1)
            aResU1=(p1>=aFirst1 && p1<=aLast1) ? p1 : resSup
        elseif in2
            aPt1=_lin2d_value(curL,resSup)
            aResU1=resSup
            p2=_lin2d_param(othL,aPt1)
            aResU2=(p2>=aFirst2 && p2<=aLast2) ? p2 : aRes2
        else
            (parCur<aFirst1-aTol11 || parCur>aLast1+aTol12 ||
             parOther<aFirst2-aTol21 || parOther>aLast2+aTol22) && return nothing
            aResU1=resSup; aResU2=aRes2
        end
    end
    q1=_lin2d_value(curL,aResU1); q2=_lin2d_value(othL,aResU2)
    aPres=((q1[1]+q2[1])*0.5,(q1[2]+q2[2])*0.5)
    return num==1 ? _ipt(aPres,aResU1,aResU2,aT1,aT2,false) :
                    _ipt(aPres,aResU2,aResU1,aT2,aT1,false)
end

# `LineCircleGeometricIntersection` — tolerance-band geometric intervals on
# the circle. Returns (CInt1, CInt2, nbsol) as PeriodicIntervals.
function _occ_lincirc_geometric(Line,Circle,Tol,TolTang)
    ax,r=_circ2d_axis(Circle)
    dO1O2=_lin2d_dist(Line,ax[1])
    R=r; RmTol=R-Tol
    binf1=0.0; bsup1=0.0; binf2=0.0; bsup2=0.0
    if dO1O2>R+Tol
        dO1O2>R+TolTang && return _pint(),_pint(),0
        nbsol=1
    else
        b2Sol=false
        if R>dO1O2+TolTang
            aX2=4.0*(R*R-dO1O2*dO1O2)
            aX2>Tol*Tol && (b2Sol=true)
        end
        if dO1O2>RmTol && !b2Sol
            dAlpha1=_gm_atan2(0.0,dO1O2)
            binf1=-dAlpha1; bsup1=dAlpha1
            nbsol=1
        else
            dx=dO1O2
            dy=R*R-dx*dx; dy=dy>=0.0 ? sqrt(dy) : 0.0
            dAlpha1=_gm_atan2(dy,dx)
            binf1=-dAlpha1; bsup2=dAlpha1
            dy=R*R-dx*dx; dy=dy>=0.0 ? sqrt(dy) : 0.0
            dAlpha1=_gm_atan2(dy,dx)
            binf2=dAlpha1; bsup1=-dAlpha1
            if dAlpha1*R<max(Tol,TolTang)
                bsup1=bsup2; nbsol=1
            else
                nbsol=2
            end
        end
    end
    dAngle1=_vec2d_angle(ax[2],Line[2])
    a,b,c=_lin2d_coefs(Line)
    d=a*ax[1][1]+b*ax[1][2]+c
    d>0.0 ? (dAngle1+=π/2) : (dAngle1-=π/2)
    if dAngle1<0.0
        dAngle1+=_PIPI
    elseif dAngle1>_PIPI
        dAngle1-=_PIPI
    end
    binf1+=dAngle1; bsup1+=dAngle1
    CInt1=_pint(binf1,bsup1)
    _pint_len(CInt1)>π && (CInt1=_pint_complement(CInt1))
    if nbsol==2
        binf2+=dAngle1; bsup2+=dAngle1
        CInt2=_pint(binf2,bsup2)
        _pint_len(CInt2)>π && (CInt2=_pint_complement(CInt2))
    else
        CInt2=_pint()
        if CInt1.bsup>_PIPI && CInt1.binf<_PIPI
            nbsol=2
            binf2=CInt1.binf; bsup2=_PIPI
            CInt1=_pint(0.0,CInt1.bsup-_PIPI)
            _pint_len(CInt1)>π && (CInt1=_pint_complement(CInt1))
            CInt2=_pint(binf2,bsup2)
            _pint_len(CInt2)>π && (CInt2=_pint_complement(CInt2))
        end
    end
    return CInt1,CInt2,nbsol
end

# `ProjectOnLAndIntersectWithLDomain` — map a circle interval to the line
# and intersect with the line domain; appends (circle, line) solution
# intervals to the output vectors and bumps the counter.
function _occ_project_on_l(Circle,Line,CDomainAndRes,LDomain,solC,solL,
                           RefLineDomain)
    CDomainAndRes.isnull && return
    ax,r=_circ2d_axis(Circle)
    Linf=_lin2d_param(Line,_circ2d_value(ax,r,CDomainAndRes.binf))
    Lsup=_lin2d_param(Line,_circ2d_value(ax,r,CDomainAndRes.bsup))
    LInter=_ivl(Linf,Lsup)
    LInterAndDomain=_ivl_intersect_bounded(LDomain,LInter)
    LInterAndDomain.isnull && return
    DomLinf=RefLineDomain.hasf ? RefLineDomain.fpar : -Inf
    DomLsup=RefLineDomain.hasl ? RefLineDomain.lpar : Inf
    Linf=LInterAndDomain.binf; Lsup=LInterAndDomain.bsup
    Linf<DomLinf && (Linf=DomLinf)
    Lsup<DomLinf && (Lsup=DomLinf)
    Linf>DomLsup && (Linf=DomLsup)
    Lsup>DomLsup && (Lsup=DomLsup)
    Cinf=CDomainAndRes.binf; Csup=CDomainAndRes.bsup
    if Cinf>=Csup
        Cinf=CDomainAndRes.binf; Csup=CDomainAndRes.bsup
    end
    cs=_pint(Cinf,Csup)
    _pint_len(cs)>π && (cs=_pint_complement(cs))
    push!(solC,cs)
    push!(solL,(binf=Linf,bsup=Lsup,hf=true,hl=true,isnull=false))
    return
end

# `IntCurve_IntConicConic::Perform(gp_Lin2d, D1, gp_Lin2d, D2, -, Tol)` —
# the lin×lin solver. `done` is always true; results are caller-ordered
# (the solver runs with the line first — `ReversedParameters` never set for
# a lin×lin GInter call).
function _occ_int_linlin(L1,Domain1,L2,Domain2,TolR)
    lpnt=Any[]; lseg=Any[]
    Tol=max(TolR,_OCC_PCONFUSION)
    U1,U2,aHalfSinL1L2,nbsol=_occ_linlin_geometric(L1,L2,Tol)
    Tan1=L1[2]; Tan2=L2[2]
    aCosT1T2=_dot2(Tan1,Tan2)
    isOpposite=aCosT1T2<0.0
    nbsol==1 && _occ_checkllcoincidence(L1,L2,Domain1,Domain2,Tol) && (nbsol=2)
    if nbsol==1
        d=0.5*Tol/aHalfSinL1L2
        U1inf=U1-d; U1sup=U1+d
        U1mU2=U1-U2; U1pU2=U1+U2
        if Domain1.hasf && _lin2d_dist(L2,Domain1.fpt)<Domain1.ftol
            U1inf>Domain1.fpar && (U1inf=Domain1.fpar)
            U1sup<Domain1.fpar && (U1sup=Domain1.fpar)
        end
        if Domain1.hasl && _lin2d_dist(L2,Domain1.lpt)<Domain1.ltol
            U1inf>Domain1.lpar && (U1inf=Domain1.lpar)
            U1sup<Domain1.lpar && (U1sup=Domain1.lpar)
        end
        if Domain2.hasf && _lin2d_dist(L1,Domain2.fpt)<Domain2.ftol
            p=_lin2d_param(L1,Domain2.fpt)
            U1inf>p && (U1inf=p)
            U1sup<p && (U1sup=p)
        end
        if Domain2.hasl && _lin2d_dist(L1,Domain2.lpt)<Domain2.ltol
            p=_lin2d_param(L1,Domain2.lpt)
            U1inf>p && (U1inf=p)
            U1sup<p && (U1sup=p)
        end
        Res1inf,Res1sup,Pos1a,Pos1b=
            _occ_domain_intersection(Domain1,U1inf,U1sup)
        if Res1sup-Res1inf<0.0
            # empty
        else
            ProdVectTan=_cross2(Tan1,Tan2)
            LongMiniSeg=Tol
            if (Res1sup-Res1inf)<=LongMiniSeg ||
               (Pos1a===Pos1b && Pos1a!==:middle)
                aCurTrans=ProdVectTan>=_OCC_TOLANG ? :out :
                          ProdVectTan<=-_OCC_TOLANG ? :in : :undecided
                np=_occ_compute_int_point(Domain1,Domain2,L1,L2,aCosT1T2,
                                          U1,U2,Res1inf,Res1sup,1,aCurTrans)
                np!==nothing && push!(lpnt,np)
            else
                if isOpposite
                    U2inf=U1pU2-Res1sup; U2sup=U1pU2-Res1inf
                else
                    U2inf=Res1inf-U1mU2; U2sup=Res1sup-U1mU2
                end
                Res2inf,Res2sup,Pos2a,Pos2b=
                    _occ_domain_intersection(Domain2,U2inf,U2sup)
                Res2sup_m_Res2inf=Res2sup-Res2inf
                if Res2sup_m_Res2inf<0.0
                    # no solution
                elseif Res2sup_m_Res2inf>LongMiniSeg ||
                       (Pos2a===Pos2b && Pos2a!==:middle)
                    if isOpposite
                        Res1inf=U1pU2-Res2sup
                        Res1sup=U1pU2-Res2inf
                        Res2inf,Res2sup=Res2sup,Res2inf
                        Pos2a,Pos2b=Pos2b,Pos2a
                    else
                        Res1inf=U1mU2+Res2inf
                        Res1sup=U1mU2+Res2sup
                    end
                    Pos1a,= _occ_findpositionll(Res1inf,Domain1)
                    Pos1b,= _occ_findpositionll(Res1sup,Domain1)
                    T1a=T2a=T1b=T2b=_itr()
                    if ProdVectTan>=_OCC_TOLANG
                        T1a=_itr(false,Pos1a,:out); T2a=_itr(false,Pos2a,:in)
                    elseif ProdVectTan<=-_OCC_TOLANG
                        T1a=_itr(false,Pos1a,:in); T2a=_itr(false,Pos2a,:out)
                    else
                        T1a=_itr(false,Pos1a,:unknown,isOpposite)
                        T2a=_itr(false,Pos2a,:unknown,isOpposite)
                    end
                    ResultIsAPoint=false
                    if (Res1sup-Res1inf)<=LongMiniSeg ||
                       abs(Res2sup-Res2inf)<=LongMiniSeg
                        ResultIsAPoint=true
                    else
                        if Pos1a===:head
                            if Pos1b!==:end && U1<Res1inf
                                ResultIsAPoint=true
                                U1=Res1inf; U2=Res2inf
                            end
                        end
                        if Pos1b===:end
                            if Pos1a!==:head && U1>Res1sup
                                ResultIsAPoint=true
                                U1=Res1sup; U2=Res2sup
                            end
                        end
                        if Pos2a===:head
                            if Pos2b!==:end && U2<Res2inf
                                ResultIsAPoint=true
                                U2=Res2inf; U1=Res1inf
                            end
                        elseif Pos2a===:end
                            if Pos2b!==:head && U2>Res2inf
                                ResultIsAPoint=true
                                U2=Res2inf; U1=Res1inf
                            end
                        end
                        if Pos2b===:head
                            if Pos2a!==:end && U2<Res2sup
                                ResultIsAPoint=true
                                U2=Res2sup; U1=Res1sup
                            end
                        elseif Pos2b===:end
                            if Pos2a!==:head && U2>Res2sup
                                ResultIsAPoint=true
                                U2=Res2sup; U1=Res1sup
                            end
                        end
                    end
                    if !ResultIsAPoint &&
                       (Pos1a!==:middle || Pos2a!==:middle)
                        if ProdVectTan>=_OCC_TOLANG
                            T1b=_itr(false,Pos1b,:out); T2b=_itr(false,Pos2b,:in)
                        elseif ProdVectTan<=-_OCC_TOLANG
                            T1b=_itr(false,Pos1b,:in); T2b=_itr(false,Pos2b,:out)
                        else
                            T1b=_itr(false,Pos1b,:unknown,isOpposite)
                            T2b=_itr(false,Pos2b,:unknown,isOpposite)
                        end
                        if Pos1a===:middle
                            t3=isOpposite ? (Pos2a===:head ? Res2sup : Res2inf) :
                                            (Pos2a===:head ? Res2inf : Res2sup)
                            Ptdebut=_lin2d_value(L2,t3)
                            Res1inf=_lin2d_param(L1,Ptdebut)
                        else
                            t4=Pos1a===:head ? Res1inf : Res1sup
                            Ptdebut=_lin2d_value(L1,t4)
                            Res2inf=_lin2d_param(L2,Ptdebut)
                        end
                        PtSeg1=(pt=Ptdebut,p1=Res1inf,p2=Res2inf,t1=T1a,t2=T2a)
                        if Pos1b!==:middle || Pos2b!==:middle
                            if Pos1b===:middle
                                Ptfin=_lin2d_value(L2,Res2sup)
                                Res1sup=_lin2d_param(L1,Ptfin)
                            else
                                Ptfin=_lin2d_value(L1,Res1sup)
                                Res2sup=_lin2d_param(L2,Ptfin)
                            end
                            PtSeg2=(pt=Ptfin,p1=Res1sup,p2=Res2sup,t1=T1b,t2=T2b)
                            push!(lseg,_iseg(PtSeg1,PtSeg2,isOpposite,false))
                        else
                            Pos1b,_=_occ_findpositionll(U1,Domain1)
                            Pos2b,_=_occ_findpositionll(U2,Domain2)
                            if ProdVectTan>=_OCC_TOLANG
                                T1b=_itr(false,Pos1b,:out)
                                T2b=_itr(false,Pos2b,:in)
                            elseif ProdVectTan<=-_OCC_TOLANG
                                T1b=_itr(false,Pos1b,:in)
                                T2b=_itr(false,Pos2b,:out)
                            else
                                T1b=_itr(false,Pos1b,:unknown,isOpposite)
                                T2b=_itr(false,Pos2b,:unknown,isOpposite)
                            end
                            PtSeg2=(pt=_lin2d_value(L2,U2),p1=U1,p2=U2,
                                    t1=T1b,t2=T2b)
                            if abs(Res1inf-U1)>LongMiniSeg &&
                               abs(Res2inf-U2)>LongMiniSeg
                                push!(lseg,_iseg(PtSeg1,PtSeg2,isOpposite,false))
                            else
                                push!(lpnt,_occ_segmenttopoint(
                                          PtSeg1,T1a,T2a,PtSeg2,T1b,T2b))
                            end
                        end
                    else
                        if Pos1b===:middle; Pos1b=Pos1a end
                        if Pos2b===:middle; Pos2b=Pos2a end
                        if ResultIsAPoint
                            if Pos1b!==:middle || Pos2b!==:middle
                                if Pos1b===:middle
                                    t2=isOpposite ?
                                        (Pos2b===:head ? Res2sup : Res2inf) :
                                        (Pos2b===:head ? Res2inf : Res2sup)
                                    Ptfin=_lin2d_value(L2,t2)
                                    Res1sup=_lin2d_param(L1,Ptfin)
                                    Pos1b,=_occ_findpositionll(Res1sup,Domain1)
                                else
                                    t1=Pos1b===:head ? Res1inf : Res1sup
                                    Ptfin=_lin2d_value(L1,t1)
                                    Res2sup=_lin2d_param(L2,Ptfin)
                                    Pos2b,=_occ_findpositionll(Res2sup,Domain2)
                                end
                                if ProdVectTan>=_OCC_TOLANG
                                    T1b=_itr(false,Pos1b,:out)
                                    T2b=_itr(false,Pos2b,:in)
                                elseif ProdVectTan<=-_OCC_TOLANG
                                    T1b=_itr(false,Pos1b,:in)
                                    T2b=_itr(false,Pos2b,:out)
                                else
                                    T1b=_itr(false,Pos1b,:unknown,isOpposite)
                                    T2b=_itr(false,Pos2b,:unknown,isOpposite)
                                end
                                PtSeg2=(pt=Ptfin,p1=Res1sup,p2=Res2sup,
                                        t1=T1b,t2=T2b)
                                push!(lpnt,PtSeg2)
                            else
                                Pos1b,=_occ_findpositionll(U1,Domain1)
                                Pos2b,=_occ_findpositionll(U2,Domain2)
                                if ProdVectTan>=_OCC_TOLANG
                                    T1b=_itr(false,Pos1b,:out)
                                    T2b=_itr(false,Pos2b,:in)
                                elseif ProdVectTan<=-_OCC_TOLANG
                                    T1b=_itr(false,Pos1b,:in)
                                    T2b=_itr(false,Pos2b,:out)
                                else
                                    T1b=_itr(false,Pos1b,:unknown,isOpposite)
                                    T2b=_itr(false,Pos2b,:unknown,isOpposite)
                                end
                                PtSeg1=(pt=_lin2d_value(L2,U2),p1=U1,p2=U2,
                                        t1=T1b,t2=T2b)
                                push!(lpnt,PtSeg1)
                            end
                        else
                            PtSeg1=(pt=_lin2d_value(L2,U2),p1=U1,p2=U2,
                                    t1=T1a,t2=T2a)
                            if Pos1b!==:middle || Pos2b!==:middle
                                if ProdVectTan>=_OCC_TOLANG
                                    T1b=_itr(false,Pos1b,:out)
                                    T2b=_itr(false,Pos2b,:in)
                                elseif ProdVectTan<=-_OCC_TOLANG
                                    T1b=_itr(false,Pos1b,:in)
                                    T2b=_itr(false,Pos2b,:out)
                                else
                                    T1b=_itr(false,Pos1b,:unknown,isOpposite)
                                    T2b=_itr(false,Pos2b,:unknown,isOpposite)
                                end
                                if Pos1b===:middle
                                    Ptfin=_lin2d_value(L2,Res2sup)
                                    Res1sup=_lin2d_param(L1,Ptfin)
                                else
                                    Ptfin=_lin2d_value(L1,Res1sup)
                                    Res2sup=_lin2d_param(L2,Ptfin)
                                end
                                PtSeg2=(pt=Ptfin,p1=Res1sup,p2=Res2sup,
                                        t1=T1b,t2=T2b)
                                if abs(U1-Res1sup)>LongMiniSeg ||
                                   abs(U2-Res2sup)>LongMiniSeg
                                    push!(lseg,_iseg(PtSeg1,PtSeg2,
                                                     isOpposite,false))
                                else
                                    push!(lpnt,_occ_segmenttopoint(
                                              PtSeg1,T1a,T2a,PtSeg2,T1b,T2b))
                                end
                            else
                                push!(lpnt,PtSeg1)
                            end
                        end
                    end
                else
                    aCurTrans=ProdVectTan>=_OCC_TOLANG ? :in :
                              ProdVectTan<=-_OCC_TOLANG ? :out : :undecided
                    np=_occ_compute_int_point(Domain2,Domain1,L2,L1,aCosT1T2,
                                              U2,U1,Res2inf,Res2sup,2,aCurTrans)
                    np!==nothing && push!(lpnt,np)
                end
            end
        end
    elseif nbsol==2
        ResHasFirstPoint=0; ResHasLastPoint=0
        ParamStart=0.0; ParamEnd=0.0; ParamStart2=0.0; ParamEnd2=0.0
        Org2SurL1=_lin2d_param(L1,L2[1])
        Domain1.hasf && (ResHasFirstPoint=1)
        Domain1.hasl && (ResHasLastPoint=1)
        if isOpposite
            Domain2.hasl && (ResHasFirstPoint+=2)
            Domain2.hasf && (ResHasLastPoint+=2)
        else
            Domain2.hasl && (ResHasLastPoint+=2)
            Domain2.hasf && (ResHasFirstPoint+=2)
        end
        if ResHasFirstPoint==0 && ResHasLastPoint==0
            push!(lseg,_iseg(isOpposite))
        else
            if ResHasFirstPoint==1
                ParamStart=Domain1.fpar
                ParamStart2=isOpposite ? Org2SurL1-ParamStart :
                                         ParamStart-Org2SurL1
            elseif ResHasFirstPoint==2
                if isOpposite
                    ParamStart2=Domain2.lpar
                    ParamStart=Org2SurL1-ParamStart2
                else
                    ParamStart2=Domain2.fpar
                    ParamStart=Org2SurL1+ParamStart2
                end
            elseif ResHasFirstPoint==3
                if isOpposite
                    ParamStart2=Domain2.lpar
                    ParamStart=Org2SurL1-ParamStart2
                    if ParamStart<Domain1.fpar
                        ParamStart=Domain1.fpar
                        ParamStart2=Org2SurL1-ParamStart
                    end
                else
                    ParamStart2=Domain2.fpar
                    ParamStart=Org2SurL1+ParamStart2
                    if ParamStart<Domain1.fpar
                        ParamStart=Domain1.fpar
                        ParamStart2=ParamStart-Org2SurL1
                    end
                end
            end
            if ResHasLastPoint==1
                ParamEnd=Domain1.lpar
                ParamEnd2=isOpposite ? Org2SurL1-ParamEnd : ParamEnd-Org2SurL1
            elseif ResHasLastPoint==2
                if isOpposite
                    ParamEnd2=Domain2.fpar
                    ParamEnd=Org2SurL1-ParamEnd2
                else
                    ParamEnd2=Domain2.lpar
                    ParamEnd=Org2SurL1+ParamEnd2
                end
            elseif ResHasLastPoint==3
                if isOpposite
                    ParamEnd2=Domain2.fpar
                    ParamEnd=Org2SurL1-ParamEnd2
                    if ParamEnd>Domain1.lpar
                        ParamEnd=Domain1.lpar
                        ParamEnd2=Org2SurL1-ParamEnd
                    end
                else
                    ParamEnd2=Domain2.lpar
                    ParamEnd=Org2SurL1+ParamEnd2
                    if ParamEnd>Domain1.lpar
                        ParamEnd=Domain1.lpar
                        ParamEnd2=ParamEnd-Org2SurL1
                    end
                end
            end
            if ResHasFirstPoint!=0
                if ResHasLastPoint!=0
                    if ParamEnd>=ParamStart-Tol
                        Pos1,_=_occ_findpositionll(ParamStart,Domain1)
                        Pos2,_=_occ_findpositionll(ParamStart2,Domain2)
                        Tinf=_itr(true,Pos1,:unknown,isOpposite)
                        Tsup=_itr(true,Pos2,:unknown,isOpposite)
                        P1=_ipt(_lin2d_value(L1,ParamStart),ParamStart,
                                ParamStart2,Tinf,Tsup,false)
                        if ParamEnd>ParamStart+Tol
                            Pos1,_=_occ_findpositionll(ParamEnd,Domain1)
                            Pos2,_=_occ_findpositionll(ParamEnd2,Domain2)
                            Tinf=_itr(true,Pos1,:unknown,isOpposite)
                            Tsup=_itr(true,Pos2,:unknown,isOpposite)
                            P2=_ipt(_lin2d_value(L1,ParamEnd),ParamEnd,
                                    ParamEnd2,Tinf,Tsup,false)
                            push!(lseg,_iseg(P1,P2,isOpposite,false))
                        else
                            push!(lpnt,P1)
                        end
                    end
                else
                    Pos1,_=_occ_findpositionll(ParamStart,Domain1)
                    Pos2,_=_occ_findpositionll(ParamStart2,Domain2)
                    Tinf=_itr(true,Pos1,:unknown,isOpposite)
                    Tsup=_itr(true,Pos2,:unknown,isOpposite)
                    P=_ipt(_lin2d_value(L1,ParamStart),ParamStart,ParamStart2,
                           Tinf,Tsup,false)
                    push!(lseg,_iseg_half(P,true,isOpposite,false))
                end
            else
                Pos1,_=_occ_findpositionll(ParamEnd,Domain1)
                Pos2,_=_occ_findpositionll(ParamEnd2,Domain2)
                Tinf=_itr(true,Pos1,:unknown,isOpposite)
                Tsup=_itr(true,Pos2,:unknown,isOpposite)
                P2=_ipt(_lin2d_value(L1,ParamEnd),ParamEnd,ParamEnd2,
                        Tinf,Tsup,false)
                push!(lseg,_iseg_half(P2,false,isOpposite,false))
            end
        end
    end
    return true,lpnt,lseg
end

# `IntCurve_IntConicConic::Perform(gp_Lin2d, LIG_Domain, gp_Circ2d,
# CIRC_Domain, TolConf, Tol)` — the lin×circ solver; `rev` plays the
# `ReversedParameters` role (result params/transitions are swapped into
# caller order at point/segment construction).
function _occ_int_lincirc(Line,LIG_Domain,Circle,CIRC_Domain,TolConf,Tol,rev)
    lpnt=Any[]; lseg=Any[]
    CInt1,CInt2,nbsol=_occ_lincirc_geometric(Line,Circle,TolConf,Tol)
    nbsol==0 && return true,lpnt,lseg
    if nbsol==2 && CInt2.bsup==CInt1.binf+_PIPI
        FirstBound=CIRC_Domain.fpar; LastBound=CIRC_Domain.lpar
        FirstTol=CIRC_Domain.ftol; LastTol=CIRC_Domain.ltol
        if CInt1.binf==0 && FirstBound-FirstTol>CInt1.bsup
            nbsol=1; CInt1=_pint(CInt2.binf,CInt2.bsup)
        elseif CInt2.bsup==_PIPI && LastBound+LastTol<CInt2.binf
            nbsol=1
        end
    end
    CDomain=_pint(CIRC_Domain)
    deltat=CDomain.bsup-CDomain.binf
    while CDomain.binf>=_PIPI; CDomain=(binf=CDomain.binf-_PIPI,bsup=CDomain.bsup,isnull=false) end
    while CDomain.binf<0.0; CDomain=(binf=CDomain.binf+_PIPI,bsup=CDomain.bsup,isnull=false) end
    CDomain=(binf=CDomain.binf,bsup=CDomain.binf+deltat,isnull=false)
    R=Circle[3]
    BinfModif=CDomain.binf-CIRC_Domain.ftol/R
    BsupModif=CDomain.bsup+CIRC_Domain.ltol/R
    deltat=BsupModif-BinfModif
    if deltat<=_PIPI
        CDomain=(binf=BinfModif,bsup=BsupModif,isnull=false)
    else
        t=(_PIPI-deltat)*0.5
        CDomain=(binf=BinfModif+t,bsup=BsupModif-t,isnull=false)
    end
    deltat=CDomain.bsup-CDomain.binf
    binf=CDomain.binf
    while binf>=_PIPI; binf-=_PIPI end
    while binf<0.0; binf+=_PIPI end
    CDomain=(binf=binf,bsup=binf+deltat,isnull=false)
    LDomain=_ivl(LIG_Domain)
    solC=Any[]; solL=Any[]
    CDomainAndRes,=_pint_first_intersection(CDomain,CInt1)
    _occ_project_on_l(Circle,Line,CDomainAndRes,LDomain,solC,solL,LIG_Domain)
    CDomainAndRes=_pint_second_intersection(CDomain,CInt1)
    _occ_project_on_l(Circle,Line,CDomainAndRes,LDomain,solC,solL,LIG_Domain)
    if nbsol==2
        CDomainAndRes,=_pint_first_intersection(CDomain,CInt2)
        _occ_project_on_l(Circle,Line,CDomainAndRes,LDomain,solC,solL,LIG_Domain)
        CDomainAndRes=_pint_second_intersection(CDomain,CInt2)
        _occ_project_on_l(Circle,Line,CDomainAndRes,LDomain,solC,solL,LIG_Domain)
    end
    NbSolTotal=length(solC)
    MaxTol=max(TolConf,Tol,1.0e-10)
    for i in 1:NbSolTotal
        if R*_pint_len(solC[i])<MaxTol && _ivl_len(solL[i])<MaxTol
            t=(solC[i].binf+solC[i].bsup)*0.5
            solC[i]=(binf=t,bsup=t,isnull=false)
            t=(solL[i].binf+solL[i].bsup)*0.5
            solL[i]=(binf=t,bsup=t,hf=true,hl=true,isnull=false)
        end
    end
    if NbSolTotal>0
        CircleAxis,=_circ2d_axis(Circle)
        LineAxis=Line
        _,Tan1,_=_circ2d_d2(CircleAxis,R,solC[1].binf)
        _,Tan2=_lin2d_d1(LineAxis,solL[1].binf)
        isOpposite=_dot2(Tan1,Tan2)<0.0
        for i in 1:NbSolTotal
            p1=solC[i].binf; p2=solC[i].bsup
            q1=CIRC_Domain.fpar; q2=CIRC_Domain.lpar
            if p1>q2
                while p1>q2; p1-=_PIPI; p2-=_PIPI end
            elseif p2<q1
                while p2<q1; p1+=_PIPI; p2+=_PIPI end
            end
            if p1<q1 && p2>q1; p1=q1 end
            if p1<q2 && p2>q2; p2=q2 end
            solC[i]=(binf=p1,bsup=p2,isnull=false)
            Linf=isOpposite ? solL[i].bsup : solL[i].binf
            Lsup=isOpposite ? solL[i].binf : solL[i].bsup
            if Linf>Lsup
                solC[i]=(binf=solC[i].bsup,bsup=solC[i].binf,isnull=false)
                Linf,Lsup=Lsup,Linf
            end
            P1a,Tan1,Norm1=_circ2d_d2(CircleAxis,R,solC[i].binf)
            P2a,Tan2=_lin2d_d1(LineAxis,Linf)
            Pos1a=_occ_determine_position(CIRC_Domain,P1a,solC[i].binf)
            Pos2a=_occ_determine_position(LIG_Domain,P2a,Linf)
            T1a,T2a=_occ_transition_lc(Pos1a,Tan1,Norm1,Pos2a,Tan2,(0.0,0.0))
            if Pos1a===:end
                Cinf=CIRC_Domain.lpar; P1a=CIRC_Domain.lpt
                Linf=_lin2d_param(Line,P1a)
                P1a,Tan1,Norm1=_circ2d_d2(CircleAxis,R,Cinf)
                P2a,Tan2=_lin2d_d1(LineAxis,Linf)
                Pos1a=_occ_determine_position(CIRC_Domain,P1a,Cinf)
                Pos2a=_occ_determine_position(LIG_Domain,P2a,Linf)
                T1a,T2a=_occ_transition_lc(Pos1a,Tan1,Norm1,Pos2a,Tan2,(0.0,0.0))
            elseif Pos1a===:head
                Cinf=CIRC_Domain.fpar; P1a=CIRC_Domain.fpt
                Linf=_lin2d_param(Line,P1a)
                P1a,Tan1,Norm1=_circ2d_d2(CircleAxis,R,Cinf)
                P2a,Tan2=_lin2d_d1(LineAxis,Linf)
                Pos1a=_occ_determine_position(CIRC_Domain,P1a,Cinf)
                Pos2a=_occ_determine_position(LIG_Domain,P2a,Linf)
                T1a,T2a=_occ_transition_lc(Pos1a,Tan1,Norm1,Pos2a,Tan2,(0.0,0.0))
            else
                Cinf=_occ_normalize_on_circle_domain(solC[i].binf,CIRC_Domain)
            end
            NewPoint1=_ipt(P1a,Linf,Cinf,T2a,T1a,rev)
            if _ivl_len(solL[i])+_pint_len(solC[i])>0.0
                P1b,Tan1,Norm1=_circ2d_d2(CircleAxis,R,solC[i].bsup)
                P2b,Tan2=_lin2d_d1(LineAxis,Lsup)
                Pos1b=_occ_determine_position(CIRC_Domain,P1b,solC[i].bsup)
                Pos2b=_occ_determine_position(LIG_Domain,P2b,Lsup)
                T1b,T2b=_occ_transition_lc(Pos1b,Tan1,Norm1,Pos2b,Tan2,(0.0,0.0))
                if Pos1b===:end
                    Csup=CIRC_Domain.lpar; P1b=CIRC_Domain.lpt
                    Lsup=_lin2d_param(Line,P1b)
                    P1b,Tan1,Norm1=_circ2d_d2(CircleAxis,R,Csup)
                    P2b,Tan2=_lin2d_d1(LineAxis,Lsup)
                    Pos1b=_occ_determine_position(CIRC_Domain,P1b,Csup)
                    Pos2b=_occ_determine_position(LIG_Domain,P2b,Lsup)
                    T1b,T2b=_occ_transition_lc(Pos1b,Tan1,Norm1,Pos2b,Tan2,(0.0,0.0))
                elseif Pos1b===:head
                    Csup=CIRC_Domain.fpar; P1b=CIRC_Domain.fpt
                    Lsup=_lin2d_param(Line,P1b)
                    P1b,Tan1,Norm1=_circ2d_d2(CircleAxis,R,Csup)
                    P2b,Tan2=_lin2d_d1(LineAxis,Lsup)
                    Pos1b=_occ_determine_position(CIRC_Domain,P1b,Csup)
                    Pos2b=_occ_determine_position(LIG_Domain,P2b,Lsup)
                    T1b,T2b=_occ_transition_lc(Pos1b,Tan1,Norm1,Pos2b,Tan2,(0.0,0.0))
                else
                    Csup=_occ_normalize_on_circle_domain(solC[i].bsup,CIRC_Domain)
                end
                NewPoint2=_ipt(P1b,Lsup,Csup,T2b,T1b,rev)
                if (abs(Csup-Cinf)*R>MaxTol && abs(Lsup-Linf)>MaxTol) ||
                   T1a.tr!==T2a.tr
                    push!(lseg,_iseg(NewPoint1,NewPoint2,isOpposite,rev))
                else
                    (Pos1a!==:middle || Pos2a!==:middle) &&
                        _ires_insert!(lpnt,NewPoint1)
                    (Pos1b!==:middle || Pos2b!==:middle) &&
                        _ires_insert!(lpnt,NewPoint2)
                end
            else
                _ires_insert!(lpnt,NewPoint1)
            end
        end
    end
    return true,lpnt,lseg
end

# `Geom2dInt_GInter` — the two-curve intersector restricted to lin2d/circ2d
# pcurves: lin×lin and lin×circ run through `IntConicConic`; circ×lin flips
# the argument order and sets `ReversedParameters`. Returns (done, points,
# segments); `nothing`-pcurve kinds fail loudly.
function _occ_ginter(L,DL,pc,DE,TolConf,Tol)
    kind=_occ_pcurve_kind(pc)
    if kind===:lin2d
        o,d=pc.lin2d
        return _occ_int_linlin((L[1],L[2]),DL,(o,d),DE,Tol)
    elseif kind===:circ2d
        o,x,r=pc.circ2d
        return _occ_int_lincirc((L[1],L[2]),DL,(o,x,r),DE,TolConf,Tol,false)
    end
    throw(ArgumentError("_occ_ginter: unsupported pcurve kind $kind"))
end

# ── Extrema_ExtPC2d on lin2d/circ2d (GExtPC analytic path) ───────────────────
# `ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)` — wraps U1 into the
# period frame anchored at UFirst, then U2 relative to U1.
function _occ_adjust_periodic(ufirst,ulast,preci)
    p=ulast-ufirst
    return function _adj(u1,u2)
        u1-=floor((u1-ufirst)/p)*p
        ulast-u1<preci && (u1-=p)
        u2-=floor((u2-u1)/p)*p
        u2-u1<preci && (u2+=p)
        return u1,u2
    end
end

# `Extrema_ExtPElC2d::Perform` for a line — returns [(param,point,sqdist)] or
# (done=false) marker. Tol is the `t3d` confusion passed by `GExtPC`.
function _occ_extpe2d_line(P,L,Tol,Uinf,Usup)
    o,d=L
    mydist=_dot2(_sub2(P,o),d)
    mydist>=Uinf-Tol && mydist<=Usup+Tol || return true,NTuple{3,Any}[]
    MyP=(fma(mydist,d[1],o[1]),fma(mydist,d[2],o[2]))
    return true,[(mydist,MyP,_sqdist2(P,MyP))]
end

# `Extrema_ExtPElC2d::Perform` for a circle — up to two extrema.
function _occ_extpe2d_circ(P,C,Tol,Uinf,Usup)
    ax,r=_circ2d_axis(C); OC=ax[1]
    _dist2(OC,P)<=_OCC_CONFUSION && return false,NTuple{3,Any}[]
    V=_dir2(_sub2(OC,P))
    P1=_add2(OC,_mul2(V,r)); U1=_circ2d_param(ax,P1)
    U2=U1+π; P2=_add2(OC,_mul2(V,-r))
    adj=_occ_adjust_periodic(Uinf,Uinf+2π,_OCC_PCONFUSION)
    myuinf,U1=adj(Uinf,U1)
    _,U2=adj(myuinf,U2)
    if U1-2π-Uinf<Tol && U1-2π-Uinf>-Tol
        U1=Uinf
        P1=_circ2d_value(ax,r,U1)
    end
    if U2-2π-Uinf<Tol && U2-2π-Uinf>-Tol
        U2=Uinf
        P2=_circ2d_value(ax,r,U2)
    end
    res=Tuple{Float64,NTuple{2,Float64},Float64}[]
    if Uinf-U1<Tol && U1-Usup<Tol
        push!(res,(U1,P1,_sqdist2(P,P1)))
    end
    if Uinf-U2<Tol && U2-Usup<Tol
        push!(res,(U2,P2,_sqdist2(P,P2)))
    end
    return true,res
end

# `Extrema_ExtPC2d`/`GExtPC2d` over a trimmed lin2d/circ2d pcurve — the
# elementary-curve path of `GExtPC::Perform`: `ExtPElC2d` candidates filtered
# by `[uinf−tolu, usup+tolu]` after `InPeriod` mapping for periodic curves.
function _occ_extpc2d(P,pc,uinf,usup)
    kind=_occ_pcurve_kind(pc)
    tolu=_occ_pcurve_resolution(pc,_OCC_CONFUSION)
    if kind===:lin2d
        done,res=_occ_extpe2d_line(P,pc.lin2d,_OCC_CONFUSION,uinf,usup)
    elseif kind===:circ2d
        done,res=_occ_extpe2d_circ(P,pc.circ2d,_OCC_CONFUSION,uinf,usup)
    else
        throw(ArgumentError("_occ_extpc2d: unsupported pcurve kind $kind"))
    end
    done || return false,Tuple{Float64,NTuple{2,Float64},Float64}[]
    per=_occ_pcurve_isperiodic(pc)
    out=Tuple{Float64,NTuple{2,Float64},Float64}[]
    for (U,pv,sd) in res
        if per
            U=_occ_inperiod(U,uinf,uinf+_PIPI)
        end
        if U>=uinf-tolu && U<=usup+tolu
            push!(out,(U,pv,sd))
        end
    end
    return true,out
end

# `ElCLib::InPeriod`.
function _occ_inperiod(u,ufirst,ulast)
    (_occ_isinf(u) || _occ_isinf(ufirst) || _occ_isinf(ulast)) && return u
    aperiod=ulast-ufirst
    aperiod<_occ_epsilon(ulast) && return u
    return max(ufirst,u+aperiod*ceil((ufirst-u)/aperiod))
end

# `BRepClass_Intersector::Perform` — CheckOn (ExtPC2d of the ray origin onto
# the trimmed pcurve) then `Geom2dInt_GInter` against the probing segment.
# `CheckSkip` is ported but unreachable here: it fires only for vertex
# tolerances above `maxTol` (0.1), while materialized vertices carry
# `Precision::Confusion()`. Returns (done, lpnt, lseg).
function _occ_intersector_perform(m::GeoModel,g,L,P,Tol,etag)
    ent=_occ_face_pcurve(m,g,etag)
    ent===nothing && return false,Any[],Any[]
    pc,deb,fin=ent
    aTolZ=Tol
    aDebTol=deb; aFinTol=fin
    if aTolZ>_OCC_CONFUSION
        aDebTol=deb-aTolZ; aFinTol=fin+aTolZ
    end
    # CheckOn — `Extrema_ExtPC2d(L.Location(), adaptor(pcurve,[aDebTol,aFinTol]))`
    done,exts=_occ_extpc2d(L[1],pc,aDebTol,aFinTol)
    aMinDist=Inf; aMinInd=0; aPar=0.0; aPntExact=(0.0,0.0)
    if done
        for (i,(U,pv,sd)) in enumerate(exts)
            if sd<aMinDist
                aMinDist=sd; aMinInd=i; aPar=U; aPntExact=pv
            end
        end
    end
    aMinInd!=0 && (aMinDist=sqrt(aMinDist))
    if aMinDist<=aTolZ
        aTolZ=_occ_refine_tolerance(g,pc,aPar,aTolZ)
        if aMinDist<=aTolZ
            aPosOnCurve=:middle
            if abs(aPar-deb)<=_OCC_CONFUSION || aPar<deb
                aPosOnCurve=:head
            elseif abs(aPar-fin)<=_OCC_CONFUSION || aPar>fin
                aPosOnCurve=:end
            end
            p=_ipt(aPntExact,0.0,aPar,_itr(:head),_itr(aPosOnCurve),false)
            return true,Any[p],Any[]
        end
    end
    pdeb=_occ_pcurve_value(pc,deb); pfin=_occ_pcurve_value(pc,fin)
    toldeb=1e-5; tolfin=1e-5
    DL=P!=floatmax(Float64) ?
        _idom(L[1],0.0,_OCC_PCONFUSION,_lin2d_value(L,P),P,_OCC_PCONFUSION) :
        _idom(L[1],0.0,_OCC_PCONFUSION,true)
    DE=_idom(pdeb,deb,toldeb,pfin,fin,tolfin)
    _occ_pcurve_isperiodic(pc) &&
        (DE=_idom_per(DE,deb,fin))
    done,lpnt,lseg=_occ_ginter(L,DL,pc,DE,_OCC_PCONFUSION,_OCC_PINTERSECTION)
    # `CheckSkip` — vertex tolerances are `Precision::Confusion()` for all
    # materialized edges, never above `myMaxTolerance` (0.1): the C++ early
    # return `!(BRep_Tool::Tolerance(aVl) > theMaxTol)` fires unconditionally.
    return done,lpnt,lseg
end

# `RefineTolerance` — only tightens on cylinder faces:
# `aTolX = |URes·dY + VRes·dX|` with `URes/VRes` the surface resolutions.
function _occ_refine_tolerance(g,pc,aT,aTolZ)
    g.occ===:cylinder || return aTolZ
    aURes=_occ_surface_uresolution(g,aTolZ)
    aVRes=_occ_surface_vresolution(g,aTolZ)
    _,aV2D=_occ_pcurve_d1(pc,aT)
    aD2D=_dir2(aV2D)
    aTolX=aURes*aD2D[2]+aVRes*aD2D[1]
    aTolX<0.0 && (aTolX=-aTolX)
    aTolX<_OCC_CONFUSION && (aTolX=_OCC_CONFUSION)
    aTolX<aTolZ && (aTolZ=aTolX)
    return aTolZ
end

# `GeomAdaptor_Surface::UResolution` (unrestricted adaptor) — periodic
# directions resolve `2·asin(R3d/(2R))` (or 2π when the radius is too small);
# a bounded cone resolves R3d / max end radius, unbounded falls back to
# `Precision::Parametric(R3d)` = R3d/100.
function _occ_surface_uresolution(g,tol)
    if g.occ===:cylinder || g.occ===:sphere
        Res=g.radius>_OCC_CONFUSION ? tol/(2.0*g.radius) : 0.0
        return Res<=1.0 ? 2.0*_gm_asin(Res) : _OCC_TWO_PI
    elseif g.occ===:torus
        Res=g.r1+g.r2>_OCC_CONFUSION ? tol/(2.0*(g.r1+g.r2)) : 0.0
        return Res<=1.0 ? 2.0*_gm_asin(Res) : _OCC_TWO_PI
    elseif g.occ===:cone
        # `GeomAbs_Cone`: unbounded V (`myVLast−myVFirst > 1e10`) falls back to
        # `Precision::Parametric`; otherwise R3d / max end radius — a plain
        # ratio, not the 2·asin tail. Materialized cones span V ∈ [0, len].
        vlen=sqrt(g.height*g.height+(g.r1-g.r2)*(g.r1-g.r2))
        vlen>1e10 && return tol*0.01
        R=max(g.r1,g.r2)
        return R>_OCC_CONFUSION ? tol/R : 0.0
    end
    return tol
end
# `GeomAdaptor_Surface::VResolution` — cylinder/cone/plane return R3d;
# sphere uses its radius and torus its minor radius.
function _occ_surface_vresolution(g,tol)
    if g.occ===:sphere
        Res=g.radius>_OCC_CONFUSION ? tol/(2.0*g.radius) : 0.0
        return Res<=1.0 ? 2.0*_gm_asin(Res) : _OCC_TWO_PI
    elseif g.occ===:torus
        Res=g.r2>_OCC_CONFUSION ? tol/(2.0*g.r2) : 0.0
        return Res<=1.0 ? 2.0*_gm_asin(Res) : _OCC_TWO_PI
    end
    return tol
end

# `Geom2dLProp_CLProps2d` over a pcurve — `IsTangentDefined`/`Tangent`/
# `Curvature`/`Normal` (orders 1..3, `myLinTol` = PConfusion). Returns
# (tangdef, tang, curv, norm) with norm=nothing when curvature is 0/∞.
function _occ_clprops2d(pc,U)
    mylinTol=_OCC_PCONFUSION; Tol=mylinTol*mylinTol
    _,d1=_occ_pcurve_d1(pc,U)
    _,_,d2=_occ_pcurve_d2(pc,U)
    order=0
    if _sqmag2(d1)>Tol
        order=1
    elseif _sqmag2(d2)>Tol
        order=2
    else
        order=0
    end
    order==0 && return false,(0.0,0.0),0.0,nothing
    tang=_dir2(d1)
    DD1=_sqmag2(d1); DD2=_sqmag2(d2)
    curv=0.0
    if DD2>Tol
        N=_cross2(d1,d2)^2
        t=N/DD1/DD2
        curv=t<=Tol ? 0.0 : sqrt(N)/DD1/sqrt(DD1)
    end
    norm=nothing
    if curv>_OCC_PCONFUSION && !_occ_isinf(curv)
        Norm=_sub2(_mul2(d2,_dot2(d1,d1)),_mul2(d1,_dot2(d1,d2)))
        norm=_dir2(Norm)
    end
    return true,tang,curv,norm
end

# `GetTangentAsChord` — fallback tangent when `IsTangentDefined` fails.
function _occ_tangent_as_chord(pc,theParam,theFirst,theLast)
    offset=0.1*(theLast-theFirst)
    if theLast-theParam<_OCC_PCONFUSION
        offset*=-1
    elseif theParam+offset>theLast
        offset=0.5*(theLast-theParam)
    end
    p1=_occ_pcurve_value(pc,theParam)
    p2=_occ_pcurve_value(pc,theParam+offset)
    chord=_sub2(p2,p1)
    offset<0.0 && (chord=(-chord[1],-chord[2]))
    _sqmag2(chord)>_OCC_SQPARTOL || return (0.0,0.0)
    return _dir2(chord)
end

# `BRepClass_Intersector::LocalGeometry` — (Tang, Norm, C).
function _occ_pcurve_localgeometry(m::GeoModel,g,etag,U)
    ent=_occ_face_pcurve(m,g,etag)
    ent===nothing && throw(ArgumentError(
        "_occ_pcurve_localgeometry: edge $etag carries no pcurve on this face"))
    pc,fpar,lpar=ent
    tangdef,Tang,C,n=_occ_clprops2d(pc,U)
    !tangdef && (Tang=_occ_tangent_as_chord(pc,U,fpar,lpar))
    Norm=n!==nothing ? n : (Tang[2],-Tang[1])
    return Tang,Norm,C
end

# ── TopTrans_CurveTransition ─────────────────────────────────────────────────
# 2-D specialisation: `myTgt` is the probing-line direction (myCurv always 0).
mutable struct _OccTopTrans
    myTgt::NTuple{2,Float64}
    TgtFirst::NTuple{2,Float64}; NormFirst::NTuple{2,Float64}
    CurvFirst::Float64; TranFirst::Symbol
    TgtLast::NTuple{2,Float64};  NormLast::NTuple{2,Float64}
    CurvLast::Float64;  TranLast::Symbol
    Init::Bool
    _OccTopTrans()=new((0.0,0.0),(0.0,0.0),(0.0,0.0),0.0,:forward,
                       (0.0,0.0),(0.0,0.0),0.0,:forward,true)
end
_occ_toptrans_reset!(tr::_OccTopTrans,T)=
    (tr.myTgt=T; tr.Init=true; return tr)

# `TopTrans::Reverse` is not needed here — `S` arrives already adjusted.
function _occ_toptrans_compare!(tr::_OccTopTrans,Tole,T,N,C,St,Or)
    S=St; O=Or
    if S===:internal
        S=_dot2(T,tr.myTgt)<0 ? _occ_reverse_orient(O) : O
    end
    if tr.Init
        tr.Init=false
        tr.TgtFirst=T; tr.NormFirst=N; tr.CurvFirst=C; tr.TranFirst=S
        tr.TgtLast=T;  tr.NormLast=N;  tr.CurvLast=C;  tr.TranLast=S
        if O===:reversed
            tr.TgtFirst=(-tr.TgtFirst[1],-tr.TgtFirst[2])
            tr.TgtLast=(-tr.TgtLast[1],-tr.TgtLast[2])
        elseif O===:internal
            if _dot2(tr.myTgt,T)>0
                tr.TgtFirst=(-tr.TgtFirst[1],-tr.TgtFirst[2])
            else
                tr.TgtLast=(-tr.TgtLast[1],-tr.TgtLast[2])
            end
        end
        return
    end
    FirstSet=false
    cosAngWithT=_dot2(tr.myTgt,T)
    O===:reversed && (cosAngWithT=-cosAngWithT)
    O===:internal && cosAngWithT>0 && (cosAngWithT=-cosAngWithT)
    cosAngWith1=_dot2(tr.myTgt,tr.TgtFirst)
    cmp=_occ_angcmp(cosAngWithT,cosAngWith1,Tole)
    if cmp<0
        FirstSet=true
        tr.TgtFirst=_occ_ori_tangent(T,O,tr.myTgt,:first)
        tr.NormFirst=N; tr.CurvFirst=C; tr.TranFirst=S
    elseif cmp==0
        if _occ_toptrans_isbefore(tr,Tole,cosAngWithT,N,C,tr.NormFirst,tr.CurvFirst)
            FirstSet=true
            tr.TgtFirst=_occ_ori_tangent(T,O,tr.myTgt,:first)
            tr.NormFirst=N; tr.CurvFirst=C; tr.TranFirst=S
        end
    end
    if !FirstSet || O===:internal
        O===:internal && (cosAngWithT=-cosAngWithT)
        cosAngWith2=_dot2(tr.myTgt,tr.TgtLast)
        cmp=_occ_angcmp(cosAngWithT,cosAngWith2,Tole)
        if cmp>0
            tr.TgtLast=_occ_ori_tangent(T,O,tr.myTgt,:last)
            tr.NormLast=N; tr.CurvLast=C; tr.TranLast=S
        elseif cmp==0
            if _occ_toptrans_isbefore(tr,Tole,cosAngWithT,tr.NormLast,tr.CurvLast,N,C)
                tr.TgtLast=_occ_ori_tangent(T,O,tr.myTgt,:last)
                tr.NormLast=N; tr.CurvLast=C; tr.TranLast=S
            end
        end
    end
    return
end
function _occ_angcmp(a1,a2,tole)
    a1-a2>tole && return 1
    a2-a1>tole && return -1
    return 0
end
_occ_reverse_orient(o)=o===:forward ? :reversed :
    o===:reversed ? :forward : o===:internal ? :external : :internal
# Orientation-aware tangent: REVERSED flips; INTERNAL flips in the direction
# of `myTgt` for first, against for last.
function _occ_ori_tangent(T,O,myTgt,which)
    if O===:reversed
        return (-T[1],-T[2])
    elseif O===:internal
        d=_dot2(myTgt,T)
        flip=which===:first ? d>0 : d<0
        flip && return (-T[1],-T[2])
    end
    return T
end
# `TopTrans_CurveTransition::IsBefore` — `myCurv` is 0 (straight reference).
function _occ_toptrans_isbefore(tr,Tole,CosAngl,N1,C1,N2,C2)
    TN1=_dot2(tr.myTgt,N1); TN2=_dot2(tr.myTgt,N2)
    OneBefore=false
    if abs(TN1)<=Tole || abs(TN2)<=Tole
        C1<C2 && (OneBefore=true)
        CosAngl>0 && (OneBefore=!OneBefore)
    elseif TN1<0
        if TN2>0
            OneBefore=true
        elseif C1>C2
            OneBefore=true
        end
    elseif TN1>0
        TN2>0 && C1<C2 && (OneBefore=true)
    end
    return OneBefore
end
_occ_toptrans_statebefore(tr)=tr.Init ? :unknown :
    (tr.TranFirst===:forward || tr.TranFirst===:external ? :out : :in)

# ── TopClass_Classifier2d ────────────────────────────────────────────────────
mutable struct _OccClassifier
    myLin::Tuple{NTuple{2,Float64},NTuple{2,Float64}}
    myParam::Float64; myTolerance::Float64
    myState::Symbol; myFirstCompare::Bool; myFirstTrans::Bool
    myClosest::Int; myIsHeadOrEnd::Bool
    lpnt::Vector{Any}; lseg::Vector{Any}
    myTrans::_OccTopTrans
    _OccClassifier()=new(((0.0,0.0),(1.0,0.0)),0.0,0.0,:unknown,
                         true,true,0,false,Any[],Any[],_OccTopTrans())
end
function _occ_classifier_reset!(c::_OccClassifier,L,P,Tol)
    c.myLin=L; c.myParam=P; c.myTolerance=Tol
    c.myState=:unknown; c.myFirstCompare=true; c.myFirstTrans=true
    c.myClosest=0; c.myIsHeadOrEnd=false
    return c
end

# `TopClass_Classifier2d::Compare` — run the intersector for edge `etag`
# (signed; the sign carries the FORWARD/REVERSED orientation).
function _occ_classifier_compare!(m::GeoModel,g,c::_OccClassifier,etag,Or)
    c.myClosest=0
    done,lpnt,lseg=_occ_intersector_perform(m,g,c.myLin,c.myParam,c.myTolerance,etag)
    !done && return
    isempty(lpnt) && isempty(lseg) && return
    dMin=floatmax(Float64); PClosest=nothing
    for (i,PInter) in enumerate(lpnt)
        if PInter.t1.pos===:head
            c.myClosest=i; c.myState=:on; return
        end
        if PInter.p1<dMin
            c.myClosest=i; PClosest=PInter; dMin=PInter.p1
        end
    end
    for (i,Seg) in enumerate(lseg)
        PInter=Seg.p1
        if PInter.t1.pos===:head
            c.myClosest=length(lpnt)+i+i-1; c.myState=:on; return
        end
        if PInter.p1<dMin
            c.myClosest=length(lpnt)+i+i-1; PClosest=PInter; dMin=PInter.p1
        end
    end
    c.myClosest==0 && return
    if Or===:internal
        c.myState=:in; return
    elseif Or===:external
        c.myState=:out; return
    end
    if !c.myFirstCompare
        dMin>c.myParam && return
    end
    c.myFirstCompare=false
    c.myParam>dMin && (c.myFirstTrans=true)
    c.myParam=dMin
    T2=PClosest.t2
    c.myIsHeadOrEnd=(T2.pos===:head || T2.pos===:end)
    SegTrans=:forward
    T1=PClosest.t1
    if T1.tr===:in
        SegTrans=Or===:reversed ? :reversed : :forward
    elseif T1.tr===:out
        SegTrans=Or===:reversed ? :forward : :reversed
    elseif T1.tr===:touch
        if T1.situ===:inside
            SegTrans=Or===:reversed ? :external : :internal
        elseif T1.situ===:outside
            SegTrans=Or===:reversed ? :internal : :external
        else
            return
        end
    else
        return
    end
    if !c.myIsHeadOrEnd
        if SegTrans===:forward || SegTrans===:external
            c.myState=:out
        else
            c.myState=:in
        end
    else
        Tang,Norm,Curv=_occ_pcurve_localgeometry(m,g,etag,PClosest.p2)
        if c.myFirstTrans
            _occ_toptrans_reset!(c.myTrans,c.myLin[2])
            c.myFirstTrans=false
        end
        Ort=T2.pos===:head ? :forward : :reversed
        _occ_toptrans_compare!(c.myTrans,eps(Float64),Tang,Norm,Curv,SegTrans,Ort)
        c.myState=_occ_toptrans_statebefore(c.myTrans)
    end
    return
end

# ── BRepTools::UVBounds ──────────────────────────────────────────────────────
# `BndLib_Add2dCurve::Add` for lin2d/circ2d: line box = the two endpoints;
# circle box = endpoints plus the four coordinate extrema that fall inside
# `[t1z, t1z + (t2−t1)]` with `t1z` the period-adjusted first parameter.
function _occ_bndlib_pcurve_box(pc,t1,t2)
    pts=NTuple{2,Float64}[]
    push!(pts,_occ_pcurve_value(pc,t1),_occ_pcurve_value(pc,t2))
    kind=_occ_pcurve_kind(pc)
    if kind===:circ2d
        o,x,r=pc.circ2d
        aCosBt,aSinBt=x
        aCosGm,aSinGm=-x[2],x[1]
        pT=Float64[]
        for i in 0:1
            aLx=i==0 ? 0.0 : 1.0; aLy=i==0 ? 1.0 : 0.0
            aBx=aLx*(-r*aSinBt)-aLy*(-r*aCosBt)
            aBy=aLx*(r*aSinGm)-aLy*(r*aCosGm)
            aB=sqrt(aBx*aBx+aBy*aBy)
            aCosFi=aBx/aB; aSinFi=aBy/aB
            aFi=_gm_acos(aCosFi)
            aSinFi<0.0 && (aFi=_PIPI-aFi)
            push!(pT,_occ_adjust_to_period(_PIPI-aFi,_PIPI))
            push!(pT,_occ_adjust_to_period(π-aFi,_PIPI))
        end
        aEps=1e-14; dT=t2-t1
        aT1z=_occ_adjust_to_period(t1,_PIPI)
        abs(aT1z)<aEps && (aT1z=0.0)
        aT2z=aT1z+dT
        abs(aT2z-_PIPI)<aEps && (aT2z=_PIPI)
        for aT in pT
            aT=(aT<aT1z ? aT+_PIPI : aT)
            aT<=aT2z && push!(pts,_occ_pcurve_value(pc,aT))
        end
    elseif kind!==:lin2d
        throw(ArgumentError(
            "_occ_bndlib_pcurve_box: unsupported pcurve kind $kind"))
    end
    umin=vmin=Inf; umax=vmax=-Inf
    for p in pts
        umin=min(umin,p[1]); vmin=min(vmin,p[2])
        umax=max(umax,p[1]); vmax=max(vmax,p[2])
    end
    return umin,vmin,umax,vmax
end

# `BndLib_Box2dCurve::AdjustToPeriod`.
function _occ_adjust_to_period(aT,aPeriod)
    if aT<0.0
        k=1+floor(Int,-aT/aPeriod)
        return aT+k*aPeriod
    elseif aT>aPeriod
        k=floor(Int,aT/aPeriod)
        return aT-k*aPeriod
    end
    return aT
end

# `BRepTools::UVBounds(F)` — the union of per-edge pcurve boxes, each clamped
# into the non-periodic surface bounds where they cross (`AddUVBounds(F,E)`);
# an empty box yields (0,0,0,0).
function _occ_uvbounds(m::GeoModel,stag::Int,g)
    umin=vmin=Inf; umax=vmax=-Inf; found=false
    su1,su2,sv1,sv2=_occ_surface_natural_bounds(g)
    for loop in get(m.surfaces,stag,Int[])
        for etag in m.loops[loop]
            ent=_occ_face_pcurve(m,g,etag)
            ent===nothing && continue
            pc,t1,t2=ent
            found=true
            aXmin,aYmin,aXmax,aYmax=_occ_bndlib_pcurve_box(pc,t1,t2)
            if !_occ_surface_isuperiodic(g)
                (aXmin<su1 && su1<aXmax) && (aXmin=su1)
                (aXmin<su2 && su2<aXmax) && (aXmax=su2)
            end
            if !_occ_surface_isvperiodic(g)
                (aYmin<sv1 && sv1<aYmax) && (aYmin=sv1)
                (aYmin<sv2 && sv2<aYmax) && (aYmax=sv2)
            end
            umin=min(umin,aXmin); vmin=min(vmin,aYmin)
            umax=max(umax,aXmax); vmax=max(vmax,aYmax)
        end
    end
    if !found
        return 0.0,0.0,0.0,0.0
    end
    return umin,umax,vmin,vmax
end

# Natural `Geom_Surface::Bounds` (U1,U2,V1,V2) for materialized surfaces —
# `±Precision::Infinite()` where OCCT leaves the direction unbounded.
function _occ_surface_natural_bounds(g)
    g.occ===:cylinder && return 0.0,_OCC_TWO_PI,-_OCC_INFINITE,_OCC_INFINITE
    g.occ===:sphere && return 0.0,_OCC_TWO_PI,-π/2,π/2
    g.occ===:cone && return 0.0,_OCC_TWO_PI,-_OCC_INFINITE,_OCC_INFINITE
    g.occ===:torus && return 0.0,_OCC_TWO_PI,0.0,_OCC_TWO_PI
    g.occ===:plane && return -_OCC_INFINITE,_OCC_INFINITE,
        -_OCC_INFINITE,_OCC_INFINITE
    throw(ArgumentError("_occ_surface_natural_bounds: $(g.occ)"))
end
_occ_surface_isuperiodic(g)=g.occ===:cylinder || g.occ===:sphere ||
    g.occ===:cone || g.occ===:torus
_occ_surface_isvperiodic(g)=g.occ===:torus

# ── BRepClass_FaceExplorer ───────────────────────────────────────────────────
# Flat edge iteration (all wires in order, signed tags); the probing
# parameter walks `Probing_Start → Probing_End` by `Probing_Step` on each edge.
mutable struct _OccFaceExplorer
    m::GeoModel; stag::Int; g
    umin::Float64; umax::Float64; vmin::Float64; vmax::Float64
    boundsdone::Bool
    curEdgeInd::Int; curEdgePar::Float64
    edges::Vector{Int}
    _OccFaceExplorer(m,stag,g)=new(m,stag,g,_OCC_INFINITE,-_OCC_INFINITE,
                                   _OCC_INFINITE,-_OCC_INFINITE,false,
                                   1,_OCC_PROBING[1],Int[])
end
# `TopExp_Explorer(myFace, TopAbs_EDGE)` — every edge occurrence of every
# wire, in wire order, with the sign carrying FORWARD/REVERSED.
function _occ_fex_edges!(e::_OccFaceExplorer)
    isempty(e.edges) || return e.edges
    for loop in get(e.m.surfaces,e.stag,Int[])
        for etag in e.m.loops[loop]
            push!(e.edges,etag)
        end
    end
    return e.edges
end
function _occ_fex_bounds!(e::_OccFaceExplorer)
    e.boundsdone && return
    e.boundsdone=true
    u1,u2,v1,v2=_occ_surface_natural_bounds(e.g)
    if _occ_isinf(u1) || _occ_isinf(u2) || _occ_isinf(v1) || _occ_isinf(v2)
        u1,u2,v1,v2=_occ_uvbounds(e.m,e.stag,e.g)
    end
    e.umin,e.umax,e.vmin,e.vmax=u1,u2,v1,v2
    return
end
# `BRepClass_FaceExplorer::CheckPoint` — mutates the point when it lies too
# far from the face box; returns false while adjusting.
function _occ_fex_checkpoint!(e::_OccFaceExplorer,P)
    !e.boundsdone && _occ_fex_bounds!(e)
    (_occ_isinf(e.umin) || _occ_isinf(e.umax) ||
     _occ_isinf(e.vmin) || _occ_isinf(e.vmax)) && return true,P
    aCenter=((e.umin+e.umax)/2,(e.vmin+e.vmax)/2)
    aDistance=_dist2(aCenter,P)
    if _occ_isinf(aDistance)
        return false,(e.umin-(e.umax-e.umin),e.vmin-(e.vmax-e.vmin))
    end
    anEpsilon=_occ_epsilon(aDistance)
    if anEpsilon>max(e.umax-e.umin,e.vmax-e.vmin)
        aLinDir=_dir2(_sub2(P,aCenter))
        return false,_add2(aCenter,_mul2(aLinDir,2.0*anEpsilon))
    end
    return true,P
end
# `BRepClass_FaceExplorer::Segment/OtherSegment` — probing-ray construction.
# Returns (valid, L, Par).
function _occ_fex_other_segment!(e::_OccFaceExplorer,P)
    edges=_occ_fex_edges!(e)
    aTolParConf2=_OCC_SQPARTOL
    n=length(edges)
    while e.curEdgeInd<=n
        i=e.curEdgeInd
        etag=edges[i]
        ent=_occ_face_pcurve(e.m,e.g,etag)
        if ent===nothing
            e.curEdgeInd+=1; e.curEdgePar=_OCC_PROBING[1]
            continue
        end
        pc,aFPar,aLPar=ent
        if _occ_isinfneg(aFPar)
            if _occ_isinfpos(aLPar)
                aFPar=-1.0; aLPar=1.0
            else
                aFPar=aLPar-1.0
            end
        elseif _occ_isinfpos(aLPar)
            aLPar=aFPar+1.0
        end
        while e.curEdgePar<_OCC_PROBING[2]
            aParamIn=e.curEdgePar*aFPar+(1.0-e.curEdgePar)*aLPar
            aPOnC,aTanVec=_occ_pcurve_d1(pc,aParamIn)
            Par=_sqdist2(aPOnC,P)
            if Par>aTolParConf2
                aLinDir=_dir2(_sub2(aPOnC,P))
                aTanMod=_sqmag2(aTanVec)
                aTanMod<aTolParConf2 &&
                    (e.curEdgePar+=_OCC_PROBING[3]; continue)
                aTanVec=_mul2(aTanVec,1.0/sqrt(aTanMod))
                aSinA=_cross2(aTanVec,aLinDir)
                SmallAngle=0.001
                isSmallAngle=false
                if abs(aSinA)<SmallAngle
                    isSmallAngle=true
                    if e.curEdgePar+_OCC_PROBING[3]<_OCC_PROBING[2]
                        e.curEdgePar+=_OCC_PROBING[3]; continue
                    end
                end
                L=(P,aLinDir)
                aFPOnC=_occ_pcurve_value(pc,aFPar)
                if _lin2d_sqdist(L,aFPOnC)>aTolParConf2
                    aLPOnC=_occ_pcurve_value(pc,aLPar)
                    if _lin2d_sqdist(L,aLPOnC)>aTolParConf2
                        if isSmallAngle
                            done,exts=_occ_extpc2d(P,pc,aFPar,aLPar)
                            if done && !isempty(exts)
                                bi=argmin([s for (_,_,s) in exts])
                                aPOnC=exts[bi][2]
                                aMinDist=exts[bi][3]
                                aFDist=_sqdist2(P,aFPOnC)
                                aLDist=_sqdist2(P,aLPOnC)
                                if aMinDist>aFDist
                                    aMinDist=aFDist; aPOnC=aFPOnC
                                end
                                if aMinDist>aLDist
                                    aMinDist=aLDist; aPOnC=aLPOnC
                                end
                                if aMinDist<Par
                                    Par=aMinDist
                                    Par<aTolParConf2 &&
                                        (e.curEdgePar+=_OCC_PROBING[3]; continue)
                                    aLinDir=_dir2(_sub2(aPOnC,P))
                                    L=(P,aLinDir)
                                end
                            end
                        end
                        e.curEdgePar+=_OCC_PROBING[3]
                        if e.curEdgePar>=_OCC_PROBING[2]
                            e.curEdgeInd+=1
                            e.curEdgePar=_OCC_PROBING[1]
                        end
                        return true,L,sqrt(Par)
                    end
                end
            end
            e.curEdgePar+=_OCC_PROBING[3]
        end
        e.curEdgeInd+=1
        e.curEdgePar=_OCC_PROBING[1]
    end
    return false,(P,(1.0,0.0)),floatmax(Float64)
end
function _occ_fex_segment!(e::_OccFaceExplorer,P)
    e.curEdgeInd=1; e.curEdgePar=_OCC_PROBING[1]
    return _occ_fex_other_segment!(e,P)
end

# ── BRepClass_FClassifier / FaceClassifier drivers ──────────────────────────
# `TopClass_FaceClassifier::Perform` — returns (rejected, nowires,
# classifier_state, edgeparam). `Reject`/`RejectWire`/`RejectEdge` are
# `Standard_False` in `BRepClass_FaceExplorer`.
function _occ_face_classify(m::GeoModel,stag::Int,g,P,Tol)
    aPoint=P
    e=_OccFaceExplorer(m,stag,g)
    aResOfPointCheck=false
    while !aResOfPointCheck
        aResOfPointCheck,aPoint=_occ_fex_checkpoint!(e,aPoint)
    end
    valid,L,aParam=_occ_fex_segment!(e,aPoint)
    nowires=true
    c=_OccClassifier()
    aState=:unknown
    while valid
        _occ_classifier_reset!(c,L,aParam,Tol)
        for loop in get(m.surfaces,stag,Int[])
            nowires=false
            for etag in m.loops[loop]
                Or=etag>0 ? :forward : :reversed
                _occ_classifier_compare!(m,g,c,etag,Or)
                aState=c.myState
                aState===:on && return false,false,aState,c.myParam
            end
            aState=c.myState
            aState===:out && return false,false,aState,c.myParam
        end
        aState=c.myState
        !c.myIsHeadOrEnd && aState!==:unknown && break
        valid,L,aParam=_occ_fex_other_segment!(e,aPoint)
    end
    return false,nowires,aState,c.myParam
end

# `BRepClass_FClassifier::State` — rejected → OUT, no wires → IN, else the
# classifier state.
function _occ_face_state(m::GeoModel,stag::Int,g,P,Tol)
    rejected,nowires,st,_=_occ_face_classify(m,stag,g,P,Tol)
    rejected && return :out
    nowires && return :in
    return st
end

# `BRepClass_FaceClassifier::Perform(F, gp_Pnt, Tol)` — `UVBounds` bounds,
# `Extrema_ExtPS` projection, min-square-distance extremum's `(u,v)`, then
# the 2-D classify. No extrema / solver failure → `:out` (rejected).
function _occ_face_classify3d(m::GeoModel,stag::Int,g,P,Tol)
    u1,u2,v1,v2=_occ_uvbounds(m,stag,g)
    acc=_occ_surface_extrema(g,P,u1,u2,v1,v2,Tol,Tol)
    (acc===nothing || isempty(acc)) && return :out
    best=acc[1]
    for cand in Iterators.drop(acc,1)
        cand[4]<best[4] && (best=cand)
    end
    return _occ_face_state(m,stag,g,(best[1],best[2]),Tol)
end

# `OCCFace::containsPoint` — `state == IN || state == ON` with
# `BRep_Tool::Tolerance` = `Precision::Confusion()` for materialized faces.
function _occ_surface_contains(m::GeoModel,stag::Int,g,P)
    st=_occ_face_classify3d(m,stag,g,P,_OCC_CONFUSION)
    return st===:in || st===:on
end

# `OCCFace::containsParam` — the 2-D classifier state IN/ON.
function _occ_surface_contains_param(m::GeoModel,stag::Int,g,uv)
    st=_occ_face_state(m,stag,g,uv,_OCC_CONFUSION)
    return st===:in || st===:on
end

# `Extrema_ExtPElC::Perform(P, gp_Lin, Tol, Uinf, Usup)` — the single chord
# extremum, accepted inside [Uinf−Tol, Usup+Tol] (Tol = `Precision::Confusion`,
# a length for the unit-direction line parameter), then `Extrema_GExtPC`'s
# `mytolu` post-check (line resolution = Tol again, non-periodic).
function _occ_extpe_line(g,P,uinf,usup)
    # gp_Vec(OR, P)·gp_Vec(L.Direction()) — the stored `gp_Dir`, no
    # renormalization.
    u=_dot3(g.dir,_arc_sub(P,g.origin))
    if u>=uinf-_OCC_CONFUSION && u<=usup+_OCC_CONFUSION
        # `OR.Translated(u·dir)` — scalar product rounded before the add.
        xyz=_add3(g.origin,_occ_mul(g.dir,u))
        return [(u,xyz,_sqdist(xyz,P))]
    end
    return Tuple{Float64,NTuple{3,Float64},Float64}[]
end

# `Extrema_ExtPElC::Perform(P, gp_Circ, Tol, Uinf, Usup)` plus the
# `Extrema_GExtPC` post-check: project P into the circle's plane, take the
# `AngleWithRef` extrema pair, `AdjustPeriodic` into [Uinf, Uinf+2π), the
# boundary snap, then `TolU = Tol/R` acceptance and the `mytolu` re-check.
# `nothing` when the projected point sits on the axis (continuum).
function _occ_extpe_circle(g,P,uinf,usup)
    O,Axe,X,R=g.center,g.n,g.X,g.r
    tolU=R>floatmin(Float64) ? _OCC_CONFUSION/R : Inf
    Pp=_add3(P,_occ_mul(Axe,-_dot3(_arc_sub(P,O),Axe)))
    OPp=_arc_sub(Pp,O)
    _occ_modulus(OPp)<_OCC_CONFUSION && return nothing
    usol1=_occ_angle_with_ref(X,OPp,Axe)
    usol1+π<_OCC_ANGULAR && (usol1=-π)
    usol1-π>-_OCC_ANGULAR && (usol1=π)
    usol2=usol1+π
    _,usol1=_occ_adjust_periodic(uinf,uinf+2π,tolU,uinf,usol1)
    _,usol2=_occ_adjust_periodic(uinf,uinf+2π,tolU,uinf,usol2)
    (usol1-2π-uinf)<tolU && (usol1-2π-uinf)>-tolU && (usol1=uinf)
    (usol2-2π-uinf)<tolU && (usol2-2π-uinf)>-tolU && (usol2=uinf)
    mytolu=R>_OCC_CONFUSION/2 ? 2.0*_gm_asin(_OCC_CONFUSION/(2.0*R)) : 2π
    out=Tuple{Float64,NTuple{3,Float64},Float64}[]
    for us in (usol1,usol2)
        if (uinf-us)<tolU && (us-usup)<tolU
            U=_occ_in_period(us,uinf,uinf+2π)
            if U>=uinf-mytolu && U<=usup+mytolu
                xyz=_occ_circle_point(g,U)
                push!(out,(U,xyz,_sqdist(xyz,P)))
            end
        end
    end
    return out
end

# `OCCEdge::_project` — the `GeomAPI_ProjectPointOnCurve` answer on the
# materialized edge kinds: padded parameter bounds (skipped when the edge
# is closed), the ExtPElC candidates, first-minimum selection. Returns
# (u, xyz) or `nothing` when no extremum is accepted (→ generic fallbacks).
function _occ_curve_project(m::GeoModel,g,tag::Int,P,caller::AbstractString)
    if g.occ===:degenerate
        return nothing                    # null `_curve` → projector fails
    end
    uinf,usup=g.t0,g.t1
    first_point,last_point=m.curves[tag]
    if first_point!=last_point            # _v0 != _v1 → pad the bounds
        du=usup-uinf
        ut=max(abs(du)*1e-8,1e-12)
        uinf-=ut; usup+=ut
    end
    if g.occ===:line
        acc=_occ_extpe_line(g,P,uinf,usup)
    elseif g.occ===:circle
        acc=_occ_extpe_circle(g,P,uinf,usup)
    else
        throw(ArgumentError(
            "_occ_curve_project: unsupported OCC curve kind $(g.occ)"))
    end
    (acc===nothing || isempty(acc)) && return nothing
    best=acc[1]
    for cand in Iterators.drop(acc,1)
        cand[3]<best[3] && (best=cand)
    end
    return (best[1],best[2])
end

# `OCCEdge::firstDer` — `BRepLProp_CLProps(prop, 1, 1e-5).D1()`: the analytic
# derivative — `ElCLib::CircleD1`, `Geom_Line::D1` (the stored `gp_Dir`, no
# renormalization), or the pcurve-on-surface chain rule for degenerate edges.
function _occ_curve_first_der(m::GeoModel,g,tag::Int,u::Float64,
                              caller::AbstractString)
    if g.occ===:circle
        return _occ_circle_derivative(g,u)
    elseif g.occ===:line
        return g.dir
    elseif g.occ===:degenerate
        return _occ_degenerate_d1(g,tag,u)
    end
    throw(ArgumentError(
        "_occ_curve_first_der: unsupported OCC curve kind $(g.occ)"))
end

# `GEdge::secondDer` — a central difference of `firstDer` at eps = 1e-3,
# one-sided at the parameter bounds (OCCEdge does not override it).
function _occ_curve_second_der(m::GeoModel,g,tag::Int,u::Float64,
                               caller::AbstractString)
    eps=1e-3
    if u-eps<=g.t0
        let x1=_occ_curve_first_der(m,g,tag,u,caller),
            x2=_occ_curve_first_der(m,g,tag,u+eps,caller)
            return ntuple(i->1000.0*(x2[i]-x1[i]),3)
        end
    elseif u+eps>=g.t1
        let x1=_occ_curve_first_der(m,g,tag,u-eps,caller),
            x2=_occ_curve_first_der(m,g,tag,u,caller)
            return ntuple(i->1000.0*(x2[i]-x1[i]),3)
        end
    end
    let x1=_occ_curve_first_der(m,g,tag,u-eps,caller),
        x2=_occ_curve_first_der(m,g,tag,u+eps,caller)
        return ntuple(i->500.0*(x2[i]-x1[i]),3)
    end
end

# `OCCEdge::curvature` — degenerate returns eps = 1e-15 directly; otherwise
# `BRepLProp_CLProps(prop, 2, eps).Curvature()` — the LProp formula
# sqrt(N)/DD1/sqrt(DD1) on the analytic D1/D2 pair, DD2 ≤ Tol → 0, floored
# to eps on exit.
function _occ_curve_curvature(m::GeoModel,g,tag::Int,u::Float64,
                              caller::AbstractString)
    eps=1e-15
    g.occ===:degenerate && return eps
    d1=_occ_curve_first_der(m,g,tag,u,caller)
    d2=if g.occ===:circle
        _occ_circle_second_derivative(g,u)
    elseif g.occ===:line
        (0.0,0.0,0.0)
    else
        throw(ArgumentError(
            "_occ_curve_curvature: unsupported OCC curve kind $(g.occ)"))
    end
    tol=eps*eps
    dd1=_sqlen(d1); dd2=_sqlen(d2)
    crv=if dd2<=tol
        0.0
    else
        n=_sqlen(_occ_cross(d1,d2))
        t=n/dd1/dd2
        t<=tol ? 0.0 : sqrt(n)/dd1/sqrt(dd1)
    end
    return crv<=eps ? eps : crv
end

# `GEdge::refineProjection` on an OCC curve — same damped Newton as the
# built-in arc port but on the edge's `[t0,t1]` range and OCC derivatives.
function _occ_curve_refine(m::GeoModel,g,tag::Int,q::NTuple{3,Float64},
                           u::Float64,relax::Float64,tol::Float64,
                           lc::Float64,caller::AbstractString)
    maxDist=tol*lc
    dPQ=_arc_sub(_model_curve_point(m,tag,u,caller),q)
    err=_occ_modulus(dPQ)
    iter=0
    while (iter+=1)<=25 && err>maxDist
        der=_occ_curve_first_der(m,g,tag,u,caller)
        du=_dot3(dPQ,der)/_dot3(der,der)
        du<tol && _occ_modulus(dPQ)>maxDist && (du=1.0)
        unew=fma(-relax,du,u)
        # `std::min/max` on a NaN step returns the lower bound (C++ ordering).
        u=isnan(unew) ? g.t0 : clamp(unew,g.t0,g.t1)
        dPQ=_arc_sub(_model_curve_point(m,tag,u,caller),q)
        err=_occ_modulus(dPQ)
    end
    return err<=maxDist,u,err
end

# `GEdge::XYZToU` on an OCC curve: 21 seeds across `[t0,t1]`, relaxed retry,
# lowest-error parameter on failure.
function _occ_curve_xyz_to_u(m::GeoModel,g,tag::Int,q::NTuple{3,Float64},
                             lc::Float64,caller::AbstractString;
                             relax::Float64=1.0)
    errors=Dict{Float64,Float64}()
    step=(g.t1-g.t0)/20.0
    u_try=g.t0; err=Inf
    for i in 0:20
        u_try=fma(step,Float64(i),g.t0)
        ok,u_try,err=_occ_curve_refine(
            m,g,tag,q,u_try,relax,1e-8,lc,caller)
        ok && return true,u_try
        errors[err]=u_try
    end
    if relax>0.1
        ok,u_try=_occ_curve_xyz_to_u(
            m,g,tag,q,lc,caller;relax=0.75*relax)
        ok && return true,u_try
        errors[_occ_modulus(_arc_sub(
            _model_curve_point(m,tag,u_try,caller),q))]=u_try
    end
    return false,errors[minimum(keys(errors))]
end

# `goldenSectionSearch`/`GEdge::closestPoint` on an OCC curve — the
# 100-sample scan plus recursive golden-section minimization on `[t0,t1]`.
function _occ_curve_golden(m::GeoModel,g,tag::Int,q::NTuple{3,Float64},
                           x1::Float64,x2::Float64,x3::Float64,
                           caller::AbstractString)
    golden2=2.0-(1.0+sqrt(5.0))/2.0
    x4=fma(golden2,x3-x2,x2)
    abs(x3-x1)<1e-9*(abs(x2)+abs(x4)) && return (x3+x1)/2
    d4=_occ_modulus(_arc_sub(q,_model_curve_point(m,tag,x4,caller)))
    d2=_occ_modulus(_arc_sub(q,_model_curve_point(m,tag,x2,caller)))
    return d4<d2 ? _occ_curve_golden(m,g,tag,q,x2,x4,x3,caller) :
                   _occ_curve_golden(m,g,tag,q,x4,x2,x1,caller)
end

function _occ_curve_closest(m::GeoModel,g,tag::Int,q::NTuple{3,Float64},
                            caller::AbstractString)
    tmin,tmax=minmax(g.t0,g.t1)
    dt=(tmax-tmin)/99.0
    dmin=1e22; topt=tmin
    for i in 0:99
        t=fma(Float64(i),dt,tmin)
        d=_occ_modulus(_arc_sub(q,_model_curve_point(m,tag,t,caller)))
        d<dmin && (topt=t; dmin=d)
    end
    t=topt==tmin ?
        _occ_curve_golden(m,g,tag,q,topt,topt+dt/2.0,topt+dt,caller) :
        topt==tmax ?
        _occ_curve_golden(m,g,tag,q,topt-dt,topt-dt/2.0,topt,caller) :
        _occ_curve_golden(m,g,tag,q,topt-dt,topt,topt+dt,caller)
    return t,_model_curve_point(m,tag,t,caller)
end

# The Moore–Penrose inverse of the 3×3 `[du; dv; 0]` Jacobian that
# `GFace::XYZtoUV` builds (`invert_singular_matrix3x3` — a 1e-16-cutoff SVD
# pinv; the pinv is mathematically unique). Returns (jac[:,1], jac[:,2]).
function _occ_jac_cols(du::NTuple{3,Float64},dv::NTuple{3,Float64})
    a=_dot3(du,du); b=_dot3(du,dv); c=_dot3(dv,dv)
    disc=sqrt(fma(a-c,a-c,4.0*b*b))
    eig2=(a+c-disc)/2.0
    if eig2>1e-32
        det=fma(a,c,-(b*b))
        j0=ntuple(i->(c*du[i]-b*dv[i])/det,3)
        j1=ntuple(i->(a*dv[i]-b*du[i])/det,3)
        return j0,j1
    elseif a+c>0.0
        inv=1.0/(a+c)
        return _occ_mul(du,inv),_occ_mul(dv,inv)
    end
    return (0.0,0.0,0.0),(0.0,0.0,0.0)
end

# `GFace::XYZtoUV` — the 9×9-seed pseudo-inverse Newton used as the OCC
# fallback (`parFromPoint` with `onSurface`, `convTestXYZ` as passed) and,
# with `onSurface=false`, the `GFace::closestPoint` catch path. Returns
# (converged, u, v); on failure the last iterate values are returned,
# matching the C++ out-parameters.
function _occ_xyz_to_uv(g,P::NTuple{3,Float64},lc::Float64,
                        onSurface::Bool,testXYZ::Bool;relax::Float64=1.0)
    precision=onSurface ? 1e-8 : 1e-3
    maxiter=onSurface ? 25 : 10
    lo,hi=_occ_surface_bounds(g)
    umin,umax,vmin,vmax=lo[1],hi[1],lo[2],hi[2]
    tol=precision*((umax-umin)^2+(vmax-vmin)^2)
    initf=(0.5,0.6,0.4,0.7,0.3,0.8,0.2,1.0,0.0)
    initu=ntuple(i->fma(initf[i],umax-umin,umin),9)
    initv=ntuple(i->fma(initf[i],vmax-vmin,vmin),9)
    U=V=Unew=Vnew=0.0
    err=1.0; err2=Inf; iter=0
    for i in 1:9, j in 1:9
        U=initu[i]; V=initv[j]
        err=1.0; iter=1
        p=_occ_surface_point(g,U,V)
        err2=_occ_modulus(_arc_sub(P,p))
        err2<1e-8*lc && return true,U,V
        while err>tol && iter<maxiter
            p,du,dv=_occ_surface_d1(g,U,V)
            j0,j1=_occ_jac_cols(du,dv)
            r=_arc_sub(P,p)
            Unew=fma(relax,fma(j0[1],r[1],fma(j0[2],r[2],j0[3]*r[3])),U)
            Vnew=fma(relax,fma(j1[1],r[1],fma(j1[2],r[2],j1[3]*r[3])),V)
            ((Unew>umax+tol || Unew<umin-tol) &&
             (Vnew>vmax+tol || Vnew<vmin-tol)) && break
            du2=Unew-U; dv2=Vnew-V
            err=du2*du2+dv2*dv2
            err2=_occ_modulus(_arc_sub(P,p))
            iter+=1; U=Unew; V=Vnew
        end
        if iter<maxiter && err<=tol && Unew<=umax && Vnew<=vmax &&
           Unew>=umin && Vnew>=vmin
            if onSurface && err2>1e-4*lc && testXYZ
                continue
            else
                return true,U,V
            end
        end
    end
    onSurface || return false,U,V
    relax<1e-3 && return false,U,V
    return _occ_xyz_to_uv(g,P,lc,onSurface,testXYZ;relax=0.75*relax)
end

# `OCCFace::closestPoint` with the `GFace::closestPoint` fallback: when the
# projector finds no accepted extremum, Gmsh minimizes the distance over the
# unrestricted surface (ALGLIB L-BFGS from the best of a 10×10 grid plus the
# initial guess). For the analytic kinds the global minimum is a raw
# `ExtPElS` candidate, so the fallback returns the period-mapped argmin over
# the unfiltered candidates; with no candidates at all (degenerate query)
# it falls back to `parFromPoint(p, false)` + evaluation like Gmsh's catch.
function _occ_surface_closest(g,P,lc::Float64)
    proj=_occ_surface_project(g,P)
    proj!==nothing && return proj
    lo,hi=_occ_surface_bounds(g)
    cands=g.occ===:cylinder ? _occ_extpe_cylinder(g,P) :
          g.occ===:sphere   ? _occ_extpe_sphere(g,P) :
          g.occ===:cone     ? _occ_extpe_cone(g,P) :
          g.occ===:torus    ? _occ_extpe_torus(g,P) :
          g.occ===:plane    ? _occ_extpe_plane(g,P) :
          throw(ArgumentError(
              "_occ_surface_closest: unsupported OCC surface kind $(g.occ)"))
    if cands!==nothing && !isempty(cands)
        best=cands[1]
        for cand in Iterators.drop(cands,1)
            cand[4]<best[4] && (best=cand)
        end
        u,v,xyz,_=best
        g.occ!==:plane && (u=_occ_in_period(u,lo[1],lo[1]+2π))
        g.occ===:torus && (v=_occ_in_period(v,lo[2],lo[2]+2π))
        return (u,v,xyz)
    end
    _,u,v=_occ_xyz_to_uv(g,P,lc,false,false)
    return (u,v,_occ_surface_point(g,u,v))
end

# `OCCFace::parFromPoint` — `_project` first; on failure
# `GFace::parFromPoint(qp, onSurface=true, convTestXYZ=true)` supplies the
# parameter even when the iteration does not converge (the forced XYZ test
# only skips candidates that converged in uv but not in space).
function _occ_surface_parameter_on_face(g,P,lc::Float64)
    proj=_occ_surface_project(g,P)
    proj!==nothing && return (proj[1],proj[2])
    # `OCCFace::parFromPoint`'s fallback forces the XYZ convergence test:
    # `GFace::parFromPoint(qp, onSurface=true, convTestXYZ=true)`.
    _,u,v=_occ_xyz_to_uv(g,P,lc,true,true)
    return (u,v)
end

# `InterpolateCurve(der=1)` on a built-in arc — `gmshEdge::firstDer`: a finite
# difference with `fd_eps=1e-8`, collapsing to a one-sided difference at the
# `[0,1]` bounds (`eps1`/`eps2` zero out at the queried end).
function _arc_first_derivative_fd(g,u::Float64)
    eps=1e-8
    eps1=u<eps ? 0.0 : eps
    eps2=u>1.0-eps ? 0.0 : eps
    p0=_arc_point(g,u-eps1); p1=_arc_point(g,u+eps2)
    inv=1.0/(eps1+eps2)
    return ((p1[1]-p0[1])*inv,(p1[2]-p0[2])*inv,(p1[3]-p0[3])*inv)
end

# `GEdge::refineProjection` on a built-in arc: damped Newton on the finite-
# difference derivative, parameter clamped to `[0,1]`, convergence when the
# residual drops under `tol·lc`. The `du < tol` comparison is signed, exactly
# like the C++ source. Returns (converged, u, err).
function _arc_refine_projection(
    g,q::NTuple{3,Float64},u::Float64,relax::Float64,tol::Float64,lc::Float64)
    maxDist=tol*lc
    dPQ=_arc_sub(_arc_point(g,u),q)
    err=sqrt(_sqlen(dPQ))
    iter=0
    while (iter+=1)<=25 && err>maxDist
        der=_arc_first_derivative_fd(g,u)
        du=_dot3(dPQ,der)/_dot3(der,der)
        du<tol && sqrt(_sqlen(dPQ))>maxDist && (du=1.0)
        u=clamp(fma(-relax,du,u),0.0,1.0)
        dPQ=_arc_sub(_arc_point(g,u),q)
        err=sqrt(_sqlen(dPQ))
    end
    return err<=maxDist,u,err
end

@inline _arc_dist(g,q::NTuple{3,Float64},t::Float64) =
    sqrt(_sqlen(_arc_sub(q,_arc_point(g,t))))

# `GEdge::XYZToU` on a built-in arc: 21 evenly spaced seeds through
# `refineProjection`, then a relaxed retry starting from the last seed's
# refined value. Returns (converged, u) — on failure `u` is the lowest-error
# parameter, matching `errorVsParameter.begin()->second`.
function _arc_xyz_to_u(
    g,q::NTuple{3,Float64},lc::Float64;relax::Float64=1.0)
    errors=Dict{Float64,Float64}()
    u_try=0.0; err=Inf
    for i in 0:20
        u_try=fma(0.05,Float64(i),0.0)
        ok,u_try,err=_arc_refine_projection(g,q,u_try,relax,1e-8,lc)
        ok && return true,u_try
        errors[err]=u_try
    end
    if relax>0.1
        ok,u_try=_arc_xyz_to_u(g,q,lc;relax=0.75*relax)
        ok && return true,u_try
        errors[_arc_dist(g,q,u_try)]=u_try
    end
    return false,errors[minimum(keys(errors))]
end

# `goldenSectionSearch` — Gmsh's recursive golden-section minimizer of
# `|x(t)−q|` on `[0,1]` (`GOLDEN2 = 2−φ`, `tau=1e-9` at the call sites).
function _arc_golden_section(
    g,q::NTuple{3,Float64},x1::Float64,x2::Float64,x3::Float64,tau::Float64)
    golden2=2.0-(1.0+sqrt(5.0))/2.0
    x4=fma(golden2,x3-x2,x2)
    abs(x3-x1)<tau*(abs(x2)+abs(x4)) && return (x3+x1)/2
    return _arc_dist(g,q,x4)<_arc_dist(g,q,x2) ?
        _arc_golden_section(g,q,x2,x4,x3,tau) :
        _arc_golden_section(g,q,x4,x2,x1,tau)
end

# `GEdge::closestPoint` on a built-in arc: a 100-sample scan for the bracket,
# then the recursive golden-section search. Returns (t, point).
function _arc_closest(g,q::NTuple{3,Float64})
    tmin,tmax=0.0,1.0
    dt=(tmax-tmin)/99.0
    dmin=1e22; topt=tmin
    for i in 0:99
        t=fma(Float64(i),dt,tmin)
        d=_arc_dist(g,q,t)
        d<dmin && (topt=t; dmin=d)
    end
    t=topt==tmin ?
        _arc_golden_section(g,q,topt,topt+dt/2.0,topt+dt,1e-9) :
        topt==tmax ?
        _arc_golden_section(g,q,topt-dt,topt-dt/2.0,topt,1e-9) :
        _arc_golden_section(g,q,topt-dt,topt,topt+dt,1e-9)
    return t,_arc_point(g,t)
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
        phi=_gm_atan2(Y[i],g.X[i])
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
    f(θ)=let (sθ,cθ)=_gm_sincos(θ); r1*rho*cθ+r2*sqrt(max(0.0,1.0-rho*rho*sθ*sθ)) end
    best=max(f(a),f(b))
    lo_k=ceil(Int,a/π-1e-12); hi_k=floor(Int,b/π+1e-12)
    for k in lo_k:hi_k
        best=max(best,f(k*π))
    end
    if r1!=r2
        s=(r1*r1-r2*r2*rho*rho)/(rho*rho*(r1*r1-r2*r2))
        if 0.0<=s<1.0
            root=_gm_acos(-sqrt(1.0-s))
            for base in (root,-root), k in -1:1
                θ=base+2k*π
                a-1e-12<=θ<=b+1e-12 || continue
                best=max(best,f(clamp(θ,a,b)))
            end
        end
    end
    return best
end
