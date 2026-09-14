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
    :line=>"Line", :circle=>"Circle", :ellipse=>"Ellipse")
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
