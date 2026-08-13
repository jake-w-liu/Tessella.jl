"""
    MeshSurface

Stage-2 surface meshing (PLAN.md §3, "mesh the flat/patch/coax faces"). Two
routes, both reducing a surface to the exact 2-D machinery of [`Mesh2D`](@ref):

* **Planar faces** — a face bounded by coplanar loops (outer + holes). Boundary
  loops are first discretized under the size field ([`Mesh1D`](@ref)), projected
  to an in-plane orthonormal frame, meshed as a PSLG (constrained Delaunay +
  size-driven Ruppert refinement), then lifted back to 3-D. Isometric ⇒ the 2-D
  quality/size guarantees hold exactly on the surface.

* **Parametric faces** `s(u,v)` — meshed in the parameter rectangle. For a
  **developable** surface (plane, cylinder, cone) the parameter map is (up to a
  constant scale per axis) isometric, so plain 2-D meshing of the scaled
  parameter domain gives exact surface sizing — this is how the coax cylinders
  are meshed. For a general surface the parameter-space size target is divided by
  the local stretch of the first fundamental form (an isotropic metric
  approximation), verified on the sphere by area/where analytic.

All meshers return a 3-D triangle [`Mesh`](@ref).
"""
module MeshSurface

using ..MeshTypes: Mesh, ntris, nnodes, node, triangle_area, validate
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..Mesh1D: mesh_segment
using ..SizeField: AbstractSizeField, size_at, ConstantSize

export mesh_planar_face, mesh_cylinder_face, mesh_parametric_face
export PlaneFrame, plane_frame, project, lift

# ── vector helpers (tuples) ─────────────────────────────────────────────────────
@inline _sub(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _norm(a) = hypot(a[1],a[2],a[3])
@inline function _unit(a)
    n = _norm(a)
    (isfinite(n) && n > 0) ||
        throw(ArgumentError("MeshSurface: vector must have finite positive length (got $a)"))
    (a[1]/n,a[2]/n,a[3]/n)
end

@inline function _point3(p, caller::AbstractString)
    length(p) >= 3 || throw(ArgumentError("$caller: each point needs three coordinates"))
    q = try
        (Float64(p[1]), Float64(p[2]), Float64(p[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: point coordinates must be Float64-convertible: $(sprint(showerror, err))"))
    end
    (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) ||
        throw(ArgumentError("$caller: point has non-finite coordinates $q"))
    return q
end

# ── planar frame ────────────────────────────────────────────────────────────────
struct PlaneFrame
    o::NTuple{3,Float64}   # origin (a loop vertex)
    u::NTuple{3,Float64}   # in-plane unit axis 1
    v::NTuple{3,Float64}   # in-plane unit axis 2
    n::NTuple{3,Float64}   # unit normal
end

function _validate_coplanar_loops(loops, fr::PlaneFrame)
    scale = 1.0
    for loop in loops, p0 in loop
        p = _point3(p0, "mesh_planar_face")
        dx=p[1]-fr.o[1];dy=p[2]-fr.o[2];dz=p[3]-fr.o[3]
        (isfinite(dx)&&isfinite(dy)&&isfinite(dz)) ||
            throw(ArgumentError("mesh_planar_face: coordinate span overflowed Float64"))
        scale = max(scale, abs(dx),abs(dy),abs(dz))
    end
    tol = 256eps(Float64)*scale
    for (li, loop) in enumerate(loops), (pi, p0) in enumerate(loop)
        p = _point3(p0, "mesh_planar_face")
        dist = abs(_dot(_sub(p, fr.o), fr.n))
        isfinite(dist) || throw(ArgumentError("mesh_planar_face: plane-distance evaluation overflowed Float64"))
        dist <= tol ||
            throw(ArgumentError("mesh_planar_face: loop $li point $pi is not coplanar (plane distance $dist > tolerance $tol)"))
    end
    return nothing
end

"""Orthonormal in-plane frame of a coplanar loop (Newell's method for the normal)."""
function plane_frame(loop)
    m = length(loop)
    m >= 3 || throw(ArgumentError("plane_frame: loop needs ≥3 points"))
    nx=ny=nz=0.0
    @inbounds for i in 1:m
        p = _point3(loop[i], "plane_frame"); q = _point3(loop[mod1(i+1,m)], "plane_frame")
        nx += (p[2]-q[2])*(p[3]+q[3])
        ny += (p[3]-q[3])*(p[1]+q[1])
        nz += (p[1]-q[1])*(p[2]+q[2])
    end
    L = hypot(nx,ny,nz)
    (isfinite(L)&&L>0) || throw(ArgumentError("plane_frame: degenerate or numerically unrepresentable loop"))
    n = (nx/L, ny/L, nz/L)
    # in-plane axis: project the world axis least aligned with n
    ax = abs(n[1]) <= abs(n[2]) ?
         (abs(n[1]) <= abs(n[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(n[2]) <= abs(n[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    u = _unit(_cross(ax, n))
    v = _cross(n, u)          # already unit (n,u orthonormal)
    return PlaneFrame(_point3(loop[1], "plane_frame"), u, v, n)
end

@inline project(fr::PlaneFrame, p) = (_dot(_sub(p, fr.o), fr.u), _dot(_sub(p, fr.o), fr.v))
@inline lift(fr::PlaneFrame, a, b) =
    (fr.o[1]+a*fr.u[1]+b*fr.v[1], fr.o[2]+a*fr.u[2]+b*fr.v[2], fr.o[3]+a*fr.u[3]+b*fr.v[3])

# ── planar face meshing ─────────────────────────────────────────────────────────
"""
    mesh_planar_face(loops, sf; min_angle_deg=25.0, max_area=Inf) -> Mesh

Mesh a planar face bounded by `loops` (a vector of loops; `loops[1]` is the outer
boundary, the rest are holes). Each loop is a vector of coplanar 3-D points in
order. Boundary edges are size-graded, then the interior is constrained-Delaunay
meshed and Ruppert-refined under `sf`, and lifted back to 3-D.
"""
function mesh_planar_face(loops::AbstractVector, sf::AbstractSizeField;
                          min_angle_deg::Real=25.0, max_area::Real=Inf)
    isempty(loops) && throw(ArgumentError("mesh_planar_face: no loops"))
    fr = plane_frame(loops[1])
    _validate_coplanar_loops(loops, fr)
    xs = Float64[]; ys = Float64[]; segs = Tuple{Int,Int}[]
    for loop in loops
        _add_loop!(xs, ys, segs, loop, sf, fr)
    end
    T = constrained_delaunay(xs, ys, segs)
    sizefn = (a, b) -> size_at(sf, lift(fr, a, b))
    interior = refine!(T; min_angle_deg=min_angle_deg, max_area=max_area, size=sizefn)
    m2 = to_mesh(T; interior=interior)
    return _lift_mesh(m2, fr)
end

# discretize one loop's edges under sf and append its cyclic boundary as segments
function _add_loop!(xs, ys, segs, loop, sf, fr)
    m = length(loop)
    m >= 3 || throw(ArgumentError("mesh_planar_face: a loop needs ≥3 points"))
    start = length(xs) + 1
    for i in 1:m
        p = _point3(loop[i], "mesh_planar_face")
        q = _point3(loop[mod1(i+1,m)], "mesh_planar_face")
        pts, _ = mesh_segment(p, q, sf)
        # append all but the last node (shared with the next edge / loop closure)
        length(pts)-1 <= typemax(Int32)-length(xs) ||
            throw(ArgumentError("mesh_planar_face: boundary node count exceeds Int32 indexing"))
        @inbounds for k in 1:length(pts)-1
            a, b = project(fr, pts[k]); push!(xs, a); push!(ys, b)
        end
    end
    stop = length(xs)
    # cyclic segments over this loop's appended nodes
    for i in start:stop
        j = i == stop ? start : i+1
        push!(segs, (i, j))
    end
    return nothing
end

function _lift_mesh(m2::Mesh, fr::PlaneFrame)
    nn = nnodes(m2)
    coords = Matrix{Float64}(undef, 3, nn)
    @inbounds for i in 1:nn
        a = m2.coords[1,i]; b = m2.coords[2,i]
        p = lift(fr, a, b)
        coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]
    end
    out=Mesh(coords; tris=copy(m2.tris), tri_tag=copy(m2.tri_tag))
    d=validate(out)
    d.ok || throw(ErrorException("mesh_planar_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    return out
end

# ── developable cylinder face ───────────────────────────────────────────────────
"""
    mesh_cylinder_face(center, axis, radius, height, sf; kwargs...) -> Mesh

Mesh the lateral surface of a cylinder: axis start `center`, direction `axis`,
`radius`, `height`. Built as a **graded structured** grid — axial levels are
size-graded (via [`Mesh1D`](@ref) along the axis), the circumference is divided
uniformly into `nθ = max(6, round(2πR/h))` sectors, and each quad is split into
two triangles with the seam wrapping (θ_{nθ} ≡ θ_0). Watertight by construction
(no seam to weld), which a Ruppert mesh of the unrolled rectangle cannot
guarantee once the two seam edges refine asymmetrically.
"""
function mesh_cylinder_face(center, axis, radius::Real, height::Real,
                            sf::AbstractSizeField; min_angle_deg::Real=25.0, max_area::Real=Inf)
    R = _surface_float(radius,"mesh_cylinder_face","radius")
    H = _surface_float(height,"mesh_cylinder_face","height")
    (isfinite(R) && isfinite(H) && R > 0 && H > 0) ||
        throw(ArgumentError("mesh_cylinder_face: radius and height must be finite and positive (got $radius, $height)"))
    angle=_surface_float(min_angle_deg,"mesh_cylinder_face","min_angle_deg")
    (0<=angle<=atand(inv(sqrt(2.0)))) ||
        throw(ArgumentError("mesh_cylinder_face: min_angle_deg must be in [0, $(atand(inv(sqrt(2.0))))] for the structured right-triangle grid"))
    area=_surface_float(max_area,"mesh_cylinder_face","max_area")
    (!isnan(area)&&area>0) ||
        throw(ArgumentError("mesh_cylinder_face: max_area must be positive or Inf"))
    c = _point3(center, "mesh_cylinder_face")
    ez = _unit(_point3(axis, "mesh_cylinder_face"))
    ax = abs(ez[1]) <= abs(ez[2]) ?
         (abs(ez[1]) <= abs(ez[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(ez[2]) <= abs(ez[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    ex = _unit(_cross(ax, ez)); ey = _cross(ez, ex)
    on(θ, z) = (c[1]+R*cos(θ)*ex[1]+R*sin(θ)*ey[1]+z*ez[1],
                c[2]+R*cos(θ)*ex[2]+R*sin(θ)*ey[2]+z*ez[2],
                c[3]+R*cos(θ)*ex[3]+R*sin(θ)*ey[3]+z*ez[3])
    # axial levels: size-graded along the axis (uses the real 3-D size field)
    apts, _ = mesh_segment(c, (c[1]+H*ez[1], c[2]+H*ez[2], c[3]+H*ez[3]), sf)
    zlev = Float64[_dot(_sub(p, c), ez) for p in apts]     # 0 .. H, graded
    nz = length(zlev)
    # Circumferential divisions use the smallest sampled target over all axial
    # levels, not only the equator (an axial size field must also refine rings).
    hmin = Inf
    for z in zlev,k in 0:11
        hmin = min(hmin, size_at(sf, on(2π*k/12,z)))
    end
    mindz=minimum(zlev[i+1]-zlev[i] for i in 1:length(zlev)-1)
    (isfinite(mindz)&&mindz>0) ||
        throw(ArgumentError("mesh_cylinder_face: axial grading produced a non-positive interval"))
    tana=tand(angle)
    chord_target=hmin
    if isfinite(area)
        chord_target=min(chord_target,angle==0 ? sqrt(2area) : sqrt(area/tana))
    end
    angle>0 && (chord_target=min(chord_target,mindz/tana))
    (isfinite(chord_target)&&chord_target>0) ||
        throw(ArgumentError("mesh_cylinder_face: requested controls are below Float64 resolution"))
    sectors = (2π)*(R/chord_target)
    sectors <= prevfloat(Float64(typemax(Int))) ||
        throw(ArgumentError("mesh_cylinder_face: requested circumferential node count exceeds the platform Int limit"))
    ntheta = max(6, ceil(Int, sectors))
    chord=R*(2sinpi(inv(ntheta)))
    (isfinite(chord)&&chord>0) ||
        throw(ArgumentError("mesh_cylinder_face: circumferential chord is below Float64 resolution"))

    # Subdivide every axial interval enough to satisfy the requested right-triangle
    # angle and area controls.  The supported angle ceiling makes adjacent feasible
    # subdivision intervals overlap, so this deterministic ceil rule has no gaps.
    refined_z=Float64[zlev[1]]
    for i in 1:length(zlev)-1
        dz=zlev[i+1]-zlev[i]
        maxdz=angle==0 ? Inf : chord/tana
        isfinite(area) && (maxdz=min(maxdz,2area/chord))
        ratio=dz/maxdz
        (!isfinite(maxdz) || (isfinite(ratio)&&ratio<=prevfloat(Float64(typemax(Int))))) ||
            throw(ArgumentError("mesh_cylinder_face: axial subdivision count exceeds the platform Int limit"))
        nsub=isfinite(maxdz) ? ceil(Int,ratio) : 1
        nsub=max(1,nsub);piece=dz/nsub
        angle>0 && piece < chord*tana*(1-128eps(Float64)) &&
            throw(ErrorException("mesh_cylinder_face: could not satisfy the requested angle bound at Float64 resolution"))
        for k in 1:nsub
            α=k/nsub
            push!(refined_z,k==nsub ? zlev[i+1] : zlev[i]*(1-α)+zlev[i+1]*α)
        end
    end
    zlev=refined_z;nz=length(zlev)
    nnout=try
        Base.checked_mul(nz, ntheta)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_cylinder_face: requested mesh dimensions overflow the platform Int limit"))
    end
    nnout<=typemax(Int32) ||
        throw(ArgumentError("mesh_cylinder_face: $nnout nodes exceed the Int32 indexing limit"))
    ntout=try
        Base.checked_mul(Base.checked_mul(2, nz-1), ntheta)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_cylinder_face: requested mesh dimensions overflow the platform Int limit"))
    end
    ntout<=typemax(Int32) ||
        throw(ArgumentError("mesh_cylinder_face: $ntout triangles exceed the Int32 topology limit"))
    # structured nodes: level j (1..nz) × sector i (0..ntheta-1)
    coords = Matrix{Float64}(undef, 3, nnout)
    idx(j, i) = (j-1)*ntheta + (mod(i, ntheta)) + 1
    @inbounds for j in 1:nz, i in 0:ntheta-1
        p = on(2π*i/ntheta, zlev[j]); n = idx(j, i)
        coords[1,n]=p[1]; coords[2,n]=p[2]; coords[3,n]=p[3]
    end
    tris = Matrix{Int32}(undef, 3, ntout); nt = 0
    @inbounds for j in 1:nz-1, i in 0:ntheta-1
        a = idx(j,i); b = idx(j,i+1); cc = idx(j+1,i+1); d = idx(j+1,i)
        nt+=1; tris[1,nt]=a; tris[2,nt]=b; tris[3,nt]=cc
        nt+=1; tris[1,nt]=a; tris[2,nt]=cc; tris[3,nt]=d
    end
    out=Mesh(coords; tris=tris[:, 1:nt])
    d=validate(out)
    d.ok || throw(ErrorException("mesh_cylinder_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    @inbounds for t in 1:ntris(out)
        a=node(out,out.tris[1,t]);b=node(out,out.tris[2,t]);c3=node(out,out.tris[3,t])
        A=triangle_area(a,b,c3);amin=_triangle_min_angle_deg(a,b,c3)
        A<=area*(1+256eps(Float64)) ||
            throw(ErrorException("mesh_cylinder_face: postcondition failed: triangle $t exceeds max_area"))
        amin+1024eps(Float64)>=angle ||
            throw(ErrorException("mesh_cylinder_face: postcondition failed: triangle $t violates min_angle_deg"))
    end
    return out
end

function _triangle_min_angle_deg(a,b,c)
    sides=sort([_norm(_sub(a,b)),_norm(_sub(b,c)),_norm(_sub(c,a))])
    sides[1]>0 || return 0.0
    q=sides./sides[3]
    return acosd(clamp((q[2]^2+1-q[1]^2)/(2q[2]),-1.0,1.0))
end

# ── general parametric face (isotropic-metric approximation) ─────────────────────
"""
    mesh_parametric_face(s, umin, umax, vmin, vmax, sf; nu=32, nv=32, kwargs...) -> Mesh

Mesh a parametric surface `s(u,v) -> (x,y,z)` over the rectangle
`[umin,umax]×[vmin,vmax]`. The parameter-space size target is the physical size
divided by the local stretch of the first fundamental form (isotropic metric
approximation via the mean of √E, √G), then the domain is 2-D meshed and mapped
to 3-D. Suitable for gently-curved patches; developable surfaces are exact via
[`mesh_cylinder_face`](@ref).
"""
function mesh_parametric_face(s, umin::Real, umax::Real, vmin::Real, vmax::Real,
                              sf::AbstractSizeField; min_angle_deg::Real=25.0, max_area::Real=Inf)
    umin=_surface_float(umin,"mesh_parametric_face","umin")
    umax=_surface_float(umax,"mesh_parametric_face","umax")
    vmin=_surface_float(vmin,"mesh_parametric_face","vmin")
    vmax=_surface_float(vmax,"mesh_parametric_face","vmax")
    (isfinite(umin) && isfinite(umax) && isfinite(vmin) && isfinite(vmax) &&
     umin < umax && vmin < vmax) ||
        throw(ArgumentError("mesh_parametric_face: require finite bounds umin < umax and vmin < vmax"))
    du = 1e-6*(umax-umin); dv = 1e-6*(vmax-vmin)
    (isfinite(du)&&isfinite(dv)&&du>0&&dv>0) ||
        throw(ArgumentError("mesh_parametric_face: parameter spans are below Float64 finite-difference resolution"))
    surf(u,v) = _point3(s(u,v), "mesh_parametric_face")
    # local stretch factors (√E, √G) from finite differences
    stretch(u,v) = begin
        ul = max(umin, u-du); ur = min(umax, u+du)
        vl = max(vmin, v-dv); vr = min(vmax, v+dv)
        su = _sub(surf(ur,v), surf(ul,v)); sv = _sub(surf(u,vr), surf(u,vl))
        lu=_norm(su);lv=_norm(sv);a=lu/(ur-ul);b=lv/(vr-vl)
        uu=(su[1]/lu,su[2]/lu,su[3]/lu);vv=(sv[1]/lv,sv[2]/lv,sv[3]/lv)
        sinang=_norm(_cross(uu,vv))
        (isfinite(a) && isfinite(b) && isfinite(sinang) && a>0 && b>0 &&
         sinang > 256eps(Float64)) ||
            throw(ArgumentError("mesh_parametric_face: singular or non-finite surface Jacobian at ($u,$v)"))
        (a,b)
    end
    # parameter-space size target: physical h divided by stretch (per axis, use min)
    sizefn = (u,v) -> begin
        p = surf(u,v); h = size_at(sf, p); (a,b) = stretch(u,v)
        hp=h/max(a,b)
        (isfinite(hp)&&hp>0) || throw(ArgumentError(
            "mesh_parametric_face: parameter-space size is not representable at ($u,$v)"))
        hp
    end
    # PSLG = parameter rectangle boundary, graded in parameter space
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    corners=[(umin,vmin),(umax,vmin),(umax,vmax),(umin,vmax)]
    start=1
    for i in 1:4
        p=corners[i]; q=corners[mod1(i+1,4)]
        # metric length along the edge in physical space
        metric = _param_edge_metric(surf, sf, p, q)
        (isfinite(metric)&&metric>0&&metric<=prevfloat(Float64(typemax(Int)))) ||
            throw(ArgumentError("mesh_parametric_face: requested boundary node count exceeds the platform Int limit"))
        n = max(1, round(Int, metric))
        n<=typemax(Int32)-length(xs) ||
            throw(ArgumentError("mesh_parametric_face: boundary node count exceeds Int32 indexing"))
        for k in 0:n-1
            t=k/n; push!(xs, p[1]+t*(q[1]-p[1])); push!(ys, p[2]+t*(q[2]-p[2]))
        end
    end
    stop=length(xs)
    for i in start:stop; push!(segs, (i, i==stop ? start : i+1)); end
    T=constrained_delaunay(xs,ys,segs)
    interior=refine!(T; min_angle_deg=min_angle_deg, max_area=max_area, size=sizefn)
    m2=to_mesh(T; interior=interior)
    # map parameter mesh to 3-D
    nn=nnodes(m2); coords=Matrix{Float64}(undef,3,nn)
    @inbounds for i in 1:nn
        p=surf(m2.coords[1,i], m2.coords[2,i]); coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]
    end
    out=Mesh(coords; tris=copy(m2.tris))
    d=validate(out)
    d.ok || throw(ErrorException("mesh_parametric_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    return out
end

@inline function _param_edge_metric(s, sf, p, q)
    N = 32; L = 0.0; prev = s(p[1], p[2]);hprev=size_at(sf,prev)
    @inbounds for k in 1:N
        t = k/N; u = p[1]+t*(q[1]-p[1]); v = p[2]+t*(q[2]-p[2])
        cur = s(u,v); h = size_at(sf, cur)
        distance=_norm(_sub(cur,prev))
        increment=distance==0 ? 0.0 : (distance/hprev)/2+(distance/h)/2
        L += increment; prev = cur;hprev=h
        isfinite(L) || throw(ArgumentError("mesh_parametric_face: boundary metric length overflowed Float64"))
    end
    return L
end

function _surface_float(x::Real,caller::AbstractString,name::AbstractString)
    try
        return Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be representable as Float64"))
    end
end

end # module MeshSurface
