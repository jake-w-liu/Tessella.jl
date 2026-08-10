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

using ..MeshTypes: Mesh, ntris, nnodes, node, triangle_area
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..Mesh1D: mesh_segment
using ..SizeField: AbstractSizeField, size_at, ConstantSize

export mesh_planar_face, mesh_cylinder_face, mesh_parametric_face
export PlaneFrame, plane_frame, project, lift

# ── vector helpers (tuples) ─────────────────────────────────────────────────────
@inline _sub(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _norm(a) = sqrt(_dot(a,a))
@inline function _unit(a)
    n = _norm(a); n == 0 && error("MeshSurface: zero-length vector"); (a[1]/n,a[2]/n,a[3]/n)
end

# ── planar frame ────────────────────────────────────────────────────────────────
struct PlaneFrame
    o::NTuple{3,Float64}   # origin (a loop vertex)
    u::NTuple{3,Float64}   # in-plane unit axis 1
    v::NTuple{3,Float64}   # in-plane unit axis 2
    n::NTuple{3,Float64}   # unit normal
end

"""Orthonormal in-plane frame of a coplanar loop (Newell's method for the normal)."""
function plane_frame(loop)
    m = length(loop)
    m >= 3 || throw(ArgumentError("plane_frame: loop needs ≥3 points"))
    nx=ny=nz=0.0
    @inbounds for i in 1:m
        p = loop[i]; q = loop[mod1(i+1,m)]
        nx += (p[2]-q[2])*(p[3]+q[3])
        ny += (p[3]-q[3])*(p[1]+q[1])
        nz += (p[1]-q[1])*(p[2]+q[2])
    end
    L = sqrt(nx^2+ny^2+nz^2)
    L == 0 && throw(ArgumentError("plane_frame: degenerate (zero-area) loop"))
    n = (nx/L, ny/L, nz/L)
    # in-plane axis: project the world axis least aligned with n
    ax = abs(n[1]) <= abs(n[2]) ?
         (abs(n[1]) <= abs(n[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(n[2]) <= abs(n[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    u = _unit(_cross(ax, n))
    v = _cross(n, u)          # already unit (n,u orthonormal)
    return PlaneFrame(Tuple(Float64.(loop[1])), u, v, n)
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
        p = Tuple(Float64.(loop[i])); q = Tuple(Float64.(loop[mod1(i+1,m)]))
        pts, _ = mesh_segment(p, q, sf)
        # append all but the last node (shared with the next edge / loop closure)
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
    return Mesh(coords; tris=copy(m2.tris), tri_tag=copy(m2.tri_tag))
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
    R = Float64(radius); H = Float64(height)
    (R > 0 && H > 0) || throw(ArgumentError("mesh_cylinder_face: radius and height must be positive"))
    c = Tuple(Float64.(center)); ez = _unit(Tuple(Float64.(axis)))
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
    # circumferential divisions from the smallest size seen around the equator
    hmin = Inf
    for k in 0:11
        hmin = min(hmin, size_at(sf, on(2π*k/12, H/2)))
    end
    ntheta = max(6, round(Int, 2π*R / hmin))
    # structured nodes: level j (1..nz) × sector i (0..ntheta-1)
    coords = Matrix{Float64}(undef, 3, nz*ntheta)
    idx(j, i) = (j-1)*ntheta + (mod(i, ntheta)) + 1
    @inbounds for j in 1:nz, i in 0:ntheta-1
        p = on(2π*i/ntheta, zlev[j]); n = idx(j, i)
        coords[1,n]=p[1]; coords[2,n]=p[2]; coords[3,n]=p[3]
    end
    tris = Matrix{Int32}(undef, 3, 2*(nz-1)*ntheta); nt = 0
    @inbounds for j in 1:nz-1, i in 0:ntheta-1
        a = idx(j,i); b = idx(j,i+1); cc = idx(j+1,i+1); d = idx(j+1,i)
        nt+=1; tris[1,nt]=a; tris[2,nt]=b; tris[3,nt]=cc
        nt+=1; tris[1,nt]=a; tris[2,nt]=cc; tris[3,nt]=d
    end
    return Mesh(coords; tris=tris[:, 1:nt])
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
    umin=Float64(umin); umax=Float64(umax); vmin=Float64(vmin); vmax=Float64(vmax)
    du = 1e-6*(umax-umin); dv = 1e-6*(vmax-vmin)
    # local stretch factors (√E, √G) from finite differences
    stretch(u,v) = begin
        su = _sub(s(u+du,v), s(u-du,v)); sv = _sub(s(u,v+dv), s(u,v-dv))
        (_norm(su)/(2du), _norm(sv)/(2dv))
    end
    # parameter-space size target: physical h divided by stretch (per axis, use min)
    sizefn = (u,v) -> begin
        p = s(u,v); h = size_at(sf, p); (a,b) = stretch(u,v)
        h / max(a, b, 1e-30)
    end
    # PSLG = parameter rectangle boundary, graded in parameter space
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    corners=[(umin,vmin),(umax,vmin),(umax,vmax),(umin,vmax)]
    start=1
    for i in 1:4
        p=corners[i]; q=corners[mod1(i+1,4)]
        # metric length along the edge in physical space
        n = max(1, round(Int, _param_edge_metric(s, sf, p, q)))
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
        p=s(m2.coords[1,i], m2.coords[2,i]); coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]
    end
    return Mesh(coords; tris=copy(m2.tris))
end

@inline function _param_edge_metric(s, sf, p, q)
    N = 32; L = 0.0; prev = s(p[1], p[2])
    @inbounds for k in 1:N
        t = k/N; u = p[1]+t*(q[1]-p[1]); v = p[2]+t*(q[2]-p[2])
        cur = s(u,v); h = size_at(sf, cur)
        L += _norm(_sub(cur, prev))/h; prev = cur
    end
    return L
end

end # module MeshSurface
