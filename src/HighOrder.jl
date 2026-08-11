"""
    HighOrder

Stage-6 high-order (curved) elements (PLAN.md §4 Stage 6). This module generates
**quadratic (10-node, P2) tetrahedra** from a linear tet [`Mesh`](@ref) by adding
one shared node at the midpoint of every edge, in gmsh's type-11 node order
(4 corners, then the 6 edge nodes on edges `(1,2),(2,3),(3,1),(1,4),(2,4),(3,4)`).

Straight-sided P2 is geometrically identical to the linear mesh (edge nodes are
exact midpoints) — the total volume is unchanged — and is what a P2 FEM solver
consumes. [`curve_to_cylinder!`](@ref) / [`curve_to_surface!`](@ref) then curve the
edge nodes onto the true boundary geometry for curved-boundary accuracy.

Curving is done safely: only mid-nodes on genuine **boundary-surface edges** are
moved (never interior chords — moving those tangles the volume), and every
projection is reverted if it would drive an incident element's Jacobian
non-positive at the P2 sample nodes, so the result never contains an inverted
element (a strong sampled validity guard, not a formal global-positivity proof).
"""
module HighOrder

using ..MeshTypes: Mesh, node, tet_volume, boundary_faces
import ..MeshTypes: nnodes, ntets     # extended for P2Mesh
using Printf: @printf

export P2Mesh, p2_tetmesh, p2_volume, write_msh_p2,
       curve_to_cylinder!, curve_to_surface!, p2_min_jacobian

"""
    P2Mesh(coords, tet10)

Quadratic tet mesh: `coords` is `3 × N` (original corners followed by edge-mid
nodes); `tet10` is `10 × ntet` in gmsh type-11 order.
"""
struct P2Mesh
    coords::Matrix{Float64}
    tet10::Matrix{Int32}
end

nnodes(p::P2Mesh) = size(p.coords, 2)
ntets(p::P2Mesh) = size(p.tet10, 2)

"""
    p2_tetmesh(m::Mesh) -> P2Mesh

Convert a linear tet mesh to straight-sided quadratic tets, sharing one mid-node
per edge (so the node count grows by exactly the number of unique edges).
"""
function p2_tetmesh(m::Mesh)
    nn = nnodes(m); nt = ntets(m)
    nt == 0 && return P2Mesh(copy(m.coords), Matrix{Int32}(undef, 10, 0))
    edgeid = Dict{Tuple{Int32,Int32}, Int32}()
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    function mid(a::Int32, b::Int32)
        key = a < b ? (a, b) : (b, a)
        get!(edgeid, key) do
            pa = node(m, a); pb = node(m, b)
            push!(xs, (pa[1]+pb[1])/2); push!(ys, (pa[2]+pb[2])/2); push!(zs, (pa[3]+pb[3])/2)
            Int32(nn + length(xs))
        end
    end
    tet10 = Matrix{Int32}(undef, 10, nt)
    @inbounds for t in 1:nt
        v1=m.tets[1,t]; v2=m.tets[2,t]; v3=m.tets[3,t]; v4=m.tets[4,t]
        tet10[1,t]=v1; tet10[2,t]=v2; tet10[3,t]=v3; tet10[4,t]=v4
        tet10[5,t]=mid(v1,v2); tet10[6,t]=mid(v2,v3); tet10[7,t]=mid(v3,v1)
        tet10[8,t]=mid(v1,v4); tet10[9,t]=mid(v2,v4); tet10[10,t]=mid(v3,v4)
    end
    N = nn + length(xs)
    coords = Matrix{Float64}(undef, 3, N)
    @inbounds for i in 1:nn; p=node(m,i); coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]; end
    @inbounds for i in eachindex(xs); coords[1,nn+i]=xs[i]; coords[2,nn+i]=ys[i]; coords[3,nn+i]=zs[i]; end
    return P2Mesh(coords, tet10)
end

"""Total volume of a straight-sided P2 mesh (from the 4 corner nodes of each tet)."""
function p2_volume(p::P2Mesh)
    v = 0.0
    @inbounds for t in 1:ntets(p)
        a=(p.coords[1,p.tet10[1,t]],p.coords[2,p.tet10[1,t]],p.coords[3,p.tet10[1,t]])
        b=(p.coords[1,p.tet10[2,t]],p.coords[2,p.tet10[2,t]],p.coords[3,p.tet10[2,t]])
        c=(p.coords[1,p.tet10[3,t]],p.coords[2,p.tet10[3,t]],p.coords[3,p.tet10[3,t]])
        d=(p.coords[1,p.tet10[4,t]],p.coords[2,p.tet10[4,t]],p.coords[3,p.tet10[4,t]])
        v += tet_volume(a,b,c,d)
    end
    return v
end

# ════════════════════════════════════════════════════════════════════════════════
# P2 element validity — Jacobian of the isoparametric map at the sample nodes.
# ════════════════════════════════════════════════════════════════════════════════
# Reference tet corners N1=(0,0,0),N2=(1,0,0),N3=(0,1,0),N4=(0,0,1);
# barycentric L1=1-r-s-t, L2=r, L3=s, L4=t.  Shape functions (type-11 order):
#   corner i : Li(2Li-1)  ->  grad = (4Li-1)·grad(Li)
#   edge(a,b): 4·La·Lb    ->  grad = 4·(La·grad(Lb) + Lb·grad(La))
# with grad(L1)=(-1,-1,-1), grad(L2)=(1,0,0), grad(L3)=(0,1,0), grad(L4)=(0,0,1).
@inline function _p2_grads(r::Float64, s::Float64, t::Float64)
    L1=1.0-r-s-t; L2=r; L3=s; L4=t
    g1=(-(4L1-1), -(4L1-1), -(4L1-1))
    g2=(4L2-1, 0.0, 0.0)
    g3=(0.0, 4L3-1, 0.0)
    g4=(0.0, 0.0, 4L4-1)
    # edge(a,b) = 4(La·grad Lb + Lb·grad La)
    g5=(4*(L1*1.0 + L2*(-1.0)), 4*(L2*(-1.0)),           4*(L2*(-1.0)))          # (1,2)
    g6=(4*(L3*1.0),             4*(L2*1.0),               0.0)                    # (2,3)
    g7=(4*(L3*(-1.0)),          4*(L1*1.0 + L3*(-1.0)),   4*(L3*(-1.0)))          # (3,1)
    g8=(4*(L4*(-1.0)),          4*(L4*(-1.0)),            4*(L1*1.0 + L4*(-1.0))) # (1,4)
    g9=(4*(L4*1.0),             0.0,                      4*(L2*1.0))             # (2,4)
    g10=(0.0,                   4*(L4*1.0),               4*(L3*1.0))             # (3,4)
    return (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)
end

# The degree-3 Lagrange lattice of the reference tet — barycentric (i,j,k,l)/3 with
# i+j+k+l=3 → 20 nodes (4 corners, 2 per edge, 1 per face). These are exactly where
# a curved P2 edge tangles an incident element, so they make a discriminating
# sampled validity check. Precomputed shape-gradients at each node (const).
const _JAC_NODES = let pts = NTuple{3,Float64}[]
    for j in 0:3, k in 0:3-j, l in 0:3-j-k
        push!(pts, (j/3, k/3, l/3))   # (r,s,t) = (L2,L3,L4)
    end
    Tuple(_p2_grads(r,s,t) for (r,s,t) in pts)
end

@inline function _detJ(coords, tet10, t::Integer, grads)
    j11=0.0;j12=0.0;j13=0.0;j21=0.0;j22=0.0;j23=0.0;j31=0.0;j32=0.0;j33=0.0
    @inbounds for k in 1:10
        v=tet10[k,t]; x=coords[1,v]; y=coords[2,v]; z=coords[3,v]; g=grads[k]
        j11+=x*g[1]; j12+=x*g[2]; j13+=x*g[3]
        j21+=y*g[1]; j22+=y*g[2]; j23+=y*g[3]
        j31+=z*g[1]; j32+=z*g[2]; j33+=z*g[3]
    end
    return j11*(j22*j33-j23*j32) - j12*(j21*j33-j23*j31) + j13*(j21*j32-j22*j31)
end

# min Jacobian determinant of P2 tet `t` over the degree-3 sample nodes.
@inline function _p2_min_detJ(coords, tet10, t::Integer)
    mn = Inf
    for grads in _JAC_NODES
        d = _detJ(coords, tet10, t, grads)
        d < mn && (mn = d)
    end
    return mn
end

"""
    p2_min_jacobian(p::P2Mesh) -> Float64

Minimum, over all elements and the degree-3 sample nodes, of the P2 isoparametric
Jacobian determinant. `> 0` ⇒ no sampled inversion (all elements valid at the
sample nodes); `≤ 0` ⇒ at least one tangled/inverted curved element. For a
straight-sided mesh this equals `6·min tet volume > 0`.
"""
function p2_min_jacobian(p::P2Mesh)
    ntets(p) == 0 && return 0.0
    mn = Inf
    for t in 1:ntets(p)
        d = _p2_min_detJ(p.coords, p.tet10, t)
        d < mn && (mn = d)
    end
    return mn
end

# ════════════════════════════════════════════════════════════════════════════════
# Curving edge nodes onto boundary geometry (safe: boundary edges only + validity)
# ════════════════════════════════════════════════════════════════════════════════
# mid-node → its two corner endpoints, and mid-node → incident tet list.
function _mid_maps(p::P2Mesh)
    endp = Dict{Int32, Tuple{Int32,Int32}}()
    midtets = Dict{Int32, Vector{Int32}}()
    slots = ((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,2,4),(10,3,4))
    @inbounds for t in 1:ntets(p), (sl,i,j) in slots
        m = p.tet10[sl,t]
        endp[m] = (p.tet10[i,t], p.tet10[j,t])
        push!(get!(() -> Int32[], midtets, m), Int32(t))
    end
    return endp, midtets
end

# undirected corner-edges that lie on a boundary surface face (of the linear tets).
function _boundary_corner_edges(p::P2Mesh)
    bf, _ = boundary_faces(@view p.tet10[1:4, :])
    edges = Set{Tuple{Int32,Int32}}()
    for f in bf
        a,b,c = f[1],f[2],f[3]
        push!(edges, minmax(a,b)); push!(edges, minmax(b,c)); push!(edges, minmax(a,c))
    end
    return edges
end

# bounding-box diagonal → relative displacement scale for the "actually moved" test.
function _bbox_diag(p::P2Mesh)
    lo1=lo2=lo3=Inf; hi1=hi2=hi3=-Inf
    @inbounds for i in 1:nnodes(p)
        x=p.coords[1,i]; y=p.coords[2,i]; z=p.coords[3,i]
        x<lo1 && (lo1=x); x>hi1 && (hi1=x)
        y<lo2 && (lo2=y); y>hi2 && (hi2=y)
        z<lo3 && (lo3=z); z>hi3 && (hi3=z)
    end
    sqrt((hi1-lo1)^2 + (hi2-lo2)^2 + (hi3-lo3)^2)
end

# Core curving. `qualifies(a,b)::Bool` gates on the two corner endpoints; `project`
# snaps a point onto the surface. Only boundary-surface edges are curved, and any
# projection that would make an incident element's min sampled Jacobian ≤ 0 is
# reverted (so the mesh never gains an inverted element). Deterministic (mid-nodes
# processed in ascending id order). Returns the number of nodes actually moved.
function _curve_boundary!(p::P2Mesh, qualifies, project; rtol::Float64=1e-6)
    ntets(p) == 0 && return 0
    bnd = _boundary_corner_edges(p)
    endp, midtets = _mid_maps(p)
    movetol2 = (rtol * _bbox_diag(p))^2
    ncurved = 0
    for mid in sort!(collect(keys(endp)))
        a, b = endp[mid]
        (minmax(a,b) in bnd) || continue        # boundary-surface edge only
        qualifies(a, b) || continue             # geometric predicate on endpoints
        @inbounds ox=p.coords[1,mid]; oy=p.coords[2,mid]; oz=p.coords[3,mid]
        q = project(ox, oy, oz)
        (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) || continue
        @inbounds begin p.coords[1,mid]=q[1]; p.coords[2,mid]=q[2]; p.coords[3,mid]=q[3] end
        ok = true
        for t in midtets[mid]
            if _p2_min_detJ(p.coords, p.tet10, t) <= 0.0; ok = false; break; end
        end
        if ok
            ((q[1]-ox)^2 + (q[2]-oy)^2 + (q[3]-oz)^2) > movetol2 && (ncurved += 1)
        else
            @inbounds begin p.coords[1,mid]=ox; p.coords[2,mid]=oy; p.coords[3,mid]=oz end
        end
    end
    return ncurved
end

"""
    curve_to_cylinder!(p::P2Mesh, center, axis, radius; rtol=1e-6) -> n_curved

Curve the quadratic mesh onto a cylinder of the given `center`, `axis`, `radius`:
every **boundary-surface** edge whose *both* corner endpoints lie on the lateral
wall (radial distance `radius` from the axis, within `rtol`) has its mid-node
projected radially onto the exact cylinder — turning the straight chord into a
true curved P2 edge. Interior chords are never touched, and any projection that
would invert an incident element is reverted (see module docstring). Flat caps and
axial spokes stay straight. Returns the number of mid-nodes moved.
"""
function curve_to_cylinder!(p::P2Mesh, center, axis, radius::Real; rtol::Real=1e-6)
    R = Float64(radius)
    c = (Float64(center[1]), Float64(center[2]), Float64(center[3]))
    ez = _unit3((Float64(axis[1]), Float64(axis[2]), Float64(axis[3])))
    tol = Float64(rtol) * R
    @inline function radial(x, y, z)
        rx=x-c[1]; ry=y-c[2]; rz=z-c[3]
        ax = rx*ez[1] + ry*ez[2] + rz*ez[3]
        vx=rx-ax*ez[1]; vy=ry-ax*ez[2]; vz=rz-ax*ez[3]
        (sqrt(vx*vx+vy*vy+vz*vz), ax, vx, vy, vz)
    end
    onwall(v::Integer) = abs(radial(p.coords[1,v], p.coords[2,v], p.coords[3,v])[1] - R) <= tol
    qualifies(a::Integer, b::Integer) = onwall(a) && onwall(b)
    function project(x, y, z)
        d, ax, vx, vy, vz = radial(x, y, z)
        d == 0 && return (Inf, Inf, Inf)          # on the axis — unprojectable, skip
        f = R/d
        (c[1] + ax*ez[1] + vx*f, c[2] + ax*ez[2] + vy*f, c[3] + ax*ez[3] + vz*f)
    end
    return _curve_boundary!(p, qualifies, project; rtol=Float64(rtol))
end

@inline function _unit3(a)
    l = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    l == 0 && error("curve_to_cylinder!: zero axis")
    (a[1]/l, a[2]/l, a[3]/l)
end

"""
    curve_to_surface!(p::P2Mesh, project, on_surface; rtol=1e-6) -> n_curved

Curve the quadratic mesh onto an **arbitrary** target surface — the general form
of [`curve_to_cylinder!`](@ref). `project(x,y,z) -> (x,y,z)` snaps a point onto the
surface; `on_surface(x,y,z) -> Bool` tests whether a *corner* lies on it. A
mid-node is curved only when it sits on a genuine **boundary-surface edge** *and*
both its corner endpoints satisfy `on_surface` — interior chords (whose endpoints
may also lie on the surface) are never moved. Any projection that would invert an
incident element, or that returns a non-finite point, is reverted, so the result
never contains an inverted element. Returns the number of mid-nodes actually moved
(displacement > `rtol` × bounding-box diagonal).
"""
function curve_to_surface!(p::P2Mesh, project, on_surface; rtol::Real=1e-6)
    onsurf(v::Integer) = @inbounds on_surface(p.coords[1,v], p.coords[2,v], p.coords[3,v])
    qualifies(a::Integer, b::Integer) = onsurf(a) && onsurf(b)
    return _curve_boundary!(p, qualifies, (x,y,z) -> project(x,y,z); rtol=Float64(rtol))
end

"""
    write_msh_p2(path, p::P2Mesh; tet_tag=zeros) -> path

Write a quadratic tet mesh as gmsh MSH v2.2 with 10-node (type-11) elements.
"""
function write_msh_p2(path::AbstractString, p::P2Mesh;
                      tet_tag::AbstractVector{<:Integer}=zeros(Int32, ntets(p)))
    length(tet_tag) == ntets(p) || throw(ArgumentError("write_msh_p2: tet_tag length mismatch"))
    open(path, "w") do io
        println(io, "\$MeshFormat"); println(io, "2.2 0 8"); println(io, "\$EndMeshFormat")
        println(io, "\$Nodes"); println(io, nnodes(p))
        @inbounds for i in 1:nnodes(p)
            @printf(io, "%d %.17g %.17g %.17g\n", i, p.coords[1,i], p.coords[2,i], p.coords[3,i])
        end
        println(io, "\$EndNodes")
        println(io, "\$Elements"); println(io, ntets(p))
        @inbounds for t in 1:ntets(p)
            tag = tet_tag[t]
            print(io, t, " 11 2 ", tag, " ", tag)
            for k in 1:10; print(io, " ", p.tet10[k,t]); end
            println(io)
        end
        println(io, "\$EndElements")
    end
    return path
end

end # module HighOrder
