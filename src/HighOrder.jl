"""
    HighOrder

Stage-6 high-order (curved) elements (PLAN.md §4 Stage 6). This module generates
**quadratic (10-node, P2) tetrahedra** from a linear tet [`Mesh`](@ref) by adding
one shared node at the midpoint of every edge, in gmsh's type-11 node order
(4 corners, then the 6 edge nodes on edges `(1,2),(2,3),(3,1),(1,4),(2,4),(3,4)`).

Straight-sided P2 is geometrically identical to the linear mesh (edge nodes are
exact midpoints) — the total volume is unchanged — and is what a P2 FEM solver
consumes. *Curving* the edge nodes onto the true geometry (for curved-boundary
accuracy) needs the geometry model and is the remaining piece.
"""
module HighOrder

using ..MeshTypes: Mesh, node, tet_volume
import ..MeshTypes: nnodes, ntets     # extended for P2Mesh
using Printf: @printf

export P2Mesh, p2_tetmesh, p2_volume, write_msh_p2, curve_to_cylinder!

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

"""
    curve_to_cylinder!(p::P2Mesh, center, axis, radius; rtol=1e-6) -> n_curved

Curve the quadratic mesh onto a cylinder: every edge-mid node whose *both*
endpoints lie on the lateral surface (distance `radius` from the axis) is
projected radially onto the exact cylinder — turning the straight chord into a
true curved P2 edge. Returns the number of nodes moved. This is genuine
curved-boundary high-order for a known primitive; the flat caps are untouched.
"""
function curve_to_cylinder!(p::P2Mesh, center, axis, radius::Real; rtol::Real=1e-6)
    R = Float64(radius)
    c = (Float64(center[1]),Float64(center[2]),Float64(center[3]))
    ez = _unit3((Float64(axis[1]),Float64(axis[2]),Float64(axis[3])))
    tol = rtol * R
    radial(v) = begin
        rel = (p.coords[1,v]-c[1], p.coords[2,v]-c[2], p.coords[3,v]-c[3])
        ax = rel[1]*ez[1]+rel[2]*ez[2]+rel[3]*ez[3]
        rv = (rel[1]-ax*ez[1], rel[2]-ax*ez[2], rel[3]-ax*ez[3])
        (sqrt(rv[1]^2+rv[2]^2+rv[3]^2), ax, rv)
    end
    # mid-node → its two corner endpoints (consistent across tets)
    endp = Dict{Int32, Tuple{Int32,Int32}}()
    slots = ((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,2,4),(10,3,4))
    @inbounds for t in 1:ntets(p), (s,i,j) in slots
        endp[p.tet10[s,t]] = (p.tet10[i,t], p.tet10[j,t])
    end
    ncurved = 0
    for (mid, (a,b)) in endp
        (abs(radial(a)[1]-R) <= tol && abs(radial(b)[1]-R) <= tol) || continue
        d, ax, rv = radial(mid)
        d == 0 && continue                       # degenerate (on the axis) — skip
        f = R/d
        @inbounds begin
            p.coords[1,mid] = c[1] + ax*ez[1] + rv[1]*f
            p.coords[2,mid] = c[2] + ax*ez[2] + rv[2]*f
            p.coords[3,mid] = c[3] + ax*ez[3] + rv[3]*f
        end
        ncurved += 1
    end
    return ncurved
end

@inline function _unit3(a)
    l=sqrt(a[1]^2+a[2]^2+a[3]^2); l==0 && error("curve_to_cylinder!: zero axis"); (a[1]/l,a[2]/l,a[3]/l)
end

"""
    write_msh_p2(path, p::P2Mesh; tags=zeros) -> path

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
