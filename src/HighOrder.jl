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

export P2Mesh, p2_tetmesh, p2_volume, write_msh_p2

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
