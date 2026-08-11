# Stage-4 measured finding (2/2): fine-surface-only tetrahedralization leaves long
# INTERIOR diagonals — so refining the boundary is not sufficient for a size bound.
#
# Mesh the box surface finely (k×k grid per face, watertight), then tetrahedralize.
# The volume (64.0) and validity are correct, but `tetrahedralize` adds NO interior
# points, so a tet spanning the bottom face to the top face keeps maxedge ≈ box
# height (4.0) even at k=4 (98 surface nodes). A genuine size bound needs interior
# Steiner points from boundary-conforming Delaunay refinement, not just a fine input.
#
# Run:  julia --project=. validation/stage4_size_refinement/fine_surface_leaves_interior_diagonals.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella.MeshTypes, Tessella.Geometry
import Tessella.Mesh3D as M3

# Watertight box [0,L]^3 subdivided k×k per face (shared grid nodes on shared edges).
function fine_box(L, k)
    coords = Dict{NTuple{3,Int},Int}(); xs = Float64[]; ys = Float64[]; zs = Float64[]
    h = L / k
    gid(i, j, l) = get!(coords, (i, j, l)) do
        push!(xs, i*h); push!(ys, j*h); push!(zs, l*h); length(xs)
    end
    tris = NTuple{3,Int}[]
    quad(a, b, c, d) = (push!(tris, (a, b, c)); push!(tris, (a, c, d)))
    for i in 0:k-1, j in 0:k-1
        quad(gid(i,j,0), gid(i,j+1,0), gid(i+1,j+1,0), gid(i+1,j,0))   # z=0
        quad(gid(i,j,k), gid(i+1,j,k), gid(i+1,j+1,k), gid(i,j+1,k))   # z=L
        quad(gid(i,0,j), gid(i+1,0,j), gid(i+1,0,j+1), gid(i,0,j+1))   # y=0
        quad(gid(i,k,j), gid(i,k,j+1), gid(i+1,k,j+1), gid(i+1,k,j))   # y=L
        quad(gid(0,i,j), gid(0,i,j+1), gid(0,i+1,j+1), gid(0,i+1,j))   # x=0
        quad(gid(k,i,j), gid(k,i+1,j), gid(k,i+1,j+1), gid(k,i,j+1))   # x=L
    end
    C = Matrix{Float64}(undef, 3, length(xs))
    for n in 1:length(xs); C[1,n]=xs[n]; C[2,n]=ys[n]; C[3,n]=zs[n]; end
    Tm = Matrix{Int32}(undef, 3, length(tris))
    for (t, f) in enumerate(tris); Tm[1,t]=f[1]; Tm[2,t]=f[2]; Tm[3,t]=f[3]; end
    return Mesh(C; tris=Tm)
end

maxedge(m) = begin
    mx = 0.0
    for t in 1:size(m.tets, 2)
        vs = (m.tets[1,t], m.tets[2,t], m.tets[3,t], m.tets[4,t])
        for i in 1:4, j in i+1:4
            p = node(m, vs[i]); q = node(m, vs[j])
            d = sqrt((p[1]-q[1])^2 + (p[2]-q[2])^2 + (p[3]-q[3])^2); d > mx && (mx = d)
        end
    end; mx
end
meshvol(m) = sum(tet_volume(node(m,m.tets[1,t]), node(m,m.tets[2,t]), node(m,m.tets[3,t]), node(m,m.tets[4,t]))
                 for t in 1:size(m.tets, 2); init=0.0)

for k in (1, 2, 3, 4)
    s = fine_box(4.0, k)
    m = M3.tetrahedralize(s); d = validate(m)
    println("k=$k: surf_nodes=$(size(s.coords,2)) surf_tris=$(size(s.tris,2)) | ntets=$(size(m.tets,2)) maxedge=$(round(maxedge(m),digits=3)) vol=$(round(meshvol(m),digits=4)) valid=$(d.ok)")
end
# Observed: maxedge stays ~4.0 (box height) even at k=4 — fine boundary, coarse interior.
