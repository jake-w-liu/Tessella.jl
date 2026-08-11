# Stage-4 measured finding (1/2): interior longest-edge midpoint refinement DIVERGES.
#
# Insert the midpoint of every edge longer than `hmax` into the 3-D Delaunay of a
# convex box, repeat. The max edge gets PINNED at 4·√2 ≈ 5.657 (the box face
# diagonal) forever: each pass adds 3 vertices / 7 tets and the max edge never
# drops. Mechanism: the face diagonal is a *boundary edge* of the input surface;
# splitting it at its midpoint (the face centre) regenerates a diagonal of equal
# length among the remaining corners. Interior insertion cannot shorten a fixed
# boundary edge — so this naive scheme never reaches the size bound.
#
# Run:  julia --project=. validation/stage4_size_refinement/diverges_interior_midpoint.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella.MeshTypes, Tessella.Geometry
import Tessella.Mesh3D as M3

function scan(T, hmax)
    best = 0.0; seen = Set{Tuple{Int32,Int32}}(); longs = Tuple{Int32,Int32}[]
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !M3._is_ghost_tet(T, t)) || continue
        vs = (M3._vert(T,t,1), M3._vert(T,t,2), M3._vert(T,t,3), M3._vert(T,t,4))
        for i in 1:4, j in i+1:4
            a = vs[i]; b = vs[j]; key = a < b ? (a, b) : (b, a)
            p = M3._pt(T, a); q = M3._pt(T, b)
            d = sqrt((p[1]-q[1])^2 + (p[2]-q[2])^2 + (p[3]-q[3])^2)
            d > best && (best = d)
            if d > hmax && !(key in seen); push!(seen, key); push!(longs, key); end
        end
    end
    return best, longs
end

s = box_surface(0, 4, 0, 4, 0, 4); nn = size(s.coords, 2)
xs = [s.coords[1,i] for i in 1:nn]; ys = [s.coords[2,i] for i in 1:nn]; zs = [s.coords[3,i] for i in 1:nn]
T = M3.delaunay3d(xs, ys, zs)
hmax = 5.0
for pass in 1:25
    best, longs = scan(T, hmax)
    println("pass=$pass nverts=$(T.nreal) nlivetets=$(M3.ntets_live(T)) maxedge=$(round(best,digits=4)) nlong=$(length(longs))")
    isempty(longs) && (println("CONVERGED"); break)
    for (a, b) in longs
        p = M3._pt(T, a); q = M3._pt(T, b)
        vid = M3._add_vertex3!(T, 0.5*(p[1]+q[1]), 0.5*(p[2]+q[2]), 0.5*(p[3]+q[3]))
        M3.insert_point3!(T, vid)
    end
end
# Observed: maxedge pinned at 5.6569 (= 4√2) for every pass — never reaches hmax=5.0.
