# Regression: the smoothers fixed only nodes on the outer boundary, so nodes on an
# internal region interface (faces between differently tagged tets) and on embedded
# triangles moved off the interface: every tet stayed positive and the total volume
# was preserved, but per-region volumes and the interface geometry silently changed.
#
# In its own module: the suite includes every test file into Main. The mesh is
# built by hand (a 4x2x2 grid of half-unit cubes, six tetrahedra each, split into
# two regions at x = 1) so the test does not pay for a mesher in the interpreted
# default suite; its interface has one node that is free to move.
module OptimizeRegionInterfaceTests
using Test
using Tessella
import Tessella.MeshTypes as MT

function two_region_grid()
    nx, ny, nz = 4, 2, 2
    id(i, j, k) = Int32(i + (nx + 1) * (j + (ny + 1) * k) + 1)
    coords = Matrix{Float64}(undef, 3, (nx + 1) * (ny + 1) * (nz + 1))
    for k in 0:nz, j in 0:ny, i in 0:nx
        coords[:, id(i, j, k)] .= (i / 2, j / 2, k / 2)
    end
    tets = Matrix{Int32}(undef, 4, 6 * nx * ny * nz)
    tags = Vector{Int32}(undef, size(tets, 2))
    column = 0
    for k in 0:nz-1, j in 0:ny-1, i in 0:nx-1
        for order in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
            corner = [i, j, k]
            path = Int32[id(corner...)]
            for axis in order
                corner[axis] += 1
                push!(path, id(corner...))
            end
            p = ntuple(n -> Tuple(coords[:, path[n]]), 4)
            MT.tet_signed_volume(p...) < 0 && ((path[3], path[4]) = (path[4], path[3]))
            tets[:, column+=1] = path
            tags[column] = i < 2 ? 1 : 2
        end
    end
    return Tessella.Mesh(coords; tets=tets, tet_tag=tags)
end

function region_volumes(m)
    v = Dict{Int32,Float64}()
    for t in 1:size(m.tets, 2)
        p = ntuple(i -> Tuple(m.coords[:, m.tets[i, t]]), 4)
        v[m.tet_tag[t]] = get(v, m.tet_tag[t], 0.0) + abs(MT.tet_signed_volume(p...))
    end
    v
end

@testset "smoothing keeps region interfaces fixed" begin
    m = two_region_grid()
    @test validate(m).ok
    # Deterministically displace interior nodes (interface nodes only within x = 1).
    c = copy(m.coords)
    bnd = Tessella.Optimize._boundary_nodes(Tessella.Mesh(m.coords; tets=m.tets))
    s = UInt64(12345)
    r() = (s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17; (Float64(s >> 11) / 2.0^53 - 0.5) * 0.04)
    for i in 1:size(c, 2)
        bnd[i] && continue
        if c[1, i] == 1.0
            c[2, i] += r(); c[3, i] += r()
        else
            c[1, i] += r(); c[2, i] += r(); c[3, i] += r()
        end
    end
    mj = Tessella.Mesh(c; tets=m.tets, tet_tag=m.tet_tag)
    @test validate(mj).ok
    onif = [i for i in 1:size(c, 2) if c[1, i] == 1.0]
    @test length(onif) == 9 && count(i -> !bnd[i], onif) == 1
    v0 = region_volumes(mj)
    for f in (x -> smooth_laplacian(x; iters=2), x -> smooth_odt(x; iters=2),
              x -> smooth_optimize(x; iters=2, sliver_deg=60.0))
        m2 = f(mj)
        @test validate(m2).ok
        @test all(m2.coords[:, i] == c[:, i] for i in onif)
        v = region_volumes(m2)
        @test isapprox(v[1], v0[1]; rtol=1e-12) && isapprox(v[2], v0[2]; rtol=1e-12)
    end
    # nodes of embedded triangles are fixed as well
    me = Tessella.Mesh(c; tets=m.tets, tris=reshape(Int32[onif[1], onif[2], onif[3]], 3, 1))
    @test all(Tessella.Optimize._boundary_nodes(me)[onif[1:3]])
end
end
