# Regression: the smoothers fixed only nodes on the outer boundary, so nodes on an
# internal region interface (faces between differently tagged tets) and on embedded
# triangles moved off the interface: every tet stayed positive and the total volume
# was preserved, but per-region volumes and the interface geometry silently changed.
using Test
using Tessella
import Tessella.MeshTypes as MT

function _region_volumes(m)
    v = Dict{Int32,Float64}()
    for t in 1:size(m.tets, 2)
        p = ntuple(i -> Tuple(m.coords[:, m.tets[i, t]]), 4)
        v[m.tet_tag[t]] = get(v, m.tet_tag[t], 0.0) + abs(MT.tet_signed_volume(p...))
    end
    v
end

@testset "smoothing keeps region interfaces fixed" begin
    m = mesh_box_regions([BoxRegion(0, 1, 0, 1, 0, 1, 1), BoxRegion(1, 2, 0, 1, 0, 1, 2)]; hmax=0.5)
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
    @test length(onif) > 3
    v0 = _region_volumes(mj)
    # Two sweeps each: the default suite runs interpreted, and a node that
    # may move leaves the interface in the first sweep.
    for f in (x -> smooth_laplacian(x; iters=2), x -> smooth_odt(x; iters=2),
              x -> smooth_optimize(x; iters=2, sliver_deg=60.0))
        m2 = f(mj)
        @test validate(m2).ok
        @test all(m2.coords[:, i] == c[:, i] for i in onif)
        v = _region_volumes(m2)
        @test isapprox(v[1], v0[1]; rtol=1e-12) && isapprox(v[2], v0[2]; rtol=1e-12)
    end
    # nodes of embedded triangles are fixed as well
    me = Tessella.Mesh(c; tets=m.tets, tris=reshape(Int32[onif[1], onif[2], onif[3]], 3, 1))
    @test all(Tessella.Optimize._boundary_nodes(me)[onif[1:3]])
end
