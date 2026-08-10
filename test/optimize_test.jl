# ── Stage-4 CRC suite: tet-mesh quality reporting + Laplacian smoothing ─────────
#
# Correctness  : quality metrics on a known-good mesh (regular-ish tets); smoothing
#                preserves total volume (boundary fixed) and validity; improves the
#                mean dihedral and does not increase the sliver count.
# Robustness   : slivery Delaunay-of-random-cloud input, boundary nodes pinned.
# Completeness  : no tet inverted (positive-volume guard), tags preserved.

using Test
using Tessella.Mesh3D
using Tessella.MeshTypes
using Tessella.Optimize

mvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)

mutable struct _RO; s::UInt64; end
_nfo(r::_RO) = (r.s ⊻= r.s<<13; r.s ⊻= r.s>>7; r.s ⊻= r.s<<17; (r.s>>11)/Float64(2^53))

@testset "Optimize (Stage 4)" begin

    @testset "mesh_quality on a single unit tet" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,4],4,1))
        q = mesh_quality(m)
        @test q.n_tets == 1
        @test q.min_volume ≈ 1/6 rtol=1e-12
        @test 0 < q.min_dihedral_deg <= q.mean_dihedral_deg <= q.max_dihedral_deg < 180
        @test q.n_slivers == 0                      # a unit right tet is not a sliver
    end

    @testset "smoothing preserves volume + validity, improves mean dihedral" begin
        r = _RO(0xABCDEF)
        n = 300
        xs=Float64[_nfo(r) for _ in 1:n]; ys=Float64[_nfo(r) for _ in 1:n]; zs=Float64[_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        q0 = mesh_quality(m)
        ms = smooth_laplacian(m; iters=10, relax=0.8)
        q1 = mesh_quality(ms)
        @test mvol(ms) ≈ mvol(m) rtol=1e-9          # boundary fixed ⇒ volume preserved
        @test validate(ms).ok                       # no inverted tets
        @test ntets(ms) == ntets(m)                 # topology unchanged
        @test q1.mean_dihedral_deg > q0.mean_dihedral_deg      # mean quality up
        @test q1.n_slivers <= q0.n_slivers          # slivers not increased
    end

    @testset "boundary nodes are pinned (cube unchanged)" begin
        # cube surface fill has only boundary nodes → smoothing is a no-op geometrically
        C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
        F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
        ct=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); ct[:,k]=Int32[f...]; end
        m = tetrahedralize(Mesh(C; tris=ct))
        ms = smooth_laplacian(m; iters=5)
        @test ms.coords ≈ m.coords                  # all nodes on the boundary ⇒ fixed
        @test mvol(ms) ≈ 1.0 rtol=1e-6
    end

    @testset "tags preserved through smoothing" begin
        r = _RO(0x1234)
        n = 120
        xs=Float64[_nfo(r) for _ in 1:n]; ys=Float64[_nfo(r) for _ in 1:n]; zs=Float64[_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        # slap arbitrary tags on
        tags = Int32[(t % 3) + 1 for t in 1:ntets(m)]
        mt = Mesh(m.coords; tets=m.tets, tet_tag=tags)
        ms = smooth_laplacian(mt; iters=3)
        @test ms.tet_tag == tags
    end
end
