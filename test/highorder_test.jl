# ── Stage-6 CRC suite: quadratic (P2) tet generation + type-11 I/O ──────────────
#
# Correctness  : exactly one shared mid-node per edge (node count = corners + edges,
#                cross-checked against MeshTypes.unique_edges); corners preserved;
#                edge nodes are exact midpoints; straight P2 volume == linear volume.
# Robustness   : single tet, filled volume mesh, empty mesh.
# Completeness : gmsh type-11 file writes and reads back to the same linear
#                connectivity (CRC) — solver-consumable.

using Test
using Tessella.MeshTypes
using Tessella.Mesh3D
using Tessella.Geometry
using Tessella.HighOrder
using Tessella.IO

lvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)

@testset "HighOrder P2 (Stage 6)" begin

    @testset "single tet → 10 nodes, midpoints, volume preserved" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,4],4,1))
        p = p2_tetmesh(m)
        @test size(p.tet10) == (10, 1)
        @test size(p.coords, 2) == 4 + 6            # 4 corners + 6 edges
        # corners unchanged
        @test p.coords[:,1:4] == m.coords
        # node 5 = midpoint of edge (1,2) = (0.5,0,0)
        @test p.coords[:, p.tet10[5,1]] ≈ [0.5,0.0,0.0]
        @test p2_volume(p) ≈ 1/6 rtol=1e-12
    end

    @testset "shared mid-nodes: node count = corners + unique edges" begin
        s = box_surface(0,1,0,1,0,1)
        m = tetrahedralize(s)
        p = p2_tetmesh(m)
        nedges = length(unique_edges(Matrix{Int32}(undef,3,0), m.tets))
        @test size(p.coords, 2) == nnodes(m) + nedges     # one shared node per edge
        @test p2_volume(p) ≈ lvol(m) rtol=1e-9            # straight P2 ⇒ same volume
    end

    @testset "gmsh type-11 round-trip (reads back to linear connectivity)" begin
        m = tetrahedralize(cylinder_surface((0.,0,0),(0.,0,1),1.0,2.0; nθ=16))
        p = p2_tetmesh(m)
        dir = mktempdir(); path = joinpath(dir, "p2.msh")
        write_msh_p2(path, p; tet_tag=fill(Int32(7), ntets(m)))
        f = read_msh(path)                                # type-11 read as 4-corner tets
        @test ntets(f.mesh) == ntets(m)
        @test mesh_crc(f.mesh).sha == mesh_crc(m).sha     # linear connectivity preserved
        @test all(f.mesh.tet_tag .== 7)
        @test validate(f.mesh).ok
    end

    @testset "empty mesh" begin
        e = Mesh(Matrix{Float64}(undef,3,0))
        p = p2_tetmesh(e)
        @test size(p.tet10, 2) == 0
        @test p2_volume(p) == 0.0
    end
end
