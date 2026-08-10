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

    @testset "curve_to_cylinder!: mid-nodes land on the true cylinder" begin
        R = 2.0
        m = tetrahedralize(cylinder_surface((0.,0,0),(0.,0,1),R,5.0; nθ=16, nz=3))
        p = p2_tetmesh(m)
        nc = curve_to_cylinder!(p, (0.,0,0), (0.,0,1), R)
        @test nc > 0
        # every mid-node whose both endpoints are on the lateral wall is now exactly
        # at radius R (a genuine curved P2 edge), not on the straight chord.
        slots=((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,2,4),(10,3,4))
        endp=Dict{Int32,Tuple{Int32,Int32}}()
        for t in 1:ntets(p), (sl,i,j) in slots; endp[p.tet10[sl,t]]=(p.tet10[i,t],p.tet10[j,t]); end
        rad(v)=sqrt(p.coords[1,v]^2 + p.coords[2,v]^2)
        maxdev = 0.0
        for (mid,(a,b)) in endp
            (abs(rad(a)-R)<1e-6 && abs(rad(b)-R)<1e-6) || continue
            maxdev = max(maxdev, abs(rad(mid)-R))
        end
        @test maxdev < 1e-12                     # curved nodes exactly on the cylinder
    end

    @testset "empty mesh" begin
        e = Mesh(Matrix{Float64}(undef,3,0))
        p = p2_tetmesh(e)
        @test size(p.tet10, 2) == 0
        @test p2_volume(p) == 0.0
    end
end
