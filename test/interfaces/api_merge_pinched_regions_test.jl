module ApiMergePinchedRegionsTest
# Two Regions may touch along an edge only (OCC boxes fragmented across a
# shared edge: Gmsh 4.15.2 meshes it, 281 nodes / 772 tets at size 0.3;
# Tessella 281 / 771). Their union has vertices whose link is two disks, so
# the merged mesh failed the whole-complex vertex-link certificate with
# "vertex 33 has a disconnected link". The certificate is now taken per
# Region; one Region that is itself pinched is still refused.
using Test, Tessella
const A = Tessella.API

tet(points) = Tessella.Mesh(reduce(hcat, [collect(p) for p in points]); tets=reshape(Int32[1, 2, 3, 4], 4, 1))

@testset "Regions pinched along an edge merge" begin
    # Two positively oriented tetrahedra sharing only the edge (0,0,0)-(0,0,1).
    first = tet([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
    second = tet([(0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)])
    @test Tessella.validate(first).ok && Tessella.validate(second).ok
    merged = A._merge_entity_meshes([(1, first), (2, second)], "pinch test")[1]
    @test size(merged.coords, 2) == 6 && size(merged.tets, 2) == 2
    # The same two cells as one Region are a pinched complex and stay refused.
    single = Tessella.Mesh(merged.coords; tets=merged.tets)
    @test !Tessella.validate(single).ok
    @test_throws ErrorException A._merge_entity_meshes([(1, single)], "pinch test")
end
end
