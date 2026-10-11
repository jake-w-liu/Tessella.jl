using Test, Tessella

@testset "constant metric keeps material tags on one shared grid" begin
    mesh=mesh_box_regions([BoxRegion(0,1,0,1,0,1,3),BoxRegion(1,2,0,1,0,1,7)],
                          Metric3(4,1,1,0,0,0);max_metric_edge=2)
    @test size(mesh.tets,2)==24 && Set(mesh.tet_tag)==Set(Int32[3,7])
    volumes=[Tessella.MeshTypes.tet_signed_volume((mesh.coords[:,vertex] for vertex in vertices)...)
             for vertices in eachcol(mesh.tets)]
    @test all(>(0),volumes) && sum(volumes[mesh.tet_tag.==3])≈1 &&
        sum(volumes[mesh.tet_tag.==7])≈1
end
