using Test, Tessella

@testset "constant metric controls the actual box grid" begin
    mesh=mesh_box(0,1,0,1,0,1,Metric3(9,1,1,0,0,0))
    @test ntuple(axis->length(unique(mesh.coords[axis,:])),3)==(7,3,3) && size(mesh.tets,2)==144
    volumes=[Tessella.MeshTypes.tet_signed_volume((mesh.coords[:,vertex] for vertex in vertices)...) for vertices in eachcol(mesh.tets)]
    @test all(>(0),volumes) && sum(volumes)≈1.0
    bound=0.0
    for vertices in eachcol(mesh.tets),i in 1:3,j in i+1:4
        difference=mesh.coords[:,vertices[j]]-mesh.coords[:,vertices[i]]
        bound=max(bound,sqrt(9difference[1]^2+difference[2]^2+difference[3]^2))
    end
    @test bound<=1.0
end
