module ChainPhysicalNamesProbe
using Tessella,Test
const API=Tessella.API
@testset "Published chain cells retain their physical names in the actual mesh" begin
    API.initialize()
    try
        API.model.add_discrete_entity(1,1)
        API.mesh.add_nodes(1,1,[11,12,13,14],[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.])
        API.mesh.add_elements_by_type(1,1,[101,102,103,104],[11,12,12,13,13,14,14,11])
        API.model.add_physical_group(1,[1];tag=10)
        API.mesh.reclassify_nodes()
        API.mesh.add_homology_request("Homology",[10],[],[1])
        @test API.mesh.compute_homology()==[(1,11)]
        @test API.model.get_physical_name(1,11)=="H_1{10}1"
        @test get(API.LAST_MESH[].physical_names,(1,11),nothing)=="H_1{10}1"
    finally
        API.finalize()
    end
end
end
