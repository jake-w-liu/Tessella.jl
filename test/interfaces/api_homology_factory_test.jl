module NativeHomologyFactoryProbe
using Tessella,Test
include("../../validation/homology_native/factory_primary_v1.jl")
const API=Tessella.API
@testset "Integral chains retain retired factory maxima" begin
    for row in HOMOLOGY_NATIVE_FACTORY_ROWS
        API.initialize()
        try
            API.option("Mesh.Renumber",0)
            API.model.add_discrete_entity(1,1)
            API.mesh.add_nodes(1,1,[11,13],[0.,0.,0.,1.,0.,0.])
            API.mesh.add_elements_by_type(1,1,[101],[11,13])
            API.model.add_physical_group(1,[1];tag=10)
            API.model.add_discrete_entity(0,99)
            API.mesh.add_elements_by_type(99,15,[9001],[11])
            API.mesh.remove_elements(0,99,[9001])
            API.mesh.add_homology_request(row["kind"],[10],[],[0])
            @test API.mesh.compute_homology()==[Tuple(pair) for pair in row["returned"]]
            entity=only(API.model.get_entities_for_physical_group(0,11))
            @test entity==row["entity"]
            types,elements,connectivity=API.mesh.get_elements(0,entity)
            @test types==row["types"]
            @test elements==row["element_tags"]
            @test connectivity==row["connectivity"]
            @test API.mesh.get_nodes(0,entity,false,false)[1]==row["node_tags"]
            @test API.model.get_physical_name(0,11)==row["physical_name"]
        finally
            API.finalize()
        end
    end
end
# The native high-tag Homology postprocessing path access-violates at this
# boundary. Julia retains safe UInt64 mesh/view records. This checks the
# integral-chain and size_t birth contract, not a successful native DLL run.
@testset "Integral chains safely retain births beyond Int32" begin
    API.initialize()
    try
        API.model.add_discrete_entity(1,1)
        API.mesh.add_nodes(1,1,[11,13],[0.,0.,0.,1.,0.,0.])
        maximum=UInt64(typemax(Int32))
        API.mesh.add_elements_by_type(1,1,[maximum],[11,13])
        API.model.add_physical_group(1,[1];tag=10)
        API.mesh.add_homology_request("Homology",[10],[],[0])
        @test API.mesh.compute_homology()==[(0,11)]
        entity=only(API.model.get_entities_for_physical_group(0,11))
        @test API.mesh.get_elements(1,1)==(Int32[1],[UInt64[maximum]],[UInt64[11,13]])
        @test API.mesh.get_elements(0,entity)==(Int32[15],[UInt64[maximum+1]],[UInt64[13]])
        @test API.mesh.get_max_element_tag()==maximum+1
        @test API.mesh.get_node(11)[1]==[0.,0.,0.]
        @test API.mesh.get_node(13)[1]==[1.,0.,0.]
        @test API.model.get_physical_groups()==[(0,11),(1,10)]
        @test API.model.get_physical_name(0,11)=="H_0{10}1"
        @test API.model.get_entities_for_physical_group(1,10)==[1]
        API.mesh.clear_homology_requests()
        @test API.mesh.compute_homology()==Tuple{Int,Int}[]
    finally
        API.finalize()
    end
end
end
