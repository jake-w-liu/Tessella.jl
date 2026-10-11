module HomologyIntegralProtocolProbe
using Tessella,Test
include("../../validation/homology_native/integral_primary_v1.jl")
const API=Tessella.API
function check_stage(expected)
    @test API.model.get_entities()==[(e["dimension"],e["tag"]) for e in expected["entities"]]
    for entity in expected["entities"]
        dim,tag=entity["dimension"],entity["tag"]
        nodes,coordinates,parameters=API.mesh.get_nodes(dim,tag,false,true)
        @test nodes==entity["node_tags"]
        @test coordinates==entity["coordinates"]
        @test parameters==entity["parameters"]
        types,cells,connectivity=API.mesh.get_elements(dim,tag)
        @test types==entity["types"]
        @test cells==entity["element_tags"]
        @test connectivity==entity["connectivity"]
    end
    for node in expected["nodes"]
        coordinates,parameters,dim,owner=API.mesh.get_node(node["tag"])
        @test coordinates==node["coordinates"]
        @test parameters==node["parameters"]
        @test (dim,owner)==(node["dimension"],node["owner"])
    end
    @test API.model.get_physical_groups()==[(p["dimension"],p["tag"]) for p in expected["physical"]]
    for physical in expected["physical"]
        dim,tag=physical["dimension"],physical["tag"]
        @test API.model.get_entities_for_physical_group(dim,tag)==physical["members"]
        @test API.model.get_physical_name(dim,tag)==physical["name"]
    end
end
@testset "Integral homology cell families, orientation and exact representatives" begin
    for row in HOMOLOGY_INTEGRAL_PROTOCOL_ROWS
        @testset "$(row["name"]) $(row["kind"])" begin
            API.initialize()
            try
                dimension=row["dimension"]
                API.option("Mesh.Renumber",0)
                API.model.add_discrete_entity(dimension,1)
                API.mesh.add_nodes(dimension,1,row["input_node_tags"],row["input_coordinates"])
                connectivity=reduce(vcat,row["input_connectivity"];init=Int[])
                API.mesh.add_elements_by_type(1,row["msh_type"],row["input_element_tags"],connectivity)
                API.model.add_physical_group(dimension,[1];tag=10)
                check_stage(row["baseline"])
                API.mesh.add_homology_request(row["kind"],[10],[],row["dimensions"])
                output=API.mesh.compute_homology()
                @test output==[Tuple(pair) for pair in row["returned_physical"]]
                check_stage(row["immediate"])
                API.mesh.reclassify_nodes()
                check_stage(row["reclassified"])
            finally
                API.finalize()
            end
        end
    end
end
end
