module HigherOriginalParametersProbe
using Test,Tessella
const API=Tessella.API
include("../artifacts/higher_original_parameters_primary_v1.jl")
function check_stage(row,stage)
    for node in stage["nodes"]
        coordinates,parameters,dimension,owner=API.mesh.get_node(node["tag"])
        @test coordinates==node["coordinates"]
        @test parameters==node["parameters"]
        @test [dimension,owner]==node["owner"]
    end
    for entity in stage["entities"]
        dimension,owner=entity["entity"]
        tags,coordinates,parameters=API.mesh.get_nodes(dimension,owner,false,true)
        @test tags==entity["node_tags"]
        @test coordinates==entity["coordinates"]
        @test parameters==entity["parameters"]
        kinds,labels,nodes=API.mesh.get_elements(dimension,owner)
        @test kinds==entity["types"]
        @test labels==entity["element_tags"]
        @test nodes==entity["connectivity"]
    end
end
@testset "Actual Higher generation retains original vertex parameters" begin
    for row in HIGHER_ORIGINAL_PARAMETERS_PRIMARY
        @testset "$(row["lower"]) order$(row["order"])" begin
            API.initialize()
            try
                API.option("Mesh.Renumber",0)
                API.option("Mesh.MeshOnlyEmpty",1)
                API.model.add_discrete_entity(0,1)
                API.model.add_discrete_entity(1,11,[1])
                API.model.add_discrete_entity(2,21,[11])
                API.mesh.add_nodes(2,21,row["source_tags"],row["source_coordinates"],row["source_parameters"])
                API.mesh.add_elements_by_type(21,2,[511],row["source_tags"])
                row["lower"]=="surface_only" || API.mesh.add_elements_by_type(11,1,[611],[102,103])
                row["lower"]=="incident_curve_point" && API.mesh.add_elements_by_type(1,15,[711],[102])
                API.mesh.reclassify_nodes()
                row["order"]==2 && API.mesh.set_order(2)
                check_stage(row,row["stages"][1])
                API.mesh.generate(2)
                check_stage(row,row["stages"][2])
            finally
                API.finalize()
            end
        end
    end
end
end
