module ReclassificationParameterLifecycleTests
using Test,Tessella
const API=Tessella.API
include(joinpath(@__DIR__,"../artifacts/reclassification_parameter_lifecycle_primary.jl"))

function check_stage(row,stage)
    for node in stage["nodes"]
        node["tag"] in row["source_tags"] || continue
        coordinates,parameters,dimension,owner=API.mesh.get_node(node["tag"])
        @test coordinates==node["coordinates"]
        @test parameters==node["parameters"]
        @test [dimension,owner]==node["owner"]
    end
    for entity in stage["entities"]
        dimension,owner=entity["entity"]
        tags,coordinates,parameters=API.mesh.get_nodes(dimension,owner,false,true)
        @test parameters==entity["parameters"]
    end
end

@testset "Order and refinement retain original reclassified node parameters" begin
    for row in _RECLASSIFICATION_PARAMETER_LIFECYCLE_PRIMARY
        @testset "$(row["lower"]) $(row["action"])" begin
            API.initialize()
            try
                API.option("Mesh.Renumber",0)
                API.model.add_discrete_entity(0,1)
                API.model.add_discrete_entity(1,11,[1])
                API.model.add_discrete_entity(2,21,[11])
                API.mesh.add_nodes(2,21,row["source_tags"],row["source_coordinates"],row["source_parameters"])
                API.mesh.add_elements_by_type(21,2,[511],row["source_tags"])
                row["lower"]=="surface_only" || API.mesh.add_elements_by_type(11,1,[611],[102,103])
                row["lower"]=="incident_curve_point" && API.mesh.add_elements_by_type(1,15,[711],[102])
                for stage in row["stages"]
                    operation=stage["stage"]
                    if operation=="reclassify"
                        API.mesh.reclassify_nodes()
                    elseif operation=="set_order2"
                        API.mesh.set_order(2)
                    elseif operation in ("refine","refine_after_order2")
                        API.mesh.refine()
                    end
                    check_stage(row,stage)
                end
            finally
                API.finalize()
            end
        end
    end
end
end
