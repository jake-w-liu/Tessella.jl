using Test,Tessella
const API=Tessella.API
const Model=Tessella.Model
include(joinpath(@__DIR__,"../artifacts/reclassification_parameter_lifecycle_primary.jl"))
@testset "Orientation adapters consume original reclassified parameters" begin
    for row in _RECLASSIFICATION_PARAMETER_LIFECYCLE_PRIMARY
        row["action"]=="set_order2" || continue
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
            API.mesh.reclassify_nodes()
            input=Model._model_orientation_record_input(API.CURRENT[])
            for node in row["stages"][2]["nodes"]
                node["tag"] in row["source_tags"] || continue
                @test input.parameters[input.rows[UInt64(node["tag"])]]==node["parameters"]
            end
            API.mesh.set_order(2)
            mesh=API.LAST_MESH[]
            input=Model._model_orientation_mixed_input(API.CURRENT[],mesh)
            for node in row["stages"][3]["nodes"]
                node["tag"] in row["source_tags"] || continue
                @test input.parameters[input.rows[UInt64(node["tag"])]]==node["parameters"]
            end
        finally
            API.finalize()
        end
    end
end
