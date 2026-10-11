using Test, Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_occ_box_binding_v1.jl"))
const _BG_MODEL = Tessella.Model

function _bg_execute(origin, extent)
    source = "SetFactory(\"OpenCASCADE\"); Box(31) = {" * join(vcat(origin, extent), ",") * "};\n"
    mktempdir() do directory
        path = joinpath(directory, "signed_box.geo")
        write(path, source)
        return Tessella.GeoExec.execute_geo(path)
    end
end

@testset "Real GEO signed Box factory bindings" begin
    for (index, row) in enumerate(BOX_SIGNED_ROWS)
        @testset "signed GEO row $index" begin
            execution = _bg_execute(row["origin"], row["extents"])
            model = execution.model
            @test execution.msg_error_count == 0
            @test collect(model.box_extents[31]) == vcat(row["lower"], row["actual_extents"])
            for point in row["points"]
                @test collect(model.points[point["tag"]]) == point["value"]
            end
            for edge in row["edges"]
                lower, upper = _BG_MODEL.model_parametrization_bounds(model, 1, edge["tag"])
                @test lower == edge["bounds"][1]
                @test upper == edge["bounds"][2]
                @test _BG_MODEL.model_value(model, 1, edge["tag"], vcat(lower, upper)) == vcat(edge["values"]...)
            end
            for face in row["faces"]
                lower, upper = _BG_MODEL.model_parametrization_bounds(model, 2, face["tag"])
                @test lower == face["bounds"][1]
                @test upper == face["bounds"][2]
                @test isapprox(_BG_MODEL.model_normal(model, face["tag"], face["center_uv"]), face["normal"]; atol=2e-11, rtol=0)
                @test isapprox(_BG_MODEL.model_value(model, 2, face["tag"], face["center_uv"]), face["value"]; atol=2e-11, rtol=0)
            end
        end
    end
end

@testset "Real GEO Box two-corner reconstruction" begin
    for (index, row) in enumerate(BOX_ENDPOINT_ROWS)
        @testset "arithmetic GEO row $index" begin
            execution = _bg_execute([row["origin"],0,0], [row["width"],1,1])
            model = execution.model
            @test execution.msg_error_count == 0
            @test model.box_extents[31] == (row["minimum"], 0.0, 0.0, row["actual_delta"], 1.0, 1.0)
            for point in row["points"]
                @test collect(model.points[point["tag"]]) == point["value"]
            end
            for (dimension, field) in ((1,"edge_domains"),(2,"face_domains"))
                for entity in row[field]
                    lower, upper = _BG_MODEL.model_parametrization_bounds(model, dimension, entity["tag"])
                    @test lower == entity["bounds"][1]
                    @test upper == entity["bounds"][2]
                end
            end
        end
    end
end
