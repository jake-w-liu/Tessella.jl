using Test, Tessella
const _BA_API = Tessella.API
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_occ_box_binding_v1.jl"))
_ba_close(found, expected) = isapprox(found, expected; atol=2e-11, rtol=0)

@testset "Pinned OCC Box admission and binding" begin
@testset "Pinned OCC Box actual span admission" begin
    for (index, row) in enumerate(BOX_ADMISSION_ROWS)
        @testset "admission row $index" begin
            _BA_API.initialize()
            try
                _BA_API.option("Geometry.Tolerance", row["geometry_tolerance"])
                model = _BA_API.CURRENT[]
                before = Dict(name => deepcopy(getfield(model, name)) for name in fieldnames(typeof(model)))
                accepted = true
                result = 0
                try
                    result = _BA_API.model.add_box(row["origin"]..., row["extents"]...; tag=31)
                catch error
                    println("BOX_ADMISSION_ERROR row=", index, " ", sprint(showerror, error))
                    accepted = false
                end
                @test accepted == row["accepted"]
                if row["accepted"] && accepted
                    @test result == 31
                    @test Set(_BA_API.model.get_entities()) == Set(Tuple(entity) for entity in row["entities"])
                    for edge in row["edges"]
                        lower, upper = _BA_API.model.get_parametrization_bounds(1, edge["tag"])
                        @test _ba_close(vcat(lower, upper), edge["bounds"])
                        @test _ba_close(_BA_API.model.get_value(1, edge["tag"], [edge["middle"]]), edge["point"])
                    end
                elseif !row["accepted"]
                    @test all(name -> getfield(model, name) == before[name], keys(before))
                    @test isempty(_BA_API.model.get_entities())
                end
            finally
                _BA_API.finalize()
            end
        end
    end
end

@testset "Pinned OCC Box signed factory and parameter bindings" begin
    for (index, row) in enumerate(BOX_SIGNED_ROWS)
        @testset "signed binding row $index" begin
            _BA_API.initialize()
            try
                accepted = true
                try
                    _BA_API.model.add_box(row["origin"]..., row["extents"]...; tag=31)
                catch error
                    println("BOX_SIGNED_ERROR row=", index, " ", sprint(showerror, error))
                    accepted = false
                end
                @test accepted
                if accepted
                    @test collect(_BA_API.CURRENT[].box_extents[31]) == vcat(row["lower"], row["actual_extents"])
                    for point in row["points"]
                        @test _BA_API.model.get_value(0, point["tag"], Float64[]) == point["value"]
                    end
                    for edge in row["edges"]
                        lower, upper = _BA_API.model.get_parametrization_bounds(1, edge["tag"])
                        @test lower == edge["bounds"][1]
                        @test upper == edge["bounds"][2]
                        @test _BA_API.model.get_value(1, edge["tag"], vcat(lower, upper)) == vcat(edge["values"]...)
                    end
                    for face in row["faces"]
                        lower, upper = _BA_API.model.get_parametrization_bounds(2, face["tag"])
                        @test lower == face["bounds"][1]
                        @test upper == face["bounds"][2]
                        @test _ba_close(_BA_API.model.get_value(2, face["tag"], face["center_uv"]), face["value"])
                        @test _ba_close(_BA_API.model.get_normal(face["tag"], face["center_uv"]), face["normal"])
                        @test Set(_BA_API.model.get_boundary([(2, face["tag"])], false, true, false)) == Set(Tuple(pair) for pair in face["curves"])
                        for sample in face["curve_reparametrization"]
                            uv = _BA_API.model.reparametrize_on_surface(1, sample["curve"], [sample["parameter"]], face["tag"])
                            @test _ba_close(uv, sample["uv"])
                            @test _ba_close(_BA_API.model.get_value(2, face["tag"], uv), sample["value"])
                        end
                    end
                    @test _BA_API.model.get_boundary([(3,31)], false, true, false) == [Tuple(pair) for pair in row["boundary"]]
                end
            finally
                _BA_API.finalize()
            end
        end
    end
end

@testset "Pinned OCC Box two-corner reconstruction arithmetic" begin
    for (index, row) in enumerate(BOX_ENDPOINT_ROWS)
        @testset "endpoint row $index" begin
            _BA_API.initialize()
            try
                accepted = true
                try
                    _BA_API.model.add_box(row["origin"], 0, 0, row["width"], 1, 1; tag=31)
                catch error
                    println("BOX_ENDPOINT_ERROR row=", index, " ", sprint(showerror, error))
                    accepted = false
                end
                @test accepted
                if accepted
                    @test _BA_API.CURRENT[].box_extents[31] == (row["minimum"], 0.0, 0.0, row["actual_delta"], 1.0, 1.0)
                    for point in row["points"]
                        @test _BA_API.model.get_value(0, point["tag"], Float64[]) == point["value"]
                    end
                    for (dimension, field) in ((1, "edge_domains"), (2, "face_domains"))
                        for entity in row[field]
                            lower, upper = _BA_API.model.get_parametrization_bounds(dimension, entity["tag"])
                            @test lower == entity["bounds"][1]
                            @test upper == entity["bounds"][2]
                        end
                    end
                end
            finally
                _BA_API.finalize()
            end
        end
    end
end

end
