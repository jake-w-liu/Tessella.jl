using Test, Tessella, SHA, TOML

const _AC_FIXTURE_DIRECTORY = joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1")
const _AC_INPUTS = (
    ("native_geometry_controls_v3", "991451ac2885f4f647b93abb65f4ff3b79948d35ad48b121c79ceac3a61e5b94"),)

@testset "CAD automatic coherence" begin
    for (stem,toml_hash) in _AC_INPUTS
        @test bytes2hex(sha256(read(joinpath(_AC_FIXTURE_DIRECTORY,stem*".toml")))) == toml_hash
        native=TOML.parsefile(joinpath(_AC_FIXTURE_DIRECTORY,stem*".toml"))
        @test native["gmsh_version"] == "4.15.2"
        @test length(native["cases"]) == 23
        for case in native["cases"]
            @testset "$(case["name"])" begin
                source=case["source"]
                @test bytes2hex(sha256(codeunits(source))) == case["source_sha256"]
                execution=mktemp() do path,io
                    write(io,source);close(io)
                    execute_geo(path;mesh_dim=0)
                end
                model=execution.model
                entities=Tessella.Model.model_entities(model)
                @test entities == [(Int(dim),Int(tag)) for (dim,tag) in zip(case["dimensions"],case["tags"])]
                @test [Tessella.Model.model_entity_type(model,dim,tag) for (dim,tag) in entities] == case["types"]
                incidence=deepcopy((model.points,model.curves,model.loops,model.surfaces,
                                    model.surface_loops,model.volumes,model.next_tag))
                for ((dim,tag),boundary) in zip(entities,case["boundaries"])
                    @test Tessella.Model.model_boundary(model,[(dim,tag)],false,true,false) ==
                        [(Int(item[1]),Int(item[2])) for item in boundary]
                end
                @test isequal(incidence,(model.points,model.curves,model.loops,model.surfaces,
                                        model.surface_loops,model.volumes,model.next_tag))
                @test sort!(collect(keys(model.points))) == case["point_tags"]
                # The pinned rotation residual is the Windows native x87
                # observation; other platforms use their own native libm.
                if Sys.iswindows() || !occursin("Pi/",source)
                    for (tag,expected) in zip(case["point_tags"],case["point_xyz_bits"])
                        xyz=Tuple(reinterpret(Float64,hex2bytes(word))[1] for word in expected)
                        @test isequal(model.points[tag],xyz)
                    end
                end
                for (name,values) in zip(case["list_names"],case["list_values"])
                    @test isequal(execution.lists[name],values)
                end
                if case["name"] in ("original_b4_free_mode0","original_b4_recombined_mode0")
                    @test length(model.points) == 44
                    @test length(model.curves) == 20
                    @test Tessella.Model._model_volume_boundary_surfaces(
                        model,only(keys(model.volumes)),"automatic coherence test")==[-1,26,13,17,21,25]
                end
            end
        end
    end
end
