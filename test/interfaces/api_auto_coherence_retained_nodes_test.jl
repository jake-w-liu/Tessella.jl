module AutoCoherenceRetainedNodesTests
using Test, Tessella, SHA, TOML
const API=Tessella.API
const FIXTURES=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1")
bits(values)=[bytes2hex(reinterpret(UInt8,[Float64(value)])) for value in values]

@testset "AutoCoherence0 retains original public B4 CAD and actual nodes" begin
    @test bytes2hex(sha256(read(joinpath(FIXTURES,"native_retained_recombined.json")))) ==
        "aa658c1de4ba4523522b0a38b967d95f417d721ab8e7564b0af6a65e83d48c40"
    @test bytes2hex(sha256(read(joinpath(FIXTURES,"native_retained_recombined.toml")))) ==
        "56e9823e8905d1f48edff0b1a0f646a9a463babbcf7e35a2397919c4ec25164b"
    native=TOML.parsefile(joinpath(FIXTURES,"native_retained_recombined.toml"))
    @test native["gmsh_version"] == "4.15.2"
    @test native["generation_order"] == [3]
    @test bytes2hex(sha256(codeunits(native["source"]))) == native["source_sha256"]
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"original_b4.geo");write(path,native["source"])
            execution=API.open_geo!(path;mesh_dim=0)
            @test length(execution.model.points) == 44
            @test length(execution.model.curves) == 20
            API.mesh.generate(3)
            entities=API.model.get_entities()
            @test entities == [(Int(row["dimension"]),Int(row["tag"])) for row in native["entities"]]
            ids,xyz,uv=API.mesh.get_nodes()
            @test length(ids) == 8052
            @test ids == UInt64.(native["global_node_tags"])
            @test bits(xyz) == native["global_xyz_bits"]
            @test bits(uv) == native["global_uv_bits"]
            @test API.mesh.get_max_node_tag() == native["max_node_tag"]
            for row in native["entities"]
                dimension,tag=Int(row["dimension"]),Int(row["tag"])
                @testset "actual entity $dimension:$tag" begin
                    ids,xyz,uv=API.mesh.get_nodes(dimension,tag,false,true)
                    @test ids == UInt64.(row["node_tags"])
                    @test bits(xyz) == row["xyz_bits"]
                    @test bits(uv) == row["uv_bits"]
                    if dimension<=1
                        # Pointer ties can change native B4 face/volume splits.
                        # Its retained Point/Curve cells and all native node
                        # identities are independent of those higher choices.
                        types,labels,cells=API.mesh.get_elements(dimension,tag)
                        @test types == row["element_types"]
                        @test labels == [UInt64.(value) for value in row["element_tags"]]
                        @test cells == [UInt64.(value) for value in row["connectivity"]]
                    end
                end
            end
            setter=native["setter_immediate"]
            @test setter["all_records_and_history_unchanged"]
            before=deepcopy((API.mesh.get_nodes(),API.mesh.get_elements(),
                    API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag()))
            API.mesh.set_size([(0,Int(setter["point"]))],setter["size"])
            @test isequal(before,(API.mesh.get_nodes(),API.mesh.get_elements(),
                API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag()))
        end
    finally
        API.finalize()
    end
end
end
