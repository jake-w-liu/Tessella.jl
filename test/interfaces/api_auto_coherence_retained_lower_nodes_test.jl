module AutoCoherenceRetainedLowerNodesTests
using Test, Tessella, SHA, TOML
const API=Tessella.API
const FIXTURES=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1")
bits(values)=[bytes2hex(reinterpret(UInt8,[Float64(value)])) for value in values]

@testset "AutoCoherence0 actual lower producers retain 8052 native nodes" begin
    @test bytes2hex(sha256(read(joinpath(FIXTURES,"native_retained_recombined.json")))) ==
        "aa658c1de4ba4523522b0a38b967d95f417d721ab8e7564b0af6a65e83d48c40"
    @test bytes2hex(sha256(read(joinpath(FIXTURES,"native_retained_lower_phase_v1.json")))) ==
        "49f8cae96d75e15186e84ba65929e99cffe2478db2dae83c96bfd286d7780e80"
    native=TOML.parsefile(joinpath(FIXTURES,"native_retained_recombined.toml"))
    @test bytes2hex(sha256(codeunits(native["source"]))) == native["source_sha256"]
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"original_b4.geo");write(path,native["source"])
            execution=API.open_geo!(path;mesh_dim=0)
            model=execution.model
            @test length(model.points)==44 && length(model.curves)==20
            volume=only(keys(model.volumes))
            literal=Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"literal retained boundary";certify_shells=false)
            @test_throws ArgumentError Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"strict retained CAD shell")
            cad=deepcopy((model.points,model.curves,model.loops,model.surfaces,
                model.surface_loops,model.volumes,model.meshing.extrude_sources,model.next_tag))
            API.mesh.generate(2)
            model=API.CURRENT[]
            @test isequal(cad,(model.points,model.curves,model.loops,model.surfaces,
                model.surface_loops,model.volumes,model.meshing.extrude_sources,model.next_tag))
            @test Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"literal retained boundary";certify_shells=false)==literal
            @test_throws ArgumentError Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"strict retained CAD shell after generation")
            @test API.model.get_entities()==[(Int(row["dimension"]),Int(row["tag"])) for row in native["entities"]]
            ids,xyz,uv=API.mesh.get_nodes()
            @test length(ids)==8052
            @test ids==UInt64.(native["global_node_tags"])
            @test bits(xyz)==native["global_xyz_bits"]
            @test bits(uv)==native["global_uv_bits"]
            @test API.mesh.get_max_node_tag()==native["max_node_tag"]
            for row in native["entities"]
                dimension,tag=Int(row["dimension"]),Int(row["tag"])
                @testset "actual lower entity $dimension:$tag" begin
                    ids,xyz,uv=API.mesh.get_nodes(dimension,tag,false,true)
                    @test ids==UInt64.(row["node_tags"])
                    @test bits(xyz)==row["xyz_bits"]
                    @test bits(uv)==row["uv_bits"]
                    if dimension<=1
                        types,labels,cells=API.mesh.get_elements(dimension,tag)
                        @test types==row["element_types"]
                        @test labels==[UInt64.(value) for value in row["element_tags"]]
                        @test cells==[UInt64.(value) for value in row["connectivity"]]
                    end
                end
            end
        end
    finally
        API.finalize()
    end
end
end
