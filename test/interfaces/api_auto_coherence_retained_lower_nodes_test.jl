module AutoCoherenceRetainedLowerNodesTests
using Test, Tessella, TOML
const API=Tessella.API
const FIXTURE=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1",
    "native_geometry_controls_v3.toml")

@testset "AutoCoherence0 preserves CAD and surface nodes" begin
    source=only(filter(row->row["name"]=="original_b4_recombined_mode0",
        TOML.parsefile(FIXTURE)["cases"]))["source"]
    source=replace(source,"Layers{500}"=>"Layers{2}")
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"original_b4.geo");write(path,source)
            execution=API.open_geo!(path;mesh_dim=0)
            model=execution.model
            @test length(model.points)==44 && length(model.curves)==20
            volume=only(keys(model.volumes))
            literal=Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"literal retained boundary";certify_shells=false)
            @test Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"certified retained operation")==literal
            cad=deepcopy((model.points,model.curves,model.loops,model.surfaces,
                model.surface_loops,model.volumes,model.meshing.extrude_sources,model.next_tag))
            API.mesh.generate(2)
            model=API.CURRENT[]
            @test isequal(cad,(model.points,model.curves,model.loops,model.surfaces,
                model.surface_loops,model.volumes,model.meshing.extrude_sources,model.next_tag))
            @test Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"literal retained boundary";certify_shells=false)==literal
            @test Tessella.Model._model_volume_boundary_surfaces(
                model,volume,"certified retained operation after generation")==literal
            ids,xyz,_=API.mesh.get_nodes()
            @test length(ids)==84
            @test allunique(ids) && all(isfinite,xyz)
        end
    finally
        API.finalize()
    end
end
end
