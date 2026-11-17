module AutoCoherenceRetainedNodesTests
using Test, Tessella, TOML
const API=Tessella.API
const FIXTURE=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1",
    "native_geometry_controls_v3.toml")

@testset "AutoCoherence0 preserves B4 nodes and mesh state" begin
    source=only(filter(row->row["name"]=="original_b4_recombined_mode0",
        TOML.parsefile(FIXTURE)["cases"]))["source"]
    source=replace(source,"Layers{500}"=>"Layers{2}")
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"original_b4.geo");write(path,source)
            execution=API.open_geo!(path;mesh_dim=0)
            @test length(execution.model.points) == 44
            @test length(execution.model.curves) == 20
            API.mesh.generate(3)
            ids,xyz,_=API.mesh.get_nodes()
            @test length(ids) == 84
            @test allunique(ids) && all(isfinite,xyz)
            before=deepcopy((API.mesh.get_nodes(),API.mesh.get_elements(),
                    API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag()))
            API.mesh.set_size([(0,1)],0.25)
            @test isequal(before,(API.mesh.get_nodes(),API.mesh.get_elements(),
                API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag()))
        end
    finally
        API.finalize()
    end
end
end
