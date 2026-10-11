@testset "current-model export without a view" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0);A.model.add("export-smoke")
        A.model.add_discrete_entity(1,7);A.mesh.add_nodes(1,7,[101,102],[0.,0,0,1,0,0])
        A.mesh.add_elements_by_type(7,1,[1001],[101,102])
        before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_tags())
        mktempdir() do directory
            path=joinpath(directory,"model.msh");A.write(path)
            mesh=Tessella.Elements.read_mixed_msh(path)
            @test mesh.entity_data.external_node_tags==[101,102] && only(mesh.blocks).nodes==reshape([1,2],2,1)
            @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_tags())==before
        end
    finally;A.finalize();end
end
