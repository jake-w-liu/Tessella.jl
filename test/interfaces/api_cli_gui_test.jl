using Test
using Tessella
using Tessella.API
using Tessella.CLI: main
using Tessella.GUI: GuiState, gui_command!
using Tessella.MeshTypes: ntets, tet_volume, node

@testset "API, CLI, and GUI" begin
    API.initialize()
    tag=API.model.add_box(0,0,0,1,1,1; tag=1)
    @test tag==1
    @test API.option("Mesh.MeshSizeFactor")==1.0
    API.option("Mesh.MeshSizeFactor", 2.0)
    @test API.option("Mesh.MeshSizeFactor")==2.0
    mesh=API.mesh.generate(3)
    @test ntets(mesh)>0
    @test_throws ArgumentError API.option("No.Such.Option")
    API.finalize()
    @test_throws ArgumentError API.mesh.generate(3)

    API.initialize()
    API.model.add_box(0,0,0,1,1,1; tag=1)
    API.model.add_box(2,0,0,1,1,1; tag=2)
    @test_throws ArgumentError API.mesh.generate(3)
    API.finalize()

    API.initialize()
    API.model.add_box(0,0,0,2,1,1; tag=1)
    API.model.add_box(0,0,0,1,1,1; tag=2)
    @test API.model.boolean_difference(1,2; tag=3)==3
    cut=API.mesh.generate(3)
    @test ntets(cut)>0
    CV=sum(tet_volume(node(cut,cut.tets[1,t]),node(cut,cut.tets[2,t]),
                      node(cut,cut.tets[3,t]),node(cut,cut.tets[4,t]))
           for t in 1:ntets(cut))
    @test CV≈1.0 atol=1e-12
    API.finalize()

    geo=mktemp() do path,io
        write(io,"Box(1) = {0, 0, 0, 1, 1, 1};\n")
        close(io)
        dest=path*".msh"
        main([path,"-3","-o",dest])
        @test isfile(dest)
        dest
    end

    s=GuiState()
    gui_command!(s,"new")
    gui_command!(s,"box 0 0 0 1 1 1 1")
    gui_command!(s,"mesh 3")
    @test s.mesh_tets>0
    gui_command!(s,"select 1 2")
    @test s.selected==[1,2]
    gui_command!(s,"camera 1 2 3")
    @test s.camera==(1.0,2.0,3.0)
    @test_throws ArgumentError gui_command!(s,"explode")

end
