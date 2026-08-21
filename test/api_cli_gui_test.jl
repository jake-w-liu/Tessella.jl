using Test
using Tessella
using Tessella.API
using Tessella.CLI: main
using Tessella.GUI: GuiState, gui_command!
using Tessella.Post: View, view_value, apply_plugin
using Tessella.MeshTypes: ntets, nnodes, Mesh, tet_volume, node
using Tessella.Geometry: box_surface

@testset "API, CLI, GUI, and post views" begin
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

    surf=box_surface(0,1,0,1,0,1)
    values=collect(Float64, 1:nnodes(surf))
    v=View("n",surf,values)
    @test view_value(v,1)==1.0
    absv=apply_plugin("Abs",v)
    @test view_value(absv,1)==1.0
    @test apply_plugin("IsosurfaceZeroCount",v)==0
    @test_throws ArgumentError apply_plugin("NoPlugin",v)
end
