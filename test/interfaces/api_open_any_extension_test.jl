module APIOpenAnyExtensionTest

using Test, Tessella

const API = Tessella.API

# Gmsh's MergeFile reads a file whose extension it does not know by its first
# line: an MSH header opens as a mesh, a post-processing view as a view, and
# anything else is parsed as a .geo script (gmsh.open("script") meshes a
# script without extension; the model is named after the file).
@testset "open a file without a known extension as Gmsh does" begin
    mktempdir() do dir
        script = joinpath(dir, "script")
        write(script, "Point(1)={0,0,0,0.5};Point(2)={1,0,0,0.5};Line(1)={1,2};\n")
        mesh = joinpath(dir, "mesh.txt")
        write(mesh, """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Nodes
1 2 1 2
1 1 0 2
1
2
0 0 0
1 0 0
\$EndNodes
\$Elements
1 1 1 1
1 1 1 1
1 1 2
\$EndElements
""")
        API.initialize()
        try
            API.option("General.Terminal", 0)
            API.open(script)
            @test API.model.get_entities() == [(0, 1), (0, 2), (1, 1)]
            @test API.model.get_current() == "script"
            API.open(mesh)
            @test API.model.get_current() == "mesh"
            @test length(API.mesh.get_nodes()[1]) == 2
            @test API.mesh.get_elements(1)[2] == [UInt64[1]]
        finally
            API.finalize()
        end
    end
end

end
