module ExtrudeFullTurnUnstructuredTest

using Test, Tessella

# Three 2*Pi/3 revolves of an unstructured face close onto the face itself;
# the closing layer must weld onto the face's actual mesh (upstream
# `pos.find`), not onto a fresh remesh of it. Gmsh 4.15.2: 122 nodes,
# 156 prisms, 78 triangles, 72 quadrangles, 60 lines, 14 points.
@testset "full-turn sweep of an unstructured face welds its seam" begin
    path = joinpath(mktempdir(), "revolve.geo")
    write(path, """
        Point(1)={1,0,0,0.4}; Point(2)={2,0,0,0.4}; Point(3)={2,0,1,0.4}; Point(4)={1,0,1,0.4};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
        s1[] = Extrude {{0,0,1},{0,0,0},2*Pi/3} { Surface{1}; Layers{2}; Recombine; };
        s2[] = Extrude {{0,0,1},{0,0,0},2*Pi/3} { Surface{s1[0]}; Layers{2}; Recombine; };
        s3[] = Extrude {{0,0,1},{0,0,0},2*Pi/3} { Surface{s2[0]}; Layers{2}; Recombine; };
        """)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        # Dimension by dimension: the volume pass reads the retained face
        # meshes. (A direct `generate(3)` also needs the higher-dimension
        # certificate to read the same retained product.)
        for dim in 1:3
            A.mesh.generate(dim)
        end
        @test length(A.mesh.get_nodes()[1]) == 122
        types, tags, _ = A.mesh.get_elements()
        @test Dict(t => length(e) for (t, e) in zip(types, tags)) ==
              Dict(1 => 60, 2 => 78, 3 => 72, 6 => 156, 15 => 14)
    finally
        A.finalize()
    end
end

end # module
