# The suite includes every test file into Main; this file keeps its names to itself.
module BackgroundMeshViewTests
using Test
using Tessella
const A = Tessella.API

# `Background Mesh View[i]` makes a PostView field of view `i` and selects it
# as the background mesh (GModel::setBackgroundMesh), as Gmsh's tutorial t7
# does. Opened from a file, the script's selection was lost and the view sized
# nothing: t7 meshed its rectangle with two triangles where Gmsh makes 3,315.

const _VIEW = """
View "sizes" {
ST(0,0,0, 1,0,0, 1,1,0){0.05, 0.2, 0.2};
ST(0,0,0, 1,1,0, 0,1,0){0.05, 0.2, 0.05};
};
"""

_square(background, new_model) = """
Merge "bg.pos";
$(new_model ? "NewModel;" : "")
Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {1,1,0}; Point(4) = {0,1,0};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
$background
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
"""

function _mesh(script)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        mktempdir() do directory
            write(joinpath(directory, "bg.pos"), _VIEW)
            path = joinpath(directory, "case.geo"); write(path, script)
            A.open(path)
            fields = [(tag, A.mesh.field.get_type(tag)) for tag in A.mesh.field.list()]
            A.mesh.generate(1)
            curves = [length(A.mesh.get_nodes(1, curve)[1]) for curve in 1:4]
            A.mesh.generate(2)
            nodes = length(A.mesh.get_nodes()[1])
            triangles = sum(length, A.mesh.get_elements(2)[2]; init=0)
            return (; fields, curves, nodes, triangles)
        end
    finally
        A.finalize()
    end
end

@testset "A view selected as background mesh sizes the mesh of an opened file" begin
    for new_model in (false, true)
        result = _mesh(_square("Background Mesh View[0];", new_model))
        @test result.fields == [(1, "PostView")]
        # Gmsh 4.15.2: 9, 4, 9 and 19 nodes inside the four sides (0.05 along
        # x = 0, 0.2 along x = 1), 148 nodes and 249 triangles.
        @test result.curves == [9, 4, 9, 19]
        if Sys.ARCH === :x86_64
            @test (result.nodes, result.triangles) == (148, 249)
        else
            @test abs(result.nodes - 148) <= 3 && abs(result.triangles - 249) <= 6
        end
    end
    # Without the statement nothing sizes the square: its corners, one node in
    # the middle and four triangles, in Gmsh as here.
    plain = _mesh(_square("", false))
    @test isempty(plain.fields)
    @test plain.curves == [0, 0, 0, 0]
    @test (plain.nodes, plain.triangles) == (5, 4)
    # The same field named in the script's own words gives the same mesh.
    explicit = _mesh(_square("Field[1] = PostView; Field[1].ViewIndex = 0; Background Field = 1;", false))
    @test explicit.curves == [9, 4, 9, 19]
    Sys.ARCH === :x86_64 && @test (explicit.nodes, explicit.triangles) == (148, 249)
end
end
