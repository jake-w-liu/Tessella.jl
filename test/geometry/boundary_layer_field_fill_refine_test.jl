module BoundaryLayerFieldFillRefineTest

using Test, Tessella

# A 2-D `BoundaryLayer Field` grows its strip and then the remaining domain is
# refined by the surface algorithm under the background size (upstream
# `meshGenerator` keeps inserting after `modifyInitialMeshForBoundaryLayers`).
# Before, the fill was a bare constrained triangulation of the boundary nodes:
# 118 nodes and triangles up to ~0.9 long for `Mesh.MeshSizeMax = 0.1`.
# Gmsh 4.15.2: 215 nodes, 232 triangles, 70 quadrangles.
@testset "BoundaryLayer field fill is refined" begin
    path = joinpath(mktempdir(), "bl.geo")
    write(path, """
        Point(1)={0,0,0,0.1}; Point(2)={1,0,0,0.1}; Point(3)={1,1,0,0.1}; Point(4)={0,1,0,0.1};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
        Mesh.MeshSizeExtendFromBoundary=0; Mesh.MeshSizeFromPoints=0; Mesh.MeshSizeFromCurvature=0;
        Field[1]=BoundaryLayer; Field[1].CurvesList={1}; Field[1].PointsList={1,2};
        Field[1].Size=0.005; Field[1].Ratio=1.3; Field[1].Thickness=0.1; Field[1].Quads=1;
        BoundaryLayer Field=1; Mesh.MeshSizeMax=0.1;
        """)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(2)
        tags, xyz, _ = A.mesh.get_nodes()
        types, _, nodes = A.mesh.get_elements(2)
        quads = length(nodes[findfirst(==(3), types)]) ÷ 4
        @test quads == 70
        @test 200 <= length(tags) <= 230
        # Every triangle edge respects the background size up to the usual
        # insertion slack.
        position = Dict(t => (xyz[3i-2], xyz[3i-1]) for (i, t) in enumerate(tags))
        tri = nodes[findfirst(==(2), types)]
        longest = 0.0
        for q in 1:3:length(tri), (a, b) in ((0, 1), (1, 2), (2, 0))
            p = position[tri[q+a]]; r = position[tri[q+b]]
            longest = max(longest, hypot(p[1] - r[1], p[2] - r[2]))
        end
        @test longest < 0.16
    finally
        A.finalize()
    end
end

end # module
