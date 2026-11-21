module GeoAttractorAnisoCurveTest

using Test, Tessella

# `Field[] = AttractorAnisoCurve` used to be rejected by both the .geo
# executor and the API session ("requires model context"). Upstream samples
# `Sampling` points per curve over its parameter range with their unit
# tangents (`AttractorAnisoCurveField::update`). Gmsh 4.15.2 meshes this
# square with Mesh.Algorithm 6 into 194 nodes: Frontal-Delaunay reads the
# field's scalar operator max(d, 0.05) (BGM_MeshSize), not its tensor.
const GEO = """
    Point(1)={0,0,0,0.1}; Point(2)={1,0,0,0.1}; Point(3)={1,1,0,0.1}; Point(4)={0,1,0,0.1};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    Point(5)={0.2,0.5,0}; Point(6)={0.8,0.5,0}; Line(5)={5,6};
    Field[1]=AttractorAnisoCurve; Field[1].CurvesList={5}; Field[1].DistMin=0.05; Field[1].DistMax=0.3;
    Field[1].SizeMinNormal=0.01; Field[1].SizeMaxNormal=0.1; Field[1].SizeMinTangent=0.1; Field[1].SizeMaxTangent=0.1;
    Background Field=1;
    """

@testset "AttractorAnisoCurve fields build from the model" begin
    path = joinpath(mktempdir(), "aniso.geo")
    write(path, GEO)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(2)
        nodes = length(A.mesh.get_nodes()[1])
        @test nodes > 100
        # Tessella's Delaunay/Frontal arms insert along the tensor and give
        # 154 nodes; Gmsh's scalar size gives 194.
        @test_broken nodes == 194
    finally
        A.finalize()
    end
end

end # module
