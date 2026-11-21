module GeoAttractorAnisoCurveTest

using Test, Tessella

# `Field[] = AttractorAnisoCurve` used to be rejected by both the .geo
# executor and the API session ("requires model context"). Upstream samples
# `Sampling` points per curve over its parameter range with their unit
# tangents (`AttractorAnisoCurveField::update`). Gmsh 4.15.2 meshes this
# square with Mesh.Algorithm 6 into 194 nodes: Frontal-Delaunay reads the
# field's scalar operator max(d, 0.05) (BGM_MeshSize), not its tensor; with
# Mesh.Algorithm 5 Gmsh gives 198 nodes.
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
        tags, xyz, _ = A.mesh.get_nodes()
        @test length(tags) == 194                    # Gmsh 4.15.2: 194
        # Mean edge length and triangle count per 0.1-wide band of distance
        # from the curve, as Gmsh 4.15.2 gives them.
        pos = Dict(tags[i] => (xyz[3i-2], xyz[3i-1]) for i in eachindex(tags))
        tri = A.mesh.get_elements(2)[3][1]
        dist(p) = (x = clamp(p[1], 0.2, 0.8); hypot(p[1] - x, p[2] - 0.5))
        sums = zeros(5); counts = zeros(Int, 5)
        for k in 1:3:length(tri)
            p = [pos[tri[k+j]] for j in 0:2]
            band = min(5, 1 + floor(Int, dist(((p[1][1]+p[2][1]+p[3][1])/3,
                                               (p[1][2]+p[2][2]+p[3][2])/3)) / 0.1))
            sums[band] += sum(hypot((p[j] .- p[mod1(j+1, 3)])...) for j in 1:3) / 3
            counts[band] += 1
        end
        @test counts == [112, 59, 48, 51, 50]
        @test round.(sums ./ counts, digits=3) == [0.058, 0.093, 0.1, 0.101, 0.095]
    finally
        A.finalize()
    end
end

end # module
