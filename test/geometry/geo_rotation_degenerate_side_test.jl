module GeoRotationDegenerateSideTest
# A rotation Extrude of a triangle with an edge on the axis
# (examples/simple_geo/cone.geo, one quarter, lc=0.5). Gmsh 4.15.2 makes the
# laterals three-curve MSH_SURF_TRIC patches and, with Geometry.ExactExtrusion
# (default 1), lifts them through the extrusion formula over the whole unit
# square: gmshEdge::reparamOnFace puts the third generatrix at V=1
# (`V = hack ? 1 : U`), so the u=0 side collapses onto the axis point.
# Tessella placed that generatrix on V=U (a triangle) while lifting through
# the square: the cone's lateral had 142 triangles against Gmsh's 160 at
# lc=0.2. Values below measured with the gmsh 4.15.2 wheel.
using Test, Tessella

const _GEO = """
lc = 0.5;
Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(6) = {0,0,3,lc};
Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,1};
Curve Loop(4) = {2,3,1};
Plane Surface(5) = {4};
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{5}; }
"""

@testset "rotation Extrude with a generatrix on the axis" begin
    api = Tessella.API
    api.initialize()
    try
        mktempdir() do dir
            path = joinpath(dir, "cone_quarter.geo")
            write(path, _GEO)
            api.open(path)
            # Gmsh: getType and oriented boundaries of the swept entities.
            @test [t for (_, t) in api.model.get_entities(2)] == [5, 12, 16, 17]
            @test api.model.get_type(2, 12) == "Surface"
            @test api.model.get_type(2, 16) == "Surface"
            @test [b[2] for b in api.model.get_boundary([(2, 12)], false, true)] == [-7, -11, 2]
            @test [b[2] for b in api.model.get_boundary([(2, 16)], false, true)] == [1, 11, -9]
            api.mesh.generate(2)
            count(t) = sum(length, api.mesh.get_elements(2, t)[2])
            @test count(5) == 17 && count(17) == 17
            # The cone lateral: Gmsh 34 triangles.
            @test count(12) == 34
            # The source Line 2, the patch's third generatrix, lies on v=1.
            line_nodes = Set(api.mesh.get_nodes(1, 2, false)[1])
            tags, _, uv = api.mesh.get_nodes(2, 12, true, true)
            on_line = [uv[2i] for (i, tag) in enumerate(tags) if tag in line_nodes]
            @test length(on_line) == length(line_nodes) > 0
            @test all(v -> abs(v - 1) <= 1e-9, on_line)
            # The flat sector: Gmsh 10 triangles, Tessella 8. Same boundary
            # parametrization and lift (Gmsh's getValue agrees), the 2-D
            # frontal mesher places one interior node where Gmsh places two.
            @test_broken count(16) == 10
        end
    finally
        api.finalize()
    end
end

end
