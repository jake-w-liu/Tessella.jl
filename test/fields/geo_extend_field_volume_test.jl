module GeoExtendFieldVolumeTest

using Test, Tessella

const API = Tessella.API

# Gmsh's Extend field with a SurfacesList (examples/boolean/extend_field.geo):
# while a volume is meshed, the size at a point is interpolated between the
# mean edge length of the nearest surface triangle (its barycentre) and
# SizeMax over DistMax (Field.cpp ExtendField::operator(), recomputeSurfaces).
# Tessella threw "SurfacesList is not supported". The geometry below is
# meshed by Gmsh 4.15.2 with 418 interior nodes (4,276 tetrahedra), and with
# 49 (1,856 tetrahedra) when SurfacesList is left out. Tessella's surface
# mesh already differs from Gmsh's, so the test checks the extension itself.
const GEO = """
SetFactory("OpenCASCADE");
Box(1) = {0,0,0, 1,1,1};
MeshSize{ : } = 0.3;
MeshSize{ PointsOf{ Surface{1}; } } = 0.08;
Mesh.MeshSizeExtendFromBoundary = 0;
Field[1] = Extend;
Field[1].SurfacesList = {Surface{:}};
Field[1].CurvesList = {Curve{:}};
Field[1].DistMax = 0.5;
Field[1].SizeMax = 0.1;
Field[1].Power = 2;
Background Field = 1;
"""

count_cells(dim) = sum(length, API.mesh.get_elements(dim)[2]; init=0)

@testset "Extend field sizes volumes from the surface triangles" begin
    mktempdir() do dir
        path = joinpath(dir, "extend.geo")
        write(path, GEO)
        interior = Dict{Bool,Int}()
        for with_surfaces in (true, false)
            write(path, with_surfaces ? GEO :
                  replace(GEO, "Field[1].SurfacesList = {Surface{:}};\n" => ""))
            API.initialize()
            try
                API.option("General.Terminal", 0)
                API.open(path)
                API.mesh.generate(3)
                @test count_cells(3) > 0
                interior[with_surfaces] = length(API.mesh.get_nodes(3, 1)[1])
            finally
                API.finalize()
            end
        end
        # Without SurfacesList the field returns MAX_LC in the volume; with
        # it, SizeMax = 0.1 and the fine triangles refine the whole volume.
        @test interior[true] > 4 * interior[false]
    end
end

end
