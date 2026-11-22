module GeoDiscreteSTLVolumeTest

using Test, Tessella

const API = Tessella.API

# examples/simple_geo/sphere-discrete.geo: a Volume bounded by a merged,
# closed STL surface, which has no boundary curves. Gmsh 4.15.2 keeps the STL
# triangles through generate(1) and generate(2) and fills the volume in
# generate(3); the counts below are Gmsh's for this cube (each face split into
# four triangles about its centre).
function cube_stl(path)
    corners = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    faces = ((1, 2, 4, 3), (5, 7, 8, 6), (1, 5, 6, 2), (3, 4, 8, 7), (1, 3, 7, 5), (2, 6, 8, 4))
    open(path, "w") do io
        println(io, "solid cube")
        for f in faces
            q = [corners[i] for i in f]
            c = ntuple(k -> sum(v[k] for v in q) / 4, 3)
            for i in 1:4
                println(io, " facet normal 0 0 0\n  outer loop")
                for p in (q[i], q[mod1(i + 1, 4)], c)
                    println(io, "   vertex ", p[1], " ", p[2], " ", p[3])
                end
                println(io, "  endloop\n endfacet")
            end
        end
        println(io, "endsolid cube")
    end
end

count_cells(dim) = (b = API.mesh.get_elements(dim)[2]; isempty(b) ? 0 : sum(length, b))

@testset "Volume bounded by a merged closed STL surface" begin
    mktempdir() do dir
        cube_stl(joinpath(dir, "cube.stl"))
        write(joinpath(dir, "cube.geo"),
              "Merge \"cube.stl\";\nSurface Loop(1) = {1};\nVolume(1) = {1};\n")
        API.initialize()
        try
            API.option("General.Terminal", 0)
            API.open(joinpath(dir, "cube.geo"))
            @test API.model.get_entities() == [(2, 1), (3, 1)]
            @test collect(API.model.get_bounding_box(3, 1)) ≈ [0, 0, 0, 1, 1, 1]
            for dim in 1:2
                API.mesh.generate(dim)
                @test length(API.mesh.get_nodes()[1]) == 14
                @test (count_cells(2), count_cells(3)) == (24, 0)
            end
            API.mesh.generate(3)
            @test length(API.mesh.get_nodes()[1]) == 14
            @test (count_cells(2), count_cells(3)) == (24, 24)
        finally
            API.finalize()
        end
    end
end


# Gmsh's readSTL parses with sscanf("%lf"); its Windows build rounds through
# the 80-bit long double, so "-0.529113" (sphere-surf.stl, node 16) reads as
# -0.5291129999999999 there and as -0.529113 elsewhere.
@testset "STL coordinates parse like Gmsh's sscanf" begin
    mktempdir() do dir
        path = joinpath(dir, "one.stl")
        write(path, """solid one
 facet normal 0 0 1
  outer loop
   vertex -0.529113 0.222735 0.818799
   vertex 1 0 0
   vertex 0 1 0
  endloop
 endfacet
endsolid one
""")
        API.initialize()
        try
            API.option("General.Terminal", 0)
            API.open(path)
            x = API.mesh.get_nodes()[2][1]
            @test x == (Sys.iswindows() ? -0.5291129999999999 : -0.529113)
        finally
            API.finalize()
        end
    end
end

end
