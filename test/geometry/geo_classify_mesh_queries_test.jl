module GeoClassifyMeshQueriesTest

using Test, Tessella

const API = Tessella.API

# classify_surfaces on a merged STL (a mesh that is already cached) and the
# mesh queries right after it, before create_geometry. Gmsh 4.15.2 on this 4x4 plate: one surface (2), one
# closed curve (3) of 16 lines whose boundary lists its point twice, one point
# (1); 25 nodes, 32 triangles. Tessella threw a BoundsError in get_elements
# (the cache kept the pre-classification ownership) and listed the point once.
@testset "mesh queries after classify_surfaces, as in Gmsh" begin
    try
        API.initialize()
        API.option("General.Terminal", 0)
        n = 4
        point(i, j) = (i / n, j / n, 0.0)
        path = joinpath(mktempdir(), "plate.stl")
        open(path, "w") do io
            println(io, "solid plate")
            for j in 0:n-1, i in 0:n-1
                for t in ((point(i, j), point(i + 1, j), point(i + 1, j + 1)),
                          (point(i, j), point(i + 1, j + 1), point(i, j + 1)))
                    println(io, " facet normal 0 0 1"); println(io, "  outer loop")
                    foreach(v -> println(io, "   vertex ", join(v, " ")), t)
                    println(io, "  endloop"); println(io, " endfacet")
                end
            end
            println(io, "endsolid plate")
        end
        API.merge(path)
        API.mesh.classify_surfaces(pi / 4, true, false, pi)
        types, tags, _ = API.mesh.get_elements()
        @test collect(types) == [1, 2, 15]
        @test length.(tags) == [16, 32, 1]
        @test length(API.mesh.get_nodes()[1]) == 25
        @test API.model.get_entities() == [(0, 1), (1, 3), (2, 2)]
        @test API.model.get_boundary([(1, 3)], false, true) == [(0, 1), (0, 1)]
        @test API.model.get_boundary([(2, 2)], false, true) == [(1, 3)]
    finally
        API.finalize()
    end
end

end
