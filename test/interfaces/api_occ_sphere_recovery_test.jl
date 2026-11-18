# The suite includes every test file into Main; this file keeps its aliases to itself.
module SphereRecoveryTests
using Test
using Tessella

# The vertices of a sphere mesh are cospherical, so their Delaunay
# tetrahedralization is decided by rounding and misses surface edges at some
# sizes. A missing edge that runs under a hull edge sits in a thin twisted
# cavity no refill can fill; it is recovered once a vertex is seeded on the
# inner side. Recovery used to seed one side only, could pick the outer one,
# and then left the whole fill to the slow float recovery.
#
# The fixture is the unit sphere's surface mesh at size 0.2 as it was when
# this was found (405 vertices, 806 triangles, two missing edges): reading it
# keeps the case independent of the surface mesher and of its cost.

function sphere_fixture()
    lines = readlines(joinpath(@__DIR__, "..", "fixtures", "sphere_surface_h02.txt"))
    nodes, triangles = parse.(Int, split(lines[1]))
    coords = Matrix{Float64}(undef, 3, nodes)
    for j in 1:nodes
        coords[:, j] = parse.(Float64, split(lines[1+j]))
    end
    tris = Matrix{Int32}(undef, 3, triangles)
    for j in 1:triangles
        tris[:, j] = parse.(Int32, split(lines[1+nodes+j]))
    end
    return coords, tris
end

@testset "A sphere surface is recovered without the float fallback" begin
    coords, tris = sphere_fixture()
    @test size(coords, 2) == 405 && size(tris, 2) == 806
    fill, reason, stats = Tessella.BoundaryRecovery.recover_boundary_flips(coords, tris)
    @test fill !== nothing
    if fill !== nothing
        @test stats.missing_edges == 2
        @test validate(fill).ok
        @test Tessella.Mesh3D._fill_boundary_is_surface(Mesh(coords; tris=tris), fill)
        # Added vertices are interior: the surface rows are kept as given.
        @test fill.coords[:, 1:size(coords, 2)] == coords
        @test size(fill.coords, 2) - size(coords, 2) <= 4
    else
        @info "sphere recovery" reason stats
    end
end
end
