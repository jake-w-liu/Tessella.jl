module GeoCreateTopologyGmshTest

using Test, Tessella

const API = Tessella.API

# mesh.create_topology is Gmsh 4.15.2's createTopologyFromMesh: the boundary
# of a region (face, curve) is looked up among the elements already in the
# model and only the rest becomes new entities, one per group of adjacent
# entities. Every expected value below was read from the Gmsh Python wheel.
# Gmsh orders new points and the lines of a split curve by pointer; those are
# compared as sets.

const CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
const TETS = [(1, 2, 3, 7), (1, 3, 4, 7), (1, 4, 8, 7), (1, 8, 5, 7), (1, 5, 6, 7), (1, 6, 2, 7)]

function cube(with_surface)
    API.model.add_discrete_entity(3, 1)
    xyz = Float64[]
    for p in CUBE
        append!(xyz, p)
    end
    connectivity = Int[]
    for t in TETS
        append!(connectivity, t)
    end
    API.mesh.add_nodes(3, 1, collect(1:8), xyz)
    API.mesh.add_elements(3, 1, [4], [collect(1:6)], [connectivity])
    if with_surface
        API.model.add_discrete_entity(2, 1)
        tris = [1, 4, 3, 1, 3, 2, 7, 8, 5, 7, 5, 6, 1, 2, 6, 1, 6, 5,
                7, 6, 2, 7, 2, 3, 7, 3, 4, 7, 4, 8, 1, 5, 8, 1, 8, 4]
        API.mesh.add_elements(2, 1, [2], [collect(10:21)], [tris])
    end
end

function plates()
    API.model.add_discrete_entity(2, 1)
    API.model.add_discrete_entity(2, 2)
    n = 2
    tags = Int[]
    xyz = Float64[]
    for j in 0:n, i in 0:2n
        push!(tags, 1 + i + (2n + 1) * j)
        append!(xyz, (i / n, j / n, 0.0))
    end
    API.mesh.add_nodes(2, 1, tags, xyz)
    for s in 1:2
        conn = Int[]
        for j in 0:n-1, i in 0:n-1
            a = 1 + (s - 1) * n + i + (2n + 1) * j
            append!(conn, (a, a + 1, a + 2n + 2, a, a + 2n + 2, a + 2n + 1))
        end
        API.mesh.add_elements(2, s, [2], [collect(100s:100s+7)], [conn])
    end
end

elements(d, t) = API.mesh.get_elements(d, t)

@testset "a region of tetrahedra gets one closed surface" begin
    try
        API.initialize(); API.option("General.Terminal", 0)
        cube(false)
        API.mesh.create_topology()
        @test API.model.get_entities() == [(2, 1), (3, 1)]
        types, tags, nodes = elements(2, 1)
        @test types == [2] && tags[1] == collect(7:18)
        @test nodes[1] == [1, 3, 2, 1, 2, 6, 1, 4, 3, 1, 8, 4, 1, 6, 5, 1, 5, 8,
                           7, 2, 3, 7, 6, 2, 7, 3, 4, 7, 4, 8, 7, 5, 6, 7, 8, 5]
        @test API.mesh.get_nodes(2, 1)[1] == collect(1:8)
        @test isempty(API.mesh.get_nodes(3, 1)[1])
        @test API.model.get_boundary([(3, 1)], false, true) == [(2, 1)]
    finally
        API.finalize()
    end
end

@testset "existing triangles bound the region; no new surface" begin
    try
        API.initialize(); API.option("General.Terminal", 0)
        cube(true)
        API.mesh.create_topology()
        @test API.model.get_entities() == [(2, 1), (3, 1)]
        @test elements(2, 1)[2][1] == collect(10:21)
        @test API.model.get_boundary([(3, 1)], false, true) == [(2, 1)]
        @test API.mesh.get_nodes(2, 1)[1] == collect(1:8)
    finally
        API.finalize()
    end
end

@testset "two plates: a curve per group of adjacent surfaces" begin
    try
        API.initialize(); API.option("General.Terminal", 0)
        plates()
        API.mesh.create_topology()
        @test API.model.get_entities() == [(0, 1), (0, 2), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2)]
        @test Set(elements(1, 1)[2][1]) == Set([208, 209, 210, 215, 218, 219])
        @test Set(elements(1, 2)[2][1]) == Set([211, 213, 214, 217, 220, 221])
        @test Set(elements(1, 3)[2][1]) == Set([212, 216])
        @test Set(vcat(elements(0, 1)[2][1], elements(0, 2)[2][1])) == Set([222, 223])
        @test Set(vcat(API.mesh.get_nodes(0, 1)[1], API.mesh.get_nodes(0, 2)[1])) == Set([3, 13])
        @test sort(API.mesh.get_nodes(1, 3)[1]) == [8]
        @test API.mesh.get_nodes(2, 1)[1] == [7]
        @test sort(API.model.get_boundary([(2, 1)], false, true)) == [(1, 1), (1, 3)]
        @test sort(API.model.get_boundary([(2, 2)], false, true)) == [(1, 2), (1, 3)]
        for c in 1:3
            @test sort(API.model.get_boundary([(1, c)], false, true)) == [(0, 1), (0, 2)]
        end
    finally
        API.finalize()
    end
end

end
