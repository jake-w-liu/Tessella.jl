using Test
using Random
using Tessella
import Tessella.MeshTypes as MT

# `_tet_topology` certifies a valid complex from global face and edge tables
# and leaves everything else to the vertex-link audit. The certificate must
# accept exactly the complexes the audit accepts: valid meshes, their damaged
# variants and small random complexes are run through both.

function ttc_agree(tets::Matrix{Int32}; pinch=Dict{Int32,Int}())
    fast = MT._tet_topology_accepts(tets, pinch)
    slow = first(MT._tet_topology_audit(tets; pinch_link_chi=pinch))
    return fast == slow, slow
end

# Kuhn subdivision of an n×n×n grid of cubes: six tetrahedra per cube.
function ttc_grid(n::Int)
    id(i, j, k) = Int32(i + (n + 1) * (j + (n + 1) * k) + 1)
    tets = Matrix{Int32}(undef, 4, 6n^3)
    column = 0
    for k in 0:n-1, j in 0:n-1, i in 0:n-1
        for (first, second, third) in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
            corner = [i, j, k]
            path = Int32[id(corner...)]
            for axis in (first, second, third)
                corner[axis] += 1
                push!(path, id(corner...))
            end
            tets[:, column+=1] = path
        end
    end
    return tets
end

function ttc_ring(n::Int)
    tets = Matrix{Int32}(undef, 4, n)
    for i in 1:n
        tets[:, i] = Int32[1, 2, 2 + i, 2 + mod1(i + 1, n)]
    end
    return tets
end

function ttc_mutate(rng, tets::Matrix{Int32})
    out = copy(tets); nt = size(out, 2); nv = Int(maximum(out))
    kind = rand(rng, 1:9)
    if kind == 1          # the same tetrahedron twice
        out = hcat(out, out[:, rand(rng, 1:nt)])
    elseif kind == 2      # missing tetrahedra
        keep = trues(nt)
        for _ in 1:rand(rng, 1:4); keep[rand(rng, 1:nt)] = false; end
        out = out[:, keep]
    elseif kind == 3      # one relabelled corner
        out[rand(rng, 1:4), rand(rng, 1:nt)] = Int32(rand(rng, 1:nv))
    elseif kind == 4      # two vertices identified
        a, b = Int32(rand(rng, 1:nv)), Int32(rand(rng, 1:nv))
        out[out .== a] .= b
    elseif kind == 5      # a tetrahedron touching at one vertex
        out = hcat(out, Int32[rand(rng, 1:nv), nv + 1, nv + 2, nv + 3])
    elseif kind == 6      # a tetrahedron touching along one edge
        t = rand(rng, 1:nt)
        out = hcat(out, Int32[out[1, t], out[2, t], nv + 1, nv + 2])
    elseif kind == 7      # a second component
        out = hcat(out, Int32[nv + 1, nv + 2, nv + 3, nv + 4])
    elseif kind == 8      # an arbitrary extra tetrahedron
        out = hcat(out, Int32.(randperm(rng, nv)[1:4]))
    else                  # relabelling of columns and slots keeps the complex
        out = out[:, randperm(rng, nt)]
        for t in 1:nt
            out[:, t] = out[randperm(rng, 4), t]
        end
    end
    return out
end

@testset "Tetrahedral topology certificate" begin
    rng = MersenneTwister(20261118)

    @testset "valid complexes are certified without the audit" begin
        for n in 1:4
            grid = ttc_grid(n)
            @test MT._tet_topology_accepts(grid, Dict{Int32,Int}())
            @test MT._tet_topology(grid) == (true, "ok")
        end
        for n in 3:9
            @test MT._tet_topology_accepts(ttc_ring(n), Dict{Int32,Int}())
            @test MT._tet_topology_accepts(ttc_ring(n)[:, 1:n-1], Dict{Int32,Int}())
        end
        @test MT._tet_topology_accepts(reshape(Int32[1, 2, 3, 4], 4, 1), Dict{Int32,Int}())
    end

    @testset "defects are left to the audit, which names them" begin
        single = reshape(Int32[1, 2, 3, 4], 4, 1)
        for (tets, reason) in (
                (hcat(single, Int32[2, 1, 3, 4]), "duplicate tetrahedron"),
                (hcat(single, Int32[1, 2, 3, 5], Int32[1, 2, 3, 6]), "non-manifold link edge"),
                (hcat(single, Int32[1, 5, 6, 7]), "disconnected link"),
                (reshape(Int32[1, 2, 3, 3], 4, 1), "repeats a vertex"))
            @test !MT._tet_topology_accepts(tets, Dict{Int32,Int}())
            ok, message = MT._tet_topology(tets)
            @test !ok
            @test occursin(reason, message)
        end
        # Labels the tables cannot index are not a verdict.
        @test !MT._tet_topology_accepts(reshape(Int32[0, 1, 2, 3], 4, 1), Dict{Int32,Int}())
        @test !MT._tet_topology_accepts(reshape(Int32[1, 2, 3, 10^6], 4, 1), Dict{Int32,Int}())
        @test MT._tet_topology(reshape(Int32[1, 2, 3, 10^6], 4, 1)) == (true, "ok")
    end

    @testset "certificate and audit agree on damaged meshes" begin
        valid = 0; invalid = 0
        for base in (ttc_grid(2), ttc_grid(3), ttc_ring(5), ttc_ring(8))
            for _ in 1:300
                mutated = ttc_mutate(rng, base)
                same, ok = ttc_agree(mutated)
                @test same
                ok ? (valid += 1) : (invalid += 1)
                same, ok = ttc_agree(ttc_mutate(rng, mutated))
                @test same
                ok ? (valid += 1) : (invalid += 1)
            end
        end
        @test valid > 100 && invalid > 100
    end

    @testset "certificate and audit agree on random small complexes" begin
        valid = 0
        for _ in 1:6000
            nv = rand(rng, 4:8); nt = rand(rng, 1:7)
            tets = Matrix{Int32}(undef, 4, nt)
            for t in 1:nt
                tets[:, t] = Int32.(randperm(rng, nv)[1:4])
            end
            same, ok = ttc_agree(tets)
            @test same
            ok && (valid += 1)
        end
        @test valid > 0
    end

    @testset "pinched boundary vertices use the caller's expected link" begin
        grid = ttc_grid(2)
        for _ in 1:100
            vertex = Int32(rand(rng, 1:Int(maximum(grid))))
            same, _ = ttc_agree(grid; pinch=Dict{Int32,Int}(vertex => rand(rng, (-1, 0, 1, 2))))
            @test same
        end
    end

    @testset "a vertex shared by every tetrahedron takes the large-bucket path" begin
        # A cone over a triangulated disc: every face of the fan starts at the apex.
        n = 60
        fan = Matrix{Int32}(undef, 4, n)
        for i in 1:n
            fan[:, i] = Int32[1, 2, 2 + i, 3 + i]
        end
        same, ok = ttc_agree(fan)
        @test same && ok
    end
end
