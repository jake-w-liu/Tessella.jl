using Test
using Tessella
import Tessella.MeshTypes as MT

# Boundary recovery fills a closed triangulated surface with tetrahedra whose
# boundary faces are exactly the given triangles. The surfaces below force
# every stage: none (convex, Delaunay-conforming), flips (a torus whose quads
# take the other diagonal), cavity refills (exactly coplanar grids), nested
# shells, and inputs that are not a valid boundary.

const BRT = Tessella.BoundaryRecovery

function brt_matrix(points, triangles)
    coords = Matrix{Float64}(undef, 3, length(points))
    for (i, p) in enumerate(points)
        coords[:, i] .= p
    end
    tris = Matrix{Int32}(undef, 3, length(triangles))
    for (i, t) in enumerate(triangles)
        tris[:, i] .= t
    end
    return coords, tris
end

function brt_icosphere(level)
    t = (1 + sqrt(5)) / 2
    unit(p) = (s = sqrt(p[1]^2 + p[2]^2 + p[3]^2); (p[1] / s, p[2] / s, p[3] / s))
    V = [unit(p) for p in ((-1.0, t, 0.0), (1.0, t, 0.0), (-1.0, -t, 0.0), (1.0, -t, 0.0), (0.0, -1.0, t),
                           (0.0, 1.0, t), (0.0, -1.0, -t), (0.0, 1.0, -t), (t, 0.0, -1.0), (t, 0.0, 1.0),
                           (-t, 0.0, -1.0), (-t, 0.0, 1.0))]
    F = [(1, 12, 6), (1, 6, 2), (1, 2, 8), (1, 8, 11), (1, 11, 12), (2, 6, 10), (6, 12, 5), (12, 11, 3),
         (11, 8, 7), (8, 2, 9), (4, 10, 5), (4, 5, 3), (4, 3, 7), (4, 7, 9), (4, 9, 10), (5, 10, 6),
         (3, 5, 12), (7, 3, 11), (9, 7, 8), (10, 9, 2)]
    for _ in 1:level
        middle = Dict{Tuple{Int,Int},Int}()
        function midpoint(a, b)
            get!(middle, minmax(a, b)) do
                push!(V, unit(((V[a][1] + V[b][1]) / 2, (V[a][2] + V[b][2]) / 2, (V[a][3] + V[b][3]) / 2)))
                length(V)
            end
        end
        G = NTuple{3,Int}[]
        for (a, b, c) in F
            ab = midpoint(a, b); bc = midpoint(b, c); ca = midpoint(c, a)
            push!(G, (a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca))
        end
        F = G
    end
    return brt_matrix(V, F)
end

# A torus grid whose quadrilaterals alternate their diagonal, so that about
# half of the surface edges are not Delaunay edges of the vertex set.
function brt_torus(nu, nv; major=2.0, minor=0.6)
    id(i, j) = mod(i, nu) * nv + mod(j, nv) + 1
    V = Vector{NTuple{3,Float64}}(undef, nu * nv)
    for i in 0:nu-1, j in 0:nv-1
        u = 2pi * i / nu; v = 2pi * j / nv
        V[id(i, j)] = ((major + minor * cos(v)) * cos(u), (major + minor * cos(v)) * sin(u), minor * sin(v))
    end
    F = NTuple{3,Int}[]
    for i in 0:nu-1, j in 0:nv-1
        a = id(i, j); b = id(i + 1, j); c = id(i + 1, j + 1); d = id(i, j + 1)
        if iseven(i + j)
            push!(F, (a, b, c), (a, c, d))
        else
            push!(F, (a, b, d), (b, c, d))
        end
    end
    return brt_matrix(V, F)
end

# The surface of the unit cube carrying an n×n grid on each face, with the
# diagonals alternating: vertices of a face are exactly coplanar.
function brt_box(n; lo=0.0, hi=1.0, flip=false)
    index = Dict{NTuple{3,Int},Int}(); V = NTuple{3,Float64}[]
    at(k) = lo + (hi - lo) * k / n
    function vertex(i, j, k)
        get!(index, (i, j, k)) do
            push!(V, (at(i), at(j), at(k))); length(V)
        end
    end
    F = NTuple{3,Int}[]
    function quad(a, b, c, d, parity)
        corners = flip ? (a, d, c, b) : (a, b, c, d)
        if parity
            push!(F, (corners[1], corners[2], corners[3]), (corners[1], corners[3], corners[4]))
        else
            push!(F, (corners[1], corners[2], corners[4]), (corners[2], corners[3], corners[4]))
        end
    end
    for i in 0:n-1, j in 0:n-1
        parity = iseven(i + j)
        quad(vertex(i, j, 0), vertex(i, j + 1, 0), vertex(i + 1, j + 1, 0), vertex(i + 1, j, 0), parity)
        quad(vertex(i, j, n), vertex(i + 1, j, n), vertex(i + 1, j + 1, n), vertex(i, j + 1, n), parity)
        quad(vertex(i, 0, j), vertex(i + 1, 0, j), vertex(i + 1, 0, j + 1), vertex(i, 0, j + 1), parity)
        quad(vertex(i, n, j), vertex(i, n, j + 1), vertex(i + 1, n, j + 1), vertex(i + 1, n, j), parity)
        quad(vertex(0, i, j), vertex(0, i, j + 1), vertex(0, i + 1, j + 1), vertex(0, i + 1, j), parity)
        quad(vertex(n, i, j), vertex(n, i + 1, j), vertex(n, i + 1, j + 1), vertex(n, i, j + 1), parity)
    end
    return brt_matrix(V, F)
end

function brt_signed_volume(coords, cells)
    total = 0.0
    for c in axes(cells, 2)
        a = coords[:, cells[1, c]]; b = coords[:, cells[2, c]]; d = coords[:, cells[3, c]]
        if size(cells, 1) == 3
            total += (a[1] * (b[2] * d[3] - b[3] * d[2]) - a[2] * (b[1] * d[3] - b[3] * d[1]) +
                      a[3] * (b[1] * d[2] - b[2] * d[1])) / 6
        else
            e = coords[:, cells[4, c]]
            u = b - a; v = d - a; w = e - a
            total += (u[1] * (v[2] * w[3] - v[3] * w[2]) - u[2] * (v[1] * w[3] - v[3] * w[1]) +
                      u[3] * (v[1] * w[2] - v[2] * w[1])) / 6
        end
    end
    return total
end

# The contract: input rows kept, every tetrahedron positive, and the faces
# owned by a single tetrahedron are exactly the input triangles.
function brt_conforms(coords, tris, mesh)
    mesh.coords[:, 1:size(coords, 2)] == coords || return false
    MT.validate(mesh).ok || return false
    counts = Dict{NTuple{3,Int32},Int}()
    for t in axes(mesh.tets, 2), skip in 1:4
        face = Tuple(sort!(Int32[mesh.tets[i, t] for i in 1:4 if i != skip]))
        counts[face] = get(counts, face, 0) + 1
    end
    boundary = Set(face for (face, count) in counts if count == 1)
    wanted = Set(Tuple(sort!(Int32[tris[1, f], tris[2, f], tris[3, f]])) for f in axes(tris, 2))
    return boundary == wanted && all(count -> count <= 2, values(counts))
end

@testset "boundary recovery" begin
    @testset "a single tetrahedron" begin
        coords, tris = brt_matrix([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
                                  [(1, 3, 2), (1, 2, 4), (2, 3, 4), (1, 4, 3)])
        mesh, reason, stats = BRT.recover_boundary_flips(coords, tris)
        @test reason == ""
        @test MT.ntets(mesh) == 1 && stats.steiner == 0
        @test brt_conforms(coords, tris, mesh)
    end

    @testset "convex surface needs no recovery" begin
        coords, tris = brt_icosphere(1)
        mesh, reason, stats = BRT.recover_boundary_flips(coords, tris)
        @test reason == ""
        @test stats.missing_edges == 0 && stats.missing_faces == 0 && stats.steiner == 0
        @test brt_conforms(coords, tris, mesh)
        @test brt_signed_volume(mesh.coords, mesh.tets) ≈ abs(brt_signed_volume(coords, tris)) rtol = 1e-12
    end

    @testset "non-convex surface with non-Delaunay edges" begin
        coords, tris = brt_torus(12, 6)
        mesh, reason, stats = BRT.recover_boundary_flips(coords, tris)
        @test reason == ""
        @test stats.missing_edges > 0
        @test brt_conforms(coords, tris, mesh)
        @test brt_signed_volume(mesh.coords, mesh.tets) ≈ abs(brt_signed_volume(coords, tris)) rtol = 1e-12
        # Added vertices are interior: they follow the input rows and no
        # boundary face uses them.
        @test size(mesh.coords, 2) == size(coords, 2) + stats.steiner
    end

    @testset "exactly coplanar vertices" begin
        coords, tris = brt_box(3)
        mesh, reason, stats = BRT.recover_boundary_flips(coords, tris)
        @test reason == ""
        @test brt_conforms(coords, tris, mesh)
        @test brt_signed_volume(mesh.coords, mesh.tets) ≈ 1.0 rtol = 1e-12
        # The result does not depend on anything but the input.
        again, _, _ = BRT.recover_boundary_flips(coords, tris)
        @test again.tets == mesh.tets && again.coords == mesh.coords
    end

    @testset "a cavity is not filled" begin
        outer, outer_tris = brt_box(2)
        inner, inner_tris = brt_box(1; lo=0.25, hi=0.75, flip=true)
        coords = hcat(outer, inner)
        tris = hcat(outer_tris, inner_tris .+ Int32(size(outer, 2)))
        mesh, reason, _ = BRT.recover_boundary_flips(coords, tris)
        @test reason == ""
        @test brt_conforms(coords, tris, mesh)
        @test brt_signed_volume(mesh.coords, mesh.tets) ≈ 1.0 - 0.125 rtol = 1e-12
    end

    @testset "a vertex inside the surface becomes an interior vertex" begin
        coords, tris = brt_box(1)
        coords = hcat(coords, [0.5, 0.5, 0.5])
        mesh, reason, _ = BRT.recover_boundary_flips(coords, tris)
        @test reason == ""
        @test brt_conforms(coords, tris, mesh)
        @test any(==(Int32(size(coords, 2))), mesh.tets)
    end

    @testset "inputs that are not a boundary are refused with a reason" begin
        coords, tris = brt_box(1)
        mesh, reason, _ = BRT.recover_boundary_flips(coords, tris[:, 1:end-1])
        @test mesh === nothing && occursin("not closed", reason)

        doubled = hcat(coords, coords[:, 1])
        mesh, reason, _ = BRT.recover_boundary_flips(doubled, tris)
        @test mesh === nothing && occursin("share coordinates", reason)

        mesh, reason, _ = BRT.recover_boundary_flips(coords, hcat(tris, tris[:, 1]))
        @test mesh === nothing && occursin("duplicate facet", reason)

        bad = copy(coords); bad[1, 1] = NaN
        mesh, reason, _ = BRT.recover_boundary_flips(bad, tris)
        @test mesh === nothing && occursin("nonfinite", reason)
    end
end
