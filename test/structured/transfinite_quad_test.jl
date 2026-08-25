using Test
using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, mixed_crc
using Tessella.MeshTypes: triangle_area

if !isdefined(Tessella, :TransfiniteQuad)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "structured",
                          "TransfiniteQuad.jl"))
end
using Tessella.TransfiniteQuad: mesh_transfinite_quad_patch

function _quad_rectangle(horizontal::Int, vertical::Int;
                         xmin=0.0, xmax=Float64(horizontal),
                         ymin=0.0, ymax=Float64(vertical))
    horizontal > 0 && vertical > 0 || throw(ArgumentError(
        "positive cell counts required"))
    bottom = [(xmin + (xmax - xmin) * i / horizontal, ymin, 0.0)
              for i in 0:horizontal]
    right = [(xmax, ymin + (ymax - ymin) * j / vertical, 0.0)
             for j in 0:vertical]
    top = [(xmax - (xmax - xmin) * i / horizontal, ymax, 0.0)
           for i in 0:horizontal]
    left = [(xmin, ymax - (ymax - ymin) * j / vertical, 0.0)
            for j in 0:vertical]
    return bottom, right, top, left
end

@inline _quad_node(i::Int, j::Int, width::Int) = Int32(i + 1 + j * width)

@inline function _mixed_node(mesh, index::Integer)
    return (mesh.coords[1, index], mesh.coords[2, index], mesh.coords[3, index])
end

function _expected_quadrangles(horizontal::Int, vertical::Int)
    width = horizontal + 1
    result = Matrix{Int32}(undef, 4, horizontal * vertical)
    cursor = 0
    for i in 0:horizontal-1, j in 0:vertical-1
        cursor += 1
        result[:, cursor] .= (_quad_node(i, j, width),
                              _quad_node(i + 1, j, width),
                              _quad_node(i + 1, j + 1, width),
                              _quad_node(i, j + 1, width))
    end
    return result
end

function _quad_area(mesh, block)
    area = 0.0
    for cell in axes(block.nodes, 2)
        p1 = _mixed_node(mesh, block.nodes[1, cell])
        p2 = _mixed_node(mesh, block.nodes[2, cell])
        p3 = _mixed_node(mesh, block.nodes[3, cell])
        p4 = _mixed_node(mesh, block.nodes[4, cell])
        area += triangle_area(p1, p2, p4) + triangle_area(p4, p2, p3)
    end
    return area
end

function _quad_polygon_area(sides)
    ring = NTuple{3,Float64}[]
    for side in sides
        append!(ring, @view side[1:end-1])
    end
    area = 0.0
    for index in eachindex(ring)
        point = ring[index]
        next = ring[mod1(index + 1, length(ring))]
        area += point[1] * next[2] - point[2] * next[1]
    end
    return abs(area) / 2
end

function _quad_edge_set(matrix)
    result = Set{NTuple{2,Int32}}()
    for index in axes(matrix, 2)
        first = matrix[1, index]
        second = matrix[2, index]
        push!(result, first < second ? (first, second) : (second, first))
    end
    return result
end

function _quad_boundary_edges(quadrangles)
    incidence = Dict{NTuple{2,Int32},Int}()
    for cell in axes(quadrangles, 2)
        nodes = quadrangles[:, cell]
        for (first, second) in ((nodes[1], nodes[2]), (nodes[2], nodes[3]),
                                (nodes[3], nodes[4]), (nodes[4], nodes[1]))
            key = first < second ? (first, second) : (second, first)
            incidence[key] = get(incidence, key, 0) + 1
        end
    end
    return Set(edge for (edge, count) in incidence if count == 1),
           maximum(values(incidence); init=0)
end

struct _CountOnlyQuadSide <: AbstractVector{NTuple{3,Float64}}
    count::Int
end
Base.size(side::_CountOnlyQuadSide) = (side.count,)
Base.getindex(::_CountOnlyQuadSide, ::Int) =
    error("resource counts were not checked before input conversion")

@noinline function _quad_allocated(sides)
    GC.gc()
    return @allocated mesh_transfinite_quad_patch(
        sides...; arrangement=:alternate_right)
end

@noinline function _quad_rejected_allocated(side)
    GC.gc()
    return @allocated try
        mesh_transfinite_quad_patch(side, side, side, side)
    catch err
        err isa ArgumentError || rethrow()
    end
end

@testset "Gmsh recombined four-sided transfinite quadrangles" begin
    @testset "type-3 topology, tags, CRC, and arrangement semantics" begin
        sides = _quad_rectangle(3, 2)
        meshes = Dict{Symbol,MixedMesh}()
        for arrangement in (:left, :right, :alternate_left, :alternate_right)
            mesh = mesh_transfinite_quad_patch(
                sides...; arrangement=arrangement, face_tag=21,
                side_tags=(11, 12, 13, 14))
            meshes[arrangement] = mesh
            @test Tessella.Elements.validate(mesh).ok
            @test size(mesh.coords) == (3, 12)
            @test length(mesh.blocks) == 2
            lines, quads = mesh.blocks
            @test lines isa ElementBlock
            @test quads isa ElementBlock
            @test (lines.msh, quads.msh) == (1, 3)
            @test size(lines.nodes) == (2, 10)
            @test size(quads.nodes) == (4, 6)
            @test quads.nodes == _expected_quadrangles(3, 2)
            @test lines.tags == Int32[11, 11, 11, 12, 12,
                                      13, 13, 13, 14, 14]
            @test quads.tags == fill(Int32(21), 6)
            boundary, max_incidence = _quad_boundary_edges(quads.nodes)
            @test max_incidence == 2
            @test boundary == _quad_edge_set(lines.nodes)
            @test _quad_area(mesh, quads) == 6.0
            @test mixed_crc(mesh) == mixed_crc(mesh_transfinite_quad_patch(
                sides...; arrangement=arrangement, face_tag=21,
                side_tags=(11, 12, 13, 14)))
        end
        crcs = Set(mixed_crc(mesh).sha for mesh in values(meshes))
        @test crcs == Set([
            "d05e9cbc57de975f1f88d8fb1da62e0636ea4bbe63a038d3ca060b16201b0ea2"])
        @test all(mesh.blocks[2].nodes == meshes[:left].blocks[2].nodes
                  for mesh in values(meshes))
    end

    @testset "average-chord Coons nodes and exact boundary preservation" begin
        bottom = [(0.0, 0.0, 0.0), (0.7, -0.25, 0.0),
                  (2.1, -0.55, 0.0), (3.2, -0.2, 0.0), (4.0, 0.0, 0.0)]
        right = [(4.0, 0.0, 0.0), (4.35, 0.8, 0.0),
                 (4.2, 2.1, 0.0), (4.0, 3.0, 0.0)]
        top = [(4.0, 3.0, 0.0), (3.1, 3.55, 0.0),
               (2.0, 3.35, 0.0), (0.8, 3.15, 0.0), (0.0, 3.0, 0.0)]
        left = [(0.0, 3.0, 0.0), (-0.3, 2.15, 0.0),
                (-0.5, 0.9, 0.0), (0.0, 0.0, 0.0)]
        sides = (bottom, right, top, left)
        mesh = mesh_transfinite_quad_patch(sides...)
        @test Tessella.Elements.validate(mesh).ok
        @test (size(mesh.coords, 2), size(mesh.blocks[1].nodes, 2),
               size(mesh.blocks[2].nodes, 2)) == (20, 14, 12)
        @test _quad_area(mesh, mesh.blocks[2]) ≈ _quad_polygon_area(sides) atol=128eps(Float64)

        width = 5
        for i in 0:4
            @test _mixed_node(mesh, _quad_node(i, 0, width)) == bottom[i + 1]
            @test _mixed_node(mesh, _quad_node(i, 3, width)) == top[5 - i]
        end
        for j in 0:3
            @test _mixed_node(mesh, _quad_node(4, j, width)) == right[j + 1]
            @test _mixed_node(mesh, _quad_node(0, j, width)) == left[4 - j]
        end

        chord(first, second) = hypot(first[1] - second[1],
                                     first[2] - second[2],
                                     first[3] - second[3])
        top_grid = reverse(top)
        left_grid = reverse(left)
        ustep = [0.5 * (chord(bottom[i + 1], bottom[i]) +
                        chord(top_grid[i + 1], top_grid[i])) for i in 1:4]
        vstep = [0.5 * (chord(right[j + 1], right[j]) +
                        chord(left_grid[j + 1], left_grid[j])) for j in 1:3]
        i = 2
        j = 1
        u = sum(ustep[1:i]) / sum(ustep)
        v = sum(vstep[1:j]) / sum(vstep)
        c1 = bottom[1]
        c2 = bottom[end]
        c3 = top_grid[end]
        c4 = top_grid[1]
        expected = ntuple(3) do coordinate
            (1 - u) * left_grid[j + 1][coordinate] +
            u * right[j + 1][coordinate] +
            (1 - v) * bottom[i + 1][coordinate] +
            v * top_grid[i + 1][coordinate] -
            ((1 - u) * (1 - v) * c1[coordinate] +
             u * (1 - v) * c2[coordinate] + u * v * c3[coordinate] +
             (1 - u) * v * c4[coordinate])
        end
        actual = _mixed_node(mesh, _quad_node(i, j, width))
        @test all(isapprox(actual[coordinate], expected[coordinate];
                           atol=32eps(Float64), rtol=32eps(Float64))
                  for coordinate in 1:3)
    end

    @testset "tilted, clockwise, and minimal patches" begin
        base = _quad_rectangle(4, 3; xmax=2.0, ymax=1.5)
        ex = (inv(sqrt(2.0)), inv(sqrt(2.0)), 0.0)
        ey = (-inv(sqrt(6.0)), inv(sqrt(6.0)), 2inv(sqrt(6.0)))
        origin = (1.0, -2.0, 3.0)
        transform(point) =
            (origin[1] + point[1] * ex[1] + point[2] * ey[1],
             origin[2] + point[1] * ex[2] + point[2] * ey[2],
             origin[3] + point[1] * ex[3] + point[2] * ey[3])
        tilted = map(side -> transform.(side), base)
        mesh = mesh_transfinite_quad_patch(tilted...)
        @test Tessella.Elements.validate(mesh).ok
        @test _quad_area(mesh, mesh.blocks[2]) ≈ 3.0 atol=512eps(Float64)

        clockwise = ([(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0, 0.0)],
                     [(0.0, 2.0, 0.0), (1.0, 2.0, 0.0),
                      (2.0, 2.0, 0.0), (3.0, 2.0, 0.0)],
                     [(3.0, 2.0, 0.0), (3.0, 1.0, 0.0), (3.0, 0.0, 0.0)],
                     [(3.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                      (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
        reversed_mesh = mesh_transfinite_quad_patch(clockwise...)
        @test Tessella.Elements.validate(reversed_mesh).ok
        @test _quad_area(reversed_mesh, reversed_mesh.blocks[2]) == 6.0

        minimal = mesh_transfinite_quad_patch(_quad_rectangle(1, 1)...)
        @test size(minimal.coords, 2) == 4
        @test size(minimal.blocks[1].nodes, 2) == 4
        @test size(minimal.blocks[2].nodes, 2) == 1
        @test Tessella.Elements.validate(minimal).ok

        origin = (0.0, 0.0, 0.0)
        u = (1.0, 1.0, 1.0)
        v = (1.0, 1.0 + 2.0^-10, 1.0)
        opposite = ntuple(d -> u[d] + v[d], 3)
        thin = ([origin, u], [u, opposite], [opposite, v], [v, origin])
        thin_mesh = mesh_transfinite_quad_patch(thin...)
        @test Tessella.Elements.validate(thin_mesh).ok
        @test size(thin_mesh.blocks[2].nodes, 2) == 1

        # Each coordinate and corner identity is exactly representable, but the
        # Float64 Newell products cancel. The rare exact projection certificate
        # must preserve this valid affine patch.
        extent = Float64(Int64(1) << 27)
        exact_u = (0.0, extent, extent + 1.0)
        exact_v = (0.0, extent + 1.0, extent + 2.0)
        exact_opposite = ntuple(d -> exact_u[d] + exact_v[d], 3)
        cancellation = ([(0.0, 0.0, 0.0), exact_u],
                        [exact_u, exact_opposite],
                        [exact_opposite, exact_v],
                        [exact_v, (0.0, 0.0, 0.0)])
        cancellation_mesh = mesh_transfinite_quad_patch(cancellation...)
        @test Tessella.Elements.validate(cancellation_mesh).ok
        @test size(cancellation_mesh.blocks[2].nodes, 2) == 1

        # The normal is nonzero in Float64 here, but normalizing these raw
        # coordinates independently destroys the exact unit projected area.
        # Geometry predicates must see the original dropped coordinates.
        raw_u = (extent, extent + 1.0, extent)
        raw_v = (extent + 1.0, extent + 2.0, extent + 1.0)
        raw_opposite = ntuple(d -> raw_u[d] + raw_v[d], 3)
        raw_projection = ([(0.0, 0.0, 0.0), raw_u],
                          [raw_u, raw_opposite],
                          [raw_opposite, raw_v],
                          [raw_v, (0.0, 0.0, 0.0)])
        raw_mesh = mesh_transfinite_quad_patch(raw_projection...)
        @test Tessella.Elements.validate(raw_mesh).ok
        @test size(raw_mesh.blocks[2].nodes, 2) == 1
    end

    @testset "validated blockers and pre-allocation resource limits" begin
        sides = _quad_rectangle(3, 2)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            tuple(sides[1]...),sides[2],sides[3],sides[4])
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            1,sides[2],sides[3],sides[4])
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            [(0.0, 0.0, 0.0)], sides[2], sides[3], sides[4])
        mismatched = copy(sides[2])
        insert!(mismatched, 2, (3.0, 0.5, 0.0))
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides[1], mismatched, sides[3], sides[4])
        broken = copy(sides[2])
        broken[1] = (3.0, nextfloat(0.0), 0.0)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides[1], broken, sides[3], sides[4])
        nonfinite = copy(sides[1])
        nonfinite[2] = (NaN, 0.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            nonfinite, sides[2], sides[3], sides[4])
        extra = Any[(point..., point == sides[1][2] ? NaN : 0.0)
                    for point in sides[1]]
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            extra, sides[2], sides[3], sides[4])
        duplicate = copy(sides[1])
        duplicate[2] = duplicate[1]
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            duplicate, sides[2], sides[3], sides[4])
        nonplanar = copy(sides[1])
        nonplanar[2] = (1.0, 0.0, 1.0e-6)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            nonplanar, sides[2], sides[3], sides[4])

        bowtie = ([(0.0, 0.0, 0.0), (1.0, 1.0, 0.0)],
                  [(1.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
                  [(0.0, 1.0, 0.0), (1.0, 0.0, 0.0)],
                  [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
        @test_throws ArgumentError mesh_transfinite_quad_patch(bowtie...)
        folded = ([(0.0, 0.0, 0.0),
                   (1.138676204796158, -0.8998793245683743, 0.0),
                   (1.0, 0.0, 0.0)],
                  [(1.0, 0.0, 0.0),
                   (0.9753676924206318, -0.6302141619964883, 0.0),
                   (1.0, 1.0, 0.0)],
                  [(1.0, 1.0, 0.0),
                   (-0.14986440109197496, 1.9513598146977102, 0.0),
                   (0.0, 1.0, 0.0)],
                  [(0.0, 1.0, 0.0),
                   (-0.9877714621697844, -0.56016161339637, 0.0),
                   (0.0, 0.0, 0.0)])
        @test_throws ArgumentError mesh_transfinite_quad_patch(folded...)

        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; arrangement=:alternate)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; arrangement="Left")
        @test_throws ArgumentError mesh_transfinite_quad_patch(sides...; face_tag=true)
        @test_throws ArgumentError mesh_transfinite_quad_patch(sides...; face_tag=-1)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; side_tags=(1, 2, 3))
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; side_tags=(1, 2, true, 4))
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; max_nodes=11)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; max_quadrangles=false)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; max_quadrangles=5)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            sides...; max_quadrangles=6.0)
        bounded = mesh_transfinite_quad_patch(
            sides...; max_nodes=12, max_quadrangles=6)
        @test (size(bounded.coords, 2), size(bounded.blocks[2].nodes, 2)) == (12, 6)

        huge = _CountOnlyQuadSide(typemax(Int))
        @test_throws ArgumentError mesh_transfinite_quad_patch(huge, huge, huge, huge)
        wide = _CountOnlyQuadSide(100_000)
        @test_throws ArgumentError mesh_transfinite_quad_patch(wide, wide, wide, wide)
        certification_limited = _CountOnlyQuadSide(40_001)
        @test_throws ArgumentError mesh_transfinite_quad_patch(
            certification_limited, certification_limited,
            certification_limited, certification_limited;
            max_nodes=typemax(Int32), max_quadrangles=typemax(Int32))
        @test _quad_rejected_allocated(wide) < 64_000

        maximum = floatmax(Float64)
        overflowing = ([(maximum, 0.0, 0.0), (-maximum, 0.0, 0.0)],
                       [(-maximum, 0.0, 0.0), (-maximum, 1.0, 0.0)],
                       [(-maximum, 1.0, 0.0), (maximum, 1.0, 0.0)],
                       [(maximum, 1.0, 0.0), (maximum, 0.0, 0.0)])
        @test_throws ArgumentError mesh_transfinite_quad_patch(overflowing...)
        @test isempty(Test.detect_ambiguities(
            Tessella.TransfiniteQuad; recursive=true))
    end

    @testset "allocation growth remains linear in output size" begin
        small_sides = _quad_rectangle(64, 64)
        large_sides = _quad_rectangle(128, 64)
        mesh_transfinite_quad_patch(small_sides...)
        mesh_transfinite_quad_patch(large_sides...)
        small = _quad_allocated(small_sides)
        large = _quad_allocated(large_sides)
        @test small > 0
        @test large > small
        @test large <= 2.30small + 524_288
        @info "transfinite quad allocation ratchet" small_bytes=small large_bytes=large
    end
end
