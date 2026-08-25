using Test
using Tessella
using Tessella.MeshTypes: boundary_edges, mesh_crc, nnodes, node, nsegs,
                          ntris, triangle_area, validate

if !isdefined(Tessella, :TransfiniteTriangle)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "structured",
                          "TransfiniteTriangle.jl"))
end
using Tessella.TransfiniteTriangle: mesh_transfinite_triangle,
                                    mesh_transfinite_triangle_patch

@inline _triangle_node(i::Int, j::Int) = Int32((i * (i + 1)) ÷ 2 + j + 1)

@inline function _triangle_lerp3(a, b, t)
    ((1 - t) * a[1] + t * b[1],
     (1 - t) * a[2] + t * b[2],
     (1 - t) * a[3] + t * b[3])
end

function _straight_triangle(divisions::Int;
                            corners=((0.0, 0.0, 0.0),
                                     (4.0, 0.0, 0.0),
                                     (0.5, 3.0, 0.0)))
    divisions > 0 || throw(ArgumentError("positive divisions required"))
    c1, c2, c3 = corners
    side1 = [_triangle_lerp3(c1, c2, i / divisions) for i in 0:divisions]
    side2 = [_triangle_lerp3(c2, c3, i / divisions) for i in 0:divisions]
    side3 = [_triangle_lerp3(c3, c1, i / divisions) for i in 0:divisions]
    return side1, side2, side3
end

function _triangle_surface_area(mesh)
    return sum(triangle_area(node(mesh, mesh.tris[1, triangle]),
                             node(mesh, mesh.tris[2, triangle]),
                             node(mesh, mesh.tris[3, triangle]))
               for triangle in 1:ntris(mesh); init=0.0)
end

function _triangle_mixed_surface_area(mesh)
    mixed_node(index) = (mesh.coords[1, Int(index)],
                         mesh.coords[2, Int(index)],
                         mesh.coords[3, Int(index)])
    area = 0.0
    for block in mesh.blocks
        if block.msh == 2
            for cell in axes(block.nodes, 2)
                area += triangle_area(
                    mixed_node(block.nodes[1, cell]),
                    mixed_node(block.nodes[2, cell]),
                    mixed_node(block.nodes[3, cell]))
            end
        elseif block.msh == 3
            for cell in axes(block.nodes, 2)
                first = mixed_node(block.nodes[1, cell])
                second = mixed_node(block.nodes[2, cell])
                third = mixed_node(block.nodes[3, cell])
                fourth = mixed_node(block.nodes[4, cell])
                area += triangle_area(first, second, third)
                area += triangle_area(first, third, fourth)
            end
        end
    end
    return area
end

function _triangle_block(mesh, msh::Int)
    matches = filter(block -> block.msh == msh, mesh.blocks)
    length(matches) == 1 || error(
        "expected one mixed block for MSH type $msh, got $(length(matches))")
    return only(matches)
end

function _triangle_canonical_block_cells(block)
    cells = Tuple{Tuple,Int32}[]
    for cell in axes(block.nodes, 2)
        nodes = Tuple(block.nodes[:, cell])
        node_count = length(nodes)
        candidates = Tuple[]
        for shift in 0:node_count-1
            push!(candidates,
                  ntuple(offset -> nodes[mod1(shift + offset, node_count)],
                         node_count))
        end
        reversed = reverse(nodes)
        for shift in 0:node_count-1
            push!(candidates,
                  ntuple(offset -> reversed[mod1(shift + offset, node_count)],
                         node_count))
        end
        push!(cells, (minimum(candidates), block.tags[cell]))
    end
    sort!(cells)
    return cells
end

function _triangle_edge_set(matrix)
    result = Set{NTuple{2,Int32}}()
    for index in axes(matrix, 2)
        a = matrix[1, index]
        b = matrix[2, index]
        push!(result, a < b ? (a, b) : (b, a))
    end
    return result
end

function _triangle_canonical_triangles(mesh)
    result = NTuple{3,Int32}[]
    for triangle in axes(mesh.tris, 2)
        values = sort(mesh.tris[:, triangle])
        push!(result, (values[1], values[2], values[3]))
    end
    sort!(result)
    return result
end

function _expected_triangles(divisions::Int)
    result = NTuple{3,Int32}[]
    for i in 0:divisions-1, j in 0:i
        v1 = _triangle_node(i, j)
        v2 = _triangle_node(i + 1, j)
        v3 = _triangle_node(i + 1, j + 1)
        if i > 0 && j < i
            values = sort(Int32[v1, v3, _triangle_node(i, j + 1)])
            push!(result, (values[1], values[2], values[3]))
        end
        values = sort(Int32[v1, v2, v3])
        push!(result, (values[1], values[2], values[3]))
    end
    sort!(result)
    return result
end

function _triangle_polygon_area(sides)
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

function _chords(side)
    result = zeros(Float64, length(side))
    for index in 1:length(side)-1
        first = side[index]
        second = side[index + 1]
        result[index + 1] = result[index] +
            hypot(second[1] - first[1], second[2] - first[2],
                  second[3] - first[3])
    end
    return result
end

function _sample_chord(side, cumulative, t)
    total = cumulative[end]
    right = findfirst(index -> t <= cumulative[index] / total,
                      2:length(cumulative))
    right === nothing && error("chord sample not bracketed")
    right = right + 1
    left = right - 1
    fraction = (t * total - cumulative[left]) /
               (cumulative[right] - cumulative[left])
    return _triangle_lerp3(side[left], side[right], fraction)
end

struct _CountOnlyTriangleSide <: AbstractVector{NTuple{3,Float64}}
    count::Int
end
Base.size(side::_CountOnlyTriangleSide) = (side.count,)
Base.getindex(::_CountOnlyTriangleSide, ::Int) =
    error("resource counts were not checked before input conversion")

@noinline function _triangle_allocated(sides)
    GC.gc()
    return @allocated mesh_transfinite_triangle(
        sides...; arrangement=:alternate_left)
end

@noinline function _triangle_patch_allocated(sides)
    GC.gc()
    return @allocated mesh_transfinite_triangle_patch(
        sides...; arrangement=:alternate_left)
end

@noinline function _triangle_rejected_allocated(side)
    GC.gc()
    return @allocated try
        mesh_transfinite_triangle(side, side, side)
    catch err
        err isa ArgumentError || rethrow()
    end
end

@testset "Gmsh three-sided planar transfinite triangles" begin
    @testset "structured topology, tags, CRC, and arrangement semantics" begin
        sides = _straight_triangle(4)
        meshes = Dict{Symbol,Tessella.MeshTypes.Mesh}()
        for arrangement in (:left, :right, :alternate_left, :alternate_right)
            mesh = mesh_transfinite_triangle(
                sides...; arrangement=arrangement, face_tag=21,
                side_tags=(11, 12, 13))
            meshes[arrangement] = mesh
            @test validate(mesh).ok
            @test (nnodes(mesh), nsegs(mesh), ntris(mesh)) == (15, 12, 16)
            @test mesh_crc(mesh).bbox == ((0.0, 0.0, 0.0), (4.0, 3.0, 0.0))
            @test _triangle_canonical_triangles(mesh) == _expected_triangles(4)
            @test mesh.tri_tag == fill(Int32(21), 16)
            @test mesh.seg_tag == Int32[fill(11, 4); fill(12, 4); fill(13, 4)]
            boundary, max_incidence = boundary_edges(mesh.tris)
            @test max_incidence == 2
            @test Set(boundary) == _triangle_edge_set(mesh.segs)
            @test _triangle_surface_area(mesh) == 6.0
            @test mesh_crc(mesh) == mesh_crc(mesh_transfinite_triangle(
                sides...; arrangement=arrangement, face_tag=21,
                side_tags=(11, 12, 13)))
        end
        crcs = Set(mesh_crc(mesh).sha for mesh in values(meshes))
        @test crcs == Set([
            "5f231f0c22f6812c247514103e21653e15281194911bf456b4c00a9d190a40df"])
        reference = _triangle_canonical_triangles(meshes[:left])
        @test all(_triangle_canonical_triangles(mesh) == reference
                  for mesh in values(meshes))
    end

    @testset "Gmsh recombined triangle/quadrangle layouts" begin
        sides = _straight_triangle(4)
        # Exact local-node records from the pinned Gmsh 4.15.2
        # Mesh.TransfiniteTri=1 + Recombine Surface oracle.
        expected_triangles = Dict(
            :left => [(1, 2, 3), (3, 5, 6), (5, 8, 9), (8, 12, 13)],
            :right => [(1, 2, 3), (3, 5, 6), (6, 9, 10), (10, 14, 15)],
            :alternate_left =>
                [(1, 2, 3), (3, 5, 6), (4, 7, 8), (10, 14, 15)],
            :alternate_right =>
                [(1, 2, 3), (3, 5, 6), (4, 7, 8), (10, 14, 15)],
        )
        expected_quadrangles = Dict(
            :left => [(2, 4, 5, 3), (4, 7, 8, 5),
                      (5, 9, 10, 6), (7, 11, 12, 8),
                      (8, 13, 14, 9), (9, 14, 15, 10)],
            :right => [(2, 4, 5, 3), (4, 7, 8, 5),
                       (5, 8, 9, 6), (7, 11, 12, 8),
                       (8, 12, 13, 9), (9, 13, 14, 10)],
            :alternate_left =>
                [(2, 4, 5, 3), (4, 8, 9, 5),
                 (5, 9, 10, 6), (7, 11, 12, 8),
                 (8, 12, 13, 9), (9, 13, 14, 10)],
            :alternate_right =>
                [(2, 4, 5, 3), (4, 8, 9, 5),
                 (5, 9, 10, 6), (7, 11, 12, 8),
                 (8, 12, 13, 9), (9, 13, 14, 10)],
        )
        expected_crcs = Dict(
            :left =>
                "09d6619152fe4f42d604c9a95e3805843825a08615d92778ae4e551d85fa2ce3",
            :right =>
                "b401b6bc71cac6dc3f7f44b20a4128ddf4bb439acec941253b7598af19800205",
            :alternate_left =>
                "5d200b76825ed699e99b125c28cef49158d30bc56a9e090d8469854cc84dabfa",
            :alternate_right =>
                "5d200b76825ed699e99b125c28cef49158d30bc56a9e090d8469854cc84dabfa",
        )
        certified = mesh_transfinite_triangle(
            sides...; face_tag=21, side_tags=(11, 12, 13))
        for arrangement in
            (:left, :right, :alternate_left, :alternate_right)
            mixed = mesh_transfinite_triangle_patch(
                sides...; arrangement, face_tag=21,
                side_tags=(11, 12, 13))
            @test Tessella.Elements.validate(mixed).ok
            @test [block.msh for block in mixed.blocks] == [1, 2, 3]
            @test mixed.coords == certified.coords
            @test mixed.blocks[1].nodes == certified.segs
            @test mixed.blocks[1].tags == certified.seg_tag
            triangle_block = _triangle_block(mixed, 2)
            quadrangle_block = _triangle_block(mixed, 3)
            @test [Tuple(triangle_block.nodes[:, cell])
                   for cell in axes(triangle_block.nodes, 2)] ==
                  expected_triangles[arrangement]
            @test [Tuple(quadrangle_block.nodes[:, cell])
                   for cell in axes(quadrangle_block.nodes, 2)] ==
                  expected_quadrangles[arrangement]
            @test triangle_block.tags == fill(Int32(21), 4)
            @test quadrangle_block.tags == fill(Int32(21), 6)
            @test _triangle_mixed_surface_area(mixed) == 6.0
            @test Tessella.Elements.mixed_crc(mixed).sha ==
                  expected_crcs[arrangement]
            @test Tessella.Elements.mixed_crc(mixed) ==
                  Tessella.Elements.mixed_crc(mesh_transfinite_triangle_patch(
                      sides...; arrangement, face_tag=21,
                      side_tags=(11, 12, 13)))
        end

        for divisions in 1:8
            mixed = mesh_transfinite_triangle_patch(
                _straight_triangle(divisions)...; arrangement=:left)
            @test size(mixed.coords, 2) ==
                  (divisions + 1) * (divisions + 2) ÷ 2
            @test size(_triangle_block(mixed, 2).nodes, 2) == divisions
            if divisions == 1
                @test [block.msh for block in mixed.blocks] == [1, 2]
            else
                @test size(_triangle_block(mixed, 3).nodes, 2) ==
                      divisions * (divisions - 1) ÷ 2
            end
            @test Tessella.Elements.validate(mixed).ok
        end
        roundtrip_source = mesh_transfinite_triangle_patch(
            sides...; arrangement=:left, face_tag=21,
            side_tags=(11, 12, 13))
        mktempdir() do directory
            for version in (2.2, 4.1), binary in (false, true)
                path = joinpath(
                    directory, "triangle-patch-$version-$binary.msh")
                Tessella.Elements.write_mixed_msh(
                    path, roundtrip_source; version, binary)
                roundtrip = Tessella.Elements.read_mixed_msh(path)
                @test Tessella.Elements.validate(roundtrip).ok
                @test roundtrip.coords == roundtrip_source.coords
                @test [block.msh for block in roundtrip.blocks] == [1, 2, 3]
                for msh in (1, 2, 3)
                    @test _triangle_canonical_block_cells(
                        _triangle_block(roundtrip, msh)) ==
                          _triangle_canonical_block_cells(
                              _triangle_block(roundtrip_source, msh))
                end
            end
        end
        @test Tessella.mesh_transfinite_triangle_patch ===
              mesh_transfinite_triangle_patch
    end

    @testset "Gmsh triangular interpolation and boundary conservation" begin
        side1 = [(0.0, 0.0, 0.0), (0.7, -0.2, 0.0),
                 (2.0, -0.5, 0.0), (3.2, -0.1, 0.0), (4.0, 0.0, 0.0)]
        side2 = [(4.0, 0.0, 0.0), (3.8, 0.8, 0.0),
                 (2.8, 1.7, 0.0), (1.7, 2.5, 0.0), (0.5, 3.0, 0.0)]
        side3 = [(0.5, 3.0, 0.0), (0.1, 2.4, 0.0),
                 (-0.2, 1.6, 0.0), (-0.3, 0.7, 0.0), (0.0, 0.0, 0.0)]
        sides = (side1, side2, side3)
        mesh = mesh_transfinite_triangle(sides...; arrangement=:right)
        @test validate(mesh).ok
        @test (nnodes(mesh), nsegs(mesh), ntris(mesh)) == (15, 12, 16)
        @test _triangle_surface_area(mesh) ≈ _triangle_polygon_area(sides) atol=128eps(Float64)

        for i in 0:4
            @test node(mesh, _triangle_node(i, 0)) == side1[i + 1]
            @test node(mesh, _triangle_node(i, i)) == side3[5 - i]
            @test node(mesh, _triangle_node(4, i)) == side2[i + 1]
        end

        # Independent direct evaluation of Gmsh 4.15.2 TRAN_TRI at row i=3,
        # ray j=1, including side-2 chord interpolation at t=j/i.
        opposite = reverse(side3)
        first_chords = _chords(side1)
        opposite_chords = _chords(opposite)
        averaged = zeros(Float64, 5)
        for index in 1:4
            averaged[index + 1] = averaged[index] +
                0.5 * (first_chords[index + 1] - first_chords[index]) +
                0.5 * (opposite_chords[index + 1] - opposite_chords[index])
        end
        i = 3
        j = 1
        u = averaged[i + 1] / averaged[end]
        v = j / i
        p2 = _sample_chord(side2, _chords(side2), v)
        expected = ntuple(3) do coordinate
            u * p2[coordinate] + (1 - v) * side1[i + 1][coordinate] +
            v * opposite[i + 1][coordinate] -
            (u * (1 - v) * side1[end][coordinate] +
             u * v * side2[end][coordinate])
        end
        actual = node(mesh, _triangle_node(i, j))
        @test all(isapprox(actual[coordinate], expected[coordinate];
                           atol=32eps(Float64), rtol=32eps(Float64))
                  for coordinate in 1:3)

        for arrangement in
            (:left, :right, :alternate_left, :alternate_right)
            mixed = mesh_transfinite_triangle_patch(
                sides...; arrangement, face_tag=17,
                side_tags=(3, 4, 5))
            @test mixed.coords == mesh.coords
            @test Tessella.Elements.validate(mixed).ok
            @test _triangle_mixed_surface_area(mixed) ≈
                  _triangle_polygon_area(sides) atol=128eps(Float64)
        end
    end

    @testset "tilted, clockwise, and long patches" begin
        base = _straight_triangle(5)
        ex = (inv(sqrt(2.0)), inv(sqrt(2.0)), 0.0)
        ey = (-inv(sqrt(6.0)), inv(sqrt(6.0)), 2inv(sqrt(6.0)))
        origin = (1.0, -2.0, 3.0)
        transform(point) =
            (origin[1] + point[1] * ex[1] + point[2] * ey[1],
             origin[2] + point[1] * ex[2] + point[2] * ey[2],
             origin[3] + point[1] * ex[3] + point[2] * ey[3])
        tilted = map(side -> transform.(side), base)
        mesh = mesh_transfinite_triangle(tilted...)
        @test validate(mesh).ok
        @test _triangle_surface_area(mesh) ≈ 6.0 atol=512eps(Float64)

        clockwise = _straight_triangle(
            4; corners=((0.0, 0.0, 0.0),
                        (0.0, 2.0, 0.0),
                        (3.0, 0.0, 0.0)))
        reversed_mesh = mesh_transfinite_triangle(
            clockwise...; arrangement=:alternate_right)
        @test validate(reversed_mesh).ok
        @test _triangle_surface_area(reversed_mesh) == 3.0
        reversed_patch = mesh_transfinite_triangle_patch(
            clockwise...; arrangement=:alternate_right)
        @test Tessella.Elements.validate(reversed_patch).ok
        @test _triangle_mixed_surface_area(reversed_patch) == 3.0

        long_mesh = mesh_transfinite_triangle(_straight_triangle(128)...)
        @test (nnodes(long_mesh), ntris(long_mesh)) == (8385, 16384)
        @test validate(long_mesh).ok
    end

    @testset "validated blockers and pre-allocation resource limits" begin
        sides = _straight_triangle(3)
        @test_throws ArgumentError mesh_transfinite_triangle(
            tuple(sides[1]...), sides[2], sides[3])
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            tuple(sides[1]...), sides[2], sides[3])
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            1, sides[2], sides[3])
        @test_throws ArgumentError mesh_transfinite_triangle(
            [(0.0, 0.0, 0.0)], sides[2], sides[3])
        @test_throws ArgumentError mesh_transfinite_triangle(
            _straight_triangle(2)[1], sides[2], sides[3])

        broken = copy(sides[2])
        broken[1] = (4.0, nextfloat(0.0), 0.0)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides[1], broken, sides[3])
        nonfinite = copy(sides[1])
        nonfinite[2] = (NaN, 0.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_triangle(
            nonfinite, sides[2], sides[3])
        extra = Any[(point..., point == sides[1][2] ? NaN : 0.0)
                    for point in sides[1]]
        @test_throws ArgumentError mesh_transfinite_triangle(
            extra, sides[2], sides[3])
        duplicate = copy(sides[1])
        duplicate[2] = duplicate[1]
        @test_throws ArgumentError mesh_transfinite_triangle(
            duplicate, sides[2], sides[3])
        nonplanar = copy(sides[1])
        nonplanar[2] = (nonplanar[2][1], nonplanar[2][2], 1.0e-6)
        @test_throws ArgumentError mesh_transfinite_triangle(
            nonplanar, sides[2], sides[3])

        bowtie = ([(0.0, 0.0, 0.0), (1.0, 1.5, 0.0), (0.0, 2.0, 0.0)],
                  [(0.0, 2.0, 0.0), (2.0, 0.0, 0.0), (1.0, 2.0, 0.0)],
                  [(1.0, 2.0, 0.0), (2.0, 2.0, 0.0), (0.0, 0.0, 0.0)])
        @test_throws ArgumentError mesh_transfinite_triangle(bowtie...)

        folded = ([(1.1338512180093638, 0.0, 0.0),
                   (1.4545499120266219, 2.5193543497749493, 0.0),
                   (-0.6291363944621956, 1.0896962000992183, 0.0)],
                  [(-0.6291363944621956, 1.0896962000992183, 0.0),
                   (-1.5694091536601567, 1.9219718965824676e-16, 0.0),
                   (-0.08730973152828207, -0.15122489100218267, 0.0)],
                  [(-0.08730973152828207, -0.15122489100218267, 0.0),
                   (0.48225203858372556, -0.835285032880679, 0.0),
                   (1.1338512180093638, 0.0, 0.0)])
        fold_error = try
            mesh_transfinite_triangle(folded...)
            nothing
        catch err
            err
        end
        @test fold_error isa ArgumentError
        @test occursin("reverses patch orientation", sprint(showerror, fold_error))
        @test_throws ArgumentError mesh_transfinite_triangle_patch(folded...)

        # The atomic triangle lattice is valid, but Gmsh's Left pairing would
        # make its final quadrangle concave. The mixed API blocks that cell.
        concave_recombined = (
            [(0.0, 0.0, 0.0),
             (0.5956665854630816, 0.337420243172981, 0.0),
             (2.2608011714552694, -1.2585914468342874, 0.0),
             (3.547736455291709, -1.1126591321601882, 0.0),
             (4.0, 0.0, 0.0)],
            [(4.0, 0.0, 0.0),
             (3.4456458160331973, 0.7367598098341835, 0.0),
             (2.9501456625323503, 1.5658637318983206, 0.0),
             (1.4731544508145447, 1.8598054426345922, 0.0),
             (0.5, 3.0, 0.0)],
            [(0.5, 3.0, 0.0),
             (1.068309161473124, 2.3031495146346597, 0.0),
             (0.33513926989785225, 2.0705405301075386, 0.0),
             (0.05902398752140388, 1.2101515997372627, 0.0),
             (0.0, 0.0, 0.0)])
        @test validate(mesh_transfinite_triangle(concave_recombined...)).ok
        concave_error = try
            mesh_transfinite_triangle_patch(
                concave_recombined...; arrangement=:left)
            nothing
        catch err
            err
        end
        @test concave_error isa ArgumentError
        @test occursin("quadrangle 6", sprint(showerror, concave_error))

        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; arrangement=:alternate)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; arrangement="Left")
        @test_throws ArgumentError mesh_transfinite_triangle(sides...; face_tag=true)
        @test_throws ArgumentError mesh_transfinite_triangle(sides...; face_tag=-1)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; face_tag=big(typemax(Int32)) + 1)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; side_tags=(1, 2))
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; side_tags=(1, true, 3))
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; side_tags=(1, -2, 3))
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_triangles=false)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_nodes=-1)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_nodes=10.0)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; arrangement=:alternate)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; max_nodes=9)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; max_triangles=2)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; max_quadrangles=2)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; max_quadrangles=true)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; side_tags=(1, 2))
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            sides...; face_tag=true)
        boolean_side = Any[(true, 0.0, 0.0), sides[1][2:end]...]
        @test_throws ArgumentError mesh_transfinite_triangle(
            boolean_side, sides[2], sides[3])
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_nodes=big(typemax(Int32)) + 1)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_nodes=9)
        @test_throws ArgumentError mesh_transfinite_triangle(
            sides...; max_triangles=8)
        bounded = mesh_transfinite_triangle(
            sides...; max_nodes=10, max_triangles=9)
        @test (nnodes(bounded), ntris(bounded)) == (10, 9)
        minimal = mesh_transfinite_triangle(_straight_triangle(1)...)
        @test (nnodes(minimal), nsegs(minimal), ntris(minimal)) == (3, 3, 1)
        @test validate(minimal).ok

        huge = _CountOnlyTriangleSide(typemax(Int))
        @test_throws ArgumentError mesh_transfinite_triangle(huge, huge, huge)
        wide = _CountOnlyTriangleSide(100_000)
        @test_throws ArgumentError mesh_transfinite_triangle(wide, wide, wide)
        @test_throws ArgumentError mesh_transfinite_triangle_patch(
            wide, wide, wide)
        @test _triangle_rejected_allocated(wide) < 64_000

        maximum = floatmax(Float64)
        overflowing = ([(maximum, 0.0, 0.0), (-maximum, 0.0, 0.0)],
                       [(-maximum, 0.0, 0.0), (-maximum, 1.0, 0.0)],
                       [(-maximum, 1.0, 0.0), (maximum, 0.0, 0.0)])
        @test_throws ArgumentError mesh_transfinite_triangle(overflowing...)
        @test isempty(Test.detect_ambiguities(
            Tessella.TransfiniteTriangle; recursive=true))
    end

    @testset "allocation growth remains linear in output size" begin
        small_sides = _straight_triangle(64)
        large_sides = _straight_triangle(128)
        mesh_transfinite_triangle(small_sides...; arrangement=:alternate_left)
        mesh_transfinite_triangle(large_sides...; arrangement=:alternate_left)
        small = _triangle_allocated(small_sides)
        large = _triangle_allocated(large_sides)
        @test small > 0
        @test large > small
        @test large <= 4.35small + 1_048_576
        @info "transfinite triangle allocation ratchet" small_bytes=small large_bytes=large

        mesh_transfinite_triangle_patch(
            small_sides...; arrangement=:alternate_left)
        mesh_transfinite_triangle_patch(
            large_sides...; arrangement=:alternate_left)
        patch_small = _triangle_patch_allocated(small_sides)
        patch_large = _triangle_patch_allocated(large_sides)
        @test patch_small > 0
        @test patch_large > patch_small
        @test patch_large <= 4.35patch_small + 1_048_576
        @info "recombined transfinite triangle allocation ratchet" small_bytes=patch_small large_bytes=patch_large
    end
end
