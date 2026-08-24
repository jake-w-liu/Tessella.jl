using Test
using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, mixed_crc
using Tessella.MeshTypes: tet_volume
using Tessella.Predicates: orient3

if !isdefined(Tessella, :TransfiniteHex)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "structured",
                          "TransfiniteHex.jl"))
end
using Tessella.TransfiniteHex: mesh_transfinite_hex

function _hex_affine_corners(origin=(0.0, 0.0, 0.0),
                             u=(2.0, 0.0, 0.0),
                             v=(0.0, 1.0, 0.0),
                             w=(0.0, 0.0, 1.0))
    add(vectors...) = ntuple(d -> sum(vector[d] for vector in vectors), 3)
    return [origin, add(origin, u), add(origin, u, v), add(origin, v),
            add(origin, w), add(origin, u, w), add(origin, u, v, w),
            add(origin, v, w)]
end

@inline _hex_node(i::Int, j::Int, k::Int, npu::Int, npv::Int) =
    Int32(i + 1 + npu * (j + npv * k))

function _expected_hexes(nu::Int, nv::Int, nw::Int)
    npu = nu + 1
    npv = nv + 1
    result = Matrix{Int32}(undef, 8, nu * nv * nw)
    cursor = 0
    for i in 0:nu-1, j in 0:nv-1, k in 0:nw-1
        cursor += 1
        result[:, cursor] .= (
            _hex_node(i, j, k, npu, npv),
            _hex_node(i + 1, j, k, npu, npv),
            _hex_node(i + 1, j + 1, k, npu, npv),
            _hex_node(i, j + 1, k, npu, npv),
            _hex_node(i, j, k + 1, npu, npv),
            _hex_node(i + 1, j, k + 1, npu, npv),
            _hex_node(i + 1, j + 1, k + 1, npu, npv),
            _hex_node(i, j + 1, k + 1, npu, npv))
    end
    return result
end

function _expected_boundary_quads(nu::Int, nv::Int, nw::Int)
    npu = nu + 1
    npv = nv + 1
    count = 2 * (nu * nv + nu * nw + nv * nw)
    result = Matrix{Int32}(undef, 4, count)
    cursor = 0
    write(nodes) = begin
        cursor += 1
        result[:, cursor] .= nodes
    end
    for i in 0:nu-1, k in 0:nw-1
        write((_hex_node(i, 0, k, npu, npv),
               _hex_node(i + 1, 0, k, npu, npv),
               _hex_node(i + 1, 0, k + 1, npu, npv),
               _hex_node(i, 0, k + 1, npu, npv)))
    end
    for j in 0:nv-1, k in 0:nw-1
        write((_hex_node(nu, j, k, npu, npv),
               _hex_node(nu, j + 1, k, npu, npv),
               _hex_node(nu, j + 1, k + 1, npu, npv),
               _hex_node(nu, j, k + 1, npu, npv)))
    end
    for i in 0:nu-1, k in 0:nw-1
        write((_hex_node(i, nv, k, npu, npv),
               _hex_node(i + 1, nv, k, npu, npv),
               _hex_node(i + 1, nv, k + 1, npu, npv),
               _hex_node(i, nv, k + 1, npu, npv)))
    end
    for j in 0:nv-1, k in 0:nw-1
        write((_hex_node(0, j, k, npu, npv),
               _hex_node(0, j + 1, k, npu, npv),
               _hex_node(0, j + 1, k + 1, npu, npv),
               _hex_node(0, j, k + 1, npu, npv)))
    end
    for i in 0:nu-1, j in 0:nv-1
        write((_hex_node(i, j, 0, npu, npv),
               _hex_node(i + 1, j, 0, npu, npv),
               _hex_node(i + 1, j + 1, 0, npu, npv),
               _hex_node(i, j + 1, 0, npu, npv)))
    end
    for i in 0:nu-1, j in 0:nv-1
        write((_hex_node(i, j, nw, npu, npv),
               _hex_node(i + 1, j, nw, npu, npv),
               _hex_node(i + 1, j + 1, nw, npu, npv),
               _hex_node(i, j + 1, nw, npu, npv)))
    end
    cursor == count || error("test boundary count invariant failed")
    return result
end

@inline function _hex_mesh_node(mesh, index::Integer)
    return (mesh.coords[1, index], mesh.coords[2, index],
            mesh.coords[3, index])
end

function _canonical_quads(matrix)
    result = NTuple{4,Int32}[]
    for cell in axes(matrix, 2)
        values = sort(matrix[:, cell])
        push!(result, (values[1], values[2], values[3], values[4]))
    end
    sort!(result)
    return result
end

function _hex_boundary(hexes)
    local_faces = ((1, 4, 3, 2), (1, 2, 6, 5), (1, 5, 8, 4),
                   (2, 3, 7, 6), (3, 4, 8, 7), (5, 6, 7, 8))
    incidence = Dict{NTuple{4,Int32},Int}()
    for cell in axes(hexes, 2), slots in local_faces
        values = sort(Int32[hexes[slots[1], cell], hexes[slots[2], cell],
                            hexes[slots[3], cell], hexes[slots[4], cell]])
        key = (values[1], values[2], values[3], values[4])
        incidence[key] = get(incidence, key, 0) + 1
    end
    boundary = sort!(collect(key for (key, count) in incidence if count == 1))
    return boundary, maximum(values(incidence); init=0)
end

function _hex_decomposed_volume(mesh, hexes)
    patterns = ((1, 2, 4, 5), (2, 4, 5, 6), (5, 6, 4, 8),
                (2, 4, 6, 3), (4, 8, 6, 3), (6, 8, 7, 3))
    volume = 0.0
    for cell in axes(hexes, 2), slots in patterns
        volume += tet_volume(
            _hex_mesh_node(mesh, hexes[slots[1], cell]),
            _hex_mesh_node(mesh, hexes[slots[2], cell]),
            _hex_mesh_node(mesh, hexes[slots[3], cell]),
            _hex_mesh_node(mesh, hexes[slots[4], cell]))
    end
    return volume
end

@inline _hex_exact_det(a, b, c) =
    a[1] * (b[2] * c[3] - b[3] * c[2]) -
    a[2] * (b[1] * c[3] - b[3] * c[1]) +
    a[3] * (b[1] * c[2] - b[2] * c[1])

function _hex_exact_integrated_volume(mesh, hexes)
    rational = Rational{BigInt}
    total = zero(rational)
    for cell in axes(hexes, 2)
        points = ntuple(8) do slot
            node = Int(hexes[slot, cell])
            ntuple(dimension -> rational(mesh.coords[dimension, node]), 3)
        end
        origin = points[1]
        subtract(a, b) = ntuple(dimension -> a[dimension] - b[dimension], 3)
        add(a, b) = ntuple(dimension -> a[dimension] + b[dimension], 3)
        deltas = ntuple(slot -> subtract(points[slot], origin), 8)
        u, v, w = deltas[2], deltas[4], deltas[5]
        uv = subtract(subtract(deltas[3], u), v)
        uw = subtract(subtract(deltas[6], u), w)
        vw = subtract(subtract(deltas[8], v), w)
        uvw = subtract(
            subtract(subtract(subtract(deltas[7], u), v), w),
            add(add(uv, uw), vw))
        du = (((0, 0, 0), u), ((0, 1, 0), uv),
              ((0, 0, 1), uw), ((0, 1, 1), uvw))
        dv = (((0, 0, 0), v), ((1, 0, 0), uv),
              ((0, 0, 1), vw), ((1, 0, 1), uvw))
        dw = (((0, 0, 0), w), ((1, 0, 0), uw),
              ((0, 1, 0), vw), ((1, 1, 0), uvw))
        coefficients = fill(zero(rational), 3, 3, 3)
        for (eu, xu) in du, (ev, xv) in dv, (ew, xw) in dw
            i = eu[1] + ev[1] + ew[1] + 1
            j = eu[2] + ev[2] + ew[2] + 1
            k = eu[3] + ev[3] + ew[3] + 1
            coefficients[i, j, k] += _hex_exact_det(xu, xv, xw)
        end
        # Integrate the power basis directly. This is intentionally independent
        # of the production Bernstein-average implementation.
        for k in 1:3, j in 1:3, i in 1:3
            total += coefficients[i, j, k] / (i * j * k)
        end
    end
    return total
end

function _hex_exact_affine_volume(corners)
    rational = Rational{BigInt}
    points = ntuple(4) do position
        point = corners[(1, 2, 4, 5)[position]]
        ntuple(dimension -> rational(point[dimension]), 3)
    end
    subtract(a, b) = ntuple(dimension -> a[dimension] - b[dimension], 3)
    return _hex_exact_det(subtract(points[2], points[1]),
                          subtract(points[3], points[1]),
                          subtract(points[4], points[1]))
end

function _hex_volume_loss_corners(power::Int; shrink=1.0)
    scale = Int64(2)^power
    u = (scale, scale, scale - 1)
    v = (scale + 1, scale + 1, scale)
    w = (scale + 1, scale, scale)
    scaled(vector) = ntuple(dimension -> shrink * vector[dimension], 3)
    return _hex_affine_corners((0.0, 0.0, 0.0),
                               scaled(u), scaled(v), scaled(w))
end

struct _UnreadHexCorners end
Base.length(::_UnreadHexCorners) = 8
Base.iterate(::_UnreadHexCorners) =
    error("corner conversion occurred before resource rejection")

@noinline function _hex_allocated(corners, cells)
    GC.gc()
    return @allocated mesh_transfinite_hex(corners, cells;
                                            arrangement=:alternate_right)
end

@noinline function _hex_rejected_allocated()
    GC.gc()
    return @allocated try
        mesh_transfinite_hex(
            _UnreadHexCorners(), (100_000, 100_000, 100_000))
    catch err
        err isa ArgumentError || rethrow()
    end
end

@testset "Gmsh recombined affine six-face transfinite hexahedra" begin
    @testset "type-5/type-3 order, tags, topology, CRC, and arrangements" begin
        corners = _hex_affine_corners()
        meshes = Dict{Symbol,MixedMesh}()
        for arrangement in
                (:left, :right, :alternate_left, :alternate_right)
            mesh = mesh_transfinite_hex(
                corners, (2, 2, 1); arrangement=arrangement, volume_tag=21,
                face_tags=(11, 12, 13, 14, 15, 16))
            meshes[arrangement] = mesh
            @test Tessella.Elements.validate(mesh).ok
            @test size(mesh.coords) == (3, 18)
            @test length(mesh.blocks) == 2
            quads, hexes = mesh.blocks
            @test quads isa ElementBlock
            @test hexes isa ElementBlock
            @test (quads.msh, hexes.msh) == (3, 5)
            @test size(quads.nodes) == (4, 16)
            @test size(hexes.nodes) == (8, 4)
            @test hexes.nodes == _expected_hexes(2, 2, 1)
            @test quads.nodes == _expected_boundary_quads(2, 2, 1)
            @test hexes.tags == fill(Int32(21), 4)
            @test quads.tags == Int32[
                11, 11, 12, 12, 13, 13, 14, 14,
                15, 15, 15, 15, 16, 16, 16, 16]
            boundary, maximum_incidence = _hex_boundary(hexes.nodes)
            @test maximum_incidence == 2
            @test boundary == _canonical_quads(quads.nodes)
            @test _hex_decomposed_volume(mesh, hexes.nodes) ≈ 2.0 rtol=64eps(Float64)
            @test mixed_crc(mesh) == mixed_crc(mesh_transfinite_hex(
                corners, (2, 2, 1); arrangement=arrangement, volume_tag=21,
                face_tags=(11, 12, 13, 14, 15, 16)))
        end
        @test Set(mixed_crc(mesh).sha for mesh in values(meshes)) == Set([
            "ed570464d1b553e8bb0850669656a0d9114be81bc516d2b6d4d6f6b8cf6a1515"])
        @test all(mesh.blocks[1].nodes == meshes[:left].blocks[1].nodes &&
                  mesh.blocks[2].nodes == meshes[:left].blocks[2].nodes
                  for mesh in values(meshes))

        one = mesh_transfinite_hex(corners)
        @test one.blocks[2].nodes[:, 1] == Int32[1, 2, 4, 3, 5, 6, 8, 7]
        @test one.blocks[1].nodes == _expected_boundary_quads(1, 1, 1)
    end

    @testset "affine interpolation, orientation, and represented volume" begin
        origin = (1.25, -2.5, 0.75)
        u = (2.0, 0.5, -0.25)
        v = (-0.4, 1.75, 0.3)
        w = (0.2, -0.35, 1.6)
        corners = _hex_affine_corners(origin, u, v, w)
        mesh = mesh_transfinite_hex(corners, (4, 3, 2))
        quads, hexes = mesh.blocks
        @test Tessella.Elements.validate(mesh).ok
        @test (size(mesh.coords, 2), size(quads.nodes, 2),
               size(hexes.nodes, 2)) == (60, 52, 24)
        node_id(i, j, k) = _hex_node(i, j, k, 5, 4)
        @test _hex_mesh_node(mesh, node_id(0, 0, 0)) == corners[1]
        @test _hex_mesh_node(mesh, node_id(4, 0, 0)) == corners[2]
        @test _hex_mesh_node(mesh, node_id(4, 3, 0)) == corners[3]
        @test _hex_mesh_node(mesh, node_id(0, 3, 0)) == corners[4]
        @test _hex_mesh_node(mesh, node_id(0, 0, 2)) == corners[5]
        @test _hex_mesh_node(mesh, node_id(4, 0, 2)) == corners[6]
        @test _hex_mesh_node(mesh, node_id(4, 3, 2)) == corners[7]
        @test _hex_mesh_node(mesh, node_id(0, 3, 2)) == corners[8]
        expected = ntuple(
            d -> origin[d] + 0.5u[d] + (1 / 3)v[d] + 0.5w[d], 3)
        actual = _hex_mesh_node(mesh, node_id(2, 1, 1))
        @test all(isapprox(actual[d], expected[d];
                           atol=32eps(Float64), rtol=32eps(Float64))
                  for d in 1:3)
        expected_volume = tet_volume(
            corners[1], corners[2], corners[4], corners[5]) * 6
        @test _hex_decomposed_volume(mesh, hexes.nodes) ≈ expected_volume rtol=256eps(Float64)

        corner_patterns = (
            (1, 2, 4, 5, -1), (2, 1, 3, 6, 1),
            (3, 4, 2, 7, -1), (4, 3, 1, 8, 1),
            (5, 6, 8, 1, 1), (6, 5, 7, 2, -1),
            (7, 8, 6, 3, 1), (8, 7, 5, 4, -1))
        @test all(begin
            nodes = hexes.nodes[:, cell]
            all(orient3(
                    _hex_mesh_node(mesh, nodes[pattern[1]]),
                    _hex_mesh_node(mesh, nodes[pattern[2]]),
                    _hex_mesh_node(mesh, nodes[pattern[3]]),
                    _hex_mesh_node(mesh, nodes[pattern[4]])) == pattern[5]
                for pattern in corner_patterns)
        end for cell in axes(hexes.nodes, 2))

        thin = _hex_affine_corners(
            (0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
            (1.0, 1.0 + 2.0^-31, 1.0),
            (1.0, 1.0, 1.0 + 2.0^-31))
        thin_mesh = mesh_transfinite_hex(thin, (2, 2, 2))
        @test Tessella.Elements.validate(thin_mesh).ok
        @test size(thin_mesh.blocks[2].nodes, 2) == 8
    end

    @testset "exact integrated-Jacobian volume conservation" begin
        adversarial = _hex_volume_loss_corners(24)
        expected = _hex_exact_affine_volume(adversarial)

        # Two and four layers happen to cancel this rounding pattern exactly.
        # They must remain valid, while a nonconserving even count still fails.
        for count in (2, 4)
            mesh = mesh_transfinite_hex(adversarial, (1, 1, count))
            @test _hex_exact_integrated_volume(mesh, mesh.blocks[2].nodes) ==
                  expected
        end
        for count in (3, 6)
            error = try
                mesh_transfinite_hex(adversarial, (1, 1, count))
                nothing
            catch err
                err
            end
            @test error isa ArgumentError
            @test occursin("exact trilinear-Jacobian integration",
                           sprint(showerror, error))
        end

        # The conservation audit is geometric, not tied to the axial loop.
        for cells in ((3, 1, 1), (1, 3, 1), (3, 3, 3))
            error = try
                mesh_transfinite_hex(adversarial, cells)
                nothing
            catch err
                err
            end
            @test error isa ArgumentError
            @test occursin("exact trilinear-Jacobian integration",
                           sprint(showerror, error))
        end

        # Neighboring binades lose or gain different represented volume; all
        # are detected rather than special-casing the original exponent.
        for power in (23, 25)
            error = try
                mesh_transfinite_hex(_hex_volume_loss_corners(power), (1, 1, 3))
                nothing
            catch err
                err
            end
            @test error isa ArgumentError
            @test occursin("exact trilinear-Jacobian integration",
                           sprint(showerror, error))
        end

        scaled_error = try
            mesh_transfinite_hex(
                _hex_volume_loss_corners(24; shrink=2.0^-500), (1, 1, 3))
            nothing
        catch err
            err
        end
        @test scaled_error isa ArgumentError
        @test occursin("exact trilinear-Jacobian integration",
                       sprint(showerror, scaled_error))

        # The affine and per-cell determinants are all below Float64's range.
        # Crossing a power-of-two carry in the balanced accumulator must not
        # create a cell-count-dependent acceptance cliff.
        delta = 2.0^-400
        underflow = _hex_affine_corners(
            (0.0, 0.0, 0.0), (delta, 0.0, 0.0),
            (0.0, delta, 0.0), (0.0, 0.0, delta))
        underflow_mesh = mesh_transfinite_hex(underflow, (1, 1, 1_025))
        @test Tessella.Elements.validate(underflow_mesh).ok
        @test (size(underflow_mesh.coords, 2),
               size(underflow_mesh.blocks[1].nodes, 2),
               size(underflow_mesh.blocks[2].nodes, 2)) ==
              (4_104, 4_102, 1_025)

        # Target the former subnormal-floor cliff without making the persistent
        # suite construct a 262,152-node mesh. The production geometry path is
        # exercised above; this isolates its balanced scaled accumulator at one
        # cell beyond the 65,536-ulp boundary.
        module_under_test = Tessella.TransfiniteHex
        count = 65_537
        tiny = Rational{BigInt}(1, big(1) << 1_200)
        expected_tiny = module_under_test._scaled_exact_positive(tiny)
        cell_tiny = module_under_test._scaled_exact_positive(tiny / count)
        bins = fill(module_under_test._Interval(0.0, 0.0), 64)
        occupied = falses(64)
        for cell in 1:count
            module_under_test._accumulate_scaled!(
                bins, occupied, cell_tiny, expected_tiny.exponent, cell)
        end
        @test isnothing(module_under_test._certify_integrated_volume(
            bins, occupied, expected_tiny))
    end

    @testset "validated blockers and pre-allocation resource limits" begin
        corners = _hex_affine_corners()
        @test_throws ArgumentError mesh_transfinite_hex(corners[1:7])
        @test_throws ArgumentError mesh_transfinite_hex([corners; corners[1]])
        short_point = Any[corners...]
        short_point[2] = (1.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_hex(short_point)
        extra_point = Any[(point..., index == 4 ? NaN : 0.0)
                          for (index, point) in pairs(corners)]
        @test_throws ArgumentError mesh_transfinite_hex(extra_point)
        nonfinite = copy(corners)
        nonfinite[7] = (NaN, 1.0, 1.0)
        @test_throws ArgumentError mesh_transfinite_hex(nonfinite)
        nonrepresentable = Any[corners...]
        nonrepresentable[2] = (big(10)^1000, 0, 0)
        @test_throws ArgumentError mesh_transfinite_hex(nonrepresentable)

        warped = copy(corners)
        warped[7] = (2.0, 1.0, 1.001)
        @test_throws ArgumentError mesh_transfinite_hex(warped)
        left_handed = _hex_affine_corners(
            (0.0, 0.0, 0.0), (0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0), (0.0, 0.0, 1.0))
        @test_throws ArgumentError mesh_transfinite_hex(left_handed)
        coplanar = _hex_affine_corners(
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0), (2.0, 0.0, 0.0))
        @test_throws ArgumentError mesh_transfinite_hex(coplanar)
        maximum = floatmax(Float64)
        overflowing_span = [(maximum, 0.0, 0.0), (-maximum, 0.0, 0.0),
                            (-maximum, 1.0, 0.0), (maximum, 1.0, 0.0),
                            (maximum, 0.0, 1.0), (-maximum, 0.0, 1.0),
                            (-maximum, 1.0, 1.0), (maximum, 1.0, 1.0)]
        @test_throws ArgumentError mesh_transfinite_hex(overflowing_span)
        overflowing_measure = _hex_affine_corners(
            (0.0, 0.0, 0.0), (1.0e150, 0.0, 0.0),
            (0.0, 1.0e150, 0.0), (0.0, 0.0, 1.0e150))
        @test_throws ArgumentError mesh_transfinite_hex(overflowing_measure)

        rounded_collapse = _hex_affine_corners(
            (2.0^53, 0.0, 0.0), (2.0, 0.0, 0.0),
            (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        @test_throws ArgumentError mesh_transfinite_hex(
            rounded_collapse, (2, 1, 1))

        @test_throws ArgumentError mesh_transfinite_hex(corners, (1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(corners, (1, 1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(corners, (0, 1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(corners, (-1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(corners, (true, 1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(corners, (1.0, 1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(
            corners, (big(typemax(Int32)) + 1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_hex(
            _UnreadHexCorners(),
            (typemax(Int32), typemax(Int32), typemax(Int32)))
        _hex_rejected_allocated()
        @test _hex_rejected_allocated() < 64_000

        @test_throws ArgumentError mesh_transfinite_hex(
            corners; arrangement=:alternate)
        @test_throws ArgumentError mesh_transfinite_hex(
            corners; arrangement="Left")
        @test_throws ArgumentError mesh_transfinite_hex(corners; volume_tag=true)
        @test_throws ArgumentError mesh_transfinite_hex(corners; volume_tag=-1)
        @test_throws ArgumentError mesh_transfinite_hex(
            corners; face_tags=(1, 2, 3, 4, 5))
        @test_throws ArgumentError mesh_transfinite_hex(
            corners; face_tags=(1, 2, 3, true, 5, 6))
        @test_throws ArgumentError mesh_transfinite_hex(
            corners, (2, 2, 1); max_nodes=17)
        @test_throws ArgumentError mesh_transfinite_hex(
            corners, (2, 2, 1); max_hexahedra=3)
        @test_throws ArgumentError mesh_transfinite_hex(
            corners, (2, 2, 1); max_boundary_quadrangles=15)
        @test_throws ArgumentError mesh_transfinite_hex(corners; max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_hex(
            corners; max_hexahedra=1.0)
        @test_throws ArgumentError mesh_transfinite_hex(
            corners; max_boundary_quadrangles=false)
        bounded = mesh_transfinite_hex(
            corners, (2, 2, 1); max_nodes=18, max_hexahedra=4,
            max_boundary_quadrangles=16)
        @test (size(bounded.coords, 2), size(bounded.blocks[1].nodes, 2),
               size(bounded.blocks[2].nodes, 2)) == (18, 16, 4)
        @test isempty(Test.detect_ambiguities(
            Tessella.TransfiniteHex; recursive=true))
    end

    @testset "allocation growth remains linear in output size" begin
        corners = _hex_affine_corners()
        small_cells = (12, 10, 8)
        large_cells = (24, 10, 8)
        mesh_transfinite_hex(corners, small_cells)
        mesh_transfinite_hex(corners, large_cells)
        small = _hex_allocated(corners, small_cells)
        large = _hex_allocated(corners, large_cells)
        @test small > 0
        @test large > small
        @test large <= 2.35small + 1_048_576
        @info "transfinite hex allocation ratchet" small_bytes=small large_bytes=large
    end
end
