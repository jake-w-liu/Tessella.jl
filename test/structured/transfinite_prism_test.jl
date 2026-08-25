using Test
using Tessella
using Tessella.MeshTypes: boundary_faces, mesh_crc, nnodes, node, ntris,
                          ntets, tet_volume, validate
using Tessella.Predicates: orient3

using Tessella.TransfinitePrism: mesh_transfinite_prism

function _affine_prism_corners(origin=(0.0, 0.0, 0.0),
                               u=(3.0, 0.0, 0.0),
                               v=(3.0, 2.0, 0.0),
                               w=(0.0, 0.0, 1.0))
    add(vectors...) = ntuple(
        dimension -> sum(vector[dimension] for vector in vectors), 3)
    return [origin, add(origin, u), add(origin, v),
            add(origin, w), add(origin, u, w), add(origin, v, w)]
end

@inline function _prism_node(nr::Int, ns::Int, i::Int, j::Int, k::Int)
    layer_nodes = 1 + nr * (ns + 1)
    return Int32(k * layer_nodes +
                 (i == 0 ? 1 : 2 + (i - 1) * (ns + 1) + j))
end

function _prism_canonical_tets(mesh)
    result = NTuple{4,Int32}[]
    for tet in axes(mesh.tets, 2)
        values = sort(mesh.tets[:, tet])
        push!(result, (values[1], values[2], values[3], values[4]))
    end
    sort!(result)
end

function _prism_canonical_triangles(matrix)
    result = NTuple{3,Int32}[]
    for triangle in axes(matrix, 2)
        values = sort(matrix[:, triangle])
        push!(result, (values[1], values[2], values[3]))
    end
    sort!(result)
end

function _prism_mesh_volume(mesh)
    return sum(tet_volume(node(mesh, mesh.tets[1, tet]),
                          node(mesh, mesh.tets[2, tet]),
                          node(mesh, mesh.tets[3, tet]),
                          node(mesh, mesh.tets[4, tet]))
               for tet in axes(mesh.tets, 2); init=0.0)
end

@inline function _prism_lerp3(a, b, t)
    return ((1 - t) * a[1] + t * b[1],
            (1 - t) * a[2] + t * b[2],
            (1 - t) * a[3] + t * b[3])
end

struct _UnreadPrismCorners end
Base.length(::_UnreadPrismCorners) = 6
Base.iterate(::_UnreadPrismCorners) =
    error("corner conversion occurred before resource rejection")

@noinline function _transfinite_prism_allocated(corners, cells)
    GC.gc()
    return @allocated mesh_transfinite_prism(corners, cells)
end

@noinline function _transfinite_prism_rejected_allocated()
    GC.gc()
    return @allocated try
        mesh_transfinite_prism(
            _UnreadPrismCorners(),
            (typemax(Int32), typemax(Int32), typemax(Int32)))
    catch err
        err isa ArgumentError || rethrow()
    end
end

@testset "five-face affine transfinite prisms" begin
    @test isdefined(Tessella, :TransfinitePrism)
    @test Tessella.mesh_transfinite_prism === mesh_transfinite_prism
    @test :mesh_transfinite_prism in names(Tessella)

    @testset "Gmsh collapsed-grid topology, tags, boundary, and CRC" begin
        corners = _affine_prism_corners()
        mesh = mesh_transfinite_prism(
            corners, (2, 3, 2); volume_tag=21,
            face_tags=(11, 12, 13, 14, 15))
        @test validate(mesh).ok
        @test (nnodes(mesh), ntris(mesh), ntets(mesh)) == (27, 46, 54)
        @test mesh.tet_tag == fill(Int32(21), 54)
        @test count(==(Int32(11)), mesh.tri_tag) == 8
        @test count(==(Int32(12)), mesh.tri_tag) == 12
        @test count(==(Int32(13)), mesh.tri_tag) == 8
        @test count(==(Int32(14)), mesh.tri_tag) == 9
        @test count(==(Int32(15)), mesh.tri_tag) == 9

        boundary, maximum_incidence = boundary_faces(mesh.tets)
        @test maximum_incidence == 2
        @test sort!(boundary) == _prism_canonical_triangles(mesh.tris)
        face_counts = (8, 12, 8, 9, 9)
        opposite = (corners[3], corners[1], corners[2], corners[4], corners[1])
        triangle = 0
        for face in 1:5, _ in 1:face_counts[face]
            triangle += 1
            @test orient3(node(mesh, mesh.tris[1, triangle]),
                          node(mesh, mesh.tris[2, triangle]),
                          node(mesh, mesh.tris[3, triangle]),
                          opposite[face]) > 0
        end
        @test triangle == ntris(mesh)
        @test mesh_crc(mesh).bbox == ((0.0, 0.0, 0.0), (3.0, 2.0, 1.0))
        @test mesh_crc(mesh).sha ==
              "a16de779890f62f8a09d928cbef67a6f13b09c6765a7d91ce8e86de78c14db6e"
        @test _prism_mesh_volume(mesh) == 3.0
        @test mesh_crc(mesh) == mesh_crc(mesh_transfinite_prism(
            corners, (2, 3, 2); volume_tag=21,
            face_tags=(11, 12, 13, 14, 15)))

        minimal = mesh_transfinite_prism(corners)
        @test (nnodes(minimal), ntris(minimal), ntets(minimal)) == (6, 8, 3)
        @test minimal.tets == Int32[1 2 4; 2 3 6; 3 4 5; 4 5 3]
        @test all(orient3(node(minimal, minimal.tets[1, tet]),
                          node(minimal, minimal.tets[2, tet]),
                          node(minimal, minimal.tets[3, tet]),
                          node(minimal, minimal.tets[4, tet])) == -1
                  for tet in axes(minimal.tets, 2))
        @test _prism_canonical_tets(minimal) == sort!(NTuple{4,Int32}[
            (1, 2, 3, 4), (2, 3, 4, 5), (3, 4, 5, 6)])
    end

    @testset "affine interpolation, exact corners, and volume conservation" begin
        origin = (1.25, -2.5, 0.75)
        u = (2.0, 0.5, -0.25)
        v = (-0.4, 1.75, 0.3)
        w = (0.2, -0.35, 1.6)
        corners = _affine_prism_corners(origin, u, v, w)
        nr, ns, nw = 4, 3, 2
        mesh = mesh_transfinite_prism(corners, (nr, ns, nw))
        @test validate(mesh).ok
        @test (nnodes(mesh), ntris(mesh), ntets(mesh)) == (51, 86, 126)

        @test node(mesh, _prism_node(nr, ns, 0, 0, 0)) == corners[1]
        @test node(mesh, _prism_node(nr, ns, nr, 0, 0)) == corners[2]
        @test node(mesh, _prism_node(nr, ns, nr, ns, 0)) == corners[3]
        @test node(mesh, _prism_node(nr, ns, 0, 0, nw)) == corners[4]
        @test node(mesh, _prism_node(nr, ns, nr, 0, nw)) == corners[5]
        @test node(mesh, _prism_node(nr, ns, nr, ns, nw)) == corners[6]

        radial = 2 / nr
        opposite = 1 / ns
        axial = 1 / nw
        lower = _prism_lerp3(
            corners[1], _prism_lerp3(corners[2], corners[3], opposite),
                       radial)
        upper = _prism_lerp3(
            corners[4], _prism_lerp3(corners[5], corners[6], opposite),
                       radial)
        expected = _prism_lerp3(lower, upper, axial)
        actual = node(mesh, _prism_node(nr, ns, 2, 1, 1))
        @test all(isapprox(actual[dimension], expected[dimension];
                           atol=32eps(Float64), rtol=32eps(Float64))
                  for dimension in 1:3)
        expected_volume = 3tet_volume(corners[1], corners[2],
                                      corners[3], corners[4])
        @test _prism_mesh_volume(mesh) ≈ expected_volume rtol=512eps(Float64)

        cancellation = _affine_prism_corners(
            (0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
            (1.0, 1.0 + 2.0^-31, 1.0),
            (1.0, 1.0, 1.0 + 2.0^-31))
        cancellation_mesh = mesh_transfinite_prism(cancellation)
        @test validate(cancellation_mesh).ok
        @test ntets(cancellation_mesh) == 3
    end

    @testset "validated blockers and pre-allocation resource limits" begin
        corners = _affine_prism_corners()
        @test_throws ArgumentError mesh_transfinite_prism(corners[1:5])
        @test_throws ArgumentError mesh_transfinite_prism([corners; [corners[1]]])
        short = Any[corners...]; short[2] = (3.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_prism(short)
        extra = Any[(point..., index == 3 ? NaN : 0.0)
                    for (index, point) in pairs(corners)]
        @test_throws ArgumentError mesh_transfinite_prism(extra)
        nonfinite = copy(corners); nonfinite[6] = (NaN, 2.0, 1.0)
        @test_throws ArgumentError mesh_transfinite_prism(nonfinite)
        nonrepresentable = Any[corners...]
        nonrepresentable[2] = (big(10)^1000, 0, 0)
        @test_throws ArgumentError mesh_transfinite_prism(nonrepresentable)

        warped = copy(corners); warped[6] = (3.0, 2.0, 1.001)
        @test_throws ArgumentError mesh_transfinite_prism(warped)
        left_handed = _affine_prism_corners(
            (0., 0., 0.), (0., 2., 0.), (3., 0., 0.), (0., 0., 1.))
        @test_throws ArgumentError mesh_transfinite_prism(left_handed)
        coplanar = _affine_prism_corners(
            (0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (2., 0., 0.))
        @test_throws ArgumentError mesh_transfinite_prism(coplanar)
        maximum = floatmax(Float64)
        overflowing = [(maximum, 0., 0.), (-maximum, 0., 0.),
                       (-maximum, 1., 0.), (maximum, 0., 1.),
                       (-maximum, 0., 1.), (-maximum, 1., 1.)]
        @test_throws ArgumentError mesh_transfinite_prism(overflowing)

        measure_base = ldexp(1.5, 551)
        measure_ulp = eps(measure_base)
        measure_point(offset) = ntuple(
            dimension -> measure_base + offset[dimension] * measure_ulp, 3)
        measure_u = (-2, -25, 58)
        measure_v = (39, -70, 43)
        measure_w = (-49, 98, -5)
        measure_add(a, b) = ntuple(
            dimension -> a[dimension] + b[dimension], 3)
        measure_overflow = measure_point.((
            (0, 0, 0), measure_u, measure_v, measure_w,
            measure_add(measure_u, measure_w),
            measure_add(measure_v, measure_w)))
        measure_error = try
            mesh_transfinite_prism(measure_overflow)
            nothing
        catch err
            err
        end
        @test measure_error isa ArgumentError
        @test occursin("must remain finite Float64 values",
                       sprint(showerror, measure_error))

        @test_throws ArgumentError mesh_transfinite_prism(corners, (1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(corners, (1, 1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(corners, (0, 1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(corners, (-1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(corners, (true, 1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(corners, (1.0, 1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (big(typemax(Int32)) + 1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_prism(
            _UnreadPrismCorners(),
            (typemax(Int32), typemax(Int32), typemax(Int32)))
        _transfinite_prism_rejected_allocated()
        @test _transfinite_prism_rejected_allocated() < 64_000

        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 3, 2); max_nodes=26)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 3, 2); max_tets=53)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 3, 2); max_boundary_triangles=45)
        bounded = mesh_transfinite_prism(
            corners, (2, 3, 2); max_nodes=27, max_tets=54,
            max_boundary_triangles=46)
        @test (nnodes(bounded), ntris(bounded), ntets(bounded)) == (27, 46, 54)
        @test_throws ArgumentError mesh_transfinite_prism(corners; max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_prism(corners; max_tets=false)
        @test_throws ArgumentError mesh_transfinite_prism(corners; max_nodes=-1)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners; max_nodes=big(typemax(Int32)) + 1)
        @test_throws ArgumentError mesh_transfinite_prism(corners; max_nodes=6.0)
        boolean_corner = Any[corners...]; boolean_corner[2] = (true, 0.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_prism(boolean_corner)

        @test_throws ArgumentError mesh_transfinite_prism(corners; volume_tag=true)
        @test_throws ArgumentError mesh_transfinite_prism(corners; volume_tag=-1)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners; volume_tag=big(typemax(Int32)) + 1)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners; face_tags=(1, 2, 3, 4))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners; face_tags=(1, 2, 3, 4, true))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners; face_tags=(1, 2, 3, 4, -1))

        thin = _affine_prism_corners(
            (1., 0., 0.), (eps(1.), 0., 0.),
            (eps(1.), 1., 0.), (0., 0., 1.))
        @test_throws ArgumentError mesh_transfinite_prism(thin, (2, 1, 1))

        # Exact affine corners can still become a folded represented grid when
        # interpolation rounds an ill-conditioned basis at intermediate axial
        # layers. Preserve Gmsh's one canonical tet order and reject the mixed
        # exact-predicate signs instead of independently flipping those tets.
        scale = Int64(2)^27
        u = (scale, scale, scale - 1)
        v = (scale + 1, scale + 1, scale)
        w = (scale + 1, scale, scale)
        add(a, b) = ntuple(dimension -> a[dimension] + b[dimension], 3)
        folded = [(0, 0, 0), u, v, w, add(u, w), add(v, w)]
        folded_error = try
            mesh_transfinite_prism(folded, (1, 1, 3))
            nothing
        catch err
            err
        end
        @test folded_error isa ArgumentError
        @test occursin("represented grid is folded", sprint(showerror, folded_error))

        volume_scale = Int64(2)^24
        volume_u = (volume_scale, volume_scale, volume_scale - 1)
        volume_v = (volume_scale + 1, volume_scale + 1, volume_scale)
        volume_w = (volume_scale + 1, volume_scale, volume_scale)
        volume_loss = [(0, 0, 0), volume_u, volume_v, volume_w,
                       add(volume_u, volume_w), add(volume_v, volume_w)]
        volume_error = try
            mesh_transfinite_prism(volume_loss, (1, 1, 3))
            nothing
        catch err
            err
        end
        @test volume_error isa ArgumentError
        @test occursin("do not conserve the affine prism volume",
                       sprint(showerror, volume_error))

        shrink = 2.0^-500
        scaled_u = ntuple(dimension -> shrink * volume_u[dimension], 3)
        scaled_v = ntuple(dimension -> shrink * volume_v[dimension], 3)
        scaled_w = ntuple(dimension -> shrink * volume_w[dimension], 3)
        scaled_loss = [(0.0, 0.0, 0.0), scaled_u, scaled_v, scaled_w,
                       add(scaled_u, scaled_w), add(scaled_v, scaled_w)]
        scaled_error = try
            mesh_transfinite_prism(scaled_loss, (1, 1, 3))
            nothing
        catch err
            err
        end
        @test scaled_error isa ArgumentError
        @test occursin("do not conserve the affine prism volume",
                       sprint(showerror, scaled_error))

        # Relative determinant accumulation must remain valid below the
        # Float64 volume range. Per-tet Float64 volume floors would make this
        # otherwise-valid result depend on whether nw is below or above 21_846.
        delta = 2.0^-1000
        underflow = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                     (0.0, delta, 0.0), (0.0, 0.0, delta),
                     (1.0, 0.0, delta), (0.0, delta, delta)]
        underflow_mesh = mesh_transfinite_prism(underflow, (1, 1, 21_847))
        @test validate(underflow_mesh).ok
        @test (nnodes(underflow_mesh), ntris(underflow_mesh),
               ntets(underflow_mesh)) == (65_544, 131_084, 65_541)

        # The exact fallback's leading numerator bits must convert to Float64
        # without rounding 2^54-1 up to the next binade.
        coefficient_a = 786_429
        coefficient_b = 262_657
        coefficient_c = 87_211
        large = Int64(2)^50
        exact_u = (coefficient_a, 0, 0)
        exact_v = (large, coefficient_b, 0)
        exact_w = (large, large, coefficient_c)
        exact_measure = [(0, 0, 0), exact_u, exact_v, exact_w,
                         add(exact_u, exact_w), add(exact_v, exact_w)]
        exact_measure_mesh = mesh_transfinite_prism(exact_measure)
        @test validate(exact_measure_mesh).ok
        @test ntets(exact_measure_mesh) == 3
        @test isempty(Test.detect_ambiguities(
            Tessella.TransfinitePrism; recursive=true))
        @test isempty(Docs.undocumented_names(
            Tessella.TransfinitePrism; private=false))
    end

    @testset "allocation growth remains linear in output size" begin
        corners = _affine_prism_corners()
        mesh_transfinite_prism(corners, (24, 12, 6))
        mesh_transfinite_prism(corners, (48, 12, 6))
        small = _transfinite_prism_allocated(corners, (24, 12, 6))
        large = _transfinite_prism_allocated(corners, (48, 12, 6))
        @test small > 0
        @test large > small
        @test large <= 2.30small + 1_048_576
        @info "transfinite prism allocation ratchet" small_bytes=small large_bytes=large
    end
end
