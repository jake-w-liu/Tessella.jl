using Test
using Tessella
using Tessella.MeshTypes: boundary_faces, mesh_crc, nnodes, node, ntris, ntets,
                          tet_volume, validate
using Tessella.Predicates: orient3

if !isdefined(Tessella, :TransfiniteVolume)
    Base.include(Tessella, joinpath(
        @__DIR__, "..", "..", "src", "structured", "TransfiniteVolume.jl"))
end
using Tessella.TransfiniteVolume: mesh_transfinite_volume

function _affine_corners(origin=(0.0, 0.0, 0.0),
                         u=(2.0, 0.0, 0.0),
                         v=(0.0, 1.0, 0.0),
                         w=(0.0, 0.0, 1.0))
    add(vectors...) = ntuple(d -> sum(vector[d] for vector in vectors), 3)
    return [origin,
            add(origin, u),
            add(origin, u, v),
            add(origin, v),
            add(origin, w),
            add(origin, u, w),
            add(origin, u, v, w),
            add(origin, v, w)]
end

function _volume_canonical_tets(mesh)
    result = NTuple{4,Int32}[]
    for tet in axes(mesh.tets, 2)
        values = sort(mesh.tets[:, tet])
        push!(result, (values[1], values[2], values[3], values[4]))
    end
    sort!(result)
end

function _volume_canonical_triangles(matrix)
    result = NTuple{3,Int32}[]
    for triangle in axes(matrix, 2)
        values = sort(matrix[:, triangle])
        push!(result, (values[1], values[2], values[3]))
    end
    sort!(result)
end

function _transfinite_mesh_volume(mesh)
    sum(tet_volume(node(mesh, mesh.tets[1, tet]),
                   node(mesh, mesh.tets[2, tet]),
                   node(mesh, mesh.tets[3, tet]),
                   node(mesh, mesh.tets[4, tet]))
        for tet in axes(mesh.tets, 2); init=0.0)
end

struct _UnreadCorners end
Base.length(::_UnreadCorners) = 8
Base.iterate(::_UnreadCorners) = error("corner conversion occurred before resource rejection")

@noinline function _transfinite_volume_allocated(corners, cells)
    GC.gc()
    return @allocated mesh_transfinite_volume(corners, cells)
end

@noinline function _transfinite_volume_rejected_allocated()
    GC.gc()
    return @allocated try
        mesh_transfinite_volume(_UnreadCorners(), (100_000, 100_000, 100_000))
    catch err
        err isa ArgumentError || rethrow()
    end
end

@testset "six-face affine transfinite volumes" begin
    @testset "Gmsh six-tet subdivision, tags, boundary, and deterministic CRC" begin
        corners = _affine_corners()
        mesh = mesh_transfinite_volume(corners, (2, 2, 1);
                                       volume_tag=21,
                                       face_tags=(11, 12, 13, 14, 15, 16))
        @test validate(mesh).ok
        @test (nnodes(mesh), ntris(mesh), ntets(mesh)) == (18, 32, 24)
        @test mesh.tet_tag == fill(Int32(21), 24)
        @test count(==(Int32(11)), mesh.tri_tag) == 4
        @test count(==(Int32(12)), mesh.tri_tag) == 4
        @test count(==(Int32(13)), mesh.tri_tag) == 4
        @test count(==(Int32(14)), mesh.tri_tag) == 4
        @test count(==(Int32(15)), mesh.tri_tag) == 8
        @test count(==(Int32(16)), mesh.tri_tag) == 8
        boundary, maximum_incidence = boundary_faces(mesh.tets)
        @test maximum_incidence == 2
        @test sort!(boundary) == _volume_canonical_triangles(mesh.tris)
        center = ntuple(d -> sum(point[d] for point in corners) / 8, 3)
        @test all(orient3(node(mesh, mesh.tris[1, triangle]),
                          node(mesh, mesh.tris[2, triangle]),
                          node(mesh, mesh.tris[3, triangle]), center) > 0
                  for triangle in axes(mesh.tris, 2))
        @test mesh_crc(mesh).bbox == ((0.0, 0.0, 0.0), (2.0, 1.0, 1.0))
        @test mesh_crc(mesh).sha ==
              "6019ca07d5659879b2d2d0bb83590bae9b3e3344c29a75abfa86b28ab55b5d75"
        @test mesh_crc(mesh) == mesh_crc(mesh_transfinite_volume(
            corners, (2, 2, 1); volume_tag=21,
            face_tags=(11, 12, 13, 14, 15, 16)))

        one = mesh_transfinite_volume(_affine_corners(), (1, 1, 1))
        @test _volume_canonical_tets(one) == sort!(NTuple{4,Int32}[
            (1, 2, 3, 5), (2, 3, 5, 6), (3, 5, 6, 7),
            (2, 3, 4, 6), (3, 4, 6, 7), (4, 6, 7, 8)])
    end

    @testset "affine interpolation, exact corners, and volume conservation" begin
        origin = (1.25, -2.5, 0.75)
        u = (2.0, 0.5, -0.25)
        v = (-0.4, 1.75, 0.3)
        w = (0.2, -0.35, 1.6)
        corners = _affine_corners(origin, u, v, w)
        mesh = mesh_transfinite_volume(corners, (4, 3, 2))
        @test validate(mesh).ok
        @test (nnodes(mesh), ntris(mesh), ntets(mesh)) == (60, 104, 144)
        node_id(i, j, k) = i + 1 + 5 * (j + 4k)
        @test node(mesh, node_id(0, 0, 0)) == corners[1]
        @test node(mesh, node_id(4, 0, 0)) == corners[2]
        @test node(mesh, node_id(4, 3, 0)) == corners[3]
        @test node(mesh, node_id(0, 3, 0)) == corners[4]
        @test node(mesh, node_id(0, 0, 2)) == corners[5]
        @test node(mesh, node_id(4, 0, 2)) == corners[6]
        @test node(mesh, node_id(4, 3, 2)) == corners[7]
        @test node(mesh, node_id(0, 3, 2)) == corners[8]
        expected = ntuple(d -> origin[d] + 0.5u[d] + (1 / 3)v[d] + 0.5w[d], 3)
        @test all(isapprox(node(mesh, node_id(2, 1, 1))[d], expected[d];
                           atol=32eps(Float64), rtol=32eps(Float64)) for d in 1:3)
        expected_volume = 6tet_volume(corners[1], corners[2], corners[4], corners[5])
        @test _transfinite_mesh_volume(mesh) ≈ expected_volume rtol=256eps(Float64)

        # Exact dyadic affine identities remain valid even when independent
        # normalization rounds a derived corner one ULP away from its sum.
        conditioned = _affine_corners(
            (0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
            (1.0, 33 / 32, 1.0), (1.0, 1.0, 33 / 32))
        conditioned_mesh = mesh_transfinite_volume(conditioned)
        @test validate(conditioned_mesh).ok
        @test (nnodes(conditioned_mesh), ntris(conditioned_mesh),
               ntets(conditioned_mesh)) == (8, 12, 6)

        cancellation = _affine_corners(
            (0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
            (1.0, 1.0 + 2.0^-31, 1.0),
            (1.0, 1.0, 1.0 + 2.0^-31))
        cancellation_mesh = mesh_transfinite_volume(cancellation)
        @test validate(cancellation_mesh).ok
        @test ntets(cancellation_mesh) == 6

        center_cancellation = _affine_corners(
            (0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
            (1.0, 1.0 + 2.0^-51, 1.0),
            (1.0, 1.0, 1.0 + 2.0^-51))
        center_mesh = mesh_transfinite_volume(center_cancellation)
        @test validate(center_mesh).ok
        @test ntets(center_mesh) == 6
    end

    @testset "validated blockers and pre-allocation resource limits" begin
        corners = _affine_corners()
        @test_throws ArgumentError mesh_transfinite_volume(corners[1:7], (1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume([corners; [corners[1]]], (1, 1, 1))
        short_point = Any[corners...]; short_point[2] = (1.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_volume(short_point, (1, 1, 1))
        extra_point = Any[(point..., index == 4 ? NaN : 0.0)
                          for (index, point) in pairs(corners)]
        @test_throws ArgumentError mesh_transfinite_volume(extra_point, (1, 1, 1))
        nonfinite = copy(corners); nonfinite[7] = (NaN, 1.0, 1.0)
        @test_throws ArgumentError mesh_transfinite_volume(nonfinite, (1, 1, 1))
        nonrepresentable = Any[corners...]; nonrepresentable[2] = (big(10)^1000, 0, 0)
        @test_throws ArgumentError mesh_transfinite_volume(nonrepresentable, (1, 1, 1))

        warped = copy(corners); warped[7] = (2.0, 1.0, 1.001)
        @test_throws ArgumentError mesh_transfinite_volume(warped, (1, 1, 1))
        left_handed = _affine_corners((0., 0., 0.), (0., 1., 0.),
                                      (1., 0., 0.), (0., 0., 1.))
        @test_throws ArgumentError mesh_transfinite_volume(left_handed, (1, 1, 1))
        coplanar = _affine_corners((0., 0., 0.), (1., 0., 0.),
                                   (0., 1., 0.), (2., 0., 0.))
        @test_throws ArgumentError mesh_transfinite_volume(coplanar, (1, 1, 1))
        maximum = floatmax(Float64)
        overflowing = [(maximum, 0., 0.), (-maximum, 0., 0.),
                       (-maximum, 1., 0.), (maximum, 1., 0.),
                       (maximum, 0., 1.), (-maximum, 0., 1.),
                       (-maximum, 1., 1.), (maximum, 1., 1.)]
        @test_throws ArgumentError mesh_transfinite_volume(overflowing, (1, 1, 1))

        @test_throws ArgumentError mesh_transfinite_volume(corners, (1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(corners, (1, 1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(corners, (0, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(corners, (-1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(corners, (true, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(corners, (1.0, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, (big(typemax(Int32)) + 1, 1, 1))
        @test_throws ArgumentError mesh_transfinite_volume(
            _UnreadCorners(), (typemax(Int32), typemax(Int32), typemax(Int32)))
        _transfinite_volume_rejected_allocated()
        @test _transfinite_volume_rejected_allocated() < 64_000

        @test_throws ArgumentError mesh_transfinite_volume(corners, (2, 2, 1); max_nodes=17)
        @test_throws ArgumentError mesh_transfinite_volume(corners, (2, 2, 1); max_tets=23)
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, (2, 2, 1); max_boundary_triangles=31)
        bounded = mesh_transfinite_volume(
            corners, (2, 2, 1); max_nodes=18, max_tets=24,
            max_boundary_triangles=32)
        @test (nnodes(bounded), ntris(bounded), ntets(bounded)) == (18, 32, 24)
        @test_throws ArgumentError mesh_transfinite_volume(corners; max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_volume(corners; max_tets=false)
        @test_throws ArgumentError mesh_transfinite_volume(corners; max_nodes=-1)
        @test_throws ArgumentError mesh_transfinite_volume(
            corners; max_nodes=big(typemax(Int32)) + 1)
        @test_throws TypeError mesh_transfinite_volume(corners; max_nodes=18.0)

        @test_throws ArgumentError mesh_transfinite_volume(corners; volume_tag=true)
        @test_throws ArgumentError mesh_transfinite_volume(corners; volume_tag=-1)
        @test_throws ArgumentError mesh_transfinite_volume(
            corners; volume_tag=big(typemax(Int32)) + 1)
        @test_throws ArgumentError mesh_transfinite_volume(corners; face_tags=(1, 2, 3))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners; face_tags=(1, 2, 3, 4, 5, true))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners; face_tags=(1, 2, 3, 4, 5, -1))

        # The affine contract can still be geometrically unrepresentable at a
        # requested resolution; reject instead of emitting collapsed cells.
        thin = _affine_corners((1., 0., 0.), (eps(1.), 0., 0.),
                               (0., 1., 0.), (0., 0., 1.))
        @test_throws ArgumentError mesh_transfinite_volume(thin, (2, 1, 1))
        @test isempty(Test.detect_ambiguities(
            Tessella.TransfiniteVolume; recursive=true))
        @test isempty(Docs.undocumented_names(
            Tessella.TransfiniteVolume; private=false))
    end

    @testset "allocation growth remains linear in output size" begin
        corners = _affine_corners()
        mesh_transfinite_volume(corners, (8, 8, 4))
        mesh_transfinite_volume(corners, (16, 8, 4))
        small = _transfinite_volume_allocated(corners, (8, 8, 4))
        large = _transfinite_volume_allocated(corners, (16, 8, 4))
        @test small > 0
        @test large > small
        @test large <= 2.30small + 1_048_576
        @info "transfinite volume allocation ratchet" small_bytes=small large_bytes=large
    end
end
