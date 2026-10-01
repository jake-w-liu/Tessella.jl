using Test
using Tessella
using Tessella.MeshTypes: boundary_faces, mesh_crc, nnodes, node, ntris, ntets,
                          tet_volume, validate
using Tessella.Predicates: orient3

if !isdefined(Tessella, :TransfiniteVolume)
    Base.include(Tessella, joinpath(
        @__DIR__, "..", "..", "src", "structured", "TransfiniteVolume.jl"))
end
using Tessella.Transfinite: mesh_transfinite_patch
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

        # Nested Float64 interpolation lost about 2% of this represented affine
        # volume near a large translation. Exact dyadic interpolation must
        # preserve the corner-defined material measure.
        base = 1e100
        ulp = eps(base)
        point(offset) = ntuple(d -> base + offset[d] * ulp, 3)
        add_offsets(a, b) = ntuple(d -> a[d] + b[d], 3)
        ru = (14, 4, -1); rv = (1, 30, 2); rw = (4, 2, 14)
        zero_offset = (0, 0, 0)
        remote = point.((zero_offset, ru, add_offsets(ru, rv), rv, rw,
                         add_offsets(ru, rw),
                         add_offsets(add_offsets(ru, rv), rw),
                         add_offsets(rv, rw)))
        remote_mesh = mesh_transfinite_volume(remote, (4, 4, 4))
        R = Rational{BigInt}
        exact_u = R.(remote[2]) .- R.(remote[1])
        exact_v = R.(remote[4]) .- R.(remote[1])
        exact_w = R.(remote[5]) .- R.(remote[1])
        exact_cross = (exact_v[2] * exact_w[3] - exact_v[3] * exact_w[2],
                       exact_v[3] * exact_w[1] - exact_v[1] * exact_w[3],
                       exact_v[1] * exact_w[2] - exact_v[2] * exact_w[1])
        exact_volume = Float64(abs(sum(exact_u .* exact_cross)))
        @test validate(remote_mesh).ok
        @test _transfinite_mesh_volume(remote_mesh) ≈ exact_volume rtol=8eps(Float64)
        remote_corner_nodes = (1, 5, 25, 21, 101, 105, 125, 121)
        @test all(node(remote_mesh, remote_corner_nodes[index]) == remote[index]
                  for index in 1:8)
        @test mesh_crc(remote_mesh).sha ==
              "cfdebd9e1af30eb255ed966e95bc3999f89d8062872d5cdb370647aa1737dfa8"

        # The exact interpolation path is reserved for exactly affine input.
        # A tolerated one-ULP corner residual must not be silently discarded.
        wide_offset = Int64(4_000_000_000_000)
        remote_origin = (base, base, base)
        remote_u = (wide_offset * ulp, 0.0, 0.0)
        remote_v = (0.0, wide_offset * ulp, 0.0)
        remote_w = (0.0, 0.0, wide_offset * ulp)
        almost_affine = _affine_corners(
            remote_origin, remote_u, remote_v, remote_w)
        almost_affine[3] = (nextfloat(almost_affine[3][1]),
                            almost_affine[3][2], almost_affine[3][3])
        almost_mesh = mesh_transfinite_volume(almost_affine, (4, 4, 4))
        almost_node_id(i, j, k) = i + 1 + 5 * (j + 5k)
        @test node(almost_mesh, almost_node_id(4, 4, 0)) == almost_affine[3]
        @test validate(almost_mesh).ok
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
        boolean_corner = Any[corners...]; boolean_corner[2] = (true, 0.0, 0.0)
        @test_throws ArgumentError mesh_transfinite_volume(boolean_corner, (1, 1, 1))

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
        @test_throws ArgumentError mesh_transfinite_volume(corners; max_nodes=18.0)

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

        # Nonzero exact predicates are insufficient when rounded grid nodes
        # fabricate material. The global determinant certificate must block it.
        remote_base = 1e100
        remote_ulp = eps(remote_base)
        remote_point(offset) = ntuple(
            dimension -> remote_base + offset[dimension] * remote_ulp, 3)
        add_offsets(a, b) = ntuple(dimension -> a[dimension] + b[dimension], 3)
        loss_u = (-24, 20, -5)
        loss_v = (69, 66, 34)
        loss_w = (-63, 70, -46)
        zero_offset = (0, 0, 0)
        material_loss = remote_point.((
            zero_offset, loss_u, add_offsets(loss_u, loss_v), loss_v, loss_w,
            add_offsets(loss_u, loss_w),
            add_offsets(add_offsets(loss_u, loss_v), loss_w),
            add_offsets(loss_v, loss_w)))
        material_error = try
            mesh_transfinite_volume(material_loss, (3, 4, 2))
            nothing
        catch err
            err
        end
        @test material_error isa ArgumentError
        @test occursin("do not conserve the affine block volume",
                       sprint(showerror, material_error))

        measure_base = ldexp(1.5, 570)
        measure_ulp = eps(measure_base)
        measure_point(offset) = ntuple(
            dimension -> measure_base + offset[dimension] * measure_ulp, 3)
        measure_u = (40, 30, 23)
        measure_v = (94, -95, -13)
        measure_w = (64, 93, -67)
        measure_overflow = measure_point.((
            zero_offset, measure_u, add_offsets(measure_u, measure_v), measure_v,
            measure_w, add_offsets(measure_u, measure_w),
            add_offsets(add_offsets(measure_u, measure_v), measure_w),
            add_offsets(measure_v, measure_w)))
        measure_error = try
            mesh_transfinite_volume(measure_overflow, (3, 1, 5))
            nothing
        catch err
            err
        end
        @test measure_error isa ArgumentError
        @test occursin("must remain finite Float64 values",
                       sprint(showerror, measure_error))
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

# Six canonical face-slot grids for `mesh_transfinite_volume(; faces=...)`.
# Slot order is Gmsh's `(vmin, umax, vmax, umin, wmin, wmax)`; each slot lists
# its four volume corners in canonical parametric order and the two logical
# axes (1=u, 2=v, 3=w) its (p1, p2) grid spans.
const _VOL_FACE_CORNERS =
    ((1, 2, 6, 5), (2, 3, 7, 6), (4, 3, 7, 8),
     (1, 4, 8, 5), (1, 2, 3, 4), (5, 6, 7, 8))
const _VOL_FACE_AXES =
    ((1, 3), (2, 3), (1, 3), (2, 3), (1, 2), (1, 2))
const _VOL_FACE_FIXED_AXIS = (2, 1, 2, 1, 3, 3)
const _VOL_EDGE_AXIS =
    Dict((1, 2) => 1, (2, 3) => 2, (3, 4) => 1, (1, 4) => 2,
         (5, 6) => 1, (6, 7) => 2, (7, 8) => 1, (5, 8) => 2,
         (1, 5) => 3, (2, 6) => 3, (3, 7) => 3, (4, 8) => 3)

# Shared boundary edges must agree bitwise between the two owning faces, so
# each of the twelve edges is discretized once in one canonical direction and
# the reversed *vector* is handed to the face going the other way — the same
# convention the model path applies to signed curves. Linear interpolation in
# the opposite direction would differ by ulps and trip the certificate.
function _volume_face_grids(corners, cells)
    linear(a, b, n) = [a .+ t .* (b .- a) for t in range(0.0, 1.0, n + 1)]
    cache = Dict{Tuple{Int,Int},Vector{NTuple{3,Float64}}}()
    function edge(a, b)
        return get!(cache, a < b ? (a, b) : (b, a)) do
            a < b ? linear(corners[a], corners[b],
                           cells[_VOL_EDGE_AXIS[(a, b)]]) :
                    linear(corners[b], corners[a],
                           cells[_VOL_EDGE_AXIS[(b, a)]])
        end
    end
    side(a, b) = a < b ? edge(a, b) : reverse(edge(a, b))
    return ntuple(6) do slot
        cs = _VOL_FACE_CORNERS[slot]
        patch = mesh_transfinite_patch(side(cs[1], cs[2]), side(cs[2], cs[3]),
                                       side(cs[3], cs[4]), side(cs[4], cs[1]);
                                       allow_warped=true)
        (patch.coords, patch.tris, patch.tri_tag)
    end
end

# Volume tab index of a canonical slot node: `(p1, p2)` on `slot` maps to
# logical `(i, j, k)` with the slot's two parametric axes taking p1/p2 and the
# third axis pinned at 0 or its maximum.
function _slot_node_id(slot, p1, p2, cells)
    nu, nv, nw = cells
    a1, a2 = _VOL_FACE_AXES[slot]
    fixed_axis = _VOL_FACE_FIXED_AXIS[slot]
    ijk = zeros(Int, 3)
    ijk[a1] = p1
    ijk[a2] = p2
    ijk[fixed_axis] = slot in (2, 3, 6) ? cells[fixed_axis] : 0
    return ijk[1] + 1 + (nu + 1) * (ijk[2] + (nv + 1) * ijk[3])
end

@testset "six canonical face grids (transfiniteHex path)" begin
    @testset "affine face grids reproduce the affine mesh" begin
        corners = _affine_corners()
        cells = (4, 3, 2)
        faces = _volume_face_grids(corners, cells)
        warped = mesh_transfinite_volume(corners, cells; faces=faces,
                                         volume_tag=21,
                                         face_tags=(11, 12, 13, 14, 15, 16))
        affine = mesh_transfinite_volume(corners, cells)
        @test validate(warped).ok
        @test (nnodes(warped), ntris(warped), ntets(warped)) ==
              (nnodes(affine), ntris(affine), ntets(affine))
        # Boundary nodes reuse the face grids bitwise; interior nodes agree
        # with the affine map to chord-ratio rounding.
        @test maximum(j -> maximum(abs, warped.coords[:, j] .-
                                        affine.coords[:, j]),
                      axes(warped.coords, 2)) <= 1e-14
        @test _volume_canonical_tets(warped) == _volume_canonical_tets(affine)
        @test _volume_canonical_triangles(warped.tris) ==
              _volume_canonical_triangles(affine.tris)
        @test warped.tet_tag == fill(Int32(21), ntets(warped))
        nu, nv, nw = cells
        @test count(==(Int32(11)), warped.tri_tag) == 2 * nu * nw
        @test count(==(Int32(12)), warped.tri_tag) == 2 * nv * nw
        @test count(==(Int32(13)), warped.tri_tag) == 2 * nu * nw
        @test count(==(Int32(14)), warped.tri_tag) == 2 * nv * nw
        @test count(==(Int32(15)), warped.tri_tag) == 2 * nu * nv
        @test count(==(Int32(16)), warped.tri_tag) == 2 * nu * nv
        for slot in 1:6
            a1, a2 = _VOL_FACE_AXES[slot]
            n1 = cells[a1]; n2 = cells[a2]
            for p2 in 0:n2, p1 in 0:n1
                id = _slot_node_id(slot, p1, p2, cells)
                @test Tuple(warped.coords[:, id]) ==
                      Tuple(faces[slot][1][:, p1 + 1 + p2 * (n1 + 1)])
            end
        end
        center = ntuple(d -> sum(point[d] for point in corners) / 8, 3)
        @test all(orient3(node(warped, warped.tris[1, triangle]),
                          node(warped, warped.tris[2, triangle]),
                          node(warped, warped.tris[3, triangle]), center) > 0
                  for triangle in axes(warped.tris, 2))
    end

    @testset "non-affine corners interpolate the face Coons interior" begin
        corners = _affine_corners((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                  (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        corners[7] = (1.2, 1.0, 1.0)
        @test_throws ArgumentError mesh_transfinite_volume(corners, (2, 2, 2))
        faces = _volume_face_grids(corners, (2, 2, 2))
        mesh = mesh_transfinite_volume(corners, (2, 2, 2); faces=faces)
        @test validate(mesh).ok
        @test (nnodes(mesh), ntris(mesh), ntets(mesh)) == (27, 48, 48)
        nid(i, j, k) = i + 1 + 3 * (j + 3k)
        @test all(isapprox(node(mesh, nid(1, 1, 1))[d],
                           (0.525, 0.5, 0.5)[d]; atol=8eps(Float64))
                  for d in 1:3)
        # Boundary nodes are the supplied face vertices bitwise.
        for slot in 1:6
            a1, a2 = _VOL_FACE_AXES[slot]
            n1, n2 = (2, 2, 2)[a1], (2, 2, 2)[a2]
            for p2 in 0:n2, p1 in 0:n1
                id = _slot_node_id(slot, p1, p2, (2, 2, 2))
                @test Tuple(mesh.coords[:, id]) ==
                      Tuple(faces[slot][1][:, p1 + 1 + p2 * (n1 + 1)])
            end
        end
        boundary, maximum_incidence = boundary_faces(mesh.tets)
        @test maximum_incidence == 2
        @test sort!(boundary) == _volume_canonical_triangles(mesh.tris)
    end

    @testset "curved boundary edges reach the interior" begin
        corners = _affine_corners((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                  (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        cells = (2, 2, 2)
        faces = collect(_volume_face_grids(corners, cells))
        # Bow edge s0s1 out of the v-min and w-min grids' shared p2=0 chain.
        # The endpoints are welded to the corner vertices, like production.
        for i in 0:2
            bowed = i == 0 ? corners[1] :
                    i == 2 ? corners[2] :
                    (0.5, 0.3 * sin(pi * 0.5i), 0.0)
            for slot in (1, 5)
                points = faces[slot][1]
                points[1, i + 1] = bowed[1]
                points[2, i + 1] = bowed[2]
                points[3, i + 1] = bowed[3]
            end
        end
        mesh = mesh_transfinite_volume(corners, cells; faces=Tuple(faces))
        @test validate(mesh).ok
        nid(i, j, k) = i + 1 + 3 * (j + 3k)
        @test node(mesh, nid(1, 0, 0)) == (0.5, 0.3, 0.0)
        # The opposite face interior is untouched by the single bowed edge.
        @test node(mesh, nid(1, 1, 0)) == (0.5, 0.5, 0.0)
        # transfiniteHex pulls the center toward the bow.
        @test all(isapprox(node(mesh, nid(1, 1, 1))[d],
                           (0.5, 0.425, 0.5)[d]; atol=16eps(Float64))
                  for d in 1:3)
    end

    @testset "face-grid and certificate validation" begin
        corners = _affine_corners()
        cells = (2, 2, 1)
        faces = _volume_face_grids(corners, cells)
        @test validate(mesh_transfinite_volume(
            corners, cells; faces=faces)).ok

        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=faces[1:5])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=(faces..., faces[1]))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=(faces[1], faces[1], faces[1],
                                   faces[1], faces[1], faces[1]))

        short = collect(faces)
        short[1] = (short[1][1][:, 1:end-1], short[1][2], short[1][3])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(short))

        badshape = collect(faces)
        badshape[1] = (badshape[1][1][1:2, :], badshape[1][2],
                       badshape[1][3])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(badshape))

        nonfinite = collect(faces)
        points = copy(nonfinite[1][1])
        points[1, 2] = NaN
        nonfinite[1] = (points, nonfinite[1][2], nonfinite[1][3])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(nonfinite))

        shorttris = collect(faces)
        shorttris[1] = (shorttris[1][1], shorttris[1][2][:, 1:end-1],
                        shorttris[1][3])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(shorttris))

        badindex = collect(faces)
        tris = copy(badindex[1][2])
        tris[1, 1] = size(badindex[1][1], 2) + 1
        badindex[1] = (badindex[1][1], tris, badindex[1][3])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(badindex))

        shorttags = collect(faces)
        shorttags[1] = (shorttags[1][1], shorttags[1][2],
                        shorttags[1][3][1:end-1])
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(shorttags))

        badtag = collect(faces)
        tags = copy(badtag[1][3])
        tags[1] = -1
        badtag[1] = (badtag[1][1], badtag[1][2], tags)
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=Tuple(badtag))

        # Perturb one interior shared-edge node in a single face: the two
        # owners of s0s1 must agree bitwise.
        drifted = collect(faces)
        points = copy(drifted[5][1])
        points[1, 2] = nextfloat(points[1, 2])
        drifted[5] = (points, drifted[5][2], drifted[5][3])
        error = try
            mesh_transfinite_volume(corners, cells; faces=Tuple(drifted))
            nothing
        catch err
            err
        end
        @test error isa ArgumentError
        @test occursin("shared edge s0s1", sprint(showerror, error))

        moved = copy(corners)
        moved[7] = moved[7] .+ (0.0, 0.0, 1e-9)
        error = try
            mesh_transfinite_volume(moved, cells; faces=faces)
            nothing
        catch err
            err
        end
        @test error isa ArgumentError
        @test occursin("corners do not match the face grids",
                       sprint(showerror, error))

        # A left-handed corner order cannot pass the corner-matching gate
        # against right-handed faces.
        left = [corners[1], corners[4], corners[3], corners[2],
                corners[5], corners[8], corners[7], corners[6]]
        @test_throws ArgumentError mesh_transfinite_volume(
            left, cells; faces=faces)

        collapsed = collect(faces)
        for slot in (1, 5)
            points = copy(collapsed[slot][1])
            points[:, 2] = points[:, 1]
            collapsed[slot] = (points, collapsed[slot][2],
                               collapsed[slot][3])
        end
        flat = copy(corners)
        flat[2] = flat[1]
        @test_throws ArgumentError mesh_transfinite_volume(
            flat, cells; faces=Tuple(collapsed))

        mismatched = _volume_face_grids(corners, (3, 2, 1))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; faces=mismatched)
    end
end

# Simplex decompositions of the emitted recombined cells, indexed in each
# cell tuple's own vertex order — the templates the kernel certifies
# internally; kept local so the partition check below is test-side.
const _VOL_SHADOW_HEX =
    ((1, 2, 4, 5), (2, 4, 5, 6), (5, 6, 4, 8),
     (2, 4, 6, 3), (4, 8, 6, 3), (6, 8, 7, 3))
const _VOL_SHADOW_PRISM_1 = ((1, 2, 3, 4), (2, 3, 4, 5), (4, 5, 3, 6))
const _VOL_SHADOW_PRISM_2 = ((3, 2, 6, 1), (2, 5, 6, 1), (6, 5, 4, 1))
const _VOL_SHADOW_PRISM_1R = ((2, 1, 5, 3), (1, 3, 4, 5), (5, 4, 6, 3))
const _VOL_SHADOW_PRISM_2R = ((3, 1, 2, 6), (1, 4, 2, 6), (6, 4, 2, 5))

function _vol_block(mesh, msh)
    found = filter(block -> block.msh == Int32(msh), mesh.blocks)
    length(found) <= 1 || error("expected at most one type-$msh block")
    return isempty(found) ? nothing : found[1]
end

# `==` on block/mesh structs is identity — compare fields explicitly.
function _vol_same_blocks(a, b)
    length(a.blocks) == length(b.blocks) || return false
    return all((ba.msh, ba.nodes, ba.tags) == (bb.msh, bb.nodes, bb.tags)
               for (ba, bb) in zip(a.blocks, b.blocks))
end

# Canonical per-block cell multisets: insensitive to emission order and
# vertex winding, so the affine and warped emitters can be compared
# topologically.
function _vol_canonical_blocks(mesh)
    return sort!(map(mesh.blocks) do block
        (block.msh,
         sort!([Tuple(sort(block.nodes[:, cell]))
                for cell in axes(block.nodes, 2)]),
         sort(block.tags))
    end)
end

@noinline function _volume_recombined_allocated(corners, cells)
    GC.gc()
    return @allocated mesh_transfinite_volume(corners, cells; recombine=true)
end

@testset "six-face recombined volumes (MixedMesh)" begin
    corners = _affine_corners()
    cells = (2, 2, 1)
    nu, nv, nw = cells
    face_cells = (nu * nw, nv * nw, nu * nw, nv * nw, nu * nv, nu * nv)
    reference = mesh_transfinite_volume(corners, cells)
    center = ntuple(d -> sum(point[d] for point in corners) / 8, 3)

    @testset "all-six recombined -> hexahedra" begin
        mesh = mesh_transfinite_volume(
            corners, cells; recombine=true, volume_tag=21,
            face_tags=(11, 12, 13, 14, 15, 16))
        @test validate(mesh).ok
        @test mesh isa Tessella.Elements.MixedMesh
        @test length(mesh.blocks) == 2
        hexes = _vol_block(mesh, 5)
        quads = _vol_block(mesh, 3)
        @test _vol_block(mesh, 2) === nothing
        @test _vol_block(mesh, 6) === nothing
        @test size(hexes.nodes, 2) == nu * nv * nw
        @test hexes.tags == fill(Int32(21), nu * nv * nw)
        @test size(quads.nodes, 2) == sum(face_cells)
        for slot in 1:6
            @test count(==(Int32(10 + slot)), quads.tags) == face_cells[slot]
        end
        # Same lattice as the simplex mesh, bitwise.
        @test mesh.coords == reference.coords
        # Ordered cell tuples are the Gmsh CREATE_HEX order.
        nid(i, j, k) = Int32(i + 1 + 3 * (j + 3k))
        @test hexes.nodes[:, 1] ==
              Int32[nid(0, 0, 0), nid(1, 0, 0), nid(1, 1, 0), nid(0, 1, 0),
                    nid(0, 0, 1), nid(1, 0, 1), nid(1, 1, 1), nid(0, 1, 1)]
        # Boundary quads cover each face's cells exactly once.
        for slot in 1:6
            a1, a2 = _VOL_FACE_AXES[slot]
            expected = NTuple{4,Int32}[]
            for p2 in 0:cells[a2]-1, p1 in 0:cells[a1]-1
                ids = sort(Int32[_slot_node_id(slot, p1, p2, cells),
                                 _slot_node_id(slot, p1 + 1, p2, cells),
                                 _slot_node_id(slot, p1 + 1, p2 + 1, cells),
                                 _slot_node_id(slot, p1, p2 + 1, cells)])
                push!(expected, Tuple(ids))
            end
            emitted = NTuple{4,Int32}[]
            for cell in axes(quads.nodes, 2)
                quads.tags[cell] == Int32(10 + slot) || continue
                push!(emitted, Tuple(sort(quads.nodes[:, cell])))
            end
            @test sort!(emitted) == sort!(expected)
        end
        # Every boundary quad is strictly outward on the affine box.
        @inbounds for cell in axes(quads.nodes, 2)
            a = Tuple(mesh.coords[:, quads.nodes[1, cell]])
            b = Tuple(mesh.coords[:, quads.nodes[2, cell]])
            c = Tuple(mesh.coords[:, quads.nodes[3, cell]])
            d = Tuple(mesh.coords[:, quads.nodes[4, cell]])
            @test orient3(a, b, c, center) > 0
            @test orient3(a, c, d, center) > 0
        end
        # The shadow partition is exactly the reference six-tet subdivision.
        shadow = NTuple{4,Int32}[]
        for cell in axes(hexes.nodes, 2)
            verts = ntuple(r -> hexes.nodes[r, cell], 8)
            for tet in _VOL_SHADOW_HEX
                push!(shadow, Tuple(sort!(
                    Int32[verts[tet[1]], verts[tet[2]],
                          verts[tet[3]], verts[tet[4]]])))
            end
        end
        @test sort!(shadow) == _volume_canonical_tets(reference)
        # Determinism.
        again = mesh_transfinite_volume(
            corners, cells; recombine=(true, true, true, true, true, true),
            volume_tag=21, face_tags=(11, 12, 13, 14, 15, 16))
        @test again.coords == mesh.coords
        @test _vol_same_blocks(again, mesh)
    end

    @testset "prism-pair masks" begin
        masks = (
            (false, true, false, true, true, true),   # v-free
            (true, false, true, false, true, true),   # u-free
            (true, true, true, true, false, false))   # w-free
        free_slots = ((1, 3), (2, 4), (5, 6))
        for (mask, frees) in zip(masks, free_slots)
            mesh = mesh_transfinite_volume(
                corners, cells; recombine=mask, volume_tag=7,
                face_tags=(11, 12, 13, 14, 15, 16))
            @test validate(mesh).ok
            @test mesh isa Tessella.Elements.MixedMesh
            prisms = _vol_block(mesh, 6)
            tris = _vol_block(mesh, 2)
            quads = _vol_block(mesh, 3)
            @test _vol_block(mesh, 5) === nothing
            @test size(prisms.nodes, 2) == 2 * nu * nv * nw
            @test prisms.tags == fill(Int32(7), 2 * nu * nv * nw)
            free_area = sum(face_cells[s] for s in frees)
            fixed_area = sum(face_cells) - free_area
            @test size(tris.nodes, 2) == 2 * free_area
            @test size(quads.nodes, 2) == fixed_area
            for slot in 1:6
                tag = Int32(10 + slot)
                if slot in frees
                    @test count(==(tag), tris.tags) == 2 * face_cells[slot]
                    @test count(==(tag), quads.tags) == 0
                else
                    @test count(==(tag), quads.tags) == face_cells[slot]
                    @test count(==(tag), tris.tags) == 0
                end
            end
            @test mesh.coords == reference.coords
            # Volume conservation: the cell shadow decompositions must sum
            # to the box volume for every valid mask. The w-free pair uses
            # the PRISM_1/2 tiling; the v/u-free pairs use the R variants
            # indexed to the orientation-fixed emitted order.
            templates = mask == (true, true, true, true, false, false) ?
                (_VOL_SHADOW_PRISM_1, _VOL_SHADOW_PRISM_2) :
                (_VOL_SHADOW_PRISM_1R, _VOL_SHADOW_PRISM_2R)
            volume = 0.0
            for cell in axes(prisms.nodes, 2)
                verts = ntuple(r -> prisms.nodes[r, cell], 6)
                template = templates[isodd(cell) ? 1 : 2]
                for tet in template
                    a = Tuple(mesh.coords[:, verts[tet[1]]])
                    b = Tuple(mesh.coords[:, verts[tet[2]]])
                    c = Tuple(mesh.coords[:, verts[tet[3]]])
                    d = Tuple(mesh.coords[:, verts[tet[4]]])
                    volume += tet_volume(a, b, c, d)
                end
            end
            @test volume ≈ 2.0 atol=1e-12
            # Outward boundary orientation on the affine box.
            for cell in axes(tris.nodes, 2)
                a = Tuple(mesh.coords[:, tris.nodes[1, cell]])
                b = Tuple(mesh.coords[:, tris.nodes[2, cell]])
                c = Tuple(mesh.coords[:, tris.nodes[3, cell]])
                @test orient3(a, b, c, center) > 0
            end
            for cell in axes(quads.nodes, 2)
                a = Tuple(mesh.coords[:, quads.nodes[1, cell]])
                b = Tuple(mesh.coords[:, quads.nodes[2, cell]])
                c = Tuple(mesh.coords[:, quads.nodes[3, cell]])
                d = Tuple(mesh.coords[:, quads.nodes[4, cell]])
                @test orient3(a, b, c, center) > 0
                @test orient3(a, c, d, center) > 0
            end
        end
    end

    @testset "ordered tuples match the Gmsh macros" begin
        unit = mesh_transfinite_volume(
            _affine_corners(), (1, 1, 1); recombine=(true, true, true, true,
                                                     false, false))
        prisms = _vol_block(unit, 6)
        @test prisms.nodes[:, 1] == Int32[1, 2, 3, 5, 6, 7]
        @test prisms.nodes[:, 2] == Int32[4, 3, 2, 8, 7, 6]
        # The w-free pair partitions the reference six-tet subdivision.
        shadow = NTuple{4,Int32}[]
        for (cell, template) in ((1, _VOL_SHADOW_PRISM_1),
                                 (2, _VOL_SHADOW_PRISM_2))
            verts = ntuple(r -> prisms.nodes[r, cell], 6)
            for tet in template
                push!(shadow, Tuple(sort!(
                    Int32[verts[tet[1]], verts[tet[2]],
                          verts[tet[3]], verts[tet[4]]])))
            end
        end
        @test sort!(shadow) == _volume_canonical_tets(
            mesh_transfinite_volume(_affine_corners(), (1, 1, 1)))

        vfree = mesh_transfinite_volume(
            _affine_corners(), (1, 1, 1); recombine=(false, true, false,
                                                     true, true, true))
        prisms = _vol_block(vfree, 6)
        # Gmsh emits the partial-mask prisms with the MPrism orientation
        # fixup applied to the literal macro order.
        @test prisms.nodes[:, 1] == Int32[2, 1, 5, 4, 3, 7]
        @test prisms.nodes[:, 2] == Int32[5, 6, 2, 7, 8, 4]

        ufree = mesh_transfinite_volume(
            _affine_corners(), (1, 1, 1); recombine=(true, false, true,
                                                     false, true, true))
        prisms = _vol_block(ufree, 6)
        @test prisms.nodes[:, 1] == Int32[4, 2, 6, 3, 1, 5]
        @test prisms.nodes[:, 2] == Int32[6, 8, 4, 5, 7, 3]
    end

    @testset "mask parsing and rejection parity" begin
        default = mesh_transfinite_volume(corners, cells)
        for disabled in (nothing, false,
                         (false, false, false, false, false, false))
            mesh = mesh_transfinite_volume(corners, cells;
                                           recombine=disabled)
            @test mesh isa Tessella.MeshTypes.Mesh
            @test mesh_crc(mesh) == mesh_crc(default)
            @test mesh.coords == default.coords
            @test mesh.tets == default.tets
            @test mesh.tris == default.tris
        end

        for bad in ((true, true, true, true, true, false),
                    (false, false, false, false, true, true),
                    (true, true, false, false, true, true),
                    (false, true, true, true, true, false),
                    (false, false, true, false, true, true))
            error = try
                mesh_transfinite_volume(corners, cells; recombine=bad)
                nothing
            catch err
                err
            end
            @test error isa ArgumentError
            @test occursin("wrong surface recombination",
                           sprint(showerror, error))
        end
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; recombine=(true, true, true, true, true))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; recombine=fill(true, 7))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; recombine=(true, true, true, true, true, 1))
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; recombine=1)
        @test_throws ArgumentError mesh_transfinite_volume(
            corners, cells; recombine="yes")
    end

    @testset "warped recombination" begin
        wcorners = _affine_corners((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                   (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        wcorners[7] = (1.2, 1.0, 1.0)
        wcells = (2, 2, 2)
        faces = _volume_face_grids(wcorners, wcells)
        for mask in ((true, true, true, true, true, true),
                     (true, true, true, true, false, false))
            mesh = mesh_transfinite_volume(
                wcorners, wcells; faces=faces, recombine=mask)
            @test validate(mesh).ok
            @test mesh isa Tessella.Elements.MixedMesh
            # The warped path shares the lattice connectivity with the
            # affine path — canonical block contents are identical; only
            # coordinates and emission order differ.
            twin = mesh_transfinite_volume(_affine_corners(), wcells;
                                           recombine=mask)
            @test _vol_canonical_blocks(mesh) == _vol_canonical_blocks(twin)
            # Boundary nodes reuse the face grids bitwise.
            for slot in 1:6
                a1, a2 = _VOL_FACE_AXES[slot]
                n1, n2 = wcells[a1], wcells[a2]
                for p2 in 0:n2, p1 in 0:n1
                    id = _slot_node_id(slot, p1, p2, wcells)
                    @test Tuple(mesh.coords[:, id]) ==
                          Tuple(faces[slot][1][:, p1 + 1 + p2 * (n1 + 1)])
                end
            end
        end
    end

    @testset "recombined allocation growth remains linear" begin
        big = _affine_corners()
        mesh_transfinite_volume(big, (8, 8, 4); recombine=true)
        mesh_transfinite_volume(big, (16, 8, 4); recombine=true)
        small = _volume_recombined_allocated(big, (8, 8, 4))
        large = _volume_recombined_allocated(big, (16, 8, 4))
        @test small > 0
        @test large > small
        @test large <= 2.40small + 1_048_576
        @info "recombined transfinite volume allocation ratchet" small_bytes=small large_bytes=large
    end
end
