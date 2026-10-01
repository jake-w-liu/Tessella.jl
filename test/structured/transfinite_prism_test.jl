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

# Canonical collapsed slot grids for an affine prism — f0=(s0,s1,s5,s4),
# f1=(s1,s2,s6,s5), f2=(s0,s2,s6,s4) are (u,w)/(v,w) bilinear fills; f4/f5
# are the collapsed triangular grids whose i=0 column repeats the apex
# bitwise.
function _affine_prism_faces(corners, nr, ns, nw)
    s0, s1, s2, s4, s5, s6 = corners
    function fill_grid(m, n, evaluate)
        points = Matrix{Float64}(undef, 3, (m + 1) * (n + 1))
        for j in 0:n, i in 0:m
            point = evaluate(i / m, j / n)
            column = i + 1 + j * (m + 1)
            points[1, column] = point[1]
            points[2, column] = point[2]
            points[3, column] = point[3]
        end
        return points
    end
    bilinear(a, b, c, d) = (u, v) -> _prism_lerp3(
        _prism_lerp3(a, b, u), _prism_lerp3(d, c, u), v)
    collapsed(a, b, c) = (u, v) -> _prism_lerp3(
        a, _prism_lerp3(b, c, v), u)
    function quad_tris(m, n)
        result = Matrix{Int32}(undef, 3, 2 * m * n)
        cursor = 0
        for j in 0:n-1, i in 0:m-1
            a = i + 1 + j * (m + 1); b = a + 1
            d = a + m + 1; e = d + 1
            result[:, cursor + 1] .= (a, b, d)
            result[:, cursor + 2] .= (b, e, d)
            cursor += 2
        end
        return result
    end
    function collapsed_tris(m, n)
        result = Matrix{Int32}(undef, 3, n * (2 * m - 1))
        nodeat(i, j) = i + 1 + j * (m + 1)
        cursor = 0
        for j in 0:n-1
            result[:, cursor + 1] .= (nodeat(0, 0),
                                    nodeat(1, j), nodeat(1, j + 1))
            cursor += 1
        end
        for i in 1:m-1, j in 0:n-1
            a = nodeat(i, j); b = nodeat(i + 1, j)
            c = nodeat(i, j + 1); g = nodeat(i + 1, j + 1)
            result[:, cursor + 1] .= (a, b, c)
            result[:, cursor + 2] .= (c, b, g)
            cursor += 2
        end
        return result
    end
    f0 = fill_grid(nr, nw, bilinear(s0, s1, s5, s4))
    f1 = fill_grid(ns, nw, bilinear(s1, s2, s6, s5))
    f2 = fill_grid(nr, nw, bilinear(s0, s2, s6, s4))
    f4 = fill_grid(nr, ns, collapsed(s0, s1, s2))
    f5 = fill_grid(nr, ns, collapsed(s4, s5, s6))
    return ((f0, quad_tris(nr, nw), fill(Int32(11), 2nr * nw)),
            (f1, quad_tris(ns, nw), fill(Int32(12), 2ns * nw)),
            (f2, quad_tris(nr, nw), fill(Int32(13), 2nr * nw)),
            (f4, collapsed_tris(nr, ns), fill(Int32(14), ns * (2nr - 1))),
            (f5, collapsed_tris(nr, ns), fill(Int32(15), ns * (2nr - 1))))
end

# Compact (`transfinite3`) slot grids: the triangular faces fill the full
# (n+1)×(n+1) square with the diagonal-aliased compact lattice.
function _compact_prism_faces(corners, n, nw)
    s0, s1, s2, s4, s5, s6 = corners
    np = n + 1
    function fill_grid(evaluate, m, nn)
        points = Matrix{Float64}(undef, 3, (m + 1) * (nn + 1))
        for j in 0:nn, i in 0:m
            point = evaluate(i, j)
            column = i + 1 + j * (m + 1)
            points[1, column] = point[1]
            points[2, column] = point[2]
            points[3, column] = point[3]
        end
        return points
    end
    bilinear(a, b, c, d) = (i, j) -> _prism_lerp3(
        _prism_lerp3(a, b, i / n), _prism_lerp3(d, c, i / n), j / nw)
    function compact_face(a, b, c)
        fill_grid(n, n) do i, j
            jp = min(j, i)
            t = i == 0 ? 0.0 : jp / i
            _prism_lerp3(a, _prism_lerp3(b, c, t), i / n)
        end
    end
    function compact_tris()
        result = Matrix{Int32}(undef, 3, n * n)
        nodeat(i, j) = i + 1 + j * np
        cursor = 0
        for i in 0:n-1, j in 0:i-1
            result[:, cursor + 1] .= (nodeat(i, j), nodeat(i + 1, j + 1),
                                    nodeat(i, j + 1))
            result[:, cursor + 2] .= (nodeat(i, j), nodeat(i + 1, j),
                                    nodeat(i + 1, j + 1))
            cursor += 2
        end
        for i in 0:n-1
            result[:, cursor + 1] .= (nodeat(i, i), nodeat(i + 1, i),
                                    nodeat(i + 1, i + 1))
            cursor += 1
        end
        return result
    end
    function quad_tris(m, nn)
        result = Matrix{Int32}(undef, 3, 2 * m * nn)
        cursor = 0
        for j in 0:nn-1, i in 0:m-1
            a = i + 1 + j * (m + 1); b = a + 1
            d = a + m + 1; e = d + 1
            result[:, cursor + 1] .= (a, b, d)
            result[:, cursor + 2] .= (b, e, d)
            cursor += 2
        end
        return result
    end
    f0 = fill_grid(bilinear(s0, s1, s5, s4), n, nw)
    f1 = fill_grid(bilinear(s1, s2, s6, s5), n, nw)
    f2 = fill_grid(bilinear(s0, s2, s6, s4), n, nw)
    f4 = compact_face(s0, s1, s2)
    f5 = compact_face(s4, s5, s6)
    return ((f0, quad_tris(n, nw), fill(Int32(11), 2n * nw)),
            (f1, quad_tris(n, nw), fill(Int32(12), 2n * nw)),
            (f2, quad_tris(n, nw), fill(Int32(13), 2n * nw)),
            (f4, compact_tris(), fill(Int32(14), n * n)),
            (f5, compact_tris(), fill(Int32(15), n * n)))
end

# Positional tetrahedron-decomposition templates for the emitted recombined
# cells — the same vertex-order conventions as the kernel's shadow
# decompositions.
const _PRISM_SHADOW_TEMPLATES = (
    ((1, 2, 3, 4), (2, 3, 4, 5), (4, 6, 5, 3)),   # wedge order (a,b,c,d,e,f)
    ((1, 2, 3, 4), (2, 3, 4, 5), (4, 5, 3, 6)),   # PRISM_1/PRISM_4 order
    ((3, 2, 6, 1), (2, 5, 6, 1), (6, 5, 4, 1)),   # PRISM_2 order
    ((2, 1, 5, 3), (1, 4, 5, 3), (5, 4, 6, 3)))   # PRISM_3 order (c,a,g,f,d,h)
const _HEX_SHADOW_TEMPLATE =
    ((1, 2, 4, 5), (2, 4, 5, 6), (5, 6, 4, 8),
     (2, 4, 6, 3), (4, 8, 6, 3), (6, 8, 7, 3))

@inline function _canon4key(a, b, c, d)
    a, b = minmax(a, b); c, d = minmax(c, d)
    a, c = minmax(a, c); b, d = minmax(b, d)
    b, c = minmax(b, c)
    return (a, b, c, d)
end

@inline function _canon3key(a, b, c)
    a <= b || ((a, b) = (b, a))
    b <= c || ((b, c) = (c, b))
    a <= b || ((a, b) = (b, a))
    return (a, b, c)
end

# Decomposes every volume cell of a recombined `MixedMesh` into tetrahedra
# of the reference simplex partition: every shadow template whose tet set is
# contained in `reference_keys` contributes — the union over all cells is
# then compared for exact multiset equality with the partition.
function _recombined_shadow_tets(mesh, reference_keys::Set)
    found = NTuple{4,Int32}[]
    for block in mesh.blocks
        block.msh in (5, 6) || continue
        templates = block.msh == 6 ? _PRISM_SHADOW_TEMPLATES :
                                    (_HEX_SHADOW_TEMPLATE,)
        for cell in axes(block.nodes, 2)
            vertices = block.nodes[:, cell]
            cell_tets = Set{NTuple{4,Int32}}()
            for template in templates
                tets = [_canon4key(vertices[t[1]], vertices[t[2]],
                                   vertices[t[3]], vertices[t[4]])
                        for t in template]
                all(tet -> tet in reference_keys, tets) ||
                    continue
                union!(cell_tets, tets)
            end
            isempty(cell_tets) &&
                error("recombined cell has no simplex-partition decomposition")
            append!(found, cell_tets)
        end
    end
    return found
end

# Audits the emitted boundary cells against the simplex mesh's boundary
# triangulation: every triangle is a boundary face and every quadrangle
# covers exactly two boundary faces — once each, with nothing left over.
function _assert_boundary_cover(mesh, reference)
    remaining = Set(_prism_canonical_triangles(reference.tris))
    for block in mesh.blocks
        block.msh in (2, 3) || continue
        for cell in axes(block.nodes, 2)
            nodes = block.nodes[:, cell]
            if block.msh == 2
                key = _canon3key(nodes[1], nodes[2], nodes[3])
                @test key in remaining
                delete!(remaining, key)
            else
                a, b, c, d = nodes[1], nodes[2], nodes[3], nodes[4]
                covered = 0
                for face in (_canon3key(a, b, c), _canon3key(a, c, d),
                             _canon3key(a, b, d), _canon3key(b, c, d))
                    if face in remaining
                        delete!(remaining, face)
                        covered += 1
                    end
                end
                @test covered == 2
            end
        end
    end
    @test isempty(remaining)
    return nothing
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

    @testset "face-grid boundary input (degenerate-hexahedron slots)" begin
        corners = _affine_prism_corners()
        nr, ns, nw = 2, 3, 2
        affine = mesh_transfinite_prism(corners, (nr, ns, nw))
        meshed = mesh_transfinite_prism(
            corners, (nr, ns, nw); faces=_affine_prism_faces(corners, nr, ns, nw),
            volume_tag=21, face_tags=(11, 12, 13, 14, 15))
        @test validate(meshed).ok
        @test (nnodes(meshed), ntris(meshed), ntets(meshed)) == (27, 46, 54)
        # Boundary tab nodes reuse the face grids bitwise.
        @test Tuple(meshed.coords[:, _prism_node(nr, ns, 0, 0, 1)]) ==
              (0.0, 0.0, 0.5)
        @test Tuple(meshed.coords[:, _prism_node(nr, ns, 2, 0, 0)]) ==
              (3.0, 0.0, 0.0)
        # Topology is canonical — identical connectivity to the affine path.
        @test _prism_canonical_tets(meshed) == _prism_canonical_tets(affine)
        @test _prism_canonical_triangles(meshed.tris) ==
              _prism_canonical_triangles(affine.tris)
        # transfiniteHex on the degenerate slot map reproduces the affine
        # interior nodes to interpolation roundoff.
        @test all(isapprox(meshed.coords[d, node_id],
                           affine.coords[d, node_id];
                           atol=64eps(3.0), rtol=64eps(3.0))
                  for node_id in 1:nnodes(affine), d in 1:3)
        @test count(==(Int32(14)), meshed.tri_tag) == 9
        @test count(==(Int32(15)), meshed.tri_tag) == 9
        @test meshed.tet_tag == fill(Int32(21), ntets(meshed))

        # Warped boundary sheets: interior rows of the f1 grid bow outward —
        # shared edges keep their bitwise values, so the weld stays exact and
        # the interior follows transfiniteHex rather than the affine map.
        warped_faces = _affine_prism_faces(corners, nr, ns, nw)
        warped = copy(warped_faces[2][1])
        for j in 1:ns-1, k in 1:nw-1
            warped[1, j + 1 + k * (ns + 1)] += 0.25
        end
        warped_mesh = mesh_transfinite_prism(
            corners, (nr, ns, nw);
            faces=(warped_faces[1],
                   (warped, warped_faces[2][2], warped_faces[2][3]),
                   warped_faces[3], warped_faces[4], warped_faces[5]))
        @test validate(warped_mesh).ok
        @test _prism_canonical_tets(warped_mesh) == _prism_canonical_tets(affine)
        warped_node = _prism_node(nr, ns, nr, 1, 1)
        @test Tuple(warped_mesh.coords[:, warped_node]) ==
              Tuple(warped[:, 2 + (ns + 1)])
        @test warped_mesh.coords[:, _prism_node(nr, ns, 1, 1, 1)] !=
              affine.coords[:, _prism_node(nr, ns, 1, 1, 1)]

        faces5 = _affine_prism_faces(corners, nr, ns, nw)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw); faces=faces5[1:4])
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw); faces=(faces5..., faces5[1]))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw);
            faces=(faces5[1][1:2], faces5[2], faces5[3], faces5[4], faces5[5]))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw);
            faces=(warped_faces[1][1][:, 1:end-1] |> p -> (p, faces5[1][2],
                   faces5[1][3]), faces5[2], faces5[3], faces5[4], faces5[5]))
        # Non-finite face coordinate.
        bad_finite = copy(faces5[1][1]); bad_finite[3, 2] = NaN
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw);
            faces=((bad_finite, faces5[1][2], faces5[1][3]),
                   faces5[2], faces5[3], faces5[4], faces5[5]))
        # Collapsed apex column must repeat the apex vertex bitwise.
        bad_apex = copy(faces5[4][1]); bad_apex[1, 1 + 1 * (nr + 1)] += 1e-9
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw);
            faces=(faces5[1], faces5[2], faces5[3],
                   (bad_apex, faces5[4][2], faces5[4][3]), faces5[5]))
        # Shared edge disagreement (f0 vs f2 on the s0-s4 collapsed edge —
        # f0's i=0 column at k=1 is a non-corner axial node).
        bad_edge = copy(faces5[1][1]); bad_edge[1, 1 + (nr + 1)] += 1e-9
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw);
            faces=((bad_edge, faces5[1][2], faces5[1][3]),
                   faces5[2], faces5[3], faces5[4], faces5[5]))
        # Corners must coincide with the face grids bitwise.
        corner_mismatch = copy(corners); corner_mismatch[5] =
            corner_mismatch[5] .+ (0.0, 0.0, 1e-9)
        @test_throws ArgumentError mesh_transfinite_prism(
            corner_mismatch, (nr, ns, nw); faces=faces5)
        # Left-handed corner order is still rejected.
        mirrored = [corners[1], corners[3], corners[2],
                    corners[4], corners[6], corners[5]]
        @test_throws ArgumentError mesh_transfinite_prism(
            mirrored, (nr, ns, nw); faces=faces5)
        # Non-affine corners are legal when the face grids carry them — the
        # corner certification only pins the six shared vertices.
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (nr, ns, nw); faces=faces5, face_tags=(1, 2, 3, 4))
    end

    @testset "compact TransfiniteTri=1 subdivision (transfinite3)" begin
        corners = _affine_prism_corners(
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0))
        # Gmsh 4.15.2 emits 6/8/3, 19/32/24, and 46/72/81
        # (nodes/boundary triangles/tetrahedra) for n=1,2,3 — the expanded
        # tab slots behind the diagonal weld onto the diagonal vertices
        # while the interior orphans stay unreferenced, exactly like
        # upstream.
        for (n, nw, counts) in ((1, 1, (6, 8, 3)),
                                (2, 2, (19, 32, 24)),
                                (3, 3, (46, 72, 81)))
            mesh = mesh_transfinite_prism(corners, (n, n, nw); compact=true)
            @test validate(mesh).ok
            @test (nnodes(mesh), ntris(mesh), ntets(mesh)) == counts
            @test all(orient3(node(mesh, mesh.tets[1, tet]),
                              node(mesh, mesh.tets[2, tet]),
                              node(mesh, mesh.tets[3, tet]),
                              node(mesh, mesh.tets[4, tet])) == -1
                      for tet in axes(mesh.tets, 2))
            boundary, maximum_incidence = boundary_faces(mesh.tets)
            @test maximum_incidence == 2
            @test sort!(boundary) == _prism_canonical_triangles(mesh.tris)
        end
        mesh3 = mesh_transfinite_prism(corners, (3, 3, 3); compact=true)
        @test _prism_mesh_volume(mesh3) ≈ 0.5 rtol = 512eps(Float64)
        # The compact lattice fills j-rows (j ≤ i vertices per row) — the
        # lower face lattice occupies nodes 1:(n+1)(n+2)/2 and the upper
        # face lattice ends the mesh, so the six canonical corners land on
        # fixed bitwise positions.
        lattice = (3 + 1) * (3 + 2) ÷ 2
        @test node(mesh3, Int32(1)) == corners[1]          # (0,0) → s0
        @test node(mesh3, Int32(4)) == corners[2]          # (3,0) → s1
        @test node(mesh3, Int32(lattice)) == corners[3]    # (3,3) → s2
        top0 = nnodes(mesh3) - lattice
        @test node(mesh3, Int32(top0 + 1)) == corners[4]   # s4
        @test node(mesh3, Int32(top0 + 4)) == corners[5]   # s5
        @test node(mesh3, Int32(top0 + lattice)) == corners[6]  # s6

        # The compact subdivision requires equal radial/opposite counts —
        # upstream's transfinite3 equal-sides audit.
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 3, 2); compact=true)
        # Bad keyword types are rejected before dispatch.
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (1, 1, 1); compact=1)
    end

    @testset "compact face-grid boundary input (expanded diagonal slots)" begin
        corners = _affine_prism_corners(
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0))
        # Expanded compact grids: slot (i,j) of the full square carries the
        # compact lattice vertex (i,min(j,i)) — upper-triangle slots repeat
        # the diagonal-edge vertex bitwise so the kernel weld sees exactly
        # (n+1)(n+2)/2 distinct nodes per triangular face.
        n, nw = 3, 3
        affine = mesh_transfinite_prism(corners, (n, n, nw); compact=true)
        meshed = mesh_transfinite_prism(
            corners, (n, n, nw);
            faces=_compact_prism_faces(corners, n, nw), compact=true,
            volume_tag=21, face_tags=(11, 12, 13, 14, 15))
        @test validate(meshed).ok
        @test (nnodes(meshed), ntris(meshed), ntets(meshed)) == (46, 72, 81)
        @test _prism_canonical_tets(meshed) == _prism_canonical_tets(affine)
        @test _prism_canonical_triangles(meshed.tris) ==
              _prism_canonical_triangles(affine.tris)
        # Face-boundary nodes reuse the supplied grids bitwise; interior
        # nodes agree to transfiniteHex interpolation roundoff.
        @test all(isapprox(meshed.coords[d, id], affine.coords[d, id];
                           atol=64eps(1.0), rtol=64eps(1.0))
                  for id in 1:nnodes(affine), d in 1:3)
        @test meshed.tet_tag == fill(Int32(21), ntets(meshed))
        @test count(==(Int32(14)), meshed.tri_tag) == n * n
        @test count(==(Int32(15)), meshed.tri_tag) == n * n

        # A face grid whose upper-triangle slots carry distinct positions
        # instead of the aliased diagonal vertex cannot weld to the compact
        # lattice — the distinct-node audit rejects it.
        bad_faces = _compact_prism_faces(corners, n, nw)
        bad4 = copy(bad_faces[4][1])
        bad4[1, 1 + 1 + 2 * (n + 1)] += 0.5   # slot (1,2) — above the diagonal
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (n, n, nw); compact=true,
            faces=(bad_faces[1], bad_faces[2], bad_faces[3],
                   (bad4, bad_faces[4][2], bad_faces[4][3]), bad_faces[5]))
        # Collapsed-style grids (n*(2n-1) triangles) are not compact inputs.
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (n, n, nw); compact=true,
            faces=(bad_faces[1], bad_faces[2], bad_faces[3],
                   (bad_faces[4][1], bad_faces[4][2][:, 1:end-1],
                    bad_faces[4][3][1:end-1]), bad_faces[5]))
    end

    @testset "recombined collapsed emission (wedge prisms + hexes)" begin
        corners = _affine_prism_corners()
        reference = mesh_transfinite_prism(corners, (3, 3, 3))
        reference_tets = _prism_canonical_tets(reference)

        # Gmsh oracle prism_r_all.geo: 18 hexes + 9 prisms, 39 quads + 6 tris.
        mesh = mesh_transfinite_prism(corners, (3, 3, 3); recombine=true)
        @test mesh isa Tessella.Elements.MixedMesh
        blocks = Dict(b.msh => b for b in mesh.blocks)
        @test sort(collect(keys(blocks))) == [2, 3, 5, 6]
        @test size(blocks[6].nodes) == (6, 9)    # wedge prisms
        @test size(blocks[5].nodes) == (8, 18)   # interior hexahedra
        @test size(blocks[2].nodes) == (3, 6)    # wedge tris only
        @test size(blocks[3].nodes) == (4, 39)   # 27 axial + 12 face quads
        @test size(mesh.coords, 2) == 52

        # Cell tetrahedron decomposition tiles the simplex partition exactly.
        @test sort(_recombined_shadow_tets(mesh, Set(reference_tets))) ==
            reference_tets

        # Every boundary cell lies on the shadow boundary, covering each
        # boundary face exactly once.
        _assert_boundary_cover(mesh, reference)

        # Tags propagate.
        tagged = mesh_transfinite_prism(corners, (3, 3, 3); recombine=true,
                                      volume_tag=7,
                                      face_tags=(11, 12, 13, 14, 15))
        tb = Dict(b.msh => b for b in tagged.blocks)
        @test all(==(Int32(7)), tb[5].tags)
        @test all(==(Int32(7)), tb[6].tags)
        @test Set(tb[2].tags) == Set(Int32[14, 15])
        @test Set(tb[3].tags) == Set(Int32[11, 12, 13, 14, 15])
    end

    @testset "recombined collapsed axial-only emission (PRISM_1/PRISM_2)" begin
        corners = _affine_prism_corners()
        # Gmsh oracle prism_r_quads.geo: 45 prisms, 27 quads + 30 tris.
        mesh = mesh_transfinite_prism(corners, (3, 3, 3);
                                      recombine=(true, true, true,
                                                 false, false))
        blocks = Dict(b.msh => b for b in mesh.blocks)
        @test sort(collect(keys(blocks))) == [2, 3, 6]
        @test size(blocks[6].nodes) == (6, 45)
        @test size(blocks[3].nodes) == (4, 27)
        @test size(blocks[2].nodes) == (3, 30)
        @test size(mesh.coords, 2) == 52

        reference = mesh_transfinite_prism(corners, (3, 3, 3))
        @test sort(_recombined_shadow_tets(
            mesh, Set(_prism_canonical_tets(reference)))) ==
            _prism_canonical_tets(reference)
        _assert_boundary_cover(mesh, reference)
    end

    @testset "recombined compact emission (PRISM_3/PRISM_4)" begin
        corners = _affine_prism_corners()
        reference = mesh_transfinite_prism(corners, (3, 3, 3); compact=true)

        # Gmsh oracle prism_rc_all.geo: 27 prisms, 33 quads + 6 tris.
        mesh = mesh_transfinite_prism(corners, (3, 3, 3);
                                      compact=true, recombine=true)
        blocks = Dict(b.msh => b for b in mesh.blocks)
        @test sort(collect(keys(blocks))) == [2, 3, 6]
        @test size(blocks[6].nodes) == (6, 27)
        @test size(blocks[3].nodes) == (4, 33)
        @test size(blocks[2].nodes) == (3, 6)
        @test size(mesh.coords, 2) == 46
        refset = Set(_prism_canonical_tets(reference))
        @test sort(_recombined_shadow_tets(mesh, refset)) ==
            _prism_canonical_tets(reference)
        _assert_boundary_cover(mesh, reference)

        # Axial-only recombination keeps both triangular faces simplex.
        quads_only = mesh_transfinite_prism(corners, (3, 3, 3);
            compact=true, recombine=(true, true, true, false, false))
        qb = Dict(b.msh => b for b in quads_only.blocks)
        @test size(qb[6].nodes) == (6, 27)
        @test size(qb[3].nodes) == (4, 27)
        @test size(qb[2].nodes) == (3, 18)
        _assert_boundary_cover(quads_only, reference)

        # Mixed triangular-face states: f4 only, f5 only.
        for mask in ((true, true, true, true, false),
                     (true, true, true, false, true))
            mixed = mesh_transfinite_prism(corners, (3, 3, 3);
                                           compact=true, recombine=mask)
            mb = Dict(b.msh => b for b in mixed.blocks)
            @test size(mb[6].nodes) == (6, 27)
            @test size(mb[3].nodes) == (4, 30)
            @test size(mb[2].nodes) == (3, 12)
            _assert_boundary_cover(mixed, reference)
        end

        # Arrangement variants emit the same counts with valid coverage.
        ref4 = mesh_transfinite_prism(corners, (4, 4, 2); compact=true)
        for arrangement in (:left, :right, :alternate_left,
                            :alternate_right, (:right, :left))
            am = mesh_transfinite_prism(corners, (4, 4, 2);
                compact=true, recombine=true, arrangement=arrangement)
            ab = Dict(b.msh => b for b in am.blocks)
            @test size(ab[6].nodes) == (6, 32)
            @test size(ab[3].nodes) == (4, 3 * 4 * 2 + 2 * 6)
            @test size(ab[2].nodes) == (3, 8)
            _assert_boundary_cover(am, ref4)
        end
    end

    @testset "recombined masks, dispatch, and error states" begin
        corners = _affine_prism_corners()
        @test mesh_transfinite_prism(corners, (2, 2, 2);
                                     recombine=nothing) isa Tessella.MeshTypes.Mesh
        @test mesh_transfinite_prism(corners, (2, 2, 2);
                                     recombine=false) isa Tessella.MeshTypes.Mesh
        @test mesh_transfinite_prism(corners, (2, 2, 2);
            recombine=(false, false, false, false, false)) isa
            Tessella.MeshTypes.Mesh

        # recombine=nothing and recombine=false are bitwise identical.
        plain = mesh_transfinite_prism(corners, (3, 2, 2))
        explicit = mesh_transfinite_prism(corners, (3, 2, 2); recombine=false)
        @test plain.coords == explicit.coords
        @test plain.tets == explicit.tets
        @test plain.tris == explicit.tris

        for mask in ((true, true, false, true, true),   # partial axial
                     (false, false, false, true, true), # triangular only
                     (true, true, true, true, false),   # mixed tri states
                     (false, true, true, true, true))
            @test_throws ArgumentError mesh_transfinite_prism(
                corners, (2, 2, 1); recombine=mask)
            # In the compact branch triangular-face states are free, but the
            # three axial faces must all be recombined.
            compact_mask = (mask[1], mask[2], mask[3], true, false)
            if mask[1] && mask[2] && mask[3]
                @test mesh_transfinite_prism(corners, (2, 2, 1);
                    compact=true, recombine=compact_mask) isa
                    Tessella.Elements.MixedMesh
            else
                @test_throws ArgumentError mesh_transfinite_prism(
                    corners, (2, 2, 1); compact=true,
                    recombine=compact_mask)
            end
        end
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 2, 1); recombine=(true, true, true))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 2, 1); recombine=(true, true, true, 1, true))
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 2, 1); compact=true, recombine=true,
            arrangement=:bogus)
        @test_throws ArgumentError mesh_transfinite_prism(
            corners, (2, 2, 1); compact=true, recombine=true,
            arrangement=(:left,))
    end

    @testset "recombined warped face-grid emission" begin
        corners = _affine_prism_corners()
        for compact in (false, true)
            simplex = mesh_transfinite_prism(corners, (3, 3, 2);
                                           compact=compact)
            faces = compact ? _compact_prism_faces(corners, 3, 2) :
                              _affine_prism_faces(corners, 3, 3, 2)
            warped = mesh_transfinite_prism(corners, (3, 3, 2);
                compact=compact, faces=faces, recombine=true)
            @test warped isa Tessella.Elements.MixedMesh
            @test validate(warped).ok
            @test sort(_recombined_shadow_tets(
                warped, Set(_prism_canonical_tets(simplex)))) ==
                _prism_canonical_tets(simplex)
            _assert_boundary_cover(warped, simplex)
        end
    end

    @testset "recombined allocation growth remains linear" begin
        small = @allocated mesh_transfinite_prism(
            _affine_prism_corners(), (3, 3, 3); recombine=true)
        large = @allocated mesh_transfinite_prism(
            _affine_prism_corners(), (6, 6, 6); recombine=true)
        tets_small, tets_large = 54, 432
        @test large / max(small, 1) < 4.0 * (tets_large / tets_small)

        csmall = @allocated mesh_transfinite_prism(
            _affine_prism_corners(), (3, 3, 3); compact=true, recombine=true)
        clarge = @allocated mesh_transfinite_prism(
            _affine_prism_corners(), (6, 6, 6); compact=true, recombine=true)
        @test clarge / max(csmall, 1) < 4.0 * (tets_large / tets_small)
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

    @testset "compact allocation growth remains linear in output size" begin
        corners = _affine_prism_corners()
        compact_allocated(cells) = (GC.gc();
            @allocated mesh_transfinite_prism(corners, cells; compact=true))
        mesh_transfinite_prism(corners, (8, 8, 8); compact=true)
        mesh_transfinite_prism(corners, (16, 16, 8); compact=true)
        small = compact_allocated((8, 8, 8))
        large = compact_allocated((16, 16, 8))
        # Doubling the radial cell count quadruples the tetrahedron and
        # interior-node totals; the weld dictionary and output matrices
        # grow linearly with them.
        @test small > 0
        @test large > small
        @test large <= 4.60small + 1_048_576
        @info "compact transfinite prism allocation ratchet" small_bytes=small large_bytes=large
    end
end
