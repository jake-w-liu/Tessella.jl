using Test
import Tessella
using Tessella.MeshTypes: Mesh, nnodes, nsegs, ntris, ntets, node,
                          tet_signed_volume, boundary_faces, validate, mesh_crc
using Tessella.Refine: refine_uniform

function _refine_unit_fixture()
    coords = Float64[0 1 0 0;
                     0 0 1 0;
                     0 0 0 1]
    return Mesh(coords;
                segs=reshape(Int32[1, 2], 2, 1),
                tris=reshape(Int32[1, 2, 3], 3, 1),
                tets=reshape(Int32[1, 2, 3, 4], 4, 1),
                seg_tag=Int32[11], tri_tag=Int32[22], tet_tag=Int32[33])
end

function _refine_tet_volumes(mesh::Mesh)
    return [tet_signed_volume(node(mesh, mesh.tets[1, cell]),
                              node(mesh, mesh.tets[2, cell]),
                              node(mesh, mesh.tets[3, cell]),
                              node(mesh, mesh.tets[4, cell]))
            for cell in axes(mesh.tets, 2)]
end

function _refine_face_incidence(tets::Matrix{Int32})
    incidence = Dict{NTuple{3,Int32},Int}()
    for cell in axes(tets, 2)
        a, b, c, d = tets[:, cell]
        for face in ((b, c, d), (a, c, d), (a, b, d), (a, b, c))
            key = Tuple(sort(collect(face)))
            incidence[key] = get(incidence, key, 0) + 1
        end
    end
    return incidence
end

function _refine_segment_chain(count::Int)
    coords = zeros(Float64, 3, count + 1)
    segments = Matrix{Int32}(undef, 2, count)
    for index in 1:count + 1
        coords[1, index] = index - 1
    end
    for cell in 1:count
        segments[1, cell] = Int32(cell)
        segments[2, cell] = Int32(cell + 1)
    end
    return Mesh(coords; segs=segments, seg_tag=fill(Int32(7), count))
end

@noinline function _refine_allocated(count::Int)
    mesh = _refine_segment_chain(count)
    GC.gc()
    return @allocated refine_uniform(mesh)
end

@testset "uniform linear-simplex refinement" begin
    @testset "Gmsh 4.15.2 red templates, tags, and checksum" begin
        source = _refine_unit_fixture()
        refined = refine_uniform(source)

        expected_coords = Float64[0 1 0 0 0.5 0 0 0.5 0.5 0;
                                  0 0 1 0 0 0.5 0 0.5 0 0.5;
                                  0 0 0 1 0 0 0.5 0 0.5 0.5]
        expected_segments = Int32[1 5;
                                  5 2]
        # Gmsh emits the first corner, center, second corner and third corner
        # children in this order for a linear triangle.
        expected_triangles = Int32[1 5 5 6;
                                   5 8 2 8;
                                   6 6 8 3]
        expected_tetrahedra = Int32[1 5 6 7 5 5 6 6;
                                    5 2 8 9 6 9 7 10;
                                    6 8 3 10 7 8 9 9;
                                    7 9 10 4 9 6 10 8]

        @test refined.coords == expected_coords
        @test refined.segs == expected_segments
        @test refined.tris == expected_triangles
        @test refined.tets == expected_tetrahedra
        @test refined.seg_tag == fill(Int32(11), 2)
        @test refined.tri_tag == fill(Int32(22), 4)
        @test refined.tet_tag == fill(Int32(33), 8)
        @test (nnodes(refined), nsegs(refined), ntris(refined), ntets(refined)) ==
              (10, 2, 4, 8)
        volumes = _refine_tet_volumes(refined)
        @test all(volume -> isapprox(volume, 1 / 48; atol=0, rtol=8eps()), volumes)
        @test sum(volumes) ≈ 1 / 6 atol=0 rtol=8eps()
        @test validate(refined).ok
        # Independently derived from the canonical cell sets above using the
        # documented little-endian Mesh-CRC stream.
        @test mesh_crc(refined).sha ==
              "db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031"

        repeated = refine_uniform(source)
        @test repeated.coords == refined.coords
        @test repeated.segs == refined.segs
        @test repeated.tris == refined.tris
        @test repeated.tets == refined.tets
        @test repeated.seg_tag == refined.seg_tag
        @test repeated.tri_tag == refined.tri_tag
        @test repeated.tet_tag == refined.tet_tag
        pristine = _refine_unit_fixture()
        @test source.coords == pristine.coords
        @test source.segs == pristine.segs
        @test source.tris == pristine.tris
        @test source.tets == pristine.tets
        @test source.seg_tag == pristine.seg_tag
        @test source.tri_tag == pristine.tri_tag
        @test source.tet_tag == pristine.tet_tag

        for field in (:coords, :segs, :tris, :tets,
                      :seg_tag, :tri_tag, :tet_tag)
            @test !Base.mightalias(getfield(source, field), getfield(refined, field))
            @test !Base.mightalias(getfield(refined, field), getfield(repeated, field))
        end
    end

    @testset "empty meshes and isolated-node compaction" begin
        isolated = Mesh(Float64[9 -3 4; 8 2 5; 7 1 6])
        empty_refined = refine_uniform(isolated; max_nodes=0, max_cells=0)
        @test size(empty_refined.coords) == (3, 0)
        @test (nsegs(empty_refined), ntris(empty_refined), ntets(empty_refined)) == (0, 0, 0)
        @test validate(empty_refined).ok

        compact_source = Mesh(Float64[0 99 2; 0 99 0; 0 99 0];
                              segs=reshape(Int32[3, 1], 2, 1),
                              seg_tag=Int32[17])
        compact = refine_uniform(compact_source)
        @test compact.coords == Float64[0 2 1; 0 0 0; 0 0 0]
        @test compact.segs == Int32[2 3; 3 1]
        @test compact.seg_tag == Int32[17, 17]
        @test validate(compact).ok
    end

    @testset "shared-face conformity and parent tags" begin
        coords = Float64[0 1 0 0 0;
                         0 0 1 0 0;
                         0 0 0 1 -1]
        source = Mesh(coords;
                      tets=Int32[1 1; 2 3; 3 2; 4 5],
                      tet_tag=Int32[41, 42])
        @test validate(source).ok
        refined = refine_uniform(source)
        @test nnodes(refined) == 14
        @test ntets(refined) == 16
        @test refined.tet_tag == vcat(fill(Int32(41), 8), fill(Int32(42), 8))
        @test validate(refined).ok
        volumes = _refine_tet_volumes(refined)
        @test all(>(0), volumes)
        @test sum(volumes[1:8]) ≈ 1 / 6 atol=0 rtol=16eps()
        @test sum(volumes[9:16]) ≈ 1 / 6 atol=0 rtol=16eps()
        boundary, maximum_incidence = boundary_faces(refined.tets)
        @test length(boundary) == 24
        @test maximum_incidence == 2
        incidence = _refine_face_incidence(refined.tets)
        # Shared coarse face (1,2,3) and its edge midpoints (6,7,10).
        for face in ((Int32(1), Int32(6), Int32(7)),
                     (Int32(6), Int32(2), Int32(10)),
                     (Int32(7), Int32(10), Int32(3)),
                     (Int32(6), Int32(10), Int32(7)))
            @test get(incidence, Tuple(sort(collect(face))), 0) == 2
        end
    end

    @testset "invalid inputs, resource bounds, and Float64 extremes" begin
        source = _refine_unit_fixture()
        inverted = Mesh(source.coords;
                        tets=reshape(Int32[1, 2, 4, 3], 4, 1))
        @test !validate(inverted).ok
        @test_throws ArgumentError refine_uniform(inverted)
        flat = Mesh(Float64[0 1 2; 0 0 0; 0 0 0];
                    tris=reshape(Int32[1, 2, 3], 3, 1))
        @test !validate(flat).ok
        @test_throws ArgumentError refine_uniform(flat)

        @test_throws ArgumentError refine_uniform(source; max_nodes=true)
        @test_throws ArgumentError refine_uniform(source; max_cells=false)
        @test_throws ArgumentError refine_uniform(source; max_nodes=10.0)
        @test_throws ArgumentError refine_uniform(source; max_cells="14")
        @test_throws ArgumentError refine_uniform(source; max_nodes=nothing)
        @test_throws ArgumentError refine_uniform(source; max_cells=missing)
        @test_throws ArgumentError refine_uniform(source; max_nodes=-1)
        @test_throws ArgumentError refine_uniform(source; max_cells=-1)
        @test_throws ArgumentError refine_uniform(
            source; max_nodes=big(typemax(Int32)) + 1)
        @test_throws ArgumentError refine_uniform(
            source; max_cells=big(typemax(Int32)) + 1)
        @test_throws ArgumentError refine_uniform(source; max_nodes=9)
        @test_throws ArgumentError refine_uniform(source; max_cells=13)
        bounded = refine_uniform(source; max_nodes=BigInt(10), max_cells=UInt(14))
        @test (nnodes(bounded), nsegs(bounded) + ntris(bounded) + ntets(bounded)) ==
              (10, 14)

        @test_throws ArgumentError Tessella.Refine._checked_mul(
            typemax(Int), 2, "test")
        @test_throws ArgumentError Tessella.Refine._checked_add(
            typemax(Int), 1, "test")

        # The former subtract-then-add formula double-rounded this midpoint by
        # one ulp. A 256-bit calculation gives the literal expected below.
        midpoint_a = 1.0567732965893495e-70
        midpoint_b = 4.960979170950746e-63
        @test Tessella.Refine._midpoint_coordinate(midpoint_a, midpoint_b) ==
              2.480489638314038e-63
        @test Tessella.Refine._midpoint_coordinate(midpoint_b, midpoint_a) ==
              2.480489638314038e-63

        upper = floatmax(Float64)
        lower = prevfloat(prevfloat(upper))
        extreme = Mesh(Float64[lower upper; 0 0; 0 0];
                       segs=reshape(Int32[1, 2], 2, 1))
        @test validate(extreme).ok
        @test isinf(lower + upper)
        extreme_refined = refine_uniform(extreme)
        @test extreme_refined.coords[1, 3] == prevfloat(upper)
        @test isfinite(extreme_refined.coords[1, 3])
        @test validate(extreme_refined).ok

        tiny = nextfloat(0.0)
        subnormal = Mesh(Float64[tiny 3tiny; 0 0; 0 0];
                         segs=reshape(Int32[1, 2], 2, 1))
        subnormal_refined = refine_uniform(subnormal)
        @test subnormal_refined.coords[1, 3] == 2tiny
        @test validate(subnormal_refined).ok

        unresolved = Mesh(Float64[1 nextfloat(1.0); 0 0; 0 0];
                          segs=reshape(Int32[1, 2], 2, 1))
        @test validate(unresolved).ok
        @test_throws ArgumentError refine_uniform(unresolved)

        coplanar = zeros(Float64, 3, 4)
        coplanar[1, :] = [0, 1, 0, 1]
        coplanar[2, :] = [0, 0, 1, 1]
        output = Matrix{Int32}(undef, 4, 1)
        @test_throws ArgumentError Tessella.Refine._write_positive_tet!(
            output, 1, (Int32(1), Int32(2), Int32(3), Int32(4)),
            coplanar, 1, 1)

        oriented_coords = Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
        Tessella.Refine._write_positive_tet!(
            output, 1, (Int32(1), Int32(3), Int32(2), Int32(4)),
            oriented_coords, 1, 1)
        coordinate_mesh = Mesh(oriented_coords)
        @test tet_signed_volume((node(coordinate_mesh, output[i, 1])
                                 for i in 1:4)...) > 0
    end

    @testset "scale and translation invariance" begin
        reference = refine_uniform(Mesh(
            Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
            tets=reshape(Int32[1, 2, 3, 4], 4, 1), tet_tag=Int32[9]))
        for scale in (1e-300, 1e-200, 1.0, 1e100)
            source = Mesh(
                Float64[0 scale 0 0; 0 0 scale 0; 0 0 0 scale];
                tets=reshape(Int32[1, 2, 3, 4], 4, 1), tet_tag=Int32[9])
            refined = refine_uniform(source)
            @test refined.tets == reference.tets
            @test refined.tet_tag == reference.tet_tag
            @test validate(refined).ok
            @test refined.coords ./ scale ≈ reference.coords atol=0 rtol=2eps()
        end

        offset = 1e100
        width = 16eps(offset)
        translated = Mesh(
            Float64[offset offset + width offset offset;
                    -offset -offset -offset + width -offset;
                    offset offset offset offset + width];
            tets=reshape(Int32[1, 2, 3, 4], 4, 1), tet_tag=Int32[9])
        translated_refined = refine_uniform(translated)
        @test translated_refined.tets == reference.tets
        @test translated_refined.tet_tag == reference.tet_tag
        @test validate(translated_refined).ok
        normalized = copy(translated_refined.coords)
        normalized[1, :] .= (normalized[1, :] .- offset) ./ width
        normalized[2, :] .= (normalized[2, :] .+ offset) ./ width
        normalized[3, :] .= (normalized[3, :] .- offset) ./ width
        @test normalized == reference.coords
    end

    @testset "allocation growth remains linear" begin
        # Warm the exact call path before measuring it.
        refine_uniform(_refine_segment_chain(32))
        small = _refine_allocated(2_000)
        large = _refine_allocated(4_000)
        @test small > 0
        @test large > small
        @test large <= 2.35small + 262_144
    end

    @testset "public documentation" begin
        @test isempty(Base.Docs.undocumented_names(Tessella.Refine; private=false))
        @test isempty(Test.detect_ambiguities(Tessella.Refine; recursive=true))
    end
end
