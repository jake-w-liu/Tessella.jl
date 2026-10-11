@testset "predicates, mesh ownership and affine geometry" begin
    predicates = Tessella.Predicates
    a, b, c = (0.0, 0.0), (1.0, 0.0), (0.0, 1.0)
    @test predicates.orient2(a, b, c) == 1
    @test predicates.orient2(b, a, c) == -1
    @test predicates.orient2(a, b, (0.5, 0.0)) == 0
    @test predicates.incircle(a, b, c, (0.25, 0.25)) == 1
    @test predicates.incircle(a, b, c, (1.0, 1.0)) == 0
    @test predicates.orient2((0.0, 0.0), (1e-200, 0.0), (0.0, 1e-200)) == 1
    @test_throws ArgumentError predicates.orient2(a, b, (NaN, 1.0))

    source = unit_tetrahedron()
    @test validate(source).ok
    @test signed_volume(source) ≈ 1 / 6
    faces, incidence = MT.boundary_faces(source.tets)
    @test length(faces) == 4
    @test incidence == 1
    @test MT.euler_characteristic(source) == 1
    @test !validate(Mesh(source.coords; tets=reshape(Int32[1, 3, 2, 4], 4, 1))).ok
    @test_throws ArgumentError Mesh(zeros(2, 4))
    @test_throws ArgumentError Mesh(source.coords; tets=reshape(Int32[1, 2, 3, 5], 4, 1))

    translated = translate_mesh(source, (2.0, -3.0, 4.0))
    @test translated.coords == source.coords .+ [2.0, -3.0, 4.0]
    @test translated.tets == source.tets && translated.tet_tag == source.tet_tag
    @test signed_volume(translated) ≈ 1 / 6
    translated.coords[1, 1] += 10
    @test MT.node(source, 1) == (0.0, 0.0, 0.0)
    scaled = dilate_mesh(source, (0.0, 0.0, 0.0), (2.0, 3.0, 4.0))
    @test validate(scaled).ok && signed_volume(scaled) ≈ 4.0
    reflected = mirror_mesh(source, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    @test validate(reflected).ok && signed_volume(reflected) ≈ 1 / 6
    @test_throws ArgumentError dilate_mesh(source, (0, 0, 0), (0, 1, 1))
end
