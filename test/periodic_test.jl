using Test
using Tessella
using Tessella.MeshTypes: Mesh, validate, nnodes

@testset "periodic identification" begin
    # Two-node pair related by translation (1,0,0).
    coords=Float64[0 1 0 1; 0 0 1 1; 0 0 0 0]
    segs=Int32[1 2 3; 2 4 4]
    m=Mesh(coords; segs=segs)
    out=periodic_identify(m,(1.0,0.0,0.0),[1,3],[2,4])
    @test validate(out).ok
    @test out.coords[1,2]≈1.0 && out.coords[1,4]≈1.0
    @test_throws ArgumentError periodic_identify(m,(0,0,0),[1],[2])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1],[3])  # 3 is not 1+(1,0,0)
end
