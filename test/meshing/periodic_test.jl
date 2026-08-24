using Test
using Tessella
using Tessella.MeshTypes: Mesh, validate, mesh_crc

@testset "periodic identification" begin
    # Two-node pair related by translation (1,0,0).
    coords=Float64[0 1 0 1; 0 0 1 1; 0 0 0 0]
    segs=Int32[1 2 3; 2 4 4]
    m=Mesh(coords; segs=segs)
    noisy=copy(coords)
    noisy[1,2]+=5e-13
    noisy[2,4]-=5e-13
    perturbed=Mesh(noisy;segs=segs)
    out=periodic_identify(perturbed,(1.0,0.0,0.0),[1,3],[2,4])
    @test validate(out).ok
    @test out.coords==coords
    @test out.segs==perturbed.segs
    @test perturbed.coords==noisy
    @test mesh_crc(out).sha=="2d4c3e493639ced1a3a2e21a948e73c36b068c18cd37781760faf32eadd8f6f0"
    @test_throws ArgumentError periodic_identify(m,(0,0,0),[1],[2])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1],[3])  # 3 is not 1+(1,0,0)
    @test_throws ArgumentError periodic_identify(m,1,[1],[2])
    @test_throws ArgumentError periodic_identify(m,(true,0,0),[1],[2])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),Bool[true],[2])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[big(2)^100],[2])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1,1],[2,4])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1,3],[2,2])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1,2],[2,4])
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1],[2];atol=true)

    overflow=Mesh(Float64[floatmax(Float64) 0;0 0;0 0])
    @test_throws ArgumentError periodic_identify(overflow,(floatmax(Float64),0,0),[1],[2];
                                                  atol=floatmax(Float64))
    invalid=Mesh(coords;segs=reshape(Int32[1,1],2,1))
    @test !validate(invalid).ok
    @test_throws ArgumentError periodic_identify(invalid,(1,0,0),[1],[2])
end
