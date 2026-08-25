using Test
using Random
using LinearAlgebra
using Tessella
using Tessella.MeshTypes: Mesh, validate, mesh_crc

function exact_affine_point(affine::AbstractMatrix{Float64},point::NTuple{3,Float64})
    R=Rational{BigInt}
    return ntuple(3) do row
        Float64(R(affine[row,4])+sum(
            R(affine[row,column])*R(point[column]) for column in 1:3))
    end
end

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
    @test_throws ArgumentError periodic_identify(m,(1,0,0),[1],[2];atol=Inf)

    overflow=Mesh(Float64[floatmax(Float64) 0;0 0;0 0])
    @test_throws ArgumentError periodic_identify(overflow,(floatmax(Float64),0,0),[1],[2];
                                                  atol=floatmax(Float64))
    invalid=Mesh(coords;segs=reshape(Int32[1,1],2,1))
    @test !validate(invalid).ok
    @test_throws ArgumentError periodic_identify(invalid,(1,0,0),[1],[2])

    @testset "general affine correspondence" begin
        # The slave chain is a noisy +90° image of the master chain. The fixed
        # checksum covers coordinate snapping, connectivity, and physical tags.
        affine=[0.0 -1.0 0.0 0.0;
                1.0  0.0 0.0 0.0;
                0.0  0.0 1.0 0.0;
                0.0  0.0 0.0 1.0]
        affine_row_major=Tuple(vec(permutedims(affine)))
        affine_coords=Float64[1 2 3 1e-13 1e-13 1e-13;
                              0 0 0 1     2     3;
                              0 0 0 0     0     0]
        affine_segs=Int32[1 2 4 5;2 3 5 6]
        affine_tags=Int32[7,7,9,9]
        affine_input=Mesh(affine_coords;segs=affine_segs,seg_tag=affine_tags)
        affine_snapshot=copy(affine_input.coords)
        matrix_output=periodic_identify_affine(
            affine_input,affine,[1,2,3],[4,5,6];atol=1e-12)
        flat_output=periodic_identify_affine(
            affine_input,affine_row_major,[1,2,3],[4,5,6];atol=1e-12)
        expected=Float64[1 2 3 0 0 0;0 0 0 1 2 3;0 0 0 0 0 0]
        @test validate(matrix_output).ok
        @test matrix_output.coords==expected
        @test flat_output.coords==expected
        @test matrix_output.segs==affine_input.segs
        @test matrix_output.seg_tag==affine_input.seg_tag
        @test affine_input.coords==affine_snapshot
        @test mesh_crc(matrix_output).sha==
              "ed4b81783f68a3bb092ac6fa2156196efe7f91fec6c7ceeef1aee154fb85261a"

        # Exact-dyadic randomized oracle: upper-triangular linear parts are
        # nonsingular by construction and include reflection, shear, scaling,
        # and translation. Rational arithmetic independently derives each snap.
        rng=MersenneTwister(0x504552494f444943)
        nonzero=(-3,-2,-1,1,2,3)
        saw_reflection=false
        for _ in 1:128
            sx=Float64(rand(rng,nonzero));sy=Float64(rand(rng,nonzero))
            sz=Float64(rand(rng,nonzero))
            saw_reflection|=sx*sy*sz<0
            xy=rand(rng,-4:4)/4; xz=rand(rng,-4:4)/4; yz=rand(rng,-4:4)/4
            tx=rand(rng,-8:8)/8;ty=rand(rng,-8:8)/8;tz=rand(rng,-8:8)/8
            transform=[sx xy xz tx;
                       0.0 sy yz ty;
                       0.0 0.0 sz tz;
                       0.0 0.0 0.0 1.0]
            points=ntuple(4) do point_index
                (Float64(8point_index+rand(rng,-3:3))/8,
                 Float64(5point_index+rand(rng,-3:3))/8,
                 Float64(3point_index+rand(rng,-3:3))/8)
            end
            oracle=map(point->exact_affine_point(transform,point),points)
            random_coords=Matrix{Float64}(undef,3,8)
            for i in 1:4
                random_coords[:,i].=points[i]
                random_coords[:,4+i].=oracle[i]
                random_coords[mod1(i,3),4+i]+=iseven(i) ? 2e-13 : -2e-13
            end
            random_input=Mesh(random_coords)
            row_major=vec(permutedims(transform))
            random_output=periodic_identify_affine(
                random_input,row_major,collect(1:4),collect(5:8);atol=1e-11)
            @test validate(random_output).ok
            @test random_output.coords[:,1:4]==random_input.coords[:,1:4]
            @test all(random_output.coords[:,4+i]==collect(oracle[i]) for i in 1:4)
            @test random_input.coords==random_coords
        end
        @test saw_reflection

        # The exact-dyadic fallback resolves a cancellation whose direct
        # Float64 addition overflows; truly unrepresentable outputs are blocked.
        huge=floatmax(Float64)
        cancellation=Mesh(Float64[huge 0;0 0;0 0])
        cancellation_affine=[1.0 0.0 0.0 -huge;
                             0.0 1.0 0.0 0.0;
                             0.0 0.0 1.0 0.0;
                             0.0 0.0 0.0 1.0]
        cancelled=periodic_identify_affine(
            cancellation,cancellation_affine,[1],[2];atol=0)
        @test cancelled.coords[:,2]==[0.0,0.0,0.0]
        unrepresentable=[2.0 0.0 0.0 0.0;
                         0.0 1.0 0.0 0.0;
                         0.0 0.0 1.0 0.0;
                         0.0 0.0 0.0 1.0]
        @test_throws ArgumentError periodic_identify_affine(
            cancellation,unrepresentable,[1],[2];atol=floatmax(Float64))

        identity=Matrix{Float64}(I,4,4)
        coincident=Mesh(Float64[0 0;0 0;0 0])
        @test periodic_identify_affine(coincident,identity,[1],[2];atol=0).coords==
              coincident.coords
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,affine,[1],[4];atol=1e-16)
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,zeros(3,3),[1],[4])
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,zeros(15),[1],[4])
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,1,[1],[4])
        boolean_affine=Matrix{Bool}(I,4,4)
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,boolean_affine,[1],[4])
        nonfinite=copy(identity);nonfinite[1,1]=Inf
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,nonfinite,[1],[4])
        projective=copy(identity);projective[4,1]=1e-20
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,projective,[1],[4])
        singular=copy(identity);singular[3,3]=0
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,singular,[1],[4])
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,identity,[1],[4];atol=NaN)
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,identity,Int[],Int[])
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,identity,[1,1],[4,5])
        @test_throws ArgumentError periodic_identify_affine(
            affine_input,identity,[1,4],[4,5])
        @test_throws ArgumentError periodic_identify_affine(
            invalid,identity,[1],[2])

        # A tolerance cannot authorize a degenerate output: post-snap mesh
        # validation remains mandatory.
        collapsing=Mesh(Float64[0 1e-3 0;0 0 0;0 0 0];
                        segs=reshape(Int32[2,3],2,1))
        @test validate(collapsing).ok
        @test_throws ErrorException periodic_identify_affine(
            collapsing,identity,[1],[2];atol=0.01)

        @test isempty(Docs.undocumented_names(Tessella.Periodic;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.Periodic;recursive=true))
    end
end
