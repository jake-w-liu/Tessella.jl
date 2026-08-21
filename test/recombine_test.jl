using Test
using Tessella

function _recombine_square(;width=1.0,tags=Int32[7,7],reverse_second=false)
    coords=Float64[0 width width 0;
                   0 0     1     1;
                   0 0     0     0]
    tris=reverse_second ? Int32[1 1;2 4;3 3] : Int32[1 1;2 3;3 4]
    segs=Int32[1 2 3 4;2 3 4 1]
    return Mesh(coords;segs=segs,tris=tris,seg_tag=Int32[2,2,3,3],tri_tag=tags)
end

function _recombine_grid(n::Int)
    n>0 || throw(ArgumentError("grid extent must be positive"))
    side=n+1;coords=Matrix{Float64}(undef,3,side^2)
    node(i,j)=Int32(i+1+j*side)
    @inbounds for j in 0:n,i in 0:n
        id=Int(node(i,j));coords[1,id]=i;coords[2,id]=j;coords[3,id]=0
    end
    tris=Matrix{Int32}(undef,3,2n^2);tags=Vector{Int32}(undef,2n^2)
    cell=0
    @inbounds for j in 0:n-1,i in 0:n-1
        a=node(i,j);b=node(i+1,j);c=node(i+1,j+1);d=node(i,j+1)
        cell+=1
        tris[1,2cell-1]=a;tris[2,2cell-1]=b;tris[3,2cell-1]=c
        tris[1,2cell]=a;tris[2,2cell]=c;tris[3,2cell]=d
        tags[2cell-1]=Int32(1);tags[2cell]=Int32(1)
    end
    return Mesh(coords;tris=tris,tri_tag=tags)
end

@noinline function _recombine_allocated(n::Int)
    mesh=_recombine_grid(n)
    recombine_triangles(mesh)
    GC.gc()
    return @allocated recombine_triangles(mesh)
end

@testset "deterministic triangle-to-quadrangle recombination" begin
    @testset "square oracle and metadata preservation" begin
        mesh=_recombine_square()
        result=recombine_triangles(mesh;
            physical_names=Dict((1,2)=>"horizontal",(1,3)=>"vertical",(2,7)=>"face"))
        @test Tessella.Elements.validate(result).ok
        @test [block.msh for block in result.blocks]==[1,3]
        @test result.blocks[1].nodes==mesh.segs
        @test result.blocks[1].tags==mesh.seg_tag
        @test result.blocks[2].nodes==reshape(Int32[1,2,3,4],4,1)
        @test result.blocks[2].tags==Int32[7]
        @test result.physical_names==
              Dict((1,2)=>"horizontal",(1,3)=>"vertical",(2,7)=>"face")
        @test Tessella.Elements.mixed_crc(result)==
              Tessella.Elements.mixed_crc(recombine_triangles(mesh;
                  physical_names=result.physical_names))

        without_segments=recombine_triangles(mesh;preserve_segments=false)
        @test [block.msh for block in without_segments.blocks]==[3]
        @test without_segments.blocks[1].nodes==reshape(Int32[1,2,3,4],4,1)
    end

    @testset "selection contracts" begin
        different_tags=recombine_triangles(
            _recombine_square(tags=Int32[7,8]);preserve_segments=false)
        @test [block.msh for block in different_tags.blocks]==[2]
        @test different_tags.blocks[1].nodes==_recombine_square().tris
        @test different_tags.blocks[1].tags==Int32[7,8]

        inconsistent=recombine_triangles(
            _recombine_square(reverse_second=true);preserve_segments=false)
        @test [block.msh for block in inconsistent.blocks]==[2]
        @test size(inconsistent.blocks[1].nodes,2)==2

        elongated=_recombine_square(width=10.0)
        accepted=recombine_triangles(elongated;min_quality=0.09,
                                     preserve_segments=false)
        rejected=recombine_triangles(elongated;min_quality=0.11,
                                     preserve_segments=false)
        @test [block.msh for block in accepted.blocks]==[3]
        @test [block.msh for block in rejected.blocks]==[2]

        concave=Mesh(Float64[0 1 .2 0;0 0 .2 1;0 0 0 0];
                     tris=Int32[1 1;2 3;3 4],tri_tag=Int32[4,4])
        concave_result=recombine_triangles(concave;preserve_segments=false)
        @test [block.msh for block in concave_result.blocks]==[2]

        vertical=Mesh(Float64[0 0 0 0;0 1 1 0;0 0 1 1];
                      tris=Int32[1 1;2 3;3 4],tri_tag=Int32[6,6])
        vertical_result=recombine_triangles(vertical;preserve_segments=false)
        @test vertical_result.blocks[1].msh==3
        @test vertical_result.blocks[1].nodes==reshape(Int32[1,2,3,4],4,1)
    end

    @testset "empty, curve-only, and error paths" begin
        empty_mesh=Mesh(zeros(3,0))
        empty_result=recombine_triangles(empty_mesh)
        @test size(empty_result.coords)==(3,0)
        @test isempty(empty_result.blocks)
        @test Tessella.Elements.validate(empty_result).ok

        curve=Mesh(Float64[0 1;0 0;0 0];segs=reshape(Int32[1,2],2,1),
                   seg_tag=Int32[5])
        curve_result=recombine_triangles(curve)
        @test curve_result.blocks[1].msh==1
        @test curve_result.blocks[1].nodes==curve.segs
        @test curve_result.blocks[1].tags==curve.seg_tag

        tetra=Mesh(Float64[0 1 0 0;0 0 1 0;0 0 0 1];
                   tets=reshape(Int32[1,2,3,4],4,1))
        @test_throws ArgumentError recombine_triangles(tetra)
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=-eps())
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=1.1)
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=NaN)
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=true)
        @test_throws ArgumentError recombine_triangles(
            _recombine_square();preserve_segments=1)
        @test_throws ArgumentError recombine_triangles(
            _recombine_square();physical_names=Dict((2,0)=>"bad"))

        duplicate=Mesh(_recombine_square().coords;
                       tris=Int32[1 1;2 2;3 3],tri_tag=Int32[1,1])
        @test !validate(duplicate).ok
        @test_throws ArgumentError recombine_triangles(duplicate)
    end

    @testset "structured grid and allocation growth" begin
        grid=_recombine_grid(12)
        result=recombine_triangles(grid;preserve_segments=false)
        @test Tessella.Elements.validate(result).ok
        @test length(result.blocks)==1
        @test result.blocks[1].msh==3
        @test size(result.blocks[1].nodes)==(4,144)
        @test all(==(Int32(1)),result.blocks[1].tags)
        @test Tessella.Elements.mixed_crc(result).sha==
              "dbb1bf17965d4e011e7f51a452c6a03e4018a628ffc7e7d9b33d9fc6b922439f"

        small=_recombine_allocated(20)
        large=_recombine_allocated(40)
        @test small>0
        @test large>small
        @test large<=5.25small+262_144
    end
end
