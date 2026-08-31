using Test
using Tessella

const _ELEMENT_QUALITY=Tessella.MeshElementQuality

function _quality_fixture()
    return Mesh(
        Float64[0 2 0 0;
                0 0 3 0;
                0 0 0 4];
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1))
end

@noinline function _quality_batch_allocation(count::Int)
    coordinates=Float64[0 1 0 0;0 0 1 0;0 0 0 1]
    tetrahedra=repeat(reshape(Int32[1,2,3,4],4,1),1,count)
    mesh=Mesh(coordinates;tets=tetrahedra)
    tags=UInt64.(1:count)
    _ELEMENT_QUALITY.mesh_element_qualities(mesh,tags,"minSICN")
    GC.gc()
    return @allocated _ELEMENT_QUALITY.mesh_element_qualities(
        mesh,tags,"minSICN")
end

@testset "scale-robust dense simplex element qualities" begin
    fixture=_quality_fixture()

    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"volume")==2.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"minSJ")==1.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"minSICN")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"gamma")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"innerRadius")==1.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"outerRadius")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"angleShape")==1.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"minEdge")==2.0
    @test _ELEMENT_QUALITY.mesh_element_quality(fixture,1,"maxEdge")==2.0
    for quality in ("minDetJac","maxDetJac","minSIGE","minIsotropy")
        @test_throws ArgumentError _ELEMENT_QUALITY.mesh_element_quality(
            fixture,1,quality)
    end

    triangle_expected=Dict(
        "minDetJac"=>6.0,
        "maxDetJac"=>6.0,
        "minSJ"=>1.0,
        "minSICN"=>0.7994080650317894,
        "minSIGE"=>0.9186606921433745,
        "gamma"=>0.7735009811261456,
        "innerRadius"=>0.6972243622680055,
        "outerRadius"=>1.802775637731994,
        "minIsotropy"=>0.7994080650317894,
        "angleShape"=>0.0029285737614364003,
        "minEdge"=>2.0,
        "maxEdge"=>sqrt(13.0),
        "volume"=>3.0)
    tetrahedron_expected=Dict(
        "minDetJac"=>24.0,
        "maxDetJac"=>24.0,
        "minSJ"=>1.0,
        "minSICN"=>0.6988644589181429,
        "minSIGE"=>0.81193794990811,
        "gamma"=>0.6424749609913188,
        "innerRadius"=>0.5766389248992606,
        "outerRadius"=>2.692582403567252,
        "minIsotropy"=>0.7229631432300404,
        "angleShape"=>1.0,
        "minEdge"=>2.0,
        "maxEdge"=>5.0,
        "volume"=>4.0)
    for (quality,expected) in triangle_expected
        @test isapprox(
            _ELEMENT_QUALITY.mesh_element_quality(fixture,2,quality),expected;
            atol=16eps(Float64),rtol=16eps(Float64))
    end
    for (quality,expected) in tetrahedron_expected
        @test isapprox(
            _ELEMENT_QUALITY.mesh_element_quality(fixture,3,quality),expected;
            atol=16eps(Float64),rtol=16eps(Float64))
    end

    ordered=_ELEMENT_QUALITY.mesh_element_qualities(
        fixture,(3,2,3),"volume")
    @test ordered==[4.0,3.0,4.0]
    ordered[1]=99.0
    @test _ELEMENT_QUALITY.mesh_element_qualities(
        fixture,(3,2,3),"volume")==[4.0,3.0,4.0]
    @test isempty(_ELEMENT_QUALITY.mesh_element_qualities(
        fixture,Int[],"minIsotropy"))

    inverted=Mesh(
        Float64[0 0 2 0;
                0 3 0 0;
                0 0 0 4];
        tets=reshape(Int32[1,2,3,4],4,1))
    @test _ELEMENT_QUALITY.mesh_element_quality(
        inverted,1,"minDetJac")==-24.0
    @test _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"minSJ")==-1.0
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"minSICN"),
        -tetrahedron_expected["minSICN"];atol=16eps(Float64),rtol=16eps(Float64))
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"minSIGE"),
        -tetrahedron_expected["minSIGE"];atol=16eps(Float64),rtol=16eps(Float64))
    @test _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"minIsotropy")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"volume")==-4.0
    @test _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"innerRadius")<0
    @test _ELEMENT_QUALITY.mesh_element_quality(inverted,1,"gamma")>0

    degenerate=Mesh(
        Float64[0 1 2 0 1 0 1 0;
                0 0 0 0 0 1 0 1;
                0 0 0 0 0 0 0 0];
        segs=reshape(Int32[1,1],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[4,5,6,7],4,1))
    for tag in 1:3,quality in ("minSJ","minSICN","gamma")
        @test _ELEMENT_QUALITY.mesh_element_quality(
            degenerate,tag,quality)==0.0
    end
    @test _ELEMENT_QUALITY.mesh_element_quality(
        degenerate,2,"minSIGE")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(
        degenerate,3,"minSIGE")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(
        degenerate,2,"minIsotropy")==0.0
    @test _ELEMENT_QUALITY.mesh_element_quality(
        degenerate,3,"minIsotropy")==0.0
    @test isinf(_ELEMENT_QUALITY.mesh_element_quality(
        degenerate,2,"outerRadius"))
    @test isinf(_ELEMENT_QUALITY.mesh_element_quality(
        degenerate,3,"outerRadius"))

    for factor in (2.0^-500,2.0^500)
        scaled=Mesh(fixture.coords.*factor;
                    segs=fixture.segs,tris=fixture.tris,tets=fixture.tets)
        for quality in ("minSJ","minSICN","minSIGE","gamma",
                        "minIsotropy","angleShape")
            for tag in 2:3
                @test isapprox(
                    _ELEMENT_QUALITY.mesh_element_quality(
                        scaled,tag,quality),
                    _ELEMENT_QUALITY.mesh_element_quality(
                        fixture,tag,quality);
                    atol=64eps(Float64),rtol=64eps(Float64))
            end
        end
        for tag in 1:3,quality in ("minEdge","maxEdge","innerRadius",
                                   "outerRadius")
            @test isapprox(
                _ELEMENT_QUALITY.mesh_element_quality(scaled,tag,quality),
                factor*_ELEMENT_QUALITY.mesh_element_quality(
                    fixture,tag,quality);
                atol=0.0,rtol=64eps(Float64))
        end
    end

    maximum=floatmax(Float64)
    overflow_frame=Mesh(
        Float64[-maximum maximum 0;
                0 0 maximum;
                0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(
            overflow_frame,1,"minSICN"),sqrt(3.0)/2;
        atol=64eps(Float64),rtol=64eps(Float64))
    @test isfinite(_ELEMENT_QUALITY.mesh_element_quality(
        overflow_frame,1,"gamma"))
    @test isinf(_ELEMENT_QUALITY.mesh_element_quality(
        overflow_frame,1,"minEdge"))
    mixed_overflow_edges=Mesh(
        Float64[-maximum maximum -maximum;
                0 0 maximum;
                0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    @test _ELEMENT_QUALITY.mesh_element_quality(
        mixed_overflow_edges,1,"minEdge")==maximum
    @test isinf(_ELEMENT_QUALITY.mesh_element_quality(
        mixed_overflow_edges,1,"maxEdge"))

    tiny=nextfloat(0.0)
    underflow_frame=Mesh(
        Float64[0 tiny 0 0;
                0 0 tiny 0;
                0 0 0 tiny];
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1))
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(
            underflow_frame,1,"minSICN"),sqrt(3.0)/2;
        atol=64eps(Float64),rtol=64eps(Float64))
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(
            underflow_frame,2,"minSICN"),sqrt(2/3);
        atol=64eps(Float64),rtol=64eps(Float64))

    extreme_aspect_triangle=Mesh(
        Float64[0 1 0;0 0 tiny;0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    extreme_aspect_tet=Mesh(
        Float64[0 1 0 0;0 0 1 0;0 0 0 tiny];
        tets=reshape(Int32[1,2,3,4],4,1))
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(
            extreme_aspect_triangle,1,"minSIGE"),
        0.769800358919501;atol=16eps(Float64),rtol=16eps(Float64))
    @test isapprox(
        _ELEMENT_QUALITY.mesh_element_quality(
            extreme_aspect_tet,1,"minSIGE"),
        0.5690355937288492;atol=16eps(Float64),rtol=16eps(Float64))
    for (mesh,tag) in ((extreme_aspect_triangle,1),(extreme_aspect_tet,1)),
        quality in ("minSICN","minSIGE","gamma","innerRadius",
                    "outerRadius","minIsotropy","angleShape","minEdge",
                    "maxEdge")
        value=_ELEMENT_QUALITY.mesh_element_quality(mesh,tag,quality)
        @test !isnan(value)
        @test setprecision(BigFloat,64) do
            _ELEMENT_QUALITY.mesh_element_quality(mesh,tag,quality)
        end==value
    end

    for call in (
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,1,"volume"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(true,),"volume"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(1.0,),"volume"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(0,),"volume"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(4,),"volume"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(
            fixture,(big(typemax(Int))+1,),"volume"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(1,),:volume),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(1,),"eta"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(fixture,(1,),"volume\0"),
        ()->_ELEMENT_QUALITY.mesh_element_qualities(
            fixture,(2,1),"minIsotropy"),
    )
        @test_throws ArgumentError call()
    end

    allocation_small=_quality_batch_allocation(10_000)
    allocation_large=_quality_batch_allocation(20_000)
    @test allocation_small<=210_000
    @test allocation_large<=2allocation_small+4_096

    @test isempty(Docs.undocumented_names(
        Tessella.MeshElementQuality;private=false))
    @test isempty(Test.detect_ambiguities(
        Tessella.MeshElementQuality;recursive=true))
end
