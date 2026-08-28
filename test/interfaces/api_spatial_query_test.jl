using Test
using Tessella

const _SPATIAL_API=Tessella.API

@testset "synchronized API spatial queries" begin
    _SPATIAL_API.finalize()
    try
        _SPATIAL_API.initialize()
        @test_throws ArgumentError _SPATIAL_API.model.get_bounding_box(-1,-1)
        @test isempty(_SPATIAL_API.model.get_entities_in_bounding_box(
            -1,-1,-1,1,1,1))

        _SPATIAL_API.model.add_box(-2,1,3,4,5,6;tag=1)
        generated=_SPATIAL_API.mesh.generate(3)
        cached=_SPATIAL_API.LAST_MESH[]
        @test _SPATIAL_API.model.get_bounding_box(3,1)==
              (-2.0,1.0,3.0,2.0,6.0,9.0)
        @test _SPATIAL_API.model.get_bounding_box(-1,-1)==
              (-2.0,1.0,3.0,2.0,6.0,9.0)
        @test _SPATIAL_API.model.get_entities_in_bounding_box(
            -2,1,3,2,6,9)==[(3,1)]
        @test isempty(_SPATIAL_API.model.get_entities_in_bounding_box(
            -2,1,3,2,6,8.9))
        @test _SPATIAL_API.LAST_MESH[]===cached

        detached=_SPATIAL_API.model.get_entities_in_bounding_box(
            -2,1,3,2,6,9)
        empty!(detached)
        @test _SPATIAL_API.model.get_entities_in_bounding_box(
            -2,1,3,2,6,9)==[(3,1)]
        @test _SPATIAL_API.LAST_MESH[]===cached

        for call in (
            ()->_SPATIAL_API.model.get_bounding_box(3,99),
            ()->_SPATIAL_API.model.get_bounding_box(4,1),
            ()->_SPATIAL_API.model.get_entities_in_bounding_box(
                NaN,1,3,2,6,9),
            ()->_SPATIAL_API.model.get_entities_in_bounding_box(
                -2,1,3,2,6,9,4),
        )
            @test_throws ArgumentError call()
            @test _SPATIAL_API.LAST_MESH[]===cached
        end

        _SPATIAL_API.model.set_tag(3,1,10)
        @test _SPATIAL_API.LAST_MESH[]===nothing
        @test _SPATIAL_API.model.get_bounding_box(3,10)==
              (-2.0,1.0,3.0,2.0,6.0,9.0)
        @test_throws ArgumentError _SPATIAL_API.model.get_bounding_box(3,1)
        @test _SPATIAL_API.model.remove_entities([(3,10)])===nothing
        @test_throws ArgumentError _SPATIAL_API.model.get_bounding_box(-1,-1)

        @test (@doc _SPATIAL_API.model.get_bounding_box)!==nothing
        @test (@doc _SPATIAL_API.model.get_entities_in_bounding_box)!==nothing
        @test isempty(Docs.undocumented_names(Tessella.API;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _SPATIAL_API.finalize()
    end
end
