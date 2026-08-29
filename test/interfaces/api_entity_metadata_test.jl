using Test
using Tessella

const _METADATA_API=Tessella.API

@testset "synchronized API entity metadata" begin
    _METADATA_API.finalize()
    try
        _METADATA_API.initialize()
        @test _METADATA_API.model.get_number_of_partitions()==0
        @test_throws ArgumentError _METADATA_API.model.get_entity_type(3,1)

        _METADATA_API.model.add_box(0,0,0,1,2,3;tag=1)
        generated=_METADATA_API.mesh.generate(3)
        cached=_METADATA_API.LAST_MESH[]
        @test _METADATA_API.model.get_entity_type(3,1)=="Volume"
        @test _METADATA_API.model.get_type(3,1)=="Volume"
        @test _METADATA_API.model.get_parent(3,1)==(-1,-1)
        partitions=_METADATA_API.model.get_partitions(3,1)
        @test isempty(partitions)
        push!(partitions,7)
        @test isempty(_METADATA_API.model.get_partitions(3,1))
        @test _METADATA_API.model.get_number_of_partitions()==0
        @test _METADATA_API.LAST_MESH[]===cached

        for call in (
            ()->_METADATA_API.model.get_entity_type(4,1),
            ()->_METADATA_API.model.get_type(3,99),
            ()->_METADATA_API.model.get_parent(3,99),
            ()->_METADATA_API.model.get_partitions(3,99),
        )
            @test_throws ArgumentError call()
            @test _METADATA_API.LAST_MESH[]===cached
        end

        _METADATA_API.model.set_tag(3,1,10)
        @test _METADATA_API.LAST_MESH[]===nothing
        @test _METADATA_API.model.get_entity_type(3,10)=="Volume"
        @test_throws ArgumentError _METADATA_API.model.get_entity_type(3,1)
        @test _METADATA_API.model.remove_entities([(3,10)])===nothing
        @test_throws ArgumentError _METADATA_API.model.get_parent(3,10)

        @test (@doc _METADATA_API.model.get_entity_type)!==nothing
        @test (@doc _METADATA_API.model.get_type)!==nothing
        @test (@doc _METADATA_API.model.get_parent)!==nothing
        @test (@doc _METADATA_API.model.get_number_of_partitions)!==nothing
        @test (@doc _METADATA_API.model.get_partitions)!==nothing
        @test isempty(Docs.undocumented_names(Tessella.API;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _METADATA_API.finalize()
    end
end
