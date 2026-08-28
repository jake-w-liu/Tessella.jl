using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _IDENTITY_API=Tessella.API

@testset "synchronized API entity identity lifecycle" begin
    _IDENTITY_API.finalize()
    try
        _IDENTITY_API.initialize()
        _IDENTITY_API.model.add_box(0,0,0,1,1,1;tag=1)
        generated=_IDENTITY_API.mesh.generate(3)
        generated_crc=mesh_crc(generated)
        cached=_IDENTITY_API.LAST_MESH[]

        @test _IDENTITY_API.model.get_entity_name(3,1)==""
        @test _IDENTITY_API.model.get_entity_name(3,99)==""
        @test _IDENTITY_API.LAST_MESH[]===cached
        @test _IDENTITY_API.model.set_entity_name(3,99,"missing")===nothing
        @test _IDENTITY_API.LAST_MESH[]===cached
        @test _IDENTITY_API.model.remove_entity_name("missing")===nothing
        @test _IDENTITY_API.LAST_MESH[]===cached
        @test_throws ArgumentError _IDENTITY_API.model.set_entity_name(3,1,7)
        @test _IDENTITY_API.LAST_MESH[]===cached

        @test _IDENTITY_API.model.set_entity_name(3,1,"domain")===nothing
        @test _IDENTITY_API.LAST_MESH[]===nothing
        @test _IDENTITY_API.model.get_entity_name(3,1)=="domain"
        regenerated=_IDENTITY_API.mesh.generate(3)
        @test mesh_crc(regenerated)==generated_crc
        cached=_IDENTITY_API.LAST_MESH[]
        @test _IDENTITY_API.model.set_entity_name(3,1,"domain")===nothing
        @test _IDENTITY_API.LAST_MESH[]===cached

        @test _IDENTITY_API.model.set_tag(3,1,10)===nothing
        @test _IDENTITY_API.LAST_MESH[]===nothing
        @test _IDENTITY_API.model.get_entities(3)==[(3,10)]
        @test _IDENTITY_API.model.get_entity_name(3,1)==""
        @test _IDENTITY_API.model.get_entity_name(3,10)=="domain"
        retagged=_IDENTITY_API.mesh.generate(3)
        @test mesh_crc(retagged)==generated_crc
        cached=_IDENTITY_API.LAST_MESH[]

        stable_model=deepcopy(_IDENTITY_API.CURRENT[])
        @test_throws ArgumentError _IDENTITY_API.model.set_tag(3,10,10)
        @test_throws ArgumentError _IDENTITY_API.model.set_tag(3,99,11)
        @test_throws ArgumentError _IDENTITY_API.model.set_tag(3,10,0)
        @test _IDENTITY_API.LAST_MESH[]===cached
        @test _IDENTITY_API.CURRENT[].volumes==stable_model.volumes
        @test _IDENTITY_API.CURRENT[].box_extents==stable_model.box_extents
        @test _IDENTITY_API.CURRENT[].entity_names==stable_model.entity_names
        @test _IDENTITY_API.CURRENT[].next_tag==stable_model.next_tag

        @test _IDENTITY_API.model.remove_entity_name("domain")===nothing
        @test _IDENTITY_API.LAST_MESH[]===nothing
        @test _IDENTITY_API.model.get_entity_name(3,10)==""
        @test mesh_crc(_IDENTITY_API.mesh.generate(3))==generated_crc
        cached=_IDENTITY_API.LAST_MESH[]
        @test _IDENTITY_API.model.remove_entity_name("")===nothing
        @test _IDENTITY_API.LAST_MESH[]===cached

        for binding in (:get_entity_name,:set_entity_name,
                        :remove_entity_name,:set_tag)
            @test (@doc getfield(_IDENTITY_API.model,binding))!==nothing
        end
        @test isempty(Docs.undocumented_names(Tessella.API;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _IDENTITY_API.finalize()
    end
end
