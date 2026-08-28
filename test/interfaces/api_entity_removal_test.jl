using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _REMOVAL_API=Tessella.API

function _removal_api_surface!()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0))
        _REMOVAL_API.model.add_point(x,y,0.0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
        _REMOVAL_API.model.add_line(first_point,last_point;tag=tag)
    end
    _REMOVAL_API.model.add_curve_loop([1,2,3,4];tag=1)
    _REMOVAL_API.model.add_plane_surface([1];tag=1)
    return nothing
end

@testset "synchronized API entity removal lifecycle" begin
    _REMOVAL_API.finalize()
    try
        _REMOVAL_API.initialize()
        _REMOVAL_API.model.add_box(0,0,0,1,1,1;tag=1)
        _REMOVAL_API.model.set_entity_name(3,1,"domain")
        generated=_REMOVAL_API.mesh.generate(3)
        generated_crc=mesh_crc(generated)
        cached=_REMOVAL_API.LAST_MESH[]

        @test _REMOVAL_API.model.remove_entities([(3,99)])===nothing
        @test _REMOVAL_API.LAST_MESH[]===cached
        @test _REMOVAL_API.model.remove_entities([])===nothing
        @test _REMOVAL_API.LAST_MESH[]===cached
        stable=deepcopy(_REMOVAL_API.CURRENT[])
        @test_throws ArgumentError _REMOVAL_API.model.remove_entities([(3,0)])
        @test_throws ArgumentError _REMOVAL_API.model.remove_entities([(4,1)])
        @test_throws ArgumentError _REMOVAL_API.model.remove_entities([(3,1)],1)
        @test _REMOVAL_API.LAST_MESH[]===cached
        @test _REMOVAL_API.CURRENT[].volumes==stable.volumes
        @test _REMOVAL_API.CURRENT[].entity_names==stable.entity_names

        @test _REMOVAL_API.model.remove_entities([(3,1)])===nothing
        @test _REMOVAL_API.LAST_MESH[]===nothing
        @test isempty(_REMOVAL_API.model.get_entities())
        @test _REMOVAL_API.model.get_entity_name(3,1)==""
        @test _REMOVAL_API.CURRENT[].next_tag[4]==1
        @test _REMOVAL_API.model.add_box(2,0,0,1,1,1;tag=0)==2
        @test mesh_crc(_REMOVAL_API.mesh.generate(3)).sha==generated_crc.sha

        _REMOVAL_API.finalize()
        _REMOVAL_API.initialize()
        _removal_api_surface!()
        surface_crc=mesh_crc(_REMOVAL_API.mesh.generate(2))
        cached=_REMOVAL_API.LAST_MESH[]
        @test _REMOVAL_API.model.remove_entities([(1,1)])===nothing
        @test _REMOVAL_API.LAST_MESH[]===cached
        @test (1,1) in _REMOVAL_API.model.get_entities()
        @test _REMOVAL_API.model.remove_entities([(2,1)],true)===nothing
        @test _REMOVAL_API.LAST_MESH[]===nothing
        @test isempty(_REMOVAL_API.model.get_entities())
        @test _REMOVAL_API.model.add_point(0,0,0;tag=1)==1
        @test _REMOVAL_API.model.add_point(1,0,0;tag=2)==2
        @test _REMOVAL_API.model.add_point(1,1,0;tag=3)==3
        @test _REMOVAL_API.model.add_point(0,1,0;tag=4)==4
        @test _REMOVAL_API.model.add_line(1,2;tag=1)==1
        @test _REMOVAL_API.model.add_line(2,3;tag=2)==2
        @test _REMOVAL_API.model.add_line(3,4;tag=3)==3
        @test _REMOVAL_API.model.add_line(4,1;tag=4)==4
        @test _REMOVAL_API.model.add_curve_loop([1,2,3,4];tag=1)==1
        @test _REMOVAL_API.model.add_plane_surface([1];tag=1)==1
        @test mesh_crc(_REMOVAL_API.mesh.generate(2))==surface_crc

        @test (@doc _REMOVAL_API.model.remove_entities)!==nothing
        @test isempty(Docs.undocumented_names(Tessella.API;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _REMOVAL_API.finalize()
    end
end
