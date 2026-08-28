using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, node, ntets, tet_volume, validate

const _BOOLEAN_API=Tessella.API

@testset "synchronized API Boolean lifecycle" begin
    _BOOLEAN_API.finalize()
    try
        _BOOLEAN_API.initialize()
        _BOOLEAN_API.model.add_box(0,0,0,2,1,1;tag=1)
        cached=_BOOLEAN_API.mesh.generate(3)
        cached_crc=mesh_crc(cached)
        cached_owner=_BOOLEAN_API.LAST_MESH[]
        stable_model=deepcopy(_BOOLEAN_API.CURRENT[])
        @test_throws ArgumentError _BOOLEAN_API.model.boolean_difference(1,999;tag=3)
        @test _BOOLEAN_API.LAST_MESH[]===cached_owner
        @test mesh_crc(_BOOLEAN_API.mesh.get())==cached_crc
        @test _BOOLEAN_API.CURRENT[].volumes==stable_model.volumes
        @test _BOOLEAN_API.CURRENT[].box_extents==stable_model.box_extents
        @test _BOOLEAN_API.CURRENT[].next_tag==stable_model.next_tag

        _BOOLEAN_API.model.add_box(0,0,0,1,1,1;tag=2)
        _BOOLEAN_API.model.add_point(0.5,0.5,0.5;tag=90)
        _BOOLEAN_API.model.embed(0,[90],3,1)
        _BOOLEAN_API.model.add_physical_group(3,[1,2];tag=10,name="operands")
        @test _BOOLEAN_API.model.boolean_difference(1,2;tag=3)==3
        model=_BOOLEAN_API.CURRENT[]
        @test sort!(collect(keys(model.volumes)))==[3]
        @test isempty(model.box_extents)
        @test isempty(model.embeds)
        @test isempty(model.physical)
        @test isempty(model.physical_names)
        @test model.booleans[3]==(op=:difference,a=1,b=2)
        @test haskey(model.boolean_operands,3)
        @test _BOOLEAN_API.LAST_MESH[]===nothing

        cut=_BOOLEAN_API.mesh.generate(3)
        @test validate(cut).ok
        cut_volume=sum(tet_volume(node(cut,cut.tets[1,t]),
                                  node(cut,cut.tets[2,t]),
                                  node(cut,cut.tets[3,t]),
                                  node(cut,cut.tets[4,t]))
                       for t in 1:ntets(cut))
        @test cut_volume≈1.0 atol=1e-12
        cut_crc=mesh_crc(cut)

        _BOOLEAN_API.model.add_sphere(30,0,0,1;tag=1)
        @test haskey(_BOOLEAN_API.CURRENT[].spheres,1)
        @test !haskey(_BOOLEAN_API.CURRENT[].box_extents,1)
        retained=mesh_model_volume(_BOOLEAN_API.CURRENT[],3)
        @test mesh_crc(retained)==cut_crc
    finally
        _BOOLEAN_API.finalize()
    end
end
