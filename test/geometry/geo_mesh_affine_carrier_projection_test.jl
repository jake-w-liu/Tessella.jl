module GeoMeshAffineCarrierProjectionTests
using Tessella,Test
include("quadtri_nonew_certificates.jl")
const C=QuadTriNoNewCertificates
const A=Tessella.API

@testset "Captured affine supports keep native public vertex identities" begin
    fixture=C.fixture("affine_carrier","Layers{1}",[0.,1.],false)
    A.initialize()
    try
        mktemp() do path,io
            write(io,fixture.source);close(io);A.open_geo!(path;mesh_dim=0)
        end
        A.mesh.generate(3)
        callbacks=Ref(0)
        A.mesh.set_size_callback((dim,tag,x,y,z,lc)->(callbacks[]+=1;lc))
        A.mesh.set_order(2)
        cache=A.LAST_MESH[];class=A.LAST_MESH_CLASS[];model=A.CURRENT[]
        authority=model.meshing.mesh_authority
        coordinates=copy(cache.coords)
        parameters=Tessella.Model._model_mesh_runtime_parameters(model,cache)
        view=A.mesh.get()
        @test C.crc_record("uniform1_false_api_p2",view)==C.artifact_pins()["uniform1_false_api_p2"]
        @test (0.5,0.,0.) in Set(Tuple(view.coords[:,row]) for row in axes(view.coords,2))
        @test (0.5,0.,0.) ∉ Set(Tuple(cache.coords[:,row]) for row in axes(cache.coords,2))
        for (row,uid) in enumerate(cache.entity_data.external_node_tags)
            point,parameter,dimension,owner=A.mesh.get_node(uid)
            @test reinterpret(UInt64,point)==reinterpret(UInt64,coordinates[:,row])
            @test (dimension,Int32(owner))==cache.entity_data.node_entities[row]
            original=parameters[row]===nothing ? Float64[] : parameters[row]
            @test reinterpret(UInt64,parameter)==reinterpret(UInt64,original)
        end
        @test callbacks[]==0
        @test cache.coords==coordinates
        @test A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===class && A.CURRENT[]===model
        @test model.meshing.mesh_authority===authority
    finally
        A.finalize()
    end
end
end
