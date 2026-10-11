module GeoMeshGeometryViewTests
using Tessella,Test
using Tessella.Elements: MixedMesh,msh_spec
include(joinpath(@__DIR__,"quadtri_nonew_certificates.jl"))
include(joinpath(@__DIR__,"quadtri_nonew_quad_patch_certificates.jl"))
const C=QuadTriNoNewCertificates
const P=QuadTriNoNewQuadPatchCertificates
const API=Tessella.API

@testset "Native geometry views and endpoint identities" begin
@testset "Native Mixed geometry getter keeps established P2 CRC and complete queries" begin
    pins=C.artifact_pins()
    for layers in (1,3),laterals in (false,true)
        fixture=C.fixture("geometry_view","Layers{$layers}",collect(range(0.,1.;length=layers+1)),laterals)
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"geometry.geo");write(path,fixture.source)
                API.open_geo!(path;mesh_dim=0)
                first=API.mesh.generate(3);signature=C.typed_signature(first)
                @test first isa MixedMesh
                API.mesh.set_order(2)
                cache=API.LAST_MESH[];class=API.LAST_MESH_CLASS[];model=API.CURRENT[]
                quadratic=API.mesh.get()
                name="uniform$(layers)_$(laterals)_api_p2"
                @test C.crc_record(name,quadratic)==pins[name]
                @test Tessella.validate(quadratic).ok
                @test C.certify_quadratic(quadratic,first).ncell==C.certify_fixture(first,fixture).ncell
                @test all(block->msh_spec(block.msh).order==2,quadratic.blocks)
                @test all(block->Tessella.Elements.msh_dimension(block.msh)==3,quadratic.blocks)
                @test !isempty(API.mesh.get_elements(0)[1])
                @test !isempty(API.mesh.get_elements(1)[1])
                @test !isempty(API.mesh.get_elements(2)[1])
                @test API.LAST_MESH[]===cache && API.LAST_MESH_CLASS[]===class && API.CURRENT[]===model
                API.mesh.set_order(1)
                @test C.typed_signature(API.mesh.get())==signature
                @test C.certify_fixture(API.mesh.get(),fixture).family_counts==C.counts(first)
            end
        finally
            API.finalize()
        end
    end
end

@testset "Native endpoint identities survive exact large offsets" begin
    for plane in (:XY,:YZ,:ZX),direction in (-1,1),laterals in (false,true)
        axes=P.FRAMES[plane]
        offset=ntuple(axis->axis==axes[3] ? 2.0^50+2. : 0.,3)
        fixture=P.fixture("endpoint_identity";plane,height=direction,laterals,layers=:four,offset)
        execution=P.execute(fixture.source;dim=3)
        source=Tessella.Model._model_mesh_product(execution.model.meshing.mesh_authority,
            2,fixture.surface;producer_order=true)
        volume=Tessella.GeoExec.geo_entity_mesh(execution,3,Int(execution.lists["sweep"][2]))
        certificate=P.certify(volume,fixture,source)
        @test abs.(diff(certificate.heights))==fill(.25,4)
        @test certificate.total==1
    end
end
end
end
