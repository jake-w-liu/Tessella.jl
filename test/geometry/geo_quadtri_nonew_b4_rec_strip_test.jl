using Test,Tessella,TOML,SHA
using Tessella.Elements:MixedMesh,ElementBlock
using Tessella.MeshTypes:nnodes
using Tessella.Model:mesh_model_surface,mesh_model_volume,model_to_mixed
if !isdefined(@__MODULE__,:QuadTriNoNewB4RecStripCertificates)
    include("quadtri_nonew_b4_rec_strip_certificates.jl")
end
module NoNewB4GeometryTests
using Test,Tessella,TOML,SHA
using Tessella.Elements:MixedMesh,ElementBlock
using Tessella.MeshTypes:nnodes
using Tessella.Model:mesh_model_surface,mesh_model_volume,model_to_mixed
using ..QuadTriNoNewB4RecStripCertificates
const B4=QuadTriNoNewB4RecStripCertificates



include("quadtri_nonew_b4_rec_strip_source_frames.jl")

function b4_product(f;entry=false)
    geometry=B4.execute(f.source;dim=entry ? 3 : 0)
    m=geometry.model;out=Int.(geometry.lists["sweep"])
    source=mesh_model_surface(m,f.surface)
    # Request all final surfaces before the standalone volume, then repeat
    # afterward. Every surface must match the emitted typed boundary.
    surfaces=[(2,tag,mesh_model_surface(m,tag)) for tag in (f.surface,out[1],out[3:end]...)]
    standalone=mesh_model_volume(m,out[2])
    volume=entry ? Tessella.geo_entity_mesh(geometry,3,out[2];geometry=true) : standalone
    # The emitted view binds the same cell/node identities but keeps its own
    # coordinate arithmetic, so it is certified against the emitted Source and
    # emitted boundary parts. The current native product and fresh standalone
    # constructor each retain their independent classified projection audits.
    emitted_source=entry ? Tessella.geo_entity_mesh(geometry,2,f.surface;geometry=true) : source
    emitted_surfaces=entry ? [(2,tag,Tessella.geo_entity_mesh(geometry,2,tag;geometry=true))
        for tag in (f.surface,out[1],out[3:end]...)] : surfaces
    cert=B4.certify(volume,f,emitted_source)
    @test Tessella.validate(volume).ok
    @test nnodes(volume)==2*(f.strip_length+1)*length(f.levels)+f.strip_length
    @test cert.families==Dict(filter(row->last(row)>0,
        [4=>2f.strip_length,5=>f.strip_length*(f.intervals-1),7=>5f.strip_length]))
    @test cert.total==cert.source.area*abs(B4.Q(last(cert.heights))-B4.Q(first(cert.heights)))
    @test B4.Quad.surface_faces(emitted_surfaces,volume)==Set(keys(cert.boundary))
    if entry
        # The cap's ordinal rule refers to the bound Source producer vector,
        # whose row order differs from the generic entity query's private rank.
        native_source=Tessella.Model._model_mesh_product(m.meshing.mesh_authority,2,f.surface;
            caller="actual B4 Source certificate",producer_order=true)
        native_volume=Tessella.geo_entity_mesh(geometry,3,out[2])
        native_surfaces=[(2,tag,Tessella.geo_entity_mesh(geometry,2,tag))
            for tag in (f.surface,out[1],out[3:end]...)]
        native_cert=B4.certify(native_volume,f,native_source)
        @test native_cert.families==cert.families
        @test native_cert.total==native_cert.source.area*abs(B4.Q(last(native_cert.heights))-B4.Q(first(native_cert.heights)))
        @test B4.Quad.surface_faces(native_surfaces,native_volume)==Set(keys(native_cert.boundary))
        native_projected=model_to_mixed(m,native_volume,3,out[2])
        @test B4.Quad.typed_signature(native_projected)==B4.Quad.typed_signature(native_volume)
    end
    projected=model_to_mixed(m,standalone,3,out[2])
    @test Tessella.validate(projected).ok
    @test B4.Quad.typed_signature(projected)==B4.Quad.typed_signature(standalone)
    if entry
        changed=deepcopy(standalone)
        coordinate=findfirst(!iszero,changed.coords)
        changed.coords[coordinate]=nextfloat(changed.coords[coordinate])
        @test_throws r"input coordinates do not match the completed sweep" model_to_mixed(m,changed,3,out[2])
    end
    own=[count(owner->owner[1]==dim,projected.entity_data.node_entities) for dim in 0:3]
    n=f.intervals;k=f.strip_length
    @test own==[8,4n+4k-8,(2k-2)*(n-1),k]
    for (dim,tag,part) in surfaces
        @test B4.Quad.typed_signature(mesh_model_surface(m,tag))==B4.Quad.typed_signature(part)
    end
    return (;geometry,source=emitted_source,volume,cert,projected)
end


function b4_atomic_standalone(f)
    long_tags=Tuple(f.curve_tags[i] for i in f.long_pair)
    malformed=replace(replace(f.source," RecombLaterals"=>""),
        "Transfinite Curve{$(join(long_tags,','))}=$(f.strip_length+1);"=>
        "Transfinite Curve{$(first(long_tags))}=$(f.strip_length+1);Transfinite Curve{$(last(long_tags))}=$(f.strip_length+2);")
    invalid=(malformed,
        replace(f.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.125,0.0,1.0}"),
        replace(f.source,"Layers{1}"=>"Layers{{1},{0.5}}"))
    for text in invalid
        @test text!=f.source
        geometry=B4.execute(text;dim=0);m=geometry.model;t=Int(geometry.lists["sweep"][2])
        before=(repr(m),m.curve_params,deepcopy(m.curve_params))
        @test_throws r"QuadTriNoNewVerts" mesh_model_volume(m,t)
        @test repr(m)==before[1] && m.curve_params===before[2] && m.curve_params==before[3]
        @test_throws r"QuadTriNoNewVerts" B4.execute(text;dim=3)
    end
end

@testset "Dynamic recombined B4 actual geometry" begin
    for (m,layers,direction,plane) in ((5,:one,1,:XY),(5,:three,-1,:YZ),
                                           (7,:graded,1,:ZX),(9,:graded,-1,:XY))
        b4_product(B4.fixture("native";coherent=true,strip_length=m,layers,height=direction,plane))
    end
    @testset "Actual sampled convex sources and signed native frames" begin
        for (m,shape,plane,long_pair,winding) in ((5,:rounded,:XY,(1,3),1),
                (7,:skew,:YZ,(2,4),-1),(9,:rounded,:XZ,(1,3),-1),
                (7,:skew,:ZX,(2,4),1)),direction in (-1,1)
            b4_product(B4.fixture("convex_native";coherent=true,strip_length=m,shape,plane,long_pair,
                winding,height=.75direction,layers=:nonbinary,tags=:sparse,
                pins=(3,4,1,2),curve_reverse=5,law="Progression",coefficient=4.))
        end
    end
    @testset "Original source identities and malformed topology" begin
        f=B4.fixture("frames";coherent=true,strip_length=5,layers=:graded)
        b4_source_frames(f);b4_malformed_sources(f);b4_node_frames(f)
        f=B4.fixture("opposite_width_frames";coherent=true,strip_length=7,layers=:one,long_pair=(2,4))
        b4_source_frames(f);b4_node_frames(f)
    end
    @testset "Standalone and GEO rejection is atomic" begin
        b4_atomic_standalone(B4.fixture("atomic";coherent=true,strip_length=5,layers=:one))
    end
end
end
