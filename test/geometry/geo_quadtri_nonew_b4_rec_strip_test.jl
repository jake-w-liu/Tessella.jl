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

function b4_saved_mesh(stage,dimension,entity)
    owner=only(row for row in stage["entities"] if (row["dim"],row["tag"])==(dimension,entity))
    ids=sort!(unique(Int(node) for cell in owner["cells"] for node in cell["nodes"]))
    position=Dict(node=>slot for (slot,node) in enumerate(ids))
    records=Dict(Int(node["tag"])=>node for node in stage["nodes"])
    coords=hcat((Float64.(records[node]["coordinates"]) for node in ids)...)
    families=sort!(unique(Int(cell["type"]) for cell in owner["cells"]))
    blocks=[ElementBlock(msh,hcat(([Int32(position[Int(node)]) for node in cell["nodes"]]
        for cell in owner["cells"] if cell["type"]==msh)...)) for msh in families]
    return MixedMesh(coords,blocks),ids
end
function b4_saved_supports(raw)
    primary=Dict(parse(Int,node)=>Int(mapped) for (node,mapped) in raw["primary_remap"])
    original=Dict(Int(node["tag"])=>node for node in raw["p1"]["nodes"])
    quadratic=Dict(Int(node["tag"])=>node for node in raw["p2"]["nodes"])
    @test Set(keys(primary))==Set(keys(original))
    for (old,new) in primary
        @test original[old]["coordinate_bits"]==quadratic[new]["coordinate_bits"]
        @test original[old]["owner"]==quadratic[new]["owner"]
    end
    owners=Dict{Tuple,Tuple{Int,Int}}()
    function add(nodes,owner)
        support=B4.key(primary[Int(node)] for node in nodes)
        previous=get(owners,support,nothing)
        if previous===nothing || owner[1]<previous[1];owners[support]=owner
        elseif owner[1]==previous[1];@test owner==previous
        end
    end
    for entity in raw["p1"]["entities"]
        dim=Int(entity["dim"]);dim==0 && continue
        owner=(dim,Int(entity["tag"]))
        for cell in entity["cells"]
            nodes=Int.(cell["nodes"]);msh=Int(cell["type"])
            if dim==1;add(nodes,owner);continue;end
            faces=dim==2 ? (Tuple(eachindex(nodes)),) : B4.FACES[msh]
            for pattern in faces
                face=Tuple(nodes[k] for k in pattern)
                for k in eachindex(face);add((face[k],face[mod1(k+1,length(face))]),owner);end
                add(face,owner)
            end
            msh==5 && add(nodes,owner)
        end
    end
    supports=Dict(parse(Int,node)=>Tuple(Int.(support)) for (node,support) in raw["supports"])
    @test Set(keys(supports))==setdiff(Set(keys(quadratic)),Set(values(primary)))
    @test length(Set(values(supports)))==length(supports)
    for (node,support) in supports
        @test Tuple(Int.(quadratic[node]["owner"]))==owners[B4.key(support)]
        expected=[sum(quadratic[source]["coordinates"][k] for source in support)/length(support) for k in 1:3]
        @test maximum(abs.(Float64.(quadratic[node]["coordinates"]).-expected))<=2e-11
    end
    # Full typed P2 face supports follow actual P1 face identity. A common
    # coordinate is not used as evidence that two interpolation nodes agree.
    basis=Dict{Int,Tuple}(node=>(node,) for node in values(primary));merge!(basis,supports)
    faces=Dict{Tuple,Vector{Tuple}}()
    for entity in raw["p2"]["entities"]
        entity["dim"]==3 || continue
        for cell in entity["cells"]
            family=Int(cell["type"])-7;nodes=Int.(cell["nodes"])
            for pattern in B4.FACES[family]
                support=B4.key(nodes[k] for k in pattern)
                actual=B4.key(node for node in nodes if Set(basis[node])⊆Set(support))
                @test length(actual)==(length(support)==3 ? 6 : 9)
                push!(get!(faces,support,Tuple[]),actual)
            end
        end
    end
    @test all(length(rows) in (1,2) for rows in values(faces))
    @test all(rows[1]==rows[2] for rows in values(faces) if length(rows)==2)
    boundary=Dict(support=>only(rows) for (support,rows) in faces if length(rows)==1)
    lower=Dict(B4.key(cell["nodes"][1:(Int(cell["type"])==9 ? 3 : 4)])=>B4.key(cell["nodes"])
        for entity in raw["p2"]["entities"] if entity["dim"]==2 for cell in entity["cells"])
    @test boundary==lower
end


const artifact=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","quadtri_nonew_b4_rec_strip_oracle.toml"))
@testset "Twelve immutable actual primary B4 captures" begin
    @test artifact["fixture_count"]==12 && artifact["gmsh_version"]=="4.15.2"
    @test artifact["actual_terminal_centers"]==84 && artifact["empty_stored_surface_uv_nodes"]==384
    for record in artifact["fixtures"]
        raw=record["payload"];f=B4.saved_fixture(record)
        @test bytes2hex(sha256(raw["input_geo"]))==record["input_literal_sha256"]
        @test raw["status"]=="certified" && !any(occursin("Error",line) for line in raw["gmsh_log"])
        @test raw["warnings_or_errors"]==record["audit"]["recorded_warnings"]
        source,_=b4_saved_mesh(raw["p1"],2,1)
        linear,_=b4_saved_mesh(raw["p1"],3,1)
        cert=B4.certify(linear,f,source;oracle=true)
        @test length(cert.centers)==f.strip_length && cert.ncell==f.strip_length*(f.intervals+6)
        quadratic,ids=b4_saved_mesh(raw["p2"],3,1)
        qcert=B4.certify_quadratic(quadratic,f,source;oracle=true)
        @test qcert.ncell==cert.ncell
        recorded=Dict(Int(row["cell"])=>row for row in raw["p2_whole_maps"] if row["entity"]==[3,1])
        owner=only(row for row in raw["p2"]["entities"] if (row["dim"],row["tag"])==(3,1))
        vertices=Dict(Int(node["tag"])=>Tuple(Float64.(node["coordinates"])) for node in raw["p2"]["nodes"])
        for cell in owner["cells"]
            map=B4.quadratic_map_certificate(Int(cell["type"]),Tuple(vertices[Int(node)] for node in cell["nodes"]))
            expected=recorded[Int(cell["tag"])]
            @test map.minimum==parse(B4.Q,expected["minimum_exact"])
            @test map.coefficients==expected["coefficients"] && collect(map.degrees)==expected["degrees"]
        end
        b4_saved_supports(raw)
    end
end
include("quadtri_nonew_b4_rec_strip_source_frames.jl")

function b4_product(f;entry=false)
    geometry=B4.execute(f.source;dim=entry ? 3 : 0)
    m=geometry.model;out=Int.(geometry.lists["sweep"])
    source=mesh_model_surface(m,f.surface)
    # Request all final surfaces before the standalone volume, then repeat
    # afterward. Every surface must match the emitted typed boundary.
    surfaces=[(2,tag,mesh_model_surface(m,tag)) for tag in (f.surface,out[1],out[3:end]...)]
    volume=entry ? Tessella.geo_entity_mesh(geometry,3,out[2]) : mesh_model_volume(m,out[2])
    cert=B4.certify(volume,f,source)
    @test Tessella.validate(volume).ok
    @test nnodes(volume)==2*(f.strip_length+1)*length(f.levels)+f.strip_length
    @test cert.families==Dict(filter(row->last(row)>0,
        [4=>2f.strip_length,5=>f.strip_length*(f.intervals-1),7=>5f.strip_length]))
    @test cert.total==cert.source.area*abs(B4.Q(last(cert.heights))-B4.Q(first(cert.heights)))
    @test B4.Quad.surface_faces(surfaces,volume)==Set(keys(cert.boundary))
    projected=model_to_mixed(m,volume,3,out[2])
    @test Tessella.validate(projected).ok
    @test B4.Quad.typed_signature(projected)==B4.Quad.typed_signature(volume)
    own=[count(owner->owner[1]==dim,projected.entity_data.node_entities) for dim in 0:3]
    n=f.intervals;k=f.strip_length
    @test own==[8,4n+4k-8,(2k-2)*(n-1),k]
    for (dim,tag,part) in surfaces
        @test B4.Quad.typed_signature(mesh_model_surface(m,tag))==B4.Quad.typed_signature(part)
    end
    return (;geometry,source,volume,cert,projected)
end

function b4_match_saved_source(record)
    f=B4.saved_fixture(record);product=b4_product(f;entry=true)
    raw=record["payload"];reference,_=b4_saved_mesh(raw["p1"],2,1)
    @test nnodes(reference)==nnodes(product.source)
    positions=Dict{Int,Int}()
    for node in 1:nnodes(product.source)
        matches=findall(i->maximum(abs.(product.source.coords[:,node].-reference.coords[:,i]))<=2e-11,
            1:nnodes(reference))
        @test length(matches)==1
        length(matches)==1 && (positions[node]=only(matches))
    end
    @test length(Set(values(positions)))==nnodes(reference)
    actual=Set(B4.cycle(Tuple(positions[Int(v)] for v in cell))
        for block in B4.blocks(product.source) for cell in eachcol(block.nodes))
    expected=Set(B4.cycle(Tuple(Int.(cell))) for block in B4.blocks(reference) for cell in eachcol(block.nodes))
    @test actual==expected
    for chain in raw["source_chains"]
        actual=product.geometry.model.curve_params[Int(chain["curve"])]
        stored=Float64.(chain["stored_owned_parameters"])
        @test first(actual)==0. && last(actual)==1. && length(actual)==length(stored)+2
        @test maximum(abs.(sort(actual[2:end-1]).-sort(stored));init=0.)<=2e-11
    end
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

get(ENV,"TESSELLA_B4_SAVED_ONLY","0")=="1" || @testset "Dynamic recombined B4 actual geometry" begin
    for m in (5,7,9),layers in (:one,:three,:graded),direction in (-1,1),plane in (:XY,:YZ,:ZX)
        b4_product(B4.fixture("native";strip_length=m,layers,height=direction,plane))
    end
    @testset "Actual sampled convex sources and signed native frames" begin
        for (m,shape,plane,long_pair,winding) in ((5,:rounded,:XY,(1,3),1),
                (7,:skew,:YZ,(2,4),-1),(9,:rounded,:XZ,(1,3),-1),
                (7,:skew,:ZX,(2,4),1)),direction in (-1,1)
            b4_product(B4.fixture("convex_native";strip_length=m,shape,plane,long_pair,
                winding,height=.75direction,layers=:nonbinary,tags=:sparse,
                pins=(3,4,1,2),curve_reverse=5,law="Progression",coefficient=4.))
        end
    end
    @testset "Exact retained primary source recipes" begin
        for record in artifact["fixtures"];b4_match_saved_source(record);end
    end
    @testset "Original source identities and malformed topology" begin
        f=B4.fixture("frames";strip_length=5,layers=:graded)
        b4_source_frames(f);b4_malformed_sources(f);b4_node_frames(f)
        f=B4.fixture("opposite_width_frames";strip_length=7,layers=:one,long_pair=(2,4))
        b4_source_frames(f);b4_node_frames(f)
    end
    @testset "Standalone and GEO rejection is atomic" begin
        b4_atomic_standalone(B4.fixture("atomic";strip_length=5,layers=:one))
    end
end
end
