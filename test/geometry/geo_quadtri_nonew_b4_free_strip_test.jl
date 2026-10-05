using Test,Tessella,TOML,SHA
if !isdefined(@__MODULE__,:QuadTriNoNewB4FreeStripCertificates)
    include("quadtri_nonew_b4_free_strip_certificates.jl")
end
module NoNewB4FreeGeometryTests
using Test,Tessella,TOML,SHA
using Tessella.Elements:MixedMesh,ElementBlock
using Tessella.MeshTypes:nnodes
using Tessella.Model:mesh_model_surface,mesh_model_volume,model_to_mixed
using ..QuadTriNoNewB4FreeStripCertificates
const B4=QuadTriNoNewB4FreeStripCertificates

function saved_mesh(stage,dimension,entity)
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

function saved_supports(raw)
    primary=Dict(parse(Int,node)=>Int(mapped) for (node,mapped) in raw["primary_remap"])
    original=Dict(Int(node["tag"])=>node for node in raw["p1"]["nodes"])
    quadratic=Dict(Int(node["tag"])=>node for node in raw["p2"]["nodes"])
    @test Set(keys(primary))==Set(keys(original))
    for phase in ("p1","p2"),node in raw[phase]["nodes"]
        for (values,bits) in ((node["coordinates"],node["coordinate_bits"]),(node["stored_parameters"],node["stored_parameter_bits"]))
            @test length(values)==length(bits)
            @test all(string(reinterpret(UInt64,Float64(value));base=16,pad=16)==bit for (value,bit) in zip(values,bits))
        end
    end
    for (old,new) in primary
        @test original[old]["coordinate_bits"]==quadratic[new]["coordinate_bits"]
        @test original[old]["owner"]==quadratic[new]["owner"]
    end
    owners=Dict{Tuple,Tuple{Int,Int}}()
    function add(nodes,owner)
        support=B4.key(primary[Int(node)] for node in nodes);previous=get(owners,support,nothing)
        if previous===nothing || owner[1]<previous[1];owners[support]=owner
        elseif owner[1]==previous[1];@test owner==previous
        end
    end
    for entity in raw["p1"]["entities"]
        dim=Int(entity["dim"]);dim==0 && continue;owner=(dim,Int(entity["tag"]))
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

function product(f;entry=false)
    geometry=B4.execute(f.source;dim=entry ? 3 : 0);model=geometry.model;out=Int.(geometry.lists["sweep"])
    source=mesh_model_surface(model,f.surface)
    surfaces=[(2,tag,mesh_model_surface(model,tag)) for tag in (f.surface,out[1],out[3:end]...)]
    volume=entry ? Tessella.geo_entity_mesh(geometry,3,out[2]) : mesh_model_volume(model,out[2])
    cert=B4.certify(volume,f,source)
    @test Tessella.validate(volume).ok
    @test length(cert.domains)==f.strip_length*f.intervals
    @test cert.total==cert.source.area*abs(B4.Q(last(cert.heights))-B4.Q(first(cert.heights)))
    @test B4.Quad.surface_faces(surfaces,volume)==Set(keys(cert.boundary))
    projected=model_to_mixed(model,volume,3,out[2])
    @test Tessella.validate(projected).ok
    @test B4.Quad.typed_signature(projected)==B4.Quad.typed_signature(volume)
    for (_,tag,part) in surfaces
        @test B4.Quad.typed_signature(mesh_model_surface(model,tag))==B4.Quad.typed_signature(part)
    end
    return (;geometry,source,volume,cert,projected)
end

const artifact=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","quadtri_nonew_b4_free_strip_oracle.toml"))
@testset "Twelve original free B4 primary whole geometry captures" begin
    @test artifact["fixture_count"]==12 && artifact["gmsh_version"]=="4.15.2"
    @test artifact["coordinate_abs_tolerance"]==2e-11
    for record in artifact["fixtures"]
        raw=record["payload"];f=B4.saved_fixture(record)
        @test bytes2hex(sha256(raw["input_geo"]))==record["input_literal_sha256"]
        @test raw["status"]=="certified" && isempty(raw["warnings_or_errors"])
        @test !any(occursin("Error",line) for line in raw["gmsh_log"])
        source,_=saved_mesh(raw["p1"],2,1);linear,_=saved_mesh(raw["p1"],3,1)
        cert=B4.certify(linear,f,source;oracle=true)
        # Zero centers is an observation about these twelve captured products.
        @test isempty(cert.centers) && cert.families==Dict(parse(Int,k)=>Int(v) for (k,v) in record["retained_capture_audit"]["families"])
        quadratic,_=saved_mesh(raw["p2"],3,1);qcert=B4.certify_quadratic(quadratic,f,source;oracle=true)
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
        saved_supports(raw)
    end
end

get(ENV,"TESSELLA_B4_FREE_SAVED_ONLY","0")=="1" || @testset "Actual native free B4 geometry and retained center sources" begin
    for record in artifact["fixtures"]
        f=B4.saved_fixture(record);actual=product(f;entry=true);raw=record["payload"]
        reference,_=saved_mesh(raw["p1"],2,1)
        @test nnodes(reference)==nnodes(actual.source)
        remap=Dict{Int,Int}()
        for node in 1:nnodes(actual.source)
            matches=findall(i->maximum(abs.(actual.source.coords[:,node].-reference.coords[:,i]))<=2e-11,1:nnodes(reference))
            @test length(matches)==1
            length(matches)==1 && (remap[node]=only(matches))
        end
        @test length(Set(values(remap)))==nnodes(reference)
        @test Set(B4.cycle(Tuple(remap[Int(node)] for node in cell)) for block in B4.blocks(actual.source) for cell in eachcol(block.nodes))==
            Set(B4.cycle(Tuple(Int.(cell))) for block in B4.blocks(reference) for cell in eachcol(block.nodes))
        for chain in raw["source_chains"]
            sampled=actual.geometry.model.curve_params[Int(chain["curve"])]
            stored=Float64.(chain["stored_owned_parameters"])
            @test length(sampled)==length(stored)+2
            @test maximum(abs.(sort(sampled[2:end-1]).-sort(stored));init=0.)<=2e-11
        end
    end
    for (m,shape,plane,long_pair,winding) in ((5,:rounded,:XY,(1,3),1),(7,:skew,:YZ,(2,4),-1),(9,:rounded,:XZ,(1,3),-1),(7,:skew,:ZX,(2,4),1)),direction in (-1,1)
        product(B4.fixture("native_variants";strip_length=m,shape,plane,long_pair,winding,
            height=.75direction,layers=:nonbinary,tags=:sparse,pins=(3,4,1,2),curve_reverse=5,law="Progression",coefficient=4.))
    end
    for actual in B4.retained_products()
        cert=actual.certificate
        @test length(cert.centers)==actual.expected_centers
        @test any(last(owner)<actual.f.intervals for owner in values(cert.center_parent))==actual.nonterminal
        @test B4.Quad.typed_signature(actual.projected)==B4.Quad.typed_signature(actual.volume)
        @test all(actual.projected.entity_data.node_entities[node]==(3,Int(actual.geometry.lists["sweep"][2])) for node in cert.centers)
        before=(repr(actual.geometry.model),copy(actual.source.coords),copy(only(B4.blocks(actual.source)).nodes))
        broken=deepcopy(actual.volume);broken.coords[actual.f.axes[3],first(cert.centers)]+=0.03125
        @test_throws r"center differs" B4.certify(broken,actual.f,actual.source)
        @test repr(actual.geometry.model)==before[1] && actual.source.coords==before[2] && only(B4.blocks(actual.source)).nodes==before[3]
    end
end
end
