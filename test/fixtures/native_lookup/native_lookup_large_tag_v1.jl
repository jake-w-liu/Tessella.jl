module NativeLookupLargeElementTagsProbeV2
using Tessella,Test,SHA,TOML,Serialization
const A=Tessella.API
const M=Tessella.Model
const E=Tessella.Elements
include(joinpath(@__DIR__,"native_lookup_large_tag_fixtures_v1.jl"))
bits(xs)=[string(reinterpret(UInt64,Float64(x));base=16,pad=16) for x in xs]
digest(value)=let stream=IOBuffer();serialize(stream,value);bytes2hex(sha256(take!(stream)));end
function seed(row)
    # Exercise the established classified-cache setter contract with literal
    # actual public UIDs. This fixture does not claim native producer narrowing.
    uid=UInt64(row["actual_tag"])
    base=row["producer"]=="add" && uid<=typemax(Int32) ? Int(uid) : 9001
    A.initialize();A.model.add_discrete_entity(2,50)
    A.mesh.add_nodes(2,50,[1,2,3],[0.,0.,0.,1.,0.,0.,0.,1.,0.])
    A.mesh.add_elements_by_type(50,2,[base],[1,2,3])
    lock(A.STATE_LOCK) do
        actual,_=A._mesh_query_pair_locked("large actual UID query fixture")
        data=deepcopy(actual.entity_data)
        data.external_element_tags[1][1]=uid
        mesh=E.MixedMesh(copy(actual.coords),deepcopy(actual.blocks);
            entity_data=data,elementary_entities=deepcopy(actual.elementary_entities))
        class=A._higher_classification_from_data(mesh,"large actual UID query fixture")
        @test only(class.public_tags.element_tags)==uid
        @test (2,Int32(50)) in class.public_tags.authority
        # Replay the captured producer's actual record label and allocation
        # history as well as its detached cache. Leaving the temporary label
        # here would deliberately retain a separate caller-owned element.
        A._raw_model_mesh_record(A.CURRENT[],2,50).element_tags[1]=uid
        A.ELEMENT_TAG_MAX[]=UInt64(row["element_max"])
        A._replace_mesh_cache_locked!(mesh,class;reconcile_records=false)
        # This is a detached classified-cache query fixture. Model producer
        # authority still has its separately recorded Int32 label limit.
        @test A.CURRENT[].meshing.mesh_authority===nothing
    end
    A.option("General.AbortOnError",0)
    return uid
end
function run(output)
    records=Dict{String,Any}[]
    @testset "Native size_t element catalogue and checked quality tags" begin
        for row in NATIVE_LOOKUP_LARGE_TAG_ROWS
            @testset "$(row["producer"])/$(row["requested_tag"])/$(row["repeat"])" begin
                try
                    uid=seed(row)
                    model=A.CURRENT[];meshing=model.meshing;authority=meshing.mesh_authority
                    mesh=A.LAST_MESH[];class=A.LAST_MESH_CLASS[];overlay=A.LAST_MESH_HIGH_ORDER[]
                    before=digest((model,A.OPTIONS,A.SESSION_VIEWS,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
                    @test bits(A.mesh.get_element_qualities(UInt64[uid],"volume"))==row["volume_bits"]
                    @test bits(A.mesh.get_element_qualities(UInt64[uid],"minDetJac"))==row["determinant_bits"]
                    @test bits(A.mesh.get_element_qualities(UInt64[uid,uid+1],"volume"))==row["mixed_volume_bits"]
                    requested=UInt64(row["requested_tag"])
                    @test bits(A.mesh.get_element_qualities(UInt64[requested],"volume"))==row["requested_volume_bits"]
                    types,tags,cells=A.mesh.get_elements()
                    @test Int.(types)==row["types"]
                    @test tags==[UInt64.(v) for v in row["elements"]]
                    @test cells==[UInt64.(v) for v in row["connectivity"]]
                    @test A.mesh.get_max_node_tag()==row["node_max"]
                    @test A.mesh.get_max_element_tag()==row["element_max"]
                    @test A.CURRENT[]===model && model.meshing===meshing && meshing.mesh_authority===authority
                    @test A.LAST_MESH[]===mesh && A.LAST_MESH_CLASS[]===class && A.LAST_MESH_HIGH_ORDER[]===overlay
                    @test digest((model,A.OPTIONS,A.SESSION_VIEWS,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
                    push!(records,Dict("producer"=>row["producer"],"requested_tag"=>string(requested),
                        "actual_tag"=>string(uid),"repeat"=>row["repeat"]))
                finally;A.finalize();end
            end
        end
    end
    open(output,"w") do stream;TOML.print(stream,Dict("julia_version"=>string(VERSION),"rows"=>records,
        "native_bits"=>true,"producer_narrowing_replayed"=>false,
        "query_fixture"=>"detached classified cache",
        "model_authority_Int32_capacity_unimplemented"=>true,
        "prior_scalar_visibility_Int_boundary_unchanged"=>true));end
end
end
