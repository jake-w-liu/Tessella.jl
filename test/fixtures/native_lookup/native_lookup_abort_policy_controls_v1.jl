module NativeLookupAbortPolicyV1
using Tessella,Test,TOML
const A=Tessella.API
const NATIVE_LOOKUP_ABORT_ROWS_V1=TOML.parsefile(
    joinpath(@__DIR__,"native_lookup_abort_policy_fixtures_v1.toml"))["rows"]
bits(values)=String[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
plain(value)=value isa Tuple || value isa AbstractVector ? Any[plain(item) for item in value] : value

function snapshot()
    tags,coordinates,parameters=A.mesh.get_nodes()
    types,elements,cells=A.mesh.get_elements()
    records=Any[]
    for (dimension,entity) in A.model.get_entities()
        nodes,xyz,uv=A.mesh.get_nodes(dimension,entity,false,true)
        families,labels,connectivity=A.mesh.get_elements(dimension,entity)
        push!(records,Dict("owner"=>[dimension,entity],"node_tags"=>nodes,
            "coordinate_bits"=>bits(xyz),"parameter_bits"=>bits(uv),"types"=>families,
            "element_tags"=>labels,"element_nodes"=>connectivity))
    end
    return Dict("node_tags"=>tags,"coordinate_bits"=>bits(coordinates),
        "parameter_bits"=>bits(parameters),"types"=>types,"element_tags"=>elements,
        "element_nodes"=>cells,"records"=>records,"node_max"=>A.mesh.get_max_node_tag(),
        "element_max"=>A.mesh.get_max_element_tag())
end
function check_snapshot(actual,expected)
    @test Set(keys(actual))==Set(keys(expected))
    for key in sort!(collect(keys(expected)))
        @test actual[key]==expected[key]
    end
end
function query(kind,tag)
    try
        result=kind=="node" ? A.mesh.get_node(tag) : A.mesh.get_element(tag)
        return Dict("threw"=>false,"result"=>plain(result),"error"=>"",
            "last_error"=>A.logger.get_last_error())
    catch failure
        failure isa InterruptException && rethrow()
        # Native backend diagnostics throw their message; ArgumentError is
        # reserved for invalid Julia frontend arguments.
        @test failure isa ErrorException
        return Dict("threw"=>true,"error"=>sprint(showerror,failure),
            "last_error"=>A.logger.get_last_error())
    end
end
function seed(row)
    A.initialize();A.mesh.generate(0)
    A.option("Mesh.Renumber",row["renumber"])
    A.option("General.AbortOnError",row["abort"])
    A.model.add("native lookup abort policy")
    A.model.add_discrete_entity(2,50)
    A.mesh.add_nodes(2,50,[101,103,107],[0.,0.,0.,1.,0.,0.,0.,1.,0.])
    A.mesh.add_elements_by_type(50,2,[9001],[101,103,107])
    return nothing
end
# generate returns a detached Tessella geometry carrier; native Msg/state
# parity is checked separately through the complete captured snapshots.
function checked_reset()
    carrier=A.mesh.generate(0)
    @test carrier isa Union{Mesh,MixedMesh}
    @test carrier!==A.LAST_MESH[]
    if !isempty(carrier.coords)
        @test !Base.mightalias(carrier.coords,A.LAST_MESH[].coords)
        carrier.coords[1,1]+=17.0
    end
    return carrier
end
function run(prefix)
    outcomes=Any[]
    @testset "Native lookup abort policy, diagnostics and order lifecycle" begin
        for (index,row) in enumerate(NATIVE_LOOKUP_ABORT_ROWS_V1)
            @testset "native row=$index kind=$(row["kind"]) mode=$(row["abort"]) R=$(row["renumber"]) UID=$(row["tag"]) repeat=$(row["repeat"])" begin
                tag=parse(UInt64,row["tag"])
                actual=Dict{String,Any}(key=>row[key] for key in ("repeat","abort","renumber","kind","tag"))
                try
                    seed(row)
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                    actual["before"]=snapshot();check_snapshot(actual["before"],row["before"])
                    refs=(A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.LAST_MESH_HIGH_ORDER[],
                        A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[])
                    actual["first"]=query(row["kind"],tag)
                    @test actual["first"]==row["first"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
                    @test all(a===b for (a,b) in zip(refs,(A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],
                        A.LAST_MESH_HIGH_ORDER[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[])))
                    actual["after_first"]=snapshot();check_snapshot(actual["after_first"],row["after_first"])
                    actual["known"]=query(row["kind"],row["kind"]=="node" ? UInt64(101) : UInt64(9001))
                    @test actual["known"]==row["known"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
                    actual["second"]=query(row["kind"],tag)
                    @test actual["second"]==row["second"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==2
                    actual["after_repeat"]=snapshot();check_snapshot(actual["after_repeat"],row["after_repeat"])
                    @test A.mesh.set_order(2)===nothing
                    actual["after_dirty_order"]=snapshot();check_snapshot(actual["after_dirty_order"],row["after_dirty_order"])
                    actual["dirty_error"]=A.logger.get_last_error();@test actual["dirty_error"]==row["dirty_error"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==2
                    checked_reset()
                    actual["after_reset"]=snapshot();check_snapshot(actual["after_reset"],row["after_reset"])
                    actual["reset_error"]=A.logger.get_last_error();@test actual["reset_error"]==row["reset_error"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                    @test A.mesh.set_order(2)===nothing
                    actual["after_clean_order"]=snapshot();check_snapshot(actual["after_clean_order"],row["after_clean_order"])
                    actual["clean_error"]=A.logger.get_last_error();@test actual["clean_error"]==row["clean_error"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                    actual["lifetime"]=query(row["kind"],tag)
                    @test actual["lifetime"]==row["lifetime"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
                finally
                    A.finalize()
                end
                try
                    A.initialize()
                    actual["persistent_error"]=A.logger.get_last_error();@test actual["persistent_error"]==row["persistent_error"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
                    checked_reset()
                    actual["cleared_error"]=A.logger.get_last_error();@test actual["cleared_error"]==row["cleared_error"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                finally
                    A.finalize()
                end
                @test actual==row
                push!(outcomes,actual)
            end
        end
    end
    open(prefix*".functional.toml","w") do stream
        TOML.print(stream,Dict("rows"=>outcomes);sorted=true)
    end
    return outcomes
end
end
