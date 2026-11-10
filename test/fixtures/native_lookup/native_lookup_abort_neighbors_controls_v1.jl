module NativeLookupAbortNeighborsV1
using Tessella,Test,TOML
const A=Tessella.API
const ROWS=TOML.parsefile(joinpath(@__DIR__,"native_lookup_abort_neighbors_fixtures_v1.toml"))["rows"]
bits(values)=String[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
plain(value)=value isa Tuple || value isa AbstractVector ? Any[plain(item) for item in value] : value
function snapshot()
    tags,coordinates,parameters=A.mesh.get_nodes()
    types,elements,cells=A.mesh.get_elements()
    return Dict("node_tags"=>tags,"coordinate_bits"=>bits(coordinates),
        "parameter_bits"=>bits(parameters),"types"=>types,"element_tags"=>elements,
        "element_nodes"=>cells,"node_max"=>A.mesh.get_max_node_tag(),"element_max"=>A.mesh.get_max_element_tag())
end
function operation(name)
    tag=UInt64(12345678)
    name=="coordinates" && return A.mesh.get_local_coordinates_in_element(tag,0.,0.,0.)
    name=="jacobian" && return A.mesh.get_jacobian(tag,[0.,0.,0.])
    name=="keys" && return A.mesh.get_keys_for_element(tag,"Lagrange",true)
    name=="quality" && return A.mesh.get_element_qualities([tag],"minSICN")
    name=="set_node" && return A.mesh.set_node(tag,[0.,0.,0.],Float64[])
    error(name)
end
function query(name)
    try
        result=operation(name)
        record=Dict{String,Any}("threw"=>false,"returned_nothing"=>result===nothing,"error"=>"","last_error"=>A.logger.get_last_error())
        result===nothing || (record["result"]=plain(result))
        return record
    catch failure
        failure isa InterruptException && rethrow()
        @test failure isa ErrorException
        return Dict("threw"=>true,"returned_nothing"=>false,"error"=>sprint(showerror,failure),"last_error"=>A.logger.get_last_error())
    end
end
function run(prefix)
    outcomes=Any[]
    @testset "Native neighboring lookup query and setter abort policy" begin
        for (index,row) in enumerate(ROWS)
            @testset "neighbor=$index operation=$(row["operation"]) mode=$(row["abort"]) R=$(row["renumber"]) repeat=$(row["repeat"])" begin
                try
                    A.initialize();A.mesh.generate(0)
                    A.option("General.AbortOnError",row["abort"])
                    A.option("Mesh.Renumber",row["renumber"])
                    A.model.add_discrete_entity(2,50)
                    A.mesh.add_nodes(2,50,[101,103,107],[0.,0.,0.,1.,0.,0.,0.,1.,0.])
                    A.mesh.add_elements_by_type(50,2,[9001],[101,103,107])
                    actual=Dict{String,Any}(key=>row[key] for key in ("repeat","abort","renumber","operation"))
                    actual["before"]=snapshot();@test actual["before"]==row["before"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                    actual["first"]=query(row["operation"]);@test actual["first"]==row["first"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
                    actual["after_first"]=snapshot();@test actual["after_first"]==row["after_first"]
                    actual["second"]=query(row["operation"]);@test actual["second"]==row["second"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==2
                    actual["after_repeat"]=snapshot();@test actual["after_repeat"]==row["after_repeat"]
                    @test actual==row
                    push!(outcomes,actual)
                finally
                    A.finalize()
                end
            end
        end
    end
    open(prefix*".neighbors.toml","w") do stream
        TOML.print(stream,Dict("rows"=>outcomes);sorted=true)
    end
end
end
