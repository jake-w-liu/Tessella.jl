module NativeLookupAbortFractionalV1
using Tessella,Test,TOML
const A=Tessella.API
const ROWS=TOML.parsefile(joinpath(@__DIR__,"native_lookup_abort_fractional_fixtures_v1.toml"))["rows"]
plain(value)=value isa Tuple || value isa AbstractVector ? Any[plain(item) for item in value] : value
function run(prefix)
    outcomes=Any[]
    @testset "Native fractional AbortOnError option normalization" begin
        for row in ROWS
            @testset "requested=$(row["requested"]) normalized=$(row["normalized"])" begin
                try
                    A.initialize();A.mesh.generate(0)
                    A.option("Mesh.Renumber",0)
                    A.option("General.AbortOnError",row["requested"])
                    @test A.option("General.AbortOnError")==row["normalized"]
                    actual=Any[]
                    for (index,expected) in enumerate(row["queries"])
                        record=Dict{String,Any}("kind"=>expected["kind"],"repeat"=>expected["repeat"])
                        try
                            value=expected["kind"]=="node" ? A.mesh.get_node(12345678) : A.mesh.get_element(12345678)
                            merge!(record,Dict("threw"=>false,"result"=>plain(value),"error"=>"","last_error"=>A.logger.get_last_error()))
                        catch failure
                            failure isa InterruptException && rethrow()
                            @test failure isa ErrorException
                            merge!(record,Dict("threw"=>true,"error"=>sprint(showerror,failure),"last_error"=>A.logger.get_last_error()))
                        end
                        @test record==expected
                        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==index
                        push!(actual,record)
                    end
                    @test A.option("General.AbortOnError")==row["after_queries"]
                    reset_carrier=A.mesh.generate(0)
                    @test reset_carrier isa Union{Mesh,MixedMesh}
                    @test reset_carrier!==A.LAST_MESH[]
                    @test A.option("General.AbortOnError")==row["after_reset"]
                    @test A.logger.get_last_error()==row["reset_error"]
                    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                    push!(outcomes,Dict("requested"=>row["requested"],"normalized"=>A.option("General.AbortOnError"),"queries"=>actual))
                finally
                    A.finalize()
                end
            end
        end
    end
    open(prefix*".fractional.toml","w") do stream
        TOML.print(stream,Dict("rows"=>outcomes);sorted=true)
    end
end
end
