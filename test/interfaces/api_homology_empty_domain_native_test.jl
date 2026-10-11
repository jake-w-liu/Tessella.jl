module APIHomologyEmptyDomainNativeTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
const NATIVE=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","homology_empty_domain_native_v1.toml"))
const BOUNDARIES=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","homology_boundaries_native_v1.toml"))
const FILE_OVERRIDE=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","homology_file_override_native_v1.toml"))
const QUEUE=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","homology_queue_native_v1.toml"))
const DISCRETE_TYPES=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","discrete_entity_types_native_v1.toml"))
const BASE="Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(1)={1,2};\nTransfinite Curve{1}=2;Physical Curve(10)={1};Physical Point(11)={1};Mesh 1;\n"
include(joinpath(@__DIR__,"native_public_state.jl"))
using .NativePublicState:public_state,assert_state

@testset "native empty homology domain API/GEO diagnostics and state" begin
    @test NATIVE["gmsh_version"]=="4.15.2" && length(NATIVE["cases"])==8
    @test !NATIVE["direct_native_Msg_counter_read"]
    for (name,expected) in sort!(collect(NATIVE["cases"]);by=first)
        @testset "$name" begin
            case=expected["case"];abort=case["abort"];subdomain=case["subdomain"]
            A.initialize()
            try
                @test A.option("General.AbortOnError")==expected["default_abort"]
                mktempdir() do directory
                    base=joinpath(directory,"interval.geo");write(base,BASE)
                    A.open_geo!(base)
                    assert_state(expected["before"])
                    before_last_error=A.logger.get_last_error()
                    before_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                    call=if case["route"]=="api"
                        A.option("General.AbortOnError",abort)
                        A.mesh.add_homology_request("Homology",[99],subdomain,[0,1])
                        ()->A.mesh.generate(1)
                    else
                        # Public GEO parses the identical complete scene in a
                        # new call; the native fragment route is kept in the
                        # fixture and checked separately through its full state.
                        # Inherit default Abort2. The file must override it to
                        # Abort0 before any homology diagnostic is emitted.
                        source="General.AbortOnError=$abort;\n"*BASE*
                            "Homology({0,1}){{99},{"*join(subdomain,',')*"}};\nMesh 1;\n"
                        path=joinpath(directory,"homology.geo");write(path,source)
                        ()->A.open_geo!(path)
                    end
                    if expected["status"]=="success"
                        result=call()
                        @test result!==nothing
                        assert_state(expected["after"])
                        @test A.logger.get_last_error()==expected["last_error"]
                        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==length(expected["error_events"])
                        case["route"]=="geo" && @test result.msg_error_count==length(expected["error_events"])
                    else
                        error=try call();nothing catch e;e isa InterruptException && rethrow();e end
                        @test error isa Exception
                        @test occursin(expected["exception"],sprint(showerror,error))
                        # Native generation has already run before the diagnostic;
                        # surviving cells remain fixed while its clock grows3→5.
                        if case["route"]=="api"
                            @test A.logger.get_last_error()==expected["last_error"]
                            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==length(expected["error_events"])
                            assert_state(expected["after"])
                        else
                            # This public frontend publishes a new scene only
                            # after execution succeeds. Its unchanged previous
                            # scene is an atomic-open contract; native merge of
                            # the fragment instead consumes the clocks above.
                            assert_state(expected["before"])
                            @test A.logger.get_last_error()==before_last_error
                            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==before_errors
                        end
                    end
                end
            finally A.finalize() end
        end
    end
end

@testset "native empty, unmeshed and invalid-dimension homology boundaries" begin
    @test BOUNDARIES["gmsh_version"]=="4.15.2" && length(BOUNDARIES["cases"])==8
    @test !BOUNDARIES["direct_native_Msg_counter_read"]
    for (name,expected) in sort!(collect(BOUNDARIES["cases"]);by=first)
        @testset "$name" begin
            case=expected["case"];scene=case["scene"]
            A.initialize()
            try
                @test A.option("General.AbortOnError")==2
                A.option("General.AbortOnError",case["abort"])
                if scene!="empty"
                    A.model.add_point(0.,0.,0.;tag=1,meshSize=1.)
                    if scene!="point"
                        A.model.add_point(1.,0.,0.;tag=2,meshSize=1.)
                        A.model.add_line(1,2;tag=1)
                        A.mesh.set_transfinite_curve(1,2)
                    end
                    A.model.add_physical_group(scene=="point" ? 0 : 1,[1];tag=11)
                    scene=="meshed_line" && A.mesh.generate(1)
                end
                assert_state(expected["before"])
                # computeHomology does not generate or reset process-global
                # Msg state. Compare the exact event delta, preserving earlier
                # cases' diagnostic history without an ambient manual reset.
                previous_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                previous_last_error=A.logger.get_last_error()
                @test A.mesh.add_homology_request("Homology",case["domain"],Int[],case["dims"])===nothing
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==previous_errors
                @test A.logger.get_last_error()==previous_last_error
                if expected["status"]=="success"
                    result=A.mesh.compute_homology()
                    @test result isa Vector{Tuple{Int,Int}}
                    @test isempty(expected["value"]) && isempty(result)
                else
                    error=try A.mesh.compute_homology();nothing catch e;e isa InterruptException && rethrow();e end
                    @test error isa ErrorException
                    @test occursin(expected["exception"],sprint(showerror,error))
                end
                assert_state(expected["after"])
                @test A.logger.get_last_error()==expected["last_error"]
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]-previous_errors==length(expected["error_events"])
            finally A.finalize() end
        end
    end
end

@testset "GEO AbortOnError override publishes only after successful import" begin
    @test FILE_OVERRIDE["gmsh_version"]=="4.15.2" && length(FILE_OVERRIDE["cases"])==2
    before=NATIVE["cases"]["abort2_api_empty"]["before"]
    for (name,expected) in sort!(collect(FILE_OVERRIDE["cases"]);by=first)
        @testset "$name" begin
            @test bytes2hex(sha256(expected["scene"]))==expected["scene_sha256"]
            A.initialize()
            try
                mktempdir() do directory
                    base=joinpath(directory,"interval.geo");write(base,BASE)
                    A.open_geo!(base);assert_state(before)
                    A.option("General.AbortOnError",expected["inherited"])
                    before_last_error=A.logger.get_last_error()
                    before_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                    source=joinpath(directory,"override.geo");write(source,expected["scene"])
                    if expected["status"]=="success"
                        result=A.open_geo!(source)
                        assert_state(expected["after"])
                        @test result.msg_error_count==length(expected["error_events"])
                    else
                        error=try A.open_geo!(source);nothing catch e;e isa InterruptException && rethrow();e end
                        @test error isa ErrorException
                        @test occursin(expected["exception"],sprint(showerror,error))
                        # Failed imports leave options and diagnostics with
                        # the same publication boundary as model/cache state.
                        assert_state(before)
                    end
                    if expected["status"]=="success"
                        @test A.option("General.AbortOnError")==expected["final_abort"]
                        @test A.logger.get_last_error()==expected["last_error"]
                        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==length(expected["error_events"])
                    else
                        @test A.option("General.AbortOnError")==expected["inherited"]
                        @test A.logger.get_last_error()==before_last_error
                        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==before_errors
                    end
                end
            finally A.finalize() end
        end
    end
end

@testset "native homology queue accepts safe dimensions before compute diagnostics" begin
    @test QUEUE["gmsh_version"]=="4.15.2" && length(QUEUE["cases"])==20
    @test !QUEUE["direct_native_Msg_counter_read"]
    for (name,expected) in sort!(collect(QUEUE["cases"]);by=first)
        @testset "$name" begin
            case=expected["case"]
            @test expected["queue_returns_none"] && isempty(expected["queue_messages"])
            @test expected["after_scalar_cell_expectations_derived_from_owned_bulk"]
            @test expected["native_stale_scalar_oracle_not_qualification"]
            A.initialize()
            try
                A.option("General.AbortOnError",case["abort"])
                A.model.add_point(0.,0.,0.;tag=1,meshSize=1.)
                A.model.add_point(1.,0.,0.;tag=2,meshSize=1.)
                A.model.add_line(1,2;tag=1);A.mesh.set_transfinite_curve(1,2)
                A.model.add_physical_group(1,[1];tag=11);A.mesh.generate(1)
                assert_state(expected["before"])
                previous_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                previous_last_error=A.logger.get_last_error()
                @test A.mesh.add_homology_request(case["kind"],[11],Int[],case["dims"])===nothing
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==previous_errors
                @test A.logger.get_last_error()==previous_last_error
                if expected["status"]=="success"
                    result=A.mesh.compute_homology()
                    @test result isa Vector{Tuple{Int,Int}}
                    @test result==[(Int(row[1]),Int(row[2])) for row in expected["value"]]
                else
                    error=try A.mesh.compute_homology();nothing catch e;e isa InterruptException && rethrow();e end
                    @test error isa ErrorException
                    @test occursin(expected["exception"],sprint(showerror,error))
                end
                # The retained warmed native scalar cache can miss newly
                # published cells. Its unsafe query records are not expected
                # values. Check every safe state field and the complete owned
                # cell topology while retaining the unsafe scalar trace in
                # the sealed artifact instead of calling that query oracle.
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]-previous_errors==length(expected["error_events"])
                @test A.logger.get_last_error()==(isempty(expected["error_events"]) ?
                    previous_last_error : expected["last_error"])
                operation_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                operation_last_error=A.logger.get_last_error()
                actual=public_state(;include_scalar_cells=false)
                safe_keys=setdiff(Set(keys(expected["after"])),Set(["cell_details"]))
                @test Set(keys(actual))==safe_keys
                for key in safe_keys
                    @test actual[key]==expected["after"][key]
                end
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==operation_errors
                @test A.logger.get_last_error()==operation_last_error
            finally A.finalize() end
        end
    end
end

@testset "homology API Cint safety and strict Model queue validation remain atomic" begin
    A.initialize()
    try
        model=A._model_locked()
        for dims in ([true],[Int64(typemax(Cint))+1],[Int64(typemin(Cint))-1],[1.0])
            @test_throws ArgumentError A.mesh.add_homology_request("Homology",Int[],Int[],dims)
            @test isempty(model.meshing.homology_requests)
        end
        @test_throws ArgumentError A.mesh.add_homology_request(:Homology,Int[],Int[],[0])
        @test isempty(model.meshing.homology_requests)
        for (kind,dims) in (("Unknown",[0]),("Homology",[-1]),("Homology",[4]))
            @test_throws ArgumentError Tessella.Model.add_homology_request!(model;kind,dims)
            @test isempty(model.meshing.homology_requests)
        end
    finally A.finalize() end
end

@testset "native discrete entity type names in all dimensions" begin
    @test DISCRETE_TYPES["gmsh_version"]=="4.15.2" && length(DISCRETE_TYPES["cases"])==4
    for expected in DISCRETE_TYPES["cases"]
        A.initialize()
        try
            dim=expected["dimension"];tag=expected["tag"]
            @test A.model.add_discrete_entity(dim,tag)==tag
            @test A.model.get_entities()==[(dim,tag)]
            @test A.model.get_type(dim,tag)==expected["type"]
            @test Tessella.Model.model_entity_type(A._model_locked(),dim,tag)==expected["type"]
        finally A.finalize() end
    end
end
end
