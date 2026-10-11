module PlaneNativeNormalQueryControlsV1
using Test,Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_plane_normal_query_fixtures_v2.jl"))
const M=Tessella.Model
_bits(values)=string.(reinterpret.(UInt64,Float64.(values));base=16,pad=16)
function _positive_box_control()
    captured=0
    probe=()->captured
    captured+=1
    return probe
end
function _has_box(value)
    value===Core.Box && return true
    value isa Core.Box && return true
    value isa Expr && return any(_has_box,value.args)
    value isa QuoteNode && return _has_box(value.value)
    value isa Core.CodeInfo && return any(_has_box,value.code)
    value isa AbstractArray && return any(_has_box,value)
    return false
end
function _lowered_body_check(function_object,signature,path)
    method=which(function_object,signature)
    @test normpath(String(method.file))==normpath(path)
    lowered=Base.code_lowered(function_object,signature)
    @test !isempty(lowered)
    @test !_has_box(lowered)
    return nothing
end
function run!()
    @test _has_box(Base.code_lowered(_positive_box_control,Tuple{}))
    reader=_positive_box_control()
    @test reader()==1
    @test any(type->type===Core.Box,fieldtypes(typeof(reader)))
    @test length(_NATIVE_PLANE_NORMAL_QUERY_ROWS_V2)==26
    expected_origin=joinpath(Base.pkgdir(Tessella),"src","geometry","ModelEntityEvaluation.jl")
    for row in _NATIVE_PLANE_NORMAL_QUERY_ROWS_V2
        @testset "$(row["name"])" begin
            execution=mktempdir() do directory
                path=joinpath(directory,"plane.geo");write(path,row["source"])
                execute_geo(path;mesh_dim=0)
            end
            model=execution.model
            saved=deepcopy(model)
            if row["degenerate"]
                # Native accepts an SVD fallback here. This preexisting Model
                # fallback is explicitly unsupported; retain its error policy
                # and native packet rather than fabricate a native normal.
                @test_throws ArgumentError M.model_normal(model,1,[0.,0.])
                @test all(name->isequal(getfield(model,name),getfield(saved,name)),fieldnames(typeof(model)))
                continue
            end
            plane=M._model_plane_frame(model,1,"normal controls")
            for query in row["queries"]
                parameters=Float64.(query["parameters"])
                warnings=String[]
                actual=M.model_normal(model,1,parameters;warn=message->push!(warnings,message))
                @test _bits(actual)==query["result"]["normal_bits"]
                @test _bits(M.model_value(model,2,1,parameters))==query["result"]["value_bits"]
                @test _bits(M.model_derivative(model,2,1,parameters))==query["result"]["derivative_bits"]
                @test _bits(M.model_second_derivative(model,2,1,parameters))==query["result"]["second_derivative_bits"]
                expected_warnings=String[message[10:end] for message in query["messages"] if startswith(message,"Warning: ")]
                @test warnings==expected_warnings
                @test M._model_plane_frame(model,1,"normal controls")==plane
                @test all(name->isequal(getfield(model,name),getfield(saved,name)),fieldnames(typeof(model)))
            end
            _lowered_body_check(M._model_native_plane_query_normal,
                Tuple{typeof(model),Int,typeof(plane),Float64,Float64,String,Nothing},expected_origin)
            _lowered_body_check(M._model_native_plane_angle,
                Tuple{NTuple{3,Float64},NTuple{3,Float64},NTuple{3,Float64},NTuple{3,Float64}},expected_origin)
            _lowered_body_check(M.model_normal,Tuple{typeof(model),Int,Vector{Float64}},expected_origin)
            _lowered_body_check(M._gmsh_native_plane_declared_triple,
                Tuple{typeof(plane.coordinates),typeof(plane.projection)},
                joinpath(Base.pkgdir(Tessella),"src","geometry","ModelGmshSurfaceProducer.jl"))
        end
    end
end
run!()
end
