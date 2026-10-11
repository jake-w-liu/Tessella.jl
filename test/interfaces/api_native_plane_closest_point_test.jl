module APINativePlaneClosestPointTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
const G=Tessella.GmshLibm
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Native built-in Plane closest point: 112 frames and twelve caller modes" begin
    path=joinpath(@__DIR__,"..","artifacts","native_plane_closest_point_primary.toml")
    raw=joinpath(@__DIR__,"..","artifacts","native_plane_closest_point_primary.json")
    @test bytes2hex(sha256(read(path)))=="98e96e3273243da24959515c2151277a04f3a4cbec55d983e1d02e291f67c257"
    @test bytes2hex(sha256(read(raw)))=="f3705045ce1d9af50254b94e9f153fc35eaf784bbb145c46112630a859206d69"
    fixture=TOML.parsefile(path)
    @test fixture["primary_json_sha256"]==bytes2hex(sha256(read(raw)))
    @test fixture["original_primary_sha256"]=="fd307e834b654ec9cf44b106873cc050c98528620d86cac2c80195deacd6480d"
    @test fixture["version"]=="4.15.2"
    @test length(fixture["cases"])==112
    @test length(fixture["modes"])==12
    context_path=joinpath(@__DIR__,"..","artifacts","native_plane_constructor_context_primary.toml")
    @test bytes2hex(sha256(read(context_path)))=="93b03d0e51d7c31a3d295df7019c84c32f0ab7b7651d86ab3b8000783dfa9cf7"
    contexts=TOML.parsefile(context_path)["cases"]
    queries=0;calls=0
    call(A.initialize)
    saved=G._win_x87_control_word()
    try
      for (case_index,control) in enumerate(fixture["cases"])
        try
            context=contexts[case_index]
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,context["lc_before_clear"])
            call(A.clear)
            call(A.model.add,"closest_point_control")
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,context["lc_after_clear"])
            for (tag,point) in enumerate(control["points"])
                call(A.model.geo.add_point,point...;tag,meshSize=.5)
            end
            for tag in 1:4;call(A.model.geo.add_line,tag,mod1(tag+1,4);tag);end
            call(A.model.geo.add_curve_loop,control["curve_loop"];tag=1)
            call(A.model.geo.add_plane_surface,[1];tag=1)
            @test isempty(call(A.model.get_entities))
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,context["lc_before_surface"])
            call(A.model.geo.synchronize)
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,context["lc_after_sync"])
            @test bits(call(A.model.get_value,2,1,[0.,0.]))==bits(control["origin"])
            @test bits(call(A.model.get_derivative,2,1,[0.,0.]))==bits(vcat(control["first"],control["second"]))
            lower,upper=call(A.model.get_parametrization_bounds,2,1)
            @test bits(lower)==bits(control["lower"])
            @test bits(upper)==bits(control["upper"])
            bound_state=digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],
                A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
            last_error=call(A.logger.get_last_error);lookup_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
            for query in control["queries"]
                queries+=1
                for mode in fixture["modes"]
                    caller=(saved&UInt16(0xf0ff))|UInt16(mode)
                    G._win_x87_set_control_word!(caller)
                    xyz,uv=call(A.model.get_closest_point,2,1,query["query"])
                    @test bits(xyz)==bits(query["xyz"])
                    @test bits(uv)==bits(query["uv"])
                    @test G._win_x87_control_word()==caller
                    calls+=1
                end
            end
            G._win_x87_set_control_word!(saved)
            @test digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],
                A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==bound_state
            @test call(A.logger.get_last_error)==last_error
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==lookup_errors
        finally
            G._win_x87_set_control_word!(saved)
        end
      end
    finally
        G._win_x87_set_control_word!(saved)
        call(A.finalize)
    end
    @test queries==fixture["unique_queries"]==2016
    @test calls==fixture["scalar_mode_calls"]==24192
end
end
end
