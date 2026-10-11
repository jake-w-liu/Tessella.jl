module APINativePlaneFreshContextTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
const G=Tessella.GmshLibm
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Native Plane constructors in 112 separate fresh sessions" begin
    path=joinpath(@__DIR__,"..","artifacts","native_plane_fresh_context_primary.toml")
    raw=joinpath(@__DIR__,"..","artifacts","native_plane_fresh_context_primary.json")
    @test bytes2hex(sha256(read(path)))=="76f37ea44c342c7a7d8a8ee57f01935e5a6a040ecbde4d51c574d976c3139d29"
    @test bytes2hex(sha256(read(raw)))=="0ae6fad2a2b9dd2230d3cd0decad4a15561a6c1159e6757e513930452f8185e2"
    fixture=TOML.parsefile(path)
    @test fixture["primary_json_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(fixture["cases"])==112
    calls=0
    for control in fixture["cases"]
        call(A.initialize)
        saved=G._win_x87_control_word()
        try
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,control["lc_before"])
            call(A.clear);call(A.model.add,"closest_point_control")
            for (tag,point) in enumerate(control["points"])
                call(A.model.geo.add_point,point...;tag,meshSize=.5)
            end
            for tag in 1:4;call(A.model.geo.add_line,tag,mod1(tag+1,4);tag);end
            call(A.model.geo.add_curve_loop,control["curve_loop"];tag=1)
            call(A.model.geo.add_plane_surface,[1];tag=1)
            @test isempty(call(A.model.get_entities))
            call(A.model.geo.synchronize)
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,control["lc_after"])
            @test bits(call(A.model.get_value,2,1,[0.,0.]))==bits(control["origin"])
            @test bits(call(A.model.get_derivative,2,1,[0.,0.]))==bits(vcat(control["first"],control["second"]))
            lower,upper=call(A.model.get_parametrization_bounds,2,1)
            @test bits(lower)==bits(control["lower"])
            @test bits(upper)==bits(control["upper"])
            before=digest((A.CURRENT[],A.LAST_MESH[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
            for query in control["queries"],mode in fixture["modes"]
                caller=(saved&UInt16(0xf0ff))|UInt16(mode)
                G._win_x87_set_control_word!(caller)
                xyz,uv=call(A.model.get_closest_point,2,1,query["query"])
                @test bits(xyz)==bits(query["xyz"])
                @test bits(uv)==bits(query["uv"])
                @test G._win_x87_control_word()==caller
                calls+=1
            end
            G._win_x87_set_control_word!(saved)
            @test digest((A.CURRENT[],A.LAST_MESH[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
        finally
            G._win_x87_set_control_word!(saved);call(A.finalize)
        end
    end
    @test calls==fixture["calls"]==24192
end
end
end
