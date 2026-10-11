module APINativePlaneResynchronizeTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
const G=Tessella.GmshLibm
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Native GEO repeat synchronization and cached Plane bounds" begin
    path=joinpath(@__DIR__,"..","artifacts","native_plane_resynchronize_primary.toml")
    raw=joinpath(@__DIR__,"..","artifacts","native_plane_resynchronize_primary.json")
    @test bytes2hex(sha256(read(path)))=="a31b1447b21e26b8b658f977bfb453fb4d313c5015500574c3931dc8164f8a84"
    @test bytes2hex(sha256(read(raw)))=="ab1b548e88153ac1a96d905ef98c1034ee8111bb377694b1809ca91310ecf539"
    fixture=TOML.parsefile(path)
    @test fixture["primary_json_sha256"]==bytes2hex(sha256(read(raw)))
    primary=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","native_plane_closest_point_primary.toml"))
    @test fixture["primary_sha256"]==primary["original_primary_sha256"]
    call(A.initialize)
    saved=G._win_x87_control_word();calls=0
    try
        for control in primary["cases"][1:fixture["prefix_cases"]]
            call(A.clear);call(A.model.add,"closest_point_control")
            for (tag,point) in enumerate(control["points"])
                call(A.model.geo.add_point,point...;tag,meshSize=.5)
            end
            for tag in 1:4;call(A.model.geo.add_line,tag,mod1(tag+1,4);tag);end
            call(A.model.geo.add_curve_loop,control["curve_loop"];tag=1)
            call(A.model.geo.add_plane_surface,[1];tag=1)
            call(A.model.geo.synchronize)
        end
        for row in fixture["rows"]
            phase=row["phase"]
            if phase=="repeat_sync" || phase=="after_point_sync"
                call(A.model.geo.synchronize)
            elseif phase=="point_set_coordinates"
                call(A.model.set_coordinates,fixture["edited_point"],fixture["edited_coordinates"]...)
            else
                @test phase=="first_sync"
            end
            control=row["state"]
            @test bits(call(A.model.get_value,2,1,[0.,0.]))==bits(control["origin"])
            @test bits(call(A.model.get_derivative,2,1,[0.,0.]))==bits(vcat(control["first"],control["second"]))
            lower,upper=call(A.model.get_parametrization_bounds,2,1)
            @test bits(lower)==bits(control["bounds"][1])
            @test bits(upper)==bits(control["bounds"][2])
            @test reinterpret(UInt64,call(A.option,"General.BoundingBoxSize"))==reinterpret(UInt64,control["lc"])
            before=digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
            last_error=call(A.logger.get_last_error);lookup_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
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
            @test digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
            @test call(A.logger.get_last_error)==last_error
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==lookup_errors
        end
    finally
        G._win_x87_set_control_word!(saved);call(A.finalize)
    end
    @test calls==fixture["calls"]==864
end
end
end
