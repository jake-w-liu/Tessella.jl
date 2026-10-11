module APINativePlaneSVDVariedContextTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
const G=Tessella.GmshLibm
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Native Plane single RHS QR reduction across 160 varied contexts" begin
    path=joinpath(@__DIR__,"..","artifacts","native_plane_svd_varied_context_v1.toml")
    raw=joinpath(@__DIR__,"..","artifacts","native_plane_svd_varied_context_v1.json")
    @test bytes2hex(sha256(read(path)))=="a4c3ec6f52767ed261401c7942ff3524b5d512226a72a96d1523ed61e87be13c"
    @test bytes2hex(sha256(read(raw)))=="ae7c3512d062700275003389681e93d698fa8649faabc56b79eb2b33fb975744"
    fixture=TOML.parsefile(path)
    @test length(fixture["cases"])==160
    call(A.initialize);saved=G._win_x87_control_word()
    try
        for (index,control) in enumerate(fixture["cases"])
            @testset "native varied Plane $index" begin
                call(A.clear);call(A.model.add,"seed")
                call(A.model.geo.add_point,-1.e15,-1.e15,-1.e15;tag=1,meshSize=1.)
                call(A.model.geo.add_point,1.e15,1.e15,1.e15;tag=2,meshSize=1.)
                call(A.model.geo.synchronize)
                call(A.clear);call(A.model.add,"fallback")
                for (tag,point) in enumerate(control["points"])
                    call(A.model.geo.add_point,point...;tag,meshSize=1.)
                end
                for tag in 1:4;call(A.model.geo.add_line,tag,mod1(tag+1,4);tag);end
                call(A.model.geo.add_curve_loop,control["loop"];tag=1)
                @test bits([call(A.option,"General.BoundingBoxSize")])==bits([control["lc_before_surface"]])
                call(A.model.geo.add_plane_surface,[1];tag=1)
                call(A.model.geo.synchronize)
                samples=[call(A.model.get_value,1,abs(curve),[t]) for curve in control["loop"] for t in (.333,.666)]
                @test bits(reduce(vcat,samples))==bits(reduce(vcat,control["samples"]))
                @test Tessella.Model._gmsh_native_plane_frame(
                    Tuple.(samples),1.e-8,control["lc_before_surface"])===nothing
                before=digest((A.CURRENT[],A.LAST_MESH[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
                @test bits(call(A.model.get_value,2,1,[0.,0.]))==bits(control["expected_origin"])
                @test bits(call(A.model.get_derivative,2,1,[0.,0.]))==
                    bits(vcat(control["expected_first"],control["expected_second"]))
                @test digest((A.CURRENT[],A.LAST_MESH[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
                @test G._win_x87_control_word()==saved
            end
        end
    finally
        G._win_x87_set_control_word!(saved);call(A.finalize)
    end
end
end
end
