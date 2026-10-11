module NativeTranslatedLineLateralTests
using Test,Tessella,TOML,SHA
const API=Tessella.API
const ARTIFACTS=joinpath(@__DIR__,"..","artifacts")
bits(values)=[bytes2hex(reinterpret(UInt8,[Float64(value)])) for value in values]

@testset "Native translated signed Line lateral point and derivatives" begin
    @test bytes2hex(sha256(read(joinpath(ARTIFACTS,"native_translated_line_lateral_values_v1.json"))))==
        "96eac72d3976f83735b0be7e642a5217a7bcb8e5be6d17aa8cda51489be00bfc"
    @test bytes2hex(sha256(read(joinpath(ARTIFACTS,"native_translated_line_lateral_values_v1.toml"))))==
        "2c25f9ed30bd8643a33eb245e64ca2552de650836d26c120a8e621ee7a913b30"
    native=TOML.parsefile(joinpath(ARTIFACTS,"native_translated_line_lateral_values_v1.toml"))
    @test native["native_json_sha256"]=="96eac72d3976f83735b0be7e642a5217a7bcb8e5be6d17aa8cda51489be00bfc"
    API.initialize()
    try
        mktempdir() do directory
            for (index,row) in pairs(native["cases"])
                API.clear()
                first=row["first"];last=row["last"];delta=row["delta"];sign=row["source_sign"]
                source="""
                Geometry.AutoCoherence=0;
                Point(1)={$(join(first,','))};Point(2)={$(join(last,','))};
                Line(1)={1,2};
                s[]=Extrude{$(join(delta,','))}{Curve{$sign};Layers{2};};
                """
                path=joinpath(directory,"lateral_$index.geo");write(path,source)
                API.open_geo!(path;mesh_dim=0)
                surface=only(entity[2] for entity in API.model.get_entities(2))
                before=deepcopy((API.CURRENT[].points,API.CURRENT[].curves,
                    API.CURRENT[].loops,API.CURRENT[].surfaces,API.CURRENT[].next_tag))
                for query in row["queries"]
                    uv=Float64.(query["uv"])
                    @test bits(API.model.get_value(2,surface,uv))==query["xyz_bits"]
                    @test bits(API.model.get_derivative(2,surface,uv))==query["d1_bits"]
                end
                @test isequal(before,(API.CURRENT[].points,API.CURRENT[].curves,
                    API.CURRENT[].loops,API.CURRENT[].surfaces,API.CURRENT[].next_tag))
                @test API.LAST_MESH[]===nothing
            end
        end
    finally
        API.finalize()
    end
end
end
