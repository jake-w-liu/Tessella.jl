module APICurvedProductViewsNativeTests
using Tessella,Test,TOML,SHA
include(joinpath(@__DIR__,"..","..","validation","native_geometry",
    "native_periodic_surface_production_controls_v1.jl"))
const C=NativePeriodicSurfaceProductionControlsV1
const A=Tessella.API
@testset "curved GEO completed products retain exact native public state" begin
    artifact=joinpath(@__DIR__,"..","artifacts","curved_product_views_native_v1.toml")
    packet=TOML.parsefile(artifact)
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "0cfbdd24d89e4047724990529c2c7a5d015c7f980f667c80ff35710cd61b9552"
    @test packet["gmsh_version"]=="4.15.2" && packet["all_three_full_states_exact"]
    @test length(packet["raw_sha256"])==3
    script=joinpath(@__DIR__,"..","..",packet["capture_script"])
    @test bytes2hex(sha256(replace(read(script,String),"\r\n"=>"\n")))==packet["capture_script_sha256"]
    A.initialize()
    try
        mktempdir() do directory
            source=joinpath(directory,"scene.geo");write(source,packet["scene"])
            @test bytes2hex(sha256(read(source)))==packet["source_sha256"]
            A.open_geo!(source;mesh_dim=0)
        end
        A.mesh.generate(2)
        # This covers every entity, directed cell corner, UID, owner, XYZ/UV
        # bit, scalar query, physical group, visibility, and factory clock.
        @test C.exact(C.state(),packet["native"])
        @test length(A.mesh.get_elements(2,1)[2][1])==104
        @test length(A.mesh.get_elements(2,3)[2][1])==11
    finally A.finalize() end
end
end
