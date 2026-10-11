module NativeInitialSurfaceBdsTests
using Tessella,Test,SHA
const fixture_path=joinpath(@__DIR__,"..","..","validation","native_geometry",
    "native_initial_surface_bds_controls_v1.jl")
include(fixture_path)
const A=Tessella.API
bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
@testset "actual initial-only BDS partition traversal and constructor births" begin
    @test bytes2hex(sha256(read(fixture_path)))==
        "8a2d71b389aa8de454f76f36d91c04b38a0d05c91cd295a086bccb313ba51f40"
    @test NATIVE_INITIAL_SURFACE_BDS_CONTROLS_V1["native_version"]=="4.15.2"
    for case in NATIVE_INITIAL_SURFACE_BDS_CONTROLS_V1["cases"]
        @testset "$(case["label"])" begin
            A.initialize()
            try
                mktempdir() do directory
                    path=joinpath(directory,"scene.geo");write(path,case["scene"])
                    @test bytes2hex(sha256(read(path)))==case["scene_sha256"]
                    A.open_geo!(path;mesh_dim=0)
                end
                A.mesh.generate(2)
                @test A.mesh.get_max_node_tag()==case["node_max"]
                @test A.mesh.get_max_element_tag()==case["element_max"]
                for owner in case["owners"]
                    nodes,xyz,uv=A.mesh.get_nodes(owner["owner"]...,false,true)
                    @test nodes==UInt64.(owner["uids"])
                    @test bits(xyz)==owner["xyz_bits"]
                    @test bits(uv)==owner["uv_bits"]
                    kinds,tags,corners=A.mesh.get_elements(owner["owner"]...)
                    @test kinds==Int32.(owner["types"])
                    @test tags==[UInt64.(block) for block in owner["tags"]]
                    @test corners==[UInt64.(block) for block in owner["corners"]]
                end
                println("INITIAL_NEIGHBOR_DONE ",case["label"]);flush(stdout)
            finally A.finalize() end
        end
    end
end

end
