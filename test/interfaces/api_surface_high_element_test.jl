module SurfaceHighElementTests
using Tessella,Test,SHA
const fixture_path=joinpath(@__DIR__,"..","..","validation","native_geometry",
    "native_surface_high_element_v1.jl")
include(fixture_path)
const A=Tessella.API
surface_high_element_bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
@testset "native Surface UInt64 constructor labels" begin
    @test bytes2hex(sha256(read(fixture_path)))==
        "d3affc8a1b4f74e095aa7a2efa769861ae6728fc1007387b682d7c0d71b68225"
    @test NATIVE_SURFACE_HIGH_ELEMENT_V1["native_version"]=="4.15.2"
    for case in NATIVE_SURFACE_HIGH_ELEMENT_V1["cases"]
        @testset "$(case["label"])" begin
            A.initialize()
            try
                mktempdir() do directory
                    path=joinpath(directory,"scene.geo")
                    # The sealed native child opened this Windows text file.
                    write(path,replace(case["scene"],"\n"=>"\r\n"))
                    @test bytes2hex(sha256(read(path)))==case["scene_sha256"]
                    A.open_geo!(path;mesh_dim=0)
                end
                A.option("Mesh.Renumber",0);A.mesh.generate(1)
                old=reduce(vcat,A.mesh.get_elements()[2])
                point=only(only(A.mesh.get_elements(0,1)[2]))
                replacement=[uid==point ? UInt64(case["before"]["element_max"]) : uid for uid in old]
                A.mesh.renumber_elements(old,replacement)
                A.option("Mesh.MeshOnlyEmpty",1)
                @test A.mesh.get_max_node_tag()==case["before"]["node_max"]
                @test A.mesh.get_max_element_tag()==case["before"]["element_max"]
                A.mesh.generate(2)
                expected=case["after"]
                @test A.mesh.get_max_node_tag()==expected["node_max"]
                @test A.mesh.get_max_element_tag()==expected["element_max"]
                for owner in expected["nodes"]
                    nodes,xyz,uv=A.mesh.get_nodes(owner["owner"]...,false,true)
                    @test nodes==UInt64.(owner["tags"])
                    @test surface_high_element_bits(xyz)==owner["xyz_bits"]
                    @test surface_high_element_bits(uv)==owner["uv_bits"]
                end
                types,tags,corners=A.mesh.get_elements(2,1)
                @test types==Int32.(expected["cells"]["types"])
                @test tags==[UInt64.(block) for block in expected["cells"]["tags"]]
                @test corners==[UInt64.(block) for block in expected["cells"]["corners"]]
                for cell in expected["scalar"]
                    kind,conn,dim,tag=A.mesh.get_element(cell["tag"])
                    @test (kind,conn,dim,tag)==(cell["type"],UInt64.(cell["corners"]),cell["owner"]...)
                end
                println("SURFACE_HIGH_ELEMENT_DONE ",case["label"])
                flush(stdout)
            finally
                A.finalize()
            end
        end
    end
end
end
