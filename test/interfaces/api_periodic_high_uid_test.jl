module PeriodicHighUidTests
using Tessella,Test
const A=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_periodic_high_uid_v1.jl"))
bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
@testset "actual transfinite and periodic Surface births retain UInt64 node labels" begin
    expected=NATIVE_PERIODIC_HIGH_UID_V1
    A.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"periodic.geo");write(path,expected["source"])
            A.open_geo!(path;mesh_dim=0)
        end
        A.option("Mesh.Renumber",0)
        A.mesh.generate(1)
        old=A.mesh.get_nodes()[1]
        tags=UInt64.(2147483748:2147483748+length(old)-1)
        A.mesh.renumber_nodes(old,tags)
        A.option("Mesh.MeshOnlyEmpty",1)
        @test A.mesh.get_max_node_tag()==expected["before"]["node_max"]
        @test A.mesh.get_max_element_tag()==expected["before"]["element_max"]
        A.mesh.generate(2)
        @test A.mesh.get_max_node_tag()==expected["after"]["node_max"]
        @test A.mesh.get_max_element_tag()==expected["after"]["element_max"]
        for surface in expected["after"]["surfaces"]
            tag=surface["tag"]
            nodes,xyz,uv=A.mesh.get_nodes(2,tag,false,true)
            types,element_tags,corners=A.mesh.get_elements(2,tag)
            @test nodes==UInt64.(surface["nodes"])
            @test bits(xyz)==surface["xyz_bits"]
            @test bits(uv)==surface["uv_bits"]
            @test types==Int32.(surface["types"])
            @test element_tags==[UInt64.(block) for block in surface["elements"]]
            @test corners==[UInt64.(block) for block in surface["corners"]]
            for uid in nodes
                point,parameters,dim,owner=A.mesh.get_node(uid)
                @test (dim,owner)==(2,tag)
                @test bits(point)==surface["xyz_bits"]
                @test bits(parameters)==surface["uv_bits"]
            end
        end
    finally
        A.finalize()
    end
end
end
