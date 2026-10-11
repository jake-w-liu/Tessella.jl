module APINativeRecombinePlaneContextTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
function check_nodes(expected)
    for (label,record) in expected
        point,parameters,dimension,entity=call(A.mesh.get_node,parse(UInt64,label))
        @test reinterpret(UInt64,point)==reinterpret(UInt64,record["xyz"])
        @test reinterpret(UInt64,parameters)==reinterpret(UInt64,Float64.(record["parameters"]))
        @test (dimension,entity)==(record["dimension"],record["entity"])
    end
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Pinned native general Plane default Surface optimization" begin
    fixture=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","native_recombine_general_plane_primary.toml"))
    @test fixture["version"]=="4.15.2"
    for control in fixture["cases"]
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        try
            for (tag,point) in enumerate(control["points"])
                call(A.model.add_point,point...;tag,meshSize=.5)
            end
            for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
            call(A.model.add_curve_loop,control["curve_loop"];tag=1)
            call(A.model.add_plane_surface,[1];tag=1)
            for (name,value) in (("Mesh.Renumber",0.),("Mesh.RecombinationAlgorithm",1.),
                                ("Mesh.RecombineNodeRepositioning",1.),("Mesh.RecombineOptimizeTopology",5.),
                                ("Mesh.Smoothing",1.),("Mesh.RecombineMinimumQuality",.01))
                call(A.option,name,value)
            end
            # Stored Surface/Curve parameter context is distinct from actual
            # mesh XYZ. CAD supplies the boundary reparametrization.
            for (tag,point) in enumerate(control["points"])
                call(A.mesh.add_nodes,0,tag,[100+tag],point)
            end
            for (curve,boundary) in enumerate(control["boundaries"])
                call(A.mesh.add_nodes,1,curve,[200+curve],boundary["point"],[.5])
                call(A.mesh.add_elements_by_type,curve,1,[2000+curve,2100+curve],
                    [100+curve,200+curve,200+curve,100+mod1(curve+1,4)])
            end
            center=control["before"]["nodes"]["401"]
            call(A.mesh.add_nodes,2,1,[401],center["xyz"],center["parameters"])
            cells=only(control["before"]["elements"])
            call(A.mesh.add_elements_by_type,1,3,cells["tags"],cells["nodes"])
            check_nodes(control["before"]["nodes"])
            @test call(A.mesh.get_max_node_tag)==401
            @test call(A.mesh.get_max_element_tag)==3004
            old=A.CURRENT[];oldcache=A.LAST_MESH[];before=digest((old,oldcache))
            word=Tessella.GmshLibm._win_x87_control_word()
            call(A.mesh.recombine)
            @test digest((old,oldcache))==before
            check_nodes(control["after"]["nodes"])
            @test call(A.mesh.get_elements,2,1)==([3],[UInt64.(cells["tags"])],[UInt64.(cells["nodes"])])
            @test call(A.mesh.get_max_node_tag)==401
            @test call(A.mesh.get_max_element_tag)==3004
            @test Tessella.GmshLibm._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
    end
end
end
end
