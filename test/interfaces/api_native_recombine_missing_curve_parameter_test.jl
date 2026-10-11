module APINativeRecombineMissingCurveParameterTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
function check_nodes(expected)
    for (label,record) in expected
        point,parameters,dimension,entity=call(A.mesh.get_node,parse(UInt64,label))
        @test reinterpret(UInt64,point)==reinterpret(UInt64,Float64.(record["xyz"]))
        @test reinterpret(UInt64,parameters)==reinterpret(UInt64,Float64.(record["parameters"]))
        @test (dimension,entity)==(record["dimension"],record["entity"])
    end
end
function insert_mesh(control)
    for (tag,point) in enumerate(control["points"])
        call(A.mesh.add_nodes,0,tag,[100+tag],point)
    end
    for curve in 1:4
        record=control["before"]["nodes"][string(200+curve)]
        call(A.mesh.add_nodes,1,curve,[200+curve],record["xyz"],Float64.(record["parameters"]))
        call(A.mesh.add_elements_by_type,curve,1,[2000+curve,2100+curve],
            [100+curve,200+curve,200+curve,100+mod1(curve+1,4)])
    end
    center=control["before"]["nodes"]["401"]
    call(A.mesh.add_nodes,2,1,[401],center["xyz"],center["parameters"])
    cells=only(control["before"]["elements"])
    call(A.mesh.add_elements_by_type,1,3,cells["tags"],cells["nodes"])
    return cells
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Pinned native missing straight-Line parameters preserve raw storage" begin
    artifact_path=joinpath(@__DIR__,"..","artifacts",
        "native_recombine_missing_curve_parameter_primary.toml")
    raw_path=joinpath(@__DIR__,"..","artifacts",
        "native_recombine_missing_curve_parameter_primary.json")
    @test bytes2hex(sha256(read(artifact_path)))=="36b9c18b9ce6538122c886031f4ab942766939564615bc43e3d563fa083a01f4"
    @test bytes2hex(sha256(read(raw_path)))=="7a4b5dfd56e5ae8f3ccaee28db638558548692cade79acc8b5dc7c92ef6f6a05"
    fixture=TOML.parsefile(artifact_path)
    @test fixture["primary_json_sha256"]=="7a4b5dfd56e5ae8f3ccaee28db638558548692cade79acc8b5dc7c92ef6f6a05"
    @test fixture["version"]=="4.15.2"
    @test length(fixture["cases"])==6
    for control in fixture["cases"]
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        try
            for (tag,point) in enumerate(control["points"])
                call(A.model.add_point,point...;tag,meshSize=.5)
            end
            for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
            call(A.model.add_curve_loop,[1,2,3,4];tag=1)
            call(A.model.add_plane_surface,[1];tag=1)
            call(A.option,"Mesh.Renumber",0.)
            for (name,expected) in (("Mesh.RecombinationAlgorithm",1.),
                ("Mesh.RecombineNodeRepositioning",1.),("Mesh.RecombineOptimizeTopology",5.),
                ("Mesh.Smoothing",1.),("Mesh.RecombineMinimumQuality",.01))
                @test call(A.option,name)==expected
            end
            cells=insert_mesh(control)
            check_nodes(control["before"]["nodes"])
            for phase in ("after","after_repeat")
                old=A.CURRENT[];oldcache=A.LAST_MESH[];before=digest((old,oldcache))
                last_error=call(A.logger.get_last_error)
                errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                word=Tessella.GmshLibm._win_x87_control_word()
                call(A.mesh.recombine)
                @test digest((old,oldcache))==before
                check_nodes(control[phase]["nodes"])
                @test call(A.mesh.get_elements,2,1)==([3],[UInt64.(cells["tags"])],[UInt64.(cells["nodes"])])
                @test call(A.mesh.get_max_node_tag)==401
                @test call(A.mesh.get_max_element_tag)==3004
                @test call(A.logger.get_last_error)==last_error
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
                @test Tessella.GmshLibm._win_x87_control_word()==word
            end
        finally
            call(A.finalize)
        end
    end
end

@testset "Unsupported curved boundary remains atomic with absent Curve parameters" begin
    control=first(filter(case->!case["curve_parameters"],
        TOML.parsefile(joinpath(@__DIR__,"..","artifacts",
            "native_recombine_missing_curve_parameter_primary.toml"))["cases"]))
    # Circle is a real GEO CAD boundary, outside this straight-Line increment.
    # The mesh has classified, positive cells; the bounded adapter must reject
    # its unsupported evaluator before publishing XYZ, UV, labels or counters.
    mktempdir() do directory
        path=joinpath(directory,"curved_boundary.geo")
        write(path,"Point(1)={0,0,0,.5}; Point(2)={1,0,0,.5}; " *
            "Point(3)={1,1,0,.5}; Point(4)={0,1,0,.5}; Point(5)={.5,.5,0,.5}; " *
            "Circle(1)={1,5,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1}; " *
            "Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};")
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        try
            call(A.open_geo!,path;mesh_dim=0)
            call(A.option,"Mesh.Renumber",0.)
            insert_mesh(control)
            before=digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],
                A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
            last_error=call(A.logger.get_last_error);errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
            caught=try call(A.mesh.recombine);nothing catch error;error end
            @test caught isa ArgumentError
            @test caught isa ArgumentError && occursin("Curve evaluation context",sprint(showerror,caught))
            @test digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],
                A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
            @test call(A.logger.get_last_error)==last_error
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
        finally
            call(A.finalize)
        end
    end
end
end
end
