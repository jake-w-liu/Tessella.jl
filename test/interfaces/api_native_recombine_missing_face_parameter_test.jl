module APINativeRecombineMissingFaceParameterTests
using Test,Tessella,TOML,SHA,Serialization
const A=Tessella.API
call(args...;kwargs...)=Base.invokelatest(args...;kwargs...)
function digest(value)
    io=IOBuffer();serialize(io,value);return bytes2hex(sha256(take!(io)))
end
function check_nodes(expected)
    for (label,record) in expected
        xyz,parameters,dimension,entity=call(A.mesh.get_node,parse(UInt64,label))
        @test reinterpret(UInt64,xyz)==reinterpret(UInt64,Float64.(record["xyz"]))
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
    call(A.mesh.add_nodes,2,1,[401],center["xyz"],Float64.(center["parameters"]))
    cells=only(control["before"]["elements"])
    call(A.mesh.add_elements_by_type,1,3,cells["tags"],cells["nodes"])
    return cells
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Native built-in Plane relocation with absent Face UV preserves empty raw storage" begin
    path=joinpath(@__DIR__,"..","artifacts","native_recombine_missing_face_parameter_primary.toml")
    raw=joinpath(@__DIR__,"..","artifacts","native_recombine_missing_face_parameter_primary.json")
    @test bytes2hex(sha256(read(path)))=="691adef81b12028174b229704f7fa3b97cc5edd3417dc420ceaf7dd462971242"
    @test bytes2hex(sha256(read(raw)))=="85988d3a6e0f8a6ada0855e0d2986dd92e30b2fae60044ff51adb2286b81dc58"
    fixture=TOML.parsefile(path)
    @test fixture["primary_json_sha256"]==bytes2hex(sha256(read(raw)))
    controls=filter(control->!control["face_parameters"],fixture["cases"])
    @test length(controls)==6
    for control in controls
        call(A.initialize)
        try
            for (tag,point) in enumerate(control["points"])
                call(A.model.add_point,point...;tag,meshSize=.5)
            end
            for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
            call(A.model.add_curve_loop,[1,2,3,4];tag=1)
            call(A.model.add_plane_surface,[1];tag=1)
            call(A.option,"Mesh.Renumber",0.)
            cells=insert_mesh(control)
            check_nodes(control["before"]["nodes"])
            old=A.CURRENT[];oldcache=A.LAST_MESH[];before=digest((old,oldcache))
            last_error=call(A.logger.get_last_error);errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
            word=Tessella.GmshLibm._win_x87_control_word()
            call(A.mesh.recombine)
            @test digest((old,oldcache))==before
            check_nodes(control["after"]["nodes"])
            @test call(A.mesh.get_elements,2,1)==([3],[UInt64.(cells["tags"])],[UInt64.(cells["nodes"])])
            @test call(A.mesh.get_max_node_tag)==401
            @test call(A.mesh.get_max_element_tag)==3004
            @test call(A.logger.get_last_error)==last_error
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
            @test Tessella.GmshLibm._win_x87_control_word()==word
            @test isempty(call(A.mesh.get_node,UInt64(401))[2])
        finally
            call(A.finalize)
        end
    end
end
@testset "Unsupported curved context with absent Face UV rejects atomically" begin
    control=first(filter(case->!case["face_parameters"],TOML.parsefile(joinpath(
        @__DIR__,"..","artifacts","native_recombine_missing_face_parameter_primary.toml"))["cases"]))
    mktempdir() do directory
        path=joinpath(directory,"missing_face_curved_boundary.geo")
        write(path,"Point(1)={0,0,0,.5};Point(2)={1,0,0,.5};Point(3)={1,1,0,.5};" *
            "Point(4)={0,1,0,.5};Point(5)={.5,.5,0,.5};Circle(1)={1,5,2};" *
            "Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};")
        call(A.initialize)
        try
            call(A.open_geo!,path;mesh_dim=0);call(A.option,"Mesh.Renumber",0.)
            insert_mesh(control)
            before=digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
            last_error=call(A.logger.get_last_error);errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
            word=Tessella.GmshLibm._win_x87_control_word()
            caught=try call(A.mesh.recombine);nothing catch error;error end
            @test caught isa ArgumentError
            @test caught isa ArgumentError && occursin("Curve evaluation context",sprint(showerror,caught))
            @test digest((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
            @test call(A.logger.get_last_error)==last_error
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
            @test Tessella.GmshLibm._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
    end
end
end
end
