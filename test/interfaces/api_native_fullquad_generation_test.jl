module APINativeFullQuadGenerationTests
using Test,Tessella,TOML,SHA,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const RAW=joinpath(@__DIR__,"../artifacts/native_fullquad_generation_primary.json")
const PROJECTION=joinpath(@__DIR__,"../artifacts/native_fullquad_generation_primary.toml")
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
function start()
    call(A.initialize)
    call(A.mesh.generate,0)
    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
end
function fixture(case)
    points=case["shape"]=="square" ? [(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)] :
        [(0.,0.,0.),(1.,0.,0.),(0.,1.,0.)]
    case["tilted"] && (points=[(x,y,x+2y) for (x,y,_) in points])
    mesh_size=get(case,"mesh_size",10.)
    for (tag,point) in enumerate(points)
        call(A.model.add_point,point...;tag,meshSize=mesh_size)
    end
    for tag in eachindex(points)
        call(A.model.add_line,tag,mod1(tag+1,length(points));tag)
    end
    call(A.model.add_curve_loop,collect(eachindex(points));tag=1)
    call(A.model.add_plane_surface,[1];tag=1)
    curve_nodes=get(case,"curve_nodes",3)
    if curve_nodes>0
        for tag in eachindex(points);call(A.mesh.set_transfinite_curve,tag,curve_nodes);end
    end
    call(A.mesh.set_recombine,2,1)
    for (name,value) in (("Mesh.Renumber",0),("Mesh.Algorithm",case["mesh_algorithm"]),
            ("Mesh.RecombinationAlgorithm",case["recombine_algorithm"]),
            ("Mesh.SecondOrderLinear",case["second_order_linear"] ? 1 : 0),
            ("Mesh.RecombineNodeRepositioning",case["reposition"]),
            ("Mesh.RecombineOptimizeTopology",case["topology"]),
            ("Mesh.Smoothing",case["smoothing"]),("Mesh.MeshSizeMin",mesh_size),
            ("Mesh.MeshSizeMax",mesh_size))
        call(A.option,name,value)
    end
end
@testset "native full-quad public generation" begin
@testset "native full-quad generation 2/3 exact directed products and raw vertices" begin
    @test bytes2hex(sha256(read(RAW)))=="4441e0362da948d65e8aaf4505a7c0fe767be039de124886558d526023879be4"
    @test bytes2hex(sha256(read(PROJECTION)))=="444c7722ac208ceeaee7c77584f134c3488445feccbfbc1a2e370f5b74e70203"
    corpus=TOML.parsefile(PROJECTION)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(RAW)))
    @test corpus["version"]=="4.15.2"
    @test length(corpus["cases"])==48
    for case in corpus["cases"]
        @testset "$(case["shape"]) mesher=$(case["mesh_algorithm"]) selector=$(case["recombine_algorithm"]) tilted=$(case["tilted"]) linear=$(case["second_order_linear"]) size=$(case["mesh_size"])" begin
        start()
        try
            fixture(case)
            model=A.CURRENT[];cache=A.LAST_MESH[]
            before=fingerprint((model,cache))
            word=Tessella.GmshLibm._win_x87_control_word()
            call(A.mesh.generate,2)
            @test fingerprint((model,cache))==before
            @test call(A.mesh.get_elements)==(case["types"],
                [UInt64.(row) for row in case["tags"]],
                [UInt64.(row) for row in case["cells"]])
            @test Set(first(call(A.mesh.get_nodes)))==Set(UInt64(row["tag"]) for row in case["nodes"])
            for row in case["nodes"]
                point,parameter,dimension,owner=call(A.mesh.get_node,row["tag"])
                @test bits(point)==bits(row["xyz"])
                @test bits(parameter)==bits(row["parameters"])
                @test (dimension,owner)==Tuple(row["owner"])
            end
            @test call(A.mesh.get_max_node_tag)==case["max_node"]
            @test call(A.mesh.get_max_element_tag)==case["max_element"]
            @test call(A.logger.get_last_error)==case["last_error"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            @test call(A.option,"Mesh.MeshSizeFactor")==1
            @test call(A.option,"Mesh.SecondOrderLinear")==Int(case["second_order_linear"])
            @test A.CURRENT[].meshing.second_order_linear==case["second_order_linear"]
            @test Tessella.GmshLibm._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
        end
    end
end

@testset "full-quad Plane mesh-periodic masters restore Curve identities before copying" begin
    raw=joinpath(@__DIR__,"../artifacts/native_fullquad_periodic_primary.json")
    projection=joinpath(@__DIR__,"../artifacts/native_fullquad_periodic_primary.toml")
    @test bytes2hex(sha256(read(raw)))=="6af7b93a8935f2e90d5f49bac160f79e4df38cd80c793723dce0574ac3465971"
    @test bytes2hex(sha256(read(projection)))=="0f7617ea22ec2afce12005c8effa006aa5b18038ef5e9c928443b36fdb8be3ff"
    corpus=TOML.parsefile(projection)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(corpus["cases"])==8
    for case in corpus["cases"]
        @testset "periodic mesher=$(case["mesh_algorithm"]) selector=$(case["recombine_algorithm"]) tilted=$(case["tilted"])" begin
        start()
        try
            points=[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)]
            case["tilted"] && (points=[(x,y,x+2y) for (x,y,_) in points])
            for index in 0:1
                offset=10index
                for (tag,(x,y,z)) in enumerate(points)
                    call(A.model.add_point,x+2index,y,z;tag=tag+offset,meshSize=10.)
                end
                for tag in 1:4
                    call(A.model.add_line,tag+offset,mod1(tag+1,4)+offset;tag=tag+offset)
                end
                call(A.model.add_curve_loop,collect(1:4).+offset;tag=index+1)
                call(A.model.add_plane_surface,[index+1];tag=index+1)
            end
            for offset in (0,10),tag in 1:4
                call(A.mesh.set_transfinite_curve,tag+offset,3)
            end
            call(A.mesh.set_periodic,2,[2],[1],case["periodic"]["affine"])
            for (name,value) in (("Mesh.Renumber",0),("Mesh.Algorithm",case["mesh_algorithm"]),
                    ("Mesh.RecombinationAlgorithm",case["recombine_algorithm"]),
                    ("Mesh.RecombineNodeRepositioning",0),("Mesh.RecombineOptimizeTopology",0),
                    ("Mesh.Smoothing",0),("Mesh.MeshSizeMin",10.),("Mesh.MeshSizeMax",10.))
                call(A.option,name,value)
            end
            for face in (1,2);call(A.mesh.set_recombine,2,face);end
            model=A.CURRENT[];cache=A.LAST_MESH[];before=fingerprint((model,cache))
            word=Tessella.GmshLibm._win_x87_control_word()
            call(A.mesh.generate,2)
            @test fingerprint((model,cache))==before
            @test call(A.mesh.get_elements)==(case["types"],
                [UInt64.(row) for row in case["tags"]],[UInt64.(row) for row in case["cells"]])
            @test Set(first(call(A.mesh.get_nodes)))==Set(UInt64(row["tag"]) for row in case["nodes"])
            for row in case["nodes"]
                point,parameter,dimension,owner=call(A.mesh.get_node,row["tag"])
                @test bits(point)==bits(row["xyz"])
                @test bits(parameter)==bits(row["parameters"])
                @test (dimension,owner)==Tuple(row["owner"])
            end
            @test call(A.mesh.get_max_node_tag)==case["max_node"]
            @test call(A.mesh.get_max_element_tag)==case["max_element"]
            @test call(A.logger.get_last_error)==case["last_error"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            @test call(A.option,"Mesh.MeshSizeFactor")==1
            mapping=call(A.mesh.get_periodic_nodes,2,2)
            @test mapping.master_entity==case["periodic"]["master"]
            @test Dict(zip(mapping.slave_nodes,mapping.master_nodes))==
                Dict(zip(UInt64.(case["periodic"]["slave_nodes"]),UInt64.(case["periodic"]["master_nodes"])))
            @test bits(mapping.affine)==bits(case["periodic"]["affine"])
            @test Tessella.GmshLibm._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
        end
    end
end

@testset "full-quad geometry option survives detached copy, retagging and GEO entry" begin
    start()
    try
        case=deepcopy(first(TOML.parsefile(PROJECTION)["cases"]))
        case["second_order_linear"]=true
        fixture(case)
        call(A.mesh.generate,2)
        model=A.CURRENT[]
        detached=call(Tessella.Model._model_mesh_transaction_copy,model)
        @test detached.meshing.second_order_linear
        @test model.meshing.second_order_linear
        call(A.model.set_tag,0,1,11)
        @test A.CURRENT[].meshing.second_order_linear
        @test call(A.option,"Mesh.SecondOrderLinear")==1
    finally
        call(A.finalize)
    end
    start()
    try
        mktempdir() do directory
            path=joinpath(directory,"linear_fullquad.geo")
            write(path,"Mesh.SecondOrderLinear=1; Point(1)={0,0,0,10}; " *
                "Point(2)={1,0,0,10}; Point(3)={0,1,0,10}; " *
                "Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1}; " *
                "Curve Loop(1)={1,2,3}; Plane Surface(1)={1}; " *
                "Transfinite Curve{1,2,3}=3; Recombine Surface{1}; " *
                "Mesh.RecombinationAlgorithm=3; Mesh.Algorithm=3; " *
                "Mesh.RecombineNodeRepositioning=0; Mesh.RecombineOptimizeTopology=0; Mesh.Smoothing=0; Mesh.Renumber=0;")
            call(A.open_geo!,path;mesh_dim=2)
        end
        @test A.CURRENT[].meshing.second_order_linear
        @test call(A.option,"Mesh.SecondOrderLinear")==1
        @test call(A.mesh.get_elements,2,1)[1]==[3]
        @test length(only(call(A.mesh.get_elements,2,1)[2]))==3
    finally
        call(A.finalize)
    end
    start()
    try
        @test call(A.option,"Mesh.SecondOrderLinear")==0
        @test !A.CURRENT[].meshing.second_order_linear
    finally
        call(A.finalize)
    end
end
end
end
end
