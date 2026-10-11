module APINativeSelector4PeriodicTests
using Test,Tessella,TOML,SHA,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const F=getfield(parentmodule(@__MODULE__),:APINativeFullQuadGenerationTests)
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
start()=call(F.start)
@testset "selector4 Plane mesh-periodic masters retain lower identities after consumed Face products" begin
    raw=joinpath(@__DIR__,"../artifacts/native_selector4_periodic_primary.json")
    projection=joinpath(@__DIR__,"../artifacts/native_selector4_periodic_primary.toml")
    @test bytes2hex(sha256(read(raw)))=="c801e0e36bf505fffe43a4e59f75c565a9e2c667379ab30b56da8e174962c12d"
    @test bytes2hex(sha256(read(projection)))=="1ec80ae78db0e3fc5bc7f0492b641a3f89027d5f0a0138248dcfdfe2b6623dba"
    corpus=TOML.parsefile(projection)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(corpus["cases"])==4
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


end
end
