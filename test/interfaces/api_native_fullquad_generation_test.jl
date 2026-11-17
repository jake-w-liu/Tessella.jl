module APINativeFullQuadGenerationTests
using Test,Tessella
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
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
@testset "full-quad geometry option survives detached copy, retagging and GEO entry" begin
    start()
    try
        case=Dict("shape"=>"triangle","tilted"=>false,"mesh_algorithm"=>3,
            "recombine_algorithm"=>3,"second_order_linear"=>true,
            "reposition"=>0,"topology"=>0,"smoothing"=>0)
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
end # module
