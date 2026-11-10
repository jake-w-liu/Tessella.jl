module APINativeRecombineDiamondTests
using Test, Tessella, Serialization, SHA
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
fingerprint(value)=begin stream=IOBuffer();serialize(stream,value);sha256(take!(stream));end
const OUTER=[(-1.,0.,0.),(-.7,-.7,0.),(0.,-1.,0.),(.7,-.7,0.),
             (1.,0.,0.),(.7,.7,0.),(0.,1.,0.),(-.7,.7,0.)]
const INNER=[(-.3,0.,0.),(0.,-.3,0.),(.3,0.,0.),(0.,.3,0.)]
const CELLS=[9 9 10 10 11 12 12; 10 1 2 4 5 6 8; 11 2 3 5 6 7 1; 12 10 4 11 12 8 9]

function setup(;inner=INNER,rotation=0,parameters=false,topology=5,offset=0.)
    call(A.initialize)
    call(A.option,"General.Terminal",0)
    call(A.option,"Mesh.Renumber",0)
    call(A.option,"Mesh.Smoothing",0)
    call(A.option,"Mesh.RecombineOptimizeTopology",topology)
    points=[(point[1]+offset,point[2]+offset,point[3]) for point in vcat(OUTER,inner)]
    for i in 1:8;call(A.model.add_point,points[i]...;tag=i);end
    for i in 1:8;call(A.model.add_line,i,mod1(i+1,8);tag=i);end
    call(A.model.add_curve_loop,collect(1:8);tag=1)
    call(A.model.add_plane_surface,[1];tag=1)
    for i in 1:8;call(A.mesh.add_nodes,0,i,[i],collect(points[i]));end
    xyz=collect(Iterators.flatten(points[9:12]))
    uv=parameters ? call(A.model.get_parametrization,2,1,xyz) : Float64[]
    call(A.mesh.add_nodes,2,1,[9,10,11,12],xyz,uv)
    cells=copy(CELLS);cells[:,1]=circshift(cells[:,1],rotation)
    call(A.mesh.add_elements_by_type,1,3,collect(1:7),vec(cells))
    return cells
end

# Independently captured with the public Gmsh 4.15.2 API, Renumber=0,
# Smoothing=0, RecombineOptimizeTopology=5. No external runtime is needed.
@testset "Native Plane diamond contraction and original cell identities" begin
    for rotation in (0,2),parameters in (false,true),offset in (0.,1e9)
        setup(;rotation,parameters,offset)
        try
            keep,drop=rotation==0 ? (9,11) : (11,9)
            before=[call(A.mesh.get_node,i) for i in 1:12]
            held=(A.CURRENT[],A.LAST_MESH[]);digest=fingerprint(held)
            call(A.mesh.recombine)
            @test fingerprint(held)==digest
            types,tags,nodes=call(A.mesh.get_elements,2,1)
            expected=copy(CELLS[:,2:end]);expected[expected.==drop].=keep
            @test types==[3]
            @test only(tags)==UInt64.(2:7)
            @test only(nodes)==UInt64.(vec(expected))
            @test sort(first(call(A.mesh.get_nodes)))==UInt64[i for i in 1:12 if i!=drop]
            @test call(A.mesh.get_node,keep)[1]==[offset,offset,0.]
            @test_throws ArgumentError call(A.mesh.get_node,drop)
            @test [call(A.mesh.get_node,i) for i in 1:8]==before[1:8]
            @test call(A.mesh.get_node,10)==before[10]
            @test call(A.mesh.get_node,12)==before[12]
            uv=call(A.mesh.get_node,keep)[2]
            if parameters
                @test uv==call(A.model.get_parametrization,2,1,[offset,offset,0.])
            else
                @test isempty(uv)
            end
            @test call(A.mesh.get_max_node_tag)==12
            @test call(A.mesh.get_max_element_tag)==7
            snapshot=(call(A.mesh.get_elements,2,1),call(A.mesh.get_nodes))
            call(A.mesh.recombine)
            @test (call(A.mesh.get_elements,2,1),call(A.mesh.get_nodes))==snapshot
        finally
            call(A.finalize)
        end
    end
end

@testset "A nonconvex trial and disabled topology preserve the input" begin
    rejected=[(-.68,-.14,0.),(-.21,-.56,0.),(.57,-.12,0.),(-.26,.08,0.)]
    for (inner,topology) in ((rejected,5),(INNER,0))
        setup(;inner,topology,parameters=true)
        try
            before=(call(A.mesh.get_elements,2,1),call(A.mesh.get_nodes),
                [call(A.mesh.get_node,i) for i in 1:12])
            call(A.mesh.recombine)
            @test (call(A.mesh.get_elements,2,1),call(A.mesh.get_nodes),
                [call(A.mesh.get_node,i) for i in 1:12])==before
        finally
            call(A.finalize)
        end
    end
end

@testset "Four degree-three corners do not form a contractible diamond" begin
    call(A.initialize)
    try
        call(A.option,"General.Terminal",0);call(A.option,"Mesh.Renumber",0)
        call(A.option,"Mesh.Smoothing",0)
        points=[(-1.,0.,0.),(0.,-1.,0.),(1.,0.,0.),(0.,1.,0.),INNER...]
        for i in 1:4;call(A.model.add_point,points[i]...;tag=i);end
        for i in 1:4;call(A.model.add_line,i,mod1(i+1,4);tag=i);end
        call(A.model.add_curve_loop,[1,2,3,4];tag=1)
        call(A.model.add_plane_surface,[1];tag=1)
        for i in 1:4;call(A.mesh.add_nodes,0,i,[i],collect(points[i]));end
        call(A.mesh.add_nodes,2,1,[5,6,7,8],collect(Iterators.flatten(INNER)))
        cells=[5,6,7,8,5,1,2,6,5,8,4,1,7,6,2,3,7,3,4,8]
        call(A.mesh.add_elements_by_type,1,3,collect(1:5),cells)
        before=(call(A.mesh.get_elements,2,1),call(A.mesh.get_nodes))
        call(A.mesh.recombine)
        @test (call(A.mesh.get_elements,2,1),call(A.mesh.get_nodes))==before
    finally
        call(A.finalize)
    end
end

@testset "Newly recombined triangle pair participates in diamond contraction" begin
    setup()
    try
        call(A.mesh.remove_elements,2,1,[1])
        call(A.mesh.add_elements_by_type,1,2,[1,8],[9,10,11,9,11,12])
        call(A.mesh.recombine)
        types,tags,nodes=call(A.mesh.get_elements,2,1)
        expected=copy(CELLS[:,2:end]);expected[expected.==9].=11
        @test types==[3]
        @test only(tags)==UInt64.(2:7)
        @test only(nodes)==UInt64.(vec(expected))
        @test sort(first(call(A.mesh.get_nodes)))==UInt64[i for i in 1:12 if i!=9]
        @test call(A.mesh.get_node,11)[1]==[0.,0.,0.]
        @test call(A.mesh.get_max_element_tag)==10
    finally
        call(A.finalize)
    end
end

@testset "Topology-only diamonds retain active boundary-layer field contexts" begin
    for parameters in (false,true)
        setup(;parameters)
        try
            field=call(A.mesh.field.add,"BoundaryLayer")
            call(A.mesh.field.set_as_boundary_layer,field)
            call(A.mesh.recombine)
            types,tags,nodes=call(A.mesh.get_elements,2,1)
            expected=copy(CELLS[:,2:end]);expected[expected.==11].=9
            @test types==[3]
            @test only(tags)==UInt64.(2:7)
            @test only(nodes)==UInt64.(vec(expected))
            @test call(A.mesh.get_node,9)[1]==[0.,0.,0.]
        finally
            call(A.finalize)
        end
    end
end
end
