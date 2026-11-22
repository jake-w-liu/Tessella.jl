module OCCBooleanSheetRimChecks
using Test,Tessella
const A=Tessella.API

# Sheets that touch only along an edge share that edge after a fragment, as
# OpenCASCADE's General Fuse glues them. Tags, edge ends and the faces'
# oriented boundaries are Gmsh 4.15.2's.
function boundary(t)
    [b[2] for b in A.model.get_boundary([(2,t)],false,true)]
end
ends()=Dict(t=>Tuple(A.model.get_adjacencies(1,t)[2]) for (_,t) in A.model.get_entities(1))
point(t)=Tuple(A.model.get_value(0,t,Float64[]))

@testset "OCC fragment of abutting rectangles" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_rectangle(0,0,0,1,1;tag=1)
        occ.add_rectangle(1,0,0,1,1;tag=2)
        out,images=occ.fragment([(2,1)],[(2,2)])
        occ.synchronize()
        @test out==[(2,1),(2,2)]
        @test images==[[(2,1)],[(2,2)]]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:2)==(6,7,2)
        @test boundary(1)==[1,2,3,4]
        @test boundary(2)==[-2,5,6,7]
        @test ends()==Dict(1=>(1,2),2=>(2,3),3=>(3,4),4=>(4,1),
                           5=>(2,5),6=>(5,6),7=>(6,3))
        @test [point(t) for t in 1:6]==[(0.0,0.0,0.0),(1.0,0.0,0.0),
            (1.0,1.0,0.0),(0.0,1.0,0.0),(2.0,0.0,0.0),(2.0,1.0,0.0)]
        @test occ.get_mass(2,1)≈1 && occ.get_mass(2,2)≈1
    finally
        A.finalize()
    end
    # a rim shared over half its length: both edges split at the
    # partner's corner and the common half is one curve
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_rectangle(0,0,0,1,1;tag=1)
        occ.add_rectangle(1,0.5,0,1,1;tag=2)
        out,images=occ.fragment([(2,1)],[(2,2)])
        occ.synchronize()
        @test images==[[(2,1)],[(2,2)]]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:2)==(8,9,2)
        @test boundary(1)==[1,2,3,4,5]
        @test boundary(2)==[-3,6,7,8,9]
        @test ends()[3]==(3,4) && ends()[9]==(8,4)
        @test point(3)==(1.0,0.5,0.0)
    finally
        A.finalize()
    end
    # three sheets around a corner
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_rectangle(0,0,0,1,1;tag=1)
        occ.add_rectangle(1,0,0,1,1;tag=2)
        occ.add_rectangle(0,1,0,1,1;tag=3)
        out,images=occ.fragment([(2,1),(2,2)],[(2,3)])
        occ.synchronize()
        @test images==[[(2,1)],[(2,2)],[(2,3)]]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:2)==(8,10,3)
        @test boundary(2)==[-2,5,6,7]
        @test boundary(3)==[-3,8,9,10]
        @test ends()[8]==(3,7) && ends()[10]==(8,4)
    finally
        A.finalize()
    end
end
end
