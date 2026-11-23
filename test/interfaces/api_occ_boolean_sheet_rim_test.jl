module OCCBooleanSheetRimChecks
using Test,Tessella
const A=Tessella.API

# Sheets that touch only along an edge share that edge after a fragment, as
# OpenCASCADE's General Fuse glues them. Tags, edge ends and the faces'
# oriented boundaries in stored order (`combined=false`) are Gmsh 4.15.2's:
# a rebuilt face keeps its source's edge order, the shared rim where the
# source had its own edge.
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
        @test boundary(2)==[5,6,7,-2]
        @test ends()==Dict(1=>(1,2),2=>(2,3),3=>(3,4),4=>(4,1),
                           5=>(2,5),6=>(5,6),7=>(6,3))
        @test [point(t) for t in 1:6]==[(0.0,0.0,0.0),(1.0,0.0,0.0),
            (1.0,1.0,0.0),(0.0,1.0,0.0),(2.0,0.0,0.0),(2.0,1.0,0.0)]
        @test occ.get_mass(2,1)≈1 && occ.get_mass(2,2)≈1
    finally
        A.finalize()
    end
    # operands swapped: the tool face keeps its source order with the
    # object's edge 4 where its own edge 2 was. Gmsh 4.15.2 lists 5,-4,6,7
    # (`combined=false`), a connected chain; the combined listing reorders.
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_rectangle(0,0,0,1,1;tag=1)
        occ.add_rectangle(1,0,0,1,1;tag=2)
        out,images=occ.fragment([(2,2)],[(2,1)])
        occ.synchronize()
        @test out==[(2,2),(2,1)]
        @test boundary(1)==[5,-4,6,7]
        @test boundary(2)==[1,2,3,4]
        @test ends()==Dict(1=>(1,2),2=>(2,3),3=>(3,4),4=>(4,1),
                           5=>(5,1),6=>(4,6),7=>(6,5))
        @test [point(t) for t in 5:6]==[(0.0,0.0,0.0),(0.0,1.0,0.0)]
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
        @test boundary(2)==[6,7,8,9,-3]
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
        @test boundary(2)==[5,6,7,-2]
        @test boundary(3)==[-3,8,9,10]
        @test ends()[8]==(3,7) && ends()[10]==(8,4)
    finally
        A.finalize()
    end
end

# extrude.geo's fragment: a rectangle, a square inside it and a disk over
# both. A split sheet lists its pieces as its own edges meet them, then as
# the other operands' edges meet them, last operand first, the piece inside
# that operand first; a coincident piece stands where its own region was.
# Gmsh 4.15.2, by piece centroid.
@testset "OCC fragment order of split sheet pieces" begin
    r2(x)=round(x;digits=2)
    for (ord,maps,cs) in (
            ((1,2,3),[[1,2,3,4],[4,3],[5,2,3]],
             [(0.45,0.43),(0.6,0.62),(0.37,0.42),(0.33,0.28),(1.04,0.6)]),
            ((1,3,2),[[1,2,3,4],[5,2,4],[3,4]],
             [(0.45,0.43),(0.6,0.62),(0.33,0.28),(0.37,0.42),(1.04,0.6)]),
            ((2,1,3),[[1,2],[3,4,2,1],[5,4,2]],
             [(0.33,0.28),(0.37,0.42),(0.45,0.43),(0.6,0.62),(1.04,0.6)]))
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            occ.add_rectangle(0,0,0,1,1;tag=1)
            occ.add_rectangle(0.2,0.2,0,0.3,0.3;tag=2)
            occ.add_disk(0.6,0.6,0,0.5,0.3;tag=3)
            occ.synchronize()
            out,images=occ.fragment([(2,ord[1])],[(2,t) for t in ord[2:end]])
            occ.synchronize()
            @test [[t for (_,t) in m] for m in images]==maps
            @test [Tuple(r2.(occ.get_center_of_mass(2,t)[1:2])) for t in 1:5]==cs
        finally
            A.finalize()
        end
    end
end
end
