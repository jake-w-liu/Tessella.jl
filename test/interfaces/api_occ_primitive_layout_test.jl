module OCCPrimitiveLayoutChecks
using Test,Tessella
const A=Tessella.API

# Scripts address the faces and edges of a primitive by tag, so the tags have
# to land where OpenCASCADE puts them. `occ.add_box` numbered its points from
# (xmin,ymin,zmin) around the bottom and made Surface 1 the zmin face; Gmsh
# 4.15.2 starts at (xmin,ymin,zmax), Surface 1 is the xmin face, and the
# three minimum faces are used reversed (their natural normal is the positive
# axis). `occ.add_wedge` is the same layout with the xmax edge at ymax
# collapsed, and an apex cone's base is Surface 2, not 3. Every list below is
# Gmsh 4.15.2's.
boundary(d,t)=[s for (_,s) in A.model.get_boundary([(d,t)],false,true,false)]
point(t)=Tuple(round.(A.model.get_bounding_box(0,t)[1:3];digits=12))

@testset "OCC primitives have OpenCASCADE's entity layout" begin
    A.initialize()
    try
        occ=A.model.occ
        box=occ.add_box(1,2,3,0.5,0.25,2)
        occ.synchronize()
        @test [point(t) for t in 1:8]==[(1.0,2.0,5.0),(1.0,2.0,3.0),
            (1.0,2.25,5.0),(1.0,2.25,3.0),(1.5,2.0,5.0),(1.5,2.0,3.0),
            (1.5,2.25,5.0),(1.5,2.25,3.0)]
        @test [boundary(1,t) for t in 1:12]==[[2,1],[1,3],[4,3],[2,4],
            [6,5],[5,7],[8,7],[6,8],[2,6],[1,5],[4,8],[3,7]]
        # Gmsh's wires as it stores them (`getBoundary(combined=false)`)
        @test [boundary(2,t) for t in 1:6]==[[-1,4,3,-2],[-5,8,7,-6],
            [-9,1,10,-5],[-11,3,12,-7],[-4,9,8,-11],[-2,10,6,-12]]
        @test boundary(3,box)==[-1,2,-3,4,-5,6]
        @test occ.get_mass(3,box)≈0.25 rtol=1e-14
        @test collect(occ.get_center_of_mass(3,box))≈[1.25,2.125,4.0] rtol=1e-14
    finally
        A.finalize()
    end

    A.initialize()
    try
        occ=A.model.occ
        wedge=occ.add_wedge(0,0,0,1,2,3)
        occ.synchronize()
        @test [point(t) for t in 1:6]==[(0.0,0.0,3.0),(0.0,0.0,0.0),
            (0.0,2.0,3.0),(0.0,2.0,0.0),(1.0,0.0,3.0),(1.0,0.0,0.0)]
        @test [boundary(1,t) for t in 1:9]==[[2,1],[1,3],[4,3],[2,4],
            [6,5],[5,3],[6,4],[2,6],[1,5]]
        @test [boundary(2,t) for t in 1:5]==[[-1,4,3,-2],[-5,7,3,-6],
            [-8,1,9,-5],[-4,8,7],[-2,9,6]]
        @test boundary(3,wedge)==[-1,2,-3,-4,5]
        @test occ.get_mass(3,wedge)≈3.0 rtol=1e-14
    finally
        A.finalize()
    end

    A.initialize()
    try
        occ=A.model.occ
        cone=occ.add_cone(0,0,0,0,0,2,1,0)
        occ.synchronize()
        @test [t for (_,t) in A.model.get_entities(2)]==[1,2]
        @test boundary(3,cone)==[1,-2]
    finally
        A.finalize()
    end
end
end
