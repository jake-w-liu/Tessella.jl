module OCCRevolveSplineChecks
using Test,Tessella
const A=Tessella.API

# `occ.revolve` of a spline threw "revolution of a tensor_bspline profile is
# not yet supported — only line curves are". The swept face is now the exact
# tensor product of the rational circle and the profile. Tags, boundary
# order, types of the curves and the areas are Gmsh 4.15.2's
# (`BRepPrimAPI_MakeRevol`).
function spline_profile(occ,points)
    occ.add_spline([occ.add_point(x,0,z) for (x,z) in points])
end

@testset "OCC revolve of a spline: full turn" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        curve=spline_profile(occ,((1,0),(1.5,0.5),(1.2,1.0)))
        @test occ.revolve([(1,curve)],0,0,0,0,0,1,2π)==[(1,1),(2,1),(1,2),(1,3)]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(3,3,1,0)
        @test [A.model.get_type(1,t) for t in 1:3]==["BSpline","Circle","Circle"]
        @test A.model.get_boundary([(2,1)],false,true)==[(1,2),(1,1),(1,-3),(1,-1)]
        @test occ.get_mass(1,2)≈2π rtol=1e-12
        @test occ.get_mass(1,3)≈2.4π rtol=1e-12
        # u = 0 is the profile itself, and the normal points off the axis
        on_profile=A.model.get_value(1,1,[0.5])
        @test collect(A.model.get_value(2,1,[0.0,0.5]))≈collect(on_profile) atol=1e-12
        normal=A.model.get_normal(1,[0.0,0.5])
        @test normal[1]>0.8 && abs(normal[2])<1e-12
        if Base.JLOptions().compile_enabled==1
            # Gmsh: 11.2982537820
            @test occ.get_mass(2,1)≈11.2982537820 atol=2e-9
            A.option("Mesh.MeshSizeMax",0.2)
            A.mesh.generate(2)
            # Gmsh: 393 nodes, 714 triangles
            @test 350<length(A.mesh.get_nodes()[1])<440
        end
    finally
        A.finalize()
    end
end

@testset "OCC revolve of a spline: quarter turn ending on the axis" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        curve=spline_profile(occ,((1,0),(0.9,0.6),(0,1.0)))
        # the profile's end on the axis sweeps a collapsed parallel, which
        # takes its tag before the rotated copy of the profile
        @test occ.revolve([(1,curve)],0,0,0,0,0,1,π/2)==[(1,4),(2,1),(1,2),(1,3)]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(4,4,1,0)
        @test [A.model.get_type(1,t) for t in 1:4]==
            ["BSpline","Circle","Unknown","BSpline"]
        @test A.model.get_boundary([(2,1)],false,true)==[(1,2),(1,4),(1,-3),(1,-1)]
        @test [A.model.get_boundary([(1,t)],false,true) for t in 1:4]==
            [[(0,1),(0,3)],[(0,1),(0,4)],[(0,3),(0,3)],[(0,4),(0,3)]]
        @test collect(A.model.get_value(0,4,Float64[]))≈[0,1,0] atol=1e-15
        if Base.JLOptions().compile_enabled==1
            # Gmsh: 1.741868289430764
            @test occ.get_mass(2,1)≈1.741868289430764 rtol=1e-10
        end
        @test_throws ArgumentError occ.revolve(
            [(1,spline_profile(occ,((1,0),(-1,0.5),(1,1.0))))],0,0,0,0,0,1,π)
    finally
        A.finalize()
    end
end
end
