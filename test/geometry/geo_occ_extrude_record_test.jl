# Extruding a surface that carries an OpenCASCADE geometry record (a `Disk`)
# moves the record with the copied top face; the translation transform used
# to be built with a missing argument and threw a `MethodError`.
using Test
using Tessella

@testset "geo Extrude of OCC-record surfaces" begin
    dir=mktempdir()
    path=joinpath(dir,"disk_extrude.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Disk(1) = {0,0,0,1,0.5};
        a() = Extrude{0,0,0.3}{ Surface{1}; };
        bb() = BoundingBox Surface{a(0)};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test length(ex.lists["a"])==3
    @test ex.lists["a"][2]==1.0
    # Gmsh 4.15.2 prints `-1 -0.5 0.3 1 0.5 0.3` for the top face's box.
    @test ex.lists["bb"]≈[-1.0,-0.5,0.3,1.0,0.5,0.3] atol=1e-9
end
