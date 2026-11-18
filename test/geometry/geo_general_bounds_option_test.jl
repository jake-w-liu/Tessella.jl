# `General.MinX`…`General.MaxZ` read the synchronized model's bounding box
# (`GModel::current()->bounds()`, 0 for an empty model); they used to stay at
# their static 0 default, so examples/boolean/hybrid_occ_builtin.geo built a
# degenerate enclosing box. Values printed by Gmsh 4.15.2 on the same script.
using Test
using Tessella

@testset "geo General.MinX..MaxZ options" begin
    dir=mktempdir()
    path=joinpath(dir,"bounds.geo")
    write(path,"""
        a() = {General.MinX, General.MaxX, General.MinY, General.MaxY, General.MinZ, General.MaxZ};
        Point(1) = {1,2,3,0.1}; Point(2) = {4,6,3,0.1};
        b() = {General.MinX, General.MaxX, General.MinY, General.MaxY, General.MinZ, General.MaxZ};
        SetFactory("Built-in");
        c() = {General.MinX, General.MaxX, General.MinY, General.MaxY, General.MinZ, General.MaxZ};
        Point(3) = {-1,0,5,0.1};
        BoundingBox {0,10,0,20,0,0};
        d() = {General.MinX, General.MaxX, General.MinY, General.MaxY, General.MinZ, General.MaxZ};
        BoundingBox;
        e() = {General.MinX, General.MaxX, General.MinY, General.MaxY, General.MinZ, General.MaxZ};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.lists["a"]==zeros(6)
    # Not yet synchronized: still the empty-model zeros.
    @test ex.lists["b"]==zeros(6)
    @test ex.lists["c"]==[1.0,4.0,2.0,6.0,3.0,3.0]
    # A forced scene box does not change the model bounds.
    @test ex.lists["d"]==[1.0,4.0,2.0,6.0,3.0,3.0]
    @test ex.lists["e"]==[-1.0,4.0,0.0,6.0,3.0,5.0]
end
