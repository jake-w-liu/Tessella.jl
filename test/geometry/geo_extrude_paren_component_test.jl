# An `Extrude {dx, dy, dz}` component starting with a parenthesized term
# followed by more expression — `(left?-1:1)*lb/2` in
# examples/simple_geo/indheat.geo, `(zmax - zmin) + 2 * l` in
# examples/boolean/hybrid_occ_builtin.geo — is a scalar FExpr, not a
# `( VExpr )` group. Coordinates printed by Gmsh 4.15.2.
using Test
using Tessella

@testset "geo Extrude parenthesized scalar components" begin
    dir=mktempdir()
    path=joinpath(dir,"paren.geo")
    write(path,"""
        left = 1; lb = 0.4; zmax = 2; zmin = 0.5; l = 0.25;
        Point(1) = {0,0,0,0.1}; Point(2) = {1,0,0,0.1};
        a[] = Extrude {(left?-1:1)*lb/2, 0, 0} { Point{1}; };
        b[] = Extrude {0, 0, (zmax - zmin) + 2 * l} { Point{2}; };
        c[] = Extrude {(1), (2)*0.5, 0} { Point{2}; };
        pa[] = Point{a[0]}; pb[] = Point{b[0]}; pc[] = Point{c[0]};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.lists["pa"]≈[-0.2,0.0,0.0]
    @test ex.lists["pb"]≈[1.0,0.0,2.0]
    @test ex.lists["pc"]≈[2.0,1.0,0.0]
end
