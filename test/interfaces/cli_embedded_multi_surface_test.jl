# `-2` on a model with an embedded entity and several surfaces used to abort
# with "classified surface projection requires exactly one selected surface"
# (tutorials/t15.geo). Gmsh 4.15.2 writes 22 nodes, 1 point and 30 triangles
# for this script.
using Test
using Tessella

@testset "CLI embedded entities with several surfaces" begin
    dir=mktempdir()
    geo=joinpath(dir,"emb2.geo")
    write(geo,"""
        Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
        Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
        Point(5) = {2,0,0,0.5}; Point(6) = {2,1,0,0.5};
        Point(7) = {0.5,0.5,0,0.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {2,5}; Line(6) = {5,6}; Line(7) = {6,3};
        Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
        Curve Loop(2) = {5,6,7,-2}; Plane Surface(2) = {2};
        Point{7} In Surface{1};
        Physical Point("p") = {7};
        Physical Surface("s") = {1,2};
        """)
    out=joinpath(dir,"emb2.msh")
    @test Tessella.CLI.main([geo,"-2","-o",out])==out
    m=read_mixed_msh(out)
    count(d)=sum(length(b.tags) for b in m.blocks
                 if Tessella.Elements.msh_dimension(b.msh)==d;init=0)
    @test size(m.coords,2)==22
    @test count(0)==1
    @test count(2)==30
    @test any(i->m.coords[:,i]≈[0.5,0.5,0.0],axes(m.coords,2))
end
