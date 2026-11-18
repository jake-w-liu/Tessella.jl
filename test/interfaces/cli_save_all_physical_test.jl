# CLI MSH output honours `Mesh.SaveAll = 0`: with physical groups present only
# elements of physical entities are written (Gmsh 4.15.2 `tutorials/t1.geo`
# -2: 70 lines on the three physical curves + 726 triangles, no point
# elements); `-save_all` writes everything. Here every node is referenced
# by a physical element, so the node count is unchanged.
using Test
using Tessella

function _msh_element_header(path)
    lines=readlines(path)
    i=findfirst(==("\$Elements"),lines)
    return parse.(Int,split(lines[i+1]))
end

@testset "CLI Mesh.SaveAll physical selection" begin
    dir=mktempdir()
    geo=joinpath(dir,"square.geo")
    write(geo,"""
        Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
        Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
        Line(1) = {1,2}; Line(2) = {3,2}; Line(3) = {3,4}; Line(4) = {4,1};
        Curve Loop(1) = {4,1,-2,3};
        Plane Surface(1) = {1};
        Physical Curve(5) = {1,2,4};
        Physical Surface("My surface") = {1};
        """)
    out=joinpath(dir,"phys.msh")
    Tessella.CLI.main([geo,"-2","-o",out])
    all_out=joinpath(dir,"all.msh")
    Tessella.CLI.main([geo,"-2","-save_all","-o",all_out])
    phys=read_mixed_msh(out)
    full=read_mixed_msh(all_out)
    dims(m)=[Tessella.Elements.msh_dimension(b.msh) for b in m.blocks]
    count(m,d)=sum(length(b.tags) for b in m.blocks
                   if Tessella.Elements.msh_dimension(b.msh)==d;init=0)
    @test count(phys,0)==0
    @test count(full,0)==4
    @test count(full,1)-count(phys,1)==count(full,1)÷4
    @test count(phys,2)==count(full,2)>0
    @test size(phys.coords,2)==size(full.coords,2)
    @test _msh_element_header(out)[1]==4

    # Nodes of a mesh-carrying entity outside every physical group (and not
    # referenced by a written element) are dropped, as in Gmsh's output for
    # tutorials/t2.geo -2.
    two=joinpath(dir,"two.geo")
    write(two,"""
        Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
        Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
        Point(5) = {3,0,0,0.5}; Point(6) = {4,0,0,0.5};
        Point(7) = {4,1,0,0.5}; Point(8) = {3,1,0,0.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
        Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
        Curve Loop(2) = {5,6,7,8}; Plane Surface(2) = {2};
        Physical Surface(1) = {1};
        """)
    two_out=joinpath(dir,"two.msh")
    Tessella.CLI.main([two,"-2","-o",two_out])
    two_all=joinpath(dir,"two_all.msh")
    Tessella.CLI.main([two,"-2","-save_all","-o",two_all])
    sel=read_mixed_msh(two_out);full2=read_mixed_msh(two_all)
    @test 2*size(sel.coords,2)==size(full2.coords,2)
    @test all(x->x<=1.0+1e-12,sel.coords[1,:])
end
