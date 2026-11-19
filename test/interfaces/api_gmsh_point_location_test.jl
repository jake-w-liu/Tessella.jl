# Point location on meshes loaded from MSH files answers what Gmsh 4.15.2 answers.
# The test/fixtures/gmsh_4152_api_*.msh files were written by Gmsh 4.15.2 and
# every expected value below is what `gmsh` returned for the same sequence of
# calls after reloading the file (the search remembers its last hit).
module GmshPointLocationParityChecks
using Test, Tessella
const A=Tessella.API
const DATA=joinpath(@__DIR__,"..","fixtures")

function load(name)
    A.finalize();A.initialize();A.option("General.Terminal",0)
    A.open(joinpath(DATA,"gmsh_4152_api_"*name*".msh"))
    return nothing
end
find(x,y,z,dim=-1,strict=false)=A.mesh.get_element_by_coordinates(x,y,z,dim,strict)
find_all(x,y,z,dim=-1,strict=false)=Int.(A.mesh.get_elements_by_coordinates(x,y,z,dim,strict))
locate(tag,x,y,z)=collect(A.mesh.get_local_coordinates_in_element(tag,x,y,z))
function same(found,tag,type,nodes,uvw)
    return found[1]==tag && found[2]==type && found[3]==nodes &&
        isapprox(collect(found[4:6]),uvw;atol=1e-12)
end

@testset "elements at a point, in MElementOctree's order" begin
    try
        load("square_o1")
        @test same(find(0.25,0.5,0.0),9,2,[1,2,4],[0.25,0.5,0.0])
        @test same(find(0.5,0.0,0.0,1),5,1,[1,2],[0.0,0.0,0.0])
        @test same(find(0.0,0.0,0.0,0),1,15,[1],[0.0,0.0,0.0])
        # A corner belongs to a triangle, two lines and a point element; the
        # list is not sorted by dimension or tag.
        @test find_all(0.0,0.0,0.0)==[9,8,5,1]
        @test find_all(0.5,0.5,0.0)==[9,10]
        # MTriangle::xyz2uvw drops the point onto the triangle along an axis
        # and reports w = 0, so a relaxed search finds the triangle below a
        # point off the mesh; a strict one does not look that far.
        @test same(find(0.25,0.5,0.3,-1,false),9,2,[1,2,4],[0.25,0.5,0.0])
        @test_throws ArgumentError find(0.25,0.5,0.3,-1,true)
        @test_throws ArgumentError find_all(0.25,0.5,0.3,-1,true)
        # Outside the mesh the relaxed search widens the reference tolerance
        # up to one and returns the first element that then contains the point.
        @test same(find(2.0,0.5,0.0,-1,false),10,2,[4,2,3],[0.5,1.5,0.0])
        @test find_all(2.0,0.5,0.0,-1,false)==[10]
        @test_throws ArgumentError find(2.0,0.5,0.0,-1,true)
        load("cube_o1")
        @test same(find(0.3,0.4,0.2),27,5,[3,1,2,4,7,5,6,8],[0.2,0.6,-0.4])
        @test same(find(0.0,0.4,0.2,2),21,3,[2,1,3,4],[-0.6,-0.2,0.0])
        @test find_all(1.0,1.0,1.0)==[7,14,15,20,22,24,26,27]
        @test same(find(0.3,0.4,1.5,-1,false),27,5,[3,1,2,4,7,5,6,8],[0.2,-2.0,-0.4])
        @test_throws ArgumentError find(0.3,0.4,1.5,-1,true)
        load("tets")
        @test same(find(0.5,0.5,0.5),34,4,[4,1,6,5],[0.0,0.0,0.5])
        @test find_all(0.5,0.5,0.5)==[34,35,36,37]
        @test same(find(3.0,0.4,0.2,-1,false),21,2,[2,1,4],[0.2,0.4,0.0])
        @test find_all(3.0,0.4,0.2,-1,false)==[21,23]
        load("pyramid")
        @test same(find(0.4,0.45,0.3),1,7,[2,3,4,1,5],[-0.1,0.2,0.375])
        load("square_o3")
        @test same(find(0.3,0.4,0.0),9,21,[1,2,4,5,6,13,14,11,12,15],[0.3,0.4,0.0])
    finally
        A.finalize()
    end
end

@testset "local coordinates are MElement::xyz2uvw" begin
    try
        load("square_o1")
        # The Newton loop of a line, point, quadrangle or curved element stops
        # on its step, so the offset from the element is added nineteen times
        # to the coordinates the element does not span.
        @test locate(5,0.25,0.5,0.3)≈[-0.5,-9.5,-5.7] atol=1e-12
        @test locate(1,0.25,0.5,0.3)≈[4.75,9.5,5.7] atol=1e-12
        @test locate(9,0.25,0.5,0.3)==[0.25,0.5,0.0]
        load("cube_o1")
        @test locate(21,0.3,0.4,0.2)≈[-0.6,-0.2,-5.7] atol=1e-12
        load("square_o3")
        @test locate(9,0.3,0.4,0.2)≈[0.3,0.4,3.8] atol=1e-12
        load("pyramid")
        @test locate(1,0.9,0.9,0.5)≈[0.8,-0.8,0.625] atol=1e-12
    finally
        A.finalize()
    end
end

end # module
