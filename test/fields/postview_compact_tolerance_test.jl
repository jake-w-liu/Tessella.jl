module PostViewCompactToleranceChecks
using Test,Tessella
const F=Tessella.SizeField

function compact_field(msh;reference_tolerance=1e-6)
    keyword=Dict(1=>:lines,2=>:triangles,3=>:quadrangles,4=>:tetrahedra,
                 5=>:hexahedra,6=>:prisms,7=>:pyramids)[msh]
    coordinates=Tessella.Elements.lagrange_nodes(msh)
    connectivity=reshape(Int32.(1:size(coordinates,2)),:,1)
    return F.PostViewField(coordinates,fill(2.,size(coordinates,2));
        crop_negative=false,use_closest=false,reference_tolerance,
        (keyword=>connectivity,)...)
end

@testset "Compact PostView boxes include the accepted reference tolerance" begin
    for (msh,family) in ((1,"Line"),(2,"Triangle"),(4,"Tetrahedron"),
                         (3,"Quadrangle"),(5,"Hexahedron"),(6,"Prism"),(7,"Pyramid"))
        @testset "$family" begin
            field=compact_field(msh;reference_tolerance=.25)
            query=msh in (2,4,6) ? (1.1,0.,0.) : (1.2,0.,0.)
            @test F.field_value(field,query...)≈2.
            @test F.field_value(field,1.3,0.,0.)==F.GMSH_MAX_SIZE
            ordinary=compact_field(msh)
            @test F.field_value(ordinary,1+5e-7,0.,0.)≈2.
            @test F.field_value(ordinary,1+2e-6,0.,0.)==F.GMSH_MAX_SIZE
        end
    end
    field=F.PostViewField([0. 1.;0. 0.;0. 0.],[1.,2.];lines=reshape(Int32[1,2],2,1),
        reference_tolerance=.25,crop_negative=false,use_closest=false)
    @test F.field_value(field,1.1,0.,0.)≈2.1
end

@testset "Compact artificial normals retain their physical scale" begin
    # The Line kernel builds normals from the full span, while a scheme
    # Jacobian uses its half-span. Its second normal also carries the
    # physical cell scale. This small real Line catches both differences.
    line=F.PostViewField([0. 10.;0. 0.;0. 0.],[2.,2.];lines=reshape(Int32[1,2],2,1),
        reference_tolerance=.25,use_closest=false)
    @test F.field_value(line,5.,1.,0.)≈2.
    @test F.field_value(line,5.,3.,0.)==F.GMSH_MAX_SIZE
    @test F.field_value(line,5.,0.,10.)≈2.
    @test F.field_value(line,5.,0.,30.)==F.GMSH_MAX_SIZE
    quad=compact_field(3;reference_tolerance=.25)
    @test F.field_value(quad,0.,0.,.05)≈2.
    @test F.field_value(quad,0.,0.,.2)==F.GMSH_MAX_SIZE
    # Compact triangles deliberately use the two-axis projection and their
    # original octree thickness. Growing an intrinsic box must not thicken
    # that distinct normal gate as a side effect.
    triangle=compact_field(2;reference_tolerance=.25)
    @test F.field_value(triangle,.2,.2,.005)≈2.
    @test F.field_value(triangle,.2,.2,.02)==F.GMSH_MAX_SIZE
end

@testset "Tilted triangle tolerance follows the lifted plane" begin
    coordinates=[0. 1. 0.;0. 0. 1.;0. 1. 1.]
    cell=reshape(Int32[1,2,3],3,1)
    for order in ((1,2,3),(2,3,1),(3,1,2))
        field=F.PostViewField(coordinates[:,collect(order)],fill(2.,3);
            triangles=cell,reference_tolerance=.25,use_closest=false)
        @test F.field_value(field,1.1,0.,1.1)≈2.
        @test F.field_value(field,.2,.2,.4)≈2.
        @test F.field_value(field,.2,.2,.405)≈2.
        @test F.field_value(field,.2,.2,1.)==F.GMSH_MAX_SIZE
        @test F.field_value(field,1.3,0.,1.3)==F.GMSH_MAX_SIZE
    end
    # The diagonal overflows, but its one-percent thickness is finite.
    extreme=F.PostViewField([-1e308 1e308 0.;0. 0. 1e308;0. 0. 0.],fill(2.,3);
        triangles=cell,use_closest=false)
    @test F.field_value(extreme,0.,1e307,1e306)≈2.
    @test F.field_value(extreme,0.,1e307,1e307)==F.GMSH_MAX_SIZE
    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        field=F.PostViewField(coordinates,fill(2.,3);triangles=cell,use_closest=false)
        query(f)=F.field_value(f,.2,.2,.4)
        for _ in 1:3;query(field);end
        @test @allocated(query(field))==0
    end
end
end
