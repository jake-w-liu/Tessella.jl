module PostViewThinGeometryChecks
using Test,Tessella
using LinearAlgebra: I
const F=Tessella.SizeField
const P=Tessella.Post

function field(msh,kind,coordinates,values)
    basis=Tessella.API._native_view_writer_basis(msh)
    scheme=P.PosScheme(basis.coefficients,basis.exponents,
        basis.coefficients,basis.exponents)
    record=P.PosElement(kind,coordinates,reshape(values,1,:);scheme)
    F.PostViewField([record];crop_negative=false,use_closest=false)
end

@testset "Scheme inversion retains thin intrinsic directions and actual normal offsets" begin
    h=1e-200
    triangle=field(2,:triangle,[0. 1. 0.;0. 0. h;0. 0. 0.],[0.,1.,2.])
    @test F.field_value(triangle,.25,.25h,0.)≈.75
    # The artificial normal is (0,0,h), so these reference offsets lie on
    # opposite sides of the existing 1e-6 normal-coordinate tolerance.
    @test F.field_value(triangle,.25,.25h,5e-7h)≈.75
    @test F.field_value(triangle,.25,.25h,2e-6h)==F.GMSH_MAX_SIZE
    @test F.field_value(triangle,.75,.75h,0.)==F.GMSH_MAX_SIZE
    tetrahedron=field(4,:tetrahedron,
        [0. 1. 0. 0.;0. 0. h 0.;0. 0. 0. h],[0.,1.,2.,3.])
    @test F.field_value(tetrahedron,.25,.25h,.25h)≈1.5
    # Column scaling alone cannot represent this determinant. The exact
    # fallback distinguishes its positive rank from a genuinely singular map.
    skew=field(4,:tetrahedron,
        [0. 1. 1. 1.;0. 0. h 0.;0. 0. 0. h],[0.,1.,2.,3.])
    @test F.field_value(skew,.75,.25h,.25h)≈1.5
    @test_throws ArgumentError field(4,:tetrahedron,
        [0. 1. 1. 1.;0. 0. h 0.;0. 0. 0. 0.],[0.,1.,2.,3.])
end

@testset "Regular curved geometry retains normals below the Float64 product range" begin
    # dX/du=3u²+alpha>0 for every real u. Only the artificial normal
    # product underflows at the center; the physical curve is regular.
    alpha=ldexp(1.,-700)
    scheme=P.PosScheme(ones(1,1),[1. 0. 0.],reshape([alpha,1.],1,2),
        [1. 0. 0.;3. 0. 0.])
    record=P.PosElement(:line,reshape([1.,0.,0.],3,1),ones(1,1);scheme)
    curve=F.PostViewField([record];crop_negative=false,use_closest=false)
    @test F.field_value(curve,0.,0.,0.)==0.
    for u in (.25,-.25,ldexp(1.,-280),-ldexp(1.,-280))
        x=u^3+alpha*u
        @test F.field_value(curve,x,0.,0.)≈u
    end
    @test F.field_value(curve,0.,0.,nextfloat(0.))==F.GMSH_MAX_SIZE
    @test F.field_value(curve,0.,alpha/2e6,0.)==0.
    @test F.field_value(curve,0.,alpha/2e5,0.)==F.GMSH_MAX_SIZE
    u=ldexp(1.,-280)
    @test F.field_value(curve,u^3+alpha*u,alpha/2e6,0.)≈u
    @test F.field_value(curve,2.,0.,0.)==F.GMSH_MAX_SIZE
    # Surface normals can underflow too, with two genuine small tangent
    # directions. The same intrinsic-rank and virtual-normal certificate
    # must distinguish a regular chart from a singular one.
    geometry=[alpha 1. 0. 0.;0. 0. alpha 1.;0. 0. 0. 0.]
    exponents=[1. 0. 0.;3. 0. 0.;0. 1. 0.;0. 3. 0.]
    record=P.PosElement(:quadrangle,Matrix{Float64}(I,3,3),ones(1,1);
        scheme=P.PosScheme(reshape([1.,2.],1,2),[1. 0. 0.;0. 1. 0.],geometry,exponents))
    surface=F.PostViewField([record];crop_negative=false,use_closest=false)
    @test F.field_value(surface,0.,0.,0.)==0.
    @test F.field_value(surface,0.,0.,nextfloat(0.))==F.GMSH_MAX_SIZE
    @test F.field_value(surface,.25^3+alpha*.25,.125^3+alpha*.125,0.)≈.5
    u=ldexp(1.,-280);v=u/2
    @test F.field_value(surface,u^3+alpha*u,v^3+alpha*v,0.)≈u+2v
end

@testset "Public view probes retain a regular virtual-normal geometry" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        alpha=ldexp(1.,-700);view=A.view.add("regular thin curved line")
        A.view.set_interpolation_matrices(view,"Line",1,[1.],[1.,0.,0.],
            2,[0.,0.,alpha,1.],[1.,0.,0.,3.,0.,0.])
        A.view.add_list_data(view,"SL",1,[0.,1.,0.,0.,0.,0.,1.])
        @test first(A.view.probe(view,0.,0.,0.))==[0.]
        @test first(A.view.probe(view,0.,0.,0.;gradient=true))≈[ldexp(1.,700),0.,0.]
        u=ldexp(1.,-280)
        @test first(A.view.probe(view,u^3+alpha*u,0.,0.))≈[u]
    finally
        A.finalize()
    end
end
end
