module PostViewThinGeometryChecks
using Test,Tessella
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
end
