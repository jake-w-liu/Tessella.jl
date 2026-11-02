module PostViewInverseConvergenceChecks
using Test, Tessella
const F=Tessella.SizeField
const P=Tessella.Post
const E=Tessella.Elements

function scheme_record(msh,coordinates,values)
    basis=Tessella.API._native_view_writer_basis(msh)
    scheme=P.PosScheme(basis.coefficients,basis.exponents,
        basis.coefficients,basis.exponents)
    family=E.msh_family(msh)==:lin ? :line : :pyramid
    return P.PosElement(family,coordinates,values;scheme)
end

function cycle_field()
    parameters=vec(E.lagrange_nodes(28)[1,:])
    coordinates=zeros(3,length(parameters))
    coordinates[1,:].=curve_x.(parameters)
    record=scheme_record(28,coordinates,reshape(parameters,1,:))
    return F.PostViewField([record];crop_negative=false,use_closest=false)
end

curve_x(u)=let t=2u-1;23t-10t^3+3t^5;end

@testset "safeguarded inversion of a globally regular curved line" begin
    # x'(u)=16+480u^2(u-1)^2 >=16 on the entire real line. Every query
    # below has one unique parameter, and the geometry has no fold or cusp.
    # Undamped Newton for x=0 nevertheless cycles forever between u=0,1.
    field=cycle_field()
    @test curve_x(.5)==0.
    @test curve_x(-1.)==-528. && curve_x(1.)==16.
    for parameter in (-1.,-.75,-.25,0.,.25,.5,.75,1.)
        @test F.field_value(field,curve_x(parameter),0.,0.)≈parameter atol=1e-10
    end
    @test F.field_value(field,17.,0.,0.)==F.GMSH_MAX_SIZE
    @test F.field_value(field,-529.,0.,0.)==F.GMSH_MAX_SIZE
    # The existing artificial-normal reference tolerance remains applicable.
    @test F.field_value(field,0.,1e-8,0.)≈.5 atol=1e-10
    @test F.field_value(field,0.,1e-2,0.)==F.GMSH_MAX_SIZE
end

@testset "the actual residual certificate retains the pyramid apex" begin
    coordinates=E.lagrange_nodes(14)
    values=reshape([2+coordinates[1,i]^2+coordinates[2,i]^2+3coordinates[3,i]^2
        for i in axes(coordinates,2)],1,:)
    field=F.PostViewField([scheme_record(14,coordinates,values)];use_closest=false)
    @test F.field_value(field,0.,0.,1.)≈5.
    @test F.field_value(field,.1,-.2,.3)≈2.32 atol=1e-10
    @test F.field_value(field,0.,0.,1.01)==F.GMSH_MAX_SIZE
end
end
