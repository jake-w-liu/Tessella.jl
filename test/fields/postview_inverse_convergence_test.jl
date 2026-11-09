module PostViewInverseConvergenceChecks
using Test, Tessella
using LinearAlgebra: I
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

function small_parameter_field(dim,alpha)
    geometry=zeros(3,2dim)
    exponents=zeros(2dim,3)
    for axis in 1:dim
        geometry[axis,2axis-1]=alpha
        geometry[axis,2axis]=1.
        exponents[2axis-1,axis]=1.
        exponents[2axis,axis]=3.
    end
    family=dim==1 ? :line : dim==2 ? :quadrangle : :tetrahedron
    scheme=P.PosScheme(ones(1,1),[1. 0. 0.],geometry,exponents)
    record=P.PosElement(family,Matrix{Float64}(I,3,3),ones(1,1);scheme)
    return F.PostViewField([record];crop_negative=false,use_closest=false)
end

@testset "small world residuals certify the actual intrinsic parameter" begin
    # Each coordinate X(t)=t^3+alpha*t is globally regular: X'(t)>0.
    # Its finite artificial normals do not require the virtual-normal path.
    # An absolute world-residual floor can nevertheless accept a parameter
    # many orders of magnitude too large near the nearly flat center.
    for (alpha,u) in ((1e-100,1e-40),(ldexp(1.,-400),ldexp(1.,-150)))
        field=small_parameter_field(1,alpha)
        @test isapprox(F.field_value(field,u^3+alpha*u,0.,0.),u;rtol=1e-8,atol=0.)
    end
    # Ordinary neighboring coordinates can converge first while the tiny
    # coordinate still needs more than 64 safeguarded Newton iterations.
    for dim in 1:3, v in (.125,.1)
        alpha=1e-100;u=1e-40
        parameters=(u,dim>=2 ? v : 0.,dim>=3 ? .25 : 0.)
        query=map(t->t^3+alpha*t,parameters)
        field=small_parameter_field(dim,alpha)
        @test isapprox(F.field_value(field,query...),u;rtol=1e-8,atol=0.)
    end
    field=small_parameter_field(1,1e-100)
    @test F.field_value(field,0.,0.,0.)==0.
    @test F.field_value(field,2.,0.,0.)==F.GMSH_MAX_SIZE
end

@testset "public probes retain small intrinsic parameters with finite normals" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        for (alpha,u) in ((1e-100,1e-40),(ldexp(1.,-400),ldexp(1.,-150)))
            view=A.view.add("regular finite-normal curved line")
            A.view.set_interpolation_matrices(view,"Line",1,[1.],[1.,0.,0.],
                2,[0.,0.,alpha,1.],[1.,0.,0.,3.,0.,0.])
            A.view.add_list_data(view,"SL",1,[0.,1.,0.,0.,0.,0.,1.])
            actual=only(first(A.view.probe(view,u^3+alpha*u,0.,0.)))
            @test isapprox(actual,u;rtol=1e-8,atol=0.)
        end
    finally
        A.finalize()
    end
end
end
