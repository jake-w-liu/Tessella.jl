module PostViewNonpartitionGeometryChecks
using Test,Tessella
using Tessella.Post: PosScheme,PosElement
using Tessella.SizeField: PostViewField,field_value,GMSH_MAX_SIZE
const A=Tessella.API
const E=[0. 0. 0.;1. 0. 0.]
const C=[.5 -.5;.5 .5]

function line_field(coordinates,geometry,exponents;native=false)
    element=PosElement(:line,coordinates,reshape([1.,2.],1,2);
        scheme=PosScheme(C,E,geometry,exponents))
    PostViewField([element];use_closest=false,_native_geometry_partition=native)
end

@testset "Nonpartition geometry uses its actual map, Jacobian and BVH" begin
    field=line_field([10. 11.;0. 0.;0. 0.],.1C,E)
    @test field_value(field,1.025,0.,0.)≈1.25
    @test field_value(field,10.25,0.,0.)==GMSH_MAX_SIZE
    @test field_value(field,.9,0.,0.)==GMSH_MAX_SIZE
    @test only(field.cell_hi)[1]<2.
    # A single control coordinate is sufficient for X(u)=10u. It does not
    # define a degenerate map simply because the control hull is a point.
    singleton=line_field(reshape([10.,0.,0.],3,1),ones(1,1),[1. 0. 0.])
    @test field_value(singleton,2.5,0.,0.)≈1.625
    @test field_value(singleton,-10.,0.,0.)≈1.
    @test field_value(singleton,10.,0.,0.)≈2.
    @test field_value(singleton,11.,0.,0.)==GMSH_MAX_SIZE
    # Every number below is binary-exact. Normalizing about the original
    # control coordinate would overflow even though the actual-map query
    # and derivative have well represented local coordinates.
    big=ldexp(1.,1000);small=ldexp(1.,-1000)
    thin=line_field([big big;0. small;0. 0.],.125C,E)
    @test field_value(thin,ldexp(1.,997),ldexp(1.,-1005),0.)≈1.25
    @test field_value(thin,big,small/4,0.)==GMSH_MAX_SIZE
    lower,upper=Tessella.SizeField._pv_scheme_geometry_bounds(only(thin.scheme_cells))
    @test lower[1]==upper[1]==ldexp(1.,997)
    @test abs(only(thin.cell_lo)[2])<small
    # Exact coefficient identity, rather than a tolerance, distinguishes an
    # explicit small defect from the native partition semantics.
    base=ldexp(1.,40);defect=ldexp(1.,-30)
    near=copy(C);near[1,1]+=defect
    shifted=line_field([base base+1;0. 0.;0. 0.],near,E)
    @test field_value(shifted,base+1024+.25,0.,0.)≈1.25
    @test field_value(shifted,base+.25,0.,0.)==GMSH_MAX_SIZE
    # Duplicate exponent columns still certify the exact unity identity.
    duplicate=line_field([base base+1;0. 0.;0. 0.],
        [.25 .25 -.5;.25 .25 .5],[0. 0. 0.;0. 0. 0.;1. 0. 0.])
    @test only(duplicate.scheme_cells).geometry===nothing
    @test field_value(duplicate,base+.25,0.,0.)≈1.25
    @test Tessella.SizeField._pv_scheme_with_values(only(field.scheme_cells),
        reshape([3.,4.],1,2)).geometry===only(field.scheme_cells).geometry

    # The producer explicitly declares its native rounded Lagrange table.
    # That provenance preserves its centered geometry at large translations.
    basis=A._native_view_writer_basis(28)
    nodes=Tessella.Elements.lagrange_nodes(28)
    coords=zeros(3,size(nodes,2));coords[1,:].=base .+(nodes[1,:].+1)./2
    values=reshape(coords[1,:].-base,1,:)
    native=PostViewField([PosElement(:line,coords,values;
        scheme=PosScheme(basis.coefficients,basis.exponents,
                         basis.coefficients,basis.exponents))];
        use_closest=false,_native_geometry_partition=true)
    @test only(native.scheme_cells).geometry===nothing
    @test field_value(native,base+.25,0.,0.)≈.25
    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        query(f,p)=field_value(f,p...)
        for _ in 1:3;query(field,(1.025,0.,0.));query(native,(base+.25,0.,0.));end
        @test @allocated(query(field,(1.025,0.,0.)))==0
        @test @allocated(query(native,(base+.25,0.,0.)))==0
    end
end

@testset "Public nonpartition probes preserve geometry, gradients and steps" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        exponents=[0.,0.,0.,1.,0.,0.]
        view=A.view.add("shifted custom map")
        A.view.set_interpolation_matrices(view,"Line",2,vec(transpose(C)),exponents,
            2,vec(transpose(.1C)),exponents)
        A.view.add_list_data(view,"SL",1,[10.,11.,0.,0.,0.,0.,1.,2.,3.,5.])
        @test first(A.view.probe(view,1.025,0.,0.;step=0))≈[1.25]
        @test first(A.view.probe(view,1.025,0.,0.;step=1))≈[3.5]
        @test first(A.view.probe(view,1.025,0.,0.;gradient=true,step=0))≈[10.,0.,0.]
        @test first(A.view.probe(view,1.025,0.,0.;gradient=true,step=1))≈[20.,0.,0.]
        @test isempty(first(A.view.probe(view,10.25,0.,0.)))
        singleton=A.view.add("one geometry control")
        A.view.set_interpolation_matrices(singleton,"Line",2,vec(transpose(C)),exponents,
            1,[1.],[1.,0.,0.])
        A.view.add_list_data(singleton,"SL",1,[10.,0.,0.,1.,2.])
        @test first(A.view.probe(singleton,2.5,0.,0.))≈[1.625]
        @test first(A.view.probe(singleton,2.5,0.,0.;gradient=true))≈[.05,0.,0.]
        thin=A.view.add("finite rebased thin map")
        A.view.set_interpolation_matrices(thin,"Line",2,vec(transpose(C)),exponents,
            2,vec(transpose(.125C)),exponents)
        big=ldexp(1.,1000);small=ldexp(1.,-1000)
        A.view.add_list_data(thin,"SL",1,[big,big,0.,small,0.,0.,1.,2.])
        @test first(A.view.probe(thin,ldexp(1.,997),ldexp(1.,-1005),0.))≈[1.25]
        @test first(A.view.probe(thin,ldexp(1.,997),ldexp(1.,-1005),0.;gradient=true))≈
            [0.,ldexp(1.,1003),0.]
    finally
        A.finalize()
    end
end
end

module PostViewSchemeCertificateBoundsChecks
using Test,Tessella
using Tessella.Post: PosElement,PosScheme
using Tessella.SizeField: PostViewField,field_value,GMSH_MAX_SIZE
@testset "Scheme BVH encloses the permitted normal tube and intrinsic extension" begin
    value_coeff=reshape([1.],1,1);value_exp=zeros(1,3)
    line=PostViewField([PosElement(:line,[0. 1.;0. 0.;0. 0.],reshape([2.],1,1);
        scheme=PosScheme(value_coeff,value_exp,[1. -1. -1.;0. 1. 1.],
            [0. 0. 0.;1. 0. 0.;1001. 0. 0.]))];use_closest=false)
    @test field_value(line,2.,0.,.5)≈2.
    @test field_value(line,2.,0.,2.)==GMSH_MAX_SIZE
    coordinates=[0. 1. 0.;0. 0. 1.;0. 0. 0.]
    geometry=[1. -1. -1. -1. -1.;0. 1. 1. 0. 0.;0. 0. 0. 1. 1.]
    exponents=[0. 0. 0.;1. 0. 0.;1001. 0. 0.;0. 1. 0.;0. 1001. 0.]
    surface=PostViewField([PosElement(:quadrangle,coordinates,reshape([2.],1,1);
        scheme=PosScheme(value_coeff,value_exp,geometry,exponents))];use_closest=false)
    @test field_value(surface,2.,2.,.5)≈2.
    @test field_value(surface,2.,2.,2.)==GMSH_MAX_SIZE
    geometry=[1. -1. -1. -1. -1.;0. 1. 1. 0. 0.;0. 0. 0. 1. 0.;0. 0. 0. 0. 1.]
    exponents=[0. 0. 0.;1. 0. 0.;100001. 0. 0.;0. 1. 0.;0. 0. 1.]
    volume=PostViewField([PosElement(:tetrahedron,
        [0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1.],reshape([2.],1,1);
        scheme=PosScheme(value_coeff,value_exp,geometry,exponents))];use_closest=false)
    u=1+5e-7
    @test field_value(volume,u+u^100001,0.,0.)≈2.
    u=1+2e-6
    @test field_value(volume,u+u^100001,0.,0.)==GMSH_MAX_SIZE
end
@testset "Rebasing retains the original artificial-normal scale" begin
    coefficients=[.5 -.5;.5 .5];exponents=[0. 0. 0.;1. 0. 0.]
    line=PostViewField([PosElement(:line,[0. 1.;0. 0.;0. 0.],reshape([1.,2.],1,2);
        scheme=PosScheme(coefficients,exponents,.1coefficients,exponents))];use_closest=false)
    @test field_value(line,.025,0.,1e-9)≈1.25
    @test field_value(line,.025,0.,1e-8)==GMSH_MAX_SIZE
    surface=PostViewField([PosElement(:triangle,[0. 1. 0.;0. 0. 1.;0. 0. 0.],reshape([2.],1,1);
        scheme=PosScheme(reshape([1.],1,1),zeros(1,3),
            .1 * [1. -1. -1.;0. 1. 0.;0. 0. 1.],[0. 0. 0.;1. 0. 0.;0. 1. 0.]))];use_closest=false)
    @test field_value(surface,.025,.025,5e-9)≈2.
    @test field_value(surface,.025,.025,5e-8)==GMSH_MAX_SIZE
end
end
