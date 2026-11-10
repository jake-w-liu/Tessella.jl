module OCCNaturalMeasureTranslationChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function natural_solid(f,name,scale,translation)
    A.initialize()
    try
        mktempdir() do directory
            source=read(joinpath(FIXTURES,name),String)
            source=replace(source,r"(?m)^\+4 0[ \t]*\r?$"=>"+2 0")
            path=joinpath(directory,"natural_solid.brep");write(path,source)
            dimension=name=="brep_native_weighted_closed.brep" ? 2 : 3
            @test A.model.occ.import_shapes(path)==[(dimension,1)]
            if dimension==2
                shell=A.model.occ.add_surface_loop([1])
                A.model.occ.add_volume([shell];tag=1)
            end
            scale==1. || A.model.occ.dilate([(3,1)],0.,0.,0.,scale,scale,scale)
            translation==(0.,0.,0.) || A.model.occ.translate([(3,1)],translation...)
            A.model.occ.synchronize()
            return f(A.model.occ)
        end
    finally
        A.finalize()
    end
end

@testset "Natural tensor Gauss rules integrate their polynomial range" begin
    for order in (8,16,32,64,128)
        nodes,weights=M._occ_natural_gauss_rule(order)
        @test all(node->0<node<1,nodes)
        @test all(>(0),weights)
        @test sum(weights)≈1. atol=3e-15
        @test all(isapprox(2sum(weight*node^degree for (node,weight) in zip(nodes,weights)),
            2/(degree+1);rtol=3e-12,atol=3e-15) for degree in 0:2:2order-2)
    end
end

@testset "Natural moment convergence requires finite component error scales" begin
    previous=(volume=1.,absolute=(volume=1.,))
    different=(volume=2.,absolute=(volume=2.,))
    unbounded=(volume=2.,absolute=(volume=Inf,))
    @test !M._occ_natural_integrals_converged(previous,different,(:volume,))
    @test !M._occ_natural_integrals_converged(previous,unbounded,(:volume,))
    near_zero=(volume=0.,absolute=(volume=1.,))
    roundoff=(volume=1e-13,absolute=(volume=1.,))
    @test M._occ_natural_integrals_converged(near_zero,roundoff,(:volume,))
end

@testset "Natural-solid divergence moments retain their local reference" begin
    sphere_I=[8π/15,0.,0.,0.,8π/15,0.,0.,0.,8π/15]
    torus_V=π^2/8
    torus_I=[torus_V*(1/2+5/128),0.,0.,0.,torus_V*(1/2+5/128),0.,
        0.,0.,torus_V*(1+3/64)]
    for (name,volume,inertia) in (("brep_native_natural_sphere.brep",4π/3,sphere_I),
                                ("brep_native_natural_torus.brep",torus_V,torus_I)),
            scale in (1.,.001),translation in ((0.,0.,0.),(1e8,-2e8,3e8),(1e9,-2e9,3e9))
        natural_solid(name,scale,translation) do occ
            @test occ.get_mass(3,1)≈volume*scale^3 rtol=2e-10 atol=0.
            absolute_tolerance=max(1e-12*scale,4eps(maximum(abs,translation)))
            @test collect(occ.get_center_of_mass(3,1))≈collect(translation) rtol=0. atol=absolute_tolerance
            expected=inertia.*scale^5
            @test occ.get_matrix_of_inertia(3,1)≈expected rtol=2e-9 atol=2e-12*maximum(abs,expected)
        end
    end
end

@testset "A translated rational natural solid matches independent Green moments" begin
    # The cross-section is the rational cubic
    # x=(1-6u²+7u³)/(1+u³), y=(3u-3u²)/(1+u³).
    # Its sections scale by (1-z²)/2. Independent 60-digit one-dimensional
    # Green integrals give the following volume, centroid and inertia.
    volume=.060643007402725003288285737111278327
    centroid=[.28014171466987301622052448664170909,
              .16470712757482137415981366018228475,0.]
    inertia=[.0090184087660961650156951462420358435,
             .00014805187961320998327960040991700162,0.,
             .00014805187961320998327960040991700162,
             .0092088865132836452484834390001992619,0.,
             0.,0.,.00090072173574409503895408892472701189]
    for (scale,translation) in ((1.,(0.,0.,0.)),(1.,(1e8,-2e8,3e8)),
                                (.001,(0.,0.,0.)))
        natural_solid("brep_native_weighted_closed.brep",scale,translation) do occ
            @test occ.get_mass(3,1)≈volume*scale^3 rtol=2e-10 atol=0.
            absolute_tolerance=max(1e-10*scale,4eps(maximum(abs,translation)))
            @test collect(occ.get_center_of_mass(3,1))≈scale.*centroid.+collect(translation) rtol=0. atol=absolute_tolerance
            expected=inertia.*scale^5
            @test occ.get_matrix_of_inertia(3,1)≈expected rtol=2e-9 atol=2e-12*maximum(abs,expected)
        end
    end
    natural_solid("brep_native_weighted_closed.brep",1.,(0.,0.,0.)) do occ
        stretch=(1e200,1e-100,1e-100)
        occ.dilate([(3,1)],0.,0.,0.,stretch...)
        occ.synchronize()
        @test occ.get_mass(3,1)≈volume rtol=2e-10 atol=0.
        actual=occ.get_center_of_mass(3,1)
        expected=centroid.*collect(stretch)
        @test all(isapprox(actual[axis],expected[axis];rtol=2e-9,atol=1e-112)
                  for axis in 1:3)
    end
end
end
