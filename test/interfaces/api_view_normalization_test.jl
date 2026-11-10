module ViewNormalizationChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const F=Tessella.SizeField
const R=Rational{BigInt}

@testset "Public probes retain original coordinates in thin maps" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        coords=Matrix{Float64}(I,3,3)
        exponents=[1. 0. 0.;0. 1. 0.;0. 0. 1.]
        for thinness in (1e-8,1e-12,1e-14), horizontal in (.125,.3), components in (1,3,9)
            geometry=[1. 0. 0.;1. thinness 0.;0. 0. 3.]
            data=collect(1.:components)
            payload=vcat(vec(coords[1,:]),vec(coords[2,:]),vec(coords[3,:]),data)
            before=copy(payload)
            view=A.view.add("thin original query")
            A.view.set_interpolation_matrices(view,"Pyramid",1,[1.],[0.,1.,0.],
                3,vec(geometry'),vec(exponents'))
            A.view.add_list_data(view,components==1 ? "SY" : components==3 ? "VY" : "TY",1,payload)
            point=(horizontal,horizontal+thinness*.25,3horizontal)
            expected=Float64((R(point[2])-R(point[1]))/R(thinness))
            values,distance=A.view.probe(view,point...)
            @test values≈expected.*data rtol=2e-13
            @test iszero(distance)
            gradients,distance=A.view.probe(view,point...;gradient=true)
            expected_gradient=vec([-inv(thinness),inv(thinness),0.].*data')
            @test gradients≈expected_gradient rtol=2e-13
            @test iszero(distance)
            @test payload==before
            @test isempty(first(A.view.probe(view,horizontal,horizontal+2thinness,3horizontal)))
        end
    finally A.finalize() end
end

@testset "Field certificate closures retain concrete locals" begin
    for fn in (F._pv_rebased_geometry_bounds,F._pv_scheme_precise_certificate)
        @test all(methods(fn)) do method
            info=Base.uncompressed_ast(method)
            !any(info.code) do statement
                statement isa Expr && statement.head==:(=) || return false
                rhs=statement.args[2]
                rhs isa Expr && rhs.head==:call && rhs.args[1]==GlobalRef(Core,:Box)
            end
        end
    end
end
end
