module NativePyramidApexChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const E=Tessella.Elements
const B=Tessella.MeshFunctionSpaces
const F=Tessella.SizeField

@testset "Native pyramid probes retain the actual collapse denominator" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("native apex")
        A.model.add_discrete_entity(3,7)
        coordinates=E.lagrange_nodes(7)
        tags=UInt64.(101:105)
        A.mesh.add_nodes(3,7,tags,vec(coordinates))
        A.mesh.add_elements_by_type(7,7,[301],tags)
        for exponent in (40,48,50,53),kind in ("NodeData","ElementNodeData")
            @testset "q=2^-$exponent $kind" begin
                q=2.0^-exponent
                reference=(q/4,q/4,1-q)
                amplitude=inv(q)
                values=[amplitude,0.,0.,0.,0.]
                native=B.mesh_basis_functions(7,collect(reference),"Lagrange")[2]
                derivatives=reshape(B.mesh_basis_functions(
                    7,collect(reference),"GradLagrange")[2],3,:)
                @test native[1]*amplitude==9/64
                view=A.view.add("native corner")
                A.view.add_model_data(view,0,"native apex",kind,
                    kind=="NodeData" ? tags : [301],
                    kind=="NodeData" ? [[value] for value in values] : [values])
                actual,distance=A.view.probe(view,reference...)
                @test actual≈[9/64] atol=1e-12 rtol=1e-12
                @test iszero(distance)
                gradient,gradient_distance=A.view.probe(view,reference...;gradient=true)
                @test gradient≈amplitude.*derivatives[:,1] rtol=1e-12
                @test iszero(gradient_distance)
                # Native interpolation has a finite canonical apex convention.
                @test A.view.probe(view,0.,0.,1.)==([0.],0.)
                apex_gradient,apex_distance=A.view.probe(view,0.,0.,1.;gradient=true)
                @test apex_gradient≈fill(-amplitude/4,3) rtol=1e-12
                @test iszero(apex_distance)
                if exponent==50
                    # The public canonical field path already uses the
                    # actual collapse denominator and must retain it.
                    field=F.PostViewField(coordinates,values;
                        pyramids=[(1,2,3,4,5)],crop_negative=false,use_closest=false)
                    @test F.field_value(field,reference...)==9/64
                end
            end
        end
        @testset "native vector and tensor components" begin
            q=2.0^-50
            point=(q/4,q/4,1-q)
            derivative=reshape(B.mesh_basis_functions(
                7,collect(point),"GradLagrange")[2],3,:)[:,1]
            for components in (3,9)
                view=A.view.add("native component corner")
                values=zeros(components,5)
                values[:,1].=collect(1:components)./q
                A.view.add_model_data(view,0,"native apex","NodeData",tags,
                    [collect(values[:,node]) for node in 1:5])
                actual,distance=A.view.probe(view,point...)
                @test actual≈(9/64).*collect(1:components) rtol=1e-12
                @test iszero(distance)
                gradient,gradient_distance=A.view.probe(view,point...;gradient=true)
                expected=vcat([component.*derivative./q for component in 1:components]...)
                @test gradient≈expected rtol=1e-12
                @test iszero(gradient_distance)
            end
        end
    finally
        A.finalize()
    end
end

@testset "Explicit POS pyramid matrices retain their denominator convention" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        coordinates=E.lagrange_nodes(7)
        basis=A._native_view_writer_basis(7)
        for exponent in (40,50,53),explicit_geometry in (false,true)
            q=2.0^-exponent
            point=(q/4,q/4,1-q)
            values=[inv(q),0.,0.,0.,0.]
            view=A.view.add("explicit matrix")
            coefficients=vec(basis.coefficients')
            exponents=vec(Float64.(basis.exponents)')
            if explicit_geometry
                # Sum=2 gives the explicit map x=2r, avoiding cancellation
                # against a translated first node while testing both matrices.
                A.view.set_interpolation_matrices(view,"Pyramid",5,
                    coefficients,exponents,5,2 .* coefficients,exponents)
            else
                A.view.set_interpolation_matrices(view,"Pyramid",5,
                    coefficients,exponents)
            end
            payload=vcat(vec(coordinates[1,:]),vec(coordinates[2,:]),
                vec(coordinates[3,:]),values)
            A.view.add_list_data(view,"SY",1,payload)
            R=Rational{BigInt}
            expected=Float64((R(q)-R(point[1])-R(point[2])+
                R(point[1])*R(point[2])/R(max(1e-14,q)))/(4R(q)))
            physical=explicit_geometry ? 2 .* point : point
            actual,distance=A.view.probe(view,physical...)
            @test actual≈[expected] atol=1e-12 rtol=1e-12
            @test iszero(distance)
        end
    finally
        A.finalize()
    end
end
end
