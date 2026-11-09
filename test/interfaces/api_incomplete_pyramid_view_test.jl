module IncompletePyramidViewChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const E=Tessella.Elements
const F=Tessella.SizeField
const B=Tessella.MeshFunctionSpaces

function scheme_values(basis,point)
    return setprecision(BigFloat,256) do
        u,v,w=BigFloat.(point)
        q=1-w
        terms=BigFloat[]
        for row in axes(basis.exponents,1)
            a,b,c=basis.exponents[row,:]
            push!(terms,iszero(q) && max(a,b)>0 ? BigFloat(0) :
                u^a*v^b*w^c/q^min(a,b))
        end
        Float64.(BigFloat.(basis.coefficients)*terms)
    end
end

@testset "Incomplete pyramid schemes preserve native edge traces at every order" begin
    for msh in 125:131
        spec=E.msh_spec(msh)
        basis=A._native_view_writer_basis(msh)
        @test size(basis.coefficients,1)==spec.nnodes
        @test size(basis.coefficients,2)==size(basis.exponents,1)
        for point in ((.1,-.2,.3),(.25,-1.,0.),(0.,0.,1.),(.15,.25,.75),
                      (.1*2.0^-20,.2*2.0^-20,1-2.0^-20),
                      (.1*2.0^-40,.2*2.0^-40,1-2.0^-40))
            expected=B.mesh_basis_functions(msh,collect(point),"Lagrange")[2]
            @test scheme_values(basis,point)≈expected atol=2e-8 rtol=2e-8
        end
        nodes=E.lagrange_nodes(msh)
        for node in (1,6,spec.nnodes)
            expected=zeros(spec.nnodes);expected[node]=1
            @test scheme_values(basis,Tuple(nodes[:,node]))≈expected atol=2e-8
        end
    end
end

@testset "Incomplete pyramids retain owned model views and sizing fields" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        for msh in (125,126,131)
            name="pyramid $msh"
            A.model.add(name)
            A.model.add_discrete_entity(3,7)
            coordinates=E.lagrange_nodes(msh)
            coordinates[3,6]+=.05
            tags=UInt64.(101:100+size(coordinates,2))
            A.mesh.add_nodes(3,7,tags,vec(coordinates))
            A.mesh.add_elements_by_type(7,msh,[301],tags)
            point=(.1,-.2,.3)
            # A single edge-node impulse distinguishes the native extension
            # from affine/corner projection and a fitted complete pyramid.
            values=zeros(length(tags));values[6]=1
            native=B.mesh_basis_functions(msh,collect(point),"Lagrange")[2]
            gradients=reshape(B.mesh_basis_functions(
                msh,collect(point),"GradLagrange")[2],3,:)
            expected=native[6]
            physical=coordinates*native
            jacobian=coordinates*transpose(gradients)
            expected_gradient=transpose(jacobian)\gradients[:,6]
            for kind in ("NodeData","ElementNodeData")
                view=A.view.add("edge impulse")
                A.view.add_model_data(view,0,name,kind,
                    kind=="NodeData" ? tags : [301],
                    kind=="NodeData" ? [[v] for v in values] : [values])
                actual,distance=A.view.probe(view,physical...)
                @test actual≈[expected] atol=2e-9
                @test iszero(distance)
                gradient,gradient_distance=A.view.probe(view,physical...;gradient=true)
                @test gradient≈expected_gradient atol=2e-8
                @test iszero(gradient_distance)
                # Disable value cropping: an interior edge basis can be negative.
                fieldid=A.mesh.field.add("PostView")
                A.mesh.field.set_number(fieldid,"ViewTag",view)
                A.mesh.field.set_number(fieldid,"CropNegativeValues",0)
                A.mesh.field.set_number(fieldid,"UseClosest",0)
                A.mesh.field.set_as_background_mesh(fieldid)
                field=lock(A.STATE_LOCK) do
                    A._session_size_field_locked(A._model_locked();raw=true)
                end
                @test F.field_value(field,physical...)≈expected atol=2e-9
                mktempdir() do directory
                    path=joinpath(directory,"native-pyramid.msh")
                    for binary in (0,1)
                        @testset "msh=$msh $kind binary=$binary" begin
                            A.option("PostProcessing.Binary",binary)
                            A.view.write(view,path)
                            imported=Tessella.read_mixed_msh(path)
                            @test only(imported.blocks).msh==msh
                            @test only(imported.data_sections).values==values
                            previous=A.view.get_tags()
                            A.open(path)
                            reopened=only(setdiff(A.view.get_tags(),previous))
                            result,reopened_distance=A.view.probe(reopened,physical...)
                            @test length(result)==1
                            if length(result)==1
                                @test only(result)≈expected atol=2e-9
                            end
                            @test iszero(reopened_distance)
                            derivative,derivative_distance=A.view.probe(
                                reopened,physical...;gradient=true)
                            @test length(derivative)==3
                            if length(derivative)==3
                                @test derivative≈expected_gradient atol=2e-8
                            end
                            @test iszero(derivative_distance)
                            A.model.set_current(name)
                        end
                    end
                end
            end
        end
    finally
        A.finalize()
    end
end
end
