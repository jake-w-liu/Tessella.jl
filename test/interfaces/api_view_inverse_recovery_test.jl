module ViewInverseRecoveryChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const F=Tessella.SizeField
const E=Tessella.Elements
const P=Tessella.Post
const B=Tessella.MeshFunctionSpaces
const retained_field=Ref{Any}(nothing)
const retained_point=Ref((0.,0.,0.))
function model_field(tag)
    lock(A.STATE_LOCK) do
        A._slot_capture_locked!()
        view=A._session_view_data(tag,"inverse regression")
        A._model_post_view_field(view,Int(tag);crop_negative=false,use_closest=false)
    end
end
@testset "Curved native pyramid data retain boundary-node queries" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("curved boundary")
        A.model.add_discrete_entity(3,7)
        coords=E.lagrange_nodes(125);coords[3,6]+=.05
        tags=UInt64.(101:121)
        A.mesh.add_nodes(3,7,tags,vec(coords))
        A.mesh.add_elements_by_type(7,125,[301],tags)
        values=zeros(21);values[6]=1
        point=Tuple(coords[:,6])
        for kind in ("NodeData","ElementNodeData")
            view=A.view.add("boundary impulse")
            A.view.add_model_data(view,0,"curved boundary",kind,
                kind=="NodeData" ? tags : [301],
                kind=="NodeData" ? [[value] for value in values] : [values])
            actual,distance=A.view.probe(view,point...)
            @test actual≈[1.] atol=1e-12
            @test iszero(distance)
            field=model_field(view)
            retained_field[]=field;retained_point[]=point
            @test F.field_value(field,point...)≈1. atol=1e-12
            mktempdir() do directory
                for binary in (0,1)
                    @testset "$kind binary=$binary" begin
                        A.option("PostProcessing.Binary",binary)
                        path=joinpath(directory,"boundary.msh")
                        A.view.write(view,path)
                        previous=A.view.get_tags();A.open(path)
                        imported=only(setdiff(A.view.get_tags(),previous))
                        actual,distance=A.view.probe(imported,point...)
                        @test actual≈[1.] atol=1e-12
                        @test iszero(distance)
                        @test F.field_value(model_field(imported),point...)≈1. atol=1e-12
                        A.model.set_current("curved boundary")
                    end
                end
            end
        end
    finally A.finalize() end
end
@testset "Explicit pyramid matrices retain precise physical queries" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        coords=E.lagrange_nodes(7);basis=A._native_view_writer_basis(7)
        for exponent in (40,50,53)
            q=2.0^-exponent;point=(q/4,q/4,1-q)
            values=[inv(q),0.,0.,0.,0.]
            view=A.view.add("explicit near apex")
            A.view.set_interpolation_matrices(view,"Pyramid",5,
                vec(basis.coefficients'),vec(Float64.(basis.exponents)'),
                5,vec(basis.coefficients'),vec(Float64.(basis.exponents)'))
            A.view.add_list_data(view,"SY",1,vcat(vec(coords[1,:]),vec(coords[2,:]),vec(coords[3,:]),values))
            R=Rational{BigInt}
            expected=Float64((R(q)-R(point[1])-R(point[2])+
                R(point[1])*R(point[2])/R(max(1e-14,q)))/(4R(q)))
            actual,distance=A.view.probe(view,point...)
            @test actual≈[expected] atol=1e-12
            @test iszero(distance)
            record=P.PosElement(:pyramid,coords,reshape(values,1,:);
                scheme=P.PosScheme(basis.coefficients,basis.exponents,basis.coefficients,basis.exponents))
            field=F.PostViewField([record];crop_negative=false,use_closest=false)
            @test F.field_value(field,point...)≈expected atol=1e-12
        end
        # A stored control point is outside this actual map z∈[0,.1].
        view=A.view.add("scaled explicit geometry")
        A.view.set_interpolation_matrices(view,"Pyramid",5,
            vec(basis.coefficients'),vec(Float64.(basis.exponents)'),
            5,vec((.1 .* basis.coefficients)'),vec(Float64.(basis.exponents)'))
        A.view.add_list_data(view,"SY",1,vcat(vec(coords[1,:]),vec(coords[2,:]),vec(coords[3,:]),ones(5)))
        @test A.view.probe(view,0.,0.,1.)==(Float64[],0.)
    finally A.finalize() end
end
@testset "Canonical pyramid fields retain tiny apex offsets" begin
    q=2.0^-53;point=(q/4,q/4,1-q)
    field=F.PostViewField(E.lagrange_nodes(7),[inv(q),0.,0.,0.,0.];
        pyramids=[(1,2,3,4,5)],crop_negative=false,use_closest=false)
    @test F.field_value(field,point...)==9/64
end
@testset "Concurrent precise fields preserve caller arithmetic scopes" begin
    field=retained_field[];point=retained_point[]
    saved_precision=precision(BigFloat);saved_rounding=rounding(BigFloat)
    before=copy(field.values)
    tasks=[Threads.@spawn begin
        bits=(128,192,256,320)[mod1(index,4)]
        mode=(RoundDown,RoundUp,RoundNearest,RoundToZero)[mod1(index,4)]
        setprecision(BigFloat,bits) do
            setrounding(BigFloat,mode) do
                good=F._pv_precise_scope() do
                    precision(BigFloat)==F._PV_PRECISE_INVERSE_BITS &&
                        rounding(BigFloat)==RoundNearest
                end
                for _ in 1:2
                    good &= isapprox(F.field_value(field,point...),1.;atol=1e-12)
                    good &= precision(BigFloat)==bits && rounding(BigFloat)==mode
                    yield()
                end
                good
            end
        end
    end for index in 1:8]
    @test all(fetch,tasks)
    @test precision(BigFloat)==saved_precision
    @test rounding(BigFloat)==saved_rounding
    @test isequal(field.values,before)
end
end
