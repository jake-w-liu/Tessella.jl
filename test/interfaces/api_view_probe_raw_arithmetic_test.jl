module ViewProbeRawArithmeticChecks
using Test,Tessella
const A=Tessella.API
const R=Rational{BigInt}

function list_record(name,h,coefficients,exponents,values;
                     geometry=[1.,-1.,0.,1.])
    coef=Float64.(coefficients);exp=Float64.(exponents);data=Float64.(values)
    geo=Float64.(geometry);geo_exp=[0.,0.,0.,1.,0.,0.]
    payload=[0.,h,0.,0.,0.,0.,data...]
    inputs=(coef,exp,data,geo,geo_exp,payload)
    snapshots=map(copy,inputs)
    tag=A.view.add(name)
    A.view.set_interpolation_matrices(tag,"Line",length(data),coef,exp,2,geo,geo_exp)
    A.view.add_list_data(tag,"SL",1,payload)
    # There is no public interpolation-matrix getter. This read-only snapshot
    # checks ownership; every evaluated value and derivative uses public probe.
    basis=deepcopy(A.SESSION_VIEWS[tag].interp)
    before=deepcopy(A.view.get_list_data(tag))
    return (;tag,inputs,snapshots,basis,before)
end

function unchanged(record)
    @test all(isequal.(record.inputs,record.snapshots))
    @test isequal(A.view.get_list_data(record.tag),record.before)
    @test isequal(A.SESSION_VIEWS[record.tag].interp,record.basis)
end

@testset "Public list probes preserve raw basis arithmetic" begin
    A.initialize()
    saved_precision=precision(BigFloat);saved_rounding=rounding(BigFloat)
    try
        A.option("General.Terminal",0)
        h=2.0^-1000;delta=nextfloat(0.0)
        @testset "mixed nodal scales and ordinary data" begin
            mixed=list_record("mixed nodal scales",h,
                [0.,-.5,.5,0.,.5,.5,1.,0.,-1.],
                [0.,0.,0.,1.,0.,0.,2.,0.,0.],[0.,delta,1.];
                geometry=[.5,-.5,.5,.5])
            @test A.view.probe(mixed.tag,h/2,0.,0.;gradient=true)==([2.0^-74,0.,0.],0.)
            @test A.view.probe(mixed.tag,h/2,0.,0.)==([1.],0.)
            unchanged(mixed)
            ordinary=list_record("ordinary scales",1.,
                [0.,-.5,.5,0.,.5,.5,1.,0.,-1.],
                [0.,0.,0.,1.,0.,0.,2.,0.,0.],[0.,1.,1.];
                geometry=[.5,-.5,.5,.5])
            @test A.view.probe(ordinary.tag,.5,0.,0.;gradient=true)==([1.,0.,0.],0.)
            unchanged(ordinary)
        end
        @testset "coefficient product before derivative folding" begin
            record=list_record("tiny derivative coefficient",h,[delta],[2.,0.,0.],[1.])
            @test A.view.probe(record.tag,h/4,0.,0.;gradient=true)==([2.0^-75,0.,0.],0.)
            @test A.view.probe(record.tag,h/4,0.,0.)==([0.],0.)
            unchanged(record)
        end
        @testset "coefficient product before a large datum" begin
            record=list_record("value rescaled by datum",1.,[delta],[1.,0.,0.],[1e308])
            expected=Float64(R(delta)*R(1e308)/2)
            @test A.view.probe(record.tag,.5,0.,0.)==([expected],0.)
            @test A.view.probe(record.tag,.5,0.,0.;gradient=true)==([2expected,0.,0.],0.)
            unchanged(record)
        end
        @testset "duplicate terms retain a small remainder" begin
            record=list_record("duplicate polynomial terms",1.,
                [1e100,1.,-1e100,0.,0.,0.,0.,0.,0.],
                [1.,0.,0.,1.,0.,0.,1.,0.,0.],[1.,0.,0.])
            @test A.view.probe(record.tag,.5,0.,0.)==([.5],0.)
            @test A.view.probe(record.tag,.5,0.,0.;gradient=true)==([1.,0.,0.],0.)
            unchanged(record)
        end
        @testset "fractional powers retain their scaled derivative" begin
            record=list_record("fractional power",h,[delta],[.5,0.,0.],[1.])
            @test A.view.probe(record.tag,h,0.,0.)==([delta],0.)
            setprecision(BigFloat,192) do
                setrounding(BigFloat,RoundDown) do
                    @test A.view.probe(record.tag,h,0.,0.;gradient=true)==([2.0^-75,0.,0.],0.)
                    @test precision(BigFloat)==192
                    @test rounding(BigFloat)==RoundDown
                end
            end
            unchanged(record)
        end
        @testset "an underflowed monomial has a finite physical derivative" begin
            record=list_record("degree1100",h,[1.],[1100.,0.,0.],[1.])
            @test A.view.probe(record.tag,h/2,0.,0.)==([0.],0.)
            @test A.view.probe(record.tag,h/2,0.,0.;gradient=true)==([ldexp(1100.,-99),0.,0.],0.)
            unchanged(record)
        end
        @testset "a huge exponent at one requires no huge rational power" begin
            exponent=1e300
            record=list_record("huge exponent endpoint",1.,[delta],[exponent,0.,0.],[1.])
            @test A.view.probe(record.tag,1.,0.,0.)==([delta],0.)
            @test A.view.probe(record.tag,1.,0.,0.;gradient=true)==([Float64(R(delta)*R(exponent)),0.,0.],0.)
            unchanged(record)
        end
        @testset "the original physical query retains small endpoint offsets" begin
            for (name,payload,point,value,gradient) in (
                    ("near last node",[-1e16,0.,0.,0.,0.,0.,1e16,0.],-.25,.25,-1.),
                    ("near first node",[0.,1e16,0.,0.,0.,0.,0.,1e16],.25,.25,1.))
                tag=A.view.add(name)
                A.view.add_list_data(tag,"SL",1,payload)
                before=deepcopy(A.view.get_list_data(tag))
                @test A.view.probe(tag,point,0.,0.)==([value],0.)
                @test A.view.probe(tag,point,0.,0.;gradient=true)==([gradient,0.,0.],0.)
                @test isequal(A.view.get_list_data(tag),before)
            end
        end
        @test precision(BigFloat)==saved_precision
        @test rounding(BigFloat)==saved_rounding
    finally
        A.finalize()
    end
end
end
