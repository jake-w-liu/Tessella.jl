module QuadTriNoNewJacobianTests

using Tessella, Test, Random
using Tessella.Model
const CUBE=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
            (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.))
const POS=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
           (0,0,1),(1,0,1),(1,1,1),(0,1,1))
cert(v,ori=1)=Model._extrude_nonew_hex_jacobian_certify(v,ori,"jac-certificate")
function exact_samples(v)
    coords=ntuple(k->ntuple(d->Rational{BigInt}(v[k][d]),3),8)
    return ntuple(27) do index
        r=big((index-1)%3)//big(2);s=big((index-1)÷3%3)//big(2)
        t=big((index-1)÷9)//big(2)
        derivatives=ntuple(3) do axis
            ntuple(3) do d
                total=big(0)//big(1)
                for k in 1:8
                    x,y,z=POS[k]
                    factor=axis==1 ? (x==0 ? -1 : 1) : (x==0 ? 1-r : r)
                    factor*=axis==2 ? (y==0 ? -1 : 1) : (y==0 ? 1-s : s)
                    factor*=axis==3 ? (z==0 ? -1 : 1) : (z==0 ? 1-t : t)
                    total+=factor*coords[k][d]
                end
                total
            end
        end
        a,b,c=derivatives
        a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+
            a[3]*(b[1]*c[2]-b[2]*c[1])
    end
end
function exact_bernstein(v)
    samples=exact_samples(v)
    function transform(values,stride)
        ntuple(27) do index
            digit=(index-1)÷stride%3
            digit==1 ? 2values[index]-(values[index-stride]+values[index+stride])/2 :
                values[index]
        end
    end
    return transform(transform(transform(samples,1),3),9)
end
function interval_bernstein(v)
    edges=Model._extrude_nonew_jac_interval_edges(v)
    samples=ntuple(i->Model._extrude_nonew_jac_sample(edges,i),27)
    first=Model._extrude_nonew_jac_bernstein_axis(samples,1)
    second=Model._extrude_nonew_jac_bernstein_axis(first,3)
    return Model._extrude_nonew_jac_bernstein_axis(second,9)
end
contains(interval,value)=
    (interval.lo==-Inf || (isfinite(interval.lo) && Rational{BigInt}(interval.lo)<=value)) &&
    (interval.hi==Inf || (isfinite(interval.hi) && Rational{BigInt}(interval.hi)>=value))

@testset "Actual Hex8 full-domain Bernstein certificate" begin
    @test cert(CUBE)===nothing
    reflected=ntuple(k->(-CUBE[k][1],CUBE[k][2],CUBE[k][3]),8)
    @test cert(reflected,-1)===nothing
    @test_throws ArgumentError cert(CUBE,-1)
    @test_throws ArgumentError cert(CUBE,0)
    @test_throws ArgumentError cert(ntuple(k->(NaN,CUBE[k][2],CUBE[k][3]),8))
    for scale in (1e-200,1e200)
        scaled=ntuple(k->ntuple(d->CUBE[k][d]*scale,3),8)
        @test cert(scaled)===nothing
        @test Model._extrude_nonew_jac_exact_positive(scaled,1)
    end
    rng=MersenneTwister(0xbef00d)
    for index in 1:40
        v=ntuple(k->ntuple(d->CUBE[k][d]+.15(rand(rng)-.5),3),8)
        exact=exact_bernstein(v);bounds=interval_bernstein(v)
        @test all(>(0),exact)
        @test Model._extrude_nonew_jac_exact_positive(v,1)
        @test cert(v)===nothing
        for k in 1:27
            @test contains(bounds[k],exact[k])
        end
    end
    # The actual public positive-volume retained Hex8, not a continuous CAD
    # approximation, has a strictly negative interior trilinear Jacobian.
    folded=((1.,1.,0.),(1.,0.,0.),(0.,0.,0.),(0.,1.,0.),
        (5.5224396333287125,-4.486401190603047,-.30524078967639356),
        (5.192952385269027,-5.156913942543361,.3594777995412404),
        (4.522439633328713,-5.486401190603047,-.30524078967639356),
        (4.8519268813883984,-4.815888438662733,-.9699593788940275))
    @test !Model._extrude_nonew_jac_exact_positive(folded,1)
    error=try cert(folded);nothing catch caught;caught end
    @test error isa ArgumentError
    @test occursin("full P1 Hex8 Jacobian positivity is not certified",error.msg)
    @test any(<=(0),exact_bernstein(folded))
    # A zero Jacobian remains unproved, even if all exact coefficients vanish.
    flat=ntuple(k->(CUBE[k][1],CUBE[k][2],0.),8)
    @test_throws ArgumentError cert(flat)
end

@testset "Outward intervals retain exact extreme arithmetic" begin
    values=(0.0,0.1,-0.1,1e-200,-1e-200,1e200,-1e200,floatmax(Float64))
    for a in values,b in values
        left=Model._ExtrudeNoNewJacInterval(a,a)
        right=Model._ExtrudeNoNewJacInterval(b,b)
        exact_left=Rational{BigInt}(a);exact_right=Rational{BigInt}(b)
        @test contains(Model._extrude_nonew_jac_add(left,right),exact_left+exact_right)
        @test contains(Model._extrude_nonew_jac_sub(left,right),exact_left-exact_right)
        @test contains(Model._extrude_nonew_jac_mul(left,right),exact_left*exact_right)
        @test contains(Model._extrude_nonew_jac_scale(left,-0.5),-exact_left/2)
    end
end

function boxes(value)
    value===Core.Box && return 1
    value isa GlobalRef && return value.mod===Core && value.name===:Box ? 1 : 0
    value isa Core.CodeInfo && return sum(boxes,value.code;init=0)
    value isa Expr && return sum(boxes,value.args;init=0)
    value isa QuoteNode && return boxes(value.value)
    return 0
end
function deliberate_box_probe(value)
    captured=value
    reader=()->captured
    captured=value+1
    return reader
end
@testset "Hex8 Jacobian helper boxing and warmed allocation" begin
    @test boxes(Base.uncompressed_ast(only(methods(deliberate_box_probe))))>0
    cert(CUBE)
    @test (@allocated cert(CUBE))==0
    for name in names(Model;all=true)
        startswith(String(name),"_extrude_nonew_jac_") ||
            name===:_extrude_nonew_hex_jacobian_certify || continue
        value=getfield(Model,name)
        value isa Function || continue
        for method in methods(value)
            @test boxes(Base.uncompressed_ast(method))==0
        end
    end
end

end # module QuadTriNoNewJacobianTests
