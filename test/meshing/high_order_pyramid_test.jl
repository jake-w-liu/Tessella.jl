module HighOrderPyramidTests

using Tessella, Test, Random
using LinearAlgebra: I
const HO=Tessella.HighOrder
const Q=Rational{BigInt}
const REFERENCE=((-1.,-1.,0.),(1.,-1.,0.),(1.,1.,0.),(-1.,1.,0.),(0.,0.,1.),
    (0.,-1.,0.),(-1.,0.,0.),(-.5,-.5,.5),(1.,0.,0.),(.5,-.5,.5),
    (0.,1.,0.),(.5,.5,.5),(-.5,.5,.5),(0.,0.,0.))
const EXPONENTS=Tuple((i,j,k) for i in 0:2 for j in 0:2 for k in 0:2-max(i,j))
cert(v,orientation=1)=HO._p2_pyramid_jacobian_certify(v,orientation,"independent-pyr14-test")

# Independent rational space from pinned pyramidalBasis.cpp/BergotBasis.cpp.
# Its cardinal interpolation is determined here by exact matrix elimination,
# without the production derivative or tensor-Bernstein weight tables.
function monomial(point,exponent)
    u,v,w=point;i,j,k=exponent
    w==1 && max(i,j)>0 && return zero(Q)
    return u^i*v^j*w^k*(1-w)^(-min(i,j))
end
function inverse_exact(input)
    n=size(input,1)
    matrix=hcat(input,Matrix{Q}(I,n,n))
    for column in 1:n
        pivot=findfirst(row->matrix[row,column]!=0,column:n)
        pivot===nothing && error("singular independent interpolation matrix")
        row=column+pivot-1
        matrix[column,:],matrix[row,:]=matrix[row,:],matrix[column,:]
        matrix[column,:]./=matrix[column,column]
        for target in 1:n
            target==column && continue
            factor=matrix[target,column]
            factor==0 || (matrix[target,:].-=factor.*matrix[column,:])
        end
    end
    return matrix[:,n+1:2n]
end
const INVERSE=inverse_exact(Q[monomial(Q.(p),e) for p in REFERENCE,e in EXPONENTS])
function physical_coefficients(v)
    values=Q[Q(v[node][d])-Q(v[1][d]) for node in 1:14,d in 1:3]
    return INVERSE*values
end
function gradient(point,exponent,axis)
    r,s,w=point;i,j,k=exponent;m=max(i,j);q=1-w
    if axis==1
        return i==0 ? zero(Q) : i*r^(i-1)*s^j*q^(m-1)*w^k
    elseif axis==2
        return j==0 ? zero(Q) : j*r^i*s^(j-1)*q^(m-1)*w^k
    end
    first=k==0 ? zero(Q) : k*r^i*s^j*q^m*w^(k-1)
    second=min(i,j)==0 ? zero(Q) : min(i,j)*r^i*s^j*q^(m-1)*w^k
    return first+second
end
function exact_jacobian(coefficients,r,s,w)
    point=(Q(r),Q(s),Q(w))
    columns=ntuple(axis->ntuple(d->sum(coefficients[k,d]*gradient(point,e,axis)
        for (k,e) in enumerate(EXPONENTS)),3),3)
    a,b,c=columns
    return a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+
           a[3]*(b[1]*c[2]-b[2]*c[1])
end
function bernstein_value(degree,index,x)
    return binomial(degree,index)*x^index*(1-x)^(degree-index)
end
const B6_INVERSE=inverse_exact(Q[bernstein_value(5,index,big(row)//big(5))
    for row in 0:5,index in 0:5])
const B4_INVERSE=inverse_exact(Q[bernstein_value(3,index,big(row)//big(3))
    for row in 0:3,index in 0:3])
function independent_bernstein(v)
    physical=physical_coefficients(v)
    values=Q[exact_jacobian(physical,2big(i)//big(5)-1,2big(j)//big(5)-1,
        big(k)//big(3)) for i in 0:5,j in 0:5,k in 0:3]
    first=Q[sum(B6_INVERSE[i,a]*values[a,j,k] for a in 1:6)
        for i in 1:6,j in 1:6,k in 1:4]
    second=Q[sum(B6_INVERSE[j,a]*first[i,a,k] for a in 1:6)
        for i in 1:6,j in 1:6,k in 1:4]
    return Q[sum(B4_INVERSE[k,a]*second[i,j,a] for a in 1:4)
        for i in 1:6,j in 1:6,k in 1:4]
end
function production_bernstein(v)
    numerator=HO._p2pyr_exact_coefficients(v)
    coordinate_denominator=foldl(lcm,(denominator(Q(value)) for point in v for value in point);init=big(1))
    return Q[numerator[i+1,j+1,k+1]//
        (64binomial(5,i)*binomial(5,j)*binomial(3,k)*coordinate_denominator^3)
        for i in 0:5,j in 0:5,k in 0:3]
end
function bernstein_evaluate(coefficients,r,s,w)
    x=(Q(r)+1)/2;y=(Q(s)+1)/2;z=Q(w)
    return sum(coefficients[i+1,j+1,k+1]*bernstein_value(5,i,x)*
        bernstein_value(5,j,y)*bernstein_value(3,k,z)
        for i in 0:5,j in 0:5,k in 0:3)
end
contains(interval,value)=
    (interval.lo==-Inf || (isfinite(interval.lo) && Q(interval.lo)<=value)) &&
    (interval.hi==Inf || (isfinite(interval.hi) && Q(interval.hi)>=value))

const WARPED=ntuple(14) do node
    u,v,w=REFERENCE[node]
    (u+u*w/32,v+v*w/64,w+(u*u+v*v)/128)
end
const EXACT_FALLBACK=ntuple(14) do node
    u,v,w=REFERENCE[node]
    (u+v+w,u-v+w,u+v-w)
end
const PRECISION_FOLD=((1e16,0.,0.),(1e16,4.,0.),(1e16,4.,1.),(1e16,0.,1.),
    (1e16+2,2.,.5),(1e16,2.,0.),(1e16,0.,.5),(1e16,1.,.25),
    (1e16,4.,.5),(1e16,3.,.25),(1e16,2.,1.),(1e16,3.,.75),
    (1e16,1.,.75),(1e16,2.,.5))

@testset "Pyramid14 whole rational reference map" begin
    @test length(EXPONENTS)==14
    @test Tessella.Elements.lagrange_nodes(14)==hcat(collect.(REFERENCE)...)
    @test cert(REFERENCE)===nothing
    @test cert(WARPED)===nothing
    reflected=ntuple(i->(-REFERENCE[i][1],REFERENCE[i][2],REFERENCE[i][3]),14)
    @test cert(reflected,-1)===nothing
    @test_throws ArgumentError cert(reflected)
    @test_throws ArgumentError cert(REFERENCE,0)
    @test_throws ArgumentError cert(ntuple(i->(NaN,0.,0.),14))
    @test_throws ArgumentError cert(ntuple(i->(Inf,0.,0.),14))
    degenerate=ntuple(i->(REFERENCE[i][1],REFERENCE[i][2],0.),14)
    @test_throws ArgumentError cert(degenerate)
    @test all(iszero,HO._p2pyr_exact_coefficients(degenerate))
    for v in (REFERENCE,WARPED,reflected,PRECISION_FOLD,EXACT_FALLBACK)
        independent=independent_bernstein(v)
        production=production_bernstein(v)
        @test size(production)==(6,6,4)
        for i in eachindex(independent)
            @test production[i]==independent[i]
        end
    end
    physical=physical_coefficients(PRECISION_FOLD)
    # The actual captured public Pyr5->Pyr14 support coordinates produce
    # X-offset=4w^2-2w, Y=2+2u, Z=1/2+v/2. P1 had constant detJ=2.
    @test exact_jacobian(physical,0,0,1//8)==-1
    @test exact_jacobian(physical,0,0,1//2)==2
    @test exact_jacobian(physical,0,0,1)==6
    @test_throws r"full P2 Pyramid14 Jacobian positivity is not certified" cert(PRECISION_FOLD)
    folded_reflection=ntuple(i->(-PRECISION_FOLD[i][1],PRECISION_FOLD[i][2],PRECISION_FOLD[i][3]),14)
    @test_throws ArgumentError cert(folded_reflection,-1)

    # Correlated derivative components of this positive affine map make the
    # interval determinant inconclusive. It exercises the actual exact path.
    @test !(HO._p2pyr_interval_bound(EXACT_FALLBACK,true).lo>0)
    @test cert(EXACT_FALLBACK)===nothing
    @test all(==(4),independent_bernstein(EXACT_FALLBACK))

    rng=MersenneTwister(0x14be551e)
    for _ in 1:12
        v=ntuple(i->ntuple(d->REFERENCE[i][d]+(rand(rng)-.5)/128,3),14)
        coefficients=production_bernstein(v)
        @test all(>(0),coefficients)
        @test cert(v)===nothing
        physical=physical_coefficients(v)
        bound=HO._p2pyr_interval_bound(v)
        normalized=HO._p2pyr_interval_bound(v,true)
        powers=ntuple(d->1-last(frexp(maximum(abs(p[d]) for p in v))),3)
        power=sum(powers)
        scale=power>=0 ? (big(1)<<power)//big(1) : big(1)//(big(1)<<(-power))
        for r in (-1//1,-1//3,1//3,1//1),s in (-1//1,1//5,1//1),w in (0//1,1//7,4//5,1//1)
            actual=exact_jacobian(physical,r,s,w)
            @test bernstein_evaluate(coefficients,r,s,w)==actual
            @test contains(bound,actual)
            @test contains(normalized,actual*scale)
        end
    end
end

@testset "Pyramid14 interval and exact extreme scales" begin
    for scale in (1e-200,1e200)
        v=ntuple(i->ntuple(d->REFERENCE[i][d]*scale,3),14)
        bound=HO._p2pyr_interval_bound(v)
        @test !(bound.lo>0)
        @test HO._p2pyr_interval_bound(v,true).lo>0
        @test cert(v)===nothing
        @test all(>(0),HO._p2pyr_exact_coefficients(v))
        reflected=ntuple(i->(-v[i][1],v[i][2],v[i][3]),14)
        @test cert(reflected,-1)===nothing
        @test_throws ArgumentError cert(v,-1)
    end
    anisotropic=ntuple(i->(REFERENCE[i][1]*1e-200,
        REFERENCE[i][2]*1e200,REFERENCE[i][3]),14)
    @test HO._p2pyr_interval_bound(anisotropic,true).lo>0
    @test cert(anisotropic)===nothing
    @test all(>(0),HO._p2pyr_exact_coefficients(anisotropic))
    values=(0.,nextfloat(0.),floatmin(Float64),.5,-2.,floatmax(Float64))
    for a in values,b in values
        difference=HO._p2_jac_interval_difference(a,b)
        @test contains(difference,Q(a)-Q(b))
        ia=HO._P2JacInterval(a,a);ib=HO._P2JacInterval(b,b)
        @test contains(HO._p2_jac_interval_add(ia,ib),Q(a)+Q(b))
        @test contains(HO._p2_jac_interval_sub(ia,ib),Q(a)-Q(b))
        @test contains(HO._p2_jac_interval_mul(ia,ib),Q(a)*Q(b))
        for k in (-4.,.25,0.,1.,2.)
            @test contains(HO._p2_jac_interval_scale(ia,k),Q(a)*Q(k))
        end
    end
    zero=HO._p2_jac_interval_mul(HO._P2_JAC_ZERO,HO._P2_JAC_UNKNOWN)
    @test zero.lo==zero.hi==0
    @test HO._p2_jac_interval_bounds(NaN,NaN)==HO._P2_JAC_UNKNOWN

    # Exact rational values validate binary normalization itself, including
    # exponent ranges that cannot be represented by a Float64 scale factor.
    powers=(-2098,-1074,-1022,-1,0,1,1022,1074,2098)
    for power in powers
        scale=power>=0 ? (big(1)<<power)//big(1) : big(1)//(big(1)<<(-power))
        for value in (values...,-nextfloat(0.),-floatmin(Float64),-floatmax(Float64))
            @test contains(HO._p2_jac_interval_scaled_value(value,power),Q(value)*scale)
            point=HO._P2JacInterval(value,value)
            @test contains(HO._p2_jac_interval_ldexp(point,power),Q(value)*scale)
        end
        for a in values,b in values
            @test contains(HO._p2pyr_interval_difference(a,b,power),(Q(a)-Q(b))*scale)
        end
        for (lo,hi) in ((-nextfloat(0.),nextfloat(0.)),
                        (-floatmax(Float64),floatmax(Float64)),(-2.,.5))
            scaled=HO._p2_jac_interval_ldexp(HO._P2JacInterval(lo,hi),power)
            @test contains(scaled,Q(lo)*scale)
            @test contains(scaled,Q(hi)*scale)
        end
    end
end

function boxed(value)
    value===Core.Box && return 1
    value isa GlobalRef && return value.mod===Core && value.name===:Box ? 1 : 0
    value isa Expr && return sum(boxed,value.args;init=0)
    value isa Core.CodeInfo && return sum(boxed,value.code;init=0)
    return 0
end
function deliberate_box(value)
    local_value=value
    capture=()->local_value
    local_value+=1
    return capture()
end
function warmed_bytes(v)
    cert(v)
    return @allocated cert(v)
end
@testset "Pyramid14 bounded certificate resources" begin
    @test isbitstype(HO._P2JacInterval)
    @test size(HO._p2pyr_exact_coefficients(REFERENCE))==(6,6,4)
    for value in (REFERENCE,WARPED,
        ntuple(i->ntuple(d->REFERENCE[i][d]*1e-200,3),14),
        ntuple(i->ntuple(d->REFERENCE[i][d]*1e200,3),14),
        ntuple(i->(REFERENCE[i][1]*1e-200,REFERENCE[i][2]*1e200,REFERENCE[i][3]),14))
        warmed_bytes(value)
        @test warmed_bytes(value)==0
    end
    @test boxed(Base.uncompressed_ast(only(methods(deliberate_box))))>0
    for symbol in names(HO;all=true)
        name=String(symbol)
        (startswith(name,"_p2pyr_") || startswith(name,"_p2_pyramid_") ||
         startswith(name,"_p2_jac_interval_")) || continue
        value=getfield(HO,symbol)
        value isa Function || continue
        for method in methods(value)
            @test boxed(Base.uncompressed_ast(method))==0
        end
    end
end

end
