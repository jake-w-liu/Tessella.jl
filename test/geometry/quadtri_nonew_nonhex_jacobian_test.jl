module QuadTriNoNewNonHexJacobianTests

using Tessella, Test, Random
using Tessella.Model
using Tessella.MeshTypes: tet_signed_volume

const PYRAMID=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(.5,.5,1.))
const PRISM=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
             (0.,0.,1.),(1.,0.,1.),(0.,1.,1.))
const FOLDED_PYRAMID=((0.,0.,0.),(1.,0.,0.),(1.,1.,1.),
                      (0.,1.,0.),(2.,2.,2.5))
const FOLDED_PRISM=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
                    (0.,0.,1.),(-2.,0.,1.),(0.,-.5,1.))
const SINGULAR_PRISM=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
                      (0.,0.,1.),(-1.,0.,1.),(0.,-1.,1.))
const TURN=2.8
const TWISTED_PRISM=(PRISM[1],PRISM[2],PRISM[3],PRISM[4],
                     (cos(TURN),sin(TURN),1.),(-sin(TURN),cos(TURN),1.))
# Its strictly positive middle Jacobian is smaller than the interval rounding
# uncertainty. This exercises exact fallback without changing the input map.
const NEAR_SINGULAR_PRISM=(PRISM[1],PRISM[2],PRISM[3],PRISM[4],
                          (-1.,1e-12,1.),(-1e-12,-1.,1.))

reflect(v)=ntuple(i->(-v[i][1],v[i][2],v[i][3]),length(v))
pyramid_cert(v)=Model._extrude_nonew_pyramid_jacobian_certify(v,"nonhex-test")
prism_cert(v)=Model._extrude_nonew_prism_jacobian_certify(v,"nonhex-test")

# Independent exact reference-gradient evaluation. This uses the rational
# five-node pyramid and six-node prism shape functions rather than production
# corner determinants, Bernstein coefficients or interval arithmetic.
function reference_jacobian(v,family,u,s,w)
    u,s,w=Rational{BigInt}.((u,s,w))
    gradients=if family===:pyramid
        q=1-w
        signs=((-1,-1),(1,-1),(1,1),(-1,1))
        base=ntuple(4) do i
            a,b=signs[i]
            (a*(q+b*s)/(4q),b*(q+a*u)/(4q),-1//4+a*b*u*s/(4q*q))
        end
        (base...,(0//1,0//1,1//1))
    else
        barycentric=(1-u-s,u,s)
        triangle=((-1,-1),(1,0),(0,1))
        ntuple(6) do i
            node=mod1(i,3)
            direction=i<=3 ? -1 : 1
            (triangle[node][1]*(1+direction*w)/2,
             triangle[node][2]*(1+direction*w)/2,
             direction*barycentric[node]/2)
        end
    end
    derivatives=ntuple(3) do axis
        ntuple(3) do coordinate
            sum(Rational{BigInt}(v[i][coordinate])*gradients[i][axis]
                for i in eachindex(v))
        end
    end
    a,b,c=derivatives
    return a[1]*(b[2]*c[3]-b[3]*c[2])-
           a[2]*(b[1]*c[3]-b[3]*c[1])+
           a[3]*(b[1]*c[2]-b[2]*c[1])
end

# In the prism the determinant is affine in the triangle coordinates. At
# each triangle vertex, recover its exact power-basis quadratic from the
# reference gradients and evaluate its endpoints and stationary minimum.
function exact_prism_positive(v)
    orientation=sign(reference_jacobian(v,:prism,0,0,-1))
    orientation==0 && return false
    for (u,s) in ((0,0),(1,0),(0,1))
        low=orientation*reference_jacobian(v,:prism,u,s,-1)
        middle=orientation*reference_jacobian(v,:prism,u,s,0)
        high=orientation*reference_jacobian(v,:prism,u,s,1)
        (low>0 && high>0) || return false
        a=(low+high)/2-middle
        b=(high-low)/2
        if a>0
            stationary=-b/(2a)
            if -1<stationary<1
                a*stationary^2+b*stationary+middle>0 || return false
            end
        end
    end
    return true
end

function exact_pyramid_positive(v)
    corners=((-1,-1),(1,-1),(1,1),(-1,1))
    values=ntuple(i->reference_jacobian(v,:pyramid,corners[i]...,0),4)
    orientation=sign(values[2]) # Same cyclic triangle as the first canonical tet.
    return orientation!=0 && all(value->value*orientation>0,values)
end

function accepts(certify,v)
    try
        certify(v)
        return true
    catch err
        err isa ArgumentError || rethrow()
        return false
    end
end

function batch_pyramid(v,n)
    for _ in 1:n;pyramid_cert(v);end
    return nothing
end

function batch_prism(v,n)
    for _ in 1:n;prism_cert(v);end
    return nothing
end

function warmed_bytes(batch::F,v::T,n::Int) where {F,T}
    batch(v,n)
    return @allocated batch(v,n)
end

@testset "Complete P1 pyramid and prism Jacobians" begin
    for v in (PYRAMID,reflect(PYRAMID))
        @test pyramid_cert(v)===nothing
        @test exact_pyramid_positive(v)
    end
    # Both canonical tetrahedra are positive, yet the rational reference map
    # has a negative determinant strictly inside the pyramid.
    @test all(tet_signed_volume(FOLDED_PYRAMID[a],FOLDED_PYRAMID[b],
        FOLDED_PYRAMID[c],FOLDED_PYRAMID[d])>0
        for (a,b,c,d) in Model._EXTRUDE_CELL_TETS[7])
    @test reference_jacobian(FOLDED_PYRAMID,:pyramid,9//20,9//20,1//2)<0
    @test_throws r"full P1 Pyr5 Jacobian" pyramid_cert(FOLDED_PYRAMID)

    # A negative middle Bernstein coefficient can still define a positive
    # quadratic. The certificate must accept the complete valid map.
    for v in (PRISM,TWISTED_PRISM,NEAR_SINGULAR_PRISM,
              reflect(PRISM),reflect(TWISTED_PRISM),reflect(NEAR_SINGULAR_PRISM))
        @test exact_prism_positive(v)
        @test prism_cert(v)===nothing
    end
    @test reference_jacobian(NEAR_SINGULAR_PRISM,:prism,1//3,1//3,0)>0
    @test reference_jacobian(FOLDED_PRISM,:prism,1//3,1//3,0)<0
    @test reference_jacobian(SINGULAR_PRISM,:prism,1//3,1//3,0)==0
    for v in (FOLDED_PRISM,SINGULAR_PRISM)
        @test !exact_prism_positive(v)
        @test_throws r"full P1 Pri6 Jacobian" prism_cert(v)
    end

    for (v,certify) in ((PYRAMID,pyramid_cert),(PRISM,prism_cert))
        @test_throws r"finite coordinates" certify(Base.setindex(v,(NaN,0.,0.),1))
        @test_throws r"Jacobian" certify(ntuple(i->(v[i][1],v[i][2],0.),length(v)))
    end

    rng=MersenneTwister(0x1516)
    decisions=Set{Bool}()
    for _ in 1:80
        v=ntuple(i->ntuple(k->PRISM[i][k]+1.5*(rand(rng)-.5),3),6)
        exact=exact_prism_positive(v)
        push!(decisions,exact)
        @test accepts(prism_cert,v)==exact
    end
    @test decisions==Set((false,true))
    for _ in 1:40
        v=ntuple(i->ntuple(k->PYRAMID[i][k]+1.5*(rand(rng)-.5),3),5)
        @test accepts(pyramid_cert,v)==exact_pyramid_positive(v)
    end

    for (v,batch) in ((PYRAMID,batch_pyramid),(PRISM,batch_prism),
                       (TWISTED_PRISM,batch_prism)),n in (1000,2000,4000)
        @test warmed_bytes(batch,v,n)==0
    end
    for name in names(Model;all=true)
        (startswith(String(name),"_extrude_nonew_prism_jac") ||
         name===:_extrude_nonew_pyramid_jacobian_certify) || continue
        f=getfield(Model,name)
        f isa Function || continue
        for method in methods(f)
            @test !occursin("Core.Box",sprint(show,Base.uncompressed_ast(method)))
        end
    end
end

end # module
