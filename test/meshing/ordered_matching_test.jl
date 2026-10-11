module OrderedMatchingTests
using Test, Tessella, Random
const O=Tessella.OrderedMatching

# Independent subset recurrence, including an unmatched first vertex. The
# objective compares cardinality before an Int128 cost sum.
function optimum(costs,mask::UInt64,memo=Dict{UInt64,Tuple{Int,Int128}}())
    mask==0 && return (0,Int128(0))
    haskey(memo,mask) && return memo[mask]
    vertex=trailing_zeros(mask)+1
    rest=mask & ~(UInt64(1)<<(vertex-1))
    best=optimum(costs,rest,memo)
    for other in vertex+1:64
        bit=UInt64(1)<<(other-1)
        rest & bit==0 && continue
        cost=get(costs,(vertex,other),nothing)
        cost===nothing && continue
        pairs,value=optimum(costs,rest & ~bit,memo)
        candidate=(pairs+1,value+cost)
        if candidate[1]>best[1] || (candidate[1]==best[1] && candidate[2]<best[2])
            best=candidate
        end
    end
    return memo[mask]=best
end

function check(n,edges)
    costs=Dict{Tuple{Int,Int},Int128}()
    for (a,b,cost) in edges
        a==b && continue
        key=minmax(a,b)
        costs[key]=min(get(costs,key,Int128(cost)),Int128(cost))
    end
    pairs,value=optimum(costs,(UInt64(1)<<n)-UInt64(1))
    result=O.ordered_perfect_matching(n,edges)
    if 2pairs!=n
        @test result===nothing
    else
        @test result!==nothing
        result===nothing && return
        @test all(result.mate[v]!=v && 1<=result.mate[v]<=n &&
            result.mate[result.mate[v]]==v for v in 1:n)
        @test sum(Int128(edges[i][3]) for i in result.selected;init=Int128(0))==value
        @test sort!([v for i in result.selected for v in edges[i][1:2]])==collect(1:n)
        @test issorted([edges[i][2] for i in result.selected])
    end
    return result
end

@testset "Ordered perfect matching: every weighted graph through five vertices" begin
    for n in 0:5
        pairs=[(a,b) for a in 1:n for b in a+1:n]
        for code in 0:3^length(pairs)-1
            digits=code;edges=Tuple{Int,Int,Int64}[]
            for (a,b) in pairs
                digit=digits%3;digits÷=3
                digit==0 || push!(edges,(a,b,digit==1 ? -7 : 3))
            end
            check(n,edges)
        end
    end
end

@testset "Ordered perfect matching: every six-vertex topology and weighted profile" begin
    pairs=[(a,b) for a in 1:6 for b in a+1:6]
    for code in UInt32(0):(UInt32(1)<<length(pairs))-UInt32(1),profile in 1:2
        edges=Tuple{Int,Int,Int64}[]
        for (index,(a,b)) in enumerate(pairs)
            code & (UInt32(1)<<(index-1))==0 && continue
            push!(edges,(a,b,profile==1 ? 0 : mod(17index+Int(code),13)-6))
        end
        check(6,edges)
    end
end

@testset "Directed parallel identity, deterministic ties and checked input" begin
    rng=MersenneTwister(610091)
    for n in 7:12,sample in 1:60
        edges=Tuple{Int,Int,Int64}[(a,b,rand(rng,-20:20))
            for a in 1:n for b in a+1:n if rand(rng)<.45]
        check(n,edges)
    end
    edges=[(1,2,3),(2,1,3),(3,4,0),(4,3,0),(1,1,-100)]
    result=check(4,edges)
    @test result.selected==[2,4]
    @test O.ordered_perfect_matching(4,edges)==result
    @test edges==[(1,2,3),(2,1,3),(3,4,0),(4,3,0),(1,1,-100)]
    @test O.ordered_perfect_matching(0,Tuple{Int,Int,Int}[])==(mate=Int[],selected=Int[])
    @test O.ordered_perfect_matching(2,[(1,1,0)])===nothing
    @test O.ordered_perfect_matching(3,[(1,2,0)])===nothing
    @test_throws ArgumentError O.ordered_perfect_matching(-1,Tuple{Int,Int,Int}[])
    @test_throws ArgumentError O.ordered_perfect_matching(2,[(0,2,1)])
    @test_throws ArgumentError O.ordered_perfect_matching(2,[(1,2,1.)])
    @test_throws ArgumentError O.ordered_perfect_matching(2,[(1,2,typemax(Int64))])
    @test_throws ArgumentError O.ordered_perfect_matching(2,[(true,2,1)])
    @test_throws ArgumentError O.ordered_perfect_matching(2,[(1,2,true)])
end
@testset "Literal large integer costs retain the exact perfect optimum" begin
    scale=Int64(10)^16
    for sign in (-1,1),reversed in (false,true),diagonals in (false,true)
        value=sign*scale
        edges=[(1,2,value),(2,3,-value),(3,4,value),(4,1,-value)]
        diagonals && append!(edges,[(1,3,value),(2,4,-value)])
        reversed && reverse!(edges)
        check(4,edges)
    end
end
end # module OrderedMatchingTests
