module WeightedMatchingTests
using Test, Tessella, Random
const W=Tessella.WeightedMatching

# Independent subset enumeration: remove the first remaining vertex either
# unmatched or with each remaining neighbor. Compare cardinality before cost.
function brute_optimum(n,costs,mask::UInt64,memo=Dict{UInt64,Tuple{Int,Int128}}())
 mask==0 && return (0,Int128(0))
 haskey(memo,mask) && return memo[mask]
 first=trailing_zeros(mask)+1;remaining=mask & ~(UInt64(1)<<(first-1))
 best=brute_optimum(n,costs,remaining,memo)
 for second in first+1:n
  bit=UInt64(1)<<(second-1)
  remaining & bit==0 && continue
  haskey(costs,(first,second)) || continue
  pairs,cost=brute_optimum(n,costs,remaining & ~bit,memo)
  candidate=(pairs+1,cost+costs[(first,second)])
  if candidate[1]>best[1] || (candidate[1]==best[1] && candidate[2]<best[2])
   best=candidate
  end
 end
 memo[mask]=best
 return best
end
function check_graph(n,edges)
 costs=Dict{Tuple{Int,Int},Int128}()
 for (a,b,cost) in edges
  a==b && continue
  edge=minmax(a,b);costs[edge]=min(get(costs,edge,Int128(cost)),Int128(cost))
 end
 mate=W.minimum_weight_matching(n,edges)
 @test all(mate[v]==0 || (mate[v]!=v && 1<=mate[v]<=n &&
   mate[mate[v]]==v && haskey(costs,minmax(v,mate[v]))) for v in 1:n)
 actual=(count(!iszero,mate)÷2,sum(costs[(v,mate[v])] for v in 1:n if mate[v]>v;init=Int128(0)))
 @test actual==brute_optimum(n,costs,(UInt64(1)<<n)-UInt64(1))
 return mate
end

@testset "All weighted general graphs through five vertices" begin
 for n in 0:5
  pairs=[(a,b) for a in 1:n for b in a+1:n]
  for code in 0:3^length(pairs)-1
   state=code;edges=Tuple{Int,Int,Int64}[]
   for (a,b) in pairs
    digit=state%3;state÷=3
    digit==0 || push!(edges,(a,b,digit==1 ? -7 : 3))
   end
   check_graph(n,edges)
  end
 end
end
@testset "Every six-vertex topology with distinct weighted profiles" begin
 pairs=[(a,b) for a in 1:6 for b in a+1:6]
 for code in UInt32(0):(UInt32(1)<<length(pairs))-UInt32(1),profile in 1:2
  edges=Tuple{Int,Int,Int64}[]
  for (index,(a,b)) in enumerate(pairs)
   code & (UInt32(1)<<(index-1))==0 && continue
   cost=profile==1 ? 0 : mod(17index+Int(code),13)-6
   push!(edges,(a,b,cost))
  end
  check_graph(6,edges)
 end
end
@testset "Larger odd cycles, nested contractions and parallel edges" begin
 rng=MersenneTwister(61009)
 for n in 7:12,sample in 1:100
  edges=Tuple{Int,Int,Int64}[(a,b,rand(rng,-20:20))
   for a in 1:n for b in a+1:n if rand(rng)<.45]
  check_graph(n,edges)
 end
 edges=[(1,2,5),(1,2,-3),(2,3,-8),(1,1,-100),(3,4,2)]
 @test check_graph(4,edges)==[2,1,4,3]
 @test edges==[(1,2,5),(1,2,-3),(2,3,-8),(1,1,-100),(3,4,2)]
 @test W.minimum_weight_matching(0,Tuple{Int,Int,Int}[])==Int[]
 @test W.minimum_weight_matching(5,Tuple{Int,Int,Int}[])==zeros(Int,5)
 @test_throws ArgumentError W.minimum_weight_matching(-1,Tuple{Int,Int,Int}[])
 @test_throws ArgumentError W.minimum_weight_matching(2,[(1,3,2)])
 @test_throws ArgumentError W.minimum_weight_matching(2,[(1,2,2.)])
 @test_throws ArgumentError W.minimum_weight_matching(2,[(1,2,typemax(Int64))])
 @test_throws ArgumentError W.minimum_weight_matching(2,[(1,2,true)])
 @test_throws ArgumentError W.minimum_weight_matching(2,[(true,2,3)])
end
@testset "Literal large integer costs retain the exact optimum" begin
 scale=Int64(10)^16
 for sign in (-1,1),reversed in (false,true),diagonals in (false,true)
  value=sign*scale
  edges=[(1,2,value),(2,3,-value),(3,4,value),(4,1,-value)]
  diagonals && append!(edges,[(1,3,value),(2,4,-value)])
  reversed && reverse!(edges)
  check_graph(4,edges)
 end
end
end # module WeightedMatchingTests
