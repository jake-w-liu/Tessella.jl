module ViewNearestResourceChecks
using Test, Tessella, Random
const A=Tessella.API
@testset "Nearest view nodes preserve filters and extreme distances" begin
 rng=MersenneTwister(718231)
 for run in 1:32
  records=NamedTuple[(kind=:point,dim=i%3,coords=randn(rng,3,4)) for i in 1:8]
  point=Tuple(randn(rng,3))
  for dim in (-1,0,1,2)
   candidates=[(hypot((point.-Tuple(r.coords[:,node]))...),i,node)
    for (i,r) in enumerate(records) for node in axes(r.coords,2)
    if dim<0 || r.dim==dim]
   expected=minimum(candidates)
   match,distance=A._view_probe_find(records,point,dim,nothing;exact_match=false)
   @test distance==expected[1]
   @test match[1]===records[expected[2]]
   @test match[2]==Tuple(records[expected[2]].coords[:,expected[3]])
  end
 end
 records=NamedTuple[(kind=:point,dim=0,coords=reshape([Float64(i),0.,0.],3,1)) for i in 1:3]
 @test A._view_probe_find(records,(0.,0.,0.),1,nothing;exact_match=false)==(nothing,Inf)
 @test A._view_probe_find(records,(0.,0.,0.),-1,records[2].coords;exact_match=false)==((records[2],(2.,0.,0.)),2.)
 @test A._view_probe_find(records,(0.,0.,0.),-1,nothing;exact_match=false,nearest_match=false)==(nothing,Inf)
 huge=floatmax(Float64)
 extreme=NamedTuple[(kind=:point,dim=0,coords=reshape([-huge,-huge,0.],3,1)),
  (kind=:point,dim=0,coords=reshape([-huge,0.,0.],3,1))]
 match,distance=A._view_probe_find(extreme,(huge,huge,0.),-1,nothing;exact_match=false)
 @test match[1]===extreme[2]
 @test isinf(distance)
 # The actual record producer uses an abstract NamedTuple vector. Warm the
 # production query before checking its allocation bound under optimization.
 if Base.JLOptions().compile_enabled!=3
  many=NamedTuple[(kind=:point,dim=0,coords=fill(Float64(129-i),3,8)) for i in 1:128]
  query()=A._view_probe_find(many,(0.,0.,0.),-1,nothing;exact_match=false)
  query();allocated=@allocated query()
  @test allocated<4096
 end
end
end
