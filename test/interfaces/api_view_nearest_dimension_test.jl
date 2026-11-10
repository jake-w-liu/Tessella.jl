module ViewNearestDimensionChecks
using Test,Tessella
const A=Tessella.API
@testset "Nearest distance spans the view before dimension selection" begin
 A.initialize()
 try
  A.option("General.Terminal",0)
  v=A.view.add("global nearest")
  A.view.add_list_data(v,"SP",1,[0.,0.,0.,11.])
  A.view.add_list_data(v,"SL",1,[2.,3.,0.,0.,0.,0.,22.,33.])
  # Neither query is on a cell. The globally nearest node belongs to a
  # different dimension in half the cases: a farther selected cell cannot
  # replace that nearest node when applying the value search filter.
  for (x,owner,value) in ((1.8,1,22.),(.2,0,11.)),dim in (-1,0,1,2),radius in (0.,-1.,.5,2.)
   values,distance=A.view.probe(v,x,0.,0.;dim,distanceMax=radius)
   expected_distance=owner==1 ? 2.0-x : x
   @test distance==expected_distance
   available=(dim<0 || dim==owner) && radius!=0 &&
    (radius<0 || expected_distance<=radius)
   @test values==(available ? [value] : Float64[])
  end
  # Private search order remains explicit: exact searches rank element
  # families first, and later records win a same-family distance tie.
  coords=reshape([1.,0.,0.],3,1)
  records=NamedTuple[(kind=:point,dim=0,coords=copy(coords)),
   (kind=:line,dim=1,coords=copy(coords)),
   (kind=:line,dim=1,coords=copy(coords))]
  selected,_=A._view_probe_find(records,(0.,0.,0.),-1,nothing;
   exact_match=false,reverse_order=true)
  @test selected[1]===records[3]
  first_selected,_=A._view_probe_find(records,(0.,0.,0.),-1,nothing;
   exact_match=false)
  @test first_selected[1]===records[1]
  empty=A.view.add("empty nearest")
  @test A.view.probe(empty,1.,0.,0.;dim=2,distanceMax=-1.)==(Float64[],-1.)
  @test A.view.probe(v,1.8,0.,0.;dim=0,numComp=1)==(Float64[],-1.)
 finally
  A.finalize()
 end
end
end
