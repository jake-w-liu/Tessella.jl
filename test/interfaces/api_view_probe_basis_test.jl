module ViewProbeBasisChecks
using Test,Tessella
const A=Tessella.API
function session(f)
 A.initialize()
 try
  A.option("General.Terminal",0)
  mktempdir(f)
 finally
  A.finalize()
 end
end
function mesh!(dir,name,msh,coords;connectivity=collect(eachindex(coords)),element=77,uid_offset=0)
 path=joinpath(dir,name*".msh");ids=UInt64[10i+1+uid_offset for i in eachindex(coords)]
 open(path,"w") do io
  println(io,raw"$MeshFormat");println(io,"2.2 0 8");println(io,raw"$EndMeshFormat")
  println(io,raw"$Nodes");println(io,length(ids))
  for (uid,p) in zip(ids,coords);println(io,uid," ",join(p," "));end
  println(io,raw"$EndNodes");println(io,raw"$Elements");println(io,1)
  println(io,element," ",msh," 2 0 1 ",join(ids[connectivity]," "));println(io,raw"$EndElements")
 end
 A.open(path)
 return A.model.get_current(),ids
end
@testset "Legacy and canonical probes retain high-order geometry, basis and ownership" begin
 session() do dir
  name,ids=mesh!(dir,"curved",8,[(0.,0.,0.),(1.,0.,0.),(.5,.5,0.)])
  for kind in ("NodeData","ElementNodeData")
   v=A.view.add(kind)
   A.view.add_model_data(v,0,name,kind,kind=="NodeData" ? ids : [77],
       kind=="NodeData" ? [[1.],[2.],[1.25]] : [[1.,2.,1.25]])
   @test A.view.probe(name,v,.25,.375,0.)≈[1.0625]
   values,distance=A.view.probe(v,.25,.375,0.)
   @test values≈[1.0625] && distance==0
   gradient,distance=A.view.probe(v,.25,.375,0.;gradient=true)
   @test gradient≈[.25,.25,0.] && distance==0
   @test A.view.probe(name,v,.5,.5,0.)≈[1.25]
   alias=A.view.add_alias(v);A.view.remove(v)
   A.model.add("unrelated target")
   @test A.view.probe(name,alias,.25,.375,0.)≈[1.0625]
   @test first(A.view.probe(alias,.25,.375,0.))≈[1.0625]
  end
 end
end
@testset "Canonical model steps, gaps, gradients and nearest-node distances" begin
 session() do dir
  name,ids=mesh!(dir,"straight",8,[(0.,0.,0.),(1.,0.,0.),(.5,0.,0.)])
  v=A.view.add("steps")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[1.],[2.],[1.25]])
  A.view.add_model_data(v,2,name,"NodeData",ids,[[4.],[6.],[4.5]])
  @test A.view.probe(v,.25,0.,0.)==([1.0625,0.,4.125],0.)
  @test A.view.probe(v,.25,0.,0.,2)==([4.125],0.)
  @test A.view.probe(v,.25,0.,0.;gradient=true)==([.5,0.,0.,0.,0.,0.,1.,0.,0.],0.)
  @test A.view.probe(v,1.5,0.,0.)==(Float64[],.5)
  @test A.view.probe(v,1.5,0.,0.;distanceMax=.49)==(Float64[],.5)
  @test A.view.probe(v,1.5,0.,0.;distanceMax=.5)==([2.,0.,6.],.5)
  @test A.view.probe(v,1.5,0.,0.;distanceMax=-1)==([2.,0.,6.],.5)
  @test A.view.probe(v,.25,0.,0.;numComp=3)==(Float64[],-1.)
  @test A.view.probe(v,.25,0.,0.;dim=2)==(Float64[],.25)
  @test_throws ArgumentError A.view.probe(v,.25,0.,0.;step=3)
  @test_throws ArgumentError A.view.probe(v,NaN,0.,0.)
 end
end
@testset "Legacy widths and allocated initial gaps retain their contracts" begin
 session() do dir
  name,ids=mesh!(dir,"components",1,[(0.,0.,0.),(1.,0.,0.)])
  vector2=A.view.add("legacy two components")
  A.view.add_model_data(vector2,0,name,"NodeData",ids,[[1.,2.],[2.,4.]])
  @test A.view.probe(name,vector2,.25,0.,0.)≈[1.25,2.5]
  @test A.view.probe(vector2,.25,0.,0.)==(Float64[],.25)
  @test A.view.probe(vector2,.25,0.,0.;numComp=2)==(Float64[],-1.)
  gap=A.view.add("allocated first gaps")
  A.view.add_model_data(gap,2,name,"NodeData",ids,[[1.],[2.]])
  @test A.view.probe(gap,.25,0.,0.)==([0.,0.,1.25],0.)
  @test A.view.probe(gap,.25,0.,0.;step=0)==([0.],0.)
  @test A.view.probe(gap,.25,0.,0.;gradient=true)==([0.,0.,0.,0.,0.,0.,1.,0.,0.],0.)
 end
end
@testset "List probes preserve components, gradients and ordered element selection" begin
 session() do dir
  v=A.view.add("vector plus scalar")
  A.view.add_list_data(v,"VL",1,[0.,1.,0.,0.,0.,0., 2.,3.,4.,2.,3.,4.])
  A.view.add_list_data(v,"SL",1,[0.,1.,0.,0.,0.,0.,1.,1.])
  @test A.view.probe(v,.25,0.,0.)==([1.],0.)
  @test A.view.probe(v,.25,0.,0.;numComp=3)==([2.,3.,4.],0.)
  @test A.view.probe(v,.25,0.,0.;numComp=9)==(Float64[],-1.)
  @test A.view.probe(v,.25,0.,0.;gradient=true)==([0.,0.,0.],0.)
  q=A.view.add("coincident triangles")
  A.view.add_list_data(q,"ST",2,[0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,1.,1.,
      0.,0.,1.,0.,1.,0.,0.,0.,0.,2.,2.,2.])
  @test A.view.probe(q,.2,.2,0.)==([2.],0.)
  @test A.view.probe(q,.2,.2,0.;xElemCoord=[0.,1.,0.],yElemCoord=[0.,0.,1.],zElemCoord=zeros(3))==([1.],0.)
  @test_throws ArgumentError A.view.probe(q,.2,.2,0.;xElemCoord=[0.],yElemCoord=[],zElemCoord=[])
 end
end
@testset "Mixed list classes select exact matches before global nearest-node fallback" begin
 session() do dir
  v=A.view.add("separated classes")
  A.view.add_list_data(v,"SL",1,[0.,1.,0.,0.,0.,0.,1.,1.])
  A.view.add_list_data(v,"VL",1,[10.,11.,0.,0.,0.,0.,2.,3.,4.,2.,3.,4.])
  @test A.view.probe(v,10.25,0.,0.)==([2.,3.,4.],0.)
  @test A.view.probe(v,10.25,0.,0.;distanceMax=100)==([2.,3.,4.],0.)
  @test A.view.probe(v,10.25,0.,0.;numComp=1)==(Float64[],-1.)
  @test A.view.probe(v,10.25,0.,0.;numComp=1,distanceMax=100)==(Float64[],.25)
  @test A.view.probe(v,9.,0.,0.;distanceMax=100)==([2.,3.,4.],1.)
  @test A.view.probe(v,9.,0.,0.;numComp=3,distanceMax=100,gradient=true)==([0.,0.,0.,0.,0.,0.,0.,0.,0.],1.)
  # Positive unsupported widths retain the documented component-filter rule.
  @test A.view.probe(v,.25,0.,0.;numComp=2)==(Float64[],-1.)
  A.plugin.set_number("Probe","View",v);A.plugin.set_number("Probe","X",10.25)
  A.plugin.set_number("Probe","Y",0.);A.plugin.set_number("Probe","Z",0.)
  A.plugin.run("Probe")
  types,data=A.view.get_list_data(last(A.view.get_tags()))
  @test types==["VP"]
  @test only(data)≈[10.25,0.,0.,2.,3.,4.]
 end
end
@testset "Probe gradients use the physical Jacobian for surfaces and vector volumes" begin
 session() do dir
  name,ids=mesh!(dir,"tilted",2,[(0.,0.,0.),(1.,0.,1.),(0.,1.,1.)])
  v=A.view.add("tangential x")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[0.],[1.],[0.]])
  @test first(A.view.probe(v,.2,.3,.5))≈[.2]
  @test first(A.view.probe(v,.2,.3,.5;gradient=true))≈[2/3,-1/3,1/3]
  coords=[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
          (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)]
  name,ids=mesh!(dir,"vector_hex",5,coords)
  v=A.view.add("affine vector")
  A.view.add_model_data(v,0,name,"NodeData",ids,
      [[p[1]+2p[2],3p[1]+4p[3],5p[2]+6p[3]] for p in coords])
  @test first(A.view.probe(v,.2,.3,.4))≈[.8,2.2,3.9]
  @test first(A.view.probe(v,.2,.3,.4;gradient=true))≈[1.,2.,0.,3.,0.,4.,0.,5.,6.]
 end
end
@testset "Pyramid gradients retain the physical map including its apex" begin
 session() do dir
  coords=[(-1.,-1.,0.),(1.,-1.,0.),(1.,1.,0.),(-1.,1.,0.),(0.,0.,1.)]
  name,ids=mesh!(dir,"pyramid",7,coords)
  v=A.view.add("linear pyramid field")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[p[1]+2p[2]+3p[3]] for p in coords])
  @test first(A.view.probe(v,.2,.1,.4))≈[1.6]
  @test first(A.view.probe(v,.2,.1,.4;gradient=true))≈[1.,2.,3.]
  @test first(A.view.probe(v,0.,0.,1.))≈[3.]
  @test first(A.view.probe(v,0.,0.,1.;gradient=true))≈[1.,2.,3.]
 end
end
@testset "Custom list value interpolation keeps its independent geometry basis" begin
 session() do dir
  v=A.view.add("quadratic on Line2 geometry")
  A.view.set_interpolation_matrices(v,"Line",3,
      [0.,-.5,.5, 0.,.5,.5, 1.,0.,-1.],
      [0.,0.,0., 1.,0.,0., 2.,0.,0.])
  A.view.add_list_data(v,"SL",1,[0.,1.,0.,0.,0.,0.,1.,2.,1.25])
  @test first(A.view.probe(v,.25,0.,0.))≈[1.0625]
  @test first(A.view.probe(v,.25,0.,0.;gradient=true))≈[.5,0.,0.]
 end
end
@testset "Custom geometry control points do not fabricate exact domain matches" begin
 session() do dir
  v=A.view.add("scaled custom geometry")
  exponents=[0.,0.,0.,1.,0.,0.]
  A.view.set_interpolation_matrices(v,"Line",2,[.5,-.5,.5,.5],exponents,
      2,[.05,-.05,.05,.05],exponents)
  A.view.add_list_data(v,"SL",1,[0.,1.,0.,0.,0.,0.,1.,2.])
  @test first(A.view.probe(v,.025,0.,0.))≈[1.25]
  @test first(A.view.probe(v,.025,0.,0.;gradient=true))≈[10.,0.,0.]
  # x=1 is a stored geometry control point outside the actual [0,0.1] map.
  @test A.view.probe(v,1.,0.,0.)==(Float64[],0.)
 end
end
@testset "Moving steps use material coordinates for all steps and their own physical frame individually" begin
 session() do dir
  name,ids=mesh!(dir,"firstframe",1,[(0.,0.,0.),(1.,0.,0.)])
  v=A.view.add("moving frames")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[1.],[2.]])
  A.model.add("another frame")
  name2,ids2=mesh!(dir,"secondframe",1,[(0.,0.,0.),(2.,0.,0.)])
  A.view.add_model_data(v,1,name2,"NodeData",ids2,[[11.],[12.]])
  @test first(A.view.probe(v,.25,0.,0.))≈[1.25,11.25]
  @test first(A.view.probe(v,.25,0.,0.;step=1))≈[11.125]
  @test first(A.view.probe(v,.25,0.,0.;gradient=true))≈[1.,0.,0.,1.,0.,0.]
  @test first(A.view.probe(v,.25,0.,0.;step=1,gradient=true))≈[.5,0.,0.]
 end
end
@testset "All-step material rows follow actual first-frame UIDs and populated history" begin
 session() do dir
  first_name,ids=mesh!(dir,"materialA",1,[(0.,0.,0.),(1.,0.,0.)])
  moving=A.view.add("reordered nodes")
  A.view.add_model_data(moving,0,first_name,"NodeData",ids,[[1.],[2.]])
  A.model.add("materialB")
  second_name,ids2=mesh!(dir,"materialB",1,[(10.,0.,0.),(11.,0.,0.)];connectivity=[2,1])
  A.view.add_model_data(moving,1,second_name,"NodeData",ids2,[[11.],[12.]])
  @test first(A.view.probe(moving,.25,0.,0.))≈[1.25,11.25]
  @test first(A.view.probe(moving,10.25,0.,0.;step=1))≈[11.25]
  @test first(A.view.probe(moving,.25,0.,0.;gradient=true))≈[1.,0.,0.,1.,0.,0.]
  A.plugin.set_number("Probe","View",moving)
  for (key,value) in (("X",.25),("Y",0.),("Z",0.));A.plugin.set_number("Probe",key,value);end
  A.plugin.run("Probe")
  @test only(last(A.view.get_list_data(last(A.view.get_tags()))))≈[.25,0.,0.,1.25,11.25]
  partial=A.view.add("partial after gap")
  A.view.add_model_data(partial,0,first_name,"NodeData",ids,[[1.],[2.]])
  A.view.add_model_data(partial,2,first_name,"NodeData",[first(ids)],[[11.]])
  @test first(A.view.probe(partial,.25,0.,0.))≈[1.25,0.,8.75]
  @test first(A.view.probe(partial,.25,0.,0.;step=2))≈[8.25]
  @test first(A.view.probe(partial,.25,0.,0.;gradient=true))≈[1.,0.,0.,0.,0.,0.,-9.,0.,0.]
  @test first(A.view.probe(partial,.25,0.,0.;step=2,gradient=true))≈[-11.,0.,0.]
  A.model.add("materialC")
  third_name,ids3=mesh!(dir,"materialC",1,[(20.,0.,0.),(21.,0.,0.)];element=88,uid_offset=100)
  for kind in ("NodeData","ElementData")
   v=A.view.add("different identities "*kind)
   A.view.add_model_data(v,0,first_name,kind,kind=="NodeData" ? ids : [77],
       kind=="NodeData" ? [[1.],[2.]] : [[1.]])
   A.view.add_model_data(v,1,third_name,kind,kind=="NodeData" ? ids3 : [88],
       kind=="NodeData" ? [[11.],[12.]] : [[11.]])
   @test first(A.view.probe(v,.25,0.,0.))≈(kind=="NodeData" ? [1.25,1.25] : [1.,1.])
   @test first(A.view.probe(v,20.25,0.,0.;step=1))≈(kind=="NodeData" ? [11.25] : [11.])
  end
 end
end
@testset "Element-node rows use their actual variable widths" begin
 session() do dir
  path=joinpath(dir,"mixed_orders.msh")
  write(path,raw"""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
5
11 0 0 0
21 1 0 0
31 2 0 0
41 3 0 0
51 2.5 0 0
$EndNodes
$Elements
2
77 1 2 0 1 11 21
88 8 2 0 1 31 41 51
$EndElements
""")
  A.open(path);name=A.model.get_current()
  v=A.view.add("ragged rows")
  A.view.add_model_data(v,0,name,"ElementNodeData",[77,88],[[1.,2.],[5.,10.,7.25]])
  @test A.view.probe(name,v,.25,0.,0.)≈[1.25]
  @test A.view.probe(name,v,2.25,0.,0.)≈[6.0625]
  @test first(A.view.probe(v,2.25,0.,0.;gradient=true))≈[4.5,0.,0.]
 end
end
@testset "Unbounded nearest-node search retains ordering when distances overflow" begin
 session() do dir
  v=A.view.add("large translated line")
  A.view.add_list_data(v,"SL",1,[-1e308,-9e307,0.,0.,0.,0.,1.,2.])
  values,distance=A.view.probe(v,1e308,0.,0.;distanceMax=-1)
  @test values≈[2.]
  @test isinf(distance) && distance>0
 end
end
@testset "Canonical model probes retain zero-valued unpopulated cells and nearest nodes" begin
 session() do dir
  path=joinpath(dir,"partial.msh")
  write(path,raw"""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
11 0 0 0
21 1 0 0
31 10 0 0
41 11 0 0
$EndNodes
$Elements
2
77 1 2 0 1 11 21
88 1 2 0 1 31 41
$EndElements
""")
  A.open(path);name=A.model.get_current()
  for kind in ("NodeData","ElementData","ElementNodeData")
   v=A.view.add(kind)
   A.view.add_model_data(v,0,name,kind,kind=="NodeData" ? [11] : [77],
       kind=="ElementNodeData" ? [[1.,2.]] : [[1.]])
   @test A.view.probe(v,10.25,0.,0.)==([0.],0.)
   @test A.view.probe(v,12.,0.,0.;distanceMax=-1)==([0.],1.)
   if kind=="NodeData"
    @test A.view.probe(name,v,.25,0.,0.)≈[.75]
    @test A.view.probe(name,v,10.25,0.,0.)==[0.]
   else
    @test_throws ArgumentError A.view.probe(name,v,10.25,0.,0.)
   end
  end
 end
end
@testset "Plugin Probe uses first-frame material coordinates for all steps" begin
 session() do dir
  name,ids=mesh!(dir,"origin",2,[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.)])
  v=A.view.add("moving triangles")
  A.view.add_model_data(v,0,name,"NodeData",ids,fill([2.],3))
  A.model.add("translated frame")
  name2,ids2=mesh!(dir,"translated",2,[(10.,0.,0.),(11.,0.,0.),(10.,1.,0.)])
  A.view.add_model_data(v,1,name2,"NodeData",ids2,fill([7.],3))
  A.plugin.set_number("Probe","View",v)
  for (key,value) in (("X",.25),("Y",.25),("Z",0.))
   A.plugin.set_number("Probe",key,value)
  end
  A.plugin.run("Probe")
  types,data=A.view.get_list_data(last(A.view.get_tags()))
  @test types==["SP"]
  @test only(data)≈[.25,.25,0.,2.,7.]
  @test A.view.probe(v,.25,.25,0.;step=1)[1]==Float64[]
  empty=A.view.add("empty")
  @test A.view.probe(empty,0.,0.,0.)==(Float64[],-1.)
  @test A.view.probe(empty,0.,0.,0.;step=0)==(Float64[],-1.)
 end
end
@testset "High-order probes preserve large constant values and zero derivatives" begin
 session() do dir
  name,ids=mesh!(dir,"constant",8,[(0.,0.,0.),(1.,0.,0.),(.5,0.,0.)])
  v=A.view.add("large constant")
  A.view.add_model_data(v,0,name,"NodeData",ids,fill([1e308],3))
  @test first(A.view.probe(v,1.,0.,0.))==[1e308]
  @test first(A.view.probe(v,1.,0.,0.;gradient=true))==[0.,0.,0.]
 end
end
@testset "Probe values and physical gradients retain a valid thin surface" begin
 session() do dir
  h=1e-200
  name,ids=mesh!(dir,"thin",2,[(0.,0.,0.),(1.,0.,0.),(0.,h,0.)])
  v=A.view.add("thin linear field")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[0.],[1.],[2.]])
  @test A.view.probe(name,v,.25,.25h,0.)≈[.75]
  @test first(A.view.probe(v,.25,.25h,0.))≈[.75]
  @test first(A.view.probe(v,.25,.25h,0.;gradient=true))≈[1.,2e200,0.]
 end
end
@testset "Finite values and gradients do not require a representable full inverse" begin
 session() do dir
  h=1e-310
  name,ids=mesh!(dir,"subnormal",2,[(0.,0.,0.),(1.,0.,0.),(0.,h,0.)])
  v=A.view.add("finite x derivative")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[0.],[1.],[0.]])
  @test A.view.probe(name,v,.25,.25h,0.)≈[.25]
  @test first(A.view.probe(v,.25,.25h,0.;gradient=true))≈[1.,0.,0.]
 end
end
@testset "Extreme value and geometry scales retain finite nonconstant gradients" begin
 for magnitude in (1e308,1e-300)
  session() do dir
   name,ids=mesh!(dir,"scaled",8,[(0.,0.,0.),(magnitude,0.,0.),(.5magnitude,0.,0.)])
   v=A.view.add("quadratic peak")
   A.view.add_model_data(v,0,name,"NodeData",ids,[[0.],[0.],[magnitude]])
   @test first(A.view.probe(v,0.,0.,0.))==[0.]
   @test first(A.view.probe(v,0.,0.,0.;gradient=true))≈[4.,0.,0.]
  end
 end
 session() do dir
  magnitude=1e-300
  name,ids=mesh!(dir,"dynamic",8,[(0.,0.,0.),(magnitude,0.,0.),(.5magnitude,0.,0.)])
  v=A.view.add("independent derivative scales")
  A.view.add_model_data(v,0,name,"NodeData",ids,[[0.],[magnitude],[1e308]])
  @test first(A.view.probe(v,.5magnitude,0.,0.;gradient=true))≈[1.,0.,0.]
 end
end
end
