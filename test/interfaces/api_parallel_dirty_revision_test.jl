module ParallelDirtyRevisionChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model
@testset "Parallel receipts observe every entity revision" begin
 dirty=A._ParallelDirty();entity=(1,Int32(7));other=(2,Int32(9))
 @test isempty(dirty)
 push!(dirty,entity);snapshot=copy(dirty)
 @test entity in dirty
 @test !A._parallel_dirty_changed(dirty,snapshot,entity)
 push!(dirty,entity)
 @test A._parallel_dirty_changed(dirty,snapshot,entity)
 @test length(dirty)==1
 @test Set(dirty)==Set([entity])
 @test !A._parallel_dirty_changed(dirty,snapshot,other)
 phase=A._ParallelDirty();union!(phase,dirty)
 @test !A._parallel_dirty_changed(phase,dirty,entity)
 push!(phase,entity)
 @test A._parallel_dirty_changed(phase,dirty,entity)
 @test !A._parallel_dirty_changed(dirty,snapshot,other)
 A.initialize()
 try
  A.option("General.Terminal",0)
  p=A.model.add_point(0.,0.,0.;meshSize=.2)
  q=A.model.add_point(1.,0.,0.;meshSize=.2)
  curve=A.model.add_line(p,q);m=A.CURRENT[]
  options=M._model_default_mesh1d_options(m,"revision regression")
  work=A._ParallelCurveWork(curve,Set([(1,Int32(curve))]),nothing)
  for mode in (:surface,:volume)
   dirty=A._ParallelDirty();m.curve_params[curve]=[0.,1.]
   mark!(refs)=mode===:surface ? A._parallel_surface_mark_dirty!(dirty,m,nothing,[],refs) :
    A._parallel_volume_mark_dirty!(dirty,m,refs)
   refs=A._parallel_surface_param_refs(m)
   m.curve_params[curve]=[0.,.5,1.];mark!(refs);snapshot=copy(dirty)
   seed=options.random_state[]
   out=A._parallel_curve_attempt(M._model_mesh_transaction_copy(m),work,m,
    A._parallel_task_options(options,seed),seed,"revision regression",false)
   @test out.ok
   @test A._parallel_curve_fresh(work,out,dirty,snapshot,options,seed)
   # This second write was invisible to a set of dirty entity identities.
   refs=A._parallel_surface_param_refs(m)
   insert!(m.curve_params[curve],2,.25);mark!(refs)
   @test !A._parallel_curve_fresh(work,out,dirty,snapshot,options,seed)
   snapshot=copy(dirty);refs=A._parallel_surface_param_refs(m)
   delete!(m.curve_params,curve);mark!(refs)
   @test A._parallel_dirty_changed(dirty,snapshot,(1,Int32(curve)))
   @test !A._parallel_curve_fresh(work,out,dirty,snapshot,options,seed)
  end
 finally
  A.finalize()
 end
end
function surface_run(workers)
 A.initialize()
 try
  A.option("General.Terminal",0);A.option("General.NumThreads",workers)
  for offset in (0.,3.,6.)
   pts=[A.model.add_point(offset+x,y,0.;meshSize=.6) for (x,y) in ((0.,0.),(1.,0.),(0.,1.))]
   lines=[A.model.add_line(pts[i],pts[mod1(i+1,3)]) for i in 1:3]
   A.model.add_plane_surface([A.model.add_curve_loop(lines)])
  end
  A._PARALLEL_ADOPTED[]=0;A.mesh.generate(2)
  mesh=A.CURRENT[].meshing.mesh_authority.mesh
  bytes=mktempdir() do dir
   file=joinpath(dir,"surfaces.msh");Tessella.Elements.write_mixed_msh(file,mesh;version=4.1);read(file)
  end
  return (bytes,Tessella.validate(mesh).ok,A._PARALLEL_ADOPTED[])
 finally
  A.finalize()
 end
end
if Threads.nthreads()>=2
 @testset "Revision transport preserves actual parallel Curve and Surface meshes" begin
  serial=surface_run(1);parallel=surface_run(4)
  @test serial[2] && parallel[2]
  @test serial[1]==parallel[1]
  @test serial[3]==0
  # Nine Curves can contribute at most nine adoptions; this also proves
 # actual Surface adoption without requiring every eligible task to detach.
 @test parallel[3]>9
 end
end
end
