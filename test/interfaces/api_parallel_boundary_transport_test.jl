module ParallelBoundaryTransportChecks
using Test,Tessella
A=Tessella.API;M=Tessella.Model
mesh=Tessella.Mesh([0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1.];tets=reshape(Int32[1,2,3,4],4,1))
constraints=(outer_faces=[(7,UInt64(21),(Int32(1),Int32(2),Int32(3)))],embedded_faces=[(8,UInt64(19),(Int32(1),Int32(3),Int32(4)))],edges=[(5,UInt64(22),(Int32(1),Int32(2)))])
nodes=Dict(Int32(2)=>(2,Int32(7)));params=Dict(Int32(2)=>[.25,.5])
surfaces=Dict(7=>[(2,Int32[1,2,3],UInt64(21))]);curves=Dict(5=>[(1,Int32[1,2],UInt64(22))])
replaced=Dict((2,Int32(7))=>Set{Tuple}([(Int32(1),Int32(2),Int32(3))]))
added=Dict((2,Int32(7))=>[(2,Int32[1,2,3],UInt64(21),Int32(9))])
wb=M._ModelVolumeBoundaryWriteback(constraints,nodes,params,surfaces,curves,replaced,added)
p=M._ModelVolumeGeometryProduction(mesh,UInt64[2,11,12,13],UInt64(13),UInt64(0),UInt64[23],UInt64(23),wb)
q=A._parallel_volume_remap_production(p,UInt64(100),UInt64(200),UInt64(10),UInt64(20))
@testset "Volume boundary receipt transport" begin
 @test A._parallel_volume_remap_supported(p)
 @test q.boundary_writeback.constraints.outer_faces[1][2]==221
 @test q.boundary_writeback.constraints.embedded_faces[1][2]==19
 @test q.boundary_writeback.constraints.edges[1][2]==222
 @test q.boundary_writeback.surfaces[7][1][3]==221
 @test q.boundary_writeback.curves[5][1][3]==222
 @test q.boundary_writeback.added[(2,Int32(7))][1][3]==221
 @test q.boundary_node_tags==UInt64[2,111,112,113]
 @test q.element_tags==UInt64[223]
 @test wb.constraints.outer_faces[1][2]==21
 @test wb.surfaces[7][1][3]==21
 @test q.boundary_writeback.nodes===wb.nodes
 @test q.boundary_writeback.node_params===wb.node_params
 @test q.boundary_writeback.replaced===wb.replaced
 @test q.boundary_writeback.added[(2,Int32(7))][1][4]==9
 @test A._parallel_volume_remap_boundary(wb,UInt64(0),UInt64(20))===wb
 A.initialize()
 try
  dirty=A._ParallelDirty();snapshot=copy(dirty)
  A._parallel_volume_mark_dirty!(dirty,A.CURRENT[],Dict{Int,Vector{Float64}}();boundary_writeback=wb)
  @test A._parallel_dirty_changed(dirty,snapshot,(2,Int32(7)))
  @test A._parallel_dirty_changed(dirty,snapshot,(1,Int32(5)))
  @test !A._parallel_dirty_changed(dirty,snapshot,(3,Int32(99)))
  snapshot=copy(dirty)
  A._parallel_volume_mark_dirty!(dirty,A.CURRENT[],Dict{Int,Vector{Float64}}();boundary_writeback=wb)
  @test A._parallel_dirty_changed(dirty,snapshot,(2,Int32(7)))
 finally
  A.finalize()
 end
end

function public_run(workers)
 A.initialize()
 try
  A.option("General.Terminal",0);A.option("General.NumThreads",workers)
  A.open(joinpath(@__DIR__,"..","fixtures","geo_parallel_recombine_rtree.geo"))
  A.mesh.generate(2)
  m=A.CURRENT[];input=m.meshing.mesh_authority
  if workers>1
   prepared=A._parallel_volume_prepare(m,[1,102],input,"WB regression";
    size_field=nothing,scope=nothing,extruded=nothing,hybrid_batch=nothing,staged=false)
   @test prepared!==nothing && all(!isnothing,prepared[1])
   data=input.mesh.entity_data
   hn=maximum(data.external_node_tags;init=UInt64(0))
   hc=maximum((tag for block in data.external_element_tags for tag in block);init=UInt64(0))
   working=M._model_mesh_transaction_copy(m)
   options=M._model_default_mesh1d_options(working,"WB regression")
   production=M.mesh_model_volume(working,1;_mesh_input=input,
    _mesh1d_options=options,_return_volume_production=true,
    _mesh_node_high_water=hn,_mesh_element_high_water=hc,
    _extrude_working_model=true,_boundary_regrade=false)
   @test production.boundary_writeback!==nothing
   @test A._parallel_volume_remap_supported(production)
  end
  A._PARALLEL_ADOPTED[]=0;A._PARALLEL_SEED_ADOPTED[]=0
  A.mesh.generate(3)
  mesh=m.meshing.mesh_authority.mesh
  bytes=mktempdir() do dir
   file=joinpath(dir,"out.msh");Tessella.Elements.write_mixed_msh(file,mesh;version=4.1);read(file)
  end
  return bytes,Tessella.validate(mesh).ok,A._PARALLEL_ADOPTED[],A._PARALLEL_SEED_ADOPTED[]
 finally
  A.finalize()
 end
end
if Threads.nthreads()>=2
 @testset "RTREE recombined volumes transport real boundary products" begin
  serial=public_run(1);parallel=public_run(4)
  @test serial[2] && parallel[2]
  @test serial[1]==parallel[1]
  @test serial[3]==0 && serial[4]==0
  @test parallel[3]>=2
  # This counter increments only for actual adopted Volume receipts that
  # consumed the RNG stream, proving both Regions ran through transport.
  @test parallel[4]==2
 end
end
end
