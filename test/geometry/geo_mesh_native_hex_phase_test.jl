module GeoMeshNativeHexPhaseTests
using Test
using Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_transfinite_hex_phase_primary_v1.jl"))
@testset "native transfinite Hex corner and cell phase" begin
 for native in NATIVE_TRANSFINITE_HEX_PHASE_PRIMARY_V1
  execution=mktemp() do path,io
   write(io,replace(native.source,";"=>";\n"));close(io);execute_geo(path;mesh_dim=3)
  end
  authority=execution.model.meshing.mesh_authority;actual=authority.mesh;data=actual.entity_data
  expected=Dict((c.owner[1],c.owner[2],uid)=>(c.msh,Tuple(nodes[:,column]))
   for c in native.cells for nodes in (reshape(c.nodes,msh_num_nodes(c.msh),:),)
   for (column,uid) in enumerate(c.tags))
  cells=Dict((msh_dimension(block.msh),data.block_entities[bi][column],uid)=>
    (block.msh,Tuple(data.external_node_tags[block.nodes[:,column]]))
   for (bi,block) in enumerate(actual.blocks) for (column,uid) in enumerate(data.external_element_tags[bi]))
  @test cells==expected
  @test Dict(zip(data.external_node_tags,data.node_entities))==Dict(node.uid=>node.owner for node in native.nodes)
  # Every captured carrier cell keeps the same actual primary support set,
  # even where the native Region reference frame changes its corner phase.
  catalog=Dict(uid=>Tuple(sort(data.external_node_tags[block.nodes[:,column]]))
   for (bi,block) in enumerate(actual.blocks) for (column,uid) in enumerate(data.external_element_tags[bi]))
  carrier=authority.carrier
  for (block,uids) in zip(carrier.blocks,authority.carrier_cells), (column,uid) in enumerate(uids)
   rows=authority.carrier_rows[block.nodes[:,column]]
   @test Tuple(sort(data.external_node_tags[rows]))==catalog[uid]
  end
 end
end
end
