# Generation runs the same unset-owner transfer rule as public selector4.
# Its subsequent deleteUnusedVertices stage removes unused vertices owned by
# this Face; public mesh.recombine does not run that deletion stage. Retained
# Point/Curve vertices and all consumed constructor clocks remain actual.
function _gmsh_native_surface_bipartite(m,tag,product,world,caller)
    mesh=product.mesh
    cells=if mesh isa Mesh
        [collect(@view(mesh.tris[:,column])) for column in axes(mesh.tris,2)]
    else
        all(block->block.msh in (2,3),mesh.blocks) || throw(ArgumentError(
            "$caller: experimental Surface producer must contain first-order Face cells"))
        [collect(@view(block.nodes[:,column])) for block in mesh.blocks for column in axes(block.nodes,2)]
    end
    isempty(cells) && throw(ArgumentError(
        "$caller: experimental Surface producer has no imported Face cells"))
    events=_primary_bipartite_face_events(cells,tag,caller)
    owner=(2,Int32(tag))
    rows=findall(entity->entity!=owner,product.node_owners)
    mapping=zeros(Int32,length(product.node_tags));mapping[rows]=Int32.(eachindex(rows))
    birth_rows=Int32[mapping[row] for row in product.birth_rows if mapping[row]!=0]
    result=MixedMesh(copy(mesh.coords[:,rows]),MixedElementBlock[])
    output=_ModelSurfaceProduction(result,product.node_tags[rows],product.node_owners[rows],
        product.node_parameters[rows],product.retained_nodes[rows],product.node_high_water,
        UInt64[],product.element_high_water,birth_rows,nothing,
        product.source_node_updates,events)
    return output,world[:,rows]
end
