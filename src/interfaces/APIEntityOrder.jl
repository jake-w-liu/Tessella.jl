# The order in which Gmsh 4.15.2 walks the elements of a model: entities by
# dimension and tag (`GModel::getEntities`), and inside an entity the element
# families in the order `GEntity::getMeshElement` concatenates them. Edge and
# face numbering (`createEdges`, `createFaces`) and the element octree follow it.

@inline function _entity_family_rank(family::Symbol)
    family===:tri && return 1
    family===:qua && return 2
    family===:tet && return 1
    family===:hex && return 2
    family===:pri && return 3
    family===:pyr && return 4
    return 5
end

# The cache catalog and one `(dimension, owner, family rank, block, column)`
# key per element, sorted into Gmsh's order. `selected`, when given, lists the
# entities to walk in the caller's own order instead of every entity.
function _entity_element_order(mesh::Union{Mesh,MixedMesh},class;selected=nothing)
    catalog=_cache_catalog(mesh,class)
    position=selected===nothing ? nothing :
        Dict((Int(dim),Int32(tag))=>index for (index,(dim,tag)) in enumerate(selected))
    keys=Tuple{Int,Int32,Int,Int,Int}[]
    for (block,(msh,dim,_,cells,owners)) in enumerate(catalog)
        family=msh_spec(Int(msh)).family
        family===:trih && continue
        rank=_entity_family_rank(family)
        classified=length(owners)==size(cells,2)
        for column in axes(cells,2)
            owner=classified ? owners[column] : Int32(0)
            if position===nothing
                push!(keys,(dim,owner,rank,block,column))
            else
                index=get(position,(dim,owner),0)
                index==0 || push!(keys,(index,owner,rank,block,column))
            end
        end
    end
    sort!(keys)
    return catalog,keys
end

# Apply `visit(block, columns)` to maximal runs of consecutive keys that stay
# in one block, so the typed per-block kernels do the element loop.
function _entity_order_runs(visit::F,keys) where F
    first=1
    while first<=length(keys)
        block=keys[first][4];last=first
        while last<length(keys) && keys[last+1][4]==block
            last+=1
        end
        visit(block,Int[keys[index][5] for index in first:last])
        first=last+1
    end
    return nothing
end

# `gmsh::model::mesh::createEdges`: every edge of every element, numbered at
# its first use in entity order.
function _mesh_edge_topology_in_entity_order(mesh::Union{Mesh,MixedMesh},class,topology;
                                             selected=nothing)
    catalog,keys=_entity_element_order(mesh,class;selected)
    replacement=MeshEntityTopology._edge_topology_copy(topology,nnodes(mesh))
    _entity_order_runs(keys) do block,columns
        msh,_,_,cells,_=catalog[block]
        MeshEntityTopology._append_generated_edges!(replacement,cells,
            MeshEntityTopology._simplex_edge_patterns(Int(msh)),Int(msh),columns)
    end
    return replacement
end

# `gmsh::model::mesh::createFaces`: every face in `MElement::getFace` order,
# triangles before quadrangles inside a prism or pyramid.
function _mesh_face_topology_in_entity_order(mesh::Union{Mesh,MixedMesh},class,topology;
                                             selected=nothing)
    catalog,keys=_entity_element_order(mesh,class;selected)
    replacement=MeshEntityTopology._face_topology_copy(topology,nnodes(mesh))
    _entity_order_runs(keys) do block,columns
        msh,dim,_,cells,_=catalog[block]
        dim>=2 || return nothing
        patterns=(MeshEntityTopology._simplex_face_patterns(Int(msh),3)...,
                  MeshEntityTopology._simplex_face_patterns(Int(msh),4)...)
        MeshEntityTopology._append_generated_faces!(replacement,cells,patterns,Int(msh),columns)
    end
    return replacement
end
