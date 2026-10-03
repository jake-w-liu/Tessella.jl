# Gmsh 4.15.2 src/mesh/meshRefine.cpp: complete order-two support nodes
# followed by family-preserving linear subdivision. The pyramid transition
# produces four pyramids and eight tetrahedra, with matching boundary faces.
const _MIXED_REFINE_TEMPLATES=Dict(
    1=>((1,((1,3),(3,2))),),
    2=>((2,((1,4,6),(4,5,6),(4,2,5),(6,5,3))),),
    3=>((3,((1,5,9,8),(5,2,6,9),(9,6,3,7),(8,9,7,4))),),
    4=>((4,((1,5,7,8),(5,2,6,10),(7,6,3,9),(8,10,9,4),
            (5,7,8,10),(5,10,6,7),(7,8,10,9),(7,9,10,6))),),
    5=>((5,((1,9,21,10,11,22,27,23),(11,22,27,23,5,17,26,18),
            (9,2,12,21,22,13,24,27),(22,13,24,27,17,6,19,26),
            (10,21,14,4,23,27,25,16),(23,27,25,16,18,26,20,8),
            (21,12,3,14,27,24,15,25),(27,24,15,25,26,19,7,20))),),
    6=>((6,((1,7,8,9,16,17),(9,16,17,4,13,14),
            (7,2,10,16,11,18),(16,11,18,13,5,15),
            (8,10,3,17,18,12),(17,18,12,14,15,6),
            (10,8,7,18,17,16),(18,17,16,15,14,13))),),
    7=>((7,((1,6,14,7,8),(6,2,9,14,10),(14,9,3,11,12),(7,14,11,4,13))),
        (4,((8,10,13,5),(10,12,13,5),(10,13,12,14),(8,13,10,14),
            (8,10,6,14),(10,12,9,14),(13,11,12,14),(8,7,13,14)))),
    15=>((15,((1,),)),))

function _mixed_refine_limit(value,name)
    (value isa Integer && !(value isa Bool)) || throw(ArgumentError(
        "API.mesh.refine: $name must be an integer other than Bool"))
    0<=value<=typemax(Int32) || throw(ArgumentError(
        "API.mesh.refine: $name must lie in 0:$(typemax(Int32))"))
    return Int(value)
end

function _mixed_refine_plan(mesh,node_limit,cell_limit,caller)
    counts=Dict{Int,Int}()
    total=0
    for (msh,cells,_) in _cache_native_blocks(mesh)
        templates=get(_MIXED_REFINE_TEMPLATES,msh,nothing)
        templates===nothing && throw(ArgumentError(
            "$caller: refinement is unavailable for element type $msh"))
        for (target,children) in templates
            count=Refine._checked_mul(length(children),size(cells,2),"output cell")
            counts[target]=Refine._checked_add(get(counts,target,0),count,"output cell")
            total=Refine._checked_add(total,count,"output cell")
        end
    end
    total<=typemax(Int32) || throw(ArgumentError("$caller: output cell count exceeds Int32 indexing"))
    total<=cell_limit || throw(ArgumentError(
        "$caller: output requires $total cells, exceeding max_cells=$cell_limit"))
    for (msh,count) in counts
        entries=Refine._checked_mul(Elements.msh_num_nodes(msh),count,"connectivity entry")
        Refine._checked_mul(entries,sizeof(Int32),"connectivity byte")
        Refine._checked_mul(count,sizeof(Int32),"cell tag byte")
    end
    nnodes(mesh)<=node_limit || throw(ArgumentError(
        "$caller: primary nodes exceed max_nodes=$node_limit"))
    # Count the exact same weighted support identities as quadratic elevation;
    # topology, rather than coordinate equality, shares all edge/face nodes.
    supports=Set{Tuple}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        msh==15 && continue
        family=msh_spec(msh).family
        reference=Elements.lagrange_nodes(msh_type(family,2))
        for slot in size(cells,1)+1:size(reference,2)
            weights=MeshFunctionSpaces._first_order_values(Val(family),
                reference[1,slot],reference[2,slot],reference[3,slot],caller,slot)
            nonzero=findall(!iszero,collect(weights))
            for column in axes(cells,2)
                key=Tuple(sort!([(cells[j,column],Float64(weights[j])) for j in nonzero]))
                push!(supports,key)
                nnodes(mesh)+length(supports)<=node_limit || throw(ArgumentError(
                    "$caller: output requires more than max_nodes=$node_limit"))
            end
        end
    end
    nodes=Refine._checked_add(nnodes(mesh),length(supports),"output node")
    nodes<=typemax(Int32) || throw(ArgumentError("$caller: output node count exceeds Int32 indexing"))
    Refine._checked_mul(Refine._checked_mul(3,nodes,"coordinate entry"),sizeof(Float64),"coordinate byte")
    return counts,nodes
end

function _mixed_refine_certify(coords,cells,column,msh,orientation,caller)
    if msh in (4,5,6,7)
        for tet in Model._EXTRUDE_CELL_TETS[msh]
            p=ntuple(j->ntuple(axis->coords[axis,cells[tet[j],column]],3),4)
            exact=Model.orient3(p...)
            volume=Model.tet_signed_volume(p...)
            exact!=0 && isfinite(volume) && volume!=0 && sign(volume)==orientation ||
                throw(ArgumentError("$caller: refined volume cell is folded or degenerate"))
        end
    elseif msh in (2,3)
        for tri in (msh==2 ? ((1,2,3),) : ((1,2,3),(1,3,4)))
            pa,pb,pc=(ntuple(axis->coords[axis,cells[k,column]],3) for k in tri)
            u=ntuple(k->pb[k]-pa[k],3);v=ntuple(k->pc[k]-pa[k],3)
            normal=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
            all(isfinite,normal) && any(!iszero,normal) || throw(ArgumentError(
                "$caller: refined surface cell is degenerate"))
        end
    elseif msh==1
        first_node,last_node=cells[1,column],cells[2,column]
        any(axis->coords[axis,first_node]!=coords[axis,last_node],1:3) || throw(ArgumentError(
            "$caller: refined segment is below Float64 coordinate resolution"))
    end
    return nothing
end

function _mixed_refine_unclassified(mesh)
    owners=Dict{Int,Vector{Int32}}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        append!(get!(()->Int32[],owners,msh),zeros(Int32,size(cells,2)))
    end
    return _mixed_classification(mesh,(0,Int32(0)),Tuple{Int,Int32}[],
        fill((0,Int32(0)),nnodes(mesh)),Dict{Tuple{Int,Int32},Vector{Int32}}(),owners)
end

function _mixed_refine_locked!(m,cached::MixedMesh,max_nodes,max_cells)
    caller="API.mesh.refine"
    node_limit=_mixed_refine_limit(max_nodes,"max_nodes")
    cell_limit=_mixed_refine_limit(max_cells,"max_cells")
    diagnostic=validate(cached)
    diagnostic.ok || throw(ArgumentError("$caller: input mesh is invalid — "*join(diagnostic.messages,"; ")))
    stored_class=_cached_classification_locked(cached)
    _mixed_require_linear_cad(m,stored_class,caller)
    class=stored_class===nothing ? _mixed_refine_unclassified(cached) : stored_class
    linear,linear_class=_mixed_linear_cache(cached,class)
    # SetOrderN calls pruneMeshVertexAssociations: unreferenced source nodes
    # disappear, while explicit point elements keep their primary vertices.
    used=falses(nnodes(linear))
    for b in linear.blocks
        used[b.nodes].=true
    end
    linear,linear_class=_mixed_select_columns(linear,linear_class,
        [collect(axes(b.nodes,2)) for b in linear.blocks];node_order=findall(used))
    tagged=stored_class!==nothing && stored_class.public_tags!==nothing
    if tagged
        ordering_class=linear_class
        order=sortperm(linear.entity_data.external_node_tags;
            by=tag->(ordering_class.node_entities[ordering_class.public_tags.node_indices[tag]],tag),
            alg=Base.Sort.MergeSort)
        linear,linear_class=_mixed_select_columns(linear,linear_class,
            [collect(axes(b.nodes,2)) for b in linear.blocks];node_order=order)
    end
    counts,node_count=_mixed_refine_plan(linear,node_limit,cell_limit,caller)
    quadratic,quadclass=_mixed_quadratic_cache(linear,linear_class,caller)
    nnodes(quadratic)==node_count || error("$caller: internal support-node count mismatch")
    connectivity=Dict(msh=>Matrix{Int32}(undef,Elements.msh_num_nodes(msh),n) for (msh,n) in counts)
    tags=Dict(msh=>Vector{Int32}(undef,n) for (msh,n) in counts)
    owners=Dict(msh=>Vector{Int32}(undef,n) for (msh,n) in counts)
    cursors=Dict(msh=>0 for msh in keys(counts))
    external=tagged ? Dict(msh=>Vector{UInt64}(undef,n) for (msh,n) in counts) : nothing
    element_max=tagged ? max(ELEMENT_TAG_MAX[],maximum(stored_class.public_tags.element_tags;
        init=UInt64(0))) : UInt64(0)
    if tagged
        converted=sum(length(block.tags) for block in linear.blocks if block.msh!=15;init=0)
        element_max<=typemax(Int32)-2converted || throw(ArgumentError(
            "$caller: refinement order-conversion tags exceed Int32"))
        element_max+=UInt64(2converted)
    end
    catalog=_cache_catalog(quadratic,quadclass)
    for (bi,(input,_,_,support,parent_owners)) in enumerate(catalog)
        original=linear.blocks[bi]
        for parent in axes(support,2)
            orientation=if original.msh in (4,5,6,7)
                volume=Model._extrude_cell_signed_volume(linear.coords,original.nodes,parent,original.msh)
                isfinite(volume) && volume!=0 || throw(ArgumentError("$caller: parent volume cell is degenerate"))
                sign(volume)
            else
                0.0
            end
            for (target,children) in _MIXED_REFINE_TEMPLATES[original.msh],child in children
                column=cursors[target]+1
                for row in eachindex(child)
                    connectivity[target][row,column]=support[child[row],parent]
                end
                _mixed_refine_certify(quadratic.coords,connectivity[target],column,target,orientation,caller)
                tags[target][column]=original.tags[parent]
                owners[target][column]=parent_owners[parent]
                if tagged
                    if target==15
                        external[target][column]=linear.entity_data.external_element_tags[bi][parent]
                    else
                        element_max<typemax(Int32) || throw(ArgumentError(
                            "$caller: refined public element tags exceed Int32"))
                        element_max+=UInt64(1)
                        external[target][column]=element_max
                    end
                end
                cursors[target]=column
            end
        end
    end
    for msh in keys(counts)
        msh in (4,5,6,7) || continue
        # RefineMesh's final setAllVolumesPositive normalizes every volume
        # child, including a mesh the caller explicitly reversed beforehand.
        Model._extrude_make_positive!(quadratic.coords,connectivity[msh],msh)
        for column in axes(connectivity[msh],2)
            _mixed_refine_certify(quadratic.coords,connectivity[msh],column,msh,1.0,caller)
        end
    end
    blocks=[ElementBlock(msh,connectivity[msh],tags[msh]) for msh in sort!(collect(keys(counts))) if counts[msh]>0]
    if tagged
        assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
            NODE_TAG_MAX[],element_max,node_limit,cell_limit,caller)
        data=quadratic.entity_data
        for position in axes(quadratic.coords,2)
            _dim01_add_node!(assembly,data.external_node_tags[position],
                ntuple(axis->quadratic.coords[axis,position],3),
                quadclass.node_entities[position],data.node_parametric[position])
        end
        for block in blocks,column in axes(block.nodes,2)
            _dim01_add_cell!(assembly,block.msh,
                (Int(msh_dimension(block.msh)),owners[block.msh][column]),
                @view(block.nodes[:,column]),external[block.msh][column];
                physical=block.tags[column])
        end
        plan=_dim01_finish!(assembly,deepcopy(m),model_entities(m),
            Set{Tuple{Int,Int32}}(),falses(length(assembly.nodes)),
            !iszero(OPTIONS["Mesh.Renumber"]),false;
            save_all=!iszero(OPTIONS["Mesh.SaveAll"]),curve_parameter_order=true)
        diagnostic=validate(plan.mesh)
        diagnostic.ok || throw(ArgumentError("$caller: refinement produced an invalid mesh — "*
            join(diagnostic.messages,"; ")))
        _commit_tagged_plan_locked!(plan;authority=stored_class.public_tags.authority)
        return _copy_mesh(plan.mesh)
    end
    refined=MixedMesh(quadratic.coords,blocks;physical_names=cached.physical_names)
    diagnostic=validate(refined)
    diagnostic.ok || throw(ArgumentError("$caller: refinement produced an invalid mesh — "*join(diagnostic.messages,"; ")))
    cache=_copy_mesh(refined)
    new_class=stored_class===nothing ? nothing : _mixed_rebind_class(quadclass,cache;
        node_entities=copy(quadclass.node_entities),owners=owners)
    _replace_mesh_cache_locked!(cache,new_class)
    return refined
end
