struct _ModelNativeHybridVolumeProduction
    mesh::MixedMesh
    boundary_node_tags::Vector{UInt64}
    node_high_water::UInt64
    temporary_elements::UInt64
    element_tags::Vector{UInt64}
    element_high_water::UInt64
    input::_ModelMeshAuthority
    constructor::NamedTuple
    registration_rows::Vector{Int32}
end

# Successful insertion registers its MVertices before split-Quad recovery
# appends Pyramid tops in the restored Region Face/Quad constructor order.
# Preserve that presentation independently of physical birth UIDs and kernel
# coordinate rows, which must keep their arithmetic/optimizer layout.
function _model_native_hybrid_registration_rows(constructor,plan,insertion,node_uids,caller)
    length(insertion.node_numbers)>=length(plan.node_uids) &&
        insertion.node_numbers[1:length(plan.node_uids)]==plan.node_uids &&
        allunique(insertion.node_numbers) && allunique(node_uids) && all(!iszero,node_uids) ||
        throw(ArgumentError("$caller: hybrid registration lacks its actual insertion/Quad birth trace"))
    current=_model_native_hybrid_constructor_input(constructor.input,constructor.boundaries,
        constructor.node_uids,constructor.coordinates,caller)
    isequal(current,constructor) || throw(ArgumentError(
        "$caller: hybrid registration changes its restored Region constructor"))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(node_uids))
    source=Set(plan.source_uids)
    registered=UInt64[uid for uid in node_uids if uid in source]
    append!(registered,(uid for uid in insertion.node_numbers[length(plan.node_uids)+1:end]
        if haskey(positions,uid)))
    births=Dict((birth.owner,birth.element_uid)=>birth for birth in plan.births)
    length(births)==length(plan.births) || throw(ArgumentError(
        "$caller: hybrid registration repeats a recorded Quad birth"))
    for ordinal in constructor.pyramid_order
        quad=constructor.quads[ordinal];birth=get(births,(quad.owner,quad.uid),nothing)
        birth!==nothing && birth.base_uids==quad.nodes && haskey(positions,birth.apex_uid) ||
            throw(ArgumentError("$caller: hybrid registration loses an actual restored Pyramid top"))
        push!(registered,birth.apex_uid)
    end
    length(registered)==length(node_uids) && allunique(registered) && Set(registered)==Set(node_uids) ||
        throw(ArgumentError("$caller: hybrid registration omits or aliases a producer-owned row"))
    return Int32[positions[uid] for uid in registered]
end

function _model_native_hybrid_volume_input(m,tag,input,caller)
    input isa _ModelMeshAuthority && input.mesh!==nothing || return nothing
    boundaries=_model_volume_boundary_surfaces(m,tag,caller)
    isempty(boundaries) && return nothing
    actual=input.mesh;data=actual.entity_data
    data===nothing && return nothing
    quads=false;used=Set{UInt64}()
    for surface in boundaries
        product=get(input.products,(2,Int32(abs(surface))),nothing)
        product===nothing && return nothing
        for block in product.blocks
            block.msh in (2,3) || return nothing
            quads|=block.msh==3 && !isempty(block.nodes)
            for local_row in block.nodes
                push!(used,data.external_node_tags[product.rows[local_row]])
            end
        end
    end
    quads || return nothing
    # meshGRegion's initial MVertex set uses current Face cell incidence.
    # Unreferenced lower rows remain global records, outside this constructor.
    ids=sort!(collect(used))
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    all(uid->haskey(rows,uid),ids) || throw(ArgumentError(
        "$caller: hybrid Source cell UID has no current authoritative vertex"))
    coords=actual.coords[:,Int32[rows[uid] for uid in ids]]
    return _model_native_hybrid_constructor_input(input,boundaries,ids,coords,caller)
end

# Admit one translated boundary pair only from copyMesh's actual UID events.
# Coordinates verify those events; neither geometry nor row order selects a pair.
function _model_native_hybrid_periodic_boundary_certificate(m,tag,constructor,pairs,caller)
    constraints=_model_volume_periodic_surface_constraints(m,tag,caller)
    isempty(constraints) && return nothing
    length(constraints)==1 || throw(ArgumentError(
        "$caller: hybrid periodic boundary requires one translated Surface pair"))
    constraint=only(constraints);affine=constraint.affine
    affine!==nothing && all(index->affine[index]==(index in (1,6,11,16) ? 1.0 : 0.0),
        (1,2,3,5,6,7,9,10,11,13,14,15,16)) || throw(ArgumentError(
        "$caller: hybrid periodic boundary supports translations only"))
    slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
    event=pairs===nothing ? nothing : get(pairs,slave,nothing)
    event!==nothing && event[1]==master && !isempty(event[2]) || throw(ArgumentError(
        "$caller: hybrid periodic boundary lacks actual Surface copy UID pairs"))
    input=constructor.input
    source=_model_periodic_surface_source(m,constraint,input,caller)
    target=get(input.surface_constructor_inputs,(2,Int32(slave)),nothing)
    source!==nothing && target!==nothing || throw(ArgumentError(
        "$caller: hybrid periodic boundary lacks owned Surface constructors"))
    _model_surface_constructor_binding(_model_mesh_product(input,2,slave;caller),target,caller)
    bindings=event[2];sd=source.entity_data;td=target.entity_data
    mapping=Dict{UInt64,UInt64}()
    seen=Set{UInt64}()
    for (slave_uid,master_uid) in bindings
        !haskey(mapping,master_uid) && !(slave_uid in seen) || throw(ArgumentError(
            "$caller: hybrid periodic copy repeats an actual UID pair"))
        mapping[master_uid]=slave_uid;push!(seen,slave_uid)
    end
    Set(keys(mapping))==Set(sd.external_node_tags) && seen==Set(td.external_node_tags) ||
        throw(ArgumentError("$caller: hybrid periodic copy omits actual constructor rows"))
    target_rows=Dict(uid=>row for (row,uid) in enumerate(td.external_node_tags))
    point_map,curve_pairs,_=_model_periodic_surface_edge_map(m,slave,master,
        affine,constraint.atol,caller)
    curve_map=Dict(master_curve=>slave_curve for (slave_curve,master_curve) in curve_pairs)
    coefficients,translation=_periodic_affine_3x4(affine)
    scale=max(1.0,maximum(abs,source.coords;init=0.0),maximum(abs,target.coords;init=0.0))
    for (row,uid) in enumerate(sd.external_node_tags)
        target_row=target_rows[mapping[uid]];owner=sd.node_entities[row]
        expected=owner[1]==0 ? (0,Int32(get(point_map,Int(owner[2]),0))) :
            owner[1]==1 ? (1,Int32(get(curve_map,Int(owner[2]),0))) : (2,Int32(slave))
        (owner[1]<2 || owner==(2,Int32(master))) && td.node_entities[target_row]==expected ||
            throw(ArgumentError("$caller: hybrid periodic copy changes bound node ownership"))
        transformed=_model_affine_point(coefficients,translation,
            Tuple(source.coords[:,row]),caller,row)
        all(axis->abs(target.coords[axis,target_row]-transformed[axis])<=constraint.atol*scale,1:3) ||
            throw(ArgumentError("$caller: hybrid periodic copy changes its bound affine coordinates"))
    end
    mapped_cells=[(Int(block.msh),Tuple(mapping[sd.external_node_tags[row]] for row in block.nodes[:,column]))
        for block in source.blocks for column in axes(block.nodes,2)]
    target_cells=[(Int(block.msh),Tuple(td.external_node_tags[block.nodes[:,column]]))
        for block in target.blocks for column in axes(block.nodes,2)]
    mapped_cells==target_cells || throw(ArgumentError(
        "$caller: hybrid periodic copy changes ordered actual Surface cell incidence"))
    return nothing
end

function _model_native_hybrid_volume_production(m,tag,constructor,options,births,caller;
                                                size_field=nothing,periodic_surface_pairs=nothing)
    _model_native_hybrid_volume_scope(m,tag,constructor,caller)
    m.meshing.algorithm_3d==1 && isempty(get(m.embeds,(3,tag),NTuple{2,Int}[])) &&
        !get(m.meshing.reverse,(3,tag),false) && get(m.meshing.smoothing,(3,tag),0)==0 || throw(ArgumentError(
        "$caller: hybrid Volume requires its owned ordinary Delaunay constructor without additional planner attributes"))
    _model_native_hybrid_periodic_boundary_certificate(m,tag,constructor,periodic_surface_pairs,caller)
    options.optimize_netgen==0 || throw(ArgumentError(
        "$caller: native hybrid Netgen optimization is not implemented"))
    node_clock=births===nothing ? maximum(constructor.input.mesh.entity_data.external_node_tags;init=UInt64(0)) : births.node_high_water
    element_clock=births===nothing ? maximum((uid for tags in constructor.input.mesh.entity_data.external_element_tags for uid in tags);init=UInt64(0)) : births.element_high_water
    work=_ModelNativeVolumeBirths(node_clock,element_clock)
    plan=_model_native_hybrid_quad_plan(constructor.node_uids,constructor.coordinates,
        constructor.quads,node_clock,caller)
    work.source=Mesh(plan.source_coordinates);work.node_tags=copy(plan.node_uids)
    work.node_high_water=plan.node_high_water
    fan=copy(plan.fan_facets)
    for triangle in constructor.triangles;fan=hcat(fan,collect(triangle.rows));end
    surface=Mesh(plan.coordinates;tris=fan)
    order=sortperm(plan.node_uids);inverse=invperm(order)
    kernel=Ref{ModelNativeInitialTet.Kernel}()
    _model_native_initial_delaunay_fill(surface,order,options,caller;births=work,kernel_output=kernel)
    facets=_model_native_hybrid_recovery_facets(plan,constructor,order,caller)
    interior=ModelNativeHybridInterior.recover(kernel[],facets)
    tets=reduce(hcat,Int32[order[row] for row in cell] for cell in interior.cells)
    mesh=Mesh(plan.coordinates;tets=tets)
    _model_native_volume_cell_births!(work,mesh)
    triangles=_model_native_hybrid_global_triangles(m,constructor.input,caller)
    seeds=_model_native_hybrid_global_size_seeds(m,constructor.input,caller)
    insertion=ModelNativeHybridInsertion.refine(mesh,work.node_tags,work.element_tags,
        triangles,_ModelAlgorithmBackgroundSize(m,3,tag,options,size_field,caller),caller;
        minimize=_resolved_extend_from_boundary(m,3,tag)==2,
        extend=_resolved_extend_from_boundary(m,3,tag)>0,
        node_high_water=work.node_high_water,element_high_water=work.element_high_water,seeds)
    smoothed=ModelNativeHybridInsertion.smooth(insertion,length(plan.node_uids),
        m.meshing.native_smoothing,caller)
    relocated=ModelNativeHybridRelocate.relocate(plan,smoothed.mesh,caller;
        birth_order=constructor.pyramid_order,node_uids=insertion.node_numbers)
    pyramids=_model_native_hybrid_pyramid_production(plan,relocated.coords,
        insertion.element_high_water,caller;birth_order=constructor.pyramid_order)
    work.cell_source=relocated;work.element_tags=insertion.element_numbers
    work.node_tags=copy(insertion.node_numbers);work.node_high_water=insertion.node_high_water
    work.element_high_water=pyramids.element_high_water
    mesh=_model_native_optimize_volume(relocated,options,work,constructor.input,
        NTuple{3,Int32}[],NTuple{2,Int32}[],caller)
    blocks=ElementBlock[ElementBlock(4,copy(mesh.tets)),ElementBlock(7,pyramids.nodes)]
    mixed=MixedMesh(mesh.coords,blocks)
    validate(mixed).ok || throw(ArgumentError("$caller: hybrid Volume produces invalid typed cells"))
    _model_native_hybrid_volume_certificate(mixed,work.node_tags,constructor,caller)
    _model_native_hybrid_periodic_boundary_certificate(m,tag,constructor,periodic_surface_pairs,caller)
    geometry=_model_explicit_volume_geometry(m,tag,caller;_mesh_input=constructor.input,_mesh1d_options=options)
    _model_certify_explicit_volume_semantics(mixed,tag,geometry.expected_volume,
        geometry.comparison_scale,caller)
    return _ModelNativeHybridVolumeProduction(mixed,copy(work.node_tags),work.node_high_water,
        work.temporary_elements,vcat(work.element_tags,pyramids.element_tags),work.element_high_water,
        constructor.input,constructor,
        _model_native_hybrid_registration_rows(constructor,plan,insertion,work.node_tags,caller))
end

# Native boundary reconstruction creates the actual GFace Tri vectors first,
# in the tag-ordered Face constructor traversal, then the split-Quad MFace
# map. Combining both into one sorted list changes physical pool cell births.
function _model_native_hybrid_recovery_facets(plan,constructor,order,caller)
    length(order)==length(plan.node_uids) && allunique(order) &&
        all(row->0<row<=length(plan.node_uids),order) &&
        plan.source_uids==constructor.node_uids &&
        isequal(plan.source_coordinates,constructor.coordinates) || throw(ArgumentError(
        "$caller: hybrid facet constructor loses its owned Source and initial row permutation"))
    inverse=invperm(order)
    triangles=NTuple{3,Int32}[]
    for triangle in constructor.triangles
        Tuple(plan.node_uids[collect(triangle.rows)])==triangle.nodes || throw(ArgumentError(
            "$caller: hybrid actual Tri constructor loses its current UID incidence"))
        push!(triangles,Tuple(Int32.(inverse[collect(triangle.rows)])))
    end
    fans=NTuple{3,Int32}[Tuple(Int32.(inverse[plan.fan_facets[:,column]])) for column in axes(plan.fan_facets,2)]
    sort!(fans;by=face->Tuple(sort!(UInt64[plan.node_uids[order[row]] for row in face])))
    return vcat(triangles,fans)
end

function _model_native_hybrid_volume_scope(m,tag,constructor,caller)
    length(m.volumes)==1 && haskey(m.volumes,tag) || throw(ArgumentError(
        "$caller: native hybrid multi-Region constructor history is not implemented"))
    _model_native_hybrid_global_triangles(m,constructor.input,caller)
    return nothing
end

# Sizing sees every current GFace Tri constructor, including a declared
# unrelated Face sharing one selected UID. Its foreign rows contribute edge
# lengths without entering the Region's actual Source or boundary registry.
function _model_native_hybrid_global_triangles(m,input,caller)
    input isa _ModelMeshAuthority && input.mesh!==nothing && input.mesh.entity_data!==nothing ||
        throw(ArgumentError("$caller: global hybrid sizing lacks current actual lower authority"))
    _model_native_hybrid_global_size_seeds(m,input,caller)
    triangles=NamedTuple[];seen=Set{UInt64}()
    for owner in sort!(collect(keys(input.products)))
        owner[1]==2 || continue
        product=input.products[owner]
        any(block->!isempty(block.nodes),product.blocks) || continue
        haskey(m.surfaces,Int(owner[2])) || throw(ArgumentError(
            "$caller: global hybrid sizing names an undeclared Surface owner"))
        actual=_model_mesh_product(input,2,Int(owner[2]);caller)
        data=actual.entity_data
        data!==nothing && length(data.external_node_tags)==nnodes(actual) &&
            allunique(data.external_node_tags) && all(!iszero,data.external_node_tags) &&
            all(isfinite,actual.coords) || throw(ArgumentError(
            "$caller: global hybrid Surface lacks unique current UID/XYZ records"))
        for (block,labels,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)
            block.msh in (2,3) && length(labels)==length(owners)==size(block.nodes,2) ||
                throw(ArgumentError("$caller: global hybrid sizing requires actual linear Surface cells"))
            for (column,uid) in enumerate(labels)
                uid>0 && !(uid in seen) && owners[column]==owner[2] || throw(ArgumentError(
                    "$caller: global hybrid sizing repeats or changes current cell UID ownership"))
                push!(seen,uid)
                block.msh==2 || continue
                size(block.nodes,1)==3 && all(row->0<row<=nnodes(actual),block.nodes[:,column]) ||
                    throw(ArgumentError("$caller: global hybrid Tri has invalid actual incidence"))
                nodes=ntuple(slot->data.external_node_tags[block.nodes[slot,column]],3)
                allunique(nodes) || throw(ArgumentError("$caller: global hybrid Tri repeats a current UID"))
                coordinates=ntuple(slot->ntuple(axis->actual.coords[axis,block.nodes[slot,column]],3),3)
                push!(triangles,(owner=owner,uid=uid,nodes=nodes,coordinates=coordinates))
            end
        end
    end
    return triangles
end

function _model_native_hybrid_volume_certificate(mesh,node_uids,constructor,caller)
    constructor.input isa _ModelMeshAuthority && length(node_uids)==nnodes(mesh) &&
        allunique(node_uids) && all(!iszero,node_uids) || throw(ArgumentError(
        "$caller: hybrid Volume lacks complete operation-owned current UID rows"))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(node_uids))
    counts=Dict{Tuple,Int}()
    for block in mesh.blocks
        block.msh in (4,7) || throw(ArgumentError("$caller: hybrid Volume has a foreign typed cell family"))
        for column in axes(block.nodes,2),face in _VOLUME_CELL_FACES[Int(block.msh)]
            key=Tuple(sort!(Int32[block.nodes[slot,column] for slot in face]))
            counts[key]=get(counts,key,0)+1
        end
    end
    all(count->count<=2,values(counts)) || throw(ArgumentError(
        "$caller: hybrid Volume has non-manifold typed face incidence"))
    wanted=Set{Tuple}()
    for face in vcat(constructor.triangles,constructor.quads)
        all(uid->haskey(positions,uid),face.nodes) || throw(ArgumentError(
            "$caller: hybrid Volume omits an actual Surface constructor vertex"))
        key=Tuple(sort!(Int32[positions[uid] for uid in face.nodes]))
        key in wanted && throw(ArgumentError("$caller: hybrid Source repeats a typed shell panel"))
        push!(wanted,key)
    end
    boundary=Set(face for (face,count) in counts if count==1)
    boundary==wanted || throw(ArgumentError(
        "$caller: hybrid Volume exterior differs from the exact actual Tri/Quad Source shell"))
    edges=_volume_edge_set(mesh)
    actual=constructor.input.mesh;data=actual.entity_data
    for curve in unique(abs(curve) for surface in constructor.boundaries for curve in
            data.entities[(2,Int(abs(surface)))].boundaries)
        product=get(constructor.input.products,(1,Int32(curve)),nothing)
        product!==nothing || throw(ArgumentError("$caller: hybrid shell lacks actual Curve[$curve] cells"))
        for block in product.blocks
            block.msh==1 || throw(ArgumentError("$caller: hybrid shell requires actual linear Curve cells"))
            for column in axes(block.nodes,2)
                ids=UInt64[data.external_node_tags[product.rows[row]] for row in block.nodes[:,column]]
                all(uid->haskey(positions,uid),ids) || throw(ArgumentError(
                    "$caller: hybrid shell omits an actual Curve cell UID"))
                minmax((positions[uid] for uid in ids)...) in edges || throw(ArgumentError(
                    "$caller: hybrid shell omits actual Line constructor incidence"))
            end
        end
    end
    return nothing
end
