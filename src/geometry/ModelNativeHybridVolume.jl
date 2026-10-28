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
        constructor.node_uids,constructor.coordinates,caller;embedded_point=constructor.embedded_point,embedded_curve=constructor.embedded_curve,
        embedded_surface=constructor.embedded_surface,embedded_constraints=constructor.embedded_constraints)
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
    embedded_point=_model_native_hybrid_embedded_point_input(m,tag,input,caller)
    embedded_point===nothing || push!(used,embedded_point.record.uid)
    embedded_curve=_model_native_hybrid_embedded_curve_input(m,tag,input,caller)
    embedded_curve===nothing || union!(used,(record.uid for record in embedded_curve.records))
    embedded_surface=_model_native_hybrid_embedded_surface_input(m,tag,input,caller)
    embedded_surface===nothing || union!(used,(record.uid for record in embedded_surface.records))
    embedded_constraints=_model_native_hybrid_embedded_constraints_input(m,tag,input,caller)
    embedded_constraints===nothing || union!(used,(record.uid for record in embedded_constraints.records))
    # Only actual selected lower cell incidence enters the native Source.
    # Unreferenced lower rows remain global.
    ids=sort!(collect(used))
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    all(uid->haskey(rows,uid),ids) || throw(ArgumentError(
        "$caller: hybrid Source cell UID has no current authoritative vertex"))
    coords=actual.coords[:,Int32[rows[uid] for uid in ids]]
    constructor=_model_native_hybrid_constructor_input(input,boundaries,ids,coords,caller;embedded_point,embedded_curve,embedded_surface,embedded_constraints)
    _model_native_hybrid_point_containment_certificate(constructor,caller)
    _model_native_hybrid_curve_containment_certificate(constructor,caller)
    _model_native_hybrid_surface_containment_certificate(constructor,caller)
    _model_native_hybrid_constraints_containment_certificate(constructor,caller)
    return constructor
end

# Admit each proper-rigid boundary pair from copyMesh's actual UID events.
# Coordinates verify those events; neither geometry nor row order selects a
# pair. Independent directions (an x/y/z-periodic box) certify per pair —
# upstream applies the same per-`GFace` `setMeshMaster` copy under each
# relation and keeps per-face correspondence records, so pairs only need to
# agree where corner nodes coincide bitwise.
function _model_native_hybrid_periodic_boundary_certificate(m,tag,constructor,pairs,caller)
    constraints=_model_volume_periodic_surface_constraints(m,tag,caller)
    isempty(constraints) && return nothing
    for constraint in constraints
        affine=constraint.affine
        # The copy map is a certified proper rigid transform: finite applied
        # entries, an orthonormal linear part and determinant +1. The stored
        # fourth row is verbatim record bytes like upstream's — the applied 3×4
        # map is what this certifies. `bound` is a Float64 rounding allowance
        # on the declared coefficients — distinct from the geometric `atol`
        # that the node transport check below applies.
        affine!==nothing && all(isfinite,affine[i] for i in 1:12) ||
            throw(ArgumentError(
                "$caller: hybrid periodic boundary requires a finite rigid Surface transform"))
        linear=(affine[1],affine[2],affine[3],affine[5],affine[6],affine[7],
                affine[9],affine[10],affine[11])
        rigid_bound=1e-12*max(1.0,maximum(abs,linear))
        for column_i in 1:3,column_j in column_i:3
            gram=linear[column_i]*linear[column_j]+linear[3+column_i]*linear[3+column_j]+
                linear[6+column_i]*linear[6+column_j]
            abs(gram-(column_i==column_j ? 1.0 : 0.0))<=rigid_bound || throw(ArgumentError(
                "$caller: hybrid periodic boundary requires an orthonormal rigid transform"))
        end
        determinant=linear[1]*(linear[5]*linear[9]-linear[6]*linear[8])-
            linear[2]*(linear[4]*linear[9]-linear[6]*linear[7])+
            linear[3]*(linear[4]*linear[8]-linear[5]*linear[7])
        abs(determinant-1.0)<=rigid_bound || throw(ArgumentError(
            "$caller: hybrid periodic boundary requires a determinant +1 rigid transform"))
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
    end
    return nothing
end

function _model_native_hybrid_volume_production(m,tag,constructor,options,births,caller;
                                                size_field=nothing,periodic_surface_pairs=nothing,
                                                optimize_pyramids::Float64=0.0)
    _model_native_hybrid_volume_scope(m,tag,constructor,caller)
    isfinite(optimize_pyramids) || throw(ArgumentError(
        "$caller: native hybrid Pyramid optimization value must be finite"))
    m.meshing.algorithm_3d==1 &&
        (isempty(get(m.embeds,(3,tag),NTuple{2,Int}[])) ||
         constructor.embedded_point!==nothing || constructor.embedded_curve!==nothing ||
         constructor.embedded_surface!==nothing || constructor.embedded_constraints!==nothing) &&
        !get(m.meshing.reverse,(3,tag),false) && get(m.meshing.smoothing,(3,tag),0)==0 || throw(ArgumentError(
        "$caller: hybrid Volume requires its owned ordinary Delaunay constructor without additional planner attributes"))
    _model_native_hybrid_periodic_boundary_certificate(m,tag,constructor,periodic_surface_pairs,caller)
    # Mesh.OptimizeNetgen on a hybrid region is an upstream no-op (the Netgen
    # optimizer skips hybrid meshes); `hybrid=true` below encodes that skip.
    node_clock=births===nothing ? maximum(constructor.input.mesh.entity_data.external_node_tags;init=UInt64(0)) : births.node_high_water
    element_clock=births===nothing ? maximum((uid for tags in constructor.input.mesh.entity_data.external_element_tags for uid in tags);init=UInt64(0)) : births.element_high_water
    work=_ModelNativeVolumeBirths(node_clock,element_clock)
    plan=_model_native_hybrid_quad_plan(constructor.node_uids,constructor.coordinates,
        constructor.quads,node_clock,caller;
        normal_factor=optimize_pyramids<0.0 ? -optimize_pyramids : 0.0)
    optimize_pyramids<0.0 && _model_native_hybrid_normal_apex_certificate(plan,constructor,caller)
    work.source=Mesh(plan.source_coordinates);work.node_tags=copy(plan.node_uids)
    work.node_high_water=plan.node_high_water
    fan=copy(plan.fan_facets)
    for triangle in constructor.triangles;fan=hcat(fan,collect(triangle.rows));end
    surface=Mesh(plan.coordinates;tris=fan)
    order=sortperm(plan.node_uids);inverse=invperm(order)
    kernel=Ref{ModelNativeInitialTet.Kernel}()
    _model_native_initial_delaunay_fill(surface,order,options,caller;births=work,kernel_output=kernel)
    facets=_model_native_hybrid_recovery_facets(plan,constructor,order,caller)
    curve_segments=vcat(_model_native_hybrid_curve_segments(constructor,plan.node_uids,caller),
        _model_native_hybrid_constraints_segments(constructor,plan.node_uids,caller))
    kernel_segments=NTuple{2,Int32}[ntuple(slot->Int32(inverse[edge[slot]]),2) for edge in curve_segments]
    sheet=_model_native_hybrid_surface_constraints(constructor,plan.node_uids,caller)
    csheet=_model_native_hybrid_constraints_facets(constructor,plan.node_uids,caller)
    all_facets=vcat(sheet.facets,csheet.facets)
    kernel_facets=NTuple{3,Int32}[ntuple(slot->Int32(inverse[face[slot]]),3) for face in all_facets]
    protected_edges=Set(minmax(edge...) for edge in vcat(curve_segments,sheet.edges,csheet.edges))
    protected_faces=Set(Tuple(sort!(collect(face))) for face in all_facets)
    interior=ModelNativeHybridInterior.recover(kernel[],facets;embedded_segments=kernel_segments,
        embedded_facets=kernel_facets)
    tets=reduce(hcat,Int32[order[row] for row in cell] for cell in interior.cells)
    mesh=Mesh(plan.coordinates;tets=tets)
    _model_native_volume_cell_births!(work,mesh)
    surface_witnesses=_model_native_hybrid_surface_witnesses(constructor)
    triangles=_model_native_hybrid_global_triangles(m,constructor.input,caller;surface_witnesses)
    seeds=_model_native_hybrid_global_size_seeds(m,constructor.input,caller;surface_witnesses)
    insertion=ModelNativeHybridInsertion.refine(mesh,work.node_tags,work.element_tags,
        triangles,_ModelAlgorithmBackgroundSize(m,3,tag,options,size_field,caller),caller;
        minimize=_resolved_extend_from_boundary(m,3,tag)==2,
        extend=_resolved_extend_from_boundary(m,3,tag)>0,
        node_high_water=work.node_high_water,element_high_water=work.element_high_water,seeds,
        protected_edges,protected_faces,
        metric=_volume_metric_insertion_field(m,size_field,options,caller),
        metric_entity=(3,tag))
    smoothed=ModelNativeHybridInsertion.smooth(insertion,length(plan.node_uids),
        m.meshing.native_smoothing,caller)
    # Native negative factors select the normal-offset births and skip the
    # Pyramid smoother; Tet insertion/optimization keep these shell rows fixed.
    relocated=optimize_pyramids<0.0 ? smoothed.mesh :
        (optimize_pyramids==1.0 ? ModelNativeHybridWinslow : ModelNativeHybridRelocate).relocate(
            plan,smoothed.mesh,caller;birth_order=constructor.pyramid_order,
            node_uids=insertion.node_numbers)
    pyramids=_model_native_hybrid_pyramid_production(plan,relocated.coords,
        insertion.element_high_water,caller;birth_order=constructor.pyramid_order)
    work.cell_source=relocated;work.element_tags=insertion.element_numbers
    work.node_tags=copy(insertion.node_numbers);work.node_high_water=insertion.node_high_water
    work.element_high_water=pyramids.element_high_water
    final_sheet=_model_native_hybrid_surface_constraints(constructor,work.node_tags,caller)
    final_csheet=_model_native_hybrid_constraints_facets(constructor,work.node_tags,caller)
    final_edges=vcat(_model_native_hybrid_curve_segments(constructor,work.node_tags,caller),
        _model_native_hybrid_constraints_segments(constructor,work.node_tags,caller),
        final_sheet.edges,final_csheet.edges)
    mesh=_model_native_optimize_volume(relocated,options,work,constructor.input,
        vcat(final_sheet.facets,final_csheet.facets),final_edges,caller;hybrid=true)
    blocks=ElementBlock[ElementBlock(4,copy(mesh.tets)),ElementBlock(7,pyramids.nodes)]
    mixed=MixedMesh(mesh.coords,blocks)
    optimize_pyramids<0.0 && _model_native_hybrid_normal_apex_fixed_certificate(
        plan,mixed,work.node_tags,caller)
    validation=validate(mixed)
    validation.ok || throw(ArgumentError(
        "$caller: hybrid Volume produces invalid typed cells: "*
        join(validation.messages,"; ")))
    _model_native_hybrid_volume_certificate(mixed,work.node_tags,constructor,caller)
    _model_native_hybrid_embedded_point_certificate(mixed,work.node_tags,constructor,caller)
    _model_native_hybrid_embedded_curve_certificate(mixed,work.node_tags,constructor,caller)
    _model_native_hybrid_embedded_surface_certificate(mixed,work.node_tags,constructor,caller)
    _model_native_hybrid_embedded_constraints_certificate(mixed,work.node_tags,constructor,caller)
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
    _model_native_hybrid_global_triangles(m,constructor.input,caller;
        surface_witnesses=_model_native_hybrid_surface_witnesses(constructor))
    return nothing
end

# Sizing sees every current GFace Tri constructor, including a declared
# unrelated Face sharing one selected UID. Its foreign rows contribute edge
# lengths without entering the Region's actual Source or boundary registry.
function _model_native_hybrid_global_triangles(m,input,caller;surface_witnesses=())
    input isa _ModelMeshAuthority && input.mesh!==nothing && input.mesh.entity_data!==nothing ||
        throw(ArgumentError("$caller: global hybrid sizing lacks current actual lower authority"))
    _model_native_hybrid_global_size_seeds(m,input,caller;surface_witnesses)
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
