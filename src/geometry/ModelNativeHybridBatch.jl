struct _ModelNativeHybridBatchProduction
    input::_ModelMeshAuthority
    regions::Vector{Int}
    groups::Vector{Vector{Int}}
    products::Dict{Int,_ModelNativeHybridVolumeProduction}
    node_high_water::UInt64
    element_high_water::UInt64
    random_state::UInt32
end

function _model_native_hybrid_groups(m,regions,input,caller)
    tags=Int.(collect(regions))
    !isempty(tags) && allunique(tags) && all(tag->haskey(m.volumes,tag),tags) ||
        throw(ArgumentError("$caller: hybrid batch lacks unique declared selected Regions"))
    sort!(tags);faces=Dict{Int,Vector{Int}}()
    for tag in tags
        boundary=abs.(_model_volume_boundary_surfaces(m,tag,caller))
        !isempty(boundary) && allunique(boundary) || throw(ArgumentError(
            "$caller: hybrid batch Region lacks a unique current declared Face closure"))
        for face in boundary;push!(get!(()->Int[],faces,face),tag);end
    end
    remaining=Set(tags);groups=Vector{Int}[]
    for seed in tags
        seed in remaining || continue
        delete!(remaining,seed);group=Int[seed];cursor=1
        while cursor<=length(group)
            tag=group[cursor];cursor+=1
            for face in abs.(_model_volume_boundary_surfaces(m,tag,caller)),neighbor in faces[face]
                if neighbor in remaining;delete!(remaining,neighbor);push!(group,neighbor);end
            end
        end
        push!(groups,sort!(group))
    end
    # Certify every group's lower input before any group consumes RNG or a
    # factory clock. Selection is from declared Face connectivity, never XYZ.
    receipts=[_model_native_hybrid_component_receipt(m,group,input,caller) for group in groups]
    return groups,receipts
end

# Disjoint periodic Regions consume the same actual Surface copy events as
# the single-Region producer. The complete lower constructor binding precedes
# RNG/factory use; coordinates only verify those explicit UID events.
function _model_native_hybrid_batch_periodic_certificate(m,receipt,pairs,optimize_pyramids,caller)
    for tag in receipt.regions
        constraints=_model_volume_periodic_surface_constraints(m,tag,caller)
        isempty(constraints) && continue
        length(receipt.regions)==1 && iszero(optimize_pyramids) || throw(ArgumentError(
            "$caller: periodic hybrid batches require separated Regions and default Pyramid optimization"))
        constructor=_model_native_hybrid_volume_input(m,tag,receipt.input,caller)
        constructor!==nothing || throw(ArgumentError(
            "$caller: periodic hybrid batch Region lacks its actual lower constructor"))
        _model_native_hybrid_periodic_boundary_certificate(m,tag,constructor,pairs,caller)
    end
    return nothing
end

function _model_native_hybrid_batch_production(m,regions,input,options,births,caller;
                                                size_field=nothing,optimize_pyramids::Float64=0.0,
                                                periodic_surface_pairs=nothing)
    isfinite(optimize_pyramids) || throw(ArgumentError("$caller: hybrid batch Pyramid setting must be finite"))
    options.optimize_netgen==0 || throw(ArgumentError("$caller: hybrid batch Netgen optimization is not implemented"))
    groups,receipts=_model_native_hybrid_groups(m,regions,input,caller)
    _model_native_hybrid_global_triangles(m,input,caller)
    _model_native_hybrid_global_size_seeds(m,input,caller)
    for receipt in receipts,tag in receipt.regions
        m.meshing.algorithm_3d==1 && isempty(get(m.embeds,(3,tag),NTuple{2,Int}[])) &&
            !get(m.meshing.reverse,(3,tag),false) && get(m.meshing.smoothing,(3,tag),0)==0 || throw(ArgumentError(
                "$caller: hybrid batch has unsupported selected Region planner attributes"))
    end
    for receipt in receipts
        _model_native_hybrid_batch_periodic_certificate(m,receipt,periodic_surface_pairs,optimize_pyramids,caller)
    end
    deferred=_ModelMesh1DOptions((name===:optimize ? 0 : getfield(options,name)
        for name in fieldnames(_ModelMesh1DOptions))...)
    node_clock=births.node_high_water;element_clock=births.element_high_water
    products=Dict{Int,_ModelNativeHybridVolumeProduction}()
    plans=Dict{Int,_ModelNativeHybridQuadPlan}()
    for receipt in receipts
        current=_ModelNativeVolumeBirths(node_clock,element_clock)
        group=_model_native_hybrid_connected_production(m,receipt,deferred,current,caller;
            size_field,optimize_pyramids,periodic_surface_pairs)
        for (tag,product) in group.products
            haskey(products,tag) && throw(ArgumentError("$caller: hybrid batch repeats a selected Region product"))
            products[tag]=product
            plans[tag]=optimize_pyramids<0.0 ? group.normal_plans[tag] : group.plan
        end
        node_clock=group.node_high_water;element_clock=group.element_high_water
    end
    tags=sort!(Int.(collect(regions)))
    # Native default optimization follows construction of every selected
    # component, including all Pyramid births. Each Region then consumes the
    # same continuing global element factory in declared Region order.
    for tag in tags
        product=products[tag]
        tet=only(filter(block->block.msh==4,product.mesh.blocks))
        pyramid=only(filter(block->block.msh==7,product.mesh.blocks))
        mesh=Mesh(product.mesh.coords;tets=tet.nodes)
        work=_ModelNativeVolumeBirths(node_clock,element_clock)
        work.source=Mesh(product.constructor.coordinates)
        work.node_tags=copy(product.boundary_node_tags)
        work.cell_source=mesh;work.element_tags=copy(product.element_tags[1:size(tet.nodes,2)])
        optimized=_model_native_optimize_volume(mesh,options,work,input,
            NTuple{3,Int32}[],NTuple{2,Int32}[],caller)
        element_clock=work.element_high_water
        blocks=ElementBlock[ElementBlock(4,copy(optimized.tets)),ElementBlock(7,copy(pyramid.nodes))]
        mixed=MixedMesh(optimized.coords,blocks)
        optimize_pyramids<0.0 && _model_native_hybrid_normal_apex_fixed_certificate(
            plans[tag],mixed,product.boundary_node_tags,caller)
        validate(mixed).ok || throw(ArgumentError("$caller: hybrid batch produces invalid typed cells"))
        _model_native_hybrid_volume_certificate(mixed,product.boundary_node_tags,product.constructor,caller)
        _model_native_hybrid_periodic_boundary_certificate(m,tag,product.constructor,periodic_surface_pairs,caller)
        geometry=_model_explicit_volume_geometry(m,tag,caller;_mesh_input=input,_mesh1d_options=options)
        _model_certify_explicit_volume_semantics(mixed,tag,geometry.expected_volume,
            geometry.comparison_scale,caller)
        products[tag]=_ModelNativeHybridVolumeProduction(mixed,copy(product.boundary_node_tags),node_clock,
            product.temporary_elements,vcat(work.element_tags,product.element_tags[size(tet.nodes,2)+1:end]),
            element_clock,input,product.constructor,copy(product.registration_rows))
    end
    # The operation publishes its final aggregate clocks on every product;
    # each retained UID still belongs to its explicit Region product once.
    for tag in tags
        product=products[tag]
        products[tag]=_ModelNativeHybridVolumeProduction(product.mesh,product.boundary_node_tags,
            node_clock,product.temporary_elements,product.element_tags,element_clock,input,product.constructor,
            copy(product.registration_rows))
    end
    return _ModelNativeHybridBatchProduction(input,tags,groups,products,node_clock,element_clock,options.random_state[])
end
