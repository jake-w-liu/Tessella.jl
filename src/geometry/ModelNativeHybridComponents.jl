# Connected Region recovery binds actual lower Face closures and keeps their
# declared Region ownership through shared-panel partitioning and publication.
# `embedded` carries the deduplicated joint embedded pool whenever any member
# Region declares `Embedded{}` entities (see the pool section of
# ModelNativeHybridEmbeddedConstraints.jl); `nothing` when no member embeds.
struct _ModelNativeHybridComponentReceipt
    input::_ModelMeshAuthority
    regions::Vector{Int}
    boundaries::Dict{Int,Vector{Int}}
    face_regions::Dict{Int,Vector{Int}}
    constructor::NamedTuple
    embedded::Union{Nothing,_ModelNativeHybridEmbeddedComponentInput}
end

struct _ModelNativeHybridRecoveryReceipt
    component::_ModelNativeHybridComponentReceipt
    plan::_ModelNativeHybridQuadPlan
    order::Vector{Int}
    kernel::ModelNativeInitialTet.Kernel
    points::Vector{NTuple{3,Float64}}
    facets::Vector{NTuple{3,Int32}}
    multiplicities::Dict{NTuple{3,Int32},Int}
end

# Project already allocated global births onto a Region's actual lower
# constructor. No vertex, owner, normal or factory UID is inferred from XYZ.
function _model_native_hybrid_region_normal_plan(m,component,plan,tag,caller)
    tag in component.regions && plan.source_uids==component.constructor.node_uids &&
        isequal(plan.source_coordinates,component.constructor.coordinates) || throw(ArgumentError(
        "$caller: hybrid normal batch loses its actual connected Source"))
    constructor=_model_native_hybrid_volume_input(m,tag,component.input,caller)
    constructor!==nothing || throw(ArgumentError("$caller: normal batch Region lacks its actual constructor"))
    global_rows=Dict(uid=>row for (row,uid) in enumerate(plan.source_uids))
    all(uid->haskey(global_rows,uid),constructor.node_uids) &&
        isequal(plan.source_coordinates[:,Int[global_rows[uid] for uid in constructor.node_uids]],
                constructor.coordinates) || throw(ArgumentError(
        "$caller: normal batch Region changes a witnessed Source UID/XYZ row"))
    source=copy(constructor.node_uids)
    local_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(source))
    originals=Dict((birth.owner,birth.element_uid)=>birth for birth in plan.births)
    length(originals)==length(plan.births) || throw(ArgumentError("$caller: normal batch repeats a global Quad birth"))
    coordinates=Matrix{Float64}(undef,3,length(source)+length(constructor.quads))
    coordinates[:,1:length(source)].=constructor.coordinates
    node_uids=copy(source);births=_ModelNativeHybridQuadBirth[]
    for (ordinal,quad) in enumerate(constructor.quads)
        original=get(originals,(quad.owner,quad.uid),nothing)
        original!==nothing && original.base_uids==quad.nodes &&
            component.face_regions[Int(quad.owner[2])]==[tag] &&
            plan.node_uids[original.apex_row]==original.apex_uid &&
            isequal(Tuple(plan.coordinates[:,original.apex_row]),original.point) || throw(ArgumentError(
            "$caller: normal batch Region has a missing or foreign actual Quad birth"))
        base=ntuple(slot->local_rows[quad.nodes[slot]],4)
        apex=Int32(length(source)+ordinal)
        coordinates[:,apex].=original.point;push!(node_uids,original.apex_uid)
        push!(births,_ModelNativeHybridQuadBirth(original.owner,original.element_uid,
            original.base_uids,base,original.apex_uid,apex,original.point))
    end
    regional_plan=_ModelNativeHybridQuadPlan(source,copy(constructor.coordinates),node_uids,
        coordinates,zeros(Int32,3,0),births,plan.initial_clock,plan.node_high_water)
    _model_native_hybrid_normal_apex_certificate(regional_plan,constructor,caller)
    return regional_plan
end

function _model_native_hybrid_component_receipt(m,region_tags,input,caller)
    tags=Int.(collect(region_tags))
    !isempty(tags) && allunique(tags) && all(tag->haskey(m.volumes,tag),tags) ||
        throw(ArgumentError("$caller: connected hybrid receipt lacks unique declared Region inputs"))
    sort!(tags)
    input isa _ModelMeshAuthority && input.mesh!==nothing && input.mesh.entity_data!==nothing ||
        throw(ArgumentError("$caller: connected hybrid receipt lacks current actual lower authority"))
    boundaries=Dict{Int,Vector{Int}}();face_regions=Dict{Int,Vector{Int}}()
    signatures=Set{Tuple}()
    for tag in tags
        faces=Int.(_model_volume_boundary_surfaces(m,tag,caller))
        !isempty(faces) && allunique(abs.(faces)) || throw(ArgumentError(
            "$caller: connected hybrid Region repeats or omits a declared Face"))
        signature=Tuple(sort(abs.(faces)))
        signature in signatures && throw(ArgumentError(
            "$caller: connected hybrid Regions have ambiguous identical Face closures"))
        push!(signatures,signature);boundaries[tag]=faces
        for face in abs.(faces)
            owners=get!(()->Int[],face_regions,face);push!(owners,tag)
            length(owners)<=2 || throw(ArgumentError(
                "$caller: connected hybrid Face has more than two Region owners"))
        end
    end
    # This receipt is one connected native meshGRegion group. Disjoint groups
    # retain their own sequential RNG/factory continuation.
    reached=Set([first(tags)]);queue=copy(tags[1:1]);cursor=1
    while cursor<=length(queue)
        tag=queue[cursor];cursor+=1
        for face in abs.(boundaries[tag]),neighbor in face_regions[face]
            if !(neighbor in reached);push!(reached,neighbor);push!(queue,neighbor);end
        end
    end
    reached==Set(tags) || throw(ArgumentError(
        "$caller: disconnected hybrid Regions require separate native constructor groups"))
    faces=sort!(collect(keys(face_regions)));actual=input.mesh;data=actual.entity_data
    used=Set{UInt64}()
    for face in faces
        product=get(input.products,(2,Int32(face)),nothing)
        product!==nothing || throw(ArgumentError("$caller: connected hybrid Face lacks actual lower cells"))
        for block in product.blocks
            block.msh in (2,3) || throw(ArgumentError(
                "$caller: connected hybrid Face requires actual linear Tri/Quad cells"))
            block.msh==3 && !isempty(block.nodes) && length(face_regions[face])==2 &&
                throw(ArgumentError("$caller: shared Quad recovery lacks a qualified native two-Region Pyramid contract"))
            for row in block.nodes;push!(used,data.external_node_tags[product.rows[row]]);end
        end
    end
    # The joint embedded pool is built in two phases: witnesses supply their
    # actual node UIDs to the shared Source before the constructor exists;
    # record classification then runs against the finished constructor's
    # literal panels (see ModelNativeHybridEmbeddedConstraints.jl).
    pool_witnesses=_model_native_hybrid_embedded_pool_witnesses(m,tags,input,caller)
    pool_witnesses===nothing || union!(used,pool_witnesses.uids)
    ids=sort!(collect(used));positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    all(uid->haskey(positions,uid),ids) || throw(ArgumentError(
        "$caller: connected hybrid Face incidence omits an actual current node UID"))
    coordinates=actual.coords[:,Int32[positions[uid] for uid in ids]]
    constructor=_model_native_hybrid_constructor_input(input,faces,ids,coordinates,caller)
    embedded=_model_native_hybrid_embedded_component_input(m,tags,input,boundaries,
        face_regions,constructor,caller;witnesses=pool_witnesses)
    return _ModelNativeHybridComponentReceipt(input,tags,boundaries,face_regions,constructor,embedded)
end

function _model_native_hybrid_component_facets(receipt,plan,caller)
    receipt.input===receipt.constructor.input && plan.source_uids==receipt.constructor.node_uids &&
        isequal(plan.source_coordinates,receipt.constructor.coordinates) || throw(ArgumentError(
        "$caller: connected hybrid facet receipt changes its owned constructor"))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(plan.node_uids))
    facets=Dict{NTuple{3,Int32},Tuple{Int,UInt64}}()
    pool=receipt.embedded
    function add(rows,owner,uid)
        key=Tuple(sort!(Int32[rows...]))
        haskey(facets,key) && throw(ArgumentError(
            "$caller: connected hybrid constructor repeats a typed Face panel"))
        haskey(receipt.face_regions,owner) || throw(ArgumentError(
            "$caller: connected hybrid panel has a foreign Face owner"))
        facets[key]=(owner,uid)
    end
    emitted=Set{Int}()
    for triangle in receipt.constructor.triangles
        face=Int(triangle.owner[2])
        subs=pool===nothing ? nothing : get(pool.splits,face,nothing)
        if subs===nothing
            add(ntuple(slot->positions[triangle.nodes[slot]],3),face,triangle.uid)
        elseif !(face in emitted)
            # A certified whole-Face split keeps the shared Face owner on every
            # published sub-facet; provenance records the Face's first literal
            # panel UID since sub-facets no longer descend from one cell.
            push!(emitted,face)
            for sub in subs
                add(ntuple(slot->positions[sub[slot]],3),face,triangle.uid)
            end
        end
    end
    for birth in plan.births,side in 1:4
        add((birth.base_rows[side],birth.base_rows[mod1(side+1,4)],birth.apex_row),
            Int(birth.owner[2]),birth.element_uid)
    end
    for tag in receipt.regions
        edges=Dict{Tuple{Int32,Int32},Int}()
        for (face,(owner,uid)) in facets
            tag in receipt.face_regions[owner] || continue
            for side in 1:3
                key=minmax(face[side],face[mod1(side+1,3)])
                edges[key]=get(edges,key,0)+1
            end
        end
        !isempty(edges) && all(==(2),values(edges)) || throw(ArgumentError(
            "$caller: connected hybrid Region[$tag] does not have a closed actual constructor shell"))
    end
    return facets
end

function _model_native_hybrid_component_partition(mesh,node_uids,cell_uids,receipt,plan,caller)
    length(node_uids)==nnodes(mesh) && length(node_uids)>=length(plan.node_uids) && allunique(node_uids) && all(!iszero,node_uids) &&
        node_uids[1:length(plan.node_uids)]==plan.node_uids &&
        isequal(mesh.coords[:,1:length(plan.source_uids)],plan.source_coordinates) &&
        length(cell_uids)==ntets(mesh) && allunique(cell_uids) && all(!iszero,cell_uids) &&
        validate(mesh).ok || throw(ArgumentError(
        "$caller: connected hybrid partition lacks complete native UID/cell constructor incidence"))
    panels=_model_native_hybrid_component_facets(receipt,plan,caller)
    incidence=Dict{NTuple{3,Int32},Vector{Int}}()
    for column in axes(mesh.tets,2),omit in 1:4
        key=Tuple(sort!(Int32[mesh.tets[slot,column] for slot in 1:4 if slot!=omit]))
        push!(get!(()->Int[],incidence,key),column)
    end
    all(rows->length(rows)<=2,values(incidence)) || throw(ArgumentError(
        "$caller: connected hybrid native cells have non-manifold face incidence"))
    all(face->haskey(incidence,face),keys(panels)) || throw(ArgumentError(
        "$caller: connected hybrid native cells omit a current actual Face panel"))
    neighbors=[Int[] for _ in axes(mesh.tets,2)];touching=[Set{Int}() for _ in axes(mesh.tets,2)]
    for (face,rows) in incidence
        if haskey(panels,face)
            owner=panels[face][1]
            length(rows)==length(receipt.face_regions[owner]) || throw(ArgumentError(
                "$caller: connected hybrid panel multiplicity changes its declared Region ownership"))
            for row in rows;push!(touching[row],owner);end
        elseif length(rows)==2
            first,last=rows;push!(neighbors[first],last);push!(neighbors[last],first)
        else
            throw(ArgumentError("$caller: connected hybrid native exterior has an unowned panel"))
        end
    end
    owners=zeros(Int32,ntets(mesh));seen=falses(ntets(mesh));used_regions=Set{Int}()
    for seed in axes(mesh.tets,2)
        seen[seed] && continue
        seen[seed]=true;queue=Int[seed];cursor=1;members=Int[];faces=Set{Int}()
        while cursor<=length(queue)
            row=queue[cursor];cursor+=1;push!(members,row);union!(faces,touching[row])
            for neighbor in neighbors[row]
                if !seen[neighbor];seen[neighbor]=true;push!(queue,neighbor);end
            end
        end
        candidates=[tag for tag in receipt.regions if Set(abs.(receipt.boundaries[tag]))==faces]
        length(candidates)==1 && !(only(candidates) in used_regions) || throw(ArgumentError(
            "$caller: connected hybrid native cell component lacks a unique declared Face closure"))
        tag=only(candidates);push!(used_regions,tag);owners[members].=Int32(tag)
    end
    used_regions==Set(receipt.regions) || throw(ArgumentError(
        "$caller: connected hybrid partition omits a declared Region"))
    return owners
end

function _model_native_hybrid_component_recovery_receipt(component,plan,order,kernel,caller)
    component isa _ModelNativeHybridComponentReceipt && kernel isa ModelNativeInitialTet.Kernel ||
        throw(ArgumentError("$caller: hybrid recovery lacks an operation-owned connected Region receipt"))
    # Rebind the selected constructor against the current lower products. A
    # stored Face list or coordinate match alone cannot supply this provenance.
    current=_model_native_hybrid_constructor_input(component.input,
        sort!(collect(keys(component.face_regions))),component.constructor.node_uids,
        component.constructor.coordinates,caller)
    isequal(current,component.constructor) || throw(ArgumentError(
        "$caller: hybrid recovery changes its actual current lower constructor"))
    length(order)==length(plan.node_uids)==kernel.count && allunique(order) &&
        all(row->0<row<=kernel.count,order) || throw(ArgumentError(
        "$caller: hybrid recovery changes its initial constructor row permutation"))
    panels=_model_native_hybrid_component_facets(component,plan,caller)
    facets=_model_native_hybrid_recovery_facets(plan,component.constructor,order,caller)
    component.embedded===nothing ||
        (facets=_model_native_hybrid_pool_split_facets(facets,component.embedded,plan,order,caller))
    inverse=invperm(order)
    multiplicities=Dict(Tuple(sort!(Int32[inverse[row] for row in panel]))=>
        length(component.face_regions[owner]) for (panel,(owner,uid)) in panels)
    Set(Tuple(sort!(collect(face))) for face in facets)==Set(keys(multiplicities)) ||
        throw(ArgumentError("$caller: hybrid recovery omits a certified current Face panel"))
    return _ModelNativeHybridRecoveryReceipt(component,plan,Int.(order),kernel,
        copy(kernel.points),copy(facets),multiplicities)
end

function _model_native_hybrid_recovery_multiplicities(receipt,kernel,faces)
    receipt isa _ModelNativeHybridRecoveryReceipt && receipt.kernel===kernel &&
        isequal(receipt.points,kernel.points) && receipt.facets==faces || throw(ArgumentError(
        "hybrid connected recovery: missing, foreign, or changed initial constructor receipt"))
    fresh=_model_native_hybrid_component_recovery_receipt(receipt.component,receipt.plan,
        receipt.order,kernel,"hybrid connected recovery")
    fresh.multiplicities==receipt.multiplicities || throw(ArgumentError(
        "hybrid connected recovery: current typed Face ownership changed"))
    return copy(receipt.multiplicities)
end

struct _ModelNativeHybridRegionBackground{T}
    fields::Dict{Int32,T}
end
function (background::_ModelNativeHybridRegionBackground)(x,y,z,tag)
    field=get(background.fields,Int32(tag),nothing)
    field===nothing && throw(ArgumentError("hybrid connected sizing: foreign current Region owner"))
    return field(x,y,z)
end

function _model_native_hybrid_connected_production(m,component,options,births,caller;
                                                   size_field=nothing,optimize_pyramids::Float64=0.0,
                                                   periodic_surface_pairs=nothing)
    # Tet optimization is deferred until every selected component has
    # completed construction with its continuing global factory clocks.
    # Mesh.OptimizeNetgen is an upstream no-op on hybrid regions, so only the
    # deferred default optimization counter remains gated.
    component isa _ModelNativeHybridComponentReceipt && !isempty(component.regions) &&
        m.meshing.algorithm_3d==1 && options.optimize==0 ||
        throw(ArgumentError("$caller: hybrid connected production requires its recorded batch and deferred default optimization"))
    isfinite(optimize_pyramids) || throw(ArgumentError("$caller: hybrid Pyramid setting must be finite"))
    fresh=_model_native_hybrid_component_receipt(m,component.regions,component.input,caller)
    fresh.boundaries==component.boundaries && fresh.face_regions==component.face_regions &&
        isequal(fresh.constructor,component.constructor) &&
        _model_native_hybrid_pool_equal(fresh.embedded,component.embedded) || throw(ArgumentError(
        "$caller: hybrid connected production changes its current declared Region or lower constructor receipt"))
    for tag in component.regions
        !get(m.meshing.reverse,(3,tag),false) && get(m.meshing.smoothing,(3,tag),0)==0 ||
            throw(ArgumentError(
                "$caller: hybrid connected production has unsupported current Region planner attributes"))
    end
    pool=component.embedded
    # The joint embedded pool supplies per-Region constructor fields; the
    # normal-apex Pyramid relocation tier keeps its own Region-local Source
    # projection and stays a loud specific reject for embedded batches.
    (pool===nothing || optimize_pyramids>=0.0) || throw(ArgumentError(
        "$caller: hybrid embedded batch requires nonnegative Pyramid optimization"))
    _model_native_hybrid_batch_periodic_certificate(m,component,periodic_surface_pairs,optimize_pyramids,caller)
    surface_witnesses=pool===nothing ? () : Tuple(pool.surfaces)
    triangles=_model_native_hybrid_global_triangles(m,component.input,caller;
        surface_witnesses)
    seeds=_model_native_hybrid_global_size_seeds(m,component.input,caller;
        surface_witnesses)
    constructor=component.constructor
    work=_ModelNativeVolumeBirths(births.node_high_water,births.element_high_water)
    plan=_model_native_hybrid_quad_plan(constructor.node_uids,constructor.coordinates,
        constructor.quads,work.node_high_water,caller;
        normal_factor=optimize_pyramids<0.0 ? -optimize_pyramids : 0.0)
    normal_plans=Dict{Int,_ModelNativeHybridQuadPlan}()
    if optimize_pyramids<0.0
        for tag in component.regions
            normal_plans[tag]=_model_native_hybrid_region_normal_plan(m,component,plan,tag,caller)
        end
    end
    work.source=Mesh(plan.source_coordinates);work.node_tags=copy(plan.node_uids)
    work.node_high_water=plan.node_high_water
    order=sortperm(plan.node_uids);kernel=Ref{ModelNativeInitialTet.Kernel}()
    _model_native_initial_delaunay_fill(Mesh(plan.coordinates;tris=plan.fan_facets),order,
        options,caller;births=work,kernel_output=kernel)
    facets=_model_native_hybrid_recovery_facets(plan,constructor,order,caller)
    embedded_segments=NTuple{2,Int32}[];embedded_facets=NTuple{3,Int32}[]
    if pool!==nothing
        # Shared interface panels carrying on-panel records are replaced by
        # their certified CDT sub-facets; strictly-interior segments and sheet
        # facets join the recovery as protected interior constraints.
        facets=_model_native_hybrid_pool_split_facets(facets,pool,plan,order,caller)
        embedded_segments,embedded_facets=
            _model_native_hybrid_pool_kernel_constraints(pool,plan.node_uids,order,caller)
    end
    recovery=_model_native_hybrid_component_recovery_receipt(component,plan,order,kernel[],caller)
    initial=ModelNativeHybridInterior.recover(kernel[],facets;component=recovery,
        embedded_segments=embedded_segments,embedded_facets=embedded_facets)
    mesh=Mesh(plan.coordinates;tets=reduce(hcat,Int32[order[row] for row in cell] for cell in initial.cells))
    _model_native_volume_cell_births!(work,mesh)
    owners=_model_native_hybrid_component_partition(mesh,plan.node_uids,work.element_tags,component,plan,caller)
    mesh=Mesh(mesh.coords;tets=mesh.tets,tet_tag=owners)
    fields=Dict(Int32(tag)=>_ModelAlgorithmBackgroundSize(m,3,tag,options,size_field,caller)
        for tag in component.regions)
    minimize=unique(_resolved_extend_from_boundary(m,3,tag)==2 for tag in component.regions)
    extend=unique(_resolved_extend_from_boundary(m,3,tag)>0 for tag in component.regions)
    length(minimize)==length(extend)==1 || throw(ArgumentError(
        "$caller: hybrid connected production has inconsistent native global sizing policy"))
    insertion=ModelNativeHybridInsertion.refine(mesh,work.node_tags,work.element_tags,
        triangles,_ModelNativeHybridRegionBackground(fields),caller;
        minimize=only(minimize),extend=only(extend),region_background=true,
        node_high_water=work.node_high_water,element_high_water=work.element_high_water,seeds,
        metric=_volume_metric_insertion_field(m,size_field,options,caller))
    smoothed=ModelNativeHybridInsertion.smooth(insertion,length(plan.node_uids),m.meshing.native_smoothing,caller)
    # Native relocation is invoked separately for each restored GRegion.
    # Its two-hop Tet neighborhood must not acquire the other Region's star
    # merely because both Regions share a fixed lower vertex.
    coordinates=copy(smoothed.mesh.coords)
    for tag in (optimize_pyramids<0.0 ? Int[] : component.regions)
        birth_indices=findall(birth->only(component.face_regions[Int(birth.owner[2])])==tag,plan.births)
        local_positions=Dict(index=>ordinal for (ordinal,index) in enumerate(birth_indices))
        local_order=Int[local_positions[index] for index in constructor.pyramid_order if haskey(local_positions,index)]
        local_plan=_ModelNativeHybridQuadPlan(plan.source_uids,plan.source_coordinates,plan.node_uids,
            coordinates,zeros(Int32,3,0),plan.births[birth_indices],plan.initial_clock,plan.node_high_water)
        columns=findall(==(tag),smoothed.mesh.tet_tag)
        local_mesh=Mesh(coordinates;tets=smoothed.mesh.tets[:,columns])
        relocated=(optimize_pyramids==1.0 ? ModelNativeHybridWinslow : ModelNativeHybridRelocate).relocate(
            local_plan,local_mesh,caller;
            birth_order=local_order,node_uids=insertion.node_numbers)
        coordinates=relocated.coords
    end
    relocated=Mesh(coordinates;tets=smoothed.mesh.tets,tet_tag=smoothed.mesh.tet_tag)
    pyramids=_model_native_hybrid_pyramid_production(plan,relocated.coords,
        insertion.element_high_water,caller;birth_order=constructor.pyramid_order)
    # The region of an apex is selected by the recorded Quad owner. Shared
    # Quads were rejected before kernel/RNG mutation by the component receipt.
    pyramid_owners=Int32[only(component.face_regions[Int(plan.births[index].owner[2])])
        for index in constructor.pyramid_order]
    products=Dict{Int,_ModelNativeHybridVolumeProduction}()
    for tag in component.regions
        tet_columns=findall(==(tag),relocated.tet_tag);pyramid_columns=findall(==(tag),pyramid_owners)
        selected=sort!(unique(vcat(vec(relocated.tets[:,tet_columns]),vec(pyramids.nodes[:,pyramid_columns]))))
        inverse=zeros(Int32,nnodes(relocated));inverse[selected].=Int32.(eachindex(selected))
        blocks=ElementBlock[ElementBlock(4,inverse[relocated.tets[:,tet_columns]]),
            ElementBlock(7,inverse[pyramids.nodes[:,pyramid_columns]])]
        mixed=MixedMesh(relocated.coords[:,selected],blocks)
        node_uids=insertion.node_numbers[selected]
        region_constructor=_model_native_hybrid_volume_input(m,tag,component.input,caller;
            embedded_pool=pool)
        region_constructor!==nothing || throw(ArgumentError("$caller: hybrid connected Region loses its current lower constructor"))
        _model_native_hybrid_volume_certificate(mixed,node_uids,region_constructor,caller;
            embedded_pool=pool)
        _model_native_hybrid_embedded_point_certificate(mixed,node_uids,region_constructor,caller)
        _model_native_hybrid_embedded_curve_certificate(mixed,node_uids,region_constructor,caller)
        _model_native_hybrid_embedded_surface_certificate(mixed,node_uids,region_constructor,caller)
        _model_native_hybrid_embedded_constraints_certificate(mixed,node_uids,region_constructor,caller)
        _model_native_hybrid_periodic_boundary_certificate(m,tag,region_constructor,periodic_surface_pairs,caller)
        products[tag]=_ModelNativeHybridVolumeProduction(mixed,node_uids,insertion.node_high_water,
            work.temporary_elements,vcat(insertion.element_numbers[tet_columns],pyramids.element_tags[pyramid_columns]),
            pyramids.element_high_water,component.input,region_constructor,
            _model_native_hybrid_registration_rows(region_constructor,plan,insertion,node_uids,caller),
            pool)
    end
    pool===nothing || _model_native_hybrid_embedded_pool_certificate(products,pool,caller)
    return (input=component.input,component=component,plan=plan,normal_plans=normal_plans,products=products,
        node_high_water=insertion.node_high_water,element_high_water=pyramids.element_high_water,
        random_state=options.random_state[])
end
