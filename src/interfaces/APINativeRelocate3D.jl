# meshRelocateVertex.cpp's tetrahedron RelocateVertices path. Public ownership,
# rather than a topological boundary mask, admits moving vertices.
struct _NativeRelocate3DQuality
    rule::_NativeCellSampleRule
    coordinates::Matrix{Float64}
    derivatives::NTuple{3,Matrix{Float64}}
end

struct _NativeRelocate3DCell
    cell::_Dim01Cell
    quality::Union{Nothing,_NativeRelocate3DQuality}
end

function _native_relocate3d_cell(cell,caller)
    # Upstream RelocateVertices(GRegion) builds vertex adjacency from the
    # region's tetrahedra, pyramids, prisms and hexahedra only; every other
    # stored cell (embedded faces, polyhedra, trihedra, uncatalogued records)
    # is excluded from relocation instead of rejecting the region.
    spec=get(Elements._MSH_CATALOG_LOOKUP,Int(cell.msh),nothing)
    (spec===nothing || spec.dim!=3 || !(spec.family in (:tet,:hex,:pri,:pyr))) &&
        return nothing
    cell.msh==4 && return _NativeRelocate3DCell(cell,nothing)
    rule=_native_cell_sample_rule(cell.msh,caller)
    coordinates=Matrix{Float64}(undef,length(cell.nodes),3)
    derivatives=ntuple(_->zeros(Float64,size(rule.ideal_gradients[1],1),3),3)
    return _NativeRelocate3DCell(cell,_NativeRelocate3DQuality(rule,coordinates,derivatives))
end

function _native_relocate3d_sicn(assembly,entry,node,trial)
    cell=entry.cell;work=entry.quality::_NativeRelocate3DQuality
    coordinates=work.coordinates;derivatives=work.derivatives
    for slot in eachindex(cell.nodes)
        position=cell.nodes[slot]
        point=position==node ? trial : assembly.nodes[position].point
        for axis in 1:3
            coordinates[slot,axis]=point[axis]
        end
    end
    for axis in 1:3
        fill!(derivatives[axis],0.0)
        _native_bezier_product!(derivatives[axis],work.rule.ideal_gradients[axis],coordinates)
    end
    minimum_value=0.0
    for sample in axes(derivatives[1],1)
        a=(derivatives[1][sample,1],derivatives[1][sample,2],derivatives[1][sample,3])
        b=(derivatives[2][sample,1],derivatives[2][sample,2],derivatives[2][sample,3])
        c=(derivatives[3][sample,1],derivatives[3][sample,2],derivatives[3][sample,3])
        value=_native_cell_sample_condition3(a,b,c)
        sample==1 && (minimum_value=value)
        value<minimum_value && (minimum_value=value)
    end
    return minimum_value
end

function _native_relocate3d_objective(assembly,node,incident,target,xi)
    original=assembly.nodes[node].point
    trial=ntuple(axis->(1.0-xi)*original[axis]+xi*target[axis],3)
    minimum_quality=1.0
    for entry in incident
        cell=entry.cell
        quality=if cell.msh==4
            points=ntuple(slot->cell.nodes[slot]==node ? trial :
                assembly.nodes[cell.nodes[slot]].point,4)
            value,_=_gopt_gamma(points...)
            Model.orient3(points...)>0 ? -value : value
        else
            abs(_native_relocate3d_sicn(assembly,entry,node,trial))*0.2
        end
        # std::min(Q,minQual) preserves Q when either comparison is unordered.
        minimum_quality=minimum_quality<quality ? minimum_quality : quality
    end
    return minimum_quality
end

function _native_relocate3d_golden(assembly,node,incident,target)
    lambda=0.5*(sqrt(5.0)-1.0)
    mu=0.5*(3.0-sqrt(5.0))
    a=0.0;b=2.0;tolerance=1e-2
    x1=b-lambda*(b-a);x2=a+lambda*(b-a)
    f1=_native_relocate3d_objective(assembly,node,incident,target,x1)
    f2=_native_relocate3d_objective(assembly,node,incident,target,x2)
    while !(abs(b-a)<tolerance)
        if f1<f2
            a=x1
            abs(b-a)<tolerance && break
            x1=x2;f1=f2;x2=b-mu*(b-a)
            f2=_native_relocate3d_objective(assembly,node,incident,target,x2)
        else
            b=x2
            abs(b-a)<tolerance && break
            x2=x1;f2=f1;x1=a+mu*(b-a)
            f1=_native_relocate3d_objective(assembly,node,incident,target,x1)
        end
    end
    return a
end

function _native_relocate3d_vertex!(assembly,node,incident)
    previous=assembly.nodes[node]
    previous.owner[1]==3 || return nothing
    x=0.0;y=0.0;z=0.0;count=0
    for entry in incident
        cell=entry.cell
        sx=0.0;sy=0.0;sz=0.0
        for position in cell.nodes
            point=assembly.nodes[position].point
            sx+=point[1];sy+=point[2];sz+=point[3]
        end
        x+=sx;y+=sy;z+=sz;count+=length(cell.nodes)
    end
    target=(x/count,y/count,z/count)
    before=_native_relocate3d_objective(assembly,node,incident,target,0.0)
    after=_native_relocate3d_objective(assembly,node,incident,target,1.0)
    after>before && (assembly.nodes[node]=_Dim01Node(previous.tag,target,
        previous.owner,previous.parameter))
    xi=_native_relocate3d_golden(assembly,node,incident,target)
    original=assembly.nodes[node].point
    # The CPP final blend multiplies the accumulated sum before dividing;
    # predividing the target changes partial relocation by a few ULPs.
    point=((1.0-xi)*original[1]+xi*x/count,
        (1.0-xi)*original[2]+xi*y/count,(1.0-xi)*original[3]+xi*z/count)
    assembly.nodes[node]=_Dim01Node(previous.tag,point,previous.owner,previous.parameter)
    return nothing
end

function _native_relocate3d_locked!(iterations,pairs,caller)
    model=_model_locked();cached=LAST_MESH[]
    class=cached===nothing ? nothing : _cached_classification_locked(cached)
    actual,actual_class=_dim01_actual_cache(cached,class;physical_names=model.physical_names)
    isempty(pairs) || @warn "Entity selection in the mesh optimizer is ignored"
    source=Model._model_mesh_transaction_copy(model)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    mapping=_dim01_ingest!(assembly,source,actual,actual_class;preserve_source_order=true)
    admitted=iszero(OPTIONS["General.AbortOnError"]) || NATIVE_LOOKUP_QUERY_ERRORS[]==0
    if admitted && iterations!=0 && iterations>=-1
        regions=Dict{Int,Vector{_NativeRelocate3DCell}}()
        for cell in assembly.cells
            cell.owner[1]==3 || continue
            entry=_native_relocate3d_cell(cell,caller)
            entry===nothing && continue
            push!(get!(()->_NativeRelocate3DCell[],regions,Int(cell.owner[2])),entry)
        end
        for (_,region) in model_entities(source,3)
            cells=get(regions,region,nothing)
            cells===nothing && continue
            adjacency=Dict{Int32,Vector{_NativeRelocate3DCell}}()
            for entry in cells,position in entry.cell.nodes
                push!(get!(()->_NativeRelocate3DCell[],adjacency,position),entry)
            end
            order=sort!(collect(keys(adjacency));by=node->assembly.nodes[node].tag)
            # The CPP's relax argument is unused by relocateVertexGolden.
            # Its loop intentionally performs niter+2 passes, including one
            # pass for -1; zero is an explicit early return before that loop.
            for _ in 1:(iterations+2),node in order
                _native_relocate3d_vertex!(assembly,node,adjacency[node])
            end
        end
    end
    isempty(assembly.nodes) && isempty(assembly.cells) && return nothing
    renumber=!iszero(OPTIONS["Mesh.Renumber"])
    plan=_dim01_finish!(assembly,source,model_entities(source),Set{Tuple{Int,Int32}}(),
        falses(length(assembly.nodes)),renumber,false;
        save_all=!iszero(OPTIONS["Mesh.SaveAll"]),source_mesh=actual,
        source_node_positions=mapping)
    if cached isa MixedMesh && (!isempty(cached.ancillary_sections) ||
        !isempty(cached.data_sections) || cached.partition_data!==nothing)
        cached.entity_data!==nothing ||
            (isempty(cached.data_sections) && cached.partition_data===nothing) ||
            throw(ArgumentError("$caller: retained file rows require explicit public-ID metadata"))
        sections=_reclassify_data_sections(cached,plan.mesh,plan,caller)
        partition=_reclassify_partition_data(cached,plan.mesh,plan,caller)
        preserved=MixedMesh(plan.mesh.coords,plan.mesh.blocks;
            physical_names=plan.mesh.physical_names,entity_data=plan.mesh.entity_data,
            elementary_entities=plan.mesh.elementary_entities,periodic_links=plan.mesh.periodic_links,
            ancillary_sections=cached.ancillary_sections,data_sections=sections,partition_data=partition)
        plan=merge(plan,(mesh=preserved,class=_mixed_rebind_class(plan.class,preserved),))
    end
    _commit_tagged_plan_locked!(plan;preserve_visibility=true,
        element_tag_map=plan.visibility_tag_map)
    renumber && _native_lookup_invalidate_locked!()
    return nothing
end
