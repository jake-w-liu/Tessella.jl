# An order-two parent can carry native periodic interpolation correspondences
# into its actual linear children. This proof never obtains identity from XYZ.
struct _SurfacePeriodicRefineReceipt
    source::MixedMesh
    source_class::Any
    source_snapshot::NamedTuple
    children::Vector{NamedTuple}
end

_surface_periodic_float_bits(values)=Tuple(reinterpret(UInt64,Float64[values...]))

function _surface_periodic_source_snapshot(source,class,caller)
    class!==nothing && hasproperty(class,:mesh) && class.mesh===source ||
        throw(ArgumentError("$caller: periodic parent classification is not bound to its actual mesh"))
    public=class===nothing ? nothing : class.public_tags
    public!==nothing || throw(ArgumentError("$caller: periodic parent has no actual UID table"))
    tags=public.node_tags
    length(tags)==nnodes(source) && length(unique(tags))==length(tags) || throw(ArgumentError(
        "$caller: periodic parent has an invalid node UID table"))
    length(class.node_entities)==length(class.node_parametric)==length(tags) || throw(ArgumentError(
        "$caller: periodic parent classification loses vertex bindings"))
    cells=NamedTuple[];used=Set{Int32}()
    for (msh,dim,offset,nodes,owners) in _cache_catalog(source,class)
        dim==2 || continue
        spec=msh_spec(msh)
        spec.family in (:tri,:qua) && spec.order==2 && !spec.serendipity || throw(ArgumentError(
            "$caller: periodic refinement requires an actual complete quadratic Surface parent"))
        for column in axes(nodes,2)
            rows=Int32[nodes[:,column]...]
            all(row->1<=row<=length(tags),rows) || throw(ArgumentError(
                "$caller: periodic parent has invalid connectivity"))
            union!(used,rows)
            push!(cells,(msh=Int(msh),owner=owners[column],
                uid=public.element_tags[offset+column],nodes=Tuple(tags[row] for row in rows)))
        end
    end
    length(unique(cell.uid for cell in cells))==length(cells) || throw(ArgumentError(
        "$caller: periodic parent repeats a Surface cell UID"))
    vertices=NamedTuple[]
    for row in sort!(collect(used))
        raw=class.node_parametric[row]
        push!(vertices,(uid=tags[row],owner=class.node_entities[row],
            xyz=_surface_periodic_float_bits(@view(source.coords[:,row])),
            parameters=raw===nothing ? nothing : _surface_periodic_float_bits(raw)))
    end
    return (;vertices,cells)
end

function _surface_periodic_expected_children(snapshot,caller)
    result=NamedTuple[]
    for parent in snapshot.cells
        family=msh_spec(parent.msh).family
        linear=family===:tri ? 2 : 3
        for (target,children) in _MIXED_REFINE_TEMPLATES[linear],child in children
            push!(result,(msh=target,owner=parent.owner,
                nodes=Tuple(parent.nodes[slot] for slot in child)))
        end
    end
    return result
end

function _surface_periodic_multiset(items)
    counts=Dict{Any,Int}()
    for item in items;counts[item]=get(counts,item,0)+1;end
    return counts
end

function _surface_periodic_refine_receipt(source,class,blocks,owners,external,caller)
    snapshot=_surface_periodic_source_snapshot(source,class,caller)
    tags=class.public_tags.node_tags
    children=NamedTuple[]
    for block in blocks
        msh_dimension(block.msh)==2 || continue
        labels=get(external,block.msh,nothing)
        labels!==nothing && length(labels)==size(block.nodes,2) || throw(ArgumentError(
            "$caller: periodic child has no actual cell UID table"))
        for column in axes(block.nodes,2)
            push!(children,(msh=Int(block.msh),owner=owners[block.msh][column],uid=labels[column],
                nodes=Tuple(tags[row] for row in block.nodes[:,column])))
        end
    end
    expected=_surface_periodic_expected_children(snapshot,caller)
    actual=[(msh=cell.msh,owner=cell.owner,nodes=cell.nodes) for cell in children]
    _surface_periodic_multiset(actual)==_surface_periodic_multiset(expected) &&
        length(unique(cell.uid for cell in children))==length(children) || throw(ArgumentError(
            "$caller: actual periodic children do not match their quadratic parent templates"))
    return _SurfacePeriodicRefineReceipt(source,class,snapshot,children)
end

function _surface_periodic_refine_validate(receipt::_SurfacePeriodicRefineReceipt,
        source,target,target_class,node_tag_map,element_tag_map,caller)
    source===receipt.source || throw(ArgumentError("$caller: periodic refinement parent identity changed"))
    _surface_periodic_source_snapshot(source,receipt.source_class,caller)==receipt.source_snapshot ||
        throw(ArgumentError("$caller: actual periodic parent bindings changed"))
    target_class!==nothing && hasproperty(target_class,:mesh) && target_class.mesh===target ||
        throw(ArgumentError("$caller: periodic child classification is not bound to its actual mesh"))
    public=target_class===nothing ? nothing : target_class.public_tags
    if public===nothing && target isa MixedMesh && target.entity_data!==nothing
        # The in-flight generation/refine plan publishes its actual UID rows
        # through `entity_data` before the classification tag table is bound.
        public=_cache_public_tags(target;caller)
    end
    public!==nothing || throw(ArgumentError("$caller: periodic child UID table is missing"))
    tags=public.node_tags
    length(tags)==nnodes(target)==length(target_class.node_entities)==
        length(target_class.node_parametric) || throw(ArgumentError(
            "$caller: periodic child classification loses vertex bindings"))
    positions=Dict(uid=>row for (row,uid) in enumerate(tags))
    length(positions)==length(tags) || throw(ArgumentError("$caller: periodic child repeats a node UID"))
    function remap_uid(mapping,uid)
        mapping===nothing && return uid
        haskey(mapping,uid) || throw(ArgumentError(
            "$caller: periodic refinement UID remap is incomplete"))
        return mapping[uid]
    end
    remap_node(uid)=remap_uid(node_tag_map,uid)
    remap_cell(uid)=remap_uid(element_tag_map,uid)
    for vertex in receipt.source_snapshot.vertices
        uid=remap_node(vertex.uid);row=get(positions,uid,0)
        row>0 || throw(ArgumentError("$caller: periodic child loses a parent vertex UID"))
        raw=target_class.node_parametric[row]
        parameters=raw===nothing ? nothing : _surface_periodic_float_bits(raw)
        target_class.node_entities[row]==vertex.owner &&
            _surface_periodic_float_bits(@view(target.coords[:,row]))==vertex.xyz &&
            parameters==vertex.parameters || throw(ArgumentError(
                "$caller: periodic child changes a retained parent vertex binding"))
    end
    expected=NamedTuple[(msh=cell.msh,owner=cell.owner,uid=remap_cell(cell.uid),
        nodes=Tuple(remap_node(uid) for uid in cell.nodes)) for cell in receipt.children]
    all(cell->all(uid->haskey(positions,uid),cell.nodes),expected) || throw(ArgumentError(
        "$caller: periodic refinement loses an actual child vertex UID"))
    actual=NamedTuple[]
    for (msh,dim,offset,nodes,owners) in _cache_catalog(target,target_class)
        dim==2 || continue
        for column in axes(nodes,2)
            push!(actual,(msh=Int(msh),owner=owners[column],uid=public.element_tags[offset+column],
                nodes=Tuple(tags[row] for row in nodes[:,column])))
        end
    end
    _surface_periodic_multiset(actual)==_surface_periodic_multiset(expected) || throw(ArgumentError(
        "$caller: actual periodic child cells or current UID bindings changed"))
    return nothing
end
