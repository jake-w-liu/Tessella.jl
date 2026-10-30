struct _ModelNativeHybridSizeSeeds
    model::GeoModel
    input::_ModelMeshAuthority
    records::Vector{NamedTuple}
    surfaces::Tuple
end
_ModelNativeHybridSizeSeeds(model,input,records)=
    _ModelNativeHybridSizeSeeds(model,input,records,())

function _model_native_hybrid_global_size_seeds(m,input,caller;surface_witnesses=())
    input isa _ModelMeshAuthority && input.mesh!==nothing && input.mesh.entity_data!==nothing ||
        throw(ArgumentError("$caller: hybrid embedded sizing lacks actual current lower authority"))
    all(owner->owner[1]!=3 || haskey(m.volumes,owner[2]),keys(m.embeds)) || throw(ArgumentError(
        "$caller: hybrid embedded sizing names an undeclared target Region"))
    any(!isempty(get(m.embeds,(3,region),NTuple{2,Int}[])) for region in keys(m.volumes)) ||
        return _ModelNativeHybridSizeSeeds(m,input,NamedTuple[])
    actual=input.mesh;data=actual.entity_data
    length(data.external_node_tags)==nnodes(actual) && allunique(data.external_node_tags) &&
        all(!iszero,data.external_node_tags) && all(isfinite,actual.coords) || throw(ArgumentError(
        "$caller: hybrid embedded sizing lacks unique actual UID/XYZ records"))
    records=NamedTuple[]
    current_cells=Dict{UInt64,Tuple{Int,Int32,Tuple}}()
    length(actual.blocks)==length(data.block_entities)==length(data.external_element_tags) ||
        throw(ArgumentError("$caller: hybrid embedded sizing has incomplete current typed cell metadata"))
    for (block,owners,labels) in zip(actual.blocks,data.block_entities,data.external_element_tags)
        length(owners)==length(labels)==size(block.nodes,2) &&
            all(row->0<row<=nnodes(actual),block.nodes) || throw(ArgumentError(
            "$caller: hybrid embedded sizing has incomplete actual typed cell incidence"))
        for (column,uid) in enumerate(labels)
            uid>0 && !haskey(current_cells,uid) || throw(ArgumentError(
                "$caller: hybrid embedded sizing repeats a current global cell UID"))
            current_cells[uid]=(Int(block.msh),owners[column],Tuple(data.external_node_tags[block.nodes[:,column]]))
        end
    end
    function record(row,owner,cell_uid,size)
        0<row<=nnodes(actual) && isfinite(size) && size>0 || throw(ArgumentError(
            "$caller: hybrid embedded size has invalid current incidence or prescribed length"))
        parameter=data.node_parametric[row]
        return (owner=owner,cell_uid=cell_uid,uid=data.external_node_tags[row],
            coordinate=ntuple(axis->actual.coords[axis,row],3),
            node_owner=data.node_entities[row],parameter=parameter===nothing ? nothing : Tuple(parameter),size=Float64(size))
    end
    entities=Set{Tuple{Int,Int}}()
    for region in sort!(collect(keys(m.volumes))),entity in get(m.embeds,(3,region),NTuple{2,Int}[])
        if entity[1]==2
            entity in entities && continue
            witnesses=_ModelNativeHybridEmbeddedSurfaceInput[witness for witness in surface_witnesses
                if witness isa _ModelNativeHybridEmbeddedSurfaceInput && witness.model===m &&
                    witness.surface==entity[2]]
            length(witnesses)==1 || throw(ArgumentError(
                "$caller: hybrid embedded Face recovery requires a separate native constructor receipt"))
            push!(entities,entity)
            embedded_surface=only(witnesses)
            _model_native_hybrid_surface_source_certificate(embedded_surface,input,
                data.external_node_tags,actual.coords,caller)
            positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
            # Every actual Face triangle contributes its literal incident
            # lengths at bound lower vertices, including the shared diagonal.
            for triangle in embedded_surface.triangles
                get(current_cells,triangle.uid,nothing)==(2,Int32(entity[2]),triangle.nodes) ||
                    throw(ArgumentError("$caller: hybrid sheet sizing changes its actual Tri3 UID/owner/corners"))
                for index in 1:3
                    first=positions[triangle.nodes[index]]
                    last=positions[triangle.nodes[mod1(index+1,3)]]
                    x=actual.coords[1,first]-actual.coords[1,last]
                    y=actual.coords[2,first]-actual.coords[2,last]
                    z=actual.coords[3,first]-actual.coords[3,last]
                    edge_length=sqrt(x*x+y*y+z*z)
                    push!(records,record(first,entity,triangle.uid,edge_length))
                    push!(records,record(last,entity,triangle.uid,edge_length))
                end
            end
            continue
        end
        entity[1] in (0,1) || throw(ArgumentError(
            "$caller: hybrid embedded Face recovery requires a separate native constructor receipt"))
        entity in entities && continue
        push!(entities,entity);dim,tag=entity
        (dim==0 ? haskey(m.points,tag) : haskey(m.curves,tag)) || throw(ArgumentError(
            "$caller: hybrid embedded sizing names an undeclared lower entity"))
        product=get(input.products,(dim,Int32(tag)),nothing)
        product!==nothing || throw(ArgumentError("$caller: hybrid embedded sizing omits current lower cells"))
        product.owner==(dim,Int32(tag)) && allunique(product.rows) &&
            all(row->0<row<=nnodes(actual),product.rows) &&
            length(product.blocks)==length(product.element_tags) || throw(ArgumentError(
                "$caller: hybrid embedded sizing changes actual lower product ownership or rows"))
        for (block,labels) in zip(product.blocks,product.element_tags)
            length(labels)==size(block.nodes,2) && all(row->0<row<=length(product.rows),block.nodes) ||
                throw(ArgumentError("$caller: hybrid embedded sizing loses current constructor incidence"))
            for (column,uid) in enumerate(labels)
                observed=(Int(block.msh),Int32(tag),Tuple(data.external_node_tags[product.rows[block.nodes[:,column]]]))
                get(current_cells,uid,nothing)==observed || throw(ArgumentError(
                    "$caller: hybrid embedded product lacks exact current cell UID/family/owner/ordered-corner binding"))
            end
        end
        if dim==0
            owned=Int32[row for row in product.rows if data.node_entities[row]==(0,Int32(tag))]
            length(owned)==1 || throw(ArgumentError(
                "$caller: hybrid embedded Point lacks its unique current native vertex constructor"))
            cells=Tuple{UInt64,Int32}[]
            for (block,labels) in zip(product.blocks,product.element_tags)
                block.msh==15 && length(labels)==size(block.nodes,2) || throw(ArgumentError(
                    "$caller: hybrid embedded Point lacks current Point-cell provenance"))
                for (column,uid) in enumerate(labels)
                    push!(cells,(uid,product.rows[block.nodes[1,column]]))
                end
            end
            length(cells)==1 && only(cells)[2]==only(owned) && only(cells)[1]>0 || throw(ArgumentError(
                "$caller: hybrid embedded Point changes its exact current Point-cell binding"))
            prescribed=get(m.point_size,tag,0.0);prescribed==0 && (prescribed=GMSH_MAX_SIZE)
            push!(records,record(only(owned),entity,only(cells)[1],prescribed))
        else
            seen=Set{UInt64}();count=0
            for (block,labels) in zip(product.blocks,product.element_tags)
                block.msh==1 && size(block.nodes,1)==2 && length(labels)==size(block.nodes,2) ||
                    throw(ArgumentError("$caller: hybrid embedded Curve requires actual linear Line-cell provenance"))
                for (column,uid) in enumerate(labels)
                    uid>0 && !(uid in seen) || throw(ArgumentError(
                        "$caller: hybrid embedded Curve repeats or omits a current cell UID"))
                    push!(seen,uid);count+=1
                    first,last=product.rows[block.nodes[:,column]]
                    first!=last || throw(ArgumentError("$caller: hybrid embedded Line repeats its current node UID"))
                    x=actual.coords[1,first]-actual.coords[1,last]
                    y=actual.coords[2,first]-actual.coords[2,last]
                    z=actual.coords[3,first]-actual.coords[3,last]
                    edge_length=sqrt(x*x+y*y+z*z)
                    push!(records,record(first,entity,uid,edge_length));push!(records,record(last,entity,uid,edge_length))
                end
            end
            count>0 || throw(ArgumentError("$caller: hybrid embedded Curve has no current Line constructor"))
        end
    end
    return _ModelNativeHybridSizeSeeds(m,input,records,Tuple(collect(surface_witnesses)))
end

function _model_native_hybrid_seed_records(receipt)
    receipt isa _ModelNativeHybridSizeSeeds || throw(ArgumentError(
        "hybrid embedded sizing: missing operation-owned lower seed receipt"))
    current=_model_native_hybrid_global_size_seeds(receipt.model,receipt.input,"hybrid embedded sizing";
        surface_witnesses=receipt.surfaces)
    isequal(current.records,receipt.records) || throw(ArgumentError(
        "hybrid embedded sizing: current lower UID/cell/XYZ/UV or prescribed size changed"))
    return receipt.records
end
