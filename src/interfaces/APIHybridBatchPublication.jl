# A native batch constructs every Region before publishing any one of them.
# Its first product therefore already carries the final shared factory clocks.
# Immutable ordered snapshots bind the remaining births to that operation;
# coordinates verify recorded UIDs and never select a UID or an owner.
struct _HigherNativeHybridBatchReceipt
    batch::Model._ModelNativeHybridBatchProduction
    node_base::UInt64
    cell_base::UInt64
    groups::Tuple
    closures::Tuple
    input_nodes::Tuple
    input_cells::Tuple
    products::Tuple
end

function _higher_hybrid_input_snapshot(input,caller)
    actual=input.mesh;data=actual.entity_data
    data!==nothing && validate(actual).ok &&
        length(data.external_node_tags)==length(data.node_entities)==length(data.node_parametric)==nnodes(actual) &&
        length(actual.blocks)==length(data.external_element_tags)==length(data.block_entities) &&
        all(length(labels)==length(owners)==size(block.nodes,2)
            for (block,labels,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)) &&
        allunique(data.external_node_tags) && all(!iszero,data.external_node_tags) ||
        throw(ArgumentError("$caller: hybrid publication lacks complete actual Source UIDs"))
    nodes=Tuple((uid=data.external_node_tags[row],owner=data.node_entities[row],
        point=ntuple(axis->actual.coords[axis,row],3),
        parameter=data.node_parametric[row]===nothing ? nothing : Tuple(data.node_parametric[row]))
        for row in axes(actual.coords,2))
    cells=Tuple((uid=uid,msh=Int(block.msh),owner=(msh_dimension(block.msh),owners[column]),physical=block.tags[column],
        nodes=Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (block,labels,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)
        for (column,uid) in enumerate(labels))
    allunique(cell.uid for cell in cells) && all(cell->cell.uid>0,cells) ||
        throw(ArgumentError("$caller: hybrid publication lacks unique actual Source cell UIDs"))
    return nodes,cells
end

function _higher_hybrid_constructor_snapshot(constructor)
    return (boundaries=Tuple(constructor.boundaries),uids=Tuple(constructor.node_uids),
        coordinates=Tuple(ntuple(axis->constructor.coordinates[axis,row],3)
            for row in axes(constructor.coordinates,2)),
        triangles=Tuple(constructor.triangles),quads=Tuple(constructor.quads),
        pyramid_order=Tuple(constructor.pyramid_order))
end

function _higher_hybrid_product_snapshot(production,tag,input,m,caller)
    production.input===input && validate(production.mesh).ok || throw(ArgumentError(
        "$caller: hybrid publication has a foreign or invalid Region product"))
    current=Model._model_native_hybrid_constructor_input(input,production.constructor.boundaries,
        production.constructor.node_uids,production.constructor.coordinates,caller)
    isequal(current,production.constructor) || throw(ArgumentError(
        "$caller: hybrid publication changed its current owned constructor presentation"))
    Model._model_native_hybrid_volume_certificate(production.mesh,production.boundary_node_tags,
        production.constructor,caller)
    ids=production.boundary_node_tags
    length(ids)==nnodes(production.mesh) && allunique(ids) && all(!iszero,ids) ||
        throw(ArgumentError("$caller: hybrid publication lacks complete product node UIDs"))
    length(production.registration_rows)==length(ids) && allunique(production.registration_rows) &&
        all(row->0<row<=length(ids),production.registration_rows) || throw(ArgumentError(
        "$caller: hybrid publication loses its actual producer registration permutation"))
    sum(size(block.nodes,2) for block in production.mesh.blocks)==length(production.element_tags) ||
        throw(ArgumentError("$caller: hybrid publication cell births do not cover typed products"))
    offset=0;cells=NamedTuple[]
    for block in production.mesh.blocks
        block.msh in (4,7) || throw(ArgumentError("$caller: hybrid publication changes cell family"))
        for column in axes(block.nodes,2)
            push!(cells,(uid=production.element_tags[offset+column],msh=Int(block.msh),
                owner=(3,Int32(tag)),physical=_dim01_physical(m,(3,Int(tag))),
                nodes=Tuple(ids[block.nodes[:,column]])))
        end
        offset+=size(block.nodes,2)
    end
    return (tag=Int(tag),production=production,mesh=production.mesh,uids=Tuple(ids),
        registration_rows=Tuple(production.registration_rows),
        physical_tags=Tuple(Tuple(block.tags) for block in production.mesh.blocks),
        coordinates=Tuple(ntuple(axis->production.mesh.coords[axis,row],3)
            for row in axes(production.mesh.coords,2)),cells=Tuple(cells),
        constructor=_higher_hybrid_constructor_snapshot(production.constructor))
end

function _higher_native_hybrid_batch_receipt(assembly,m,batch,input,regions,caller)
    batch isa Model._ModelNativeHybridBatchProduction && batch.input===input &&
        batch.regions==sort!(Int.(collect(regions))) && allunique(batch.regions) &&
        Set(keys(batch.products))==Set(batch.regions) || throw(ArgumentError(
        "$caller: hybrid publication lacks its complete operation-owned Region batch"))
    nodes,cells=_higher_hybrid_input_snapshot(input,caller)
    maximum((node.uid for node in nodes);init=UInt64(0))<=assembly.node_max<=batch.node_high_water &&
        maximum((cell.uid for cell in cells);init=UInt64(0))<=assembly.cell_max<=batch.element_high_water ||
        throw(ArgumentError("$caller: hybrid publication changes its pre-operation factory clocks"))
    source=Dict(node.uid=>node for node in nodes)
    born_nodes=Set{UInt64}();born_cells=Set{UInt64}();records=NamedTuple[]
    closures=Tuple((tag,Tuple(Model._model_volume_boundary_surfaces(m,tag,caller))) for tag in batch.regions)
    for tag in batch.regions
        production=batch.products[tag]
        production.node_high_water==batch.node_high_water &&
            production.element_high_water==batch.element_high_water || throw(ArgumentError(
            "$caller: hybrid Region product loses the final operation factory clocks"))
        record=_higher_hybrid_product_snapshot(production,tag,input,m,caller)
        closure=Model._model_mesh_projection_closure(m,3,tag,caller)
        for (row,uid) in enumerate(record.uids)
            node=get(source,uid,nothing)
            if node===nothing
                assembly.node_max<uid<=batch.node_high_water && !(uid in born_nodes) ||
                    throw(ArgumentError("$caller: hybrid batch repeats or aliases a Region-owned node birth"))
                push!(born_nodes,uid)
            else
                node.owner[1]<3 && node.owner in closure && isequal(node.point,record.coordinates[row]) ||
                    throw(ArgumentError("$caller: hybrid batch consumes a foreign or changed actual lower UID"))
            end
        end
        for cell in record.cells
            assembly.cell_max<cell.uid<=batch.element_high_water && !(cell.uid in born_cells) ||
                throw(ArgumentError("$caller: hybrid batch repeats or collides with an actual cell UID"))
            push!(born_cells,cell.uid)
        end
        push!(records,record)
    end
    receipt=_HigherNativeHybridBatchReceipt(batch,assembly.node_max,assembly.cell_max,
        Tuple(Tuple(group) for group in batch.groups),closures,nodes,cells,Tuple(records))
    _higher_native_hybrid_batch_current(receipt,assembly,m,caller)
    return receipt
end

function _higher_native_hybrid_batch_current(receipt,assembly,m,caller)
    receipt isa _HigherNativeHybridBatchReceipt || throw(ArgumentError(
        "$caller: hybrid publication lacks a native operation-owned birth receipt"))
    batch=receipt.batch
    Tuple(Tuple(group) for group in batch.groups)==receipt.groups &&
        Tuple((tag,Tuple(Model._model_volume_boundary_surfaces(m,tag,caller))) for tag in batch.regions)==receipt.closures &&
        Set(keys(batch.products))==Set(record.tag for record in receipt.products) ||
        throw(ArgumentError("$caller: hybrid publication changed its selected Region topology"))
    current_nodes,current_cells=_higher_hybrid_input_snapshot(batch.input,caller)
    isequal(current_nodes,receipt.input_nodes) && isequal(current_cells,receipt.input_cells) ||
        throw(ArgumentError("$caller: hybrid publication changed actual Source XYZ, UV, owner, or ordered cells"))
    # Lower records can have only Point UV tuples and short cell supports.
    # Region births add no-parameter nodes and Tet/Pyramid support tuples.
    nodes=Dict{UInt64,NamedTuple}(node.uid=>node for node in receipt.input_nodes)
    cells=Dict{UInt64,NamedTuple}(cell.uid=>cell for cell in receipt.input_cells)
    for record in receipt.products
        production=get(batch.products,record.tag,nothing)
        production===record.production && production.mesh===record.mesh &&
            isequal(_higher_hybrid_product_snapshot(production,record.tag,batch.input,m,caller),record) ||
            throw(ArgumentError("$caller: hybrid publication changed its recorded constructor or product"))
        for (row,uid) in enumerate(record.uids)
            haskey(nodes,uid) && continue
            nodes[uid]=(uid=uid,owner=(3,Int32(record.tag)),point=record.coordinates[row],parameter=nothing)
        end
        for cell in record.cells
            haskey(cells,cell.uid) && throw(ArgumentError("$caller: hybrid publication aliases cell ownership"))
            cells[cell.uid]=cell
        end
    end
    allunique(node.tag for node in assembly.nodes) && allunique(cell.tag for cell in assembly.cells) ||
        throw(ArgumentError("$caller: hybrid publication assembly repeats an actual UID"))
    for (row,node) in enumerate(assembly.nodes)
        expected=get(nodes,node.tag,nothing)
        expected!==nothing && get(assembly.node_positions,node.tag,Int32(0))==row &&
            node.owner==expected.owner && isequal(node.point,expected.point) &&
            isequal(node.parameter===nothing ? nothing : Tuple(node.parameter),expected.parameter) ||
            throw(ArgumentError("$caller: hybrid publication assembly changed an actual node receipt"))
    end
    for cell in assembly.cells
        expected=get(cells,cell.tag,nothing)
        expected!==nothing && cell.msh==expected.msh && cell.owner==expected.owner && cell.physical==expected.physical &&
            all(row->0<row<=length(assembly.nodes),cell.nodes) &&
            Tuple(assembly.nodes[row].tag for row in cell.nodes)==expected.nodes ||
            throw(ArgumentError("$caller: hybrid publication assembly changed an ordered cell receipt"))
    end
    receipt.node_base<=assembly.node_max<=batch.node_high_water &&
        receipt.cell_base<=assembly.cell_max<=batch.element_high_water || throw(ArgumentError(
        "$caller: hybrid publication assembly lost the shared operation factory clocks"))
    return nothing
end

function _higher_native_hybrid_batch_binding(receipt,assembly,m,mesh,class,dimension,tag,
        source_row_tags,source_node_base,source_node_high_water,source_cell_tags,source_cell_high_water,caller;
        source_row_order=nothing)
    _higher_native_hybrid_batch_current(receipt,assembly,m,caller)
    records=filter(record->record.tag==tag,receipt.products)
    length(records)==1 || throw(ArgumentError("$caller: hybrid publication uses a foreign Region"))
    record=only(records)
    dimension==3 && mesh===record.mesh && class.mesh===mesh && class.entity==(3,Int32(tag)) &&
        length(class.node_entities)==length(class.node_parametric)==nnodes(mesh) &&
        source_row_tags!==nothing && source_cell_tags!==nothing &&
        source_row_order!==nothing && Tuple(source_row_order)==record.registration_rows &&
        Tuple(source_row_tags)==record.uids && source_node_base==receipt.node_base &&
        source_node_high_water==receipt.batch.node_high_water &&
        Tuple(source_cell_tags)==Tuple(cell.uid for cell in record.cells) &&
        source_cell_high_water==receipt.batch.element_high_water || throw(ArgumentError(
        "$caller: hybrid publication mixes missing, foreign, or reordered producer births"))
    lower=Dict(node.uid=>node for node in receipt.input_nodes)
    for (row,uid) in enumerate(record.uids)
        source=get(lower,uid,nothing)
        owner=source===nothing ? (3,Int32(tag)) : source.owner
        parameter=source===nothing ? nothing : source.parameter
        class.node_entities[row]==owner &&
            isequal(class.node_parametric[row]===nothing ? nothing : Tuple(class.node_parametric[row]),parameter) ||
            throw(ArgumentError("$caller: hybrid publication changes a recorded node owner or raw UV"))
    end
    for block in mesh.blocks
        get(class.cell_entities,Int(block.msh),nothing)==fill(Int32(tag),size(block.nodes,2)) ||
            throw(ArgumentError("$caller: hybrid publication changes its typed Region ownership"))
    end
    _higher_native_hybrid_batch_cell_births(receipt,assembly,tag,source_cell_tags,source_cell_high_water,caller)
    return nothing
end

function _higher_native_hybrid_batch_cell_births(receipt,assembly,tag,ids,high_water,caller)
    receipt isa _HigherNativeHybridBatchReceipt || throw(ArgumentError(
        "$caller: hybrid cell publication lacks its operation receipt"))
    records=filter(record->record.tag==tag,receipt.products)
    length(records)==1 && ids!==nothing && high_water==receipt.batch.element_high_water ||
        throw(ArgumentError("$caller: hybrid cell publication lacks a unique producer"))
    expected=Tuple(cell.uid for cell in only(records).cells)
    Tuple(ids)==expected && allunique(ids) && all(id->receipt.cell_base<id<=high_water,ids) &&
        isempty(intersect(Set(ids),Set(cell.tag for cell in assembly.cells))) || throw(ArgumentError(
        "$caller: hybrid cell publication repeats, reorders, or collides with actual births"))
    return nothing
end
