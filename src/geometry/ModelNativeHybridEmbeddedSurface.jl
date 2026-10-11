# One literal interior Tri3 sheet, with its actual Point/Line boundary products.
struct _ModelNativeHybridEmbeddedSurfaceInput
    model::GeoModel
    region::Int
    surface::Int
    input::_ModelMeshAuthority
    constructor::MixedMesh
    loops::Dict{Int,Vector{Int}}
    endpoints::Dict{Int,NTuple{2,Int}}
    triangle::NamedTuple
    records::Vector{NamedTuple}
    lower_cells::Vector{NamedTuple}
end

function _model_native_hybrid_embedded_surface_input(m,region,input,caller)
    embeds=get(m.embeds,(3,region),NTuple{2,Int}[])
    length(m.volumes)==1 && length(embeds)==1 && only(embeds)[1]==2 || return nothing
    surface=only(embeds)[2]
    haskey(m.surfaces,surface) && length(m.surfaces[surface])==1 &&
        get(m.surface_types,surface,:plane)==:plane &&
        isempty(get(m.embeds,(2,surface),NTuple{2,Int}[])) || throw(ArgumentError(
        "$caller: embedded Surface requires one declared triangle boundary without nested constraints"))
    loop=only(m.surfaces[surface]);curves=_model_projection_surface_curves(m,surface)
    length(curves)==3 && length(m.loops[loop])==3 || throw(ArgumentError(
        "$caller: embedded Surface requires three actual boundary Curves"))
    all(curve->get(m.curve_types,curve,:line)==:line,curves) || throw(ArgumentError(
        "$caller: embedded Surface requires its declared straight Line boundary"))
    product=_model_mesh_product(input,2,surface;caller)
    product isa MixedMesh && product.entity_data!==nothing && length(product.blocks)==1 &&
        only(product.blocks).msh==2 && size(only(product.blocks).nodes)==(3,1) || throw(ArgumentError(
        "$caller: embedded Surface needs its one actual linear Tri3 cell"))
    constructor=get(input.surface_constructor_inputs,(2,Int32(surface)),nothing)
    constructor!==nothing || throw(ArgumentError(
        "$caller: embedded Surface lacks its owned actual Surface constructor"))
    _model_surface_constructor_binding(product,constructor,caller)
    data=product.entity_data;block=only(product.blocks)
    labels=only(data.external_element_tags);owners=only(data.block_entities)
    length(labels)==length(owners)==1 && only(labels)>0 && only(owners)==surface ||
        throw(ArgumentError("$caller: embedded Surface changes its actual Tri3 UID/owner"))
    nodes=ntuple(slot->data.external_node_tags[block.nodes[slot,1]],3)
    allunique(nodes) && all(!iszero,nodes) || throw(ArgumentError(
        "$caller: embedded Surface has duplicate or missing actual vertex UIDs"))
    triangle=(owner=(2,Int32(surface)),uid=only(labels),nodes=nodes)
    endpoints=Dict(curve=>m.curves[curve] for curve in curves)
    points=sort!(unique!(Int[point for pair in values(endpoints) for point in pair]))
    length(points)==3 || throw(ArgumentError(
        "$caller: embedded Surface requires three distinct declared endpoint Points"))
    point_uids=Dict{Int,UInt64}();lower_cells=NamedTuple[]
    for point in points
        actual=_model_mesh_product(input,0,point;caller)
        actual isa MixedMesh && actual.entity_data!==nothing && length(actual.blocks)==1 &&
            only(actual.blocks).msh==15 && size(only(actual.blocks).nodes)==(1,1) || throw(ArgumentError(
            "$caller: embedded Surface endpoint lacks its actual Point cell"))
        pdata=actual.entity_data;row=only(only(actual.blocks).nodes)
        labels=only(pdata.external_element_tags)
        length(labels)==1 && only(labels)>0 && pdata.node_entities[row]==(0,Int32(point)) ||
            throw(ArgumentError("$caller: embedded Surface endpoint changes its Point UID/owner"))
        uid=pdata.external_node_tags[row];point_uids[point]=uid
        push!(lower_cells,(owner=(0,Int32(point)),uid=only(labels),nodes=(uid,)))
    end
    Set(values(point_uids))==Set(nodes) || throw(ArgumentError(
        "$caller: embedded Surface vertices omit or alias declared endpoint Point UIDs"))
    expected_edges=Set(minmax(nodes[index],nodes[mod1(index+1,3)]) for index in 1:3)
    actual_edges=Set{NTuple{2,UInt64}}()
    for curve in curves
        actual=_model_mesh_product(input,1,curve;caller)
        actual isa MixedMesh && actual.entity_data!==nothing && length(actual.blocks)==1 &&
            only(actual.blocks).msh==1 && size(only(actual.blocks).nodes)==(2,1) || throw(ArgumentError(
            "$caller: embedded Surface boundary lacks its actual Line2 cell"))
        cdata=actual.entity_data;labels=only(cdata.external_element_tags)
        line=ntuple(slot->cdata.external_node_tags[only(actual.blocks).nodes[slot,1]],2)
        length(labels)==1 && only(labels)>0 && only(only(cdata.block_entities))==curve &&
            minmax(line...)==minmax((point_uids[point] for point in endpoints[curve])...) &&
            !(minmax(line...) in actual_edges) || throw(ArgumentError(
            "$caller: embedded Surface boundary changes actual Line2 UID/owner/incidence"))
        push!(actual_edges,minmax(line...))
        push!(lower_cells,(owner=(1,Int32(curve)),uid=only(labels),nodes=line))
    end
    actual_edges==expected_edges || throw(ArgumentError(
        "$caller: embedded Surface triangle does not match its actual boundary Line2 cells"))
    positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    records=NamedTuple[]
    for point in points
        uid=point_uids[point];row=positions[uid]
        data.node_entities[row]==(0,Int32(point)) || throw(ArgumentError(
            "$caller: embedded Surface loses a current lower node owner"))
        push!(records,(uid=uid,coordinate=Tuple(product.coords[:,row]),node_owner=data.node_entities[row]))
    end
    return _ModelNativeHybridEmbeddedSurfaceInput(m,region,surface,input,constructor,
        Dict(loop=>copy(m.loops[loop])),endpoints,triangle,records,lower_cells)
end

function _model_native_hybrid_surface_source_certificate(witness,input,uids,coordinates,caller)
    witness===nothing && return nothing
    witness isa _ModelNativeHybridEmbeddedSurfaceInput && witness.input===input &&
        get(witness.model.embeds,(3,witness.region),NTuple{2,Int}[])==[(2,witness.surface)] ||
        throw(ArgumentError("$caller: embedded Surface has a foreign or changed constructor"))
    fresh=_model_native_hybrid_embedded_surface_input(witness.model,witness.region,input,caller)
    fresh!==nothing && fresh.constructor===witness.constructor && fresh.loops==witness.loops &&
        fresh.endpoints==witness.endpoints && fresh.triangle==witness.triangle &&
        fresh.records==witness.records && fresh.lower_cells==witness.lower_cells || throw(ArgumentError(
        "$caller: embedded Surface changes its actual current UID/owner/cell witness"))
    positions=Dict(uid=>row for (row,uid) in enumerate(uids))
    for record in witness.records
        row=get(positions,record.uid,0)
        row>0 && isequal(Tuple(coordinates[:,row]),record.coordinate) || throw(ArgumentError(
            "$caller: embedded Surface UID/XYZ is missing or changed in the native Source"))
    end
    return nothing
end

function _model_native_hybrid_surface_containment_certificate(constructor,caller)
    witness=constructor.embedded_surface;witness===nothing && return nothing
    _model_native_hybrid_surface_source_certificate(witness,constructor.input,
        constructor.node_uids,constructor.coordinates,caller)
    positions=Dict(uid=>row for (row,uid) in enumerate(constructor.node_uids))
    panels=vcat(constructor.triangles,constructor.quads)
    boundary=Set(uid for panel in panels for uid in panel.nodes)
    coordinate(uid)=Tuple(constructor.coordinates[:,positions[uid]])
    all(record->!(record.uid in boundary),witness.records) || throw(ArgumentError(
        "$caller: embedded Surface must be strictly inside the shell"))
    for panel in panels
        triangles=length(panel.nodes)==3 ? (panel.nodes,) :
            ((panel.nodes[1],panel.nodes[2],panel.nodes[3]),(panel.nodes[1],panel.nodes[3],panel.nodes[4]))
        for triangle in triangles
            a,b,c=map(coordinate,triangle)
            sides=[orient3(a,b,c,coordinate(uid)) for uid in boundary]
            all(isfinite,sides) || throw(ArgumentError("$caller: embedded Surface shell halfspaces are not finite"))
            signs=Set(sign(side) for side in sides if side!=0)
            length(signs)==1 || throw(ArgumentError("$caller: embedded Surface needs a certified convex shell"))
            all(record->begin
                side=orient3(a,b,c,record.coordinate)
                isfinite(side) && side!=0 && sign(side)==only(signs)
            end,witness.records) || throw(ArgumentError(
                "$caller: embedded Surface is outside or on the shell boundary"))
        end
    end
    return nothing
end

function _model_native_hybrid_surface_constraints(constructor,uids,caller)
    witness=constructor.embedded_surface
    witness===nothing && return (facets=NTuple{3,Int32}[],edges=NTuple{2,Int32}[])
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(uids))
    all(record->haskey(positions,record.uid),witness.records) || throw(ArgumentError(
        "$caller: embedded Surface constraint omits an actual node UID"))
    face=ntuple(slot->positions[witness.triangle.nodes[slot]],3)
    edges=NTuple{2,Int32}[ntuple(slot->positions[cell.nodes[slot]],2)
                         for cell in witness.lower_cells if cell.owner[1]==1]
    return (facets=[face],edges=edges)
end

function _model_native_hybrid_embedded_surface_certificate(mesh,uids,constructor,caller)
    witness=constructor.embedded_surface;witness===nothing && return nothing
    _model_native_hybrid_surface_source_certificate(witness,constructor.input,uids,mesh.coords,caller)
    _model_native_hybrid_surface_containment_certificate(constructor,caller)
    constraints=_model_native_hybrid_surface_constraints(constructor,uids,caller)
    face=only(constraints.facets);sides=Float64[]
    for row in face
        _model_native_hybrid_interior_node_star_certificate(mesh,row,caller,"embedded Surface vertex")
    end
    for block in mesh.blocks,column in axes(block.nodes,2)
        nodes=block.nodes[:,column];all(row->row in nodes,face) || continue
        block.msh==4 || throw(ArgumentError("$caller: embedded Surface facet has a foreign typed cell role"))
        opposite=only(row for row in nodes if !(row in face))
        push!(sides,orient3((Tuple(mesh.coords[:,row]) for row in (face...,opposite))...))
    end
    length(sides)==2 && all(isfinite,sides) && Set(sign.(sides))==Set((-1.0,1.0)) ||
        throw(ArgumentError("$caller: embedded Surface facet lacks two-sided positive Tet incidence"))
    for edge in constraints.edges
        links=Set{NTuple{2,Int32}}();degrees=Dict{Int32,Int}()
        for block in mesh.blocks,column in axes(block.nodes,2)
            nodes=block.nodes[:,column];all(row->row in nodes,edge) || continue
            block.msh==4 || throw(ArgumentError("$caller: embedded Surface edge has a foreign typed cell role"))
            link=minmax((row for row in nodes if !(row in edge))...)
            link in links && throw(ArgumentError("$caller: embedded Surface repeats a Tet edge-link panel"))
            push!(links,link)
            for row in link;degrees[row]=get(degrees,row,0)+1;end
        end
        length(links)>=3 && all(==(2),values(degrees)) || throw(ArgumentError(
            "$caller: embedded Surface edge lacks its closed interior Tet link"))
        reached=Set([first(first(links))]);changed=true
        while changed
            changed=false
            for link in links
                any(row->row in reached,link) || continue
                for row in link
                    row in reached || (push!(reached,row);changed=true)
                end
            end
        end
        reached==Set(keys(degrees)) || throw(ArgumentError(
            "$caller: embedded Surface has disconnected Tet edge links"))
    end
    return nothing
end
