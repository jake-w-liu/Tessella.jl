# A literal one- or two-Tri3 planar sheet with actual Point/Line boundary products.
struct _ModelNativeHybridEmbeddedSurfaceInput
    model::GeoModel
    region::Int
    surface::Int
    input::_ModelMeshAuthority
    constructor::MixedMesh
    loops::Dict{Int,Vector{Int}}
    endpoints::Dict{Int,NTuple{2,Int}}
    triangles::Vector{NamedTuple}
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
        "$caller: embedded Surface requires one declared planar boundary without nested constraints"))
    loop=only(m.surfaces[surface]);curves=_model_projection_surface_curves(m,surface)
    length(curves) in (3,4) && length(m.loops[loop])==length(curves) || throw(ArgumentError(
        "$caller: embedded Surface requires three or four actual boundary Curves"))
    all(curve->get(m.curve_types,curve,:line)==:line,curves) || throw(ArgumentError(
        "$caller: embedded Surface requires its declared straight Line boundary"))
    product=_model_mesh_product(input,2,surface;caller)
    product isa MixedMesh && product.entity_data!==nothing && length(product.blocks)==1 &&
        only(product.blocks).msh==2 && size(only(product.blocks).nodes)==(3,length(curves)-2) || throw(ArgumentError(
        "$caller: embedded Surface needs its actual one- or two-cell linear Tri3 sheet"))
    constructor=get(input.surface_constructor_inputs,(2,Int32(surface)),nothing)
    constructor!==nothing || throw(ArgumentError(
        "$caller: embedded Surface lacks its owned actual Surface constructor"))
    _model_surface_constructor_binding(product,constructor,caller)
    data=product.entity_data;block=only(product.blocks)
    labels=only(data.external_element_tags);owners=only(data.block_entities)
    length(labels)==length(owners)==size(block.nodes,2) && all(>(0),labels) &&
        allunique(labels) && all(==(surface),owners) || throw(ArgumentError(
        "$caller: embedded Surface changes its actual Tri3 UID/owner"))
    triangles=NamedTuple[];used=Set{UInt64}()
    for (column,uid) in enumerate(labels)
        nodes=ntuple(slot->data.external_node_tags[block.nodes[slot,column]],3)
        allunique(nodes) && all(!iszero,nodes) || throw(ArgumentError(
            "$caller: embedded Surface has duplicate or missing actual vertex UIDs"))
        union!(used,nodes)
        push!(triangles,(owner=(2,Int32(surface)),uid=uid,nodes=nodes))
    end
    endpoints=Dict(curve=>m.curves[curve] for curve in curves)
    points=sort!(unique!(Int[point for pair in values(endpoints) for point in pair]))
    length(points)==length(curves) || throw(ArgumentError(
        "$caller: embedded Surface requires distinct declared boundary endpoint Points"))
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
    Set(values(point_uids))==used || throw(ArgumentError(
        "$caller: embedded Surface vertices omit or alias declared endpoint Point UIDs"))
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
    positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    records=NamedTuple[]
    for point in points
        uid=point_uids[point];row=positions[uid]
        data.node_entities[row]==(0,Int32(point)) || throw(ArgumentError(
            "$caller: embedded Surface loses a current lower node owner"))
        push!(records,(uid=uid,coordinate=Tuple(product.coords[:,row]),node_owner=data.node_entities[row]))
    end
    _model_native_hybrid_surface_disk_certificate(triangles,records,actual_edges,caller)
    return _ModelNativeHybridEmbeddedSurfaceInput(m,region,surface,input,constructor,
        Dict(loop=>copy(m.loops[loop])),endpoints,triangles,records,lower_cells)
end

# The initial extension admits a single triangle or two coherently oriented
# triangles meeting at one actual edge. Their literal outer edges must equal
# the complete lower Line-cell boundary; no triangulation is reconstructed.
function _model_native_hybrid_surface_disk_certificate(triangles,records,boundary,caller)
    incidence=Dict{NTuple{2,UInt64},Vector{NTuple{2,UInt64}}}()
    supports=Set{NTuple{3,UInt64}}()
    for triangle in triangles
        support=Tuple(sort!(collect(triangle.nodes)))
        support in supports && throw(ArgumentError("$caller: embedded Surface repeats a literal Tri3 support"))
        push!(supports,support)
        for index in 1:3
            edge=(triangle.nodes[index],triangle.nodes[mod1(index+1,3)])
            push!(get!(incidence,minmax(edge...),NTuple{2,UInt64}[]),edge)
        end
    end
    all(edges->length(edges) in (1,2),values(incidence)) &&
        all(edges->length(edges)==1 || first(edges)==reverse(last(edges)),values(incidence)) ||
        throw(ArgumentError("$caller: embedded Surface has incoherent or nonmanifold actual triangle edges"))
    observed=Set(edge for (edge,entries) in incidence if length(entries)==1)
    observed==boundary && count(entries->length(entries)==2,values(incidence))==length(triangles)-1 &&
        length(records)-length(incidence)+length(triangles)==1 || throw(ArgumentError(
        "$caller: embedded Surface triangles do not form their complete actual Line-boundary disk"))
    coordinates=Dict(record.uid=>record.coordinate for record in records)
    a,b,c=map(uid->coordinates[uid],first(triangles).nodes)
    normal=LinearAlgebra.cross(collect(b).-collect(a),collect(c).-collect(a))
    magnitude=LinearAlgebra.norm(normal)
    isfinite(magnitude) && magnitude>0 || throw(ArgumentError(
        "$caller: embedded Surface has a degenerate actual triangle"))
    scale=max(1.0,maximum(LinearAlgebra.norm(collect(record.coordinate).-collect(a)) for record in records))
    bound=128eps(Float64)*scale*magnitude
    # orient3 is an exact sign predicate; plane distance needs the adaptive
    # determinant magnitude on the same literal stored coordinates.
    determinants=[parentmodule(@__MODULE__).Predicates._orient3_value(a,b,c,record.coordinate)
                  for record in records]
    isfinite(bound) && all(value->isfinite(value) && abs(value)<=bound,determinants) ||
        throw(ArgumentError("$caller: embedded Surface is not a certified planar sheet"))
    drop=argmax(abs.(normal));axes=filter(!=(drop),collect(1:3))
    project(uid)=ntuple(index->coordinates[uid][axes[index]],2)
    allunique(project(record.uid) for record in records) || throw(ArgumentError(
        "$caller: embedded Surface repeats a projected vertex"))
    areas=[orient2((project(uid) for uid in triangle.nodes)...) for triangle in triangles]
    all(isfinite,areas) && all(!iszero,areas) && length(Set(sign.(areas)))==1 ||
        throw(ArgumentError("$caller: embedded Surface has degenerate or opposing actual triangle winding"))
    if length(triangles)==2
        diagonal=only(edge for (edge,entries) in incidence if length(entries)==2)
        opposite=[only(uid for uid in triangle.nodes if !(uid in diagonal)) for triangle in triangles]
        sides=[orient2(project(diagonal[1]),project(diagonal[2]),project(uid)) for uid in opposite]
        all(isfinite,sides) && Set(sign.(sides))==Set((-1.0,1.0)) || throw(ArgumentError(
            "$caller: embedded Surface triangles overlap or fold across their actual shared edge"))
    end
    return nothing
end

function _model_native_hybrid_surface_source_certificate(witness,input,uids,coordinates,caller)
    witness===nothing && return nothing
    witness isa _ModelNativeHybridEmbeddedSurfaceInput && witness.input===input &&
        get(witness.model.embeds,(3,witness.region),NTuple{2,Int}[])==[(2,witness.surface)] ||
        throw(ArgumentError("$caller: embedded Surface has a foreign or changed constructor"))
    fresh=_model_native_hybrid_embedded_surface_input(witness.model,witness.region,input,caller)
    fresh!==nothing && fresh.constructor===witness.constructor && fresh.loops==witness.loops &&
        fresh.endpoints==witness.endpoints && fresh.triangles==witness.triangles &&
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
    facets=NTuple{3,Int32}[ntuple(slot->positions[triangle.nodes[slot]],3)
                          for triangle in witness.triangles]
    # Every literal sheet edge, including the shared diagonal, stays protected.
    edges=NTuple{2,Int32}[ntuple(slot->positions[cell.nodes[slot]],2)
                         for cell in witness.lower_cells if cell.owner[1]==1]
    seen=Set(minmax(edge...) for edge in edges)
    for face in facets,index in 1:3
        edge=(face[index],face[mod1(index+1,3)]);support=minmax(edge...)
        support in seen && continue
        push!(seen,support);push!(edges,edge)
    end
    return (;facets,edges)
end

function _model_native_hybrid_embedded_surface_certificate(mesh,uids,constructor,caller)
    witness=constructor.embedded_surface;witness===nothing && return nothing
    _model_native_hybrid_surface_source_certificate(witness,constructor.input,uids,mesh.coords,caller)
    _model_native_hybrid_surface_containment_certificate(constructor,caller)
    constraints=_model_native_hybrid_surface_constraints(constructor,uids,caller)
    for row in unique!(Int32[row for face in constraints.facets for row in face])
        _model_native_hybrid_interior_node_star_certificate(mesh,row,caller,"embedded Surface vertex")
    end
    for face in constraints.facets
        sides=Float64[]
        for block in mesh.blocks,column in axes(block.nodes,2)
            nodes=block.nodes[:,column];all(row->row in nodes,face) || continue
            block.msh==4 || throw(ArgumentError("$caller: embedded Surface facet has a foreign typed cell role"))
            opposite=only(row for row in nodes if !(row in face))
            push!(sides,orient3((Tuple(mesh.coords[:,row]) for row in (face...,opposite))...))
        end
        length(sides)==2 && all(isfinite,sides) && Set(sign.(sides))==Set((-1.0,1.0)) ||
            throw(ArgumentError("$caller: embedded Surface facet lacks two-sided positive Tet incidence"))
    end
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
