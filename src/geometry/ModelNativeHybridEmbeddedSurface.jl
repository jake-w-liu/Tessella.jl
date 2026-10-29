# A literal Tri3 sheet with actual Point/Line boundary products, optional
# interior loops, subdivision chains, interior vertices and nested declared
# Point/Curve constraints. Planar sheets certify as a dominant-axis planar
# disk; warped ruled patches certify the same embedded-disk incidence in the
# declared (u,v) parameter domain.
struct _ModelNativeHybridEmbeddedSurfaceInput
    model::GeoModel
    region::Int
    surface::Int
    input::_ModelMeshAuthority
    constructor::Union{Nothing,MixedMesh}
    loops::Dict{Int,Vector{Int}}
    endpoints::Dict{Int,NTuple{2,Int}}
    triangles::Vector{NamedTuple}
    records::Vector{NamedTuple}
    lower_cells::Vector{NamedTuple}
    inventory::Vector{NTuple{2,Int}}
end

# The unique actual Point-cell witness of one declared vertex entity.
function _model_native_hybrid_embedded_sheet_point(input,point,caller)
    actual=_model_mesh_product(input,0,point;caller)
    actual isa MixedMesh && actual.entity_data!==nothing && length(actual.blocks)==1 &&
        only(actual.blocks).msh==15 && size(only(actual.blocks).nodes)==(1,1) || throw(ArgumentError(
        "$caller: embedded Surface vertex lacks its actual Point cell"))
    pdata=actual.entity_data;row=only(only(actual.blocks).nodes)
    labels=only(pdata.external_element_tags)
    length(labels)==1 && only(labels)>0 && pdata.node_entities[row]==(0,Int32(point)) ||
        throw(ArgumentError("$caller: embedded Surface vertex changes its Point UID/owner"))
    uid=pdata.external_node_tags[row]
    return uid,(owner=(0,Int32(point)),uid=only(labels),nodes=(uid,))
end

# The complete actual Line2 cell chain of one declared Curve.
function _model_native_hybrid_embedded_sheet_curve(input,curve,caller)
    product=_model_mesh_product(input,1,curve;caller)
    product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
        "$caller: embedded Surface boundary lacks its actual Curve discretization"))
    data=product.entity_data;cells=NamedTuple[];seen=Set{UInt64}()
    supports=Set{Tuple{UInt64,UInt64}}()
    for (block,labels,owners) in zip(product.blocks,data.external_element_tags,data.block_entities)
        block.msh==1 && size(block.nodes,1)==2 && length(labels)==length(owners)==size(block.nodes,2) ||
            throw(ArgumentError("$caller: embedded Surface boundary needs actual linear Line2 cells"))
        for (column,uid) in enumerate(labels)
            nodes=ntuple(slot->data.external_node_tags[block.nodes[slot,column]],2)
            uid>0 && !(uid in seen) && owners[column]==curve && all(!iszero,nodes) && allunique(nodes) ||
                throw(ArgumentError("$caller: embedded Surface boundary repeats or changes a current cell UID/owner"))
            support=minmax(nodes...)
            support in supports && throw(ArgumentError("$caller: embedded Surface boundary repeats a current segment"))
            push!(seen,uid);push!(supports,support)
            push!(cells,(owner=(1,Int32(curve)),uid=uid,nodes=nodes))
        end
    end
    isempty(cells) && throw(ArgumentError("$caller: embedded Surface boundary has no actual segment cells"))
    return cells
end

# The literal cells form one connected chain ending at the declared endpoint
# UIDs; subdivision nodes keep degree two like the embedded-Curve inventory.
function _model_native_hybrid_embedded_sheet_chain(cells,endpoint_uids,caller,label)
    used=Set{UInt64}()
    for cell in cells,uid in cell.nodes;push!(used,uid);end
    degrees=Dict(uid=>0 for uid in used)
    for cell in cells,uid in cell.nodes;degrees[uid]+=1;end
    expected_ends=first(endpoint_uids)==last(endpoint_uids) ? Set{UInt64}() : Set(endpoint_uids)
    Set(uid for (uid,degree) in degrees if degree==1)==expected_ends &&
        all(degree->degree in (1,2),values(degrees)) && all(uid->uid in used,endpoint_uids) ||
        throw(ArgumentError("$caller: $label cells do not form their declared endpoint chain"))
    reached=Set([first(endpoint_uids)]);changed=true
    while changed
        changed=false
        for cell in cells
            any(uid->uid in reached,cell.nodes) || continue
            for uid in cell.nodes
                uid in reached || (push!(reached,uid);changed=true)
            end
        end
    end
    reached==used || throw(ArgumentError("$caller: $label has disconnected actual segments"))
    return used
end

function _model_native_hybrid_embedded_surface_input(m,region,input,caller;surface=nothing)
    embeds=copy(get(m.embeds,(3,region),NTuple{2,Int}[]))
    if surface===nothing
        # The original single-entity route remains separate.
        haskey(m.volumes,region) && length(embeds)==1 && only(embeds)[1]==2 || return nothing
        surface=only(embeds)[2]
    else
        haskey(m.volumes,region) && (2,surface) in embeds &&
            all(owner->owner[1] in (0,1,2),embeds) || throw(ArgumentError(
            "$caller: selected embedded Surface lacks its actual Region declaration"))
    end
    # Planar sheets certify in a dominant-axis projection; genuinely warped
    # boundary-loop patches (`:ruled`, the upstream `RuledSurface` GFace)
    # certify in the patch's own (u,v) domain instead — the same literal
    # Tri3 cells and lower products feed the kernel either way.
    haskey(m.surfaces,surface) && !isempty(m.surfaces[surface]) &&
        get(m.surface_types,surface,:plane) in (:plane,:ruled) || throw(ArgumentError(
        "$caller: embedded Surface requires a declared planar or ruled boundary loop sheet"))
    loops=m.surfaces[surface]
    nested=copy(get(m.embeds,(2,surface),NTuple{2,Int}[]))
    allunique(nested) && all(owner->owner[1] in (0,1),nested) || throw(ArgumentError(
        "$caller: embedded Surface supports only declared nested Point/Curve constraints"))
    nested_points=Int[tag for (dim,tag) in nested if dim==0]
    nested_curves=Int[tag for (dim,tag) in nested if dim==1]
    curves=_model_projection_surface_curves(m,surface)
    declared=Int[abs(curve) for loop in loops for curve in m.loops[loop]]
    all(loop->length(m.loops[loop])>=3,loops) && length(declared)==length(curves) &&
        allunique(declared) || throw(ArgumentError(
        "$caller: embedded Surface requires three or more distinct boundary Curves per loop"))
    # Any declared Curve whose actual product is a linear Line2 chain is a
    # legal sheet boundary: nonlinear Circles/Splines discretize to the same
    # cell family upstream inserts into the volume triangulation.
    all(curve->!(curve in curves),nested_curves) || throw(ArgumentError(
        "$caller: embedded Surface names a boundary Curve as a nested constraint"))
    product=_model_mesh_product(input,2,surface;caller)
    product isa MixedMesh && product.entity_data!==nothing && length(product.blocks)==1 &&
        only(product.blocks).msh==2 && size(only(product.blocks).nodes,1)==3 &&
        !isempty(only(product.blocks).nodes) || throw(ArgumentError(
        "$caller: embedded Surface needs its actual linear Tri3 sheet"))
    # A planar sheet binds its unlifted PSLG constructor input; a warped
    # ruled patch has no such registered input — its evidence is the lifted
    # actual product plus the declared (u,v) reparameterization below.
    constructor=get(input.surface_constructor_inputs,(2,Int32(surface)),nothing)
    if constructor===nothing
        get(m.surface_types,surface,:plane)==:ruled || throw(ArgumentError(
            "$caller: embedded Surface lacks its owned actual Surface constructor"))
    else
        _model_surface_constructor_binding(product,constructor,caller)
    end
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
    Set(data.external_node_tags)==used || throw(ArgumentError(
        "$caller: embedded Surface has an actual vertex outside its Tri3 cells"))
    endpoints=Dict(curve=>m.curves[curve] for curve in curves)
    points=sort!(unique!(Int[point for pair in values(endpoints) for point in pair]))
    length(points)==length(declared) || throw(ArgumentError(
        "$caller: embedded Surface requires distinct declared boundary endpoint Points"))
    nested_endpoints=sort!(unique!(Int[point for curve in nested_curves
                                       for point in m.curves[curve]]))
    point_uids=Dict{Int,UInt64}();lower_cells=NamedTuple[]
    for point in unique!(sort!(vcat(points,nested_points,nested_endpoints)))
        uid,cell=_model_native_hybrid_embedded_sheet_point(input,point,caller)
        point_uids[point]=uid;push!(lower_cells,cell)
    end
    sheet_edges=Set(minmax(triangle.nodes[index],triangle.nodes[mod1(index+1,3)])
                    for triangle in triangles for index in 1:3)
    actual_edges=Set{NTuple{2,UInt64}}()
    for curve in curves
        cells=_model_native_hybrid_embedded_sheet_curve(input,curve,caller)
        expected=Tuple(point_uids[point] for point in endpoints[curve])
        _model_native_hybrid_embedded_sheet_chain(cells,expected,caller,"embedded Surface boundary")
        for cell in cells
            support=minmax(cell.nodes...)
            support in actual_edges && throw(ArgumentError(
                "$caller: embedded Surface repeats an actual boundary segment"))
            support in sheet_edges || throw(ArgumentError(
                "$caller: embedded Surface drops an actual boundary segment"))
            push!(actual_edges,support);push!(lower_cells,cell)
        end
    end
    for curve in nested_curves
        cells=_model_native_hybrid_embedded_sheet_curve(input,curve,caller)
        expected=Tuple(point_uids[point] for point in m.curves[curve])
        _model_native_hybrid_embedded_sheet_chain(cells,expected,caller,"embedded Surface nested Curve")
        for cell in cells
            minmax(cell.nodes...) in sheet_edges || throw(ArgumentError(
                "$caller: embedded Surface drops a nested Curve segment"))
            push!(lower_cells,cell)
        end
    end
    Set(values(point_uids)) ⊆ used || throw(ArgumentError(
        "$caller: embedded Surface drops a declared Point vertex"))
    positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    point_set=Set(points);nested_point_set=Set(nested_points)
    nested_endpoint_set=Set(nested_endpoints);curve_set=Set(curves);nested_curve_set=Set(nested_curves)
    records=NamedTuple[]
    for uid in sort!(collect(used))
        row=positions[uid];owner=data.node_entities[row]
        admissible=owner==(2,Int32(surface)) ||
            (owner[1]==0 && (Int(owner[2]) in point_set || Int(owner[2]) in nested_point_set ||
                             Int(owner[2]) in nested_endpoint_set)) ||
            (owner[1]==1 && (Int(owner[2]) in curve_set || Int(owner[2]) in nested_curve_set))
        admissible || throw(ArgumentError(
            "$caller: embedded Surface vertex $uid keeps a foreign lower owner $owner"))
        for (point,point_uid) in point_uids
            point_uid==uid && owner!=(0,Int32(point)) && throw(ArgumentError(
                "$caller: embedded Surface loses a declared Point vertex owner"))
        end
        push!(records,(uid=uid,coordinate=Tuple(product.coords[:,row]),node_owner=owner))
    end
    parameters=nothing
    if get(m.surface_types,surface,:plane)==:ruled
        # Upstream's `RuledSurface` GFace triangulates inside `parBounds`, so
        # the literal lifted facets form their embedded disk in the patch's
        # own (u,v) domain. Each stored vertex is reparameterized through the
        # same `XYZtoUV` solve the dim-2 ruled mesher uses.
        silent=Returns(nothing)
        span=maximum(LinearAlgebra.norm(collect(record.coordinate).-
            collect(first(records).coordinate)) for record in records)
        lc=max(1.0,span)
        parameters=Dict{UInt64,NTuple{2,Float64}}()
        for record in records
            u,v=_ruled_xyz_to_uv(m,surface,record.coordinate,1.0,lc,false,
                                 silent,silent,caller)
            (isfinite(u) && isfinite(v)) || throw(ArgumentError(
                "$caller: embedded Surface vertex has no finite ruled parameter"))
            parameters[record.uid]=(u,v)
        end
    end
    _model_native_hybrid_surface_disk_certificate(triangles,records,actual_edges,
        length(loops),caller;parameters)
    return _ModelNativeHybridEmbeddedSurfaceInput(m,region,surface,input,constructor,
        Dict(loop=>copy(m.loops[loop]) for loop in loops),endpoints,triangles,records,
        lower_cells,embeds)
end

# The literal boundary must be `cycles` connected cycles and the triangle dual
# graph connected. Exact projected incidence checks certify adjacent and
# nonadjacent triangles; no triangulation, source coordinates or constructor
# order is reconstructed. `parameters` supplies the certifying 2D embedding:
# nothing means the sheet is planar and is certified in its dominant-axis
# projection (the planarity block below also runs), while a ruled sheet
# supplies each vertex's own (u,v) parameter and skips only that block.
function _model_native_hybrid_surface_disk_certificate(triangles,records,boundary,cycles,caller;
                                                       parameters=nothing)
    incidence=Dict{NTuple{2,UInt64},Vector{NTuple{2,UInt64}}}()
    adjacent=Dict{NTuple{2,UInt64},Vector{Int}}()
    supports=Set{NTuple{3,UInt64}}()
    for (column,triangle) in enumerate(triangles)
        support=Tuple(sort!(collect(triangle.nodes)))
        support in supports && throw(ArgumentError("$caller: embedded Surface repeats a literal Tri3 support"))
        push!(supports,support)
        for index in 1:3
            edge=(triangle.nodes[index],triangle.nodes[mod1(index+1,3)])
            push!(get!(incidence,minmax(edge...),NTuple{2,UInt64}[]),edge)
            push!(get!(adjacent,minmax(edge...),Int[]),column)
        end
    end
    all(edges->length(edges) in (1,2),values(incidence)) &&
        all(edges->length(edges)==1 || first(edges)==reverse(last(edges)),values(incidence)) ||
        throw(ArgumentError("$caller: embedded Surface has incoherent or nonmanifold actual triangle edges"))
    observed=Set(edge for (edge,entries) in incidence if length(entries)==1)
    # A connected planar sheet with `cycles` boundary loops has Euler value 2-cycles.
    observed==boundary && length(records)-length(incidence)+length(triangles)==2-cycles ||
        throw(ArgumentError(
        "$caller: embedded Surface triangles do not form their complete actual Line-boundary sheet"))
    boundary_graph=Dict{UInt64,Vector{UInt64}}()
    for (vertex_a,vertex_b) in boundary
        push!(get!(boundary_graph,vertex_a,UInt64[]),vertex_b)
        push!(get!(boundary_graph,vertex_b,UInt64[]),vertex_a)
    end
    all(neighbors->length(neighbors)==2,values(boundary_graph)) || throw(ArgumentError(
        "$caller: embedded Surface lower boundary is not an endpoint cycle set"))
    reached=Set{UInt64}();components=0
    for start in keys(boundary_graph)
        start in reached && continue
        components+=1;pending=UInt64[start]
        while !isempty(pending)
            vertex=pop!(pending);vertex in reached && continue
            push!(reached,vertex);append!(pending,boundary_graph[vertex])
        end
    end
    components==cycles || throw(ArgumentError(
        "$caller: embedded Surface lower boundary loses or merges declared cycles"))
    triangle_graph=[Int[] for _ in triangles]
    for columns in values(adjacent)
        length(columns)==2 || continue
        triangle_a,triangle_b=columns
        push!(triangle_graph[triangle_a],triangle_b);push!(triangle_graph[triangle_b],triangle_a)
    end
    pending_triangles=Int[1];reached_triangles=Set{Int}()
    while !isempty(pending_triangles)
        column=pop!(pending_triangles);column in reached_triangles && continue
        push!(reached_triangles,column);append!(pending_triangles,triangle_graph[column])
    end
    length(reached_triangles)==length(triangles) || throw(ArgumentError(
        "$caller: embedded Surface has disconnected actual triangles"))
    coordinates=Dict(record.uid=>record.coordinate for record in records)
    project=if parameters===nothing
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
        uid->ntuple(index->coordinates[uid][axes[index]],2)
    else
        # A warped patch keeps its literal 3D lift; the declared (u,v)
        # domain is where the embedded-disk incidence certificates hold.
        uid->parameters[uid]
    end
    allunique(project(record.uid) for record in records) || throw(ArgumentError(
        "$caller: embedded Surface repeats a projected vertex"))
    areas=[orient2((project(uid) for uid in triangle.nodes)...) for triangle in triangles]
    all(isfinite,areas) && all(!iszero,areas) && length(Set(sign.(areas)))==1 ||
        throw(ArgumentError("$caller: embedded Surface has degenerate or opposing actual triangle winding"))
    for (diagonal,columns) in adjacent
        length(columns)==2 || continue
        opposite=[only(uid for uid in triangles[column].nodes if !(uid in diagonal)) for column in columns]
        sides=[orient2(project(diagonal[1]),project(diagonal[2]),project(uid)) for uid in opposite]
        all(isfinite,sides) && Set(sign.(sides))==Set((-1.0,1.0)) || throw(ArgumentError(
            "$caller: embedded Surface triangles overlap or fold across their actual shared edge"))
    end
    on_segment(point,first,last)=all(axis->min(first[axis],last[axis])<=point[axis]<=
        max(first[axis],last[axis]),1:2)
    for first_column in 1:length(triangles)-1,last_column in first_column+1:length(triangles)
        first_triangle=triangles[first_column].nodes;last_triangle=triangles[last_column].nodes
        for i in 1:3,j in 1:3
            first_edge=(first_triangle[i],first_triangle[mod1(i+1,3)])
            last_edge=(last_triangle[j],last_triangle[mod1(j+1,3)])
            minmax(first_edge...)==minmax(last_edge...) && continue
            a,b=map(project,first_edge);c,d=map(project,last_edge)
            ab_c=orient2(a,b,c);ab_d=orient2(a,b,d)
            cd_a=orient2(c,d,a);cd_b=orient2(c,d,b)
            ab_c!=0 && ab_d!=0 && sign(ab_c)!=sign(ab_d) &&
                cd_a!=0 && cd_b!=0 && sign(cd_a)!=sign(cd_b) && throw(ArgumentError(
                "$caller: embedded Surface actual triangle edges cross"))
            for (uid,point,side) in ((last_edge[1],c,ab_c),(last_edge[2],d,ab_d))
                side==0 && on_segment(point,a,b) && !(uid in first_edge) && throw(ArgumentError(
                    "$caller: embedded Surface has an overlapping edge or unbound edge vertex"))
            end
            for (uid,point,side) in ((first_edge[1],a,cd_a),(first_edge[2],b,cd_b))
                side==0 && on_segment(point,c,d) && !(uid in last_edge) && throw(ArgumentError(
                    "$caller: embedded Surface has an overlapping edge or unbound edge vertex"))
            end
        end
        for (selected,other) in ((first_triangle,last_triangle),(last_triangle,first_triangle)),uid in selected
            uid in other && continue
            point=project(uid)
            sides=ntuple(i->orient2(project(other[i]),project(other[mod1(i+1,3)]),point),3)
            all(!iszero,sides) && length(Set(sign.(sides)))==1 && throw(ArgumentError(
                "$caller: embedded Surface actual triangle interiors overlap"))
        end
    end
    return nothing
end

function _model_native_hybrid_surface_source_certificate(witness,input,uids,coordinates,caller)
    witness===nothing && return nothing
    witness isa _ModelNativeHybridEmbeddedSurfaceInput && witness.input===input &&
        get(witness.model.embeds,(3,witness.region),NTuple{2,Int}[])==witness.inventory &&
        (2,witness.surface) in witness.inventory || throw(ArgumentError(
        "$caller: embedded Surface has a foreign or changed constructor"))
    fresh=_model_native_hybrid_embedded_surface_input(witness.model,witness.region,input,caller;
        surface=witness.surface)
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
    entries=_model_native_hybrid_convex_shell(panels,coordinate,boundary,caller,"embedded Surface")
    if entries!==nothing
        _model_native_hybrid_convex_inside(witness.records,entries,caller,"embedded Surface")
    else
        for record in witness.records
            _model_native_hybrid_shell_inside(record.coordinate,panels,coordinate,
                                              caller,"embedded Surface")
        end
        lookup=Dict(record.uid=>record.coordinate for record in witness.records)
        # Every literal sheet edge stays clear of the shell panels, and no
        # literal panel edge pierces a sheet facet: strictly-inside nodes
        # plus pairwise non-intersection certify the lifted disk inside.
        sheet_edges=Set(minmax(triangle.nodes[index],triangle.nodes[mod1(index+1,3)])
                        for triangle in witness.triangles for index in 1:3)
        _model_native_hybrid_shell_segments_inside(
            sheet_edges,uid->lookup[uid],panels,coordinate,caller,"embedded Surface")
        for panel in panels,panel_triangle in _model_native_hybrid_panel_triangles(panel)
            for index in 1:3
                p,q=map(coordinate,(panel_triangle[index],panel_triangle[mod1(index+1,3)]))
                for triangle in witness.triangles
                    a,b,c=map(uid->lookup[uid],triangle.nodes)
                    _model_native_hybrid_segment_triangle_relation(p,q,a,b,c)==0 ||
                        throw(ArgumentError(
                        "$caller: embedded Surface is outside or on the hybrid shell boundary"))
                end
            end
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
