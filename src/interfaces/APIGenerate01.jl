# Detached lower-dimensional generation. Raw vertex identity is its source tag;
# coordinates are never used to merge distinct source vertices.
struct _Dim01Node
    tag::UInt64
    point::NTuple{3,Float64}
    owner::Tuple{Int,Int32}
    parameter::Union{Nothing,Vector{Float64}}
end

mutable struct _Dim01Cell
    msh::Int
    owner::Tuple{Int,Int32}
    nodes::Vector{Int32}
    tag::UInt64
    source_tag::UInt64
    physical::Int32
end

mutable struct _Dim01Assembly
    nodes::Vector{_Dim01Node}
    cells::Vector{_Dim01Cell}
    node_positions::Dict{UInt64,Int32}
    node_max::UInt64
    cell_max::UInt64
    max_nodes::Int
    max_cells::Int
    caller::String
end

function _dim01_next_tag!(assembly::_Dim01Assembly,node::Bool)
    previous=node ? assembly.node_max : assembly.cell_max
    previous<typemax(UInt64) || throw(ArgumentError(
        "$(assembly.caller): mesh tag allocation exceeds UInt64"))
    value=previous+UInt64(1)
    node ? (assembly.node_max=value) : (assembly.cell_max=value)
    return value
end

function _dim01_add_node!(assembly::_Dim01Assembly,tag::UInt64,point,owner,
                          parameter=nothing)
    position=get(assembly.node_positions,tag,Int32(0))
    if position!=0
        previous=assembly.nodes[position]
        previous.owner==owner && previous.point==point || throw(ArgumentError(
            "$(assembly.caller): source node $tag has conflicting ownership or coordinates"))
        # Raw records can omit parameters that are retained in the classified
        # cache. Fill that absence without replacing caller-provided values.
        if previous.parameter===nothing && parameter!==nothing
            assembly.nodes[position]=_Dim01Node(tag,point,owner,copy(parameter))
        end
        return position
    end
    length(assembly.nodes)<assembly.max_nodes || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_nodes=$(assembly.max_nodes)"))
    all(isfinite,point) || throw(ArgumentError(
        "$(assembly.caller): source node $tag has non-finite coordinates"))
    push!(assembly.nodes,_Dim01Node(tag,point,owner,parameter))
    position=Int32(length(assembly.nodes))
    assembly.node_positions[tag]=position
    assembly.node_max=max(assembly.node_max,tag)
    return position
end

function _dim01_add_cell!(assembly::_Dim01Assembly,msh,owner,nodes,tag;
                          source_tag=tag,physical=Int32(0))
    length(assembly.cells)<assembly.max_cells || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_cells=$(assembly.max_cells)"))
    spec=msh_spec(msh)
    length(nodes)==spec.nnodes && msh_dimension(msh)==owner[1] ||
        throw(ArgumentError("$(assembly.caller): invalid source MSH $msh connectivity"))
    spec.family===:pnt || spec.order in (1,2) || throw(ArgumentError(
        "$(assembly.caller): lower-dimensional generation supports orders one and two; " *
        "source MSH $msh has order $(spec.order)"))
    push!(assembly.cells,_Dim01Cell(Int(msh),owner,Vector{Int32}(nodes),
                                   UInt64(tag),UInt64(source_tag),physical))
    assembly.cell_max=max(assembly.cell_max,UInt64(tag))
    return nothing
end

_dim01_record(m,entity)=get(m.discrete,entity,get(m.meshing.attached,entity,nothing))

function _dim01_cell_order(cell)
    family=msh_spec(cell.msh).family
    rank=family===:pnt ? 0 : family===:lin ? 1 : family===:tri ? 2 :
         family===:qua ? 3 : family===:tet ? 4 : family===:hex ? 5 :
         family===:pri ? 6 : family===:pyr ? 7 : 8
    return (cell.owner,rank)
end

function _dim01_fully_discrete(m,entity)
    record=get(m.discrete,entity,nothing)
    record===nothing && return false
    entity[1]==1 && return !Model._model_discrete_curve_parametrized(record)
    return true
end

function _dim01_physical(m,entity)
    return Model._model_projection_legacy_tag(
        Model._model_projection_physical_tags(m,entity[1],entity[2]))
end

function _dim01_ingest_record_cells!(assembly,m,records)
    cell_tags=Dict{UInt64,Int}()
    raw_cells=Set{Tuple{Int,Tuple{Int,Int32},Tuple}}()
    for ((dim,tag),record) in records
        owner=(dim,Int32(tag))
        for column in eachindex(record.element_tags)
            nodes=Int32[]
            for source in record.element_nodes[column]
                position=get(assembly.node_positions,UInt64(source),Int32(0))
                position!=0 || throw(ArgumentError(
                    "$(assembly.caller): element $(record.element_tags[column]) " *
                    "references unknown source node $source"))
                push!(nodes,position)
            end
            external=UInt64(record.element_tags[column])
            haskey(cell_tags,external) && throw(ArgumentError(
                "$(assembly.caller): duplicate source element tag $external"))
            _dim01_add_cell!(assembly,record.element_types[column],owner,nodes,external;
                             physical=_dim01_physical(m,(dim,tag)))
            push!(raw_cells,(Int(record.element_types[column]),owner,Tuple(nodes)))
            cell_tags[external]=length(assembly.cells)
        end
    end
    return cell_tags,raw_cells
end

function _dim01_ingest!(assembly,m,cached,cached_class)
    records=sort!(collect(Model._discrete_mesh_records_model(m));by=first)
    for ((dim,tag),record) in records
        owner=(dim,Int32(tag))
        for column in eachindex(record.node_tags)
            parameter=size(record.node_params,1)==dim &&
                      size(record.node_params,2)==length(record.node_tags) ?
                      collect(@view record.node_params[:,column]) : nothing
            point=ntuple(axis->record.node_coords[axis,column],3)
            _dim01_add_node!(assembly,UInt64(record.node_tags[column]),point,owner,parameter)
        end
    end
    if cached===nothing
        _dim01_ingest_record_cells!(assembly,m,records)
        return Int32[]
    end
    cached_class===nothing && throw(ArgumentError(
        "$(assembly.caller): existing mesh needs entity classification"))
    data=cached isa MixedMesh ? cached.entity_data : nothing
    authoritative=cached_class.public_tags!==nothing
    mapping=Vector{Int32}(undef,nnodes(cached))
    # A generated cache without an external-tag table predates this planner.
    # Match only a unique vertex of the same owner; ambiguous raw identity is
    # rejected instead of guessing which coincident source tag it represented.
    raw_positions=Dict{Tuple{Tuple{Int,Int32},NTuple{3,Float64}},Vector{Int32}}()
    for (position,node) in enumerate(assembly.nodes)
        push!(get!(()->Int32[],raw_positions,(node.owner,node.point)),Int32(position))
    end
    for column in axes(cached.coords,2)
        point=ntuple(axis->cached.coords[axis,column],3)
        owner=cached_class.node_entities[column]
        parameter=data===nothing ? nothing : data.node_parametric[column]
        if authoritative
            mapping[column]=_dim01_add_node!(assembly,cached_class.public_tags.node_tags[column],
                point,owner,parameter)
        else
            matches=get(raw_positions,(owner,point),Int32[])
            length(matches)<=1 || throw(ArgumentError(
                "$(assembly.caller): an older cache has ambiguous coincident " *
                "raw source vertices on entity $owner"))
            if length(matches)==1
                mapping[column]=matches[1]
                node=assembly.nodes[mapping[column]]
                if node.parameter===nothing && parameter!==nothing
                    assembly.nodes[mapping[column]]=_Dim01Node(
                        node.tag,node.point,node.owner,parameter)
                end
            else
                external=UInt64(column)
                haskey(assembly.node_positions,external) &&
                    (external=_dim01_next_tag!(assembly,true))
                mapping[column]=_dim01_add_node!(assembly,external,point,owner,parameter)
            end
        end
    end
    # A raw cell may reference a published cache-owned primary/support node.
    # Resolve those exact labels before parsing record connectivity; raw nodes
    # still establish ownership first and parameter enrichment remains intact.
    cell_tags,raw_cells=_dim01_ingest_record_cells!(assembly,m,records)
    for (bi,(msh,dim,offset,cells,owners)) in enumerate(_cache_catalog(cached,cached_class))
        physical=_cache_native_blocks(cached)[bi][3]
        for column in axes(cells,2)
            nodes=mapping[cells[:,column]]
            owner=(dim,owners[column])
            # Legacy native projection labels are synthetic. A lower boundary
            # cell can already exist in an attached sparse source record.
            # Preserve that actual source tag rather than mirror it twice.
            !authoritative && (Int(msh),owner,Tuple(nodes)) in raw_cells && continue
            external=authoritative ? cached_class.public_tags.element_tags[offset+column] :
                     data===nothing ? UInt64(offset+column) :
                                      data.external_element_tags[bi][column]
            previous=get(cell_tags,external,0)
            if previous!=0
                old=assembly.cells[previous]
                if old.owner==owner && old.msh==msh && old.nodes==nodes
                    continue
                end
                !authoritative || throw(ArgumentError(
                    "$(assembly.caller): cached element $external conflicts with its source record"))
                external=_dim01_next_tag!(assembly,false)
            end
            _dim01_add_cell!(assembly,msh,owner,nodes,external;physical=physical[column])
            cell_tags[external]=length(assembly.cells)
        end
    end
    return mapping
end

function _dim01_first_order!(assembly,m;skip_discrete=true)
    high_order=falses(length(assembly.nodes))
    for cell in assembly.cells
        cell.owner[1]==0 && continue
        skip_discrete && _dim01_fully_discrete(m,(cell.owner[1],Int(cell.owner[2]))) && continue
        spec=msh_spec(cell.msh)
        target=msh_type(spec.family,1)
        primary=Elements.msh_num_nodes(target)
        for node in @view cell.nodes[primary+1:end]
            high_order[node]=true
        end
        resize!(cell.nodes,primary)
        cell.msh=target
        # Gmsh's SetOrder1 reconstructs even an already linear native cell.
        cell.tag=_dim01_next_tag!(assembly,false)
    end
    return high_order
end

function _dim01_point_phase!(assembly,m,entities)
    first_nodes=Dict{Int,Int32}()
    owned=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (position,node) in enumerate(assembly.nodes)
        push!(get!(()->Int32[],owned,node.owner),Int32(position))
    end
    have_cells=Set(cell.owner for cell in assembly.cells)
    for (dim,tag) in entities
        dim==0 || continue
        owner=(0,Int32(tag))
        nodes=get(owned,owner,Int32[])
        if isempty(nodes)
            point=get(m.points,tag,nothing)
            point===nothing && continue
            node=_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),point,owner,Float64[])
            push!(nodes,node)
        end
        first_nodes[tag]=nodes[1]
        owner in have_cells && continue
        _dim01_add_cell!(assembly,15,owner,Int32[nodes[end]],
            _dim01_next_tag!(assembly,false);source_tag=UInt64(0),
            physical=_dim01_physical(m,(dim,tag)))
    end
    return first_nodes
end

function _dim01_emit_curve!(assembly,m,curve,params,first_nodes;
                            reuse=false,high_order=BitVector(),owned_nodes=nothing)
    owner=(1,Int32(curve))
    endpoints=get(m.curves,curve,nothing)
    if endpoints===nothing
        record=_dim01_record(m,(1,curve))
        points=record===nothing ? Int[] :
            Int[tag for (dim,tag) in record.boundary if dim==0]
        length(points) in (1,2) || throw(ArgumentError(
            "$(assembly.caller): parametrized discrete Curve[$curve] " *
            "needs one or two boundary points"))
        endpoints=(points[1],points[end])
    end
    begin_point,end_point=endpoints
    begin_node=get(first_nodes,begin_point,Int32(0))
    end_node=get(first_nodes,end_point,Int32(0))
    begin_node!=0 && end_node!=0 || throw(ArgumentError(
        "$(assembly.caller): Curve[$curve] has no boundary mesh vertices"))
    length(params)>=2 || return nothing
    chain=Int32[begin_node]
    sample_positions=Int32[]
    sample_values=Float64[]
    sample_axis=1
    if reuse
        positions=owned_nodes===nothing ? eachindex(assembly.nodes) :
                  get(owned_nodes,owner,Int32[])
        for position in positions
            source=assembly.nodes[position]
            source.owner==owner || continue
            position<=length(high_order) && high_order[position] && continue
            push!(sample_positions,Int32(position))
        end
        if !isempty(sample_positions)
            spans=ntuple(3) do axis
                values=(assembly.nodes[position].point[axis] for position in sample_positions)
                maximum(values)-minimum(values)
            end
            sample_axis=argmax(spans)
            sample_values=let axis=sample_axis
                sort!(sample_positions;by=position->assembly.nodes[position].point[axis])
                Float64[assembly.nodes[position].point[axis]
                        for position in sample_positions]
            end
        end
    end
    for slot in 2:length(params)-1
        parameter=params[slot]
        point=Model._model_curve_part_point(m,curve,parameter,assembly.caller)
        node=Int32(0)
        if reuse
            # Older surface/volume caches do not store every boundary line,
            # but the model retains its 1D parameters and classified nodes.
            # Match a unique actual sample; never merge ambiguous raw tags.
            candidates=Int32[]
            tolerance=64eps(Float64)*max(1.0,maximum(abs,point))
            lo=searchsortedfirst(sample_values,prevfloat(point[sample_axis]-tolerance))
            hi=searchsortedlast(sample_values,nextfloat(point[sample_axis]+tolerance))
            for sample in lo:hi
                position=sample_positions[sample]
                source=assembly.nodes[position]
                maximum(abs(source.point[axis]-point[axis]) for axis in 1:3)<=tolerance || continue
                push!(candidates,Int32(position))
            end
            length(candidates)<=1 || throw(ArgumentError(
                "$(assembly.caller): retained Curve[$curve] sample has ambiguous source identity"))
            isempty(candidates) || (node=candidates[1])
        end
        if node==0
            node=_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),
                                 point,owner,Float64[parameter])
        end
        push!(chain,node)
    end
    push!(chain,end_node)
    for slot in 1:length(chain)-1
        _dim01_add_cell!(assembly,1,owner,chain[slot:slot+1],
            _dim01_next_tag!(assembly,false);source_tag=UInt64(0),
            physical=_dim01_physical(m,(1,curve)))
    end
    return nothing
end

function _dim01_curve_phase!(assembly,m,options,entities,first_nodes,high_order)
    pending=Set{Int}(tag for (dim,tag) in entities if dim==1)
    removed=Set{Tuple{Int,Int32}}()
    have_cells=Set(cell.owner for cell in assembly.cells)
    source_cells=length(assembly.cells)
    owned_nodes=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (position,node) in enumerate(assembly.nodes)
        push!(get!(()->Int32[],owned_nodes,node.owner),Int32(position))
    end
    iterations=0
    while !isempty(pending)
        for curve in sort!(collect(pending))
            owner=(1,Int32(curve))
            # meshGEdge's skip gates preserve the stored lines. The grading
            # helper's raw-record lookup covers fully discrete entities, while
            # this gate also handles native attached records and cached lines.
            keep=_dim01_fully_discrete(m,(1,curve)) ||
                 Model._model_curve_degenerate_kind(m,curve)==1 ||
                 (options.mesh_only_visible && model_entity_visibility(m,1,curve)==0) ||
                 (options.mesh_only_empty && owner in have_cells)
            if keep
                if !(owner in have_cells) &&
                   haskey(m.curve_params,curve)
                    _dim01_emit_curve!(assembly,m,curve,m.curve_params[curve],first_nodes;
                                       reuse=true,high_order=high_order,owned_nodes=owned_nodes)
                    push!(have_cells,owner)
                end
                delete!(pending,curve)
                continue
            end
            result=Model._model_mesh_curve!(m,curve,options,assembly.caller)
            result===:pending && continue
            if result===:keep
                if !(owner in have_cells) &&
                   haskey(m.curve_params,curve)
                    _dim01_emit_curve!(assembly,m,curve,m.curve_params[curve],first_nodes;
                                       reuse=true,high_order=high_order,owned_nodes=owned_nodes)
                    push!(have_cells,owner)
                end
                delete!(pending,curve)
                continue
            end
            params=result
            m.curve_params[curve]=params
            push!(removed,owner)
            _dim01_emit_curve!(assembly,m,curve,params,first_nodes)
            push!(have_cells,owner)
            delete!(pending,curve)
        end
        isempty(pending) && break
        iterations+=1
        iterations>options.max_retries && break
    end
    # Keep generated cells in their allocation order and remove stale source
    # cells in one pass, instead of scanning all cells for every meshed curve.
    if !isempty(removed)
        write=1
        for read in eachindex(assembly.cells)
            read<=source_cells && assembly.cells[read].owner in removed && continue
            assembly.cells[write]=assembly.cells[read]
            write+=1
        end
        resize!(assembly.cells,write-1)
    end
    return removed
end

function _dim01_support_key(primary,weights)
    terms=Tuple{Int32,Float64}[(primary[i],Float64(weights[i]))
                              for i in eachindex(primary) if !iszero(weights[i])]
    sort!(terms;by=first)
    return Tuple(terms)
end

function _dim01_entity_closures(m,entities,assembly)
    owned=Dict(entity=>Set{Int32}() for entity in entities)
    for (position,node) in enumerate(assembly.nodes)
        push!(get!(()->Set{Int32}(),owned,(node.owner[1],Int(node.owner[2]))),Int32(position))
    end
    for cell in assembly.cells
        union!(get!(()->Set{Int32}(),owned,(cell.owner[1],Int(cell.owner[2]))),cell.nodes)
    end
    # Include the actual boundary closure, even when no explicit boundary
    # elements were present in the source cache.
    for entity in entities
        for (_,boundary) in model_boundary(m,[entity],false,false,true)
            union!(owned[entity],get(owned,(0,boundary),Set{Int32}()))
        end
        dim,tag=entity
        dim==0 && continue
        for child in model_boundary(m,[entity],false,false,false)
            union!(owned[entity],get(owned,child,Set{Int32}()))
        end
    end
    return owned
end

function _dim01_quadratic!(assembly,m,entities;only_visible=false)
    entity_orders=Dict{Tuple{Int,Int32},Int}()
    for cell in assembly.cells
        cell.owner[1]>0 || continue
        only_visible && model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))==0 && continue
        order=msh_spec(cell.msh).order
        previous=get(entity_orders,cell.owner,order)
        previous==order || throw(ArgumentError(
            "$(assembly.caller): mixed polynomial orders on entity $(cell.owner); " *
            "Gmsh's first-element high-order transition is not implemented"))
        entity_orders[cell.owner]=order
    end
    any(cell->cell.owner[1]>0 && msh_spec(cell.msh).order==1 &&
        (!only_visible || model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))!=0),
        assembly.cells) || return false
    # The general affine full-P2 stencils share topological supports. Native
    # non-affine CAD placement remains explicit until its geometry stencil is
    # available; already quadratic discrete cells retain their actual supports.
    native_entities=Set{Tuple{Int,Int32}}(cell.owner for cell in assembly.cells
        if cell.owner[1]>0 && msh_spec(cell.msh).order==1 &&
           (!only_visible || model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))!=0) &&
           !_dim01_fully_discrete(m,(cell.owner[1],Int(cell.owner[2]))))
    if !isempty(native_entities)
        boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}()
        for (dim,tag) in entities
            key=(dim,Int32(tag))
            key in native_entities || continue
            boundaries[key]=Int32[abs(t) for (_,t) in
                model_boundary(m,[(dim,tag)],false,false,false)]
            for child in model_boundary(m,[(dim,tag)],false,false,false)
                boundaries[(child[1],Int32(child[2]))]=Int32[]
            end
        end
        guard_class=(boundaries=boundaries,node_entities=collect(native_entities))
        _mixed_require_linear_cad(m,guard_class,assembly.caller)
    end
    closures=_dim01_entity_closures(m,entities,assembly)
    incident_entities=Dict{Int32,Vector{Tuple{Int,Int}}}()
    for entity in entities,node in closures[entity]
        push!(get!(()->Tuple{Int,Int}[],incident_entities,node),entity)
    end
    support_nodes=Dict{Tuple,Int32}()
    # Existing discrete P2 nodes must be used as supports by adjacent linear
    # cells, rather than being recreated at straight interpolation positions.
    for cell in assembly.cells
        spec=msh_spec(cell.msh)
        cell.owner[1]>0 && spec.order==2 || continue
        only_visible && model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))==0 && continue
        primary=Elements.msh_num_nodes(msh_type(spec.family,1))
        reference=Elements.lagrange_nodes(cell.msh)
        for slot in primary+1:spec.nnodes
            weights=MeshFunctionSpaces._first_order_values(Val(spec.family),
                reference[1,slot],reference[2,slot],reference[3,slot],assembly.caller,slot)
            key=_dim01_support_key(@view(cell.nodes[1:primary]),weights)
            previous=get(support_nodes,key,cell.nodes[slot])
            previous==cell.nodes[slot] || throw(ArgumentError(
                "$(assembly.caller): existing quadratic cells have conflicting support-node identities"))
            support_nodes[key]=cell.nodes[slot]
        end
    end
    for cell in assembly.cells
        spec=msh_spec(cell.msh)
        cell.owner[1]>0 && spec.order==1 || continue
        only_visible && model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))==0 && continue
        target=msh_type(spec.family,2)
        reference=Elements.lagrange_nodes(target)
        primary=copy(cell.nodes)
        for slot in length(primary)+1:size(reference,2)
            weights=MeshFunctionSpaces._first_order_values(Val(spec.family),
                reference[1,slot],reference[2,slot],reference[3,slot],assembly.caller,slot)
            key=_dim01_support_key(primary,weights)
            node=get(support_nodes,key,Int32(0))
            if node==0
                point=ntuple(axis->sum(weight*assembly.nodes[index].point[axis]
                                     for (index,weight) in key),3)
                owner=cell.owner
                for candidate in get(incident_entities,key[1][1],Tuple{Int,Int}[])
                    candidate[1]>cell.owner[1] && break
                    all(term->term[1] in closures[candidate],key) || continue
                    owner=(candidate[1],Int32(candidate[2]))
                    break
                end
                parameter=nothing
                if owner[1]==1 && haskey(m.curves,Int(owner[2]))
                    parameters=Model.model_parametrization(m,1,Int(owner[2]),collect(point))
                    parameter=Float64[parameters[1]]
                    point=Model._model_curve_part_point(m,Int(owner[2]),parameters[1],assembly.caller)
                elseif owner[1]==2 && haskey(m.surfaces,Int(owner[2]))
                    parameter=Model.model_parametrization(m,2,Int(owner[2]),collect(point))
                end
                node=_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),point,owner,parameter)
                support_nodes[key]=node
            end
            push!(cell.nodes,node)
        end
        cell.msh=target
        cell.tag=_dim01_next_tag!(assembly,false)
        _api_p2_constructed_cell_certify(assembly.nodes,cell,assembly.caller)
    end
    return true
end

function _dim01_entities(m,entities,assembly,kept)
    metadata=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}()
    owned=Dict{Tuple{Int,Int32},Vector{Int}}()
    for position in kept
        push!(get!(()->Int[],owned,assembly.nodes[position].owner),position)
    end
    for entity in entities
        dim,tag=entity
        boundary=dim==0 ? Int32[] : Int32[t for (_,t) in
            model_boundary(m,[entity],false,true,false)]
        boundaries[(dim,Int32(tag))]=abs.(boundary)
        record=_dim01_record(m,entity)
        bounds=if haskey(m.discrete,entity)
            positions=get(owned,(dim,Int32(tag)),Int[])
            if isempty(positions)
                (0.0,0.0,0.0,0.0,0.0,0.0)
            elseif dim==0
                point=assembly.nodes[positions[1]].point
                (point...,point...)
            else
                ntuple(6) do axis
                    coordinate=mod1(axis,3)
                    values=(assembly.nodes[i].point[coordinate] for i in positions)
                    axis<=3 ? minimum(values) : maximum(values)
                end
            end
        else
            model_bounding_box(m,dim,tag)
        end
        embedded=dim==2 ? Int32[t for (d,t) in get(m.embeds,entity,NTuple{2,Int}[]) if d==1] : Int32[]
        metadata[entity]=Elements.MixedEntity(dim,tag,bounds;
            physical_tags=Model._model_projection_physical_tags(m,dim,tag),
            boundaries=boundary,embedded_curves=embedded)
    end
    return metadata,boundaries
end

function _dim01_finish!(assembly,m,entities,discarded,high_order,renumber,p2_pass;
                        save_all=false,node_limit=assembly.max_nodes,
                        cell_limit=assembly.max_cells,source_mesh=nothing,
                        source_node_positions=Int32[],curve_parameter_order::Bool=false)
    referenced=falses(length(assembly.nodes))
    for cell in assembly.cells,node in cell.nodes
        referenced[node]=true
    end
    if p2_pass
        # pruneMeshVertexAssociations clears associations and then visits
        # volume, surface, curve and Point cells. The first entity at the
        # lowest incident dimension owns the vertex, including old primaries.
        incident=fill((4,Int32(0)),length(assembly.nodes))
        cell_order=sortperm(assembly.cells;by=cell->(-cell.owner[1],cell.owner[2]),
                            alg=Base.Sort.MergeSort)
        for position in cell_order
            cell=assembly.cells[position]
            for node in cell.nodes
                incident[node][1]>cell.owner[1] && (incident[node]=cell.owner)
            end
        end
        for position in eachindex(assembly.nodes)
            referenced[position] || continue
            node=assembly.nodes[position]
            owner=incident[position]
            parameter=node.owner==owner ? node.parameter : nothing
            assembly.nodes[position]=_Dim01Node(node.tag,node.point,owner,parameter)
        end
    end
    # Initial SetOrder1 removes old native interpolation vertices. Remeshing
    # deletes native entity-owned vertices; retained discrete references survive.
    kept=Int[i for (i,node) in enumerate(assembly.nodes) if referenced[i] ||
        (!p2_pass && !(node.owner in discarded) &&
         !(i<=length(high_order) && high_order[i]))]
    length(kept)<=node_limit || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_nodes=$node_limit"))
    length(assembly.cells)<=cell_limit || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_cells=$cell_limit"))
    if curve_parameter_order
        # meshRefine.cpp: Subdivide(GEdge) restores mesh_vertices to native
        # parameter order after the high-order support pass. Ordinary order
        # changes retain their own source-tag ordering and do not enter here.
        for position in kept
            node=assembly.nodes[position]
            node.owner[1]==1 && haskey(m.curves,Int(node.owner[2])) || continue
            parameter=node.parameter
            if parameter===nothing || length(parameter)!=1
                parameter=Model.model_parametrization(m,1,Int(node.owner[2]),collect(node.point))
                assembly.nodes[position]=_Dim01Node(node.tag,node.point,node.owner,parameter)
            end
        end
        sort!(kept;by=i->begin
            node=assembly.nodes[i]
            native=node.owner[1]==1 && haskey(m.curves,Int(node.owner[2]))
            (node.owner,native ? node.parameter[1] : 0.0,
             p2_pass && !native ? node.tag : UInt64(0))
        end,alg=Base.Sort.MergeSort)
    elseif p2_pass
        sort!(kept;by=i->(assembly.nodes[i].owner,assembly.nodes[i].tag),alg=Base.Sort.MergeSort)
    else
        sort!(kept;by=i->assembly.nodes[i].owner,alg=Base.Sort.MergeSort)
    end
    mapping=zeros(Int32,length(assembly.nodes))
    mapping[kept]=Int32.(1:length(kept))
    save_physical=!save_all && any(!isempty,values(m.physical))
    if renumber && save_physical && !isempty(kept)
        # The partial-save renumber pass first forceNums every vertex to the
        # sentinel nv+1. forceNum advances Gmsh's allocation high-water mark.
        assembly.node_max=max(assembly.node_max,UInt64(length(kept))+UInt64(1))
    end
    physical_support=falses(length(assembly.nodes))
    if save_physical
        for cell in assembly.cells
            _dim01_physical(m,(cell.owner[1],Int(cell.owner[2])))!=0 || continue
            physical_support[cell.nodes].=true
        end
    end
    node_tags=UInt64[assembly.nodes[i].tag for i in kept]
    if renumber
        next=UInt64(0)
        for selected in (true,false), (slot,position) in enumerate(kept)
            save_physical && physical_support[position]!=selected && continue
            !save_physical && !selected && continue
            next+=UInt64(1)
            node_tags[slot]=next
        end
    end
    node_tag_map=Dict(assembly.nodes[i].tag=>node_tags[slot] for (slot,i) in enumerate(kept))
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    # Display flags belong to the surviving MElement, not to a reused numeric
    # tag. Factory order conversions and remeshing allocate a new identity.
    visibility_tag_map=Dict{UInt64,UInt64}()
    if renumber
        next=UInt64(0)
        for selected in (true,false),cell in assembly.cells
            saved=_dim01_physical(m,(cell.owner[1],Int(cell.owner[2])))!=0
            save_physical && saved!=selected && continue
            !save_physical && !selected && continue
            next+=UInt64(1)
            cell.source_tag>0 && cell.tag==cell.source_tag &&
                (visibility_tag_map[cell.source_tag]=next)
            cell.tag=next
        end
    end
    order=Int[];groups=Dict{Int,Vector{Int}}()
    element_tag_map=Dict{UInt64,UInt64}()
    for (slot,cell) in enumerate(assembly.cells)
        cell.source_tag!=0 && (element_tag_map[cell.source_tag]=cell.tag)
        !renumber && cell.source_tag>0 && cell.tag==cell.source_tag &&
            (visibility_tag_map[cell.source_tag]=cell.tag)
        haskey(groups,cell.msh) || push!(order,cell.msh)
        push!(get!(()->Int[],groups,cell.msh),slot)
    end
    coords=Matrix{Float64}(undef,3,length(kept))
    node_entities=Tuple{Int,Int32}[]
    node_parameters=Union{Nothing,Vector{Float64}}[]
    for (slot,i) in enumerate(kept)
        source=assembly.nodes[i]
        for axis in 1:3
            coords[axis,slot]=source.point[axis]
        end
        push!(node_entities,source.owner)
        push!(node_parameters,source.parameter)
    end
    metadata,boundaries=_dim01_entities(m,entities,assembly,kept)
    blocks=ElementBlock[];block_owners=Vector{Int32}[];external_cells=Vector{UInt64}[]
    owners=Dict{Int,Vector{Int32}}()
    for msh in order
        positions=groups[msh]
        nodes=Matrix{Int32}(undef,Elements.msh_num_nodes(msh),length(positions))
        physical=Vector{Int32}(undef,length(positions))
        classified=Vector{Int32}(undef,length(positions))
        external=Vector{UInt64}(undef,length(positions))
        for (column,position) in enumerate(positions)
            cell=assembly.cells[position]
            nodes[:,column]=mapping[cell.nodes]
            physical[column]=cell.physical
            classified[column]=cell.owner[2]
            external[column]=cell.tag
        end
        push!(blocks,ElementBlock(msh,nodes,physical))
        push!(block_owners,classified);push!(external_cells,external)
        owners[msh]=classified
    end
    data=Elements.MixedEntityData(metadata;node_entities=node_entities,
        node_parametric=node_parameters,external_node_tags=node_tags,
        block_entities=block_owners,external_element_tags=external_cells)
    mesh=MixedMesh(coords,blocks;physical_names=m.physical_names,entity_data=data,
                   elementary_entities=block_owners)
    classified_entities=Tuple{Int,Int32}[(dim,Int32(tag)) for (dim,tag) in entities]
    primary=isempty(classified_entities) ? (0,Int32(0)) : classified_entities[end]
    class=_mixed_classification(mesh,primary,classified_entities,node_entities,boundaries,owners)
    links=Elements.MixedPeriodicLink[]
    for ((dim,slave),_) in sort!(collect(m.periodic);by=first)
        dim==1 || continue
        correspondence=_dim01_periodic_nodes(m,mesh,class,slave,assembly.caller;
                                             include_high_order=true)
        push!(links,Elements.MixedPeriodicLink(1,slave,correspondence.master_entity,
            correspondence.slave_nodes,correspondence.master_nodes;
            affine=correspondence.affine))
    end
    if source_mesh isa MixedMesh
        for link in source_mesh.periodic_links
            link.dim==1 && haskey(m.periodic,(1,Int(link.slave_entity))) && continue
            slaves=mapping[source_node_positions[link.slave_nodes]]
            masters=mapping[source_node_positions[link.master_nodes]]
            all(!iszero,slaves) && all(!iszero,masters) || throw(ArgumentError(
                "$(assembly.caller): existing periodic dimension $(link.dim) " *
                "entity $(link.slave_entity) loses interpolation vertices; " *
                "higher-dimensional correspondence reconstruction is not implemented"))
            push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,
                link.master_entity,slaves,masters;affine=link.affine))
        end
    end
    if !isempty(links)
        mesh=MixedMesh(coords,blocks;physical_names=m.physical_names,entity_data=data,
                       elementary_entities=block_owners,periodic_links=links)
        class=_mixed_rebind_class(class,mesh)
    end
    _dim01_reconcile_records!(m,entities,mesh,node_tag_map)
    return (model=m,mesh=mesh,class=class,authority=Set(entities),
        node_tag_map=node_tag_map,element_tag_map=element_tag_map,
        visibility_tag_map=visibility_tag_map,
        max_node_tag=assembly.node_max,max_element_tag=assembly.cell_max)
end

function _dim01_reconcile_records!(m,entities,mesh,node_tag_map)
    data=mesh.entity_data
    maximum(data.external_node_tags;init=UInt64(0))<=typemax(Int32) &&
        maximum((maximum(tags;init=UInt64(0)) for tags in data.external_element_tags);
                init=UInt64(0))<=typemax(Int32) || throw(ArgumentError(
            "API.mesh.generate: reconciled model record tags exceed Int32"))
    positions_by_owner=Dict{Tuple{Int,Int},Vector{Int}}()
    for (position,owner) in enumerate(data.node_entities)
        push!(get!(()->Int[],positions_by_owner,(owner[1],Int(owner[2]))),position)
    end
    records=Dict{Tuple{Int,Int},DiscreteEntity}()
    for entity in entities
        dim,tag=entity
        previous=_dim01_record(m,entity)
        record=previous===nothing ? DiscreteEntity() : deepcopy(previous)
        positions=get(positions_by_owner,entity,Int[])
        record.node_tags=Int32.(data.external_node_tags[positions])
        record.node_coords=mesh.coords[:,positions]
        parameters=data.node_parametric[positions]
        record.node_params=all(value->value!==nothing && length(value)==dim,parameters) ?
            (isempty(parameters) ? zeros(dim,0) : hcat(parameters...)) : zeros(0,length(positions))
        empty!(record.element_types);empty!(record.element_tags);empty!(record.element_nodes)
        record.aux_params=Dict{Int32,Vector{Float64}}(
            Int32(node_tag_map[UInt64(source)])=>values for (source,values) in record.aux_params
            if haskey(node_tag_map,UInt64(source)))
        records[entity]=record
    end
    for (bi,block) in enumerate(mesh.blocks)
        dim=msh_dimension(block.msh)
        for column in axes(block.nodes,2)
            record=get(records,(dim,Int(data.block_entities[bi][column])),nothing)
            record===nothing && continue
            push!(record.element_types,Int32(block.msh))
            push!(record.element_tags,Int32(data.external_element_tags[bi][column]))
            push!(record.element_nodes,Int32.(data.external_node_tags[block.nodes[:,column]]))
        end
    end
    for (entity,record) in records
        if haskey(m.discrete,entity)
            m.discrete[entity]=record
        else
            m.meshing.attached[entity]=record
        end
    end
    return nothing
end

function _generate_dim01_plan(source_model,cached,cached_class,dimension::Int,
                              options::Model._ModelMesh1DOptions,caller;
                              element_order=1,renumber=true,
                              save_all=false,
                              max_nodes=typemax(Int32),max_cells=typemax(Int32),
                              initial_max_node_tag=UInt64(0),
                              initial_max_element_tag=UInt64(0),
                              generation_meshing=nothing)
    dimension in (0,1) || throw(ArgumentError("$caller: dim must be zero or one"))
    element_order in (1,2) || throw(ArgumentError(
        "$caller: mesh order $element_order is not supported (supported: 1, 2)"))
    max_nodes isa Integer && !(max_nodes isa Bool) && 0<=max_nodes<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_nodes must be a nonnegative Int32-representable integer"))
    max_cells isa Integer && !(max_cells isa Bool) && 0<=max_cells<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_cells must be a nonnegative Int32-representable integer"))
    m=deepcopy(source_model)
    # API options describe this generation, while source-cache completion above
    # it still uses the attributes of the mesh being retained. Apply only after
    # staging so a rejected request leaves the live model/session untouched.
    generation_meshing===nothing ||
        _apply_generation_meshing_options!(m,generation_meshing)
    entities=model_entities(m)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    source_node_positions=_dim01_ingest!(assembly,m,cached,cached_class)
    # Staging may retain vertices that the detached pass later deletes. Bound
    # its growth by source size plus requested result capacity; certify the
    # exact result counts before allocating final coordinate/connectivity arrays.
    assembly.max_nodes=min(Int(typemax(Int32)),Base.checked_add(length(assembly.nodes),Int(max_nodes)))
    assembly.max_cells=min(Int(typemax(Int32)),Base.checked_add(length(assembly.cells),Int(max_cells)))
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    high_order=_dim01_first_order!(assembly,m)
    discarded=Set{Tuple{Int,Int32}}()
    if dimension==1
        Model._model_require_point_mesh_geometry(m,caller)
        for entity in entities
            entity[1]>=2 && !_dim01_fully_discrete(m,entity) &&
                push!(discarded,(entity[1],Int32(entity[2])))
        end
        filter!(cell->!(cell.owner in discarded),assembly.cells)
        first_nodes=_dim01_point_phase!(assembly,m,entities)
        union!(discarded,_dim01_curve_phase!(assembly,m,options,entities,first_nodes,high_order))
    end
    # orientMeshGEdge runs on every generation, including dimension zero and
    # retained discrete lines. It reverses the current mesh when requested.
    for cell in assembly.cells
        cell.owner[1]==1 && get(m.meshing.reverse,(1,Int(cell.owner[2])),false) || continue
        cell.nodes=cell.nodes[_mixed_reverse_permutation(cell.msh)]
    end
    # The post-generation high-order pass does not run for a point-only mesh.
    p2_pass=element_order==2 && any(cell->cell.owner[1]>0,assembly.cells)
    p2_pass && _dim01_quadratic!(assembly,m,entities;only_visible=options.mesh_only_visible)
    return _dim01_finish!(assembly,m,entities,discarded,high_order,renumber,p2_pass;
                         save_all=save_all,node_limit=Int(max_nodes),cell_limit=Int(max_cells),
                         source_mesh=cached,source_node_positions=source_node_positions)
end

# Immediate setOrder is distinct from the global generation option: it does
# not run Mesh0D/Mesh1D or the line orientation constraints.
function _set_order_dim01_plan(source_model,cached,cached_class,order::Int,caller;
                               renumber=true,save_all=false,only_visible=false,
                               max_nodes=typemax(Int32),max_cells=typemax(Int32),
                               initial_max_node_tag=UInt64(0),
                               initial_max_element_tag=UInt64(0))
    order in (1,2) || throw(ArgumentError(
        "$caller: mesh order $order is not supported (supported: 1, 2)"))
    max_nodes isa Integer && !(max_nodes isa Bool) && 0<=max_nodes<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_nodes must be a nonnegative Int32-representable integer"))
    max_cells isa Integer && !(max_cells isa Bool) && 0<=max_cells<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_cells must be a nonnegative Int32-representable integer"))
    m=deepcopy(source_model)
    entities=model_entities(m)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    source_node_positions=_dim01_ingest!(assembly,m,cached,cached_class)
    assembly.max_nodes=min(Int(typemax(Int32)),Base.checked_add(length(assembly.nodes),Int(max_nodes)))
    assembly.max_cells=min(Int(typemax(Int32)),Base.checked_add(length(assembly.cells),Int(max_cells)))
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    high_order=order==1 ? _dim01_first_order!(assembly,m;skip_discrete=false) :
                         falses(length(assembly.nodes))
    order==2 && _dim01_quadratic!(assembly,m,entities;only_visible=only_visible)
    return _dim01_finish!(assembly,m,entities,Set{Tuple{Int,Int32}}(),high_order,
                         renumber,order==2;save_all=save_all,
                         node_limit=Int(max_nodes),cell_limit=Int(max_cells),
                         source_mesh=cached,source_node_positions=source_node_positions)
end

# The surface projector's curve-edge catalog deliberately excludes line-only
# caches. This mapper consumes actual classified Line connectivity instead.
# It returns dense indices; the public API translates them through its tag table.
function _dim01_periodic_nodes(m,mesh,class,slave::Int,caller;
                               include_high_order=false)
    constraint=get(m.periodic,(1,slave),nothing)
    constraint===nothing && return (master_entity=slave,slave_nodes=Int32[],
                                    master_nodes=Int32[],affine=nothing)
    master=abs(Int(constraint.master_entity))
    master_set=Set{Int32}();slave_set=Set{Int32}()
    for (msh,dim,_,cells,owners) in _cache_catalog(mesh,class)
        dim==1 || continue
        width=include_high_order ? size(cells,1) :
              Elements.msh_num_nodes(msh_type(msh_spec(msh).family,1))
        for column in axes(cells,2)
            target=owners[column]==master ? master_set :
                   owners[column]==slave ? slave_set : nothing
            target===nothing && continue
            union!(target,@view cells[1:width,column])
        end
    end
    (isempty(master_set) || isempty(slave_set)) && return (
        master_entity=master,slave_nodes=Int32[],master_nodes=Int32[],
        affine=constraint.affine)
    length(master_set)==length(slave_set) || throw(ArgumentError(
        "$caller: periodic Curve[$slave]/Curve[$master] node counts differ"))
    master_nodes=sort!(collect(master_set))
    slave_nodes=Int32[]
    if constraint.affine===nothing
        # Orientation-only relations share native parameter distributions.
        master_bounds=Model._model_curve_param_bounds(m,master,caller)
        slave_bounds=Model._model_curve_param_bounds(m,slave,caller)
        masters=sort!([(Model.model_parametrization(m,1,master,
                       collect(@view mesh.coords[:,node]))[1],node) for node in master_set])
        slaves=sort!([(Model.model_parametrization(m,1,slave,
                      collect(@view mesh.coords[:,node]))[1],node) for node in slave_set])
        constraint.reversed && reverse!(slaves)
        empty!(master_nodes)
        for ((mu,mnode),(su,snode)) in zip(masters,slaves)
            mapped=constraint.reversed ? slave_bounds[2]-su+master_bounds[1] :
                                         su-slave_bounds[1]+master_bounds[1]
            abs(mu-mapped)<=max(constraint.atol,64eps(Float64)) ||
                throw(ArgumentError("$caller: periodic curve parameters do not match"))
            push!(master_nodes,mnode);push!(slave_nodes,snode)
        end
    else
        positions=Dict{NTuple{3,Float64},Vector{Int32}}()
        for node in slave_set
            point=ntuple(axis->mesh.coords[axis,node]+0.0,3)
            push!(get!(()->Int32[],positions,point),node)
        end
        used=Set{Int32}()
        coefficients,translation=_periodic_affine_3x4(constraint.affine)
        for node in master_nodes
            point=ntuple(axis->mesh.coords[axis,node],3)
            expected=Model._model_affine_point(coefficients,translation,point,caller,Int(node))
            key=ntuple(axis->expected[axis]+0.0,3)
            candidates=Int32[candidate for candidate in get(positions,key,Int32[])
                              if !(candidate in used)]
            if isempty(candidates)
                # Geometry endpoints may satisfy the declared affine relation
                # within its tolerance rather than by identical Float64 bits.
                for candidate in slave_set
                    candidate in used && continue
                    actual=ntuple(axis->mesh.coords[axis,candidate],3)
                    Model._model_point_distance(actual,expected)<=constraint.atol || continue
                    push!(candidates,candidate)
                end
            end
            length(candidates)==1 || throw(ArgumentError(
                "$caller: periodic Curve[$slave] has missing or ambiguous node correspondence"))
            push!(slave_nodes,candidates[1]);push!(used,candidates[1])
        end
    end
    return (master_entity=master,slave_nodes=slave_nodes,
            master_nodes=master_nodes,affine=constraint.affine)
end
