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
    source_tag::Union{Nothing,UInt64}
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
    node_counter_wrapped::Bool
    cell_counter_wrapped::Bool
end

# Detached planning records its own constructor events; it does not mutate session clocks.
_Dim01Assembly(nodes,cells,positions,node_max,cell_max,max_nodes,max_cells,caller)=
    _Dim01Assembly(nodes,cells,positions,node_max,cell_max,max_nodes,max_cells,caller,false,false)

function _dim01_next_tag!(assembly::_Dim01Assembly,node::Bool)
    previous=node ? assembly.node_max : assembly.cell_max
    value=previous+UInt64(1)
    if node
        assembly.node_max=value
        iszero(value) && (assembly.node_counter_wrapped=true)
    else
        assembly.cell_max=value
        iszero(value) && (assembly.cell_counter_wrapped=true)
    end
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
    assembly.node_counter_wrapped || (assembly.node_max=max(assembly.node_max,tag))
    return position
end

function _dim01_add_cell!(assembly::_Dim01Assembly,msh,owner,nodes,tag;
                          source_tag=tag,physical=Int32(0),
                          preserve_source_order::Bool=false)
    length(assembly.cells)<assembly.max_cells || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_cells=$(assembly.max_cells)"))
    spec=msh_spec(msh)
    length(nodes)==spec.nnodes && msh_dimension(msh)==owner[1] ||
        throw(ArgumentError("$(assembly.caller): invalid source MSH $msh connectivity"))
    preserve_source_order || spec.family===:pnt || spec.order in (1,2) ||
        (spec.family===:lin && 3<=spec.order<=10) || throw(ArgumentError(
        "$(assembly.caller): lower-dimensional generation supports orders one and two; " *
        "source MSH $msh has order $(spec.order)"))
    push!(assembly.cells,_Dim01Cell(Int(msh),owner,Vector{Int32}(nodes),
                                   UInt64(tag),source_tag===nothing ? nothing : UInt64(source_tag),physical))
    assembly.cell_counter_wrapped || (assembly.cell_max=max(assembly.cell_max,UInt64(tag)))
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

function _dim01_ingest_record_cells!(assembly,m,records;preserve_source_order::Bool=false,
                                    retired_nodes=nothing)
    cell_tags=Dict{UInt64,Int}()
    raw_cells=Set{Tuple{Int,Tuple{Int,Int32},Tuple}}()
    count=sum(length(record.element_tags) for (_,record) in records;init=0)
    sizehint!(cell_tags,count);sizehint!(raw_cells,count)
    for ((dim,tag),record) in records
        owner=(dim,Int32(tag))
        retired=retired_nodes===nothing ? nothing : get(retired_nodes,owner,nothing)
        for column in eachindex(record.element_tags)
            nodes=Int32[]
            for source in record.element_nodes[column]
                position=get(assembly.node_positions,UInt64(source),Int32(0))
                if retired!==nothing && UInt64(source) in retired
                    # A generation-only placeholder keeps this MElement's
                    # SetOrder1 allocation without reading the retired Point.
                    position=Int32(0)
                end
                position!=0 || (retired!==nothing && UInt64(source) in retired) ||
                    throw(ArgumentError(
                        "$(assembly.caller): element $(record.element_tags[column]) " *
                        "references unknown source node $source"))
                push!(nodes,position)
            end
            external=UInt64(record.element_tags[column])
            haskey(cell_tags,external) && throw(ArgumentError(
                "$(assembly.caller): duplicate source element tag $external"))
            _dim01_add_cell!(assembly,record.element_types[column],owner,nodes,external;
                             physical=_dim01_physical(m,(dim,tag)),preserve_source_order)
            push!(raw_cells,(Int(record.element_types[column]),owner,Tuple(nodes)))
            cell_tags[external]=length(assembly.cells)
        end
    end
    return cell_tags,raw_cells
end

function _dim01_ingest!(assembly,m,cached,cached_class;preserve_source_order::Bool=false,
                       retired_nodes=nothing)
    records=sort!(collect(Model._discrete_mesh_records_model(m));by=first)
    for ((dim,tag),record) in records
        owner=(dim,Int32(tag))
        for column in eachindex(record.node_tags)
            parameter=Model._model_mesh_node_parameters(record,column)
            point=ntuple(axis->record.node_coords[axis,column],3)
            _dim01_add_node!(assembly,UInt64(record.node_tags[column]),point,owner,parameter)
        end
    end
    if cached===nothing
        _dim01_ingest_record_cells!(assembly,m,records;preserve_source_order=preserve_source_order,retired_nodes)
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
    if !authoritative
        sizehint!(raw_positions,length(assembly.nodes))
        for (position,node) in enumerate(assembly.nodes)
            push!(get!(()->Int32[],raw_positions,(node.owner,node.point)),Int32(position))
        end
    end
    for column in axes(cached.coords,2)
        point=ntuple(axis->cached.coords[axis,column],3)
        owner=cached_class.node_entities[column]
        parameter=cached_class.node_parametric[column]
        if authoritative
            tag_uid=cached_class.public_tags.node_tags[column]
            existing=get(assembly.node_positions,tag_uid,Int32(0))
            if existing!=0
                # A classification (ClassifySurfaces) re-owns nodes on the
                # live records without republishing this snapshot, so the
                # record claim is current and the cached owner can lag.
                # Only coordinate divergence is a genuine source conflict.
                known=assembly.nodes[existing]
                known.point==point || throw(ArgumentError(
                    "$(assembly.caller): source node $tag_uid has conflicting " *
                    "ownership or coordinates"))
                if known.parameter===nothing && parameter!==nothing
                    assembly.nodes[existing]=_Dim01Node(
                        known.tag,known.point,known.owner,copy(parameter))
                end
                mapping[column]=existing
            else
                mapping[column]=_dim01_add_node!(assembly,tag_uid,point,owner,parameter)
            end
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
    cell_tags,raw_cells=_dim01_ingest_record_cells!(assembly,m,records;preserve_source_order=preserve_source_order,retired_nodes)
    for (bi,(msh,dim,offset,cells,owners)) in enumerate(_cache_catalog(cached,cached_class))
        physical=_cache_native_blocks(cached)[bi][3]
        for column in axes(cells,2)
            owner=(dim,owners[column])
            nodes=mapping[cells[:,column]]
            retired=retired_nodes===nothing ? nothing : get(retired_nodes,owner,nothing)
            if retired!==nothing
                for slot in eachindex(nodes)
                    assembly.nodes[nodes[slot]].tag in retired && (nodes[slot]=Int32(0))
                end
            end
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
                # Classification re-owns source elements without republishing
                # this snapshot: identical type/connectivity under a new
                # record owner is the same logical element, not a conflict.
                authoritative && old.msh==Int(msh) && old.nodes==nodes && continue
                !authoritative || throw(ArgumentError(
                    "$(assembly.caller): cached element $external conflicts with its source record"))
                external=_dim01_next_tag!(assembly,false)
            end
            _dim01_add_cell!(assembly,msh,owner,nodes,external;
                             physical=physical[column],preserve_source_order)
            cell_tags[external]=length(assembly.cells)
        end
    end
    return mapping
end

function _dim01_first_order!(assembly,m;skip_discrete=true,retired_nodes=nothing)
    high_order=falses(length(assembly.nodes))
    for cell in assembly.cells
        cell.owner[1]==0 && continue
        skip_discrete && _dim01_fully_discrete(m,(cell.owner[1],Int(cell.owner[2]))) && continue
        spec=msh_spec(cell.msh)
        target=msh_type(spec.family,1)
        primary=Elements.msh_num_nodes(target)
        for node in @view cell.nodes[primary+1:end]
            node==0 || (high_order[node]=true)
        end
        resize!(cell.nodes,primary)
        cell.msh=target
        # Gmsh's SetOrder1 reconstructs even an already linear native cell.
        cell.tag=_dim01_next_tag!(assembly,false)
    end
    # Allocation follows the complete native factory traversal above. Only
    # now may the detached transaction remove explicit retired-Point cells.
    retired_nodes===nothing || filter!(cell->!any(iszero,cell.nodes),assembly.cells)
    return high_order
end

function _dim01_point_phase!(assembly,m,entities;create_cells=true,create_nodes=true)
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
            create_nodes || continue
            point=get(m.points,tag,nothing)
            point===nothing && haskey(m.discrete,(0,tag)) &&
                (point=(0.0,0.0,0.0))
            point===nothing && continue
            node=_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),point,owner,Float64[])
            push!(nodes,node)
        end
        first_nodes[tag]=nodes[1]
        create_cells || continue
        owner in have_cells && continue
        _dim01_add_cell!(assembly,15,owner,Int32[nodes[end]],
            _dim01_next_tag!(assembly,false);source_tag=nothing,
            physical=_dim01_physical(m,(dim,tag)))
    end
    return first_nodes
end

function _dim01_curve_endpoint_counts(m)
    counts=Dict{Int,Int}()
    for (curve,(first,last)) in m.curves
        counts[first]=get(counts,first,0)+1
        first==last || (counts[last]=get(counts,last,0)+1)
    end
    for ((dimension,curve),record) in m.discrete
        dimension==1 && !haskey(m.curves,curve) || continue
        points=Set(abs(tag) for (dim,tag) in record.boundary if dim==0)
        for point in points
            counts[point]=get(counts,point,0)+1
        end
    end
    return counts
end

function _dim01_extruded_curve_points(assembly,m,curve,parameters,curve_births)
    curve_births===nothing && return nothing
    attributes=Model._extrude_gate(m,1,curve)
    attributes===nothing && return nothing
    link=get(m.meshing.extrude_sources,(1,curve),nothing)
    (link===nothing || link[1]!=1) && return nothing
    spec=get(m.meshing.extrude_specs,(1,curve),nothing)
    spec===nothing && return nothing
    source_curve=abs(link[2])
    source=get(curve_births,(1,Int32(source_curve)),nothing)
    # A freshly emitted Curve supplies its actual oriented Line-node chain.
    # Retained/imported products without that receipt keep their existing path.
    source===nothing && return nothing
    count=length(parameters)
    length(source)==count || throw(ArgumentError(
        "$(assembly.caller): copied Curve[$curve] source births do not align with its parameters"))
    isempty(attributes.heights) && throw(ArgumentError(
        "$(assembly.caller): copied Curve[$curve] lacks its extrusion height"))
    level=Float64(attributes.heights[end])
    points=Vector{NTuple{3,Float64}}(undef,max(0,count-2))
    for slot in 2:count-1
        source_slot=link[2]>0 ? slot : count-slot+1
        position=get(assembly.node_positions,source[source_slot],Int32(0))
        position>0 || throw(ArgumentError(
            "$(assembly.caller): copied Curve[$curve] loses an actual source vertex"))
        node=assembly.nodes[position]
        node.owner==(1,Int32(source_curve)) || throw(ArgumentError(
            "$(assembly.caller): copied Curve[$curve] source vertex has conflicting ownership"))
        # meshGEdgeExtruded copies XYZ through Extrude and stores the mirrored
        # parameter separately; evaluating the target at that U can change XYZ.
        points[slot-1]=Model._extrude_at(spec,node.point,level,assembly.caller)
    end
    return points
end

function _dim01_emit_curve!(assembly,m,curve,params,first_nodes;
                             reuse=false,high_order=BitVector(),owned_nodes=nothing,
                             geometry_birth=nothing,geometry_parameters=params,
                             curve_births=nothing,endpoint_counts=nothing,birth=nothing,
                             copied_points=nothing)
    owner=(1,Int32(curve))
    endpoints=get(m.curves,curve,nothing)
    if endpoints===nothing
        record=_dim01_record(m,(1,curve))
        points=record===nothing ? Int[] :
            Int[abs(tag) for (dim,tag) in record.boundary if dim==0]
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
    native_geometry=Model._model_native_curve_mesh_geometry(m,curve,assembly.caller)
    if !reuse && begin_point==end_point && native_geometry!==nothing &&
            endpoint_counts!==nothing && get(endpoint_counts,begin_point,0)==1
        source=assembly.nodes[begin_node]
        point=Model._model_high_order_curve_point(
            native_geometry,native_geometry.lower,assembly.caller)
        assembly.nodes[begin_node]=_Dim01Node(source.tag,point,source.owner,source.parameter)
        m.meshing.native_curve_endpoints[begin_point]=
            Model._ModelNativeCurveEndpoint(curve,source.tag,point)
    end
    next_birth=1
    for slot in 2:length(params)-1
        if !reuse && birth!==nothing
            source_row=birth.rows[slot]-1
            while next_birth<source_row
                _dim01_next_tag!(assembly,true)
                next_birth+=1
            end
        end
        parameter=params[slot]
        point=copied_points===nothing ? (native_geometry===nothing ?
            Model._model_curve_part_point(m,curve,parameter,assembly.caller;
                _native_producer=true) :
            Model._model_high_order_curve_point(native_geometry,parameter,assembly.caller)) :
            copied_points[slot-1]
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
                distance=0.0
                for axis in 1:3
                    distance=max(distance,abs(source.point[axis]-point[axis]))
                end
                distance<=tolerance || continue
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
        geometry_birth===nothing || _dim01_geometry_curve_node!(
            geometry_birth,assembly,m,curve,node,geometry_parameters[slot])
        !reuse && birth!==nothing && (next_birth+=1)
        push!(chain,node)
    end
    if !reuse && birth!==nothing
        while next_birth<=birth.count
            _dim01_next_tag!(assembly,true)
            next_birth+=1
        end
    end
    push!(chain,end_node)
    # Store the actual oriented emitter identities before high-order pruning
    # can change endpoint ownership. Point-cell construction is independent.
    curve_births===nothing || (curve_births[owner]=UInt64[assembly.nodes[row].tag for row in chain])
    for slot in 1:length(chain)-1
        _dim01_add_cell!(assembly,1,owner,chain[slot:slot+1],
            _dim01_next_tag!(assembly,false);source_tag=nothing,
            physical=_dim01_physical(m,(1,curve)))
    end
    return nothing
end

function _dim01_curve_phase!(assembly,m,options,entities,first_nodes,high_order;
                              done=nothing,geometry_birth=nothing,curve_births=nothing,
                              reuse_curve_params::Bool=false)
    endpoint_counts=_dim01_curve_endpoint_counts(m)
    pending=Set{Int}(tag for (dim,tag) in entities if dim==1)
    # Generator::Mesh1D resets every Curve to PENDING before its skip gates.
    done===nothing || filter!(owner->owner[1]!=1,done)
    removed=Set{Tuple{Int,Int32}}()
    have_cells=Set(cell.owner for cell in assembly.cells)
    source_cells=length(assembly.cells)
    owned_nodes=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (position,node) in enumerate(assembly.nodes)
        push!(get!(()->Int32[],owned_nodes,node.owner),Int32(position))
    end
    iterations=0
    while !isempty(pending)
        stalled=length(pending)
        ordered=sort!(collect(pending))
        # PAR01 — `General.NumThreads` detaches runs of independent singleton
        # Curves: each eligible produce runs the identical
        # `_model_mesh_curve!` call on a transaction copy under
        # `Threads.@threads`, and the unchanged serial merge below replays or
        # re-runs each slot in tag order.  `nothing` keeps the exact serial
        # path; emission stays serial, so public UIDs are minted in the same
        # order under any worker count.
        prepared=_parallel_curve_prepare(m,ordered,options;
            have_cells,reuse_curve_params)
        parallel_results=prepared===nothing ? nothing :
            fill!(Vector{Any}(undef,length(ordered)),nothing)
        parallel_dirty=Set{Tuple{Int,Int32}}()
        parallel_last=prepared===nothing ? 0 :
            something(findlast(!isnothing,prepared),0)
        parallel_precise=prepared===nothing ? false :
            GmshLibm._win_cli_precision()
        for (index,curve) in enumerate(ordered)
            owner=(1,Int32(curve))
            # meshGEdge's skip gates preserve the stored lines. The grading
            # helper's raw-record lookup covers fully discrete entities, while
            # this gate also handles native attached records and cached lines.
            # A consumer plan (a Surface producer's lower-dimensional input)
            # reuses the stored `GEdge::mesh_vertices` discretization
            # verbatim: `meshGFace` upstream never re-runs `meshGEdge`, so a
            # retained `curve_params` entry is emitted, not regraded.
            keep=_dim01_fully_discrete(m,(1,curve)) ||
                 Model._model_curve_degenerate_kind(m,curve)==1 ||
                 (options.mesh_only_visible && model_entity_visibility(m,1,curve)==0) ||
                 (options.mesh_only_empty && owner in have_cells) ||
                 (reuse_curve_params && haskey(m.curve_params,curve))
            if keep
                if done!==nothing && !_dim01_fully_discrete(m,(1,curve)) &&
                   Model._model_curve_degenerate_kind(m,curve)==1
                    push!(done,owner)
                end
                if !(owner in have_cells) &&
                   haskey(m.curve_params,curve)
                    _dim01_emit_curve!(assembly,m,curve,m.curve_params[curve],first_nodes;
                                        reuse=true,high_order=high_order,owned_nodes=owned_nodes,curve_births)
                    push!(have_cells,owner)
                end
                delete!(pending,curve)
                continue
            end
            work=prepared===nothing ? nothing : prepared[index]
            if work!==nothing && parallel_results[index]===nothing
                vend=index
                while vend<length(ordered) && prepared[vend+1]!==nothing
                    vend+=1
                end
                _parallel_curve_fill!(parallel_results,index:vend,m,prepared,
                    options,parallel_dirty,parallel_precise,assembly.caller)
            end
            refs=index<parallel_last ? _parallel_surface_param_refs(m) : nothing
            result=work!==nothing && parallel_results[index]!==nothing ?
                _parallel_curve_result(m,work,parallel_results[index],
                    parallel_dirty,options,assembly.caller) :
                Model._model_mesh_curve!(m,curve,options,assembly.caller;
                    _native_producer=true)
            if result===:pending
                refs===nothing ||
                    _parallel_volume_mark_dirty!(parallel_dirty,m,refs)
                continue
            end
            if result===:keep
                if !(owner in have_cells) &&
                   haskey(m.curve_params,curve)
                    _dim01_emit_curve!(assembly,m,curve,m.curve_params[curve],first_nodes;
                                        reuse=true,high_order=high_order,owned_nodes=owned_nodes,curve_births)
                    push!(have_cells,owner)
                end
                delete!(pending,curve)
                refs===nothing ||
                    _parallel_volume_mark_dirty!(parallel_dirty,m,refs)
                continue
            end
            birth=result isa Model._ModelNativeCurveGrade ? result : nothing
            params=birth===nothing ? result : birth.parameters
            params=_dim01_native_periodic_parameters(assembly,m,curve,params,curve_births)
            m.curve_params[curve]=params
            done===nothing || push!(done,owner)
            push!(removed,owner)
            geometry_parameters=geometry_birth===nothing ? params :
                _dim01_geometry_curve_parameters!(geometry_birth,m,curve,params,options,assembly.caller)
            copied_points=_dim01_extruded_curve_points(assembly,m,curve,params,curve_births)
            _dim01_emit_curve!(assembly,m,curve,params,first_nodes;
                geometry_birth,geometry_parameters,curve_births,endpoint_counts,birth,copied_points)
            push!(have_cells,owner)
            delete!(pending,curve)
            refs===nothing || _parallel_volume_mark_dirty!(parallel_dirty,m,refs)
        end
        isempty(pending) && break
        if length(pending)==stalled
            # A stalled pass whose remainder contains a periodic dependency
            # cycle can never resolve by copying — upstream's retry starves
            # it. Mesh one seed member natively per cycle so the ordinary
            # deferred copy converges the rest; the seed keeps its own
            # grading rather than `_dim01_native_periodic_parameters`'s
            # master-node projection (no master nodes exist yet).
            for seed in Model._model_periodic_cycle_seeds(m,pending)
                owner=(1,Int32(seed))
                result=Model._model_mesh_curve!(m,seed,options,assembly.caller;
                    _native_producer=true,_force_native=true)
                result===:pending && continue
                if result===:keep
                    if !(owner in have_cells) && haskey(m.curve_params,seed)
                        _dim01_emit_curve!(assembly,m,seed,m.curve_params[seed],first_nodes;
                            reuse=true,high_order=high_order,owned_nodes=owned_nodes,curve_births)
                        push!(have_cells,owner)
                    end
                    delete!(pending,seed)
                    continue
                end
                birth=result isa Model._ModelNativeCurveGrade ? result : nothing
                params=birth===nothing ? result : birth.parameters
                m.curve_params[seed]=params
                done===nothing || push!(done,owner)
                push!(removed,owner)
                geometry_parameters=geometry_birth===nothing ? params :
                    _dim01_geometry_curve_parameters!(geometry_birth,m,seed,params,options,assembly.caller)
                copied_points=_dim01_extruded_curve_points(assembly,m,seed,params,curve_births)
                _dim01_emit_curve!(assembly,m,seed,params,first_nodes;
                    geometry_birth,geometry_parameters,curve_births,endpoint_counts,birth,copied_points)
                push!(have_cells,owner)
                delete!(pending,seed)
            end
            length(pending)<stalled && continue
        end
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

function _dim01_discrete_curve_parameter(assembly,m,curve,geometry,position,point_positions;
                                       periodic_end=false)
    node=assembly.nodes[position]
    record=m.discrete[(1,curve)]
    boundary=record.boundary
    if !isempty(boundary)
        begin_point=first(boundary)
        end_point=last(boundary)
        begin_position=get(point_positions,(begin_point[1],Int32(abs(begin_point[2]))),Int32(0))
        end_position=get(point_positions,(end_point[1],Int32(abs(end_point[2]))),Int32(0))
        periodic=begin_point[1]==end_point[1] && abs(begin_point[2])==abs(end_point[2])
        periodic_end && periodic && position==end_position && return geometry.t1
        position==begin_position && return geometry.t0
        position==end_position && return geometry.t1
    end
    if node.parameter!==nothing && !isempty(node.parameter) && node.parameter[1]!=1e6
        parameter=node.parameter[1]
    else
        _,parameter=Model._model_discrete_curve_xyz_to_u(
            geometry,node.point,Model._model_lc(m),assembly.caller)
    end
    parameter<1e6 || throw(ArgumentError(
        "$(assembly.caller): cannot reparametrize node $(node.tag) on discrete Curve[$curve]"))
    return parameter
end

# Endpoint identity uses the first mesh vertex on the geometric Point.
# Original MEdge/MFace parameters remain independent of the current mesh XYZ.
function _api_p2_curve_parameter(m,curve,geometry,position,point,parameter,point_positions;
                               periodic_end=false)
    endpoints=get(m.curves,curve,nothing)
    if endpoints===nothing
        record=get(m.discrete,(1,curve),nothing)
        boundary=record===nothing ? Tuple{Int,Int}[] : record.boundary
        endpoints=length(boundary)>=2 ? (abs(first(boundary)[2]),abs(last(boundary)[2])) : nothing
    end
    if endpoints!==nothing
        first_position=get(point_positions,(0,Int32(endpoints[1])),Int32(0))
        last_position=get(point_positions,(0,Int32(endpoints[2])),Int32(0))
        periodic_end && endpoints[1]==endpoints[2] && position==last_position && return geometry.upper
        position==first_position && return geometry.lower
        position==last_position && return geometry.upper
    end
    u=parameter!==nothing && !isempty(parameter) && parameter[1]!=1e6 ? parameter[1] :
        Model._model_high_order_curve_inverse(m,curve,geometry,point,"API.mesh high-order Curve")
    u<1e6 || throw(ArgumentError("API.mesh: cannot reparametrize support node on Curve[$curve]"))
    return u
end

# Captured discrete CAD can disagree with edited mesh vertices or their original
# parameters. HighOrder.cpp still constructs those exact supports, including a
# folded Line map. This admission concerns that bound provider; other generated
# P2 cells retain the FEM whole-domain certificate.
function _dim01_captured_curve_p2(assembly,m,cell)
    cell.msh==8 && cell.owner[1]==1 || return false
    entity=(1,Int(cell.owner[2]))
    record=get(m.discrete,entity,nothing)
    record!==nothing && record.geometry_parametrized &&
        Model._model_discrete_curve_geometry(m,entity[2])!==nothing || return false
    for position in cell.nodes
        node=assembly.nodes[position]
        get(assembly.node_positions,node.tag,Int32(0))==position || throw(ArgumentError(
            "$(assembly.caller): captured Curve[$(entity[2])] loses its current node identity"))
        (node.owner==cell.owner || node.owner[1]==0 &&
            any(entry->entry[1]==0 && abs(entry[2])==node.owner[2],
                record.boundary)) || throw(ArgumentError(
            "$(assembly.caller): captured Curve[$(entity[2])] references an unrelated node owner"))
        all(isfinite,node.point) || throw(ArgumentError(
            "$(assembly.caller): captured Curve[$(entity[2])] has non-finite coordinates"))
    end
    support=assembly.nodes[last(cell.nodes)]
    support.owner==cell.owner && support.parameter!==nothing &&
        length(support.parameter)==1 && all(isfinite,support.parameter) ||
        throw(ArgumentError("$(assembly.caller): captured Curve[$(entity[2])] has an invalid support parameter"))
    return true
end

# Receipt for a support in this detached native CAD order operation. The
# public UID and current assembly row establish identity; saved coordinates
# and parameters validate that the recorded support has not been replaced.
struct _Dim01NativeSupportReceipt
    tag::UInt64
    owner::Tuple{Int,Int32}
    point::NTuple{3,Float64}
    parameter::Union{Nothing,Vector{Float64}}
    carriers::Vector{Int32}
end

function _dim01_native_support_receipt!(receipts,assembly,position,carriers)
    node=assembly.nodes[position]
    get(assembly.node_positions,node.tag,Int32(0))==position || throw(ArgumentError(
        "$(assembly.caller): native CAD support loses its current node identity"))
    all(isfinite,node.point) && (node.parameter===nothing || all(isfinite,node.parameter)) ||
        throw(ArgumentError("$(assembly.caller): native CAD support is non-finite"))
    receipts[position]=_Dim01NativeSupportReceipt(node.tag,node.owner,node.point,
        node.parameter===nothing ? nothing : copy(node.parameter),Int32[carriers...])
    return nothing
end

# HighOrder.cpp uses plain MVertex interpolation when SecondOrderLinear is
# selected. A native CAD entity still owns that construction even if rounding
# folds its map; absence of stored parameters is part of the native result.
function _dim01_native_interpolated_owner(m,owner)
    dimension,tag=owner
    dimension in (1,2) || return false
    _dim01_fully_discrete(m,(dimension,Int(tag))) && return false
    if dimension==1
        return haskey(m.curves,Int(tag)) ||
            Model._model_discrete_curve_geometry(m,Int(tag))!==nothing
    end
    return haskey(m.surfaces,Int(tag))
end

function _dim01_native_support_seeds!(receipts,assembly,m,only_visible)
    for cell in assembly.cells
        spec=msh_spec(cell.msh)
        spec.order==2 && cell.owner[1] in (1,2) || continue
        only_visible && model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))==0 && continue
        _dim01_fully_discrete(m,(cell.owner[1],Int(cell.owner[2]))) && continue
        primary=Elements.msh_num_nodes(msh_type(spec.family,1))
        reference=Elements.lagrange_nodes(cell.msh)
        for slot in primary+1:spec.nnodes
            weights=MeshFunctionSpaces._first_order_values(Val(spec.family),
                reference[1,slot],reference[2,slot],reference[3,slot],assembly.caller,slot)
            key=_dim01_support_key(@view(cell.nodes[1:primary]),weights)
            carriers=Int32[first(term) for term in key]
            _dim01_native_support_receipt!(receipts,assembly,cell.nodes[slot],carriers)
        end
    end
    return nothing
end

function _dim01_native_constructed_cell_receipt(assembly,receipts,cell,closures)
    cell.owner[1] in (1,2) || return false
    spec=msh_spec(cell.msh)
    primary=Elements.msh_num_nodes(msh_type(spec.family,1))
    original=get(closures,(cell.owner[1],Int(cell.owner[2])),nothing)
    original===nothing && return false
    pending=Int32[cell.nodes[primary+1:end]...]
    visited=Set{Int32}()
    while !isempty(pending)
        position=pop!(pending)
        position in visited && continue
        push!(visited,position)
        receipt=get(receipts,position,nothing)
        receipt===nothing && return false
        node=assembly.nodes[position]
        get(assembly.node_positions,receipt.tag,Int32(0))==position &&
            node.tag==receipt.tag && node.owner==receipt.owner &&
            isequal(node.point,receipt.point) && isequal(node.parameter,receipt.parameter) ||
            throw(ArgumentError("$(assembly.caller): native CAD support receipt changed before publication"))
        for carrier in receipt.carriers
            carrier in original || push!(pending,carrier)
        end
    end
    return true
end

# All CAD and source records are fixed during this detached P2 operation.
# Reuse its already established Plane frame, preserving public validation and
# the distinct EDGE CAD-coordinate / FACE edited-native-Point semantics.
function _dim01_native_plane_frame_reusable(m,surface)
    haskey(m.surfaces,surface) || return false
    haskey(m.discrete,(2,surface)) && return false
    Model._surface_type(m,surface)==:plane || return false
    geometry=get(m.surface_geometry,surface,nothing)
    return !(geometry!==nothing && hasproperty(geometry,:occ))
end

function _dim01_native_plane_parametrization_cached!(frames,m,surface,point)
    caller="model_parametrization"
    _,entity_tag=Model._model_evaluation_entity(m,2,surface,(0,1,2),caller)
    values=Model._model_evaluation_values(collect(point),3,caller,"coordinates")
    coordinate=(values[1],values[2],values[3])
    return _dim01_surface_support_parameters!(frames,m,entity_tag,coordinate,caller)
end

function _dim01_native_plane_surface_uv_cached!(frames,m,assembly,node,surface)
    source=assembly.nodes[node]
    dim,entity=source.owner
    parameter=source.parameter
    dim==2 && entity==surface && parameter!==nothing && length(parameter)>=2 &&
        return parameter[1:2]
    tag=Int(surface)
    _dim01_native_plane_frame_reusable(m,tag) ||
        return _native_order_surface_uv(m,assembly,node,surface)
    if dim==0 && haskey(m.points,Int(entity))
        input=Float64[]
    elseif dim==1 && haskey(m.curves,Int(entity))
        parameter!==nothing && !isempty(parameter) || return nothing
        input=parameter[1:1]
    else
        return _dim01_native_plane_parametrization_cached!(frames,m,tag,source.point)
    end
    caller="model_reparametrize_on_surface"
    dimension,entity_tag=Model._model_evaluation_entity(m,dim,Int(entity),(0,1),caller)
    _,plane_tag=Model._model_metadata_entity(m,2,tag,caller)
    Model._model_reparametrization_selector(0,caller)
    values=Model._model_evaluation_values(input,1,caller,"parametric coordinates")
    # Public reparametrization establishes its frame after validation and
    # before evaluating the CAD Point/Curve coordinate, even on a cold cache.
    if !haskey(frames,plane_tag)
        frames[plane_tag]=Model._model_plane_frame(m,plane_tag,caller)
    end
    if dimension==0
        isempty(values) || throw(ArgumentError(
            "$caller: Point parametric coordinates must be empty"))
        coordinate=m.points[entity_tag]
    else
        coordinate=Model._model_curve_point(m,entity_tag,only(values),caller,1)
    end
    return _dim01_surface_support_parameters!(frames,m,plane_tag,coordinate,caller)
end

function _dim01_native_plane_surface_face_uv_cached!(frames,m,assembly,node,surface)
    tag=Int(surface)
    _dim01_native_plane_frame_reusable(m,tag) ||
        return _native_order_surface_face_uv(m,assembly,node,surface)
    source=assembly.nodes[node]
    dimension,entity=source.owner
    if dimension==0 && haskey(m.points,Int(entity)) &&
            Model.model_entity_type(m,2,tag)=="Plane"
        kernel=m.meshing.native_geo_kernel
        origin=get(kernel.point_aliases,Int(entity),Int(entity))
        if haskey(kernel.points,origin)
            return _dim01_native_plane_parametrization_cached!(frames,m,tag,source.point)
        end
    end
    return _dim01_native_plane_surface_uv_cached!(frames,m,assembly,node,surface)
end

function _dim01_native_plane_support!(frames,m,assembly,cell,primary,slot,owner,point)
    surface=Int(owner[2])
    fallback=_dim01_surface_support_parameters!(frames,m,surface,point,assembly.caller)
    cell.owner==owner && Model.model_entity_type(m,2,surface)=="Plane" ||
        return point,fallback,Int32[],false
    family=msh_spec(cell.msh).family
    edges=_native_order_edges(family)
    local_slot=slot-length(primary)
    if local_slot<=length(edges)
        first,last=edges[local_slot]
        first=primary[first+1];last=primary[last+1]
        u0=_dim01_native_plane_surface_uv_cached!(frames,m,assembly,first,owner[2])
        u1=_dim01_native_plane_surface_uv_cached!(frames,m,assembly,last,owner[2])
        if u0===nothing || u1===nothing
            # reparamMeshEdgeOnFace falls back for both endpoints together.
            u0=_dim01_surface_support_parameters!(frames,m,surface,assembly.nodes[first].point,assembly.caller)
            u1=_dim01_surface_support_parameters!(frames,m,surface,assembly.nodes[last].point,assembly.caller)
        end
        parameter=Float64[u0[axis]+(u1[axis]-u0[axis])*0.5 for axis in 1:2]
        carriers=Int32[first,last]
    elseif family===:qua
        row=only(_native_inner_rows(:qua,2))
        boundary=cell.nodes
        uv=[_dim01_native_plane_surface_face_uv_cached!(frames,m,assembly,node,owner[2]) for node in boundary]
        if all(value->value!==nothing,uv)
            parameter=zeros(Float64,2)
            for (column,weight) in row
                parameter[1]+=weight*uv[column][1]
                parameter[2]+=weight*uv[column][2]
            end
        else
            # A failed vertex reparameterization uses the physical native
            # boundary stencil before the Plane closest-point projection.
            physical=(0.0,0.0,0.0)
            for (column,weight) in row
                p=assembly.nodes[boundary[column]].point
                physical=(physical[1]+weight*p[1],physical[2]+weight*p[2],physical[3]+weight*p[3])
            end
            parameter=_dim01_surface_support_parameters!(frames,m,surface,physical,assembly.caller)
        end
        carriers=copy(boundary)
    else
        return point,fallback,Int32[],false
    end
    geometry=get(m.surface_geometry,surface,nothing)
    placed=if geometry!==nothing && hasproperty(geometry,:occ)
        values=model_value(m,2,surface,parameter)
        (values[1],values[2],values[3])
    else
        Model._model_plane_point(frames[surface],parameter[1],parameter[2],assembly.caller,1)
    end
    return placed,parameter,carriers,true
end

function _dim01_surface_support_parameters!(frames,m,surface,point,caller)
    geometry=get(m.surface_geometry,surface,nothing)
    if !haskey(m.discrete,(2,surface)) && Model._surface_type(m,surface)==:plane &&
            !(geometry!==nothing && hasproperty(geometry,:occ))
        # All supports in this detached order operation see the same CAD and
        # source records. Compute its native Plane frame once; evaluating the
        # same frame retains the public parametrization's literal arithmetic.
        frame=if haskey(frames,surface)
            frames[surface]
        else
            value=Model._model_plane_frame(m,surface,caller)
            frames[surface]=value
            value
        end
        return collect(Model._model_plane_parameters(frame,point,caller,1))
    end
    return Model.model_parametrization(m,2,surface,collect(point))
end

# A Surface factory consumes the actual EDGE endpoints or the eight produced
# FACE boundary slots. CAD boundary incidence alone does not create a lower
# mesh factory: shared supports are reused through their existing UID key.
function _dim01_native_ruled_support!(m,assembly,cell,primary,slot,owner)
    surface=Int(owner[2])
    owner==cell.owner && owner[1]==2 &&
        Model._model_high_order_ruled_surface(m,surface,assembly.caller) ||
        return nothing
    family=msh_spec(cell.msh).family
    edges=_native_order_edges(family)
    local_slot=slot-length(primary)
    if local_slot<=length(edges)
        a,b=edges[local_slot]
        first,last=primary[a+1],primary[b+1]
        u0=_native_order_surface_uv(m,assembly,first,owner[2])
        u1=_native_order_surface_uv(m,assembly,last,owner[2])
        initial=u0===nothing || u1===nothing ? nothing :
            Float64[u0[axis]+(u1[axis]-u0[axis])*0.5 for axis in 1:2]
        a=assembly.nodes[first].point;b=assembly.nodes[last].point
        physical=ntuple(axis->0.5*b[axis]+0.5*a[axis],3)
        carriers=Int32[first,last]
    elseif family===:qua && slot==9
        length(cell.nodes)==8 || throw(ArgumentError(
            "$(assembly.caller): native Quad2 center lacks its produced boundary slots"))
        physical=(0.0,0.0,0.0)
        uv=[_native_order_surface_face_uv(m,assembly,node,owner[2]) for node in cell.nodes]
        initial=all(value->value!==nothing,uv) ? zeros(Float64,2) : nothing
        for (column,weight) in only(_native_inner_rows(:qua,2))
            p=assembly.nodes[cell.nodes[column]].point
            physical=(physical[1]+weight*p[1],physical[2]+weight*p[2],physical[3]+weight*p[3])
            if initial!==nothing
                initial[1]+=weight*uv[column][1]
                initial[2]+=weight*uv[column][2]
            end
        end
        carriers=copy(cell.nodes)
    else
        return nothing
    end
    point,parameter,_=Model._model_closest_ruled_native(
        m,surface,physical,assembly.caller;initial)
    return point,collect(parameter),carriers,true
end

# HighOrder.cpp's Curve edge container rejects repeated non-degenerate
# endpoint identities when constructing new edge supports. Keep that native
# admission separate from coordinate coincidence and the FEM Jacobian check.
function _dim01_p2_curve_edge_guard(assembly,m,only_visible)
    selected=0
    for cell in assembly.cells
        cell.owner[1]==1 && msh_spec(cell.msh).family===:lin &&
            msh_spec(cell.msh).order==1 || continue
        only_visible && model_entity_visibility(m,1,Int(cell.owner[2]))==0 && continue
        cell.nodes[1]!=cell.nodes[2] && (selected+=1)
    end
    selected<2 && return nothing
    occupied=Dict{Tuple{Int32,Int32},Int32}()
    sizehint!(occupied,selected)
    for cell in assembly.cells
        cell.owner[1]==1 && msh_spec(cell.msh).family===:lin &&
            msh_spec(cell.msh).order==1 || continue
        only_visible && model_entity_visibility(m,1,Int(cell.owner[2]))==0 && continue
        first,last=cell.nodes[1],cell.nodes[2]
        first==last && continue
        edge=minmax(first,last)
        haskey(occupied,edge) && throw(ArgumentError(
            "$(assembly.caller): Mesh edges from different curves share nodes: "*
            "create a finer mesh (curve involved: $(cell.owner[2]))"))
        occupied[edge]=cell.owner[2]
    end
    return nothing
end

# Builtin Plane Quad2 centers retain their independent geometry-carrier
# parameter receipt while actual placement shares native FACE semantics.
function _dim01_native_quad2_center(assembly,m,cell,slot,linear,point,frames)
    linear && return nothing
    cell.owner[1]==2 && msh_family(cell.msh)===:qua && slot==9 || return nothing
    surface=Int(cell.owner[2])
    haskey(m.surfaces,surface) && !haskey(m.discrete,(2,surface)) || return nothing
    geometry=get(m.surface_geometry,surface,nothing)
    Model._surface_type(m,surface)===:plane &&
        !(geometry!==nothing && hasproperty(geometry,:occ)) || return nothing
    length(cell.nodes)==8 || throw(ArgumentError(
        "$(assembly.caller): native Quad2 center lacks its eight produced boundary slots"))
    return _dim01_native_plane_support!(frames,m,assembly,cell,
        @view(cell.nodes[1:4]),slot,cell.owner,point)
end

# OCC analytic revolution faces (Cylinder/Sphere/Cone/Torus records) project
# through the same `model_closest_point`/`model_value` evaluators the order>=3
# native face factories already consume. Plane records keep the Plane path and
# NURBS-carrying records keep the NURBS arm below.
_occ_p2_curved_surface(geometry)=geometry!==nothing && hasproperty(geometry,:occ) &&
    geometry.occ in (:cylinder,:sphere,:cone,:torus)

function _dim01_quadratic!(assembly,m,entities;only_visible=false,linear=false,geometry_birth=nothing)
    _dim01_p2_curve_edge_guard(assembly,m,only_visible)
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
    if !linear && !isempty(native_entities)
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
    curve_geometry=Dict{Int,Any}()
    surface_frames=Dict{Int,Any}()
    native_supports=Dict{Int32,_Dim01NativeSupportReceipt}()
    receipt_capacity=0
    for original_cell in assembly.cells
        original_spec=msh_spec(original_cell.msh)
        original_cell.owner[1] in (1,2) && original_spec.order in (1,2) || continue
        target_spec=msh_spec(msh_type(original_spec.family,2))
        receipt_capacity+=target_spec.nnodes-Elements.msh_num_nodes(msh_type(original_spec.family,1))
    end
    sizehint!(native_supports,receipt_capacity)
    _dim01_native_support_seeds!(native_supports,assembly,m,only_visible)
    point_positions=Dict{Tuple{Int,Int32},Int32}()
    for (position,node) in enumerate(assembly.nodes)
        node.owner[1]==0 && get!(point_positions,node.owner,Int32(position))
    end
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
                if !(cell.owner[1]==2 &&
                        Model._model_high_order_ruled_surface(m,Int(cell.owner[2]),assembly.caller))
                    for candidate in get(incident_entities,key[1][1],Tuple{Int,Int}[])
                        candidate[1]>cell.owner[1] && break
                        all(term->term[1] in closures[candidate],key) || continue
                        owner=(candidate[1],Int32(candidate[2]))
                        break
                    end
                end
                parameter=nothing
                legacy_surface_parameter=nothing
                center=nothing
                native_placed=linear && _dim01_native_interpolated_owner(m,owner)
                receipt_carriers=Int32[first(term) for term in key]
                if !linear && owner[1]==1
                    curve=Int(owner[2])
                    geometry=get!(curve_geometry,curve) do
                        Model._model_high_order_curve_geometry(m,curve,assembly.caller)
                    end
                    if geometry!==nothing
                        endpoints=spec.family===:lin && cell.owner==owner ? primary :
                            Int32[key[1][1],key[2][1]]
                        first_node=assembly.nodes[endpoints[1]];last_node=assembly.nodes[endpoints[2]]
                        u0=_api_p2_curve_parameter(m,curve,geometry,endpoints[1],
                            first_node.point,first_node.owner==owner ? first_node.parameter : nothing,
                            point_positions)
                        u1=_api_p2_curve_parameter(m,curve,geometry,endpoints[2],
                            last_node.point,last_node.owner==owner ? last_node.parameter : nothing,
                            point_positions;periodic_end=true)
                        u=Model._model_high_order_curve_p2_parameter(geometry,u0,u1,assembly.caller)
                        parameter=Float64[u]
                        point=Model._model_high_order_curve_point(geometry,u,assembly.caller)
                        native_placed=true
                    end
                elseif !linear && owner[1]==2 && haskey(m.surfaces,Int(owner[2]))
                    ruled=_dim01_native_ruled_support!(m,assembly,cell,primary,slot,owner)
                    if ruled!==nothing
                        point,parameter,receipt_carriers,native_placed=ruled
                    elseif Model._model_nurbs_surface(get(m.surface_geometry,Int(owner[2]),nothing))
                        point,parameter=Model._model_nurbs_surface_support(
                            m,Int(owner[2]),assembly.nodes,key,assembly.caller;
                            quadrangle_boundary=spec.family===:qua && slot==9 ? cell.nodes : nothing)
                    elseif _occ_p2_curved_surface(get(m.surface_geometry,Int(owner[2]),nothing))
                        # Gmsh's high-order pass owns GFace supports by
                        # evaluating the face at the support's parameter;
                        # `closestPoint` projects the reference interpolation
                        # onto the analytic face and reports its (u,v).
                        closest,projected=Model.model_closest_point(
                            m,2,Int(owner[2]),collect(point))
                        point=(closest[1],closest[2],closest[3])
                        parameter=projected
                        native_placed=true
                    else
                        center=_dim01_native_quad2_center(
                            assembly,m,cell,slot,linear,point,surface_frames)
                        if center===nothing
                            point,parameter,receipt_carriers,native_placed=_dim01_native_plane_support!(
                                surface_frames,m,assembly,cell,primary,slot,owner,point)
                        else
                            # Geometry carrier keeps its original first-order parameter;
                            # actual Quad2 placement consumes native FACE boundary values.
                            legacy_surface_parameter=_dim01_surface_support_parameters!(
                                surface_frames,m,Int(owner[2]),point,assembly.caller)
                            point,parameter,receipt_carriers,native_placed=center
                        end
                    end
                end
                node=_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),point,owner,parameter)
                native_placed && _dim01_native_support_receipt!(
                    native_supports,assembly,node,receipt_carriers)
                geometry_birth===nothing || _dim01_geometry_support!(geometry_birth,assembly,m,node,key;linear)
                if center!==nothing && geometry_birth!==nothing &&
                        haskey(geometry_birth.nodes,assembly.nodes[node].tag)
                    retained=geometry_birth.nodes[assembly.nodes[node].tag]
                    geometry_birth.nodes[retained.tag]=_Dim01Node(retained.tag,
                        retained.point,retained.owner,legacy_surface_parameter)
                end
                support_nodes[key]=node
            end
            push!(cell.nodes,node)
        end
        cell.msh=target
        cell.tag=_dim01_next_tag!(assembly,false)
        (!linear && _dim01_captured_curve_p2(assembly,m,cell)) ||
            _dim01_native_constructed_cell_receipt(assembly,native_supports,cell,closures) ||
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
        boundary=if dim==0
            Int32[]
        elseif dim==1 && haskey(m.curves,tag)
            first_point,last_point=m.curves[tag]
            Int32[first_point,-last_point]
        else
            Int32[t for (_,t) in model_boundary(m,[entity],false,true,false)]
        end
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
        metadata[entity]=Model._model_mesh_entity_metadata(m,dim,tag,bounds;
            physical_tags=Model._model_projection_physical_tags(m,dim,tag),
            boundaries=boundary,embedded_curves=embedded)
    end
    return metadata,boundaries
end

function _dim01_finish!(assembly,m,entities,discarded,high_order,renumber,p2_pass;
                        save_all=false,node_limit=assembly.max_nodes,
                        cell_limit=assembly.max_cells,source_mesh=nothing,
                        source_node_positions=Int32[],curve_parameter_order::Bool=false,
                        node_tag_order::Bool=false,record_view::Bool=false,curve_births=nothing,
                        native_orphan_births=nothing,periodic_parent_receipt=nothing,
                        periodic_uid_operation::Symbol=:preserve)
    periodic_uid_operation in (:preserve,:order,:refine,:defer) || throw(ArgumentError(
        "$(assembly.caller): invalid periodic UID publication operation"))
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
            assembly.nodes[position]=_Dim01Node(node.tag,node.point,owner,node.parameter)
        end
    end
    # Initial SetOrder1 removes old native interpolation vertices. Remeshing
    # deletes native entity-owned vertices; retained discrete references survive.
    kept=Int[]
    sizehint!(kept,count(referenced))
    for (i,node) in enumerate(assembly.nodes)
        if referenced[i] || (!p2_pass && ((native_orphan_births!==nothing && node.tag in native_orphan_births) ||
            (!(node.owner in discarded) && !(i<=length(high_order) && high_order[i]))))
            push!(kept,i)
        end
    end
    length(kept)<=node_limit || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_nodes=$node_limit"))
    length(assembly.cells)<=cell_limit || throw(ArgumentError(
        "$(assembly.caller): mesh exceeds max_cells=$cell_limit"))
    if curve_parameter_order
        # Subdivide(GEdge) compares MVertex::getParameter(0). A plain
        # interpolation vertex returns zero and keeps its empty native data.
        sort!(kept;by=i->begin
            node=assembly.nodes[i]
            parameter=node.parameter
            u=node.owner[1]==1 && parameter!==nothing && !isempty(parameter) ? parameter[1] : 0.0
            (node.owner,u)
        end,alg=Base.Sort.MergeSort)
    elseif p2_pass || node_tag_order
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
        if _native_has_ghost_entities(m)
            ghost_owners=_imported_ghost_uid_owners(m)
            ghost_ranks=_imported_ghost_uid_ranks(m)
            occurrences=Tuple{Tuple{Int,Int32},Int,Int,Int}[]
            unchanged=Bool[cell.source_tag!==nothing && cell.tag==cell.source_tag for cell in assembly.cells]
            for (position,cell) in enumerate(assembly.cells)
                receipt=get(m.meshing.imported_entities,(cell.owner[1],Int(cell.owner[2])),nothing)
                if receipt===nothing || !startswith(receipt.kind,"Ghost ")
                    push!(occurrences,(cell.owner,_native_mesh_element_bucket(cell.msh),position,position))
                end
                for owner in get(ghost_owners,cell.tag,Tuple{Int,Int32}[])
                    owner[1]==cell.owner[1] || continue
                    rank=ghost_ranks[(owner[1],Int(owner[2]))][cell.tag]
                    push!(occurrences,(owner,_native_mesh_element_bucket(cell.msh),rank,position))
                end
            end
            sort!(occurrences;alg=Base.Sort.MergeSort)
            for selected in (true,false),(owner,_,_,position) in occurrences
                saved=_dim01_physical(m,(owner[1],Int(owner[2])))!=0
                save_physical && saved!=selected && continue
                !save_physical && !selected && continue
                next+=UInt64(1)
                cell=assembly.cells[position]
                unchanged[position] && (visibility_tag_map[cell.source_tag]=next)
                cell.tag=next
            end
        else
            for selected in (true,false),cell in assembly.cells
                saved=_dim01_physical(m,(cell.owner[1],Int(cell.owner[2])))!=0
                save_physical && saved!=selected && continue
                !save_physical && !selected && continue
                next+=UInt64(1)
                cell.source_tag!==nothing && cell.tag==cell.source_tag &&
                    (visibility_tag_map[cell.source_tag]=next)
                cell.tag=next
            end
        end
    end
    order=Int[];groups=Dict{Int,Vector{Int}}()
    element_tag_map=Dict{UInt64,UInt64}()
    for (slot,cell) in enumerate(assembly.cells)
        cell.source_tag!==nothing && (element_tag_map[cell.source_tag]=cell.tag)
        !renumber && cell.source_tag!==nothing && cell.tag==cell.source_tag &&
            (visibility_tag_map[cell.source_tag]=cell.tag)
        haskey(groups,cell.msh) || push!(order,cell.msh)
        push!(get!(()->Int[],groups,cell.msh),slot)
    end
    Model._model_imported_ghost_remap!(m,visibility_tag_map;drop_removed=true)
    coords=Matrix{Float64}(undef,3,length(kept))
    node_entities=Tuple{Int,Int32}[]
    node_parameters=Union{Nothing,Vector{Float64}}[]
    sizehint!(node_entities,length(kept))
    sizehint!(node_parameters,length(kept))
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
    projected_parameters=Union{Nothing,Vector{Float64}}[
        Model._model_mesh_projected_parameter(owner,value)
        for (owner,value) in zip(node_entities,node_parameters)]
    data=Elements._native_mixed_entity_data(metadata;node_entities=node_entities,
        node_parametric=projected_parameters,external_node_tags=node_tags,
        block_entities=block_owners,external_element_tags=external_cells)
    mesh=MixedMesh(coords,blocks;physical_names=m.physical_names,entity_data=data,
                   elementary_entities=block_owners)
    classified_entities=Tuple{Int,Int32}[(dim,Int32(tag)) for (dim,tag) in entities]
    primary=isempty(classified_entities) ? (0,Int32(0)) : classified_entities[end]
    class=_mixed_classification(mesh,primary,classified_entities,node_entities,boundaries,owners;
        node_parametric=node_parameters)
    if record_view
        return (model=m,mesh=mesh,class=class,authority=Set(entities),
            node_tag_map=node_tag_map,element_tag_map=element_tag_map,
            visibility_tag_map=visibility_tag_map,
            max_node_tag=assembly.node_max,max_element_tag=assembly.cell_max,
            node_counter_wrapped=assembly.node_counter_wrapped,
            element_counter_wrapped=assembly.cell_counter_wrapped)
    end
    links=Elements.MixedPeriodicLink[]
    endpoint_constraints=Model.ModelPeriodicConstraint[]
    for ((dim,slave),constraint) in sort!(collect(m.periodic);by=first)
        dim==1 || continue
        push!(endpoint_constraints,constraint)
        correspondence=_dim01_periodic_nodes(m,mesh,class,slave,assembly.caller;
                                             include_high_order=true,source_mesh,node_tag_map,curve_births)
        push!(links,Elements.MixedPeriodicLink(1,slave,correspondence.master_entity,
            correspondence.slave_nodes,correspondence.master_nodes;
            affine=correspondence.affine))
    end
    emitted_endpoints=Set{Int}()
    if !isempty(endpoint_constraints)
        # The surface projection's periodic endpoint relations are part of the
        # actual merged product — without them a relation's corner vertices
        # serialize without their dim-0 master links.
        point_nodes=Dict{Int,Int32}()
        for (row,owner) in enumerate(node_entities)
            owner[1]==0 || continue
            point=Int(owner[2])
            haskey(point_nodes,point) && throw(ErrorException(
                "$(assembly.caller): periodic Point[$point] has multiple actual nodes"))
            point_nodes[point]=Int32(row)
        end
        parents=Model._model_periodic_endpoint_forest(
            m,endpoint_constraints,assembly.caller)
        for slave_point in sort!(collect(keys(parents)))
            master_point,affine=parents[slave_point]
            # A point whose entity was not generated has no actual node — like
            # a starved curve relation, its endpoint link is mesh-inert here.
            haskey(point_nodes,slave_point) && haskey(point_nodes,master_point) ||
                continue
            push!(links,Elements.MixedPeriodicLink(0,slave_point,master_point,
                Int32[point_nodes[slave_point]],Int32[point_nodes[master_point]];
                affine=affine))
            # MSH stores one master per slave vertex — like the dim-1 dedupe,
            # the emitted forest link supersedes a stored link for the same
            # slave endpoint even when the stored master differs.
            push!(emitted_endpoints,slave_point)
        end
    end
    if source_mesh isa MixedMesh
        for link in source_mesh.periodic_links
            link.dim==1 && haskey(m.periodic,(1,Int(link.slave_entity))) && continue
            link.dim==2 && continue
            link.dim==0 && Int(link.slave_entity) in emitted_endpoints &&
                continue
            slaves=mapping[source_node_positions[link.slave_nodes]]
            masters=mapping[source_node_positions[link.master_nodes]]
            all(!iszero,slaves) && all(!iszero,masters) || throw(ArgumentError(
                "$(assembly.caller): existing periodic dimension $(link.dim) " *
                "entity $(link.slave_entity) loses interpolation vertices; " *
                "higher-dimensional correspondence reconstruction is not implemented"))
            push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,
                link.master_entity,slaves,masters;affine=link.affine))
        end
        append!(links,filter(link->link.dim==2,
            _periodic_surface_order_links(source_mesh,mesh,class,assembly.caller;
                node_tag_map,element_tag_map,parent_receipt=periodic_parent_receipt)))
    end
    # A periodic slave that meshes through its own Transfinite run stores only
    # the relation marker — the node map materializes in the classified
    # geometry projection. A rebuild over a link-free carrier must re-emit that
    # marker the same way the dim-0/dim-1 relations above are re-derived.
    emitted_surfaces=Set{Int}(Int(link.slave_entity) for link in links if link.dim==2)
    for ((dim,slave),constraint) in sort!(collect(m.periodic);by=first)
        dim==2 || continue
        slave in emitted_surfaces && continue
        haskey(m.meshing.transfinite_surfaces,slave) || continue
        any(owner->owner==(2,Int32(slave)),node_entities) || continue
        push!(links,Elements.MixedPeriodicLink(2,constraint.slave_entity,
            constraint.master_entity,Int32[],Int32[];affine=constraint.affine))
    end
    if !isempty(links)
        mesh=MixedMesh(coords,blocks;physical_names=m.physical_names,entity_data=data,
                       elementary_entities=block_owners,periodic_links=links)
        class=_mixed_rebind_class(class,mesh)
    end
    if periodic_uid_operation===:preserve
        Model._model_periodic_uid_remap!(m,node_tag_map)
    elseif periodic_uid_operation===:order
        _periodic_uid_history_order!(m,mesh,class,node_tag_map,assembly.caller)
    elseif periodic_uid_operation===:refine
        _periodic_uid_history_refine!(m,source_mesh,mesh,class,node_tag_map,
            element_tag_map,periodic_parent_receipt,assembly.caller)
    end
    Model._model_native_curve_endpoint_remap!(m,node_tag_map)
    _dim01_reconcile_records!(m,entities,mesh,node_tag_map;runtime_parameters=node_parameters)
    return (model=m,mesh=mesh,class=class,authority=Set(entities),
        node_tag_map=node_tag_map,element_tag_map=element_tag_map,
        visibility_tag_map=visibility_tag_map,
        max_node_tag=assembly.node_max,max_element_tag=assembly.cell_max,
            node_counter_wrapped=assembly.node_counter_wrapped,
            element_counter_wrapped=assembly.cell_counter_wrapped)
end

function _dim01_reconcile_records!(m,entities,mesh,node_tag_map;runtime_parameters=nothing)
    data=mesh.entity_data
    positions_by_owner=Dict{Tuple{Int,Int},Vector{Int}}()
    for (position,owner) in enumerate(data.node_entities)
        push!(get!(()->Int[],positions_by_owner,(owner[1],Int(owner[2]))),position)
    end
    records=Dict{Tuple{Int,Int},DiscreteEntity}()
    for entity in entities
        dim,tag=entity
        previous=_dim01_record(m,entity)
        # Every mesh payload below is rebuilt from actual identity rows.
        # Copy the boundary and retained parameter registries without cloning
        # the node and cell arrays that this reconciliation replaces.
        record=previous===nothing ? DiscreteEntity() : DiscreteEntity(
            copy(previous.boundary),UInt64[],zeros(3,0),zeros(0,0),
            Int32[],UInt64[],Vector{UInt64}[],deepcopy(previous.aux_params),
            runtime_parameters===nothing ? deepcopy(previous.node_native_params) :
                Dict{UInt64,Vector{Float64}}(),previous.geometry_parametrized)
        positions=get(positions_by_owner,entity,Int[])
        record.node_tags=copy(data.external_node_tags[positions])
        record.node_coords=mesh.coords[:,positions]
        parameters=data.node_parametric[positions]
        record.node_params=all(value->value!==nothing && length(value)==dim,parameters) ?
            (isempty(parameters) ? zeros(dim,0) : hcat(parameters...)) : zeros(0,length(positions))
        if runtime_parameters!==nothing
            record.node_native_params=Dict(data.external_node_tags[position]=>
                (runtime_parameters[position]===nothing ? Float64[] : copy(runtime_parameters[position]))
                for position in positions)
        else
            rebound=Dict{UInt64,Vector{Float64}}()
            owned=Set(record.node_tags)
            for (tag,value) in record.node_native_params
                haskey(node_tag_map,tag) || continue
                mapped=node_tag_map[tag]
                mapped in owned || continue
                rebound[mapped]=copy(value)
            end
            record.node_native_params=rebound
        end
        empty!(record.element_types);empty!(record.element_tags);empty!(record.element_nodes)
        record.aux_params=Dict{UInt64,Vector{Float64}}(
            node_tag_map[UInt64(source)]=>values for (source,values) in record.aux_params
            if haskey(node_tag_map,UInt64(source)))
        records[entity]=record
    end
    for (bi,block) in enumerate(mesh.blocks)
        dim=msh_dimension(block.msh)
        for column in axes(block.nodes,2)
            record=get(records,(dim,Int(data.block_entities[bi][column])),nothing)
            record===nothing && continue
            push!(record.element_types,Int32(block.msh))
            push!(record.element_tags,data.external_element_tags[bi][column])
            push!(record.element_nodes,copy(data.external_node_tags[block.nodes[:,column]]))
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

# Synchronization preserves a surviving Curve's old element incidence even
# after it retires a Point MVertex. Admit only the recorded retired Point UIDs
# in a dimension-one transaction that will replace this Curve, retaining gates.
function _dim01_retired_curve_nodes(m,dimension,options)
    dimension==1 || return nothing
    retired=m.meshing.native_geo_kernel.retired_point_curves
    isempty(retired) && return nothing
    bound=Dict{Tuple{Int,Int32},Set{UInt64}}()
    for curve in keys(retired)
        haskey(m.curves,curve) || continue
        _dim01_fully_discrete(m,(1,curve)) && continue
        Model._model_curve_degenerate_kind(m,curve)==1 && continue
        options.mesh_only_visible && model_entity_visibility(m,1,curve)==0 && continue
        record=_dim01_record(m,(1,curve))
        if options.mesh_only_empty
            params=get(m.curve_params,curve,nothing)
            params!==nothing && !isempty(params) && continue
            record!==nothing && !isempty(record.element_tags) && continue
        end
        bound[(1,Int32(curve))]=retired[curve]
    end
    return bound
end

function _generate_dim01_plan(source_model,cached,cached_class,dimension::Int,
                              options::Model._ModelMesh1DOptions,caller;
                              element_order=1,renumber=true,linear=false,incomplete=false,
                              save_all=false,
                              max_nodes=typemax(Int32),max_cells=typemax(Int32),
                              initial_max_node_tag=UInt64(0),
                              initial_max_element_tag=UInt64(0),
                              initial_node_counter_wrapped::Bool=false,
                              initial_element_counter_wrapped::Bool=false,
                              generation_meshing=nothing,
                              reuse_curve_params::Bool=false)
    dimension in (0,1) || throw(ArgumentError("$caller: dim must be zero or one"))
    1<=element_order<=10 || throw(ArgumentError(
        "$caller: Curve mesh order $element_order is not supported (supported: 1 through 10)"))
    max_nodes isa Integer && !(max_nodes isa Bool) && 0<=max_nodes<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_nodes must be a nonnegative Int32-representable integer"))
    max_cells isa Integer && !(max_cells isa Bool) && 0<=max_cells<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_cells must be a nonnegative Int32-representable integer"))
    m=Model._model_mesh_transaction_copy(source_model)
    # API options describe this generation, while source-cache completion above
    # it still uses the attributes of the mesh being retained. Apply only after
    # staging so a rejected request leaves the live model/session untouched.
    generation_meshing===nothing ||
        _apply_generation_meshing_options!(m,generation_meshing)
    entities=model_entities(m)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller),
        initial_node_counter_wrapped,initial_element_counter_wrapped)
    retired_nodes=_dim01_retired_curve_nodes(m,dimension,options)
    source_node_positions=_dim01_ingest!(assembly,m,cached,cached_class;preserve_source_order=true,retired_nodes)
    previous=m.meshing.mesh_authority
    geometry_birth=_dim01_geometry_birth(previous)
    curve_births=Dict{Tuple{Int,Int32},Vector{UInt64}}()
    done=previous===nothing ? Set{Tuple{Int,Int32}}() : copy(previous.done)
    # Staging may retain vertices that the detached pass later deletes. Bound
    # its growth by source size plus requested result capacity; certify the
    # exact result counts before allocating final coordinate/connectivity arrays.
    assembly.max_nodes=min(Int(typemax(Int32)),Base.checked_add(length(assembly.nodes),Int(max_nodes)))
    assembly.max_cells=min(Int(typemax(Int32)),Base.checked_add(length(assembly.cells),Int(max_cells)))
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    high_order=_dim01_first_order!(assembly,m;retired_nodes)
    discarded=Set{Tuple{Int,Int32}}()
    if dimension==1
        for entity in entities
            entity[1]>=2 && !_dim01_fully_discrete(m,entity) &&
                push!(discarded,(entity[1],Int32(entity[2])))
        end
        filter!(cell->!(cell.owner in discarded),assembly.cells)
        filter!(owner->!(owner in discarded),done)
        first_nodes=_dim01_point_phase!(assembly,m,entities)
        union!(discarded,_dim01_curve_phase!(assembly,m,options,entities,first_nodes,high_order;done,geometry_birth,curve_births,reuse_curve_params))
        for (dim,tag) in discarded
            dim==1 && delete!(m.meshing.native_geo_kernel.retired_point_curves,Int(tag))
        end
    end
    # orientMeshGEdge runs on every generation, including dimension zero and
    # retained discrete lines. It reverses the current mesh when requested.
    for cell in assembly.cells
        cell.owner[1]==1 && get(m.meshing.reverse,(1,Int(cell.owner[2])),false) || continue
        cell.nodes=cell.nodes[_mixed_reverse_permutation(cell.msh)]
    end
    # The post-generation high-order pass does not run for a point-only mesh.
    p2_pass=element_order>1 && any(cell->cell.owner[1]>0,assembly.cells)
    if element_order>2
        _dim01_higher_order!(assembly,m,entities,element_order;
            only_visible=options.mesh_only_visible,linear,incomplete)
    elseif p2_pass
        _dim01_quadratic!(assembly,m,entities;only_visible=options.mesh_only_visible,linear,geometry_birth)
    end
    plan=_dim01_finish!(assembly,m,entities,discarded,high_order,renumber,p2_pass;
                         save_all=save_all,node_limit=Int(max_nodes),cell_limit=Int(max_cells),
                         source_mesh=cached,source_node_positions=source_node_positions,curve_births)
    return merge(plan,(mesh_done=done,mesh_dimension=dimension,
        geometry_carrier=_dim01_geometry_carrier(plan,geometry_birth),curve_births=curve_births))
end

# Immediate setOrder is distinct from the global generation option: it does
# not run Mesh0D/Mesh1D or the line orientation constraints.
function _set_order_dim01_plan(source_model,cached,cached_class,order::Int,caller;
                               renumber=true,save_all=false,only_visible=false,linear=false,incomplete=false,
                               max_nodes=typemax(Int32),max_cells=typemax(Int32),
                               initial_max_node_tag=UInt64(0),
                               initial_max_element_tag=UInt64(0))
    1<=order<=10 || throw(ArgumentError(
        "$caller: Curve mesh order $order is not supported (supported: 1 through 10)"))
    max_nodes isa Integer && !(max_nodes isa Bool) && 0<=max_nodes<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_nodes must be a nonnegative Int32-representable integer"))
    max_cells isa Integer && !(max_cells isa Bool) && 0<=max_cells<=typemax(Int32) ||
        throw(ArgumentError("$caller: max_cells must be a nonnegative Int32-representable integer"))
    m=Model._model_mesh_transaction_copy(source_model)
    entities=model_entities(m)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    source_node_positions=_dim01_ingest!(assembly,m,cached,cached_class;preserve_source_order=true)
    periodic_source=cached
    if !isempty(m.meshing.periodic_uid_history)
        if !(cached isa MixedMesh)
            source_view=_dim01_finish!(deepcopy(assembly),m,entities,
                Set{Tuple{Int,Int32}}(),falses(length(assembly.nodes)),false,false;
                record_view=true)
            periodic_source,source_class=_periodic_uid_history_source(
                m,source_view.mesh,source_view.class,caller)
            source_node_positions=Int32[assembly.node_positions[uid]
                for uid in periodic_source.entity_data.external_node_tags]
        else
            periodic_source,_=_periodic_uid_history_source(m,cached,cached_class,caller)
        end
    end
    assembly.max_nodes=min(Int(typemax(Int32)),Base.checked_add(length(assembly.nodes),Int(max_nodes)))
    assembly.max_cells=min(Int(typemax(Int32)),Base.checked_add(length(assembly.cells),Int(max_cells)))
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    high_order=order==1 ? _dim01_first_order!(assembly,m;skip_discrete=false) :
                         falses(length(assembly.nodes))
    geometry_birth=_dim01_geometry_birth(m.meshing.mesh_authority)
    if order>2 || (order==2 && (incomplete || any(cell->msh_spec(cell.msh).order>2,assembly.cells)))
        _dim01_higher_order!(assembly,m,entities,order;only_visible,linear,incomplete)
    elseif order==2
        _dim01_quadratic!(assembly,m,entities;only_visible,linear,geometry_birth)
    end
    plan=_dim01_finish!(assembly,m,entities,Set{Tuple{Int,Int32}}(),high_order,
                         renumber,order>1;save_all=save_all,
                         node_limit=Int(max_nodes),cell_limit=Int(max_cells),
                         source_mesh=periodic_source,source_node_positions=source_node_positions,
                         periodic_uid_operation=:order)
    previous=m.meshing.mesh_authority
    return previous!==nothing && previous.carrier isa MixedMesh ?
        merge(plan,(geometry_carrier=_dim01_geometry_carrier(plan,geometry_birth),)) : plan
end

# The surface projector's curve-edge catalog deliberately excludes line-only
# caches. This mapper consumes actual classified Line connectivity instead.
# It returns dense indices; the public API translates them through its tag table.
include("APIPeriodicUIDCorrespondence.jl")
include("APIPeriodicSurfaceOrder.jl")
include("APIPeriodicSurfaceRefineReceipt.jl")
include("APIPeriodicUIDHistory.jl")

function _dim01_periodic_nodes(m,mesh,class,slave::Int,caller;
                               include_high_order=false,source_mesh=nothing,node_tag_map=nothing,curve_births=nothing)
    constraint=get(m.periodic,(1,slave),nothing)
    constraint===nothing && return (master_entity=slave,slave_nodes=Int32[],
                                    master_nodes=Int32[],affine=nothing)
    bound=_dim01_periodic_uid_nodes(m,mesh,class,slave,caller;
        include_high_order,source_mesh,node_tag_map,curve_births)
    bound===nothing || return bound
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
