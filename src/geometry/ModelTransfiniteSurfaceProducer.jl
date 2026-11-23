# Transfinite boundary identity comes from the actual Curve Line chains and
# the kernel's tagged side events. Geometry interpolation shares those birth
# rows; coordinate coincidence never establishes a source identity.
struct _ModelSurfaceGeometryProduction{P,G}
    production::P
    geometry::G
end

function _model_transfinite_source_admitted(m,t,input,caller)
    input===nothing && return _model_transfinite_declared_counts_match(m,t,caller)
    loops=m.surfaces[t]
    length(loops)==1 || return true # The regular TF dispatcher owns this diagnostic.
    curves=Int[abs(signed) for signed in m.loops[only(loops)]
        if !(abs(signed) in m.meshing.degenerated)]
    length(curves)==4 || return true
    counts=Int[]
    for curve in curves
        product=_model_mesh_product(input,1,curve;caller)
        product===nothing && return true
        # computeEdgeLoops visits the current first CAD Point vertex followed
        # by every GEdge-owned vertex. P2 association pruning can move an old
        # endpoint onto the Curve without changing its Line-cell identity.
        push!(counts,2+count(==((1,Int32(curve))),product.entity_data.node_entities))
    end
    return counts[1]==counts[3] && counts[2]==counts[4]
end

# Without produced curves: Gmsh's MeshTransfiniteSurface gives up ("cannot
# be meshed using the transfinite algorithm") when opposite sides of a
# four-curve loop carry different node counts, and the surface is meshed by
# the unstructured algorithm instead. Counts come from the transfinite
# records (or a curve's own parameters); anything else keeps the transfinite
# path and its diagnostics.
function _model_transfinite_declared_counts_match(m,t,caller)
    loops=m.surfaces[t]
    length(loops)==1 || return true
    curves=Int[abs(signed) for signed in m.loops[only(loops)]
        if !(abs(signed) in m.meshing.degenerated)]
    length(curves)==4 || return true
    counts=Int[]
    for curve in curves
        spec=get(m.meshing.transfinite_curves,curve,nothing)
        if spec!==nothing
            push!(counts,_flexible_transfinite_nodes(m,spec.num_nodes,curve,caller))
        elseif haskey(m.curve_params,curve)
            push!(counts,length(m.curve_params[curve]))
        else
            return true
        end
    end
    return counts[1]==counts[3] && counts[2]==counts[4]
end

function _model_transfinite_side_keys(m,curves,junctions,caller)
    n=length(junctions)
    return ntuple(n) do side
        first=junctions[side];last=junctions[mod1(side+1,n)]
        matches=Int[]
        for signed in curves
            curve=abs(signed);a,b=m.curves[curve]
            a==first && b==last && push!(matches,curve)
            b==first && a==last && push!(matches,-curve)
        end
        length(matches)==1 || throw(ArgumentError(
            "$caller: transfinite side does not identify exactly one boundary Curve"))
        only(matches)
    end
end
_model_transfinite_side_tags(m,curves,junctions,caller)=
    map(abs,_model_transfinite_side_keys(m,curves,junctions,caller))

function _model_transfinite_boundary_events(kernel)
    if kernel isa Mesh
        return [(kernel.segs,kernel.seg_tag)]
    end
    return [(block.nodes,block.tags) for block in kernel.blocks if block.msh==1]
end

# Whether a boundary node of transfinite Surface `t` has a geometry birth
# position other than its actual coordinates: the only input in which the
# geometry layout of the face differs from its actual layout.
function _model_transfinite_geometry_moved(m,t,input,geometry_points,caller)
    for signed in m.loops[only(m.surfaces[t])]
        curve=abs(signed)
        curve in m.meshing.degenerated && continue
        source=_model_mesh_curve_births(input,m,curve,caller)
        source===nothing && continue
        data=source.mesh.entity_data;coords=source.mesh.coords
        for row in source.rows
            point=get(geometry_points,data.external_node_tags[row],nothing)
            point===nothing && continue
            point===(coords[1,row],coords[2,row],coords[3,row]) || return true
        end
    end
    return false
end

function _model_transfinite_surface_production(m,t,kernel,mesh,input,high_water,
                                              param_sizes,caller;
                                              side_keys,size_field=nothing,
                                              _geometry_points=nothing,
                                              birth_rows=Int32.(1:nnodes(mesh)),producer_parameters=nothing)
    limit=high_water===nothing ?
        maximum(input.mesh.entity_data.external_node_tags;init=UInt64(0)) : UInt64(high_water)
    maximum(input.mesh.entity_data.external_node_tags;init=UInt64(0))<=limit ||
        throw(ArgumentError("$caller: transfinite high water omits actual source nodes"))
    tags=zeros(UInt64,nnodes(mesh));owners=fill((2,Int32(t)),nnodes(mesh))
    parameters=Vector{Union{Nothing,Vector{Float64}}}(undef,nnodes(mesh))
    fill!(parameters,nothing);retained=falses(nnodes(mesh))
    for signed in side_keys
        curve=abs(signed)
        source=_model_mesh_curve_births(input,m,curve,caller)
        source===nothing && throw(ArgumentError(
            "$caller: transfinite Curve[$curve] has no actual Line chain"))
        rows=Int32[]
        for (nodes,physical) in _model_transfinite_boundary_events(kernel)
            for column in axes(nodes,2)
                physical[column]==curve || continue
                a,b=nodes[:,column]
                isempty(rows) ? push!(rows,a) : rows[end]==a || throw(ArgumentError(
                    "$caller: transfinite Curve[$curve] side events are not consecutive"))
                push!(rows,b)
            end
        end
        source_rows=signed>0 ? source.rows : reverse(source.rows)
        length(rows)==length(source_rows) || throw(ArgumentError(
            "$caller: transfinite side births disagree with actual Curve[$curve]"))
        data=source.mesh.entity_data
        for (row,source_row) in zip(rows,source_rows)
            tag=data.external_node_tags[source_row]
            tags[row] in (0,tag) || throw(ArgumentError(
                "$caller: transfinite junction has conflicting actual node identities"))
            all(axis->isequal(mesh.coords[axis,row],source.mesh.coords[axis,source_row]),1:3) || throw(ArgumentError(
                "$caller: transfinite source node $tag changed its actual coordinates"))
            tags[row]=tag;owners[row]=data.node_entities[source_row]
            parameters[row]=data.node_parametric[source_row]
            retained[row]=true
        end
    end
    length(birth_rows)==nnodes(mesh) && allunique(birth_rows) &&
        all(row->1<=row<=nnodes(mesh),birth_rows) || throw(ArgumentError(
            "$caller: transfinite node birth rows do not cover its actual grid"))
    fresh=Int[]
    for row in birth_rows
        retained[row] && continue
        limit<typemax(UInt64) || throw(ArgumentError(
            "$caller: transfinite node allocation exhausts UInt64 public labels"))
        limit+=1;tags[row]=limit
        if producer_parameters===nothing
            push!(fresh,row)
        else
            parameters[row]=copy(get(producer_parameters,row) do
                throw(ArgumentError("$caller: emitted Surface node lacks its actual UV binding"))
            end)
        end
    end
    # One batched parametrization: its per-call setup (surface frame, model
    # bounding box) is independent of the queried points.
    if !isempty(fresh)
        uv=model_parametrization(m,2,t,vec(mesh.coords[:,fresh]))
        length(uv)==2*length(fresh) || throw(ArgumentError(
            "$caller: transfinite Surface parametrization lost a node"))
        for (k,row) in enumerate(fresh)
            parameters[row]=uv[2k-1:2k]
        end
    end
    allunique(tags) || throw(ArgumentError(
        "$caller: transfinite emitted rows repeat an actual node identity"))
    product=_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,limit,
        UInt64[],UInt64(0),copy(birth_rows))
    _geometry_points===nothing && return product
    isempty(_geometry_points) && return product
    # The geometry layout differs from the actual one only through the birth
    # positions of boundary nodes (`_transfinite_surface_sides`). When every
    # boundary node is born where it now is, meshing the face a second time
    # gives the same grid: copy it instead.
    geometry=_model_transfinite_geometry_moved(m,t,input,_geometry_points,caller) ?
        _transfinite_surface_mesh(m,t,param_sizes,caller;size_field,
            _mesh_input=input,_geometry_points,_return_surface_production=false) :
        deepcopy(mesh)
    nnodes(geometry)==nnodes(mesh) && _surface_dim2_blocks(geometry)==_surface_dim2_blocks(mesh) ||
        throw(ArgumentError("$caller: transfinite geometry does not share actual emitter events"))
    return _ModelSurfaceGeometryProduction(product,geometry)
end
