using ..SizeField: AbstractAnisoField, metric_edge_length
using ..PipelineSupport: pipeline_float, pipeline_limit, PIPELINE_DEFAULT_MAX_NODES
using ..Mesh3D: _EHeap, _hpush!, _hpop!

@inline _metric_surface_edges(cell)=((cell[1],cell[2]),(cell[2],cell[3]),(cell[3],cell[1]))

function _metric_surface_incidence(cells,caller)
    incidence=Dict{NTuple{2,Int32},Vector{Int}}()
    for (column,cell) in enumerate(cells),edge in _metric_surface_edges(cell)
        push!(get!(()->Int[],incidence,minmax(edge...)),column)
    end
    all(hits->length(hits) in (1,2),values(incidence)) || throw(ArgumentError(
        "$caller: initial Surface has nonmanifold edge incidence"))
    return incidence
end

function _metric_surface_birth(high,count,caller,label)
    UInt64(count)<=typemax(UInt64)-high || throw(ArgumentError(
        "$caller: $label births exhaust UInt64 labels"))
    return high+UInt64(count)
end

function _metric_surface_admit(m,tag,caller)
    haskey(m.surfaces,tag) || throw(ArgumentError("$caller: unknown Surface[$tag]"))
    _surface_type(m,tag)===:plane || throw(ArgumentError("$caller: requires a planar Surface"))
    haskey(m.meshing.transfinite_surfaces,tag) && throw(ArgumentError(
        "$caller: transfinite Surface construction is not supported"))
    _extrude_gate(m,2,tag)===nothing || throw(ArgumentError(
        "$caller: extruded Surface construction is not supported"))
    _model_surface_recombined(m,tag) && throw(ArgumentError(
        "$caller: recombination is not supported"))
    isempty(get(m.embeds,(2,tag),NTuple{2,Int}[])) || throw(ArgumentError(
        "$caller: embedded Surface constraints are not supported"))
    get(m.meshing.smoothing,(2,tag),0)==0 || throw(ArgumentError(
        "$caller: Surface smoothing is not supported"))
    curves=_gmsh_native_surface_curves(m,tag)
    all(curve->_curve_type(m,curve)===:line,curves) || throw(ArgumentError(
        "$caller: requires straight boundary Curves"))
    curve_set=Set(curves)
    any(constraint->constraint.dim==2 &&
        (abs(Int(constraint.slave_entity))==tag || abs(Int(constraint.master_entity))==tag) ||
        constraint.dim==1 && (abs(Int(constraint.slave_entity)) in curve_set ||
        abs(Int(constraint.master_entity)) in curve_set),model_periodic_constraints(m)) &&
        throw(ArgumentError("$caller: periodic Surface or boundary Curve construction is not supported"))
    _model_surface_plane(m,tag,caller)
    return nothing
end

function _metric_surface_refine(m,tag,field,target,node_limit,triangle_limit,
                                input,options,node_high_water,element_high_water,caller)
    authority,limit=_gmsh_native_surface_lower_input(m,options,input,node_high_water,caller)
    frame=_gmsh_native_surface_frame(m,tag,options,caller)
    source=_gmsh_native_surface_source(m,tag,authority,frame,caller;ctx_lc=options.ctx_lc)
    length(source.tags)<=node_limit || throw(ArgumentError(
        "$caller: retained lower vertices exceed max_nodes"))
    # Without interior constraints, Euler's disk-with-holes relation fixes
    # the initial triangle count before its triangulation is allocated.
    initial_count=length(source.tags)+2(length(m.surfaces[tag])-1)-2
    initial_count<=triangle_limit || throw(ArgumentError(
        "$caller: initial Surface exceeds max_triangles"))
    initial=_gmsh_native_initial_surface(m,tag,options,authority,limit,caller;
        algorithm=3,element_high_water)
    mesh=initial.mesh
    nnodes(mesh)<=node_limit && ntris(mesh)<=triangle_limit || throw(ArgumentError(
        "$caller: initial lower mesh exceeds max_nodes or max_triangles"))
    validate(mesh).ok || throw(ArgumentError("$caller: initial Surface is invalid"))
    coordinates=NTuple{3,Float64}[Tuple(mesh.coords[:,row]) for row in axes(mesh.coords,2)]
    parameters=NTuple{2,Float64}[_gmsh_native_surface_parameters(frame,point) for point in coordinates]
    cells=NTuple{3,Int32}[Tuple(cell) for cell in eachcol(mesh.tris)]
    incidence=_metric_surface_incidence(cells,caller)
    boundary=Set(edge for (edge,hits) in incidence if length(hits)==1)
    original_boundary=Set(minmax(initial.node_tags[edge[1]],initial.node_tags[edge[2]]) for edge in boundary)
    wanted_boundary=Set(minmax(source.tags[edge[1]],source.tags[edge[2]]) for edge in source.boundaries)
    original_boundary==wanted_boundary || throw(ArgumentError(
        "$caller: initial Surface loses its actual lower Curve-cell incidence"))
    score(edge)=metric_edge_length(field,coordinates[edge[1]],coordinates[edge[2]];entity=(2,tag))
    all(edge->score(edge)<=target,boundary) || throw(ArgumentError(
        "$caller: retained boundary exceeds max_metric_edge; provide finer lower Curve cells"))
    labels=copy(initial.element_tags);numbers=copy(initial.node_tags)
    owners=copy(initial.node_owners);node_parameters=deepcopy(initial.node_parameters)
    retained=copy(initial.retained_nodes);birth_rows=copy(initial.birth_rows)
    alive=trues(length(cells));live_count=length(cells)
    node_max=initial.node_high_water;cell_max=initial.element_high_water
    heap=_EHeap();queued=Set{NTuple{2,Int32}}()
    function enqueue(edge)
        haskey(incidence,edge) && !(edge in queued) || return
        violation=score(edge)
        isfinite(violation) && violation>=0 || throw(ArgumentError(
            "$caller: edge metric must be finite and nonnegative"))
        if violation>target
            edge in boundary && throw(ArgumentError("$caller: retained boundary metric changed"))
            push!(queued,edge);_hpush!(heap,(violation,edge[1],edge[2]))
        end
    end
    for edge in sort!(collect(keys(incidence)));enqueue(edge);end
    locations=Set(coordinates)
    while !isempty(heap.d)
        _,a,b=_hpop!(heap);edge=(a,b);delete!(queued,edge)
        haskey(incidence,edge) || continue
        score(edge)>target || continue
        hits=copy(incidence[edge]);length(hits)==2 || throw(ArgumentError(
            "$caller: metric refinement would split an owned lower boundary"))
        length(coordinates)<node_limit && live_count+length(hits)<=triangle_limit ||
            throw(ArgumentError("$caller: metric refinement exceeds max_nodes or max_triangles"))
        next_node=_metric_surface_birth(node_max,1,caller,"node")
        next_cell=_metric_surface_birth(cell_max,2length(hits),caller,"cell")
        midpoint=ntuple(axis->_convex_coordinate(coordinates[a][axis],coordinates[b][axis],.5),3)
        midpoint in locations && throw(ArgumentError(
            "$caller: edge midpoint is not a distinct representable Surface vertex"))
        uv=ntuple(axis->_convex_coordinate(parameters[a][axis],parameters[b][axis],.5),2)
        row=Int32(length(coordinates)+1)
        children=NTuple{3,Int32}[]
        for column in hits
            cell=cells[column]
            slot=only(index for index in 1:3 if minmax(cell[index],cell[mod1(index+1,3)])==edge)
            first,last,other=cell[slot],cell[mod1(slot+1,3)],cell[mod1(slot+2,3)]
            winding=orient2(parameters[first],parameters[last],parameters[other])
            first_child=orient2(parameters[first],uv,parameters[other])
            second_child=orient2(uv,parameters[last],parameters[other])
            isfinite(winding) && isfinite(first_child) && isfinite(second_child) &&
                winding!=0 && sign(first_child)==sign(winding) && sign(second_child)==sign(winding) ||
                throw(ArgumentError("$caller: midpoint creates a flat or inverted Surface triangle"))
            child_areas=(triangle_area(coordinates[first],midpoint,coordinates[other]),
                         triangle_area(midpoint,coordinates[last],coordinates[other]))
            all(area->isfinite(area) && area>0,child_areas) || throw(ArgumentError(
                "$caller: midpoint creates an unrepresentable physical Surface triangle"))
            push!(children,(first,row,other),(row,last,other))
        end
        push!(coordinates,midpoint);push!(parameters,uv);push!(locations,midpoint)
        push!(numbers,next_node);push!(owners,(2,Int32(tag)));push!(node_parameters,collect(uv))
        push!(retained,false);push!(birth_rows,row);node_max=next_node
        for column in hits
            alive[column]=false
            for old_edge in _metric_surface_edges(cells[column])
                key=minmax(old_edge...);filter!(!=(column),incidence[key])
                isempty(incidence[key]) && delete!(incidence,key)
            end
        end
        new_edges=Set{NTuple{2,Int32}}()
        for child in children
            push!(cells,child);push!(alive,true);cell_max+=UInt64(1);push!(labels,cell_max)
            for child_edge in _metric_surface_edges(child)
                key=minmax(child_edge...);push!(get!(()->Int[],incidence,key),length(cells));push!(new_edges,key)
            end
        end
        cell_max==next_cell || error("$caller: internal cell birth count mismatch")
        live_count+=length(hits)
        for new_edge in sort!(collect(new_edges));enqueue(new_edge);end
    end
    live=findall(alive);connectivity=Matrix{Int32}(undef,3,length(live))
    for (column,index) in enumerate(live);connectivity[:,column].=cells[index];end
    output_coordinates=Matrix{Float64}(undef,3,length(coordinates))
    for (row,point) in enumerate(coordinates);output_coordinates[:,row].=point;end
    output=Mesh(output_coordinates;tris=connectivity)
    validate(output).ok || throw(ErrorException("$caller: refined Surface is invalid"))
    final_incidence=_metric_surface_incidence(NTuple{3,Int32}[Tuple(cell) for cell in eachcol(output.tris)],caller)
    Set(edge for (edge,hits) in final_incidence if length(hits)==1)==boundary || throw(ErrorException(
        "$caller: refined Surface changes its retained boundary"))
    all(edge->score(edge)<=target,keys(final_incidence)) || throw(ErrorException(
        "$caller: final edge exceeds max_metric_edge"))
    area(mesh)=sum(triangle_area((Tuple(mesh.coords[:,row]) for row in cell)...) for cell in eachcol(mesh.tris))
    initial_area=area(mesh);final_area=area(output)
    isfinite(initial_area) && initial_area>0 && isfinite(final_area) &&
        isapprox(final_area,initial_area;rtol=256eps(Float64)*max(1,ntris(output)),atol=0) ||
        throw(ErrorException("$caller: refinement changes the Surface area"))
    isequal(output.coords[:,1:nnodes(mesh)],mesh.coords) && numbers[1:nnodes(mesh)]==initial.node_tags ||
        throw(ErrorException("$caller: refinement changes a retained lower vertex"))
    return _ModelSurfaceProduction(output,numbers,owners,node_parameters,retained,node_max,
        labels[live],cell_max,birth_rows)
end

"""
    mesh_model_surface(model, tag, field::AbstractAnisoField;
                       max_metric_edge=1.0, max_nodes=10_000_000,
                       max_triangles=20_000_000) -> Mesh

Generate an ordinary planar triangular Surface using a spatially varying tensor
field. Conforming interior edge bisection measures the full world-coordinate
metric, equivalently its tangential pullback, at both endpoints and the midpoint
of every edge. The output satisfies `metric_edge_length(field,a,b)<=max_metric_edge`.

Actual lower Point/Line products retain their UID, XYZ and boundary incidence.
They must already satisfy the requested bound: coarse boundary edges are rejected
instead of replacing owned Curve cells. Straight planar boundaries, including
holes, are supported. Transfinite/extruded/periodic/recombined Surfaces, embeddings
and Surface smoothing are rejected. New vertices and cells have owned constructor
births; resource and representability failures return no partial mesh or model
mutation. The supplied field, including any explicit wrapper, is used as given.

This explicit tensor overload leaves the existing native scalar generation path
unchanged. It does not claim BAMG ordering, native anisotropic API generation, a
minimum triangle angle, or pointwise bounds between the three metric samples.
"""
function mesh_model_surface(m::GeoModel,tag::Integer,field::AbstractAnisoField;
        max_metric_edge=1.0,max_nodes=PIPELINE_DEFAULT_MAX_NODES,max_triangles=20_000_000,
        _mesh_input=nothing,_mesh1d_options=nothing,_mesh_node_high_water=nothing,
        _mesh_element_high_water=UInt64(0),_return_surface_production=false)
    caller="mesh_model_surface tensor metric"
    target=pipeline_float(max_metric_edge,caller,"max_metric_edge")
    target>0 || throw(ArgumentError("$caller: max_metric_edge must be positive"))
    node_limit=pipeline_limit(max_nodes,caller,"max_nodes",typemax(Int32))
    triangle_limit=pipeline_limit(max_triangles,caller,"max_triangles",typemax(Int32)÷2)
    t=_tag(tag,caller,2);_metric_surface_admit(m,t,caller)
    working=_model_mesh_transaction_copy(m)
    options=_model_algorithm_size_options(working,_mesh1d_options,caller)
    product=_metric_surface_refine(working,t,field,target,node_limit,triangle_limit,
        _mesh_input,options,_mesh_node_high_water,_mesh_element_high_water,caller)
    return _return_surface_production ? product : product.mesh
end
