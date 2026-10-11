# Native recombination stages, independently expressed from their numerical
# and topological contracts. The SDK controls are kept outside this module;
# no Gmsh or restricted matching implementation is embedded here.
struct _PrimarySurfaceContext{F,P}
    entity::Int32
    node_entities::Vector{Tuple{Int,Int32}}
    node_tags::Vector{UInt64}
    parameters::Matrix{Float64}
    evaluate::F
    parameter_present::BitVector
    project::P
    reposition::Bool
    topology_passes::Int
    smoothing::Int
    topology_admitted::Bool
end

function _PrimarySurfaceContext(entity,node_entities,node_tags,parameters,evaluate;
        reposition=true,topology_passes=5,smoothing=1,topology_admitted=true,
        parameter_present=trues(length(node_tags)),project=nothing)
    length(node_entities)==length(node_tags)==size(parameters,2) ||
        throw(ArgumentError("native Surface context rows do not align"))
    size(parameters,1)==2 || throw(ArgumentError("native Surface context requires two parameter rows"))
    length(parameter_present)==length(node_tags) && all(value->value isa Bool,parameter_present) ||
        throw(ArgumentError("native Surface context parameter presence does not bind input rows"))
    all(owner->owner[1] in 0:3 && owner[2]>0,node_entities) ||
        throw(ArgumentError("native Surface context requires classified positive entity identities"))
    all(tag->tag isa Integer && !(tag isa Bool) && 0<tag<=typemax(UInt64),node_tags) &&
        allunique(node_tags) || throw(ArgumentError("native Surface context requires unique positive public node identities"))
    topology_passes isa Integer && !(topology_passes isa Bool) &&
        smoothing isa Integer && !(smoothing isa Bool) ||
        throw(ArgumentError("native optimization pass counts must be integers"))
    return _PrimarySurfaceContext(Int32(entity),Tuple{Int,Int32}[(Int(d),Int32(t)) for (d,t) in node_entities],
        UInt64.(node_tags),copy(parameters),evaluate,BitVector(parameter_present),project,Bool(reposition),Int(topology_passes),
        Int(smoothing),Bool(topology_admitted))
end

# This is the affine Plane operation order used by native public evaluation
# and by relocation's successive Float64 additions.
struct _PrimaryPlaneEvaluation
    origin::NTuple{3,Float64}
    first::NTuple{3,Float64}
    second::NTuple{3,Float64}
end
@inline (plane::_PrimaryPlaneEvaluation)(u,v)=ntuple(axis->
    plane.origin[axis]+plane.first[axis]*u+plane.second[axis]*v,3)
@inline function _primary_plane_parameters(plane::_PrimaryPlaneEvaluation,point)
    delta=ntuple(axis->point[axis]-plane.origin[axis],3)
    return (_dot3(delta,plane.first),_dot3(delta,plane.second))
end
struct _PrimaryBlockedSurfaceEvaluation
    message::String
end
(evaluation::_PrimaryBlockedSurfaceEvaluation)(u,v)=throw(ArgumentError(evaluation.message))

function _primary_relocation_fraction(quality_at)
    target=quality_at(0.0)
    target>0.5 && return 0.0
    root=sqrt(5.0);ratio=(root-1.0)*0.5;complement=(3.0-root)*0.5
    lower=0.0;upper=1.0
    left=upper-ratio*(upper-lower);right=lower+ratio*(upper-lower)
    left_score=quality_at(left);right_score=quality_at(right)
    while abs(upper-lower)>=1e-6
        if left_score<right_score
            lower=left
            abs(upper-lower)<1e-6 && break
            left=right;left_score=right_score
            right=upper-complement*(upper-lower);right_score=quality_at(right)
        else
            upper=right
            abs(upper-lower)<1e-6 && break
            right=left;right_score=left_score
            left=lower+complement*(upper-lower);left_score=quality_at(left)
        end
    end
    return quality_at(lower)<target ? 0.0 : lower
end

function _primary_relocate_vertex!(coords,parameters,triangles,quadrangles,vertex,evaluate,
                                   triangle_columns,quad_columns)
    total_u=0.0;total_v=0.0;count=0
    for column in triangle_columns,node in triangles[column]
        total_u+=parameters[1,node];total_v+=parameters[2,node];count+=1
    end
    for column in quad_columns,node in quadrangles[column]
        total_u+=parameters[1,node];total_v+=parameters[2,node];count+=1
    end
    count==0 && return nothing
    reciprocal=1.0/count
    target=(total_u*reciprocal,total_v*reciprocal)
    before=(parameters[1,vertex],parameters[2,vertex])
    all(isfinite,target) && all(isfinite,before) || throw(ArgumentError(
        "native relocation requires stored finite Surface parameters"))
    original=ntuple(axis->coords[axis,vertex],3)
    function quality_at(fraction)
        trial_u=target[1]*(1.0-fraction)+before[1]*fraction
        trial_v=target[2]*(1.0-fraction)+before[2]*fraction
        trial_point=evaluate(trial_u,trial_v)
        all(isfinite,trial_point) || throw(ArgumentError("native Surface evaluation is nonfinite"))
        for axis in 1:3;coords[axis,vertex]=trial_point[axis];end
        score=1.0
        try
            for column in triangle_columns
                score=_primary_min(abs(_primary_triangle_gamma(coords,triangles[column])),score)
            end
            for column in quad_columns
                score=_primary_min(_gmsh_recombine_pair_measure(coords,quadrangles[column]),score)
            end
            return score
        finally
            for axis in 1:3;coords[axis,vertex]=original[axis];end
        end
    end
    fraction=_primary_relocation_fraction(quality_at)
    u=target[1]*(1.0-fraction)+before[1]*fraction
    v=target[2]*(1.0-fraction)+before[2]*fraction
    point=evaluate(u,v)
    all(isfinite,point) || throw(ArgumentError("native Surface evaluation is nonfinite"))
    for axis in 1:3;coords[axis,vertex]=point[axis];end
    parameters[1,vertex]=u;parameters[2,vertex]=v
    return nothing
end

function _primary_relocate_world_vertex!(coords,triangles,quadrangles,vertex,project,
                                        triangle_columns,quad_columns)
    project===nothing && throw(ArgumentError("native relocation requires a supported closest-point context"))
    total_x=0.0;total_y=0.0;total_z=0.0;count=0
    for column in triangle_columns,node in triangles[column]
        total_x+=coords[1,node];total_y+=coords[2,node];total_z+=coords[3,node];count+=1
    end
    for column in quad_columns,node in quadrangles[column]
        total_x+=coords[1,node];total_y+=coords[2,node];total_z+=coords[3,node];count+=1
    end
    count==0 && return nothing
    reciprocal=1.0/count;target=(total_x*reciprocal,total_y*reciprocal,total_z*reciprocal)
    original=ntuple(axis->coords[axis,vertex],3)
    function point_at(fraction)
        query=ntuple(axis->target[axis]*(1.0-fraction)+original[axis]*fraction,3)
        projected=project(query)
        all(isfinite,projected) || throw(ArgumentError("native Surface closest point is nonfinite"))
        return projected
    end
    function quality_at(fraction)
        trial_point=point_at(fraction)
        for axis in 1:3;coords[axis,vertex]=trial_point[axis];end
        score=1.0
        try
            for column in triangle_columns
                score=_primary_min(abs(_primary_triangle_gamma(coords,triangles[column])),score)
            end
            for column in quad_columns
                score=_primary_min(_gmsh_recombine_pair_measure(coords,quadrangles[column]),score)
            end
            return score
        finally
            for axis in 1:3;coords[axis,vertex]=original[axis];end
        end
    end
    point=point_at(_primary_relocation_fraction(quality_at))
    for axis in 1:3;coords[axis,vertex]=point[axis];end
    # An unparameterized native vertex retains empty raw parameter storage.
    return nothing
end

function _primary_relocate!(coords,triangles,quadrangles,context::_PrimarySurfaceContext)
    (!context.reposition || context.evaluate===nothing || context.smoothing<=0) && return nothing
    triangles_at=[Int[] for _ in axes(coords,2)]
    quads_at=[Int[] for _ in axes(coords,2)]
    for (column,cell) in enumerate(triangles),node in cell;push!(triangles_at[node],column);end
    for (column,cell) in enumerate(quadrangles),node in cell;push!(quads_at[node],column);end
    order=sortperm(context.node_tags)
    for _ in 1:context.smoothing,vertex in order
        context.node_entities[vertex][1]==2 || continue
        isempty(triangles_at[vertex]) && isempty(quads_at[vertex]) && continue
        if context.parameter_present[vertex]
            _primary_relocate_vertex!(coords,context.parameters,triangles,quadrangles,vertex,
                context.evaluate,triangles_at[vertex],quads_at[vertex])
        else
            _primary_relocate_world_vertex!(coords,triangles,quadrangles,vertex,
                context.project,triangles_at[vertex],quads_at[vertex])
            # A subsequently moved parameterized neighbor reparametrizes this
            # vertex from its current XYZ; this local cache is never raw storage.
            uv=_primary_plane_parameters(context.evaluate,ntuple(axis->coords[axis,vertex],3))
            context.parameters[1,vertex]=uv[1];context.parameters[2,vertex]=uv[2]
        end
    end
    return nothing
end

# One complete degree-two pass. Repeating passes is controlled by the native
# outer stage schedule, including its final empty pass and relocation.
function _primary_collapse_degree_two(coords,triangles,quadrangles,tags,births,
                                     context::_PrimarySurfaceContext,counter)
    adjacency=[Tuple{Bool,Int}[] for _ in axes(coords,2)]
    for (column,cell) in enumerate(triangles),node in cell;push!(adjacency[node],(false,column));end
    for (column,cell) in enumerate(quadrangles),node in cell;push!(adjacency[node],(true,column));end
    touched=falses(length(quadrangles));removed=falses(size(coords,2))
    added=NTuple{4,Int32}[];added_tags=Int32[];added_births=UInt64[]
    for vertex in sortperm(context.node_tags)
        context.node_entities[vertex]==(2,context.entity) || continue
        incident=adjacency[vertex]
        length(incident)==2 && all(first,incident) || continue
        first_quad=incident[1][2];second_quad=incident[2][2]
        (touched[first_quad] || touched[second_quad]) && continue
        first_cell=quadrangles[first_quad];second_cell=quadrangles[second_quad]
        position=findfirst(==(Int32(vertex)),first_cell)
        a=first_cell[mod1(position+1,4)];b=first_cell[mod1(position+2,4)]
        c=first_cell[mod1(position+3,4)]
        remaining=Int32[node for node in second_cell if node!=vertex && node!=a && node!=c]
        length(remaining)==1 && a in second_cell && c in second_cell || throw(ArgumentError(
            "native degree-two topology has unsupported incidence"))
        tags[first_quad]==tags[second_quad] || throw(ArgumentError(
            "native degree-two topology crosses a physical-tag boundary"))
        merged=(a,b,c,only(remaining))
        allunique(merged) || throw(ArgumentError("native degree-two topology repeats a vertex"))
        counter=Base.checked_add(counter,UInt64(1))
        push!(added,merged);push!(added_tags,tags[first_quad]);push!(added_births,counter)
        touched[first_quad]=true;touched[second_quad]=true;removed[vertex]=true
    end
    return (quadrangles=append!(quadrangles[.!touched],added),
        tags=append!(tags[.!touched],added_tags),births=append!(births[.!touched],added_births),
        removed,counter,changes=length(added))
end

# Only a qualifying diamond requires the still-unimplemented collapse
# policy. An ordinary mesh with no such configuration proceeds normally.
function _primary_require_no_diamonds(triangles,quadrangles,context::_PrimarySurfaceContext)
    context.evaluate===nothing && return nothing
    touched=Set{Int32}(node for cell in triangles for node in cell)
    degrees=zeros(Int,length(context.node_tags))
    for cell in quadrangles,node in cell;degrees[node]+=1;end
    for cell in quadrangles
        any(node->node in touched,cell) && continue
        all(node->context.node_entities[node]==(2,context.entity),cell) || continue
        ((degrees[cell[1]]==degrees[cell[3]]==3) ||
         (degrees[cell[2]]==degrees[cell[4]]==3)) || continue
        throw(ArgumentError("native recombination diamond optimization is not implemented for this configuration"))
    end
    return nothing
end

function _primary_optimized_emission(coords,source_triangles,source_tags,
                                     emission::_PrimaryWeightedEmission,context::_PrimarySurfaceContext,
                                     minimum_quality::Float64)
    size(coords,2)==length(context.node_tags) || throw(ArgumentError("native Surface context does not bind input rows"))
    triangles=NTuple{3,Int32}[];triangle_tags=Int32[];triangle_births=UInt64[]
    for column in axes(source_triangles,2)
        emission.used[column] && continue
        push!(triangles,ntuple(row->source_triangles[row,column],3))
        push!(triangle_tags,source_tags[column]);push!(triangle_births,UInt64(0))
    end
    append!(triangles,emission.triangles);append!(triangle_tags,emission.triangle_tags)
    append!(triangle_births,emission.triangle_births)
    quadrangles=copy(emission.quadrangles);quad_tags=copy(emission.quadrangle_tags)
    quad_births=copy(emission.quadrangle_births);counter=emission.factory_births
    removed=falses(size(coords,2))
    _primary_relocate!(coords,triangles,quadrangles,context)
    if context.topology_passes>0 && context.topology_admitted
        for _ in 0:context.topology_passes
            changes=0
            while true
                pass=_primary_collapse_degree_two(coords,triangles,quadrangles,quad_tags,quad_births,context,counter)
                quadrangles=pass.quadrangles;quad_tags=pass.tags;quad_births=pass.births
                counter=pass.counter;removed .|= pass.removed;changes+=pass.changes
                pass.changes==0 && break
            end
            _primary_require_no_diamonds(triangles,quadrangles,context)
            _primary_relocate!(coords,triangles,quadrangles,context)
            changes==0 && break
        end
    end
    kept=NTuple{4,Int32}[];kept_tags=Int32[];kept_births=UInt64[]
    for index in eachindex(quadrangles)
        nodes=quadrangles[index]
        split=_primary_postquality_split(coords,nodes,minimum_quality)
        if split===nothing
            push!(kept,nodes);push!(kept_tags,quad_tags[index]);push!(kept_births,quad_births[index])
        else
            _primary_supported_quad_split(coords,split.triangles) || throw(ArgumentError(
                "native postquality materialization has no supported consistent projection"))
            for (cell,offset) in zip(split.triangles,split.birth_offsets)
                push!(triangles,cell);push!(triangle_tags,quad_tags[index])
                push!(triangle_births,Base.checked_add(counter,UInt64(offset)))
            end
            counter=Base.checked_add(counter,UInt64(split.reserved))
        end
    end
    _primary_relocate!(coords,triangles,kept,context)
    for nodes in kept
        _strict_convex_projection(coords,nodes) && last(_quad_quality_range_result(coords,nodes)) ||
            throw(ArgumentError("native optimized quadrangle is outside the supported convex geometry certificate"))
    end
    result=_PrimaryWeightedEmission(trues(size(source_triangles,2)),triangles,triangle_tags,triangle_births,
        kept,kept_tags,kept_births,counter,emission.pair_count,emission.pair_passes)
    return (emission=result,removed_nodes=removed,node_parameters=context.parameters)
end
