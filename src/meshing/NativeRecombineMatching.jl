# The graph interpretation below is independently written from the pinned
# Gmsh 4.15.2 public behavior and meshGFaceOptimize.cpp graph specification.
# It does not adapt the restricted Concorde implementation. Integer matching
# uses the separately documented BSD kernel and independent ordered
# singleton/half-integral policy. General graphs retaining half-cycles still
# have an explicit primary tie-parity gap.

struct _PrimaryWeightedGraph
    edges::Vector{Tuple{Int,Int,Int64}}
    physical::Vector{_QuadCandidate}
end

function _primary_weighted_graph(coords,triangles;
                                 protected=nothing,embedded=nothing,
                                 node_owner_dims=nothing)
    size(coords,1)==3 && size(triangles,1)==3 || throw(ArgumentError(
        "recombine_triangles: weighted graph requires three-coordinate triangles"))
    count=size(triangles,2)
    count<=typemax(Int32) || throw(ArgumentError(
        "recombine_triangles: weighted triangle identities exceed Int32"))
    if node_owner_dims!==nothing
        length(node_owner_dims)==size(coords,2) || throw(ArgumentError(
            "recombine_triangles: node classification count does not match coordinates"))
        all(dim->dim isa Integer && !(dim isa Bool) && 0<=dim<=3,node_owner_dims) ||
            throw(ArgumentError("recombine_triangles: node classifications must lie in 0:3"))
    end
    # Callers supply columns in increasing public-tag order. The first
    # triangle's directed edge is retained independently of that key order.
    owners=Dict{NTuple{2,Int32},Tuple{Int,Int}}()
    for column in 1:count
        for (i,j) in ((1,2),(2,3),(3,1))
            first_node=triangles[i,column];second_node=triangles[j,column]
            1<=first_node<=size(coords,2) && 1<=second_node<=size(coords,2) ||
                throw(ArgumentError("recombine_triangles: weighted graph node is out of range"))
            edge=_edge_key(first_node,second_node)
            first,second=get(owners,edge,(0,0))
            second==0 || throw(ArgumentError(
                "recombine_triangles: non-manifold edge $edge has more than two incident triangles"))
            owners[edge]=first==0 ? (column,0) : (first,column)
        end
    end
    physical=_QuadCandidate[]
    periodic=Dict{Int32,Tuple{Int,Int}}()
    for edge in sort!(collect(keys(owners)))
        first,second=owners[edge]
        if second!=0
            ((protected!==nothing && edge in protected) ||
             (embedded!==nothing && edge in embedded)) && continue
            triangle1=(triangles[1,first],triangles[2,first],triangles[3,first])
            triangle2=(triangles[1,second],triangles[2,second],triangles[3,second])
            a,b=_gmsh_first_candidate_edge(triangle1,edge)
            c=_third_vertex(triangle1,a,b);d=_third_vertex(triangle2,a,b)
            # This is the raw constructor order used to evaluate eta. It is
            # not a declaration that this physical quadrangle is admissible.
            # Geometry is certified only AFTER solving the full objective.
            push!(physical,_QuadCandidate(Int32(first),Int32(second),(a,c,b,d),0.0,edge))
        elseif embedded===nothing || !(edge in embedded)
            for vertex in edge
                previous=get(periodic,vertex,(0,0))
                if previous[1]==0
                    periodic[vertex]=(first,0)
                elseif previous[1]!=first
                    periodic[vertex]=(previous[1],first)
                else
                    delete!(periodic,vertex)
                end
            end
        end
    end
    _gmsh_sort_candidates!(physical,coords,triangles)
    edges=Tuple{Int,Int,Int64}[]
    sizehint!(edges,length(physical)+length(periodic))
    penalty_cut=floor(Int64,1000.0*_gm87_exp(0.1))
    for pair in physical
        eta=_gmsh_recombine_pair_measure(coords,pair.nodes)
        raw=floor(Int64,1000.0*_gm87_exp(-eta))
        lower_owned=node_owner_dims===nothing ? 0 :
            sum(node_owner_dims[Int(node)]<2 for node in pair.nodes)
        cost=raw>penalty_cut && lower_owned>2 ? Int64(5000) :
            raw>=1000 && lower_owned>2 ? Int64(10000) : raw
        push!(edges,(Int(pair.first_triangle),Int(pair.second_triangle),cost))
    end
    for vertex in sort!(collect(keys(periodic)))
        first,second=periodic[vertex]
        # Preserve literal directed/parallel edges and the native map's
        # null value initialization. Solvers validate and ignore self loops.
        push!(edges,(first,second==0 ? 1 : second,Int64(100000)))
    end
    return _PrimaryWeightedGraph(edges,physical)
end

function _primary_weighted_matching(coords,triangles;kwargs...)
    graph=_primary_weighted_graph(coords,triangles;kwargs...)
    matching=ordered_perfect_matching(size(triangles,2),graph.edges)
    return (graph=graph,matching=matching)
end

@inline function _primary_unit3(vector)
    length=sqrt(_dot3(vector,vector))
    length==0.0 && return vector
    reciprocal=1.0/length
    return (vector[1]*reciprocal,vector[2]*reciprocal,vector[3]*reciprocal)
end

@inline function _primary_cross3(a,b)
    return (a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
end

# This compatibility gamma keeps the native raw norme/product operation
# order. It is separate from the scale-robust public quality report.
function _primary_triangle_gamma(coords,nodes::NTuple{3,Int32})
    points=ntuple(i->ntuple(d->coords[d,Int(nodes[i])],3),3)
    a=_primary_unit3(_sub3(points[3],points[2]))
    b=_primary_unit3(_sub3(points[1],points[3]))
    c=_primary_unit3(_sub3(points[2],points[1]))
    first=_primary_cross3(b,c);second=_primary_cross3(c,a);third=_primary_cross3(a,b)
    sina=sqrt(_dot3(first,first));sinb=sqrt(_dot3(second,second));sinc=sqrt(_dot3(third,third))
    sina==0.0 && sinb==0.0 && sinc==0.0 && return 0.0
    return 2.0*(2.0*sina*sinb*sinc/(sina+sinb+sinc))
end

@inline function _primary_triangle_normal(coords,nodes::NTuple{3,Int32})
    points=ntuple(i->ntuple(d->coords[d,Int(nodes[i])],3),3)
    return _primary_unit3(_primary_cross3(_sub3(points[2],points[1]),_sub3(points[3],points[1])))
end

@inline _primary_min(first::Float64,second::Float64)=second<first ? second : first

# A split is chosen only AFTER raw weighted matching and quad publication.
# Native constructs all four candidate triangles, including the two that
# are discarded. Offsets identify their original factory birth order.
function _primary_postquality_split(coords,nodes::NTuple{4,Int32},minimum_quality::Float64)
    minimum_quality>0.0 || return nothing
    _gmsh_recombine_pair_measure(coords,nodes)<minimum_quality+1e-12 || return nothing
    a,b,c,d=nodes
    first=(a,b,c);second=(c,d,a);third=(b,c,d);fourth=(d,a,b)
    quality1=_primary_min(_primary_triangle_gamma(coords,first),_primary_triangle_gamma(coords,second))
    quality2=_primary_min(_primary_triangle_gamma(coords,third),_primary_triangle_gamma(coords,fourth))
    orientation2=_dot3(_primary_triangle_normal(coords,third),_primary_triangle_normal(coords,fourth))
    return quality1>quality2 || orientation2<0.0 ?
        (triangles=(first,second),birth_offsets=(1,2),reserved=4) :
        (triangles=(third,fourth),birth_offsets=(3,4),reserved=4)
end

struct _PrimaryWeightedEmission
    used::BitVector
    triangles::Vector{NTuple{3,Int32}}
    triangle_tags::Vector{Int32}
    triangle_births::Vector{UInt64}
    quadrangles::Vector{NTuple{4,Int32}}
    quadrangle_tags::Vector{Int32}
    quadrangle_births::Vector{UInt64}
    factory_births::UInt64
    pair_count::Int
    pair_passes::Int
end

function _primary_matched_cycle(triangles,pair::_QuadCandidate)
    first=Int(pair.first_triangle);second=Int(pair.second_triangle)
    a=(triangles[1,first],triangles[2,first],triangles[3,first])
    b=(triangles[1,second],triangles[2,second],triangles[3,second])
    opposite=_third_vertex(b,pair.edge...)
    start=1
    for index in 1:3
        a[index] in b || (start=index;break)
    end
    return (a[start],a[mod1(start+1,3)],opposite,a[mod1(start+2,3)])
end

@inline function _primary_greedy_cycle(pair::_QuadCandidate)
    first,opposite_first,second,opposite_second=pair.nodes
    return (first,opposite_second,second,opposite_first)
end

function _primary_supported_triangle_split(coords,triangles,pair,cells)
    first_column=Int(pair.first_triangle);second_column=Int(pair.second_triangle)
    first=(triangles[1,first_column],triangles[2,first_column],triangles[3,first_column])
    second=(triangles[1,second_column],triangles[2,second_column],triangles[3,second_column])
    for axes in ((1,2),(2,3),(3,1))
        orientation=orient2(_projected_point(coords,first[1],axes),
            _projected_point(coords,first[2],axes),_projected_point(coords,first[3],axes))
        orientation==0 && continue
        orient2(_projected_point(coords,second[1],axes),
            _projected_point(coords,second[2],axes),_projected_point(coords,second[3],axes))==orientation || continue
        reference_orientation=orientation;projection_axes=axes
        all(cell->orient2(_projected_point(coords,cell[1],projection_axes),
            _projected_point(coords,cell[2],projection_axes),_projected_point(coords,cell[3],projection_axes))==reference_orientation,cells) && return true
    end
    return false
end

# For a preexisting quadrangle the selected two triangles themselves certify
# a simple consistently oriented boundary projection. Both possible splits
# retain the same four directed boundary edges and shoelace area.
function _primary_supported_quad_split(coords,cells)
    for projection_axes in ((1,2),(2,3),(3,1))
        first=cells[1]
        orientation=orient2(_projected_point(coords,first[1],projection_axes),
            _projected_point(coords,first[2],projection_axes),
            _projected_point(coords,first[3],projection_axes))
        orientation==0 && continue
        reference_orientation=orientation
        all(cell->orient2(_projected_point(coords,cell[1],projection_axes),
            _projected_point(coords,cell[2],projection_axes),
            _projected_point(coords,cell[3],projection_axes))==reference_orientation,cells) && return true
    end
    return false
end

# New triangles and quadrangles carry their exact relative factory birth,
# including discarded temporary cells. The caller can check its global tag
# limit before atomically publishing any changed model, cache or counter.
function _primary_weighted_emission(coords,triangles,triangle_tags;
                                    protected=nothing,embedded=nothing,
                                    node_owner_dims=nothing,minimum_quality::Float64=0.0,
                                    algorithm::Symbol=:blossom,recombine_angle::Float64=1.0,
                                    existing_quadrangles=nothing,existing_quadrangle_tags=nothing,
                                    defer_postquality::Bool=false)
    isfinite(minimum_quality) || throw(ArgumentError(
        "recombine_triangles: primary minimum quality must be finite"))
    length(triangle_tags)==size(triangles,2) || throw(ArgumentError(
        "recombine_triangles: weighted triangle physical-tag count mismatch"))
    algorithm in (:greedy,:blossom) || throw(ArgumentError(
        "recombine_triangles: primary algorithm must be :greedy or :blossom"))
    isfinite(recombine_angle) && 0.0<=recombine_angle<=90.0 || throw(ArgumentError(
        "recombine_triangles: primary recombine angle must lie in 0:90"))
    if existing_quadrangles!==nothing
        size(existing_quadrangles,1)==4 || throw(ArgumentError(
            "recombine_triangles: existing quadrangles require four nodes"))
        existing_quadrangle_tags!==nothing &&
            length(existing_quadrangle_tags)==size(existing_quadrangles,2) || throw(ArgumentError(
            "recombine_triangles: existing quadrangle physical-tag count mismatch"))
        for column in axes(existing_quadrangles,2)
            nodes=ntuple(row->existing_quadrangles[row,column],4)
            allunique(nodes) && all(node->node isa Integer && !(node isa Bool) &&
                1<=node<=size(coords,2),nodes) || throw(ArgumentError(
                "recombine_triangles: existing quadrangle has invalid node identities"))
        end
    elseif existing_quadrangle_tags!==nothing
        throw(ArgumentError("recombine_triangles: existing quadrangle tags require quadrangles"))
    end
    graph=_primary_weighted_graph(coords,triangles;protected,embedded,node_owner_dims)
    even_blossom=algorithm===:blossom && iseven(size(triangles,2))
    matching=even_blossom ? ordered_perfect_matching(size(triangles,2),graph.edges) : nothing
    native_retry=even_blossom && matching===nothing
    selected=matching===nothing ? Int[] : matching.selected
    if matching===nothing
        greedy_used=falses(size(triangles,2))
        for (index,pair) in enumerate(graph.physical)
            first=Int(pair.first_triangle);second=Int(pair.second_triangle)
            if !greedy_used[first] && !greedy_used[second] &&
               _gmsh_recombine_pair_measure(coords,pair.nodes)<recombine_angle
                push!(selected,index);greedy_used[first]=true;greedy_used[second]=true
            end
        end
    end
    used=falses(size(triangles,2))
    new_triangles=NTuple{3,Int32}[];new_triangle_tags=Int32[];triangle_births=UInt64[]
    raw_quads=NTuple{4,Int32}[];raw_tags=Int32[];raw_births=UInt64[]
    raw_pairs=Union{Nothing,_QuadCandidate}[]
    if existing_quadrangles!==nothing
        for column in axes(existing_quadrangles,2)
            push!(raw_quads,ntuple(row->Int32(existing_quadrangles[row,column]),4))
            push!(raw_tags,existing_quadrangle_tags[column])
            push!(raw_births,UInt64(0));push!(raw_pairs,nothing)
        end
    end
    pair_passes=native_retry ? 2 : 1
    counter=Base.checked_mul(UInt64(length(graph.physical)),UInt64(pair_passes))
    for index in selected
        index<=length(graph.physical) || continue
        pair=graph.physical[index]
        first=Int(pair.first_triangle);second=Int(pair.second_triangle)
        triangle_tags[first]==triangle_tags[second] || throw(ArgumentError(
            "recombine_triangles: primary optimal matching crosses a physical-tag boundary"))
        physical=triangle_tags[first]
        nodes=matching===nothing ? _primary_greedy_cycle(pair) :
            _primary_matched_cycle(triangles,pair)
        a,b,c,d=nodes
        cells=if matching!==nothing && embedded!==nothing && _edge_key(a,c) in embedded
            ((a,b,c),(c,d,a))
        elseif matching!==nothing && embedded!==nothing && _edge_key(b,d) in embedded
            ((b,c,d),(d,a,b))
        else
            nothing
        end
        if cells===nothing
            counter=Base.checked_add(counter,UInt64(1))
            push!(raw_quads,nodes);push!(raw_tags,physical);push!(raw_births,counter);push!(raw_pairs,pair)
        else
            _primary_supported_triangle_split(coords,triangles,pair,cells) || throw(ArgumentError(
                "recombine_triangles: primary embedded-diagonal materialization has no supported consistent projection"))
            for cell in cells
                counter=Base.checked_add(counter,UInt64(1))
                push!(new_triangles,cell);push!(new_triangle_tags,physical);push!(triangle_births,counter)
            end
        end
        used[first]=true;used[second]=true
    end
    if defer_postquality
        return _PrimaryWeightedEmission(used,new_triangles,new_triangle_tags,triangle_births,
            raw_quads,raw_tags,raw_births,counter,length(graph.physical),pair_passes)
    end
    quads=NTuple{4,Int32}[];quad_tags=Int32[];quad_births=UInt64[]
    for index in eachindex(raw_quads)
        nodes=raw_quads[index];physical=raw_tags[index];pair=raw_pairs[index]
        split=_primary_postquality_split(coords,nodes,minimum_quality)
        if split===nothing
            (_strict_convex_projection(coords,nodes) &&
             last(_quad_quality_range_result(coords,nodes))) || throw(ArgumentError(
                "recombine_triangles: primary optimal matching selected a quadrangle outside the supported convex geometry certificate"))
            push!(quads,nodes);push!(quad_tags,physical);push!(quad_births,raw_births[index])
        else
            supported=pair===nothing ? _primary_supported_quad_split(coords,split.triangles) :
                _primary_supported_triangle_split(coords,triangles,pair,split.triangles)
            supported || throw(ArgumentError(
                "recombine_triangles: primary postquality materialization has no supported consistent projection"))
            for (cell,offset) in zip(split.triangles,split.birth_offsets)
                push!(new_triangles,cell);push!(new_triangle_tags,physical)
                push!(triangle_births,Base.checked_add(counter,UInt64(offset)))
            end
            counter=Base.checked_add(counter,UInt64(split.reserved))
        end
    end
    return _PrimaryWeightedEmission(used,new_triangles,new_triangle_tags,triangle_births,
        quads,quad_tags,quad_births,counter,length(graph.physical),pair_passes)
end

function _primary_emission_product(mesh,emission::_PrimaryWeightedEmission;
                                   preserve_segments,physical_names,full_quad)
    retained=count(!,emission.used)
    triangle_count=retained+length(emission.triangles)
    full_quad && triangle_count!=0 && throw(ArgumentError(
        "recombine_triangles: full_quad requested but $triangle_count triangles remain in primary materialization"))
    blocks=Elements.ElementBlock[]
    births=Vector{UInt64}[]
    if preserve_segments && size(mesh.segs,2)>0
        push!(blocks,Elements.ElementBlock(1,mesh.segs,mesh.seg_tag))
        push!(births,zeros(UInt64,size(mesh.segs,2)))
    end
    if triangle_count>0
        cells=Matrix{Int32}(undef,3,triangle_count)
        tags=Vector{Int32}(undef,triangle_count)
        cell_births=zeros(UInt64,triangle_count)
        destination=0
        for column in axes(mesh.tris,2)
            emission.used[column] && continue
            destination+=1
            for row in 1:3;cells[row,destination]=mesh.tris[row,column];end
            tags[destination]=mesh.tri_tag[column]
        end
        for index in eachindex(emission.triangles)
            destination+=1
            for row in 1:3;cells[row,destination]=emission.triangles[index][row];end
            tags[destination]=emission.triangle_tags[index]
            cell_births[destination]=emission.triangle_births[index]
        end
        push!(blocks,Elements.ElementBlock(2,cells,tags));push!(births,cell_births)
    end
    if !isempty(emission.quadrangles)
        cells=Matrix{Int32}(undef,4,length(emission.quadrangles))
        for column in eachindex(emission.quadrangles),row in 1:4
            cells[row,column]=emission.quadrangles[column][row]
        end
        push!(blocks,Elements.ElementBlock(3,cells,emission.quadrangle_tags))
        push!(births,copy(emission.quadrangle_births))
    end
    result=Elements.MixedMesh(mesh.coords,blocks;physical_names)
    diagnostic=Elements.validate(result)
    diagnostic.ok || throw(ErrorException(
        "recombine_triangles: primary output validation failed — "*join(diagnostic.messages,"; ")))
    return (mesh=result,primary_pair_passes=emission.pair_passes,
        primary_factory_births=emission.factory_births,element_births=births)
end
