"""
    Recombine

Deterministic, topology-preserving triangle-to-quadrangle recombination for
validated surface meshes. The result uses Gmsh's fixed-node mixed-element model.
"""
module Recombine

using ..Predicates: orient2
import ..MeshTypes
import ..Elements

export recombine_triangles

struct _QuadCandidate
    first_triangle::Int32
    second_triangle::Int32
    nodes::NTuple{4,Int32}
    quality::Float64
    edge::NTuple{2,Int32}
end

@inline _edge_key(a::Int32,b::Int32)=a<b ? (a,b) : (b,a)

@inline function _third_vertex(triangle::NTuple{3,Int32},u::Int32,v::Int32)
    @inbounds for node in triangle
        node!=u && node!=v && return node
    end
    throw(ArgumentError("recombine_triangles: malformed shared triangle edge"))
end

@inline function _rotate_quad_minimum(nodes::NTuple{4,Int32})
    position=1
    @inbounds for i in 2:4
        nodes[i]<nodes[position] && (position=i)
    end
    return ntuple(i->nodes[mod1(position+i-1,4)],4)
end

@inline function _projected_point(coords,node::Int32,axes::NTuple{2,Int})
    i=Int(node)
    return (coords[axes[1],i],coords[axes[2],i])
end

function _strict_convex_projection(coords,nodes::NTuple{4,Int32})
    choices=((1,2),(2,3),(3,1))
    for axes in choices
        a=_projected_point(coords,nodes[1],axes)
        b=_projected_point(coords,nodes[2],axes)
        c=_projected_point(coords,nodes[3],axes)
        reference=orient2(a,b,c);reference==0 && continue
        convex=true
        @inbounds for i in 2:4
            a=_projected_point(coords,nodes[i],axes)
            b=_projected_point(coords,nodes[mod1(i+1,4)],axes)
            c=_projected_point(coords,nodes[mod1(i+2,4)],axes)
            if orient2(a,b,c)!=reference
                convex=false;break
            end
        end
        convex && return true
    end
    return false
end

@inline _sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _cross3(a,b)=(a[2]*b[3]-a[3]*b[2],
                      a[3]*b[1]-a[1]*b[3],
                      a[1]*b[2]-a[2]*b[1])
@inline _dot3(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _norm3(a)=hypot(a[1],a[2],a[3])

function _quad_quality(coords,nodes::NTuple{4,Int32})
    scale=0.0
    @inbounds for node in nodes,d in 1:3
        scale=max(scale,abs(coords[d,Int(node)]))
    end
    scale>0 || return 0.0
    points=ntuple(4) do i
        node=Int(nodes[i])
        (coords[1,node]/scale,coords[2,node]/scale,coords[3,node]/scale)
    end
    edges=ntuple(i->_sub3(points[mod1(i+1,4)],points[i]),4)
    lengths=ntuple(i->_norm3(edges[i]),4)
    minimum_length=minimum(lengths);maximum_length=maximum(lengths)
    minimum_length>0 && isfinite(maximum_length) || return 0.0
    minimum_sine=1.0
    @inbounds for i in 1:4
        previous=edges[mod1(i-1,4)]
        current=edges[i]
        sine=_norm3(_cross3(previous,current))/(lengths[mod1(i-1,4)]*lengths[i])
        minimum_sine=min(minimum_sine,sine)
    end
    normal1=_cross3(_sub3(points[2],points[1]),_sub3(points[3],points[1]))
    normal2=_cross3(_sub3(points[3],points[1]),_sub3(points[4],points[1]))
    norm1=_norm3(normal1);norm2=_norm3(normal2)
    norm1>0 && norm2>0 || return 0.0
    alignment=_dot3(normal1,normal2)/(norm1*norm2)
    alignment>0 || return 0.0
    quality=min(minimum_length/maximum_length,minimum_sine,min(alignment,1.0))
    return isfinite(quality) ? clamp(quality,0.0,1.0) : 0.0
end

function _candidate(coords,triangles,first_triangle::Int32,second_triangle::Int32,
                    first_u::Int32,first_v::Int32,
                    second_u::Int32,second_v::Int32,edge)
    # Consistently oriented neighboring triangles traverse their shared edge in
    # opposite directions. Leaving an inconsistent pair uncombined preserves
    # the input rather than silently emitting an inside-out quadrangle.
    (second_u==first_v && second_v==first_u) || return nothing
    first=Int(first_triangle);second=Int(second_triangle)
    triangle1=(triangles[1,first],triangles[2,first],triangles[3,first])
    triangle2=(triangles[1,second],triangles[2,second],triangles[3,second])
    opposite1=_third_vertex(triangle1,first_u,first_v)
    opposite2=_third_vertex(triangle2,first_u,first_v)
    nodes=_rotate_quad_minimum((first_u,opposite2,first_v,opposite1))
    _strict_convex_projection(coords,nodes) || return nothing
    quality=_quad_quality(coords,nodes)
    quality>0 || return nothing
    return _QuadCandidate(first_triangle,second_triangle,nodes,quality,edge)
end

function _recombine_quality(value)
    value isa Bool && throw(ArgumentError(
        "recombine_triangles: min_quality must not be Bool"))
    value isa Real || throw(ArgumentError(
        "recombine_triangles: min_quality must be real"))
    quality=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "recombine_triangles: min_quality must be Float64-representable"))
    end
    isfinite(quality) || throw(ArgumentError(
        "recombine_triangles: min_quality must be finite"))
    0<=quality<=1 || throw(ArgumentError(
        "recombine_triangles: min_quality must lie in 0:1"))
    return quality
end

"""
    recombine_triangles(mesh; min_quality=0, preserve_segments=true,
                        physical_names=Dict(), algorithm=:greedy,
                        full_quad=false) -> MixedMesh

Pair adjacent, consistently oriented triangles into four-node Gmsh quadrangles.
`:greedy` accepts candidates in deterministic shape-score order. `:blossom` uses
Edmonds maximum-cardinality matching on the dual of eligible same-tag pairs,
trying neighbors in that same quality order. `full_quad=true` requires
`:blossom` and throws if any triangle remains unmatched. Only strictly convex
projected quadrangles with matching physical tags and quality at least
`min_quality` are accepted. Unpaired triangles are retained unless `full_quad`
is requested. Segment connectivity and all per-cell physical tags are preserved
by default.

The input must be a validated surface/curve mesh without tetrahedra. This is a
surface recombination operation; it does not generate a structured grid or modify
the geometry and it never pairs triangles across a physical-tag boundary.
"""
function recombine_triangles(mesh::MeshTypes.Mesh;min_quality=0.0,
                             preserve_segments=true,
                             physical_names=Dict{Tuple{Int,Int},String}(),
                             algorithm=:greedy,
                             full_quad=false)
    preserve_segments isa Bool || throw(ArgumentError(
        "recombine_triangles: preserve_segments must be Bool"))
    full_quad isa Bool || throw(ArgumentError(
        "recombine_triangles: full_quad must be Bool"))
    algorithm isa Symbol || throw(ArgumentError(
        "recombine_triangles: algorithm must be :greedy or :blossom"))
    algorithm in (:greedy,:blossom) || throw(ArgumentError(
        "recombine_triangles: algorithm must be :greedy or :blossom"))
    full_quad && algorithm!==:blossom && throw(ArgumentError(
        "recombine_triangles: full_quad requires algorithm=:blossom"))
    threshold=_recombine_quality(min_quality)
    size(mesh.tets,2)==0 || throw(ArgumentError(
        "recombine_triangles: input must not contain tetrahedra"))
    diagnostic=MeshTypes.validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "recombine_triangles: input mesh is invalid — "*
        join(diagnostic.messages,"; ")))

    owners=Dict{NTuple{2,Int32},NTuple{3,Int32}}()
    completed=Set{NTuple{2,Int32}}()
    candidates=_QuadCandidate[]
    triangles=mesh.tris
    @inbounds for triangle_index in axes(triangles,2)
        triangle=Int32(triangle_index)
        for local_edge in ((1,2),(2,3),(3,1))
            u=triangles[local_edge[1],triangle_index]
            v=triangles[local_edge[2],triangle_index]
            edge=_edge_key(u,v)
            if edge in completed
                throw(ArgumentError(
                    "recombine_triangles: non-manifold edge $edge has more than two incident triangles"))
            elseif haskey(owners,edge)
                first_triangle,first_u,first_v=owners[edge]
                candidate=_candidate(mesh.coords,triangles,first_triangle,triangle,
                                     first_u,first_v,u,v,edge)
                candidate===nothing || push!(candidates,candidate)
                push!(completed,edge)
            else
                owners[edge]=(triangle,u,v)
            end
        end
    end
    sort!(candidates;by=c->(-c.quality,c.edge,c.first_triangle,c.second_triangle),
          alg=MergeSort)

    used=falses(size(triangles,2));accepted=_QuadCandidate[]
    if algorithm===:blossom
        n=size(triangles,2)
        adj=[Int[] for _ in 1:n]
        bypair=Dict{Tuple{Int,Int},_QuadCandidate}()
        for candidate in candidates
            a=Int(candidate.first_triangle);b=Int(candidate.second_triangle)
            mesh.tri_tag[a]==mesh.tri_tag[b] || continue
            candidate.quality>=threshold || continue
            push!(adj[a],b); push!(adj[b],a)
            bypair[(min(a,b),max(a,b))]=candidate
        end
        for v in 1:n
            sort!(adj[v]; by=u->begin
                c=bypair[(min(v,u),max(v,u))]
                (-c.quality,c.edge,c.first_triangle,c.second_triangle)
            end)
        end
        mate=_edmonds_matching(n,adj)
        for v in 1:n
            u=mate[v]
            u>v || continue
            candidate=bypair[(v,u)]
            used[v]=true; used[u]=true
            push!(accepted,candidate)
        end
        if full_quad && any(!,used)
            leftover=count(!,used)
            throw(ArgumentError(
                "recombine_triangles: full_quad requested but $leftover triangles remain unmatched"))
        end
    else
        for candidate in candidates
            first=Int(candidate.first_triangle);second=Int(candidate.second_triangle)
            if !used[first] && !used[second] &&
               mesh.tri_tag[first]==mesh.tri_tag[second] &&
               candidate.quality>=threshold
                used[first]=true;used[second]=true;push!(accepted,candidate)
            end
        end
    end

    blocks=Elements.ElementBlock[]
    if preserve_segments && size(mesh.segs,2)>0
        push!(blocks,Elements.ElementBlock(1,mesh.segs,mesh.seg_tag))
    end
    remaining=count(!,used)
    if remaining>0
        nodes=Matrix{Int32}(undef,3,remaining);tags=Vector{Int32}(undef,remaining)
        destination=0
        @inbounds for triangle in axes(triangles,2)
            used[triangle] && continue
            destination+=1
            for local_node in 1:3
                nodes[local_node,destination]=triangles[local_node,triangle]
            end
            tags[destination]=mesh.tri_tag[triangle]
        end
        push!(blocks,Elements.ElementBlock(2,nodes,tags))
    end
    if !isempty(accepted)
        nodes=Matrix{Int32}(undef,4,length(accepted))
        tags=Vector{Int32}(undef,length(accepted))
        @inbounds for (i,candidate) in pairs(accepted)
            for local_node in 1:4
                nodes[local_node,i]=candidate.nodes[local_node]
            end
            tags[i]=mesh.tri_tag[Int(candidate.first_triangle)]
        end
        push!(blocks,Elements.ElementBlock(3,nodes,tags))
    end
    result=Elements.MixedMesh(mesh.coords,blocks;physical_names=physical_names)
    output_diagnostic=Elements.validate(result)
    output_diagnostic.ok || throw(ErrorException(
        "recombine_triangles: internal output validation failed — "*
        join(output_diagnostic.messages,"; ")))
    return result
end

# Edmonds' blossom algorithm: maximum-cardinality matching on a general graph.
# Neighbors are tried in the given order so the matching is deterministic.
function _edmonds_matching(n::Int, adj::Vector{Vector{Int}})
    mate=zeros(Int,n)
    function lca(a::Int,b::Int,base::Vector{Int},parent::Vector{Int})
        seen=falses(n)
        while true
            a=base[a]
            seen[a]=true
            mate[a]==0 && break
            a=base[parent[mate[a]]]
        end
        while true
            b=base[b]
            seen[b] && return b
            b=base[parent[mate[b]]]
        end
    end
    function mark_path!(blossom,v::Int,b::Int,children::Int,
                        base::Vector{Int},parent::Vector{Int})
        while base[v]!=b
            blossom[base[v]]=true
            blossom[base[mate[v]]]=true
            parent[v]=children
            children=mate[v]
            v=parent[mate[v]]
        end
    end
    function augment(root::Int)
        used=falses(n)
        parent=zeros(Int,n)
        base=collect(1:n)
        q=Int[root]; used[root]=true; head=1
        while head<=length(q)
            v=q[head]; head+=1
            for u in adj[v]
                if base[v]==base[u] || mate[v]==u
                    continue
                elseif u==root || (mate[u]!=0 && parent[mate[u]]!=0)
                    b=lca(v,u,base,parent)
                    blossom=falses(n)
                    mark_path!(blossom,v,b,u,base,parent)
                    mark_path!(blossom,u,b,v,base,parent)
                    for i in 1:n
                        blossom[base[i]] || continue
                        base[i]=b
                        used[i] && continue
                        used[i]=true
                        push!(q,i)
                    end
                elseif parent[u]==0
                    parent[u]=v
                    if mate[u]==0
                        # reconstruct augmenting path
                        while u!=0
                            pv=parent[u]; nu=mate[pv]
                            mate[u]=pv; mate[pv]=u
                            u=nu
                        end
                        return true
                    end
                    matched=mate[u]
                    push!(q,matched); used[matched]=true
                end
            end
        end
        return false
    end
    grew=true
    while grew
        grew=false
        for v in 1:n
            mate[v]==0 && augment(v) && (grew=true)
        end
    end
    return mate
end

end # module Recombine
