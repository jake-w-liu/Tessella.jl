"""
    Recombine

Deterministic, topology-preserving triangle-to-quadrangle recombination for
validated surface meshes. The result uses Gmsh's fixed-node mixed-element model.
"""
module Recombine

using ..Predicates: orient2
import ..MeshTypes
import ..Elements
using ..GmshLibm: _gm87_atan2

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
    # explicit rotation: a closure over the reassigned `position` would box it
    @inbounds return (nodes[position],nodes[mod1(position+1,4)],
                      nodes[mod1(position+2,4)],nodes[mod1(position+3,4)])
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

@inline _div3(a,denominator)=(a[1]/denominator,a[2]/denominator,a[3]/denominator)
@inline _maxabs3(a)=max(abs(a[1]),abs(a[2]),abs(a[3]))

# Compare divide-before-subtract spans with differences of the represented
# input coordinates. Keep healthy old arithmetic within sixteen rounding
# units; cancellation beyond that bound uses the local geometry workspace.
@inline function _quad_span_cancellation(coords,a::Int32,b::Int32,span,scale)
    first=Int(a);second=Int(b)
    reference=((coords[1,second]-coords[1,first])/scale,
               (coords[2,second]-coords[2,first])/scale,
               (coords[3,second]-coords[3,first])/scale)
    span==reference && return false
    length=_norm3(reference)
    return !isfinite(length) || _norm3(_sub3(span,reference))>16eps(Float64)*length
end

# Exceptional-range geometry uses each raw edge independently. A quality score
# may round to zero while represented geometry remains nondegenerate and valid.
function _quad_local_vectors(coords,nodes::NTuple{4,Int32})
    n1=Int(nodes[1]);n2=Int(nodes[2]);n3=Int(nodes[3]);n4=Int(nodes[4])
    p1=(coords[1,n1],coords[2,n1],coords[3,n1])
    p2=(coords[1,n2],coords[2,n2],coords[3,n2])
    p3=(coords[1,n3],coords[2,n3],coords[3,n3])
    p4=(coords[1,n4],coords[2,n4],coords[3,n4])
    edges=(_sub3(p2,p1),_sub3(p3,p2),_sub3(p4,p3),_sub3(p1,p4))
    if !all(edge->all(isfinite,edge),edges)
        coordinate_scale=max(_maxabs3(p1),_maxabs3(p2),_maxabs3(p3),_maxabs3(p4))
        coordinate_scale>0 && isfinite(coordinate_scale) || return nothing
        q1=_div3(p1,coordinate_scale);q2=_div3(p2,coordinate_scale)
        q3=_div3(p3,coordinate_scale);q4=_div3(p4,coordinate_scale)
        return (_sub3(q2,q1),_sub3(q3,q2),_sub3(q4,q3),_sub3(q1,q4))
    end
    return edges
end

@inline function _quad_edge_measure(edge)
    scale=_maxabs3(edge)
    scale>0 && isfinite(scale) || return ((0.0,0.0,0.0),(0,0.0),false)
    scaled=_div3(edge,scale);length=_norm3(scaled)
    unit=_div3(scaled,length)
    mantissa,exponent=frexp(scale)
    mantissa*=length
    if mantissa>=1.0
        mantissa*=0.5;exponent+=1
    end
    return (unit,(exponent,mantissa),true)
end

function _quad_quality_exact_range(coords,nodes::NTuple{4,Int32})
    R=Rational{BigInt}
    points=ntuple(i->ntuple(d->R(coords[d,Int(nodes[i])]),3),4)
    edges=ntuple(i->_sub3(points[mod1(i+1,4)],points[i]),4)
    squared=ntuple(i->_dot3(edges[i],edges[i]),4)
    minimum(squared)>0 || return (0.0,false)
    corners=ntuple(i->_cross3(edges[i],edges[mod1(i+1,4)]),4)
    normal_squared=ntuple(i->_dot3(corners[i],corners[i]),4)
    minimum(normal_squared)>0 || return (0.0,false)
    alignment=_dot3(corners[1],corners[3])
    alignment>0 || return (0.0,false)
    squared_quality=min(minimum(squared)/maximum(squared),
        minimum(ntuple(i->normal_squared[i]/(squared[i]*squared[mod1(i+1,4)]),4)),
        alignment^2/(normal_squared[1]*normal_squared[3]),one(R))
    score=setrounding(BigFloat,RoundNearest) do
        setprecision(BigFloat,256) do
            Float64(sqrt(BigFloat(squared_quality)),RoundNearest)
        end
    end
    return (score,true)
end

function _quad_quality_range_result(coords,nodes::NTuple{4,Int32})
    edges=_quad_local_vectors(coords,nodes)
    edges===nothing && return (0.0,false)
    measures=ntuple(i->_quad_edge_measure(edges[i]),4)
    all(measure->measure[3],measures) || return _quad_quality_exact_range(coords,nodes)
    units=ntuple(i->measures[i][1],4)
    lengths=ntuple(i->measures[i][2],4)
    corners=ntuple(i->_cross3(units[i],units[mod1(i+1,4)]),4)
    sines=ntuple(i->_norm3(corners[i]),4)
    # A zero or cancellation-dominated unit cross needs exact represented
    # geometry. Its positive sine can itself round to zero, independently of
    # the validity of the two triangles and their alignment.
    for i in 1:4
        first=units[i];second=units[mod1(i+1,4)]
        permanent=(abs(first[2]*second[3])+abs(first[3]*second[2]))+
            (abs(first[3]*second[1])+abs(first[1]*second[3]))+
            (abs(first[1]*second[2])+abs(first[2]*second[1]))
        sines[i]>64eps(Float64)*permanent || return _quad_quality_exact_range(coords,nodes)
    end
    normal1=_div3(corners[1],sines[1]);normal2=_div3(corners[3],sines[3])
    alignment=_dot3(normal1,normal2)
    permanent=sum(abs(normal1[i]*normal2[i]) for i in 1:3)
    abs(alignment)>64eps(Float64)*permanent || return _quad_quality_exact_range(coords,nodes)
    alignment>0 || return (0.0,false)
    low=minimum(lengths);high=maximum(lengths)
    ratio=ldexp(low[2]/high[2],low[1]-high[1])
    quality=min(ratio,minimum(sines),min(alignment,1.0))
    return (clamp(quality,0.0,1.0),true)
end

_quad_quality_range(coords,nodes::NTuple{4,Int32})=
    first(_quad_quality_range_result(coords,nodes))


function _quad_quality(coords,nodes::NTuple{4,Int32})
    scale=0.0
    @inbounds for node in nodes,d in 1:3
        scale=max(scale,abs(coords[d,Int(node)]))
    end
    scale>0 || return 0.0
    # explicit points: a closure over the reassigned `scale` would box it
    n1=Int(nodes[1]);n2=Int(nodes[2]);n3=Int(nodes[3]);n4=Int(nodes[4])
    points=((coords[1,n1]/scale,coords[2,n1]/scale,coords[3,n1]/scale),
            (coords[1,n2]/scale,coords[2,n2]/scale,coords[3,n2]/scale),
            (coords[1,n3]/scale,coords[2,n3]/scale,coords[3,n3]/scale),
            (coords[1,n4]/scale,coords[2,n4]/scale,coords[3,n4]/scale))
    edges=ntuple(i->_sub3(points[mod1(i+1,4)],points[i]),4)
    @inbounds for i in 1:4
        _quad_span_cancellation(coords,nodes[i],nodes[mod1(i+1,4)],edges[i],scale) &&
            return _quad_quality_range(coords,nodes)
    end
    (_quad_span_cancellation(coords,nodes[1],nodes[3],_sub3(points[3],points[1]),scale) ||
     _quad_span_cancellation(coords,nodes[1],nodes[4],_sub3(points[4],points[1]),scale)) &&
        return _quad_quality_range(coords,nodes)
    lengths=ntuple(i->_norm3(edges[i]),4)
    minimum_length=minimum(lengths);maximum_length=maximum(lengths)
    minimum_length>0 && isfinite(maximum_length) || return _quad_quality_range(coords,nodes)
    minimum_sine=1.0
    @inbounds for i in 1:4
        previous=edges[mod1(i-1,4)]
        current=edges[i]
        denominator=lengths[mod1(i-1,4)]*lengths[i]
        denominator>=floatmin(Float64) || return _quad_quality_range(coords,nodes)
        sine=_norm3(_cross3(previous,current))/denominator
        minimum_sine=min(minimum_sine,sine)
    end
    normal1=_cross3(_sub3(points[2],points[1]),_sub3(points[3],points[1]))
    normal2=_cross3(_sub3(points[3],points[1]),_sub3(points[4],points[1]))
    norm1=_norm3(normal1);norm2=_norm3(normal2)
    norm1>0 && norm2>0 || return _quad_quality_range(coords,nodes)
    denominator=norm1*norm2
    denominator>=floatmin(Float64) || return _quad_quality_range(coords,nodes)
    alignment=_dot3(normal1,normal2)/denominator
    alignment>0 || return 0.0
    quality=min(minimum_length/maximum_length,minimum_sine,min(alignment,1.0))
    return isfinite(quality) && quality>0 ? clamp(quality,0.0,1.0) : _quad_quality_range(coords,nodes)
end

# Pinned Gmsh qualityMeasures.cpp qmQuadrangle::eta; RecombineTriangle
# uses this signed corner-angle measure for the strict greedy angle admission.
function _gmsh_recombine_pair_measure(coords,nodes::NTuple{4,Int32})
    # The compatibility measure deliberately uses raw represented coordinates.
    # Gmsh's cross-product norm can underflow or overflow; normalizing these
    # vectors would change its strict angle-admission decision.
    n1=Int(nodes[1]);n2=Int(nodes[2]);n3=Int(nodes[3]);n4=Int(nodes[4])
    points=((coords[1,n1],coords[2,n1],coords[3,n1]),
            (coords[1,n2],coords[2,n2],coords[3,n2]),
            (coords[1,n3],coords[2,n3],coords[3,n3]),
            (coords[1,n4],coords[2,n4],coords[3,n4]))
    edges=ntuple(i->_sub3(points[mod1(i+1,4)],points[i]),4)
    corners=ntuple(i->_cross3(edges[i],edges[mod1(i+1,4)]),4)
    deviation=0.0
    for i in 1:4
        first=_sub3(points[mod1(i-1,4)],points[i])
        second=_sub3(points[mod1(i+1,4)],points[i])
        normal=_cross3(first,second)
        angle=180.0*_gm87_atan2(sqrt(_dot3(normal,normal)),_dot3(first,second))/pi
        # std::min(180., angle) keeps its first operand when angle is NaN.
        angle=angle<180.0 ? angle : 180.0
        deviation=max(deviation,abs(90.0-angle))
    end
    sign=1.0
    for i in 2:4
        _dot3(corners[1],corners[i])<0 && (sign=-1.0)
    end
    return sign*(1.0-deviation*(1.0/90.0))
end

@inline function _gmsh_priority_corner(coords,a::Int32,b::Int32,c::Int32)
    point=(coords[1,Int(b)],coords[2,Int(b)],coords[3,Int(b)])
    first=_sub3((coords[1,Int(a)],coords[2,Int(a)],coords[3,Int(a)]),point)
    second=_sub3((coords[1,Int(c)],coords[2,Int(c)],coords[3,Int(c)]),point)
    normal=_cross3(first,second)
    return abs(90.0-180.0*_gm87_atan2(sqrt(_dot3(normal,normal)),_dot3(first,second))/pi)
end

function _gmsh_recombine_pair_priority(coords,triangles,candidate::_QuadCandidate)
    first=Int(candidate.first_triangle);second=Int(candidate.second_triangle)
    triangle1=(triangles[1,first],triangles[2,first],triangles[3,first])
    triangle2=(triangles[1,second],triangles[2,second],triangles[3,second])
    directed=_gmsh_first_candidate_edge(triangle1,candidate.edge)
    n1,n2=directed
    n3=_third_vertex(triangle1,n1,n2);n4=_third_vertex(triangle2,n1,n2)
    # MEdge retains the first triangle's directed edge. The constructor visits
    # n4,n2,n3,n1, and std::max(new,old) returns new when either is NaN.
    quality=_gmsh_priority_corner(coords,n1,n4,n2)
    for deviation in (_gmsh_priority_corner(coords,n4,n2,n3),
            _gmsh_priority_corner(coords,n2,n3,n1),
            _gmsh_priority_corner(coords,n3,n1,n4))
        quality=deviation<quality ? quality : deviation
    end
    return quality
end

@inline function _gmsh_first_candidate_edge(triangle,edge)
    for (i,j) in ((1,2),(2,3),(3,1))
        _edge_key(triangle[i],triangle[j])==edge && return (triangle[i],triangle[j])
    end
    throw(ErrorException("recombine_triangles: candidate edge is absent from its first triangle"))
end

struct _GmshPairOrder
    priority::Float64
    candidate::_QuadCandidate
end
@inline _gmsh_pair_less(a,b)=a.priority<b.priority

include("GmshPairSort.jl")

function _gmsh_sort_candidates!(candidates,coords,triangles)
    sort!(candidates;by=c->c.edge,alg=MergeSort)
    ordered=Vector{_GmshPairOrder}(undef,length(candidates))
    for index in eachindex(candidates)
        candidate=candidates[index]
        ordered[index]=_GmshPairOrder(_gmsh_recombine_pair_priority(coords,triangles,candidate),candidate)
    end
    _gmsh_pair_sort!(ordered)
    for index in eachindex(candidates)
        candidates[index]=ordered[index].candidate
    end
    return candidates
end

function _recombine_angle(value)
    value===nothing && return nothing
    value isa Real && !(value isa Bool) || throw(ArgumentError(
        "recombine_triangles: recombine_angle must be real"))
    angle=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("recombine_triangles: recombine_angle must be Float64-representable"))
    end
    isfinite(angle) && 0<=angle<=90 || throw(ArgumentError(
        "recombine_triangles: recombine_angle must lie in 0:90"))
    return angle
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
    if quality==0.0
        quality,valid=_quad_quality_range_result(coords,nodes)
        valid || return nothing
    end
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
by default. `recombine_angle` optionally applies the pinned Gmsh strict
signed corner-angle measure admission to greedy pairing. With an explicit angle,
an odd triangle count makes `:blossom` fall back to that greedy pass, as in pinned
Gmsh. Even-count Blossom matching ignores the angle when its native boundary
closure graph admits a perfect matching; otherwise it retries greedy pairing.
Without an angle, `:blossom`
retains the standalone maximum-cardinality matching contract.

The input must be a validated surface/curve mesh without tetrahedra. This is a
surface recombination operation; it does not generate a structured grid or modify
the geometry and it never pairs triangles across a physical-tag boundary.
"""
function recombine_triangles(mesh::MeshTypes.Mesh;min_quality=0.0,
                             preserve_segments=true,
                             physical_names=Dict{Tuple{Int,Int},String}(),
                             algorithm=:greedy,
                             full_quad=false,
                             protected_edges=nothing,recombine_angle=nothing)
    return _recombine_triangles(mesh;min_quality,preserve_segments,physical_names,
        algorithm,full_quad,protected_edges,recombine_angle).mesh
end

# Keep primary pair-pass history available to the API label allocator while
# preserving the public function's detached MixedMesh return contract.
function _recombine_triangles(mesh::MeshTypes.Mesh;min_quality=0.0,
                              preserve_segments=true,
                              physical_names=Dict{Tuple{Int,Int},String}(),
                              algorithm=:greedy,full_quad=false,
                              protected_edges=nothing,recombine_angle=nothing,
                              embedded_edges=nothing)
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
    protected = protected_edges===nothing ? nothing :
        begin
            protected_edges isa AbstractSet || throw(ArgumentError(
                "recombine_triangles: protected_edges must be a set of " *
                "two-node edge tuples"))
            protected_edges
        end
    threshold=_recombine_quality(min_quality)
    angle=_recombine_angle(recombine_angle)
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
                (candidate===nothing ||
                 (protected!==nothing && edge in protected)) ||
                    push!(candidates,candidate)
                push!(completed,edge)
            else
                owners[edge]=(triangle,u,v)
            end
        end
    end
    if angle===nothing
        sort!(candidates;by=c->(-c.quality,c.edge,c.first_triangle,c.second_triangle),alg=MergeSort)
    else
        _gmsh_sort_candidates!(candidates,mesh.coords,triangles)
    end

    used=falses(size(triangles,2));accepted=_QuadCandidate[]
    primary_retry=false
    # Pinned Gmsh falls back to its angle-filtered greedy pass when the
    # triangle count is odd. With no angle, keep the standalone matching API.
    if algorithm===:blossom && (angle===nothing || iseven(size(triangles,2)))
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
        if angle!==nothing && any(iszero,mate)
            # An internal-only unmatched triangle does not imply primary
            # failure: native Blossom also contains high-cost boundary links.
            closure=_primary_closure_graph(triangles,protected,embedded_edges)
            primary_retry=any(iszero,_edmonds_matching(n,closure))
        end
        if !primary_retry
            for v in 1:n
                u=mate[v]
                if angle===nothing
                    u>v || continue
                    candidate=bypair[(v,u)]
                else
                    # Native match publication visits the larger triangle
                    # endpoint. Mate selection and the standalone contract
                    # retain their original ordering and behavior.
                    0<u<v || continue
                    candidate=bypair[(u,v)]
                end
                used[v]=true; used[u]=true
                push!(accepted,candidate)
            end
        end
    end
    if algorithm!==:blossom ||
       (angle!==nothing && isodd(size(triangles,2))) || primary_retry
        for candidate in candidates
            first=Int(candidate.first_triangle);second=Int(candidate.second_triangle)
            if !used[first] && !used[second] &&
               mesh.tri_tag[first]==mesh.tri_tag[second] &&
               candidate.quality>=threshold &&
               (angle===nothing || _gmsh_recombine_pair_measure(mesh.coords,candidate.nodes)<angle)
                used[first]=true;used[second]=true;push!(accepted,candidate)
            end
        end
    end

    if full_quad && any(!,used)
        leftover=count(!,used)
        throw(ArgumentError(
            "recombine_triangles: full_quad requested but $leftover triangles remain unmatched"))
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
    return (mesh=result,primary_pair_passes=primary_retry ? 2 : 1)
end

# The pinned graph includes all unprotected internal triangle adjacencies,
# including pairs screened out by the safer geometric/min_quality contract.
# Its additional boundary links join the two boundary triangles incident at
# each vertex; they are matching edges only, never physical quadrangles.
function _primary_closure_graph(triangles,protected,embedded)
    n=size(triangles,2)
    adjacency=[Int[] for _ in 1:n]
    owners=Dict{NTuple{2,Int32},Tuple{Int,Int}}()
    @inbounds for triangle in axes(triangles,2)
        for (i,j) in ((1,2),(2,3),(3,1))
            edge=_edge_key(triangles[i,triangle],triangles[j,triangle])
            first,_=get(owners,edge,(0,0))
            owners[edge]=first==0 ? (triangle,0) : (first,triangle)
        end
    end
    periodic=Dict{Int32,Tuple{Int,Int}}()
    for edge in sort!(collect(keys(owners)))
        first,second=owners[edge]
        if second!=0
            protected!==nothing && edge in protected && continue
            push!(adjacency[first],second);push!(adjacency[second],first)
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
    for (first,second) in values(periodic)
        # Native t2n[nullptr] value-initializes to triangle index zero. Preserve
        # that pinned graph convention for a lone boundary incidence.
        other=second==0 ? 1 : second
        first==other && continue
        push!(adjacency[first],other);push!(adjacency[other],first)
    end
    return adjacency
end

# Edmonds' blossom algorithm: maximum-cardinality matching on a general graph.
# Neighbors are tried in the given order so the matching is deterministic.
function _edmonds_matching(n::Int, adj::Vector{Vector{Int}})
    mate=zeros(Int,n)
    # A matching search uses linear scratch storage. Reuse it across roots
    # and contractions rather than allocating full-n arrays on every search.
    seen=falses(n)
    used=falses(n)
    parent=zeros(Int,n)
    base=collect(1:n)
    blossom=falses(n)
    q=Int[]
    sizehint!(q,n)
    function lca(a::Int,b::Int,base::Vector{Int},parent::Vector{Int})
        fill!(seen,false)
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
        fill!(used,false)
        fill!(parent,0)
        for vertex in 1:n
            base[vertex]=vertex
        end
        empty!(q);push!(q,root);used[root]=true;head=1
        while head<=length(q)
            v=q[head]; head+=1
            for u in adj[v]
                if base[v]==base[u] || mate[v]==u
                    continue
                elseif u==root || (mate[u]!=0 && parent[mate[u]]!=0)
                    b=lca(v,u,base,parent)
                    fill!(blossom,false)
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
