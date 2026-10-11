# Ordered integer matching for bipartite graphs. This is an independent
# primal-dual implementation: feasible singleton duals, tight-edge alternating
# trees, and a full primal/dual certificate. It does not adapt Concorde code.
# General graphs use the separately licensed BSD weighted Blossom kernel.
module OrderedMatching

using ..WeightedMatching: minimum_weight_matching, _SeededInitial

struct _OrderedGraph
    edges::Vector{Tuple{Int,Int,Int64}}
    incident::Vector{Vector{Int}}
end

@inline function _other(g::_OrderedGraph, edge::Int, vertex::Int)
    a,b,_=g.edges[edge]
    return vertex==a ? b : a
end

@inline function _slack(g::_OrderedGraph, dual, edge::Int)
    a,b,cost=g.edges[edge]
    return Base.checked_sub(Base.checked_mul(Int64(2),cost),
        Base.checked_add(dual[a],dual[b]))
end

function _graph(count::Int, input)
    count>=0 || throw(ArgumentError("ordered matching: count must be nonnegative"))
    count<=typemax(Int)÷2 || throw(ArgumentError("ordered matching: workspace count overflows Int"))
    edges=Tuple{Int,Int,Int64}[]
    incident=[Int[] for _ in 1:count]
    if Base.IteratorSize(typeof(input)) isa Union{Base.HasLength,Base.HasShape}
        sizehint!(edges,length(input))
    end
    limit=typemax(Int64)÷(16max(Int128(1),Int128(count)))
    for (a,b,cost) in input
        a isa Integer && b isa Integer && !(a isa Bool) && !(b isa Bool) ||
            throw(ArgumentError("ordered matching: endpoints must be integer vertices"))
        1<=a<=count && 1<=b<=count || throw(ArgumentError("ordered matching: vertex out of range"))
        cost isa Integer && !(cost isa Bool) || throw(ArgumentError("ordered matching: costs must be integers"))
        -limit<=cost<=limit || throw(ArgumentError("ordered matching: integer cost exceeds checked dual range"))
        a=Int(a);b=Int(b)
        # Preserve every input edge and its direction, including equal-cost
        # parallel edges. A self-loop cannot cover two matching vertices.
        push!(edges,(a,b,Int64(cost)))
        a==b && continue
        push!(incident[a],length(edges));push!(incident[b],length(edges))
    end
    return _OrderedGraph(edges,incident)
end

function _bipartite(g::_OrderedGraph)
    count=length(g.incident)
    color=zeros(Int8,count);queue=Int[];sizehint!(queue,count)
    for root in 1:count
        color[root]==0 || continue
        empty!(queue);push!(queue,root);color[root]=1;head=1
        while head<=length(queue)
            vertex=queue[head];head+=1
            for edge in g.incident[vertex]
                other=_other(g,edge,vertex)
                if color[other]==0
                    color[other]=-color[vertex];push!(queue,other)
                elseif color[other]==color[vertex]
                    return false
                end
            end
        end
    end
    return true
end

function _certificate(g::_OrderedGraph, mate_edge, dual)
    primal=Int128(0)
    for vertex in eachindex(mate_edge)
        edge=mate_edge[vertex]
        @assert edge!=0
        other=_other(g,edge,vertex)
        @assert other!=vertex && mate_edge[other]==edge && _slack(g,dual,edge)==0
        a,b,cost=g.edges[edge]
        vertex==b && (primal+=cost)
    end
    for edge in eachindex(g.edges)
        a,b,_=g.edges[edge]
        a==b || @assert _slack(g,dual,edge)>=0
    end
    @assert 2primal==sum(Int128(value) for value in dual;init=Int128(0))
    return nothing
end

function _bipartite_perfect(g::_OrderedGraph; reject_odd_cycles::Bool=false)
    count=length(g.incident)
    any(isempty,g.incident) && return nothing
    dual=Vector{Int64}(undef,count)
    mate_edge=zeros(Int,count)
    for vertex in 1:count
        dual[vertex]=minimum(g.edges[edge][3] for edge in g.incident[vertex])
    end
    # Visit vertices in input order and incident edges in reverse insertion
    # order. Tightening an unmatched singleton preserves dual feasibility.
    for vertex in 1:count
        mate_edge[vertex]==0 || continue
        delta=typemax(Int64)
        for edge in Iterators.reverse(g.incident[vertex])
            slack=_slack(g,dual,edge)
            delta=min(delta,slack)
            other=_other(g,edge,vertex)
            if slack==0 && mate_edge[vertex]==0 && mate_edge[other]==0
                mate_edge[vertex]=edge;mate_edge[other]=edge
            end
        end
        dual[vertex]=Base.checked_add(dual[vertex],delta)
    end
    label=zeros(Int8,count);parent=zeros(Int,count);parent_edge=zeros(Int,count)
    children=[Int[] for _ in 1:count]
    stack=Int[];visited=Int[];sizehint!(stack,count);sizehint!(visited,count)
    for root in 1:count
        mate_edge[root]==0 || continue
        # Only vertices reached by the previous alternating tree can own
        # dirty labels, parent edges or children. Clear them before reuse.
        for vertex in visited
            label[vertex]=0;parent[vertex]=0;parent_edge[vertex]=0
            empty!(children[vertex])
        end
        empty!(visited)
        label[root]=1;push!(visited,root)
        endpoint=0
        while endpoint==0
            empty!(stack);push!(stack,root)
            while !isempty(stack) && endpoint==0
                vertex=pop!(stack)
                if label[vertex]==1
                    for edge in Iterators.reverse(g.incident[vertex])
                        _slack(g,dual,edge)==0 || continue
                        other=_other(g,edge,vertex)
                        # A tight edge between two outer vertices closes an
                        # odd alternating cycle. The singleton tree cannot
                        # certify its contraction; let the Blossom kernel
                        # handle that graph from its original inputs.
                        reject_odd_cycles && label[other]==1 && return nothing
                        label[other]==0 || continue
                        parent[other]=vertex;parent_edge[other]=edge
                        label[other]=2;push!(visited,other)
                        push!(children[vertex],other)
                        if mate_edge[other]==0
                            endpoint=other;break
                        end
                        partner=_other(g,mate_edge[other],other)
                        @assert label[partner]==0
                        parent[partner]=other;parent_edge[partner]=mate_edge[other]
                        label[partner]=1;push!(visited,partner)
                        push!(children[other],partner)
                    end
                end
                append!(stack,children[vertex])
            end
            endpoint!=0 && break
            delta=typemax(Int64)
            for vertex in visited
                label[vertex]==1 || continue
                for edge in g.incident[vertex]
                    other=_other(g,edge,vertex)
                    label[other]==0 || continue
                    delta=min(delta,_slack(g,dual,edge))
                end
            end
            delta==typemax(Int64) && return nothing
            @assert delta>0
            if reject_odd_cycles
                # Both outer singleton duals grow by delta. Keep all such
                # edges feasible, including non-tree edges. Encountering
                # their event first requires an odd-cycle contraction.
                for vertex in visited
                    label[vertex]==1 || continue
                    for edge in g.incident[vertex]
                        other=_other(g,edge,vertex)
                        label[other]==1 || continue
                        Int128(_slack(g,dual,edge))<2Int128(delta) && return nothing
                    end
                end
            end
            for vertex in visited
                dual[vertex]=label[vertex]==1 ? Base.checked_add(dual[vertex],delta) :
                    Base.checked_sub(dual[vertex],delta)
            end
        end
        # The final free vertex closes an alternating path to the free root.
        # Preserve the directed input edge recorded at discovery.
        while endpoint!=0
            predecessor=parent[endpoint]
            @assert predecessor!=0
            edge=parent_edge[endpoint]
            @assert edge!=0
            previous=mate_edge[predecessor]
            mate_edge[predecessor]=edge;mate_edge[endpoint]=edge
            endpoint=previous==0 ? 0 : _other(g,previous,predecessor)
        end
    end
    _certificate(g,mate_edge,dual)
    return mate_edge
end

include("OrderedFractionalMatching.jl")

"""
    ordered_perfect_matching(count, edges) -> NamedTuple or nothing

Find an integer minimum-cost perfect matching. Return `nothing` if one cannot
cover every vertex. `mate` contains reciprocal vertex identities and `selected`
contains directed input-edge indices, visited by their original second endpoint.
Bipartite graphs follow the independently tested input-order alternating-tree
policy. General graphs first run certified fractional singleton searches.
Integral completions retain that input-order policy; remaining half-cycles
seed the BSD Blossom kernel with their directed rings, rounded mates, and
certified singleton duals instead of restarting from an empty forest.
"""
function ordered_perfect_matching(count::Int,input)
    graph=_graph(count,input)
    isodd(count) && return nothing
    count==0 && return (mate=Int[],selected=Int[])
    bipartite=_bipartite(graph)
    mate_edge=nothing
    if bipartite
        mate_edge=_bipartite_perfect(graph)
        mate_edge===nothing && return nothing
    else
        state=_fractional_singleton_state(graph)
        if state!==nothing && !any(!iszero,state.cycle_owner)
            _certificate(graph,state.mate,state.dual)
            mate_edge=state.mate
        end
        if mate_edge===nothing
            mate=state===nothing ? minimum_weight_matching(count,graph.edges) :
                minimum_weight_matching(count,graph.edges,
                    _fractional_singleton_seed(graph,state))
            any(iszero,mate) && return nothing
            cheapest=Dict{Tuple{Int,Int},Int}()
            for (index,(a,b,cost)) in enumerate(graph.edges)
                a==b && continue
                key=minmax(a,b);previous=get(cheapest,key,0)
                if previous==0 || cost<=graph.edges[previous][3]
                    cheapest[key]=index
                end
            end
            mate_edge=Int[cheapest[minmax(vertex,mate[vertex])] for vertex in 1:count]
        end
    end
    selected=Int[];sizehint!(selected,count÷2)
    mate=Vector{Int}(undef,count)
    for vertex in 1:count
        edge=mate_edge[vertex]
        mate[vertex]=_other(graph,edge,vertex)
        vertex==graph.edges[edge][2] && push!(selected,edge)
    end
    @assert length(selected)==count÷2
    return (mate=mate,selected=selected)
end

end # module OrderedMatching
