# Independent degree-one fractional matching with singleton duals. Extreme
# solutions contain matched edges and disjoint odd cycles of half-edges.
# A completed integral solution is certified by the same feasible singleton
# duals. This implementation does not translate the restricted native solver.
struct _FractionalSingletonState
    graph::_OrderedGraph
    dual::Vector{Int64}
    mate::Vector{Int}
    cycle_owner::Vector{Int}
    cycles::Vector{Tuple{Vector{Int},Vector{Int}}}
    label::Vector{Int8}
    parent::Vector{Int}
    incoming::Vector{Int}
    children::Vector{Vector{Int}}
    stack::Vector{Int}
    visited::Vector{Int}
end

function _fractional_singleton_augment!(s,endpoint)
    while endpoint!=0
        predecessor=s.parent[endpoint]
        @assert predecessor!=0
        edge=s.incoming[endpoint]
        @assert edge!=0
        previous=s.mate[predecessor]
        s.mate[endpoint]=edge;s.mate[predecessor]=edge
        endpoint=previous==0 ? 0 : _other(s.graph,previous,predecessor)
    end
    return nothing
end

function _fractional_singleton_open_cycle!(s,vertex)
    id=s.cycle_owner[vertex]
    vertices,edges=s.cycles[id]
    at=only(findall(==(vertex),vertices));count=length(vertices)
    @assert isodd(count)
    for node in vertices
        s.cycle_owner[node]=0;s.mate[node]=0
    end
    # Removing the entering vertex leaves an even path with one complete
    # alternating matching. The incoming tree then covers the entering vertex.
    for step in 1:2:count-2
        edge=edges[mod1(at+step,count)]
        first,second,_=s.graph.edges[edge]
        s.mate[first]=edge;s.mate[second]=edge
    end
    @assert s.mate[vertex]==0 && all(node->node==vertex || s.mate[node]!=0,vertices)
    return nothing
end

function _fractional_singleton_close_cycle!(s,first,second,joining_edge)
    first_path=Int[first];second_path=Int[second]
    while s.parent[last(first_path)]!=0
        push!(first_path,s.parent[last(first_path)])
    end
    while s.parent[last(second_path)]!=0
        push!(second_path,s.parent[last(second_path)])
    end
    second_set=Set(second_path)
    common=findfirst(in(second_set),first_path)
    @assert common!==nothing
    base=first_path[common]
    resize!(first_path,common)
    resize!(second_path,only(findall(==(base),second_path)))
    vertices=vcat(reverse(first_path),second_path[1:end-1])
    edges=vcat(reverse(s.incoming[first_path[1:end-1]]),[joining_edge],
        s.incoming[second_path[1:end-1]])
    if s.graph.edges[joining_edge][1]==first
        # The observed construction publishes the ring against the stored
        # direction of the closing edge: keep the base first, then the second
        # tree side, and traverse the joining edge second-to-first.
        vertices=vcat([base],reverse(vertices[2:end]));reverse!(edges)
    end
    @assert length(vertices)==length(edges) && isodd(length(vertices)) && allunique(vertices)
    previous=s.mate[base]
    prefix=previous==0 ? 0 : _other(s.graph,previous,base)
    push!(s.cycles,(vertices,edges));id=length(s.cycles)
    for vertex in vertices
        @assert s.cycle_owner[vertex]==0
        s.cycle_owner[vertex]=id;s.mate[vertex]=0
    end
    # The half-cycle covers its base; reverse the even stem to cover the free
    # root and the remaining stem vertices, releasing the old base mate.
    if prefix!=0
        s.mate[prefix]=0
        _fractional_singleton_augment!(s,prefix)
    end
    return nothing
end

function _fractional_singleton_certificate(s)
    g=s.graph;count=length(s.mate)
    @assert all(vertex->(s.mate[vertex]!=0)!=(s.cycle_owner[vertex]!=0),1:count)
    for edge in eachindex(g.edges)
        first_vertex,second_vertex,_=g.edges[edge]
        first_vertex==second_vertex || @assert _slack(g,s.dual,edge)>=0
    end
    twice_primal=Int128(0)
    for vertex in 1:count
        edge=s.mate[vertex];edge==0 && continue
        other=_other(g,edge,vertex)
        @assert other!=vertex && s.mate[other]==edge && _slack(g,s.dual,edge)==0
        vertex==g.edges[edge][2] && (twice_primal+=2Int128(g.edges[edge][3]))
    end
    for (id,(vertices,edges)) in enumerate(s.cycles)
        s.cycle_owner[first(vertices)]==id || continue
        @assert length(vertices)==length(edges)>=3 && isodd(length(vertices)) && allunique(vertices)
        @assert all(vertex->s.cycle_owner[vertex]==id,vertices)
        for index in eachindex(edges)
            edge=edges[index];vertex=vertices[index]
            @assert _other(g,edge,vertex)==vertices[mod1(index+1,length(vertices))]
            @assert _slack(g,s.dual,edge)==0
            twice_primal+=g.edges[edge][3]
        end
    end
    # Every covered vertex contributes its dual exactly once, including the
    # two half-edges at each odd-cycle vertex. Weak duality proves optimality.
    @assert sum(Int128,s.dual)==twice_primal
    return nothing
end

function _fractional_singleton_state(g::_OrderedGraph)
    count=length(g.incident)
    any(isempty,g.incident) && return nothing
    dual=Int64[minimum(g.edges[edge][3] for edge in g.incident[vertex]) for vertex in 1:count]
    s=_FractionalSingletonState(g,dual,zeros(Int,count),zeros(Int,count),
        Tuple{Vector{Int},Vector{Int}}[],zeros(Int8,count),zeros(Int,count),zeros(Int,count),
        [Int[] for _ in 1:count],Int[],Int[])
    sizehint!(s.stack,count);sizehint!(s.visited,count)
    for vertex in 1:count
        s.mate[vertex]==0 || continue
        delta=typemax(Int64)
        for edge in Iterators.reverse(g.incident[vertex])
            slack=_slack(g,dual,edge)
            @assert slack>=0
            delta=min(delta,slack)
            other=_other(g,edge,vertex)
            if slack==0 && s.mate[vertex]==0 && s.mate[other]==0
                s.mate[vertex]=edge;s.mate[other]=edge
            end
        end
        dual[vertex]=Base.checked_add(dual[vertex],delta)
    end
    for root in 1:count
        (s.mate[root]==0 && s.cycle_owner[root]==0) || continue
        for vertex in s.visited
            s.label[vertex]=0;s.parent[vertex]=0;s.incoming[vertex]=0
            empty!(s.children[vertex])
        end
        empty!(s.visited);s.label[root]=1;push!(s.visited,root)
        completed=false
        while !completed
            empty!(s.stack);push!(s.stack,root)
            while !isempty(s.stack) && !completed
                vertex=pop!(s.stack)
                if s.label[vertex]==1
                    for edge in Iterators.reverse(g.incident[vertex])
                        _slack(g,dual,edge)==0 || continue
                        other=_other(g,edge,vertex)
                        s.label[other]==2 && continue
                        if s.label[other]==1
                            _fractional_singleton_close_cycle!(s,vertex,other,edge)
                            completed=true;break
                        end
                        s.parent[other]=vertex;s.incoming[other]=edge;s.label[other]=2
                        push!(s.visited,other);push!(s.children[vertex],other)
                        if s.cycle_owner[other]!=0
                            _fractional_singleton_open_cycle!(s,other)
                            _fractional_singleton_augment!(s,other)
                            completed=true;break
                        elseif s.mate[other]==0
                            _fractional_singleton_augment!(s,other)
                            completed=true;break
                        end
                        partner=_other(g,s.mate[other],other)
                        @assert s.label[partner]==0
                        s.parent[partner]=other;s.incoming[partner]=s.mate[other];s.label[partner]=1
                        push!(s.visited,partner);push!(s.children[other],partner)
                    end
                end
                append!(s.stack,s.children[vertex])
            end
            completed && break
            delta=typemax(Int64)
            for vertex in s.visited
                s.label[vertex]==1 || continue
                for edge in g.incident[vertex]
                    other=_other(g,edge,vertex);slack=_slack(g,dual,edge)
                    @assert slack>=0
                    if s.label[other]==0
                        delta=min(delta,slack)
                    elseif s.label[other]==1
                        # Tight tree paths give the outer vertices equal dual
                        # parity, so their doubled cost slack is even.
                        @assert iseven(slack)
                        delta=min(delta,slack÷2)
                    end
                end
            end
            delta==typemax(Int64) && return nothing
            @assert delta>0
            for vertex in s.visited
                dual[vertex]=s.label[vertex]==1 ? Base.checked_add(dual[vertex],delta) :
                    Base.checked_sub(dual[vertex],delta)
            end
        end
    end
    _fractional_singleton_certificate(s)
    return s
end

function _fractional_singleton_perfect(g::_OrderedGraph)
    s=_fractional_singleton_state(g)
    s===nothing && return nothing
    any(!iszero,s.cycle_owner) && return nothing
    _certificate(g,s.mate,s.dual)
    return s.mate
end

# Round each retained half-cycle onto the alternating parity its directed
# ring stores: the vertex just before the smallest member in the published
# orientation stays free and the remaining positions pair through the stored
# ring edges. The weighted kernel consumes the same ring in the opposite
# traversal order, so the seed hands it the anchored reversal. The kernel
# seed carries the resulting vertex mates, the certified singleton duals,
# and the consumption-direction rings.
function _fractional_singleton_seed(g::_OrderedGraph,s)
    mate=Int[s.mate[v]==0 ? 0 : _other(g,s.mate[v],v) for v in eachindex(s.mate)]
    cycles=Vector{Int}[]
    for (id,(vertices,edges)) in enumerate(s.cycles)
        s.cycle_owner[first(vertices)]==id || continue
        count=length(vertices);free=mod1(argmin(vertices)-1,count)
        for step in 1:2:count-2
            a,b,_=g.edges[edges[mod1(free+step,count)]]
            mate[a]=b;mate[b]=a
        end
        push!(cycles,vcat(vertices[1],reverse(vertices[2:end])))
    end
    return _SeededInitial(mate,s.dual,cycles)
end
