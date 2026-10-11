# Integer primal-dual weighted Blossom matching, based on Edmonds/Galil.
# Adapted from NetworkX 3.6.1 (NetworkX Developers, BSD-3-Clause).
# Full license and pinned source provenance: LICENSES/networkx-BSD-3-Clause.txt
# and LICENSES/weighted-matching-provenance.md. No foreign production backend.
module WeightedMatching

export minimum_weight_matching

mutable struct _WeightedBlossom
    children::Vector{Int}
    edges::Vector{Tuple{Int,Int}}
    best_edges::Vector{Tuple{Int,Int}}
    best_valid::Bool
end
mutable struct _WeightedMatching
    count::Int
    neighbors::Vector{Vector{Int}}
    weight::Vector{Int64}
    edges::Vector{Tuple{Int,Int}}
    edge_index::Dict{Tuple{Int,Int},Int}
    mate::Vector{Int}
    inside::Vector{Int}
    parent::Vector{Int}
    base::Vector{Int}
    label::Vector{Int}
    label_edge::Vector{Tuple{Int,Int}}
    best::Vector{Tuple{Int,Int}}
    vertex_dual::Vector{Int64}
    blossom_dual::Vector{Int64}
    blossoms::Vector{Union{Nothing,_WeightedBlossom}}
    active::BitVector
    allowed::BitVector
    queue::Vector{Int}
    free::Vector{Int}
    leaves::Vector{Int}
    leaf_stack::Vector{Int}
    scan_path::Vector{Int}
    ancestors_first::Vector{Int}
    ancestors_second::Vector{Int}
    best_to::Vector{Tuple{Int,Int}}
    best_order::Vector{Int}
end
@inline _weighted_slack(s,v,w)=
    Base.checked_sub(Base.checked_add(s.vertex_dual[v],s.vertex_dual[w]),
        Base.checked_mul(Int64(2),s.weight[s.edge_index[minmax(v,w)]]))
@inline _weighted_allow!(s,v,w)=(s.allowed[s.edge_index[minmax(v,w)]]=true)
@inline _weighted_at(v,i)=v[mod(i,length(v))+1]
@inline _weighted_blossom(s,b)=s.blossoms[b]::_WeightedBlossom
function _weighted_leaves(s,b)
    out=s.leaves;stack=s.leaf_stack
    empty!(out);empty!(stack);push!(stack,b)
    while !isempty(stack)
        child=pop!(stack)
        if child<=s.count
            push!(out,child)
        else
            append!(stack,_weighted_blossom(s,child).children)
        end
    end
    return out
end
function _weighted_label!(s,w,label,v)
    b=s.inside[w]
    @assert s.label[w]==0 && s.label[b]==0
    s.label[w]=label;s.label[b]=label
    edge=v==0 ? (0,0) : (v,w)
    s.label_edge[w]=edge;s.label_edge[b]=edge
    s.best[w]=(0,0);s.best[b]=(0,0)
    if label==1
        if b<=s.count
            push!(s.queue,b)
        else
            append!(s.queue,_weighted_leaves(s,b))
        end
    elseif label==2
        base=s.base[b]
        _weighted_label!(s,s.mate[base],1,base)
    end
    return nothing
end
function _weighted_scan!(s,v,w)
    path=s.scan_path;empty!(path);base=0
    while v!=0
        b=s.inside[v]
        if s.label[b]&4!=0
            base=s.base[b];break
        end
        @assert s.label[b]==1
        push!(path,b);s.label[b]=5
        if s.label_edge[b]==(0,0)
            @assert s.mate[s.base[b]]==0
            v=0
        else
            @assert s.label_edge[b][1]==s.mate[s.base[b]]
            v=s.label_edge[b][1];b=s.inside[v]
            @assert s.label[b]==2
            v=s.label_edge[b][1]
        end
        w!=0 && ((v,w)=(w,v))
    end
    for b in path
        s.label[b]=1
    end
    return base
end
# Rebuild a contracted group's retained per-target best-edge cache from its
# children's valid lists or a fresh leaf scan, then its own minimum entry.
function _weighted_collect_best!(s,b)
    storage=_weighted_blossom(s,b);order=s.best_order
    for target in order;s.best_to[target]=(0,0);end
    empty!(order)
    for child in storage.children
        if child>s.count && _weighted_blossom(s,child).best_valid
            child_storage=_weighted_blossom(s,child)
            for edge in child_storage.best_edges
                _weighted_consider_best!(s,b,edge)
            end
            child_storage.best_valid=false
        else
            for i in _weighted_leaves(s,child),j in s.neighbors[i]
                _weighted_consider_best!(s,b,(i,j))
            end
        end
        s.best[child]=(0,0)
    end
    empty!(storage.best_edges)
    for target in order;push!(storage.best_edges,s.best_to[target]);end
    storage.best_valid=true
    for edge in storage.best_edges
        if s.best[b]==(0,0) || _weighted_slack(s,edge...)<_weighted_slack(s,s.best[b]...)
            s.best[b]=edge
        end
    end
    return nothing
end
function _weighted_add!(s,base,v,w)
    bb=s.inside[base];bv=s.inside[v];bw=s.inside[w]
    isempty(s.free) && throw(ErrorException("weighted matching: blossom workspace exhausted"))
    b=pop!(s.free)
    s.base[b]=base;s.parent[b]=0;s.label[b]=1
    s.label_edge[b]=s.label_edge[bb];s.best[b]=(0,0)
    s.blossom_dual[b]=0;s.active[b]=true
    if s.blossoms[b]===nothing
        s.blossoms[b]=_WeightedBlossom(Int[],Tuple{Int,Int}[],Tuple{Int,Int}[],false)
    end
    storage=_weighted_blossom(s,b)
    path=storage.children;edges=storage.edges
    empty!(path);empty!(edges);push!(edges,(v,w));storage.best_valid=false
    s.parent[bb]=b
    while bv!=bb
        s.parent[bv]=b;push!(path,bv);push!(edges,s.label_edge[bv])
        v=s.label_edge[bv][1];bv=s.inside[v]
    end
    push!(path,bb);reverse!(path);reverse!(edges)
    while bw!=bb
        s.parent[bw]=b;push!(path,bw);push!(edges,reverse(s.label_edge[bw]))
        w=s.label_edge[bw][1];bw=s.inside[w]
    end
    for vertex in _weighted_leaves(s,b)
        s.label[s.inside[vertex]]==2 && push!(s.queue,vertex)
        s.inside[vertex]=b
    end
    _weighted_collect_best!(s,b)
    return nothing
end
function _weighted_consider_best!(s,b,edge)
    i,j=edge
    s.inside[j]==b && ((i,j)=(j,i))
    target=s.inside[j]
    if target!=b && s.label[target]==1
        old=s.best_to[target]
        if old==(0,0) || _weighted_slack(s,i,j)<_weighted_slack(s,old...)
            old==(0,0) && push!(s.best_order,target)
            s.best_to[target]=edge
        end
    end
    return nothing
end
function _weighted_expand!(s,b,endstage)
    storage=_weighted_blossom(s,b);children=storage.children;edges=storage.edges
    for child in children
        s.parent[child]=0
        if child>s.count && endstage && s.blossom_dual[child]==0
            _weighted_expand!(s,child,endstage)
        else
            for vertex in _weighted_leaves(s,child)
                s.inside[vertex]=child
            end
        end
    end
    if !endstage && s.label[b]==2
        entry=s.inside[s.label_edge[b][2]]
        j=findfirst(==(entry),children)-1
        if isodd(j)
            j-=length(children);step=1
        else
            step=-1
        end
        v,w=s.label_edge[b]
        while j!=0
            p,q=step==1 ? _weighted_at(edges,j) : reverse(_weighted_at(edges,j-1))
            s.label[w]=0;s.label[q]=0
            _weighted_label!(s,w,2,v)
            _weighted_allow!(s,p,q);j+=step
            v,w=step==1 ? _weighted_at(edges,j) : reverse(_weighted_at(edges,j-1))
            _weighted_allow!(s,v,w);j+=step
        end
        bw=_weighted_at(children,j)
        s.label[w]=2;s.label[bw]=2
        s.label_edge[w]=(v,w);s.label_edge[bw]=(v,w);s.best[bw]=(0,0)
        j+=step
        while _weighted_at(children,j)!=entry
            bv=_weighted_at(children,j)
            if s.label[bv]==1
                j+=step;continue
            end
            vertex=bv
            if bv>s.count
                vertex=0
                for leaf in _weighted_leaves(s,bv)
                    if s.label[leaf]!=0
                        vertex=leaf;break
                    end
                end
                if vertex==0
                    j+=step;continue
                end
            end
            if s.label[vertex]!=0
                @assert s.label[vertex]==2 && s.inside[vertex]==bv
                s.label[vertex]=0;s.label[s.mate[s.base[bv]]]=0
                _weighted_label!(s,vertex,2,s.label_edge[vertex][1])
            end
            j+=step
        end
    end
    s.label[b]=0;s.label_edge[b]=(0,0);s.best[b]=(0,0);s.active[b]=false
    s.parent[b]=0;s.base[b]=0;s.blossom_dual[b]=0
    empty!(children);empty!(edges);empty!(storage.best_edges);storage.best_valid=false
    push!(s.free,b)
    return nothing
end
function _weighted_augment_blossom!(s,b,vertex)
    child=vertex
    while s.parent[child]!=b
        child=s.parent[child]
    end
    child>s.count && _weighted_augment_blossom!(s,child,vertex)
    children=_weighted_blossom(s,b).children;edges=_weighted_blossom(s,b).edges
    i=findfirst(==(child),children)-1;j=i
    if isodd(i)
        j-=length(children);step=1
    else
        step=-1
    end
    while j!=0
        j+=step;child=_weighted_at(children,j)
        w,x=step==1 ? _weighted_at(edges,j) : reverse(_weighted_at(edges,j-1))
        child>s.count && _weighted_augment_blossom!(s,child,w)
        j+=step;child=_weighted_at(children,j)
        child>s.count && _weighted_augment_blossom!(s,child,x)
        s.mate[w]=x;s.mate[x]=w
    end
    if i>0
        reverse!(children,1,i);reverse!(children,i+1,length(children));reverse!(children)
        reverse!(edges,1,i);reverse!(edges,i+1,length(edges));reverse!(edges)
    end
    s.base[b]=s.base[first(children)]
    @assert s.base[b]==vertex
    return nothing
end
function _weighted_augment!(s,v,w)
    for (vertex,partner) in ((v,w),(w,v))
        while true
            b=s.inside[vertex]
            @assert s.label[b]==1
            b>s.count && _weighted_augment_blossom!(s,b,vertex)
            s.mate[vertex]=partner
            s.label_edge[b]==(0,0) && break
            t=s.label_edge[b][1];bt=s.inside[t]
            @assert s.label[bt]==2
            vertex,partner=s.label_edge[bt]
            bt>s.count && _weighted_augment_blossom!(s,bt,partner)
            s.mate[partner]=vertex
        end
    end
    return nothing
end
function _weighted_verify(s)
    offset=max(0,-minimum(s.vertex_dual))
    for (i,j) in s.edges
        slack=_weighted_slack(s,i,j)
        ib=s.ancestors_first;jb=s.ancestors_second
        empty!(ib);empty!(jb);push!(ib,i);push!(jb,j)
        while s.parent[last(ib)]!=0
            push!(ib,s.parent[last(ib)])
        end
        while s.parent[last(jb)]!=0
            push!(jb,s.parent[last(jb)])
        end
        reverse!(ib);reverse!(jb)
        for (a,b) in zip(ib,jb)
            a==b || break
            slack=Base.checked_add(slack,Base.checked_mul(Int64(2),s.blossom_dual[a]))
        end
        @assert slack>=0
        s.mate[i]==j && @assert s.mate[j]==i && slack==0
    end
    for v in 1:s.count
        @assert s.vertex_dual[v]+offset>=0
        s.mate[v]==0 && @assert s.vertex_dual[v]+offset==0
    end
    for b in s.count+1:length(s.base)
        s.active[b] || continue
        @assert s.blossom_dual[b]>=0
        if s.blossom_dual[b]>0
            edges=_weighted_blossom(s,b).edges
            @assert isodd(length(edges))
            for index in 2:2:length(edges)
                i,j=edges[index]
                @assert s.mate[i]==j && s.mate[j]==i
            end
        end
    end
    return nothing
end
# Consume the validated directed edge input once. Parallel edges collapse to
# one undirected pair at minimum cost; `directed` additionally retains the
# last minimum-cost directed form of each pair, which the seeded continuation
# uses to publish a contracted ring in the observed orientation.
function _weighted_input(count::Int,input;directed::Bool=false)
    count>=0 || throw(ArgumentError("weighted matching: count must be nonnegative"))
    count<=typemax(Int)÷2 || throw(ArgumentError("weighted matching: workspace count overflows Int"))
    neighbors=[Int[] for _ in 1:count]
    costs=Int64[];edges=Tuple{Int,Int}[];edge_index=Dict{Tuple{Int,Int},Int}()
    directed_edges=directed ? Tuple{Int,Int}[] : nothing
    if Base.IteratorSize(typeof(input)) isa Union{Base.HasLength,Base.HasShape}
        entries=length(input)
        entries<=typemax(Int)÷2 || throw(ArgumentError("weighted matching: edge workspace overflows Int"))
        # Parallel edges collapse into one undirected pair. Bound reservation
        # by that number while still validating every input edge in order.
        possible=Int128(count)*(Int128(count)-1)÷2
        entries=Int(min(Int128(entries),possible))
        sizehint!(costs,entries);sizehint!(edges,entries)
        # Tuple-key hash clustering otherwise causes a second table growth
        # at some shifted input sizes. Reserve a low-load sparse index once.
        sizehint!(edge_index,2entries)
        directed && sizehint!(directed_edges,entries)
    end
    # This conservative range keeps doubled slacks and accumulated blossom
    # duals inside Int64; checked arithmetic remains active in the kernel.
    cost_limit=typemax(Int64)÷(16*max(Int128(1),Int128(count)))
    for (v,w,cost) in input
        v isa Integer && w isa Integer && !(v isa Bool) && !(w isa Bool) ||
            throw(ArgumentError("weighted matching: endpoints must be integer vertices"))
        1<=v<=count && 1<=w<=count || throw(ArgumentError("weighted matching: vertex out of range"))
        cost isa Integer && !(cost isa Bool) || throw(ArgumentError("weighted matching: costs must be integers"))
        -cost_limit<=cost<=cost_limit || throw(ArgumentError("weighted matching: integer cost exceeds checked dual range"))
        v==w && continue
        v=Int(v);w=Int(w);edge=minmax(v,w);value=Int64(cost)
        index=get(edge_index,edge,0)
        if index==0
            push!(neighbors[v],w);push!(neighbors[w],v)
            push!(edges,edge);push!(costs,value);edge_index[edge]=length(edges)
            directed && push!(directed_edges,(v,w))
        else
            if directed && value<=costs[index]
                directed_edges[index]=(v,w)
            end
            costs[index]=min(costs[index],value)
        end
    end
    return neighbors,costs,edges,edge_index,directed_edges
end
"""
    minimum_weight_matching(count, edges) -> Vector{Int}

Return a minimum-cost maximum-cardinality matching on an undirected general
graph with vertices `1:count`. Each edge is `(first, second, integer_cost)`.
Unmatched vertices have mate zero; parallel edges retain their minimum cost,
and self-loops are ignored. Integer primal-dual arithmetic also certifies the
returned optimum. This function does not impose a primary engine's tie rule.
"""
function minimum_weight_matching(count::Int,input)
    neighbors,costs,edges,edge_index,_=_weighted_input(count,input)
    isempty(costs) && return zeros(Int,count)
    maximum_cost=maximum(costs)
    offset=Base.checked_add(maximum_cost,Int64(1))
    for index in eachindex(costs);costs[index]=Base.checked_sub(offset,costs[index]);end
    maximum_weight=maximum(costs)
    slots=2count
    storage=Vector{Union{Nothing,_WeightedBlossom}}(nothing,slots)
    s=_WeightedMatching(count,neighbors,costs,edges,edge_index,zeros(Int,count),collect(1:count),
        zeros(Int,slots),vcat(collect(1:count),zeros(Int,count)),zeros(Int,slots),fill((0,0),slots),
        fill((0,0),slots),fill(maximum_weight,count),zeros(Int64,slots),
        storage,vcat(trues(count),falses(count)),falses(length(edges)),Int[],
        collect(slots:-1:count+1),Int[],Int[],Int[],Int[],Int[],fill((0,0),slots),Int[])
    for workspace in (s.queue,s.leaves,s.leaf_stack,s.scan_path,s.ancestors_first,s.ancestors_second,s.best_order)
        sizehint!(workspace,count)
    end
    while true
        fill!(s.label,0);fill!(s.label_edge,(0,0));fill!(s.best,(0,0))
        for b in count+1:length(s.base)
            s.active[b] && (_weighted_blossom(s,b).best_valid=false)
        end
        fill!(s.allowed,false);empty!(s.queue)
        for v in 1:count
            s.mate[v]==0 && s.label[s.inside[v]]==0 && _weighted_label!(s,v,1,0)
        end
        augmented=false
        while true
            while !isempty(s.queue) && !augmented
                v=pop!(s.queue)
                @assert s.label[s.inside[v]]==1
                for w in s.neighbors[v]
                    bv=s.inside[v];bw=s.inside[w];bv==bw && continue
                    edge_id=s.edge_index[minmax(v,w)]
                    slack=_weighted_slack(s,v,w)
                    slack<=0 && (s.allowed[edge_id]=true)
                    if s.allowed[edge_id]
                        if s.label[bw]==0
                            _weighted_label!(s,w,2,v)
                        elseif s.label[bw]==1
                            base=_weighted_scan!(s,v,w)
                            if base!=0
                                _weighted_add!(s,base,v,w)
                            else
                                _weighted_augment!(s,v,w);augmented=true;break
                            end
                        elseif s.label[w]==0
                            s.label[w]=2;s.label_edge[w]=(v,w)
                        end
                    elseif s.label[bw]==1
                        (s.best[bv]==(0,0) || slack<_weighted_slack(s,s.best[bv]...)) && (s.best[bv]=(v,w))
                    elseif s.label[w]==0
                        (s.best[w]==(0,0) || slack<_weighted_slack(s,s.best[w]...)) && (s.best[w]=(v,w))
                    end
                end
            end
            augmented && break
            kind=0;delta=Int64(0);edge=(0,0);blossom=0
            for v in 1:count
                if s.label[s.inside[v]]==0 && s.best[v]!=(0,0)
                    value=_weighted_slack(s,s.best[v]...)
                    if kind==0 || value<delta
                        kind=2;delta=value;edge=s.best[v]
                    end
                end
            end
            for b in eachindex(s.base)
                if s.active[b] && s.parent[b]==0 && s.label[b]==1 && s.best[b]!=(0,0)
                    slack=_weighted_slack(s,s.best[b]...)
                    @assert iseven(slack)
                    value=slack÷2
                    if kind==0 || value<delta
                        kind=3;delta=value;edge=s.best[b]
                    end
                end
            end
            for b in count+1:length(s.base)
                if s.active[b] && s.parent[b]==0 && s.label[b]==2 &&
                   (kind==0 || s.blossom_dual[b]<delta)
                    kind=4;delta=s.blossom_dual[b];blossom=b
                end
            end
            if kind==0
                kind=1;delta=max(0,minimum(s.vertex_dual))
            end
            for v in 1:count
                if s.label[s.inside[v]]==1
                    s.vertex_dual[v]=Base.checked_sub(s.vertex_dual[v],delta)
                elseif s.label[s.inside[v]]==2
                    s.vertex_dual[v]=Base.checked_add(s.vertex_dual[v],delta)
                end
            end
            for b in count+1:length(s.base)
                s.active[b] && s.parent[b]==0 || continue
                if s.label[b]==1
                    s.blossom_dual[b]=Base.checked_add(s.blossom_dual[b],delta)
                elseif s.label[b]==2
                    s.blossom_dual[b]=Base.checked_sub(s.blossom_dual[b],delta)
                end
            end
            kind==1 && break
            if kind==2 || kind==3
                _weighted_allow!(s,edge...);push!(s.queue,edge[1])
            else
                _weighted_expand!(s,blossom,false)
            end
        end
        for v in 1:count
            s.mate[v]!=0 && @assert s.mate[s.mate[v]]==v
        end
        augmented || break
        for b in count+1:length(s.base)
            s.active[b] && s.parent[b]==0 && s.label[b]==1 && s.blossom_dual[b]==0 &&
                _weighted_expand!(s,b,true)
        end
    end
    _weighted_verify(s)
    return s.mate
end

include("SeededWeightedMatching.jl")

end # module WeightedMatching
