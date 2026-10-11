# Polynomial primal-dual weighted blossom, integer weights. The algorithm is
# Edmonds/Galil; bookkeeping follows the independently published NetworkX
# reference implementation (BSD-3-Clause, NetworkX developers).
# https://networkx.org/documentation/stable/_modules/networkx/algorithms/matching.html
# Independently verified against exhaustive small graphs and the integer
# NetworkX oracle; production uses only this Julia implementation.
mutable struct _NativeWeightedBlossom
    children::Vector{Int}
    edges::Vector{Tuple{Int,Int}}
    best_edges::Union{Nothing,Vector{Tuple{Int,Int}}}
end
mutable struct _NativeWeightedMatching
    count::Int
    neighbors::Vector{Vector{Int}}
    weight::Dict{Tuple{Int,Int},Int}
    mate::Vector{Int}
    inside::Vector{Int}
    parent::Vector{Int}
    base::Vector{Int}
    label::Vector{Int}
    label_edge::Vector{Tuple{Int,Int}}
    best::Vector{Tuple{Int,Int}}
    vertex_dual::Vector{Int}
    blossom_dual::Vector{Int}
    blossoms::Vector{_NativeWeightedBlossom}
    active::BitVector
    allowed::Set{Tuple{Int,Int}}
    queue::Vector{Int}
end
@inline _native_weighted_slack(s,v,w)=
    s.vertex_dual[v]+s.vertex_dual[w]-2*s.weight[minmax(v,w)]
@inline _native_weighted_at(v,i)=v[mod(i,length(v))+1]
function _native_weighted_leaves!(out,s,b)
    if b<=s.count
        push!(out,b)
    else
        for child in s.blossoms[b].children
            _native_weighted_leaves!(out,s,child)
        end
    end
    return out
end
_native_weighted_leaves(s,b)=_native_weighted_leaves!(Int[],s,b)
function _native_weighted_label!(s,w,label,v)
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
            append!(s.queue,_native_weighted_leaves(s,b))
        end
    elseif label==2
        base=s.base[b]
        _native_weighted_label!(s,s.mate[base],1,base)
    end
    return nothing
end
function _native_weighted_scan!(s,v,w)
    path=Int[];base=0
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
function _native_weighted_add!(s,base,v,w)
    bb=s.inside[base];bv=s.inside[v];bw=s.inside[w]
    b=length(s.base)+1
    push!(s.base,base);push!(s.parent,0);push!(s.label,1)
    push!(s.label_edge,s.label_edge[bb]);push!(s.best,(0,0))
    push!(s.blossom_dual,0);push!(s.active,true)
    path=Int[];edges=Tuple{Int,Int}[(v,w)]
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
    push!(s.blossoms,_NativeWeightedBlossom(path,edges,nothing))
    for vertex in _native_weighted_leaves(s,b)
        s.label[s.inside[vertex]]==2 && push!(s.queue,vertex)
        s.inside[vertex]=b
    end
    best_to=Dict{Int,Tuple{Int,Int}}();order=Int[]
    for child in path
        list=Tuple{Int,Int}[]
        if child>s.count && s.blossoms[child].best_edges!==nothing
            list=s.blossoms[child].best_edges
            s.blossoms[child].best_edges=nothing
        else
            for i in _native_weighted_leaves(s,child),j in s.neighbors[i]
                push!(list,(i,j))
            end
        end
        for edge in list
            i,j=edge
            s.inside[j]==b && ((i,j)=(j,i))
            bj=s.inside[j]
            if bj!=b && s.label[bj]==1
                old=get(best_to,bj,(0,0))
                if old==(0,0) || _native_weighted_slack(s,i,j)<_native_weighted_slack(s,old...)
                    old==(0,0) && push!(order,bj)
                    best_to[bj]=edge
                end
            end
        end
        s.best[child]=(0,0)
    end
    s.blossoms[b].best_edges=[best_to[j] for j in order]
    for edge in s.blossoms[b].best_edges
        if s.best[b]==(0,0) || _native_weighted_slack(s,edge...)<_native_weighted_slack(s,s.best[b]...)
            s.best[b]=edge
        end
    end
    return nothing
end
function _native_weighted_expand!(s,b,endstage)
    children=s.blossoms[b].children;edges=s.blossoms[b].edges
    for child in children
        s.parent[child]=0
        if child>s.count && endstage && s.blossom_dual[child]==0
            _native_weighted_expand!(s,child,endstage)
        else
            for vertex in _native_weighted_leaves(s,child)
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
            p,q=step==1 ? _native_weighted_at(edges,j) : reverse(_native_weighted_at(edges,j-1))
            s.label[w]=0;s.label[q]=0
            _native_weighted_label!(s,w,2,v)
            push!(s.allowed,minmax(p,q));j+=step
            v,w=step==1 ? _native_weighted_at(edges,j) : reverse(_native_weighted_at(edges,j-1))
            push!(s.allowed,minmax(v,w));j+=step
        end
        bw=_native_weighted_at(children,j)
        s.label[w]=2;s.label[bw]=2
        s.label_edge[w]=(v,w);s.label_edge[bw]=(v,w);s.best[bw]=(0,0)
        j+=step
        while _native_weighted_at(children,j)!=entry
            bv=_native_weighted_at(children,j)
            if s.label[bv]==1
                j+=step;continue
            end
            vertex=bv
            if bv>s.count
                reachable=findfirst(v->s.label[v]!=0,_native_weighted_leaves(s,bv))
                if reachable===nothing
                    j+=step;continue
                end
                vertex=_native_weighted_leaves(s,bv)[reachable]
            end
            if s.label[vertex]!=0
                @assert s.label[vertex]==2 && s.inside[vertex]==bv
                s.label[vertex]=0;s.label[s.mate[s.base[bv]]]=0
                _native_weighted_label!(s,vertex,2,s.label_edge[vertex][1])
            end
            j+=step
        end
    end
    s.label[b]=0;s.label_edge[b]=(0,0);s.best[b]=(0,0);s.active[b]=false
    return nothing
end
function _native_weighted_augment_blossom!(s,b,vertex)
    child=vertex
    while s.parent[child]!=b
        child=s.parent[child]
    end
    child>s.count && _native_weighted_augment_blossom!(s,child,vertex)
    children=s.blossoms[b].children;edges=s.blossoms[b].edges
    i=findfirst(==(child),children)-1;j=i
    if isodd(i)
        j-=length(children);step=1
    else
        step=-1
    end
    while j!=0
        j+=step;child=_native_weighted_at(children,j)
        w,x=step==1 ? _native_weighted_at(edges,j) : reverse(_native_weighted_at(edges,j-1))
        child>s.count && _native_weighted_augment_blossom!(s,child,w)
        j+=step;child=_native_weighted_at(children,j)
        child>s.count && _native_weighted_augment_blossom!(s,child,x)
        s.mate[w]=x;s.mate[x]=w
    end
    s.blossoms[b].children=vcat(children[i+1:end],children[1:i])
    s.blossoms[b].edges=vcat(edges[i+1:end],edges[1:i])
    s.base[b]=s.base[first(s.blossoms[b].children)]
    @assert s.base[b]==vertex
    return nothing
end
function _native_weighted_augment!(s,v,w)
    for (vertex,partner) in ((v,w),(w,v))
        while true
            b=s.inside[vertex]
            @assert s.label[b]==1
            b>s.count && _native_weighted_augment_blossom!(s,b,vertex)
            s.mate[vertex]=partner
            s.label_edge[b]==(0,0) && break
            t=s.label_edge[b][1];bt=s.inside[t]
            @assert s.label[bt]==2
            vertex,partner=s.label_edge[bt]
            bt>s.count && _native_weighted_augment_blossom!(s,bt,partner)
            s.mate[partner]=vertex
        end
    end
    return nothing
end
function _native_weighted_verify(s)
    offset=max(0,-minimum(s.vertex_dual))
    for ((i,j),weight) in s.weight
        slack=s.vertex_dual[i]+s.vertex_dual[j]-2*weight
        ib=Int[i];jb=Int[j]
        while s.parent[last(ib)]!=0
            push!(ib,s.parent[last(ib)])
        end
        while s.parent[last(jb)]!=0
            push!(jb,s.parent[last(jb)])
        end
        reverse!(ib);reverse!(jb)
        for (a,b) in zip(ib,jb)
            a==b || break
            slack+=2*s.blossom_dual[a]
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
            edges=s.blossoms[b].edges
            @assert isodd(length(edges))
            for (i,j) in edges[2:2:end]
                @assert s.mate[i]==j && s.mate[j]==i
            end
        end
    end
    return nothing
end
function _native_weighted_perfect_matching(count,input)
    count==0 && return Tuple{Int,Int}[]
    neighbors=[Int[] for _ in 1:count]
    weights=Dict{Tuple{Int,Int},Int}()
    maximum_cost=maximum((cost for (_,_,cost) in input);init=0)
    for (v,w,cost) in input
        v==w && continue
        edge=minmax(v,w);value=maximum_cost+1-cost
        if !haskey(weights,edge)
            push!(neighbors[v],w);push!(neighbors[w],v)
            weights[edge]=value
        else
            weights[edge]=max(weights[edge],value)
        end
    end
    maximum_weight=maximum(values(weights);init=0)
    s=_NativeWeightedMatching(count,neighbors,weights,zeros(Int,count),collect(1:count),
        zeros(Int,count),collect(1:count),zeros(Int,count),fill((0,0),count),
        fill((0,0),count),fill(maximum_weight,count),zeros(Int,count),
        [_NativeWeightedBlossom(Int[],Tuple{Int,Int}[],nothing) for _ in 1:count],
        trues(count),Set{Tuple{Int,Int}}(),Int[])
    while true
        fill!(s.label,0);fill!(s.label_edge,(0,0));fill!(s.best,(0,0))
        for b in count+1:length(s.base)
            s.active[b] && (s.blossoms[b].best_edges=nothing)
        end
        empty!(s.allowed);empty!(s.queue)
        for v in 1:count
            s.mate[v]==0 && s.label[s.inside[v]]==0 && _native_weighted_label!(s,v,1,0)
        end
        augmented=false
        while true
            while !isempty(s.queue) && !augmented
                v=pop!(s.queue)
                @assert s.label[s.inside[v]]==1
                for w in s.neighbors[v]
                    bv=s.inside[v];bw=s.inside[w];bv==bw && continue
                    edge=minmax(v,w)
                    slack=_native_weighted_slack(s,v,w)
                    slack<=0 && push!(s.allowed,edge)
                    if edge in s.allowed
                        if s.label[bw]==0
                            _native_weighted_label!(s,w,2,v)
                        elseif s.label[bw]==1
                            base=_native_weighted_scan!(s,v,w)
                            if base!=0
                                _native_weighted_add!(s,base,v,w)
                            else
                                _native_weighted_augment!(s,v,w);augmented=true;break
                            end
                        elseif s.label[w]==0
                            s.label[w]=2;s.label_edge[w]=(v,w)
                        end
                    elseif s.label[bw]==1
                        (s.best[bv]==(0,0) || slack<_native_weighted_slack(s,s.best[bv]...)) && (s.best[bv]=(v,w))
                    elseif s.label[w]==0
                        (s.best[w]==(0,0) || slack<_native_weighted_slack(s,s.best[w]...)) && (s.best[w]=(v,w))
                    end
                end
            end
            augmented && break
            kind=0;delta=0;edge=(0,0);blossom=0
            for v in 1:count
                if s.label[s.inside[v]]==0 && s.best[v]!=(0,0)
                    value=_native_weighted_slack(s,s.best[v]...)
                    if kind==0 || value<delta
                        kind=2;delta=value;edge=s.best[v]
                    end
                end
            end
            for b in eachindex(s.base)
                if s.active[b] && s.parent[b]==0 && s.label[b]==1 && s.best[b]!=(0,0)
                    slack=_native_weighted_slack(s,s.best[b]...)
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
                    s.vertex_dual[v]-=delta
                elseif s.label[s.inside[v]]==2
                    s.vertex_dual[v]+=delta
                end
            end
            for b in count+1:length(s.base)
                s.active[b] && s.parent[b]==0 || continue
                if s.label[b]==1
                    s.blossom_dual[b]+=delta
                elseif s.label[b]==2
                    s.blossom_dual[b]-=delta
                end
            end
            kind==1 && break
            if kind==2 || kind==3
                push!(s.allowed,minmax(edge...));push!(s.queue,edge[1])
            else
                _native_weighted_expand!(s,blossom,false)
            end
        end
        for v in 1:count
            s.mate[v]!=0 && @assert s.mate[s.mate[v]]==v
        end
        augmented || break
        for b in count+1:length(s.base)
            s.active[b] && s.parent[b]==0 && s.label[b]==1 && s.blossom_dual[b]==0 &&
                _native_weighted_expand!(s,b,true)
        end
    end
    _native_weighted_verify(s)
    any(==(0),s.mate) && return nothing
    return [(v,s.mate[v]) for v in 1:count if s.mate[v]>v]
end
