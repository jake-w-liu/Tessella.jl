# Half-integral seeded continuation of the weighted Blossom search, based on
# the same Edmonds/Galil primal-dual bookkeeping as the BSD kernel above.
# When the certified degree-one relaxation retains tight odd cycles, each ring
# contracts as a pre-formed blossom and the same search resumes on that seed
# instead of restarting from an empty forest. The observed native construction
# order is retained on the remaining ties: the forest grows one actual free
# base at a time, each newly labeled matched subtree is scanned depth-first,
# every scalar dual event rebuilds the scan queue in label-tree order with
# leaves rotated to their anchor, and a fresh outer contraction that covers
# the current free base after a positive dual step defers the forest to the
# next actual free root. Retained cycles anchor their leaf order just before
# the smallest leaf; formed blossoms anchor at their free base; retained
# cycles are not expanded at stage end even at zero dual. All dual arithmetic
# and the final primal/dual certificate are the checked kernel's own.

# A certified half-integral partial solution: vertex mates for the matched
# part, the certified feasible singleton duals, and the retained directed
# odd-cycle vertex rings (each ring leaves exactly one vertex unmatched).
struct _SeededInitial
    mate::Vector{Int}
    dual::Vector{Int64}
    cycles::Vector{Vector{Int}}
end

mutable struct _SeededWeightedMatching
    count::Int
    neighbors::Vector{Vector{Int}}
    weight::Vector{Int64}
    edges::Vector{Tuple{Int,Int}}
    directed::Vector{Tuple{Int,Int}}
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
    initial_cycle::BitVector
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
    offset::Int64
    batches::Vector{Vector{Int}}
    batch_depth::Int
    roots::Vector{Int}
    outer_members::Vector{Int}
    outer_children::Vector{Union{Nothing,Vector{Int}}}
    outer_incoming::Vector{Int}
    outer_marks::Vector{Int8}
    outer_order::Vector{Int}
    outer_roots::Vector{Int}
    outer_stack::Vector{Int}
    positive_tree_dual::Bool
    root_rotation_pending::Bool
    next_free_root::Int
    allow_dual::Bool
end

# Forward constructor-order leaf walk; the plain leaf walk yields stack order.
function _seeded_ring_leaves!(s,b)
    out=s.leaves;stack=s.leaf_stack
    empty!(out);empty!(stack);push!(stack,b)
    while !isempty(stack)
        child=pop!(stack)
        if child<=s.count
            push!(out,child)
        else
            append!(stack,Iterators.reverse(_weighted_blossom(s,child).children))
        end
    end
    return out
end

function _seeded_label!(s,w,label,v)
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
            leaves=_seeded_ring_leaves!(s,b)
            # Retained cycles re-enter just before their smallest leaf; a
            # formed blossom re-enters at its actual free base.
            at=s.initial_cycle[b] ? mod1(argmin(leaves)-1,length(leaves)) :
                findfirst(==(s.base[b]),leaves)
            @assert at!==nothing
            for i in at:length(leaves);push!(s.queue,leaves[i]);end
            for i in 1:at-1;push!(s.queue,leaves[i]);end
        end
    elseif label==2
        base=s.base[b]
        _seeded_label!(s,s.mate[base],1,base)
    end
    return nothing
end

function _seeded_add!(s,base,v,w)
    # The retained directed form of the tight edge fixes ring publication
    # orientation: a closing edge stored opposite the scanned direction walks
    # the constructor ring in the opposite sense.
    reverse_constructor=s.directed[s.edge_index[minmax(v,w)]][1]!=v
    bb=s.inside[base];bv=s.inside[v];bw=s.inside[w]
    isempty(s.free) && throw(ErrorException("seeded weighted matching: blossom workspace exhausted"))
    b=pop!(s.free)
    s.base[b]=base;s.parent[b]=0;s.label[b]=1
    s.label_edge[b]=s.label_edge[bb];s.best[b]=(0,0)
    s.blossom_dual[b]=0;s.active[b]=true;s.initial_cycle[b]=false
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
    if reverse_constructor
        reverse!(path,2);reverse!(edges)
        for index in eachindex(edges)
            first,second=edges[index];edges[index]=(second,first)
        end
    end
    for vertex in _weighted_leaves(s,b)
        s.label[s.inside[vertex]]==2 && push!(s.queue,vertex)
        s.inside[vertex]=b
    end
    _weighted_collect_best!(s,b)
    # A new outer contraction after a positive scalar dual event can finish
    # its tight queue before the next actual free root is searched. Existing
    # blossoms do not set this flag when a forest is restarted.
    if s.positive_tree_dual && s.label_edge[b]==(0,0) && s.mate[s.base[b]]==0
        s.root_rotation_pending=true
    end
    return nothing
end

# Rebuild the actual label forest for a rescan: preorder over the top-level
# outer groups, roots ascending, each group's children ascending by the index
# of the edge that labeled their inner link.
function _seeded_tree_order!(s,outer)
    marks=s.outer_marks;children=s.outer_children;incoming=s.outer_incoming
    roots=s.outer_roots;empty!(roots)
    order=s.outer_order;empty!(order)
    for node in outer
        marks[node]=1
        kids=children[node]
        kids===nothing && (children[node]=kids=Int[])
        empty!(kids)
    end
    for node in outer
        edge=s.label_edge[node]
        if edge==(0,0)
            push!(roots,node);continue
        end
        inner=s.inside[edge[1]]
        @assert s.label[inner]==2
        link=s.label_edge[inner]
        @assert link!=(0,0)
        parent=s.inside[link[1]]
        @assert marks[parent]==1 && parent!=node
        incoming[node]=s.edge_index[minmax(link[1],link[2])]
        push!(children[parent]::Vector{Int},node)
    end
    @assert !isempty(roots)
    for node in outer
        sort!(children[node]::Vector{Int};by=child->(incoming[child],child))
    end
    sort!(roots)
    stack=s.outer_stack;empty!(stack)
    for i in length(roots):-1:1;push!(stack,roots[i]);end
    while !isempty(stack)
        node=pop!(stack)
        @assert marks[node]==1
        marks[node]=2
        push!(order,node)
        kids=children[node]::Vector{Int}
        for i in length(kids):-1:1;push!(stack,kids[i]);end
    end
    for node in outer;marks[node]=0;end
    return order
end

# A scalar dual step tightens every attaining edge at once. Revisit the whole
# existing forest in label-tree order instead of prioritizing one cached edge.
function _seeded_rescan_queue!(s)
    empty!(s.queue)
    outer=s.outer_members;empty!(outer)
    for b in eachindex(s.base)
        s.active[b] && s.parent[b]==0 && s.label[b]==1 && push!(outer,b)
    end
    for b in _seeded_tree_order!(s,outer)
        leaves=_seeded_ring_leaves!(s,b)
        at=s.initial_cycle[b] ? mod1(argmin(leaves)-1,length(leaves)) :
            findfirst(==(s.base[b]),leaves)
        at===nothing && error("seeded weighted matching: outer rescan lost a free base")
        for i in at:length(leaves);push!(s.queue,leaves[i]);end
        for i in 1:at-1;push!(s.queue,leaves[i]);end
    end
    return nothing
end

# Scan exactly the entries appended to the queue after `first`, in arrival
# order; the region is consumed whether or not an augment results. Buffers
# are pooled by recursion depth so deeper descents cannot clobber this frame.
function _seeded_scan_batch!(s,first)
    depth=s.batch_depth+=1
    depth>length(s.batches) && push!(s.batches,Int[])
    batch=s.batches[depth];empty!(batch)
    queue=s.queue
    for i in first+1:length(queue);push!(batch,queue[i]);end
    resize!(queue,first)
    found=false
    for vertex in batch
        _seeded_depth_scan!(s,vertex) && (found=true;break)
    end
    s.batch_depth-=1
    return found
end

function _seeded_depth_scan!(s,v)
    @assert s.label[s.inside[v]]==1
    # Prefer an already tight connection to an actual free base before
    # descending any matched-neighbor subtree. This check does not label
    # siblings or change the existing adjacency order and dual arithmetic.
    for w in s.neighbors[v]
        bv=s.inside[v];bw=s.inside[w]
        bv==bw && continue
        s.label[bw]==0 && s.mate[s.base[bw]]==0 || continue
        _weighted_slack(s,v,w)<=0 || continue
        _weighted_allow!(s,v,w)
        _seeded_label!(s,w,1,0);_weighted_augment!(s,v,w)
        return true
    end
    for w in s.neighbors[v]
        bv=s.inside[v];bw=s.inside[w];bv==bw && continue
        edge_id=s.edge_index[minmax(v,w)];slack=_weighted_slack(s,v,w)
        slack<=0 && (s.allowed[edge_id]=true)
        if s.allowed[edge_id]
            if s.label[bw]==0
                if s.mate[s.base[bw]]==0
                    _seeded_label!(s,w,1,0);_weighted_augment!(s,v,w)
                    return true
                end
                pending=length(s.queue);_seeded_label!(s,w,2,v)
                # Complete the newly discovered alternating subtree before
                # returning to the remaining incident edges of this vertex.
                _seeded_scan_batch!(s,pending) && return true
            elseif s.label[bw]==1
                base=_weighted_scan!(s,v,w)
                if base!=0
                    pending=length(s.queue);_seeded_add!(s,base,v,w)
                    # The new blossom's own leaves replace the inner entries
                    # the contraction just queued, anchored at its free base.
                    resize!(s.queue,pending)
                    blossom=s.inside[base]
                    leaves=_weighted_leaves(s,blossom)
                    at=findfirst(==(s.base[blossom]),leaves)
                    @assert at!==nothing
                    for i in at:length(leaves);push!(s.queue,leaves[i]);end
                    for i in 1:at-1;push!(s.queue,leaves[i]);end
                    _seeded_scan_batch!(s,pending) && return true
                else
                    _weighted_augment!(s,v,w);return true
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
    return false
end

function _seeded_reset_forest!(s)
    fill!(s.label,0);fill!(s.label_edge,(0,0));fill!(s.best,(0,0))
    for b in s.count+1:length(s.base)
        s.active[b] && (_weighted_blossom(s,b).best_valid=false)
    end
    fill!(s.allowed,false);empty!(s.queue)
    return nothing
end

function _seeded_dual_step!(s,delta)
    for vertex in 1:s.count
        if s.label[s.inside[vertex]]==1
            s.vertex_dual[vertex]=Base.checked_sub(s.vertex_dual[vertex],delta)
        elseif s.label[s.inside[vertex]]==2
            s.vertex_dual[vertex]=Base.checked_add(s.vertex_dual[vertex],delta)
        end
    end
    for b in s.count+1:length(s.base)
        s.active[b] && s.parent[b]==0 || continue
        if s.label[b]==1
            s.blossom_dual[b]=Base.checked_add(s.blossom_dual[b],delta)
        elseif s.label[b]==2
            s.blossom_dual[b]=Base.checked_sub(s.blossom_dual[b],delta)
        end
    end
    return nothing
end

# Normalize an initially half-integral free group through its actual tight
# alternating tree. Unreached matched pairs remain untouched. Tight
# connections can augment before any positive parity dual step.
function _seeded_normalize_cycles!(s,initial::_SeededInitial)
    count=s.count
    for vertices in initial.cycles
        at=0
        for i in eachindex(vertices)
            initial.mate[vertices[i]]==0 && (at=i;break)
        end
        @assert at!=0
        root=vertices[at]
        (s.mate[root]==0 && isodd(initial.dual[root])) || continue
        _seeded_reset_forest!(s)
        _seeded_label!(s,root,1,0)
        augmented=false
        while true
            while !isempty(s.queue) && !augmented
                vertex=popfirst!(s.queue)
                augmented=_seeded_depth_scan!(s,vertex)
            end
            augmented && break
            zero_inner=0
            for b in count+1:length(s.base)
                if s.active[b] && s.parent[b]==0 && s.label[b]==2 &&
                   s.blossom_dual[b]==0
                    zero_inner=b;break
                end
            end
            if zero_inner!=0
                _weighted_expand!(s,zero_inner,false)
                continue
            end
            _seeded_dual_step!(s,Int64(1))
            break
        end
    end
    return nothing
end

function _seeded_matching_run(count::Int,input,initial::_SeededInitial;
                              reverse_neighbors::Bool,allow_dual::Bool)
    neighbors,costs,edges,edge_index,directed=_weighted_input(count,input;directed=true)
    isempty(costs) && return zeros(Int,count)
    maximum_cost=maximum(costs)
    offset=Base.checked_add(maximum_cost,Int64(1))
    for index in eachindex(costs);costs[index]=Base.checked_sub(offset,costs[index]);end
    maximum_weight=maximum(costs)
    slots=2count
    storage=Vector{Union{Nothing,_WeightedBlossom}}(nothing,slots)
    s=_SeededWeightedMatching(count,neighbors,costs,edges,directed,edge_index,
        zeros(Int,count),collect(1:count),
        zeros(Int,slots),vcat(collect(1:count),zeros(Int,count)),zeros(Int,slots),
        fill((0,0),slots),fill((0,0),slots),fill(maximum_weight,count),zeros(Int64,slots),
        storage,vcat(trues(count),falses(count)),falses(slots),falses(length(edges)),
        Int[],collect(slots:-1:count+1),Int[],Int[],Int[],Int[],Int[],fill((0,0),slots),Int[],
        offset,Vector{Int}[],0,Int[],Int[],Vector{Union{Nothing,Vector{Int}}}(nothing,slots),
        zeros(Int,slots),
        zeros(Int8,slots),Int[],Int[],Int[],false,false,0,allow_dual)
    sizehint!(s.queue,count)
    s.mate[:]=initial.mate
    for vertex in 1:count
        s.vertex_dual[vertex]=Base.checked_sub(offset,initial.dual[vertex])
    end
    if allow_dual
        # Contract each retained tight odd ring at its actual free base; the
        # internal alternating pairs come from the verified seed mates.
        for vertices in initial.cycles
            at=0
            for i in eachindex(vertices)
                initial.mate[vertices[i]]==0 && (at=i;break)
            end
            @assert at!=0
            children=[vertices[mod1(at+i,length(vertices))]
                      for i in 0:length(vertices)-1]
            links=[(children[i],children[mod1(i+1,length(children))])
                   for i in eachindex(children)]
            isempty(s.free) && throw(ErrorException(
                "seeded weighted matching: blossom workspace exhausted"))
            b=pop!(s.free)
            s.base[b]=vertices[at];s.active[b]=true;s.initial_cycle[b]=true
            s.blossoms[b]=_WeightedBlossom(children,links,Tuple{Int,Int}[],false)
            for vertex in children;s.parent[vertex]=b;s.inside[vertex]=b;end
        end
    end
    reverse_neighbors && foreach(reverse!,s.neighbors)
    @assert all(edge->_weighted_slack(s,edge...)>=0,s.edges)
    @assert all(v->s.mate[v]==0 || (s.mate[s.mate[v]]==v &&
        _weighted_slack(s,v,s.mate[v])==0),1:count)
    allow_dual && _seeded_normalize_cycles!(s,initial)
    while true
        stage_root=0;restart_forest=false
        s.positive_tree_dual=false;s.root_rotation_pending=false
        _seeded_reset_forest!(s)
        roots=s.roots;empty!(roots)
        for vertex in count:-1:1
            s.mate[vertex]==0 && push!(roots,vertex)
        end
        for v in roots
            if s.mate[v]==0 && s.label[s.inside[v]]==0
                s.next_free_root!=0 && s.next_free_root!=v && continue
                stage_root=v;s.next_free_root=0
                _seeded_label!(s,v,1,0);break
            end
        end
        augmented=false
        while true
            while !isempty(s.queue) && !augmented
                v=popfirst!(s.queue)
                augmented=_seeded_depth_scan!(s,v)
            end
            augmented && break
            if s.root_rotation_pending && allow_dual && stage_root!=0
                # Cycle through the current root order. Contraction changes
                # group membership, so require another current free base.
                at=findfirst(==(stage_root),roots)
                for step in 1:length(roots)-1
                    next=roots[mod1(at+step,length(roots))]
                    s.mate[next]==0 && s.inside[next]!=s.inside[stage_root] || continue
                    s.next_free_root=next;restart_forest=true
                    break
                end
                s.root_rotation_pending=false
                restart_forest && break
            end
            allow_dual || break
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
            _seeded_dual_step!(s,delta)
            delta>0 && kind!=1 && (s.positive_tree_dual=true)
            kind==1 && break
            if kind==2 || kind==3
                _weighted_allow!(s,edge...)
                _seeded_rescan_queue!(s)
            else
                _weighted_expand!(s,blossom,false)
            end
        end
        restart_forest && continue
        for v in 1:count
            s.mate[v]!=0 && @assert s.mate[s.mate[v]]==v
        end
        augmented || break
        for b in count+1:length(s.base)
            !s.initial_cycle[b] && s.active[b] && s.parent[b]==0 && s.label[b]==1 &&
                s.blossom_dual[b]==0 && _weighted_expand!(s,b,true)
        end
    end
    !allow_dual && any(iszero,s.mate) && return nothing
    # An incomplete result means no perfect matching exists; the calling
    # contract reruns the cold kernel rather than certify a partial state.
    any(iszero,s.mate) && return s.mate
    _weighted_verify(s)
    if !allow_dual
        # No dual changed and every selected edge is tight, so twice the
        # matched cost equals the certified fractional singleton bound.
        twice=Int128(0)
        for v in 1:count
            twice+=Int128(Base.checked_sub(offset,
                s.weight[s.edge_index[minmax(v,s.mate[v])]]))
        end
        @assert twice==sum(Int128,initial.dual)
    end
    return s.mate
end

function minimum_weight_matching(count::Int,input,initial::_SeededInitial)
    owned=input isa AbstractVector ? input : collect(input)
    if !isempty(initial.cycles)
        # A pure tight-edge descent in reversed adjacency order can already
        # close every free base; when it does, the retained certificate makes
        # the result a certified integral optimum with no dual change.
        trial=_seeded_matching_run(count,owned,initial;
            reverse_neighbors=true,allow_dual=false)
        trial!==nothing && return trial
    end
    result=_seeded_matching_run(count,owned,initial;
        reverse_neighbors=false,allow_dual=true)
    # The seeded traversal promises the same maximum-cardinality contract as
    # the cold search; if it ever stops short, recover through the cold run.
    any(iszero,result) && return minimum_weight_matching(count,owned)
    return result
end
