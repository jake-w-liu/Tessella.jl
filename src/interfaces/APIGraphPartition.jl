# Native multilevel graph partitioning and renumbering.
#
# This file provides the quality partition tier that upstream Gmsh implements
# through METIS, plus the native renumbering tiers that replace METIS-based
# and space-filling orderings upstream exposes through
# `gmsh.model.mesh.computeRenumbering`. Everything here is plain Julia with no
# external library dependency.
#
# Element-dual graph semantics follow the public gmsh SDK `Graph` builder:
# every element is a vertex; two elements are adjacent when they share at
# least `numCommonNodesInDualGraph` primary nodes — the per-family threshold
# table below mirrors the upstream values (a facet's worth of nodes for
# same-dimension pairs, less across dimensions so lower-dimensional cells
# attach to the higher-dimensional cells they bound). Higher-order elements
# only contribute their primary (corner) vertices.
#
# The multilevel pipeline mirrors the documented METIS behaviour: sorted
# heavy-edge coarsening, greedy region-growing initial partition, and
# Fiduccia–Mattheyses / greedy boundary refinement while uncoarsening.
# All tie-breaking is index-based so results are deterministic.

# ---------------------------------------------------------------------------
# Element-dual graph construction
# ---------------------------------------------------------------------------

const _GP_FAMILIES=(:lin,:tri,:qua,:tet,:hex,:pri,:pyr,:trih)

# Upstream `MElement::numCommonNodesInDualGraph` table: minimum shared node
# count for two elements of families `fi`/`fj` to be dual-graph adjacent.
# Symmetric for every catalogued pair; the generic fallback matches
# upstream's polygon/polyhedron behaviour (any real geometric contact).
function _gp_dual_threshold(fi::Symbol,fj::Symbol)
    if fi===:pnt; return 1
    elseif fi===:lin; return fj in (:pnt,:lin) ? 1 : 2
    elseif fi===:tri; return fj===:pnt ? 1 : fj===:lin ? 2 : fj in (:tri,:qua) ? 2 : 3
    elseif fi===:qua; return fj===:pnt ? 1 : fj===:lin ? 2 : fj in (:tri,:qua) ? 2 : 4
    elseif fi===:tet; return fj===:pnt ? 1 : fj===:lin ? 2 : fj===:qua ? 4 : 3
    elseif fi===:hex; return fj===:pnt ? 1 : fj===:lin ? 2 : fj===:tri ? 3 : 4
    elseif fi in (:pri,:pyr); return fj===:pnt ? 1 : fj===:lin ? 2 : fj in (:qua,:hex) ? 4 : 3
    elseif fi===:trih; return fj===:pnt ? 1 : fj===:lin ? 2 : fj in (:qua,:hex) ? 4 : fj===:trih ? 2 : 3
    end
    return fj===:pnt ? 1 : fj===:lin ? 2 : 2
end

@inline function _gp_family_dim(fi::Symbol)
    fi===:pnt && return 0
    fi===:lin && return 1
    fi in (:tri,:qua) && return 2
    return 3
end

# Compact CSR dual graph plus per-vertex element weights.
struct _GPGraph
    xadj::Vector{Int32}
    adjncy::Vector{Int32}
    adjwgt::Vector{Int32}
    vwgt::Vector{Int32}
end
@inline _gp_nvertices(g::_GPGraph)=length(g.vwgt)

# Upstream `Graph::fillDefaultWeights`: the all-1 short-circuit and the
# -1 auto test inspect only the seven listed families — the trihedron weight
# is ignored upstream, and in the configured branch unlisted families
# (points, trihedra, polygons) weigh 1.
function _gp_element_weights(families::Vector{Symbol},weights::NTuple{8,Int})
    if all(==(1),weights[1:7])
        return ones(Int32,length(families))
    elseif any(<(0),weights[1:7])
        topdim=maximum(_gp_family_dim,families;init=0)
        return Int32[_gp_family_dim(f)==topdim ? 1 : 0 for f in families]
    else
        wmap=(lin=weights[7],tri=weights[1],qua=weights[2],tet=weights[3],
              hex=weights[4],pri=weights[5],pyr=weights[6])
        return Int32[get(wmap,f,1) for f in families]
    end
end

# Partition solver options collected from the session option dictionary or
# from an explicit NamedTuple supplied by the .geo executor (where the user
# may have assigned `Mesh.Metis*` options inside the script itself).
struct _GPOptions
    algorithm::Int        # Mesh.MetisAlgorithm: 1=recursive bisection
                          # (upstream default), 2=k-way; other values fall
                          # back to recursive. -1 = Tessella extension
                          # selecting the legacy BFS tier.
    edge_matching::Int    # Mesh.MetisEdgeMatching: 1=random 2=sorted
                          # heavy-edge matching
    ufactor::Int          # Mesh.MetisMaxLoadImbalance, 1e-3 fraction;
                          # <0 uses the algorithm default
    objective::Int        # Mesh.MetisObjective: 1=cut 2=volume
    refinement::Int       # Mesh.MetisRefinementAlgorithm: 1=FM
                          # 2=greedy boundary
    minconn::Int          # Mesh.MetisMinConn
    weights::NTuple{8,Int}
end

# Upstream balance defaults (documented under `Mesh.MetisMaxLoadImbalance`):
# recursive bisection allows 1 per-mille load imbalance, k-way 30.
_gp_default_ufactor(algorithm::Int)=algorithm==2 ? 30 : 1

function _gp_weights_ntuple(getval)
    return (getval("Mesh.PartitionTriWeight",-1),getval("Mesh.PartitionQuadWeight",-1),
            getval("Mesh.PartitionTetWeight",-1),getval("Mesh.PartitionHexWeight",-1),
            getval("Mesh.PartitionPrismWeight",-1),getval("Mesh.PartitionPyramidWeight",-1),
            getval("Mesh.PartitionLineWeight",-1),getval("Mesh.PartitionTrihedronWeight",0))
end

function _gp_options(opts)
    if opts===nothing
        getval=(k,d)->Int(round(get(OPTIONS,k,Float64(d))))
    else
        getval=(k,d)->Int(round(get(opts,Symbol(replace(k,"Mesh."=>"")),d)))
    end
    return _GPOptions(getval("Mesh.MetisAlgorithm",1),
                      getval("Mesh.MetisEdgeMatching",2),
                      getval("Mesh.MetisMaxLoadImbalance",-1),
                      getval("Mesh.MetisObjective",1),
                      getval("Mesh.MetisRefinementAlgorithm",2),
                      getval("Mesh.MetisMinConn",-1),
                      _gp_weights_ntuple(getval))
end

# Iterate (family-symbol, msh-type, cells) triples in the same catalog order
# `_partition_elements_simple` numbers elements, so vertex ids stay aligned
# with the dense assignment indices the partition producers consume.
_partition_dual_blocks(mesh)=
    Tuple((msh_family(msh),Int(msh),cells)
          for (msh,cells,_) in _cache_native_blocks(mesh))

# Build the undirected element-dual graph. Returns `(graph, families)` with
# vertices numbered in mesh element order (block order, then column).
function _gp_build_dual(mesh)
    blocks=_partition_dual_blocks(mesh)
    total=sum((size(cells,2) for (_,_,cells) in blocks);init=0)
    families=Vector{Symbol}(undef,total)
    nn=nnodes(mesh)
    incidence=[Int32[] for _ in 1:nn]
    column=0
    # Primary (corner) vertices only, matching upstream `numPrimaryVertices`.
    for (family,_,cells) in blocks
        nprim=Elements._msh_primary_nodes(family)
        for c in axes(cells,2)
            column+=1
            families[column]=family
            for r in 1:nprim
                push!(incidence[Int(cells[r,c])],Int32(column))
            end
        end
    end
    # For each element, walk its nodes' incidence lists and count shared
    # nodes per candidate neighbour; emit an undirected edge when either
    # endpoint's threshold rule accepts (the table is symmetric, so this is
    # exactly `count >= threshold`).
    counts=zeros(Int32,total)
    touched=Int32[]
    xadj=Vector{Int32}(undef,total+1);xadj[1]=1
    adjncy=Int32[]
    column=0
    for (family,_,cells) in blocks
        nprim=Elements._msh_primary_nodes(family)
        for c in axes(cells,2)
            column+=1
            i=column
            fi=family
            empty!(touched)
            for r in 1:nprim
                node=Int(cells[r,c])
                for j in incidence[node]
                    j==i && continue
                    if counts[j]==0; push!(touched,j); end
                    counts[j]+=1
                end
            end
            start=length(adjncy)+1
            for j in touched
                cj=counts[j];counts[j]=0
                thr=min(_gp_dual_threshold(fi,families[j]),
                        _gp_dual_threshold(families[j],fi))
                cj>=thr || continue
                push!(adjncy,j)
            end
            start<=length(adjncy) && sort!(@view adjncy[start:end])
            xadj[i+1]=length(adjncy)+1
        end
    end
    return _GPGraph(xadj,adjncy,ones(Int32,length(adjncy)),
                    ones(Int32,total)),families
end

# ---------------------------------------------------------------------------
# Deterministic binary max-heap (key, tiebreak) — used by all greedy loops.
# Stored as (priority, -index); comparisons are lexicographic so equal
# priorities resolve to the smallest vertex index.
# ---------------------------------------------------------------------------

function _gp_heap_push!(heap::Vector{Tuple{Int32,Int32}},item)
    push!(heap,item)
    i=length(heap)
    while i>1
        p=i>>>1
        heap[p]>=heap[i] && break
        heap[p],heap[i]=heap[i],heap[p];i=p
    end
    return heap
end

function _gp_heap_pop!(heap::Vector{Tuple{Int32,Int32}})
    top=heap[1]
    last=pop!(heap)
    n=length(heap)
    if n>0
        heap[1]=last;i=1
        while true
            l=2i;r=l+1;c=i
            l<=n && heap[l]>heap[c] && (c=l)
            r<=n && heap[r]>heap[c] && (c=r)
            c==i && break
            heap[c],heap[i]=heap[i],heap[c];i=c
        end
    end
    return top
end

# ---------------------------------------------------------------------------
# Coarsening: sorted heavy-edge matching (SHEM)
# ---------------------------------------------------------------------------

# Vertices are visited in ascending (degree, index) order for determinism
# (`Mesh.MetisEdgeMatching=3`, the upstream default). Each vertex pairs with
# its heaviest unmatched neighbour (lowest index on ties); unmatched vertices
# stay singleton.
function _gp_match(g::_GPGraph,sorted::Bool)
    n=_gp_nvertices(g)
    match=zeros(Int32,n)
    if sorted
        degree=Int32[g.xadj[i+1]-g.xadj[i] for i in 1:n]
        order=sortperm(degree;alg=Base.Sort.DEFAULT_STABLE)
    else
        order=collect(1:n)
    end
    ncoarse=0
    for v in order
        match[v]!=0 && continue
        best=Int32(0);bestw=Int32(-1)
        for e in g.xadj[v]:g.xadj[v+1]-1
            u=g.adjncy[e]
            match[u]!=0 && continue
            w=g.adjwgt[e]
            if w>bestw || (w==bestw && best!=0 && u<best)
                best=u;bestw=w
            end
        end
        ncoarse+=1
        best==0 ? match[v]=Int32(ncoarse) :
                  (match[v]=match[best]=Int32(ncoarse))
    end
    return ncoarse,match
end

# Build the coarse CSR graph induced by `match` (fine vertex → coarse id).
function _gp_coarsen(g::_GPGraph,ncoarse::Int,match::Vector{Int32})
    members=[Int32[] for _ in 1:ncoarse]
    for v in 1:_gp_nvertices(g);push!(members[match[v]],Int32(v));end
    cvwgt=zeros(Int32,ncoarse)
    wacc=zeros(Int32,ncoarse)
    seen=Int32[]
    xadj=Vector{Int32}(undef,ncoarse+1);xadj[1]=1
    adjncy=Int32[];adjwgt=Int32[]
    for c in 1:ncoarse
        empty!(seen)
        for v in members[c]
            cvwgt[c]+=g.vwgt[v]
            for e in g.xadj[v]:g.xadj[v+1]-1
                cu=match[g.adjncy[e]]
                cu==c && continue
                wacc[cu]==0 && push!(seen,cu)
                wacc[cu]+=g.adjwgt[e]
            end
        end
        sort!(seen)
        for cu in seen
            push!(adjncy,cu);push!(adjwgt,wacc[cu]);wacc[cu]=0
        end
        xadj[c+1]=length(adjncy)+1
    end
    return _GPGraph(xadj,adjncy,adjwgt,cvwgt)
end

# Full coarsening chain. Returns (levels, cmaps); levels[1] is the input
# graph, cmaps[i] maps level-i vertices to level-(i+1) ids.
function _gp_coarsen_until(g::_GPGraph,limit::Int,edge_matching::Int)
    levels=_GPGraph[g]
    cmaps=Vector{Int32}[]
    # MetisEdgeMatching: 1=random (index order here — deterministic), the
    # default 2/METIS default is sorted heavy-edge matching.
    sorted=edge_matching!=1
    while _gp_nvertices(levels[end])>limit
        ncoarse,match=_gp_match(levels[end],sorted)
        # Stop when matching can no longer shrink the graph meaningfully.
        ncoarse>=ceil(Int,_gp_nvertices(levels[end])*0.9) && break
        coarse=_gp_coarsen(levels[end],ncoarse,match)
        push!(cmaps,match);push!(levels,coarse)
    end
    return levels,cmaps
end

# ---------------------------------------------------------------------------
# Initial bisection: pseudo-peripheral seed + greedy gain growth
# ---------------------------------------------------------------------------

# Farthest vertex from `seeds` under BFS distance (deterministic ties).
function _gp_bfs_distances(g::_GPGraph,sources)
    n=_gp_nvertices(g)
    dist=fill(Int32(-1),n)
    queue=Int32[]
    for s in sources;dist[s]=0;push!(queue,Int32(s));end
    head=1
    while head<=length(queue)
        v=queue[head];head+=1
        for e in g.xadj[v]:g.xadj[v+1]-1
            u=g.adjncy[e]
            dist[u]<0 || continue
            dist[u]=dist[v]+1;push!(queue,u)
        end
    end
    return dist
end

function _gp_pseudo_peripheral(g::_GPGraph)
    n=_gp_nvertices(g)
    degree=Int32[g.xadj[i+1]-g.xadj[i] for i in 1:n]
    root=argmin(degree)
    dist=_gp_bfs_distances(g,(root,))
    best=maximum(dist)
    best<0 && return root
    seed=findfirst(==(best),dist)
    dist=_gp_bfs_distances(g,(seed,))
    best=maximum(dist)
    best<=0 && return seed
    return findfirst(==(best),dist)
end

# Greedy region growing: assign weight `target` to side A, expanding from
# `seed` by maximum gain (edge weight into A; ties → smallest index).
function _gp_grow_bisect(g::_GPGraph,seed::Int,target::Int)
    n=_gp_nvertices(g)
    side=falses(n)
    gain=zeros(Int32,n)
    heap=Tuple{Int32,Int32}[]
    side[seed]=true;wA=Int(g.vwgt[seed])
    function push_neighbors(v)
        for e in g.xadj[v]:g.xadj[v+1]-1
            u=g.adjncy[e]
            side[u] && continue
            gain[u]+=g.adjwgt[e]
            _gp_heap_push!(heap,(gain[u],Int32(-u)))
        end
    end
    push_neighbors(seed)
    while wA<target && !isempty(heap)
        gv,nu=_gp_heap_pop!(heap);u=-nu
        (side[u] || gain[u]!=gv) && continue   # stale
        side[u]=true;wA+=g.vwgt[u]
        push_neighbors(u)
    end
    # Any vertices left unreachable (disconnected) join B; callers handle
    # connectivity afterwards.
    return side,wA
end

# ---------------------------------------------------------------------------
# FM boundary refinement for a bisection
# ---------------------------------------------------------------------------

# Compute the current edge cut of a bisection.
function _gp_bisect_cut(g::_GPGraph,side)
    cut=Int64(0)
    for v in 1:_gp_nvertices(g)
        for e in g.xadj[v]:g.xadj[v+1]-1
            u=g.adjncy[e]
            u>v && side[u]!=side[v] && (cut+=g.adjwgt[e])
        end
    end
    return cut
end

# One FM pass over a bisection. Any boundary vertex may move; the pass rolls
# back to the best prefix measured as (balance excess, cut) — so restoring an
# out-of-window side always outranks a marginal cut improvement, matching
# METIS's behaviour of first satisfying the balance constraint.
function _gp_fm_bisect_pass!(g::_GPGraph,side::BitVector,wA::Int,
                             loA::Int,hiA::Int)
    n=_gp_nvertices(g)
    gain=zeros(Int32,n)
    heap=Tuple{Int32,Int32}[]
    for v in 1:n
        inw=Int32(0);exw=Int32(0)
        for e in g.xadj[v]:g.xadj[v+1]-1
            w=g.adjwgt[e]
            side[g.adjncy[e]]==side[v] ? inw+=w : exw+=w
        end
        gain[v]=exw-inw
        exw>0 && _gp_heap_push!(heap,(gain[v],Int32(-v)))
    end
    locked=falses(n)
    cut=_gp_bisect_cut(g,side)
    excess(w)=max(loA-w,w-hiA,0)
    bestex=excess(wA);bestcut=cut;bestlen=0
    moved_v=Int32[]
    while !isempty(heap)
        gv,nu=_gp_heap_pop!(heap);v=-nu
        (locked[v] || gv!=gain[v]) && continue   # stale entry
        delta=side[v] ? -Int(g.vwgt[v]) : Int(g.vwgt[v])
        locked[v]=true
        side[v]=!side[v];wA+=delta;cut-=gv
        push!(moved_v,Int32(v))
        # Re-key unlocked neighbours: their gains shift by ±2w as v flips.
        for e in g.xadj[v]:g.xadj[v+1]-1
            u=g.adjncy[e]
            locked[u] && continue
            w=g.adjwgt[e]
            side[u]==side[v] ? (gain[u]-=2w) : (gain[u]+=2w)
            _gp_heap_push!(heap,(gain[u],Int32(-u)))
        end
        ex=excess(wA)
        if ex<bestex || (ex==bestex && cut<bestcut)
            bestex=ex;bestcut=cut;bestlen=length(moved_v)
        end
    end
    for k in length(moved_v):-1:bestlen+1
        v=moved_v[k];side[v]=!side[v]
    end
    wA=0
    for v in 1:n; side[v] && (wA+=g.vwgt[v]); end
    return bestcut,Int(wA)
end

# Refine a bisection until no pass improves (bounded passes).
function _gp_refine_bisect!(g::_GPGraph,side::BitVector,wA::Int,
                          loA::Int,hiA::Int;max_passes::Int=10)
    cut=_gp_bisect_cut(g,side)
    for _ in 1:max_passes
        newcut,wA=_gp_fm_bisect_pass!(g,side,wA,loA,hiA)
        newcut<cut || break
        cut=newcut
    end
    return wA
end

# ---------------------------------------------------------------------------
# Multilevel bisection
# ---------------------------------------------------------------------------

# A few deterministic starting seeds (METIS's `niter` role): pseudo-
# peripheral, min/max-degree, first and median vertices.
function _gp_seed_set(g::_GPGraph)
    n=_gp_nvertices(g)
    degree=Int32[g.xadj[i+1]-g.xadj[i] for i in 1:n]
    seeds=Int[_gp_pseudo_peripheral(g),argmin(degree),argmax(degree),
              1,(n+1)÷2]
    return unique!(seeds)
end

# Bisect `g` so side A carries about `targetA` weight within `slack`.
# Returns the side bitvector (true = A).
function _gp_bisect_ml(g::_GPGraph,targetA::Int,slack::Int,
                       edge_matching::Int,refinement::Int)
    n=_gp_nvertices(g)
    n<=1 && return trues(n)
    levels,cmaps=_gp_coarsen_until(g,max(32,2),edge_matching)
    coarsest=levels[end]
    # Multi-seed initial partition on the (small) coarsest graph: grow+refine
    # from each seed, keep the best cut — deterministic `niter` role.
    side=falses(_gp_nvertices(coarsest));wA=0;bestcut=typemax(Int64)
    for seed in _gp_seed_set(coarsest)
        cand,cw=_gp_grow_bisect(coarsest,seed,targetA)
        refinement>0 &&
            (cw=_gp_refine_bisect!(coarsest,cand,cw,
                                   targetA-slack,targetA+slack))
        c=_gp_bisect_cut(coarsest,cand)
        if c<bestcut
            bestcut=c;side=cand;wA=cw
        end
    end
    # Uncoarsen: project the coarse bisection down and refine each level.
    for level in length(levels)-1:-1:1
        cmap=cmaps[level]
        fine=levels[level]
        fineside=falses(_gp_nvertices(fine))
        fwA=0
        for v in 1:_gp_nvertices(fine)
            fineside[v]=side[cmap[v]]
            fineside[v] && (fwA+=fine.vwgt[v])
        end
        side=fineside;wA=Int(fwA)
        refinement>0 &&
            (wA=_gp_refine_bisect!(fine,side,wA,targetA-slack,targetA+slack))
    end
    return side
end

# ---------------------------------------------------------------------------
# Multilevel k-way partitioning
# ---------------------------------------------------------------------------

# Initial k-way partition of the coarsest graph: farthest-first seeds, then
# per-part greedy growth to weight targets.
function _gp_initial_kway(g::_GPGraph,nparts::Int,targets::Vector{Int})
    n=_gp_nvertices(g)
    part=zeros(Int32,n)
    degree=Int32[g.xadj[i+1]-g.xadj[i] for i in 1:n]
    seeds=Int32[argmin(degree)]
    for p in 2:nparts
        dist=_gp_bfs_distances(g,Tuple(Int(s) for s in seeds))
        dist2=copy(dist)
        # farthest reachable vertex from all existing seeds; fall back to the
        # lowest unassigned vertex when the graph is disconnected.
        best=-1;seed=0
        for v in 1:n
            part[v]!=0 && continue
            d=dist2[v]<0 ? typemax(Int32)-1 : dist2[v]
            d>best && (best=d;seed=v)
        end
        seed==0 && break
        push!(seeds,Int32(seed))
    end
    # Grow each part round-robin by gain until targets are met.
    gain=zeros(Int32,n)
    heaps=[Tuple{Int32,Int32}[] for _ in 1:nparts]
    partw=zeros(Int,nparts)
    for (p,s) in enumerate(seeds)
        part[s]==0 || continue
        part[s]=p;partw[p]=g.vwgt[s]
        for e in g.xadj[s]:g.xadj[s+1]-1
            u=g.adjncy[e]
            part[u]==0 && _gp_heap_push!(heaps[p],(g.adjwgt[e],Int32(-u)))
        end
    end
    changed=true
    while changed && !all(p->partw[p]>=targets[p],1:nparts)
        changed=false
        for p in 1:nparts
            partw[p]>=targets[p] && continue
            while !isempty(heaps[p])
                gv,nu=_gp_heap_pop!(heaps[p]);u=-nu
                part[u]!=0 && continue
                # recompute gain into part p for freshness
                gw=Int32(0)
                for e in g.xadj[u]:g.xadj[u+1]-1
                    part[g.adjncy[e]]==p && (gw+=g.adjwgt[e])
                end
                gw!=gv && (_gp_heap_push!(heaps[p],(gw,Int32(-u)));continue)
                part[u]=p;partw[p]+=g.vwgt[u];changed=true
                for e in g.xadj[u]:g.xadj[u+1]-1
                    w=g.adjncy[e]
                    part[w]==0 &&
                        _gp_heap_push!(heaps[p],(g.adjwgt[e],Int32(-w)))
                end
                break
            end
        end
    end
    # Leftover/unreachable vertices join the adjacent part sharing the most
    # weight, else the lightest part.
    for v in 1:n
        part[v]!=0 && continue
        wsum=zeros(Int32,nparts)
        for e in g.xadj[v]:g.xadj[v+1]-1
            p=part[g.adjncy[e]]
            p!=0 && (wsum[p]+=g.adjwgt[e])
        end
        if any(>(0),wsum)
            part[v]=argmax(wsum)
        else
            part[v]=argmin(partw)
        end
        partw[part[v]]+=g.vwgt[v]
    end
    return part
end

# Greedy k-way boundary refinement pass (upstream REFINEMENT_GREEDY /
# FM style): vertices move to the part they share the most weight with,
# subject to part bounds.
function _gp_refine_kway!(g::_GPGraph,part::Vector{Int32},nparts::Int,
                        lo::Vector{Int},hi::Vector{Int};max_passes::Int=10)
    n=_gp_nvertices(g)
    partw=zeros(Int,nparts)
    for v in 1:n
        partw[part[v]]+=Int(g.vwgt[v])
    end
    excess(w,p)=max(lo[p]-w,w-hi[p],0)
    for _ in 1:max_passes
        moved=false
        wsum=zeros(Int32,nparts)
        for v in 1:n
            own=part[v]
            fill!(wsum,Int32(0))
            for e in g.xadj[v]:g.xadj[v+1]-1
                wsum[part[g.adjncy[e]]]+=g.adjwgt[e]
            end
            wv=Int(g.vwgt[v])
            own_after=excess(partw[own]-wv,own)
            own_before=excess(partw[own],own)
            best=0;bestgain=typemin(Int)
            for p in 1:nparts
                p==own && continue
                gain_p=Int(wsum[p])-Int(wsum[own])
                before=own_before+excess(partw[p],p)
                after=own_after+excess(partw[p]+wv,p)
                # Cut-improving moves must not worsen balance; balance-
                # repairing moves run regardless of gain.
                ((gain_p>0 && after<=before) || after<before) || continue
                (gain_p>bestgain || (gain_p==bestgain && p<best)) ||
                    continue
                bestgain=gain_p;best=p
            end
            best==0 && continue
            part[v]=best;partw[own]-=wv;partw[best]+=wv
            moved=true
        end
        moved || break
    end
    return part
end

# Multilevel k-way partition. Returns 1-based part labels.
function _gp_kway_ml(g::_GPGraph,nparts::Int,ufactor::Int,
                     edge_matching::Int,refinement::Int)
    n=_gp_nvertices(g)
    total=Int(sum(g.vwgt))
    slack=max(1,total*ufactor÷1000)
    targets=[round(Int,total*p/nparts)-round(Int,total*(p-1)/nparts)
             for p in 1:nparts]
    lo=[max(0,t-slack) for t in targets]
    hi=[t+slack for t in targets]
    levels,cmaps=_gp_coarsen_until(g,max(32,nparts*4),edge_matching)
    part=_gp_initial_kway(levels[end],nparts,targets)
    refinement>0 && _gp_refine_kway!(levels[end],part,nparts,lo,hi)
    for level in length(levels)-1:-1:1
        cmap=cmaps[level]
        fine=levels[level]
        fpart=Vector{Int32}(undef,_gp_nvertices(fine))
        for v in 1:_gp_nvertices(fine);fpart[v]=part[cmap[v]];end
        part=fpart
        refinement>0 && _gp_refine_kway!(fine,part,nparts,lo,hi)
    end
    return part
end

# ---------------------------------------------------------------------------
# Recursive bisection driver
# ---------------------------------------------------------------------------

# Partition the induced subgraph on `verts` into `nparts` via recursive
# bisection. Returns a length-nparts Vector{Vector{Int32}} of fine vertices.
function _gp_recursive_split(g::_GPGraph,verts::Vector{Int32},nparts::Int,
                             ufactor::Int,edge_matching::Int,refinement::Int)
    nparts<=1 && return [verts]
    total=sum(Int(g.vwgt[v]) for v in verts)
    nleft=nparts÷2;nright=nparts-nleft
    targetA=round(Int,total*nleft/nparts)
    slack=max(1,total*ufactor÷1000)
    # Build the induced subgraph on `verts` (compact CSR, same weights).
    index=zeros(Int32,_gp_nvertices(g))
    for (i,v) in enumerate(verts);index[v]=i;end
    xadj=Vector{Int32}(undef,length(verts)+1);xadj[1]=1
    adjncy=Int32[];adjwgt=Int32[];vwgt=Vector{Int32}(undef,length(verts))
    for (i,v) in enumerate(verts)
        vwgt[i]=g.vwgt[v]
        for e in g.xadj[v]:g.xadj[v+1]-1
            j=index[g.adjncy[e]]
            j==0 && continue
            push!(adjncy,j);push!(adjwgt,g.adjwgt[e])
        end
        xadj[i+1]=length(adjncy)+1
    end
    sub=_GPGraph(xadj,adjncy,adjwgt,vwgt)
    side=_gp_bisect_ml(sub,targetA,slack,edge_matching,refinement)
    left=Int32[verts[i] for i in 1:length(verts) if side[i]]
    right=Int32[verts[i] for i in 1:length(verts) if !side[i]]
    out=_gp_recursive_split(g,left,nleft,ufactor,edge_matching,refinement)
    append!(out,_gp_recursive_split(g,right,nright,ufactor,
                                    edge_matching,refinement))
    return out
end

# ---------------------------------------------------------------------------
# Connectivity, components and relabelling
# ---------------------------------------------------------------------------

# Connected components restricted to an index set.
function _gp_components_in(g::_GPGraph,members::Vector{Int32})
    inset=falses(_gp_nvertices(g))
    for v in members;inset[v]=true;end
    comps=Vector{Int32}[]
    seen=falses(_gp_nvertices(g))
    for s in members
        seen[s] && continue
        push!(comps,Int32[]);c=length(comps)
        queue=Int32[s];seen[s]=true
        while !isempty(queue)
            v=pop!(queue);push!(comps[c],v)
            for e in g.xadj[v]:g.xadj[v+1]-1
                u=g.adjncy[e]
                (inset[u] && !seen[u]) || continue
                seen[u]=true;push!(queue,u)
            end
        end
        sort!(comps[c])
    end
    return comps
end

# Repair connectivity: for each part, keep its largest component and donate
# the others to the adjacent part sharing the most edge weight. Attaching a
# donated component to a neighbouring region preserves the receiver's
# connectivity, so the total component count strictly decreases.
function _gp_repair_connectivity!(g::_GPGraph,part::Vector{Int32},nparts::Int)
    for _ in 1:_gp_nvertices(g)+nparts
        dirty=false
        members=[Int32[] for _ in 1:nparts]
        for v in 1:_gp_nvertices(g);push!(members[part[v]],v);end
        for p in 1:nparts
            comps=_gp_components_in(g,members[p])
            length(comps)<=1 && continue
            dirty=true
            sort!(comps;by=c->(-length(c),c[1]))
            keep=comps[1]
            for c in 2:length(comps)
                wsum=zeros(Int32,nparts)
                for v in comps[c]
                    for e in g.xadj[v]:g.xadj[v+1]-1
                        wsum[part[g.adjncy[e]]]+=g.adjwgt[e]
                    end
                end
                wsum[p]=0
                dest=any(>(0),wsum) ? argmax(wsum) :
                    argmin([length(members[q]) for q in 1:nparts])
                for v in comps[c];part[v]=dest;end
            end
        end
        dirty || break
    end
    return part
end

# Re-label parts so the part containing vertex 1 is 1, and each successive
# new part first met scanning vertex order takes the next free label. This
# canonicalises label ordering across paths (bisection vs k-way vs simple).
function _gp_relabel!(part::Vector{Int32})
    map=zeros(Int32,maximum(part;init=0))
    next=0
    for v in 1:length(part)
        p=part[v]
        if map[p]==0;next+=1;map[p]=next;end
        part[v]=map[p]
    end
    return part
end

# ---------------------------------------------------------------------------
# Top-level partition entry
# ---------------------------------------------------------------------------

# Legacy contiguous quota split used for edge-less duals and degenerate
# cases — identical to the BFS tier's assignment for those inputs.
function _gp_contiguous_split(total::Int,nparts::Int)
    assignment=Vector{Int32}(undef,total)
    quota=cld(total,nparts)
    for e in 1:total
        assignment[e]=Int32(min(cld(e,quota),nparts))
    end
    return assignment
end

# Partition `g` into `nparts`. Vertex labels are 1-based, deterministic, and
# connected whenever the dual allows it.
function _gp_partition(g::_GPGraph,nparts::Int,opts::_GPOptions)
    n=_gp_nvertices(g)
    nparts<=1 && return ones(Int32,n)
    n==0 && return Int32[]
    # More parts than vertices — or an edge-less dual such as a soup of
    # isolated elements — has no cut to optimise: keep the legacy contiguous
    # quota convention (identical to the BFS tier's assignment).
    if nparts>=n || isempty(g.adjncy)
        return _gp_contiguous_split(n,nparts)
    end
    ufactor=opts.ufactor<0 ? _gp_default_ufactor(opts.algorithm) : opts.ufactor
    comps=_gp_components_in(g,collect(Int32,1:n))
    part=zeros(Int32,n)
    if length(comps)>1
        # Disconnected dual: allocate whole components to part ranges
        # proportional to weight so each part stays inside one component
        # whenever possible.
        total=sum(Int,v for v in g.vwgt)
        compw=[sum(Int(g.vwgt[v]) for v in c) for c in comps]
        shares=[max(1,round(Int,nparts*w/total)) for w in compw]
        # Normalise shares to nparts.
        while sum(shares)>nparts
            i=argmax(shares);shares[i]-=1
        end
        rest=nparts-sum(shares)
        while rest>0
            i=argmax(compw./shares);shares[i]+=1;rest-=1
        end
        offset=0
        for (c,verts) in enumerate(comps)
            k=shares[c]
            if k<=1
                for v in verts;part[v]=offset+1;end
            else
                pieces=_gp_recursive_split(g,verts,k,ufactor,
                                           opts.edge_matching,opts.refinement)
                for (q,piece) in enumerate(pieces)
                    for v in piece;part[v]=offset+q;end
                end
            end
            offset+=k
        end
        if offset<nparts
            for v in 1:n
                part[v]==0 && (part[v]=nparts)
            end
        end
    elseif opts.algorithm==2
        part=_gp_kway_ml(g,nparts,ufactor,opts.edge_matching,opts.refinement)
    else
        pieces=_gp_recursive_split(g,collect(Int32,1:n),nparts,ufactor,
                                   opts.edge_matching,opts.refinement)
        for (q,piece) in enumerate(pieces)
            for v in piece;part[v]=q;end
        end
    end
    # Final boundary polish on the fine graph regardless of path: greedy
    # k-way refinement with the configured balance bounds (the upstream
    # refinement stage runs at full resolution too).
    if opts.refinement>0 && nparts>1
        total=sum(Int,g.vwgt)
        slack=max(1,total*ufactor÷1000)
        targets=[round(Int,total*p/nparts)-round(Int,total*(p-1)/nparts)
                 for p in 1:nparts]
        lo=[max(0,t-slack) for t in targets]
        hi=[t+slack for t in targets]
        _gp_refine_kway!(g,part,nparts,lo,hi;max_passes=4)
    end
    _gp_repair_connectivity!(g,part,nparts)
    return _gp_relabel!(part)
end

# Public entry used by `_partition_elements`: build the dual, weight it per
# the configured `Mesh.Partition*Weight` options and partition.
function _gp_partition_mesh(mesh,nparts::Int,opts::_GPOptions)
    # A single partition needs no dual graph — identical to the BFS tier.
    nparts<=1 && return ones(Int32,
        sum((size(cells,2) for (_,_,cells) in _partition_dual_blocks(mesh));
            init=0))
    graph,families=_gp_build_dual(mesh)
    graph=_GPGraph(graph.xadj,graph.adjncy,graph.adjwgt,
                   _gp_element_weights(families,opts.weights))
    return _gp_partition(graph,nparts,opts)
end

# `_partition_elements` — the quality partition tier. `Mesh.MetisAlgorithm`
# selects the solver: 1 (upstream default) runs recursive multilevel
# bisection, 2 runs the multilevel k-way path; the Tessella extension value
# -1 selects the legacy linear-time BFS tier (`_partition_elements_simple`,
# which remains what `SimplePartition`-style geometric partitioning
# approximates for an upstream no-METIS build).
function _partition_elements(mesh::Union{Mesh,MixedMesh},num_partitions::Int;
                             options=nothing)
    opts=_gp_options(options)
    opts.algorithm<0 &&
        return _partition_elements_simple(mesh,num_partitions)
    return _gp_partition_mesh(mesh,num_partitions,opts)
end

# ---------------------------------------------------------------------------
# Renumbering: Hilbert space-filling order + 2-opt polish + nested dissection
# ---------------------------------------------------------------------------

# Skilling axes-to-transpose: converts `X` (n coordinates of `b` bits each)
# to the transposed Hilbert integer bits in-place. Standard public-domain
# algorithm used for Hilbert curve indexing.
function _gp_hilbert_transpose!(X::Vector{UInt64},b::Int)
    n=length(X)
    M=UInt64(1)<<(b-1)
    Q=M
    while Q>1
        P=Q-1
        for i in 1:n
            if (X[i]&Q)!=0
                X[1]⊻=P
            else
                t=(X[1]⊻X[i])&P;X[1]⊻=t;X[i]⊻=t
            end
        end
        Q>>=1
    end
    for i in 2:n;X[i]⊻=X[i-1];end
    t=UInt64(0);Q=M
    while Q>1
        (X[n]&Q)!=0 && (t⊻=Q-1)
        Q>>=1
    end
    for i in 1:n;X[i]⊻=t;end
    return X
end

# Hilbert index for a point quantized to `b` bits per axis (n dims ≤ 4).
function _gp_hilbert_index(coords::Vector{UInt64},b::Int)
    X=_gp_hilbert_transpose!(copy(coords),b)
    # Interleave bits, most significant first.
    n=length(X);index=zero(UInt64)
    for bit in b-1:-1:0
        for i in 1:n
            index=(index<<1)|((X[i]>>bit)&1)
        end
    end
    return index
end

# Sort node positions by Hilbert curve order. `coords` is 3×n. Returns the
# permutation (order[new_position] = old node index).
function _mesh_hilbert_node_order(coords::AbstractMatrix)
    n=size(coords,2)
    n<=1 && return collect(1:n)
    b=21  # 63-bit index for 3-D
    lo=vec(minimum(coords,dims=2))
    hi=vec(maximum(coords,dims=2))
    scale=[hi[d]>lo[d] ? (UInt64(1)<<b-1)/(hi[d]-lo[d]) : 0.0 for d in 1:3]
    keys=Vector{UInt64}(undef,n)
    X=Vector{UInt64}(undef,3)
    for i in 1:n
        for d in 1:3
            t=hi[d]>lo[d] ? (coords[d,i]-lo[d])*scale[d] : 0.0
            X[d]=UInt64(clamp(floor(Int64,t),0,Int64((UInt64(1)<<b)-1)))
        end
        keys[i]=_gp_hilbert_index(X,b)
    end
    order=sortperm(keys;alg=Base.Sort.DEFAULT_STABLE)
    # 2-opt polish on the ordered path through physical space.
    _gp_two_opt!(order,coords)
    return order
end

# Deterministic bounded 2-opt pass: scan adjacent windows, reverse any
# segment that shortens the path. Improves cache-locality of the ordering.
function _gp_two_opt!(order::Vector{Int},coords::AbstractMatrix;
                      passes::Int=4,window::Int=32)
    n=length(order)
    n<4 && return order
    dist(a,b)=((coords[1,a]-coords[1,b])^2+(coords[2,a]-coords[2,b])^2+
               (coords[3,a]-coords[3,b])^2)
    for _ in 1:passes
        improved=false
        for i in 1:n-1
            a=order[i];a2=order[i+1]
            da=dist(a,a2)
            for j in i+2:min(n-1,i+window)
                b=order[j];b2=order[j+1]
                if dist(a,b)+dist(a2,b2)<da+dist(b,b2)-1e-18
                    reverse!(@view order[i+1:j])
                    improved=true
                    a2=order[i+1];da=dist(a,a2)
                end
            end
        end
        improved || break
    end
    return order
end

# Nested-dissection ordering on the node adjacency (list-of-sets like the
# RCM path uses). Recursively bisects the vertex-induced graph, ordering
# interior(A), then B, then the boundary vertex separator covering the cut.
function _mesh_node_nd_order(adjacency)
    n=length(adjacency)
    # Convert to CSR for the partitioner machinery.
    xadj=Vector{Int32}(undef,n+1);xadj[1]=1
    adjncy=Int32[]
    for v in 1:n
        for u in sort!(collect(adjacency[v]))
            push!(adjncy,Int32(u))
        end
        xadj[v+1]=length(adjncy)+1
    end
    g=_GPGraph(xadj,adjncy,ones(Int32,length(adjncy)),ones(Int32,n))
    order=Int32[]
    _gp_nd_recurse!(g,collect(Int32,1:n),order)
    return Int.(order)
end

function _gp_nd_recurse!(g::_GPGraph,verts::Vector{Int32},order::Vector{Int32})
    length(verts)<=4 && (append!(order,verts);return)
    index=zeros(Int32,_gp_nvertices(g))
    for (i,v) in enumerate(verts);index[v]=i;end
    xadj=Vector{Int32}(undef,length(verts)+1);xadj[1]=1
    adjncy=Int32[];adjwgt=Int32[]
    for (i,v) in enumerate(verts)
        for e in g.xadj[v]:g.xadj[v+1]-1
            j=index[g.adjncy[e]]
            j==0 && continue
            push!(adjncy,j);push!(adjwgt,g.adjwgt[e])
        end
        xadj[i+1]=length(adjncy)+1
    end
    sub=_GPGraph(xadj,adjncy,adjwgt,ones(Int32,length(verts)))
    isempty(adjncy) && (append!(order,verts);return)
    total=length(verts)
    target=total÷2
    side=_gp_bisect_ml(sub,target,max(1,total*30÷1000),3,1)
    # Vertex separator = the smaller side's boundary vertices.
    bndA=Int32[];bndB=Int32[]
    for i in 1:length(verts)
        boundary=false
        for e in sub.xadj[i]:sub.xadj[i+1]-1
            side[sub.adjncy[e]]!=side[i] && (boundary=true;break)
        end
        boundary || continue
        side[i] ? push!(bndA,i) : push!(bndB,i)
    end
    sep=length(bndA)<=length(bndB) ? bndA : bndB
    insep=falses(length(verts));for i in sep;insep[i]=true;end
    left=Int32[verts[i] for i in 1:length(verts) if side[i] && !insep[i]]
    right=Int32[verts[i] for i in 1:length(verts) if !side[i] && !insep[i]]
    isempty(left)||isempty(right) &&
        (append!(order,left);append!(order,right);
         append!(order,Int32[verts[i] for i in sep]);return)
    _gp_nd_recurse!(g,left,order)
    _gp_nd_recurse!(g,right,order)
    append!(order,Int32[verts[i] for i in sep])
    return order
end
