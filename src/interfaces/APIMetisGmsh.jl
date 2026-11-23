# Gmsh 4.15.2's default mesh partitioner (Mesh.MetisAlgorithm = 1), ported from
# src/mesh/meshPartition.cpp (makeGraph, Graph::createDualGraph(false),
# fillDefaultWeights, partitionGraph and its topology correction) and from the
# METIS 5 copy in contrib/metis that Gmsh builds: METIS_PartGraphRecursive with
# the PMETIS settings of SetupCtrl (rtype FM, ctype SHEM, iptype GROW, ncuts 1,
# niter 10, ufactor 1, CoarsenTo 20, no2hop 0, objtype CUT), single constraint.
#
# Build facts reproduced here:
# - idx_t is 64-bit (HAVE_64BIT_SIZE_T), real_t is 32-bit float; every
#   float/double promotion follows the C expression it comes from.
# - USE_GKRAND is not defined, so GKlib draws from the C runtime's rand(); on
#   the Windows build that is the MSVCRT linear congruential generator, seeded
#   by InitRandom(-1) -> srand(4321). irand() (idx_t) is gk_randint64:
#   `(uint64_t)rand() << 32 | (uint64_t)rand()`. The order in which the two
#   rand() calls of that expression run is unspecified in C; the order the
#   Windows build uses is `_MG_RAND_LEFT_FIRST` (checked against Gmsh).
# - iabs is labs, which is 32-bit on Windows; every value passed to it here is
#   far below 2^31, so plain abs is equivalent.
#
# Internally vertex ids are 1-based; values whose numeric value matters to
# METIS (random positions, the 2-hop hash keys, the 0-initialised mark array,
# bndptr positions) keep their C 0-based meaning.

# ---------------------------------------------------------------------------
# C runtime rand() and GKlib's idx_t random routines
# ---------------------------------------------------------------------------

const _MG_RAND_LEFT_FIRST=Ref(true)

mutable struct _MGRand
    state::UInt32
end

@inline function _mg_crand!(r::_MGRand)
    r.state=r.state*UInt32(214013)+UInt32(2531011)
    return Int64((r.state>>16)&0x7fff)
end

@inline function _mg_irand!(r::_MGRand)
    a=_mg_crand!(r); b=_mg_crand!(r)
    return _MG_RAND_LEFT_FIRST[] ? (a<<32)|b : (b<<32)|a
end

@inline _mg_irand_in_range!(r::_MGRand,max::Int64)=_mg_irand!(r)%max

# irandArrayPermute(n, p, nshuffles, 1); p holds 1-based values.
function _mg_rand_permute!(r::_MGRand,p::Vector{Int64},n::Int64,nshuffles::Int64)
    @inbounds for i in 1:n
        p[i]=i
    end
    if n<10
        @inbounds for _ in 1:n
            v=_mg_irand_in_range!(r,n)+1
            u=_mg_irand_in_range!(r,n)+1
            p[v],p[u]=p[u],p[v]
        end
    else
        @inbounds for _ in 1:nshuffles
            v=_mg_irand_in_range!(r,n-3)+1
            u=_mg_irand_in_range!(r,n-3)+1
            p[v],p[u+2]=p[u+2],p[v]
            p[v+1],p[u+3]=p[u+3],p[v+1]
            p[v+2],p[u]=p[u],p[v+2]
            p[v+3],p[u+1]=p[u+1],p[v+3]
        end
    end
    return p
end

# ---------------------------------------------------------------------------
# GKlib max priority queue (gk_mkpqueue.h, rpq: real_t keys, KEY_LT = '>')
# ---------------------------------------------------------------------------

mutable struct _MGPQ
    nnodes::Int
    key::Vector{Float32}
    val::Vector{Int}
    locator::Vector{Int}
end
_MGPQ(n::Int)=_MGPQ(0,Vector{Float32}(undef,n),Vector{Int}(undef,n),fill(-1,n))

function _mg_pq_reset!(q::_MGPQ)
    @inbounds for i in q.nnodes-1:-1:0
        q.locator[q.val[i+1]]=-1
    end
    q.nnodes=0
end

function _mg_pq_insert!(q::_MGPQ,node::Int,key::Float32)
    heap_k=q.key; heap_v=q.val; loc=q.locator
    i=q.nnodes; q.nnodes+=1
    @inbounds while i>0
        j=(i-1)>>1
        if key>heap_k[j+1]
            heap_k[i+1]=heap_k[j+1]; heap_v[i+1]=heap_v[j+1]
            loc[heap_v[i+1]]=i
            i=j
        else
            break
        end
    end
    @inbounds heap_k[i+1]=key; heap_v[i+1]=node; loc[node]=i
    return nothing
end

@inline function _mg_pq_sift_down!(q::_MGPQ,i::Int,newkey::Float32,nnodes::Int)
    heap_k=q.key; heap_v=q.val; loc=q.locator
    @inbounds while (j=2i+1)<nnodes
        if heap_k[j+1]>newkey
            if j+1<nnodes && heap_k[j+2]>heap_k[j+1]
                j+=1
            end
            heap_k[i+1]=heap_k[j+1]; heap_v[i+1]=heap_v[j+1]; loc[heap_v[i+1]]=i
            i=j
        elseif j+1<nnodes && heap_k[j+2]>newkey
            j+=1
            heap_k[i+1]=heap_k[j+1]; heap_v[i+1]=heap_v[j+1]; loc[heap_v[i+1]]=i
            i=j
        else
            break
        end
    end
    return i
end

@inline function _mg_pq_sift_up!(q::_MGPQ,i::Int,newkey::Float32)
    heap_k=q.key; heap_v=q.val; loc=q.locator
    @inbounds while i>0
        j=(i-1)>>1
        if newkey>heap_k[j+1]
            heap_k[i+1]=heap_k[j+1]; heap_v[i+1]=heap_v[j+1]; loc[heap_v[i+1]]=i
            i=j
        else
            break
        end
    end
    return i
end

function _mg_pq_delete!(q::_MGPQ,node::Int)
    loc=q.locator
    @inbounds i=loc[node]
    @inbounds loc[node]=-1
    q.nnodes-=1
    @inbounds if q.nnodes>0 && q.val[q.nnodes+1]!=node
        node=q.val[q.nnodes+1]
        newkey=q.key[q.nnodes+1]
        oldkey=q.key[i+1]
        i=newkey>oldkey ? _mg_pq_sift_up!(q,i,newkey) : _mg_pq_sift_down!(q,i,newkey,q.nnodes)
        q.key[i+1]=newkey; q.val[i+1]=node; loc[node]=i
    end
    return nothing
end

function _mg_pq_update!(q::_MGPQ,node::Int,newkey::Float32)
    @inbounds i=q.locator[node]
    @inbounds oldkey=q.key[i+1]
    i=newkey>oldkey ? _mg_pq_sift_up!(q,i,newkey) : _mg_pq_sift_down!(q,i,newkey,q.nnodes)
    @inbounds q.key[i+1]=newkey; q.val[i+1]=node; q.locator[node]=i
    return nothing
end

function _mg_pq_gettop!(q::_MGPQ)
    q.nnodes==0 && return -1
    q.nnodes-=1
    @inbounds vtx=q.val[1]
    @inbounds q.locator[vtx]=-1
    if (i=q.nnodes)>0
        @inbounds key=q.key[i+1]; node=q.val[i+1]
        i=_mg_pq_sift_down!(q,0,key,q.nnodes)
        @inbounds q.key[i+1]=key; q.val[i+1]=node; q.locator[node]=i
    end
    return vtx
end

# ---------------------------------------------------------------------------
# GKlib's GK_MKQSORT (glibc-derived quicksort + insertion sort), as ikvsorti
# ---------------------------------------------------------------------------

function _mg_ikvsorti!(keys::Vector{Int64},vals::Vector{Int64},n::Int)
    n==0 && return
    @inline lt(a,b)=@inbounds keys[a]<keys[b]
    @inline function sw(a,b)
        @inbounds keys[a],keys[b]=keys[b],keys[a]
        @inbounds vals[a],vals[b]=vals[b],vals[a]
    end
    TH=4
    if n>TH
        lo=1; hi=n
        slo=Int[0]; shi=Int[0]; top=1          # one dummy entry, as `_top = _stack + 1`
        while top>0
            mid=lo+((hi-lo)>>1)
            lt(mid,lo) && sw(mid,lo)
            if lt(hi,mid)
                sw(mid,hi)
                lt(mid,lo) && sw(mid,lo)
            end
            left=lo+1; right=hi-1
            while true
                while lt(left,mid); left+=1; end
                while lt(mid,right); right-=1; end
                if left<right
                    sw(left,right)
                    if mid==left
                        mid=right
                    elseif mid==right
                        mid=left
                    end
                    left+=1; right-=1
                elseif left==right
                    left+=1; right-=1
                    break
                end
                left<=right || break
            end
            if right-lo<=TH
                if hi-left<=TH
                    top-=1
                    if top>0
                        lo=slo[top+1]; hi=shi[top+1]
                    end
                else
                    lo=left
                end
            elseif hi-left<=TH
                hi=right
            elseif right-lo>hi-left
                top+=1
                if length(slo)<top; push!(slo,0); push!(shi,0); end
                slo[top]=lo; shi[top]=right
                lo=left
            else
                top+=1
                if length(slo)<top; push!(slo,0); push!(shi,0); end
                slo[top]=left; shi[top]=hi
                hi=right
            end
        end
    end
    endp=n; tmp=1
    thresh=min(1+TH,endp)
    for run in 2:thresh
        lt(run,tmp) && (tmp=run)
    end
    tmp!=1 && sw(tmp,1)
    run=2
    while (run+=1)<=endp
        tmp=run-1
        while lt(run,tmp); tmp-=1; end
        tmp+=1
        if tmp!=run
            hk=keys[run]; hv=vals[run]
            @inbounds for h in run:-1:tmp+1
                keys[h]=keys[h-1]; vals[h]=vals[h-1]
            end
            keys[tmp]=hk; vals[tmp]=hv
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# METIS graph and control structures
# ---------------------------------------------------------------------------

mutable struct _MGGraph
    nvtxs::Int
    nedges::Int
    xadj::Vector{Int}       # C offsets: edges of vertex v are xadj[v]+1:xadj[v+1]
    vwgt::Vector{Int}
    adjncy::Vector{Int}     # 1-based vertex ids
    adjwgt::Vector{Int}
    tvwgt::Int
    invtvwgt::Float32
    label::Vector{Int}
    cmap::Vector{Int}
    where::Vector{Int}      # 0/1
    pwgts::Vector{Int}
    bndptr::Vector{Int}     # -1 or 0-based position in bndind
    bndind::Vector{Int}
    id::Vector{Int}
    ed::Vector{Int}
    mincut::Int
    nbnd::Int
    coarser::Union{Nothing,_MGGraph}
    finer::Union{Nothing,_MGGraph}
end

function _MGGraph(nvtxs,xadj,vwgt,adjncy,adjwgt)
    g=_MGGraph(nvtxs,xadj[nvtxs+1],xadj,vwgt,adjncy,adjwgt,0,0f0,Int[],Int[],
               Int[],Int[],Int[],Int[],Int[],Int[],-1,-1,nothing,nothing)
    _mg_setup_tvwgt!(g)
    return g
end

function _mg_setup_tvwgt!(g::_MGGraph)
    s=0
    @inbounds for i in 1:g.nvtxs; s+=g.vwgt[i]; end
    g.tvwgt=s
    g.invtvwgt=Float32(1.0/(s>0 ? Float64(s) : 1.0))
end

mutable struct _MGCtrl
    rng::_MGRand
    niter::Int
    CoarsenTo::Int
    maxvwgt::Int
    ubfactor::Float32
    pijbm::Vector{Float32}
    no2hop::Bool
end

# ---------------------------------------------------------------------------
# Coarsening (coarsen.c)
# ---------------------------------------------------------------------------

@inline _mg_deg(g::_MGGraph,v::Int)=@inbounds g.xadj[v+1]-g.xadj[v]

function _mg_coarsen(ctrl::_MGCtrl,graph::_MGGraph)
    eqewgts=true
    @inbounds for i in 2:graph.nedges
        if graph.adjwgt[1]!=graph.adjwgt[i]
            eqewgts=false; break
        end
    end
    ctrl.maxvwgt=trunc(Int,1.5*Float64(graph.tvwgt)/Float64(ctrl.CoarsenTo))
    while true
        if eqewgts || graph.nedges==0
            _mg_match!(ctrl,graph,false)
        else
            _mg_match!(ctrl,graph,true)
        end
        graph=graph.coarser
        eqewgts=false
        (graph.nvtxs>ctrl.CoarsenTo && graph.nvtxs<0.85*graph.finer.nvtxs &&
         graph.nedges>graph.nvtxs÷2) || break
    end
    return graph
end

# BucketSortKeysInc
function _mg_bucket_sort_inc(n::Int,max::Int,keys::Vector{Int},tperm::Vector{Int})
    counts=zeros(Int,max+2)
    @inbounds for i in 1:n; counts[keys[i]+1]+=1; end
    @inbounds for i in 2:max+1; counts[i]+=counts[i-1]; end
    @inbounds for i in max+2:-1:2; counts[i]=counts[i-1]; end
    counts[1]=0
    perm=Vector{Int}(undef,n)
    @inbounds for ii in 1:n
        i=tperm[ii]
        perm[counts[keys[i]+1]+1]=i
        counts[keys[i]+1]+=1
    end
    return perm
end

# Match_RM (shem=false) and Match_SHEM (shem=true), single constraint.
function _mg_match!(ctrl::_MGCtrl,graph::_MGGraph,shem::Bool)
    nvtxs=graph.nvtxs
    xadj=graph.xadj; vwgt=graph.vwgt; adjncy=graph.adjncy; adjwgt=graph.adjwgt
    maxvwgt=ctrl.maxvwgt
    match=fill(-1,nvtxs)
    if shem
        tperm=_mg_rand_permute!(ctrl.rng,Vector{Int}(undef,nvtxs),nvtxs,nvtxs÷8)
        avgdegree=trunc(Int,0.7*Float64(xadj[nvtxs+1]÷nvtxs))
        degrees=Vector{Int}(undef,nvtxs)
        @inbounds for i in 1:nvtxs
            d=xadj[i+1]-xadj[i]
            degrees[i]=d>avgdegree ? avgdegree : d
        end
        perm=_mg_bucket_sort_inc(nvtxs,avgdegree,degrees,tperm)
    else
        perm=_mg_rand_permute!(ctrl.rng,Vector{Int}(undef,nvtxs),nvtxs,nvtxs÷8)
    end
    nunmatched=0
    last_unmatched=0
    @inbounds for pi in 0:nvtxs-1
        i=perm[pi+1]
        match[i]==-1 || continue
        maxidx=i
        maxwgt=-1
        if vwgt[i]<maxvwgt
            if xadj[i]==xadj[i+1]
                last_unmatched=max(pi,last_unmatched)+1
                while last_unmatched<nvtxs
                    j=perm[last_unmatched+1]
                    if match[j]==-1
                        maxidx=j; break
                    end
                    last_unmatched+=1
                end
            else
                if shem
                    for j in xadj[i]+1:xadj[i+1]
                        k=adjncy[j]
                        if match[k]==-1 && maxwgt<adjwgt[j] && vwgt[i]+vwgt[k]<=maxvwgt
                            maxidx=k; maxwgt=adjwgt[j]
                        end
                    end
                else
                    for j in xadj[i]+1:xadj[i+1]
                        k=adjncy[j]
                        if match[k]==-1 && vwgt[i]+vwgt[k]<=maxvwgt
                            maxidx=k; break
                        end
                    end
                end
                if maxidx==i && 3*vwgt[i]<maxvwgt
                    nunmatched+=1
                    maxidx=-1
                end
            end
        end
        if maxidx!=-1
            match[i]=maxidx
            match[maxidx]=i
        end
    end
    if !ctrl.no2hop && Float64(nunmatched)>0.10*Float64(nvtxs)
        _mg_match_2hop!(ctrl,graph,perm,match,nunmatched)
    end
    cmap=Vector{Int}(undef,nvtxs)
    cnvtxs=0
    @inbounds for i in 1:nvtxs
        if match[i]==-1
            match[i]=i
            cnvtxs+=1; cmap[i]=cnvtxs
        elseif i<=match[i]
            cnvtxs+=1; cmap[i]=cnvtxs; cmap[match[i]]=cnvtxs
        end
    end
    graph.cmap=cmap
    _mg_create_coarse!(graph,cnvtxs,match)
    return cnvtxs
end

function _mg_match_2hop!(ctrl,graph,perm,match,nunmatched)
    n=graph.nvtxs
    nunmatched=_mg_match_2hop_any!(graph,perm,match,nunmatched,2)
    nunmatched=_mg_match_2hop_all!(graph,perm,match,nunmatched,64)
    if Float64(nunmatched)>1.5*0.10*Float64(n)
        nunmatched=_mg_match_2hop_any!(graph,perm,match,nunmatched,3)
    end
    if Float64(nunmatched)>2.0*0.10*Float64(n)
        nunmatched=_mg_match_2hop_any!(graph,perm,match,nunmatched,n)
    end
    return nunmatched
end

function _mg_match_2hop_any!(graph,perm,match,nunmatched,maxdegree)
    nvtxs=graph.nvtxs; xadj=graph.xadj; adjncy=graph.adjncy
    colptr=zeros(Int,nvtxs+1)
    @inbounds for i in 1:nvtxs
        if match[i]==-1 && xadj[i+1]-xadj[i]<maxdegree
            for j in xadj[i]+1:xadj[i+1]
                colptr[adjncy[j]]+=1
            end
        end
    end
    @inbounds for i in 2:nvtxs; colptr[i]+=colptr[i-1]; end      # MAKECSR(i, nvtxs, colptr)
    @inbounds for i in nvtxs+1:-1:2; colptr[i]=colptr[i-1]; end
    colptr[1]=0
    rowind=Vector{Int}(undef,colptr[nvtxs+1])
    @inbounds for pi in 1:nvtxs
        i=perm[pi]
        if match[i]==-1 && xadj[i+1]-xadj[i]<maxdegree
            for j in xadj[i]+1:xadj[i+1]
                c=adjncy[j]
                rowind[colptr[c]+1]=i
                colptr[c]+=1
            end
        end
    end
    @inbounds for i in nvtxs+1:-1:2; colptr[i]=colptr[i-1]; end  # SHIFTCSR
    colptr[1]=0
    @inbounds for pi in 1:nvtxs
        i=perm[pi]
        colptr[i+1]-colptr[i]<2 && continue
        jj=colptr[i+1]; j=colptr[i]
        while j<jj
            if match[rowind[j+1]]==-1
                jj-=1
                while jj>j
                    if match[rowind[jj+1]]==-1
                        match[rowind[j+1]]=rowind[jj+1]
                        match[rowind[jj+1]]=rowind[j+1]
                        nunmatched-=2
                        break
                    end
                    jj-=1
                end
            end
            j+=1
        end
    end
    return nunmatched
end

function _mg_match_2hop_all!(graph,perm,match,nunmatched,maxdegree)
    nvtxs=graph.nvtxs; xadj=graph.xadj; adjncy=graph.adjncy
    mask=typemax(Int64)÷maxdegree
    keys=Int64[]; vals=Int64[]
    @inbounds for pi in 1:nvtxs
        i=perm[pi]
        idegree=xadj[i+1]-xadj[i]
        if match[i]==-1 && idegree>1 && idegree<maxdegree
            k=0
            for j in xadj[i]+1:xadj[i+1]
                k+=(adjncy[j]-1)%mask
            end
            push!(vals,i)
            push!(keys,(k%mask)*maxdegree+idegree)
        end
    end
    ncand=length(keys)
    _mg_ikvsorti!(keys,vals,ncand)
    mark=zeros(Int,nvtxs)                    # holds 0-based ids, as in C
    @inbounds for pi in 1:ncand
        i=vals[pi]
        match[i]==-1 || continue
        for j in xadj[i]+1:xadj[i+1]
            mark[adjncy[j]]=i-1
        end
        for pk in pi+1:ncand
            k=vals[pk]
            match[k]==-1 || continue
            keys[pi]!=keys[pk] && break
            xadj[i+1]-xadj[i]!=xadj[k+1]-xadj[k] && break
            jj=xadj[k]+1
            while jj<=xadj[k+1]
                mark[adjncy[jj]]!=i-1 && break
                jj+=1
            end
            if jj==xadj[k+1]+1
                match[i]=k; match[k]=i
                nunmatched-=2
                break
            end
        end
    end
    return nunmatched
end

# CreateCoarseGraph / CreateCoarseGraphNoMask (identical adjacency order).
function _mg_create_coarse!(graph::_MGGraph,cnvtxs::Int,match::Vector{Int})
    nvtxs=graph.nvtxs; xadj=graph.xadj; vwgt=graph.vwgt
    adjncy=graph.adjncy; adjwgt=graph.adjwgt; cmap=graph.cmap
    cxadj=Vector{Int}(undef,cnvtxs+1); cvwgt=Vector{Int}(undef,cnvtxs)
    cadjncy=Vector{Int}(undef,graph.nedges); cadjwgt=Vector{Int}(undef,graph.nedges)
    htable=fill(-1,cnvtxs)
    cxadj[1]=0; c=0; cnedges=0
    @inbounds for v in 1:nvtxs
        u=match[v]
        u<v && continue
        c+=1
        cvwgt[c]=vwgt[v]
        nedges=0
        for j in xadj[v]+1:xadj[v+1]
            k=cmap[adjncy[j]]
            m=htable[k]
            if m==-1
                cadjncy[cnedges+nedges+1]=k; cadjwgt[cnedges+nedges+1]=adjwgt[j]
                htable[k]=nedges; nedges+=1
            else
                cadjwgt[cnedges+m+1]+=adjwgt[j]
            end
        end
        if v!=u
            cvwgt[c]+=vwgt[u]
            for j in xadj[u]+1:xadj[u+1]
                k=cmap[adjncy[j]]
                m=htable[k]
                if m==-1
                    cadjncy[cnedges+nedges+1]=k; cadjwgt[cnedges+nedges+1]=adjwgt[j]
                    htable[k]=nedges; nedges+=1
                else
                    cadjwgt[cnedges+m+1]+=adjwgt[j]
                end
            end
            if (j=htable[c])!=-1
                nedges-=1
                cadjncy[cnedges+j+1]=cadjncy[cnedges+nedges+1]
                cadjwgt[cnedges+j+1]=cadjwgt[cnedges+nedges+1]
                htable[c]=-1
            end
        end
        for j in 1:nedges
            htable[cadjncy[cnedges+j]]=-1
        end
        cnedges+=nedges
        cxadj[c+1]=cnedges
    end
    resize!(cadjncy,cnedges); resize!(cadjwgt,cnedges)
    cg=_MGGraph(cnvtxs,cxadj,cvwgt,cadjncy,cadjwgt)
    cg.finer=graph
    graph.coarser=cg
    return cg
end

# ---------------------------------------------------------------------------
# 2-way partition parameters, balance and FM refinement
# ---------------------------------------------------------------------------

@inline function _mg_bnd_insert!(g::_MGGraph,nbnd::Int,i::Int)
    @inbounds g.bndind[nbnd+1]=i
    @inbounds g.bndptr[i]=nbnd
    return nbnd+1
end
@inline function _mg_bnd_delete!(g::_MGGraph,nbnd::Int,i::Int)
    bndind=g.bndind; bndptr=g.bndptr
    nbnd-=1
    @inbounds bndind[bndptr[i]+1]=bndind[nbnd+1]
    @inbounds bndptr[bndind[nbnd+1]]=bndptr[i]
    @inbounds bndptr[i]=-1
    return nbnd
end

function _mg_alloc_2way!(g::_MGGraph)
    n=g.nvtxs
    g.pwgts=zeros(Int,2); g.where=Vector{Int}(undef,n); g.bndptr=Vector{Int}(undef,n)
    g.bndind=Vector{Int}(undef,n); g.id=Vector{Int}(undef,n); g.ed=Vector{Int}(undef,n)
end

function _mg_compute_2way_params!(g::_MGGraph)
    n=g.nvtxs; xadj=g.xadj; vwgt=g.vwgt; adjncy=g.adjncy; adjwgt=g.adjwgt
    where=g.where; id=g.id; ed=g.ed
    g.pwgts[1]=0; g.pwgts[2]=0
    fill!(g.bndptr,-1)
    @inbounds for i in 1:n
        g.pwgts[where[i]+1]+=vwgt[i]
    end
    nbnd=0; mincut=0
    @inbounds for i in 1:n
        istart=xadj[i]; iend=xadj[i+1]
        me=where[i]; tid=0; ted=0
        for j in istart+1:iend
            if me==where[adjncy[j]]
                tid+=adjwgt[j]
            else
                ted+=adjwgt[j]
            end
        end
        id[i]=tid; ed[i]=ted
        if ted>0 || istart==iend
            nbnd=_mg_bnd_insert!(g,nbnd,i)
            mincut+=ted
        end
    end
    g.mincut=mincut÷2
    g.nbnd=nbnd
end

function _mg_load_imbalance_diff(g::_MGGraph,pijbm::Vector{Float32},ubf::Float32)
    mx=-1.0f0
    @inbounds for j in 1:2
        cur=Float32(g.pwgts[j])*pijbm[j]-ubf
        cur>mx && (mx=cur)
    end
    return mx
end

function _mg_balance2way!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32})
    _mg_load_imbalance_diff(g,ctrl.pijbm,ctrl.ubfactor)<=0 && return
    d=ntp[1]*Float32(g.tvwgt)-Float32(g.pwgts[1])
    abs(trunc(Int,d))<3*g.tvwgt÷g.nvtxs && return
    if g.nbnd>0
        _mg_bnd2way_balance!(ctrl,g,ntp)
    else
        _mg_general2way_balance!(ctrl,g,ntp)
    end
end

function _mg_bnd2way_balance!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32})
    n=g.nvtxs; xadj=g.xadj; vwgt=g.vwgt; adjncy=g.adjncy; adjwgt=g.adjwgt
    where=g.where; id=g.id; ed=g.ed; pwgts=g.pwgts; bndptr=g.bndptr; bndind=g.bndind
    tpw0=trunc(Int,Float32(g.tvwgt)*ntp[1])
    tpw=(tpw0,g.tvwgt-tpw0)
    mindiff=abs(tpw0-pwgts[1])
    from=pwgts[1]<tpw0 ? 1 : 0
    to=(from+1)%2
    queue=_MGPQ(n)
    moved=fill(-1,n)
    nbnd=g.nbnd
    perm=_mg_rand_permute!(ctrl.rng,Vector{Int}(undef,max(nbnd,1)),nbnd,nbnd÷5)
    @inbounds for ii in 1:nbnd
        v=bndind[perm[ii]]
        if where[v]==from && vwgt[v]<=mindiff
            _mg_pq_insert!(queue,v,Float32(ed[v]-id[v]))
        end
    end
    mincut=g.mincut
    @inbounds for nswaps in 0:n-1
        (higain=_mg_pq_gettop!(queue))==-1 && break
        pwgts[to+1]+vwgt[higain]>tpw[to+1] && break
        mincut-=(ed[higain]-id[higain])
        pwgts[to+1]+=vwgt[higain]; pwgts[from+1]-=vwgt[higain]
        where[higain]=to
        moved[higain]=nswaps
        id[higain],ed[higain]=ed[higain],id[higain]
        if ed[higain]==0 && xadj[higain]<xadj[higain+1]
            nbnd=_mg_bnd_delete!(g,nbnd,higain)
        end
        for j in xadj[higain]+1:xadj[higain+1]
            k=adjncy[j]
            kwgt=to==where[k] ? adjwgt[j] : -adjwgt[j]
            id[k]+=kwgt; ed[k]-=kwgt
            if bndptr[k]!=-1
                if ed[k]==0
                    nbnd=_mg_bnd_delete!(g,nbnd,k)
                    if moved[k]==-1 && where[k]==from && vwgt[k]<=mindiff
                        _mg_pq_delete!(queue,k)
                    end
                else
                    if moved[k]==-1 && where[k]==from && vwgt[k]<=mindiff
                        _mg_pq_update!(queue,k,Float32(ed[k]-id[k]))
                    end
                end
            else
                if ed[k]>0
                    nbnd=_mg_bnd_insert!(g,nbnd,k)
                    if moved[k]==-1 && where[k]==from && vwgt[k]<=mindiff
                        _mg_pq_insert!(queue,k,Float32(ed[k]-id[k]))
                    end
                end
            end
        end
    end
    g.mincut=mincut
    g.nbnd=nbnd
end

function _mg_general2way_balance!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32})
    n=g.nvtxs; xadj=g.xadj; vwgt=g.vwgt; adjncy=g.adjncy; adjwgt=g.adjwgt
    where=g.where; id=g.id; ed=g.ed; pwgts=g.pwgts; bndptr=g.bndptr
    tpw0=trunc(Int,Float32(g.tvwgt)*ntp[1])
    tpw=(tpw0,g.tvwgt-tpw0)
    mindiff=abs(tpw0-pwgts[1])
    from=pwgts[1]<tpw0 ? 1 : 0
    to=(from+1)%2
    queue=_MGPQ(n)
    moved=fill(-1,n)
    perm=_mg_rand_permute!(ctrl.rng,Vector{Int}(undef,n),n,n÷5)
    @inbounds for ii in 1:n
        i=perm[ii]
        if where[i]==from && vwgt[i]<=mindiff
            _mg_pq_insert!(queue,i,Float32(ed[i]-id[i]))
        end
    end
    mincut=g.mincut
    nbnd=g.nbnd
    @inbounds for nswaps in 0:n-1
        (higain=_mg_pq_gettop!(queue))==-1 && break
        pwgts[to+1]+vwgt[higain]>tpw[to+1] && break
        mincut-=(ed[higain]-id[higain])
        pwgts[to+1]+=vwgt[higain]; pwgts[from+1]-=vwgt[higain]
        where[higain]=to
        moved[higain]=nswaps
        id[higain],ed[higain]=ed[higain],id[higain]
        if ed[higain]==0 && bndptr[higain]!=-1 && xadj[higain]<xadj[higain+1]
            nbnd=_mg_bnd_delete!(g,nbnd,higain)
        end
        if ed[higain]>0 && bndptr[higain]==-1
            nbnd=_mg_bnd_insert!(g,nbnd,higain)
        end
        for j in xadj[higain]+1:xadj[higain+1]
            k=adjncy[j]
            kwgt=to==where[k] ? adjwgt[j] : -adjwgt[j]
            id[k]+=kwgt; ed[k]-=kwgt
            if moved[k]==-1 && where[k]==from && vwgt[k]<=mindiff
                _mg_pq_update!(queue,k,Float32(ed[k]-id[k]))
            end
            if ed[k]==0 && bndptr[k]!=-1
                nbnd=_mg_bnd_delete!(g,nbnd,k)
            elseif ed[k]>0 && bndptr[k]==-1
                nbnd=_mg_bnd_insert!(g,nbnd,k)
            end
        end
    end
    g.mincut=mincut
    g.nbnd=nbnd
end

# FM_2WayCutRefine
function _mg_fm2way!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32},niter::Int)
    n=g.nvtxs; xadj=g.xadj; vwgt=g.vwgt; adjncy=g.adjncy; adjwgt=g.adjwgt
    where=g.where; id=g.id; ed=g.ed; pwgts=g.pwgts; bndptr=g.bndptr; bndind=g.bndind
    moved=fill(-1,n); swaps=Vector{Int}(undef,n); perm=Vector{Int}(undef,n)
    tpw0=trunc(Int,Float32(g.tvwgt)*ntp[1])
    tpw1=g.tvwgt-tpw0
    limit=trunc(Int,min(max(0.01*Float64(n),15.0),100.0))
    avgvwgt=min((pwgts[1]+pwgts[2])÷20,2*(pwgts[1]+pwgts[2])÷n)
    queues=(_MGPQ(n),_MGPQ(n))
    origdiff=abs(tpw0-pwgts[1])
    @inbounds for pass in 0:niter-1
        _mg_pq_reset!(queues[1]); _mg_pq_reset!(queues[2])
        mincutorder=-1
        newcut=mincut=initcut=g.mincut
        mindiff=abs(tpw0-pwgts[1])
        nbnd=g.nbnd
        _mg_rand_permute!(ctrl.rng,perm,nbnd,nbnd)
        for ii in 1:nbnd
            v=bndind[perm[ii]]
            _mg_pq_insert!(queues[where[v]+1],v,Float32(ed[v]-id[v]))
        end
        nswaps=0
        while nswaps<n
            from=tpw0-pwgts[1]<tpw1-pwgts[2] ? 0 : 1
            to=(from+1)%2
            (higain=_mg_pq_gettop!(queues[from+1]))==-1 && break
            newcut-=(ed[higain]-id[higain])
            pwgts[to+1]+=vwgt[higain]; pwgts[from+1]-=vwgt[higain]
            if (newcut<mincut && abs(tpw0-pwgts[1])<=origdiff+avgvwgt) ||
               (newcut==mincut && abs(tpw0-pwgts[1])<mindiff)
                mincut=newcut
                mindiff=abs(tpw0-pwgts[1])
                mincutorder=nswaps
            elseif nswaps-mincutorder>limit
                newcut+=(ed[higain]-id[higain])
                pwgts[from+1]+=vwgt[higain]; pwgts[to+1]-=vwgt[higain]
                break
            end
            where[higain]=to
            moved[higain]=nswaps
            swaps[nswaps+1]=higain
            id[higain],ed[higain]=ed[higain],id[higain]
            if ed[higain]==0 && xadj[higain]<xadj[higain+1]
                nbnd=_mg_bnd_delete!(g,nbnd,higain)
            end
            for j in xadj[higain]+1:xadj[higain+1]
                k=adjncy[j]
                kwgt=to==where[k] ? adjwgt[j] : -adjwgt[j]
                id[k]+=kwgt; ed[k]-=kwgt
                if bndptr[k]!=-1
                    if ed[k]==0
                        nbnd=_mg_bnd_delete!(g,nbnd,k)
                        moved[k]==-1 && _mg_pq_delete!(queues[where[k]+1],k)
                    else
                        moved[k]==-1 && _mg_pq_update!(queues[where[k]+1],k,Float32(ed[k]-id[k]))
                    end
                else
                    if ed[k]>0
                        nbnd=_mg_bnd_insert!(g,nbnd,k)
                        moved[k]==-1 && _mg_pq_insert!(queues[where[k]+1],k,Float32(ed[k]-id[k]))
                    end
                end
            end
            nswaps+=1
        end
        for i in 1:nswaps
            moved[swaps[i]]=-1
        end
        nswaps-=1
        while nswaps>mincutorder
            higain=swaps[nswaps+1]
            to=where[higain]=(where[higain]+1)%2
            id[higain],ed[higain]=ed[higain],id[higain]
            if ed[higain]==0 && bndptr[higain]!=-1 && xadj[higain]<xadj[higain+1]
                nbnd=_mg_bnd_delete!(g,nbnd,higain)
            elseif ed[higain]>0 && bndptr[higain]==-1
                nbnd=_mg_bnd_insert!(g,nbnd,higain)
            end
            pwgts[to+1]+=vwgt[higain]; pwgts[(to+1)%2+1]-=vwgt[higain]
            for j in xadj[higain]+1:xadj[higain+1]
                k=adjncy[j]
                kwgt=to==where[k] ? adjwgt[j] : -adjwgt[j]
                id[k]+=kwgt; ed[k]-=kwgt
                if bndptr[k]!=-1 && ed[k]==0
                    nbnd=_mg_bnd_delete!(g,nbnd,k)
                end
                if bndptr[k]==-1 && ed[k]>0
                    nbnd=_mg_bnd_insert!(g,nbnd,k)
                end
            end
            nswaps-=1
        end
        g.mincut=mincut
        g.nbnd=nbnd
        (mincutorder<=0 || mincut==initcut) && break
    end
end

# ---------------------------------------------------------------------------
# Initial partition (initpart.c)
# ---------------------------------------------------------------------------

function _mg_init2way!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32},niparts::Int)
    if g.nedges==0
        _mg_random_bisection!(ctrl,g,ntp,niparts)
    else
        _mg_grow_bisection!(ctrl,g,ntp,niparts)
    end
end

function _mg_random_bisection!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32},niparts::Int)
    n=g.nvtxs; vwgt=g.vwgt
    _mg_alloc_2way!(g)
    where=g.where
    bestwhere=Vector{Int}(undef,n); perm=Vector{Int}(undef,n)
    zeromaxpwgt=trunc(Int,ctrl.ubfactor*Float32(g.tvwgt)*ntp[1])
    bestcut=0
    for inbfs in 0:niparts-1
        fill!(where,1)
        if inbfs>0
            _mg_rand_permute!(ctrl.rng,perm,n,n÷2)
            pw1=g.tvwgt; pw0=0
            @inbounds for ii in 1:n
                i=perm[ii]
                if pw0+vwgt[i]<zeromaxpwgt
                    where[i]=0
                    pw0+=vwgt[i]; pw1-=vwgt[i]
                    pw0>zeromaxpwgt && break
                end
            end
        end
        _mg_compute_2way_params!(g)
        _mg_balance2way!(ctrl,g,ntp)
        _mg_fm2way!(ctrl,g,ntp,4)
        if inbfs==0 || bestcut>g.mincut
            bestcut=g.mincut
            copyto!(bestwhere,where)
            bestcut==0 && break
        end
    end
    g.mincut=bestcut
    copyto!(where,bestwhere)
end

function _mg_grow_bisection!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32},niparts::Int)
    n=g.nvtxs; xadj=g.xadj; vwgt=g.vwgt; adjncy=g.adjncy
    _mg_alloc_2way!(g)
    where=g.where
    bestwhere=Vector{Int}(undef,n); queue=Vector{Int}(undef,n); touched=Vector{Int}(undef,n)
    ubf=ctrl.ubfactor
    onemaxpwgt=trunc(Int,ubf*Float32(g.tvwgt)*ntp[2])
    oneminpwgt=trunc(Int,(1.0/Float64(ubf))*Float64(g.tvwgt)*Float64(ntp[2]))
    rng=ctrl.rng
    bestcut=0
    for inbfs in 0:niparts-1
        fill!(where,1); fill!(touched,0)
        pw1=g.tvwgt; pw0=0
        queue[1]=_mg_irand_in_range!(rng,n)+1
        touched[queue[1]]=1
        first=0; last=1; nleft=n-1; drain=false
        @inbounds while true
            if first==last
                (nleft==0 || drain) && break
                k=_mg_irand_in_range!(rng,nleft)
                i=1
                while i<=n
                    if touched[i]==0
                        k==0 && break
                        k-=1
                    end
                    i+=1
                end
                queue[1]=i; touched[i]=1
                first=0; last=1
                nleft-=1
            end
            i=queue[first+1]; first+=1
            if pw0>0 && pw1-vwgt[i]<oneminpwgt
                drain=true
                continue
            end
            where[i]=0
            pw0+=vwgt[i]; pw1-=vwgt[i]
            pw1<=onemaxpwgt && break
            drain=false
            for j in xadj[i]+1:xadj[i+1]
                k=adjncy[j]
                if touched[k]==0
                    last+=1; queue[last]=k
                    touched[k]=1
                    nleft-=1
                end
            end
        end
        pw1==0 && (where[_mg_irand_in_range!(rng,n)+1]=1)
        pw0==0 && (where[_mg_irand_in_range!(rng,n)+1]=0)
        _mg_compute_2way_params!(g)
        _mg_balance2way!(ctrl,g,ntp)
        _mg_fm2way!(ctrl,g,ntp,ctrl.niter)
        if inbfs==0 || bestcut>g.mincut
            bestcut=g.mincut
            copyto!(bestwhere,where)
            bestcut==0 && break
        end
    end
    g.mincut=bestcut
    copyto!(where,bestwhere)
end

# ---------------------------------------------------------------------------
# Uncoarsening (refine.c)
# ---------------------------------------------------------------------------

function _mg_project2way!(g::_MGGraph)
    cg=g.coarser
    _mg_alloc_2way!(g)
    cwhere=cg.where; cbndptr=cg.bndptr
    n=g.nvtxs; cmap=g.cmap; xadj=g.xadj; adjncy=g.adjncy; adjwgt=g.adjwgt
    where=g.where; id=g.id; ed=g.ed
    fill!(g.bndptr,-1)
    @inbounds for i in 1:n
        j=cmap[i]
        where[i]=cwhere[j]
        cmap[i]=cbndptr[j]
    end
    nbnd=0
    @inbounds for i in 1:n
        istart=xadj[i]; iend=xadj[i+1]
        tid=0; ted=0
        if cmap[i]==-1
            for j in istart+1:iend; tid+=adjwgt[j]; end
        else
            me=where[i]
            for j in istart+1:iend
                if me==where[adjncy[j]]
                    tid+=adjwgt[j]
                else
                    ted+=adjwgt[j]
                end
            end
        end
        id[i]=tid; ed[i]=ted
        if ted>0 || istart==iend
            nbnd=_mg_bnd_insert!(g,nbnd,i)
        end
    end
    g.mincut=cg.mincut
    g.nbnd=nbnd
    g.pwgts[1]=cg.pwgts[1]; g.pwgts[2]=cg.pwgts[2]
    g.coarser=nothing
end

function _mg_refine2way!(ctrl::_MGCtrl,orggraph::_MGGraph,g::_MGGraph,ntp::NTuple{2,Float32})
    _mg_compute_2way_params!(g)
    while true
        _mg_balance2way!(ctrl,g,ntp)
        _mg_fm2way!(ctrl,g,ntp,ctrl.niter)
        g===orggraph && break
        g=g.finer
        _mg_project2way!(g)
    end
end

# ---------------------------------------------------------------------------
# Recursive bisection driver (pmetis.c)
# ---------------------------------------------------------------------------

function _mg_multilevel_bisect!(ctrl::_MGCtrl,g::_MGGraph,ntp::NTuple{2,Float32})
    ctrl.pijbm[1]=g.invtvwgt/ntp[1]
    ctrl.pijbm[2]=g.invtvwgt/ntp[2]
    cg=_mg_coarsen(ctrl,g)
    niparts=cg.nvtxs<=ctrl.CoarsenTo ? 5 : 7
    _mg_init2way!(ctrl,cg,ntp,niparts)
    _mg_refine2way!(ctrl,g,cg,ntp)
    return g.mincut
end

function _mg_split!(g::_MGGraph)
    n=g.nvtxs; xadj=g.xadj; vwgt=g.vwgt; adjncy=g.adjncy; adjwgt=g.adjwgt
    label=g.label; where=g.where; bndptr=g.bndptr
    rename=Vector{Int}(undef,n)
    snv=[0,0]
    @inbounds for i in 1:n
        k=where[i]+1
        snv[k]+=1; rename[i]=snv[k]
    end
    sx=(Int[0],Int[0]); sw=(Int[],Int[]); sa=(Int[],Int[]); saw=(Int[],Int[]); sl=(Int[],Int[])
    @inbounds for i in 1:n
        p=where[i]+1
        if bndptr[i]==-1
            for j in xadj[i]+1:xadj[i+1]
                push!(sa[p],adjncy[j]); push!(saw[p],adjwgt[j])
            end
        else
            for j in xadj[i]+1:xadj[i+1]
                k=adjncy[j]
                if where[k]==p-1
                    push!(sa[p],k); push!(saw[p],adjwgt[j])
                end
            end
        end
        push!(sw[p],vwgt[i])
        push!(sl[p],label[i])
        push!(sx[p],length(sa[p]))
    end
    out=ntuple(2) do p
        a=sa[p]
        @inbounds for i in eachindex(a); a[i]=rename[a[i]]; end
        sg=_MGGraph(snv[p],sx[p],sw[p],a,saw[p])
        sg.label=sl[p]
        sg
    end
    return out
end

function _mg_rb!(ctrl::_MGCtrl,g::_MGGraph,nparts::Int,part::Vector{Int},
                 tpwgts::Vector{Float32},toff::Int,fpart::Int)
    g.nvtxs==0 && return 0
    half=nparts>>1
    s=0f0
    for k in 1:half; s+=tpwgts[toff+k]; end
    ntp=(s,Float32(1.0-Float64(s)))
    objval=_mg_multilevel_bisect!(ctrl,g,ntp)
    @inbounds for i in 1:g.nvtxs
        part[g.label[i]]=g.where[i]+fpart
    end
    lg=rg=g
    nparts>2 && ((lg,rg)=_mg_split!(g))
    wsum=0f0
    for k in 1:half; wsum+=tpwgts[toff+k]; end
    a=Float32(1.0/Float64(wsum))
    for k in 1:half; tpwgts[toff+k]*=a; end
    b=Float32(1.0/(1.0-Float64(wsum)))
    for k in half+1:nparts; tpwgts[toff+k]*=b; end
    if nparts>3
        objval+=_mg_rb!(ctrl,lg,half,part,tpwgts,toff,fpart)
        objval+=_mg_rb!(ctrl,rg,nparts-half,part,tpwgts,toff+half,fpart+half)
    elseif nparts==3
        objval+=_mg_rb!(ctrl,rg,nparts-half,part,tpwgts,toff+half,fpart+half)
    end
    return objval
end

"""
    _mg_part_graph_recursive(xadj, adjncy, vwgt, nparts) -> Vector{Int64}

METIS_PartGraphRecursive as Gmsh calls it (single constraint, unit edge weights,
uniform target weights, default ubvec, C numbering). `xadj`/`adjncy` are the
0-based CSR arrays; the result holds 0-based parts.
"""
function _mg_part_graph_recursive(xadj::Vector{Int64},adjncy::Vector{Int64},
                                  vwgt::Vector{Int64},nparts::Int)::Vector{Int64}
    nvtxs=length(xadj)-1
    rng=_MGRand(UInt32(4321))
    ubf=Float32(Float64(Float32(1.0+0.001*1))+0.0000499)
    ctrl=_MGCtrl(rng,10,20,0,ubf,zeros(Float32,2),false)
    tpwgts=fill(Float32(1.0/nparts),nparts)
    g=_MGGraph(nvtxs,copy(xadj),copy(vwgt),adjncy.+1,ones(Int,length(adjncy)))
    g.label=collect(1:nvtxs)
    part=zeros(Int,nvtxs)
    _mg_rb!(ctrl,g,nparts,part,tpwgts,0,0)
    return part
end

# ---------------------------------------------------------------------------
# Gmsh's element graph (meshPartition.cpp)
# ---------------------------------------------------------------------------

# Families: 0 point, 1 line, 2 triangle, 3 quadrangle, 4 tetrahedron,
# 5 hexahedron, 6 prism, 7 pyramid, 8 trihedron.
function _mg_family(etype::Integer)
    name=lowercase(first(mesh.get_element_properties(etype)))
    for (p,f) in (("point",0),("line",1),("triangle",2),("quadrilateral",3),
                  ("tetrahedron",4),("hexahedron",5),("prism",6),("pyramid",7),
                  ("trihedron",8))
        startswith(name,p) && return f
    end
    error("Gmsh's partitioner graph does not include element type $etype ($name)")
end
const _MG_NPRIMARY=(1,2,3,4,4,8,6,5,4)
const _MG_DIM=(0,1,2,2,3,3,3,3,3)

# MElement::numCommonNodesInDualGraph(this = fi, other = fj)
function _mg_common_nodes(fi::Int,fj::Int)
    fi==0 && return 1
    fj==0 && return 1
    fi==1 && return fj==1 ? 1 : 2
    fj==1 && return 2
    if fi==2 || fi==3
        (fj==2 || fj==3) && return 2
        return fi==2 ? 3 : 4
    end
    fi==4 && return fj==3 ? 4 : 3
    fi==5 && return fj==2 ? 3 : 4
    (fi==6 || fi==7) && return (fj==3 || fj==5) ? 4 : 3
    return (fj==3 || fj==5) ? 4 : (fj==8 ? 2 : 3)          # trihedron
end

"""
    _mg_gmsh_partition(nparts) -> Dict{Int,Int}

Partition the current session mesh as `gmsh.model.mesh.partition(nparts)` does
with the default options, and return each element tag's 1-based partition.
"""
function _mg_gmsh_partition(nparts::Int)
    # makeGraph(model, graph, -1): volumes (tet, hex, prism, pyramid, trihedron),
    # surfaces (triangle, quadrangle), curves, points; entities by tag.
    order_in_dim=((0,),(1,),(2,3),(4,5,6,7,8))
    etags=Int[]; fams=Int[]; eptr=Int[0]; eind=Int[]
    nodeid=Dict{Int,Int}()
    for d in 3:-1:0
        for (_,tag) in model.get_entities(d)
            types,tags,ntags=mesh.get_elements(d,tag)
            fam=[_mg_family(t) for t in types]
            for f in order_in_dim[d+1]
                for b in eachindex(types)
                    fam[b]==f || continue
                    np=_MG_NPRIMARY[f+1]
                    nn=length(ntags[b])÷max(length(tags[b]),1)
                    for e in eachindex(tags[b])
                        push!(etags,Int(tags[b][e])); push!(fams,f)
                        for q in 1:np
                            nt=Int(ntags[b][(e-1)*nn+q])
                            push!(eind,get!(nodeid,nt,length(nodeid)+1))
                        end
                        push!(eptr,length(eind))
                    end
                end
            end
        end
    end
    return _mg_partition_graph_elements(etags,fams,eptr,eind,length(nodeid),nparts)
end

# makeGraph's element list (tags, families, CSR primary nodes) → Gmsh's
# partitions, METIS plus the topology correction.
function _mg_partition_graph_elements(etags,fams,eptr,eind,nnodes::Int,nparts::Int)
    ne=length(etags)
    ne==0 && return Dict{Int,Int}()
    meshdim=maximum(f->_MG_DIM[f+1],fams)
    nparts<=1 && return Dict(t=>1 for t in etags)
    # Graph::createDualGraph(false)
    nn=nnodes
    nptr=zeros(Int,nn+1)
    for j in 1:eptr[end]; nptr[eind[j]+1]+=1; end
    for i in 2:nn+1; nptr[i]+=nptr[i-1]; end
    nind=Vector{Int}(undef,eptr[end])
    fillp=copy(nptr)
    for i in 1:ne, j in eptr[i]+1:eptr[i+1]
        n=eind[j]; fillp[n]+=1; nind[fillp[n]]=i
    end
    xadj=zeros(Int64,ne+1); adj=Int64[]
    marker=zeros(Int,ne); nbrs=Int[]
    for i in 1:ne
        empty!(nbrs)
        for j in eptr[i]+1:eptr[i+1]
            n=eind[j]
            for k in nptr[n]+1:nptr[n+1]
                e=nind[k]
                e==i && continue
                marker[e]==0 && push!(nbrs,e)
                marker[e]+=1
            end
        end
        for e in nbrs
            marker[e]>=_mg_common_nodes(fams[i],fams[e]) && push!(adj,e-1)
            marker[e]=0
        end
        xadj[i+1]=length(adj)
    end
    # fillDefaultWeights with the default -1 weights
    vwgt=Int64[_MG_DIM[f+1]==meshdim ? 1 : 0 for f in fams]
    epart=_mg_part_graph_recursive(xadj,adj,vwgt,nparts)
    # Check and correct the topology
    for i in 1:3, j in 1:ne
        dj=_MG_DIM[fams[j]+1]
        dj==meshdim && continue
        for k in xadj[j]+1:xadj[j+1]
            e=adj[k]+1
            if dj==_MG_DIM[fams[e]+1]-i && epart[j]!=epart[e]
                epart[j]=epart[e]
                break
            end
        end
    end
    return Dict(etags[i]=>Int(epart[i])+1 for i in 1:ne)
end

# The partitioner settings this port covers: Gmsh's defaults (recursive
# bisection, sorted heavy-edge matching, default imbalance, min-cut).
function _mg_gmsh_default_options()
    get(OPTIONS,"Mesh.MetisAlgorithm",1.0)==1 || return false
    get(OPTIONS,"Mesh.MetisEdgeMatching",2.0)==2 || return false
    get(OPTIONS,"Mesh.MetisMaxLoadImbalance",-1.0)<0 || return false
    get(OPTIONS,"Mesh.MetisObjective",1.0)==1 || return false
    return true
end

# The same partition for a classified model mesh (a `.geo` PartitionMesh):
# makeGraph's order from the classification's cell owners, returned as the
# dense assignment (`_cache_catalog` element order) the producers consume.
function _mg_gmsh_partition_assignment(mesh,class,nparts::Int)
    order_in_dim=((0,),(1,),(2,3),(4,5,6,7,8))
    catalog=_cache_catalog(mesh,class)
    total=sum((size(cells,2) for (_,_,_,cells,_) in catalog);init=0)
    etags=Int[]; fams=Int[]; eptr=Int[0]; eind=Int[]
    nodeid=Dict{Int,Int}()
    for d in 3:-1:0
        owners=Set{Int32}()
        for (msh,dim,_,_,own) in catalog
            dim==d && union!(owners,own)
        end
        for entity in sort!(collect(owners)), f in order_in_dim[d+1]
            for (msh,dim,offset,cells,own) in catalog
                (dim==d && _mg_family(msh)==f) || continue
                np=_MG_NPRIMARY[f+1]
                for column in axes(cells,2)
                    own[column]==entity || continue
                    push!(etags,offset+column); push!(fams,f)
                    for q in 1:np
                        push!(eind,get!(nodeid,Int(cells[q,column]),length(nodeid)+1))
                    end
                    push!(eptr,length(eind))
                end
            end
        end
    end
    parts=_mg_partition_graph_elements(etags,fams,eptr,eind,length(nodeid),nparts)
    assignment=zeros(Int32,total)
    for (index,part) in parts; assignment[index]=Int32(part); end
    return assignment
end
