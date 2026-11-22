# Gmsh 4.15.2's partition entities and partition topology, ported from
# src/mesh/meshPartition.cpp: createNewEntities, divideNonConnectedEntities,
# createPartitionTopology (with assignNewEntityBRep, getReferenceElement,
# assignPartitionBoundary, assignBrep) and assignMeshVertices.
#
# Entity tags, element numbers and node ownership follow Gmsh's construction
# order. Gmsh walks its boundary maps (std::unordered_map keyed by MEdge, MFace
# and MVertex*) in the iteration order of libstdc++'s hashtable, which decides
# the tags of interface entities and the numbers of their new elements; that
# order is reproduced by `_GPHashOrder` (GCC 11's prime rehash policy, the
# 32-bit `__prime_list` of the Windows build, FNV-1a hashes of MEdgeHash and
# MFaceHash, MVertexPtrHash = node number). On a build whose standard library
# differs (libc++ on macOS) Gmsh itself numbers these entities differently.
#
# Supported cells: first-order points, lines, triangles, quadrangles and
# tetrahedra; anything else falls back to the earlier producers.

const _GP_PRIMES=UInt64[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,103,109,
    113,127,137,139,149,157,167,179,193,199,211,227,241,257,277,293,313,337,359,383,409,439,467,503,541,577,
    619,661,709,761,823,887,953,1031,1109,1193,1289,1381,1493,1613,1741,1879,2029,2179,2357,2549,2753,2971,
    3209,3469,3739,4027,4349,4703,5087,5503,5953,6427,6949,7517,8123,8783,9497,10273,11113,12011,12983,14033,
    15173,16411,17749,19183,20753,22447,24281,26267,28411,30727,33223,35933,38873,42043,45481,49201,53201,
    57557,62233,67307,72817,78779,85229,92203,99733,107897,116731,126271,136607,147793,159871,172933,187091,
    202409,218971,236897,256279,277261,299951,324503,351061,379787,410857,444487,480881,520241,562841,608903,
    658753,712697,771049,834181,902483,976369,1056323,1142821,1236397,1337629,1447153,1565659,1693859,1832561,
    1982627,2144977,2320627,2510653,2716249,2938679,3179303,3439651,3721303,4026031,4355707,4712381,5098259,
    5515729,5967347,6456007,6984629,7556579,8175383,8844859,9569143,10352717,11200489,12117689,13109983,
    14183539,15345007,16601593,17961079,19431899,21023161,22744717,24607243,26622317,28802401,31160981,
    33712729,36473443,39460231,42691603,46187573,49969847,54061849,58488943,63278561,68460391,74066549,
    80131819,86693767,93793069,101473717,109783337,118773397,128499677,139022417,150406843,162723577,
    176048909,190465427,206062531,222936881,241193053,260944219,282312799,305431229,330442829,357502601,
    386778277,418451333,452718089,489790921,529899637,573292817,620239453,671030513,725980837,785430967,
    849749479,919334987,994618837,1076067617,1164186217,1259520799,1362662261,1474249943,1594975441,
    1725587117,1866894511,2019773507,2185171673,2364114217,2557710269,2767159799,2993761039,3238918481,
    3504151727,3791104843,4101556399,4294967291]
const _GP_FAST_BKT=(2,2,2,3,5,5,7,7,11,11,11,11,13,13)

# libstdc++ `_Prime_rehash_policy::_M_next_bkt` (max load factor 1).
function _gp_next_bkt(n::UInt64)
    n<UInt64(length(_GP_FAST_BKT)) && return n==0 ? (UInt64(1),UInt64(0)) :
        (UInt64(_GP_FAST_BKT[n+1]),UInt64(_GP_FAST_BKT[n+1]))
    last=length(_GP_PRIMES)-1      # lower_bound over [primes+6, last prime)
    i=searchsortedfirst(Base.view(_GP_PRIMES,7:last),n)+6
    i>=last+1 && return (_GP_PRIMES[end],typemax(UInt64))
    return (_GP_PRIMES[i],_GP_PRIMES[i])
end

# The iteration order of a libstdc++ unordered_map with unique keys: a single
# forward list threaded through buckets, a new key inserted at the front of
# its bucket (or of the whole list when the bucket is empty), rehashing by
# `_M_rehash_aux`.
mutable struct _GPHashOrder{K}
    keys::Vector{K}
    codes::Vector{UInt64}
    next::Vector{Int}                # 0 ends the list
    buckets::Vector{Int}             # -1 empty, 0 before-begin, else node
    head::Int
    next_resize::UInt64
    index::Dict{K,Int}
end
_GPHashOrder{K}() where {K}=_GPHashOrder{K}(K[],UInt64[],Int[],Int[-1],0,UInt64(0),Dict{K,Int}())

_gp_next_of(m::_GPHashOrder,p::Int)=p==0 ? m.head : m.next[p]
_gp_set_next!(m::_GPHashOrder,p::Int,v::Int)=(p==0 ? (m.head=v) : (m.next[p]=v);nothing)

function _gp_rehash!(m::_GPHashOrder,n::UInt64)
    buckets=fill(-1,Int(n))
    p=m.head;m.head=0;bbegin=0
    while p!=0
        nxt=m.next[p];bkt=Int(m.codes[p]%n)+1
        if buckets[bkt]==-1
            m.next[p]=m.head;m.head=p;buckets[bkt]=0
            m.next[p]!=0 && (buckets[bbegin]=p)
            bbegin=bkt
        else
            m.next[p]=_gp_next_of(m,buckets[bkt]);_gp_set_next!(m,buckets[bkt],p)
        end
        p=nxt
    end
    m.buckets=buckets
    return nothing
end

# `operator[]`: the node of `key`, inserted if absent.
function _gp_node!(m::_GPHashOrder{K},key::K,code::UInt64) where {K}
    found=get(m.index,key,0)
    found!=0 && return found
    nbkt=UInt64(length(m.buckets));nelt=UInt64(length(m.keys))
    if nelt+1>m.next_resize
        minb=Float64(max(nelt+1,m.next_resize==0 ? UInt64(11) : UInt64(0)))
        if minb>=nbkt
            target=max(UInt64(floor(minb))+1,nbkt*2)
            n,resize=_gp_next_bkt(target)
            m.next_resize=resize
            _gp_rehash!(m,n)
        else
            m.next_resize=nbkt
        end
    end
    push!(m.keys,key);push!(m.codes,code);push!(m.next,0)
    node=length(m.keys);m.index[key]=node
    n=UInt64(length(m.buckets));bkt=Int(code%n)+1
    if m.buckets[bkt]!=-1
        before=m.buckets[bkt]
        m.next[node]=_gp_next_of(m,before);_gp_set_next!(m,before,node)
    else
        m.next[node]=m.head;m.head=node
        if m.next[node]!=0
            m.buckets[Int(m.codes[m.next[node]]%n)+1]=node
        end
        m.buckets[bkt]=0
    end
    return node
end

function _gp_iteration(m::_GPHashOrder)
    order=Int[];p=m.head
    while p!=0;push!(order,p);p=m.next[p];end
    return order
end

# Gmsh's FNV-1a over the bytes of a size_t array (64-bit constants).
function _gp_fnv1a(values::NTuple{N,UInt64}) where {N}
    h=0xcbf29ce484222325
    for v in values,b in 0:7
        h=(h⊻((v>>(8b))&0xff))*0x00000100000001b3
    end
    return h
end

# ── cells ────────────────────────────────────────────────────────────────────
# Gmsh TYPE_* parents: :pnt :lin :tri :qua :tet.
const _GP_TYPE=Dict(15=>:pnt,1=>:lin,2=>:tri,3=>:qua,4=>:tet)
const _GP_NPRIM=Dict(:pnt=>1,:lin=>2,:tri=>3,:qua=>4,:tet=>4)
const _GP_EDGES=Dict(:lin=>((1,2),),:tri=>((1,2),(2,3),(3,1)),:qua=>((1,2),(2,3),(3,4),(4,1)),
    :tet=>((1,2),(2,3),(3,1),(4,1),(4,3),(4,2)),:pnt=>())
const _GP_FACES=Dict(:tri=>((1,2,3),),:qua=>((1,2,3,4),),
    :tet=>((1,3,2),(1,2,4),(1,4,3),(4,2,3)),:lin=>(),:pnt=>())
_gp_dim(t::Symbol)=t===:pnt ? 0 : t===:lin ? 1 : (t===:tri || t===:qua) ? 2 : 3
# Element storage order inside one entity (GFace: triangles, quadrangles;
# GRegion: tetrahedra first).
_gp_rank(t::Symbol)=t===:qua ? 2 : 1

function _gp_common(a::Symbol,b::Symbol)
    a===:pnt && return 1
    a===:lin && return (b===:lin || b===:pnt) ? 1 : 2
    if a===:tri || a===:qua
        b===:pnt && return 1
        (b===:lin || b===:tri || b===:qua) && return 2
        return a===:tri ? 3 : 4
    end
    # tetrahedron
    b===:pnt && return 1
    b===:lin && return 2
    b===:qua && return 4
    return 3
end

mutable struct _GPCells
    type::Vector{Symbol}
    verts::Vector{Vector{Int}}      # node rows
    num::Vector{UInt64}             # element number
end

mutable struct _GPEnt
    dim::Int
    tag::Int
    parent::Tuple{Int,Int}
    partitions::Vector{Int}
    partition::Bool
    cells::Vector{Int}
    bnd::Vector{Int}                # GFace l_edges / GRegion l_faces (tags)
    dirs::Vector{Int}
    v0::Int                          # GEdge begin / end vertex tags (0 none)
    v1::Int
    up::Vector{Int}                  # vertex edges / edge faces / face regions
    alive::Bool
end
_GPEnt(dim,tag,parent,partitions,partition)=_GPEnt(dim,tag,parent,partitions,partition,Int[],
    Int[],Int[],0,0,Int[],true)

struct _GPModel
    ents::Vector{Dict{Int,_GPEnt}}   # by dim+1
    cells::_GPCells
    owner::Vector{Tuple{Int,Int}}    # cell -> (dim, tag)
    vnum::Vector{UInt64}             # node row -> node number
    counter::Base.RefValue{UInt64}   # MElement global number
end

_gp_sorted(g::_GPModel,dim)=sort!([e for e in values(g.ents[dim+1]) if e.alive];by=e->e.tag)
_gp_maxtag(g::_GPModel,dim)=maximum((e.tag for e in values(g.ents[dim+1]) if e.alive);init=0)
_gp_ordered(g::_GPModel,e::_GPEnt)=length(e.cells)<2 ? e.cells :
    sort(e.cells;by=c->_gp_rank(g.cells.type[c]),alg=Base.Sort.MergeSort)
function _gp_add_cell!(g::_GPModel,e::_GPEnt,c::Int)
    push!(e.cells,c);g.owner[c]=(e.dim,e.tag);return nothing
end
function _gp_new_cell!(g::_GPModel,type::Symbol,verts::Vector{Int})
    push!(g.cells.type,type);push!(g.cells.verts,verts)
    g.counter[]+=1;push!(g.cells.num,g.counter[]);push!(g.owner,(-1,-1))
    return length(g.cells.type)
end
_gp_prim(g::_GPModel,c::Int)=Base.view(g.cells.verts[c],1:_GP_NPRIM[g.cells.type[c]])
function _gp_new_entity!(g::_GPModel,dim,tag,parent,partitions)
    e=_GPEnt(dim,tag,parent,partitions,true);g.ents[dim+1][tag]=e;return e
end
_gp_ent(g::_GPModel,key)=g.ents[key[1]+1][key[2]]

# Graph::createDualGraph(false) over `cells` (fillElementsToNodesMap order).
function _gp_dual(g::_GPModel,cells::Vector{Int})
    ne=length(cells);local_index=Dict{Int,Int}()
    eptr=Vector{Int}(undef,ne+1);eptr[1]=0;eind=Int[]
    for (i,c) in enumerate(cells)
        for row in _gp_prim(g,c)
            push!(eind,get!(local_index,row,length(local_index)))
        end
        eptr[i+1]=length(eind)
    end
    nn=length(local_index)
    incident=[Int[] for _ in 1:nn]
    for i in 1:ne,j in eptr[i]+1:eptr[i+1]
        push!(incident[eind[j]+1],i-1)
    end
    marker=zeros(Int,ne);nbrs=Int[]
    adj=Vector{Vector{Int}}(undef,ne)
    for i in 1:ne
        empty!(nbrs)
        for j in eptr[i]+1:eptr[i+1],k in incident[eind[j]+1]
            k==i-1 && continue
            marker[k+1]==0 && push!(nbrs,k)
            marker[k+1]+=1
        end
        list=Int[]
        for k in nbrs
            marker[k+1]>=_gp_common(g.cells.type[cells[i]],g.cells.type[cells[k+1]]) && push!(list,k)
            marker[k+1]=0
        end
        adj[i]=list
    end
    return adj
end

# fillConnectedElements, with Gmsh's use of 0 as the visited mark.
function _gp_components(g::_GPModel,cells::Vector{Int},adj::Vector{Vector{Int}})
    ne=length(cells);ne<2 && return Vector{Vector{Int}}()
    adj=[copy(a) for a in adj]
    components=Vector{Vector{Int}}();start=0;size=0;isolated=0
    while true
        stack=[start];members=Set{Int}([start])
        while !isempty(stack)
            top=pop!(stack);push!(members,top)
            for (slot,k) in enumerate(adj[top+1])
                if k!=0;push!(stack,k);adj[top+1][slot]=0;end
            end
        end
        component=sort!([cells[i+1] for i in members];by=c->g.cells.num[c])
        push!(components,component);size+=length(members)
        size==ne && break
        start=0
        for i in 0:ne-1
            if any(!=(0),adj[i+1]);start=i;break;end
        end
        if start==0
            skip=0
            for i in 1:ne-1
                if isempty(adj[i+1])
                    if skip==isolated;start=i;isolated+=1;break;end
                    skip+=1
                end
            end
        end
    end
    return components
end

function _gp_connected_by_count(adj::Vector{Vector{Int}})
    n=UInt64(length(adj));edges=UInt64(sum(length,adj;init=0)÷2)
    return ((n-1)*(n-2))÷2<edges
end

# divideNonConnectedEntities(model, dim, ...).
function _gp_divide!(g::_GPModel,dim::Int)
    if dim<0 || dim==0
        number=_gp_maxtag(g,0)
        for v in _gp_sorted(g,0)
            v.partition && length(v.cells)>1 || continue
            for c in copy(v.cells)
                pv=_gp_new_entity!(g,0,number+=1,v.parent,copy(v.partitions))
                _gp_add_cell!(g,pv,c)
                for etag in v.up
                    e=get(g.ents[2],etag,nothing);(e===nothing || !e.alive) && continue
                    if e.v0==v.tag;e.v0=pv.tag;_gp_push_up!(pv,etag);end
                    if e.v1==v.tag;e.v1=pv.tag;_gp_push_up!(pv,etag);end
                end
            end
            v.alive=false;empty!(v.cells)
        end
    end
    for d in (1,2,3)
        (dim<0 || dim==d) || continue
        number=_gp_maxtag(g,d)
        for e in _gp_sorted(g,d)
            e.partition || continue
            cells=_gp_ordered(g,e)
            adj=_gp_dual(g,cells)
            _gp_connected_by_count(adj) && continue
            components=_gp_components(g,cells,adj)
            length(components)>1 || continue
            uppers=d<3 ? [u for u in e.up if haskey(g.ents[d+2],u) && g.ents[d+2][u].alive] : Int[]
            old=Int[]
            for u in uppers
                push!(old,_gp_del_bnd!(_gp_ent(g,(d+1,u)),e.tag))
            end
            for component in components
                pe=_gp_new_entity!(g,d,number+=1,e.parent,copy(e.partitions))
                for c in component;_gp_add_cell!(g,pe,c);end
                for (j,u) in enumerate(uppers)
                    up=_gp_ent(g,(d+1,u));push!(up.bnd,pe.tag);push!(up.dirs,old[j])
                    _gp_push_up!(pe,u)
                end
            end
            e.alive=false;empty!(e.cells)
        end
    end
    return nothing
end
_gp_push_up!(e::_GPEnt,u::Int)=(u in e.up || push!(e.up,u);nothing)
function _gp_del_bnd!(e::_GPEnt,tag::Int)
    i=findfirst(==(tag),e.bnd);i===nothing && return 0
    deleteat!(e.bnd,i);isempty(e.dirs) && return 0
    o=e.dirs[i];deleteat!(e.dirs,i);return o
end

# computeOrientation(reference, element).
function _gp_orientation(g::_GPModel,reference::Int,element::Int)
    t=g.cells.type[element];ev=_gp_prim(g,element);rt=g.cells.type[reference]
    rv=_gp_prim(g,reference)
    if _gp_dim(t)==2
        key=sort(collect(ev))
        for f in _GP_FACES[rt]
            fv=[rv[k] for k in f]
            sort(fv)==key && return fv==collect(ev) ? 1 : -1
        end
    elseif _gp_dim(t)==1
        key=minmax(ev[1],ev[2])
        for (a,b) in _GP_EDGES[rt]
            minmax(rv[a],rv[b])==key && return (rv[a],rv[b])==(ev[1],ev[2]) ? 1 : -1
        end
    elseif _gp_dim(t)==0
        rv[1]==ev[1] && return 1
        length(rv)>1 && rv[2]==ev[1] && return -1
    end
    return 0
end

# GFace::setEdge / GRegion::setFace / GEdge::setVertex, with the reverse link.
function _gp_set_bnd!(g::_GPModel,upper::_GPEnt,lower::_GPEnt,ori::Int)
    if upper.dim==1
        ori>0 ? (upper.v0=lower.tag) : ori<0 ? (upper.v1=lower.tag) : nothing
    else
        push!(upper.bnd,lower.tag);push!(upper.dirs,ori)
    end
    if upper.dim==3
        # GFace::addRegion keeps two regions.
        length(lower.up)<2 && !(upper.tag in lower.up) && push!(lower.up,upper.tag)
    else
        _gp_push_up!(lower,upper.tag)
    end
    return nothing
end

# getReferenceElement.
function _gp_reference(g::_GPModel,pairs::Vector{Tuple{Int,Vector{Int}}})
    minimum_size=minimum(length(p[2]) for p in pairs)
    candidates=[(p[1],copy(p[2])) for p in pairs if length(p[2])==minimum_size]
    if length(candidates)==length(pairs) && all(p->p[2]==candidates[1][2],candidates)
        return 0
    end
    while length(candidates)>1
        lowest=typemax(Int)
        for p in candidates
            isempty(p[2]) && return candidates[1][1]
            lowest=min(lowest,p[2][1])
        end
        kept=Tuple{Int,Vector{Int}}[]
        for p in candidates
            p[2][1]==lowest && (popfirst!(p[2]);push!(kept,p))
        end
        candidates=kept
    end
    return candidates[1][1]
end

# One boundary level of createPartitionTopology: `entries[key]` lists
# (element, partitions); `make(reference, key)` returns the new element.
function _gp_boundary_level!(g::_GPModel,dim::Int,map::_GPHashOrder,entries,make)
    number=_gp_maxtag(g,dim)
    created=Dict{Tuple{Vector{Int},Tuple{Int,Int}},Int}()
    for node in _gp_iteration(map)
        pairs=entries[node]
        partitions=sort!(unique!(reduce(vcat,(p[2] for p in pairs);init=Int[])))
        length(partitions)<2 && continue
        reference=_gp_reference(g,pairs)
        reference==0 && continue
        refent=g.owner[reference]
        key=(partitions,refent)
        tag=get(created,key,0);fresh=tag==0
        if fresh
            tag=(number+=1);created[key]=tag
            _gp_new_entity!(g,dim,tag,_gp_ent(g,refent).parent,partitions)
        end
        pe=_gp_ent(g,(dim,tag))
        c=make(reference,map.keys[node]);_gp_add_cell!(g,pe,c)
        if fresh
            first_ref=Dict{Tuple{Int,Int},Int}()
            for (element,_) in pairs
                haskey(first_ref,g.owner[element]) || (first_ref[g.owner[element]]=element)
            end
            for owner in sort!(collect(keys(first_ref)))
                upper=_gp_ent(g,owner)
                _gp_set_bnd!(g,upper,pe,_gp_orientation(g,first_ref[owner],pe.cells[1]))
            end
        end
    end
    return nothing
end

function _gp_boundary_elements(g::_GPModel,cells::Vector{Int},adj,part::Vector{Int},
                               nparts::Int,dim::Int)
    sets=[Set{Int}() for _ in 1:nparts]
    for i in eachindex(cells),k in adj[i]
        part[i]!=part[k+1] && _gp_dim(g.cells.type[cells[i]])==dim && push!(sets[part[i]+1],cells[i])
    end
    return [sort!(collect(s);by=c->g.cells.num[c]) for s in sets]
end

# makeGraph(model, graph, selectDim) element order.
function _gp_graph_cells(g::_GPModel,dim::Int)
    cells=Int[]
    for e in _gp_sorted(g,dim);append!(cells,_gp_ordered(g,e));end
    return cells
end

function _gp_topology!(g::_GPModel,graph_cells::Vector{Int},graph_adj,graph_part::Vector{Int},
                       nparts::Int,meshdim::Int)
    boundary=_gp_boundary_elements(g,graph_cells,graph_adj,graph_part,nparts,meshdim)
    # assignNewEntityBRep over the full mesh graph.
    seen=Set{Tuple{Tuple{Int,Int},Tuple{Int,Int}}}()
    brep=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int,Int}}}()
    for (i,c) in enumerate(graph_cells),k in graph_adj[i]
        other=graph_cells[k+1]
        _gp_dim(g.cells.type[c])==_gp_dim(g.cells.type[other])+1 || continue
        g1=g.owner[c];g2=g.owner[other]
        (g1,g2) in seen && continue
        push!(seen,(g1,g2))
        push!(get!(()->Tuple{Int,Int,Int}[],brep,g1),(_gp_orientation(g,c,other),g2[1],g2[2]))
    end
    for g1 in sort!(collect(keys(brep)))
        for (ori,d2,t2) in sort!(unique(brep[g1]))
            _gp_set_bnd!(g,_gp_ent(g,g1),_gp_ent(g,(d2,t2)),ori)
        end
    end
    vnum=g.vnum
    if meshdim>=3
        map=_GPHashOrder{NTuple{4,UInt64}}();entries=Vector{Vector{Tuple{Int,Vector{Int}}}}()
        for p in 1:nparts,c in boundary[p]
            v=_gp_prim(g,c)
            for f in _GP_FACES[g.cells.type[c]]
                nums=sort([vnum[v[k]] for k in f])
                key=ntuple(i->i<=length(nums) ? nums[i] : UInt64(0),4)
                node=_gp_node!(map,key,_gp_fnv1a(key))
                node>length(entries) && push!(entries,Tuple{Int,Vector{Int}}[])
                push!(entries[node],(c,[p]))
            end
        end
        make=(reference,key)->begin
            rv=_gp_prim(g,reference)
            for f in _GP_FACES[g.cells.type[reference]]
                nums=sort([vnum[rv[k]] for k in f])
                if ntuple(i->i<=length(nums) ? nums[i] : UInt64(0),4)==key
                    return _gp_new_cell!(g,length(f)==3 ? :tri : :qua,[rv[k] for k in f])
                end
            end
            error("partition face not found on its reference element")
        end
        _gp_boundary_level!(g,2,map,entries,make)
        _gp_divide!(g,2)
    end
    if meshdim>=2
        map=_GPHashOrder{NTuple{2,UInt64}}();entries=Vector{Vector{Tuple{Int,Vector{Int}}}}()
        function add_edges(c,partitions)
            v=_gp_prim(g,c)
            for (a,b) in _GP_EDGES[g.cells.type[c]]
                key=minmax(vnum[v[a]],vnum[v[b]])
                node=_gp_node!(map,key,_gp_fnv1a(key))
                node>length(entries) && push!(entries,Tuple{Int,Vector{Int}}[])
                push!(entries[node],(c,partitions))
            end
        end
        if meshdim==2
            for p in 1:nparts,c in boundary[p];add_edges(c,[p]);end
        else
            cells=_gp_graph_cells(g,2);adj=_gp_dual(g,cells)
            faces=[e for e in _gp_sorted(g,2) if e.partition && !isempty(e.cells)]
            part=Int[];parts=Vector{Vector{Int}}()
            for (i,e) in enumerate(faces)
                push!(parts,e.partitions);append!(part,fill(i-1,length(e.cells)))
            end
            sub=_gp_boundary_elements(g,cells,adj,part,length(faces),2)
            for i in eachindex(faces),c in sub[i];add_edges(c,parts[i]);end
        end
        make=(reference,key)->begin
            rv=_gp_prim(g,reference)
            for (a,b) in _GP_EDGES[g.cells.type[reference]]
                minmax(vnum[rv[a]],vnum[rv[b]])==key && return _gp_new_cell!(g,:lin,[rv[a],rv[b]])
            end
            error("partition edge not found on its reference element")
        end
        _gp_boundary_level!(g,1,map,entries,make)
        _gp_divide!(g,1)
    end
    if meshdim>=1
        map=_GPHashOrder{UInt64}();entries=Vector{Vector{Tuple{Int,Vector{Int}}}}()
        function add_vertices(c,partitions)
            for row in _gp_prim(g,c)
                node=_gp_node!(map,vnum[row],vnum[row])
                node>length(entries) && push!(entries,Tuple{Int,Vector{Int}}[])
                push!(entries[node],(c,partitions))
            end
        end
        if meshdim==1
            for p in 1:nparts,c in boundary[p];add_vertices(c,[p]);end
        else
            cells=_gp_graph_cells(g,1);adj=_gp_dual(g,cells)
            curves=[e for e in _gp_sorted(g,1) if e.partition && !isempty(e.cells)]
            part=Int[];parts=Vector{Vector{Int}}()
            for (i,e) in enumerate(curves)
                push!(parts,e.partitions);append!(part,fill(i-1,length(e.cells)))
            end
            sub=_gp_boundary_elements(g,cells,adj,part,length(curves),1)
            for i in eachindex(curves),c in sub[i];add_vertices(c,parts[i]);end
        end
        rows=Dict(vnum[row]=>row for row in eachindex(vnum))
        make=(reference,key)->_gp_new_cell!(g,:pnt,[rows[key]])
        _gp_boundary_level!(g,0,map,entries,make)
        _gp_divide!(g,0)
    end
    return nothing
end

# assignMeshVertices: each node belongs to the first entity (points, curves,
# surfaces, volumes; by tag; elements in storage order) that uses it.
function _gp_mesh_vertices(g::_GPModel)
    owner=fill((-1,-1),length(g.vnum));order=Int[]
    for dim in 0:3,e in _gp_sorted(g,dim),c in _gp_ordered(g,e),row in g.cells.verts[c]
        owner[row]==(-1,-1) || continue
        owner[row]=(dim,e.tag);push!(order,row)
    end
    return owner,order
end

# The partition plan: Gmsh's PartitionUsingThisSplit (and the entity layer of
# PartitionMesh) for an explicit per-cell assignment in catalog order.
function _partition_gmsh_supported(cached,class)
    class!==nothing && class.public_tags!==nothing || return false
    _high_order_overlay(cached)===nothing || return false
    cached isa MixedMesh && cached.partition_data!==nothing && return false
    for row in _cache_catalog(cached,class)
        isempty(row[4]) && continue
        haskey(_GP_TYPE,Int(row[1])) || return false
    end
    return true
end

function _partition_gmsh_plan(model,cached,class,count,assignment,caller;topology::Bool)
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    public=class.public_tags
    total=sum(size(row[4],2) for row in catalog;init=0)
    length(assignment)==total || throw(ArgumentError(
        "$caller: partition assignment rows do not align with the cell catalog"))
    for (_,record) in Model._discrete_mesh_records_model(model)
        all(uid->haskey(public.node_indices,uid),record.node_tags) &&
        all(uid->haskey(public.element_indices,uid),record.element_tags) || throw(ArgumentError(
            "$caller: partitioning cannot omit caller-owned records outside the cache"))
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    _dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    global_entities,_=_dim01_entities(model,model_entities(model),assembly,eachindex(assembly.nodes))
    # The Gmsh-side model: every entity with its cells in storage order.
    cells=_GPCells(Symbol[],Vector{Int}[],UInt64[])
    owner=Tuple{Int,Int}[];part=Int[]
    maxnum=UInt64(0)
    for (msh,dim,offset,block,owners) in catalog
        t=_GP_TYPE[Int(msh)]
        for column in axes(block,2)
            push!(cells.type,t);push!(cells.verts,Int.(block[:,column]))
            num=public.element_tags[offset+column];push!(cells.num,num);maxnum=max(maxnum,num)
            push!(owner,(dim,Int(owners[column])));push!(part,Int(assignment[offset+column]))
        end
    end
    g=_GPModel([Dict{Int,_GPEnt}() for _ in 1:4],cells,owner,copy(public.node_tags),
        Ref(max(maxnum,UInt64(ELEMENT_TAG_MAX[]))))
    for (dim,tag) in model_entities(model)
        g.ents[dim+1][tag]=_GPEnt(dim,tag,(-1,-1),Int[],false)
    end
    for c in eachindex(owner)
        e=get(g.ents[owner[c][1]+1],owner[c][2],nothing)
        e===nothing && (e=g.ents[owner[c][1]+1][owner[c][2]]=_GPEnt(owner[c][1],owner[c][2],(-1,-1),Int[],false))
        push!(e.cells,c)
    end
    # makeGraph(model, graph, -1): regions, faces, curves, points.
    graph_cells=Int[]
    for dim in (3,2,1,0);append!(graph_cells,_gp_graph_cells(g,dim));end
    meshdim=maximum(_gp_dim(t) for t in cells.type)
    graph_adj=_gp_dual(g,graph_cells)
    graph_part=[part[c]-1 for c in graph_cells]
    # createNewEntities.
    for dim in 0:3
        number=_gp_maxtag(g,dim)
        for e in _gp_sorted(g,dim)
            e.partition && continue
            children=Dict{Int,_GPEnt}()
            for c in _gp_ordered(g,e)
                p=part[c]
                child=get(children,p,nothing)
                child===nothing && (child=children[p]=_gp_new_entity!(g,dim,number+=1,(dim,e.tag),[p]))
                _gp_add_cell!(g,child,c)
            end
            empty!(e.cells)
        end
    end
    if topology
        _gp_divide!(g,-1)
        _gp_topology!(g,graph_cells,graph_adj,graph_part,count,meshdim)
    end
    node_owner,node_order=_gp_mesh_vertices(g)
    return _partition_gmsh_publish(model,cached,class,count,g,node_owner,node_order,
        global_entities,caller)
end

function _partition_gmsh_publish(model,cached,class,count,g::_GPModel,node_owner,node_order,
                                 global_entities,caller)
    public=class.public_tags
    staged=Model._model_mesh_transaction_copy(model)
    # Nodes: Gmsh's getNodes order; nodes no entity uses keep their place
    # after them with their previous owner.
    used=falses(length(g.vnum));for row in node_order;used[row]=true;end
    rows=vcat(node_order,[row for row in eachindex(g.vnum) if !used[row]])
    inverse=zeros(Int32,length(rows));for (i,row) in enumerate(rows);inverse[row]=Int32(i);end
    coords=cached.coords[:,rows]
    node_entities=Tuple{Int,Int32}[used[row] ? (node_owner[row][1],Int32(node_owner[row][2])) :
        class.node_entities[row] for row in rows]
    node_parametric=Vector{Union{Nothing,Vector{Float64}}}(undef,length(rows))
    for (i,row) in enumerate(rows)
        p=class.node_parametric[row]
        # An MVertex keeps its parameters when assignMeshVertices moves it to
        # a partition entity; Tessella's node records hold them only on an
        # entity of the same dimension.
        moved=used[row] && length(something(p,Float64[]))!=node_owner[row][1]
        node_parametric[i]=p===nothing || moved ? nothing : copy(p)
    end
    # Partition entities, their physical groups (the parent's, Gmsh's
    # assignPhysicals without PartitionCreatePhysicals) and boundaries.
    part_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    for dim in 0:3,e in _gp_sorted(g,dim)
        e.partition || continue
        Model.add_discrete_entity!(staged,dim,e.tag)
        groups=e.parent[1]==dim ? Model._model_projection_physical_tags(staged,dim,e.parent[2]) : Int32[]
        staged.entity_physicals[(dim,e.tag)]=copy(groups)
        for group in groups
            members=staged.physical[(dim,Int(group))]
            e.tag in members || push!(members,e.tag)
        end
        boundary=dim==1 ? Int32[t for t in (e.v0,-e.v1) if t!=0] :
            Int32[e.dirs[i]<0 ? -e.bnd[i] : e.bnd[i] for i in eachindex(e.bnd)
                  if (dim==0 ? false : get(g.ents[dim],e.bnd[i],nothing)!==nothing && g.ents[dim][e.bnd[i]].alive)]
        staged.discrete[(dim,e.tag)].boundary=Tuple{Int,Int}[(dim-1,Int(t)) for t in boundary]
        referenced=unique!(Int32[inverse[row] for c in e.cells for row in g.cells.verts[c]])
        part_entities[(dim,e.tag)]=Elements.MixedEntity(dim,e.tag,_partition_bounds(coords,referenced);
            physical_tags=groups,boundaries=boundary,partitions=Int32.(e.partitions),
            parent=(e.parent[1],Int32(e.parent[2])))
    end
    for dim in 0:3
        staged.next_tag[dim+1]=max(staged.next_tag[dim+1],_gp_maxtag(g,dim))
    end
    # Blocks: one per Gmsh type, cells by entity then storage order.
    msh_of=Dict(:pnt=>15,:lin=>1,:tri=>2,:qua=>3,:tet=>4)
    blocks=ElementBlock[];owners=Vector{Int32}[];uids=Vector{UInt64}[]
    for t in (:pnt,:lin,:tri,:qua,:tet)
        selected=Int[];selowner=Int32[]
        for dim in 0:3,e in _gp_sorted(g,dim),c in _gp_ordered(g,e)
            g.cells.type[c]===t || continue
            push!(selected,c);push!(selowner,Int32(e.tag))
        end
        isempty(selected) && continue
        nodes=Matrix{Int32}(undef,_GP_NPRIM[t],length(selected))
        physical=Vector{Int32}(undef,length(selected))
        for (j,c) in enumerate(selected)
            for (k,row) in enumerate(g.cells.verts[c]);nodes[k,j]=inverse[row];end
            entity=get(part_entities,(_gp_dim(t),Int(selowner[j])),nothing)
            physical[j]=entity===nothing || isempty(entity.physical_tags) ? Int32(0) :
                first(entity.physical_tags)
        end
        push!(blocks,ElementBlock(msh_of[t],nodes,physical));push!(owners,selowner)
        push!(uids,g.cells.num[selected])
    end
    data=Elements._native_mixed_entity_data(global_entities;node_entities,node_parametric,
        external_node_tags=public.node_tags[rows],block_entities=owners,external_element_tags=uids)
    partition=Elements.MixedPartitionData(count;entities=part_entities)
    target=MixedMesh(coords,blocks;entity_data=data,elementary_entities=owners,
        physical_names=staged.physical_names,partition_data=partition,
        periodic_links=cached isa MixedMesh ? cached.periodic_links : Elements.MixedPeriodicLink[],
        ancillary_sections=cached isa MixedMesh ? cached.ancillary_sections : [],
        data_sections=cached isa MixedMesh ? cached.data_sections : Elements.MshDataSection[])
    _geo_import_entity_receipts!(staged,target)
    identities=Dict(uid=>uid for uid in public.node_tags)
    declarations=union(Set(keys(global_entities)),Set(keys(part_entities)))
    _dim01_reconcile_records!(staged,declarations,target,identities;runtime_parameters=data.node_parametric)
    staged.meshing.mesh_authority=Model._model_mesh_authority(target;caller,
        private_rank=collect(UInt64,1:nnodes(target)))
    cell_map=Dict(uid=>uid for uid in public.element_tags)
    return (;model=staged,mesh=target,class=_higher_classification_from_data(target,caller),
        authority=Set{Tuple{Int,Int32}}((dim,Int32(tag)) for (dim,tag) in declarations),
        max_node_tag=NODE_TAG_MAX[],max_element_tag=max(ELEMENT_TAG_MAX[],g.counter[]),
        node_tag_map=identities,element_tag_map=cell_map,visibility_tag_map=cell_map)
end
