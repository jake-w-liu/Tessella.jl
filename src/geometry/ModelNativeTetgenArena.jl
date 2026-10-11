module ModelNativeTetgenArena
using ...Predicates: orient3

# Private native interior arena: recursive link and two-to-three flips. Logical identities
# are native point rows and stable pool slots; coordinates are predicates only.
struct Handle
    cell::Int32
    version::UInt8
end
Handle(cell::Integer,version::Integer)=Handle(Int32(cell),UInt8(version))
const NULL=Handle(0,0)
const ORG=(4,4,2,2,3,1,1,3,2,3,4,1)
const DEST=(3,1,1,3,2,3,4,1,4,4,2,2)
const APEX=(2,3,4,1,4,4,2,2,3,1,1,3)
const OPPO=(1,2,3,4,1,2,3,4,1,2,3,4)
const ESYM=(9,6,11,4,3,7,1,5,10,0,8,2)
mutable struct Cell
    vertices::NTuple{4,Int32}
    links::NTuple{4,Handle}
    alive::Bool
end
mutable struct Arena
    points::Vector{NTuple{3,Float64}}
    cells::Vector{Cell}
    nreal::Int
    hints::Vector{Int32}
    dead::Vector{Int32}
    seed::UInt32
    protected_faces::Set{NTuple{3,Int32}}
    protected_segments::Set{Tuple{Int32,Int32}}
    trace::Vector{NamedTuple}
    pool_events::Vector{NamedTuple}
end
struct Unsupported <: Exception
    reason::String
end
Base.showerror(io::IO,error::Unsupported)=print(io,"Native Tetgen arena prototype: ",error.reason)
facekey(face)=Tuple(sort!(Int32[face...]))
edgekey(a,b)=minmax(Int32(a),Int32(b))
@inline vertex(a,h,table)=a.cells[h.cell].vertices[table[Int(h.version)+1]]
@inline org(a,h)=vertex(a,h,ORG)
@inline dest(a,h)=vertex(a,h,DEST)
@inline apex(a,h)=vertex(a,h,APEX)
@inline oppo(a,h)=vertex(a,h,OPPO)
@inline ordered(a,h)=(org(a,h),dest(a,h),apex(a,h),oppo(a,h))
@inline enext(h)=Handle(h.cell,mod(Int(h.version)+4,12))
@inline eprev(h)=Handle(h.cell,mod(Int(h.version)+8,12))
@inline esym(h)=Handle(h.cell,ESYM[Int(h.version)+1])
@inline eprevesym(h)=esym(eprev(h))
@inline eorgoppo(h)=eprev(esym(enext(h)))
@inline edestoppo(h)=enext(esym(eprev(h)))
@inline orient(a,i,j,k,l)=orient3(a.points[i],a.points[j],a.points[k],a.points[l])
function findhandle(a,cell,first,second,third)
    for version in 0:11
        h=Handle(cell,version)
        (org(a,h),dest(a,h),apex(a,h))==(first,second,third) && return h
    end
    throw(Unsupported("inconsistent oriented face incidence"))
end
function bond!(a,first,second)
    (org(a,first),dest(a,first),apex(a,first))==
        (dest(a,second),org(a,second),apex(a,second)) ||
        throw(Unsupported("bond does not reverse its common directed edge"))
    v1=Int(first.version);v2=Int(second.version)
    a.cells[first.cell].links=Base.setindex(a.cells[first.cell].links,
        Handle(second.cell,(v2&3)+mod((v1&12)+(v2&12),12)),(v1&3)+1)
    a.cells[second.cell].links=Base.setindex(a.cells[second.cell].links,
        Handle(first.cell,(v1&3)+mod((v2&12)+(v1&12),12)),(v2&3)+1)
    return nothing
end
function fsym(a,h)
    encoded=a.cells[h.cell].links[(Int(h.version)&3)+1]
    encoded.cell==0 && throw(Unsupported("dummy/hull traversal is outside the interior slice"))
    return Handle(encoded.cell,mod(Int(encoded.version)+12-(Int(h.version)&12),12))
end
@inline fnext(a,h)=fsym(a,esym(h))
function Arena(kernel;facets=NTuple{3,Int32}[],segments=Tuple{Int32,Int32}[],seed=1)
    points=copy(kernel.points)
    cells=Cell[]
    for source in kernel.cells
        source.vertices[1]==0 && continue
        tuple=source.vertices
        sign=orient3((points[node] for node in tuple)...)
        sign>0 && (tuple=(tuple[2],tuple[1],tuple[3],tuple[4]))
        push!(cells,Cell(tuple,(NULL,NULL,NULL,NULL),true))
    end
    a=Arena(points,cells,kernel.count,zeros(Int32,length(points)),Int32[],UInt32(seed),
        Set{NTuple{3,Int32}}(facekey(face) for face in facets),
        Set{Tuple{Int32,Int32}}(edgekey(edge...) for edge in segments),NamedTuple[],NamedTuple[])
    pending=Dict{NTuple{3,Int32},Handle}()
    for cell in eachindex(cells),version in 0:3
        h=Handle(cell,version)
        key=facekey((org(a,h),dest(a,h),apex(a,h)))
        previous=pop!(pending,key,NULL)
        if previous.cell==0
            pending[key]=h
        else
            second=findhandle(a,cell,dest(a,previous),org(a,previous),apex(a,previous))
            bond!(a,previous,second)
        end
    end
    # The real-pool traversal fills point2tet in physical order; native hull
    # slots are outside this prototype and are never chosen as point hints.
    for cell in eachindex(cells),node in cells[cell].vertices
        a.hints[node]=Int32(cell)
    end
    return a
end
function point2tetorg(a,node)
    cell=a.hints[node]
    cell!=0 && a.cells[cell].alive || throw(Unsupported("missing/stale point2tet hint"))
    position=findfirst(==(node),a.cells[cell].vertices)
    position===nothing && throw(Unsupported("point2tet hint no longer owns the node"))
    return Handle(cell,(11,3,7,0)[position])
end
function randomnation!(a,choices::Integer)
    0<choices<=typemax(UInt32) || throw(ArgumentError("randomnation choices"))
    choices=UInt32(choices)
    if choices>=714025
        newrandom=(a.seed*UInt32(1366)+UInt32(150889))%UInt32(714025)
        a.seed=(newrandom*UInt32(1366)+UInt32(150889))%UInt32(714025)
        newrandom=newrandom*(choices÷UInt32(714025))+a.seed
        return newrandom>=choices ? newrandom-choices : newrandom
    end
    a.seed=(a.seed*UInt32(1366)+UInt32(150889))%UInt32(714025)
    return a.seed%choices
end
function finddirection!(a,h,endpoint)
    origin=org(a,h)
    dest(a,h)==endpoint && return h,:across_vertex
    apex(a,h)==endpoint && return eprevesym(h),:across_vertex
    for _ in 1:max(100,20length(a.cells))
        p=oppo(a,h)
        p==endpoint && return enext(esym(h)),:across_vertex
        b=dest(a,h);c=apex(a,h)
        hori=orient(a,origin,b,c,endpoint)
        rori=orient(a,b,origin,p,endpoint)
        lori=orient(a,origin,c,p,endpoint)
        move=if hori>0
            if rori>0
                if lori>0
                    draw=randomnation!(a,3)
                    draw==0 ? :h : draw==1 ? :r : :l
                else
                    randomnation!(a,2)!=0 ? :h : :r
                end
            elseif lori>0
                randomnation!(a,2)!=0 ? :h : :l
            else
                :h
            end
        elseif rori>0
            lori>0 ? (randomnation!(a,2)!=0 ? :r : :l) : :r
        elseif lori>0
            :l
        else
            if hori==0
                rori==0 && return h,:across_vertex
                lori==0 && return eprevesym(h),:across_vertex
                return h,:across_edge
            elseif rori==0
                lori==0 && return enext(esym(h)),:across_vertex
                return enext(esym(h)),:across_edge
            elseif lori==0
                return eprevesym(h),:across_edge
            end
            return h,:across_face
        end
        h=move===:r ? fnext(a,h) : move===:l ? enext(fnext(a,eprev(h))) : enext(fsym(a,h))
        org(a,h)==origin || throw(Unsupported("finddirection changed its fixed origin"))
    end
    throw(Unsupported("finddirection exceeded prototype work limit"))
end
function spin(a,start)
    ring=Handle[];h=start
    while true
        push!(ring,h);h=fnext(a,h)
        h.cell==start.cell && break
        length(ring)<=length(a.cells) || throw(Unsupported("edge ring does not close"))
    end
    return ring
end
function crossing(a,facet,first,last)
    (first in facet || last in facet) && return :none
    p,q,r=facet
    s0=orient(a,p,q,r,first);s1=orient(a,p,q,r,last)
    (s0>0 && s1>0 || s0<0 && s1<0) && return :none
    (s0==0 || s1==0) && throw(Unsupported("touching/coplanar facet intersection requires native classifier"))
    values=(orient(a,first,last,p,q),orient(a,first,last,q,r),orient(a,first,last,r,p))
    all(>=(0),values) || all(<=(0),values) || return :none
    return any(iszero,values) ? :across_edge : :across_face
end
function flip32!(a,ring)
    length(ring)==3 || return false
    first,second,third=ring
    e=org(a,first);d=dest(a,first)
    edgekey(e,d) in a.protected_segments && return false
    p=apex(a,first);q=apex(a,second);r=apex(a,third)
    ordered(a,first)==(e,d,p,q) && ordered(a,second)==(e,d,q,r) &&
        ordered(a,third)==(e,d,r,p) || throw(Unsupported("three-cell star is not in native oriented ring order"))
    orient(a,q,p,r,e)<0 && orient(a,p,q,r,d)<0 || return false
    for h in ring
        facekey((org(a,h),dest(a,h),apex(a,h))) in a.protected_faces && return false
    end
    top=Handle[fsym(a,eorgoppo(h)) for h in ring]
    bottom=Handle[fsym(a,edestoppo(h)) for h in ring]
    before=Tuple(ordered(a,h) for h in ring)
    a.cells[first.cell].vertices=(p,q,r,d)
    a.cells[first.cell].links=(NULL,NULL,NULL,NULL)
    a.cells[second.cell].vertices=(q,p,r,e)
    a.cells[second.cell].links=(NULL,NULL,NULL,NULL)
    a.cells[third.cell].vertices=(0,0,0,0)
    a.cells[third.cell].links=(NULL,NULL,NULL,NULL)
    a.cells[third.cell].alive=false
    previous=isempty(a.dead) ? Int32(0) : last(a.dead)
    push!(a.dead,third.cell)
    push!(a.pool_events,(kind=:deallocate,slot=third.cell,previous_free_top=previous,
        free_count=length(a.dead)))
    h0=Handle(first.cell,11);h1=Handle(second.cell,11)
    bond!(a,h0,h1)
    for index in 1:3
        bond!(a,esym(h0),top[index]);h0=enext(h0)
    end
    for index in 1:3
        bond!(a,esym(h1),bottom[index]);h1=eprev(h1)
    end
    # Native flip32's five setpoint2tet calls, in literal order.
    for node in (p,q,r,d)
        a.hints[node]=first.cell
    end
    a.hints[e]=second.cell
    push!(a.trace,(edge=(e,d),pool=(first.cell,second.cell,third.cell),before=before,
        after=(a.cells[first.cell].vertices,a.cells[second.cell].vertices)))
    return true
end
# Bounded literal flipnm branch: interior reflex-link reduction by three-cell
# flips and interior face 2-to-3. Hull/dummy, rollback and deeper link levels remain explicit
# unsupported branches. Physical pool rows and the parent ring phase survive.
function flip23!(a,first,second,counters)
    p,q,r,d=ordered(a,first);e=oppo(a,second)
    ordered(a,second)==(q,p,r,e) || throw(Unsupported("face23 inputs do not share the native reversed face"))
    facekey((p,q,r)) in a.protected_faces && return nothing
    top=Handle[];bottom=Handle[]
    h=first
    for _ in 1:3
        push!(top,fnext(a,h));h=enext(h)
    end
    h=second
    for _ in 1:3
        push!(bottom,fnext(a,h));h=eprev(h)
    end
    created=if isempty(a.dead)
        push!(a.cells,Cell((0,0,0,0),(NULL,NULL,NULL,NULL),false))
        push!(counters,0)
        slot=Int32(length(a.cells))
        push!(a.pool_events,(kind=:allocate,slot=slot,previous_free_top=Int32(0),free_count=0))
        slot
    else
        slot=pop!(a.dead)
        push!(a.pool_events,(kind=:reuse,slot=slot,previous_free_top=slot,free_count=length(a.dead)))
        slot
    end
    result=Handle[Handle(first.cell,11),Handle(second.cell,11),Handle(created,11)]
    tuples=((e,d,p,q),(e,d,q,r),(e,d,r,p))
    for i in 1:3
        a.cells[result[i].cell].vertices=tuples[i]
        a.cells[result[i].cell].links=(NULL,NULL,NULL,NULL)
        a.cells[result[i].cell].alive=true
        counters[result[i].cell]=0
    end
    for i in 1:3
        bond!(a,esym(result[i]),result[mod(i,3)+1])
    end
    for i in 1:3
        bond!(a,eorgoppo(result[i]),top[i])
    end
    for i in 1:3
        bond!(a,edestoppo(result[i]),bottom[i])
    end
    a.hints[p]=result[1].cell
    a.hints[q]=result[1].cell
    a.hints[r]=result[2].cell
    a.hints[d]=result[1].cell
    a.hints[e]=result[1].cell
    return result
end

function recursive_star!(a,ring,counters,facet;level=0,edgepivot=0,maxlevel=1)
    n=length(ring)
    n>=3 || throw(Unsupported("an edge star has fewer than three cells"))
    pa=org(a,ring[1]);pb=dest(a,ring[1])
    if n==3
        first,second=ring[1],ring[2]
        flip32!(a,ring) || return n
        for h in ring
            counters[h.cell]=0
        end
        ring[1]=Handle(first.cell,11);ring[2]=Handle(second.cell,11)
        return 2
    end
    reflex=0
    for i in 1:n
        previous=mod(i-2,n)+1;following=mod(i,n)+1
        facekey((org(a,ring[i]),dest(a,ring[i]),apex(a,ring[i]))) in a.protected_faces && continue
        (counters[ring[i].cell]>1 || counters[ring[previous].cell]>1) && continue
        pc=apex(a,ring[i]);pd=apex(a,ring[following]);pe=apex(a,ring[previous])
        reducible=false
        if orient(a,pb,pc,pd,pe)>0 && orient(a,pc,pa,pd,pe)>0
            sign=orient(a,pa,pb,pd,pe)
            reducible=sign>0 || (sign==0 && n==4)
        end
        if reducible
            if crossing(a,facet,pe,pd)===:none
                result=flip23!(a,ring[i],fsym(a,ring[i]),counters)
                result===nothing && continue
                replacement=edestoppo(result[1])
                counters[replacement.cell]+=1
                ring[previous]=replacement
                deleteat!(ring,i)
                nn=recursive_star!(a,ring,counters,facet;level,edgepivot,maxlevel)
                nn>2 && sign==0 && throw(Unsupported("flat face23 requires the native rollback branch"))
                return nn
            end
        else
            reflex+=1
        end
    end
    reflex>0 && level<maxlevel || return n
    for i in 1:n
        previous=mod(i-2,n)+1;following=mod(i,n)+1
        (counters[ring[i].cell]>1 || counters[ring[previous].cell]>1) && continue
        pc=apex(a,ring[i]);pd=apex(a,ring[following]);pe=apex(a,ring[previous])
        pivot=0;selected=NULL
        if orient(a,pb,pc,pd,pe)<=0
            selected=enext(ring[i]);pivot=1
        elseif orient(a,pc,pa,pd,pe)<=0
            selected=eprev(ring[i]);pivot=2
        end
        pivot==0 && continue
        edgekey(org(a,selected),dest(a,selected)) in a.protected_segments && continue
        selected=esym(selected)
        child=spin(a,selected)
        sum(counters[h.cell] for h in child)>2 && continue
        for h in child
            counters[h.cell]+=1
        end
        nn=recursive_star!(a,child,counters,facet;level=level+1,edgepivot=pivot,maxlevel)
        if nn==2
            replacement=if pivot==1
                # tmpabtets[0]=[d,a,e,b] => [a,b,e,d].
                enext(esym(enext(child[1])))
            else
                # tmpabtets[1]=[b,d,e,a] => [a,b,e,d].
                eprev(esym(eprev(child[2])))
            end
            (org(a,replacement),dest(a,replacement))==(pa,pb) ||
                throw(Unsupported("reflex-link replacement lost the parent edge phase"))
            counters[replacement.cell]+=1
            ring[previous]=replacement
            deleteat!(ring,i)
            return recursive_star!(a,ring,counters,facet;level,edgepivot,maxlevel)
        end
        for j in 1:nn
            counters[child[j].cell]-=1
        end
    end
    return n
end
function removeedge_recursive!(a,h,facet)
    edgekey(org(a,h),dest(a,h)) in a.protected_segments && return false
    ring=spin(a,h)
    counters=zeros(Int,length(a.cells))
    for current in ring
        counters[current.cell]=1
    end
    nn=recursive_star!(a,ring,counters,facet)
    if nn>2
        for index in 1:nn
            counters[ring[index].cell]=0
        end
    end
    return nn==2
end

function recoverface!(a,facet)
    facet=Tuple(Int32.(facet))
    for index in 0:2
        for _ in 1:max(100,10length(a.cells))
            first=facet[index+1];last=facet[mod(index+1,3)+1];target=facet[mod(index+2,3)+1]
            start,kind=finddirection!(a,point2tetorg(a,first),last)
            kind===:across_vertex && dest(a,start)==last ||
                throw(Unsupported("the target facet edge has not been recovered"))
            ring=spin(a,start)
            for h in ring
                if apex(a,h)==target
                    push!(a.protected_faces,facekey(facet))
                    result=h
                    for _ in 1:index
                        result=eprev(result)
                    end
                    return result
                end
            end
            edge=NULL
            for h in ring
                if crossing(a,facet,apex(a,h),oppo(a,h))!==:none
                    edge=edestoppo(h)
                    break
                end
            end
            edge.cell!=0 || throw(Unsupported("no crossing edge found in the ordered facet-edge spin"))
            removeedge_recursive!(a,edge,facet) && continue
            break
        end
    end
    throw(Unsupported("the missing facet requires a larger-star/hull/Steiner recovery branch"))
end

function native_stack!(a,objects)
    stack=eltype(objects)[]
    for (index,item) in enumerate(objects)
        position=Int(randomnation!(a,index))+1
        if position<=length(stack)
            push!(stack,stack[position]);stack[position]=item
        else
            push!(stack,item)
        end
    end
    return stack
end
function recoverboundary!(a,facets)
    # unifysegments retains the first physical subface-edge occurrence and
    # its directed vertices; explicit Curve marking does not reverse it.
    segments=Tuple{Int32,Int32}[];seen=Set{Tuple{Int32,Int32}}()
    for face in facets,index in 1:3
        edge=(Int32(face[index]),Int32(face[mod(index,3)+1]))
        key=edgekey(edge...)
        if !(key in seen)
            push!(seen,key);push!(segments,edge)
        end
    end
    segment_stack=native_stack!(a,segments)
    segment_order=Tuple{Int32,Int32}[]
    while !isempty(segment_stack)
        edge=pop!(segment_stack);push!(segment_order,edge)
        h,kind=finddirection!(a,point2tetorg(a,edge[1]),edge[2])
        kind===:across_vertex && dest(a,h)==edge[2] ||
            throw(Unsupported("missing boundary segment requires native edge recovery"))
        push!(a.protected_segments,edgekey(edge...))
    end
    face_stack=native_stack!(a,NTuple{3,Int32}[Tuple(Int32.(face)) for face in facets])
    face_order=NTuple{3,Int32}[]
    while !isempty(face_stack)
        facet=pop!(face_stack);push!(face_order,facet)
        recoverface!(a,facet)
    end
    return (physical_segments=segments,segment_order=segment_order,
            facet_order=face_order,seed=a.seed)
end

function real_cells(a)
    result=NTuple{4,Int32}[]
    slots=Int32[]
    for (index,cell) in enumerate(a.cells)
        cell.alive && all(node->0<node<=a.nreal,cell.vertices) || continue
        push!(result,cell.vertices);push!(slots,Int32(index))
    end
    return result,slots
end
function consistency(a)
    for (index,cell) in enumerate(a.cells)
        cell.alive || continue
        for version in 0:11
            h=Handle(index,version)
            encoded=cell.links[(version&3)+1]
            encoded.cell==0 && continue
            other=fsym(a,h)
            a.cells[other.cell].alive && fsym(a,other)==h || return false
        end
    end
    return true
end

"""
    recover_convex(kernel, facets; seed=1)

Recover a closed convex PLC from explicit primary-node row facets. This bounded
arena admits existing boundary segments, strict interior three-cell flips and
two-to-three reductions through one reflex-link level. It retains physical pool
slot order and returns native negative-orient3/positive-MTet corner tuples.

Hull/dummy flips, missing segments, coplanar intersections, deeper link levels,
rollback, Steiner insertion, nonconvex carving and later Lawson/optimization
stages require separate native paths. They are explicit Unsupported outcomes.
The input kernel is read-only; all working points, pool rows and marks are owned
by the returned arena. Coordinate predicates do not define node identity.
"""
function recover_convex(kernel,facets;seed=1)
    faces=NTuple{3,Int32}[Tuple(Int32.(face)) for face in facets]
    isempty(faces) && throw(Unsupported("a closed PLC facet list is required"))
    edge_counts=Dict{Tuple{Int32,Int32},Int}()
    wanted=Set{NTuple{3,Int32}}()
    for face in faces
        length(unique(face))==3 && all(node->0<node<=kernel.count,face) ||
            throw(Unsupported("PLC facets require three distinct current primary rows"))
        key=facekey(face)
        key in wanted && throw(Unsupported("duplicate PLC facets require an explicit native subface binding"))
        push!(wanted,key)
        for index in 1:3
            edge=edgekey(face[index],face[mod(index,3)+1])
            edge_counts[edge]=get(edge_counts,edge,0)+1
        end
        positive=false;negative=false
        for node in 1:kernel.count
            sign=orient3(kernel.points[face[1]],kernel.points[face[2]],
                kernel.points[face[3]],kernel.points[node])
            positive|=sign>0;negative|=sign<0
            positive && negative && throw(Unsupported("nonconvex/outside PLC points require native carving"))
        end
        (positive || negative) || throw(Unsupported("a flat PLC requires the native degenerate-domain path"))
    end
    all(==(2),values(edge_counts)) || throw(Unsupported("PLC facet edges must form a closed two-face manifold"))
    arena=Arena(kernel;seed)
    protocol=recoverboundary!(arena,faces)
    consistency(arena) || throw(Unsupported("recovered physical arena has inconsistent face bonds"))
    cells,slots=real_cells(arena)
    incidence=Dict{NTuple{3,Int32},Int}()
    for cell in cells
        orient3((arena.points[node] for node in cell)...)<0 ||
            throw(Unsupported("recovered interior cells require native positive orientation"))
        for omitted in 1:4
            face=Tuple(cell[index] for index in 1:4 if index!=omitted)
            key=facekey(face)
            incidence[key]=get(incidence,key,0)+1
        end
    end
    all(count->count<=2,values(incidence)) || throw(Unsupported("recovered cells are not a manifold"))
    observed=Set(key for (key,count) in incidence if count==1)
    observed==wanted || throw(Unsupported("real-cell boundary differs from the explicit PLC facet rows"))
    return (cells=cells,pool_slots=slots,arena=arena,protocol=protocol)
end
include("ModelNativeTetgenQuality.jl")
end
