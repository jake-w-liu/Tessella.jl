# Operation-owned native hull slots and carving traversal. Dummy points are
# topology only; interior predicates never receive the sentinel row zero.
module ModelNativeHybridCarve
const A=parentmodule(@__MODULE__).ModelNativeTetgenArena
rawneighbor(a,h,index=Int(h.version)&3)=a.cells[h.cell].links[index+1]
function allocate!(a,vertices)
    slot=if isempty(a.dead)
        push!(a.cells,A.Cell(vertices,(A.NULL,A.NULL,A.NULL,A.NULL),true))
        slot=Int32(length(a.cells))
        push!(a.pool_events,(kind=:allocate,slot=slot,previous_free_top=Int32(0),free_count=0))
        slot
    else
        slot=pop!(a.dead)
        a.cells[slot]=A.Cell(vertices,(A.NULL,A.NULL,A.NULL,A.NULL),true)
        push!(a.pool_events,(kind=:reuse,slot=slot,previous_free_top=slot,free_count=length(a.dead)))
        slot
    end
    return slot
end
function deallocate!(a,slot)
    a.cells[slot].alive || throw(A.Unsupported("hybrid carving deallocates a dead pool slot"))
    a.cells[slot].vertices=(0,0,0,0);a.cells[slot].links=(A.NULL,A.NULL,A.NULL,A.NULL);a.cells[slot].alive=false
    previous=isempty(a.dead) ? Int32(0) : last(a.dead)
    push!(a.dead,slot)
    push!(a.pool_events,(kind=:deallocate,slot=slot,previous_free_top=previous,free_count=length(a.dead)))
end
function connect_hulls!(a,hulls)
    pending=Dict{NTuple{3,Int32},A.Handle}()
    for slot in hulls,version in 0:2
        h=A.Handle(slot,version);face=(A.org(a,h),A.dest(a,h),A.apex(a,h))
        0 in face || throw(A.Unsupported("hybrid hull face lost its dummy point"))
        key=A.facekey(face);other=pop!(pending,key,A.NULL)
        if other.cell==0
            pending[key]=h
        else
            match=A.findhandle(a,slot,A.dest(a,other),A.org(a,other),A.apex(a,other))
            A.bond!(a,other,match)
        end
    end
    isempty(pending) || throw(A.Unsupported("hybrid hull faces do not form a closed manifold"))
end
function initial_hulls!(a)
    isempty(a.dead) || throw(A.Unsupported("hybrid initial hull reconstruction expects untouched native cells"))
    original=length(a.cells);hulls=Int32[]
    for slot in 1:original,version in 0:3
        h=A.Handle(slot,version)
        rawneighbor(a,h).cell==0 || continue
        p,q,r=A.org(a,h),A.dest(a,h),A.apex(a,h)
        new=allocate!(a,(q,p,r,Int32(0)));push!(hulls,new)
        A.bond!(a,h,A.Handle(new,11))
    end
    connect_hulls!(a,hulls)
    return hulls
end
function carve!(a,original_hulls;exterior_facets=nothing)
    frontier=exterior_facets===nothing ? a.protected_faces : exterior_facets
    all(face->face in a.protected_faces,frontier) || throw(A.Unsupported(
        "hybrid carving frontier contains an unrecovered exterior facet"))
    infected=Set{Int32}();queue=A.Handle[];oldhulls=Int32[]
    for slot in original_hulls
        h=A.Handle(slot,11)
        A.facekey((A.org(a,h),A.dest(a,h),A.apex(a,h))) in frontier && continue
        push!(infected,slot);push!(oldhulls,slot)
        other=rawneighbor(a,h)
        other.cell>0 || throw(A.Unsupported("hybrid initial hull has no real neighbor"))
        if !(other.cell in infected)
            push!(infected,other.cell);push!(queue,other)
        end
    end
    position=1
    while position<=length(queue)
        h=queue[position];position+=1;phase=Int(h.version)&3
        for increment in 1:3
            other=rawneighbor(a,h,mod(phase+increment,4))
            other.cell>0 || throw(A.Unsupported("hybrid exterior flood lost its native hull neighbor"))
            other.cell in infected && continue
            key=A.facekey((A.org(a,other),A.dest(a,other),A.apex(a,other)))
            if !(key in frontier)
                0 in a.cells[other.cell].vertices && throw(A.Unsupported("hybrid unprotected hull was not seeded"))
                push!(infected,other.cell);push!(queue,other)
            elseif 0 in a.cells[other.cell].vertices
                throw(A.Unsupported("hybrid carving needs exterior protected-face deletion"))
            end
        end
    end
    hulls=Int32[]
    for exterior in queue,version in 0:3
        interior=rawneighbor(a,exterior,version)
        interior.cell in infected && continue
        key=A.facekey((A.org(a,interior),A.dest(a,interior),A.apex(a,interior)))
        key in frontier || throw(A.Unsupported("hybrid carved frontier lacks its literal subface"))
        p,q,r=A.org(a,interior),A.dest(a,interior),A.apex(a,interior)
        new=allocate!(a,(q,p,r,Int32(0)));push!(hulls,new)
        A.bond!(a,interior,A.Handle(new,11))
        a.hints[p]=interior.cell;a.hints[q]=interior.cell;a.hints[r]=interior.cell
    end
    connect_hulls!(a,hulls)
    for h in queue;deallocate!(a,h.cell);end
    for slot in oldhulls;deallocate!(a,slot);end
    return (hulls=hulls,exterior_order=Int32[h.cell for h in queue],old_hulls=oldhulls)
end
end
