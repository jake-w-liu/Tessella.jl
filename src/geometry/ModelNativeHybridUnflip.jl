# Native interior Delaunay recovery and exact recorded inverse handle phases.
module ModelNativeHybridUnflip
const A=parentmodule(@__MODULE__).ModelNativeTetgenArena
const R=parentmodule(@__MODULE__).ModelNativeHybridRecovery
function metric(a,cell)
    # meshGRegionBoundaryRecovery.cpp defines its local orient4d as zero
    # before including tetgenBR.cxx. Every embedded tetprismvol is zero.
    return 0.0
end
function flip23!(a,first,second,counters,delta)
    before=metric(a,a.cells[first.cell].vertices)+metric(a,a.cells[second.cell].vertices)
    result=A.flip23!(a,first,second,counters)
    result===nothing && return nothing
    after=sum(metric(a,a.cells[h.cell].vertices) for h in result)
    delta[]+=after-before
    return result
end
function flip32!(a,ring,delta)
    before=sum(metric(a,a.cells[h.cell].vertices) for h in ring)
    A.flip32!(a,ring) || return false
    after=metric(a,a.cells[ring[1].cell].vertices)+metric(a,a.cells[ring[2].cell].vertices)
    delta[]+=after-before
    return true
end
function expand!(ring,index,first,second)
    n=length(ring)+1
    insert!(ring,index,second)
    ring[mod(index-2,n)+1]=first
end
function flip3!(a,ring,counters,delta,cutoff,level,pivot,collected)
    if cutoff isa R.QualityConstraint
        face=ntuple(i->A.apex(a,ring[i]),3)
        R.eligible32(a,cutoff,face,A.org(a,ring[1]),A.dest(a,ring[1]),level,pivot) || return false
    end
    flip32!(a,ring,delta) || return false
    first=A.Handle(ring[1].cell,11);second=A.Handle(ring[2].cell,11)
    for h in ring;counters[h.cell]=0;end
    if cutoff isa Float64 && level==0 && (delta[]>=0 || abs(delta[])<cutoff)
        result=flip23!(a,first,second,counters,delta)
        result===nothing && throw(A.Unsupported("hybrid objective unflip could not restore its directed face"))
        empty!(ring);append!(ring,result)
        for h in ring;counters[h.cell]+=1;end
        return false
    end
    empty!(ring);append!(ring,(first,second))
    level==0 ? append!(collected,(first,second)) : push!(collected,pivot==1 ? second : first)
    return true
end
function inverse!(a,ring,history,counters,delta,collected,pivot)
    length(ring)==2 || throw(A.Unsupported("hybrid child inverse lacks its two directed output faces"))
    restored=flip23!(a,ring[1],ring[2],counters,delta)
    restored===nothing && throw(A.Unsupported("hybrid child inverse lost its exact directed face"))
    empty!(ring);append!(ring,restored)
    resize!(collected,length(collected)-(pivot==0 ? 2 : 1))
    for record in history
        index=record.index
        replacement=ring[mod(index-2,length(ring))+1]
        if record.kind===:face
            edge=A.enext(A.esym(A.eprev(replacement)))
            cavity=A.spin(a,edge)
            length(cavity)==3 && flip32!(a,cavity,delta) || throw(A.Unsupported(
                "hybrid recorded face inverse lost its exact three-cell cavity"))
            first=A.Handle(cavity[1].cell,11);second=A.Handle(cavity[2].cell,11)
            counters[first.cell]=0;counters[second.cell]=0
            expand!(ring,index,A.esym(second),first)
            resize!(collected,length(collected)-2)
        elseif record.kind===:child
            child=A.Handle[A.NULL,A.NULL]
            if record.pivot==1
                child[1]=A.eprev(A.esym(A.eprev(replacement)));child[2]=A.fsym(a,child[1])
            else
                child[2]=A.enext(A.esym(A.enext(replacement)));child[1]=A.fsym(a,child[2])
            end
            inverse!(a,child,record.history,counters,delta,collected,record.pivot)
            first=record.pivot==1 ? A.esym(A.enext(child[2])) : A.esym(A.eprev(child[2]))
            second=record.pivot==1 ? A.eprev(A.esym(child[1])) : A.enext(A.esym(child[1]))
            expand!(ring,index,first,second)
        else
            throw(A.Unsupported("hybrid child inverse has an unknown recorded reduction"))
        end
    end
    return nothing
end
function recursive!(a,ring,counters,delta,cutoff,collected;level=0,pivot=0,history=Any[])
    n=length(ring)
    n==3 && return flip3!(a,ring,counters,delta,cutoff,level,pivot,collected)
    n>=4 || throw(A.Unsupported("hybrid Delaunay unflip has an invalid interior star"))
    pa=A.org(a,ring[1]);pb=A.dest(a,ring[1]);reflex=0
    for index in 1:n
        previous=mod(index-2,n)+1;following=mod(index,n)+1
        A.facekey((A.org(a,ring[index]),A.dest(a,ring[index]),A.apex(a,ring[index]))) in a.protected_faces && continue
        (counters[ring[index].cell]>1 || counters[ring[previous].cell]>1) && continue
        pc=A.apex(a,ring[index]);pd=A.apex(a,ring[following]);pe=A.apex(a,ring[previous])
        all(>(0),(pa,pb,pc,pd,pe)) || throw(A.Unsupported("hybrid Delaunay unflip requires the native hull branch"))
        sign=A.orient(a,pa,pb,pd,pe)
        if A.orient(a,pb,pc,pd,pe)>0 && A.orient(a,pc,pa,pd,pe)>0 && sign>=0
            cutoff isa R.QualityConstraint && !R.eligible23(a,cutoff,(pa,pb,pc),pe,pd) && continue
            result=flip23!(a,ring[index],A.fsym(a,ring[index]),counters,delta)
            result===nothing && continue
            append!(collected,(result[2],result[3]))
            replacement=A.edestoppo(result[1]);counters[replacement.cell]+=1
            ring[previous]=replacement;deleteat!(ring,index)
            if recursive!(a,ring,counters,delta,cutoff,collected;level,pivot,history)
                push!(history,(kind=:face,index=index))
                return true
            end
            replacement=ring[mod(index-2,length(ring))+1]
            edge=A.edestoppo(replacement);cavity=A.spin(a,edge)
            length(cavity)==3 && flip32!(a,cavity,delta) || throw(A.Unsupported("hybrid face23 inverse lost its exact three-cell cavity"))
            first=A.Handle(cavity[1].cell,11);second=A.Handle(cavity[2].cell,11)
            counters[first.cell]+=1;counters[second.cell]+=1
            expand!(ring,index,A.esym(second),first)
            resize!(collected,length(collected)-2)
        else
            reflex+=1
        end
    end
    reflex>0 && level<1 || return false
    for index in 1:n
        previous=mod(index-2,n)+1;following=mod(index,n)+1
        (counters[ring[index].cell]>1 || counters[ring[previous].cell]>1) && continue
        pc=A.apex(a,ring[index]);pd=A.apex(a,ring[following]);pe=A.apex(a,ring[previous])
        pivot=0;selected=A.NULL
        if A.orient(a,pb,pc,pd,pe)<=0
            selected=A.enext(ring[index]);pivot=1
        elseif A.orient(a,pc,pa,pd,pe)<=0
            selected=A.eprev(ring[index]);pivot=2
        end
        pivot==0 && continue
        A.edgekey(A.org(a,selected),A.dest(a,selected)) in a.protected_segments && continue
        selected=A.esym(selected);child=A.spin(a,selected)
        length(child)>=3 || throw(A.Unsupported("hybrid Delaunay inverse has an invalid child star"))
        sum(counters[h.cell] for h in child)>2 && continue
        for h in child;counters[h.cell]+=1;end
        child_history=Any[]
        if recursive!(a,child,counters,delta,cutoff,collected;level=level+1,pivot,history=child_history)
            replacement=pivot==1 ? A.enext(A.esym(A.enext(child[1]))) : A.eprev(A.esym(A.eprev(child[2])))
            (A.org(a,replacement),A.dest(a,replacement))==(pa,pb) || throw(A.Unsupported("hybrid child inverse lost its parent edge"))
            counters[replacement.cell]+=1;ring[previous]=replacement;deleteat!(ring,index)
            if recursive!(a,ring,counters,delta,cutoff,collected;level,pivot,history)
                push!(history,(kind=:child,index=index,pivot=pivot,history=child_history))
                return true
            end
            replacement=ring[mod(index-2,length(ring))+1]
            if pivot==1
                child[1]=A.eprev(A.esym(A.eprev(replacement)));child[2]=A.fsym(a,child[1])
            else
                child[2]=A.enext(A.esym(A.enext(replacement)));child[1]=A.fsym(a,child[2])
            end
            inverse!(a,child,child_history,counters,delta,collected,pivot)
            first=pivot==1 ? A.esym(A.enext(child[2])) : A.esym(A.eprev(child[2]))
            second=pivot==1 ? A.eprev(A.esym(child[1])) : A.enext(A.esym(child[1]))
            counters[first.cell]+=1;counters[second.cell]+=1
            expand!(ring,index,first,second)
        else
            for h in child;counters[h.cell]-=1;end
        end
    end
    return false
end
function removeedge!(a,current,cutoff)
    A.edgekey(A.org(a,current),A.dest(a,current)) in a.protected_segments && return (removed=false,delta=0.0)
    ring=A.spin(a,current);counters=zeros(Int,length(a.cells));delta=Ref(0.0)
    for h in ring;counters[h.cell]=1;end
    collected=A.Handle[]
    removed=recursive!(a,ring,counters,delta,cutoff,collected)
    return (removed=removed,delta=delta[],collected=collected)
end
end
