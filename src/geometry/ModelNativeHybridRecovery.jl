# Operation-specific segment recovery for the hybrid Quad fan. Generic

# recover_convex keeps its conservative admission guards. This interior first

# pass follows native recoveredgebyflips/removefacebyflips; larger edge stars,

# full-search, hull and Steiner branches remain explicit Unsupported outcomes.

module ModelNativeHybridRecovery

const A=parentmodule(@__MODULE__).ModelNativeTetgenArena



# Native checkflipeligibility distinguishes the newly created edge from

# the other two edges of each new face. Endpoint sharing is topological.

function intersection(arena,face,segment)

    any(node->node in face,segment) && return (:share_vertex,0)

    p,q,r=face;first,last=segment

    a=A.orient(arena,p,q,r,first);b=A.orient(arena,p,q,r,last)

    (a>0 && b>0 || a<0 && b<0) && return (:none,0)

    (a==0 || b==0) && throw(A.Unsupported(

        "hybrid segment flip eligibility needs native touching/coplanar classification"))

    values=(A.orient(arena,first,last,p,q),A.orient(arena,first,last,q,r),A.orient(arena,first,last,r,p))

    all(>=(0),values) || all(<=(0),values) || return (:none,0)

    position=findfirst(iszero,values)

    return position===nothing ? (:across_face,0) : (:across_edge,position)

end

function eligible23(arena,segment,face,e,d)
    segment===nothing && return true

    for corner in face

        kind,edge=intersection(arena,(e,d,corner),segment)

        (kind===:across_face || kind===:across_edge && edge==1) && return false

    end

    return true

end

function eligible32(arena,segment,face)
    segment===nothing && return true

    kind,_=intersection(arena,face,segment)

    return kind!==:across_face

end

struct QualityConstraint
    input::Float64
    output::Base.RefValue{Float64}
end
function quality_allowed(arena,constraint::QualityConstraint,candidates)
    for cell in candidates
        value=minimum(A.quality_dihedrals(ntuple(i->arena.points[cell[i]],4),A.QualityWorkspace()))
        difference=value-constraint.input
        abs(difference/constraint.input)<1e-12 && (difference=0.0)
        difference>0 || return false
        constraint.output[]=min(constraint.output[],value)
    end
    return true
end
eligible23(arena,constraint::QualityConstraint,face,e,d)=
    quality_allowed(arena,constraint,((e,d,face[2],face[3]),(e,d,face[3],face[1])))
function eligible32(arena,constraint::QualityConstraint,face,e,d,level,edgepivot)
    p,q,r=face
    candidates=level==0 ? ((p,q,r,d),(q,p,r,e)) : edgepivot==1 ? ((q,p,r,e),) : ((p,q,r,d),)
    return quality_allowed(arena,constraint,candidates)
end
eligible32(arena,segment,face,e,d,level,edgepivot)=eligible32(arena,segment,face)
function removeedge!(arena,current,segment)

    A.edgekey(A.org(arena,current),A.dest(arena,current)) in arena.protected_segments && return false

    ring=A.spin(arena,current);counters=zeros(Int,length(arena.cells))

    for cell in ring;counters[cell.cell]=1;end

    result=hybrid_recursive_star!(arena,ring,counters,segment)

    return result==2

end



function hybrid_recursive_star!(a,ring,counters,facet;level=0,edgepivot=0,maxlevel=1)

    n=length(ring)

    n>=3 || throw(A.Unsupported("an edge star has fewer than three cells"))

    pa=A.org(a,ring[1]);pb=A.dest(a,ring[1])

    if n==3

        first,second=ring[1],ring[2]

        eligible32(a,facet,(A.apex(a,ring[1]),A.apex(a,ring[2]),A.apex(a,ring[3])),pa,pb,level,edgepivot) || return n

        A.flip32!(a,ring) || return n

        for h in ring

            counters[h.cell]=0

        end

        ring[1]=A.Handle(first.cell,11);ring[2]=A.Handle(second.cell,11)

        return 2

    end

    reflex=0

    for i in 1:n

        previous=mod(i-2,n)+1;following=mod(i,n)+1

        A.facekey((A.org(a,ring[i]),A.dest(a,ring[i]),A.apex(a,ring[i]))) in a.protected_faces && continue

        (counters[ring[i].cell]>1 || counters[ring[previous].cell]>1) && continue

        pc=A.apex(a,ring[i]);pd=A.apex(a,ring[following]);pe=A.apex(a,ring[previous])

        reducible=false

        if A.orient(a,pb,pc,pd,pe)>0 && A.orient(a,pc,pa,pd,pe)>0

            sign=A.orient(a,pa,pb,pd,pe)

            reducible=sign>0 || (sign==0 && n==4)

        end

        if reducible

            if eligible23(a,facet,(pa,pb,pc),pe,pd)

                result=A.flip23!(a,ring[i],A.fsym(a,ring[i]),counters)

                result===nothing && continue

                replacement=A.edestoppo(result[1])

                counters[replacement.cell]+=1

                ring[previous]=replacement

                deleteat!(ring,i)

                nn=hybrid_recursive_star!(a,ring,counters,facet;level,edgepivot,maxlevel)

                nn>2 && sign==0 && throw(A.Unsupported("flat face23 requires the native rollback branch"))

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

        pc=A.apex(a,ring[i]);pd=A.apex(a,ring[following]);pe=A.apex(a,ring[previous])

        pivot=0;selected=A.NULL

        if A.orient(a,pb,pc,pd,pe)<=0

            selected=A.enext(ring[i]);pivot=1

        elseif A.orient(a,pc,pa,pd,pe)<=0

            selected=A.eprev(ring[i]);pivot=2

        end

        pivot==0 && continue

        A.edgekey(A.org(a,selected),A.dest(a,selected)) in a.protected_segments && continue

        selected=A.esym(selected)

        child=A.spin(a,selected)

        sum(counters[h.cell] for h in child)>2 && continue

        for h in child

            counters[h.cell]+=1

        end

        nn=hybrid_recursive_star!(a,child,counters,facet;level=level+1,edgepivot=pivot,maxlevel)

        if nn==2

            replacement=if pivot==1

                # tmpabtets[0]=[d,a,e,b] => [a,b,e,d].

                A.enext(A.esym(A.enext(child[1])))

            else

                # tmpabtets[1]=[b,d,e,a] => [a,b,e,d].

                A.eprev(A.esym(A.eprev(child[2])))

            end

            (A.org(a,replacement),A.dest(a,replacement))==(pa,pb) ||

                throw(A.Unsupported("reflex-link replacement lost the parent edge phase"))

            counters[replacement.cell]+=1

            ring[previous]=replacement

            deleteat!(ring,i)

            return hybrid_recursive_star!(a,ring,counters,facet;level,edgepivot,maxlevel)

        end

        for j in 1:nn

            counters[child[j].cell]-=1

        end

    end

    return n

end

function removeface!(arena,first,segment)

    face=A.facekey((A.org(arena,first),A.dest(arena,first),A.apex(arena,first)))

    face in arena.protected_faces && return false

    second=A.fsym(arena,first)

    p,q,r,d=A.ordered(arena,first);e=A.oppo(arena,second)

    selected=if A.orient(arena,p,q,d,e)<=0

        first

    elseif A.orient(arena,q,r,d,e)<=0

        A.enext(first)

    elseif A.orient(arena,r,p,d,e)<=0

        A.eprev(first)

    else

        counters=zeros(Int,length(arena.cells))

        A.flip23!(arena,first,second,counters)===nothing && return false

        return true

    end

    return removeedge!(arena,selected,segment)

end



function recoversegment!(arena,start,endpoint)

    for _ in 1:max(100,10length(arena.cells))

        current,kind=A.finddirection!(arena,A.point2tetorg(arena,start),endpoint)

        if kind===:across_vertex

            A.dest(arena,current)==endpoint || throw(A.Unsupported(

                "hybrid segment crosses another existing constructor vertex"))

            push!(arena.protected_segments,A.edgekey(start,endpoint))

            return current

        end

        # Native enextesym goes from [start,b,c,d] to its opposite face.

        crossing=A.esym(A.enext(current))

        if kind===:across_face

            removeface!(arena,crossing,(start,endpoint)) && continue

        elseif kind===:across_edge

            removeedge!(arena,crossing,(start,endpoint)) && continue

        else

            throw(A.Unsupported("hybrid segment has an unsupported native crossing class"))

        end

        return nothing

    end

    throw(A.Unsupported("hybrid segment recovery exceeded the native first-pass work bound"))

end



function recoverboundary!(arena,facets;embedded_segments=NTuple{2,Int32}[],
                          embedded_facets=NTuple{3,Int32}[])

    segments=Tuple{Int32,Int32}[];seen=Set{Tuple{Int32,Int32}}()

    for face in vcat(facets,embedded_facets),index in 1:3

        edge=(Int32(face[index]),Int32(face[mod1(index+1,3)]))

        key=A.edgekey(edge...)

        if !(key in seen)

            push!(seen,key);push!(segments,edge)

        end

    end

    # Interior constraints use their actual lower cell endpoints. They
    # share real segment recovery/protection with the boundary, but never
    # become exterior facets.
    for segment in embedded_segments
        edge=Tuple(Int32.(segment))
        allunique(edge) && all(row->0<row<=length(arena.points),edge) || throw(A.Unsupported(
            "hybrid embedded segment has missing or duplicate current rows"))
        key=A.edgekey(edge...)
        key in seen && throw(A.Unsupported("hybrid embedded segment repeats a boundary or constraint edge"))
        push!(seen,key);push!(segments,edge)
    end

    stack=A.native_stack!(arena,segments);order=Tuple{Int32,Int32}[]

    while !isempty(stack)

        edge=pop!(stack);push!(order,edge)

        result=recoversegment!(arena,edge...)

        result===nothing && (result=recoversegment!(arena,reverse(edge)...))

        result!==nothing || throw(A.Unsupported(

            "hybrid boundary segment requires native full-search or Steiner recovery"))

    end

    face_stack=A.native_stack!(arena,NTuple{3,Int32}[Tuple(Int32.(face)) for face in vcat(facets,embedded_facets)])

    face_order=NTuple{3,Int32}[]

    while !isempty(face_stack)

        facet=pop!(face_stack);push!(face_order,facet)

        A.recoverface!(arena,facet)

    end

    return (physical_segments=segments,segment_order=order,facet_order=face_order,seed=arena.seed)

end

end
