# Experimental operation-owned interior Lawson recovery. Boundary carving,

# native hull pool slots and recursive unflips are separate explicit stages.

module ModelNativeHybridLawson

const A=parentmodule(@__MODULE__).ModelNativeTetgenArena

using ...Predicates: insphere,orient3

function lawson!(a,outside;queue_cells=nothing)

    stack=A.Handle[];marked=Set{Tuple{Int32,Int}}();unflippable=NamedTuple[]

    counters=zeros(Int,length(a.cells));flips23=0;flips32=0

    pivot(h)=(h.cell,Int(h.version)&3)

    function pushface(h)

        a.cells[h.cell].alive && !(h.cell in outside) || return

        key=pivot(h)

        if !(key in marked);push!(marked,key);push!(stack,h);end

    end

    function clearcells(slots)

        for slot in slots,index in 0:3;delete!(marked,(Int32(slot),index));end

    end

    function after23(result)

        for h in result;pushface(A.esym(A.eprev(h)));end

        for h in result;pushface(A.esym(A.enext(h)));end

    end

    function after32(head,tail)

        pushface(A.esym(A.enext(head)));pushface(A.esym(A.eprev(tail)))

        pushface(A.esym(A.eprev(head)));pushface(A.esym(A.enext(tail)))

        pushface(A.esym(head));pushface(A.esym(tail))

    end

    selected=queue_cells===nothing ? eachindex(a.cells) : [h.cell for h in queue_cells]
    for slot in selected,version in 0:3

        a.cells[slot].alive && !(Int32(slot) in outside) || continue

        h=A.Handle(slot,version)

        encoded=a.cells[slot].links[version+1]

        encoded.cell==0 && continue

        !(pivot(A.fsym(a,h)) in marked) && pushface(h)

    end

    for round in 1:10001

        roundflips=0

        while !isempty(stack)

            h=pop!(stack)

            a.cells[h.cell].alive && pivot(h) in marked || continue

            delete!(marked,pivot(h))

            h.cell in outside && continue

            first=h;second=A.fsym(a,first)

            second.cell in outside && continue # explicit current boundary

            A.facekey((A.org(a,first),A.dest(a,first),A.apex(a,first))) in a.protected_faces && continue

            second_vertices=a.cells[second.cell].vertices
            points=ntuple(i->a.points[second_vertices[i]],4)

            sign=insphere(points...,a.points[A.oppo(a,first)])

            sign<0 || continue

            d=A.oppo(a,first);e=A.oppo(a,second)

            delta=ntuple(i->a.points[d][i]-a.points[e][i],3)

            length3=sqrt(delta[1]^2+delta[2]^2+delta[3]^2)^3

            orientation=0.0

            for _ in 1:3

                orientation=A.orient(a,A.org(a,first),A.dest(a,first),d,e)

                if orientation>0

                    firstboundary=A.facekey((A.org(a,A.esym(first)),A.dest(a,A.esym(first)),A.apex(a,A.esym(first)))) in a.protected_faces

                    secondboundary=A.facekey((A.org(a,A.esym(second)),A.dest(a,A.esym(second)),A.apex(a,A.esym(second)))) in a.protected_faces

                    if firstboundary || secondboundary

                        p=a.points[A.org(a,first)];q=a.points[A.dest(a,first)];r=a.points[d];s=a.points[e]

                        ax=p[1]-s[1];bx=q[1]-s[1];cx=r[1]-s[1]

                        ay=p[2]-s[2];by=q[2]-s[2];cy=r[2]-s[2]

                        az=p[3]-s[3];bz=q[3]-s[3];cz=r[3]-s[3]

                        volume=ax*(by*cz-bz*cy)+bx*(cy*az-cz*ay)+cx*(ay*bz-az*by)

                        abs(volume)/length3<1e-12 && (orientation=0.0)

                    end

                end

                orientation<=0 && break

                first=A.enext(first);second=A.eprev(second)

            end

            if orientation>0

                clearcells((first.cell,second.cell))

                result=A.flip23!(a,first,second,counters)

                result===nothing && throw(A.Unsupported("owned Lawson face23 unexpectedly rejected"))

                clearcells(h.cell for h in result);after23(result)

                flips23+=1;roundflips+=1;continue

            end

            A.edgekey(A.org(a,first),A.dest(a,first)) in a.protected_segments && continue

            edge=A.esym(first);ring=A.spin(a,edge)

            if length(ring)==3 && all(h->!(h.cell in outside),ring)

                clearcells(h.cell for h in ring)

                if A.flip32!(a,ring)

                    after32(A.Handle(ring[1].cell,11),A.Handle(ring[2].cell,11))

                    flips32+=1;roundflips+=1;continue

                end

            elseif length(ring)==4 && orientation==0

                throw(A.Unsupported("owned Lawson requires native flat 4-to-4/hull branch"))

            end

            push!(unflippable,(handle=edge,vertices=A.ordered(a,edge)))

        end

        isempty(unflippable) && return (flips23=flips23,flips32=flips32,unflippable=unflippable)

        roundflips==0 && return (flips23=flips23,flips32=flips32,unflippable=unflippable)

        for record in unflippable

            h=record.handle

            a.cells[h.cell].alive && A.ordered(a,h)==record.vertices && pushface(h)

        end

        empty!(unflippable)

    end

    throw(A.Unsupported("owned Lawson exceeded native 10000 round limit"))

end

end
