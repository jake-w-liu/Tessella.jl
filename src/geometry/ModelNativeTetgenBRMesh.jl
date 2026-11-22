# tetgenmesh::reconstructmesh (meshGRegionBoundaryRecovery.cpp) up to and
# including recoverdelaunay(): the pool of tetrahedra in the Delaunay
# kernel's allocation order, hull tetrahedra, subfaces and segments from the
# surface triangles, unifysegments, recoverboundary, carveholes and the
# Lawson pass.

function _emptytg(points,ninput,seed,epsilon)
    n=length(points)
    m=TG(copy(points),Int32.(1:n),zeros(Int32,n),zeros(Int,n),0,n,ninput,
        Int[],Int32[],Int[],Int[],Int32[],0,0,0,ShellPool(),ShellPool(),UInt32(seed),
        false,false,false,1,-1,-1,epsilon,0.0,0.0,0.0,
        Tri[],BadFace[],Tri[],Tri[],Tri[],Sh[],Sh[],Sh[],Sh[],Tri[],Sh[],Sh[],Int32[],Int32[],
        Vector{Tri}[],0,Vector{Tri}(undef,3),Vector{Tri}[])
    sizehint!(m.nb,32n);sizehint!(m.tv,32n);sizehint!(m.tseg,48n);sizehint!(m.tsub,32n)
    sizehint!(m.tmark,8n)
    for p in 1:n;setpointtype!(m,p,VOLVERTEX);end
    return m
end

# The tetrahedra of the Delaunay kernel, in its allocation order, become the
# first pool slots; faces are bonded, hull tetrahedra appended in traversal
# order, and point2tet points at the last real tetrahedron of each vertex.
function _reconstruct!(m,kernel)
    # Live kernel cells become slots 1, 2, ... in order. A face is bonded
    # when its second tetrahedron is made; the kernel's own adjacency names
    # that neighbour, so no face table is built.
    slotof=zeros(Int32,length(kernel.cells))
    # Kernel face k is opposite vertex OPPOSITE[k]; tetgen face f is opposite
    # vertex f+1.
    facefor=(2,3,4,1)          # kernel face opposite kernel vertex q
    for (ci,source) in enumerate(kernel.cells)
        source.vertices[1]==0 && continue
        p=source.vertices
        ori=orient3(m.xyz[p[1]],m.xyz[p[2]],m.xyz[p[3]],m.xyz[p[4]])
        ori==0 && throw(Unsupported("a degenerate tetrahedron in the initial Delaunay mesh"))
        swapped=ori>0
        swapped && (p=(p[2],p[1],p[3],p[4]))
        t=maketetrahedron!(m)
        slotof[ci]=Int32(t.tet)
        setvertices!(m,t,p[1],p[2],p[3],p[4])
        for f in 0:3
            q=f+1
            swapped && q<=2 && (q=3-q)
            other=source.neighbors[facefor[q]]
            other==0 && continue
            o=Int(slotof[other])
            o==0 && continue
            h=Tri(t.tet,f)
            a=org(m,h);b=dest(m,h);c=apex(m,h)
            bonded=false
            for g0 in 0:3
                x=tvert(m,o,4+g0)
                (x==a || x==b || x==c) && continue
                # The neighbour's face opposite x; its matching version
                # reverses the directed edge [a,b].
                for v in g0:4:11
                    g=Tri(o,v)
                    if org(m,g)==b && dest(m,g)==a && apex(m,g)==c
                        bond!(m,h,g);bonded=true;break
                    end
                end
                break
            end
            bonded || throw(Unsupported("the initial Delaunay adjacency is inconsistent"))
        end
    end
    real=m.tmax
    slot=1
    while slot<=m.tmax
        if tvert(m,slot,4)!=0 && tvert(m,slot,7)!=DUMMY
            for v in 0:3
                tetloop=Tri(slot,v)
                if m.nb[_nbi(tetloop)]==0
                    hulltet=maketetrahedron!(m)
                    p0=org(m,tetloop);p1=dest(m,tetloop);p2=apex(m,tetloop)
                    setvertices!(m,hulltet,p1,p0,p2,DUMMY)
                    bond!(m,tetloop,hulltet)
                    for _ in 1:3
                        face2=fsym(m,hulltet)
                        while true
                            face2.tet==0 && break
                            face2=esym(face2)
                            apex(m,face2)==DUMMY && break
                            x=m.nb[_nbi(face2)]
                            if x==0
                                face2=TNULL;break
                            end
                            face2=fsym(m,face2)
                        end
                        face2.tet!=0 && bond!(m,esym(hulltet),face2)
                        hulltet=enext(hulltet)
                    end
                end
                setpoint2tet!(m,tvert(m,slot,4+v),Tri(slot,0))
            end
        end
        slot+=1
    end
    return real
end

# Subfaces and segments of the surface triangles (in Gmsh's GFace and triangle
# order), then unifysegments().
function _surface!(m,facets,ang_tol)
    n=length(facets)
    sizehint!(m.subfaces.w,11n);sizehint!(m.subfaces.mark,n);sizehint!(m.subfaces.flag,n)
    sizehint!(m.subsegs.w,33n);sizehint!(m.subsegs.mark,3n);sizehint!(m.subsegs.flag,3n)
    for face in facets
        for j in 1:3
            pointtype(m,face[j])==VOLVERTEX && setpointtype!(m,face[j],FACETVERTEX)
        end
        sh=makeshellface!(m.subfaces)
        setshvertices!(m.subfaces,sh,face[1],face[2],face[3])
        for _ in 1:3
            seg=makeshellface!(m.subsegs)
            setshvertices!(m.subsegs,seg,sorg(m.subfaces,sh),sdest(m.subfaces,sh),0)
            m.subsegs.mark[seg.sh]=-1
            ssbond!(m,sh,seg)
            sh=senext(sh)
        end
    end
    _unifysegments!(m,ang_tol)
end

function _unifysegments!(m,ang_tol)
    F=m.subfaces
    n=length(m.xyz)
    counts=zeros(Int,n+1)
    for slot in 1:F.maxitems
        s=Sh(slot,0);shdead(F,s) && continue
        counts[sw(F,s,3)]+=1;counts[sw(F,s,4)]+=1;sw(F,s,5)!=0 && (counts[sw(F,s,5)]+=1)
    end
    first=zeros(Int,n+1);acc=0
    for i in 1:n
        first[i]=acc;acc+=counts[i]
    end
    first[n+1]=acc
    list=Vector{Sh}(undef,acc);fill=copy(first)
    for slot in 1:F.maxitems
        s=Sh(slot,0);shdead(F,s) && continue
        j=sw(F,s,3);fill[j]+=1;list[fill[j]]=Sh(slot,0)
        j=sw(F,s,4);fill[j]+=1;list[fill[j]]=Sh(slot,2)
        j=sw(F,s,5);fill[j]+=1;list[fill[j]]=Sh(slot,4)
    end
    link=Sh[]
    S=m.subsegs
    for slot in 1:S.maxitems
        seg=Sh(slot,0);shdead(S,seg) && continue
        torg=sorg(S,seg);tdest=sdest(S,seg)
        empty!(link)
        for k in first[torg]+1:first[torg+1]
            sface=list[k]
            shdead(F,sface) && continue
            sdest(F,sface)!=tdest && (sface=sesym(senext2(sface)))
            sdest(F,sface)!=tdest && continue
            if length(link)>=2
                at=length(link)
                for q in 1:length(link)-1
                    f1=link[q];f2=link[q+1]
                    ori1=orient3d(m,torg,tdest,sapex(F,f1),sapex(F,f2))
                    ori2=orient3d(m,torg,tdest,sapex(F,f1),sapex(F,sface))
                    stop=false
                    if ori1>0
                        if ori2>0
                            ori3=orient3d(m,torg,tdest,sapex(F,f2),sapex(F,sface))
                            ori3==0 && throw(Unsupported("overlapping facets"))
                            stop=ori3>0
                        elseif ori2<0
                            stop=true
                        else
                            ori3=orient3d(m,torg,tdest,sapex(F,f2),sapex(F,sface))
                            ori3<=0 && throw(Unsupported("overlapping facets"))
                            stop=true
                        end
                    elseif ori1<0
                        if ori2<0
                            ori3=orient3d(m,torg,tdest,sapex(F,f2),sapex(F,sface))
                            ori3==0 && throw(Unsupported("overlapping facets"))
                            stop=ori3>0
                        elseif ori2==0
                            ori3=orient3d(m,torg,tdest,sapex(F,f2),sapex(F,sface))
                            ori3>0 && throw(Unsupported("overlapping facets"))
                        end
                    else
                        if ori2<0
                            stop=true
                        elseif ori2==0
                            throw(Unsupported("overlapping facets"))
                        end
                    end
                    if stop
                        at=q;break
                    end
                end
                insert!(link,at+1,sface)
            elseif length(link)==1
                ori1=orient3d(m,torg,tdest,sapex(F,link[1]),sapex(F,sface))
                if ori1==0
                    n1=facenormal(m.xyz[torg],m.xyz[tdest],m.xyz[sapex(F,link[1])])
                    n2=facenormal(m.xyz[torg],m.xyz[tdest],m.xyz[sapex(F,sface)])
                    dot3(n1,n2)>0 && throw(Unsupported("overlapping facets"))
                end
                push!(link,sface)
            else
                push!(link,sface)
            end
        end
        for f in link
            testseg=sspivot(F,f)
            if testseg.sh!=0 && testseg.sh!=seg.sh && !shdead(S,testseg)
                shellfacedealloc!(S,testseg.sh)
            end
            ssbond!(m,f,seg)
        end
        if length(link)>1
            for k in 1:length(link)
                f1=link[k];f2=k<length(link) ? link[k+1] : link[1]
                n1=facenormal(m.xyz[torg],m.xyz[tdest],m.xyz[sapex(F,f1)])
                n2=facenormal(m.xyz[torg],m.xyz[tdest],m.xyz[sapex(F,f2)])
                cosang=dot3(n1,n2)/(sqrt(dot3(n1,n1))*sqrt(dot3(n2,n2)))
                cosang=clamp(cosang,-1.0,1.0)
                acos(cosang)<ang_tol && throw(Unsupported("overlapping facets"))
                sbond1!(F,f1,f2)
            end
        end
    end
    return nothing
end

function _recoverboundary!(m)
    m.checksubsegflag=true
    S=m.subsegs
    stack=m.subsegstack
    i=0
    for slot in 1:S.maxitems
        seg=Sh(slot,0);shdead(S,seg) && continue
        s=randomnation!(m,i+1)
        push!(stack,s+1<=length(stack) ? stack[s+1] : seg)
        stack[s+1]=seg
        i+=1
    end
    ms=length(stack);nit=0
    m.autofliplinklevel=1
    misseglist=Sh[]
    while true
        recoversegments!(m,misseglist,false)
        isempty(misseglist) && break
        m.fliplinklevel>=0 && break
        if length(misseglist)>=ms
            nit+=1
            nit>=3 && (m.fliplinklevel=100000)
        else
            ms=length(misseglist)
            nit>0 && (nit-=1)
        end
        append!(stack,misseglist);empty!(misseglist)
        m.autofliplinklevel+=1
    end
    while !isempty(misseglist)
        ms=length(misseglist)
        append!(stack,misseglist);empty!(misseglist)
        recoversegments!(m,misseglist,true)
        length(misseglist)<ms || break
    end
    isempty(misseglist) ||
        throw(Unsupported("segment recovery requires Steiner points"))
    m.checksubfaceflag=true
    F=m.subfaces
    fstack=m.subfacstack
    i=0
    for slot in 1:F.maxitems
        sh=Sh(slot,0);shdead(F,sh) && continue
        s=randomnation!(m,i+1)
        push!(fstack,s+1<=length(fstack) ? fstack[s+1] : sh)
        fstack[s+1]=sh
        i+=1
    end
    ms=length(fstack);nit=0
    m.fliplinklevel=-1;m.autofliplinklevel=1
    misshlist=Sh[]
    while true
        recoversubfaces!(m,misshlist)
        isempty(misshlist) && break
        m.fliplinklevel>=0 && break
        if length(misshlist)>=ms
            nit+=1
            nit>=3 && (m.fliplinklevel=100000)
        else
            ms=length(misshlist)
            nit>0 && (nit-=1)
        end
        append!(fstack,misshlist);empty!(misshlist)
        m.autofliplinklevel+=1
    end
    isempty(misshlist) || throw(Unsupported("facet recovery requires Steiner points"))
    return nothing
end

function _carveholes!(m)
    tetarray=Tri[];hullarray=Tri[]
    for slot in 1:m.tmax
        tvert(m,slot,4)==0 && continue
        tetloop=Tri(slot,11)
        if ishulltet(m,tetloop) && !issubface(m,tetloop)
            infect!(m,tetloop);push!(hullarray,tetloop)
            neightet=rawnb(m,slot,3)
            if !infected(m,neightet)
                infect!(m,neightet);push!(tetarray,neightet)
            end
        end
    end
    i=1
    while i<=length(tetarray)
        parytet=tetarray[i]
        j=parytet.ver&3
        for k in 1:3
            neightet=rawnb(m,parytet.tet,(j+k)%4)
            if !infected(m,neightet)
                if !issubface(m,neightet)
                    infect!(m,neightet);push!(tetarray,neightet)
                elseif ishulltet(m,neightet)
                    infect!(m,neightet);push!(hullarray,neightet)
                    checksh=tspivot(m,neightet)
                    sinfect!(m.subfaces,checksh);push!(m.subfacstack,checksh)
                end
            elseif issubface(m,neightet)
                checksh=tspivot(m,neightet)
                if !sinfected(m.subfaces,checksh)
                    sinfect!(m.subfaces,checksh);push!(m.subfacstack,checksh)
                end
            end
        end
        i+=1
    end
    for p in 1:length(m.xyz)
        pointtype(m,p)==DEADVERTEX && continue
        if pointtype(m,p)!=UNUSEDVERTEX && pointtype(m,p)!=DUPLICATEDVERTEX
            neightet=tdec(m.p2tet[p])
            neightet.tet!=0 && infected(m,neightet) && push!(m.cavetetvertlist,Int32(p))
            push!(m.subvertstack,Int32(p))
        end
    end
    if !isempty(tetarray)
        S=m.subsegs;F=m.subfaces
        for slot in 1:S.maxitems
            seg=Sh(slot,0);shdead(S,seg) && continue
            neightet=sstpivot1(m,seg)
            neightet.tet!=0 && infected(m,neightet) && push!(m.subsegstack,seg)
        end
        newhull=Tri[]
        for parytet in tetarray,j in 0:3
            tetloop=rawnb(m,parytet.tet,j)
            if !infected(m,tetloop)
                checksh=tspivot(m,tetloop)
                checksh.sh==0 && throw(Unsupported("carveholes: an unprotected new hull face"))
                hulltet=maketetrahedron!(m)
                pa=org(m,tetloop);pb=dest(m,tetloop);pc=apex(m,tetloop)
                setvertices!(m,hulltet,pb,pa,pc,DUMMY)
                bond!(m,tetloop,hulltet)
                tsbond!(m,hulltet,sesym(checksh))
                for _ in 1:3
                    if issubseg(m,tetloop)
                        seg=tsspivot1(m,tetloop)
                        tssbond1!(m,hulltet,seg);sstbond1!(m,seg,hulltet)
                    end
                    tetloop=enext(tetloop);hulltet=eprev(hulltet)
                end
                setpoint2tet!(m,pa,Tri(tetloop.tet,0))
                setpoint2tet!(m,pb,Tri(tetloop.tet,0))
                setpoint2tet!(m,pc,Tri(tetloop.tet,0))
                push!(newhull,Tri(parytet.tet,j))
            end
        end
        for parytet in newhull
            neightet=fsym(m,parytet)
            hulltet=fsym(m,neightet)
            for _ in 1:3
                casface=esym(hulltet)
                if m.nb[_nbi(casface)]==0
                    neightet=parytet
                    while true
                        neightet=fnext(m,neightet)
                        infected(m,neightet) || break
                    end
                    ishulltet(m,neightet) || (neightet=esym(fsym(m,neightet)))
                    bond!(m,casface,neightet)
                end
                hulltet=enext(hulltet);parytet=enext(parytet)
            end
        end
        if !isempty(m.subfacstack)
            for parysh in m.subfacstack
                cur=parysh
                for _ in 1:3
                    casingout=spivot(F,cur)
                    checkseg=sspivot(F,cur)
                    if casingout.sh!=0
                        casingin=casingout
                        while true
                            checksh=spivot(F,casingin)
                            checksh.sh==cur.sh && break
                            casingin=checksh
                        end
                        if casingin.sh!=casingout.sh
                            sbond1!(F,casingin,casingout)
                        else
                            sdissolve!(F,casingout)
                        end
                        checkseg.sh!=0 && ssbond!(m,casingout,checkseg)
                    elseif checkseg.sh!=0
                        shellfacedealloc!(S,checkseg.sh)
                    end
                    cur=senext(cur)
                end
                shellfacedealloc!(F,cur.sh)
            end
            empty!(m.subfacstack)
        end
        for seg in m.subsegstack
            if !shdead(S,seg)
                neightet=sstpivot1(m,seg)
                infected(m,neightet) && shellfacedealloc!(S,seg.sh)
            end
        end
        empty!(m.subsegstack)
        for p in m.cavetetvertlist
            neightet=tdec(m.p2tet[p])
            infected(m,neightet) && setpointtype!(m,p,UNUSEDVERTEX)
        end
        empty!(m.cavetetvertlist)
        for t in tetarray;tetrahedrondealloc!(m,t.tet);end
        for t in hullarray;tetrahedrondealloc!(m,t.tet);end
    end
    m.nonconvex=true
    for slot in 1:m.tmax
        tvert(m,slot,4)==0 && continue
        tvert(m,slot,7)==DUMMY && flippush!(m,fsym(m,Tri(slot,11)))
    end
    fc=Flip();fc.enqflag=2
    lawsonflip3d!(m,fc)
    empty!(m.unflipqueue)
    return nothing
end

# suppresssteinerpoints(): every vertex is queued by carveholes() under -Y/2;
# only Steiner points (none without boundary Steiner insertion) are acted on.
function _suppresssteinerpoints!(m)
    for p in m.subvertstack
        t=pointtype(m,p)
        t in (FREESEGVERTEX,FREEFACETVERTEX,FREEVOLVERTEX) &&
            throw(Unsupported("suppresssteinerpoints with boundary-recovery Steiner points"))
    end
    empty!(m.subvertstack)
end

function _recoverdelaunay!(m)
    for slot in 1:m.tmax
        (tvert(m,slot,4)==0 || tvert(m,slot,7)==DUMMY) && continue
        for v in 0:3
            neightet=rawnb(m,slot,v)
            facemarked(m,neightet) || flippush!(m,Tri(slot,v))
        end
    end
    fc=Flip()
    fc.bak_tetprism_vol=0.0
    fc.remove_ndelaunay_edge=true
    fc.enqflag=2
    lawsonflip3d!(m,fc)
    isempty(m.unflipqueue) && return nothing
    fc.unflip=true;fc.collectnewtets=true;fc.enqflag=0
    m.autofliplinklevel=1;m.fliplinklevel=-1
    bak=m.flipstarsize;m.flipstarsize=10
    flipqueue=m.unflipqueue;m.unflipqueue=BadFace[]
    nextflipqueue=BadFace[]
    while !isempty(flipqueue)
        i=1
        while i<=length(flipqueue)
            bface=flipqueue[i]
            ok,tt=getedge!(m,bface.forg,bface.fdest,bface.tt)
            bface=BadFace(tt,bface.key,bface.cent,bface.forg,bface.fdest,bface.fapex,bface.foppo)
            flipqueue[i]=bface
            if ok
                nn,tt=removeedgebyflips!(m,tt,fc)
                bface=BadFace(tt,bface.key,bface.cent,bface.forg,bface.fdest,bface.fapex,bface.foppo)
                flipqueue[i]=bface
                if nn==2
                    throw(Unsupported("recoverdelaunay: an improving edge flip (tetprismvol is zero in Gmsh)"))
                else
                    push!(nextflipqueue,bface)
                    fc.tetprism_vol_sum=0.0
                end
            end
            i+=1
        end
        empty!(flipqueue)
        flipqueue,nextflipqueue=nextflipqueue,flipqueue
        if !isempty(flipqueue)
            m.autofliplinklevel>=1 && break
            m.autofliplinklevel+=1
        end
    end
    m.flipstarsize=bak
    empty!(m.unflipqueue)
    return nothing
end
