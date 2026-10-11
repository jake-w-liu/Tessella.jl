# Detached native initial-mesh implementation. The caller supplies the actual
# producer insertion order; this core does not infer identity from coordinates.
# PolyMesh's face/edge/vertex allocation order follows meshPolyMesh.h, and the
# walk, edge recovery and region coloring follow meshTriangulation.cpp.

struct _GmshOrientScratch
    B::Vector{Float64}
    u::Vector{Float64}
    C1::Vector{Float64}
    C2::Vector{Float64}
    D::Vector{Float64}
end
_GmshOrientScratch()=_GmshOrientScratch(zeros(4),zeros(4),zeros(8),zeros(12),zeros(16))

function _gmsh_orient2_value(a,b,c,s::_GmshOrientScratch)
    ax,ay=a;bx,by=b;cx,cy=c
    left=(ax-cx)*(by-cy);right=(ay-cy)*(bx-cx)
    det=left-right
    if left>0
        right<=0 && return det
        total=left+right
    elseif left<0
        right>=0 && return det
        total=-left-right
    else
        return det
    end
    bound=Predicates.CCWERRBOUND_A*total
    (det>=bound || -det>=bound) && return det
    acx=ax-cx;bcx=bx-cx;acy=ay-cy;bcy=by-cy
    lh,ll=Predicates._two_product_pair(acx,bcy)
    rh,rl=Predicates._two_product_pair(acy,bcx)
    b3,b2,b1,b0=Predicates._two_two_diff(lh,ll,rh,rl)
    Predicates._store4!(s.B,b3,b2,b1,b0)
    det=Predicates._estimate(4,s.B)
    bound=Predicates.CCWERRBOUND_B*total
    (det>=bound || -det>=bound) && return det
    axt=Predicates._two_diff_tail(ax,cx,acx)
    bxt=Predicates._two_diff_tail(bx,cx,bcx)
    ayt=Predicates._two_diff_tail(ay,cy,acy)
    byt=Predicates._two_diff_tail(by,cy,bcy)
    iszero(axt) && iszero(bxt) && iszero(ayt) && iszero(byt) && return det
    bound=Predicates.CCWERRBOUND_C*total+Predicates.RESULTERRBOUND*abs(det)
    det+=(acx*byt+bcy*axt)-(acy*bxt+bcx*ayt)
    (det>=bound || -det>=bound) && return det
    l1,l0=Predicates._two_product_pair(axt,bcy)
    r1,r0=Predicates._two_product_pair(ayt,bcx)
    u3,u2,u1,u0=Predicates._two_two_diff(l1,l0,r1,r0)
    Predicates._store4!(s.u,u3,u2,u1,u0)
    n1=Predicates._fast_expansion_sum_zeroelim!(4,s.B,4,s.u,s.C1)
    l1,l0=Predicates._two_product_pair(acx,byt)
    r1,r0=Predicates._two_product_pair(acy,bxt)
    u3,u2,u1,u0=Predicates._two_two_diff(l1,l0,r1,r0)
    Predicates._store4!(s.u,u3,u2,u1,u0)
    n2=Predicates._fast_expansion_sum_zeroelim!(n1,s.C1,4,s.u,s.C2)
    l1,l0=Predicates._two_product_pair(axt,byt)
    r1,r0=Predicates._two_product_pair(ayt,bxt)
    u3,u2,u1,u0=Predicates._two_two_diff(l1,l0,r1,r0)
    Predicates._store4!(s.u,u3,u2,u1,u0)
    nd=Predicates._fast_expansion_sum_zeroelim!(n2,s.C2,4,s.u,s.D)
    return s.D[nd]
end

struct _GmshPolyMesh
    x::Vector{Float64}
    y::Vector{Float64}
    vertex_edge::Vector{Int32}
    vertex_label::Vector{Int64}
    origin::Vector{Int32}
    face::Vector{Int32}
    previous::Vector{Int32}
    next::Vector{Int32}
    opposite::Vector{Int32}
    edge_data::Vector{Int32}
    face_edge::Vector{Int32}
    face_data::Vector{Int32}
    orient::_GmshOrientScratch
end
_GmshPolyMesh()=_GmshPolyMesh(Float64[],Float64[],Int32[],Int64[],
    Int32[],Int32[],Int32[],Int32[],Int32[],Int32[],Int32[],Int32[],_GmshOrientScratch())

function _gmsh_poly_vertex!(p,x,y,label=-1)
    length(p.x)<typemax(Int32) || throw(ArgumentError("native initial mesh: vertex capacity exceeds Int32"))
    push!(p.x,x);push!(p.y,y);push!(p.vertex_edge,0);push!(p.vertex_label,label)
    return Int32(length(p.x))
end

function _gmsh_poly_edge!(p,vertex)
    length(p.origin)<typemax(Int32) || throw(ArgumentError("native initial mesh: edge capacity exceeds Int32"))
    push!(p.origin,vertex);push!(p.face,0);push!(p.previous,0)
    push!(p.next,0);push!(p.opposite,0);push!(p.edge_data,-1)
    return Int32(length(p.origin))
end

function _gmsh_poly_face!(p,edge;data=-1)
    length(p.face_edge)<typemax(Int32) || throw(ArgumentError("native initial mesh: face capacity exceeds Int32"))
    push!(p.face_edge,edge);push!(p.face_data,data)
    return Int32(length(p.face_edge))
end

function _gmsh_poly_create_face!(p,f,a,b,c,h0,h1,h2)
    p.origin[h0]=a;p.origin[h1]=b;p.origin[h2]=c
    p.vertex_edge[a]=h0;p.vertex_edge[b]=h1;p.vertex_edge[c]=h2
    p.next[h0]=h1;p.previous[h1]=h0
    p.next[h1]=h2;p.previous[h2]=h1
    p.next[h2]=h0;p.previous[h0]=h2
    p.face[h0]=f;p.face[h1]=f;p.face[h2]=f;p.face_edge[f]=h0
    return nothing
end

function _gmsh_poly_rectangle(xmin,xmax,ymin,ymax)
    p=_GmshPolyMesh()
    mm=_gmsh_poly_vertex!(p,xmin,ymin)
    mM=_gmsh_poly_vertex!(p,xmin,ymax)
    MM=_gmsh_poly_vertex!(p,xmax,ymax)
    Mm=_gmsh_poly_vertex!(p,xmax,ymin)
    h0=_gmsh_poly_edge!(p,mm);h1=_gmsh_poly_edge!(p,MM);h2=_gmsh_poly_edge!(p,Mm)
    f0=_gmsh_poly_face!(p,h0)
    _gmsh_poly_create_face!(p,f0,mm,MM,Mm,h0,h1,h2)
    o0=_gmsh_poly_edge!(p,MM);o1=_gmsh_poly_edge!(p,mm);o2=_gmsh_poly_edge!(p,mM)
    f1=_gmsh_poly_face!(p,o0)
    _gmsh_poly_create_face!(p,f1,MM,mm,mM,o0,o1,o2)
    p.opposite[h0]=o0;p.opposite[o0]=h0
    return p
end

@inline _gmsh_poly_point(p,v)=(p.x[v],p.y[v])
@inline _gmsh_poly_end(p,h)=p.origin[p.next[h]]

function _gmsh_poly_swap!(p,h0)
    o0=p.opposite[h0]
    o0==0 && return -1
    h1=p.next[h0];h2=p.next[h1]
    o1=p.next[o0];o2=p.next[o1]
    v0=p.origin[o1];v1=p.origin[o2];v2=p.origin[o0];v3=p.origin[h2]
    f0=p.face[h0];f1=p.face[o2]
    _gmsh_poly_create_face!(p,f0,v0,v1,v3,o1,o0,h2)
    _gmsh_poly_create_face!(p,f1,v1,v2,v3,o2,h1,h0)
    return 0
end

function _gmsh_poly_delaunay(p,h)
    opposite=p.opposite[h]
    opposite==0 && return false
    a=p.origin[h];b=p.origin[p.next[h]];c=p.origin[p.next[p.next[h]]]
    d=p.origin[p.next[p.next[opposite]]]
    return Predicates.incircle(_gmsh_poly_point(p,a),_gmsh_poly_point(p,b),
        _gmsh_poly_point(p,c),_gmsh_poly_point(p,d))<0
end

function _gmsh_poly_walk(p,f,x,y)
    h=p.face_edge[f];point=(x,y)
    for _ in 1:Base.checked_add(Base.checked_mul(4,length(p.face_edge)),8)
        a=p.origin[h];b=p.origin[p.next[h]];c=p.origin[p.next[p.next[h]]]
        s0=-_gmsh_orient2_value(_gmsh_poly_point(p,a),_gmsh_poly_point(p,b),point,p.orient)
        s1=-_gmsh_orient2_value(_gmsh_poly_point(p,b),_gmsh_poly_point(p,c),point,p.orient)
        s2=-_gmsh_orient2_value(_gmsh_poly_point(p,c),_gmsh_poly_point(p,a),point,p.orient)
        if s0>=0 && s1>=0 && s2>=0
            return p.face[h]
        elseif s0<=0 && s1>=0 && s2>=0
            h=p.opposite[h]
        elseif s1<=0 && s0>=0 && s2>=0
            h=p.opposite[p.next[h]]
        elseif s2<=0 && s0>=0 && s1>=0
            h=p.opposite[p.next[p.next[h]]]
        elseif s0<=0 && s1<=0
            h=s0>s1 ? p.opposite[h] : p.opposite[p.next[h]]
        elseif s0<=0 && s2<=0
            h=s0>s2 ? p.opposite[h] : p.opposite[p.next[p.next[h]]]
        elseif s1<=0 && s2<=0
            h=s1>s2 ? p.opposite[p.next[h]] : p.opposite[p.next[p.next[h]]]
        else
            throw(ArgumentError("native initial mesh: invalid orientation in face walk"))
        end
        h==0 && return Int32(0)
    end
    throw(ArgumentError("native initial mesh: face walk did not terminate"))
end

function _gmsh_poly_split!(p,f,x,y,label)
    v=_gmsh_poly_vertex!(p,x,y,label)
    h0=p.face_edge[f];h1=p.next[h0];h2=p.next[h1]
    v0=p.origin[h0];v1=p.origin[h1];v2=p.origin[h2]
    r0=_gmsh_poly_edge!(p,v);r1=_gmsh_poly_edge!(p,v);r2=_gmsh_poly_edge!(p,v)
    s0=_gmsh_poly_edge!(p,v0);s1=_gmsh_poly_edge!(p,v1);s2=_gmsh_poly_edge!(p,v2)
    for (a,b) in ((r0,s0),(r1,s1),(r2,s2))
        p.opposite[a]=b;p.opposite[b]=a
    end
    f1=_gmsh_poly_face!(p,r1;data=p.face_data[f])
    f2=_gmsh_poly_face!(p,r2;data=p.face_data[f])
    _gmsh_poly_create_face!(p,f,v0,v1,v,h0,s1,r0)
    _gmsh_poly_create_face!(p,f1,v1,v2,v,h1,s2,r1)
    _gmsh_poly_create_face!(p,f2,v2,v0,v,h2,s0,r2)
    stack=Int32[h0,h1,h2];touched=Set{Int32}()
    while !isempty(stack)
        h=pop!(stack);push!(touched,h)
        _gmsh_poly_delaunay(p,h) || continue
        _gmsh_poly_swap!(p,h)
        for k in (h,p.opposite[h])
            k==0 && continue
            b=p.next[k];c=p.next[b]
            !(b in touched) && !(p.opposite[b] in touched) && push!(stack,b)
            !(c in touched) && !(p.opposite[c] in touched) && push!(stack,c)
        end
    end
    return v
end

function _gmsh_poly_get_edge(p,a,b)
    start=p.vertex_edge[a];h=start
    for _ in eachindex(p.origin)
        _gmsh_poly_end(p,h)==b && return h
        h=p.opposite[h]
        h==0 && return Int32(0)
        h=p.next[h]
        h==start && return Int32(0)
    end
    throw(ArgumentError("native initial mesh: invalid vertex edge ring"))
end

function _gmsh_poly_intersect(p,a,b,c,d)
    pa=_gmsh_poly_point(p,a);pb=_gmsh_poly_point(p,b)
    pc=_gmsh_poly_point(p,c);pd=_gmsh_poly_point(p,d)
    s0=_gmsh_orient2_value(pa,pb,pc,p.orient)
    s1=_gmsh_orient2_value(pa,pb,pd,p.orient)
    s0*s1>=0 && return false
    t0=_gmsh_orient2_value(pc,pd,pa,p.orient)
    t1=_gmsh_orient2_value(pc,pd,pb,p.orient)
    return t0*t1<0
end

function _gmsh_poly_recover!(p,a,b)
    start=p.vertex_edge[a];h=start;queue=Int32[]
    for _ in eachindex(p.origin)
        v1=p.origin[p.next[h]];v2=p.origin[p.next[p.next[h]]]
        (v1==b || v2==b) && return 0
        if _gmsh_poly_intersect(p,a,b,v1,v2)
            push!(queue,p.next[h]);break
        end
        h=p.opposite[p.next[p.next[h]]]
        (h==0 || h==start) && break
    end
    isempty(queue) && return -1
    reached=false
    for _ in eachindex(p.origin)
        h=p.opposite[last(queue)]
        h==0 && return -2
        h=p.next[h];v1=p.origin[h];v2=p.origin[p.next[h]]
        if v2==b
            reached=true;break
        elseif _gmsh_poly_intersect(p,a,b,v1,v2)
            push!(queue,h)
        else
            h=p.next[h];v1=p.origin[h];v2=p.origin[p.next[h]]
            if v2==b
                reached=true;break
            elseif _gmsh_poly_intersect(p,a,b,v1,v2)
                push!(queue,h)
            else
                return -3
            end
        end
    end
    reached || throw(ArgumentError("native initial mesh: edge walk did not terminate"))
    crossings=length(queue);position=1
    bound=Base.checked_add(Base.checked_mul(length(p.origin),length(p.origin)),8)
    steps=0
    while position<=length(queue)
        steps+=1
        steps<=bound || throw(ArgumentError("native initial mesh: edge recovery did not terminate"))
        h=queue[position];position+=1;o=p.opposite[h]
        v0=p.origin[h];v1=p.origin[p.next[h]]
        v2=p.origin[p.next[p.next[h]]];v3=p.origin[p.next[p.next[o]]]
        if _gmsh_poly_intersect(p,v0,v1,v2,v3)
            remains=_gmsh_poly_intersect(p,a,b,v2,v3)
            _gmsh_poly_swap!(p,h)
            remains && push!(queue,h)
        else
            push!(queue,h)
        end
    end
    return crossings
end

function _gmsh_poly_color!(p,h,color)
    h!=0 || throw(ArgumentError("native initial mesh: coloring requires a face seed"))
    stack=Int32[p.face[h]];other=Int32(0)
    while !isempty(stack)
        f=pop!(stack);p.face_data[f]=color;edge=p.face_edge[f]
        for _ in 1:3
            opposite=p.opposite[edge]
            if p.edge_data[edge]==-1 && opposite!=0 && p.face_data[p.face[opposite]]==-1
                push!(stack,p.face[opposite])
            elseif p.edge_data[edge]!=-1 && opposite!=0
                other=opposite
            end
            edge=p.next[edge]
        end
    end
    return other
end

function _gmsh_poly_initial(xs,ys,segments,insertion_order;vertex_labels=1:length(xs),surface_tag=1)
    length(xs)==length(ys)==length(vertex_labels) || throw(ArgumentError("native initial mesh: input rows differ"))
    n=length(xs)
    n>=3 && all(isfinite,xs) && all(isfinite,ys) || throw(ArgumentError("native initial mesh: finite planar vertices required"))
    sort!(collect(Int,insertion_order))==collect(1:n) || throw(ArgumentError("native initial mesh: insertion order must be a row permutation"))
    allunique(vertex_labels) || throw(ArgumentError("native initial mesh: source identities must be unique"))
    xmin,xmax=extrema(xs);ymin,ymax=extrema(ys)
    cx=(xmin+xmax)*.5;cy=(ymin+ymax)*.5
    xmin=(xmin-cx)*1.1+cx;xmax=(xmax-cx)*1.1+cx
    ymin=(ymin-cy)*1.1+cy;ymax=(ymax-cy)*1.1+cy
    all(isfinite,(xmin,xmax,ymin,ymax)) && xmin<xmax && ymin<ymax ||
        throw(ArgumentError("native initial mesh: invalid enclosing rectangle"))
    p=_gmsh_poly_rectangle(xmin,xmax,ymin,ymax);rows=zeros(Int32,n);f=Int32(1)
    for row in insertion_order
        f=_gmsh_poly_walk(p,f,xs[row],ys[row])
        f!=0 || throw(ArgumentError("native initial mesh: source vertex lies outside its rectangle"))
        rows[row]=_gmsh_poly_split!(p,f,xs[row],ys[row],vertex_labels[row])
    end
    for (a,b,owner) in segments
        1<=a<=n && 1<=b<=n && a!=b && owner>0 || throw(ArgumentError("native initial mesh: invalid source edge"))
        _gmsh_poly_recover!(p,rows[a],rows[b])>=0 || throw(ArgumentError("native initial mesh: cannot recover source edge"))
        h=_gmsh_poly_get_edge(p,rows[a],rows[b])
        if h!=0
            p.edge_data[h]=owner
            p.opposite[h]!=0 && (p.edge_data[p.opposite[h]]=owner)
        end
    end
    other=_gmsh_poly_color!(p,p.vertex_edge[1],-2)
    _gmsh_poly_color!(p,other,surface_tag)
    for _ in 1:100
        count=0
        for h in eachindex(p.origin)
            o=p.opposite[h]
            o!=0 && p.face_data[p.face[h]]==surface_tag && p.face_data[p.face[o]]==surface_tag || continue
            _gmsh_poly_delaunay(p,h) || continue
            a=p.origin[h];b=p.origin[p.next[h]]
            c=p.origin[p.next[p.next[h]]];d=p.origin[p.next[p.next[o]]]
            if _gmsh_poly_intersect(p,a,b,c,d)
                _gmsh_poly_swap!(p,h);count+=1
            end
        end
        count==0 && break
    end
    return (mesh=p,source_rows=rows)
end
