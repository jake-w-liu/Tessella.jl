# Normalized circumradius insertion for the ordinary native 3-D Delaunay
# mesher. Existing edge-subdivision and point-cloud kernel entry points keep
# their own algorithms. This route preserves the caller's original node IDs.
import ..Predicates: insphere, _orient3_value

struct _GmshTetEntry
    radius::Float64
    key::NTuple{4,Int32}
    tet::Int32
    number::UInt64
end

mutable struct _GmshTetQueue
    heap::Vector{_GmshTetEntry}
    pending::Set{NTuple{4,Int32}}
end

struct _GmshDelaunay3Product
    mesh::Mesh
    element_numbers::Vector{UInt64}
    element_high_water::UInt64
    node_numbers::Vector{UInt64}
    node_high_water::UInt64
end

_GmshTetQueue()=_GmshTetQueue(_GmshTetEntry[],Set{NTuple{4,Int32}}())
@inline _gmsh3_key(T,t)=_sort4t(_vert(T,t,1),_vert(T,t,2),_vert(T,t,3),_vert(T,t,4))
@inline _gmsh3_before(a,b)=a.radius>b.radius ||
    (a.radius==b.radius && a.number<b.number)

_gmsh3_push!(queue::_GmshTetQueue,entry)=_gmsh3_push!(queue.heap,entry)

function _gmsh3_push!(heap::Vector{_GmshTetEntry},entry::_GmshTetEntry)
    push!(heap,entry)
    slot=length(heap)
    @inbounds while slot>1
        parent=slot>>1
        _gmsh3_before(entry,heap[parent]) || break
        heap[slot]=heap[parent]
        slot=parent
    end
    @inbounds heap[slot]=entry
    return nothing
end

function _gmsh3_pop!(queue)
    first=_gmsh3_pop_entry!(queue.heap)
    first===nothing || delete!(queue.pending,first.key)
    return first
end

# The heap's best entry, removed; `nothing` when it is empty.
function _gmsh3_pop_entry!(heap::Vector{_GmshTetEntry})
    isempty(heap) && return nothing
    first=heap[1]
    last=pop!(heap)
    if !isempty(heap)
        slot=1
        @inbounds while 2slot<=length(heap)
            child=2slot
            child<length(heap) && _gmsh3_before(heap[child+1],heap[child]) && (child+=1)
            _gmsh3_before(heap[child],last) || break
            heap[slot]=heap[child]
            slot=child
        end
        @inbounds heap[slot]=last
    end
    return first
end

@inline _gmsh3_sub(a,b)=ntuple(k->a[k]-b[k],3)
@inline _gmsh3_dot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _gmsh3_cross(a,b)=(a[2]*b[3]-a[3]*b[2],
    a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
@inline _gmsh3_det(a,b,c,d)=_gmsh3_dot(_gmsh3_sub(b,a),
    _gmsh3_cross(_gmsh3_sub(c,a),_gmsh3_sub(d,a)))

function _gmsh3_exact_center(a,b,c,d)
    R=Rational{BigInt}
    origin=ntuple(k->R(a[k]),3)
    u=ntuple(k->R(b[k])-origin[k],3)
    v=ntuple(k->R(c[k])-origin[k],3)
    w=ntuple(k->R(d[k])-origin[k],3)
    vw=_gmsh3_cross(v,w);wu=_gmsh3_cross(w,u);uv=_gmsh3_cross(u,v)
    denominator=2*_gmsh3_dot(u,vw)
    iszero(denominator) && throw(ArgumentError(
        "native 3-D Delaunay insertion: degenerate tetrahedron"))
    uu=_gmsh3_dot(u,u);vv=_gmsh3_dot(v,v);ww=_gmsh3_dot(w,w)
    return ntuple(k->Float64(origin[k]+(uu*vw[k]+vv*wu[k]+ww*uv[k])/denominator),3)
end

function _gmsh3_center(T,t)
    a=_pt(T,_vert(T,t,1));b=_pt(T,_vert(T,t,2))
    c=_pt(T,_vert(T,t,3));d=_pt(T,_vert(T,t,4))
    raw_u=_gmsh3_sub(b,a);raw_v=_gmsh3_sub(c,a);raw_w=_gmsh3_sub(d,a)
    magnitude=max(maximum(abs,raw_u),maximum(abs,raw_v),maximum(abs,raw_w))
    if !(isfinite(magnitude) && magnitude>0)
        return _gmsh3_exact_center(a,b,c,d)
    end
    scale=0x1p-200<=magnitude<=0x1p200 ? 1.0 : exp2(exponent(magnitude))
    u=ntuple(k->raw_u[k]/scale,3)
    v=ntuple(k->raw_v[k]/scale,3)
    w=ntuple(k->raw_w[k]/scale,3)
    vw=_gmsh3_cross(v,w);wu=_gmsh3_cross(w,u);uv=_gmsh3_cross(u,v)
    determinant=_gmsh3_dot(u,vw)
    bound=sum(abs(u[k]*vw[k]) for k in 1:3)
    abs(determinant)>64eps(Float64)*bound || return _gmsh3_exact_center(a,b,c,d)
    uu=_gmsh3_dot(u,u);vv=_gmsh3_dot(v,v);ww=_gmsh3_dot(w,w)
    center=ntuple(k->a[k]+scale*(uu*vw[k]+vv*wu[k]+ww*uv[k])/(2determinant),3)
    all(isfinite,center) || return _gmsh3_exact_center(a,b,c,d)
    return center
end

function _gmsh3_native_center(a,b,c,d)
    u=ntuple(axis->b[axis]-a[axis],3)
    v=ntuple(axis->c[axis]-a[axis],3)
    w=ntuple(axis->d[axis]-a[axis],3)
    uu=u[1]*u[1]+u[2]*u[2]+u[3]*u[3]
    vv=v[1]*v[1]+v[2]*v[2]+v[3]*v[3]
    ww=w[1]*w[1]+w[2]*w[2]+w[3]*w[3]
    vw=(v[2]*w[3]-w[2]*v[3],v[3]*w[1]-w[3]*v[1],v[1]*w[2]-w[1]*v[2])
    wu=(w[2]*u[3]-u[2]*w[3],w[3]*u[1]-u[3]*w[1],w[1]*u[2]-u[1]*w[2])
    uv=(u[2]*v[3]-v[2]*u[3],u[3]*v[1]-v[3]*u[1],u[1]*v[2]-v[1]*u[2])
    orientation=_orient3_value(b,c,d,a)
    if isfinite(orientation) && iszero(orientation)
        # An exactly-flat constructor cell has no finite circumcentre; its
        # barycentre is the interior candidate the cavity refinement needs
        # to consume it (upstream's insertion admits the same transient
        # cells rather than rejecting them).
        return ntuple(axis->(a[axis]+b[axis]+c[axis]+d[axis])/4,3)
    end
    (isfinite(orientation) && !iszero(orientation)) || return _gmsh3_exact_center(a,b,c,d)
    denominator=.5/orientation
    result=ntuple(axis->(uu*vw[axis]+vv*wu[axis]+ww*uv[axis])*denominator+a[axis],3)
    return all(isfinite,result) ? result : _gmsh3_exact_center(a,b,c,d)
end

@inline _gmsh3_native_center(T,vertices)=_gmsh3_native_center(
    _pt(T,vertices[1]),_pt(T,vertices[2]),_pt(T,vertices[3]),_pt(T,vertices[4]))

@inline function _gmsh3_length(pa,pb,native_arithmetic)
    dx=pa[1]-pb[1];dy=pa[2]-pb[2];dz=pa[3]-pb[3]
    return native_arithmetic ? sqrt(dx*dx+dy*dy+dz*dz) : hypot(dx,dy,dz)
end

@inline function _gmsh3_native_radius(T,vertices,sizes,background,extend)
    points=ntuple(slot->_pt(T,vertices[slot]),4)
    # A co-spherical coplanar constructor cell has an infinite circumradius —
    # it is the queue's highest-priority refinement target upstream too.
    _orient3_value(points[2],points[3],points[4],points[1])==0.0 && return floatmax(Float64)
    center=_gmsh3_native_center(points...)
    dx=points[1][1]-center[1];dy=points[1][2]-center[2];dz=points[1][3]-center[3]
    radius=sqrt(dx*dx+dy*dy+dz*dz)
    lc=.25*(background[vertices[1]]+background[vertices[2]]+
              background[vertices[3]]+background[vertices[4]])
    if extend
        lc=min(lc,.25*(sizes[vertices[1]]+sizes[vertices[2]]+
                        sizes[vertices[3]]+sizes[vertices[4]]))
    end
    return radius/lc
end

function _gmsh3_birth_numbers(values,count,caller,kind)
    values===nothing && return UInt64.(1:count)
    length(values)==count && all(value->value isa Integer && !(value isa Bool) &&
        0<value<=typemax(UInt64),values) || throw(ArgumentError(
        "$caller: $kind numbers must bind every row with positive UInt64 integers"))
    result=UInt64.(values)
    allunique(result) || throw(ArgumentError("$caller: $kind numbers must be unique"))
    return result
end

function _gmsh3_birth_maximum(value,numbers,caller,kind)
    maximum_number=maximum(numbers;init=UInt64(0))
    value===nothing && return maximum_number
    value isa Integer && !(value isa Bool) && maximum_number<=value<=typemax(UInt64) ||
        throw(ArgumentError("$caller: $kind high-water mark does not cover its born numbers"))
    return UInt64(value)
end

function _gmsh3_triangulation(mesh)
    T=Triangulation3(collect(mesh.coords[1,:]),collect(mesh.coords[2,:]),
                     collect(mesh.coords[3,:]))
    tags=Int32[]
    for column in axes(mesh.tets,2)
        a,b,c,d=mesh.tets[:,column]
        orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d))<0 && ((c,d)=(d,c))
        _newtet!(T,a,b,c,d)
        push!(tags,isempty(mesh.tet_tag) ? Int32(0) : mesh.tet_tag[column])
    end
    faces=Dict{NTuple{3,Int32},Tuple{Int32,Int}}()
    complete=Set{NTuple{3,Int32}}()
    for t in eachindex(T.alive),slot in 1:4
        face=_face(T,t,slot);key=_sort3t(face...)
        key in complete && throw(ArgumentError(
            "native 3-D Delaunay insertion: a face has more than two incident tetrahedra"))
        previous=pop!(faces,key,nothing)
        if previous===nothing
            faces[key]=(Int32(t),slot)
        else
            neighbour,other_slot=previous
            _setn!(T,t,slot,neighbour);_setn!(T,neighbour,other_slot,t)
            push!(complete,key)
        end
    end
    isempty(T.alive) || (T.last=Int32(1))
    return T,tags,faces
end

@inline function _gmsh3_update_size!(sizes,a,b,T,minimize,native_arithmetic=false)
    pa=_pt(T,a);pb=_pt(T,b)
    length=_gmsh3_length(pa,pb,native_arithmetic)
    sizes[a]=minimize ? min(sizes[a],length) : max(sizes[a],length)
    sizes[b]=minimize ? min(sizes[b],length) : max(sizes[b],length)
    return nothing
end

function _gmsh3_initial_sizes(T,mesh,exterior,minimize,embedded_sizes;native_arithmetic=false,
                              sheet_faces=())
    sizes=fill(minimize ? Inf : 0.0,T.nreal)
    boundary=falses(T.nreal)
    for (node,value) in embedded_sizes
        1<=node<=T.nreal || throw(ArgumentError(
            "native 3-D Delaunay insertion: embedded size references an unknown node"))
        (isfinite(value) && value>0) || throw(ArgumentError(
            "native 3-D Delaunay insertion: embedded size must be finite and positive"))
        sizes[node]=value
    end
    for face in keys(exterior)
        for node in face;boundary[node]=true;end
        for (i,j) in ((1,2),(2,3),(3,1))
            _gmsh3_update_size!(sizes,face[i],face[j],T,minimize,native_arithmetic)
        end
    end
    # Embedded surface triangles also contribute sizes and boundary markers.
    for column in axes(mesh.tris,2)
        face=(mesh.tris[1,column],mesh.tris[2,column],mesh.tris[3,column])
        for node in face;boundary[node]=true;end
        for (i,j) in ((1,2),(2,3),(3,1))
            _gmsh3_update_size!(sizes,face[i],face[j],T,minimize,native_arithmetic)
        end
    end
    # setLcs runs over the triangles of every GFace, embedded ones included:
    # the protected sheet faces of an embedded Surface are boundary vertices
    # sized by their own triangle edges, never by the tetrahedra around them.
    # (A face also in `exterior` or `mesh.tris` changes nothing: the update
    # is idempotent.)
    for face in sheet_faces
        all(node->1<=node<=T.nreal,face) || continue
        for node in face;boundary[node]=true;end
        for (i,j) in ((1,2),(2,3),(3,1))
            _gmsh3_update_size!(sizes,face[i],face[j],T,minimize,native_arithmetic)
        end
    end
    for t in eachindex(T.alive),i in 1:4,j in i+1:4
        a=_vert(T,t,i);b=_vert(T,t,j)
        pa=_pt(T,a);pb=_pt(T,b)
        length=_gmsh3_length(pa,pb,native_arithmetic)
        boundary[a] || (sizes[a]=minimize ? min(sizes[a],length) : max(sizes[a],length))
        boundary[b] || (sizes[b]=minimize ? min(sizes[b],length) : max(sizes[b],length))
    end
    all(value->isfinite(value) && value>0,sizes) || throw(ArgumentError(
        "native 3-D Delaunay insertion: an initial node has no usable incident size"))
    return sizes,copy(sizes)
end

function _gmsh3_queue_tet!(queue,T,t,sizes,background,extend,skipped,numbers,storage=nothing)
    T.alive[t] || return nothing
    key=_gmsh3_key(T,t)
    (key in skipped || key in queue.pending) && return nothing
    normalized=_gmsh3_normalized_radius(T,t,sizes,background,extend,storage)
    normalized>=1 || return nothing
    push!(queue.pending,key)
    _gmsh3_push!(queue,_GmshTetEntry(normalized,key,Int32(t),numbers[t]))
    return nothing
end

# The queue of `_refine_gmsh_delaunay3`: a cell is offered once, when it is
# made, and every cell it is offered with contains the vertex just inserted,
# so no offered cell can be pending or refused already. A hash set of the
# pending and refused vertex sets would answer "no" to every question.
function _gmsh3_offer_tet!(heap::Vector{_GmshTetEntry},T,t,sizes,background,extend,numbers,storage)
    @inbounds T.alive[t] || return nothing
    normalized=_gmsh3_normalized_radius(T,t,sizes,background,extend,storage)
    normalized>=1 || return nothing
    _gmsh3_push!(heap,_GmshTetEntry(normalized,_gmsh3_key(T,t),Int32(t),numbers[t]))
    return nothing
end

function _gmsh3_normalized_radius(T,t,sizes,background,extend,storage)
    normalized=if storage===nothing
        center=_gmsh3_center(T,t);point=_pt(T,_vert(T,t,1))
        radius=hypot(point[1]-center[1],point[2]-center[2],point[3]-center[3])
        lc=sum(background[_vert(T,t,k)]/4 for k in 1:4)
        extend && (lc=min(lc,sum(sizes[_vert(T,t,k)]/4 for k in 1:4)))
        radius/lc
    else
        _gmsh3_native_radius(T,storage[t],sizes,background,extend)
    end
    isfinite(normalized) || throw(ArgumentError(
        "native 3-D Delaunay insertion: normalized circumradius is not finite"))
    return normalized
end

function _gmsh3_begin_cavity!(T,seed,point,tags,protected_faces)
    slots=length(T.alive)
    previous=length(T.mark)
    resize!(T.mark,slots)
    previous<slots && fill!(@view(T.mark[previous+1:slots]),Int32(0))
    if T.epoch==typemax(Int32)
        fill!(T.mark,Int32(0));T.epoch=Int32(0)
    end
    T.epoch+=Int32(1)
    empty!(T.cavity);empty!(T.stack)
    push!(T.cavity,Int32(seed));push!(T.stack,Int32(seed))
    T.mark[seed]=T.epoch
    while !isempty(T.stack)
        t=pop!(T.stack)
        for slot in 1:4
            nb=_nbr(T,t,slot)
            (nb!=0 && T.mark[nb]!=T.epoch && tags[nb]==tags[t]) || continue
            _sort3t(_face(T,t,slot)...) in protected_faces && continue
            a=_pt(T,_vert(T,nb,1));b=_pt(T,_vert(T,nb,2))
            c=_pt(T,_vert(T,nb,3));d=_pt(T,_vert(T,nb,4))
            insphere(a,b,c,d,point)>0 || continue
            T.mark[nb]=T.epoch
            push!(T.cavity,nb);push!(T.stack,nb)
        end
    end
    return nothing
end

function _gmsh3_shell!(T)
    empty!(T.boundary)
    for t in T.cavity,slot in 1:4
        nb=_nbr(T,t,slot)
        (nb!=0 && T.mark[nb]==T.epoch) && continue
        face=_face(T,t,slot)
        push!(T.boundary,(face[1],face[2],face[3],nb))
    end
    return nothing
end

# Grow the cavity across every boundary face the point does not see, taking
# the first such face in cavity-and-slot order each time. A boundary face
# stops being one only when its neighbour joins the cavity and never becomes
# one again, so after an addition the next such face is never before the one
# just used: the scan resumes there rather than restarting (same cavity, same
# order, linear instead of quadratic in the cavity size).
function _gmsh3_star!(T,point,tags,protected_faces)
    i=1
    while i<=length(T.cavity)
        t=T.cavity[i]
        for slot in 1:4
            nb=_nbr(T,t,slot)
            (nb!=0 && T.mark[nb]==T.epoch) && continue
            face=_face(T,t,slot)
            orient3(_pt(T,face[1]),_pt(T,face[2]),_pt(T,face[3]),point)<0 && continue
            (nb!=0 && tags[nb]==tags[t] &&
             !(_sort3t(face...) in protected_faces)) || return false
            T.mark[nb]=T.epoch;push!(T.cavity,nb)
        end
        i+=1
    end
    _gmsh3_shell!(T)
    return true
end

# --- Gmsh's own order (meshGRegionDelaunayInsertion.cpp), on cells stored in
# MTetrahedron vertex order (`storage`). Face i of a cell is faces[i] of
# that file; it is opposite storage vertex _GMSH_OPPOSITE[i]. Every face so
# taken is oriented like `_face` (the opposite vertex on its negative side).
const _GMSH_FACES=((1,2,3),(1,3,4),(1,4,2),(2,4,3))
const _GMSH_OPPOSITE=(4,2,3,1)
@inline function _gmsh3_slot(T,t,x)
    @inbounds for k in 1:4
        _vert(T,t,k)==x && return k
    end
    throw(ArgumentError("native 3-D Delaunay insertion: stored cell order lost a vertex"))
end
@inline function _gmsh3_gface(s,i)
    f=_GMSH_FACES[i]
    return (s[f[1]],s[f[2]],s[f[3]])
end
@inline _gmsh3_gnbr(T,storage,t,i)=_nbr(T,t,_gmsh3_slot(T,t,storage[t][_GMSH_OPPOSITE[i]]))

# findCavity(): breadth first from the worst cell; a face whose neighbour is
# missing, across an embedded face, of another region or not in conflict
# joins the shell when it is met.
function _gmsh3_gmsh_cavity!(T,storage,seed,point,tags,protected_faces,shell)
    slots=length(T.alive)
    previous=length(T.mark)
    resize!(T.mark,slots)
    previous<slots && fill!(@view(T.mark[previous+1:slots]),Int32(0))
    if T.epoch==typemax(Int32)
        fill!(T.mark,Int32(0));T.epoch=Int32(0)
    end
    T.epoch+=Int32(1)
    empty!(T.cavity);empty!(shell)
    push!(T.cavity,Int32(seed));T.mark[seed]=T.epoch
    head=1
    while head<=length(T.cavity)
        t=T.cavity[head];head+=1
        for i in 1:4
            nb=_gmsh3_gnbr(T,storage,t,i)
            if nb!=0 && _sort3t(_gmsh3_gface(storage[t],i)...) in protected_faces
                nb=0
            end
            if nb==0
                push!(shell,(Int32(t),Int8(i)))
            elseif T.mark[nb]!=T.epoch
                a=_pt(T,_vert(T,nb,1));b=_pt(T,_vert(T,nb,2))
                c=_pt(T,_vert(T,nb,3));d=_pt(T,_vert(T,nb,4))
                if insphere(a,b,c,d,point)>0 && tags[nb]==tags[t]
                    T.mark[nb]=T.epoch;push!(T.cavity,nb)
                else
                    push!(shell,(Int32(t),Int8(i)))
                end
            end
        end
    end
    return nothing
end

# MTetrahedron::xyz2uvw (sys3x3 on getMat) and isInside with the reference
# element tolerance, over the cavity in order: the first cell containing the
# point and its barycentric weights in stored vertex order.
function _gmsh3_gmsh_containing(T,storage,point;tolerance=1e-6)
    for t in T.cavity
        s=storage[t]
        p0=_pt(T,s[1]);p1=_pt(T,s[2]);p2=_pt(T,s[3]);p3=_pt(T,s[4])
        m00=p1[1]-p0[1];m01=p2[1]-p0[1];m02=p3[1]-p0[1]
        m10=p1[2]-p0[2];m11=p2[2]-p0[2];m12=p3[2]-p0[2]
        m20=p1[3]-p0[3];m21=p2[3]-p0[3];m22=p3[3]-p0[3]
        b0=point[1]-p0[1];b1=point[2]-p0[2];b2=point[3]-p0[3]
        det=m00*(m11*m22-m12*m21)-m01*(m10*m22-m12*m20)+m02*(m10*m21-m11*m20)
        if det==0.0
            u=v=w=0.0
        else
            ud=1.0/det
            u=b0*(m11*m22-m12*m21)-m01*(b1*m22-m12*b2)+m02*(b1*m21-m11*b2)
            v=m00*(b1*m22-m12*b2)-b0*(m10*m22-m12*m20)+m02*(m10*b2-b1*m20)
            w=m00*(m11*b2-b1*m21)-m01*(m10*b2-b1*m20)+b0*(m10*m21-m11*m20)
            u*=ud;v*=ud;w*=ud
        end
        (u< -tolerance || v< -tolerance || w< -tolerance || u>((1.0+tolerance)-v-w)) && continue
        return t,(1.0-u-v-w,u,v,w)
    end
    return nothing
end

@inline function _gmsh3_shell_find(shell,storage,key)
    for (position,(t,i)) in enumerate(shell)
        _sort3t(_gmsh3_gface(storage[t],i)...)==key && return position
    end
    return 0
end

# makeCavityStarShaped() until no face is wrong: collect every shell face the
# point does not see, then extend across each one still on the shell, the
# neighbour's four faces toggling in or out of the shell (extendCavity).
function _gmsh3_gmsh_star!(T,storage,point,tags,protected_faces,shell,wrong)
    while true
        empty!(wrong)
        for (t,i) in shell
            f=_gmsh3_gface(storage[t],i)
            orient3(_pt(T,f[1]),_pt(T,f[2]),_pt(T,f[3]),point)<0 || push!(wrong,(t,i))
        end
        isempty(wrong) && return true
        for (t,i) in wrong
            _gmsh3_shell_find(shell,storage,_sort3t(_gmsh3_gface(storage[t],i)...))==0 && continue
            nb=_gmsh3_gnbr(T,storage,t,i)
            (nb!=0 && tags[nb]==tags[t] &&
             !(_sort3t(_gmsh3_gface(storage[t],i)...) in protected_faces)) || return false
            for j in 1:4
                position=_gmsh3_shell_find(shell,storage,_sort3t(_gmsh3_gface(storage[nb],j)...))
                if position==0
                    push!(shell,(Int32(nb),Int8(j)))
                else
                    deleteat!(shell,position)
                end
            end
            push!(T.cavity,nb);T.mark[nb]=T.epoch
        end
    end
end

function _gmsh3_gmsh_boundary!(T,storage,shell)
    empty!(T.boundary)
    for (t,i) in shell
        f=_gmsh3_gface(storage[t],i)
        push!(T.boundary,(f[1],f[2],f[3],_gmsh3_gnbr(T,storage,t,i)))
    end
    return nothing
end

function _gmsh3_weights(a,b,c,d,point)
    determinant=_gmsh3_det(a,b,c,d)
    weights=(_gmsh3_det(point,b,c,d)/determinant,
        _gmsh3_det(a,point,c,d)/determinant,
        _gmsh3_det(a,b,point,d)/determinant,
        _gmsh3_det(a,b,c,point)/determinant)
    if all(isfinite,weights)
        return weights
    end
    R=Rational{BigInt}
    ar=ntuple(k->R(a[k]),3);br=ntuple(k->R(b[k]),3)
    cr=ntuple(k->R(c[k]),3);dr=ntuple(k->R(d[k]),3)
    pr=ntuple(k->R(point[k]),3)
    denominator=_gmsh3_det(ar,br,cr,dr)
    return (Float64(_gmsh3_det(pr,br,cr,dr)/denominator),
        Float64(_gmsh3_det(ar,pr,cr,dr)/denominator),
        Float64(_gmsh3_det(ar,br,pr,dr)/denominator),
        Float64(_gmsh3_det(ar,br,cr,pr)/denominator))
end

function _gmsh3_containing_weights(T,point)
    for t in T.cavity
        inside=true
        for slot in 1:4
            face=_face(T,t,slot)
            if orient3(_pt(T,face[1]),_pt(T,face[2]),_pt(T,face[3]),point)>0
                inside=false;break
            end
        end
        inside || continue
        a=_pt(T,_vert(T,t,1));b=_pt(T,_vert(T,t,2))
        c=_pt(T,_vert(T,t,3));d=_pt(T,_vert(T,t,4))
        return t,_gmsh3_weights(a,b,c,d,point)
    end
    return nothing
end

# Reusable containers of `_gmsh3_compatible`: one cavity shell is examined
# per insertion, so its vertex and edge tables are emptied rather than rebuilt.
# The shell holds a few dozen triangles, so its vertices and edges live in two
# small open-addressing tables whose slots are valid for one call (`epoch`):
# nothing is hashed through `Set`/`Dict` and nothing is cleared.
mutable struct _GmshShellScratch
    epoch::UInt32
    vertex::Vector{Int32}
    vertex_epoch::Vector{UInt32}
    edge::Vector{UInt64}
    edge_count::Vector{UInt8}
    edge_epoch::Vector{UInt32}
end
_GmshShellScratch()=_GmshShellScratch(UInt32(0),zeros(Int32,256),zeros(UInt32,256),
    zeros(UInt64,256),zeros(UInt8,256),zeros(UInt32,256))

function _gmsh3_shell_reset!(scratch::_GmshShellScratch,triangles::Int)
    # Both tables stay at most a quarter full: 3 corners and 3 edge uses a triangle.
    capacity=length(scratch.vertex)
    if capacity<12triangles
        capacity=nextpow(2,12triangles)
        scratch.vertex=zeros(Int32,capacity);scratch.vertex_epoch=zeros(UInt32,capacity)
        scratch.edge=zeros(UInt64,capacity);scratch.edge_count=zeros(UInt8,capacity)
        scratch.edge_epoch=zeros(UInt32,capacity);scratch.epoch=UInt32(0)
    elseif scratch.epoch==typemax(UInt32)
        fill!(scratch.vertex_epoch,UInt32(0));fill!(scratch.edge_epoch,UInt32(0))
        scratch.epoch=UInt32(0)
    end
    scratch.epoch+=UInt32(1)
    return nothing
end

@inline _gmsh3_edge_code(u::Int32,v::Int32)=u<=v ? (UInt64(reinterpret(UInt32,u))<<32)|UInt64(reinterpret(UInt32,v)) :
    (UInt64(reinterpret(UInt32,v))<<32)|UInt64(reinterpret(UInt32,u))

# Slot of `node` in the vertex table; `found` tells whether it was recorded this call.
@inline function _gmsh3_shell_vertex_slot(scratch::_GmshShellScratch,node::Int32)
    mask=length(scratch.vertex)-1
    slot=Int((reinterpret(UInt32,node)*0x9e3779b1)>>7)&mask
    @inbounds while scratch.vertex_epoch[slot+1]==scratch.epoch
        scratch.vertex[slot+1]==node && return slot+1,true
        slot=(slot+1)&mask
    end
    return slot+1,false
end

@inline function _gmsh3_shell_edge_slot(scratch::_GmshShellScratch,code::UInt64)
    mask=length(scratch.edge)-1
    slot=Int((code*0x9e3779b97f4a7c15)>>40)&mask
    @inbounds while scratch.edge_epoch[slot+1]==scratch.epoch
        scratch.edge[slot+1]==code && return slot+1,true
        slot=(slot+1)&mask
    end
    return slot+1,false
end

# The end nodes of the protected edges, as a mask over node numbers. A cavity
# is asked about every edge of every one of its cells at every insertion,
# and nearly all of those edges join interior nodes, which no protected edge
# touches: one bit test per end spares them the two hashed lookups.
function _gmsh3_protected_nodes(protected_edges,nodes::Int)
    mask=falses(nodes)
    for (u,v) in protected_edges
        0<u<=nodes && (mask[u]=true)
        0<v<=nodes && (mask[v]=true)
    end
    return mask
end

function _gmsh3_compatible(T,point,lc,protected_faces,protected_edges,
                           scratch::_GmshShellScratch=_GmshShellScratch(),
                           protected_nodes::Union{Nothing,BitVector}=nothing;
                           keep_vertices::Bool=true)
    _gmsh3_shell_reset!(scratch,length(T.boundary))
    epoch=scratch.epoch
    nvertices=0;nedges=0
    @inbounds for (a,b,c,_) in T.boundary
        for node in (a,b,c)
            slot,found=_gmsh3_shell_vertex_slot(scratch,node)
            if !found
                scratch.vertex[slot]=node;scratch.vertex_epoch[slot]=epoch
                nvertices+=1
            end
            p=_pt(T,node)
            hypot(p[1]-point[1],p[2]-point[2],p[3]-point[3])<.05lc && return false
        end
        for (u,v) in ((a,b),(b,c),(c,a))
            slot,found=_gmsh3_shell_edge_slot(scratch,_gmsh3_edge_code(u,v))
            if found
                # A third triangle on an edge: the shell is not a surface.
                scratch.edge_count[slot]==0x02 && return false
                scratch.edge_count[slot]=0x02
            else
                scratch.edge[slot]=_gmsh3_edge_code(u,v);scratch.edge_epoch[slot]=epoch
                scratch.edge_count[slot]=0x01
                nedges+=1
            end
        end
    end
    # Every edge on exactly two triangles, and the Euler characteristic of a sphere.
    3length(T.boundary)==2nedges || return false
    nvertices-nedges+length(T.boundary)==2 || return false
    guarded=!isempty(protected_faces)
    for t in T.cavity,slot in 1:4
        # A replacement must keep every preexisting vertex and every typed
        # embedded face. The input's original node labels remain usable.
        # Gmsh does not ask this: a vertex inside the cavity is dropped with
        # its cells (_deleteUnusedVertices removes it afterwards).
        keep_vertices && !last(_gmsh3_shell_vertex_slot(scratch,_vert(T,t,slot))) && return false
        guarded || continue
        # Only a face between two cavity cells disappears; look that up
        # before hashing the face.
        nb=_nbr(T,t,slot)
        (nb!=0 && T.mark[nb]==T.epoch &&
         _sort3t(_face(T,t,slot)...) in protected_faces) && return false
    end
    # A protected edge of a cavity cell must survive on the shell. The cavity
    # holds a few dozen cells while a model protects every curve segment, so
    # the cells' own edges are looked up rather than every protected edge
    # tested against every cell.
    isempty(protected_edges) && return true
    masked=protected_nodes!==nothing
    known=masked ? length(protected_nodes) : 0
    for t in T.cavity
        vertices=(_vert(T,t,1),_vert(T,t,2),_vert(T,t,3),_vert(T,t,4))
        for i in 1:3,j in i+1:4
            u=vertices[i];v=vertices[j]
            if masked
                # A node beyond the mask was born after it was built and is
                # no end of a protected edge.
                (0<u<=known && 0<v<=known && @inbounds(protected_nodes[u]) &&
                 @inbounds(protected_nodes[v])) || continue
            end
            ((u,v) in protected_edges || (v,u) in protected_edges) || continue
            last(_gmsh3_shell_edge_slot(scratch,_gmsh3_edge_code(u,v))) || return false
        end
    end
    return true
end

function _refine_gmsh_delaunay3(mesh::Mesh,background_size;
        extend_boundary::Bool=true,minimize_boundary::Bool=false,
        embedded_sizes=Dict{Int32,Float64}(),
        protected_faces=Set{NTuple{3,Int32}}(),
        protected_edges=Set{NTuple{2,Int32}}(),
        element_numbers=nothing,element_high_water=nothing,
        node_numbers=nothing,node_high_water=nothing,
        native_storage_order::Bool=false,return_numbers::Bool=false,
        audit::Bool=true,validated::Bool=false,
        max_nodes::Int=PIPELINE_DEFAULT_MAX_NODES,
        max_tets::Int=PIPELINE_DEFAULT_MAX_TETS,max_iterations::Int=1_000_000)
    caller="native 3-D Delaunay insertion"
    max_nodes=_refine_limit3(max_nodes,"max_nodes";minimum=1)
    max_tets=_refine_limit3(max_tets,"max_tets";minimum=1)
    max_iterations=_refine_limit3(max_iterations,"max_iterations";minimum=1)
    nn=size(mesh.coords,2);nt=size(mesh.tets,2)
    (nn<=max_nodes && nt<=max_tets) || throw(ArgumentError(
        "$caller: input exceeds the construction limits"))
    # `validated`: the caller has run this same `validate(mesh)` on this mesh
    # and it passed (the volume driver's certified fill); it is not repeated.
    if !validated
        diagnostic=validate(mesh)
        diagnostic.ok || throw(ArgumentError("$caller: invalid input — "*join(diagnostic.messages,"; ")))
    end
    T,tags,exterior=_gmsh3_triangulation(mesh)
    numbers=_gmsh3_birth_numbers(element_numbers,nt,caller,"element")
    number=_gmsh3_birth_maximum(element_high_water,numbers,caller,"element")
    vertex_numbers=_gmsh3_birth_numbers(node_numbers,nn,caller,"node")
    vertex_number=_gmsh3_birth_maximum(node_high_water,vertex_numbers,caller,"node")
    storage=native_storage_order ? NTuple{4,Int32}[
        ntuple(slot->mesh.tets[slot,column],4) for column in axes(mesh.tets,2)] : nothing
    sizes,background=_gmsh3_initial_sizes(T,mesh,exterior,minimize_boundary,embedded_sizes;native_arithmetic=native_storage_order,
        sheet_faces=protected_faces)
    # Every cell is offered to the queue once, when it is made, and the cells
    # made by an insertion all contain the new vertex: an offered cell is never
    # pending or refused already, and an entry is popped at most once. So the
    # queue needs no set of pending cells, and a refused cell is marked by its
    # slot (cleared when the slot is reused) instead of hashed by its vertices.
    heap=_GmshTetEntry[];skipped=falses(length(T.alive))
    shell_scratch=_GmshShellScratch()
    # `protected_edges` is not edited while this loop runs.
    protected_nodes=_gmsh3_protected_nodes(protected_edges,nn)
    for t in eachindex(T.alive)
        _gmsh3_offer_tet!(heap,T,t,sizes,background,extend_boundary,numbers,storage)
    end
    live=nt;iteration=0;newtets=Int32[]
    gshell=Tuple{Int32,Int8}[];gwrong=Tuple{Int32,Int8}[]
    while true
        entry=_gmsh3_pop_entry!(heap)
        entry===nothing && break
        t=entry.tet
        (T.alive[t] && _gmsh3_key(T,t)==entry.key) || continue
        iteration<max_iterations || throw(ArgumentError("$caller: iteration limit reached"))
        iteration+=1
        point=storage===nothing ? _gmsh3_center(T,t) :
            _gmsh3_native_center(T,storage[t])
        all(isfinite,point) || throw(ArgumentError("$caller: circumcenter overflows Float64"))
        # The stack search MVertex has a birth even when no containing cell is found.
        vertex_number=Base.checked_add(vertex_number,UInt64(1))
        if storage===nothing
            _gmsh3_begin_cavity!(T,t,point,tags,protected_faces)
            containing=_gmsh3_containing_weights(T,point)
        else
            _gmsh3_gmsh_cavity!(T,storage,t,point,tags,protected_faces,gshell)
            containing=_gmsh3_gmsh_containing(T,storage,point)
        end
        if containing===nothing
            skipped[t]=true;continue
        end
        # The candidate MVertex is allocated before star and distance rejection.
        vertex_number=Base.checked_add(vertex_number,UInt64(1))
        if storage===nothing
            if !_gmsh3_star!(T,point,tags,protected_faces)
                skipped[t]=true;continue
            end
        else
            if !_gmsh3_gmsh_star!(T,storage,point,tags,protected_faces,gshell,gwrong)
                skipped[t]=true;continue
            end
            _gmsh3_gmsh_boundary!(T,storage,gshell)
        end
        owner,weights=containing
        propagated=0.0
        for k in 1:4
            propagated+=weights[k]*sizes[storage===nothing ? _vert(T,owner,k) : storage[owner][k]]
        end
        target=Float64(background_size(point...))
        (isfinite(propagated) && propagated>0 && isfinite(target) && target>0) ||
            throw(ArgumentError("$caller: insertion size must be finite and positive"))
        lc=extend_boundary ? min(propagated,target) : target
        # MTetrahedron objects for every shell facet are constructed before
        # insertVertexB's distance rejection. Their operation identities are
        # therefore consumed even when the candidate cavity is discarded.
        candidates=UInt64(length(T.boundary))
        candidates<=typemax(UInt64)-number || throw(ArgumentError(
            "$caller: candidate element number allocation overflows"))
        first_number=number
        number+=candidates
        if !_gmsh3_compatible(T,point,lc,protected_faces,protected_edges,shell_scratch,
                              protected_nodes;keep_vertices=storage===nothing)
            skipped[t]=true;continue
        end
        T.nreal<max_nodes || throw(ArgumentError("$caller: node limit reached"))
        remaining=live-length(T.cavity)
        length(T.boundary)<=max_tets-remaining || throw(ArgumentError("$caller: tetrahedron limit reached"))
        region=tags[owner]
        push!(T.x,point[1]);push!(T.y,point[2]);push!(T.z,point[3])
        push!(T.vtet,Int32(0));T.nreal+=1
        push!(sizes,propagated);push!(background,target)
        push!(vertex_numbers,vertex_number)
        for old in T.cavity;_killtet!(T,old);end
        empty!(newtets)
        _retriangulate3!(T,T.boundary,Int32(T.nreal),newtets)
        resize!(tags,length(T.alive))
        resize!(numbers,length(T.alive))
        resize!(skipped,length(T.alive))
        storage===nothing || resize!(storage,length(T.alive))
        for (ordinal,new) in enumerate(newtets)
            tags[new]=region
            numbers[new]=first_number+UInt64(ordinal)
            # insertVertexB: MTetrahedron(face vertices in stored order, v).
            if storage!==nothing
                face=T.boundary[ordinal]
                storage[new]=(face[1],face[2],face[3],Int32(T.nreal))
            end
            skipped[new]=false
            _gmsh3_offer_tet!(heap,T,new,sizes,background,extend_boundary,numbers,storage)
        end
        live=remaining+length(newtets)
    end
    # Rows no cell uses any more (swallowed by a cavity, as Gmsh allows)
    # are dropped like _deleteUnusedVertices; the others keep their order.
    used=falses(T.nreal)
    for t in eachindex(T.alive)
        T.alive[t] || continue
        for k in 1:4;used[_vert(T,t,k)]=true;end
    end
    for column in axes(mesh.tris,2),k in 1:3;used[mesh.tris[k,column]]=true;end
    for column in axes(mesh.segs,2),k in 1:2;used[mesh.segs[k,column]]=true;end
    for (row,_) in embedded_sizes;used[row]=true;end
    storage===nothing && fill!(@view(used[1:nn]),true)
    rowmap=zeros(Int32,T.nreal);kept=0
    for node in 1:T.nreal
        used[node] && (kept+=1;rowmap[node]=Int32(kept))
    end
    coordinates=Matrix{Float64}(undef,3,kept)
    for node in 1:T.nreal
        used[node] || continue
        r=rowmap[node]
        coordinates[1,r]=T.x[node];coordinates[2,r]=T.y[node];coordinates[3,r]=T.z[node]
    end
    kept<T.nreal && (vertex_numbers=vertex_numbers[used])
    tets=Matrix{Int32}(undef,4,live)
    output_tags=isempty(mesh.tet_tag) ? Int32[] : Vector{Int32}(undef,live)
    output_numbers=return_numbers ? Vector{UInt64}(undef,live) : UInt64[]
    output_rows=findall(T.alive)
    if storage!==nothing
        # Largest radius first, then the element number: numbers are unique, so
        # sorting (-radius, number, row) triples gives the same order as the
        # comparison sort of the rows, without a closure per comparison.
        order=Vector{Tuple{Float64,UInt64,Int32}}(undef,length(output_rows))
        for (index,t) in enumerate(output_rows)
            radius=skipped[t] ? 0. :
                _gmsh3_native_radius(T,storage[t],sizes,background,extend_boundary)
            order[index]=(-radius,numbers[t],Int32(t))
        end
        sort!(order)
        for index in eachindex(order)
            output_rows[index]=Int(order[index][3])
        end
    end
    column=0
    for t in output_rows
        column+=1
        vertices=storage===nothing ? (_vert(T,t,1),_vert(T,t,2),_vert(T,t,4),_vert(T,t,3)) : storage[t]
        for slot in 1:4;tets[slot,column]=rowmap[vertices[slot]];end
        isempty(output_tags) || (output_tags[column]=tags[t])
        return_numbers && (output_numbers[column]=numbers[t])
    end
    output=Mesh(coordinates;segs=Int32.(rowmap[mesh.segs]),tris=Int32.(rowmap[mesh.tris]),tets=tets,
        seg_tag=copy(mesh.seg_tag),tri_tag=copy(mesh.tri_tag),tet_tag=output_tags)
    audit && (validate(output).ok || throw(ErrorException("$caller: invalid replacement mesh")))
    _single_tet_faces(output.tets)==_single_tet_faces(mesh.tets) ||
        throw(ErrorException("$caller: exterior boundary changed"))
    return return_numbers ? _GmshDelaunay3Product(output,output_numbers,number,vertex_numbers,vertex_number) : output
end
