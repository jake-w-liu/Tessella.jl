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

function _gmsh3_push!(queue,entry)
    heap=queue.heap
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
    heap=queue.heap
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
    delete!(queue.pending,first.key)
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

function _gmsh3_initial_sizes(T,mesh,exterior,minimize,embedded_sizes;native_arithmetic=false)
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
    normalized>=1 || return nothing
    push!(queue.pending,key)
    _gmsh3_push!(queue,_GmshTetEntry(normalized,key,Int32(t),numbers[t]))
    return nothing
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

function _gmsh3_star!(T,point,tags,protected_faces)
    while true
        _gmsh3_shell!(T)
        added=false
        for t in T.cavity,slot in 1:4
            nb=_nbr(T,t,slot)
            (nb!=0 && T.mark[nb]==T.epoch) && continue
            face=_face(T,t,slot)
            orient3(_pt(T,face[1]),_pt(T,face[2]),_pt(T,face[3]),point)<0 && continue
            (nb!=0 && tags[nb]==tags[t] &&
             !(_sort3t(face...) in protected_faces)) || return false
            T.mark[nb]=T.epoch;push!(T.cavity,nb)
            added=true
            break
        end
        added || return true
    end
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

function _gmsh3_compatible(T,point,lc,protected_faces,protected_edges)
    shell_vertices=Set{Int32}()
    shell_edges=Dict{NTuple{2,Int32},Int}()
    for (a,b,c,_) in T.boundary
        for node in (a,b,c)
            push!(shell_vertices,node)
            p=_pt(T,node)
            hypot(p[1]-point[1],p[2]-point[2],p[3]-point[3])<.05lc && return false
        end
        for (u,v) in ((a,b),(b,c),(c,a))
            key=_edgekey(u,v)
            shell_edges[key]=get(shell_edges,key,0)+1
        end
    end
    all(==(2),values(shell_edges)) || return false
    length(shell_vertices)-length(shell_edges)+length(T.boundary)==2 || return false
    for t in T.cavity,slot in 1:4
        # A replacement must keep every preexisting vertex and every typed
        # embedded face. The input's original node labels remain usable.
        _vert(T,t,slot) in shell_vertices || return false
        face=_face(T,t,slot);nb=_nbr(T,t,slot)
        (_sort3t(face...) in protected_faces && nb!=0 &&
         T.mark[nb]==T.epoch) && return false
    end
    for edge in protected_edges
        touched=false
        for t in T.cavity
            vertices=ntuple(k->_vert(T,t,k),4)
            if edge[1] in vertices && edge[2] in vertices
                touched=true;break
            end
        end
        touched && !haskey(shell_edges,_edgekey(edge...)) && return false
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
        max_nodes::Int=PIPELINE_DEFAULT_MAX_NODES,
        max_tets::Int=PIPELINE_DEFAULT_MAX_TETS,max_iterations::Int=1_000_000)
    caller="native 3-D Delaunay insertion"
    max_nodes=_refine_limit3(max_nodes,"max_nodes";minimum=1)
    max_tets=_refine_limit3(max_tets,"max_tets";minimum=1)
    max_iterations=_refine_limit3(max_iterations,"max_iterations";minimum=1)
    nn=size(mesh.coords,2);nt=size(mesh.tets,2)
    (nn<=max_nodes && nt<=max_tets) || throw(ArgumentError(
        "$caller: input exceeds the construction limits"))
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError("$caller: invalid input — "*join(diagnostic.messages,"; ")))
    T,tags,exterior=_gmsh3_triangulation(mesh)
    numbers=_gmsh3_birth_numbers(element_numbers,nt,caller,"element")
    number=_gmsh3_birth_maximum(element_high_water,numbers,caller,"element")
    vertex_numbers=_gmsh3_birth_numbers(node_numbers,nn,caller,"node")
    vertex_number=_gmsh3_birth_maximum(node_high_water,vertex_numbers,caller,"node")
    storage=native_storage_order ? NTuple{4,Int32}[
        ntuple(slot->mesh.tets[slot,column],4) for column in axes(mesh.tets,2)] : nothing
    sizes,background=_gmsh3_initial_sizes(T,mesh,exterior,minimize_boundary,embedded_sizes;native_arithmetic=native_storage_order)
    queue=_GmshTetQueue();skipped=Set{NTuple{4,Int32}}()
    for t in eachindex(T.alive)
        _gmsh3_queue_tet!(queue,T,t,sizes,background,extend_boundary,skipped,numbers,storage)
    end
    live=nt;iteration=0;newtets=Int32[]
    while true
        entry=_gmsh3_pop!(queue)
        entry===nothing && break
        t=entry.tet
        (T.alive[t] && _gmsh3_key(T,t)==entry.key && !(entry.key in skipped)) || continue
        iteration<max_iterations || throw(ArgumentError("$caller: iteration limit reached"))
        iteration+=1
        point=storage===nothing ? _gmsh3_center(T,t) :
            _gmsh3_native_center(T,storage[t])
        all(isfinite,point) || throw(ArgumentError("$caller: circumcenter overflows Float64"))
        # The stack search MVertex has a birth even when no containing cell is found.
        vertex_number=Base.checked_add(vertex_number,UInt64(1))
        _gmsh3_begin_cavity!(T,t,point,tags,protected_faces)
        containing=_gmsh3_containing_weights(T,point)
        if containing===nothing
            push!(skipped,entry.key);continue
        end
        # The candidate MVertex is allocated before star and distance rejection.
        vertex_number=Base.checked_add(vertex_number,UInt64(1))
        if !_gmsh3_star!(T,point,tags,protected_faces)
            push!(skipped,entry.key);continue
        end
        owner,weights=containing
        propagated=0.0
        for k in 1:4
            propagated+=weights[k]*sizes[_vert(T,owner,k)]
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
        if !_gmsh3_compatible(T,point,lc,protected_faces,protected_edges)
            push!(skipped,entry.key);continue
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
        storage===nothing || resize!(storage,length(T.alive))
        for (ordinal,new) in enumerate(newtets)
            tags[new]=region
            numbers[new]=first_number+UInt64(ordinal)
            storage===nothing || (storage[new]=(_vert(T,new,1),_vert(T,new,2),
                _vert(T,new,4),_vert(T,new,3)))
            _gmsh3_queue_tet!(queue,T,new,sizes,background,extend_boundary,skipped,numbers,storage)
        end
        live=remaining+length(newtets)
    end
    coordinates=Matrix{Float64}(undef,3,T.nreal)
    for node in 1:T.nreal
        coordinates[1,node]=T.x[node];coordinates[2,node]=T.y[node];coordinates[3,node]=T.z[node]
    end
    tets=Matrix{Int32}(undef,4,live)
    output_tags=isempty(mesh.tet_tag) ? Int32[] : Vector{Int32}(undef,live)
    output_numbers=return_numbers ? Vector{UInt64}(undef,live) : UInt64[]
    output_rows=findall(T.alive)
    if storage!==nothing
        radii=Float64[_gmsh3_key(T,t) in skipped ? 0. :
            _gmsh3_native_radius(T,storage[t],sizes,background,extend_boundary) for t in output_rows]
        pool_rows=output_rows
        permutation=sortperm(eachindex(pool_rows);lt=(a,b)->radii[a]>radii[b] ? true :
            radii[a]<radii[b] ? false : numbers[pool_rows[a]]<numbers[pool_rows[b]])
        output_rows=output_rows[permutation]
    end
    column=0
    for t in output_rows
        column+=1
        vertices=storage===nothing ? (_vert(T,t,1),_vert(T,t,2),_vert(T,t,4),_vert(T,t,3)) : storage[t]
        for slot in 1:4;tets[slot,column]=vertices[slot];end
        isempty(output_tags) || (output_tags[column]=tags[t])
        return_numbers && (output_numbers[column]=numbers[t])
    end
    output=Mesh(coordinates;segs=copy(mesh.segs),tris=copy(mesh.tris),tets=tets,
        seg_tag=copy(mesh.seg_tag),tri_tag=copy(mesh.tri_tag),tet_tag=output_tags)
    validate(output).ok || throw(ErrorException("$caller: invalid replacement mesh"))
    Set(first(boundary_faces(output.tets)))==Set(first(boundary_faces(mesh.tets))) ||
        throw(ErrorException("$caller: exterior boundary changed"))
    return return_numbers ? _GmshDelaunay3Product(output,output_numbers,number,vertex_numbers,vertex_number) : output
end
