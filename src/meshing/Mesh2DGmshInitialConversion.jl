# Convert the complete non-box PolyMesh topology to the native insertion
# storage. Region membership is a separate boundary-parity operation: an
# embedded closed curve does not remove the faces enclosed by that curve.
function _gmsh_poly_triangulation(product,boundaries,internals;cleanup_surface=nothing)
    p=product.mesh;rows=product.source_rows;n=length(rows)
    source=zeros(Int32,length(p.x))
    source[rows]=Int32.(1:n)
    T=Mesh2D.Triangulation(p.x[rows],p.y[rows])
    face_order=collect(eachindex(p.face_edge))
    if cleanup_surface!==nothing
        # BDS_Mesh::cleanup uses libstdc++'s bidirectional std::partition.
        # Deleted box/exterior faces are swapped with retained faces from the
        # tail. Keeping this allocation traversal matters to initial-only
        # MTriangle births; a stable filter changes public cell identities.
        kept=Bool[p.face_data[f]==cleanup_surface for f in face_order]
        first=1;last=length(face_order)
        while first<=last
            while first<=last && kept[face_order[first]]
                first+=1
            end
            while first<=last && !kept[face_order[last]]
                last-=1
            end
            first>=last && break
            face_order[first],face_order[last]=face_order[last],face_order[first]
            first+=1;last-=1
        end
    end
    for f in face_order
        h=p.face_edge[f]
        a=source[p.origin[h]]
        b=source[p.origin[p.next[h]]]
        c=source[p.origin[p.next[p.next[h]]]]
        (a==0 || b==0 || c==0) && continue
        # PolyMesh stores real faces clockwise. Insertion uses CCW storage.
        Mesh2D._newtri!(T,a,c,b)
    end
    isempty(T.alive) && throw(ArgumentError(
        "native initial mesh: no real face remains after box removal"))
    unmatched=Dict{NTuple{2,Int32},Tuple{Int32,Int32}}()
    present=Set{NTuple{2,Int32}}()
    real_count=length(T.alive)
    for t in 1:real_count,k in 1:3
        a,b=Mesh2D._edge(T,t,k);key=Mesh2D._skey(a,b)
        push!(present,key)
        if haskey(unmatched,key)
            previous,slot=pop!(unmatched,key)
            Mesh2D._setn!(T,t,k,previous)
            Mesh2D._setn!(T,previous,slot,t)
        else
            unmatched[key]=(Int32(t),Int32(k))
        end
    end
    # Follow actual real-face storage order for the hull ghost allocation.
    for t in 1:real_count,k in 1:3
        Mesh2D._nbr(T,t,k)==0 || continue
        a,b=Mesh2D._edge(T,t,k)
        ghost=Mesh2D._newtri!(T,b,a,Mesh2D.GHOST)
        Mesh2D._setn!(T,t,k,ghost)
        Mesh2D._setn!(T,ghost,3,t)
    end
    empty!(unmatched)
    for t in real_count+1:length(T.alive),k in 1:2
        a,b=Mesh2D._edge(T,t,k);key=Mesh2D._skey(a,b)
        if haskey(unmatched,key)
            previous,slot=pop!(unmatched,key)
            Mesh2D._setn!(T,t,k,previous)
            Mesh2D._setn!(T,previous,slot,t)
        else
            unmatched[key]=(Int32(t),Int32(k))
        end
    end
    isempty(unmatched) || throw(ArgumentError(
        "native initial mesh: real faces do not form a closed hull"))
    for (records,destination) in ((boundaries,T.seg),(internals,T.internal))
        for (a,b) in records
            1<=a<=n && 1<=b<=n && a!=b || throw(ArgumentError(
                "native initial mesh: invalid classified constraint"))
            key=Mesh2D._skey(Int32(a),Int32(b))
            key in present || throw(ArgumentError(
                "native initial mesh: a classified constraint was not recovered"))
            push!(destination,key)
        end
    end
    T.last=Int32(1)
    Mesh2D._assert_consistent(T,"native initial mesh conversion")
    return T
end
