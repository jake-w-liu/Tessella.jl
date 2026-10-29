# A selected Point is an actual lower constructor, not an unreferenced row.
struct _ModelNativeHybridEmbeddedPointInput
    model::GeoModel
    region::Int
    input::_ModelMeshAuthority
    seeds::_ModelNativeHybridSizeSeeds
    record::NamedTuple
    inventory::Vector{NTuple{2,Int}}
end

function _model_native_hybrid_embedded_point_input(m,region,input,caller;point=nothing,
                                                    surface_witnesses=())
    embeds=copy(get(m.embeds,(3,region),NTuple{2,Int}[]))
    if point===nothing
        # The original single-entity route remains separate.
        haskey(m.volumes,region) && length(embeds)==1 && only(embeds)[1]==0 || return nothing
        point=only(embeds)[2]
    else
        haskey(m.volumes,region) && (0,point) in embeds && all(owner->owner[1] in (0,1,2),embeds) ||
            throw(ArgumentError("$caller: selected embedded Point lacks its actual Region declaration"))
    end
    seeds=_model_native_hybrid_global_size_seeds(m,input,caller;surface_witnesses)
    records=[record for record in seeds.records if record.owner==(0,point)]
    length(records)==1 || throw(ArgumentError(
        "$caller: embedded Point lacks its unique actual Point-cell witness"))
    return _ModelNativeHybridEmbeddedPointInput(m,region,input,seeds,only(records),embeds)
end

function _model_native_hybrid_point_source_certificate(witness,input,uids,coordinates,caller)
    witness===nothing && return nothing
    witness isa _ModelNativeHybridEmbeddedPointInput && witness.input===input &&
        get(witness.model.embeds,(3,witness.region),NTuple{2,Int}[])==witness.inventory &&
        witness.record.owner in witness.inventory && all(owner->owner[1] in (0,1,2),witness.inventory) ||
        throw(ArgumentError("$caller: embedded Point has a foreign or changed constructor"))
    record=witness.record
    record in _model_native_hybrid_seed_records(witness.seeds) &&
        record.node_owner==(0,Int32(record.owner[2])) || throw(ArgumentError(
        "$caller: embedded Point loses its exact current Point-cell/UID/owner binding"))
    row=findfirst(==(record.uid),uids)
    row!==nothing && isequal(Tuple(coordinates[:,row]),record.coordinate) ||
        throw(ArgumentError("$caller: embedded Point UID is missing or changed in the native Source"))
    return nothing
end

# The literal fan triangles of one shell panel (a Quad splits on its kept
# diagonal, exactly like the constructor stores it).
_model_native_hybrid_panel_triangles(panel)=
    length(panel.nodes)==3 ? (panel.nodes,) :
        ((panel.nodes[1],panel.nodes[2],panel.nodes[3]),
         (panel.nodes[1],panel.nodes[3],panel.nodes[4]))

# Coplanar point-on-triangle test in the panel's own plane (dominant normal
# axis dropped, exact orient2 edge signs).
function _model_native_hybrid_coplanar_on_triangle(p,a,b,c)
    u=(b[1]-a[1],b[2]-a[2],b[3]-a[3]);v=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
    normal=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    drop=argmax(abs.(collect(normal)));axes=filter(!=(drop),(1,2,3))
    project(x)=(x[axes[1]],x[axes[2]])
    pa,pb,pc,pp=map(project,(a,b,c,p))
    s1=orient2(pa,pb,pp);s2=orient2(pb,pc,pp);s3=orient2(pc,pa,pp)
    return (all(side->side>=0,(s1,s2,s3)) || all(side->side<=0,(s1,s2,s3)))
end

# A coplanar segment properly intersects the triangle when an endpoint lies
# inside it, a triangle edge is straddled, or the segment lies on an edge
# line (collinear overlap is contact, not a miss).
function _model_native_hybrid_coplanar_segment_triangle(p,q,a,b,c)
    u=(b[1]-a[1],b[2]-a[2],b[3]-a[3]);v=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
    normal=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    drop=argmax(abs.(collect(normal)));axes=filter(!=(drop),(1,2,3))
    project(x)=(x[axes[1]],x[axes[2]])
    pa,pb,pc,pp,pq=map(project,(a,b,c,p,q))
    inside(x)=begin
        s1=orient2(pa,pb,x);s2=orient2(pb,pc,x);s3=orient2(pc,pa,x)
        all(side->side>=0,(s1,s2,s3)) || all(side->side<=0,(s1,s2,s3))
    end
    (inside(pp) || inside(pq)) && return true
    for edge in ((pa,pb),(pb,pc),(pc,pa))
        ea,eb=edge
        s1=orient2(pp,pq,ea);s2=orient2(pp,pq,eb)
        t1=orient2(ea,eb,pp);t2=orient2(ea,eb,pq)
        s1==0 && s2==0 && return true
        s1!=0 && s2!=0 && t1!=0 && t2!=0 &&
            sign(s1)!=sign(s2) && sign(t1)!=sign(t2) && return true
    end
    return false
end

# Exact open-segment/literal-triangle relation: 1 = proper crossing through
# the triangle interior, -1 = degenerate contact (an endpoint in the plane on
# the triangle, a grazing edge hit, or a coplanar intersection), 0 = no
# intersection.
function _model_native_hybrid_segment_triangle_relation(p,q,a,b,c)
    sp=orient3(a,b,c,p);sq=orient3(a,b,c,q)
    (isfinite(sp) && isfinite(sq)) || return -1
    if sp==0 || sq==0
        sp==0 && sq==0 &&
            return _model_native_hybrid_coplanar_segment_triangle(p,q,a,b,c) ? -1 : 0
        on=sp==0 ? p : q
        return _model_native_hybrid_coplanar_on_triangle(on,a,b,c) ? -1 : 0
    end
    sign(sp)==sign(sq) && return 0
    s1=orient3(p,q,a,b);s2=orient3(p,q,b,c);s3=orient3(p,q,c,a)
    (isfinite(s1) && isfinite(s2) && isfinite(s3)) || return -1
    (s1==0 || s2==0 || s3==0) && return -1
    sign(s1)==sign(s2)==sign(s3) && return 1
    return 0
end

# The shell is certified convex exactly when every literal panel plane keeps
# the whole boundary strictly on one side — the classic halfspace test the
# original certificates used. The return is one `(a,b,c,side)` entry per
# literal fan triangle, or `nothing` for a nonconvex shell, which falls back
# to the exact ray-parity and segment-clearance certificates below.
function _model_native_hybrid_convex_shell(panels,coordinate,boundary,caller,label)
    entries=NamedTuple[]
    for panel in panels,triangle in _model_native_hybrid_panel_triangles(panel)
        a,b,c=map(coordinate,triangle)
        sides=[orient3(a,b,c,coordinate(uid)) for uid in boundary]
        all(isfinite,sides) || throw(ArgumentError(
            "$caller: $label shell halfspaces are not finite"))
        signs=Set(sign(side) for side in sides if side!=0)
        length(signs)==1 || return nothing
        push!(entries,(a=a,b=b,c=c,side=only(signs)))
    end
    return entries
end

# Strict interior on a nonconvex shell: cast an exact ray to a point far past
# the axis-aligned bounding box and count literal triangle crossings. A ray
# that grazes a vertex, edge or panel plane retries with the next
# deterministic direction; an entity coplanar with every tried direction is
# genuinely on the boundary and rejects there.
function _model_native_hybrid_shell_inside(point,panels,coordinate,caller,label)
    bounds=[coordinate(uid) for panel in panels for uid in panel.nodes]
    span=max(1.0,maximum(maximum(b[axis] for b in bounds)-
                         minimum(b[axis] for b in bounds) for axis in 1:3))
    directions=((1.0,0.37139067635410372,0.19518001458970673),
                (0.41421356237309503,1.0,0.27182818284590454),
                (0.23606797749978981,0.57735026918962584,1.0),
                (1.0,1.0,1.0),
                (1.0,-0.70710678118654752,0.61803398874989490),
                (-0.61803398874989490,1.0,-0.41421356237309515))
    for direction in directions
        far=ntuple(axis->point[axis]+direction[axis]*span*4,3)
        crossings=0;degenerate=false
        for panel in panels,triangle in _model_native_hybrid_panel_triangles(panel)
            a,b,c=map(coordinate,triangle)
            relation=_model_native_hybrid_segment_triangle_relation(point,far,a,b,c)
            if relation<0
                degenerate=true;break
            end
            crossings+=relation
        end
        degenerate && continue
        isodd(crossings) || throw(ArgumentError(
            "$caller: $label is outside or on the hybrid shell boundary"))
        return nothing
    end
    throw(ArgumentError("$caller: $label is outside or on the hybrid shell boundary"))
end

# On a certified convex shell a strictly-interior entity stays on the inward
# halfspace side of every literal panel.
function _model_native_hybrid_convex_inside(records,entries,caller,label)
    for entry in entries
        all(record->begin
            value=orient3(entry.a,entry.b,entry.c,record.coordinate)
            isfinite(value) && value!=0 && sign(value)==entry.side
        end,records) || throw(ArgumentError(
            "$caller: $label is outside or on the hybrid shell boundary"))
    end
    return nothing
end

# On a nonconvex shell each literal segment must miss every panel interior:
# strictly-inside endpoints plus no proper crossing and no degenerate contact
# certify the open chain inside.
function _model_native_hybrid_shell_segments_inside(segments,lookup,panels,coordinate,caller,label)
    for segment in segments
        p,q=map(lookup,segment)
        for panel in panels,triangle in _model_native_hybrid_panel_triangles(panel)
            a,b,c=map(coordinate,triangle)
            _model_native_hybrid_segment_triangle_relation(p,q,a,b,c)==0 || throw(
                ArgumentError("$caller: $label is outside or on the hybrid shell boundary"))
        end
    end
    return nothing
end

_model_native_hybrid_shell_boundary(panels)=
    Set(uid for panel in panels for uid in panel.nodes)

function _model_native_hybrid_point_containment_certificate(constructor,caller)
    witness=constructor.embedded_point;witness===nothing && return nothing
    _model_native_hybrid_point_source_certificate(witness,constructor.input,
        constructor.node_uids,constructor.coordinates,caller)
    positions=Dict(uid=>row for (row,uid) in enumerate(constructor.node_uids))
    panels=vcat(constructor.triangles,constructor.quads)
    boundary=_model_native_hybrid_shell_boundary(panels)
    witness.record.uid in boundary && throw(ArgumentError(
        "$caller: embedded Point must be strictly inside the hybrid shell"))
    coordinate(uid)=Tuple(constructor.coordinates[:,positions[uid]])
    entries=_model_native_hybrid_convex_shell(panels,coordinate,boundary,caller,"embedded Point")
    if entries!==nothing
        _model_native_hybrid_convex_inside((witness.record,),entries,caller,"embedded Point")
    else
        _model_native_hybrid_shell_inside(witness.record.coordinate,panels,coordinate,
                                          caller,"embedded Point")
    end
    return nothing
end

function _model_native_hybrid_embedded_point_certificate(mesh,uids,constructor,caller)
    witness=constructor.embedded_point;witness===nothing && return nothing
    _model_native_hybrid_point_source_certificate(witness,constructor.input,uids,mesh.coords,caller)
    _model_native_hybrid_point_containment_certificate(constructor,caller)
    point=Int32(findfirst(==(witness.record.uid),uids))
    return _model_native_hybrid_interior_node_star_certificate(mesh,point,caller,"embedded Point")
end

# A referenced interior lower vertex must have one closed spherical Tet link.
function _model_native_hybrid_interior_node_star_certificate(mesh,point,caller,label)
    links=Set{NTuple{3,Int32}}()
    for block in mesh.blocks,column in axes(block.nodes,2)
        nodes=block.nodes[:,column];point in nodes || continue
        block.msh==4 || throw(ArgumentError(
            "$caller: $label has a foreign hybrid cell role"))
        face=Tuple(sort!(Int32[node for node in nodes if node!=point]))
        face in links && throw(ArgumentError("$caller: $label repeats a Tet star face"))
        push!(links,face)
    end
    isempty(links) && throw(ArgumentError(
        "$caller: $label UID is unreferenced by the generated volume cells"))
    edges=Dict{Tuple{Int32,Int32},Vector{Int}}();vertices=Set{Int32}();faces=collect(links)
    for (index,face) in enumerate(faces)
        union!(vertices,face)
        for edge in ((face[1],face[2]),(face[1],face[3]),(face[2],face[3]))
            push!(get!(()->Int[],edges,edge),index)
        end
    end
    all(indices->length(indices)==2,values(edges)) &&
        length(vertices)-length(edges)+length(faces)==2 || throw(ArgumentError(
        "$caller: $label lacks a closed interior Tet star"))
    reached=Set([1]);queue=[1];cursor=1
    while cursor<=length(queue)
        current=queue[cursor];cursor+=1
        for indices in values(edges)
            current in indices || continue
            for neighbor in indices
                neighbor in reached || (push!(reached,neighbor);push!(queue,neighbor))
            end
        end
    end
    length(reached)==length(faces) || throw(ArgumentError(
        "$caller: $label has disconnected Tet star components"))
    return nothing
end
