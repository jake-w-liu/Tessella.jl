"""
    BoundaryLayer

First-order boundary-layer element topology: prismatic extrusion of a
triangle surface along area-weighted vertex normals (type-6 prisms), planar
planar-polyline extrusion to type-3 quads with optional convex-corner fans
(type-2 triangles in the first layer, type-3 quads in subsequent layers),
and filled extrusion where the remaining core after a layer extrusion is
tetrahedralized behind a certified conforming prism/tetrahedron interface.
This is the element-topology counterpart of [`BoundaryLayerField`](@ref).
"""
module BoundaryLayer

using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets, triangle_area, tet_volume,
                   boundary_faces
using ..Elements: ElementBlock, MixedMesh, validate
using ..Predicates: orient2, orient3
using ..Mesh3D: delaunay3d, to_mesh3, recover_segment3, recover_triangle3,
                _raygrid, _inside_grid
using ..RecoverCDT: recover_boundary_cdt

export mesh_boundary_layer, mesh_boundary_layer_2d, mesh_boundary_layer_filled

function _finite(value, caller, name)
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    v=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    isfinite(v) || throw(ArgumentError("$caller: $name must be finite"))
    return v
end

function _bounded_int(value::Integer, caller, name; minimum::Int=0)
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    result=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name exceeds the platform Int range"))
    end
    result>=minimum || throw(ArgumentError(
        "$caller: $name must be ≥ $minimum"))
    return result
end

function _checked_add(a::Int,b::Int,caller,name)
    (a>=0 && b>=0) || throw(ArgumentError(
        "$caller: internal negative operand while computing $name"))
    a<=typemax(Int)-b || throw(ArgumentError(
        "$caller: $name exceeds the platform Int range"))
    return a+b
end

function _checked_mul(a::Int,b::Int,caller,name)
    (a>=0 && b>=0) || throw(ArgumentError(
        "$caller: internal negative operand while computing $name"))
    (a==0 || b<=typemax(Int)÷a) || throw(ArgumentError(
        "$caller: $name exceeds the platform Int range"))
    return a*b
end

function _validate_input(mesh::Mesh,caller,kind)
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: input $kind is invalid — "*join(diagnostic.messages,"; ")))
    ntets(mesh)==0 || throw(ArgumentError(
        "$caller: input $kind must not contain tetrahedra"))
    return nothing
end

function _require_all_referenced(cells,nn::Int,caller,kind)
    referenced=falses(nn)
    @inbounds for cell in axes(cells,2),slot in axes(cells,1)
        referenced[Int(cells[slot,cell])]=true
    end
    missing=findfirst(!,referenced)
    missing===nothing || throw(ArgumentError(
        "$caller: input $kind node $missing is not referenced"))
    return nothing
end

function _bounded_index_set(values,upper::Int,caller,name)
    applicable(iterate,values) || throw(ArgumentError(
        "$caller: $name must be an iterable of integer indices"))
    result=Set{Int}()
    count=0
    for raw in values
        count=_checked_add(count,1,caller,"$name count")
        count<=upper || throw(ArgumentError(
            "$caller: $name contains more than $upper entries"))
        raw isa Bool && throw(ArgumentError(
            "$caller: $name entry $count must not be Bool"))
        raw isa Integer || throw(ArgumentError(
            "$caller: $name entry $count must be an integer"))
        index=try
            Int(raw)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: $name entry $count exceeds the platform Int range"))
        end
        1<=index<=upper || throw(ArgumentError(
            "$caller: $name index $index is out of range 1:$upper"))
        index in result && throw(ArgumentError(
            "$caller: duplicate $name index $index"))
        push!(result,index)
    end
    return result
end

function _layer_offsets(hw, ra, nl, caller)
    offsets=Vector{Float64}(undef,nl)
    width=hw
    offset=0.0
    @inbounds for k in 1:nl
        offset+=width
        (isfinite(offset) && offset>0) || throw(ArgumentError(
            "$caller: layer offset $k is not finite and positive"))
        offsets[k]=offset
        if k<nl
            width*=ra
            (isfinite(width) && width>0) || throw(ArgumentError(
                "$caller: layer width $(k+1) is not finite and positive"))
        end
    end
    return offsets
end

"""
    mesh_boundary_layer(surface; hwall, ratio, nlayers,
                        max_prisms=10_000_000) -> MixedMesh

Extrude every triangle of a valid surface mesh through `nlayers` first-order
type-6 prisms along area-weighted vertex normals. The first-layer width is
`hwall`; each subsequent width is multiplied by `ratio`, which must be greater
than one. All surface nodes must be referenced by triangles, and the input must
not contain tetrahedra.

The operation checks all numeric conversions and output counts before
allocation, leaves `surface` unchanged, and returns a structurally validated
mixed mesh. Use [`mesh_boundary_layer_filled`](@ref) when the remaining enclosed
volume must also be tetrahedralized.
"""
function mesh_boundary_layer(surface::Mesh; hwall::Real, ratio::Real, nlayers::Integer,
                             max_prisms::Integer=10_000_000)
    caller="mesh_boundary_layer"
    ntris(surface)>0 || throw(ArgumentError("$caller: surface has no triangles"))
    _validate_input(surface,caller,"surface")
    _require_all_referenced(surface.tris,nnodes(surface),caller,"surface")
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=_bounded_int(nlayers,caller,"nlayers";minimum=1)
    prism_limit=_bounded_int(max_prisms,caller,"max_prisms")
    npr=_checked_mul(nl,ntris(surface),caller,"prism count")
    npr<=prism_limit || throw(ArgumentError(
        "$caller: $npr prisms exceed max_prisms=$prism_limit"))
    npr<=typemax(Int32) || throw(ArgumentError("$caller: prism count exceeds Int32"))

    nv=nnodes(surface)
    normals=zeros(Float64,3,nv)
    @inbounds for t in 1:ntris(surface)
        i,j,k=Int(surface.tris[1,t]),Int(surface.tris[2,t]),Int(surface.tris[3,t])
        a=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
        b=(surface.coords[1,j],surface.coords[2,j],surface.coords[3,j])
        c=(surface.coords[1,k],surface.coords[2,k],surface.coords[3,k])
        ab=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
        ac=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
        n=(ab[2]*ac[3]-ab[3]*ac[2], ab[3]*ac[1]-ab[1]*ac[3], ab[1]*ac[2]-ab[2]*ac[1])
        area=triangle_area(a,b,c)
        for id in (i,j,k)
            normals[1,id]+=n[1]*area; normals[2,id]+=n[2]*area; normals[3,id]+=n[3]*area
        end
    end
    @inbounds for i in 1:nv
        L=hypot(normals[1,i],normals[2,i],normals[3,i])
        L>0 || throw(ArgumentError("$caller: vertex $i has a zero normal"))
        normals[1,i]/=L; normals[2,i]/=L; normals[3,i]/=L
    end

    offsets=_layer_offsets(hw,ra,nl,caller)

    layer_count=_checked_add(nl,1,caller,"layer-node multiplier")
    nout=_checked_mul(nv,layer_count,caller,"node count")
    nout<=typemax(Int32) || throw(ArgumentError("$caller: node count exceeds Int32"))
    coords=Matrix{Float64}(undef,3,nout)
    @inbounds for i in 1:nv
        coords[1,i]=surface.coords[1,i]; coords[2,i]=surface.coords[2,i]; coords[3,i]=surface.coords[3,i]
    end
    @inbounds for k in 1:nl, i in 1:nv
        id=k*nv+i
        coords[1,id]=surface.coords[1,i]+offsets[k]*normals[1,i]
        coords[2,id]=surface.coords[2,i]+offsets[k]*normals[2,i]
        coords[3,id]=surface.coords[3,i]+offsets[k]*normals[3,i]
        all(isfinite, (coords[1,id],coords[2,id],coords[3,id])) ||
            throw(ArgumentError("$caller: extruded node is non-finite"))
    end

    prisms=Matrix{Int32}(undef,6,npr)
    cursor=0
    @inbounds for k in 0:nl-1, t in 1:ntris(surface)
        cursor+=1
        b1=Int32(k*nv+Int(surface.tris[1,t]))
        b2=Int32(k*nv+Int(surface.tris[2,t]))
        b3=Int32(k*nv+Int(surface.tris[3,t]))
        t1=Int32((k+1)*nv+Int(surface.tris[1,t]))
        t2=Int32((k+1)*nv+Int(surface.tris[2,t]))
        t3=Int32((k+1)*nv+Int(surface.tris[3,t]))
        prisms[:,cursor].=(b1,b2,b3,t1,t2,t3)
    end
    mesh=MixedMesh(coords,[ElementBlock(6,prisms)])
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return mesh
end

function _left_unit(ax, ay, bx, by, caller, seg)
    tx,ty=bx-ax,by-ay
    (isfinite(tx) && isfinite(ty)) || throw(ArgumentError(
        "$caller: segment $seg coordinate span overflows Float64"))
    L=hypot(tx,ty)
    (isfinite(L) && L>0) || throw(ArgumentError(
        "$caller: segment $seg has zero or non-finite length"))
    return (-ty/L, tx/L)
end

function _plane_unit_normal(raw, caller)
    (raw isa Tuple || raw isa AbstractVector) || throw(ArgumentError(
        "$caller: plane_normal must be a three-component tuple or vector"))
    component_count=try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: plane_normal must have a finite declared length"))
    end
    component_count==3 || throw(ArgumentError(
        "$caller: plane_normal must have exactly three components"))
    values=try
        ntuple(i->raw[i],3)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: plane_normal components must be indexable"))
    end
    normal=ntuple(i->_finite(values[i],caller,"plane_normal[$i]"),3)
    normal_scale=max(abs(normal[1]),abs(normal[2]),abs(normal[3]))
    normal_scale>0 || throw(ArgumentError(
        "$caller: plane_normal must have finite positive length"))
    scaled=ntuple(i->normal[i]/normal_scale,3)
    length_scaled=hypot(scaled...)
    (isfinite(length_scaled) && length_scaled>0) || throw(ArgumentError(
        "$caller: plane_normal cannot be normalized"))
    unit=ntuple(i->scaled[i]/length_scaled,3)
    all(isfinite,unit) || throw(ArgumentError(
        "$caller: normalized plane_normal is not finite"))
    return unit
end

function _plane_basis(normal,caller)
    ax,ay,az=abs(normal[1]),abs(normal[2]),abs(normal[3])
    # Prefer y on ties so the default +z normal retains the historical x/y
    # coordinate frame and therefore its deterministic output bytes.
    reference=ay<=ax && ay<=az ? (0.0,1.0,0.0) :
              ax<=az ? (1.0,0.0,0.0) : (0.0,0.0,1.0)
    raw_u=_cross3(reference,normal)
    length_u=hypot(raw_u...)
    (isfinite(length_u) && length_u>0) || throw(ArgumentError(
        "$caller: could not construct a stable plane frame"))
    u=ntuple(i->raw_u[i]/length_u,3)
    raw_v=_cross3(normal,u)
    length_v=hypot(raw_v...)
    (isfinite(length_v) && length_v>0) || throw(ArgumentError(
        "$caller: could not complete a stable plane frame"))
    v=ntuple(i->raw_v[i]/length_v,3)
    all(isfinite,u) && all(isfinite,v) || throw(ArgumentError(
        "$caller: plane frame is not finite"))
    return u,v
end

function _validate_plane(points,normal,caller)
    origin=points[1]
    scale=0.0
    maximum_distance=0.0
    maximum_node=1
    @inbounds for (index,point) in pairs(points)
        delta=(point[1]-origin[1],point[2]-origin[2],point[3]-origin[3])
        all(isfinite,delta) || throw(ArgumentError(
            "$caller: coordinate span overflows Float64 at node $index"))
        span=hypot(delta...)
        isfinite(span) || throw(ArgumentError(
            "$caller: coordinate span is not finite at node $index"))
        scale=max(scale,span)
        distance=abs(_dot3(delta,normal))
        isfinite(distance) || throw(ArgumentError(
            "$caller: plane distance is not finite at node $index"))
        if distance>maximum_distance
            maximum_distance=distance
            maximum_node=index
        end
    end
    tolerance=8192eps(Float64)*max(scale,1.0)
    maximum_distance<=tolerance || throw(ArgumentError(
        "$caller: node $maximum_node is outside plane_normal's plane "*
        "(distance $maximum_distance exceeds $tolerance)"))
    return origin
end

@inline function _plane_point(point,origin,u,v)
    delta=(point[1]-origin[1],point[2]-origin[2],point[3]-origin[3])
    return (_dot3(delta,u),_dot3(delta,v))
end

@inline function _offset_point(point,u,v,dx,dy,offset)
    direction=(dx*u[1]+dy*v[1],dx*u[2]+dy*v[2],dx*u[3]+dy*v[3])
    return (point[1]+offset*direction[1],
            point[2]+offset*direction[2],
            point[3]+offset*direction[3])
end

function _polyline_vertices(curve::Mesh, caller)
    n=nsegs(curve)
    n>0 || throw(ArgumentError("$caller: curve has no segments"))
    nv=nnodes(curve)
    next_vertex=Dict{Int,Int}()
    indegree=zeros(Int,nv)
    referenced=falses(nv)
    @inbounds for s in 1:n
        a=Int(curve.segs[1,s]); b=Int(curve.segs[2,s])
        (1<=a<=nv && 1<=b<=nv) || throw(ArgumentError("$caller: segment $s is out of range"))
        a==b && throw(ArgumentError("$caller: segment $s is degenerate"))
        haskey(next_vertex,a) && throw(ArgumentError(
            "$caller: vertex $a has more than one outgoing segment"))
        indegree[b]=_checked_add(indegree[b],1,caller,"vertex $b indegree")
        indegree[b]<=1 || throw(ArgumentError(
            "$caller: vertex $b has more than one incoming segment"))
        next_vertex[a]=b
        referenced[a]=true;referenced[b]=true
    end
    missing=findfirst(!,referenced)
    missing===nothing || throw(ArgumentError(
        "$caller: input curve node $missing is not referenced"))
    starts=Int[];ends=Int[]
    @inbounds for v in 1:nv
        outgoing=haskey(next_vertex,v)
        indegree[v]==0 && outgoing && push!(starts,v)
        indegree[v]==1 && !outgoing && push!(ends,v)
        (indegree[v]<=1 && (outgoing || indegree[v]==1)) || throw(ArgumentError(
            "$caller: segments are not a coherently directed single polyline"))
    end
    if length(starts)==1 && length(ends)==1
        closed=false
        start=only(starts)
    elseif isempty(starts) && isempty(ends) && length(next_vertex)==nv
        closed=true
        start=minimum(keys(next_vertex))
    else
        throw(ArgumentError(
            "$caller: segments are not a coherently directed single polyline"))
    end
    verts=Int[]
    seen=falses(nv)
    current=start
    while true
        seen[current] && throw(ArgumentError(
            "$caller: segments are not a single directed polyline"))
        push!(verts,current);seen[current]=true
        nxt=get(next_vertex,current,0)
        nxt==0 && break
        if closed && nxt==start
            break
        end
        current=nxt
        length(verts)>n+1 && throw(ArgumentError("$caller: segments are not a single polyline"))
    end
    expected=closed ? n : n+1
    (length(verts)==expected && all(seen)) || throw(ArgumentError(
        "$caller: segments are not a single directed polyline"))
    return verts, closed
end

@inline function _positive_quad(a,b,c,d)
    return orient2(a,b,c)>0 && orient2(b,c,d)>0 &&
           orient2(c,d,a)>0 && orient2(d,a,b)>0
end

function _id_layer(id_of, v, layer, ray, fan)
    return layer==0 ? id_of[(v,0,0)] : id_of[(v,layer, fan ? ray : 0)]
end

"""
    mesh_boundary_layer_2d(curve; hwall, ratio, nlayers, fans=(),
                           fan_elements=5, plane_normal=(0.0,0.0,1.0),
                           max_cells=10_000_000) -> MixedMesh

Extrude a valid, coherently directed planar polyline along its left normals into
first-order type-3 quadrangles. `plane_normal` gives the oriented normal
direction used to define "left"; it defaults to positive `z` for backward
compatibility. The curve must lie in the plane through its first node orthogonal
to that normal, within a scale-aware floating tolerance.

`fans` lists convex interior vertex indices that receive `fan_elements`
first-layer triangles and matching ring quadrangles instead of a single
averaged-normal column. Every curve node must belong to the chain,
`fan_elements` must be at least two, and `max_cells` is a nonnegative allocation
bound. The input is left unchanged.
"""
function mesh_boundary_layer_2d(curve::Mesh; hwall::Real, ratio::Real, nlayers::Integer,
                                fans=(), fan_elements::Integer=5,
                                plane_normal=(0.0,0.0,1.0),
                                max_cells::Integer=10_000_000)
    caller="mesh_boundary_layer_2d"
    nsegs(curve)>0 || throw(ArgumentError("$caller: curve has no segments"))
    (ntris(curve)==0 && ntets(curve)==0) || throw(ArgumentError(
        "$caller: input curve must contain only segment cells"))
    _validate_input(curve,caller,"curve")
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=_bounded_int(nlayers,caller,"nlayers";minimum=1)
    nfan=_bounded_int(fan_elements,caller,"fan_elements";minimum=2)
    normal=_plane_unit_normal(plane_normal,caller)
    cell_limit=_bounded_int(max_cells,caller,"max_cells")
    nv=nnodes(curve)
    nv>=2 || throw(ArgumentError("$caller: curve has too few nodes"))
    verts, closed=_polyline_vertices(curve,caller)
    nchain=length(verts)
    nchain_segs=closed ? nchain : nchain-1
    fan_set=_bounded_index_set(fans,nv,caller,"fan")
    pos=Dict{Int,Int}()
    for (idx,v) in enumerate(verts)
        pos[v]=idx
    end
    for v in fan_set
        haskey(pos,v) || throw(ArgumentError("$caller: fan vertex $v is not on the polyline"))
        idx=pos[v]
        interior=closed || (idx>1 && idx<nchain)
        interior || throw(ArgumentError("$caller: fan vertex $v is not an interior vertex"))
    end
    n_fans=length(fan_set)
    ntris_out=_checked_mul(n_fans,nfan,caller,"fan triangle count")
    main_quads=_checked_mul(nchain_segs,nl,caller,"strip quadrangle count")
    ring_quads=_checked_mul(ntris_out,nl-1,caller,"fan ring quadrangle count")
    nquads=_checked_add(main_quads,ring_quads,caller,"quadrangle count")
    ncells=_checked_add(ntris_out,nquads,caller,"cell count")
    ncells<=cell_limit || throw(ArgumentError(
        "$caller: $ncells cells exceed max_cells=$cell_limit"))
    ncells<=typemax(Int32) || throw(ArgumentError("$caller: cell count exceeds Int32"))
    ncells>0 || throw(ArgumentError("$caller: no cells to emit"))
    points_per_layer=_checked_add(nchain,ntris_out,caller,"points per layer")
    new_points=_checked_mul(nl,points_per_layer,caller,"extruded node count")
    npoints=_checked_add(nv,new_points,caller,"node count")
    npoints<=typemax(Int32) || throw(ArgumentError("$caller: node count exceeds Int32"))

    points=Vector{NTuple{3,Float64}}(undef,nv)
    @inbounds for i in 1:nv
        points[i]=(curve.coords[1,i],curve.coords[2,i],curve.coords[3,i])
    end
    plane_origin=_validate_plane(points,normal,caller)
    basis_u,basis_v=_plane_basis(normal,caller)
    plane_points=Vector{NTuple{2,Float64}}(undef,nv)
    @inbounds for i in 1:nv
        plane_points[i]=_plane_point(points[i],plane_origin,basis_u,basis_v)
        all(isfinite,plane_points[i]) || throw(ArgumentError(
            "$caller: projected node $i is not finite"))
    end
    seg_left=Vector{NTuple{2,Float64}}(undef,nchain_segs)
    @inbounds for i in 1:nchain_segs
        a=verts[i]; b=verts[i==nchain ? 1 : i+1]
        seg_left[i]=_left_unit(
            plane_points[a][1],plane_points[a][2],
            plane_points[b][1],plane_points[b][2],caller,i)
    end
    n_in=Vector{NTuple{2,Float64}}(undef,nchain)
    n_out=Vector{NTuple{2,Float64}}(undef,nchain)
    n_avg=Vector{NTuple{2,Float64}}(undef,nchain)
    @inbounds for i in 1:nchain
        if closed
            nin=seg_left[i==1 ? nchain_segs : i-1]
            nout=seg_left[i]
        elseif i==1
            nin=nout=seg_left[1]
        elseif i==nchain
            nin=nout=seg_left[nchain_segs]
        else
            nin=seg_left[i-1]; nout=seg_left[i]
        end
        n_in[i]=nin; n_out[i]=nout
        sx,sy=nin[1]+nout[1], nin[2]+nout[2]
        L=hypot(sx,sy)
        L>0 || throw(ArgumentError("$caller: vertex $(verts[i]) has a zero normal"))
        n_avg[i]=(sx/L,sy/L)
    end
    for v in fan_set
        i=pos[v]
        nin=n_in[i]; nout=n_out[i]
        cr=nin[1]*nout[2]-nin[2]*nout[1]
        dt=nin[1]*nout[1]+nin[2]*nout[2]
        θ=atan(cr,dt)
        θ>0 || throw(ArgumentError(
            "$caller: fan vertex $v does not have a positive CCW sector"))
    end

    offsets=_layer_offsets(hw,ra,nl,caller)
    id_of=Dict{Tuple{Int,Int,Int},Int}()
    @inbounds for v in verts
        id_of[(v,0,0)]=v
    end
    @inbounds for k in 1:nl, (idx,v) in enumerate(verts)
        if v in fan_set
            nin=n_in[idx]; nout=n_out[idx]
            θ=atan(nin[1]*nout[2]-nin[2]*nout[1], nin[1]*nout[1]+nin[2]*nout[2])
            for ray in 0:nfan
                α=ray/nfan*θ
                c,s=cos(α),sin(α)
                dx=nin[1]*c-nin[2]*s
                dy=nin[1]*s+nin[2]*c
                new_point=_offset_point(
                    points[v],basis_u,basis_v,dx,dy,offsets[k])
                all(isfinite,new_point) || throw(ArgumentError(
                    "$caller: extruded node is non-finite"))
                new_plane_point=_plane_point(
                    new_point,plane_origin,basis_u,basis_v)
                all(isfinite,new_plane_point) || throw(ArgumentError(
                    "$caller: projected extruded node is non-finite"))
                push!(points,new_point)
                push!(plane_points,new_plane_point)
                id_of[(v,k,ray)]=length(points)
            end
        else
            nx,ny=n_avg[idx]
            new_point=_offset_point(
                points[v],basis_u,basis_v,nx,ny,offsets[k])
            all(isfinite,new_point) || throw(ArgumentError(
                "$caller: extruded node is non-finite"))
            new_plane_point=_plane_point(
                new_point,plane_origin,basis_u,basis_v)
            all(isfinite,new_plane_point) || throw(ArgumentError(
                "$caller: projected extruded node is non-finite"))
            push!(points,new_point)
            push!(plane_points,new_plane_point)
            id_of[(v,k,0)]=length(points)
        end
    end
    length(points)==npoints || throw(ErrorException("$caller: node count mismatch"))
    length(plane_points)==npoints || throw(ErrorException(
        "$caller: projected node count mismatch"))
    _validate_plane(points,normal,caller)

    quads=Matrix{Int32}(undef,4,nquads)
    tris=ntris_out==0 ? Matrix{Int32}(undef,3,0) : Matrix{Int32}(undef,3,ntris_out)
    qcursor=0; tcursor=0
    @inbounds for i in 1:nchain_segs
        a=verts[i]; b=verts[i==nchain ? 1 : i+1]
        a_fan=a in fan_set; b_fan=b in fan_set
        ray_a=a_fan ? nfan : 0
        ray_b=0
        for k in 1:nl
            b1=_id_layer(id_of,a,k-1,ray_a,a_fan)
            b2=_id_layer(id_of,b,k-1,ray_b,b_fan)
            t2=_id_layer(id_of,b,k,ray_b,b_fan)
            t1=_id_layer(id_of,a,k,ray_a,a_fan)
            p1=plane_points[b1]; p2=plane_points[b2]
            p3=plane_points[t2]; p4=plane_points[t1]
            _positive_quad(p1,p2,p3,p4) || throw(ArgumentError(
                "$caller: quadrangle on segment $i layer $k has a zero or "*
                "reversed corner Jacobian"))
            qcursor+=1
            quads[:,qcursor].=(Int32(b1),Int32(b2),Int32(t2),Int32(t1))
        end
    end
    @inbounds for (idx,v) in enumerate(verts)
        v in fan_set || continue
        origin=id_of[(v,0,0)]
        for ray in 0:nfan-1
            i1=id_of[(v,1,ray)]; i2=id_of[(v,1,ray+1)]
            p0=plane_points[origin]; p1=plane_points[i1]; p2=plane_points[i2]
            orient2(p0,p1,p2)>0 ||
                throw(ArgumentError("$caller: inverted fan triangle at vertex $v"))
            tcursor+=1
            tris[:,tcursor].=(Int32(origin),Int32(i1),Int32(i2))
        end
        for k in 2:nl, ray in 0:nfan-1
            b1=id_of[(v,k-1,ray)]; b2=id_of[(v,k-1,ray+1)]
            t2=id_of[(v,k,ray+1)]; t1=id_of[(v,k,ray)]
            p1=plane_points[b1]; p2=plane_points[t1]
            p3=plane_points[t2]; p4=plane_points[b2]
            _positive_quad(p1,p2,p3,p4) || throw(ArgumentError(
                "$caller: fan quadrangle at vertex $v layer $k has a zero or "*
                "reversed corner Jacobian"))
            qcursor+=1
            quads[:,qcursor].=(Int32(b1),Int32(t1),Int32(t2),Int32(b2))
        end
    end
    qcursor==nquads || throw(ErrorException("$caller: quadrangle count mismatch"))
    tcursor==ntris_out || throw(ErrorException("$caller: triangle count mismatch"))

    coords=Matrix{Float64}(undef,3,length(points))
    @inbounds for i in eachindex(points)
        coords[1,i]=points[i][1]; coords[2,i]=points[i][2]; coords[3,i]=points[i][3]
    end
    blocks=ElementBlock[]
    ntris_out>0 && push!(blocks,ElementBlock(2,tris))
    push!(blocks,ElementBlock(3,quads))
    mesh=MixedMesh(coords,blocks)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return mesh
end

# Closed-manifold wall decomposition: union-find over triangles joined by shared
# edges. Returns (comp_of, comp_tris, nc) where comp_of[v] is the 1-based wall
# index of vertex v and comp_tris[c] lists that wall's triangles. Blocks on
# unreferenced nodes, non-closed edges, and pinched vertices (a vertex shared by
# edge-disjoint walls), which have no single well-defined extrusion normal.
function _wall_components(surface::Mesh, caller)
    nv=nnodes(surface); nt=ntris(surface)
    nt>0 || throw(ArgumentError("$caller: surface has no triangles"))
    referenced=falses(nv)
    inc=Dict{NTuple{2,Int32},Int}()
    edge_tri=Dict{NTuple{2,Int32},Int32}()
    tparent=Vector{Int32}(undef,nt)
    @inbounds for t in 1:nt; tparent[t]=Int32(t); end
    tfind!(a::Int32) = begin
        while tparent[a]!=a
            tparent[a]=tparent[tparent[a]]
            a=tparent[a]
        end
        return a
    end
    @inbounds for f in 1:nt
        a=surface.tris[1,f]; b=surface.tris[2,f]; c=surface.tris[3,f]
        (1<=a<=nv && 1<=b<=nv && 1<=c<=nv) ||
            throw(ArgumentError("$caller: triangle $f references a node out of range"))
        (a==b || b==c || a==c) &&
            throw(ArgumentError("$caller: triangle $f is degenerate"))
        referenced[a]=true; referenced[b]=true; referenced[c]=true
        for (u,v) in ((a,b),(b,c),(c,a))
            key=u<v ? (u,v) : (v,u)
            inc[key]=get(inc,key,0)+1
            prev=get(edge_tri,key,Int32(0))
            if prev==0
                edge_tri[key]=Int32(f)
            else
                ra=tfind!(prev); rb=tfind!(Int32(f))
                ra!=rb && (tparent[ra]=rb)
            end
        end
    end
    @inbounds for v in 1:nv
        referenced[v] || throw(ArgumentError(
            "$caller: node $v is not referenced by any triangle"))
    end
    for (e,n) in inc
        n==2 || throw(ArgumentError(
            "$caller: wall edge $e has incidence $n; every wall must be closed and manifold"))
    end
    vwall=Dict{Int32,Int32}()
    @inbounds for f in 1:nt
        r=tfind!(Int32(f))
        for v in (surface.tris[1,f],surface.tris[2,f],surface.tris[3,f])
            prev=get(vwall,v,Int32(0))
            if prev==0
                vwall[v]=r
            elseif prev!=r
                throw(ArgumentError(
                    "$caller: node $v is pinched by edge-disjoint walls; " *
                    "extrusion has no single normal there"))
            end
        end
    end
    roots=Dict{Int32,Int}()
    comp_of=Vector{Int}(undef,nv)
    for (v,r) in sort!(collect(vwall); by=p->p[1])
        idx=get!(roots,r,length(roots)+1)
        comp_of[Int(v)]=idx
    end
    nc=length(roots)
    comp_tris=[Int[] for _ in 1:nc]
    @inbounds for f in 1:nt
        push!(comp_tris[comp_of[Int(surface.tris[1,f])]],f)
    end
    return comp_of, comp_tris, nc
end

# Signed enclosed volume by the divergence theorem over triangles as wound.
function _div_volume(coords, faces)
    total=0.0
    @inbounds for (i,j,k) in faces
        ax,ay,az=coords[1,i],coords[2,i],coords[3,i]
        bx,by,bz=coords[1,j],coords[2,j],coords[3,j]
        cx,cy,cz=coords[1,k],coords[2,k],coords[3,k]
        total += (ax*(by*cz-bz*cy)+ay*(bz*cx-bx*cz)+az*(bx*cy-by*cx))/6.0
    end
    return total
end

@inline _cross3(a,b)=(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _dot3(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

# Volume of one prism wedge from its six nodes (bottom tri then top tri). Side
# quads of a per-vertex-normal sweep are generally skew, so each is split along
# the canonical diagonal (smaller global id to smaller global id) — the same
# choice either neighbour wedge makes for their shared face, which makes the
# layer stack an exact partition. Subfaces are oriented away from the centroid
# with the exact orient3 predicate, so folded or flat wedges report ≤ 0.
function _prism_volume6(coords, v::NTuple{6,Int})
    pts=NTuple{3,Float64}[(coords[1,v[i]],coords[2,v[i]],coords[3,v[i]]) for i in 1:6]
    cen=(sum(p[1] for p in pts)/6.0,sum(p[2] for p in pts)/6.0,sum(p[3] for p in pts)/6.0)
    vol_face(a,b,c)=begin
        s=_dot3(a,_cross3(b,c))/6.0
        # Tessella's orient3 is negative on the right-hand side of (a,b,c), so an
        # outward-facing winding (centroid on the opposite side) reads positive.
        orient3(a,b,c,cen)>0 ? s : -s
    end
    # bottom cap reversed, top cap as wound
    total=vol_face(pts[1],pts[3],pts[2])+vol_face(pts[4],pts[5],pts[6])
    for (i,j) in ((1,2),(2,3),(3,1))
        bi=pts[i]; bj=pts[j]; ti=pts[i+3]; tj=pts[j+3]
        if v[i]<v[j]
            total+=vol_face(bi,bj,tj)+vol_face(bi,tj,ti)
        else
            total+=vol_face(bj,bi,ti)+vol_face(bj,ti,tj)
        end
    end
    return total
end

_norm_key(x,y,z)=((x==0 ? 0.0 : x),(y==0 ? 0.0 : y),(z==0 ? 0.0 : z))

# Ray-parity point-in-closed-surface test along one fixed direction.
function _ray_inside(p, dir, coords, faces)
    crossings=0
    @inbounds for (a,b,c) in faces
        pa=(coords[1,a],coords[2,a],coords[3,a])
        pb=(coords[1,b],coords[2,b],coords[3,b])
        pc=(coords[1,c],coords[2,c],coords[3,c])
        e1=(pb[1]-pa[1],pb[2]-pa[2],pb[3]-pa[3])
        e2=(pc[1]-pa[1],pc[2]-pa[2],pc[3]-pa[3])
        pv=_cross3(dir,e2)
        det=_dot3(e1,pv)
        abs(det)>1e-14 || continue
        inv=1.0/det
        s=(p[1]-pa[1],p[2]-pa[2],p[3]-pa[3])
        u=_dot3(s,pv)*inv
        (u<-1e-12 || u>1+1e-12) && continue
        q=_cross3(s,e1)
        v=_dot3(dir,q)*inv
        (v<-1e-12 || u+v>1+1e-12) && continue
        t=_dot3(e2,q)*inv
        t>1e-10 && (crossings+=1)
    end
    return isodd(crossings)
end

# Classify two walls as :nested_a_in_b, :nested_b_in_a, :disjoint, or :ambiguous.
# Three generic directions must agree per test; any disagreement is ambiguous.
function _wall_relation(acoords, afaces, bcoords, bfaces)
    dirs=((0.8131,-0.3742,0.4459),(-0.2667,0.8089,0.5241),(0.3559,0.4527,-0.8196))
    pa=(acoords[1,afaces[1][1]],acoords[2,afaces[1][1]],acoords[3,afaces[1][1]])
    pb=(bcoords[1,bfaces[1][1]],bcoords[2,bfaces[1][1]],bcoords[3,bfaces[1][1]])
    ain=sum(_ray_inside(pa,d,bcoords,bfaces) for d in dirs)
    bin=sum(_ray_inside(pb,d,acoords,afaces) for d in dirs)
    (ain==0||ain==3) && (bin==0||bin==3) || return :ambiguous
    if ain==3 && bin==0; return :nested_a_in_b; end
    if bin==3 && ain==0; return :nested_b_in_a; end
    ain==0 && bin==0 && return :disjoint
    return :ambiguous
end

"""
    mesh_boundary_layer_filled(surface; hwall, ratio, nlayers, cavities=(),
                               max_prisms=10_000_000, max_tets=10_000_000)
                               -> MixedMesh

Extrude every closed manifold wall of `surface` into first-order type-6 prism
layers along area-weighted vertex normals — **into** each solid wall's enclosed
interior, and **away from** it for walls listed in `cavities` (1-based wall
indices, i.e. interior holes) — then tetrahedralize the remaining core behind
the last layer and merge interface nodes onto the prism-stack numbering.

The core is filled by a two-stage ladder: an exact-coordinate Float64 Delaunay
path for Delaunay-friendly caps (flat/mostly-planar walls) under a hard Steiner
budget, then an exact-rational conforming-recovery pass for harder caps. Both
stages are held to the same independent certificates; whichever succeeds first
is returned.

The returned mesh is certified before return:

- every wall vertex survives recovery and merges into one shared node per cap
  position, so prisms and core tets conform with no cracks;
- the recovered tetrahedron boundary tiles the offset caps exactly (equal total
  area to 1e-9 relative);
- per-wall shell identities `|V(S)| ∓ V(prisms) == |V(cap)|` (minus for solids,
  plus for cavities) and the global fill identity
  `V(tets) == Σ V(solid caps) − Σ V(cavity caps)` hold to 1e-9 relative;
- every prism and tetrahedron has strictly positive exact-predicate volume.

Too-large `hwall` (self-intersecting layers), interpenetrating walls, pinched
or open walls, or an unmeshable core are explicit blockers — never defective
meshes. The input must be a valid triangle surface without tetrahedra; resource
limits and all integer controls are validated before output allocation.
"""
function mesh_boundary_layer_filled(surface::Mesh; hwall::Real, ratio::Real,
                                    nlayers::Integer, cavities=(),
                                    max_prisms::Integer=10_000_000,
                                    max_tets::Integer=10_000_000)
    caller="mesh_boundary_layer_filled"
    ntris(surface)>0 || throw(ArgumentError("$caller: surface has no triangles"))
    _validate_input(surface,caller,"surface")
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=_bounded_int(nlayers,caller,"nlayers";minimum=1)
    prism_limit=_bounded_int(max_prisms,caller,"max_prisms")
    tet_limit=_bounded_int(max_tets,caller,"max_tets")
    nv=nnodes(surface); nt=ntris(surface)

    comp_of,comp_tris,nc=_wall_components(surface,caller)
    cavity_set=_bounded_index_set(cavities,nc,caller,"cavity")
    nc>length(cavity_set) || throw(ArgumentError(
        "$caller: every wall is a cavity; there is no solid core to fill"))

    npr=_checked_mul(nl,nt,caller,"prism count")
    npr<=prism_limit || throw(ArgumentError(
        "$caller: $npr prisms exceed max_prisms=$prism_limit"))
    npr<=typemax(Int32) || throw(ArgumentError("$caller: prism count exceeds Int32"))
    layer_count=_checked_add(nl,1,caller,"layer-node multiplier")
    nout=_checked_mul(nv,layer_count,caller,"node count")
    nout<=typemax(Int32) || throw(ArgumentError("$caller: node count exceeds Int32"))

    offsets=_layer_offsets(hw,ra,nl,caller)

    # Per-wall enclosed volume with inherited winding (positive ⇔ outward).
    wall_faces=[NTuple{3,Int}[] for _ in 1:nc]
    @inbounds for f in 1:nt
        push!(wall_faces[comp_of[Int(surface.tris[1,f])]],
              (Int(surface.tris[1,f]),Int(surface.tris[2,f]),Int(surface.tris[3,f])))
    end
    wall_vol=[_div_volume(surface.coords,wall_faces[c]) for c in 1:nc]
    for c in 1:nc
        (isfinite(wall_vol[c]) && abs(wall_vol[c])>0) || throw(ArgumentError(
            "$caller: wall $c does not have a finite nonzero enclosed volume"))
    end

    # Solid walls grow against their winding normal (into the material);
    # cavity walls grow along it (into the surrounding material).
    dirs=[(c in cavity_set ? 1.0 : -1.0)*sign(wall_vol[c]) for c in 1:nc]

    # Area-weighted vertex normals scaled into each wall's growth direction.
    normals=zeros(Float64,3,nv)
    @inbounds for t in 1:nt
        i,j,k=Int(surface.tris[1,t]),Int(surface.tris[2,t]),Int(surface.tris[3,t])
        a=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
        b=(surface.coords[1,j],surface.coords[2,j],surface.coords[3,j])
        cc=(surface.coords[1,k],surface.coords[2,k],surface.coords[3,k])
        ab=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
        ac=(cc[1]-a[1],cc[2]-a[2],cc[3]-a[3])
        n=(ab[2]*ac[3]-ab[3]*ac[2], ab[3]*ac[1]-ab[1]*ac[3], ab[1]*ac[2]-ab[2]*ac[1])
        area=triangle_area(a,b,cc)
        for id in (i,j,k)
            normals[1,id]+=n[1]*area; normals[2,id]+=n[2]*area; normals[3,id]+=n[3]*area
        end
    end
    @inbounds for i in 1:nv
        L=hypot(normals[1,i],normals[2,i],normals[3,i])
        L>0 || throw(ArgumentError("$caller: vertex $i has a zero normal"))
        s=dirs[comp_of[i]]/L
        normals[1,i]*=s; normals[2,i]*=s; normals[3,i]*=s
    end

    coords=Matrix{Float64}(undef,3,nout)
    @inbounds for i in 1:nv
        coords[1,i]=surface.coords[1,i]; coords[2,i]=surface.coords[2,i]; coords[3,i]=surface.coords[3,i]
    end
    @inbounds for k in 1:nl, i in 1:nv
        id=k*nv+i
        d=offsets[k]
        coords[1,id]=surface.coords[1,i]+d*normals[1,i]
        coords[2,id]=surface.coords[2,i]+d*normals[2,i]
        coords[3,id]=surface.coords[3,i]+d*normals[3,i]
        all(isfinite,(coords[1,id],coords[2,id],coords[3,id])) ||
            throw(ArgumentError("$caller: extruded node is non-finite"))
    end

    # Prism layers share the plain mesh_boundary_layer layout and orientation.
    prisms=Matrix{Int32}(undef,6,npr)
    cursor=0
    @inbounds for k in 0:nl-1, t in 1:nt
        cursor+=1
        prisms[:,cursor].=(Int32(k*nv+Int(surface.tris[1,t])),
                           Int32(k*nv+Int(surface.tris[2,t])),
                           Int32(k*nv+Int(surface.tris[3,t])),
                           Int32((k+1)*nv+Int(surface.tris[1,t])),
                           Int32((k+1)*nv+Int(surface.tris[2,t])),
                           Int32((k+1)*nv+Int(surface.tris[3,t])))
    end
    cursor==npr || throw(ErrorException("$caller: prism count mismatch"))

    # Multi-wall separation gate. AABBs farther apart than both depths prove
    # the shells can never meet. Otherwise the pair is classified by ray parity:
    # genuinely nested walls (cavities) are allowed — cap crossings are caught
    # by the exact downstream certificates — while disjoint-but-close walls are
    # blocked, and ambiguous classifications are blocked conservatively.
    if nc>1
        boxes=Vector{NTuple{6,Float64}}(undef,nc)
        @inbounds for c in 1:nc
            lo=(Inf,Inf,Inf); hi=(-Inf,-Inf,-Inf)
            for f in wall_faces[c], i in f
                lo=(min(lo[1],surface.coords[1,i]),min(lo[2],surface.coords[2,i]),min(lo[3],surface.coords[3,i]))
                hi=(max(hi[1],surface.coords[1,i]),max(hi[2],surface.coords[2,i]),max(hi[3],surface.coords[3,i]))
            end
            boxes[c]=(lo[1],lo[2],lo[3],hi[1],hi[2],hi[3])
        end
        depth=offsets[end]
        for a in 1:nc, b in a+1:nc
            A=boxes[a]; B=boxes[b]
            gap=max(A[1]-B[4], B[1]-A[4], A[2]-B[5], B[2]-A[5], A[3]-B[6], B[3]-A[6])
            gap>2*depth && continue
            rel=_wall_relation(surface.coords,wall_faces[a],surface.coords,wall_faces[b])
            rel===:disjoint && throw(ArgumentError(
                "$caller: walls $a and $b approach within the layer depth; " *
                "reduce hwall so distinct shells never meet"))
            rel===:ambiguous && throw(ArgumentError(
                "$caller: walls $a and $b are too close to certify separation; " *
                "reduce hwall so distinct shells never meet"))
        end
    end

    # Cap surface at the last layer: an offset copy with inherited connectivity.
    cap_coords=coords[:,nl*nv+1:nl*nv+nv]
    cap=Mesh(cap_coords; tris=surface.tris)
    cap_off=nout-nv
    cap_keys=Dict{NTuple{3,Float64},Int}()
    @inbounds for i in 1:nv
        key=_norm_key(cap_coords[1,i],cap_coords[2,i],cap_coords[3,i])
        previous=get(cap_keys,key,0)
        previous==0 || throw(ArgumentError(
            "$caller: offset cap nodes $(previous-cap_off) and $i coincide"))
        cap_keys[key]=cap_off+i
    end

    # Stage 1 — exact-coordinate Float64 Delaunay with facet recovery. No SoS
    # perturbation, so wall vertices keep bit-exact coordinates and the
    # interface merge is key-exact. Recovery inserts Steiner points only on the
    # cap features themselves; any failure falls through to the rational stage.
    errors=String[]
    result=nothing
    reason=Ref("")
    try
        cand=_float_fill(cap)
        m=_merge_fill(cand,cap_keys,nout,caller,reason)
        if m===nothing
            push!(errors,"float fill: "*reason[])
        else
            result=_assemble_and_certify(m...,cand,coords,cap_coords,prisms,
                surface.tris,cap_off,nt,tet_limit,wall_faces,comp_tris,
                cavity_set,wall_vol,caller,reason)
            result===nothing && push!(errors,"float fill: "*reason[])
        end
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        push!(errors,"float fill: "*sprint(showerror,err))
    end

    # Stage 2 — exact-rational conforming recovery (bit-exact interface).
    if result===nothing
        reason[]=""
        try
            cand=recover_boundary_cdt(cap)
            m=_merge_fill(cand,cap_keys,nout,caller,reason)
            if m===nothing
                push!(errors,"exact recovery: "*reason[])
            else
                result=_assemble_and_certify(m...,cand,coords,cap_coords,prisms,
                    surface.tris,cap_off,nt,tet_limit,wall_faces,comp_tris,
                    cavity_set,wall_vol,caller,reason)
                result===nothing && push!(errors,"exact recovery: "*reason[])
            end
        catch err
            err isa InterruptException && rethrow()
            (err isa ArgumentError || err isa ErrorException) || rethrow()
            push!(errors,"exact recovery: "*sprint(showerror,err))
        end
    end

    result===nothing && throw(ErrorException(
        "$caller: remaining core could not be tetrahedralized and certified " *
        "(reduce hwall/nlayers or coarsen the wall). Attempts: "*join(errors," | ")))

    all_coords,tets=result
    blocks=[ElementBlock(6,prisms),ElementBlock(4,tets)]
    mesh=MixedMesh(all_coords,blocks)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return mesh
end

# Sorted-node face/edge sets of a tet mesh, used to skip recovery for features
# that are already present exactly.
@inline function _sort3i(a::Int32,b::Int32,c::Int32)
    b<c || ((b,c)=(c,b))
    a<b && return (a,b,c)
    a<c && return (b,a,c)
    return (b,c,a)
end

function _tet_face_set(m)
    s=Set{NTuple{3,Int32}}()
    @inbounds for t in axes(m.tets,2)
        v=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
        push!(s,_sort3i(v[2],v[3],v[4]))
        push!(s,_sort3i(v[1],v[3],v[4]))
        push!(s,_sort3i(v[1],v[2],v[4]))
        push!(s,_sort3i(v[1],v[2],v[3]))
    end
    return s
end

function _tet_edge_set(m)
    s=Set{NTuple{2,Int32}}()
    @inbounds for t in axes(m.tets,2)
        v=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
        for (a,b) in ((v[1],v[2]),(v[1],v[3]),(v[1],v[4]),(v[2],v[3]),(v[2],v[4]),(v[3],v[4]))
            push!(s,a<b ? (a,b) : (b,a))
        end
    end
    return s
end

# Stage-1 core engine for Delaunay-friendly caps: exact-coordinate Float64
# Delaunay (no SoS perturbation, so wall vertices stay bit-exact) followed by
# targeted recovery of whichever cap edges are missing. A cheap pre-check
# defers fundamentally non-Delaunay caps (smooth/near-cospherical surfaces lose
# almost every crease edge) to the exact-rational stage, and a hard Steiner
# budget bounds any recovery work instead of grinding.
function _float_fill(cap)
    nn=size(cap.coords,2); ntri=size(cap.tris,2)
    xs=Vector{Float64}(undef,nn); ys=similar(xs); zs=similar(xs)
    @inbounds for i in 1:nn
        xs[i]=cap.coords[1,i]; ys[i]=cap.coords[2,i]; zs[i]=cap.coords[3,i]
    end
    T=delaunay3d(xs,ys,zs; perturb=false)
    m=to_mesh3(T)
    pt(i)=(cap.coords[1,i],cap.coords[2,i],cap.coords[3,i])
    # Hardness gate: smooth/near-cospherical caps lose almost every crease edge
    # in the Delaunay tessellation, while flat caps lose only ambiguous quad
    # diagonals. Defer only in the former case.
    es=_tet_edge_set(m)
    nedge=0; miss_e=0
    @inbounds for t in 1:ntri
        for (i,j) in ((1,2),(2,3),(3,1))
            a=Int(cap.tris[i,t]); b=Int(cap.tris[j,t])
            k=a<b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))
            nedge+=1
            k in es || (miss_e+=1)
        end
    end
    3*miss_e>nedge && throw(ErrorException(
        "float fill: $miss_e/$nedge crease edges are non-Delaunay; deferring to exact recovery"))
    budget=max(64,2nn)
    for e in sort!(collect(Set{NTuple{2,Int}}(
            Int(cap.tris[i,t])<Int(cap.tris[j,t]) ? (Int(cap.tris[i,t]),Int(cap.tris[j,t])) :
            (Int(cap.tris[j,t]),Int(cap.tris[i,t]))
            for t in 1:ntri for (i,j) in ((1,2),(2,3),(3,1)))))
        (Int32(e[1]),Int32(e[2])) in es && continue
        size(m.coords,2)-nn>=budget && throw(ErrorException(
            "float fill: Steiner budget exhausted; deferring to exact recovery"))
        m=recover_segment3(m,pt(e[1]),pt(e[2]))
        es=_tet_edge_set(m)
    end
    fs=_tet_face_set(m)
    @inbounds for t in 1:ntri
        a=Int32(cap.tris[1,t]); b=Int32(cap.tris[2,t]); c=Int32(cap.tris[3,t])
        _sort3i(a,b,c) in fs && continue
        size(m.coords,2)-nn>=budget && throw(ErrorException(
            "float fill: Steiner budget exhausted; deferring to exact recovery"))
        m=recover_triangle3(m,pt(Int(a)),pt(Int(b)),pt(Int(c)))
        fs=_tet_face_set(m)
    end
    keep=falses(size(m.tets,2))
    g=_raygrid(cap)
    @inbounds for t in axes(m.tets,2)
        cx=(m.coords[1,m.tets[1,t]]+m.coords[1,m.tets[2,t]]+
            m.coords[1,m.tets[3,t]]+m.coords[1,m.tets[4,t]])/4
        cy=(m.coords[2,m.tets[1,t]]+m.coords[2,m.tets[2,t]]+
            m.coords[2,m.tets[3,t]]+m.coords[2,m.tets[4,t]])/4
        cz=(m.coords[3,m.tets[1,t]]+m.coords[3,m.tets[2,t]]+
            m.coords[3,m.tets[3,t]]+m.coords[3,m.tets[4,t]])/4
        keep[t]=_inside_grid((cx,cy,cz),g)
    end
    used=collect(Set{Int32}(v for t in findall(keep) for v in view(m.tets,:,t)))
    sort!(used; by=v->(m.coords[1,v],m.coords[2,v],m.coords[3,v]))
    nid=Dict{Int32,Int32}()
    coords=Matrix{Float64}(undef,3,length(used))
    @inbounds for (k,v) in enumerate(used)
        nid[v]=Int32(k)
        coords[1,k]=m.coords[1,v]; coords[2,k]=m.coords[2,v]; coords[3,k]=m.coords[3,v]
    end
    kept=findall(keep)
    tets=Matrix{Int32}(undef,4,length(kept))
    @inbounds for (j,t) in enumerate(kept), r in 1:4
        tets[r,j]=nid[m.tets[r,t]]
    end
    return Mesh(coords; tets=tets)
end

# Merge a candidate core onto the prism-stack numbering through bit-exact cap
# keys. Both engines preserve the wall's input coordinates exactly — the float
# stage runs unperturbed and recovery inserts only points constructed from
# those same exact coordinates — so a core node either *is* a cap node or
# becomes a new Steiner node. Every cap node must be claimed exactly once, or
# the merge is rejected (returns nothing with `reason` set).
function _merge_fill(fm, cap_keys, nout, caller, reason::Ref{String})
    nv=length(cap_keys); nf=size(fm.coords,2)
    remap=Vector{Int32}(undef,nf)
    newpts=NTuple{3,Float64}[]
    total=nout
    used=Set{Int}()
    claimed=Set{Int}()
    @inbounds for j in 1:nf
        pj=(fm.coords[1,j],fm.coords[2,j],fm.coords[3,j])
        gid=get(cap_keys,_norm_key(pj[1],pj[2],pj[3]),0)
        if gid!=0
            gid in claimed && (reason[]="two core nodes claim one cap position";
                               return nothing)
            push!(claimed,gid)
        else
            total+1<=typemax(Int32) ||
                throw(ArgumentError("$caller: combined node count exceeds Int32"))
            push!(newpts,_norm_key(pj[1],pj[2],pj[3]))
            total+=1
            gid=total
        end
        gid in used && (reason[]="collapsed node binding"; return nothing)
        push!(used,gid)
        remap[j]=Int32(gid)
    end
    length(claimed)==nv || (reason[]="core does not recover every wall node";
                            return nothing)
    return remap,newpts,total
end

# Build the combined coordinate/connectivity arrays and run the full
# certification stack. Returns (all_coords, tets) on success, nothing otherwise
# with `reason` set.
function _assemble_and_certify(remap, newpts, total, fm, stack_coords,
                               cap_coords, prisms, surface_tris, cap_off, nt,
                               max_tets, wall_faces, comp_tris, cavity_set,
                               wall_vol, caller, reason::Ref{String})
    ntet=size(fm.tets,2)
    (ntet>0 && ntet<=max_tets) ||
        (reason[]="core tetrahedron count $ntet out of contract"; return nothing)
    nnew=length(newpts)
    all_coords=Matrix{Float64}(undef,3,total)
    @inbounds for i in axes(stack_coords,2)
        all_coords[1,i]=stack_coords[1,i]
        all_coords[2,i]=stack_coords[2,i]
        all_coords[3,i]=stack_coords[3,i]
    end
    @inbounds for j in 1:nnew
        i=total-nnew+j
        p=newpts[j]
        all_coords[1,i]=p[1]; all_coords[2,i]=p[2]; all_coords[3,i]=p[3]
    end
    tets=Matrix{Int32}(undef,4,ntet)
    @inbounds for t in 1:ntet, r in 1:4
        tets[r,t]=remap[Int(fm.tets[r,t])]
    end

    # Interface tiling: the core boundary must reproduce the cap area exactly.
    bnd,maxinc=boundary_faces(tets)
    maxinc<=2 || (reason[]="non-manifold core (face incidence $maxinc)"; return nothing)
    area_bnd=0.0
    @inbounds for f in bnd
        area_bnd+=triangle_area((all_coords[1,f[1]],all_coords[2,f[1]],all_coords[3,f[1]]),
                                (all_coords[1,f[2]],all_coords[2,f[2]],all_coords[3,f[2]]),
                                (all_coords[1,f[3]],all_coords[2,f[3]],all_coords[3,f[3]]))
    end
    area_cap=0.0
    @inbounds for t in 1:nt
        area_cap+=triangle_area((cap_coords[1,surface_tris[1,t]],cap_coords[2,surface_tris[1,t]],cap_coords[3,surface_tris[1,t]]),
                                (cap_coords[1,surface_tris[2,t]],cap_coords[2,surface_tris[2,t]],cap_coords[3,surface_tris[2,t]]),
                                (cap_coords[1,surface_tris[3,t]],cap_coords[2,surface_tris[3,t]],cap_coords[3,surface_tris[3,t]]))
    end
    (area_cap>0 && abs(area_bnd-area_cap)<=1e-9*area_cap) ||
        (reason[]="core boundary area $area_bnd does not tile cap area $area_cap";
         return nothing)

    # Strictly positive exact-predicate cell volumes everywhere.
    @inbounds for p in axes(prisms,2)
        _prism_volume6(all_coords,(Int(prisms[1,p]),Int(prisms[2,p]),Int(prisms[3,p]),
                                   Int(prisms[4,p]),Int(prisms[5,p]),Int(prisms[6,p])))>0 ||
            (reason[]="prism $p is degenerate or inverted; hwall exceeds this wall's feature size";
             return nothing)
    end
    @inbounds for t in 1:ntet
        pa=(all_coords[1,tets[1,t]],all_coords[2,tets[1,t]],all_coords[3,tets[1,t]])
        pb=(all_coords[1,tets[2,t]],all_coords[2,tets[2,t]],all_coords[3,tets[2,t]])
        pc=(all_coords[1,tets[3,t]],all_coords[2,tets[3,t]],all_coords[3,tets[3,t]])
        pd=(all_coords[1,tets[4,t]],all_coords[2,tets[4,t]],all_coords[3,tets[4,t]])
        -orient3(pa,pb,pc,pd)>0 ||
            (reason[]="tetrahedron $t is degenerate or inverted"; return nothing)
    end

    # Per-wall shell identities and the global fill identity.
    nc=length(wall_faces)
    cap_vols=Vector{Float64}(undef,nc)
    for c in 1:nc
        cap_vols[c]=abs(_div_volume(all_coords,
            [(cap_off+f[1],cap_off+f[2],cap_off+f[3]) for f in wall_faces[c]]))
        pv=0.0
        for t in comp_tris[c], k in 0:(size(prisms,2)÷nt)-1
            col=k*nt+t
            pv+=_prism_volume6(all_coords,(Int(prisms[1,col]),Int(prisms[2,col]),
                Int(prisms[3,col]),Int(prisms[4,col]),Int(prisms[5,col]),Int(prisms[6,col])))
        end
        expected=c in cavity_set ? cap_vols[c]-abs(wall_vol[c]) : abs(wall_vol[c])-cap_vols[c]
        abs(pv-expected)<=1e-9*abs(wall_vol[c]) ||
            (reason[]="wall $c shell volume $pv violates its identity (expected $expected)";
             return nothing)
    end
    expected_fill=sum(c in cavity_set ? -cap_vols[c] : cap_vols[c] for c in 1:nc)
    vfill=0.0
    @inbounds for t in 1:ntet
        vfill+=tet_volume((all_coords[1,tets[1,t]],all_coords[2,tets[1,t]],all_coords[3,tets[1,t]]),
                          (all_coords[1,tets[2,t]],all_coords[2,tets[2,t]],all_coords[3,tets[2,t]]),
                          (all_coords[1,tets[3,t]],all_coords[2,tets[3,t]],all_coords[3,tets[3,t]]),
                          (all_coords[1,tets[4,t]],all_coords[2,tets[4,t]],all_coords[3,tets[4,t]]))
    end
    expected_fill>0 || (reason[]="no positive fill region"; return nothing)
    abs(vfill-expected_fill)<=1e-9*expected_fill ||
        (reason[]="filled volume $vfill violates the global identity (expected $expected_fill)";
         return nothing)

    return all_coords,tets
end

end # module
