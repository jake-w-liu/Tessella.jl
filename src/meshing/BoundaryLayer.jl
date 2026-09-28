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
                _raygrid, _inside_grid, _rb_fan_steiner, _rb_orient_facets,
                _signed_vol6
using ..RecoverCDT: recover_boundary_cdt

export mesh_boundary_layer, mesh_boundary_layer_2d, mesh_boundary_layer_filled,
       mesh_boundary_layer_fan

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
        # The cross product is already twice the area-weighted unit normal.
        for id in (i,j,k)
            normals[1,id]+=n[1]; normals[2,id]+=n[2]; normals[3,id]+=n[3]
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

The core is filled by a bounded ladder: an exact-coordinate Float64 Delaunay
path for Delaunay-friendly caps (flat/mostly-planar walls) under a hard Steiner
budget, an interior-kernel fan for star-shaped caps, then an exact-rational
conforming-recovery pass for harder caps. All
stages are held to the same independent certificates; whichever succeeds first
is returned.

The returned mesh is certified before return:

- every wall vertex survives recovery and merges into one shared node per cap
  position, so prisms and core tets conform with no cracks;
- the recovered tetrahedron boundary has exactly the prism cap's triangular
  faces and equal total area to 1e-9 relative;
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
        # The cross product is already twice the area-weighted unit normal.
        for id in (i,j,k)
            normals[1,id]+=n[1]; normals[2,id]+=n[2]; normals[3,id]+=n[3]
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

    # A fan uses only one interior Steiner point and preserves every cap face.
    # It is available only when the existing kernel test certifies a connected
    # star-shaped cap; cavities and non-star-shaped caps proceed to recovery.
    if result===nothing && nt<=tet_limit
        reason[]=""
        facets=NTuple{3,Int32}[(cap.tris[1,t],cap.tris[2,t],cap.tris[3,t]) for t in 1:nt]
        cand=_rb_fan_steiner(cap.coords[1,:],cap.coords[2,:],cap.coords[3,:],facets)
        if cand===nothing
            push!(errors,"star fan: cap has no certified interior kernel")
        else
            m=_merge_fill(cand,cap_keys,nout,caller,reason)
            if m===nothing
                push!(errors,"star fan: "*reason[])
            else
                result=_assemble_and_certify(m...,cand,coords,cap_coords,prisms,
                    surface.tris,cap_off,nt,tet_limit,wall_faces,comp_tris,
                    cavity_set,wall_vol,caller,reason)
                result===nothing && push!(errors,"star fan: "*reason[])
            end
        end
    end

    # Ear-clipping fill for non-star-shaped caps: preserves every cap face by
    # vertex identity while carving off one vertex's cone at a time.
    if result===nothing && nt<=tet_limit
        reason[]=""
        try
            cand=_vclip_fill(cap)
            cand===nothing && throw(ErrorException("no removable ear"))
            m=_merge_fill(cand,cap_keys,nout,caller,reason)
            if m===nothing
                push!(errors,"ear clip: "*reason[])
            else
                result=_assemble_and_certify(m...,cand,coords,cap_coords,prisms,
                    surface.tris,cap_off,nt,tet_limit,wall_faces,comp_tris,
                    cavity_set,wall_vol,caller,reason)
                result===nothing && push!(errors,"ear clip: "*reason[])
            end
        catch err
            err isa InterruptException && rethrow()
            (err isa ArgumentError || err isa ErrorException) || rethrow()
            push!(errors,"ear clip: "*sprint(showerror,err))
        end
    end

    # Exact-rational conforming recovery (bit-exact interface).
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

function _cap_node_ids(m,cap)
    ids=Dict(_norm_key(cap.coords[1,i],cap.coords[2,i],cap.coords[3,i])=>Int32(i)
             for i in axes(cap.coords,2))
    return [get(ids,_norm_key(m.coords[1,i],m.coords[2,i],m.coords[3,i]),Int32(0))
            for i in axes(m.coords,2)]
end

function _tet_face_set(m,cap)
    ids=_cap_node_ids(m,cap)
    s=Set{NTuple{3,Int32}}()
    @inbounds for t in axes(m.tets,2)
        v=(ids[m.tets[1,t]],ids[m.tets[2,t]],ids[m.tets[3,t]],ids[m.tets[4,t]])
        for (a,b,c) in ((v[2],v[3],v[4]),(v[1],v[3],v[4]),
                        (v[1],v[2],v[4]),(v[1],v[2],v[3]))
            min(a,b,c)>0 && push!(s,_sort3i(a,b,c))
        end
    end
    return s
end

function _tet_edge_set(m,cap)
    ids=_cap_node_ids(m,cap)
    s=Set{NTuple{2,Int32}}()
    @inbounds for t in axes(m.tets,2)
        v=(ids[m.tets[1,t]],ids[m.tets[2,t]],ids[m.tets[3,t]],ids[m.tets[4,t]])
        for (a,b) in ((v[1],v[2]),(v[1],v[3]),(v[1],v[4]),(v[2],v[3]),(v[2],v[4]),(v[3],v[4]))
            min(a,b)>0 && push!(s,a<b ? (a,b) : (b,a))
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
    es=_tet_edge_set(m,cap)
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
        es=_tet_edge_set(m,cap)
    end
    fs=_tet_face_set(m,cap)
    @inbounds for t in 1:ntri
        a=Int32(cap.tris[1,t]); b=Int32(cap.tris[2,t]); c=Int32(cap.tris[3,t])
        _sort3i(a,b,c) in fs && continue
        size(m.coords,2)-nn>=budget && throw(ErrorException(
            "float fill: Steiner budget exhausted; deferring to exact recovery"))
        m=recover_triangle3(m,pt(Int(a)),pt(Int(b)),pt(Int(c)))
        fs=_tet_face_set(m,cap)
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
    # single-assignment copies: closures over the reassigned `m` box it
    final_tets=m.tets;final_coords=m.coords
    used=collect(Set{Int32}(v for t in findall(keep) for v in view(final_tets,:,t)))
    sort!(used; by=v->(final_coords[1,v],final_coords[2,v],final_coords[3,v]))
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

# Segment (p,q) piercing the interior of triangle (a,b,c) — a float heuristic
# for the ear-clipping fill; exactness is enforced by the downstream certify.
function _seg_pierces_tri(p,q,a,b,c)
    ab,ac,ap=(b[1]-a[1],b[2]-a[2],b[3]-a[3]),(c[1]-a[1],c[2]-a[2],c[3]-a[3]),(a[1]-p[1],a[2]-p[2],a[3]-p[3])
    n=_cross3(ab,ac); d=(q[1]-p[1],q[2]-p[2],q[3]-p[3]); den=_dot3(n,d)
    abs(den)<1e-14 && return false
    t=_dot3(n,ap)/den
    (t<=1e-9 || t>=1-1e-9) && return false
    h=(p[1]+d[1]*t,p[2]+d[2]*t,p[3]+d[3]*t)
    v0,v1,v2=ab,ac,(h[1]-a[1],h[2]-a[2],h[3]-a[3])
    d00,d01,d11=_dot3(v0,v0),_dot3(v0,v1),_dot3(v1,v1)
    d20,d21=_dot3(v2,v0),_dot3(v2,v1)
    den2=d00*d11-d01*d01
    abs(den2)<1e-30 && return false
    v=(d11*d20-d01*d21)/den2; w=(d00*d21-d01*d20)/den2
    return 1-v-w>1e-9 && v>1e-9 && w>1e-9
end
function _tris_cross3(p1,p2,p3,a,b,c)
    _seg_pierces_tri(p1,p2,a,b,c) && return true
    _seg_pierces_tri(p2,p3,a,b,c) && return true
    _seg_pierces_tri(p3,p1,a,b,c) && return true
    _seg_pierces_tri(a,b,p1,p2,p3) && return true
    _seg_pierces_tri(b,c,p1,p2,p3) && return true
    _seg_pierces_tri(c,a,p1,p2,p3) && return true
    return false
end

# Strict-interior point-in-tetrahedron on the exact orient3 predicate: the
# point must lie on the opposite vertex's side of every face — the expected
# sign of each face test is derived from the tetrahedron's own orientation,
# so either winding is accepted.
function _pt_in_tet(p,a,b,c,d)
    s=orient3(a,b,c,d)
    s==0 && return false
    orient3(b,c,d,p)==-s || return false
    orient3(a,d,c,p)==-s || return false
    orient3(a,b,d,p)==-s || return false
    orient3(a,c,b,p)==-s || return false
    return true
end

# Tetrahedron interiors overlapping — vertex containment either way or a proper
# face-face crossing. Shared faces and shared boundary edges are not overlaps.
function _tets_overlap(t1,t2)
    @inbounds for q in t1
        _pt_in_tet(q,t2[1],t2[2],t2[3],t2[4]) && return true
    end
    @inbounds for q in t2
        _pt_in_tet(q,t1[1],t1[2],t1[3],t1[4]) && return true
    end
    @inbounds for f1 in ((t1[2],t1[3],t1[4]),(t1[1],t1[4],t1[3]),
                        (t1[1],t1[2],t1[4]),(t1[1],t1[3],t1[2])),
                f2 in ((t2[2],t2[3],t2[4]),(t2[1],t2[4],t2[3]),
                       (t2[1],t2[2],t2[4]),(t2[1],t2[3],t2[2]))
        _tris_cross3(f1[1],f1[2],f1[3],f2[1],f2[2],f2[3]) && return true
    end
    return false
end

# Ear-clipping conforming fill of a closed triangle surface: each level picks a
# vertex whose link polygon admits a fan triangulation whose coned-off ear
# region stays inside the domain, emits the ear tets (v, link-tri), and recurses
# on the reduced surface. Preserves every input facet by vertex identity — the
# strict contract recover_boundary_cdt cannot give — so it is the last resort
# before the exact engine when the cap is not star-shaped. Returns a Mesh with
# only tets, or nothing when no removable ear exists (a Schoenhardt-type cap,
# which the exact recovery stage then handles or rejects explicitly).
function _vclip_fill(cap::Mesh)
    res=_vclip_rec(cap,NTuple{4,NTuple{3,Float64}}[])
    res===nothing && return nothing
    C,tets=res
    out=Matrix{Int32}(undef,4,length(tets))
    @inbounds for (j,tt) in enumerate(tets), r in 1:4
        out[r,j]=tt[r]
    end
    return Mesh(C;tets=out)
end
function _vclip_rec(cap::Mesh,anc::Vector{NTuple{4,NTuple{3,Float64}}})
    nf=size(cap.tris,2); nn=size(cap.coords,2)
    facets=NTuple{3,Int32}[(cap.tris[1,t],cap.tris[2,t],cap.tris[3,t]) for t in 1:nf]
    oriented=_rb_orient_facets(facets)
    oriented===nothing && return nothing
    vt(i)=(cap.coords[1,i],cap.coords[2,i],cap.coords[3,i])
    if nn==4
        nf==4 || return nothing
        v6=_signed_vol6(vt(1),vt(2),vt(3),vt(4))
        v6==0 && return nothing
        return (cap.coords,[v6>0 ? (Int32(1),Int32(2),Int32(3),Int32(4)) :
                                   (Int32(2),Int32(1),Int32(3),Int32(4))])
    end
    nn<4 && return nothing
    g=_raygrid(cap)
    vfaces=[Int[] for _ in 1:nn]
    @inbounds for (t,(a,b,c)) in enumerate(oriented), v in (a,b,c)
        push!(vfaces[Int(v)],t)
    end
    tets=NTuple{4,Int32}[]
    for v in sortperm(length.(vfaces))           # fewest faces first, stable
        fs=vfaces[v]; k=length(fs)
        k<3 && continue
        ein=Dict{NTuple{2,Int32},Int}()
        @inbounds for t in fs
            a,b,c=oriented[t]
            for e in (a<b ? (a,b) : (b,a), b<c ? (b,c) : (c,b), a<c ? (a,c) : (c,a))
                ein[e]=get(ein,e,0)+1
            end
        end
        link=NTuple{2,Int32}[]
        for (e,cnt) in ein
            cnt==1 && push!(link,e)
            cnt>2 && (k=-1; break)
        end
        (k==-1 || length(link)!=k) && continue
        nxt=Dict{Int32,Vector{Int32}}()
        for (x,y) in link
            push!(get!(nxt,x,Int32[]),y); push!(get!(nxt,y,Int32[]),x)
        end
        all(==(2),length.(values(nxt))) || continue
        start=minimum(e->min(e[1],e[2]),link)
        cyc=Int32[start]; prev=Int32(-1); cur=start
        for _ in 1:k+2
            w=nxt[cur][1]==prev ? nxt[cur][2] : nxt[cur][1]
            w==start && break
            push!(cyc,w); prev,cur=cur,w
        end
        length(cyc)==k || continue
        cyc=cyc[2]<cyc[end] ? cyc : vcat(cyc[1],reverse(cyc[2:end]))
        # try every fan triangulation of the link polygon (apex = each vertex);
        # accept the first whose ear region is inside and whose new link
        # triangles do not properly cross any surface face outside the fan
        notfan=[t for t in 1:nf if !(t in fs)]
        ptris=nothing
        for ap in 1:k
            A=cyc[ap]; ok=true; cand=NTuple{3,Int32}[]
            for i in 2:k-1
                v2=cyc[mod1(ap+i-1,k)]; v3=cyc[mod1(ap+i,k)]
                push!(cand,(A,v2,v3))
            end
            for pt in cand
                q=((vt(v)[1]+vt(pt[1])[1]+vt(pt[2])[1]+vt(pt[3])[1])/4,
                   (vt(v)[2]+vt(pt[1])[2]+vt(pt[2])[2]+vt(pt[3])[2])/4,
                   (vt(v)[3]+vt(pt[1])[3]+vt(pt[2])[3]+vt(pt[3])[3])/4)
                _inside_grid(q,g) || (ok=false; break)
            end
            ok || continue
            # the ear region must be empty of surface material: no link
            # triangle may properly cross a remaining face and no remaining
            # vertex may lie strictly inside any ear tet (exact orient3);
            # ear-tet centroids inside R plus boundary conformance then give
            # E ⊆ R, so the residual surface bounds R ∖ E
            pv=vt(v)
            for pt in cand, t in notfan
                a,b,c=oriented[t]
                _tris_cross3(vt(pt[1]),vt(pt[2]),vt(pt[3]),
                             vt(a),vt(b),vt(c)) && (ok=false; break)
            end
            ok || continue
            nvset=Set{Int32}()
            for t in notfan, r in 1:3
                push!(nvset,oriented[t][r])
            end
            for pt in cand
                p1,p2,p3=vt(pt[1]),vt(pt[2]),vt(pt[3])
                for u in nvset
                    _pt_in_tet(vt(u),pv,p1,p2,p3) && (ok=false; break)
                end
                ok || break
            end
            ok || continue
            # the ear region must leave the remaining domain: with the link
            # triangles replacing the v-fan, every ear tet must sit OUTSIDE the
            # residual surface — otherwise later ears re-carve its space
            rfaces=NTuple{3,Int32}[oriented[t] for t in notfan]
            append!(rfaces,cand)
            rg=_raygrid(Mesh(cap.coords;tris=reshape(
                [rfaces[t][r] for t in 1:length(rfaces),r in 1:3],3,:)))
            for pt in cand
                q=((vt(v)[1]+vt(pt[1])[1]+vt(pt[2])[1]+vt(pt[3])[1])/4,
                   (vt(v)[2]+vt(pt[1])[2]+vt(pt[2])[2]+vt(pt[3])[2])/4,
                   (vt(v)[3]+vt(pt[1])[3]+vt(pt[2])[3]+vt(pt[3])[3])/4)
                _inside_grid(q,rg) && (ok=false; break)
            end
            ok || continue
            # no new ear tet may overlap a tetrahedron already emitted at a
            # shallower level — deeper residual surfaces can legitimately
            # re-enter earlier ears' space, so this is checked pairwise rather
            # than through the surface
            et=NTuple{4,NTuple{3,Float64}}[]
            for pt in cand
                p1,p2,p3=vt(pt[1]),vt(pt[2]),vt(pt[3])
                v6=_signed_vol6(pv,p1,p2,p3)
                v6==0 && (ok=false; break)
                push!(et,v6>0 ? (pv,p1,p2,p3) : (p1,pv,p2,p3))
            end
            ok || continue
            for i in eachindex(et), j in 1:i-1
                _tets_overlap(et[i],et[j]) && (ok=false; break)
            end
            ok || continue
            for e in et
                for a in anc
                    _tets_overlap(e,a) && (ok=false; break)
                end
                ok || break
            end
            ok && (ptris=cand; break)
        end
        ptris===nothing && continue
        pv=vt(v)
        n0=length(tets)
        for pt in ptris
            v6=_signed_vol6(pv,vt(pt[1]),vt(pt[2]),vt(pt[3]))
            v6==0 && break
            push!(tets, v6>0 ? (Int32(v),pt[1],pt[2],pt[3]) :
                               (pt[1],Int32(v),pt[2],pt[3]))
        end
        length(tets)-n0==length(ptris) || (resize!(tets,n0); continue)
        # residual surface: S - fan + link tris, with v removed from numbering
        keep=sort!(setdiff(1:nn,[v]))
        remap=Dict{Int32,Int32}(Int32(u)=>Int32(i) for (i,u) in enumerate(keep))
        S2=Int32[]
        @inbounds for (t,(a,b,c)) in enumerate(oriented)
            t in fs && continue
            push!(S2,remap[a],remap[b],remap[c])
        end
        for pt in ptris
            push!(S2,remap[pt[1]],remap[pt[2]],remap[pt[3]])
        end
        anc2=copy(anc)
        for pt in ptris
            p1,p2,p3=vt(pt[1]),vt(pt[2]),vt(pt[3])
            v6=_signed_vol6(pv,p1,p2,p3)
            push!(anc2,v6>0 ? (pv,p1,p2,p3) : (p1,pv,p2,p3))
        end
        sub=_vclip_rec(Mesh(cap.coords[:,keep];tris=reshape(S2,3,:)),anc2)
        sub===nothing && (resize!(tets,n0); continue)
        scoords,stets=sub
        nkeep=length(keep)
        extra=scoords[:,nkeep+1:end]
        C=hcat(cap.coords,extra)
        okall=true
        for (a,b,c,d) in stets
            A=a<=nkeep ? Int32(keep[Int(a)]) : Int32(nn+(Int(a)-nkeep))
            B=b<=nkeep ? Int32(keep[Int(b)]) : Int32(nn+(Int(b)-nkeep))
            C2=c<=nkeep ? Int32(keep[Int(c)]) : Int32(nn+(Int(c)-nkeep))
            D=d<=nkeep ? Int32(keep[Int(d)]) : Int32(nn+(Int(d)-nkeep))
            v6=_signed_vol6((C[1,A],C[2,A],C[3,A]),(C[1,B],C[2,B],C[3,B]),
                            (C[1,C2],C[2,C2],C[3,C2]),(C[1,D],C[2,D],C[3,D]))
            v6==0 && (okall=false; break)
            push!(tets, v6>0 ? (A,B,C2,D) : (B,A,C2,D))
        end
        okall || (resize!(tets,n0); continue)
        return (C,tets)
    end
    return nothing
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

    # A first-order prism's triangular cap cannot meet subdivided tet faces:
    # equal area and shared original vertices still permit hanging nodes.
    bnd,maxinc=boundary_faces(tets)
    maxinc<=2 || (reason[]="non-manifold core (face incidence $maxinc)"; return nothing)
    cap_faces=Set(_sort3i(Int32(cap_off+surface_tris[1,t]),
                         Int32(cap_off+surface_tris[2,t]),
                         Int32(cap_off+surface_tris[3,t])) for t in 1:nt)
    Set(_sort3i(f...) for f in bnd)==cap_faces ||
        (reason[]="core boundary does not match the prism cap facets"; return nothing)
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

# ═══ joined multi-region boundary layers: 3-D edge fans and corner blocks ══════

# Canonical split of a quad given cyclically (a,b,c,d): the diagonal whose
# endpoints hold the smaller pair minimum, so any two cells sharing the face
# triangulate it identically regardless of their own winding conventions.
@inline function _fan_quad_tris(a::Int32,b::Int32,c::Int32,d::Int32)
    if min(a,c)<min(b,d)
        return (_sort3i(a,b,c),_sort3i(a,c,d))
    end
    return (_sort3i(b,c,d),_sort3i(b,d,a))
end

# The eight sorted-vertex triangle faces of a type-6 prism (bottom cap, top
# cap, three side quads split on the canonical diagonal).
function _fan_prism_faces(v::NTuple{6,Int32}, sink)
    push!(sink,_sort3i(v[1],v[2],v[3]))
    push!(sink,_sort3i(v[4],v[5],v[6]))
    for (i,j) in ((1,2),(2,3),(3,1))
        f1,f2=_fan_quad_tris(v[i],v[j],v[j+3],v[i+3])
        push!(sink,f1); push!(sink,f2)
    end
    return sink
end

# Signed prism volume under the canonical side triangulation — identical in
# spirit to _prism_volume6 but the skew-side diagonal is the shared rule, so
# joined cells on both sides of an interface measure the same surface.
function _fan_prism_volume(coords, v::NTuple{6,<:Integer})
    pts=NTuple{3,Float64}[(coords[1,v[i]],coords[2,v[i]],coords[3,v[i]]) for i in 1:6]
    cen=(sum(p[1] for p in pts)/6.0,sum(p[2] for p in pts)/6.0,sum(p[3] for p in pts)/6.0)
    function vf(a,b,c)
        s=_dot3(a,_cross3(b,c))/6.0
        return orient3(a,b,c,cen)>0 ? s : -s
    end
    total=vf(pts[1],pts[2],pts[3])+vf(pts[4],pts[5],pts[6])
    for (i,j) in ((1,2),(2,3),(3,1))
        bi=pts[i]; bj=pts[j]; ti=pts[i+3]; tj=pts[j+3]
        if min(v[i],v[j+3])<min(v[j],v[i+3])
            total+=vf(bi,bj,tj)+vf(bi,tj,ti)
        else
            total+=vf(bj,bi,ti)+vf(bj,ti,tj)
        end
    end
    return total
end

# The two vertices of surface triangle t other than v.
@inline function _fan_others(tris,t,v)
    a=Int(tris[1,t]); b=Int(tris[2,t]); c=Int(tris[3,t])
    a==v && return b,c
    b==v && return a,c
    return a,b
end

function _fan_edge_map(tris,nt)
    em=Dict{NTuple{2,Int32},Vector{Int32}}()
    sizehint!(em,nt*3÷2)
    @inbounds for t in 1:nt
        a=tris[1,t]; b=tris[2,t]; c=tris[3,t]
        for (x,y) in ((a,b),(b,c),(c,a))
            key=x<y ? (x,y) : (y,x)
            push!(get!(em,key,Int32[]),Int32(t))
        end
    end
    return em
end

# Breadth-first check that one region's triangles are connected through shared
# edges; disconnected regions have no single coherent extrusion.
function _fan_check_connected(em,tris,region_of,r,trs,seen,inseen,caller)
    length(trs)<=1 && return nothing
    for t in trs; inseen[t]=true; end
    first=trs[1]; seen[first]=true; stack=Int[first]; nseen=1
    while !isempty(stack)
        t=pop!(stack)
        a=tris[1,t]; b=tris[2,t]; c=tris[3,t]
        for (x,y) in ((a,b),(b,c),(c,a))
            key=x<y ? (x,y) : (y,x)
            for t2 in em[key]
                u=Int(t2)
                if inseen[u] && !seen[u]
                    seen[u]=true; nseen+=1; push!(stack,u)
                end
            end
        end
    end
    for t in trs; inseen[t]=false; seen[t]=false; end
    nseen==length(trs) || throw(ArgumentError(
        "$caller: region $r is not edge-connected"))
    return nothing
end

# Cyclic walk through the triangle fan around one manifold vertex: returns the
# v-edges crossed between consecutive fan triangles in walk order as
# (other_endpoint, region_in, region_out) triples. A vertex whose incident
# triangles do not form a single cycle is pinched — no coherent extrusion.
function _fan_vertex_walk(tris, inc, v, region_of, em, caller)
    t0=inc[1]; n_inc=length(inc)
    o1,o2=_fan_others(tris,t0,v)
    crossings=Vector{NTuple{3,Int32}}(undef,n_inc)
    seen=falses(length(region_of)); seen[t0]=true
    cur=t0; cure=Int32(min(o1,o2)); ncross=0
    while true
        key=Int32(v)<cure ? (Int32(v),cure) : (cure,Int32(v))
        t12=em[key]
        tnext=Int(t12[1]==cur ? t12[2] : t12[1])
        ncross+=1
        crossings[ncross]=(cure,Int32(region_of[cur]),Int32(region_of[tnext]))
        tnext==t0 && break
        (seen[tnext] || ncross>=n_inc) && throw(ArgumentError(
            "$caller: vertex $v triangle fan does not form a single cycle"))
        seen[tnext]=true
        x,y=_fan_others(tris,tnext,v)
        (x==cure || y==cure) || throw(ArgumentError(
            "$caller: vertex $v triangle fan is inconsistent"))
        cure=Int32(x==cure ? y : x)
        cur=tnext
    end
    ncross==n_inc || throw(ArgumentError(
        "$caller: vertex $v is pinched by disjoint triangle fans"))
    return crossings
end

# Interior (non-endpoint) directions of one fan arc: `arctype` 1 rotates n_lo
# onto n_hi through the short dihedral sector, 2 through the reflex sector.
# Directions depend only on the arc key so distinct edges sharing the key share
# the same node columns.
function _fan_arc_dirs(n_lo,n_hi,arctype,nfan,caller,tag)
    g=_cross3(n_lo,n_hi); gl=hypot(g[1],g[2],g[3])
    co=_dot3(n_lo,n_hi)
    if gl<=64eps(Float64)
        co>0 && return fill(n_lo,nfan+1)
        ax,ay,az=abs(n_lo[1]),abs(n_lo[2]),abs(n_lo[3])
        ref=ay<=ax && ay<=az ? (0.0,1.0,0.0) :
            ax<=az ? (1.0,0.0,0.0) : (0.0,0.0,1.0)
        g=_cross3(ref,n_lo); gl=hypot(g[1],g[2],g[3])
        (isfinite(gl) && gl>0) || throw(ArgumentError(
            "$caller: cannot determine a fan axis at $tag"))
    end
    u=(g[1]/gl,g[2]/gl,g[3]/gl)
    θ=atan(gl,co)
    rays=Vector{NTuple{3,Float64}}(undef,nfan+1)
    if arctype==1
        for j in 0:nfan
            φ=j*θ/nfan
            rays[j+1]=_fan_rot(n_lo,u,φ)
        end
    else
        for j in 0:nfan
            φ=-j*(2π-θ)/nfan
            rays[j+1]=_fan_rot(n_lo,u,φ)
        end
    end
    return rays
end

@inline function _fan_rot(v,u,φ)
    c,s=cos(φ),sin(φ)
    cr=_cross3(u,v); d=_dot3(u,v)*(1-c)
    return (v[1]*c+cr[1]*s+u[1]*d,
            v[2]*c+cr[2]*s+u[2]*d,
            v[3]*c+cr[3]*s+u[3]*d)
end

@inline function _fan_fid(idN,idF,v,rlo,rhi,at,nfan,j,k)
    k==0 && return Int32(v)
    j==0 && return idN[(Int32(rlo),Int32(v),Int32(k))]
    j==nfan && return idN[(Int32(rhi),Int32(v),Int32(k))]
    return idF[(Int32(v),Int32(rlo),Int32(rhi),Int32(at),Int32(j),Int32(k))]
end

"""
    mesh_boundary_layer_fan(surface; regions, hwall, ratio, nlayers,
                            unlayered=(), cavities=(), fan_elements=3,
                            max_prisms=10_000_000, max_tets=10_000_000)
                            -> MixedMesh

Extrude the closed manifold wall(s) of `surface` into joined first-order
boundary layers: `regions` partitions the surface triangles into connected
patches (an iterable of triangle-index iterables covering every triangle
exactly once), each layered region growing its own per-vertex, area-weighted
normal columns by `nlayers` widths `hwall·ratio^(k-1)` — **into** each solid
wall's interior and **away from** walls listed in `cavities`. Along a boundary
edge between two layered regions the sector between the two extrusion wedges
is swept by `fan_elements` strip columns — one prism at the surface level and
two prisms per subsequent layer — and at a vertex where layered region arcs
meet, the residual spherical polygon is filled with corner cells sharing every
incident fan strip's end face. Regions in `unlayered` (1-based region indices)
emit no cells; their triangles stay on the core boundary and layered
neighbours end against them.

The remaining core behind the last layer is tetrahedralized by the same
certified ladder as [`mesh_boundary_layer_filled`](@ref) and the whole
assembly is certified before return: every cell has strictly positive
exact-predicate volume, the union boundary equals the input surface exactly,
the core boundary matches the stack's exposed faces node-for-node, each wall's
layer-cell volume equals `|V(wall)|−|V(cap)|` (sign reversed for cavities),
and the global fill identity `V(tets) == Σ V(solid caps) − Σ V(cavity caps)`
holds to 1e-9 relative.

Pinched or open walls, non-connected or non-contiguous regions, coincident
wedge normals on one side of an edge only, oversized `hwall`, and unmeshable
caps are explicit blockers — never defective meshes. Layer parameters are
uniform across regions; per-region grading is not expressible.
"""
function mesh_boundary_layer_fan(surface::Mesh; regions, hwall::Real,
                                 ratio::Real, nlayers::Integer, unlayered=(),
                                 cavities=(), fan_elements::Integer=3,
                                 max_prisms::Integer=10_000_000,
                                 max_tets::Integer=10_000_000)
    caller="mesh_boundary_layer_fan"
    ntris(surface)>0 || throw(ArgumentError("$caller: surface has no triangles"))
    _validate_input(surface,caller,"surface")
    applicable(iterate,regions) || throw(ArgumentError(
        "$caller: regions must be an iterable of triangle index sets"))
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=_bounded_int(nlayers,caller,"nlayers";minimum=1)
    nfan=_bounded_int(fan_elements,caller,"fan_elements";minimum=1)
    prism_limit=_bounded_int(max_prisms,caller,"max_prisms")
    tet_limit=_bounded_int(max_tets,caller,"max_tets")
    nv=nnodes(surface); nt=ntris(surface)

    comp_of,comp_tris,nc=_wall_components(surface,caller)
    cavity_set=_bounded_index_set(cavities,nc,caller,"cavity")
    nc>length(cavity_set) || throw(ArgumentError(
        "$caller: every wall is a cavity; there is no solid core to fill"))

    # Region partition: every triangle in exactly one edge-connected region.
    region_of=zeros(Int32,nt); nreg=0; region_tris=Vector{Int}[]
    for rs in regions
        nreg+=1
        s=_bounded_index_set(rs,nt,caller,"region $nreg")
        isempty(s) && throw(ArgumentError("$caller: region $nreg is empty"))
        trs=sort!(collect(s))
        for t in trs
            region_of[t]==0 || throw(ArgumentError(
                "$caller: triangle $t belongs to multiple regions"))
            region_of[t]=Int32(nreg)
        end
        push!(region_tris,trs)
        nreg<=typemax(Int32) || throw(ArgumentError(
            "$caller: region count exceeds Int32"))
    end
    nreg>0 || throw(ArgumentError("$caller: regions must not be empty"))
    missing=findfirst(==(0),region_of)
    missing===nothing || throw(ArgumentError(
        "$caller: regions do not cover triangle $missing"))
    unlayered_set=_bounded_index_set(unlayered,nreg,caller,"unlayered")
    layered=[!(r in unlayered_set) for r in 1:nreg]
    any(layered) || throw(ArgumentError(
        "$caller: every region is unlayered; there is nothing to extrude"))

    em=_fan_edge_map(surface.tris,nt)
    seen=falses(nt); inseen=falses(nt)
    for r in 1:nreg
        _fan_check_connected(em,surface.tris,region_of,r,region_tris[r],
                             seen,inseen,caller)
    end

    vinc=[Int[] for _ in 1:nv]
    @inbounds for t in 1:nt, i in 1:3
        push!(vinc[Int(surface.tris[i,t])],t)
    end
    bedges=NTuple{2,Int32}[]
    for e in sort!(collect(keys(em)))
        t12=em[e]
        region_of[t12[1]]!=region_of[t12[2]] && push!(bedges,e)
    end

    # Per-wall enclosed volume, growth direction, material-side orient3 sign.
    wall_faces=[NTuple{3,Int}[] for _ in 1:nc]
    @inbounds for f in 1:nt
        push!(wall_faces[comp_of[Int(surface.tris[1,f])]],
              (Int(surface.tris[1,f]),Int(surface.tris[2,f]),Int(surface.tris[3,f])))
    end
    wall_vol=[_div_volume(surface.coords,wall_faces[c]) for c in 1:nc]
    for c in 1:nc
        (isfinite(wall_vol[c]) && wall_vol[c]!=0) || throw(ArgumentError(
            "$caller: wall $c does not have a finite nonzero enclosed volume"))
    end
    dirs=[(c in cavity_set ? 1.0 : -1.0)*sign(wall_vol[c]) for c in 1:nc]
    sgn=[-dirs[c] for c in 1:nc]

    # Per-(region,vertex) inward normals for every region. Unlayered regions
    # emit no columns, but their normals take part in endpoint merging: a
    # layered neighbour must turn away from the unlayered wall so its side
    # face bounds the core instead of covering input triangles.
    nacc=Dict{Tuple{Int32,Int32},NTuple{3,Float64}}()
    @inbounds for r in 1:nreg
        for t in region_tris[r]
            i,j,k=Int(surface.tris[1,t]),Int(surface.tris[2,t]),Int(surface.tris[3,t])
            a=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
            b=(surface.coords[1,j],surface.coords[2,j],surface.coords[3,j])
            cc=(surface.coords[1,k],surface.coords[2,k],surface.coords[3,k])
            ab=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
            ac=(cc[1]-a[1],cc[2]-a[2],cc[3]-a[3])
            n=(ab[2]*ac[3]-ab[3]*ac[2],ab[3]*ac[1]-ab[1]*ac[3],ab[1]*ac[2]-ab[2]*ac[1])
            for id in (i,j,k)
                key=(Int32(r),Int32(id))
                old=get(nacc,key,(0.0,0.0,0.0))
                nacc[key]=(old[1]+n[1],old[2]+n[2],old[3]+n[3])
            end
        end
    end
    ndir=Dict{Tuple{Int32,Int32},NTuple{3,Float64}}()
    for ((r,v),n) in nacc
        L=hypot(n[1],n[2],n[3])
        L>0 || throw(ArgumentError("$caller: vertex $v has a zero normal in region $r"))
        s=dirs[comp_of[Int(v)]]/L
        ndir[(r,v)]=(n[1]*s,n[2]*s,n[3]*s)
    end

    # Vertex walks for every boundary vertex; region arcs must be contiguous.
    bverts=Int[]
    for e in bedges
        v,w=Int(e[1]),Int(e[2])
        (isempty(bverts) || bverts[end]!=v) && push!(bverts,v)
        (isempty(bverts) || bverts[end]!=w) && push!(bverts,w)
    end
    bverts=sort!(unique!(bverts))
    walk_of=Dict{Int,Vector{NTuple{3,Int32}}}()
    for v in bverts
        w=_fan_vertex_walk(surface.tris,vinc[v],v,region_of,em,caller)
        ntrans=0
        for x in w
            x[2]!=x[3] && (ntrans+=1)
        end
        distinct=length(Set(x[2] for x in w))
        ntrans==distinct || throw(ArgumentError(
            "$caller: regions are not contiguous around vertex $v"))
        walk_of[v]=w
    end

    # Fan-arc classification per layered-layered boundary edge. On a convex
    # (short-sector) dihedral the sector swept between the two layer normals
    # lies inside the union of the two extrusion wedges, so wedge cells would
    # double-cover the slabs: the two regions must instead share one merged
    # offset node at each endpoint, which makes their side faces on the edge
    # coincide into interior faces. Only reflex (long-sector) dihedrals leave
    # a genuine gap and receive `fan_elements` strip columns. Merging can tilt
    # a third edge's normals into a new convex or coincident pair, so the
    # classification iterates to a fixed point; merges only ever make endpoint
    # normals more coincident, bounding the number of passes.
    earc=Dict{NTuple{2,Int32},Int8}()
    eflat=Set{NTuple{2,Int32}}()
    converged=false
    for _pass in 1:4nreg+8
        merged=false
        empty!(earc); empty!(eflat)
        for e in bedges
            t1,t2=em[e]
            r1=Int(region_of[t1]); r2=Int(region_of[t2])
            L1=layered[r1]; L2=layered[r2]
            (L1 || L2) || continue
            v,w=Int(e[1]),Int(e[2])
            rlo,rhi=min(r1,r2),max(r1,r2)
            nlv=ndir[(Int32(rlo),Int32(v))]; nhv=ndir[(Int32(rhi),Int32(v))]
            nlw=ndir[(Int32(rlo),Int32(w))]; nhw=ndir[(Int32(rhi),Int32(w))]
            fv=(nlv==nhv); fw=(nlw==nhw)
            fv && fw && (push!(eflat,e); continue)
            o1,o2=_fan_edge_dihedral(surface,e,t1,t2)
            s=sgn[comp_of[v]]
            (o1==0 || o2==0) && throw(ArgumentError(
                "$caller: boundary edge $e has a degenerate dihedral angle"))
            sign(o1)!=sign(o2) && throw(ArgumentError(
                "$caller: boundary edge $e has an inconsistent dihedral angle"))
            at=sign(o1)==s ? Int8(1) : Int8(2)
            if at==1
                for u in (v,w)
                    (u==v && fv) && continue
                    (u==w && fw) && continue
                    s1=ndir[(Int32(rlo),Int32(u))]
                    s2=ndir[(Int32(rhi),Int32(u))]
                    m=(s1[1]+s2[1],s1[2]+s2[2],s1[3]+s2[3])
                    L=hypot(m[1],m[2],m[3])
                    L>0 || throw(ArgumentError(
                        "$caller: boundary edge $e forces opposing layer "*
                        "normals at vertex $u"))
                    mn=(m[1]/L,m[2]/L,m[3]/L)
                    ndir[(Int32(rlo),Int32(u))]=mn
                    ndir[(Int32(rhi),Int32(u))]=mn
                end
                merged=true
            elseif L1 && L2
                (fv || fw) && throw(ArgumentError(
                    "$caller: boundary edge $e has coincident layer normals "*
                    "at one endpoint only; the fan sector tapers to zero width"))
                earc[e]=at
            end
        end
        merged || (converged=true; break)
    end
    converged || throw(ArgumentError(
        "$caller: convex boundary sectors did not converge to merged layer normals"))

    # Per-vertex fan-arc list in walk order, and corner classification.
    arc_of_v=Dict{Int,Vector{NTuple{4,Int32}}}()
    corner_v=Int[]
    for v in bverts
        w=walk_of[v]
        all_layered=true
        for (other,rin,rout) in w
            layered[Int(rin)] || (all_layered=false)
            layered[Int(rout)] || (all_layered=false)
        end
        # fan-arc crossings (other endpoint, arctype, rin, rout) in walk order
        list=NTuple{4,Int32}[]
        for (other,rin,rout) in w
            rin==rout && continue
            e=Int32(v)<other ? (Int32(v),other) : (other,Int32(v))
            e in eflat && continue
            at=get(earc,e,Int8(0))
            at==0 && continue
            push!(list,(other,at,rin,rout))
        end
        arc_of_v[v]=list
        if all_layered && length(list)>=2
            same2=length(list)==2 &&
                sort!([list[1][3],list[1][4]])==sort!([list[2][3],list[2][4]]) &&
                list[1][2]==list[2][2]
            same2 || push!(corner_v,v)
        end
    end
    sort!(corner_v)

    # Corner polygon spec: cyclic node entries (kind,a,b,c,j):
    # kind 0 → region node (a=region); kind 1 → fan interior node
    # (a=rlo,b=rhi,c=arctype,j=canonical ray).
    polyspec_of=Dict{Int,Vector{NTuple{5,Int32}}}()
    mdir=Dict{Int,NTuple{3,Float64}}()
    for v in corner_v
        P=NTuple{5,Int32}[]
        for (other,at,rin,rout) in arc_of_v[v]
            rlo,rhi=min(Int(rin),Int(rout)),max(Int(rin),Int(rout))
            push!(P,(0,rin,0,0,0))
            if rin==rlo
                for j in 1:nfan-1
                    push!(P,(1,Int32(rlo),Int32(rhi),at,Int32(j)))
                end
            else
                for j in nfan-1:-1:1
                    push!(P,(1,Int32(rlo),Int32(rhi),at,Int32(j)))
                end
            end
        end
        polyspec_of[v]=P
        m=(0.0,0.0,0.0)
        for (other,at,rin,rout) in arc_of_v[v]
            n=ndir[(rin,Int32(v))]
            m=(m[1]+n[1],m[2]+n[2],m[3]+n[3])
        end
        L=hypot(m[1],m[2],m[3])
        L>0 || throw(ArgumentError(
            "$caller: corner vertex $v has a zero interior direction"))
        mdir[v]=(m[1]/L,m[2]/L,m[3]/L)
    end

    # Ray families for every distinct arc key.
    arcrays=Dict{NTuple{4,Int32},Vector{NTuple{3,Float64}}}()
    for v in bverts
        for (other,at,rin,rout) in arc_of_v[v]
            rlo,rhi=min(Int(rin),Int(rout)),max(Int(rin),Int(rout))
            key=(Int32(v),Int32(rlo),Int32(rhi),at)
            haskey(arcrays,key) && continue
            arcrays[key]=_fan_arc_dirs(ndir[(Int32(rlo),Int32(v))],
                ndir[(Int32(rhi),Int32(v))],Int(at),nfan,caller,
                "vertex $v regions ($rlo,$rhi)")
        end
    end

    offsets=_layer_offsets(hw,ra,nl,caller)

    # Node emission: input nodes first, then per layer — deduplicated region
    # columns (regions sharing a bitwise normal share the column), fan
    # interior rays, corner centers.
    points=Vector{NTuple{3,Float64}}(undef,nv)
    @inbounds for i in 1:nv
        points[i]=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
    end
    idN=Dict{NTuple{3,Int32},Int32}()
    idF=Dict{NTuple{6,Int32},Int32}()
    idC=Dict{NTuple{2,Int32},Int32}()
    linc=[Int[] for _ in 1:nv]
    @inbounds for ((r,v),n) in ndir
        layered[Int(r)] && push!(linc[Int(v)],Int(r))
    end
    for v in 1:nv
        sort!(linc[v]); unique!(linc[v])
    end
    dedup=Dict{Tuple{Int32,NTuple{3,Float64}},Int32}()
    for k in 1:nl
        d=offsets[k]
        empty!(dedup)
        for v in 1:nv
            pv=points[v]
            for r in linc[v]
                n=ndir[(Int32(r),Int32(v))]
                id=get(dedup,(Int32(v),n),Int32(0))
                if id==0
                    push!(points,(pv[1]+d*n[1],pv[2]+d*n[2],pv[3]+d*n[3]))
                    id=Int32(length(points)); dedup[(Int32(v),n)]=id
                end
                idN[(Int32(r),Int32(v),Int32(k))]=id
            end
        end
        for v in bverts
            pv=points[v]
            for (other,at,rin,rout) in arc_of_v[v]
                rlo,rhi=min(Int(rin),Int(rout)),max(Int(rin),Int(rout))
                rays=arcrays[(Int32(v),Int32(rlo),Int32(rhi),at)]
                for j in 1:nfan-1
                    key=(Int32(v),Int32(rlo),Int32(rhi),at,Int32(j),Int32(k))
                    haskey(idF,key) && continue
                    rj=rays[j+1]
                    push!(points,(pv[1]+d*rj[1],pv[2]+d*rj[2],pv[3]+d*rj[3]))
                    idF[key]=Int32(length(points))
                end
            end
        end
        for v in corner_v
            m=mdir[v]
            pv=points[v]
            push!(points,(pv[1]+d*m[1],pv[2]+d*m[2],pv[3]+d*m[3]))
            idC[(Int32(v),Int32(k))]=Int32(length(points))
        end
    end
    nout=length(points)
    nout<=typemax(Int32) || throw(ArgumentError("$caller: node count exceeds Int32"))
    @inbounds for p in points
        all(isfinite,p) || throw(ArgumentError(
            "$caller: extruded node is non-finite"))
    end

    # Cell counts before allocation.
    ntL=0
    for r in 1:nreg
        layered[r] && (ntL+=length(region_tris[r]))
    end
    nLL=length(earc)
    ncorn=length(corner_v)
    corner_tris=0; corner_prisms=0
    for v in corner_v
        Lp=length(polyspec_of[v])
        corner_tris=_checked_add(corner_tris,Lp,caller,"corner triangle count")
        corner_prisms=_checked_add(corner_prisms,Lp*(nl-1),caller,"corner prism count")
    end
    npr=_checked_add(_checked_mul(nl,ntL,caller,"region prism count"),
        _checked_add(_checked_mul(nLL,_checked_mul(nfan,2nl-1,caller,"fan cell count"),caller,"fan prism count"),
            corner_prisms,caller,"prism count"),caller,"prism count")
    npr<=prism_limit || throw(ArgumentError(
        "$caller: $npr prisms exceed max_prisms=$prism_limit"))
    npr<=typemax(Int32) || throw(ArgumentError("$caller: prism count exceeds Int32"))
    corner_tris=_checked_mul(corner_tris,1,caller,"corner tet count")
    corner_tris<=tet_limit || throw(ArgumentError(
        "$caller: $corner_tris stack tets exceed max_tets=$tet_limit"))

    prisms=Matrix{Int32}(undef,6,npr)
    pw=Vector{Int32}(undef,npr)
    stets=Matrix{Int32}(undef,4,corner_tris)
    tw=Vector{Int32}(undef,corner_tris)
    pc=0; tc=0
    @inbounds for k in 0:nl-1, r in 1:nreg
        layered[r] || continue
        for t in region_tris[r]
            i,j,k2=Int(surface.tris[1,t]),Int(surface.tris[2,t]),Int(surface.tris[3,t])
            b1=k==0 ? Int32(i) : idN[(Int32(r),Int32(i),Int32(k))]
            b2=k==0 ? Int32(j) : idN[(Int32(r),Int32(j),Int32(k))]
            b3=k==0 ? Int32(k2) : idN[(Int32(r),Int32(k2),Int32(k))]
            t1=idN[(Int32(r),Int32(i),Int32(k+1))]
            t2=idN[(Int32(r),Int32(j),Int32(k+1))]
            t3=idN[(Int32(r),Int32(k2),Int32(k+1))]
            pc+=1; prisms[:,pc].=(b1,b2,b3,t1,t2,t3); pw[pc]=Int32(comp_of[i])
        end
    end
    for e in sort!(collect(keys(earc)))
        v,w=Int(e[1]),Int(e[2])
        t1,t2=em[e]
        r1,r2=Int(region_of[t1]),Int(region_of[t2])
        rlo,rhi=min(r1,r2),max(r1,r2)
        at=earc[e]
        for k in 1:nl, j in 0:nfan-1
            a=_fan_fid(idN,idF,v,rlo,rhi,at,nfan,j,k-1)
            b=_fan_fid(idN,idF,w,rlo,rhi,at,nfan,j,k-1)
            c=_fan_fid(idN,idF,w,rlo,rhi,at,nfan,j+1,k-1)
            d=_fan_fid(idN,idF,v,rlo,rhi,at,nfan,j+1,k-1)
            A=_fan_fid(idN,idF,v,rlo,rhi,at,nfan,j,k)
            B=_fan_fid(idN,idF,w,rlo,rhi,at,nfan,j,k)
            C=_fan_fid(idN,idF,w,rlo,rhi,at,nfan,j+1,k)
            D=_fan_fid(idN,idF,v,rlo,rhi,at,nfan,j+1,k)
            if k==1
                pc+=1; prisms[:,pc].=(a,A,D,b,B,C); pw[pc]=Int32(comp_of[v])
            else
                if min(a,c)<min(b,d)
                    pc+=1; prisms[:,pc].=(a,b,c,A,B,C); pw[pc]=Int32(comp_of[v])
                    pc+=1; prisms[:,pc].=(a,c,d,A,C,D); pw[pc]=Int32(comp_of[v])
                else
                    pc+=1; prisms[:,pc].=(b,c,d,B,C,D); pw[pc]=Int32(comp_of[v])
                    pc+=1; prisms[:,pc].=(b,d,a,B,D,A); pw[pc]=Int32(comp_of[v])
                end
            end
        end
    end
    for v in corner_v
        P=polyspec_of[v]; m=length(P)
        for k in 1:nl
            pk=Vector{Int32}(undef,m)
            for mi in 1:m
                e=P[mi]
                pk[mi]=e[1]==0 ? idN[(e[2],Int32(v),Int32(k))] :
                    idF[(Int32(v),e[2],e[3],e[4],e[5],Int32(k))]
            end
            if k==1
                pv=points[v]; c1=idC[(Int32(v),Int32(1))]
                cp=points[Int(c1)]
                for mi in 1:m
                    m2=mi==m ? 1 : mi+1
                    pa=points[Int(pk[mi])]; pb=points[Int(pk[m2])]
                    tc+=1
                    if _signed_vol6(pv,pa,pb,cp)>=0
                        stets[:,tc].=(Int32(v),pk[mi],pk[m2],c1)
                    else
                        stets[:,tc].=(Int32(v),pk[m2],pk[mi],c1)
                    end
                    tw[tc]=Int32(comp_of[v])
                end
            else
                pk1=Vector{Int32}(undef,m)
                for mi in 1:m
                    e=P[mi]
                    pk1[mi]=e[1]==0 ? idN[(e[2],Int32(v),Int32(k-1))] :
                        idF[(Int32(v),e[2],e[3],e[4],e[5],Int32(k-1))]
                end
                c0=idC[(Int32(v),Int32(k-1))]; c1=idC[(Int32(v),Int32(k))]
                for mi in 1:m
                    m2=mi==m ? 1 : mi+1
                    pc+=1
                    prisms[:,pc].=(c0,pk1[mi],pk1[m2],c1,pk[mi],pk[m2])
                    pw[pc]=Int32(comp_of[v])
                end
            end
        end
    end
    pc==npr || throw(ErrorException("$caller: prism count mismatch"))
    tc==corner_tris || throw(ErrorException("$caller: corner tet count mismatch"))

    # Derive the core boundary: every stack face with incidence one that is not
    # an input triangle, plus every unlayered input triangle.
    inc=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    facebuf=NTuple{3,Int32}[]
    @inbounds for p in 1:npr
        v6=Tuple(prisms[:,p])
        empty!(facebuf)
        _fan_prism_faces(v6,facebuf)
        for f in facebuf
            cnt,_=get(inc,f,(0,Int32(0)))
            inc[f]=(cnt+1,pw[p])
        end
    end
    @inbounds for t in 1:corner_tris
        v4=(stets[1,t],stets[2,t],stets[3,t],stets[4,t])
        for f in (_sort3i(v4[2],v4[3],v4[4]),_sort3i(v4[1],v4[3],v4[4]),
                  _sort3i(v4[1],v4[2],v4[4]),_sort3i(v4[1],v4[2],v4[3]))
            cnt,_=get(inc,f,(0,Int32(0)))
            inc[f]=(cnt+1,tw[t])
        end
    end
    inputset=Set{NTuple{3,Int32}}()
    @inbounds for t in 1:nt
        push!(inputset,_sort3i(surface.tris[1,t],surface.tris[2,t],surface.tris[3,t]))
    end
    cap_faces=NTuple{3,Int32}[]
    cap_walls=Int32[]
    for (f,(cnt,w)) in inc
        cnt==1 && !(f in inputset) && (push!(cap_faces,f); push!(cap_walls,w))
    end
    @inbounds for t in 1:nt
        layered[Int(region_of[t])] && continue
        f=(surface.tris[1,t],surface.tris[2,t],surface.tris[3,t])
        push!(cap_faces,f); push!(cap_walls,Int32(comp_of[Int(f[1])]))
    end
    perm=sortperm(cap_faces)
    cap_faces=cap_faces[perm]; cap_walls=cap_walls[perm]
    capnodes=Int32[]
    seen_node=falses(nout)
    for f in cap_faces, i in f
        seen_node[Int(i)] || (push!(capnodes,i); seen_node[Int(i)]=true)
    end
    sort!(capnodes)
    caplocal=Dict{Int32,Int32}(g=>Int32(i) for (i,g) in enumerate(capnodes))
    cap_coords=Matrix{Float64}(undef,3,length(capnodes))
    cap_keys=Dict{NTuple{3,Float64},Int}()
    @inbounds for (i,g) in enumerate(capnodes)
        p=points[Int(g)]
        cap_coords[1,i]=p[1]; cap_coords[2,i]=p[2]; cap_coords[3,i]=p[3]
        key=_norm_key(p[1],p[2],p[3])
        previous=get(cap_keys,key,0)
        previous==0 || throw(ArgumentError(
            "$caller: cap nodes $(previous) and $g coincide"))
        cap_keys[key]=Int(g)
    end
    cap_tris=Matrix{Int32}(undef,3,length(cap_faces))
    @inbounds for (i,f) in enumerate(cap_faces)
        cap_tris[1,i]=caplocal[f[1]]; cap_tris[2,i]=caplocal[f[2]]; cap_tris[3,i]=caplocal[f[3]]
    end
    cap=Mesh(cap_coords;tris=cap_tris)

    # Same bounded fill ladder as mesh_boundary_layer_filled: exact-coordinate
    # Delaunay, interior-kernel fan, ear clipping, exact-rational conforming
    # recovery.
    errors=String[]
    result=nothing
    reason=Ref("")
    try
        cand=_float_fill(cap)
        m=_merge_fill(cand,cap_keys,nout,caller,reason)
        if m===nothing
            push!(errors,"float fill: "*reason[])
        else
            result=_fan_assemble_certify(m...,cand,points,prisms,pw,stets,tw,
                cap_faces,cap_walls,surface.tris,wall_faces,comp_tris,
                cavity_set,wall_vol,dirs,tet_limit,caller,reason)
            result===nothing && push!(errors,"float fill: "*reason[])
        end
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        push!(errors,"float fill: "*sprint(showerror,err))
    end
    if result===nothing && size(cap.tris,2)<=tet_limit
        reason[]=""
        ncf=size(cap.tris,2)
        facets=NTuple{3,Int32}[(cap.tris[1,t],cap.tris[2,t],cap.tris[3,t]) for t in 1:ncf]
        cand=_rb_fan_steiner(cap.coords[1,:],cap.coords[2,:],cap.coords[3,:],facets)
        if cand===nothing
            push!(errors,"star fan: cap has no certified interior kernel")
        else
            m=_merge_fill(cand,cap_keys,nout,caller,reason)
            if m===nothing
                push!(errors,"star fan: "*reason[])
            else
                result=_fan_assemble_certify(m...,cand,points,prisms,pw,stets,tw,
                    cap_faces,cap_walls,surface.tris,wall_faces,comp_tris,
                    cavity_set,wall_vol,dirs,tet_limit,caller,reason)
                result===nothing && push!(errors,"star fan: "*reason[])
            end
        end
    end
    if result===nothing && size(cap.tris,2)<=tet_limit
        reason[]=""
        try
            cand=_vclip_fill(cap)
            cand===nothing && throw(ErrorException("no removable ear"))
            m=_merge_fill(cand,cap_keys,nout,caller,reason)
            if m===nothing
                push!(errors,"ear clip: "*reason[])
            else
                result=_fan_assemble_certify(m...,cand,points,prisms,pw,stets,tw,
                    cap_faces,cap_walls,surface.tris,wall_faces,comp_tris,
                    cavity_set,wall_vol,dirs,tet_limit,caller,reason)
                result===nothing && push!(errors,"ear clip: "*reason[])
            end
        catch err
            err isa InterruptException && rethrow()
            (err isa ArgumentError || err isa ErrorException) || rethrow()
            push!(errors,"ear clip: "*sprint(showerror,err))
        end
    end
    if result===nothing
        reason[]=""
        try
            cand=recover_boundary_cdt(cap)
            m=_merge_fill(cand,cap_keys,nout,caller,reason)
            if m===nothing
                push!(errors,"exact recovery: "*reason[])
            else
                result=_fan_assemble_certify(m...,cand,points,prisms,pw,stets,tw,
                    cap_faces,cap_walls,surface.tris,wall_faces,comp_tris,
                    cavity_set,wall_vol,dirs,tet_limit,caller,reason)
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

    all_coords,all_tets=result
    blocks=ElementBlock[]
    npr>0 && push!(blocks,ElementBlock(6,prisms))
    push!(blocks,ElementBlock(4,all_tets))
    mesh=MixedMesh(all_coords,blocks)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return mesh
end

# Signed orient3-based dihedral votes for a boundary edge's two incident
# triangles: the opposite vertex of each triangle tested against the other
# triangle's plane. Both must agree for the sector to classify.
function _fan_edge_dihedral(surface,e,t1,t2)
    p(i)=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
    tris=surface.tris
    a=Int(e[1]); b=Int(e[2])
    x,y,z=Int(tris[1,t2]),Int(tris[2,t2]),Int(tris[3,t2])
    apex1= x!=a && x!=b ? x : y!=a && y!=b ? y : z
    x,y,z=Int(tris[1,t1]),Int(tris[2,t1]),Int(tris[3,t1])
    apex2= x!=a && x!=b ? x : y!=a && y!=b ? y : z
    o1=orient3(p(Int(tris[1,t1])),p(Int(tris[2,t1])),p(Int(tris[3,t1])),p(apex1))
    o2=orient3(p(Int(tris[1,t2])),p(Int(tris[2,t2])),p(Int(tris[3,t2])),p(apex2))
    return o1,o2
end

# Certification for the joined-stack assembly: builds the merged coordinate
# array, remaps the core candidate, then checks (a) the union of stack cells
# and core tets is a manifold whose only boundary faces are exactly the input
# triangles, (b) every prism and tet has strictly positive exact-predicate
# volume, (c) each wall's layer-cell volume equals |V(wall)|−|V(cap)| (sign
# reversed for cavities) and the global fill identity holds, all to 1e-9
# relative. Returns (all_coords, all_tets) or nothing with `reason` set.
function _fan_assemble_certify(remap,newpts,total,fm,points,prisms,pw,stets,tw,
                               cap_faces,cap_walls,surface_tris,wall_faces,
                               comp_tris,cavity_set,wall_vol,dirs,max_tets,
                               caller,reason::Ref{String})
    ncore=size(fm.tets,2); nstack=size(stets,2)
    ntet=ncore+nstack
    (ntet>0 && ntet<=max_tets) ||
        (reason[]="core tetrahedron count $ntet out of contract"; return nothing)
    all_coords=Matrix{Float64}(undef,3,total)
    @inbounds for i in eachindex(points)
        all_coords[1,i]=points[i][1]; all_coords[2,i]=points[i][2]; all_coords[3,i]=points[i][3]
    end
    nnew=length(newpts)
    @inbounds for j in 1:nnew
        i=total-nnew+j; p=newpts[j]
        all_coords[1,i]=p[1]; all_coords[2,i]=p[2]; all_coords[3,i]=p[3]
    end
    all_tets=Matrix{Int32}(undef,4,ntet)
    @inbounds for t in 1:nstack, r in 1:4
        all_tets[r,t]=stets[r,t]
    end
    @inbounds for t in 1:ncore, r in 1:4
        all_tets[r,nstack+t]=remap[Int(fm.tets[r,t])]
    end

    # Global watertightness: incidence over stack cells plus core tets; the
    # only count-1 faces allowed are the input triangles themselves.
    inc=Dict{NTuple{3,Int32},Int}()
    facebuf=NTuple{3,Int32}[]
    @inbounds for p in axes(prisms,2)
        v6=(prisms[1,p],prisms[2,p],prisms[3,p],prisms[4,p],prisms[5,p],prisms[6,p])
        empty!(facebuf)
        _fan_prism_faces(v6,facebuf)
        for f in facebuf
            inc[f]=get(inc,f,0)+1
        end
    end
    @inbounds for t in axes(all_tets,2)
        a,b,c,d=all_tets[1,t],all_tets[2,t],all_tets[3,t],all_tets[4,t]
        for f in (_sort3i(b,c,d),_sort3i(a,c,d),_sort3i(a,b,d),_sort3i(a,b,c))
            inc[f]=get(inc,f,0)+1
        end
    end
    inputset=Set{NTuple{3,Int32}}()
    @inbounds for t in axes(surface_tris,2)
        push!(inputset,_sort3i(surface_tris[1,t],surface_tris[2,t],surface_tris[3,t]))
    end
    nbnd=0; badfaces=NTuple{3,Int32}[]; overfaces=NTuple{3,Int32}[]
    for (f,cnt) in inc
        if cnt>2
            push!(overfaces,f)
        elseif cnt==1
            nbnd+=1
            f in inputset || push!(badfaces,f)
        end
    end
    isempty(overfaces) ||
        (reason[]="non-manifold union (face $(minimum(overfaces)) has incidence >2)";
         return nothing)
    isempty(badfaces) ||
        (reason[]="union boundary face $(minimum(badfaces)) is not part of the input surface";
         return nothing)
    nbnd==length(inputset) ||
        (reason[]="union boundary has $nbnd faces, expected $(length(inputset)) input triangles";
         return nothing)

    # Strictly positive exact-predicate cell volumes; layer volume per wall.
    nc=length(wall_faces)
    vlayer=zeros(Float64,nc)
    @inbounds for p in axes(prisms,2)
        v=_fan_prism_volume(all_coords,(Int(prisms[1,p]),Int(prisms[2,p]),
            Int(prisms[3,p]),Int(prisms[4,p]),Int(prisms[5,p]),Int(prisms[6,p])))
        v>0 || (reason[]="prism $p is degenerate or inverted; hwall exceeds this wall's feature size";
                return nothing)
        vlayer[Int(pw[p])]+=v
    end
    @inbounds for t in 1:nstack
        pa=(all_coords[1,all_tets[1,t]],all_coords[2,all_tets[1,t]],all_coords[3,all_tets[1,t]])
        pb=(all_coords[1,all_tets[2,t]],all_coords[2,all_tets[2,t]],all_coords[3,all_tets[2,t]])
        pc2=(all_coords[1,all_tets[3,t]],all_coords[2,all_tets[3,t]],all_coords[3,all_tets[3,t]])
        pd=(all_coords[1,all_tets[4,t]],all_coords[2,all_tets[4,t]],all_coords[3,all_tets[4,t]])
        -orient3(pa,pb,pc2,pd)>0 ||
            (reason[]="corner tetrahedron $t is degenerate or inverted"; return nothing)
        vlayer[Int(tw[t])]+=abs(_signed_vol6(pa,pb,pc2,pd))/6.0
    end
    vfill=0.0
    @inbounds for t in nstack+1:ntet
        pa=(all_coords[1,all_tets[1,t]],all_coords[2,all_tets[1,t]],all_coords[3,all_tets[1,t]])
        pb=(all_coords[1,all_tets[2,t]],all_coords[2,all_tets[2,t]],all_coords[3,all_tets[2,t]])
        pc2=(all_coords[1,all_tets[3,t]],all_coords[2,all_tets[3,t]],all_coords[3,all_tets[3,t]])
        pd=(all_coords[1,all_tets[4,t]],all_coords[2,all_tets[4,t]],all_coords[3,all_tets[4,t]])
        -orient3(pa,pb,pc2,pd)>0 ||
            (reason[]="tetrahedron $t is degenerate or inverted"; return nothing)
        vfill+=tet_volume(pa,pb,pc2,pd)
    end

    # Per-wall cap volumes: each cap face is oriented so the owning core tet's
    # apex lies on its positive orient3 side — consistently away from the
    # core for solids and into the void for cavities.
    apex=Dict{NTuple{3,Int32},Int32}()
    @inbounds for t in nstack+1:ntet
        a,b,c,d=all_tets[1,t],all_tets[2,t],all_tets[3,t],all_tets[4,t]
        apex[_sort3i(b,c,d)]=a; apex[_sort3i(a,c,d)]=b
        apex[_sort3i(a,b,d)]=c; apex[_sort3i(a,b,c)]=d
    end
    cap_signed=zeros(Float64,nc)
    @inbounds for (i,f) in enumerate(cap_faces)
        skey=_sort3i(f[1],f[2],f[3])
        ap=get(apex,skey,Int32(0))
        ap==0 && (reason[]="cap face $skey has no owning core tetrahedron";
                  return nothing)
        pa=(all_coords[1,f[1]],all_coords[2,f[1]],all_coords[3,f[1]])
        pb=(all_coords[1,f[2]],all_coords[2,f[2]],all_coords[3,f[2]])
        pc2=(all_coords[1,f[3]],all_coords[2,f[3]],all_coords[3,f[3]])
        pd=(all_coords[1,ap],all_coords[2,ap],all_coords[3,ap])
        q=_dot3(pa,_cross3(pb,pc2))/6.0
        cap_signed[Int(cap_walls[i])]+= orient3(pa,pb,pc2,pd)>0 ? q : -q
    end
    @inbounds for c in 1:nc
        capv=abs(cap_signed[c])
        expected=c in cavity_set ? capv-abs(wall_vol[c]) : abs(wall_vol[c])-capv
        abs(vlayer[c]-expected)<=1e-9*abs(wall_vol[c]) ||
            (reason[]="wall $c layer volume $(vlayer[c]) violates its identity (expected $expected)";
             return nothing)
    end
    expected_fill=sum(c in cavity_set ? -abs(cap_signed[c]) : abs(cap_signed[c]) for c in 1:nc)
    expected_fill>0 || (reason[]="no positive fill region"; return nothing)
    abs(vfill-expected_fill)<=1e-9*expected_fill ||
        (reason[]="filled volume $vfill violates the global identity (expected $expected_fill)";
         return nothing)
    return all_coords,all_tets
end

end # module
