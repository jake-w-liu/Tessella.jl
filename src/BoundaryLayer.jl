"""
    BoundaryLayer

First-order boundary-layer element topology: prismatic extrusion of a
triangle surface along area-weighted vertex normals (type-6 prisms), and
planar polyline extrusion to type-3 quads with optional convex-corner fans
(type-2 triangles in the first layer, type-3 quads in subsequent layers).
This is the element-topology counterpart of [`BoundaryLayerField`](@ref).
"""
module BoundaryLayer

using ..MeshTypes: Mesh, nnodes, nsegs, ntris, triangle_area, tet_volume
using ..Elements: ElementBlock, MixedMesh, validate
using ..Predicates: orient3

export mesh_boundary_layer, mesh_boundary_layer_2d

function _finite(value, caller, name)
    v=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    isfinite(v) || throw(ArgumentError("$caller: $name must be finite"))
    return v
end

function _layer_offsets(hw, ra, nl, caller)
    offsets=Vector{Float64}(undef,nl)
    @inbounds for k in 1:nl
        offsets[k]=hw*(ra^k-1)/(ra-1)
        isfinite(offsets[k]) || throw(ArgumentError("$caller: layer offset overflowed"))
        offsets[k]>0 || throw(ArgumentError("$caller: layer offset must be positive"))
    end
    return offsets
end

function mesh_boundary_layer(surface::Mesh; hwall::Real, ratio::Real, nlayers::Integer,
                             max_prisms::Integer=10_000_000)
    caller="mesh_boundary_layer"
    ntris(surface)>0 || throw(ArgumentError("$caller: surface has no triangles"))
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=Int(nlayers); nl>=1 || throw(ArgumentError("$caller: nlayers must be ≥ 1"))
    npr=nl*ntris(surface)
    npr<=max_prisms || throw(ArgumentError("$caller: $npr prisms exceed max_prisms=$max_prisms"))
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

    nout=nv*(nl+1)
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
    L=hypot(tx,ty)
    L>0 || throw(ArgumentError("$caller: segment $seg has zero length"))
    return (-ty/L, tx/L)
end

function _polyline_vertices(curve::Mesh, caller)
    n=nsegs(curve)
    n>0 || throw(ArgumentError("$caller: curve has no segments"))
    nv=nnodes(curve)
    adj=Dict{Int,Vector{Int}}()
    @inbounds for s in 1:n
        a=Int(curve.segs[1,s]); b=Int(curve.segs[2,s])
        (1<=a<=nv && 1<=b<=nv) || throw(ArgumentError("$caller: segment $s is out of range"))
        a==b && throw(ArgumentError("$caller: segment $s is degenerate"))
        push!(get!(()->Int[], adj, a), b)
        push!(get!(()->Int[], adj, b), a)
    end
    for (v,nbrs) in adj
        unique!(sort!(nbrs))
        1<=length(nbrs)<=2 || throw(ArgumentError(
            "$caller: vertex $v is not on a single polyline"))
    end
    endpoints=sort!(collect(v for (v,nbrs) in adj if length(nbrs)==1))
    if length(endpoints)==2
        closed=false
        start=endpoints[1]
    elseif isempty(endpoints) && length(adj)==n
        closed=true
        start=minimum(keys(adj))
    else
        throw(ArgumentError("$caller: segments are not a single polyline"))
    end
    verts=Int[start]
    prev=0
    while true
        nbrs=adj[verts[end]]
        nxt=if prev==0
            nbrs[1]
        else
            nbrs[1]==prev ? nbrs[end] : nbrs[1]
        end
        if closed && nxt==start
            break
        end
        nxt in verts && throw(ArgumentError("$caller: segments are not a single polyline"))
        push!(verts,nxt)
        prev=verts[end-1]
        if !closed && length(adj[nxt])==1
            break
        end
        length(verts)>n+1 && throw(ArgumentError("$caller: segments are not a single polyline"))
    end
    expected=closed ? n : n+1
    length(verts)==expected || throw(ArgumentError("$caller: segments are not a single polyline"))
    return verts, closed
end

function _signed_quad(ax,ay,bx,by,cx,cy,dx,dy)
    return ax*by-ay*bx + bx*cy-by*cx + cx*dy-cy*dx + dx*ay-dy*ax
end

function _id_layer(id_of, v, layer, ray, fan)
    return layer==0 ? id_of[(v,0,0)] : id_of[(v,layer, fan ? ray : 0)]
end

"""
    mesh_boundary_layer_2d(curve; hwall, ratio, nlayers, fans=(), fan_elements=5)

Extrude a constant-`z` polyline along its left-normals (CCW of the directed
chain) into first-order type-3 quadrangles. `fans` lists convex interior
vertices that receive `fan_elements` first-layer triangles and matching
ring quadrangles instead of a single averaged-normal column.
"""
function mesh_boundary_layer_2d(curve::Mesh; hwall::Real, ratio::Real, nlayers::Integer,
                                fans=(), fan_elements::Integer=5,
                                max_cells::Integer=10_000_000)
    caller="mesh_boundary_layer_2d"
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=Int(nlayers); nl>=1 || throw(ArgumentError("$caller: nlayers must be ≥ 1"))
    nfan=Int(fan_elements)
    nv=nnodes(curve)
    nv>=2 || throw(ArgumentError("$caller: curve has too few nodes"))
    z0=curve.coords[3,1]
    @inbounds for i in 1:nv
        abs(curve.coords[3,i]-z0)<=1e-12 || throw(ArgumentError(
            "$caller: curve is not constant-z; general-plane extrusion is not implemented"))
    end
    verts, closed=_polyline_vertices(curve,caller)
    nchain=length(verts)
    nchain_segs=closed ? nchain : nchain-1
    fan_set=Set{Int}()
    for raw in fans
        v=Int(raw)
        1<=v<=nv || throw(ArgumentError("$caller: fan vertex $v is out of range"))
        v in fan_set && throw(ArgumentError("$caller: duplicate fan vertex $v"))
        push!(fan_set,v)
    end
    if !isempty(fan_set)
        nfan>=2 || throw(ArgumentError("$caller: fan_elements must be ≥ 2"))
    end
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
    ntris_out=n_fans*nfan
    nquads=nchain_segs*nl + n_fans*nfan*max(nl-1,0)
    ncells=ntris_out+nquads
    ncells<=max_cells || throw(ArgumentError("$caller: $ncells cells exceed max_cells=$max_cells"))
    ncells<=typemax(Int32) || throw(ArgumentError("$caller: cell count exceeds Int32"))
    ncells>0 || throw(ArgumentError("$caller: no cells to emit"))

    points=Vector{NTuple{3,Float64}}(undef,nv)
    @inbounds for i in 1:nv
        points[i]=(curve.coords[1,i],curve.coords[2,i],z0)
    end
    seg_left=Vector{NTuple{2,Float64}}(undef,nchain_segs)
    @inbounds for i in 1:nchain_segs
        a=verts[i]; b=verts[i==nchain ? 1 : i+1]
        seg_left[i]=_left_unit(points[a][1],points[a][2],points[b][1],points[b][2],caller,i)
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
        px,py,_=points[v]
        if v in fan_set
            nin=n_in[idx]; nout=n_out[idx]
            θ=atan(nin[1]*nout[2]-nin[2]*nout[1], nin[1]*nout[1]+nin[2]*nout[2])
            for ray in 0:nfan
                α=ray/nfan*θ
                c,s=cos(α),sin(α)
                dx=nin[1]*c-nin[2]*s
                dy=nin[1]*s+nin[2]*c
                push!(points,(px+offsets[k]*dx, py+offsets[k]*dy, z0))
                all(isfinite, points[end]) || throw(ArgumentError("$caller: extruded node is non-finite"))
                id_of[(v,k,ray)]=length(points)
            end
        else
            nx,ny=n_avg[idx]
            push!(points,(px+offsets[k]*nx, py+offsets[k]*ny, z0))
            all(isfinite, points[end]) || throw(ArgumentError("$caller: extruded node is non-finite"))
            id_of[(v,k,0)]=length(points)
        end
    end
    length(points)<=typemax(Int32) || throw(ArgumentError("$caller: node count exceeds Int32"))

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
            p1=points[b1]; p2=points[b2]; p3=points[t2]; p4=points[t1]
            _signed_quad(p1[1],p1[2],p2[1],p2[2],p3[1],p3[2],p4[1],p4[2])>0 ||
                throw(ArgumentError("$caller: inverted quadrangle on segment $i layer $k"))
            qcursor+=1
            quads[:,qcursor].=(Int32(b1),Int32(b2),Int32(t2),Int32(t1))
        end
    end
    @inbounds for (idx,v) in enumerate(verts)
        v in fan_set || continue
        origin=id_of[(v,0,0)]
        for ray in 0:nfan-1
            i1=id_of[(v,1,ray)]; i2=id_of[(v,1,ray+1)]
            p0=points[origin]; p1=points[i1]; p2=points[i2]
            (p1[1]-p0[1])*(p2[2]-p0[2])-(p1[2]-p0[2])*(p2[1]-p0[1])>0 ||
                throw(ArgumentError("$caller: inverted fan triangle at vertex $v"))
            tcursor+=1
            tris[:,tcursor].=(Int32(origin),Int32(i1),Int32(i2))
        end
        for k in 2:nl, ray in 0:nfan-1
            b1=id_of[(v,k-1,ray)]; b2=id_of[(v,k-1,ray+1)]
            t2=id_of[(v,k,ray+1)]; t1=id_of[(v,k,ray)]
            p1=points[b1]; p2=points[t1]; p3=points[t2]; p4=points[b2]
            _signed_quad(p1[1],p1[2],p2[1],p2[2],p3[1],p3[2],p4[1],p4[2])>0 ||
                throw(ArgumentError("$caller: inverted fan quadrangle at vertex $v layer $k"))
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

end # module
