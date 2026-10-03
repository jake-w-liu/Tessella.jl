# ═══════════════════════════════════════════════════════════════════════════
# `Layers` extruded-entity meshing — the sweep kernels of
# meshGEdgeExtruded.cpp / meshGFaceExtruded.cpp / meshGRegionExtruded.cpp.
#
# A `Layers`-attributed extrusion records, per generated entity, the
# `_GeoExtrudeParams` (`m.meshing.extrude`), its upstream source link
# (`m.meshing.extrude_sources` — `(1, ±c)` generatrix for a lateral surface,
# `(2, ±s)` source surface for an extruded volume or its top copy) and the
# `_ExtrudeSpec` transform (`m.meshing.extrude_specs`). The kernels below
# re-evaluate the transform at each layer level's normalized parameter and
# emit the swept element layout upstream produces: lateral surfaces become
# quadrangle strips (recombined) or triangle pairs, top surfaces verbatim
# copies, and volumes prisms/hexahedra (recombined) or the global
# diagonal-compatible three-tetrahedron subdivision upstream's
# `SubdivideExtrudedMesh` computes.
# ═══════════════════════════════════════════════════════════════════════════

# `ExtrudeParams::u(iLayer, iElemLayer)` flattened — the normalized sweep
# parameter at every level boundary. Length `sum(layers)+1`, starting at 0
# and ending at `heights[end]` (upstream does not force it to 1 either).
function _extrude_level_us(params::_GeoExtrudeParams)
    us=Float64[0.0]
    h0=0.0
    for j in eachindex(params.layers)
        n=params.layers[j]
        h1=params.heights[j]
        for k in 1:n
            push!(us,h0+k/n*(h1-h0))
        end
        h0=h1
    end
    return us
end

# `ExtrudeParams::Extrude(t)` — the rigid sweep transform at the normalized
# parameter: `x += t·T` for `:translate`, `ProtudeXYZ` (rotate about the
# axis point by `angle·t`) for `:rotate`, and `ProtudeXYZ` followed by
# `x += t·T` for `:translate_rotate`.
function _extrude_at(spec,p::NTuple{3,Float64},t::Float64,
                     caller::AbstractString)
    if spec.type===:translate
        T=spec.T
        return (p[1]+t*T[1],p[2]+t*T[2],p[3]+t*T[3])
    elseif spec.type===:rotate
        return _extrude_rotate_point(spec.axis,spec.origin,spec.angle*t,p)
    end
    T=spec.T
    p=_extrude_rotate_point(spec.axis,spec.origin,spec.angle*t,p)
    return _gmsh_matvec4x4(
        _gmsh_translation_step((t*T[1],t*T[2],t*T[3])),p)
end

# Rigid transforms depend on the layer level, so build O(levels) records
# and reuse them at every node instead of allocating per point evaluation.
function _extrude_level_transforms(spec,us,caller::AbstractString)
    spec.type===:translate && return nothing
    if spec.type===:rotate
        transforms=Vector{NTuple{3,NTuple{16,Float64}}}(undef,length(us))
        for l in eachindex(us)
            transforms[l]=_extrude_rotation_steps(
                spec.axis,spec.origin,spec.angle*us[l])
        end
        return transforms
    end
    transforms=Vector{NTuple{4,NTuple{16,Float64}}}(undef,length(us))
    for l in eachindex(us)
        u=us[l]
        steps=_extrude_rotation_steps(spec.axis,spec.origin,spec.angle*u)
        transforms[l]=(steps...,_gmsh_translation_step(
            (u*spec.T[1],u*spec.T[2],u*spec.T[3])))
    end
    return transforms
end

@inline function _extrude_level_at(spec,p::NTuple{3,Float64},u::Float64,
                                   transforms,l::Int)
    if transforms===nothing
        return (p[1]+u*spec.T[1],p[2]+u*spec.T[2],p[3]+u*spec.T[3])
    end
    for step in transforms[l]
        p=_gmsh_matvec4x4(step,p)
    end
    return p
end

# The evaluated node chain `_model_curve_mesh_parts` emits for `curve` —
# identical params, evaluator, and vertex snaps, so swept columns weld
# bitwise onto the curve's own line elements. Lazily grades on demand like
# `_model_surface_mesh_curves!` (upstream's global 1-D pass always ran).
function _extrude_curve_nodes(m::GeoModel,curve::Int,caller::AbstractString)
    haskey(m.curves,curve) || throw(ArgumentError(
        "$caller: unknown Curve[$curve]"))
    params=get(m.curve_params,curve,nothing)
    if params===nothing
        _model_surface_mesh_curves!(m,Int[curve],caller)
        params=get(m.curve_params,curve,nothing)
        params===nothing && throw(ArgumentError(
            "$caller: Curve[$curve] could not be meshed"))
    end
    a,b=m.curves[curve]
    closed=a==b
    us=Float64.(params)
    if closed && length(us)>1
        _,hi=_model_curve_param_bounds(m,curve,caller)
        us[end]>=hi && (us=us[1:end-1])
    end
    nodes=[_model_curve_part_point(m,curve,u,caller) for u in us]
    if !isempty(nodes)
        haskey(m.points,a) && (nodes[1]=m.points[a])
        !closed && haskey(m.points,b) && (nodes[end]=m.points[b])
    end
    return nodes
end

# `pos.find` — one entry per bitwise-unique swept coordinate, in
# first-referenced order (upstream's rtree vertex identity).
@inline function _extrude_node!(lookup::Dict{NTuple{3,Float64},Int32},
                                coords::Vector{NTuple{3,Float64}},
                                p::NTuple{3,Float64})
    i=get(lookup,p,Int32(0))
    i!=0 && return i
    i=Int32(length(coords)+1)
    lookup[p]=i
    push!(coords,p)
    return i
end

# The sorted position pair keying upstream's `std::pair<MVertex*,MVertex*>`
# edge sets — lexicographic order stands in for pointer order (the only
# ordering requirement is determinism within the subdivision).
@inline _extrude_ekey(a::NTuple{3,Float64},b::NTuple{3,Float64})=
    isless(a,b) ? (a,b) : (b,a)
@inline _extrude_ein(E,a::NTuple{3,Float64},b::NTuple{3,Float64})=
    _extrude_ekey(a,b) in E

# `createQuaTri` verbatim — degenerate collapses first, then the recombine /
# constrained-diagonal emission order. `edges` is the global subdivision
# diagonal set (`constrainedEdges` upstream); vertex identity is the
# bitwise position tuple.
function _extrude_quatri!(tris,quads,v0,v1,v2,v3,recombine::Bool,
                          edges,t::Int,caller::AbstractString)
    if v0==v1 || v1==v3
        push!(tris,(v0,v3,v2))
        return nothing
    elseif v0==v2 || v2==v3
        push!(tris,(v0,v1,v3))
        return nothing
    elseif v0==v3 || v1==v2
        throw(ErrorException(
            "$caller: Incoherent extruded quadrangle in surface $t"))
    end
    if recombine
        if edges===nothing
            push!(quads,(v0,v1,v3,v2))
        elseif _extrude_ein(edges,v1,v2)
            push!(tris,(v2,v1,v0))
            push!(tris,(v2,v3,v1))
        elseif _extrude_ein(edges,v0,v3)
            push!(tris,(v2,v3,v0))
            push!(tris,(v0,v3,v1))
        else
            push!(quads,(v0,v1,v3,v2))
        end
        return nothing
    end
    if edges===nothing
        push!(tris,(v0,v1,v3))
        push!(tris,(v0,v3,v2))
        return nothing
    end
    if _extrude_ein(edges,v1,v2)
        push!(tris,(v2,v1,v0))
        push!(tris,(v2,v3,v1))
    else
        push!(tris,(v2,v3,v0))
        push!(tris,(v0,v3,v1))
    end
    return nothing
end

# Emit the swept surface cell lists as a `Mesh` (triangles only) or
# `MixedMesh` (with quadrangles). Node order is generatrix-node × level,
# matching upstream's vertex-creation order.
function _extrude_surface_part(tris,quads,coords,lookup)
    out=Matrix{Float64}(undef,3,length(coords))
    for (i,(x,y,z)) in enumerate(coords)
        out[1,i]=x;out[2,i]=y;out[3,i]=z
    end
    # Parts carry zero cell tags — physical ownership is re-derived from the
    # model downstream, the same contract `mesh_model_surface` obeys.
    if !isempty(quads)
        blocks=ElementBlock[]
        if !isempty(tris)
            nodes=Matrix{Int32}(undef,3,length(tris))
            for (i,cell) in enumerate(tris)
                nodes[1,i]=lookup[cell[1]]
                nodes[2,i]=lookup[cell[2]]
                nodes[3,i]=lookup[cell[3]]
            end
            push!(blocks,ElementBlock(2,nodes))
        end
        qnodes=Matrix{Int32}(undef,4,length(quads))
        for (i,cell) in enumerate(quads)
            for k in 1:4
                qnodes[k,i]=lookup[cell[k]]
            end
        end
        push!(blocks,ElementBlock(3,qnodes))
        return MixedMesh(out,blocks)
    end
    tnodes=Matrix{Int32}(undef,3,length(tris))
    for (i,cell) in enumerate(tris)
        tnodes[1,i]=lookup[cell[1]]
        tnodes[2,i]=lookup[cell[2]]
        tnodes[3,i]=lookup[cell[3]]
    end
    return Mesh(out;tris=tnodes)
end

# `extrudeMesh(GEdge *from, GFace *to, ...)` — sweep the generatrix's stored
# node chain into the lateral grid. Upstream queries every cell corner
# through `pos.find(Extrude(u, p))` — the transform evaluated at each level
# — so the grid is the pure evaluation image of the generatrix chain:
# interior-level nodes, connector-curve nodes and the chapeau row all share
# the same eval positions and weld bitwise when the parts merge.
# `edges` carries the global subdivision constraints when
# `SubdivideExtrudedMesh` remeshes the surface.
function _extrude_lateral_mesh(m::GeoModel,t::Int,params::_GeoExtrudeParams,
                               spec,gen_signed::Int,caller::AbstractString;
                               edges=nothing)
    gen=abs(gen_signed)
    chain=_extrude_curve_nodes(m,gen,caller)
    a,b=m.curves[gen]
    closed=a==b
    nnode=length(chain)
    nseg=closed ? nnode : nnode-1
    nseg>=1 || throw(ArgumentError(
        "$caller: Surface[$t] generatrix Curve[$gen] has no segments"))
    us=_extrude_level_us(params)
    nlev=length(us)
    transforms=_extrude_level_transforms(spec,us,caller)
    grid=Matrix{NTuple{3,Float64}}(undef,nnode,nlev)
    for i in 1:nnode
        grid[i,1]=chain[i]
        for l in 2:nlev
            grid[i,l]=_extrude_level_at(spec,chain[i],us[l],transforms,l)
        end
    end
    for r in _extrude_quadtri_regions(m,t,caller)
        seam=_extrude_quadtri_seam_index(m,r,caller)
        seam===nothing && continue
        for i in 1:nnode
            grid[i,nlev]=_extrude_quadtri_snap_seam(seam,grid[i,nlev],caller)
        end
        break
    end
    tris=NTuple{3,NTuple{3,Float64}}[]
    quads=NTuple{4,NTuple{3,Float64}}[]
    recombine=_extrude_quadtri_lateral_recombine(m,t,params,caller)
    for i in 1:nseg,l in 1:nlev-1
        i2=closed ? mod1(i+1,nnode) : i+1
        _extrude_quatri!(tris,quads,grid[i,l],grid[i2,l],
                         grid[i,l+1],grid[i2,l+1],recombine,edges,t,caller)
    end
    # `pos` registration order: per generatrix node, per level.
    lookup=Dict{NTuple{3,Float64},Int32}()
    coords=NTuple{3,Float64}[]
    for i in 1:nnode,l in 1:nlev
        _extrude_node!(lookup,coords,grid[i,l])
    end
    return _extrude_surface_part(tris,quads,coords,lookup)
end

# `copyMesh(GFace *from, GFace *to, ...)` — the top surface duplicates the
# source's mesh verbatim: every node is `pos.find(Extrude(u_top, p))`, the
# transform evaluated at the last layer level, which welds bitwise onto the
# chapeau curves' nodes (themselves transform evaluations).
function _extrude_top_mesh(m::GeoModel,t::Int,params::_GeoExtrudeParams,
                           spec,src_signed::Int,caller::AbstractString;
                           min_angle_deg::Real=25.0,
                           max_periodic_passes=8,
                           size_field::Union{Nothing,AbstractSizeField}=nothing)
    src=abs(src_signed)
    srcmesh=mesh_model_surface(m,src;min_angle_deg=min_angle_deg,
        max_periodic_passes=max_periodic_passes,size_field=size_field)
    u_top=_extrude_level_us(params)[end]
    transforms=_extrude_level_transforms(spec,(u_top,),caller)
    scoords=srcmesh.coords
    lookup=Dict{NTuple{3,Float64},Int32}()
    coords=NTuple{3,Float64}[]
    remap=Vector{Int32}(undef,nnodes(srcmesh))
    for i in 1:nnodes(srcmesh)
        p=(scoords[1,i],scoords[2,i],scoords[3,i])
        remap[i]=_extrude_node!(lookup,coords,
            _extrude_level_at(spec,p,u_top,transforms,1))
    end
    if srcmesh isa Mesh
        ntri=ntris(srcmesh)
        nodes=Matrix{Int32}(undef,3,ntri)
        for i in 1:ntri
            nodes[1,i]=remap[srcmesh.tris[1,i]]
            nodes[2,i]=remap[srcmesh.tris[2,i]]
            nodes[3,i]=remap[srcmesh.tris[3,i]]
        end
        out=Matrix{Float64}(undef,3,length(coords))
        for (i,(x,y,z)) in enumerate(coords)
            out[1,i]=x;out[2,i]=y;out[3,i]=z
        end
        return Mesh(out;tris=nodes)
    end
    if _extrude_quadtri_top(m,t,params,caller)
        return _extrude_quadtri_top_part(m,src,srcmesh,coords,remap,caller)
    end
    blocks=ElementBlock[]
    for block in srcmesh.blocks
        block isa ElementBlock && msh_dimension(block.msh)==2 || continue
        width=msh_num_nodes(block.msh)
        width in (3,4) || continue
        nodes=Matrix{Int32}(undef,width,size(block.nodes,2))
        for i in axes(block.nodes,2),k in 1:width
            nodes[k,i]=remap[block.nodes[k,i]]
        end
        push!(blocks,ElementBlock(block.msh,nodes))
    end
    out=Matrix{Float64}(undef,3,length(coords))
    for (i,(x,y,z)) in enumerate(coords)
        out[1,i]=x;out[2,i]=y;out[3,i]=z
    end
    return MixedMesh(out,blocks)
end

# Per-source-node level columns for the volume sweep — upstream's
# `getExtrudedVertices`: every column entry is `pos.find(Extrude(u, p))`,
# i.e. the transform evaluated at each level boundary. Boundary welding
# falls out of position equality: connector/chapeau mesh nodes and the top
# copy vertices are themselves transform evaluations (`extrudeMesh` /
# `copyMesh`), so the level rows land bitwise on the entity nodes. A vertex
# on the extrusion-fixed set (e.g. a collapsed rotation connector) keeps a
# constant column because its evals coincide, exactly like upstream's
# `pos.find` folding the duplicate position onto the same vertex.
function _extrude_volume_columns(m::GeoModel,t::Int,src::Int,srcmesh,
                                 params::_GeoExtrudeParams,spec,us,
                                 caller::AbstractString)
    nlev=length(us)
    nn=nnodes(srcmesh)
    scoords=srcmesh.coords
    transforms=_extrude_level_transforms(spec,us,caller)
    cols=Matrix{NTuple{3,Float64}}(undef,nlev,nn)
    for i in 1:nn
        p=(scoords[1,i],scoords[2,i],scoords[3,i])
        cols[1,i]=p
        for l in 2:nlev
            cols[l,i]=_extrude_level_at(spec,p,us[l],transforms,l)
        end
    end
    seam=_extrude_quadtri_seam_index(m,t,caller)
    if seam!==nothing
        for i in 1:nn
            cols[nlev,i]=_extrude_quadtri_snap_seam(seam,cols[nlev,i],caller)
        end
    end
    return cols
end

# `createPriPyrTet` verbatim — vertex identity is the bitwise position
# tuple, so `v[i]==v[i+3]` marks a collapsed column.
function _extrude_pripyrtet!(tets,pyramids,prisms,
                             v::NTuple{6,NTuple{3,Float64}})
    dup=Int[]
    for i in 1:3
        v[i]==v[i+3] && push!(dup,i)
    end
    if length(dup)==2
        if dup[1]==1 && dup[2]==2
            push!(tets,(v[1],v[2],v[3],v[6]))
        elseif dup[1]==2 && dup[2]==3
            push!(tets,(v[1],v[2],v[3],v[4]))
        else
            push!(tets,(v[1],v[2],v[3],v[5]))
        end
    elseif length(dup)==1
        if dup[1]==1
            push!(pyramids,(v[2],v[5],v[6],v[3],v[1]))
        elseif dup[1]==2
            push!(pyramids,(v[1],v[3],v[6],v[4],v[2]))
        else
            push!(pyramids,(v[1],v[2],v[5],v[4],v[3]))
        end
    else
        push!(prisms,v)
    end
    return nothing
end

# `createHexPri` verbatim — `v[i]==v[i+4]` marks a collapsed column; an
# unmatched two-collapse pattern emits upstream's error and no cell.
function _extrude_hexpri!(hexes,prisms,degenerate,
                          v::NTuple{8,NTuple{3,Float64}})
    dup=Int[]
    for i in 1:4
        v[i]==v[i+4] && push!(dup,i)
    end
    if length(dup)==2
        if dup[1]==1 && dup[2]==2
            push!(prisms,(v[1],v[4],v[8],v[2],v[3],v[7]))
        elseif dup[1]==2 && dup[2]==3
            push!(prisms,(v[1],v[2],v[5],v[4],v[3],v[8]))
        elseif dup[1]==3 && dup[2]==4
            push!(prisms,(v[1],v[4],v[5],v[2],v[3],v[6]))
        elseif dup[1]==1 && dup[2]==4
            push!(prisms,(v[1],v[2],v[6],v[4],v[3],v[7]))
        else
            # "Wrong hexahedron in extrusion" — no cell, like upstream.
            degenerate[]=true
        end
    else
        push!(hexes,v)
    end
    return nothing
end

# One swept volume: node columns plus the emitted/pending cell records in
# upstream's storage split (tetrahedra, hexahedra, prisms, pyramids).
mutable struct _ExtrudeVolumeSweep
    tag::Int
    recombine::Bool
    cols::Matrix{NTuple{3,Float64}}
    prism6::Vector{NTuple{6,NTuple{3,Float64}}}
    tets::Vector{NTuple{4,NTuple{3,Float64}}}
    hexes::Vector{NTuple{8,NTuple{3,Float64}}}
    prisms::Vector{NTuple{6,NTuple{3,Float64}}}
    pyramids::Vector{NTuple{5,NTuple{3,Float64}}}
    interior::Vector{NTuple{3,Float64}}
    degenerate::Vector{Bool}
end

# `extrudeMesh(GFace *from, GRegion *to, ...)` — sweep the source surface's
# cells into per-level records. Recombined volumes emit through
# `createPriPyrTet`/`createHexPri` directly; non-recombined sweeps keep the
# 6-vertex prism records for the global diagonal pass upstream runs in
# `SubdivideExtrudedMesh` (its initial `createPriPyrTet` output is deleted
# before `phase3`, so nothing is emitted here).
function _extrude_volume_sweep(m::GeoModel,t::Int,params::_GeoExtrudeParams,
                               spec,src_signed::Int,caller::AbstractString)
    src=abs(src_signed)
    srcmesh=mesh_model_surface(m,src)
    us=_extrude_level_us(params)
    nlev=length(us)
    cols=_extrude_volume_columns(m,t,src,srcmesh,params,spec,us,caller)
    sweep=_ExtrudeVolumeSweep(t,params.recombine,cols,
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false])
    tricells=NTuple{3,Int}[]
    quadcells=NTuple{4,Int}[]
    if srcmesh isa Mesh
        for i in 1:ntris(srcmesh)
            push!(tricells,(srcmesh.tris[1,i],srcmesh.tris[2,i],
                            srcmesh.tris[3,i]))
        end
    else
        for block in srcmesh.blocks
            block isa ElementBlock && msh_dimension(block.msh)==2 ||
                continue
            for i in axes(block.nodes,2)
                if block.msh==2
                    push!(tricells,(block.nodes[1,i],block.nodes[2,i],
                                    block.nodes[3,i]))
                elseif block.msh==3
                    push!(quadcells,(block.nodes[1,i],block.nodes[2,i],
                                     block.nodes[3,i],block.nodes[4,i]))
                end
            end
        end
    end
    recombine=params.recombine
    for cell in tricells,l in 1:nlev-1
        c1,c2,c3=cell
        v=(cols[l,c1],cols[l,c2],cols[l,c3],
           cols[l+1,c1],cols[l+1,c2],cols[l+1,c3])
        if recombine && params.quad_to_tri===:add_verts
            push!(sweep.prisms,v)
        elseif recombine
            _extrude_pripyrtet!(sweep.tets,sweep.pyramids,sweep.prisms,v)
        else
            push!(sweep.prism6,v)
        end
    end
    for cell in quadcells,l in 1:nlev-1
        c1,c2,c3,c4=cell
        v=(cols[l,c1],cols[l,c2],cols[l,c3],cols[l,c4],
           cols[l+1,c1],cols[l+1,c2],cols[l+1,c3],cols[l+1,c4])
        if recombine && params.quad_to_tri===:add_verts
            push!(sweep.hexes,v)
        elseif recombine
            _extrude_hexpri!(sweep.hexes,sweep.prisms,sweep.degenerate,v)
        else
            throw(ArgumentError(
                "$caller: Cannot extrude quadrangles without Recombine"))
        end
    end
    return sweep
end

# `phase1` — pick one diagonal per prism lateral face. The first pass only
# takes a face's "other" diagonal when the default is unclaimed; the retry
# (`ntry == 2`, the Michel Benhamou split) picks by vertex order — here the
# lexicographic position order, the deterministic counterpart of `getNum`.
function _extrude_phase1!(records,edges,ntry::Int)
    for v in records
        p1,p2,p3,p4,p5,p6=v
        if ntry==1
            _extrude_ein(edges,p1,p5) ||
                push!(edges,_extrude_ekey(p2,p4))
            _extrude_ein(edges,p5,p3) ||
                push!(edges,_extrude_ekey(p2,p6))
            _extrude_ein(edges,p4,p3) ||
                push!(edges,_extrude_ekey(p1,p6))
        else
            push!(edges,isless(p2,p1) ? _extrude_ekey(p2,p4) :
                                        _extrude_ekey(p1,p5))
            push!(edges,isless(p2,p3) ? _extrude_ekey(p2,p6) :
                                        _extrude_ekey(p5,p3))
            push!(edges,isless(p1,p3) ? _extrude_ekey(p1,p6) :
                                        _extrude_ekey(p4,p3))
        end
    end
    return nothing
end

# `phase2` — a prism whose three lateral diagonals form either "twisted"
# combination cannot split into three tets; rotate one diagonal to its
# partner until no prism remains twisted. `swapset` keeps already-rotated
# diagonals from oscillating back, exactly like `edges_swap`.
function _extrude_phase2!(records,edges,swapset)
    swap=0
    for v in records
        p1,p2,p3,p4,p5,p6=v
        if _extrude_ein(edges,p4,p2) && _extrude_ein(edges,p5,p3) &&
           _extrude_ein(edges,p1,p6)
            swap+=1
            if !_extrude_ein(swapset,p4,p2)
                delete!(edges,_extrude_ekey(p4,p2))
                push!(edges,_extrude_ekey(p1,p5))
                push!(swapset,_extrude_ekey(p4,p2))
                push!(swapset,_extrude_ekey(p1,p5))
            elseif !_extrude_ein(swapset,p5,p3)
                delete!(edges,_extrude_ekey(p5,p3))
                push!(edges,_extrude_ekey(p2,p6))
                push!(swapset,_extrude_ekey(p5,p3))
                push!(swapset,_extrude_ekey(p2,p6))
            elseif !_extrude_ein(swapset,p1,p6)
                delete!(edges,_extrude_ekey(p1,p6))
                push!(edges,_extrude_ekey(p4,p3))
                push!(swapset,_extrude_ekey(p1,p6))
                push!(swapset,_extrude_ekey(p4,p3))
            end
        elseif _extrude_ein(edges,p1,p5) && _extrude_ein(edges,p2,p6) &&
               _extrude_ein(edges,p4,p3)
            swap+=1
            if !_extrude_ein(swapset,p1,p5)
                delete!(edges,_extrude_ekey(p1,p5))
                push!(edges,_extrude_ekey(p4,p2))
                push!(swapset,_extrude_ekey(p1,p5))
                push!(swapset,_extrude_ekey(p4,p2))
            elseif !_extrude_ein(swapset,p2,p6)
                delete!(edges,_extrude_ekey(p2,p6))
                push!(edges,_extrude_ekey(p5,p3))
                push!(swapset,_extrude_ekey(p2,p6))
                push!(swapset,_extrude_ekey(p5,p3))
            elseif !_extrude_ein(swapset,p4,p3)
                delete!(edges,_extrude_ekey(p4,p3))
                push!(edges,_extrude_ekey(p1,p6))
                push!(swapset,_extrude_ekey(p4,p3))
                push!(swapset,_extrude_ekey(p1,p6))
            end
        end
    end
    return swap
end

# `createTet` — emit only cells on four distinct vertices.
@inline function _extrude_tet!(tets,a,b,c,d)
    (a==b || a==c || a==d || b==c || b==d || c==d) ||
        push!(tets,(a,b,c,d))
    return nothing
end

# `phase3` — emit the three tetrahedra per prism for whichever of the six
# valid lateral-diagonal combinations the edge set holds.
function _extrude_phase3!(sweep::_ExtrudeVolumeSweep,edges)
    tets=sweep.tets
    for v in sweep.prism6
        p1,p2,p3,p4,p5,p6=v
        if _extrude_ein(edges,p4,p2) && _extrude_ein(edges,p5,p3) &&
           _extrude_ein(edges,p4,p3)
            _extrude_tet!(tets,p1,p2,p3,p4)
            _extrude_tet!(tets,p4,p5,p6,p3)
            _extrude_tet!(tets,p2,p4,p5,p3)
        elseif _extrude_ein(edges,p4,p2) && _extrude_ein(edges,p2,p6) &&
               _extrude_ein(edges,p4,p3)
            _extrude_tet!(tets,p1,p2,p3,p4)
            _extrude_tet!(tets,p4,p5,p6,p2)
            _extrude_tet!(tets,p4,p2,p6,p3)
        elseif _extrude_ein(edges,p4,p2) && _extrude_ein(edges,p2,p6) &&
               _extrude_ein(edges,p6,p1)
            _extrude_tet!(tets,p1,p2,p3,p6)
            _extrude_tet!(tets,p4,p5,p6,p2)
            _extrude_tet!(tets,p2,p4,p6,p1)
        elseif _extrude_ein(edges,p5,p1) && _extrude_ein(edges,p5,p3) &&
               _extrude_ein(edges,p4,p3)
            _extrude_tet!(tets,p1,p2,p3,p5)
            _extrude_tet!(tets,p4,p5,p6,p3)
            _extrude_tet!(tets,p1,p4,p5,p3)
        elseif _extrude_ein(edges,p5,p1) && _extrude_ein(edges,p5,p3) &&
               _extrude_ein(edges,p6,p1)
            _extrude_tet!(tets,p1,p2,p3,p5)
            _extrude_tet!(tets,p4,p5,p6,p1)
            _extrude_tet!(tets,p1,p3,p5,p6)
        elseif _extrude_ein(edges,p5,p1) && _extrude_ein(edges,p2,p6) &&
               _extrude_ein(edges,p6,p1)
            _extrude_tet!(tets,p1,p2,p3,p6)
            _extrude_tet!(tets,p4,p5,p6,p1)
            _extrude_tet!(tets,p1,p2,p5,p6)
        end
    end
    return nothing
end

# `SubdivideExtrudedMesh` — global diagonal selection over every swept,
# non-recombined volume: `phase1` seeds laterals, `phase2` untwists to a
# fixpoint (one retry on a stalled cycle, exactly like `ntry`), then
# `phase3` emits the tets. Returns the shared diagonal set the lateral
# surfaces remesh against.
function _extrude_subdivide!(sweeps::Vector{_ExtrudeVolumeSweep},
                             caller::AbstractString)
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    solved=false
    for ntry in 1:2
        empty!(edges)
        for sweep in sweeps
            _extrude_phase1!(sweep.prism6,edges,ntry)
        end
        swapset=Set{NTuple{2,NTuple{3,Float64}}}()
        prev=0
        stall=false
        while true
            swap=0
            for sweep in sweeps
                swap+=_extrude_phase2!(sweep.prism6,edges,swapset)
            end
            swap==0 && (solved=true;break)
            if prev!=0 && prev==swap
                ntry==1 && (stall=true;break)
                throw(ErrorException(
                    "$caller: unable to subdivide extruded mesh: change " *
                    "surface mesh or recombine extrusion instead"))
            end
            prev=swap
        end
        solved && break
        stall || break
    end
    solved || throw(ErrorException(
        "$caller: unable to subdivide extruded mesh: change surface " *
        "mesh or recombine extrusion instead"))
    for sweep in sweeps
        _extrude_phase3!(sweep,edges)
    end
    return edges
end

# `GModel::setAllVolumesPositive` — `MElement::setVolumePositive` reverses
# any negative-Jacobian volume cell. The swaps are each element's
# `reverse()`: tet and prism exchange the first two corners of both
# triangle/quad caps; hexahedron and pyramid exchange corner 1 with corner 3
# (plus the matching top pair for the hex).
const _EXTRUDE_REVERSE_SWAPS=Dict{Int,NTuple{2,NTuple{2,Int}}}(
    4=>((1,2),(0,0)),5=>((1,3),(5,7)),
    6=>((1,2),(4,5)),7=>((1,3),(0,0)))

# An orientation-consistent tetrahedralization per volume cell type — the
# hex fan along the (1,7) body diagonal, the prism fan along (1,6), and the
# pyramid base split on (1,3). Signed `tet_signed_volume` sums give the
# cell's Jacobian sign, the quantity `setVolumePositive` tests.
const _EXTRUDE_CELL_TETS=Dict{Int,Vector{NTuple{4,Int}}}(
    4=>[(1,2,3,4)],
    5=>[(1,2,3,7),(1,3,4,7),(1,4,8,7),(1,8,5,7),(1,5,6,7),(1,6,2,7)],
    6=>[(1,2,3,6),(1,2,6,5),(1,5,6,4)],
    7=>[(1,2,3,5),(1,3,4,5)])

@inline function _extrude_cell_pt(coords::Matrix{Float64},
                                  nodes::Matrix{Int32},i::Int,cell::Int)
    v=nodes[i,cell]
    return (coords[1,v],coords[2,v],coords[3,v])
end

function _extrude_cell_signed_volume(coords::Matrix{Float64},
                                     nodes::Matrix{Int32},cell::Int,msh::Int)
    total=0.0
    @inbounds for (i1,i2,i3,i4) in _EXTRUDE_CELL_TETS[msh]
        total+=tet_signed_volume(
            _extrude_cell_pt(coords,nodes,i1,cell),
            _extrude_cell_pt(coords,nodes,i2,cell),
            _extrude_cell_pt(coords,nodes,i3,cell),
            _extrude_cell_pt(coords,nodes,i4,cell))
    end
    return total
end

function _extrude_make_positive!(coords::Matrix{Float64},
                                 nodes::Matrix{Int32},msh::Int)
    (a1,b1),(a2,b2)=_EXTRUDE_REVERSE_SWAPS[msh]
    @inbounds for cell in axes(nodes,2)
        volume=_extrude_cell_signed_volume(coords,nodes,cell,msh)
        isfinite(volume) && volume!=0 || throw(ArgumentError(
            "extruded mesh: degenerate or nonfinite volume in element type $msh"))
        volume>0 && continue
        nodes[a1,cell],nodes[b1,cell]=nodes[b1,cell],nodes[a1,cell]
        a2==0 && continue
        nodes[a2,cell],nodes[b2,cell]=nodes[b2,cell],nodes[a2,cell]
    end
    return nothing
end

# Emit a swept volume as a `Mesh` (tetrahedra only) or `MixedMesh`. Nodes
# are the deduplicated column positions in source-node × level order, and
# blocks follow ascending MSH type like the writer's serialization order.
function _extrude_volume_part(m::GeoModel,sweep::_ExtrudeVolumeSweep,
                              caller::AbstractString)
    sweep.degenerate[] && throw(ArgumentError(
        "$caller: Wrong hexahedron in extrusion"))
    lookup=Dict{NTuple{3,Float64},Int32}()
    coords=NTuple{3,Float64}[]
    cols=sweep.cols
    for i in axes(cols,2),l in axes(cols,1)
        _extrude_node!(lookup,coords,cols[l,i])
    end
    for p in sweep.interior
        _extrude_node!(lookup,coords,p)
    end
    # Transition fans introduce body vertices after the column sweep.
    for cells in (sweep.tets,sweep.hexes,sweep.prisms,sweep.pyramids),cell in cells,p in cell
        _extrude_node!(lookup,coords,p)
    end
    out=Matrix{Float64}(undef,3,length(coords))
    for (i,(x,y,z)) in enumerate(coords)
        out[1,i]=x;out[2,i]=y;out[3,i]=z
    end
    if isempty(sweep.hexes) && isempty(sweep.prisms) &&
       isempty(sweep.pyramids)
        ncells=length(sweep.tets)
        nodes=Matrix{Int32}(undef,4,ncells)
        for (i,cell) in enumerate(sweep.tets)
            nodes[1,i]=lookup[cell[1]]
            nodes[2,i]=lookup[cell[2]]
            nodes[3,i]=lookup[cell[3]]
            nodes[4,i]=lookup[cell[4]]
        end
        _extrude_make_positive!(out,nodes,4)
        return Mesh(out;tets=nodes)
    end
    blocks=ElementBlock[]
    if !isempty(sweep.tets)
        nodes=Matrix{Int32}(undef,4,length(sweep.tets))
        for (i,cell) in enumerate(sweep.tets),k in 1:4
            nodes[k,i]=lookup[cell[k]]
        end
        _extrude_make_positive!(out,nodes,4)
        push!(blocks,ElementBlock(4,nodes))
    end
    if !isempty(sweep.hexes)
        nodes=Matrix{Int32}(undef,8,length(sweep.hexes))
        for (i,cell) in enumerate(sweep.hexes),k in 1:8
            nodes[k,i]=lookup[cell[k]]
        end
        _extrude_make_positive!(out,nodes,5)
        push!(blocks,ElementBlock(5,nodes))
    end
    if !isempty(sweep.prisms)
        nodes=Matrix{Int32}(undef,6,length(sweep.prisms))
        for (i,cell) in enumerate(sweep.prisms),k in 1:6
            nodes[k,i]=lookup[cell[k]]
        end
        _extrude_make_positive!(out,nodes,6)
        push!(blocks,ElementBlock(6,nodes))
    end
    if !isempty(sweep.pyramids)
        nodes=Matrix{Int32}(undef,5,length(sweep.pyramids))
        for (i,cell) in enumerate(sweep.pyramids),k in 1:5
            nodes[k,i]=lookup[cell[k]]
        end
        _extrude_make_positive!(out,nodes,7)
        push!(blocks,ElementBlock(7,nodes))
    end
    return MixedMesh(out,blocks)
end

# The `(dim,tag)` gate for a swept entity — `ep->mesh.ExtrudeMesh` is set
# only by `Layers`, so an extrusion without layers keeps ordinary meshing.
function _extrude_gate(m::GeoModel,dim::Int,tag::Int)
    params=get(m.meshing.extrude,(dim,tag),nothing)
    params===nothing && return nothing
    isempty(params.layers) && return nothing
    return params
end

function _extrude_entity_params(m::GeoModel,dim::Int,tag::Int,
                                caller::AbstractString)
    params=_extrude_gate(m,dim,tag)
    params===nothing && return nothing
    if params.quad_to_tri===:no_new_verts && params.recombine &&
       (dim==3 || any(r->_extrude_is_quadtri(_extrude_gate(m,3,r)),
                     _extrude_quadtri_regions(m,tag,caller)))
        kind="QuadTriNoNewVerts"
        throw(ArgumentError(
            "$caller: $kind on $(dim==2 ? "Surface" : "Volume")[$tag] " *
            "requires the QuadToTri extrusion kernel, which Tessella " *
            "does not implement"))
    end
    spec=get(m.meshing.extrude_specs,(dim,tag),nothing)
    link=get(m.meshing.extrude_sources,(dim,tag),nothing)
    (spec===nothing || link===nothing) && throw(ArgumentError(
        "$caller: $(dim==2 ? "Surface" : "Volume")[$tag] Layers " *
        "extrusion is missing recorded extrusion parameters"))
    return params,spec,link
end

# `MeshExtrudedSurface` — dispatch on `geo.Mode`: `(1, ·)` sources are
# EXTRUDED_ENTITY laterals, `(2, ·)` sources are COPIED_ENTITY tops.
function _extrude_surface_mesh(m::GeoModel,t::Int,caller::AbstractString;
                               min_angle_deg::Real=25.0,
                               max_periodic_passes=8,
                               size_field::Union{Nothing,AbstractSizeField}=nothing,
                               edges=nothing)
    entry=_extrude_entity_params(m,2,t,caller)
    entry===nothing && return nothing
    params,spec,link=entry
    if link[1]==1
        return _extrude_lateral_mesh(m,t,params,spec,link[2],caller;
                                   edges=edges)
    end
    return _extrude_top_mesh(m,t,params,spec,link[2],caller;
        min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
        size_field=size_field)
end

# `meshGRegionExtruded` + the region-local half of `SubdivideExtrudedMesh` —
# the standalone `mesh_model_volume` arm runs the subdivision over this one
# volume's prism records (identical to the global result when no second
# extruded volume shares its lateral faces).
function _extrude_volume_mesh(m::GeoModel,t::Int,caller::AbstractString)
    entry=_extrude_entity_params(m,3,t,caller)
    entry===nothing && return nothing
    params,spec,link=entry
    sweep=_extrude_volume_sweep(m,t,params,spec,link[2],caller)
    if params.recombine && params.quad_to_tri===:add_verts
        _extrude_quadtri_addverts!(m,sweep,caller)
    elseif !params.recombine
        _extrude_subdivide!(_ExtrudeVolumeSweep[sweep],caller)
    end
    return _extrude_volume_part(m,sweep,caller)
end

# The model-wide extrusion pass upstream interleaves into `Mesh 3`:
# `meshGRegionExtruded` for every swept volume first, then the shared
# `SubdivideExtrudedMesh` diagonal selection, then lateral-surface remesh
# against the global edge set. Returns the per-volume parts plus the
# remesh bookkeeping `_geo_mesh_model` applies to the surface parts.
function _extrude_volume_pass(m::GeoModel,caller::AbstractString)
    tags=Int[]
    for tag in keys(m.volumes)
        _extrude_gate(m,3,tag)===nothing || push!(tags,tag)
    end
    isempty(tags) && return nothing
    sort!(tags)
    sweeps=_ExtrudeVolumeSweep[]
    for tag in tags
        params,spec,link=_extrude_entity_params(m,3,tag,caller)
        push!(sweeps,_extrude_volume_sweep(m,tag,params,spec,link[2],caller))
    end
    subdivided=[s for s in sweeps if !s.recombine]
    edges=nothing
    isempty(subdivided) || (edges=_extrude_subdivide!(subdivided,caller))
    for sweep in sweeps
        params=_extrude_gate(m,3,sweep.tag)
        params.recombine && params.quad_to_tri===:add_verts || continue
        _extrude_quadtri_addverts!(m,sweep,caller;edges=edges)
    end
    parts=Dict{Int,Union{Mesh,MixedMesh}}()
    for sweep in sweeps
        parts[sweep.tag]=_extrude_volume_part(m,sweep,caller)
    end
    remesh=Set{Int}()
    for sweep in sweeps
        volume_params=_extrude_gate(m,3,sweep.tag)
        quadtri=sweep.recombine && volume_params.quad_to_tri===:add_verts
        sweep.recombine && !(quadtri && edges!==nothing) && continue
        for sl in m.volumes[sweep.tag],s in m.surface_loops[sl]
            ep=get(m.meshing.extrude,(2,abs(s)),nothing)
            ep===nothing && continue
            isempty(ep.layers) && continue
            ep.recombine && !quadtri && continue
            link=get(m.meshing.extrude_sources,(2,abs(s)),nothing)
            link===nothing && continue
            link[1]==1 && push!(remesh,abs(s))
        end
    end
    return (parts=parts,edges=edges,remesh=remesh)
end

include("ModelExtrudeQuadTri.jl")
