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
    if spec.type===:boundary_layer
        # The chapeau of a `BOUNDARY_LAYER` curve extrusion has no rigid
        # transform: each generatrix node follows its own field direction to
        # the strip's top row. `spec.tops` maps the evaluated generatrix
        # position to that exact emitted top node; mirror-evaluated keys can
        # drift an ulp off the stored chain, so fall back to the nearest row
        # within a scale-aware tolerance.
        tops=spec.tops
        top=get(tops,p,nothing)
        top===nothing || return top
        isempty(tops) && throw(ArgumentError(
            "$caller: boundary-layer copy has no stored top positions"))
        scale=1.0
        for pos in keys(tops)
            scale=max(scale,abs(pos[1]),abs(pos[2]),abs(pos[3]))
        end
        best=nothing;bd=Inf;bd2=Inf
        for (pos,v) in tops
            d=hypot(p[1]-pos[1],p[2]-pos[2],p[3]-pos[3])
            if d<bd
                bd2=bd;bd=d;best=v
            elseif d<bd2
                bd2=d
            end
        end
        bd<=1e-9*scale && return best
        # The stored keys were evaluated at the generatrix's refresh-time
        # parameters; a later native re-grading can shift the same logical
        # node a few e-9 along the curve (a curved generatrix's u drifts
        # when the 1-D producer re-integrates it). The row correspondence is
        # still unambiguous — the drift stays far inside the key's Voronoi
        # cell — so accept the nearest row when it dominates the runner-up
        # by an order of magnitude and the displacement stays well under
        # any reachable layer offset. A genuinely foreign point sits between
        # keys with best≈second distances and still rejects.
        bd<=1e-7*scale && bd<=0.1*bd2 && return best
        throw(ArgumentError(
            "$caller: boundary-layer copy evaluated $(p) off the wall's " *
            "stored nodes"))
    elseif spec.type===:translate
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
    (spec.type===:translate || spec.type===:boundary_layer) && return nothing
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
    # Exact full turns are the identity — the closing layer row lands
    # bitwise on the source row so the swept seam welds.
    spec.type===:rotate && mod(spec.angle*u,2π)==0.0 && return p
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

# The operation retains its own position registry and indexed source columns.
# Coordinates verify these events; they do not rediscover producer identities.
struct _ExtrudeMeshProduction{M}
    mesh::M
    position_rows::Dict{NTuple{3,Float64},Int32}
    source_node_tags::Vector{UInt64}
    column_rows::Matrix{Int32}
    cell_events::Vector{Tuple{Int,Int}}
end

_ExtrudeMeshProduction(mesh,positions,source_tags,columns)=
    _ExtrudeMeshProduction(mesh,positions,source_tags,columns,Tuple{Int,Int}[])

function _extrude_curve_source(m,curve,input,caller)
    input===nothing && return (_extrude_curve_nodes(m,curve,caller),UInt64[],
        NTuple{3,Float64}[],UInt64[])
    rows,mesh=if m.curves[curve][1]!=m.curves[curve][2]
        source=_model_mesh_curve_births(input,m,curve,caller;used_only=true)
        source!==nothing || throw(ArgumentError("$caller: extrusion lacks actual Curve[$curve] cells"))
        source.rows,source.mesh
    else
        product_mesh=_model_mesh_product(input,1,curve;caller)
        product_mesh!==nothing || throw(ArgumentError("$caller: extrusion lacks actual closed Curve[$curve] cells"))
        # A closed producer starts at its actual Point vertex and lists each
        # subsequent Line vertex once. Check its literal oriented cell cycle.
        cycle_used=Set(row for block in product_mesh.blocks for row in block.nodes)
        produced=get(input.producer_rows,(1,Int32(curve)),nothing)
        cycle_rows=if produced===nothing
            Int32[row for row in 1:nnodes(product_mesh) if row in cycle_used]
        else
            local_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(product_mesh.entity_data.external_node_tags))
            unique(Int32[local_rows[input.mesh.entity_data.external_node_tags[row]] for row in produced
                if get(local_rows,input.mesh.entity_data.external_node_tags[row],Int32(0)) in cycle_used])
        end
        length(cycle_rows)==length(cycle_used) || throw(ArgumentError(
            "$caller: closed Curve[$curve] lost actual producer incidence"))
        all(block->block.msh==1,product_mesh.blocks) || throw(ArgumentError(
            "$caller: extrusion needs linear closed Curve[$curve] cells"))
        edges=Set(minmax(block.nodes[1,c],block.nodes[2,c]) for block in product_mesh.blocks for c in axes(block.nodes,2))
        expected=Set(minmax(cycle_rows[i],cycle_rows[mod1(i+1,length(cycle_rows))]) for i in eachindex(cycle_rows))
        length(edges)==sum(size(block.nodes,2) for block in product_mesh.blocks)==length(cycle_rows) && edges==expected ||
            throw(ArgumentError("$caller: closed Curve[$curve] lacks a complete actual producer cycle"))
        product_mesh.entity_data.node_entities[first(cycle_rows)]==(0,Int32(m.curves[curve][1])) || throw(ArgumentError(
            "$caller: closed Curve[$curve] producer does not start at its actual Point"))
        cycle_rows,product_mesh
    end
    used=Set(rows)
    orphans=Int32[row for row in 1:nnodes(mesh) if !(row in used)]
    all(row->mesh.entity_data.node_entities[row]==(1,Int32(curve)),orphans) || throw(ArgumentError(
        "$caller: Curve[$curve] has foreign unreferenced source rows"))
    return (NTuple{3,Float64}[ntuple(a->mesh.coords[a,row],3) for row in rows],
        mesh.entity_data.external_node_tags[rows],
        NTuple{3,Float64}[ntuple(a->mesh.coords[a,row],3) for row in orphans],
        mesh.entity_data.external_node_tags[orphans])
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

# ── Boundary-layer curve extrusion (`Extrude {Curve{c}; Layers{...}}`) ────────
# Gmsh's `BOUNDARY_LAYER` mode sweeps each generatrix node along the smoothed
# normal field (`boundaryLayerField` / `smooth_data`) rather than a rigid
# transform: source-edge element normals accumulate by position, so a vertex
# shared between sibling walls grows one common column. Tessella delegates the
# strip emission to `BoundaryLayer._bl_2d_offsets` — the same kernel that backs
# `mesh_boundary_layer_2d` — and keeps the sibling entities welding bitwise by
# feeding the kernel's exact per-node directions back into the connector and
# chapeau `extrude_specs` records.

# The wall polygon's plane normal via Newell's method, sign-aligned toward +z
# so the positive-generatrix growth side reproduces upstream's `t × ẑ`
# convention in either winding. Straight (degenerate-winding) walls fall back
# to the tangent-perpendicular component of ẑ, then x̂/ŷ.
function _extrude_bl_wall_normal(chain,closed::Bool,caller::AbstractString)
    nnode=length(chain)
    nseg=closed ? nnode : nnode-1
    cx=cy=cz=0.0
    for p in chain
        cx+=p[1];cy+=p[2];cz+=p[3]
    end
    inv=1.0/nnode;cx*=inv;cy*=inv;cz*=inv
    nx=ny=nz=0.0
    for i in 1:nseg
        j=i==nnode ? 1 : i+1
        ux=chain[i][1]-cx;uy=chain[i][2]-cy;uz=chain[i][3]-cz
        wx=chain[j][1]-cx;wy=chain[j][2]-cy;wz=chain[j][3]-cz
        nx+=uy*wz-uz*wy;ny+=uz*wx-ux*wz;nz+=ux*wy-uy*wx
    end
    L=hypot(nx,ny,nz)
    if L>0
        nx/=L;ny/=L;nz/=L
        nz<0 && return (-nx,-ny,-nz)
        return (nx,ny,nz)
    end
    nseg>=1 || throw(ArgumentError(
        "$caller: boundary-layer wall has no segments"))
    t=(chain[2][1]-chain[1][1],chain[2][2]-chain[1][2],chain[2][3]-chain[1][3])
    Lt=hypot(t...)
    Lt>0 || throw(ArgumentError(
        "$caller: boundary-layer wall has a zero-length first segment"))
    t=(t[1]/Lt,t[2]/Lt,t[3]/Lt)
    for ref in ((0.0,0.0,1.0),(1.0,0.0,0.0),(0.0,1.0,0.0))
        d=t[1]*ref[1]+t[2]*ref[2]+t[3]*ref[3]
        nx,ny,nz=ref[1]-d*t[1],ref[2]-d*t[2],ref[3]-d*t[3]
        Ln=hypot(nx,ny,nz)
        Ln>0.1 && return (nx/Ln,ny/Ln,nz/Ln)
    end
    throw(ArgumentError(
        "$caller: boundary-layer wall has no stable plane normal"))
end

# The kernel's `plane_normal` for a signed generatrix: a reversed record
# traverses the wall backwards, so its growth side flips.
_extrude_bl_plane_normal(chain,closed::Bool,positive::Bool,
                         caller::AbstractString) =
    (pn=_extrude_bl_wall_normal(chain,closed,caller);
     positive ? (-pn[1],-pn[2],-pn[3]) : pn)

# Coherently directed polyline `Mesh` — the `_bl_2d_offsets` input.
function _extrude_bl_polyline(chain,closed::Bool)
    nnode=length(chain)
    nseg=closed ? nnode : nnode-1
    coords=Matrix{Float64}(undef,3,nnode)
    for (i,p) in enumerate(chain)
        coords[1,i]=p[1];coords[2,i]=p[2];coords[3,i]=p[3]
    end
    segs=Matrix{Int32}(undef,2,nseg)
    for i in 1:nseg
        segs[1,i]=Int32(i);segs[2,i]=Int32(i==nnode ? 1 : i+1)
    end
    return Mesh(coords;segs=segs)
end

# `smooth_data::add(x,y,z,3,nn)` — every wall segment's left normal
# (`cross(peff, t)`: `t × ẑ` upstream, generalized to each wall's own plane)
# accumulates on both endpoint positions; the merged per-position sum
# normalizes into the growth direction. Positions whose incident segments
# cancel (opposing walls) stay unset and keep the kernel's own average.
function _extrude_bl_normal_field(walls)
    sums=Dict{NTuple{3,Float64},NTuple{3,Float64}}()
    for wall in walls
        chain,closed,peff=wall
        nnode=length(chain)
        nseg=closed ? nnode : nnode-1
        for i in 1:nseg
            j=i==nnode ? 1 : i+1
            p,q=chain[i],chain[j]
            t=(q[1]-p[1],q[2]-p[2],q[3]-p[3])
            Lt=hypot(t...)
            Lt>0 || continue
            t=(t[1]/Lt,t[2]/Lt,t[3]/Lt)
            n=(peff[2]*t[3]-peff[3]*t[2],
               peff[3]*t[1]-peff[1]*t[3],
               peff[1]*t[2]-peff[2]*t[1])
            for pos in (p,q)
                s=get(sums,pos,(0.0,0.0,0.0))
                sums[pos]=(s[1]+n[1],s[2]+n[2],s[3]+n[3])
            end
        end
    end
    field=Dict{NTuple{3,Float64},NTuple{3,Float64}}()
    for (pos,s) in sums
        L=hypot(s[1],s[2],s[3])
        L>0 || continue
        field[pos]=(s[1]/L,s[2]/L,s[3]/L)
    end
    return field
end

# `(chain, closed, peff)` for every boundary-layer surface's generatrix — the
# mesh-time equivalent of the exec-time entity list, so the field merges
# across `Extrude` statements just as upstream accumulates every
# BoundaryLayerCurve's edges.
function _extrude_bl_walls(m::GeoModel,caller::AbstractString;
                           _mesh_input=nothing)
    walls=Tuple{Vector{NTuple{3,Float64}},Bool,NTuple{3,Float64}}[]
    for (d,tag) in sort!(collect(keys(m.meshing.extrude_specs)))
        d==2 || continue
        spec=m.meshing.extrude_specs[(d,tag)]
        spec isa NamedTuple || continue
        get(spec,:type,nothing)===:boundary_layer || continue
        link=get(m.meshing.extrude_sources,(2,tag),nothing)
        link===nothing && continue
        chain,_,orphans,_=_extrude_curve_source(m,abs(link[2]),_mesh_input,caller)
        isempty(orphans) || throw(ArgumentError(
            "$caller: Surface[$tag] boundary-layer generatrix " *
            "Curve[$(abs(link[2]))] has unreferenced mesh nodes"))
        a,b=m.curves[abs(link[2])]
        closed=a==b
        peff=_extrude_bl_plane_normal(chain,closed,link[2]>0,caller)
        push!(walls,(chain,closed,peff))
    end
    return walls
end

# One wall's strip emission: explicit cumulative offsets plus the merged-field
# direction per shared wall vertex.
function _extrude_bl_kernel(chain,closed::Bool,peff,field,offsets,caller)
    poly=_extrude_bl_polyline(chain,closed)
    vd=Dict{Int,NTuple{3,Float64}}()
    for (i,pos) in enumerate(chain)
        d=get(field,pos,nothing)
        d===nothing || (vd[i]=d)
    end
    # `BoundaryLayer` is a sibling submodule included after `Model` — resolve
    # it lazily like `.geo`'s `SizeField` lookup.
    return getfield(parentmodule(@__MODULE__),:BoundaryLayer)._bl_2d_offsets(
        poly,offsets;plane_normal=peff,vertex_dirs=vd,caller=caller)
end

# Exact-match lookup, else the nearest kernel row — mirrored chapeau evals can
# land an ulp off the stored chain position. `unbounded` skips the tolerance
# for off-wall lookups (the chapeau's interior control copies).
function _extrude_bl_top_at(top_of::AbstractDict,p,caller;unbounded::Bool=false)
    top=get(top_of,p,nothing)
    top===nothing || return top
    isempty(top_of) && throw(ArgumentError(
        "$caller: boundary-layer wall produced no top positions"))
    scale=1.0
    for pos in keys(top_of)
        scale=max(scale,abs(pos[1]),abs(pos[2]),abs(pos[3]))
    end
    best=nothing;bd=Inf
    for (pos,v) in top_of
        d=(p[1]-pos[1])^2+(p[2]-pos[2])^2+(p[3]-pos[3])^2
        d<bd && (bd=d;best=v)
    end
    unbounded || bd<=(1e-9*scale)^2 || throw(ArgumentError(
        "$caller: boundary-layer copy evaluated $(p) off the wall's " *
        "stored nodes"))
    return best
end

# `copyMesh`'s position→top map for the chapeau: every stored chapeau
# parameter evals (`_model_curve_part_point`) to a generatrix position —
# direct for a forward copy, `t1+u_min-v` mirrored for a reversed one — which
# keys the kernel's exact top-row node.
function _extrude_bl_tops_dict(m::GeoModel,gen::Int,chapeau::Int,
                               signed_c::Int,top_of::AbstractDict,
                               caller::AbstractString)
    params=get(m.meshing.extrude,(1,chapeau),nothing)
    params===nothing && return Dict{NTuple{3,Float64},NTuple{3,Float64}}()
    vals=_model_curve_extrude_params(m,chapeau,params,
        _model_default_mesh1d_options(m,caller),caller)
    vals===:pending && throw(ArgumentError(
        "$caller: Curve[$chapeau] boundary-layer source has no stored " *
        "mesh parameters"))
    _,t1=_model_curve_param_bounds(m,chapeau,caller)
    u_min=_model_curve_param_bounds(m,gen,caller)[1]
    tops=Dict{NTuple{3,Float64},NTuple{3,Float64}}()
    for v in vals
        usrc=signed_c<0 ? (t1+u_min)-v : v
        p=_periodic_curve_point(m,gen,usrc,caller)
        tops[p]=_extrude_bl_top_at(top_of,p,caller)
    end
    return tops
end

# Write the kernel's exact directions/positions back onto the sibling
# records: connector curves keep `:translate` specs whose `T` is the unit
# direction (the stored level params are absolute heights, so `p + h·T`
# lands bitwise on the strip's rows), the chapeau keeps the position→top map
# its `_extrude_at` branch reads, and every generated vertex takes the
# kernel's top-row coordinate.
function _extrude_bl_refresh!(m::GeoModel,spec,gen::Int,signed_c::Int,
                              closed::Bool,kernel,
                              caller::AbstractString;
                              _params_model::GeoModel=m)
    a,b=m.curves[gen]
    sbeg,send=signed_c>0 ? (a,b) : (b,a)
    beg_pos=m.points[sbeg];end_pos=m.points[send]
    cbeg,cend=spec.connectors
    # Sibling walls share coincident connectors — the operation's coherence
    # pass keeps one tag; refresh only rewrites survivors.
    if haskey(m.curves,cbeg)
        m.points[m.curves[cbeg][2]]=kernel.tops[beg_pos]
        m.meshing.extrude_specs[(1,cbeg)]=
            (type=:translate,T=kernel.dirs[beg_pos])
    end
    if cend>0 && haskey(m.curves,cend)
        m.points[m.curves[cend][2]]=kernel.tops[end_pos]
        m.meshing.extrude_specs[(1,cend)]=
            (type=:translate,T=kernel.dirs[end_pos])
    end
    chapeau=spec.chapeau
    if haskey(m.curves,chapeau)
        tops=_extrude_bl_tops_dict(_params_model,gen,chapeau,signed_c,
                                   kernel.tops,caller)
        m.meshing.extrude_specs[(1,chapeau)]=(type=:boundary_layer,tops=tops)
        m.points[m.curves[chapeau][1]]=kernel.tops[beg_pos]
        closed || (m.points[m.curves[chapeau][2]]=kernel.tops[end_pos])
        # The chapeau's materialized control copies follow it to the top
        # row — otherwise coincident-with-source copies merge into the
        # source's endpoints and `Coherence` collapses the record back into
        # its generatrix. Interior copies of non-line sources bind to the
        # nearest wall row.
        for cp in get(m.curve_control_points,chapeau,Int[])
            haskey(m.points,cp) || continue
            m.points[cp]=_extrude_bl_top_at(kernel.tops,m.points[cp],caller;
                                            unbounded=true)
        end
    end
    return nothing
end

# Non-recombined product: split each strip quadrangle on the `(b1,b2,t2)` /
# `(b1,t2,t1)` diagonal — the same fan `_extrude_quatri!` emits.
function _extrude_bl_unrecombined(mesh::MixedMesh)
    ntri=0
    for block in mesh.blocks
        ntri+=block.msh==3 ? 2*size(block.nodes,2) :
              block.msh==2 ? size(block.nodes,2) : 0
    end
    tris=Matrix{Int32}(undef,3,ntri)
    c=0
    for block in mesh.blocks
        if block.msh==3
            for i in axes(block.nodes,2)
                b1,b2,t2,t1=block.nodes[:,i]
                c+=1;tris[:,c].=(b1,b2,t2)
                c+=1;tris[:,c].=(b1,t2,t1)
            end
        elseif block.msh==2
            for i in axes(block.nodes,2)
                c+=1;tris[:,c].=block.nodes[:,i]
            end
        end
    end
    return Mesh(mesh.coords;tris=tris)
end

# `extrudeMesh(GEdge *from, GFace *to)` under `BOUNDARY_LAYER`: the wall's
# stored node chain grown along the merged normal field. The kernel returns
# the certified strip; refreshing the sibling specs keeps the connector and
# chapeau evals bitwise-aligned with the emitted columns.
function _extrude_bl_lateral_mesh(m::GeoModel,t::Int,params::_GeoExtrudeParams,
                                  spec,gen_signed::Int,caller::AbstractString;
                                  _mesh_input=nothing)
    gen=abs(gen_signed)
    chain,_,orphans,_=_extrude_curve_source(m,gen,_mesh_input,caller)
    isempty(orphans) || throw(ArgumentError(
        "$caller: Surface[$t] boundary-layer generatrix Curve[$gen] has " *
        "unreferenced mesh nodes"))
    a,b=m.curves[gen]
    closed=a==b
    peff=_extrude_bl_plane_normal(chain,closed,gen_signed>0,caller)
    walls=_extrude_bl_walls(m,caller;_mesh_input)
    field=_extrude_bl_normal_field(walls)
    kernel=_extrude_bl_kernel(chain,closed,peff,field,
        _extrude_level_us(params)[2:end],caller)
    _extrude_bl_refresh!(m,spec,gen,gen_signed,closed,kernel,caller)
    return params.recombine ? kernel.mesh : _extrude_bl_unrecombined(kernel.mesh)
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
                               edges=nothing,_mesh_input=nothing,_return_production=false,
                               _register_orphans=true)
    gen=abs(gen_signed)
    chain,source_tags,orphan_points,orphan_tags=_extrude_curve_source(m,gen,_mesh_input,caller)
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
    columns=_return_production ? zeros(Int32,nlev,nnode+length(orphan_points)) : nothing
    for i in 1:nnode,l in 1:nlev
        row=_extrude_node!(lookup,coords,grid[i,l])
        columns===nothing || (columns[l,i]=row)
    end
    # Native extrudeMesh registers owned Curve orphans only on interior
    # levels. They create Surface orphans, never additional boundary Lines or
    # a fabricated chapeau vertex. Literal Line cells still own the strip.
    for (i,point) in enumerate(_register_orphans ? orphan_points : NTuple{3,Float64}[]),l in 2:nlev-1
        row=_extrude_node!(lookup,coords,_extrude_level_at(spec,point,us[l],transforms,l))
        columns===nothing || (columns[l,nnode+i]=row)
    end
    mesh=_extrude_surface_part(tris,quads,coords,lookup)
    return _return_production ? _ExtrudeMeshProduction(mesh,lookup,vcat(source_tags,orphan_tags),columns) : mesh
end

# `copyMesh(GFace *from, GFace *to, ...)` — the top surface duplicates the
# source's mesh verbatim: every node is `pos.find(Extrude(u_top, p))`, the
# transform evaluated at the last layer level, which welds bitwise onto the
# chapeau curves' nodes (themselves transform evaluations).
function _extrude_top_mesh(m::GeoModel,t::Int,params::_GeoExtrudeParams,
                           spec,src_signed::Int,caller::AbstractString;
                           min_angle_deg::Real=25.0,
                           max_periodic_passes=8,
                           size_field::Union{Nothing,AbstractSizeField}=nothing,
                           source_mesh=nothing,_mesh_input=nothing,_return_production=false)
    src=abs(src_signed)
    srcmesh=source_mesh===nothing ?
        mesh_model_surface(m,src;min_angle_deg=min_angle_deg,
            max_periodic_passes=max_periodic_passes,size_field=size_field,
            _mesh_input,_boundary_regrade=false) :
        source_mesh
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
        mesh=Mesh(out;tris=nodes)
        source_tags=srcmesh isa MixedMesh && srcmesh.entity_data!==nothing ?
            srcmesh.entity_data.external_node_tags : UInt64[]
        return _return_production ? _ExtrudeMeshProduction(mesh,lookup,source_tags,reshape(remap,1,:)) : mesh
    end
    if _extrude_quadtri_top(m,t,params,caller)
        return _extrude_quadtri_top_part(m,src,srcmesh,coords,remap,caller;
                                        quad_to_tri=params.quad_to_tri)
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
    mesh=MixedMesh(out,blocks)
    source_tags=srcmesh.entity_data===nothing ? UInt64[] : srcmesh.entity_data.external_node_tags
    return _return_production ? _ExtrudeMeshProduction(mesh,lookup,source_tags,reshape(remap,1,:)) : mesh
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
    source_node_tags::Vector{UInt64}
    source_used::BitVector
    prism6::Vector{NTuple{6,NTuple{3,Float64}}}
    tets::Vector{NTuple{4,NTuple{3,Float64}}}
    hexes::Vector{NTuple{8,NTuple{3,Float64}}}
    prisms::Vector{NTuple{6,NTuple{3,Float64}}}
    pyramids::Vector{NTuple{5,NTuple{3,Float64}}}
    interior::Vector{NTuple{3,Float64}}
    degenerate::Vector{Bool}
    cell_events::Vector{Tuple{Int,Int}}
end

_ExtrudeVolumeSweep(tag,recombine,cols,prism6,tets,hexes,prisms,pyramids,interior,degenerate)=
    _ExtrudeVolumeSweep(tag,recombine,cols,UInt64[],trues(size(cols,2)),prism6,tets,hexes,prisms,pyramids,interior,degenerate,Tuple{Int,Int}[])

# `extrudeMesh(GFace *from, GRegion *to, ...)` — sweep the source surface's
# cells into per-level records. Recombined volumes emit through
# `createPriPyrTet`/`createHexPri` directly; non-recombined sweeps keep the
# 6-vertex prism records for the global diagonal pass upstream runs in
# `SubdivideExtrudedMesh` (its initial `createPriPyrTet` output is deleted
# before `phase3`, so nothing is emitted here).
function _extrude_volume_sweep(m::GeoModel,t::Int,params::_GeoExtrudeParams,
                               spec,src_signed::Int,caller::AbstractString;_mesh_input=nothing,
                               _record_cell_events=false)
    src=abs(src_signed)
    srcmesh=mesh_model_surface(m,src;_mesh_input,_boundary_regrade=false)
    us=_extrude_level_us(params)
    nlev=length(us)
    cols=_extrude_volume_columns(m,t,src,srcmesh,params,spec,us,caller)
    source_tags=srcmesh isa MixedMesh && srcmesh.entity_data!==nothing ?
        srcmesh.entity_data.external_node_tags : UInt64[]
    source_used=falses(nnodes(srcmesh))
    if srcmesh isa Mesh
        source_used[srcmesh.tris].=true
    else
        for block in srcmesh.blocks
            block isa ElementBlock && msh_dimension(block.msh)==2 || continue
            source_used[block.nodes].=true
        end
    end
    sweep=_ExtrudeVolumeSweep(t,params.recombine,cols,source_tags,source_used,
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false],Tuple{Int,Int}[])
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
            before=(length(sweep.tets),length(sweep.pyramids),length(sweep.prisms))
            _extrude_pripyrtet!(sweep.tets,sweep.pyramids,sweep.prisms,v)
            for (msh,cells,count) in ((4,sweep.tets,before[1]),(7,sweep.pyramids,before[2]),
                                      (6,sweep.prisms,before[3]))
                !_record_cell_events || length(cells)==count || push!(sweep.cell_events,(msh,length(cells)))
            end
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
            before=(length(sweep.hexes),length(sweep.prisms))
            _extrude_hexpri!(sweep.hexes,sweep.prisms,sweep.degenerate,v)
            !_record_cell_events || length(sweep.hexes)==before[1] || push!(sweep.cell_events,(5,length(sweep.hexes)))
            !_record_cell_events || length(sweep.prisms)==before[2] || push!(sweep.cell_events,(6,length(sweep.prisms)))
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
                              caller::AbstractString;_return_column_rows::Bool=false,_return_production::Bool=false)
    sweep.degenerate[] && throw(ArgumentError(
        "$caller: Wrong hexahedron in extrusion"))
    events=!_return_production ? Tuple{Int,Int}[] : sweep.recombine ? sweep.cell_events :
        Tuple{Int,Int}[(4,column) for column in eachindex(sweep.tets)]
    lookup=Dict{NTuple{3,Float64},Int32}()
    coords=NTuple{3,Float64}[]
    cols=sweep.cols
    column_rows=(_return_column_rows || _return_production) ? zeros(Int32,size(cols)) : nothing
    for i in axes(cols,2),l in axes(cols,1)
        !sweep.source_used[i] && l in (firstindex(cols,1),lastindex(cols,1)) && continue
        row=_extrude_node!(lookup,coords,cols[l,i])
        column_rows===nothing || (column_rows[l,i]=row)
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
        mesh=Mesh(out;tets=nodes)
        return _return_production ? _ExtrudeMeshProduction(mesh,lookup,sweep.source_node_tags,column_rows,events) :
            _return_column_rows ? (mesh=mesh,column_rows=column_rows) : mesh
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
    mesh=MixedMesh(out,blocks)
    return _return_production ? _ExtrudeMeshProduction(mesh,lookup,sweep.source_node_tags,column_rows,events) :
        _return_column_rows ? (mesh=mesh,column_rows=column_rows) : mesh
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
                               edges=nothing,scope=nothing,_mesh_input=nothing,_return_production=false)
    planned=_extrude_nonew_surface_mesh(m,t,caller;scope=scope,
        min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
        size_field=size_field)
    planned===nothing || return planned
    entry=_extrude_entity_params(m,2,t,caller)
    entry===nothing && return nothing
    params,spec,link=entry
    if spec isa NamedTuple && get(spec,:type,nothing)===:boundary_layer
        link[1]==1 || throw(ArgumentError(
            "$caller: Surface[$t] boundary-layer extrusion needs a curve " *
            "generatrix"))
        return _extrude_bl_lateral_mesh(m,t,params,spec,link[2],caller;
                                        _mesh_input)
    end
    if link[1]==1
        return _extrude_lateral_mesh(m,t,params,spec,link[2],caller;
                                   edges=edges,_mesh_input,
                                   _return_production=_return_production && params.quad_to_tri===:none)
    end
    return _extrude_top_mesh(m,t,params,spec,link[2],caller;
        min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
        size_field=size_field,_mesh_input,
        _return_production=_return_production && params.quad_to_tri===:none)
end

# `meshGRegionExtruded` + the region-local half of `SubdivideExtrudedMesh` —
# the standalone `mesh_model_volume` arm runs the subdivision over this one
# volume's prism records (identical to the global result when no second
# extruded volume shares its lateral faces).
function _extrude_curve_params_stage(m::GeoModel)
    # Lazy curve grading and surface boundary writeback replace dictionary
    # entries; they do not mutate the stored parameter vectors or geometry.
    # Keep those writes local without copying the size callback or geometry.
    working=GeoModel(ntuple(i->getfield(m,i),Val(fieldcount(GeoModel)))...)
    working.curve_params=copy(m.curve_params)
    return working
end

function _extrude_volume_mesh(m::GeoModel,t::Int,caller::AbstractString;
                              scope=nothing,
                              size_field::Union{Nothing,AbstractSizeField}=nothing,
                              _working_model::Bool=false,_mesh_input=nothing)
    planned=_extrude_nonew_volume_mesh(m,t,caller;scope=scope,size_field=size_field)
    planned===nothing || return planned
    entry=_extrude_entity_params(m,3,t,caller)
    entry===nothing && return nothing
    params,spec,link=entry
    staged=params.recombine && params.quad_to_tri===:add_verts && !_working_model
    working=staged ? _extrude_curve_params_stage(m) : m
    sweep=_extrude_volume_sweep(working,t,params,spec,link[2],caller;_mesh_input)
    if params.recombine && params.quad_to_tri===:add_verts
        _extrude_quadtri_addverts!(working,sweep,caller)
    elseif !params.recombine
        _extrude_subdivide!(_ExtrudeVolumeSweep[sweep],caller)
    end
    part=_extrude_volume_part(working,sweep,caller)
    staged && (m.curve_params=working.curve_params)
    return part
end

# The model-wide extrusion pass upstream interleaves into `Mesh 3`:
# `meshGRegionExtruded` for every swept volume first, then the shared
# `SubdivideExtrudedMesh` diagonal selection, then lateral-surface remesh
# against the global edge set. Returns the per-volume parts plus the
# remesh bookkeeping `_geo_mesh_model` applies to the surface parts.
function _extrude_volume_pass(m::GeoModel,caller::AbstractString;
                             scope=nothing,_mesh_input=nothing,volume_tags=keys(m.volumes))
    tags=Int[]
    for tag in volume_tags
        _extrude_gate(m,3,tag)===nothing || push!(tags,tag)
    end
    isempty(tags) && return nothing
    sort!(tags)
    scope===nothing && (scope=_extrude_nonew_scope(m,caller;tags=tags,_mesh_input))
    sweeps=_ExtrudeVolumeSweep[]
    parts=Dict{Int,Union{Mesh,MixedMesh}}()
    productions=Dict{Int,_ExtrudeMeshProduction}()
    temporary_elements=Dict{Int,Int}()
    for tag in tags
        params,spec,link=_extrude_entity_params(m,3,tag,caller)
        if _extrude_is_nonew(params)
            parts[tag]=_extrude_nonew_volume_mesh(m,tag,caller;scope=scope)
            continue
        end
        push!(sweeps,_extrude_volume_sweep(m,tag,params,spec,link[2],caller;_mesh_input,
            _record_cell_events=_mesh_input!==nothing && params.quad_to_tri===:none))
    end
    subdivided=[s for s in sweeps if !s.recombine]
    for sweep in subdivided
        # meshGRegionExtruded constructs these initial cells before the
        # shared subdivision deletes them. Their labels remain consumed.
        temporary_elements[sweep.tag]=length(sweep.prism6)
    end
    edges=nothing
    isempty(subdivided) || (edges=_extrude_subdivide!(subdivided,caller))
    for sweep in sweeps
        params=_extrude_gate(m,3,sweep.tag)
        params.recombine && params.quad_to_tri===:add_verts || continue
        _extrude_quadtri_addverts!(m,sweep,caller;edges=edges)
    end
    for sweep in sweeps
        params=_extrude_gate(m,3,sweep.tag)
        production=_extrude_volume_part(m,sweep,caller;
            _return_production=_mesh_input!==nothing && params.quad_to_tri===:none)
        if production isa _ExtrudeMeshProduction
            parts[sweep.tag]=production.mesh;productions[sweep.tag]=production
        else
            parts[sweep.tag]=production
        end
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
    return (parts=parts,productions=productions,temporary_elements=temporary_elements,
            edges=edges,remesh=remesh,scope=scope)
end

include("ModelExtrudeQuadTri.jl")
include("ModelExtrudeQuadTriNoNew.jl")
include("ModelExtrudeNoNewScope.jl")
include("ModelExtrudeNoNewRetainedPlan.jl")
include("ModelExtrudeNoNewRetainedBoundary.jl")
include("ModelExtrudeNoNewProjection.jl")
