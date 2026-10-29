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

# Position-keyed row lookup shared by the boundary-layer arms of
# `_extrude_at`/`_extrude_level_at`: exact hit first, then the nearest row
# within a scale-aware tolerance — mirror-evaluated keys can drift an ulp
# off the stored positions, and a later native re-grading shifts the same
# logical node a few e-9 along its curve. The row correspondence stays
# unambiguous — the drift sits far inside the key's Voronoi cell — so the
# nearest row is accepted when it dominates the runner-up by an order of
# magnitude and the displacement stays well under any reachable layer
# offset. A genuinely foreign point sits between keys with best≈second
# distances and still rejects.
function _extrude_bl_row_lookup(rows::AbstractDict{NTuple{3,Float64},V},
                                p::NTuple{3,Float64},
                                caller::AbstractString) where V
    top=get(rows,p,nothing)
    top===nothing || return top
    isempty(rows) && throw(ArgumentError(
        "$caller: boundary-layer copy has no stored top positions"))
    scale=1.0
    for pos in keys(rows)
        scale=max(scale,abs(pos[1]),abs(pos[2]),abs(pos[3]))
    end
    best=nothing;bd=Inf;bd2=Inf
    for (pos,v) in rows
        d=hypot(p[1]-pos[1],p[2]-pos[2],p[3]-pos[3])
        if d<bd
            bd2=bd;bd=d;best=v
        elseif d<bd2
            bd2=d
        end
    end
    bd<=1e-9*scale && return best
    bd<=1e-7*scale && bd<=0.1*bd2 && return best
    throw(ArgumentError(
        "$caller: boundary-layer copy evaluated $(p) off the wall's " *
        "stored nodes"))
end

# `ExtrudeParams::Extrude(t)` — the rigid sweep transform at the normalized
# parameter: `x += t·T` for `:translate`, `ProtudeXYZ` (rotate about the
# axis point by `angle·t`) for `:rotate`, and `ProtudeXYZ` followed by
# `x += t·T` for `:translate_rotate`.
function _extrude_at(spec,p::NTuple{3,Float64},t::Float64,
                     caller::AbstractString)
    if spec.type===:boundary_layer
        # The chapeau of a `BOUNDARY_LAYER` extrusion has no rigid transform:
        # each generatrix node follows its own field direction to the copy's
        # top row. `spec.tops` maps the evaluated source position to that
        # exact emitted top node.
        return _extrude_bl_row_lookup(spec.tops,p,caller)
    elseif spec.type===:translate
        # A boundary-layer ScaleLast connector remaps the last layer's `u`
        # interval through its source vertex's accumulated scale
        # (`ExtrudeParams::Extrude` — `t = (t-h[N-2])·scale + h[N-2]`).
        if haskey(spec,:scale) && spec.scale!==nothing
            hprev,by=spec.scale
            t>hprev && (t=hprev+(t-hprev)*by)
        end
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
    if spec.type===:boundary_layer && haskey(spec,:levels)
        # A resolved `bl_surface` spec stores each source position's whole
        # displaced column (`x + t·nn` per level, ScaleLast already baked
        # into the last row). `u` is the absolute level height — callers pass
        # the exact stored `us` entries, so the row index resolves by the
        # first `uvals` hit; `u <= 0` is the source level itself.
        u<=0 && return p
        row=_extrude_bl_row_lookup(spec.levels,p,
            "boundary-layer level evaluation")
        k=searchsortedfirst(spec.uvals,u)
        return row[min(k,length(row))]
    end
    if spec.type===:boundary_layer
        # Curve-arm `boundary_layer` records carry only the top row;
        # intermediate-level queries land on it like `_extrude_at`.
        return u<=0 ? p : _extrude_at(spec,p,u,
            "boundary-layer level evaluation")
    end
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
function _extrude_curve_nodes(m::GeoModel,curve::Int,caller::AbstractString;
                              _mesh1d_options=nothing)
    haskey(m.curves,curve) || throw(ArgumentError(
        "$caller: unknown Curve[$curve]"))
    params=get(m.curve_params,curve,nothing)
    if params===nothing
        _model_surface_mesh_curves!(m,Int[curve],caller;
            _mesh1d_options=_mesh1d_options)
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

function _extrude_curve_source(m,curve,input,caller;_mesh1d_options=nothing)
    input===nothing && return (
        _extrude_curve_nodes(m,curve,caller;_mesh1d_options=_mesh1d_options),
        UInt64[],NTuple{3,Float64}[],UInt64[])
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

# One boundary-layer wall record: the generatrix's stored node chain, whether
# the source curve is closed, its plane-aligned growth normal, and the
# statement-level `mesh.BoundaryLayerIndex`/`mesh.ViewIndex`/`mesh.ScaleLast`
# the entity records carry. `surf` is the signed source face of a 3-D
# boundary-layer strip (0 for the ordinary 2-D curve arm): the field then
# comes from the face's elements, and the strip kernel cannot rely on a
# planar frame.
const _ExtrudeBlWall=NamedTuple{(:chain,:closed,:peff,:index,:view,:scale_last,:surf),
    Tuple{Vector{NTuple{3,Float64}},Bool,NTuple{3,Float64},Int,Int,Bool,Int}}

# `smooth_data::add(x,y,z,3,nn)` — every wall segment's left normal
# (`cross(peff, t)`: `t × ẑ` upstream, generalized to each wall's own plane)
# accumulates on both endpoint positions; the merged per-position mean
# normalizes into the growth direction. Positions whose incident segments
# cancel (opposing walls) stay unset and keep the kernel's own average.
#
# Upstream keeps two independent `normals[index]` fields selected by each
# copy's `mesh.BoundaryLayerIndex` (`Using Index[i]`): walls under different
# indices never share a column, and `normalsCoherence` only bridges the two
# fields at explicitly merged points. Walls whose index falls outside 0..1
# contribute nothing (upstream errors and skips them), leaving their strip on
# the kernel's own vertex normals.
#
# ScaleLast walls additionally accumulate each source element's average edge
# length per vertex (`smooth_data::add_scale`); non-ScaleLast walls at an
# index where any wall is ScaleLast contribute 1.0 instead
# (`skipScaleCalc`), so a shared vertex "averages" between scaled and
# unscaled. The mean is upstream's running-average update.
# `samplers` maps a nonnegative `ViewIndex` to the resolved
# `OctreePost`-equivalent vector sampler (`PostViewField` restricted to the
# view's vector records, `use_closest=false`) or `nothing` when that view
# carries no vector data — `addExtrudeNormals`'s `octree && !gouraud` gate.
#
# The source-face mesh `addExtrudeNormals` accumulates over: with a
# `_mesh_input`, the carried dim-2 product (the pipeline's produced face
# mesh); otherwise the native producer's own output — requested as a
# production so a plane face takes the same `_gmsh_native_initial_surface`
# path the dim-2 phase's `_return_surface_production` call uses instead of
# the lazy CDT fallback. `production.mesh` and a retained product are both
# returned as the plain mesh.
function _extrude_bl_srcmesh(m::GeoModel,stag::Int;
                             size_field=nothing,_mesh_input=nothing,
                             _mesh1d_options=nothing)
    mesh=mesh_model_surface(m,stag;
        size_field=size_field,_mesh_input=_mesh_input,
        _boundary_regrade=false,_mesh1d_options=_mesh1d_options,
        _return_surface_production=true)
    return mesh isa Mesh || mesh isa MixedMesh ? mesh : mesh.mesh
end

# `_model` enables the 3-D arm: pending/resolved `bl_surface` records
# (COPIED_ENTITY chapeaux of `Extrude{Surface{...};Layers{...}}`) name their
# source faces, whose mesh elements contribute the way upstream's
# `addExtrudeNormals(sourceFaces, sourceFaceInfo, faceSkipScaleCalc)` does —
# and when any source face exists upstream skips the edge arm entirely.
# `_mesh_input` reuses the produced face mesh the volume sweep consumes, so
# field keys land bitwise on the sweep's source-node columns.
function _extrude_bl_normal_field(walls,samplers::AbstractDict=Dict{Int,Any}();
                                  _model=nothing,_mesh_input=nothing,
                                  _mesh1d_options=nothing,_size_field=nothing)
    sums=Dict{Int,Dict{NTuple{3,Float64},NTuple{3,Float64}}}()
    counts=Dict{Int,Dict{NTuple{3,Float64},Int}}()
    scale_sums=Dict{Int,Dict{NTuple{3,Float64},Float64}}()
    scale_counts=Dict{Int,Dict{NTuple{3,Float64},Int}}()
    # `calcLayerScaleFactor[index]` — a global gate upstream sets when any
    # copy at that index is ScaleLast; only then does anything add_scale.
    # Face sources arm it through their chapeau's `mesh.ScaleLast` the same
    # way edge walls do.
    scale_gate=Dict{Int,Bool}()
    for wall in walls
        0<=wall.index<=1 || continue
        wall.scale_last && (scale_gate[wall.index]=true)
    end
    faces=_model===nothing ? nothing : _extrude_bl_face_infos(_model)
    if faces!==nothing && !isempty(faces.infos)
        for (stag,infos) in faces.infos
            for (_,index,_) in infos
                faces.scale_last[stag] &&
                    0<=index<=1 && (scale_gate[index]=true)
            end
        end
    end
    sf_mod=getfield(parentmodule(@__MODULE__),:SizeField)
    if faces!==nothing && !isempty(faces.infos)
        # `sourceFaces` non-empty: `addExtrudeNormals` iterates the source
        # faces only — every source *edge* contribution is skipped, even for
        # sibling curve walls. Each `(invert,index,view)` info entry is one
        # boundary-layer copy's contribution set.
        for (stag,infos) in faces.infos
            # `_return_surface_production` selects the native producer the
            # dim-2 phase uses — `addExtrudeNormals` accumulates over the
            # face's *produced* vertices, and the returned production's mesh
            # is the discretization the consumers' `_mesh_input` products
            # carry (a lazy CDT would miss interior nodes like the
            # recombined face's center fan).
            srcmesh=_extrude_bl_srcmesh(_model,stag;
                size_field=_size_field,_mesh_input=_mesh_input,
                _mesh1d_options=_mesh1d_options)
            for (invert,index,view) in infos
                0<=index<=1 || continue
                sf=get!(Dict{NTuple{3,Float64},NTuple{3,Float64}},
                        sums,index)
                cf=get!(Dict{NTuple{3,Float64},Int},counts,index)
                gated=get(scale_gate,index,false)
                ss=gated ? get!(Dict{NTuple{3,Float64},Float64},
                                scale_sums,index) : nothing
                sc=gated ? get!(Dict{NTuple{3,Float64},Int},
                                scale_counts,index) : nothing
                skip_scale=!faces.scale_last[stag]
                sampler=view>=0 ? get(samplers,view,nothing) : nothing
                _extrude_bl_face_accumulate!(srcmesh,sf,cf,ss,sc;
                    invert=invert,sampler=sampler,skip_scale=skip_scale,
                    sf_mod=sf_mod)
            end
        end
    else
        for wall in walls
            0<=wall.index<=1 || continue
            chain=wall.chain;closed=wall.closed;peff=wall.peff
            nnode=length(chain)
            nseg=closed ? nnode : nnode-1
            sf=get!(Dict{NTuple{3,Float64},NTuple{3,Float64}},
                    sums,wall.index)
            cf=get!(Dict{NTuple{3,Float64},Int},counts,wall.index)
            gated=get(scale_gate,wall.index,false)
            ss=gated ? get!(Dict{NTuple{3,Float64},Float64},
                            scale_sums,wall.index) : nothing
            sc=gated ? get!(Dict{NTuple{3,Float64},Int},
                            scale_counts,wall.index) : nothing
            sampler=wall.view>=0 ? get(samplers,wall.view,nothing) : nothing
            if sampler!==nothing
                # Post-view direction: upstream dedupes the wall's vertices
                # (`verts` set) and adds each `searchVector` result once —
                # the raw interpolated vector, magnitude included, with no
                # `invert` applied to signed generatrices. A miss contributes
                # the zero vector the C++ call leaves in `nn`.
                seen=Set{NTuple{3,Float64}}()
                for pos in chain
                    pos in seen && continue
                    push!(seen,pos)
                    v=sf_mod.vector_raw(sampler,pos[1],pos[2],pos[3])
                    n=v===nothing ? (0.0,0.0,0.0) : v
                    s=get(sf,pos,(0.0,0.0,0.0))
                    sf[pos]=(s[1]+n[1],s[2]+n[2],s[3]+n[3])
                    cf[pos]=get(cf,pos,0)+1
                end
            else
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
                        s=get(sf,pos,(0.0,0.0,0.0))
                        sf[pos]=(s[1]+n[1],s[2]+n[2],s[3]+n[3])
                        cf[pos]=get(cf,pos,0)+1
                    end
                end
            end
            # `add_scale` stays per element-vertex in both modes — a shared
            # vertex averages its incident segment lengths either way
            # (`aveLength = GetAveEdgeLength(elem)`; `skipScaleCalc` → 1.0).
            ss===nothing && continue
            for i in 1:nseg
                j=i==nnode ? 1 : i+1
                p,q=chain[i],chain[j]
                Lt=hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3])
                Lt>0 || continue
                seg_scale=wall.scale_last ? Lt : 1.0
                for pos in (p,q)
                    ss[pos]=get(ss,pos,0.0)+seg_scale
                    sc[pos]=get(sc,pos,0)+1
                end
            end
        end
    end
    fields=Dict{Int,Dict{NTuple{3,Float64},NTuple{3,Float64}}}()
    for (index,sf) in sums
        cf=counts[index]
        field=Dict{NTuple{3,Float64},NTuple{3,Float64}}()
        for (pos,s) in sf
            # `smooth_data` keeps a running mean per position — kept raw
            # here; unit normalization is `_extrude_bl_view_fields`'s
            # global `normalize` pass so a single vector-view wall leaves
            # every index's field un-normalized like upstream.
            c=cf[pos]
            field[pos]=(s[1]/c,s[2]/c,s[3]/c)
        end
        fields[index]=field
    end
    scales=Dict{Int,Dict{NTuple{3,Float64},Float64}}()
    for (index,ss) in scale_sums
        sc=scale_counts[index]
        field=Dict{NTuple{3,Float64},Float64}()
        for (pos,s) in ss
            field[pos]=s/sc[pos]
        end
        scales[index]=field
    end
    return (fields=fields,scales=scales)
end

# The `(invert,index,view)` info set `Mesh2DWithBoundaryLayers` collects per
# source face from every COPIED_ENTITY boundary-layer record — pending
# `bl_surface` specs plus their resolved `boundary_layer` forms keep the
# same bookkeeping fields (`bl3d` marks both). `scale_last` folds
# `faceSkipScaleCalc[tag]`: the last copy record's flag wins upstream, so
# iterating sorted spec tags reproduces the same per-source flag.
function _extrude_bl_face_infos(m::GeoModel)
    infos=Dict{Int,Set{Tuple{Bool,Int,Int}}}()
    scale_last=Dict{Int,Bool}()
    for (d,tag) in sort!(collect(keys(m.meshing.extrude_specs)))
        d==2 || continue
        spec=m.meshing.extrude_specs[(d,tag)]
        spec isa NamedTuple || continue
        get(spec,:bl3d,false) || continue
        src=Int(get(spec,:src,0))
        src==0 && continue
        stag=abs(src)
        push!(get!(infos,stag,Set{Tuple{Bool,Int,Int}}()),
              (src<0,Int(get(spec,:index,0)),Int(get(spec,:view,-1))))
        params=get(m.meshing.extrude,(2,tag),nothing)
        scale_last[stag]=params===nothing ? false : params.scale_last
    end
    return (infos=infos,scale_last=scale_last)
end

# `addExtrudeNormals` over one source face's triangle/quadrangle elements:
# gouraud arm — unit element normals angle-weighted per corner
# (`getAngleAtVertex·2/π`), or raw normals plus `aveLength` scale
# accumulation when the index's scale gate is armed; vector-view arm —
# deduped `searchVector` samples (`invert` never applies). `invert` follows
# the chapeau's signed `geo.Source`.
function _extrude_bl_face_accumulate!(srcmesh,sf,cf,ss,sc;
                                    invert::Bool,sampler,skip_scale::Bool,
                                    sf_mod)
    scoords=srcmesh.coords
    posat(row)=(scoords[1,row],scoords[2,row],scoords[3,row])
    cells=_extrude_bl_face_cells(srcmesh)
    if sampler!==nothing
        seen=Set{NTuple{3,Float64}}()
        for cell in cells
            if ss!==nothing
                Lt=_extrude_bl_ave_edge_length(map(posat,cell))
                seg_scale=skip_scale ? 1.0 : Lt
                for row in cell
                    pos=posat(row)
                    push!(seen,pos)
                    seg_scale==0.0 && continue
                    ss[pos]=get(ss,pos,0.0)+seg_scale
                    sc[pos]=get(sc,pos,0)+1
                end
            else
                for row in cell
                    push!(seen,posat(row))
                end
            end
        end
        for pos in seen
            v=sf_mod.vector_raw(sampler,pos[1],pos[2],pos[3])
            n=v===nothing ? (0.0,0.0,0.0) : v
            s=get(sf,pos,(0.0,0.0,0.0))
            sf[pos]=(s[1]+n[1],s[2]+n[2],s[3]+n[3])
            cf[pos]=get(cf,pos,0)+1
        end
        return nothing
    end
    gated=ss!==nothing
    for cell in cells
        pts=map(posat,cell)
        n=if length(cell)==3
            _extrude_bl_cell_normal(pts[1],pts[2],pts[3])
        else
            _extrude_bl_quad_normal(pts[1],pts[2],pts[3],pts[4])
        end
        n===nothing && continue
        invert && (n=(-n[1],-n[2],-n[3]))
        if gated
            Lt=skip_scale ? 1.0 : _extrude_bl_ave_edge_length(pts)
            for k in eachindex(cell)
                pos=pts[k]
                s=get(sf,pos,(0.0,0.0,0.0))
                sf[pos]=(s[1]+n[1],s[2]+n[2],s[3]+n[3])
                cf[pos]=get(cf,pos,0)+1
                Lt==0.0 && continue
                ss[pos]=get(ss,pos,0.0)+Lt
                sc[pos]=get(sc,pos,0)+1
            end
        else
            for k in eachindex(cell)
                pos=pts[k]
                fact=_extrude_bl_corner_angle(pts,k)*2/π
                w=fact>0 ? fact : 1.0
                s=get(sf,pos,(0.0,0.0,0.0))
                sf[pos]=(s[1]+n[1]*w,s[2]+n[2]*w,s[3]+n[3]*w)
                cf[pos]=get(cf,pos,0)+1
            end
        end
    end
    return nothing
end

# The triangle/quadrangle node-row tuples a face mesh contributes to the
# `addExtrudeNormals` accumulation — `Mesh` tris plus a `MixedMesh`'s
# dim-2 element blocks.
function _extrude_bl_face_cells(srcmesh)
    cells=Tuple{Vararg{Int}}[]
    if srcmesh isa Mesh
        for i in 1:ntris(srcmesh)
            push!(cells,(srcmesh.tris[1,i],srcmesh.tris[2,i],
                         srcmesh.tris[3,i]))
        end
    else
        for block in srcmesh.blocks
            block isa ElementBlock && msh_dimension(block.msh)==2 ||
                continue
            for i in axes(block.nodes,2)
                block.msh==2 &&
                    push!(cells,(block.nodes[1,i],block.nodes[2,i],
                                 block.nodes[3,i]))
                block.msh==3 &&
                    push!(cells,(block.nodes[1,i],block.nodes[2,i],
                                 block.nodes[3,i],block.nodes[4,i]))
            end
        end
    end
    return cells
end

# `MElement::getFace(0).normal()` — unit normal of a triangle/quadrangle
# element (diagonal cross for the quad). `nothing` when degenerate.
function _extrude_bl_cell_normal(a,b,c)
    ux,uy,uz=b[1]-a[1],b[2]-a[2],b[3]-a[3]
    vx,vy,vz=c[1]-a[1],c[2]-a[2],c[3]-a[3]
    nx,ny,nz=uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx
    L=hypot(nx,ny,nz)
    L>0 || return nothing
    return (nx/L,ny/L,nz/L)
end
_extrude_bl_quad_normal(a,b,c,d)=
    _extrude_bl_cell_normal(a,c,b,d)

# `GetAveEdgeLength` — mean of the element's edge lengths.
function _extrude_bl_ave_edge_length(pts)
    n=length(pts)
    total=0.0
    for i in 1:n
        p,q=pts[i],pts[mod1(i+1,n)]
        total+=hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3])
    end
    return total/n
end

# `MElement::getAngleAtVertex` — the interior angle at corner `k` in
# radians, between the two incident polygon edges.
function _extrude_bl_corner_angle(pts,k::Int)
    n=length(pts)
    p=pts[k]
    q1=pts[mod1(k-1,n)];q2=pts[mod1(k+1,n)]
    ux,uy,uz=q1[1]-p[1],q1[2]-p[2],q1[3]-p[3]
    vx,vy,vz=q2[1]-p[1],q2[2]-p[2],q2[3]-p[3]
    cx,cy,cz=uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx
    return atan(hypot(cx,cy,cz),ux*vx+uy*vy+uz*vz)
end

# Upstream `smooth_data::normalize` — unit-normalize each field entry in
# place; a zero mean drops out like the previous sum-normalized path did.
function _extrude_bl_unitize!(field::Dict{NTuple{3,Float64},NTuple{3,Float64}})
    for (pos,s) in collect(field)
        L=hypot(s[1],s[2],s[3])
        if L>0
            field[pos]=(s[1]/L,s[2]/L,s[3]/L)
        else
            delete!(field,pos)
        end
    end
    return field
end

# `View[-3]` — upstream's `special3dbox` post-normalize snap: every field
# entry's component becomes ±1/0 by a 0.1 deadband, turning merged diagonals
# into axis-aligned growth.
function _extrude_bl_snap_axis!(field::Dict{NTuple{3,Float64},NTuple{3,Float64}})
    for (pos,n) in field
        field[pos]=(n[1]<-0.1 ? -1.0 : n[1]>0.1 ? 1.0 : 0.0,
                    n[2]<-0.1 ? -1.0 : n[2]>0.1 ? 1.0 : 0.0,
                    n[3]<-0.1 ? -1.0 : n[3]>0.1 ? 1.0 : 0.0)
    end
    return field
end

# `View[-5]` — upstream's `extrudeField`: scale every normalized entry by the
# background field's value at its position (no background → unchanged).
function _extrude_bl_scale_background!(field,bg)
    bg===nothing && return field
    evalbg=getfield(parentmodule(@__MODULE__),:SizeField).field_value
    for (pos,n) in field
        r=evalbg(bg,pos[1],pos[2],pos[3])
        (isfinite(r) && r!=0) || continue
        field[pos]=(n[1]*r,n[2]*r,n[3]*r)
    end
    return field
end

# One node position's mean normal direction and incident scale over the
# *live* source-face mesh — the same arms `_extrude_bl_face_accumulate!`
# contributes for every `(invert,index,view)` copy, plus the
# `_extrude_bl_view_fields` post passes (`unitize` covers
# `smooth_data::normalize`; `snap`/`scale_bg` only fire when no resolved
# vector view exists, mirroring the global `normalize` flag).
function _extrude_bl_live_node(srcmesh,pos,cells,normalmap,ctx)
    sampler=ctx.sampler
    sf_mod=getfield(parentmodule(@__MODULE__),:SizeField)
    if sampler!==nothing
        v=sf_mod.vector_raw(sampler,pos[1],pos[2],pos[3])
        return v===nothing ? (0.0,0.0,0.0) : v
    end
    contrib=_extrude_bl_live_contrib(srcmesh,pos,cells,normalmap,ctx)
    c=length(contrib)
    c==0 && return (0.0,0.0,0.0)
    # `smooth_data`'s running mean per position (upstream also divides a
    # vertex's accumulated normals by its contribution count).
    sx=sy=sz=0.0
    for v in contrib
        sx+=v[1];sy+=v[2];sz+=v[3]
    end
    n=(sx/c,sy/c,sz/c)
    if ctx.unitize
        L=hypot(n[1],n[2],n[3])
        n=L>0 ? (n[1]/L,n[2]/L,n[3]/L) : (0.0,0.0,0.0)
        if ctx.snap
            n=(n[1]<-0.1 ? -1.0 : n[1]>0.1 ? 1.0 : 0.0,
               n[2]<-0.1 ? -1.0 : n[2]>0.1 ? 1.0 : 0.0,
               n[3]<-0.1 ? -1.0 : n[3]>0.1 ? 1.0 : 0.0)
        end
        if ctx.scale_bg && ctx.bg!==nothing
            r=sf_mod.field_value(ctx.bg,pos[1],pos[2],pos[3])
            (isfinite(r) && r!=0) && (n=(n[1]*r,n[2]*r,n[3]*r))
        end
    end
    return n
end

# Per-cell normal (or angle-weighted) contributions incident to `pos` —
# the `_extrude_bl_face_accumulate!` gouraud/gated arms restricted to the
# elements touching one vertex. `normalmap` caches element normals.
function _extrude_bl_live_contrib(srcmesh,pos,cells,normalmap,ctx)
    posat(row)=(srcmesh.coords[1,row],srcmesh.coords[2,row],
                srcmesh.coords[3,row])
    out=NTuple{3,Float64}[]
    for (ci,cell) in enumerate(cells)
        hit=0
        for k in eachindex(cell)
            posat(cell[k])==pos && (hit=k;break)
        end
        hit==0 && continue
        n=get!(normalmap,ci) do
            pts=map(posat,cell)
            en=length(cell)==3 ?
                _extrude_bl_cell_normal(pts[1],pts[2],pts[3]) :
                _extrude_bl_quad_normal(pts[1],pts[2],pts[3],pts[4])
            en===nothing ? (0.0,0.0,0.0) :
                ctx.invert ? (-en[1],-en[2],-en[3]) : en
        end
        n==(0.0,0.0,0.0) && continue
        if ctx.gated
            push!(out,n)
        else
            pts=map(posat,cell)
            fact=_extrude_bl_corner_angle(pts,hit)*2/π
            w=fact>0 ? fact : 1.0
            push!(out,(n[1]*w,n[2]*w,n[3]*w))
        end
    end
    return out
end

# `boundaryLayerVertices` covers `face->mesh_vertices` — the *live* face
# mesh — but a resolved `bl_surface` spec's `levels`/`tops` were keyed by
# the mesh `prepare` measured. The detached dim-0/1 phases can resolve the
# spec before the native surface producer's exact discretization exists;
# when the consumer's mesh is a different discretization of the same face,
# extend the row maps for its uncovered positions with the same
# incident-element accumulation, keyed once per mesh object.
function _extrude_bl_ensure_srcmesh!(m::GeoModel,dim::Int,tag::Int,
                                   caller::AbstractString,srcmesh)
    spec=get(m.meshing.extrude_specs,(dim,tag),nothing)
    spec isa NamedTuple || return spec
    ctx=get(spec,:bl_ctx,nothing)
    ctx===nothing && return spec
    oid=objectid(srcmesh)
    oid in ctx.seen && return spec
    push!(ctx.seen,oid)
    levels=spec.levels;tops=spec.tops
    uvals=spec.uvals;nl=length(uvals)
    cells=_extrude_bl_face_cells(srcmesh)
    normalmap=Dict{Int,NTuple{3,Float64}}()
    changed=false
    for i in 1:nnodes(srcmesh)
        pos=(srcmesh.coords[1,i],srcmesh.coords[2,i],srcmesh.coords[3,i])
        # The live mesh is authoritative: a position the prepare-time mesh
        # already covered still gets recomputed — its incident elements (and
        # therefore its angle-weighted normal and `GetAveEdgeLength` scale)
        # belong to the discretization the consumers actually sweep.
        n=_extrude_bl_live_node(srcmesh,pos,cells,normalmap,ctx)
        row=n==(0.0,0.0,0.0) ? fill(pos,nl) :
            [(pos[1]+t*n[1],pos[2]+t*n[2],pos[3]+t*n[3]) for t in uvals]
        if ctx.scale_active && n!=(0.0,0.0,0.0)
            s=_extrude_bl_live_scale(srcmesh,pos,cells,ctx)
            if s!=1.0
                for l in eachindex(uvals)
                    t=uvals[l]
                    t>ctx.hprev || continue
                    t2=ctx.hprev+(t-ctx.hprev)*s
                    row[l]=(pos[1]+t2*n[1],pos[2]+t2*n[2],
                            pos[3]+t2*n[3])
                end
            end
        end
        levels[pos]=row;tops[pos]=row[end];changed=true
    end
    changed || return spec
    return spec
end

# `GetAveEdgeLength` mean over the elements incident to `pos` — the
# `add_scale` accumulation for a vertex `prepare`'s mesh never covered.
function _extrude_bl_live_scale(srcmesh,pos,cells,ctx)
    posat(row)=(srcmesh.coords[1,row],srcmesh.coords[2,row],
                srcmesh.coords[3,row])
    total=0.0;count=0
    for cell in cells
        any(row->posat(row)==pos,cell) || continue
        pts=map(posat,cell)
        Lt=ctx.skip_scale ? 1.0 : _extrude_bl_ave_edge_length(pts)
        Lt==0.0 && continue
        total+=Lt;count+=1
    end
    return count==0 ? 1.0 : total/count
end

# `(chain, closed, peff, index, view, scale_last)` for every boundary-layer
# surface's generatrix — the mesh-time equivalent of the exec-time entity
# list, so the field merges across `Extrude` statements just as upstream
# accumulates every BoundaryLayerCurve's edges.
function _extrude_bl_walls(m::GeoModel,caller::AbstractString;
                           _mesh_input=nothing,_mesh1d_options=nothing)
    walls=_ExtrudeBlWall[]
    for (d,tag) in sort!(collect(keys(m.meshing.extrude_specs)))
        d==2 || continue
        spec=m.meshing.extrude_specs[(d,tag)]
        spec isa NamedTuple || continue
        get(spec,:type,nothing)===:boundary_layer || continue
        link=get(m.meshing.extrude_sources,(2,tag),nothing)
        link===nothing && continue
        chain,_,orphans,_=_extrude_curve_source(m,abs(link[2]),_mesh_input,
            caller;_mesh1d_options=_mesh1d_options)
        isempty(orphans) || throw(ArgumentError(
            "$caller: Surface[$tag] boundary-layer generatrix " *
            "Curve[$(abs(link[2]))] has unreferenced mesh nodes"))
        a,b=m.curves[abs(link[2])]
        closed=a==b
        peff=_extrude_bl_plane_normal(chain,closed,link[2]>0,caller)
        params=get(m.meshing.extrude,(2,tag),nothing)
        push!(walls,(chain=chain,closed=closed,peff=peff,
            index=Int(get(spec,:index,0)),view=Int(get(spec,:view,-1)),
            scale_last=params===nothing ? false : params.scale_last,
            surf=Int(get(spec,:surf,0))))
    end
    return walls
end

# Upstream's `addExtrudeNormals` view semantics are global over the whole
# pass: `special3dbox`/`extrudeField` snap/scale every built field once any
# wall selects them. Returns `(snap_axis, scale_background)` flags.
function _extrude_bl_view_flags(walls)
    snap=false;scale=false
    for wall in walls
        wall.view==-3 && (snap=true)
        wall.view==-5 && (scale=true)
    end
    return (snap,scale)
end

# One wall's strip emission: explicit cumulative offsets plus the merged-field
# direction per shared wall vertex. `scale_last` carries the wall's per-vertex
# last-layer scale (`Dict` node-index→factor) for upstream's `ScaleLastLayer`
# rescale of the last `u` interval.
function _extrude_bl_kernel(chain,closed::Bool,peff,field,offsets,caller;
                            scale_last=nothing,scale_hprev=0.0,
                            normalize_dirs::Bool=true)
    poly=_extrude_bl_polyline(chain,closed)
    vd=Dict{Int,NTuple{3,Float64}}()
    for (i,pos) in enumerate(chain)
        d=get(field,pos,nothing)
        d===nothing || (vd[i]=d)
    end
    # `BoundaryLayer` is a sibling submodule included after `Model` — resolve
    # it lazily like `.geo`'s `SizeField` lookup.
    return getfield(parentmodule(@__MODULE__),:BoundaryLayer)._bl_2d_offsets(
        poly,offsets;plane_normal=peff,vertex_dirs=vd,caller=caller,
        vertex_scale_last=scale_last,scale_hprev=scale_hprev,
        normalize_dirs=normalize_dirs)
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
# keys the kernel's exact top-row node. The ordered `(usrc, top)` row backs
# the same discrete-edge eval upstream's `boundaryLayerEdge` runs — a query
# between stored nodes interpolates the top polyline instead of rejecting
# an off-key position.
function _extrude_bl_tops_dict(m::GeoModel,gen::Int,chapeau::Int,
                               signed_c::Int,top_of::AbstractDict,
                               caller::AbstractString)
    params=get(m.meshing.extrude,(1,chapeau),nothing)
    params===nothing && return (Dict{NTuple{3,Float64},NTuple{3,Float64}}(),
                                Tuple{Float64,NTuple{3,Float64}}[])
    vals=_model_curve_extrude_params(m,chapeau,params,
        _model_default_mesh1d_options(m,caller),caller)
    vals===:pending && throw(ArgumentError(
        "$caller: Curve[$chapeau] boundary-layer source has no stored " *
        "mesh parameters"))
    _,t1=_model_curve_param_bounds(m,chapeau,caller)
    u_min=_model_curve_param_bounds(m,gen,caller)[1]
    tops=Dict{NTuple{3,Float64},NTuple{3,Float64}}()
    row=Tuple{Float64,NTuple{3,Float64}}[]
    for v in vals
        usrc=signed_c<0 ? (t1+u_min)-v : v
        p=_periodic_curve_point(m,gen,usrc,caller)
        top=_extrude_bl_top_at(top_of,p,caller)
        tops[p]=top
        push!(row,(usrc,top))
    end
    sort!(row;by=first)
    # A closed source can store the seam parameter twice; duplicate
    # parameters describe the same top vertex.
    keep=ones(Bool,length(row))
    for i in 2:length(row)
        row[i][1]==row[i-1][1] && (keep[i]=false)
    end
    all(keep) || (row=row[keep])
    return (tops,row)
end

# Continuous `boundaryLayerEdge` eval: bitwise key hits return the stored
# top node; a re-graded query drifting an ulp off the stored key still binds
# by node identity — upstream's `pos.find` answers with the Voronoi cell's
# stored vertex — and only a genuinely interior probe interpolates the
# ordered top row. A spec without row data (built before this field
# existed) keeps the strict position-lookup behavior.
function _extrude_bl_row_eval(spec,p,usrc::Float64,caller::AbstractString)
    tops=get(spec,:tops,nothing)
    if tops!==nothing
        top=get(tops,p,nothing)
        top===nothing || return top
        if !isempty(tops)
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
            bd<=1e-7*scale && bd<=0.1*bd2 && return best
        end
    end
    row=get(spec,:top_row,nothing)
    if row!==nothing && !isempty(row)
        length(row)==1 && return row[1][2]
        usrc<=row[1][1] && return row[1][2]
        usrc>=row[end][1] && return row[end][2]
        i=searchsortedlast(row,usrc;by=first)
        i=min(max(i,1),length(row)-1)
        u0,t0=row[i];u1,t1=row[i+1]
        u1==u0 && return t0
        f=(usrc-u0)/(u1-u0)
        abs(f)<=1e-9 && return t0
        abs(1-f)<=1e-9 && return t1
        return (t0[1]+f*(t1[1]-t0[1]),
                t0[2]+f*(t1[2]-t0[2]),
                t0[3]+f*(t1[3]-t0[3]))
    end
    return _extrude_at(spec,p,1.0,caller)
end

# Write the kernel's exact directions/positions back onto the sibling
# records: connector curves keep `:translate` specs whose `T` is the unit
# direction (the stored level params are absolute heights, so `p + h·T`
# lands bitwise on the strip's rows), the chapeau keeps the position→top map
# its `_extrude_at` branch reads, and every generated vertex takes the
# kernel's top-row coordinate. A ScaleLast connector additionally carries
# `scale=(h_prev, s)` — the `Extrude(iLayer=last)` rescale evaluated at the
# connector's own source vertex.
function _extrude_bl_refresh!(m::GeoModel,spec,gen::Int,signed_c::Int,
                              closed::Bool,kernel,
                              caller::AbstractString;
                              _params_model::GeoModel=m)
    a,b=m.curves[gen]
    sbeg,send=signed_c>0 ? (a,b) : (b,a)
    beg_pos=m.points[sbeg];end_pos=m.points[send]
    cbeg,cend=spec.connectors
    scale_last=get(spec,:scale_last,nothing)
    hprev=get(spec,:scale_hprev,0.0)
    # Sibling walls share coincident connectors — the operation's coherence
    # pass keeps one tag; refresh only rewrites survivors.
    if haskey(m.curves,cbeg)
        m.points[m.curves[cbeg][2]]=kernel.tops[beg_pos]
        s=scale_last===nothing ? 1.0 : get(scale_last,beg_pos,1.0)
        m.meshing.extrude_specs[(1,cbeg)]=
            (type=:translate,T=kernel.dirs[beg_pos],
             scale=s==1.0 ? nothing : (hprev,s))
    end
    if cend>0 && haskey(m.curves,cend)
        m.points[m.curves[cend][2]]=kernel.tops[end_pos]
        s=scale_last===nothing ? 1.0 : get(scale_last,end_pos,1.0)
        m.meshing.extrude_specs[(1,cend)]=
            (type=:translate,T=kernel.dirs[end_pos],
             scale=s==1.0 ? nothing : (hprev,s))
    end
    chapeau=spec.chapeau
    if haskey(m.curves,chapeau)
        tops,row=_extrude_bl_tops_dict(_params_model,gen,chapeau,signed_c,
                                       kernel.tops,caller)
        m.meshing.extrude_specs[(1,chapeau)]=
            (type=:boundary_layer,tops=tops,top_row=row)
        m.points[m.curves[chapeau][1]]=kernel.tops[beg_pos]
        closed || (m.points[m.curves[chapeau][2]]=kernel.tops[end_pos])
        # The chapeau's materialized control copies keep `MSH_POINT` type
        # upstream, so `boundaryLayerVertices` never displaces them: index-0
        # copies rewired into their sources at the coherence merge and
        # nonzero-index orphans simply stay coincident with the generatrix.
        # Only a boundary-layer-typed control record (a welded beg/end that
        # also serves as a control) takes the top row.
        for cp in get(m.curve_control_points,chapeau,Int[])
            haskey(m.points,cp) || continue
            _point_bl_class(m,cp)[1] || continue
            m.points[cp]=_extrude_bl_top_at(kernel.tops,m.points[cp],caller;
                                            unbounded=true)
        end
    end
    return nothing
end

# The 3-D counterpart of `_extrude_bl_refresh!`: the source-face field — not
# a planar kernel — produced the strip's `dirs`/`tops`, so the connector
# `:translate` records, the A-curve's `tops`/`top_row` map, and the
# `MSH_POINT_BND_LAYER` vertices are written straight from that field.
# `boundaryLayerVertices` semantics are the same: every boundary-layer copy
# vertex moves to `x + t_top·nn` (with the ScaleLast `get_scale` remap).
function _extrude_bl_surf_refresh!(m::GeoModel,spec,gen::Int,signed_c::Int,
                                   closed::Bool,dirs,tops,scale_field,hprev,
                                   caller::AbstractString;
                                   _params_model::GeoModel=m)
    a,b=m.curves[gen]
    sbeg,send=signed_c>0 ? (a,b) : (b,a)
    beg_pos=m.points[sbeg];end_pos=m.points[send]
    cbeg,cend=spec.connectors
    if haskey(m.curves,cbeg)
        m.points[m.curves[cbeg][2]]=tops[beg_pos]
        s=scale_field===nothing ? 1.0 : get(scale_field,beg_pos,1.0)
        (isfinite(s) && abs(s-1.0)<=1e-12) && (s=1.0)
        m.meshing.extrude_specs[(1,cbeg)]=
            (type=:translate,T=dirs[beg_pos],
             scale=s==1.0 ? nothing : (hprev,s))
    end
    if cend>0 && haskey(m.curves,cend)
        m.points[m.curves[cend][2]]=tops[end_pos]
        s=scale_field===nothing ? 1.0 : get(scale_field,end_pos,1.0)
        (isfinite(s) && abs(s-1.0)<=1e-12) && (s=1.0)
        m.meshing.extrude_specs[(1,cend)]=
            (type=:translate,T=dirs[end_pos],
             scale=s==1.0 ? nothing : (hprev,s))
    end
    chapeau=spec.chapeau
    if haskey(m.curves,chapeau)
        ctops,row=_extrude_bl_tops_dict(_params_model,gen,chapeau,signed_c,
                                       tops,caller)
        m.meshing.extrude_specs[(1,chapeau)]=
            (type=:boundary_layer,tops=ctops,top_row=row)
        m.points[m.curves[chapeau][1]]=tops[beg_pos]
        closed || (m.points[m.curves[chapeau][2]]=tops[end_pos])
        for cp in get(m.curve_control_points,chapeau,Int[])
            haskey(m.points,cp) || continue
            _point_bl_class(m,cp)[1] || continue
            m.points[cp]=_extrude_bl_top_at(tops,m.points[cp],caller;
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
                                  _mesh_input=nothing,_size_field=nothing,
                                  _options=nothing)
    gen=abs(gen_signed)
    chain,_,orphans,_=_extrude_curve_source(m,gen,_mesh_input,caller;
        _mesh1d_options=_options)
    isempty(orphans) || throw(ArgumentError(
        "$caller: Surface[$t] boundary-layer generatrix Curve[$gen] has " *
        "unreferenced mesh nodes"))
    a,b=m.curves[gen]
    closed=a==b
    peff=_extrude_bl_plane_normal(chain,closed,gen_signed>0,caller)
    index=Int(get(spec,:index,0))
    walls=_extrude_bl_walls(m,caller;_mesh_input,_mesh1d_options=_options)
    # `extrudeField` reads `BGM_MeshSize` — the mesh pass's background field,
    # carried on the 1-D options bundle when no explicit `size_field` rides
    # this call.
    background=_size_field===nothing ?
        (_options===nothing ? nothing : getproperty(_options,:inner_field)) :
        _size_field
    fields=_extrude_bl_view_fields(walls,caller;
                                   _background_field=background,
                                   _options=_options,
                                   _views=m.meshing.post_views,
                                   _model=m,_mesh_input=_mesh_input,
                                   _mesh1d_options=_options,
                                   _size_field=_size_field)
    # Resolve any standalone point-generatrix connectors against the same
    # field — upstream displaces them inside this pass.
    _extrude_bl_prepare!(m,caller;_fields=fields,_mesh_input,
        _options=_options,_mesh1d_options=_options,_size_field=_size_field)
    field=get(fields.fields,index,_EMPTY_BL_FIELD)
    offsets=_extrude_level_us(params)[2:end]
    # The `calcLayerScaleFactor[index]` gate: only when some copy at this
    # index is ScaleLast does the last `u` interval rescale per vertex —
    # upstream's `Extrude(iLayer=NbLayer-1)` consults the accumulated scale
    # for every BL entity at the index, ScaleLast or not.
    scale_field=get(fields.scales,index,nothing)
    hprev=length(params.layers)>=2 ? params.heights[end-1] : 0.0
    vertex_scale=nothing
    if scale_field!==nothing
        vertex_scale=Dict{Int,Float64}()
        for (i,pos) in enumerate(chain)
            s=get(scale_field,pos,1.0)
            abs(s-1.0)<=1e-12 && (s=1.0)
            vertex_scale[i]=s
        end
    end
    spec=(;spec...,scale_last=scale_field,scale_hprev=hprev)
    kernel=_extrude_bl_kernel(chain,closed,peff,field,offsets,caller;
                              scale_last=vertex_scale,scale_hprev=hprev,
                              normalize_dirs=fields.normalize)
    _extrude_bl_refresh!(m,spec,gen,gen_signed,closed,kernel,caller)
    return params.recombine ? kernel.mesh : _extrude_bl_unrecombined(kernel.mesh)
end

const _EMPTY_BL_FIELD=Dict{NTuple{3,Float64},NTuple{3,Float64}}()

# Upstream resolves boundary-layer metadata inside `Mesh2DWithBoundaryLayers`:
# an out-of-range `BoundaryLayerIndex` is a non-fatal `Msg::Error` and the
# wall's copy keeps the default field; a `ViewIndex` of -3/-5 selects the
# global special modes, a nonnegative index must name a live post view
# (otherwise `Unknown View[i]: using normals instead`), and any other
# negative value silently keeps plain normals.
function _extrude_bl_check_meta!(walls,options,views=nothing)
    options===nothing && return nothing
    nviews=views===nothing ? 0 : length(views)
    for wall in walls
        wall.index in (0,1) || _model_mesh_error(options,
            "Invalid boundary layer index $(wall.index)")
        v=wall.view
        (v==-1 || v==-3 || v==-5) && continue
        (v>=0 && v<nviews) && continue
        v>=0 && _model_mesh_error(options,
            "Unknown View[$v]: using normals instead")
    end
    return nothing
end

# `_SessionView`/`PosView` → element records for the vector sampler.
# `PosView` exposes `.elements` directly; a session view's flat list data
# decodes through the API writer's record builder. Anything else has no
# vector data and takes the gouraud fallback.
function _extrude_bl_view_records(view,caller)
    if hasproperty(view,:elements)
        return view.elements
    elseif hasproperty(view,:list_types)
        api=getfield(parentmodule(@__MODULE__),:API)
        return api._native_list_view_records(view,caller)
    end
    return nothing
end

# `OctreePost(PView::list[view])` restricted to the view's vector records —
# `getNumVectors() > 0` is the gouraud→octree switch upstream. Returns
# `nothing` when the view carries no vector data.
function _extrude_bl_view_sampler(view,caller)
    records=_extrude_bl_view_records(view,caller)
    records===nothing && return nothing
    vectors=[e for e in records
        if hasproperty(e,:values) && e.values isa AbstractArray &&
           size(e.values,1)==3]
    isempty(vectors) && return nothing
    sf_mod=getfield(parentmodule(@__MODULE__),:SizeField)
    return sf_mod.PostViewField(vectors;use_closest=false)
end

# Per-`ViewIndex` sampler table for the wall set — keyed on the view index
# so sibling walls sharing a `Using View[i]` share one octree. 3-D
# boundary-layer copies carry `Using View[i]` on their `bl_surface` specs
# rather than a wall record, so `face_infos` entries join the same table.
function _extrude_bl_view_samplers(walls,views,caller;face_infos=nothing)
    samplers=Dict{Int,Any}()
    views===nothing && return samplers
    for wall in walls
        v=wall.view
        (v>=0 && v<length(views)) || continue
        haskey(samplers,v) && continue
        samplers[v]=_extrude_bl_view_sampler(views[v+1],caller)
    end
    if face_infos!==nothing
        for (_,infos) in face_infos
            for (_,_,v) in infos
                (v>=0 && v<length(views)) || continue
                haskey(samplers,v) && continue
                samplers[v]=_extrude_bl_view_sampler(views[v+1],caller)
            end
        end
    end
    return samplers
end

# The `bl3d` copy-record analogue of `_extrude_bl_view_flags`: the
# `special3dbox` snap (`View[-3]`) and `extrudeField` background scale
# (`View[-5]`) modes arm on a surface copy's `Using View[...]` the same way
# they arm on a wall's.
function _extrude_bl_face_view_flags(face_infos)
    snap=scale=false
    face_infos===nothing && return (snap=snap,scale=scale)
    for (_,infos) in face_infos
        for (_,_,v) in infos
            v==-3 && (snap=true)
            v==-5 && (scale=true)
        end
    end
    return (snap=snap,scale=scale)
end

# `_extrude_bl_check_meta!` for the `bl3d` records — same recoverable
# diagnostics on `Using Index[...]`/`Using View[...]` metadata a wall
# carries (`addExtrudeNormals` runs its `Msg::Error` switch per copied
# entity the same way it runs it per boundary-layer wall).
function _extrude_bl_check_meta_faces!(face_infos,options,views=nothing)
    (options===nothing || face_infos===nothing) && return nothing
    nviews=views===nothing ? 0 : length(views)
    for (_,infos) in face_infos
        for (_,index,v) in infos
            index in (0,1) || _model_mesh_error(options,
                "Invalid boundary layer index $index")
            (v==-1 || v==-3 || v==-5) && continue
            (v>=0 && v<nviews) && continue
            v>=0 && _model_mesh_error(options,
                "Unknown View[$v]: using normals instead")
        end
    end
    return nothing
end

# Build the per-index smooth-data fields for the current wall set, applying
# upstream's global `special3dbox`/`extrudeField` post-normalize passes.
# `normalize` follows upstream's single flag: one resolved vector view
# (`octree && !gouraud`) skips normalization AND the snap/field post passes
# for every index — the view vectors and gouraud means stay raw, and the
# extrude kernel consumes their magnitudes (`x += n[i]*t`).
function _extrude_bl_view_fields(walls,caller;
                                 _background_field=nothing,_options=nothing,
                                 _views=nothing,_model=nothing,
                                 _mesh_input=nothing,_mesh1d_options=nothing,
                                 _size_field=nothing)
    _extrude_bl_check_meta!(walls,_options,_views)
    face_infos=_model===nothing ? nothing :
        _extrude_bl_face_infos(_model).infos
    _extrude_bl_check_meta_faces!(face_infos,_options,_views)
    samplers=_extrude_bl_view_samplers(walls,_views,caller;
                                       face_infos=face_infos)
    normalize=!any(s->s!==nothing,values(samplers))
    built=_extrude_bl_normal_field(walls,samplers;
        _model=_model,_mesh_input=_mesh_input,
        _mesh1d_options=_mesh1d_options,_size_field=_size_field)
    snap,scale_bg=_extrude_bl_view_flags(walls)
    face_flags=_extrude_bl_face_view_flags(face_infos)
    snap = snap || face_flags.snap
    scale_bg = scale_bg || face_flags.scale
    for (index,field) in built.fields
        if normalize
            _extrude_bl_unitize!(field)
            snap && _extrude_bl_snap_axis!(field)
            scale_bg && _extrude_bl_scale_background!(field,_background_field)
        end
    end
    # The kernel consumes the field as upstream's `Extrude` reads `nn`: raw,
    # never re-normalized. `normalize=false` (`_bl_2d_offsets`'s raw mode) is
    # therefore required not only for vector views but also for the
    # `special3dbox`/`extrudeField` post passes — a snapped (1,-1,0) or
    # background-scaled (0,-2,0) entry loses its magnitude to a unitizing
    # re-pass, and either can leave the wall plane.
    return (;built...,normalize=normalize && !snap && !scale_bg,
            samplers=samplers,bg=_background_field)
end

# Every curve whose stored discretization a pending boundary-layer
# resolution is about to read: the wall generatrices (2-D arms and 3-D
# lateral strips alike — `Mesh2DWithBoundaryLayers` merges the smooth-data
# field globally) plus the boundary curves of each `bl_surface`/lateral
# source face, whose `meshGFace` chains feed `addExtrudeNormals`.
function _extrude_bl_seed_curves(m::GeoModel)
    curves=Set{Int}()
    for ((d,tag),spec) in m.meshing.extrude_specs
        spec isa NamedTuple || continue
        st=get(spec,:type,nothing)
        if d==2 && st===:boundary_layer
            link=get(m.meshing.extrude_sources,(2,tag),nothing)
            link!==nothing && link[1]==1 &&
                push!(curves,abs(link[2]))
            stag=abs(Int(get(spec,:surf,0)))
            stag==0 && continue
        elseif st===:bl_surface
            stag=abs(Int(get(spec,:src,0)))
            stag==0 && continue
        else
            continue
        end
        haskey(m.surfaces,stag) || continue
        for l in m.surfaces[stag],c in m.loops[l]
            push!(curves,abs(c))
        end
    end
    return curves
end

# `Mesh2DWithBoundaryLayers` resolves its displaced-position maps off the
# native 1-D producer's stored grading — the parameters the curve phase
# publishes on the actual products. A lazy grade drifts ~1e-9 from it on
# curved walls (the `_native_producer` lc clamps differ in the last ulp),
# leaving every position-keyed `levels`/`tops` map split into stale keys
# and stale values the welded consumers cannot reconcile bitwise. Re-grade
# the seed set natively on `sm` — the same transport the exec-time stage
# copy already performs for `working` — so the maps land on the final
# coordinates. Entries are always re-graded: an existing lazy entry is
# exactly the drift being removed, and a native one re-grades identically.
# Curves the producer cannot grade yet (`:pending` periodic copies, deferred
# sources) keep the lazy chain, exactly `_extrude_curve_nodes`'s fallback.
function _extrude_bl_seed_native!(sm::GeoModel,curves,
                                  options::_ModelMesh1DOptions,
                                  caller::AbstractString)
    pending=Set{Int}(curves)
    for _ in 1:max(1,length(curves))
        isempty(pending) && break
        progressed=false
        for c in sort!(collect(pending))
            graded=try
                _model_mesh_curve!(sm,c,options,caller;_native_producer=true)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                nothing
            end
            graded===:pending && continue
            if graded isa Vector{Float64}
                sm.curve_params[c]=graded
            elseif graded isa _ModelNativeCurveGrade
                sm.curve_params[c]=Float64.(graded.parameters)
            else
                # `:keep` — the stored discretization stands untouched.
            end
            progressed=true;delete!(pending,c)
        end
        progressed || break
    end
    return nothing
end

# `Mesh2DWithBoundaryLayers`'s `ExtrudePoint` arm: a standalone
# point-generatrix connector (`type=:bl_point` until resolved) sweeps along
# the smooth field's normal at its source position — zero when the position
# has no accumulated entry, which leaves the `MSH_POINT_BND_LAYER` copy
# coincident with its source exactly like upstream's degenerate output. The
# resolved `:translate` spec (plus its ScaleLast remap, when both the copy
# selects ScaleLast and the index's scale gate is armed) is idempotent:
# re-entry finds no pending specs.
# Re-entrancy gate: meshing a source face through the native producer runs a
# detached dim-01 plan whose `_model_point_mesh_parts` re-enters this
# function on its transaction copy — that inner resolution is incidental to
# face meshing (upstream's `(*it)->mesh(true)` never re-enters the BL pass).
# Task-local so detached parallel surface units keep their own gate.
function _extrude_bl_prepare!(m::GeoModel,caller::AbstractString;
                              _fields=nothing,_mesh_input=nothing,
                              _options=nothing,_model=nothing,
                              _mesh1d_options=nothing,_size_field=nothing)
    get(task_local_storage(),:TESS_BL_PREPARE,false) && return nothing
    task_local_storage(:TESS_BL_PREPARE,true)
    try
        _extrude_bl_prepare_inner!(m,caller;_fields,_mesh_input,_options,
            _model,_mesh1d_options,_size_field)
    finally
        delete!(task_local_storage(),:TESS_BL_PREPARE)
    end
    return nothing
end

function _extrude_bl_prepare_inner!(m::GeoModel,caller::AbstractString;
                              _fields=nothing,_mesh_input=nothing,
                              _options=nothing,_model=nothing,
                              _mesh1d_options=nothing,_size_field=nothing)
    pending=Tuple{Int,Any}[]
    pending_surf=Tuple{Int,Int,Any}[]
    pending_lat=Tuple{Int,Any}[]
    for ((d,tag),spec) in m.meshing.extrude_specs
        spec isa NamedTuple || continue
        if d==1 && get(spec,:type,nothing)===:bl_point
            push!(pending,(tag,spec))
        elseif d in (2,3) && get(spec,:type,nothing)===:bl_surface
            push!(pending_surf,(d,tag,spec))
        elseif d==2 && get(spec,:type,nothing)===:boundary_layer &&
            Int(get(spec,:surf,0))!=0 && !haskey(spec,:levels)
            # A 3-D boundary-layer strip: the source-face field resolves its
            # generatrix columns and sibling connector/chapeau records.
            push!(pending_lat,(tag,spec))
        end
    end
    isempty(pending) && isempty(pending_surf) && isempty(pending_lat) &&
        return nothing
    field_model=_model===nothing ? m : _model
    # The chains and source-face meshes these maps key on must be the
    # *native* 1-D grading — the parameters the curve phase will publish —
    # not a lazy grade (whose `_native_producer`-off clamps drift ~1e-9 in
    # the last ulp). Re-grade the wall generatrices and face boundaries
    # before any lazy `_extrude_curve_source`/`_extrude_bl_srcmesh` read can
    # freeze drifted entries into the stored spec maps. `_mesh1d_options`
    # is the caller's own producer options so the seeded grading is the
    # same one the products carry.
    seed_options=_mesh1d_options===nothing ?
        _model_default_mesh1d_options(m,caller) : _mesh1d_options
    seed_curves=_extrude_bl_seed_curves(m)
    if !isempty(seed_curves)
        _extrude_bl_seed_native!(m,seed_curves,seed_options,caller)
        field_model===m ||
            _extrude_bl_seed_native!(field_model,seed_curves,
                                     seed_options,caller)
    end
    if _fields===nothing
        walls=_extrude_bl_walls(m,caller;_mesh_input,
            _mesh1d_options=_mesh1d_options)
        _fields=_extrude_bl_view_fields(walls,caller;_options=_options,
                                        _views=m.meshing.post_views,
                                        _model=field_model,
                                        _mesh_input=_mesh_input,
                                        _mesh1d_options=_mesh1d_options,
                                        _size_field=_size_field)
    end
    # `Mesh2DWithBoundaryLayers`'s `ExtrudeSurface` arm: each COPIED_ENTITY
    # chapeau/volume record resolves into a column field — the displaced
    # position `x + t·nn` for every level `t`, keyed by the source face's
    # mesh-node positions, exactly like `boundaryLayerVertices`. ScaleLast
    # remaps the top level through the index's scale field the way
    # `ExtrudeParams::Extrude(iLayer, iElemLayer, ..., scaleFactor)` does.
    # `bl_ctx` lets a later consumer re-extend the row maps against the
    # *actual* discretization it sweeps when the producer's face mesh
    # differs from the one resolved here (`face->mesh_vertices` upstream).
    if !isempty(pending_surf)
        samplers_ctx=get(_fields,:samplers,nothing)
        rawnorm_ctx=samplers_ctx===nothing ||
            !any(s->s!==nothing,values(samplers_ctx))
        snap_ctx,scale_bg_ctx=_extrude_bl_view_flags(
            _extrude_bl_walls(m,caller;_mesh_input,
                _mesh1d_options=_mesh1d_options))
        face_flags_ctx=_extrude_bl_face_view_flags(
            _extrude_bl_face_infos(m).infos)
        snap_ctx = snap_ctx || face_flags_ctx.snap
        scale_bg_ctx = scale_bg_ctx || face_flags_ctx.scale
        faces_scale_ctx=_extrude_bl_face_infos(m).scale_last
        bg_ctx=get(_fields,:bg,nothing)
    end
    for (d,tag,spec) in pending_surf
        src=Int(get(spec,:src,0))
        src==0 && continue
        stag=abs(src)
        haskey(m.surfaces,stag) || continue
        srcmesh=_extrude_bl_srcmesh(field_model,stag;
            size_field=_size_field,_mesh_input=_mesh_input,
            _mesh1d_options=_mesh1d_options)
        index=Int(get(spec,:index,0))
        field=get(_fields.fields,index,_EMPTY_BL_FIELD)
        scale_field=get(_fields.scales,index,nothing)
        scale_active=scale_field!==nothing && get(spec,:scale_last,false)
        params=get(m.meshing.extrude,(d,tag),nothing)
        params===nothing && continue
        offsets=_extrude_level_us(params)[2:end]
        nl=length(offsets)
        levels=Dict{NTuple{3,Float64},Vector{NTuple{3,Float64}}}()
        scales=scale_active ?
            Dict{NTuple{3,Float64},Float64}() : nothing
        hprev=length(params.layers)>=2 ? params.heights[end-1] : 0.0
        t_top=params.heights[end]
        for i in 1:size(srcmesh.coords,2)
            pos=(srcmesh.coords[1,i],srcmesh.coords[2,i],srcmesh.coords[3,i])
            n=get(field,pos,(0.0,0.0,0.0))
            # `nn[0] = 0` — a position the field never covered leaves every
            # level coincident with the source vertex.
            row=n==(0.0,0.0,0.0) ? fill(pos,nl) :
                [(pos[1]+t*n[1],pos[2]+t*n[2],pos[3]+t*n[3]) for t in offsets]
            if scale_active
                s=get(scale_field,pos,1.0)
                (isfinite(s) && abs(s-1.0)<=1e-12) && (s=1.0)
                scales[pos]=s
                # `Extrude(iLayer==NbLayer-1, scale)` — every element level
                # inside the last layer group remaps its `u` as
                # `u = h[N-2] + (u-h[N-2])·scale` (`u·scale` when NbLayer≤1,
                # i.e. `hprev==0`), not just the final row.
                if s!=1.0 && n!=(0.0,0.0,0.0)
                    for l in eachindex(offsets)
                        t=offsets[l]
                        t>hprev || continue
                        t2=hprev+(t-hprev)*s
                        row[l]=(pos[1]+t2*n[1],pos[2]+t2*n[2],
                                pos[3]+t2*n[3])
                    end
                end
            end
            levels[pos]=row
        end
        # `tops` keeps the top-row mapping `_extrude_at` and the volume
        # sweep already read; `levels` serves per-level evaluation.
        tops=Dict{NTuple{3,Float64},NTuple{3,Float64}}(
            pos=>row[end] for (pos,row) in levels)
        view_ctx=Int(get(spec,:view,-1))
        bl_ctx=(seen=Set{UInt}((objectid(srcmesh),)),
            view=view_ctx,
            sampler=view_ctx>=0 && samplers_ctx!==nothing ?
                get(samplers_ctx,view_ctx,nothing) : nothing,
            gated=haskey(_fields.scales,index),
            skip_scale=!get(faces_scale_ctx,stag,false),
            invert=src<0,
            unitize=rawnorm_ctx,
            snap=rawnorm_ctx && snap_ctx,
            scale_bg=rawnorm_ctx && scale_bg_ctx,
            bg=bg_ctx,
            scale_active=scale_active,
            hprev=hprev,t_top=t_top)
        m.meshing.extrude_specs[(d,tag)]=(type=:boundary_layer,bl3d=true,
            src=src,tops=tops,levels=levels,uvals=offsets,
            scales=scales,hprev=hprev,
            index=index,view=view_ctx,
            scale_last=get(spec,:scale_last,false),bl_ctx=bl_ctx)
    end
    # The `MSH_SURF_BND_LAYER` laterals of a 3-D boundary layer resolve the
    # same way, keyed on their source *edge's* stored chain — a subset of
    # the source face's mesh nodes, so the strip columns weld bitwise onto
    # the chapeau/volume rows. `levels` marks the spec resolved; the generic
    # `_extrude_lateral_mesh` product grid (quads or the
    # `SubdivideExtrudedMesh` triangle pairs) sweeps it.
    for (tag,spec) in pending_lat
        link=get(m.meshing.extrude_sources,(2,tag),nothing)
        link===nothing && continue
        signed_c=link[2];gen=abs(signed_c)
        haskey(m.curves,gen) || continue
        chain,_,orphans,_=_extrude_curve_source(m,gen,_mesh_input,caller;
            _mesh1d_options=_mesh1d_options)
        a,b=m.curves[gen];closed=a==b
        index=Int(get(spec,:index,0))
        field=get(_fields.fields,index,_EMPTY_BL_FIELD)
        scale_field=get(_fields.scales,index,nothing)
        params=get(m.meshing.extrude,(2,tag),nothing)
        params===nothing && continue
        scale_active=scale_field!==nothing && params.scale_last
        offsets=_extrude_level_us(params)[2:end]
        nl=length(offsets)
        hprev=length(params.layers)>=2 ? params.heights[end-1] : 0.0
        t_top=params.heights[end]
        dirs=Dict{NTuple{3,Float64},NTuple{3,Float64}}()
        levels=Dict{NTuple{3,Float64},Vector{NTuple{3,Float64}}}()
        for pos in chain
            n=get(field,pos,(0.0,0.0,0.0))
            dirs[pos]=n
            row=n==(0.0,0.0,0.0) ? fill(pos,nl) :
                [(pos[1]+t*n[1],pos[2]+t*n[2],pos[3]+t*n[3])
                 for t in offsets]
            if scale_active
                s=get(scale_field,pos,1.0)
                (isfinite(s) && abs(s-1.0)<=1e-12) && (s=1.0)
                if s!=1.0 && n!=(0.0,0.0,0.0)
                    for l in eachindex(offsets)
                        t=offsets[l]
                        t>hprev || continue
                        t2=hprev+(t-hprev)*s
                        row[l]=(pos[1]+t2*n[1],pos[2]+t2*n[2],
                                pos[3]+t2*n[3])
                    end
                end
            end
            levels[pos]=row
        end
        tops=Dict{NTuple{3,Float64},NTuple{3,Float64}}(
            pos=>row[end] for (pos,row) in levels)
        m.meshing.extrude_specs[(2,tag)]=(;spec...,levels=levels,tops=tops,
                                          uvals=offsets)
        _extrude_bl_surf_refresh!(m,spec,gen,signed_c,closed,dirs,tops,
            scale_active ? scale_field : nothing,hprev,caller)
    end
    for (tag,spec) in pending
        haskey(m.curves,tag) || continue
        a,b=m.curves[tag]
        src_pt=Int(get(spec,:src,a))
        haskey(m.points,src_pt) || continue
        pos=m.points[src_pt]
        index=Int(get(spec,:index,0))
        field=get(_fields.fields,index,_EMPTY_BL_FIELD)
        n=get(field,pos,(0.0,0.0,0.0))
        params=get(m.meshing.extrude,(1,tag),nothing)
        # `scaleLast && scaleFactor` — both the copy's own flag and the
        # index's armed scale gate must hold before the last interval scales.
        scale_field=get(_fields.scales,index,nothing)
        s=(scale_field===nothing || !get(spec,:scale_last,false)) ? 1.0 :
            get(scale_field,pos,1.0)
        (isfinite(s) && abs(s-1.0)<=1e-12) && (s=1.0)
        hprev=params===nothing || length(params.layers)<2 ? 0.0 :
              params.heights[end-1]
        m.meshing.extrude_specs[(1,tag)]=(type=:translate,T=n,
            scale=s==1.0 ? nothing : (hprev,s))
        # The `ExtrudePoint` copy moves to the connector's top level along
        # the resolved normal; a zero normal leaves it coincident.
        dest=a==src_pt ? b : a
        haskey(m.points,dest) || continue
        t_top=params===nothing ? 0.0 : params.heights[end]
        s!=1.0 && (t_top=hprev+(t_top-hprev)*s)
        m.points[dest]=(pos[1]+t_top*n[1],pos[2]+t_top*n[2],pos[3]+t_top*n[3])
    end
    return nothing
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
    blspec=spec isa NamedTuple && get(spec,:bl3d,false)
    srcmesh=source_mesh===nothing ?
        (blspec ?
            _extrude_bl_srcmesh(m,src;size_field=size_field,
                _mesh_input=_mesh_input) :
            mesh_model_surface(m,src;min_angle_deg=min_angle_deg,
                max_periodic_passes=max_periodic_passes,
                size_field=size_field,_mesh_input,
                _boundary_regrade=false)) :
        source_mesh
    ensured=_extrude_bl_ensure_srcmesh!(m,2,t,caller,srcmesh)
    ensured===nothing || (spec=ensured)
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
    # Replaying an already-meshed model sweeps the retained actual Source
    # product — upstream `extrudeMesh` reads `from->mesh_vertices`, the
    # mesh the source GFace currently carries — so a `mesh_model_volume`
    # on a generated model reproduces the pass's column order bitwise.
    # When no product is retained this falls through to a fresh mesh, the
    # same as an explicit `_mesh_input` miss.
    input=_mesh_input===nothing ? m.meshing.mesh_authority : _mesh_input
    srcmesh=spec isa NamedTuple && get(spec,:bl3d,false) ?
        _extrude_bl_srcmesh(m,src;_mesh_input=input) :
        mesh_model_surface(m,src;_mesh_input=input,_boundary_regrade=false)
    ensured=_extrude_bl_ensure_srcmesh!(m,3,t,caller,srcmesh)
    ensured===nothing || (spec=ensured)
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
                                caller::AbstractString;_mesh_input=nothing,
                                _mesh1d_options=nothing,_size_field=nothing)
    params=_extrude_gate(m,dim,tag)
    params===nothing && return nothing
    spec=get(m.meshing.extrude_specs,(dim,tag),nothing)
    link=get(m.meshing.extrude_sources,(dim,tag),nothing)
    (spec===nothing || link===nothing) && throw(ArgumentError(
        "$caller: $(dim==2 ? "Surface" : "Volume")[$tag] Layers " *
        "extrusion is missing recorded extrusion parameters"))
    if spec isa NamedTuple
        st=get(spec,:type,nothing)
        if st===:bl_surface ||
            (dim==2 && st===:boundary_layer &&
             Int(get(spec,:surf,0))!=0 && !haskey(spec,:levels))
            # Upstream resolves every boundary-layer copy's displaced rows
            # inside `Mesh2DWithBoundaryLayers` before the swept mesh reads
            # them — pending records finish their column fields here.
            _extrude_bl_prepare!(m,caller;_mesh_input,
                _mesh1d_options=_mesh1d_options,_size_field=_size_field)
            spec=get(m.meshing.extrude_specs,(dim,tag),nothing)
        end
    end
    return params,spec,link
end

# `MeshExtrudedSurface` — dispatch on `geo.Mode`: `(1, ·)` sources are
# EXTRUDED_ENTITY laterals, `(2, ·)` sources are COPIED_ENTITY tops.
function _extrude_surface_mesh(m::GeoModel,t::Int,caller::AbstractString;
                               min_angle_deg::Real=25.0,
                               max_periodic_passes=8,
                               size_field::Union{Nothing,AbstractSizeField}=nothing,
                               edges=nothing,scope=nothing,_mesh_input=nothing,_return_production=false,
                               _options=nothing)
    planned=_extrude_nonew_surface_mesh(m,t,caller;scope=scope,
        min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
        size_field=size_field)
    planned===nothing || return planned
    entry=_extrude_entity_params(m,2,t,caller;_mesh_input,
        _mesh1d_options=_options,_size_field=size_field)
    entry===nothing && return nothing
    params,spec,link=entry
    if spec isa NamedTuple && get(spec,:type,nothing)===:boundary_layer &&
        link[1]==1 && !haskey(spec,:levels)
        return _extrude_bl_lateral_mesh(m,t,params,spec,link[2],caller;
                                        _mesh_input,_size_field=size_field,
                                        _options=_options)
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
    entry=_extrude_entity_params(m,3,t,caller;_mesh_input,
        _size_field=size_field)
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
        params,spec,link=_extrude_entity_params(m,3,tag,caller;_mesh_input)
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
