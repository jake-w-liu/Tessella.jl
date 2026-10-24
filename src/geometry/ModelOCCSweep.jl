# `gmsh.model.occ.extrude` — `OCC_Internals::extrude`: a
# `BRepPrimAPI_MakePrism`-style translation sweep per input dimension.
# The bounded implementation covers the all-linear profiles exactly:
# points sweep to line edges, line edges to ruled parallelogram plane
# faces, and line-bounded plane faces to prism solids; a single-circle
# disk swept along its own normal reuses the compact cylinder
# publication (upstream shares the source face as the bottom cap — the
# producer builds a geometrically identical cap instead, so the source
# stays bound as an independent entity). `numElements`/`heights`/
# `recombine` engage the native Layers machinery through
# `_model_occ_sweep_transport!`, which mirrors upstream's
# `ExtrudeParams` attachment onto every created entity. Curved profiles
# other than disks, solid inputs and in-plane directions reject loudly.
# Output order follows upstream's `getReturnedShapes` reorder — per
# input `[top (dim-1), body (dim), laterals …]` — with dim-0 inputs
# emitting just the swept edge, all grouped by ascending dimension like
# `_extrude`'s per-dim compound pass.

function _model_occ_extrude_direction(dx,dy,dz,caller)
    d=_finite3(dx,dy,dz,caller)
    _occ_modulus(d)>0 || throw(ArgumentError(
        "$caller: extrusion direction must have positive length"))
    return d
end

# Per-call sweep bookkeeping: translated vertices/edges and swept
# lateral edges deduplicate across inputs (one compound like upstream),
# `cursor` hands out per-dim tags and resyncs when an eager producer
# (the cylinder path) advances `kernel.maxima` mid-call.
mutable struct _ModelOCCSweep
    d::NTuple{3,Float64}
    caller::String
    pending::Vector{Tuple{Int,Int,UInt64,_ModelOCCShape}}
    shapes::Dict{UInt64,_ModelOCCShape}
    tags::Dict{UInt64,Int}
    cursor::Vector{Int}
    next_id::UInt64
    vmap::Dict{UInt64,UInt64}
    emap::Dict{UInt64,UInt64}
    lmap::Dict{UInt64,UInt64}
    rmap::Dict{UInt64,UInt64}
    fmap::Dict{UInt64,UInt64}
    # `smap` carries the extrusion's swept faces (face id -> generatrix edge
    # id; `fmap` plays the same role for the revolve); `links` holds the
    # cap/solid/compact-publication source records that the maps cannot
    # express — `(dim, entity_tag) => (source_dim, source_tag)` mesh links
    # consumed by `_model_occ_sweep_transport!` when `numElements` engages
    # the Layers machinery.
    smap::Dict{UInt64,UInt64}
    links::Vector{Tuple{Int,Int,Tuple{Int,Int}}}
end

function _model_occ_sweep_setup(d,caller,kernel)
    return _ModelOCCSweep(d,caller,
        Tuple{Int,Int,UInt64,_ModelOCCShape}[],
        Dict{UInt64,_ModelOCCShape}(),Dict{UInt64,Int}(),
        collect(kernel.maxima),kernel.next_shape,
        Dict{UInt64,UInt64}(),Dict{UInt64,UInt64}(),
        Dict{UInt64,UInt64}(),Dict{UInt64,UInt64}(),
        Dict{UInt64,UInt64}(),Dict{UInt64,UInt64}(),
        Tuple{Int,Int,Tuple{Int,Int}}[])
end

# `next_shape` only advances at the flush, so identity allocation is
# sweep-scoped; it still resyncs defensively if an eager producer bumped
# `kernel.next_shape` mid-call, and the flush takes `max` against those
# eagerly-bound shapes.
function _model_occ_sweep_id!(sw::_ModelOCCSweep,kernel)
    sw.next_id=max(sw.next_id,kernel.next_shape)
    sw.next_id<=typemax(UInt64)-1 || throw(ArgumentError(
        "$(sw.caller): OCC shape identity allocation overflow"))
    id=sw.next_id
    sw.next_id+=1
    return id
end

function _model_occ_sweep_tag!(sw::_ModelOCCSweep,kernel,dim)
    index=dim+3
    tag=sw.cursor[index]=max(sw.cursor[index],kernel.maxima[index])+1
    tag<=typemax(Int32) || throw(ArgumentError(
        "$(sw.caller): OCC entity tag allocation exceeds Int32"))
    return tag
end

function _model_occ_sweep_bind!(sw::_ModelOCCSweep,dim,tag,id,shape)
    push!(sw.pending,(dim,tag,id,shape))
    sw.shapes[id]=shape
    sw.tags[id]=tag
    return nothing
end

function _model_occ_sweep_vertex!(sw::_ModelOCCSweep,kernel,src_id)
    haskey(sw.vmap,src_id) && return sw.vmap[src_id]
    source=kernel.shapes[src_id]
    point=_occ3_add(source.geometry.point,sw.d)
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,0)
    _model_occ_sweep_bind!(sw,0,tag,id,
        _ModelOCCShape(0,(point=point,
                          mesh_size=source.geometry.mesh_size),UInt64[]))
    sw.vmap[src_id]=id
    return id
end

# Swept vertex edge src → src+d.
function _model_occ_sweep_lateral_edge!(sw::_ModelOCCSweep,kernel,src_id)
    haskey(sw.lmap,src_id) && return sw.lmap[src_id]
    top=_model_occ_sweep_vertex!(sw,kernel,src_id)
    a=kernel.shapes[src_id].geometry.point
    b=sw.shapes[top].geometry.point
    diff=_arc_sub(b,a)
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:line,geometry=
            (occ=:line,origin=a,dir=_occ_dir(diff),
             t0=0.0,t1=_occ_modulus(diff))),
            UInt64[src_id,top]))
    sw.lmap[src_id]=id
    return id
end

# Translated copy of a line edge — the prism "top" for curve inputs.
function _model_occ_sweep_top_edge!(sw::_ModelOCCSweep,kernel,src_id,caller)
    haskey(sw.emap,src_id) && return sw.emap[src_id]
    source=kernel.shapes[src_id]
    g=source.geometry.geometry
    g.occ===:line || throw(ArgumentError(
        "$caller: extrusion of OCC curve kind $(g.occ) is not yet " *
        "supported"))
    a_id,b_id=source.children
    a=_model_occ_sweep_vertex!(sw,kernel,a_id)
    b=_model_occ_sweep_vertex!(sw,kernel,b_id)
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:line,geometry=
            (occ=:line,origin=_occ3_add(g.origin,sw.d),dir=g.dir,
             t0=g.t0,t1=g.t1)),
            UInt64[a,b]))
    sw.emap[src_id]=id
    return id
end

_model_occ_sweep_shape_of(sw::_ModelOCCSweep,kernel,id)=
    get(sw.shapes,id) do
        kernel.shapes[id]
    end

_model_occ_sweep_tag_of(sw::_ModelOCCSweep,kernel,id,caller)=
    get(sw.tags,id) do
        _model_occ_shape_tag(kernel,id,caller)
    end

# dim-0 input: one new top vertex and the swept line edge. Upstream's
# reorder gate needs dim>=1, so only the edge is reported.
function _model_occ_extrude_point!(sw::_ModelOCCSweep,kernel,src_id)
    src=_model_occ_sweep_lateral_edge!(sw,kernel,src_id)
    return sw.tags[src]
end

# dim-1 input: swept parallelogram face sharing the source edge as its
# bottom boundary, the translated copy as top, and two swept lateral
# edges. Reports (top edge, face, laterals) in profile order.
function _model_occ_extrude_line!(sw::_ModelOCCSweep,kernel,src_id,caller)
    source=kernel.shapes[src_id]
    a_id,b_id=source.children
    top=_model_occ_sweep_top_edge!(sw,kernel,src_id,caller)
    la=_model_occ_sweep_lateral_edge!(sw,kernel,a_id)
    lb=_model_occ_sweep_lateral_edge!(sw,kernel,b_id)
    uses=Tuple{UInt64,Bool}[(src_id,false),(lb,false),
                           (top,true),(la,true)]
    shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    geometry=_model_occ_plane_face_geometry(shape_of,tag_of,[uses],caller)
    wire_id=_model_occ_sweep_id!(sw,kernel)
    wire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
    _model_occ_sweep_bind!(sw,-1,wire_tag,wire_id,
        _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses]))
    face_id=_model_occ_sweep_id!(sw,kernel)
    face_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,face_tag,face_id,
        _ModelOCCShape(2,(kind=:plane,geometry=geometry),
                       UInt64[wire_id]))
    sw.smap[face_id]=src_id
    return sw.tags[top],face_tag,[sw.tags[la],sw.tags[lb]]
end

# Single-circle disk swept along its normal → the compact cylinder
# publication. Upstream shares the source face as the bottom cap; here
# the shell reuses it directly and the lateral wire reuses the source
# circle — the native Layers machinery's bitwise row welding then makes
# the swept strip close on the source entities' own cells (the seam
# runs through the profile's vertex, the rims grade at the source's
# phase). `sw` records carry the mesh-transport links.
function _model_occ_extrude_disk!(m::GeoModel,kernel,sw,face_id,d,caller)
    face=kernel.shapes[face_id]
    wire=kernel.shapes[only(face.children)]
    edge_id=only(wire.geometry.uses)[1]
    g=kernel.shapes[edge_id].geometry.geometry
    # Phase 1 already certified the full-circle, axis-aligned profile.
    radius=g.occ===:circle ? g.r : g.rx
    h=_occ_modulus(d);nd=_occ_dir(d)
    src_vid=kernel.shapes[edge_id].children[1]
    src_vtx=kernel.shapes[src_vid].geometry.point
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    src_tag=tag_of(face_id);circle_tag=tag_of(edge_id)
    vtx_tag=tag_of(src_vid)
    # `flips` marks an anti-normal sweep: the profile circle's parameter
    # runs backward around +d̂, so wire flags and pcurve legs mirror.
    flips=_dot3(nd,g.n)<0
    Xf,Yf=_occ_ax2_setx(nd,g.X)
    top_vtx=_occ3_add(src_vtx,d)
    top_vid=_model_occ_sweep_id!(sw,kernel)
    tv_tag=_model_occ_sweep_tag!(sw,kernel,0)
    _model_occ_sweep_bind!(sw,0,tv_tag,top_vid,
        _ModelOCCShape(0,(point=top_vtx,mesh_size=
            kernel.shapes[src_vid].geometry.mesh_size),UInt64[]))
    rim_id=_model_occ_sweep_id!(sw,kernel)
    rim_tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,rim_tag,rim_id,
        _ModelOCCShape(1,(kind=:circle,
            geometry=(occ=:circle,center=_occ3_add(g.center,d),n=nd,
                      X=Xf,Y=Yf,r=radius,t0=0.0,t1=_OCC_TWO_PI)),
            UInt64[top_vid,top_vid]))
    seam_id=_model_occ_sweep_id!(sw,kernel)
    seam_tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,seam_tag,seam_id,
        _ModelOCCShape(1,(kind=:line,
            geometry=(occ=:line,origin=src_vtx,dir=nd,t0=0.0,t1=h)),
            UInt64[src_vid,top_vid]))
    lat_pc=Dict{Int,NamedTuple}(
        rim_tag=>(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing),
        seam_tag=>(fwd=(lin2d=((_OCC_TWO_PI,0.0),(0.0,1.0)),),
                   rev=(lin2d=((0.0,0.0),(0.0,1.0)),)))
    lat_pc[circle_tag]=flips ?
        (fwd=(lin2d=((_OCC_TWO_PI,0.0),(-1.0,0.0)),),rev=nothing) :
        (fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing)
    lwire_id=_model_occ_sweep_id!(sw,kernel)
    lwire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
    _model_occ_sweep_bind!(sw,-1,lwire_tag,lwire_id,
        _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(rim_id,true),
            (seam_id,true),(edge_id,flips),(seam_id,false)],),
            UInt64[rim_id,seam_id,edge_id,seam_id]))
    lat_id=_model_occ_sweep_id!(sw,kernel)
    lat_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,lat_tag,lat_id,
        _ModelOCCShape(2,(kind=:cylinder,
            geometry=(occ=:cylinder,center=g.center,axis=nd,X=Xf,Y=Yf,
                      radius=radius,height=h,pcurves=lat_pc)),
            UInt64[lwire_id]))
    top_center=_occ3_add(g.center,d)
    twire_id=_model_occ_sweep_id!(sw,kernel)
    twire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
    _model_occ_sweep_bind!(sw,-1,twire_tag,twire_id,
        _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(rim_id,false)],),
                       UInt64[rim_id]))
    top_pc=Dict{Int,NamedTuple}(
        rim_tag=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),radius),),rev=nothing))
    cap_id=_model_occ_sweep_id!(sw,kernel)
    cap_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,cap_tag,cap_id,
        _ModelOCCShape(2,(kind=:plane,
            geometry=(occ=:plane,center=top_center,axis=nd,X=Xf,Y=Yf,
                      reversed=false,pcurves=top_pc,
                      uvb=_model_occ_prim_uvb(top_pc,
                          Dict(rim_tag=>(0.0,_OCC_TWO_PI))))),
            UInt64[twire_id]))
    shell_id=_model_occ_sweep_id!(sw,kernel)
    shell_tag=_model_occ_sweep_tag!(sw,kernel,-2)
    src_axis=face.geometry.geometry.axis
    _model_occ_sweep_bind!(sw,-2,shell_tag,shell_id,
        _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[(lat_id,false),
            (cap_id,false),(face_id,_dot3(src_axis,d)>0)],),
            UInt64[lat_id,cap_id,face_id]))
    solid_id=_model_occ_sweep_id!(sw,kernel)
    solid_tag=_model_occ_sweep_tag!(sw,kernel,3)
    _model_occ_sweep_bind!(sw,3,solid_tag,solid_id,
        _ModelOCCShape(3,(kind=:solid,
            primitive=(occ=:cylinder,
                record=(center=g.center,axis=d,radius=radius,height=h))),
            UInt64[shell_id]))
    push!(sw.links,(1,rim_tag,(1,flips ? -circle_tag : circle_tag)))
    push!(sw.links,(1,seam_tag,(0,vtx_tag)))
    push!(sw.links,(2,lat_tag,(1,circle_tag)))
    push!(sw.links,(2,cap_tag,(2,src_tag)))
    push!(sw.links,(3,solid_tag,(2,src_tag)))
    return cap_tag,solid_tag,[lat_tag]
end

# dim-2 input: all-line plane faces sweep to a prism — translated caps,
# one ruled lateral per wire edge, one exterior shell (hole laterals
# join the same shell since the swept void is open at the caps). Each
# new face record bakes its intended outward normal into `reversed`
# like `add_box`, so the shell enumerates all-forward uses except the
# shared source cap.
function _model_occ_extrude_planar!(sw::_ModelOCCSweep,kernel,src_id,caller)
    source=kernel.shapes[src_id]
    g=source.geometry.geometry
    # Solid consumers read `σ·axis` from the shell use flag, so the
    # source cap's outward flip keys on the raw stored axis.
    src_normal=g.axis
    shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    outer_normal=_model_occ_newell_normal(_model_occ_wire_samples(
        shape_of,kernel.shapes[first(source.children)].geometry.uses,
        caller))
    _occ_modulus(outer_normal)>_OCC_CONFUSION || throw(ArgumentError(
        "$caller: face boundary has no consistent normal"))
    outer_normal=_occ_dir(outer_normal)
    laterals=Tuple{UInt64,Int}[]
    cap_uses=Vector{Tuple{UInt64,Bool}}[]
    for (windex,wire_id) in enumerate(source.children)
        wire=kernel.shapes[wire_id]
        hole=windex>1
        cap_wire=Tuple{UInt64,Bool}[]
        for (edge_id,rev) in wire.geometry.uses
            a_id,b_id=kernel.shapes[edge_id].children
            a=kernel.shapes[a_id].geometry.point
            b=kernel.shapes[b_id].geometry.point
            tangent=_arc_sub(b,a)
            top=_model_occ_sweep_top_edge!(sw,kernel,edge_id,caller)
            push!(cap_wire,(top,rev))
            la=_model_occ_sweep_lateral_edge!(sw,kernel,a_id)
            lb=_model_occ_sweep_lateral_edge!(sw,kernel,b_id)
            # The quad traversal is built on the wire's traversal flag;
            # hole walls traverse opposite so their Newell lands on the
            # cavity's outward (into-material) side directly — `axis`
            # carries the effective normal for every solid consumer.
            f=hole ? !rev : rev
            luses=f ? Tuple{UInt64,Bool}[(edge_id,true),(la,false),
                                         (top,false),(lb,true)] :
                      Tuple{UInt64,Bool}[(edge_id,false),(lb,false),
                                         (top,true),(la,true)]
            lgeo=_model_occ_plane_face_geometry(
                shape_of,tag_of,[luses],caller)
            udir=_occ_dir(f ? _arc_sub(a,b) : tangent)
            intended=_occ_crossed(udir,outer_normal)
            # The traversal's own Newell is û_f×d̂ — it lands on the
            # intended side only when d̂ agrees with the profile normal;
            # the mirror traversal corrects the anti-normal sweep.
            if _dot3(lgeo.axis,intended)<0
                luses=Tuple{UInt64,Bool}[(e,!r) for (e,r) in
                    Iterators.reverse(luses)]
                lgeo=_model_occ_plane_face_geometry(
                    shape_of,tag_of,[luses],caller)
            end
            _dot3(lgeo.axis,intended)>0.999 || throw(ArgumentError(
                "$caller: internal lateral face winding is inconsistent"))
            lwire_id=_model_occ_sweep_id!(sw,kernel)
            lwire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
            _model_occ_sweep_bind!(sw,-1,lwire_tag,lwire_id,
                _ModelOCCShape(-1,(uses=luses,),
                              UInt64[u[1] for u in luses]))
            lface_id=_model_occ_sweep_id!(sw,kernel)
            lface_tag=_model_occ_sweep_tag!(sw,kernel,2)
            _model_occ_sweep_bind!(sw,2,lface_tag,lface_id,
                _ModelOCCShape(2,(kind=:plane,geometry=lgeo),
                              UInt64[lwire_id]))
            sw.smap[lface_id]=edge_id
            push!(laterals,(lface_id,lface_tag))
        end
        push!(cap_uses,cap_wire)
    end
    # When the sweep runs against the profile normal the translated cap
    # must wind opposite so its Newell — hence `axis` — points along d̂.
    cap_geometry=_model_occ_plane_face_geometry(
        shape_of,tag_of,cap_uses,caller)
    if _dot3(cap_geometry.axis,sw.d)<0
        cap_uses=[Tuple{UInt64,Bool}[(e,!r) for (e,r) in
                    Iterators.reverse(cap)] for cap in cap_uses]
        cap_geometry=_model_occ_plane_face_geometry(
            shape_of,tag_of,cap_uses,caller)
    end
    cap_wire_ids=UInt64[]
    for cap in cap_uses
        wire_id=_model_occ_sweep_id!(sw,kernel)
        wire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
        _model_occ_sweep_bind!(sw,-1,wire_tag,wire_id,
            _ModelOCCShape(-1,(uses=cap,),UInt64[u[1] for u in cap]))
        push!(cap_wire_ids,wire_id)
    end
    cap_id=_model_occ_sweep_id!(sw,kernel)
    cap_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,cap_tag,cap_id,
        _ModelOCCShape(2,(kind=:plane,geometry=cap_geometry),
                       cap_wire_ids))
    shell_uses=Tuple{UInt64,Bool}[(cap_id,false),
        (src_id,_dot3(src_normal,sw.d)>0)]
    for (fid,_) in laterals
        push!(shell_uses,(fid,false))
    end
    shell_id=_model_occ_sweep_id!(sw,kernel)
    shell_tag=_model_occ_sweep_tag!(sw,kernel,-2)
    _model_occ_sweep_bind!(sw,-2,shell_tag,shell_id,
        _ModelOCCShape(-2,(uses=shell_uses,),
                      UInt64[u[1] for u in shell_uses]))
    solid_id=_model_occ_sweep_id!(sw,kernel)
    solid_tag=_model_occ_sweep_tag!(sw,kernel,3)
    _model_occ_sweep_bind!(sw,3,solid_tag,solid_id,
        _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id]))
    src_tag=tag_of(src_id)
    push!(sw.links,(2,cap_tag,(2,src_tag)))
    push!(sw.links,(3,solid_tag,(2,src_tag)))
    return cap_tag,solid_tag,[tag for (_,tag) in laterals]
end

# `numElements`/`heights`/`recombine` — the OCC sweep meshing options.
# Upstream attaches `ExtrudeParams` (`mesh.ExtrudeMesh`) to every entity the
# sweep creates; here the equivalent `_GeoExtrudeParams` record plus the
# `(dim,tag)`-keyed source/spec links feed the existing Julia-native Layers
# machinery once the shapes materialize. `heights` are the upstream
# normalized cumulative layer boundaries — verbatim into `params.heights`,
# defaulting to uniform spacing when omitted. Without `numElements` the
# options are inert (upstream's `ExtrudeMesh` gate stays closed) except a
# bare `heights` list, which is a clear argument error.
function _model_occ_sweep_mesh_params(num_elements,heights,recombine,caller)
    isempty(num_elements) && !isempty(heights) && throw(ArgumentError(
        "$caller: heights requires numElements"))
    isempty(num_elements) && return nothing
    layers=Int[]
    for e in num_elements
        (e isa Integer || (e isa Real && isinteger(e))) && e>0 ||
            throw(ArgumentError(
                "$caller: numElements entries must be positive integers"))
        push!(layers,Int(e))
    end
    hs=isempty(heights) ?
        Float64[j/length(layers) for j in 1:length(layers)] :
        Float64[_finite_scalar(h,caller,"heights entry") for h in heights]
    length(hs)==length(layers) || throw(ArgumentError(
        "$caller: heights must have one cumulative entry per numElements " *
        "layer ($(length(hs)) heights for $(length(layers)) layers)"))
    for j in eachindex(hs)
        hs[j]>(j==1 ? 0.0 : hs[j-1]) || throw(ArgumentError(
            "$caller: heights must be strictly increasing positive " *
            "fractions normalized to 1"))
    end
    hs[end]==1.0 || throw(ArgumentError(
        "$caller: heights must be normalized to 1 — the last cumulative " *
        "layer boundary is the swept cap"))
    return (layers=layers,heights=hs,scale_last=false,
            recombine=Bool(recombine),quad_to_tri=:none,
            recomb_laterals=false)::_GeoExtrudeParams
end

# Attach the `ExtrudeParams` records to every created entity. The maps carry
# shape ids resolved through `sw.tags`; `sw.links` holds tag-level records
# for caps, solids and the compact cylinder/torus publications.
function _model_occ_sweep_transport!(m::GeoModel,kernel,
                                     sw::_ModelOCCSweep,
                                     params::_GeoExtrudeParams,spec)
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,sw.caller)
    function attach(dim::Int,tag::Int,src::NTuple{2,Int})
        m.meshing.extrude[(dim,tag)]=params
        m.meshing.extrude_sources[(dim,tag)]=src
        m.meshing.extrude_specs[(dim,tag)]=spec
    end
    # Vertex-swept edges: extrusion connectors, revolve rim circles, and
    # apex degenerate edges all take the point-generatrix link.
    for (vid,eid) in sw.lmap
        eid==0 && continue
        attach(1,sw.tags[eid],(0,tag_of(vid)))
    end
    for (vid,cid) in sw.rmap
        cid==0 && continue
        attach(1,sw.tags[cid],(0,tag_of(vid)))
    end
    # Translated top-edge copies take the copyMesh link on the source edge.
    for (eid,teid) in sw.emap
        attach(1,sw.tags[teid],(1,tag_of(eid)))
    end
    # Swept surfaces take the EXTRUDED_ENTITY link on their generatrix
    # edge — `smap` keys the face id to its generatrix edge id.
    for (fid,eid) in sw.smap
        attach(2,sw.tags[fid],(1,tag_of(eid)))
    end
    for (eid,fid) in sw.fmap
        attach(2,sw.tags[fid],(1,tag_of(eid)))
    end
    for (dim,tag,src) in sw.links
        attach(dim,tag,src)
    end
    return nothing
end

function _model_occ_extrude(m::GeoModel,dim_tags,dx,dy,dz;
                            num_elements=Int[],heights=Float64[],
                            recombine=false,caller="OCC.extrude")
    params=_model_occ_sweep_mesh_params(num_elements,heights,recombine,caller)
    d=_model_occ_extrude_direction(dx,dy,dz,caller)
    kernel=m.occ_kernel
    # Phase 1: validate and classify every input before any bind so a
    # rejection leaves the registry untouched.
    classified=Tuple{Int,UInt64,Symbol}[]
    for entry in dim_tags
        dim=_model_occ_dimension(first(entry),caller)
        0<=dim<=3 || throw(ArgumentError(
            "$caller: Wrong input dimension in extrusion"))
        tag=_model_occ_requested_tag(last(entry),caller)
        id=get(kernel.bindings,(dim,tag),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $dim " *
            "with tag $tag"))
        kind=if dim==0
            :point
        elseif dim==1
            g=kernel.shapes[id].geometry
            g.kind===:line || throw(ArgumentError(
                "$caller: extrusion of OCC curve kind $(g.kind) is " *
                "not yet supported"))
            a_id,b_id=kernel.shapes[id].children
            tangent=_arc_sub(
                kernel.shapes[b_id].geometry.point,
                kernel.shapes[a_id].geometry.point)
            cross=_occ_cross(tangent,d)
            _occ_modulus(cross)>_OCC_CONFUSION*_occ_modulus(tangent)*
                _occ_modulus(d) || throw(ArgumentError(
                "$caller: extrusion direction is parallel to the curve"))
            :line
        elseif dim==2
            g=kernel.shapes[id].geometry
            g.kind===:plane || throw(ArgumentError(
                "$caller: extrusion of OCC surface kind $(g.kind) is " *
                "not yet supported"))
            all_lines=true;edges=0
            for wire_id in kernel.shapes[id].children
                for (edge_id,_) in
                    kernel.shapes[wire_id].geometry.uses
                    edges+=1
                    eshape=kernel.shapes[edge_id]
                    if eshape.geometry.kind===:line
                        a_id,b_id=eshape.children
                        tangent=_arc_sub(
                            kernel.shapes[b_id].geometry.point,
                            kernel.shapes[a_id].geometry.point)
                        cross=_occ_cross(tangent,d)
                        _occ_modulus(cross)>_OCC_CONFUSION*
                            _occ_modulus(tangent)*_occ_modulus(d) ||
                            throw(ArgumentError(
                            "$caller: extrusion direction is " *
                            "parallel to a boundary curve"))
                    else
                        all_lines=false
                    end
                end
            end
            if all_lines
                abs(_dot3(g.geometry.axis,d))>
                    _OCC_CONFUSION*_occ_modulus(d) || throw(
                    ArgumentError(
                    "$caller: extrusion direction is coplanar with " *
                    "the surface"))
                :planar
            elseif edges==1
                eid=only(kernel.shapes[only(
                    kernel.shapes[id].children)].geometry.uses)[1]
                eg=kernel.shapes[eid].geometry.geometry
                if eg.occ===:circle ||
                        (eg.occ===:ellipse &&
                         abs(eg.rx-eg.ry)<=_OCC_CONFUSION*eg.rx)
                    eg.t1-eg.t0==_OCC_TWO_PI || throw(ArgumentError(
                        "$caller: extrusion of partial-circle " *
                        "profiles is not yet supported"))
                    _occ_modulus(_occ_cross(d,eg.n))<
                        _OCC_CONFUSION*_occ_modulus(d) || throw(
                        ArgumentError(
                        "$caller: extrusion of a disk must follow " *
                        "its surface normal"))
                    :disk
                elseif eg.occ===:ellipse
                    throw(ArgumentError(
                        "$caller: extrusion of elliptical profiles " *
                        "is not yet supported"))
                else
                    :unsupported
                end
            else
                :unsupported
            end
        else
            throw(ArgumentError(
                "$caller: extrusion of OCC volumes is not yet supported"))
        end
        kind===:unsupported && throw(ArgumentError(
            "$caller: extrusion of this OCC entity is not yet supported"))
        push!(classified,(dim,id,kind))
    end
    sw=_model_occ_sweep_setup(d,caller,kernel)
    out=Tuple{Int,Int}[]
    # Phase 2: construct per ascending dimension like upstream's
    # `_extrude` per-dim compound pass.
    for dim in (0,1,2)
        for (edim,id,kind) in classified
            edim==dim || continue
            if dim==0
                push!(out,(1,_model_occ_extrude_point!(sw,kernel,id)))
            elseif dim==1
                top,body,laterals=_model_occ_extrude_line!(
                    sw,kernel,id,caller)
                push!(out,(1,top));push!(out,(2,body))
                for tag in laterals
                    push!(out,(1,tag))
                end
            else
                top,body,laterals=kind===:disk ?
                    _model_occ_extrude_disk!(m,kernel,sw,id,d,caller) :
                    _model_occ_extrude_planar!(sw,kernel,id,caller)
                push!(out,(2,top));push!(out,(3,body))
                for tag in laterals
                    push!(out,(2,tag))
                end
            end
        end
    end
    for (dimension,tag,id,shape) in sw.pending
        _model_occ_bind_shape!(kernel,dimension,tag,id,shape)
    end
    isempty(sw.pending) ||
        (kernel.next_shape=max(kernel.next_shape,sw.next_id))
    params===nothing || _model_occ_sweep_transport!(
        m,kernel,sw,params,_extrude_spec_translate(d))
    return out
end

# ---------------------------------------------------------------------------
# `gmsh.model.occ.revolve` — `OCC_Internals::revolve`: a
# `BRepPrimAPI_MakeRevol`-style rotational sweep per input dimension.
# The bounded implementation covers sweeps of axial-plane profiles
# exactly: points sweep to rim circles (arcs on partial turns), line
# edges to cylinder/cone laterals or annulus/disk plane faces,
# single-wire all-line plane faces to revolved solids (cylinder/cone
# compact encodings when the profile produces a primitive, `sector`
# records on partial turns), and a full-circle disk pierced by the
# axis to the compact torus publication (a torus segment with two
# meridian disk caps on partial turns). At 2·Pi the rotated "top"
# coincides with the profile — the source entity answers as `top` and
# both meridian legs reuse the source edge; at partial angles the top
# is the rotated copy and the two walls are the source face and its
# rotated image. Non-axial-plane profiles, axis-crossing profiles,
# curved profiles other than disks, solid inputs and multi-wire faces
# reject loudly. `numElements`/`heights`/`recombine` engage the same
# native Layers transport as the extrusion. Rim curves, meridian
# faces, rotated copies and degenerate apex edges deduplicate across
# shared profile entities through `rmap`/`fmap`/`vmap`/`emap` like the
# prism sweep's maps.

# (rho,h,radial vector) of a point in the revolve frame: h is the axial
# coordinate along unit axis n, rho the distance to the axis line.
function _model_occ_revol_rh(pt,p,n)
    w=_arc_sub(pt,p);h=_dot3(w,n)
    r=_occ3_add(w,_occ3_scale(n,-h))
    return _occ_modulus(r),h,r
end

# Rotated copy of a profile vertex at the sweep end `a` (partial turns
# only — at 2·Pi the source vertex answers). An on-axis vertex rotates
# to itself, so the source id answers directly and stays shared with the
# profile edges that reference it. Deduplicated through `vmap` like the
# prism sweep's translated copies. The position takes the spec's step
# chain — not the gp_Trsf compact form — so the vertex welds bitwise
# onto the Layers emitter's closing row (the two formulas differ in the
# last ulp through the 1-cos assembly).
function _model_occ_revol_vertex!(sw,kernel,src_id,p,n,Ma,loca,spec)
    haskey(sw.vmap,src_id) && return sw.vmap[src_id]
    source=kernel.shapes[src_id]
    pt=source.geometry.point
    rho,_,_=_model_occ_revol_rh(pt,p,n)
    if rho<=_OCC_CONFUSION
        sw.vmap[src_id]=src_id
        return src_id
    end
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,0)
    _model_occ_sweep_bind!(sw,0,tag,id,
        _ModelOCCShape(0,(point=_extrude_rotate_point(
                              spec.axis,spec.origin,spec.angle,pt),
                          mesh_size=source.geometry.mesh_size),UInt64[]))
    sw.vmap[src_id]=id
    return id
end

# Rotated copy of a profile line edge — the revolve "top" edge at the
# sweep end (partial turns), deduplicated through `emap` like the prism
# sweep's translated tops.
function _model_occ_revol_top_edge!(sw,kernel,src_id,p,n,Ma,loca,spec,
                                    caller)
    haskey(sw.emap,src_id) && return sw.emap[src_id]
    source=kernel.shapes[src_id]
    g=source.geometry.geometry
    g.occ===:line || throw(ArgumentError(
        "$caller: revolve of OCC curve kind $(g.occ) is not yet " *
        "supported"))
    a_id,b_id=source.children
    a=_model_occ_revol_vertex!(sw,kernel,a_id,p,n,Ma,loca,spec)
    b=_model_occ_revol_vertex!(sw,kernel,b_id,p,n,Ma,loca,spec)
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:line,geometry=
            (occ=:line,origin=_occ_xform_pnt(Ma,loca,g.origin),
             dir=_occ_xform_dir(Ma,g.dir),t0=g.t0,t1=g.t1)),
            UInt64[a,b]))
    sw.emap[src_id]=id
    return id
end

# Rim circle (full turn) or rim arc (partial turn over `[0,a]`) swept by
# a profile vertex, deduplicated through `rmap`; an on-axis vertex maps
# to the zero sentinel and is skipped by the callers (it sweeps to
# itself, not a curve). Partial arcs end at the rotated copy vertex.
function _model_occ_revol_rim!(sw,kernel,vid,p,n,a,Ma,loca,spec,caller)
    haskey(sw.rmap,vid) && return sw.rmap[vid]
    v=kernel.shapes[vid].geometry.point
    rho,h,r=_model_occ_revol_rh(v,p,n)
    if rho<=_OCC_CONFUSION
        sw.rmap[vid]=zero(UInt64)
        return zero(UInt64)
    end
    Xc,Yc=_occ_ax2_setx(n,_occ_dir(r))
    end_vid=a<_OCC_TWO_PI ?
        _model_occ_revol_vertex!(sw,kernel,vid,p,n,Ma,loca,spec) : vid
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:circle,
            geometry=(occ=:circle,
                      center=_occ3_add(p,_occ3_scale(n,h)),
                      n=_occ_crossed(Xc,Yc),X=Xc,Y=Yc,r=rho,
                      t0=0.0,t1=a)),UInt64[vid,end_vid]))
    sw.rmap[vid]=id
    return id
end

# Degenerate rim edge for an apex vertex (rho==0) on a cone lateral —
# the producer convention binds the collapsed meridian limit to the
# (shared) profile vertex. Partial sweeps collapse over `[0,a]`.
function _model_occ_revol_apex!(sw,kernel,vid,face_geometry,a,caller)
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:degenerate,
            geometry=(occ=:degenerate,t0=0.0,t1=a,
                      face=face_geometry)),UInt64[vid,vid]))
    # The degenerate collapsed rim is the on-axis vertex's swept edge — the
    # point-generatrix link applies like the rim circles' (`lmap` carries
    # vertex->edge for the revolve; the extrusion never sees it).
    sw.lmap[vid]=id
    return id
end

# lin2d origin shifted so the pcurve evaluates correctly on an edge
# whose own parameter starts at t0 (lin2d evals at absolute t).
_occ_pcurve_shift(o,d,t0)=(o[1]-d[1]*t0,o[2]-d[2]*t0)

# Pappus moments of the revolved profile polygon: over fan triangles
# from the first vertex, ∫∫ρ dA drives the volume (a·∫∫ρ dA over the
# sweep angle `a`) and the ρ-weighted axial moment ∫∫h·ρ dA drives the
# centroid's axial coordinate — for a linear function the triangle mean
# is the vertex average, and for the bilinear h·ρ the mean is
# [(Σh)(Σρ)+Σhᵢρᵢ]/12. For partial turns the centroid also rides the
# sector bisector at radius 2·sin(a/2)·∫∫ρ² dA/(a·J0) — the quadratic
# ρ² integrates exactly through the edge-midpoint rule
# A/3·Σ ρ(mᵢⱼ)² with ρ(mᵢⱼ)=(ρᵢ+ρⱼ)/2 for linear ρ. Signed terms
# cancel in the ratios; V takes |a·J0|.
function _model_occ_revol_profile_moments(verts,a)
    r0,h0=verts[1]
    J0=J1=J2=0.0
    for i in 2:length(verts)-1
        r1,h1=verts[i];r2,h2=verts[i+1]
        A2=(r1-r0)*(h2-h0)-(r2-r0)*(h1-h0)
        sr=r0+r1+r2;sh=h0+h1+h2
        J0+=A2*sr/6
        J1+=A2*(sh*sr+r0*h0+r1*h1+r2*h2)/24
        J2+=A2*((r0+r1)^2+(r1+r2)^2+(r2+r0)^2)/24
    end
    s2,c2=_gm_sincos(a/2)
    arm=J0==0.0 ? 0.0 : 2*s2*J2/(a*J0)
    return abs(a*J0),J0==0.0 ? 0.0 : J1/J0,arm
end

# Meridian face swept by a profile line edge: the canonical
# cylinder/cone lateral record (ascending meridian, X at the seam's
# radial direction) with the source edge reused as seam and rim
# circles (or a degenerate apex edge) per endpoint. Partial turns carry
# `angle=a` on the lateral record, trim the rim circles to arcs, and
# split the seam into the source edge at u=0 plus its rotated copy at
# u=a. `fmap` deduplicates per profile edge. Returns the swept face's
# shape id.
function _model_occ_revol_meridian!(sw,kernel,src_id,p,n,a,Ma,loca,spec,
                                    caller)
    haskey(sw.fmap,src_id) && return sw.fmap[src_id]
    source=kernel.shapes[src_id]
    g=source.geometry.geometry
    a_id,b_id=source.children
    pa=kernel.shapes[a_id].geometry.point
    pb=kernel.shapes[b_id].geometry.point
    ra,ha,rd_a=_model_occ_revol_rh(pa,p,n)
    rb,hb,rd_b=_model_occ_revol_rh(pb,p,n)
    len=_occ_modulus(_arc_sub(pb,pa))
    rhat=_occ_dir(ra>=rb ? rd_a : rd_b)
    lo_vid,hi_vid=ha<=hb ? (a_id,b_id) : (b_id,a_id)
    rho_lo,rho_hi=ha<=hb ? (ra,rb) : (rb,ra)
    h_lo=min(ha,hb);span=abs(hb-ha)
    Xc,Yc=_occ_ax2_setx(n,rhat)
    cen=_occ3_add(p,_occ3_scale(n,h_lo))
    cylinder=abs(rho_hi-rho_lo)<=_OCC_CONFUSION*max(rho_hi,1.0)
    full=a==_OCC_TWO_PI
    seam_tag=_model_occ_sweep_tag_of(sw,kernel,src_id,caller)
    rot_id=full ? zero(UInt64) :
        _model_occ_revol_top_edge!(sw,kernel,src_id,p,n,Ma,loca,spec,
                                   caller)
    rot_tag=full ? 0 :
        _model_occ_sweep_tag_of(sw,kernel,rot_id,caller)
    # The meridian-leg pcurves map the edge's own parameter to the
    # meridian v: ascending profiles run 0->len, descending len->0. A
    # full turn reuses the source edge on both seam columns (u=0 rev,
    # u=2*Pi fwd); a partial turn puts the source edge at u=0 and its
    # rotated copy at u=a — each column map serves both directions
    # since each edge is used once.
    gspan=g.t1-g.t0
    sc=(a_id==lo_vid ? 1.0 : -1.0)*len/gspan
    o_seam=_occ_pcurve_shift(a_id==lo_vid ? (0.0,0.0) : (0.0,len),
                           (0.0,sc),g.t0)
    lat_pcurves=Dict{Int,NamedTuple}()
    if full
        lat_pcurves[seam_tag]=
            (fwd=(lin2d=((_OCC_TWO_PI,o_seam[2]),(0.0,sc)),),
             rev=(lin2d=((0.0,o_seam[2]),(0.0,sc)),))
    else
        lat_pcurves[seam_tag]=
            (fwd=(lin2d=((0.0,o_seam[2]),(0.0,sc)),),
             rev=(lin2d=((0.0,o_seam[2]),(0.0,sc)),))
        lat_pcurves[rot_tag]=
            (fwd=(lin2d=((a,o_seam[2]),(0.0,sc)),),
             rev=(lin2d=((a,o_seam[2]),(0.0,sc)),))
    end
    c_lo=c_hi=zero(UInt64)
    lat_geometry=nothing
    if cylinder
        lat_geometry=(occ=:cylinder,center=cen,axis=n,X=Xc,Y=Yc,
                      radius=rho_lo,height=span,angle=a,
                      pcurves=lat_pcurves)
    else
        lat_geometry=(occ=:cone,center=cen,axis=n,X=Xc,Y=Yc,
                      r1=rho_lo,r2=rho_hi,height=span,angle=a,
                      pcurves=lat_pcurves)
    end
    for (vid,vpos,is_hi) in ((lo_vid,0.0,false),(hi_vid,len,true))
        rho=vid==lo_vid ? rho_lo : rho_hi
        if rho<=_OCC_CONFUSION
            eid=_model_occ_revol_apex!(sw,kernel,vid,lat_geometry,a,caller)
        else
            eid=_model_occ_revol_rim!(sw,kernel,vid,p,n,a,Ma,loca,spec,
                                      caller)
        end
        etag=_model_occ_sweep_tag_of(sw,kernel,eid,caller)
        lat_pcurves[etag]=(fwd=(lin2d=((0.0,vpos),(1.0,0.0)),),
                           rev=nothing)
        is_hi ? (c_hi=eid) : (c_lo=eid)
    end
    # Producer lateral wire: [rim_hi rev, u=0 leg down, rim_lo fwd,
    # u=u_max leg up]. Full turns reuse the source edge on both legs;
    # partial turns give the u=a leg the rotated copy. Leg flags adapt
    # to the source edge's own child order (children run lo->hi like
    # the producers, or hi->lo).
    wire=Tuple{UInt64,Bool}[(c_hi,true)]
    push!(wire,(src_id,a_id==lo_vid))
    push!(wire,(c_lo,false))
    push!(wire,(full ? src_id : rot_id,!(a_id==lo_vid)))
    wire_id=_model_occ_sweep_id!(sw,kernel)
    wire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
    _model_occ_sweep_bind!(sw,-1,wire_tag,wire_id,
        _ModelOCCShape(-1,(uses=wire,),UInt64[u[1] for u in wire]))
    face_id=_model_occ_sweep_id!(sw,kernel)
    face_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,face_tag,face_id,
        _ModelOCCShape(2,(kind=cylinder ? :cylinder : :cone,
                         geometry=lat_geometry),UInt64[wire_id]))
    sw.fmap[src_id]=face_id
    return face_id
end

# Annulus (or disk) plane face swept by a radial profile edge at axial
# coordinate h — two concentric circles for a full turn, one open
# sector wire (outer arc, rotated radial leg, inner arc, source radial
# leg) for a partial turn. `axis_sign` chooses the stored plane axis =
# ±n — the caller derives it from the profile direction (standalone
# faces) or the required outward normal (solid caps). Returns the face
# id.
function _model_occ_revol_annulus!(sw,kernel,src_id,p,n,a,Ma,loca,spec,
                                   axis_sign,caller)
    haskey(sw.fmap,src_id) && return sw.fmap[src_id]
    source=kernel.shapes[src_id]
    a_id,b_id=source.children
    pa=kernel.shapes[a_id].geometry.point
    pb=kernel.shapes[b_id].geometry.point
    ra,ha,_=_model_occ_revol_rh(pa,p,n)
    rb,hb,_=_model_occ_revol_rh(pb,p,n)
    c_in=c_out=zero(UInt64)
    r_lo=min(ra,rb);r_hi=max(ra,rb)
    c_out=_model_occ_revol_rim!(
        sw,kernel,ra>=rb ? a_id : b_id,p,n,a,Ma,loca,spec,caller)
    r_lo>_OCC_CONFUSION &&
        (c_in=_model_occ_revol_rim!(
            sw,kernel,ra>=rb ? b_id : a_id,p,n,a,Ma,loca,spec,caller))
    shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    full=a==_OCC_TWO_PI
    uses=Vector{Tuple{UInt64,Bool}}[]
    if full
        # Wind the outer circle so the face's Newell axis lands on
        # axis_sign·n; hole wires keep the same stored orientation (the
        # children index marks them, like `add_plane_surface` users).
        flag=axis_sign<0
        push!(uses,Tuple{UInt64,Bool}[(c_out,flag)])
        c_in!=0 && push!(uses,Tuple{UInt64,Bool}[(c_in,flag)])
    else
        # Sector wire traversed CCW about +n (the plane fit normalizes
        # conic orientations onto its axis): the outer rim arc spans
        # u=0->a, the rotated radial leg descends to the inner arc
        # (which spans u=a->0 reversed), and the source radial leg
        # closes back out. Child order decides each leg's flag.
        rot_id=_model_occ_revol_top_edge!(
            sw,kernel,src_id,p,n,Ma,loca,spec,caller)
        wire=Tuple{UInt64,Bool}[(c_out,false),(rot_id,ra<rb)]
        c_in!=0 && push!(wire,(c_in,true))
        push!(wire,(src_id,ra>=rb))
        push!(uses,wire)
    end
    geometry=_model_occ_plane_face_geometry(shape_of,tag_of,uses,caller)
    wire_ids=UInt64[]
    for use in uses
        wire_id=_model_occ_sweep_id!(sw,kernel)
        wire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
        _model_occ_sweep_bind!(sw,-1,wire_tag,wire_id,
            _ModelOCCShape(-1,(uses=use,),UInt64[u[1] for u in use]))
        push!(wire_ids,wire_id)
    end
    face_id=_model_occ_sweep_id!(sw,kernel)
    face_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,face_tag,face_id,
        _ModelOCCShape(2,(kind=:plane,geometry=geometry),wire_ids))
    sw.fmap[src_id]=face_id
    return face_id
end

# dim-0 input: the vertex sweeps a rim circle (full turn) or rim arc
# (partial turn). Upstream emits [top vertex, swept edge]; at 2·Pi the
# top is the source vertex, at partial angles the rotated copy.
function _model_occ_revol_point!(sw,kernel,src_id,p,n,a,Ma,loca,spec,
                                 caller)
    eid=_model_occ_revol_rim!(sw,kernel,src_id,p,n,a,Ma,loca,spec,caller)
    top_vid=a<_OCC_TWO_PI ?
        _model_occ_revol_vertex!(sw,kernel,src_id,p,n,Ma,loca,spec) :
        src_id
    top_tag=_model_occ_sweep_tag_of(sw,kernel,top_vid,caller)
    return top_tag,_model_occ_sweep_tag_of(sw,kernel,eid,caller)
end

# dim-1 input: the line sweeps a meridian lateral or an annulus/disk
# plane face; laterals are the rim curves swept by its vertices. The
# reported top edge is the source at 2·Pi, its rotated copy otherwise.
function _model_occ_revol_line!(sw,kernel,src_id,p,n,a,Ma,loca,spec,
                                kind,caller)
    a_id,b_id=kernel.shapes[src_id].children
    pa=kernel.shapes[a_id].geometry.point
    pb=kernel.shapes[b_id].geometry.point
    ra,_,_=_model_occ_revol_rh(pa,p,n)
    rb,_,_=_model_occ_revol_rh(pb,p,n)
    full=a==_OCC_TWO_PI
    if kind===:annulus
        # Standalone swept-face normal follows the profile's radial
        # direction (upstream's sweep convention): outward radial
        # profile -> +n, inward -> -n.
        face_id=_model_occ_revol_annulus!(
            sw,kernel,src_id,p,n,a,Ma,loca,spec,sign(rb-ra),caller)
    else
        face_id=_model_occ_revol_meridian!(
            sw,kernel,src_id,p,n,a,Ma,loca,spec,caller)
    end
    top_eid=full ? src_id :
        _model_occ_revol_top_edge!(sw,kernel,src_id,p,n,Ma,loca,spec,
                                   caller)
    top_tag=_model_occ_sweep_tag_of(sw,kernel,top_eid,caller)
    face_tag=_model_occ_sweep_tag_of(sw,kernel,face_id,caller)
    laterals=Int[]
    for vid in (a_id,b_id)
        eid=_model_occ_revol_rim!(sw,kernel,vid,p,n,a,Ma,loca,spec,caller)
        eid!=0 && push!(laterals,
            _model_occ_sweep_tag_of(sw,kernel,eid,caller))
    end
    return top_tag,face_tag,laterals
end

# dim-2 disk input: a full-circle profile in an axial plane at radius
# rho>r_disk sweeps the compact torus publication. The producer's solid
# answers as body, its single torus face as the one lateral. The rim
# vertex and the seam meridian reuse the profile's own vertex and circle
# — the native Layers strip then welds bitwise onto the source entities'
# cells, the same convention as the compact disk prism. `sw` records
# carry the mesh-transport links: the rim parallel through the seam
# takes the point-generatrix link (via `rmap`), the torus face links the
# profile circle and the solid links the source face.
function _model_occ_revol_disk!(m::GeoModel,kernel,sw,face_id,p,n,a,
                                Ma,loca,spec,caller)
    face=kernel.shapes[face_id]
    wire=kernel.shapes[only(face.children)]
    edge_id=only(wire.geometry.uses)[1]
    g=kernel.shapes[edge_id].geometry.geometry
    radius=g.occ===:circle ? g.r : g.rx
    _,h_c,rd=_model_occ_revol_rh(g.center,p,n)
    r1=_occ_modulus(rd)
    cen=_occ3_add(p,_occ3_scale(n,h_c))
    src_vid=kernel.shapes[edge_id].children[1]
    src_vtx=kernel.shapes[src_vid].geometry.point
    # Seam radial direction and parallel: the torus frame X points at the
    # seam azimuth so the wire legs sit at u=0/u=u_max like the
    # producer's.
    rs,hs,rd_s=_model_occ_revol_rh(src_vtx,p,n)
    Xs=_occ_dir(rd_s)
    Ys=_occ_crossed(n,Xs)
    # Poloidal angle of the seam vertex on the profile circle, measured
    # from the outer radial toward +n (the producer's meridian winds
    # about Dm=-Y at the seam azimuth).
    v_s=atan(hs-h_c,rs-r1)
    flips=_dot3(g.n,Ys)>0.0
    full=a==_OCC_TWO_PI
    eq_id=_model_occ_revol_rim!(sw,kernel,src_vid,p,n,a,Ma,loca,spec,
                                caller)
    eq_id==zero(UInt64) && throw(ArgumentError(
        "$caller: cannot revolve a disk seam vertex on the revolution axis"))
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    src_tag=tag_of(face_id);circle_tag=tag_of(edge_id)
    eq_tag=tag_of(eq_id)
    σ=flips ? -1.0 : 1.0
    v2=flips ? v_s-_OCC_TWO_PI : v_s+_OCC_TWO_PI
    shell_uses=Tuple{UInt64,Bool}[]
    tpc=Dict{Int,NamedTuple}(
        eq_tag=>(fwd=(lin2d=((0.0,v_s),(1.0,0.0)),),
                 rev=(lin2d=((0.0,v2),(1.0,0.0)),)))
    if full
        # The source circle answers both meridian legs — u=0 under one
        # traversal flag, u=2*Pi under the other.
        tpc[circle_tag]=flips ?
            (fwd=(lin2d=((0.0,v_s),(0.0,-1.0)),),
             rev=(lin2d=((_OCC_TWO_PI,v_s),(0.0,-1.0)),)) :
            (fwd=(lin2d=((_OCC_TWO_PI,v_s),(0.0,1.0)),),
             rev=(lin2d=((0.0,v_s),(0.0,1.0)),))
        torus_uses=Tuple{UInt64,Bool}[(eq_id,true),
            (edge_id,flips),(eq_id,false),(edge_id,!flips)]
        end_cap_id=zero(UInt64);end_cap_tag=0
    else
        # Partial turn: the profile circle bounds the u=0 cut alone and
        # a rotated copy bounds the u=a cut (the producer's
        # start/end-meridian pair); both columns map through the same
        # poloidal direction σ.
        rot_vid=_model_occ_revol_vertex!(sw,kernel,src_vid,p,n,Ma,loca,
                                         spec)
        rot_id=_model_occ_sweep_id!(sw,kernel)
        rot_tag=_model_occ_sweep_tag!(sw,kernel,1)
        _model_occ_sweep_bind!(sw,1,rot_tag,rot_id,
            _ModelOCCShape(1,(kind=g.occ,
                geometry=(occ=g.occ,
                    center=_occ_xform_pnt(Ma,loca,g.center),
                    n=_occ_xform_dir(Ma,g.n),
                    X=_occ_xform_dir(Ma,g.X),Y=_occ_xform_dir(Ma,g.Y),
                    (g.occ===:circle ? (;r=g.r) : (;rx=g.rx,ry=g.ry))...,
                    t0=g.t0,t1=g.t1)),UInt64[rot_vid,rot_vid]))
        tpc[circle_tag]=
            (fwd=(lin2d=((0.0,v_s),(0.0,σ)),),
             rev=(lin2d=((0.0,v_s),(0.0,σ)),))
        tpc[rot_tag]=
            (fwd=(lin2d=((a,v_s),(0.0,σ)),),
             rev=(lin2d=((a,v_s),(0.0,σ)),))
        torus_uses=Tuple{UInt64,Bool}[(eq_id,true),
            (edge_id,true),(eq_id,false),(rot_id,false)]
        # Meridian caps: the source disk bounds the u=0 cut and its
        # rotated copy the u=a cut. Cap outward normals are −(n×Xs) and
        # +rot(n×Xs) — the shell flags key on the record axes.
        M=_occ_crossed(n,Xs)
        cap_wire_id=_model_occ_sweep_id!(sw,kernel)
        cap_wire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
        _model_occ_sweep_bind!(sw,-1,cap_wire_tag,cap_wire_id,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(rot_id,false)],),
                UInt64[rot_id]))
        shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
        cap_geometry=_model_occ_plane_face_geometry(
            shape_of,tag_of,
            Vector{Tuple{UInt64,Bool}}[[(rot_id,false)]],caller)
        end_cap_id=_model_occ_sweep_id!(sw,kernel)
        end_cap_tag=_model_occ_sweep_tag!(sw,kernel,2)
        _model_occ_sweep_bind!(sw,2,end_cap_tag,end_cap_id,
            _ModelOCCShape(2,(kind=:plane,geometry=cap_geometry),
                UInt64[cap_wire_id]))
        src_axis=face.geometry.geometry.axis
        rot_M=_occ_xform_dir(Ma,M)
        push!(shell_uses,(face_id,_dot3(src_axis,M)>0.0))
        push!(shell_uses,(end_cap_id,_dot3(cap_geometry.axis,rot_M)<0.0))
        push!(sw.links,(2,end_cap_tag,(2,src_tag)))
        push!(sw.links,(1,rot_tag,(1,circle_tag)))
    end
    lwire_id=_model_occ_sweep_id!(sw,kernel)
    lwire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
    _model_occ_sweep_bind!(sw,-1,lwire_tag,lwire_id,
        _ModelOCCShape(-1,(uses=torus_uses,),
            UInt64[u[1] for u in torus_uses]))
    face_id2=_model_occ_sweep_id!(sw,kernel)
    face_tag=_model_occ_sweep_tag!(sw,kernel,2)
    _model_occ_sweep_bind!(sw,2,face_tag,face_id2,
        _ModelOCCShape(2,(kind=:torus,
            geometry=(occ=:torus,center=cen,axis=n,X=Xs,Y=Ys,
                      r1=r1,r2=radius,angle=a,pcurves=tpc)),
            UInt64[lwire_id]))
    shell_id=_model_occ_sweep_id!(sw,kernel)
    shell_tag=_model_occ_sweep_tag!(sw,kernel,-2)
    _model_occ_sweep_bind!(sw,-2,shell_tag,shell_id,
        _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[
            (face_id2,false),shell_uses...],),
            UInt64[face_id2,first.(shell_uses)...]))
    solid_id=_model_occ_sweep_id!(sw,kernel)
    solid_tag=_model_occ_sweep_tag!(sw,kernel,3)
    _model_occ_sweep_bind!(sw,3,solid_tag,solid_id,
        _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id]))
    push!(sw.links,(2,face_tag,(1,circle_tag)))
    push!(sw.links,(3,solid_tag,(2,src_tag)))
    top_tag=full ? src_tag : end_cap_tag
    return top_tag,solid_tag,[face_tag]
end

# dim-2 input: a single-wire all-line plane face lying in an axial
# plane sweeps a revolved solid — every boundary edge contributes its
# swept surface (on-axis edges collapse), the shell's use flags carry
# the outward convention measured by `σ·axis`/`_occ_surface_normal`,
# and profiles that land exactly on a primitive layout carry the
# compact encoding for the measure/sync consumers.
function _model_occ_revol_planar!(sw,kernel,src_id,p,n,a,Ma,loca,spec,
                                  caller)
    source=kernel.shapes[src_id]
    g=source.geometry.geometry
    face_axis=g.axis
    wire=kernel.shapes[only(source.children)]
    full=a==_OCC_TWO_PI
    laterals=UInt64[]
    shell_uses=Tuple{UInt64,Bool}[]
    cyl_faces=cone_faces=0
    caps=Tuple{Float64,Float64}[]
    verts=Tuple{Float64,Float64}[]
    rot_wire=Tuple{UInt64,Bool}[]
    rhat=nothing
    for (edge_id,rev) in wire.geometry.uses
        a_id,b_id=kernel.shapes[edge_id].children
        pa=kernel.shapes[a_id].geometry.point
        pb=kernel.shapes[b_id].geometry.point
        ra,ha,rd_a=_model_occ_revol_rh(pa,p,n)
        rb,hb,rd_b=_model_occ_revol_rh(pb,p,n)
        push!(verts,rev ? (rb,hb) : (ra,ha))
        rhat===nothing && ra>_OCC_CONFUSION &&
            (rhat=_occ_dir(rd_a))
        rhat===nothing && rb>_OCC_CONFUSION &&
            (rhat=_occ_dir(rd_b))
        len=_occ_modulus(_arc_sub(pb,pa))
        # Profile-edge outward in the meridian plane: the wire's CCW
        # traversal (material left, like the extrusion caps) puts the
        # outward on the right of the traversal direction.
        tangent=_arc_sub(pb,pa)
        udir=_occ_dir(rev ? _occ3_scale(tangent,-1.0) : tangent)
        o3=_occ_crossed(udir,face_axis)
        dh=abs(hb-ha)
        full || push!(rot_wire,(_model_occ_revol_top_edge!(
            sw,kernel,edge_id,p,n,Ma,loca,spec,caller),rev))
        if ra<=_OCC_CONFUSION && rb<=_OCC_CONFUSION
            continue  # on-axis profile edge sweeps nothing
        elseif dh<=_OCC_CONFUSION*len
            fid=_model_occ_revol_annulus!(
                sw,kernel,edge_id,p,n,a,Ma,loca,spec,
                sign(_dot3(o3,n)),caller)
            rec=_model_occ_sweep_shape_of(sw,kernel,fid).geometry
            rev_flag=_dot3(rec.geometry.axis,o3)<0
            push!(shell_uses,(fid,rev_flag))
            push!(caps,((ha+hb)/2,min(ra,rb)))
        else
            fid=_model_occ_revol_meridian!(
                sw,kernel,edge_id,p,n,a,Ma,loca,spec,caller)
            rec=_model_occ_sweep_shape_of(sw,kernel,fid).geometry
            nv=_occ_surface_normal(rec.geometry,0.0,
                                   _occ_modulus(_arc_sub(pb,pa))/2)
            push!(shell_uses,(fid,_dot3(nv,o3)<0))
            rec.kind===:cylinder ? (cyl_faces+=1) : (cone_faces+=1)
        end
        push!(laterals,fid)
    end
    rhat===nothing && throw(ArgumentError(
        "$caller: revolve profile lies on the revolution axis"))
    top_tag=src_tag=
        _model_occ_sweep_tag_of(sw,kernel,src_id,caller)
    if !full
        # Meridian walls: the source face bounds the u=0 cut, its
        # rotated copy bounds u=a. The wall outward normals are
        # ∓(n×rhat) — flag keys on the record axes like the other
        # shell uses.
        shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
        tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
        rot_geometry=_model_occ_plane_face_geometry(
            shape_of,tag_of,[rot_wire],caller)
        rot_wire_id=_model_occ_sweep_id!(sw,kernel)
        rot_wire_tag=_model_occ_sweep_tag!(sw,kernel,-1)
        _model_occ_sweep_bind!(sw,-1,rot_wire_tag,rot_wire_id,
            _ModelOCCShape(-1,(uses=rot_wire,),
                UInt64[u[1] for u in rot_wire]))
        rot_face_id=_model_occ_sweep_id!(sw,kernel)
        rot_face_tag=_model_occ_sweep_tag!(sw,kernel,2)
        _model_occ_sweep_bind!(sw,2,rot_face_tag,rot_face_id,
            _ModelOCCShape(2,(kind=:plane,geometry=rot_geometry),
                UInt64[rot_wire_id]))
        M=_occ_crossed(n,rhat)
        push!(shell_uses,(src_id,_dot3(face_axis,M)>0.0))
        push!(shell_uses,(rot_face_id,
            _dot3(rot_geometry.axis,_occ_xform_dir(Ma,M))<0.0))
        push!(sw.links,(2,rot_face_tag,(2,src_tag)))
        top_tag=rot_face_tag
    end
    shell_id=_model_occ_sweep_id!(sw,kernel)
    shell_tag=_model_occ_sweep_tag!(sw,kernel,-2)
    _model_occ_sweep_bind!(sw,-2,shell_tag,shell_id,
        _ModelOCCShape(-2,(uses=shell_uses,),
                      UInt64[u[1] for u in shell_uses]))
    solid_geometry=(kind=:solid,)
    # Compact encodings when the swept boundary is exactly a primitive:
    # one cylinder + two axis-centered disk caps, or one cone + its
    # disk cap(s). Full turns carry the `primitive` record; partial
    # turns carry the `sector` wedge record the partial-sweep
    # primitives already measure through.
    # Each entry in `caps` is (axial h, inner radius); a disk cap's
    # inner radius is ~0 and its h must sit at a quadric end.
    cap_is_disk(c)=c[2]<=_OCC_CONFUSION
    cap_at(rec,c)=begin
        h_lo=_dot3(_arc_sub(rec.center,p),n)
        abs(c[1]-h_lo)<=_OCC_CONFUSION ||
            abs(c[1]-(h_lo+rec.height))<=_OCC_CONFUSION
    end
    shape_kind(fid)=
        _model_occ_sweep_shape_of(sw,kernel,fid).geometry.kind
    if cyl_faces==1 && cone_faces==0 && length(caps)==2 &&
            all(cap_is_disk,caps)
        rec=_model_occ_sweep_shape_of(sw,kernel,only(
            fid for (fid,_) in shell_uses
            if shape_kind(fid)===:cylinder)).geometry.geometry
        all(c->cap_at(rec,c),caps) &&
            (solid_geometry=full ?
                (kind=:solid,
                 primitive=(occ=:cylinder,
                     record=(center=rec.center,
                             axis=_occ3_scale(rec.axis,rec.height),
                             radius=rec.radius,height=rec.height))) :
                (kind=:solid,
                 sector=(occ=:cylinder,center=rec.center,
                         axis=_occ3_scale(rec.axis,rec.height),
                         X=rec.X,Y=rec.Y,radius=rec.radius,
                         height=rec.height,angle=a)))
    elseif cone_faces==1 && cyl_faces==0 &&
            1<=length(caps)<=2 && all(cap_is_disk,caps)
        rec=_model_occ_sweep_shape_of(sw,kernel,only(
            fid for (fid,_) in shell_uses
            if shape_kind(fid)===:cone)).geometry.geometry
        all(c->cap_at(rec,c),caps) &&
            (solid_geometry=full ?
                (kind=:solid,
                 primitive=(occ=:cone,
                     record=(center=rec.center,
                             axis=_occ3_scale(rec.axis,rec.height),
                             r1=rec.r1,r2=rec.r2,height=rec.height))) :
                (kind=:solid,
                 sector=(occ=:cone,center=rec.center,
                         axis=_occ3_scale(rec.axis,rec.height),
                         X=rec.X,Y=rec.Y,r1=rec.r1,r2=rec.r2,
                         height=rec.height,angle=a)))
    end
    # Every revolved solid carries Pappus integrals as a fallback
    # encoding: the profile moments give the closed-form volume and
    # the centroid's axial coordinate for the measure path even when
    # the boundary isn't a compact primitive. Partial turns add the
    # sector bisector and its radial arm — the centroid rides the
    # u=a/2 meridian off axis.
    if !haskey(solid_geometry,:primitive) &&
            !haskey(solid_geometry,:sector)
        V,h_c,arm=_model_occ_revol_profile_moments(verts,a)
        V>_OCC_CONFUSION || throw(ArgumentError(
            "$caller: revolved profile has zero swept volume"))
        solid_geometry=(kind=:solid,
            primitive=(occ=:revolved,
                record=full ? (point=p,axis=n,V=V,h_c=h_c) :
                    (point=p,axis=n,V=V,h_c=h_c,angle=a,
                     bisector=_occ_dir(_occ3_add(
                         _occ3_scale(rhat,cos(a/2)),
                         _occ3_scale(_occ_crossed(n,rhat),sin(a/2)))),
                     arm=arm)))
    end
    solid_id=_model_occ_sweep_id!(sw,kernel)
    solid_tag=_model_occ_sweep_tag!(sw,kernel,3)
    _model_occ_sweep_bind!(sw,3,solid_tag,solid_id,
        _ModelOCCShape(3,solid_geometry,UInt64[shell_id]))
    push!(sw.links,(3,solid_tag,(2,src_tag)))
    lat_tags=Int[_model_occ_sweep_tag_of(sw,kernel,fid,caller)
                 for fid in laterals]
    return top_tag,solid_tag,lat_tags
end

# Phase-1 classifier for a line input: rejects non-:line records,
# fully on-axis segments and segments whose endpoints sit in different
# axial half-planes (profile crossing or leaving the revolution's
# meridian plane). Returns :annulus for radial profiles, :meridian
# otherwise (cylinder/cone decided at build time).
function _model_occ_revol_check_line(kernel,src_id,p,n,caller)
    source=kernel.shapes[src_id]
    g=source.geometry
    g.kind===:line || throw(ArgumentError(
        "$caller: revolution of a $(g.kind) profile is not yet " *
        "supported — only line curves are"))
    a_id,b_id=source.children
    pa=kernel.shapes[a_id].geometry.point
    pb=kernel.shapes[b_id].geometry.point
    ra,ha,rd_a=_model_occ_revol_rh(pa,p,n)
    rb,hb,rd_b=_model_occ_revol_rh(pb,p,n)
    len=_occ_modulus(_arc_sub(pb,pa))
    ra>_OCC_CONFUSION || rb>_OCC_CONFUSION || throw(ArgumentError(
        "$caller: cannot revolve a line lying on the revolution axis"))
    if ra>_OCC_CONFUSION && rb>_OCC_CONFUSION
        _dot3(_occ_dir(rd_a),_occ_dir(rd_b))>1.0-1.0e-10 ||
            throw(ArgumentError(
                "$caller: revolve profile is not planar with the " *
                "revolution axis (endpoints in different meridian " *
                "planes or the profile crosses the axis)"))
    end
    return abs(hb-ha)<=_OCC_CONFUSION*len ? :annulus : :meridian
end

# Phase-1 classifier for a face input: a single-wire all-line plane
# face whose plane contains the revolution axis (:planar), or a
# full-circle disk pierced by an axis lying in its plane (:disk).
# Everything else — curved wires, multi-wire faces, non-axial planes —
# rejects here so phase 2 cannot fail after binding.
function _model_occ_revol_check_face(kernel,face_id,p,n,caller)
    face=kernel.shapes[face_id]
    g=face.geometry
    g.kind===:plane || throw(ArgumentError(
        "$caller: revolution of a $(g.kind) face is not yet " *
        "supported — only planar profiles are"))
    length(face.children)==1 || throw(ArgumentError(
        "$caller: revolution of a face with interior wires is not " *
        "yet supported"))
    wire=kernel.shapes[only(face.children)]
    uses=wire.geometry.uses
    kinds=[kernel.shapes[eid].geometry.kind
           for (eid,_) in uses]
    if length(uses)==1 &&
            (kinds[1]===:circle ||
             (kinds[1]===:ellipse &&
              let e=kernel.shapes[uses[1][1]].geometry.geometry
                  abs(e.rx-e.ry)<=_OCC_CONFUSION*e.rx
              end))
        # Disk profile: the axis must lie in the disk's own plane —
        # its normal is perpendicular to n and the center's offset
        # from the axis line has no n-normal component — and the
        # disk must clear the axis (rho > radius) or the swept torus
        # self-intersects at the spindle.
        gg=kernel.shapes[uses[1][1]].geometry.geometry
        radius=gg.occ===:circle ? gg.r : gg.rx
        _,_,rd=_model_occ_revol_rh(gg.center,p,n)
        abs(_dot3(g.geometry.axis,n))<=1.0e-10 || throw(ArgumentError(
            "$caller: cannot revolve a disk whose plane does not " *
            "contain the revolution axis"))
        _occ_modulus(rd)>radius+_OCC_CONFUSION || throw(ArgumentError(
            "$caller: cannot revolve a disk that reaches or crosses " *
            "the revolution axis (spindle torus)"))
        return :disk
    end
    all(k->k===:line,kinds) || throw(ArgumentError(
        "$caller: revolution of a face bounded by non-line curves " *
        "is not yet supported"))
    # The face plane must be an axial plane: its axis is perpendicular
    # to n and every vertex's radial direction lies in the plane.
    abs(_dot3(g.geometry.axis,n))<=1.0e-10 || throw(ArgumentError(
        "$caller: revolve profile plane does not contain the " *
        "revolution axis"))
    for (edge_id,_) in uses
        edge=kernel.shapes[edge_id]
        a_id,b_id=edge.children
        ra,_,rd_a=_model_occ_revol_rh(
            kernel.shapes[a_id].geometry.point,p,n)
        rb,_,rd_b=_model_occ_revol_rh(
            kernel.shapes[b_id].geometry.point,p,n)
        (ra>_OCC_CONFUSION ?
            abs(_dot3(g.geometry.axis,_occ_dir(rd_a))) : 0.0)<=1.0e-10 &&
        (rb>_OCC_CONFUSION ?
            abs(_dot3(g.geometry.axis,_occ_dir(rd_b))) : 0.0)<=1.0e-10 ||
            throw(ArgumentError(
                "$caller: revolve profile is not planar with the " *
                "revolution axis"))
        if ra>_OCC_CONFUSION && rb>_OCC_CONFUSION
            _dot3(_occ_dir(rd_a),_occ_dir(rd_b))>1.0-1.0e-10 ||
                throw(ArgumentError(
                    "$caller: revolve profile crosses the " *
                    "revolution axis"))
        end
    end
    return :planar
end

# Descendant vertex ids of a sweep input (dim 0..2) — used by the Layers
# transport gate: an on-axis profile vertex collapses its swept column,
# which the shared prism-collapse emission does not yet support.
function _model_occ_sweep_vertex_ids(kernel,id,dim)
    dim==0 && return UInt64[id]
    dim==1 && return collect(kernel.shapes[id].children)
    out=UInt64[]
    for wire_id in kernel.shapes[id].children,
        (edge_id,_) in kernel.shapes[wire_id].geometry.uses
        append!(out,kernel.shapes[edge_id].children)
    end
    return unique!(out)
end

function _model_occ_revolve(m::GeoModel,dim_tags,x,y,z,ax,ay,az,angle;
                            num_elements=Int[],heights=Float64[],
                            recombine=false)
    caller="OCC.revolve"
    params=_model_occ_sweep_mesh_params(num_elements,heights,recombine,caller)
    p=_finite3(x,y,z,caller)
    d=_finite3(ax,ay,az,caller)
    _occ_modulus(d)>0 || throw(ArgumentError(
        "$caller: revolution axis direction must have positive length"))
    n=_occ_dir(d)
    a=_finite_scalar(angle,caller,"angle")
    a==0.0 && throw(ArgumentError("$caller: angle must be nonzero"))
    # Revolving by -a about n is the same as a about -n.
    a<0.0 && (a=-a; n=_occ3_scale(n,-1.0))
    a<=_OCC_TWO_PI || throw(ArgumentError(
        "$caller: angle magnitude must not exceed 2*Pi"))
    full=a==_OCC_TWO_PI
    # The sweep-end rotation — partial turns place a rotated copy of
    # every profile entity at u=a (rim arcs, rotated edges and faces);
    # at 2*Pi the copies coincide with the sources and are skipped.
    # `Ma/loca` is the gp_Trsf for records; `spec` is the Layers
    # transport's step-chain transform — vertex positions take the
    # latter so they weld bitwise onto the emitted mesh rows.
    Ma,loca=_occ_trsf_rot(p,n,a)
    spec=_extrude_spec_rotate(n,p,a,caller)
    kernel=m.occ_kernel
    # Phase 1: validate and classify every input before any bind, like
    # the prism sweep — a rejection leaves the registry untouched.
    classified=Tuple{Int,UInt64,Symbol}[]
    for entry in dim_tags
        dim=_model_occ_dimension(first(entry),caller)
        0<=dim<=3 || throw(ArgumentError(
            "$caller: Wrong input dimension in revolve"))
        tag=_model_occ_requested_tag(last(entry),caller)
        id=get(kernel.bindings,(dim,tag),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $dim " *
            "with tag $tag"))
        kind=if dim==0
            rho,_,_=_model_occ_revol_rh(
                kernel.shapes[id].geometry.point,p,n)
            rho>_OCC_CONFUSION || throw(ArgumentError(
                "$caller: cannot revolve a point on the revolution axis"))
            :point
        elseif dim==1
            _model_occ_revol_check_line(kernel,id,p,n,caller)
        elseif dim==2
            _model_occ_revol_check_face(kernel,id,p,n,caller)
        else
            throw(ArgumentError(
                "$caller: revolution of OCC volumes is not yet supported"))
        end
        if params!==nothing && dim>0
            for vid in _model_occ_sweep_vertex_ids(kernel,id,dim)
                rho,_,_=_model_occ_revol_rh(
                    kernel.shapes[vid].geometry.point,p,n)
                rho>_OCC_CONFUSION || throw(ArgumentError(
                    "$caller: extrusion meshing of axis-touching " *
                    "profiles is not yet supported"))
            end
        end
        push!(classified,(dim,id,kind))
    end
    # A full-turn strip whose last level welds onto the source row needs
    # at least three subdivisions: with two, each arc's second gap lands
    # on the same four rows as the first — coincident cells that no mesh
    # can represent. A pure point generatrix only sweeps a ring curve —
    # two arcs still tile it, but one is a degenerate self-loop. Partial
    # turns end on the rotated copy instead of welding — one layer
    # already tiles the strip.
    if params!==nothing
        need=full ? (any(dim>0 for (dim,_,_) in classified) ? 3 : 2) : 1
        sum(params.layers)<need && throw(ArgumentError(
            "$caller: a $(full ? "full-turn" : "partial") revolve " *
            "needs at least $need total numElements subdivisions"))
    end
    sw=_model_occ_sweep_setup(n,caller,kernel)
    out=Tuple{Int,Int}[]
    for dim in (0,1,2)
        for (edim,id,kind) in classified
            edim==dim || continue
            if dim==0
                top,body=_model_occ_revol_point!(
                    sw,kernel,id,p,n,a,Ma,loca,spec,caller)
                push!(out,(0,top));push!(out,(1,body))
            elseif dim==1
                top,body,laterals=_model_occ_revol_line!(
                    sw,kernel,id,p,n,a,Ma,loca,spec,kind,caller)
                push!(out,(1,top));push!(out,(2,body))
                for tag in laterals
                    push!(out,(1,tag))
                end
            else
                top,body,laterals=kind===:disk ?
                    _model_occ_revol_disk!(m,kernel,sw,id,p,n,a,Ma,loca,
                        spec,caller) :
                    _model_occ_revol_planar!(sw,kernel,id,p,n,a,Ma,loca,
                        spec,caller)
                push!(out,(2,top));push!(out,(3,body))
                for tag in laterals
                    push!(out,(2,tag))
                end
            end
        end
    end
    for (dimension,tag,id,shape) in sw.pending
        _model_occ_bind_shape!(kernel,dimension,tag,id,shape)
    end
    isempty(sw.pending) ||
        (kernel.next_shape=max(kernel.next_shape,sw.next_id))
    params===nothing || _model_occ_sweep_transport!(
        m,kernel,sw,params,spec)
    return out
end

# ---------------------------------------------------------------------------
# `model.occ.addThickSections` — thicken a bound OCC surface into solids
# along its normal. Upstream gmsh 4.15.2 publishes no endpoint of this
# name (its OCC set ends at `addThruSections`/`addThickSolid`); the
# contract here follows the OCCT `BRepOffset` thickening recipe — a face
# offset copy at each listed distance bounds one *thick section*, so
# `offsets` reads as the signed section boundaries measured along the
# face normal with an implicit leading 0 at the source surface. Each of
# the `n` slabs [o_{i-1}, o_i] is the exact prism sweep of the lower
# section face along ±n̂ — the same records `_model_occ_extrude_planar!`
# certifies — and adjacent slabs share the intermediate offset face,
# wires, edges and vertices, matching the Layers-style conforming
# decomposition. The offset list must be strictly monotone so the
# sections cannot overlap or degenerate. The bounded scope mirrors the
# prism sweep's two exact profile classes: all-line boundaries (ruled
# side walls) and the single full-circle profile (compact cylinder
# wall); every other surface family rejects loudly. `tag` binds the sole
# resulting solid like upstream's `_multiBind` tag rule — it applies
# only when a single entity is produced.

# Retag a bound shape in place: upstream's `_multiBind` honors `tag`
# only for a sole result entity, which the sweep cursor already
# allocated by the time the caller knows the answer. The entity is a
# leaf (no bound parent references it), so moving its tag touches only
# the three lookup tables plus the per-dim maximum.
function _model_occ_retag_shape!(kernel,dimension,from,to,caller)
    from==to && return nothing
    id=get(kernel.bindings,(dimension,from),nothing)
    id===nothing && throw(ArgumentError(
        "$caller: OCC entity ($dimension,$from) is unbound"))
    haskey(kernel.bindings,(dimension,to)) && throw(ArgumentError(
        "$caller: OCC entity ($dimension,$to) already exists"))
    _model_occ_entity_map_unbind!(kernel.entity_maps[dimension+3],from)
    delete!(kernel.bindings,(dimension,from))
    kernel.bindings[(dimension,to)]=id
    _model_occ_entity_map_bind!(kernel.entity_maps[dimension+3],to)
    kernel.shape_tags[id]=(dimension,to)
    kernel.maxima[dimension+3]=maximum(
        (t for (d,t) in keys(kernel.bindings) if d==dimension);init=0)
    kernel.changed=true
    return nothing
end

function _model_occ_thick_sections(m::GeoModel,tag,surface_tag,offsets;
                                   caller="OCC.addThickSections")
    kernel=m.occ_kernel
    requested=_model_occ_requested_tag(tag,caller)
    # A negative tag flips the face normal — the same reversed-use
    # convention `add_wire` applies to negative curve tags.
    st=_model_occ_requested_tag(surface_tag,caller)
    st==0 && throw(ArgumentError("$caller: OCC Surface[0] is invalid"))
    flip=st<0
    fid=get(kernel.bindings,(2,abs(st)),nothing)
    fid===nothing && throw(ArgumentError(
        "$caller: unknown OCC Surface[$(abs(st))]"))
    fshape=kernel.shapes[fid]
    fg=fshape.geometry
    fg.kind===:plane || throw(ArgumentError(
        "$caller: addThickSections needs a planar surface " *
        "(got kind=$(fg.kind))"))
    offs=Float64[_finite_scalar(o,caller,"offset") for o in offsets]
    isempty(offs) && throw(ArgumentError(
        "$caller: addThickSections requires at least one offset"))
    # The source surface is the implicit 0-boundary; every listed offset
    # bounds one slab, so the sequence must march strictly in one
    # direction along the normal.
    sign_of=offs[1]>0.0 ? 1.0 : offs[1]<0.0 ? -1.0 : 0.0
    sign_of==0.0 && throw(ArgumentError(
        "$caller: the first offset must be non-zero — the source " *
        "surface is the implicit zero-thickness section"))
    bounds=Vector{Float64}(undef,length(offs)+1)
    bounds[1]=0.0
    for i in eachindex(offs)
        bounds[i+1]=offs[i]
        step=bounds[i+1]-bounds[i]
        step*sign_of>0.0 || throw(ArgumentError(
            "$caller: offsets must be strictly monotone — the thick " *
            "sections would overlap or degenerate"))
    end
    # Profile admissibility — the same two classes the prism sweep
    # certifies exactly (the sweep direction below is always the face
    # normal, so the direction gates cannot fire).
    edges=0;all_lines=true
    for wid in fshape.children,
            (eid,_) in kernel.shapes[wid].geometry.uses
        edges+=1
        kernel.shapes[eid].geometry.kind===:line || (all_lines=false)
    end
    if !all_lines
        edges==1 || throw(ArgumentError(
            "$caller: addThickSections supports all-line boundaries " *
            "or a single full-circle profile only"))
        eid=only(kernel.shapes[only(fshape.children)].geometry.uses)[1]
        eg=kernel.shapes[eid].geometry.geometry
        if eg.occ===:circle || (eg.occ===:ellipse &&
                abs(eg.rx-eg.ry)<=_OCC_CONFUSION*eg.rx)
            eg.t1-eg.t0==_OCC_TWO_PI || throw(ArgumentError(
                "$caller: addThickSections supports full-circle " *
                "profiles only"))
        else
            throw(ArgumentError(
                "$caller: addThickSections supports circular " *
                "profiles only"))
        end
    end
    # The face's canonical normal — `reversed` marks records baked with
    # the opposite intent, the same convention `_occ_surface_normal`
    # applies.
    n_eff=get(fg.geometry,:reversed,false) ? (.-fg.geometry.axis) :
                                           fg.geometry.axis
    flip && (n_eff=(.-n_eff))
    ndir=_occ_dir(n_eff)
    nsolids=length(offs)
    # Upstream `_multiBind` binds `tag` to the first entity only and
    # errors when a second entity of that dimension would share it —
    # reject before any bind so the failure stays atomic.
    requested>0 && nsolids>1 && throw(ArgumentError(
        "$caller: cannot bind $nsolids volumes to a single tag " *
        "$requested — use tag=-1 for multi-slab thickening"))
    requested>0 && haskey(kernel.bindings,(3,requested)) &&
        throw(ArgumentError(
            "$caller: OCC Volume[$requested] already exists"))
    # Each interval is the lower face swept by the signed step along
    # n̂; the previous slab's cap is the next section's source face, so
    # the intermediate offset faces stay shared entities.
    solid_tags=Int[]
    prev_tag=abs(st)
    for i in eachindex(offs)
        step=bounds[i+1]-bounds[i]
        d=_occ3_scale(ndir,step)
        seg=_model_occ_extrude(m,[(2,prev_tag)],d[1],d[2],d[3];
                               caller=caller)
        # dim-2 extrusion reports [(2,cap),(3,solid),(2,laterals…)].
        cap_tag=seg[1][2]; solid_tag=seg[2][2]
        push!(solid_tags,solid_tag)
        prev_tag=cap_tag
    end
    if requested>0 && nsolids==1
        _model_occ_retag_shape!(kernel,3,solid_tags[1],requested,caller)
        solid_tags[1]=requested
    end
    return Tuple{Int,Int}[(3,t) for t in solid_tags]
end
