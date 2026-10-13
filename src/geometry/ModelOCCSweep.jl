# `gmsh.model.occ.extrude` — `OCC_Internals::extrude`: a
# `BRepPrimAPI_MakePrism`-style translation sweep per input dimension.
# The bounded implementation covers the all-linear profiles exactly:
# points sweep to line edges, line edges to ruled parallelogram plane
# faces, and line-bounded plane faces to prism solids; a single-circle
# disk swept along its own normal reuses the compact cylinder
# publication (upstream shares the source face as the bottom cap — the
# producer builds a geometrically identical cap instead, so the source
# stays bound as an independent entity). Curved profiles, solid inputs,
# in-plane directions and the extrusion-meshing options
# (`numElements`/`heights`/`recombine`) reject loudly. Output order
# follows upstream's `getReturnedShapes` reorder — per input
# `[top (dim-1), body (dim), laterals …]` — with dim-0 inputs emitting
# just the swept edge, all grouped by ascending dimension like
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
end

function _model_occ_sweep_setup(d,caller,kernel)
    return _ModelOCCSweep(d,caller,
        Tuple{Int,Int,UInt64,_ModelOCCShape}[],
        Dict{UInt64,_ModelOCCShape}(),Dict{UInt64,Int}(),
        collect(kernel.maxima),kernel.next_shape,
        Dict{UInt64,UInt64}(),Dict{UInt64,UInt64}(),
        Dict{UInt64,UInt64}(),Dict{UInt64,UInt64}(),
        Dict{UInt64,UInt64}())
end

# `next_shape` only advances at the flush, so identity allocation is
# sweep-scoped; it still resyncs when an eager producer (the cylinder
# path) bumped `kernel.next_shape` mid-call, and the flush takes `max`
# against those eagerly-bound shapes.
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
    return sw.tags[top],face_tag,[sw.tags[la],sw.tags[lb]]
end

# Single-circle disk swept along its normal → the compact cylinder
# publication. The producer's +d̂ cap answers as "top", its lateral as
# the single lateral face.
function _model_occ_extrude_disk!(m::GeoModel,kernel,face_id,d,caller)
    face=kernel.shapes[face_id]
    wire=kernel.shapes[only(face.children)]
    edge_id=only(wire.geometry.uses)[1]
    g=kernel.shapes[edge_id].geometry.geometry
    # Phase 1 already certified the full-circle, axis-aligned profile.
    radius=g.occ===:circle ? g.r : g.rx
    solid_tag=_model_occ_add_cylinder!(m,g.center[1],g.center[2],
        g.center[3],d[1],d[2],d[3],radius)
    solid_id=kernel.bindings[3,solid_tag]
    top_center=_occ3_add(g.center,d)
    top_tag=0;lateral_tag=0
    for (fid,_) in kernel.shapes[only(
            kernel.shapes[solid_id].children)].geometry.uses
        fg=kernel.shapes[fid].geometry
        if fg.kind===:cylinder
            lateral_tag=_model_occ_shape_tag(kernel,fid,caller)
        elseif fg.kind===:plane &&
                _occ_modulus(_arc_sub(fg.geometry.center,top_center))<
                _OCC_CONFUSION*_occ_modulus(d)
            top_tag=_model_occ_shape_tag(kernel,fid,caller)
        end
    end
    (top_tag>0 && lateral_tag>0) || throw(ArgumentError(
        "$caller: internal cylinder sweep publication is inconsistent"))
    return top_tag,solid_tag,[lateral_tag]
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
    return cap_tag,solid_tag,[tag for (_,tag) in laterals]
end

function _model_occ_extrude(m::GeoModel,dim_tags,dx,dy,dz;
                            num_elements=Int[],heights=Float64[],
                            recombine=false)
    caller="OCC.extrude"
    isempty(num_elements) && isempty(heights) && !recombine ||
        throw(ArgumentError(
            "$caller: extrusion meshing (numElements, heights, " *
            "recombine) is not yet supported"))
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
                    _model_occ_extrude_disk!(m,kernel,id,d,caller) :
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
    return out
end

# ---------------------------------------------------------------------------
# `gmsh.model.occ.revolve` — `OCC_Internals::revolve`: a
# `BRepPrimAPI_MakeRevol`-style rotational sweep per input dimension.
# The bounded implementation covers full-turn (angle == 2·Pi) sweeps of
# axial-plane profiles exactly: points sweep to rim circles, line edges
# to cylinder/cone laterals or annulus/disk plane faces, single-wire
# all-line plane faces to revolved solids (cylinder/cone compact
# encodings when the profile produces a primitive), and a full-circle
# disk pierced by the axis to the compact torus publication. Partial
# angles, non-axial-plane profiles, axis-crossing profiles, curved
# profiles other than disks, solid inputs, multi-wire faces and the
# extrusion-meshing options reject loudly. At 2·Pi the rotated "top"
# coincides with the profile — the source entity answers as `top` —
# and each boundary edge's swept surface answers as a lateral, matching
# upstream's `getReturnedShapes` reorder. Rim circles, meridian faces
# and degenerate apex edges deduplicate across shared profile entities
# through `rmap`/`fmap` like the prism sweep's maps.

# (rho,h,radial vector) of a point in the revolve frame: h is the axial
# coordinate along unit axis n, rho the distance to the axis line.
function _model_occ_revol_rh(pt,p,n)
    w=_arc_sub(pt,p);h=_dot3(w,n)
    r=_occ3_add(w,_occ3_scale(n,-h))
    return _occ_modulus(r),h,r
end

# Rim circle swept by a profile vertex, deduplicated through `rmap`;
# an on-axis vertex maps to the zero sentinel and is skipped by the
# callers (it sweeps to itself, not a curve).
function _model_occ_revol_rim!(sw,kernel,vid,p,n,caller)
    haskey(sw.rmap,vid) && return sw.rmap[vid]
    v=kernel.shapes[vid].geometry.point
    rho,h,r=_model_occ_revol_rh(v,p,n)
    if rho<=_OCC_CONFUSION
        sw.rmap[vid]=zero(UInt64)
        return zero(UInt64)
    end
    Xc,Yc=_occ_ax2_setx(n,_occ_dir(r))
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:circle,
            geometry=(occ=:circle,
                      center=_occ3_add(p,_occ3_scale(n,h)),
                      n=_occ_crossed(Xc,Yc),X=Xc,Y=Yc,r=rho,
                      t0=0.0,t1=_OCC_TWO_PI)),UInt64[vid,vid]))
    sw.rmap[vid]=id
    return id
end

# Degenerate rim edge for an apex vertex (rho==0) on a cone lateral —
# the producer convention binds the collapsed meridian limit to the
# (shared) profile vertex.
function _model_occ_revol_apex!(sw,kernel,vid,face_geometry,caller)
    id=_model_occ_sweep_id!(sw,kernel)
    tag=_model_occ_sweep_tag!(sw,kernel,1)
    _model_occ_sweep_bind!(sw,1,tag,id,
        _ModelOCCShape(1,(kind=:degenerate,
            geometry=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI,
                      face=face_geometry)),UInt64[vid,vid]))
    return id
end

# lin2d origin shifted so the pcurve evaluates correctly on an edge
# whose own parameter starts at t0 (lin2d evals at absolute t).
_occ_pcurve_shift(o,d,t0)=(o[1]-d[1]*t0,o[2]-d[2]*t0)

# Pappus moments of the revolved profile polygon: over fan triangles
# from the first vertex, ∫∫ρ dA drives the volume (2π·∫∫ρ dA) and
# the ρ-weighted axial moment ∫∫h·ρ dA drives the centroid's axial
# coordinate — for a linear function the triangle mean is the vertex
# average, and for the bilinear h·ρ the mean is
# [(Σh)(Σρ)+Σhᵢρᵢ]/12. Signed terms cancel in the ratio; V takes |J0|.
function _model_occ_revol_profile_moments(verts)
    r0,h0=verts[1]
    J0=J1=0.0
    for i in 2:length(verts)-1
        r1,h1=verts[i];r2,h2=verts[i+1]
        A2=(r1-r0)*(h2-h0)-(r2-r0)*(h1-h0)
        sr=r0+r1+r2;sh=h0+h1+h2
        J0+=A2*sr/6
        J1+=A2*(sh*sr+r0*h0+r1*h1+r2*h2)/24
    end
    return abs(2π*J0),J0==0.0 ? 0.0 : J1/J0
end

# Meridian face swept by a profile line edge: the canonical
# cylinder/cone lateral record (ascending meridian, X at the seam's
# radial direction) with the source edge reused as seam and rim
# circles (or a degenerate apex edge) per endpoint. `fmap` deduplicates
# per profile edge. Returns the swept face's shape id.
function _model_occ_revol_meridian!(sw,kernel,src_id,p,n,caller)
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
    seam_tag=_model_occ_sweep_tag_of(sw,kernel,src_id,caller)
    # The seam pcurve maps the source edge's own parameter to the
    # meridian v: ascending profiles run 0->len, descending len->0.
    gspan=g.t1-g.t0
    sc=(a_id==lo_vid ? 1.0 : -1.0)*len/gspan
    o_seam=_occ_pcurve_shift(a_id==lo_vid ? (0.0,0.0) : (0.0,len),
                           (0.0,sc),g.t0)
    c_lo=c_hi=zero(UInt64)
    lat_pcurves=Dict{Int,NamedTuple}(
        seam_tag=>(fwd=(lin2d=((_OCC_TWO_PI,o_seam[2]),(0.0,sc)),),
                   rev=(lin2d=((0.0,o_seam[2]),(0.0,sc)),)))
    lat_geometry=nothing
    if cylinder
        lat_geometry=(occ=:cylinder,center=cen,axis=n,X=Xc,Y=Yc,
                      radius=rho_lo,height=span,pcurves=lat_pcurves)
    else
        lat_geometry=(occ=:cone,center=cen,axis=n,X=Xc,Y=Yc,
                      r1=rho_lo,r2=rho_hi,height=span,pcurves=lat_pcurves)
    end
    for (vid,vpos,is_hi) in ((lo_vid,0.0,false),(hi_vid,len,true))
        rho=vid==lo_vid ? rho_lo : rho_hi
        if rho<=_OCC_CONFUSION
            eid=_model_occ_revol_apex!(sw,kernel,vid,lat_geometry,caller)
        else
            eid=_model_occ_revol_rim!(sw,kernel,vid,p,n,caller)
        end
        etag=_model_occ_sweep_tag_of(sw,kernel,eid,caller)
        lat_pcurves[etag]=(fwd=(lin2d=((0.0,vpos),(1.0,0.0)),),
                           rev=nothing)
        is_hi ? (c_hi=eid) : (c_lo=eid)
    end
    # Producer lateral wire: [rim_hi rev, seam leg, rim_lo fwd, seam
    # leg]. The seam legs' flags adapt to the source edge's own child
    # order (children run lo->hi like the producers, or hi->lo).
    wire=Tuple{UInt64,Bool}[(c_hi,true)]
    push!(wire,(src_id,a_id==lo_vid))
    push!(wire,(c_lo,false))
    push!(wire,(src_id,!(a_id==lo_vid)))
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
# coordinate h. `axis_sign` chooses the stored plane axis = ±n — the
# caller derives it from the profile direction (standalone faces) or
# the required outward normal (solid caps). Returns the face id.
function _model_occ_revol_annulus!(sw,kernel,src_id,p,n,axis_sign,caller)
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
        sw,kernel,ra>=rb ? a_id : b_id,p,n,caller)
    r_lo>_OCC_CONFUSION &&
        (c_in=_model_occ_revol_rim!(
            sw,kernel,ra>=rb ? b_id : a_id,p,n,caller))
    shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    # Wind the outer circle so the face's Newell axis lands on
    # axis_sign·n; hole wires keep the same stored orientation (the
    # children index marks them, like `add_plane_surface` users).
    flag=axis_sign<0
    uses=[Tuple{UInt64,Bool}[(c_out,flag)]]
    c_in!=0 && push!(uses,Tuple{UInt64,Bool}[(c_in,flag)])
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

# dim-0 input: the vertex sweeps a rim circle. Upstream emits
# [top vertex, swept edge]; at 2·Pi the top is the source vertex.
function _model_occ_revol_point!(sw,kernel,src_id,p,n,caller)
    eid=_model_occ_revol_rim!(sw,kernel,src_id,p,n,caller)
    src_tag=_model_occ_sweep_tag_of(sw,kernel,src_id,caller)
    return src_tag,_model_occ_sweep_tag_of(sw,kernel,eid,caller)
end

# dim-1 input: the line sweeps a meridian lateral or an annulus/disk
# plane face; laterals are the rim circles swept by its vertices.
function _model_occ_revol_line!(sw,kernel,src_id,p,n,kind,caller)
    a_id,b_id=kernel.shapes[src_id].children
    pa=kernel.shapes[a_id].geometry.point
    pb=kernel.shapes[b_id].geometry.point
    ra,_,_=_model_occ_revol_rh(pa,p,n)
    rb,_,_=_model_occ_revol_rh(pb,p,n)
    if kind===:annulus
        # Standalone swept-face normal follows the profile's radial
        # direction (upstream's sweep convention): outward radial
        # profile -> +n, inward -> -n.
        face_id=_model_occ_revol_annulus!(
            sw,kernel,src_id,p,n,sign(rb-ra),caller)
    else
        face_id=_model_occ_revol_meridian!(sw,kernel,src_id,p,n,caller)
    end
    src_tag=_model_occ_sweep_tag_of(sw,kernel,src_id,caller)
    face_tag=_model_occ_sweep_tag_of(sw,kernel,face_id,caller)
    laterals=Int[]
    for vid in (a_id,b_id)
        eid=_model_occ_revol_rim!(sw,kernel,vid,p,n,caller)
        eid!=0 && push!(laterals,
            _model_occ_sweep_tag_of(sw,kernel,eid,caller))
    end
    return src_tag,face_tag,laterals
end

# dim-2 disk input: a full-circle profile in an axial plane at radius
# rho>r_disk sweeps the compact torus publication. The producer's solid
# answers as body, its single torus face as the one lateral.
function _model_occ_revol_disk!(m::GeoModel,kernel,face_id,p,n,caller)
    face=kernel.shapes[face_id]
    wire=kernel.shapes[only(face.children)]
    edge_id=only(wire.geometry.uses)[1]
    g=kernel.shapes[edge_id].geometry.geometry
    radius=g.occ===:circle ? g.r : g.rx
    _,h_c,rd=_model_occ_revol_rh(g.center,p,n)
    r1=_occ_modulus(rd)
    cen=_occ3_add(p,_occ3_scale(n,h_c))
    solid_tag=_model_occ_add_torus!(
        m,cen[1],cen[2],cen[3],r1,radius;zAxis=collect(n))
    solid_id=kernel.bindings[3,solid_tag]
    face_tag=0
    for (fid,_) in kernel.shapes[only(
            kernel.shapes[solid_id].children)].geometry.uses
        kernel.shapes[fid].geometry.kind===:torus &&
            (face_tag=_model_occ_shape_tag(kernel,fid,caller))
    end
    face_tag>0 || throw(ArgumentError(
        "$caller: internal torus sweep publication is inconsistent"))
    src_tag=_model_occ_shape_tag(kernel,face_id,caller)
    return src_tag,solid_tag,[face_tag]
end

# dim-2 input: a single-wire all-line plane face lying in an axial
# plane sweeps a revolved solid — every boundary edge contributes its
# swept surface (on-axis edges collapse), the shell's use flags carry
# the outward convention measured by `σ·axis`/`_occ_surface_normal`,
# and profiles that land exactly on a primitive layout carry the
# compact encoding for the measure/sync consumers.
function _model_occ_revol_planar!(sw,kernel,src_id,p,n,caller)
    source=kernel.shapes[src_id]
    g=source.geometry.geometry
    face_axis=g.axis
    wire=kernel.shapes[only(source.children)]
    laterals=UInt64[]
    shell_uses=Tuple{UInt64,Bool}[]
    cyl_faces=cone_faces=0
    caps=Tuple{Float64,Float64}[]
    verts=Tuple{Float64,Float64}[]
    for (edge_id,rev) in wire.geometry.uses
        a_id,b_id=kernel.shapes[edge_id].children
        pa=kernel.shapes[a_id].geometry.point
        pb=kernel.shapes[b_id].geometry.point
        ra,ha,_=_model_occ_revol_rh(pa,p,n)
        rb,hb,_=_model_occ_revol_rh(pb,p,n)
        push!(verts,rev ? (rb,hb) : (ra,ha))
        len=_occ_modulus(_arc_sub(pb,pa))
        # Profile-edge outward in the meridian plane: the wire's CCW
        # traversal (material left, like the extrusion caps) puts the
        # outward on the right of the traversal direction.
        tangent=_arc_sub(pb,pa)
        udir=_occ_dir(rev ? _occ3_scale(tangent,-1.0) : tangent)
        o3=_occ_crossed(udir,face_axis)
        dh=abs(hb-ha)
        if ra<=_OCC_CONFUSION && rb<=_OCC_CONFUSION
            continue  # on-axis profile edge sweeps nothing
        elseif dh<=_OCC_CONFUSION*len
            fid=_model_occ_revol_annulus!(
                sw,kernel,edge_id,p,n,sign(_dot3(o3,n)),caller)
            rec=_model_occ_sweep_shape_of(sw,kernel,fid).geometry
            rev_flag=_dot3(rec.geometry.axis,o3)<0
            push!(shell_uses,(fid,rev_flag))
            push!(caps,((ha+hb)/2,min(ra,rb)))
        else
            fid=_model_occ_revol_meridian!(sw,kernel,edge_id,p,n,caller)
            rec=_model_occ_sweep_shape_of(sw,kernel,fid).geometry
            nv=_occ_surface_normal(rec.geometry,0.0,
                                   _occ_modulus(_arc_sub(pb,pa))/2)
            push!(shell_uses,(fid,_dot3(nv,o3)<0))
            rec.kind===:cylinder ? (cyl_faces+=1) : (cone_faces+=1)
        end
        push!(laterals,fid)
    end
    shell_id=_model_occ_sweep_id!(sw,kernel)
    shell_tag=_model_occ_sweep_tag!(sw,kernel,-2)
    _model_occ_sweep_bind!(sw,-2,shell_tag,shell_id,
        _ModelOCCShape(-2,(uses=shell_uses,),
                      UInt64[u[1] for u in shell_uses]))
    solid_geometry=(kind=:solid,)
    # Compact encodings when the swept boundary is exactly a primitive:
    # one cylinder + two axis-centered disk caps, or one cone + its
    # disk cap(s).
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
            (solid_geometry=(kind=:solid,
                primitive=(occ=:cylinder,
                    record=(center=rec.center,
                            axis=_occ3_scale(rec.axis,rec.height),
                            radius=rec.radius,height=rec.height))))
    elseif cone_faces==1 && cyl_faces==0 &&
            1<=length(caps)<=2 && all(cap_is_disk,caps)
        rec=_model_occ_sweep_shape_of(sw,kernel,only(
            fid for (fid,_) in shell_uses
            if shape_kind(fid)===:cone)).geometry.geometry
        all(c->cap_at(rec,c),caps) &&
            (solid_geometry=(kind=:solid,
                primitive=(occ=:cone,
                    record=(center=rec.center,
                            axis=_occ3_scale(rec.axis,rec.height),
                            r1=rec.r1,r2=rec.r2,height=rec.height))))
    end
    # Every revolved solid carries Pappus integrals as a fallback
    # encoding: the profile moments give the closed-form volume and
    # the centroid's axial coordinate for the measure path even when
    # the boundary isn't a compact primitive.
    if !haskey(solid_geometry,:primitive)
        V,h_c=_model_occ_revol_profile_moments(verts)
        V>_OCC_CONFUSION || throw(ArgumentError(
            "$caller: revolved profile has zero swept volume"))
        solid_geometry=(kind=:solid,
            primitive=(occ=:revolved,
                record=(point=p,axis=n,V=V,h_c=h_c)))
    end
    solid_id=_model_occ_sweep_id!(sw,kernel)
    solid_tag=_model_occ_sweep_tag!(sw,kernel,3)
    _model_occ_sweep_bind!(sw,3,solid_tag,solid_id,
        _ModelOCCShape(3,solid_geometry,UInt64[shell_id]))
    src_tag=_model_occ_sweep_tag_of(sw,kernel,src_id,caller)
    lat_tags=Int[_model_occ_sweep_tag_of(sw,kernel,fid,caller)
                 for fid in laterals]
    return src_tag,solid_tag,lat_tags
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

function _model_occ_revolve(m::GeoModel,dim_tags,x,y,z,ax,ay,az,angle;
                            num_elements=Int[],heights=Float64[],
                            recombine=false)
    caller="OCC.revolve"
    isempty(num_elements) && isempty(heights) && !recombine ||
        throw(ArgumentError(
            "$caller: revolve meshing (numElements, heights, " *
            "recombine) is not yet supported"))
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
    a==_OCC_TWO_PI || throw(ArgumentError(
        "$caller: partial revolve angles are not yet supported"))
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
        push!(classified,(dim,id,kind))
    end
    sw=_model_occ_sweep_setup(n,caller,kernel)
    out=Tuple{Int,Int}[]
    for dim in (0,1,2)
        for (edim,id,kind) in classified
            edim==dim || continue
            if dim==0
                top,body=_model_occ_revol_point!(sw,kernel,id,p,n,caller)
                push!(out,(0,top));push!(out,(1,body))
            elseif dim==1
                top,body,laterals=_model_occ_revol_line!(
                    sw,kernel,id,p,n,kind,caller)
                push!(out,(1,top));push!(out,(2,body))
                for tag in laterals
                    push!(out,(1,tag))
                end
            else
                top,body,laterals=kind===:disk ?
                    _model_occ_revol_disk!(m,kernel,id,p,n,caller) :
                    _model_occ_revol_planar!(sw,kernel,id,p,n,caller)
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
    return out
end
