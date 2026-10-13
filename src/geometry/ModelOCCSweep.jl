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
end

function _model_occ_sweep_setup(d,caller,kernel)
    return _ModelOCCSweep(d,caller,
        Tuple{Int,Int,UInt64,_ModelOCCShape}[],
        Dict{UInt64,_ModelOCCShape}(),Dict{UInt64,Int}(),
        collect(kernel.maxima),kernel.next_shape,
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
