# `gmsh.model.occ.addSurfaceLoop(surfaceTags, tag=-1, sewing=false)` — the
# `BRepPrim_Builder::MakeShell` + `AddShellFace` path: a shell over the
# bound faces, all FORWARD, with no closure check — upstream defers
# manifold validity to ShapeFix/meshing, and the emitted `m.surface_loops`
# entry follows the same deferred rule. With `sewing=true` upstream runs
# `BRepBuilderAPI_Sewing` and recursively rebinds the result: the sewed
# shell's faces are new TopoDS objects on a welded graph, so every member
# re-publishes under fresh tags. `_model_occ_sew_surface_loop!` mirrors
# that — coincident boundary vertices weld at `Precision::Confusion`,
# coincident edges (same endpoint pair + same curve locus) merge into one
# shared edge, and each member face re-publishes as a copy whose wires
# reference the welded graph. Member edges stored against the weld
# direction flip their wire flags and reparametrize their `lin2d`/`circ2d`
# pcurves onto the shared edge's parameter map; records that cannot carry
# the reversed trace refuse the weld instead of corrupting the face.
#
# `gmsh.model.occ.addVolume(shellTags, tag=-1)` —
# `BRepBuilderAPI_MakeSolid`: the first shell is the exterior, further
# shells are cavities.
#
# `gmsh.model.occ.addBox(x,y,z,dx,dy,dz,tag=-1)` — `BRepPrimAPI_MakeBox`
# between the two corner points; a zero extent rejects as "Degenerate
# box". The composed shell shares edges across faces like the sewn
# primitive, and the whole graph — 8 vertices, 12 edges, 6 wires, 6
# plane faces with `lin2d` pcurves, shell, solid — publishes atomically
# like `makeRectangle`/`makeDisk`.
#
# `gmsh.model.occ.addWedge(x,y,z,dx,dy,dz,tag=-1,ltx=0,zAxis)` —
# `BRepPrimAPI_MakeWedge`: the xz trapezoid (0,0)→(dx,0)→(ltx,dz)→(0,dz)
# in the `gp_Ax2(p,zAxis)` frame, extruded `dy` along +y. `ltx==0`
# collapses the top edge into the ridge vertex — a triangular profile,
# 5 faces rather than 6 — and `ltx==dx` reproduces the box layout. A
# zero extent rejects as "Degenerate wedge"; `zAxis` follows the same
# strict 3-vector contract as `add_torus`.

function _model_occ_add_surface_loop!(m::GeoModel,face_tags;tag=-1,
                                      sewing=false)
    caller="OCC.add_surface_loop"
    kernel=m.occ_kernel
    sewing isa Bool || throw(ArgumentError("$caller: sewing must be Bool"))
    allocated=_model_occ_next_tag(kernel,-2,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(-2,allocated)) && throw(ArgumentError(
        "$caller: OCC Surface Loop[$allocated] already exists"))
    face_ids=UInt64[]
    for raw in face_tags
        current=_model_occ_requested_tag(raw,caller)
        id=get(kernel.bindings,(2,current),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: unknown OCC Surface[$current]"))
        push!(face_ids,id)
    end
    sewing && !isempty(face_ids) &&
        return _model_occ_sew_surface_loop!(m,face_ids,allocated,caller)
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    uses=Tuple{UInt64,Bool}[(fid,false) for fid in face_ids]
    _model_occ_bind_shape!(kernel,-2,allocated,id,
        _ModelOCCShape(-2,(uses=uses,),face_ids))
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# Point on a stored edge record — used to probe the traversal direction of
# closed-edge candidates whose endpoints coincide. `nothing` for records
# with no parametric evaluator (they weld only by identical records).
function _model_occ_sew_edge_point(g,t)
    g.occ===:line && return _occ_line_point(g,t)
    g.occ===:circle && return _occ_circle_point(g,t)
    g.occ===:ellipse && return _occ_ellipse_point(g,t)
    return nothing
end

# `BRepBuilderAPI_Sewing` — weld the member faces' boundary graph:
# vertices coincide within `Precision::Confusion`; edges weld when their
# endpoint vertex groups match (either order) and their stored curves are
# the same locus — lines by endpoints, circle arcs by frame/radius/span/
# midpoint, and anything else only by identical records (so a welded
# member never needs a parameter flip the record kinds cannot express).
function _model_occ_sew_edge_coincident(g1,g2,tol)
    g1.occ===g2.occ || return false
    if g1.occ===:line
        return true
    elseif g1.occ===:circle
        _model_occ_vertices_coincide(g1.center,g2.center) || return false
        abs(g1.r-g2.r)<=tol*max(1.0,g1.r) || return false
        abs(_dot3(g1.n,g2.n))>1.0-1e-12 || return false
        abs((g1.t1-g1.t0)-(g2.t1-g2.t0))<=tol*max(1.0,g1.t1-g1.t0) ||
            return false
        m1=_occ_circle_point(g1,0.5(g1.t0+g1.t1))
        m2=_occ_circle_point(g2,0.5(g2.t0+g2.t1))
        return _model_occ_vertices_coincide(m1,m2)
    elseif g1.occ===:degenerate
        return true                     # single-point locus
    end
    # Identical records parametrize the same locus with the same
    # parameters — they never need a direction flip.
    return _model_occ_geometry_equal(g1,g2)
end

# A stored pcurve record can be reparametrized onto a welded edge's
# parameter map only through `_brep_pc_reparam`, which supports `lin2d`,
# `circ2d`, and `neg2d` wrappers of those. Anything else refuses the weld.
_occ_sew_pc_reparamable(pc)=
    pc===nothing ? true : hasproperty(pc,:neg2d) ?
        _occ_sew_pc_reparamable(pc.neg2d[1]) :
        (hasproperty(pc,:lin2d) || hasproperty(pc,:circ2d))

# Reparametrize a stored pcurve by the member→weld map `s = α·t + β`
# (`α=-1` for direction-flipped members); `_brep_pc_reparam` returns the
# inner record plus the handedness flag the `neg2d` wrapper restores.
function _model_occ_sew_shift_pc(pc,α,β,tol,caller)
    pc===nothing && return nothing
    pc2,neg=_brep_pc_reparam(pc,α,β,tol,caller)
    return neg ? (neg2d=(pc2,)) : pc2
end

# `BRepBuilderAPI_Sewing` + recursive `bind`: emit a fresh, welded
# boundary graph — one vertex per coincidence group, one edge per
# coincident-edge group, a new wire and a new face copy per member —
# bound under fresh tags in upstream's face→wires→edges→vertices
# exploration order, then the shell itself under `allocated`.
function _model_occ_sew_surface_loop!(m::GeoModel,face_ids,allocated,caller)
    kernel=m.occ_kernel; shapes=kernel.shapes; tol=_OCC_CONFUSION
    # Per-member wire use lists, in member order.
    member_wires=[Vector{Tuple{UInt64,Bool}}[
        collect(shapes[wid].geometry.uses) for wid in shapes[fid].children]
        for fid in face_ids]
    # Vertex weld groups over every boundary endpoint.
    vparent=Dict{UInt64,UInt64}()
    function vfind(v)
        r=get(vparent,v,v)
        while haskey(vparent,r); r=vparent[r]; end
        while r!=v; nxt=get(vparent,v,v); vparent[v]=r; v=nxt; end
        return r
    end
    vertices=UInt64[]
    for wires in member_wires, uses in wires, (eid,_) in uses
        a,b=shapes[eid].children
        for v in (a,b)
            vfind(v); v in vertices || push!(vertices,v)
        end
    end
    for i in eachindex(vertices), j in 1:i-1
        a,b=vertices[i],vertices[j]
        vfind(a)==vfind(b) && continue
        _model_occ_vertices_coincide(shapes[a].geometry.point,
            shapes[b].geometry.point) || continue
        vparent[vfind(a)]=vfind(b)
    end
    # Faces that reference each boundary edge (for the flip support gate).
    edge_faces=Dict{UInt64,Vector{Int}}()
    for (fi,wires) in enumerate(member_wires), uses in wires, (eid,_) in uses
        push!(get!(edge_faces,eid,Int[]),fi)
    end
    # Edge weld: greedy groups in first-use (upstream exploration) order.
    # `erep[eid]` is the group's canonical member; `eflip[eid]` records
    # whether the member is stored opposite the canonical direction.
    erep=Dict{UInt64,UInt64}(); eflip=Dict{UInt64,Bool}()
    reps=UInt64[]; emembers=Dict{UInt64,Vector{UInt64}}()
    for wires in member_wires, uses in wires, (eid,_) in uses
        haskey(erep,eid) && continue
        a,b=shapes[eid].children
        ga,gb=vfind(a),vfind(b)
        ga>gb && ((ga,gb)=(gb,ga))
        g=shapes[eid].geometry.geometry
        match=UInt64(0); flip=false
        for r in reps
            ra,rb=shapes[r].children
            rga,rgb=vfind(ra),vfind(rb)
            rga>rgb && ((rga,rgb)=(rgb,rga))
            (rga,rgb)==(ga,gb) || continue
            rg=shapes[r].geometry.geometry
            _model_occ_sew_edge_coincident(g,rg,tol) || continue
            flipped=false
            if ga==gb && g.occ!==:degenerate &&
                    hasproperty(g,:t0) && hasproperty(rg,:t0)
                # Closed edge — both ends weld to one vertex group, so the
                # traversal direction is invisible to endpoints: probe the
                # quarter point against the rep's forward/reversed quarter.
                # Neither coincidence means the curves only share a chord.
                qm=_model_occ_sew_edge_point(g,g.t0+0.25(g.t1-g.t0))
                if qm!==nothing
                    qs=rg.t1-rg.t0
                    qf=_model_occ_sew_edge_point(rg,rg.t0+0.25qs)
                    qr=_model_occ_sew_edge_point(rg,rg.t1-0.25qs)
                    if _model_occ_vertices_coincide(qm,qf)
                        flipped=false
                    elseif qf!==nothing && qr!==nothing &&
                            _model_occ_vertices_coincide(qm,qr)
                        flipped=true
                    else
                        continue
                    end
                end
            else
                flipped=_model_occ_vertices_coincide(shapes[a].geometry.point,
                    shapes[rb].geometry.point) &&
                    !_model_occ_vertices_coincide(shapes[a].geometry.point,
                        shapes[ra].geometry.point)
            end
            if flipped
                g.occ in (:line,:circle,:ellipse) || continue
            end
            # A member joins the weld only when every owner face's stored
            # pcurve entry survives the member→rep parameter map (α=-1 for
            # flipped members, a nonzero shift for off-phase aligned ones —
            # both run through `_brep_pc_reparam`, which supports only
            # `lin2d`/`circ2d`/`neg2d`). Otherwise the member keeps its own
            # copy rather than corrupting a trace it cannot express.
            needs_map=flipped ||
                (hasproperty(g,:t0) && hasproperty(rg,:t0) &&
                 g.t0!=rg.t0)
            if needs_map
                all(edge_faces[eid]) do fi
                    geo=shapes[face_ids[fi]].geometry.geometry
                    hasproperty(geo,:pcurves) || return true
                    ent=get(geo.pcurves,
                        _model_occ_shape_tag(kernel,eid,caller),nothing)
                    ent===nothing && return true
                    _occ_sew_pc_reparamable(get(ent,:fwd,nothing)) &&
                        _occ_sew_pc_reparamable(get(ent,:rev,nothing))
                end || continue
            end
            match=r; flip=flipped; break
        end
        erep[eid]=match==0 ? eid : match
        eflip[eid]=match==0 ? false : flip
        r=erep[eid]
        members=get!(emembers,r,UInt64[])
        eid in members || push!(members,eid)
        match==0 && push!(reps,eid)
    end
    # Upstream recursive-bind order: per member face — face tag, wire
    # tags, then its unique edges (each followed by its own vertices).
    nfaces=length(face_ids)
    nwires=sum(length,member_wires)
    ngroups_v=length(Set(vfind(v) for v in vertices))
    count=nfaces+nwires+length(reps)+ngroups_v+1
    kernel.maxima[5]<=typemax(Int32)-nfaces || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-nwires || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-length(reps) || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[3]<=typemax(Int32)-ngroups_v || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    ids=_model_occ_reserve_shapes(kernel,count,caller)
    cursor=0
    take!()=(cursor+=1; ids[cursor])
    vnew=Dict{UInt64,Tuple{UInt64,Int}}()   # vertex group rep → (id, tag)
    enew=Dict{UInt64,Tuple{UInt64,Int}}()   # edge group rep → (id, tag)
    epos=Dict{UInt64,Int}()                 # edge group rep → pending index
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    # Per-face plan entries: (member index, face tag, new face id,
    # [(wire tag, new wire id) per wire]).
    plans=Tuple{Int,Int,UInt64,Vector{Tuple{Int,UInt64}}}[]
    new_face_ids=UInt64[]
    c0=kernel.maxima[3]; c1=kernel.maxima[4]
    cw=kernel.maxima[2]; cf=kernel.maxima[5]
    for (fi,fid) in enumerate(face_ids)
        cf+=1
        fid_new=take!(); push!(new_face_ids,fid_new)
        wplan=Tuple{Int,UInt64}[]
        for _ in member_wires[fi]
            cw+=1; push!(wplan,(cw,take!()))
        end
        # Unique edge groups in first-use order across this face's wires;
        # each edge binds its welded endpoint vertices inside its own
        # recursion, exactly like upstream `bind(edge)`.
        seen=Set{UInt64}()
        for uses in member_wires[fi], (eid,_) in uses
            r=erep[eid]
            r in seen || push!(seen,r)
            haskey(enew,r) && continue
            c1+=1
            enew[r]=(take!(),c1)
            cids=UInt64[]
            for cv in shapes[r].children
                g=vfind(cv)
                if !haskey(vnew,g)
                    c0+=1
                    vg=shapes[g].geometry
                    vnew[g]=(take!(),c0)
                    push!(pending,(0,c0,vnew[g][1],_ModelOCCShape(0,
                        (point=vg.point,
                         mesh_size=get(vg,:mesh_size,0.0)),UInt64[])))
                end
                push!(cids,vnew[g][1])
            end
            push!(pending,(1,c1,enew[r][1],_ModelOCCShape(1,
                (kind=shapes[r].geometry.kind,
                 geometry=shapes[r].geometry.geometry),cids)))
            epos[r]=length(pending)
        end
        push!(plans,(fi,cf,fid_new,wplan))
    end
    # Wire copies with welded uses, then face copies whose pcurve dicts
    # re-key onto the welded edge tags under the member's parameter map
    # (shift for aligned members, reversal for flipped ones).
    new_geos=Dict{Int,Any}()
    for (fi,ftag,fid_new,wplan) in plans
        new_wire_ids=UInt64[]
        for (wi,(wtag,wid)) in enumerate(wplan)
            uses=member_wires[fi][wi]
            nuses=Tuple{UInt64,Bool}[(enew[erep[e]][1],xor(rev,eflip[e]))
                                     for (e,rev) in uses]
            push!(pending,(-1,wtag,wid,_ModelOCCShape(-1,(uses=nuses,),
                UInt64[first(u) for u in nuses])))
            push!(new_wire_ids,wid)
        end
        fg=shapes[face_ids[fi]].geometry
        geo=fg.geometry
        if hasproperty(geo,:pcurves)
            pcnew=Dict{Int,NamedTuple}()
            for uses in member_wires[fi], (eid,_) in uses
                r=erep[eid]
                key=enew[r][2]
                haskey(pcnew,key) && continue
                old_tag=_model_occ_shape_tag(kernel,eid,caller)
                ent=get(geo.pcurves,old_tag,nothing)
                ent===nothing && continue
                fwd=get(ent,:fwd,nothing); rev=get(ent,:rev,nothing)
                mg=shapes[eid].geometry.geometry
                rg=shapes[r].geometry.geometry
                if eflip[eid]
                    β=mg.t0+rg.t1
                    # The welded edge walks opposite the member: a FORWARD
                    # use of the weld is the member's reversed occurrence,
                    # so the stored fwd/rev trace columns swap (a member
                    # with only `fwd` keeps `rev=nothing` on the weld —
                    # the lookup falls back to the shared column either
                    # way).
                    fsrc=rev===nothing ? fwd : rev
                    ent=(fwd=_model_occ_sew_shift_pc(fsrc,-1.0,β,tol,caller),
                         rev=rev===nothing ? nothing :
                             _model_occ_sew_shift_pc(fwd,-1.0,β,tol,caller))
                elseif hasproperty(mg,:t0) && hasproperty(rg,:t0)
                    δ=mg.t0-rg.t0
                    δ==0.0 || (ent=(fwd=_model_occ_sew_shift_pc(
                                         fwd,1.0,δ,tol,caller),
                                    rev=_model_occ_sew_shift_pc(
                                         rev,1.0,δ,tol,caller)))
                end
                pcnew[key]=ent
            end
            geo=merge(geo,(pcurves=pcnew,))
        end
        new_geos[fi]=geo
        push!(pending,(2,ftag,fid_new,
            _ModelOCCShape(2,(kind=fg.kind,geometry=geo),new_wire_ids)))
    end
    # Edges whose records embed a face (degenerate apex/pole edges, surface
    # isos) rebind `face` to a sewn owner geometry so their pcurve lookups
    # resolve under the welded tags — owners with a pcurve entry for the
    # member are preferred so the degenerate trace stays available.
    for (r,idx) in epos
        d,t,i,s=pending[idx]
        g=s.geometry.geometry
        hasproperty(g,:face) || continue
        owner=typemax(Int)
        for e in emembers[r], fi in edge_faces[e]
            fg=shapes[face_ids[fi]].geometry.geometry
            hasproperty(fg,:pcurves) || continue
            haskey(fg.pcurves,_model_occ_shape_tag(kernel,e,caller)) ||
                continue
            fi<owner && (owner=fi)
        end
        owner==typemax(Int) && (owner=minimum(
            fi for e in emembers[r] for fi in edge_faces[e]))
        pending[idx]=(d,t,i,_ModelOCCShape(1,
            (kind=s.geometry.kind,geometry=merge(g,(face=new_geos[owner],))),
            s.children))
    end
    sid=take!()
    suse=Tuple{UInt64,Bool}[(f,false) for f in new_face_ids]
    push!(pending,(-2,allocated,sid,
        _ModelOCCShape(-2,(uses=suse,),new_face_ids)))
    cursor==count || throw(ArgumentError(
        "$caller: internal sewn-graph accounting mismatch"))
    for (d,t,i,s) in pending
        _model_occ_bind_shape!(kernel,d,t,i,s)
    end
    kernel.next_shape=ids[end]+UInt64(1)
    return allocated
end

function _model_occ_add_volume!(m::GeoModel,shell_tags;tag=-1)
    caller="OCC.add_volume"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    shell_ids=UInt64[]
    for raw in shell_tags
        current=_model_occ_requested_tag(raw,caller)
        id=get(kernel.bindings,(-2,current),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: unknown OCC Surface Loop[$current]"))
        push!(shell_ids,id)
    end
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,3,allocated,id,
        _ModelOCCShape(3,(kind=:solid,),shell_ids))
    kernel.next_shape=id+UInt64(1)
    return allocated
end

function _model_occ_add_box!(m::GeoModel,x,y,z,dx,dy,dz;tag=-1)
    caller="OCC.add_box"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    corner=(_finite_scalar(x,caller,"x"),_finite_scalar(y,caller,"y"),
            _finite_scalar(z,caller,"z"))
    ext=(_finite_scalar(dx,caller,"dx"),_finite_scalar(dy,caller,"dy"),
         _finite_scalar(dz,caller,"dz"))
    any(==(0.0),ext) && throw(ArgumentError("$caller: Degenerate box"))
    # BRepPrimAPI_MakeBox(P1,P2): the box between opposite corners —
    # negative extents simply place P2 below P1 on that axis.
    lo=ntuple(i->min(corner[i],corner[i]+ext[i]),3)
    hi=ntuple(i->max(corner[i],corner[i]+ext[i]),3)
    v=((lo[1],lo[2],lo[3]),(hi[1],lo[2],lo[3]),(hi[1],hi[2],lo[3]),
       (lo[1],hi[2],lo[3]),(lo[1],lo[2],hi[3]),(hi[1],lo[2],hi[3]),
       (hi[1],hi[2],hi[3]),(lo[1],hi[2],hi[3]))
    # Edge list: 4 bottom, 4 top, 4 vertical.
    ep=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
        (1,5),(2,6),(3,7),(4,8))
    # Per-face wire traversal, wound counterclockwise about the outward
    # normal — bottom/top/front/back/left/right; -k walks edge k against
    # its stored direction.
    face_loops=(((-4,-3,-2,-1),(0.0,0.0,-1.0)),
                ((5,6,7,8),(0.0,0.0,1.0)),
                ((1,10,-5,-9),(0.0,-1.0,0.0)),
                ((12,-7,-11,3),(0.0,1.0,0.0)),
                ((9,-8,-12,4),(-1.0,0.0,0.0)),
                ((2,11,-6,-10),(1.0,0.0,0.0)))
    count=8+12+6+6+1+1
    kernel.maxima[3]<=typemax(Int32)-8 || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-12 || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-6 || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-6 || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    ids=_model_occ_reserve_shapes(kernel,count,caller)
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    vertex_ids=ids[1:8]
    for i in 1:8
        push!(pending,(0,kernel.maxima[3]+i,vertex_ids[i],
            _ModelOCCShape(0,(point=v[i],mesh_size=0.0),UInt64[])))
    end
    edge_ids=ids[9:20]
    for i in 1:12
        a,b=ep[i]
        pa,pb=v[a],v[b]
        diff=_arc_sub(pb,pa)
        geometry=(occ=:line,origin=pa,dir=_occ_dir(diff),
                  t0=0.0,t1=_occ_modulus(diff))
        push!(pending,(1,kernel.maxima[4]+i,edge_ids[i],
            _ModelOCCShape(1,(kind=:line,geometry=geometry),
                           UInt64[vertex_ids[a],vertex_ids[b]])))
    end
    wire_ids=ids[21:26]
    face_ids=ids[27:32]
    shape_index=Dict{UInt64,_ModelOCCShape}()
    for (_,_,id,shape) in pending
        shape_index[id]=shape
    end
    edge_tag_of=Dict{UInt64,Int}(edge_ids[i]=>kernel.maxima[4]+i for i in 1:12)
    for (index,(loop,_)) in enumerate(face_loops)
        uses=Tuple{UInt64,Bool}[(edge_ids[abs(k)],k<0) for k in loop]
        push!(pending,(-1,kernel.maxima[2]+index,wire_ids[index],
            _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses])))
    end
    shape_of=id->get(shape_index,id) do
        kernel.shapes[id]
    end
    tag_of=id->get(()->_model_occ_shape_tag(kernel,id,caller),
                   edge_tag_of,id)
    for (index,(loop,normal)) in enumerate(face_loops)
        uses=Tuple{UInt64,Bool}[(edge_ids[abs(k)],k<0) for k in loop]
        geometry=_model_occ_plane_face_geometry(
            shape_of,tag_of,[uses],caller)
        # The wire is wound outward: the Newell fit must point along the
        # intended normal — a mismatch means a winding bug above.
        _dot3(geometry.axis,normal)>0.999 || throw(ArgumentError(
            "$caller: internal box face winding is inconsistent"))
        push!(pending,(2,kernel.maxima[5]+index,face_ids[index],
            _ModelOCCShape(2,(kind=:plane,geometry=geometry),
                           UInt64[wire_ids[index]])))
    end
    shell_id=ids[33]
    uses=Tuple{UInt64,Bool}[(fid,false) for fid in face_ids]
    push!(pending,(-2,kernel.maxima[1]+1,shell_id,
        _ModelOCCShape(-2,(uses=uses,),face_ids)))
    solid_id=ids[34]
    push!(pending,(3,allocated,solid_id,
        _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id])))
    for (dimension,tag,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,tag,id,shape)
    end
    kernel.next_shape=solid_id+UInt64(1)
    return allocated
end

function _model_occ_add_wedge!(m::GeoModel,x,y,z,dx,dy,dz;
                               tag=-1,ltx=0.0,zAxis=Float64[])
    caller="OCC.add_wedge"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    corner=_finite3(x,y,z,caller)
    ex=_finite_scalar(dx,caller,"dx")
    ey=_finite_scalar(dy,caller,"dy")
    ez=_finite_scalar(dz,caller,"dz")
    lt=_finite_scalar(ltx,caller,"ltx")
    (ex==0.0 || ey==0.0 || ez==0.0) && throw(ArgumentError(
        "$caller: Degenerate wedge"))
    axis_dir=if isempty(zAxis)
        (0.0,0.0,1.0)
    else
        length(zAxis)==3 || throw(ArgumentError(
            "$caller: zAxis must be a 3-vector"))
        dir=_finite3(zAxis[1],zAxis[2],zAxis[3],caller)
        _occ_modulus(dir)>0 || throw(ArgumentError(
            "$caller: zAxis must have positive length"))
        _occ_dir(dir)
    end
    Xf,Yf=_occ_ax2_frame(axis_dir)
    world(lx,ly,lz)=_occ3_add(corner,_occ3_add(
        _occ3_scale(Xf,lx),_occ3_add(_occ3_scale(Yf,ly),
                                     _occ3_scale(axis_dir,lz))))
    # BRepPrim_Wedge profile in the local xz plane, extruded along +y.
    # ltx==0 merges the top vertices into the ridge — a triangle.
    profile=lt==0.0 ? ((0.0,0.0),(ex,0.0),(0.0,ez)) :
                      ((0.0,0.0),(ex,0.0),(lt,ez),(0.0,ez))
    k=length(profile)
    shoelace=sum(profile[i][1]*profile[mod1(i+1,k)][2]-
                 profile[mod1(i+1,k)][1]*profile[i][2] for i in 1:k)
    # Traversal sign: wound as (a0→b0→b1→a1) a side quad's Newell is
    # left-normal·dy — outward exactly when sign(dy·shoelace)<0.
    side_forward=sign(ey*shoelace)<0
    cap0_forward=sign(shoelace)==sign(ey)
    nverts,nedges,nfaces=2k,3k,k+2
    kernel.maxima[3]<=typemax(Int32)-nverts || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-nedges || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-nfaces || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-nfaces || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    ids=_model_occ_reserve_shapes(kernel,nverts+nedges+2nfaces+2,caller)
    vertex_ids=ids[1:nverts]
    edge_ids=ids[nverts+1:nverts+nedges]
    wire_ids=ids[nverts+nedges+1:nverts+nedges+nfaces]
    face_ids=ids[nverts+nedges+nfaces+1:nverts+nedges+2nfaces]
    shell_id=ids[end-1]
    solid_id=ids[end]
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    for i in 1:k
        for (ring,ly) in ((0,0.0),(k,ey))
            px,pz=profile[i]
            push!(pending,(0,kernel.maxima[3]+ring+i,vertex_ids[ring+i],
                _ModelOCCShape(0,(point=world(px,ly,pz),mesh_size=0.0),
                               UInt64[])))
        end
    end
    for j in 1:k
        a,b=j,mod1(j+1,k)
        pa,pb=profile[a],profile[b]
        for (ring,ly) in ((0,0.0),(k,ey))
            w0=world(pa[1],ly,pa[2]); w1=world(pb[1],ly,pb[2])
            diff=_arc_sub(w1,w0)
            push!(pending,(1,kernel.maxima[4]+ring+j,edge_ids[ring+j],
                _ModelOCCShape(1,(kind=:line,geometry=
                    (occ=:line,origin=w0,dir=_occ_dir(diff),
                     t0=0.0,t1=_occ_modulus(diff))),
                    UInt64[vertex_ids[ring+a],vertex_ids[ring+b]])))
        end
        w0=world(pa[1],0.0,pa[2]); w1=world(pa[1],ey,pa[2])
        diff=_arc_sub(w1,w0)
        push!(pending,(1,kernel.maxima[4]+2k+j,edge_ids[2k+j],
            _ModelOCCShape(1,(kind=:line,geometry=
                (occ=:line,origin=w0,dir=_occ_dir(diff),
                 t0=0.0,t1=_occ_modulus(diff))),
                UInt64[vertex_ids[a],vertex_ids[k+a]])))
    end
    # Wire order lists, signed edge ids: k profile edges per ring, then k
    # extrusion edges; -k walks edge k against its stored direction.
    wires=Vector{Int}[]
    for j in 1:k
        a,b=j,mod1(j+1,k)
        push!(wires,side_forward ?
            [j,2k+b,-(k+j),-(2k+a)] : [2k+a,k+j,-(2k+b),-j])
    end
    order=collect(1:k)
    cap0=cap0_forward ? order : [-e for e in Iterators.reverse(order)]
    push!(wires,cap0)
    ring=order.+k
    cap1=!cap0_forward ? ring : [-e for e in Iterators.reverse(ring)]
    push!(wires,cap1)
    shape_index=Dict{UInt64,_ModelOCCShape}()
    for (_,_,id,shape) in pending
        shape_index[id]=shape
    end
    edge_tag_of=Dict{UInt64,Int}(edge_ids[i]=>kernel.maxima[4]+i
                                 for i in 1:nedges)
    for (index,loop) in enumerate(wires)
        uses=Tuple{UInt64,Bool}[(edge_ids[abs(e)],e<0) for e in loop]
        push!(pending,(-1,kernel.maxima[2]+index,wire_ids[index],
            _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses])))
    end
    shape_of=id->get(shape_index,id) do
        kernel.shapes[id]
    end
    tag_of=id->get(()->_model_occ_shape_tag(kernel,id,caller),
                   edge_tag_of,id)
    for (index,loop) in enumerate(wires)
        uses=Tuple{UInt64,Bool}[(edge_ids[abs(e)],e<0) for e in loop]
        geometry=_model_occ_plane_face_geometry(
            shape_of,tag_of,[uses],caller)
        # Outward check like add_box: side faces read the profile-edge
        # left normal with the traversal sign; caps follow ∓sign(dy)·Y.
        if index<=k
            a,b=index,mod1(index+1,k)
            p=profile[b][1]-profile[a][1]; q=profile[b][2]-profile[a][2]
            s=sign(shoelace)
            s==0.0 && continue
            nlocal=(q*s,0.0,-p*s)
            m=_occ_modulus(nlocal)
            m==0.0 && continue
            intended=_occ3_add(_occ3_scale(Xf,nlocal[1]/m),
                _occ3_scale(axis_dir,nlocal[3]/m))
        else
            intended=_occ3_scale(Yf,index==k+1 ? -sign(ey) : sign(ey))
        end
        _dot3(geometry.axis,intended)>0.999 || throw(ArgumentError(
            "$caller: internal wedge face winding is inconsistent"))
        push!(pending,(2,kernel.maxima[5]+index,face_ids[index],
            _ModelOCCShape(2,(kind=:plane,geometry=geometry),
                           UInt64[wire_ids[index]])))
    end
    uses=Tuple{UInt64,Bool}[(fid,false) for fid in face_ids]
    push!(pending,(-2,kernel.maxima[1]+1,shell_id,
        _ModelOCCShape(-2,(uses=uses,),face_ids)))
    push!(pending,(3,allocated,solid_id,
        _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id])))
    for (dimension,tag,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,tag,id,shape)
    end
    kernel.next_shape=solid_id+UInt64(1)
    return allocated
end
