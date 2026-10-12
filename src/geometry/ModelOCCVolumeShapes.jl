# `gmsh.model.occ.addSurfaceLoop(surfaceTags, tag=-1, sewing=false)` — the
# `BRepPrim_Builder::MakeShell` + `AddShellFace` path: a shell over the
# bound faces, all FORWARD, with no closure check — upstream defers
# manifold validity to ShapeFix/meshing, and the emitted `m.surface_loops`
# entry follows the same deferred rule. `sewing=true`
# (`BRepBuilderAPI_Sewing`, welding geometrically coincident but
# topologically distinct edges into shell-local face copies) is not
# implemented yet and rejects loudly rather than approximating.
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
    sewing && throw(ArgumentError(
        "$caller: sewing is not supported on OCC surface loops"))
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
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    uses=Tuple{UInt64,Bool}[(fid,false) for fid in face_ids]
    _model_occ_bind_shape!(kernel,-2,allocated,id,
        _ModelOCCShape(-2,(uses=uses,),face_ids))
    kernel.next_shape=id+UInt64(1)
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
