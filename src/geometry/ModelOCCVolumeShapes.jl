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
