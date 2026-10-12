# `BRepPrimAPI_MakeCylinder`/`MakeCone`/`MakeSphere`/`MakeTorus` — the
# primitive OneAxis/revolution solids bound into the OCC shape registry. The
# emitted graphs replicate the layouts `BRepPrim_OneAxis` produces upstream
# (the same layouts `_materialize_*!` materialize for the geo pipeline):
# cylinder/cone carry two rim vertices, two rim circles and one seam line,
# lateral + cap faces and a `[lateral,+top,-bottom]` shell; a sphere carries
# two pole vertices with degenerate edges and a trimmed meridian; a torus
# carries rim vertices, an equator and meridian circles, and Plane caps for
# partial sweeps. All sub-entities — vertices, edges, wires, faces, shell,
# solid — publish atomically like `add_box`.
#
# Partial sweeps are genuine backlog: `BRepPrim_OneAxis` produces extra
# rectangular side faces whose layout has no verified reference yet, so
# `angle`/`angle1`/`angle2`/`angle3` away from the full-rotation defaults
# reject loudly rather than approximating an unverified topology.

# UV bounds for a `:plane` record built against not-yet-bound edges:
# `trange` maps a pcurve's edge tag to its (t0,t1) span.
function _model_occ_prim_uvb(pcurves,trange)
    umin=vmin=Inf; umax=vmax=-Inf
    for (ct,e) in pcurves
        t0,t1=trange[ct]
        for pc in (e.fwd,e.rev)
            pc===nothing && continue
            lo,hi=_occ_pcurve_fillbndbox(pc,t0,t1)
            umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
            umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
        end
    end
    return ([umin,vmin],[umax,vmax])
end

function _model_occ_add_cylinder!(m::GeoModel,x,y,z,dx,dy,dz,r;
                                  tag=-1,angle=2π)
    caller="OCC.add_cylinder"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    base=_finite3(x,y,z,caller)
    axis=_finite3(dx,dy,dz,caller)
    h=_occ_modulus(axis)
    h>0 || throw(ArgumentError(
        "$caller: Cannot build cylinder of zero height"))
    radius=_finite_scalar(r,caller,"radius")
    radius>0 || throw(ArgumentError("$caller: radius must be positive"))
    a=_finite_scalar(angle,caller,"angle")
    (a>0.0 && a<=_OCC_TWO_PI) || throw(ArgumentError(
        "$caller: Cannot build cylinder with angle <= 0 or angle > 2*Pi"))
    a==_OCC_TWO_PI || throw(ArgumentError(
        "$caller: partial cylinder angles are not implemented"))
    # Same allocation discipline as add_box: every sub-tag is checked
    # against Int32 before any shape is bound.
    kernel.maxima[3]<=typemax(Int32)-2 || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    n=_occ_dir(axis)
    X,Y=_occ_ax2_frame(n)
    M,loc=_occ_end_trsf(base,n)
    Xc,Yc=_occ_ax2_setx(n,X)
    nc=_occ_crossed(Xc,Yc)
    w_top=_occ_vertex(base,n,X,(radius,h),M,loc)
    w_bot=_occ_vertex(base,n,X,(radius,0.0),M,loc)
    top_center=_occ_translated(base,n,h)
    seam_origin=_occ_xform_pnt(M,loc,_occ_translated(base,X,radius))
    seam_dir=_occ_xform_dir(M,n)
    ids=_model_occ_reserve_shapes(kernel,13,caller)
    vid_top,vid_bot=ids[1:2]
    cid_top,cid_seam,cid_bot=ids[3:5]
    wid_lat,wid_top,wid_bot=ids[6:8]
    fid_lat,fid_top,fid_bot=ids[9:11]
    shell_id,solid_id=ids[12:13]
    t_top=kernel.maxima[4]+1; t_seam=kernel.maxima[4]+2
    t_bot=kernel.maxima[4]+3
    s_lat=kernel.maxima[5]+1; s_top=kernel.maxima[5]+2
    s_bot=kernel.maxima[5]+3
    top_pc=Dict{Int,NamedTuple}(
        t_top=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),radius),),rev=nothing))
    bot_pc=Dict{Int,NamedTuple}(
        t_bot=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),radius),),rev=nothing))
    trange=Dict(t_top=>(0.0,_OCC_TWO_PI),t_bot=>(0.0,_OCC_TWO_PI))
    pending=[
        (0,kernel.maxima[3]+1,vid_top,
            _ModelOCCShape(0,(point=w_top,mesh_size=0.0),UInt64[])),
        (0,kernel.maxima[3]+2,vid_bot,
            _ModelOCCShape(0,(point=w_bot,mesh_size=0.0),UInt64[])),
        (1,t_top,cid_top,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=top_center,n=nc,X=Xc,Y=Yc,
                          r=radius,t0=0.0,t1=_OCC_TWO_PI)),
                           UInt64[vid_top,vid_top])),
        (1,t_seam,cid_seam,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=seam_origin,dir=seam_dir,
                          t0=0.0,t1=h)),
                           UInt64[vid_bot,vid_top])),
        (1,t_bot,cid_bot,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=base,n=nc,X=Xc,Y=Yc,
                          r=radius,t0=0.0,t1=_OCC_TWO_PI)),
                           UInt64[vid_bot,vid_bot])),
        (-1,kernel.maxima[2]+1,wid_lat,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_top,true),
                (cid_seam,true),(cid_bot,false),(cid_seam,false)],),
                           UInt64[cid_top,cid_seam,cid_bot,cid_seam])),
        (-1,kernel.maxima[2]+2,wid_top,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_top,false)],),
                           UInt64[cid_top])),
        (-1,kernel.maxima[2]+3,wid_bot,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_bot,false)],),
                           UInt64[cid_bot])),
        (2,s_lat,fid_lat,
            _ModelOCCShape(2,(kind=:cylinder,
                geometry=(occ=:cylinder,center=base,axis=n,X=X,Y=Y,
                          radius=radius,height=h,
                          pcurves=Dict{Int,NamedTuple}(
                              t_top=>(fwd=(lin2d=((0.0,h),(1.0,0.0)),),
                                      rev=nothing),
                              t_bot=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),
                                      rev=nothing),
                              t_seam=>(fwd=(lin2d=((_OCC_TWO_PI,0.0),
                                                   (0.0,1.0)),),
                                       rev=(lin2d=((0.0,0.0),(0.0,1.0)),))))),
                           UInt64[wid_lat])),
        (2,s_top,fid_top,
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=top_center,axis=n,X=Xc,Y=Yc,
                          reversed=false,pcurves=top_pc,
                          uvb=_model_occ_prim_uvb(top_pc,trange))),
                           UInt64[wid_top])),
        (2,s_bot,fid_bot,
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=base,axis=n,X=Xc,Y=Yc,
                          reversed=true,pcurves=bot_pc,
                          uvb=_model_occ_prim_uvb(bot_pc,trange))),
                           UInt64[wid_bot])),
        (-2,kernel.maxima[1]+1,shell_id,
            _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[(fid_lat,false),
                (fid_top,false),(fid_bot,true)],),
                           UInt64[fid_lat,fid_top,fid_bot])),
        (3,allocated,solid_id,
            _ModelOCCShape(3,(kind=:solid,
                primitive=(occ=:cylinder,
                           record=(center=base,axis=axis,radius=radius,
                                   height=h))),UInt64[shell_id]))]
    for (dimension,ptag,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,ptag,id,shape)
    end
    kernel.next_shape=solid_id+UInt64(1)
    return allocated
end

function _model_occ_add_cone!(m::GeoModel,x,y,z,dx,dy,dz,r1,r2;
                              tag=-1,angle=2π)
    caller="OCC.add_cone"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    base=_finite3(x,y,z,caller)
    axis=_finite3(dx,dy,dz,caller)
    h=_occ_modulus(axis)
    h>0 || throw(ArgumentError(
        "$caller: Cannot build cone of zero height"))
    ra=_finite_scalar(r1,caller,"r1")
    rb=_finite_scalar(r2,caller,"r2")
    (ra>=0.0 && rb>=0.0 && (ra>0.0 || rb>0.0)) || throw(ArgumentError(
        "$caller: radii must be non-negative with at least one positive"))
    a=_finite_scalar(angle,caller,"angle")
    a>0.0 || throw(ArgumentError(
        "$caller: Cone angle should be positive"))
    a==_OCC_TWO_PI || throw(ArgumentError(
        "$caller: partial cone angles are not implemented"))
    apex_top=rb==0.0; apex_bot=ra==0.0
    kernel.maxima[3]<=typemax(Int32)-2 || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    n=_occ_dir(axis)
    X,Y=_occ_ax2_frame(n)
    M,loc=_occ_end_trsf(base,n)
    dr=rb-ra
    sa=_gm_atan(dr/h)
    len=sqrt(fma(h,h,dr*dr))
    dir2d=_occ_dir2(_gm_sincos(sa))
    mp(v)=_occ_lin2d_eval((ra,0.0),dir2d,v)
    mpT=mp(len); mpB=mp(0.0)
    Xc,Yc=_occ_ax2_setx(n,X)
    nc=_occ_crossed(Xc,Yc)
    MY,_=_occ_trsf_rot(base,Y,sa)
    d0=_occ_rot_dir(MY,n)
    w_top=_occ_vertex(base,n,X,mpT,M,loc)
    w_bot=_occ_vertex(base,n,X,mpB,M,loc)
    top_center=_occ_translated(base,n,mpT[2])
    seam_origin=_occ_xform_pnt(M,loc,_occ_translated(base,X,ra))
    seam_dir=_occ_xform_dir(M,d0)
    ids=_model_occ_reserve_shapes(kernel,13,caller)
    vid_top,vid_bot=ids[1:2]
    cid_top,cid_seam,cid_bot=ids[3:5]
    wid_lat,wid_top,wid_bot=ids[6:8]
    fid_lat,fid_top,fid_bot=ids[9:11]
    shell_id,solid_id=ids[12:13]
    t_top=kernel.maxima[4]+1; t_seam=kernel.maxima[4]+2
    t_bot=kernel.maxima[4]+3
    s_lat=kernel.maxima[5]+1; s_top=kernel.maxima[5]+2
    s_bot=kernel.maxima[5]+3
    lat_pcurves=Dict{Int,NamedTuple}(
        t_top=>(fwd=(lin2d=((0.0,len),(1.0,0.0)),),rev=nothing),
        t_bot=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing),
        t_seam=>(fwd=(lin2d=((_OCC_TWO_PI,-0.0),(0.0,1.0)),),
                 rev=(lin2d=((0.0,-0.0),(0.0,1.0)),)))
    lat_geometry=(occ=:cone,center=base,axis=n,X=Xc,Y=Yc,r1=ra,r2=rb,
                  height=h,pcurves=lat_pcurves)
    # A zero-radius end collapses to a degenerate apex edge bound to the
    # lateral face (`setTrimmed`); a positive end keeps its rim circle and
    # Plane cap. Tag bookkeeping still reserves both ends' tags so an apex
    # case is just the same record family minus the cap.
    trange=Dict{Int,NTuple{2,Float64}}()
    apex_top || (trange[t_top]=(0.0,_OCC_TWO_PI))
    apex_bot || (trange[t_bot]=(0.0,_OCC_TWO_PI))
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    push!(pending,(0,kernel.maxima[3]+1,vid_top,
        _ModelOCCShape(0,(point=w_top,mesh_size=0.0),UInt64[])))
    push!(pending,(0,kernel.maxima[3]+2,vid_bot,
        _ModelOCCShape(0,(point=w_bot,mesh_size=0.0),UInt64[])))
    if apex_top
        push!(pending,(1,t_top,cid_top,
            _ModelOCCShape(1,(kind=:degenerate,
                geometry=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI,
                          face=lat_geometry)),
                           UInt64[vid_top,vid_top])))
    else
        push!(pending,(1,t_top,cid_top,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=top_center,n=nc,X=Xc,Y=Yc,
                          r=mpT[1],t0=0.0,t1=_OCC_TWO_PI)),
                           UInt64[vid_top,vid_top])))
    end
    push!(pending,(1,t_seam,cid_seam,
        _ModelOCCShape(1,(kind=:line,
            geometry=(occ=:line,origin=seam_origin,dir=seam_dir,
                      t0=0.0,t1=len)),
                       UInt64[vid_bot,vid_top])))
    if apex_bot
        push!(pending,(1,t_bot,cid_bot,
            _ModelOCCShape(1,(kind=:degenerate,
                geometry=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI,
                          face=lat_geometry)),
                           UInt64[vid_bot,vid_bot])))
    else
        push!(pending,(1,t_bot,cid_bot,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=base,n=nc,X=Xc,Y=Yc,
                          r=mpB[1],t0=0.0,t1=_OCC_TWO_PI)),
                           UInt64[vid_bot,vid_bot])))
    end
    push!(pending,(-1,kernel.maxima[2]+1,wid_lat,
        _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_top,true),
            (cid_seam,true),(cid_bot,false),(cid_seam,false)],),
                       UInt64[cid_top,cid_seam,cid_bot,cid_seam])))
    apex_top || push!(pending,(-1,kernel.maxima[2]+2,wid_top,
        _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_top,false)],),
                       UInt64[cid_top])))
    apex_bot || push!(pending,(-1,kernel.maxima[2]+3,wid_bot,
        _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_bot,false)],),
                       UInt64[cid_bot])))
    push!(pending,(2,s_lat,fid_lat,
        _ModelOCCShape(2,(kind=:cone,geometry=lat_geometry),
                       UInt64[wid_lat])))
    shell_uses=Tuple{UInt64,Bool}[(fid_lat,false)]
    if !apex_top
        top_pc=Dict{Int,NamedTuple}(
            t_top=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),mpT[1]),),
                    rev=nothing))
        push!(pending,(2,s_top,fid_top,
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=top_center,axis=n,X=Xc,Y=Yc,
                          reversed=false,pcurves=top_pc,
                          uvb=_model_occ_prim_uvb(top_pc,trange))),
                           UInt64[wid_top])))
        push!(shell_uses,(fid_top,false))
    end
    if !apex_bot
        bot_pc=Dict{Int,NamedTuple}(
            t_bot=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),mpB[1]),),
                    rev=nothing))
        push!(pending,(2,s_bot,fid_bot,
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=base,axis=n,X=Xc,Y=Yc,
                          reversed=true,pcurves=bot_pc,
                          uvb=_model_occ_prim_uvb(bot_pc,trange))),
                           UInt64[wid_bot])))
        push!(shell_uses,(fid_bot,true))
    end
    push!(pending,(-2,kernel.maxima[1]+1,shell_id,
        _ModelOCCShape(-2,(uses=shell_uses,),
                       UInt64[u[1] for u in shell_uses])))
    push!(pending,(3,allocated,solid_id,
        _ModelOCCShape(3,(kind=:solid,
            primitive=(occ=:cone,
                       record=(center=base,axis=axis,r1=ra,r2=rb,
                               height=h))),UInt64[shell_id])))
    for (dimension,ptag,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,ptag,id,shape)
    end
    kernel.next_shape=solid_id+UInt64(1)
    return allocated
end

function _model_occ_add_sphere!(m::GeoModel,xc,yc,zc,radius;
                                tag=-1,angle1=-π/2,angle2=π/2,angle3=2π)
    caller="OCC.add_sphere"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    center=_finite3(xc,yc,zc,caller)
    r=_finite_scalar(radius,caller,"radius")
    r>0 || throw(ArgumentError(
        "$caller: Sphere radius should be positive"))
    a1=_finite_scalar(angle1,caller,"angle1")
    a2=_finite_scalar(angle2,caller,"angle2")
    a3=_finite_scalar(angle3,caller,"angle3")
    (a3>0.0 && a3<=_OCC_TWO_PI) || throw(ArgumentError(
        "$caller: Cannot build sphere with angle <= 0 or angle > 2*Pi"))
    (a1==-π/2 && a2==π/2 && a3==_OCC_TWO_PI) || throw(ArgumentError(
        "$caller: partial sphere angles are not implemented"))
    kernel.maxima[3]<=typemax(Int32)-2 || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-3 || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    n=(0.0,0.0,1.0)
    X,Y=_occ_ax2_setx(n,(1.0,0.0,0.0))
    M,loc=_occ_end_trsf(center,n)
    mp(v)=_occ_circ2d_eval((0.0,0.0),(1.0,0.0),r,v)
    w_n=_occ_vertex(center,n,X,mp(π/2),M,loc)
    w_s=_occ_vertex(center,n,X,mp(-π/2),M,loc)
    Dm=(-Y[1],-Y[2],-Y[3])
    Xm,Ym=_occ_ax2_setx(Dm,X)
    mer=_occ_xform_ax2(M,loc,(loc=center,axis=Dm,X=Xm,Y=Ym))
    ids=_model_occ_reserve_shapes(kernel,9,caller)
    vid_n,vid_s=ids[1:2]
    cid_n,cid_mer,cid_s=ids[3:5]
    wid=ids[6]
    fid,shell_id,solid_id=ids[7:9]
    t_n=kernel.maxima[4]+1; t_mer=kernel.maxima[4]+2
    t_s=kernel.maxima[4]+3
    s_face=kernel.maxima[5]+1
    face_geometry=(occ=:sphere,center=center,radius=r,
                   axis=n,X=X,Y=Y,
                   pcurves=Dict{Int,NamedTuple}(
                       t_n=>(fwd=(lin2d=((0.0,π/2),(1.0,0.0)),),
                             rev=nothing),
                       t_s=>(fwd=(lin2d=((0.0,-π/2),(1.0,0.0)),),
                             rev=nothing),
                       t_mer=>(fwd=(lin2d=((_OCC_TWO_PI,-_OCC_TWO_PI),
                                           (0.0,1.0)),),
                               rev=(lin2d=((0.0,-_OCC_TWO_PI),(0.0,1.0)),))))
    pending=[
        (0,kernel.maxima[3]+1,vid_n,
            _ModelOCCShape(0,(point=w_n,mesh_size=0.0),UInt64[])),
        (0,kernel.maxima[3]+2,vid_s,
            _ModelOCCShape(0,(point=w_s,mesh_size=0.0),UInt64[])),
        (1,t_n,cid_n,
            _ModelOCCShape(1,(kind=:degenerate,
                geometry=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI,
                          face=face_geometry)),
                           UInt64[vid_n,vid_n])),
        (1,t_mer,cid_mer,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=mer.loc,n=mer.axis,X=mer.X,
                          Y=mer.Y,r=r,t0=1.5π,t1=2.5π)),
                           UInt64[vid_s,vid_n])),
        (1,t_s,cid_s,
            _ModelOCCShape(1,(kind=:degenerate,
                geometry=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI,
                          face=face_geometry)),
                           UInt64[vid_s,vid_s])),
        (-1,kernel.maxima[2]+1,wid,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_n,true),
                (cid_mer,true),(cid_s,false),(cid_mer,false)],),
                           UInt64[cid_n,cid_mer,cid_s,cid_mer])),
        (2,s_face,fid,
            _ModelOCCShape(2,(kind=:sphere,geometry=face_geometry),
                           UInt64[wid])),
        (-2,kernel.maxima[1]+1,shell_id,
            _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[(fid,false)],),
                           UInt64[fid])),
        (3,allocated,solid_id,
            _ModelOCCShape(3,(kind=:solid,
                primitive=(occ=:sphere,
                           record=(center=center,radius=r))),
                UInt64[shell_id]))]
    for (dimension,ptag,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,ptag,id,shape)
    end
    kernel.next_shape=solid_id+UInt64(1)
    return allocated
end

function _model_occ_add_torus!(m::GeoModel,x,y,z,r1,r2;
                               tag=-1,angle=2π,zAxis=Float64[])
    caller="OCC.add_torus"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,3,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
        "$caller: OCC Volume[$allocated] already exists"))
    center=_finite3(x,y,z,caller)
    ra=_finite_scalar(r1,caller,"r1")
    rb=_finite_scalar(r2,caller,"r2")
    (ra>0.0 && rb>0.0) || throw(ArgumentError(
        "$caller: Torus radii should be positive"))
    a=_finite_scalar(angle,caller,"angle")
    (a>0.0 && a<=_OCC_TWO_PI) || throw(ArgumentError(
        "$caller: angle must lie in (0, 2*Pi]"))
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
    full=a>=_OCC_TWO_PI
    # Full torus: rim vertex, equator + meridian circles, one Torus face,
    # shell [face]. Partial: two rim vertices, a trimmed equator arc, two
    # meridians, the Torus face, and two Plane caps.
    nvertices=full ? 1 : 2
    ncurves=full ? 2 : 3
    nwires=full ? 1 : 3
    nfaces=full ? 1 : 3
    kernel.maxima[3]<=typemax(Int32)-nvertices || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-ncurves || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-nwires || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-nfaces || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    n=axis_dir
    X,Y=_occ_ax2_frame(n)
    M,loc=_occ_end_trsf(center,n)
    mp(v)=_occ_circ2d_eval((ra,0.0),(1.0,0.0),rb,v)
    mpE=mp(_OCC_TWO_PI)
    w_rim=_occ_vertex(center,n,X,mpE,M,loc)
    Dm=(-Y[1],-Y[2],-Y[3])
    Xm,Ym=_occ_ax2_setx(Dm,X)
    locm=_occ_translated(center,X,ra)
    Xe,Ye=_occ_ax2_setx(n,X)
    ne=_occ_crossed(Xe,Ye)
    ce=_occ_translated(center,n,mpE[2])
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    if full
        ids=_model_occ_reserve_shapes(kernel,7,caller)
        vid_rim,cid_eq,cid_mer,wid,fid,shell_id,solid_id=ids
        t_eq=kernel.maxima[4]+1; t_mer=kernel.maxima[4]+2
        s_face=kernel.maxima[5]+1
        mer=_occ_xform_ax2(M,loc,(loc=locm,axis=Dm,X=Xm,Y=Ym))
        pending=[
            (0,kernel.maxima[3]+1,vid_rim,
                _ModelOCCShape(0,(point=w_rim,mesh_size=0.0),UInt64[])),
            (1,t_eq,cid_eq,
                _ModelOCCShape(1,(kind=:circle,
                    geometry=(occ=:circle,center=ce,n=ne,X=Xe,Y=Ye,
                              r=mpE[1],t0=0.0,t1=_OCC_TWO_PI)),
                               UInt64[vid_rim,vid_rim])),
            (1,t_mer,cid_mer,
                _ModelOCCShape(1,(kind=:circle,
                    geometry=(occ=:circle,center=mer.loc,n=mer.axis,
                              X=mer.X,Y=mer.Y,r=rb,t0=0.0,t1=_OCC_TWO_PI)),
                               UInt64[vid_rim,vid_rim])),
            (-1,kernel.maxima[2]+1,wid,
                _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_eq,true),
                    (cid_mer,false),(cid_eq,false),(cid_mer,true)],),
                               UInt64[cid_eq,cid_mer,cid_eq,cid_mer])),
            (2,s_face,fid,
                _ModelOCCShape(2,(kind=:torus,
                    geometry=(occ=:torus,center=center,axis=n,X=X,Y=Y,
                              r1=ra,r2=rb,angle=a,
                              pcurves=Dict{Int,NamedTuple}(
                                  t_eq=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),
                                         rev=(lin2d=((0.0,_OCC_TWO_PI),
                                                     (1.0,0.0)),)),
                                  t_mer=>(fwd=(lin2d=((_OCC_TWO_PI,-0.0),
                                                      (0.0,1.0)),),
                                          rev=(lin2d=((0.0,-0.0),
                                                      (0.0,1.0)),))))),
                               UInt64[wid])),
            (-2,kernel.maxima[1]+1,shell_id,
                _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[(fid,false)],),
                               UInt64[fid])),
            (3,allocated,solid_id,
                _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id]))]
    else
        ids=_model_occ_reserve_shapes(kernel,13,caller)
        vid_end,vid_start,cid_arc,cid_end,cid_start,
            wid_torus,wid_start,wid_end,
            fid_torus,fid_start,fid_end,shell_id,solid_id=ids
        t_arc=kernel.maxima[4]+1; t_end=kernel.maxima[4]+2
        t_start=kernel.maxima[4]+3
        s_torus=kernel.maxima[5]+1; s_start=kernel.maxima[5]+2
        s_end=kernel.maxima[5]+3
        Ma,loca=_occ_trsf_rot(center,n,a)
        w_end=_occ_vertex(center,n,X,mpE,Ma,loca)
        w_start=ntuple(i->(center[i]+n[i]*mpE[2])+X[i]*mpE[1],3)
        mer_e=_occ_xform_ax2(Ma,loca,(loc=locm,axis=Dm,X=Xm,Y=Ym))
        mer_s=_occ_xform_ax2(
            ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)),(0.0,0.0,0.0),
            (loc=locm,axis=Dm,X=Xm,Y=Ym))
        cap_e=_occ_rotate_ax2(Ma,loca,(loc=center,axis=Dm,X=Xm,Y=Ym))
        start_pc=Dict{Int,NamedTuple}(
            t_start=>(fwd=(circ2d=((ra,0.0),(1.0,0.0),rb),),rev=nothing))
        end_pc=Dict{Int,NamedTuple}(
            t_end=>(fwd=(circ2d=((ra,0.0),(1.0,0.0),rb),),rev=nothing))
        trange=Dict(t_start=>(0.0,_OCC_TWO_PI),t_end=>(0.0,_OCC_TWO_PI))
        pending=[
            (0,kernel.maxima[3]+1,vid_end,
                _ModelOCCShape(0,(point=w_end,mesh_size=0.0),UInt64[])),
            (0,kernel.maxima[3]+2,vid_start,
                _ModelOCCShape(0,(point=w_start,mesh_size=0.0),UInt64[])),
            (1,t_arc,cid_arc,
                _ModelOCCShape(1,(kind=:circle,
                    geometry=(occ=:circle,center=ce,n=ne,X=Xe,Y=Ye,
                              r=mpE[1],t0=0.0,t1=a)),
                               UInt64[vid_start,vid_end])),
            (1,t_end,cid_end,
                _ModelOCCShape(1,(kind=:circle,
                    geometry=(occ=:circle,center=mer_e.loc,n=mer_e.axis,
                              X=mer_e.X,Y=mer_e.Y,r=rb,t0=0.0,
                              t1=_OCC_TWO_PI)),
                               UInt64[vid_end,vid_end])),
            (1,t_start,cid_start,
                _ModelOCCShape(1,(kind=:circle,
                    geometry=(occ=:circle,center=mer_s.loc,n=mer_s.axis,
                              X=mer_s.X,Y=mer_s.Y,r=rb,t0=0.0,
                              t1=_OCC_TWO_PI)),
                               UInt64[vid_start,vid_start])),
            (-1,kernel.maxima[2]+1,wid_torus,
                _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_arc,true),
                    (cid_start,true),(cid_arc,false),(cid_end,false)],),
                               UInt64[cid_arc,cid_start,cid_arc,cid_end])),
            (-1,kernel.maxima[2]+2,wid_start,
                _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_start,false)],),
                               UInt64[cid_start])),
            (-1,kernel.maxima[2]+3,wid_end,
                _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_end,false)],),
                               UInt64[cid_end])),
            (2,s_torus,fid_torus,
                _ModelOCCShape(2,(kind=:torus,
                    geometry=(occ=:torus,center=center,axis=n,X=X,Y=Y,
                              r1=ra,r2=rb,angle=a,
                              pcurves=Dict{Int,NamedTuple}(
                                  t_arc=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),
                                          rev=(lin2d=((0.0,_OCC_TWO_PI),
                                                      (1.0,0.0)),)),
                                  t_start=>(fwd=(lin2d=((0.0,-0.0),
                                                        (0.0,1.0)),),
                                            rev=nothing),
                                  t_end=>(fwd=(lin2d=((a,-0.0),(0.0,1.0)),),
                                          rev=nothing)))),
                               UInt64[wid_torus])),
            (2,s_start,fid_start,
                _ModelOCCShape(2,(kind=:plane,
                    geometry=(occ=:plane,center=center,axis=Dm,X=Xm,Y=Ym,
                              reversed=false,pcurves=start_pc,
                              uvb=_model_occ_prim_uvb(start_pc,trange))),
                               UInt64[wid_start])),
            (2,s_end,fid_end,
                _ModelOCCShape(2,(kind=:plane,
                    geometry=(occ=:plane,center=cap_e.loc,axis=cap_e.axis,
                              X=cap_e.X,Y=cap_e.Y,reversed=true,
                              pcurves=end_pc,
                              uvb=_model_occ_prim_uvb(end_pc,trange))),
                               UInt64[wid_end])),
            (-2,kernel.maxima[1]+1,shell_id,
                _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[(fid_torus,false),
                    (fid_start,false),(fid_end,true)],),
                               UInt64[fid_torus,fid_start,fid_end])),
            (3,allocated,solid_id,
                _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id]))]
    end
    for (dimension,ptag,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,ptag,id,shape)
    end
    kernel.next_shape=solid_id+UInt64(1)
    return allocated
end
