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
# Partial sweeps (`angle`/`angle1`/`angle2`/`angle3` off the full-rotation
# defaults) follow `BRepPrim_OneAxis`'s sector layout: the closed rim circles
# become trimmed arcs over `[0,a]`, each arc endpoint gains a radial edge to
# a cap-axis vertex, an axis edge joins the two axis vertices, and two
# meridian Plane faces close the sector. Profile endpoints sitting on the
# axis (cone apices, sphere poles) collapse to a shared vertex plus a
# degenerate edge spanning the sweep — apex caps and their radials vanish,
# and the side faces become triangles (both-axis-end spheres give bigons).
# Sector solids carry a `sector` measure record instead of `primitive`: they
# publish no compact `m.cylinders`/`m.spheres`/`m.cones` encoding, so volume
# meshing flows through the explicit boundary like Boolean outputs.

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
    a==_OCC_TWO_PI || return _model_occ_add_cylinder_sector!(
        m,kernel,allocated,base,axis,h,radius,a,caller)
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

# Partial `BRepPrim_OneAxis` cylinder — the sector layout `BRepPrim_Wedge`
# emits for `angle < 2*Pi`: the closed rim circles become arcs over [0,a],
# each arc endpoint gains a radial edge to a cap-axis vertex, an axis edge
# joins the cap centers, two sector Plane caps close the rims, and two
# meridian Plane faces close the sweep. Vertex order mirrors upstream (end
# rim, start rim, then axis points); the solid carries a `sector` measure
# record rather than the compact `primitive` encoding, so the boundary
# flows through the explicit-graph meshers like Boolean outputs.
function _model_occ_add_cylinder_sector!(m::GeoModel,kernel,allocated,base,
                                        axis,h,radius,a,caller)
    kernel.maxima[3]<=typemax(Int32)-6 || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-9 || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<=typemax(Int32)-5 || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-5 || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    n=_occ_dir(axis)
    X,Y=_occ_ax2_frame(n)
    Xc,Yc=_occ_ax2_setx(n,X)
    nc=_occ_crossed(Xc,Yc)
    Ma,loca=_occ_trsf_rot(base,n,a)
    top_center=_occ_translated(base,n,h)
    Xa=_occ_xform_dir(Ma,Xc)
    sa,ca=_gm_sincos(a)
    # Start-frame points take the raw build-frame arithmetic like the
    # torus's `w_start`; end points rotate it by `a`.
    w_ts=ntuple(i->(base[i]+n[i]*h)+X[i]*radius,3)
    w_bs=ntuple(i->base[i]+X[i]*radius,3)
    w_te=_occ_vertex(base,n,X,(radius,h),Ma,loca)
    w_be=_occ_vertex(base,n,X,(radius,0.0),Ma,loca)
    ids=_model_occ_reserve_shapes(kernel,27,caller)
    vid_te,vid_ts,vid_be,vid_bs,vid_tc,vid_bc=ids[1:6]
    cid_top,cid_seam_e,cid_bot,cid_seam_s,cid_rt0,cid_rta,cid_rba,cid_rb0,
        cid_axis=ids[7:15]
    wid_lat,wid_top,wid_bot,wid_s0,wid_sa=ids[16:20]
    fid_lat,fid_top,fid_bot,fid_s0,fid_sa=ids[21:25]
    shell_id,solid_id=ids[26:27]
    vb=kernel.maxima[3]; cb=kernel.maxima[4]
    wb=kernel.maxima[2]; sb=kernel.maxima[5]
    t_top=cb+1; t_seam_e=cb+2; t_bot=cb+3; t_seam_s=cb+4; t_rt0=cb+5
    t_rta=cb+6; t_rba=cb+7; t_rb0=cb+8; t_axis=cb+9
    s_lat=sb+1; s_top=sb+2; s_bot=sb+3; s_s0=sb+4; s_sa=sb+5
    lat_pcurves=Dict{Int,NamedTuple}(
        t_top=>(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing),
        t_seam_e=>(fwd=(lin2d=((a,0.0),(0.0,1.0)),),rev=nothing),
        t_bot=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing),
        t_seam_s=>(fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),rev=nothing))
    top_pc=Dict{Int,NamedTuple}(
        t_top=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),radius),),rev=nothing),
        t_rta=>(fwd=(lin2d=((0.0,0.0),(ca,sa)),),rev=nothing),
        t_rt0=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing))
    bot_pc=Dict{Int,NamedTuple}(
        t_bot=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),radius),),rev=nothing),
        t_rba=>(fwd=(lin2d=((0.0,0.0),(ca,sa)),),rev=nothing),
        t_rb0=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing))
    s0_pc=Dict{Int,NamedTuple}(
        t_axis=>(fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),rev=nothing),
        t_rt0=>(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing),
        t_seam_s=>(fwd=(lin2d=((radius,0.0),(0.0,1.0)),),rev=nothing),
        t_rb0=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing))
    sa_pc=Dict{Int,NamedTuple}(
        t_axis=>(fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),rev=nothing),
        t_rta=>(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing),
        t_seam_e=>(fwd=(lin2d=((radius,0.0),(0.0,1.0)),),rev=nothing),
        t_rba=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing))
    trange=Dict{Int,NTuple{2,Float64}}(
        t_top=>(0.0,a),t_bot=>(0.0,a),t_axis=>(0.0,h),
        t_seam_e=>(0.0,h),t_seam_s=>(0.0,h),
        t_rt0=>(0.0,radius),t_rta=>(0.0,radius),
        t_rba=>(0.0,radius),t_rb0=>(0.0,radius))
    pending=[
        (0,vb+1,vid_te,
            _ModelOCCShape(0,(point=w_te,mesh_size=0.0),UInt64[])),
        (0,vb+2,vid_ts,
            _ModelOCCShape(0,(point=w_ts,mesh_size=0.0),UInt64[])),
        (0,vb+3,vid_be,
            _ModelOCCShape(0,(point=w_be,mesh_size=0.0),UInt64[])),
        (0,vb+4,vid_bs,
            _ModelOCCShape(0,(point=w_bs,mesh_size=0.0),UInt64[])),
        (0,vb+5,vid_tc,
            _ModelOCCShape(0,(point=top_center,mesh_size=0.0),UInt64[])),
        (0,vb+6,vid_bc,
            _ModelOCCShape(0,(point=base,mesh_size=0.0),UInt64[])),
        (1,t_top,cid_top,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=top_center,n=nc,X=Xc,Y=Yc,
                          r=radius,t0=0.0,t1=a)),
                           UInt64[vid_ts,vid_te])),
        (1,t_seam_e,cid_seam_e,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=w_be,dir=n,t0=0.0,t1=h)),
                           UInt64[vid_be,vid_te])),
        (1,t_bot,cid_bot,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=base,n=nc,X=Xc,Y=Yc,
                          r=radius,t0=0.0,t1=a)),
                           UInt64[vid_bs,vid_be])),
        (1,t_seam_s,cid_seam_s,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=w_bs,dir=n,t0=0.0,t1=h)),
                           UInt64[vid_bs,vid_ts])),
        (1,t_rt0,cid_rt0,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=top_center,dir=Xc,t0=0.0,
                          t1=radius)),
                           UInt64[vid_tc,vid_ts])),
        (1,t_rta,cid_rta,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=top_center,dir=Xa,t0=0.0,
                          t1=radius)),
                           UInt64[vid_tc,vid_te])),
        (1,t_rba,cid_rba,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=base,dir=Xa,t0=0.0,
                          t1=radius)),
                           UInt64[vid_bc,vid_be])),
        (1,t_rb0,cid_rb0,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=base,dir=Xc,t0=0.0,
                          t1=radius)),
                           UInt64[vid_bc,vid_bs])),
        (1,t_axis,cid_axis,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=base,dir=n,t0=0.0,t1=h)),
                           UInt64[vid_bc,vid_tc])),
        (-1,wb+1,wid_lat,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_top,true),
                (cid_seam_s,true),(cid_bot,false),(cid_seam_e,false)],),
                           UInt64[cid_top,cid_seam_s,cid_bot,cid_seam_e])),
        (-1,wb+2,wid_top,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_top,false),
                (cid_rta,true),(cid_rt0,false)],),
                           UInt64[cid_top,cid_rta,cid_rt0])),
        (-1,wb+3,wid_bot,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_bot,false),
                (cid_rba,true),(cid_rb0,false)],),
                           UInt64[cid_bot,cid_rba,cid_rb0])),
        (-1,wb+4,wid_s0,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_axis,true),
                (cid_rt0,true),(cid_seam_s,false),(cid_rb0,false)],),
                           UInt64[cid_axis,cid_rt0,cid_seam_s,cid_rb0])),
        (-1,wb+5,wid_sa,
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(cid_axis,true),
                (cid_rta,true),(cid_seam_e,false),(cid_rba,false)],),
                           UInt64[cid_axis,cid_rta,cid_seam_e,cid_rba])),
        (2,s_lat,fid_lat,
            _ModelOCCShape(2,(kind=:cylinder,
                geometry=(occ=:cylinder,center=base,axis=n,X=X,Y=Y,
                          radius=radius,height=h,angle=a,
                          pcurves=lat_pcurves)),
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
        (2,s_s0,fid_s0,
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=base,axis=_occ_crossed(Xc,n),
                          X=Xc,Y=n,reversed=false,pcurves=s0_pc,
                          uvb=_model_occ_prim_uvb(s0_pc,trange))),
                           UInt64[wid_s0])),
        (2,s_sa,fid_sa,
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=base,axis=_occ_crossed(Xa,n),
                          X=Xa,Y=n,reversed=true,pcurves=sa_pc,
                          uvb=_model_occ_prim_uvb(sa_pc,trange))),
                           UInt64[wid_sa])),
        (-2,kernel.maxima[1]+1,shell_id,
            _ModelOCCShape(-2,(uses=Tuple{UInt64,Bool}[(fid_lat,false),
                (fid_top,false),(fid_bot,true),(fid_s0,false),
                (fid_sa,true)],),
                           UInt64[fid_lat,fid_top,fid_bot,fid_s0,fid_sa])),
        (3,allocated,solid_id,
            _ModelOCCShape(3,(kind=:solid,
                sector=(occ=:cylinder,center=base,axis=axis,X=Xc,Y=Yc,
                        radius=radius,height=h,angle=a)),
                UInt64[shell_id]))]
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
    a>_OCC_TWO_PI && throw(ArgumentError(
        "$caller: Cone angle should be at most 2*Pi"))
    apex_top=rb==0.0; apex_bot=ra==0.0
    a==_OCC_TWO_PI || return _model_occ_add_cone_sector!(
        m,kernel,allocated,base,axis,h,ra,rb,apex_top,apex_bot,a,caller)
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

# Partial `BRepPrim_OneAxis` cone — the sector layout for `angle < 2*Pi`.
# A frustum end keeps its rim arc plus radials and a sector cap; an apex end
# collapses to a shared vertex plus a degenerate edge spanning the sweep, so
# the side faces become triangles and no cap binds on that end. Like the
# cylinder sector the solid carries a `sector` measure record.
function _model_occ_add_cone_sector!(m::GeoModel,kernel,allocated,base,axis,
                                     h,ra,rb,apex_top,apex_bot,a,caller)
    nverts=apex_top||apex_bot ? 4 : 6
    nedges=apex_top||apex_bot ? 7 : 9
    ncaps=(apex_top ? 0 : 1)+(apex_bot ? 0 : 1)
    nfaces=1+ncaps+2
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
    n=_occ_dir(axis)
    X,Y=_occ_ax2_frame(n)
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
    Ma,loca=_occ_trsf_rot(base,n,a)
    sa3,ca3=_gm_sincos(a)
    Xa=_occ_xform_dir(Ma,Xc)
    da=_occ_xform_dir(Ma,d0)
    top_center=_occ_translated(base,n,mpT[2])
    # Start-frame points take the raw build-frame arithmetic like the
    # torus's `w_start`; end points rotate it by `a`. An apex collapses to
    # its axis point.
    w_te=_occ_vertex(base,n,X,mpT,Ma,loca)
    w_ts=ntuple(i->(base[i]+n[i]*mpT[2])+X[i]*mpT[1],3)
    w_be=_occ_vertex(base,n,X,mpB,Ma,loca)
    w_bs=ntuple(i->(base[i]+n[i]*mpB[2])+X[i]*mpB[1],3)
    vb=kernel.maxima[3]; cb=kernel.maxima[4]
    wb=kernel.maxima[2]; sb=kernel.maxima[5]
    # Edge tags in upstream order: top rim (or apex-degenerate), end seam,
    # bottom rim (or apex-degenerate), start seam, top radials, bottom
    # radials — each family present only for non-apex ends — then the axis.
    t_top=cb+1; t_seam_e=cb+2; t_bot=cb+3; t_seam_s=cb+4
    t_rt0=cb+5; t_rta=cb+6
    t_rba=cb+5+(apex_top ? 0 : 2); t_rb0=cb+6+(apex_top ? 0 : 2)
    t_axis=cb+nedges
    lat_pcurves=Dict{Int,NamedTuple}(
        t_top=>(fwd=(lin2d=((0.0,len),(1.0,0.0)),),rev=nothing),
        t_seam_e=>(fwd=(lin2d=((a,0.0),(0.0,1.0)),),rev=nothing),
        t_bot=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing),
        t_seam_s=>(fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),rev=nothing))
    lat_geometry=(occ=:cone,center=base,axis=n,X=Xc,Y=Yc,r1=ra,r2=rb,
                  height=h,angle=a,pcurves=lat_pcurves)
    trange=Dict{Int,NTuple{2,Float64}}(
        t_top=>(0.0,a),t_bot=>(0.0,a),t_axis=>(0.0,h),
        t_seam_e=>(0.0,len),t_seam_s=>(0.0,len))
    if !apex_top
        trange[t_rt0]=(0.0,mpT[1]); trange[t_rta]=(0.0,mpT[1])
    end
    if !apex_bot
        trange[t_rba]=(0.0,mpB[1]); trange[t_rb0]=(0.0,mpB[1])
    end
    ids=_model_occ_reserve_shapes(kernel,nverts+nedges+2*nfaces+2,caller)
    vid(k)=ids[k]; eid(k)=ids[nverts+k]
    wid(k)=ids[nverts+nedges+k]; fid(k)=ids[nverts+nedges+nfaces+k]
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    # Vertices: end/start rim verts per end (a shared apex vertex on an
    # apex end), then the axis points of non-apex ends.
    if apex_top
        v_apex=vid(1); v_be=vid(2); v_bs=vid(3); v_bc=vid(4)
        v_te=v_ts=v_apex; v_tc=v_apex
        push!(pending,(0,vb+1,v_apex,
            _ModelOCCShape(0,(point=w_te,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+2,v_be,
            _ModelOCCShape(0,(point=w_be,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+3,v_bs,
            _ModelOCCShape(0,(point=w_bs,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+4,v_bc,
            _ModelOCCShape(0,(point=base,mesh_size=0.0),UInt64[])))
    elseif apex_bot
        v_apex=vid(1); v_te=vid(2); v_ts=vid(3); v_tc=vid(4)
        v_be=v_bs=v_apex; v_bc=v_apex
        push!(pending,(0,vb+1,v_apex,
            _ModelOCCShape(0,(point=w_be,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+2,v_te,
            _ModelOCCShape(0,(point=w_te,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+3,v_ts,
            _ModelOCCShape(0,(point=w_ts,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+4,v_tc,
            _ModelOCCShape(0,(point=top_center,mesh_size=0.0),UInt64[])))
    else
        v_te,v_ts,v_be,v_bs,v_tc,v_bc=vid.([1,2,3,4,5,6])
        push!(pending,(0,vb+1,v_te,
            _ModelOCCShape(0,(point=w_te,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+2,v_ts,
            _ModelOCCShape(0,(point=w_ts,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+3,v_be,
            _ModelOCCShape(0,(point=w_be,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+4,v_bs,
            _ModelOCCShape(0,(point=w_bs,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+5,v_tc,
            _ModelOCCShape(0,(point=top_center,mesh_size=0.0),UInt64[])))
        push!(pending,(0,vb+6,v_bc,
            _ModelOCCShape(0,(point=base,mesh_size=0.0),UInt64[])))
    end
    c_top=eid(1); c_seam_e=eid(2); c_bot=eid(3); c_seam_s=eid(4)
    c_rt0=apex_top ? UInt64(0) : eid(5); c_rta=apex_top ? UInt64(0) : eid(6)
    c_rba=apex_top ? eid(5) : eid(7); c_rb0=apex_top ? eid(6) : eid(8)
    c_axis=eid(nedges)
    if apex_top
        push!(pending,(1,t_top,c_top,
            _ModelOCCShape(1,(kind=:degenerate,
                geometry=(occ=:degenerate,t0=0.0,t1=a,face=lat_geometry)),
                           UInt64[v_apex,v_apex])))
    else
        push!(pending,(1,t_top,c_top,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=top_center,n=nc,X=Xc,Y=Yc,
                          r=mpT[1],t0=0.0,t1=a)),
                           UInt64[v_ts,v_te])))
    end
    push!(pending,(1,t_seam_e,c_seam_e,
        _ModelOCCShape(1,(kind=:line,
            geometry=(occ=:line,origin=w_be,dir=da,t0=0.0,t1=len)),
                       UInt64[v_be,v_te])))
    if apex_bot
        push!(pending,(1,t_bot,c_bot,
            _ModelOCCShape(1,(kind=:degenerate,
                geometry=(occ=:degenerate,t0=0.0,t1=a,face=lat_geometry)),
                           UInt64[v_apex,v_apex])))
    else
        push!(pending,(1,t_bot,c_bot,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=base,n=nc,X=Xc,Y=Yc,
                          r=mpB[1],t0=0.0,t1=a)),
                           UInt64[v_bs,v_be])))
    end
    push!(pending,(1,t_seam_s,c_seam_s,
        _ModelOCCShape(1,(kind=:line,
            geometry=(occ=:line,origin=w_bs,dir=d0,t0=0.0,t1=len)),
                       UInt64[v_bs,v_ts])))
    if !apex_top
        push!(pending,(1,t_rt0,c_rt0,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=top_center,dir=Xc,t0=0.0,
                          t1=mpT[1])),
                           UInt64[v_tc,v_ts])))
        push!(pending,(1,t_rta,c_rta,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=top_center,dir=Xa,t0=0.0,
                          t1=mpT[1])),
                           UInt64[v_tc,v_te])))
    end
    if !apex_bot
        push!(pending,(1,t_rba,c_rba,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=base,dir=Xa,t0=0.0,
                          t1=mpB[1])),
                           UInt64[v_bc,v_be])))
        push!(pending,(1,t_rb0,c_rb0,
            _ModelOCCShape(1,(kind=:line,
                geometry=(occ=:line,origin=base,dir=Xc,t0=0.0,
                          t1=mpB[1])),
                           UInt64[v_bc,v_bs])))
    end
    push!(pending,(1,t_axis,c_axis,
        _ModelOCCShape(1,(kind=:line,
            geometry=(occ=:line,origin=base,dir=n,t0=0.0,t1=h)),
                       UInt64[v_bc,v_tc])))
    # Faces: lateral, top cap, bottom cap, start side, end side — apex ends
    # drop their cap. Wire/face ordinal tracks the face position.
    push!(pending,(-1,wb+1,wid(1),
        _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(c_top,true),
            (c_seam_s,true),(c_bot,false),(c_seam_e,false)],),
                       UInt64[c_top,c_seam_s,c_bot,c_seam_e])))
    push!(pending,(2,sb+1,fid(1),
        _ModelOCCShape(2,(kind=:cone,geometry=lat_geometry),
                       UInt64[wid(1)])))
    shell_uses=Tuple{UInt64,Bool}[(fid(1),false)]
    nf=1
    if !apex_top
        top_pc=Dict{Int,NamedTuple}(
            t_top=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),mpT[1]),),rev=nothing),
            t_rta=>(fwd=(lin2d=((0.0,0.0),(ca3,sa3)),),rev=nothing),
            t_rt0=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing))
        nf+=1
        push!(pending,(-1,wb+nf,wid(nf),
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(c_top,false),
                (c_rta,true),(c_rt0,false)],),
                           UInt64[c_top,c_rta,c_rt0])))
        push!(pending,(2,sb+nf,fid(nf),
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=top_center,axis=n,X=Xc,Y=Yc,
                          reversed=false,pcurves=top_pc,
                          uvb=_model_occ_prim_uvb(top_pc,trange))),
                           UInt64[wid(nf)])))
        push!(shell_uses,(fid(nf),false))
    end
    if !apex_bot
        bot_pc=Dict{Int,NamedTuple}(
            t_bot=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),mpB[1]),),rev=nothing),
            t_rba=>(fwd=(lin2d=((0.0,0.0),(ca3,sa3)),),rev=nothing),
            t_rb0=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing))
        nf+=1
        push!(pending,(-1,wb+nf,wid(nf),
            _ModelOCCShape(-1,(uses=Tuple{UInt64,Bool}[(c_bot,false),
                (c_rba,true),(c_rb0,false)],),
                           UInt64[c_bot,c_rba,c_rb0])))
        push!(pending,(2,sb+nf,fid(nf),
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=base,axis=n,X=Xc,Y=Yc,
                          reversed=true,pcurves=bot_pc,
                          uvb=_model_occ_prim_uvb(bot_pc,trange))),
                           UInt64[wid(nf)])))
        push!(shell_uses,(fid(nf),true))
    end
    # Meridian side faces — a quad for two rim ends, a triangle through the
    # apex for an apex end. The frame is (radial direction, axis).
    for u_end in (false,true)
        Xu=u_end ? Xa : Xc
        c_seam_u=u_end ? c_seam_e : c_seam_s
        t_seam_u=u_end ? t_seam_e : t_seam_s
        c_rt_u=u_end ? c_rta : c_rt0
        t_rt_u=u_end ? t_rta : t_rt0
        c_rb_u=u_end ? c_rba : c_rb0
        t_rb_u=u_end ? t_rba : t_rb0
        rev=u_end
        side_pc=Dict{Int,NamedTuple}(
            t_axis=>(fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),rev=nothing),
            t_seam_u=>(fwd=(lin2d=((mpB[1],0.0),dir2d),),rev=nothing))
        side_uses=Tuple{UInt64,Bool}[(c_axis,true)]
        side_children=UInt64[c_axis]
        if apex_top
            side_pc[t_rb_u]=(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing)
            append!(side_uses,
                Tuple{UInt64,Bool}[(c_seam_u,false),(c_rb_u,false)])
            append!(side_children,UInt64[c_seam_u,c_rb_u])
        elseif apex_bot
            side_pc[t_rt_u]=(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing)
            append!(side_uses,
                Tuple{UInt64,Bool}[(c_rt_u,true),(c_seam_u,false)])
            append!(side_children,UInt64[c_rt_u,c_seam_u])
        else
            side_pc[t_rt_u]=(fwd=(lin2d=((0.0,h),(1.0,0.0)),),rev=nothing)
            side_pc[t_rb_u]=(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing)
            append!(side_uses,
                Tuple{UInt64,Bool}[(c_rt_u,true),(c_seam_u,false),
                    (c_rb_u,false)])
            append!(side_children,UInt64[c_rt_u,c_seam_u,c_rb_u])
        end
        nf+=1
        push!(pending,(-1,wb+nf,wid(nf),
            _ModelOCCShape(-1,(uses=side_uses,),side_children)))
        push!(pending,(2,sb+nf,fid(nf),
            _ModelOCCShape(2,(kind=:plane,
                geometry=(occ=:plane,center=base,axis=_occ_crossed(Xu,n),
                          X=Xu,Y=n,reversed=rev,pcurves=side_pc,
                          uvb=_model_occ_prim_uvb(side_pc,trange))),
                           UInt64[wid(nf)])))
        push!(shell_uses,(fid(nf),rev))
    end
    shell_id=ids[nverts+nedges+2*nfaces+1]
    solid_id=ids[nverts+nedges+2*nfaces+2]
    push!(pending,(-2,kernel.maxima[1]+1,shell_id,
        _ModelOCCShape(-2,(uses=shell_uses,),
                       UInt64[u[1] for u in shell_uses])))
    push!(pending,(3,allocated,solid_id,
        _ModelOCCShape(3,(kind=:solid,
            sector=(occ=:cone,center=base,axis=axis,X=Xc,Y=Yc,r1=ra,r2=rb,
                    height=h,angle=a)),
            UInt64[shell_id])))
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
    (a1>=-π/2 && a1<=π/2 && a2>=-π/2 && a2<=π/2) || throw(ArgumentError(
        "$caller: angle1 and angle2 should be in [-Pi/2, Pi/2]"))
    a1<a2 || throw(ArgumentError(
        "$caller: angle1 should be smaller than angle2"))
    (a1==-π/2 && a2==π/2 && a3==_OCC_TWO_PI) || return _model_occ_add_sphere_sector!(
        m,kernel,allocated,center,r,a1,a2,a3,caller)
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
# Partial `BRepPrim_OneAxis` sphere — `angle1`/`angle2` trim the latitude
# band, `angle3` sweeps the azimuth over `[0,a3]`. A non-pole latitude
# bound carries a latitude arc (a closed circle on a seam vertex for a
# full azimuth; a trimmed arc plus a cap-sector Plane and radial edges to
# a cap-axis vertex for a partial azimuth), while a pole bound collapses
# to a vertex plus a degenerate edge spanning the sweep. The axis edge
# and the two meridian Plane faces exist only for `angle3 < 2*Pi`. The
# solid carries a `sector` measure record — no compact `m.spheres`
# encoding — so volume meshing uses the explicit boundary like Boolean
# outputs.
function _model_occ_add_sphere_sector!(m,kernel,allocated,center,r,
                                       a1,a2,a3,caller)
    n=(0.0,0.0,1.0)
    X,Y=_occ_ax2_setx(n,(1.0,0.0,0.0))
    pole_top=a2>=π/2; pole_bot=a1<=-π/2
    azim=a3<_OCC_TWO_PI
    s1,c1=_gm_sincos(a1); s2,c2=_gm_sincos(a2)
    z1=r*s1; z2=r*s2; rr1=r*c1; rr2=r*c2
    M0,loc0=_occ_end_trsf(center,n)
    Ma,loca=_occ_trsf_rot(center,n,a3)
    mp(v)=_occ_circ2d_eval((0.0,0.0),(1.0,0.0),r,v)
    ws(v)=ntuple(i->(center[i]+n[i]*mp(v)[2])+X[i]*mp(v)[1],3)
    we(v)=_occ_vertex(center,n,X,mp(v),Ma,loca)
    Da=_occ_xform_dir(Ma,X)
    nc0=_occ_crossed(X,n); nca=_occ_crossed(Da,n)
    sa3,ca3=_gm_sincos(a3)
    nverts=azim ? 6-2*Int(pole_top)-2*Int(pole_bot) : 2
    nedges=azim ? 9-2*Int(pole_top)-2*Int(pole_bot) : 3
    ncaps=Int(!pole_top)+Int(!pole_bot)
    nfaces=1+ncaps+(azim ? 2 : 0)
    kernel.maxima[3]<=typemax(Int32)-nverts || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-nedges || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[5]<=typemax(Int32)-nfaces || throw(ArgumentError(
        "$caller: OCC surface tag allocation exceeds Int32"))
    kernel.maxima[1]<=typemax(Int32)-1 || throw(ArgumentError(
        "$caller: OCC surface loop tag allocation exceeds Int32"))
    ids=_model_occ_reserve_shapes(kernel,nverts+nedges+2*nfaces+2,caller)
    vb=kernel.maxima[3]; cb=kernel.maxima[4]
    wb=kernel.maxima[2]; sb=kernel.maxima[5]
    vid(k)=ids[k]; eid(k)=ids[nverts+k]
    wid(k)=ids[nverts+nedges+k]; fid(k)=ids[nverts+nedges+nfaces+k]
    vids=Dict{Symbol,UInt64}(); cids=Dict{Symbol,UInt64}()
    etags=Dict{Symbol,Int}(); trange=Dict{Int,NTuple{2,Float64}}()
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    # Vertices — non-pole rim endpoints (start before end), poles, then
    # cap-axis points.
    nv=0
    function emit_vertex!(key,p)
        nv+=1; vids[key]=vid(nv)
        push!(pending,(0,vb+nv,vids[key],
            _ModelOCCShape(0,(point=p,mesh_size=0.0),UInt64[])))
    end
    pole_bot || emit_vertex!(:u0b,ws(a1))
    (azim && !pole_bot) && emit_vertex!(:uab,we(a1))
    pole_top || emit_vertex!(:u0t,ws(a2))
    (azim && !pole_top) && emit_vertex!(:uat,we(a2))
    pole_bot && emit_vertex!(:ps,_occ_translated(center,n,-r))
    pole_top && emit_vertex!(:pn,_occ_translated(center,n,r))
    (azim && !pole_bot) && emit_vertex!(:cb,_occ_translated(center,n,z1))
    (azim && !pole_top) && emit_vertex!(:ct,_occ_translated(center,n,z2))
    vrim0(bot)=bot ? (pole_bot ? vids[:ps] : vids[:u0b]) :
                     (pole_top ? vids[:pn] : vids[:u0t])
    vrima(bot)=azim ?
        (bot ? (pole_bot ? vids[:ps] : vids[:uab]) :
               (pole_top ? vids[:pn] : vids[:uat])) : vrim0(bot)
    vax(bot)=bot ? (pole_bot ? vids[:ps] : vids[:cb]) :
                   (pole_top ? vids[:pn] : vids[:ct])
    zax(bot)=bot ? (pole_bot ? -r : z1) : (pole_top ? r : z2)
    # Edge keys and tag assignment — tags precede records so degenerate
    # edges can reference the finished sphere face.
    ekeys=Symbol[pole_bot ? :degb : :latb, azim ? :mer0 : :seam,
                 pole_top ? :degt : :latt]
    azim && push!(ekeys,:mera)
    if azim
        !pole_bot && append!(ekeys,(:rb0,:rba))
        !pole_top && append!(ekeys,(:rt0,:rta))
        push!(ekeys,:axis)
    end
    for (i,key) in enumerate(ekeys)
        etags[key]=cb+i; cids[key]=eid(i)
    end
    id2tag=Dict{UInt64,Int}(cids[k]=>etags[k] for k in ekeys)
    sph_pc=Dict{Int,NamedTuple}()
    for (bot,pole) in ((true,pole_bot),(false,pole_top))
        key=bot ? (pole ? :degb : :latb) : (pole ? :degt : :latt)
        sph_pc[etags[key]]=
            (fwd=(lin2d=((0.0,bot ? a1 : a2),(1.0,0.0)),),rev=nothing)
    end
    if azim
        sph_pc[etags[:mer0]]=
            (fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),rev=nothing)
        sph_pc[etags[:mera]]=
            (fwd=(lin2d=((a3,0.0),(0.0,1.0)),),rev=nothing)
    else
        sph_pc[etags[:seam]]=
            (fwd=(lin2d=((_OCC_TWO_PI,0.0),(0.0,1.0)),),
             rev=(lin2d=((0.0,0.0),(0.0,1.0)),))
    end
    sph_geometry=(occ=:sphere,center=center,radius=r,axis=n,X=X,Y=Y,
                  angle=a3,lat1=a1,lat2=a2,pcurves=sph_pc)
    function edge_geometry(key)
        if key===:latb || key===:latt
            bot=key===:latb
            (occ=:circle,
             center=_occ_translated(center,n,bot ? z1 : z2),n=n,X=X,Y=Y,
             r=bot ? rr1 : rr2,t0=0.0,t1=azim ? a3 : _OCC_TWO_PI)
        elseif key===:degb || key===:degt
            (occ=:degenerate,t0=0.0,t1=azim ? a3 : _OCC_TWO_PI,
             face=sph_geometry)
        elseif key===:mer0 || key===:seam
            (occ=:circle,center=center,n=nc0,X=X,Y=n,r=r,t0=a1,t1=a2)
        elseif key===:mera
            (occ=:circle,center=center,n=nca,X=Da,Y=n,r=r,t0=a1,t1=a2)
        elseif key===:axis
            (occ=:line,origin=center,dir=n,t0=zax(true),t1=zax(false))
        else
            bot=key===:rb0 || key===:rba
            (occ=:line,
             origin=_occ_translated(center,n,bot ? z1 : z2),
             dir=key===:rb0 || key===:rt0 ? X : Da,t0=0.0,
             t1=bot ? rr1 : rr2)
        end
    end
    edge_children(key)=
        key===:latb ? (vrim0(true),vrima(true)) :
        key===:latt ? (vrim0(false),vrima(false)) :
        key===:degb ? (vrim0(true),vrim0(true)) :
        key===:degt ? (vrim0(false),vrim0(false)) :
        key===:mera ? (vrima(true),vrima(false)) :
        key∈(:mer0,:seam) ? (vrim0(true),vrim0(false)) :
        key===:axis ? (vax(true),vax(false)) :
        key===:rb0 ? (vax(true),vrim0(true)) :
        key===:rba ? (vax(true),vrima(true)) :
        key===:rt0 ? (vax(false),vrim0(false)) :
                     (vax(false),vrima(false))
    for key in ekeys
        g=edge_geometry(key); a,b=edge_children(key)
        trange[etags[key]]=(get(g,:t0,0.0),get(g,:t1,0.0))
        kind=g.occ===:degenerate ? :degenerate :
            g.occ===:circle ? :circle : :line
        push!(pending,(1,etags[key],cids[key],
            _ModelOCCShape(1,(kind=kind,geometry=g),UInt64[a,b])))
    end
    # Wires — vertex-chained (edge,forward) multisets matching the public
    # layouts — then faces over them, in sphere/caps/sides order.
    nf=0
    shell_uses=Tuple{UInt64,Bool}[]
    function emit_face!(wire_uses,face_geometry,shell_rev)
        nf+=1
        push!(pending,(-1,wb+nf,wid(nf),
            _ModelOCCShape(-1,(uses=wire_uses,),
                           UInt64[u[1] for u in wire_uses])))
        push!(pending,(2,sb+nf,fid(nf),
            _ModelOCCShape(2,(kind=face_geometry.occ,
                geometry=face_geometry),UInt64[wid(nf)])))
        push!(shell_uses,(fid(nf),shell_rev))
    end
    fwd(k)=(cids[k],false); rev(k)=(cids[k],true)
    if azim
        sph_uses=if pole_top&&pole_bot
            Tuple{UInt64,Bool}[fwd(:degb),fwd(:mera),rev(:degt),
                               rev(:mer0)]
        elseif pole_top
            Tuple{UInt64,Bool}[rev(:degt),rev(:mer0),fwd(:latb),
                               fwd(:mera)]
        elseif pole_bot
            Tuple{UInt64,Bool}[rev(:latt),rev(:mer0),fwd(:degb),
                               fwd(:mera)]
        else
            Tuple{UInt64,Bool}[fwd(:latb),fwd(:mera),rev(:latt),
                               rev(:mer0)]
        end
    else
        sph_uses=Tuple{UInt64,Bool}[
            rev(pole_top ? :degt : :latt),rev(:seam),
            fwd(pole_bot ? :degb : :latb),fwd(:seam)]
    end
    emit_face!(sph_uses,sph_geometry,false)
    # Cap sector/disks at non-pole bounds.
    for bot in (true,false)
        (bot ? pole_bot : pole_top) && continue
        rr=bot ? rr1 : rr2; zi=bot ? z1 : z2
        cc=_occ_translated(center,n,zi)
        latk=bot ? :latb : :latt
        if azim
            rk0=bot ? :rb0 : :rt0; rka=bot ? :rba : :rta
            cap_pc=Dict{Int,NamedTuple}(
                etags[latk]=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),rr),),
                              rev=nothing),
                etags[rka]=>(fwd=(lin2d=((0.0,0.0),(ca3,sa3)),),
                             rev=nothing),
                etags[rk0]=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),
                             rev=nothing))
            emit_face!(Tuple{UInt64,Bool}[fwd(latk),rev(rka),fwd(rk0)],
                (occ=:plane,center=cc,axis=n,X=X,Y=Y,reversed=bot,
                 pcurves=cap_pc,uvb=_model_occ_prim_uvb(cap_pc,trange)),
                bot)
        else
            cap_pc=Dict{Int,NamedTuple}(
                etags[latk]=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),rr),),
                              rev=nothing))
            emit_face!(Tuple{UInt64,Bool}[fwd(latk)],
                (occ=:plane,center=cc,axis=n,X=X,Y=Y,reversed=bot,
                 pcurves=cap_pc,uvb=_model_occ_prim_uvb(cap_pc,trange)),
                bot)
        end
    end
    # Meridian side planes — quads for two rim bounds, triangles through
    # a pole, bigons for two poles. The u=0 plane's outward normal is
    # X×n = -Y; at u=a it is n×Da.
    if azim
        for u_end in (false,true)
            Xu=u_end ? Da : X; mkey=u_end ? :mera : :mer0
            side_pc=Dict{Int,NamedTuple}(
                etags[mkey]=>(fwd=(circ2d=((0.0,0.0),(1.0,0.0),r),),
                              rev=nothing),
                etags[:axis]=>(fwd=(lin2d=((0.0,0.0),(0.0,1.0)),),
                               rev=nothing))
            pole_bot || (side_pc[etags[u_end ? :rba : :rb0]]=
                (fwd=(lin2d=((0.0,z1),(1.0,0.0)),),rev=nothing))
            pole_top || (side_pc[etags[u_end ? :rta : :rt0]]=
                (fwd=(lin2d=((0.0,z2),(1.0,0.0)),),rev=nothing))
            uses=if pole_top&&pole_bot
                Tuple{UInt64,Bool}[fwd(mkey),rev(:axis)]
            elseif pole_top
                u_end ?
                    Tuple{UInt64,Bool}[fwd(:mera),rev(:axis),fwd(:rba)] :
                    Tuple{UInt64,Bool}[fwd(:mer0),rev(:axis),fwd(:rb0)]
            elseif pole_bot
                u_end ?
                    Tuple{UInt64,Bool}[fwd(:mera),rev(:rta),rev(:axis)] :
                    Tuple{UInt64,Bool}[fwd(:mer0),rev(:rt0),rev(:axis)]
            else
                u_end ?
                    Tuple{UInt64,Bool}[fwd(:mera),rev(:rta),rev(:axis),
                        fwd(:rba)] :
                    Tuple{UInt64,Bool}[fwd(:mer0),rev(:rt0),rev(:axis),
                        fwd(:rb0)]
            end
            emit_face!(uses,
                (occ=:plane,center=center,axis=_occ_crossed(Xu,n),X=Xu,
                 Y=n,reversed=u_end,pcurves=side_pc,
                 uvb=_model_occ_prim_uvb(side_pc,trange)),u_end)
        end
    end
    shell_id=ids[nverts+nedges+2*nfaces+1]
    solid_id=ids[nverts+nedges+2*nfaces+2]
    push!(pending,(-2,kernel.maxima[1]+1,shell_id,
        _ModelOCCShape(-2,(uses=shell_uses,),
                       UInt64[u[1] for u in shell_uses])))
    push!(pending,(3,allocated,solid_id,
        _ModelOCCShape(3,(kind=:solid,
            sector=(occ=:sphere,center=center,axis=n,X=X,Y=Y,radius=r,
                    angle=a3,lat1=a1,lat2=a2)),
            UInt64[shell_id])))
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
