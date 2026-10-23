# OCC circle-edge fillets (rolling-ball torus blends on coaxial revolution
# faces), curved wireTags trims on tensor surfaces, and the geometric
# seam-period/topology queries used by trim and BREP unwrapping.
using Test
using Tessella
using Tessella.Model: _model_occ_fillet, _model_occ_add_cylinder!,
                      _model_occ_add_cone!, _model_occ_add_torus!,
                      _model_occ_add_point!, _model_occ_add_circle!,
                      _model_occ_add_wire!, _model_occ_add_bspline_surface!,
                      _model_occ_shape_tag, _occ_pcurve_d0,
                      _model_nurbs_surface_d2, _model_nurbs_surface_bounds,
                      _model_nurbs_patch_seam_period, _brep_uperiod,
                      model_value, model_parametrization_bounds, model_boundary

@testset "circle-edge fillet smoke" begin
    # OCC cylinder r=1, h=2: curves 1(top),2(seam),3(bottom). Fillet the
    # bottom rim r=0.25: expected rolling-ball torus centred at r=0.75,z=0.25
    # with tangent rings (r=1,z=0.25) on the lateral and (r=0.75,z=0) on cap.
    m=GeoModel()
    _model_occ_add_cylinder!(m,0,0,0,0,0,2,1)
    out=_model_occ_fillet(m,[1],[3],[0.25])
    @test out==[(3,1)]
    k=m.occ_kernel
    surfs=[g.geometry.occ for (d,t) in sort!(collect(keys(k.bindings)))
           if d==2 for g in (k.shapes[k.bindings[(d,t)]].geometry,)]
    @test sort!(surfs)==[:cylinder,:plane,:plane,:torus]
    sid=k.shapes[k.bindings[(3,1)]].children[1]
    uses=k.shapes[sid].geometry.uses
    @test length(uses)==4
    fid_t=only(fid for (fid,_) in uses
               if get(k.shapes[fid].geometry.geometry,:occ,nothing)===:torus)
    tg=k.shapes[fid_t].geometry.geometry
    @test tg.r1≈0.75 atol=1e-12
    @test tg.r2≈0.25
    @test tg.center[3]≈0.25 atol=1e-12
    rings=Tuple{Float64,Float64}[]
    for ((d,t),id) in k.bindings
        d==1 || continue
        g=k.shapes[id].geometry.geometry
        get(g,:occ,nothing)===:circle || continue
        push!(rings,(g.r,g.center[3]))
    end
    @test any(q->abs(q[1]-1.0)<1e-9 && abs(q[2]-0.25)<1e-9,rings)
    @test any(q->abs(q[1]-0.75)<1e-9 && abs(q[2])<1e-9,rings)
    @test !haskey(k.bindings,(1,3))
    fid_l=only(fid for (fid,_) in uses
               if get(k.shapes[fid].geometry.geometry,:occ,nothing)===:cylinder)
    lg=k.shapes[fid_l].geometry.geometry
    wl=k.shapes[fid_l].children[1]
    wuses=k.shapes[wl].geometry.uses
    @test length(wuses)==4
    occs=[get(k.shapes[u].geometry.geometry,:occ,nothing) for (u,_) in wuses]
    @test occs==[:circle,:line,:circle,:line]
    seam_ids=[u for (u,_) in wuses if
              get(k.shapes[u].geometry.geometry,:occ,nothing)===:line]
    sg=k.shapes[seam_ids[1]].geometry.geometry
    @test sg.t1-sg.t0≈1.75 atol=1e-9
    etags=Set(keys(lg.pcurves))
    @test length(etags)>=3
    @test_throws ArgumentError _model_occ_fillet(m,[1],[3],[-0.1])
end

@testset "circle-edge fillet — both rims / rejection" begin
    m=GeoModel()
    _model_occ_add_cylinder!(m,0,0,0,0,0,2,1)
    out=_model_occ_fillet(m,[1],[1,3],[0.25,0.25])
    @test out==[(3,1)]
    k=m.occ_kernel
    ntor=count(id->get(k.shapes[id].geometry.geometry,:occ,nothing)===:torus,
               (id for ((d,_),id) in k.bindings if d==2))
    @test ntor==2
    fid_l=only(id for ((d,_),id) in k.bindings if d==2 &&
               get(k.shapes[id].geometry.geometry,:occ,nothing)===:cylinder)
    wl=k.shapes[fid_l].children[1]
    wuses=k.shapes[wl].geometry.uses
    @test length(wuses)==4
    seam_ids=[u for (u,_) in wuses if
              get(k.shapes[u].geometry.geometry,:occ,nothing)===:line]
    sg=k.shapes[seam_ids[1]].geometry.geometry
    @test sg.t1-sg.t0≈1.5 atol=1e-9

    # Overrun: ρ=2.5 on a tall-radius short cylinder — the strip leg covers
    # the top rim, a genuine boundary-overrun rejection.
    m=GeoModel()
    _model_occ_add_cylinder!(m,0,0,0,0,0,2,3)
    @test_throws ArgumentError _model_occ_fillet(m,[1],[3],[2.5])

    # ρ > R puts the tube centre across the axis — spindle rejection.
    m=GeoModel()
    _model_occ_add_cylinder!(m,0,0,0,0,0,2,1)
    @test_throws ArgumentError _model_occ_fillet(m,[1],[3],[1.5])

    m=GeoModel()
    _model_occ_add_cylinder!(m,0,0,0,0,0,2,1)
    @test_throws ArgumentError _model_occ_fillet(m,[1],[99],[0.1])
end

@testset "circle-edge fillet — cone rim / self-shared rim" begin
    # Truncated cone r1=2 (z=0) → r2=1 (z=2). Fillet the wide bottom rim.
    m=GeoModel()
    _model_occ_add_cone!(m,0,0,0,0,0,2,2,1)
    out=_model_occ_fillet(m,[1],[3],[0.3])
    @test out==[(3,1)]
    k=m.occ_kernel
    fid_t=only(id for ((d,_),id) in k.bindings if d==2 &&
               get(k.shapes[id].geometry.geometry,:occ,nothing)===:torus)
    tg=k.shapes[fid_t].geometry.geometry
    @test tg.r2≈0.3
    @test 0.0<tg.r1<2.0
    @test tg.center[3]>0.0
    @test haskey(k.bindings,(1,1)) && !haskey(k.bindings,(1,3))

    # A lone torus' equator is shared by the face with itself — reject.
    m=GeoModel()
    _model_occ_add_torus!(m,0,0,0,3,1)
    @test_throws ArgumentError _model_occ_fillet(m,[1],[1],[0.2])
end

@testset "OCC wireTags curved trim on bspline surface" begin
    # Flat 3x3 quadratic patch over x,y ∈ [-1,2] at z=0; trim it with a
    # circle wire — the classic bounded-face curved trim (IGES 144 style).
    m=GeoModel()
    pts=Int[]
    for y in (-1.0,0.5,2.0), x in (-1.0,0.5,2.0)
        push!(pts,_model_occ_add_point!(m,x,y,0.0))
    end
    ct=_model_occ_add_circle!(m,0.5,0.5,0.0,0.4)
    wt=_model_occ_add_wire!(m,[ct])
    st=_model_occ_add_bspline_surface!(m,pts,3;wire_tags=[wt])
    k=m.occ_kernel
    fsh=k.shapes[k.bindings[(2,st)]]
    g=fsh.geometry.geometry
    @test get(g,:trimmed,false)===true
    @test fsh.children==[k.bindings[(-1,wt)]]
    @test haskey(g.pcurves,ct)
    # The projected pcurve must re-evaluate onto the 3D circle: the patch
    # is flat, so S(uv) sits in z=0 at radius 0.4 from (0.5,0.5).
    eg=k.shapes[k.bindings[(1,ct)]].geometry.geometry
    for i in 0:4
        t=eg.t0+(eg.t1-eg.t0)*i/4
        uv=_occ_pcurve_d0(g.pcurves[ct].fwd,t)
        P=first(_model_nurbs_surface_d2(g,uv...))
        @test P[3]≈0.0 atol=1e-9
        @test hypot(P[1]-0.5,P[2]-0.5)≈0.4 atol=1e-7
    end

    # Outer wire + hole: two children in declaration order.
    m=GeoModel()
    pts=Int[]
    for y in (-1.0,0.5,2.0), x in (-1.0,0.5,2.0)
        push!(pts,_model_occ_add_point!(m,x,y,0.0))
    end
    c1=_model_occ_add_circle!(m,0.5,0.5,0.0,0.4)
    c2=_model_occ_add_circle!(m,0.5,0.5,0.0,0.15)
    w1=_model_occ_add_wire!(m,[c1])
    w2=_model_occ_add_wire!(m,[c2])
    st=_model_occ_add_bspline_surface!(m,pts,3;wire_tags=[w1,w2])
    k=m.occ_kernel
    fsh=k.shapes[k.bindings[(2,st)]]
    @test fsh.children==[k.bindings[(-1,w1)],k.bindings[(-1,w2)]]
    @test haskey(fsh.geometry.geometry.pcurves,c2)

    # Rejections: unknown wire, duplicate wire, wire outside the patch.
    m=GeoModel()
    pts=Int[]
    for y in (-1.0,0.5,2.0), x in (-1.0,0.5,2.0)
        push!(pts,_model_occ_add_point!(m,x,y,0.0))
    end
    @test_throws ArgumentError _model_occ_add_bspline_surface!(
        m,pts,3;wire_tags=[77])
    c3=_model_occ_add_circle!(m,10.0,10.0,0.0,0.4)
    w3=_model_occ_add_wire!(m,[c3])
    @test_throws ArgumentError _model_occ_add_bspline_surface!(
        m,pts,3;wire_tags=[w3,w3])
    @test_throws ArgumentError _model_occ_add_bspline_surface!(
        m,pts,3;wire_tags=[w3])
end

@testset "NURBS tensor seam period / topology queries" begin
    # Geometrically closed-in-u clamped patch: the first and last pole
    # columns coincide (distinct tags) so the u-boundaries meet — the
    # geometric seam test must report the knot-domain width as period.
    m=GeoModel()
    cols=[(1.0,0.0),(0.0,1.0),(-1.0,0.0),(1.0,0.0)]
    pts=Int[]
    for z in (0.0,1.0), (x,y) in cols
        push!(pts,_model_occ_add_point!(m,x,y,z))
    end
    st=_model_occ_add_bspline_surface!(m,pts,4;degree_u=3,degree_v=1)
    k=m.occ_kernel
    g=k.shapes[k.bindings[(2,st)]].geometry.geometry
    lo,hi=_model_nurbs_surface_bounds(g)
    @test lo==[0.0,0.0] && hi==[1.0,1.0]
    @test _model_nurbs_patch_seam_period(g,1)≈1.0
    @test _model_nurbs_patch_seam_period(g,2)==0.0
    @test _brep_uperiod(g)≈1.0

    # An open patch has no seam in either direction.
    m=GeoModel()
    pts=Int[]
    for y in (0.0,0.5,1.0), x in (0.0,0.5,1.0)
        push!(pts,_model_occ_add_point!(m,x,y,0.0))
    end
    st=_model_occ_add_bspline_surface!(m,pts,3)
    k=m.occ_kernel
    g=k.shapes[k.bindings[(2,st)]].geometry.geometry
    @test _model_nurbs_patch_seam_period(g,1)==0.0
    @test _model_nurbs_patch_seam_period(g,2)==0.0
    @test _brep_uperiod(g)==0.0
end
