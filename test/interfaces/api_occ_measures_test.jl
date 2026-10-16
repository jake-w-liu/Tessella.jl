# `gmsh.model.occ` kernel-side BRepGProp measures — `get_mass` and
# `get_center_of_mass` on the bound kernel graph (unsynchronized, like
# upstream `OCC_Internals::getMass`/`getCenterOfMass`). Curve measures
# are per-record closed forms or Gauss-Legendre panel quadrature; faces
# use the plane shoelace with conic-segment corrections, quadric
# closed-form patch integrals, or NURBS tensor quadrature; solids use
# the compact `primitive` record, torus-sector Pappus integrals, or the
# divergence theorem over all-plane shells.
using Test

using Tessella
const A=Tessella.API

function _occ_setup()
    A.initialize()
    return A.model.occ
end

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "occ.get_mass / get_center_of_mass — curves" begin
    occ=_occ_setup()
    p1=occ.add_point(0.0,0.0,0.0); p2=occ.add_point(3.0,4.0,0.0)
    l=occ.add_line(p1,p2)
    @test occ.get_mass(1,l)≈5.0
    @test collect(occ.get_center_of_mass(1,l))≈[1.5,2.0,0.0]
    # Quarter arc: L=r·Δθ; centroid rides the bisector at r·sinc(Δθ/2).
    c0=occ.add_point(0.0,0.0,0.0)
    a=occ.add_point(1.0,0.0,0.0); b=occ.add_point(0.0,1.0,0.0)
    arc=occ.add_circle_arc(a,c0,b)
    @test occ.get_mass(1,arc)≈π/2
    @test collect(occ.get_center_of_mass(1,arc))≈[2/π,2/π,0.0] atol=1e-9
    # Full ellipse: quadrature recovers the elliptic-integral length and
    # the center centroid to panel accuracy.
    e=occ.add_ellipse(0.0,0.0,0.0,2.0,1.0)
    @test occ.get_mass(1,e)≈9.6884482245 atol=1e-6
    @test collect(occ.get_center_of_mass(1,e))≈[0.0,0.0,0.0] atol=1e-9
    # Half-ellipse: symmetric about the minor axis.
    ea=occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=0.0,angle2=Float64(π))
    @test occ.get_mass(1,ea)≈4.8442241103 atol=1e-6
    cx,cy,_=occ.get_center_of_mass(1,ea)
    @test cx≈0.0 atol=1e-9
    @test 0.0<cy<1.0
    # Spline curves measure by per-knot-span quadrature; the symmetric
    # control polygon puts the centroid on x=1.5.
    pts=[occ.add_point(x,y,0.0) for (x,y) in ((0,0),(1,2),(2,2),(3,0))]
    bs=occ.add_bspline(pts)
    @test occ.get_mass(1,bs)>4.0
    cx2,_,_=occ.get_center_of_mass(1,bs)
    @test cx2≈1.5 atol=1e-9
end

@testset "occ.get_mass / get_center_of_mass — faces" begin
    occ=_occ_setup()
    # Unit square: shoelace area and polygon centroid.
    p=[occ.add_point(x,y,0.0) for (x,y) in ((0,0),(2,0),(2,2),(0,2))]
    ls=[occ.add_line(p[i],p[mod1(i+1,4)]) for i in 1:4]
    w=occ.add_wire(ls)
    f=occ.add_plane_surface([w])
    @test occ.get_mass(2,f)≈4.0
    @test collect(occ.get_center_of_mass(2,f))≈[1.0,1.0,0.0]
    # Square minus centered circular hole — the conic-segment shoelace
    # path keeps area |outer|−|hole| and the symmetric centroid.
    cp=[occ.add_point(1.5,1.0,0.0),occ.add_point(1.0,1.0,0.0),
        occ.add_point(0.5,1.0,0.0)]
    c1=occ.add_circle_arc(cp[1],cp[2],cp[3])
    c2=occ.add_circle_arc(cp[3],cp[2],cp[1])
    wh=occ.add_wire([c1,c2])
    fh=occ.add_plane_surface([w,wh])
    @test occ.get_mass(2,fh)≈4.0-π/4 atol=1e-9
    @test collect(occ.get_center_of_mass(2,fh))≈[1.0,1.0,0.0] atol=1e-9
end

@testset "occ.get_mass / get_center_of_mass — quadric faces" begin
    occ=_occ_setup()
    occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,2.0,1.0)
    @test occ.get_mass(2,1)≈4.0π atol=1e-9   # lateral face
    @test collect(occ.get_center_of_mass(2,1))≈[0.0,0.0,1.0] atol=1e-9
    @test occ.get_mass(2,2)≈π atol=1e-9     # top cap
    @test collect(occ.get_center_of_mass(2,2))≈[0.0,0.0,2.0] atol=1e-9
    # Cone lateral A=π(r1+r2)·slant; centroid at h·(r1+2r2)/(3(r1+r2))
    # from the r1 base — the apex limit is h/3, cylinder limit h/2.
    occ.add_cone(0.0,0.0,0.0,0.0,0.0,3.0,2.0,1.0)
    @test occ.get_mass(2,4)≈3.0π*sqrt(10.0) atol=1e-9
    @test collect(occ.get_center_of_mass(2,4))≈
        [0.0,0.0,4.0/3] atol=1e-9
    occ.add_sphere(0.0,0.0,0.0,2.0)
    @test occ.get_mass(2,7)≈16.0π atol=1e-9
    @test collect(occ.get_center_of_mass(2,7))≈[0.0,0.0,0.0] atol=1e-9
    occ.add_torus(0.0,0.0,0.0,3.0,1.0)
    @test occ.get_mass(2,8)≈12.0π*π atol=1e-8
    @test collect(occ.get_center_of_mass(2,8))≈[0.0,0.0,0.0] atol=1e-9
end

@testset "occ.get_mass / get_center_of_mass — solids" begin
    occ=_occ_setup()
    b=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    @test occ.get_mass(3,b)≈1.0
    @test collect(occ.get_center_of_mass(3,b))≈[0.5,0.5,0.5]
    b2=occ.add_box(0.0,0.0,0.0,2.0,1.0,1.0)
    @test occ.get_mass(3,b2)≈2.0
    @test collect(occ.get_center_of_mass(3,b2))≈[1.0,0.5,0.5]
    c=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,2.0,1.0)
    @test occ.get_mass(3,c)≈2.0π atol=1e-12
    @test collect(occ.get_center_of_mass(3,c))≈[0.0,0.0,1.0] atol=1e-12
    s=occ.add_sphere(1.0,2.0,3.0,2.0)
    @test occ.get_mass(3,s)≈32.0π/3 atol=1e-12
    @test collect(occ.get_center_of_mass(3,s))≈[1.0,2.0,3.0] atol=1e-12
    # Frustum r1=2,r2=1,h=3: V=7π; z̄=h(r1²+2r1r2+3r2²)/(4(r1²+r1r2+r2²)).
    k=occ.add_cone(0.0,0.0,0.0,0.0,0.0,3.0,2.0,1.0)
    @test occ.get_mass(3,k)≈7.0π atol=1e-12
    @test collect(occ.get_center_of_mass(3,k))≈[0.0,0.0,33.0/28] atol=1e-12
    t=occ.add_torus(0.0,0.0,0.0,3.0,1.0)
    @test occ.get_mass(3,t)≈6.0π*π atol=1e-10
    @test collect(occ.get_center_of_mass(3,t))≈[0.0,0.0,0.0] atol=1e-9
end

@testset "occ.get_mass / get_center_of_mass — partial torus" begin
    occ=_occ_setup()
    # Sector φ=π/2, R=3, r=1: V=φ·πr²·R=3π²/2; centroid rides the
    # u=π/4 bisector at (R+r²/4R)·sinc(φ/2)=37/12·2√2/π.
    tp=occ.add_torus(0.0,0.0,0.0,3.0,1.0;angle=Float64(π)/2)
    @test occ.get_mass(3,tp)≈3.0π*π/2 atol=1e-10
    arm=(37/12)*2*sqrt(2)/π
    @test collect(occ.get_center_of_mass(3,tp))≈
        [arm*sqrt(2)/2,arm*sqrt(2)/2,0.0] atol=1e-9
end

@testset "occ.get_mass / get_center_of_mass — misc dims + lifecycle" begin
    occ=_occ_setup()
    # Points, wires and shells hit no GProp arm upstream → 0 + origin.
    p=occ.add_point(1.0,2.0,3.0)
    @test occ.get_mass(0,p)≈0.0
    @test collect(occ.get_center_of_mass(0,p))≈[0.0,0.0,0.0]
    p2=occ.add_point(2.0,0.0,0.0)
    l=occ.add_line(p,p2); w=occ.add_wire([l])
    @test occ.get_mass(-1,w)≈0.0
    @test collect(occ.get_center_of_mass(-1,w))≈[0.0,0.0,0.0]
    b=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    @test occ.get_mass(-2,1)≈0.0
    @test collect(occ.get_center_of_mass(-2,1))≈[0.0,0.0,0.0]
    # Empty volume: no shells → the empty-GProp result upstream reports.
    v=occ.add_volume(Int[])
    @test occ.get_mass(3,v)≈0.0
    @test collect(occ.get_center_of_mass(3,v))≈[0.0,0.0,0.0]
    # Degenerate pole edges: zero length, centroid on the pole.
    occ.add_sphere(0.0,0.0,0.0,2.0)
    m=A._with_model(identity)
    deg=sort([t for ((d,t),id) in m.occ_kernel.bindings
              if d==1 && m.occ_kernel.shapes[id].geometry.geometry.occ===
                  :degenerate])
    @test length(deg)==2
    @test occ.get_mass(1,deg[1])≈0.0
    cx,cy,cz=occ.get_center_of_mass(1,deg[1])
    @test abs(cz)≈2.0
    # Unknown entity and out-of-range dims reject.
    @test occursin("Unknown OpenCASCADE entity of dimension 3 with tag 99",
                   _errtext(()->occ.get_mass(3,99)))
    @test occursin("dimension",_errtext(()->occ.get_mass(4,1)))
end

@testset "occ.get_mass — kernel-side view across transforms" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,2.0,1.0)
    occ.translate([(3,v)],5.0,0.0,0.0)
    # Transform preserves volume and shifts the centroid — visible
    # kernel-side before synchronization, like upstream's registry.
    @test occ.get_mass(3,v)≈2.0π atol=1e-12
    @test collect(occ.get_center_of_mass(3,v))≈[5.0,0.0,1.0] atol=1e-12
    occ.synchronize()
    @test collect(occ.get_center_of_mass(3,v))≈[5.0,0.0,1.0] atol=1e-12
    # Mirror preserves volume.
    occ2=_occ_setup()
    b=occ2.add_box(2.0,0.0,0.0,1.0,1.0,1.0)
    occ2.mirror([(3,b)],1.0,0.0,0.0,0.0)
    @test occ2.get_mass(3,b)≈1.0
    @test collect(occ2.get_center_of_mass(3,b))≈[-2.5,0.5,0.5] atol=1e-12
end

@testset "occ.get_mass — conic-faced solid split" begin
    occ=_occ_setup()
    # Circular-boundary faces reach the divergence path for get_mass
    # (their shoelace handles conics); the centroid uses the closed-form
    # plane-face moment (n̂·c̄)·A·c̄ so conic wires compute too.
    d1=occ.add_disk(0.0,0.0,0.0,1.0,1.0)
    d2=occ.add_disk(0.0,0.0,3.0,1.0,1.0)
    sh=occ.add_surface_loop([d1,d2])
    v=occ.add_volume([sh])
    @test occ.get_mass(3,v)≈π atol=1e-9
    cx,cy,cz=occ.get_center_of_mass(3,v)
    @test cx≈0.0 atol=1e-12
    @test cy≈0.0 atol=1e-12
    @test cz≈2.25 atol=1e-12   # M=π·3·(0,0,3); V=π → z̄=9/4
end
