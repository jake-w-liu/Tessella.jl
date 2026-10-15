# `model.occ.revolve` — `OCC_Internals::revolve` rotational sweeps
# (`BRepPrimAPI_MakeRevol`) on the pending-shape kernel. Bounded
# coverage: at angle 2·Pi, points sweep rim circles, line edges sweep
# cylinder/cone lateral records or annulus/disk plane faces, a
# single-wire all-line plane face in an axial plane sweeps a revolved
# solid (compact cylinder/cone encodings when the boundary matches a
# primitive layout, the `:revolved` Pappus encoding otherwise), and a
# full-circle disk pierced by an in-plane axis reuses the torus
# publication. Output ordering follows upstream's `getReturnedShapes`
# reorder — `[top (dim-1), body (dim), laterals …]` per input, grouped
# by ascending dimension; at the closed angle the rotated top
# coincides with the profile so the source entity answers as top.
# Rim circles/meridian faces deduplicate across shared profile
# sub-entities; swept rim circles bind to the source vertices.
# `numElements`/`heights`/`recombine` transport onto the Julia-native
# Layers machinery (the upstream `ExtrudeMesh` attachment). Partial
# angles sweep sector topology: rim arcs, rotated "top" copies and the
# two meridian caps bound the wedge. Non-axial profiles, axis-crossing
# or multi-wire faces, curved boundaries and solid inputs reject
# loudly in phase 1 before any bind; under Layers an axis-touching
# profile and a sub-three-subdivision full turn reject too — the
# former collapses the column machinery, the latter welds coincident
# cells onto the source row.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

_current_model()=A._with_model(identity)
_kernel()=A._with_model(m->m.occ_kernel)
_errtext(f)=(try;f();catch err;return sprint(showerror,err);end;"")

function _occ_counts()
    d=Dict{Int,Int}()
    for (dim,tag) in A.model.occ.get_entities()
        d[dim]=get(d,dim,0)+1
    end
    return d
end

# A rectangle in the xz meridian plane spanning rho_lo..rho_hi,
# h_lo..h_hi (all four corner radii >= 0).
function _meridian_rect(occ,rho_lo,rho_hi,h_lo,h_hi)
    ps=[occ.add_point(rho_lo,0.0,h_lo),occ.add_point(rho_hi,0.0,h_lo),
        occ.add_point(rho_hi,0.0,h_hi),occ.add_point(rho_lo,0.0,h_hi)]
    es=[occ.add_line(ps[1],ps[2]),occ.add_line(ps[2],ps[3]),
        occ.add_line(ps[3],ps[4]),occ.add_line(ps[4],ps[1])]
    return ps,es,occ.add_plane_surface([occ.add_wire(es)])
end

@testset "occ.revolve point and line inputs" begin
    occ=_occ_setup()
    p=occ.add_point(3.0,0.0,1.0)
    out=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    @test length(out)==2
    @test out[1]==(0,p)
    @test out[2][1]==1
    kernel=_kernel()
    eid=kernel.bindings[1,out[2][2]]
    eg=kernel.shapes[eid].geometry.geometry
    @test eg.occ===:circle
    @test eg.r≈3.0
    @test collect(eg.center)≈[0.0,0.0,1.0]
    @test collect(eg.n)≈[0.0,0.0,1.0]
    @test eg.t1-eg.t0==2pi
    # The rim circle binds both ends to the source vertex.
    a_id,b_id=kernel.shapes[eid].children
    @test a_id==b_id==kernel.bindings[0,p]

    # Axial line at rho=2, h in [0,3]: cylinder lateral + two rims.
    pa=occ.add_point(2.0,0.0,0.0)
    pb=occ.add_point(2.0,0.0,3.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    @test length(out)==4
    @test out[1]==(1,l)
    @test out[2][1]==2
    @test out[3][1]==1 && out[4][1]==1
    fid=kernel.bindings[2,out[2][2]]
    fg=kernel.shapes[fid].geometry
    @test fg.kind===:cylinder
    @test fg.geometry.radius≈2.0
    @test fg.geometry.height≈3.0
    @test collect(fg.geometry.center)≈[0.0,0.0,0.0]
    # The seam is the source edge; the wire has the canonical
    # [rim_hi, seam, rim_lo, seam] shape.
    uses=kernel.shapes[only(kernel.shapes[fid].children)].geometry.uses
    @test length(uses)==4
    @test uses[2][1]==kernel.bindings[1,l]
    @test uses[4][1]==kernel.bindings[1,l]
    occ.synchronize()
    @test occ.get_mass(out[2]...)≈2pi*2*3
    # Rim circles bind to the source vertices.
    for lat in out[3:4]
        rid=kernel.bindings[1,lat[2]]
        ra_id,rb_id=kernel.shapes[rid].children
        @test ra_id==rb_id
        @test kernel.shapes[rid].geometry.geometry.r≈2.0
    end

    # Radial line rho 1..4 at h=0: annulus face + two rims.
    pa=occ.add_point(1.0,0.0,0.0)
    pb=occ.add_point(4.0,0.0,0.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    @test length(out)==4
    fid=kernel.bindings[2,out[2][2]]
    fg=kernel.shapes[fid].geometry
    @test fg.kind===:plane
    @test length(kernel.shapes[fid].children)==2  # outer + hole wires
    @test occ.get_mass(out[2]...)≈pi*15

    # Slanted line (rho 1->3, h 0->2): cone lateral.
    pa=occ.add_point(1.0,0.0,0.0)
    pb=occ.add_point(3.0,0.0,2.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    fid=kernel.bindings[2,out[2][2]]
    fg=kernel.shapes[fid].geometry
    @test fg.kind===:cone
    @test fg.geometry.r1≈1.0 && fg.geometry.r2≈3.0
    @test occ.get_mass(out[2]...)≈pi*4*sqrt(8.0)

    # Descending child order (a=hi, b=lo): the seam's `lin2d` pcurve
    # maps the edge's own parameter to v: slant->0 on the u=2π leg.
    ld=occ.add_line(pb,pa)
    out=occ.revolve([(1,ld)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    fid=kernel.bindings[2,out[2][2]]
    fg=kernel.shapes[fid].geometry
    @test fg.kind===:cone
    seid=kernel.bindings[1,ld]
    sg=kernel.shapes[seid].geometry.geometry
    pc=fg.geometry.pcurves[ld]
    @test M._occ_pcurve_d0(pc.fwd,sg.t0)[2]≈sqrt(8.0)
    @test M._occ_pcurve_d0(pc.fwd,sg.t1)[2]≈0.0 atol=1e-12
    @test M._occ_pcurve_d0(pc.rev,sg.t0)[2]≈sqrt(8.0)
    @test M._occ_pcurve_d0(pc.rev,sg.t1)[2]≈0.0 atol=1e-12

    # On-axis line and cross-meridian line reject.
    a=occ.add_point(0.0,0.0,0.0)
    b=occ.add_point(0.0,0.0,1.0)
    laxis=occ.add_line(a,b)
    @test_throws ArgumentError occ.revolve(
        [(1,laxis)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    a=occ.add_point(1.0,0.0,0.0)
    b=occ.add_point(0.0,1.0,0.0)
    lcross=occ.add_line(a,b)
    @test_throws ArgumentError occ.revolve(
        [(1,lcross)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
end

@testset "occ.revolve planar profiles to solids" begin
    occ=_occ_setup()
    # Rectangle on the axis -> compact cylinder encoding.
    _,_,f=_meridian_rect(occ,0.0,2.0,0.0,5.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    @test out[1]==(2,f)
    @test out[2][1]==3
    @test length(out)==5  # top, solid, 3 laterals (axis edge skips)
    kernel=_kernel()
    sid=kernel.bindings[3,out[2][2]]
    prim=kernel.shapes[sid].geometry.primitive
    @test prim.occ===:cylinder
    @test prim.record.radius≈2.0
    @test prim.record.height≈5.0
    @test occ.get_mass(out[2]...)≈pi*4*5
    @test collect(occ.get_center_of_mass(out[2]...))≈[0.0,0.0,2.5]
    # Shell: two disk caps + one cylinder, axis edge contributes none.
    shell=kernel.shapes[only(kernel.shapes[sid].children)]
    kinds=sort!([kernel.shapes[fid].geometry.kind
                 for (fid,_) in shell.geometry.uses])
    @test kinds==[:cylinder,:plane,:plane]

    # Triangle on the axis -> compact cone encoding.
    p0=occ.add_point(0.0,0.0,0.0)
    pa=occ.add_point(3.0,0.0,0.0)
    pb=occ.add_point(0.0,0.0,4.0)
    es=[occ.add_line(p0,pa),occ.add_line(pa,pb),occ.add_line(pb,p0)]
    f2=occ.add_plane_surface([occ.add_wire(es)])
    out=occ.revolve([(2,f2)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    sid=kernel.bindings[3,out[2][2]]
    prim=kernel.shapes[sid].geometry.primitive
    @test prim.occ===:cone
    @test prim.record.r1≈3.0 && prim.record.r2≈0.0
    @test occ.get_mass(out[2]...)≈pi*9*4/3

    # Offset rectangle -> revolved annular solid, Pappus encoding.
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    n1=count(e->e[1]==1,occ.get_entities())
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    n2=count(e->e[1]==1,occ.get_entities())
    @test n2-n1==4  # four deduplicated rim circles
    sid=kernel.bindings[3,out[2][2]]
    prim=kernel.shapes[sid].geometry.primitive
    @test prim.occ===:revolved
    @test prim.record.V≈2pi*1.5*1
    @test prim.record.h_c≈0.5
    @test occ.get_mass(out[2]...)≈2pi*1.5*1
    @test collect(occ.get_center_of_mass(out[2]...))≈[0.0,0.0,0.5]
    @test collect(occ.get_bounding_box(3,out[2][2]))≈
        [-2.0,-2.0,0.0,2.0,2.0,1.0]
    shell=kernel.shapes[only(kernel.shapes[sid].children)]
    kinds=sort!([kernel.shapes[fid].geometry.kind
                 for (fid,_) in shell.geometry.uses])
    @test kinds==[:cylinder,:cylinder,:plane,:plane]
    # Every face's `σ·axis` (or implied normal for the quadrics)
    # points away from the material: the bottom cap reads -z at h=0,
    # the top +z, and exactly one of the two cylinders is reversed.
    for (fid,rev) in shell.geometry.uses
        g=kernel.shapes[fid].geometry
        g.kind===:plane || continue
        n=rev ? M._occ3_scale(g.geometry.axis,-1.0) : g.geometry.axis
        if g.geometry.center[3]≈0.0
            @test n[3]<0
        else
            @test n[3]>0
        end
    end
    @test Set([rev for (fid,rev) in shell.geometry.uses
               if kernel.shapes[fid].geometry.kind===:cylinder])==
        Set([true,false])

    # Disk in an axial plane -> compact torus publication.
    d=occ.add_disk(5.0,0.0,0.0,1.0,1.0;zAxis=[0.0,1.0,0.0])
    out=occ.revolve([(2,d)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    @test length(out)==3
    @test out[1]==(2,d)
    @test out[2][1]==3
    @test out[3][1]==2
    @test occ.get_mass(out[2]...)≈2*pi*pi*5
    @test collect(occ.get_center_of_mass(out[2]...))≈[0.0,0.0,0.0] atol=1e-12
end

@testset "occ.revolve partial angles" begin
    occ=_occ_setup()
    kernel=_kernel()
    a=pi/2

    # Point -> rim arc over [0,a]; the rotated copy is the top.
    p=occ.add_point(2.0,0.0,0.0)
    out=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    @test length(out)==2
    @test out[1][1]==0 && out[1][2]!=p
    @test out[2][1]==1
    eid=kernel.bindings[1,out[2][2]]
    eg=kernel.shapes[eid].geometry.geometry
    @test eg.occ===:circle
    @test eg.r≈2.0
    @test eg.t0==0.0 && eg.t1==a
    a_id,b_id=kernel.shapes[eid].children
    @test a_id==kernel.bindings[0,p]
    @test b_id==kernel.bindings[0,out[1][2]]
    @test collect(kernel.shapes[b_id].geometry.point)≈[0.0,2.0,0.0] atol=1e-14

    # Axial line -> cylinder-sector lateral; rim arcs are the laterals.
    pa=occ.add_point(2.0,0.0,0.0)
    pb=occ.add_point(2.0,0.0,3.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    @test length(out)==4
    @test out[1][1]==1 && out[1][2]!=l   # rotated copy
    @test out[2][1]==2
    fid=kernel.bindings[2,out[2][2]]
    fg=kernel.shapes[fid].geometry
    @test fg.kind===:cylinder
    @test fg.geometry.angle≈a
    @test fg.geometry.radius≈2.0
    @test fg.geometry.height≈3.0
    occ.synchronize()
    @test occ.get_mass(out[2]...)≈a*2*3
    for lat in out[3:4]
        rg=kernel.shapes[kernel.bindings[1,lat[2]]].geometry.geometry
        @test rg.occ===:circle && rg.t1≈a
    end

    # Radial line -> annular-sector plane face.
    pa=occ.add_point(1.0,0.0,0.0)
    pb=occ.add_point(4.0,0.0,0.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    fid=kernel.bindings[2,out[2][2]]
    @test kernel.shapes[fid].geometry.kind===:plane
    @test occ.get_mass(out[2]...)≈(a/2)*(4*4-1*1)

    # Meridian rectangle -> solid sector: walls + two meridian caps;
    # the rotated copy answers as top.
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    @test out[1][1]==2 && out[1][2]!=f
    @test out[2][1]==3
    @test length(out)==6
    sid=kernel.bindings[3,out[2][2]]
    prim=kernel.shapes[sid].geometry.primitive
    @test prim.occ===:revolved
    @test prim.record.V≈a*1.5*1
    @test prim.record.angle≈a
    @test occ.get_mass(out[2]...)≈a*1.5*1
    com=occ.get_center_of_mass(out[2]...)
    # Bisector direction at a/2; arm = J2-based radial offset.
    arm=prim.record.arm
    @test collect(com)≈[arm*sqrt(0.5),arm*sqrt(0.5),0.5] atol=1e-12
    shell=kernel.shapes[only(kernel.shapes[sid].children)]
    @test length(shell.geometry.uses)==6  # 4 walls + 2 caps

    # Axis-touching rectangle -> compact cylinder sector record.
    _,_,f=_meridian_rect(occ,0.0,2.0,0.0,5.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    sid=kernel.bindings[3,out[2][2]]
    sec=kernel.shapes[sid].geometry.sector
    @test sec.occ===:cylinder
    @test sec.radius≈2.0 && sec.height≈5.0 && sec.angle≈a
    @test occ.get_mass(out[2]...)≈a*2*2*5/2
    @test collect(occ.get_bounding_box(3,out[2][2]))≈
        [0.0,0.0,0.0,2.0,2.0,5.0]

    # Triangle on the axis -> cone sector.
    p0=occ.add_point(0.0,0.0,0.0)
    pa=occ.add_point(3.0,0.0,0.0)
    pb=occ.add_point(0.0,0.0,4.0)
    es=[occ.add_line(p0,pa),occ.add_line(pa,pb),occ.add_line(pb,p0)]
    f2=occ.add_plane_surface([occ.add_wire(es)])
    out=occ.revolve([(2,f2)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    sid=kernel.bindings[3,out[2][2]]
    sec=kernel.shapes[sid].geometry.sector
    @test sec.occ===:cone
    @test occ.get_mass(out[2]...)≈a*3*4/2*1.0  # a·A·rho_bar
        atol=1e-12

    # Disk -> torus segment; the rotated disk copy is the top.
    d=occ.add_disk(5.0,0.0,0.0,1.0,1.0;zAxis=[0.0,1.0,0.0])
    out=occ.revolve([(2,d)],0.0,0.0,0.0,0.0,0.0,1.0,a)
    @test out[1][1]==2 && out[1][2]!=d
    @test out[2][1]==3 && out[3][1]==2
    @test occ.get_mass(out[2]...)≈(a/(2pi))*2*pi*pi*5
    bb=occ.get_bounding_box(3,out[2][2])
    @test bb[1]>=-1e-12 && bb[2]>=-1e-12 && bb[4]≈6.0

    # Negative partial angle is the flipped-axis positive turn.
    p=occ.add_point(2.0,0.0,0.0)
    o1=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,-1.0,-a)
    e1=kernel.shapes[kernel.bindings[1,o1[2][2]]].geometry.geometry
    @test collect(e1.n)≈[0.0,0.0,1.0]
    @test e1.t1≈a
    # The endpoint lands as a +z/+a rotation would place it.
    b_id=kernel.shapes[kernel.bindings[1,o1[2][2]]].children[2]
    @test collect(kernel.shapes[b_id].geometry.point)≈[0.0,2.0,0.0] atol=1e-14
end

@testset "occ.revolve axes, angles and ordering" begin
    occ=_occ_setup()
    # Negative angle is the flipped-axis full turn.
    p=occ.add_point(2.0,0.0,0.0)
    o1=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    o2=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,-1.0,-2pi)
    kernel=_kernel()
    e1=kernel.shapes[kernel.bindings[1,o1[2][2]]].geometry.geometry
    e2=kernel.shapes[kernel.bindings[1,o2[2][2]]].geometry.geometry
    # (axis -z, angle -2pi) is the same turn as (axis +z, angle +2pi):
    # the swept circles are identical.
    @test collect(e1.center)≈collect(e2.center)
    @test e1.r≈e2.r≈2.0
    @test collect(e1.n)≈collect(e2.n)

    # Translated axis point: line at (2,2,*) about (1,2,*) -> r=1 cyl.
    pa=occ.add_point(2.0,2.0,0.0)
    pb=occ.add_point(2.0,2.0,4.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],1.0,2.0,0.0,0.0,0.0,1.0,2pi)
    occ.synchronize()
    @test occ.get_mass(out[2]...)≈2pi*1*4

    # Tilted axis: an x-directed line at rho=3 about the x-axis.
    pa=occ.add_point(0.0,3.0,0.0)
    pb=occ.add_point(2.0,3.0,0.0)
    l=occ.add_line(pa,pb)
    out=occ.revolve([(1,l)],0.0,0.0,0.0,1.0,0.0,0.0,2pi)
    fid=_kernel().bindings[2,out[2][2]]
    fg=_kernel().shapes[fid].geometry.geometry
    @test fg.occ===:cylinder
    @test collect(fg.axis)≈[1.0,0.0,0.0]

    # Ordering: dims group 0,1,2 regardless of input order; per input
    # [top, body, laterals].
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    p=occ.add_point(3.0,0.0,0.0)
    q1=occ.add_point(4.0,0.0,0.0)
    q2=occ.add_point(4.0,0.0,1.0)
    ln=occ.add_line(q1,q2)
    out=occ.revolve([(2,f),(0,p),(1,ln)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    @test out[1]==(0,p)
    @test out[2][1]==1
    @test out[3]==(1,ln)
    @test out[4][1]==2
    @test out[7]==(2,f)
    @test out[8][1]==3
    @test all(o[1]==2 for o in out[9:end])
end

@testset "occ.revolve rejections are atomic" begin
    occ=_occ_setup()
    function check(f)
        before=_occ_counts()
        @test_throws ArgumentError f()
        @test _occ_counts()==before
    end
    p_on=occ.add_point(0.0,0.0,1.0)
    p_off=occ.add_point(1.0,0.0,0.0)
    check(()->occ.revolve([(0,p_on)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    check(()->occ.revolve([(0,p_off)],0.0,0.0,0.0,0.0,0.0,0.0,2pi))
    check(()->occ.revolve([(0,p_off)],0.0,0.0,0.0,0.0,0.0,1.0,0.0))
    check(()->occ.revolve([(0,p_off)],0.0,0.0,0.0,0.0,0.0,1.0,3pi))
    check(()->occ.revolve([(0,p_off)],0.0,0.0,0.0,0.0,0.0,NaN,2pi))
    check(()->occ.revolve([(0,999)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    check(()->occ.revolve([(4,1)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    check(()->occ.revolve([(0,p_off)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
        numElements=[1]))
    check(()->occ.revolve([(0,p_off)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
        heights=[0.5]))
    c=occ.add_circle(2.0,0.0,0.0,0.5;zAxis=[0.0,1.0,0.0])
    check(()->occ.revolve([(1,c)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    # A face in the xy plane does not contain the z axis.
    pf=[occ.add_point(1.0,0.0,0.0),occ.add_point(2.0,0.0,0.0),
        occ.add_point(2.0,1.0,0.0),occ.add_point(1.0,1.0,0.0)]
    ef=[occ.add_line(pf[1],pf[2]),occ.add_line(pf[2],pf[3]),
        occ.add_line(pf[3],pf[4]),occ.add_line(pf[4],pf[1])]
    fnp=occ.add_plane_surface([occ.add_wire(ef)])
    check(()->occ.revolve([(2,fnp)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    # Hole-bounded face rejects.
    po=[occ.add_point(1.0,0.0,0.0),occ.add_point(2.0,0.0,0.0),
        occ.add_point(2.0,0.0,1.0),occ.add_point(1.0,0.0,1.0)]
    eo=[occ.add_line(po[1],po[2]),occ.add_line(po[2],po[3]),
        occ.add_line(po[3],po[4]),occ.add_line(po[4],po[1])]
    pi_=[occ.add_point(1.2,0.0,0.2),occ.add_point(1.8,0.0,0.2),
         occ.add_point(1.8,0.0,0.8),occ.add_point(1.2,0.0,0.8)]
    ei=[occ.add_line(pi_[1],pi_[2]),occ.add_line(pi_[2],pi_[3]),
        occ.add_line(pi_[3],pi_[4]),occ.add_line(pi_[4],pi_[1])]
    fh2=occ.add_plane_surface([occ.add_wire(eo),occ.add_wire(ei)])
    check(()->occ.revolve([(2,fh2)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    # Under Layers a two-subdivision full turn welds coincident cells
    # onto the source row, and an axis-touching profile collapses the
    # column machinery — both reject before any bind.
    check(()->occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
        numElements=[2]))
    _,_,fax=_meridian_rect(occ,0.0,2.0,0.0,1.0)
    check(()->occ.revolve([(2,fax)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
        numElements=[4]))
    # Spindle disk and axis-parallel disk reject.
    dbad=occ.add_disk(0.5,0.0,0.0,1.0,1.0;zAxis=[0.0,1.0,0.0])
    check(()->occ.revolve([(2,dbad)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    dpar=occ.add_disk(5.0,0.0,0.0,1.0,1.0;zAxis=[0.0,0.0,1.0])
    check(()->occ.revolve([(2,dpar)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    bx=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    check(()->occ.revolve([(3,bx)],0.0,0.0,0.0,0.0,0.0,1.0,2pi))
    # A mixed valid+invalid call still rolls back entirely.
    lcross=occ.add_line(occ.add_point(1.0,0.0,0.0),
                        occ.add_point(0.0,1.0,0.0))
    check(()->occ.revolve([(0,p_off),(1,lcross)],
        0.0,0.0,0.0,0.0,0.0,1.0,2pi))
end

@testset "occ.revolve cylinder solid meshes" begin
    occ=_occ_setup()
    _,_,f=_meridian_rect(occ,0.0,1.5,0.0,2.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    occ.synchronize()
    solid_tag=out[2][2]
    m=_current_model()
    @test haskey(m.volumes,solid_tag)
    @test haskey(m.cylinders,solid_tag)
    A.mesh.generate(3)
    node_tags,coords=A.mesh.get_nodes(3,solid_tag,true)
    @test length(node_tags)>0
    ok=true
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        ok &= sqrt(x*x+y*y)<=1.5+1e-6 &&
              -1e-6<=z<=2.0+1e-6
    end
    @test ok
    element_types,element_tags,_=A.mesh.get_elements(3,solid_tag)
    @test element_types==Int32[4]
end

@testset "occ.revolve Layers transport" begin
    # A point generatrix under Layers sweeps a graded rim circle; two
    # subdivisions tile it as a two-arc ring.
    occ=_occ_setup()
    p=occ.add_point(3.0,0.0,1.0)
    out=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
                    numElements=[2])
    rim_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(1)
    _,coords=A.mesh.get_nodes(1,rim_tag,true)
    @test length(coords)>=2*3
    element_types,_,_=A.mesh.get_elements(1,rim_tag)
    @test element_types==Int32[1]

    # Annular profile -> revolved solid: the swept volume grades at the
    # normalized heights; recombination emits prisms (type 6).
    occ=_occ_setup()
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
                    numElements=[4])
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    _,coords=A.mesh.get_nodes(3,solid_tag,true)
    phis=sort!(unique!(round.(
        [let a=mod(atan(coords[i+1],coords[i]),2pi)
             a>2pi-1e-6 ? 0.0 : a
         end for i in 1:3:length(coords)],digits=9)))
    @test phis≈collect(range(0.0,2pi,length=5)[1:4]) atol=1e-9

    occ=_occ_setup()
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
                    numElements=[4],recombine=true)
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test 6 in element_types
    @test !(4 in element_types)

    # Disk -> torus with graded heights: layers [2,2] at cumulative
    # boundaries 0.4/1.0 grade the turn at 0, 0.2, 0.4, 0.7 of a turn —
    # azimuths 0, 0.4pi, 0.8pi, 1.4pi (the 2pi level welds onto 0).
    occ=_occ_setup()
    d=occ.add_disk(5.0,0.0,0.0,1.0,1.0;zAxis=[0.0,1.0,0.0])
    out=occ.revolve([(2,d)],0.0,0.0,0.0,0.0,0.0,1.0,2pi;
                    numElements=[2,2],heights=[0.4,1.0],
                    recombine=true)
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    _,coords=A.mesh.get_nodes(3,solid_tag,true)
    phis=sort!(unique!(round.(
        [let a=mod(atan(coords[i+1],coords[i]),2pi)
             a>2pi-1e-6 ? 0.0 : a
         end for i in 1:3:length(coords)],digits=9)))
    @test phis≈[0.0,0.4pi,0.8pi,1.4pi] atol=1e-9
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test 6 in element_types

    # Partial turn: the same transport grades the open sector — the
    # closing column lands on the rotated copy, not the source row, so
    # a single subdivision already tiles the strip.
    occ=_occ_setup()
    _,_,f=_meridian_rect(occ,1.0,2.0,0.0,1.0)
    out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,pi/2;
                    numElements=[3],recombine=true)
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    _,coords=A.mesh.get_nodes(3,solid_tag,true)
    phis=sort!(unique!(round.(
        [mod(atan(coords[i+1],coords[i]),2pi)
         for i in 1:3:length(coords)],digits=9)))
    @test phis≈collect(range(0.0,pi/2,length=4)) atol=1e-9
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test 6 in element_types
    @test !(4 in element_types)
end
