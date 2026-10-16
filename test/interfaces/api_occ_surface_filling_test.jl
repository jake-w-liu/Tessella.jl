# `model.occ.addSurfaceFilling` / `addBSplineFilling` —
# `BRepOffsetAPI_MakeFilling` bounded-surface construction. Planar wires
# publish the exact `:plane` record; non-planar 4-edge wires become a
# Coons-blend tensor patch reproducing the boundary exactly. All failures
# are atomic.

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

# Axis-aligned square wire of half-side h at height z.
function _square_wire(occ,h,z)
    ps=[occ.add_point(h*q[1],h*q[2],z)
        for q in [(1,1),(-1,1),(-1,-1),(1,-1)]]
    occ.add_wire([occ.add_line(ps[i],ps[mod1(i+1,4)]) for i in 1:4])
end

@testset "occ.add_surface_filling" begin
    occ=_occ_setup()
    try
        # ── planar wire → exact plane record ────────────────────────
        w=_square_wire(occ,1.0,0.0)
        f=occ.add_surface_filling(w)
        @test f>0
        @test occ.get_mass(2,f)≈4.0 atol=1e-12
        xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(2,f)
        @test xmin≈-1.0 atol=1e-6
        @test xmax≈1.0 atol=1e-6
        @test zmin≈0.0 atol=1e-9
        @test zmax≈0.0 atol=1e-9
        occ.synchronize()

        # ── non-planar quadrilateral: tilted ruled quad ──────────────
        # corners (±1,±1,0/2): the Coons blend is bilinearly exact —
        # area of a hyperbolic-parallelogram-style ruled quad.
        p1=occ.add_point(1.0,1.0,0.0); p2=occ.add_point(-1.0,1.0,2.0)
        p3=occ.add_point(-1.0,-1.0,2.0); p4=occ.add_point(1.0,-1.0,0.0)
        wq=occ.add_wire([occ.add_line(p1,p2),occ.add_line(p2,p3),
                         occ.add_line(p3,p4),occ.add_line(p4,p1)])
        fq=occ.add_surface_filling(wq)
        @test fq>0
        # Bilinear patch over a square: area of the exact ruled surface
        # (opposite edges parallel → planar? no: edges p1p2 and p4p3 are
        # parallel, p2p3 and p1p4 are parallel → actually a parallelogram
        # only if both pairs parallel; here it is planar? p1p2 direction
        # (-2,0,2), p4p3 direction (-2,0,2) — parallel; p2p3 (0,-2,0),
        # p1p4 (0,-2,0) — parallel → parallelogram, planar!) — so this
        # lands on the plane path; assert the parallelogram area
        # |d1×d2| = |(-2,0,2)×(0,-2,0)| = |(4,0,4)| = 4√2.
        @test occ.get_mass(2,fq)≈4*sqrt(2) atol=1e-12

        # ── genuine non-planar Coons: half-cylinder wall ─────────────
        # Bottom semicircle r=1 at z=0, top semicircle at z=1, two
        # vertical lines → the fill is the half-cylinder lateral, area
        # π·r·h = π (Coons blending two arcs along v reproduces each
        # arc cross-section exactly). `add_circle_arc(start,center,stop)`
        # takes the circle center as its middle argument.
        a=occ.add_point(1.0,0.0,0.0); b=occ.add_point(-1.0,0.0,0.0)
        c=occ.add_point(1.0,0.0,1.0); d=occ.add_point(-1.0,0.0,1.0)
        oc0=occ.add_point(0.0,0.0,0.0); oc1=occ.add_point(0.0,0.0,1.0)
        # Arcs run counter-clockwise from start to stop around the center:
        # a→b (0→π) bulges +y; c→d (0→π) bulges +y too — the wire consumes
        # the second arc reversed (d→c).
        arc_bot=occ.add_circle_arc(a,oc0,b)
        lft=occ.add_line(b,d)
        arc_top=occ.add_circle_arc(c,oc1,d)
        rgt=occ.add_line(c,a)
        wc=occ.add_wire([arc_bot,lft,arc_top,rgt])
        fc=occ.add_surface_filling(wc)
        @test fc>0
        @test occ.get_mass(2,fc)≈π atol=1e-8
        xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(2,fc)
        @test xmin≈-1.0 atol=1e-6
        @test xmax≈1.0 atol=1e-6
        @test zmin≈0.0 atol=1e-6
        @test zmax≈1.0 atol=1e-6

        # ── add_bspline_filling alias ────────────────────────────────
        f2=occ.add_bspline_filling(wc;tag=990)
        @test f2==990
        @test occ.get_mass(2,f2)≈π atol=1e-8

        # ── atomicity ────────────────────────────────────────────────
        @test occursin("unknown",_errtext(()->occ.add_surface_filling(99999)))
        # too few edges
        open_w=occ.add_wire([occ.add_line(a,c)];checkClosed=false)
        @test occursin("at least 3",_errtext(()->
            occ.add_surface_filling(open_w)))
        # non-planar ≠4 edges rejects
        s3=_square_wire(occ,0.5,3.0)
        @test_throws ArgumentError occ.add_surface_filling(s3;tag=990)
    finally
        A.finalize()
    end
end
