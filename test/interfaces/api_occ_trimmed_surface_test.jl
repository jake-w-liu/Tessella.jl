# `model.occ.addTrimmedSurface` — `BRepBuilderAPI_MakeFace(S,W)` bounded
# faces on a parametric base. Edges project to `lin2d`/`spl2d` pcurves;
# no surfaceTag falls back to the plane path. Failures are atomic.
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

# Face of `vol` with the largest mass — picks the barrel over the caps.
function _biggest_face(occ)
    best=(0.0,-1)
    for (d,t) in occ.get_entities(2)
        m=occ.get_mass(2,t)
        m>best[1] && (best=(m,t))
    end
    return best[2]
end

@testset "occ.add_trimmed_surface" begin
    occ=_occ_setup()
    try
        # ── plane fallback: no surfaceTag on a planar wire ──────────
        ps=[occ.add_point(q[1],q[2],0.0)
            for q in [(2.,0.),(2.,3.),(0.,3.),(0.,0.)]]
        w=occ.add_wire([occ.add_line(ps[i],ps[mod1(i+1,4)]) for i in 1:4])
        f=occ.add_trimmed_surface(w)
        @test f>0
        @test occ.get_mass(2,f)≈6.0 atol=1e-9

        # ── explicit plane base with the same wire ──────────────────
        base=occ.get_entities(2)   # the face just made — any plane works
        f2=occ.add_trimmed_surface(w,first(t for (d,t) in base))
        @test f2>0
        @test occ.get_mass(2,f2)≈6.0 atol=1e-9

        # ── cylinder barrel trim: half-wrap strip ───────────────────
        # Barrel r=1, z∈[0,4]. Wire = two semicircle arcs on the barrel
        # (θ∈[0,π]) at z=2 and z=3 plus two vertical lines: a closed
        # on-surface loop bounding the θ∈[0,π], z∈[2,3] strip.
        occ.add_cylinder(0.,0.,0.,0.,0.,4.,1.)
        barrel=_biggest_face(occ)
        a=occ.add_point(1.,0.,2.); b=occ.add_point(-1.,0.,2.)
        c=occ.add_point(1.,0.,3.); d=occ.add_point(-1.,0.,3.)
        oc2=occ.add_point(0.,0.,2.); oc3=occ.add_point(0.,0.,3.)
        # The 180° arcs normalize to n=−z so both bulge −y; the wire
        # consumes the second arc reversed (d→c) for a coherent loop.
        arc_lo=occ.add_circle_arc(a,oc2,b)
        lft=occ.add_line(b,d)
        arc_hi=occ.add_circle_arc(c,oc3,d)
        rgt=occ.add_line(c,a)
        wt=occ.add_wire([arc_lo,lft,arc_hi,rgt])
        ft=occ.add_trimmed_surface(wt,barrel)
        @test ft>0
        # strip area = (arc length π·r)·(height 1) = π
        @test occ.get_mass(2,ft)≈π atol=1e-6
        # analytic centroid of the half-cylinder strip: `add_circle_arc`
        # normalizes the 180° a→b arc about n=−z, so the strip sits on
        # the −y side: centroid (0, −2r/π, mid-z).
        cx,cy,cz=occ.get_center_of_mass(2,ft)
        @test cx≈0.0 atol=1e-8
        @test cy≈-2/π atol=1e-6
        @test cz≈2.5 atol=1e-8
        xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(2,ft)
        @test xmin≈-1.0 atol=1e-4
        @test xmax≈1.0 atol=1e-4
        @test zmin≈2.0 atol=1e-4
        @test zmax≈3.0 atol=1e-4

        # ── chord edges project to spl2d (non-linear uv traces) ─────
        # Three lines between on-barrel points at different heights: the
        # straight 3-D chords project as curved uv splines bounding a
        # genuine curvilinear patch.
        q1=occ.add_point(1.,0.,1.); q2=occ.add_point(0.,1.,2.)
        q3=occ.add_point(-0.5,-sqrt(3)/2,1.5)
        wc=occ.add_wire([occ.add_line(q1,q2),occ.add_line(q2,q3),
                         occ.add_line(q3,q1)])
        fc=occ.add_trimmed_surface(wc,barrel)
        @test fc>0
        mc=occ.get_mass(2,fc)
        @test mc>0.1
        @test mc<8π   # strictly inside the barrel
        xc,yc,zc=occ.get_center_of_mass(2,fc)
        @test 1.0<zc<2.0                  # centroid within the height span
        @test xc*xc+yc*yc<1.0             # centroid inside the bore

        # ── atomicity ───────────────────────────────────────────────
        nfaces=length(occ.get_entities(2))
        @test occursin("unknown",_errtext(()->occ.add_trimmed_surface(99999)))
        @test occursin("unknown",_errtext(()->occ.add_trimmed_surface(wt,99999)))
        @test length(occ.get_entities(2))==nfaces
    finally
        A.finalize()
    end
end
