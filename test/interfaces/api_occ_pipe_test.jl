# `model.occ.addPipe` — `BRepOffsetAPI_MakePipe` semantics via sampled
# parallel-transport section lofting: the profile transports rigidly
# along spine stations (conic densify at π/16, free-form 16 spans) and
# the transported wires feed `addThruSections`. A straight spine with a
# circle profile yields an exact rational cylinder; an arc spine
# approximates the torus segment.

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

@testset "occ.add_pipe" begin
    occ=_occ_setup()
    try
        # ── straight spine, circle profile → exact cylinder ──────────
        # Circle r=0.5 in the xy-plane at the spine start (0,0,0);
        # spine = line z:0→1. Two-station loft = exact rational
        # cylinder, lateral area 2πr·h = π.
        wc=occ.add_wire([occ.add_circle(0.0,0.0,0.0,0.5)])
        pa=occ.add_point(0.0,0.0,0.0); pb=occ.add_point(0.0,0.0,1.0)
        ws=occ.add_wire([occ.add_line(pa,pb)])
        out=occ.add_pipe(wc,ws)
        @test !isempty(out)
        f=first(dt for dt in out if dt[1]==2)[2]
        @test occ.get_mass(2,f)≈π atol=1e-8
        xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(2,f)
        @test xmin≈-0.5 atol=1e-6
        @test xmax≈0.5 atol=1e-6
        @test zmax≈1.0 atol=1e-6

        # ── face profile → solid pipe (exact cylinder volume) ────────
        dk=occ.add_disk(0.0,0.0,2.0,0.5,0.5)
        wsp=occ.add_wire([occ.add_line(occ.add_point(0.0,0.0,2.0),
                                       occ.add_point(0.0,0.0,3.0))])
        v=occ.add_pipe((2,dk),wsp)
        vt=v isa Integer ? v : first(dt for dt in v if dt[1]==3)[2]
        @test vt>0
        @test occ.get_mass(3,vt)≈π*0.25 atol=1e-8
        c=occ.get_center_of_mass(3,vt)
        @test c[1]≈0.0 atol=1e-9
        @test c[2]≈0.0 atol=1e-9
        @test c[3]≈2.5 atol=1e-9

        # ── arc spine → quarter-torus tube ───────────────────────────
        # Spine: 90° circle arc r=1 in the xy plane from (1,0,0) to
        # (0,1,0); profile circle r=0.3 centered at (1,0,0) in the
        # xz-plane (normal +y, matching the start tangent).
        sc=occ.add_point(0.0,0.0,0.0)
        sp1=occ.add_point(1.0,0.0,0.0); sp2=occ.add_point(0.0,1.0,0.0)
        wsa=occ.add_wire([occ.add_circle_arc(sp1,sc,sp2)])
        prc=occ.add_wire([occ.add_circle(1.0,0.0,0.0,0.3;
                                         zAxis=[0.0,1.0,0.0])])
        tout=occ.add_pipe(prc,wsa)
        @test !isempty(tout)
        # Quarter-torus lateral area: (2πr)(arc len) = 2π·0.3·(π/2) =
        # 0.3π² ≈ 2.96; the lofted patch is a dense-section
        # approximation (within a few percent).
        At=sum(dt[1]==2 ? occ.get_mass(2,dt[2]) : 0.0 for dt in tout)
        @test At≈0.3*π*π rtol=3e-2
        occ.synchronize()

        # ── atomicity ────────────────────────────────────────────────
        @test occursin("unknown",_errtext(()->occ.add_pipe(99999,ws)))
        @test occursin("unknown",_errtext(()->occ.add_pipe(wc,99999)))
        @test occursin("unknown",_errtext(()->occ.add_pipe((2,99999),ws)))
    finally
        A.finalize()
    end
end
