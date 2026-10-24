# `model.occ.addThickSections` — Tessella extension (gmsh 4.15.2 has no
# endpoint of this name; verified against the installed Python module,
# the generated api/gen.py, GModelIO_OCC.h and the .geo grammar). Each
# listed offset bounds one solid slab along the face normal; adjacent
# slabs share the intermediate section face/wires/edges/vertices —
# the OCCT `BRepOffset` thick-sections recipe built on the certified
# prism/disk sweep. Bounded to the sweep's exact profile classes:
# all-line boundaries and single full-circle profiles.
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

@testset "occ.add_thick_sections" begin
    occ=_occ_setup()
    try
        # ── single slab off a 3×2 rectangle at z=0 ──────────────────
        f=occ.add_rectangle(0.,0.,0.,3.,2.)
        out=occ.add_thick_sections(f,[1.5])
        @test length(out)==1 && out[1][1]==3
        v=out[1][2]
        @test occ.get_mass(3,v)≈9.0 atol=1e-6          # 3·2·1.5
        xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(3,v)
        @test collect((xmin,ymin,zmin))≈[0.0,0.0,0.0] atol=1e-6
        @test collect((xmax,ymax,zmax))≈[3.0,2.0,1.5] atol=1e-6
        @test collect(occ.get_center_of_mass(3,v))≈[1.5,1.0,0.75] atol=1e-6

        # ── two offsets: stacked slabs sharing the middle section ────
        nfaces_before=length(occ.get_entities(2))
        out2=occ.add_thick_sections(f,[0.5,1.5])
        @test length(out2)==2 && all(d==3 for (d,_) in out2)
        va,vb=out2[1][2],out2[2][2]
        @test va!=vb
        @test occ.get_mass(3,va)≈3.0 atol=1e-6         # 3·2·0.5
        @test occ.get_mass(3,vb)≈6.0 atol=1e-6         # 3·2·1.0
        _,_,zmin_a,_,_,zmax_a=occ.get_bounding_box(3,va)
        _,_,zmin_b,_,_,zmax_b=occ.get_bounding_box(3,vb)
        @test zmin_a≈0.0 atol=1e-6
        @test zmax_a≈0.5 atol=1e-6
        @test zmin_b≈0.5 atol=1e-6
        @test zmax_b≈1.5 atol=1e-6
        @test collect(occ.get_center_of_mass(3,vb))≈[1.5,1.0,1.0] atol=1e-6
        # Each shell has 6 faces and the middle section is ONE shared
        # entity at z=0.5 — so only 10 NEW faces bind (2 outer caps +
        # 8 laterals; the source face `f` was already bound).
        loops_a,faces_a=occ.get_surface_loops(va)
        loops_b,faces_b=occ.get_surface_loops(vb)
        fa=only(faces_a); fb=only(faces_b)
        @test length(fa)==6 && length(fb)==6
        shared=intersect(fa,fb)
        @test length(shared)==1
        _,_,zc=occ.get_center_of_mass(2,only(shared))
        @test zc≈0.5 atol=1e-6
        @test length(occ.get_entities(2))==nfaces_before+10

        # ── negative offsets march along −n̂ ─────────────────────────
        f3=occ.add_rectangle(0.,10.,0.,1.,1.)
        vn=only(occ.add_thick_sections(f3,[-1.0]))[2]
        @test occ.get_mass(3,vn)≈1.0 atol=1e-6
        _,_,zmin_n,_,_,zmax_n=occ.get_bounding_box(3,vn)
        @test zmin_n≈-1.0 atol=1e-6
        @test zmax_n≈0.0 atol=1e-6

        # ── negative surface tag flips the direction ────────────────
        f4=occ.add_rectangle(0.,20.,0.,1.,1.)
        vr=only(occ.add_thick_sections(-f4,[2.0]))[2]
        _,_,zmin_r,_,_,zmax_r=occ.get_bounding_box(3,vr)
        @test zmin_r≈-2.0 atol=1e-6
        @test zmax_r≈0.0 atol=1e-6

        # ── disk profile → compact cylinder ──────────────────────────
        df=occ.add_disk(0.,0.,30.,1.,1.)
        vd=only(occ.add_thick_sections(df,[2.0]))[2]
        @test occ.get_mass(3,vd)≈2π atol=1e-5
        _,_,zmin_d,_,_,zmax_d=occ.get_bounding_box(3,vd)
        @test zmin_d≈30.0 atol=1e-6
        @test zmax_d≈32.0 atol=1e-6

        # ── tag binds the sole result (upstream `_multiBind` rule) ──
        f5=occ.add_rectangle(0.,0.,50.,1.,1.)
        @test occ.add_thick_sections(f5,[0.25];tag=77)==[(3,77)]
        @test occ.get_mass(3,77)≈0.25 atol=1e-6

        # ── error classification, atomic ─────────────────────────────
        nvol=length(occ.get_entities(3))
        nface=length(occ.get_entities(2))
        @test occursin("unknown OCC Surface",_errtext(()->
            occ.add_thick_sections(99999,[1.0])))
        @test occursin("at least one offset",_errtext(()->
            occ.add_thick_sections(f,Float64[])))
        @test occursin("non-zero",_errtext(()->
            occ.add_thick_sections(f,[0.0])))
        @test occursin("strictly monotone",_errtext(()->
            occ.add_thick_sections(f,[1.0,0.5])))
        @test occursin("strictly monotone",_errtext(()->
            occ.add_thick_sections(f,[1.0,1.0])))
        @test occursin("cannot bind",_errtext(()->
            occ.add_thick_sections(f,[1.0,2.0];tag=88)))
        @test_throws ArgumentError occ.add_thick_sections(f,[NaN])
        # non-planar surface (sphere face) rejects loudly
        s=occ.add_sphere(0.,0.,100.,1.)
        sph_face=first(t for (d,t) in occ.get_entities(2)
                       if abs(last(occ.get_center_of_mass(2,t))-100)<0.9)
        @test occursin("planar",_errtext(()->
            occ.add_thick_sections(sph_face,[0.1])))
        # elliptical disk profile is outside the certified classes
        ef=occ.add_disk(0.,0.,200.,2.,1.)
        @test occursin("circular",_errtext(()->
            occ.add_thick_sections(ef,[1.0])))
        # failures bind nothing — only the sphere (1 volume, 1 face)
        # and the ellipse disk (1 face) were legitimately added
        @test length(occ.get_entities(3))==nvol+1
        @test length(occ.get_entities(2))==nface+2
    finally
        A.finalize()
    end
end
