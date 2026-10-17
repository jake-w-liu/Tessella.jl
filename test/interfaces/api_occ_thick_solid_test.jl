# `model.occ.addThickSolid` — `BRepOffsetAPI_MakeThickSolid` on planar
# solids: the offset body is the half-space intersection of shifted
# kept-face planes + unshifted excluded-face planes; the shell wall is a
# boolean difference. Curved solids reject loudly.
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

@testset "occ.add_thick_solid" begin
    occ=_occ_setup()
    try
        # ── inward shell of a box, all faces kept → closed hollow box ─
        # Box [0,2]³, offset −0.2: inner body [0.2,1.8]³, wall volume
        # = 8 − 1.6³.
        b=occ.add_box(0.,0.,0.,2.,2.,2.)
        t=occ.add_thick_solid(b,Int[],-0.2)
        @test t>0
        @test occ.get_mass(3,t)≈8.0-1.6^3 atol=1e-6
        xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(3,t)
        @test xmin≈0.0 atol=1e-6
        @test xmax≈2.0 atol=1e-6

        # ── open-top box: exclude the +z face ───────────────────────
        gmsh_faces=[t for (d,t) in occ.get_entities(2)]
        # identify the top face: its centroid z == 2
        top=-1
        for f in gmsh_faces
            _,_,z=occ.get_center_of_mass(2,f)
            z≈2.0 && (top=f)
        end
        @test top>0
        t2=occ.add_thick_solid(b,[top],-0.2)
        @test t2>0
        # inner body: sides shrink to [0.2,1.8]² but the open top keeps
        # the cavity running to z=2 → inner volume = 1.6·1.6·1.8 = 4.608
        @test occ.get_mass(3,t2)≈8.0-1.6*1.6*1.8 atol=1e-6

        # ── outward offset: wall grows outside the original ─────────
        t3=occ.add_thick_solid(b,Int[],0.2)
        @test t3>0
        # dilated body [−0.2,2.2]³ minus original [0,2]³
        @test occ.get_mass(3,t3)≈2.4^3-8.0 atol=1e-6

        # ── non-box planar solid: triangular prism wedge ────────────
        # Wedge: right triangle base (legs 3 along x, 2 along y),
        # height 1. Faces: bottom/top triangles + 3 side quads.
        p1=occ.add_point(0.,0.,0.); p2=occ.add_point(3.,0.,0.)
        p3=occ.add_point(0.,2.,0.)
        p4=occ.add_point(0.,0.,1.); p5=occ.add_point(3.,0.,1.)
        p6=occ.add_point(0.,2.,1.)
        l1=occ.add_line(p1,p2); l2=occ.add_line(p2,p3); l3=occ.add_line(p3,p1)
        l4=occ.add_line(p4,p5); l5=occ.add_line(p5,p6); l6=occ.add_line(p6,p4)
        l7=occ.add_line(p1,p4); l8=occ.add_line(p2,p5); l9=occ.add_line(p3,p6)
        # Wires wound so every face normal is outward (negative first
        # tag binds w.Reversed()).
        wb=occ.add_wire([-l2,l3,l1])      # bottom: −z
        wt=occ.add_wire([l4,l5,l6])       # top: +z
        ws1=occ.add_wire([l1,l8,l4,l7])   # y=0 face: −y
        ws2=occ.add_wire([l2,l9,l5,l8])   # hypotenuse: +(2,3,0)/√13
        ws3=occ.add_wire([-l9,l6,l7,l3])  # x=0 face: −x
        fb=occ.add_plane_surface([wb]); ft=occ.add_plane_surface([wt])
        fs1=occ.add_plane_surface([ws1]); fs2=occ.add_plane_surface([ws2])
        fs3=occ.add_plane_surface([ws3])
        sl=occ.add_surface_loop([fb,ft,fs1,fs2,fs3])
        wv=occ.add_volume([sl])
        @test occ.get_mass(3,wv)≈3.0 atol=1e-6     # V = ½·3·2·1
        wt2=occ.add_thick_solid(wv,Int[],-0.2)
        @test wt2>0
        inner_mass=occ.get_mass(3,wv)-occ.get_mass(3,wt2)
        # inner body = same-shape prism shrunk by 0.2 on every face
        @test inner_mass>0.0 && inner_mass<3.0

        # ── curved solid rejects loudly ─────────────────────────────
        cyl=occ.add_cylinder(0.,0.,0.,0.,0.,2.,1.)
        @test occursin("planar",_errtext(()->occ.add_thick_solid(cyl,Int[],-0.1)))

        # ── atomicity ───────────────────────────────────────────────
        nvol=length(occ.get_entities(3))
        @test occursin("unknown",_errtext(()->occ.add_thick_solid(99999,Int[],-0.1)))
        @test_throws ArgumentError occ.add_thick_solid(b,Int[],0.0)
        @test occursin("not a face",_errtext(()->
            occ.add_thick_solid(b,[99999],-0.1)))
        @test length(occ.get_entities(3))==nvol
    finally
        A.finalize()
    end
end
