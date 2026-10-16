# `model.occ.addThruSections` — `BRepOffsetAPI_ThruSections` lofting over
# consecutive wire sections. Smooth mode interpolates section centroids
# along a spine (approximating cones/frustums); `makeRuled` interpolates
# degree-1 between matching vertices (exact for prisms/pyramid frustums).
# `makeSolid` adds planar caps and binds shell+solid; surface mode binds
# one face per lateral. All failures are atomic.

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

# Axis-aligned circular wire at height z with radius r.
_circle_wire(occ,r,z)=occ.add_wire([occ.add_circle(0.0,0.0,z,r)])

# Axis-aligned square wire of half-side h at height z.
function _square_wire(occ,h,z)
    ps=[occ.add_point(h*q[1],h*q[2],z)
        for q in [(1,1),(-1,1),(-1,-1),(1,-1)]]
    occ.add_wire([occ.add_line(ps[i],ps[mod1(i+1,4)]) for i in 1:4])
end

@testset "occ.add_thru_sections" begin
    occ=_occ_setup()
    try
        # ── two-section circular frustum, makeSolid ──────────────────
        w1=_circle_wire(occ,1.0,0.0)
        w2=_circle_wire(occ,0.5,2.0)
        out=occ.add_thru_sections([w1,w2])
        @test out==[(3,out[1][2])]
        @test out[1][1]==3
        solid=out[1][2]
        # Analytic frustum: V = πh/3 (R²+Rr+r²) = π·2/3·1.75.
        @test occ.get_mass(3,solid)≈3.665191429188092 rtol=1e-6
        cx,cy,cz=occ.get_center_of_mass(3,solid)
        @test cx≈0.0 atol=1e-6
        @test cy≈0.0 atol=1e-6
        # Frustum centroid measured from the wider base:
        #   z̄ = h(R²+2Rr+3r²)/(4(R²+Rr+r²)) = 2·2.75/7.
        @test cz≈2*2.75/7 rtol=1e-4
        bb=occ.get_bounding_box(3,solid)
        @test bb[1]≈-1.0 atol=1e-6   # Confusion-margin padded like upstream
        @test bb[4]≈1.0 atol=1e-6
        @test bb[3]≈0.0 atol=1e-6
        @test bb[6]≈2.0 atol=1e-6
        # Caps + one lateral per section pair = 3 surfaces.
        @test length(occ.get_entities(2))==3

        # ── surface mode ─────────────────────────────────────────────
        w3=_circle_wire(occ,1.0,4.0)
        w4=_circle_wire(occ,0.5,6.0)
        out=occ.add_thru_sections([w3,w4];makeSolid=false)
        @test all(d==2 for (d,_) in out)
        # Lateral area of a circular frustum: π(R+r)·slant,
        # slant=√(h²+(R−r)²)=√4.25; smooth spine approximates the cone.
        @test occ.get_mass(2,out[1][2])≈pi*1.5*sqrt(4.25) rtol=1e-3

        # ── three-section loft through a waist ───────────────────────
        w5=_circle_wire(occ,1.0,8.0)
        w6=_circle_wire(occ,0.6,9.0)
        w7=_circle_wire(occ,0.9,10.0)
        out=occ.add_thru_sections([w5,w6,w7])
        @test out[1][1]==3
        m=occ.get_mass(3,out[1][2])
        @test 1.0<m<8.0   # smooth waist — between cylinders of the extremes

        # ── ruled square frustum (exact planar laterals) ─────────────
        s1=_square_wire(occ,1.0,12.0)
        s2=_square_wire(occ,0.5,14.0)
        out=occ.add_thru_sections([s1,s2];makeRuled=true)
        # Pyramid frustum: V = h/3·(A1+√A1A2+A2) = 2/3·(4+2+1) = 14/3.
        @test occ.get_mass(3,out[1][2])≈14/3 rtol=1e-9
        cx,cy,cz=occ.get_center_of_mass(3,out[1][2])
        @test cx≈0.0 atol=1e-9
        @test cy≈0.0 atol=1e-9
        @test 12.0<cz<14.0

        # ── ruled square→near-apex pyramid ───────────────────────────
        s3=_square_wire(occ,1.0,16.0)
        apex=_square_wire(occ,1e-3,18.0)
        out=occ.add_thru_sections([s3,apex];makeRuled=true)
        m=occ.get_mass(3,out[1][2])
        @test m>2.6 && m<2.68   # → 8/3 as apex half-side → 0

        # ── atomic rejection ─────────────────────────────────────────
        before=length(occ.get_entities(-1))
        @test occursin("at least 2 wires",
                       _errtext(()->occ.add_thru_sections([w1])))
        @test occursin("Unknown OpenCASCADE",
                       _errtext(()->occ.add_thru_sections([w1,99999])))
        @test occursin("smoothing",
                       _errtext(()->occ.add_thru_sections([w1,w2];
                                                         smoothing=true)))
        @test occursin("continuity",
                       _errtext(()->occ.add_thru_sections([w1,w2];
                                                         continuity="bogus")))
        @test occursin("parametrization",
                       _errtext(()->occ.add_thru_sections([w1,w2];
                                                         parametrization="bogus")))
        @test length(occ.get_entities(-1))==before
    finally
        A.finalize()
    end
end
