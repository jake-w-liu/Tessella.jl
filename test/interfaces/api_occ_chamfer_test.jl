# `model.occ.chamfer` — BRepFilletAPI_MakeChamfer for planar-faced
# solids: each chamfered edge adds one half-space plane through the
# lines at distances d1/d2 inside its adjacent faces; the chamfered
# body is the (exact) extended half-space intersection.
using Test
using Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

# Edges of the box rooted at (x0,0,0) — bounding-box query.
function _edges_in_box(occ,x0)
    [t for (d,t) in occ.get_entities_in_bounding_box(
         x0-1e-9,-1e-9,-1e-9,x0+1+1e-9,1+1e-9,1+1e-9;dim=1)]
end

@testset "occ.chamfer" begin
    A.initialize(); occ=A.model.occ
    try
        # Single symmetric chamfer on a unit-cube edge: removes a
        # triangular prism of legs d — V = 1 - d²/2.
        v1=occ.add_box(0,0,0,1,1,1)
        e1=first(_edges_in_box(occ,0))
        out=occ.chamfer([v1],[e1],0.25,0.25)
        @test out==[(3,v1)]
        @test occ.get_mass(3,v1)≈1.0-0.25^2/2 atol=1e-9
        @test length(occ.get_entities(2))==7   # 6 faces + chamfer
        bb=occ.get_bounding_box(3,v1)
        @test all(collect(bb).≈(0.0,0.0,0.0,1.0,1.0,1.0))

        # Asymmetric distances: legs d1,d2 → prism area d1*d2/2.
        v2=occ.add_box(2,0,0,1,1,1)
        e2=first(_edges_in_box(occ,2))
        occ.chamfer([v2],[e2],0.2,0.4)
        @test occ.get_mass(3,v2)≈1.0-0.2*0.4/2 atol=1e-9

        # The four vertical edges of a third cube — disjoint prisms,
        # no corner interaction: V = 1 - 4·d²/2.
        v3=occ.add_box(4,0,0,1,1,1)
        verts3=[t for t in _edges_in_box(occ,4)
                if (b=occ.get_bounding_box(1,t);
                    b[6]-b[3]>0.5 && b[4]-b[1]<1e-6 &&
                    b[5]-b[2]<1e-6)]
        @test length(verts3)==4
        occ.chamfer([v3],verts3,0.2,0.2)
        @test occ.get_mass(3,v3)≈1.0-4*0.2^2/2 atol=1e-9

        # Overrun rejection: a leg longer than the face reaches the
        # opposite corner.
        v4=occ.add_box(6,0,0,1,1,1)
        e4=first(_edges_in_box(occ,6))
        @test occursin("overruns",
                       _errtext(()->occ.chamfer([v4],[e4],1.1,1.1)))
        # Failed chamfer leaves the volume intact (atomicity).
        @test occ.get_mass(3,v4)≈1.0 atol=1e-9

        # Unknown entities and edges of other volumes reject.
        @test occursin("unknown OCC Volume",
                       _errtext(()->occ.chamfer([99999],[1],0.1,0.1)))
        @test occursin("unknown OCC Curve",
                       _errtext(()->occ.chamfer([v1],[99999],0.1,0.1)))
        e_v1=first(_edges_in_box(occ,0))
        @test occursin("not an edge",
                       _errtext(()->occ.chamfer([v4],[e_v1],0.1,0.1)))
    finally
        A.finalize()
    end
end
