# `model.occ.fillet` — BRepFilletAPI_MakeFillet for planar-faced
# solids: each filleted edge becomes an exact cylindrical sector
# tangent to both adjacent faces; end vertices become circular arcs
# on the (orthogonal) corner faces.
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

function _edges_in_box(occ,x0)
    [t for (d,t) in occ.get_entities_in_bounding_box(
         x0-1e-9,-1e-9,-1e-9,x0+1+1e-9,1+1e-9,1+1e-9;dim=1)]
end

@testset "occ.fillet" begin
    A.initialize(); occ=A.model.occ
    try
        # Single-edge fillet on a unit cube: removes the prism
        # between the edge and the tangent lines minus the cylinder
        # sector — V = 1 - r²(1-π/4)·L.
        r=0.25
        v1=occ.add_box(0,0,0,1,1,1)
        e1=first(_edges_in_box(occ,0))
        out=occ.fillet([v1],[e1],r)
        @test out==[(3,v1)]
        @test occ.get_mass(3,v1)≈1.0-r^2*(1-π/4) atol=1e-8

        # Four parallel edges — disjoint fillets.
        v2=occ.add_box(2,0,0,1,1,1)
        ves=[t for t in _edges_in_box(occ,2)
             if (b=occ.get_bounding_box(1,t);
                 b[6]-b[3]>0.5 && b[4]-b[1]<1e-6 && b[5]-b[2]<1e-6)]
        @test length(ves)==4
        occ.fillet([v2],ves,fill(0.15,4))
        @test occ.get_mass(3,v2)≈1.0-4*0.15^2*(1-π/4) atol=1e-8

        # Filleted edges sharing a vertex reject (corner blend).
        v3=occ.add_box(4,0,0,1,1,1)
        es3=_edges_in_box(occ,4)
        # two edges incident at the same corner: pick any pair that
        # shares an endpoint — corner (4,0,0): the three edges whose
        # bbox contains that vertex
        corner=[t for t in es3
                if (b=occ.get_bounding_box(1,t);
                    abs(b[1]-4)<1e-6 && abs(b[2])<1e-6 &&
                    (abs(b[3])<1e-6 || abs(b[6])<1e-6))]
        @test length(corner)>=2
        @test occursin("share a vertex",
                       _errtext(()->occ.fillet([v3],corner[1:2],0.1)))

        # Rejections.
        @test occursin("unknown OCC Volume",
                       _errtext(()->occ.fillet([99999],[1],0.1)))
        v4=occ.add_box(6,0,0,1,1,1)
        @test occursin("unknown OCC Curve",
                       _errtext(()->occ.fillet([v4],[99999],0.1)))
        @test occursin("radii must be positive",
                       _errtext(()->occ.fillet([v2],[1],0.0)))
        e_v1=first(_edges_in_box(occ,0))
        @test occursin("not an edge",
                       _errtext(()->occ.fillet([v3],[e_v1],0.1)))
    finally
        A.finalize()
    end
end
