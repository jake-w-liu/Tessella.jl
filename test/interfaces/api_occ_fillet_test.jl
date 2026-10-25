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

# `occ.fillet` on B-rep lofts (the t19 case): a full-circle rim shared
# by a cap face and the periodic tensor-B-spline lateral takes the
# swept rolling-ball path — each rim becomes a periodic tensor-NURBS
# fillet face, the kept faces are rebuilt as trimmed strips and the
# periodic seam edge (the same face on both uses) is skipped per-edge
# inside a list, while a seam-only request reports the upstream
# "no suitable edges" failure.
# Partial arcs and non-circle rims keep the upstream-style
# rejections, and oversized radii carry the spindle/overrun
# certificates.
@testset "occ.fillet loft rims" begin
    A.initialize(); occ=A.model.occ
    try
        # Three-circle smooth loft — the t19 volume.
        c1=occ.add_circle(0,0,0,0.5); l1=occ.add_curve_loop([c1])
        c2=occ.add_circle(0.1,0.05,1,0.1); l2=occ.add_curve_loop([c2])
        c3=occ.add_circle(-0.1,-0.1,2,0.3); l3=occ.add_curve_loop([c3])
        occ.add_thru_sections([l1,l2,l3]); occ.synchronize()
        m0=occ.get_mass(3,1)
        @test m0>0.31

        # The periodic seam spans the full height; the lateral face
        # owns both uses. It is skipped inside a mixed list, but a
        # seam-only request has no suitable edges — matching the
        # OCC "no suitable edges for chamfer or fillet" rejection.
        seam=[t for (d,t) in occ.get_entities(1)
              if (b=occ.get_bounding_box(1,t); b[6]-b[3]>1.0)]
        @test length(seam)==1
        @test occursin("no valid edges",
                       _errtext(()->occ.fillet([1],seam,0.1)))

        # The verbatim t19 edge list — `Unique(Boundary{Surface{f}})`
        # — carries both rims plus the seam; the rims take the swept
        # path and the seam is skipped inline.
        f=A.model.get_boundary([(3,1)],true,false,false)
        etags=sort!(unique!(abs.([t for (d,t) in
            A.model.get_boundary(f,false,false,false)])))
        @test length(etags)==3
        occ.fillet([1],etags,fill(0.1,length(etags))); occ.synchronize()
        m1=occ.get_mass(3,1)
        @test 0.28<m1<m0
        # Rebuilt boundary: trimmed lateral + both caps + 2 fillet
        # pipe faces (the volume's model tables track the rebuild).
        @test length(A.model.get_boundary([(3,1)],true,false,false))==5
        bb=occ.get_bounding_box(3,1)
        @test bb[1]>-0.8 && bb[4]<0.8 && bb[3]>-1e-6 && bb[6]<2.0+1e-6

        # Copy/translate: the filleted B-rep is invariant under
        # rigid copies (trimmed-NURBS CDT quadrature level).
        dup=occ.copy([(3,1)]); occ.translate(dup,5,0,0); occ.synchronize()
        vc=dup[1][2]
        @test occ.get_mass(3,vc)≈m1 atol=2e-4
        g0=occ.get_center_of_mass(3,1); gc=occ.get_center_of_mass(3,vc)
        @test gc[1]≈g0[1]+5 atol=5e-4
        @test gc[2]≈g0[2] atol=5e-4
        @test gc[3]≈g0[3] atol=5e-4

        # Partial-circle rims stay rejected — the sector cylinder's
        # top/bottom rims are half arcs (z-span≈0, wide in x and y).
        v2=occ.add_cylinder(10,0,0,0,0,2,0.5;angle=pi)
        occ.synchronize()
        arcs=[t for (d,t) in occ.get_entities_in_bounding_box(
                9.4,-0.6,-0.1,10.6,0.6,2.1;dim=1)
              if (b=occ.get_bounding_box(1,t);
                  b[6]-b[3]<1e-6 && min(b[4]-b[1],b[5]-b[2])>0.3)]
        @test !isempty(arcs)
        @test occursin("full 2π circle",
                       _errtext(()->occ.fillet([v2],arcs[1:1],0.1)))

        # Non-circle rims keep the upstream rejection text.
        ec=occ.add_circle(20,0,0,0.5); lc=occ.add_curve_loop([ec])
        ee=occ.add_ellipse(20.05,0.02,1.0,0.4,0.25)
        le=occ.add_curve_loop([ee])
        ve=occ.add_thru_sections([lc,le])[1][2]; occ.synchronize()
        @test occursin("must be a line or a full circle",
                       _errtext(()->occ.fillet([ve],[ee],0.05)))

        # Oversized radii keep their certificates: the rolling tube
        # centre collapsing onto the rim axis (spindle) and the
        # tangent trace overrunning a face domain.
        d1=occ.add_circle(30,0,0,0.5); m4=occ.add_curve_loop([d1])
        d2=occ.add_circle(30.1,0.05,1,0.1); m5=occ.add_curve_loop([d2])
        d3=occ.add_circle(29.9,-0.1,2,0.3); m6=occ.add_curve_loop([d3])
        v4=occ.add_thru_sections([m4,m5,m6])[1][2]; occ.synchronize()
        @test occursin("tube centre",_errtext(()->occ.fillet([v4],[d1],0.4)))
        @test occursin("overruns",_errtext(()->occ.fillet([v4],[d1],0.9)))
    finally
        A.finalize()
    end
end
