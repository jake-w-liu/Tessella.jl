# `model.occ.addCurveLoops` — `OCC_Internals::addCurveLoops` searches the
# given curves for closed loops (`ShapeAnalysis_FreeBounds` semantics:
# endpoints connect on `Precision::Confusion` coincidence, not just shared
# vertex entities). Each closed walk publishes a curve loop under a fresh
# tag; unconsumed curves come back as leftovers.

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

function _poly_edges(occ,pts;closed=true)
    n=length(pts)
    ps=[occ.add_point(p...) for p in pts]
    [occ.add_line(ps[i],ps[mod1(i+1,n)]) for i in 1:(closed ? n : n-1)]
end

@testset "occ.add_curve_loops" begin
    occ=_occ_setup()
    try
        # Plain closed square -> one loop, no leftovers.
        es=_poly_edges(occ,[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)])
        loops,leftover=occ.add_curve_loops(es)
        @test length(loops)==1
        @test isempty(leftover)
        @test loops[1] isa Int
        # The found loop is a real closed wire usable by add_plane_surface.
        f=occ.add_plane_surface([loops[1]])
        @test occ.get_mass(2,f)≈1.0

        # Square plus a dangling edge -> loop + one leftover.
        es2=_poly_edges(occ,[(0.,0.,2.),(1.,0.,2.),(1.,1.,2.),(0.,1.,2.)])
        tail=occ.add_line(occ.add_point(2.0,0.0,2.0),
                          occ.add_point(3.0,0.0,2.0))
        loops,leftover=occ.add_curve_loops(vcat(es2,tail))
        @test length(loops)==1
        @test leftover==[tail]

        # Two disjoint triangles -> two loops.
        t1=_poly_edges(occ,[(0.,0.,4.),(1.,0.,4.),(0.,1.,4.)])
        t2=_poly_edges(occ,[(2.,0.,4.),(3.,0.,4.),(2.,1.,4.)])
        loops,leftover=occ.add_curve_loops(vcat(t1,t2))
        @test length(loops)==2
        @test isempty(leftover)

        # Purely open chain -> no loops, everything returns.
        open_es=_poly_edges(occ,[(0.,0.,6.),(1.,0.,6.),(2.,0.,6.)];
                            closed=false)
        loops,leftover=occ.add_curve_loops(open_es)
        @test isempty(loops)
        @test sort(leftover)==sort(open_es)

        # Single closed circle edge is a loop by itself.
        c=occ.add_circle(0.0,0.0,8.0,1.0)
        loops,leftover=occ.add_curve_loops([c])
        @test length(loops)==1
        @test isempty(leftover)

        # Figure-8: two squares sharing only vertex v8 — discovery must not
        # fuse them into one self-touching wire.
        v8=occ.add_point(4.0,0.0,0.0)
        sq1=[occ.add_line(v8,occ.add_point(5.0,0.0,0.0)),
             occ.add_line(occ.add_point(5.0,0.0,0.0),
                          occ.add_point(5.0,1.0,0.0)),
             occ.add_line(occ.add_point(5.0,1.0,0.0),
                          occ.add_point(4.0,1.0,0.0)),
             occ.add_line(occ.add_point(4.0,1.0,0.0),v8)]
        sq2=[occ.add_line(v8,occ.add_point(3.0,0.0,0.0)),
             occ.add_line(occ.add_point(3.0,0.0,0.0),
                          occ.add_point(3.0,-1.0,0.0)),
             occ.add_line(occ.add_point(3.0,-1.0,0.0),
                          occ.add_point(4.0,-1.0,0.0)),
             occ.add_line(occ.add_point(4.0,-1.0,0.0),v8)]
        # Interleave the edge order so the walk is tempted to mix loops.
        loops,leftover=occ.add_curve_loops(
            [sq1[1],sq2[1],sq1[2],sq2[2],sq1[3],sq2[3],sq1[4],sq2[4]])
        @test length(loops)==2
        @test isempty(leftover)
        for loop in loops
            f=occ.add_plane_surface([loop])
            @test occ.get_mass(2,f)≈1.0
        end

        # Coincident-but-distinct vertices still connect (FreeBounds).
        q1=[occ.add_line(occ.add_point(0.0,0.0,10.0),
                         occ.add_point(1.0,0.0,10.0)),
            occ.add_line(occ.add_point(1.0,0.0,10.0),
                         occ.add_point(1.0,1.0,10.0)),
            occ.add_line(occ.add_point(1.0,1.0,10.0),
                         occ.add_point(0.0,1.0,10.0)),
            occ.add_line(occ.add_point(0.0,1.0,10.0),
                         occ.add_point(0.0,0.0,10.0))]
        loops,leftover=occ.add_curve_loops(q1)
        @test length(loops)==1
        @test isempty(leftover)

        # Unknown tags reject atomically; empty input returns empty.
        @test occursin("unknown OCC Curve",
                       _errtext(()->occ.add_curve_loops([99999])))
        loops,leftover=occ.add_curve_loops(Int[])
        @test isempty(loops) && isempty(leftover)
    finally
        A.finalize()
    end
end
