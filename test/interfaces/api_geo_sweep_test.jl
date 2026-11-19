module ApiGeoSweepTest

using Test, Tessella

# Built-in `geo.extrude`/`revolve`/`twist` with structured layers, the
# `geo.mesh.*` attribute setters, and positional `add_point`. Expected
# `outDimTags` and mesh counts are Gmsh 4.15.2's for the same calls.
@testset "API built-in geo sweeps" begin
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        g = A.model.geo
        p = [g.add_point(x, y, 0, 0.5) for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1))]
        l = [g.add_line(p[i], p[mod1(i + 1, 4)]) for i in 1:4]
        s = g.add_plane_surface([g.add_curve_loop(l)])
        for c in l
            g.mesh.set_transfinite_curve(c, 4)
        end
        g.mesh.set_transfinite_surface(s)
        g.mesh.set_recombine(2, s)
        @test g.extrude([(2, s)], 0, 0, 1; numElements=[2, 3], heights=[0.4, 1],
                        recombine=true) ==
              [(2, 26), (3, 1), (2, 13), (2, 17), (2, 21), (2, 25)]
        # A curve on the rotation axis is degenerate: Gmsh returns the curve.
        @test g.revolve([(1, l[2])], 1, 0, 0, 0, 1, 0, pi / 2; numElements=[4],
                        recombine=true) == [(1, 2)]
        @test g.twist([(0, p[1])], 0, 0, 0, 0, 0, 1, 0, 0, 1, pi / 4;
                      numElements=[3]) == [(0, 5), (1, 28)]
        g.synchronize()
        @test [t for (_, t) in A.model.get_entities(0)] ==
              [1, 2, 3, 4, 5, 6, 10, 14, 21, 22, 23, 24]
        A.mesh.generate(3)
        @test length(A.mesh.get_nodes()[1]) == 102
        types, tags, _ = A.mesh.get_elements()
        counts = Dict(t => length(e) for (t, e) in zip(types, tags))
        @test counts == Dict(1 => 47, 3 => 78, 5 => 45, 15 => 12)
    finally
        A.finalize()
    end
end

end # module
