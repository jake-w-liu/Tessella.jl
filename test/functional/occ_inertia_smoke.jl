module OccInertiaSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC kernel-side matrix of inertia — registry read only; the
        # interface test exercises per-dim curves, quadric faces, plane
        # tensors, primitive and divergence-path solids, the wedge fan,
        # transforms, degenerate edges and the unsupported-boundary and
        # empty-GProp fall-throughs.
        v=occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
        @test occ.get_matrix_of_inertia(3,v)≈
            [2/3,-0.25,-0.25, -0.25,2/3,-0.25, -0.25,-0.25,2/3] atol=1e-12
        s=occ.add_sphere(0.0,0.0,0.0,1.0)
        @test occ.get_matrix_of_inertia(3,s)≈
            [8π/15,0.0,0.0, 0.0,8π/15,0.0, 0.0,0.0,8π/15] atol=1e-12
        p1=occ.add_point(1.0,0.0,0.0); p2=occ.add_point(3.0,0.0,0.0)
        l=occ.add_line(p1,p2)
        @test occ.get_matrix_of_inertia(1,l)≈
            [0.0,0.0,0.0, 0.0,26/3,0.0, 0.0,0.0,26/3]
        @test occ.get_matrix_of_inertia(0,p1)≈zeros(9)
        @test_throws ArgumentError occ.get_matrix_of_inertia(3,99)
        occ.synchronize()
        @test occ.get_matrix_of_inertia(3,v)[9]≈2/3 atol=1e-12
    finally
        api.finalize()
    end
end
end
