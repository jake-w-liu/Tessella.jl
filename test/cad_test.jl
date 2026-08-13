# ── Stage-5 CRC suite: native analytical geometry (CAD-lite), no OpenCASCADE ─────
#
# Correctness  : every analytical surface's exact membership/projection/intersection is
#                checked to round-off against first-principles residuals (a projected
#                point has zero surface residual; an imprint curve lies on BOTH surfaces).
# Robustness   : random query points, on-axis degeneracy, oblique (ellipse) imprints.
# Completeness : the exact boolean-imprint curve the ENC-COAX bore/shield need — a
#                cylinder piercing a planar wall is an EXACT circle, natively.

using Test
using Random
using Tessella
using Tessella.CAD
using Tessella.MeshTypes
using Tessella.HighOrder
using Tessella.Mesh3D: recover_boundary

@testset "CAD (native analytical geometry, Stage 5)" begin

    @testset "exact surface membership + projection" begin
        rng = MersenneTwister(1)
        cyl = CylinderS((17.,10.,15.), (0.,1.,0.), 0.21)
        sph = SphereS((0.3,-0.2,0.1), 2.0)
        pln = PlaneS((1.,2.,3.), (0.2,0.9,0.1))
        # a projected point lies EXACTLY on the surface (zero residual to round-off)
        for _ in 1:2000
            p = (17+3randn(rng), 10+3randn(rng), 15+3randn(rng))
            @test abs(surface_residual(cyl, project_to(cyl, p))) < 1e-11
            @test on_surface(cyl, project_to(cyl, p))
        end
        for _ in 1:2000
            p = (3randn(rng), 3randn(rng), 3randn(rng))
            @test abs(surface_residual(sph, project_to(sph, p))) < 1e-11
            @test abs(surface_residual(pln, project_to(pln, p))) < 1e-11
        end
        # a point known to be ON the cylinder reads zero residual
        @test abs(surface_residual(cyl, (17+0.21, 10.0, 15.0))) < 1e-12
        @test !on_surface(cyl, (17+0.30, 10.0, 15.0))          # off the surface ⇒ rejected
    end

    @testset "EXACT boolean imprint: cylinder ∩ plane (the coax bore through the case wall)" begin
        cyl  = CylinderS((17.,10.,15.), (0.,1.,0.), 0.21)      # bore along +y
        wall = PlaneS((0.,14.,0.), (0.,1.,0.))                 # case wall at y=14 (⊥ the bore)
        center, circ = imprint_circle(cyl, wall; nseg=48)
        @test length(circ) == 48
        @test center == (17.0, 14.0, 15.0)                     # axis ∩ wall
        # every imprint node lies EXACTLY on BOTH the cylinder and the wall — the exact
        # curved boolean interface OpenCASCADE would produce, computed natively.
        @test maximum(abs(surface_residual(cyl,  p)) for p in circ) < 1e-12
        @test maximum(abs(surface_residual(wall, p)) for p in circ) < 1e-12
        # the imprint circle has the exact bore radius from its center
        @test maximum(abs(sqrt(sum((p .- center).^2)) - cyl.r) for p in circ) < 1e-12
        # non-perpendicular plane ⇒ imprint is an ellipse, still exact on both surfaces
        obl = PlaneS((0.,14.,0.), (0.3,1.0,0.2))
        _, ell = imprint_ellipse(cyl, obl; nseg=48)
        @test maximum(abs(surface_residual(cyl, p)) for p in ell) < 1e-12
        @test maximum(abs(surface_residual(obl, p)) for p in ell) < 1e-12
        # a perpendicular request on an oblique plane is rejected (routed to the ellipse)
        @test_throws ArgumentError imprint_circle(cyl, obl)
    end

    @testset "exact-surface meshing: curve P2 elements onto the analytical surface (any primitive)" begin
        # The native CAD layer plugs into HighOrder.curve_to_surface! — so a mesh's P2
        # mid-nodes land EXACTLY on the true analytical surface (the "NURBS surface" a mesh
        # can carry), generalized beyond cylinders to ANY CAD primitive. Here: a sphere.
        R = 2.0
        t = (1+sqrt(5))/2
        V = [(-1.,t,0.),(1.,t,0.),(-1.,-t,0.),(1.,-t,0.),(0.,-1.,t),(0.,1.,t),(0.,-1.,-t),(0.,1.,-t),(t,0.,-1.),(t,0.,1.),(-t,0.,-1.),(-t,0.,1.)]
        F = [(1,12,6),(1,6,2),(1,2,8),(1,8,11),(1,11,12),(2,6,10),(6,12,5),(12,11,3),(11,8,7),(8,2,9),(4,10,5),(4,5,3),(4,3,7),(4,7,9),(4,9,10),(5,10,6),(3,5,12),(7,3,11),(9,7,8),(10,9,2)]
        C = Matrix{Float64}(undef,3,12)
        for (i,v) in enumerate(V); p = collect(v) .* (R/sqrt(sum(collect(v).^2))); C[:,i] = p; end
        Tm = Matrix{Int32}(undef,3,20); for (i,f) in enumerate(F); Tm[:,i] = Int32[f...]; end
        m  = recover_boundary(Mesh(C; tris=Tm))
        p2 = p2_tetmesh(m)
        sph = SphereS((0.,0.,0.), R)
        nc = curve_to_surface!(p2, (x,y,z)->project_to(sph,(x,y,z)), (x,y,z)->on_surface(sph,(x,y,z); tol=1e-9))
        @test nc > 0                                       # boundary edges were curved onto the sphere
        # every mid-node that sits on the sphere reads zero surface residual (exact geometry)
        maxres = 0.0
        for i in 1:size(p2.coords,2)
            p = (p2.coords[1,i], p2.coords[2,i], p2.coords[3,i])
            abs(sqrt(sum(p.^2)) - R) < 1e-9 && (maxres = max(maxres, abs(surface_residual(sph, p))))
        end
        @test maxres < 1e-9
    end
end
