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
using Tessella.CAD

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
end
