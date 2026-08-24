# ── Stage-5 CRC suite: surface diagnostics ("heal, don't fail" detection) ───────
#
# Correctness  : a clean closed cube surface passes every check; each injected
#                defect (open edge, non-manifold, degenerate, duplicate, flipped
#                orientation, coincident vertex) is detected exactly.
# Robustness   : coincident-vertex hash grid; tolerance scaling.
# Completeness : is_meshable gates volume meshing (go/no-go with a precise report).

using Test
using Tessella.MeshTypes
using Tessella.Heal

# clean unit-cube surface (12 outward-consistent triangles)
function clean_cube()
    C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
    t=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    Mesh(C; tris=t)
end

@testset "Heal surface diagnostics (Stage 5)" begin

    @testset "clean cube passes every check" begin
        m = clean_cube()
        ok, r = is_meshable(m)
        @test ok
        @test r.closed && r.manifold && r.oriented
        @test r.n_open_edges == 0 && r.n_nonmanifold_edges == 0
        @test r.n_degenerate_tris == 0 && r.n_duplicate_tris == 0 && r.n_coincident_pairs == 0
        @test r.min_edge_length ≈ 1.0
        @test isempty(r.messages)
    end

    @testset "open surface (missing face) detected" begin
        m = clean_cube()
        open_tris = m.tris[:, 1:end-1]              # drop one triangle
        mo = Mesh(m.coords; tris=open_tris)
        ok, r = is_meshable(mo)
        @test !ok
        @test !r.closed
        @test r.n_open_edges > 0
        @test any(occursin("not closed", s) for s in r.messages)
    end

    @testset "flipped-orientation triangle detected" begin
        m = clean_cube()
        tris = copy(m.tris)
        tris[2,1], tris[3,1] = tris[3,1], tris[2,1]   # reverse triangle 1
        mf = Mesh(m.coords; tris=tris)
        ok, r = is_meshable(mf)
        @test !ok
        @test r.closed && r.manifold                  # still closed & manifold…
        @test !r.oriented                             # …but orientation inconsistent
    end

    @testset "degenerate triangle detected" begin
        # add a zero-area triangle (3 collinear points)
        C = Float64[0 1 2; 0 0 0; 0 0 0]
        m = Mesh(C; tris=reshape(Int32[1,2,3],3,1))
        r = surface_diagnostics(m)
        @test r.n_degenerate_tris == 1
    end

    @testset "duplicate triangle detected" begin
        C = Float64[0 1 0; 0 0 1; 0 0 0]
        m = Mesh(C; tris=Int32[1 1; 2 2; 3 3])       # same triangle twice
        r = surface_diagnostics(m)
        @test r.n_duplicate_tris == 1
        @test !r.oriented                            # both traverse edges the same way
    end

    @testset "non-manifold edge detected" begin
        # edge (1,2) shared by three triangles
        C = Float64[0 1 0 1 0.5; 0 0 1 1 -1; 0 0 0 0 0]
        tris = Int32[1 1 1; 2 2 2; 3 4 5]
        m = Mesh(C; tris=tris)
        r = surface_diagnostics(m)
        @test r.n_nonmanifold_edges >= 1
        @test !r.manifold
    end

    @testset "vertex-pinch topology detected" begin
        C=Float64[0 1 0 -1 0;0 0 1 0 -1;0 0 0 0 0]
        m=Mesh(C;tris=Int32[1 1;2 4;3 5])
        r=surface_diagnostics(m)
        @test !r.manifold
        @test r.n_nonmanifold_edges==0
        @test any(occursin("disconnected triangle link",s) for s in r.messages)
    end

    @testset "coincident vertices detected" begin
        # vertex 9 coincides with vertex 1
        m = clean_cube()
        C = hcat(m.coords, m.coords[:,1])            # duplicate node 1 as node 9
        m2 = Mesh(C; tris=m.tris)
        r = surface_diagnostics(m2; tol=1e-9)
        @test r.n_coincident_pairs >= 1
    end

    @testset "empty, non-finite, invalid-tolerance, and far-origin inputs block safely" begin
        empty_surface = Mesh(Matrix{Float64}(undef, 3, 0))
        ok, r = is_meshable(empty_surface)
        @test !ok && !r.closed
        @test any(occursin("no triangles", s) for s in r.messages)

        for badtol in (-1.0, Inf, NaN)
            @test_throws ArgumentError surface_diagnostics(empty_surface; tol=badtol)
        end

        cube = clean_cube()
        far = Mesh(cube.coords .+ 1.0e15; tris=cube.tris)
        okfar, rfar = is_meshable(far)
        @test okfar
        @test rfar.n_coincident_pairs == 0
        @test rfar.min_edge_length ≈ 1.0
        far_dup = Mesh(hcat(far.coords, far.coords[:,1]); tris=far.tris)
        @test surface_diagnostics(far_dup).n_coincident_pairs == 1

        C = copy(cube.coords); C[1,1] = NaN
        nanmesh=Mesh(cube.coords;tris=cube.tris);nanmesh.coords[1,1]=NaN
        oknan, rnan = is_meshable(nanmesh)
        @test !oknan
        @test any(occursin("non-finite coordinates", s) for s in rnan.messages)

        infmesh=Mesh(cube.coords;tris=cube.tris);infmesh.coords[1,1]=Inf
        okinf, rinf = is_meshable(infmesh)
        @test !okinf
        @test any(occursin("non-finite coordinates", s) for s in rinf.messages)
    end
end
