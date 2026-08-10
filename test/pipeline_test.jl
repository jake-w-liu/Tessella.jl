# ── Integration CRC suite: the top-level mesh_volume pipeline ───────────────────
#
# Correctness  : end-to-end surface → volume mesh at the exact expected volume.
# Robustness   : the "validated or explicit blocker" contract — a defective input
#                raises a precise error, never a silent bad mesh.
# Completeness : output validates; smoothing preserves volume.

using Test
using Tessella
using Tessella.MeshTypes

mvpvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                for t in 1:ntets(m); init=0.0)

function _cube_surface()
    C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
    t=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    Mesh(C; tris=t)
end

@testset "mesh_volume pipeline (integration)" begin
    @testset "clean surface → validated volume mesh" begin
        m = mesh_volume(_cube_surface())
        @test validate(m).ok
        @test mvpvol(m) ≈ 1.0 rtol=1e-6
        @test mesh_quality(m).n_tets == ntets(m)
    end

    @testset "defective surface → explicit blocker (not silent)" begin
        cube = _cube_surface()
        openm = Mesh(cube.coords; tris=cube.tris[:, 1:end-1])   # missing a face
        @test_throws ArgumentError mesh_volume(openm)
        # check=false bypasses the gate (caller takes responsibility)
        @test mesh_volume(openm; check=false) isa Mesh
    end

    @testset "smoothing preserves volume through the pipeline" begin
        m1 = mesh_volume(_cube_surface(); smooth=false)
        m2 = mesh_volume(_cube_surface(); smooth=true, smooth_iters=8)
        @test mvpvol(m1) ≈ mvpvol(m2) rtol=1e-6
    end
end
