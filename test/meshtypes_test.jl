# ── Stage-0 CRC suite: mesh container, topology, quality, checksum ──────────────
#
# Correctness  : geometry checked against independent oracles — Heron's-formula
#                area, a linear-solve (Cramer) circumcenter, analytic volumes,
#                hand-computed regular-tetrahedron dihedral/radius-edge, and Euler
#                characteristics of canonical complexes (single tet, unit cube).
# Robustness   : degenerate (flat/inverted) cells, empty meshes, isolated nodes.
# Completeness : validate() flags inverted/degenerate/non-manifold; mesh_crc is
#                deterministic, order-invariant, and mutation-sensitive.

using Test
using Tessella.MeshTypes

norm(v) = sqrt(sum(x^2 for x in v))

# ── independent oracles ─────────────────────────────────────────────────────────
heron(a, b, c) = begin
    A = norm(collect(b) .- collect(a)); B = norm(collect(c) .- collect(b)); C = norm(collect(a) .- collect(c))
    s = (A + B + C) / 2
    sqrt(max(s*(s-A)*(s-B)*(s-C), 0.0))
end

# Circumcenter by solving 2(v_j - v_1)·O = |v_j|²-|v_1|² with Cramer's rule.
function oracle_circumradius(a, b, c, d)
    v1 = collect(Float64, a); v2 = collect(Float64, b); v3 = collect(Float64, c); v4 = collect(Float64, d)
    A = [ (v2 .- v1)'; (v3 .- v1)'; (v4 .- v1)' ] .* 2
    rhs = [sum(v2.^2)-sum(v1.^2), sum(v3.^2)-sum(v1.^2), sum(v4.^2)-sum(v1.^2)]
    det(M) = M[1,1]*(M[2,2]*M[3,3]-M[2,3]*M[3,2]) -
             M[1,2]*(M[2,1]*M[3,3]-M[2,3]*M[3,1]) +
             M[1,3]*(M[2,1]*M[3,2]-M[2,2]*M[3,1])
    D = det(A)
    abs(D) < 1e-300 && return Inf
    cx = det([rhs A[:,2] A[:,3]]) / D
    cy = det([A[:,1] rhs A[:,3]]) / D
    cz = det([A[:,1] A[:,2] rhs]) / D
    return norm([cx,cy,cz] .- v1)
end

# Positively-oriented 6-tet Kuhn decomposition of the unit cube.
function cube_mesh()
    coords = Float64[0 1 1 0 0 1 1 0;   # x
                     0 0 1 1 0 0 1 1;   # y
                     0 0 0 0 1 1 1 1]   # z
    raw = [ (1,2,3,7), (1,2,6,7), (1,4,3,7), (1,4,8,7), (1,5,6,7), (1,5,8,7) ]
    tets = Matrix{Int32}(undef, 4, length(raw))
    for (k, t) in enumerate(raw)
        a = (coords[1,t[1]],coords[2,t[1]],coords[3,t[1]])
        b = (coords[1,t[2]],coords[2,t[2]],coords[3,t[2]])
        c = (coords[1,t[3]],coords[2,t[3]],coords[3,t[3]])
        d = (coords[1,t[4]],coords[2,t[4]],coords[3,t[4]])
        v = tet_signed_volume(a,b,c,d)
        tt = v < 0 ? (t[1],t[2],t[4],t[3]) : t   # flip to positive orientation
        tets[:,k] = Int32[tt...]
    end
    return Mesh(coords; tets=tets)
end

@testset "MeshTypes (Stage 0)" begin

    @testset "accessors + allocation-free node()" begin
        m = cube_mesh()
        @test nnodes(m) == 8
        @test ntets(m) == 6
        @test node(m, 7) == (1.0, 1.0, 1.0)
        @test node2(m, 3) == (1.0, 1.0)
        node(m, 1); node2(m, 1)
        @test (@allocated node(m, 1)) == 0
        @test (@allocated node2(m, 1)) == 0
    end

    @testset "triangle_area vs Heron" begin
        for (a,b,c) in [ ((0,0,0),(3,0,0),(0,4,0)),        # 3-4-5 → area 6
                         ((0,0,0),(1,0,0),(0,1,0)),        # area 0.5
                         ((1,2,3),(4,0,1),(2,2,5)),
                         ((0,0,0),(2,0,0),(1,0,0)) ]       # degenerate → 0
            @test triangle_area(a,b,c) ≈ heron(a,b,c) atol=1e-12
        end
        @test triangle_area((0,0,0),(3,0,0),(0,4,0)) ≈ 6.0
        @test triangle_area((0,0,0),(2,0,0),(1,0,0)) == 0.0
    end

    @testset "tet volume analytic + sign" begin
        @test tet_volume((0,0,0),(1,0,0),(0,1,0),(0,0,1)) ≈ 1/6
        @test tet_signed_volume((0,0,0),(1,0,0),(0,1,0),(0,0,1)) ≈ 1/6   # positive
        @test tet_signed_volume((0,0,0),(1,0,0),(0,0,1),(0,1,0)) ≈ -1/6  # swapped → negative
        @test tet_signed_volume((0,0,0),(1,0,0),(2,0,0),(3,1,0)) == 0.0  # coplanar → 0
    end

    @testset "circumradius vs Cramer oracle + known values" begin
        @test tet_circumradius((0,0,0),(1,0,0),(0,1,0),(0,0,1)) ≈ sqrt(3)/2
        @test tet_circumradius((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)) ≈ sqrt(3)
        # random non-degenerate tets vs the linear-solve oracle
        st = UInt64(0x9E3779B97F4A7C15)
        nf() = (st ⊻= st<<13; st ⊻= st>>7; st ⊻= st<<17; (st>>11)/Float64(2^53) - 0.5)
        tested = 0
        for _ in 1:2000
            a=(nf(),nf(),nf()); b=(nf(),nf(),nf()); c=(nf(),nf(),nf()); d=(nf(),nf(),nf())
            abs(tet_signed_volume(a,b,c,d)) < 1e-3 && continue
            @test tet_circumradius(a,b,c,d) ≈ oracle_circumradius(a,b,c,d) rtol=1e-9
            tested += 1
        end
        @test tested > 500
    end

    @testset "regular-tet quality (hand values)" begin
        a,b,c,d = (1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)
        mn, mx = tet_dihedral_extrema(a,b,c,d)
        @test mn ≈ acos(1/3) atol=1e-12        # all six dihedrals equal
        @test mx ≈ acos(1/3) atol=1e-12
        @test tet_radius_edge(a,b,c,d) ≈ sqrt(3)/(2*sqrt(2)) rtol=1e-12
        # right-corner tet has a π/2 dihedral present
        mn2, mx2 = tet_dihedral_extrema((0,0,0),(1,0,0),(0,1,0),(0,0,1))
        @test mn2 > 0 && mx2 < pi
    end

    @testset "boundary extraction" begin
        # single tet: 4 boundary faces, manifold
        one = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,4],4,1))
        bf, maxinc = boundary_faces(one.tets)
        @test length(bf) == 4
        @test maxinc == 1
        # unit cube: surface = 12 triangles, every internal face shared by 2 tets
        m = cube_mesh()
        bf2, maxinc2 = boundary_faces(m.tets)
        @test length(bf2) == 12
        @test maxinc2 == 2
        # triangle surface boundary edges: single triangle → 3 boundary edges
        tsurf = Mesh(Float64[0 1 0; 0 0 1; 0 0 0]; tris=reshape(Int32[1,2,3],3,1))
        be, _ = boundary_edges(tsurf.tris)
        @test length(be) == 3
    end

    @testset "Euler characteristic" begin
        one = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,4],4,1))
        @test euler_characteristic(one) == 1     # ball
        @test boundary_euler(one) == 2           # boundary sphere
        m = cube_mesh()
        @test euler_characteristic(m) == 1       # cube ≅ ball
        @test boundary_euler(m) == 2             # cube surface ≅ sphere
        @test is_closed_manifold(m)
        # segment nodes are NOT counted in V (E carries no segment edges): a tet plus
        # a disjoint segment still reads χ=1, not 3 (regression).
        seg = Mesh(Float64[0 1 0 0 5 6; 0 0 1 0 0 0; 0 0 0 1 0 0];
                   tets=reshape(Int32[1,2,3,4],4,1), segs=reshape(Int32[5,6],2,1))
        @test euler_characteristic(seg) == 1
    end

    @testset "bounding_box / mesh_crc exclude unreferenced nodes" begin
        # node 5 is isolated far away; the contract is "over referenced nodes", and
        # the bbox feeds mesh_crc — a stray node must not enlarge either.
        m = Mesh(Float64[0 1 0 0 1000; 0 0 1 0 0; 0 0 0 1 0]; tets=reshape(Int32[1,2,3,4],4,1))
        lo, hi = bounding_box(m)
        @test lo == (0.0, 0.0, 0.0)
        @test hi == (1.0, 1.0, 1.0)                       # not (1000,1,1)
        @test mesh_crc(m).bbox == ((0.0,0.0,0.0),(1.0,1.0,1.0))
        # referenced segment nodes DO count toward the geometric bbox
        ms = Mesh(Float64[0 1 0 0 2; 0 0 1 0 3; 0 0 0 1 0];
                  tets=reshape(Int32[1,2,3,4],4,1), segs=reshape(Int32[1,5],2,1))
        @test bounding_box(ms)[2] == (2.0, 3.0, 1.0)
    end

    @testset "validate contract" begin
        m = cube_mesh()
        d = validate(m)
        @test d.ok
        @test isempty(d.messages)
        # inverted tet → flagged
        bad = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,4,3],4,1))
        db = validate(bad)
        @test !db.ok
        @test any(occursin("negatively-oriented", s) for s in db.messages)
        # degenerate (flat) tet → flagged
        flat = Mesh(Float64[0 1 2 3; 0 0 0 1; 0 0 0 0]; tets=reshape(Int32[1,2,3,4],4,1))
        @test !validate(flat).ok
    end

    @testset "constructor rejects bad input" begin
        @test_throws ArgumentError Mesh(Float64[0 1; 0 0]; )                 # coords not 3×n
        @test_throws ArgumentError Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,9],4,1))  # id out of range
        @test_throws ArgumentError Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tris=reshape(Int32[1,2,3],3,1), tri_tag=Int32[1,2])  # tag mismatch
    end

    @testset "mesh_crc determinism + order-invariance + mutation-sensitivity" begin
        m = cube_mesh()
        crc1 = mesh_crc(m)
        crc2 = mesh_crc(cube_mesh())
        @test crc1.sha == crc2.sha
        @test crc1.n_tets == 6 && crc1.n_nodes == 8
        @test crc1.n_boundary_faces == 12

        # permuting vertices within tets and reordering tets must NOT change the sha
        tets = copy(m.tets)
        perm = [3,1,2,4]                       # even perm keeps orientation class
        tets2 = tets[perm, [6,5,4,3,2,1]]      # reorder cells + permute within
        m2 = Mesh(m.coords; tets=tets2)
        @test mesh_crc(m2).sha == crc1.sha

        # changing one connectivity entry MUST change the sha (mutation sensitivity)
        tets3 = copy(m.tets); tets3[4, 1] = tets3[4,1] == 8 ? Int32(6) : Int32(8)
        m3 = Mesh(m.coords; tets=tets3)
        @test mesh_crc(m3).sha != crc1.sha

        # crc_string is printable
        @test occursin("tets=6", crc_string(crc1))
    end

    @testset "empty mesh is well-defined" begin
        e = Mesh(Matrix{Float64}(undef, 3, 0))
        @test nnodes(e) == 0 && ntets(e) == 0
        c = mesh_crc(e)
        @test c.n_tets == 0
        @test validate(e).ok           # vacuously valid; empty-region detection is caller's job
    end
end
