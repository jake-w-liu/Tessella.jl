# ── Stage-0 CRC suite: .msh v2/v4 round-trip, STL ingest, .geo scan ─────────────
#
# Correctness  : hand-written gmsh-format samples (an independent transcription of
#                the documented format) are parsed and checked node-by-node /
#                element-by-element; round-trips are checked by mesh_crc equality.
# Robustness   : non-contiguous node tags, physical groups, empty mesh, both
#                format versions, ASCII+binary STL, welding of coincident verts.
# Completeness : cross-format (v2→v4→v2) preserves connectivity CRC exactly.

using Test
using Tessella.MeshTypes
using Tessella.IO

# Positively-oriented 6-tet Kuhn cube (same helper as meshtypes_test, redefined
# locally so this file runs standalone).
function _cube()
    coords = Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    raw = [ (1,2,3,7), (1,2,6,7), (1,4,3,7), (1,4,8,7), (1,5,6,7), (1,5,8,7) ]
    tets = Matrix{Int32}(undef, 4, length(raw))
    tags = Int32[]
    for (k, t) in enumerate(raw)
        f(i) = (coords[1,i],coords[2,i],coords[3,i])
        v = tet_signed_volume(f(t[1]),f(t[2]),f(t[3]),f(t[4]))
        tt = v < 0 ? (t[1],t[2],t[4],t[3]) : t
        tets[:,k] = Int32[tt...]
        push!(tags, Int32(4))     # physical tag 4 ("case")
    end
    return Mesh(coords; tets=tets, tet_tag=tags)
end

@testset "IO (Stage 0)" begin
    dir = mktempdir()

    @testset "v2 round-trip preserves connectivity CRC" begin
        m = _cube()
        p = joinpath(dir, "cube_v2.msh")
        write_msh(p, m; version=2.2, physical_names=Dict((3,4)=>"case"))
        f = read_msh(p)
        @test mesh_crc(f.mesh).sha == mesh_crc(m).sha
        @test ntets(f.mesh) == 6
        @test all(f.mesh.tet_tag .== 4)
        @test f.physical_names[(3,4)] == "case"
    end

    @testset "v4.1 round-trip preserves connectivity CRC" begin
        m = _cube()
        p = joinpath(dir, "cube_v4.msh")
        write_msh(p, m; version=4.1, physical_names=Dict((3,4)=>"case"))
        f = read_msh(p)
        @test mesh_crc(f.mesh).sha == mesh_crc(m).sha
        @test ntets(f.mesh) == 6
        @test all(f.mesh.tet_tag .== 4)
        @test f.physical_names[(3,4)] == "case"
    end

    @testset "cross-format v2 → v4 → v2 keeps CRC" begin
        m = _cube()
        p2 = joinpath(dir, "x_v2.msh"); p4 = joinpath(dir, "x_v4.msh"); p2b = joinpath(dir, "x_v2b.msh")
        write_msh(p2, m; version=2.2)
        m2 = read_msh(p2).mesh
        write_msh(p4, m2; version=4.1)
        m4 = read_msh(p4).mesh
        write_msh(p2b, m4; version=2.2)
        m2b = read_msh(p2b).mesh
        @test mesh_crc(m).sha == mesh_crc(m2b).sha
    end

    @testset "mixed-dim mesh (segs+tris+tets) round-trip" begin
        coords = Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
        m = Mesh(coords;
                 segs=reshape(Int32[1,2],2,1), seg_tag=Int32[6],
                 tris=reshape(Int32[1,2,3],3,1), tri_tag=Int32[3],
                 tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[1])
        for ver in (2.2, 4.1)
            p = joinpath(dir, "mixed_$(ver).msh")
            write_msh(p, m; version=ver)
            f = read_msh(p)
            @test mesh_crc(f.mesh).sha == mesh_crc(m).sha
            @test nsegs(f.mesh)==1 && ntris(f.mesh)==1 && ntets(f.mesh)==1
            @test f.mesh.seg_tag[1]==6 && f.mesh.tri_tag[1]==3 && f.mesh.tet_tag[1]==1
        end
    end

    @testset "parse hand-written gmsh v2 sample (format oracle)" begin
        sample = """
        \$MeshFormat
        2.2 0 8
        \$EndMeshFormat
        \$PhysicalNames
        1
        2 5 "surf"
        \$EndPhysicalNames
        \$Nodes
        3
        10 0.0 0.0 0.0
        20 2.0 0.0 0.0
        30 0.0 3.0 0.0
        \$EndNodes
        \$Elements
        1
        1 2 2 5 1 10 20 30
        \$EndElements
        """
        p = joinpath(dir, "sample_v2.msh"); write(p, sample)
        f = read_msh(p)
        @test nnodes(f.mesh) == 3
        @test ntris(f.mesh) == 1
        @test f.mesh.tri_tag[1] == 5
        @test f.physical_names[(2,5)] == "surf"
        # non-contiguous tags (10,20,30) relabelled but geometry intact:
        @test triangle_area(node(f.mesh,f.mesh.tris[1,1]), node(f.mesh,f.mesh.tris[2,1]), node(f.mesh,f.mesh.tris[3,1])) ≈ 3.0
    end

    @testset "parse hand-written gmsh v4.1 sample (format oracle)" begin
        sample = """
        \$MeshFormat
        4.1 0 8
        \$EndMeshFormat
        \$PhysicalNames
        1
        3 1 "vol"
        \$EndPhysicalNames
        \$Entities
        0 0 0 1
        1 0 0 0 1 1 1 1 1 0
        \$EndEntities
        \$Nodes
        1 4 1 4
        3 1 0 4
        1
        2
        3
        4
        0 0 0
        1 0 0
        0 1 0
        0 0 1
        \$EndNodes
        \$Elements
        1 1 1 1
        3 1 4 1
        1 1 2 3 4
        \$EndElements
        """
        p = joinpath(dir, "sample_v4.msh"); write(p, sample)
        f = read_msh(p)
        @test nnodes(f.mesh) == 4
        @test ntets(f.mesh) == 1
        @test f.mesh.tet_tag[1] == 1
        @test f.physical_names[(3,1)] == "vol"
        @test tet_volume(node(f.mesh,1),node(f.mesh,2),node(f.mesh,3),node(f.mesh,4)) ≈ 1/6
        @test validate(f.mesh).ok
    end

    @testset "STL ingest welds coincident vertices (ASCII)" begin
        # tetrahedron surface: 4 triangles, 4 unique vertices (each shared by 3 faces)
        A=(0.0,0.0,0.0); B=(1.0,0.0,0.0); C=(0.0,1.0,0.0); D=(0.0,0.0,1.0)
        faces = [(A,B,C),(A,B,D),(A,C,D),(B,C,D)]
        buf = IOBuffer()
        println(buf, "solid tet")
        for (p,q,r) in faces
            println(buf, "facet normal 0 0 0"); println(buf, "  outer loop")
            for v in (p,q,r); println(buf, "    vertex $(v[1]) $(v[2]) $(v[3])"); end
            println(buf, "  endloop"); println(buf, "endfacet")
        end
        println(buf, "endsolid tet")
        pstl = joinpath(dir, "tet.stl"); write(pstl, String(take!(buf)))
        m = read_stl(pstl)
        @test nnodes(m) == 4          # welded from 12 raw vertices
        @test ntris(m) == 4
        be, maxinc = boundary_edges(m.tris)
        @test isempty(be)             # closed surface: no boundary edges
        @test maxinc == 2             # manifold
    end

    @testset "STL ingest (binary)" begin
        A=(0.0f0,0.0f0,0.0f0); B=(1.0f0,0.0f0,0.0f0); C=(0.0f0,1.0f0,0.0f0); D=(0.0f0,0.0f0,1.0f0)
        faces = [(A,B,C),(A,B,D),(A,C,D),(B,C,D)]
        buf = IOBuffer()
        write(buf, zeros(UInt8, 80))            # header
        write(buf, UInt32(length(faces)))       # triangle count
        for (p,q,r) in faces
            write(buf, 0.0f0, 0.0f0, 0.0f0)     # normal
            for v in (p,q,r); write(buf, v[1], v[2], v[3]); end
            write(buf, UInt16(0))               # attribute
        end
        pstl = joinpath(dir, "tet_bin.stl"); write(pstl, take!(buf))
        m = read_stl(pstl)
        @test nnodes(m) == 4
        @test ntris(m) == 4
    end

    @testset "read_geo_params on the enclosure fixture" begin
        gp = read_geo_params(joinpath(@__DIR__, "fixtures", "enclosure_coax_junction.geo"))
        @test gp.mesh_size_min ≈ 0.00026669999999999933 rtol=1e-9
        @test gp.random_seed == 1
        # physical groups declared in the .geo (dim, tag) => name
        @test gp.physical_groups[(3,1)] == "air"
        @test gp.physical_groups[(3,2)] == "coax_pin"
        @test gp.physical_groups[(3,4)] == "case"
        @test gp.physical_groups[(2,3)] == "resistor"
        @test gp.physical_groups[(2,7)] == "radiation"
        @test gp.physical_groups[(1,6)] == "p1_line"
        # the three volumes the mesher must fill are all present
        vols = sort([name for ((d,t),name) in gp.physical_groups if d==3])
        @test vols == ["air", "case", "coax_pin"]
    end

    @testset "empty mesh round-trips" begin
        e = Mesh(Matrix{Float64}(undef, 3, 0))
        for ver in (2.2, 4.1)
            p = joinpath(dir, "empty_$(ver).msh")
            write_msh(p, e; version=ver)
            f = read_msh(p)
            @test nnodes(f.mesh) == 0 && ntets(f.mesh) == 0
        end
    end
end
