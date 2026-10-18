# ── .p3d (PLOT3D structured grid) reader + .mail/.ir3 (CEA/Iridium) writers ──
#
# Correctness  : hand-written `.p3d` fixtures follow the pinned
#                `GModelIO_P3D.cpp` token-stream layout (block count, `ni nj
#                nk` per block — `nk == 1` is the 2-D quad-sheet variant —
#                then three `ni*nj*nk` coordinate passes, `i` fastest), and
#                block cells/entities are checked literally. `.mail`/`.ir3`
#                writers are checked against the pinned
#                `GModelIO_MAIL.cpp`/`GModelIO_IR3.cpp` record layouts —
#                `.mail` and `.ir3` are write-only formats upstream (no
#                reader exists), and `.p3d` writes only structured
#                transfinite grids upstream, so `write_p3d` rejects like the
#                "No structured grids to save" path for every representable
#                mesh.
# Robustness   : non-positive block counts/dims, truncated token streams,
#                non-numeric data, count limits, special blocks, partition
#                tag type and session-level dispatch guards.
# Completeness : `API.open`/`API.merge` accept `.p3d`; `API.write` emits
#                `.mail`/`.ir3` (elementary tag per the upstream default)
#                and rejects `.p3d` loudly; `IO.write_ir3` honors the
#                `element_tag_type` keyword directly.

using Test
import Tessella
using Printf: @sprintf
using Tessella.Elements: MixedMesh, ElementBlock, SpecialElementBlock
using Tessella.IO

const _P3D_MAIL_API = Tessella.API

# Single 2×2×2 block: x runs fastest, then j, then k; one unit-cube hex.
const _P3D_HEX = """
1
2 2 2
0 1 0 1 0 1 0 1
0 0 1 1 0 0 1 1
0 0 0 0 1 1 1 1
"""

# Two blocks on one deck: a 3×2×1 quad sheet (entity 1 of dimension 2) and a
# 2×2×2 hex block (entity 1 of dimension 3 — counters are per-dimension).
const _P3D_MIXED = """
2
3 2 1
2 2 2
1 2 3 1 2 3
0 0 0 1 1 1
0 0 0 0 0 0
0 1 0 1 0 1 0 1
0 0 1 1 0 0 1 1
0 0 0 0 1 1 1 1
"""

@testset "P3D/MAIL/IR3 mesh I/O" begin

    @testset "P3D reader: single hex block" begin
        dir = mktempdir()
        path = joinpath(dir, "cube.p3d")
        Base.write(path, _P3D_HEX)
        f = read_p3d(path)
        @test f isa P3dFile
        @test f.node_tags == UInt64.(1:8)
        @test size(f.mesh.coords, 2) == 8
        # node row = i + (j-1)*ni + (k-1)*ni*nj
        @test f.mesh.coords[:, 1] == [0.0, 0.0, 0.0]
        @test f.mesh.coords[:, 4] == [1.0, 1.0, 0.0]
        @test f.mesh.coords[:, 8] == [1.0, 1.0, 1.0]
        @test [b.msh for b in f.mesh.blocks] == [5]
        # hex: (i,j,k),(i+1,j,k),(i+1,j+1,k),(i,j+1,k) then k+1
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 4, 3, 5, 6, 8, 7], 8, 1)
        @test f.mesh.elementary_entities == [Int32[1]]
        @test f.mesh.blocks[1].tags == Int32[0]      # no physicals in .p3d
        @test f.element_tags == [UInt64[1]]
    end

    @testset "P3D reader: mixed quad/hex blocks, per-dim entities" begin
        dir = mktempdir()
        path = joinpath(dir, "mixed.p3d")
        Base.write(path, _P3D_MIXED)
        f = read_p3d(path)
        @test size(f.mesh.coords, 2) == 14
        @test [b.msh for b in f.mesh.blocks] == [3, 5]
        # 3×2 sheet: quads [1 2 5 4] and [2 3 6 5]
        @test f.mesh.blocks[1].nodes == Int32[1 2; 2 3; 5 6; 4 5]
        @test f.mesh.blocks[2].nodes ==
            reshape(Int32[7, 8, 10, 9, 11, 12, 14, 13], 8, 1)
        @test f.mesh.elementary_entities == [Int32[1, 1], Int32[1]]
        @test f.mesh.coords[:, 6] == [3.0, 1.0, 0.0]
    end

    @testset "P3D reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".p3d"))
            Base.write(p, text)
            return p
        end
        @test_throws ArgumentError read_p3d(_w("0\n"))           # no blocks
        @test_throws ArgumentError read_p3d(_w("-2\n"))          # negative
        @test_throws ArgumentError read_p3d(_w("x\n"))           # bad int
        @test_throws ArgumentError read_p3d(_w("1\n"))           # missing dims
        @test_throws ArgumentError read_p3d(_w("1\n2 2\n"))      # short dims
        @test_throws ArgumentError read_p3d(_w("1\n2 2 x\n"))    # bad dim
        @test_throws ArgumentError read_p3d(_w("1\n0 2 2\n"))    # zero dim
        @test_throws ArgumentError read_p3d(_w("1\n2 2 -1\n"))   # negative dim
        @test_throws ArgumentError read_p3d(_w(                  # truncated x
            "1\n2 2 2\n0 1 0 1\n"))
        @test_throws ArgumentError read_p3d(_w(                  # truncated z
            "1\n2 2 2\n0 1 0 1 0 1 0 1\n0 0 1 1 0 0 1 1\n0 0\n"))
        @test_throws ArgumentError read_p3d(_w(                  # bad float
            "1\n2 2 2\n0 1 0 1 0 1 zz 1\n0 0 1 1 0 0 1 1\n" *
            "0 0 0 0 1 1 1 1\n"))
        @test_throws ArgumentError read_p3d(_w(                  # non-finite
            "1\n2 2 2\n0 1 0 1 0 1 NaN 1\n0 0 1 1 0 0 1 1\n" *
            "0 0 0 0 1 1 1 1\n"))
        @test_throws ArgumentError read_p3d(joinpath(dir, "missing.p3d"))
        @test_throws ArgumentError read_p3d(_w(_P3D_HEX); max_nodes=7)
        @test_throws ArgumentError read_p3d(_w(_P3D_HEX); max_elements=0)
        @test_throws ArgumentError read_p3d(_w(_P3D_HEX); max_file_bytes=4)
        # a 1×N degenerate block is legal: nodes but no cells
        f = read_p3d(_w("1\n1 3 1\n0 0 0 1 1 1 0 0 0\n"))
        @test size(f.mesh.coords, 2) == 3
        @test isempty(f.mesh.blocks)
        # trailing tokens after the declared data are ignored like upstream
        f = read_p3d(_w(_P3D_HEX * "9 9 9\n"))
        @test size(f.mesh.coords, 2) == 8
    end

    @testset "P3D writer: no structured grids to save" begin
        dir = mktempdir()
        hex = MixedMesh(zeros(3, 8),
            [ElementBlock(5, reshape(Int32[1:8;], 8, 1), Int32[1])])
        out = joinpath(dir, "x.p3d")
        @test_throws ArgumentError write_p3d(out, hex)
        @test !isfile(out)                       # no partial file left behind
        @test_throws ArgumentError write_p3d(out, 42)
        @test_throws ArgumentError write_p3d(joinpath(dir, "nodir", "x.p3d"),
            hex)
    end

    @testset "MAIL writer: CEA triangulation deck" begin
        dir = mktempdir()
        coords = Float64[0 1 1 0; 0 0 1 0; 0 0 0 2.5]
        blocks = ElementBlock[
            ElementBlock(2, Int32[1 2; 2 4; 3 3], Int32[7, 7]),   # 2 tris
            ElementBlock(3, reshape(Int32[1, 2, 4, 3], 4, 1), Int32[0]), # quad
            ElementBlock(4, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0]), # tet
            ElementBlock(1, reshape(Int32[1, 4], 2, 1), Int32[0]),       # line
        ]
        mesh = MixedMesh(coords, blocks)
        path = joinpath(dir, "tri.mail")
        @test write_mail(path, mesh) == path
        lines = readlines(path)
        @test lines[1] == " 4 2"
        @test length(lines) == 1 + 4 + 2 + 2   # header + nodes + tris + edges
        @test all(l -> startswith(l, " "), lines)
        @test lines[2] == @sprintf(" %19.10E %19.10E %19.10E", 0.0, 0.0, 0.0)
        @test lines[5] == @sprintf(" %19.10E %19.10E %19.10E", 0.0, 0.0, 2.5)
        @test lines[6] == " 1 2 3"
        @test lines[7] == " 2 4 3"
        @test lines[8:9] == [" 0 0 0", " 0 0 0"]
        # higher-order triangles emit their three corner vertices
        tri6 = MixedMesh(zeros(3, 6),
            [ElementBlock(9, reshape(Int32[1:6;], 6, 1), Int32[1])])
        p6 = joinpath(dir, "tri6.mail")
        write_mail(p6, tri6)
        l6 = readlines(p6)
        @test l6[1] == " 6 1"
        @test l6[8] == " 1 2 3"
        # no triangles → an empty-but-valid deck
        pe = joinpath(dir, "empty.mail")
        write_mail(pe, MixedMesh(zeros(3, 2), ElementBlock[]))
        @test readlines(pe)[1] == " 2 0"
        special = MixedMesh(zeros(3, 2),
            [SpecialElementBlock(67, reshape(Int32[1, 2], 2, 1))])
        @test_throws ArgumentError write_mail(joinpath(dir, "s.mail"), special)
        @test_throws ArgumentError write_mail(
            joinpath(dir, "nodir", "x.mail"), tri6)
    end

    @testset "IR3 writer: record layout and tag types" begin
        dir = mktempdir()
        coords = Float64[0 1 1 0 0.5; 0 0 1 0 0.5; 0 0 0 1 0.25]
        blocks = ElementBlock[
            ElementBlock(1, reshape(Int32[1, 5], 2, 1), Int32[0]),   # skipped
            ElementBlock(2, Int32[1 2; 2 4; 3 3], Int32[7, 9]),      # 2 tris
            ElementBlock(3, reshape(Int32[1, 2, 4, 3], 4, 1), Int32[7]),
            ElementBlock(4, reshape(Int32[1, 2, 3, 5], 4, 1), Int32[4]),
            ElementBlock(15, reshape(Int32[5], 1, 1), Int32[0]),     # skipped
        ]
        mesh = MixedMesh(coords, blocks;
            elementary_entities=[Int32[1], Int32[11, 11], Int32[12], Int32[21],
                                 Int32[31]])
        path = joinpath(dir, "m.ir3")
        @test write_ir3(path, mesh) == path
        lines = readlines(path)
        @test lines[1] == "33"
        @test lines[2] == "5 3 1"                    # nv num2D num3D
        @test length(lines) == 2 + 5 + 3 + 1
        @test lines[3] == "1 0 0 0"
        @test lines[7] == "5 0.5 0.5 0.25"
        # dim-2 records numbered 1..3 with elementary tags, then dim-3 restarts
        @test lines[8] == "1 11 3 1 2 3"
        @test lines[9] == "2 11 3 2 4 3"
        @test lines[10] == "3 12 4 1 2 4 3"
        @test lines[11] == "1 21 4 1 2 3 5"
        # Mesh.SaveElementTagType=2 → abs(physical)
        p2 = joinpath(dir, "phys.ir3")
        write_ir3(p2, mesh; element_tag_type=2)
        l2 = readlines(p2)
        @test l2[8] == "1 7 3 1 2 3"
        @test l2[9] == "2 9 3 2 4 3"
        @test l2[11] == "1 4 4 1 2 3 5"
        # (a negative physical tag would reverse connectivity upstream, but
        # `ElementBlock` validates tags non-negative so the path is
        # unreachable through construction — kept as defensive parity.)
        # partition tags (type 3) are not representable; other values fall
        # back to elementary like upstream
        @test_throws ArgumentError write_ir3(joinpath(dir, "p3.ir3"), mesh;
            element_tag_type=3)
        p9 = joinpath(dir, "nine.ir3")
        write_ir3(p9, mesh; element_tag_type=9)
        @test readlines(p9)[8] == "1 11 3 1 2 3"
        # single-dimension header form: only dim-2 cells
        tris = MixedMesh(coords[:, 1:3],
            [ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[1])])
        pt = joinpath(dir, "t.ir3")
        write_ir3(pt, tris)
        @test readlines(pt)[2] == "3 1"
        # no dim-2/3 cells → `nv 0`
        pe = joinpath(dir, "e.ir3")
        write_ir3(pe, MixedMesh(zeros(3, 2), ElementBlock[]))
        @test readlines(pe)[2] == "2 0"
        special = MixedMesh(zeros(3, 2),
            [SpecialElementBlock(67, reshape(Int32[1, 2], 2, 1))])
        @test_throws ArgumentError write_ir3(joinpath(dir, "s.ir3"), special)
        @test_throws ArgumentError write_ir3(
            joinpath(dir, "nodir", "x.ir3"), tris)
    end

    @testset "session open/merge/write dispatch" begin
        _P3D_MAIL_API.initialize()
        try
            _P3D_MAIL_API.option("General.Terminal", 0)
            mktempdir() do dir
                p3d = joinpath(dir, "cube.p3d")
                Base.write(p3d, _P3D_MIXED)
                _P3D_MAIL_API.model.add("p3dhost")
                @test _P3D_MAIL_API.open(p3d) === nothing
                entities = _P3D_MAIL_API.model.get_entities()
                @test (2, 1) in entities && (3, 1) in entities
                types = _P3D_MAIL_API.mesh.get_elements()[1]
                @test Int32(3) in types && Int32(5) in types
                # merge keeps the host model and adds another block's cells
                host = copy(entities)
                @test _P3D_MAIL_API.merge(p3d) === nothing
                merged = _P3D_MAIL_API.model.get_entities()
                @test issubset(host, merged)
                @test length(_P3D_MAIL_API.mesh.get_elements()[1]) >=
                    length(types)
                # .mail/.ir3 export through API.write; .p3d rejects
                mail = joinpath(dir, "out.mail")
                ir3 = joinpath(dir, "out.ir3")
                @test _P3D_MAIL_API.write(mail) === nothing
                @test isfile(mail)
                @test _P3D_MAIL_API.write(ir3) === nothing
                @test isfile(ir3)
                @test readlines(ir3)[1] == "33"
                @test_throws ArgumentError _P3D_MAIL_API.write(
                    joinpath(dir, "out.p3d"))
                @test !isfile(joinpath(dir, "out.p3d"))
                # .mail/.ir3 are write-only formats upstream
                @test_throws ArgumentError _P3D_MAIL_API.open(mail)
                @test_throws ArgumentError _P3D_MAIL_API.merge(ir3)
            end
        finally
            _P3D_MAIL_API.finalize()
        end
    end
end
