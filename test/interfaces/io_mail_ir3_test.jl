# ── .mail (CEA triangulation) and .ir3 (Iridium) readers ──
#
# Correctness  : `.mail`/`.ir3` are write-only formats upstream, so the
#                readers are literal inverses of `write_mail`/`write_ir3`.
#                Round-trip fixtures check node rows, connectivity, per-cell
#                entity tags and the `nv nt`/`33` headers; hand-written decks
#                cover the tolerated edge-block tail, sparse vertex indices,
#                the per-section dimension split and the single-count
#                header's arity/planarity dimension resolution.
# Robustness   : bad format lines/counts, truncated token streams, dangling
#                or duplicate node references, negative tags, unsupported
#                vertex counts and the count limits.
# Completeness : `API.open`/`API.merge` route `.mail`/`.ir3` through the
#                same foreign-mesh staging channel as `.unv`/`.p3d`.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock
using Tessella.IO

const _MAIL_IR3_API = Tessella.API

@testset "MAIL/IR3 mesh readers" begin

    @testset "MAIL reader: writer round-trip" begin
        dir = mktempdir()
        coords = Float64[0 1 1 0.5; 0 0 1 0.5; 0 0 0 2.5]
        mesh = MixedMesh(coords,
            [ElementBlock(2, Int32[1 2; 2 4; 3 3], Int32[0, 0])])
        path = joinpath(dir, "t.mail")
        @test write_mail(path, mesh) == path
        f = read_mail(path)
        @test f isa MailFile
        @test f.node_tags == UInt64.(1:4)
        @test f.mesh.coords == coords
        @test [b.msh for b in f.mesh.blocks] == [2]
        @test f.mesh.blocks[1].nodes == Int32[1 2; 2 4; 3 3]
        @test f.mesh.elementary_entities == [Int32[1, 1]]   # implicit surface
        @test f.mesh.blocks[1].tags == Int32[0, 0]
        @test f.element_tags == [UInt64[1, 2]]
        # quads/tets/lines are skipped by the writer like upstream
        mixed = MixedMesh(coords,
            [ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0]),
             ElementBlock(4, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0]),
             ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[0]),
             ElementBlock(1, reshape(Int32[1, 4], 2, 1), Int32[0])])
        pm = joinpath(dir, "m.mail")
        write_mail(pm, mixed)
        fm = read_mail(pm)
        @test [b.msh for b in fm.mesh.blocks] == [2]
        @test size(fm.mesh.blocks[1].nodes, 2) == 1
        # a deck that ends cleanly after the triangle block is tolerated
        pe = joinpath(dir, "noedges.mail")
        write_mail(pe, mesh)
        lines = readlines(pe)
        Base.write(pe, join(lines[1:end-2], "\n") * "\n")     # drop edge rows
        fe = read_mail(pe)
        @test f.mesh.coords == fe.mesh.coords
        @test fe.mesh.blocks[1].nodes == f.mesh.blocks[1].nodes
        # higher-order triangles read back their three corner vertices
        p6 = joinpath(dir, "tri6.mail")
        write_mail(p6, MixedMesh(zeros(3, 6),
            [ElementBlock(9, reshape(Int32[1:6;], 6, 1), Int32[1])]))
        f6 = read_mail(p6)
        @test [b.msh for b in f6.mesh.blocks] == [2]
        @test f6.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3], 3, 1)
        @test f6.mesh.elementary_entities == [Int32[1]]
    end

    @testset "MAIL reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".mail"))
            Base.write(p, text)
            return p
        end
        @test_throws ArgumentError read_mail(_w("-1 2\n"))        # negative nv
        @test_throws ArgumentError read_mail(_w("2 -1\n"))        # negative nt
        @test_throws ArgumentError read_mail(_w("x 2\n"))         # bad int
        @test_throws ArgumentError read_mail(_w("2 1\n"))         # no vertices
        @test_throws ArgumentError read_mail(_w(
            "1 1\n0 0 zz\n1 1 1\n"))                            # bad float
        @test_throws ArgumentError read_mail(_w(
            "1 1\n0 0 NaN\n1 1 1\n"))                           # non-finite
        @test_throws ArgumentError read_mail(_w(
            "3 1\n0 0 0\n1 0 0\n0 1 0\n1 2\n"))                 # short record
        @test_throws ArgumentError read_mail(_w(
            "3 1\n0 0 0\n1 0 0\n0 1 0\n1 2 x\n"))               # bad index
        @test_throws ArgumentError read_mail(_w(
            "3 1\n0 0 0\n1 0 0\n0 1 0\n1 2 0\n"))               # index 0
        @test_throws ArgumentError read_mail(_w(
            "3 1\n0 0 0\n1 0 0\n0 1 0\n1 2 4\n"))               # dangling
        @test_throws ArgumentError read_mail(_w(
            "3 1\n0 0 0\n1 0 0\n0 1 0\n1 2 3\n0 0\n"))          # half edge
        @test_throws ArgumentError read_mail(_w(
            "3 1\n0 0 0\n1 0 0\n0 1 0\n1 2 3\n0 x 0\n"))        # bad edge
        @test_throws ArgumentError read_mail(joinpath(dir, "missing.mail"))
        good = _w("4 2\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n" *
                  "1 2 3\n2 3 4\n0 0 0\n0 0 0\n")
        @test_throws ArgumentError read_mail(good; max_nodes=3)
        @test_throws ArgumentError read_mail(good; max_elements=1)
        @test_throws ArgumentError read_mail(good; max_file_bytes=4)
        f = read_mail(good)
        @test size(f.mesh.coords, 2) == 4
        @test size(f.mesh.blocks[1].nodes, 2) == 2
        # an empty-but-valid deck
        fe = read_mail(_w("2 0\n0 0 0\n1 0 0\n"))
        @test size(fe.mesh.coords, 2) == 2
        @test isempty(fe.mesh.blocks)
    end

    @testset "IR3 reader: three-count header round-trip" begin
        dir = mktempdir()
        coords = Float64[0 1 1 0 0.5; 0 0 1 0 0.5; 0 0 0 1 0.25]
        mesh = MixedMesh(coords,
            [ElementBlock(1, reshape(Int32[1, 5], 2, 1), Int32[0]),   # skipped
             ElementBlock(2, Int32[1 2; 2 4; 3 3], Int32[7, 9]),
             ElementBlock(3, reshape(Int32[1, 2, 4, 3], 4, 1), Int32[7]),
             ElementBlock(4, reshape(Int32[1, 2, 3, 5], 4, 1), Int32[4]),
             ElementBlock(15, reshape(Int32[5], 1, 1), Int32[0])];  # skipped
            elementary_entities=[Int32[1], Int32[11, 11], Int32[12],
                                 Int32[21], Int32[31]])
        path = joinpath(dir, "m.ir3")
        @test write_ir3(path, mesh) == path
        @test readlines(path)[2] == "5 3 1"
        f = read_ir3(path)
        @test f isa Ir3File
        @test f.node_tags == UInt64.(1:5)
        @test f.mesh.coords == coords
        # dim-2 records keep file order, then dim-3 restarts its own section
        @test [b.msh for b in f.mesh.blocks] == [2, 3, 4]
        @test f.mesh.blocks[1].nodes == Int32[1 2; 2 4; 3 3]
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 2, 4, 3], 4, 1)
        @test f.mesh.blocks[3].nodes == reshape(Int32[1, 2, 3, 5], 4, 1)
        # the tag column round-trips as the elementary entity tag
        @test f.mesh.elementary_entities ==
            [Int32[11, 11], Int32[12], Int32[21]]
        @test all(b -> all(==(0), b.tags), f.mesh.blocks)
        # element tags are dense file-order counters across the dim split
        @test f.element_tags == [UInt64[1, 2], UInt64[3], UInt64[4]]
        # same arity in both sections decodes per section: tri6 vs pri6
        ho = MixedMesh(zeros(3, 9),
            [ElementBlock(9, reshape(Int32[1:6;], 6, 1), Int32[2]),
             ElementBlock(6, reshape(Int32[1, 2, 3, 7, 8, 9], 6, 1), Int32[3])];
            elementary_entities=[Int32[5], Int32[6]])
        ph = joinpath(dir, "ho.ir3")
        write_ir3(ph, ho)
        @test readlines(ph)[2] == "9 1 1"
        fh = read_ir3(ph)
        @test [b.msh for b in fh.mesh.blocks] == [9, 6]
        @test fh.mesh.blocks[1].nodes == reshape(Int32[1:6;], 6, 1)
        @test fh.mesh.blocks[2].nodes == reshape(Int32[1, 2, 3, 7, 8, 9], 6, 1)
        @test fh.mesh.elementary_entities == [Int32[5], Int32[6]]
    end

    @testset "IR3 reader: single-count header dimension resolution" begin
        dir = mktempdir()
        # tri3 arity exists only in dim 2 — resolved without geometry
        tris = MixedMesh(Float64[0 1 0; 0 0 1; 0 0 1],
            [ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[4])];
            elementary_entities=[Int32[8]])
        pt = joinpath(dir, "t.ir3")
        write_ir3(pt, tris)
        @test readlines(pt)[2] == "3 1"
        ft = read_ir3(pt)
        @test [b.msh for b in ft.mesh.blocks] == [2]
        @test ft.mesh.elementary_entities == [Int32[8]]
        # pyr5 arity exists only in dim 3
        pyrs = MixedMesh(Float64[0 1 1 0 0.5; 0 0 1 1 0.5; 0 0 0 0 1],
            [ElementBlock(7, reshape(Int32[1:5;], 5, 1), Int32[0])])
        pp = joinpath(dir, "p.ir3")
        write_ir3(pp, pyrs)
        @test readlines(pp)[2] == "5 1"
        fp = read_ir3(pp)
        @test [b.msh for b in fp.mesh.blocks] == [7]
        @test fp.mesh.blocks[1].nodes == reshape(Int32[1:5;], 5, 1)
        # arity-4 is ambiguous (qua4/tet4): the non-planar tet reads dim-3
        tets = MixedMesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1],
            [ElementBlock(4, reshape(Int32[1:4;], 4, 1), Int32[0])])
        pv = joinpath(dir, "v.ir3")
        write_ir3(pv, tets)
        @test readlines(pv)[2] == "4 1"
        fv = read_ir3(pv)
        @test [b.msh for b in fv.mesh.blocks] == [4]
        @test fv.mesh.blocks[1].nodes == reshape(Int32[1:4;], 4, 1)
        # the same record over a planar node set reads dim-2
        quads = MixedMesh(Float64[0 1 1 0; 0 0 1 1; 0 0 0 0],
            [ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0])])
        pq = joinpath(dir, "q.ir3")
        write_ir3(pq, quads)
        fq = read_ir3(pq)
        @test [b.msh for b in fq.mesh.blocks] == [3]
        @test fq.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        # sparse vertex indices become the file node labels
        sparse = joinpath(dir, "s.ir3")
        Base.write(sparse, "33\n3 1\n10 0 0 0\n20 1 0 0\n30 0 1 0\n" *
                           "1 7 3 10 20 30\n")
        fs = read_ir3(sparse)
        @test fs.node_tags == UInt64[10, 20, 30]
        @test fs.mesh.coords == Float64[0 1 0; 0 0 1; 0 0 0]
        @test fs.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3], 3, 1)
        @test fs.mesh.elementary_entities == [Int32[7]]
        # no dim-2/3 cells → header `nv 0`, an empty block list
        pe = joinpath(dir, "e.ir3")
        write_ir3(pe, MixedMesh(zeros(3, 2), ElementBlock[]))
        @test readlines(pe)[2] == "2 0"
        fe = read_ir3(pe)
        @test size(fe.mesh.coords, 2) == 2
        @test isempty(fe.mesh.blocks)
    end

    @testset "IR3 reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".ir3"))
            Base.write(p, text)
            return p
        end
        @test_throws ArgumentError read_ir3(_w(""))               # empty
        @test_throws ArgumentError read_ir3(_w("34\n3 1\n"))      # bad format
        @test_throws ArgumentError read_ir3(_w("33\n"))           # no header
        @test_throws ArgumentError read_ir3(_w("33\n4\n"))        # one int
        @test_throws ArgumentError read_ir3(_w("33\n4 1 0 9\n"))  # four ints
        @test_throws ArgumentError read_ir3(_w("33\n-1 0\n"))     # negative nv
        @test_throws ArgumentError read_ir3(_w("33\n3 -1\n"))     # negative ne
        @test_throws ArgumentError read_ir3(_w("33\n3 1 -1\n"))   # negative n3
        @test_throws ArgumentError read_ir3(_w("33\nx 1\n"))      # bad int
        @test_throws ArgumentError read_ir3(_w(
            "33\n1 0\n"))                                       # missing row
        @test_throws ArgumentError read_ir3(_w(
            "33\n1 0\n0 0 0 0\n"))                              # index 0
        @test_throws ArgumentError read_ir3(_w(
            "33\n2 0\n1 0 0 0\n1 1 0 0\n"))                     # dup index
        @test_throws ArgumentError read_ir3(_w(
            "33\n1 0\n1 0 0 zz\n"))                             # bad float
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "0 0 3 1 2 3\n"))                                   # element num 0
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "1 -2 3 1 2 3\n"))                                  # negative tag
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "1 0 0\n"))                                         # zero count
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "1 0 2 1 2\n"))                                     # arity 2
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "1 0 3 1 2 9\n"))                                   # dangling node
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "1 0 3 1 2\n"))                                     # short record
        # arity 3 has no dim-3 type even inside a declared 3-int section
        @test_throws ArgumentError read_ir3(_w(
            "33\n3 0 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
            "1 0 3 1 2 3\n"))
        # a dim-3 record inside a declared dim-2 section fails the same way
        @test_throws ArgumentError read_ir3(_w(
            "33\n4 1 0\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n" *
            "1 0 5 1 2 3 4 1\n"))
        @test_throws ArgumentError read_ir3(joinpath(dir, "missing.ir3"))
        good = _w("33\n4 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n" *
                  "1 0 4 1 2 3 4\n")
        @test_throws ArgumentError read_ir3(good; max_nodes=3)
        @test_throws ArgumentError read_ir3(good; max_elements=0)
        @test_throws ArgumentError read_ir3(good; max_file_bytes=4)
        f = read_ir3(good)
        @test [b.msh for b in f.mesh.blocks] == [4]   # non-planar → tet
        # trailing tokens after the declared records are ignored like .p3d
        ft = read_ir3(_w("33\n3 1\n1 0 0 0\n2 1 0 0\n3 0 1 0\n" *
                         "1 0 3 1 2 3\n9 9 9\n"))
        @test [b.msh for b in ft.mesh.blocks] == [2]
    end

    @testset "session open/merge dispatch" begin
        _MAIL_IR3_API.initialize()
        try
            _MAIL_IR3_API.option("General.Terminal", 0)
            mktempdir() do dir
                mail = joinpath(dir, "tri.mail")
                write_mail(mail, MixedMesh(Float64[0 1 0; 0 0 1; 0 0 0],
                    [ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[0])]))
                _MAIL_IR3_API.model.add("mailhost")
                @test _MAIL_IR3_API.open(mail) === nothing
                entities = _MAIL_IR3_API.model.get_entities()
                @test (2, 1) in entities                      # implicit surface
                @test Int32(2) in _MAIL_IR3_API.mesh.get_elements()[1]
                # merge keeps the host and adds the incoming cells
                host = copy(entities)
                @test _MAIL_IR3_API.merge(mail) === nothing
                @test issubset(host, _MAIL_IR3_API.model.get_entities())
                # an .ir3 deck publishes its records on the tag-column tags
                ir3 = joinpath(dir, "m.ir3")
                write_ir3(ir3, MixedMesh(
                    Float64[0 1 1 0 0.5; 0 0 1 0 0.5; 0 0 0 1 0.25],
                    [ElementBlock(2, Int32[1 2; 2 4; 3 3], Int32[0, 0]),
                     ElementBlock(4, reshape(Int32[1, 2, 3, 5], 4, 1),
                                  Int32[0])];
                    elementary_entities=[Int32[11, 11], Int32[21]]))
                @test _MAIL_IR3_API.open(ir3) === nothing
                ents = _MAIL_IR3_API.model.get_entities()
                @test (2, 11) in ents && (3, 21) in ents
                types = _MAIL_IR3_API.mesh.get_elements()[1]
                @test Int32(2) in types && Int32(4) in types
            end
        finally
            _MAIL_IR3_API.finalize()
        end
    end
end
