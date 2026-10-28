# ── .su2 (SU2 native) and .bdf/.nas (Nastran bulk data) mesh I/O ─────────────
#
# Correctness  : hand-written fixtures transcribe the documented file layouts
#                (SU2 NDIME/NPOIN/NELEM/NMARK keyword sections; Nastran GRID
#                and element cards in fixed 8-column, comma-delimited free and
#                +-continuation forms). A real box mesh generated through
#                `API.mesh.generate` is exported with `API.write` and
#                re-imported with `API.open` for all three extensions,
#                checking node/element/physical parity where representable.
# Robustness   : `D`/`d` and implied exponents, `%`/`$` comments, mixed
#                index/non-index point rows, dangling node references,
#                duplicate labels, multi-zone and INCLUDE rejection, skipped
#                `cp` frames, malformed cards and unrepresentable element
#                types on write.
# Completeness : SU2 markers ↔ physical groups+names, BDF property ids ↔
#                physical tags; `API.merge` renumbers file-local labels into
#                fresh session UIDs and preserves the existing model.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock, SpecialElementBlock,
    msh_num_nodes
using Tessella.IO

const _SU2_BDF_API=Tessella.API

# 3-D SU2 fixture: 0-based indexed point rows, a `D` exponent, NELEM records
# (tet + pyramid) and two named boundary markers (triangles + a quad).
const _SU2_FIXTURE = """
% a comment line is skipped
NDIME= 3
NPOIN= 5
0.0 0.0 0.0 0
1.0D+00 0.0 0.0 1
0.0 1.0 0.0 2
0.0 0.0 1.0 3
0.5 0.5 0.5 4
NELEM= 2
10 0 1 2 3 0
14 0 1 2 4 3 1
NMARK= 2
MARKER_TAG= wall
MARKER_ELEMS= 2
5 0 1 2
5 0 1 3
MARKER_TAG= inlet
MARKER_ELEMS= 1
9 0 2 3 1
"""

# 2-D SU2 fixture without the trailing point-index column; markers are lines.
const _SU2_2D_FIXTURE = """
NDIME= 2
NPOIN= 3
0.0 0.0
1.0 0.0
0.0 1.0
NELEM= 1
5 0 1 2 0
NMARK= 1
MARKER_TAG= edge
MARKER_ELEMS= 2
3 0 1
3 1 2
"""

# Right-justified 8-column Nastran field helpers for hand-built cards.
_bdf8(fields...) = join(lpad(f, 8) for f in fields)

# Sorted connectivity sets per MSH family, keyed by node UID (as in the
# .unv/.mesh suite) plus a coordinate-keyed variant for formats whose node
# numbering is positional rather than UID-preserving.
function _cells_by_type(get_elements)
    types, _, nodes = get_elements
    out = Dict{Int32, Vector{Vector{UInt64}}}()
    for i in eachindex(types)
        arity = msh_num_nodes(types[i])
        flat = nodes[i]
        out[types[i]] = sort!([sort!(flat[(k-1)*arity+1:k*arity])
                               for k in 1:length(flat)÷arity])
    end
    return out
end

function _cells_coord_sets(get_elements, uids, coords)
    types, _, conn = get_elements
    pos = Dict{UInt64, Int}(u => i for (i, u) in enumerate(uids))
    out = Dict{Int32, Vector{Vector{NTuple{3, Float64}}}}()
    for i in eachindex(types)
        arity = msh_num_nodes(types[i])
        flat = conn[i]
        cells = [flat[(k-1)*arity+1:k*arity] for k in 1:length(flat)÷arity]
        out[types[i]] = sort!([sort!(
            [NTuple{3, Float64}((coords[3pos[u]-2], coords[3pos[u]-1],
                                coords[3pos[u]])) for u in c])
            for c in cells])
    end
    return out
end

@testset "SU2/Nastran mesh I/O" begin

    @testset "SU2 reader: sections, markers, names" begin
        dir = mktempdir()
        path = joinpath(dir, "fixture.su2")
        Base.write(path, _SU2_FIXTURE)
        f = read_su2(path)
        @test f isa Su2File
        @test f.node_tags == UInt64[1, 2, 3, 4, 5]
        @test size(f.mesh.coords, 2) == 5
        @test f.mesh.coords[:, 2] == [1.0, 0.0, 0.0]   # 1.0D+00
        @test [b.msh for b in f.mesh.blocks] == [4, 7, 2, 3]
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 2, 3, 5, 4], 5, 1)
        @test f.mesh.blocks[3].nodes == Int32[1 1; 2 2; 3 4]
        @test f.mesh.blocks[4].nodes == reshape(Int32[1, 3, 4, 2], 4, 1)
        # NELEM cells are unmarked; marker i carries physical/entity i.
        @test [b.tags for b in f.mesh.blocks] ==
              [Int32[0], Int32[0], Int32[1, 1], Int32[2]]
        @test f.mesh.elementary_entities ==
              [Int32[1], Int32[1], Int32[1, 1], Int32[2]]
        @test f.mesh.physical_names[(2, 1)] == "wall"
        @test f.mesh.physical_names[(2, 2)] == "inlet"
        @test f.element_tags == [UInt64[1], UInt64[2], UInt64[3, 4], UInt64[5]]
    end

    @testset "SU2 reader: 2-D, non-indexed points" begin
        dir = mktempdir()
        path = joinpath(dir, "flat.su2")
        Base.write(path, _SU2_2D_FIXTURE)
        f = read_su2(path)
        @test f.node_tags == UInt64[1, 2, 3]
        @test f.mesh.coords[:, 3] == [0.0, 1.0, 0.0]
        @test [b.msh for b in f.mesh.blocks] == [2, 1]
        @test [b.tags for b in f.mesh.blocks] == [Int32[0], Int32[1, 1]]
        @test f.mesh.physical_names[(1, 1)] == "edge"
    end

    @testset "SU2 writer round-trip preserves markers and names" begin
        dir = mktempdir()
        path = joinpath(dir, "fixture.su2")
        Base.write(path, _SU2_FIXTURE)
        f = read_su2(path)
        rt = joinpath(dir, "rt.su2")
        @test write_su2(rt, f.mesh) == rt
        text = read(rt, String)
        @test occursin("NDIME= 3", text)
        @test occursin("MARKER_TAG= wall", text)
        @test occursin("MARKER_TAG= inlet", text)
        f2 = read_su2(rt)
        @test f2.mesh.coords == f.mesh.coords
        @test [b.msh for b in f2.mesh.blocks] == [4, 7, 2, 3]
        @test [b.nodes for b in f2.mesh.blocks] ==
              [b.nodes for b in f.mesh.blocks]
        @test [b.tags for b in f2.mesh.blocks] ==
              [b.tags for b in f.mesh.blocks]
        @test f2.mesh.physical_names == f.mesh.physical_names
    end

    @testset "SU2 reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".su2"))
            Base.write(p, text)
            return p
        end
        @test_throws ArgumentError read_su2(_w(           # multi-zone
            "NZONE= 2\nNDIME= 3\n"))
        @test_throws ArgumentError read_su2(_w(           # bad NDIME
            "NDIME= 4\n"))
        @test_throws ArgumentError read_su2(_w(           # NPOIN before NDIME
            "NPOIN= 1\n0.0 0.0 0.0 0\n"))
        @test_throws ArgumentError read_su2(_w(           # truncated NPOIN
            "NDIME= 3\nNPOIN= 2\n0.0 0.0 0.0 0\n"))
        @test_throws ArgumentError read_su2(_w(           # mixed index columns
            "NDIME= 3\nNPOIN= 2\n0.0 0.0 0.0 0\n1.0 0.0 0.0\n"))
        @test_throws ArgumentError read_su2(_w(           # duplicate index
            "NDIME= 3\nNPOIN= 2\n0.0 0.0 0.0 0\n1.0 0.0 0.0 0\n"))
        @test_throws ArgumentError read_su2(_w(           # unsupported type
            "NDIME= 3\nNPOIN= 1\n0.0 0.0 0.0 0\nNELEM= 1\n4 0 0\n"))
        @test_throws ArgumentError read_su2(_w(           # dangling point
            "NDIME= 3\nNPOIN= 1\n0.0 0.0 0.0 0\nNELEM= 1\n10 0 1 2 3 0\n"))
        @test_throws ArgumentError read_su2(_w(           # wrong arity
            "NDIME= 3\nNPOIN= 4\n" *
            "0.0 0.0 0.0 0\n1.0 0.0 0.0 1\n0.0 1.0 0.0 2\n0.0 0.0 1.0 3\n" *
            "NELEM= 1\n10 0 1 2\n"))
        @test_throws ArgumentError read_su2(_w(           # marker before NMARK
            "NDIME= 3\nNPOIN= 1\n0.0 0.0 0.0 0\nMARKER_TAG= x\n"))
        @test_throws ArgumentError read_su2(_w(           # missing MARKER_ELEMS
            "NDIME= 3\nNPOIN= 1\n0.0 0.0 0.0 0\nNMARK= 1\n" *
            "MARKER_TAG= x\nNELEM= 0\n"))
        @test_throws ArgumentError read_su2(_w(           # NMARK count mismatch
            "NDIME= 3\nNPOIN= 1\n0.0 0.0 0.0 0\nNMARK= 2\n" *
            "MARKER_TAG= x\nMARKER_ELEMS= 0\n"))
        @test_throws ArgumentError read_su2(joinpath(dir, "missing.su2"))
        @test_throws ArgumentError read_su2(_w(
            "NDIME= 3\nNPOIN= 1\n0.0 0.0 0.0 0\n"); max_nodes=0)
    end

    @testset "BDF reader: fixed fields, all families" begin
        dir = mktempdir()
        grids = [_bdf8("GRID", string(i), "", string(1.0i), "0.", "0.")
                 for i in 1:20]
        cards = [
            "\$ a comment",
            "BEGIN BULK",
            grids...,
            _bdf8("CBAR", "1", "5", "1", "2"),
            _bdf8("CROD", "2", "5", "2", "3"),
            _bdf8("CTUBE", "3", "5", "3", "4"),
            _bdf8("CBEAM", "4", "5", "4", "5"),
            _bdf8("CTRIA3", "5", "7", "1", "2", "3"),
            # CTRIA6 in the MSC layout: corners, blank THETA/ZOFFS/TFLAG,
            # continuation T1..T3 blanks then midside nodes.
            _bdf8("CTRIA6", "6", "7", "1", "2", "3", "", "", "", "+") *
                "\n" * _bdf8("+A", "", "", "", "4", "5", "6"),
            _bdf8("CQUAD4", "7", "8", "1", "2", "3", "4"),
            _bdf8("CQUAD8", "8", "8", "1", "2", "3", "4", "5", "6", "+") *
                "\n" * _bdf8("+B", "7", "8"),
            _bdf8("CTETRA", "9", "9", "1", "2", "3", "4"),
            # CTETRA10 in the MSC layout: G1..G4, two blanks, cont id,
            # continuation G5..G10.
            _bdf8("CTETRA", "10", "9", "1", "2", "3", "4", "", "", "+") *
                "\n" * _bdf8("+C", "5", "6", "7", "8", "9", "10"),
            _bdf8("CPENTA", "11", "10", "1", "2", "3", "4", "5", "6"),
            _bdf8("CPENTA15", "12", "10", "1", "2", "3", "4", "5", "6",
                  "+") * "\n" *
                _bdf8("+D", "7", "8", "9", "10", "11", "12", "13", "14") *
                "\n" * _bdf8("+E", "15"),
            _bdf8("CHEXA", "13", "11", "1", "2", "3", "4", "5", "6", "+") *
                "\n" * _bdf8("+F", "7", "8"),
            _bdf8("CHEXA20", "14", "11", "1", "2", "3", "4", "5", "6",
                  "+") * "\n" *
                _bdf8("+G", "7", "8", "9", "10", "11", "12", "13", "14") *
                "\n" * _bdf8("+H", "15", "16", "17", "18", "19", "20"),
            _bdf8("CPYRA", "15", "12", "1", "2", "3", "4", "5"),
            _bdf8("CPYRA13", "16", "12", "1", "2", "3", "4", "5", "6",
                  "+") * "\n" *
                _bdf8("+I", "7", "8", "9", "10", "11", "12", "13"),
            "PSHELL     1       1       0.1                     1",   # skipped
            "CORD2R        1       0       0.      0.      0.",        # skipped
            "ENDDATA",
            "GRID   900               9.      9.      9.",           # ignored
        ]
        path = joinpath(dir, "fixture.bdf")
        Base.write(path, join(cards, "\n") * "\n")
        f = read_bdf(path)
        @test f isa BdfFile
        @test f.node_tags == UInt64.(1:20)
        @test f.mesh.coords[:, 7] == [7.0, 0.0, 0.0]
        @test [b.msh for b in f.mesh.blocks] ==
              [1, 2, 9, 3, 16, 4, 11, 6, 18, 5, 17, 7, 19]
        @test f.mesh.blocks[1].nodes == Int32[1 2 3 4; 2 3 4 5]
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 2, 3], 3, 1)
        @test f.mesh.blocks[3].nodes == reshape(Int32[1:6;], 6, 1)
        @test f.mesh.blocks[4].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        @test f.mesh.blocks[5].nodes == reshape(Int32[1:8;], 8, 1)
        @test f.mesh.blocks[7].nodes == reshape(Int32[1:10;], 10, 1)
        @test f.mesh.blocks[8].nodes == reshape(Int32[1:6;], 6, 1)
        @test f.mesh.blocks[9].nodes == reshape(Int32[1:15;], 15, 1)
        @test f.mesh.blocks[10].nodes == reshape(Int32[1:8;], 8, 1)
        @test f.mesh.blocks[11].nodes == reshape(Int32[1:20;], 20, 1)
        @test f.mesh.blocks[12].nodes == reshape(Int32[1:5;], 5, 1)
        @test f.mesh.blocks[13].nodes == reshape(Int32[1:13;], 13, 1)
        @test [b.tags for b in f.mesh.blocks][1] == Int32[5, 5, 5, 5]
        @test f.mesh.elementary_entities[2] == Int32[7]
        @test f.mesh.blocks[3].tags == Int32[7]
        # property ids surface on every block
        @test [first(b.tags) for b in f.mesh.blocks] ==
              Int32[5, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12]
        @test f.element_tags[1] == UInt64[1, 2, 3, 4]
        @test f.element_tags[end] == UInt64[16]
    end

    @testset "BDF reader: free format, exponents, cp frames" begin
        dir = mktempdir()
        path = joinpath(dir, "free.bdf")
        Base.write(path, join([
            "\$ free-format deck",
            "BEGIN BULK",
            "GRID,1,,0.,0.,0.",
            "GRID,2,,1.d0,0.,0.",
            "GRID,3,,0.,1.0+0,0.",
            "GRID,4,,0.,0.,2.5-1",
            "GRID 5 2.0 0.5 0.25",                # whitespace free form
            "GRID,9,5,1.,1.,1.",                  # alternate frame: skipped
            "CTRIA3,1,7,1,2,3",
            "CTRIA3,2,7,1,2,5",
            "CTETRA,3,9,1,2,3,4",
            "ENDDATA",
            "GRID,900,,9.,9.,9.",                 # ignored after ENDDATA
        ], "\n") * "\n")
        f = @test_logs (:warn,) match_mode=:any read_bdf(path)
        @test f.node_tags == UInt64[1, 2, 3, 4, 5]
        @test f.mesh.coords[:, 2] == [1.0, 0.0, 0.0]
        @test f.mesh.coords[:, 3] == [0.0, 1.0, 0.0]
        @test f.mesh.coords[:, 4] == [0.0, 0.0, 0.25]
        @test f.mesh.coords[:, 5] == [2.0, 0.5, 0.25]
        @test [b.msh for b in f.mesh.blocks] == [2, 4]
        @test f.mesh.blocks[1].nodes == Int32[1 1; 2 2; 3 5]
        @test [b.tags for b in f.mesh.blocks] == [Int32[7, 7], Int32[9]]
        # an element that references a skipped cp-frame node is a hard error
        bad = joinpath(dir, "dangling_cp.bdf")
        Base.write(bad, "GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\n" *
                        "GRID,9,5,0.,1.,0.\nCTRIA3,1,7,1,2,9\nENDDATA\n")
        @test_throws ArgumentError read_bdf(bad)
    end

    @testset "BDF writer round-trips" begin
        dir = mktempdir()
        # linear mesh: line + tri + quad + tet + prism + hex + pyr
        coords = Float64[0 1 1 0 0 1 1 0 0.5 0.5 0 1 1 0 0 1 1 0 0.5 0.5;
                         0 0 1 1 0 0 1 1 0.5 0.5 0 0 1 1 0 0 1 1 0.5 0.5;
                         0 0 0 0 1 1 1 1 0.5 0.5 0 0 0 0 1 1 1 1 2.0 2.0]
        blocks = ElementBlock[
            ElementBlock(1, reshape(Int32[1, 2], 2, 1), Int32[4]),
            ElementBlock(2, Int32[1 1; 2 3; 4 4], Int32[7, 7]),
            ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[8]),
            ElementBlock(4, reshape(Int32[1, 2, 3, 5], 4, 1), Int32[9]),
            ElementBlock(6, reshape(Int32[1, 2, 3, 5, 6, 7], 6, 1), Int32[10]),
            ElementBlock(5, reshape(Int32[1, 2, 3, 4, 5, 6, 7, 8], 8, 1),
                         Int32[11]),
            ElementBlock(7, reshape(Int32[1, 2, 3, 4, 9], 5, 1), Int32[12])]
        mesh = MixedMesh(coords, blocks)
        path = joinpath(dir, "lin.bdf")
        @test write_bdf(path, mesh) == path
        text = read(path, String)
        @test occursin("BEGIN BULK", text) && occursin("ENDDATA", text)
        @test occursin("CTETRA  ", text) && occursin("CPYRA   ", text)
        f = read_bdf(path)
        @test f.mesh.coords == coords
        @test [b.msh for b in f.mesh.blocks] == [1, 2, 3, 4, 6, 5, 7]
        @test [b.nodes for b in f.mesh.blocks] ==
              [b.nodes for b in blocks]
        @test [b.tags for b in f.mesh.blocks] ==
              [b.tags for b in blocks]
        # second-order mesh: tri6 + quad8 + tet10 + penta15 + hex20 + pyra13
        c2 = zeros(3, 20)
        c2[1, :] = 1:20
        blocks2 = ElementBlock[
            ElementBlock(9, reshape(Int32[1:6;], 6, 1), Int32[3]),
            ElementBlock(16, reshape(Int32[1:8;], 8, 1), Int32[4]),
            ElementBlock(11, reshape(Int32[1:10;], 10, 1), Int32[5]),
            ElementBlock(18, reshape(Int32[1:15;], 15, 1), Int32[6]),
            ElementBlock(17, reshape(Int32[1:20;], 20, 1), Int32[7]),
            ElementBlock(19, reshape(Int32[1:13;], 13, 1), Int32[8])]
        mesh2 = MixedMesh(c2, blocks2)
        p2 = joinpath(dir, "ord2.bdf")
        write_bdf(p2, mesh2)
        f2 = read_bdf(p2)
        @test [b.msh for b in f2.mesh.blocks] == [9, 16, 11, 18, 17, 19]
        @test [b.nodes for b in f2.mesh.blocks] ==
              [b.nodes for b in blocks2]
        @test [b.tags for b in f2.mesh.blocks] ==
              [b.tags for b in blocks2]
        # `.nas` shares the writer/reader path
        nas = joinpath(dir, "lin.nas")
        write_bdf(nas, mesh)
        @test [b.msh for b in read_bdf(nas).mesh.blocks] == [1, 2, 3, 4, 6, 5, 7]
        # empty mesh writes an empty deck
        ep = joinpath(dir, "empty.bdf")
        write_bdf(ep, MixedMesh(zeros(3, 0), ElementBlock[]))
        @test size(read_bdf(ep).mesh.coords, 2) == 0
    end

    @testset "BDF reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".bdf"))
            Base.write(p, text)
            return p
        end
        header = "GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nGRID,3,,0.,1.,0.\n"
        @test_throws ArgumentError read_bdf(_w(
            "INCLUDE 'other.bdf'\n" * header))
        @test_throws ArgumentError read_bdf(_w(           # dangling node
            header * "CTRIA3,1,7,1,2,9\nENDDATA\n"))
        @test_throws ArgumentError read_bdf(_w(           # duplicate GRID
            "GRID,1,,0.,0.,0.\nGRID,1,,1.,0.,0.\nENDDATA\n"))
        @test_throws ArgumentError read_bdf(_w(           # duplicate EID
            header * "CBAR,1,5,1,2\nCTRIA3,1,7,1,2,3\nENDDATA\n"))
        @test_throws ArgumentError read_bdf(_w(           # bad integer
            "GRID,xx,,0.,0.,0.\nENDDATA\n"))
        @test_throws ArgumentError read_bdf(_w(           # bad real
            "GRID,1,,zz,0.,0.\nENDDATA\n"))
        @test_throws ArgumentError read_bdf(_w(           # negative pid
            header * "CTRIA3,1,-2,1,2,3\nENDDATA\n"))
        @test_throws ArgumentError read_bdf(_w(           # fixed-field cp junk
            _bdf8("GRID", "1", "2.5", "0.", "0.", "0.")))
        @test_throws ArgumentError read_bdf(_w(
            header * "CTRIA3,1,7,1,2,3\nENDDATA\n"); max_elements=0)
        @test_throws ArgumentError read_bdf(joinpath(dir, "missing.bdf"))
        # an unrecognized C-card warns once but is skipped
        warn_path = joinpath(dir, "warn.bdf")
        Base.write(warn_path, header * "CFOOCARD,1,2,3\n" *
                   "CTRIA3,1,7,1,2,3\nENDDATA\n")
        w = @test_logs (:warn,) match_mode=:any read_bdf(warn_path)
        @test only(w.mesh.blocks).msh == 2
    end

    @testset "writers reject unrepresentable blocks" begin
        dir = mktempdir()
        quad9 = MixedMesh(zeros(3, 9),
            [ElementBlock(10, reshape(Int32[1:9;], 9, 1), Int32[1])])
        @test_throws ArgumentError write_su2(joinpath(dir, "q.su2"), quad9)
        @test_throws ArgumentError write_bdf(joinpath(dir, "q.bdf"), quad9)
        special = MixedMesh(zeros(3, 2),
            [SpecialElementBlock(67, reshape(Int32[1, 2], 2, 1))])
        @test_throws ArgumentError write_su2(joinpath(dir, "s.su2"), special)
        @test_throws ArgumentError write_bdf(joinpath(dir, "s.bdf"), special)
        # SU2 has no 1-D or 0-D domain
        line_only = MixedMesh(zeros(3, 2),
            [ElementBlock(1, reshape(Int32[1, 2], 2, 1), Int32[1])])
        @test_throws ArgumentError write_su2(joinpath(dir, "l.su2"), line_only)
        points_only = MixedMesh(zeros(3, 2),
            [ElementBlock(15, reshape(Int32[1, 2], 1, 2), Int32[1, 1])])
        @test_throws ArgumentError write_su2(joinpath(dir, "p.su2"), points_only)
        # BDF tolerates (skips) point cells: the same mesh writes an empty deck
        pp = joinpath(dir, "p.bdf")
        write_bdf(pp, points_only)
        @test isempty(read_bdf(pp).mesh.blocks)
        # a 3-node line is not representable in either format
        seg3 = MixedMesh(zeros(3, 3),
            [ElementBlock(8, reshape(Int32[1, 2, 3], 3, 1), Int32[1])])
        @test_throws ArgumentError write_su2(joinpath(dir, "t.su2"), seg3)
        @test_throws ArgumentError write_bdf(joinpath(dir, "t.bdf"), seg3)
        tris = MixedMesh(zeros(3, 3),
            [ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[1])])
        @test_throws ArgumentError write_su2(joinpath(dir, "nodir", "x.su2"),
            tris)
        @test_throws ArgumentError write_bdf(joinpath(dir, "nodir", "x.bdf"),
            points_only)
    end

    @testset "session round-trip through API.write/API.open" begin
        _SU2_BDF_API.initialize()
        try
            _SU2_BDF_API.option("General.Terminal", 0)
            _SU2_BDF_API.option("Mesh.MeshSizeFactor", 8.0)
            _SU2_BDF_API.option("Mesh.SaveAll", 1.0)
            _SU2_BDF_API.model.add("su2box")
            _SU2_BDF_API.model.add_box(0, 0, 0, 1, 1, 1; tag=1)
            _SU2_BDF_API.model.add_physical_group(3, [1]; tag=5, name="vol")
            _SU2_BDF_API.model.add_physical_group(2, [1, 2]; tag=6, name="side")
            generated = _SU2_BDF_API.mesh.generate(3)
            @test size(generated.coords, 2) > 0
            source_nodes = _SU2_BDF_API.mesh.get_nodes()
            source_cells = _cells_by_type(_SU2_BDF_API.mesh.get_elements())
            source_coord_cells = _cells_coord_sets(
                _SU2_BDF_API.mesh.get_elements(), source_nodes[1],
                source_nodes[2])
            @test 4 in keys(source_cells)
            @test 2 in keys(source_cells)
            @test 1 in keys(source_cells)

            mktempdir() do dir
                # Write every format while the generated model is current —
                # `API.write` always exports the current model.
                for ext in ("su2", "bdf", "nas")
                    @test _SU2_BDF_API.write(joinpath(dir, "box." * ext)) ===
                        nothing
                    @test isfile(joinpath(dir, "box." * ext))
                end
                # BDF keeps node/element UIDs and physical tags; lines and
                # boundary tris survive; point cells are dropped.
                for ext in ("bdf", "nas")
                    _SU2_BDF_API.model.add("re" * ext)
                    @test _SU2_BDF_API.open(joinpath(dir, "box." * ext)) ===
                        nothing
                    @test "box" in _SU2_BDF_API.model.list()
                    re_nodes = _SU2_BDF_API.mesh.get_nodes()
                    @test re_nodes[1] == source_nodes[1]          # node UIDs
                    # 8-column fields round-trip to ~8 significant digits
                    @test isapprox(re_nodes[2], source_nodes[2]; atol=1e-7)
                    re_cells = _cells_by_type(_SU2_BDF_API.mesh.get_elements())
                    for msh in (Int32(1), Int32(2), Int32(4))
                        @test re_cells[msh] == source_cells[msh]
                    end
                    @test !haskey(re_cells, Int32(15))
                    groups = _SU2_BDF_API.model.get_physical_groups()
                    @test (2, 6) in groups && (3, 5) in groups
                    # property ids carry no names
                    @test _SU2_BDF_API.model.get_physical_name(3, 5) == ""
                end
                # SU2 keeps tets and boundary tris only; physical groups come
                # back as markers (sorted by tag: the unmarked boundary first,
                # then `side`), with marker names preserved.
                _SU2_BDF_API.model.add("resu2")
                @test _SU2_BDF_API.open(joinpath(dir, "box.su2")) === nothing
                re_nodes = _SU2_BDF_API.mesh.get_nodes()
                @test re_nodes[2] == source_nodes[2]              # coordinates
                re_coord_cells = _cells_coord_sets(
                    _SU2_BDF_API.mesh.get_elements(), re_nodes[1], re_nodes[2])
                @test re_coord_cells[Int32(2)] == source_coord_cells[Int32(2)]
                @test re_coord_cells[Int32(4)] == source_coord_cells[Int32(4)]
                @test !haskey(re_coord_cells, Int32(1))
                @test !haskey(re_coord_cells, Int32(15))
                groups = _SU2_BDF_API.model.get_physical_groups()
                @test (2, 1) in groups && (2, 2) in groups
                @test !((3, 5) in groups)
                @test _SU2_BDF_API.model.get_physical_name(2, 1) == "unmarked"
                @test _SU2_BDF_API.model.get_physical_name(2, 2) == "side"
            end
        finally
            _SU2_BDF_API.finalize()
        end
    end

    @testset "session merge keeps the host model intact" begin
        _SU2_BDF_API.initialize()
        try
            _SU2_BDF_API.option("General.Terminal", 0)
            _SU2_BDF_API.model.add("host")
            _SU2_BDF_API.model.add_discrete_entity(1, 7)
            _SU2_BDF_API.model.add_discrete_entity(2, 9)
            _SU2_BDF_API.mesh.add_nodes(1, 7, [1, 2], [0.0, 0, 0, 1, 0, 0])
            _SU2_BDF_API.mesh.add_nodes(2, 9, [3, 4], [0.5, 0, 0, 0.5, 0.5, 0])
            _SU2_BDF_API.mesh.add_elements_by_type(7, 1, [11], [1, 2])
            _SU2_BDF_API.mesh.add_elements_by_type(9, 2, [12], [2, 3, 4])
            host_entities = _SU2_BDF_API.model.get_entities()
            @test (1, 7) in host_entities && (2, 9) in host_entities

            mktempdir() do dir
                su2 = joinpath(dir, "add.su2")
                Base.write(su2, _SU2_FIXTURE)
                @test _SU2_BDF_API.merge(su2) === nothing
                after_su2 = _SU2_BDF_API.model.get_entities()
                @test issubset(host_entities, after_su2)
                @test (3, 1) in after_su2             # NELEM domain entity
                @test (2, 1) in after_su2 && (2, 2) in after_su2
                groups = _SU2_BDF_API.model.get_physical_groups()
                @test (2, 1) in groups && (2, 2) in groups
                @test _SU2_BDF_API.model.get_physical_name(2, 1) == "wall"
                @test _SU2_BDF_API.model.get_physical_name(2, 2) == "inlet"

                bdf = joinpath(dir, "add.bdf")
                Base.write(bdf, join([
                    "BEGIN BULK",
                    "GRID,1,,0.,0.,0.",
                    "GRID,2,,1.,0.,0.",
                    "GRID,3,,0.,1.,0.",
                    "GRID,4,,0.,0.,1.",
                    "CTRIA3,1,42,1,2,3",
                    "CTETRA,2,43,1,2,3,4",
                    "ENDDATA",
                ], "\n") * "\n")
                @test _SU2_BDF_API.merge(bdf) === nothing
                after_bdf = _SU2_BDF_API.model.get_entities()
                @test issubset(after_su2, after_bdf)
                @test (2, 42) in after_bdf            # BDF property → entity
                @test (3, 43) in after_bdf
                groups = _SU2_BDF_API.model.get_physical_groups()
                @test (2, 42) in groups && (3, 43) in groups
            end
        finally
            _SU2_BDF_API.finalize()
        end
    end

    @testset "session-level dispatch guards" begin
        _SU2_BDF_API.initialize()
        try
            _SU2_BDF_API.option("General.Terminal", 0)
            _SU2_BDF_API.model.add("guards")
            mktempdir() do dir
                for ext in ("su2", "bdf", "nas")
                    missing = joinpath(dir, "missing." * ext)
                    @test_throws ArgumentError _SU2_BDF_API.open(missing)
                    @test_throws ArgumentError _SU2_BDF_API.merge(missing)
                end
                bad = joinpath(dir, "x.foo")
                @test_throws ArgumentError _SU2_BDF_API.write(bad)
            end
        finally
            _SU2_BDF_API.finalize()
        end
    end

    # ── SU2 keyword variants ──────────────────────────────────────────────────
    # Upstream Gmsh 4.15.2 has no SU2 reader at all (only writeSU2): the
    # reader scope here is a Tessella extension. Multi-zone and periodic
    # decks cannot be expressed through the single-domain bridge and must
    # reject loudly; ancillary SU2 solver keys (MARKER_PERIODIC etc.) are
    # skipped like upstream's unknown-keyword scan.
    @testset "su2 zone/marker variants" begin
        _SU2_BDF_API.initialize()
        try
            _SU2_BDF_API.option("General.Terminal", 0)
            mktempdir() do dir
                # Multi-zone and periodic forms reject with the key named.
                for (key, deck) in (
                    ("NZONE", "NDIME= 2\nNZONE= 2\nIZONE= 1\nNPOIN= 0\n"),
                    ("IZONE", "NDIME= 2\nIZONE= 1\nNPOIN= 0\n"),
                    ("NPERIODIC", "NDIME= 2\nNPOIN= 0\nNPERIODIC= 1\n"))
                    p = joinpath(dir, key * ".su2")
                    write(p, deck)
                    _SU2_BDF_API.model.add("z" * key)
                    err = try
                        _SU2_BDF_API.open(p); nothing
                    catch e; e; end
                    @test err isa ArgumentError
                    @test occursin(key, sprint(showerror, err))
                end
                # Marker elements carrying the optional trailing element index
                # plus an ignored solver key between marker sections.
                p = joinpath(dir, "marks.su2")
                write(p, """
NDIME= 2
NPOIN= 3
0.0 0.0 0
1.0 0.0 1
0.0 1.0 2
NELEM= 1
5 0 1 2 0
NMARK= 2
MARKER_TAG= wall
MARKER_ELEMS= 1
3 0 1 7
MARKER_PERIODIC= NO
MARKER_TAG= inlet
MARKER_ELEMS= 1
3 1 2 8
""")
                _SU2_BDF_API.model.add("marks")
                @test _SU2_BDF_API.open(p) === nothing
                names = Dict(_SU2_BDF_API.model.get_physical_groups()
                    .=> [_SU2_BDF_API.model.get_physical_name(d, t)
                         for (d, t) in _SU2_BDF_API.model.get_physical_groups()])
                @test names[(1, 1)] == "wall" && names[(1, 2)] == "inlet"
            end
        finally
            _SU2_BDF_API.finalize()
        end
    end
end
