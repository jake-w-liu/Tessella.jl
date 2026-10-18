# ── .inp (Abaqus input deck) and .neu (Gambit neutral) mesh I/O ───────────────
#
# Correctness  : hand-written fixtures transcribe the documented layouts
#                (Abaqus `*NODE`/`*ELEMENT` keyword blocks with `*ELSET`/NSET
#                set membership and `**` comments; Gambit `CONTROL INFO` /
#                `NODAL COORDINATES` / `ELEMENTS/CELLS` / `ELEMENT GROUP` /
#                `BOUNDARY CONDITIONS` sections ended by `ENDOFSECTION`).
#                Session round-trips export with `API.write` and re-import
#                with `API.open`, checking physical tags, names and entity
#                recovery; `API.merge` offsets file-local labels.
# Robustness   : `D` exponents, trailing-comma continuations and suffix
#                modifiers on read; unknown element types warn+skip; dangling
#                node/set/boundary references, duplicate labels, bad arities,
#                `*INCLUDE`, missing counts and truncated sections reject;
#                unrepresentable element types reject loudly on write.
# Completeness : `.inp` `Volume|Surface|Line|Point<n>` element sets →
#                elementary entities, `Physical<Dim><tag>` names pin physical
#                tags, other set names become named physicals; `.neu` group
#                ordinals → elementary entities, `GROUP:` numbers → physical
#                tags, boundary sets → face/edge cells with their own
#                physical; names round-trip through `mesh.physical_names`.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock, SpecialElementBlock,
    msh_num_nodes
using Tessella.IO

const _INP_NEU_API=Tessella.API

# Hand-written Abaqus deck: comments/`*Heading` skip, `D` exponents, a
# continuation record (element 2), an entity-convention ELSET (`Volume3`,
# `Surface2`), a `Physical<Dim><tag>` set, a custom-named set, `GENERATE`
# ranges, an unknown element type (warn+skip) and a skipped `*Step` block.
const _INP_FIXTURE = """
** Tessella test deck
*Heading
 hand-written fixture
*Node, Nset=ALLNODES
      1,    0.0,    0.0,    0.0
      2,  1.0D0,    0.0,    0.0
      3,    0.0,    1.0,    0.0
      4,    0.0,    0.0,    1.0
      5,    0.5,    0.5,    0.5
*Element, type=C3D4, ELSET=Volume3
 1, 1, 2, 3, 4
*Element, type=C3D5,
 2, 1, 2, 3, 5,
 4
*Element, type=CPS3, ELSET=Surface2
 3, 1, 2, 3
*Element, type=B21
 4, 1, 2
*Element, type=FAKETYPE42
 5, 1, 2
*Elset, elset=PhysicalVolume7
 1, 2
*Elset, elset=mywall
 3
*Elset, elset=PhysicalLine4, GENERATE
 4, 4, 1
*Nset, Nset=corners, generate
 1, 3, 1
*Step
*Static
1.0, 1.0
*End Step
"""

# Hand-written Gambit deck: hex (type 4, wrapped 8-node record), pyramid
# (type 7), two groups (one with NFLAGS 0), one boundary quad face on the
# hex's bottom side.
const _NEU_FIXTURE = """
        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Tessella test deck
PROGRAM:                Gambit     VERSION:  2.4.6
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         9         2         2         1         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   1.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         5   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         6   1.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         7   0.00000000000e+00   1.00000000000e+00   1.00000000000e+00
         8   1.00000000000e+00   1.00000000000e+00   1.00000000000e+00
         9   5.00000000000e-01   5.00000000000e-01   2.00000000000e+00
ENDOFSECTION
      ELEMENTS/CELLS 2.0.0
       1  4  8        1       2       4       3
               5       6       8       7
       2  7  5        5       6       8       7       9
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:          5 ELEMENTS:          1 MATERIAL:          0 NFLAGS:          1
                            fluid
       0
       1
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:          6 ELEMENTS:          1 MATERIAL:          0 NFLAGS:          0
                            roof
       2
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
                            wall         1         1         0        51
         1    4    5
ENDOFSECTION
"""

# Small 2-D Gambit deck: two triangles, one group, one edge boundary record.
const _NEU_2D_FIXTURE = """
        CONTROL INFO 2.0.0
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         2         1         1         2         2
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00
         4   1.00000000000e+00   1.00000000000e+00
ENDOFSECTION
      ELEMENTS/CELLS 2.0.0
       1  3  3        1       2       3
       2  3  3        2       4       3
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:          9 ELEMENTS:          2 MATERIAL:          0 NFLAGS:          0
                          sheet
       1       2
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
                        edgebc         1         1         0        44
         1    3    1
ENDOFSECTION
"""

@testset "INP/NEU mesh I/O" begin

    @testset "INP reader: keywords, sets, names" begin
        dir = mktempdir()
        path = joinpath(dir, "fixture.inp")
        Base.write(path, _INP_FIXTURE)
        f = @test_logs (:warn,) match_mode=:any read_inp(path)
        @test f isa InpFile
        @test f.node_tags == UInt64[1, 2, 3, 4, 5]
        @test f.mesh.coords[:, 2] == [1.0, 0.0, 0.0]    # 1.0D0
        @test f.mesh.coords[:, 5] == [0.5, 0.5, 0.5]
        # FAKETYPE42 was skipped; blocks keep file order by type.
        @test [b.msh for b in f.mesh.blocks] == [4, 7, 2, 1]
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 2, 3, 5, 4], 5, 1)
        @test f.mesh.blocks[3].nodes == reshape(Int32[1, 2, 3], 3, 1)
        @test f.mesh.blocks[4].nodes == reshape(Int32[1, 2], 2, 1)
        # Volume3 / Surface2 ELSET conventions → elementary entity tags.
        @test f.mesh.elementary_entities ==
              [Int32[3], Int32[0], Int32[2], Int32[0]]
        # PhysicalVolume7 pins tag 7; mywall gets the lowest free dim-2 tag;
        # PhysicalLine4 pins tag 4 (GENERATE range resolved).
        @test [b.tags for b in f.mesh.blocks] ==
              [Int32[7], Int32[7], Int32[1], Int32[4]]
        @test f.mesh.physical_names == Dict((2, 1) => "mywall")
        @test f.element_tags ==
              [UInt64[1], UInt64[2], UInt64[3], UInt64[4]]
    end

    @testset "INP reader: suffix modifiers and set precedence" begin
        dir = mktempdir()
        path = joinpath(dir, "suffix.inp")
        Base.write(path, """
            *Node
            1, 0., 0., 0.
            2, 1., 0., 0.
            3, 0., 1., 0.
            4, 0., 0., 1.
            *Element, type=C3D4R
            1, 1, 2, 3, 4
            *Element, type=S4R
            2, 1, 2, 3, 4
            *Elset, elset=first
            1
            *Elset, elset=second
            1
            """)
        f = read_inp(path)
        @test [b.msh for b in f.mesh.blocks] == [4, 3]   # R-suffixes stripped
        # first physical set wins on overlap
        @test f.mesh.blocks[1].tags == Int32[1]
        @test f.mesh.physical_names[(3, 1)] == "first"
    end

    @testset "INP writer round-trips" begin
        dir = mktempdir()
        coords = Float64[0 1 0 0 0.5 0 1 1 0 0 1 1 0 0.5;
                         0 0 1 0 0.5 0 0 1 1 0 0 1 1 0.5;
                         0 0 0 1 0.5 0 0 0 0 1 1 1 1 2.0]
        blocks = ElementBlock[
            ElementBlock(1, reshape(Int32[1, 2], 2, 1), Int32[4]),
            ElementBlock(2, Int32[1 1; 2 2; 3 5], Int32[6, 6]),
            ElementBlock(4, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[5]),
            ElementBlock(7, reshape(Int32[10, 11, 12, 13, 14], 5, 1),
                         Int32[0])]
        mesh = MixedMesh(coords, blocks;
            physical_names=Dict((3, 5) => "vol", (1, 4) => "rod"),
            elementary_entities=[Int32[1], Int32[9, 9], Int32[3], Int32[0]])
        path = joinpath(dir, "deck.inp")
        @test write_inp(path, mesh) == path
        text = read(path, String)
        @test occursin("*Heading", text) && occursin("*NODE", text)
        @test occursin("*ELEMENT, type=C3D4, ELSET=Volume3", text)
        @test occursin("*ELEMENT, type=CPS3, ELSET=Surface9", text)
        @test occursin("*ELEMENT, type=T3D2, ELSET=Line1", text)
        @test occursin("*ELEMENT, type=C3D5, ELSET=Volume0", text)
        @test occursin("*ELSET,ELSET=vol", text)
        @test occursin("*ELSET,ELSET=PhysicalSurface6", text)
        @test occursin("*NSET,NSET=rod", text)
        f = read_inp(path)
        @test f.mesh.coords == coords
        # *ELEMENT blocks emit sorted by (dim, entity, type)
        @test [b.msh for b in f.mesh.blocks] == [1, 2, 4, 7]
        @test f.mesh.elementary_entities ==
              [Int32[1], Int32[9, 9], Int32[3], Int32[0]]
        # named sets keep their name but take the lowest free tag; the
        # unnamed dim-2 physical round-trips its tag 6 via the
        # `PhysicalSurface6` fallback name.
        @test [b.tags for b in f.mesh.blocks] ==
              [Int32[1], Int32[6, 6], Int32[1], Int32[0]]
        @test f.mesh.physical_names ==
              Dict((3, 1) => "vol", (1, 1) => "rod")
        # second-order connectivity permutations invert exactly
        c2 = zeros(3, 27); c2[1, :] = 1:27
        blocks2 = ElementBlock[
            ElementBlock(8, reshape(Int32[1, 2, 3], 3, 1), Int32[1]),
            ElementBlock(11, reshape(Int32[1:10;], 10, 1), Int32[2]),
            ElementBlock(17, reshape(Int32[1:20;], 20, 1), Int32[3]),
            ElementBlock(18, reshape(Int32[1:15;], 15, 1), Int32[4]),
            ElementBlock(19, reshape(Int32[1:13;], 13, 1), Int32[5]),
            ElementBlock(12, reshape(Int32[1:27;], 27, 1), Int32[6]),
            ElementBlock(10, reshape(Int32[1:9;], 9, 1), Int32[7])]
        mesh2 = MixedMesh(c2, blocks2)
        p2 = joinpath(dir, "ord2.inp")
        write_inp(p2, mesh2)
        text2 = read(p2, String)
        @test occursin("type=C3D10", text2) && occursin("type=C3D20", text2)
        @test occursin("type=C3D15", text2) && occursin("type=C3D13", text2)
        @test occursin("type=C3D27", text2) && occursin("type=M3D9", text2)
        @test occursin("type=T3D3", text2)
        f2 = read_inp(p2)
        bymsh = Dict(b.msh => b for b in f2.mesh.blocks)
        @test sort(collect(keys(bymsh))) == [8, 10, 11, 12, 17, 18, 19]
        for b in blocks2
            @test bymsh[b.msh].nodes == b.nodes
            @test bymsh[b.msh].tags == b.tags
        end
    end

    @testset "INP reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".inp"))
            Base.write(p, text)
            return p
        end
        nodes = "*Node\n1, 0., 0., 0.\n2, 1., 0., 0.\n3, 0., 1., 0.\n"
        @test_throws ArgumentError read_inp(_w(           # missing file body
            "*INCLUDE, INPUT=other.inp\n" * nodes))
        @test_throws ArgumentError read_inp(_w(           # INPUT= deck ref
            "*Node, INPUT=other.inp\n"))
        @test_throws ArgumentError read_inp(_w(           # *ELEMENT w/o TYPE
            nodes * "*Element\n1, 1, 2, 3\n"))
        @test_throws ArgumentError read_inp(_w(           # dangling node
            nodes * "*Element, type=CPS3\n1, 1, 2, 9\n"))
        @test_throws ArgumentError read_inp(_w(           # duplicate node
            "*Node\n1, 0., 0.\n1, 0., 0.\n"))
        @test_throws ArgumentError read_inp(_w(           # duplicate element
            nodes * "*Element, type=CPS3\n1, 1, 2, 3\n" *
                    "*Element, type=CPS3\n1, 2, 3, 1\n"))
        @test_throws ArgumentError read_inp(_w(           # bad arity
            nodes * "*Element, type=CPS3\n1, 1, 2, 3, 1\n"))
        @test_throws ArgumentError read_inp(_w(           # bad integer
            nodes * "*Element, type=CPS3\nx, 1, 2, 3\n"))
        @test_throws ArgumentError read_inp(_w(           # bad coordinate
            "*Node\n1, q, 0., 0.\n"))
        @test_throws ArgumentError read_inp(_w(           # unnamed *ELSET
            nodes * "*Element, type=CPS3\n1, 1, 2, 3\n*Elset\n1\n"))
        @test_throws ArgumentError read_inp(_w(           # bad GENERATE range
            nodes * "*Elset, elset=s, generate\n3, 1\n"))
        @test_throws ArgumentError read_inp(_w(           # exploding GENERATE
            nodes * "*Elset, elset=s, generate\n1, 999999999\n"))
        @test_throws ArgumentError read_inp(joinpath(dir, "missing.inp"))
        @test_throws ArgumentError read_inp(_w(nodes); max_nodes=2)
        @test_throws ArgumentError read_inp(
            _w(nodes * "*Element, type=CPS3\n1, 1, 2, 3\n"); max_elements=0)
        # unknown element type warns once and is skipped
        warn_path = joinpath(dir, "warn.inp")
        Base.write(warn_path, nodes *
            "*Element, type=UNKNOWN9\n1, 1, 2\n*Element, type=CPS3\n1, 1, 2, 3\n")
        w = @test_logs (:warn,) match_mode=:any read_inp(warn_path)
        @test only(w.mesh.blocks).msh == 2
        # a dangling set membership warns but survives
        dang = joinpath(dir, "dang.inp")
        Base.write(dang, nodes * "*Element, type=CPS3\n1, 1, 2, 3\n" *
            "*Elset, elset=s\n1, 99\n")
        d = @test_logs (:warn,) match_mode=:any read_inp(dang)
        @test d.mesh.blocks[1].tags == Int32[1]
        @test d.mesh.physical_names[(2, 1)] == "s"
    end

    @testset "NEU reader: sections, groups, boundaries" begin
        dir = mktempdir()
        path = joinpath(dir, "fixture.neu")
        Base.write(path, _NEU_FIXTURE)
        f = read_neu(path)
        @test f isa NeuFile
        @test f.node_tags == UInt64.(1:9)
        @test f.mesh.coords[:, 9] == [0.5, 0.5, 2.0]
        # hex + pyramid cells plus the generated boundary quad
        @test [b.msh for b in f.mesh.blocks] == [5, 7, 3]
        # Gambit 3↔4 base-corner swap recovered to Gmsh order
        @test f.mesh.blocks[1].nodes == reshape(Int32[1:8;], 8, 1)
        @test f.mesh.blocks[2].nodes == reshape(Int32[5, 6, 7, 8, 9], 5, 1)
        # boundary record `1 4 5`: hex element 1, side 5 → bottom quad
        @test f.mesh.blocks[3].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        @test f.mesh.elementary_entities ==
              [Int32[1], Int32[2], Int32[1]]
        @test [b.tags for b in f.mesh.blocks] ==
              [Int32[5], Int32[6], Int32[1]]
        @test f.mesh.physical_names[(3, 5)] == "fluid"
        @test f.mesh.physical_names[(3, 6)] == "roof"
        @test f.mesh.physical_names[(2, 1)] == "wall"
    end

    @testset "NEU reader: 2-D deck" begin
        dir = mktempdir()
        path = joinpath(dir, "flat.neu")
        Base.write(path, _NEU_2D_FIXTURE)
        f = read_neu(path)
        @test size(f.mesh.coords, 2) == 4
        @test all(iszero, f.mesh.coords[3, :])
        @test [b.msh for b in f.mesh.blocks] == [2, 1]
        @test f.mesh.blocks[1].nodes == Int32[1 2; 2 4; 3 3]
        # boundary `1 3 1`: tri element 1, side 1 → edge {1,2}
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 2], 2, 1)
        @test f.mesh.blocks[1].tags == Int32[9, 9]
        @test f.mesh.blocks[2].tags == Int32[1]
        @test f.mesh.physical_names[(2, 9)] == "sheet"
        @test f.mesh.physical_names[(1, 1)] == "edgebc"
    end

    @testset "NEU writer round-trips" begin
        dir = mktempdir()
        # tet with a boundary tri on its {1,2,4} face (side table resolves it)
        coords = Float64[0 1 0 0;
                         0 0 1 0;
                         0 0 0 1]
        blocks = ElementBlock[
            ElementBlock(4, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[5]),
            ElementBlock(2, reshape(Int32[1, 2, 4], 3, 1), Int32[9])]
        mesh = MixedMesh(coords, blocks;
            physical_names=Dict((3, 5) => "vol", (2, 9) => "inlet"))
        path = joinpath(dir, "m.neu")
        @test write_neu(path, mesh) == path
        text = read(path, String)
        @test occursin("CONTROL INFO", text)
        @test occursin("NUMNP     NELEM", text)
        @test occursin("NODAL COORDINATES", text)
        @test occursin("ELEMENTS/CELLS", text)
        @test occursin("ELEMENT GROUP", text)
        @test occursin("BOUNDARY CONDITIONS", text)
        @test occursin("GROUP:", text)
        @test occursin("inlet", text)
        f = read_neu(path)
        @test f.mesh.coords == coords
        @test [b.msh for b in f.mesh.blocks] == [4, 2]
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        # the boundary tri regenerated with physical 9 and its name
        @test sort(vec(f.mesh.blocks[2].nodes)) == Int32[1, 2, 4]
        @test f.mesh.blocks[1].tags == Int32[5]
        @test f.mesh.blocks[2].tags == Int32[1]
        @test f.mesh.physical_names[(3, 5)] == "Material group 5" ||
              f.mesh.physical_names[(3, 5)] == "vol"
        @test f.mesh.physical_names[(2, 1)] == "inlet"
        # a physical name matching the fallback recovers the pinned tag
        mesh2 = MixedMesh(coords, blocks)
        p2 = joinpath(dir, "m2.neu")
        write_neu(p2, mesh2)
        f2 = read_neu(p2)
        @test f2.mesh.blocks[2].tags == Int32[9]   # PhysicalSurface9 fallback
        # 2-D deck: quad + edge boundary cell
        flat = MixedMesh(
            Float64[0 1 0 1; 0 0 1 1; 0 0 0 0],
            [ElementBlock(3, reshape(Int32[1, 2, 4, 3], 4, 1), Int32[2]),
             ElementBlock(1, reshape(Int32[1, 2], 2, 1), Int32[7])])
        p3 = joinpath(dir, "flat.neu")
        write_neu(p3, flat)
        f3 = read_neu(p3)
        @test [b.msh for b in f3.mesh.blocks] == [3, 1]
        @test f3.mesh.blocks[1].tags == Int32[2]
        @test f3.mesh.blocks[2].tags == Int32[7]   # PhysicalLine7 fallback
    end

    @testset "NEU reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".neu"))
            Base.write(p, text)
            return p
        end
        head = "        CONTROL INFO 2.0.0\n" *
               "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n" *
               "         4         1         0         0         3         3\n" *
               "ENDOFSECTION\n"
        coord = "   NODAL COORDINATES 2.0.0\n" *
                "1 0.0 0.0 0.0\n2 1.0 0.0 0.0\n3 0.0 1.0 0.0\n4 0.0 0.0 1.0\n" *
                "ENDOFSECTION\n"
        cells = "      ELEMENTS/CELLS 2.0.0\n1 6 4 1 2 3 4\nENDOFSECTION\n"
        @test_throws ArgumentError read_neu(_w(            # no CONTROL INFO
            coord * cells))
        @test_throws ArgumentError read_neu(_w(            # truncated nodes
            head * "   NODAL COORDINATES 2.0.0\n1 0.0 0.0 0.0\n"))
        @test_throws ArgumentError read_neu(_w(            # missing elements
            head * coord))
        @test_throws ArgumentError read_neu(_w(            # unknown type
            head * coord *
            "      ELEMENTS/CELLS 2.0.0\n1 42 4 1 2 3 4\nENDOFSECTION\n"))
        @test_throws ArgumentError read_neu(_w(            # ndp mismatch
            head * coord *
            "      ELEMENTS/CELLS 2.0.0\n1 6 5 1 2 3 4 2\nENDOFSECTION\n"))
        @test_throws ArgumentError read_neu(_w(            # dangling node
            head * coord *
            "      ELEMENTS/CELLS 2.0.0\n1 6 4 1 2 3 9\nENDOFSECTION\n"))
        @test_throws ArgumentError read_neu(_w(            # dangling group ref
            head * coord * cells *
            "       ELEMENT GROUP 2.0.0\n" *
            "GROUP:          1 ELEMENTS:          1 MATERIAL:          0 " *
            "NFLAGS:          0\ng\n9\nENDOFSECTION\n"))
        @test_throws ArgumentError read_neu(_w(            # bad side number
            head * coord * cells *
            "       BOUNDARY CONDITIONS 2.0.0\nb 1 1 0 0\n1 6 9\n" *
            "ENDOFSECTION\n"))
        @test_throws ArgumentError read_neu(_w(            # wrong parent type
            head * coord * cells *
            "       BOUNDARY CONDITIONS 2.0.0\nb 1 1 0 0\n1 4 1\n" *
            "ENDOFSECTION\n"))
        @test_throws ArgumentError read_neu(joinpath(dir, "missing.neu"))
        @test_throws ArgumentError read_neu(_w(head * coord * cells);
            max_elements=0)
    end

    @testset "writers reject unrepresentable blocks" begin
        dir = mktempdir()
        quad9 = MixedMesh(zeros(3, 9),
            [ElementBlock(10, reshape(Int32[1:9;], 9, 1), Int32[1])])
        @test_throws ArgumentError write_neu(joinpath(dir, "q.neu"), quad9)
        # order-3 tri / full-P2 prism / special records have no Abaqus name
        tri10 = MixedMesh(zeros(3, 10),
            [ElementBlock(21, reshape(Int32[1:10;], 10, 1), Int32[1])])
        @test_throws ArgumentError write_inp(joinpath(dir, "t.inp"), tri10)
        pri18 = MixedMesh(zeros(3, 18),
            [ElementBlock(13, reshape(Int32[1:18;], 18, 1), Int32[1])])
        @test_throws ArgumentError write_inp(joinpath(dir, "p.inp"), pri18)
        special = MixedMesh(zeros(3, 2),
            [SpecialElementBlock(67, reshape(Int32[1, 2], 2, 1))])
        @test_throws ArgumentError write_inp(joinpath(dir, "s.inp"), special)
        @test_throws ArgumentError write_neu(joinpath(dir, "s.neu"), special)
        # .neu carries only linear cells: tri6 has no Gambit code
        tri6 = MixedMesh(zeros(3, 6),
            [ElementBlock(9, reshape(Int32[1:6;], 6, 1), Int32[1])])
        @test_throws ArgumentError write_neu(joinpath(dir, "t6.neu"), tri6)
        # .neu has no point-only meshes
        points = MixedMesh(zeros(3, 2),
            [ElementBlock(15, reshape(Int32[1, 2], 1, 2), Int32[1, 1])])
        @test_throws ArgumentError write_neu(joinpath(dir, "p.neu"), points)
        tris = MixedMesh(zeros(3, 3),
            [ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[1])])
        @test_throws ArgumentError write_inp(
            joinpath(dir, "nodir", "x.inp"), tris)
        @test_throws ArgumentError write_neu(
            joinpath(dir, "nodir", "x.neu"), tris)
    end

    @testset "session round-trip through API.write/API.open" begin
        _INP_NEU_API.initialize()
        try
            _INP_NEU_API.option("General.Terminal", 0)
            _INP_NEU_API.option("Mesh.SaveAll", 1.0)
            _INP_NEU_API.model.add("deck")
            _INP_NEU_API.model.add_discrete_entity(3, 1)
            _INP_NEU_API.model.add_discrete_entity(2, 9)
            _INP_NEU_API.mesh.add_nodes(3, 1, [1, 2, 3, 4],
                [0.0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
            _INP_NEU_API.mesh.add_elements_by_type(1, 4, [10], [1, 2, 3, 4])
            # the boundary tri shares the tet's {1,2,3} face so .neu resolves
            # it as a boundary record
            _INP_NEU_API.mesh.add_elements_by_type(9, 2, [20], [1, 2, 3])
            _INP_NEU_API.model.add_physical_group(3, [1]; tag=5, name="vol")
            _INP_NEU_API.model.add_physical_group(2, [9]; tag=6)
            mktempdir() do dir
                inp = joinpath(dir, "deck.inp")
                neu = joinpath(dir, "deck.neu")
                @test _INP_NEU_API.write(inp) === nothing
                @test _INP_NEU_API.write(neu) === nothing
                text = read(inp, String)
                @test occursin("*ELEMENT, type=C3D4, ELSET=Volume1", text)
                @test occursin("*ELEMENT, type=CPS3, ELSET=Surface9", text)
                @test occursin("*ELSET,ELSET=vol", text)
                @test occursin("*ELSET,ELSET=PhysicalSurface6", text)
                @test occursin("*NSET,NSET=vol", text)

                _INP_NEU_API.model.add("reinp")
                @test _INP_NEU_API.open(inp) === nothing
                groups = _INP_NEU_API.model.get_physical_groups()
                # named "vol" takes the lowest free tag; the unnamed group
                # round-trips tag 6 through the PhysicalSurface6 set name
                @test (3, 1) in groups && (2, 6) in groups
                @test _INP_NEU_API.model.get_physical_name(3, 1) == "vol"
                entities = _INP_NEU_API.model.get_entities()
                @test (3, 1) in entities && (2, 9) in entities
                types, _, _ = _INP_NEU_API.mesh.get_elements()
                @test Int32(4) in types && Int32(2) in types

                _INP_NEU_API.model.add("reneu")
                @test _INP_NEU_API.open(neu) === nothing
                groups = _INP_NEU_API.model.get_physical_groups()
                # GROUP: 5 carries the tag and name; the boundary tri
                # regenerates with the PhysicalSurface6-pinned tag
                @test (3, 5) in groups && (2, 6) in groups
                @test _INP_NEU_API.model.get_physical_name(3, 5) == "vol"
                types, _, _ = _INP_NEU_API.mesh.get_elements()
                @test Int32(4) in types && Int32(2) in types
            end
        finally
            _INP_NEU_API.finalize()
        end
    end

    @testset "session merge keeps the host model intact" begin
        _INP_NEU_API.initialize()
        try
            _INP_NEU_API.option("General.Terminal", 0)
            _INP_NEU_API.model.add("host")
            _INP_NEU_API.model.add_discrete_entity(1, 7)
            _INP_NEU_API.mesh.add_nodes(1, 7, [1, 2], [0.0, 0, 0, 1, 0, 0])
            _INP_NEU_API.mesh.add_elements_by_type(7, 1, [11], [1, 2])
            host_entities = _INP_NEU_API.model.get_entities()
            @test (1, 7) in host_entities

            mktempdir() do dir
                inp = joinpath(dir, "add.inp")
                Base.write(inp, """
                    *Node
                    1, 0., 0., 0.
                    2, 1., 0., 0.
                    3, 0., 1., 0.
                    *Element, type=CPS3, ELSET=Surface4
                    1, 1, 2, 3
                    *Elset, elset=side
                    1
                    """)
                @test _INP_NEU_API.merge(inp) === nothing
                after = _INP_NEU_API.model.get_entities()
                @test issubset(host_entities, after)
                @test (2, 4) in after                  # ELSET=Surface4 entity
                groups = _INP_NEU_API.model.get_physical_groups()
                @test any(g -> g[1] == 2, groups)
                @test "side" in (
                    _INP_NEU_API.model.get_physical_name(g[1], g[2])
                    for g in groups if g[1] == 2)

                neu = joinpath(dir, "add.neu")
                Base.write(neu, _NEU_2D_FIXTURE)
                @test _INP_NEU_API.merge(neu) === nothing
                after2 = _INP_NEU_API.model.get_entities()
                @test issubset(after, after2)
                groups = _INP_NEU_API.model.get_physical_groups()
                @test (2, 9) in groups                 # GROUP: 9 → physical
                @test _INP_NEU_API.model.get_physical_name(2, 9) == "sheet"
            end
        finally
            _INP_NEU_API.finalize()
        end
    end

    @testset "session-level dispatch guards" begin
        _INP_NEU_API.initialize()
        try
            _INP_NEU_API.option("General.Terminal", 0)
            _INP_NEU_API.model.add("guards")
            mktempdir() do dir
                for ext in ("inp", "neu")
                    missing = joinpath(dir, "missing." * ext)
                    @test_throws ArgumentError _INP_NEU_API.open(missing)
                    @test_throws ArgumentError _INP_NEU_API.merge(missing)
                    bad = joinpath(dir, "bad." * ext)
                    Base.write(bad, "garbage\n")
                    @test_throws ArgumentError _INP_NEU_API.open(bad)
                end
            end
        finally
            _INP_NEU_API.finalize()
        end
    end
end
