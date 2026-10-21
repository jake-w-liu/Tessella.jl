# ── .diff (Diffpack GridFE) reader/writer ──
#
# Correctness  : `.diff` is Diffpack's GridFE deck — bidirectional upstream
#                (`GModelIO_DIFF.cpp`). Hand-written decks pin the counted
#                banner headers (`Number of space dim.`/`elements`/`nodes`,
#                `Max number of nodes in an element`, `N Boundary
#                indicators:`), the `num ( x , y , z ) [k] b1 ..` node rows
#                and the `num Elm* subdomain n1 ..` element records with the
#                upstream read `mapping` tables. Writer output pins the GridFE
#                text, the `getVertexDIFF` connectivity placements and the
#                boundary-indicator propagation from dim-1 entities.
# Robustness   : missing/malformed headers and `#` separators, negative or
#                over-limit counts, bad node records (bad int/float, negative
#                node number, short indicator tail), duplicate nodes (first
#                wins, warning), dangling node references, unknown `Elm*`
#                types, negative subdomains, truncated files and the
#                read-time limits.
# Completeness : `API.open`/`API.merge`/`API.write` route `.diff` through the
#                same foreign-mesh staging channel as `.mail`/`.ir3`.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock, MixedEntity, MixedEntityData
using Tessella.IO

const _DIFF_API = Tessella.API

const _F64V = Union{Nothing,Vector{Float64}}

# A full GridFE deck for `nv` node record strings in `node_lines` and element
# record strings `elt_lines`; `bis` lists the banner boundary indicators.
function _diff_deck(nv, node_lines, elt_lines; bis=Int[], maxn=4)
    io = IOBuffer()
    println(io, "\n")
    println(io, " Finite element mesh (GridFE):\n")
    println(io, " Number of space dim. =   3")
    println(io, " Number of elements   =  ", length(elt_lines))
    println(io, " Number of nodes      =  ", nv, "\n")
    println(io, " All elements are of the same type: dpTRUE")
    println(io, " Max number of nodes in an element: ", maxn, " ")
    println(io, " Only one subdomain               : dpFALSE")
    println(io, " Lattice data                     ? 0\n\n\n")
    print(io, " ", length(bis), " Boundary indicators:  ")
    for b in bis
        print(io, " ", b)
    end
    println(io, "\n\n")
    println(io, "#")
    for line in node_lines
        println(io, line)
    end
    println(io, "\n#")
    for line in elt_lines
        println(io, line)
    end
    println(io)
    return String(take!(io))
end

@testset "DIFF mesh reader/writer" begin

    @testset "DIFF reader: hand-written deck" begin
        dir = mktempdir()
        function _w(text, name="t$(abs(hash(text))).diff")
            p = joinpath(dir, name)
            Base.write(p, text)
            return p
        end
        # sparse labels, [k] boundary-indicator tails, material column,
        # element-number column ignored (identity is file order)
        deck = _diff_deck(4,
            [" 10 ( 0 , 0 , 0 ) [1] 7",
             " 20 ( 1 , 0 , 0 ) [2] 7 8",
             " 30 ( 1 , 1 , 0 ) [0]",
             " 40 ( 0 , 1 , 0 ) [1] 8"],
            [" 77 ElmB4n2D 9 10 20 40 30"]; bis=[7, 8])
        f = read_diff(_w(deck))
        @test f isa DiffFile
        @test f.node_tags == UInt64[10, 20, 30, 40]
        @test f.mesh.coords == Float64[0 1 1 0; 0 0 1 1; 0 0 0 0]
        @test [b.msh for b in f.mesh.blocks] == [3]
        # file order 10 20 40 30 → Gmsh slots (1,2,4,3): rows stay dense
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        @test f.mesh.elementary_entities == [Int32[9]]   # material column
        @test f.mesh.blocks[1].tags == Int32[0]          # no physicals
        @test f.element_tags == [UInt64[1]]              # dense file order
        # a mixed-type deck keeps per-record arity and order of appearance
        mixed = _diff_deck(4,
            [" 1 ( 0 , 0 , 0 ) [0]",
             " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 1 , 1 , 0 ) [0]",
             " 4 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 2 1 2 3",
             " 2 ElmB4n2D 3 1 2 3 4"])
        fm = read_diff(_w(mixed))
        @test [b.msh for b in fm.mesh.blocks] == [2, 3]
        @test fm.mesh.elementary_entities == [Int32[2], Int32[3]]
        @test fm.element_tags == [UInt64[1], UInt64[2]]
        # the linear-hex read mapping: file slots feed Gmsh slots
        # (5,6,1,2,8,7,4,3) over sparse node labels
        hexdeck = _diff_deck(8,
            [" $(11i) ( $i , 0 , 0 ) [0]" for i in 1:8],
            [" 5 ElmB8n3D 4 11 22 33 44 55 66 77 88"]; maxn=8)
        fh = read_diff(_w(hexdeck))
        @test [b.msh for b in fh.mesh.blocks] == [5]
        @test fh.node_tags == UInt64[11, 22, 33, 44, 55, 66, 77, 88]
        @test vec(fh.mesh.blocks[1].nodes) ==
              Int32[3, 4, 8, 7, 1, 2, 6, 5]
        @test fh.mesh.elementary_entities == [Int32[4]]
        # tetrahedra read identity
        tet = _diff_deck(4,
            [" $i ( $(i-1) , 0 , 0 ) [0]" for i in 1:4],
            [" 1 ElmT4n3D 6 1 2 3 4"])
        @test read_diff(_w(tet)).mesh.blocks[1].msh == 4
        # a zero-element deck is legal and yields an empty block list
        f0 = read_diff(_w(_diff_deck(3,
            [" 1 ( 0 , 0 , 0 ) [0]", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"], String[])))
        @test isempty(f0.mesh.blocks)
        # duplicate node number warns and keeps the first coordinates
        dup = _diff_deck(4,
            [" 1 ( 0 , 0 , 0 ) [0]",
             " 2 ( 1 , 0 , 0 ) [0]",
             " 2 ( 9 , 9 , 9 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"]; maxn=3)
        fd = @test_logs (:warn, r"skipping duplicate") read_diff(_w(dup))
        @test fd.node_tags == UInt64[1, 2, 3]
        @test fd.mesh.coords[:, 2] == [1.0, 0.0, 0.0]
        @test fd.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3], 3, 1)
    end

    @testset "DIFF reader: failure paths" begin
        dir = mktempdir()
        function _w(text, name="b$(abs(hash(text))).diff")
            p = joinpath(dir, name)
            Base.write(p, text)
            return p
        end
        nodes3 = [" 1 ( 0 , 0 , 0 ) [0]", " 2 ( 1 , 0 , 0 ) [0]",
                  " 3 ( 0 , 1 , 0 ) [0]"]
        good = _diff_deck(3, nodes3, [" 1 ElmT3n2D 1 1 2 3"]; maxn=3)
        @test_throws ArgumentError read_diff(_w("nothing here\n"))
        @test_throws ArgumentError read_diff(_w(replace(good,
            "Number of elements   =" => "Number of elemxnts   =")))
        @test_throws ArgumentError read_diff(_w(replace(good,
            "Number of nodes      =" => "Number of nudes      =")))
        @test_throws ArgumentError read_diff(_w(replace(good,
            "Number of elements   =  1" => "Number of elements   =  -1")))
        @test_throws ArgumentError read_diff(_w(replace(good,
            "Number of nodes      =  3" => "Number of nodes      =  -3")))
        @test_throws ArgumentError read_diff(_w(replace(good,
            "Max number of nodes in an element: 3" =>
            "Max number of nodes in an element: -1")))
        # no `Boundary indicators:` anchor before the node section
        @test_throws ArgumentError read_diff(_w(replace(good,
            "Boundary indicators:" => "indicators:")))
        # malformed node rows
        @test_throws ArgumentError read_diff(_w(_diff_deck(3,
            [" 1 ( 0 , 0 , 0 )", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"])))                       # missing [k]
        @test_throws ArgumentError read_diff(_w(_diff_deck(3,
            [" 0 ( 0 , 0 , 0 ) [0]", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"])))                       # node number 0
        @test_throws ArgumentError read_diff(_w(_diff_deck(3,
            [" 1 ( zz , 0 , 0 ) [0]", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"])))                       # bad float
        @test_throws ArgumentError read_diff(_w(_diff_deck(3,
            [" 1 ( 0 , 0 , 0 ) [-1]", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"])))                       # k < 0
        @test_throws ArgumentError read_diff(_w(_diff_deck(3,
            [" 1 ( 0 , 0 , 0 ) [2] 7", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"])))                       # short tail
        @test_throws ArgumentError read_diff(_w(_diff_deck(3,
            [" 1 ( 0 , 0 , 0 ) [1] x", " 2 ( 1 , 0 , 0 ) [0]",
             " 3 ( 0 , 1 , 0 ) [0]"],
            [" 1 ElmT3n2D 1 1 2 3"])))                       # bad indicator
        # element records
        @test_throws ArgumentError read_diff(_w(_diff_deck(3, nodes3,
            [" 1 ElmZ9n3D 1 1 2 3"])))                       # unknown type
        @test_throws ArgumentError read_diff(_w(_diff_deck(3, nodes3,
            [" 1 ElmT3n2D -2 1 2 3"])))                      # negative subdomain
        @test_throws ArgumentError read_diff(_w(_diff_deck(3, nodes3,
            [" 1 ElmT3n2D 1 1 2 7"])))                       # dangling node
        @test_throws ArgumentError read_diff(_w(_diff_deck(3, nodes3,
            [" 1 ElmT3n2D 1 1 2"])))                         # short record
        @test_throws ArgumentError read_diff(_w(_diff_deck(3, nodes3,
            [" x ElmT3n2D 1 1 2 3"])))                       # bad element num
        # truncated decks: banner stops early, no node rows, no element rows
        head = split(good, "0 Boundary indicators:")[1]
        @test_throws ArgumentError read_diff(_w(head))
        nonodes = split(good, "#\n 1 ( ")[1] * "#\n"
        @test_throws ArgumentError read_diff(_w(nonodes))
        noelts = split(good, "#\n 1 ElmT3n2D")[1] * "#\n"
        @test_throws ArgumentError read_diff(_w(noelts))
        @test_throws ArgumentError read_diff(joinpath(dir, "missing.diff"))
        goodp = _w(good)
        @test_throws ArgumentError read_diff(goodp; max_nodes=2)
        @test_throws ArgumentError read_diff(goodp; max_elements=0)
        @test_throws ArgumentError read_diff(goodp; max_file_bytes=4)
        @test_throws ArgumentError read_diff(goodp; max_nodes=-1)
        @test_throws ArgumentError read_diff(3.0)
    end

    @testset "DIFF writer: GridFE deck" begin
        dir = mktempdir()
        coords = Float64[0 1 1 0 0.5; 0 0 1 1 0.5; 0 0 0 0 0]
        # MixedMesh validation: a cell's legacy `tags` entry must equal the
        # entity's first physical tag (0 when the entity has none).
        blocks = [ElementBlock(1, Int32[1 2; 2 3], Int32[7, 7]),
                  ElementBlock(2, Int32[1 2; 2 5; 4 1], Int32[5, 5]),
                  ElementBlock(1, reshape(Int32[4, 3], 2, 1), Int32[8])]
        data = MixedEntityData(Dict{Tuple{Int,Int},MixedEntity}(
            (1, 20) => MixedEntity(1, 20, (0., 0., 0., 1., 0., 0.);
                physical_tags=Int32[7]),
            (1, 21) => MixedEntity(1, 21, (0., 1., 0., 1., 1., 0.);
                physical_tags=Int32[8]),
            (2, 10) => MixedEntity(2, 10, (0., 0., 0., 1., 1., 0.);
                physical_tags=Int32[5]));
            node_entities=Tuple{Int,Int32}[(1, 20), (1, 20), (1, 21),
                                           (1, 21), (2, 10)],
            node_parametric=_F64V[nothing for _ in 1:5],
            external_node_tags=UInt64.(1:5),
            block_entities=[Int32[20, 20], Int32[10, 10], Int32[21]],
            external_element_tags=[UInt64[101, 102], UInt64[201, 202],
                                   UInt64[301]])
        mesh = MixedMesh(coords, blocks; entity_data=data,
            elementary_entities=[Int32[20, 20], Int32[10, 10], Int32[21]])
        path = joinpath(dir, "m.diff")
        @test write_diff(path, mesh) == path
        lines = readlines(path)
        @test any(==(" Finite element mesh (GridFE):"), lines)
        @test any(==(" Number of elements   =  2"), lines)   # top-dim only
        @test any(==(" Number of nodes      =  5"), lines)
        @test any(l -> occursin("All elements are of the same type: dpTRUE", l),
                  lines)
        @test any(l -> occursin("Max number of nodes in an element: 3", l),
                  lines)
        @test in(" 2 Boundary indicators:   7 8", lines)
        # curve 20 (physical 7) owns nodes 1,2,3; curve 21 (physical 8) owns
        # 3,4 — node 3 sits on both → [2] 7 8; node 5 is interior → [0].
        node_line(i) = lines[findfirst(
            l -> startswith(l, " $i ("), lines)]
        @test occursin(r"\[1\]  7\s*$", node_line(1))
        @test occursin(r"\[1\]  7\s*$", node_line(2))
        @test occursin(r"\[2\]  7  8\s*$", node_line(3))
        @test occursin(r"\[1\]  8\s*$", node_line(4))
        @test occursin(r"\[0\]\s*$", node_line(5))
        @test occursin(r"1\.0000000000000000E\+00", node_line(2))
        @test any(==("1 ElmT3n2D 5  1 2 4"), lines)
        @test any(==("2 ElmT3n2D 5  2 5 1"), lines)
        # reading the emitted deck reproduces connectivity, coords and the
        # subdomain column as the elementary tag (physicals are not carried)
        f = read_diff(path)
        @test f.mesh.coords == coords
        @test [b.msh for b in f.mesh.blocks] == [2]
        @test f.mesh.blocks[1].nodes == Int32[1 2; 2 5; 4 1]
        @test f.mesh.elementary_entities == [Int32[5, 5]]
        @test f.node_tags == UInt64.(1:5)
        @test f.element_tags == [UInt64[1, 2]]
        # no physical groups anywhere → subdomain/indicators are entity tags
        flat = joinpath(dir, "flat.diff")
        write_diff(flat, MixedMesh(coords, blocks;
            elementary_entities=[Int32[20, 20], Int32[10, 10], Int32[21]]))
        flines = readlines(flat)
        @test in(" 2 Boundary indicators:   20 21", flines)
        @test any(==("1 ElmT3n2D 10  1 2 4"), flines)
        @test read_diff(flat).mesh.elementary_entities == [Int32[10, 10]]
        # an entity with several physicals emits one record per membership
        multi = MixedEntityData(Dict{Tuple{Int,Int},MixedEntity}(
            (2, 10) => MixedEntity(2, 10, (0., 0., 0., 1., 1., 0.);
                physical_tags=Int32[5, 6]));
            node_entities=fill((2, 10), 5),
            node_parametric=_F64V[nothing for _ in 1:5],
            external_node_tags=UInt64.(1:5),
            block_entities=[Int32[10, 10]],
            external_element_tags=[UInt64[1, 2]])
        mp = joinpath(dir, "multi.diff")
        write_diff(mp, MixedMesh(coords,
            [ElementBlock(2, Int32[1 2; 2 5; 4 1], Int32[5, 5])];
            entity_data=multi))
        mlines = readlines(mp)
        @test in("1 ElmT3n2D 5  1 2 4", mlines)
        @test in("2 ElmT3n2D 6  1 2 4", mlines)
        @test in("3 ElmT3n2D 5  2 5 1", mlines)
        @test in("4 ElmT3n2D 6  2 5 1", mlines)
        # header still counts one cell per element (upstream `numElements`);
        # the reader therefore consumes the first two records only — the
        # write→read asymmetry is inherited from upstream
        @test any(==(" Number of elements   =  2"), mlines)
        fm = read_diff(mp)
        @test fm.mesh.elementary_entities == [Int32[5, 6]]
        # node 3 is used by no cell, so the writer omits it (upstream
        # `indexMeshVertices`) and the read back compacts labels 4/5
        @test fm.node_tags == UInt64[1, 2, 4, 5]
        @test fm.mesh.blocks[1].nodes == Int32[1 1; 2 2; 3 3]
        @test fm.mesh.coords == coords[:, [1, 2, 4, 5]]
    end

    @testset "DIFF writer: connectivity permutations and coverage" begin
        dir = mktempdir()
        coords = repeat(reshape(Float64.(1:8), 1, 8), 3, 1)
        # hex8 write placement (3,4,8,7,1,2,6,5) on dense labels
        hexpath = joinpath(dir, "hex.diff")
        write_diff(hexpath, MixedMesh(coords,
            [ElementBlock(5, reshape(Int32(1):Int32(8), 8, 1), Int32[0])];
            elementary_entities=[Int32[4]]))
        @test in("1 ElmB8n3D 4  3 4 8 7 1 2 6 5", readlines(hexpath))
        # quad4 write placement swaps slots 3/4 like upstream `ElmB4n2D`
        qpath = joinpath(dir, "quad.diff")
        write_diff(qpath, MixedMesh(coords[:, 1:4],
            [ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0])];
            elementary_entities=[Int32[7]]))
        @test in("1 ElmB4n2D 7  1 2 4 3", readlines(qpath))
        # only the top dimension is written; dim-1 cells mark boundary nodes
        mixed = MixedMesh(coords,
            [ElementBlock(5, reshape(Int32(1):Int32(8), 8, 1), Int32[0]),
             ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0]),
             ElementBlock(1, reshape(Int32[1, 2], 2, 1), Int32[0])];
            elementary_entities=[Int32[4], Int32[3], Int32[2]])
        mpath = joinpath(dir, "mixed.dim.diff")
        write_diff(mpath, mixed)
        mlines = readlines(mpath)
        @test in("1 ElmB8n3D 4  3 4 8 7 1 2 6 5", mlines)
        @test !any(l -> occursin("ElmB4n2D", l), mlines)
        @test in(" 1 Boundary indicators:   3", mlines)   # dim-1 entity tag
        @test any(l -> startswith(l, " 4 ("), mlines)
        # scaling_factor multiplies coordinates like Mesh.ScalingFactor
        sp = joinpath(dir, "scaled.diff")
        write_diff(sp, MixedMesh(coords[:, 1:4],
            [ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1), Int32[0])];
            elementary_entities=[Int32[7]]); scaling_factor=2.0)
        @test any(l -> occursin("2.0000000000000000E+00", l), readlines(sp))
        @test_throws ArgumentError write_diff(
            sp, MixedMesh(coords[:, 1:4],
                [ElementBlock(3, reshape(Int32[1, 2, 3, 4], 4, 1),
                                Int32[0])]); scaling_factor=Inf)
        # prisms have no DIFF type: the element section is empty like
        # upstream, but their nodes are still written
        prism = MixedMesh(Float64[0 1 0 0 1 0; 0 0 1 0 0 1; 0 0 0 1 1 1],
            [ElementBlock(6, reshape(Int32(1):Int32(6), 6, 1), Int32[0])];
            elementary_entities=[Int32[2]])
        pp = joinpath(dir, "prism.diff")
        write_diff(pp, prism)
        plines = readlines(pp)
        @test any(==(" Number of elements   =  0"), plines)
        @test any(==(" Number of nodes      =  6"), plines)
        @test in(" 0 Boundary indicators:  ", plines)
        @test_throws ArgumentError write_diff(
            joinpath(dir, "nodir", "x.diff"), prism)
    end

    @testset "DIFF writer round-trip: tet deck with physicals" begin
        dir = mktempdir()
        # tet4 + boundary tris: dim=3, indicators come from the dim-1 surface
        coords = Float64[0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0]
        data = MixedEntityData(Dict{Tuple{Int,Int},MixedEntity}(
            (3, 1) => MixedEntity(3, 1, (0., 0., 0., 1., 1., 1.);
                physical_tags=Int32[4]),
            (2, 7) => MixedEntity(2, 7, (0., 0., 0., 1., 1., 0.);
                physical_tags=Int32[9]));
            node_entities=Tuple{Int,Int32}[(2, 7), (2, 7), (2, 7), (2, 7),
                                           (3, 1), (3, 1)],
            node_parametric=_F64V[nothing for _ in 1:6],
            external_node_tags=UInt64.(1:6),
            block_entities=[Int32[1, 1], Int32[7, 7]],
            external_element_tags=[UInt64[1, 2], UInt64[3, 4]])
        mesh = MixedMesh(coords,
            [ElementBlock(4, Int32[1 1; 2 2; 3 5; 4 6], Int32[4, 4]),
             ElementBlock(2, Int32[1 1; 2 2; 3 4], Int32[9, 9])];
            entity_data=data,
            elementary_entities=[Int32[1, 1], Int32[7, 7]])
        path = joinpath(dir, "tet.diff")
        write_diff(path, mesh)
        lines = readlines(path)
        @test in(" 1 Boundary indicators:   9", lines)
        @test in("1 ElmT4n3D 4  1 2 3 4", lines)
        @test in("2 ElmT4n3D 4  1 2 5 6", lines)
        @test !any(l -> occursin("ElmT3n2D", l), lines)
        f = read_diff(path)
        @test [b.msh for b in f.mesh.blocks] == [4]
        @test f.mesh.elementary_entities == [Int32[4, 4]]
        @test f.mesh.coords == coords
    end

    @testset "session open/merge/write dispatch" begin
        _DIFF_API.initialize()
        try
            _DIFF_API.option("General.Terminal", 0)
            mktempdir() do dir
                deck = joinpath(dir, "q.diff")
                Base.write(deck, _diff_deck(4,
                    [" 1 ( 0 , 0 , 0 ) [0]", " 2 ( 1 , 0 , 0 ) [0]",
                     " 3 ( 1 , 1 , 0 ) [0]", " 4 ( 0 , 1 , 0 ) [0]"],
                    [" 1 ElmB4n2D 5 1 2 4 3"]))
                _DIFF_API.model.add("diffhost")
                @test _DIFF_API.open(deck) === nothing
                entities = _DIFF_API.model.get_entities()
                @test (2, 5) in entities          # subdomain column → entity
                @test Int32(3) in _DIFF_API.mesh.get_elements()[1]
                host = copy(entities)
                @test _DIFF_API.merge(deck) === nothing
                @test issubset(host, _DIFF_API.model.get_entities())
                # export routes through write_diff and re-imports
                out = joinpath(dir, "out.diff")
                @test _DIFF_API.write(out) === nothing
                f = read_diff(out)
                @test [b.msh for b in f.mesh.blocks] == [3]
                @test f.mesh.elementary_entities == [Int32[5, 5]]
                @test_throws ArgumentError _DIFF_API.open(
                    joinpath(dir, "bad.ext"))
            end
        finally
            _DIFF_API.finalize()
        end
    end
end
