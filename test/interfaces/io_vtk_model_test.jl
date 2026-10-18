# ── .vtk (legacy VTK) model-level mesh I/O ──
#
# Correctness  : hand-written fixtures follow the pinned `GModelIO_VTK.cpp`
#                token-stream layout (`# vtk DataFile Version 2.0`, title,
#                `ASCII`, `DATASET UNSTRUCTURED_GRID`, `POINTS n float`,
#                `CELLS m size` count+0-based-index records, `CELL_TYPES m`).
#                VTK type → MSH type decoding covers the fixed linear and
#                quadratic records plus complete-order Lagrange cells (68–74)
#                whose order is inferred from the node count; node ordering is
#                checked literal-for-literal against the inverse of the shared
#                `API._native_model_view_vtk_cell` write permutation (e.g.
#                VTK_QUADRATIC_TETRA swaps Gmsh tet10 nodes 9/10). `POLYDATA`
#                `POLYGONS` arity decode and `LINES` polyline→`lin2` splitting
#                with the upstream `centerline` physical group are covered.
# Robustness   : XML/BINARY dialects, truncated token streams, bad keywords,
#                count/arity mismatches, out-of-range connectivity, limits,
#                special blocks and unsupported (incomplete) cell types.
# Completeness : `API.open`/`API.merge` accept `.vtk`; `API.write` emits a
#                legacy ASCII `UNSTRUCTURED_GRID` whose physical tags land in
#                the upstream `SCALARS CellEntityIds` section (ignored on
#                read-back like upstream). Write→read is a literal round trip.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock, SpecialElementBlock
using Tessella.IO

const _VTK_API = Tessella.API

# One tri3 + one lin2 on four points.
const _VTK_TRI_LINE = """
# vtk DataFile Version 2.0
triangle and line
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0  1 0 0  1 1 0  0 1 0
CELLS 2 7
3 0 1 2
2 0 3
CELL_TYPES 2
5
3
"""

# Quadratic tet (VTK 24): connectivity 0..9 exercises the 9/10 node swap.
const _VTK_TET10 = """
# vtk DataFile Version 2.0
quadratic tet
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 10 float
0 0 0  1 0 0  0 1 0  0 0 1  .1 0 0  .6 .3 0  .2 0 .1  0 0 .6  .5 0 .5  0 .3 .3
CELLS 1 11
10 0 1 2 3 4 5 6 7 8 9
CELL_TYPES 1
24
"""

# Lagrange cells: curve(4 nodes), triangle(10), tetrahedron(20).
const _VTK_LAGRANGE = """
# vtk DataFile Version 2.0
lagrange cells
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 20 double
0 0 0 1 0 0 0 1 0 0 0 1 .33 0 0 .67 0 0 .67 .33 0 .33 .67 0 0 .67 0 0 0 .33
.33 0 .67 .67 0 .33 .25 0 .25 .33 .25 .25 .25 .33 .33 0 .25 .75 0 .25 .25
0 .75 .25 0 .5 .5 .5 .5 .25
CELLS 3 34
4 0 9 10 11
10 0 1 2 4 5 6 7 8 11 12
20 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
CELL_TYPES 3
68
69
71
"""

# POLYDATA with two POLYGONS (arity decode). Upstream reads a single section
# keyword per dataset, so a file carrying both POLYGONS and LINES keeps only
# the first — the two paths get separate fixtures below.
const _VTK_POLYDATA = """
# vtk DataFile Version 2.0
polydata
ASCII
DATASET POLYDATA
POINTS 5 float
0 0 0  1 0 0  1 1 0  0 1 0  .5 .5 0
POLYGONS 2 8
3 0 1 2
4 0 1 2 3
"""

# POLYDATA LINES: two polylines, each split into lin2 segments on its own
# curve entity under the upstream "centerline" physical group.
const _VTK_LINES = """
# vtk DataFile Version 2.0
polylines
ASCII
DATASET POLYDATA
POINTS 5 float
0 0 0  1 0 0  1 1 0  0 1 0  .5 .5 0
LINES 2 8
3 0 1 4
2 1 3
"""

@testset "VTK model mesh I/O" begin

    @testset "reader: UNSTRUCTURED_GRID basics" begin
        dir = mktempdir()
        path = joinpath(dir, "tri_line.vtk")
        Base.write(path, _VTK_TRI_LINE)
        f = read_vtk(path)
        @test f isa VtkFile
        @test f.node_tags == UInt64.(1:4)
        @test size(f.mesh.coords, 2) == 4
        @test f.mesh.coords[:, 2] == [1.0, 0.0, 0.0]
        @test [b.msh for b in f.mesh.blocks] == [2, 1]
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3], 3, 1)
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 4], 2, 1)
        # per-dimension entities: surface 1, curve 1
        @test f.mesh.elementary_entities == [Int32[1], Int32[1]]
        @test f.mesh.blocks[1].tags == Int32[0]
        @test f.element_tags == [UInt64[1], UInt64[2]]
        # point cell + vertex entity counter
        pf = joinpath(dir, "point.vtk")
        Base.write(pf, "# vtk DataFile Version 2.0\np\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 2 float\n0 0 0 1 0 0\n" *
            "CELLS 2 2\n1 0\n1 1\nCELL_TYPES 2\n1\n1\n")
        fp = read_vtk(pf)
        @test [b.msh for b in fp.mesh.blocks] == [15]
        @test fp.mesh.elementary_entities == [Int32[1, 2]]
    end

    @testset "reader: VTK→MSH node permutations" begin
        dir = mktempdir()
        path = joinpath(dir, "tet10.vtk")
        Base.write(path, _VTK_TET10)
        f = read_vtk(path)
        @test [b.msh for b in f.mesh.blocks] == [11]
        # inverse of the writer's tet10 order [1..8,10,9]
        @test f.mesh.blocks[1].nodes == reshape(
            Int32[1, 2, 3, 4, 5, 6, 7, 8, 10, 9], 10, 1)

        api_layout = Tessella.API._native_model_view_vtk_cell
        # every fixed read type uses the inverse of the shared write order
        for (vtk, msh) in (3 => 1, 21 => 8, 5 => 2, 22 => 9, 9 => 3, 23 => 16,
                           28 => 10, 10 => 4, 24 => 11, 12 => 5, 25 => 17,
                           29 => 12, 13 => 6, 26 => 18, 32 => 13, 14 => 7,
                           27 => 19, 1 => 15)
            order = api_layout(msh, "t").order
            conn = collect(Int32(0):Int32(length(order) - 1))
            points3 = join(("$i 0 0" for i in 0:length(order) - 1), " ")
            cell_line = string(length(order), " ", join(conn, " "))
            text = "# vtk DataFile Version 2.0\nx\nASCII\n" *
                "DATASET UNSTRUCTURED_GRID\n" *
                "POINTS $(length(order)) float\n$points3\n" *
                "CELLS 1 $(length(order) + 1)\n$cell_line\n" *
                "CELL_TYPES 1\n$vtk\n"
            p = joinpath(dir, "cell_$(vtk).vtk")
            Base.write(p, text)
            fr = read_vtk(p)
            @test [b.msh for b in fr.mesh.blocks] == [msh]
            nodes = fr.mesh.blocks[1].nodes[:, 1]
            @test all(nodes[order[k]] == conn[k] + 1 for k in eachindex(order))
        end
    end

    @testset "reader: Lagrange cell types 68-74" begin
        dir = mktempdir()
        path = joinpath(dir, "lagrange.vtk")
        Base.write(path, _VTK_LAGRANGE)
        f = read_vtk(path)
        @test [b.msh for b in f.mesh.blocks] == [26, 21, 29]
        api_layout = Tessella.API._native_model_view_vtk_cell
        conn_lin = Int32[1, 10, 11, 12]
        nodes = f.mesh.blocks[1].nodes[:, 1]
        order = api_layout(26, "t").order
        @test all(nodes[order[k]] == conn_lin[k] for k in eachindex(order))
        # a 9-node Lagrange triangle has no complete order -> loud reject
        bad = joinpath(dir, "lag_bad.vtk")
        Base.write(bad, "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 9 float\n" *
            join(("$i 0 0" for i in 0:8), " ") * "\n" *
            "CELLS 1 10\n9 0 1 2 3 4 5 6 7 8\nCELL_TYPES 1\n69\n")
        @test_throws ArgumentError read_vtk(bad)
    end

    @testset "reader: POLYDATA polygons and polylines" begin
        dir = mktempdir()
        path = joinpath(dir, "poly.vtk")
        Base.write(path, _VTK_POLYDATA)
        f = read_vtk(path)
        @test [b.msh for b in f.mesh.blocks] == [2, 3]
        @test f.mesh.blocks[1].nodes == reshape(Int32[1, 2, 3], 3, 1)
        @test f.mesh.blocks[2].nodes == reshape(Int32[1, 2, 3, 4], 4, 1)
        @test f.mesh.elementary_entities == [Int32[1], Int32[1]]

        lp = joinpath(dir, "lines.vtk")
        Base.write(lp, _VTK_LINES)
        fl = read_vtk(lp)
        @test [b.msh for b in fl.mesh.blocks] == [1]
        @test fl.mesh.blocks[1].nodes == Int32[1 2 2; 2 5 4]
        # one curve entity per polyline, physical 1 named "centerline"
        @test fl.mesh.elementary_entities[1] == Int32[1, 1, 2]
        @test fl.mesh.blocks[1].tags == Int32[1, 1, 1]
        @test fl.mesh.physical_names[(1, 1)] == "centerline"
    end

    @testset "reader: failure paths" begin
        dir = mktempdir()
        function _w(text, name = string(hash(text)))
            p = joinpath(dir, name * ".vtk")
            Base.write(p, text)
            return p
        end
        head = "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 4 float\n" *
            "0 0 0 1 0 0 1 1 0 0 1 0\n"
        @test_throws ArgumentError read_vtk(_w(""))
        @test_throws ArgumentError read_vtk(_w("<VTKFile type=\"UnstructuredGrid\">\n"))
        @test_throws ArgumentError read_vtk(_w("<?xml version=\"1.0\"?>\n"))
        @test_throws ArgumentError read_vtk(_w(
            "# vtk DataFile Version 2.0\nx\nBINARY\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 1 float\n0 0 0\n"))
        @test_throws ArgumentError read_vtk(_w(
            "# vtk DataFile Version 2.0\nx\nGZIP\n"))
        @test_throws ArgumentError read_vtk(_w(
            "# vtk DataFile Version 2.0\nx\nASCII\nDATASET STRUCTURED_GRID\n"))
        @test_throws ArgumentError read_vtk(_w(
            "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nVERTICES 1\n"))
        @test_throws ArgumentError read_vtk(_w(                       # no points
            "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 0 float\n"))
        @test_throws ArgumentError read_vtk(_w(                       # bad point type
            "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 1 int\n0 0 0\n"))
        @test_throws ArgumentError read_vtk(_w(                       # truncated coords
            "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 4 float\n0 0 0\n"))
        @test_throws ArgumentError read_vtk(_w(                       # bad float
            "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 1 float\n0 zz 0\n"))
        @test_throws ArgumentError read_vtk(_w(                       # non-finite
            "# vtk DataFile Version 2.0\nx\nASCII\n" *
            "DATASET UNSTRUCTURED_GRID\nPOINTS 1 float\n0 NaN 0\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # no cells
            "CELL_TYPES 0\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # truncated cell
            "CELLS 1 3\n3 0 1\nCELL_TYPES 1\n5\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # index range
            "CELLS 1 4\n3 0 1 9\nCELL_TYPES 1\n5\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # negative index
            "CELLS 1 4\n3 0 -1 2\nCELL_TYPES 1\n5\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # missing CELL_TYPES
            "CELLS 1 4\n3 0 1 2\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # type count
            "CELLS 1 4\n3 0 1 2\nCELL_TYPES 2\n5\n5\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # unsupported type
            "CELLS 1 3\n2 0 1\nCELL_TYPES 1\n2\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # arity mismatch
            "CELLS 1 5\n4 0 1 2 3\nCELL_TYPES 1\n5\n"))
        @test_throws ArgumentError read_vtk(_w(head *                  # lagrange unknown
            "CELLS 1 3\n2 0 1\nCELL_TYPES 1\n70\n"))  # quad lagrange w/ 2 nodes
        @test_throws ArgumentError read_vtk(joinpath(dir, "missing.vtk"))
        @test_throws ArgumentError read_vtk(_w(_VTK_TRI_LINE); max_nodes = 3)
        @test_throws ArgumentError read_vtk(_w(_VTK_TRI_LINE); max_elements = 1)
        @test_throws ArgumentError read_vtk(_w(_VTK_TRI_LINE); max_file_bytes = 4)
        # trailing CELL_DATA/POINT_DATA sections are ignored like upstream
        tagged = _w(head * "CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n" *
            "CELL_DATA 1\nSCALARS CellEntityIds int 1\nLOOKUP_TABLE default\n7\n" *
            "POINT_DATA 4\nSCALARS junk float 1\nLOOKUP_TABLE default\n" *
            "0 1 2 3\n")
        ft = read_vtk(tagged)
        @test ft.mesh.blocks[1].tags == Int32[0]
    end

    @testset "writer: file layout and round trip" begin
        dir = mktempdir()
        coords = zeros(3, 30)
        coords[1, :] = collect(0.0:29.0)
        blocks = ElementBlock[
            ElementBlock(1, reshape(Int32[1, 2], 2, 1)),
            ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[7]),
            ElementBlock(11, reshape(Int32[1:10;], 10, 1)),
            ElementBlock(17, reshape(Int32[1:20;], 20, 1)),
            ElementBlock(13, reshape(Int32[1:18;], 18, 1)),
            ElementBlock(19, reshape(Int32[1:13;], 13, 1)),
            ElementBlock(21, reshape(Int32[1:10;], 10, 1)),
            ElementBlock(26, reshape(Int32[1:4;], 4, 1))]
        mesh = MixedMesh(coords, blocks)
        path = joinpath(dir, "out.vtk")
        @test write_vtk(path, mesh) == path
        lines = readlines(path)
        @test lines[1] == "# vtk DataFile Version 2.0"
        @test lines[3] == "ASCII"
        @test lines[4] == "DATASET UNSTRUCTURED_GRID"
        @test lines[5] == "POINTS 30 double"
        i = findfirst(startswith("CELLS "), lines)
        @test lines[i] == "CELLS 8 88"
        # tet10 emitted in VTK order: Gmsh nodes 9/10 swapped
        tet_line = lines[i + 3]
        @test tet_line == "10 0 1 2 3 4 5 6 7 9 8"
        j = findfirst(startswith("CELL_TYPES "), lines)
        @test parse.(Int, lines[j + 1:j + 8]) == [3, 5, 24, 25, 32, 27, 69, 68]
        # physical tag -> CellEntityIds section
        @test "SCALARS CellEntityIds int 1" in lines
        k = findfirst(==("LOOKUP_TABLE default"), lines)
        @test parse.(Int, lines[k + 1:k + 8]) == [-1, 7, -1, -1, -1, -1, -1, -1]

        # write->read round trip is literal for connectivity and types
        f2 = read_vtk(path)
        @test f2.mesh.coords == coords
        @test [b.msh for b in f2.mesh.blocks] == [b.msh for b in blocks]
        for (b1, b2) in zip(blocks, f2.mesh.blocks)
            @test b1.nodes == b2.nodes
        end
        # tags are reset: CellEntityIds is ignored on read like upstream
        @test all(all(==(0), b.tags) for b in f2.mesh.blocks)

        # writer rejects unsupported inputs before touching the target
        @test_throws ArgumentError write_vtk(joinpath(dir, "s.vtk"),
            MixedMesh(coords,
                [SpecialElementBlock(67, reshape(Int32[1, 2], 2, 1))]))
        @test_throws ArgumentError write_vtk(joinpath(dir, "i.vtk"),
            MixedMesh(coords,
                [ElementBlock(20, reshape(Int32[1:9;], 9, 1))]))
        @test !isfile(joinpath(dir, "i.vtk"))
        @test_throws ArgumentError write_vtk(
            joinpath(dir, "nodir", "x.vtk"), mesh)
        @test_throws ArgumentError write_vtk(joinpath(dir, "t.vtk"), mesh;
            title = "bad\ntitle")
    end

    @testset "session open/merge/write dispatch" begin
        _VTK_API.initialize()
        try
            _VTK_API.option("General.Terminal", 0)
            mktempdir() do dir
                vtk = joinpath(dir, "tri_line.vtk")
                Base.write(vtk, _VTK_TRI_LINE)
                _VTK_API.model.add("vtkhost")
                @test _VTK_API.open(vtk) === nothing
                entities = _VTK_API.model.get_entities()
                @test (1, 1) in entities && (2, 1) in entities
                types = _VTK_API.mesh.get_elements()[1]
                @test Int32(2) in types && Int32(1) in types
                host = copy(entities)
                @test _VTK_API.merge(vtk) === nothing
                merged = _VTK_API.model.get_entities()
                @test issubset(host, merged)
                out = joinpath(dir, "session.vtk")
                @test _VTK_API.write(out) === nothing
                @test isfile(out)
                @test readlines(out)[1] == "# vtk DataFile Version 2.0"
                # the exported file opens again through the same dispatch
                @test _VTK_API.open(out) === nothing
                # XML .vtu is out of scope for the model filename APIs
                @test_throws ArgumentError _VTK_API.open(
                    joinpath(dir, "x.vtu"))
                @test_throws ArgumentError _VTK_API.write(
                    joinpath(dir, "x.vtu"))
                @test !isfile(joinpath(dir, "x.vtu"))
            end
        finally
            _VTK_API.finalize()
        end
    end
end
