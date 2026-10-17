# ── .unv (I-DEAS Universal) and .mesh (INRIA/Medit) mesh I/O ──────────────────
#
# Correctness  : hand-written fixtures transcribe the documented file layouts
#                (UNV datasets 2411/2412; Medit MeshVersionFormatted sections);
#                a real box mesh generated through `API.mesh.generate` is
#                exported with `API.write` and re-imported with `API.open` for
#                both extensions, checking node/element/physical parity.
# Robustness   : `D` exponents, malformed records, unknown descriptors
#                (warn + skip, upstream behaviour), dangling node references,
#                empty sections, binary-version rejection, duplicate labels,
#                and unrepresentable element types on write.
# Completeness : UNV elementary/physical fields and Medit `ref` tags are
#                carried end-to-end; `API.merge` renumbers file-local labels
#                into fresh session UIDs and preserves the existing model.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock, SpecialElementBlock
using Tessella.IO

const _UNV_MESH_API=Tessella.API

# Hand fixture: non-contiguous labels, `D` exponents, a beam-descriptor record
# (21 → 2-node line) with its orientation line, a triangle and a tetrahedron.
const _UNV_FIXTURE = """
    -1
  2411
         5         1         1        11
   0.0000000000000000D+00   0.0D+00   0.0D+00
         7         1         1        11
   1.0D+00   0.0D+00   0.0D+00
         9         1         1        11
   0.0D+00   1.0D+00   0.0D+00
        11         1         1        11
   0.0D+00   0.0D+00   1.0D+00
    -1
    -1
  2412
        20        21        42         7         3         2
         0         0         0
         5         7
        30        41        43         8         3         3
         5         7         9
        40       111        44         9         3         4
         5         7         9        11
    -1
"""

const _MESH_FIXTURE = """
MeshVersionFormatted 1
 Dimension
 3
 Vertices
 4
 0. 0. 0. 0
 1. 0. 0. 0
 0. 1. 0. 0
 0. 0. 1e0 0
 Edges
 1
 1 4 7
 Triangles
 2
 1 2 3 8
 1 3 4 8
 Quadrilaterals
 0
 Tetrahedra
 1
 1 2 3 4 9
 Hexahedra
 0
 Corners
 1
 1
 RequiredVertices
 2
 1
 4
 End
"""

# Sorted connectivity sets per MSH family, keyed by node UID.
function _cells_by_type(get_elements)
    types, tags, nodes = get_elements
    out = Dict{Int32, Vector{Vector{UInt64}}}()
    for i in eachindex(types)
        out[types[i]] = sort!([sort!(vec(collect(c))) for c in nodes[i]])
    end
    return out
end

@testset "UNV/Medit mesh I/O" begin

    @testset "UNV reader: datasets, D exponents, tags" begin
        dir = mktempdir()
        path = joinpath(dir, "fixture.unv")
        Base.write(path, _UNV_FIXTURE)
        f = read_unv(path)
        @test f isa UnvFile
        @test f.node_tags == UInt64[5, 7, 9, 11]
        @test size(f.mesh.coords, 2) == 4
        @test f.mesh.coords[:, 2] == [1.0, 0.0, 0.0]   # 1.0D+00
        @test [b.msh for b in f.mesh.blocks] == [1, 2, 4]
        @test [b.nodes for b in f.mesh.blocks] ==
              [reshape(Int32[1, 2], 2, 1), reshape(Int32[1, 2, 3], 3, 1),
               reshape(Int32[1, 2, 3, 4], 4, 1)]
        @test [b.tags for b in f.mesh.blocks] == [Int32[7], Int32[8], Int32[9]]
        @test f.mesh.elementary_entities == [Int32[42], Int32[43], Int32[44]]
        @test f.element_tags == [UInt64[20], UInt64[30], UInt64[40]]
        # Round-trip preserves every channel.
        rt = joinpath(dir, "rt.unv")
        write_unv(rt, f.mesh)
        f2 = read_unv(rt)
        @test f2.mesh.coords == f.mesh.coords
        @test [b.msh for b in f2.mesh.blocks] == [1, 2, 4]
        @test [b.nodes for b in f2.mesh.blocks] == [b.nodes for b in f.mesh.blocks]
        @test [b.tags for b in f2.mesh.blocks] == [b.tags for b in f.mesh.blocks]
        @test f2.mesh.elementary_entities == f.mesh.elementary_entities
    end

    @testset "UNV reader: descriptor aliases and orders" begin
        dir = mktempdir()
        # 6-node triangle via descriptor 92 (file order differs from Gmsh order)
        path = joinpath(dir, "tri6.unv")
        Base.write(path, """
    -1
  2411
$(join(["        $i         1         1        11\n   $(Float64(i)) 0.0 0.0" for i in 1:6], "\n"))
    -1
    -1
  2412
         1        92        10         5         7         6
         1         2         3         4         5         6
    -1
""")
        f = read_unv(path)
        @test only(f.mesh.blocks).msh == 9
        # descriptor 92 permutes (1,3,5,2,4,6): file labels 1..6 map to
        # rows (1,3,5,2,4,6)
        @test only(f.mesh.blocks).nodes == reshape(Int32[1, 3, 5, 2, 4, 6], 6, 1)
        # unknown descriptor warns once and skips; the file still parses
        warn_path = joinpath(dir, "warn.unv")
        Base.write(warn_path, replace(_UNV_FIXTURE, "       111" => "       999"))
        w = read_unv(warn_path)
        @test [b.msh for b in w.mesh.blocks] == [1, 2]
        @test length(w.element_tags) == 2
    end

    @testset "UNV reader: failure paths" begin
        dir = mktempdir()
        function _w(text)
            p = joinpath(dir, string(hash(text), ".unv"))
            Base.write(p, text)
            return p
        end
        # malformed 2411 record (missing a field)
        @test_throws ArgumentError read_unv(_w("""
    -1
  2411
         1         1         1
   0.0 0.0 0.0
    -1
"""))
        # dangling node reference
        @test_throws ArgumentError read_unv(_w("""
    -1
  2411
         1         1         1        11
   0.0 0.0 0.0
    -1
    -1
  2412
         1        21         1         0         7         2
         0         0         0
         1         9
    -1
"""))
        # truncated element node list (file ends mid-record)
        @test_throws ArgumentError read_unv(_w("""
    -1
  2411
         1         1         1        11
   0.0 0.0 0.0
         2         1         1        11
   1.0 0.0 0.0
    -1
    -1
  2412
         1        21         1         0         7         2
         0         0         0
         1
"""))
        # duplicate node label
        @test_throws ArgumentError read_unv(_w("""
    -1
  2411
         1         1         1        11
   0.0 0.0 0.0
         1         1         1        11
   1.0 0.0 0.0
    -1
"""))
        # node count mismatch for the descriptor
        @test_throws ArgumentError read_unv(_w("""
    -1
  2411
         1         1         1        11
   0.0 0.0 0.0
         2         1         1        11
   1.0 0.0 0.0
         3         1         1        11
   0.0 1.0 0.0
         4         1         1        11
   1.0 1.0 0.0
    -1
    -1
  2412
         1       111         1         0         7         4
         1         2         3
    -1
"""))
        @test_throws ArgumentError read_unv(joinpath(dir, "missing.unv"))
        @test_throws ArgumentError read_unv(_w("""
    -1
  2411
         1         1         1        11
   0.0 0.0 0.0
    -1
"""); max_nodes=0)
    end

    @testset "Medit reader: sections, refs, leniency" begin
        dir = mktempdir()
        path = joinpath(dir, "fixture.mesh")
        Base.write(path, _MESH_FIXTURE)
        f = read_mesh(path)
        @test f isa MeshFile
        @test size(f.mesh.coords, 2) == 4
        @test f.mesh.coords[:, 4] == [0.0, 0.0, 1.0]   # 1e0
        @test [b.msh for b in f.mesh.blocks] == [1, 2, 4]
        @test [b.nodes for b in f.mesh.blocks] ==
              [reshape(Int32[1, 4], 2, 1), Int32[1 1; 2 3; 3 4],
               reshape(Int32[1, 2, 3, 4], 4, 1)]
        # each `ref` becomes both the elementary entity tag and physical tag
        @test [b.tags for b in f.mesh.blocks] == [Int32[7], Int32[8, 8], Int32[9]]
        @test f.mesh.elementary_entities == [Int32[7], Int32[8, 8], Int32[9]]
        # round-trip is symmetric
        rt = joinpath(dir, "rt.mesh")
        write_mesh(rt, f.mesh)
        f2 = read_mesh(rt)
        @test f2.mesh.coords == f.mesh.coords
        @test [b.msh for b in f2.mesh.blocks] == [1, 2, 4]
        @test [b.nodes for b in f2.mesh.blocks] == [b.nodes for b in f.mesh.blocks]
        @test [b.tags for b in f2.mesh.blocks] == [b.tags for b in f.mesh.blocks]
        # inline `Dimension 3` and inline section counts are accepted
        inl = joinpath(dir, "inline.mesh")
        Base.write(inl, "MeshVersionFormatted 2\nDimension 3\n" *
            "Vertices 2\n0 0 0 1\n1 0 0 2\nEdges 1\n1 2 3\nEnd\n")
        fi = read_mesh(inl)
        @test size(fi.mesh.coords, 2) == 2
        @test only(fi.mesh.blocks).msh == 1
        @test only(fi.mesh.blocks).tags == Int32[3]
    end

    @testset "Medit reader: failure paths" begin
        dir = mktempdir()
        function _wm(text)
            p = joinpath(dir, string(hash(text), ".mesh"))
            Base.write(p, text)
            return p
        end
        @test_throws ArgumentError read_mesh(_wm(
            "Vertices\n 1\n 0 0 0 0\n End\n"))                       # no header
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 3\n End\n"))                      # binary marker
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 9\n End\n"))                      # unknown version
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 2\n Dimension\n 4\n End\n"))      # bad dimension
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 2\n Dimension 3\n Vertices\n 1\n" *
            " 0 0 0 0\n Edges\n 1\n 1 9 3\n"))                      # dangling ref
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 2\n Dimension 3\n Vertices\n 1\n" *
            " 0 0 zz 0\n"))                                         # bad float
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 2\n Dimension 3\n Vertices\n 1\n" *
            " 0 0 0 0\n Triangles\n 1\n 1 1 1 -2\n"))               # negative ref
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 2\n Dimension 3\n Vertices\n 0\n" *
            " Edges\n 2\n 1 1 1\n 1 1 1\n"))                        # before Vertices
        @test_throws ArgumentError read_mesh(_wm(
            "MeshVersionFormatted 2\n Dimension 3\n Vertices\n 1\n" *
            " 0 0 0 0\n Edges\n 1\n 1 1 3\n"); max_elements=0)
        @test_throws ArgumentError read_mesh(joinpath(dir, "missing.mesh"))
    end

    @testset "writers reject unrepresentable blocks" begin
        dir = mktempdir()
        coords = zeros(3, 9)
        coords[1, :] = [0, 1, 1, 0, 0.5, 1, 0.5, 0, 0.5]
        quad9 = MixedMesh(coords,
            [ElementBlock(10, reshape(Int32[1:9;], 9, 1), Int32[1])])
        @test_throws ArgumentError write_unv(joinpath(dir, "q.unv"), quad9)
        @test_throws ArgumentError write_mesh(joinpath(dir, "q.mesh"), quad9)
        special = MixedMesh(zeros(3, 2),
            [SpecialElementBlock(67, reshape(Int32[1, 2], 2, 1))])
        @test_throws ArgumentError write_unv(joinpath(dir, "s.unv"), special)
        @test_throws ArgumentError write_mesh(joinpath(dir, "s.mesh"), special)
        empty = MixedMesh(zeros(3, 0), ElementBlock[])
        ep = joinpath(dir, "empty.unv")
        @test write_unv(ep, empty) == ep
        @test read_unv(ep).node_tags == UInt64[]
        em = joinpath(dir, "empty.mesh")
        @test write_mesh(em, empty) == em
        @test size(read_mesh(em).mesh.coords, 2) == 0
        @test_throws ArgumentError write_unv(
            joinpath(dir, "nodir", "x.unv"), empty)
        @test_throws ArgumentError write_mesh(
            joinpath(dir, "nodir", "x.mesh"), empty)
    end

    @testset "session round-trip through API.write/API.open" begin
        _UNV_MESH_API.initialize()
        try
            _UNV_MESH_API.option("General.Terminal", 0)
            _UNV_MESH_API.option("Mesh.MeshSizeFactor", 8.0)
            _UNV_MESH_API.option("Mesh.SaveAll", 1.0)
            _UNV_MESH_API.model.add("unvbox")
            _UNV_MESH_API.model.add_box(0, 0, 0, 1, 1, 1; tag=1)
            _UNV_MESH_API.model.add_physical_group(3, [1]; tag=5, name="vol")
            _UNV_MESH_API.model.add_physical_group(2, [1, 2]; tag=6, name="side")
            generated = _UNV_MESH_API.mesh.generate(3)
            @test size(generated.coords, 2) > 0
            source_nodes = _UNV_MESH_API.mesh.get_nodes()
            source_cells = _cells_by_type(_UNV_MESH_API.mesh.get_elements())
            @test 4 in keys(source_cells)        # tets present
            @test 2 in keys(source_cells)        # boundary triangles present
            @test 1 in keys(source_cells)        # boundary lines present
            source_groups = _UNV_MESH_API.model.get_physical_groups()
            @test (2, 6) in source_groups && (3, 5) in source_groups

            mktempdir() do dir
                for ext in ("unv", "mesh")
                    path = joinpath(dir, "box." * ext)
                    @test _UNV_MESH_API.write(path) === nothing
                    @test isfile(path)
                    _UNV_MESH_API.model.add("re" * ext)
                    @test _UNV_MESH_API.open(path) === nothing
                    @test "box" in _UNV_MESH_API.model.list()
                    re_nodes = _UNV_MESH_API.mesh.get_nodes()
                    @test re_nodes[1] == source_nodes[1]          # node UIDs
                    @test re_nodes[2] == source_nodes[2]          # coordinates
                    re_cells = _cells_by_type(_UNV_MESH_API.mesh.get_elements())
                    for msh in (Int32(1), Int32(2), Int32(4))
                        @test re_cells[msh] == source_cells[msh]
                    end
                    # point cells are not representable in either format and
                    # are dropped (descriptor-0 warn+skip / Medit omission)
                    @test !haskey(re_cells, Int32(15))
                    @test (2, 6) in _UNV_MESH_API.model.get_physical_groups()
                    @test (3, 5) in _UNV_MESH_API.model.get_physical_groups()
                end
            end
        finally
            _UNV_MESH_API.finalize()
        end
    end

    @testset "session merge keeps the host model intact" begin
        _UNV_MESH_API.initialize()
        try
            _UNV_MESH_API.option("General.Terminal", 0)
            _UNV_MESH_API.model.add("host")
            _UNV_MESH_API.model.add_discrete_entity(1, 7)
            _UNV_MESH_API.model.add_discrete_entity(2, 9)
            _UNV_MESH_API.mesh.add_nodes(1, 7, [1, 2], [0.0, 0, 0, 1, 0, 0])
            _UNV_MESH_API.mesh.add_nodes(2, 9, [3, 4], [0.5, 0, 0, 0.5, 0.5, 0])
            _UNV_MESH_API.mesh.add_elements_by_type(7, 1, [11], [1, 2])
            _UNV_MESH_API.mesh.add_elements_by_type(9, 2, [12], [2, 3, 4])
            host_entities = _UNV_MESH_API.model.get_entities()
            @test (1, 7) in host_entities && (2, 9) in host_entities

            mktempdir() do dir
                unv = joinpath(dir, "add.unv")
                Base.write(unv, """
    -1
  2411
       101         1         1        11
   0.0D+00   0.0   0.0
       102         1         1        11
   2.0D+00   0.0   0.0
    -1
    -1
  2412
        55        21        42         3         7         2
         0         0         0
       101       102
    -1
""")
                @test _UNV_MESH_API.merge(unv) === nothing
                after_unv = _UNV_MESH_API.model.get_entities()
                @test issubset(host_entities, after_unv)
                @test (1, 42) in after_unv          # UNV elementary field

                meshp = joinpath(dir, "add.mesh")
                Base.write(meshp, """
MeshVersionFormatted 2
 Dimension
 3
 Vertices
 3
 0 0 0 30
 1 0 0 30
 0 1 0 30
 Triangles
 1
 1 2 3 30
 End
""")
                @test _UNV_MESH_API.merge(meshp) === nothing
                after_mesh = _UNV_MESH_API.model.get_entities()
                @test issubset(after_unv, after_mesh)
                @test (2, 30) in after_mesh         # Medit ref → entity tag
                @test (2, 30) in _UNV_MESH_API.model.get_physical_groups()
                types = _UNV_MESH_API.mesh.get_elements()[1]
                @test Int32(1) in types && Int32(2) in types
            end
        finally
            _UNV_MESH_API.finalize()
        end
    end

    @testset "session-level dispatch guards" begin
        _UNV_MESH_API.initialize()
        try
            _UNV_MESH_API.option("General.Terminal", 0)
            _UNV_MESH_API.model.add("guards")
            mktempdir() do dir
                for ext in ("unv", "mesh")
                    missing = joinpath(dir, "missing." * ext)
                    @test_throws ArgumentError _UNV_MESH_API.open(missing)
                    @test_throws ArgumentError _UNV_MESH_API.merge(missing)
                end
                bad = joinpath(dir, "x.vtk")
                @test_throws ArgumentError _UNV_MESH_API.write(bad)
            end
        finally
            _UNV_MESH_API.finalize()
        end
    end
end
