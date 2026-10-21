# ── .meshb (INRIA/Medit binary GMF) mesh I/O ─────────────────────────────
#
# Correctness  : embedded base64 fixtures follow the libMeshb `GmfKwdCod`
#                layout (keyword code, next-keyword offset, count, records)
#                and decode under every storage version — v1: f32/i32/i32,
#                v2: f64/i32/i32, v3: f64/i32/i64-offsets, v4: i64 fields —
#                in little- and big-endian form; a meshio-produced file and
#                Dimension 2 (z filled with 0) are covered, and the same
#                mesh round-trips as ASCII `.mesh` without drift.
# Robustness   : bad magic/version/dimension, a missing `Dimension`
#                keyword, elements before `Vertices`, dangling vertex
#                references, negative refs, inconsistent next-keyword
#                offsets, unskippable keywords, truncation, file/node/
#                element limits, and unrepresentable types or widths on
#                write.
# Completeness : trailing `ref` fields become elementary+physical tags,
#                TetrahedraP2 permutes like the ASCII reader, unrecognized
#                keywords (`Corners`, `RequiredVertices`, ...) skip through
#                their offsets, and `API.open`/`API.merge`/`API.write` route
#                `.meshb` through the same foreign-mesh channel as `.mesh`.

using Test
import Tessella
using Tessella.Elements: MixedMesh, ElementBlock, SpecialElementBlock
using Tessella.IO
using Base64: base64decode

const _MESHB_API=Tessella.API

# libMeshb-layout cube fixture (verified against meshio's GMF reader): 8
# vertices, one record each of Edges/Triangles/Quadrilaterals/Tetrahedra/
# Hexahedra with refs 7/8/9/11/12, plus `Corners` and `RequiredVertices`
# sections that must be skipped through their next-keyword offsets.
const _FX_V3 = base64decode(
    "AQAAAAMAAAADAAAAGAAAAAAAAAADAAAABAAAAAgBAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAADwPwAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAA8D8AAAAA" *
    "AADwPwAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAAAIAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8EAAAAAAAAAAAA8D8AAAAAAAAAAAAAAAAAAPA/BAAAAAAAAAAAAPA/" *
    "AAAAAAAA8D8AAAAAAADwPwQAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAA8D8EAAAABQAAACQB" *
    "AAAAAAAAAQAAAAEAAAACAAAABwAAAAYAAABEAQAAAAAAAAEAAAABAAAAAgAAAAMAAAAIAAAA" *
    "BwAAAGgBAAAAAAAAAQAAAAEAAAACAAAAAwAAAAQAAAAJAAAACAAAAIwBAAAAAAAAAQAAAAEA" *
    "AAACAAAAAwAAAAUAAAALAAAACgAAAMABAAAAAAAAAQAAAAEAAAACAAAAAwAAAAQAAAAFAAAA" *
    "BgAAAAcAAAAIAAAADAAAAA0AAADYAQAAAAAAAAIAAAABAAAABAAAAA8AAADsAQAAAAAAAAEA" *
    "AAACAAAANgAAAAAAAAAAAAAA")

const _FX_V1 = base64decode(
    "AQAAAAEAAAADAAAAFAAAAAMAAAAEAAAAoAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAACAPwAA" *
    "AAAAAAAAAgAAAAAAgD8AAIA/AAAAAAIAAAAAAAAAAACAPwAAAAACAAAAAAAAAAAAAAAAAIA/" *
    "BAAAAAAAgD8AAAAAAACAPwQAAAAAAIA/AACAPwAAgD8EAAAAAAAAAAAAgD8AAIA/BAAAAAUA" *
    "AAC4AAAAAQAAAAEAAAACAAAABwAAAAYAAADUAAAAAQAAAAEAAAACAAAAAwAAAAgAAAAHAAAA" *
    "9AAAAAEAAAABAAAAAgAAAAMAAAAEAAAACQAAAAgAAAAUAQAAAQAAAAEAAAACAAAAAwAAAAUA" *
    "AAALAAAACgAAAEQBAAABAAAAAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgAAAAMAAAA" *
    "DQAAAFgBAAACAAAAAQAAAAQAAAAPAAAAaAEAAAEAAAACAAAANgAAAAAAAAA=")

const _FX_V2 = base64decode(
    "AQAAAAIAAAADAAAAFAAAAAMAAAAEAAAAAAEAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8AAAAAAAAAAAAAAAAAAAAAAgAAAAAAAAAAAPA/AAAAAAAA8D8AAAAA" *
    "AAAAAAIAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAPA/BAAAAAAAAAAAAPA/AAAAAAAAAAAAAAAAAADwPwQAAAAAAAAAAADwPwAAAAAAAPA/" *
    "AAAAAAAA8D8EAAAAAAAAAAAAAAAAAAAAAADwPwAAAAAAAPA/BAAAAAUAAAAYAQAAAQAAAAEA" *
    "AAACAAAABwAAAAYAAAA0AQAAAQAAAAEAAAACAAAAAwAAAAgAAAAHAAAAVAEAAAEAAAABAAAA" *
    "AgAAAAMAAAAEAAAACQAAAAgAAAB0AQAAAQAAAAEAAAACAAAAAwAAAAUAAAALAAAACgAAAKQB" *
    "AAABAAAAAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgAAAAMAAAANgAAAAAAAAA=")

const _FX_V4 = base64decode(
    "AQAAAAQAAAADAAAAGAAAAAAAAAADAAAABAAAACwBAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAA" *
    "AAAAAAAA8D8AAAAAAADwPwAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAA" *
    "AAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADwPwQAAAAAAAAAAAAAAAAA8D8AAAAA" *
    "AAAAAAAAAAAAAPA/BAAAAAAAAAAAAAAAAADwPwAAAAAAAPA/AAAAAAAA8D8EAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8AAAAAAADwPwQAAAAAAAAABQAAAFgBAAAAAAAAAQAAAAAAAAABAAAA" *
    "AAAAAAIAAAAAAAAABwAAAAAAAAAGAAAAjAEAAAAAAAABAAAAAAAAAAEAAAAAAAAAAgAAAAAA" *
    "AAADAAAAAAAAAAgAAAAAAAAABwAAAMgBAAAAAAAAAQAAAAAAAAABAAAAAAAAAAIAAAAAAAAA" *
    "AwAAAAAAAAAEAAAAAAAAAAkAAAAAAAAACAAAAAQCAAAAAAAAAQAAAAAAAAABAAAAAAAAAAIA" *
    "AAAAAAAAAwAAAAAAAAAFAAAAAAAAAAsAAAAAAAAACgAAAGACAAAAAAAAAQAAAAAAAAABAAAA" *
    "AAAAAAIAAAAAAAAAAwAAAAAAAAAEAAAAAAAAAAUAAAAAAAAABgAAAAAAAAAHAAAAAAAAAAgA" *
    "AAAAAAAADAAAAAAAAAANAAAAhAIAAAAAAAACAAAAAAAAAAEAAAAAAAAABAAAAAAAAAAPAAAA" *
    "oAIAAAAAAAABAAAAAAAAAAIAAAAAAAAANgAAAAAAAAAAAAAA")

const _FX_V3BE = base64decode(
    "AAAAAQAAAAMAAAADAAAAAAAAABgAAAADAAAABAAAAAAAAAEIAAAACAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAA/8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACP/AAAAAAAAA/8AAA" *
    "AAAAAAAAAAAAAAAAAAAAAgAAAAAAAAAAP/AAAAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAAAA" *
    "AAAAAAAAP/AAAAAAAAAAAAAEP/AAAAAAAAAAAAAAAAAAAD/wAAAAAAAAAAAABD/wAAAAAAAA" *
    "P/AAAAAAAAA/8AAAAAAAAAAAAAQAAAAAAAAAAD/wAAAAAAAAP/AAAAAAAAAAAAAEAAAABQAA" *
    "AAAAAAEkAAAAAQAAAAEAAAACAAAABwAAAAYAAAAAAAABRAAAAAEAAAABAAAAAgAAAAMAAAAI" *
    "AAAABwAAAAAAAAFoAAAAAQAAAAEAAAACAAAAAwAAAAQAAAAJAAAACAAAAAAAAAGMAAAAAQAA" *
    "AAEAAAACAAAAAwAAAAUAAAALAAAACgAAAAAAAAHAAAAAAQAAAAEAAAACAAAAAwAAAAQAAAAF" *
    "AAAABgAAAAcAAAAIAAAADAAAAA0AAAAAAAAB2AAAAAIAAAABAAAABAAAAA8AAAAAAAAB7AAA" *
    "AAEAAAACAAAANgAAAAAAAAAA")

# Dimension 2, v3: same connectivity, 2 reals per vertex, no skip sections.
const _FX_2D = base64decode(
    "AQAAAAMAAAADAAAAGAAAAAAAAAACAAAABAAAAMgAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8AAAAAAAAAAAIAAAAAAAAAAADwPwAAAAAAAPA/AgAAAAAAAAAAAAAA" *
    "AAAAAAAA8D8CAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAADwPwAAAAAAAAAABAAAAAAA" *
    "AAAAAPA/AAAAAAAA8D8EAAAAAAAAAAAAAAAAAAAAAADwPwQAAAAFAAAA5AAAAAAAAAABAAAA" *
    "AQAAAAIAAAAHAAAABgAAAAQBAAAAAAAAAQAAAAEAAAACAAAAAwAAAAgAAAAHAAAAKAEAAAAA" *
    "AAABAAAAAQAAAAIAAAADAAAABAAAAAkAAAAIAAAATAEAAAAAAAABAAAAAQAAAAIAAAADAAAA" *
    "BQAAAAsAAAAKAAAAgAEAAAAAAAABAAAAAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgA" *
    "AAAMAAAANgAAAAAAAAAAAAAA")

# TetrahedraP2 (code 30) + EdgesP2 (code 25): file-order connectivity
# 1..10 must permute to (1..8,10,9) in the Gmsh-order block.
const _FX_TET10 = base64decode(
    "AQAAAAMAAAADAAAAGAAAAAAAAAADAAAABAAAAEABAAAAAAAACgAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAADwPwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAhAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQQAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABhA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAcQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" *
    "IEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACJAAAAAAAAAAAAAAAAAAAAAAAAAAAAeAAAA" *
    "fAEAAAAAAAABAAAAAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgAAAAJAAAACgAAACoA" *
    "AAAZAAAAnAEAAAAAAAABAAAAAQAAAAIAAAADAAAANwAAADYAAAAAAAAAAAAAAA==")

# Structural rejects: element keyword before `Vertices`; unrecognized
# keyword 77 whose next-keyword offset is 0; a `Vertices` offset that does
# not match the stream position; `Vertices` where `Dimension` was expected.
const _FX_BAD_ORDER = base64decode(
    "AQAAAAMAAAADAAAAGAAAAAAAAAADAAAABQAAADQAAAAAAAAAAQAAAAEAAAACAAAABwAAAAQA" *
    "AAAkAQAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA8D8AAAAA" *
    "AAAAAAAAAAAAAAAAAgAAAAAAAAAAAPA/AAAAAAAA8D8AAAAAAAAAAAIAAAAAAAAAAAAAAAAA" *
    "AAAAAPA/AAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAPA/BAAAAAAAAAAAAPA/" *
    "AAAAAAAAAAAAAAAAAADwPwQAAAAAAAAAAADwPwAAAAAAAPA/AAAAAAAA8D8EAAAAAAAAAAAA" *
    "AAAAAAAAAADwPwAAAAAAAPA/BAAAADYAAAAAAAAAAAAAAA==")
const _FX_BAD_KWD = base64decode(
    "AQAAAAMAAAADAAAAGAAAAAAAAAADAAAATQAAAAAAAAAAAAAAKgAAADYAAAAAAAAAAAAAAA==")
const _FX_BAD_POS = base64decode(
    "AQAAAAMAAAADAAAAGAAAAAAAAAADAAAABAAAAOcDAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAADwPwAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAA8D8AAAAA" *
    "AADwPwAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAAAIAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8EAAAAAAAAAAAA8D8AAAAAAAAAAAAAAAAAAPA/BAAAAAAAAAAAAPA/" *
    "AAAAAAAA8D8AAAAAAADwPwQAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAA8D8EAAAABgAAACgB" *
    "AAAAAAAAAQAAAAEAAAACAAAAAwAAAAgAAAA2AAAAAAAAAAAAAAA=")
const _FX_BAD_DIMKWD = base64decode(
    "AQAAAAMAAAAEAAAAGAAAAAAAAAADAAAABAAAAAgBAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAAAAAAAAADwPwAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAA8D8AAAAA" *
    "AADwPwAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAAAIAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8EAAAAAAAAAAAA8D8AAAAAAAAAAAAAAAAAAPA/BAAAAAAAAAAAAPA/" *
    "AAAAAAAA8D8AAAAAAADwPwQAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAA8D8EAAAABgAAACgB" *
    "AAAAAAAAAQAAAAEAAAACAAAAAwAAAAgAAAA2AAAAAAAAAAAAAAA=")

# meshio-produced (third-party, libMeshb-layout) file: v4, triangle/quad/
# tet/hex records with refs 5/9/11/3.
const _FX_MESHIO = base64decode(
    "AQAAAAQAAAADAAAAGAAAAAAAAAADAAAABAAAACwBAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAPA/AAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAA" *
    "AAAAAAAA8D8AAAAAAADwPwAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAA" *
    "AAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADwPwEAAAAAAAAAAAAAAAAA8D8AAAAA" *
    "AAAAAAAAAAAAAPA/AQAAAAAAAAAAAAAAAADwPwAAAAAAAPA/AAAAAAAA8D8BAAAAAAAAAAAA" *
    "AAAAAAAAAAAAAAAA8D8AAAAAAADwPwEAAAAAAAAABgAAAIABAAAAAAAAAgAAAAAAAAABAAAA" *
    "AAAAAAIAAAAAAAAAAwAAAAAAAAAFAAAAAAAAAAEAAAAAAAAAAwAAAAAAAAAEAAAAAAAAAAUA" *
    "AAAAAAAABwAAALwBAAAAAAAAAQAAAAAAAAAFAAAAAAAAAAYAAAAAAAAABwAAAAAAAAAIAAAA" *
    "AAAAAAkAAAAAAAAACAAAAPgBAAAAAAAAAQAAAAAAAAABAAAAAAAAAAIAAAAAAAAAAwAAAAAA" *
    "AAAFAAAAAAAAAAsAAAAAAAAACgAAAFQCAAAAAAAAAQAAAAAAAAABAAAAAAAAAAIAAAAAAAAA" *
    "AwAAAAAAAAAEAAAAAAAAAAUAAAAAAAAABgAAAAAAAAAHAAAAAAAAAAgAAAAAAAAAAwAAAAAA" *
    "AAA2AAAAAAAAAAAAAAA=")

const _CUBE_COORDS = [0.0 1.0 1.0 0.0 0.0 1.0 1.0 0.0;
                      0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
                      0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0]
const _CUBE_BLOCKS = [1, 2, 3, 4, 5]
const _CUBE_NODES = [Int32[1; 2;;], Int32[1; 2; 3;;], Int32[1; 2; 3; 4;;],
                     Int32[1; 2; 3; 5;;], reshape(Int32(1):Int32(8), 8, 1)]
const _CUBE_REFS = [Int32[7], Int32[8], Int32[9], Int32[11], Int32[12]]

function _fx(dir,name,bytes)
    path=joinpath(dir,name)
    Base.write(path,bytes)
    return path
end

function _check_cube(f)
    @test f isa MeshFile
    @test size(f.mesh.coords,2)==8
    @test f.mesh.coords==_CUBE_COORDS
    @test [b.msh for b in f.mesh.blocks]==_CUBE_BLOCKS
    @test [b.nodes for b in f.mesh.blocks]==_CUBE_NODES
    @test [b.tags for b in f.mesh.blocks]==_CUBE_REFS
    @test f.mesh.elementary_entities==_CUBE_REFS
    @test f.node_tags==UInt64.(1:8)
end

@testset "binary Medit (.meshb) I/O" begin

    @testset "reader: versions, widths, endianness, skips" begin
        mktempdir() do dir
            # v3 little-endian: i64 offsets, f64 coords, i32 records; the
            # Corners/RequiredVertices keywords are hopped by offset.
            _check_cube(read_meshb(_fx(dir,"v3.meshb",_FX_V3)))
            # v1: f32 coords and i32 offsets/records.
            _check_cube(read_meshb(_fx(dir,"v1.meshb",_FX_V1)))
            # v2: f64 coords with i32 offsets/records.
            _check_cube(read_meshb(_fx(dir,"v2.meshb",_FX_V2)))
            # v4: i64 offsets and i64 records.
            _check_cube(read_meshb(_fx(dir,"v4.meshb",_FX_V4)))
            # big-endian v3 detected through the byte-swapped magic tag.
            _check_cube(read_meshb(_fx(dir,"v3be.meshb",_FX_V3BE)))
            # Dimension 2 stores two reals per vertex; z fills with 0.
            f2=read_meshb(_fx(dir,"d2.meshb",_FX_2D))
            @test size(f2.mesh.coords,2)==8
            @test all(==(0.0),f2.mesh.coords[3,:])
            @test [b.msh for b in f2.mesh.blocks]==_CUBE_BLOCKS
            # A real meshio-produced file decodes with the same contract.
            fm=read_meshb(_fx(dir,"meshio.meshb",_FX_MESHIO))
            @test [b.msh for b in fm.mesh.blocks]==[2,3,4,5]
            @test [b.tags for b in fm.mesh.blocks]==
                  [Int32[5,5],Int32[9],Int32[11],Int32[3]]
            # TetrahedraP2/EdgesP2 keywords permute like the ASCII reader.
            ft=read_meshb(_fx(dir,"t10.meshb",_FX_TET10))
            @test [b.msh for b in ft.mesh.blocks]==[11,8]
            @test ft.mesh.blocks[1].nodes==reshape(Int32[1,2,3,4,5,6,7,8,10,9],10,1)
            @test ft.mesh.blocks[1].tags==Int32[42]
            @test ft.mesh.blocks[2].nodes==reshape(Int32[1,2,3],3,1)
            @test ft.mesh.blocks[2].tags==Int32[55]
        end
    end

    @testset "reader: failure paths" begin
        mktempdir() do dir
            base=_fx(dir,"ok.meshb",_FX_V3)
            function _mut(name,f)
                bytes=copy(_FX_V3);f(bytes)
                return _fx(dir,name,bytes)
            end
            # bad MeshVersionFormatted code (and an ASCII .mesh byte stream)
            @test_throws ArgumentError read_meshb(_mut("magic.meshb",
                b->(b[1]=2)))
            ascii=_fx(dir,"ascii.meshb",codeunits(
                "MeshVersionFormatted 2\n Dimension 3\n Vertices 0\n End\n"))
            @test_throws ArgumentError read_meshb(ascii)
            # unsupported version word
            @test_throws ArgumentError read_meshb(_mut("ver.meshb",
                b->(b[5:8].=reinterpret(UInt8,UInt32[9]))))
            # Dimension keyword missing / bad dimension value
            @test_throws ArgumentError read_meshb(_fx(
                dir,"dimkwd.meshb",_FX_BAD_DIMKWD))
            @test_throws ArgumentError read_meshb(_mut("dim.meshb",
                b->(b[21:24].=reinterpret(UInt8,UInt32[4]))))
            # elements before Vertices, dangling vertex, negative ref
            @test_throws ArgumentError read_meshb(_fx(
                dir,"order.meshb",_FX_BAD_ORDER))
            @test_throws ArgumentError read_meshb(_mut("dangle.meshb",
                b->(b[309:312].=reinterpret(UInt8,UInt32[99]))))
            @test_throws ArgumentError read_meshb(_mut("negref.meshb",
                b->(b[321:324].=reinterpret(UInt8,UInt32[0xFFFFFFFF]))))
            # inconsistent / unskippable offsets
            @test_throws ArgumentError read_meshb(_fx(
                dir,"badpos.meshb",_FX_BAD_POS))
            @test_throws ArgumentError read_meshb(_fx(
                dir,"badkwd.meshb",_FX_BAD_KWD))
            # truncated inside a Triangles record
            @test_throws ArgumentError read_meshb(_fx(
                dir,"trunc.meshb",_FX_V3[1:320]))
            # empty file and missing file
            @test_throws ArgumentError read_meshb(_fx(dir,"empty.meshb",UInt8[]))
            @test_throws ArgumentError read_meshb(joinpath(dir,"missing.meshb"))
            # limits
            @test_throws ArgumentError read_meshb(base;max_nodes=7)
            @test_throws ArgumentError read_meshb(base;max_elements=4)
            @test_throws ArgumentError read_meshb(base;max_file_bytes=10)
            @test_throws ArgumentError read_meshb(42)
            # non-finite vertex coordinate
            @test_throws ArgumentError read_meshb(_mut("nan.meshb",
                b->(b[41:48].=reinterpret(UInt8,UInt64[0x7FF8000000000000]))))
            # the text reader must not silently accept binary bytes
            @test_throws ArgumentError read_mesh(base)
        end
    end

    @testset "writer: all versions round-trip" begin
        mktempdir() do dir
            source=read_meshb(_fx(dir,"src.meshb",_FX_V3)).mesh
            for version in 1:4
                path=joinpath(dir,"rt$version.meshb")
                @test write_meshb(path,source;version=version)==path
                f=read_meshb(path)
                @test f.mesh.coords==source.coords          # 0/1: exact in f32
                @test [b.msh for b in f.mesh.blocks]==_CUBE_BLOCKS
                @test [b.nodes for b in f.mesh.blocks]==_CUBE_NODES
                @test [b.tags for b in f.mesh.blocks]==_CUBE_REFS
                @test f.mesh.elementary_entities==_CUBE_REFS
            end
            # default version picks 3 (Int32 counts): i64 offsets, i32 records
            path=joinpath(dir,"default.meshb")
            write_meshb(path,source)
            bytes=read(path)
            @test reinterpret(Int32,bytes[1:4])[]==Int32(1)
            @test reinterpret(Int32,bytes[5:8])[]==Int32(3)
            _check_cube(read_meshb(path))
            # version=1 stores f32 reals; exact fractions survive the trip.
            fp=joinpath(dir,"f32.meshb")
            write_meshb(fp,MixedMesh([0.0 1 0.5; 0 0 0.25; 0 0 0],
                [ElementBlock(2,Int32[1;2;3;;],Int32[9])]);version=1)
            g1=read_meshb(fp)
            @test g1.mesh.coords==[0.0 1 0.5; 0 0 0.25; 0 0 0]
            @test only(g1.mesh.blocks).tags==Int32[9]
            # a ref at the Int32 ceiling still fits the v3 field widths
            big=MixedMesh(zeros(3,4),[ElementBlock(4,Int32[1;2;3;4;;],
                Int32[typemax(Int32)])])
            bp2=joinpath(dir,"ceil.meshb")
            write_meshb(bp2,big;version=3)
            @test only(read_meshb(bp2).mesh.blocks).tags==Int32[typemax(Int32)]
            # version must stay in the 1–4 libMeshb range
            @test_throws ArgumentError write_meshb(
                joinpath(dir,"w9.meshb"),source;version=9)
            @test_throws ArgumentError write_meshb(
                joinpath(dir,"w5.meshb"),source;version=5)
            @test_throws ArgumentError write_meshb(
                joinpath(dir,"w0.meshb"),source;version=-1)
            # unrepresentable families reject before touching the target
            quad9=MixedMesh(zeros(3,9),
                [ElementBlock(10,reshape(Int32[1:9;],9,1),Int32[1])])
            @test_throws ArgumentError write_meshb(
                joinpath(dir,"q.meshb"),quad9)
            @test !isfile(joinpath(dir,"q.meshb"))
            special=MixedMesh(zeros(3,2),
                [SpecialElementBlock(67,reshape(Int32[1,2],2,1))])
            @test_throws ArgumentError write_meshb(
                joinpath(dir,"s.meshb"),special)
            @test_throws ArgumentError write_meshb(
                joinpath(dir,"nodir","x.meshb"),source)
            # empty meshes still emit the Vertices-0 skeleton
            empty=MixedMesh(zeros(3,0),ElementBlock[])
            ep=joinpath(dir,"empty.meshb")
            @test write_meshb(ep,empty)==ep
            @test size(read_meshb(ep).mesh.coords,2)==0
        end
    end

    @testset "binary/ASCII parity and P2 permutation on write" begin
        mktempdir() do dir
            # The P2-tet block exercises the same swap the ASCII writer uses.
            mixed=MixedMesh(
                reduce(hcat,([Float64(i),0.0,0.0] for i in 1:10)),
                ElementBlock[ElementBlock(11,
                    reshape(Int32(1):Int32(10),10,1),Int32[42]),
                    ElementBlock(8,Int32[1;2;3;;],Int32[55])])
            mp=joinpath(dir,"p2.mesh");bp=joinpath(dir,"p2.meshb")
            write_mesh(mp,mixed);write_meshb(bp,mixed)
            ta=read_mesh(mp);tb=read_meshb(bp)
            @test ta.mesh.coords==tb.mesh.coords
            # EdgesP2 precedes TetrahedraP2 in the Medit section order
            @test [b.msh for b in ta.mesh.blocks]==[8,11]
            @test [b.msh for b in tb.mesh.blocks]==[8,11]
            @test [b.nodes for b in ta.mesh.blocks]==
                  [b.nodes for b in tb.mesh.blocks]
            @test [b.tags for b in ta.mesh.blocks]==
                  [b.tags for b in tb.mesh.blocks]
            # A binary write→read of the cube matches the ASCII channel too.
            source=read_meshb(_fx(dir,"src.meshb",_FX_V3)).mesh
            write_mesh(joinpath(dir,"x.mesh"),source)
            write_meshb(joinpath(dir,"x.meshb"),source)
            a=read_mesh(joinpath(dir,"x.mesh"))
            b=read_meshb(joinpath(dir,"x.meshb"))
            @test a.mesh.coords==b.mesh.coords
            @test [bl.nodes for bl in a.mesh.blocks]==
                  [bl.nodes for bl in b.mesh.blocks]
            @test [bl.tags for bl in a.mesh.blocks]==
                  [bl.tags for bl in b.mesh.blocks]
        end
    end

    @testset "session open/write/merge route" begin
        _MESHB_API.initialize()
        try
            _MESHB_API.option("General.Terminal",0)
            # like upstream, non-physical entities export only with SaveAll
            _MESHB_API.option("Mesh.SaveAll",1.0)
            _MESHB_API.model.add("host")
            _MESHB_API.model.add_discrete_entity(1,7)
            _MESHB_API.model.add_discrete_entity(2,8)
            _MESHB_API.model.add_discrete_entity(3,9)
            _MESHB_API.mesh.add_nodes(1,7,[1,2],[0.0,0,0,1,0,0])
            _MESHB_API.mesh.add_nodes(2,8,[3,4,5],
                [0.0,1,0, 1,1,0, 0,1,1])
            _MESHB_API.mesh.add_elements_by_type(7,1,[11],[1,2])
            _MESHB_API.mesh.add_elements_by_type(8,2,[12],[3,4,5])
            _MESHB_API.mesh.add_elements_by_type(9,4,[13],[1,2,3,5])
            _MESHB_API.model.add_physical_group(2,[8];tag=60,name="surf")
            host=_MESHB_API.model.get_entities()
            @test (1,7) in host && (2,8) in host && (3,9) in host

            mktempdir() do dir
                # API.write emits the binary GMF stream; API.open re-imports it.
                out=joinpath(dir,"host.meshb")
                @test _MESHB_API.write(out)===nothing
                bytes=read(out)
                @test reinterpret(Int32,bytes[1:4])[]==Int32(1)
                _MESHB_API.model.add("re")
                @test _MESHB_API.open(out)===nothing
                @test "host" in _MESHB_API.model.list()
                entities=_MESHB_API.model.get_entities()
                # ref = physical (60) for the surfaced tri, entity tags else
                @test (1,7) in entities && (2,60) in entities &&
                      (3,9) in entities
                @test (2,60) in _MESHB_API.model.get_physical_groups()
                types=_MESHB_API.mesh.get_elements()[1]
                @test Int32(1) in types && Int32(2) in types &&
                      Int32(4) in types
                # a file fixture merges into the live model like `.mesh`.
                fx=_fx(dir,"fx.meshb",_FX_V3)
                before=_MESHB_API.model.get_entities()
                @test _MESHB_API.merge(fx)===nothing
                after=_MESHB_API.model.get_entities()
                @test issubset(before,after)
                @test (1,7) in after && (2,8) in after && (2,9) in after &&
                      (3,11) in after && (3,12) in after
                merged=_MESHB_API.mesh.get_elements()[1]
                @test Int32(3) in merged && Int32(5) in merged
                # dispatch guards stay consistent with the other formats.
                @test_throws ArgumentError _MESHB_API.open(
                    joinpath(dir,"missing.meshb"))
                @test_throws ArgumentError _MESHB_API.merge(
                    joinpath(dir,"missing.meshb"))
                # ASCII content under a .meshb name fails the binary reader.
                bad=_fx(dir,"bad.meshb",codeunits(
                    "MeshVersionFormatted 2\n Dimension 3\n End\n"))
                @test_throws ArgumentError _MESHB_API.open(bad)
            end
        finally
            _MESHB_API.finalize()
        end
    end
end
