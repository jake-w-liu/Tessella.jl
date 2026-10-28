# ── Long-tail surface/mesh decks: .off/.obj/.ply2 (write+read round trip) and
#    .dat/.m/.wrl/.key/.rad/.x3d (write) plus .celum rejection ────────────────
#
# Correctness: a discrete-surface fixture (unit square: 2 triangles + 1 quad
#   face on surface 1, a boundary line on curve 9, a point cell) is exported
#   through `API.write` for every new extension; OFF/OBJ/PLY2 re-import through
#   `API.open` and are compared on coordinates and face connectivity counts.
#   Write-only decks are checked structurally against the upstream layouts
#   captured from Gmsh 4.15.2 GModelIO_{TOCHNOG,MATLAB,VRML,KEY,RAD,X3D}.cpp.
# Robustness : .celum must reject truthfully (no parametric-surface normals/
#   curvatures on the foreign projection).

using Test
import Tessella

const _SD_API = Tessella.API

function _sd_counts()
    els = _SD_API.mesh.get_elements()
    ntri = 0; nquad = 0; nline = 0; npnt = 0
    types = els[1]; tags = els[2]
    for i in eachindex(types)
        t = Int(types[i]); n = length(tags[i])
        t == 2 && (ntri += n); t == 3 && (nquad += n)
        t == 1 && (nline += n); t == 15 && (npnt += n)
    end
    return (ntri=ntri, nquad=nquad, n2=ntri+nquad, nline=nline, npnt=npnt,
        nnodes=length(_SD_API.mesh.get_nodes()[1]))
end

@testset "long-tail mesh decks" begin
    _SD_API.initialize()
    try
        _SD_API.option("General.Terminal", 0)
        _SD_API.option("Mesh.SaveAll", 1.0)
        _SD_API.model.add("deck")
        _SD_API.model.add_discrete_entity(2, 1)
        _SD_API.model.add_discrete_entity(1, 9)
        _SD_API.model.add_discrete_entity(0, 5)
        # Square quad (1-4) plus a detached triangle (5-7); a 2-node line on
        # curve 9 and a point cell — no collinear faces.
        _SD_API.mesh.add_nodes(2, 1, [1, 2, 3, 4, 5, 6, 7],
            [0.0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
             2.0, 0, 0, 3, 0, 0, 2.5, 1, 0])
        _SD_API.mesh.add_nodes(1, 9, [8, 9], [0.0, 0.5, 0, 1.0, 0.5, 0])
        _SD_API.mesh.add_nodes(0, 5, [10], [0.0, 0.0, 0])
        _SD_API.mesh.add_elements_by_type(1, 3, [12], [1, 2, 3, 4])
        _SD_API.mesh.add_elements_by_type(1, 2, [10], [5, 6, 7])
        _SD_API.mesh.add_elements_by_type(9, 1, [13], [8, 9])
        _SD_API.mesh.add_elements_by_type(5, 15, [14], [10])
        _SD_API.model.add_physical_group(2, [1]; tag=7, name="plate")
        counts = _sd_counts()
        @test counts.ntri == 1 && counts.nquad == 1 && counts.nline == 1 &&
            counts.nnodes == 10
        mktempdir() do dir
            # Round-trippable decks.
            for ext in ("off", "obj", "ply2")
                p = joinpath(dir, "mesh." * ext)
                @test _SD_API.write(p) === nothing
                @test isfile(p)
            end
            off = readlines(joinpath(dir, "mesh.off"))
            hdr = split(off[1])
            @test hdr[1] == "OFF"
            @test parse(Int, hdr[2]) == counts.nnodes
            @test parse(Int, hdr[3]) == counts.n2
            # OFF re-import: nodes and 3/4 face records survive (0-based).
            _SD_API.model.add("reoff")
            @test _SD_API.open(joinpath(dir, "mesh.off")) === nothing
            re = _sd_counts()
            # upstream `_storeVerticesInEntities` drops vertices no face
            # references on read-back: the 7 surface nodes survive, the
            # line/point nodes are dropped.
            @test re.nnodes == 7
            @test re.ntri == counts.ntri && re.nquad == counts.nquad
            obj = readlines(joinpath(dir, "mesh.obj"))
            @test obj[1] == "# Gmsh OBJ export"
            @test count(l->startswith(l, "v "), obj) == counts.nnodes
            @test count(l->startswith(l, "f "), obj) == counts.n2
            _SD_API.model.add("reobj")
            @test _SD_API.open(joinpath(dir, "mesh.obj")) === nothing
            @test _sd_counts().n2 == counts.n2
            ply = readlines(joinpath(dir, "mesh.ply2"))
            @test parse(Int, ply[1]) == counts.nnodes
            @test parse(Int, ply[2]) == counts.ntri
            _SD_API.model.add("reply")
            @test _SD_API.open(joinpath(dir, "mesh.ply2")) === nothing
            @test _sd_counts().ntri == counts.ntri

            # Write-only decks — structural checks against upstream layouts.
            # Re-select the original model: `API.open` above switched the
            # current model to the last re-import.
            _SD_API.model.set_current("deck")
            for ext in ("dat", "m", "wrl", "key", "rad", "x3d")
                p = joinpath(dir, "mesh." * ext)
                @test _SD_API.write(p) === nothing
                @test isfile(p)
                @test filesize(p) > 0
            end
            dat = read(joinpath(dir, "mesh.dat"), String)
            @test occursin("Tochnog Input File Created by Gmsh", dat)
            @test occursin(r"node \d+ [\d\.eE+-]+ [\d\.eE+-]+", dat)
            @test occursin(r"element \d+ -tria3", dat)
            @test occursin("element_group", dat)
            mm = read(joinpath(dir, "mesh.m"), String)
            @test occursin("%% Matlab mesh", mm)
            @test occursin("msh.nbNod = $(counts.nnodes);", mm)
            @test occursin("msh.TRIANGLES =[", mm)
            @test occursin("msh.QUADS =[", mm)
            wrl = read(joinpath(dir, "mesh.wrl"), String)
            @test occursin("#VRML V1.0 ascii", wrl)
            @test occursin("Coordinate3 {", wrl)
            @test occursin("IndexedFaceSet {", wrl)
            @test occursin("IndexedLineSet {", wrl)
            key = read(joinpath(dir, "mesh.key"), String)
            @test occursin("*KEYWORD", key)
            @test occursin("*NODE", key)
            @test occursin("*ELEMENT_SHELL", key)
            @test occursin("*END", key)
            rad = read(joinpath(dir, "mesh.rad"), String)
            @test occursin("#RADIOSS STARTER", rad)
            @test occursin("/BEGIN", rad)
            @test occursin("/NODE", rad)
            @test occursin("/SH3N/", rad)
            x3d = read(joinpath(dir, "mesh.x3d"), String)
            @test occursin("<X3D profile='Interchange'", x3d)
            @test occursin("<TriangleSet", x3d)
            @test isfile(joinpath(dir, "mesh_x3d.html"))
            # CELUM rejects truthfully — no parametric normals upstream of the
            # mesh projection.
            @test_throws ArgumentError _SD_API.write(joinpath(dir, "mesh.celum"))
        end
    finally
        _SD_API.finalize()
    end
end

# ── Stanford .ply reader: ASCII + binary little/big-endian ───────────────────
# Upstream `GModel::readPLY` accepts all three encodings (binary marked
# experimental); `.ply` is read-only there — the writable dialect is `.ply2`.
# Faces are triangles built from the first three indices regardless of the
# declared count (upstream `sscanf` semantics); unreferenced vertices are
# dropped exactly like `_storeVerticesInEntities`.
function _sd_write_ply_le(path, coords, faces; be=false, extra=false)
    open(path, "w") do io
        endian = be ? "binary_big_endian" : "binary_little_endian"
        props = extra ? "property uchar red\n" : ""
        print(io, "ply\nformat $endian 1.0\ncomment devin ply fixture\n",
            "element vertex $(length(coords))\n",
            "property float x\nproperty float y\nproperty float z\n",
            props,
            "element face $(length(faces))\n",
            "property list uchar int vertex_indices\nend_header\n")
        pack = be ? (v -> bswap(v)) : identity
        for (i, xyz) in enumerate(coords)
            for x in xyz; write(io, pack(Float32(x))); end
            extra && write(io, UInt8(i))
        end
        for f in faces
            write(io, UInt8(length(f)))
            for i in f; write(io, pack(Int32(i - 1))); end
        end
    end
    return path
end

@testset "ply reader (ascii/binary)" begin
    _SD_API.initialize()
    try
        _SD_API.option("General.Terminal", 0)
        mktempdir() do dir
            # ASCII: square quad face + detached triangle + one unreferenced
            # vertex. Upstream takes the first three quad indices → 2 tris.
            asc = joinpath(dir, "tri.ply")
            write(asc, """
ply
format ascii 1.0
comment ascii fixture
element vertex 6
property float x
property float y
property float z
element face 2
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
1 1 0
0 1 0
3 0 0
4 0 0
4 0 1 2 3
3 4 5 3
""")
            _SD_API.model.add("ply_a")
            @test _SD_API.open(asc) === nothing
            c = _sd_counts()
            @test c.ntri == 2 && c.nnodes == 6  # all 6 verts referenced

            le = joinpath(dir, "tri_le.ply")
            _sd_write_ply_le(le, [(0.0, 0, 0), (1, 0, 0), (1, 1, 0),
                (0, 1, 0), (9, 9, 9)], [[1, 2, 3], [1, 3, 4]])
            _SD_API.model.add("ply_le")
            @test _SD_API.open(le) === nothing
            c = _sd_counts()
            @test c.ntri == 2 && c.nnodes == 4  # vertex 5 unreferenced → dropped
            nodes = _SD_API.mesh.get_nodes()
            xs = nodes[2]
            @test maximum(xs[1:3:end]) ≈ 1.0  # coords decoded little-endian

            be = joinpath(dir, "tri_be.ply")
            _sd_write_ply_le(be, [(0.0, 0, 0), (1, 0, 0), (1, 1, 0),
                (0, 1, 0)], [[1, 2, 3], [1, 3, 4]]; be=true)
            _SD_API.model.add("ply_be")
            @test _SD_API.open(be) === nothing
            c = _sd_counts()
            @test c.ntri == 2 && c.nnodes == 4
            xs = _SD_API.mesh.get_nodes()[2]
            @test maximum(xs[1:3:end]) ≈ 1.0  # big-endian swap honored

            # Extra scalar vertex property: readable but dropped (upstream
            # would build a per-node PView; the mesh bridge has no channel).
            xp = joinpath(dir, "tri_prop.ply")
            _sd_write_ply_le(xp, [(0.0, 0, 0), (1, 0, 0), (1, 1, 0)],
                [[1, 2, 3]]; extra=true)
            _SD_API.model.add("ply_prop")
            @test _SD_API.open(xp) === nothing
            @test _sd_counts().ntri == 1

            # Rejections.
            write(joinpath(dir, "bad_fmt.ply"), "ply\nformat garbage 1.0\nelement vertex 0\nproperty float x\nproperty float y\nproperty float z\nend_header\n")
            @test_throws ArgumentError begin
                _SD_API.model.add("r1"); _SD_API.open(joinpath(dir, "bad_fmt.ply"))
            end
            write(joinpath(dir, "bad_elem.ply"), "ply\nformat ascii 1.0\nelement edge 1\nproperty int x\nend_header\n")
            @test_throws ArgumentError begin
                _SD_API.model.add("r2"); _SD_API.open(joinpath(dir, "bad_elem.ply"))
            end
            # binary quad face: upstream reads exactly 3 indices per record,
            # desyncing on 4 — Tessella rejects instead of fabricating.
            q = joinpath(dir, "quad_be.ply")
            _sd_write_ply_le(q, [(0.0, 0, 0), (1, 0, 0), (1, 1, 0),
                (0, 1, 0)], [[1, 2, 3, 4]])
            @test_throws ArgumentError begin
                _SD_API.model.add("r3"); _SD_API.open(q)
            end
            # `.ply` is read-only upstream: writes still reject.
            @test_throws ArgumentError _SD_API.write(joinpath(dir, "out.ply"))
        end
    finally
        _SD_API.finalize()
    end
end

# ── Dependency-named rejections (.med/.cgns) and .cdf fallthrough ────────────
# Upstream `.med`/`.mmed`/`.omm` and `.cgns` reads require optional libraries
# this project does not link; the diagnostic must name the missing dependency.
# Upstream has no `.cdf` dispatch at all — unknown extensions fall through to
# the `.geo` parser, which Tessella reproduces literally.
@testset "format rejection diagnostics" begin
    _SD_API.initialize()
    try
        _SD_API.option("General.Terminal", 0)
        _SD_API.model.add("rej")
        mktempdir() do dir
            for ext in ("med", "mmed", "omm", "cgns")
                p = joinpath(dir, "x." * ext)
                write(p, "garbage\n")
                err = try
                    _SD_API.open(p); nothing
                catch e; e; end
                @test err isa ArgumentError
                msg = sprint(showerror, err)
                @test occursin(
                    ext == "cgns" ? "CGNS" : "MED", msg)
                @test occursin("compiled with", msg)
                err = try
                    _SD_API.merge(p); nothing
                catch e; e; end
                @test err isa ArgumentError
            end
            # `.cdf` is routed to the .geo executor like upstream's readGEO
            # fallthrough: a valid script opens, garbage fails as a parse
            # error rather than an extension complaint.
            good = joinpath(dir, "pt.cdf")
            write(good, "Point(1) = {0,0,0,1};\n")
            _SD_API.model.add("cdfsrc")
            @test _SD_API.open(good) === nothing
            bad = joinpath(dir, "raw.cdf")
            write(bad, "\x00\x01\x02binary blob\n")
            err = try
                _SD_API.model.add("cdfbad"); _SD_API.open(bad); nothing
            catch e; e; end
            @test err isa ArgumentError
            @test !occursin("supported file extensions", sprint(showerror, err))
        end
    finally
        _SD_API.finalize()
    end
end
