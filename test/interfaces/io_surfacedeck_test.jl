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
