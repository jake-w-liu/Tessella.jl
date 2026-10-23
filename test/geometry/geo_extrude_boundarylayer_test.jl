using Test
using Tessella
import Tessella.MeshTypes as MT

# `.geo` boundary-layer extrusion — `Extrude {Curve{...}; Layers{...}}`:
# the no-displacement form grows a lateral strip along the smoothed wall
# normal field (upstream's `BOUNDARY_LAYER` Extrude mode) instead of a rigid
# sweep. The generated entities keep the ordinary extrusion record set —
# `DuplicateCurve` chapeau, `ExtrudePoint` connectors, `Layers` params — so
# `Mesh 2` sweeps them through `mesh_boundary_layer_2d` coherently.

function _execute_bl_source(source::AbstractString; mesh_dim=0)
    return mktemp() do path, io
        write(io, source)
        close(io)
        execute_geo(path; mesh_dim=mesh_dim)
    end
end

function _bl_error(source::AbstractString)
    try
        ex=_execute_bl_source(source)
        # Recoverable `Error   :` diagnostics surface through
        # `exec.errors`, not a throw — re-materialize them so the
        # assertions below keep checking the upstream text.
        !isempty(ex.errors) ?
            ArgumentError(join(ex.errors,"\n")) : nothing
    catch err
        return err
    end
end

function _bl_cells(mesh)
    tris = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 2; init = 0)
    quads = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 3; init = 0)
    return (tris=tris, quads=quads)
end

# MixedMesh surface area via the (1,k,k+1) triangle fan — covers tri and
# quad blocks alike.
function _bl_area(mesh)
    area = 0.0
    for b in mesh.blocks, col in eachcol(b.nodes)
        pts = [Tuple(mesh.coords[:, Int(v)]) for v in col]
        for k in 2:length(pts)-1
            area += MT.triangle_area(pts[1], pts[k], pts[k+1])
        end
    end
    return area
end

@testset ".geo boundary-layer curve Extrude" begin
    # `Extrude {Curve{1}; Layers{3}}` on the unit x-line grows the strip on
    # the `t × ẑ` side: tops at y=-1 with the ordinary curve-extrude tag
    # layout out=[chapeau, surface, connector_end, -connector_beg].
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{3}; Recombine; };
        """)
    @test r.lists["out"] == [2.0, 5.0, 4.0, -3.0]
    @test r.model.curves[2] == (3, 4)       # chapeau copy on the top row
    @test r.model.curves[3] == (1, 3)       # connector at the signed begin
    @test r.model.curves[4] == (2, 4)       # connector at the signed end
    @test r.model.loops[only(r.model.surfaces[5])] == [1, 4, -2, -3]
    @test r.model.points[3] == (0.0, -1.0, 0.0)
    @test r.model.points[4] == (1.0, -1.0, 0.0)
    # The generated entities keep their `Layers` sweep records.
    @test r.model.meshing.extrude[(2, 5)].layers == [3]
    @test r.model.meshing.extrude_sources[(2, 5)] == (1, 1)
    @test r.model.meshing.extrude_sources[(1, 2)] == (1, 1)
    @test r.model.meshing.extrude_sources[(1, 3)] == (0, 1)
    @test r.model.meshing.extrude_sources[(1, 4)] == (0, 2)

    # `Mesh 2`: the recombined strip emits one quadrangle per layer row; the
    # connectors carry the same level columns, and the chapeau meshes like
    # the generatrix copy — all welded onto the strip's boundary nodes.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{3}; Recombine; };
        Mesh 2;
        """)
    surface = geo_entity_mesh(r, 2, 5)
    @test _bl_cells(surface) == (tris=0, quads=3)
    ys = sort!(unique!(round.(vec(surface.coords[2, :]); digits=12)))
    @test ys ≈ [-1.0, -2/3, -1/3, 0.0]
    @test _bl_cells(geo_entity_mesh(r, 1, 2)) == (tris=0, quads=0)
    conn = geo_entity_mesh(r, 1, 3)
    @test sum(size(b.nodes, 2) for b in conn.blocks) == 3
    @test sort!(vec(conn.coords[2, :])) ≈ [-1.0, -2/3, -1/3, 0.0]

    # Without Recombine the kernel's quadrangles split on the lateral
    # diagonal, matching the ordinary non-recombined strip layout.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{3}; };
        Mesh 2;
        """)
    @test _bl_cells(geo_entity_mesh(r, 2, 5)) == (tris=6, quads=0)

    # `Layers{3,5}` — counts and an absolute cumulative height of 5.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{3,5}; Recombine; };
        Mesh 2;
        """)
    surface = geo_entity_mesh(r, 2, 5)
    @test _bl_cells(surface) == (tris=0, quads=3)
    @test r.model.points[3] == (0.0, -5.0, 0.0)

    # Documented nested `Layers{{c1,c2},{h1,h2}}` — cumulative heights.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{{2,4},{0.2,0.8}}; Recombine; };
        Mesh 2;
        """)
    surface = geo_entity_mesh(r, 2, 5)
    @test _bl_cells(surface) == (tris=0, quads=6)
    @test r.model.points[3] == (0.0, -0.8, 0.0)

    # A negative generatrix extrudes the reversed record and keeps its
    # signed tag in the lateral tail, like the ordinary path.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{-1}; Layers{3}; Recombine; };
        """)
    @test r.lists["out"] == [2.0, 5.0, -1.0, 4.0, -3.0]
    @test r.model.loops[only(r.model.surfaces[5])] == [-1, 4, -2, -3]
    @test r.model.points[3] == (1.0, 1.0, 0.0)
    @test r.model.points[4] == (0.0, 1.0, 0.0)

    # Sibling walls in one statement share the smoothed normal field: the
    # two coincident connectors off the shared endpoint merge through the
    # operation's coherence pass, and each strip keeps welding columns.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {2,0,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        out[] = Extrude { Curve{1}; Curve{2}; Layers{3}; Recombine; };
        Mesh 2;
        """)
    @test length(r.lists["out"]) == 8
    s1 = geo_entity_mesh(r, 2, Int(r.lists["out"][2]))
    s2 = geo_entity_mesh(r, 2, Int(r.lists["out"][6]))
    @test _bl_cells(s1) == (tris=0, quads=3)
    @test _bl_cells(s2) == (tris=0, quads=3)
    @test _bl_area(s1) + _bl_area(s2) ≈ 2.0

    # An L-shaped wall pair: the corner vertex's merged growth direction is
    # the normalized sum of both segments' normals (upstream's smooth_data
    # position accumulation), so both strips weld onto one corner column.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {1,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        out[] = Extrude { Curve{1}; Curve{2}; Layers{1}; Recombine; };
        """)
    corner_tops = Set(p for p in values(r.model.points) if p != (0.0,0.0,0.0) &&
        p != (1.0,0.0,0.0) && p != (1.0,1.0,0.0))
    @test corner_tops == Set([(0.0,-1.0,0.0),(1.0+sqrt(0.5),-sqrt(0.5),0.0),
                              (2.0,1.0,0.0)])

    # Boundary-layer extrusion still requires `Layers` and rejects the
    # parameters upstream does not honor on this mode.
    for src in ("out[] = Extrude { Curve{1}; Recombine; };",
                "out[] = Extrude { Curve{1}; Layers{3}; ScaleLastLayer; };",
                "out[] = Extrude { Curve{1}; Layers{3}; QuadTriAddVerts; };",
                "out[] = Extrude { Surface{1}; Layers{3}; };")
        @test _bl_error("""
            Point(1) = {0,0,0,1};
            Point(2) = {1,0,0,1};
            Line(1) = {1,2};
            Plane Surface(1) = {};
            $src
            """) isa ArgumentError
    end

    # Under the OpenCASCADE factory upstream rejects boundary-layer
    # extrusion in the parser action.
    @test _bl_error("""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{3}; };
        """) isa ArgumentError

    # Ordinary extrusions are untouched: a displacement still sweeps.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude {0,0,1} { Curve{1}; Layers{3}; Recombine; };
        Mesh 2;
        """)
    @test r.lists["out"] == [2.0, 5.0, 4.0, -3.0]
    surface = geo_entity_mesh(r, 2, 5)
    @test _bl_cells(surface) == (tris=0, quads=3)
    zs = extrema(vec(surface.coords[3, :]))
    @test zs == (0.0, 1.0)
end
