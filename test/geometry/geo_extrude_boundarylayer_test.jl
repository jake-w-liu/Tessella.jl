using Test
using Tessella
import Tessella.MeshTypes as MT

# `.geo` boundary-layer extrusion — `Extrude {Curve{...}; Layers{...}}`:
# the no-displacement form grows a lateral strip along the smoothed wall
# normal field (upstream's `BOUNDARY_LAYER` Extrude mode) instead of a rigid
# sweep. The generated entities keep the ordinary extrusion record set —
# `DuplicateCurve` chapeau, `ExtrudePoint` connectors, `Layers` params — so
# `Mesh 2` sweeps them through `mesh_boundary_layer_2d` coherently.

function _execute_bl_source(source::AbstractString; mesh_dim=0,
                            files=Dict{String,String}())
    dir=mktempdir()
    for (name,content) in files
        write(joinpath(dir,name),content)
    end
    path=joinpath(dir,"script.geo")
    write(path,source)
    return execute_geo(path; mesh_dim=mesh_dim)
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
    # `CompareTwoPoints` vertex classes upstream: the chapeau's ordinary
    # control copies (tags 3,4) fold into the coincident sources while the
    # `MSH_POINT_BND_LAYER` beg/end copies (5,6) weld the connector dests —
    # the surviving pair carries the displaced positions.
    @test r.model.curves[2] == (5, 6)       # chapeau copy on the top row
    @test r.model.curves[3] == (1, 5)       # connector at the signed begin
    @test r.model.curves[4] == (2, 6)       # connector at the signed end
    @test r.model.loops[only(r.model.surfaces[5])] == [1, 4, -2, -3]
    @test r.model.points[5] == (0.0, -1.0, 0.0)
    @test r.model.points[6] == (1.0, -1.0, 0.0)
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
    @test r.model.points[5] == (0.0, -5.0, 0.0)

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
    @test r.model.points[5] == (0.0, -0.8, 0.0)

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
    @test r.model.points[5] == (1.0, 1.0, 0.0)
    @test r.model.points[6] == (0.0, 1.0, 0.0)

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

    # `Using Index[i]` partitions the vertex merge classes upstream's
    # `CompareTwoPoints` computes from `Typ`/`boundaryLayerIndex`. The
    # index-0 wall's ordinary control copies fold into the coincident
    # sources, while the index-1 wall's identical copies stay coincident
    # orphans; both walls' `MSH_POINT_BND_LAYER` beg/end copies weld their
    # own connector dests and carry the displaced positions. Gmsh 4.15.2's
    # entity store for this fixture is byte-identical: 9 points, 8 curves,
    # 2 surfaces with the same tags, bounds and loops.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {1,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        out[] = Extrude { Curve{1}; Layers{2}; Recombine; };
        out2[] = Extrude { Curve{2}; Layers{2}; Using Index[1]; Recombine; };
        """)
    @test r.errors == String[]
    @test r.lists["out"] == [3.0, 6.0, 5.0, -4.0]
    @test r.lists["out2"] == [7.0, 10.0, 9.0, -8.0]
    @test sort!(collect(keys(r.model.points))) == [1,2,3,6,7,8,9,10,11]
    @test r.model.points[6] == (0.0,-1.0,0.0)   # wall-1 survivors, displaced
    @test r.model.points[7] == (1.0,-1.0,0.0)
    @test r.model.points[8] == (1.0,0.0,0.0)    # index-1 orphans, source pos
    @test r.model.points[9] == (1.0,1.0,0.0)
    @test r.model.points[10] == (2.0,0.0,0.0)   # wall-2 survivors, displaced
    @test r.model.points[11] == (2.0,1.0,0.0)
    @test sort!(collect(keys(r.model.curves))) == [1,2,3,4,5,7,8,9]
    @test r.model.curves[3] == (6,7)
    @test r.model.curves[7] == (10,11)
    @test r.model.curves[4] == (1,6)
    @test r.model.curves[5] == (2,7)
    @test r.model.curves[8] == (2,10)
    @test r.model.curves[9] == (3,11)
    @test r.model.loops[only(r.model.surfaces[6])] == [1,5,-3,-4]
    @test r.model.loops[only(r.model.surfaces[10])] == [2,9,-7,-8]
    @test sort!(collect(keys(r.model.surfaces))) == [6,10]
    # Merge classes: index-0 boundary-layer copies weld only with each
    # other; the index-1 orphans keep `MSH_POINT` type (no boundary-layer
    # class) at their source positions.
    bl=r.model.meshing.point_bl
    @test bl[6] == (true,0) && bl[7] == (true,0)
    @test bl[8] == (false,1) && bl[9] == (false,1)
    @test bl[10] == (true,1) && bl[11] == (true,1)

    # A lone `Using Index[1]` wall shows the same orphan pair — upstream's
    # `DuplicateCurve` control copies carry the index but keep ordinary
    # type, so nothing coincident survives to fold them into.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {1,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        out[] = Extrude { Curve{2}; Layers{2}; Using Index[1]; Recombine; };
        """)
    @test r.errors == String[]
    @test r.lists["out"] == [3.0, 6.0, 5.0, -4.0]
    @test sort!(collect(keys(r.model.points))) == 1:7
    @test r.model.points[4] == (1.0,0.0,0.0)    # coincident orphans
    @test r.model.points[5] == (1.0,1.0,0.0)
    @test r.model.points[6] == (2.0,0.0,0.0)    # displaced survivors
    @test r.model.points[7] == (2.0,1.0,0.0)

    # Boundary-layer extrusion still requires `Layers` and rejects the
    # parameters upstream does not honor on this mode. `Surface{1}` itself
    # is supported now — the empty `Plane Surface(1) = {}` keeps this row on
    # the diagnostic path it always exercised.
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

@testset ".geo boundary-layer curved generatrix" begin
    # A `Circle` generatrix is re-graded by the native 1-D producer between
    # the boundary layer's refresh-time parameter sweep and the chapeau's
    # own emission: the same logical wall node re-evaluates a few e-9 along
    # the curve. The chapeau's top-row lookup has to recover that
    # correspondence (the drift stays far inside the key's Voronoi cell) —
    # upstream's `copyMesh` maps the stored source *vertices* themselves,
    # so the row binding is by node identity, not bitwise position.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {1,1,0,1};
        Circle(1) = {1,2,3};
        out[] = Extrude { Curve{1}; Layers{3}; Recombine; };
        Mesh 1;
        """)
    @test r.errors == String[]
    @test r.lists["out"] == [2.0, 5.0, 4.0, -3.0]

    # The generatrix keeps two Line cells; its regraded interior node is
    # the arc midpoint Tessella's native grading produces (upstream emits
    # the same point).
    gen = geo_entity_mesh(r, 1, 1)
    @test sum(size(b.nodes, 2) for b in gen.blocks) == 2
    @test Set(Tuple(round.(gen.coords[:, i]; digits=9))
              for i in axes(gen.coords, 2)) ==
          Set([(0.0, 0.0, 0.0), (0.292893219, 0.707106781, 0.0),
               (1.0, 1.0, 0.0)])

    # Chapeau: one Line cell per generatrix element on the welded top row.
    # The unique positions match Gmsh 4.15.2's `.msh` emission — the oracle
    # writes `0.9238795320849413 -0.382683433394379` for the begin top and
    # `1.000000000787784 7.877838e-10` for the interior top (its smoothed
    # endpoint normals accumulate in the same order).
    chapeau = geo_entity_mesh(r, 1, 2)
    @test sum(size(b.nodes, 2) for b in chapeau.blocks) == 2
    @test Set(Tuple(round.(chapeau.coords[:, i]; digits=9))
              for i in axes(chapeau.coords, 2)) ==
          Set([(0.923879532, -0.382683433, 0.0), (1.000000001, 1e-9, 0.0),
               (1.382683431, 0.076120467, 0.0)])

    # Each connector emits one segment per layer along the smoothed end
    # normal, sharing the strip's column nodes.
    c3 = geo_entity_mesh(r, 1, 3)
    @test sum(size(b.nodes, 2) for b in c3.blocks) == 3
    @test Set(Tuple(round.(c3.coords[:, i]; digits=9))
              for i in axes(c3.coords, 2)) ==
          Set([(0.0, 0.0, 0.0), (0.307959844, -0.127561144, 0.0),
               (0.615919688, -0.255122289, 0.0),
               (0.923879532, -0.382683433, 0.0)])
    c4 = geo_entity_mesh(r, 1, 4)
    @test sum(size(b.nodes, 2) for b in c4.blocks) == 3
    @test Set(Tuple(round.(c4.coords[:, i]; digits=9))
              for i in axes(c4.coords, 2)) ==
          Set([(1.0, 1.0, 0.0), (1.127561144, 0.692040156, 0.0),
               (1.255122288, 0.384080311, 0.0),
               (1.382683431, 0.076120467, 0.0)])
end

@testset "boundary-layer top-row lookup tolerance" begin
    # `_extrude_at`'s `:boundary_layer` arm binds an evaluated generatrix
    # position to its stored top node. Exact keys hit bitwise; a few e-9
    # of re-evaluation drift is recovered only when the nearest row is
    # unambiguous (an order of magnitude closer than the runner-up); a
    # genuinely foreign point between rows still rejects.
    spec = (type=:boundary_layer,
            tops=Dict((0.0, 0.0, 0.0) => (0.0, -1.0, 0.0),
                      (1.0, 0.0, 0.0) => (1.0, -1.0, 0.0)))
    @test Tessella.Model._extrude_at(spec, (0.0, 0.0, 0.0), 0.5, "test") ==
          (0.0, -1.0, 0.0)
    # The strict 1e-9·scale arm.
    @test Tessella.Model._extrude_at(spec, (1e-10, 0.0, 0.0), 0.5, "test") ==
          (0.0, -1.0, 0.0)
    # The observed ~2.2e-9 re-grading drift recovers through the
    # Voronoi-dominance arm (nearest key ~5e5× closer than the runner-up).
    @test Tessella.Model._extrude_at(spec, (2.2e-9, 0.0, 0.0), 0.5, "test") ==
          (0.0, -1.0, 0.0)
    # Drift past the 1e-7·scale cap still rejects even when unambiguous.
    @test_throws ArgumentError Tessella.Model._extrude_at(
        spec, (1e-6, 0.0, 0.0), 0.5, "test")
    # A point between keys (best ≈ runner-up distance) rejects.
    @test_throws ArgumentError Tessella.Model._extrude_at(
        spec, (0.5, 0.0, 0.0), 0.5, "test")
    @test_throws ArgumentError Tessella.Model._extrude_at(
        spec, (0.5, 5.0, 0.0), 0.5, "test")
    @test_throws ArgumentError Tessella.Model._extrude_at(
        (type=:boundary_layer, tops=Dict{NTuple{3,Float64},
                                       NTuple{3,Float64}}()),
        (0.0, 0.0, 0.0), 0.5, "test")
end

const _BL_VIEW_POS = """
View "wall" {
VP(0.0,0.0,0.0) {0.0,2.0,0.0};
VP(1.0,0.0,0.0) {0.0,2.0,0.0};
};
"""

@testset "boundary-layer Extrude Using View[i]" begin
    # `Using View[0]` resolves the merged `.pos` view's raw interpolated
    # vectors as the growth direction — upstream `OctreePost::searchVector`
    # feeds `x += n[i]*t` with the magnitude included, so a (0,2,0) field
    # grows the strip to y=+2 (the gouraud side is y=-1, and the raw vector
    # points opposite the wall normal — the ordinary positive-Jacobian gate
    # does not apply to view-driven strips).
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Merge "wall_view.pos";
        out[] = Extrude { Curve{1}; Layers{2}; Using View[0]; Recombine; };
        Mesh 2;
        """, files=Dict("wall_view.pos"=>_BL_VIEW_POS))
    @test r.errors == String[]
    @test r.lists["out"] == [2.0, 5.0, 4.0, -3.0]
    @test r.model.points[5] == (0.0, 2.0, 0.0)
    @test r.model.points[6] == (1.0, 2.0, 0.0)
    surface = geo_entity_mesh(r, 2, 5)
    @test _bl_cells(surface) == (tris=0, quads=2)
    ys = sort!(unique!(round.(vec(surface.coords[2, :]); digits=12)))
    @test ys ≈ [0.0, 1.0, 2.0]
    conn = geo_entity_mesh(r, 1, 3)
    @test sum(size(b.nodes, 2) for b in conn.blocks) == 2
    @test sort!(vec(conn.coords[2, :])) ≈ [0.0, 1.0, 2.0]

    # Out-of-plane view vectors displace the strip off the wall plane —
    # upstream emits the warped row without a planarity gate (oracle tops
    # (0,1,1),(1,1,1) for a (0,1,1) field).
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Merge "warp.pos";
        out[] = Extrude { Curve{1}; Layers{1}; Using View[0]; Recombine; };
        """,
        files=Dict("warp.pos"=>"""
        View "warp" {
        VP(0.0,0.0,0.0) {0.0,1.0,1.0};
        VP(1.0,0.0,0.0) {0.0,1.0,1.0};
        };
        """))
    @test r.errors == String[]
    @test r.model.points[5] == (0.0, 1.0, 1.0)
    @test r.model.points[6] == (1.0, 1.0, 1.0)

    # A view carrying no vector records (scalar payload) leaves the wall on
    # the gouraud field — upstream's `getNumVectors()>0` switch.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Merge "scalar.pos";
        out[] = Extrude { Curve{1}; Layers{1}; Using View[0]; Recombine; };
        """,
        files=Dict("scalar.pos"=>"""
        View "scalar" {
        SP(0.0,0.0,0.0) {5.0};
        SP(1.0,0.0,0.0) {7.0};
        };
        """))
    @test r.errors == String[]
    @test r.model.points[5] == (0.0, -1.0, 0.0)
    @test r.model.points[6] == (1.0, -1.0, 0.0)

    # `Using View[9]` is a recoverable upstream `Msg::Error`; the wall falls
    # back to the gouraud field and `Mesh 2` still produces.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{1}; Using View[9]; Recombine; };
        Mesh 2;
        """)
    @test any(occursin("Unknown View[9]: using normals instead", e)
              for e in r.errors)
    @test r.model.points[5] == (0.0, -1.0, 0.0)
    @test r.model.points[6] == (1.0, -1.0, 0.0)

    # `View[-3]` snaps each normalized component with a 0.1 deadband — the
    # 45-degree wall's (±0.707,∓0.707,0) normals land on (±1,∓1,0), matching
    # the oracle tops (1,-1,0) and (2,0,0).
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,1,0,1};
        Line(1) = {1,2};
        out[] = Extrude { Curve{1}; Layers{1}; Using View[-3]; Recombine; };
        """)
    @test r.errors == String[]
    @test all(r.model.points[5] .≈ (1.0, -1.0, 0.0))
    @test all(r.model.points[6] .≈ (2.0, 0.0, 0.0))

    # `View[-5]` scales the normalized field by the background field value —
    # `F=2` doubles the (0,-1,0) normal, oracle tops (0,-2),(1,-2).
    r = _execute_bl_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Field[1] = MathEval;
        Field[1].F = "2";
        Background Field = 1;
        out[] = Extrude { Curve{1}; Layers{1}; Using View[-5]; Recombine; };
        """)
    @test r.errors == String[]
    @test all(r.model.points[5] .≈ (0.0, -2.0, 0.0))
    @test all(r.model.points[6] .≈ (1.0, -2.0, 0.0))
end

@testset ".geo boundary-layer 3-D surface Extrude" begin
    src = """
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {1,1,0,1};
        Point(4) = {0,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,4};
        Line(4) = {4,1};
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        """
    # `Extrude {Surface{1}; Layers{3}; Recombine}` — the 3-D boundary layer
    # clones the face's edges (`DuplicateCurve` top-row copies 6..9), sweeps
    # per-point connector columns (11,12,16,20), one `MSH_SURF_BND_LAYER`
    # strip per source edge (13,17,21,25), a delayed-rename chapeau (26) and
    # the volume (1) — the exact NEWREG order the oracle emits.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{3}; Recombine; };
        """)
    @test r.errors == String[]
    @test r.lists["out"] == [26.0, 1.0, 13.0, 17.0, 21.0, 25.0]
    m = r.model
    # Upstream leaves the chapeau's geometry unextruded until mesh time —
    # `-0 -o *.geo_unrolled` writes `Point(7) = {0,0,0,…}` and fails to
    # serialize the internal boundary-layer curve type — so the stored
    # coordinates stay on the source plane.
    @test m.points[7] == (0.0, 0.0, 0.0)
    @test m.points[8] == (1.0, 0.0, 0.0)
    @test m.points[12] == (1.0, 1.0, 0.0)
    @test m.points[16] == (0.0, 1.0, 0.0)
    @test m.curves[6] == (7, 8)
    @test m.curves[7] == (8, 12)
    @test m.curves[8] == (12, 16)
    @test m.curves[9] == (16, 7)
    @test m.curves[11] == (1, 7)
    @test m.curves[12] == (2, 8)
    @test m.curves[16] == (3, 12)
    @test m.curves[20] == (4, 16)
    @test m.loops[only(m.surfaces[13])] == [1, 12, -6, -11]
    @test m.loops[only(m.surfaces[17])] == [2, 16, -7, -12]
    @test m.loops[only(m.surfaces[21])] == [3, 20, -8, -16]
    @test m.loops[only(m.surfaces[25])] == [4, 11, -9, -20]
    @test m.loops[only(m.surfaces[26])] == [6, 7, 8, 9]
    @test m.surface_loops[only(m.volumes[1])] == [-1, 26, 13, 17, 21, 25]
    # `geo.Source` records — signed source on the chapeau/volume, the
    # generatrix edge on each lateral.
    @test m.meshing.extrude_sources[(3, 1)] == (2, 1)
    @test m.meshing.extrude_sources[(2, 26)] == (2, 1)
    @test m.meshing.extrude_sources[(2, 13)] == (1, 1)
    @test get(m.meshing.extrude_specs[(2, 26)], :type, nothing) in
        (:bl_surface, :boundary_layer)
    @test Int(get(m.meshing.extrude_specs[(2, 13)], :surf, 0)) == 1

    # `Mesh 3`: the recombined product matches the oracle layout — the
    # source face's native mesh carries a center fan (its interior node
    # grows its own displaced column), each lateral emits one quadrangle
    # per layer row and the swept volume is twelve prisms.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{3}; Recombine; };
        """; mesh_dim=3)
    mesh = r.mesh
    @test size(mesh.coords, 2) == 20
    cells = (tris = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 2; init = 0),
             quads = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 3; init = 0),
             prisms = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 6; init = 0),
             lines = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 1; init = 0))
    @test cells == (tris = 8, quads = 12, prisms = 12, lines = 20)
    want = Set{NTuple{3,Float64}}()
    for corner in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5))
        for z in (0.0, 1/3, 2/3, 1.0)
            push!(want, (corner[1], corner[2], z))
        end
    end
    got = Set(Tuple(round.(mesh.coords[:, i]; digits=12))
              for i in axes(mesh.coords, 2))
    want = Set(Tuple(round.(collect(p); digits=12)) for p in want)
    @test got == want
    for t in (13, 17, 21, 25)
        @test _bl_cells(geo_entity_mesh(r, 2, t)) == (tris=0, quads=3)
    end
    @test _bl_cells(geo_entity_mesh(r, 2, 26)) == (tris=4, quads=0)

    # Without `Recombine` the strips split on the lateral diagonal and the
    # volume subdivides into tetrahedra — the oracle's 20 nodes, 20 line
    # segments, 32 triangles and 36 tetrahedra.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{3}; };
        """; mesh_dim=3)
    mesh = r.mesh
    @test size(mesh.coords, 2) == 20
    @test size(mesh.segs, 2) == 20
    @test size(mesh.tris, 2) == 32
    @test size(mesh.tets, 2) == 36

    # The signed source inverts the growth direction; upstream's lateral
    # `out` filter compares `abs(tag)` against the *signed* input, so the
    # source face itself lists between the body and the side surfaces.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{-1}; Layers{3}; Recombine; };
        """; mesh_dim=3)
    @test r.lists["out"] == [26.0, 1.0, 1.0, 13.0, 17.0, 21.0, 25.0]
    @test r.model.surface_loops[only(r.model.volumes[1])] ==
        [-1, 26, 13, 17, 21, 25]
    mesh = r.mesh
    @test size(mesh.coords, 2) == 20
    @test sort!(unique!(round.(mesh.coords[3, :]; digits=12))) ≈
        [-1.0, -2/3, -1/3, 0.0]
    cells = (quads = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 3; init = 0),
             prisms = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 6; init = 0))
    @test cells == (quads = 12, prisms = 12)

    # `ScaleLastLayer` rescales every element level inside the last layer
    # group by the vertex's incident mean element edge length — for the
    # center-fan face mesh every vertex measures ~0.8047379, the oracle's
    # exact level set.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{3}; ScaleLastLayer; Recombine; };
        """; mesh_dim=3)
    zs = sort!(unique!(round.(r.mesh.coords[3, :]; digits=9)))
    @test zs ≈ [0.0, 0.268245951, 0.536491903, 0.804737854]

    # `Using Index[i]` keeps the same layout — the index only partitions
    # which contributions accumulate into the normal field.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{3}; Using Index[0]; Recombine; };
        """; mesh_dim=3)
    @test r.lists["out"] == [26.0, 1.0, 13.0, 17.0, 21.0, 25.0]
    @test size(r.mesh.coords, 2) == 20
    @test sort!(unique!(round.(r.mesh.coords[3, :]; digits=12))) ≈
        [0.0, 1/3, 2/3, 1.0]

    # `Using View[i]` samples the post view's vector records instead of the
    # gouraud field — a uniform (0,0,2) field doubles the growth direction.
    # Two `VT` vector triangles tile the source face so every sample point
    # lies inside an element; upstream's `searchVector` interpolates the
    # (constant) nodal vectors rather than falling back to the zero-vector
    # miss that keeps a vertex on its source row.
    r = _execute_bl_source(src * """
        Merge "vec.pos";
        out[] = Extrude { Surface{1}; Layers{3}; Using View[0]; Recombine; };
        """;
        files=Dict("vec.pos"=>"""
        View "vec" {
        VT(0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0)
            {0.0,0.0,2.0, 0.0,0.0,2.0, 0.0,0.0,2.0};
        VT(0.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0)
            {0.0,0.0,2.0, 0.0,0.0,2.0, 0.0,0.0,2.0};
        };
        """), mesh_dim=3)
    @test r.errors == String[]
    @test sort!(unique!(round.(r.mesh.coords[3, :]; digits=12))) ≈
        [0.0, 2/3, 4/3, 2.0]
end

# The lateral strips' plane from the executed mesh — used to pin each
# canonical tag to its geometric edge.
function _bl_strip_plane(r, tag)
    part = geo_entity_mesh(r, 2, tag)
    isempty(part.coords) && return nothing
    return extrema(part.coords, dims = 2)
end

@testset ".geo boundary-layer holed surface Extrude" begin
    # A holed source face canonicalizes its wires once in
    # `SetSurfaceGeneratrices`: a positive hole reference enumerates its
    # curves reversed and negated, so `Extrude{Surface}` allocates the
    # hole's chapeau copies, `geo.Source` links and lateral tags in the
    # canonical `-8,-7,-6,-5` order rather than the declared `{5,6,7,8}`.
    # Gmsh 4.15.2 `-0` on this fixture emits exactly:
    #   Curve Loop(37) = {-8, 36, -14, -35};  Curve Loop(41) = {-7, 40, -15, -36};
    #   Curve Loop(45) = {-6, 44, -16, -40};  Curve Loop(49) = {-5, 35, -17, -44};
    #   Surface Loop(1) = {1, 50, 21, 25, 29, 33, 37, 41, 45, 49};
    # with `SurfacesOrientations` -1 on the four hole laterals.
    src = """
        Point(1) = {0,0,0,1};
        Point(2) = {4,0,0,1};
        Point(3) = {4,4,0,1};
        Point(4) = {0,4,0,1};
        Point(5) = {1,1,0,0.5};
        Point(6) = {3,1,0,0.5};
        Point(7) = {3,3,0,0.5};
        Point(8) = {1,3,0,0.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
        Curve Loop(1) = {1,2,3,4};
        Curve Loop(2) = {5,6,7,8};
        Plane Surface(1) = {1,2};
        """
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{2}; Recombine; };
        """)
    @test r.errors == String[]
    @test r.lists["out"] == [50.0, 1.0, 21.0, 25.0, 29.0, 33.0, 37.0, 41.0,
                             45.0, 49.0]
    m = r.model
    # Canonical hole laterals — the signed generatrix record, its end
    # connector, the negated chapeau copy, and the negated beg connector.
    @test m.loops[only(m.surfaces[21])] == [1, 20, -10, -19]
    @test m.loops[only(m.surfaces[25])] == [2, 24, -11, -20]
    @test m.loops[only(m.surfaces[29])] == [3, 28, -12, -24]
    @test m.loops[only(m.surfaces[33])] == [4, 19, -13, -28]
    @test m.loops[only(m.surfaces[37])] == [-8, 36, -14, -35]
    @test m.loops[only(m.surfaces[41])] == [-7, 40, -15, -36]
    @test m.loops[only(m.surfaces[45])] == [-6, 44, -16, -40]
    @test m.loops[only(m.surfaces[49])] == [-5, 35, -17, -44]
    # The copy set itself lands on the canonical edges: curves 14..17 are
    # the reversed records' duplicates — copies of -8,-7,-6,-5 traversing
    # the hole clockwise. Their wired endpoints are the `MSH_POINT`-class
    # orphans' siblings — the beg/end `MSH_POINT_BND_LAYER` copies weld
    # per vertex under `CompareTwoPoints`'s type-plus-index merge key, so
    # 27,28,32,36 are the survivors (the control orphans stay at
    # 25,26,29,30,33,34,37,38 and never join the curve's wired pair).
    @test m.curves[14] == (27, 28)
    @test m.curves[15] == (28, 32)
    @test m.curves[16] == (32, 36)
    @test m.curves[17] == (36, 27)
    @test m.points[27] == (1.0, 1.0, 0.0)
    @test m.points[28] == (1.0, 3.0, 0.0)
    @test m.points[32] == (3.0, 3.0, 0.0)
    @test m.points[36] == (3.0, 1.0, 0.0)
    # Signed `geo.Source` links follow the canonical entries, and the
    # chapeau/volume carry the signed source face.
    @test m.meshing.extrude_sources[(2, 37)] == (1, -8)
    @test m.meshing.extrude_sources[(2, 41)] == (1, -7)
    @test m.meshing.extrude_sources[(2, 45)] == (1, -6)
    @test m.meshing.extrude_sources[(2, 49)] == (1, -5)
    @test m.meshing.extrude_sources[(3, 1)] == (2, 1)
    @test m.meshing.extrude_sources[(2, 50)] == (2, 1)
    # The volume's shell lists every lateral in canonical order — the
    # hole's `SurfacesOrientations` entries are the negative tags.
    @test m.surface_loops[only(m.volumes[1])] ==
        [-1, 50, 21, 25, 29, 33, -37, -41, -45, -49]
    # The chapeau's hole wire stores its declared traversal — the copies
    # reversed and negated — so a later generatrix pass reconstructs the
    # upstream `Curve Loop(50) = {10,...,17}` flat order.
    @test length(m.surfaces[50]) == 2
    @test m.loops[m.surfaces[50][2]] == [-17, -16, -15, -14]

    # `Mesh 3`: the oracle's 156 nodes, 144 triangles (72+72 on the source
    # and displaced chapeau), 64 lateral quadrangles and 144 swept prisms.
    # Hole strips point their normals into the solid annulus while the
    # outer strips point outward, matching the oracle block-for-block.
    r = _execute_bl_source(src * """
        out[] = Extrude { Surface{1}; Layers{2}; Recombine; };
        """; mesh_dim=3)
    mesh = r.mesh
    @test r.errors == String[]
    @test size(mesh.coords, 2) == 156
    cells = (tris = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 2; init = 0),
             quads = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 3; init = 0),
             prisms = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 6; init = 0))
    @test cells == (tris = 144, quads = 64, prisms = 144)
    @test sort!(unique!(round.(mesh.coords[3, :]; digits=12))) ≈ [0.0, 0.5, 1.0]
    for t in (21, 25, 29, 33, 37, 41, 45, 49)
        @test _bl_cells(geo_entity_mesh(r, 2, t)) == (tris=0, quads=8)
    end
    # Each canonical tag binds to its oracle edge plane: the hole strips
    # sit on the inner square's left (x=1), top (y=3), right (x=3) and
    # bottom (y=1) edges in tag order — the declared order would instead
    # put them on y=1, x=3, y=3, x=1.
    plane = _bl_strip_plane(r, 37)
    @test abs(plane[1][1] - 1.0) < 1e-9 && abs(plane[1][2] - 1.0) < 1e-9
    plane = _bl_strip_plane(r, 41)
    @test abs(plane[2][1] - 3.0) < 1e-9 && abs(plane[2][2] - 3.0) < 1e-9
    plane = _bl_strip_plane(r, 45)
    @test abs(plane[1][1] - 3.0) < 1e-9 && abs(plane[1][2] - 3.0) < 1e-9
    plane = _bl_strip_plane(r, 49)
    @test abs(plane[2][1] - 1.0) < 1e-9 && abs(plane[2][2] - 1.0) < 1e-9
end

@testset ".geo boundary-layer native grading at small lc" begin
    # BL preparation seeds each generatrix and source boundary curve with
    # the native 1-D producer — under the caller's `_ModelMesh1DOptions` —
    # so the position-keyed `levels`/`tops` maps store exactly the
    # coordinates the curve producer later emits. Before the seeding the
    # lazy grading fallback drifted a few ulps at this `lc`: ownership then
    # rejected a lower cell whose stored row no longer matched the emitted
    # node, and `Mesh 3` failed inside certificate validation. The unit
    # square at `lc = 0.5` exercises the drift; Gmsh 4.15.2 meshes it to
    # the counts asserted below.
    r = _execute_bl_source("""
        Point(1) = {0,0,0,0.5};
        Point(2) = {1,0,0,0.5};
        Point(3) = {1,1,0,0.5};
        Point(4) = {0,1,0,0.5};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,4};
        Line(4) = {4,1};
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        out[] = Extrude { Surface{1}; Layers{3}; Recombine; };
        """; mesh_dim=3)
    @test r.errors == String[]
    mesh = r.mesh
    @test size(mesh.coords, 2) == 48
    cells = (tris = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 2; init = 0),
             quads = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 3; init = 0),
             prisms = sum(size(b.nodes, 2) for b in mesh.blocks if b.msh == 6; init = 0))
    @test cells == (tris = 28, quads = 24, prisms = 42)
    @test sort!(unique!(round.(mesh.coords[3, :]; digits=12))) ≈
        [0.0, 1/3, 2/3, 1.0]
end
