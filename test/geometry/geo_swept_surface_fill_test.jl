# The suite includes every test file into Main; this file keeps its names to itself.
module SweptSurfaceFillTests
using Test
using Tessella
import Tessella.Model as M
import Tessella.Mesh2D as K
const A = Tessella.API

# A disk swept by one quarter turn of a screw, without Layers: the first step
# of the coil of Gmsh's examples/simple_geo/indheat.geo. Its four lateral
# faces are long strips (0.016 by 0.16) between a quarter arc, its swept copy
# and two helices, which Gmsh stores as splines.
const _SCREW = """
lc = 0.02; r = 0.11; rc = 0.01; hc = 0.25; turns = 5;
Point(1) = {0, -r, -hc/2, lc};
Point(2) = {0, -r+rc, -hc/2, lc};
Point(3) = {0, -r, -hc/2+rc, lc};
Point(4) = {0, -r-rc, -hc/2, lc};
Point(5) = {0, -r, -hc/2-rc, lc};
Circle(1) = {2,1,3}; Circle(2) = {3,1,4}; Circle(3) = {4,1,5}; Circle(4) = {5,1,2};
Curve Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
tmp[] = Extrude { {0,0,hc/Floor(turns)/4}, {0,0,1}, {0,0,0}, Pi/2 } { Surface{6}; };
"""

# A quarter cylinder whose arc ends 1e-5 beside its last point, as written
# coordinates do (examples/simple_geo/piece.geo): Gmsh draws the arc with the
# radius of its first point and accepts the other up to 10% off.
const _LOOSE_ARC = """
Point(1) = {0, 0, 0, 0.2}; Point(2) = {1, 0, 0, 0.2}; Point(3) = {0, 1.00001, 0, 0.2};
Circle(1) = {2,1,3};
out[] = Extrude {0, 0, 0.5} { Curve{1}; };
"""

_count(dim, tag) = (blocks = A.mesh.get_elements(dim, tag)[2]; isempty(blocks) ? 0 : sum(length, blocks))
_count(dim) = (blocks = A.mesh.get_elements(dim)[2]; isempty(blocks) ? 0 : sum(length, blocks))

function _with_geo(body, text)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        mktempdir() do directory
            path = joinpath(directory, "case.geo"); write(path, text)
            A.open(path)
            body()
        end
    finally
        A.finalize()
    end
end

@testset "A screw-swept solid without Layers meshes as in Gmsh" begin
    _with_geo(_SCREW) do
        laterals = [tag for (_, tag) in A.model.get_entities(2) if A.model.get_type(2, tag) != "Plane"]
        caps = [tag for (_, tag) in A.model.get_entities(2) if A.model.get_type(2, tag) == "Plane"]
        @test length(laterals) == 4 && length(caps) == 2
        A.mesh.generate(2)
        # The curves carry Gmsh's nodes: one on each quarter arc, 7 to 9 on the helices.
        interior = Dict(tag => length(A.mesh.get_nodes(1, tag)[1]) for (_, tag) in A.model.get_entities(1))
        @test sort!(collect(values(interior))) == [1, 1, 1, 1, 1, 1, 1, 1, 7, 8, 8, 9]
        # A helix is a spline through points of the true helix: its nodes lie
        # beside the swept surface, so the surface cannot be asked to invert
        # them, and Gmsh does not ask.
        helix = first(tag for (tag, n) in interior if n == 7)
        _, xyz, _ = A.mesh.get_nodes(1, helix)
        @test maximum(abs(hypot(xyz[3k-2], xyz[3k-1]) - 0.1) for k in 1:7) > 1e-5
        # The caps are Gmsh's; every strip gets interior points (Gmsh: 29, 31, 29, 29 triangles, 93 nodes).
        @test all(_count(2, tag) == 14 for tag in caps)
        # The strips are meshed as Gmsh meshes them (PolyMesh initial mesh,
        # delaunayizeBDS, frontal points on the surface): the same counts.
        @test all(!isempty(A.mesh.get_nodes(2, tag)[1]) for tag in laterals)
        if Sys.ARCH === :x86_64
            @test [_count(2, tag) for tag in laterals] == [29, 31, 29, 29]
            @test length(A.mesh.get_nodes()[1]) == 93
        else
            @test all(27 <= _count(2, tag) <= 33 for tag in laterals)
            @test abs(length(A.mesh.get_nodes()[1]) - 93) <= 3
        end
        surface_cells = _count(2)
        A.mesh.generate(3)
        # The volume fills from that boundary and leaves it alone (Gmsh: 151 tetrahedra).
        @test _count(2) == surface_cells
        Sys.ARCH === :x86_64 ? (@test _count(3) == 151) : (@test 130 <= _count(3) <= 170)
    end
end

@testset "A swept surface resolved once gives the points it gives one by one" begin
    _with_geo(_SCREW) do
        A._with_model() do m
            for (_, tag) in A.model.get_entities(2)
                A.model.get_type(2, tag) == "Plane" && continue
                lift = M._ruled_surface_evaluator(m, tag, "test")
                @test all(lift(u, v) === M._ruled_surface_point(m, tag, u, v, "test")
                          for u in 0.0:0.125:1.0, v in 0.0:0.125:1.0)
                if Base.JLOptions().compile_enabled == 1
                    probe(f) = @allocated f(0.3, 0.4)
                    probe(lift)
                    @test probe(lift) == 0
                end
            end
        end
    end
    # A patch that is not a rigid sweep (its corner moved) keeps the blend.
    _with_geo(_LOOSE_ARC) do
        A._with_model() do m
            tag = only(tag for (_, tag) in A.model.get_entities(2))
            @test M._ruled_exact_extrusion(m, tag, m.loops[first(m.surfaces[tag])], "test") === nothing
            lift = M._ruled_surface_evaluator(m, tag, "test")
            @test all(lift(u, v) === M._ruled_surface_point(m, tag, u, v, "test")
                      for u in 0.0:0.25:1.0, v in 0.0:0.25:1.0)
        end
    end
end

# Every built-in curve kind the 1-D mesher asks points of.
const _CURVES = """
Point(1) = {0, 0, 0}; Point(2) = {1, 0, 0}; Point(3) = {0, 1, 0}; Point(4) = {-1, 0, 0};
Point(5) = {2, 0.5, 0.1}; Point(6) = {3, -0.5, 0.3}; Point(7) = {4, 0.2, 0}; Point(8) = {5, 0, 0.4};
Point(9) = {2, 0, 0}; Point(10) = {0, 0.5, 0};
Line(1) = {1, 2};
Circle(2) = {2, 1, 3};
Ellipse(3) = {9, 1, 9, 10};
Spline(4) = {2, 5, 6, 7, 8};
BSpline(5) = {2, 5, 6, 7, 8};
Bezier(6) = {2, 5, 6, 7, 8};
"""

@testset "A curve resolved once gives the points it gives one by one" begin
    _with_geo(_CURVES) do
        A._with_model() do m
            for tag in 1:6
                point = M._model_curve_point_evaluator(m, tag, "test")
                t0, t1 = M._model_curve_param_bounds(m, tag, "test")
                @test all(point(t0 + k * (t1 - t0) / 16) === M._model_curve_point(m, tag, t0 + k * (t1 - t0) / 16, "test")
                          for k in 0:16)
            end
            if Base.JLOptions().compile_enabled == 1
                # A spline asked point by point collects its control points each time.
                probe(f) = @allocated f(0.37)
                resolved = M._model_curve_point_evaluator(m, 4, "test")
                general = u -> M._model_curve_point(m, 4, u, "test")
                probe(resolved); probe(general)
                @test probe(resolved) == 0
                @test probe(general) > 0
            end
        end
    end
end

@testset "An arc that ends beside its last point bounds a meshable surface" begin
    _with_geo(_LOOSE_ARC) do
        tag = only(tag for (_, tag) in A.model.get_entities(2))
        A._with_model() do m
            # The arc's end and its end point differ by the radius mismatch.
            t0, t1 = M._model_curve_param_bounds(m, 1, "test")
            stop = M._model_curve_point(m, 1, t1, "test")
            @test 5e-6 < hypot((stop .- m.points[3])...) < 2e-5
        end
        A.mesh.generate(2)
        @test _count(2, tag) > 0
        # The mesh node of Point 3 is the point, not the arc's end.
        _, xyz, _ = A.mesh.get_nodes(0, 3)
        @test xyz == [0.0, 1.00001, 0.0]
    end
end

# The sphere of Gmsh's examples/simple_geo/sphere.geo: eight three-sided ruled
# patches between quarter arcs.
const _SPHERE = """
Mesh.Algorithm = 6;
lc = .55;
Point(1) = {0,0,0,lc}; Point(2) = {1,0,0,lc}; Point(3) = {0,1,0,lc};
Point(4) = {-1,0,0,lc}; Point(5) = {0,-1,0,lc}; Point(6) = {0,0,-1,lc}; Point(7) = {0,0,1,lc};
Circle(1) = {2,1,3}; Circle(2) = {3,1,4}; Circle(3) = {4,1,5}; Circle(4) = {5,1,2};
Circle(5) = {3,1,6}; Circle(6) = {6,1,5}; Circle(7) = {5,1,7}; Circle(8) = {7,1,3};
Circle(9) = {2,1,7}; Circle(10) = {7,1,4}; Circle(11) = {4,1,6}; Circle(12) = {6,1,2};
Curve Loop(13) = {2,8,-10}; Surface(14) = {13};
Curve Loop(15) = {10,3,7}; Surface(16) = {15};
Curve Loop(17) = {-8,-9,1}; Surface(18) = {17};
Curve Loop(19) = {-11,-2,5}; Surface(20) = {19};
Curve Loop(21) = {-5,-12,-1}; Surface(22) = {21};
Curve Loop(23) = {-3,11,6}; Surface(24) = {23};
Curve Loop(25) = {-7,4,9}; Surface(26) = {25};
Curve Loop(27) = {-4,12,-6}; Surface(28) = {27};
"""

@testset "A scalar size on a ruled patch is measured on the surface" begin
    _with_geo(_SPHERE) do
        A.mesh.generate(2)
        # Gmsh: 112 triangles and 59 nodes, 13 or 15 triangles a patch. Reading
        # lengths through the patch's tensor at one point of each cell gave 152.
        # The reference is Gmsh's Windows build; Tessella follows each
        # platform's own libm (GmshLibm), and with glibc one patch takes 15.
        if Sys.iswindows() && Sys.ARCH === :x86_64
            @test [_count(2, tag) for (_, tag) in A.model.get_entities(2)] ==
                  [13, 15, 15, 15, 13, 13, 15, 13]
            @test length(A.mesh.get_nodes()[1]) == 59
        else
            @test 104 <= _count(2) <= 120
            @test abs(length(A.mesh.get_nodes()[1]) - 59) <= 4
            @test all(12 <= _count(2, tag) <= 16 for (_, tag) in A.model.get_entities(2))
        end
        # Every node is on the unit sphere, but for the arcs' center.
        _, xyz, _ = A.mesh.get_nodes()
        radii = [hypot(xyz[3k-2], xyz[3k-1], xyz[3k]) for k in 1:length(xyz)÷3]
        @test count(iszero, radii) == 1
        @test all(r -> iszero(r) || abs(r - 1) < 1e-12, radii)
    end
end

@testset "The initial triangulation is made Delaunay in the surface's metric" begin
    # The parameters of a strip ten times longer than wide, mapped to the unit
    # square: one node in the middle of each short side, seven along each long one.
    xs = Float64[]; ys = Float64[]
    for u in (0.0, 0.5, 1.0); push!(xs, u); push!(ys, 0.0); end
    for k in 1:7; push!(xs, 1.0); push!(ys, k / 8); end
    for u in (1.0, 0.5, 0.0); push!(xs, u); push!(ys, 1.0); end
    for k in 7:-1:1; push!(xs, 0.0); push!(ys, k / 8); end
    n = length(xs)
    segments = [(k, mod1(k + 1, n)) for k in 1:n]
    strip(u, v) = (1.0, 0.0, 100.0)            # E, F, G of (u, v) -> (u, 10v)
    function longest(T, interior)
        # the longest cell edge, measured on the strip
        worst = 0.0
        for t in eachindex(T.alive)
            (T.alive[t] && interior[t] && !K._is_ghost_tri(T, t)) || continue
            for k in 1:3
                a, b = K._edge(T, t, k)
                pa = K._pt(T, a); pb = K._pt(T, b)
                worst = max(worst, hypot(pa[1] - pb[1], 10 * (pa[2] - pb[2])))
            end
        end
        return worst
    end
    function violations(T, interior, ground)
        bad = 0
        for t in eachindex(T.alive)
            (T.alive[t] && interior[t] && !K._is_ghost_tri(T, t)) || continue
            for k in 1:3
                other = K._nbr(T, t, k)
                (other == 0 || !interior[other] || K._is_ghost_tri(T, other)) && continue
                a, b = K._edge(T, t, k)
                K._constrained_edge(T, K._skey(a, b)) && continue
                p = K._vert(T, t, k); s = K._vert(T, other, K._nslot(T, other, t))
                pa = K._pt(T, a); pb = K._pt(T, b); pp = K._pt(T, p); ps = K._pt(T, s)
                metric = ground(0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2]))
                center, radius = K._gmsh_circumcenter_metric(pa, pb, pp, metric)
                inside = K._gmsh_metric_distance2(center, ps, metric) < radius * (1 - 1e-9)
                convex = Tessella.Predicates.orient2(pp, pa, ps) * Tessella.Predicates.orient2(pp, ps, pb) > 0
                inside && convex && (bad += 1)
            end
        end
        return bad
    end
    T = K.constrained_delaunay(xs, ys, segments)
    interior = Vector{Bool}(collect(K.classify_interior(T)))
    cells = count(t -> T.alive[t] && interior[t] && !K._is_ghost_tri(T, t), eachindex(T.alive))
    @test cells == n - 2
    @test longest(T, interior) > 3        # flat cells reach along the strip
    @test violations(T, interior, strip) > 0
    @test K._gmsh_delaunayize_initial!(T, strip, interior) > 0
    @test longest(T, interior) < 2        # now they cross it
    @test violations(T, interior, strip) == 0
    @test count(t -> T.alive[t] && interior[t] && !K._is_ghost_tri(T, t), eachindex(T.alive)) == cells
    @test K._gmsh_delaunayize_initial!(T, strip, interior) == 0
    # Where parameters are lengths the triangulation already is Delaunay: nothing moves.
    flat = K.constrained_delaunay(xs, ys, segments)
    inside = Vector{Bool}(collect(K.classify_interior(flat)))
    before = copy(flat.tv)
    @test K._gmsh_delaunayize_initial!(flat, (u, v) -> (1.0, 0.0, 1.0), inside) == 0
    @test flat.tv == before
    # Four vertex ids key a configuration whatever their order.
    @test K._gmsh_swapquad(7, 3, 9, 1) === (Int32(1), Int32(3), Int32(7), Int32(9))
    @test K._gmsh_swapquad(1, 3, 7, 9) === K._gmsh_swapquad(9, 7, 3, 1)
end
end
