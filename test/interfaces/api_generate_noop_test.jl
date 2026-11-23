# The suite includes every test file into Main; this file keeps its aliases to itself.
module APIGenerateNoopTests
using Test
using Tessella
const A = Tessella.API

# generate(3) of a model without Volumes has nothing to mesh. It takes a
# shortcut (`_higher_noop_generate_locked!`) instead of the full plan; both
# must leave the same session: nodes, cells and labels, parametric
# coordinates, maximum labels, element display flags, the returned mesh and
# what a later operation (set_order) builds on it.

function build()
    A.option("General.Terminal", 0)
    geo = A.model.geo
    lc = 0.45
    geo.add_point(0, 0, 0, lc, 1); geo.add_point(1, 0, 0, lc, 2)
    geo.add_point(1, 1, 0, lc, 3); geo.add_point(0, 1, 0, lc, 4)
    geo.add_point(2, 0, 0, lc, 5); geo.add_point(2, 1, 0, lc, 6)
    geo.add_line(1, 2, 1); geo.add_line(2, 3, 2); geo.add_line(3, 4, 3); geo.add_line(4, 1, 4)
    geo.add_line(2, 5, 5); geo.add_line(5, 6, 6); geo.add_line(6, 3, 7)
    geo.add_curve_loop([1, 2, 3, 4], 1); geo.add_plane_surface([1], 1)
    geo.add_curve_loop([5, 6, 7, -2], 2); geo.add_plane_surface([2], 2)
    geo.synchronize()
    A.model.add_physical_group(2, [1]; tag=7)
    A.model.add_physical_group(1, [1, 2]; tag=8)
end

function state()
    nodes = A.mesh.get_nodes(-1, -1, false, true)
    elements = [A.mesh.get_elements(d, t) for (d, t) in A.model.get_entities()]
    tags = reduce(vcat, (reduce(vcat, e[2]; init=UInt64[]) for e in elements); init=UInt64[])
    (nodes, elements, A.mesh.get_max_node_tag(), A.mesh.get_max_element_tag(),
     A.mesh.get_visibility(tags))
end

function run_sequence(shortcut::Bool)
    previous = Tessella.API._HIGHER_NOOP_SHORTCUT[]
    Tessella.API._HIGHER_NOOP_SHORTCUT[] = shortcut
    A.initialize()
    try
        build()
        A.mesh.generate(2)
        _, tags, _ = A.mesh.get_elements(2, 1)
        A.mesh.set_visibility(tags[1][1:2], 0)
        _, points, _ = A.mesh.get_elements(0, 1)
        A.mesh.set_visibility(points[1], 0)
        first = A.mesh.generate(3); s1 = state()
        second = A.mesh.generate(3); s2 = state()
        # Labels in another order than generation's: the plan renumbers them.
        A.mesh.renumber_nodes(); A.mesh.renumber_elements()
        A.mesh.generate(3); s3 = state()
        A.mesh.set_order(2); s4 = state()
        return (first.coords, second.coords, typeof(first), typeof(second), s1, s2, s3, s4)
    finally
        A.finalize()
        Tessella.API._HIGHER_NOOP_SHORTCUT[] = previous
    end
end

@testset "no-op generate(3) equals the full plan" begin
    fast = run_sequence(true)
    slow = run_sequence(false)
    @test fast[1] == slow[1] && fast[2] == slow[2]
    @test fast[3] === slow[3] && fast[4] === slow[4]
    for k in 5:8
        @test fast[k] == slow[k]
    end
end

# Allocation bounds, compiled runs only (the suite runs interpreted). The
# no-op costs about its returned copy of the mesh; a generation of a box
# about a dozen heap objects per cell (13.0 at 0.1.3, measured on a
# 179,000-tetrahedron box).
objects(f) = (GC.gc(); a = Base.gc_num(); f(); d = Base.GC_Diff(Base.gc_num(), a);
    (d.allocd, d.malloc + d.realloc + d.poolalloc + d.bigalloc))
if Base.JLOptions().compile_enabled == 1
    @testset "allocation bounds" begin
        A.initialize()
        try
            build()
            A.mesh.generate(2); A.mesh.generate(3); A.mesh.generate(3)
            bytes, _ = objects(() -> A.mesh.generate(3))
            mesh = A.mesh.generate(3)
            @test bytes <= 4 * Base.summarysize(mesh) + 200_000
        finally
            A.finalize()
        end
        A.initialize()
        try
            A.option("General.Terminal", 0)
            A.model.occ.add_box(0, 0, 0, 1, 1, 1); A.model.occ.synchronize()
            A.option("Mesh.MeshSizeMin", 0.08); A.option("Mesh.MeshSizeMax", 0.08)
            A.mesh.generate(3)
            A.mesh.clear()
            _, count = objects(() -> A.mesh.generate(3))
            cells = sum(length, A.mesh.get_elements()[2])
            @test count / cells <= 16
        finally
            A.finalize()
        end
    end
end

end
