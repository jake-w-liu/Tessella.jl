# `plugin("CurvedBndDist")` — Gmsh 4.15.2 `CurvedBndDist.cpp`: a bare
# `GMSH_PostPlugin` with no options (only the implicit `Run` action exists —
# the oracle rejects every named option). It scans surface element edge
# closures; for each edge whose first high-order vertex is classified on a
# model curve it emits one `SL` record (values 0, 1) between the farthest
# point of one curve and its closest point on the other — the Hausdorff
# distance between the CAD arc and the finite-element nodal edge. Bounded
# subset: order-2 overlays only (a P2 mid-edge node classified on a curve
# anchors the comparison); lower-order meshes produce upstream's empty view.
# Upstream's adaptive `parametricLine::discretize` recursion does not
# terminate on curved geometry, so the comparison samples both curves at a
# fixed 64-point resolution — the documented semantics, not the literal
# (and defective) upstream max-pair accumulation.
module PluginCurvedBndDistChecks
using Test,Tessella
const A=Tessella.API

function _disk(name;size=0.8)
    A.model.add(name)
    A.model.occ.add_disk(0.,0.,0.,1.,1.)
    A.model.occ.synchronize()
    A.option("Mesh.MeshSizeMin",size)
    A.option("Mesh.MeshSizeMax",size)
    A.mesh.generate(2)
end

@testset "Plugin(CurvedBndDist) option table" begin
    A.initialize()
    try
        _disk("cbd0")
        # Upstream declares no options: only the implicit `Run` action
        # exists and every named option is unknown.
        @test_throws ArgumentError A.plugin.set_number(
            "CurvedBndDist","View",-1.)
        @test_throws ArgumentError A.plugin.set_number(
            "CurvedBndDist","Thickness",0.1)
        @test_throws ArgumentError A.plugin.set_string(
            "CurvedBndDist","Curves","1")
        # `Run` via set_number dispatches like upstream's action column.
        A.plugin.set_number("CurvedBndDist","Run",1.0)
        @test length(A.view.get_tags())==1
    finally
        A.finalize()
    end
end

@testset "Plugin(CurvedBndDist) empty on linear mesh" begin
    A.initialize()
    try
        _disk("cbd1")
        A.plugin.run("CurvedBndDist")
        tags=A.view.get_tags()
        @test length(tags)==1
        # First-order edges have no high-order vertex on a curve: the
        # upstream view is created empty.
        types,data=A.view.get_list_data(only(tags))
        @test isempty(types)
        @test isempty(data)
    finally
        A.finalize()
    end
end

@testset "Plugin(CurvedBndDist) SL records on order-2 disk" begin
    A.initialize()
    try
        _disk("cbd2")
        A.mesh.set_order(2)
        A.plugin.run("CurvedBndDist")
        tags=A.view.get_tags()
        @test length(tags)==1
        types,data=A.view.get_list_data(only(tags))
        @test types==["SL"]
        payload=only(data)
        nrec=length(payload)÷8
        @test nrec>0
        # One record per boundary edge whose mid node sits on the circle.
        for i in 1:nrec
            base=(i-1)*8
            p=(payload[base+1],payload[base+3],payload[base+5])
            q=(payload[base+2],payload[base+4],payload[base+6])
            @test payload[base+7]==0.0
            @test payload[base+8]==1.0
            # The Hausdorff gap along a coarse P2 edge is bounded by the
            # element size — far below it for a smooth circle.
            @test hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3])<0.8
        end
    finally
        A.finalize()
    end
end
end
