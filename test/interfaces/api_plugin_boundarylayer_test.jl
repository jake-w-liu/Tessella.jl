# `plugin("BoundaryLayer")` — Gmsh 4.15.2 `BoundaryLayer.cpp`: parses the
# space/comma-separated entity selections (non-numeric words and unknown
# entity tags are silently skipped), computes geometric layer widths, then
# mutates the mesh in place (new layer nodes on adjacent curves/faces,
# zero-thickness quads, column splits, Winslow untangling). That in-place
# mesh surgery is not expressible through Tessella's discrete-entity
# records, so the registration keeps the exact upstream option table and
# the selection parsing, then rejects the run explicitly instead of
# pretending to grow a boundary layer.
module PluginBoundaryLayerChecks
using Test,Tessella
const A=Tessella.API

function _square(name)
    A.model.add(name)
    A.model.add_discrete_entity(1,1)
    A.model.add_discrete_entity(2,2)
    A.mesh.add_nodes(2,2,[1,2,3],
        [0.,0.,0.,  1.,0.,0.,  0.,1.,0.])
    A.mesh.add_elements_by_type(2,2,[10],[1,2,3])
end

@testset "Plugin(BoundaryLayer) exact option table" begin
    A.initialize()
    try
        _square("bl")
        # Upstream numeric options (BoundaryLayer.cpp).
        A.plugin.set_number("BoundaryLayer","Thickness",0.05)
        A.plugin.set_number("BoundaryLayer","Size",0.005)
        A.plugin.set_number("BoundaryLayer","Ratio",1.1)
        A.plugin.set_number("BoundaryLayer","Run",0.0)
        # Upstream string options holding entity selections.
        for opt in ("Volumes","Surfaces","Curves","Points")
            A.plugin.set_string("BoundaryLayer",opt,"")
        end
        A.plugin.set_string("BoundaryLayer","Curves","1, 2,3")
        @test_throws ArgumentError A.plugin.set_number(
            "BoundaryLayer","NoSuchOption",1.)
        @test_throws ArgumentError A.plugin.set_string(
            "BoundaryLayer","NoSuchOption","x")
    finally
        A.finalize()
    end
end

@testset "Plugin(BoundaryLayer) rejects the mesh mutation" begin
    A.initialize()
    try
        _square("bl2")
        # Unknown entity tags are silently skipped by upstream's
        # getXByTag resolution — the run still fails with the documented
        # not-implemented error (not a tag error).
        A.plugin.set_string("BoundaryLayer","Curves","1, 999")
        err=try
            A.plugin.run("BoundaryLayer");nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("BoundaryLayer",sprint(showerror,err))
        @test occursin("not implemented",sprint(showerror,err))
    finally
        A.finalize()
    end
end
end
