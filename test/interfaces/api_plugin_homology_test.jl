# `plugin("HomologyComputation")` — Gmsh 4.15.2 `HomologyComputation.cpp`
# bridged to Tessella's native `compute_homology!` solver. Upstream parses
# the physical-group lists, computes the (relative) (co)homology basis,
# publishes each non-zero chain as a discrete entity plus physical group
# (`H_<dim>{<domain>}<i>` naming), and — `CreatePostProcessingViews != 0` —
# attaches one `"<phys>=<name>"` view per chain. On a single closed curve
# loop H_1 = 1, so the run must publish exactly one 1-chain covering the
# loop; a repeated run recomputes fresh requests instead of accumulating.
module PluginHomologyChecks
using Test,Tessella
const A=Tessella.API

function _loop()
    A.model.add("loop")
    A.model.add_discrete_entity(1,1)
    A.mesh.add_nodes(1,1,[11,12,13,14],
        [0.,0.,0., 1.,0.,0., 1.,1.,0., 0.,1.,0.])
    A.mesh.add_elements_by_type(1,1,[101,102,103,104],
        [11,12,12,13,13,14,14,11])
    A.model.add_physical_group(1,[1];tag=10)
    A.mesh.reclassify_nodes()
end

@testset "Plugin(HomologyComputation)" begin
    A.initialize()
    try
        _loop()
        A.plugin.set_string("HomologyComputation","DomainPhysicalGroups","10")
        A.plugin.set_string("HomologyComputation","DimensionOfChainsToSave","1")
        A.plugin.run("HomologyComputation")

        # One new physical group carrying the H_1 generator.
        groups=A.model.get_physical_groups()
        @test (1,10) in groups
        new=setdiff(groups,[(1,10)])
        @test length(new)==1
        (dim,phys)=only(new)
        @test dim==1
        @test A.model.get_physical_name(dim,phys)=="H_1{10}1"
        ents=A.model.get_entities_for_physical_group(dim,phys)
        @test length(ents)==1

        # The chain entity holds four LINE cells — the whole loop.
        etypes,etags,enodes=A.mesh.get_elements(1,only(ents))
        @test etypes==[1]
        @test length(only(etags))==4
        chain_edges=Set{Tuple{UInt64,UInt64}}()
        flat=only(enodes)
        for i in 1:2:length(flat)
            push!(chain_edges,Tuple(sort!(collect(flat[i:i+1]))))
        end
        @test chain_edges==Set([(11,12),(12,13),(13,14),(11,14)])

        # `CreatePostProcessingViews=1` (default) publishes "<phys>=<name>".
        tags=A.view.get_tags()
        @test length(tags)==1
        dt,data=A.view.get_list_data(only(tags))
        @test dt==["SL"]
        # 4 records × (2 xyz + 2 values); every coefficient expands to a
        # unit-valued cell.
        d=only(data)
        @test length(d)==4*8
        @test all(≈(1.0),d[7:8:end]) && all(≈(1.0),d[8:8:end])

        # A rerun recomputes fresh requests: no duplicated basis groups.
        A.plugin.run("HomologyComputation")
        @test length(setdiff(A.model.get_physical_groups(),[(1,10)]))==2
    finally
        A.finalize()
    end
end

@testset "Plugin(HomologyComputation) option validation" begin
    A.initialize()
    try
        _loop()
        # Upstream `parseStringOpt` errors on unexpected characters at run.
        A.plugin.set_string("HomologyComputation","DomainPhysicalGroups","10,x")
        @test_throws ArgumentError A.plugin.run("HomologyComputation")
        A.plugin.set_string("HomologyComputation","DomainPhysicalGroups","10")
        # Both compute arms off → no-op run, no views.
        A.plugin.set_number("HomologyComputation","ComputeHomology",0.0)
        A.plugin.set_number("HomologyComputation","ComputeCohomology",0.0)
        A.plugin.run("HomologyComputation")
        @test isempty(A.view.get_tags())
        # Cohomology arm publishes H^ chains.
        A.plugin.set_number("HomologyComputation","ComputeCohomology",1.0)
        A.plugin.set_string("HomologyComputation","DomainPhysicalGroups","10")
        A.plugin.set_string("HomologyComputation","DimensionOfChainsToSave","1")
        A.plugin.run("HomologyComputation")
        names=[A.model.get_physical_name(d,p) for (d,p) in
            setdiff(A.model.get_physical_groups(),[(1,10)])]
        @test any(startswith("H^1"),names)
    finally
        A.finalize()
    end
end
end
