# `.geo` `Plugin(Name).Opt = v` / `Plugin(Name).Run` — upstream executes
# plugins against the model view list, which `_GEO_EXEC_VIEWS` carries here
# (`Merge "*.pos"` order). Options share the API plugin spec; runs swap the
# geo store into `SESSION_VIEWS` under `STATE_LOCK` and restore it.
using Test
using Tessella

function _execute_geo_source(source::AbstractString;files=Dict{String,String}(),mesh_dim=0)
    dir=mktempdir()
    path=joinpath(dir,"script.geo")
    write(path,source)
    for (name,content) in files
        write(joinpath(dir,name),content)
    end
    return execute_geo(path;mesh_dim=mesh_dim)
end

const _TRI_POS="""
View "v" {
ST(0,0,0, 1,0,0, 0,1,0){1.0, 2.0, 3.0};
};
"""

@testset ".geo Plugin(...) dispatch" begin
    # `Run` publishes plugin views onto the execution's own view list —
    # `View = i` addresses them by merge order like upstream `View[i]`.
    ex=_execute_geo_source("""
        Merge "v.pos";
        Plugin(Integrate).View = 0;
        Plugin(Integrate).Run;
        Plugin(MinMax).View = 0;
        Plugin(MinMax).Run;""";
        files=Dict("v.pos"=>_TRI_POS))
    @test length(ex.views)==4
    @test ex.views[1].name=="v"
    @test ex.views[2].name=="v_Integrate"
    @test ex.views[3].name=="v_Min"
    @test ex.views[4].name=="v_Max"
    # ∫v dA over the unit triangle: area 0.5 × mean 2.0 = 1.0.
    integral=ex.views[2].elements[1]
    @test integral.kind==:point
    @test integral.values[1]≈1.0
    @test ex.views[3].elements[1].values[1]≈1.0
    @test ex.views[4].elements[1].values[1]≈3.0

    # `View = -1` (default) resolves the last view, and option expressions
    # evaluate through the `.geo` numeric context.
    ex=_execute_geo_source("""
        Merge "v.pos";
        Plugin(Integrate).View = 1 - 1;
        Plugin(Integrate).Run;
        Plugin(Integrate).View = -1;
        Plugin(Integrate).Run;""";
        files=Dict("v.pos"=>_TRI_POS))
    @test length(ex.views)==3
    # The second Integrate ran on view 1 — the first `v_Integrate` SP record.
    @test ex.views[3].name=="v_Integrate_Integrate"
    @test ex.views[3].elements[1].values[1]≈1.0

    # A view→view plugin: CutPlane slices the triangle on x = 0.5
    # (upstream's plane is `A·x + B·y + C·z + D = 0`, so `D = -0.5`).
    ex=_execute_geo_source("""
        Merge "v.pos";
        Plugin(CutPlane).A = 1;
        Plugin(CutPlane).D = -0.5;
        Plugin(CutPlane).View = 0;
        Plugin(CutPlane).Run;""";
        files=Dict("v.pos"=>_TRI_POS))
    @test length(ex.views)==2
    @test ex.views[2].name=="v_CutPlane"
    @test !isempty(ex.views[2].elements)
    for element in ex.views[2].elements
        for node in axes(element.coords,2)
            @test element.coords[1,node]≈0.5 atol=1e-6
        end
    end

    # `Run = 0` stores the option without executing (upstream's Run button
    # fires only on a nonzero write).
    ex=_execute_geo_source("""
        Merge "v.pos";
        Plugin(Integrate).Run = 0;""";
        files=Dict("v.pos"=>_TRI_POS))
    @test length(ex.views)==1

    # Option sets without `Run` leave the view list alone.
    ex=_execute_geo_source("""
        Merge "v.pos";
        Plugin(Integrate).View = 0;""";
        files=Dict("v.pos"=>_TRI_POS))
    @test length(ex.views)==1

    # Malformed / unknown names and options reject like the API surface.
    @test_throws ArgumentError _execute_geo_source(
        "Plugin(NoSuchPlugin).Run;")
    @test_throws ArgumentError _execute_geo_source(
        "Plugin(Integrate).NoSuchOption = 1;")
    @test_throws ArgumentError _execute_geo_source(
        "Plugin(Integrate).View;")

    # Mesh plugins bind the mid-file mesh as the current slot: SimplePartition
    # slices it and the produced `_MeshPartition` republishes through the
    # `PartitionMesh` route — `Save` emits `$PartitionedEntities`.
    dir=mktempdir()
    geo=joinpath(dir,"part.geo")
    write(geo,"""
        Point(1) = {0,0,0,.5}; Point(2) = {1,0,0,.5};
        Point(3) = {1,1,0,.5}; Point(4) = {0,1,0,.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        Physical Surface("s") = {1};
        Mesh 2;
        Plugin(SimplePartition).NumSlicesX = 2;
        Plugin(SimplePartition).NumSlicesY = 2;
        Plugin(SimplePartition).Run;
        Save "part.msh";
        """)
    ex=execute_geo(geo)
    out=read(joinpath(dir,"part.msh"),String)
    @test occursin("\$PartitionedEntities",out)
    # The session view store is untouched by `.geo` plugin runs.
    @test isempty(Tessella.API.SESSION_VIEWS)
end
