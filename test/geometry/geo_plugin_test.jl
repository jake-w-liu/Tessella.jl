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
    @test ex.views[2].name=="v_Levelset"
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

    # Unknown plugin names/options are recoverable `Error   :` diagnostics
    # upstream — execution continues and the run records them.
    ex=_execute_geo_source("Plugin(NoSuchPlugin).Run;")
    @test ex.errors==["Failed action 'Run' or unknown plugin 'NoSuchPlugin'"]
    ex=_execute_geo_source("Plugin(Integrate).NoSuchOption = 1;")
    @test ex.errors==["Unknown option 'NoSuchOption' or plugin 'Integrate'"]
    ex=_execute_geo_source("Plugin(Integrate).View;")
    @test ex.errors==["Failed action 'View' or unknown plugin 'Integrate'"]

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

@testset ".geo Save/Print/Delete View" begin
    dir=mktempdir()
    write(joinpath(dir,"in.pos"),"""
        View "v0" {
        SP(0,0,0){1.0};
        SP(1,0,0){2.0};
        };
        View "v1" {
        };
        """)
    geo=joinpath(dir,"t.geo")
    write(geo,"""
        Merge "in.pos";
        Save View[0] "out0.pos";
        Save View[1] "out1.pos";
        Delete View[0];
        """)
    ex=execute_geo(geo)
    @test length(ex.views)==1
    # `Save View` exports `.pos` records from the exec view list — including
    # an empty view's name-only table (upstream writes the same block).
    saved=read(joinpath(dir,"out0.pos"),String)
    @test occursin("View \"v0\"",saved)
    @test occursin("SP(0,0,0){1};",saved)
    @test occursin("SP(1,0,0){2};",saved)
    printed=read(joinpath(dir,"out1.pos"),String)
    @test printed=="View \"v1\" {\n};\n" || occursin("View \"v1\"",printed)
    # `Delete View[0]` left only the (now index-0) empty view.
    @test ex.views[1].name=="v1"
    @test isempty(ex.views[1].elements)
    # Out-of-range and non-verb view operations are recoverable upstream
    # `yymsg(0)` diagnostics, not rejections.
    ex=_execute_geo_source(
        "Merge \"in.pos\";\nSave View[9] \"x.pos\";",
        files=Dict("in.pos"=>"View \"v\" {\n};\n"))
    @test ex.errors==["Unknown view 9"]
    ex=_execute_geo_source("Save View[0] \"x.pos\";")
    @test ex.errors==["Unknown view 0"]
    ex=_execute_geo_source(
        "Merge \"in.pos\";\nDelete View[5];",
        files=Dict("in.pos"=>"View \"v\" {\n};\n"))
    @test ex.errors==["Unknown view 5"]
    # A non-`.pos` extension takes upstream's flat-text view export.
    write(joinpath(dir,"in2.pos"),"View \"v\" {\nSP(0,0,0){1.0};\n};\n")
    geo=joinpath(dir,"vt.geo")
    write(geo,"Merge \"in2.pos\";\nSave View[0] \"x.txt\";")
    ex=execute_geo(geo)
    @test isempty(ex.errors)
    @test read(joinpath(dir,"x.txt"),String)=="0 0 0 0 0 0 0 1 \n" ||
          read(joinpath(dir,"x.txt"),String)=="0 0 0 0 0 0 0 1 \r\n"
    # Upstream: `Print` is not a view verb — "Unknown operation 'Print' on
    # view N" — but only after the index resolves to a real view.
    ex=_execute_geo_source(
        "Merge \"in.pos\";\nPrint View[0] \"x.pos\";",
        files=Dict("in.pos"=>"View \"v\" {\n};\n"))
    @test ex.errors==["Unknown operation 'Print' on view 0"]
    ex=_execute_geo_source("Print View[9] \"x.pos\";")
    @test ex.errors==["Unknown view 9"]
    # `Delete Empty Views` drops only the record-less entries.
    geo2=joinpath(dir,"t2.geo")
    write(geo2,"""
        Merge "in.pos";
        Delete Empty Views;
        """)
    ex2=execute_geo(geo2)
    @test length(ex2.views)==1
    @test ex2.views[1].name=="v0"
end

@testset ".geo Alias/Combine views" begin
    dir=mktempdir()
    nl=Sys.iswindows() ? "\r\n" : "\n"
    # Two views named "a" plus singleton "b" — mirrors the upstream probe:
    # by-name combine keeps the singleton first, then "a_Combine".
    write(joinpath(dir,"in.pos"),"""
        View "a" {
        SP(0,0,0){1.0};
        };
        View "b" {
        SP(1,1,1){2.0};
        };
        View "a" {
        SP(2,2,2){3.0};
        };
        """)

    # `Alias View[i]`/`AliasWithOptions View[i]` append an independent copy.
    geo=joinpath(dir,"alias.geo")
    write(geo,"""
        Merge "in.pos";
        Alias View[0];
        AliasWithOptions View[0];
        Save View[4] "alias.pos";
        """)
    ex=execute_geo(geo)
    @test length(ex.views)==5
    @test ex.views[4].name=="a"
    @test length(ex.views[4].elements)==1
    @test read(joinpath(dir,"alias.pos"),String)==
        "View \"a\" {"*nl*"SP(0,0,0){1};"*nl*"};"*nl
    ex=_execute_geo_source("Alias View[0];")
    @test ex.errors==["Unknown view 0"]
    ex=_execute_geo_source(
        "Merge \"in.pos\";\nAlias View[9];",
        files=Dict("in.pos"=>"View \"v\" {\n};\n"))
    @test ex.errors==["Unknown view 9"]

    # `Combine ElementsByViewName`: singletons stay, multi-view names merge.
    geo=joinpath(dir,"byname.geo")
    write(geo,"""
        Merge "in.pos";
        Combine ElementsByViewName;
        Save View[0] "s0.pos";
        Save View[1] "s1.pos";
        """)
    ex=execute_geo(geo)
    @test [v.name for v in ex.views]==["b","a_Combine"]
    @test read(joinpath(dir,"s1.pos"),String)==
        "View \"a_Combine\" {"*nl*"SP(0,0,0){1};"*nl*"SP(2,2,2){3};"*nl*"};"*nl

    # `Combine ElementsFromAllViews` collapses everything into "all_Combine".
    geo=joinpath(dir,"all.geo")
    write(geo,"""
        Merge "in.pos";
        Combine ElementsFromAllViews;
        Save View[0] "all.pos";
        """)
    ex=execute_geo(geo)
    @test length(ex.views)==1
    @test read(joinpath(dir,"all.pos"),String)==
        "View \"all_Combine\" {"*nl*"SP(0,0,0){1};"*nl*"SP(1,1,1){2};"*nl*
        "SP(2,2,2){3};"*nl*"};"*nl

    # `Combine TimeStepsByViewName` zips aligned rows into step columns —
    # byte-identical to the upstream `a_Combine` table.
    geo=joinpath(dir,"ts.geo")
    write(geo,"""
        Merge "in.pos";
        Combine TimeStepsByViewName;
        Save View[1] "ts.pos";
        """)
    ex=execute_geo(geo)
    @test [v.name for v in ex.views]==["b","a_Combine"]
    @test read(joinpath(dir,"ts.pos"),String)==
        "View \"a_Combine\" {"*nl*"TIME{0,1};"*nl*"SP(0,0,0){1,3};"*nl*"};"*nl

    ex=_execute_geo_source(
        "Merge \"in.pos\";\nCombine BogusMode;",
        files=Dict("in.pos"=>"View \"v\" {\n};\n"))
    @test ex.errors==["Unknown 'Combine' command"]
    # Combining zero views is a no-op upstream.
    ex0=_execute_geo_source("Combine ElementsFromAllViews;")
    @test isempty(ex0.views)
end
