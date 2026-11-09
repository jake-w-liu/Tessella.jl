module CLIMeshAdaptViewFrontendChecks
using Test,Tessella,LinearAlgebra
const C=Tessella.CLI
const G=Tessella.GeoExec
const E=Tessella.Elements

const PATCH="""
General.Terminal=0;Mesh.MeshSizeMin=.45;Mesh.MeshSizeMax=.45;Mesh.Smoothing=0;
Point(1)={0,0,0,.45};Point(2)={1,0,0,.45};
Point(3)={1,1,.2,.45};Point(4)={0,1,0,.45};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Surface(1)={1};Transfinite Curve{1,2,3,4}=4;
"""

function patch_checks(mesh)
    @test Tessella.validate(mesh).ok
    triangles=[block for block in mesh.blocks if E.msh_dimension(block.msh)==2]
    @test !isempty(triangles) && all(block->block.msh==2,triangles)
    @test all(p->abs(p[3]-.2p[1]*p[2])<=1e-10,eachcol(mesh.coords))
    area=0.
    for block in triangles,cell in eachcol(block.nodes)
        p,q,r=(mesh.coords[:,row] for row in cell)
        area+=abs((q[1]-p[1])*(r[2]-p[2])-(q[2]-p[2])*(r[1]-p[1]))/2
    end
    @test area≈1. atol=1e-10
end

function error_text(f)
    try;f();return "";catch error;return sprint(showerror,error);end
end

@testset "CLI and GEO options reach native MeshAdapt on an actual ruled patch" begin
    mktempdir() do dir
        source=joinpath(dir,"patch.geo");output=joinpath(dir,"patch.msh")
        write(source,PATCH*"\nMesh.Algorithm=6;")
        @test C.main([source,"-2","-algo","meshadapt","-o",output])==output
        patch_checks(E.read_mixed_msh(output))
        # The in-script selector reaches the same actual producer without
        # a CLI override. Export remains part of the public frontend path.
        write(source,PATCH*"\nMesh.Algorithm=1;")
        @test C.main([source,"-2","-o",output])==output
        patch_checks(E.read_mixed_msh(output))
    end
end

@testset "GEO view-combination aliases retain real elements and time steps" begin
    mktempdir() do dir
        views=joinpath(dir,"views.pos")
        write(views,"""
View "same" { SP(0,0,0){1}; };
View "same" { SP(0,0,0){2}; };
View "other" { SP(0,0,0){3}; };
""")
        before=read(views)
        source=joinpath(dir,"combine.geo")
        for command in ("TimeSteps","TimeStepsByViewName")
            write(source,"Merge \"views.pos\";Combine $command;")
            execution=G.execute_geo(source)
            @test isempty(execution.errors)
            @test [view.name for view in execution.views]==["other","same_Combine"]
            @test vec(only(last(execution.views).elements).values)==[1.,2.]
        end
        write(source,"Merge \"views.pos\";Combine Views;")
        execution=G.execute_geo(source)
        @test isempty(execution.errors)
        @test length(execution.views)==1 && length(only(execution.views).elements)==3
        write(source,"Merge \"views.pos\";Combine TimeStepsFromAllViews;")
        execution=G.execute_geo(source)
        @test isempty(execution.errors)
        @test vec(only(only(execution.views).elements).values)==[1.,2.,3.]
        @test read(views)==before
        @test C.main([views,"-combine"])==views
        @test occursin("duplicate -combine",error_text(()->C.main([views,"-combine","-combine"])))
        model=joinpath(dir,"hybrid.msh")
        write(model,raw"""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1
11 0 0 0
$EndNodes
$Elements
1
77 15 2 0 1 11
$EndElements
$NodeData
1
"same"
1
0
3
0
1
1
11 7
$EndNodeData
""")
        # The backend's hybrid-data diagnostic proves the flag executes
        # the real combination after import rather than only storing a flag.
        @test occursin("Cannot combine hybrid",error_text(()->C.main([model,"-merge",views,"-combine"])))
    end
end

@testset "CLI combination feeds actual background sizing and preserves its source" begin
    mktempdir() do dir
        views=joinpath(dir,"sizes.pos")
        record(h)="ST(-1,-1,0,3,-1,0,-1,3,0){$h,$h,$h};"
        write(views,"View \"same\" {$(record(.25))};View \"same\" {$(record(.5))};View \"other\" {$(record(.4))};")
        source=joinpath(dir,"square.geo")
        square=replace(PATCH,"{1,1,.2,.45}"=>"{1,1,0,.45}","Surface(1)"=>"Plane Surface(1)",
            "Transfinite Curve{1,2,3,4}=4;"=>"", "Mesh.MeshSizeMin=.45;Mesh.MeshSizeMax=.45;"=>"")
        write(source,square*"\nMesh.Algorithm=1;Merge \"sizes.pos\";Field[1]=PostView;Field[1].ViewIndex=1;Background Field=1;")
        coarse=joinpath(dir,"coarse.msh");combined=joinpath(dir,"combined.msh")
        C.main([source,"-2","-o",coarse])
        @test C.main([source,"-combine","-2","-o",combined])==combined
        count2(mesh)=sum(size(block.nodes,2) for block in mesh.blocks if E.msh_dimension(block.msh)==2;init=0)
        coarse_mesh=E.read_mixed_msh(coarse);combined_mesh=E.read_mixed_msh(combined)
        @test Tessella.validate(combined_mesh).ok
        @test count2(combined_mesh)>count2(coarse_mesh)
        # A BGM file merged after the input can lose its old tag during
        # combination; its actual combined source must remain usable.
        write(source,square*"\nMesh.Algorithm=1;")
        bgm_output=joinpath(dir,"bgm.msh")
        @test C.main([source,"-merge",views,"-combine","-bgm",views,"-2","-o",bgm_output])==bgm_output
        @test Tessella.validate(E.read_mixed_msh(bgm_output)).ok
    end
end

@testset "CLI combination preserves GEO PostView selectors" begin
    mktempdir() do dir
        views=joinpath(dir,"sizes.pos")
        write(views,"View \"same\" { SL(0,0,0,1,0,0){.5,.5}; };" *
            "View \"same\" { SL(0,0,0,1,0,0){.5,.5}; };" *
            "View \"other\" { SL(0,0,0,1,0,0){.9,.9}; };")
        before=read(views)
        source=joinpath(dir,"line.geo");output=joinpath(dir,"line.msh")
        prefix="General.Terminal=0;Mesh.MeshSizeFromPoints=0;" *
            "Mesh.MeshSizeExtendFromBoundary=0;Merge \"sizes.pos\";" *
            "Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(1)={1,2};Field[1]=PostView;"
        for selection in ("ViewTag=1;","ViewIndex=0;","IView=0;",
                "IView=2;Field[1].ViewIndex=0;","ViewIndex=2;Field[1].IView=0;","")
            field=isempty(selection) ? "" : "Field[1]."*selection
            write(source,prefix*field*"Background Field=1;")
            @test C.main([source,"-combine","-1","-o",output])==output
            @test size(E.read_mixed_msh(output).coords,2)==3
        end
        @test read(views)==before
    end
end
end
