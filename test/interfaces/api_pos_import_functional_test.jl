module PosImportFunctionalChecks
using Test,Tessella
const A=Tessella.API
const P=Tessella.Post
const E=Tessella.Elements

snapshot()=(models=A.model.list(),name=A.model.get_current(),file=A.model.get_file_name(),
    entities=A.model.get_entities(),nodes=A.mesh.get_nodes(),cells=A.mesh.get_elements(),
    views=A.view.get_tags(),factor=A.option("Mesh.MeshSizeFactor"),log=A.logger.get())
vtk_time(path)=let lines=readlines(path),row=findfirst(line->startswith(line,"TIME 1 "),lines)
    parse.(Float64,split(lines[row+1]))
end

function exercise(directory)
    path=joinpath(directory,"lists.pos")
    Base.write(path,"View \"scalar time\" {TIME{0.25,1.5};SL(0,0,0,1,0,0){1,2,3,4};SP(2,0,0){9,10};};\n"*
        "View \"vector\" {VT(0,0,0,1,0,0,0,1,0){1,0,0,0,1,0,0,0,1};};\nView \"empty\" {};\n")
    empty_path=joinpath(directory,"empty.pos");Base.write(empty_path,"// no views\n")
    @testset "actual list open and literal multistep times" begin
        A.initialize()
        try
            A.option("General.Terminal",0);prior=A.view.add("prior")
            @test A.open(path)===nothing && A.model.list()==["lists"] &&
                A.model.get_file_name()==path && isempty(A.model.get_entities())
            tags=setdiff(A.view.get_tags(),[prior]);scalar,vector,empty=tags
            @test [A.view.get_string(tag,"Name") for tag in tags]==["scalar time","vector","empty"] &&
                [A.view.get_number(tag,"NbTimeStep") for tag in tags]==[2,1,0]
            expected=(String["SP","SL"],[[2.,0,0,9,10],[0.,1,0,0,0,0,1,2,3,4]])
            @test A.view.get_list_data(scalar)==expected &&
                A.view.get_list_data(vector)==(["VT"],[[0.,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]])
            detached=A.view.get_list_data(scalar);detached[2][1][1]=99
            @test A.view.get_list_data(scalar)==expected
            A.view.set_string(scalar,"Name","renamed");A.view.set_string(scalar,"custom","retained")
            @test A.view.get_string(scalar,"Name")=="renamed" && A.view.get_string(scalar,"custom")=="retained"
            pos=joinpath(directory,"roundtrip.pos");csv=joinpath(directory,"roundtrip.csv")
            vtk=joinpath(directory,"roundtrip.vtk")
            for output in (pos,csv,vtk);A.view.write(scalar,output);end
            parsed=only(P.read_pos(pos))
            @test parsed.name=="renamed" && parsed.time==[.25,1.5] &&
                parsed.elements[2].values==reshape([1.,2,3,4],1,2,2)
            @test [parse(Float64,split(line,',')[2]) for line in readlines(csv)]==[.25,.25,1.5,1.5]
            @test vtk_time(vtk)==[.25,1.5]
            # Existing fields select a step index; physical TIME selection is separate.
            @test Tessella.SizeField.field_value(P.postview_field(parsed),.5,0,0)≈1.5
        finally;A.finalize();end
    end
    @testset "single and empty times retain declared metadata" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            timed=joinpath(directory,"timed.pos")
            Base.write(timed,"View \"single\" {TIME{0.25};SP(0,0,0){5};};\n"*
                "View \"plain empty\" {};\nView \"timed empty\" {TIME{0};};\n")
            A.model.add("named empty");names=A.model.list();names[end]="timed"
            A.open(timed);tags=A.view.get_tags()
            @test A.model.list()==names &&
                [A.view.get_number(tag,"NbTimeStep") for tag in tags]==[1,0,0]
            times=Vector{Float64}[]
            for (ordinal,tag) in enumerate(tags)
                output=joinpath(directory,"time_$ordinal.pos");A.view.write(tag,output)
                push!(times,only(P.read_pos(output)).time)
            end
            @test times==[[.25],Float64[],[0.]]
            old=A.view.add("old implicit");A.view.add_list_data(old,"SP",1,[0.,0,0,7])
            output=joinpath(directory,"old-implicit.pos");A.view.write(old,output)
            @test !occursin("TIME",read(output,String)) && only(P.read_pos(output)).time==[0.]
        finally;A.finalize();end
    end
    @testset "held model views and views-only merge" begin
        A.initialize()
        try
            A.option("General.Terminal",0);A.model.add("kept");A.model.set_file_name("kept.msh")
            A.model.add_discrete_entity(0,17);A.mesh.add_nodes(0,17,[501],[2.,3,4])
            A.mesh.add_elements_by_type(17,15,[701],[501]);A.mesh.rebuild_node_cache(false)
            held=A.view.add("held");A.view.add_model_data(held,0,"kept","NodeData",[501],[[7.25]])
            models=A.model.list();A.open(path)
            @test A.model.list()==vcat(models,["lists"]) && A.model.get_current()=="lists" &&
                isempty(first(A.mesh.get_nodes())) && A.view.get_model_data(held,0)[2]==UInt64[501]
            A.option("Mesh.Binary",1);A.option("Mesh.MshFileVersion",4.1)
            output=joinpath(directory,"held-binary.msh");A.view.write(held,output)
            mesh=E.read_mixed_msh(output)
            @test mesh.entity_data.external_node_tags==UInt64[501] &&
                mesh.entity_data.node_entities==[(0,Int32(17))] && only(mesh.data_sections).values==[7.25]
            A.model.set_current("kept");before=snapshot();A.merge(path);after=snapshot()
            @test Base.merge(after,(views=before.views,))==before && length(after.views)==length(before.views)+3
            before=snapshot();A.merge(empty_path)
            @test snapshot()==before
            models=A.model.list();views=A.view.get_tags();A.open(empty_path)
            @test A.model.list()==vcat(models,["empty"]) && A.view.get_tags()==views && A.model.get_file_name()==empty_path
        finally;A.finalize();end
    end
    @testset "all 24 witnessed linear families" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.merge(joinpath(@__DIR__,"../data/session_pos_import/all-linear-types.pos"))
            tag=only(A.view.get_tags());types,data=A.view.get_list_data(tag)
            @test types==split("SP VP TP SL VL TL ST VT TT SQ VQ TQ SS VS TS SH VH TH SI VI TI SY VY TY")
            coords=(reshape([0.,0,0],3,1),[0. 1;0 0;0 0],[0. 1 0;0 0 1;0 0 0],
                [0. 1 1 0;0 0 1 1;0 0 0 0],[0. 1 0 0;0 0 1 0;0 0 0 1],
                [0. 1 1 0 0 1 1 0;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1],
                [0. 1 0 0 1 0;0 0 1 0 0 1;0 0 0 1 1 1],
                [0. 1 1 0 .5;0 0 1 1 .5;0 0 0 0 1])
            # The point's literal XYZ are all zero; the remaining corners are distinct source coordinates.
            flat=Vector{Float64}[]
            for (shape,corners) in enumerate(coords),components in (1,3,9)
                xyz=shape==1 ? [0.,0,0] : vec(permutedims(corners))
                push!(flat,vcat(xyz,collect(1.:components*size(corners,2))))
            end
            @test data==flat
        finally;A.finalize();end
    end
    @testset "whole-file rejection preserves session and tag clock" begin
        A.initialize()
        try
            A.option("General.Terminal",0);A.logger.start()
            A.model.add("retained");A.model.add_discrete_entity(0,9)
            prior=A.view.add("prior");A.view.add_list_data(prior,"SP",1,[0.,0,0,7])
            bodies=("View \"bad\" {SL(0,0){1};};",
                "View \"text\" {T2(1,2,10){\"a\"};};",
                "View \"high\" {SL2(0,0,0,1,0,0,.5,0,0){1,2,3};};",
                "View \"scheme\" {INTERPOLATION_SCHEME{{1}}{{0}};};")
            for (ordinal,body) in enumerate(bodies)
                invalid=joinpath(directory,"invalid_$ordinal.pos")
                Base.write(invalid,"View \"valid prefix\" {SP(0,0,0){1};};\n"*body)
                before=snapshot()
                @test_throws ArgumentError (isodd(ordinal) ? A.open(invalid) : A.merge(invalid))
                @test snapshot()==before && A.view.get_list_data(prior)==(["SP"],[[0.,0,0,7]])
            end
            @test A.view.add("clock witness")==prior+1
        finally;A.finalize();end
    end
end

directory=get(ENV,"TESSELLA_POS_IMPORT_DIR","")
if isempty(directory);mktempdir(exercise);else;mkpath(directory);exercise(directory);end
end
