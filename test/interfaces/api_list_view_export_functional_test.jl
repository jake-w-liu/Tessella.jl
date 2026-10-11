module ListViewExportFunctionalChecks
using Test, Tessella
const A=Tessella.API

# Parse the public legacy VTK grammar, independently of the writer helpers.
function read_vtk(path)
    lines=readlines(path)
    @test first(lines)=="# vtk DataFile Version 2.0"
    @test lines[3]=="ASCII"
    words=split(join(lines[4:end],' '));position=Ref(1)
    word()=let value=words[position[]];position[]+=1;value;end
    integer()=parse(Int,word())
    numbers(n)=[parse(Float64,word()) for _ in 1:n]
    @test word()=="DATASET" && word()=="UNSTRUCTURED_GRID"
    points=zeros(3,0);cells=Vector{Int}[];types=Int[];arrays=Dict{String,Matrix{Float64}}()
    point_count=0
    while position[]<=length(words)
        key=word()
        if key=="POINTS"
            count=integer();@test word()=="double"
            points=reshape(numbers(3count),3,count)
        elseif key=="CELLS"
            count=integer();extent=integer();consumed=0
            for _ in 1:count
                width=integer();push!(cells,[integer() for _ in 1:width]);consumed+=width+1
            end
            @test consumed==extent
        elseif key=="CELL_TYPES"
            types=[integer() for _ in 1:integer()]
        elseif key=="POINT_DATA"
            point_count=integer();@test point_count==size(points,2)
        elseif key=="FIELD"
            word();count=integer()
            for _ in 1:count
                name=word();components=integer();tuples=integer();@test word()=="double"
                arrays[name]=reshape(numbers(components*tuples),components,tuples)
            end
        elseif key=="SCALARS"
            name=word();@test word()=="double";components=integer()
            @test word()=="LOOKUP_TABLE" && word()=="default"
            arrays[name]=reshape(numbers(components*point_count),components,point_count)
        elseif key=="VECTORS"
            name=word();@test word()=="double"
            arrays[name]=reshape(numbers(3point_count),3,point_count)
        else
            error("unexpected VTK token $key")
        end
    end
    return (;points,cells,types,arrays)
end

raw_rows(path,separator)=map(readlines(path)) do line
    parse.(Float64,separator==' ' ? split(line) : split(line,separator))
end

@testset "list TXT, CSV and VTK retain coordinates, independent corners and steps" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        view=A.view.add("scalar list")
        point=[0.,0.,1.,nextfloat(7.),17.]
        triangle=[0.,2.,0.,0.,0.,3.,1.,1.,1.,1.,3.,4.,11.,13.,14.]
        A.view.add_list_data(view,"SP",1,point)
        A.view.add_list_data(view,"ST",1,triangle)
        before=A.view.get_list_data(view)
        point[end]=999.;triangle[1]=999.
        mktempdir() do directory
            txt=joinpath(directory,"values.txt");csv=joinpath(directory,"values.csv")
            vtk=joinpath(directory,"values.vtk")
            for path in (txt,csv,vtk);A.view.write(view,path);end
            rows=raw_rows(txt,' ')
            @test length(rows)==4
            @test rows[1]==[0,0,0,0,0,0,1,nextfloat(7.)]
            @test rows[2]==[0,0,0,1,0,0,1,1,0,0,0,1,2,0,1,3,0,0,0,1,0,3,1,4]
            @test rows[3]==[1,1,0,0,0,0,1,17]
            @test rows[4][8:8:end]==[11,13,14]
            @test raw_rows(csv,',')==rows
            @test occursin(',',read(csv,String))
            parsed=read_vtk(vtk)
            @test parsed.points==[0 0 2 0;0 0 0 3;1 1 1 1]
            @test parsed.cells==[[0],[1,2,3]] && parsed.types==[1,5]
            @test parsed.arrays["TIME"]==[0. 1.]
            @test parsed.arrays["step_0"]==reshape([nextfloat(7.),1.,3.,4.],1,4)
            @test parsed.arrays["step_1"]==[17. 11. 13. 14.]
            for path in (txt,csv)
                bytes=read(path);A.view.write(view,path,true)
                @test read(path)==vcat(bytes,bytes)
            end
            bytes=read(vtk)
            @test_throws ArgumentError A.view.write(view,vtk,true)
            @test read(vtk)==bytes
            # The do-block closes its stream: immediate rename/read succeeds.
            moved=joinpath(directory,"moved.vtk");mv(vtk,moved)
            @test read(moved)==bytes
        end
        @test A.view.get_list_data(view)==before
    finally
        A.finalize()
    end
end

@testset "list vector, nonsymmetric tensor and volume cell exports" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        vector=A.view.add("vector list")
        A.view.add_list_data(vector,"VL",1,
            [2.,4.,3.,5.,7.,9.,1.,2.,3.,4.,5.,6.,11.,12.,13.,14.,15.,16.])
        tensor=A.view.add("tensor \"list\"")
        A.view.add_list_data(tensor,"TP",1,vcat([2.,3.,4.],collect(1.:18.)))
        volume=A.view.add("volume shapes")
        shapes=[("SS",[0 1 0 0;0 0 1 0;0 0 0 1]),
            ("SH",[0 1 1 0 0 1 1 0;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1]),
            ("SI",[0 1 0 0 1 0;0 0 1 0 0 1;0 0 0 1 1 1]),
            ("SY",[0 1 1 0 .5;0 0 1 1 .5;0 0 0 0 1])]
        for (kind,coords) in shapes
            A.view.add_list_data(volume,kind,1,vcat(vec(permutedims(coords)),collect(1.:size(coords,2))))
        end
        mktempdir() do directory
            for (tag,name,components,expected,second) in (
                (vector,"vector",3,reshape(collect(1.:6.),3,2),reshape(collect(11.:16.),3,2)),
                (tensor,"tensor",9,reshape(collect(1.:9.),9,1),reshape(collect(10.:18.),9,1)))
                path=joinpath(directory,name*".vtk");A.view.write(tag,path)
                parsed=read_vtk(path)
                @test parsed.arrays["step_0"]==expected
                @test parsed.arrays["step_1"]==second
                @test size(parsed.arrays["step_0"],1)==components
                csv=joinpath(directory,name*".csv");txt=joinpath(directory,name*".txt")
                A.view.write(tag,csv);A.view.write(tag,txt)
                @test raw_rows(csv,',')==raw_rows(txt,' ')
                @test length(first(raw_rows(csv,',')))==(7+components)*size(expected,2)
            end
            path=joinpath(directory,"volume.vtk");A.view.write(volume,path)
            parsed=read_vtk(path)
            @test parsed.types==[10,12,13,14]
            @test length.(parsed.cells)==[4,8,6,5]
            @test parsed.points==hcat((coords for (_,coords) in shapes)...)
        end
    finally
        A.finalize()
    end
end

@testset "list export validates complete supported payloads before opening" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        malformed=A.view.add("incomplete")
        A.view.add_list_data(malformed,"ST",1,zeros(10))
        unsupported=A.view.add("unsupported cell")
        A.view.add_list_data(unsupported,"ST2",1,zeros(12))
        mismatch=A.view.add("different time counts")
        A.view.add_list_data(mismatch,"SP",1,[0.,0.,0.,1.,2.])
        A.view.add_list_data(mismatch,"SL",1,[0.,1.,0.,0.,0.,0.,1.,2.])
        mixed=A.view.add("mixed component counts")
        A.view.add_list_data(mixed,"SP",1,[0.,0.,0.,1.])
        A.view.add_list_data(mixed,"VP",1,[1.,2.,3.,4.,5.,6.])
        A.model.add("node model");A.model.add_discrete_entity(0,1)
        A.mesh.add_nodes(0,1,[101],[0.,0.,0.])
        mixed_storage=A.view.add("mixed storage")
        A.view.add_list_data(mixed_storage,"SP",1,[0.,0.,0.,1.])
        A.view.add_model_data(mixed_storage,0,"node model","NodeData",[101],[[2.]])
        mktempdir() do directory
            for extension in (".txt",".csv",".vtk"),tag in (malformed,unsupported,mismatch,mixed_storage)
                path=joinpath(directory,"preserved"*extension);write(path,"sentinel")
                @test_throws ArgumentError A.view.write(tag,path)
                @test read(path,String)=="sentinel"
            end
            path=joinpath(directory,"mixed.vtk");write(path,"sentinel")
            @test_throws ArgumentError A.view.write(mixed,path)
            @test read(path,String)=="sentinel"
            # Mixed components have an unambiguous raw text/CSV representation.
            for extension in (".txt",".csv")
                path=joinpath(directory,"mixed"*extension);A.view.write(mixed,path)
                @test length(readlines(path))==2
            end
        end
    finally
        A.finalize()
    end
end
end
