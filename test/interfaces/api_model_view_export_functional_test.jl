module ModelViewFormatFunctionalChecks
using Test, Tessella
const A=Tessella.API

function vtk(path)
    lines=readlines(path);@test first(lines)=="# vtk DataFile Version 2.0"
    @test lines[3]=="ASCII"
    tokens=split(join(lines[4:end],' '));index=Ref(1)
    take()=let value=tokens[index[]];index[]+=1;value;end
    integer()=parse(Int,take())
    @test take()=="DATASET" && take()=="UNSTRUCTURED_GRID"
    points=zeros(3,0);cells=Vector{Int}[];types=Int[];fields=Dict{Tuple{String,String},Any}()
    association="global"
    while index[]<=length(tokens)
        key=take()
        if key=="POINTS"
            n=integer();@test take()=="double"
            points=reshape([parse(Float64,take()) for _ in 1:3n],3,n)
        elseif key=="CELLS"
            n=integer();extent=integer();consumed=0
            for _ in 1:n
                width=integer();push!(cells,[integer() for _ in 1:width]);consumed+=width+1
            end
            @test consumed==extent
        elseif key=="CELL_TYPES"
            n=integer();types=[integer() for _ in 1:n]
        elseif key in ("POINT_DATA","CELL_DATA")
            association=key;n=integer()
            @test n==(key=="POINT_DATA" ? size(points,2) : length(cells))
        elseif key=="FIELD"
            take();n=integer()
            for _ in 1:n
                name=take();components=integer();tuples=integer();type=take()
                T=type=="double" ? Float64 : type=="vtktypeuint64" ? UInt64 : Int
                fields[(association,name)]=reshape([parse(T,take()) for _ in 1:components*tuples],components,tuples)
            end
        else
            error("unsupported VTK token $key")
        end
    end
    return (;points,cells,types,fields)
end

function csv_rows(path)
    function fields(line)
        chars=collect(line);row=String[];io=IOBuffer();quoted=false;i=1
        while i<=length(chars)
            char=chars[i]
            if char=='"'
                if quoted && i<length(chars) && chars[i+1]=='"'
                    print(io,'"');i+=1
                else
                    quoted=!quoted
                end
            elseif char==',' && !quoted
                push!(row,String(take!(io)))
            else
                print(io,char)
            end
            i+=1
        end
        @test !quoted
        push!(row,String(take!(io)));return row
    end
    lines=readlines(path);header=fields(first(lines))
    return [Dict(zip(header,fields(line))) for line in lines[2:end]]
end

function tetrahedron(name)
    A.model.add(name);A.model.add_discrete_entity(3,7)
    base=(UInt64(1)<<53)+UInt64(101);tags=base.+UInt64.(0:9)
    # Gmsh Tet10 edge order: 12,23,31,41,43,42.
    coords=[0 2 0 0 1 1 0 0 0 1;0 0 3 0 0 1.5 1.5 0 1.5 0;0 0 0 4 0 0 0 2 2 2]
    A.mesh.add_nodes(3,7,tags,vec(coords))
    A.mesh.add_nodes(3,7,[base+100],[9.,8.,7.])
    A.mesh.add_elements_by_type(7,11,[301],tags)
    return tags,base+100,coords
end

@testset "model sparse data exports exact UIDs, complete Tet10 topology and times" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        name="p2,\"model\"";tags,orphan,coords=tetrahedron(name)
        view=A.view.add("partial nodal values")
        A.view.add_model_data(view,0,name,"NodeData",[tags[1],tags[5],orphan],[[1.25],[nextfloat(2.5)],[9.]],.25)
        A.view.add_model_data(view,2,name,"NodeData",[tags[10]],[[3.75]],1.75)
        before=(A.view.get_model_data(view,0),A.view.get_model_data(view,2),A.mesh.get_nodes())
        mktempdir() do directory
            path=joinpath(directory,"sparse.vtk");A.view.write(view,path)
            parsed=vtk(path);uids=vec(parsed.fields[("POINT_DATA","node_uid")]);positions=Dict(uid=>i for (i,uid) in enumerate(uids))
            @test Set(uids)==Set(vcat(tags,[orphan]))
            @test parsed.types==[24]
            @test uids[only(parsed.cells).+1]==tags[[1,2,3,4,5,6,7,8,10,9]]
            @test parsed.points[:,positions[tags[9]]]==[0,1.5,2]
            @test parsed.fields[("global","STEP")]==[0 2]
            @test parsed.fields[("global","TIME")]==[.25 1.75]
            @test parsed.fields[("global","ALLOCATED_RANGES")]==reshape([0,2],2,1)
            first=vec(parsed.fields[("POINT_DATA","step_0")]);present=vec(parsed.fields[("POINT_DATA","present_0")])
            @test first[positions[tags[1]]]==1.25 && first[positions[tags[5]]]==nextfloat(2.5) && first[positions[orphan]]==9.
            @test sum(present)==3 && all(isnan(first[i]) for i in eachindex(present) if present[i]==0)
            @test parsed.fields[("POINT_DATA","step_2")][positions[tags[10]]]==3.75
            @test vec(parsed.fields[("CELL_DATA","element_uid")])==UInt64[301]
            @test vec(parsed.fields[("CELL_DATA","element_entity_dim")])==[3] && vec(parsed.fields[("CELL_DATA","element_entity_tag")])==[7]
            path=joinpath(directory,"sparse.csv");A.view.write(view,path);rows=csv_rows(path)
            data=[row for row in rows if row["record"]=="data"]
            @test length(data)==4 && all(row["model"]==name for row in rows)
            @test Set(parse.(UInt64,[row["uid"] for row in data]))==Set([tags[1],tags[5],tags[10],orphan])
            @test only(row["value_1"] for row in data if row["uid"]==string(tags[5]))==string(nextfloat(2.5))
            @test parse.(UInt64,split(only(row["connectivity"] for row in rows if row["record"]=="cell"),';'))==tags
            @test count(row->row["record"]=="node",rows)==11
            @test only(row["step_last"] for row in rows if row["record"]=="range")=="2"
            tensor=A.view.add("tensor per cell")
            A.view.add_model_data(tensor,0,name,"ElementData",[301],[collect(1.:9.)],.5,9)
            path=joinpath(directory,"tensor.vtk");A.view.write(tensor,path);parsed=vtk(path)
            @test parsed.fields[("CELL_DATA","step_0")]==reshape(collect(1.:9.),9,1)
            @test parsed.fields[("CELL_DATA","present_0")]==reshape([1],1,1)
        end
        @test (A.view.get_model_data(view,0),A.view.get_model_data(view,2),A.mesh.get_nodes())==before
    finally
        A.finalize()
    end
end

@testset "ElementNodeData preserves discontinuous values and unknown rows fail atomically" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("triangles");A.model.add_discrete_entity(2,7)
        A.mesh.add_nodes(2,7,[101,102,103,104,105],[0.,0.,0.,2.,0.,0.,0.,3.,0.,2.,3.,0.,9.,8.,7.])
        A.mesh.add_elements_by_type(7,2,[301,302],[101,102,103,102,104,103])
        view=A.view.add("discontinuous")
        A.view.add_model_data(view,0,"triangles","ElementNodeData",[301,302],[[1.,2.,3.],[20.,40.,30.]])
        mktempdir() do directory
            path=joinpath(directory,"element-nodes.vtk");A.view.write(view,path);parsed=vtk(path)
            @test parsed.cells==[[0,1,2],[3,4,5]] && parsed.types==[5,5]
            @test vec(parsed.fields[("POINT_DATA","node_uid")])==UInt64[101,102,103,102,104,103,105]
            @test vec(parsed.fields[("POINT_DATA","step_0")])[1:6]==[1,2,3,20,40,30]
            @test isnan(last(parsed.fields[("POINT_DATA","step_0")])) && last(parsed.fields[("POINT_DATA","present_0")])==0
            path=joinpath(directory,"element-nodes.csv");A.view.write(view,path);rows=csv_rows(path)
            @test count(row->row["record"]=="data",rows)==6
            @test [parse(Float64,row["value_1"]) for row in rows if row["record"]=="data" && row["node_uid"]=="102"]==[2.,20.]
            bad=A.view.add("orphan view UID")
            A.view.add_model_data(bad,0,"triangles","NodeData",[999],[[1.]])
            for extension in (".vtk",".csv")
                path=joinpath(directory,"kept"*extension);write(path,"sentinel")
                @test_throws ArgumentError A.view.write(bad,path)
                @test read(path,String)=="sentinel"
                @test_throws ArgumentError A.view.write(view,path,true)
                @test read(path,String)=="sentinel"
            end
        end
    finally
        A.finalize()
    end
end

@testset "model VTK prevalidates unsupported higher order and held multi-model steps" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("cubic");A.model.add_discrete_entity(2,7)
        _,dimension,_,count,reference,_=A.mesh.get_element_properties(21)
        coords=vcat(reshape(reference,dimension,count),zeros(3-dimension,count))
        tags=collect(101:100+count);A.mesh.add_nodes(2,7,tags,vec(coords))
        A.mesh.add_elements_by_type(7,21,[301],tags)
        cubic=A.view.add("cubic values");A.view.add_model_data(cubic,0,"cubic","NodeData",[101],[[1.]])
        mktempdir() do directory
            path=joinpath(directory,"cubic.vtk");write(path,"sentinel")
            @test_throws ArgumentError A.view.write(cubic,path)
            @test read(path,String)=="sentinel"
            path=joinpath(directory,"cubic.csv");A.view.write(cubic,path);rows=csv_rows(path)
            @test length(split(only(row["connectivity"] for row in rows if row["record"]=="cell"),';'))==count
            multi=A.view.add("two models")
            for (step,name,uid,value) in ((0,"first",201,2.5),(1,"second",202,7.5))
                A.model.add(name);A.model.add_discrete_entity(0,1)
                A.mesh.add_nodes(0,1,[uid],[Float64(step),0.,0.])
                A.view.add_model_data(multi,step,name,"NodeData",[uid],[[value]],.25+step)
            end
            path=joinpath(directory,"multi.vtk");A.view.write(multi,path)
            @test !isfile(path)
            @test vtk(joinpath(directory,"multi_0.vtk")).fields[("POINT_DATA","step_0")]==reshape([2.5],1,1)
            @test vtk(joinpath(directory,"multi_1.vtk")).fields[("POINT_DATA","step_1")]==reshape([7.5],1,1)
            A.view.add_model_data(multi,2,"second","NodeData",[999],[[9.]])
            first=joinpath(directory,"kept_0.vtk");second=joinpath(directory,"kept_1.vtk")
            write(first,"first sentinel");write(second,"second sentinel")
            @test_throws ArgumentError A.view.write(multi,joinpath(directory,"kept.vtk"))
            @test read(first,String)=="first sentinel" && read(second,String)=="second sentinel"
            duplicate=A.view.add("duplicate held names")
            for (step,value) in ((0,1.25),(1,7.5))
                A.model.add("duplicate");A.model.add_discrete_entity(0,1)
                A.mesh.add_nodes(0,1,[701],[0.,0.,0.])
                A.mesh.add_elements_by_type(1,15,[801],[701])
                # Empty model_name binds the actual current slot, even when
                # another held model has the same public name and node UID.
                A.view.add_model_data(duplicate,step,"","NodeData",[701],[[value]])
            end
            path=joinpath(directory,"duplicates.csv");A.view.write(duplicate,path);rows=csv_rows(path)
            data=[row for row in rows if row["record"]=="data"]
            @test [(row["dataset"],row["model"],row["uid"],row["value_1"]) for row in data]==
                [("0","duplicate","701","1.25"),("1","duplicate","701","7.5")]
            for record in ("range","node","cell","data")
                @test Set(row["dataset"] for row in rows if row["record"]==record)==Set(["0","1"])
            end
        end
    finally
        A.finalize()
    end
end
end
