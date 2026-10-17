module VTKLagrangeCellChecks
using Test,Tessella,TOML
const A=Tessella.API
const CASES=TOML.parsefile(joinpath(@__DIR__,"..","fixtures","vtk_lagrange_public_cases.toml"))

function read_vtk(path)
    lines=readlines(path)
    modern=first(lines)=="# vtk DataFile Version 5.1"
    (modern || first(lines)=="# vtk DataFile Version 2.0") && lines[3]=="ASCII" || error("invalid VTK header")
    tokens=split(join(lines[4:end],' '));index=Ref(1)
    take()=let value=tokens[index[]];index[]+=1;value;end
    integer()=parse(Int,take())
    take()=="DATASET" && take()=="UNSTRUCTURED_GRID" || error("invalid VTK dataset")
    points=zeros(3,0);cells=Vector{Int}[];types=Int[]
    fields=Dict{Tuple{String,String},Any}();association="global"
    while index[]<=length(tokens)
        key=take()
        if key=="POINTS"
            count=integer();take()=="double" || error("invalid point type")
            points=reshape([parse(Float64,take()) for _ in 1:3count],3,count)
        elseif key=="CELLS"
            count=integer();extent=integer()
            if modern
                take()=="OFFSETS" && take()=="vtktypeint64" || error("invalid cell offsets")
                offsets=[integer() for _ in 1:count]
                count>=1 && first(offsets)==0 && last(offsets)==extent && issorted(offsets) || error("invalid cell extent")
                take()=="CONNECTIVITY" && take()=="vtktypeint64" || error("invalid cell connectivity")
                nodes=[integer() for _ in 1:extent]
                for cell in 1:count-1;push!(cells,nodes[offsets[cell]+1:offsets[cell+1]]);end
            else
                consumed=0
                for _ in 1:count
                    width=integer();push!(cells,[integer() for _ in 1:width]);consumed+=width+1
                end
                consumed==extent || error("invalid cell extent")
            end
        elseif key=="CELL_TYPES"
            count=integer();types=[integer() for _ in 1:count]
        elseif key in ("POINT_DATA","CELL_DATA")
            association=key;count=integer()
            count==(key=="POINT_DATA" ? size(points,2) : length(cells)) || error("invalid data count")
        elseif key=="FIELD"
            take();count=integer()
            for _ in 1:count
                name=take();components=integer();tuples=integer();type=take()
                T=type=="double" ? Float64 : type=="vtktypeuint64" ? UInt64 : Int
                fields[(association,name)]=reshape([parse(T,take()) for _ in 1:components*tuples],components,tuples)
            end
        else
            error("unexpected VTK token $key")
        end
    end
    (;points,cells,types,fields,modern)
end

world(reference)=[2. 1. .5;.5 3. 1.;1. .5 4.]*reference .+ [7.,-2.,5.]

function fixture(msh,name)
    A.model.add(name)
    _,dimension,_,count,local_coordinates,_=A.mesh.get_element_properties(msh)
    reference=vcat(reshape(local_coordinates,dimension,count),zeros(3-dimension,count))
    family=Tessella.Elements.msh_family(msh)
    if family===:qua || family===:pyr
        reference[1:2,:].=(reference[1:2,:].+1)./2
    elseif family===:hex
        reference.=(reference.+1)./2
    elseif family===:pri
        reference[3,:].=(reference[3,:].+1)./2
    end
    coordinates=world(reference)
    base=(UInt64(1)<<53)+UInt64(10000msh)
    tags=base.+UInt64.(count:-1:1)
    A.model.add_discrete_entity(dimension,7)
    A.mesh.add_nodes(dimension,7,tags,vec(coordinates))
    A.mesh.add_elements_by_type(7,msh,[5000+msh],tags)
    _,actual,_=A.mesh.get_elements(dimension,7);cell=only(only(actual))
    values=[[Float64(node),Float64(node^2),-Float64(node)] for node in 1:count]
    (;name,tags,cell,coordinates,reference,values)
end

function check_export(parsed,model,case)
    slots=only(parsed.cells).+1
    expected=world(reduce(hcat,case["nodes"])./case["order"])
    uids=vec(parsed.fields[("POINT_DATA","node_uid")])
    positions=Dict(uid=>index for (index,uid) in enumerate(model.tags))
    @test parsed.modern && parsed.types==[case["vtk_type"]] && length(slots)==length(model.tags)
    @test isapprox(parsed.points[:,slots],expected;atol=1e-12)
    @test Set(uids)==Set(model.tags) && all(isapprox(parsed.points[:,slot],model.coordinates[:,positions[uids[slot]]];atol=1e-12) for slot in slots)
    @test parsed.fields[("POINT_DATA","step_0")][:,slots]==hcat((model.values[positions[uids[slot]]] for slot in slots)...)
    @test parsed.fields[("POINT_DATA","step_2")][:,slots]==hcat((2 .* model.values[positions[uids[slot]]] .+ [3.,5.,7.] for slot in slots)...)
    @test all(==(1),parsed.fields[("POINT_DATA","present_0")]) && all(==(1),parsed.fields[("POINT_DATA","present_2")])
    @test parsed.fields[("global","STEP")]==[0 2] && parsed.fields[("global","TIME")]==[.25 1.75] &&
          vec(parsed.fields[("CELL_DATA","element_uid")])==[model.cell]
end

function exercise(directory)
    A.initialize()
    try
        A.option("General.Terminal",0)
        for (name,case) in sort!(collect(CASES);by=first)
            println("VTK Lagrange ",name);flush(stdout)
            model=fixture(case["msh"],"lagrange-"*name)
            second=[2 .* values .+ [3.,5.,7.] for values in model.values]
            for kind in ("NodeData","ElementNodeData")
                view=A.view.add(kind)
                tags=kind=="NodeData" ? model.tags : [model.cell]
                first_values=kind=="NodeData" ? model.values : [vcat(model.values...)]
                second_values=kind=="NodeData" ? second : [vcat(second...)]
                A.view.add_model_data(view,0,model.name,kind,tags,first_values,.25,3)
                A.view.add_model_data(view,2,model.name,kind,tags,second_values,1.75,3)
                before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0),A.view.get_model_data(view,2))
                path=joinpath(directory,name*"-"*kind*".vtk")
                A.view.write(view,path)
                check_export(read_vtk(path),model,case)
                @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0),A.view.get_model_data(view,2))==before
            end
        end
    finally
        A.finalize()
    end
end

@testset "complete cubic and catalog maximum VTK cells retain every node and value" begin
    retained=get(ENV,"TESSELLA_VTK_LAGRANGE_RETAIN_DIR",nothing)
    if retained===nothing
        mktempdir(exercise)
    else
        mkpath(retained);exercise(retained)
    end
end

@testset "incomplete higher cells and Pyramid14 reject before output" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            for msh in (20,)
                model=fixture(msh,"unsupported-"*string(msh));view=A.view.add("unsupported")
                A.view.add_model_data(view,0,model.name,"NodeData",model.tags,model.values,0.,3)
                before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0))
                path=joinpath(directory,"unsupported.vtk");write(path,"sentinel")
                @test_throws ArgumentError A.view.write(view,path)
                @test read(path,String)=="sentinel"
                @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0))==before
            end
        end
    finally
        A.finalize()
    end
end
end
