module VTKQuadraticFamilyChecks
using Test, Tessella
const A=Tessella.API

function read_vtk(path)
    lines=readlines(path)
    @test first(lines)=="# vtk DataFile Version 2.0" && lines[3]=="ASCII"
    tokens=split(join(lines[4:end],' '));index=Ref(1)
    take()=let value=tokens[index[]];index[]+=1;value;end
    integer()=parse(Int,take())
    @test take()=="DATASET" && take()=="UNSTRUCTURED_GRID"
    points=zeros(3,0);cells=Vector{Int}[];types=Int[]
    fields=Dict{Tuple{String,String},Any}();association="global"
    while index[]<=length(tokens)
        key=take()
        if key=="POINTS"
            count=integer();@test take()=="double"
            points=reshape([parse(Float64,take()) for _ in 1:3count],3,count)
        elseif key=="CELLS"
            count=integer();extent=integer();consumed=0
            for _ in 1:count
                width=integer();push!(cells,[integer() for _ in 1:width]);consumed+=width+1
            end
            @test consumed==extent
        elseif key=="CELL_TYPES"
            count=integer();types=[integer() for _ in 1:count]
        elseif key in ("POINT_DATA","CELL_DATA")
            association=key;count=integer()
            @test count==(key=="POINT_DATA" ? size(points,2) : length(cells))
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
    return (;points,cells,types,fields)
end

# Independent topology expectations: corner vertices, edge midpoints, face
# centers and the volume center. No production permutation is consulted.
function topology(msh)
    if msh==10
        vertices=[0. 1. 1. 0.;0. 0. 1. 1.;0. 0. 0. 0.]
        edges=((1,2),(2,3),(3,4),(4,1));faces=((1,2,3,4),)
        type=28;interior=false
    elseif msh==12
        vertices=[0. 1. 1. 0. 0. 1. 1. 0.;0. 0. 1. 1. 0. 0. 1. 1.;0. 0. 0. 0. 1. 1. 1. 1.]
        edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
            (1,5),(2,6),(3,7),(4,8))
        faces=((1,4,8,5),(2,3,7,6),(1,2,6,5),(4,3,7,8),(1,2,3,4),(5,6,7,8))
        type=29;interior=true
    else
        vertices=[0. 1. 0. 0. 1. 0.;0. 0. 1. 0. 0. 1.;0. 0. 0. 1. 1. 1.]
        edges=((1,2),(2,3),(3,1),(4,5),(5,6),(6,4),(1,4),(2,5),(3,6))
        faces=((1,2,5,4),(2,3,6,5),(3,1,4,6))
        type=32;interior=false
    end
    center(indices)=vec(sum(vertices[:,collect(indices)];dims=2))/length(indices)
    reference=hcat(vertices,(center(edge) for edge in edges)...,(center(face) for face in faces)...)
    interior && (reference=hcat(reference,center(axes(vertices,2))))
    return type,reference
end

world(reference)=[2. 1. .5;.5 3. 1.;1. .5 4.]*reference .+ [7.,-2.,5.]
value(point)=[1+point[1]+2point[2]+3point[3],point[1]*point[2]-point[3],point[3]^2+2point[1]]

function fixture(msh)
    name="quadratic-$msh";A.model.add(name)
    _,dimension,_,count,local_coordinates,_=A.mesh.get_element_properties(msh)
    reference=vcat(reshape(local_coordinates,dimension,count),zeros(3-dimension,count))
    if msh==13
        reference[3,:].=(reference[3,:].+1)./2
    else
        reference[1:dimension,:].=(reference[1:dimension,:].+1)./2
    end
    coordinates=world(reference)
    # Descending, exactly retained labels above Float64's integer precision.
    base=(UInt64(1)<<53)+UInt64(100msh)
    tags=base.+UInt64.(count:-1:1)
    requested=base+(UInt64(1)<<31)+UInt64(1000)
    A.model.add_discrete_entity(dimension,7)
    A.mesh.add_nodes(dimension,7,tags,vec(coordinates))
    A.mesh.add_elements_by_type(7,msh,[requested],tags)
    _,element_tags,_=A.mesh.get_elements(dimension,7);cell=only(only(element_tags))
    # Opaque Gmsh 4.15.2 observations: element constructors normalize requested
    # tags through signed Int32. Bind view rows to the actual public identity.
    native_actual=Dict(10=>UInt64(18446744071562069968),
        12=>UInt64(18446744071562070168),13=>UInt64(18446744071562070268))
    @test cell==native_actual[msh] && cell>(UInt64(1)<<53) && cell!=requested
    type,vtk_reference=topology(msh);expected=world(vtk_reference)
    lookup=Dict(Tuple(coordinates[:,node])=>tags[node] for node in axes(coordinates,2))
    expected_tags=[lookup[Tuple(point)] for point in eachcol(expected)]
    values=[value(point) for point in eachcol(coordinates)]
    return (;name,tags,cell,coordinates,type,expected,expected_tags,values)
end

@testset "complete P2 Quad9, Hex27 and Prism18 preserve topology, UIDs and every vector step" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            for msh in (10,12,13)
                model=fixture(msh);view=A.view.add("node vectors")
                A.view.add_model_data(view,0,model.name,"NodeData",model.tags,model.values,.25,3)
                second=[2 .* row .+ [3.,5.,7.] for row in model.values]
                A.view.add_model_data(view,2,model.name,"NodeData",model.tags,second,1.75,3)
                before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0),A.view.get_model_data(view,2))
                path=joinpath(directory,"nodes-$msh.vtk");A.view.write(view,path);parsed=read_vtk(path)
                slots=only(parsed.cells).+1;uids=vec(parsed.fields[("POINT_DATA","node_uid")])
                @test parsed.types==[model.type] && length(slots)==length(model.tags)
                @test parsed.points[:,slots]==model.expected
                @test uids[slots]==model.expected_tags && Set(uids)==Set(model.tags)
                @test parsed.fields[("POINT_DATA","step_0")][:,slots]==hcat((value(point) for point in eachcol(model.expected))...)
                @test parsed.fields[("POINT_DATA","step_2")][:,slots]==hcat((2 .* value(point) .+ [3.,5.,7.] for point in eachcol(model.expected))...)
                @test all(==(1),parsed.fields[("POINT_DATA","present_0")]) && all(==(1),parsed.fields[("POINT_DATA","present_2")])
                @test parsed.fields[("global","STEP")]==[0 2] && parsed.fields[("global","TIME")]==[.25 1.75]
                @test vec(parsed.fields[("CELL_DATA","element_uid")])==[model.cell]
                @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0),A.view.get_model_data(view,2))==before
                nodal=A.view.add("element node vectors")
                A.view.add_model_data(nodal,0,model.name,"ElementNodeData",[model.cell],[vcat(model.values...)],.5,3)
                A.view.add_model_data(nodal,1,model.name,"ElementNodeData",[model.cell],[vcat(second...)],2.5,3)
                @test A.view.get_model_data(nodal,0)==(["ElementNodeData"],[model.cell],vcat(model.values...))
                path=joinpath(directory,"element-nodes-$msh.vtk");A.view.write(nodal,path);parsed=read_vtk(path)
                slots=only(parsed.cells).+1
                @test parsed.types==[model.type] && parsed.points[:,slots]==model.expected
                @test vec(parsed.fields[("POINT_DATA","node_uid")])[slots]==model.expected_tags
                @test parsed.fields[("POINT_DATA","step_0")][:,slots]==hcat((value(point) for point in eachcol(model.expected))...)
                @test parsed.fields[("POINT_DATA","step_1")][:,slots]==hcat((2 .* value(point) .+ [3.,5.,7.] for point in eachcol(model.expected))...)
                @test parsed.fields[("global","TIME")]==[.5 2.5]
            end
        end
    finally
        A.finalize()
    end
end

@testset "view row UIDs accept full UInt64 and reject invalid tags atomically" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("wide point");A.model.add_discrete_entity(0,7)
        uid=(UInt64(1)<<63)+UInt64(123)
        A.mesh.add_nodes(0,7,[uid],[7.,-2.,5.]);A.mesh.add_elements_by_type(7,15,[301],[uid])
        view=A.view.add("wide point value")
        A.view.add_model_data(view,0,"wide point","NodeData",[uid],[[1.25]])
        before=A.view.get_model_data(view,0)
        @test before==(["NodeData"],[uid],[1.25])
        for invalid in (true,-1,big(1)<<64,1.0)
            @test_throws ArgumentError A.view.add_model_data(view,0,"wide point","NodeData",[invalid],[[9.]])
        end
        @test A.view.get_model_data(view,0)==before
        mktempdir() do directory
            path=joinpath(directory,"wide.vtk");A.view.write(view,path);parsed=read_vtk(path)
            @test vec(parsed.fields[("POINT_DATA","node_uid")])==[uid]
            @test parsed.fields[("POINT_DATA","step_0")]==reshape([1.25],1,1)
            element=A.view.add("wide element tensor")
            requested=(UInt64(1)<<31)+UInt64(301)
            A.mesh.add_elements_by_type(7,15,[requested],[uid])
            actual=only(tag for tag in only(A.mesh.get_elements(0,7)[2]) if tag>(UInt64(1)<<53))
            A.view.add_model_data(element,0,"wide point","ElementData",[actual],[collect(1.:9.)],.5,9)
            @test A.view.get_model_data(element,0)==(["ElementData"],[actual],collect(1.:9.))
            path=joinpath(directory,"wide-tensor.vtk");A.view.write(element,path);parsed=read_vtk(path)
            position=only(findall(==(actual),vec(parsed.fields[("CELL_DATA","element_uid")])))
            @test parsed.fields[("CELL_DATA","step_0")][:,position]==collect(1.:9.)
        end
    finally
        A.finalize()
    end
end

@testset "complete Pyramid14 exports as VTK type 74 and CSV retains all nodes" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("pyramid");A.model.add_discrete_entity(3,7)
        _,dimension,_,count,reference,_=A.mesh.get_element_properties(14)
        tags=collect(101:100+count)
        A.mesh.add_nodes(3,7,tags,reference);A.mesh.add_elements_by_type(7,14,[301],tags)
        view=A.view.add("pyramid values")
        A.view.add_model_data(view,0,"pyramid","NodeData",tags,[[Float64(tag)] for tag in tags])
        mktempdir() do directory
            path=joinpath(directory,"pyramid.vtk");A.view.write(view,path)
            lines=readlines(path)
            @test occursin("CELL_TYPES",join(lines,' ')) && occursin(r"\b74\b",join(lines,' '))
            @test any(line->occursin("5.1",line)||occursin("Lagrange",line),lines[1:4])
            path=joinpath(directory,"pyramid.csv");A.view.write(view,path)
            rows=split.(readlines(path),',');header=first(rows);column=Dict(name=>index for (index,name) in enumerate(header))
            cell=only(row for row in rows[2:end] if row[column["record"]]=="cell")
            @test parse.(Int,split(cell[column["connectivity"]],';'))==tags && length(tags)==14
        end
    finally
        A.finalize()
    end
end
end
