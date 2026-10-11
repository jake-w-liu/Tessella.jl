module VTKSerendipityChecks
using Test, Tessella
const A=Tessella.API

function expected_topology(msh)
    if msh==16
        vertices=[0. 1. 1. 0.;0. 0. 1. 1.;0. 0. 0. 0.]
        edges=((1,2),(2,3),(3,4),(4,1));type=23
    elseif msh==17
        vertices=[0. 1. 1. 0. 0. 1. 1. 0.;0. 0. 1. 1. 0. 0. 1. 1.;0. 0. 0. 0. 1. 1. 1. 1.]
        edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8));type=25
    elseif msh==18
        vertices=[0. 1. 0. 0. 1. 0.;0. 0. 1. 0. 0. 1.;0. 0. 0. 1. 1. 1.]
        edges=((1,2),(2,3),(3,1),(4,5),(5,6),(6,4),(1,4),(2,5),(3,6));type=26
    else
        vertices=[0. 1. 1. 0. .5;0. 0. 1. 1. .5;0. 0. 0. 0. 1.]
        edges=((1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5));type=27
    end
    return type,hcat(vertices,((vertices[:,a]+vertices[:,b])/2 for (a,b) in edges)...)
end

function read_cell(path,count)
    lines=readlines(path);start=findfirst(startswith("POINTS "),lines)
    points=hcat((parse.(Float64,split(lines[start+row])) for row in 1:count)...)
    cell=split(lines[findfirst(startswith("CELLS "),lines)+1]);slots=parse.(Int,cell[2:end]).+1
    field(name,T;components=1)=let row=findfirst(startswith(name*" "),lines)
        reshape(parse.(T,split(join(lines[row+1:row+count],' '))),components,count)
    end
    return (;points,slots,field,width=parse(Int,first(cell)),
        type=parse(Int,lines[findfirst(startswith("CELL_TYPES "),lines)+1]))
end

world(reference)=[2. 1. .5;.5 3. 1.;1. .5 4.]*reference .+ [7.,-2.,5.]
value(point)=[point[1]+2point[2]-point[3],point[3]^2+point[1]]

@testset "quadratic serendipity VTK keeps corner and edge geometry and values" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            for msh in (16,17,18,19)
                name="serendipity-$msh";A.model.add(name)
                _,dimension,_,count,local_coordinates,_=A.mesh.get_element_properties(msh)
                reference=vcat(reshape(local_coordinates,dimension,count),zeros(3-dimension,count))
                if msh in (16,19)
                    reference[1:2,:].=(reference[1:2,:].+1)./2
                elseif msh==17
                    reference.=(reference.+1)./2
                else
                    reference[3,:].=(reference[3,:].+1)./2
                end
                coordinates=world(reference);tags=collect(101:100+count)
                A.model.add_discrete_entity(dimension,7)
                A.mesh.add_nodes(dimension,7,tags,vec(coordinates))
                A.mesh.add_elements_by_type(7,msh,[301],tags)
                type,target=expected_topology(msh);expected=world(target)
                lookup=Dict(Tuple(coordinates[:,node])=>tags[node] for node in 1:count)
                expected_tags=[lookup[Tuple(point)] for point in eachcol(expected)]
                values=[value(point) for point in eachcol(coordinates)]
                view=A.view.add("serendipity vectors")
                A.view.add_model_data(view,0,name,"NodeData",tags,values,.25,2)
                before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0))
                path=joinpath(directory,"nodes-$msh.vtk");A.view.write(view,path);parsed=read_cell(path,count)
                @test parsed.type==type && parsed.width==count
                @test parsed.points[:,parsed.slots]==expected
                @test vec(parsed.field("node_uid",UInt64))[parsed.slots]==UInt64.(expected_tags)
                @test parsed.field("step_0",Float64;components=2)[:,parsed.slots]==hcat((value(point) for point in eachcol(expected))...)
                @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0))==before
                nodal=A.view.add("serendipity element vectors")
                A.view.add_model_data(nodal,0,name,"ElementNodeData",[301],[vcat(values...)],.5,2)
                path=joinpath(directory,"nodal-$msh.vtk");A.view.write(nodal,path);parsed=read_cell(path,count)
                @test parsed.type==type && parsed.width==count && parsed.points[:,parsed.slots]==expected
                @test vec(parsed.field("node_uid",UInt64))[parsed.slots]==UInt64.(expected_tags)
                @test parsed.field("step_0",Float64;components=2)[:,parsed.slots]==hcat((value(point) for point in eachcol(expected))...)
            end
        end
    finally
        A.finalize()
    end
end
end
