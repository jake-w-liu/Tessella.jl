module VTKLagrangeCurveChecks
using Test, Tessella
const A=Tessella.API

function read_curve(path,count)
    lines=readlines(path)
    start=findfirst(startswith("POINTS "),lines)
    points=hcat((parse.(Float64,split(lines[start+row])) for row in 1:count)...)
    cell=split(lines[findfirst(startswith("CELLS "),lines)+1])
    slots=parse.(Int,cell[2:end]).+1
    field(name,T;components=1)=let row=findfirst(startswith(name*" "),lines)
        reshape(parse.(T,split(join(lines[row+1:row+count],' '))),components,count)
    end
    return (;lines,points,slots,width=parse(Int,first(cell)),field,
        type=parse(Int,lines[findfirst(startswith("CELL_TYPES "),lines)+1]))
end

@testset "complete Lagrange curves retain curved geometry and nodal values" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            for (degree,msh) in ((3,26),(10,66))
                name="curve-$degree";A.model.add(name);A.model.add_discrete_entity(1,7)
                parameters=vcat([0.,1.],[index/degree for index in 1:degree-1])
                coordinates=hcat(([2+3t,-1+t^2,.5+t^3] for t in parameters)...)
                tags=collect(101:101+degree)
                A.mesh.add_nodes(1,7,tags,vec(coordinates))
                A.mesh.add_elements_by_type(7,msh,[301],tags)
                values=[[1+t,2-t^2] for t in parameters]
                later=[row.+[3.,5.] for row in values]
                view=A.view.add("curve vectors")
                A.view.add_model_data(view,0,name,"NodeData",tags,values,.25,2)
                A.view.add_model_data(view,2,name,"NodeData",tags,later,1.5,2)
                before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0))
                path=joinpath(directory,"curve-$degree.vtk");A.view.write(view,path)
                parsed=read_curve(path,degree+1)
                @test parsed.type==68 && parsed.width==degree+1
                @test parsed.points[:,parsed.slots]==coordinates
                @test vec(parsed.field("node_uid",UInt64))[parsed.slots]==UInt64.(tags)
                @test parsed.field("step_0",Float64;components=2)[:,parsed.slots]==hcat(values...)
                @test parsed.field("step_2",Float64;components=2)[:,parsed.slots]==hcat(later...)
                @test all(==(1),parsed.field("present_0",Int)) && all(==(1),parsed.field("present_2",Int))
                @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(view,0))==before
                nodal=A.view.add("discontinuous curve vectors")
                A.view.add_model_data(nodal,0,name,"ElementNodeData",[301],[vcat(later...)],.5,2)
                path=joinpath(directory,"nodal-$degree.vtk");A.view.write(nodal,path)
                parsed=read_curve(path,degree+1)
                @test parsed.type==68 && parsed.width==degree+1 && parsed.points[:,parsed.slots]==coordinates
                @test vec(parsed.field("node_uid",UInt64))[parsed.slots]==UInt64.(tags)
                @test parsed.field("step_0",Float64;components=2)[:,parsed.slots]==hcat(later...)
            end
        end
    finally
        A.finalize()
    end
end
end
