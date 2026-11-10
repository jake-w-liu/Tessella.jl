module PluginDistanceThin3DTests
using Test,Tessella
const A=Tessella.API

# Eight Hex8 cells cover box(1,h,h), all six boundary Quad4 panels fixed.
function thin_box(h;offset=0.)
    id(i,j,k)=1+i+3j+9k
    A.model.add("thin diffusion box")
    for tag in 1:6;A.model.add_discrete_entity(2,tag);end
    A.model.add_discrete_entity(3,8)
    points=[(offset+i/2,j*h/2,k*h/2) for k in 0:2 for j in 0:2 for i in 0:2]
    A.mesh.add_nodes(3,8,collect(1:27),collect(Iterators.flatten(points)))
    cells=Vector{Int}[]
    for k in 0:1,j in 0:1,i in 0:1
        push!(cells,[id(i,j,k),id(i+1,j,k),id(i+1,j+1,k),id(i,j+1,k),
            id(i,j,k+1),id(i+1,j,k+1),id(i+1,j+1,k+1),id(i,j+1,k+1)])
    end
    faces=Vector{Vector{Int}}[]
    for axis in 1:3,side in (0,2)
        panel=Vector{Int}[]
        for a in 0:1,b in 0:1
            ids=Int[]
            for (u,v) in ((a,b),(a+1,b),(a+1,b+1),(a,b+1))
                xyz=axis==1 ? (side,u,v) : axis==2 ? (u,side,v) : (u,v,side)
                push!(ids,id(xyz...))
            end
            push!(panel,ids)
        end
        push!(faces,panel)
    end
    for (tag,panel) in enumerate(faces)
        A.mesh.add_elements_by_type(tag,3,collect(200+4tag:203+4tag),
            collect(Iterators.flatten(panel)))
    end
    A.mesh.add_elements_by_type(8,5,collect(301:308),collect(Iterators.flatten(cells)))
    return points
end

function expected_distance(h)
    return setprecision(BigFloat,512) do
        H=BigFloat(h);mu=H*sqrt(1+2H^2)
        s=(H^2/8)/(mu^2*(4/9)*(H^2+2)+H^2/27)
        Float64(-mu*log1p(-s))
    end
end

@testset "Thin 3D diffusion retains representable distances" begin
    for h in (1e-20,1e-200)
        for offset in (0.,10.)
            A.initialize()
            try
                A.option("General.Terminal",0)
                thin_box(h;offset=offset)
                nodes=A.mesh.get_nodes()
                elements=A.mesh.get_elements()
                model=A.model.get_current()
                A.plugin.set_number("Distance","DistanceType",h)
                A.plugin.run("Distance")
                types,data=A.view.get_list_data(only(A.view.get_tags()))
                found=Float64[]
                center=(offset+.5,.5h,.5h)
                for (type,flat) in zip(types,data)
                    @test type=="SH"
                    for base in 0:32:length(flat)-1,node in 1:8
                        point=ntuple(axis->flat[base+8(axis-1)+node],3)
                        point==center && push!(found,flat[base+24+node])
                    end
                end
                expected=expected_distance(h)
                @test length(found)==8
                @test all(value->isapprox(value,expected;rtol=2e-11),found)
                @test expected>0.
                @test A.mesh.get_nodes()==nodes
                @test A.mesh.get_elements()==elements
                @test A.model.get_current()==model
            finally
                A.finalize()
            end
        end
    end
end
end
