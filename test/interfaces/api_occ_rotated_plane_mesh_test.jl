module OCCRotatedPlaneMeshChecks
using Test,Tessella
const A=Tessella.API

function surface_mesh_area()
    tags,xyz,_=A.mesh.get_nodes()
    index=Dict(tag=>i for (i,tag) in enumerate(tags))
    area=0.0
    types,_,nodes=A.mesh.get_elements(2)
    for (k,type) in enumerate(types)
        type==2 || continue
        list=nodes[k]
        for j in 1:3:length(list)
            a,b,c=(index[list[j]],index[list[j+1]],index[list[j+2]])
            p=ntuple(i->xyz[3a-3+i],3);q=ntuple(i->xyz[3b-3+i],3)
            r=ntuple(i->xyz[3c-3+i],3)
            u=q.-p;w=r.-p
            area+=0.5*sqrt((u[2]*w[3]-u[3]*w[2])^2+(u[3]*w[1]-u[1]*w[3])^2+
                           (u[1]*w[2]-u[2]*w[1])^2)
        end
    end
    return area
end

# The vertices of a rotated planar face are stored a rounding residual off
# their common plane. The planar mesher's exact coplanarity predicate used to
# reject them ("Point[3] is not coplanar with Surface[1]"), so a generically
# rotated box — or the caps of a cylinder with a tilted axis — could not be
# meshed at all.
@testset "OCC planar faces mesh after a generic rotation" begin
    cases=(
        ("rotated box",10.0,occ->begin
            box=occ.add_box(0,0,0,2,1,1)
            occ.rotate([(3,box)],0.1,0.2,0.3,1,1,0.5,0.7)
        end),
        ("obliquely mirrored box",10.0,occ->begin
            box=occ.add_box(0,0,0,2,1,1)
            occ.mirror([(3,box)],1,1,0.5,0.25)
        end),
        ("rotated wedge",3+sqrt(2),occ->begin
            wedge=occ.add_wedge(0,0,0,1,1,1)
            occ.rotate([(3,wedge)],0.1,0.2,0.3,1,1,0.5,0.7)
        end),
        ("rotated rectangle",2.0,occ->begin
            rectangle=occ.add_rectangle(0,0,0,2,1)
            occ.rotate([(2,rectangle)],0.1,0.2,0.3,1,1,0.5,0.7)
        end),
    )
    for (name,area,build) in cases
        @testset "$name" begin
            A.initialize()
            try
                A.option("General.Terminal",0)
                build(A.model.occ)
                A.model.occ.synchronize()
                A.option("Mesh.MeshSizeMax",0.5)
                A.mesh.generate(2)
                @test surface_mesh_area()≈area rtol=1e-9
            finally
                A.finalize()
            end
        end
    end
    # a point genuinely off the plane is still rejected
    M=Tessella.Model
    square=NTuple{3,Float64}[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,1e-9)]
    @test_throws ArgumentError M._model_surface_projection(
        square,[1,2,3,4],1,"test")
    tilted=NTuple{3,Float64}[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,3e-17)]
    @test M._model_surface_projection(tilted,[1,2,3,4],1,"test")[4]==(1,2)
end
end
