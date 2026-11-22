# A cylinder converted by `occ.convert_to_nurbs` meshes in 3D (it does in
# Gmsh 4.15.2). Its converted cap keeps its plane record beside the patch
# and has no vertex table; the surface producer read `vertex_uv` from it
# ("type NamedTuple has no field vertex_uv").
module ApiOccNurbsMeshTest
using Test, Tessella

@testset "occ.convert_to_nurbs cylinder meshes" begin
    A=Tessella.API; occ=A.model.occ
    A.initialize()
    try
        A.option("General.Terminal",0)
        cyl=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,4.0,1.5)
        occ.convert_to_nurbs([(3,cyl)])
        occ.synchronize()
        A.option("Mesh.MeshSizeMax",1.0)
        A.mesh.generate(3)
        ets,etags,_=A.mesh.get_elements(3,cyl)
        @test ets==Int32[4]
        @test length(etags[1])>0
    finally
        A.finalize()
    end
end
end
