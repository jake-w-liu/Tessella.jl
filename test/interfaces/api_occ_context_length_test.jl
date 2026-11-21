module OCCContextLengthChecks
using Test, Tessella
const A=Tessella.API

# Gmsh sizes a Curve whose end Points carry no size at CTX::lc/10
# (`LC_MVertex_PNTS`), so with only `Mesh.MeshSizeMax` set every coarse mesh
# depends on CTX::lc, the diagonal of the model's bounding box. For an
# OpenCASCADE model that box is the union of each edge's and face's own
# `BRepBndLib` box. The model box used here was the vertex cloud plus the
# primitive records (read with the cylinder's base centre as its mid-height
# centre), so a torus (one vertex) had lc 2.83 for 7.72 and a cylinder 3.32
# for 2.45, and their coarse meshes had up to twice Gmsh's triangles or half
# of them. Each `lc` below is `General.BoundingBoxSize` of Gmsh 4.15.2 after
# `occ.synchronize()`; the segment counts are its `generate(1)` at
# `Mesh.MeshSizeMax=0.4`. (A cylinder cut by a cylinder, 2.449490069, and two
# fused spheres, 4.123105965, agree too; a Boolean is too slow to run here
# without compilation.)
function with_model(f,build)
    A.initialize()
    try
        A.option("General.Terminal",0)
        build(A.model.occ)
        A.model.occ.synchronize()
        return f()
    finally
        A.finalize()
    end
end

@testset "CTX::lc of an OpenCASCADE model is Gmsh's" begin
    converted(build)=o->o.convert_to_nurbs([(3,build(o))])
    for (name,build,lc) in (
            ("cylinder",o->o.add_cylinder(0,0,0,0,0,2,0.5),2.449490069),
            ("cone",o->o.add_cone(0,0,0,0,0,2,1,0.5),3.464101962),
            ("torus",o->o.add_torus(0,0,0,2,0.5),7.718720651),
            ("sphere",o->o.add_sphere(0,0,0,1),3.464101962),
            ("box",o->o.add_box(0,0,0,1,2,3),3.741657707),
            ("disk",o->o.add_disk(0,0,0,1,0.5),2.236068246),
            ("partial sphere",o->o.add_sphere(0,0,0,1;angle1=-0.5,angle2=1.0,angle3=2.0),2.179504756),
            ("partial torus",o->o.add_torus(0,0,0,2,0.5;angle=2.0),4.796567649),
            ("tilted cylinder",o->o.add_cylinder(1,2,3,1,1,1,0.5),3.146264716),
            ("NURBS cylinder",converted(o->o.add_cylinder(0,0,0,0,0,2,0.5)),3.041381609),
            ("NURBS torus",converted(o->o.add_torus(0,0,0,2,0.5)),11.586630535))
        @testset "$name" begin
            @test with_model(()->A.NATIVE_CONTEXT_LC[],build)≈lc atol=2e-9
        end
    end
end

@testset "coarse Curves of the curved primitives have Gmsh's segments" begin
    segments()=[sum(length,A.mesh.get_elements(1,tag)[2];init=0)
                for (_,tag) in A.model.get_entities(1)]
    for (name,build,expected) in (
            ("cylinder",o->o.add_cylinder(0,0,0,0,0,2,0.5),[13,9,13]),
            ("cone",o->o.add_cone(0,0,0,0,0,2,1,0.5),[10,6,19]),
            ("torus",o->o.add_torus(0,0,0,2,0.5),[40,8]))
        @testset "$name" begin
            counts=with_model(build) do
                A.option("Mesh.MeshSizeMax",0.4)
                A.mesh.generate(1)
                segments()
            end
            @test counts==expected
        end
    end
end

end
