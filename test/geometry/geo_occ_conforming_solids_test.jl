# `.geo` meshing fills OpenCASCADE primitive and Boolean solids from their own
# boundary Surface meshes, like the API session generator, so `Mesh.MeshSizeMax`
# is honoured on an OCC sphere (Gmsh 4.15.2: 380 surface triangles at 0.3)
# instead of the fixed analytic tessellation.
module GeoOccConformingSolidsTest
using Test, Tessella

@testset "geo OCC solids mesh from their boundary surfaces" begin
    dir=mktempdir()
    path=joinpath(dir,"sphere.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Sphere(1) = {0,0,0,1};
        Mesh.MeshSizeMax = 0.3;
        """)
    ex=execute_geo(path;mesh_dim=2)
    @test isempty(ex.errors)
    @test ex.model.meshing.conforming_solids
    ntri=sum(size(part[3].tris,2) for part in ex.mesh_parts
             if part[1]==2 && part[3] isa Tessella.MeshTypes.Mesh;init=0)
    @test 300<=ntri<=500
end
end
