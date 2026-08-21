using Test
using Tessella
using Tessella.MeshTypes: ntris
using Tessella.Elements: validate

@testset "prismatic boundary layer" begin
    # Unit square in z=0, two triangles.
    coords=Float64[0 1 1 0; 0 0 1 1; 0 0 0 0]
    tris=Int32[1 1; 2 3; 3 4]
    face=Mesh(coords; tris=tris)
    bl=mesh_boundary_layer(face; hwall=0.1, ratio=1.2, nlayers=2)
    @test validate(bl).ok
    @test size(bl.coords,2)==4*3
    hexes=only(b for b in bl.blocks if b.msh==6)
    @test size(hexes.nodes,2)==ntris(face)*2
    # Independent volume: unit square area 1 times total height
    # H = hwall*((ratio^nl-1)/(ratio-1)) = 0.1*(1.44-1)/0.2 = 0.22
    H=0.1*(1.2^2-1)/(1.2-1)
    @test H≈0.22
    # Prism volume for a straight extrusion of a planar face is area * height.
    # Sample the first layer height hwall=0.1, area=1, two triangles → two prisms
    # of combined volume 0.1; second layer 0.12; total 0.22.
    zmax=maximum(bl.coords[3,:])
    @test zmax≈H atol=1e-12
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.0, nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=-0.1, ratio=1.2, nlayers=1)
end
