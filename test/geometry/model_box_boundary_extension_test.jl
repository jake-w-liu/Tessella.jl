module ModelBoxBoundaryExtensionChecks
using Test, Tessella
const M=Tessella.Model
const N=Tessella.MeshTypes
const F=Tessella.SizeField

# A Box is filled from its six face meshes by Gmsh's Delaunay insertion, like
# any explicit Volume. Until 0.1.3 it was filled by bisecting every edge
# longer than twice the size, which this file asserted as `edge <= 2h` on
# every edge; Gmsh itself does not satisfy that (its longest edge on this
# model is 2.1 h) and the bisection gave several times its tetrahedra.
@testset "A Box extends boundary sizing without a selected field" begin
    model=GeoModel()
    M.add_box!(model,0.,0.,0.,1.,1.,1.;tag=1)
    model.meshing.algorithm_2d=6
    for point in keys(model.points)
        model.point_size[point]=.25
    end
    surface=M._model_explicit_volume_geometry(model,1,"boundary extension regression").surface
    field=M._volume_boundary_size_field(model,1,surface,"boundary extension regression")
    mesh=mesh_model_volume(model,1)
    @test validate(mesh).ok
    volume=sum(N.tet_volume(N.node(mesh,cell[1]),N.node(mesh,cell[2]),
                            N.node(mesh,cell[3]),N.node(mesh,cell[4]))
               for cell in eachcol(mesh.tets))
    @test volume≈1.0 atol=1e-12
    @test N.nnodes(mesh)>N.nnodes(surface)+1
    edges=Set(minmax(Int(cell[i]),Int(cell[j])) for cell in eachcol(mesh.tets)
              for i in 1:3 for j in i+1:4)
    ratios=map(collect(edges)) do (a,b)
        p=N.node(mesh,a);q=N.node(mesh,b)
        midpoint=ntuple(i->(p[i]+q[i])/2,3)
        h=minimum(F.size_at(field,r...) for r in (p,midpoint,q))
        return hypot((q[1]-p[1],q[2]-p[2],q[3]-p[3])...)/h
    end
    sort!(ratios)
    # The interior follows the size of the boundary: the typical edge is the
    # prescribed size and none is far from it (Gmsh 4.15.2 on the same Box:
    # median 1.03 h, longest 2.1 h, 386 tetrahedra on its own 264 triangles).
    @test 0.8<=ratios[end÷2]<=1.25
    @test ratios[end]<=3
    @test size(surface.tris,2)<size(mesh.tets,2)<2.5*size(surface.tris,2)
    actual=Set(Tuple(sort([N.node(mesh,v) for v in cell]))
               for cell in first(N.boundary_faces(mesh.tets)))
    expected=Set(Tuple(sort([N.node(surface,v) for v in cell]))
                 for cell in eachcol(surface.tris))
    @test actual==expected
end
end
