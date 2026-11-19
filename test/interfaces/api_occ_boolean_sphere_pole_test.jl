module OCCBooleanSpherePoleChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A box corner at the sphere center: the three box faces through the corner
# cut the sphere along its equator, its seam meridian and a second meridian
# that runs into the pole. The kept spherical piece is the octant bounded by
# those three arcs and a quarter of the pole's degenerate edge. Counts are
# Gmsh 4.15.2's; volume, area and centroid are closed forms.
#
# The one Boolean serves three regressions: the section through the pole
# (the shell used not to close), the measure of the trimmed sphere face (it
# reported the whole sphere, 28.27) and its bounding box (the whole ball).
@testset "OCC Boolean sections through a sphere pole" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(0,0,0,1.5)
        out,images=occ.intersect([(3,box)],[(3,ball)])
        @test out==[(3,1)]
        @test images==[[(3,1)],[(3,1)]]
        @test Tuple(length(occ.get_entities(d)) for d in 0:3)==(4,7,4,1)
        @test occ.get_mass(3,1)≈π*1.5^3/6 rtol=1e-8
        kernel=A.CURRENT[].occ_kernel
        spherical=[id for ((d,_),id) in kernel.bindings if d==2 &&
            get(kernel.shapes[id].geometry.geometry,:occ,nothing)===:sphere]
        @test length(spherical)==1
        face=kernel.shapes[only(spherical)]
        # three arcs and a quarter of the pole edge
        edges=kernel.shapes[only(face.children)].children
        @test length(edges)==4
        spans=[g.t1-g.t0 for g in (kernel.shapes[e].geometry.geometry
                                   for e in edges) if g.occ===:degenerate]
        @test spans≈[π/2] atol=1e-12
        area,centroid=M._model_occ_face_measure(kernel,only(spherical))
        @test area≈π*1.5^2/2 rtol=1e-6
        @test collect(centroid)≈fill(0.75,3) atol=1e-5
        @test collect(M._model_occ_shape_bbox(kernel,only(spherical)))≈
            [0,0,0,1.5,1.5,1.5] atol=1e-9
    finally
        A.finalize()
    end
end
end
