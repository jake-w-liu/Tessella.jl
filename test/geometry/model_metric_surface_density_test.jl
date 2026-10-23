module ModelMetricSurfaceDensityChecks
using Test,Tessella,LinearAlgebra
const K=Tessella.Model
const N=Tessella.MeshTypes
const F=Tessella.SizeField

# The conforming tensor producer must fill a planar Surface at the density
# the metric actually prescribes: a pure midpoint-bisection fill can satisfy
# every edge certificate while collapsing all directions toward the smallest
# prescribed axis, overshooting the anisotropic count by tens of times. The
# metric-Delaunay Steiner insertions plus the certified flip/removal
# coarsening sweeps exist to keep the count inside the density band the
# field implies. This set pins that band on constant, graded and sheared
# metrics — the regression is the collapse, not the precise constant.
function square(;extent=2.0,size=0.045)
    model=GeoModel()
    a=-extent;b=extent
    points=[add_point!(model,x,y,0.0;mesh_size=size)
            for (x,y) in ((a,a),(b,a),(b,b),(a,b))]
    curves=[add_line!(model,points[index],points[mod1(index+1,4)])
            for index in 1:4]
    tag=add_plane_surface!(model,[add_curve_loop!(model,curves)])
    return model,tag
end

function edges_of(mesh)
    seen=Set{NTuple{2,Int32}}()
    for cell in eachcol(mesh.tris),pair in
            ((cell[1],cell[2]),(cell[2],cell[3]),(cell[3],cell[1]))
        push!(seen,minmax(pair...))
    end
    return seen
end

function metric_bound(mesh,field,tag,target)
    all(edges_of(mesh)) do (a,b)
        F.metric_edge_length(field,Tuple(mesh.coords[:,a]),
            Tuple(mesh.coords[:,b]);entity=(2,tag))<=target
    end
end

# Optimal unit-metric triangles carry metric area sqrt(3)/4, so the count a
# constant-or-slowly-varying field prescribes is 4*area*sqrt(detM)/sqrt(3).
# A certified-but-collapsed fill overshoots by the aspect ratio.
theory(area,sqrtdet)=4area*sqrtdet/sqrt(3.0)

@testset "anisotropic metric surface keeps prescribed density" begin
    model,tag=square()
    field=MathEvalAnisoField(m11="1",m22="400",m33="1")
    mesh=mesh_model_surface(model,tag,field;max_metric_edge=1.0)
    @test validate(mesh).ok
    @test metric_bound(mesh,field,tag,1.0)
    # sqrt(det)=20 over area 16 -> ~739 optimal triangles; the certified
    # collapse historically emitted ~15x that. Stay inside the parity band.
    expected=theory(16.0,20.0)
    @test 0.3expected<N.ntris(mesh)<5.0expected
    @test sum(N.triangle_area((Tuple(mesh.coords[:,row]) for row in cell)...)
              for cell in eachcol(mesh.tris))≈16.0 rtol=1e-9
end

@testset "graded and sheared metrics stay inside the band" begin
    model,tag=square()
    field=MathEvalAnisoField(m11="1",m22="(8*x+1)^2 + 1",m33="1")
    mesh=mesh_model_surface(model,tag,field;max_metric_edge=1.0)
    @test validate(mesh).ok
    @test metric_bound(mesh,field,tag,1.0)
    # Numerical density integral for m22=(8x+1)^2+1 on [-2,2]^2 is ~130
    # metric area units -> ~301 optimal triangles.
    @test 0.3*301.0<N.ntris(mesh)<5.0*301.0

    model,tag=square()
    field=MathEvalAnisoField(m11="1",m22="25",m33="1",m12="(x+y)/4")
    mesh=mesh_model_surface(model,tag,field;max_metric_edge=1.0)
    @test validate(mesh).ok
    @test metric_bound(mesh,field,tag,1.0)
    @test N.ntris(mesh)<1200
end
end
