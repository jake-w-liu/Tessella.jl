module ModelMetricDelaunayArmChecks
using Test,Tessella,LinearAlgebra,Statistics
const K=Tessella.Model
const N=Tessella.MeshTypes

# Unit square with straight boundary Lines; transfinite N=6 fixes the
# boundary sampling so the only free output is the interior triangulation.
function rectangle(;counts=6)
    model=GeoModel()
    points=[add_point!(model,x,y,0.;mesh_size=.5) for (x,y) in
            ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    curves=[add_line!(model,points[i],points[mod1(i+1,4)]) for i in 1:4]
    for curve in curves
        K.set_transfinite_curve!(model,curve,counts)
    end
    loop=add_curve_loop!(model,curves)
    tag=add_plane_surface!(model,[loop])
    return model,tag
end
function area(mesh)
    sum(norm(cross(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
                   mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]]))/2 for cell in eachcol(mesh.tris))
end
# Mean |dx|/|dy| over unique interior edges: under m11=16 an x-aligned chord
# is metric-expensive, so a tensor-aware Delaunay arm elongates its cells
# along the weak direction (y) — the anisotropy ratio drops below 1.
function edge_anisotropy(mesh)
    seen=Set{Tuple{Int32,Int32}}()
    dx=Float64[];dy=Float64[]
    for cell in eachcol(mesh.tris),(a,b) in ((1,2),(2,3),(3,1))
        key=minmax(cell[a],cell[b])
        key in seen && continue
        push!(seen,key)
        pa,pb=mesh.coords[:,key[1]],mesh.coords[:,key[2]]
        push!(dx,abs(pa[1]-pb[1]));push!(dy,abs(pa[2]-pb[2]))
    end
    return mean(dx)/mean(dy)
end

@testset "algorithm 5 Delaunay arm consumes a tensor field" begin
    # xfield and yfield scalarize to the same background size (min-direction
    # h=0.25); only the transported tensor differs, so any output difference
    # is attributable to the metric-Delaunay cavity and circumcenter tests.
    xmodel,xtag=rectangle();xmodel.meshing.algorithm_2d=5
    ymodel,ytag=rectangle();ymodel.meshing.algorithm_2d=5
    xmesh=mesh_model_surface(xmodel,xtag;
        size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
    ymesh=mesh_model_surface(ymodel,ytag;
        size_field=MathEvalAnisoField(m11="1",m22="16",m33="1"))
    @test validate(xmesh).ok && validate(ymesh).ok
    @test N.ntris(xmesh)>0 && N.ntris(ymesh)>0
    @test isapprox(area(xmesh),1.;atol=1e-12) && isapprox(area(ymesh),1.;atol=1e-12)
    @test xmesh.tris!=ymesh.tris || !isequal(xmesh.coords,ymesh.coords)
    # Directional elongation follows the transported tensor, not a constant:
    # m11=16 suppresses x-extent while m22=16 suppresses y-extent.
    xr=edge_anisotropy(xmesh);yr=edge_anisotropy(ymesh)
    @test xr<1<yr
    @test xr<0.75 && yr>1.5

    # The scalar arm on the same fixture is unchanged by the tensor path.
    model,tag=rectangle();model.meshing.algorithm_2d=5
    scalar=mesh_model_surface(model,tag)
    @test validate(scalar).ok && isapprox(area(scalar),1.;atol=1e-12)
end

@testset "algorithm 6 frontal arm consumes the same tensor" begin
    xmodel,xtag=rectangle();xmodel.meshing.algorithm_2d=6
    ymodel,ytag=rectangle();ymodel.meshing.algorithm_2d=6
    xmesh=mesh_model_surface(xmodel,xtag;
        size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
    ymesh=mesh_model_surface(ymodel,ytag;
        size_field=MathEvalAnisoField(m11="1",m22="16",m33="1"))
    @test validate(xmesh).ok && validate(ymesh).ok
    @test isapprox(area(xmesh),1.;atol=1e-12) && isapprox(area(ymesh),1.;atol=1e-12)
    @test xmesh.tris!=ymesh.tris || !isequal(xmesh.coords,ymesh.coords)
    @test edge_anisotropy(xmesh)<edge_anisotropy(ymesh)
end

end
