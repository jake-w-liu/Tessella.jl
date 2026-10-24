# `Mesh.Algorithm=2` (AUTO) per-face resolution.
#
# Verified on the Gmsh 4.15.2 binary: AUTO picks Delaunay on planar faces
# and MeshAdapt on every non-planar class — Cylinder, Sphere, ruled
# Surface and BSpline all report "Meshing surface N (<kind>, MeshAdapt)".
# No native MeshAdapt producer exists, so the implemented Frontal-Delaunay
# construction serves the off-plane admission: upstream's outcome (a
# valid mesh) through the closest native producer. Selectors that name
# machinery without a native producer (1/4/8/9/10/11) keep their loud
# rejection.

using Test
using Tessella
const M=Tessella.Model
const MT=Tessella.MeshTypes

_error(f)=try;f();nothing catch err;err end

function _square_model(;ruled=false)
    m=GeoModel()
    for (tag,(x,y)) in enumerate(((0.0,0.0),(2.0,0.0),
                                  (2.0,2.0),(0.0,2.0)))
        add_point!(m,x,y,0;tag=tag)
    end
    for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        add_line!(m,a,b;tag=tag)
    end
    add_curve_loop!(m,[1,2,3,4];tag=1)
    ruled ? add_ruled_surface!(m,[1];tag=2) : add_plane_surface!(m,[1];tag=1)
    return m
end

function _sector_model()
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0;tag=1,angle=pi)
    M._model_occ_synchronize!(m)
    return m,findfirst(==(:cylinder),m.surface_types)
end

_tri_area(a,b,c)=0.5*hypot(
    (b[2]-a[2])*(c[3]-a[3])-(b[3]-a[3])*(c[2]-a[2]),
    (b[3]-a[3])*(c[1]-a[1])-(b[1]-a[1])*(c[3]-a[3]),
    (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))

function _mesh_area(mesh)
    total=0.0
    for cell in 1:MT.ntris(mesh)
        a,b,c=mesh.tris[:,cell]
        total+=_tri_area(mesh.coords[:,a],mesh.coords[:,b],mesh.coords[:,c])
    end
    return total
end

@testset "AUTO on a planar face selects the Delaunay arm" begin
    # Upstream: "Meshing surface 1 (Plane, Delaunay)". Tessella resolves
    # `:plane` to algorithm 5 — the PSLG refine fills a valid 4.0 mesh.
    m=_square_model()
    m.meshing.algorithm_2d=2
    mesh=M.mesh_model_surface(m,1)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    @test _mesh_area(mesh)≈4.0 atol=1e-6
end

@testset "AUTO on a ruled face selects the Frontal arm" begin
    # Upstream's AUTO picks MeshAdapt on a ruled "Surface" — the native
    # stand-in is Frontal-Delaunay, which fills a coplanar-bounded ruled
    # patch through the `:frontal` refine arm.
    m=_square_model(ruled=true)
    m.meshing.algorithm_2d=2
    mesh=M.mesh_model_surface(m,2)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    @test _mesh_area(mesh)≈4.0 atol=1e-6
end

@testset "AUTO on a curved OCC face meshes the (u,v) domain" begin
    # Upstream: "(Cylinder, MeshAdapt)" and a valid tessellation. The
    # parametric kernel is the native producer for the Frontal admission —
    # a half-sector lateral (r=1, h=2, θ∈[0,π]) yields r·h·θ = 2π area,
    # with the flat facets under-approximating the curved band by ~2%.
    m,lat=_sector_model()
    m.meshing.algorithm_2d=2
    mesh=M.mesh_model_surface(m,lat)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    @test _mesh_area(mesh)≈2π atol=0.15
end

@testset "explicit producer-less selectors still reject under AUTO" begin
    # AUTO itself is covered above; the unchanged explicit-algorithm
    # contract keeps its loud rejection on the same curved face.
    for algorithm in (1,4,8,9,10,11)
        m,lat=_sector_model()
        m.meshing.algorithm_2d=algorithm
        err=_error(()->M.mesh_model_surface(m,lat))
        @test err isa ArgumentError
        @test occursin("not implemented",sprint(showerror,err))
    end
end
