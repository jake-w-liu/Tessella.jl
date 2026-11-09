# `Mesh.Algorithm` selector dispatch on OCC-curved volume boundary faces.
#
# Upstream `meshGFace` resolves `meshAttributes.algorithm`/`Mesh.Algorithm` on
# parametric OCC faces through the same admission as planar faces, so neither
# the `_model_occ_curved_surface_mesh` kernel (explicit shells) nor the
# primitive factories' analytic boundary (implicit volumes) may silently
# substitute their own construction for a selector naming different
# machinery. The Delaunay-family policies 3/5/6 are admitted; AUTO resolves
# per face to Delaunay on planes and local MeshAdapt on curved charts,
# including their periodic and collapsed-pole covers. Every selector
# without a native producer keeps the loud not-implemented rejection.
#
# Before the admission was wired, every `Mesh.Algorithm` value produced the
# identical mesh on a cylinder volume — a silent degradation.

module ModelOCCCurvedAlgorithmChecks
using Test
using Tessella
const M=Tessella.Model
const MT=Tessella.MeshTypes

function _sector(algorithm;per_face=nothing)
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0;tag=1,angle=pi)
    M._model_occ_synchronize!(m)
    m.meshing.algorithm_2d=algorithm
    per_face===nothing || M.set_algorithm!(m,2,per_face[1],per_face[2])
    return m
end

function _primitive_cylinder(algorithm)
    m=GeoModel()
    add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0;tag=1)
    m.meshing.algorithm_2d=algorithm
    return m
end

_error(f)=try;f();nothing catch err;err end

@testset "explicit OCC sector: curved boundary honors Mesh.Algorithm" begin
    # The sector's lateral is a `:cylinder` face meshed by the parametric
    # OCC kernel; its four caps are ordinary planes. 1 runs local adaptation,
    # while 5 and 6 fill a valid solid through the Delaunay-family admission.
    for algorithm in (1,5,6)
        mesh=M.mesh_model_volume(_sector(algorithm),1)
        @test MT.validate(mesh).ok
        @test MT.ntets(mesh)>0
    end

    # `set_algorithm!(2, lateral, 3)` suppresses the kernel's interior (u,v)
    # lattice seeds — upstream `ALGO_2D_INITIAL_ONLY` bounds the face to its
    # recovered 1-D boundary vertices — so the solid refines everywhere else
    # but loses exactly the lateral's interior nodes.
    lateral=findfirst(==(:cylinder),_sector(6).surface_types)
    refined=M.mesh_model_volume(_sector(6),1)
    initial=M.mesh_model_volume(_sector(6;per_face=(lateral,3)),1)
    @test MT.validate(initial).ok
    @test MT.ntets(initial)>0
    @test MT.nnodes(initial)<MT.nnodes(refined)

    # AUTO resolves per face like upstream: plane caps take Delaunay and
    # the curved lateral takes local MeshAdapt — the whole solid
    # fills instead of rejecting.
    auto=M.mesh_model_volume(_sector(2),1)
    @test MT.validate(auto).ok
    @test MT.ntets(auto)>0

    # Selectors without a native producer keep the loud rejection — the
    # lateral rejects 8, 9, 10 and 11 even though the caps would accept
    # several of them.
    for algorithm in (8,9,10,11)
        err=_error(()->M.mesh_model_volume(_sector(algorithm),1))
        @test err isa ArgumentError
        @test occursin("not implemented",sprint(showerror,err))
    end
    # Upstream's BAMG selector is now the certified OCC metric producer on
    # the curved lateral: it requires an anisotropic background field, so a
    # global algorithm-7 request with no field keeps a loud rejection —
    # naming the requirement rather than silently substituting the scalar
    # context — while the same request with a tensor field fills the solid.
    err=_error(()->M.mesh_model_volume(_sector(7),1))
    @test err isa ArgumentError
    @test occursin("anisotropic",sprint(showerror,err))
    aniso=MathEvalAnisoField(m11="1",m22="1",m33="1")
    mesh=M.mesh_model_volume(_sector(6;per_face=(lateral,7)),1;
        size_field=aniso)
    @test MT.validate(mesh).ok
    @test MT.ntets(mesh)>0
    # A per-face request names the rejecting surface in the diagnostic.
    err=_error(()->M.mesh_model_volume(_sector(6;per_face=(lateral,8)),1))
    @test err isa ArgumentError
    @test occursin("Surface[$lateral]",sprint(showerror,err))
end

@testset "primitive volumes route MeshAdapt through their actual faces" begin
    # The canned `cylinder_surface` boundary is the native analytic product
    # serving 3/5/6. MeshAdapt and non-planar AUTO consume generated faces;
    # every selector that
    # names machinery without a native producer now rejects loudly instead
    # of silently returning the same tessellation.
    # AUTO on this cylinder resolves to the same algorithm-1 producer; the
    # dedicated AUTO checks below and in the sector testset cover that route.
    for algorithm in (1,3,5,6)
        mesh=M.mesh_model_volume(_primitive_cylinder(algorithm),1)
        @test MT.validate(mesh).ok
        @test MT.ntets(mesh)>0
    end
    for algorithm in (7,8,9,10,11)
        err=_error(()->M.mesh_model_volume(_primitive_cylinder(algorithm),1))
        @test err isa ArgumentError
        @test occursin("not implemented",sprint(showerror,err))
    end
    # A per-face request on the materialized cap records rejects too.
    m=_primitive_cylinder(6)
    cap=findfirst(==(:plane),m.surface_types)
    M.set_algorithm!(m,2,cap,8)
    err=_error(()->M.mesh_model_volume(m,1))
    @test err isa ArgumentError
    @test occursin("Surface[$cap]",sprint(showerror,err))

    # AUTO still resolves per face: an all-plane primitive box takes the
    # Delaunay branch, while closed sphere/cone/torus charts adapt their
    # periodic covers and retain the actual lower-dimensional identities.
    box=GeoModel()
    add_box!(box,0.0,0.0,0.0,1.0,1.0,1.0;tag=1)
    box.meshing.algorithm_2d=2
    box_face=first(sort!(collect(keys(box.surfaces))))
    mesh=M.mesh_model_surface(box,box_face)
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    # An all-plane primitive Box consumes its actual algorithm-1 face PLC.
    local_box=GeoModel()
    add_box!(local_box,0.0,0.0,0.0,1.0,1.0,1.0;tag=1)
    local_box.meshing.algorithm_2d=1
    local_mesh=M.mesh_model_volume(local_box,1)
    @test MT.validate(local_mesh).ok && MT.ntets(local_mesh)>0
    sphere=GeoModel()
    add_sphere!(sphere,0.0,0.0,0.0,1.0;tag=1)
    sphere.meshing.algorithm_2d=2
    smesh=M.mesh_model_volume(sphere,1)
    @test MT.validate(smesh).ok && MT.ntets(smesh)>0
    explicit=GeoModel()
    add_sphere!(explicit,0.0,0.0,0.0,1.0;tag=1)
    explicit.meshing.algorithm_2d=1
    explicit_mesh=M.mesh_model_volume(explicit,1)
    @test MT.validate(explicit_mesh).ok && MT.ntets(explicit_mesh)>0
end

@testset "Mesh.Algorithm3D ordinary-volume admission" begin
    # Delaunay(1), initial-only(3) and R-tree(9) are the implemented
    # ordinary policies — 9 aliases the Delaunay fill unless Recombine3D
    # asks for the hex-combining arm.  The remaining documented selectors
    # (Frontal 4 = netgen, MMG3D 7, HXT 10) have no native kernel and keep
    # their loud rejection.
    delaunay=M.mesh_model_volume(_primitive_cylinder(6),1)
    @test MT.ntets(delaunay)>0
    m=_primitive_cylinder(6)
    m.meshing.algorithm_3d=3
    initial=M.mesh_model_volume(m,1)
    @test MT.ntets(initial)>0
    m=_primitive_cylinder(6)
    m.meshing.algorithm_3d=9
    rtree=M.mesh_model_volume(m,1)
    @test MT.ntets(rtree)>0
    for algorithm3d in (4,7,10)
        model=_primitive_cylinder(6)
        model.meshing.algorithm_3d=algorithm3d
        err=_error(()->M.mesh_model_volume(model,1))
        @test err isa ArgumentError
        @test occursin("not implemented",sprint(showerror,err))
    end
end
end
