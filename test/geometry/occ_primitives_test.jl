# Materialized Cylinder/Sphere/Cone boundaries — Gmsh 4.15.2 OCC layout parity.
#
# `add_cylinder!`/`add_sphere!`/`add_cone!` build real boundary representations:
# rim/pole Points, closed-circle and degenerate edges, Cylinder/Sphere/Cone and
# Plane faces, and a Surface Loop — while retaining the compact encoding that
# drives native volume tessellation, exactly the `add_box!` dual-registration
# pattern.

using Test
using Tessella
using Tessella.Model: model_boundary, model_entities, model_entity_type,
                      model_value, model_derivative, model_curvature,
                      model_parametrization_bounds, model_bounding_box,
                      model_set_tag!, transform_entities!, remove_entities!,
                      duplicate_entities!, mesh_model_volume,
                      _affine_translation, _affine_rotation, _affine_dilation,
                      _materialized_curved_consistent
using Tessella.MeshTypes: ntets, validate

@testset "OCC cylinder materialization" begin
    m=GeoModel()
    @test add_cylinder!(m,0,0,0,0,0,2,1)==1
    # Gmsh 4.15.2 `addCylinder(0,0,0,0,0,2,1)`: 2 rim points (top then bottom),
    # top circle, seam line, bottom circle; lateral [-top,-seam,bottom,seam];
    # caps Plane{top}, Plane{bottom}; shell [lateral,+top,-bottom].
    @test sort!(collect(keys(m.points)))==[1,2]
    @test m.points[1]==(1.0,0.0,2.0) && m.points[2]==(1.0,0.0,0.0)
    @test sort!(collect(keys(m.curves)))==[1,2,3]
    @test m.curves[1]==(1,1) && m.curves[2]==(2,1) && m.curves[3]==(2,2)
    @test sort!(collect(keys(m.surfaces)))==[1,2,3]
    @test sort!(collect(keys(m.surface_loops)))==[1]
    @test m.volumes[1]==[1]
    @test m.surface_loops[1]==[1,2,-3]
    @test model_entity_type(m,1,1)=="Circle"
    @test model_entity_type(m,1,2)=="Line"
    @test model_entity_type(m,1,3)=="Circle"
    @test model_entity_type(m,2,1)=="Cylinder"
    @test model_entity_type(m,2,2)=="Plane"
    @test model_entity_type(m,2,3)=="Plane"
    # Signed boundaries through the public query API.
    @test last.(model_boundary(m,[(3,1)],true,true,false))==[1,2,-3]
    @test last.(model_boundary(m,[(2,1)],false,true,false))==[-1,-2,3,2]
    @test last.(model_boundary(m,[(2,2)],false,true,false))==[1]
    @test last.(model_boundary(m,[(2,3)],false,true,false))==[3]
    @test last.(model_boundary(m,[(1,1)],false,false,false))==[1,1]
    @test last.(model_boundary(m,[(1,2)],false,false,false))==[2,1]
    @test last.(model_boundary(m,[(1,3)],false,false,false))==[2,2]
    # OCC parametrization: circles span [0,2π], the seam is arc-length.
    @test model_parametrization_bounds(m,1,1)==([0.0],[2π])
    @test model_parametrization_bounds(m,1,2)==([0.0],[2.0])
    @test model_value(m,1,1,[0.0])≈[1.0,0.0,2.0]
    @test model_value(m,1,1,[π])≈[-1.0,0.0,2.0] atol=1e-15
    @test model_value(m,1,3,[π/2])≈[0.0,1.0,0.0] atol=1e-15
    @test model_value(m,1,2,[1.0])≈[1.0,0.0,1.0] atol=1e-15
    @test model_derivative(m,1,2,[0.5])≈[0.0,0.0,1.0] atol=1e-15
    @test model_curvature(m,1,1,[0.3])≈[1.0]
    @test model_curvature(m,1,2,[0.5])≈[0.0] atol=1e-15
    # Analytic bounds from the encoding and the OCC geometry records.
    @test model_bounding_box(m,3,1)==
          (-1.0,-1.0,0.0,1.0,1.0,2.0)
    @test model_bounding_box(m,2,1)==
          (-1.0,-1.0,0.0,1.0,1.0,2.0)
    @test model_bounding_box(m,1,1)==
          (-1.0,-1.0,2.0,1.0,1.0,2.0)
    # The compact encoding stays valid and drives tessellation.
    @test _materialized_curved_consistent(m,1)
    mesh=mesh_model_volume(m,1)
    @test ntets(mesh)>0 && validate(mesh).ok
end

@testset "OCC sphere materialization" begin
    m=GeoModel()
    @test add_sphere!(m,0,0,0,1)==1
    # `addSphere(0,0,0,1)`: north then south pole; degenerate north edge,
    # meridian circle, degenerate south edge; one Sphere face; shell [face].
    @test sort!(collect(keys(m.points)))==[1,2]
    @test m.points[1]==(0.0,0.0,1.0) && m.points[2]==(0.0,0.0,-1.0)
    @test sort!(collect(keys(m.curves)))==[1,2,3]
    @test m.curves[1]==(1,1) && m.curves[2]==(2,1) && m.curves[3]==(2,2)
    @test sort!(collect(keys(m.surfaces)))==[1]
    @test m.volumes[1]==[1]
    @test m.surface_loops[1]==[1]
    @test model_entity_type(m,2,1)=="Sphere"
    @test model_entity_type(m,1,1)=="Unknown"
    @test model_entity_type(m,1,2)=="Circle"
    @test model_entity_type(m,1,3)=="Unknown"
    @test last.(model_boundary(m,[(3,1)],true,true,false))==[1]
    @test last.(model_boundary(m,[(2,1)],false,true,false))==[-1,-2,3,2]
    # Degenerate edges evaluate to their pole with zero derivatives.
    @test model_value(m,1,1,[0.0])≈[0.0,0.0,1.0]
    @test model_value(m,1,1,[2.0])≈[0.0,0.0,1.0]
    @test model_derivative(m,1,1,[1.0])==[0.0,0.0,0.0]
    @test model_curvature(m,1,1,[0.0])≈[0.0]
    # Meridian: OCC trims to [3π/2,5π/2], sweeping through +x̂.
    t0,t1=model_parametrization_bounds(m,1,2)
    @test t0==[1.5π] && t1==[2.5π]
    @test model_value(m,1,2,[2π])≈[1.0,0.0,0.0] atol=1e-15
    @test model_value(m,1,2,[1.5π])≈[0.0,0.0,-1.0] atol=1e-15
    @test model_curvature(m,1,2,[2.0])≈[1.0]
    # The single meridian wire must not bound the face — the surface bbox
    # comes from the spherical geometry record.
    @test model_bounding_box(m,2,1)==
          (-1.0,-1.0,-1.0,1.0,1.0,1.0)
    @test model_bounding_box(m,3,1)==
          (-1.0,-1.0,-1.0,1.0,1.0,1.0)
    @test _materialized_curved_consistent(m,1)
    mesh=mesh_model_volume(m,1)
    @test ntets(mesh)>0 && validate(mesh).ok
end

@testset "OCC cone materialization" begin
    # Truncated cone mirrors the cylinder layout.
    m=GeoModel()
    @test add_cone!(m,0,0,0,0,0,2,2,1)==1
    @test sort!(collect(keys(m.points)))==[1,2]
    @test m.points[1]==(1.0,0.0,2.0) && m.points[2]==(2.0,0.0,0.0)
    @test sort!(collect(keys(m.curves)))==[1,2,3]
    @test sort!(collect(keys(m.surfaces)))==[1,2,3]
    @test m.surface_loops[1]==[1,2,-3]
    @test model_entity_type(m,2,1)=="Cone"
    @test model_entity_type(m,1,1)=="Circle"
    @test model_entity_type(m,1,3)=="Circle"
    @test _materialized_curved_consistent(m,1)
    @test model_bounding_box(m,3,1)==
          (-2.0,-2.0,0.0,2.0,2.0,2.0)
    @test ntets(mesh_model_volume(m,1))>0

    # r2 == 0: degenerate apex edge replaces the top circle; no top cap.
    m=GeoModel()
    @test add_cone!(m,0,0,0,0,0,2,2,0)==1
    @test m.points[1]==(0.0,0.0,2.0) && m.points[2]==(2.0,0.0,0.0)
    @test sort!(collect(keys(m.curves)))==[1,2,3]
    @test sort!(collect(keys(m.surfaces)))==[1,2]
    @test m.surface_loops[1]==[1,-2]
    @test model_entity_type(m,1,1)=="Unknown"
    @test model_entity_type(m,1,3)=="Circle"
    @test last.(model_boundary(m,[(2,1)],false,true,false))==[-1,-2,3,2]
    @test model_value(m,1,1,[0.0])≈[0.0,0.0,2.0]
    @test _materialized_curved_consistent(m,1)
    @test ntets(mesh_model_volume(m,1))>0

    # r1 == 0: apex at the base; only the top cap.
    m=GeoModel()
    @test add_cone!(m,0,0,0,0,0,2,0,2)==1
    @test m.points[1]==(2.0,0.0,2.0) && m.points[2]==(0.0,0.0,0.0)
    @test sort!(collect(keys(m.surfaces)))==[1,2]
    @test m.surface_loops[1]==[1,2]
    @test model_entity_type(m,1,1)=="Circle"
    @test model_entity_type(m,1,3)=="Unknown"
    @test model_value(m,1,3,[0.0])≈[0.0,0.0,0.0]
    @test _materialized_curved_consistent(m,1)
    @test ntets(mesh_model_volume(m,1))>0

    # Both radii zero is still rejected before any entity exists.
    m=GeoModel()
    @test_throws ArgumentError add_cone!(m,0,0,0,0,0,2,0,0)
    @test isempty(m.points) && isempty(m.curves) && isempty(m.volumes)
end

@testset "materialized primitive consistency and lifecycle" begin
    m=GeoModel()
    add_cylinder!(m,0,0,0,0,0,2,1)
    # Whole-volume affine updates the encoding and every OCC record.
    transform_entities!(m,_affine_translation((1.0,0.0,0.0),"test"),[(3,1)])
    @test m.cylinders[1].center==(1.0,0.0,0.0)
    @test m.curve_geometry[1].center==(1.0,0.0,2.0)
    @test m.surface_geometry[1].center==(1.0,0.0,0.0)
    @test _materialized_curved_consistent(m,1)
    @test ntets(mesh_model_volume(m,1))>0
    # Independent vertex moves drop the stale encoding instead of leaving
    # inconsistent dual state; the explicit shell remains.
    transform_entities!(m,_affine_translation((0.5,0.0,0.0),"test"),[(0,2)])
    @test !haskey(m.cylinders,1)
    @test !_materialized_curved_consistent(m,1)
    # Curved shells cannot ride the planar explicit-shell mesher.
    @test_throws ArgumentError mesh_model_volume(m,1)

    # Retagging migrates encodings and OCC state with the entity.
    m=GeoModel()
    add_sphere!(m,0,0,0,1)
    model_set_tag!(m,3,1,7)
    @test haskey(m.spheres,7) && !haskey(m.spheres,1)
    @test _materialized_curved_consistent(m,7)
    model_set_tag!(m,1,2,9)
    @test m.curves[9]==(2,1)
    @test haskey(m.curve_geometry,9) && !haskey(m.curve_geometry,2)
    @test _materialized_curved_consistent(m,7)

    # Non-recursive volume removal keeps the shell but drops the encoding.
    m=GeoModel()
    add_cone!(m,0,0,0,0,0,2,2,1)
    remove_entities!(m,[(3,1)])
    @test !haskey(m.volumes,1) && !haskey(m.cones,1)
    @test sort!(collect(keys(m.surfaces)))==[1,2,3]

    # OCC-copy duplication shares boundary copies across the whole call —
    # each source entity maps to exactly one copy (Gmsh occ.copy layout).
    m=GeoModel()
    add_cylinder!(m,0,0,0,0,0,2,1)
    out=duplicate_entities!(m,[(3,1)])
    @test out==[(3,4)]
    @test length(m.points)==4 && length(m.curves)==6 &&
          length(m.surfaces)==6 && length(m.volumes)==2
    @test haskey(m.cylinders,4)
    @test _materialized_curved_consistent(m,4)
    @test ntets(mesh_model_volume(m,4))>0
end

@testset ".geo primitive materialization" begin
    r=mktemp() do path,io
        write(io,raw"""
            Cylinder(1) = {0,0,0,0,0,2,1};
            Sphere(2) = {5,0,0,1};
            Cone(3) = {0,5,0,0,0,1,1,0.5};
            """)
        close(io)
        execute_geo(path)
    end
    model=r.model
    @test sort!(collect(keys(model.volumes)))==[1,2,3]
    @test length(model.points)==6 && length(model.curves)==9 &&
          length(model.surfaces)==7
    @test last.(model_boundary(model,[(3,1)],true,true,false))==[1,2,-3]
    @test last.(model_boundary(model,[(3,2)],true,true,false))==[4]
    @test last.(model_boundary(model,[(3,3)],true,true,false))==[5,6,-7]
    @test model_entity_type(model,2,4)=="Sphere"
    @test model_entity_type(model,2,5)=="Cone"
    for v in 1:3
        @test _materialized_curved_consistent(model,v)
        @test ntets(mesh_model_volume(model,v))>0
    end
end
