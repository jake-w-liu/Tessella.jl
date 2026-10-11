module ExplicitEmittedProjectionChecks
using Test,Tessella
using Tessella.MeshTypes: nnodes
const M=Tessella.Model
const SOURCE=raw"""
Point(1) = {0,0,0,0.4}; Point(2) = {1,0,0,0.4};
Point(3) = {1,1,0,0.4}; Point(4) = {0,1,0,0.4};
Point(5) = {0,0,1,0.4}; Point(6) = {1,0,1,0.4};
Point(7) = {1,1,1,0.4}; Point(8) = {0,1,1,0.4};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
Curve Loop(3) = {1,10,-5,-9}; Curve Loop(4) = {2,11,-6,-10};
Curve Loop(5) = {3,12,-7,-11}; Curve Loop(6) = {4,9,-8,-12};
Plane Surface(1) = {1}; Plane Surface(2) = {2};
Plane Surface(3) = {3}; Plane Surface(4) = {4};
Plane Surface(5) = {5}; Plane Surface(6) = {6};
Surface Loop(1) = {1,2,3,4,5,6}; Volume(1) = {1};
Mesh 3;
"""
@testset "plain explicit Volume projection binds its emitted constructor" begin
    execution=mktemp() do path,io
        write(io,SOURCE);close(io);execute_geo(path;mesh_dim=3)
    end
    @test isempty(execution.model.meshing.transfinite_volumes)
    @test isempty(execution.model.meshing.extrude)
    part=geo_entity_mesh(execution,3,1;geometry=true)
    @test part isa Tessella.Mesh
    authority=execution.model.meshing.mesh_authority
    state=M._model_mesh_emitted_state(authority.mesh)
    @test authority.emitted_geometry!==nothing
    projected=model_to_mixed(execution.model,part,3,1)
    @test validate(projected).ok
    @test projected.coords[:,1:nnodes(part)]==part.coords
    @test all(isnothing,projected.entity_data.node_parametric)
    @test !projected.entity_data.native_mesh_tags
    @test M._model_mesh_emitted_state_matches(authority.mesh,state)
    @test_throws ArgumentError M._model_mesh_projection_input(
        execution.model,projected,3,1,"plain emitted geometry role")
    changed=deepcopy(part);changed.coords[1,1]=nextfloat(changed.coords[1,1])
    @test_throws ArgumentError model_to_mixed(execution.model,changed,3,1)
    @test M._model_mesh_emitted_state_matches(authority.mesh,state)
    tagged=geo_entity_mesh(execution,3,1)
    actual=model_to_mixed(execution.model,tagged,3,1)
    @test validate(actual;reject_duplicate_point_cells=false).ok
    @test actual.entity_data.node_parametric[1:nnodes(tagged)]==tagged.entity_data.node_parametric
end
end
