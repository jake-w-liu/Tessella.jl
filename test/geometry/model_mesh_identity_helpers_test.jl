using Test
using Tessella

@testset "Raw discrete curve meshes retain stored node classification" begin
    model=GeoModel()
    add_point!(model,0.,0,0;tag=1)
    add_point!(model,1.,0,0;tag=2)
    Tessella.Model.add_discrete_entity!(model,1,1,[(0,1),(0,2)])
    # A stored mesh can enumerate its vertices independently of curve order.
    Tessella.Model.add_discrete_nodes!(model,1,1,[11,12,13],
                                      [0.5,0,0,1,0,0,0,0,0])
    Tessella.Model.add_discrete_elements!(model,1,1,[1],[[101,102]],
                                         [[13,11,11,12]])
    part=only(Tessella.Model._model_curve_mesh_parts(model,"test"))[3]
    owners=Tessella.Model._model_mesh_part_node_entities(model,1,1,part,"test")
    @test owners==fill((1,Int32(1)),3)
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[(1,1,part)]
    append!(parts,Tessella.Model._model_point_mesh_parts(model,"test"))
    part_owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    merged,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=part_owners)
    @test size(merged.coords,2)==5
    @test length(merged.segs)==4
    @test count(==((1,1)),values(node_owners))==3
    @test count(==((0,1)),values(node_owners))==1
    @test count(==((0,2)),values(node_owners))==1
    # A graded discrete curve replaces its raw records and uses the actual
    # model vertices at the emitted chain's two boundary positions.
    model.curve_params[1]=[0.,0.5,1.]
    remeshed=Mesh(Float64[0 0.5 1;0 0 0;0 0 0];segs=Int32[1 2;2 3])
    @test Tessella.Model._model_mesh_part_node_entities(model,1,1,remeshed,"test")==
        [(0,Int32(1)),(1,Int32(1)),(0,Int32(2))]
    # Closed native curve parts omit the duplicated endpoint. Their last
    # coordinate column remains a curve vertex, rather than a second point.
    model.curves[2]=(1,1)
    closed=Mesh(Float64[0 1 1;0 0 1;0 0 0];
                segs=Int32[1 2 3;2 3 1])
    closed_owners=Tessella.Model._model_mesh_part_node_entities(model,1,2,closed,"test")
    @test closed_owners[1]==(0,Int32(1))
    @test closed_owners[2:3]==[(1,Int32(2)),(1,Int32(2))]
end

@testset "Native point attachments replace the generated geometry vertex" begin
    model=GeoModel()
    add_point!(model,0.,0,0;tag=1)
    Tessella.Model.add_discrete_nodes!(model,0,1,[11],[2.,0,0])
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
        Tessella.Model._model_point_mesh_parts(model,"test")...]
    @test length(parts)==1
    @test only(parts)[3].coords==reshape([2.,0,0],3,1)
    owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    mesh,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=owners)
    @test size(mesh.coords,2)==1
    @test node_owners[1]==(0,1)
    cells=Tessella.GeoExec._geo_homology_cells(model,mesh,parts,"test";
        node_owner=node_owners,part_node_entities=owners)
    @test length(cells[(0,1)])==1
    @test only(only(cells[(0,1)])[2])==1
    Tessella.Model.add_physical_group!(model,0,[1];tag=7)
    Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0])
    Tessella.Model.compute_homology!(model,cells)
    @test Set(values(model.physical_names))==Set(["H_0{7}1"])
    # The first stored vertex is the endpoint reused by incident curves.
    # Additional point vertices are retained when that endpoint is on geometry.
    on_geometry=GeoModel()
    add_point!(on_geometry,0.,0,0;tag=1)
    add_point!(on_geometry,1.,0,0;tag=2)
    add_line!(on_geometry,1,2;tag=1)
    Tessella.Model.add_discrete_nodes!(on_geometry,0,1,[11,12],[0.,0,0,2.,0,0])
    @test sum(part->size(part[3].coords,2),
        Tessella.Model._model_point_mesh_parts(on_geometry,"test"))==3
    add_point!(model,1.,0,0;tag=3)
    add_line!(model,1,3;tag=1)
    before=copy(model.meshing.attached[(0,1)].node_coords)
    @test_throws r"Point\[1\].*off its geometry.*Curve\[1\]" Tessella.Model._model_point_mesh_parts(model,"test")
    @test model.meshing.attached[(0,1)].node_coords==before
    @test isempty(model.curve_params)
    volume=GeoModel()
    add_box!(volume,0.,0,0,1,1,1)
    Tessella.Model.add_discrete_nodes!(volume,0,1,[11],[2.,0,1])
    @test_throws r"Point\[1\].*off its geometry.*Curve\[1\]" mesh_model_surface(volume,1)
    @test_throws r"Point\[1\].*off its geometry.*Curve\[1\]" mesh_model_volume(volume,1)
    @test isempty(volume.curve_params)
end

@testset "Raw coincident node tags fail before topology is welded" begin
    # Gmsh 4.15.2 MeshOnlyEmpty retains these six curve vertices and the
    # two separate model points; H_0 on the curve has rank two. The native
    # coordinate merger must not collapse the disconnected chains to one.
    for native in (false,true),signed_zero in (false,true)
        model=GeoModel()
        add_point!(model,0.,0,0;tag=1)
        add_point!(model,1.,0,0;tag=2)
        native ? add_line!(model,1,2;tag=1) :
            Tessella.Model.add_discrete_entity!(model,1,1,[(0,1),(0,2)])
        x0=signed_zero ? -0.0 : 0.0
        Tessella.Model.add_discrete_nodes!(model,1,1,collect(11:16),
            [0.,0,0,0.5,0,0,1,0,0,x0,0,0,0.5,0,0,1,0,0])
        Tessella.Model.add_discrete_elements!(model,1,1,[1],[[101,102,103,104]],
            [[11,12,12,13,14,15,15,16]])
        record=get(model.discrete,(1,1),get(model.meshing.attached,(1,1),nothing))
        before_coords=copy(record.node_coords)
        before_nodes=deepcopy(record.element_nodes)
        @test_throws r"raw Curve\[1\].*tags 11 and 14.*identical coordinates" Tessella.Model._model_curve_mesh_parts(
            model,"test")
        @test isequal(record.node_coords,before_coords)
        @test record.element_nodes==before_nodes
        @test !haskey(model.curve_params,1)
    end
    # Point records pass through the same coordinate merger. Signed zero
    # is one position; distinct stored tags at it still have distinct identity.
    for signed_zero in (false,true)
        model=GeoModel()
        Tessella.Model.add_discrete_entity!(model,0,1)
        Tessella.Model.add_discrete_nodes!(model,0,1,[11,12],
            [0.,0,0,signed_zero ? -0.0 : 0.0,0,0])
        Tessella.Model.add_discrete_elements!(model,0,1,[15],[[101,102]],[[11,12]])
        @test_throws r"raw Point\[1\].*tags 11 and 12.*identical coordinates" Tessella.Model._model_point_mesh_parts(
            model,"test")
        @test model.discrete[(0,1)].node_tags==Int32[11,12]
    end
end

@testset "Native curve attachments retain raw connectivity ownership" begin
    model=GeoModel()
    add_point!(model,0.,0,0;tag=1)
    add_point!(model,1.,0,0;tag=2)
    add_line!(model,1,2;tag=1)
    Tessella.Model.add_discrete_nodes!(model,1,1,[11,12,13],
                                      [0.5,0,0,1,0,0,0,0,0])
    Tessella.Model.add_discrete_elements!(model,1,1,[1],[[101,102]],
                                         [[13,11,11,12]])
    @test !haskey(model.curve_params,1)
    part=only(Tessella.Model._model_curve_mesh_parts(model,"test"))[3]
    owners=Tessella.Model._model_mesh_part_node_entities(model,1,1,part,"test")
    @test owners==fill((1,Int32(1)),3)
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[(1,1,part)]
    append!(parts,Tessella.Model._model_point_mesh_parts(model,"test"))
    part_owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    merged,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=part_owners)
    @test size(merged.coords,2)==5
    @test Set(values(node_owners))==Set([(1,1),(0,1),(0,2)])
end
