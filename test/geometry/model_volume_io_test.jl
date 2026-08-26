using Test
using Tessella
using Tessella.MeshTypes: Mesh, nnodes, ntets, tet_volume, node
using Tessella.Mesh3D: insert_steiner3, recover_segment3

function _classified_volume_fixture()
    model=GeoModel()
    add_box!(model,0,0,0,1,1,1;tag=1)

    add_point!(model,0.3,0.7,0.8;tag=10)
    add_point!(model,0.2,0.2,0.2;tag=11)
    add_point!(model,0.8,0.2,0.2;tag=12)
    add_line!(model,11,12;tag=20)

    add_point!(model,0.2,0.2,0.5;tag=21)
    add_point!(model,0.8,0.2,0.5;tag=22)
    add_point!(model,0.5,0.8,0.5;tag=23)
    add_line!(model,21,22;tag=31)
    add_line!(model,22,23;tag=32)
    add_line!(model,23,21;tag=33)
    add_curve_loop!(model,[31,32,33];tag=30)
    add_plane_surface!(model,[30];tag=30)

    add_point!(model,0.3,0.35,0.5;tag=24)
    add_point!(model,0.7,0.35,0.5;tag=25)
    add_point!(model,0.5,0.35,0.5;tag=26)
    add_line!(model,24,25;tag=34)
    embed!(model,0,[26],2,30)
    embed!(model,1,[34],2,30)

    embed!(model,0,[10],3,1)
    embed!(model,1,[20],3,1)
    embed!(model,2,[30],3,1)
    add_physical_group!(
        model,0,[10,11,12,21,22,23,24,25,26];
        tag=41,name="embedded points")
    add_physical_group!(
        model,1,[20,31,32,33,34];tag=42,name="embedded curves")
    add_physical_group!(model,2,[30];tag=43,name="embedded sheet")
    add_physical_group!(model,3,[1];tag=44,name="domain")
    return model,mesh_model_volume(model,1)
end

function _add_explicit_cube_shell!(model,offset::Int,lo::Float64,hi::Float64)
    coordinates=((lo,lo,lo),(hi,lo,lo),(hi,hi,lo),(lo,hi,lo),
                 (lo,lo,hi),(hi,lo,hi),(hi,hi,hi),(lo,hi,hi))
    for (index,coordinate) in pairs(coordinates)
        add_point!(model,coordinate...;tag=offset+index)
    end
    endpoints=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
               (1,5),(2,6),(3,7),(4,8))
    for (index,(first_point,last_point)) in pairs(endpoints)
        add_line!(model,offset+first_point,offset+last_point;tag=offset+index)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),
           (3,12,-7,-11),(4,9,-8,-12))
    signed_tag(tag)=sign(tag)*(offset+abs(tag))
    for (index,curves) in pairs(loops)
        tag=offset+index
        add_curve_loop!(model,signed_tag.(curves);tag=tag)
        add_plane_surface!(model,[tag];tag=tag)
    end
    surfaces=collect((offset+1):(offset+6))
    return add_surface_loop!(model,surfaces;tag=offset+1)
end

function _mesh_volume_value(mesh)
    return sum(tet_volume(
        node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
        node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
        for cell in 1:ntets(mesh))
end

@testset "classified native volume projection" begin
    model,mesh=_classified_volume_fixture()
    projected=model_to_mixed(model,mesh,3,1)
    @test validate(projected).ok
    @test isempty(projected.periodic_links)
    @test [block.msh for block in projected.blocks]==[15,1,2,4]

    point_block=only(findall(block->block.msh==15,projected.blocks))
    line_block=only(findall(block->block.msh==1,projected.blocks))
    surface_block=only(findall(block->block.msh==2,projected.blocks))
    volume_block=only(findall(block->block.msh==4,projected.blocks))
    expected_points=Set(Int32[10,11,12,21,22,23,24,25,26])
    expected_curves=Set(Int32[20,31,32,33,34])
    @test Set(projected.entity_data.block_entities[point_block])==expected_points
    @test Set(projected.entity_data.block_entities[line_block])==expected_curves
    @test projected.entity_data.block_entities[surface_block]==
          fill(Int32(30),length(projected.blocks[surface_block].tags))
    @test projected.entity_data.block_entities[volume_block]==
          fill(Int32(1),ntets(mesh))
    @test projected.entity_data.entities[(2,30)].boundaries==Int32[31,32,33]
    @test projected.entity_data.entities[(2,30)].embedded_curves==Int32[34]
    @test isempty(projected.entity_data.entities[(3,1)].boundaries)
    @test projected.physical_names==Dict(
        (0,41)=>"embedded points",(1,42)=>"embedded curves",
        (2,43)=>"embedded sheet",(3,44)=>"domain")

    classified_points=Dict(Int(entity[2])=>Int32(node) for (node,entity) in
        enumerate(projected.entity_data.node_entities) if entity[1]==0)
    @test Set(keys(classified_points))==Set(Int.(expected_points))
    @test all(==(Int32(41)),projected.blocks[point_block].tags)
    @test all(==(Int32(42)),projected.blocks[line_block].tags)
    @test all(==(Int32(43)),projected.blocks[surface_block].tags)
    @test all(==(Int32(44)),projected.blocks[volume_block].tags)
    @test size(projected.blocks[surface_block].nodes,2)>0
    @test any(entity->entity==(2,Int32(30)),
              projected.entity_data.node_entities)
    nested_point=only(node for (node,entity) in
        enumerate(projected.entity_data.node_entities)
        if entity==(0,Int32(26)))
    nested_line_cells=only(block for block in projected.blocks if block.msh==1)
    @test nested_point in nested_line_cells.nodes
    @test any(nested_point in projected.blocks[surface_block].nodes[:,cell]
              for cell in axes(projected.blocks[surface_block].nodes,2))
    @test nnodes(mesh)==size(projected.coords,2)

    plain_model=deepcopy(model)
    empty!(plain_model.embeds)
    plain=model_to_mixed(plain_model,mesh,3,1)
    @test [block.msh for block in plain.blocks]==[4]
    @test Set(keys(plain.entity_data.entities))==Set([(3,1)])
    @test plain.physical_names==Dict((3,44)=>"domain")

    crc=mixed_crc(projected)
    @test crc.sha==
          "e6a1a6de65b65987c543553d6456e4607b43fd3f3127294d926237888c9b5453"
    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"classified-volume-$version-$binary.msh")
            @test write_mixed_msh(
                path,projected;version=version,binary=binary)==path
            reread=read_mixed_msh(path)
            @test validate(reread).ok
            @test reread.physical_names==projected.physical_names
            if version==4.1
                @test mixed_crc(reread)==crc
                @test haskey(reread.entity_data.entities,(2,30))
                @test haskey(reread.entity_data.entities,(3,1))
                @test reread.entity_data.entities[(2,30)].embedded_curves==
                      Int32[34]
                @test isempty(reread.entity_data.entities[(3,1)].boundaries)
            else
                @test reread.entity_data===nothing
                ownership=Dict(block.msh=>Set(reread.elementary_entities[index])
                    for (index,block) in pairs(reread.blocks))
                @test ownership[15]==expected_points
                @test ownership[1]==expected_curves
                @test ownership[2]==Set(Int32[30])
                @test ownership[4]==Set(Int32[1])
            end
        end
    end

    @test_throws ArgumentError model_to_mixed(model,mesh,1,1)
    @test_throws ArgumentError model_to_mixed(model,mesh,true,1)
    @test_throws ArgumentError model_to_mixed(model,mesh,3,99)
    with_triangles=Mesh(
        mesh.coords;tris=reshape(mesh.tets[1:3,1],3,1),tets=mesh.tets)
    @test_throws ArgumentError model_to_mixed(model,with_triangles,3,1)
    tagged=Mesh(mesh.coords;tets=mesh.tets,tet_tag=fill(Int32(1),ntets(mesh)))
    @test_throws ArgumentError model_to_mixed(model,tagged,3,1)
    empty=Mesh(zeros(3,0))
    @test_throws ArgumentError model_to_mixed(model,empty,3,1)

    unrelated=Mesh(
        [2.0 3.0 2.0 2.0; 2.0 2.0 3.0 2.0; 2.0 2.0 2.0 3.0];
        tets=reshape(Int32[1,2,3,4],4,1))
    @test validate(unrelated).ok
    @test_throws ArgumentError model_to_mixed(plain_model,unrelated,3,1)

    explicit_shell=deepcopy(model)
    explicit_shell.volumes[1]=[30]
    @test_throws ArgumentError model_to_mixed(explicit_shell,mesh,3,1)

    off_sheet=deepcopy(model)
    add_point!(off_sheet,0.5,0.35,0.6;tag=27)
    embed!(off_sheet,0,[27],2,30)
    off_sheet_mesh,_=insert_steiner3(mesh,(0.5,0.35,0.6))
    @test_throws ArgumentError model_to_mixed(off_sheet,off_sheet_mesh,3,1)

    off_sheet_curve=deepcopy(model)
    add_point!(off_sheet_curve,0.3,0.35,0.6;tag=28)
    add_point!(off_sheet_curve,0.7,0.35,0.6;tag=29)
    add_line!(off_sheet_curve,28,29;tag=35)
    embed!(off_sheet_curve,1,[35],2,30)
    off_sheet_curve_mesh=recover_segment3(
        mesh,(0.3,0.35,0.6),(0.7,0.35,0.6))
    @test_throws ArgumentError model_to_mixed(
        off_sheet_curve,off_sheet_curve_mesh,3,1)

    bounding_embed=deepcopy(model)
    embed!(bounding_embed,1,[31],2,30)
    @test_throws ArgumentError mesh_model_volume(bounding_embed,1)

    overlapping_curve=deepcopy(model)
    add_line!(overlapping_curve,11,12;tag=40)
    embed!(overlapping_curve,1,[40],3,1)
    @test_throws ArgumentError model_to_mixed(overlapping_curve,mesh,3,1)

    overlapping_surface=deepcopy(model)
    add_plane_surface!(overlapping_surface,[30];tag=31)
    embed!(overlapping_surface,2,[31],3,1)
    @test_throws ArgumentError model_to_mixed(overlapping_surface,mesh,3,1)

    periodic=deepcopy(model)
    add_point!(periodic,0.2,0.3,0.2;tag=13)
    add_point!(periodic,0.8,0.3,0.2;tag=14)
    add_line!(periodic,13,14;tag=21)
    translate=(1.0,0.0,0.0,0.0,
               0.0,1.0,0.0,0.1,
               0.0,0.0,1.0,0.0,
               0.0,0.0,0.0,1.0)
    set_periodic!(periodic,1,[21],[20],translate)
    @test_throws ArgumentError model_to_mixed(periodic,mesh,3,1)
end

@testset "explicit modeled volume shells" begin
    model=GeoModel()
    @test _add_explicit_cube_shell!(model,0,0.0,1.0)==1
    @test add_volume!(model,[1];tag=1)==1
    @test Tessella.Model.model_entity(model,3,1)==[1]
    add_physical_group!(model,0,collect(1:8);tag=61,name="corners")
    add_physical_group!(model,1,collect(1:12);tag=62,name="edges")
    add_physical_group!(model,2,collect(1:6);tag=63,name="boundary")
    add_physical_group!(model,3,[1];tag=64,name="domain")

    mesh=mesh_model_volume(model,1)
    @test validate(mesh).ok
    @test ntets(mesh)>0
    @test _mesh_volume_value(mesh)≈1.0 atol=1e-12
    projected=model_to_mixed(model,mesh,3,1)
    @test validate(projected).ok
    @test [block.msh for block in projected.blocks]==[15,1,2,4]
    point_block=only(findall(block->block.msh==15,projected.blocks))
    curve_block=only(findall(block->block.msh==1,projected.blocks))
    surface_block=only(findall(block->block.msh==2,projected.blocks))
    volume_block=only(findall(block->block.msh==4,projected.blocks))
    @test Set(projected.entity_data.block_entities[point_block])==Set(Int32.(1:8))
    @test Set(projected.entity_data.block_entities[curve_block])==Set(Int32.(1:12))
    @test Set(projected.entity_data.block_entities[surface_block])==Set(Int32.(1:6))
    @test projected.entity_data.block_entities[volume_block]==
          fill(Int32(1),ntets(mesh))
    @test projected.entity_data.entities[(3,1)].boundaries==Int32.(1:6)
    @test projected.entity_data.entities[(2,3)].boundaries==Int32[1,10,-5,-9]
    @test projected.physical_names==Dict(
        (0,61)=>"corners",(1,62)=>"edges",(2,63)=>"boundary",(3,64)=>"domain")
    @test all(==(Int32(61)),projected.blocks[point_block].tags)
    @test all(==(Int32(62)),projected.blocks[curve_block].tags)
    @test all(==(Int32(63)),projected.blocks[surface_block].tags)
    @test all(==(Int32(64)),projected.blocks[volume_block].tags)
    crc=mixed_crc(projected)
    @test crc.sha==
          "9bce88e319c67236317df64b876739a62f80982ed86eca028bd1e7bda022bcb6"

    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"explicit-shell-$version-$binary.msh")
            @test write_mixed_msh(path,projected;version=version,binary=binary)==path
            reread=read_mixed_msh(path)
            @test validate(reread).ok
            @test reread.physical_names==projected.physical_names
            if version==4.1
                @test mixed_crc(reread)==crc
                @test reread.entity_data.entities[(3,1)].boundaries==Int32.(1:6)
            else
                @test reread.entity_data===nothing
                ownership=Dict(block.msh=>Set(reread.elementary_entities[index])
                    for (index,block) in pairs(reread.blocks))
                @test ownership[15]==Set(Int32.(1:8))
                @test ownership[1]==Set(Int32.(1:12))
                @test ownership[2]==Set(Int32.(1:6))
                @test ownership[4]==Set(Int32[1])
            end
        end
    end

    hollow=GeoModel()
    @test _add_explicit_cube_shell!(hollow,0,0.0,2.0)==1
    @test _add_explicit_cube_shell!(hollow,100,0.5,1.5)==101
    @test add_volume!(hollow,[1,101];tag=1)==1
    hollow_mesh=mesh_model_volume(hollow,1)
    @test validate(hollow_mesh).ok
    @test _mesh_volume_value(hollow_mesh)≈7.0 atol=1e-12
    hollow_projected=model_to_mixed(hollow,hollow_mesh,3,1)
    @test validate(hollow_projected).ok
    @test hollow_projected.entity_data.entities[(3,1)].boundaries==
          Int32[1,2,3,4,5,6,-101,-102,-103,-104,-105,-106]
    @test mixed_crc(hollow_projected).sha==
          "83721952195b78f4d18b9e5ff862ef7629b9fbe6b642f33d6649d85e83d0c8b2"

    signed=GeoModel()
    _add_explicit_cube_shell!(signed,0,0.0,1.0)
    @test add_surface_loop!(signed,[-1,2,3,4,5,6];tag=20)==20
    add_point!(signed,0.25,0.5,1.0;tag=20)
    add_point!(signed,0.75,0.5,1.0;tag=21)
    add_point!(signed,0.5,0.25,1.0;tag=22)
    add_line!(signed,20,21;tag=20)
    embed!(signed,0,[22],2,2)
    embed!(signed,1,[20],2,2)
    @test add_volume!(signed,[20];tag=1)==1
    signed_mesh=mesh_model_volume(signed,1)
    signed_projected=model_to_mixed(signed,signed_mesh,3,1)
    @test signed_projected.entity_data.entities[(3,1)].boundaries==
          Int32[-1,2,3,4,5,6]
    @test signed_projected.entity_data.entities[(2,2)].embedded_curves==Int32[20]

    @test_throws ArgumentError add_surface_loop!(model,[1];tag=20)
    @test_throws ArgumentError add_surface_loop!(model,[1,1];tag=20)
    @test_throws ArgumentError add_surface_loop!(model,[1,2,3,4,5,99];tag=20)
    @test_throws ArgumentError add_volume!(model,[99];tag=2)
    @test_throws ArgumentError add_volume!(model,[0];tag=2)
    @test_throws ArgumentError add_volume!(model,[1,1];tag=2)
    @test_throws ArgumentError add_volume!(model,[1];tag=1)
    @test_throws ArgumentError add_surface_loop!(model,[true];tag=20)
    @test_throws ArgumentError add_volume!(model,[true];tag=2)

    duplicate_shell=deepcopy(model)
    @test add_surface_loop!(duplicate_shell,collect(1:6);tag=20)==20
    @test add_surface_loop!(duplicate_shell,collect(1:6))==21
    @test add_plane_surface!(duplicate_shell,[1])==7
    @test_throws ArgumentError add_volume!(duplicate_shell,[1,20];tag=2)

    disconnected=GeoModel()
    _add_explicit_cube_shell!(disconnected,0,0.0,1.0)
    _add_explicit_cube_shell!(disconnected,100,2.0,3.0)
    @test_throws ArgumentError add_surface_loop!(
        disconnected,vcat(collect(1:6),collect(101:106));tag=200)

    outside_cavity=GeoModel()
    _add_explicit_cube_shell!(outside_cavity,0,0.0,2.0)
    _add_explicit_cube_shell!(outside_cavity,100,3.0,3.5)
    add_volume!(outside_cavity,[1,101];tag=1)
    @test_throws ArgumentError mesh_model_volume(outside_cavity,1)
    @test_throws ArgumentError model_to_mixed(outside_cavity,mesh,3,1)

    overlapping_cavities=GeoModel()
    _add_explicit_cube_shell!(overlapping_cavities,0,0.0,4.0)
    _add_explicit_cube_shell!(overlapping_cavities,100,0.5,2.5)
    _add_explicit_cube_shell!(overlapping_cavities,200,1.5,3.5)
    add_volume!(overlapping_cavities,[1,101,201];tag=1)
    @test_throws ArgumentError mesh_model_volume(overlapping_cavities,1)

    nonplanar=deepcopy(model)
    nonplanar.points[8]=(0.0,1.0,1.1)
    @test_throws ArgumentError mesh_model_volume(nonplanar,1)

    damaged=Mesh(mesh.coords;tets=mesh.tets[:,1:(end-1)])
    @test_throws ArgumentError model_to_mixed(model,damaged,3,1)
    missing_loop=deepcopy(model)
    missing_loop.volumes[1]=[999]
    @test_throws ArgumentError model_to_mixed(missing_loop,mesh,3,1)
end
