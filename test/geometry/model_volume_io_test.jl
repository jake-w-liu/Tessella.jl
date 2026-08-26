using Test
using Tessella
using Tessella.MeshTypes: Mesh, nnodes, ntets

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

    embed!(model,0,[10],3,1)
    embed!(model,1,[20],3,1)
    embed!(model,2,[30],3,1)
    add_physical_group!(
        model,0,[10,11,12,21,22,23];tag=41,name="embedded points")
    add_physical_group!(
        model,1,[20,31,32,33];tag=42,name="embedded curves")
    add_physical_group!(model,2,[30];tag=43,name="embedded sheet")
    add_physical_group!(model,3,[1];tag=44,name="domain")
    return model,mesh_model_volume(model,1)
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
    expected_points=Set(Int32[10,11,12,21,22,23])
    expected_curves=Set(Int32[20,31,32,33])
    @test Set(projected.entity_data.block_entities[point_block])==expected_points
    @test Set(projected.entity_data.block_entities[line_block])==expected_curves
    @test projected.entity_data.block_entities[surface_block]==
          fill(Int32(30),length(projected.blocks[surface_block].tags))
    @test projected.entity_data.block_entities[volume_block]==
          fill(Int32(1),ntets(mesh))
    @test projected.entity_data.entities[(2,30)].boundaries==Int32[31,32,33]
    @test isempty(projected.entity_data.entities[(2,30)].embedded_curves)
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
    @test nnodes(mesh)==size(projected.coords,2)

    plain_model=deepcopy(model)
    empty!(plain_model.embeds)
    plain=model_to_mixed(plain_model,mesh,3,1)
    @test [block.msh for block in plain.blocks]==[4]
    @test Set(keys(plain.entity_data.entities))==Set([(3,1)])
    @test plain.physical_names==Dict((3,44)=>"domain")

    crc=mixed_crc(projected)
    @test crc.sha==
          "d12446ff4f9c24254d02ae3938fc51b1010063399ab1a4d42b8947bd6637c1f8"
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

    nested=deepcopy(model)
    embed!(nested,0,[10],2,30)
    @test_throws ArgumentError model_to_mixed(nested,mesh,3,1)

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
