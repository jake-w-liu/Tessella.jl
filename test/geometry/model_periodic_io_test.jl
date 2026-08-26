using Test
using Tessella
using Tessella.MeshTypes: Mesh, nnodes, ntris

function _periodic_projection_fixture()
    model=GeoModel()
    for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),
                                  (1.0,1.0),(0.0,1.0)))
        add_point!(model,x,y,0;tag=tag,mesh_size=0.5)
    end
    for (tag,(first_point,last_point)) in
        enumerate(((1,2),(2,3),(3,4),(4,1)))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    add_physical_group!(model,0,[1,2,3,4];tag=11,name="corners")
    add_physical_group!(model,1,[1,2,3,4];tag=12,name="walls")
    add_physical_group!(model,2,[1];tag=13,name="domain")
    add_physical_group!(model,2,[1];tag=14,name="domain_alias")
    affine=(1.0,0.0,0.0,1.0,
            0.0,1.0,0.0,0.0,
            0.0,0.0,1.0,0.0,
            0.0,0.0,0.0,1.0)
    set_periodic!(model,1,[2],[4],affine)
    return model,mesh_model_surface(model,1)
end

function _holed_projection_fixture()
    model=GeoModel()
    coordinates=((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0),
                 (0.25,0.25),(0.75,0.25),(0.75,0.75),(0.25,0.75))
    for (tag,(x,y)) in enumerate(coordinates)
        add_point!(model,x,y,0;tag=tag,mesh_size=0.5)
    end
    endpoints=((1,2),(2,3),(3,4),(4,1),
               (5,6),(6,7),(7,8),(8,5))
    for (tag,(first_point,last_point)) in enumerate(endpoints)
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_curve_loop!(model,[5,6,7,8];tag=2)
    add_plane_surface!(model,[1,2];tag=1)
    add_physical_group!(model,1,[5,6,7,8];tag=22,name="hole")
    add_physical_group!(model,2,[1];tag=23,name="domain")
    return model,mesh_model_surface(model,1)
end

function _embedded_projection_fixture()
    model=GeoModel()
    coordinates=((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0),
                 (0.25,0.5),(0.75,0.5),(0.5,0.5))
    for (tag,(x,y)) in enumerate(coordinates)
        add_point!(model,x,y,0;tag=tag,mesh_size=0.5)
    end
    for (tag,(first_point,last_point)) in enumerate(
            ((1,2),(2,3),(3,4),(4,1),(5,6)))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    embed!(model,0,[7],2,1)
    embed!(model,1,[5],2,1)
    add_physical_group!(model,0,[5,6,7];tag=31,name="embedded points")
    add_physical_group!(model,1,[5];tag=32,name="embedded line")
    add_physical_group!(model,2,[1];tag=33,name="domain")
    return model,mesh_model_surface(model,1)
end

@testset "native model periodic projection to mixed MSH metadata" begin
    model,mesh=_periodic_projection_fixture()
    projected=model_to_mixed(model,mesh,1)
    @test validate(projected).ok
    projected_crc=mixed_crc(projected)
    @test projected_crc.sha==
          "d9aa0af0ed218f321adea7b7276312583ee31f4747e3771ca410b87be3b628b7"
    @test [block.msh for block in projected.blocks]==[15,1,2]
    @test [length(block.tags) for block in projected.blocks]==[4,16,22]
    @test all(==(11),projected.blocks[1].tags)
    @test all(==(12),projected.blocks[2].tags)
    @test all(==(13),projected.blocks[3].tags)
    @test projected.physical_names==Dict(
        (0,11)=>"corners",(1,12)=>"walls",(2,13)=>"domain",
        (2,14)=>"domain_alias")

    data=projected.entity_data
    @test data!==nothing
    @test Set(keys(data.entities))==
          Set([(0,tag) for tag in 1:4] ∪
              [(1,tag) for tag in 1:4] ∪ [(2,1)])
    @test data.entities[(1,2)].boundaries==Int32[2,-3]
    @test data.entities[(2,1)].boundaries==Int32[1,2,3,4]
    @test data.entities[(2,1)].physical_tags==Int32[13,14]
    @test count(entity->entity[1]==0,data.node_entities)==4
    @test count(entity->entity[1]==1,data.node_entities)==12
    @test count(entity->entity[1]==2,data.node_entities)==nnodes(mesh)-16
    @test projected.elementary_entities==data.block_entities
    @test count(==(Int32(2)),data.block_entities[2])==4
    @test all(==(Int32(1)),data.block_entities[3])

    @test length(projected.periodic_links)==3
    @test sort([(Int(link.slave_entity),Int(link.master_entity))
                for link in projected.periodic_links if link.dim==0])==
          [(2,1),(3,4)]
    curve_link=only(filter(link->link.dim==1,projected.periodic_links))
    mapping=model_periodic_nodes(model,mesh,1,2)
    @test curve_link.slave_entity==2
    @test curve_link.master_entity==4
    @test curve_link.slave_nodes==mapping.slave_nodes
    @test curve_link.master_nodes==mapping.master_nodes
    @test curve_link.affine==mapping.affine

    format_crcs=Dict(2.2=>Set{String}(),4.1=>Set{String}())
    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"projected-$version-$binary.msh")
            @test write_mixed_msh(
                path,projected;version=version,binary=binary)==path
            reread=read_mixed_msh(path)
            @test validate(reread).ok
            @test length(reread.periodic_links)==3
            reread_curve=only(filter(link->link.dim==1,reread.periodic_links))
            @test reread_curve.slave_entity==2
            @test reread_curve.master_entity==4
            @test length(reread_curve.slave_nodes)==5
            @test reread_curve.affine==mapping.affine
            @test reread.physical_names==projected.physical_names
            if version==2.2
                @test reread.entity_data===nothing
                @test reread.elementary_entities!==nothing
            else
                @test reread.entity_data!==nothing
                @test reread.elementary_entities===nothing
                @test mixed_crc(reread)==projected_crc
            end
            push!(format_crcs[version],mixed_crc(reread).sha)
        end
    end
    @test all(length(crcs)==1 for crcs in values(format_crcs))
    @test only(format_crcs[2.2])==
          "12a1eb50575a3af08273b1a0fdefca49d7e4b01b4898573e6346b6b61b4978c3"
    @test only(format_crcs[4.1])==projected_crc.sha

    owned_crc=mixed_crc(projected)
    mesh.coords[1,1]+=10
    model.physical_names[(2,13)]="changed"
    @test mixed_crc(projected)==owned_crc

    fresh_model,fresh_mesh=_periodic_projection_fixture()
    @test_throws ArgumentError model_to_mixed(fresh_model,fresh_mesh,99)
    with_segments=Mesh(
        fresh_mesh.coords;segs=reshape(Int32[1,2],2,1),tris=fresh_mesh.tris)
    @test_throws ArgumentError model_to_mixed(fresh_model,with_segments,1)
    tagged=Mesh(
        fresh_mesh.coords;tris=fresh_mesh.tris,
        tri_tag=fill(Int32(7),ntris(fresh_mesh)))
    @test_throws ArgumentError model_to_mixed(fresh_model,tagged,1)
    invalid=deepcopy(fresh_mesh);invalid.coords[1,1]=NaN
    @test_throws ArgumentError model_to_mixed(fresh_model,invalid,1)
    damaged=Mesh(fresh_mesh.coords;tris=fresh_mesh.tris[:,1:end-1])
    @test_throws ArgumentError model_to_mixed(fresh_model,damaged,1)
    volumetric=mesh_box(0,1,0,1,0,1;hmax=2)
    @test_throws ArgumentError model_to_mixed(fresh_model,volumetric,1)
    @test_throws ArgumentError model_to_mixed(fresh_model,fresh_mesh,true)
    @test_throws ArgumentError model_to_mixed(
        fresh_model,Mesh(zeros(3,0)),1)

    nonplanar=deepcopy(fresh_mesh);nonplanar.coords[3,1]=1e-6
    @test_throws ArgumentError model_to_mixed(fresh_model,nonplanar,1)
    nonmanifold=Mesh(
        [0.0 1.0 0.5 0.5 0.5;
         0.0 0.0 0.5 -0.5 0.25;
         0.0 0.0 0.0 0.0 0.0];
        tris=Int32[1 2 1;2 1 2;3 4 5])
    @test_throws ArgumentError model_to_mixed(fresh_model,nonmanifold,1)

    unsynchronized=Mesh(
        [0.0 1.0 1.0 0.0 0.0;
         0.0 0.0 1.0 1.0 0.5;
         0.0 0.0 0.0 0.0 0.0];
        tris=Int32[1 2 3;2 3 4;5 5 5])
    incompatibility=try
        model_to_mixed(fresh_model,unsynchronized,1)
        nothing
    catch err
        err
    end
    @test incompatibility isa ArgumentError
    @test occursin(
        "model_to_mixed: periodic Curve[2] mapping is incompatible",
        sprint(showerror,incompatibility))

    unsnapped=deepcopy(fresh_mesh)
    slave=first(model_periodic_nodes(
        fresh_model,fresh_mesh,1,2).slave_nodes)
    unsnapped.coords[1,slave]+=1e-13
    @test_throws ArgumentError model_to_mixed(fresh_model,unsnapped,1)

    embedded=deepcopy(fresh_model)
    add_point!(embedded,0.5,0.5,0;tag=5,mesh_size=0.5)
    embed!(embedded,0,[5],2,1)
    @test_throws ArgumentError model_to_mixed(embedded,fresh_mesh,1)

    double_periodic=deepcopy(fresh_model)
    translate_y=(1.0,0.0,0.0,0.0,
                 0.0,1.0,0.0,1.0,
                 0.0,0.0,1.0,0.0,
                 0.0,0.0,0.0,1.0)
    set_periodic!(double_periodic,1,[3],[1],translate_y)
    double_mesh=mesh_model_surface(double_periodic,1)
    double_projection=model_to_mixed(double_periodic,double_mesh,1)
    @test validate(double_projection).ok
    @test length(double_projection.periodic_links)==5
    @test sort([(Int(link.slave_entity),Int(link.master_entity))
                for link in double_projection.periodic_links if link.dim==0])==
          [(2,1),(3,2),(4,1)]
    @test sort([(Int(link.slave_entity),Int(link.master_entity))
                for link in double_projection.periodic_links if link.dim==1])==
          [(2,4),(3,1)]
    @test mixed_crc(double_projection).sha==
          "231b6d20877b7427735f96454a7f2ddcbe7d1d5368db8ab7cb9e05bbc5ffcf34"

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end

@testset "embedded native surface projection" begin
    model,mesh=_embedded_projection_fixture()
    projected=model_to_mixed(model,mesh,1)
    @test validate(projected).ok
    @test isempty(projected.periodic_links)
    @test [block.msh for block in projected.blocks]==[15,1,2]
    @test size(projected.blocks[1].nodes,2)==7
    @test Set(projected.entity_data.block_entities[1])==Set(Int32.(1:7))
    @test Set(projected.entity_data.block_entities[2])==Set(Int32.(1:5))
    @test projected.entity_data.entities[(2,1)].embedded_curves==Int32[5]
    @test projected.entity_data.entities[(2,1)].boundaries==Int32[1,2,3,4]
    @test projected.entity_data.entities[(1,5)].boundaries==Int32[5,-6]
    @test projected.entity_data.entities[(1,5)].physical_tags==Int32[32]
    @test projected.entity_data.entities[(0,7)].physical_tags==Int32[31]
    @test projected.physical_names==Dict(
        (0,31)=>"embedded points",(1,32)=>"embedded line",(2,33)=>"domain")

    point_nodes=Dict(Int(entity[2])=>node for (node,entity) in
        enumerate(projected.entity_data.node_entities) if entity[1]==0)
    @test Set(keys(point_nodes))==Set(1:7)
    @test projected.entity_data.node_entities[point_nodes[7]]==(0,Int32(7))
    embedded_line_cells=findall(==(Int32(5)),
        projected.entity_data.block_entities[2])
    @test length(embedded_line_cells)==2
    @test all(cell->point_nodes[7] in projected.blocks[2].nodes[:,cell],
              embedded_line_cells)
    @test all(cell->all(node->node in Set(values(point_nodes)) ||
            projected.entity_data.node_entities[node][1]==1,
            projected.blocks[2].nodes[:,cell]),embedded_line_cells)

    crc=mixed_crc(projected)
    @test crc.sha==
          "e762c7c566f1e5768ad1e2849302815dbfd9d19a14c1b3840abcefa4aedcaf43"
    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"embedded-$version-$binary.msh")
            @test write_mixed_msh(
                path,projected;version=version,binary=binary)==path
            reread=read_mixed_msh(path)
            @test validate(reread).ok
            if version==4.1
                @test reread.entity_data.entities[(2,1)].embedded_curves==Int32[5]
                @test mixed_crc(reread)==crc
            else
                @test reread.entity_data===nothing
                point_block=only(findall(block->block.msh==15,reread.blocks))
                line_block=only(findall(block->block.msh==1,reread.blocks))
                @test Set(reread.elementary_entities[point_block])==
                      Set(Int32.(1:7))
                @test Set(reread.elementary_entities[line_block])==
                      Set(Int32.(1:5))
            end
        end
    end

    periodic_embedded=deepcopy(model)
    add_point!(periodic_embedded,0.25,0.75,0;tag=8,mesh_size=0.5)
    add_point!(periodic_embedded,0.75,0.75,0;tag=9,mesh_size=0.5)
    add_line!(periodic_embedded,8,9;tag=6)
    embed!(periodic_embedded,1,[6],2,1)
    translate=(1.0,0.0,0.0,0.0,
               0.0,1.0,0.0,0.25,
               0.0,0.0,1.0,0.0,
               0.0,0.0,0.0,1.0)
    set_periodic!(periodic_embedded,1,[6],[5],translate)
    periodic_error=try
        model_to_mixed(periodic_embedded,mesh,1)
        nothing
    catch err
        err
    end
    @test periodic_error isa ArgumentError
    @test occursin("periodic embedded curves are not supported",
                   sprint(showerror,periodic_error))
end

@testset "holed native surface projection" begin
    model,mesh=_holed_projection_fixture()
    projected=model_to_mixed(model,mesh,1)
    @test validate(projected).ok
    @test isempty(projected.periodic_links)
    @test projected.entity_data!==nothing
    @test projected.entity_data.entities[(2,1)].boundaries==
          Int32[1,2,3,4,-8,-7,-6,-5]
    @test Set(projected.entity_data.block_entities[2])==Set(Int32.(1:8))
    @test projected.entity_data.entities[(1,5)].physical_tags==Int32[22]
    @test projected.entity_data.entities[(1,4)].physical_tags==Int32[]
    @test projected.physical_names==Dict((1,22)=>"hole",(2,23)=>"domain")
    @test mixed_crc(projected).sha==
          "654d305c58cfe9db2f87ac1424a31912863a21f75196beabc3f05c37d8f6e73f"
end
