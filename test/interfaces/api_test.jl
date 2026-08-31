using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, node, ntets, tet_volume, validate

const _API=Tessella.API

@testset "owned and validated API session" begin
    _API.finalize()
    @test_throws ArgumentError _API.option("Mesh.MeshSizeFactor")
    @test_throws ArgumentError _API.option("Mesh.MeshSizeFactor",2.0)
    @test_throws ArgumentError _API.model.add_box(0,0,0,1,1,1;tag=1)
    @test_throws ArgumentError _API.mesh.get()
    @test_throws ArgumentError _API.mesh.set_size([(0,1)],0.5)
    @test_throws ArgumentError _API.open_geo!("missing.geo")

    try
        @test _API.initialize()===nothing
        @test _API.option("Mesh.MeshSizeMin")==0.0
        @test _API.option("Mesh.MeshSizeMax")==1.0e22
        @test _API.option("Mesh.MeshSizeFactor")==1.0
        @test _API.option("Mesh.MeshSizeFactor",2)==2.0
        @test_throws ArgumentError _API.option("No.Such.Option")
        @test_throws ArgumentError _API.option("No.Such.Option",1.0)
        @test_throws ArgumentError _API.option("Mesh.MeshSizeFactor",true)
        @test_throws ArgumentError _API.option("Mesh.MeshSizeFactor",Inf)
        @test_throws ArgumentError _API.option("Mesh.MeshSizeFactor",0.0)
        @test_throws ArgumentError _API.option("Mesh.MeshSizeMin",-1.0)
        @test_throws ArgumentError _API.option("Mesh.MeshSizeMax",0.0)

        @test _API.option("Mesh.MeshSizeMax",1.0)==1.0
        @test_throws ArgumentError _API.option("Mesh.MeshSizeMin",2.0)
        @test _API.option("Mesh.MeshSizeMin")==0.0
        @test _API.option("Mesh.MeshSizeMax")==1.0
        @test _API.option("Mesh.MeshSizeMin",0.5)==0.5
        @test_throws ArgumentError _API.option("Mesh.MeshSizeMax",0.25)
        @test _API.option("Mesh.MeshSizeMin")==0.5
        @test _API.option("Mesh.MeshSizeMax")==1.0

        _API.initialize()
        @test _API.option("Mesh.MeshSizeMin")==0.0
        @test _API.option("Mesh.MeshSizeMax")==1.0e22
        @test _API.option("Mesh.MeshSizeFactor")==1.0

        @test _API.model.add_box(0,0,0,1,1,1;tag=1)==1
        generated=_API.mesh.generate(3)
        @test validate(generated).ok
        expected_crc=mesh_crc(generated)
        @test expected_crc.sha==
              "e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83"
        node_tags,node_coordinates,node_parameters=_API.mesh.get_nodes()
        @test node_tags==UInt64.(1:9)
        @test reshape(node_coordinates,3,:)==generated.coords
        @test isempty(node_parameters)
        element_types,element_tags,element_nodes=_API.mesh.get_elements()
        @test element_types==Int32[4]
        @test element_tags==[UInt64.(1:12)]
        @test element_nodes==[UInt64.(vec(generated.tets))]
        @test _API.mesh.get_element_types()==Int32[4]
        @test _API.mesh.get_elements_by_type(4)==
              (UInt64.(1:12),UInt64.(vec(generated.tets)))
        type_node_tags,type_node_coordinates,type_node_parameters=
            _API.mesh.get_nodes_by_element_type(4)
        @test type_node_tags==UInt64.(vec(generated.tets))
        @test reshape(type_node_coordinates,3,:)==
              generated.coords[:,Int.(type_node_tags)]
        @test isempty(type_node_parameters)
        @test length(_API.mesh.get_barycenters(4,-1,false,false))==36
        @test length(_API.mesh.get_element_edge_nodes(4))==144
        @test length(_API.mesh.get_element_face_nodes(4,3))==144
        @test _API.mesh.get_max_node_tag()==UInt64(9)
        @test _API.mesh.get_max_element_tag()==UInt64(12)

        cached=_API.mesh.get()
        @test cached!==generated && cached.coords!==generated.coords
        generated.coords[1,1]+=17.0
        @test mesh_crc(_API.mesh.get())==expected_crc
        cached.coords[2,1]+=19.0
        @test mesh_crc(_API.mesh.get())==expected_crc

        @test_throws ArgumentError _API.model.add_box(0,0,0,1,1,1;tag=1)
        @test mesh_crc(_API.mesh.get())==expected_crc
        @test_throws ArgumentError _API.mesh.generate(false)
        @test_throws ArgumentError _API.mesh.generate(4)
        @test_throws ArgumentError _API.mesh.generate(big(typemax(Int))+1)
        @test mesh_crc(_API.mesh.get())==expected_crc

        @test _API.model.add_box(2,0,0,1,1,1;tag=2)==2
        @test_throws ArgumentError _API.mesh.get()
        @test_throws ArgumentError _API.mesh.generate(3)
    finally
        _API.finalize()
    end

    try
        _API.initialize()
        _API.model.add_box(0,0,0,2,1,1;tag=1)
        _API.model.add_box(0,0,0,1,1,1;tag=2)
        @test _API.model.boolean_difference(1,2;tag=3)==3
        cut=_API.mesh.generate(3)
        @test validate(cut).ok && ntets(cut)>0
        cut_volume=sum(tet_volume(node(cut,cut.tets[1,t]),node(cut,cut.tets[2,t]),
                                  node(cut,cut.tets[3,t]),node(cut,cut.tets[4,t]))
                       for t in 1:ntets(cut))
        @test cut_volume≈1.0 atol=1e-12
    finally
        _API.finalize()
    end

    mktempdir() do directory
        valid=joinpath(directory,"box.geo")
        invalid=joinpath(directory,"invalid.geo")
        periodic=joinpath(directory,"periodic.geo")
        write(valid,"Box(1) = {0, 0, 0, 1, 1, 1};\n")
        write(invalid,"Extrude {0, 0, 1} { Volume{1}; }\n")
        write(periodic,"""
            Point(1) = {0, 0, 0, 0.5};
            Point(2) = {1, 0, 0, 0.5};
            Point(3) = {1, 1, 0, 0.5};
            Point(4) = {0, 1, 0, 0.5};
            Line(1) = {1, 2};
            Line(2) = {2, 3};
            Line(3) = {3, 4};
            Line(4) = {4, 1};
            Curve Loop(1) = {1, 2, 3, 4};
            Plane Surface(1) = {1};
            periodicSlave = Sqrt(4);
            periodicMaster = 8 / 2;
            periodicShift = Cos(0);
            periodicZero = Atan2(0, 1);
            Periodic Curve {periodicSlave} = {periodicMaster}
              Translate {periodicShift, Sin(0), periodicZero};
            """)

        @test_throws ArgumentError _API.open_geo!(valid)
        try
            _API.initialize()
            @test _API.model.add_box(0,0,0,1,1,1;tag=1)==1
            @test_throws ArgumentError _API.open_geo!(invalid)
            @test_throws ArgumentError _API.model.add_box(0,0,0,1,1,1;tag=1)
            @test _API.model.add_box(2,0,0,1,1,1;tag=2)==2

            result=_API.open_geo!(valid)
            @test result.mesh===nothing
            @test add_box!(result.model,2,0,0,1,1,1;tag=2)==2
            @test _API.model.add_box(2,0,0,1,1,1;tag=2)==2
            @test_throws ArgumentError _API.mesh.get()

            periodic_result=_API.open_geo!(periodic;mesh_dim=2)
            @test periodic_result.mesh!==nothing
            periodic_crc=mesh_crc(periodic_result.mesh)
            @test periodic_crc.sha==
                  "3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5"
            periodic_mapping=_API.mesh.get_periodic_nodes(1,2)
            @test periodic_mapping.master_entity==4
            @test length(periodic_mapping.slave_nodes)==5
            periodic_result.mesh.coords[1,1]+=10
            @test mesh_crc(_API.mesh.get())==periodic_crc
        finally
            _API.finalize()
        end
    end

    @test _API.finalize()===nothing
    @test _API.finalize()===nothing
    @test isempty(Docs.undocumented_names(Tessella.API;private=false))
end

@testset "owned Physical-group lifecycle through API" begin
    _API.finalize()
    try
        _API.initialize()
        @test _API.model.add_point(0,0,0;tag=1)==1
        @test _API.model.add_point(1,0,0;tag=2)==2
        @test _API.model.add_point(0,1,0;tag=3)==3
        @test _API.model.add_line(1,2;tag=1)==1
        @test _API.model.add_line(2,3;tag=2)==2
        @test _API.model.add_line(3,1;tag=3)==3
        @test _API.model.add_curve_loop([1,2,3];tag=10)==10
        @test _API.model.add_plane_surface([10];tag=1)==1
        @test _API.model.add_physical_group(0,[1];name="first")==1
        @test _API.model.add_physical_group(1,[1];name="second")==2
        @test _API.model.add_physical_group(0,[2];name="first")==3
        @test _API.model.add_physical_group(2,[1];name="first")==4
        @test _API.CURRENT[].physical_names==Dict(
            (0,1)=>"first",(1,2)=>"second",(2,4)=>"first")
        @test _API.CURRENT[].physical_tag_max==4

        expected_groups=[(0,1),(0,3),(1,2),(2,4)]
        @test _API.model.get_physical_groups()==expected_groups
        @test _API.model.get_physical_groups(0)==[(0,1),(0,3)]
        @test _API.model.get_entities_for_physical_group(0,1)==[1]
        @test _API.model.get_physical_groups_for_entity(0,1)==[1]
        @test _API.model.get_physical_name(0,1)=="first"
        @test _API.model.get_physical_name(0,3)==""
        @test _API.model.get_physical_name(0,99)==""
        @test _API.model.get_entities_for_physical_name("first")==
              [(0,1),(2,1)]
        group_pairs,entity_pairs=_API.model.get_physical_groups_entities()
        @test group_pairs==expected_groups
        @test entity_pairs==[[(0,1)],[(0,2)],[(1,1)],[(2,1)]]

        detached_groups=_API.model.get_physical_groups()
        detached_members=_API.model.get_entities_for_physical_group(0,1)
        detached_pairs=_API.model.get_physical_groups_entities()[2]
        push!(detached_groups,(3,99))
        push!(detached_members,99)
        push!(detached_pairs[1],(0,99))
        @test _API.model.get_physical_groups()==expected_groups
        @test _API.model.get_entities_for_physical_group(0,1)==[1]
        @test _API.model.get_physical_groups_entities()[2][1]==[(0,1)]

        generated=_API.mesh.generate(2)
        @test validate(generated).ok
        cached=_API.LAST_MESH[]
        @test _API.model.get_physical_groups()==expected_groups
        @test _API.LAST_MESH[]===cached
        @test _API.model.set_physical_name(0,3,"first")===nothing
        @test _API.model.set_physical_name(0,99,"ghost")===nothing
        @test _API.model.remove_physical_name("missing")===nothing
        @test _API.model.remove_physical_groups([(0,99)])===nothing
        @test _API.LAST_MESH[]===cached

        stable_groups=_API.model.get_physical_groups()
        @test_throws ArgumentError _API.model.remove_physical_groups(
            [(0,1),(4,1)])
        @test _API.model.get_physical_groups()==stable_groups
        @test _API.LAST_MESH[]===cached
        @test_throws ArgumentError _API.model.remove_physical_groups([(0,1),1])
        @test_throws ArgumentError _API.model.get_physical_groups(4)
        @test_throws ArgumentError _API.model.get_entities_for_physical_group(0,99)
        @test_throws ArgumentError _API.model.get_physical_groups_for_entity(0,99)
        @test_throws ArgumentError _API.model.get_entities_for_physical_name("missing")
        @test_throws ArgumentError _API.model.set_physical_name(0,1,1)
        @test_throws ArgumentError _API.model.remove_physical_name(1)

        @test _API.model.set_physical_name(0,3,"third")===nothing
        @test _API.model.get_physical_name(0,3)=="third"
        @test _API.LAST_MESH[]===nothing
        _API.mesh.generate(2)
        @test _API.model.remove_physical_name("first")===nothing
        @test _API.model.get_physical_name(0,1)==""
        @test _API.model.get_physical_name(2,4)==""
        @test _API.model.get_physical_name(0,3)=="third"
        @test _API.LAST_MESH[]===nothing
        _API.mesh.generate(2)
        @test _API.model.remove_physical_groups([(1,2)])===nothing
        @test _API.model.get_physical_groups()==[(0,1),(0,3),(2,4)]
        @test _API.CURRENT[].physical_tag_max==4
        @test _API.LAST_MESH[]===nothing
        @test _API.model.add_physical_group(1,[1];name="replacement")==5
        _API.mesh.generate(2)
        @test _API.model.remove_physical_groups()===nothing
        @test isempty(_API.model.get_physical_groups())
        @test _API.LAST_MESH[]===nothing
        @test _API.model.add_physical_group(2,[1];name="fresh")==6
    finally
        _API.finalize()
    end
end

@testset "owned point mesh-size API" begin
    _API.finalize()
    try
        _API.initialize()
        @test _API.model.add_box(0,0,0,1,1,1;tag=1)==1
        @test _API.model.add_point(0.25,0.25,0.25;tag=101)==101
        @test _API.model.add_point(0.75,0.75,0.75;tag=102)==102
        initial=_API.mesh.generate(3)
        @test validate(initial).ok
        @test mesh_crc(initial).sha==
              "e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83"

        @test _API.mesh.set_size((0=>101,0=>102),0.25)===nothing
        @test _API.CURRENT[].point_size[101]==0.25
        @test _API.CURRENT[].point_size[102]==0.25
        @test_throws ArgumentError _API.mesh.get()

        refreshed=_API.mesh.generate(3)
        @test validate(refreshed).ok
        refreshed_crc=mesh_crc(refreshed)
        @test refreshed_crc==mesh_crc(initial)
        stable_sizes=copy(_API.CURRENT[].point_size)
        for (dim_tags,size) in (
                ([(1,101)],0.5),([(false,101)],0.5),
                ([(0,101),(0,999)],0.5),([(0,101.0)],0.5),
                ([(0,101)],0.0),((),0.5),([101],0.5),
                ("(0, 101)",0.5))
            @test_throws ArgumentError _API.mesh.set_size(dim_tags,size)
            @test _API.CURRENT[].point_size==stable_sizes
            @test mesh_crc(_API.mesh.get())==refreshed_crc
        end
    finally
        _API.finalize()
    end
end

@testset "explicit volume-shell API" begin
    _API.finalize()
    try
        _API.initialize()
        points=((0.0,0.0,0.0),(1.0,0.0,0.0),
                (0.0,1.0,0.0),(0.0,0.0,1.0))
        for (tag,point) in pairs(points)
            @test _API.model.add_point(point...;tag=tag)==tag
        end
        edges=((1,2),(2,3),(3,1),(1,4),(2,4),(3,4))
        for (tag,(first_point,last_point)) in pairs(edges)
            @test _API.model.add_line(first_point,last_point;tag=tag)==tag
        end
        loops=((1,2,3),(1,5,-4),(2,6,-5),(3,4,-6))
        for (tag,curves) in pairs(loops)
            @test _API.model.add_curve_loop(curves;tag=tag)==tag
            @test _API.model.add_plane_surface([tag];tag=tag)==tag
        end
        @test _API.model.add_surface_loop([1,2,3,4];tag=1)==1
        @test _API.model.add_volume([1];tag=1)==1
        generated=_API.mesh.generate(3)
        @test validate(generated).ok
        @test ntets(generated)>0
        volume=sum(tet_volume(
            node(generated,generated.tets[1,cell]),
            node(generated,generated.tets[2,cell]),
            node(generated,generated.tets[3,cell]),
            node(generated,generated.tets[4,cell])) for cell in 1:ntets(generated))
        @test volume≈1/6 atol=1e-12
        @test mesh_crc(generated).sha==
              "71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b"
        expected=mesh_crc(generated)
        @test_throws ArgumentError _API.model.add_surface_loop([1];tag=2)
        @test mesh_crc(_API.mesh.get())==expected
    finally
        _API.finalize()
    end
end

@testset "periodic volume-boundary API" begin
    fixture=normpath(joinpath(
        @__DIR__,"..","fixtures","periodic_surface_volume.geo"))
    source=replace(
        read(fixture,String),
        "Periodic Surface {4} = {6} Translate {1, 0, 0};\n"=>"",
        "Periodic Surface {5} = {3} Translate {0, 1, 0};\n"=>"")
    mktemp() do path,io
        write(io,source)
        close(io)
        _API.finalize()
        try
            _API.initialize()
            built=_API.open_geo!(path)
            @test built.mesh===nothing
            translate_x=(1.0,0.0,0.0,1.0,
                         0.0,1.0,0.0,0.0,
                         0.0,0.0,1.0,0.0,
                         0.0,0.0,0.0,1.0)
            translate_y=(1.0,0.0,0.0,0.0,
                         0.0,1.0,0.0,1.0,
                         0.0,0.0,1.0,0.0,
                         0.0,0.0,0.0,1.0)
            @test _API.mesh.set_periodic(
                2,[4],[6],translate_x)===nothing
            @test _API.mesh.set_periodic(
                2,[5],[3],translate_y)===nothing
            @test_throws ArgumentError _API.mesh.get()

            generated=_API.mesh.generate(3)
            @test validate(generated).ok
            expected=mesh_crc(generated)
            @test expected.sha==
                  "2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808"
            for (slave,master,pairs) in ((4,6,5),(5,3,4))
                mapping=_API.mesh.get_periodic_nodes(2,slave)
                @test mapping.master_entity==master
                @test length(mapping.slave_nodes)==
                      length(mapping.master_nodes)==pairs
            end
            first_mapping=_API.mesh.get_periodic_nodes(2,4)
            first_mapping.master_nodes[1]=0
            @test _API.mesh.get_periodic_nodes(2,4).master_nodes!=
                  first_mapping.master_nodes
            generated.coords[1,1]+=10
            @test mesh_crc(_API.mesh.get())==expected

            refined=_API.mesh.refine()
            @test mesh_crc(refined).sha==
                  "5db960c919e54384ecf25fe27144a8636d9b8557d2afe052c77c7c545ffe8722"
            for (slave,master,pairs) in ((4,6,13),(5,3,9))
                mapping=_API.mesh.get_periodic_nodes(2,slave)
                @test mapping.master_entity==master
                @test length(mapping.slave_nodes)==
                      length(mapping.master_nodes)==pairs
            end
        finally
            _API.finalize()
        end
    end
end

@testset "periodic API ownership and cache invalidation" begin
    _API.finalize()
    @test_throws ArgumentError _API.mesh.set_periodic(
        1,[2],[4],(1.0,0.0,0.0,0.0,
                    0.0,1.0,0.0,0.0,
                    0.0,0.0,1.0,0.0,
                    0.0,0.0,0.0,1.0))
    @test_throws ArgumentError _API.mesh.get_periodic_nodes(1,2)

    try
        _API.initialize()
        for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),
                                      (1.0,1.0),(0.0,1.0)))
            @test _API.model.add_point(
                x,y,0;tag=tag,meshSize=0.5)==tag
        end
        for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            @test _API.model.add_line(first,last;tag=tag)==tag
        end
        @test _API.model.add_curve_loop([1,2,3,4];tag=1)==1
        @test _API.model.add_plane_surface([1];tag=1)==1

        translation=Float64[
            1,0,0,1,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1,
        ]
        @test _API.mesh.set_periodic(
            1,[2],[4],translation)===nothing
        translation[4]=99
        @test_throws ArgumentError _API.mesh.get_periodic_nodes(1,2)

        generated=_API.mesh.generate(2)
        expected=mesh_crc(generated)
        @test expected.sha==
              "3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5"
        mapping=_API.mesh.get_periodic_nodes(1,2)
        @test mapping.master_entity==4
        @test mapping.affine[4]==1
        @test length(mapping.slave_nodes)==length(mapping.master_nodes)==5
        cached=_API.mesh.get()
        for (slave,master) in zip(mapping.slave_nodes,mapping.master_nodes)
            @test Tuple(cached.coords[:,slave])==
                  (cached.coords[1,master]+1,cached.coords[2,master],
                   cached.coords[3,master])
        end
        mapping.master_nodes[1]=1
        @test first(_API.mesh.get_periodic_nodes(1,2).master_nodes)!=1
        generated.coords[1,1]+=10
        @test mesh_crc(_API.mesh.get())==expected

        identity=(1.0,0.0,0.0,0.0,
                  0.0,1.0,0.0,0.0,
                  0.0,0.0,1.0,0.0,
                  0.0,0.0,0.0,1.0)
        @test_throws ArgumentError _API.mesh.set_periodic(
            2,[3],[1],identity)
        @test mesh_crc(_API.mesh.get())==expected
        @test_throws ArgumentError _API.mesh.get_periodic_nodes(1,4)

        translate_y=(1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,1.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0)
        @test _API.mesh.set_periodic(
            1,[3],[1],translate_y)===nothing
        @test_throws ArgumentError _API.mesh.get()
        @test_throws ArgumentError _API.mesh.get_periodic_nodes(1,2)

        double_periodic=_API.mesh.generate(2)
        @test validate(double_periodic).ok
        @test mesh_crc(double_periodic).sha==
              "95ef6d0db94505d4f35ff870af09e952d74a32508a338b3994af347b406e9d05"
        @test length(_API.mesh.get_periodic_nodes(1,2).slave_nodes)==5
        @test length(_API.mesh.get_periodic_nodes(1,3).slave_nodes)==5
    finally
        _API.finalize()
    end
end

@testset "embedded periodic-curve API" begin
    _API.finalize()
    try
        _API.initialize()
        coordinates=((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0),
                     (0.25,0.25),(0.75,0.25),
                     (0.25,0.75),(0.75,0.75),(0.5,0.25))
        for (tag,(x,y)) in pairs(coordinates)
            @test _API.model.add_point(
                x,y,0;tag=tag,meshSize=0.5)==tag
        end
        endpoints=((1,2),(2,3),(3,4),(4,1),(5,6),(7,8))
        for (tag,(first_point,last_point)) in pairs(endpoints)
            @test _API.model.add_line(
                first_point,last_point;tag=tag)==tag
        end
        @test _API.model.add_curve_loop([1,2,3,4];tag=1)==1
        @test _API.model.add_plane_surface([1];tag=1)==1
        @test _API.model.embed(0,[9],2,1)==1
        @test _API.model.embed(1,[5,6],2,1)==1
        translation=(1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,0.5,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0)
        @test _API.mesh.set_periodic(
            1,[6],[5],translation)===nothing
        generated=_API.mesh.generate(2)
        @test validate(generated).ok
        @test mesh_crc(generated).sha==
              "9794a65ea5402683d0d50612522c2f71f7c98ec2a9f6b9e6b49a61e62cd85cf2"
        mapping=_API.mesh.get_periodic_nodes(1,6)
        @test mapping.master_entity==5
        @test length(mapping.slave_nodes)==3
        cached=_API.mesh.get()
        for (slave,master) in zip(mapping.slave_nodes,mapping.master_nodes)
            @test Tuple(cached.coords[:,slave])==
                  (cached.coords[1,master],cached.coords[2,master]+0.5,
                   cached.coords[3,master])
        end
    finally
        _API.finalize()
    end
end

@testset "periodic-curve dependency graph API" begin
    _API.finalize()
    try
        _API.initialize()
        coordinates=(
            (1,0.0,0.0),(2,1.0,0.0),(3,1.0,1.0),(4,0.0,1.0),
            (101,0.2,0.2),(102,0.8,0.2),
            (103,0.2,0.5),(104,0.8,0.5),
            (105,0.2,0.8),(106,0.8,0.8),(107,0.425,0.8))
        for (tag,x,y) in coordinates
            @test _API.model.add_point(
                x,y,0;tag=tag,meshSize=0.4)==tag
        end
        endpoints=((1,1,2),(2,2,3),(3,3,4),(4,4,1),
                   (30,101,102),(20,103,104),(10,106,105))
        for (tag,first_point,last_point) in endpoints
            @test _API.model.add_line(
                first_point,last_point;tag=tag)==tag
        end
        @test _API.model.add_curve_loop([1,2,3,4];tag=1)==1
        @test _API.model.add_plane_surface([1];tag=1)==1
        @test _API.model.embed(0,[107],2,1)==1
        @test _API.model.embed(1,[30,20,10],2,1)==1
        translation=(1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,0.3,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0)
        @test _API.mesh.set_periodic(
            1,[20,10],[30,20],translation)===nothing
        generated=_API.mesh.generate(2)
        @test validate(generated).ok
        @test mesh_crc(generated).sha==
              "dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338"
        cached=_API.mesh.get()
        for (slave_entity,master_entity) in ((10,20),(20,30))
            mapping=_API.mesh.get_periodic_nodes(1,slave_entity)
            @test mapping.master_entity==master_entity
            @test length(mapping.slave_nodes)==9
            for (slave,master) in zip(mapping.slave_nodes,
                                      mapping.master_nodes)
                @test Tuple(cached.coords[:,slave])==
                      (cached.coords[1,master],cached.coords[2,master]+0.3,
                       cached.coords[3,master])
            end
        end
        @test_throws ArgumentError _API.mesh.set_periodic(1,[30],[10],(
            1.0,0.0,0.0,0.0,
            0.0,1.0,0.0,-0.6,
            0.0,0.0,1.0,0.0,
            0.0,0.0,0.0,1.0))
        @test mesh_crc(_API.mesh.get())==mesh_crc(cached)
    finally
        _API.finalize()
    end
end
