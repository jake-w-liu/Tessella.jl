using Test
using LinearAlgebra
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume, validate

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
              "eae8751b0dad3b89f2d7a4416ea079a352a3bd6b8eff31ef7eb7d8bb78d8509a"
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
        # Adding a volume invalidates the cache; generating again meshes both
        # volumes into the merged cache with per-entity classification.
        _API.mesh.generate(3)
        @test ntets(_API.mesh.get())==24
        @test _API.mesh.get_elements(3,1)[1]==Int32[4]
        @test _API.mesh.get_elements(3,2)[1]==Int32[4]
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
        write(valid,"SetFactory(\"OpenCASCADE\");\nBox(1) = {0, 0, 0, 1, 1, 1};\n")
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
        @test _API.model.remove_physical_name("missing")===nothing
        @test _API.LAST_MESH[]===cached
        # `setPhysicalName` binds names to groupless tags unconditionally, and
        # `removePhysicalGroup` erases name bindings unconditionally — both are
        # real mutations that invalidate the cached mesh.
        @test _API.model.set_physical_name(0,99,"ghost")===nothing
        @test _API.model.get_physical_name(0,99)=="ghost"
        @test _API.LAST_MESH[]===nothing
        @test _API.model.remove_physical_groups([(0,99)])===nothing
        @test _API.model.get_physical_name(0,99)==""
        @test _API.LAST_MESH[]===nothing

        stable_groups=_API.model.get_physical_groups()
        @test_throws ArgumentError _API.model.remove_physical_groups(
            [(0,1),(4,1)])
        @test _API.model.get_physical_groups()==stable_groups
        @test _API.LAST_MESH[]===nothing
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
              "eae8751b0dad3b89f2d7a4416ea079a352a3bd6b8eff31ef7eb7d8bb78d8509a"

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
              "03cfdc7130ae46c251a59237671e3bb37dcba83b690accde3770ac4a78d4cbb4"
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
                  "a58374071a4c485a339e1c5b48b8b0f3e69bf362ff0e41f57ca1a665139e81df"
            for (slave,master,pairs) in ((4,6,25),(5,3,25))
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
                  "97cc7537053d447ca2c9bff1be0be82812c6ce1b7b384d54af3f50c5ee038408"
            for (slave,master,pairs) in ((4,6,81),(5,3,81))
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
        # A cycle-closing declaration stores like upstream (`setMeshMaster`
        # has no cycle check); the mesh cache invalidates on the constraint
        # change and the next `generate` converges the cyclic parameter
        # fixpoint, including the 30 -> 10 closing link at -0.6.
        @test _API.mesh.set_periodic(1,[30],[10],(
            1.0,0.0,0.0,0.0,
            0.0,1.0,0.0,-0.6,
            0.0,0.0,1.0,0.0,
            0.0,0.0,0.0,1.0))===nothing
        regen=_API.mesh.generate(2)
        @test validate(regen).ok
        recached=_API.mesh.get()
        closing=_API.mesh.get_periodic_nodes(1,30)
        @test closing.master_entity==10
        @test length(closing.slave_nodes)==9
        for (slave,master) in zip(closing.slave_nodes,
                                  closing.master_nodes)
            @test Tuple(recached.coords[:,slave])==
                  (recached.coords[1,master],recached.coords[2,master]-0.6,
                   recached.coords[3,master])
        end
    finally
        _API.finalize()
    end
end

@testset "model name/file and periodic/renumbering parity" begin
    _API.finalize()
    @test_throws ArgumentError _API.model.get_current()
    @test_throws ArgumentError _API.model.set_current("unnamed")
    @test_throws ArgumentError _API.model.get_file_name()
    @test_throws ArgumentError _API.model.set_file_name("x.msh")
    try
        _API.initialize()
        @test _API.model.get_current()==""
        @test _API.model.get_file_name()==""
        @test _API.model.set_current("")===nothing
        @test_throws ArgumentError _API.model.set_current("other")
        @test _API.model.set_file_name("fixture.msh")===nothing
        @test _API.model.get_file_name()=="fixture.msh"
        for call in (()->_API.model.set_current(3),
                     ()->_API.model.set_current(true),
                     ()->_API.model.set_file_name(3),
                     ()->_API.model.set_file_name(nothing))
            @test_throws ArgumentError call()
        end
        # initialize() resets session model state.
        _API.initialize()
        @test _API.model.get_current()==""
        @test _API.model.get_file_name()==""

        _API.model.add_point(0,0,0;tag=1);_API.model.add_point(1,0,0;tag=2)
        _API.model.add_point(0,1,0;tag=3);_API.model.add_point(1,1,0;tag=4)
        _API.model.add_line(1,2;tag=1);_API.model.add_line(3,4;tag=2)
        _API.model.add_line(2,4;tag=3);_API.model.add_line(1,3;tag=4)
        _API.model.add_curve_loop([1,3,-2,-4];tag=1)
        _API.model.add_plane_surface([1];tag=1)
        _API.mesh.set_periodic(1,[2],[1],(
            1.0,0.0,0.0,0.0,
            0.0,1.0,0.0,1.0,
            0.0,0.0,1.0,0.0,
            0.0,0.0,0.0,1.0))
        # Entities with no periodic master map to themselves, like Gmsh 4.15.2.
        @test _API.mesh.get_periodic(1,[1,2,3,4])==Int32[1,1,3,4]
        @test _API.mesh.get_periodic(2,[1])==Int32[1]
        @test _API.mesh.get_periodic(0,[1])==Int32[1]
        @test _API.mesh.get_periodic(1,Int[])==Int32[]
        for call in (()->_API.mesh.get_periodic(1,[99]),
                     ()->_API.mesh.get_periodic(1,[2.5]),
                     ()->_API.mesh.get_periodic(1,[true]),
                     ()->_API.mesh.get_periodic(-1,[1]),
                     ()->_API.mesh.get_periodic(4,[1]),
                     ()->_API.mesh.get_periodic(1,3))
            @test_throws ArgumentError call()
        end
        # removeConstraints retains periodic relations, embeddings, and Point
        # sizes in Gmsh 4.15.2; the native model stores none of the cleared
        # attribute kinds, so this is a validated no-op.
        @test _API.mesh.remove_constraints()===nothing
        @test _API.mesh.remove_constraints([(0,1),(1,2),(2,1)])===nothing
        @test _API.mesh.get_periodic(1,[2])==Int32[1]
        for call in (()->_API.mesh.remove_constraints([(1,99)]),
                     ()->_API.mesh.remove_constraints([(4,1)]),
                     ()->_API.mesh.remove_constraints([(1,)]),
                     ()->_API.mesh.remove_constraints(3))
            @test_throws ArgumentError call()
        end

        generated=_API.mesh.generate(2)
        @test validate(generated).ok
        old_tags,new_tags=_API.mesh.compute_renumbering()
        @test old_tags==UInt64.(1:nnodes(generated))
        @test sort(Int.(new_tags))==collect(1:nnodes(generated))
        triangle_tags,_=_API.mesh.get_elements_by_type(2)
        restricted_old,restricted_new=_API.mesh.compute_renumbering(
            "RCMK",triangle_tags[1:2])
        @test length(restricted_old)==length(restricted_new)
        @test issorted(restricted_old)
        @test sort(Int.(restricted_new))==collect(1:length(restricted_new))
        for call in (()->_API.mesh.compute_renumbering("Hilbert"),
                     ()->_API.mesh.compute_renumbering(3),
                     ()->_API.mesh.compute_renumbering("RCMK",[999]),
                     ()->_API.mesh.compute_renumbering("RCMK",[-1]),
                     ()->_API.mesh.compute_renumbering("RCMK",3))
            @test_throws ArgumentError call()
        end
        # optimize() runs the boundary-preserving tet optimizer on the cache;
        # a 2-D cache is unchanged, entity scoping freezes unselected nodes,
        # and invalid scopes/methods fail explicitly.
        @test _API.mesh.optimize()===nothing
        @test validate(_API.mesh.get()).ok
        @test _API.mesh.optimize("",false,1,[(2,1)])===nothing
        for call in (()->_API.mesh.optimize("Netgen"),
                     ()->_API.mesh.optimize(3),
                     ()->_API.mesh.optimize("",3,1),
                     ()->_API.mesh.optimize("",false,-1),
                     ()->_API.mesh.optimize("",false,1,[(3,99)]))
            @test_throws ArgumentError call()
        end
        # Per-element visibility is raw display state like Gmsh 4.15.2: default
        # 1, stored values pass through, unknown tags silently report 0, and
        # the state resets when the cache is replaced.
        triangle_tags,_=_API.mesh.get_elements_by_type(2)
        @test _API.mesh.get_visibility(triangle_tags)==
              fill(Int32(1),length(triangle_tags))
        @test _API.mesh.set_visibility(triangle_tags[1:1],0)===nothing
        @test _API.mesh.get_visibility(triangle_tags[1:2])==Int32[0,1]
        @test _API.mesh.set_visibility([99999],5)===nothing
        @test _API.mesh.get_visibility([99999,triangle_tags[1]])==Int32[0,0]
        for call in (()->_API.mesh.set_visibility([1],"x"),
                     ()->_API.mesh.set_visibility([1.5],1),
                     ()->_API.mesh.set_visibility(1,0),
                     ()->_API.mesh.get_visibility(3))
            @test_throws ArgumentError call()
        end
        # Per-window visibility is validated display state: window indices are
        # non-negative integers and unknown element tags are dropped silently.
        @test _API.model.set_visibility_per_window(0)===nothing
        @test _API.model.set_visibility_per_window(1,2)===nothing
        @test _API.mesh.set_visibility_per_window(1,0)===nothing
        @test _API.mesh.set_visibility_per_window(99999,0,3)===nothing
        for call in (()->_API.model.set_visibility_per_window(1,-1),
                     ()->_API.model.set_visibility_per_window("x"),
                     ()->_API.mesh.set_visibility_per_window(1.5,0),
                     ()->_API.mesh.set_visibility_per_window(1,0,-2),
                     ()->_API.mesh.set_visibility_per_window(1,"x"))
            @test_throws ArgumentError call()
        end
        _API.mesh.clear()
        @test _API.mesh.get_visibility([1])==Int32[0]
        _API.mesh.generate(2)
        metadata=Docs.meta(Tessella.API.mesh)
        for name in (:get_periodic,:remove_constraints,:compute_renumbering,
                     :optimize,:set_visibility,:get_visibility,
                     :set_visibility_per_window)
            @test haskey(metadata,Docs.Binding(Tessella.API.mesh,name))
        end
        # Multi-model lifecycle matches Gmsh 4.15.2: `initialize` creates one
        # unnamed model (""); `add` always appends a fresh model and selects
        # it; `set_current` switches between slots, preserving each model's
        # geometry, mesh, and visibility state; `remove` deletes the current
        # slot and selects the last remaining one.
        @test _API.model.list()==[""]
        @test _API.model.add("second")===nothing
        @test _API.model.list()==["","second"]
        @test _API.model.get_current()=="second"
        # The new model starts empty; the first model keeps its state.
        @test _API.model.get_entities()==Tuple{Int,Int}[]
        @test_throws ArgumentError _API.mesh.get()
        @test _API.model.add_point(9,9,9;tag=8)==8
        _API.model.set_current("")
        @test _API.model.get_entities(0)==[(0,1),(0,2),(0,3),(0,4)]
        @test size(_API.mesh.get().tris,2)>0
        _API.model.set_current("second")
        @test _API.model.get_entities()==[(0,8)]
        @test_throws ArgumentError _API.mesh.get()
        # `remove` drops the current slot and selects the last remaining model.
        _API.model.remove()
        @test _API.model.list()==[""]
        @test _API.model.get_current()==""
        @test size(_API.mesh.get().tris,2)>0
        _API.model.remove()
        @test_throws ArgumentError _API.model.remove()
        @test _API.model.list()==String[]
        @test_throws ArgumentError _API.model.get_entities()
        @test _API.model.add("fresh")===nothing
        @test _API.model.list()==["fresh"]
        @test _API.model.get_current()=="fresh"
        @test _API.model.set_current("fresh")===nothing
        @test_throws ArgumentError _API.model.set_current("missing")
        # Duplicate names are allowed; `set_current` selects the first match.
        @test _API.model.add("fresh")===nothing
        @test _API.model.list()==["fresh","fresh"]
        @test _API.model.set_current("fresh")===nothing
        @test _API.model.add_point(1,1,1;tag=7)==7
        @test _API.model.get_entities()==[(0,7)]
        _API.model.set_file_name("fresh.msh")
        @test _API.model.get_file_name()=="fresh.msh"
        metadata_model=Docs.meta(Tessella.API.model)
        for name in (:get_current,:set_current,:get_file_name,:set_file_name,
                     :add,:remove,:list,:set_visibility_per_window)
            @test haskey(metadata_model,
                         Docs.Binding(Tessella.API.model,name))
        end
    finally
        _API.finalize()
    end
    @test_throws ArgumentError _API.mesh.compute_renumbering()
end

@testset "volume attribute consumption and empty-cache parity" begin
    _API.initialize()
    try
        # Transfinite volume: an explicit cube with 12 transfinite edges and 6
        # transfinite faces fills through mesh_transfinite_volume.
        m=_API.model
        for (x,y,z) in [(0,0,0),(1,0,0),(1,1,0),(0,1,0),
                        (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
            m.add_point(x,y,z)
        end
        for (a,b) in [(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
                      (1,5),(2,6),(3,7),(4,8)]
            m.add_line(a,b)
        end
        faces=[[1,2,3,4],[5,6,7,8],[1,10,-5,-9],
               [2,11,-6,-10],[3,12,-7,-11],[4,9,-8,-12]]
        loop_tags=[m.add_curve_loop(f) for f in faces]
        surf_tags=[m.add_plane_surface([l]) for l in loop_tags]
        shell=m.add_surface_loop(surf_tags)
        m.add_volume([shell])
        for c in 1:12
            _API.mesh.set_transfinite_curve(c,4)
        end
        for s in surf_tags
            _API.mesh.set_transfinite_surface(s)
        end
        _API.mesh.set_transfinite_volume(1)
        _API.mesh.generate(3)
        node_tags,_=_API.mesh.get_nodes()
        element_types,element_tags,_=_API.mesh.get_elements()
        @test length(node_tags)==64        # 4x4x4 structured grid
        @test element_types==Int32[4]
        @test length(element_tags[1])==162 # 3x3x3 cells x 6 tets
        # A mismatched edge family fails explicitly.
        _API.mesh.set_transfinite_curve(1,6)
        @test_throws ArgumentError _API.mesh.generate(3)
    finally
        _API.finalize()
    end
    _API.initialize()
    try
        # set_reverse on a volume flips every tetrahedron; classification and
        # entity queries keep working on the intentionally inverted complex.
        _API.model.add_box(0,0,0,1,1,1)
        _API.mesh.set_reverse(3,1)
        _API.mesh.generate(3)
        node_tags,node_coords,_=_API.mesh.get_nodes()
        index=Dict(t=>i for (i,t) in enumerate(node_tags))
        coords=reshape(node_coords,3,:)
        _,_,blocks=_API.mesh.get_elements()
        flat=blocks[1]
        negative=count(1:div(length(flat),4)) do cell
            a,b,c,d=(coords[:,index[flat[4*(cell-1)+slot]]] for slot in 1:4)
            dot(c-a,cross(b-a,d-a))>0
        end
        # dot(c-a, cross(b-a,d-a)) = -det[b-a,c-a,d-a]: >0 means kernel-
        # negative, so a fully reversed volume reports every tet negative.
        @test negative==div(length(flat),4)
        element_types,_,_=_API.mesh.get_elements(3,1)
        @test element_types==Int32[4]
        # set_smoothing runs boundary-preserving Laplacian iterations.
        _API.mesh.clear()
        _API.mesh.set_reverse(3,1,false)
        _API.mesh.set_smoothing(3,1,2)
        _API.mesh.generate(3)
        @test validate(_API.mesh.get()).ok
        # A volume size callback densifies the interior.
        _API.mesh.set_size_callback((dim,tag,x,y,z,lc)->0.15)
        _API.mesh.generate(3)
        callback_nodes,_=_API.mesh.get_nodes()
        @test length(callback_nodes)>100
        _API.mesh.set_size_callback(nothing)
        # optimize accepts entity scopes and the 2-D method aliases; unknown
        # entities and unimplemented methods still fail.
        _API.mesh.optimize("",false,1,[(3,1)])
        @test validate(_API.mesh.get()).ok
        @test_throws ArgumentError _API.mesh.optimize("",false,1,[(3,99)])
        @test_throws ArgumentError _API.mesh.optimize("Netgen")
        # After clear(), queries answer empty arrays instead of throwing —
        # Gmsh parity for an unmeshed-but-known model.
        _API.mesh.clear()
        @test _API.mesh.get_nodes()[1]==UInt64[]
        @test _API.mesh.get_elements()[1]==Int32[]
        @test _API.mesh.get_element_types()==Int32[]
        @test _API.mesh.get_max_node_tag()==0
        @test _API.mesh.get_max_element_tag()==0
        @test _API.mesh.get_nodes(3,1)[1]==UInt64[]
        @test_throws ArgumentError _API.mesh.get_nodes(3,99)
    finally
        _API.finalize()
    end
    _API.initialize()
    try
        # Boundary Point sizes propagate into an explicit volume's interior
        # (Gmsh MeshSizeFromBoundary semantics for dim-3 generation).
        m=_API.model
        for (x,y,z) in [(0,0,0),(1,0,0),(1,1,0),(0,1,0),
                        (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
            m.add_point(x,y,z)
        end
        for (a,b) in [(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
                      (1,5),(2,6),(3,7),(4,8)]
            m.add_line(a,b)
        end
        faces=[[1,2,3,4],[5,6,7,8],[1,10,-5,-9],
               [2,11,-6,-10],[3,12,-7,-11],[4,9,-8,-12]]
        loop_tags=[m.add_curve_loop(f) for f in faces]
        surf_tags=[m.add_plane_surface([l]) for l in loop_tags]
        shell=m.add_surface_loop(surf_tags)
        m.add_volume([shell])
        for i in 1:8
            _API.mesh.set_size([(0,i)],0.08)
        end
        _API.mesh.generate(3)
        node_tags,_=_API.mesh.get_nodes()
        @test length(node_tags)>1000
        @test validate(_API.mesh.get()).ok
    finally
        _API.finalize()
    end
    _API.initialize()
    try
        # A 3-sided transfinite surface defaults to Gmsh's legacy
        # Mesh.TransfiniteTri=0 collapsed-quadrilateral algorithm: n=4
        # divisions give 1+n(n+1)=21 nodes and n(2n-1)=28 triangles.
        m=_API.model
        pa=m.add_point(0,0,0);pb=m.add_point(1,0,0);pc=m.add_point(0,1,0)
        ea=m.add_line(pa,pb);eb=m.add_line(pb,pc);ec=m.add_line(pc,pa)
        tri_loop=m.add_curve_loop([ea,eb,ec])
        tri_face=m.add_plane_surface([tri_loop])
        for e in (ea,eb,ec)
            _API.mesh.set_transfinite_curve(e,5)
        end
        _API.mesh.set_transfinite_surface(tri_face)
        _API.mesh.generate(2)
        all_nodes,all_coords,_=_API.mesh.get_nodes()
        @test length(all_nodes)==21
        lattice=sort!([(all_coords[3i-2],all_coords[3i-1]) for i in 1:21])
        @test lattice==sort!(vcat([(0.0,0.0)],
            [(i*(4-j)/16,i*j/16) for i in 1:4 for j in 0:4]))
        # Boundary nodes classify on the curves/points, so the surface entity
        # reports only its interior nodes — matching Gmsh getNodes(2,tag).
        tri_nodes,_=_API.mesh.get_nodes(2,tri_face)
        @test length(tri_nodes)==9
        _,_,tri_blocks=_API.mesh.get_elements(2,tri_face)
        @test div(length(tri_blocks[1]),3)==28
        @test validate(_API.mesh.get()).ok
        # The collapsed algorithm accepts unequal counts when the two sides
        # incident to the rotated collapsed corner match: (5,5,6) succeeds.
        _API.mesh.set_transfinite_curve(ec,6)
        _API.mesh.generate(2)
        un_nodes,_=_API.mesh.get_nodes()
        @test length(un_nodes)==25
        # (5,6,7) has no valid collapsed corner and still fails explicitly.
        _API.mesh.set_transfinite_curve(eb,7)
        @test_throws ArgumentError _API.mesh.generate(2)
        _API.mesh.set_transfinite_curve(eb,5)
        _API.mesh.set_transfinite_curve(ec,5)
        # Mesh.TransfiniteTri=1 selects the compact triangular-lattice
        # algorithm: n=4 divisions give (n+1)(n+2)/2 nodes and n^2 triangles.
        @test _API.option("Mesh.TransfiniteTri")==0.0
        _API.option("Mesh.TransfiniteTri",1)
        @test _API.option("Mesh.TransfiniteTri")==1.0
        @test_throws ArgumentError _API.option("Mesh.TransfiniteTri",2)
        _API.mesh.generate(2)
        all_nodes,all_coords,_=_API.mesh.get_nodes()
        @test length(all_nodes)==15
        lattice=sort!([(all_coords[3i-2],all_coords[3i-1]) for i in 1:15])
        @test lattice==sort!([(i/4,j/4) for i in 0:4 for j in 0:4-i])
        tri_nodes,_=_API.mesh.get_nodes(2,tri_face)
        @test length(tri_nodes)==3
        _,_,tri_blocks=_API.mesh.get_elements(2,tri_face)
        @test div(length(tri_blocks[1]),3)==16
        @test validate(_API.mesh.get()).ok
        # Mismatched transfinite curve counts fail explicitly under Tri=1.
        _API.mesh.set_transfinite_curve(ec,6)
        @test_throws ArgumentError _API.mesh.generate(2)
        _API.mesh.set_transfinite_curve(ec,5)
        # A loop built from reversed curve signs still meshes.
        _API.mesh.clear()
        _API.model.remove_entities([(2,tri_face)])
        rev_loop=m.add_curve_loop([-ec,-eb,-ea])
        rev_face=m.add_plane_surface([rev_loop])
        _API.mesh.set_transfinite_surface(rev_face)
        _API.mesh.generate(2)
        rev_nodes,_=_API.mesh.get_nodes(2,rev_face)
        @test length(rev_nodes)==3
        # Explicitly pinned corners select the collapsed corner under the
        # default Tri=0 algorithm: pinning (pc,pa,pb) collapses at
        # pc=(0,1,0), unreachable by auto-detection.
        _API.option("Mesh.TransfiniteTri",0)
        _API.mesh.clear()
        _API.model.remove_entities([(2,rev_face)])
        pin_loop=m.add_curve_loop([ea,eb,ec])
        pin_face=m.add_plane_surface([pin_loop])
        _API.mesh.set_transfinite_surface(pin_face,"Left",[pc,pa,pb])
        _API.mesh.generate(2)
        pin_nodes,pin_coords,_=_API.mesh.get_nodes()
        @test length(pin_nodes)==21
        pin_lattice=sort!([(pin_coords[3i-2],pin_coords[3i-1]) for i in 1:21])
        @test pin_lattice==sort!(vcat([(0.0,1.0)],
            [(i*j/16,(4-i)/4) for i in 1:4 for j in 0:4]))
        @test validate(_API.mesh.get()).ok
        # Corner tags that are not the surface junctions fail explicitly.
        free_point=m.add_point(0.25,0.25,0)
        _API.mesh.set_transfinite_surface(pin_face,"Left",[pa,pb,free_point])
        @test_throws ArgumentError _API.mesh.generate(2)
    finally
        _API.finalize()
    end
end

@testset "API option Mesh.FlexibleTransfinite plumbing" begin
    _API.finalize()
    try
        _API.initialize()
        # Option surface: defaults, truncation, recombination clamp, and the
        # `MeshSizeFactor`/`CharacteristicLengthFactor` single-field alias.
        @test _API.option("Mesh.FlexibleTransfinite")==0.0
        @test _API.option("Mesh.RecombineAll")==0.0
        @test _API.option("Mesh.RecombinationAlgorithm")==1.0
        @test _API.option("Mesh.CharacteristicLengthFactor")==1.0
        @test _API.option("Mesh.FlexibleTransfinite",1)==1.0
        @test _API.option("Mesh.FlexibleTransfinite")==1.0
        @test _API.option("Mesh.RecombinationAlgorithm",9)==0.0
        @test _API.option("Mesh.CharacteristicLengthFactor",2)==2.0
        @test _API.option("Mesh.MeshSizeFactor")==2.0
        @test _API.option("Mesh.MeshSizeFactor",3)==3.0
        @test _API.option("Mesh.CharacteristicLengthFactor")==3.0
        @test_throws ArgumentError _API.option("Mesh.CharacteristicLengthFactor",0.0)

        m=_API.model
        m.add_point(0,0,0);m.add_point(1,0,0);m.add_point(1,1,0);m.add_point(0,1,0)
        m.add_line(1,2);m.add_line(2,3);m.add_line(3,4);m.add_line(4,1)
        loop=m.add_curve_loop([1,2,3,4])
        m.add_plane_surface([loop])
        for c in 1:4
            _API.mesh.set_transfinite_curve(c,11)
        end
        _API.mesh.set_transfinite_surface(1)
        # FlexibleTransfinite=1 with lcFactor=3: 11 -> 3 nodes per side on
        # the structured grid.
        _API.mesh.generate(2)
        nodes,_,_=_API.mesh.get_nodes()
        @test length(nodes)==9
        _API.option("Mesh.FlexibleTransfinite",0)
        _API.mesh.generate(2)
        nodes,_,_=_API.mesh.get_nodes()
        @test length(nodes)==121
    finally
        _API.finalize()
    end
end
