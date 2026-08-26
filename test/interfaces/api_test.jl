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
            Periodic Curve {2} = {4} Translate {1, 0, 0};
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
