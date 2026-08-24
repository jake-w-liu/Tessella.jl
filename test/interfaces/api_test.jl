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
        write(valid,"Box(1) = {0, 0, 0, 1, 1, 1};\n")
        write(invalid,"Extrude {0, 0, 1} { Volume{1}; }\n")

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
        finally
            _API.finalize()
        end
    end

    @test _API.finalize()===nothing
    @test _API.finalize()===nothing
    @test isempty(Docs.undocumented_names(Tessella.API;private=false))
end
