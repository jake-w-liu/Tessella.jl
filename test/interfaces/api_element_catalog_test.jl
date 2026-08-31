using Test
using SHA
using Tessella
using Tessella.MeshTypes: mesh_crc

const _CATALOG_API=Tessella.API
const _CATALOG_FAMILY_NAMES=Dict(
    :pnt=>"Point",:lin=>"Line",:tri=>"Triangle",:qua=>"Quadrangle",
    :tet=>"Tetrahedron",:hex=>"Hexahedron",:pri=>"Prism",
    :pyr=>"Pyramid",:trih=>"Trihedron")
const _CATALOG_SPECIAL_PROPERTIES=(34,35,69,133,134,135,136)

function _catalog_write_i32!(stream,value)
    write(stream,htol(reinterpret(UInt32,Int32(value))))
end

function _catalog_api_sha()
    stream=IOBuffer()
    tags=vcat(sort!(collect(keys(MSH_CATALOG))),
              collect(_CATALOG_SPECIAL_PROPERTIES))
    for tag in tags
        name,dimension,order,num_nodes,coordinates,num_primary_nodes=
            _CATALOG_API.mesh.get_element_properties(tag)
        write(stream,htol(UInt64(tag)))
        bytes=codeunits(name)
        write(stream,htol(UInt64(length(bytes))))
        write(stream,bytes)
        for value in (dimension,order,num_nodes,num_primary_nodes)
            _catalog_write_i32!(stream,value)
        end
        write(stream,htol(UInt64(length(coordinates))))
        for value in coordinates
            write(stream,htol(reinterpret(UInt64,value)))
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

@testset "element catalog queries through API" begin
    _CATALOG_API.finalize()

    for tag in sort!(collect(keys(MSH_CATALOG)))
        spec=msh_spec(tag)
        family_name=_CATALOG_FAMILY_NAMES[spec.family]
        @test _CATALOG_API.mesh.get_element_type(
            family_name,spec.order,spec.serendipity)==Int32(tag)

        direct=msh_properties(tag)
        queried=_CATALOG_API.mesh.get_element_properties(tag)
        @test queried==(
            direct.name,Int32(direct.dim),Int32(direct.order),
            Int32(direct.num_nodes),direct.local_node_coordinates,
            Int32(direct.num_primary_nodes))
    end

    @test _CATALOG_API.mesh.get_element_type("triangle",2,true)==Int32(9)
    @test _CATALOG_API.mesh.get_element_type("LINE",10,true)==Int32(66)
    @test _CATALOG_API.mesh.get_element_type("pyramid",1,true)==Int32(7)
    @test _CATALOG_API.mesh.get_element_type("Point",-1,true)==Int32(15)
    @test _CATALOG_API.mesh.get_element_type(
        "Trihedron",typemax(Int32),false)==Int32(140)

    expected_special=Dict(
        34=>("Polygon",2),35=>("Polyhedron",3),
        69=>("Polygon Border",2),133=>("Point Xfem",0),
        134=>("Line Xfem",1),135=>("Triangle Xfem",2),
        136=>("Tetrahedron Xfem",3))
    for (tag,(name,dimension)) in expected_special
        @test _CATALOG_API.mesh.get_element_properties(tag)==
              (name,Int32(dimension),Int32(1),Int32(0),Float64[],Int32(0))
    end

    @test _catalog_api_sha()==
          "4c8e4086febda2bc3f482087d589eb157dbe57394f30781449a17a64e1b2a65c"

    detached=_CATALOG_API.mesh.get_element_properties(4)
    detached[5][1]=99
    @test _CATALOG_API.mesh.get_element_properties(4)[5][1]==0

    for call in (
        ()->_CATALOG_API.mesh.get_element_type(:Line,1),
        ()->_CATALOG_API.mesh.get_element_type("",1),
        ()->_CATALOG_API.mesh.get_element_type("Quadrilateral",1),
        ()->_CATALOG_API.mesh.get_element_type("Line\0",1),
        ()->_CATALOG_API.mesh.get_element_type("Line",true),
        ()->_CATALOG_API.mesh.get_element_type("Triangle",-1),
        ()->_CATALOG_API.mesh.get_element_type("Triangle",11),
        ()->_CATALOG_API.mesh.get_element_type(
            "Triangle",big(typemax(Int32))+1),
        ()->_CATALOG_API.mesh.get_element_type("Triangle",1,0),
        ()->_CATALOG_API.mesh.get_element_properties(true),
        ()->_CATALOG_API.mesh.get_element_properties(0),
        ()->_CATALOG_API.mesh.get_element_properties(67),
        ()->_CATALOG_API.mesh.get_element_properties(138),
        ()->_CATALOG_API.mesh.get_element_properties(999),
    )
        @test_throws ArgumentError call()
    end

    try
        _CATALOG_API.initialize()
        _CATALOG_API.model.add_box(0,0,0,1,1,1;tag=1)
        generated=_CATALOG_API.mesh.generate(3)
        expected_crc=mesh_crc(generated)
        @test _CATALOG_API.mesh.get_element_type("Tetrahedron",1)==Int32(4)
        @test _CATALOG_API.mesh.get_element_properties(4)[1]=="Tetrahedron 4"
        @test mesh_crc(_CATALOG_API.mesh.get())==expected_crc
    finally
        _CATALOG_API.finalize()
    end
end
