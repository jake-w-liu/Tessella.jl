module GeoHigherDefaultTests
using Test,Tessella,SHA,TOML
using Tessella.MeshTypes: Mesh,mesh_crc
const API=Tessella.API
const GEO=Tessella.GeoExec

@testset "Ordinary GEO generation publishes literal complete owner products" begin
    path=joinpath(@__DIR__,"../artifacts/onlyempty_geo_default_square_primary.toml")
    @test bytes2hex(sha256(replace(read(path,String),"\r\n"=>"\n")))==
        "f4393385a11e8211886f1f4979e73c1e6b26c82eea3659a16edd4a26e2e87640"
    primary=TOML.parsefile(path)
    @test primary["gmsh_version"]=="4.15.2"
    @test primary["primary_json_sha256"]=="edc61727b74b987d1ea4bacf8e771524208945696ae03b0e4a53768b8979b0cb"
    current=API.CURRENT[];cache=API.LAST_MESH[];classification=API.LAST_MESH_CLASS[]
    for row in primary["rows"]
        mktempdir() do directory
            source=joinpath(directory,"default.geo")
            write(source,"""
            Mesh.Renumber=$(row["renumber"]); Mesh.ElementOrder=$(row["order"]);
            Mesh.MeshOnlyEmpty=0;
            Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
            Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
            Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
            Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
            Transfinite Curve{1:4}=2; Transfinite Surface{1}; Mesh 2;
            """)
            result=GEO.execute_geo(source)
            authority=result.model.meshing.mesh_authority
            actual=authority.mesh;data=actual.entity_data
            @test result.mesh isa Mesh
            @test authority.carrier===result.mesh
            @test actual.coords[:,authority.carrier_rows]==result.mesh.coords
            @test size(result.mesh.segs,2)==4 && size(result.mesh.tris,2)==2
            @test data.external_node_tags==row["all_node_tags"]
            @test Dict(string(block.msh)=>size(block.nodes,2) for block in actual.blocks)==row["counts"]
            @test allunique(data.external_node_tags)
            @test allunique(vcat(data.external_element_tags...))
            for entity in row["entities"]
                dimension,tag=entity["owner"]
                part=GEO.geo_entity_mesh(result,dimension,tag)
                metadata=part.entity_data
                owned=Int[position for (position,owner) in enumerate(metadata.node_entities)
                          if owner==(dimension,Int32(tag))]
                @test metadata.external_node_tags[owned]==entity["nodes"]
                @test vec(part.coords[:,owned])≈Float64.(entity["xyz"]) atol=2e-11 rtol=0
                parameters=reduce(vcat,(something(metadata.node_parametric[position],Float64[])
                                       for position in owned);init=Float64[])
                @test parameters≈Float64.(entity["parameters"]) atol=2e-11 rtol=0
                @test Int[block.msh for block in part.blocks]==entity["types"]
                @test metadata.external_element_tags==entity["elements"]
                @test [vec(metadata.external_node_tags[block.nodes]) for block in part.blocks]==entity["connections"]
            end
            @test API.CURRENT[]===current && API.LAST_MESH[]===cache && API.LAST_MESH_CLASS[]===classification
        end
    end
end

function prepare_point_vertices()
    API.initialize()
    for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        API.model.add_point(x,y,0.;tag,meshSize=1.)
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        API.model.add_line(first,last;tag)
        API.mesh.set_transfinite_curve(tag,3)
    end
    API.model.add_curve_loop([1,2,3,4];tag=1)
    API.model.add_plane_surface([1];tag=1)
    API.mesh.set_transfinite_surface(1)
    API.mesh.add_nodes(0,1,[101,102,103],zeros(9))
end

@testset "GEO native carrier binds actual coincident Point vertices after order" begin
    # Independent Gmsh4.15.2 Point rows used by the API companion test:
    # 75161be150ba8924b0a1910568adbba199309f8fcf02366420bee9d10321b141.
    # P1 retains all three rows; P2 retains referenced rows and reclassifies
    # the first Curve endpoint. A coordinate cannot stand in for its UID.
    for order in (1,2),renumber in (0,1),supplied in (false,true)
        try
            prepare_point_vertices()
            supplied && API.mesh.add_elements_by_type(1,15,[777,778],[101,102])
            model=Tessella.Model._model_mesh_transaction_copy(API.CURRENT[])
            before=(API.CURRENT[],API.LAST_MESH[],API.LAST_MESH_CLASS[])
            context=Tessella.IO._GeoNumericContext()
            context.option_numbers[("Mesh",0,"ElementOrder")]=Float64(order)
            context.option_numbers[("Mesh",0,"Renumber")]=Float64(renumber)
            context.option_numbers[("Mesh",0,"MeshOnlyEmpty")]=0.
            carrier,_=GEO._geo_mesh_model(model,2,context)
            authority=model.meshing.mesh_authority
            point=Tessella.Model._model_mesh_product(authority,0,1)
            metadata=point.entity_data
            owned=findall(==((0,Int32(1))),metadata.node_entities)
            expected=order==1 ? (renumber==0 ? [101,102,103] : [1,2,3]) :
                supplied ? (renumber==0 ? [101,102] : [1,2]) :
                (renumber==0 ? [103] : [1])
            @test metadata.external_node_tags[owned]==UInt64.(expected)
            @test point.coords[:,owned]==zeros(3,length(expected))
            cell=only(point.blocks)
            expected_cells=supplied ? (renumber==0 ? [777,778] : [1,2]) : [1]
            expected_connections=supplied ? (renumber==0 ? [101,102] : [1,2]) :
                (renumber==0 ? [103] : [order==1 ? 3 : 1])
            @test only(metadata.external_element_tags)==UInt64.(expected_cells)
            @test vec(metadata.external_node_tags[cell.nodes])==UInt64.(expected_connections)
            @test carrier isa Mesh && authority.carrier===carrier
            @test authority.mesh.coords[:,authority.carrier_rows]==carrier.coords
            @test allunique(authority.carrier_rows)
            @test (API.CURRENT[],API.LAST_MESH[],API.LAST_MESH_CLASS[])===before
        finally
            API.finalize()
        end
    end
end
end
