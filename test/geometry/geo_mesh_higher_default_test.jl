module GeoHigherDefaultTests
using Test,Tessella
using Tessella.MeshTypes: Mesh
const API=Tessella.API
const GEO=Tessella.GeoExec

function prepare_point_vertices(;curve_nodes=3,point_vertices=true)
    API.initialize()
    for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        API.model.add_point(x,y,0.;tag,meshSize=1.)
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        API.model.add_line(first,last;tag)
        API.mesh.set_transfinite_curve(tag,curve_nodes)
    end
    API.model.add_curve_loop([1,2,3,4];tag=1)
    API.model.add_plane_surface([1];tag=1)
    API.mesh.set_transfinite_surface(1)
    point_vertices && API.mesh.add_nodes(0,1,[101,102,103],zeros(9))
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

@testset "Surface producer rows retain Point-owned constructor corners" begin
    for order in (1,2)
        try
            prepare_point_vertices(;curve_nodes=2,point_vertices=false)
            model=Tessella.Model._model_mesh_transaction_copy(API.CURRENT[])
            context=Tessella.IO._GeoNumericContext()
            context.option_numbers[("Mesh",0,"ElementOrder")]=Float64(order)
            context.option_numbers[("Mesh",0,"Renumber")]=0.
            GEO._geo_mesh_model(model,2,context)
            authority=model.meshing.mesh_authority
            surface=(2,Int32(1))
            @test haskey(authority.producer_rows,surface)
            @test all(owner->owner[1]!=0,keys(authority.producer_rows))
            rows=authority.producer_rows[surface]
            @test length(rows)==4
            @test Set(authority.mesh.entity_data.node_entities[rows])==
                Set((0,Int32(tag)) for tag in 1:4)
            for tag in 1:4
                point=Tessella.Model._model_mesh_product(authority,0,tag;producer_order=true)
                @test size(point.coords,2)==1
                @test point.entity_data.node_entities==[(0,Int32(tag))]
            end
        finally
            API.finalize()
        end
    end
end
end
