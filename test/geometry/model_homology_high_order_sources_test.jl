module HomologySourceOrderProbe
using Tessella,Test
include("../../validation/homology_native/high_order_source_primary_v1.jl")
const API=Tessella.API
function check_stage(expected)
    @test API.model.get_entities()==[(e["dimension"],e["tag"]) for e in expected["entities"]]
    for entity in expected["entities"]
        dim,tag=entity["dimension"],entity["tag"]
        nodes,coordinates,parameters=API.mesh.get_nodes(dim,tag,false,true)
        @test nodes==entity["node_tags"]
        @test coordinates==entity["coordinates"]
        @test parameters==entity["parameters"]
        types,cells,connectivity=API.mesh.get_elements(dim,tag)
        @test types==entity["types"]
        @test cells==entity["element_tags"]
        @test connectivity==entity["connectivity"]
    end
    for node in expected["nodes"]
        coordinates,parameters,dimension,owner=API.mesh.get_node(node["tag"])
        @test coordinates==node["coordinates"]
        @test parameters==node["parameters"]
        @test (dimension,owner)==(node["dimension"],node["owner"])
    end
    @test API.model.get_physical_groups()==[(p["dimension"],p["tag"]) for p in expected["physical"]]
    for physical in expected["physical"]
        dim,tag=physical["dimension"],physical["tag"]
        @test API.model.get_entities_for_physical_group(dim,tag)==physical["members"]
        @test API.model.get_physical_name(dim,tag)==physical["name"]
    end
end
@testset "Actual higher-order source identity during integral chain publication" begin
    for row in HOMOLOGY_SOURCE_ORDER_ROWS
        @testset "$(row["name"]) $(row["kind"]) edited=$(row["edited"])" begin
            API.initialize()
            try
                dim=row["dimension"]
                API.option("Mesh.Renumber",0)
                API.model.add_discrete_entity(dim,1)
                API.mesh.add_nodes(dim,1,row["input_nodes"],row["input_coordinates"],row["input_parameters"])
                API.mesh.add_elements_by_type(1,row["msh"],[1001],row["input_nodes"])
                API.model.add_physical_group(dim,[1];tag=10)
                check_stage(row["baseline"])
                points=reshape(row["input_coordinates"],3,:)
                coordinates=Dict(Int32(tag)=>ntuple(axis->points[axis,slot],3)
                    for (slot,tag) in enumerate(row["input_nodes"]))
                volume=Tessella.Model._native_homology_measure(row["msh"],Int32.(row["input_nodes"]),coordinates)
                @test isapprox(volume,row["volume"];rtol=0,atol=2e-11)
                API.mesh.add_homology_request(row["kind"],[10],[],collect(0:dim))
                @test API.mesh.compute_homology()==[Tuple(pair) for pair in row["returned"]]
                check_stage(row["immediate"])
                API.mesh.reclassify_nodes()
                check_stage(row["reclassified"])
            finally
                API.finalize()
            end
        end
    end
end
@testset "Higher-order raw association without basis or map reconstruction" begin
    for row in HOMOLOGY_SOURCE_ASSOCIATION_ROWS
        @testset "$(row["name"]) edited=$(row["edited"]) $(row["kind"]) source" begin
            API.initialize()
            try
                dim=row["dimension"]
                API.option("Mesh.Renumber",0)
                API.model.add_discrete_entity(dim,1)
                API.mesh.add_nodes(dim,1,row["input_nodes"],row["input_coordinates"],row["input_parameters"])
                API.mesh.add_elements_by_type(1,row["msh"],[1001],row["input_nodes"])
                API.model.add_physical_group(dim,[1];tag=10)
                check_stage(row["baseline"])
                API.mesh.reclassify_nodes()
                check_stage(row["association"])
                API.mesh.reclassify_nodes()
                check_stage(row["association"])
            finally
                API.finalize()
            end
        end
    end
end
end
