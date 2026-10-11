module OriginalVertexTypeProbe
using Tessella,Test
include("../artifacts/original_vertex_type_primary_v1.jl")
const API=Tessella.API
function check_stage(expected,dimension)
    tags,coordinates,parameters=API.mesh.get_nodes(dimension,7,false,true)
    @test tags==expected["node_tags"]
    @test coordinates==expected["coordinates"]
    @test parameters==expected["parameters"]
    types,cells,connectivity=API.mesh.get_elements(dimension,7)
    @test types==expected["types"]
    @test cells==expected["element_tags"]
    @test connectivity==expected["connectivity"]
    for node in expected["nodes"]
        xyz,parameters,dim,owner=API.mesh.get_node(node["tag"])
        @test xyz==node["coordinates"]
        @test parameters==node["parameters"]
        @test (dim,owner)==(node["dimension"],node["owner"])
    end
end
@testset "Original vertex type preserves addNodes parameter semantics" begin
    for row in ORIGINAL_VERTEX_TYPE_ROWS
        dimension=row["dimension"]
        @testset "dimension$dimension supplied$(row["supplied"])" begin
            API.initialize()
            try
                API.option("Mesh.Renumber",0)
                API.model.add_discrete_entity(dimension,7)
                API.mesh.add_nodes(dimension,7,[11,12,13,14],
                    [0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],row["input_parameters"])
                kind=(15,1,2,4)[dimension+1]
                labels=dimension==0 ? [201,202,203,204] : dimension==1 ? [201,202] : [201]
                connectivity=dimension==2 ? [11,12,13] : [11,12,13,14]
                API.mesh.add_elements_by_type(7,kind,labels,connectivity)
                check_stage(row["stages"][1],dimension)
                API.mesh.set_node(11,[.02,.03,.04],[.9+.1*k for k in 0:dimension-1])
                check_stage(row["stages"][2],dimension)
                API.mesh.reclassify_nodes()
                check_stage(row["stages"][3],dimension)
            finally
                API.finalize()
            end
        end
    end
end
end
