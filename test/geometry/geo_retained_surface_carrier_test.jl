module GeoRetainedSurfaceCarrierTests
using Tessella,Test,LinearAlgebra
const M=Tessella.Model
const E=Tessella.Elements

function scene(points,only_empty,phases)
    return """
    SetFactory("OpenCASCADE");
    Box(1)={0,0,0,1,1,1};
    Point(11)={3,0,0,1};Point(12)={4,0,0,1};Point(13)={4,1,0,1};Point(14)={3,1,0,1};
    Line(21)={11,12};Line(22)={12,13};Line(23)={13,14};Line(24)={14,11};
    Curve Loop(7)={21,22,23,24};Plane Surface(7)={7};
    Transfinite Curve{1:12,21:24}=$points;Transfinite Surface{1:7};Transfinite Volume{1};
    Mesh.Renumber=0;Mesh.MeshOnlyEmpty=$only_empty;$phases
    """
end

function surface_state(mesh)
    data=mesh.entity_data
    return (nodes=copy(data.external_node_tags),owners=copy(data.node_entities),
        coordinates=copy(mesh.coords),parameters=deepcopy(data.node_parametric),
        types=[block.msh for block in mesh.blocks],
        corners=[reshape(data.external_node_tags[vec(block.nodes)],size(block.nodes)) for block in mesh.blocks])
end

function signed_volume(mesh)
    return sum(det(hcat(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
        mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]],
        mesh.coords[:,cell[4]]-mesh.coords[:,cell[1]]))/6
        for block in mesh.blocks for cell in eachcol(block.nodes);init=0.)
end

@testset "GEO triangular geometry carrier retains disjoint actual Surface products" begin
    for points in (2,3),only_empty in (0,1)
        @testset "TF$points-OnlyEmpty$only_empty" begin
            mktempdir() do directory
                executions=Any[]
                for (index,phases) in enumerate(("Mesh 2;","Mesh 2;Mesh 3;","Mesh 2;Mesh 3;Mesh 3;"))
                    path=joinpath(directory,"retained$index.geo")
                    write(path,scene(points,only_empty,phases))
                    push!(executions,Tessella.execute_geo(path))
                end
                products=[M._model_mesh_product(execution.model.meshing.mesh_authority,2,7)
                    for execution in executions]
                @test all(mesh->mesh isa E.MixedMesh,products)
                before=surface_state(products[1])
                for (execution,product) in zip(executions[2:end],products[2:end])
                    @test execution.mesh isa (only_empty==0 ? Tessella.MeshTypes.Mesh : E.MixedMesh)
                    @test isequal(surface_state(product),before)
                    @test sum(size(block.nodes,2) for block in product.blocks)==2*(points-1)^2
                    volume=M._model_mesh_product(execution.model.meshing.mesh_authority,3,1)
                    @test volume isa E.MixedMesh
                    @test all(block->block.msh==4,volume.blocks)
                    @test signed_volume(volume)≈1.
                    @test all(det(hcat(volume.coords[:,cell[2]]-volume.coords[:,cell[1]],
                        volume.coords[:,cell[3]]-volume.coords[:,cell[1]],
                        volume.coords[:,cell[4]]-volume.coords[:,cell[1]]))>0
                        for block in volume.blocks for cell in eachcol(block.nodes))
                    actual=execution.model.meshing.mesh_authority.mesh
                    if only_empty==0
                        @test execution.mesh.coords!==actual.coords
                        snapshot=copy(actual.coords)
                        execution.mesh.coords[1,1]+=0.125
                        @test isequal(actual.coords,snapshot)
                    end
                end
            end
        end
    end
end
end
