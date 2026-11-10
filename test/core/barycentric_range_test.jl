module BarycentricRangeChecks
using Test,Tessella
const E=Tessella.Elements
const M=Tessella.MeshTypes

function classified(mesh;copies=1)
    mixed=E.simplex_to_mixed(mesh);dim=isempty(mesh.tets) ? 2 : 3
    count=size(mixed.coords,2);block=mixed.blocks[1]
    cells=repeat(block.nodes,1,copies)
    bounds=(minimum(mixed.coords[1,:]),minimum(mixed.coords[2,:]),minimum(mixed.coords[3,:]),
            maximum(mixed.coords[1,:]),maximum(mixed.coords[2,:]),maximum(mixed.coords[3,:]))
    entity=E.MixedEntity(dim,1,bounds;physical_tags=[7])
    parameters=[fill(1e308,dim) for _ in 1:count]
    data=E.MixedEntityData(Dict((dim,1)=>entity);
        node_entities=fill((dim,Int32(1)),count),node_parametric=parameters,
        external_node_tags=UInt64.(1:count),block_entities=[fill(Int32(1),copies)],
        external_element_tags=[UInt64.(1:copies)])
    return E.MixedMesh(mixed.coords,[E.ElementBlock(block.msh,cells,fill(Int32(7),copies))];entity_data=data)
end

function expected_center(coords)
    return setprecision(BigFloat,4096) do
        [Float64(sum(BigFloat,coords[d,:])/size(coords,2)) for d in 1:3]
    end
end

@testset "Barycentric subdivision retains finite extreme means" begin
    for degree in (3,4),xs in ((8e307,1e308,8e307,8e307),
                              (floatmax(Float64),floatmax(Float64),-floatmax(Float64),-floatmax(Float64)),
                              (-8e307,-1e308,-8e307,-8e307))
        coords=zeros(3,degree);coords[1,:].=xs[1:degree]
        coords[2,2]=1.;coords[3,3]=1.
        degree==4 && (coords[2,4]=1.;coords[3,4]=2.)
        nodes=reshape(Int32.(1:degree),degree,1)
        mesh=degree==3 ? M.Mesh(coords;tris=nodes,tri_tag=[7]) : M.Mesh(coords;tets=nodes,tet_tag=[7])
        for source in (mesh,E.simplex_to_mixed(mesh),classified(mesh))
            original=copy(source.coords)
            refined=E.barycentric_refine(source)
            @test all(isfinite,refined.coords)
            @test refined.coords[:,end]==expected_center(coords)
            @test refined.coords[:,1:degree]==coords
            @test source.coords==original
            if refined isa E.MixedMesh
                @test size(refined.blocks[1].nodes,2)==degree
                @test all(==(7),refined.blocks[1].tags)
                if refined.entity_data!==nothing
                    @test refined.entity_data.node_parametric[end]==fill(1e308,degree==3 ? 2 : 3)
                    @test refined.entity_data.node_entities[end]==(degree==3 ? 2 : 3,Int32(1))
                    @test refined.entity_data.external_node_tags[end]==degree+1
                end
            else
                @test degree==3 ? M.ntris(refined)==3 : M.ntets(refined)==4
            end
        end
    end
end
end
