module DuplicatePointRefineTests
using Test,Tessella
const API=Tessella.API

# Independent pinned Gmsh 4.15.2 result: duplicate Point cells 777/778
# retain separate identities on node1 through order2, order1 and refinement.
# Primary JSON SHA256 59dc9aa6ac07bc7d6029a02fb6d32b25bdb2f7ffc21d9fc6338b40cbe730a2b2.
function setup(;native=false)
    API.initialize()
    API.option("Mesh.Renumber",0)
    API.model.add_discrete_entity(0,91)
    API.model.add_discrete_entity(2,92)
    if native
        E=Tessella.Elements
        mesh=E.MixedMesh(Float64[3 0 1 0;3 0 0 1;0 0 0 0],
            [E.ElementBlock(15,reshape(Int32[1,1],1,2)),
             E.ElementBlock(2,reshape(Int32[2,3,4],3,1))])
        class=API._mixed_classification(mesh,(2,Int32(92)),
            [(0,Int32(91)),(2,Int32(92))],
            [(0,Int32(91)),(2,Int32(92)),(2,Int32(92)),(2,Int32(92))],
            Dict{Tuple{Int,Int32},Vector{Int32}}(),
            Dict(15=>Int32[91,91],2=>Int32[92]))
        table=API._cache_public_tags(mesh;node_tags=1:4,element_tags=[777,778,888])
        API._replace_mesh_cache_locked!(mesh,API._classification_with_public_tags(class,table))
        @test API.LAST_MESH[].entity_data===nothing
        return nothing
    end
    API.mesh.add_nodes(0,91,[1],[3.,3.,0.])
    API.mesh.add_nodes(2,92,[2,3,4],[0.,0,0,1,0,0,0,1,0])
    API.mesh.add_elements_by_type(91,15,[777,778],[1,1])
    API.mesh.add_elements_by_type(92,2,[888],[2,3,4])
    return nothing
end

@testset "Point exception preserves strict structural validation" begin
    E=Tessella.Elements
    point=E.MixedMesh(reshape([3.,3.,0.],3,1),
        [E.ElementBlock(15,reshape(Int32[1,1],1,2))])
    @test !E.validate(point).ok
    @test E.validate(point;reject_duplicate_point_cells=false).ok
    @test_throws ArgumentError E.validate(point;reject_duplicate_point_cells=1)
    duplicate=E.MixedMesh(Float64[0 1;0 0;0 0],
        [E.ElementBlock(1,Int32[1 2;2 1])])
    @test !E.validate(duplicate;reject_duplicate_point_cells=false).ok
    repeated=E.MixedMesh(Float64[0 1;0 0;0 0],
        [E.ElementBlock(1,reshape(Int32[1,1],2,1))])
    @test !E.validate(repeated;reject_duplicate_point_cells=false).ok
end

@testset "Distinct Point cells round-trip through supported MSH formats" begin
    E=Tessella.Elements
    try
        setup();API.mesh.set_order(1)
        mesh=API.LAST_MESH[]
        @test mesh.entity_data!==nothing
        mktempdir() do directory
            for version in (2.2,4.1),binary in (false,true)
                path=joinpath(directory,"point_$(version)_$(binary).msh")
                E.write_mixed_msh(path,mesh;version,binary)
                result=E.read_mixed_msh(path)
                points=[(index,block) for (index,block) in enumerate(result.blocks) if block.msh==15]
                @test length(points)==1
                index,block=only(points)
                @test size(block.nodes)==(1,2) && block.nodes[1,1]==block.nodes[1,2]
                @test result.coords[:,block.nodes[1,1]]==[3.,3.,0.]
                @test E.validate(result;reject_duplicate_point_cells=false).ok
                if version==4.1
                    @test result.entity_data.external_element_tags[index]==UInt64[777,778]
                    @test result.entity_data.node_entities[block.nodes[1,1]]==(0,Int32(91))
                end
            end
            duplicate=E.MixedMesh(Float64[0 1;0 0;0 0],
                [E.ElementBlock(1,Int32[1 2;2 1])])
            invalid=joinpath(directory,"invalid.msh")
            write(invalid,"existing file")
            @test_throws r"duplicate type 1" E.write_mixed_msh(invalid,duplicate)
            @test read(invalid,String)=="existing file"
        end
    finally
        API.finalize()
    end
end

function check_points()
    @test API.mesh.get_elements_by_type(15,91)==(UInt64[777,778],UInt64[1,1])
    @test API.mesh.get_nodes(0,91)==(UInt64[1],[3.,3.,0.],Float64[])
    @test API.mesh.get_element(777)==(15,UInt64[1],0,91)
    @test API.mesh.get_element(778)==(15,UInt64[1],0,91)
end

@testset "Separate Point cell identities survive actual refinement" begin
    for native in (false,true)
      try
        setup(;native)
        for order in (1,2,1)
            API.mesh.set_order(order)
            check_points()
        end
        before=(model=API.CURRENT[],cache=API.LAST_MESH[],class=API.LAST_MESH_CLASS[],
                nodes=API.mesh.get_nodes(),elements=API.mesh.get_elements(),
                history=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]))
        @test_throws r"max_cells" API.mesh.refine(;max_cells=5)
        @test API.CURRENT[]===before.model && API.LAST_MESH[]===before.cache
        @test API.LAST_MESH_CLASS[]===before.class
        @test API.mesh.get_nodes()==before.nodes && API.mesh.get_elements()==before.elements
        @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.history
        API.mesh.refine()
        check_points()
        triangles,connections=API.mesh.get_elements_by_type(2,92)
        @test length(triangles)==4 && length(connections)==12
        @test sum(API.mesh.get_element_qualities(triangles,"volume"))≈.5 atol=2e-11
        @test allunique(API.mesh.get_nodes()[1])
        @test allunique(vcat(API.mesh.get_elements()[2]...))
      finally
        API.finalize()
      end
    end
end
end
