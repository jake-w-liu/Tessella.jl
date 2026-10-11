module ElementFamilyGroupingTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
const E=Tessella.Elements
const PRIMARY=joinpath(@__DIR__,"..","artifacts","same_family_query_native_v1.toml")
const NATIVE=TOML.parsefile(PRIMARY)
bits(values)=string.(reinterpret(UInt64,Float64.(values));base=16,pad=16)
function check_native_queries(expected;msh=1)
    for entity in (-1,7)
        nt,xyz,uv=A.mesh.get_nodes_by_element_type(msh,entity)
        j,d,p=A.mesh.get_jacobians(msh,msh==1 ? [0.,0.,0.] : [0.2,0.3,0.],entity)
        kt,ke,kc=A.mesh.get_keys(msh,"Lagrange",entity)
        ht,he,hc=A.mesh.get_keys(msh,"H1Legendre3",entity)
        tags,nodes=A.mesh.get_elements_by_type(msh,entity)
        actual=Dict("nodes"=>string.(nt),"xyz_bits"=>bits(xyz),"uv_bits"=>bits(uv),
            "element_tags"=>string.(tags),"element_nodes"=>string.(nodes),
            "jacobian_bits"=>bits(j),"determinant_bits"=>bits(d),"mapped_bits"=>bits(p),
            "barycenter_bits"=>bits(A.mesh.get_barycenters(msh,entity,false,false)),
            "edge_nodes"=>string.(A.mesh.get_element_edge_nodes(msh,entity)),
            "orientations"=>A.mesh.get_basis_functions_orientation(msh,"Lagrange",entity),
            "key_types"=>Int.(kt),"key_entities"=>string.(ke),"key_coordinate_bits"=>bits(kc),
            "hierarchical_key_types"=>Int.(ht),"hierarchical_key_entities"=>string.(he),
            "hierarchical_key_coordinate_bits"=>bits(hc),
            "hierarchical_orientations"=>A.mesh.get_basis_functions_orientation(msh,"H1Legendre3",entity))
        msh==2 && (actual["face_nodes"]=string.(A.mesh.get_element_face_nodes(msh,3,entity)))
        for key in keys(expected[string(entity)])
            @test actual[key]==expected[string(entity)][key]
        end
    end
end

@testset "Same-owner imported blocks form one public MSH family bucket" begin
    @test NATIVE["gmsh_version"]=="4.15.2"
    @test bytes2hex(sha256(replace(read(PRIMARY,String),"\r\n"=>"\n")))==
        "ba6897f10ec2534db5174498c610573186d04d3c1f935636f021806e2aaffc72"
    data=joinpath(@__DIR__,"..","data","native_msh_import_uint64_v1")
    # Independent Gmsh4.15.2 retained_same_owner packet in native_controls_v2.toml.
    mktempdir() do directory
        source=joinpath(directory,"same_owner.geo")
        files=("io_route_native_line_orphan2_input_v1.msh","disjoint_curve7_ascii4.msh")
        write(source,join("Merge \"$(replace(joinpath(data,name),'\\'=>'/'))\";\n" for name in files))
        A.initialize()
        try
            A.open_geo!(source)
            expected=(Int32[1],[UInt64[701,1201]],[UInt64[17,29,41,43]])
            @test A.mesh.get_elements(1,7)==expected
            @test A.mesh.get_elements()==expected
            @test A.mesh.get_element_types(1,7)==expected[1]
            @test A.mesh.get_elements_by_type(1,7)==(expected[2][1],expected[3][1])
            # The independent retained_same_owner packet explicitly rebuilds
            # scalar catalogues after its MSH2 then MSH4 imports.
            A.mesh.rebuild_node_cache(false)
            A.mesh.rebuild_element_cache(false)
            @test A.mesh.get_element(701)==(1,UInt64[17,29],1,7)
            @test A.mesh.get_element(1201)==(1,UInt64[41,43],1,7)
            check_native_queries(NATIVE["before"])
            cache=A.LAST_MESH[];snapshot=deepcopy(cache)
            @test_throws ArgumentError A.mesh.reorder_elements(1,7,[0,0])
            @test A.LAST_MESH[]===cache
            @test cache.entity_data.external_element_tags==snapshot.entity_data.external_element_tags
            A.mesh.reorder_elements(1,7,[1,0])
            @test A.mesh.get_elements(1,7)==(Int32[1],[UInt64[1201,701]],[UInt64[41,43,17,29]])
            check_native_queries(NATIVE["reordered"])
            @test A.mesh.get_element(701)==(1,UInt64[17,29],1,7)
            @test A.mesh.get_element(1201)==(1,UInt64[41,43],1,7)
        finally
            A.finalize()
        end
    end
end

@testset "Parent-family node queries retain interleaved polynomial order" begin
    A.initialize()
    try
        mktempdir() do directory
            data=joinpath(@__DIR__,"..","data","native_msh_import_uint64_v1")
            names=("io_route_native_line_orphan2_input_v1.msh",
                "same_owner_quadratic_curve7_ascii4.msh","disjoint_curve7_ascii4.msh")
            path=joinpath(directory,"mixed_order.geo")
            write(path,join("Merge \"$(replace(joinpath(data,name),'\\'=>'/'))\";\n" for name in names))
            A.open_geo!(path)
            # Isolated Gmsh4.15.2 returns the same actual parent-family incidences
            # for requested MSH1/MSH8 and the aggregate/entity7 selectors.
            expected_nodes=UInt64[17,29,53,61,71,41,43]
            expected_xyz=Float64[0,-0.,0,1,0,0,4,0,0,5,0,0,4.5,0,0,2,-0.,0,3,0,0]
            cache=A.LAST_MESH[];before=deepcopy(cache)
            for msh in (1,8),entity in (-1,7),parametric in (false,true)
                nodes,xyz,uv=A.mesh.get_nodes_by_element_type(msh,entity,parametric)
                @test nodes==expected_nodes
                @test isequal(xyz,expected_xyz)
                @test isempty(uv)
            end
            @test A.mesh.get_elements_by_type(1,7)==(UInt64[701,1201],UInt64[17,29,41,43])
            @test A.mesh.get_elements_by_type(8,7)==(UInt64[901],UInt64[53,61,71])
            @test A.LAST_MESH[]===cache
            @test cache.entity_data.external_element_tags==before.entity_data.external_element_tags
            @test [block.nodes for block in cache.blocks]==[block.nodes for block in before.blocks]
        end
    finally A.finalize() end
end

@testset "Repeated Triangle family keeps face incidences and public bubble identities" begin
    primary=joinpath(@__DIR__,"..","artifacts","triangle_family_query_native_v1.toml")
    native=TOML.parsefile(primary)
    @test bytes2hex(sha256(replace(read(primary,String),"\r\n"=>"\n")))==
        "5e50932e92392ea9283f8a177c94eb9aab559bedd6eeb83527f8cca2a4e4f6e3"
    A.initialize()
    try
        mktempdir() do directory
            input=joinpath(@__DIR__,"..","data","native_msh_import_uint64_v1","duplicate_triangle_family_ascii4.msh")
            path=joinpath(directory,"triangles.geo");write(path,"Merge \"$(replace(input,'\\'=>'/'))\";\n")
            A.open_geo!(path)
            check_native_queries(native["before"];msh=2)
            # A valid detached cache may store the two Triangle blocks around
            # an unrelated Line block, while representing the same native mesh.
            base=A.LAST_MESH[];data=base.entity_data
            triangle=only(block for block in base.blocks if block.msh==2)
            line=only(block for block in base.blocks if block.msh==1)
            blocks=[E.ElementBlock(2,triangle.nodes[:,1:1],triangle.tags[1:1]),line,
                E.ElementBlock(2,triangle.nodes[:,2:2],triangle.tags[2:2])]
            metadata=E._native_mixed_entity_data(data.entities;
                node_entities=data.node_entities,node_parametric=data.node_parametric,
                external_node_tags=data.external_node_tags,
                block_entities=[Int32[7],Int32[8],Int32[7]],
                external_element_tags=[UInt64[701],UInt64[901],UInt64[1201]])
            interleaved=E.MixedMesh(base.coords,blocks;entity_data=metadata)
            A._replace_mesh_cache_locked!(interleaved,A._higher_classification_from_data(interleaved,"interleaved Triangle fixture"))
            A._native_lookup_invalidate_locked!()
            @test [block.msh for block in A.LAST_MESH[].blocks]==[2,1,2]
            check_native_queries(native["before"];msh=2)
            A.mesh.reorder_elements(2,7,[1,0])
            check_native_queries(native["reordered"];msh=2)
            @test A.mesh.get_element(701)==(2,UInt64[17,29,41],2,7)
            @test A.mesh.get_element(1201)==(2,UInt64[43,53,61],2,7)
            @test A.mesh.get_element(901)==(1,UInt64[43,53],1,8)
            @test A.mesh.get_elements_by_type(1,8)==(UInt64[901],UInt64[43,53])
            @test [block.msh for block in A.LAST_MESH[].blocks]==[2,1,2]
        end
    finally A.finalize() end
end

@testset "Family grouping preserves sparse labels, dense slots and task slices" begin
    A.initialize()
    try
        A.model.add_discrete_entity(1,7)
        coords=Float64[0 1 2 3;-0. 0. -0. 0.;0 0 0 0]
        blocks=[E.ElementBlock(1,reshape(Int32[1,2],2,1)),
            E.ElementBlock(1,reshape(Int32[3,4],2,1))]
        metadata=E.MixedEntityData(Dict((1,7)=>E.MixedEntity(1,7,[0.,0.,0.,3.,0.,0.]));
            node_entities=fill((1,Int32(7)),4),node_parametric=fill(nothing,4),
            external_node_tags=UInt64[17,29,41,43],
            block_entities=[Int32[7],Int32[7]],
            external_element_tags=[UInt64[701],UInt64[1201]])
        cache=E.MixedMesh(coords,blocks;entity_data=metadata)
        class=A._higher_classification_from_data(cache,"family grouping fixture")
        A._replace_mesh_cache_locked!(cache,class)
        before=(copy(cache.coords),[copy(block.nodes) for block in cache.blocks],
            deepcopy(cache.entity_data.external_element_tags))
        expected=(Int32[1],[UInt64[701,1201]],[UInt64[17,29,41,43]])
        @test A.mesh.get_elements(1,7)==expected
        @test A.mesh.get_elements(1)==expected
        @test A.mesh.get_elements_by_type(1,7)==(expected[2][1],expected[3][1])
        @test A.mesh.get_elements_by_type(1)==(expected[2][1],expected[3][1])
        for tasks in (2,3,4),task in 0:tasks-1
            chosen=A._mesh_task_range(2,task,tasks)
            @test A.mesh.get_elements_by_type(1,7,task,tasks)==
                (expected[2][1][chosen],vec(reshape(expected[3][1],2,2)[:,chosen]))
        end
        @test A.mesh.get_element(701)==(1,UInt64[17,29],1,7)
        @test A.mesh.get_element(1201)==(1,UInt64[41,43],1,7)
        @test A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===class
        @test isequal(cache.coords,before[1])
        @test [block.nodes for block in cache.blocks]==before[2]
        @test cache.entity_data.external_element_tags==before[3]
        result=A.mesh.get_elements(1,7)
        result[2][1][1]=999;result[3][1][1]=999
        @test A.mesh.get_elements(1,7)==expected

        # Another family occupies dense slot2 between the two matching blocks;
        # its public label is independent of the family-relative positions.
        A.model.add_discrete_entity(2,9)
        interleaved=[blocks[1],E.ElementBlock(2,reshape(Int32[1,2,3],3,1)),blocks[2]]
        metadata2=E.MixedEntityData(merge(metadata.entities,
            Dict((2,9)=>E.MixedEntity(2,9,[0.,0.,0.,3.,0.,0.])));
            node_entities=metadata.node_entities,node_parametric=metadata.node_parametric,
            external_node_tags=metadata.external_node_tags,
            block_entities=[Int32[7],Int32[9],Int32[7]],
            external_element_tags=[UInt64[701],UInt64[901],UInt64[1201]])
        cache2=E.MixedMesh(coords,interleaved;entity_data=metadata2)
        A._replace_mesh_cache_locked!(cache2,A._higher_classification_from_data(cache2,"interleaved fixture"))
        A._native_lookup_invalidate_locked!()
        @test A.mesh.get_elements_by_type(1,7)==(expected[2][1],expected[3][1])
        @test A.mesh.get_elements_by_type(1)==(expected[2][1],expected[3][1])
        @test A.mesh.get_elements_by_type(2,9)==(UInt64[901],UInt64[17,29,41])
        @test A.mesh.get_elements_by_type(1,7,1,2)==(UInt64[1201],UInt64[41,43])
        @test A.mesh.get_element(901)==(2,UInt64[17,29,41],2,9)
        A.model.add_discrete_entity(1,8)
        @test A.mesh.get_elements_by_type(1,8)==(UInt64[],UInt64[])
    finally
        A.finalize()
    end
end

@testset "Repeated family query allocation growth stays linear" begin
    A.initialize()
    allocations=Int[]
    try
        A.model.add_discrete_entity(1,7)
        for count in (1000,2000,4000)
            block=E.ElementBlock(1,repeat(reshape(Int32[1,2],2,1),1,count))
            labels=[collect(UInt64,1:count),collect(UInt64,count+1:2count)]
            data=E._native_mixed_entity_data(Dict((1,7)=>E.MixedEntity(1,7,[0.,0.,0.,1.,0.,0.]));
                node_entities=fill((1,Int32(7)),2),node_parametric=fill(nothing,2),
                external_node_tags=UInt64[17,29],block_entities=[fill(Int32(7),count),fill(Int32(7),count)],
                external_element_tags=labels)
            cache=E.MixedMesh(Float64[0 1;0 0;0 0],[block,block];entity_data=data)
            A._replace_mesh_cache_locked!(cache,A._higher_classification_from_data(cache,"repeated family resource fixture"))
            getter=()->A.mesh.get_elements_by_type(1,7)
            getter()
            push!(allocations,minimum(@allocated(getter()) for _ in 1:3))
            tags,nodes=getter()
            @test tags==vcat(labels...)
            @test nodes==repeat(UInt64[17,29],2count)
            @test A.mesh.get_elements_by_type(1,7,2count-1,2count)==(UInt64[2count],UInt64[17,29])
            @test A.LAST_MESH[]===cache
        end
        @test allocations[2]<=2.15*allocations[1]+65536
        @test allocations[3]<=2.15*allocations[2]+65536
        println("Repeated family public query allocation sweep: ",allocations)
    finally A.finalize() end
end
end
