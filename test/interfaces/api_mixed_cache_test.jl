using Test
using Tessella

function _api_mixed_square_fixture()
    api=Tessella.API
    api.initialize()
    for (x,y) in ((0,0),(1,0),(1,1),(0,1))
        api.model.add_point(x,y,0)
    end
    for (a,b) in ((1,2),(2,3),(3,4),(4,1))
        api.model.add_line(a,b)
    end
    api.model.add_curve_loop([1,2,3,4])
    api.model.add_plane_surface([1])
    for curve in 1:4
        api.mesh.set_transfinite_curve(curve,3)
    end
    api.mesh.set_transfinite_surface(1)
    api.mesh.set_recombine(2,1)
    return api.mesh.generate(2)
end

@testset "Mixed and record reflections use the exact determinant sign" begin
    api=Tessella.API
    api.initialize()
    try
        coords=1e-150 .* Float64[0 1 0 0;0 0 1 0;0 0 0 1]
        api.model.add_discrete_entity(3,1)
        api.mesh.add_nodes(3,1,[11,12,13,14],vec(coords))
        api.mesh.add_elements_by_type(1,4,[101],[11,12,13,14])
        mixed=MixedMesh(coords,[ElementBlock(4,reshape(Int32[1,2,3,4],4,1))])
        api._replace_mesh_cache_locked!(mixed,nothing)
        # det([-3 2 1;4 1 2;1 3 5]) = -22. Its scaled Float64
        # cofactor sum is Inf-Inf, while transformed coordinates stay finite.
        affine=zeros(4,4);affine[4,4]=1
        affine[1:3,1:3]=1e150 .* Float64[-3 2 1;4 1 2;1 3 5]
        result=api.mesh.affine_transform(affine)
        nodes=only(result.blocks).nodes[:,1]
        points=ntuple(i->Tuple(result.coords[:,nodes[i]]),4)
        @test all(isfinite,result.coords)
        @test Tessella.MeshTypes.tet_signed_volume(points...)≈22/6
        record=api.CURRENT[].discrete[(3,1)]
        @test all(isfinite,record.node_coords)
        @test record.element_nodes[1]!=Int32[11,12,13,14]
        indices=Dict(tag=>i for (i,tag) in enumerate(record.node_tags))
        record_points=ntuple(i->Tuple(record.node_coords[:,indices[record.element_nodes[1][i]]]),4)
        @test Tessella.MeshTypes.tet_signed_volume(record_points...)≈22/6
    finally
        api.finalize()
    end
end

@testset "API generation rejects displaced incident native vertex meshes" begin
    api=Tessella.API
    try
        before=_api_mixed_square_fixture()
        api.mesh.add_nodes(0,1,[100],[2.,0,0])
        @test_throws r"Point\[1\].*off its geometry.*Curve\[1\]" api.mesh.generate(2)
        @test api.mesh.get().coords==before.coords
        @test api.mesh.get().blocks[1].nodes==before.blocks[1].nodes
    finally
        api.finalize()
    end
end

@testset "API simplex merge preserves independent coincident entities" begin
    api=Tessella.API
    api.initialize()
    try
        for offset in (0,4)
            for (index,(x,y)) in enumerate(((0,0),(1,0),(1,1),(0,1)))
                api.model.add_point(x,y,0;tag=offset+index)
            end
            for (index,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
                curve=offset+index
                api.model.add_line(offset+a,offset+b;tag=curve)
                api.mesh.set_transfinite_curve(curve,3)
            end
            surface=offset==0 ? 1 : 2
            api.model.add_curve_loop(collect(offset.+(1:4));tag=surface)
            api.model.add_plane_surface([surface];tag=surface)
            api.mesh.set_transfinite_surface(surface)
        end
        generated=api.mesh.generate(2)
        @test generated isa Mesh
        @test size(generated.coords,2)==18
        @test size(generated.tris,2)==16
        @test validate(generated).ok
        one=api.mesh.get_nodes(2,1,true)[1]
        two=api.mesh.get_nodes(2,2,true)[1]
        @test length(one)==length(two)==9
        @test isempty(intersect(one,two))
        @test length(api.mesh.get_elements(2,1)[2][1])==8
        @test length(api.mesh.get_elements(2,2)[2][1])==8
        api.mesh.set_order(2)
        @test api.mesh.get_elements(2,1)[1]==Int32[9]
        @test length(api.mesh.get_elements(2,1)[3][1])==48
        @test length(api.mesh.get_nodes()[1])==50
        @test isempty(intersect(api.mesh.get_nodes(2,1,true)[1],
                                api.mesh.get_nodes(2,2,true)[1]))
    finally
        api.finalize()
    end
end

@testset "Unclassified simplex overlays retain dimension-specific elements" begin
    api=Tessella.API
    api.initialize()
    try
        mesh=Mesh(Float64[0 1 0 0;0 0 1 0;0 0 0 1];
            tets=reshape(Int32[1,2,3,4],4,1))
        api._replace_mesh_cache_locked!(mesh,nothing)
        api.mesh.set_order(2)
        global_elements=api.mesh.get_elements()
        volume_elements=api.mesh.get_elements(3)
        @test volume_elements==global_elements
        @test volume_elements[1]==Int32[11]
        @test length(only(volume_elements[3]))==10
    finally
        api.finalize()
    end
end

@testset "API curved CAD elevation and refinement fail before cache mutation" begin
    api=Tessella.API
    source="""
        Point(1)={1,0,0};Point(2)={2,0,0};Point(3)={0,2,0};
        Point(4)={0,1,0};Point(5)={0,0,0};
        Line(1)={1,2};Circle(2)={2,5,3};Line(3)={3,4};Circle(4)={4,5,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{1,2,3,4}=2;Transfinite Surface{1};Recombine Surface{1};
        """
    api.initialize()
    try
        mktemp() do path,io
            write(io,source);close(io)
            api.open_geo!(path)
        end
        before=api.mesh.generate(2)
        @test before isa MixedMesh
        @test_throws ArgumentError api.mesh.set_order(2)
        @test api.CURRENT[].meshing.order==1
        @test api.mesh.get().coords==before.coords
        @test api.mesh.get().blocks[1].nodes==before.blocks[1].nodes
        @test_throws ArgumentError api.mesh.refine()
        @test api.mesh.get().coords==before.coords
        @test api.mesh.get_element_types()==Int32[3]
    finally
        api.finalize()
    end
end

@testset "API native hexahedron cache generation, classification, and order" begin
    api=Tessella.API
    api.initialize()
    try
        api.model.add_box(0,0,0,1,1,1)
        for (_,curve) in api.model.get_entities(1)
            api.mesh.set_transfinite_curve(curve,3)
        end
        for (_,surface) in api.model.get_entities(2)
            api.mesh.set_transfinite_surface(surface)
            api.mesh.set_recombine(2,surface)
        end
        api.mesh.set_transfinite_volume(1)
        generated=api.mesh.generate(3)
        @test generated isa MixedMesh
        @test validate(generated).ok
        @test size(generated.coords,2)==27
        @test api.mesh.get_element_types(3,1)==Int32[5]
        @test length(api.mesh.get_elements(3,1)[2][1])==8
        @test length(api.mesh.get_nodes(3,1)[1])==1
        @test length(api.mesh.get_nodes(2,1)[1])==1
        api.mesh.set_order(2)
        @test size(api.mesh.get().coords,2)==125
        @test api.mesh.get_element_types(3,1)==Int32[12]
        @test length(api.mesh.get_elements(3,1)[3][1])==216
        @test length(api.mesh.get_nodes(3,1)[1])==27
        @test length(api.mesh.get_nodes(2,1)[1])==9
        @test length(api.mesh.get_nodes(1,1)[1])==3
        api.mesh.set_order(1)
        @test api.mesh.get().coords==generated.coords
        api.mesh.clear([(3,1)])
        @test isempty(api.mesh.get_element_types())
        @test size(api.mesh.get().coords,2)==26
    finally
        api.finalize()
    end
end

@testset "API mixed merge preserves coincident independent entity identity" begin
    api=Tessella.API
    try
        _api_mixed_square_fixture()
        for (x,y) in ((0,0),(1,0),(1,1),(0,1))
            api.model.add_point(x,y,0)
        end
        for (a,b) in ((5,6),(6,7),(7,8),(8,5))
            api.model.add_line(a,b)
        end
        api.model.add_curve_loop([5,6,7,8])
        api.model.add_plane_surface([2])
        for curve in 5:8
            api.mesh.set_transfinite_curve(curve,3)
        end
        api.mesh.set_transfinite_surface(2)
        api.mesh.set_recombine(2,2)
        generated=api.mesh.generate(2)
        @test generated isa MixedMesh
        @test size(generated.coords,2)==18
        @test validate(generated).ok
        @test length(api.mesh.get_elements(2,1)[2][1])==4
        @test length(api.mesh.get_elements(2,2)[2][1])==4
        one=api.mesh.get_nodes(2,1,true)[1]
        two=api.mesh.get_nodes(2,2,true)[1]
        @test isempty(intersect(one,two))
        api.mesh.set_order(2)
        @test size(api.mesh.get().coords,2)==50
        @test isempty(intersect(api.mesh.get_nodes(2,1,true)[1],
                                api.mesh.get_nodes(2,2,true)[1]))
    finally
        api.finalize()
    end
end

@testset "Native quadratic cache elevations preserve every linear reference family" begin
    api=Tessella.API
    elements=Tessella.Elements
    matrix=Float64[2 .25 .5; .5 3 .25; .25 .5 4]
    shift=[10.,-5.,3.]
    for msh in (1,2,3,4,5,6,7,15)
        reference=elements.lagrange_nodes(msh)
        dim=elements.msh_dimension(msh)
        linear=MixedMesh(matrix*reference.+shift,
            [ElementBlock(msh,reshape(Int32.(1:size(reference,2)),:,1),Int32[9])])
        owner=(dim,Int32(1))
        class=api._mixed_classification(linear,owner,[owner],
            fill(owner,size(reference,2)),Dict(owner=>Int32[]),Dict(msh=>Int32[1]))
        quadratic,quadratic_class=api._mixed_quadratic_cache(linear,class,"test")
        target=msh==15 ? 15 : elements.msh_type(elements.msh_family(msh),2)
        @test quadratic.blocks[1].msh==target
        @test quadratic.coords≈matrix*elements.lagrange_nodes(target).+shift
        @test validate(quadratic).ok
        @test quadratic_class.cell_entities[target]==Int32[1]
        reversed,reversed_class=api._mixed_reverse_cache(quadratic,quadratic_class)
        restored,_=api._mixed_reverse_cache(reversed,reversed_class)
        @test restored.blocks[1].nodes==quadratic.blocks[1].nodes
        roundtrip,_=api._mixed_linear_cache(quadratic,quadratic_class)
        @test roundtrip.coords==linear.coords
        @test roundtrip.blocks[1].nodes==linear.blocks[1].nodes
    end
end

@testset "Native mixed quadratic elevation prunes unreferenced vertices" begin
    api=Tessella.API
    coords=Float64[0 1 1 0 0.5;0 0 1 1 0.5;0 0 0 0 3]
    linear=MixedMesh(coords,[ElementBlock(3,reshape(Int32[1,2,3,4],4,1))])
    owner=(2,Int32(1))
    class=api._mixed_classification(linear,owner,[owner],fill(owner,5),
        Dict(owner=>Int32[]),Dict(3=>Int32[1]))
    quadratic,quadratic_class=api._mixed_quadratic_cache(linear,class,"test")
    @test size(quadratic.coords,2)==9
    @test !any(node->quadratic.coords[:,node]==coords[:,5],axes(quadratic.coords,2))
    restored,restored_class=api._mixed_linear_cache(quadratic,quadratic_class)
    @test restored.coords==coords[:,1:4]
    @test restored.blocks[1].nodes==linear.blocks[1].nodes
    @test restored_class.node_entities==class.node_entities[1:4]
end

@testset "API preserves native quadrangle caches and their lifecycle" begin
    api=Tessella.API
    try
        generated=_api_mixed_square_fixture()
        @test generated isa MixedMesh
        @test validate(generated).ok
        @test api.mesh.get_element_types(2,1)==Int32[3]
        types,tags,nodes=api.mesh.get_elements(2,1)
        @test types==Int32[3]
        @test length(tags[1])==4
        @test length(nodes[1])==16
        @test_throws ArgumentError api.mesh.optimize()
        @test api.mesh.optimize("",false,0)===nothing
        @test api.mesh.recombine()===nothing
        @test_throws ArgumentError api.mesh.split_quadrangles()
        @test_throws ArgumentError api.mesh.compute_cross_field()
        @test length(api.mesh.get_nodes(2,1)[1])==1
        @test length(api.mesh.get_nodes(2,1,true)[1])==9
        first_tag=tags[1][1]
        kind,connectivity,dim,owner=api.mesh.get_element(first_tag)
        @test (kind,dim,owner)==(3,2,1)
        @test connectivity==nodes[1][1:4]
        generated.coords.=17
        @test maximum(api.mesh.get().coords)==1
        before=api.mesh.get_elements(2,1)[3][1]
        api.mesh.reorder_elements(3,1,[3,2,1,0])
        @test reshape(api.mesh.get_elements(2,1)[3][1],4,:)==
              reverse(reshape(before,4,:);dims=2)
        api.mesh.reverse_elements([first_tag])
        reversed=api.mesh.get_element(first_tag)[2]
        @test reversed==reshape(api.mesh.get_elements(2,1)[3][1],4,:)[:,1]
        api.mesh.reverse_elements([first_tag])
        api.mesh.reverse([(2,1)])
        api.mesh.reverse([(2,1)])
        @test validate(api.mesh.get()).ok
        snapshot=api.mesh.get().coords
        api.mesh.affine_transform([1,0,0,2,0,1,0,3,0,0,1,4])
        @test api.mesh.get().coords==snapshot.+[2,3,4]
        @test api.mesh.get_element(first_tag)[3:4]==(2,1)
        api.mesh.set_node(1,[7,8,9])
        @test api.mesh.get_node(1)[1]==[7,8,9]
        api.mesh.set_node(1,api.mesh.get().coords[:,1].-[5,5,5])
        count=size(api.mesh.get().coords,2)
        api.mesh.renumber_nodes([1,count],[count,1])
        @test api.mesh.get_node(count)[1]==[2,3,4]
        api.mesh.renumber_nodes([1,count],[count,1])
        api.mesh.renumber_elements([Int(tags[1][1]),Int(tags[1][end])],
                                  [Int(tags[1][end]),Int(tags[1][1])])
        @test api.mesh.get_element_types(2,1)==Int32[3]
        api.mesh.remove_elements(2,1,[first_tag])
        @test length(api.mesh.get_elements(2,1)[2][1])==3
        @test size(api.mesh.get().coords,2)==count
        api.mesh.clear([(2,1)])
        @test isempty(api.mesh.get_element_types(2,1))
        @test isempty(api.mesh.get_element_types(1))
        api.model.add("other")
        @test_throws ArgumentError api.mesh.get()
        api.model.set_current("")
        @test api.mesh.get() isa MixedMesh
    finally
        api.finalize()
    end
end

@testset "API native quadratic quadrangles retain shared nodes and owners" begin
    api=Tessella.API
    try
        _api_mixed_square_fixture()
        linear=api.mesh.get()
        api.mesh.set_order(2)
        elevated=api.mesh.get()
        @test elevated isa MixedMesh
        @test api.mesh.get_element_types(2,1)==Int32[10]
        @test isempty(api.mesh.get_element_types(1))
        @test size(elevated.coords,2)==25
        @test length(api.mesh.get_nodes(2,1)[1])==9
        @test length(api.mesh.get_nodes(1,1)[1])==3
        @test length(api.mesh.get_elements(2,1)[3][1])==36
        @test validate(elevated).ok
        api.mesh.reverse([(2,1)])
        api.mesh.reverse([(2,1)])
        @test api.mesh.get().coords==elevated.coords
        @test api.mesh.get_elements(2,1)[3][1]==UInt64.(vec(elevated.blocks[1].nodes))
        api.mesh.set_order(1)
        @test api.mesh.get_element_types(2,1)==Int32[3]
        @test api.mesh.get().coords==linear.coords
        api.mesh.set_order(2)
        # Generation reads the global option; set_order changes the current mesh.
        api.option("Mesh.ElementOrder",2)
        @test api.mesh.generate(2) isa MixedMesh
        @test api.mesh.get_element_types(2,1)==Int32[10]
    finally
        api.finalize()
    end
end
