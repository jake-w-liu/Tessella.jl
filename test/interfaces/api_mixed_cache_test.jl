using Test
using Tessella

function _api_mixed_square_fixture(;renumber=1)
    api=Tessella.API
    api.initialize()
    api.option("Mesh.Renumber",renumber)
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

include("../artifacts/displaced_point_generation_primary_v1.jl")
@testset "API generation retains extra Point nodes by actual entity identity" begin
    api=Tessella.API
    for row in DISPLACED_POINT_GENERATION_PRIMARY
        try
            _api_mixed_square_fixture(;renumber=row["renumber"])
            api.mesh.add_nodes(0,1,[100],[2.,0,0])
            @test row["error"]===nothing
            api.mesh.generate(2)
            for entity in row["after_generation"]
                dimension,owner=entity["entity"]
                tags,coordinates,parameters=api.mesh.get_nodes(dimension,owner,false,true)
                @test tags==entity["tags"]
                @test coordinates==entity["coordinates"]
                @test parameters==entity["parameters"]
                types,elements,nodes=api.mesh.get_elements(dimension,owner)
                @test types==entity["types"]
                @test elements==entity["elements"]
                @test nodes==entity["nodes"]
            end
        finally
            api.finalize()
        end
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

@testset "API curved CAD elevation and refinement update the cache" begin
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
        mktempdir() do directory
            path=joinpath(directory,"curved.geo")
            write(path,source)
            api.open_geo!(path)
        end
        before=api.mesh.generate(2)
        @test before isa MixedMesh
        api.mesh.set_order(2)
        @test api.CURRENT[].meshing.order==2
        elevated=api.mesh.get()
        @test api.mesh.get_element_types()==Int32[8,10,15]
        @test size(elevated.coords,2)==9
        # Circle-arc midside nodes sit on the arcs, not on the chords.
        @test elevated.coords[:,6]≈[sqrt(2),sqrt(2),0.0]
        @test elevated.coords[:,8]≈[sqrt(2)/2,sqrt(2)/2,0.0]
        api.mesh.refine()
        @test api.mesh.get().coords!=elevated.coords
        @test size(api.mesh.get().coords,2)==10
        @test api.mesh.get_element_types()==Int32[1,3,15]
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
        # `clear` is per-entity upstream — the volume record is removed but
        # the boundary curve, surface and point records survive.
        @test api.mesh.get_element_types()==Int32[1,3,15]
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
        # Re-derived 2026-10-26 against Gmsh 4.15.2: `optimize` runs the
        # native tagged pipeline which, like upstream OptimizeMesh, only
        # optimizes Region cells — a 2-D quadrangle cache passes through
        # unchanged (tags [13..16] kept) instead of throwing.
        @test api.mesh.optimize()===nothing
        @test api.mesh.optimize("",false,0)===nothing
        @test api.mesh.recombine()===nothing
        # Cross fields still require a linear Tri3 surface cache; upstream's
        # quad attempt is an unreproducible PETSc crash, not a clean reject.
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
        pre_reverse=api.mesh.get_element(first_tag)[2]
        api.mesh.reverse_elements([first_tag])
        reversed=api.mesh.get_element(first_tag)[2]
        # Gmsh 4.15.2 reverseElements emits [a,d,c,b] for a quadrangle, and the
        # reordered block stores the element at its new column, not column 1.
        @test reversed==pre_reverse[[1,4,3,2]]
        block_tags,block_nodes=api.mesh.get_elements(2,1)[2:3]
        @test reversed==
              reshape(block_nodes[1],4,:)[:,findfirst(==(first_tag),block_tags[1])]
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
        # Per-entity upstream: clearing Surface 1 drops only its quadrangles;
        # the boundary curve Line2 records survive (dim-1 reports [1]).
        @test api.mesh.get_element_types(1)==Int32[1]
        api.model.add("other")
        @test_throws ArgumentError api.mesh.get()
        api.model.set_current("")
        @test api.mesh.get() isa MixedMesh
        # Re-derived 2026-10-26 against Gmsh 4.15.2 splitQuadrangles: the four
        # quadrangles become eight Tri3 under fresh tags (upstream emits
        # 19..22,25,26,31,32 here) and the retired quad tags report
        # "Unknown element" on lookup. The Tri3 cache then satisfies the
        # cross-field Tri3 admission gate.
        api.mesh.generate(2)
        @test api.mesh.get_element_types(2,1)==Int32[3]
        quad_tags=api.mesh.get_elements(2,1)[2][1]
        @test api.mesh.split_quadrangles()===nothing
        @test api.mesh.get_element_types(2,1)==Int32[2]
        @test length(api.mesh.get_elements(2,1)[2][1])==8
        @test_throws ErrorException api.mesh.get_element(quad_tags[1])
        @test length(api.mesh.compute_cross_field())==3
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
        # Re-derived 2026-10-26 against Gmsh 4.15.2: order-2 elevation keeps
        # the boundary curves as Line3 elements — upstream reports [8] too.
        @test api.mesh.get_element_types(1)==Int32[8]
        @test size(elevated.coords,2)==25
        @test length(api.mesh.get_nodes(2,1)[1])==9
        @test length(api.mesh.get_nodes(1,1)[1])==3
        @test length(api.mesh.get_elements(2,1)[3][1])==36
        @test validate(elevated).ok
        api.mesh.reverse([(2,1)])
        api.mesh.reverse([(2,1)])
        @test api.mesh.get().coords==elevated.coords
        # Re-derived 2026-10-26 against Gmsh 4.15.2: `get_elements` reports
        # public node tags (entity-birth order, byte-identical to upstream's
        # Quad9 block) while `blocks[].nodes` are internal indices ordered
        # differently. Verify the published tag and stored index bind the same
        # node — coordinates agree within the ~1e-13 transfinite jitter.
        node_tags,node_coords,_=api.mesh.get_nodes()
        public_point=Dict(node_tags[k]=>
            Tuple(node_coords[3k-2:3k]) for k in eachindex(node_tags))
        published=api.mesh.get_elements(2,1)[3][1]
        @test all(eachindex(published)) do j
            published_point=public_point[published[j]]
            stored_point=Tuple(elevated.coords[:,elevated.blocks[1].nodes[j]])
            all(k->abs(published_point[k]-stored_point[k])<=1e-9,1:3)
        end
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
