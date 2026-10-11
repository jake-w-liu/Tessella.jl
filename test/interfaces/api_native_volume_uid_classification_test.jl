# Additional UID admission tests use real source products and production meshes.
# The six SDK lifecycle controls are in api_native_generic_volume_source_test.jl.
module ActualVolumeUIDContractChecks
using Tessella,Test
const A=Tessella.API
function tetrahedron!()
    for (tag,x,y,z) in ((10,0.,0.,0.),(2,1.,0.,0.),(7,0.,1.,0.),(5,0.,0.,1.))
        A.model.add_point(x,y,z;tag)
    end
    for (tag,p,q) in ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5));A.model.add_line(p,q;tag);end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),(23,[3,4,-12]),(24,[11,6,-4]))
        A.model.add_curve_loop(curves;tag);A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop([21,-22,23,-24];tag=30);A.model.add_volume([30];tag=40)
end
function assembled(model,input)
    assembly=A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(0),UInt64(0),typemax(Int32),typemax(Int32),"actual UID contract")
    parameters=Tessella.Model._model_mesh_runtime_parameters(model,input.mesh)
    data=input.mesh.entity_data
    for (row,uid) in enumerate(data.external_node_tags)
        A._dim01_add_node!(assembly,uid,Tuple(input.mesh.coords[:,row]),data.node_entities[row],parameters[row])
    end
    return assembly
end
function contracts()
    A.initialize()
    try
        A.option("Mesh.Renumber",0);tetrahedron!()
        A.mesh.add_nodes(0,10,[101,102],[.125,0.,0.,.2,0.,0.])
        A.mesh.add_elements_by_type(10,15,[401,402],[101,102]);A.mesh.generate(2)
        model=A._with_model() do current;deepcopy(current);end
        input=model.meshing.mesh_authority;assembly=assembled(model,input)
        production=Tessella.Model.mesh_model_volume(model,40;_mesh_input=input,_return_volume_production=true,
            _mesh_node_high_water=A.mesh.get_max_node_tag(),_mesh_element_high_water=A.mesh.get_max_element_tag())
        mesh=production.mesh;tags=production.boundary_node_tags
        classify(m,authority,geometry,ids)=A._higher_volume_input_classification(assembly,m,geometry,authority,40,ids,"actual UID contract")
        class=classify(model,input,mesh,tags)
        @test 101 in tags && !(102 in tags)
        @test A.mesh.get_nodes(0,10)[1]==UInt64[101,102]
        @test class.node_entities[only(findall(==(UInt64(101)),tags))]==(0,Int32(10))
        # Exact ownership follows each actual UID and raw parameters follow its
        # real assembly record, including retained parameters after owner edits.
        for (row,uid) in enumerate(tags)
            position=get(assembly.node_positions,uid,Int32(0));position==0 && continue
            @test class.node_entities[row]==assembly.nodes[position].owner
            @test isequal(class.node_parametric[row],assembly.nodes[position].parameter)
        end
        bad=copy(tags);bad[only(findall(==(UInt64(101)),tags))]=typemax(UInt64)-UInt64(100)
        @test_throws ArgumentError classify(model,input,mesh,bad)
        duplicated=copy(tags);duplicated[end]=duplicated[1]
        @test_throws ArgumentError classify(model,input,mesh,duplicated)
        foreign=deepcopy(input);source=only(findall(==(UInt64(101)),foreign.mesh.entity_data.external_node_tags))
        foreign.mesh.entity_data.node_entities[source]=(3,Int32(41))
        @test_throws ArgumentError classify(model,foreign,mesh,tags)
        changed=copy(mesh.coords);changed[1,only(findall(==(UInt64(101)),tags))]=nextfloat(.125)
        altered=Tessella.MeshTypes.Mesh(changed;tets=mesh.tets)
        @test_throws ArgumentError classify(model,input,altered,tags)
        reversed=Tessella.MeshTypes.Mesh(mesh.coords;tets=mesh.tets[[2,1,3,4],:])
        old_connectivity=copy(reversed.tets)
        negative=classify(model,input,reversed,tags)
        @test reversed.tets==old_connectivity
        @test negative.node_entities==class.node_entities
        @test negative.node_parametric==class.node_parametric
        @test negative.edge_entities==class.edge_entities && negative.face_entities==class.face_entities
        # Two different Surface owners claiming one actual facet must retain
        # the existing same-dimension ambiguity diagnostic.
        ambiguous=deepcopy(input)
        first_product=ambiguous.products[(2,Int32(21))];second_product=ambiguous.products[(2,Int32(22))]
        rows=sort!(unique(vcat(first_product.rows,second_product.rows)))
        inverse=Dict(row=>Int32(index) for (index,row) in enumerate(rows))
        remapped=map(product->Tessella.Elements.ElementBlock[Tessella.Elements.ElementBlock(block.msh,
            map(slot->inverse[product.rows[slot]],block.nodes),block.tags) for block in product.blocks],(second_product,first_product))
        labels=vcat(deepcopy(second_product.element_tags),deepcopy(first_product.element_tags))
        ambiguous.products[(2,Int32(22))]=Tessella.Model._ModelMeshProduct((2,Int32(22)),rows,vcat(remapped...),labels)
        ambiguous.producer_rows[(2,Int32(22))]=rows
        @test_throws ArgumentError classify(model,ambiguous,mesh,tags)
    finally
        A.finalize()
    end
end
function interior_births()
    A.initialize()
    try
        A.option("Mesh.Renumber",0);tetrahedron!()
        A.mesh.set_size([(0,tag) for tag in (2,5,7,10)],.3);A.mesh.generate(2)
        model=A._with_model() do current;deepcopy(current);end
        input=model.meshing.mesh_authority;assembly=assembled(model,input)
        production=Tessella.Model.mesh_model_volume(model,40;_mesh_input=input,_return_volume_production=true,
            _mesh_node_high_water=A.mesh.get_max_node_tag(),_mesh_element_high_water=A.mesh.get_max_element_tag())
        tags=production.boundary_node_tags
        class=A._higher_volume_input_classification(assembly,model,production.mesh,input,40,tags,"actual interior UID contract")
        fresh=[row for (row,uid) in enumerate(tags) if !haskey(assembly.node_positions,uid)]
        @test !isempty(fresh)
        @test all(row->class.node_entities[row]==(3,Int32(40)),fresh)
        @test all(row->class.node_parametric[row]===nothing,fresh)
        @test allunique(tags) && maximum(tags)<=production.node_high_water
    finally
        A.finalize()
    end
end
function transfinite_existing_only()
    A.initialize()
    try
        A.model.add_box(0,0,0,1,1,1;tag=1)
        for (_,tag) in A.model.get_entities(1);A.mesh.set_transfinite_curve(tag,3);end
        for (_,tag) in A.model.get_entities(2);A.mesh.set_transfinite_surface(tag);end
        A.mesh.set_transfinite_volume(1);A.mesh.generate(2)
        old=A.mesh.get_nodes()[1]
        A.mesh.generate(3)
        @test Set(old)⊆Set(A.mesh.get_nodes()[1])
        @test length(A.mesh.get_elements_by_type(4,1)[1])==48
        model=A._with_model() do current;current;end
        @test model.meshing.mesh_authority!==nothing
        # This path keeps its native existing-only mask; the new ordinary UID
        # classification branch explicitly excludes that mask in source.
        @test length(model.meshing.mesh_authority.mesh.entity_data.node_entities)==27
    finally
        A.finalize()
    end
end
@testset "Actual native Volume UID contract negative and reverse checks" begin
    contracts();interior_births();transfinite_existing_only()
end
end
