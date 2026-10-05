using Test, Tessella
using Tessella.Elements: ElementBlock, MixedMesh

const _PEL=Tessella.API

# The independent Gmsh 4.15.2 fixture is retained in test/tmp/
# b4_cross_family_renumber_primary_v1.json, SHA256
# D46F930B4EF99004CCCE33D473FBCBA2469C5B4C3DA836AA64C43AFEACB87D29.
# Its separated exact unit Tet4/Hex8/Pyr5 maps establish the simultaneous
# labels below; stored geometry is independently checked against these maps.
const _PEL_COORDS=(Float64[0 1 0 0;0 0 1 0;0 0 0 1],
    Float64[3 4 4 3 3 4 4 3;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1],
    Float64[6 7 7 6 6.5;0 0 1 1 .5;0 0 0 0 1])
const _PEL_NODES=([101,105,109,113],[211,215,219,223,227,231,235,239],
                  [311,315,319,323,327])

function _pel_mixed(;tagged=false,quadratic=false,same_owner=false)
    _PEL.initialize();_PEL.option("Mesh.Renumber",0)
    owners=same_owner ? [91,91,91] : [11,22,33]
    for owner in unique(owners);_PEL.model.add_discrete_entity(3,owner);end
    if tagged
        for (msh,owner,nodes,coords,label) in zip((4,5,7),owners,_PEL_NODES,_PEL_COORDS,(1001,7007,9009))
            _PEL.mesh.add_nodes(3,owner,nodes,vec(coords))
            _PEL.mesh.add_elements_by_type(owner,msh,[label],nodes)
        end
        _PEL.mesh.generate(0)
    else
        coordinates=hcat(_PEL_COORDS...)
        blocks=ElementBlock[];node_owners=Tuple{Int,Int32}[];offset=0
        for (msh,owner,coords) in zip((4,5,7),owners,_PEL_COORDS)
            width=size(coords,2)
            push!(blocks,ElementBlock(msh,reshape(Int32.(offset.+(1:width)),width,1)))
            append!(node_owners,fill((3,Int32(owner)),width));offset+=width
        end
        mesh=MixedMesh(coordinates,blocks)
        entities=[(3,Int32(owner)) for owner in unique(owners)]
        class=_PEL._mixed_classification(mesh,first(entities),entities,node_owners,
            Dict{Tuple{Int,Int32},Vector{Int32}}(),
            Dict(msh=>Int32[owner] for (msh,owner) in zip((4,5,7),owners)))
        _PEL._replace_mesh_cache_locked!(mesh,class)
        _PEL.mesh.renumber_elements([1,2,3],[1001,7007,9009])
    end
    quadratic && _PEL.mesh.set_order(2)
    return owners
end

function _pel_native_square(;quadratic=false)
    _PEL.initialize();_PEL.option("Mesh.Renumber",0)
    for (x,y) in ((0,0),(1,0),(1,1),(0,1));_PEL.model.add_point(x,y,0);end
    for (a,b) in ((1,2),(2,3),(3,4),(4,1));_PEL.model.add_line(a,b);end
    _PEL.model.add_curve_loop([1,2,3,4]);_PEL.model.add_plane_surface([1])
    _PEL.model.add_physical_group(2,[1];tag=99)
    for curve in 1:4;_PEL.mesh.set_transfinite_curve(curve,3);end
    _PEL.mesh.set_transfinite_surface(1);_PEL.mesh.generate(2)
    quadratic && _PEL.mesh.set_order(2)
    @test _PEL.LAST_MESH[] isa Mesh
end

function _pel_state()
    slot=_PEL.MODEL_SLOTS[_PEL._current_slot_index_locked()]
    return (model=_PEL.CURRENT[],value=repr(_PEL.CURRENT[]),mesh=_PEL.LAST_MESH[],
        class=_PEL.LAST_MESH_CLASS[],overlay=_PEL.LAST_MESH_HIGH_ORDER[],
        parent=_PEL.LAST_MESH_HIGH_ORDER_MESH[],mids=_PEL.LAST_MESH_HIGH_ORDER_MIDS[],
        nodes=_PEL.mesh.get_nodes(),elements=_PEL.mesh.get_elements(),
        maxima=(_PEL.NODE_TAG_MAX[],_PEL.ELEMENT_TAG_MAX[]),
        flags=copy(_PEL.ELEMENT_VISIBILITY[]),windows=deepcopy(slot.window_element_visibility))
end

function _pel_unchanged(before)
    after=_pel_state()
    for key in (:model,:mesh,:class,:overlay,:parent,:mids);@test getproperty(after,key)===getproperty(before,key);end
    for key in (:value,:nodes,:elements,:maxima,:flags,:windows);@test getproperty(after,key)==getproperty(before,key);end
end

@testset "Public element labels preserve native and tagged mixed slots" begin
    cases=((Int[],Int[],[1,2,3]),([1001],[55],[55,56,57]),
        ([666666],[44],[45,46,47]),([1001,1001],[17,18],[18,19,20]),
        ([1001,7007,9009],[7007,9009,1001],[7007,9009,1001]),
        ([1001,7007,9009],[1,2,3],[1,2,3]),
        ([1001,7007,9009],[12001,15001,17001],[12001,15001,17001]))
    for tagged in (false,true),same_owner in (false,true),(old,new,expected) in cases
        try
            owners=_pel_mixed(;tagged,same_owner)
            cache=_PEL.LAST_MESH[];coords=copy(cache.coords)
            connectivity=[copy(block.nodes) for block in cache.blocks]
            node_payload=_PEL.mesh.get_nodes();class=_PEL.LAST_MESH_CLASS[]
            _PEL.mesh.set_visibility([1001,9009],0)
            _PEL.mesh.renumber_elements(old,new)
            @test _PEL.mesh.get_elements()[1]==Int32[4,5,7]
            @test only.(_PEL.mesh.get_elements()[2])==UInt64.(expected)
            @test _PEL.LAST_MESH[].coords==coords
            @test [block.nodes for block in _PEL.LAST_MESH[].blocks]==connectivity
            @test _PEL.mesh.get_nodes()==node_payload
            @test _PEL.mesh.get_visibility(expected)==Int32[0,1,0]
            @test _PEL.LAST_MESH_CLASS[].edge_entities==class.edge_entities
            @test _PEL.mesh.get_max_element_tag()==max(9009,maximum(expected))
            for (msh,owner,label,nodes) in zip((4,5,7),owners,expected,connectivity)
                @test _PEL.mesh.get_element(label)[1]==msh
                @test _PEL.mesh.get_element(label)[3:4]==(3,owner)
                @test _PEL.mesh.get_elements_by_type(msh)[1]==UInt64[label]
                @test only(_PEL.mesh.get_elements(3,owner)[2][findfirst(==(msh),_PEL.mesh.get_elements(3,owner)[1])])==label
                @test _PEL.mesh.get_element_qualities([label],"minEdge")[1]≈1 atol=2e-11
            end
            next=_PEL.mesh.get_max_element_tag()+1
            _PEL.mesh.add_nodes(3,owners[1],[],vec(_PEL_COORDS[1].+[10,0,0]))
            added=_PEL.mesh.get_nodes(3,owners[1])[1][end-3:end]
            _PEL.mesh.add_elements_by_type(owners[1],4,[],added)
            @test _PEL.mesh.get_element(next)[1]==4
        finally
            _PEL.finalize()
        end
    end
end

@testset "Native simplex P1/P2 labels keep actual maps and query identities" begin
    for quadratic in (false,true),edited in (false,true)
        try
            _pel_native_square(;quadratic)
            if quadratic && edited
                overlay=_PEL.LAST_MESH_HIGH_ORDER[]
                # A finite singular map is still stored/queryable. Label edits
                # must not newly reject or rebuild its interpolation rows.
                overlay.coords[:,overlay.tri6[4,1]]=overlay.coords[:,overlay.tri6[1,1]]
            end
            before=_pel_state();types,tags,nodes=before.elements
            old=vcat(tags...);new=reverse(old).+1000
            records=Dict(tag=>_PEL.mesh.get_element(tag) for tag in old)
            jacobians=Dict(tag=>_PEL.mesh.get_jacobian(tag,[.2,.2,0]) for tag in old)
            located=edited ? UInt64[] : _PEL.mesh.get_elements_by_coordinates(.5,.5,0,2,true)
            quality=_PEL.mesh.get_element_qualities(old,"volume")
            _PEL.mesh.set_visibility(old[1:2:end],0)
            _PEL.mesh.renumber_elements(old,new)
            @test _PEL.LAST_MESH[]===before.mesh
            @test _PEL.LAST_MESH_HIGH_ORDER[]===before.overlay
            @test _PEL.LAST_MESH_HIGH_ORDER_MESH[]===before.parent
            @test _PEL.LAST_MESH_HIGH_ORDER_MIDS[]===before.mids
            @test _PEL.mesh.get_nodes()==before.nodes
            @test _PEL.mesh.get_elements()[1]==types
            @test _PEL.mesh.get_elements()[3]==nodes
            @test _PEL.mesh.get_element_qualities(new,"volume")==quality
            for (source,target) in zip(old,new)
                @test _PEL.mesh.get_element(target)==records[source]
                @test _PEL.mesh.get_jacobian(target,[.2,.2,0])==jacobians[source]
            end
            if !edited
                mapping=Dict(zip(old,new))
                expected=sort!([mapping[tag] for tag in located])
                @test _PEL.mesh.get_elements_by_coordinates(.5,.5,0,2,true)==expected
                answer=_PEL.mesh.get_element_by_coordinates(.5,.5,0,2,true)
                source=old[findfirst(==(first(expected)),new)]
                @test answer[1]==first(expected)
                @test answer[2:3]==records[source][1:2]
            end
            table=_PEL.LAST_MESH_CLASS[].public_tags
            @test table.mesh===before.mesh
            @test length(table.node_tags)==length(before.nodes[1])
            for (badold,badnew) in ((old,[0;new[2:end]]),(new,fill(77,length(new))),
                    (new,new[1:end-1]),([new[1]],[typemax(Int32)+1]),([new[1]],[typemax(Int32)]))
                snapshot=_pel_state()
                @test_throws ArgumentError _PEL.mesh.renumber_elements(badold,badnew)
                _pel_unchanged(snapshot)
            end
            snapshot=_pel_state();_PEL.mesh.set_order(quadratic ? 2 : 1)
            @test _PEL.LAST_MESH[]===snapshot.mesh
            @test _PEL.LAST_MESH_HIGH_ORDER[]===snapshot.overlay
            @test _PEL.mesh.get_elements()==snapshot.elements
            _PEL.mesh.reverse_elements([new[1]])
            @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===_PEL.LAST_MESH[]
            @test _PEL.mesh.get_element(new[1])[1]==records[old[1]][1]
        finally
            _PEL.finalize()
        end
    end
end

@testset "Actual quadratic mixed labels retain geometry, owners and visibility" begin
    for tagged in (false,true),same_owner in (false,true)
        try
            owners=_pel_mixed(;tagged,same_owner,quadratic=true)
            before=_pel_state()
            @test before.elements[1]==Int32[11,12,14]
            old=vcat(before.elements[2]...)
            @test old==UInt64[9010,9011,9012]
            new=old[[2,3,1]]
            scalar=Dict(tag=>_PEL.mesh.get_element(tag) for tag in old)
            jac=Dict(tag=>_PEL.mesh.get_jacobian(tag,[.2,.2,.1]) for tag in old)
            _PEL.mesh.set_visibility(old[[1,3]],0)
            _PEL.mesh.renumber_elements(old,new)
            @test _PEL.mesh.get_elements()[2]==[UInt64[tag] for tag in new]
            @test _PEL.mesh.get_elements()[3]==before.elements[3]
            @test _PEL.mesh.get_nodes()==before.nodes
            @test _PEL.mesh.get_visibility(new)==Int32[0,1,0]
            @test _PEL.mesh.get_max_element_tag()==9012
            for (source,target) in zip(old,new)
                @test _PEL.mesh.get_element(target)==scalar[source]
                @test _PEL.mesh.get_jacobian(target,[.2,.2,.1])==jac[source]
                @test _PEL.mesh.get_keys_for_element(target,"Lagrange")[2]==scalar[source][2]
            end
        finally
            _PEL.finalize()
        end
    end
end

@testset "Native label tables rebind through actual P2 lifecycle" begin
    try
        _pel_native_square()
        old_nodes=_PEL.mesh.get_nodes();old_elements=vcat(_PEL.mesh.get_elements()[2]...)
        labels=old_elements.+9000
        _PEL.mesh.renumber_elements(old_elements,labels)
        renamed_nodes=reverse(old_nodes[1]).+2000
        _PEL.mesh.renumber_nodes(old_nodes[1],renamed_nodes)
        @test _PEL.mesh.get_nodes()[1]==renamed_nodes
        for (index,label) in enumerate(renamed_nodes)
            @test _PEL.mesh.get_node(label)[1]==old_nodes[2][3index-2:3index]
        end
        linear_tags,linear_nodes=_PEL.mesh.get_elements_by_type(2)
        orientations=Int32[]
        for (column,tag) in enumerate(linear_tags)
            vertices=linear_nodes[3column-2:3column]
            rank=sum(sum(vertices[j]<vertices[i] for j in i+1:3)*factorial(3-i) for i in 1:2)
            push!(orientations,Int32(rank))
            @test _PEL.mesh.get_basis_functions_orientation_for_element(tag,"H1Legendre3")==rank
        end
        @test _PEL.mesh.get_basis_functions_orientation(2,"H1Legendre3")==orientations
        previous_node=_PEL.mesh.get_max_node_tag();previous_element=_PEL.mesh.get_max_element_tag()
        _PEL.mesh.set_order(2)
        @test _PEL.LAST_MESH[] isa Mesh
        @test _PEL.LAST_MESH_HIGH_ORDER[] isa P2TriMesh
        @test _PEL.mesh.get_nodes()[1][1:length(renamed_nodes)]==renamed_nodes
        @test minimum(_PEL.mesh.get_nodes()[1][length(renamed_nodes)+1:end])==previous_node+1
        tags,nodes=_PEL.mesh.get_elements_by_type(9)
        @test _PEL.mesh.get_nodes_by_element_type(9)[1]==nodes
        physical=_PEL.mesh.get_nodes_for_physical_group(2,99)
        expected=sort!(unique(_PEL.mesh.get_nodes(2,1,true)[1]))
        @test physical[1]==expected
        for (index,label) in enumerate(physical[1])
            @test physical[2][3index-2:3index]==_PEL.mesh.get_node(label)[1]
        end
        @test tags==previous_element .+ UInt64.(1:length(tags))
        for (column,tag) in enumerate(tags)
            connectivity=nodes[6column-5:6column]
            @test _PEL.mesh.get_element(tag)[2]==connectivity
            @test _PEL.mesh.get_keys_for_element(tag,"Lagrange")[2]==connectivity
            space=Tessella.MeshFunctionSpaces._function_space("H1Legendre5","label test")
            counts=Tessella.MeshFunctionSpaces._h1_counts(:tri,space.key_order)
            @test length(_PEL.mesh.get_keys_for_element(tag,"H1Legendre5")[2])==counts.vertex+counts.edge+counts.face+counts.bubble
        end
        scalar=Dict(tag=>_PEL.mesh.get_element(tag) for tag in tags)
        _PEL.mesh.set_visibility(tags[1:2:end],0)
        visibility=_PEL.mesh.get_visibility(tags)
        _PEL.mesh.reorder_elements(9,1,reverse(collect(0:length(tags)-1)))
        @test _PEL.mesh.get_elements_by_type(9)[1]==reverse(tags)
        @test _PEL.mesh.get_visibility(tags)==visibility
        for tag in tags;@test _PEL.mesh.get_element(tag)==scalar[tag];end
        @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===_PEL.LAST_MESH[]
        _PEL.mesh.remove_elements(2,1,[tags[1]])
        @test Set(_PEL.mesh.get_elements_by_type(9)[1])==Set(tags[2:end])
        for tag in tags[2:end];@test _PEL.mesh.get_element(tag)==scalar[tag];end
        @test _PEL.mesh.get_visibility(tags[2:end])==visibility[2:end]
        coordinates=_PEL.mesh.get_nodes()[2]
        _PEL.mesh.affine_transform([1.,0,0,2,0,1,0,3,0,0,1,4])
        moved=_PEL.mesh.get_nodes()[2]
        @test reshape(moved,3,:)≈reshape(coordinates,3,:).+[2,3,4] atol=2e-11
        @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===_PEL.LAST_MESH[]
        _PEL.mesh.clear([(2,1)])
        @test isempty(_PEL.mesh.get_elements()[1])
        @test _PEL.mesh.get_max_element_tag()>=maximum(tags)
    finally
        _PEL.finalize()
    end
    try
        _pel_native_square(;quadratic=true)
        old=vcat(_PEL.mesh.get_elements()[2]...)
        _PEL.mesh.renumber_elements(old,old.+9000)
        highwater=_PEL.mesh.get_max_element_tag()
        before=_PEL.mesh.get_nodes()
        result=_PEL.mesh.refine()
        @test result isa Mesh && _PEL.LAST_MESH[] isa Mesh
        @test _PEL.LAST_MESH_HIGH_ORDER[] isa P2TriMesh
        @test _PEL.mesh.get_elements()[1]==Int32[9]
        tags=vcat(_PEL.mesh.get_elements()[2]...)
        @test minimum(tags)==highwater+1
        @test sum(_PEL.mesh.get_element_qualities(tags,"volume"))≈1 atol=2e-11
        nodes=_PEL.mesh.get_nodes()
        @test length(unique(nodes[1]))==length(nodes[1])
        @test Set(before[1])⊆Set(nodes[1])
        for tag in before[1]
            index=findfirst(==(tag),before[1])
            @test _PEL.mesh.get_node(tag)[1]==before[2][3index-2:3index]
        end
        node_labels=nodes[1][1:Tessella.MeshTypes.nnodes(_PEL.LAST_MESH[])]
        _PEL.mesh.set_order(1)
        @test _PEL.mesh.get_nodes()[1]==node_labels
        @test _PEL.LAST_MESH_HIGH_ORDER[]===nothing
        node_max=_PEL.mesh.get_max_node_tag()
        _PEL.mesh.set_order(2)
        @test minimum(_PEL.mesh.get_nodes()[1][length(node_labels)+1:end])==node_max+1
        @test _PEL.LAST_MESH[] isa Mesh
    finally
        _PEL.finalize()
    end
end

@testset "Actual native tetrahedral bubble keys use public element labels" begin
    try
        _PEL.initialize();_PEL.option("Mesh.Renumber",0)
        _PEL.model.add_discrete_entity(3,1)
        mesh=Mesh(_PEL_COORDS[1];tets=reshape(Int32[1,2,3,4],4,1))
        class=_PEL._MeshClassification(mesh,(3,Int32(1)),[(3,Int32(1))],
            fill((3,Int32(1)),4),Dict{Tuple{Int,Int32},Vector{Int32}}(),
            Int32[],Int32[],Int32[1])
        _PEL._replace_mesh_cache_locked!(mesh,class)
        _PEL.mesh.renumber_elements([1],[1001])
        _PEL.mesh.renumber_nodes([1,2,3,4],[9,7,5,3])
        @test _PEL.mesh.get_keys_for_element(1001,"H1Legendre4")[2][end]==1001
        @test _PEL.mesh.get_keys(4,"H1Legendre4")[2][end]==1001
        @test _PEL.mesh.get_basis_functions_orientation_for_element(1001,"H1Legendre4")==23
        _PEL.mesh.set_order(2)
        @test _PEL.mesh.get_element(1002)[1]==11
        @test _PEL.mesh.get_keys_for_element(1002,"H1Legendre4")[2][end]==1002
        @test _PEL.mesh.get_keys(11,"H1Legendre4")[2][end]==1002
        @test _PEL.mesh.get_keys_for_element(1002,"Lagrange")[2]==_PEL.mesh.get_element(1002)[2]
    finally
        _PEL.finalize()
    end
end

@testset "Native Mesh and Tet10 overlays rename across stored dimensions" begin
    for quadratic in (false,true)
        try
            _PEL.initialize();_PEL.option("Mesh.Renumber",0)
            for (dim,tag) in ((1,11),(2,22),(3,33));_PEL.model.add_discrete_entity(dim,tag);end
            mesh=Mesh(_PEL_COORDS[1];segs=reshape(Int32[1,2],2,1),
                tris=reshape(Int32[1,2,3],3,1),tets=reshape(Int32[1,2,3,4],4,1))
            entities=[(1,Int32(11)),(2,Int32(22)),(3,Int32(33))]
            owners=[(1,Int32(11)),(1,Int32(11)),(2,Int32(22)),(3,Int32(33))]
            class=_PEL._MeshClassification(mesh,last(entities),entities,owners,
                Dict{Tuple{Int,Int32},Vector{Int32}}(),Int32[11],Int32[22],Int32[33])
            _PEL._replace_mesh_cache_locked!(mesh,class)
            quadratic && _PEL.mesh.set_order(2)
            _PEL.mesh.renumber_elements([1,2,3],[1001,7007,9009])
            before=_pel_state();scalar=Dict(tag=>_PEL.mesh.get_element(tag) for tag in (1001,7007,9009))
            _PEL.mesh.set_visibility([1001,9009],0)
            _PEL.mesh.renumber_elements([1001,7007,9009],[7007,9009,1001])
            @test _PEL.LAST_MESH[]===before.mesh
            @test _PEL.LAST_MESH_HIGH_ORDER[]===before.overlay
            @test _PEL.LAST_MESH_HIGH_ORDER_MIDS[]===before.mids
            @test _PEL.mesh.get_nodes()==before.nodes
            @test _PEL.mesh.get_elements()[1]==Int32[1,2,quadratic ? 11 : 4]
            @test _PEL.mesh.get_elements()[2]==[UInt64[7007],UInt64[9009],UInt64[1001]]
            @test _PEL.mesh.get_elements()[3]==before.elements[3]
            @test _PEL.mesh.get_visibility([7007,9009,1001])==Int32[0,1,0]
            for (source,target) in zip((1001,7007,9009),(7007,9009,1001))
                @test _PEL.mesh.get_element(target)==scalar[source]
            end
            @test _PEL.mesh.get_element_by_coordinates(.1,.1,.1,3,true)[1]==1001
            @test _PEL.mesh.get_element_qualities([1001],"volume")[1]≈1/6 atol=2e-11
            @test _PEL.mesh.get_max_element_tag()==9009
        finally
            _PEL.finalize()
        end
    end
end

@testset "Native label record additions preserve overlays and reject collisions atomically" begin
    for quadratic in (false,true)
        try
            _pel_native_square(;quadratic)
            old=vcat(_PEL.mesh.get_elements()[2]...)
            _PEL.mesh.renumber_elements(old,old.+9000)
            cache=_PEL.LAST_MESH[];overlay=_PEL.LAST_MESH_HIGH_ORDER[]
            node_max=_PEL.mesh.get_max_node_tag();element_max=_PEL.mesh.get_max_element_tag()
            _PEL.mesh.add_nodes(2,1,[],[0.,0,0,1,0,0,0,1,0])
            @test _PEL.LAST_MESH[]===cache && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
            added=UInt64.(node_max.+(1:3))
            _PEL.mesh.add_elements_by_type(1,2,[],added)
            @test _PEL.LAST_MESH[]===cache && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
            @test _PEL.mesh.get_element(element_max+1)==(Int32(2),added,2,1)
            @test element_max+1 in vcat(_PEL.mesh.get_elements(2,1)[2]...)
            @test isempty(_PEL.LAST_MESH_CLASS[].public_tags.authority)
            snapshot=_pel_state()
            @test_throws ArgumentError _PEL.mesh.add_nodes(2,1,[_PEL.mesh.get_nodes()[1][1]],[0.,0,0])
            _pel_unchanged(snapshot)
            @test_throws ArgumentError _PEL.mesh.add_elements_by_type(1,2,[old[1]+9000],added)
            _pel_unchanged(snapshot)
            # The explicitly attached cell participates in the same detached
            # partial/default mapping, while its coordinates stay separate.
            _PEL.mesh.renumber_elements([element_max+1],[77])
            @test _PEL.mesh.get_element(77)==(Int32(2),added,2,1)
            @test _PEL.LAST_MESH[]===cache && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
            @test length(_PEL.mesh.get_nodes()[1])==node_max+3
        finally
            _PEL.finalize()
        end
    end
end

@testset "Uncovered records survive native label order and refinement allocations" begin
    try
        _pel_native_square()
        old=vcat(_PEL.mesh.get_elements()[2]...)
        _PEL.mesh.renumber_elements(old,old.+3000)
        _PEL.mesh.add_nodes(2,1,[80,81,82],[0.,0,0,1,0,0,0,1,0])
        _PEL.mesh.add_elements_by_type(1,2,[700],[80,81,82])
        record=deepcopy(_PEL.CURRENT[].meshing.attached[(2,1)])
        maximum_node=_PEL.mesh.get_max_node_tag()
        _PEL.mesh.set_order(2)
        native=_PEL.LAST_MESH_CLASS[].public_tags
        @test minimum(native.node_tags[10:end])==maximum_node+1
        @test _PEL.mesh.get_node(80)[1]==[0.,0,0]
        @test _PEL.mesh.get_element(700)==(Int32(2),UInt64[80,81,82],2,1)
        @test _PEL.CURRENT[].meshing.attached[(2,1)].node_tags==record.node_tags
        @test _PEL.CURRENT[].meshing.attached[(2,1)].element_nodes==record.element_nodes
        @test _PEL.mesh.get_element_types(2,1)==Int32[2,9]
        element_max=_PEL.mesh.get_max_element_tag()
        _PEL.mesh.refine()
        @test _PEL.LAST_MESH[] isa Mesh && _PEL.LAST_MESH_HIGH_ORDER[] isa P2TriMesh
        @test minimum(_PEL.LAST_MESH_CLASS[].public_tags.element_tags)==element_max+1
        @test _PEL.mesh.get_element(700)==(Int32(2),UInt64[80,81,82],2,1)
        @test _PEL.mesh.get_node(80)[1]==[0.,0,0]
        snapshot=_pel_state()
        @test_throws ArgumentError _PEL.mesh.remove_elements(2,1,[700,999999])
        _pel_unchanged(snapshot)
        cache_tag=first(_PEL.LAST_MESH_CLASS[].public_tags.element_tags)
        _PEL.mesh.remove_elements(2,1,[cache_tag,700])
        @test !(700 in vcat(_PEL.mesh.get_elements()[2]...))
        @test !(cache_tag in vcat(_PEL.mesh.get_elements()[2]...))
        @test _PEL.mesh.get_node(80)[1]==[0.,0,0]
    finally
        _PEL.finalize()
    end
    try
        _pel_mixed()
        _PEL.mesh.add_nodes(3,11,[80,81,82,83],vec(_PEL_COORDS[1].+[10,0,0]))
        _PEL.mesh.add_elements_by_type(11,4,[800],[80,81,82,83])
        before=_PEL.mesh.get_element(800)
        element_max=_PEL.mesh.get_max_element_tag()
        _PEL.mesh.refine()
        @test _PEL.LAST_MESH[] isa MixedMesh
        @test all(block->Tessella.Elements.msh_spec(block.msh).order==1,_PEL.LAST_MESH[].blocks)
        @test minimum(_PEL.LAST_MESH_CLASS[].public_tags.element_tags)==element_max+1
        @test _PEL.mesh.get_element(800)==before
        @test isempty(_PEL.LAST_MESH_CLASS[].public_tags.authority)
    finally
        _PEL.finalize()
    end
end

function _pel_scale_fixture(count)
    _PEL.initialize()
    _PEL.model.add_discrete_entity(2,1)
    coordinates=zeros(3,3count)
    for column in 1:count
        coordinates[:,3column-2:3column]=Float64[3column 3column+1 3column;0 0 1;0 0 0]
    end
    mesh=Mesh(coordinates;tris=reshape(Int32.(1:3count),3,count))
    class=_PEL._MeshClassification(mesh,(2,Int32(1)),[(2,Int32(1))],
        fill((2,Int32(1)),3count),Dict{Tuple{Int,Int32},Vector{Int32}}(),
        Int32[],fill(Int32(1),count),Int32[])
    _PEL._replace_mesh_cache_locked!(mesh,class)
    return mesh
end

@testset "Native periodic nodes and keys retain sparse public labels" begin
    try
        _PEL.initialize()
        for (x,y) in ((0,0),(1,0),(0,1),(1,1));_PEL.model.add_point(x,y,0);end
        for (a,b) in ((1,2),(3,4),(2,4),(1,3));_PEL.model.add_line(a,b);end
        _PEL.model.add_curve_loop([1,3,-2,-4]);_PEL.model.add_plane_surface([1])
        for curve in 1:4;_PEL.mesh.set_transfinite_curve(curve,3);end
        _PEL.mesh.set_transfinite_surface(1)
        affine=[1.,0,0,0,0,1,0,1,0,0,1,0]
        _PEL.mesh.set_periodic(1,[2],[1],affine)
        _PEL.mesh.generate(2)
        before=_PEL.mesh.get_periodic_nodes(1,2)
        @test length(before.slave_nodes)==3
        for (slave,master) in zip(before.slave_nodes,before.master_nodes)
            @test _PEL.mesh.get_node(slave)[1]≈_PEL.mesh.get_node(master)[1].+[0,1,0] atol=2e-11
        end
        slave_nodes=sort!(copy(before.slave_nodes);by=node->_PEL.mesh.get_node(node)[1][1])
        _PEL.mesh.add_elements_by_type(2,1,[501,502],slave_nodes[[1,2,2,3]])
        keys=_PEL.mesh.get_periodic_keys(1,"Lagrange",2)
        @test !isempty(keys[4])
        elements=vcat(_PEL.mesh.get_elements()[2]...)
        _PEL.mesh.renumber_elements(elements,elements.+1000)
        nodes=_PEL.mesh.get_nodes()[1];labels=reverse(nodes).+2000
        mapping=Dict(zip(nodes,labels))
        _PEL.mesh.renumber_nodes(nodes,labels)
        after=_PEL.mesh.get_periodic_nodes(1,2)
        @test after.slave_nodes==[mapping[tag] for tag in before.slave_nodes]
        @test after.master_nodes==[mapping[tag] for tag in before.master_nodes]
        renamed=_PEL.mesh.get_periodic_keys(1,"Lagrange",2)
        @test renamed[4]==[mapping[tag] for tag in keys[4]]
        @test renamed[5]==[mapping[tag] for tag in keys[5]]
        @test renamed[6:7]==keys[6:7]
    finally
        _PEL.finalize()
    end
end

@testset "Native actual midpoint labels survive node and element deduplication" begin
    try
        _PEL.initialize();_PEL.option("Mesh.Renumber",0)
        for owner in (1,2);_PEL.model.add_discrete_entity(3,owner);end
        coordinates=hcat(_PEL_COORDS[1],_PEL_COORDS[1][:,1])
        mesh=Mesh(coordinates;tets=Int32[1 5;2 2;3 3;4 4])
        class=_PEL._MeshClassification(mesh,(3,Int32(1)),[(3,Int32(1))],
            fill((3,Int32(1)),5),Dict{Tuple{Int,Int32},Vector{Int32}}(),
            Int32[],Int32[],Int32[1,1])
        _PEL._replace_mesh_cache_locked!(mesh,class)
        _PEL.mesh.renumber_elements([1,2],[700,800]);_PEL.mesh.set_order(2)
        original_nodes=_PEL.mesh.get_nodes()[1]
        labels=reverse(original_nodes).+1000
        _PEL.mesh.renumber_nodes(original_nodes,labels)
        before=_PEL.LAST_MESH_HIGH_ORDER[]
        slots=Tessella.HighOrder._P2_EDGE_SLOTS
        edge_slot=only(slot for (slot,i,j) in slots if minmax(i,j)==(1,2))
        first_mid=labels[before.tet10[edge_slot,1]]
        duplicate_mid=labels[before.tet10[edge_slot,2]]
        raw=[duplicate_mid,labels[3],labels[4],labels[1]]
        _PEL.mesh.add_elements_by_type(2,4,[900],raw)
        duplicate_tags=_PEL.mesh.get_duplicate_nodes()
        @test labels[1] in duplicate_tags && labels[5] in duplicate_tags
        @test first_mid in duplicate_tags && duplicate_mid in duplicate_tags
        coordinate=copy(before.coords[:,before.tet10[edge_slot,2]])
        before.coords[3,before.tet10[edge_slot,2]]+=.125
        snapshot=_pel_state()
        @test_throws ArgumentError _PEL.mesh.remove_duplicate_nodes()
        _pel_unchanged(snapshot)
        before.coords[:,before.tet10[edge_slot,2]]=coordinate
        before.coords[:,before.tet10[edge_slot,2]]=before.coords[:,1]
        snapshot=_pel_state()
        @test_throws ArgumentError _PEL.mesh.remove_duplicate_nodes()
        _pel_unchanged(snapshot)
        before.coords[:,before.tet10[edge_slot,2]]=coordinate
        # Gmsh 4.15.2 GModel.cpp removeDuplicateMeshVertices deliberately
        # ignores getNum(): the first stored node survives each coincidence.
        _PEL.mesh.remove_duplicate_nodes()
        @test _PEL.mesh.get_element(900)[2]==UInt64[first_mid,labels[3],labels[4],labels[1]]
        @test _PEL.mesh.get_nodes()[1][1:4]==labels[1:4]
        @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===_PEL.LAST_MESH[]
        @test length(_PEL.LAST_MESH_CLASS[].public_tags.node_tags)==10
        @test length(unique(_PEL.mesh.get_nodes()[1]))==10
        cells=_PEL.mesh.get_elements_by_type(11)
        @test cells[1]==UInt64[801,802]
        @test reshape(cells[2],10,:)[:,1]==reshape(cells[2],10,:)[:,2]
        _PEL.mesh.set_visibility([801],0)
        _PEL.mesh.remove_duplicate_elements()
        @test _PEL.mesh.get_elements_by_type(11)[1]==UInt64[801]
        @test _PEL.mesh.get_visibility([801])==Int32[0]
        @test _PEL.mesh.get_element_qualities([801],"volume")[1]≈1/6 atol=2e-11
        @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===_PEL.LAST_MESH[]
        _PEL.mesh.add_nodes(3,2,[3000],_PEL.mesh.get_node(first_mid)[1])
        _PEL.mesh.add_elements_by_type(2,4,[901],[3000,labels[3],labels[4],labels[1]])
        @test 3000 in _PEL.mesh.get_duplicate_nodes()
        _PEL.mesh.remove_duplicate_nodes()
        @test _PEL.mesh.get_element(901)[2]==UInt64[first_mid,labels[3],labels[4],labels[1]]
        @test !(3000 in _PEL.mesh.get_nodes()[1])
        @test _PEL.mesh.get_max_node_tag()==3000
    finally
        _PEL.finalize()
    end
end

@testset "Native regeneration rejects unresolved cache references atomically" begin
    try
        _pel_native_square(;quadratic=true)
        old=vcat(_PEL.mesh.get_elements()[2]...)
        _PEL.mesh.renumber_elements(old,old.+1000)
        _PEL.option("Mesh.ElementOrder",2)
        _PEL.model.add_discrete_entity(1,90)
        nodes=_PEL.mesh.get_nodes()[1]
        _PEL.mesh.add_elements_by_type(90,1,[900],nodes[[1,10]])
        for renumber in (0,1)
            _PEL.option("Mesh.Renumber",renumber)
            snapshot=_pel_state()
            @test_throws ArgumentError _PEL.mesh.generate(2)
            _pel_unchanged(snapshot)
            @test _PEL.mesh.get_element(900)[2]==nodes[[1,10]]
        end
        _PEL.mesh.remove_elements(1,90,[900])
        _PEL.option("Mesh.Renumber",0)
        # Independent Gmsh 4.15.2 Point-record allocator proof:
        # public_labels_square_point_records_primary_v2.json, SHA256
        # 5128641775becb0b2d07f0e850a1908f93edd25118b16027a350393c7c4d8802.
        # Native Curve phase retention and global order remain a separate parity fixture.
        point_labels=Dict(point=>_PEL.mesh.get_nodes(0,point)[1] for point in 1:4)
        for (point,node,cell,x) in ((91,6001,901,10.),(92,6002,902,11.))
            _PEL.model.add_discrete_entity(0,point)
            _PEL.mesh.add_nodes(0,point,[node],[x,0,0])
            _PEL.mesh.add_elements_by_type(point,15,[cell],[node])
        end
        records=Dict(cell=>_PEL.mesh.get_element(cell) for cell in (901,902))
        node_highwater=_PEL.mesh.get_max_node_tag();element_highwater=_PEL.mesh.get_max_element_tag()
        _PEL.mesh.generate(2)
        @test _PEL.LAST_MESH[] isa Mesh && _PEL.LAST_MESH_HIGH_ORDER[] isa P2TriMesh
        public=_PEL.LAST_MESH_CLASS[].public_tags
        @test public.mesh===_PEL.LAST_MESH[]
        @test Dict(point=>_PEL.mesh.get_nodes(0,point)[1] for point in 1:4)==point_labels
        @test minimum(filter(tag->tag>node_highwater,public.node_tags))==node_highwater+1
        @test all(tag->tag>element_highwater,public.element_tags)
        @test Dict(cell=>_PEL.mesh.get_element(cell) for cell in (901,902))==records
        @test _PEL.mesh.get_node(6001)[1]==[10.,0,0]
        @test _PEL.mesh.get_node(6002)[1]==[11.,0,0]
        all_nodes=_PEL.mesh.get_nodes()[1]
        _,all_cells,all_connections=_PEL.mesh.get_elements()
        @test length(unique(all_nodes))==length(all_nodes)
        @test length(unique(vcat(all_cells...)))==sum(length,all_cells)
        known_nodes=Set(all_nodes)
        @test all(node->node in known_nodes,Iterators.flatten(all_connections))
        @test sum(_PEL.mesh.get_element_qualities(public.element_tags,"volume"))≈1 atol=2e-11
    finally
        _PEL.finalize()
    end
end

@testset "Virgin native caches accept node-first public labels" begin
    for quadratic in (false,true)
        try
            _pel_native_square(;quadratic)
            mesh=_PEL.LAST_MESH[];overlay=_PEL.LAST_MESH_HIGH_ORDER[]
            mids=_PEL.LAST_MESH_HIGH_ORDER_MIDS[]
            @test _PEL.LAST_MESH_CLASS[].public_tags===nothing
            nodes,coords,params=_PEL.mesh.get_nodes()
            elements=_PEL.mesh.get_elements()
            labels=UInt64.(10003 .+ 7 .* (0:length(nodes)-1))
            records=[_PEL.mesh.get_node(node) for node in nodes]
            _PEL.mesh.renumber_nodes(nodes,labels)
            @test _PEL.LAST_MESH[]===mesh && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
            @test _PEL.LAST_MESH_HIGH_ORDER_MIDS[]===mids
            @test _PEL.mesh.get_nodes()==(labels,coords,params)
            @test _PEL.mesh.get_elements()[2]==elements[2]
            @test [_PEL.mesh.get_node(label) for label in labels]==records
            @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===mesh
            _PEL.mesh.renumber_nodes(labels,reverse(labels))
            @test _PEL.LAST_MESH[]===mesh && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
            @test _PEL.mesh.get_nodes()==(reverse(labels),coords,params)
            @test [_PEL.mesh.get_node(label) for label in reverse(labels)]==records
            _PEL.mesh.set_order(quadratic ? 2 : 1)
            @test _PEL.LAST_MESH[]===mesh && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
        finally
            _PEL.finalize()
        end
    end
end

@testset "Native Point source labels survive higher generation and order" begin
    # Independent Gmsh 4.15.2 authority: b4_free_api01_point_primary_v2.json,
    # SHA256 68a54a0696fb4116befdfe615eb4285d7e383b3f151f46d15f6cb909101f53c1.
    for renumber in (0,1)
        _PEL.initialize()
        try
            _PEL.option("Mesh.Renumber",renumber);_PEL.option("Mesh.ElementOrder",1)
            _PEL.model.add_box(0,0,0,1,1,1;tag=1)
            point=_PEL.model.get_value(0,1,Float64[])
            @test point==[0.,0.,1.]
            _PEL.mesh.add_nodes(0,1,[111],point)
            _PEL.mesh.add_elements_by_type(1,15,[777],[111])
            for (_,curve) in _PEL.model.get_entities(1);_PEL.mesh.set_transfinite_curve(curve,2);end
            for (_,surface) in _PEL.model.get_entities(2)
                _PEL.mesh.set_transfinite_surface(surface);_PEL.mesh.set_recombine(2,surface)
            end
            _PEL.mesh.set_transfinite_volume(1)
            node=UInt64(renumber==0 ? 111 : 1);element=UInt64(renumber==0 ? 777 : 1)
            for (step,volume_type,expected_nodes) in ((:generate3,5,8),(:order2,12,27),
                    (:order1,5,8),(:generate3,5,8),(:generate2,0,8),(:generate3,5,8))
                if step===:order2
                    _PEL.mesh.set_order(2)
                elseif step===:order1
                    _PEL.mesh.set_order(1)
                else
                    _PEL.mesh.generate(step===:generate2 ? 2 : 3)
                end
                @test _PEL.LAST_MESH[] isa MixedMesh
                @test _PEL.LAST_MESH_CLASS[].public_tags.mesh===_PEL.LAST_MESH[]
                @test _PEL.mesh.get_nodes(0,1,true,true)==([node],point,Float64[])
                @test _PEL.mesh.get_elements(0,1)==(Int32[15],[[element]],[[node]])
                @test _PEL.mesh.get_node(node)==(point,Float64[],0,1)
                @test _PEL.mesh.get_element(element)==(Int32(15),[node],0,1)
                @test _PEL.mesh.get_nodes_by_element_type(15,1,true)==([node],point,Float64[])
                @test _PEL.mesh.get_barycenters(15,1,false,true)==point
                published=_PEL.LAST_MESH[];overlay=_PEL.LAST_MESH_HIGH_ORDER[]
                @test _PEL.mesh.get_local_coordinates_in_element(element,point...)==(0.,0.,0.)
                @test _PEL.mesh.get_element_by_coordinates(point...,0,true)==
                    (element,Int32(15),[node],0.,0.,0.)
                @test _PEL.mesh.get_elements_by_coordinates(point...,0,true)==[element]
                jac=_PEL.mesh.get_jacobian(element,[0.,0.,0.])
                @test jac[1]==[1.,0,0,0,1,0,0,0,1] && jac[2]==[1.] && jac[3]==point
                @test _PEL.mesh.get_jacobians(15,[0.,0.,0.],1)==jac
                @test _PEL.mesh.get_basis_functions_orientation_for_element(element,"Lagrange")==0
                @test _PEL.mesh.get_basis_functions_orientation(15,"Lagrange",1)==Int32[0]
                @test _PEL.mesh.get_keys_for_element(element,"Lagrange",true)==(Int32[0],[node],point)
                @test _PEL.mesh.get_keys(15,"Lagrange",1,true)==(Int32[0],[node],point)
                @test _PEL.LAST_MESH[]===published && _PEL.LAST_MESH_HIGH_ORDER[]===overlay
                tags,xyz,params=_PEL.mesh.get_nodes()
                @test length(tags)==expected_nodes && allunique(tags) && count(==(node),tags)==1
                @test isempty(params) && length(xyz)==3expected_nodes
                types,ids,connectivity=_PEL.mesh.get_elements()
                @test allunique(vcat(ids...)) && all(n in Set(tags) for nodes in connectivity for n in nodes)
                @test _PEL.mesh.get_element_types(3,1)==(volume_type==0 ? Int32[] : Int32[volume_type])
            end
            for (_,curve) in _PEL.model.get_entities(1);_PEL.mesh.set_transfinite_curve(curve,3);end
            for dimension in (1,0,3)
                _PEL.mesh.generate(dimension)
                @test _PEL.mesh.get_nodes(0,1,true,true)==([node],point,Float64[])
                @test _PEL.mesh.get_elements(0,1)==(Int32[15],[[element]],[[node]])
                @test _PEL.mesh.get_nodes_by_element_type(15,1,true)==([node],point,Float64[])
                @test _PEL.mesh.get_keys_for_element(element,"Lagrange",true)==(Int32[0],[node],point)
                if dimension<=1
                    @test isempty(_PEL.mesh.get_elements(3,1)[1])
                    @test length(_PEL.mesh.get_nodes()[1])==20
                else
                    @test _PEL.mesh.get_element_types(3,1)==Int32[5]
                end
            end
        finally
            _PEL.finalize()
        end
    end
end

@testset "Native labeled node removal keeps surviving record references valid" begin
    try
        _pel_native_square(;quadratic=true)
        old=vcat(_PEL.mesh.get_elements()[2]...);_PEL.mesh.renumber_elements(old,old.+1000)
        triangle,connectivity=_PEL.mesh.get_elements_by_type(9,1)
        _PEL.mesh.add_elements_by_type(1,9,[7000],connectivity[1:6])
        @test _PEL.mesh.get_nodes_by_element_type(9,1,false)[2][end-17:end]==
            vcat((_PEL.mesh.get_node(tag)[1] for tag in connectivity[1:6])...)
        @test _PEL.mesh.get_barycenters(9,1,false,false)[end-2:end]≈
            sum((_PEL.mesh.get_node(tag)[1] for tag in connectivity[1:6]))/6 atol=2e-11
        state=_pel_state()
        @test_throws ArgumentError _PEL.mesh.set_order(1)
        _pel_unchanged(state)
        _PEL.mesh.remove_elements(2,1,[7000]);_PEL.mesh.set_order(1)
        @test _PEL.LAST_MESH_HIGH_ORDER[]===nothing
    finally
        _PEL.finalize()
    end
    try
        _pel_native_square()
        triangle,connectivity=_PEL.mesh.get_elements_by_type(2,1)
        state=_pel_state()
        @test_throws ArgumentError _PEL.mesh.add_elements_by_type(1,2,[triangle[1]],connectivity[1:3])
        _pel_unchanged(state)
    finally
        _PEL.finalize()
    end
    try
        _pel_native_square(;quadratic=true)
        _,nodes=_PEL.mesh.get_elements_by_type(9,1)
        @test _PEL.LAST_MESH_CLASS[].public_tags===nothing
        _PEL.mesh.add_elements_by_type(1,9,[7000],nodes[1:6])
        @test _PEL.LAST_MESH_CLASS[].public_tags!==nothing
        state=_pel_state()
        @test_throws ArgumentError _PEL.mesh.set_order(1)
        _pel_unchanged(state)
    finally
        _PEL.finalize()
    end
end

@testset "Native cached Point labels survive order conversion and refinement" begin
    _PEL.initialize()
    try
        _PEL.option("Mesh.Renumber",0)
        _PEL.model.add_discrete_entity(0,91);_PEL.model.add_discrete_entity(0,93)
        _PEL.model.add_discrete_entity(2,92)
        mesh=MixedMesh(Float64[3 4 0 1 0;3 3 0 0 1;0 0 0 0 0],
            [ElementBlock(15,reshape(Int32[1,2],1,2)),ElementBlock(2,reshape(Int32[3,4,5],3,1))])
        class=_PEL._mixed_classification(mesh,(2,Int32(92)),[(2,Int32(92)),(0,Int32(91)),(0,Int32(93))],
            [(0,Int32(91)),(0,Int32(93)),(2,Int32(92)),(2,Int32(92)),(2,Int32(92))],
            Dict{Tuple{Int,Int32},Vector{Int32}}(),Dict(15=>Int32[91,93],2=>Int32[92]))
        _PEL._replace_mesh_cache_locked!(mesh,class)
        _PEL.mesh.renumber_elements([1,2,3],[777,778,888])
        for operation in (()->_PEL.mesh.set_order(2),()->_PEL.mesh.set_order(1),()->_PEL.mesh.refine())
            operation()
            @test _PEL.LAST_MESH[] isa MixedMesh && _PEL.LAST_MESH[].entity_data===nothing
            @test _PEL.mesh.get_element(777)==(Int32(15),UInt64[1],0,91)
            @test _PEL.mesh.get_element(778)==(Int32(15),UInt64[2],0,93)
            @test _PEL.mesh.get_elements_by_type(15)==(UInt64[777,778],UInt64[1,2])
            @test _PEL.mesh.get_nodes(0,91)[2]==[3.,3,0]
        end
    finally
        _PEL.finalize()
    end
end

@testset "Generation high-water includes retained independent records" begin
    _PEL.initialize()
    try
        _PEL.option("Mesh.Renumber",1);_PEL.model.add_box(0,0,0,1,1,1;tag=1)
        _PEL.mesh.add_nodes(0,1,[1],[0.,0,1])
        _PEL.model.add_discrete_entity(3,99)
        _PEL.mesh.add_nodes(3,99,[2,3,4,5],[10.,0,0,11,0,0,10,1,0,10,0,1])
        _PEL.mesh.add_elements_by_type(99,4,[1],[2,3,4,5])
        for (_,curve) in _PEL.model.get_entities(1);_PEL.mesh.set_transfinite_curve(curve,2);end
        for (_,surface) in _PEL.model.get_entities(2)
            _PEL.mesh.set_transfinite_surface(surface);_PEL.mesh.set_recombine(2,surface)
        end
        _PEL.mesh.set_transfinite_volume(1);_PEL.mesh.generate(3)
        node_max=_PEL.mesh.get_max_node_tag();element_max=_PEL.mesh.get_max_element_tag()
        @test _PEL.NODE_TAG_MAX[]==node_max && _PEL.ELEMENT_TAG_MAX[]==element_max
        @test node_max==12 && element_max==2
        _PEL.mesh.add_nodes(3,99,UInt64[],[12.,2,2])
        @test _PEL.mesh.get_node(node_max+1)[1]==[12.,2,2]
        _,connectivity=_PEL.mesh.get_elements_by_type(4,99)
        _PEL.mesh.add_elements_by_type(99,4,UInt64[],connectivity)
        @test _PEL.mesh.get_element(element_max+1)[2]==connectivity
        _PEL.mesh.set_order(2)
        tags=_PEL.mesh.get_nodes()[1]
        @test allunique(tags)
        @test _PEL.mesh.get_node(node_max+1)[1]==[12.,2,2]
    finally
        _PEL.finalize()
    end
end

@testset "Pure labels preserve populated metadata and empty blocks" begin
    try
        _pel_mixed(;tagged=true)
        original=_PEL.LAST_MESH[];data=original.entity_data
        blocks=original.blocks
        metadata=Tessella.Elements.MixedEntityData(data.entities;
            node_entities=data.node_entities,node_parametric=[collect(original.coords[:,i]) for i in axes(original.coords,2)],
            external_node_tags=data.external_node_tags,
            block_entities=data.block_entities,external_element_tags=data.external_element_tags)
        ancillary=Tessella.Elements.MshAncillarySection("Comments",false,UInt8[65,10],Int32(1))
        view=Tessella.Elements.MshDataSection("NodeData",["field"],[0.],[0,1,1],
            Int32[],Int32[1],Int32[],false,[3.],Int32(1))
        periodic=Tessella.Elements.MixedPeriodicLink(3,22,11,[5],[1])
        mesh=MixedMesh(original.coords,blocks;entity_data=metadata,
            physical_names=Dict((3,42)=>"material"),periodic_links=[periodic],
            ancillary_sections=[ancillary],data_sections=[view])
        # Existing storage can be edited after construction. A label operation
        # retains that zero-width block without reconstructing or validating it.
        push!(mesh.blocks,ElementBlock(2,zeros(Int32,3,0)))
        push!(mesh.entity_data.block_entities,Int32[])
        push!(mesh.entity_data.external_element_tags,UInt64[])
        mesh.elementary_entities===nothing || length(mesh.elementary_entities)==4 || push!(mesh.elementary_entities,Int32[])
        class=_PEL._mixed_rebind_class(_PEL.LAST_MESH_CLASS[],mesh)
        _PEL._replace_mesh_cache_locked!(mesh,class)
        before=(coords=copy(mesh.coords),blocks=[(b.msh,copy(b.nodes),copy(b.tags)) for b in mesh.blocks],
            owners=deepcopy(mesh.entity_data),names=copy(mesh.physical_names),
            periodic=repr(mesh.periodic_links),ancillary=repr(mesh.ancillary_sections),views=repr(mesh.data_sections))
        _PEL.mesh.renumber_elements([1001,7007,9009],[7007,9009,1001])
        result=_PEL.LAST_MESH[]
        @test result.coords==before.coords
        @test [(b.msh,b.nodes,b.tags) for b in result.blocks]==before.blocks
        @test result.entity_data.external_node_tags==before.owners.external_node_tags
        @test result.entity_data.node_entities==before.owners.node_entities
        @test result.entity_data.node_parametric==before.owners.node_parametric
        @test result.entity_data.block_entities==before.owners.block_entities
        @test repr(result.entity_data.entities)==repr(before.owners.entities)
        @test result.physical_names==before.names
        @test repr(result.periodic_links)==before.periodic
        @test repr(result.ancillary_sections)==before.ancillary
        @test repr(result.data_sections)==before.views
        @test result.entity_data.external_element_tags==[UInt64[7007],UInt64[9009],UInt64[1001],UInt64[]]
    finally
        _PEL.finalize()
    end
end

@testset "Public element label plans have bounded linear retained storage" begin
    allocations=Int[];retained=Int[]
    try
        for count in (1000,2000,4000)
            mesh=_pel_scale_fixture(count)
            _PEL.mesh.renumber_elements()
            _PEL.mesh.renumber_elements([1],[10001])
            _PEL.mesh.renumber_elements()
            push!(allocations,minimum(@allocated(_PEL.mesh.renumber_elements()) for _ in 1:3))
            push!(retained,Base.summarysize(_PEL.LAST_MESH_CLASS[].public_tags))
            @test _PEL.LAST_MESH[]===mesh
            @test _PEL.LAST_MESH_CLASS[].public_tags.element_tags==UInt64.(1:count)
            @test _PEL.mesh.get_element_qualities([1,count],"volume")==[.5,.5]
            @test _PEL.mesh.get_max_element_tag()==10000+count
        end
        for index in 2:3
            @test allocations[index]<=2.5allocations[index-1]+262144
            @test retained[index]<=2.5retained[index-1]+262144
        end
        @info "Public element label storage" cells=(1000,2000,4000) allocations=Tuple(allocations) retained=Tuple(retained)
    finally
        _PEL.finalize()
    end
end
