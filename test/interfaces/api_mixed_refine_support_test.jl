using Test,Tessella

function _mrs_native_prism(api,offset,x,z)
    points=((x,0.,z),(x+1,0.,z),(x,1.,z),
            (x,0.,z+1),(x+1,0.,z+1),(x,1.,z+1))
    for (i,p) in enumerate(points)
        api.model.add_point(p...;tag=offset+i)
    end
    for (i,(a,b)) in enumerate(((1,2),(2,3),(1,3),(4,5),(5,6),(4,6),
                                (1,4),(2,5),(3,6)))
        api.model.add_line(offset+a,offset+b;tag=offset+i)
        api.mesh.set_transfinite_curve(offset+i,2)
    end
    for (i,cycle) in enumerate(((1,8,-4,-7),(2,9,-5,-8),(3,9,-6,-7),
                                (1,2,-3),(4,5,-6)))
        api.model.add_curve_loop([sign(c)*(offset+abs(c)) for c in cycle];tag=offset+i)
        api.model.add_plane_surface([offset+i];tag=offset+i)
    end
    api.model.add_surface_loop([offset+1,offset+2,-offset-3,-offset-4,offset+5];tag=offset+1)
    api.model.add_volume([offset+1];tag=offset+1)
    for (i,corners) in enumerate(((1,2,5,4),(2,3,6,5),(1,3,6,4),(1,2,3),(4,5,6)))
        api.mesh.set_transfinite_surface(offset+i,"Left",[offset+p for p in corners])
        i<=3 && api.mesh.set_recombine(2,offset+i)
    end
    api.mesh.set_transfinite_volume(offset+1,collect(offset+1:offset+6))
end

function _mrs_native_fixture(api,family)
    if family==6
        _mrs_native_prism(api,0,0.,0.)
        _mrs_native_prism(api,100,1.,1.)
    else
        api.model.add_box(0,0,0,1,1,1)
        api.model.add_box(1,1,1,1,1,1)
        for (_,curve) in api.model.get_entities(1)
            api.mesh.set_transfinite_curve(curve,2)
        end
        for (_,surface) in api.model.get_entities(2)
            api.mesh.set_transfinite_surface(surface)
            api.mesh.set_recombine(2,surface)
        end
        for (_,volume) in api.model.get_entities(3)
            api.mesh.set_transfinite_volume(volume)
        end
    end
    api.mesh.generate(3)
end

# Complete imported cells provide actual lower carriers even when an older
# classification has no support catalog. Sparse node labels are deliberately
# opposite to source storage order, and the two bodies share only one point.
function _mrs_tagged_fixture(api,family)
    reference=Tessella.Elements.lagrange_nodes(family)
    edges=Tessella.API._api_support_edges(Tessella.Elements.msh_family(family))
    faces=Tessella.Elements._VOLUME_CELL_FACES[family]
    for region in 0:1
        offset=100region
        points=family==6 ? reference.*[1.,1.,.5] .+[region,0.,region+.5] :
                          reference.*.5 .+(region+.5)
        n=size(points,2)
        labels=UInt64[1000+17*(2n-(region*n+i)) for i in 1:n]
        for i in 1:n
            api.model.add_discrete_entity(0,offset+i)
            api.mesh.add_nodes(0,offset+i,[labels[i]],points[:,i])
            # Only one Point cell will reference the merged vertex. Duplicate
            # Point15 cells would make the imported refinement input invalid.
            if !(region==1 && i==1)
                api.mesh.add_elements_by_type(offset+i,15,[10000+offset+i],[labels[i]])
            end
        end
        for (i,(a,b)) in enumerate(edges)
            api.model.add_discrete_entity(1,offset+i,[offset+a+1,offset+b+1])
            api.mesh.add_elements_by_type(offset+i,1,[20000+offset+i],[labels[a+1],labels[b+1]])
        end
        for (i,face) in enumerate(faces)
            boundary=Int[]
            for j in eachindex(face)
                pair=minmax(face[j]-1,face[mod1(j+1,length(face))]-1)
                edge=findfirst(e->minmax(e...)==pair,edges)
                @assert edge!==nothing
                push!(boundary,offset+edge)
            end
            api.model.add_discrete_entity(2,offset+i,boundary)
            api.mesh.add_elements_by_type(offset+i,length(face)==3 ? 2 : 3,
                [30000+offset+i],labels[collect(face)])
        end
        api.model.add_discrete_entity(3,offset+1,collect(offset+1:offset+length(faces)))
        api.mesh.add_elements_by_type(offset+1,family,[40000+offset+1],labels)
    end
    api.option("Mesh.Renumber",0)
    api.mesh.generate(0)
    @test api.LAST_MESH_CLASS[].public_tags!==nothing
end

function _mrs_assert_actual_supports(api)
    mesh=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
    edge_keys=Set{NTuple{2,Int32}}()
    face_keys=Set{NTuple{3,Int32}}()
    quad_keys=Set{NTuple{4,Int32}}()
    for block in mesh.blocks,column in axes(block.nodes,2)
        family=Tessella.Elements.msh_family(block.msh)
        for (i,j) in api._api_support_edges(family)
            push!(edge_keys,minmax(block.nodes[i+1,column],block.nodes[j+1,column]))
        end
        for (i,j,k) in api._api_support_faces(family)
            push!(face_keys,api._api_support_face(block.nodes[i,column],
                block.nodes[j,column],block.nodes[k,column]))
        end
        for (i,j,k,l) in api._api_support_quads(family)
            push!(quad_keys,api._api_support_quad(block.nodes[i,column],block.nodes[j,column],
                block.nodes[k,column],block.nodes[l,column]))
        end
    end
    @test Set(keys(class.edge_entities)) ⊆ edge_keys
    @test Set(keys(class.face_entities)) ⊆ face_keys
    @test Set(keys(class.quad_entities)) ⊆ quad_keys
    @test count(owner->owner[1]==1,values(class.edge_entities))>0
    @test all(owner->owner[1]==2,values(class.face_entities))
    @test all(owner->owner[1]==2,values(class.quad_entities))
end

function _mrs_point_owner(api,target,expected)
    nodes=api.mesh.get_nodes()[1]
    selected=[node for node in nodes if Tuple(api.mesh.get_node(node)[1])==target]
    @test length(selected)==1
    isempty(selected) && return
    @test api.mesh.get_node(only(selected))[3:4]==expected
end

@testset "Mixed refinement inherits actual child support carriers" begin
    api=Tessella.API
    for family in (5,6),tagged in (false,true)
        try
            api.initialize()
            tagged ? _mrs_tagged_fixture(api,family) : _mrs_native_fixture(api,family)
            before=length(api.mesh.get_nodes()[1])
            api.mesh.remove_duplicate_nodes()
            @test length(api.mesh.get_nodes()[1])==before-1
            api.mesh.set_order(2)
            for round in 1:2
                api.mesh.refine()
                @test api.mesh.get_element_types(3)==Int32[family]
                _mrs_assert_actual_supports(api)
                api.mesh.set_order(2)
                @test api.mesh.get_element_types(3)==Int32[family==5 ? 12 : 13]
                step=1/2.0^(round+1)
                if !tagged && family==6
                    _mrs_point_owner(api,(1+step,0.,1.),(1,101))
                    _mrs_point_owner(api,(1.,step,1.),(1,103))
                    _mrs_point_owner(api,(1+step,0.,1+step),(2,101))
                    _mrs_point_owner(api,(1.,step,1+step),(2,103))
                elseif !tagged
                    for (target,surface) in (((1.,1+step,1+step),7),
                            ((1+step,1.,1+step),9),((1+step,1+step,1.),11))
                        _mrs_point_owner(api,target,(2,surface))
                    end
                end
                for (_,curve) in api.model.get_entities(1)
                    @test length(api.mesh.get_nodes(1,curve)[1])==2^(round+1)-1
                end
                for (_,surface) in api.model.get_entities(2)
                    corners=length(api.LAST_MESH_CLASS[].boundaries[(2,Int32(surface))])
                    expected=corners==4 ? (2^(round+1)-1)^2 :
                        (2^(round+1)-1)*(2^(round+1)-2)÷2
                    @test length(api.mesh.get_nodes(2,surface)[1])==expected
                end
                tagged && @test api.LAST_MESH_CLASS[].public_tags!==nothing
            end
        finally
            api.finalize()
        end
    end
    @testset "Imported body ties keep their existing supported contract" begin
        for dimension in (2,3)
            try
                api.initialize()
                api.model.add_discrete_entity(dimension,1)
                api.model.add_discrete_entity(dimension,2)
                if dimension==3
                    api.mesh.add_nodes(3,1,[1,2,3,4],[0.,0,0,1,0,0,0,1,0,0,0,1])
                    api.mesh.add_nodes(3,2,[5],[0.,0,-1])
                    api.mesh.add_elements_by_type(1,4,[1],[1,2,3,4])
                    api.mesh.add_elements_by_type(2,4,[2],[1,3,2,5])
                else
                    api.mesh.add_nodes(2,1,[1,2,3],[0.,0,0,1,0,0,0,1,0])
                    api.mesh.add_nodes(2,2,[4],[0.,-1,0])
                    api.mesh.add_elements_by_type(1,2,[1],[1,2,3])
                    api.mesh.add_elements_by_type(2,2,[2],[2,1,4])
                end
                api.mesh.generate(0)
                api.mesh.refine()
                @test sum(length,api.mesh.get_elements(dimension)[2])==(dimension==2 ? 8 : 16)
                @test length(api.mesh.get_nodes()[1])==(dimension==2 ? 9 : 14)
                @test all(node->api.mesh.get_node(node)[3] == dimension,api.mesh.get_nodes()[1])
            finally
                api.finalize()
            end
        end
    end
end
