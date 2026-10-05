using Tessella.Elements:msh_dimension,msh_spec,msh_type

# Match one independently certified volume to its classified boundary. This
# correspondence preserves public node identity; it is never used to weld
# coordinates across independent, coincident entities.
function free_actual_carriers(carriers;source=nothing,volume=nothing,projected=nothing,public_tags=nothing)
    model=API.CURRENT[]
    source===nothing && (source=Tessella.Model.mesh_model_surface(model,carriers.source))
    volume===nothing && (volume=Tessella.Model.mesh_model_volume(model,carriers.volume))
    projected===nothing && (projected=Tessella.Model.model_to_mixed(model,volume,3,carriers.volume))
    tags,xyz,_=API.mesh.get_nodes(3,carriers.volume,true)
    remap=if public_tags===nothing
        positions=Dict(Tuple(p)=>tag for (tag,p) in zip(tags,eachcol(reshape(xyz,3,:))))
        @test length(positions)==length(tags)
        [positions[Tuple(projected.coords[:,i])] for i in axes(projected.coords,2)]
    else
        # The public setters receive this exact projection-index -> external-ID
        # map. Verify it directly; coordinates never determine seeded identity.
        @test length(public_tags)==size(projected.coords,2) && allunique(public_tags)
        @test Set(public_tags)==Set(tags)
        for i in eachindex(public_tags)
            @test API.mesh.get_node(public_tags[i])[1]==projected.coords[:,i]
        end
        public_tags
    end
    edges=Dict{Tuple,Tuple{Int,Int}}()
    faces=Dict{Tuple,Tuple{Int,Int}}()
    quads=Dict{Tuple,Tuple{Int,Int}}()
    closures=Dict{Tuple{Int,Int},Set{UInt64}}()
    owners=Dict(remap[i]=>(Int(owner[1]),Int(owner[2]))
        for (i,owner) in enumerate(projected.entity_data.node_entities))
    function insert!(table,support,owner)
        old=get(table,support,nothing)
        if old===nothing || owner[1]<old[1]
            table[support]=owner
        elseif owner[1]==old[1]
            @test old==owner
        end
    end
    for (index,block) in enumerate(projected.blocks)
        dim=msh_dimension(block.msh);dim<=2 || continue
        for (column,cell) in enumerate(eachcol(block.nodes))
            owner=(dim,Int(projected.elementary_entities[index][column]))
            public=Tuple(remap[node] for node in cell)
            union!(get!(closures,owner,Set{UInt64}()),public)
            if dim==1
                insert!(edges,Cert.key(public),owner)
            elseif dim==2
                insert!(faces,Cert.key(public),owner)
                for row in eachindex(public)
                    insert!(edges,Cert.key((public[row],public[mod1(row+1,length(public))])),owner)
                end
                length(public)==4 && insert!(quads,Cert.key(public),owner)
            end
        end
    end
    @test Set(keys(owners))==Set(tags)
    for (node,owner) in owners
        @test API.mesh.get_node(node)[3:4]==owner
    end
    return (;edges,faces,quads,closures,owners,source,projected,
        primary=Set(tags),volume=carriers.volume)
end

function free_public_support_nodes(volume)
    types,_,connectivity=API.mesh.get_elements(3,volume)
    result=Dict{UInt64,Tuple}()
    function insert(node,support)
        @test !haskey(result,node) || result[node]==support
        result[node]=support
    end
    for (msh,nodes) in zip(types,connectivity)
        for edge in eachcol(reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:))
            insert(edge[3],Cert.key(edge[1:2]))
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(API.mesh.get_element_face_nodes(msh,4,volume,false),9,:))
                insert(face[9],Cert.key(face[1:4]))
            end
        end
        if msh==12
            for cell in eachcol(reshape(nodes,27,:))
                insert(cell[27],Cert.key(cell[1:8]))
            end
        end
    end
    @test length(Set(values(result)))==length(result)
    return result
end

free_support_carrier(support,expected)=length(support)==2 ?
    get(expected.edges,support,(3,expected.volume)) : length(support)==4 ?
    get(expected.quads,support,(3,expected.volume)) : (3,expected.volume)

function free_check_parameters(entity,closure,owned)
    dim,tag=entity
    own,ownxyz,ownuv=API.mesh.get_nodes(dim,tag,false,true)
    alltags,xyz,uv=API.mesh.get_nodes(dim,tag,true,true)
    @test Set(own)==owned && Set(alltags)==closure
    @test allunique(own) && allunique(alltags)
    @test length(ownxyz)==3length(owned) && length(ownuv)==dim*length(owned)
    @test length(xyz)==3length(closure) && length(uv)==dim*length(closure)
    @test isempty(xyz) || maximum(abs.(API.model.get_value(dim,tag,uv).-xyz))<=2e-11
    for node in alltags
        point,parameters,owner_dim,owner=API.mesh.get_node(node)
        @test length(parameters)==owner_dim
        @test maximum(abs.(API.model.get_value(owner_dim,owner,parameters).-point))<=2e-11
    end
end

function free_check_typed_faces(volume,expected)
    occurrences=Dict{Tuple,Vector{Tuple}}()
    for msh in API.mesh.get_element_types(3,volume)
        family=Dict(11=>4,12=>5,13=>6,14=>7)[Int(msh)]
        for width in (3,4)
            any(length(pattern)==width for pattern in Cert.FACES[family]) || continue
            nodes=API.mesh.get_element_face_nodes(msh,width,volume,false)
            for face in eachcol(reshape(nodes,width==3 ? 6 : 9,:))
                support=Cert.key(face[1:width])
                push!(get!(occurrences,support,Tuple[]),Cert.key(face))
            end
        end
    end
    @test all(length(rows) in (1,2) for rows in values(occurrences))
    @test all(first(rows)==last(rows) for rows in values(occurrences) if length(rows)==2)
    boundary=Set(support for (support,rows) in occurrences if length(rows)==1)
    @test boundary==Set(keys(expected.faces))
    return occurrences
end

function free_check_actual_p2(carriers,f,expected;certify_geometry=true)
    supports=free_public_support_nodes(carriers.volume)
    alltags=Set(API.mesh.get_nodes(3,carriers.volume,true)[1])
    @test Set(keys(supports))==setdiff(alltags,expected.primary)
    @test isempty(intersect(Set(keys(supports)),expected.primary))
    derived=Dict(expected.owners)
    for (node,support) in supports
        point,_,dim,owner=API.mesh.get_node(node)
        carrier=free_support_carrier(support,expected)
        @test (dim,owner)==carrier
        derived[node]=carrier
        corners=[API.mesh.get_node(primary)[1] for primary in support]
        mean=[sum(c[d] for c in corners)/length(corners) for d in 1:3]
        @test maximum(abs.(point.-mean))<=2e-11
    end
    for (entity,primary) in expected.closures
        entity[1] in (1,2) || continue
        closure=union(primary,Set(node for (node,support) in supports if all(in(primary),support)))
        owned=Set(node for (node,owner) in derived if owner==entity)
        free_check_parameters(entity,closure,owned)
    end
    faces=free_check_typed_faces(carriers.volume,expected)
    certificate=certify_geometry ? Cert.certify_quadratic(Cert.public_volume(API,carriers.volume),f,expected.source) : nothing
    return (;supports,certificate,faces)
end
