using Tessella.Elements:msh_dimension

function b4_actual_carriers(carriers)
    model=API.CURRENT[]
    source=Tessella.Model.mesh_model_surface(model,carriers.source)
    volume=Tessella.Model.mesh_model_volume(model,carriers.volume)
    projected=Tessella.Model.model_to_mixed(model,volume,3,carriers.volume)
    tags,xyz,_=API.mesh.get_nodes(3,carriers.volume,true)
    points=reshape(xyz,3,:)
    # The independent product certificate establishes distinct coordinates
    # within this one volume. This is only a query correspondence; it does
    # not identify coincident nodes in unrelated regions or change a cache.
    positions=Dict(Tuple(p)=>tag for (tag,p) in zip(tags,eachcol(points)))
    @test length(positions)==length(tags)
    remap=[positions[Tuple(projected.coords[:,i])] for i in axes(projected.coords,2)]
    edges=Dict{Tuple,Tuple{Int,Int}}();quads=Dict{Tuple,Tuple{Int,Int}}()
    function insert!(table,support,owner)
        old=get(table,support,nothing)
        if old===nothing || owner[1]<old[1]
            table[support]=owner
        elseif owner[1]==old[1]
            @test old==owner
        end
    end
    for (index,block) in enumerate(projected.blocks)
        dim=msh_dimension(block.msh);dim in (1,2) || continue
        owners=projected.elementary_entities[index]
        for (column,cell) in enumerate(eachcol(block.nodes))
            owner=(dim,Int(owners[column]))
            if block.msh==1
                insert!(edges,Cert.key(remap[v] for v in cell),owner)
            else
                for row in eachindex(cell)
                    insert!(edges,Cert.key((remap[cell[row]],remap[cell[mod1(row+1,length(cell))]])),owner)
                end
                block.msh==3 && insert!(quads,Cert.key(remap[v] for v in cell),owner)
            end
        end
    end
    primary_keys=Dict(node=>(Tuple(point),dim,owner) for node in tags
        for (point,_,dim,owner) in (API.mesh.get_node(node),))
    return (;edges,quads,source,primary=Set(tags),primary_keys,volume=carriers.volume)
end

function b4_public_support_nodes(volume)
    types,tags,connectivity=API.mesh.get_elements(3,volume)
    result=Dict{UInt64,Tuple}()
    function insert(node,support)
        @test !haskey(result,node) || result[node]==support
        result[node]=support
    end
    for (msh,ids,nodes) in zip(types,tags,connectivity)
        for edge in eachcol(reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:))
            insert(edge[3],Cert.key(edge[1:2]))
        end
        # Pri18 includes Quad9 face centers.
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
    return result
end

function b4_check_parameters(dim,entity,boundary,expected)
    tags,xyz,uv=API.mesh.get_nodes(dim,entity,boundary,true)
    @test length(tags)==expected && allunique(tags)
    @test length(xyz)==3expected && length(uv)==dim*expected
    @test maximum(abs.(API.model.get_value(dim,entity,uv).-xyz))<=2e-11
    for node in tags
        p,parameters,owner_dim,owner=API.mesh.get_node(node)
        @test length(parameters)==owner_dim
        @test maximum(abs.(API.model.get_value(owner_dim,owner,parameters).-p))<=2e-11
    end
    return Set(tags)
end

function b4_check_actual_p2(carriers,f,expected)
    m,n=f.strip_length,f.intervals
    check_owner_totals(carriers.volume,m,n,2)
    supports=b4_public_support_nodes(carriers.volume)
    @test length(supports)==expected_quadratic(m,n)-expected_primary(m,n)
    # Order elevation can renumber primaries; retain their exact position and owner binding.
    primary=Set{UInt64}()
    types,_,connectivity=API.mesh.get_elements(3,carriers.volume)
    for (msh,nodes) in zip(types,connectivity)
        _,_,_,width,_,corners=API.mesh.get_element_properties(msh)
        for cell in eachcol(reshape(nodes,width,:)),node in cell[1:corners]
            push!(primary,node)
        end
    end
    current_keys=Dict((Tuple(point),dim,owner)=>node for node in primary
        for (point,_,dim,owner) in (API.mesh.get_node(node),))
    @test length(current_keys)==length(primary)==length(expected.primary)
    binding=Dict(node=>current_keys[key] for (node,key) in expected.primary_keys)
    edges=Dict(Cert.key(binding[node] for node in support)=>owner for (support,owner) in expected.edges)
    quads=Dict(Cert.key(binding[node] for node in support)=>owner for (support,owner) in expected.quads)
    @test Set(keys(supports))==setdiff(Set(API.mesh.get_nodes(3,carriers.volume,true)[1]),Set(values(binding)))
    for (node,support) in supports
        p,uv,dim,owner=API.mesh.get_node(node)
        carrier=length(support)==2 ? get(edges,support,(3,carriers.volume)) :
            length(support)==4 ? get(quads,support,(3,carriers.volume)) : (3,carriers.volume)
        @test (dim,owner)==carrier
        corners=[API.mesh.get_node(v)[1] for v in support]
        mean=[sum(c[d] for c in corners)/length(corners) for d in 1:3]
        @test maximum(abs.(p.-mean))<=2e-11
    end
    for surface in (carriers.source,carriers.top)
        own,closure=expected_cap(m)
        @test b4_check_parameters(2,surface,false,own)⊆b4_check_parameters(2,surface,true,closure)
    end
    for surface in carriers.laterals
        link=API.CURRENT[].meshing.extrude_sources[(2,surface)]
        curve=abs(link[2]);slot=findfirst(==(curve),f.curve_tags)
        @test slot!==nothing
        own,closure=slot in f.long_pair ? expected_long(m,n) : expected_short(n)
        @test b4_check_parameters(2,surface,false,own)⊆b4_check_parameters(2,surface,true,closure)
    end
    for (slot,curve) in pairs(f.curve_tags)
        segments=slot in f.long_pair ? m : 1
        b4_check_parameters(1,curve,false,2segments-1)
        b4_check_parameters(1,curve,true,2segments+1)
    end
    actual=Cert.public_volume(API,carriers.volume)
    certificate=Cert.certify_quadratic(actual,f,expected.source)
    @test certificate.ncell==m*(n+6) && length(certificate.centers)==m
    return (;supports,certificate)
end

function b4_native_construction_cases()
    for (m,profile,direction) in ((5,:one,1),(7,:nonbinary,-1),(9,:one,1))
        f=Cert.fixture("actual_p2";coherent=true,strip_length=m,layers=profile,height=direction,pins=(1,2,3,4))
        with_model(f,carriers->begin
            API.mesh.generate(3)
            check_owner_totals(carriers.volume,m,f.intervals,1)
            @test all(!isempty(API.mesh.get_elements(dim,tag)[1])
                for dim in 0:2 for (_,tag) in API.model.get_entities(dim))
            expected=b4_actual_carriers(carriers)
            p1=Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source)
            @test p1.total==p1.source.area*abs(Cert.Q(last(p1.heights))-Cert.Q(first(p1.heights)))
            API.mesh.set_order(2)
            b4_check_actual_p2(carriers,f,expected)
            actual_before=(API.mesh.get_nodes(),API.mesh.get_elements())
            API.mesh.set_order(2)
            @test actual_before==(API.mesh.get_nodes(),API.mesh.get_elements())
        end)
    end
end
