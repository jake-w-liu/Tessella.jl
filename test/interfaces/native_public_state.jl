module NativePublicState
using Tessella,Test
const A=Tessella.API
bits(v)=string.(reinterpret(UInt64,Float64.(v));base=16,pad=16)
cell_rows(types,tags,nodes)=[Dict("type"=>Int(t),"tags"=>string.(labels),"nodes"=>string.(conn))
    for (t,labels,conn) in zip(types,tags,nodes)]
function public_state(;include_scalar_cells=true)
    nt,xyz,uv=A.mesh.get_nodes();types,tags,nodes=A.mesh.get_elements()
    entities=Any[]
    for (d,t) in A.model.get_entities()
        en,co,pa=A.mesh.get_nodes(d,t,false,true)
        push!(entities,Dict("owner"=>[Int(d),Int(t)],"kind"=>A.model.get_type(d,t),
            "boundary"=>[collect(p) for p in A.model.get_boundary([(d,t)],false,true,false)],
            "physical"=>Int.(A.model.get_physical_groups_for_entity(d,t)),
            "parent"=>collect(A.model.get_parent(d,t)),"partitions"=>Int.(A.model.get_partitions(d,t)),
            "nodes"=>string.(en),"xyz_bits"=>bits(co),"uv_bits"=>bits(pa),
            "cells"=>cell_rows(A.mesh.get_elements(d,t)...)))
    end
    node_details=Any[]
    for uid in nt
        co,pa,d,t=A.mesh.get_node(uid)
        push!(node_details,Dict("uid"=>string(uid),"owner"=>[Int(d),Int(t)],
            "xyz_bits"=>bits(co),"uv_bits"=>bits(pa)))
    end
    cell_details=Any[]
    if include_scalar_cells
        for bucket in tags,uid in bucket
            kind,connected,d,t=A.mesh.get_element(uid)
            push!(cell_details,Dict("uid"=>string(uid),"type"=>Int(kind),"nodes"=>string.(connected),
                "owner"=>[Int(d),Int(t)]))
        end
    end
    groups=[Dict("dim"=>Int(d),"tag"=>Int(t),"name"=>A.model.get_physical_name(d,t),
        "members"=>Int.(A.model.get_entities_for_physical_group(d,t))) for (d,t) in A.model.get_physical_groups()]
    result=Dict("nodes"=>string.(nt),"xyz_bits"=>bits(xyz),"uv_bits"=>bits(uv),
        "cells"=>cell_rows(types,tags,nodes),"entities"=>entities,"node_details"=>node_details,
        "groups"=>groups,"max_node"=>string(A.mesh.get_max_node_tag()),
        "max_element"=>string(A.mesh.get_max_element_tag()))
    include_scalar_cells && (result["cell_details"]=cell_details)
    return result
end
function assert_state(expected)
    actual=public_state()
    for key in keys(expected)
        @test actual[key]==expected[key]
    end
end
end
