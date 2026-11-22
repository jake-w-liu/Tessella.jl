# The entity that owns each order-two support node, by Gmsh 4.15.2's rule:
# SetOrderN creates an edge or face node with the lowest-dimensional element
# that has that edge or face — a Curve line, then a Surface cell — and the
# Volume owns the rest. Read from the public P1 mesh of the Volume's closure
# (getNodes(dim,tag) per entity in Gmsh agrees node for node with Tessella on
# these fixtures, checked at orders 1 and 2).
function gmsh_carriers(api,volume)
    key(nodes)=Tuple(sort!(Int.(collect(nodes))))
    edges=Dict{Tuple,Tuple{Int,Int}}();quads=Dict{Tuple,Tuple{Int,Int}}()
    closure=Set{Tuple{Int,Int}}()
    pending=[(3,volume)]
    while !isempty(pending)
        entity=pop!(pending)
        entity in closure && continue
        push!(closure,entity)
        entity[1]>1 && append!(pending,[(d,abs(t)) for (d,t) in
            api.model.get_boundary([entity],false,false,false)])
    end
    for dim in (1,2),(d,tag) in sort!(collect(closure))
        d==dim || continue
        types,_,connectivity=api.mesh.get_elements(dim,tag)
        for (msh,nodes) in zip(types,connectivity)
            width=Tessella.Elements.msh_spec(msh).nnodes
            for cell in eachcol(reshape(nodes,width,:))
                if dim==1
                    get!(edges,key(cell[1:2]),(1,tag))
                else
                    corners=msh==2 ? 3 : 4
                    for i in 1:corners
                        get!(edges,key((cell[i],cell[mod1(i+1,corners)])),(2,tag))
                    end
                    corners==4 && get!(quads,key(cell[1:4]),(2,tag))
                end
            end
        end
    end
    return edges,quads
end
