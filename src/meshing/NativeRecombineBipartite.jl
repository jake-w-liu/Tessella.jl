# Public experimental quadrangulation has a separate transfer admission rule.
# Imported Face owner labels are unset. Labelling and topology only mutate
# vertex labels and preserve these unset Face labels, so the public exporter
# admits no Face for a positive CAD tag. Derive the surviving cells from that
# owner invariant; do not route this selector through a pairing algorithm.
const _PRIMARY_BIPARTITE_TYPE_ERROR =
    "GFace2PolyMesh only support quads (element type 3) and triangles (element type 2)"
const _PRIMARY_BIPARTITE_EMPTY_ERROR =
    "exception: access violation reading 0x0000000000000000"

function _primary_bipartite_face_events(cells,owner,caller)
    events=Tuple{Symbol,String}[
        (:trace,"$owner  Quadrangulation of face"),
        (:trace,"$owner  Transfer to Half Edge done")]
    for count in _primary_bipartite_boundary_lengths(cells,caller)
        isodd(count) || continue
        push!(events,(:warning,
            "Bipartite algorithm requires an even numer of boundary vertices on every edge loop that bounds face $owner. This loop has $count nodes and some triangles will remain."))
    end
    for message in ("Initial labelling done","internal labelling done",
            "optimum labelling done","steiner points added","quads created","2 gmsh done")
        push!(events,(:trace,"$owner  $message"))
    end
    return events
end

function _primary_bipartite_boundary_lengths(cells,caller)
    origin=Int32[];following=Int32[];opposite=Int32[]
    first_seen=Int32[];last_edge=Dict{Int32,Int32}()
    groups=Dict{Tuple{Int32,Int32},Vector{Int32}}()
    for cell in cells
        width=length(cell)
        start=Int32(length(origin)+1)
        for index in 1:width
            a=cell[index];b=cell[mod1(index+1,width)]
            length(origin)<typemax(Int32) || throw(ArgumentError(
                "$caller: experimental half-edge capacity exceeds Int32"))
            edge=Int32(length(origin)+1)
            haskey(last_edge,a) || push!(first_seen,a)
            last_edge[a]=edge
            push!(origin,a);push!(following,start+Int32(mod(index,width)))
            push!(opposite,0)
            push!(get!(()->Int32[],groups,minmax(a,b)),edge)
        end
    end
    for edges in values(groups)
        length(edges)<=2 || throw(ArgumentError(
            "$caller: experimental nonmanifold boundary traversal is not implemented"))
        length(edges)==2 || continue
        a,b=edges
        origin[a]==origin[following[b]] && origin[b]==origin[following[a]] ||
            throw(ArgumentError("$caller: experimental inconsistent half-edge orientation is not implemented"))
        opposite[a]=b;opposite[b]=a
    end
    visited=Set{Int32}();lengths=Int[]
    for vertex in first_seen
        vertex in visited && continue
        first=last_edge[vertex];edge=first
        boundary=false
        for _ in 1:length(origin)+1
            if opposite[edge]==0
                boundary=true;break
            end
            edge=following[opposite[edge]]
            edge==first && break
        end
        boundary || continue
        push!(visited,vertex)
        edge=first;count=0;closed=false
        for _ in 1:2length(origin)+1
            if opposite[edge]==0
                edge=following[edge];count+=1
                target=origin[edge]
                if target in visited
                    closed=true;break
                end
                push!(visited,target)
            else
                edge=following[opposite[edge]]
            end
        end
        closed || throw(ArgumentError("$caller: experimental boundary traversal failed to close"))
        push!(lengths,count)
    end
    return lengths
end

# Detached raw result. The caller publishes the completed prefix before a
# native Msg error throws; automatic numbering follows only a finished raw
# stage. Non-throwing Msg errors still allow later faces to be processed.
function _primary_bipartite_raw(mesh,face_tags,caller;abort=2)
    data=mesh.entity_data
    consumed=Set{Int32}();events=Tuple{Symbol,String}[]
    finished=true
    groups=Dict{Int32,Vector{Tuple{Int,Int}}}()
    for (bi,block) in enumerate(mesh.blocks)
        Elements.msh_dimension(block.msh)==2 || continue
        for column in axes(block.nodes,2)
            owner=data.block_entities[bi][column]
            push!(get!(()->Tuple{Int,Int}[],groups,owner),(bi,column))
        end
    end
    for owner in sort!(unique!(Int32.(face_tags)))
        selected=get(groups,owner,Tuple{Int,Int}[])
        if isempty(selected)
            push!(events,(:trace,"$owner  Quadrangulation of face"))
            push!(events,(:fault,_PRIMARY_BIPARTITE_EMPTY_ERROR))
            finished=false;break
        end
        if any(entry->mesh.blocks[first(entry)].msh ∉ (2,3),selected)
            push!(events,(:trace,"$owner  Quadrangulation of face"))
            push!(events,(:error,_PRIMARY_BIPARTITE_TYPE_ERROR))
            if abort in (2,3,4)
                finished=false;break
            end
            push!(events,(:trace,"$owner  Transfer to Half Edge done"))
            continue
        end
        cells=[collect(@view(mesh.blocks[bi].nodes[:,column])) for (bi,column) in selected]
        append!(events,_primary_bipartite_face_events(cells,owner,caller))
        push!(consumed,owner)
    end
    keep_columns=Vector{Int}[]
    for (bi,block) in enumerate(mesh.blocks)
        push!(keep_columns,[column for column in axes(block.nodes,2) if
            !(Elements.msh_dimension(block.msh)==2 && data.block_entities[bi][column] in consumed)])
    end
    return (;keep_columns,events,finished,consumed)
end
