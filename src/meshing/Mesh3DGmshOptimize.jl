# Native tetrahedron quality optimization. Storage order and element births
# belong to this operation; the caller retains its public node identities.

const _GOPT_EDGES=((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
const _GOPT_EDGE_FACES=((1,3),(1,2),(2,3),(1,4),(3,4),(2,4))
const _GOPT_FACES=((1,2,3),(1,3,4),(1,2,4),(2,3,4))
const _GOPT_VERTEX_FACES=((1,2,3),(1,3,4),(1,2,4),(2,3,4))
const _GOPT_EDGE_PERMUTATION=(1,4,2,3,6,5)


# Literal native ring3–7 pattern order; pinned CPP SHA-256 abfb4c8cbda1356aeaf8cf017eea2ebf103acbbd2a22664a49664c0157bec658.
const _GOPT_SWAP_PATTERNS=Dict{Int,Tuple{Tuple,Tuple}}(
    3=>(((1,2,3),),((1,),)),
    4=>(((1,2,3),(3,4,1),(2,3,4),(4,1,2)),((1,2),(3,4))),
    5=>(((1,2,3),(1,3,4),(1,4,5),(1,2,5),(2,4,5),(2,3,4),(3,4,5),(1,3,5),(1,2,4),(2,3,5)),((1,2,3),(4,5,6),(1,7,8),(3,6,9),(4,7,10))),
    6=>(((1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,3,6),(3,5,6),(3,4,5),(1,4,6),(4,5,6),(1,3,5),(3,4,6),(2,3,4),(1,2,4),(1,2,6),(2,5,6),(2,4,5),(1,2,5),(2,4,6),(2,3,5),(2,3,6)),((1,2,3,4),(1,5,6,7),(1,2,8,9),(1,4,7,10),(1,5,9,11),(3,4,12,13),(12,14,15,16),(8,9,12,13),(4,12,16,17),(9,12,14,18),(7,14,15,19),(4,7,17,19),(6,7,14,20),(9,11,14,20))),
    7=>(((1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,6,7),(1,4,7),(4,6,7),(4,5,6),(1,5,7),(5,6,7),(1,4,6),(4,5,7),(1,3,5),(3,4,5),(1,3,7),(3,6,7),(3,5,6),(1,3,6),(3,5,7),(3,4,6),(3,4,7),(1,2,4),(2,3,4),(1,2,5),(2,4,5),(1,2,7),(2,6,7),(2,5,6),(1,2,6),(2,5,7),(2,4,6),(2,4,7),(2,3,5),(2,3,6),(2,3,7)),((1,2,3,4,5),(1,2,6,7,8),(1,2,3,9,10),(1,2,5,8,11),(1,2,6,10,12),(1,4,5,13,14),(1,14,15,16,17),(1,9,10,13,14),(1,5,14,17,18),(1,10,14,15,19),(1,8,15,16,20),(1,5,8,18,20),(1,7,8,15,21),(1,10,12,15,21),(3,4,5,22,23),(6,7,8,22,23),(3,9,10,22,23),(5,8,11,22,23),(6,10,12,22,23),(4,5,23,24,25),(23,25,26,27,28),(9,10,23,24,25),(5,23,25,28,29),(10,23,25,26,30),(8,23,26,27,31),(5,8,23,29,31),(7,8,23,26,32),(10,12,23,26,32),(4,5,14,24,33),(14,26,27,28,33),(9,10,14,24,33),(5,14,28,29,33),(10,14,26,30,33),(14,17,26,27,34),(5,14,17,29,34),(14,16,17,26,35),(10,14,19,26,35),(8,20,26,27,34),(5,8,20,29,34),(8,16,20,26,35),(7,8,21,26,35),(10,12,21,26,35))),
)

mutable struct _GmshOptimizeTet
    vertices::NTuple{4,Int32}
    neighbors::Vector{Int32}
    number::UInt64
    physical::Int32
    quality::Float64
    deleted::Bool
end

mutable struct _GmshOptimizeState
    coords::Matrix{Float64}
    cells::Vector{_GmshOptimizeTet}
    element_max::UInt64
    embedded_faces::Set{NTuple{3,Int32}}
    embedded_edges::Set{Tuple{Int32,Int32}}
    movable::BitVector
    # Live-cell carrier counts for each embedded entity. Embedded faces are
    # adjacency cuts (`_gopt_connect!`), but a swap can still delete every
    # carrier of a face or edge the ring walk never crosses; the counts let a
    # swap detect that case and only accept a reconnection that reproduces it.
    embedded_face_carriers::Dict{NTuple{3,Int32},Int}
    embedded_edge_carriers::Dict{Tuple{Int32,Int32},Int}
    # Vertices of embedded faces. A reconnection that creates a non-embedded
    # face entirely from these vertices adds a second covering layer over the
    # embedded sheet, so such triangulations are rejected outright.
    embedded_verts::Set{Int32}
end

# `delta` is +1 on cell creation and -1 on deletion.
function _gopt_note_carriers!(state,vertices,delta)
    if !isempty(state.embedded_face_carriers)
        for face in 1:4
            key=_gopt_face(vertices,face)
            haskey(state.embedded_face_carriers,key) &&
                (state.embedded_face_carriers[key]+=delta)
        end
    end
    if !isempty(state.embedded_edge_carriers)
        for e in _GOPT_EDGES
            key=_gopt_edge(vertices[e[1]],vertices[e[2]])
            haskey(state.embedded_edge_carriers,key) &&
                (state.embedded_edge_carriers[key]+=delta)
        end
    end
    return nothing
end

@inline _gopt_point(state,row)=
    (state.coords[1,row],state.coords[2,row],state.coords[3,row])
@inline _gopt_subtract(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _gopt_dot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _gopt_cross(a,b)=(a[2]*b[3]-a[3]*b[2],
    -a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
@inline _gopt_norm(a)=sqrt(_gopt_dot(a,a))
@inline _gopt_area(a,b,c)=0.5*_gopt_norm(
    _gopt_cross(_gopt_subtract(c,b),_gopt_subtract(a,b)))

function _gopt_gamma(a,b,c,d)
    volume=_orient3_value(a,b,c,d)/6.0
    abs(volume)==0.0 && return 0.0,abs(volume)
    ab=_gopt_subtract(b,a);ac=_gopt_subtract(c,a);ad=_gopt_subtract(d,a)
    dc=_gopt_subtract(d,c);db=_gopt_subtract(d,b);cb=_gopt_subtract(c,b)
    la=_gopt_dot(ab,ab);lb=_gopt_dot(ac,ac);lc=_gopt_dot(ad,ad)
    lA=_gopt_dot(dc,dc);lB=_gopt_dot(db,db);lC=_gopt_dot(cb,cb)
    lalA=sqrt(la*lA);lblB=sqrt(lb*lB);lclC=sqrt(lc*lC)
    inside=(lalA+lblB+lclC)*(lalA+lblB-lclC)*
        (lalA-lblB+lclC)*(-lalA+lblB+lclC)
    inside<=0.0 && return 0.0,abs(volume)
    partR=sqrt(inside)*(1.0/24.0)
    s1=abs(_gopt_area(a,b,c));s2=abs(_gopt_area(a,c,d))
    s3=abs(_gopt_area(a,b,d));s4=abs(_gopt_area(b,c,d))
    rho=3*3*volume/(s1+s2+s3+s4)
    return rho*volume/partR,abs(volume)
end

function _gopt_quality(state,vertices)
    a=_gopt_point(state,vertices[1]);b=_gopt_point(state,vertices[2])
    c=_gopt_point(state,vertices[3]);d=_gopt_point(state,vertices[4])
    return _gopt_gamma(a,b,c,d)
end

function _gopt_volume(state,vertices)
    a=_gopt_point(state,vertices[1]);b=_gopt_point(state,vertices[2])
    c=_gopt_point(state,vertices[3]);d=_gopt_point(state,vertices[4])
    e1=_gopt_subtract(b,a);e2=_gopt_subtract(c,a);e3=_gopt_subtract(d,a)
    determinant=e1[1]*(e2[2]*e3[3]-e3[2]*e2[3])-
        e2[1]*(e1[2]*e3[3]-e3[2]*e1[3])+
        e3[1]*(e1[2]*e2[3]-e2[2]*e1[3])
    return determinant*0.166666666666666666
end

@inline _gopt_face(vertices,face)=_sorted_face3(
    vertices[_GOPT_FACES[face][1]],vertices[_GOPT_FACES[face][2]],
    vertices[_GOPT_FACES[face][3]])
@inline _gopt_edge(first,last)=first<last ? (first,last) : (last,first)

# Reconnection retains neighbors outside the supplied local cavity. Embedded
# faces remain cuts in the graph, as in native connectTets.
function _gopt_connect!(state,indices)
    faces=Dict{NTuple{3,Int32},Tuple{Int32,Int}}()
    for index in indices
        cell=state.cells[index]
        cell.deleted && continue
        for face in 1:4
            key=_gopt_face(cell.vertices,face)
            key in state.embedded_faces && continue
            previous=get(faces,key,nothing)
            if previous===nothing
                faces[key]=(Int32(index),face)
            elseif previous[1]!=index
                state.cells[previous[1]].neighbors[previous[2]]=Int32(index)
                cell.neighbors[face]=previous[1]
            end
        end
    end
    return nothing
end

function _gopt_edge_cavity(state,start,edge)
    cell=state.cells[start]
    local_edge=_GOPT_EDGE_PERMUTATION[edge]
    first=cell.vertices[_GOPT_EDGES[local_edge][1]]
    last=cell.vertices[_GOPT_EDGES[local_edge][2]]
    last_in_ring=cell.vertices[_GOPT_EDGES[7-local_edge][1]]
    cavity=Int32[Int32(start)];ring=Int32[last_in_ring]
    current=Int32(start)
    while true
        vertices=state.cells[current].vertices
        opposite=_GOPT_EDGES[7-local_edge]
        first_opposite=vertices[opposite[1]]
        last_opposite=vertices[opposite[2]]
        role=first_opposite==last_in_ring ? 2 : 1
        last_in_ring=first_opposite==last_in_ring ? last_opposite : first_opposite
        target_slot=opposite[role]
        candidate_faces=_GOPT_EDGE_FACES[local_edge]
        face=target_slot in _GOPT_FACES[candidate_faces[1]] ? candidate_faces[1] :
            target_slot in _GOPT_FACES[candidate_faces[2]] ? candidate_faces[2] : 0
        face==0 && return nothing
        current=state.cells[current].neighbors[face]
        current==0 && return nothing
        state.cells[current].deleted && return nothing
        current==start && break
        push!(ring,last_in_ring);push!(cavity,current)
        local_edge=0
        vertices=state.cells[current].vertices
        for index in 1:6
            a=vertices[_GOPT_EDGES[index][1]];b=vertices[_GOPT_EDGES[index][2]]
            if (a==first && b==last) || (a==last && b==first)
                local_edge=index
                break
            end
        end
        local_edge==0 && return nothing
        length(cavity)>1000 && return nothing
    end
    outside=Int32[]
    for index in cavity, neighbor in state.cells[index].neighbors
        neighbor==0 || neighbor in outside || neighbor in cavity || push!(outside,neighbor)
    end
    return first,last,cavity,outside,ring
end

function _gopt_new_cell!(state,vertices,quality,physical)
    state.element_max=Base.checked_add(state.element_max,UInt64(1))
    length(state.cells)<typemax(Int32) || throw(ArgumentError(
        "native tetrahedron optimization exceeds Int32 storage capacity"))
    push!(state.cells,_GmshOptimizeTet(vertices,zeros(Int32,4),
        state.element_max,physical,quality,false))
    _gopt_note_carriers!(state,vertices,1)
    return Int32(length(state.cells))
end

function _gopt_swap_edge!(state,new_cells,start,edge)
    found=_gopt_edge_cavity(state,start,edge)
    found===nothing && return false
    first,last,cavity,outside,ring=found
    pattern=get(_GOPT_SWAP_PATTERNS,length(ring),nothing)
    pattern===nothing && return false
    # Embedded entities are hard constraints. Faces/edges carried by a cell
    # outside the cavity survive regardless; an entity whose every carrier is
    # deleted must be reproduced by the chosen reconnection.
    constrained_faces=Dict{NTuple{3,Int32},Int}()
    constrained_edges=Dict{Tuple{Int32,Int32},Int}()
    for index in cavity
        vertices=state.cells[index].vertices
        for face in 1:4
            key=_gopt_face(vertices,face)
            haskey(state.embedded_face_carriers,key) &&
                (constrained_faces[key]=get(constrained_faces,key,0)+1)
        end
        for e in _GOPT_EDGES
            key=_gopt_edge(vertices[e[1]],vertices[e[2]])
            haskey(state.embedded_edge_carriers,key) &&
                (constrained_edges[key]=get(constrained_edges,key,0)+1)
        end
    end
    filter!(pair->state.embedded_face_carriers[pair[1]]==pair[2],
            constrained_faces)
    filter!(pair->state.embedded_edge_carriers[pair[1]]==pair[2],
            constrained_edges)
    volume_ref=0.0;quality_ref=1.0
    for index in cavity
        cell=state.cells[index]
        volume_ref+=abs(_gopt_volume(state,cell.vertices))
        quality_ref=min(quality_ref,cell.quality)
    end
    triangles,triangulations=pattern
    count=length(triangles)
    qualities_first=Vector{Float64}(undef,count)
    qualities_last=similar(qualities_first)
    volumes_first=similar(qualities_first);volumes_last=similar(qualities_first)
    for index in eachindex(triangles)
        slots=triangles[index]
        a=ring[slots[1]];b=ring[slots[2]];c=ring[slots[3]]
        qualities_first[index],volumes_first[index]=_gopt_quality(state,(a,b,c,first))
        qualities_last[index],volumes_last[index]=_gopt_quality(state,(a,b,c,last))
    end
    best=-1.0;best_index=1
    for (index,triangulation) in enumerate(triangulations)
        quality=1.0;volume=0.0
        for triangle in triangulation
            quality=min(quality,qualities_first[triangle])
            quality=min(quality,qualities_last[triangle])
            volume+=volumes_first[triangle]+volumes_last[triangle]
        end
        abs(volume-volume_ref)>1e-10*(volume+volume_ref) && (quality=-1.0)
        if quality!=-1.0 &&
                !(isempty(constrained_faces) && isempty(constrained_edges) &&
                  isempty(state.embedded_verts))
            uncovered_faces=Set(keys(constrained_faces))
            uncovered_edges=Set(keys(constrained_edges))
            violation=false
            for triangle in triangulation
                slots=triangles[triangle]
                a=ring[slots[1]];b=ring[slots[2]];c=ring[slots[3]]
                for vertices in ((a,b,c,first),(c,b,a,last))
                    for face in 1:4
                        key=_gopt_face(vertices,face)
                        delete!(uncovered_faces,key)
                        if !violation &&
                                !haskey(state.embedded_face_carriers,key) &&
                                key[1] in state.embedded_verts &&
                                key[2] in state.embedded_verts &&
                                key[3] in state.embedded_verts
                            violation=true
                        end
                    end
                    for e in _GOPT_EDGES
                        delete!(uncovered_edges,
                            _gopt_edge(vertices[e[1]],vertices[e[2]]))
                    end
                end
            end
            isempty(uncovered_faces) && isempty(uncovered_edges) &&
                !violation || (quality=-1.0)
        end
        if quality>best
            best=quality;best_index=index
        end
    end
    best<=quality_ref+1e-20 && return false
    for triangle in triangulations[best_index]
        slots=triangles[triangle]
        a=ring[slots[1]];b=ring[slots[2]];c=ring[slots[3]]
        physical=state.cells[cavity[1]].physical
        first_cell=_gopt_new_cell!(state,(a,b,c,first),qualities_first[triangle],physical)
        last_cell=_gopt_new_cell!(state,(c,b,a,last),qualities_last[triangle],physical)
        push!(outside,first_cell);push!(outside,last_cell)
        push!(new_cells,first_cell);push!(new_cells,last_cell)
    end
    for index in cavity
        cell=state.cells[index]
        cell.deleted=true
        _gopt_note_carriers!(state,cell.vertices,-1)
    end
    _gopt_connect!(state,outside)
    return true
end

# Use a stack of continuation frames to retain native depth-first visitation
# order without imposing the native recursion's process-stack limit.
function _gopt_vertex_cavity(state,start,vertex)
    cavity=Int32[Int32(start)];seen=Set(cavity)
    frames=Tuple{Int32,Int}[(Int32(start),1)]
    while !isempty(frames)
        index,position=pop!(frames)
        cell=state.cells[index]
        cell.deleted && return nothing
        role=findfirst(==(vertex),cell.vertices)
        role===nothing && return nothing
        position>3 && continue
        push!(frames,(index,position+1))
        neighbor=cell.neighbors[_GOPT_VERTEX_FACES[role][position]]
        neighbor==0 || neighbor in seen || begin
            push!(seen,neighbor);push!(cavity,neighbor);push!(frames,(neighbor,1))
        end
    end
    return cavity
end

function _gopt_smooth_vertex!(state,index,role)
    cell=state.cells[index]
    cell.deleted && return false
    vertex=cell.vertices[role]
    state.movable[vertex] || return false
    cavity=_gopt_vertex_cavity(state,index,vertex)
    cavity===nothing && return false
    x=0.0;y=0.0;z=0.0;volume_total=0.0;worst=1.0
    for index in cavity
        cell=state.cells[index];v=cell.vertices
        volume=abs(_gopt_volume(state,v))
        worst=min(worst,cell.quality)
        x+=0.25*(state.coords[1,v[1]]+state.coords[1,v[2]]+
            state.coords[1,v[3]]+state.coords[1,v[4]])*volume
        y+=0.25*(state.coords[2,v[1]]+state.coords[2,v[2]]+
            state.coords[2,v[3]]+state.coords[2,v[4]])*volume
        z+=0.25*(state.coords[3,v[1]]+state.coords[3,v[2]]+
            state.coords[3,v[3]]+state.coords[3,v[4]])*volume
        volume_total+=volume
    end
    old=_gopt_point(state,vertex)
    state.coords[1,vertex]=x/volume_total
    state.coords[2,vertex]=y/volume_total
    state.coords[3,vertex]=z/volume_total
    worst_after=1.0;volume_after=0.0
    qualities=Vector{Float64}(undef,length(cavity))
    for (position,index) in enumerate(cavity)
        quality,volume=_gopt_quality(state,state.cells[index].vertices)
        qualities[position]=quality;volume_after+=volume
        worst_after=min(worst_after,quality)
    end
    if abs(volume_after-volume_total)>1e-10*volume_total || worst_after<worst
        state.coords[1,vertex]=old[1];state.coords[2,vertex]=old[2];state.coords[3,vertex]=old[3]
        return false
    end
    for (position,index) in enumerate(cavity)
        state.cells[index].quality=qualities[position]
    end
    return true
end

function _gmsh_optimize_tetrahedra(mesh::Mesh;
        element_numbers=collect(UInt64,1:size(mesh.tets,2)),
        element_high_water=maximum(element_numbers;init=UInt64(0)),
        movable=nothing,threshold=0.3,embedded_faces=Set{NTuple{3,Int32}}(),
        embedded_edges=Set{Tuple{Int32,Int32}}(),positive=true)
    isfinite(threshold) || throw(ArgumentError("native optimizer threshold must be finite"))
    numbers=collect(UInt64,element_numbers)
    length(numbers)==size(mesh.tets,2) && allunique(numbers) || throw(ArgumentError(
        "native optimizer element identities must bind every input cell uniquely"))
    maximum(numbers;init=UInt64(0))<=element_high_water || throw(ArgumentError(
        "native optimizer element high-water mark does not cover its inputs"))
    movement=if movable===nothing
        fixed=falses(size(mesh.coords,2))
        for face in first(boundary_faces(mesh.tets)),node in face
            fixed[node]=true
        end
        .!fixed
    else
        length(movable)==size(mesh.coords,2) || throw(ArgumentError(
            "native optimizer movement mask must bind every input node"))
        BitVector(movable)
    end
    cells=_GmshOptimizeTet[]
    face_carriers=Dict(key=>0 for key in embedded_faces)
    edge_carriers=Dict(key=>0 for key in embedded_edges)
    embedded_verts=Set{Int32}(v for f in embedded_faces for v in f)
    state=_GmshOptimizeState(copy(mesh.coords),cells,UInt64(element_high_water),
        Set{NTuple{3,Int32}}(embedded_faces),Set{Tuple{Int32,Int32}}(embedded_edges),
        movement,face_carriers,edge_carriers,embedded_verts)
    for column in axes(mesh.tets,2)
        vertices=ntuple(slot->mesh.tets[slot,column],4)
        quality,_=_gopt_quality(state,vertices)
        push!(cells,_GmshOptimizeTet(vertices,zeros(Int32,4),numbers[column],
            mesh.tet_tag[column],quality,false))
        _gopt_note_carriers!(state,vertices,1)
    end
    _gopt_connect!(state,eachindex(cells))
    active=collect(Int32,eachindex(cells))
    swaps=0;relocations=0;previous_worst=0.0
    if threshold>0.0
        while true
            new_cells=Int32[]
            for index in active
                cell=cells[index]
                cell.deleted && continue
                if cell.quality<threshold
                    for edge in 1:6
                        slots=_GOPT_EDGES[_GOPT_EDGE_PERMUTATION[edge]]
                        key=_gopt_edge(cell.vertices[slots[1]],cell.vertices[slots[2]])
                        key in state.embedded_edges && continue
                        if _gopt_swap_edge!(state,new_cells,index,edge)
                            swaps+=1
                            break
                        end
                    end
                end
            end
            isempty(new_cells) && break
            for index in new_cells
                cells[index].deleted || push!(active,index)
            end
            for index in active
                cell=cells[index]
                cell.deleted && continue
                if cell.quality<threshold
                    for role in 1:4
                        _gopt_smooth_vertex!(state,index,role) && (relocations+=1)
                    end
                end
            end
            worst=1.0
            for index in active
                cells[index].deleted || (worst=min(worst,cells[index].quality))
            end
            previous_worst!=0.0 && worst-previous_worst<1e-6 && break
            previous_worst=worst
        end
    end
    kept=Int32[index for index in active if !cells[index].deleted]
    output=Matrix{Int32}(undef,4,length(kept));tags=UInt64[];physical=Int32[]
    for (column,index) in enumerate(kept)
        vertices=cells[index].vertices
        if positive && _gopt_volume(state,vertices)<0.0
            vertices=(vertices[2],vertices[1],vertices[3],vertices[4])
        end
        for slot in 1:4
            output[slot,column]=vertices[slot]
        end
        push!(tags,cells[index].number)
        push!(physical,cells[index].physical)
    end
    product=Mesh(state.coords;segs=mesh.segs,seg_tag=mesh.seg_tag,
        tris=mesh.tris,tri_tag=mesh.tri_tag,tets=output,
        tet_tag=physical)
    return (mesh=product,element_numbers=tags,element_high_water=state.element_max,
        edge_swaps=swaps,vertex_relocations=relocations)
end
