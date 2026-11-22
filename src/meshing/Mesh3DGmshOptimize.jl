# SPDX-FileCopyrightText: 1997-2025 C. Geuzaine, J.-F. Remacle
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: MIT AND GPL-2.0-or-later
# Contains Gmsh 4.15.2 adaptations: mesh/meshGRegionLocalMeshMod.cpp; ring swap pattern order.
# Julia adaptation/changes: Tessella contributors, 2026; notices corrected 2026-11-14.
# Licenses apply to their respective portions; see LICENSES/README.md and NOTICES.md#gmsh.

# Native tetrahedron quality optimization. Storage order and element births
# belong to this operation; the caller retains its public node identities.

const _GOPT_EDGES=((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
const _GOPT_EDGE_FACES=((1,3),(1,2),(2,3),(1,4),(3,4),(2,4))
const _GOPT_FACES=((1,2,3),(1,3,4),(1,2,4),(2,3,4))
const _GOPT_VERTEX_FACES=((1,2,3),(1,3,4),(1,2,4),(2,3,4))
const _GOPT_FACE_APEX=(4,2,3,1)
const _GOPT_EDGE_PERMUTATION=(1,4,2,3,6,5)


# Literal native ring3–7 pattern order; pinned CPP SHA-256 abfb4c8cbda1356aeaf8cf017eea2ebf103acbbd2a22664a49664c0157bec658.
const _GOPT_SWAP_PATTERNS=Dict{Int,Tuple{Tuple,Tuple}}(
    3=>(((1,2,3),),((1,),)),
    4=>(((1,2,3),(3,4,1),(2,3,4),(4,1,2)),((1,2),(3,4))),
    5=>(((1,2,3),(1,3,4),(1,4,5),(1,2,5),(2,4,5),(2,3,4),(3,4,5),(1,3,5),(1,2,4),(2,3,5)),((1,2,3),(4,5,6),(1,7,8),(3,6,9),(4,7,10))),
    6=>(((1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,3,6),(3,5,6),(3,4,5),(1,4,6),(4,5,6),(1,3,5),(3,4,6),(2,3,4),(1,2,4),(1,2,6),(2,5,6),(2,4,5),(1,2,5),(2,4,6),(2,3,5),(2,3,6)),((1,2,3,4),(1,5,6,7),(1,2,8,9),(1,4,7,10),(1,5,9,11),(3,4,12,13),(12,14,15,16),(8,9,12,13),(4,12,16,17),(9,12,14,18),(7,14,15,19),(4,7,17,19),(6,7,14,20),(9,11,14,20))),
    7=>(((1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,6,7),(1,4,7),(4,6,7),(4,5,6),(1,5,7),(5,6,7),(1,4,6),(4,5,7),(1,3,5),(3,4,5),(1,3,7),(3,6,7),(3,5,6),(1,3,6),(3,5,7),(3,4,6),(3,4,7),(1,2,4),(2,3,4),(1,2,5),(2,4,5),(1,2,7),(2,6,7),(2,5,6),(1,2,6),(2,5,7),(2,4,6),(2,4,7),(2,3,5),(2,3,6),(2,3,7)),((1,2,3,4,5),(1,2,6,7,8),(1,2,3,9,10),(1,2,5,8,11),(1,2,6,10,12),(1,4,5,13,14),(1,14,15,16,17),(1,9,10,13,14),(1,5,14,17,18),(1,10,14,15,19),(1,8,15,16,20),(1,5,8,18,20),(1,7,8,15,21),(1,10,12,15,21),(3,4,5,22,23),(6,7,8,22,23),(3,9,10,22,23),(5,8,11,22,23),(6,10,12,22,23),(4,5,23,24,25),(23,25,26,27,28),(9,10,23,24,25),(5,23,25,28,29),(10,23,25,26,30),(8,23,26,27,31),(5,8,23,29,31),(7,8,23,26,32),(10,12,23,26,32),(4,5,14,24,33),(14,26,27,28,33),(9,10,14,24,33),(5,14,28,29,33),(10,14,26,30,33),(14,17,26,27,34),(5,14,17,29,34),(14,16,17,26,35),(10,14,19,26,35),(8,20,26,27,34),(5,8,20,29,34),(8,16,20,26,35),(7,8,21,26,35),(10,12,21,26,35))),
)

# The same table with concrete element types. The literal above nests tuples
# of different lengths, so iterating it dispatches dynamically on every step
# of the innermost swap loop.
const _GOPT_SWAP_TABLE=Dict{Int,Tuple{Vector{NTuple{3,Int}},Vector{Vector{Int}}}}(
    size=>(NTuple{3,Int}[triangle for triangle in triangles],
           Vector{Int}[collect(Int,triangulation) for triangulation in triangulations])
    for (size,(triangles,triangulations)) in _GOPT_SWAP_PATTERNS)

# `neighbors` is stored inline: a vector per cell was two more heap objects
# for every tetrahedron of the mesh and for every cell a swap creates. The
# record is immutable, so the cell vector holds it inline instead of one heap
# object per cell; the setters below replace a record in place.
struct _GmshOptimizeTet
    vertices::NTuple{4,Int32}
    neighbors::NTuple{4,Int32}
    number::UInt64
    physical::Int32
    quality::Float64
    deleted::Bool
end

@inline function _gopt_set_neighbor!(cells,index,slot,value)
    c=cells[index]
    cells[index]=_GmshOptimizeTet(c.vertices,Base.setindex(c.neighbors,Int32(value),slot),
        c.number,c.physical,c.quality,c.deleted)
    return nothing
end
@inline function _gopt_set_quality!(cells,index,quality)
    c=cells[index]
    cells[index]=_GmshOptimizeTet(c.vertices,c.neighbors,c.number,c.physical,quality,c.deleted)
    return nothing
end
@inline function _gopt_delete!(cells,index)
    c=cells[index]
    cells[index]=_GmshOptimizeTet(c.vertices,c.neighbors,c.number,c.physical,c.quality,true)
    return nothing
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

# Whether the cell has negative volume, by the filtered exact predicate.
function _gopt_negative(state,vertices)
    a=_gopt_point(state,vertices[1]);b=_gopt_point(state,vertices[2])
    c=_gopt_point(state,vertices[3]);d=_gopt_point(state,vertices[4])
    return _orient3_value(a,b,c,d)>0.0
end

# +1 or -1 when the exact orientation and the Float64 volume `_gopt_finish`
# orients by agree on the sign of the cell, 0 otherwise.
function _gopt_orientation(state,vertices)
    exact=_orient3_value(_gopt_point(state,vertices[1]),_gopt_point(state,vertices[2]),
        _gopt_point(state,vertices[3]),_gopt_point(state,vertices[4]))
    float=_gopt_volume(state,vertices)
    # `_orient3_value` is det(a-d,b-d,c-d), the negative of 6 × the volume.
    exact<0.0 && float>0.0 && return Int8(1)
    exact>0.0 && float<0.0 && return Int8(-1)
    return Int8(0)
end

@inline _gopt_face(vertices,face)=_sorted_face3(
    vertices[_GOPT_FACES[face][1]],vertices[_GOPT_FACES[face][2]],
    vertices[_GOPT_FACES[face][3]])
@inline _gopt_edge(first,last)=first<last ? (first,last) : (last,first)

# Reconnection retains neighbors outside the supplied local cavity. Embedded
# faces remain cuts in the graph, as in native connectTets.
function _gopt_connect!(state,indices)
    # Ingesting a whole mesh bonds every cell at once: sort the faces rather
    # than hash them.
    indices isa AbstractUnitRange && length(indices)>256 &&
        return _gopt_connect_sorted!(state,indices)
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
                _gopt_set_neighbor!(state.cells,previous[1],previous[2],index)
                _gopt_set_neighbor!(state.cells,index,face,previous[1])
            end
        end
    end
    return nothing
end

# The same bonds as the hashed pass above, from sorted faces. The owners of a
# face must stay in the order the cells were given (each later owner bonds
# with the first, exactly as the hashed pass leaves them), so the sort has to
# be stable. Faces are dealt into one bucket per lowest vertex by counting,
# which keeps that order, and each bucket (about twenty faces) is put in order
# by a stable insertion sort on its two remaining vertices: no comparison of
# triples and no merge buffer for the whole mesh.
function _gopt_connect_sorted!(state,indices::AbstractUnitRange)
    nv=size(state.coords,2)
    (isempty(indices) || last(indices)>=(1<<29)) &&
        return _gopt_connect_merge_sorted!(state,indices)
    embedded=!isempty(state.embedded_faces)
    stop=zeros(Int32,nv+2)
    @inbounds for index in indices
        cell=state.cells[index]
        cell.deleted && continue
        for face in 1:4
            key=_gopt_face(cell.vertices,face)
            embedded && key in state.embedded_faces && continue
            (0<key[1]<=nv) || return _gopt_connect_merge_sorted!(state,indices)
            stop[key[1]+2]+=Int32(1)
        end
    end
    @inbounds for vertex in 1:nv
        stop[vertex+2]+=stop[vertex+1]
    end
    total=Int(stop[nv+2])
    rest=Vector{UInt64}(undef,total)     # (middle, largest) vertex of a face
    who=Vector{Int32}(undef,total)       # 4*cell + face - 1
    @inbounds for index in indices
        cell=state.cells[index]
        cell.deleted && continue
        for face in 1:4
            key=_gopt_face(cell.vertices,face)
            embedded && key in state.embedded_faces && continue
            position=(stop[key[1]+1]+=Int32(1))
            rest[position]=(UInt64(key[2])<<32)|UInt64(key[3])
            who[position]=(Int32(index)<<2)|Int32(face-1)
        end
    end
    # After the fill, bucket `vertex` spans `stop[vertex]+1 : stop[vertex+1]`.
    @inbounds for vertex in 1:nv
        first=Int(stop[vertex])+1;last_=Int(stop[vertex+1])
        first<last_ || continue
        _sort_face_bucket!(rest,who,first,last_)
        position=first
        while position<=last_
            next=position+1
            head=who[position]
            while next<=last_ && rest[next]==rest[position]
                other=who[next]
                head_cell=head>>2;other_cell=other>>2
                _gopt_set_neighbor!(state.cells,head_cell,Int(head&Int32(3))+1,other_cell)
                _gopt_set_neighbor!(state.cells,other_cell,Int(other&Int32(3))+1,head_cell)
                next+=1
            end
            position=next
        end
    end
    return nothing
end

# The plain stable sort, for complexes the bucket pass does not take.
function _gopt_connect_merge_sorted!(state,indices::AbstractUnitRange)
    records=Vector{Tuple{NTuple{3,Int32},Int32,Int8}}()
    sizehint!(records,4*length(indices))
    for index in indices
        cell=state.cells[index]
        cell.deleted && continue
        for face in 1:4
            key=_gopt_face(cell.vertices,face)
            key in state.embedded_faces && continue
            push!(records,(key,Int32(index),Int8(face)))
        end
    end
    sort!(records;by=first,alg=MergeSort)
    position=1;total=length(records)
    @inbounds while position<=total
        next=position+1
        head=records[position]
        while next<=total && records[next][1]==head[1]
            other=records[next]
            _gopt_set_neighbor!(state.cells,head[2],Int(head[3]),other[2])
            _gopt_set_neighbor!(state.cells,other[2],Int(other[3]),head[2])
            next+=1
        end
        position=next
    end
    return nothing
end

# Containers of one edge swap, emptied and refilled by every attempt instead
# of allocated by it: a swap looks at a ring of three to seven cells, so its
# face tables are short vectors searched in order rather than hash tables.
struct _GoptSwapScratch
    cavity::Vector{Int32}
    ring::Vector{Int32}
    outside::Vector{Int32}
    faces::Vector{NTuple{3,Int32}}      # distinct faces of the cavity cells
    face_count::Vector{Int}
    face_apex::Vector{Int32}
    candidate::Vector{NTuple{3,Int32}}  # distinct faces of one reconnection
    candidate_count::Vector{Int}
    qualities_first::Vector{Float64}
    qualities_last::Vector{Float64}
    volumes_first::Vector{Float64}
    volumes_last::Vector{Float64}
    signs_first::Vector{Int8}
    signs_last::Vector{Int8}
    records::Vector{Tuple{NTuple{3,Int32},Int32,Int8}}
    # vertex relocation: the cells around a vertex, the walk's frames and
    # the trial qualities
    star::Vector{Int32}
    frames::Vector{Tuple{Int32,Int}}
    qualities::Vector{Float64}
end
_GoptSwapScratch()=_GoptSwapScratch(Int32[],Int32[],Int32[],NTuple{3,Int32}[],Int[],Int32[],
    NTuple{3,Int32}[],Int[],Float64[],Float64[],Float64[],Float64[],Int8[],Int8[],
    Tuple{NTuple{3,Int32},Int32,Int8}[],Int32[],Tuple{Int32,Int}[],Float64[])

# Position of `key` in `keys`, 0 when absent.
@inline function _gopt_find(keys::Vector{NTuple{3,Int32}},key::NTuple{3,Int32})
    @inbounds for index in eachindex(keys)
        keys[index]==key && return index
    end
    return 0
end

# `_gopt_connect!` for the few cells around one swap, without a hash table:
# the faces are put in order by a stable insertion sort, so every later owner
# of a face bonds with its first owner exactly as the hashed pass does.
function _gopt_connect_small!(state,indices,scratch::_GoptSwapScratch)
    records=scratch.records
    empty!(records)
    embedded=!isempty(state.embedded_faces)
    @inbounds for index in indices
        cell=state.cells[index]
        cell.deleted && continue
        for face in 1:4
            key=_gopt_face(cell.vertices,face)
            embedded && key in state.embedded_faces && continue
            push!(records,(key,Int32(index),Int8(face)))
        end
    end
    @inbounds for i in 2:length(records)
        record=records[i];j=i-1
        while j>=1 && records[j][1]>record[1]
            records[j+1]=records[j];j-=1
        end
        records[j+1]=record
    end
    position=1;total=length(records)
    @inbounds while position<=total
        head=records[position];next=position+1
        while next<=total && records[next][1]==head[1]
            other=records[next]
            if other[2]!=head[2]
                _gopt_set_neighbor!(state.cells,head[2],Int(head[3]),other[2])
                _gopt_set_neighbor!(state.cells,other[2],Int(other[3]),head[2])
            end
            next+=1
        end
        position=next
    end
    return nothing
end

function _gopt_edge_cavity(state,start,edge)
    scratch=_GoptSwapScratch()
    found=_gopt_edge_cavity!(scratch,state,start,edge)
    found===nothing && return nothing
    return found[1],found[2],scratch.cavity,scratch.outside,scratch.ring
end

# The cells around the edge `edge` of cell `start`, the ring of vertices
# opposite it and the cells next to that cavity, into `scratch`; returns the
# edge's two vertices, or `nothing` when the ring does not close.
function _gopt_edge_cavity!(scratch::_GoptSwapScratch,state,start,edge)
    cell=state.cells[start]
    local_edge=_GOPT_EDGE_PERMUTATION[edge]
    first=cell.vertices[_GOPT_EDGES[local_edge][1]]
    last=cell.vertices[_GOPT_EDGES[local_edge][2]]
    last_in_ring=cell.vertices[_GOPT_EDGES[7-local_edge][1]]
    cavity=scratch.cavity;ring=scratch.ring
    empty!(cavity);empty!(ring);empty!(scratch.outside)
    push!(cavity,Int32(start));push!(ring,last_in_ring)
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
    outside=scratch.outside
    for index in cavity, neighbor in state.cells[index].neighbors
        neighbor==0 || neighbor in outside || neighbor in cavity || push!(outside,neighbor)
    end
    return first,last
end

function _gopt_new_cell!(state,vertices,quality,physical)
    state.element_max=Base.checked_add(state.element_max,UInt64(1))
    length(state.cells)<typemax(Int32) || throw(ArgumentError(
        "native tetrahedron optimization exceeds Int32 storage capacity"))
    push!(state.cells,_GmshOptimizeTet(vertices,(Int32(0),Int32(0),Int32(0),Int32(0)),
        state.element_max,physical,quality,false))
    _gopt_note_carriers!(state,vertices,1)
    return Int32(length(state.cells))
end

function _gopt_swap_edge!(state,new_cells,start,edge,scratch::_GoptSwapScratch=_GoptSwapScratch())
    found=_gopt_edge_cavity!(scratch,state,start,edge)
    found===nothing && return false
    first,last=found
    cavity=scratch.cavity;outside=scratch.outside;ring=scratch.ring
    pattern=get(_GOPT_SWAP_TABLE,length(ring),nothing)
    pattern===nothing && return false
    # Embedded entities are hard constraints. Faces/edges carried by a cell
    # outside the cavity survive regardless; an entity whose every carrier is
    # deleted must be reproduced by the chosen reconnection. (Without any
    # embedded entity there is nothing to collect.)
    constrained=!isempty(state.embedded_face_carriers) ||
        !isempty(state.embedded_edge_carriers)
    constrained_faces=constrained ? Dict{NTuple{3,Int32},Int}() : nothing
    constrained_edges=constrained ? Dict{Tuple{Int32,Int32},Int}() : nothing
    if constrained
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
    end
    volume_ref=0.0;quality_ref=1.0
    # Face incidence of the cavity: each distinct face, how many cavity cells
    # carry it and the apex opposite it in the last of them.
    faces=scratch.faces;incidence=scratch.face_count;apexof=scratch.face_apex
    empty!(faces);empty!(incidence);empty!(apexof)
    for index in cavity
        cell=state.cells[index]
        volume_ref+=abs(_gopt_volume(state,cell.vertices))
        quality_ref=min(quality_ref,cell.quality)
        for face in 1:4
            key=_gopt_face(cell.vertices,face)
            apex=cell.vertices[_GOPT_FACE_APEX[face]]
            position=_gopt_find(faces,key)
            if position==0
                push!(faces,key);push!(incidence,1);push!(apexof,apex)
            else
                incidence[position]+=1;apexof[position]=apex
            end
        end
    end
    # Face incidence inside the cavity defines the certificate a reconnection
    # must satisfy: faces carried once form the cavity boundary and must be
    # reproduced exactly once by the new cells; every other new-cell face must
    # pair with another new cell; and the new apex on a boundary face must
    # stay on the cavity side of that face. The sidedness clause is what
    # keeps a cone tet of a non-convex ring ear from re-entering a region
    # already covered by a live cell outside the cavity — vertex-set
    # duplicates and folded-over output both enter through that path.
    # Upstream patterns on a complete convex star always satisfy the whole
    # certificate, so a refusal here only ever rejects a corrupting pattern.
    any(count->count>2,incidence) && return false
    # The boundary faces are the faces carried once (incidence 1).
    no_boundary=!any(==(1),incidence)
    triangles,triangulations=pattern
    count=length(triangles)
    qualities_first=resize!(scratch.qualities_first,count)
    qualities_last=resize!(scratch.qualities_last,count)
    volumes_first=resize!(scratch.volumes_first,count)
    volumes_last=resize!(scratch.volumes_last,count)
    # Orientation of every candidate cell, as the sign shared by its exact
    # determinant and its Float64 one (0 when they disagree or vanish).
    signs_first=resize!(scratch.signs_first,count);signs_last=resize!(scratch.signs_last,count)
    for index in eachindex(triangles)
        slots=triangles[index]
        a=ring[slots[1]];b=ring[slots[2]];c=ring[slots[3]]
        qualities_first[index],volumes_first[index]=_gopt_quality(state,(a,b,c,first))
        qualities_last[index],volumes_last[index]=_gopt_quality(state,(a,b,c,last))
        signs_first[index]=_gopt_orientation(state,(a,b,c,first))
        signs_last[index]=_gopt_orientation(state,(c,b,a,last))
    end
    best=-1.0;best_index=1
    for (index,triangulation) in enumerate(triangulations)
        quality=1.0;volume=0.0
        # Gmsh's gamma and volume test take absolute volumes, so a pattern
        # holding one flat or inverted cell of near-zero volume passes them
        # whenever the cavity was already that poor (cone.geo's Volume 2: four
        # coplanar boundary vertices). A reconnection of a valid cavity orients
        # every new cell alike; anything else is refused.
        orientation=signs_first[triangulation[1]]
        for triangle in triangulation
            quality=min(quality,qualities_first[triangle])
            quality=min(quality,qualities_last[triangle])
            volume+=volumes_first[triangle]+volumes_last[triangle]
            (signs_first[triangle]==orientation && signs_last[triangle]==orientation) ||
                (orientation=Int8(0))
        end
        orientation==0 && (quality=-1.0)
        abs(volume-volume_ref)>1e-10*(volume+volume_ref) && (quality=-1.0)
        if quality!=-1.0
            counts=scratch.candidate_count;keys_=scratch.candidate
            empty!(counts);empty!(keys_)
            for triangle in triangulation
                slots=triangles[triangle]
                a=ring[slots[1]];b=ring[slots[2]];c=ring[slots[3]]
                for vertices in ((a,b,c,first),(c,b,a,last)),face in 1:4
                    key=_gopt_face(vertices,face)
                    position=_gopt_find(keys_,key)
                    if position==0
                        push!(keys_,key);push!(counts,1)
                    else
                        counts[position]+=1
                    end
                end
            end
            # A tiling of the cavity: each face of the new cells carried once
            # is a cavity boundary face, each carried twice is not, and every
            # boundary face is carried exactly once.
            tiling=isempty(keys_) ? no_boundary : true
            for (position,key) in enumerate(keys_)
                where=_gopt_find(faces,key)
                boundary=where!=0 && incidence[where]==1
                tiling=counts[position]==1 ? boundary :
                    counts[position]==2 ? !boundary : false
                tiling || break
            end
            if tiling
                for (where,key) in enumerate(faces)
                    incidence[where]==1 || continue
                    position=_gopt_find(keys_,key)
                    (position!=0 && counts[position]==1) || (tiling=false;break)
                end
            end
            if tiling
                for triangle in triangulation
                    slots=triangles[triangle]
                    a=ring[slots[1]];b=ring[slots[2]];c=ring[slots[3]]
                    for vertices in ((a,b,c,first),(c,b,a,last)),face in 1:4
                        key=_gopt_face(vertices,face)
                        where=_gopt_find(faces,key)
                        (where!=0 && incidence[where]==1) || continue
                        interior=apexof[where]
                        fa=_gopt_point(state,key[1]);fb=_gopt_point(state,key[2])
                        fc=_gopt_point(state,key[3])
                        s_new=_orient3_value(fa,fb,fc,
                            _gopt_point(state,vertices[_GOPT_FACE_APEX[face]]))
                        s_old=_orient3_value(fa,fb,fc,_gopt_point(state,interior))
                        s_new*s_old<0 && (tiling=false)
                    end
                    tiling || break
                end
            end
            tiling || (quality=-1.0)
        end
        if quality!=-1.0 &&
                !((constrained_faces===nothing || isempty(constrained_faces)) &&
                  (constrained_edges===nothing || isempty(constrained_edges)) &&
                  isempty(state.embedded_verts))
            uncovered_faces=constrained_faces===nothing ? Set{NTuple{3,Int32}}() :
                Set(keys(constrained_faces))
            uncovered_edges=constrained_edges===nothing ? Set{Tuple{Int32,Int32}}() :
                Set(keys(constrained_edges))
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
        _gopt_delete!(state.cells,index)
        _gopt_note_carriers!(state,state.cells[index].vertices,-1)
    end
    _gopt_connect_small!(state,outside,scratch)
    return true
end

# Use a stack of continuation frames to retain native depth-first visitation
# order without imposing the native recursion's process-stack limit.
function _gopt_vertex_cavity(state,start,vertex)
    cavity=Int32[];frames=Tuple{Int32,Int}[]
    return _gopt_vertex_cavity!(cavity,frames,state,start,vertex)
end

# The same walk into reused vectors. The cells around one vertex are a few
# dozen, so "already seen" is a search of `cavity` rather than a hash set.
function _gopt_vertex_cavity!(cavity::Vector{Int32},frames::Vector{Tuple{Int32,Int}},
                              state,start,vertex)
    empty!(cavity);empty!(frames)
    push!(cavity,Int32(start));push!(frames,(Int32(start),1))
    while !isempty(frames)
        index,position=pop!(frames)
        cell=state.cells[index]
        cell.deleted && return nothing
        role=findfirst(==(vertex),cell.vertices)
        role===nothing && return nothing
        position>3 && continue
        push!(frames,(index,position+1))
        neighbor=cell.neighbors[_GOPT_VERTEX_FACES[role][position]]
        neighbor==0 || neighbor in cavity || begin
            push!(cavity,neighbor);push!(frames,(neighbor,1))
        end
    end
    return cavity
end

function _gopt_smooth_vertex!(state,index,role,scratch::_GoptSwapScratch=_GoptSwapScratch())
    cell=state.cells[index]
    cell.deleted && return false
    vertex=cell.vertices[role]
    state.movable[vertex] || return false
    cavity=_gopt_vertex_cavity!(scratch.star,scratch.frames,state,index,vertex)
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
    qualities=resize!(scratch.qualities,length(cavity))
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
        _gopt_set_quality!(state.cells,index,qualities[position])
    end
    return true
end

function _gopt_ingest(mesh::Mesh,element_numbers,element_high_water,movable,
                      embedded_faces,embedded_edges)
    numbers=collect(UInt64,element_numbers)
    length(numbers)==size(mesh.tets,2) && allunique(numbers) || throw(ArgumentError(
        "native optimizer element identities must bind every input cell uniquely"))
    maximum(numbers;init=UInt64(0))<=element_high_water || throw(ArgumentError(
        "native optimizer element high-water mark does not cover its inputs"))
    all(face->length(face)==3 && allunique(face) &&
        all(row->row isa Integer && !(row isa Bool) && 0<row<=size(mesh.coords,2),face),
        embedded_faces) || throw(ArgumentError("native optimizer embedded faces have invalid current rows"))
    all(edge->length(edge)==2 && allunique(edge) &&
        all(row->row isa Integer && !(row isa Bool) && 0<row<=size(mesh.coords,2),edge),
        embedded_edges) || throw(ArgumentError("native optimizer embedded edges have invalid current rows"))
    # The owned constructor supplies cyclic faces and directed edges; only
    # internal incidence/cut keys forget that presentation.
    embedded_faces=Set{NTuple{3,Int32}}(_sorted_face3(face...) for face in embedded_faces)
    embedded_edges=Set{Tuple{Int32,Int32}}(_gopt_edge(edge...) for edge in embedded_edges)
    movement=if movable===nothing
        fixed=falses(size(mesh.coords,2))
        for face in _single_tet_faces(mesh.tets),node in face
            fixed[node]=true
        end
        .!fixed
    else
        length(movable)==size(mesh.coords,2) || throw(ArgumentError(
            "native optimizer movement mask must bind every input node"))
        BitVector(movable)
    end
    cells=_GmshOptimizeTet[]
    # Swaps append the cells they make; room for them saves copying the
    # whole vector when it would otherwise double.
    sizehint!(cells,size(mesh.tets,2)+size(mesh.tets,2)>>2)
    face_carriers=Dict(key=>0 for key in embedded_faces)
    edge_carriers=Dict(key=>0 for key in embedded_edges)
    embedded_verts=Set{Int32}(v for f in embedded_faces for v in f)
    state=_GmshOptimizeState(copy(mesh.coords),cells,UInt64(element_high_water),
        Set{NTuple{3,Int32}}(embedded_faces),Set{Tuple{Int32,Int32}}(embedded_edges),
        movement,face_carriers,edge_carriers,embedded_verts)
    for column in axes(mesh.tets,2)
        vertices=ntuple(slot->mesh.tets[slot,column],4)
        quality,_=_gopt_quality(state,vertices)
        push!(cells,_GmshOptimizeTet(vertices,(Int32(0),Int32(0),Int32(0),Int32(0)),numbers[column],
            mesh.tet_tag[column],quality,false))
        _gopt_note_carriers!(state,vertices,1)
    end
    _gopt_connect!(state,eachindex(cells))
    return state,collect(Int32,eachindex(cells))
end

function _gopt_finish(state,mesh,active,swaps,relocations;positive=true)
    kept=Int32[]
    for index in active
        state.cells[index].deleted || push!(kept,index)
    end
    output=Matrix{Int32}(undef,4,length(kept));tags=UInt64[];physical=Int32[]
    for (column,index) in enumerate(kept)
        vertices=state.cells[index].vertices
        # The exact sign decides: on a cell of near-zero volume the Float64
        # determinant can have the other sign, and a swap by it inverts the
        # cell (cone.geo's Volume 2, four almost coplanar boundary vertices).
        if positive && _gopt_negative(state,vertices)
            vertices=(vertices[2],vertices[1],vertices[3],vertices[4])
        end
        for slot in 1:4
            output[slot,column]=vertices[slot]
        end
        push!(tags,state.cells[index].number)
        push!(physical,state.cells[index].physical)
    end
    product=Mesh(state.coords;segs=mesh.segs,seg_tag=mesh.seg_tag,
        tris=mesh.tris,tri_tag=mesh.tri_tag,tets=output,
        tet_tag=physical)
    return (mesh=product,element_numbers=tags,element_high_water=state.element_max,
        edge_swaps=swaps,vertex_relocations=relocations)
end

function _gmsh_optimize_tetrahedra(mesh::Mesh;
        element_numbers=collect(UInt64,1:size(mesh.tets,2)),
        element_high_water=maximum(element_numbers;init=UInt64(0)),
        movable=nothing,threshold=0.3,embedded_faces=Set{NTuple{3,Int32}}(),
        embedded_edges=Set{Tuple{Int32,Int32}}(),positive=true,relocate=true)
    isfinite(threshold) || throw(ArgumentError("native optimizer threshold must be finite"))
    state,active=_gopt_ingest(mesh,element_numbers,element_high_water,movable,
        embedded_faces,embedded_edges)
    cells=state.cells
    swaps=0;relocations=0;previous_worst=0.0
    scratch=_GoptSwapScratch()
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
                        if _gopt_swap_edge!(state,new_cells,index,edge,scratch)
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
            # optimizeMesh's relocate pass only runs when the region holds no
            # hexahedra, prisms or pyramids; hybrid regions get swaps alone.
            if relocate
                for index in active
                    cell=cells[index]
                    cell.deleted && continue
                    if cell.quality<threshold
                        for role in 1:4
                            _gopt_smooth_vertex!(state,index,role,scratch) && (relocations+=1)
                        end
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
    return _gopt_finish(state,mesh,active,swaps,relocations;positive)
end

# Independent Netgen-role volume pass. Ng_OptimizeVolumeMesh drives Netgen's
# optsteps3d rounds of "cmdmustm" (CombineImprove, ImproveMesh, SplitImprove,
# SwapImprove, SwapImprove2, ImproveMesh). This pass reproduces the certified
# swap sweep and per-point quality line search roles; CombineImprove node
# removal, SplitImprove midpoint insertion and the upstream-disabled surface
# swap are not implemented, and it never drops element coverage.
function _ng_vertex_quality(state,cavity,vertex,trial)
    worst=1.0
    for index in cavity
        cell=state.cells[index]
        cell.deleted && continue
        vertices=cell.vertices
        points=ntuple(slot->vertices[slot]==vertex ? trial :
            _gopt_point(state,vertices[slot]),4)
        value,_=_gopt_gamma(points...)
        orient3(points...)>0 && (value=-value)
        value<worst && (worst=value)
    end
    return worst
end

function _ng_vertex_objective(state,cavity,vertex,old,target,xi)
    trial=((1.0-xi)*old[1]+xi*target[1],
           (1.0-xi)*old[2]+xi*target[2],
           (1.0-xi)*old[3]+xi*target[3])
    return _ng_vertex_quality(state,cavity,vertex,trial)
end

function _ng_smooth_vertex!(state,index,role)
    cell=state.cells[index]
    cell.deleted && return false
    vertex=cell.vertices[role]
    state.movable[vertex] || return false
    cavity=_gopt_vertex_cavity(state,index,vertex)
    cavity===nothing && return false
    x=0.0;y=0.0;z=0.0;count=0
    for index in cavity
        vertices=state.cells[index].vertices
        for slot in 1:4
            node=vertices[slot]
            node==vertex && continue
            point=_gopt_point(state,node)
            x+=point[1];y+=point[2];z+=point[3];count+=1
        end
    end
    count==0 && return false
    target=(x/count,y/count,z/count)
    old=_gopt_point(state,vertex)
    lambda=0.5*(sqrt(5.0)-1.0);mu=0.5*(3.0-sqrt(5.0))
    a=0.0;b=2.0;tolerance=1e-2
    x1=b-lambda*(b-a);x2=a+lambda*(b-a)
    f1=_ng_vertex_objective(state,cavity,vertex,old,target,x1)
    f2=_ng_vertex_objective(state,cavity,vertex,old,target,x2)
    while !(abs(b-a)<tolerance)
        if f1<f2
            a=x1
            abs(b-a)<tolerance && break
            x1=x2;f1=f2;x2=b-mu*(b-a)
            f2=_ng_vertex_objective(state,cavity,vertex,old,target,x2)
        else
            b=x2
            abs(b-a)<tolerance && break
            x2=x1;f2=f1;x1=a+mu*(b-a)
            f1=_ng_vertex_objective(state,cavity,vertex,old,target,x1)
        end
    end
    # Strict improvement only: a plateau or zero-length move reports no work,
    # which keeps the caller's round loop a bounded strict-descent sequence.
    before=_ng_vertex_objective(state,cavity,vertex,old,target,0.0)
    after=_ng_vertex_objective(state,cavity,vertex,old,target,a)
    after<=before && return false
    point=((1.0-a)*old[1]+a*target[1],
           (1.0-a)*old[2]+a*target[2],
           (1.0-a)*old[3]+a*target[3])
    state.coords[1,vertex]=point[1];state.coords[2,vertex]=point[2]
    state.coords[3,vertex]=point[3]
    for index in cavity
        cell=state.cells[index]
        cell.deleted && continue
        quality,_=_gopt_quality(state,cell.vertices)
        _gopt_set_quality!(state.cells,index,quality)
    end
    return true
end

function _netgen_optimize_tetrahedra(mesh::Mesh;
        element_numbers=collect(UInt64,1:size(mesh.tets,2)),
        element_high_water=maximum(element_numbers;init=UInt64(0)),
        movable=nothing,embedded_faces=Set{NTuple{3,Int32}}(),
        embedded_edges=Set{Tuple{Int32,Int32}}(),positive=true,rounds=3)
    rounds>=0 || throw(ArgumentError("native Netgen-role rounds must be non-negative"))
    state,active=_gopt_ingest(mesh,element_numbers,element_high_water,movable,
        embedded_faces,embedded_edges)
    cells=state.cells
    swaps=0;relocations=0
    scratch=_GoptSwapScratch()
    for _ in 1:rounds
        changed=false
        new_cells=Int32[]
        for index in active
            cells[index].deleted && continue
            for role in 1:4
                _ng_smooth_vertex!(state,index,role) &&
                    (relocations+=1;changed=true)
            end
        end
        for index in active
            cell=cells[index]
            cell.deleted && continue
            for edge in 1:6
                slots=_GOPT_EDGES[_GOPT_EDGE_PERMUTATION[edge]]
                key=_gopt_edge(cell.vertices[slots[1]],cell.vertices[slots[2]])
                key in state.embedded_edges && continue
                _gopt_swap_edge!(state,new_cells,index,edge,scratch) &&
                    (swaps+=1;changed=true;break)
            end
        end
        for index in new_cells
            cells[index].deleted || push!(active,index)
        end
        for index in active
            cells[index].deleted && continue
            for role in 1:4
                _ng_smooth_vertex!(state,index,role) &&
                    (relocations+=1;changed=true)
            end
        end
        changed || break
    end
    return _gopt_finish(state,mesh,active,swaps,relocations;positive)
end
