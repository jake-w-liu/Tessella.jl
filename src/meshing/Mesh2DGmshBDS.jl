# `delaunayizeBDS` on the BDS mesh `meshGenerator` builds from the initial
# PolyMesh (meshGFace.cpp, meshGFaceBDS.cpp, BDS.cpp). On a plane the initial
# mesh is already Delaunay in the face metric and the pass swaps nothing; on a
# curved chart it decides the cells the frontal insertion starts from, and
# Gmsh's swaps depend on the order of its edge vector, on which face of an
# edge is listed first and on the node order each swap gives its two new
# faces. All of that is kept here:
#
#   * every PolyMesh face, box faces included, becomes a BDS face in face
#     order (`add_triangle`), which fixes the edge vector;
#   * faces outside the domain are deleted and the triangle and edge vectors
#     are compacted with libstdc++'s bidirectional `std::partition`;
#   * the swap pass walks the edge vector, appended edges included, asks
#     `edgeSwapTestDelaunayAniso` (first fundamental form at the edge's
#     parametric middle, each four-vertex configuration once) and swaps
#     through `BDS_Mesh::swap_edge` with `BDS_SwapEdgeTestQuality(false)`.
#
# `ids` are the BDS point numbers (`all_vertices` rank, box points -1..-4);
# the result lists the surviving faces in BDS2GMSH order by their nodes as
# `getNodes` returns them, and the swaps in the order they were made.
struct _GmshBDSResult
    faces::Vector{NTuple{3,Int32}}
    swaps::Vector{NTuple{4,Int32}}
end

mutable struct _GmshBDS
    u::Vector{Float64}
    v::Vector{Float64}
    id::Vector{Int}
    degree::Vector{Int32}
    edge_p1::Vector{Int32}
    edge_p2::Vector{Int32}
    edge_faces::Vector{Vector{Int32}}
    edge_deleted::BitVector
    edge_constrained::BitVector
    face_edges::Vector{NTuple{3,Int32}}
    face_deleted::BitVector
    edges::Vector{Int32}
    triangles::Vector{Int32}
    lookup::Dict{NTuple{2,Int32},Int32}
end

@inline _gmsh_bds_key(a::Int32,b::Int32)=a<b ? (a,b) : (b,a)

function _gmsh_bds_add_edge!(m::_GmshBDS,a::Int32,b::Int32)
    key=_gmsh_bds_key(a,b)
    found=get(m.lookup,key,Int32(0))
    found!=0 && !m.edge_deleted[found] && return found
    return _gmsh_bds_new_edge!(m,a,b)
end

function _gmsh_bds_new_edge!(m::_GmshBDS,a::Int32,b::Int32)
    m.id[a]<m.id[b] || ((a,b)=(b,a))
    push!(m.edge_p1,a);push!(m.edge_p2,b);push!(m.edge_faces,Int32[])
    push!(m.edge_deleted,false);push!(m.edge_constrained,false)
    e=Int32(length(m.edge_p1))
    push!(m.edges,e)
    m.degree[a]+=1;m.degree[b]+=1
    m.lookup[_gmsh_bds_key(a,b)]=e
    return e
end

function _gmsh_bds_new_face!(m::_GmshBDS,e1::Int32,e2::Int32,e3::Int32)
    push!(m.face_edges,(e1,e2,e3));push!(m.face_deleted,false)
    f=Int32(length(m.face_edges))
    push!(m.edge_faces[e1],f);push!(m.edge_faces[e2],f);push!(m.edge_faces[e3],f)
    push!(m.triangles,f)
    return f
end

function _gmsh_bds_del_face!(m::_GmshBDS,f::Int32)
    for e in m.face_edges[f]
        list=m.edge_faces[e]
        index=findfirst(==(f),list)
        index===nothing || deleteat!(list,index)
    end
    m.face_deleted[f]=true
    return nothing
end

function _gmsh_bds_del_edge!(m::_GmshBDS,e::Int32)
    m.degree[m.edge_p1[e]]-=1;m.degree[m.edge_p2[e]]-=1
    m.edge_deleted[e]=true
    key=_gmsh_bds_key(m.edge_p1[e],m.edge_p2[e])
    get(m.lookup,key,Int32(0))==e && delete!(m.lookup,key)
    return nothing
end

@inline function _gmsh_bds_common(m::_GmshBDS,a::Int32,b::Int32)
    p1=m.edge_p1[a];p2=m.edge_p2[a]
    (p1==m.edge_p1[b] || p1==m.edge_p2[b]) && return p1
    return p2
end

@inline function _gmsh_bds_nodes(m::_GmshBDS,f::Int32)
    e1,e2,e3=m.face_edges[f]
    return (_gmsh_bds_common(m,e1,e3),_gmsh_bds_common(m,e1,e2),
            _gmsh_bds_common(m,e2,e3))
end

@inline function _gmsh_bds_opposite(m::_GmshBDS,f::Int32,p1::Int32,p2::Int32)
    nodes=_gmsh_bds_nodes(m,f)
    nodes[1]!=p1 && nodes[1]!=p2 && return nodes[1]
    nodes[2]!=p1 && nodes[2]!=p2 && return nodes[2]
    return nodes[3]
end

function _gmsh_bds_face_edge(m::_GmshBDS,f::Int32,a::Int32,b::Int32)
    for e in m.face_edges[f]
        p1=m.edge_p1[e];p2=m.edge_p2[e]
        ((p1==a && p2==b) || (p1==b && p2==a)) && return e
    end
    return Int32(0)
end

# libstdc++ `std::partition` for bidirectional iterators: retained items from
# the back fill the first deleted slots.
function _gmsh_bds_partition!(items::Vector{Int32},deleted::BitVector)
    first=1;last=length(items)
    while true
        while first<=last && !deleted[items[first]]
            first+=1
        end
        first>last && break
        while first<=last && deleted[items[last]]
            last-=1
        end
        first>last && break
        items[first],items[last]=items[last],items[first]
        first+=1;last-=1
    end
    resize!(items,first-1)
    return items
end

@inline function _gmsh_bds_area(m::_GmshBDS,a::Int32,b::Int32,c::Int32)
    x1=m.u[a]-m.u[b];y1=m.v[a]-m.v[b]
    x2=m.u[a]-m.u[c];y2=m.v[a]-m.v[c]
    return 0.5*(x1*y2-y1*x2)
end

function _gmsh_bds_swap!(m::_GmshBDS,e::Int32)
    m.edge_deleted[e] && return false
    faces=m.edge_faces[e]
    length(faces)==2 || return false
    m.edge_constrained[e] && return false
    p1=m.edge_p1[e];p2=m.edge_p2[e]
    f0=faces[1];f1=faces[2]
    pts1=_gmsh_bds_nodes(m,f0)
    op0=_gmsh_bds_opposite(m,f0,p1,p2);op1=_gmsh_bds_opposite(m,f1,p1,p2)
    orientation=0
    for i in 1:3
        if pts1[i]==p1
            orientation=pts1[mod1(i+1,3)]==p2 ? 1 : -1
            break
        end
    end
    # BDS_SwapEdgeTestQuality(false, true): no quality test, but the small
    # triangle guard and the parametric convexity test.
    (m.degree[p1]<=3 || m.degree[p2]<=3) && return false
    s1=abs(_gmsh_bds_area(m,p1,p2,op0));s2=abs(_gmsh_bds_area(m,p1,p2,op1))
    s3=abs(_gmsh_bds_area(m,p1,op0,op1));s4=abs(_gmsh_bds_area(m,p2,op0,op1))
    abs(s1+s2-s3-s4)>1e-12*(s3+s4) && return false
    p1_op1=_gmsh_bds_face_edge(m,f0,p1,op0)
    op1_p2=_gmsh_bds_face_edge(m,f0,op0,p2)
    p1_op2=_gmsh_bds_face_edge(m,f1,p1,op1)
    op2_p2=_gmsh_bds_face_edge(m,f1,op1,p2)
    (p1_op1==p1_op2 || op2_p2==op1_p2) && return false
    _gmsh_bds_del_face!(m,f0)
    _gmsh_bds_del_face!(m,m.edge_faces[e][1])
    _gmsh_bds_del_edge!(m,e)
    new=_gmsh_bds_new_edge!(m,op0,op1)
    if orientation==1
        _gmsh_bds_new_face!(m,p1_op1,p1_op2,new)
        _gmsh_bds_new_face!(m,new,op2_p2,op1_p2)
    else
        _gmsh_bds_new_face!(m,p1_op2,p1_op1,new)
        _gmsh_bds_new_face!(m,op2_p2,new,op1_p2)
    end
    return true
end

# `inCircumCircleAniso(gf, p1, p2, p3, p4, metric)`: the four-point overload,
# with its radius-dependent tolerance.
function _gmsh_bds_in_circle(pa,pb,pc,pd,metric)
    center,radius2=_gmsh_circumcenter_metric(pa,pb,pc,metric)
    distance=_gmsh_metric_distance2(center,pd,metric)
    tolerance=radius2<=1e3 ? 1e-12 : radius2<=1e5 ? 1e-11 : 1e-9
    return distance<radius2-tolerance
end

"""
    _gmsh_bds_delaunayize(p, ids, interior_faces, constrained, ground)

Replay `meshGenerator`'s BDS construction from the PolyMesh `p` and run
`delaunayizeBDS`. `ids[v]` is the BDS number of PolyMesh vertex `v`,
`interior_faces[f]` whether face `f` is in the domain, `constrained` the
recovered model edges as PolyMesh vertex pairs and `ground(u,v)` the first
fundamental form `(E,F,G)`.
"""
function _gmsh_bds_delaunayize(p,ids::AbstractVector{<:Integer},
                               interior_faces::AbstractVector{Bool},
                               constrained,ground::G) where {G}
    n=length(p.x)
    m=_GmshBDS(copy(p.x),copy(p.y),collect(Int,ids),zeros(Int32,n),
        Int32[],Int32[],Vector{Int32}[],falses(0),falses(0),
        NTuple{3,Int32}[],falses(0),Int32[],Int32[],Dict{NTuple{2,Int32},Int32}())
    for f in eachindex(p.face_edge)
        h=p.face_edge[f]
        a=p.origin[h];b=p.origin[p.next[h]];c=p.origin[p.next[p.next[h]]]
        e1=_gmsh_bds_add_edge!(m,a,b)
        e2=_gmsh_bds_add_edge!(m,b,c)
        e3=_gmsh_bds_add_edge!(m,c,a)
        _gmsh_bds_new_face!(m,e1,e2,e3)
    end
    for (a,b) in constrained
        e=get(m.lookup,_gmsh_bds_key(Int32(a),Int32(b)),Int32(0))
        e==0 && throw(ArgumentError("native initial mesh: a recovered model edge is missing"))
        m.edge_constrained[e]=true
    end
    for f in m.triangles
        interior_faces[f] || _gmsh_bds_del_face!(m,f)
    end
    _gmsh_bds_partition!(m.triangles,m.face_deleted)
    for e in m.edges
        isempty(m.edge_faces[e]) && _gmsh_bds_del_edge!(m,e)
    end
    _gmsh_bds_partition!(m.edges,m.edge_deleted)
    examined=Set{NTuple{4,Int}}()
    swaps=NTuple{4,Int32}[]
    while true
        count=0
        index=1
        while index<=length(m.edges)
            e=m.edges[index];index+=1
            m.edge_deleted[e] && continue
            faces=m.edge_faces[e]
            length(faces)==2 || continue
            p1=m.edge_p1[e];p2=m.edge_p2[e]
            op0=_gmsh_bds_opposite(m,faces[1],p1,p2)
            op1=_gmsh_bds_opposite(m,faces[2],p1,p2)
            quad=Tuple(sort!(Int[m.id[p1],m.id[p2],m.id[op0],m.id[op1]]))
            quad in examined && continue
            push!(examined,quad)
            metric=ground(0.5*(m.u[p1]+m.u[p2]),0.5*(m.v[p1]+m.v[p2]))
            _gmsh_bds_in_circle((m.u[p1],m.v[p1]),(m.u[p2],m.v[p2]),
                (m.u[op0],m.v[op0]),(m.u[op1],m.v[op1]),metric) || continue
            if _gmsh_bds_swap!(m,e)
                count+=1
                push!(swaps,(p1,p2,op0,op1))
            end
        end
        count==0 && break
    end
    faces=NTuple{3,Int32}[]
    for f in m.triangles
        m.face_deleted[f] || push!(faces,_gmsh_bds_nodes(m,f))
    end
    return _GmshBDSResult(faces,swaps)
end
