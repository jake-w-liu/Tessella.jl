# The BDS mesh `meshGenerator` builds from the initial PolyMesh (meshGFace.cpp,
# meshGFaceBDS.cpp, BDS.cpp), with the two passes that run on it:
#
#   * `delaunayizeBDS`, which on a curved chart decides the cells the frontal
#     insertion starts from (on a plane the initial mesh is already Delaunay
#     in the face metric and it swaps nothing);
#   * `refineMeshBDS`, the MeshAdapt algorithm: split, swap, collapse and
#     smooth passes until every edge is between 0.7 and 1.4 sizes long.
#
# Gmsh's results depend on the order of its edge and triangle vectors, on the
# order of each point's edge list and each edge's face list, on the node order
# `getNodes` derives from a face's edges, and on libstdc++'s `std::partition`
# when the vectors are compacted. All of that is kept: every PolyMesh face, box
# faces included, becomes a BDS face in face order (`add_triangle`), which
# fixes the edge vector; points are numbered by their rank in `all_vertices`
# (node number order), the box points -1..-4, inserted points from
# `MAXPOINTNUMBER` up.

const _GMSH_BDS_MAX_LC=1e22

mutable struct _GmshBDS
    u::Vector{Float64}
    v::Vector{Float64}
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    lc::Vector{Float64}
    lcbgm::Vector{Float64}
    id::Vector{Int}
    # BDS_GeomEntity of each point and edge: classification degree (0 vertex,
    # 1 curve, 2 face; -1 none) and tag.
    pdeg::Vector{Int8}
    ptag::Vector{Int}
    pedges::Vector{Vector{Int32}}
    pdeleted::BitVector
    edge_p1::Vector{Int32}
    edge_p2::Vector{Int32}
    edge_faces::Vector{Vector{Int32}}
    edge_deleted::BitVector
    edge_deg::Vector{Int8}
    face_edges::Vector{NTuple{3,Int32}}
    face_deleted::BitVector
    face_g::BitVector
    edges::Vector{Int32}
    triangles::Vector{Int32}
    by_id::Dict{Int,Int32}
    maxpoint::Int
end

struct _GmshBDSResult
    faces::Vector{NTuple{3,Int32}}
    swaps::Vector{NTuple{4,Int32}}
end

function _gmsh_bds_add_point!(m::_GmshBDS,id::Int,u,v,X,Y,Z;deg=Int8(-1),tag=0)
    push!(m.u,u);push!(m.v,v);push!(m.X,X);push!(m.Y,Y);push!(m.Z,Z)
    push!(m.lc,0.0);push!(m.lcbgm,0.0);push!(m.id,id)
    push!(m.pdeg,Int8(deg));push!(m.ptag,tag);push!(m.pedges,Int32[])
    push!(m.pdeleted,false)
    p=Int32(length(m.u));m.by_id[id]=p
    m.maxpoint=max(m.maxpoint,id)
    return p
end

function _gmsh_bds_find_edge(m::_GmshBDS,a::Int32,b::Int32)
    for e in m.pedges[a]
        (m.edge_p1[e]==b || m.edge_p2[e]==b) && return e
    end
    return Int32(0)
end

function _gmsh_bds_new_edge!(m::_GmshBDS,a::Int32,b::Int32)
    m.id[a]<m.id[b] || ((a,b)=(b,a))
    push!(m.edge_p1,a);push!(m.edge_p2,b);push!(m.edge_faces,Int32[])
    push!(m.edge_deleted,false);push!(m.edge_deg,Int8(-1))
    e=Int32(length(m.edge_p1))
    push!(m.edges,e)
    push!(m.pedges[a],e);push!(m.pedges[b],e)
    return e
end

function _gmsh_bds_add_edge!(m::_GmshBDS,a::Int32,b::Int32)
    found=_gmsh_bds_find_edge(m,a,b)
    found!=0 && return found
    return _gmsh_bds_new_edge!(m,a,b)
end

function _gmsh_bds_new_face!(m::_GmshBDS,e1::Int32,e2::Int32,e3::Int32,g::Bool=false)
    push!(m.face_edges,(e1,e2,e3));push!(m.face_deleted,false);push!(m.face_g,g)
    f=Int32(length(m.face_edges))
    push!(m.edge_faces[e1],f);push!(m.edge_faces[e2],f);push!(m.edge_faces[e3],f)
    push!(m.triangles,f)
    return f
end

_gmsh_bds_add_triangle!(m::_GmshBDS,a::Int32,b::Int32,c::Int32)=
    _gmsh_bds_new_face!(m,_gmsh_bds_add_edge!(m,a,b),_gmsh_bds_add_edge!(m,b,c),
                        _gmsh_bds_add_edge!(m,c,a))

function _gmsh_bds_del_face!(m::_GmshBDS,f::Int32)
    for e in m.face_edges[f]
        filter!(!=(f),m.edge_faces[e])
    end
    m.face_deleted[f]=true
    return nothing
end

function _gmsh_bds_del_edge!(m::_GmshBDS,e::Int32)
    filter!(!=(e),m.pedges[m.edge_p1[e]]);filter!(!=(e),m.pedges[m.edge_p2[e]])
    m.edge_deleted[e]=true
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

@inline _gmsh_bds_other(m::_GmshBDS,e::Int32,p::Int32)=
    m.edge_p1[e]==p ? m.edge_p2[e] : m.edge_p1[e]

function _gmsh_bds_face_edge(m::_GmshBDS,f::Int32,a::Int32,b::Int32)
    for e in m.face_edges[f]
        p1=m.edge_p1[e];p2=m.edge_p2[e]
        ((p1==a && p2==b) || (p1==b && p2==a)) && return e
    end
    return Int32(0)
end

# `BDS_Point::getTriangles`: the faces of each edge in edge-list order.
function _gmsh_bds_point_faces(m::_GmshBDS,p::Int32)
    faces=Int32[]
    for e in m.pedges[p],f in m.edge_faces[e]
        f in faces || push!(faces,f)
    end
    return faces
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

_gmsh_bds_cleanup!(m::_GmshBDS)=(_gmsh_bds_partition!(m.triangles,m.face_deleted);
    _gmsh_bds_partition!(m.edges,m.edge_deleted);nothing)

# `surface_triangle_param` without degenerated points.
@inline function _gmsh_bds_area(m::_GmshBDS,a::Int32,b::Int32,c::Int32)
    x1=m.u[a]-m.u[b];y1=m.v[a]-m.v[b]
    x2=m.u[a]-m.u[c];y2=m.v[a]-m.v[c]
    return 0.5*(x1*y2-y1*x2)
end

@inline function _gmsh_bds_norme(a)
    modulus=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
    modulus==0.0 && return a
    inverse=1.0/modulus
    return (a[1]*inverse,a[2]*inverse,a[3]*inverse)
end

@inline _gmsh_bds_prodve(a,b)=(a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],
                               a[1]*b[2]-a[2]*b[1])

# `qmTriangle::gamma` on BDS coordinates.
function _gmsh_bds_gamma(m::_GmshBDS,pa::Int32,pb::Int32,pc::Int32)
    a=_gmsh_bds_norme((m.X[pc]-m.X[pb],m.Y[pc]-m.Y[pb],m.Z[pc]-m.Z[pb]))
    b=_gmsh_bds_norme((m.X[pa]-m.X[pc],m.Y[pa]-m.Y[pc],m.Z[pa]-m.Z[pc]))
    c=_gmsh_bds_norme((m.X[pb]-m.X[pa],m.Y[pb]-m.Y[pa],m.Z[pb]-m.Z[pa]))
    va=_gmsh_bds_prodve(b,c);sina=sqrt(va[1]*va[1]+va[2]*va[2]+va[3]*va[3])
    vb=_gmsh_bds_prodve(c,a);sinb=sqrt(vb[1]*vb[1]+vb[2]*vb[2]+vb[3]*vb[3])
    vc=_gmsh_bds_prodve(a,b);sinc=sqrt(vc[1]*vc[1]+vc[2]*vc[2]+vc[3]*vc[3])
    (sina==0.0 && sinb==0.0 && sinc==0.0) && return 0.0
    return 2*(2*sina*sinb*sinc/(sina+sinb+sinc))
end

# `normal_triangle`: `vector_triangle` then `norme`.
function _gmsh_bds_normal(m::_GmshBDS,p1::Int32,p2::Int32,p3::Int32)
    a=(m.X[p1]-m.X[p2],m.Y[p1]-m.Y[p2],m.Z[p1]-m.Z[p2])
    b=(m.X[p1]-m.X[p3],m.Y[p1]-m.Y[p3],m.Z[p1]-m.Z[p3])
    c=(a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
    return _gmsh_bds_norme(c)
end

# `BDS_Mesh::swap_edge`. `quality` selects `BDS_SwapEdgeTestQuality(true)`
# (delaunayizeBDS asks with `false`); the small-triangle guard and the
# parametric convexity test always apply.
function _gmsh_bds_swap!(m::_GmshBDS,e::Int32,quality::Bool)
    m.edge_deleted[e] && return false
    faces=m.edge_faces[e]
    length(faces)==2 || return false
    m.edge_deg[e]==1 && return false
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
    if quality
        qa1,qa2,qb1,qb2=orientation==1 ?
            (_gmsh_bds_gamma(m,p1,p2,op0),_gmsh_bds_gamma(m,p2,p1,op1),
             _gmsh_bds_gamma(m,p1,op1,op0),_gmsh_bds_gamma(m,op1,p2,op0)) :
            (_gmsh_bds_gamma(m,p2,p1,op0),_gmsh_bds_gamma(m,p1,p2,op1),
             _gmsh_bds_gamma(m,p1,op0,op1),_gmsh_bds_gamma(m,op1,op0,p2))
        min(qb1,qb2)>min(qa1,qa2) || return false
    end
    for p in (p1,p2)
        d=m.pdeg[p];n=length(m.pedges[p])
        d==2 && n<=4 && return false
        0<=d<2 && n<=3 && return false
    end
    s1=abs(_gmsh_bds_area(m,p1,p2,op0));s2=abs(_gmsh_bds_area(m,p1,p2,op1))
    s3=abs(_gmsh_bds_area(m,p1,op0,op1));s4=abs(_gmsh_bds_area(m,p2,op0,op1))
    abs(s1+s2-s3-s4)>1e-12*(s3+s4) && return false
    p1_op1=_gmsh_bds_face_edge(m,f0,p1,op0)
    op1_p2=_gmsh_bds_face_edge(m,f0,op0,p2)
    p1_op2=_gmsh_bds_face_edge(m,f1,p1,op1)
    op2_p2=_gmsh_bds_face_edge(m,f1,op1,p2)
    (p1_op1==p1_op2 || op2_p2==op1_p2) && return false
    g1=m.face_g[f0];g2=m.face_g[f1]
    _gmsh_bds_del_face!(m,f0)
    _gmsh_bds_del_face!(m,m.edge_faces[e][1])
    ge=m.edge_deg[e]
    _gmsh_bds_del_edge!(m,e)
    new=_gmsh_bds_new_edge!(m,op0,op1)
    if orientation==1
        _gmsh_bds_new_face!(m,p1_op1,p1_op2,new,g1)
        _gmsh_bds_new_face!(m,new,op2_p2,op1_p2,g2)
    else
        _gmsh_bds_new_face!(m,p1_op2,p1_op1,new,g1)
        _gmsh_bds_new_face!(m,op2_p2,new,op1_p2,g2)
    end
    m.edge_deg[new]=ge
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
    _gmsh_bds_build(p, ids, interior_faces, constrained; lift, owners)

`meshGenerator`'s BDS mesh of the PolyMesh `p`: `ids[v]` is the BDS number
of PolyMesh vertex `v`, `interior_faces[f]` whether face `f` is in the
domain, `constrained` the recovered model edges as PolyMesh vertex pairs,
`lift(u,v)` the face point and `owners[v]` the `(dim, tag)` of the model
entity of each non-box vertex.
"""
function _gmsh_bds_build(p,ids::AbstractVector{<:Integer},
                         interior_faces::AbstractVector{Bool},constrained;
                         lift=nothing,owners=nothing)
    m=_GmshBDS(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],
        Float64[],Int[],Int8[],Int[],Vector{Int32}[],falses(0),Int32[],Int32[],
        Vector{Int32}[],falses(0),Int8[],NTuple{3,Int32}[],falses(0),falses(0),
        Int32[],Int32[],Dict{Int,Int32}(),typemin(Int))
    for v in eachindex(p.x)
        X,Y,Z=lift===nothing || ids[v]<0 ? (0.0,0.0,0.0) : lift(p.x[v],p.y[v])
        owner=owners===nothing || ids[v]<0 ? (2,0) : owners[v]
        _gmsh_bds_add_point!(m,Int(ids[v]),p.x[v],p.y[v],X,Y,Z;
            deg=ids[v]<0 ? 2 : owner[1],tag=ids[v]<0 ? 0 : Int(owner[2]))
    end
    for f in eachindex(p.face_edge)
        h=p.face_edge[f]
        _gmsh_bds_add_triangle!(m,p.origin[h],p.origin[p.next[h]],
            p.origin[p.next[p.next[h]]])
    end
    for (a,b) in constrained
        e=_gmsh_bds_find_edge(m,Int32(a),Int32(b))
        e==0 && throw(ArgumentError("native initial mesh: a recovered model edge is missing"))
        m.edge_deg[e]=1
    end
    for f in m.triangles
        interior_faces[f] ? (m.face_g[f]=true) : _gmsh_bds_del_face!(m,f)
    end
    _gmsh_bds_partition!(m.triangles,m.face_deleted)
    for e in m.edges
        if isempty(m.edge_faces[e])
            _gmsh_bds_del_edge!(m,e)
        else
            m.edge_deg[e]<0 && (m.edge_deg[e]=2)
            for q in (m.edge_p1[e],m.edge_p2[e])
                (m.pdeg[q]<0 || m.pdeg[q]>m.edge_deg[e]) && (m.pdeg[q]=m.edge_deg[e];
                    m.ptag[q]=m.edge_deg[e]==2 ? 1 : m.ptag[q])
            end
        end
    end
    _gmsh_bds_partition!(m.edges,m.edge_deleted)
    for v in 1:min(4,length(p.x))
        m.id[v]<0 && (m.pdeleted[v]=true)
    end
    return m
end

# `delaunayizeBDS`; returns the swaps made, as (p1, p2, op0, op1).
function _gmsh_bds_delaunayize!(m::_GmshBDS,ground::G) where {G}
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
            if _gmsh_bds_swap!(m,e,false)
                count+=1
                push!(swaps,(p1,p2,op0,op1))
            end
        end
        count==0 && break
    end
    return swaps
end

_gmsh_bds_faces(m::_GmshBDS)=NTuple{3,Int32}[_gmsh_bds_nodes(m,f) for f in m.triangles
                                            if !m.face_deleted[f]]

"""
    _gmsh_bds_delaunayize(p, ids, interior_faces, constrained, ground)

Replay `meshGenerator`'s BDS construction from the PolyMesh `p` and run
`delaunayizeBDS` with the first fundamental form `ground(u,v) = (E,F,G)`.
"""
function _gmsh_bds_delaunayize(p,ids::AbstractVector{<:Integer},
                               interior_faces::AbstractVector{Bool},
                               constrained,ground::G;owners=nothing) where {G}
    m=_gmsh_bds_build(p,ids,interior_faces,constrained;owners)
    swaps=_gmsh_bds_delaunayize!(m,ground)
    return _GmshBDSResult(_gmsh_bds_faces(m),swaps)
end

# ── MeshAdapt: refineMeshBDS ────────────────────────────────────────────────

# The face `refineMeshBDS` works on: `lift(u,v)`, `bgm(u,v,x,y,z)` (BGM_MeshSize),
# whether it is a plane (edge lengths are chords) and `extend` (Mesh.
# MeshSizeExtendFromBoundary).
struct _GmshMeshAdaptFace{L,B}
    lift::L
    bgm::B
    plane::Bool
    extend::Bool
end

@inline _gmsh_bds_point_lc(m::_GmshBDS,face,p::Int32)=
    face.extend ? min(m.lc[p],m.lcbgm[p]) : m.lcbgm[p]

@inline function _gmsh_bds_length(m::_GmshBDS,a::Int32,b::Int32)
    dx=m.X[a]-m.X[b];dy=m.Y[a]-m.Y[b];dz=m.Z[a]-m.Z[b]
    return sqrt(dx*dx+dy*dy+dz*dz)
end

function _gmsh_bds_linear_length(m::_GmshBDS,face,a::Int32,b::Int32)
    face.plane && return _gmsh_bds_length(m,a,b)
    g=face.lift(0.5*(m.u[a]+m.u[b]),0.5*(m.v[a]+m.v[b]))
    dx1=m.X[a]-g[1];dy1=m.Y[a]-g[2];dz1=m.Z[a]-g[3]
    dx2=m.X[b]-g[1];dy2=m.Y[b]-g[2];dz2=m.Z[b]-g[3]
    return sqrt(dx1*dx1+dy1*dy1+dz1*dz1)+sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
end

# `NewGetLc(edge)`: the edge length over `correctLC_`.
function _gmsh_bds_edge_size(m::_GmshBDS,face,a::Int32,b::Int32)
    l1=_gmsh_bds_point_lc(m,face,a);l2=_gmsh_bds_point_lc(m,face,b)
    l=0.5*(l1+l2)
    U=0.5*m.u[a]+(1-0.5)*m.u[b];V=0.5*m.v[a]+(1-0.5)*m.v[b]
    g=face.lift(U,V)
    l=min(l,face.bgm(U,V,g[1],g[2],g[3]))
    return _gmsh_bds_linear_length(m,face,a,b)/l
end

function _gmsh_bds_edges_sort!(list::Vector{Tuple{Float64,Int32}},m::_GmshBDS)
    sort!(list;lt=(x,y)->begin
        if abs(x[1]-y[1])<1e-10
            a1=m.id[m.edge_p1[x[2]]];b1=m.id[m.edge_p1[y[2]]]
            a1==b1 && return m.id[m.edge_p2[x[2]]]<m.id[m.edge_p2[y[2]]]
            return a1<b1
        end
        return x[1]<y[1]
    end)
    return list
end

function _gmsh_bds_split_edge!(m::_GmshBDS,e::Int32,mid::Int32)
    p1=m.edge_p1[e];p2=m.edge_p2[e]
    faces=m.edge_faces[e]
    length(faces)==2 || return false
    f0=faces[1];f1=faces[2]
    op0=_gmsh_bds_opposite(m,f0,p1,p2);op1=_gmsh_bds_opposite(m,f1,p1,p2)
    pts1=_gmsh_bds_nodes(m,f0)
    orientation=0
    for i in 1:3
        if pts1[i]==p1
            orientation=pts1[mod1(i+1,3)]==p2 ? 1 : -1
            break
        end
    end
    p1_op1=_gmsh_bds_face_edge(m,f0,p1,op0)
    op1_p2=_gmsh_bds_face_edge(m,f0,op0,p2)
    p1_op2=_gmsh_bds_face_edge(m,f1,p1,op1)
    op2_p2=_gmsh_bds_face_edge(m,f1,op1,p2)
    g1=m.face_g[f0];g2=m.face_g[f1];ge=m.edge_deg[e]
    _gmsh_bds_del_face!(m,f0)
    _gmsh_bds_del_face!(m,m.edge_faces[e][1])
    _gmsh_bds_del_edge!(m,e)
    p1_mid=_gmsh_bds_new_edge!(m,p1,mid)
    mid_p2=_gmsh_bds_new_edge!(m,mid,p2)
    op1_mid=_gmsh_bds_new_edge!(m,op0,mid)
    mid_op2=_gmsh_bds_new_edge!(m,mid,op1)
    if orientation==1
        _gmsh_bds_new_face!(m,op1_mid,p1_op1,p1_mid,g1)
        _gmsh_bds_new_face!(m,mid_op2,op2_p2,mid_p2,g2)
        _gmsh_bds_new_face!(m,op1_p2,op1_mid,mid_p2,g1)
        _gmsh_bds_new_face!(m,p1_op2,mid_op2,p1_mid,g2)
    else
        _gmsh_bds_new_face!(m,p1_op1,op1_mid,p1_mid,g1)
        _gmsh_bds_new_face!(m,op2_p2,mid_op2,mid_p2,g2)
        _gmsh_bds_new_face!(m,op1_mid,op1_p2,mid_p2,g1)
        _gmsh_bds_new_face!(m,mid_op2,p1_op2,p1_mid,g2)
    end
    m.edge_deg[p1_mid]=ge;m.edge_deg[mid_p2]=ge
    m.edge_deg[op1_mid]=g1 ? 2 : -1;m.edge_deg[mid_op2]=g2 ? 2 : -1
    m.pdeg[mid]=ge;m.ptag[mid]=ge==2 ? 1 : 0
    return true
end

function _gmsh_bds_split_pass!(m::_GmshBDS,face,MAXE)
    list=Tuple{Float64,Int32}[]
    for e in m.edges
        (!m.edge_deleted[e] && length(m.edge_faces[e])==2 && m.edge_deg[e]==2) || continue
        lone=_gmsh_bds_edge_size(m,face,m.edge_p1[e],m.edge_p2[e])
        lone>MAXE && push!(list,(-lone,e))
    end
    _gmsh_bds_edges_sort!(list,m)
    mids=zeros(Int32,length(list))
    for (i,(_,e)) in enumerate(list)
        m.edge_deleted[e] && continue
        a=m.edge_p1[e];b=m.edge_p2[e]
        U=0.5*(m.u[a]+m.u[b]);V=0.5*(m.v[a]+m.v[b])
        g=face.lift(U,V)
        mid=_gmsh_bds_add_point!(m,m.maxpoint+1,U,V,g[1],g[2],g[3])
        m.lc[mid]=0.5*(m.lc[a]+m.lc[b])
        m.lcbgm[mid]=face.bgm(U,V,g[1],g[2],g[3])
        mids[i]=mid
    end
    count=0
    for (i,(_,e)) in enumerate(list)
        m.edge_deleted[e] && continue
        mid=mids[i];mid==0 && continue
        if _gmsh_bds_split_edge!(m,e,mid)
            count+=1
        else
            m.pdeleted[mid]=true;delete!(m.by_id,m.id[mid])
        end
    end
    return count
end

function _gmsh_bds_swap_pass!(m::_GmshBDS,face)
    threshold=cos(10.0*pi/180.0)
    count=0
    original=length(m.edges)
    index=1
    while index<=2*original && index<=length(m.edges)
        e=m.edges[index];index+=1
        (!m.edge_deleted[e] && length(m.edge_faces[e])==2) || continue
        p1=m.edge_p1[e];p2=m.edge_p2[e]
        f0,f1=m.edge_faces[e]
        op0=_gmsh_bds_opposite(m,f0,p1,p2);op1=_gmsh_bds_opposite(m,f1,p1,p2)
        qa=min(_gmsh_bds_gamma(m,p1,p2,op0),_gmsh_bds_gamma(m,p1,p2,op1))
        qb=min(_gmsh_bds_gamma(m,p1,op0,op1),_gmsh_bds_gamma(m,p2,op0,op1))
        n1=_gmsh_bds_normal(m,_gmsh_bds_nodes(m,f0)...)
        n2=_gmsh_bds_normal(m,_gmsh_bds_nodes(m,f1)...)
        n1[1]*n2[1]+n1[2]*n2[2]+n1[3]*n2[3]>threshold || continue
        qb>qa || continue
        _gmsh_bds_swap!(m,e,true) && (count+=1)
    end
    _gmsh_bds_cleanup!(m)
    return count
end

function _gmsh_bds_collapse_max_size(m::_GmshBDS,face,e::Int32,p::Int32)
    o=_gmsh_bds_other(m,e,p)
    largest=0.0
    for edge in m.pedges[p]
        # `BDS_Edge(newP1, newP2)` orders its ends by number.
        a=m.edge_p1[edge]==p ? o : m.edge_p1[edge]
        b=m.edge_p2[edge]==p ? o : m.edge_p2[edge]
        m.id[a]<m.id[b] || ((a,b)=(b,a))
        largest=max(largest,_gmsh_bds_edge_size(m,face,a,b))
    end
    return largest
end

function _gmsh_bds_collapse!(m::_GmshBDS,e::Int32,p::Int32)
    length(m.edge_faces[e])==2 || return false
    (m.pdeg[p]==0 || m.pdeg[p]==1) && return false
    m.edge_deg[e]==2 && (m.pdeg[p]!=2 || m.ptag[p]!=1) && return false
    f0,f1=m.edge_faces[e]
    a=m.edge_p1[e];b=m.edge_p2[e]
    o0=_gmsh_bds_opposite(m,f0,a,b);o1=_gmsh_bds_opposite(m,f1,a,b)
    for edge in m.pedges[o0]
        _gmsh_bds_other(m,edge,o0)==o1 && return false
    end
    for q in (o0,o1)
        d=m.pdeg[q];n=length(m.pedges[q])
        d==2 && n<=4 && return false
        0<=d<2 && n<=3 && return false
    end
    faces=_gmsh_bds_point_faces(m,p)
    o=_gmsh_bds_other(m,e,p)
    kept=NTuple{3,Int32}[];kept_g=Bool[]
    area_old=0.0;area_new=0.0
    for f in faces
        pts=_gmsh_bds_nodes(m,f)
        sold=abs(_gmsh_bds_area(m,pts...))
        area_old+=sold
        e in m.face_edges[f] && continue
        moved=map(q->q==p ? o : q,pts)
        snew=abs(_gmsh_bds_area(m,moved...))
        snew<0.02*sold && return false
        area_new+=snew
        push!(kept,moved);push!(kept_g,m.face_g[f])
    end
    abs(area_old-area_new)>1e-12*(area_old+area_new) && return false
    for f in faces
        _gmsh_bds_del_face!(m,f)
    end
    old=Tuple{Int32,Int32,Int8}[]
    for edge in copy(m.pedges[p])
        q1=m.edge_p1[edge]==p ? o : m.edge_p1[edge]
        q2=m.edge_p2[edge]==p ? o : m.edge_p2[edge]
        push!(old,(q1,q2,m.edge_deg[edge]))
        _gmsh_bds_del_edge!(m,edge)
    end
    for (k,nodes) in enumerate(kept)
        f=_gmsh_bds_add_triangle!(m,nodes...)
        m.face_g[f]=kept_g[k]
    end
    for (q1,q2,g) in old
        edge=_gmsh_bds_find_edge(m,q1,q2)
        edge!=0 && m.edge_deg[edge]<0 && (m.edge_deg[edge]=g)
    end
    return true
end

function _gmsh_bds_collapse_pass!(m::_GmshBDS,face,MINE,MAXNP)
    list=Tuple{Float64,Int32}[]
    for e in m.edges
        (!m.edge_deleted[e] && length(m.edge_faces[e])==2 && m.edge_deg[e]==2) || continue
        lone=_gmsh_bds_edge_size(m,face,m.edge_p1[e],m.edge_p2[e])
        lone<MINE && push!(list,(lone,e))
    end
    _gmsh_bds_edges_sort!(list,m)
    count=0
    for (value,e) in list
        m.edge_deleted[e] && continue
        a=m.edge_p1[e];b=m.edge_p2[e]
        lone1=0.0;allow1=false
        if m.id[a]>MAXNP
            lone1=_gmsh_bds_collapse_max_size(m,face,e,a)
            allow1=abs(lone1-1.0)<abs(value-1.0)
        end
        lone2=0.0;allow2=false
        if m.id[b]>MAXNP
            lone2=_gmsh_bds_collapse_max_size(m,face,e,b)
            allow2=abs(lone2-1.0)<abs(value-1.0)
        end
        p=if allow1 && allow2
            abs(lone1-lone2)<1e-12 ? (m.id[a]<m.id[b] ? a : b) :
                (abs(lone1-1.0)<abs(lone2-1.0) ? a : b)
        elseif allow1
            a
        elseif allow2
            b
        else
            Int32(0)
        end
        p!=0 && _gmsh_bds_collapse!(m,e,p) && (count+=1)
    end
    return count
end

# `getOrderedNeighboringVertices`.
function _gmsh_bds_ring(m::_GmshBDS,p::Int32,faces)
    ring=Int32[]
    isempty(faces) && return nothing
    while true
        found=false
        for f in faces
            pts=_gmsh_bds_nodes(m,f)
            pp=pts[1]==p ? (pts[2],pts[3]) : pts[2]==p ? (pts[1],pts[3]) :
                (pts[1],pts[2])
            if isempty(ring)
                push!(ring,pp[1],pp[2]);found=true;break
            else
                q0=ring[end-1];q1=ring[end]
                if q1==pp[1] && q0!=pp[2]
                    push!(ring,pp[2]);found=true;break
                elseif q1==pp[2] && q0!=pp[1]
                    push!(ring,pp[1]);found=true;break
                end
            end
        end
        length(ring)==length(faces) && break
        found || return nothing
    end
    return ring
end

function _gmsh_bds_tutte(m::_GmshBDS,p::Int32,ring)
    isempty(ring) && return 1e22,0.0
    E=0.0;largest=0.0;smallest=0.0
    for (i,q) in enumerate(ring)
        dx=m.X[p]-m.X[q];dy=m.Y[p]-m.Y[q];dz=m.Z[p]-m.Z[q]
        l2=dx*dx+dy*dy+dz*dz
        largest=i>1 ? max(largest,l2) : l2
        smallest=i>1 ? min(smallest,l2) : l2
        E+=l2
    end
    largest==0 && return 1e22,0.0
    return E,smallest/largest
end

function _gmsh_bds_cavity_valid(m::_GmshBDS,p::Int32,ring)
    here=(m.u[p],m.v[p])
    sign=orient2(here,(m.u[ring[1]],m.v[ring[1]]),(m.u[ring[2]],m.v[ring[2]]))
    n=length(ring)
    for i in 2:n
        q0=ring[i];q1=ring[mod1(i+1,n)]
        s=orient2(here,(m.u[q0],m.v[q0]),(m.u[q1],m.v[q1]))
        sign*s<=0 && return false
    end
    return true
end

function _gmsh_bds_kernel(m::_GmshBDS,p::Int32,ring)
    pp=(m.u[p],m.v[p]);ll=m.lc[p]
    kernel=NTuple{2,Float64}[(m.u[q],m.v[q]) for q in ring]
    sizes=Float64[m.lc[q] for q in ring]
    n=length(kernel)
    for i in 1:n
        now=kernel[i];lc_now=sizes[i]
        for j in 1:n
            (i!=j && i!=mod1(j+1,n)) || continue
            q0=kernel[j];q1=kernel[mod1(j+1,n)]
            x=_gmsh_sys2x2(now[1]-pp[1],q0[1]-q1[1],now[2]-pp[2],q0[2]-q1[2],
                           q0[1]-pp[1],q0[2]-pp[2])
            if x[1]>0 && x[1]<1.0
                now=(pp[1]*(1.0-x[1])+now[1]*x[1],pp[2]*(1.0-x[1])+now[2]*x[1])
                lc_now=ll*(1.0-x[1])+lc_now*x[1]
            end
        end
        kernel[i]=now;sizes[i]=lc_now
    end
    return kernel,sizes
end

# `smooth_point_centroid` with `minimizeTutteEnergyParam` and, on a plane,
# `minimizeTutteEnergyProj` through the plane's projection.
function _gmsh_bds_smooth!(m::_GmshBDS,face,p::Int32,threshold,project)
    m.pdeg[p]<=1 && return false
    m.ptag[p]<0 && return true
    faces=_gmsh_bds_point_faces(m,p)
    ring=_gmsh_bds_ring(m,p,faces)
    ring===nothing && return false
    E_unmoved,ratio=_gmsh_bds_tutte(m,p,ring)
    ratio>threshold && return false
    kernel,sizes=_gmsh_bds_kernel(m,p,ring)
    oldX=m.X[p];oldY=m.Y[p];oldZ=m.Z[p];oldU=m.u[p];oldV=m.v[p]
    U=0.0;V=0.0;LC=0.0;sum=0.0
    for (i,k) in enumerate(kernel)
        g=face.lift(k[1],k[2])
        du=m.u[p]-k[1];dv=m.v[p]-k[2]
        denominator=du*du+dv*dv
        if denominator!=0
            dx=m.X[p]-g[1];dy=m.Y[p]-g[2];dz=m.Z[p]-g[3]
            fact=sqrt((dx*dx+dy*dy+dz*dz)/denominator)
            sum+=fact;U+=k[1]*fact;V+=k[2]*fact;LC+=sizes[i]*fact
        end
    end
    if sum!=0
        U/=sum;V/=sum;LC/=sum
    end
    g=face.lift(U,V)
    m.X[p]=g[1];m.Y[p]=g[2];m.Z[p]=g[3]
    E_moved,ratio2=_gmsh_bds_tutte(m,p,ring)
    param=false
    if E_moved<E_unmoved
        m.u[p]=U;m.v[p]=V
        if _gmsh_bds_cavity_valid(m,p,ring)
            m.lc[p]=LC
            param=ratio2>0.25
        else
            m.X[p]=oldX;m.Y[p]=oldY;m.Z[p]=oldZ;m.u[p]=oldU;m.v[p]=oldV
        end
    else
        m.X[p]=oldX;m.Y[p]=oldY;m.Z[p]=oldZ
    end
    param && return true
    project===nothing && return false
    # minimizeTutteEnergyProj, from wherever the Param step left the point.
    _gmsh_bds_project!(m,p,ring,E_unmoved,project) || return false
    E_unmoved,_=_gmsh_bds_tutte(m,p,ring)
    _gmsh_bds_project!(m,p,ring,E_unmoved,project)
    return true
end

function _gmsh_bds_project!(m::_GmshBDS,p::Int32,ring,E_unmoved,project)
    oldX=m.X[p];oldY=m.Y[p];oldZ=m.Z[p];oldU=m.u[p];oldV=m.v[p]
    x=(0.0,0.0,0.0);total=0.0
    n=length(ring)
    for i in 1:n
        a=ring[i];b=ring[mod1(i+1,n)]
        v1=(m.X[a]-oldX,m.Y[a]-oldY,m.Z[a]-oldZ)
        v2=(m.X[b]-oldX,m.Y[b]-oldY,m.Z[b]-oldZ)
        c=(v1[2]*v2[3]-v2[2]*v1[3],-(v1[1]*v2[3]-v2[1]*v1[3]),v1[1]*v2[2]-v2[1]*v1[2])
        nrm=sqrt(c[1]*c[1]+c[2]*c[2]+c[3]*c[3])
        w=nrm/3.0
        x=(x[1]+(m.X[a]+oldX+m.X[b])*w,x[2]+(m.Y[a]+oldY+m.Y[b])*w,
           x[3]+(m.Z[a]+oldZ+m.Z[b])*w)
        total+=nrm
    end
    x=(x[1]/total,x[2]/total,x[3]/total)
    uv=project(x)
    m.u[p]=uv[1];m.v[p]=uv[2]
    if _gmsh_bds_cavity_valid(m,p,ring)
        m.X[p]=x[1];m.Y[p]=x[2];m.Z[p]=x[3]
        E_moved,_=_gmsh_bds_tutte(m,p,ring)
        E_moved<E_unmoved && return true
    end
    m.X[p]=oldX;m.Y[p]=oldY;m.Z[p]=oldZ;m.u[p]=oldU;m.v[p]=oldV
    return false
end

function _gmsh_bds_smooth_pass!(m::_GmshBDS,face,project)
    for p in sort!(Int32[q for q in eachindex(m.id) if !m.pdeleted[q]];by=q->m.id[q])
        _gmsh_bds_smooth!(m,face,p,0.5,project)
    end
    return nothing
end

# `computeNodalSizes`: a model point's size is its longest model edge.
function _gmsh_bds_nodal_sizes!(m::_GmshBDS)
    for p in eachindex(m.id)
        m.pdeleted[p] && continue
        L=0.0;ne=0
        for e in m.pedges[p]
            if m.edge_deg[e]==1
                l=_gmsh_bds_length(m,m.edge_p1[e],m.edge_p2[e])
                L=ne>0 ? max(L,l) : l
                ne+=1
            end
        end
        if m.pdeg[p]>=0 && m.ptag[p]>0
            ne==0 && (L=_GMSH_BDS_MAX_LC)
            m.lc[p]=L;m.lcbgm[p]=L
        end
    end
    return nothing
end

function _gmsh_bds_refine!(m::_GmshBDS,face,NIT::Int,nodal::Bool,project)
    MAXNP=m.maxpoint
    nodal && _gmsh_bds_nodal_sizes!(m)
    MINE=0.7;MAXE=1.4
    IT=0
    while true
        nb_split=_gmsh_bds_split_pass!(m,face,MAXE)
        IT==0 && (nb_split+=_gmsh_bds_split_pass!(m,face,MAXE))
        _gmsh_bds_smooth_pass!(m,face,project)
        _gmsh_bds_swap_pass!(m,face)
        _gmsh_bds_smooth_pass!(m,face,project)
        nb_collapse=_gmsh_bds_collapse_pass!(m,face,MINE,MAXNP)
        _gmsh_bds_smooth_pass!(m,face,project)
        _gmsh_bds_swap_pass!(m,face)
        _gmsh_bds_smooth_pass!(m,face,project)
        if IT==abs(NIT)
            nb_collapse+=_gmsh_bds_collapse_pass!(m,face,0.45,MAXNP)
            _gmsh_bds_smooth_pass!(m,face,project)
        end
        _gmsh_bds_cleanup!(m)
        minL=1e22;maxL=0.0
        for e in m.edges
            m.edge_deleted[e] && continue
            a=m.edge_p1[e];b=m.edge_p2[e]
            lone=2*_gmsh_bds_length(m,a,b)/
                (_gmsh_bds_point_lc(m,face,a)+_gmsh_bds_point_lc(m,face,b))
            maxL=max(maxL,lone);minL=min(minL,lone)
        end
        ((minL>MINE && maxL<MAXE) || IT>abs(NIT)) && break
        IT+=1
        (nb_split==0 && nb_collapse==0) && break
    end
    return nothing
end
