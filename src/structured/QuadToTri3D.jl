# A port of Gmsh 4.15.2's general QuadToTri region mesher
# (contrib/QuadTri/QuadTriExtruded3D.cpp, QuadTriUtils.cpp, the top-surface
# part of QuadTriExtruded2D.cpp), by T. S. Strickler.
#
# Gmsh decides diagonals by comparing MVertex pointers. Here a vertex is an
# integer whose order is Gmsh's creation order of the mesh vertices, which is
# the order those pointers follow when Gmsh allocates them one after another.
# Indices into element vertex lists are kept 0-based as in the C++ (`V(i)` reads
# entry i), so each branch can be checked against the original line by line.
module QuadToTri3D


const Edge=Tuple{Int,Int}

@inline _ekey(a::Int,b::Int)=a<b ? (a,b) : (b,a)
edge_exists(a,b,edges::Set{Edge})=(_ekey(a,b) in edges)
create_edge!(edges::Set{Edge},a,b)=(push!(edges,_ekey(a,b));nothing)
delete_edge!(edges::Set{Edge},a,b)=(delete!(edges,_ekey(a,b));nothing)

# A forbidden quadrangle is stored by its diagonal from the lowest vertex.
function _forbidden_key(v)
    low=0
    for p in 1:3
        v[p+1]<v[low+1] && (low=p)
    end
    return (v[low+1],v[(low+2)%4+1])
end
create_forbidden!(forbidden::Set{Edge},v)=(push!(forbidden,_forbidden_key(v));nothing)
forbidden_exists(v,forbidden::Set{Edge})=(_forbidden_key(v) in forbidden)

# Element identity: (0,i) for the i-th source triangle, (1,i) for a quadrangle
# (0-based, Gmsh's face->triangles / face->quadrangles order).
const ElemId=Tuple{Int,Int}

# The swept region as Gmsh's QuadToTri code sees it (single region, its
# laterals bounded by no other QuadToTri region).
#  - `tris`, `quads`: source element vertices (Gmsh element vertex order).
#  - `layers`: NbElmLayer per layer group; `vid(v,level)` the mesh vertex of
#    source vertex v at global level (0 = source, sum(layers) = top); a vertex
#    on a rotation axis maps to itself at every level. `coord(id)` its point.
#  - `boundary`: source vertices on the source face's bounding curves.
#  - `rotate`: a Rotate or Translate+Rotate extrusion (ep->geo.Type).
#  - `lateral_lines[s]`: the source curve segments (v0,v1) of lateral s;
#    `lateral_flag[s]`: createQuaTri's tri_quad_flag (0, 1 or 2) and
#    `lateral_recombine[s]` the lateral's own Recombine.
#  - `addverts`: QuadTriAddVerts; `structured_root`: the source is
#    transfinite or extruded (MeshQuadToTriTopSurface for added vertices).
#  - `next_id`: the first id free for body-centred vertices.
mutable struct QTRegion
    tris::Vector{NTuple{3,Int}}
    quads::Vector{NTuple{4,Int}}
    layers::Vector{Int}
    vid::Function
    coord::Function
    boundary::Set{Int}
    rotate::Bool
    lateral_lines::Vector{Vector{Edge}}
    lateral_flag::Vector{Int}
    lateral_recombine::Vector{Bool}
    addverts::Bool
    structured_root::Bool
    next_id::Int
end

_level(r::QTRegion,j::Int,k::Int)=sum(r.layers[1:j];init=0)+k

# getExtrudedVertices: the element's vertices at (j,k), then at (j,k+1).
function extruded_vertices(r::QTRegion,elem::ElemId,j::Int,k::Int)
    src=elem[1]==0 ? r.tris[elem[2]+1] : r.quads[elem[2]+1]
    l=_level(r,j,k)
    return Int[(r.vid(v,l) for v in src)...,(r.vid(v,l+1) for v in src)...]
end

function element_vertices(r::QTRegion,elem::ElemId)
    return elem[1]==0 ? collect(r.tris[elem[2]+1]) : collect(r.quads[elem[2]+1])
end

mutable struct QTState
    quadToTri_edges::Set{Edge}
    forbidden_edges::Set{Edge}
    lat_tri_diags::Set{Edge}
    problems::Dict{ElemId,Set{Tuple{Int,Int}}}
end
QTState()=QTState(Set{Edge}(),Set{Edge}(),Set{Edge}(),Dict{ElemId,Set{Tuple{Int,Int}}}())

add_problem!(problems,elem,j,k)=(push!(get!(()->Set{Tuple{Int,Int}}(),problems,elem),(j,k));nothing)

# fill_touch_bnd
function fill_touch_bnd(vert_bnd,n_lat)
    touch=zeros(Int,n_lat)
    for i in 0:n_lat-1
        if vert_bnd[i+1]
            touch[i+1]=1
            touch[(i+n_lat-1)%n_lat+1]=1
        end
    end
    return touch
end

# Face types of one swept element (getFaceTypes): faces 0..n_lat-1 are the
# laterals starting at vertex 0, n_lat the bottom, n_lat+1 the top.
mutable struct FaceTypes
    degen::Vector{Int}
    single_tri::Vector{Int}
    recomb::Vector{Int}
    fixed_diag::Vector{Int}
    adj_diag::Vector{Int}
    free::Vector{Int}
end
FaceTypes()=FaceTypes(Int[],Int[],Int[],Int[],Int[],Int[])

function get_face_types(r::QTRegion,j,k,verts,st::QTState,vert_bnd)
    ft=FaceTypes()
    j_top_start=length(r.layers)-1
    k_top_start=r.layers[end]-1
    n_lat=length(verts)==6 ? 3 : 4
    V=i->verts[i+1]
    nfix1=[-p*p-p-1 for p in 0:n_lat+1];nfix2=[-p*p-p-2 for p in 0:n_lat+1]
    nadj1=copy(nfix1);nadj2=copy(nfix2)
    free_flag=zeros(Int,n_lat+2)
    touch_bnd=fill_touch_bnd(vert_bnd,n_lat)
    E=st.quadToTri_edges;F=st.forbidden_edges;L=st.lat_tri_diags
    for p in 0:n_lat-1
        p2=(p+1)%n_lat
        v_face=(V(p),V(p2),V(p2+n_lat),V(p+n_lat))
        if V(p)==V(n_lat+p) && V(p2)==V(n_lat+p2)
            push!(ft.degen,p)
        elseif (V(p)==V(n_lat+p) && V(p2)!=V(n_lat+p2)) ||
               (V(p)!=V(n_lat+p) && V(p2)==V(n_lat+p2))
            push!(ft.single_tri,p)
        elseif forbidden_exists(v_face,F)
            push!(ft.recomb,p)
        elseif edge_exists(V(p),V(p2+n_lat),E)
            nfix1[p+1]=p;nfix2[p+1]=p2+n_lat;push!(ft.fixed_diag,p)
        elseif edge_exists(V(n_lat+p),V(p2),E)
            nfix1[p+1]=p+n_lat;nfix2[p+1]=p2;push!(ft.fixed_diag,p)
        elseif edge_exists(V(p),V(p2+n_lat),L)
            nadj1[p+1]=p;nadj2[p+1]=p2+n_lat;push!(ft.adj_diag,p)
        elseif edge_exists(V(n_lat+p),V(p2),L)
            nadj1[p+1]=p+n_lat;nadj2[p+1]=p2;push!(ft.adj_diag,p)
        elseif touch_bnd[p+1]==0 && (j<j_top_start || (j==j_top_start && k<k_top_start))
            push!(ft.recomb,p)
        else
            push!(ft.free,p);free_flag[p+1]=1
        end
    end
    if n_lat==3
        push!(ft.single_tri,3)
    else
        v_bot=(V(0),V(1),V(2),V(3))
        if (j==0 && k==0) || forbidden_exists(v_bot,F)
            push!(ft.recomb,4)
        elseif edge_exists(V(0),V(2),E)
            nfix1[5]=0;nfix2[5]=2;push!(ft.fixed_diag,4)
        elseif edge_exists(V(1),V(3),E)
            nfix1[5]=1;nfix2[5]=3;push!(ft.fixed_diag,4)
        elseif edge_exists(V(0),V(2),L)
            nadj1[5]=0;nadj2[5]=2;push!(ft.adj_diag,4)
        elseif edge_exists(V(1),V(3),L)
            nadj1[5]=1;nadj2[5]=3;push!(ft.adj_diag,4)
        else
            push!(ft.free,4);free_flag[5]=1
        end
    end
    if n_lat==3
        push!(ft.single_tri,4)
    else
        v_top=(V(4),V(5),V(6),V(7))
        if forbidden_exists(v_top,F)
            push!(ft.recomb,5)
        elseif edge_exists(V(4),V(6),E)
            nfix1[6]=4;nfix2[6]=6;push!(ft.fixed_diag,5)
        elseif edge_exists(V(5),V(7),E)
            nfix1[6]=5;nfix2[6]=7;push!(ft.fixed_diag,5)
        elseif edge_exists(V(4),V(6),L)
            nadj1[6]=4;nadj2[6]=6;push!(ft.adj_diag,5)
        elseif edge_exists(V(5),V(7),L)
            nadj1[6]=5;nadj2[6]=7;push!(ft.adj_diag,5)
        else
            push!(ft.free,5);free_flag[6]=1
        end
    end
    return ft,nfix1,nfix2,nadj1,nadj2,free_flag
end

# bruteForceEdgeQuadToTriPrism
function brute_force_prism!(st::QTState,elem,j,k,verts,ft::FaceTypes,
        edges_new,forbidden_new,problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
    V=i->verts[i+1]
    E=st.quadToTri_edges
    num_degen=length(ft.degen);num_single_tri=length(ft.single_tri)
    num_recomb=length(ft.recomb);num_fixed_diag=length(ft.fixed_diag)
    num_adj_diag=length(ft.adj_diag);num_free=length(ft.free)
    for ind in ft.recomb
        v=(V(ind),V((ind+1)%3),V((ind+1)%3+3),V(ind+3))
        create_forbidden!(forbidden_new,v);create_forbidden!(st.forbidden_edges,v)
    end
    (3 in ft.single_tri && 4 in ft.single_tri) ||
        error("bruteForceEdgeQuadToTriPrism(): invalid face code for top and/or bottom")
    if num_single_tri==4 && num_degen==1
        return
    elseif num_single_tri==4
        if num_adj_diag>0
            ind=ft.adj_diag[1]
            create_edge!(edges_new,V(nadj1[ind+1]),V(nadj2[ind+1]))
            create_edge!(E,V(nadj1[ind+1]),V(nadj2[ind+1]))
        end
        return
    end
    if num_recomb==2
        if num_free>0
            for ind in ft.free
                v=(V(ind),V((ind+1)%3),V((ind+1)%3+3),V(ind+3))
                create_forbidden!(forbidden_new,v);create_forbidden!(st.forbidden_edges,v)
            end
        else
            if num_adj_diag>0
                ind=ft.adj_diag[1]
                create_edge!(E,V(nadj1[ind+1]),V(nadj2[ind+1]))
                create_edge!(edges_new,V(nadj1[ind+1]),V(nadj2[ind+1]))
            end
            add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
        end
        return
    end
    if num_recomb==1 && num_free==2
        p_recomb=ft.recomb[1]
        p_free=(p_recomb+2)%3
        ind_low=V(p_free)<V(p_free+3) ? p_free : p_free+3
        add=ind_low<3 ? 0 : 3
        for d in (1,2)
            other=(ind_low-add+d)%3+3-add
            create_edge!(E,V(ind_low),V(other));create_edge!(edges_new,V(ind_low),V(other))
        end
        return
    end
    if num_free>=3
        ind_low=0
        for p in 1:5
            V(ind_low)>V(p) && (ind_low=p)
        end
        add=ind_low<3 ? 0 : 3
        for d in (1,2)
            other=(ind_low-add+d)%3+3-add
            create_edge!(E,V(ind_low),V(other));create_edge!(edges_new,V(ind_low),V(other))
        end
        return
    end
    valid_division=false
    face_done=[0,0]
    for t in 0:3
        n1=[-p*p-p-1 for p in 0:2];n2=[-p*p-p-2 for p in 0:2]
        face_is_free=falses(3)
        for p in 0:2
            if nfix1[p+1]>=0
                n1[p+1]=nfix1[p+1];n2[p+1]=nfix2[p+1]
            elseif t==0 && nadj1[p+1]>=0
                n1[p+1]=nadj1[p+1];n2[p+1]=nadj2[p+1]
            elseif (t>=1 && t<3 && nadj1[p+1]>=0) || (t==2 && free_flag[p+1]!=0)
                if (V(p)<V((p+1)%3) && V(p)<V(p+3)) ||
                   (V((p+1)%3+3)<V((p+1)%3) && V((p+1)%3+3)<V(p+3))
                    n1[p+1]=p;n2[p+1]=(p+1)%3+3
                else
                    n1[p+1]=p+3;n2[p+1]=(p+1)%3
                end
            elseif (t==3 && nadj1[p+1]>=0) || free_flag[p+1]!=0
                face_is_free[p+1]=true
            end
        end
        for p in 0:2
            q=(p+1)%3
            if face_is_free[p+1] && face_is_free[q+1]
                face_is_free[p+1]=false
                if V(q)<V(q+3)
                    n1[p+1]=p+3;n2[p+1]=q
                else
                    n1[p+1]=p;n2[p+1]=q+3
                end
            end
            if face_is_free[p+1] && n1[q+1]>=0
                valid_division=true;face_done[1]=p;face_done[2]=q
                if n1[q+1]==q || n2[q+1]==q
                    n1[p+1]=p+3;n2[p+1]=q
                else
                    n1[p+1]=p;n2[p+1]=q+3
                end
            elseif face_is_free[q+1] && n1[p+1]>=0
                valid_division=true;face_done[1]=p;face_done[2]=q
                if n2[p+1]==q || n1[p+1]==q
                    n1[q+1]=q;n2[q+1]=(p+2)%3+3
                else
                    n1[q+1]=q+3;n2[q+1]=(p+2)%3
                end
            elseif n2[p+1]==n1[q+1] || n2[p+1]==n2[q+1] || n1[p+1]==n1[q+1] || n1[p+1]==n2[q+1]
                valid_division=true;face_done[1]=p;face_done[2]=q
            end
            if valid_division
                create_edge!(E,V(n1[p+1]),V(n2[p+1]));create_edge!(edges_new,V(n1[p+1]),V(n2[p+1]))
                create_edge!(E,V(n1[q+1]),V(n2[q+1]));create_edge!(edges_new,V(n1[q+1]),V(n2[q+1]))
                break
            end
        end
        (valid_division || num_fixed_diag+num_recomb==3) && break
    end
    for ind in ft.adj_diag
        if ind!=face_done[1] && ind!=face_done[2]
            create_edge!(edges_new,V(nadj1[ind+1]),V(nadj2[ind+1]))
            create_edge!(E,V(nadj1[ind+1]),V(nadj2[ind+1]))
        end
    end
    if !valid_division
        add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
    end
    return
end

include("QuadToTri3DDegenerate.jl")
include("QuadToTri3DHexa.jl")
include("QuadToTri3DEdges.jl")
include("QuadToTri3DElements.jl")
include("QuadToTri3DHexElements.jl")
include("QuadToTri3DDriver.jl")

end
