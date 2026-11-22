# QuadToTriCreateElements and the element builders it calls
# (QuadTriExtruded3D.cpp, QtMakeCentroidVertex / QtFindVertsCentroid of
# QuadTriUtils.cpp).

mutable struct QTMesh
    tets::Vector{NTuple{4,Int}}
    pyramids::Vector{NTuple{5,Int}}
    prisms::Vector{NTuple{6,Int}}
    hexes::Vector{NTuple{8,Int}}
    # Body-centred vertices, in creation order: id => coordinates.
    new_ids::Vector{Int}
    new_coords::Vector{NTuple{3,Float64}}
    lookup::Dict{NTuple{3,Float64},Int}
    next_id::Int
    coord::Function
end

QTMesh(coord,next_id)=QTMesh(NTuple{4,Int}[],NTuple{5,Int}[],NTuple{6,Int}[],NTuple{8,Int}[],
    Int[],NTuple{3,Float64}[],Dict{NTuple{3,Float64},Int}(),next_id,coord)

# QtFindVertsCentroid: mean of the distinct swept vertices.
function verts_centroid(out::QTMesh,v)
    n_lat=length(v)==6 ? 3 : length(v)==8 ? 4 : 0
    x=0.0;y=0.0;z=0.0;count=0
    for p in 0:length(v)-1
        if (length(v)==6 || length(v)==8) && p>n_lat-1 && v[p+1]==v[p-n_lat+1]
            continue
        end
        count+=1
        c=out.coord(v[p+1])
        x+=c[1];y+=c[2];z+=c[3]
    end
    return (x/count,y/count,z/count)
end

function make_centroid_vertex!(out::QTMesh,v)
    c=verts_centroid(out,v)
    id=get(out.lookup,c,0)
    id!=0 && return id
    id=out.next_id;out.next_id+=1
    out.lookup[c]=id;push!(out.new_ids,id);push!(out.new_coords,c)
    return id
end

add_tet!(out,a,b,c,d)=push!(out.tets,(a,b,c,d))
add_pyr!(out,a,b,c,d,e)=push!(out.pyramids,(a,b,c,d,e))
add_pri!(out,a,b,c,d,e,f)=push!(out.prisms,(a,b,c,d,e,f))
add_hex!(out,a,b,c,d,e,f,g,h)=push!(out.hexes,(a,b,c,d,e,f,g,h))

function mesh_with_internal_vertex!(out::QTMesh,v,n1,n2)
    n_lat=length(v)==6 ? 3 : 4
    V=i->v[i+1]
    c=verts_centroid(out,v)
    v_int=get(out.lookup,c,0)
    v_int==0 && return   # "Could not find centroid vertex"
    for p in 0:n_lat-1
        p2=(p+1)%n_lat
        if V(p)==V(p+n_lat) && V(p2)==V(p2+n_lat)
            continue
        elseif V(p)==V(p+n_lat) || V(p2)==V(p2+n_lat)
            v_dup=V(p)==V(p+n_lat) ? V(p) : V(p2)
            v_non_dup=v_dup==V(p) ? V(p2) : V(p)
            v_non_dup2=v_non_dup==V(p) ? V(p+n_lat) : V(p2+n_lat)
            add_tet!(out,v_dup,v_non_dup,v_non_dup2,v_int)
        elseif n1[p+1]==p || n2[p+1]==p
            add_tet!(out,V(p),V(p2),V(p2+n_lat),v_int)
            add_tet!(out,V(p),V(p2+n_lat),V(p+n_lat),v_int)
        elseif n1[p+1]==p+n_lat || n2[p+1]==p+n_lat
            add_tet!(out,V(p),V(p2),V(p+n_lat),v_int)
            add_tet!(out,V(p2),V(p2+n_lat),V(p+n_lat),v_int)
        else
            add_pyr!(out,V(p),V(p2),V(p2+n_lat),V(p+n_lat),v_int)
        end
    end
    if n_lat==3
        add_tet!(out,V(0),V(1),V(2),v_int)
        add_tet!(out,V(3),V(5),V(4),v_int)
    else
        for p in 4:5
            add=p==4 ? 0 : 4
            if n1[p+1]==0+add || n2[p+1]==0+add
                add_tet!(out,V(add),V(1+add),V(2+add),v_int)
                add_tet!(out,V(add),V(2+add),V(3+add),v_int)
            elseif n1[p+1]==1+add || n2[p+1]==1+add
                add_tet!(out,V(1+add),V(2+add),V(3+add),v_int)
                add_tet!(out,V(1+add),V(3+add),V(add),v_int)
            else
                add_pyr!(out,V(add),V(1+add),V(2+add),V(3+add),v_int)
            end
        end
    end
    return
end

_is_problem(problems,source,j,k)=haskey(problems,source) && (j,k) in problems[source]

function pri_pyr_tet!(out::QTMesh,st::QTState,v,j,k,source,problems_new,bnd_elem,is_addverts)
    V=i->v[i+1]
    dup=Int[]
    for i in 0:2
        V(i)==V(i+3) && push!(dup,i)
    end
    m=length(dup)
    is_problem=!is_addverts && _is_problem(st.problems,source,j,k)
    (m>0 && is_problem) && return
    m>2 && return
    n1=fill(-1,3);n2=fill(-2,3)
    found_diags=false
    E=st.quadToTri_edges
    for p in 0:2
        n1[p+1]=-p*p-p-1;n2[p+1]=-p*p-p-2
        if V(p)==V(p+3) || V((p+1)%3)==V((p+1)%3+3)
            continue
        elseif edge_exists(V(p),V((p+1)%3+3),E)
            n1[p+1]=p;n2[p+1]=(p+1)%3+3;found_diags=true
        elseif edge_exists(V(p+3),V((p+1)%3),E)
            n1[p+1]=p+3;n2[p+1]=(p+1)%3;found_diags=true
        end
    end
    if is_addverts && bnd_elem && found_diags
        mesh_with_internal_vertex!(out,v,n1,n2)
        return
    end
    if m==2 && !is_problem
        dup_face=(dup[1]+1)%3==dup[2] ? dup[1] : dup[2]
        add_tet!(out,V(dup_face),V((dup_face+1)%3),V((dup_face+2)%3),V((dup_face+2)%3+3))
        return
    elseif m==1 && !is_problem
        p=(dup[1]+1)%3
        if n1[p+1]>=0
            add=n1[p+1]<3 ? 3 : -3
            add_tet!(out,V(n1[p+1]),V(n2[p+1]),V(n1[p+1]+add),V(dup[1]))
            add_tet!(out,V(n1[p+1]),V(n2[p+1]-add),V(n2[p+1]),V(dup[1]))
        else
            add_pyr!(out,V(p),V(p+3),V((p+1)%3+3),V((p+1)%3),V(dup[1]))
        end
        return
    elseif m==0 && !is_problem
        if !found_diags
            add_pri!(out,V(0),V(1),V(2),V(3),V(4),V(5))
            return
        else
            for p in 0:2
                if n1[p+1]>=0 && n2[p+1]==n1[(p+1)%3+1]
                    add=n2[p+1]<3 ? 3 : -3
                    tet_top=n2[p+1]+add;pyr_top=n2[p+1];base=(p+2)%3
                    add_tet!(out,V(n1[p+1]),V(n2[p+1]),V(n2[(p+1)%3+1]),V(tet_top))
                    if n1[base+1]>=0
                        add_2=n1[base+1]<3 ? 3 : -3
                        add_tet!(out,V(n1[base+1]),V(n2[base+1]),V(n1[base+1]+add_2),V(pyr_top))
                        add_tet!(out,V(n1[base+1]),V(n2[base+1]),V(n2[base+1]-add_2),V(pyr_top))
                    else
                        add_pyr!(out,V(base),V(base+3),V((base+1)%3+3),V((base+1)%3),V(pyr_top))
                    end
                    return
                end
            end
        end
    end
    if !is_problem
        make_centroid_vertex!(out,v)
        add_problem!(problems_new,source,j,k)
        is_problem=true
    end
    mesh_with_internal_vertex!(out,v,n1,n2)
    return
end
