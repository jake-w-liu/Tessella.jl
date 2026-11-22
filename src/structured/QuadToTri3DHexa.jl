# addEdgesForQuadToTriFullHexa and bruteForceEdgeQuadToTriHexa
# (QuadTriExtruded3D.cpp).

# The lowest-pointer diagonal of face p (lateral p<4, bottom 4, top 5).
function _hexa_lowest_diag(V,p)
    if p<4
        if (V(p)<V(p+4) && V(p)<V((p+1)%4)) ||
           (V((p+1)%4+4)<V(p+4) && V((p+1)%4+4)<V((p+1)%4))
            return p,(p+1)%4+4
        end
        return p+4,(p+1)%4
    end
    add=p==4 ? 0 : 4
    if (V(0+add)<V(1+add) && V(0+add)<V(3+add)) || (V(2+add)<V(1+add) && V(2+add)<V(3+add))
        return 0+add,2+add
    end
    return 1+add,3+add
end

# Opposite faces p1,p2 with definite diagonals that cannot both cut the hexa.
function _hexa_pair_incompatible(n1,n2,p1,p2)
    if p1<4
        return ((n1[p1+1]==p1 || n2[p1+1]==p1) && n1[p2+1]!=p2+4 && n2[p2+1]!=p2+4) ||
               ((n1[p1+1]==p1+4 || n2[p1+1]==p1+4) && n1[p2+1]!=p2 && n2[p2+1]!=p2)
    end
    add=p1==4 ? 4 : -4
    return n1[p1+1]+add!=n1[p2+1] && n2[p1+1]+add!=n1[p2+1]
end

# Fix the free members of the opposite pair (p1,p2) from each other.
function _hexa_pair_complete!(V,n1,n2,face_is_free,p1,p2)
    if face_is_free[p1+1] && face_is_free[p2+1]
        face_is_free[p1+1]=false
        n1[p1+1],n2[p1+1]=_hexa_lowest_diag(V,p1)
    end
    if n1[p1+1]>=0 && face_is_free[p2+1]
        if p1<4
            if n1[p1+1]==p1 || n2[p1+1]==p1
                n1[p2+1]=p2+4;n2[p2+1]=(p2+1)%4
            else
                n1[p2+1]=p2;n2[p2+1]=(p2+1)%4+4
            end
        else
            add=p1==4 ? 4 : -4
            n1[p2+1]=n1[p1+1]+add;n2[p2+1]=n2[p1+1]+add
        end
        face_is_free[p2+1]=false
    elseif n1[p2+1]>=0 && face_is_free[p1+1]
        if p1<4
            if n1[p2+1]==p2 || n2[p2+1]==p2
                n1[p1+1]=p1+4;n2[p1+1]=(p1+1)%4
            else
                n1[p1+1]=p1;n2[p1+1]=(p1+1)%4+4
            end
        else
            add=p2==4 ? 4 : -4
            n1[p1+1]=n1[p2+1]+add;n2[p1+1]=n2[p2+1]+add
        end
        face_is_free[p1+1]=false
    end
    return nothing
end

function _hexa_face_vertices(V,p)
    p<4 && return (V(p),V((p+1)%4),V((p+1)%4+4),V(p+4))
    add=p==4 ? 0 : 4
    return (V(add),V(add+1),V(add+2),V(add+3))
end

function full_hexa_edges!(st::QTState,elem,j,k,verts,ft::FaceTypes,
        edges_new,forbidden_new,problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
    V=i->verts[i+1]
    E=st.quadToTri_edges
    edge!(a,b)=(create_edge!(E,V(a),V(b));create_edge!(edges_new,V(a),V(b)))
    num_fixed_diag=length(ft.fixed_diag);num_adj_diag=length(ft.adj_diag)
    num_free=length(ft.free)
    if num_free>=6
        ind_low=0
        for p in 1:7
            V(p)<V(ind_low) && (ind_low=p)
        end
        add=ind_low<4 ? 0 : 4
        edge!(ind_low,(ind_low-add+2)%4+add)
        edge!(ind_low,(ind_low-add+1)%4+4-add)
        edge!(ind_low,(ind_low-add+3)%4+4-add)
        return
    end
    face_done=[-1,-2,-3,-4]
    p1_hold=-1;p2_hold=-2
    n1_hold=[0,0];n2_hold=[0,0]
    valid_division=false
    for t in 0:3
        n1=[-p*p-p-1 for p in 0:5];n2=[-p*p-p-2 for p in 0:5]
        face_is_free=falses(6)
        for p in 0:5
            if nfix1[p+1]>=0
                n1[p+1]=nfix1[p+1];n2[p+1]=nfix2[p+1]
            elseif t==0 && nadj1[p+1]>=0
                n1[p+1]=nadj1[p+1];n2[p+1]=nadj2[p+1]
            elseif (t>=1 && t<3 && nadj1[p+1]>=0) || (t==2 && free_flag[p+1]!=0)
                if p<4
                    if (V(p)<V((p+1)%4) && V(p)<V(p+4)) ||
                       (V((p+1)%4+4)<V((p+1)%4) && V((p+1)%4+4)<V(p+4))
                        n1[p+1]=p;n2[p+1]=(p+1)%4+4
                    else
                        n1[p+1]=p+4;n2[p+1]=(p+1)%4
                    end
                else
                    add=p==4 ? 0 : 4
                    if (V(0+add)<V(1+add) && V(0+add)<V(3+add)) ||
                       (V(2+add)<V(1+add) && V(2+add)<V(3+add))
                        n1[p+1]=0+add;n2[p+1]=2+add
                    else
                        n1[p+1]=1+add;n2[p+1]=3+add
                    end
                end
            elseif t==3 && (nadj1[p+1]>=0 || free_flag[p+1]!=0)
                face_is_free[p+1]=true
            end
        end
        if !valid_division
            for p in 0:2
                p1,p2=p>1 ? (4,5) : (p,(p+2)%4)
                ((n1[p1+1]<0 && !face_is_free[p1+1]) || (n1[p2+1]<0 && !face_is_free[p2+1])) && continue
                n1[p1+1]>=0 && n1[p2+1]>=0 && _hexa_pair_incompatible(n1,n2,p1,p2) && continue
                _hexa_pair_complete!(V,n1,n2,face_is_free,p1,p2)
                if p1_hold<0 || p2_hold<0
                    required=0
                    (nfix1[p1+1]>=0 || nadj1[p1+1]>=0) && (required+=1)
                    (nfix1[p2+1]>=0 || nadj1[p2+1]>=0) && (required+=1)
                    if p1_hold<0 && num_fixed_diag+num_adj_diag<=required
                        p1_hold=p1;p2_hold=p2
                        n1_hold[1]=n1[p1+1];n2_hold[1]=n2[p1+1]
                        n1_hold[2]=n1[p2+1];n2_hold[2]=n2[p2+1]
                    end
                end
                for s in 0:2
                    (s==p1 || (s>1 && (p1==4 || p1==5))) && continue
                    s1,s2=s>1 ? (4,5) : (s,(s+2)%4)
                    ((n1[s1+1]<0 && !face_is_free[s1+1]) || (n1[s2+1]<0 && !face_is_free[s2+1])) && continue
                    n1[s1+1]>=0 && n1[s2+1]>=0 && _hexa_pair_incompatible(n1,n2,s1,s2) && continue
                    valid_division=true
                    face_done[1]=p1;face_done[2]=p2;face_done[3]=s1;face_done[4]=s2
                    _hexa_pair_complete!(V,n1,n2,face_is_free,s1,s2)
                    valid_division && break
                end
                if !valid_division
                    for s in 0:7
                        add=s<4 ? 0 : 4
                        faces=[-1,-1]
                        third_face=s<4 ? 4 : 5
                        count=0
                        if s-add!=p1 && s-add!=p2 &&
                           (n1[s-add+1]==s || n2[s-add+1]==s || face_is_free[s-add+1])
                            faces[count+1]=s-add;count+=1
                        end
                        q=(s-add+3)%4
                        if q!=p1 && q!=p2 && (n1[q+1]==s || n2[q+1]==s || face_is_free[q+1])
                            faces[count+1]=q;count+=1
                        end
                        if count<2 && third_face!=p1 && third_face!=p2 &&
                           (n1[third_face+1]==s || n2[third_face+1]==s || face_is_free[third_face+1])
                            faces[count+1]=third_face;count+=1
                        end
                        count<2 && continue
                        valid_division=true
                        face_done[1]=p1;face_done[2]=p2;face_done[3]=faces[1];face_done[4]=faces[2]
                        if face_is_free[s-add+1]
                            n1[s-add+1]=s;n2[s-add+1]=(s-add+1)%4+4-add
                        end
                        if face_is_free[q+1]
                            n1[q+1]=s;n2[q+1]=(s-add+3)%4+4-add
                        end
                        if face_is_free[third_face+1]
                            n1[third_face+1]=s;n2[third_face+1]=(s-add+2)%4+add
                        end
                        valid_division && break
                    end
                end
                if valid_division
                    for f in face_done
                        edge!(n1[f+1],n2[f+1])
                    end
                    break
                end
            end
        end
        if !valid_division
            for p in 0:7
                add=p<4 ? 0 : 4
                third_face=p<4 ? 4 : 5
                q=(p-add+3)%4
                if (n1[p-add+1]==p || n2[p-add+1]==p || face_is_free[p-add+1]) &&
                   (n1[q+1]==p || n2[q+1]==p || face_is_free[q+1]) &&
                   (n1[third_face+1]==p || n2[third_face+1]==p || face_is_free[third_face+1])
                    valid_division=true
                    if face_is_free[p-add+1]
                        n1[p-add+1]=p;n2[p-add+1]=(p-add+1)%4+4-add
                    end
                    if face_is_free[q+1]
                        n1[q+1]=p;n2[q+1]=(p-add+3)%4+4-add
                    end
                    if face_is_free[third_face+1]
                        n1[third_face+1]=p;n2[third_face+1]=(p-add+2)%4+add
                    end
                    face_done[1]=p-add;face_done[2]=q;face_done[3]=third_face
                    edge!(n1[p-add+1],n2[p-add+1])
                    edge!(n1[q+1],n2[q+1])
                    edge!(n1[third_face+1],n2[third_face+1])
                    break
                end
            end
        end
        valid_division && break
    end
    if !valid_division && p1_hold>=0
        valid_division=true
        face_done[1]=p1_hold;face_done[2]=p2_hold
        edge!(n1_hold[1],n2_hold[1]);edge!(n1_hold[2],n2_hold[2])
        for s_tmp in ft.free
            (s_tmp==p1_hold || s_tmp==p2_hold) && continue
            v_free=_hexa_face_vertices(V,s_tmp)
            create_forbidden!(st.forbidden_edges,v_free);create_forbidden!(forbidden_new,v_free)
        end
    end
    for ind in ft.adj_diag
        ind in face_done && continue
        edge!(nadj1[ind+1],nadj2[ind+1])
    end
    if !valid_division
        add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
    end
    return
end

function brute_force_hexa!(st::QTState,elem,j,k,verts,ft::FaceTypes,
        edges_new,forbidden_new,problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
    V=i->verts[i+1]
    E=st.quadToTri_edges
    for ind in ft.recomb
        v=_hexa_face_vertices(V,ind)
        create_forbidden!(forbidden_new,v);create_forbidden!(st.forbidden_edges,v)
    end
    if haskey(st.problems,elem) && (j,k) in st.problems[elem]
        for ind in ft.adj_diag
            create_edge!(E,V(nadj1[ind+1]),V(nadj2[ind+1]))
            create_edge!(edges_new,V(nadj1[ind+1]),V(nadj2[ind+1]))
        end
        return
    end
    nsingle=length(ft.single_tri)
    if (nsingle!=0 && nsingle!=2) || (!isempty(ft.degen) && nsingle!=2)
        # "bad degenerated extrusion encountered"
        add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
        if !isempty(ft.adj_diag)
            ind=ft.adj_diag[1]
            if nadj1[ind+1]>=0
                create_edge!(edges_new,V(nadj1[ind+1]),V(nadj2[ind+1]))
                create_edge!(E,V(nadj1[ind+1]),V(nadj2[ind+1]))
            end
        end
        return
    end
    if 6-nsingle-length(ft.recomb)<2
        if length(ft.adj_diag)+length(ft.fixed_diag)>0
            add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
            for ind in ft.adj_diag
                create_edge!(E,V(nadj1[ind+1]),V(nadj2[ind+1]))
                create_edge!(edges_new,V(nadj1[ind+1]),V(nadj2[ind+1]))
            end
        elseif !isempty(ft.free)
            v_free=_hexa_face_vertices(V,ft.free[1])
            create_forbidden!(st.forbidden_edges,v_free);create_forbidden!(forbidden_new,v_free)
        end
        return
    end
    if !isempty(ft.degen)
        two_pt_degen_hexa_edges!(st,elem,j,k,verts,ft,edges_new,problems_new,
            nfix1,nfix2,nadj1,nadj2,free_flag)
    elseif nsingle==2
        one_pt_degen_hexa_edges!(st,elem,j,k,verts,ft,edges_new,problems_new,
            nfix1,nfix2,nadj1,nadj2,free_flag)
    else
        full_hexa_edges!(st,elem,j,k,verts,ft,edges_new,forbidden_new,problems_new,
            nfix1,nfix2,nadj1,nadj2,free_flag)
    end
    return
end
