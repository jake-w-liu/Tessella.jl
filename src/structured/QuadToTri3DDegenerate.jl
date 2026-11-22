# Brute-force diagonals of degenerate swept hexahedra (QuadTriExtruded3D.cpp:
# addEdgesForQuadToTriTwoPtDegenHexa, addEdgesForQuadToTriOnePtDegenHexa).

# The diagonal choice shared by the degenerate-hexahedron loops for face
# p_tmp at pass t: (n1,n2,is_free) after the fixed / preferred /
# lowest-pointer / free rules; n1,n2 keep their sentinels otherwise.
function _hexa_face_choice(V,t,p_tmp,nfix1,nfix2,nadj1,nadj2,free_flag,n1,n2)
    is_free=false
    if nfix1[p_tmp+1]>=0
        n1=nfix1[p_tmp+1];n2=nfix2[p_tmp+1]
    elseif t==0 && nadj1[p_tmp+1]>=0
        n1=nadj1[p_tmp+1];n2=nadj2[p_tmp+1]
    elseif (t>=1 && t<3 && nadj1[p_tmp+1]>=0) || (t==2 && free_flag[p_tmp+1]!=0)
        if p_tmp<4
            if (V(p_tmp)<V((p_tmp+1)%4) && V(p_tmp)<V(p_tmp+4)) ||
               (V((p_tmp+1)%4+4)<V((p_tmp+1)%4) && V((p_tmp+1)%4+4)<V(p_tmp+4))
                n1=p_tmp;n2=(p_tmp+1)%4+4
            else
                n1=p_tmp+4;n2=(p_tmp+1)%4
            end
        else
            add=p_tmp==4 ? 0 : 4
            if (V(0+add)<V(1+add) && V(0+add)<V(3+add)) ||
               (V(2+add)<V(1+add) && V(2+add)<V(3+add))
                n1=0+add;n2=2+add
            else
                n1=1+add;n2=3+add
            end
        end
    elseif t==3 && (nadj1[p_tmp+1]>=0 || free_flag[p_tmp+1]!=0)
        is_free=true
    end
    return n1,n2,is_free
end

function two_pt_degen_hexa_edges!(st::QTState,elem,j,k,verts,ft::FaceTypes,
        edges_new,problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
    V=i->verts[i+1]
    E=st.quadToTri_edges
    edge!(a,b)=(create_edge!(E,V(a),V(b));create_edge!(edges_new,V(a),V(b)))
    num_free=length(ft.free)
    degen_face=ft.degen[1]
    if num_free>=3
        ind_low=0
        for p in 1:7
            V(ind_low)>V(p) && (ind_low=p)
        end
        add=ind_low<4 ? 0 : 4
        if V(ind_low-add)==V(ind_low+4-add)
            edge!(ind_low,(ind_low-add+2)%4)
            edge!(ind_low,(ind_low-add+2)%4+4)
        else
            edge!(ind_low,(ind_low-add+2)%4+4-add)
            ind2=V((ind_low-add+1)%4)==V((ind_low-add+1)%4+4) ?
                (ind_low-add+3)%4+4-add : (ind_low-add+1)%4+4-add
            edge!(ind_low,ind2)
        end
        return
    end
    face_done=[-1,-2]
    valid_division=false
    for t in 0:3
        p_top=5;p_bot=4;p_lat=(degen_face+2)%4
        n1_bot,n2_bot,p_bot_is_free=_hexa_face_choice(V,t,p_bot,nfix1,nfix2,nadj1,nadj2,free_flag,-4,-5)
        n1_top,n2_top,p_top_is_free=_hexa_face_choice(V,t,p_top,nfix1,nfix2,nadj1,nadj2,free_flag,-2,-3)
        n1_lat,n2_lat,p_lat_is_free=_hexa_face_choice(V,t,p_lat,nfix1,nfix2,nadj1,nadj2,free_flag,-10,-11)
        if (n1_top>=0 || p_top_is_free) && (n1_bot>=0 || p_bot_is_free)
            if p_top_is_free && p_bot_is_free
                if (V(4)<V(5) && V(4)<V(7)) || (V(6)<V(5) && V(6)<V(7))
                    n1_top=4;n2_top=6
                else
                    n1_top=5;n2_top=7
                end
                n1_bot=n1_top-4;n2_bot=n2_top-4
            elseif n1_top>=0 && p_bot_is_free
                n1_bot=n1_top-4;n2_bot=n2_top-4
            elseif n1_bot>=0 && p_top_is_free
                n1_top=n1_bot+4;n2_top=n2_bot+4
            end
            if V(n1_top)==V(n1_bot) || V(n1_top)==V(n2_bot) ||
               V(n2_top)==V(n1_bot) || V(n2_top)==V(n2_bot)
                valid_division=true;face_done[1]=p_top;face_done[2]=p_bot
                edge!(n1_top,n2_top);edge!(n1_bot,n2_bot)
            end
        end
        if !valid_division && (n1_top>=0 || p_top_is_free) && (n1_lat>=0 || p_lat_is_free)
            a=(degen_face+2)%4+4
            if p_top_is_free && p_lat_is_free
                if (V(4)<V(5) && V(4)<V(7)) || (V(6)<V(5) && V(6)<V(7))
                    n1_top=4;n2_top=6
                else
                    n1_top=5;n2_top=7
                end
                if n1_top==a || n2_top==a
                    n1_lat=a;n2_lat=(n1_lat-4+1)%4
                else
                    n1_lat=(degen_face+3)%4+4;n2_lat=(n1_lat-4+3)%4
                end
            elseif n1_top>=0 && p_lat_is_free
                if n1_top==a || n2_top==a
                    n1_lat=a;n2_lat=(n1_lat-4+1)%4
                else
                    n1_lat=(degen_face+3)%4+4;n2_lat=(n1_lat-4+3)%4
                end
            elseif n1_lat>=0 && p_top_is_free
                if n1_lat==a || n2_lat==a
                    n1_top=a;n2_top=(n1_top-4+2)%4+4
                else
                    n1_top=(degen_face+3)%4+4;n2_top=(n1_top-4+2)%4+4
                end
            end
            if V(n1_top)==V(n1_lat) || V(n1_top)==V(n2_lat) ||
               V(n2_top)==V(n1_lat) || V(n2_top)==V(n2_lat)
                valid_division=true;face_done[1]=p_top;face_done[2]=p_lat
                edge!(n1_top,n2_top);edge!(n1_lat,n2_lat)
            end
        end
        if !valid_division && (n1_bot>=0 || p_bot_is_free) && (n1_lat>=0 || p_lat_is_free)
            a=(degen_face+2)%4
            if p_bot_is_free && p_lat_is_free
                if (V(0)<V(1) && V(0)<V(3)) || (V(2)<V(1) && V(2)<V(3))
                    n1_bot=0;n2_bot=2
                else
                    n1_bot=1;n2_bot=3
                end
                if n1_bot==a || n2_bot==a
                    n1_lat=a;n2_lat=(n1_lat+1)%4+4
                else
                    n1_lat=(degen_face+3)%4;n2_lat=(n1_lat+3)%4+4
                end
            elseif n1_bot>=0 && p_lat_is_free
                if n1_bot==a || n2_bot==a
                    n1_lat=a;n2_lat=(n1_lat+1)%4+4
                else
                    n1_lat=(degen_face+3)%4;n2_lat=(n1_lat+3)%4+4
                end
            elseif n1_lat>=0 && p_bot_is_free
                if n1_lat==a || n2_lat==a
                    n1_bot=a;n2_bot=(n1_bot+2)%4
                else
                    n1_bot=(degen_face+3)%4;n2_bot=(n1_bot+2)%4
                end
            end
            if V(n1_bot)==V(n1_lat) || V(n1_bot)==V(n2_lat) ||
               V(n2_bot)==V(n1_lat) || V(n2_bot)==V(n2_lat)
                valid_division=true;face_done[1]=p_bot;face_done[2]=p_lat
                edge!(n1_bot,n2_bot);edge!(n1_lat,n2_lat)
            end
        end
        valid_division && break
    end
    for ind in ft.adj_diag
        (ind!=face_done[1] && ind!=face_done[2]) && edge!(nadj1[ind+1],nadj2[ind+1])
    end
    if !valid_division
        add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
    end
    return
end

function one_pt_degen_hexa_edges!(st::QTState,elem,j,k,verts,ft::FaceTypes,
        edges_new,problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
    V=i->verts[i+1]
    E=st.quadToTri_edges
    edge!(a,b)=(create_edge!(E,V(a),V(b));create_edge!(edges_new,V(a),V(b)))
    num_free=length(ft.free)
    s0=ft.single_tri[1]
    degen_ind=V(s0)==V(s0+4) ? s0 : ft.single_tri[2]
    if num_free>=4
        edge!(degen_ind,(degen_ind+2)%4)
        edge!(degen_ind,(degen_ind+2)%4+4)
        return
    end
    face_done=[-1,-2,-3]
    valid_division=false
    d1=(degen_ind+1)%4;d2=(degen_ind+2)%4;d3=(degen_ind+3)%4
    for t in 0:3
        n1_top,n2_top,p_top_is_free=_hexa_face_choice(V,t,5,nfix1,nfix2,nadj1,nadj2,free_flag,-2,-3)
        n1_bot,n2_bot,p_bot_is_free=_hexa_face_choice(V,t,4,nfix1,nfix2,nadj1,nadj2,free_flag,-4,-5)
        n1_lat1,n2_lat1,p_lat1_is_free=_hexa_face_choice(V,t,d1,nfix1,nfix2,nadj1,nadj2,free_flag,-10,-11)
        n1_lat2,n2_lat2,p_lat2_is_free=_hexa_face_choice(V,t,d2,nfix1,nfix2,nadj1,nadj2,free_flag,-12,-13)
        if !valid_division
            if p_top_is_free && p_bot_is_free
                if (V(0)<V(1) && V(0)<V(3)) || (V(2)<V(1) && V(2)<V(3))
                    n1_bot=0;n2_bot=2
                else
                    n1_bot=1;n2_bot=3
                end
                n1_top=n1_bot+4;n2_top=n2_bot+4
            elseif p_top_is_free && n1_bot>=0
                n1_top=n1_bot+4;n2_top=n2_bot+4
            elseif p_bot_is_free && n1_top>=0
                n1_bot=n1_top-4;n2_bot=n2_top-4
            end
            if n1_top-4==n1_bot || n2_top-4==n1_bot || n1_top-4==n2_bot || n2_top-4==n2_bot
                valid_division=true;face_done[1]=4;face_done[2]=5
                edge!(n1_top,n2_top);edge!(n1_bot,n2_bot)
            end
        end
        if !valid_division
            if (n1_top==d2+4 || n2_top==d2+4 || p_top_is_free) &&
               (n2_lat1==d2+4 || p_lat1_is_free) && (n1_lat2==d2+4 || p_lat2_is_free)
                valid_division=true
                p_top_is_free && (n1_top=degen_ind+4;n2_top=d2+4)
                p_lat1_is_free && (n1_lat1=d1;n2_lat1=d2+4)
                p_lat2_is_free && (n1_lat2=d2+4;n2_lat2=d3)
                face_done[1]=5;face_done[2]=d1;face_done[3]=d2
                edge!(n1_top,n2_top)
            elseif (n1_bot==d2 || n2_bot==d2 || p_bot_is_free) &&
                   (n2_lat1==d2 || p_lat1_is_free) && (n1_lat2==d2 || p_lat2_is_free)
                valid_division=true
                p_bot_is_free && (n1_bot=degen_ind;n2_bot=d2)
                p_lat1_is_free && (n1_lat1=d1+4;n2_lat1=d2)
                p_lat2_is_free && (n1_lat2=d2;n2_lat2=d3+4)
                face_done[1]=4;face_done[2]=d1;face_done[3]=d2
                edge!(n1_bot,n2_bot)
            end
            if valid_division
                edge!(n1_lat1,n2_lat1);edge!(n1_lat2,n2_lat2)
            end
        end
        if !valid_division
            if (n1_top>=0 && n1_top!=d2+4 && n2_top!=d2+4) || p_top_is_free
                if n1_lat1==d1+4 || p_lat1_is_free
                    valid_division=true;face_done[1]=d1
                    p_lat1_is_free && (n1_lat1=d1+4;n2_lat1=d2)
                    edge!(n1_lat1,n2_lat1)
                elseif n2_lat2==d3+4 || p_lat2_is_free
                    valid_division=true;face_done[1]=d2
                    p_lat2_is_free && (n1_lat2=d2;n2_lat2=d3+4)
                    edge!(n1_lat2,n2_lat2)
                end
                if valid_division
                    face_done[2]=5
                    p_top_is_free && (n1_top=d1+4;n2_top=d3+4)
                    edge!(n1_top,n2_top)
                end
            end
            if (!valid_division && n1_bot>=0 && n1_bot!=d2 && n2_bot!=d2) || p_bot_is_free
                if n1_lat1==d1 || p_lat1_is_free
                    valid_division=true;face_done[1]=d1
                    p_lat1_is_free && (n1_lat1=d1;n2_lat1=d2+4)
                    edge!(n1_lat1,n2_lat1)
                elseif n2_lat2==d3 || p_lat2_is_free
                    valid_division=true;face_done[1]=d2
                    p_lat2_is_free && (n1_lat2=d2+4;n2_lat2=d3)
                    edge!(n1_lat2,n2_lat2)
                end
                if valid_division
                    face_done[2]=4
                    p_bot_is_free && (n1_bot=d1;n2_bot=d3)
                    edge!(n1_bot,n2_bot)
                end
            end
        end
        if !valid_division
            if (n1_lat1==d1 || p_lat1_is_free) && (n1_lat2==d2 || p_lat2_is_free)
                valid_division=true
                n1_lat1=d1;n2_lat1=d2+4;n1_lat2=d2;n2_lat2=d3+4
            elseif (n1_lat1==d1+4 || p_lat1_is_free) && (n1_lat2==d2+4 || p_lat2_is_free)
                valid_division=true
                n1_lat1=d1+4;n2_lat1=d2;n1_lat2=d2+4;n2_lat2=d3
            end
            if valid_division
                face_done[1]=d1;face_done[2]=d2
                edge!(n1_lat1,n2_lat1);edge!(n1_lat2,n2_lat2)
            end
        end
        valid_division && break
    end
    for ind in ft.adj_diag
        (ind!=face_done[1] && ind!=face_done[2] && ind!=face_done[3]) &&
            edge!(nadj1[ind+1],nadj2[ind+1])
    end
    if !valid_division
        add_problem!(st.problems,elem,j,k);add_problem!(problems_new,elem,j,k)
    end
    return
end
