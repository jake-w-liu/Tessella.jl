# QuadToTriEdgeGenerator and the edging passes it calls
# (QuadTriExtruded3D.cpp), with CategorizedSourceElements (QuadTriUtils.cpp).

mutable struct CategorizedSourceElements
    three_bnd_pt_tri::Set{Int}
    four_bnd_pt_quad::Set{Int}
    other_bnd_tri::Set{Int}
    other_bnd_quad::Set{Int}
    internal_tri_touch_one_bnd_pt_quad::Set{Int}
    internal_quad_touch_one_bnd_pt_quad::Set{Int}
    two_bnd_pt_quad_touch_one_bnd_pt_quad::Set{Int}
    internal_tri::Set{Int}
    internal_quad::Set{Int}
    tri_bool::Vector{Bool}
    quad_bool::Vector{Bool}
end

function categorize(r::QTRegion)
    c=CategorizedSourceElements(Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),
        Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),
        falses(4length(r.tris)),falses(5length(r.quads)))
    pivots=Int[]
    for t in 0:1
        n=t==0 ? length(r.tris) : length(r.quads)
        for i in 0:n-1
            ev=t==0 ? r.tris[i+1] : r.quads[i+1]
            bnd_count=0;bnd_vert=0
            for kk in 0:length(ev)-1
                if ev[kk+1] in r.boundary
                    t==0 ? (c.tri_bool[4i+kk+2]=true) : (c.quad_bool[5i+kk+2]=true)
                    bnd_vert=kk;bnd_count+=1
                end
            end
            bnd_count>0 && (t==0 ? (c.tri_bool[4i+1]=true) : (c.quad_bool[5i+1]=true))
            if bnd_count==0
                push!(t==0 ? c.internal_tri : c.internal_quad,i)
            elseif bnd_count==1 || bnd_count==2 || (bnd_count==3 && t==1)
                push!(t==0 ? c.other_bnd_tri : c.other_bnd_quad,i)
            elseif (bnd_count==3 && t==0) || bnd_count==4
                push!(t==0 ? c.three_bnd_pt_tri : c.four_bnd_pt_quad,i)
            end
            t==1 && bnd_count==1 && push!(pivots,ev[(bnd_vert+2)%4+1])
        end
    end
    isempty(pivots) && return c
    pivots=sort!(unique!(pivots))   # std::set<MVertex*> order
    for s in 0:1,q in 0:1
        if s==0
            q==0 || continue
            elems=c.internal_tri
        else
            elems=q==0 ? c.internal_quad : c.other_bnd_quad
        end
        for i in sort!(collect(elems))
            vs=s==0 ? r.tris[i+1] : r.quads[i+1]
            if s==1 && q==1
                count(kk->c.quad_bool[5i+kk+2],0:length(vs)-1)==2 || continue
            end
            found=false
            for piv in pivots
                for tt in 0:length(vs)-1
                    if piv==vs[tt+1]
                        found=true
                        if s==0
                            push!(c.internal_tri_touch_one_bnd_pt_quad,i)
                            c.tri_bool[4i+1]=true;c.tri_bool[4i+tt+2]=true
                        elseif q==0
                            push!(c.internal_quad_touch_one_bnd_pt_quad,i)
                            c.quad_bool[5i+1]=true;c.quad_bool[5i+tt+2]=true
                        else
                            push!(c.two_bnd_pt_quad_touch_one_bnd_pt_quad,i)
                        end
                        break
                    end
                end
                found && break
            end
        end
    end
    return c
end

_tri_bnd(c,i)=c.tri_bool[4i+2:4i+4]
_quad_bnd(c,i)=c.quad_bool[5i+2:5i+5]

# getExtrudedLateralVertices for source vertices a,b at (j,k).
function lateral_vertices(r::QTRegion,a,b,j,k)
    l=_level(r,j,k)
    return [r.vid(a,l),r.vid(b,l),r.vid(a,l+1),r.vid(b,l+1)]
end

function extrude_diags!(st::QTState,r::QTRegion,v,j_start,k_start,j_top,k_top,elem,
        edges_new,forbidden_new,problems_new)
    src=element_vertices(r,elem);elem_size=length(src)
    for w in 0:1
        for (v1,v2) in sort!(collect(w==0 ? forbidden_new : edges_new))
            ind1=something(findlast(==(v1),v),0)-1
            ind2=something(findlast(==(v2),v),0)-1
            (ind1<0 || ind2<0) && error("ExtrudeDiags(): could not find vertex indices")
            sv1=ind1<elem_size ? src[ind1+1] : src[ind1-elem_size+1]
            sv2=ind2<elem_size ? src[ind2+1] : src[ind2-elem_size+1]
            for j in j_start:j_top
                k_first=j==j_start ? k_start+1 : 0
                k_stop=j==j_top ? k_top : r.layers[j+1]
                for k in k_first:k_stop-1
                    e=lateral_vertices(r,sv1,sv2,j,k)
                    if w==0
                        create_forbidden!(st.forbidden_edges,(e[1],e[2],e[4],e[3]))
                    else
                        create_edge!(st.quadToTri_edges,ind1<elem_size ? e[1] : e[3],
                            ind2<elem_size ? e[2] : e[4])
                    end
                end
            end
        end
    end
    for key in keys(problems_new)
        for j in j_start:j_top
            k_first=j==j_start ? k_start+1 : 0
            k_stop=j==j_top ? k_top : r.layers[j+1]
            for k in k_first:k_stop-1
                add_problem!(st.problems,key,j,k)
            end
        end
    end
    return
end

_all_layers(r::QTRegion,st::QTState)=!isempty(st.lat_tri_diags) || (!isempty(r.quads) && r.rotate)

function edges_all_verts_on_bnd!(st::QTState,r::QTRegion,c::CategorizedSourceElements)
    finish_all_tri=false
    tri_tmp=copy(c.three_bnd_pt_tri)
    quad_tmp=r.addverts ? Set{Int}() : copy(c.four_bnd_pt_quad)
    nlayer=length(r.layers)
    while !isempty(tri_tmp) || !isempty(quad_tmp)
        s=0
        while s<2
            set_elems=s==0 ? tri_tmp : quad_tmp
            edges_new=Set{Edge}();forbidden_new=Set{Edge}()
            problems_new=Dict{ElemId,Set{Tuple{Int,Int}}}()
            done=Int[]
            for i in sort!(collect(set_elems))
                elem=(s,i)
                vert_bnd=s==0 ? _tri_bnd(c,i) : _quad_bnd(c,i)
                all_levels=_all_layers(r,st)
                j_start,k_start=all_levels ? (0,0) : (nlayer-1,r.layers[end]-1)
                num_levels=all_levels ? sum(r.layers) : 1
                num_edged=0
                ntop1=-1;ntop2=-2
                if s==1
                    verts_top=extruded_vertices(r,elem,nlayer-1,r.layers[end]-1)
                    if edge_exists(verts_top[5],verts_top[7],st.quadToTri_edges)
                        ntop1=4;ntop2=6
                    elseif edge_exists(verts_top[6],verts_top[8],st.quadToTri_edges)
                        ntop1=5;ntop2=7
                    end
                end
                for j in j_start:nlayer-1
                    for k in (j==j_start ? k_start : 0):r.layers[j+1]-1
                        verts=extruded_vertices(r,elem,j,k)
                        V=x->verts[x+1]
                        if length(verts)==8 && j==j_start && k==k_start &&
                           !edge_exists(V(0),V(2),st.quadToTri_edges) &&
                           !edge_exists(V(1),V(3),st.quadToTri_edges)
                            v_bot=(V(0),V(1),V(2),V(3))
                            create_forbidden!(st.forbidden_edges,v_bot)
                            create_forbidden!(forbidden_new,v_bot)
                        end
                        if s==0
                            ft,nfix1,nfix2,nadj1,nadj2,free_flag=get_face_types(r,j,k,verts,st,vert_bnd)
                            if finish_all_tri || length(ft.free)+length(ft.adj_diag)<2
                                brute_force_prism!(st,elem,j,k,verts,ft,edges_new,forbidden_new,
                                    problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
                                num_edged+=1
                            end
                        else
                            changed_diag=false
                            for g in 0:1
                                if g==0 && ntop1>=0 && ntop2>=0 &&
                                   !edge_exists(V(4),V(6),st.quadToTri_edges) &&
                                   !edge_exists(V(5),V(7),st.quadToTri_edges)
                                    v_top=(V(4),V(5),V(6),V(7))
                                    if !forbidden_exists(v_top,st.forbidden_edges)
                                        changed_diag=true
                                        create_edge!(st.quadToTri_edges,V(ntop1),V(ntop2))
                                    end
                                end
                                if g==1 && changed_diag && haskey(st.problems,elem) &&
                                   (j,k) in st.problems[elem]
                                    delete!(st.problems[elem],(j,k))
                                    isempty(st.problems[elem]) && delete!(st.problems,elem)
                                    delete_edge!(st.quadToTri_edges,V(ntop1),V(ntop2))
                                elseif g==1
                                    haskey(st.problems,elem) && isempty(st.problems[elem]) &&
                                        delete!(st.problems,elem)
                                    break
                                end
                                ft,nfix1,nfix2,nadj1,nadj2,free_flag=get_face_types(r,j,k,verts,st,vert_bnd)
                                brute_force_hexa!(st,elem,j,k,verts,ft,edges_new,forbidden_new,
                                    problems_new,nfix1,nfix2,nadj1,nadj2,free_flag)
                            end
                            num_edged+=1
                        end
                    end
                end
                num_edged==num_levels && push!(done,i)
                s==1 && break
            end
            for i in done
                delete!(set_elems,i)
            end
            if s==0 && (!isempty(edges_new) || !isempty(problems_new))
                s=-1
            end
            isempty(quad_tmp) && (finish_all_tri=true)
            s+=1
        end
    end
    return
end

function edges_other_bnd_prisms!(st::QTState,r::QTRegion,c::CategorizedSourceElements)
    (isempty(st.lat_tri_diags) && isempty(r.quads)) && return
    nlayer=length(r.layers)
    E=st.quadToTri_edges;L=st.lat_tri_diags
    for i in sort!(collect(c.other_bnd_tri))
        elem=(0,i);elem_size=3
        vert_bnd=_tri_bnd(c,i)
        j_start,k_start=_all_layers(r,st) ? (0,0) : (nlayer-1,r.layers[end]-1)
        for j in j_start:nlayer-1, k in (j==j_start ? k_start : 0):r.layers[j+1]-1
            verts=extruded_vertices(r,elem,j,k);V=x->verts[x+1]
            for p in 0:elem_size-1
                vert_bnd[p+1] || continue
                p2=(p+1)%elem_size;p3=(p+elem_size-1)%elem_size
                degen=V(p)==V(p+elem_size)
                !degen && !vert_bnd[p2+1] && create_edge!(E,V(p),V(p2+elem_size))
                !degen && !vert_bnd[p3+1] && create_edge!(E,V(p),V(p3+elem_size))
                if edge_exists(V(p),V(p2+elem_size),L) && !edge_exists(V(p+elem_size),V(p2),E)
                    create_edge!(E,V(p),V(p2+elem_size))
                elseif edge_exists(V(p+elem_size),V(p2),L) && !edge_exists(V(p),V(p2+elem_size),E)
                    create_edge!(E,V(p+elem_size),V(p2))
                end
            end
        end
    end
    return
end

function _two_bnd_pair(vert_bnd,elem_size)
    p1=-1;p2=-1
    for p in 0:elem_size-1
        if vert_bnd[p+1]
            p1=p
            p2_tmp=(p+1)%elem_size;p3_tmp=(p+elem_size-1)%elem_size
            if vert_bnd[p2_tmp+1]
                p1=p;p2=p2_tmp
            elseif vert_bnd[p3_tmp+1]
                p1=p3_tmp;p2=p
            end
            break
        end
    end
    return p1,p2
end

function edges_other_bnd_hexa!(st::QTState,r::QTRegion,c::CategorizedSourceElements)
    nlayer=length(r.layers)
    E=st.quadToTri_edges;L=st.lat_tri_diags
    j_top_start=nlayer-1;k_top_start=r.layers[end]-1
    for i in sort!(collect(c.other_bnd_quad))
        elem=(1,i);elem_size=4
        vert_bnd=_quad_bnd(c,i)
        bnd_count=0;skip=0;one_point_ind=0
        for s in 0:3
            if vert_bnd[s+1]
                bnd_count+=1;one_point_ind=s
            else
                skip=s
            end
        end
        j_start,k_start=_all_layers(r,st) ? (0,0) : (nlayer-1,r.layers[end]-1)
        for j in j_start:nlayer-1, k in (j==j_start ? k_start : 0):r.layers[j+1]-1
            verts=extruded_vertices(r,elem,j,k);V=x->verts[x+1]
            for p in 0:elem_size-1
                p2=(p+1)%elem_size;p3=(p+elem_size-1)%elem_size
                if vert_bnd[p+1]
                    degen=V(p)==V(p+elem_size)
                    !degen && !vert_bnd[p2+1] && create_edge!(E,V(p),V(p2+elem_size))
                    !degen && !vert_bnd[p3+1] && create_edge!(E,V(p),V(p3+elem_size))
                    if edge_exists(V(p),V(p2+elem_size),L) && !edge_exists(V(p+elem_size),V(p2),E)
                        create_edge!(E,V(p),V(p2+elem_size))
                    elseif edge_exists(V(p+elem_size),V(p2),L) && !edge_exists(V(p),V(p2+elem_size),E)
                        create_edge!(E,V(p+elem_size),V(p2))
                    end
                elseif j==j_top_start && k==k_top_start && !vert_bnd[p2+1]
                    if V(p+elem_size)<V(p2+elem_size)
                        create_edge!(E,V(p+elem_size),V(p2))
                    else
                        create_edge!(E,V(p2+elem_size),V(p))
                    end
                end
            end
        end
        if bnd_count==2
            done=false
            for j in j_start:nlayer-1
                k_stop=r.layers[j+1]-(j==nlayer-1 ? 1 : 0)
                for k in (j==j_start ? k_start : 0):k_stop-1
                    verts=extruded_vertices(r,elem,j,k);V=x->verts[x+1]
                    p1,p2=_two_bnd_pair(vert_bnd,elem_size)
                    if p2<0 && p1<0
                        break
                    end
                    if edge_exists(V(0),V(2),E)
                        create_edge!(E,V(4),V(6))
                    elseif edge_exists(V(1),V(3),E)
                        create_edge!(E,V(5),V(7))
                    elseif p2<0 && p1>=0
                        create_edge!(E,V((p1+1)%elem_size+elem_size),V((p1+elem_size-1)%elem_size+elem_size))
                    elseif p1>=0 && edge_exists(V(p1),V(p2+elem_size),E)
                        create_edge!(E,V(p2+elem_size),V((p2+2)%elem_size+elem_size))
                    elseif p1>=0 && edge_exists(V(p1+elem_size),V(p2),E)
                        create_edge!(E,V(p1+elem_size),V((p1+2)%elem_size+elem_size))
                    end
                end
            end
        elseif bnd_count==3 && (j_start<j_top_start || k_start<k_top_start)
            verts=extruded_vertices(r,elem,j_start,k_start);V=x->verts[x+1]
            edges_new=Set{Edge}();forbidden_new=Set{Edge}()
            problems_new=Dict{ElemId,Set{Tuple{Int,Int}}}()
            create_edge!(E,V(skip+elem_size),V((skip+2)%elem_size+elem_size))
            create_edge!(edges_new,V(skip+elem_size),V((skip+2)%elem_size+elem_size))
            extrude_diags!(st,r,verts,j_start,k_start,j_top_start,k_top_start,elem,
                edges_new,forbidden_new,problems_new)
        elseif bnd_count==1 && (j_start<j_top_start || k_start<k_top_start)
            verts=extruded_vertices(r,elem,j_start,k_start);V=x->verts[x+1]
            p=one_point_ind
            if r.layers[1]>1
                j_next=0;k_next=1
            else
                j_next=min(nlayer-1,1);k_next=0
            end
            divide_next=false;verts_next=Int[]
            if j_next<j_top_start || k_next<k_top_start
                divide_next=true
                verts_next=extruded_vertices(r,elem,j_next,k_next)
            end
            W=x->verts_next[x+1]
            edges_new=Set{Edge}();forbidden_new=Set{Edge}()
            problems_new=Dict{ElemId,Set{Tuple{Int,Int}}}()
            p2=(p+1)%elem_size;p3=(p+2)%elem_size;p4=(p+elem_size-1)%elem_size
            edge_exists(V(p2+elem_size),V(p3),E) || create_edge!(E,V(p2),V(p3+elem_size))
            edge_exists(V(p4+elem_size),V(p3),E) || create_edge!(E,V(p4),V(p3+elem_size))
            divide_next && !edge_exists(W(p2),W(p3+elem_size),E) &&
                create_edge!(E,W(p2+elem_size),W(p3))
            divide_next && !edge_exists(W(p4),W(p3+elem_size),E) &&
                create_edge!(E,W(p4+elem_size),W(p3))
            if V(p)==V(p+elem_size) &&
               (edge_exists(V(p2),V(p3+elem_size),E) || V(p2)==V(p2+elem_size)) &&
               (edge_exists(V(p4),V(p3+elem_size),E) || V(p4)==V(p4+elem_size)) &&
               !edge_exists(V(p2+elem_size),V(p4+elem_size),E)
                create_edge!(E,V(p+elem_size),V(p3+elem_size))
                create_edge!(edges_new,V(p+elem_size),V(p3+elem_size))
            else
                create_edge!(E,V(p2+elem_size),V(p4+elem_size))
                create_edge!(edges_new,V(p2+elem_size),V(p4+elem_size))
            end
            extrude_diags!(st,r,verts,0,0,j_top_start,k_top_start,elem,
                edges_new,forbidden_new,problems_new)
        end
    end
    return
end

function edges_touch_pivot_vert!(st::QTState,r::QTRegion,c::CategorizedSourceElements)
    (isempty(st.lat_tri_diags) && (isempty(r.quads) || !r.rotate)) && return
    nlayer=length(r.layers)
    (nlayer<=1 && r.layers[1]<=1) && return
    E=st.quadToTri_edges
    j_top_start=nlayer-1;k_top_start=r.layers[end]-1
    for s in 0:1
        set_elems=s==0 ? c.internal_tri_touch_one_bnd_pt_quad : c.internal_quad_touch_one_bnd_pt_quad
        for i in sort!(collect(set_elems))
            elem=(s,i);elem_size=s==0 ? 3 : 4
            v_index=-1
            for p in 0:elem_size-1
                if (s==0 && c.tri_bool[4i+2+p]) || (s==1 && c.quad_bool[5i+2+p])
                    v_index=p;break
                end
            end
            v_index<0 && continue
            verts=extruded_vertices(r,elem,0,0);V=x->verts[x+1]
            if r.layers[1]>1
                j_next=0;k_next=1
            else
                j_next=min(nlayer-1,1);k_next=0
            end
            divide_next=false;verts_next=Int[]
            if j_next<j_top_start || k_next<k_top_start
                divide_next=true
                verts_next=extruded_vertices(r,elem,j_next,k_next)
            end
            W=x->verts_next[x+1]
            p2=(v_index+1)%elem_size;p3=(v_index+elem_size-1)%elem_size
            if !edge_exists(V(v_index),V(p2+elem_size),E)
                create_edge!(E,V(v_index+elem_size),V(p2))
                divide_next && create_edge!(E,W(v_index),W(p2+elem_size))
            end
            if !edge_exists(V(v_index),V(p3+elem_size),E)
                create_edge!(E,V(v_index+elem_size),V(p3))
                divide_next && create_edge!(E,W(v_index),W(p3+elem_size))
            end
            if s==1 && V(v_index+elem_size)!=V(v_index) &&
               !edge_exists(V(p2+elem_size),V(p3+elem_size),E)
                create_edge!(E,V(v_index+elem_size),V((v_index+2)%elem_size+elem_size))
            end
        end
    end
    for i in sort!(collect(c.two_bnd_pt_quad_touch_one_bnd_pt_quad))
        elem=(1,i);elem_size=4
        vert_bnd=_quad_bnd(c,i)
        count(vert_bnd)==2 || continue
        verts=extruded_vertices(r,elem,0,0);V=x->verts[x+1]
        p1,p2=_two_bnd_pair(vert_bnd,elem_size)
        (p2<0 && p1<0) && break
        edges_new=Set{Edge}();forbidden_new=Set{Edge}()
        if p1>=0 && p2>=0 &&
           edge_exists(V((p1+2)%elem_size),V((p1+elem_size-1)%elem_size+elem_size),E) &&
           !edge_exists(V(p1+elem_size),V((p1+2)%elem_size+elem_size),E)
            a=V((p1+1)%elem_size+elem_size);b=V((p1+elem_size-1)%elem_size+elem_size)
            create_edge!(E,a,b);create_edge!(edges_new,a,b)
        elseif p1>=0 && p2>=0 &&
           edge_exists(V((p1+2)%elem_size+elem_size),V((p1+elem_size-1)%elem_size),E) &&
           !edge_exists(V((p1+1)%elem_size+elem_size),V((p1+elem_size-1)%elem_size+elem_size),E)
            a=V(p1+elem_size);b=V((p1+2)%elem_size+elem_size)
            create_edge!(E,a,b);create_edge!(edges_new,a,b)
        end
        if !isempty(edges_new)
            # Gmsh passes a local, discarded problem map here.
            local_state=QTState(st.quadToTri_edges,st.forbidden_edges,st.lat_tri_diags,
                Dict{ElemId,Set{Tuple{Int,Int}}}())
            extrude_diags!(local_state,r,verts,0,0,j_top_start,k_top_start,elem,
                edges_new,forbidden_new,Dict{ElemId,Set{Tuple{Int,Int}}}())
        end
    end
    return
end

function edges_internal_top_layer!(st::QTState,r::QTRegion,c::CategorizedSourceElements)
    isempty(r.quads) && return
    E=st.quadToTri_edges
    for s in 0:1
        set_elems=s==0 ? c.internal_tri : c.internal_quad
        for i in sort!(collect(set_elems))
            elem=(s,i);elem_size=s==0 ? 3 : 4
            verts=extruded_vertices(r,elem,length(r.layers)-1,r.layers[end]-1);V=x->verts[x+1]
            for p in 0:elem_size-1
                ind1=p+elem_size;ind2=(p+1)%elem_size+elem_size
                if V(ind1)<V(ind2)
                    ind_low=ind1;ind_low_2=ind2-elem_size
                else
                    ind_low=ind2;ind_low_2=ind1-elem_size
                end
                edge_exists(V(ind_low-elem_size),V(ind_low_2+elem_size),E) ||
                    create_edge!(E,V(ind_low),V(ind_low_2))
            end
        end
    end
    return
end

# QuadToTriGetRegionDiags on the laterals' first meshes (`lateral_cells[s]`:
# per line, per level, the cells createQuaTri made) and the top diagonals,
# then the edging passes of QuadToTriEdgeGenerator. A lone region's laterals
# are not shared with another QuadToTri region, so their quadrangles are
# forbidden, and they are extruded entities, so their diagonals adjustable.
function edge_generator!(st::QTState,r::QTRegion,c::CategorizedSourceElements,
        lateral_cells,top_diagonals)
    for (s,lines) in enumerate(r.lateral_lines)
        for (li,(v0,v1)) in enumerate(lines)
            for cell in lateral_cells[s][li]
                length(cell)==4 && create_forbidden!(st.forbidden_edges,cell)
            end
            # A line whose end lies on the axis gives single triangles only.
            (r.vid(v0,0)==r.vid(v0,1) || r.vid(v1,0)==r.vid(v1,1)) && continue
            any(cell->length(cell)==3,lateral_cells[s][li]) || continue
            for level in 0:sum(r.layers)-1
                verts=(r.vid(v0,level),r.vid(v1,level),r.vid(v0,level+1),r.vid(v1,level+1))
                tris=[cell for cell in lateral_cells[s][li] if length(cell)==3 &&
                      all(x->x in verts,cell)]
                isempty(tris) && continue
                t=first(tris)
                if (verts[1] in t) && (verts[4] in t)
                    create_edge!(st.lat_tri_diags,verts[1],verts[4])
                elseif (verts[2] in t) && (verts[3] in t)
                    create_edge!(st.lat_tri_diags,verts[2],verts[3])
                end
            end
        end
    end
    for d in top_diagonals
        create_edge!(st.quadToTri_edges,d...)
    end
    r.addverts && union!(st.quadToTri_edges,st.lat_tri_diags)
    (isempty(st.lat_tri_diags) && isempty(r.quads)) && return
    r.addverts && return
    edges_all_verts_on_bnd!(st,r,c)
    edges_other_bnd_prisms!(st,r,c)
    edges_other_bnd_hexa!(st,r,c)
    edges_touch_pivot_vert!(st,r,c)
    edges_internal_top_layer!(st,r,c)
    return
end
