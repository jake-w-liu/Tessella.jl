# createTwoPtDegenHexElems, createOnePtDegenHexElems, createFullHexElems,
# QuadToTriHexPri, addBodyCenteredVertices and QuadToTriCreateElements
# (QuadTriExtruded3D.cpp).

function two_pt_degen_hex_elems!(out::QTMesh,v,dup,n1,n2,found_diags)
    V=i->v[i+1];N1=i->n1[i+1];N2=i->n2[i+1]
    degen_face=(dup[1]+1)%4==dup[2] ? dup[1] : dup[2]
    if !found_diags
        add_pri!(out,V(degen_face),V((degen_face+3)%4),V((degen_face+3)%4+4),
            V((degen_face+1)%4),V((degen_face+2)%4),V((degen_face+2)%4+4))
        return true
    end
    tet_top=-1;pyr_top=-2
    if N1(4)>=0 && (N1(4)+4==N1(5) || N2(4)+4==N1(5))
        if N1(4)==degen_face || N2(4)==degen_face
            pyr_top=degen_face;tet_top=(degen_face+1)%4
        else
            pyr_top=(degen_face+1)%4;tet_top=degen_face
        end
        add_tet!(out,V(pyr_top),V((pyr_top+2)%4),V((pyr_top+2)%4+4),V(tet_top))
        base=(degen_face+2)%4
        if N1(base)>=0
            add=N1(base)>3 ? -4 : 4
            add_tet!(out,V(N1(base)),V(N2(base)),V(N1(base)+add),V(pyr_top))
            add_tet!(out,V(N1(base)),V(N2(base)-add),V(N2(base)),V(pyr_top))
        else
            add_pyr!(out,V(base),V(base+4),V((base+1)%4+4),V((base+1)%4),V(pyr_top))
        end
        return true
    end
    d2=(degen_face+2)%4
    if tet_top<0 && N1(4)>=0 && (N1(4)==N1(d2) || N2(4)==N1(d2) || N1(4)==N2(d2) || N2(4)==N2(d2))
        pyr_top=N1(d2)<4 ? N1(d2) : N2(d2)
        tet_top=pyr_top==d2 ? (degen_face+3)%4 : d2
        add_tet!(out,V(pyr_top),V(tet_top+4),V((pyr_top+2)%4),V(tet_top))
        base=5
        if N1(base)>=0
            add_tet!(out,V(N1(base)),V(N2(base)),V((N1(base)-4+1)%4+4),V(pyr_top))
            add_tet!(out,V(N1(base)),V(N2(base)),V((N1(base)-4+3)%4+4),V(pyr_top))
        else
            add_pyr!(out,V(4),V(5),V(6),V(7),V(pyr_top))
        end
        return true
    end
    if tet_top<0 && N1(5)>=0 && (N1(5)==N1(d2) || N2(5)==N1(d2) || N1(5)==N2(d2) || N2(5)==N2(d2))
        pyr_top=N1(d2)<4 ? N2(d2) : N1(d2)
        tet_top=pyr_top==d2+4 ? (degen_face+3)%4+4 : d2+4
        add_tet!(out,V(pyr_top),V((pyr_top-4+2)%4+4),V(tet_top-4),V(tet_top))
        base=4
        if N1(base)>=0
            add_tet!(out,V(N1(base)),V(N2(base)),V((N1(base)+1)%4),V(pyr_top))
            add_tet!(out,V(N1(base)),V(N2(base)),V((N1(base)+3)%4),V(pyr_top))
        else
            add_pyr!(out,V(0),V(1),V(2),V(3),V(pyr_top))
        end
        return true
    end
    return false
end

function one_pt_degen_hex_elems!(out::QTMesh,v,dup,n1,n2,found_diags)
    V=i->v[i+1];N1=i->n1[i+1];N2=i->n2[i+1]
    if !found_diags
        add_hex!(out,V(0),V(1),V(2),V(3),V(4),V(5),V(6),V(7))
        return true
    end
    d0=dup[1];d1=(d0+1)%4;d2=(d0+2)%4;d3=(d0+3)%4
    if (N1(4)==d0 || N2(4)==d0) && (N1(5)==d0+4 || N2(5)==d0+4)
        if N1(d1)>=0
            if N1(d1)==d1
                add_tet!(out,V(d1),V(d1+4),V(d2+4),V(d0))
                add_tet!(out,V(d1),V(d2+4),V(d2),V(d0))
            elseif N1(d1)==d1+4
                add_tet!(out,V(d1+4),V(d2),V(d1),V(d0))
                add_tet!(out,V(d1+4),V(d2+4),V(d2),V(d0))
            end
        else
            add_pyr!(out,V(d1),V(d1+4),V(d2+4),V(d2),V(d0))
        end
        if N1(d2)>=0
            if N1(d2)==d2
                add_tet!(out,V(d2),V(d2+4),V(d3+4),V(d0))
                add_tet!(out,V(d2),V(d3+4),V(d3),V(d0))
            elseif N1(d2)==d2+4
                add_tet!(out,V(d2+4),V(d3),V(d2),V(d0))
                add_tet!(out,V(d2+4),V(d3+4),V(d3),V(d0))
            end
        else
            add_pyr!(out,V(d2),V(d2+4),V(d3+4),V(d3),V(d0))
        end
        return true
    elseif (N1(4)==d1 || N2(4)==d1) && (N1(5)==d1+4 || N2(5)==d1+4)
        tet_top=-1;n_inner=[-1,-2]
        if N1(d1)>=0 && N1(d2)!=N2(d1)
            tet_top=N1(d1);n_inner[1]=N1(d1)
            add=n_inner[1]>3 ? 0 : 4
            n_inner[2]=d3+add
        elseif N1(d2)>=0 && N1(d2)!=N2(d1)
            tet_top=N2(d2);n_inner[2]=N2(d2)
            add=n_inner[2]>3 ? 0 : 4
            n_inner[1]=d1+add
        elseif N1(d1)>=0
            tet_top=N2(d1)
        end
        if tet_top>=0
            add=tet_top>3 ? 0 : 4
            add_tet!(out,V(d1+add),V(d2+add),V(d3+add),V(tet_top))
            if n_inner[1]<0
                add_pyr!(out,V(d1),V(d1+4),V(d3+4),V(d3),V(tet_top))
            else
                base=tet_top-4+add==d1 ? d2 : d1
                if N1(base)<0
                    add_pyr!(out,V(base),V(base+4),V((base+1)%4+4),V((base+1)%4),V(tet_top))
                else
                    add_2=N1(base)>3 ? -4 : 4
                    add_tet!(out,V(N1(base)),V(N1(base)+add_2),V(N2(base)),V(tet_top))
                    add_tet!(out,V(N1(base)),V(N2(base)),V(N2(base)-add_2),V(tet_top))
                end
            end
        else
            add_pri!(out,V(d1),V(d2),V(d3),V(d1+4),V(d2+4),V(d3+4))
        end
        if n_inner[1]<0
            add_pyr!(out,V(d1),V(d1+4),V(d3+4),V(d3),V(d0))
        else
            add=n_inner[1]>3 ? -4 : 4
            add_tet!(out,V(n_inner[1]),V(n_inner[1]+add),V(n_inner[2]),V(d0))
            add_tet!(out,V(n_inner[1]),V(n_inner[2]),V(n_inner[2]-add),V(d0))
        end
        return true
    elseif N1(d1)>=0 || N1(d2)>=0
        base=-1;top=-2
        if N2(d1)>=0 && N1(d2)==N2(d1)
            if N1(4)==N1(d2) || N2(4)==N1(d2)
                top=N1(d2);base=5
            elseif N1(5)==N1(d2) || N2(5)==N1(d2)
                top=N1(d2);base=4
            end
        end
        if base<0 && N2(d2)>=0
            if N2(d2)==d3+4 && (N1(5)==d3+4 || N2(5)==d3+4)
                top=d3+4;base=4
            elseif N2(d2)==d3 && (N1(4)==d3 || N2(4)==d3)
                top=d3;base=5
            end
        end
        if base<0 && N1(d1)>=0
            if N1(d1)==d1+4 && (N1(5)==d1+4 || N2(5)==d1+4)
                top=d1+4;base=4
            elseif N1(d1)==d1 && (N1(4)==d1 || N2(4)==d1)
                top=d1;base=5
            end
        end
        if base<0 && N1(d1)>=0 && N1(d2)>=0 && N1(d2)!=N2(d1)
            top=N1(d1);base=top<4 ? 5 : 4
        end
        if base>=0 && top>=0
            if top==N1(d1) || top==N2(d2)
                add=base==4 ? 0 : 4
                if N1(base)<0
                    add_pyr!(out,V(0+add),V(1+add),V(2+add),V(3+add),V(top))
                else
                    add_tet!(out,V(N1(base)),V((N1(base)-add+1)%4+add),V(N2(base)),V(top))
                    add_tet!(out,V(N1(base)),V(N2(base)),V((N1(base)-add+3)%4+add),V(top))
                end
                base2=-1
                if base==4 && (N1(5)==top || N2(5)==top)
                    base2=top==d1+4 ? d2 : d1
                elseif base==5 && (N1(4)==top || N2(4)==top)
                    base2=top==d1 ? d2 : d1
                else
                    base2=base==5 ? 4 : 5
                end
                if base2!=4 && base2!=5
                    add_tet!(out,V(top),V((top-4+add+2)%4),V((top-4+add+2)%4+4),V(d0))
                    if N1(base2)>=0
                        add_base2=N1(base2)<4 ? 4 : -4
                        add_tet!(out,V(N1(base2)),V(N1(base2)+add_base2),V(N2(base2)),V(top))
                        add_tet!(out,V(N1(base2)),V(N2(base2)),V(N2(base2)-add_base2),V(top))
                    else
                        add_pyr!(out,V(base2),V(base2+4),V((base2+1)%4+4),V((base2+1)%4),V(top))
                    end
                else
                    top2=(top-4+add+2)%4+add
                    if N1(base2)>=0
                        add_tet!(out,V(N1(base2)),V((N1(base2)-4+add+1)%4+4-add),V(N2(base2)),V(top2))
                        add_tet!(out,V(N1(base2)),V((N1(base2)-4+add+3)%4+4-add),V(N2(base2)),V(top2))
                    else
                        add_pyr!(out,V(0+4-add),V(1+4-add),V(2+4-add),V(3+4-add),V(top2))
                    end
                    tet_base=top2-add==d1 ? d1 : d2
                    if tet_base==d1
                        add_tet_base=N2(tet_base)>3 ? -4 : 4
                        add_tet!(out,V(N1(tet_base)),V(N2(tet_base)),V(N2(tet_base)+add_tet_base),V(top))
                    else
                        add_tet_base=N1(tet_base)>3 ? -4 : 4
                        add_tet!(out,V(N1(tet_base)),V(N2(tet_base)),V(N1(tet_base)+add_tet_base),V(top))
                    end
                end
                return true
            elseif top==d2 || top==d2+4
                add=base==4 ? 0 : 4
                if N1(base)>=0
                    add_tet!(out,V(N1(base)),V((N1(base)-add+1)%4+add),V(N2(base)),V(top))
                    add_tet!(out,V(N1(base)),V((N1(base)-add+3)%4+add),V(N2(base)),V(top))
                else
                    add_pyr!(out,V(0+add),V(1+add),V(2+add),V(3+add),V(top))
                end
                add_tet!(out,V(d0),V(d1),V(d1+4),V(top))
                add_tet!(out,V(d0),V(d3),V(d3+4),V(top))
                return true
            end
        end
    end
    return false
end

function full_hex_elems!(out::QTMesh,v,n1,n2,found_diags)
    V=i->v[i+1];N1=i->n1[i+1];N2=i->n2[i+1]
    if !found_diags
        add_hex!(out,V(0),V(1),V(2),V(3),V(4),V(5),V(6),V(7))
        return true
    end
    pyramid_top=-1;base_face_1=-1
    for s in 0:1
        p=s==0 ? 4 : 5
        if N1(p)>=0
            add=s==0 ? 0 : 4
            p1=N1(p)-add;p2=(N1(p)-add+3)%4
            if N1(p1)>=0 && N1(p2)>=0 && (N1(p)==N1(p1) || N1(p)==N2(p1)) &&
               (N1(p)==N1(p2) || N1(p)==N2(p2))
                pyramid_top=N1(p);base_face_1=s==0 ? 5 : 4
                break
            end
            p1=N2(p)-add;p2=(N2(p)-add+3)%4
            if N1(p1)>=0 && N1(p2)>=0 && (N2(p)==N1(p1) || N2(p)==N2(p1)) &&
               (N2(p)==N1(p2) || N2(p)==N2(p2))
                pyramid_top=N2(p);base_face_1=s==0 ? 5 : 4
                break
            end
        end
    end
    if pyramid_top>=0
        add=pyramid_top<4 ? 0 : 4
        base_face_2=(pyramid_top-add+1)%4
        base_face_3=(pyramid_top-add+2)%4
        for p in 0:2
            if p==0
                a=base_face_1==4 ? 0 : 4
                b=(0+a,1+a,2+a,3+a);bn1=N1(base_face_1);bn2=N2(base_face_1)
            elseif p==1
                b0=base_face_2;b=(b0,b0+4,(b0+1)%4+4,(b0+1)%4)
                # Gmsh reads n2 of base_face_3 here.
                bn1=N1(base_face_2);bn2=N2(base_face_3)
            else
                b0=base_face_3;b=(b0,b0+4,(b0+1)%4+4,(b0+1)%4)
                bn1=N1(base_face_3);bn2=N2(base_face_3)
            end
            if bn1<0
                add_pyr!(out,V(b[1]),V(b[2]),V(b[3]),V(b[4]),V(pyramid_top))
            elseif bn1==b[1] || bn2==b[1]
                add_tet!(out,V(b[1]),V(b[2]),V(b[3]),V(pyramid_top))
                add_tet!(out,V(b[1]),V(b[3]),V(b[4]),V(pyramid_top))
            else
                add_tet!(out,V(b[1]),V(b[2]),V(b[4]),V(pyramid_top))
                add_tet!(out,V(b[2]),V(b[3]),V(b[4]),V(pyramid_top))
            end
        end
        return true
    end
    prism_v=zeros(Int,2,6);prism_face=zeros(Int,2,2)
    t_prism=zeros(Int,2,3,4);p_prism=zeros(Int,2,5)
    for s in 0:1
        p_prism[s+1,:].=-2s-2
        for t in 0:2,w in 0:3
            t_prism[s+1,t+1,w+1]=-2*s*t*w-2
        end
    end
    for p_ind in 0:2
        p=p_ind==2 ? 4 : p_ind
        n_div1=-1;n_div2=-2
        if p<4
            if (N1(p)==p || N2(p)==p) && (N1((p+2)%4)==(p+2)%4+4 || N2((p+2)%4)==(p+2)%4+4)
                n_div1=p;n_div2=(p+1)%4+4
            end
            if (N1(p)==p+4 || N2(p)==p+4) && (N1((p+2)%4)==(p+2)%4 || N2((p+2)%4)==(p+2)%4)
                n_div1=p+4;n_div2=(p+1)%4
            end
        else
            if (N1(4)==0 || N2(4)==0) && (N1(5)==4 || N2(5)==4)
                n_div1=0;n_div2=2
            end
            if (N1(4)==1 || N2(4)==1) && (N1(5)==5 || N2(5)==5)
                n_div1=1;n_div2=3
            end
        end
        n_div1<0 && continue
        ext=zeros(Int,2,2,2);intern=zeros(Int,2)
        for s in 0:1
            intern[s+1]=-5-s*s-s
            for t in 0:1
                ext[s+1,t+1,1]=-1-t*t-t;ext[s+1,t+1,2]=-2-t*t-t
            end
        end
        if (n_div1>3 && n_div2<4) || (n_div1<4 && n_div2>3)
            p0=n_div1<4 ? p : (p+2)%4
            prism_v[1,:]=[p0,(p0+1)%4,(p0+1)%4+4,(p0+3)%4,(p0+2)%4,(p0+2)%4+4]
            prism_v[2,:]=[p0,p0+4,(p0+1)%4+4,(p0+3)%4,(p0+3)%4+4,(p0+2)%4+4]
            prism_face[1,:]=[4,(p0+1)%4];prism_face[2,:]=[(p0+3)%4,5]
        else
            p0=n_div1
            prism_v[1,:]=[p0,(p0+1)%4,(p0+2)%4,p0+4,(p0+1)%4+4,(p0+2)%4+4]
            prism_v[2,:]=[p0,(p0+3)%4,(p0+2)%4,p0+4,(p0+3)%4+4,(p0+2)%4+4]
            prism_face[1,:]=[p0,(p0+1)%4];prism_face[2,:]=[(p0+3)%4,(p0+2)%4]
        end
        for s in 1:2
            f0=prism_face[s,1];f1=prism_face[s,2]
            if N1(f0)==prism_v[s,2] || N2(f0)==prism_v[s,2]
                ext[s,1,1]=1;ext[s,1,2]=3
            elseif N1(f0)==prism_v[s,5] || N2(f0)==prism_v[s,5]
                ext[s,1,1]=4;ext[s,1,2]=0
            end
            if N1(f1)==prism_v[s,2] || N2(f1)==prism_v[s,2]
                ext[s,2,1]=1;ext[s,2,2]=5
            elseif N1(f1)==prism_v[s,5] || N2(f1)==prism_v[s,5]
                ext[s,2,1]=4;ext[s,2,2]=2
            end
        end
        if ext[1,1,1]>=0 && ext[1,2,1]!=ext[1,1,1]
            intern[1]=ext[1,1,2];intern[2]=ext[1,1,2]==0 ? 5 : 2
        elseif ext[1,2,1]>=0 && ext[1,1,1]!=ext[1,2,1]
            intern[1]=ext[1,2,2];intern[2]=ext[1,2,2]==2 ? 3 : 0
        end
        if (ext[2,1,1]!=ext[2,2,1] || (ext[2,1,1]<0 && ext[2,2,1]<0)) &&
           intern[1]>=0 && intern[1]!=ext[2,1,2] && intern[1]!=ext[2,2,2] &&
           intern[2]!=ext[2,1,2] && intern[2]!=ext[2,2,2]
            continue
        elseif intern[1]<0 && ext[2,1,1]>=0 && ext[2,2,1]!=ext[2,1,1]
            intern[1]=ext[2,1,2];intern[2]=ext[2,1,2]==0 ? 5 : 2
        elseif intern[1]<0 && ext[2,2,1]>=0 && ext[2,1,1]!=ext[2,2,1]
            intern[1]=ext[2,2,2];intern[2]=ext[2,2,2]==2 ? 3 : 0
        end
        if intern[1]>=0 && ((ext[1,1,1]<0 && ext[1,2,1]<0) || (ext[2,1,1]<0 && ext[2,2,1]<0))
            continue
        end
        for s in 1:2
            if ext[s,1,1]<0 && ext[s,2,1]<0
                continue
            elseif ext[s,1,1]>=0 && ext[s,1,1]==ext[s,2,1]
                add=ext[s,1,1]<3 ? 3 : -3
                t_prism[s,1,:]=[ext[s,1,1],ext[s,1,2],ext[s,2,2],ext[s,1,1]+add]
                if intern[1]<0
                    p_prism[s,:]=[0,2,5,3,ext[s,1,1]]
                else
                    v_tmp1=(intern[1]==0 || intern[2]==0) ? 2 : 0
                    v_tmp2=(intern[1]==0 || intern[2]==0) ? 3 : 5
                    t_prism[s,2,:]=[intern[1],intern[2],v_tmp1,ext[s,1,1]]
                    t_prism[s,3,:]=[intern[1],intern[2],v_tmp2,ext[s,1,1]]
                end
            elseif intern[1]>=0
                p_tet_start=ext[s,1,1]>=0 ? ext[s,1,1] : ext[s,2,1]
                p_pyr_top=p_tet_start==ext[s,1,1] ? ext[s,1,2] : ext[s,2,2]
                p_tet_top=p_tet_start==ext[s,1,1] ? (ext[s,1,2]==3 ? 0 : 3) :
                    (ext[s,2,2]==5 ? 2 : 5)
                t_prism[s,1,:]=[p_tet_start,intern[1],intern[2],p_tet_top]
                if p_tet_start==ext[s,1,1]
                    if ext[s,2,1]<0
                        p_prism[s,:]=[1,2,5,4,p_pyr_top]
                    else
                        t_prism[s,2,:]=[ext[s,2,1],ext[s,2,2],ext[s,2,1]==4 ? 1 : 4,p_pyr_top]
                        t_prism[s,3,:]=[ext[s,2,1],ext[s,2,2],ext[s,2,1]==4 ? 5 : 2,p_pyr_top]
                    end
                else
                    if ext[s,1,1]<0
                        p_prism[s,:]=[0,1,4,3,p_pyr_top]
                    else
                        t_prism[s,2,:]=[ext[s,1,1],ext[s,1,2],ext[s,1,1]==4 ? 3 : 4,p_pyr_top]
                        t_prism[s,3,:]=[ext[s,1,1],ext[s,1,2],ext[s,1,1]==4 ? 1 : 0,p_pyr_top]
                    end
                end
            end
            for t in 1:3,w in 1:4
                t_prism[s,t,w]>=0 && (t_prism[s,t,w]=prism_v[s,t_prism[s,t,w]+1])
            end
            for t in 1:5
                p_prism[s,t]>=0 && (p_prism[s,t]=prism_v[s,p_prism[s,t]+1])
            end
        end
        for s in 1:2
            if t_prism[s,1,1]<0 && p_prism[s,1]<0
                add_pri!(out,(V(prism_v[s,i]) for i in 1:6)...)
            else
                for t in 1:3
                    t_prism[s,t,1]>=0 && add_tet!(out,(V(t_prism[s,t,w]) for w in 1:4)...)
                end
                p_prism[s,1]>=0 && add_pyr!(out,(V(p_prism[s,i]) for i in 1:5)...)
            end
        end
        return true
    end
    return false
end

function hex_pri!(out::QTMesh,st::QTState,v,j,k,source,problems_new,bnd_elem,is_addverts)
    V=i->v[i+1]
    dup=Int[]
    for i in 0:3
        V(i)==V(i+4) && push!(dup,i)
    end
    m=length(dup)
    is_problem=!is_addverts && _is_problem(st.problems,source,j,k)
    found_diags=false
    n1=fill(-3,6);n2=fill(-4,6)
    E=st.quadToTri_edges
    for p in 0:5
        n1[p+1]=-p*p-p-1;n2[p+1]=-p*p-p-2
        if p<4
            if edge_exists(V(p),V((p+1)%4+4),E)
                n1[p+1]=p;n2[p+1]=(p+1)%4+4;found_diags=true
            elseif edge_exists(V((p+1)%4),V(p+4),E)
                n1[p+1]=p+4;n2[p+1]=(p+1)%4;found_diags=true
            end
        else
            add=p==4 ? 0 : 4
            if edge_exists(V(0+add),V(2+add),E)
                n1[p+1]=0+add;n2[p+1]=2+add;found_diags=true
            elseif edge_exists(V(1+add),V(3+add),E)
                n1[p+1]=1+add;n2[p+1]=3+add;found_diags=true
            end
        end
    end
    m>2 && return
    if is_addverts && (found_diags || m==1)
        mesh_with_internal_vertex!(out,v,n1,n2)
        return
    elseif m==2 && !is_problem
        two_pt_degen_hex_elems!(out,v,dup,n1,n2,found_diags) && return
    elseif m==1 && !is_problem
        one_pt_degen_hex_elems!(out,v,dup,n1,n2,found_diags) && return
    elseif !is_problem
        full_hex_elems!(out,v,n1,n2,found_diags) && return
    end
    if !is_problem
        push!(out.errors,"In QuadToTriHexPri(), Extruded hexahedron needs subdivision, but cannot  be divided without internal vertex, and was not previously detected as such.  This is a bug. Please Report.")
        make_centroid_vertex!(out,v)
        add_problem!(problems_new,source,j,k)
    end
    mesh_with_internal_vertex!(out,v,n1,n2)
    return
end

function add_body_centered_vertices!(out::QTMesh,st::QTState,r::QTRegion,c::CategorizedSourceElements)
    nlayer=length(r.layers)
    if !r.addverts
        for elem in sort!(collect(keys(st.problems)))
            for (j,k) in sort!(collect(st.problems[elem]))
                make_centroid_vertex!(out,extruded_vertices(r,elem,j,k))
            end
        end
        return
    end
    lat_size=length(st.lat_tri_diags)
    E=st.quadToTri_edges
    for t in 0:1,s in 0:2
        set_elems=t==0 ? (s==0 ? c.three_bnd_pt_tri : s==1 ? c.other_bnd_tri : c.internal_tri) :
                         (s==0 ? c.four_bnd_pt_quad : s==1 ? c.other_bnd_quad : c.internal_quad)
        for i in sort!(collect(set_elems))
            elem=(t,i);elem_size=t==0 ? 3 : 4
            found_diags=false;degen_hex=false
            v3=extruded_vertices(r,elem,0,0);V=x->v3[x+1]
            for p in 0:elem_size-1
                if t==1 && V(p)==V(p+elem_size) &&
                   V((p+1)%elem_size)!=V((p+1)%elem_size+elem_size) &&
                   V((p+elem_size-1)%elem_size)!=V((p+elem_size-1)%elem_size+elem_size)
                    degen_hex=true;break
                end
                (lat_size==0 || s>=2) && continue
                if V(p)==V(p+elem_size) || V((p+1)%elem_size)==V((p+1)%elem_size+elem_size)
                    continue
                elseif edge_exists(V(p),V((p+1)%elem_size+elem_size),E)
                    found_diags=true
                elseif edge_exists(V(p+elem_size),V((p+1)%elem_size),E)
                    found_diags=true
                end
            end
            (t==0 && (!found_diags || s==2)) && continue
            if (s<2 && found_diags) || degen_hex
                j_start=0;k_start=0
            else
                j_start=nlayer-1;k_start=r.layers[j_start+1]-1
            end
            # Gmsh restarts k at k_start in every later layer too.
            for j in j_start:nlayer-1,k in k_start:r.layers[j+1]-1
                make_centroid_vertex!(out,extruded_vertices(r,elem,j,k))
            end
        end
    end
    return
end

function create_elements!(out::QTMesh,st::QTState,r::QTRegion,c::CategorizedSourceElements)
    problems_new=Dict{ElemId,Set{Tuple{Int,Int}}}()
    add_body_centered_vertices!(out,st,r,c)
    nlayer=length(r.layers)
    for s in 0:2
        set_elems=s==0 ? c.three_bnd_pt_tri : s==1 ? c.other_bnd_tri : c.internal_tri
        for i in sort!(collect(set_elems)),j in 0:nlayer-1,k in 0:r.layers[j+1]-1
            pri_pyr_tet!(out,st,extruded_vertices(r,(0,i),j,k),j,k,(0,i),problems_new,s<2,r.addverts)
        end
    end
    for s in 0:2
        set_elems=s==0 ? c.four_bnd_pt_quad : s==1 ? c.other_bnd_quad : c.internal_quad
        for i in sort!(collect(set_elems)),j in 0:nlayer-1,k in 0:r.layers[j+1]-1
            hex_pri!(out,st,extruded_vertices(r,(1,i),j,k),j,k,(1,i),problems_new,s<2,r.addverts)
        end
    end
    for (elem,pairs) in problems_new
        union!(get!(()->Set{Tuple{Int,Int}}(),st.problems,elem),pairs)
    end
    return
end
