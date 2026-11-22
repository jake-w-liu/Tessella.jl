# The surfaces around a QuadToTri region and meshQuadToTriRegion
# (meshGFaceExtruded.cpp createQuaTri, QuadTriExtruded2D.cpp
# MeshQuadToTriTopSurface / MeshQuadToTriTopUnstructured,
# QuadTriExtruded3D.cpp QuadToTriLateralRemesh / meshQuadToTriRegion).

# createQuaTri on v = (v0@l, v1@l, v0@l+1, v1@l+1).
function create_qua_tri(v,flag::Int,recombine::Bool,constrained)
    cells=Vector{Tuple{Vararg{Int}}}()
    if v[1]==v[2] || v[2]==v[4]
        push!(cells,(v[1],v[4],v[3]))
    elseif v[1]==v[3] || v[3]==v[4]
        push!(cells,(v[1],v[2],v[4]))
    elseif v[1]==v[4] || v[2]==v[3]
        # "Incoherent extruded quadrangle"
    elseif (recombine && flag!=2) || flag==1
        if constrained===nothing
            push!(cells,(v[1],v[2],v[4],v[3]))
        elseif edge_exists(v[2],v[3],constrained)
            push!(cells,(v[3],v[2],v[1]),(v[3],v[4],v[2]))
        elseif edge_exists(v[1],v[4],constrained)
            push!(cells,(v[3],v[4],v[1]),(v[1],v[4],v[2]))
        else
            push!(cells,(v[1],v[2],v[4],v[3]))
        end
    elseif constrained===nothing
        push!(cells,(v[1],v[2],v[4]),(v[1],v[4],v[3]))
    elseif edge_exists(v[2],v[3],constrained)
        push!(cells,(v[3],v[2],v[1]),(v[3],v[4],v[2]))
    else
        push!(cells,(v[3],v[4],v[1]),(v[1],v[4],v[2]))
    end
    return cells
end

# MeshExtrudedSurface of every lateral, cells per source line.
function lateral_meshes(r::QTRegion,constrained)
    out=Vector{Vector{Vector{Tuple{Vararg{Int}}}}}()
    for (s,lines) in enumerate(r.lateral_lines)
        per_line=Vector{Vector{Tuple{Vararg{Int}}}}()
        for (v0,v1) in lines
            cells=Vector{Tuple{Vararg{Int}}}()
            for level in 0:sum(r.layers)-1
                v=(r.vid(v0,level),r.vid(v1,level),r.vid(v0,level+1),r.vid(v1,level+1))
                append!(cells,create_qua_tri(v,r.lateral_flag[s],r.lateral_recombine[s],constrained))
            end
            push!(per_line,cells)
        end
        push!(out,per_line)
    end
    return out
end

# MeshQuadToTriTopSurface: the top triangles of the source quadrangles and
# their diagonals.
function top_surface(r::QTRegion)
    top=sum(r.layers)
    cells=NTuple{3,Int}[];diagonals=Edge[]
    # copyMesh: the source triangles first, then the divided quadrangles.
    for t in r.tris
        push!(cells,map(v->r.vid(v,top),t))
    end
    for q in r.quads
        verts=[r.vid(v,top) for v in q]
        V=i->verts[i+1]
        if !r.addverts
            bnd=[q[p+1] in r.boundary for p in 0:3]
            count_bnd=count(bnd)
            low_index=-1
            if count_bnd in (1,2,3)
                for p in 0:3
                    if !bnd[p+1] && V(p)!=q[p+1]
                        if low_index<0 || V(p)<V(low_index)
                            low_index=p
                        end
                    end
                end
                low_index<0 && (low_index=0)
            else
                low_index=0
                for p in 1:3
                    V(p)<V(low_index) && (low_index=p)
                end
            end
            push!(cells,(V(low_index),V((low_index+1)%4),V((low_index+2)%4)),
                        (V(low_index),V((low_index+2)%4),V((low_index+3)%4)))
            push!(diagonals,_ekey(V(low_index),V((low_index+2)%4)))
        elseif !r.structured_root
            # MeshQuadToTriTopUnstructured, with Gmsh's int edge lengths.
            dist(a,b)=(ca=r.coord(a);cb=r.coord(b);hypot(ca[1]-cb[1],ca[2]-cb[2],ca[3]-cb[3]))
            mag_sq_ave=0.0
            for p in 0:3
                d_leg=trunc(Int,dist(V(p),V((p+1)%4)))
                mag_sq_ave+=d_leg*d_leg
            end
            mag_sq_ave/=4
            d1=dist(V(0),V(2));d2=dist(V(1),V(3))
            if abs(d1*d1-mag_sq_ave)<=abs(d2*d2-mag_sq_ave)
                push!(cells,(V(0),V(1),V(2)),(V(0),V(2),V(3)))
                push!(diagonals,_ekey(V(0),V(2)))
            else
                push!(cells,(V(1),V(2),V(3)),(V(1),V(3),V(0)))
                push!(diagonals,_ekey(V(1),V(3)))
            end
        else
            push!(cells,(V(0),V(2),V(3)),(V(0),V(1),V(2)))
            push!(diagonals,_ekey(V(0),V(2)))
        end
    end
    return cells,diagonals
end

# meshQuadToTriRegion for a lone QuadToTri region: its volume cells, the
# remeshed laterals, the top triangles and the body-centred vertices.
function qt_mesh_region(r::QTRegion)
    st=QTState()
    c=categorize(r)
    first_laterals=lateral_meshes(r,nothing)
    top_cells,top_diagonals=top_surface(r)
    edge_generator!(st,r,c,first_laterals,top_diagonals)
    out=QTMesh(r.coord,r.next_id)
    create_elements!(out,st,r,c)
    laterals=lateral_meshes(r,st.quadToTri_edges)
    return (mesh=out,laterals=laterals,top=top_cells,state=st)
end
