# ── MeshAdapt (Mesh.Algorithm 1) on a built-in plane ───────────────────────
#
# `meshGenerator` with algorithm 1: the same PolyMesh initial mesh and BDS
# mesh as Delaunay and Frontal-Delaunay, then `refineMeshBDS` twice (first
# with the nodal size field computed from the 1-D mesh) and BDS2GMSH. There is
# no Bowyer-Watson insertion and no Laplace smoothing afterwards.

function _gmsh_native_meshadapt_admitted(m::GeoModel,tag::Int)
    _surface_type(m,tag)===:plane || return false
    geometry=get(m.surface_geometry,tag,nothing)
    geometry!==nothing && hasproperty(geometry,:occ) && return false
    _model_surface_recombined(m,tag) && return false
    any(pair->first(pair)[1]==2 && (first(pair)[2]==tag ||
        abs(Int(last(pair).master_entity))==tag),m.periodic) && return false
    return true
end

struct _GmshMeshAdaptBackground{B}
    background::B
end
(background::_GmshMeshAdaptBackground)(u,v,x,y,z)=background.background(x,y,z)

# `gmshFace::closestPoint` of a plane: the normal equations of the mean
# plane's basis.
struct _GmshPlaneProjection
    frame::_GmshNativePlaneFrame
end
function (projection::_GmshPlaneProjection)(x)
    f=projection.frame
    B=(x[1]-f.origin[1],x[2]-f.origin[2],x[3]-f.origin[3])
    M=((f.first[1],f.second[1]),(f.first[2],f.second[2]),(f.first[3],f.second[3]))
    BN=ntuple(i->begin
        value=0.0
        for k in 1:3
            value+=B[k]*M[k][i]
        end
        value
    end,2)
    MN=ntuple(i->ntuple(j->begin
        value=0.0
        for k in 1:3
            value+=M[k][i]*M[k][j]
        end
        value
    end,2),2)
    return Mesh2D._gmsh_sys2x2(MN[1][1],MN[1][2],MN[2][1],MN[2][2],BN[1],BN[2])
end

function _gmsh_native_meshadapt_surface(m::GeoModel,tag::Int,options,input,
                                        high_water,caller;size_field=nothing,
                                        element_high_water=UInt64(0),forced=nothing,
                                        _periodic_positions=nothing)
    authority,limit=_gmsh_native_surface_lower_input(m,options,input,high_water,caller)
    frame=_gmsh_native_surface_frame(m,tag,options,caller)
    next_tag=Ref{UInt64}(limit)
    source=_gmsh_native_surface_source(m,tag,authority,frame,caller;
        ctx_lc=options.ctx_lc,forced,next_tag,_periodic_positions)
    limit=next_tag[]
    (;inside,p,row_of,ids,interior_faces,constrained,owners)=
        _gmsh_native_bds_input(source,tag)
    lift=_GmshNativeSurfaceLift(frame)
    bds=Mesh2D._gmsh_bds_build(p,ids,interior_faces,constrained;lift,owners)
    temporary=count(!,bds.face_deleted)
    background=_ModelAlgorithmBackgroundSize(m,2,tag,options,size_field,caller)
    embedded_points,_=_model_surface_embedding_tags(m,tag,caller)
    embedded_rows=Set{Int}(source.point_positions[point] for point in embedded_points)
    # The characteristic lengths `meshGenerator` gives the model points; the
    # nodal size field below replaces them but on embedded points.
    for v in eachindex(p.x)
        ids[v]<0 && continue
        row=row_of[v];world=source.world[row];owner=source.owners[row]
        value=_ModelAlgorithmBackgroundSize(m,Int(owner[1]),Int(owner[2]),
            options,size_field,caller)(world...)
        bds.lcbgm[v]=value;bds.lc[v]=value
    end
    plane_metric=(Mesh2D._gmsh_dot3(frame.first,frame.first),
                  Mesh2D._gmsh_dot3(frame.second,frame.first),
                  Mesh2D._gmsh_dot3(frame.second,frame.second))
    Mesh2D._gmsh_bds_delaunayize!(bds,(u,v)->plane_metric)
    for v in eachindex(p.x)
        ids[v]>=0 && Int(row_of[v]) in embedded_rows || continue
        bds.pdeg[v]=2;bds.ptag[v]=-1
    end
    extend=_resolved_extend_from_boundary(m,2,tag)
    face=Mesh2D._GmshMeshAdaptFace(lift,_GmshMeshAdaptBackground(background),true,
        extend>0 || extend==-2)
    projection=_GmshPlaneProjection(frame)
    steps=10
    Mesh2D._gmsh_bds_refine!(bds,face,steps,true,projection)
    Mesh2D._gmsh_bds_refine!(bds,face,steps,false,projection)
    # BDS2GMSH: cells in triangle-vector order; a BDS point without a model
    # node becomes a face node when a cell first names it.
    npoly=length(p.x)
    output=zeros(Int32,length(bds.id))
    coords=NTuple{3,Float64}[];tags=UInt64[];node_owners=Tuple{Int,Int32}[]
    parameters=Union{Nothing,Vector{Float64}}[];retained=Bool[]
    for v in 1:npoly
        ids[v]<0 && continue
        row=row_of[v]
        push!(coords,source.world[row]);push!(tags,source.tags[row])
        push!(node_owners,source.owners[row]);push!(parameters,source.parameters[row])
        push!(retained,!(row in source.forced_rows))
        output[v]=Int32(length(coords))
    end
    cells=NTuple{3,Int32}[]
    for nodes in Mesh2D._gmsh_bds_faces(bds)
        for q in nodes
            output[q]!=0 && continue
            push!(coords,(bds.X[q],bds.Y[q],bds.Z[q]))
            limit+=UInt64(1);push!(tags,limit)
            push!(node_owners,(2,Int32(tag)));push!(parameters,[bds.u[q],bds.v[q]])
            push!(retained,false)
            output[q]=Int32(length(coords))
        end
        push!(cells,(output[nodes[1]],output[nodes[2]],output[nodes[3]]))
    end
    # Nodes no cell uses (an unused model node cannot occur on a recovered
    # boundary) are dropped as `_gmsh_native_surface_export` drops them.
    used=falses(length(coords))
    for c in cells,q in c
        used[q]=true
    end
    keep=findall(used);mapping=zeros(Int32,length(coords))
    mapping[keep]=Int32.(1:length(keep))
    matrix=Matrix{Float64}(undef,3,length(keep))
    for (k,q) in enumerate(keep)
        matrix[1,k]=coords[q][1];matrix[2,k]=coords[q][2];matrix[3,k]=coords[q][3]
    end
    triangles=Matrix{Int32}(undef,3,length(cells))
    for (column,c) in enumerate(cells),slot in 1:3
        triangles[slot,column]=mapping[c[slot]]
    end
    births_start=element_high_water+UInt64(temporary)
    births=UInt64[births_start+UInt64(k) for k in 1:length(cells)]
    mesh=Mesh(matrix;tris=triangles)
    get(m.meshing.reverse,(2,tag),false) &&
        (mesh=_reversed_surface_mesh(mesh,caller,tag))
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid MeshAdapt Surface[$tag] — "*join(diagnostic.messages,"; ")))
    return _ModelSurfaceProduction(mesh,tags[keep],node_owners[keep],parameters[keep],
        BitVector(retained[keep]),limit,births,births_start+UInt64(length(cells)),
        Int32.(1:length(keep)),nothing,Dict{UInt64,NTuple{3,Float64}}(),
        Tuple{Symbol,String}[],_gmsh_native_surface_curve_chains(source))
end
