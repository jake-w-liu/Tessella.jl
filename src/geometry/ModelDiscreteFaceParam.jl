# SPDX-FileCopyrightText: 1997-2025 C. Geuzaine, J.-F. Remacle
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: MIT AND GPL-2.0-or-later
# Contains Gmsh 4.15.2 adaptations: geo/discreteFace.cpp (createGeometry,
# _createGeometryFromSTL, param::checkPlanar, point, closestPoint, firstDer,
# normal, intersectionWithCircle), geo/GModelParametrize.cpp
# (computeParametrization), geo/MEdge.cpp (SortEdgeConsecutive),
# numeric/Numeric.cpp (signedDistancePointTriangle, inv3x3, sys2x2).
# Julia adaptation/changes: Tessella contributors, 2026.
# Licenses apply to their respective portions; see LICENSES/README.md and NOTICES.md#gmsh.

# ---- discreteFace parametrization ------------------------------------------
#
# `discreteFace::createGeometry` maps the surface's triangulation harmonically
# onto the unit disk (`computeParametrization`, Floater's mean-value weights,
# the longest boundary loop fixed on the circle) and keeps the result apart
# from the mesh in `discreteFace::_param`: its own copy of the triangles in
# parameter space (`t2d`) and in space (`t3d`), a parameter-space locator
# (`oct`) and a spatial one (`rtree3d`). A planar patch replaces the disk by
# its mean-plane coordinates (`param::checkPlanar`). Every geometric query of
# the face — point, derivatives, normal, closest point, bounding box — reads
# this record, so the face mesh itself can be deleted and regenerated.
#
# Tessella stores the record as `m.surface_geometry[tag]` with the property
# `discrete_surface=true`, beside the discrete Curve snapshots that live in
# `m.curve_geometry`.

# Uniform bucket grid over triangle bounding boxes: the parameter-space
# `MElementOctree` (2-D, `dims[3]==1`) and the spatial `rtree3d`.
struct _DiscreteFaceGrid
    lo::NTuple{3,Float64}
    scale::NTuple{3,Float64}
    dims::NTuple{3,Int}
    start::Vector{Int32}
    items::Vector{Int32}
end

struct _DiscreteFaceParam
    node_tags::Vector{UInt64}
    # `stl_vertices_uv`: the harmonic map before the planar simplification.
    stl_uv::Vector{NTuple{2,Float64}}
    # `_param.v2d` (parameter vertices actually meshed) and `_param.v3d`.
    uv::Vector{NTuple{2,Float64}}
    xyz::Vector{NTuple{3,Float64}}
    triangles::Vector{NTuple{3,Int32}}
    bbox::NTuple{6,Float64}
    # `parBounds`: `param()`'s (-1,1) defaults unless `checkPlanar` set them.
    bounds::NTuple{4,Float64}
    planar::Bool
    grid2::_DiscreteFaceGrid
    grid3::_DiscreteFaceGrid
end

function _discrete_face_param(m::GeoModel,tag::Integer)
    geometry=get(m.surface_geometry,Int(tag),nothing)
    geometry!==nothing && hasproperty(geometry,:discrete_surface) &&
        geometry.discrete_surface || return nothing
    haskey(m.discrete,(2,Int(tag))) || return nothing
    return geometry.param::_DiscreteFaceParam
end

# ---- bucket grids ----------------------------------------------------------

function _discrete_face_grid(boxes::Vector{NTuple{6,Float64}},dimension::Int)
    lo=[Inf,Inf,Inf];hi=[-Inf,-Inf,-Inf]
    for box in boxes,axis in 1:3
        lo[axis]=min(lo[axis],box[axis]);hi[axis]=max(hi[axis],box[axis+3])
    end
    count=max(length(boxes),1)
    per_axis=max(1,round(Int,dimension==2 ? sqrt(count) : cbrt(count)))
    dims=ntuple(axis->axis<=dimension ? per_axis : 1,3)
    scale=ntuple(axis->begin
        span=hi[axis]-lo[axis]
        isfinite(span) && span>0 ? dims[axis]/span : 0.0
    end,3)
    origin=ntuple(axis->isfinite(lo[axis]) ? lo[axis] : 0.0,3)
    cell(axis,value)=clamp(floor(Int,(value-origin[axis])*scale[axis])+1,1,dims[axis])
    ncell=dims[1]*dims[2]*dims[3]
    counts=zeros(Int32,ncell+1)
    linear(i,j,k)=i+(j-1)*dims[1]+(k-1)*dims[1]*dims[2]
    for box in boxes
        for k in cell(3,box[3]):cell(3,box[6]),j in cell(2,box[2]):cell(2,box[5]),
                i in cell(1,box[1]):cell(1,box[4])
            counts[linear(i,j,k)+1]+=1
        end
    end
    start=Vector{Int32}(undef,ncell+1);start[1]=1
    for c in 1:ncell
        start[c+1]=start[c]+counts[c+1]
    end
    items=Vector{Int32}(undef,start[end]-1)
    fill_at=copy(start)
    for (index,box) in enumerate(boxes)
        for k in cell(3,box[3]):cell(3,box[6]),j in cell(2,box[2]):cell(2,box[5]),
                i in cell(1,box[1]):cell(1,box[4])
            c=linear(i,j,k)
            items[fill_at[c]]=Int32(index);fill_at[c]+=1
        end
    end
    return _DiscreteFaceGrid(origin,scale,dims,start,items)
end

@inline _discrete_face_grid_cell(grid::_DiscreteFaceGrid,axis,value)=
    clamp(floor(Int,(value-grid.lo[axis])*grid.scale[axis])+1,1,grid.dims[axis])

@inline function _discrete_face_grid_linear(grid::_DiscreteFaceGrid,i,j,k)
    return i+(j-1)*grid.dims[1]+(k-1)*grid.dims[1]*grid.dims[2]
end

# ---- computeParametrization -----------------------------------------------

# `SortEdgeConsecutive` over `MEdge`s given as (getVertex(0), getVertex(1));
# the adjacency map is ordered by node number (`MVertexPtrLessThan`).
function _discrete_face_sort_edge_consecutive(edges::Vector{NTuple{2,UInt64}})
    chains=Vector{UInt64}[]
    isempty(edges) && return chains
    adjacency=Dict{UInt64,Tuple{UInt64,UInt64}}()   # 0 is "nullptr"
    for (v0,v1) in edges
        for (a,b) in ((v0,v1),(v1,v0))
            pair=get(adjacency,a,nothing)
            if pair===nothing
                adjacency[a]=(b,UInt64(0))
            elseif pair[2]==0
                adjacency[a]=(pair[1],b)
            else
                return nothing
            end
        end
    end
    order=sort!(collect(keys(adjacency)))
    while !isempty(adjacency)
        start=UInt64(0)
        for key in order
            pair=get(adjacency,key,nothing)
            pair===nothing && continue
            start==0 && (start=key)
            if pair[2]==0
                start=key
                break
            end
        end
        pair=adjacency[start]
        previous=pair[2]==start ? pair[1] : pair[2]
        current=start
        chain=UInt64[]
        while true
            push!(chain,current)
            pair=get(adjacency,current,nothing)
            pair===nothing && return nothing
            delete!(adjacency,current)
            temp=current
            if pair[1]==previous
                current=pair[2]
            elseif pair[2]==previous
                current=pair[1]
            else
                break
            end
            previous=temp
            current==start && push!(chain,current)
            (current!=start && current!=0) || break
        end
        length(chain)>2 && chain[end-1]==chain[end] && pop!(chain)
        push!(chains,chain)
    end
    return chains
end

@inline _discrete_face_distance(a,b)=
    sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2)

# `computeParametrization(triangles, ...)`: returns (nodes, uv, stl triangles)
# or `nothing` with the reason, as upstream's `Msg::Error` branches.
function _discrete_face_compute_parametrization(
        triangles::Vector{NTuple{3,UInt64}},coords::Dict{UInt64,NTuple{3,Float64}})
    isempty(triangles) && return nothing,"no triangles"
    node_index=Dict{UInt64,Int}();nodes=UInt64[]
    # std::map<MEdge, std::vector<MTriangle*>, MEdgeLessThan>: the key keeps
    # the vertex order of the first triangle that inserted it.
    edge_slot=Dict{NTuple{2,UInt64},Int}()
    edge_vertices=NTuple{2,UInt64}[];edge_tris=Vector{Int}[]
    for (t,tri) in enumerate(triangles)
        for j in 1:3
            v=tri[j]
            haskey(node_index,v) || (push!(nodes,v);node_index[v]=length(nodes))
            a=tri[j];b=tri[mod1(j+1,3)]
            key=a<b ? (a,b) : (b,a)
            slot=get(edge_slot,key,0)
            if slot==0
                push!(edge_vertices,(a,b));push!(edge_tris,Int[t])
                edge_slot[key]=length(edge_vertices)
            else
                push!(edge_tris[slot],t)
            end
        end
    end
    order=sortperm(edge_vertices;by=e->e[1]<e[2] ? e : (e[2],e[1]))
    boundary=NTuple{2,UInt64}[]
    for slot in order
        count=length(edge_tris[slot])
        if count==1
            push!(boundary,edge_vertices[slot])
        elseif count!=2
            return nothing,"Wrong topology of triangulation for parametrization: "*
                "one edge is incident to $count triangles"
        end
    end
    loops=_discrete_face_sort_edge_consecutive(boundary)
    loops===nothing && return nothing,"Wrong topology of boundary mesh for parametrization"
    (isempty(loops) || length(loops[1])<2) &&
        return nothing,"Invalid exterior boundary mesh for parametrization"
    point(v)=coords[v]
    loop=1;longest=0.0
    for (i,chain) in enumerate(loops)
        l=0.0
        for j in 2:length(chain)
            l+=_discrete_face_distance(point(chain[j]),point(chain[j-1]))
        end
        l>longest && (longest=l;loop=i)
    end
    reverse=true
    ref=(loops[loop][1],loops[loop][2])
    for tri in triangles
        for j in 1:3
            if tri[j]==ref[1] && tri[mod1(j+1,3)]==ref[2]
                reverse=false
                break
            end
        end
        reverse || break
    end
    # Upstream reverses `vs[0]`, not `vs[loop]`.
    reverse && reverse!(loops[1])
    n=length(nodes)
    u=zeros(n);v=zeros(n);bc=falses(n)
    outer=loops[loop]
    index=node_index[outer[1]]
    bc[index]=true;u[index]=1.0;v[index]=0.0
    current_length=0.0
    for i in 2:length(outer)-1
        current_length+=_discrete_face_distance(point(outer[i]),point(outer[i-1]))
        angle=2*pi*current_length/longest
        index=node_index[outer[i]]
        bc[index]=true;u[index]=cos(angle);v[index]=sin(angle)
    end
    rows=Int[];columns=Int[];values=Float64[]
    third(t,a,b)=begin
        tri=triangles[t]
        w=tri[1]
        (w==a || w==b) && (w=tri[2])
        (w==a || w==b) && (w=tri[3])
        w
    end
    nnz=0
    for slot in order
        e0,e1=edge_vertices[slot]
        for ij in 0:1
            v0=ij==0 ? e0 : e1;v1=ij==0 ? e1 : e0
            index0=node_index[v0]
            bc[index0] && continue
            index1=node_index[v1]
            p0=point(v0);p1=point(v1)
            e=(p1[1]-p0[1],p1[2]-p0[2],p1[3]-p0[3])
            ne=sqrt(e[1]*e[1]+e[2]*e[2]+e[3]*e[3])
            left=point(third(edge_tris[slot][1],v0,v1))
            a=(left[1]-p0[1],left[2]-p0[2],left[3]-p0[3])
            na=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
            thetaL=acos((a[1]*e[1]+a[2]*e[2]+a[3]*e[3])/(na*ne))
            thetaR=0.0
            if length(edge_tris[slot])==2
                right=point(third(edge_tris[slot][2],v0,v1))
                b=(right[1]-p0[1],right[2]-p0[2],right[3]-p0[3])
                nb=sqrt(b[1]*b[1]+b[2]*b[2]+b[3]*b[3])
                thetaR=acos((b[1]*e[1]+b[2]*e[2]+b[3]*e[3])/(nb*ne))
            end
            c=(tan(0.5*thetaL)+tan(0.5*thetaR))/ne
            push!(rows,index0);push!(columns,index1);push!(values,-c)
            push!(rows,index0);push!(columns,index0);push!(values,c)
            nnz+=2
        end
    end
    if nnz>0
        for i in 1:length(outer)-1
            row=node_index[outer[i]]
            push!(rows,row);push!(columns,row);push!(values,1.0)
        end
        rhs_u=zeros(n);rhs_v=zeros(n)
        for i in 1:length(outer)-1
            row=node_index[outer[i]]
            rhs_u[row]+=u[row];rhs_v[row]+=v[row]
        end
        solved=_discrete_face_linear_solve(n,rows,columns,values,rhs_u,rhs_v)
        solved===nothing && return nothing,"parametrization linear solve failed"
        u,v=solved
    end
    return (nodes,[(u[i],v[i]) for i in 1:n],
            [NTuple{3,Int32}(Int32(node_index[w]) for w in tri) for tri in triangles]),""
end

# Upstream solves with PETSc's direct LU (`-ksp_type preonly -pc_type lu`).
# Dense LU reproduces that solution to rounding for the patch sizes STL
# classification produces; large patches use restarted GMRES with a Jacobi
# preconditioner, iterated to a tolerance well below the mesher's.
function _discrete_face_linear_solve(n,rows,columns,values,rhs_u,rhs_v)
    if n<=2000
        A=zeros(n,n)
        for k in eachindex(values)
            A[rows[k],columns[k]]+=values[k]
        end
        factor=try
            LinearAlgebra.lu(A)
        catch error
            error isa InterruptException && rethrow()
            return nothing
        end
        LinearAlgebra.issuccess(factor) || return nothing
        return factor\rhs_u,factor\rhs_v
    end
    # CSR assembly with duplicate summation.
    perm=sortperm(collect(zip(rows,columns)))
    rowptr=zeros(Int,n+1);cols=Int[];vals=Float64[]
    last=(0,0)
    for k in perm
        key=(rows[k],columns[k])
        if key==last
            vals[end]+=values[k]
        else
            push!(cols,key[2]);push!(vals,values[k]);rowptr[key[1]+1]+=1
            last=key
        end
    end
    rowptr[1]=1
    for i in 1:n
        rowptr[i+1]+=rowptr[i]
    end
    diagonal=ones(n)
    for i in 1:n,k in rowptr[i]:rowptr[i+1]-1
        cols[k]==i && vals[k]!=0 && (diagonal[i]=vals[k])
    end
    multiply!(out,x)=begin
        for i in 1:n
            s=0.0
            for k in rowptr[i]:rowptr[i+1]-1
                s+=vals[k]*x[cols[k]]
            end
            out[i]=s
        end
        out
    end
    results=Vector{Float64}[]
    for b in (rhs_u,rhs_v)
        x=_discrete_face_gmres(multiply!,diagonal,b,n)
        x===nothing && return nothing
        push!(results,x)
    end
    return results[1],results[2]
end

function _discrete_face_gmres(multiply!,diagonal,b,n;restart=60,maxiter=4000,tol=1e-14)
    x=zeros(n);r=similar(x);w=similar(x)
    bnorm=LinearAlgebra.norm(b);bnorm==0 && return x
    V=zeros(n,restart+1);H=zeros(restart+1,restart);cs=zeros(restart);sn=zeros(restart)
    g=zeros(restart+1)
    iterations=0
    while iterations<maxiter
        multiply!(r,x);r.=(b.-r)./diagonal
        beta=LinearAlgebra.norm(r)
        beta<=tol*LinearAlgebra.norm(b./diagonal) && return x
        V[:,1].=r./beta;fill!(g,0.0);g[1]=beta
        k_used=0
        for k in 1:restart
            iterations+=1
            multiply!(w,view(V,:,k));w./=diagonal
            for j in 1:k
                H[j,k]=LinearAlgebra.dot(view(V,:,j),w);w.-=H[j,k].*view(V,:,j)
            end
            H[k+1,k]=LinearAlgebra.norm(w)
            H[k+1,k]>0 && (V[:,k+1].=w./H[k+1,k])
            for j in 1:k-1
                t=cs[j]*H[j,k]+sn[j]*H[j+1,k]
                H[j+1,k]=-sn[j]*H[j,k]+cs[j]*H[j+1,k];H[j,k]=t
            end
            d=hypot(H[k,k],H[k+1,k])
            cs[k]=H[k,k]/d;sn[k]=H[k+1,k]/d
            H[k,k]=d;H[k+1,k]=0.0
            g[k+1]=-sn[k]*g[k];g[k]=cs[k]*g[k]
            k_used=k
            abs(g[k+1])<=tol*LinearAlgebra.norm(b./diagonal) && break
        end
        y=zeros(k_used)
        for i in k_used:-1:1
            s=g[i]
            for j in i+1:k_used
                s-=H[i,j]*y[j]
            end
            y[i]=s/H[i,i]
        end
        for j in 1:k_used
            x.+=y[j].*view(V,:,j)
        end
    end
    multiply!(r,x)
    LinearAlgebra.norm(b.-r)<=1e-10*bnorm || return nothing
    return x
end

# ---- _createGeometryFromSTL ------------------------------------------------

# `discreteFace::param::checkPlanar`: the mean plane of the t3d corners (one
# sample per corner, as upstream collects them), the `matchMeshTolerance`
# test, then plane coordinates for every parameter vertex.
function _discrete_face_check_planar(xyz,triangles,tolerance)
    samples=NTuple{3,Float64}[]
    sizehint!(samples,3*length(triangles))
    lo=(Inf,Inf,Inf);hi=(-Inf,-Inf,-Inf)
    for tri in triangles,j in 1:3
        p=xyz[tri[j]]
        push!(samples,p)
        lo=min.(lo,p);hi=max.(hi,p)
    end
    length(samples)>=3 || return nothing
    frame=try
        _gmsh_native_plane_small_svd_frame(samples)
    catch error
        error isa InterruptException && rethrow()
        return nothing
    end
    center=zeros(3)
    for point in samples,axis in 1:3
        center[axis]+=point[axis]
    end
    for axis in 1:3;center[axis]/=length(samples);end
    a,b,c=frame.normal
    d=Mesh2D._gmsh_dot3((center[1],center[2],center[3]),frame.normal)
    diagonal=_discrete_face_distance(lo,hi)
    for p in samples
        F=a*p[1]+b*p[2]+c*p[3]-d
        abs(F)>tolerance*diagonal && return nothing
    end
    uv=Vector{NTuple{2,Float64}}(undef,length(xyz))
    origin=frame.origin
    for (index,p) in enumerate(xyz)
        t0=-(a*p[1]+b*p[2]+c*p[3]+d)/(a*a+b*b+c*c)
        projected=(p[1]+a*t0,p[2]+b*t0,p[3]+c*t0)
        offset=Mesh2D._gmsh_sub3(projected,origin)
        uv[index]=(Mesh2D._gmsh_dot3(offset,frame.first),
                   Mesh2D._gmsh_dot3(offset,frame.second))
    end
    return uv
end

function _discrete_face_build_param(nodes,stl_uv,xyz,triangles;tolerance=1e-6)
    lo=(Inf,Inf,Inf);hi=(-Inf,-Inf,-Inf)
    for p in xyz
        lo=min.(lo,p);hi=max.(hi,p)
    end
    planar_uv=_discrete_face_check_planar(xyz,triangles,tolerance)
    planar=planar_uv!==nothing
    uv=planar ? planar_uv : copy(stl_uv)
    bounds=(-1.0,1.0,-1.0,1.0)
    if planar
        umin=Inf;umax=-Inf;vmin=Inf;vmax=-Inf
        for tri in triangles,j in 1:3
            q=uv[tri[j]]
            umin=min(umin,q[1]);vmin=min(vmin,q[2])
            umax=max(umax,q[1]);vmax=max(vmax,q[2])
        end
        bounds=(umin,umax,vmin,vmax)
    end
    boxes2=NTuple{6,Float64}[]
    boxes3=NTuple{6,Float64}[]
    for tri in triangles
        a,b,c=uv[tri[1]],uv[tri[2]],uv[tri[3]]
        push!(boxes2,(min(a[1],b[1],c[1]),min(a[2],b[2],c[2]),0.0,
                      max(a[1],b[1],c[1]),max(a[2],b[2],c[2]),0.0))
        p,q,r=xyz[tri[1]],xyz[tri[2]],xyz[tri[3]]
        push!(boxes3,(min(p[1],q[1],r[1]),min(p[2],q[2],r[2]),min(p[3],q[3],r[3]),
                      max(p[1],q[1],r[1]),max(p[2],q[2],r[2]),max(p[3],q[3],r[3])))
    end
    return _DiscreteFaceParam(nodes,stl_uv,uv,xyz,triangles,(lo...,hi...),bounds,
        planar,_discrete_face_grid(boxes2,2),_discrete_face_grid(boxes3,3))
end

# `discreteFace::createGeometry` for one record; `nothing` when upstream
# leaves `_param` empty.
function _discrete_face_create_geometry(record::DiscreteEntity,
        coords::Dict{UInt64,NTuple{3,Float64}};tolerance=1e-6)
    triangles=NTuple{3,UInt64}[]
    for (index,msh) in enumerate(record.element_types)
        msh==2 || continue
        nodes=record.element_nodes[index]
        push!(triangles,(nodes[1],nodes[2],nodes[3]))
    end
    isempty(triangles) && return nothing
    for tri in triangles,node in tri
        haskey(coords,node) || return nothing
    end
    result,_=_discrete_face_compute_parametrization(triangles,coords)
    result===nothing && return nothing
    nodes,stl_uv,stl_triangles=result
    xyz=NTuple{3,Float64}[coords[node] for node in nodes]
    return _discrete_face_build_param(nodes,stl_uv,xyz,stl_triangles;tolerance)
end

# ---- queries -----------------------------------------------------------------

# `MYxyz2uvw` on a parameter triangle.
@inline function _discrete_face_local_uv(P::_DiscreteFaceParam,t,u,v)
    tri=P.triangles[t]
    p0=P.uv[tri[1]];p1=P.uv[tri[2]];p2=P.uv[tri[3]]
    m00=p1[1]-p0[1];m01=p2[1]-p0[1];m10=p1[2]-p0[2];m11=p2[2]-p0[2]
    r0=u-p0[1];r1=v-p0[2]
    det=m00*m11-m10*m01
    return (r0*m11-m01*r1)/det,(m00*r1-m10*r0)/det
end

# `_param.oct->find(u, v, 0, -1, true)`: a parameter triangle containing
# `(u,v)` within the reference-element tolerance, or 0.
function _discrete_face_locate(P::_DiscreteFaceParam,u::Float64,v::Float64)
    grid=P.grid2
    (isfinite(u) && isfinite(v)) || return 0
    i=_discrete_face_grid_cell(grid,1,u);j=_discrete_face_grid_cell(grid,2,v)
    c=_discrete_face_grid_linear(grid,i,j,1)
    tolerance=1e-6
    best=0;best_value=-Inf
    @inbounds for k in grid.start[c]:grid.start[c+1]-1
        t=grid.items[k]
        a,b=_discrete_face_local_uv(P,t,u,v)
        (isfinite(a) && isfinite(b)) || continue
        lowest=min(a,b,1.0-a-b)
        lowest>=0.0 && return Int(t)
        if a>=-tolerance && b>=-tolerance && a<=(1.0+tolerance)-b && lowest>best_value
            best=Int(t);best_value=lowest
        end
    end
    return best
end

# `discreteFace::point`: (x,y,z) and whether the triangle was found.
function _discrete_face_point(P::_DiscreteFaceParam,u::Float64,v::Float64)
    t=_discrete_face_locate(P,u,v)
    t==0 && return (1.e21,1.e21,1.e21),false
    a,b=_discrete_face_local_uv(P,t,u,v)
    tri=P.triangles[t]
    eval=(1.0-a-b,a,b)
    X=0.0;Y=0.0;Z=0.0
    for io in 1:3
        p=P.xyz[tri[io]]
        X+=p[1]*eval[io];Y+=p[2]*eval[io];Z+=p[3]*eval[io]
    end
    return (X,Y,Z),true
end

# `discreteFace::firstDer`.
function _discrete_face_first_derivatives(P::_DiscreteFaceParam,u::Float64,v::Float64)
    t=_discrete_face_locate(P,u,v)
    t==0 && return (1.0,0.0,0.0),(0.0,1.0,0.0)
    tri=P.triangles[t]
    x1=P.xyz[tri[1]];x2=P.xyz[tri[2]];x3=P.xyz[tri[3]]
    M3D=((x2[1]-x1[1],x3[1]-x1[1]),(x2[2]-x1[2],x3[2]-x1[2]),(x2[3]-x1[3],x3[3]-x1[3]))
    q1=P.uv[tri[1]];q2=P.uv[tri[2]];q3=P.uv[tri[3]]
    M2D=(((q3[2]-q1[2]),-(q3[1]-q1[1])),(-(q2[2]-q1[2]),(q2[1]-q1[1])))
    det=1.0/(M2D[1][1]*M2D[2][2]-M2D[2][1]*M2D[1][2])
    d(i,j)=begin
        s=0.0
        for k in 1:2
            s+=det*M3D[i][k]*M2D[k][j]
        end
        s
    end
    return (d(1,1),d(2,1),d(3,1)),(d(1,2),d(2,2),d(3,2))
end

# `discreteFace::normal`.
function _discrete_face_normal(P::_DiscreteFaceParam,u::Float64,v::Float64)
    t=_discrete_face_locate(P,u,v)
    t==0 && return (0.0,0.0,1.0)
    tri=P.triangles[t]
    x1=P.xyz[tri[1]];x2=P.xyz[tri[2]];x3=P.xyz[tri[3]]
    return Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(
        Mesh2D._gmsh_sub3(x2,x1),Mesh2D._gmsh_sub3(x3,x1)))
end

# `signedDistancePointTriangle`: (distance, closest point), or `nothing` for
# a degenerate triangle (upstream returns with `d` unset).
function _discrete_face_signed_distance(p1,p2,p3,p)
    t1=Mesh2D._gmsh_sub3(p2,p1);t2=Mesh2D._gmsh_sub3(p3,p1);t3=Mesh2D._gmsh_sub3(p3,p2)
    n=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(t1,t2))
    n2t1=Mesh2D._gmsh_dot3(t1,t1);n2t2=Mesh2D._gmsh_dot3(t2,t2)
    n2t3=Mesh2D._gmsh_dot3(t3,t3)
    m=((t1[1],t2[1],-n[1]),(t1[2],t2[2],-n[2]),(t1[3],t2[3],-n[3]))
    det=m[1][1]*(m[2][2]*m[3][3]-m[2][3]*m[3][2])-
        m[1][2]*(m[2][1]*m[3][3]-m[2][3]*m[3][1])+
        m[1][3]*(m[2][1]*m[3][2]-m[2][2]*m[3][1])
    det==0.0 && return nothing
    ud=1.0/det
    i00=(m[2][2]*m[3][3]-m[2][3]*m[3][2])*ud
    i10=-(m[2][1]*m[3][3]-m[2][3]*m[3][1])*ud
    i20=(m[2][1]*m[3][2]-m[2][2]*m[3][1])*ud
    i01=-(m[1][2]*m[3][3]-m[1][3]*m[3][2])*ud
    i11=(m[1][1]*m[3][3]-m[1][3]*m[3][1])*ud
    i21=-(m[1][1]*m[3][2]-m[1][2]*m[3][1])*ud
    i02=(m[1][2]*m[2][3]-m[1][3]*m[2][2])*ud
    i12=-(m[1][1]*m[2][3]-m[1][3]*m[2][1])*ud
    i22=(m[1][1]*m[2][2]-m[1][2]*m[2][1])*ud
    pp1=Mesh2D._gmsh_sub3(p,p1)
    u=i00*pp1[1]+i01*pp1[2]+i02*pp1[3]
    v=i10*pp1[1]+i11*pp1[2]+i12*pp1[3]
    d=i20*pp1[1]+i21*pp1[2]+i22*pp1[3]
    sign=d>0 ? 1.0 : -1.0
    d==0.0 && (sign=1.0)
    along(a,e,t)=(a[1]+e[1]*t,a[2]+e[2]*t,a[3]+e[3]*t)
    if u>=0.0 && v>=0.0 && 1.0-u-v>=0.0
        close=(p1[1]+t1[1]*u+t2[1]*v,p1[2]+t1[2]*u+t2[2]*v,p1[3]+t1[3]*u+t2[3]*v)
        return d,close
    end
    t12=Mesh2D._gmsh_dot3(pp1,t1)/n2t1
    t13=Mesh2D._gmsh_dot3(pp1,t2)/n2t2
    pp2=Mesh2D._gmsh_sub3(p,p2)
    t23=Mesh2D._gmsh_dot3(pp2,t3)/n2t3
    d=1.e10;close=p1
    if t12>=0 && t12<=1.0
        q=along(p1,t1,t12)
        d=sign*min(abs(d),_discrete_face_distance(p,q));close=q
    end
    if t13>=0 && t13<=1.0
        q=along(p1,t2,t13)
        _discrete_face_distance(p,q)<abs(d) && (close=q)
        d=sign*min(abs(d),_discrete_face_distance(p,q))
    end
    if t23>=0 && t23<=1.0
        q=along(p2,t3,t23)
        _discrete_face_distance(p,q)<abs(d) && (close=q)
        d=sign*min(abs(d),_discrete_face_distance(p,q))
    end
    for q in (p1,p2,p3)
        if _discrete_face_distance(p,q)<abs(d)
            close=q
            d=sign*min(abs(d),_discrete_face_distance(p,q))
        end
    end
    return d,close
end

# `discreteFace::closestPoint(p, maxDistance)`: (xyz, uv, triangle).
function _discrete_face_closest_point(P::_DiscreteFaceParam,p::NTuple{3,Float64},
                                      max_distance::Float64=1e-6)
    grid=P.grid3
    candidates=Int32[]
    best=0;best_distance=1.e22;best_point=p
    radius=max_distance
    for _ in 1:2000
        lo=(p[1]-radius,p[2]-radius,p[3]-radius)
        hi=(p[1]+radius,p[2]+radius,p[3]+radius)
        empty!(candidates)
        for k in _discrete_face_grid_cell(grid,3,lo[3]):_discrete_face_grid_cell(grid,3,hi[3]),
                j in _discrete_face_grid_cell(grid,2,lo[2]):_discrete_face_grid_cell(grid,2,hi[2]),
                i in _discrete_face_grid_cell(grid,1,lo[1]):_discrete_face_grid_cell(grid,1,hi[1])
            c=_discrete_face_grid_linear(grid,i,j,k)
            for s in grid.start[c]:grid.start[c+1]-1
                push!(candidates,grid.items[s])
            end
        end
        sort!(candidates);unique!(candidates)
        best_distance=1.e22
        for t in candidates
            tri=P.triangles[t]
            a=P.xyz[tri[1]];b=P.xyz[tri[2]];c=P.xyz[tri[3]]
            # RTree overlap test on the triangle's bounding box.
            (min(a[1],b[1],c[1])<=hi[1] && max(a[1],b[1],c[1])>=lo[1] &&
             min(a[2],b[2],c[2])<=hi[2] && max(a[2],b[2],c[2])>=lo[2] &&
             min(a[3],b[3],c[3])<=hi[3] && max(a[3],b[3],c[3])>=lo[3]) || continue
            result=_discrete_face_signed_distance(a,b,c,p)
            result===nothing && continue
            d,close=result
            if abs(d)<best_distance
                best_distance=abs(d);best=Int(t);best_point=close
            end
        end
        best!=0 && break
        radius*=2.0
    end
    best==0 && return p,(0.0,0.0),0
    tri=P.triangles[best]
    x0=P.xyz[tri[1]];x1=P.xyz[tri[2]];x2=P.xyz[tri[3]]
    # `MElement::xyz2uvw` of the closest point on the linear triangle.
    e1=Mesh2D._gmsh_sub3(x1,x0);e2=Mesh2D._gmsh_sub3(x2,x0)
    r=Mesh2D._gmsh_sub3(best_point,x0)
    a11=Mesh2D._gmsh_dot3(e1,e1);a12=Mesh2D._gmsh_dot3(e1,e2);a22=Mesh2D._gmsh_dot3(e2,e2)
    b1=Mesh2D._gmsh_dot3(r,e1);b2=Mesh2D._gmsh_dot3(r,e2)
    det=a11*a22-a12*a12
    V=det==0 ? 0.0 : (b1*a22-a12*b2)/det
    W=det==0 ? 0.0 : (a11*b2-a12*b1)/det
    U=1.0-V-W
    q0=P.uv[tri[1]];q1=P.uv[tri[2]];q2=P.uv[tri[3]]
    uv=(U*q0[1]+V*q1[1]+W*q2[1],U*q0[2]+V*q1[2]+W*q2[2])
    xyz=(U*x0[1]+V*x1[1]+W*x2[1],U*x0[2]+V*x1[2]+W*x2[2],U*x0[3]+V*x1[3]+W*x2[3])
    return xyz,uv,best
end

# `discreteFace::parFromPoint`.
_discrete_face_par_from_point(P::_DiscreteFaceParam,p)=
    _discrete_face_closest_point(P,(Float64(p[1]),Float64(p[2]),Float64(p[3])),1e-6)[2]

# `discreteFace::intersectionWithCircle`: the frontal point on the circle of
# radius `R` about `p` in the plane normal to `n1 x n2`, searched in the
# parameter triangle at `uv` first and then in every triangle in order.
function _discrete_face_intersection_with_circle(P::_DiscreteFaceParam,n1,n2,p,R,uv)
    start=_discrete_face_locate(P,uv[1],uv[2])
    n=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(n1,n2))
    N=length(P.triangles)
    first=start==0 ? 1 : 0
    for i in first:N
        t=i==0 ? start : i
        tri=P.triangles[t]
        v0=P.xyz[tri[1]];v1=P.xyz[tri[2]];v2=P.xyz[tri[3]]
        q0=P.uv[tri[1]];q1=P.uv[tri[2]];q2=P.uv[tri[3]]
        t1=Mesh2D._gmsh_sub3(v1,v0);t2=Mesh2D._gmsh_sub3(v2,v0)
        tn=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(t1,t2))
        d=Mesh2D._gmsh_cross3(n,tn)
        Mesh2D._gmsh_norm3(d)<1.e-12 && continue
        d=Mesh2D._gmsh_normalize3(d)
        rhs=(Mesh2D._gmsh_dot3(n,p),Mesh2D._gmsh_dot3(v0,tn))
        x0=(0.0,0.0,0.0)
        m=((n[2],n[3]),(tn[2],tn[3]))
        if abs(m[1][1]*m[2][2]-m[2][1]*m[1][2])>1.e-12
            r=something(_discrete_face_sys2x2(m,rhs),(0.0,0.0))
            x0=(0.0,r[1],r[2])
        else
            m=((n[1],n[3]),(tn[1],tn[3]))
            if abs(m[1][1]*m[2][2]-m[2][1]*m[1][2])>1.e-12
                r=something(_discrete_face_sys2x2(m,rhs),(0.0,0.0))
                x0=(r[1],0.0,r[2])
            else
                m=((n[1],n[2]),(tn[1],tn[2]))
                r=_discrete_face_sys2x2(m,rhs)
                r===nothing && continue
                x0=(r[1],r[2],0.0)
            end
        end
        px=Mesh2D._gmsh_sub3(p,x0)
        a=1.0
        b=-2*Mesh2D._gmsh_dot3(d,px)
        c=Mesh2D._gmsh_dot3(px,px)-R*R
        delta=b*b-4*a*c
        delta>=0 || continue
        sign=Mesh2D._gmsh_dot3(n2,d)>0 ? 1.0 : -1.0
        ta=(-b+sign*sqrt(delta))/(2.0*a)
        tb=(-b-sign*sqrt(delta))/(2.0*a)
        for s in ((x0[1]+d[1]*ta,x0[2]+d[2]*ta,x0[3]+d[3]*ta),
                  (x0[1]+d[1]*tb,x0[2]+d[2]*tb,x0[3]+d[3]*tb))
            mat=((Mesh2D._gmsh_dot3(t1,t1),Mesh2D._gmsh_dot3(t1,t2)),
                 (Mesh2D._gmsh_dot3(t1,t2),Mesh2D._gmsh_dot3(t2,t2)))
            sv=Mesh2D._gmsh_sub3(s,v0)
            local_uv=_discrete_face_sys2x2(mat,(Mesh2D._gmsh_dot3(sv,t1),
                                                 Mesh2D._gmsh_dot3(sv,t2)))
            local_uv===nothing && (local_uv=(0.0,0.0))
            lu1,lu2=local_uv
            if lu1>=-1.e-6 && lu2>=-1.e-6 && 1.0-lu1-lu2>=-1.e-6
                w=1.0-lu1-lu2
                return (q0[1]*w+q1[1]*lu1+q2[1]*lu2,q0[2]*w+q1[2]*lu1+q2[2]*lu2)
            end
        end
    end
    return nothing
end

# `sys2x2`; `nothing` for the singular branch (upstream zeroes `res`).
function _discrete_face_sys2x2(m,b)
    norm=m[1][1]^2+m[2][2]^2+m[1][2]^2+m[2][1]^2
    det=m[1][1]*m[2][2]-m[2][1]*m[1][2]
    (norm==0.0 || abs(det)/norm<1.e-16) && return nothing
    ud=1.0/det
    return ((b[1]*m[2][2]-m[1][2]*b[2])*ud,(m[1][1]*b[2]-m[2][1]*b[1])*ud)
end

# Evaluation frames for `ModelEntityEvaluation`'s discrete-surface queries:
# the stored parameter triangles with their spatial corners.
function _discrete_face_param_frames(P::_DiscreteFaceParam)
    frames=Tuple{Vector{Vector{Float64}},Vector{NTuple{3,Float64}}}[]
    for tri in P.triangles
        push!(frames,([Float64[P.uv[k][1],P.uv[k][2]] for k in tri],
                      NTuple{3,Float64}[P.xyz[k] for k in tri]))
    end
    return frames
end

# ---- the chart the native surface producer meshes on -------------------------

struct _DiscreteFaceLift
    param::_DiscreteFaceParam
end
(lift::_DiscreteFaceLift)(u,v)=_discrete_face_point(lift.param,Float64(u),Float64(v))[1]
Mesh2D._gmsh_surface_point_succeeded(lift::_DiscreteFaceLift,u,v)=
    _discrete_face_locate(lift.param,Float64(u),Float64(v))!=0

# `buildMetric`: `firstDer` dot products.
struct _DiscreteFaceForm
    param::_DiscreteFaceParam
end
function (form::_DiscreteFaceForm)(u,v)
    du,dv=_discrete_face_first_derivatives(form.param,Float64(u),Float64(v))
    return (Mesh2D._gmsh_dot3(du,du),Mesh2D._gmsh_dot3(dv,du),Mesh2D._gmsh_dot3(dv,dv))
end

const _DiscreteFaceChart=_GmshRuledChart{_DiscreteFaceLift,_DiscreteFaceForm}

_discrete_face_chart(m::GeoModel,tag::Int,P::_DiscreteFaceParam)=
    _GmshRuledChart(_DiscreteFaceLift(P),_DiscreteFaceForm(P),tag,_model_lc(m),false)

# `reparamMeshVertexOnFace` on a DiscreteSurface: `parFromPoint` of the
# node's coordinates, whatever entity holds it.
function _gmsh_ruled_source_parameters(m::GeoModel,chart::_DiscreteFaceChart,
                                       owner,parameter,point,caller)
    return _discrete_face_par_from_point(chart.lift.param,point)
end

# `optimalPointFrontalB` on a DiscreteSurface (HXT build):
# `intersectionWithCircle` replaces `intersectCurveSurface`. Upstream's
# solution lands in a local `uv[2]` that shadows the caller's `newPoint`
# (discreteFace.cpp), so a successful intersection inserts the fast point
# itself and a failed one inserts nothing. A degenerate seed cell keeps the
# fast point.
function Mesh2D._gmsh_frontal_surface_point(fast,lift::_DiscreteFaceLift,v1,v2,v3,d)
    middle=((v1[1]+v2[1])*0.5,(v1[2]+v2[2])*0.5,(v1[3]+v2[3])*0.5)
    edge=Mesh2D._gmsh_sub3(v2,v1)
    toward=Mesh2D._gmsh_sub3(v3,middle)
    n1=Mesh2D._gmsh_cross3(edge,toward)
    Mesh2D._gmsh_norm3(n1)<1e-12 && return fast
    n2=Mesh2D._gmsh_cross3(n1,edge)
    n1=Mesh2D._gmsh_normalize3(n1);n2=Mesh2D._gmsh_normalize3(n2)
    found=_discrete_face_intersection_with_circle(lift.param,n1,n2,middle,d,fast)
    return found===nothing ? nothing : (fast[1],fast[2])
end

# ---- discreteFace::mesh --------------------------------------------------------

# Boundary Curves in `gf->edges()` order: the classification's bound edges.
function _discrete_face_curves(m::GeoModel,tag::Int)
    record=m.discrete[(2,tag)]
    curves=Int[];seen=Set{Int}()
    for (dimension,curve) in record.boundary
        dimension==1 || continue
        curve=abs(curve)
        curve in seen && continue
        push!(curves,curve);push!(seen,curve)
    end
    return curves
end

function _discrete_face_source(m::GeoModel,tag::Int,input,chart,caller)
    source=_ModelSurfaceSource()
    for curve in _discrete_face_curves(m,tag)
        _gmsh_native_surface_curve!(source,m,input,curve,chart,false,caller;
            surface_tag=tag,ctx_lc=0.0)
    end
    return source
end

# `orientMeshGFace` for a parametrized discrete face: the first cell with a
# vertex the face owns compares its normal with `discreteFace::normal` there.
function _discrete_face_orientation(P::_DiscreteFaceParam,mesh,owners,parameters,tag)
    owner=(2,Int32(tag))
    for column in axes(mesh.tris,2)
        for slot in 1:3
            row=mesh.tris[slot,column]
            owners[row]==owner || continue
            uv=parameters[row]
            (uv===nothing || length(uv)<2) && continue
            nf=_discrete_face_normal(P,uv[1],uv[2])
            a=ntuple(k->mesh.coords[k,mesh.tris[1,column]],3)
            b=ntuple(k->mesh.coords[k,mesh.tris[2,column]],3)
            c=ntuple(k->mesh.coords[k,mesh.tris[3,column]],3)
            ne=Mesh2D._gmsh_cross3(Mesh2D._gmsh_sub3(b,a),Mesh2D._gmsh_sub3(c,a))
            return Mesh2D._gmsh_dot3(ne,nf)>0 ? 1 : -1
        end
    end
    return 0
end

# `discreteFace::mesh` → `GFace::mesh` → `meshGenerator` on the stored
# parametrization: the ruled-chart producer with the discrete lift, metric,
# boundary reparametrization and frontal point.
function _discrete_face_surface_production(m::GeoModel,tag::Int,options,input,
        high_water,caller;algorithm=6,size_field=nothing,element_high_water=UInt64(0))
    P=_discrete_face_param(m,tag)
    P===nothing && throw(ArgumentError(
        "$caller: discrete Surface[$tag] has no parametrization"))
    algorithm in (3,5,6) || throw(ArgumentError(
        "$caller: mesh algorithm $algorithm on parametrized discrete Surface[$tag] " *
        "is not implemented; supported: 3 (initial mesh only), 5 (Delaunay), " *
        "6 (Frontal-Delaunay)"))
    authority,limit=_gmsh_native_surface_lower_input(m,options,input,high_water,caller)
    chart=_discrete_face_chart(m,tag,P)
    source=_discrete_face_source(m,tag,authority,chart,caller)
    isempty(source.tags) && throw(ArgumentError(
        "$caller: discrete Surface[$tag] has no boundary mesh"))
    ground=chart.form
    T,inside,initial=_gmsh_native_ruled_initial(source,tag,ground)
    curve_chains=_gmsh_native_surface_curve_chains(source)
    numbers=source.tags
    trace=nothing
    if length(source.tags)!=3 && algorithm!=3
        sizes=_gmsh_native_surface_vertex_sizes(T,inside,source,m,Int[])
        background=_GmshRuledSurfaceBackground(chart,
            _ModelAlgorithmBackgroundSize(m,2,tag,options,size_field,caller))
        background_sizes=copy(sizes)
        extend=_resolved_extend_from_boundary(m,2,tag)
        trace=Mesh2D._gmsh_surface_trace(T,inside,source.world,chart.lift,
            element_high_water,sizes,background_sizes,extend>0 || extend==-2;
            temporary_cells=count(inside),
            metric=Mesh2D._GmshFirstFundamental(ground),initial)
        insertion=refine!(T;min_angle_deg=0,gmsh_insertion=true,
            gmsh_algorithm=algorithm==5 ? :delaunay : :frontal,
            gmsh_extend_boundary=extend>0 || extend==-2,
            vertex_sizes=sizes,vertex_bgm=background_sizes,bgm_size=background,
            model_vertices=Int32[],
            gmsh_vertex_numbers=source.tags,gmsh_vertex_high_water=limit,
            return_vertex_numbers=true,_gmsh_surface_trace=trace)
        inside=insertion.interior;numbers=insertion.vertex_numbers
        limit=insertion.vertex_high_water
    end
    product=_gmsh_native_surface_export(T,inside,numbers,source,chart,tag,limit;
        projected=true,trace)
    if trace===nothing
        cells=ntris(product.mesh)
        births=UInt64[element_high_water+UInt64(column) for column in 1:cells]
        product=_ModelSurfaceProduction(product.mesh,product.node_tags,
            product.node_owners,product.node_parameters,product.retained_nodes,
            product.node_high_water,births,element_high_water+UInt64(cells),
            product.birth_rows)
    end
    mesh=product.mesh
    physical=Matrix{Float64}(undef,3,nnodes(mesh))
    for row in axes(physical,2)
        original=get(source.positions,product.node_tags[row],0)
        physical[:,row].=original>0 ? source.world[original] :
            _surface_metric_lift(chart,mesh.coords[1,row],mesh.coords[2,row])
    end
    if trace!==nothing
        product,physical=_gmsh_native_surface_smooth(m,tag,product,physical,chart,caller)
        mesh=product.mesh
    end
    nnodes(mesh)==length(product.node_tags) || throw(ArgumentError(
        "$caller: Surface[$tag] postprocessing changed the actual producer node set"))
    for row in axes(mesh.coords,2)
        mesh.coords[1,row]=physical[1,row];mesh.coords[2,row]=physical[2,row]
        mesh.coords[3,row]=physical[3,row]
    end
    # Not ported: meshGFace::operator() re-meshes a face whose Delaunay mesh
    # fails `isMeshValid` with MeshAdapt (Mesh.AlgorithmSwitchOnFailure).
    _discrete_face_orientation(P,mesh,product.node_owners,product.node_parameters,tag)<0 &&
        (mesh=_reversed_surface_mesh(mesh,caller,tag))
    get(m.meshing.reverse,(2,tag),false) &&
        (mesh=_reversed_surface_mesh(mesh,caller,tag))
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid discrete Surface[$tag] — "*join(diagnostic.messages,"; ")))
    return _ModelSurfaceProduction(mesh,product.node_tags,product.node_owners,
        product.node_parameters,product.retained_nodes,product.node_high_water,
        product.element_tags,product.element_high_water,product.birth_rows,nothing,
        product.source_node_updates,product.native_recombine_events,curve_chains)
end

function _discrete_face_mesh_model_surface(m::GeoModel,tag::Int,caller;
        size_field=nothing,_mesh_input=nothing,_mesh_node_high_water=nothing,
        _mesh_element_high_water=UInt64(0),_mesh1d_options=nothing,
        _return_surface_production=false)
    options=_model_algorithm_size_options(m,_mesh1d_options,caller)
    algorithm=_surface_mesh_algorithm(m,tag)
    production=_discrete_face_surface_production(m,tag,options,_mesh_input,
        _mesh_node_high_water,caller;algorithm,size_field,
        element_high_water=_mesh_element_high_water)
    return _return_surface_production ? production : production.mesh
end
