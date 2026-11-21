# One candidate of the seed search. Its parameters are a pair of numbers, not
# a list: a projection makes several thousand of them.
struct _ModelNURBSProjectionSample
    uv::NTuple{2,Float64}
    xyz::NTuple{3,Float64}
    squared_distance::Float64
    kind::Symbol
    i::Int
    j::Int
end
function _model_nurbs_projection_grid_axis(lower,upper,count)
    width=upper-lower;offset=width/count/100.0
    step=(width-offset)/(count-1);value=offset/2.0+lower
    axis=Vector{Float64}(undef,count)
    for i in 1:count
        axis[i]=value;value+=step
    end
    return axis
end
function _model_nurbs_projection_knot_axis(knots,degree,lower,upper,count)
    unique_knots=unique(knots);parameters=Float64[lower];previous=lower
    for i in 1:length(unique_knots)-1
        unique_knots[i]>=upper-1e-9 && break
        unique_knots[i+1]<lower+1e-9 && continue
        step=(unique_knots[i+1]-unique_knots[i])/max(degree,2)
        for k in 1:degree
            parameter=unique_knots[i]+k*step
            parameter>upper-1e-9 && break
            if parameter>previous+1e-9
                push!(parameters,parameter);previous=parameter
            end
        end
    end
    push!(parameters,upper)
    return length(parameters)<count ?
        _model_nurbs_projection_grid_axis(lower,upper,count) : parameters
end
function _model_nurbs_projection_sample(cache,target,uv::NTuple{2,Float64},kind,i,j)
    xyz=_model_nurbs_cache_evaluate(cache,uv[1],uv[2],0)
    return _ModelNURBSProjectionSample(uv,xyz,_model_nurbs_root_distance(xyz,target),kind,i,j)
end
function _model_nurbs_projection_edge(cache,target,a,b,axis)
    chord=_model_nurbs_root_distance(a.xyz,b.xyz)
    difference=abs(a.squared_distance-b.squared_distance)
    chord<=2e-9 && return a
    difference>=chord-2e-9 && return a.squared_distance>b.squared_distance ? b : a
    ratio=_model_nurbs_root_dot(target.-a.xyz,b.xyz.-a.xyz)/chord
    uv=axis==1 ? (a.uv[1]+ratio*(b.uv[1]-a.uv[1]),a.uv[2]) :
                 (a.uv[1],a.uv[2]+ratio*(b.uv[2]-a.uv[2]))
    return _model_nurbs_projection_sample(cache,target,uv,axis==1 ? :u : :v,a.i,a.j)
end
_model_nurbs_projection_same_sample(a,b)=
    a.kind==b.kind && a.i==b.i && a.j==b.j

# The sample grid of one patch: its parameters and the surface's points there.
# `Extrema_GenExtPS` builds this grid when it is handed the surface
# (`Initialize` -> `BuildGrid`) and only measures it for each point
# (`Perform`). Building it for every projection was 1,936 or more surface
# evaluations per point: opening Gmsh's tutorials/t20.geo, which slices an
# imported STEP solid, projected 32,000 points and spent 90% of a minute and
# 105 GB of allocation on it.
struct _ModelNURBSProjectionGrid
    # What the grid was built from besides the control net that keys it.
    degrees::NTuple{2,Int}
    knots_u::Vector{Float64}
    knots_v::Vector{Float64}
    weights::Matrix{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    us::Vector{Float64}
    vs::Vector{Float64}
    points::Matrix{NTuple{3,Float64}}
end

# Keyed weakly by the control net, which a grid does not reference: a grid
# lives as long as its surface does. Arrays compare by content, so a hit is
# checked against the rest of the surface before it is used.
const _MODEL_NURBS_PROJECTION_GRIDS=
    WeakKeyDict{Matrix{NTuple{3,Float64}},_ModelNURBSProjectionGrid}()

function _model_nurbs_projection_grid_matches(grid::_ModelNURBSProjectionGrid,surface)
    return grid.degrees==(surface.degree_u,surface.degree_v) &&
        (grid.knots_u===surface.knots_u || grid.knots_u==surface.knots_u) &&
        (grid.knots_v===surface.knots_v || grid.knots_v==surface.knots_v) &&
        (grid.weights===surface.weights || grid.weights==surface.weights)
end

function _model_nurbs_projection_build_grid(cache)
    surface=cache.patch.surface
    lower=[surface.knots_u[surface.degree_u+1],surface.knots_v[surface.degree_v+1]]
    upper=[surface.knots_u[end-surface.degree_u],surface.knots_v[end-surface.degree_v]]
    for i in 1:2
        allowance=max(abs(upper[i]-lower[i])*1e-8,1e-12)
        lower[i]-=allowance;upper[i]+=allowance
    end
    us=_model_nurbs_projection_knot_axis(surface.knots_u,surface.degree_u,lower[1],upper[1],44)
    vs=_model_nurbs_projection_knot_axis(surface.knots_v,surface.degree_v,lower[2],upper[2],44)
    points=Matrix{NTuple{3,Float64}}(undef,length(us),length(vs))
    for i in eachindex(us),j in eachindex(vs)
        points[i,j]=_model_nurbs_cache_evaluate(cache,us[i],vs[j],0)
    end
    return _ModelNURBSProjectionGrid((surface.degree_u,surface.degree_v),
        surface.knots_u,surface.knots_v,surface.weights,lower,upper,us,vs,points)
end

function _model_nurbs_projection_grid(cache)
    surface=cache.patch.surface
    grid=get(_MODEL_NURBS_PROJECTION_GRIDS,surface.controls,nothing)
    grid!==nothing && _model_nurbs_projection_grid_matches(grid,surface) && return grid
    grid=_model_nurbs_projection_build_grid(cache)
    _MODEL_NURBS_PROJECTION_GRIDS[surface.controls]=grid
    return grid
end

# The candidates of one projection: grid nodes, grid edges and grid cells.
# One set per task, grown to the largest grid met.
struct _ModelNURBSProjectionScratch
    nodes::Vector{_ModelNURBSProjectionSample}
    ue::Vector{_ModelNURBSProjectionSample}
    ve::Vector{_ModelNURBSProjectionSample}
    faces::Vector{_ModelNURBSProjectionSample}
end
const _MODEL_NURBS_PROJECTION_SCRATCH_KEY=:_tessella_nurbs_projection_scratch

function _model_nurbs_projection_scratch(nu::Int,nv::Int)
    scratch=get!(task_local_storage(),_MODEL_NURBS_PROJECTION_SCRATCH_KEY) do
        _ModelNURBSProjectionScratch(_ModelNURBSProjectionSample[],
            _ModelNURBSProjectionSample[],_ModelNURBSProjectionSample[],
            _ModelNURBSProjectionSample[])
    end::_ModelNURBSProjectionScratch
    length(scratch.nodes)<nu*nv && resize!(scratch.nodes,nu*nv)
    length(scratch.ue)<(nu-1)*nv && resize!(scratch.ue,(nu-1)*nv)
    length(scratch.ve)<nu*(nv-1) && resize!(scratch.ve,nu*(nv-1))
    length(scratch.faces)<(nu-1)*(nv-1) && resize!(scratch.faces,(nu-1)*(nv-1))
    return scratch
end

function _model_nurbs_projection_seeds(cache,grid::_ModelNURBSProjectionGrid,target)
    us,vs,points=grid.us,grid.vs,grid.points
    nu,nv=length(us),length(vs)
    scratch=_model_nurbs_projection_scratch(nu,nv)
    nodes,ue,ve,faces=scratch.nodes,scratch.ue,scratch.ve,scratch.faces
    # Column-major positions in the four lists.
    node(i,j)=(j-1)*nu+i
    uedge(i,j)=(j-1)*(nu-1)+i
    vedge(i,j)=(j-1)*nu+i
    face(i,j)=(j-1)*(nu-1)+i
    @inbounds for i in 1:nu,j in 1:nv
        xyz=points[i,j]
        nodes[node(i,j)]=_ModelNURBSProjectionSample((us[i],vs[j]),xyz,
            _model_nurbs_root_distance(xyz,target),:node,i,j)
    end
    @inbounds for i in 1:nu,j in 1:nv
        i<nu && (ue[uedge(i,j)]=_model_nurbs_projection_edge(
            cache,target,nodes[node(i,j)],nodes[node(i+1,j)],1))
        j<nv && (ve[vedge(i,j)]=_model_nurbs_projection_edge(
            cache,target,nodes[node(i,j)],nodes[node(i,j+1)],2))
    end
    @inbounds for i in 1:nu-1,j in 1:nv-1
        u0,u1=ue[uedge(i,j)],ue[uedge(i,j+1)];v0,v1=ve[vedge(i,j)],ve[vedge(i+1,j)]
        out=abs(u0.squared_distance-u1.squared_distance)>=_model_nurbs_root_distance(u0.xyz,u1.xyz)-2e-9 ||
            abs(v0.squared_distance-v1.squared_distance)>=_model_nurbs_root_distance(v0.xyz,v1.xyz)-2e-9
        if out
            umin=u0.squared_distance<u1.squared_distance ? u0 : u1
            vmin=v0.squared_distance<v1.squared_distance ? v0 : v1
            faces[face(i,j)]=umin.squared_distance<vmin.squared_distance ? umin : vmin
        else
            uv=(0.5*(u0.uv[1]+u1.uv[1]),0.5*(v0.uv[2]+v1.uv[2]))
            faces[face(i,j)]=_model_nurbs_projection_sample(cache,target,uv,:face,i,j)
        end
    end
    starts=Vector{Float64}[]
    @inbounds for i in 1:nu-1,j in 1:nv-1
        sample=faces[face(i,j)];kind,iu,iv=sample.kind,sample.i,sample.j
        admit=kind==:face
        if !admit
            if kind==:u
                admit=iv==1 || iv==nv
            elseif kind==:v
                admit=iu==1 || iu==nu
            elseif kind==:node
                admit=(iu==1 || iu==nu) && (iv==1 || iv==nv)
            end
            if !admit
                if kind==:u || kind==:node && (iu==1 || iu==nu)
                    admit=j>1 && _model_nurbs_projection_same_sample(faces[face(i,j-1)],sample)
                elseif kind==:v || kind==:node && (iv==1 || iv==nv)
                    admit=i>1 && _model_nurbs_projection_same_sample(faces[face(i-1,j)],sample)
                elseif iu==i && iv==j
                    admit=i>1 && j>1 &&
                        _model_nurbs_projection_same_sample(faces[face(i,j-1)],sample) &&
                        _model_nurbs_projection_same_sample(faces[face(i-1,j-1)],sample) &&
                        _model_nurbs_projection_same_sample(faces[face(i-1,j)],sample)
                end
            end
        end
        admit && push!(starts,[sample.uv[1],sample.uv[2]])
    end
    neighbors=((-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1))
    @inbounds for i in 1:nu,j in 1:nv
        here=nodes[node(i,j)].squared_distance
        admit=true
        for (di,dj) in neighbors
            ii,jj=i+di,j+dj
            if 1<=ii<=nu && 1<=jj<=nv && nodes[node(ii,jj)].squared_distance>here
                admit=false;break
            end
        end
        admit && push!(starts,[us[i],vs[j]])
    end
    return starts
end

function _model_nurbs_projection_cached(patch,target)
    cache=_ModelNURBSProjectionCache(patch)
    grid=_model_nurbs_projection_grid(cache)
    lower,upper=grid.lower,grid.upper
    starts=_model_nurbs_projection_seeds(cache,grid,target)
    function_value=_ModelNURBSProjectionFunction(cache,target)
    states=Tuple{Vector{Float64},NTuple{3,Float64}}[]
    for start in starts
        root=_model_nurbs_projection_root!(function_value,start,lower,upper)
        root===nothing && continue
        recorded_uv=copy(root.last_uv)
        any(state->_model_nurbs_root_distance(state[1],recorded_uv)<=1e-18,states) && continue
        push!(states,(recorded_uv,root.last_xyz))
    end
    isempty(states) && throw(ArgumentError("NURBS projection: no normal-equation state was recorded"))
    best=states[1]
    for state in states[2:end]
        if _model_nurbs_root_distance(state[2],target)<_model_nurbs_root_distance(best[2],target)
            best=state
        end
    end
    return (uv=best[1],xyz=best[2],starts=starts,states=states)
end
