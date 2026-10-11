struct _ModelNURBSProjectionSample
    uv::Vector{Float64}
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
function _model_nurbs_projection_sample(cache,target,uv,kind,i,j)
    xyz=_model_nurbs_cache_evaluate(cache,uv[1],uv[2],0)
    return _ModelNURBSProjectionSample(uv,xyz,_model_nurbs_root_distance(xyz,target),kind,i,j)
end
function _model_nurbs_projection_edge(cache,target,a,b,axis)
    chord=_model_nurbs_root_distance(a.xyz,b.xyz)
    difference=abs(a.squared_distance-b.squared_distance)
    chord<=2e-9 && return a
    difference>=chord-2e-9 && return a.squared_distance>b.squared_distance ? b : a
    ratio=_model_nurbs_root_dot(target.-a.xyz,b.xyz.-a.xyz)/chord
    uv=copy(a.uv);uv[axis]+=ratio*(b.uv[axis]-a.uv[axis])
    return _model_nurbs_projection_sample(cache,target,uv,axis==1 ? :u : :v,a.i,a.j)
end
_model_nurbs_projection_same_sample(a,b)=a!==nothing &&
    a.kind==b.kind && a.i==b.i && a.j==b.j

function _model_nurbs_projection_seeds(cache,target,lower,upper)
    surface=cache.patch.surface
    us=_model_nurbs_projection_knot_axis(surface.knots_u,surface.degree_u,lower[1],upper[1],44)
    vs=_model_nurbs_projection_knot_axis(surface.knots_v,surface.degree_v,lower[2],upper[2],44)
    nu,nv=length(us),length(vs)
    nodes=Matrix{_ModelNURBSProjectionSample}(undef,nu,nv)
    for i in 1:nu,j in 1:nv
        nodes[i,j]=_model_nurbs_projection_sample(cache,target,[us[i],vs[j]],:node,i,j)
    end
    ue=Matrix{_ModelNURBSProjectionSample}(undef,nu-1,nv)
    ve=Matrix{_ModelNURBSProjectionSample}(undef,nu,nv-1)
    for i in 1:nu,j in 1:nv
        i<nu && (ue[i,j]=_model_nurbs_projection_edge(cache,target,nodes[i,j],nodes[i+1,j],1))
        j<nv && (ve[i,j]=_model_nurbs_projection_edge(cache,target,nodes[i,j],nodes[i,j+1],2))
    end
    faces=Matrix{_ModelNURBSProjectionSample}(undef,nu-1,nv-1)
    for i in 1:nu-1,j in 1:nv-1
        u0,u1=ue[i,j],ue[i,j+1];v0,v1=ve[i,j],ve[i+1,j]
        out=abs(u0.squared_distance-u1.squared_distance)>=_model_nurbs_root_distance(u0.xyz,u1.xyz)-2e-9 ||
            abs(v0.squared_distance-v1.squared_distance)>=_model_nurbs_root_distance(v0.xyz,v1.xyz)-2e-9
        if out
            umin=u0.squared_distance<u1.squared_distance ? u0 : u1
            vmin=v0.squared_distance<v1.squared_distance ? v0 : v1
            faces[i,j]=umin.squared_distance<vmin.squared_distance ? umin : vmin
        else
            uv=[0.5*(u0.uv[1]+u1.uv[1]),0.5*(v0.uv[2]+v1.uv[2])]
            faces[i,j]=_model_nurbs_projection_sample(cache,target,uv,:face,i,j)
        end
    end
    starts=Vector{Float64}[]
    for i in 1:nu-1,j in 1:nv-1
        sample=faces[i,j];kind,iu,iv=sample.kind,sample.i,sample.j
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
                    admit=j>1 && _model_nurbs_projection_same_sample(faces[i,j-1],sample)
                elseif kind==:v || kind==:node && (iv==1 || iv==nv)
                    admit=i>1 && _model_nurbs_projection_same_sample(faces[i-1,j],sample)
                elseif iu==i && iv==j
                    admit=i>1 && j>1 &&
                        _model_nurbs_projection_same_sample(faces[i,j-1],sample) &&
                        _model_nurbs_projection_same_sample(faces[i-1,j-1],sample) &&
                        _model_nurbs_projection_same_sample(faces[i-1,j],sample)
                end
            end
        end
        admit && push!(starts,copy(sample.uv))
    end
    neighbors=((-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1))
    for i in 1:nu,j in 1:nv
        admit=true
        for (di,dj) in neighbors
            ii,jj=i+di,j+dj
            if 1<=ii<=nu && 1<=jj<=nv && nodes[ii,jj].squared_distance>nodes[i,j].squared_distance
                admit=false;break
            end
        end
        admit && push!(starts,copy(nodes[i,j].uv))
    end
    return starts
end

function _model_nurbs_projection_cached(patch,target)
    surface=patch.surface
    lower=[surface.knots_u[surface.degree_u+1],surface.knots_v[surface.degree_v+1]]
    upper=[surface.knots_u[end-surface.degree_u],surface.knots_v[end-surface.degree_v]]
    for i in 1:2
        allowance=max(abs(upper[i]-lower[i])*1e-8,1e-12)
        lower[i]-=allowance;upper[i]+=allowance
    end
    cache=_ModelNURBSProjectionCache(patch)
    starts=_model_nurbs_projection_seeds(cache,target,lower,upper)
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
