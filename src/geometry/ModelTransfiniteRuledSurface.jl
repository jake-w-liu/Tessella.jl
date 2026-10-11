# A ruled transfinite Surface interpolates boundary parameters before its CAD
# evaluation. Its logical grid order is separate from public node birth order.
struct _ModelTransfiniteRuledGrid
    model::GeoModel
    surface::Int
    sides::NTuple{4,Vector{NTuple{2,Float64}}}
    parameters::Dict{Tuple{Int,Int},NTuple{2,Float64}}
    u::Vector{Float64}
    v::Vector{Float64}
    caller::String
end

function _model_transfinite_ruled_chords(first,opposite,caller)
    length(first)==length(opposite) || throw(ArgumentError("$caller: ruled TF opposite side counts differ"))
    result=zeros(Float64,length(first));total=0.0
    for i in 1:length(first)-1
        a=ntuple(axis->first[i+1][axis]-first[i][axis],3)
        b=ntuple(axis->opposite[i+1][axis]-opposite[i][axis],3)
        d1=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
        d2=sqrt(b[1]*b[1]+b[2]*b[2]+b[3]*b[3])
        step=0.5*(d1+d2)
        next=total+step
        isfinite(next) && next>total || throw(ArgumentError("$caller: ruled TF chord is not representable"))
        total=next;result[i+1]=total
    end
    total>0 || throw(ArgumentError("$caller: ruled TF side is empty"))
    result./=total
    all(i->isfinite(result[i]) && result[i-1]<result[i],2:length(result)) ||
        throw(ArgumentError("$caller: ruled TF chord parameters do not increase"))
    return result
end

function _model_transfinite_ruled_grid(m,t,curves,frame,points,input,caller)
    _surface_type(m,t)===:ruled && !_occ_surface_entity(m,t) || return nothing
    length(curves)==4 && all(curve->_curve_type(m,abs(curve))===:line,curves) || return nothing
    gens=m.loops[only(m.surfaces[t])]
    length(gens)==4 && allunique(abs.(gens)) || return nothing
    keys=_model_transfinite_side_keys(m,curves,frame,caller)
    sides=ntuple(4) do side
        signed=keys[side];curve=abs(signed)
        natural=findfirst(gen->abs(gen)==curve,gens)
        natural===nothing && throw(ArgumentError(
            "$caller: ruled Surface[$t] grid loses its CAD generatrix"))
        source=_model_mesh_curve_births(input,m,curve,caller)
        params=if source===nothing
            t0,t1=_model_curve_param_bounds(m,curve,caller)
            _model_curve_transfinite_native_params(m,curve,t0,t1,
                m.meshing.transfinite_curves[curve],caller)
        else
            data=source.mesh.entity_data
            map(enumerate(source.rows)) do entry
                index,row=entry
                index==1 && return 0.0
                index==length(source.rows) && return 1.0
                raw=data.node_parametric[row]
                return raw!==nothing && length(raw)==1 ? only(raw) :
                    only(model_parametrization(m,1,curve,source.mesh.coords[:,row]))
            end
        end
        all(isfinite,params) || throw(ArgumentError(
            "$caller: ruled Surface[$t] boundary has non-finite Curve parameters"))
        uv=NTuple{2,Float64}[]
        for raw in params
            p=gens[natural]>0 ? raw : 1.0-raw
            push!(uv,natural==1 ? (p,0.0) : natural==2 ? (1.0,p) :
                natural==3 ? (1.0-p,1.0) : (0.0,1.0-p))
        end
        signed<0 && reverse!(uv)
        uv
    end
    u=_model_transfinite_ruled_chords(points[1],reverse(points[3]),caller)
    v=_model_transfinite_ruled_chords(points[2],reverse(points[4]),caller)
    return _ModelTransfiniteRuledGrid(m,t,sides,
        Dict{Tuple{Int,Int},NTuple{2,Float64}}(),u,v,String(caller))
end

function Transfinite._interpolate_point(grid::_ModelTransfiniteRuledGrid,u,v,i::Int,j::Int)
    bottom,right,back_top,back_left=grid.sides
    L=length(bottom)-1;H=length(right)-1
    length(back_top)==L+1 && length(back_left)==H+1 && 0<i<L && 0<j<H ||
        throw(ArgumentError("$(grid.caller): ruled Surface grid interpolation loses its side rows"))
    u=grid.u[i+1];v=grid.v[j+1]
    left=back_left[H-j+1];top=back_top[L-i+1]
    s1=bottom[1];s2=bottom[end];s3=right[end];s4=back_top[end]
    uv=ntuple(2) do axis
        (1.0-u)*left[axis]+u*right[j+1][axis]+(1.0-v)*bottom[i+1][axis]+v*top[axis]-
        ((1.0-u)*(1.0-v)*s1[axis]+u*(1.0-v)*s2[axis]+
         u*v*s3[axis]+(1.0-u)*v*s4[axis])
    end
    grid.parameters[(i,j)]=uv
    return _ruled_surface_point(grid.model,grid.surface,uv[1],uv[2],grid.caller)
end

_model_transfinite_ruled_parameters(grid,mesh,L,H,caller)=nothing
function _model_transfinite_ruled_parameters(grid::_ModelTransfiniteRuledGrid,mesh,L,H,caller)
    length(grid.parameters)==(L-1)*(H-1) || throw(ArgumentError(
        "$caller: ruled Surface[$(grid.surface)] lacks its emitted UV rows"))
    parameters=Dict{Int32,Vector{Float64}}()
    for ((i,j),uv) in grid.parameters
        0<i<L && 0<j<H || throw(ArgumentError("$caller: ruled Surface UV row is outside its interior"))
        row=Int32(i+1+j*(L+1))
        point=_ruled_surface_point(grid.model,grid.surface,uv[1],uv[2],caller)
        all(axis->isequal(point[axis],mesh.coords[axis,row]),1:3) || throw(ArgumentError(
            "$caller: ruled Surface[$(grid.surface)] emitted UV no longer binds its actual coordinate"))
        parameters[row]=collect(uv)
    end
    return parameters
end

# A triangular ruled Surface retains the UV returned by its projection. Public
# metadata uses that operation's parameters, rather than inverting its XYZ again.
struct _ModelTransfiniteRuledTriProjection
    model::GeoModel
    surface::Int
    lc::Float64
    caller::String
    parameters::Dict{Int32,NTuple{2,Float64}}
end

function _model_transfinite_ruled_tri_uv(project::_ModelTransfiniteRuledTriProjection,point)
    silent=Returns(nothing)
    u,v=_ruled_xyz_to_uv(project.model,project.surface,point,1.0,project.lc,
        false,silent,silent,project.caller;on_surface=false)
    return (u,v)
end

function (project::_ModelTransfiniteRuledTriProjection)(point)
    u,v=_model_transfinite_ruled_tri_uv(project,point)
    return _ruled_surface_point(project.model,project.surface,u,v,project.caller)
end

function TransfiniteTriangle._project_point(project::_ModelTransfiniteRuledTriProjection,
                                          point,row::Int32)
    row>0 || throw(ArgumentError("$(project.caller): ruled triangular projection lacks its grid row"))
    uv=_model_transfinite_ruled_tri_uv(project,point)
    haskey(project.parameters,row) && throw(ArgumentError(
        "$(project.caller): ruled triangular projection repeats its grid row"))
    project.parameters[row]=uv
    return _ruled_surface_point(project.model,project.surface,uv[1],uv[2],project.caller)
end

_model_transfinite_projection_parameters(project,mesh,caller)=nothing
function _model_transfinite_projection_parameters(project::_ModelTransfiniteRuledTriProjection,
                                                  mesh,caller)
    parameters=Dict{Int32,Vector{Float64}}()
    for (row,uv) in project.parameters
        1<=row<=nnodes(mesh) || throw(ArgumentError("$caller: ruled triangular UV row is outside its grid"))
        point=_ruled_surface_point(project.model,project.surface,uv[1],uv[2],caller)
        all(axis->isequal(point[axis],mesh.coords[axis,row]),1:3) || throw(ArgumentError(
            "$caller: ruled Surface[$(project.surface)] projected UV no longer binds its actual coordinate"))
        parameters[row]=collect(uv)
    end
    return parameters
end
