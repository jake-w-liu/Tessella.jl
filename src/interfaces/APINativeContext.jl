# CTX::lc is owned by the API session, independently of any model or mesh.
const NATIVE_CONTEXT_LC=Ref{Float64}(1.0)
const NATIVE_CONTEXT_FORCED_BBOX=Ref{Bool}(false)

function _reset_native_context_locked!()
    NATIVE_CONTEXT_LC[]=1.0
    NATIVE_CONTEXT_FORCED_BBOX[]=false
    return nothing
end

function _synchronize_native_context_locked!(m::GeoModel)
    NATIVE_CONTEXT_FORCED_BBOX[] && return nothing
    NATIVE_CONTEXT_LC[]=Model._model_mesh_lc(
        m,OPTIONS["Geometry.Tolerance"],"API.model.synchronize")
    return nothing
end

function _synchronize_native_context()
    return lock(STATE_LOCK) do
        model=_model_locked()
        Model._native_geo_kernel_synchronize!(model;
            context_lc=NATIVE_CONTEXT_LC[],geometry_tolerance=OPTIONS["Geometry.Tolerance"])
        _synchronize_native_context_locked!(model)
        _replace_mesh_cache_locked!(nothing;preserve_source_inputs=true)
    end
end

function _set_native_geo_transfinite_curve(tag,num_nodes,mesh_type="Progression",coef=1.0)
    return lock(STATE_LOCK) do
        Model._native_geo_kernel_set_transfinite_curve!(
            _model_locked(),tag,num_nodes,mesh_type,coef)
    end
end

function _add_native_geo_point(x,y,z;meshSize=0.0,tag=-1)
    return lock(STATE_LOCK) do
        Model._native_geo_kernel_add_point!(_model_locked(),x,y,z,meshSize,tag)
    end
end
function _add_native_geo_line(first,last;tag=-1)
    return lock(STATE_LOCK) do
        Model._native_geo_kernel_add_line!(_model_locked(),first,last,tag)
    end
end
function _add_native_geo_curve_loop(curves;tag=-1)
    return lock(STATE_LOCK) do
        Model._native_geo_kernel_add_curve_loop!(_model_locked(),curves,tag)
    end
end
function _add_native_geo_plane_surface(loops;tag=-1)
    return lock(STATE_LOCK) do
        Model._native_geo_kernel_add_plane_surface!(_model_locked(),loops,tag)
    end
end

# Built-in GEO sweeps (`geo.extrude`/`revolve`/`twist`/`extrudeBoundaryLayer`)
# act on GEO internals. The buffered kernel is published to the model first,
# the sweep runs on the model exactly like the `.geo` `Extrude` statement, and
# the created Points/Curves join the kernel so later `geo.add_*` calls can
# reference them. Returns Gmsh's `outDimTags`.
function _native_geo_dim_tags(dim_tags,caller)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dimTags must be a vector of (dim, tag) pairs"))
    out=NTuple{2,Int}[]
    for entry in dim_tags
        pair=entry isa Pair ? (first(entry),last(entry)) :
            ((entry isa Tuple || entry isa AbstractVector) && length(entry)==2) ?
            (entry[1],entry[2]) : throw(ArgumentError(
                "$caller: each dimTags entry must be a (dim, tag) pair"))
        all(value->value isa Integer && !(value isa Bool),pair) || throw(ArgumentError(
            "$caller: dimTags entries must be integers"))
        push!(out,(Int(pair[1]),Int(pair[2])))
    end
    return out
end

function _native_geo_publish_locked!(m::GeoModel)
    Model._native_geo_kernel_synchronize!(m;
        context_lc=NATIVE_CONTEXT_LC[],geometry_tolerance=OPTIONS["Geometry.Tolerance"])
    return nothing
end

function _native_geo_after_sweep_locked!(m::GeoModel,points_before,curves_before)
    kernel=m.meshing.native_geo_kernel
    # Coherence merges inside the sweep delete GEO records; drop them from the
    # kernel too so the next synchronization does not resurrect them.
    removed=Tuple{Int,Int}[]
    for tag in keys(kernel.points); haskey(m.points,tag) || push!(removed,(0,tag)); end
    for tag in keys(kernel.curves); haskey(m.curves,tag) || push!(removed,(1,tag)); end
    for tag in keys(kernel.plane_surfaces)
        haskey(m.surfaces,tag) || push!(removed,(2,tag))
    end
    isempty(removed) || Model._native_geo_kernel_forget!(m,removed)
    for tag in sort!(collect(setdiff(keys(m.points),points_before)))
        haskey(kernel.points,tag) || Model._native_geo_kernel_admit_point!(m,tag)
    end
    for tag in sort!(collect(setdiff(keys(m.curves),curves_before)))
        haskey(kernel.curves,tag) || Model._native_geo_kernel_admit_curve!(m,tag)
    end
    _synchronize_native_context_locked!(m)
    _replace_mesh_cache_locked!(nothing;preserve_source_inputs=true)
    return nothing
end

function _native_geo_sweep(kind::Symbol,dim_tags,arguments;numElements=Int[],
                           heights=Float64[],recombine=false,caller)
    entities=_native_geo_dim_tags(dim_tags,caller)
    return lock(STATE_LOCK) do
        m=_model_locked()
        _native_geo_publish_locked!(m)
        params=Model._model_occ_sweep_mesh_params(numElements,heights,recombine,caller)
        points_before=Set(keys(m.points));curves_before=Set(keys(m.curves))
        dims=Int[]
        tags=if kind===:translate
            Model.extrude_entities!(m,entities,Float64.(Tuple(arguments));
                params,caller,out_dims=dims)
        elseif kind===:rotate
            x,y,z,ax,ay,az,angle=arguments
            Model.revolve_entities!(m,entities,(ax,ay,az),(x,y,z),angle;
                params,caller,out_dims=dims)
        else
            x,y,z,dx,dy,dz,ax,ay,az,angle=arguments
            Model.twist_entities!(m,entities,(ax,ay,az),(x,y,z),(dx,dy,dz),angle;
                params,caller,out_dims=dims)
        end
        _native_geo_after_sweep_locked!(m,points_before,curves_before)
        NTuple{2,Int}[(dims[k],Int(tags[k])) for k in eachindex(tags)]
    end
end

function _native_geo_extrude_boundary_layer(dim_tags,numElements=[1],heights=Float64[],
                                            recombine=false,second=false,viewIndex=-1)
    caller="API.model.geo.extrude_boundary_layer"
    entities=_native_geo_dim_tags(dim_tags,caller)
    (second isa Bool && recombine isa Bool) || throw(ArgumentError(
        "$caller: recombine and second must be Bool"))
    viewIndex isa Integer || throw(ArgumentError("$caller: viewIndex must be an integer"))
    layers=Int[]
    for e in numElements
        (e isa Integer && !(e isa Bool) && e>0) || throw(ArgumentError(
            "$caller: numElements entries must be positive integers"))
        push!(layers,Int(e))
    end
    isempty(layers) && throw(ArgumentError("$caller: numElements must not be empty"))
    # Boundary-layer heights are absolute layer thicknesses (upstream keeps
    # them verbatim, unlike the normalized sweep heights).
    hs=isempty(heights) ? Float64[j/length(layers) for j in 1:length(layers)] :
        Float64[Model._finite_scalar(h,caller,"heights entry") for h in heights]
    length(hs)==length(layers) || throw(ArgumentError(
        "$caller: heights must have one entry per numElements layer"))
    params=(layers=layers,heights=hs,scale_last=false,recombine=recombine,
            quad_to_tri=:none,recomb_laterals=false)
    return lock(STATE_LOCK) do
        m=_model_locked()
        _native_geo_publish_locked!(m)
        points_before=Set(keys(m.points));curves_before=Set(keys(m.curves))
        context=GeoExec._GeoNumericContext()
        context.native_context_lc=NATIVE_CONTEXT_LC[]
        blmeta=(boundary_layer_index=second ? 1 : 0,view_index=Int(viewIndex))
        dims=Int[]
        tags=GeoExec._geo_boundary_layer_extrude!(m,entities,params,blmeta,context,caller;
            out_dims=dims)
        out=NTuple{2,Int}[(dims[k],Int(tags[k])) for k in eachindex(tags)]
        _native_geo_after_sweep_locked!(m,points_before,curves_before)
        out
    end
end

# `geo.mesh.*` attributes: built-in entities become model entities at
# synchronization, so the kernel is published first and the attribute is then
# recorded on the model like `mesh.*`.
function _native_geo_mesh_attribute(f)
    return lock(STATE_LOCK) do
        m=_model_locked()
        _native_geo_publish_locked!(m)
        f(m)
        _synchronize_native_context_locked!(m)
        nothing
    end
end

function _native_geo_set_size(dim_tags,size)
    caller="API.model.geo.mesh.set_size"
    entities=_native_geo_dim_tags(dim_tags,caller)
    value=Model._mesh_attr_finite(size,caller,"size")
    lock(STATE_LOCK) do
        kernel=_model_locked().meshing.native_geo_kernel
        for (dim,tag) in entities
            dim==0 || continue
            record=get(kernel.points,tag,nothing)
            record===nothing && throw(ArgumentError(
                "$caller: unknown built-in GEO Point[$tag]"))
            kernel.points[tag]=(record[1],record[2],record[3],value)
        end
    end
    return nothing
end
