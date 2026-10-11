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
