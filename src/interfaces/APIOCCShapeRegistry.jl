"""Pending OCC shape operations, separate from the public geometry model."""
module occ
import ....Model
using ...API: _with_model, OPTIONS

add_point(x,y,z;tag=-1,meshSize=0.0)=_with_model() do m
    Model._model_occ_add_point!(m,x,y,z;tag,mesh_size=meshSize)
end

add_bezier(points;kwargs...)=_with_model() do m
    Model._model_occ_add_bezier!(m,points;kwargs...)
end

add_bspline(points;kwargs...)=_with_model() do m
    Model._model_occ_add_bspline!(m,points;tolerance=OPTIONS["Geometry.Tolerance"],kwargs...)
end

add_line(start,stop;kwargs...)=_with_model() do m
    Model._model_occ_add_line!(m,start,stop;kwargs...)
end

add_bezier_surface(points,num_points_u;kwargs...)=_with_model() do m
    Model._model_occ_add_bezier_surface!(m,points,num_points_u;kwargs...)
end

add_bspline_surface(points,num_points_u;kwargs...)=_with_model() do m
    Model._model_occ_add_bspline_surface!(m,points,num_points_u;kwargs...)
end

get_entities(dim=-1)=_with_model(m->Model._model_occ_entities(m,dim))
get_max_tag(dim)=_with_model(m->Model._model_occ_max_tag(m,dim))
set_max_tag(dim,tag)=_with_model(m->Model._model_occ_set_max_tag!(m,dim,tag))
remove(dim_tags,recursive=false)=_with_model(m->Model._model_occ_remove!(m,dim_tags,recursive))

function synchronize()
    _with_model(synchronize_geometry=true) do m
        Model._model_occ_synchronize!(m)
    end
    return nothing
end
end
