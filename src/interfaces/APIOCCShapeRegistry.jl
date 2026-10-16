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

add_spline(points;kwargs...)=_with_model() do m
    Model._model_occ_add_spline!(m,points;tolerance=OPTIONS["Geometry.Tolerance"],kwargs...)
end

add_line(start,stop;kwargs...)=_with_model() do m
    Model._model_occ_add_line!(m,start,stop;kwargs...)
end

add_circle_arc(start,middle,stop;tag=-1,center=true)=_with_model() do m
    Model._model_occ_add_circle_arc!(m,start,middle,stop;tag,center)
end

add_circle(x,y,z,r;tag=-1,angle1=0.0,angle2=2*pi,
           zAxis=Float64[],xAxis=Float64[])=_with_model() do m
    Model._model_occ_add_circle!(m,x,y,z,r;tag,angle1,angle2,
                                 zaxis=zAxis,xaxis=xAxis)
end

add_ellipse_arc(start,center,major,stop;tag=-1)=_with_model() do m
    Model._model_occ_add_ellipse_arc!(m,start,center,major,stop;tag)
end

add_ellipse(x,y,z,rx,ry;tag=-1,angle1=0.0,angle2=2*pi,
            zAxis=Float64[],xAxis=Float64[])=_with_model() do m
    Model._model_occ_add_ellipse!(m,x,y,z,rx,ry;tag,angle1,angle2,
                                  zaxis=zAxis,xaxis=xAxis)
end

add_disk(x,y,z,rx,ry;tag=-1,
         zAxis=Float64[],xAxis=Float64[])=_with_model() do m
    Model._model_occ_add_disk!(m,x,y,z,rx,ry;tag,zaxis=zAxis,xaxis=xAxis)
end

add_wire(curve_tags;tag=-1,checkClosed=false)=_with_model() do m
    Model._model_occ_add_wire!(m,curve_tags;tag,check_closed=checkClosed)
end

add_curve_loop(curve_tags;tag=-1)=_with_model() do m
    Model._model_occ_add_curve_loop!(m,curve_tags;tag)
end

add_curve_loops(curve_tags)=_with_model() do m
    Model._model_occ_curve_loops(m,curve_tags)
end

add_plane_surface(wire_tags;tag=-1)=_with_model() do m
    Model._model_occ_add_plane_surface!(m,wire_tags;tag)
end

add_surface_loop(surface_tags;tag=-1,sewing=false)=_with_model() do m
    Model._model_occ_add_surface_loop!(m,surface_tags;tag,sewing)
end

add_volume(shell_tags;tag=-1)=_with_model() do m
    Model._model_occ_add_volume!(m,shell_tags;tag)
end

add_box(x,y,z,dx,dy,dz;tag=-1)=_with_model() do m
    Model._model_occ_add_box!(m,x,y,z,dx,dy,dz;tag)
end

add_cylinder(x,y,z,dx,dy,dz,r;tag=-1,angle=2*pi)=_with_model() do m
    Model._model_occ_add_cylinder!(m,x,y,z,dx,dy,dz,r;tag,angle)
end

add_cone(x,y,z,dx,dy,dz,r1,r2;tag=-1,angle=2*pi)=_with_model() do m
    Model._model_occ_add_cone!(m,x,y,z,dx,dy,dz,r1,r2;tag,angle)
end

add_sphere(xc,yc,zc,radius;tag=-1,angle1=-pi/2,angle2=pi/2,
           angle3=2*pi)=_with_model() do m
    Model._model_occ_add_sphere!(m,xc,yc,zc,radius;tag,angle1,angle2,angle3)
end

add_torus(x,y,z,r1,r2;tag=-1,angle=2*pi,zAxis=Float64[])=_with_model() do m
    Model._model_occ_add_torus!(m,x,y,z,r1,r2;tag,angle,zAxis)
end

add_wedge(x,y,z,dx,dy,dz;tag=-1,ltx=0.,zAxis=Float64[])=_with_model() do m
    Model._model_occ_add_wedge!(m,x,y,z,dx,dy,dz;tag,ltx,zAxis)
end

add_rectangle(x,y,z,dx,dy;tag=-1,roundedRadius=0.)=_with_model() do m
    Model._model_occ_add_rectangle!(m,x,y,z,dx,dy;tag,
                                    rounded_radius=roundedRadius)
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

translate(dim_tags,dx,dy,dz)=_with_model(m->Model._model_occ_translate!(m,dim_tags,dx,dy,dz))
rotate(dim_tags,x,y,z,ax,ay,az,angle)=_with_model(m->Model._model_occ_rotate!(m,dim_tags,x,y,z,ax,ay,az,angle))
dilate(dim_tags,x,y,z,a,b,c)=_with_model(m->Model._model_occ_dilate!(m,dim_tags,x,y,z,a,b,c))
mirror(dim_tags,a,b,c,d)=_with_model(m->Model._model_occ_mirror!(m,dim_tags,a,b,c,d))
symmetrize(dim_tags,a,b,c,d)=mirror(dim_tags,a,b,c,d)
affine_transform(dim_tags,mat)=_with_model(m->Model._model_occ_affine_transform!(m,dim_tags,mat))
copy(dim_tags)=_with_model(m->Model._model_occ_copy_entities!(m,dim_tags))

get_bounding_box(dim,tag)=_with_model(m->Model._model_occ_get_bounding_box(m,dim,tag))
get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax;dim=-1)=
    _with_model(m->Model._model_occ_entities_in_bounding_box(
        m,xmin,ymin,zmin,xmax,ymax,zmax,dim))
get_curve_loops(surface_tag)=_with_model(m->Model._model_occ_get_curve_loops(m,surface_tag))
get_surface_loops(volume_tag)=_with_model(m->Model._model_occ_get_surface_loops(m,volume_tag))
get_mass(dim,tag)=_with_model(m->Model._model_occ_get_mass(m,dim,tag))
get_center_of_mass(dim,tag)=
    _with_model(m->Model._model_occ_get_center_of_mass(m,dim,tag))
get_matrix_of_inertia(dim,tag)=
    _with_model(m->Model._model_occ_get_matrix_of_inertia(m,dim,tag))
get_closest_entities(x,y,z,dim_tags,n=1)=
    _with_model(m->Model._model_occ_get_closest_entities(m,x,y,z,dim_tags,n))
get_distance(dim1,tag1,dim2,tag2)=
    _with_model(m->Model._model_occ_get_distance(m,dim1,tag1,dim2,tag2))
extrude(dim_tags,dx,dy,dz;numElements=Int[],heights=Float64[],
        recombine=false)=_with_model(
    m->Model._model_occ_extrude(m,dim_tags,dx,dy,dz;
        num_elements=numElements,heights=heights,recombine=recombine))
revolve(dim_tags,x,y,z,ax,ay,az,angle;numElements=Int[],
        heights=Float64[],recombine=false)=_with_model(
    m->Model._model_occ_revolve(m,dim_tags,x,y,z,ax,ay,az,angle;
        num_elements=numElements,heights=heights,recombine=recombine))
fuse(objectDimTags,toolDimTags;tag=-1,removeObject=true,removeTool=true)=
    _with_model(m->Model._model_occ_boolean(
        m,:union,objectDimTags,toolDimTags,tag,removeObject,removeTool,
        "OCC.fuse"))
cut(objectDimTags,toolDimTags;tag=-1,removeObject=true,removeTool=true)=
    _with_model(m->Model._model_occ_boolean(
        m,:difference,objectDimTags,toolDimTags,tag,removeObject,removeTool,
        "OCC.cut"))
intersect(objectDimTags,toolDimTags;tag=-1,removeObject=true,removeTool=true)=
    _with_model(m->Model._model_occ_boolean(
        m,:intersection,objectDimTags,toolDimTags,tag,removeObject,
        removeTool,"OCC.intersect"))
fragment(objectDimTags,toolDimTags;tag=-1,removeObject=true,removeTool=true)=
    _with_model(m->Model._model_occ_boolean(
        m,:fragments,objectDimTags,toolDimTags,tag,removeObject,removeTool,
        "OCC.fragment"))
remove_all_duplicates()=_with_model(
    m->Model._model_occ_remove_all_duplicates(m,"OCC.removeAllDuplicates"))
fillet_2d(edge_tag1,edge_tag2,radius;tag=-1,pointTag=-1,reverse=false)=
    _with_model(m->Model._model_occ_fillet_2d!(
        m,edge_tag1,edge_tag2,radius;tag,point_tag=pointTag,reverse))
chamfer_2d(edge_tag1,edge_tag2,distance1,distance2;tag=-1)=
    _with_model(m->Model._model_occ_chamfer_2d!(
        m,edge_tag1,edge_tag2,distance1,distance2;tag))
offset_curve(curveLoopTag,offset)=_with_model(
    m->Model._model_occ_offset_wire!(m,curveLoopTag,offset;
        caller="OCC.offsetCurve"))
convert_to_nurbs(dimTags)=_with_model(
    m->Model._model_occ_convert_to_nurbs!(m,dimTags,"OCC.convertToNURBS"))
add_thru_sections(wireTags;tag=-1,makeSolid=true,makeRuled=false,
    maxDegree=-1,continuity="",parametrization="",smoothing=false)=
    _with_model(m->Model._model_occ_thru_sections!(m,tag,wireTags;
        make_solid=makeSolid,make_ruled=makeRuled,max_degree=maxDegree,
        continuity=continuity,parametrization=parametrization,
        smoothing=smoothing))

module mesh
    import .....Model
    using ....API: _with_model
    set_size(dim_tags,size)=_with_model(
        m->Model._model_occ_set_size!(m,dim_tags,size))
end

function synchronize()
    _with_model(synchronize_geometry=true) do m
        Model._model_occ_synchronize!(m)
    end
    return nothing
end
end
