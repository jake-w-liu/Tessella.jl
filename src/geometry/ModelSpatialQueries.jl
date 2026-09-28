const _MODEL_EMPTY_BOUNDS =
    (Inf,Inf,Inf,-Inf,-Inf,-Inf)

@inline function _model_bounds_include(
    bounds::NTuple{6,Float64},point::NTuple{3,Float64})
    return (min(bounds[1],point[1]),min(bounds[2],point[2]),
            min(bounds[3],point[3]),max(bounds[4],point[1]),
            max(bounds[5],point[2]),max(bounds[6],point[3]))
end

@inline function _model_bounds_union(
    first::NTuple{6,Float64},second::NTuple{6,Float64})
    return (min(first[1],second[1]),min(first[2],second[2]),
            min(first[3],second[3]),max(first[4],second[4]),
            max(first[5],second[5]),max(first[6],second[6]))
end

function _model_bounds_checked(
    values,caller::AbstractString,entity::AbstractString)
    bounds=NTuple{6,Float64}(values)
    all(isfinite,bounds) || throw(ArgumentError(
        "$caller: $entity has a non-finite bounding box"))
    (bounds[1]<=bounds[4] && bounds[2]<=bounds[5] &&
     bounds[3]<=bounds[6]) || throw(ErrorException(
        "$caller: $entity has an invalid bounding box; rebuild the model"))
    return bounds
end

function _model_bounds_from_points(
    m::GeoModel,points,caller::AbstractString,entity::AbstractString)
    isempty(points) && throw(ErrorException(
        "$caller: $entity has no boundary Points; rebuild the model"))
    bounds=_MODEL_EMPTY_BOUNDS
    for point in points
        coordinate=get(m.points,point,nothing)
        coordinate===nothing && throw(ErrorException(
            "$caller: $entity references missing Point[$point]; rebuild the model"))
        bounds=_model_bounds_include(bounds,coordinate)
    end
    return _model_bounds_checked(bounds,caller,entity)
end

function _model_bounds_from_surface_mesh(
    mesh::Mesh,caller::AbstractString,entity::AbstractString)
    ntris(mesh)>0 || throw(ErrorException(
        "$caller: $entity has no boundary triangles"))
    bounds=_MODEL_EMPTY_BOUNDS
    node_count=size(mesh.coords,2)
    @inbounds for triangle in axes(mesh.tris,2), local_node in axes(mesh.tris,1)
        node=Int(mesh.tris[local_node,triangle])
        1<=node<=node_count || throw(ErrorException(
            "$caller: $entity has an invalid boundary node; rebuild the model"))
        point=(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
        bounds=_model_bounds_include(bounds,point)
    end
    return _model_bounds_checked(bounds,caller,entity)
end

function _model_bounds_unit_axis(axis,caller::AbstractString,entity::AbstractString)
    scale=max(abs(axis[1]),abs(axis[2]),abs(axis[3]))
    (isfinite(scale) && scale>0) || throw(ErrorException(
        "$caller: $entity has an invalid axis; rebuild the model"))
    scaled=(axis[1]/scale,axis[2]/scale,axis[3]/scale)
    magnitude=hypot(scaled...)
    (isfinite(magnitude) && magnitude>0) || throw(ErrorException(
        "$caller: $entity has an invalid axis; rebuild the model"))
    return (scaled[1]/magnitude,scaled[2]/magnitude,scaled[3]/magnitude)
end

function _model_cylinder_bounds(spec,caller::AbstractString,entity::AbstractString)
    center=spec.center
    height=spec.height
    radius=spec.radius
    (isfinite(height) && height>0 && isfinite(radius) && radius>0) ||
        throw(ErrorException(
            "$caller: $entity has invalid cylinder parameters; rebuild the model"))
    axis=_model_bounds_unit_axis(spec.axis,caller,entity)
    endpoint=ntuple(index->center[index]+height*axis[index],3)
    # The transverse norm avoids cancellation in sqrt(1-axis[index]^2)
    # when the axis is nearly parallel to a coordinate direction.
    radial=ntuple(index->
        radius*hypot(axis[mod1(index+1,3)],axis[mod1(index+2,3)]),3)
    return _model_bounds_checked((
        ntuple(index->min(center[index],endpoint[index])-radial[index],3)...,
        ntuple(index->max(center[index],endpoint[index])+radial[index],3)...),
        caller,entity)
end

function _model_cone_bounds(spec,caller::AbstractString,entity::AbstractString)
    center=spec.center
    height=spec.height
    r1=spec.r1
    r2=spec.r2
    (isfinite(height) && height>0 && isfinite(r1) && isfinite(r2) &&
     r1>=0 && r2>=0 && (r1>0 || r2>0)) || throw(ErrorException(
        "$caller: $entity has invalid cone parameters; rebuild the model"))
    axis=_model_bounds_unit_axis(spec.axis,caller,entity)
    endpoint=ntuple(index->center[index]+height*axis[index],3)
    orthogonal=ntuple(index->
        hypot(axis[mod1(index+1,3)],axis[mod1(index+2,3)]),3)
    return _model_bounds_checked((
        ntuple(index->min(center[index]-r1*orthogonal[index],
                          endpoint[index]-r2*orthogonal[index]),3)...,
        ntuple(index->max(center[index]+r1*orthogonal[index],
                          endpoint[index]+r2*orthogonal[index]),3)...),
        caller,entity)
end

function _model_volume_bounds(
    m::GeoModel,tag::Int,caller::AbstractString)
    entity="Volume[$tag]"
    encodings=(haskey(m.box_extents,tag),haskey(m.cylinders,tag),
               haskey(m.spheres,tag),haskey(m.cones,tag),
               haskey(m.booleans,tag))
    # Materialized primitives keep their compact encoding alongside the shell
    # topology — legal only while the materialized entities still satisfy it
    # (corners for a box, rim/pole satisfaction for curved solids). A Boolean
    # encoding likewise coexists with its materialized result shell: the
    # encoding drives the mesh path while the entities serve topology queries.
    shelled=!isempty(m.volumes[tag])
    materialized=any(encodings[1:4]) && shelled
    count(identity,encodings)+Int(shelled && !any(encodings))<=1 ||
        throw(ErrorException(
            "$caller: $entity has multiple native encodings; rebuild the model"))
    if materialized && encodings[1]
        x0,y0,z0,dx,dy,dz=m.box_extents[tag]
        scale=max(1.0,abs(x0),abs(y0),abs(z0),abs(dx),abs(dy),abs(dz))
        corners=NTuple{3,Float64}[
            (x0+ix*dx,y0+iy*dy,z0+iz*dz) for ix in (0,1),iy in (0,1),iz in (0,1)]
        owned=[m.points[point] for point in _model_volume_owned_points(m,tag)]
        (length(owned)==8 && all(
            corner->any(p->_points_close(p,corner,1e-9*scale),owned),
            corners)) || throw(ErrorException(
                "$caller: $entity has inconsistent box encoding and shell " *
                "topology; rebuild the model"))
    elseif materialized
        _materialized_curved_consistent(m,tag) || throw(ErrorException(
            "$caller: $entity has inconsistent primitive encoding and shell " *
            "topology; rebuild the model"))
    end
    haskey(m.booleans,tag)==haskey(m.boolean_operands,tag) ||
        throw(ErrorException(
            "$caller: $entity has inconsistent Boolean snapshot ownership; " *
            "rebuild the model"))

    if encodings[1]
        x,y,z,dx,dy,dz=m.box_extents[tag]
        (isfinite(dx) && isfinite(dy) && isfinite(dz) &&
         dx>0 && dy>0 && dz>0) || throw(ErrorException(
            "$caller: $entity has invalid box parameters; rebuild the model"))
        return _model_bounds_checked(
            (x,y,z,x+dx,y+dy,z+dz),caller,entity)
    elseif encodings[2]
        return _model_cylinder_bounds(m.cylinders[tag],caller,entity)
    elseif encodings[3]
        sphere=m.spheres[tag]
        center=sphere.center
        radius=sphere.radius
        (isfinite(radius) && radius>0) || throw(ErrorException(
            "$caller: $entity has invalid sphere parameters; rebuild the model"))
        return _model_bounds_checked((
            center[1]-radius,center[2]-radius,center[3]-radius,
            center[1]+radius,center[2]+radius,center[3]+radius),caller,entity)
    elseif encodings[4]
        return _model_cone_bounds(m.cones[tag],caller,entity)
    elseif encodings[5]
        surface=_boolean_result_surface(m,tag,caller)
        return _model_bounds_from_surface_mesh(surface,caller,entity)
    end

    shelled || throw(ArgumentError(
        "$caller: $entity has no native solid encoding"))
    bounds=_MODEL_EMPTY_BOUNDS
    for curve in _model_boundary_curves(m,3,tag,caller,entity)
        bounds=_model_bounds_union(
            bounds,_model_curve_bounding_box(m,curve,caller))
    end
    return _model_bounds_checked(bounds,caller,entity)
end

# Exact bounds for one explicit curve: endpoint bounds for a Line, the true arc
# bounding box (endpoints plus axis-aligned extrema) for curved kinds.
function _model_curve_bounding_box(
    m::GeoModel,curve::Int,caller::AbstractString)
    occ=_occ_geometry_checked(m,curve,caller)
    occ===nothing || return _model_bounds_checked(
        occ.occ===:circle ? _occ_circle_bounding_box(occ) :
        occ.occ===:line ? _model_bounds_from_points(
            m,m.curves[curve],caller,"Curve[$curve]") :
        begin
            point=m.points[m.curves[curve][1]]
            (point[1],point[2],point[3],point[1],point[2],point[3])
        end,caller,"Curve[$curve]")
    if _curve_type(m,curve)==:line
        return _model_bounds_from_points(
            m,m.curves[curve],caller,"Curve[$curve]")
    end
    _curve_type(m,curve) in _SPLINE_CURVE_TYPES && return _model_bounds_checked(
        _spline_bounding_box(_spline_geometry(m,curve,caller)),
        caller,"Curve[$curve]")
    return _model_bounds_checked(
        _arc_bounding_box(_arc_geometry(m,curve,caller)),
        caller,"Curve[$curve]")
end

# Boundary curves of an explicit surface or volume, unsigned, deduplicated.
function _model_boundary_curves(m::GeoModel,dimension::Int,tag::Int,
                                caller::AbstractString,entity::AbstractString)
    curves=Set{Int}()
    if dimension==2
        for curve in _model_surface_curves(m,tag)
            push!(curves,curve)
        end
    else
        shells=m.volumes[tag]
        isempty(shells) && throw(ArgumentError(
            "$caller: $entity has no explicit surface-loop topology; " *
            "define it from Surface Loop entities or select explicit Points"))
        for shell in shells
            haskey(m.surface_loops,shell) || throw(ArgumentError(
                "$caller: $entity references unknown Surface Loop[$shell]"))
            for signed_surface in m.surface_loops[shell]
                surface=abs(signed_surface)
                haskey(m.surfaces,surface) || throw(ArgumentError(
                    "$caller: $entity references unknown Surface[$surface]"))
                for curve in _model_surface_curves(m,surface)
                    push!(curves,curve)
                end
            end
        end
    end
    isempty(curves) && throw(ErrorException(
        "$caller: $entity has no boundary curves; rebuild the model"))
    return curves
end

function _model_entity_bounding_box(
    m::GeoModel,dimension::Int,tag::Int,caller::AbstractString)
    entities=_model_entity_dictionary(m,dimension)
    haskey(entities,tag) || throw(ArgumentError(
        "$caller: unknown entity ($dimension,$tag)"))
    if dimension==0
        point=m.points[tag]
        return _model_bounds_checked(
            (point[1],point[2],point[3],point[1],point[2],point[3]),
            caller,"Point[$tag]")
    elseif dimension==1
        return _model_curve_bounding_box(m,tag,caller)
    elseif dimension==2
        surface_geometry=get(m.surface_geometry,tag,nothing)
        # A sphere's wire only spans a single meridian; its face bbox is the
        # analytic ball box. Cylinder/cone faces reduce to their end circles,
        # which the boundary-curve union already covers exactly.
        if surface_geometry!==nothing && hasproperty(surface_geometry,:occ) &&
                surface_geometry.occ===:sphere
            C=surface_geometry.center; r=surface_geometry.radius
            return _model_bounds_checked(
                (C[1]-r,C[2]-r,C[3]-r,C[1]+r,C[2]+r,C[3]+r),
                caller,"Surface[$tag]")
        end
        # A partial torus's extrema can sit on interior meridians no wire
        # covers; the face box is the exact one-variable sweep maximization.
        if surface_geometry!==nothing && hasproperty(surface_geometry,:occ) &&
                surface_geometry.occ===:torus
            return _model_bounds_checked(
                _occ_torus_bounding_box(surface_geometry),
                caller,"Surface[$tag]")
        end
        bounds=_MODEL_EMPTY_BOUNDS
        for curve in _model_boundary_curves(
            m,2,tag,caller,"Surface[$tag]")
            bounds=_model_bounds_union(
                bounds,_model_curve_bounding_box(m,curve,caller))
        end
        return _model_bounds_checked(bounds,caller,"Surface[$tag]")
    end
    return _model_volume_bounds(m,tag,caller)
end

"""
    model_bounding_box(model, dim, tag) -> NTuple{6,Float64}

Return `(xmin, ymin, zmin, xmax, ymax, zmax)` for one existing positive-tag
entity. Passing `dim=-1, tag=-1` returns the union over the whole nonempty model.
Explicit straight-edge topology and native primitive bounds are analytical;
Boolean bounds belong to the operation-time result snapshot.
"""
function model_bounding_box(m::GeoModel,dim,tag)
    caller="model_bounding_box"
    if dim isa Integer && !(dim isa Bool) && tag isa Integer &&
            !(tag isa Bool) && dim==-1 && tag==-1
        entities=model_entities(m)
        isempty(entities) && throw(ArgumentError(
            "$caller: the model has no entities"))
        bounds=_MODEL_EMPTY_BOUNDS
        for (dimension,entity_tag) in entities
            bounds=_model_bounds_union(bounds,_model_entity_bounding_box(
                m,dimension,entity_tag,caller))
        end
        return _model_bounds_checked(bounds,caller,"model")
    end
    dimension=_dimension(dim,caller)
    entity_tag=_tag(tag,caller,dimension)
    entity_tag>=0 || throw(ArgumentError(
        "$caller: entity tag must be non-negative"))
    return _model_entity_bounding_box(m,dimension,entity_tag,caller)
end

function _model_spatial_coordinate(value,caller::AbstractString,name::AbstractString)
    value isa Real || throw(ArgumentError(
        "$caller: $name must be a real number"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    coordinate=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $name must be Float64-representable"))
    end
    isfinite(coordinate) || throw(ArgumentError(
        "$caller: $name must be finite"))
    return coordinate
end

"""
    model_entities_in_bounding_box(model, xmin, ymin, zmin,
                                   xmax, ymax, zmax, dim=-1)

Return detached, sorted entities whose complete analytical bounding box is inside
the supplied finite box. `dim=-1` selects all dimensions; reversed bounds simply
select no entities. Embeddings do not enlarge their target entity.
"""
function model_entities_in_bounding_box(
    m::GeoModel,xmin,ymin,zmin,xmax,ymax,zmax,dim=-1)
    caller="model_entities_in_bounding_box"
    dimension=_query_dimension(dim,caller)
    values=(xmin,ymin,zmin,xmax,ymax,zmax)
    names=("xmin","ymin","zmin","xmax","ymax","zmax")
    box=ntuple(index->_model_spatial_coordinate(
        values[index],caller,names[index]),6)
    (box[1]<=box[4] && box[2]<=box[5] && box[3]<=box[6]) ||
        return Tuple{Int,Int}[]
    entities=model_entities(m,dimension)
    result=Tuple{Int,Int}[]
    sizehint!(result,length(entities))
    for (entity_dimension,tag) in entities
        bounds=_model_entity_bounding_box(
            m,entity_dimension,tag,caller)
        bounds[1]>=box[1] && bounds[2]>=box[2] && bounds[3]>=box[3] &&
        bounds[4]<=box[4] && bounds[5]<=box[5] && bounds[6]<=box[6] &&
            push!(result,(entity_dimension,tag))
    end
    return result
end

# ── Entity distance ──────────────────────────────────────────────────────────
#
# `model_distance` mirrors `gmsh.model.getDistance`: the minimum Euclidean
# distance between two entities plus one attaining point on each. Distance to
# a Volume is the distance to the SOLID — zero whenever the other entity
# meets or lies inside it — matching OCC `BRepExtrema_DistShapeShape` solid
# semantics (`getDistance` on OCC volumes returns the inner point for
# containment). Boundary distances reuse the verified projectors
# (`_occ_curve_project`/`_occ_curve_closest`, `_spline_closest`,
# `_arc_closest`, the OCC face extrema with the `BRepClass` trim classifier,
# exact rational line/plane projection, discrete element closests). Generic
# pairs minimize through deterministic multi-seed scans plus alternating-
# projection refinement; every accepted pair is a mutual-projection fixed
# point, i.e. a true local minimum, and the seed grids cover every basin for
# the supported kinds (line↔line additionally uses the segment closed form).
# Solid containment uses the compact primitive encodings while they satisfy
# the materialized shell and a `_BREP_RAY_DIRS` ray-parity argument over the
# shell faces otherwise, in the spirit of `_brep_point_in_solid`. Ambiguous
# directions retry a fixed direction set; a query no direction resolves
# raises rather than guessing.

const _DIST_SCAN_CURVE=33
const _DIST_SCAN_FACE=9
const _DIST_ALT_ITERS=32
const _DIST_GN_ITERS=30
const _DIST_SEED_CAP=16
const _DIST_RAY_SAMPLE=192

# Surface tags bounding a volume: the declared boundary of a discrete volume,
# else the signed shell faces of a shelled volume — sorted, deduplicated.
function _dist_volume_faces(m::GeoModel,tag::Int,caller::AbstractString)
    record=get(m.discrete,(3,tag),nothing)
    if record!==nothing
        faces=Int[]
        for (bdim,btag) in record.boundary
            bdim==2 || throw(ArgumentError(
                "$caller: discrete Volume[$tag] boundary ($bdim,$btag) is " *
                "not a surface"))
            (haskey(m.surfaces,btag) ||
             haskey(m.discrete,(2,btag))) || throw(ArgumentError(
                "$caller: discrete Volume[$tag] references unknown " *
                "Surface[$btag]"))
            push!(faces,btag)
        end
        isempty(faces) && throw(ArgumentError(
            "$caller: discrete Volume[$tag] has no declared boundary to " *
            "measure against"))
        return sort!(unique!(faces))
    end
    shells=get(m.volumes,tag,Int[])
    isempty(shells) && return Int[]
    faces=Set{Int}()
    for shell in shells
        haskey(m.surface_loops,shell) || throw(ArgumentError(
            "$caller: Volume[$tag] references unknown Surface Loop[$shell]"))
        for signed_surface in m.surface_loops[shell]
            surface=abs(signed_surface)
            (haskey(m.surfaces,surface) ||
             haskey(m.discrete,(2,surface))) || throw(ArgumentError(
                "$caller: Volume[$tag] references unknown Surface[$surface]"))
            push!(faces,surface)
        end
    end
    return sort!(collect(faces))
end

# 3-D coordinate of an existing Point entity, finite-checked.
function _dist_point_coord(m::GeoModel,tag::Int,caller::AbstractString)
    record=get(m.discrete,(0,tag),nothing)
    coordinate=if record!==nothing
        isempty(record.node_coords) && throw(ArgumentError(
            "$caller: discrete Point[$tag] has no node"))
        NTuple{3,Float64}(record.node_coords[:,1])
    else
        m.points[tag]
    end
    all(isfinite,coordinate) || throw(ArgumentError(
        "$caller: Point[$tag] has a non-finite coordinate"))
    return coordinate
end

# Corner coordinates of a discrete entity's dim-`dim` elements (lin2/tri3/quad
# corners; quads split like `_discrete_element_frames`).
function _dist_discrete_frames(m::GeoModel,record::DiscreteEntity,dim::Int,
                               caller::AbstractString,entity::AbstractString)
    coords=_model_record_node_coords(m)
    frames=Vector{NTuple{3,Float64}}[]
    for (index,msh_type) in enumerate(record.element_types)
        msh_dimension(msh_type)!=dim && continue
        nodes=record.element_nodes[index]
        count=dim==1 ? 2 : msh_family(msh_type)===:qua ? 4 : 3
        length(nodes)<count && continue
        points=Vector{NTuple{3,Float64}}(undef,count)
        for i in 1:count
            coordinate=get(coords,nodes[i],nothing)
            coordinate===nothing && throw(ErrorException(
                "$caller: $entity references missing node $(nodes[i]); " *
                "rebuild the model"))
            points[i]=coordinate
        end
        if dim==2 && msh_family(msh_type)===:qua
            push!(frames,[points[1],points[2],points[3]])
            push!(frames,[points[1],points[3],points[4]])
        else
            push!(frames,points)
        end
    end
    isempty(frames) && throw(ArgumentError(
        "$caller: $entity has no elements"))
    return frames
end

# Parameter interval of a native (non-discrete) curve.
function _dist_curve_bounds(m::GeoModel,tag::Int,caller::AbstractString)
    occ=_occ_geometry_checked(m,tag,caller)
    occ!==nothing && return occ.t0,occ.t1
    kind=_curve_type(m,tag)
    if kind in _SPLINE_CURVE_TYPES
        spline=_spline_geometry(m,tag,caller)
        return spline.ubeg,spline.uend
    end
    return 0.0,1.0
end

# Closest point on a native or discrete curve → (sqd, xyz).
function _dist_curve_project(m::GeoModel,tag::Int,q::NTuple{3,Float64},
                             lc::Float64,caller::AbstractString)
    record=get(m.discrete,(1,tag),nothing)
    if record!==nothing
        best=Inf;best_point=(0.0,0.0,0.0)
        for points in _dist_discrete_frames(
            m,record,1,caller,"discrete Curve[$tag]")
            t,dist=_segment_closest(q,points[1],points[2])
            if dist<best
                best=dist
                a,b=points[1],points[2]
                best_point=(a[1]+t*(b[1]-a[1]),a[2]+t*(b[2]-a[2]),
                            a[3]+t*(b[3]-a[3]))
            end
        end
        return best,best_point
    end
    occ=_occ_geometry_checked(m,tag,caller)
    if occ!==nothing
        proj=_occ_curve_project(m,occ,tag,q,caller)
        t,pt=proj===nothing ? _occ_curve_closest(m,occ,tag,q,caller) : proj
        return _sqdist(q,pt),pt
    end
    kind=_curve_type(m,tag)
    if kind==:line
        line=_model_line_geometry(m,tag,caller)
        exact,_=_model_line_parameter_exact(line,q)
        exact=clamp(exact,zero(exact),one(exact))
        pt=_model_line_point(line,exact,caller,1)
        return _sqdist(q,pt),pt
    elseif kind in _SPLINE_CURVE_TYPES
        _,pt=_spline_closest(_spline_geometry(m,tag,caller),q)
        return _sqdist(q,pt),pt
    end
    _,pt=_arc_closest(_arc_geometry(m,tag,caller),q)
    return _sqdist(q,pt),pt
end

# Boundary-curve closest distances of a surface's wires → (sqd, xyz).
function _dist_face_boundary_project(m::GeoModel,stag::Int,
                                     q::NTuple{3,Float64},lc::Float64,
                                     caller::AbstractString)
    curves=sort!(unique!(_model_surface_curves(m,stag)))
    isempty(curves) && return Inf,(0.0,0.0,0.0)
    best=Inf;best_point=(0.0,0.0,0.0)
    for curve in curves
        dist,pt=_dist_curve_project(m,curve,q,lc,caller)
        dist<best && (best=dist;best_point=pt)
    end
    return best,best_point
end

# Sampled-polyline boundary polygons of an explicit Plane in its 2-D
# projection axes — the same flattening `_model_plane_boundary_polygons`
# performs, except curved boundary curves are sampled rather than rejected.
function _dist_plane_polygons(m::GeoModel,stag::Int,plane,
                              caller::AbstractString)
    first_axis,second_axis=plane.projection
    out=Vector{NTuple{2,Float64}}[]
    for loop in m.surfaces[stag]
        polygon=NTuple{2,Float64}[]
        for signed_curve in m.loops[loop]
            curve=abs(signed_curve)
            t0,t1=_dist_curve_bounds(m,curve,caller)
            kind=_curve_type(m,curve)
            occ=_occ_geometry_checked(m,curve,caller)
            samples=(occ===nothing && kind==:line) ? 2 : 17
            for i in 0:samples-1
                t=fma(t1-t0,signed_curve>0 ? i/(samples-1) :
                          1.0-i/(samples-1),t0)
                point=_model_curve_point(m,curve,t,caller)
                push!(polygon,(point[first_axis],point[second_axis]))
            end
        end
        push!(out,polygon)
    end
    return out
end

# :in/:out/:on state of a 2-D point against sampled loop polygons:
# inside the outer loop but outside every hole loop → :in; on any loop → :on.
function _dist_position_in_polygons(point::NTuple{2,Float64},
                                    polygons)
    isempty(polygons) && return :out
    outer=_model_loop_position2(point,first(polygons))
    outer==0 && return :out
    outer==2 && return :on
    for polygon in Iterators.drop(polygons,1)
        position=_model_loop_position2(point,polygon)
        position==2 && return :on
        position==1 && return :out
    end
    return :in
end

# :in/:out/:on state of an on-plane point against the sampled polygons.
function _dist_plane_state(m::GeoModel,stag::Int,plane,
                           xyz::NTuple{3,Float64},caller::AbstractString)
    first_axis,second_axis=plane.projection
    point=(xyz[first_axis],xyz[second_axis])
    polygons=_dist_plane_polygons(m,stag,plane,caller)
    isempty(polygons) && throw(ErrorException(
        "$caller: Plane[$stag] has no boundary loops; rebuild the model"))
    return _dist_position_in_polygons(point,polygons)
end

# Boundary-curve samples of a face — seeds that always lie on the trim.
function _dist_face_edge_samples(m::GeoModel,stag::Int,
                                 caller::AbstractString)
    points=NTuple{3,Float64}[]
    for curve in sort!(unique!(_model_surface_curves(m,stag)))
        t0,t1=_dist_curve_bounds(m,curve,caller)
        for i in 0:4
            push!(points,_model_curve_point(
                m,curve,fma(t1-t0,i/4.0,t0),caller))
        end
    end
    return points
end

# Multi-seed Gauss–Newton closest point on a :ruled/:tric patch over [0,1]².
function _dist_ruled_closest(m::GeoModel,stag::Int,q::NTuple{3,Float64},
                             caller::AbstractString)
    best=Inf;best_point=(0.0,0.0,0.0);found=false
    for init_u in 1:5, init_v in 1:5
        u=(init_u-0.5)/5.0;v=(init_v-0.5)/5.0
        for _ in 1:_DIST_GN_ITERS
            point=_ruled_surface_point(m,stag,u,v,caller)
            du=_ruled_surface_d1(m,stag,u,v,1,caller)
            dv=_ruled_surface_d1(m,stag,u,v,2,caller)
            r=_arc_sub(point,q)
            a11=_dot3(du,du);a12=_dot3(du,dv);a22=_dot3(dv,dv)
            determinant=a11*a22-a12*a12
            abs(determinant)<=1e-30 && break
            b1=-_dot3(du,r);b2=-_dot3(dv,r)
            step_u=(a22*b1-a12*b2)/determinant
            step_v=(a11*b2-a12*b1)/determinant
            next_u=clamp(u+step_u,0.0,1.0)
            next_v=clamp(v+step_v,0.0,1.0)
            step=(next_u-u)^2+(next_v-v)^2
            u=next_u;v=next_v
            step<=1e-24 && break
        end
        point=_ruled_surface_point(m,stag,u,v,caller)
        dist=_sqdist(point,q)
        if dist<best
            best=dist;best_point=point;found=true
        end
    end
    found || return Inf,(0.0,0.0,0.0)
    return best,best_point
end

# Closest point on a native or discrete trimmed surface → (sqd, xyz).
function _dist_face_project(m::GeoModel,stag::Int,q::NTuple{3,Float64},
                            lc::Float64,caller::AbstractString)
    record=get(m.discrete,(2,stag),nothing)
    if record!==nothing
        best=Inf;best_point=(0.0,0.0,0.0)
        for points in _dist_discrete_frames(
            m,record,2,caller,"discrete Surface[$stag]")
            weights,dist=_triangle_closest(
                q,points[1],points[2],points[3])
            if dist<best
                best=dist
                best_point=ntuple(axis->sum(
                    k->weights[k]*points[k][axis],1:3),3)
            end
        end
        return best,best_point
    end
    geometry=get(m.surface_geometry,stag,nothing)
    if geometry!==nothing && hasproperty(geometry,:occ)
        # The unrestricted argmin inside the trim is the trimmed-face answer
        # outright; when it falls outside, in-trim extrema (interior local
        # minima) compete with the boundary curves for the constrained min.
        u0,v0,xyz0=_occ_surface_closest(geometry,q,lc)
        state=_occ_face_state(m,stag,geometry,(u0,v0),_OCC_CONFUSION)
        (state===:in || state===:on) && return _sqdist(q,xyz0),xyz0
        umin,umax,vmin,vmax=_occ_uvbounds(m,stag,geometry)
        acc=_occ_surface_extrema(
            geometry,q,umin,umax,vmin,vmax,_OCC_CONFUSION,_OCC_CONFUSION)
        best=Inf;best_point=(0.0,0.0,0.0)
        if acc!==nothing
            for (u,v,xyz,sqd) in acc
                state=_occ_face_state(
                    m,stag,geometry,(u,v),_OCC_CONFUSION)
                (state===:in || state===:on) || continue
                sqd<best && (best=sqd;best_point=xyz)
            end
        end
        boundary,boundary_point=_dist_face_boundary_project(
            m,stag,q,lc,caller)
        boundary<best && (best=boundary;best_point=boundary_point)
        isfinite(best) || throw(ErrorException(
            "$caller: Surface[$stag] admits no closest point"))
        return best,best_point
    end
    kind=_surface_type(m,stag)
    if kind===:plane
        plane=_model_plane_frame(m,stag,caller)
        exact=_model_plane_parameters_exact(plane,q)
        xyz=_model_plane_point(plane,exact[1],exact[2],caller,1)
        state=_dist_plane_state(m,stag,plane,xyz,caller)
        if state!==:out
            return _sqdist(q,xyz),xyz
        end
        best,best_point=_dist_face_boundary_project(m,stag,q,lc,caller)
        isfinite(best) || throw(ErrorException(
            "$caller: Plane[$stag] admits no closest point"))
        return best,best_point
    elseif kind in (:ruled,:tric)
        best,best_point=_dist_ruled_closest(m,stag,q,caller)
        boundary,boundary_point=_dist_face_boundary_project(
            m,stag,q,lc,caller)
        boundary<best && (best=boundary;best_point=boundary_point)
        isfinite(best) || throw(ErrorException(
            "$caller: Surface[$stag] admits no closest point"))
        return best,best_point
    end
    throw(ArgumentError(
        "$caller: Surface[$stag] has unsupported kind $kind"))
end

# Unified entity projector: (sqd, closest point on the entity).
# dim 3 resolves to solid distance (0, q) when q is inside.
function _dist_entity_project(m::GeoModel,dim::Int,tag::Int,
                              q::NTuple{3,Float64},lc::Float64,
                              caller::AbstractString)
    dim==0 && return _sqdist(q,_dist_point_coord(m,tag,caller)),
                     _dist_point_coord(m,tag,caller)
    dim==1 && return _dist_curve_project(m,tag,q,lc,caller)
    dim==2 && return _dist_face_project(m,tag,q,lc,caller)
    return _dist_volume_project(m,tag,q,lc,caller)
end

# Parameter-space sample points of a curve or face → (points, grid_dims).
# grid_dims is (nu,nv) for parametrized faces, nothing otherwise.
function _dist_scan(m::GeoModel,dim::Int,tag::Int,caller::AbstractString)
    dim==0 && return [_dist_point_coord(m,tag,caller)],nothing
    record=get(m.discrete,(dim,tag),nothing)
    if record!==nothing
        points=NTuple{3,Float64}[]
        if dim==1
            for frame in _dist_discrete_frames(
                m,record,1,caller,"discrete Curve[$tag]")
                for i in 0:4
                    t=i/4
                    push!(points,ntuple(axis->
                        frame[1][axis]+t*(frame[2][axis]-frame[1][axis]),3))
                end
            end
        else
            for frame in _dist_discrete_frames(
                m,record,2,caller,"discrete Surface[$tag]")
                a,b,c=frame
                push!(points,a);push!(points,b);push!(points,c)
                push!(points,ntuple(axis->(a[axis]+b[axis]+c[axis])/3,3))
            end
        end
        return points,nothing
    end
    if dim==1
        t0,t1=_dist_curve_bounds(m,tag,caller)
        points=NTuple{3,Float64}[
            _model_curve_point(m,tag,fma(
                t1-t0,i/(_DIST_SCAN_CURVE-1),t0),caller)
            for i in 0:_DIST_SCAN_CURVE-1]
        return points,nothing
    end
    geometry=get(m.surface_geometry,tag,nothing)
    if geometry!==nothing && hasproperty(geometry,:occ)
        lo,hi=_occ_surface_bounds(geometry)
        umin,vmin=lo[1],lo[2];umax,vmax=hi[1],hi[2]
        points=NTuple{3,Float64}[];all_in=true
        for j in 0:_DIST_SCAN_FACE-1, i in 0:_DIST_SCAN_FACE-1
            uv=(fma(umax-umin,i/(_DIST_SCAN_FACE-1),umin),
                fma(vmax-vmin,j/(_DIST_SCAN_FACE-1),vmin))
            state=_occ_face_state(m,tag,geometry,uv,_OCC_CONFUSION)
            if state===:in || state===:on
                push!(points,_occ_surface_point(geometry,uv[1],uv[2]))
            else
                all_in=false
            end
        end
        if !all_in
            append!(points,_dist_face_edge_samples(m,tag,caller))
        end
        isempty(points) && throw(ErrorException(
            "$caller: Surface[$tag] admits no in-trim samples"))
        return points,all_in ? (_DIST_SCAN_FACE,_DIST_SCAN_FACE) : nothing
    end
    kind=_surface_type(m,tag)
    if kind in (:ruled,:tric)
        points=NTuple{3,Float64}[
            _ruled_surface_point(m,tag,
                i/(_DIST_SCAN_FACE-1),j/(_DIST_SCAN_FACE-1),caller)
            for j in 0:_DIST_SCAN_FACE-1 for i in 0:_DIST_SCAN_FACE-1]
        return points,(_DIST_SCAN_FACE,_DIST_SCAN_FACE)
    end
    plane=_model_plane_frame(m,tag,caller)
    (umin,vmin),(umax,vmax)=_model_plane_parameter_bounds(plane,caller)
    polygons=_dist_plane_polygons(m,tag,plane,caller)
    isempty(polygons) && throw(ErrorException(
        "$caller: Plane[$tag] has no boundary loops; rebuild the model"))
    first_axis,second_axis=plane.projection
    points=NTuple{3,Float64}[];all_in=true
    for j in 0:_DIST_SCAN_FACE-1, i in 0:_DIST_SCAN_FACE-1
        xyz=_model_plane_point(plane,
            fma(umax-umin,i/(_DIST_SCAN_FACE-1),umin),
            fma(vmax-vmin,j/(_DIST_SCAN_FACE-1),vmin),caller,
            i+1+j*_DIST_SCAN_FACE)
        state=_dist_position_in_polygons(
            (xyz[first_axis],xyz[second_axis]),polygons)
        if state===:in || state===:on
            push!(points,xyz)
        else
            all_in=false
        end
    end
    if !all_in
        append!(points,_dist_face_edge_samples(m,tag,caller))
    end
    isempty(points) && throw(ErrorException(
        "$caller: Plane[$tag] admits no in-trim samples"))
    return points,all_in ? (_DIST_SCAN_FACE,_DIST_SCAN_FACE) : nothing
end

# Local-minimum indices of a scan-distance list (1-D parameter order).
function _dist_local_minima_1d(distances::Vector{Float64})
    out=Int[1]
    for i in 2:length(distances)-1
        distances[i]<=min(distances[i-1],distances[i+1]) && push!(out,i)
    end
    length(distances)>=2 && push!(out,length(distances))
    return out
end

# Local-minimum linear indices of a nu×nv scan-distance matrix.
function _dist_local_minima_2d(distances::Matrix{Float64})
    nu,nv=size(distances)
    out=Int[]
    for j in 1:nv, i in 1:nu
        value=distances[i,j]
        local_min=true
        for (ii,jj) in ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
            (1<=ii<=nu && 1<=jj<=nv) || continue
            value>distances[ii,jj] && (local_min=false;break)
        end
        local_min && push!(out,i+(j-1)*nu)
    end
    return out
end

# Alternating-projection refinement of a candidate pair. Each projection
# step cannot increase the pair distance, so the sequence is monotone;
# iteration stops on stall or the iteration cap.
function _dist_alternate(m::GeoModel,dim_a::Int,tag_a::Int,
                         dim_b::Int,tag_b::Int,
                         p_a::NTuple{3,Float64},p_b::NTuple{3,Float64},
                         lc::Float64,caller::AbstractString)
    best_sqd=_sqdist(p_a,p_b);best_a=p_a;best_b=p_b
    for _ in 1:_DIST_ALT_ITERS
        _,next_b=_dist_entity_project(m,dim_b,tag_b,best_a,lc,caller)
        _,next_a=_dist_entity_project(m,dim_a,tag_a,next_b,lc,caller)
        sqd=_sqdist(next_a,next_b)
        best_sqd-sqd<=1e-18*max(1.0,best_sqd) && break
        best_sqd=sqd;best_a=next_a;best_b=next_b
    end
    return best_sqd,best_a,best_b
end

# Scan + alternate minimization between dim>=1 entities → (sqd, pa, pb).
function _dist_minimize(m::GeoModel,dim_a::Int,tag_a::Int,
                        dim_b::Int,tag_b::Int,lc::Float64,
                        caller::AbstractString)
    samples_a,dims_a=_dist_scan(m,dim_a,tag_a,caller)
    vals_a=Tuple{Float64,NTuple{3,Float64},NTuple{3,Float64}}[
        let (sqd,pb)=_dist_entity_project(m,dim_b,tag_b,pa,lc,caller)
            (sqd,pa,pb)
        end for pa in samples_a]
    samples_b,dims_b=_dist_scan(m,dim_b,tag_b,caller)
    vals_b=Tuple{Float64,NTuple{3,Float64},NTuple{3,Float64}}[
        let (sqd,pa)=_dist_entity_project(m,dim_a,tag_a,pb,lc,caller)
            (sqd,pa,pb)
        end for pb in samples_b]
    dist_a=[v[1] for v in vals_a];dist_b=[v[1] for v in vals_b]
    indices_a=dims_a===nothing ? _dist_local_minima_1d(dist_a) :
        _dist_local_minima_2d(reshape(dist_a,dims_a))
    indices_b=dims_b===nothing ? _dist_local_minima_1d(dist_b) :
        _dist_local_minima_2d(reshape(dist_b,dims_b))
    candidates=Tuple{Float64,NTuple{3,Float64},NTuple{3,Float64}}[
        vals_a[i] for i in indices_a]
    append!(candidates,(vals_b[i] for i in indices_b))
    sort!(candidates;by=v->v[1])
    length(candidates)>_DIST_SEED_CAP &&
        (candidates=candidates[1:_DIST_SEED_CAP])
    best=Inf;best_a=(0.0,0.0,0.0);best_b=(0.0,0.0,0.0)
    for (_,p_a,p_b) in candidates
        sqd,ref_a,ref_b=_dist_alternate(
            m,dim_a,tag_a,dim_b,tag_b,p_a,p_b,lc,caller)
        sqd<best && (best=sqd;best_a=ref_a;best_b=ref_b)
    end
    isfinite(best) || throw(ErrorException(
        "$caller: distance minimization found no candidate pairs"))
    return best,best_a,best_b
end

# Ericson RTCD §5.1.9 — squared distance and closest points of two segments.
function _dist_seg_seg(p1::NTuple{3,Float64},q1::NTuple{3,Float64},
                       p2::NTuple{3,Float64},q2::NTuple{3,Float64})
    d1=_arc_sub(q1,p1);d2=_arc_sub(q2,p2);r=_arc_sub(p1,p2)
    a=_dot3(d1,d1);e=_dot3(d2,d2);f=_dot3(d2,r)
    if a<=eps() && e<=eps()
        return _sqdist(p1,p2),p1,p2
    end
    if a<=eps()
        s=0.0;t=clamp(f/e,0.0,1.0)
    else
        c=_dot3(d1,r)
        if e<=eps()
            t=0.0;s=clamp(-c/a,0.0,1.0)
        else
            b=_dot3(d1,d2);denominator=a*e-b*b
            s=denominator!=0.0 ?
                clamp((b*f-c*e)/denominator,0.0,1.0) : 0.0
            t_numerator=b*s+f
            if t_numerator<0.0
                t=0.0;s=clamp(-c/a,0.0,1.0)
            elseif t_numerator>e
                t=1.0;s=clamp((b-c)/a,0.0,1.0)
            else
                t=t_numerator/e
            end
        end
    end
    c1=_add3(p1,_occ_mul(d1,s));c2=_add3(p2,_occ_mul(d2,t))
    return _sqdist(c1,c2),c1,c2
end

# Endpoints of a segment-parametrized curve (built-in or OCC Line), else
# nothing.
function _dist_curve_segment(m::GeoModel,tag::Int,caller::AbstractString)
    haskey(m.discrete,(1,tag)) && return nothing
    occ=_occ_geometry_checked(m,tag,caller)
    if occ!==nothing
        occ.occ===:line || return nothing
    else
        _curve_type(m,tag)===:line || return nothing
    end
    a,b=m.curves[tag]
    haskey(m.points,a) || throw(ArgumentError(
        "$caller: Curve[$tag] references unknown Point[$a]"))
    haskey(m.points,b) || throw(ArgumentError(
        "$caller: Curve[$tag] references unknown Point[$b]"))
    return m.points[a],m.points[b]
end

# face↔face: interior scans plus every boundary curve ↔ other face.
function _dist_face_face(m::GeoModel,tag_a::Int,tag_b::Int,
                       lc::Float64,caller::AbstractString)
    best,best_a,best_b=_dist_minimize(m,2,tag_a,2,tag_b,lc,caller)
    for curve in sort!(unique!(_model_surface_curves(m,tag_a)))
        sqd,p_c,p_f=_dist_minimize(m,1,curve,2,tag_b,lc,caller)
        sqd<best && (best=sqd;best_a=p_c;best_b=p_f)
    end
    for curve in sort!(unique!(_model_surface_curves(m,tag_b)))
        sqd,p_c,p_f=_dist_minimize(m,1,curve,2,tag_a,lc,caller)
        sqd<best && (best=sqd;best_a=p_f;best_b=p_c)
    end
    return best,best_a,best_b
end

# ── solid containment ────────────────────────────────────────────────────────

# Analytic containment under a compact encoding that — when the volume is
# also shelled — still satisfies its materialized topology (the same
# consistency gate `_model_volume_bounds` enforces). Returns `nothing` when
# the volume has no primitive encoding to test against.
function _dist_encoded_contains(m::GeoModel,tag::Int,q::NTuple{3,Float64},
                                lc::Float64,caller::AbstractString)
    entity="Volume[$tag]"
    encodings=(haskey(m.box_extents,tag),haskey(m.cylinders,tag),
               haskey(m.spheres,tag),haskey(m.cones,tag),
               haskey(m.booleans,tag))
    shelled=!isempty(get(m.volumes,tag,Int[]))
    materialized=any(encodings[1:4]) && shelled
    count(identity,encodings)+Int(shelled && !any(encodings))<=1 ||
        throw(ErrorException(
            "$caller: $entity has multiple native encodings; rebuild the " *
            "model"))
    tol=1e-9*max(1.0,lc)
    if encodings[1]
        if materialized
            x0,y0,z0,dx,dy,dz=m.box_extents[tag]
            scale=max(1.0,abs(x0),abs(y0),abs(z0),abs(dx),abs(dy),abs(dz))
            corners=NTuple{3,Float64}[
                (x0+ix*dx,y0+iy*dy,z0+iz*dz)
                for ix in (0,1),iy in (0,1),iz in (0,1)]
            owned=[m.points[point]
                   for point in _model_volume_owned_points(m,tag)]
            (length(owned)==8 && all(
                corner->any(p->_points_close(p,corner,1e-9*scale),owned),
                corners)) || throw(ErrorException(
                    "$caller: $entity has inconsistent box encoding and " *
                    "shell topology; rebuild the model"))
        end
        x,y,z,dx,dy,dz=m.box_extents[tag]
        return x-tol<=q[1]<=x+dx+tol && y-tol<=q[2]<=y+dy+tol &&
               z-tol<=q[3]<=z+dz+tol
    elseif any(encodings[2:4])
        materialized && !_materialized_curved_consistent(m,tag) &&
            throw(ErrorException(
                "$caller: $entity has inconsistent primitive encoding and " *
                "shell topology; rebuild the model"))
        if encodings[2]
            rec=m.cylinders[tag]
            axis=_b3mul(rec.axis,1.0/rec.height)
            rel=_b3(q,rec.center)
            axial=_b3dot(rel,axis)
            radial=_b3(rel,_b3mul(axis,axial))
            return -tol<=axial<=rec.height+tol &&
                   _b3dot(radial,radial)<=(rec.radius+tol)^2
        elseif encodings[3]
            rec=m.spheres[tag]
            return _b3dot(_b3(q,rec.center),_b3(q,rec.center))<=
                (rec.radius+tol)^2
        end
        rec=m.cones[tag]
        axis=_b3mul(rec.axis,1.0/rec.height)
        rel=_b3(q,rec.center)
        axial=_b3dot(rel,axis)
        axial<-tol && return false
        axial>rec.height+tol && return false
        radial=_b3(rel,_b3mul(axis,axial))
        radius=rec.r1+axial*(rec.r2-rec.r1)/rec.height
        return _b3dot(radial,radial)<=(radius+tol)^2
    end
    return nothing
end

# Ray parameter interval where the ray meets an entity's bounding box.
function _dist_ray_bbox(m::GeoModel,dim::Int,tag::Int,
                        P::NTuple{3,Float64},D::NTuple{3,Float64},
                        caller::AbstractString)
    bounds=_model_entity_bounding_box(m,dim,tag,caller)
    t_lo=-Inf;t_hi=Inf
    for axis in 1:3
        d=D[axis]
        if abs(d)<=eps()
            (bounds[axis]<=P[axis]<=bounds[axis+3]) || return nothing
            continue
        end
        a=(bounds[axis]-P[axis])/d;b=(bounds[axis+3]-P[axis])/d
        lo,hi=minmax(a,b)
        t_lo=max(t_lo,lo);t_hi=min(t_hi,hi)
    end
    (t_hi>=t_lo && t_hi>0.0) || return nothing
    return max(t_lo,0.0),t_hi
end

# Ray∩quadric roots with a tangency band → (ts, ambiguous).
function _dist_ray_quadric(A::Float64,B::Float64,C::Float64)
    epsd=8.0*eps()*max(B*B,abs(4.0*A*C),1.0)
    if abs(A)<=1e-30
        abs(B)<=1e-30 && return Float64[],C==0.0
        return [-C/B],false
    end
    disc=B*B-4.0*A*C
    disc<0.0 && return (Float64[],abs(disc)<=epsd)
    abs(disc)<=epsd && return Float64[],true
    sq=sqrt(disc)
    return [(-B+sq)/(2.0*A),(-B-sq)/(2.0*A)],false
end

# Ray-parameter roots of ray∩analytic-surface (unbounded kind), +ambiguity.
function _dist_ray_surface_params(g,P::NTuple{3,Float64},
                                  D::NTuple{3,Float64})
    if g.occ===:plane
        den=_dot3(D,g.axis)
        if abs(den)<=1e-12
            return Float64[],
                abs(_dot3(_arc_sub(P,g.center),g.axis))<=_BREP_TOL
        end
        return [_dot3(_arc_sub(g.center,P),g.axis)/den],false
    elseif g.occ===:sphere
        rel=_arc_sub(P,g.center)
        B=2.0*_dot3(D,rel);C=_dot3(rel,rel)-g.radius^2
        return _dist_ray_quadric(1.0,B,C)
    elseif g.occ===:cylinder || g.occ===:cone
        rel=_arc_sub(P,g.center)
        da=_dot3(D,g.axis);oa=_dot3(rel,g.axis)
        dp=_arc_sub(D,_occ_mul(g.axis,da))
        rp=_arc_sub(rel,_occ_mul(g.axis,oa))
        if g.occ===:cylinder
            A=_dot3(dp,dp);B=2.0*_dot3(dp,rp)
            C=_dot3(rp,rp)-g.radius^2
        else
            k=(g.r2-g.r1)/g.height
            A=_dot3(dp,dp)-(k*da)^2
            B=2.0*(_dot3(dp,rp)-k*k*oa*da-k*g.r1*da)
            C=_dot3(rp,rp)-(g.r1+k*oa)^2
        end
        if abs(A)<=_BREP_TOL
            abs(B)>_BREP_TOL && return [-C/B],false
            return Float64[],abs(C)<=1e-7*max(1.0,g.occ===:cylinder ?
                g.radius^2 : (g.r1+g.r2)^2)
        end
        return _dist_ray_quadric(A,B,C)
    end
    throw(ArgumentError(
        "_dist_ray_surface_params: unsupported OCC surface kind $(g.occ)"))
end

# Ray∩torus roots by signed-distance bracketing + bisection on the ray∩bbox
# interval; tangential grazes report ambiguous.
function _dist_ray_torus(m::GeoModel,stag::Int,g,
                         P::NTuple{3,Float64},D::NTuple{3,Float64},
                         caller::AbstractString)
    interval=_dist_ray_bbox(m,2,stag,P,D,caller)
    interval===nothing && return Float64[],false
    t_lo,t_hi=interval
    function implicit(t)
        point=_add3(P,_occ_mul(D,t))
        rel=_arc_sub(point,g.center)
        x=_dot3(rel,g.X);y=_dot3(rel,_occ_frame_y(g));z=_dot3(rel,g.axis)
        rho=sqrt(x*x+y*y)
        return (rho-g.r1)^2+z*z-g.r2^2
    end
    roots=Float64[]
    ambiguous=false
    scale=max(1.0,g.r1^2,g.r2^2)
    band=1e-12*scale
    previous_t=t_lo;previous_g=implicit(t_lo)
    abs(previous_g)<=band && (ambiguous=true)
    for i in 1:_DIST_RAY_SAMPLE
        t=fma(t_hi-t_lo,i/_DIST_RAY_SAMPLE,t_lo)
        value=implicit(t)
        abs(value)<=band && (ambiguous=true)
        if (previous_g<0.0)!=(value<0.0) && previous_g!=value
            lo_t,hi_t=previous_t,t
            for _ in 1:60
                mid=(lo_t+hi_t)/2
                mid_g=implicit(mid)
                (mid_g<0.0)==(previous_g<0.0) ? (lo_t=mid) : (hi_t=mid)
            end
            push!(roots,(lo_t+hi_t)/2)
        end
        previous_t=t;previous_g=value
    end
    return roots,ambiguous
end

# Ray∩face crossing count with an ambiguity flag → (count, ambiguous).
function _dist_ray_face(m::GeoModel,stag::Int,P::NTuple{3,Float64},
                        D::NTuple{3,Float64},lc::Float64,
                        caller::AbstractString)
    tol=1e-9*max(1.0,lc)
    record=get(m.discrete,(2,stag),nothing)
    if record!==nothing
        hits=Float64[];ambiguous=false
        for points in _dist_discrete_frames(
            m,record,2,caller,"discrete Surface[$stag]")
            hit=_dist_ray_triangle(P,D,points[1],points[2],points[3],tol)
            hit===nothing && continue
            t,u,v,edge=hit
            edge && (ambiguous=true;continue)
            t>tol && push!(hits,t)
        end
        sort!(hits)
        count=0;previous=-Inf
        for t in hits
            t-previous<=tol && continue
            count+=1;previous=t
        end
        return count,ambiguous
    end
    geometry=get(m.surface_geometry,stag,nothing)
    if geometry!==nothing && hasproperty(geometry,:occ)
        if geometry.occ===:torus
            params,ambiguous=_dist_ray_torus(m,stag,geometry,P,D,caller)
        else
            params,ambiguous=_dist_ray_surface_params(geometry,P,D)
        end
        count=0
        for t in params
            if abs(t)<=tol
                ambiguous=true;continue
            end
            t>tol || continue
            point=_add3(P,_occ_mul(D,t))
            u,v,_=_occ_surface_closest(geometry,point,lc)
            state=_occ_face_state(m,stag,geometry,(u,v),_OCC_CONFUSION)
            (state===:on || state===:unknown) && (ambiguous=true)
            state===:in && (count+=1)
        end
        return count,ambiguous
    end
    kind=_surface_type(m,stag)
    if kind===:plane
        plane=_model_plane_frame(m,stag,caller)
        rhs=plane.properties[4]
        den=_dot3(D,plane.normal)
        if abs(den)<=1e-12
            return 0,abs(_dot3(plane.normal,P)-rhs)<=tol
        end
        t=(rhs-_dot3(plane.normal,P))/den
        if abs(t)<=tol
            return 0,true
        end
        t>tol || return 0,false
        point=_add3(P,_occ_mul(D,t))
        state=_dist_plane_state(m,stag,plane,point,caller)
        state===:in && return 1,false
        state===:on && return 0,true
        return 0,false
    elseif kind in (:ruled,:tric)
        return _dist_ray_ruled(m,stag,P,D,tol,caller)
    end
    throw(ArgumentError(
        "$caller: Surface[$stag] has unsupported kind $kind"))
end

# Möller–Trumbore ray∩triangle → (t,u,v,edge_hit) or nothing.
function _dist_ray_triangle(P::NTuple{3,Float64},D::NTuple{3,Float64},
                            a,b,c,tol::Float64)
    e1=_arc_sub(b,a);e2=_arc_sub(c,a)
    pvec=_model_cross3(D,e2)
    determinant=_dot3(e1,pvec)
    if abs(determinant)<=1e-30
        n=_model_cross3(e1,e2)
        nl=sqrt(_sqlen(n))
        nl<=0.0 && return nothing
        offset=_dot3(_arc_sub(P,a),n)/nl
        abs(offset)<=tol && return (-1.0,0.0,0.0,true)
        return nothing
    end
    inv=1.0/determinant
    tvec=_arc_sub(P,a)
    u=_dot3(tvec,pvec)*inv
    (u<-tol || u>1.0+tol) && return nothing
    qvec=_model_cross3(tvec,e1)
    v=_dot3(D,qvec)*inv
    (v<-tol || u+v>1.0+tol) && return nothing
    t=_dot3(e2,qvec)*inv
    band=1e-9
    edge=u<band || v<band || u+v>1.0-band
    return (t,u,v,edge)
end

# (u,v,t) Newton solve S(u,v) = P + t·D for a ruled patch; converged roots are
# classified by their param-box position.
function _dist_ray_ruled(m::GeoModel,stag::Int,P::NTuple{3,Float64},
                         D::NTuple{3,Float64},tol::Float64,
                         caller::AbstractString)
    interval=_dist_ray_bbox(m,2,stag,P,D,caller)
    interval===nothing && return 0,false
    t_lo,t_hi=interval
    roots=Tuple{Float64,Float64,Float64}[]
    ambiguous=false
    for init_u in 1:7, init_v in 1:7, init_t in 1:3
        u=(init_u-0.5)/7.0;v=(init_v-0.5)/7.0
        t=fma(t_hi-t_lo,init_t/4.0,t_lo)
        converged=false
        for _ in 1:15
            point=_ruled_surface_point(m,stag,u,v,caller)
            du=_ruled_surface_d1(m,stag,u,v,1,caller)
            dv=_ruled_surface_d1(m,stag,u,v,2,caller)
            residual=_arc_sub(point,_add3(P,_occ_mul(D,t)))
            jac=_invert_singular3x3(
                [du[1] dv[1] -D[1];du[2] dv[2] -D[2];du[3] dv[3] -D[3]])
            step=(jac[1,1]*residual[1]+jac[1,2]*residual[2]+
                  jac[1,3]*residual[3],
                  jac[2,1]*residual[1]+jac[2,2]*residual[2]+
                  jac[2,3]*residual[3],
                  jac[3,1]*residual[1]+jac[3,2]*residual[2]+
                  jac[3,3]*residual[3])
            u-=step[1];v-=step[2];t-=step[3]
            if step[1]^2+step[2]^2+step[3]^2<=1e-24*max(1.0,t_lo+t_hi)^2
                converged=true;break
            end
        end
        converged || continue
        point=_ruled_surface_point(m,stag,u,v,caller)
        resid=_occ_modulus(_arc_sub(point,_add3(P,_occ_mul(D,t))))
        resid>100.0*tol && (ambiguous=true;continue)
        (t>tol && -1e-9<=u<=1.0+1e-9 && -1e-9<=v<=1.0+1e-9) ||
            continue
        t_end=t
        duplicate=any(root->abs(root[3]-t_end)<=tol,roots)
        duplicate || push!(roots,(u,v,t_end))
    end
    count=0
    for (u,v,t) in roots
        if u<1e-9 || u>1.0-1e-9 || v<1e-9 || v>1.0-1e-9
            ambiguous=true;continue
        end
        count+=1
    end
    # Tangency probe: a ray nearly tangent at a root grazes the face.
    for (u,v,t) in roots
        normal=_ruled_surface_normal(m,stag,u,v,caller)
        abs(_dot3(normal,D))<1e-6 && (ambiguous=true)
    end
    return count,ambiguous
end

# Ray∩tri-mesh parity count → (count, ambiguous).
function _dist_ray_mesh(mesh::Mesh,P::NTuple{3,Float64},D::NTuple{3,Float64},
                        tol::Float64)
    hits=Float64[];ambiguous=false
    for triangle in axes(mesh.tris,2)
        points=ntuple(k->_model_mesh_coordinate(mesh,mesh.tris[k,triangle]),3)
        hit=_dist_ray_triangle(P,D,points[1],points[2],points[3],tol)
        hit===nothing && continue
        t,u,v,edge=hit
        edge && (ambiguous=true;continue)
        t>tol && push!(hits,t)
    end
    sort!(hits)
    count=0;previous=-Inf
    for t in hits
        t-previous<=tol && continue
        count+=1;previous=t
    end
    return count,ambiguous
end

# Point-in-solid by ray parity over the boundary faces or a snapshot mesh.
function _dist_point_in_shell(m::GeoModel,tag::Int,q::NTuple{3,Float64},
                              lc::Float64,caller::AbstractString)
    tol=1e-9*max(1.0,lc)
    faces=_dist_volume_faces(m,tag,caller)
    if isempty(faces)
        haskey(m.booleans,tag) || throw(ArgumentError(
            "$caller: Volume[$tag] has no boundary to contain against"))
        mesh=_boolean_result_surface(m,tag,caller)
        ntris(mesh)>0 || throw(ErrorException(
            "$caller: Volume[$tag] has no boundary triangles"))
        for D in _BREP_RAY_DIRS
            count,ambiguous=_dist_ray_mesh(mesh,q,D,tol)
            ambiguous && continue
            return isodd(count)
        end
        return nothing
    end
    for D in _BREP_RAY_DIRS
        count=0;ambiguous=false
        for stag in faces
            c,amb=_dist_ray_face(m,stag,q,D,lc,caller)
            count+=c;ambiguous|=amb
        end
        ambiguous && continue
        return isodd(count)
    end
    return nothing
end

# Solid containment: encodings → shell parity → boolean snapshot mesh.
# `nothing` when the geometry cannot decide (all rays degenerate).
function _dist_point_in_volume(m::GeoModel,tag::Int,q::NTuple{3,Float64},
                               lc::Float64,caller::AbstractString)
    encoded=_dist_encoded_contains(m,tag,q,lc,caller)
    encoded!==nothing && return encoded
    return _dist_point_in_shell(m,tag,q,lc,caller)
end

# Deterministic interior boundary sample of a volume: a face param-center of
# a shelled/discrete boundary, else a snapshot triangle centroid.
function _dist_volume_rep(m::GeoModel,tag::Int,which::Int,
                          lc::Float64,caller::AbstractString)
    faces=_dist_volume_faces(m,tag,caller)
    if isempty(faces)
        haskey(m.booleans,tag) || throw(ArgumentError(
            "$caller: Volume[$tag] has no materialized boundary"))
        mesh=_boolean_result_surface(m,tag,caller)
        ntris(mesh)>0 || throw(ErrorException(
            "$caller: Volume[$tag] has no boundary triangles"))
        index=mod(which-1,ntris(mesh))+1
        corners=ntuple(k->_model_mesh_coordinate(mesh,mesh.tris[k,index]),3)
        return ntuple(axis->(corners[1][axis]+corners[2][axis]+
                             corners[3][axis])/3.0,3)
    end
    stag=faces[mod(which-1,length(faces))+1]
    points,_=_dist_scan(m,2,stag,caller)
    isempty(points) && throw(ErrorException(
        "$caller: Surface[$stag] admits no samples; rebuild the model"))
    return points[(length(points)-1)÷2+1]
end

# Deterministic rep point on a curve (mid-parameter).
function _dist_curve_rep(m::GeoModel,tag::Int,which::Int,
                         caller::AbstractString)
    record=get(m.discrete,(1,tag),nothing)
    if record!==nothing
        frames=_dist_discrete_frames(m,record,1,caller,"discrete Curve[$tag]")
        frame=frames[mod(which-1,length(frames))+1]
        return ntuple(axis->(frame[1][axis]+frame[2][axis])/2,3)
    end
    t0,t1=_dist_curve_bounds(m,tag,caller)
    frac=(0.5,0.25,0.75,0.125,0.875)[mod(which-1,5)+1]
    return _model_curve_point(m,tag,fma(t1-t0,frac,t0),caller)
end

# Deterministic rep point on a face (param-center sample).
function _dist_face_rep(m::GeoModel,tag::Int,which::Int,
                        caller::AbstractString)
    record=get(m.discrete,(2,tag),nothing)
    if record!==nothing
        frames=_dist_discrete_frames(m,record,2,caller,
                                     "discrete Surface[$tag]")
        frame=frames[mod(which-1,length(frames))+1]
        return ntuple(axis->(frame[1][axis]+frame[2][axis]+frame[3][axis])/3,3)
    end
    points,_=_dist_scan(m,2,tag,caller)
    isempty(points) && throw(ErrorException(
        "$caller: Surface[$tag] admits no samples; rebuild the model"))
    return points[mod(which-1,length(points))+1]
end

# Boundary distance of an entity to a shelled/mesh-shell volume.
function _dist_volume_boundary_min(m::GeoModel,vtag::Int,dim::Int,tag::Int,
                                   lc::Float64,caller::AbstractString)
    faces=_dist_volume_faces(m,vtag,caller)
    best=Inf;best_a=(0.0,0.0,0.0);best_b=(0.0,0.0,0.0)
    if isempty(faces)
        haskey(m.booleans,vtag) || throw(ArgumentError(
            "$caller: Volume[$vtag] has no boundary to measure against"))
        mesh=_boolean_result_surface(m,vtag,caller)
        best,best_a,best_b=_dist_entity_mesh_min(
            m,dim,tag,mesh,lc,caller)
        return best,best_a,best_b
    end
    for stag in faces
        sqd,p_a,p_b=_dist_minimize(m,dim,tag,2,stag,lc,caller)
        sqd<best && (best=sqd;best_a=p_a;best_b=p_b)
    end
    isfinite(best) || throw(ErrorException(
        "$caller: Volume[$vtag] admits no boundary closest point"))
    return best,best_a,best_b
end

# Point↔volume solid distance: inside → (0, q), else the boundary min.
function _dist_volume_project(m::GeoModel,tag::Int,q::NTuple{3,Float64},
                              lc::Float64,caller::AbstractString)
    inside=_dist_point_in_volume(m,tag,q,lc,caller)
    inside===true && return 0.0,q
    best=Inf;best_point=(0.0,0.0,0.0)
    faces=_dist_volume_faces(m,tag,caller)
    if isempty(faces)
        haskey(m.booleans,tag) || throw(ArgumentError(
            "$caller: Volume[$tag] has no boundary to measure against"))
        mesh=_boolean_result_surface(m,tag,caller)
        for triangle in axes(mesh.tris,2)
            points=ntuple(
                k->_model_mesh_coordinate(mesh,mesh.tris[k,triangle]),3)
            weights,dist=_triangle_closest(
                q,points[1],points[2],points[3])
            if dist<best
                best=dist
                best_point=ntuple(axis->sum(
                    k->weights[k]*points[k][axis],1:3),3)
            end
        end
    else
        for stag in faces
            dist,pt=_dist_face_project(m,stag,q,lc,caller)
            dist<best && (best=dist;best_point=pt)
        end
    end
    isfinite(best) || throw(ErrorException(
        "$caller: Volume[$tag] admits no boundary closest point"))
    tol=1e-9*max(1.0,lc)
    inside===nothing && best>tol*tol && throw(ErrorException(
        "$caller: solid containment query is degenerate on all ray " *
        "directions"))
    return best,best_point
end

# Entity↔tri-mesh boundary min — the unshelled-boolean fallback.
function _dist_entity_mesh_min(m::GeoModel,dim::Int,tag::Int,mesh::Mesh,
                               lc::Float64,caller::AbstractString)
    best=Inf;best_a=(0.0,0.0,0.0);best_b=(0.0,0.0,0.0)
    if dim==0
        q=_dist_point_coord(m,tag,caller)
        for triangle in axes(mesh.tris,2)
            points=ntuple(
                k->_model_mesh_coordinate(mesh,mesh.tris[k,triangle]),3)
            weights,dist=_triangle_closest(
                q,points[1],points[2],points[3])
            if dist<best
                best=dist;best_a=q
                best_b=ntuple(axis->sum(
                    k->weights[k]*points[k][axis],1:3),3)
            end
        end
    else
        samples,_=_dist_scan(m,dim,tag,caller)
        for triangle in axes(mesh.tris,2)
            points=ntuple(
                k->_model_mesh_coordinate(mesh,mesh.tris[k,triangle]),3)
            for sample in samples
                a=sample
                for _ in 1:8
                    _,b=_dist_mesh_tri_project(points,a)
                    _,next_a=_dist_entity_project(
                        m,dim,tag,b,lc,caller)
                    sqd=_sqdist(next_a,b)
                    sqd<best && (best=sqd;best_a=next_a;best_b=b)
                    next_a==a && break
                    a=next_a
                end
            end
        end
    end
    isfinite(best) || throw(ErrorException(
        "$caller: Volume boundary mesh has no triangles"))
    return best,best_a,best_b
end

# Closest point on a single triangle → (sqd, xyz).
function _dist_mesh_tri_project(points,q::NTuple{3,Float64})
    weights,dist=_triangle_closest(q,points[1],points[2],points[3])
    return dist,ntuple(axis->sum(k->weights[k]*points[k][axis],1:3),3)
end

# Closest point on a whole tri-mesh → (sqd, xyz).
function _dist_mesh_project(mesh::Mesh,q::NTuple{3,Float64})
    best=Inf;best_point=(0.0,0.0,0.0)
    for triangle in axes(mesh.tris,2)
        points=ntuple(
            k->_model_mesh_coordinate(mesh,mesh.tris[k,triangle]),3)
        weights,dist=_triangle_closest(q,points[1],points[2],points[3])
        if dist<best
            best=dist
            best_point=ntuple(axis->sum(
                k->weights[k]*points[k][axis],1:3),3)
        end
    end
    return best,best_point
end

# Tri-mesh↔tri-mesh min via alternating triangle projection.
function _dist_mesh_mesh_min(mesh_a::Mesh,mesh_b::Mesh)
    best=Inf;best_a=(0.0,0.0,0.0);best_b=(0.0,0.0,0.0)
    ntris(mesh_a)>0 && ntris(mesh_b)>0 || throw(ErrorException(
        "mesh↔mesh distance requires non-empty meshes"))
    for triangle in axes(mesh_a.tris,2)
        points=ntuple(
            k->_model_mesh_coordinate(mesh_a,mesh_a.tris[k,triangle]),3)
        seeds=(points[1],points[2],points[3],
               ntuple(axis->(points[1][axis]+points[2][axis]+
                             points[3][axis])/3.0,3))
        for seed in seeds
            a=seed
            for _ in 1:16
                _,b=_dist_mesh_project(mesh_b,a)
                _,next_a=_dist_mesh_project(mesh_a,b)
                sqd=_sqdist(next_a,b)
                sqd<best && (best=sqd;best_a=next_a;best_b=b)
                next_a==a && break
                a=next_a
            end
        end
    end
    return best,best_a,best_b
end

# Boundary minimum between two volumes, mixing shells and snapshot meshes.
function _dist_shells_min(m::GeoModel,tag_a::Int,tag_b::Int,
                          lc::Float64,caller::AbstractString)
    faces_a=_dist_volume_faces(m,tag_a,caller)
    faces_b=_dist_volume_faces(m,tag_b,caller)
    isempty(faces_a) && !haskey(m.booleans,tag_a) && throw(ArgumentError(
        "$caller: Volume[$tag_a] has no boundary to measure against"))
    isempty(faces_b) && !haskey(m.booleans,tag_b) && throw(ArgumentError(
        "$caller: Volume[$tag_b] has no boundary to measure against"))
    best=Inf;best_a=(0.0,0.0,0.0);best_b=(0.0,0.0,0.0)
    if isempty(faces_a) && isempty(faces_b)
        mesh_a=_boolean_result_surface(m,tag_a,caller)
        mesh_b=_boolean_result_surface(m,tag_b,caller)
        return _dist_mesh_mesh_min(mesh_a,mesh_b)
    elseif isempty(faces_a)
        mesh_a=_boolean_result_surface(m,tag_a,caller)
        for fb in faces_b
            sqd,p_fb,p_mesh=_dist_entity_mesh_min(
                m,2,fb,mesh_a,lc,caller)
            sqd<best && (best=sqd;best_a=p_mesh;best_b=p_fb)
        end
    elseif isempty(faces_b)
        mesh_b=_boolean_result_surface(m,tag_b,caller)
        for fa in faces_a
            sqd,p_fa,p_mesh=_dist_entity_mesh_min(
                m,2,fa,mesh_b,lc,caller)
            sqd<best && (best=sqd;best_a=p_fa;best_b=p_mesh)
        end
    else
        for fa in faces_a, fb in faces_b
            cand,p_fa,p_fb=_dist_face_face(m,fa,fb,lc,caller)
            cand<best && (best=cand;best_a=p_fa;best_b=p_fb)
        end
    end
    isfinite(best) || throw(ErrorException(
        "$caller: distance minimization found no candidate pairs"))
    return best,best_a,best_b
end

# Pair dispatch on canonical (da <= db) entities → (sqd, pa, pb).
function _dist_pair(m::GeoModel,dim_a::Int,tag_a::Int,
                    dim_b::Int,tag_b::Int,lc::Float64,
                    caller::AbstractString)
    # Contact scale: an entity pair whose boundary residual fits inside the
    # OCC Confusion band counts as touching (distance 0 up to tolerance).
    contact=1e-7*max(1.0,lc)
    if dim_a==0 && dim_b==0
        pa=_dist_point_coord(m,tag_a,caller)
        pb=_dist_point_coord(m,tag_b,caller)
        return _sqdist(pa,pb),pa,pb
    elseif dim_a==0
        q=_dist_point_coord(m,tag_a,caller)
        sqd,pb=_dist_entity_project(m,dim_b,tag_b,q,lc,caller)
        return sqd,q,pb
    elseif dim_a==1 && dim_b==1
        seg_a=_dist_curve_segment(m,tag_a,caller)
        seg_b=_dist_curve_segment(m,tag_b,caller)
        if seg_a!==nothing && seg_b!==nothing
            return _dist_seg_seg(seg_a[1],seg_a[2],seg_b[1],seg_b[2])
        end
        return _dist_minimize(m,1,tag_a,1,tag_b,lc,caller)
    elseif dim_a==1 && dim_b==2
        return _dist_minimize(m,1,tag_a,2,tag_b,lc,caller)
    elseif dim_a==2 && dim_b==2
        return _dist_face_face(m,tag_a,tag_b,lc,caller)
    elseif dim_a<3 && dim_b==3
        # curve/face ↔ volume: boundary minimum plus a containment rep test.
        sqd,p_a,p_b=_dist_volume_boundary_min(
            m,tag_b,dim_a,tag_a,lc,caller)
        sqd<=contact*contact && return sqd,p_a,p_b
        resolved=false;inner_rep=nothing
        for which in 1:4
            rep=dim_a==1 ? _dist_curve_rep(m,tag_a,which,caller) :
                           _dist_face_rep(m,tag_a,which,caller)
            state=_dist_point_in_volume(m,tag_b,rep,lc,caller)
            state===nothing && continue
            resolved=true
            state && (inner_rep=rep;break)
        end
        inner_rep!==nothing && return 0.0,inner_rep,inner_rep
        resolved || throw(ErrorException(
            "$caller: solid containment query is degenerate on all ray " *
            "directions"))
        return sqd,p_a,p_b
    end
    # (3,3): 0 when shells touch, when a rep of A lies in B, or B in A.
    resolved=false;inner_rep=nothing
    for which in 1:4
        rep_a=_dist_volume_rep(m,tag_a,which,lc,caller)
        state=_dist_point_in_volume(m,tag_b,rep_a,lc,caller)
        state===nothing && continue
        resolved=true
        state && (inner_rep=rep_a;break)
    end
    inner_rep!==nothing && return 0.0,inner_rep,inner_rep
    for which in 1:4
        rep_b=_dist_volume_rep(m,tag_b,which,lc,caller)
        state=_dist_point_in_volume(m,tag_a,rep_b,lc,caller)
        state===nothing && continue
        resolved=true
        state && (inner_rep=rep_b;break)
    end
    inner_rep!==nothing && return 0.0,inner_rep,inner_rep
    sqd,p_a,p_b=_dist_shells_min(m,tag_a,tag_b,lc,caller)
    (resolved || sqd<=contact*contact) || throw(ErrorException(
        "$caller: solid containment query is degenerate on all ray " *
        "directions"))
    return sqd,p_a,p_b
end

"""
    model_distance(model, dim1, tag1, dim2, tag2)
        -> (Float64, NTuple{3,Float64}, NTuple{3,Float64})

Minimum Euclidean distance between two entities plus one attaining point on
each — `(distance, point1, point2)` — matching `gmsh.model.getDistance`
semantics for dimensions 0–3. Distance to a Volume is the distance to the
solid itself: it is zero whenever the other entity meets or lies inside the
volume, and `point1 == point2` is then an interior contact point rather than
a boundary projection (OCC `BRepExtrema_DistShapeShape` convention).
When the contact set is not unique the returned pair is deterministic.

Analytic primitives (point/Line/arc/spline curves, Plane and OCC
Plane/Cylinder/Sphere/Cone/Torus faces) resolve through the same projectors
as [`model_closest_point`](@ref); explicit Plane trims use sampled-boundary
classification. Curve↔curve, curve↔face, and face↔face pairs minimize by
deterministic multi-seed scans plus alternating-projection refinement —
line↔line uses the segment closed form. Volume queries combine the
boundary distances with a solid-containment test: the compact primitive
encodings while they satisfy their materialized shells, else ray parity
over the materialized shell (or the operation-time snapshot mesh of a
Boolean result). Queries no fixed ray direction resolves — and geometry
records that cannot supply a boundary — raise rather than approximate.
"""
function model_distance(m::GeoModel,dim1,tag1,dim2,tag2)
    caller="model_distance"
    d1,t1=_model_metadata_entity(m,dim1,tag1,caller)
    d2,t2=_model_metadata_entity(m,dim2,tag2,caller)
    # Canonicalize the pair order by (dimension,tag) so equal-dimension
    # queries also return bitwise-symmetric closest-point pairs.
    swapped=(d1,t1)>(d2,t2)
    da,ta,db,tb=swapped ? (d2,t2,d1,t1) : (d1,t1,d2,t2)
    lc=_model_lc(m)
    sqd,pa,pb=_dist_pair(m,da,ta,db,tb,lc,caller)
    distance=sqrt(sqd)
    (isfinite(distance) && all(isfinite,pa) && all(isfinite,pb)) ||
        throw(ErrorException(
            "$caller: non-finite geometry in the distance query"))
    return swapped ? (distance,pb,pa) : (distance,pa,pb)
end
