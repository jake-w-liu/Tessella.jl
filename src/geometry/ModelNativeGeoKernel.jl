# Built-in GEO owns CAD identities and declared Curve constraints independently
# of the synchronized GModel. API retags and live mesh setters do not rewrite them.
function _native_geo_kernel_admit_point!(m::GeoModel,tag::Int)
    coordinate=m.points[tag]
    m.meshing.native_geo_kernel.points[tag]=(
        coordinate[1],coordinate[2],coordinate[3],get(m.point_size,tag,0.0))
    m.meshing.native_geo_kernel.point_max=max(m.meshing.native_geo_kernel.point_max,tag)
    delete!(m.meshing.native_geo_kernel.point_aliases,tag)
    return nothing
end

function _native_geo_kernel_admit_curve!(m::GeoModel,tag::Int)
    kernel=m.meshing.native_geo_kernel
    endpoints=m.curves[tag]
    controls=get(m.curve_control_points,tag,nothing)
    ids=controls===nothing ? collect(endpoints) : controls
    # A Curve in another kernel has no built-in Vertex records.
    all(point->haskey(kernel.points,get(kernel.point_aliases,point,point)),ids) ||
        return nothing
    geometry=get(m.curve_geometry,tag,nothing)
    geometry!==nothing && haskey(geometry,:occ) && return nothing
    ends=ntuple(axis->get(kernel.point_aliases,endpoints[axis],endpoints[axis]),2)
    bound_controls=controls===nothing ? nothing :
        Int[get(kernel.point_aliases,point,point) for point in controls]
    kernel.curves[tag]=_ModelNativeGeoCurve(
        ends,bound_controls,get(m.curve_types,tag,:line),deepcopy(geometry))
    kernel.curve_max=max(kernel.curve_max,tag)
    delete!(kernel.curve_aliases,tag)
    return nothing
end

function _native_geo_kernel_declare_transfinite!(m::GeoModel,tag::Int,spec)
    kernel=m.meshing.native_geo_kernel
    haskey(kernel.curves,tag) || return nothing
    kernel.transfinite_curves[tag]=spec
    return nothing
end

function _native_geo_kernel_set_transfinite_curve!(
    m::GeoModel,tag,num_nodes,mesh_type="Progression",coef=1.0)
    caller="API.model.geo.mesh.set_transfinite_curve"
    requested=_mesh_attr_positive_int(tag,caller,"tag")
    count=max(2,_mesh_attr_positive_int(num_nodes,caller,"num_nodes"))
    mesh_type isa AbstractString || throw(ArgumentError(
        "$caller: mesh_type must be a string"))
    kind=mesh_type=="Bump" ? :bump : mesh_type=="Beta" ? :beta : :progression
    coefficient=_mesh_attr_finite(coef,caller,"coef")
    spec=(num_nodes=count,kind=kind,coef=abs(coefficient),reversed=coefficient<0)
    kernel=m.meshing.native_geo_kernel
    if requested==0
        for curve in keys(kernel.curves)
            kernel.transfinite_curves[curve]=spec
        end
    else
        curve=abs(requested)
        haskey(kernel.curves,curve) || throw(ArgumentError(
            "$caller: unknown built-in GEO Curve[$curve]"))
        kernel.transfinite_curves[curve]=spec
    end
    return nothing
end

function _native_geo_kernel_retag!(kernel,dimension::Int,old_tag::Int,new_tag::Int)
    records=dimension==0 ? kernel.points : dimension==1 ? kernel.curves :
        dimension==2 ? kernel.plane_surfaces : nothing
    records===nothing && return nothing
    aliases=dimension==0 ? kernel.point_aliases : dimension==1 ? kernel.curve_aliases : kernel.surface_aliases
    origin=get(aliases,old_tag,old_tag)
    haskey(records,origin) || return nothing
    delete!(aliases,old_tag)
    origin==new_tag || (aliases[new_tag]=origin)
    if dimension==2 && haskey(kernel.plane_frames,old_tag)
        kernel.plane_frames[new_tag]=pop!(kernel.plane_frames,old_tag)
    end
    if dimension==1 && haskey(kernel.retired_point_curves,old_tag)
        kernel.retired_point_curves[new_tag]=pop!(kernel.retired_point_curves,old_tag)
    end
    return nothing
end

function _native_geo_kernel_update_point!(m::GeoModel,tag::Int)
    kernel=m.meshing.native_geo_kernel
    origin=get(kernel.point_aliases,tag,tag)
    haskey(kernel.points,origin) || return nothing
    coordinate=m.points[tag]
    kernel.points[origin]=(
        coordinate[1],coordinate[2],coordinate[3],get(m.point_size,tag,0.0))
    return nothing
end

function _native_geo_kernel_synchronize!(m::GeoModel;context_lc=nothing,geometry_tolerance=1.0e-8)
    kernel=m.meshing.native_geo_kernel
    isempty(kernel.points) && isempty(kernel.curves) && isempty(kernel.plane_surfaces) && return nothing
    surface_removals=Tuple{Int,Int}[(2,tag) for (tag,origin) in kernel.surface_aliases
        if tag!=origin && haskey(m.surfaces,tag)]
    isempty(surface_removals) || remove_entities!(m,surface_removals,false)
    # A retagged Curve is absent from GEO's tree, so native synchronization
    # removes it and creates a fresh GEdge under the original Curve tag.
    curve_removals=Tuple{Int,Int}[(1,tag) for (tag,origin) in kernel.curve_aliases
        if tag!=origin && haskey(m.curves,tag)]
    isempty(curve_removals) || remove_entities!(m,curve_removals,false)
    for (tag,record) in kernel.points
        m.points[tag]=(record[1],record[2],record[3])
        m.point_size[tag]=record[4]
        m.next_tag[1]=max(m.next_tag[1],tag)
    end
    # Rebind CAD dependencies before the second removal pass. Existing Curve
    # discretizations remain attached when their GEdge identity survives.
    for (tag,record) in kernel.curves
        m.curves[tag]=record.endpoints
        m.next_tag[2]=max(m.next_tag[2],tag)
        if record.controls===nothing
            delete!(m.curve_control_points,tag)
        else
            m.curve_control_points[tag]=copy(record.controls)
        end
        record.kind===:line ? delete!(m.curve_types,tag) :
            (m.curve_types[tag]=record.kind)
        record.geometry===nothing ? delete!(m.curve_geometry,tag) :
            (m.curve_geometry[tag]=deepcopy(record.geometry))
    end
    for tag in sort!(collect(keys(kernel.curve_loops)))
        definition=kernel.curve_loops[tag]
        if haskey(m.loops,tag)
            m.loops[tag]==definition || throw(ArgumentError(
                "native GEO synchronize: conflicting Curve Loop[$tag] declaration"))
        else
            add_curve_loop!(m,definition;tag)
        end
    end
    for tag in sort!(collect(keys(kernel.plane_surfaces)))
        definition=kernel.plane_surfaces[tag]
        if haskey(m.surfaces,tag)
            m.surfaces[tag]==definition || throw(ArgumentError(
                "native GEO synchronize: conflicting Plane Surface[$tag] declaration"))
        else
            add_plane_surface!(m,definition;tag)
        end
        # Every explicit GEO synchronization resets its native CAD frame.
        # This happens before the session bounding-box update, including
        # repeat synchronization of an already published Plane.
        lc=context_lc===nothing ? _model_mesh_lc(m,Float64(geometry_tolerance),"native GEO synchronize") : Float64(context_lc)
        frame=_gmsh_native_plane_birth_frame(m,tag,Float64(geometry_tolerance),lc)
        bounds=_gmsh_native_plane_boundary_bounds(m,tag,frame)
        kernel.plane_frames[tag]=(origin=frame.origin,first=frame.first,
            second=frame.second,normal=frame.normal,bounds=bounds)
    end
    point_removals=Tuple{Int,Int}[(0,tag) for (tag,origin) in kernel.point_aliases
        if tag!=origin && haskey(m.points,tag)]
    # Capture exact retained MLine incidence before the second removal pass
    # destroys the retagged Point's MVertex. Public Curve records survive sync;
    # only a later detached remeshing transaction can discard those old cells.
    retired=Tuple{Int,Dict{Int,Set{UInt64}}}[]
    for (_,point) in point_removals
        vertex=get(m.discrete,(0,point),get(m.meshing.attached,(0,point),nothing))
        vertex===nothing && continue
        labels=Set(vertex.node_tags)
        affected=Dict{Int,Set{UInt64}}()
        for curve in keys(kernel.curves)
            record=get(m.discrete,(1,curve),get(m.meshing.attached,(1,curve),nothing))
            record===nothing && continue
            bound=Set{UInt64}()
            for nodes in record.element_nodes
                for label in nodes
                    label in labels && push!(bound,UInt64(label))
                end
            end
            isempty(bound) || (affected[curve]=bound)
        end
        push!(retired,(point,affected))
    end
    isempty(point_removals) || remove_entities!(m,point_removals,false)
    for (point,curves) in retired
        haskey(m.points,point) && continue
        for (curve,labels) in curves
            union!(get!(()->Set{UInt64}(),kernel.retired_point_curves,curve),labels)
        end
    end
    empty!(kernel.point_aliases)
    empty!(kernel.curve_aliases)
    empty!(kernel.surface_aliases)
    # resetMeshAttributes copies the internal declaration even on a pure sync;
    # API-only transfinite laws have no corresponding internal declaration.
    for tag in keys(kernel.curves)
        spec=get(kernel.transfinite_curves,tag,nothing)
        if spec===nothing
            delete!(m.meshing.transfinite_curves,tag)
        else
            m.meshing.transfinite_curves[tag]=spec
        end
    end
    return nothing
end

# Native geo.addPoint/addLine buffer internal records until synchronization.
# Keep the existing model.add_* convenience factories separate.
function _native_geo_kernel_add_point!(m::GeoModel,x,y,z,mesh_size,tag)
    caller="API.model.geo.add_point"
    coordinate=_finite3(x,y,z,caller)
    size=_mesh_attr_finite(mesh_size,caller,"mesh_size")
    requested=_mesh_attr_positive_int(tag,caller,"tag")
    kernel=m.meshing.native_geo_kernel
    assigned=requested<0 ? kernel.point_max+1 : requested
    haskey(kernel.points,assigned) && throw(ArgumentError(
        "$caller: built-in GEO Point[$assigned] already exists"))
    kernel.points[assigned]=(coordinate[1],coordinate[2],coordinate[3],size)
    kernel.point_max=max(kernel.point_max,assigned)
    return assigned
end

function _native_geo_kernel_add_line!(m::GeoModel,first,last,tag)
    caller="API.model.geo.add_line"
    first_tag=_mesh_attr_positive_int(first,caller,"start tag")
    last_tag=_mesh_attr_positive_int(last,caller,"end tag")
    requested=_mesh_attr_positive_int(tag,caller,"tag")
    kernel=m.meshing.native_geo_kernel
    assigned=requested<0 ? kernel.curve_max+1 : requested
    haskey(kernel.curves,assigned) && throw(ArgumentError(
        "$caller: built-in GEO Curve[$assigned] already exists"))
    haskey(kernel.points,first_tag) && haskey(kernel.points,last_tag) ||
        throw(ArgumentError("$caller: unknown built-in GEO endpoint"))
    kernel.curves[assigned]=_ModelNativeGeoCurve((first_tag,last_tag),nothing,:line,nothing)
    kernel.curve_max=max(kernel.curve_max,assigned)
    return assigned
end

function _native_geo_kernel_add_curve_loop!(m::GeoModel,curves,tag)
    caller="API.model.geo.add_curve_loop"
    requested=_mesh_attr_positive_int(tag,caller,"tag")
    (curves isa AbstractVector || curves isa Tuple) && !isempty(curves) ||
        throw(ArgumentError("$caller: curves must be a nonempty vector or tuple"))
    definition=Int[_mesh_attr_positive_int(value,caller,"curve tag") for value in curves]
    kernel=m.meshing.native_geo_kernel
    all(value->value!=0 && haskey(kernel.curves,abs(value)),definition) ||
        throw(ArgumentError("$caller: unknown built-in GEO Curve"))
    endpoints=[value>0 ? kernel.curves[value].endpoints : reverse(kernel.curves[-value].endpoints)
        for value in definition]
    all(index->endpoints[index][2]==endpoints[mod1(index+1,length(endpoints))][1],eachindex(endpoints)) ||
        throw(ArgumentError("$caller: Curve Loop is not an ordered closed wire"))
    assigned=requested<0 ? kernel.loop_max+1 : requested
    haskey(kernel.curve_loops,assigned) && throw(ArgumentError("$caller: built-in GEO Curve Loop[$assigned] already exists"))
    kernel.curve_loops[assigned]=definition;kernel.loop_max=max(kernel.loop_max,assigned)
    return assigned
end

function _native_geo_kernel_add_plane_surface!(m::GeoModel,loops,tag)
    caller="API.model.geo.add_plane_surface"
    requested=_mesh_attr_positive_int(tag,caller,"tag")
    (loops isa AbstractVector || loops isa Tuple) && !isempty(loops) ||
        throw(ArgumentError("$caller: loops must be a nonempty vector or tuple"))
    definition=Int[_mesh_attr_positive_int(value,caller,"loop tag") for value in loops]
    kernel=m.meshing.native_geo_kernel
    all(value->value>0 && haskey(kernel.curve_loops,value),definition) ||
        throw(ArgumentError("$caller: unknown built-in GEO Curve Loop"))
    assigned=requested<0 ? kernel.surface_max+1 : requested
    haskey(kernel.plane_surfaces,assigned) && throw(ArgumentError("$caller: built-in GEO Plane Surface[$assigned] already exists"))
    kernel.plane_surfaces[assigned]=definition;kernel.surface_max=max(kernel.surface_max,assigned)
    return assigned
end

function _native_geo_kernel_forget!(m::GeoModel,removed)
    kernel=m.meshing.native_geo_kernel
    for (dimension,tag) in removed
        if dimension==0
            delete!(kernel.points,tag)
            delete!(kernel.point_aliases,tag)
            tag==kernel.point_max && (kernel.point_max-=1)
        elseif dimension==1
            delete!(kernel.curves,tag)
            delete!(kernel.transfinite_curves,tag)
            delete!(kernel.retired_point_curves,tag)
            delete!(kernel.curve_aliases,tag)
            tag==kernel.curve_max && (kernel.curve_max-=1)
        elseif dimension==2
            delete!(kernel.plane_surfaces,tag)
            delete!(kernel.plane_frames,tag)
            delete!(kernel.surface_aliases,tag)
            tag==kernel.surface_max && (kernel.surface_max-=1)
        end
    end
    return nothing
end
