# A completed NoNew plan lives only within the mesh operation that requested
# it. In particular, no model/session cache survives edits to Layers, source
# discretization or geometry.
struct _ExtrudeNoNewCompletePlan
    volume::Union{Mesh,MixedMesh}
    surfaces::Dict{Int,Union{Mesh,MixedMesh}}
    columns::Matrix{NTuple{3,Float64}}
    catalog::Union{_ExtrudeNoNewCatalog,_ExtrudeNoNewTwoTriCatalog,_ExtrudeNoNewQuadPatchCatalog,_ExtrudeNoNewQuadStripCatalog}
end

struct _ExtrudeNoNewScope
    model::GeoModel
    volumes::Dict{Int,_ExtrudeNoNewCompletePlan}
    surfaces::Dict{Int,Union{Mesh,MixedMesh}}
end

@inline _extrude_is_nonew(params)=params!==nothing &&
    !isempty(params.layers) && params.recombine &&
    params.quad_to_tri===:no_new_verts

function _extrude_nonew_scope_contract(m::GeoModel,t::Int,
                                       caller::AbstractString)
    link=get(m.meshing.extrude_sources,(3,t),nothing)
    link===nothing && return nothing # The planner diagnoses missing provenance.
    source=abs(link[2])
    surfaces=Set(abs.(_model_volume_boundary_surfaces(m,t,caller)))
    curves=Set{Int}()
    points=Set{Int}()
    for surface in surfaces
        surface==source || !haskey(m.meshing.transfinite_surfaces,surface) ||
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts transfinite cap/lateral overrides " *
                "require the constrained-boundary planner"))
        get(m.meshing.reverse,(2,surface),false) && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts reversed boundary surfaces " *
            "require the oriented-boundary planner"))
        for loop in m.surfaces[surface],signed in m.loops[loop]
            curve=abs(signed)
            push!(curves,curve)
            a,b=m.curves[curve]
            push!(points,a);push!(points,b)
        end
    end
    get(m.meshing.reverse,(3,t),false) && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts reversed volumes " *
        "require the oriented-boundary planner"))
    for curve in curves
        get(m.meshing.reverse,(1,curve),false) && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts reversed boundary curves " *
            "require the oriented-boundary planner"))
    end
    for constraint in values(m.periodic)
        entities=constraint.dim==2 ? surfaces :
                 constraint.dim==1 ? curves :
                 constraint.dim==0 ? points : nothing
        entities===nothing && continue
        (abs(constraint.slave_entity) in entities ||
         abs(constraint.master_entity) in entities) && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts periodic boundary constraints " *
            "require the neighboring-region planner"))
    end
    for (dim,members) in m.meshing.compounds
        entities=dim==2 ? surfaces : dim==1 ? curves : nothing
        entities===nothing && continue
        for tag in members
            tag>0 && tag in entities && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts compound boundary entities " *
                "require the neighboring-region planner"))
        end
    end
    return nothing
end

@inline function _extrude_nonew_coordinate(mesh,node)
    return ntuple(k->_model_projection_coordinate_key(mesh.coords[k,node]),3)
end

@inline function _extrude_nonew_face_key(face::NTuple{3,Int32})
    a,b,c=_model_projection_face_key(face)
    return (a,b,c,Int32(0))
end

@inline _extrude_nonew_face_key(face::NTuple{4,Int32})=
    _model_projection_face_key(face)

struct _ExtrudeNoNewFaceIncidence
    count::UInt8
    orientation::UInt8
end

@inline function _extrude_nonew_face_cycle(face::NTuple{N,Int32}) where N
    first=1
    for index in 2:N
        face[index]<face[first] && (first=index)
    end
    origin=first
    return ntuple(k->k<=N ? face[mod1(origin+k-1,N)] : Int32(0),4)
end

@inline function _extrude_nonew_face_orientation(cycle::NTuple{4,Int32},width::Int)
    # Canonical rotation starts at the minimum vertex. For a quadrangle the
    # middle vertex's rank identifies the opposite pair, and one bit records
    # the order of the two remaining vertices. Reversal flips only that bit.
    width==3 && return UInt8(cycle[2]<cycle[3])
    middle_rank=Int(cycle[2]<cycle[3])+Int(cycle[4]<cycle[3])
    return UInt8(middle_rank<<1) | UInt8(cycle[2]<cycle[4])
end

function _extrude_nonew_count_faces!(
        counts::Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence},
        nodes,msh::Int,caller::AbstractString)
    # The shared topology catalog records vertex incidence, and does not orient
    # all faces outwards. This bounded check uses outward cycles for positive
    # Tet4, Hex8, Pri6 and Pyr5 maps instead.
    faces=msh==4 ? ((1,3,2),(1,2,4),(2,3,4),(3,1,4)) :
          msh==5 ? ((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)) :
          msh==6 ? ((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)) :
          msh==7 ? ((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)) : nothing
    faces===nothing && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts oriented face check does not support volume type $msh"))
    for column in axes(nodes,2),face in faces
        vertices=ntuple(k->nodes[face[k],column],length(face))
        key=_extrude_nonew_face_key(vertices)
        orientation=_extrude_nonew_face_orientation(
            _extrude_nonew_face_cycle(vertices),length(face))
        previous=get(counts,key,nothing)
        if previous===nothing
            counts[key]=_ExtrudeNoNewFaceIncidence(0x01,orientation)
        else
            previous.count==1 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts volume face $key has incidence $(Int(previous.count)+1)"))
            orientation==xor(previous.orientation,0x01) ||
                throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts internal volume face $key has equal or inconsistent orientation"))
            counts[key]=_ExtrudeNoNewFaceIncidence(0x02,previous.orientation)
        end
    end
    return nothing
end

@inline _extrude_nonew_face_count(count::Int)=count
@inline _extrude_nonew_face_count(incidence::_ExtrudeNoNewFaceIncidence)=incidence.count

function _extrude_nonew_count_faces!(counts,nodes,msh::Int,
                                     caller::AbstractString)
    faces=get(_VOLUME_CELL_FACES,msh,nothing)
    faces===nothing && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts emitted unsupported volume type $msh"))
    for column in axes(nodes,2),face in faces
        vertices=ntuple(k->nodes[face[k],column],length(face))
        key=_extrude_nonew_face_key(vertices)
        count=get(counts,key,0)+1
        count<=2 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts volume face $key has incidence $count"))
        counts[key]=count
    end
    return nothing
end

function _extrude_nonew_surface_faces!(actual,part,lookup,tag::Int,
                                       caller::AbstractString)
    remap=Vector{Int32}(undef,nnodes(part))
    for node in eachindex(remap)
        target=get(lookup,_extrude_nonew_coordinate(part,node),Int32(0))
        target>0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts Surface[$tag] node $node " *
            "is absent from the completed volume"))
        remap[node]=target
    end
    if part isa Mesh
        for column in axes(part.tris,2)
            face=ntuple(k->remap[part.tris[k,column]],3)
            key=_extrude_nonew_face_key(face)
            key in actual && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts boundary face $key is repeated"))
            push!(actual,key)
        end
    else
        for block in part.blocks
            block isa ElementBlock || continue
            msh_dimension(block.msh)==2 || continue
            block.msh in (2,3) || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts Surface[$tag] has " *
                "unsupported type $(block.msh)"))
            width=size(block.nodes,1)
            for column in axes(block.nodes,2)
                face=ntuple(k->remap[block.nodes[k,column]],width)
                key=_extrude_nonew_face_key(face)
                key in actual && throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts boundary face $key is repeated"))
                push!(actual,key)
            end
        end
    end
    return nothing
end

# Compare typed faces, so a boundary quadrangle and two triangles are never
# treated as equivalent. Coordinate lookup is confined to this one certified
# isolated region; independent coincident entities never share this registry.
function _extrude_nonew_certify_boundary(volume,surfaces,
        caller::AbstractString;oriented_internal::Bool=false,face_capacity::Int=0)
    counts=oriented_internal ?
        Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence}() :
        Dict{NTuple{4,Int32},Int}()
    face_capacity>0 && sizehint!(counts,face_capacity)
    return _extrude_nonew_certify_boundary(volume,surfaces,counts,caller)
end

function _extrude_nonew_certify_boundary(volume,surfaces,counts,
                                         caller::AbstractString)
    if volume isa Mesh
        _extrude_nonew_count_faces!(counts,volume.tets,4,caller)
    else
        for block in volume.blocks
            block isa ElementBlock || continue
            msh_dimension(block.msh)==3 || continue
            _extrude_nonew_count_faces!(counts,block.nodes,block.msh,caller)
        end
    end
    isempty(counts) && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts completed volume has no cells"))
    boundary=Set{NTuple{4,Int32}}()
    sizehint!(boundary,count(value->_extrude_nonew_face_count(value)==1,values(counts)))
    for (key,incidence) in counts
        _extrude_nonew_face_count(incidence)==1 && push!(boundary,key)
    end
    lookup=Dict{NTuple{3,Float64},Int32}()
    sizehint!(lookup,nnodes(volume))
    for node in 1:nnodes(volume)
        coordinate=_extrude_nonew_coordinate(volume,node)
        haskey(lookup,coordinate) && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts completed volume has " *
            "distinct nodes at one coordinate"))
        lookup[coordinate]=Int32(node)
    end
    actual=Set{NTuple{4,Int32}}()
    sizehint!(actual,length(boundary))
    for tag in sort!(collect(keys(surfaces)))
        _extrude_nonew_surface_faces!(actual,surfaces[tag],lookup,tag,caller)
    end
    actual==boundary || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts finalized surfaces do not exactly " *
        "cover the typed volume boundary"))
    return boundary
end

function _extrude_nonew_complete_plan(m::GeoModel,t::Int,
                                      caller::AbstractString;
                                      min_angle_deg::Real=25.0,
                                      max_periodic_passes=8,
                                      size_field::Union{Nothing,AbstractSizeField}=nothing,
                                      _working_model::Bool=false)
    working=_working_model ? m : deepcopy(m)
    _extrude_nonew_scope_contract(working,t,caller)
    plan=_extrude_nonew_plan(working,t,caller;min_angle_deg=min_angle_deg,
        max_periodic_passes=max_periodic_passes,size_field=size_field)
    volume=_extrude_volume_part(working,plan.sweep,caller)
    surfaces=Dict{Int,Union{Mesh,MixedMesh}}(plan.source_tag=>plan.source_mesh)
    quad_patch=plan.catalog isa _ExtrudeNoNewQuadPatchCatalog
    quad_strip=plan.catalog isa _ExtrudeNoNewQuadStripCatalog
    params,spec,link=_extrude_entity_params(working,2,plan.top_tag,caller)
    surfaces[plan.top_tag]=quad_patch ?
        _extrude_nonew_quad_patch_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : quad_strip ?
        _extrude_nonew_quad_strip_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : _extrude_top_mesh(
        working,plan.top_tag,params,spec,link[2],caller;
        min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
        size_field=size_field,source_mesh=plan.source_mesh)
    two_tri=plan.catalog isa _ExtrudeNoNewTwoTriCatalog
    for tag in plan.lateral_tags
        params,spec,link=_extrude_entity_params(working,2,tag,caller)
        surfaces[tag]=quad_patch ?
            _extrude_nonew_quad_patch_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) : quad_strip ?
            _extrude_nonew_quad_strip_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) : two_tri ?
            _extrude_nonew_two_tri_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) :
            _extrude_lateral_mesh(
                working,tag,params,spec,link[2],caller;edges=plan.edges)
    end
    face_capacity=quad_patch ? _extrude_nonew_quad_patch_face_capacity(plan.catalog) : quad_strip ?
        _extrude_nonew_quad_strip_face_capacity(plan.catalog) : two_tri ?
        (volume isa Mesh ? 16 : 7)*length(plan.catalog.layer_refs)+2 : 0
    _extrude_nonew_certify_boundary(volume,surfaces,caller;
        oriented_internal=two_tri || quad_patch || quad_strip,face_capacity=face_capacity)
    _working_model || (m.curve_params=working.curve_params)
    return _ExtrudeNoNewCompletePlan(volume,surfaces,plan.sweep.cols,plan.catalog)
end

function _extrude_nonew_scope(m::GeoModel,caller::AbstractString;
                               tags=keys(m.volumes),
                               min_angle_deg::Real=25.0,
                               max_periodic_passes=8,
                               size_field::Union{Nothing,AbstractSizeField}=nothing,
                               _working_model::Bool=false)
    selected=Int[]
    for tag in tags
        _extrude_is_nonew(_extrude_gate(m,3,tag)) && push!(selected,tag)
    end
    isempty(selected) && return nothing
    sort!(selected)
    working=_working_model ? m : deepcopy(m)
    volumes=Dict{Int,_ExtrudeNoNewCompletePlan}()
    surfaces=Dict{Int,Union{Mesh,MixedMesh}}()
    for tag in selected
        completed=_extrude_nonew_complete_plan(working,tag,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field,_working_model=true)
        volumes[tag]=completed
        for (surface,part) in completed.surfaces
            haskey(surfaces,surface) && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared Surface[$surface] " *
                "requires the shared-region planner"))
            surfaces[surface]=part
        end
    end
    m.curve_params=working.curve_params
    return _ExtrudeNoNewScope(m,volumes,surfaces)
end

@inline function _extrude_nonew_require_scope(scope::_ExtrudeNoNewScope,
                                              m::GeoModel,caller::AbstractString)
    scope.model===m || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts mesh scope belongs to another model"))
    return nothing
end

function _extrude_nonew_surface_mesh(m::GeoModel,t::Int,
                                      caller::AbstractString;
                                      scope=nothing,
                                      min_angle_deg::Real=25.0,
                                      max_periodic_passes=8,
                                      size_field::Union{Nothing,AbstractSizeField}=nothing)
    if scope!==nothing
        _extrude_nonew_require_scope(scope,m,caller)
        part=get(scope.surfaces,t,nothing)
        part===nothing || return part
    end
    # Do not plan a volume while meshing its plain source. The builder obtains
    # that source through the ordinary surface path before a scope exists.
    _extrude_is_nonew(_extrude_gate(m,2,t)) || return nothing
    owner=0
    for region in _extrude_quadtri_regions(m,t,caller)
        _extrude_is_nonew(_extrude_gate(m,3,region)) || continue
        owner==0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared Surface[$t] " *
            "requires the shared-region planner"))
        owner=region
    end
    owner==0 && return nothing
    if scope===nothing
        completed=_extrude_nonew_complete_plan(m,owner,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field)
        return completed.surfaces[t]
    end
    throw(ArgumentError(
        "$caller: QuadTriNoNewVerts Surface[$t] is absent from this mesh scope"))
end

function _extrude_nonew_volume_mesh(m::GeoModel,t::Int,
                                     caller::AbstractString;
                                     scope=nothing,
                                     min_angle_deg::Real=25.0,
                                     max_periodic_passes=8,
                                     size_field::Union{Nothing,AbstractSizeField}=nothing)
    _extrude_is_nonew(_extrude_gate(m,3,t)) || return nothing
    if scope===nothing
        return _extrude_nonew_complete_plan(m,t,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field).volume
    end
    _extrude_nonew_require_scope(scope,m,caller)
    completed=get(scope.volumes,t,nothing)
    completed===nothing && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts Volume[$t] is absent from this mesh scope"))
    return completed.volume
end
