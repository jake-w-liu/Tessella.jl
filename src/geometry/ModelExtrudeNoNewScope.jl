# A completed NoNew plan lives only within the mesh operation that requested
# it. In particular, no model/session cache survives edits to Layers, source
# discretization or geometry.
struct _ExtrudeNoNewCompletePlan
    volume::Union{Mesh,MixedMesh}
    surfaces::Dict{Int,Union{Mesh,MixedMesh}}
    columns::Matrix{NTuple{3,Float64}}
    catalog::Union{_ExtrudeNoNewCatalog,_ExtrudeNoNewTwoTriCatalog,_ExtrudeNoNewQuadPatchCatalog,_ExtrudeNoNewQuadStripCatalog,_ExtrudeNoNewThreeQuadStripCatalog,_ExtrudeNoNewFourQuadStripCatalog,_ExtrudeNoNewRectGridCatalog,_ExtrudeNoNewB4StripCatalog,_ExtrudeNoNewB4Free.Catalog}
    column_rows::Matrix{Int32}
    # The generic single-quad catalog retains its 6×interval face-state
    # picks; shared-lateral partners read their bound wall diagonal from it.
    faces::Matrix{UInt8}
end

_ExtrudeNoNewCompletePlan(volume,surfaces,columns,catalog)=
    _ExtrudeNoNewCompletePlan(volume,surfaces,columns,catalog,
        Matrix{Int32}(undef,0,0),Matrix{UInt8}(undef,0,0))
_ExtrudeNoNewCompletePlan(volume,surfaces,columns,catalog,column_rows)=
    _ExtrudeNoNewCompletePlan(volume,surfaces,columns,catalog,column_rows,
        Matrix{UInt8}(undef,0,0))

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
    surfaces=Set(abs.(_extrude_nonew_boundary_surfaces(m,t,caller)))
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

@inline function _extrude_nonew_count_face!(
        counts::Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence},
        nodes,column,face::NTuple{N,Int},caller::AbstractString) where N
    vertices=ntuple(k->nodes[face[k],column],Val(N))
    key=_extrude_nonew_face_key(vertices)
    orientation=_extrude_nonew_face_orientation(
        _extrude_nonew_face_cycle(vertices),N)
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
    return nothing
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
        _extrude_nonew_count_face!(counts,nodes,column,face,caller)
    end
    return nothing
end

@inline _extrude_nonew_face_count(count::Int)=count
@inline _extrude_nonew_face_count(incidence::_ExtrudeNoNewFaceIncidence)=incidence.count

@inline function _extrude_nonew_count_face!(counts,nodes,column,
        face::NTuple{N,Int},caller::AbstractString) where N
    vertices=ntuple(k->nodes[face[k],column],Val(N))
    key=_extrude_nonew_face_key(vertices)
    count=get(counts,key,0)+1
    count<=2 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts volume face $key has incidence $count"))
    counts[key]=count
    return nothing
end

function _extrude_nonew_count_faces!(counts,nodes,msh::Int,
                                     caller::AbstractString)
    faces=get(_VOLUME_CELL_FACES,msh,nothing)
    faces===nothing && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts emitted unsupported volume type $msh"))
    for column in axes(nodes,2),face in faces
        _extrude_nonew_count_face!(counts,nodes,column,face,caller)
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
                                      _working_model::Bool=false,_mesh_input=nothing,
                                      _trace_columns::Bool=false,
                                      _shared_diagonals=nothing)
    working=_working_model ? m : deepcopy(m)
    _extrude_nonew_scope_contract(working,t,caller)
    plan=_extrude_nonew_plan(working,t,caller;min_angle_deg=min_angle_deg,
        max_periodic_passes=max_periodic_passes,size_field=size_field,_mesh_input=_mesh_input,
        _shared_diagonals=_shared_diagonals)
    # The insertion renderer and indexed B4 constructors bind their actual
    # column rows explicitly. Other indexed families retain their existing
    # dispatch until they provide the corresponding constructor provenance.
    traced=_trace_columns && plan.sweep isa _ExtrudeVolumeSweep
    emitted=traced ? _extrude_volume_part(working,plan.sweep,caller;
                                          _return_column_rows=true) :
                    _extrude_volume_part(working,plan.sweep,caller)
    volume=traced ? emitted.mesh : emitted
    column_rows=traced ? emitted.column_rows : Matrix{Int32}(undef,0,0)
    if _trace_columns && plan.sweep isa _ExtrudeNoNewIndexedSweep &&
            plan.catalog isa Union{_ExtrudeNoNewB4StripCatalog,_ExtrudeNoNewB4Free.Catalog}
        column_rows=_extrude_nonew_b4_column_rows(plan.sweep,plan.catalog,caller)
    end
    surfaces=Dict{Int,Union{Mesh,MixedMesh}}(plan.source_tag=>plan.source_mesh)
    quad_patch=plan.catalog isa _ExtrudeNoNewQuadPatchCatalog
    quad_strip=plan.catalog isa _ExtrudeNoNewQuadStripCatalog
    three_quad_strip=plan.catalog isa _ExtrudeNoNewThreeQuadStripCatalog
    four_quad_strip=plan.catalog isa _ExtrudeNoNewFourQuadStripCatalog
    dynamic_grid=plan.catalog isa _ExtrudeNoNewDynamicGridCatalog
    retained_b4=_extrude_nonew_retained_b4_catalog(plan.catalog)
    params,spec,link=_extrude_entity_params(working,2,plan.top_tag,caller)
    surfaces[plan.top_tag]=quad_patch ?
        _extrude_nonew_quad_patch_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : quad_strip ?
        _extrude_nonew_quad_strip_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : three_quad_strip ?
        _extrude_nonew_three_quad_strip_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : four_quad_strip ?
        _extrude_nonew_four_quad_strip_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : dynamic_grid ?
        _extrude_nonew_rect_grid_top(working,plan.top_tag,params,
            plan.sweep.cols,plan.catalog,caller) : _extrude_top_mesh(
        working,plan.top_tag,params,spec,link[2],caller;
        min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
        size_field=size_field,source_mesh=plan.source_mesh)
    if _mesh_input!==nothing
        # Every completed family consumes the actual copied-cap vector. The
        # dynamic B4 planner already consumed its paired Quad diagonals; other
        # admitted families verify their actual cap against the typed shell
        # below rather than silently replacing caller topology with a new cap.
        actual_cap=_model_mesh_product(_mesh_input,2,plan.top_tag;caller)
        if actual_cap!==nothing && any(block->!isempty(block.nodes),actual_cap.blocks)
            surfaces[plan.top_tag]=actual_cap
        end
    end
    two_tri=plan.catalog isa _ExtrudeNoNewTwoTriCatalog
    for tag in plan.lateral_tags
        params,spec,link=_extrude_entity_params(working,2,tag,caller)
        surfaces[tag]=quad_patch ?
            _extrude_nonew_quad_patch_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) : quad_strip ?
            _extrude_nonew_quad_strip_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) : three_quad_strip ?
            _extrude_nonew_three_quad_strip_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) : four_quad_strip ?
            _extrude_nonew_four_quad_strip_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) : dynamic_grid ?
            _extrude_nonew_rect_grid_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,retained_b4 ? empty(plan.edges) : plan.edges,caller) : two_tri ?
            _extrude_nonew_two_tri_lateral(working,tag,params,link[2],
                plan.sweep.cols,plan.catalog,plan.edges,caller) :
            _extrude_lateral_mesh(
                working,tag,params,spec,link[2],caller;edges=plan.edges)
    end
    face_capacity=quad_patch ? _extrude_nonew_quad_patch_face_capacity(plan.catalog) : quad_strip ?
        _extrude_nonew_quad_strip_face_capacity(plan.catalog) : three_quad_strip ?
        _extrude_nonew_three_quad_strip_face_capacity(plan.catalog) : four_quad_strip ?
        _extrude_nonew_four_quad_strip_face_capacity(plan.catalog) : dynamic_grid ?
        _extrude_nonew_dynamic_grid_face_capacity(plan.catalog) : two_tri ?
        (volume isa Mesh ? 16 : 7)*length(plan.catalog.layer_refs)+2 : 0
    if retained_b4
        _extrude_nonew_retained_b4_certificate(working,t,volume,surfaces,plan.catalog,caller;
            face_capacity)
    else
        _extrude_nonew_certify_boundary(volume,surfaces,caller;
            oriented_internal=two_tri || quad_patch || quad_strip || three_quad_strip || four_quad_strip || dynamic_grid,
            face_capacity=face_capacity)
    end
    _working_model || (m.curve_params=working.curve_params)
    return _ExtrudeNoNewCompletePlan(volume,surfaces,plan.sweep.cols,plan.catalog,column_rows,
        hasproperty(plan,:faces) ? plan.faces : Matrix{UInt8}(undef,0,0))
end

# Both indexed B4 emitters write each primary row through this integer map.
# Check that map against the completed constructor before another view uses it.
function _extrude_nonew_b4_column_rows(sweep::_ExtrudeNoNewIndexedSweep,
        catalog::Union{_ExtrudeNoNewB4StripCatalog,_ExtrudeNoNewB4Free.Catalog},caller)
    levels,nodes=size(sweep.cols)
    nodes==length(catalog.source.source_coordinates) || throw(ArgumentError(
        "$caller: indexed B4 Source column count differs from its constructor"))
    rows=Matrix{Int32}(undef,levels,nodes)
    for node in 1:nodes,level in 1:levels
        row=_extrude_nonew_rect_grid_node(Int32(node),level,levels)
        0<row<=nnodes(sweep.volume) &&
            isequal(_extrude_nonew_coordinate(sweep.volume,row),sweep.cols[level,node]) ||
            throw(ArgumentError("$caller: indexed B4 column loses its constructor row"))
        rows[level,node]=row
    end
    return rows
end

# Native Source coordinates and emitted Source birth values are separate
# inputs. Reuse the completed topology decisions while deriving and certifying
# the emitted columns and body means from the latter, through explicit UIDs.
function _extrude_nonew_b4_geometry(m::GeoModel,t::Int,completed::_ExtrudeNoNewCompletePlan,
                                  source_point::F,caller) where F
    catalog=completed.catalog
    catalog isa Union{_ExtrudeNoNewB4StripCatalog,_ExtrudeNoNewB4Free.Catalog} || return nothing
    link=get(m.meshing.extrude_sources,(3,t),nothing)
    link!==nothing && link[1]==2 || throw(ArgumentError(
        "$caller: emitted B4 sweep lacks its Source provenance"))
    source_tag=abs(link[2]);source=completed.surfaces[source_tag]
    source isa MixedMesh && source.entity_data!==nothing || throw(ArgumentError(
        "$caller: emitted B4 Source lacks its actual column UIDs"))
    tags=source.entity_data.external_node_tags
    length(tags)==size(completed.columns,2) && allunique(tags) || throw(ArgumentError(
        "$caller: emitted B4 Source UIDs do not align with its columns"))
    coordinates=copy(source.coords);changed=false
    for (row,uid) in enumerate(tags)
        point=source_point(uid)
        point===nothing && continue
        old=(coordinates[1,row],coordinates[2,row],coordinates[3,row])
        isequal(old,point) && continue
        coordinates[1,row]=point[1];coordinates[2,row]=point[2];coordinates[3,row]=point[3]
        changed=true
    end
    changed || return nothing
    geometry_source=MixedMesh(coordinates,source.blocks)
    params,spec,_=_extrude_entity_params(m,3,t,caller)
    original=catalog.source
    geometry_data=_extrude_nonew_rect_grid_source(m,source_tag,geometry_source,spec,caller;
                                                mode=:b4_strip,bound=original)
    geometry_data.source_cells==original.source_cells &&
        geometry_data.edges==original.edges && geometry_data.edge_cells==original.edge_cells &&
        geometry_data.normal_axis==original.normal_axis &&
        geometry_data.normal_direction==original.normal_direction &&
        geometry_data.source_orientation==original.source_orientation || throw(ArgumentError(
            "$caller: emitted B4 Source changes its bound topology or orientation"))
    geometry_catalog=if catalog isa _ExtrudeNoNewB4StripCatalog
        _ExtrudeNoNewB4StripCatalog(geometry_data,catalog.top_states,catalog.cell_counts,
            catalog.face_capacity,catalog.diagonal_count,catalog.recombined,catalog.levels,
            catalog.layer_refs)
    else
        _ExtrudeNoNewB4Free.Catalog(geometry_data,catalog.top_states,catalog.template_indices,
            catalog.problem_positions,catalog.problem_masks,catalog.cell_counts,
            catalog.face_capacity,catalog.diagonal_count,catalog.recombined,catalog.levels,
            catalog.layer_refs)
    end
    columns=_extrude_volume_columns(m,t,source_tag,geometry_source,params,spec,catalog.levels,caller)
    _extrude_nonew_rect_grid_product_certify(columns,geometry_catalog,spec,caller)
    geometry=if catalog isa _ExtrudeNoNewB4StripCatalog
        top,laterals=_extrude_nonew_region_boundary(m,t,source_tag,caller;sides=4)
        _extrude_nonew_b4_strip_finish(m,t,params,source_tag,geometry_source,
            geometry_catalog,columns,top,laterals,caller).sweep.volume
    else
        _ExtrudeNoNewB4Free.emit(geometry_catalog,columns,spec,caller).volume
    end
    actual=completed.volume
    nnodes(geometry)==nnodes(actual) && length(geometry.blocks)==length(actual.blocks) &&
        all(a.msh==b.msh && a.nodes==b.nodes && a.tags==b.tags
            for (a,b) in zip(geometry.blocks,actual.blocks)) || throw(ArgumentError(
                "$caller: emitted B4 geometry changes its bound cell identities"))
    return geometry
end

# ─── Bounded shared-lateral slice (LAY01) ───────────────────────────────────
# Two adjacent QuadTriNoNewVerts sweeps may share exactly one literal lateral
# Surface when both record the identical linear translation operation. The
# shared Surface is emitted once — the first partner's plan binds its actual
# cells — and its corner/column incidence is transported into the second
# plan as a per-layer diagonal contract, so each Region's lateral boundary
# sees identical wall sub-edges. A joint oriented-incidence check then
# certifies opposite child faces, matching side diagonals, positive typed
# volumes and the closed exterior boundary. Everything outside this slice
# (third sharer, multiple or non-literal sharing, non-lateral partners,
# unequal or nonlinear sweeps, non-generic catalogs) keeps the original
# loud diagnostic.

@noinline _extrude_nonew_shared_planner_error(caller,surface)=
    throw(ArgumentError(
        "$caller: QuadTriNoNewVerts shared Surface[$surface] " *
        "requires the shared-region planner"))

function _extrude_nonew_shell_sign(m::GeoModel,t::Int,s::Int)
    for shell in m.volumes[t],signed in m.surface_loops[abs(shell)]
        abs(signed)==s && return sign(signed)
    end
    return 0
end

# One boundary Surface of `t` is bounded by more than `t`. Admit only the
# literal pair: `t` plus exactly one other layered QuadTriNoNewVerts Region
# created by the same translation, with `s` a lateral generatrix product
# owned by both and the shared Curve a boundary of both Sources.
function _extrude_nonew_shared_validate(m::GeoModel,t::Int,s::Int,
                                        adjacent,caller::AbstractString)
    length(adjacent)==2 && t in adjacent ||
        _extrude_nonew_shared_planner_error(caller,s)
    partner=adjacent[1]==t ? adjacent[2] : adjacent[1]
    _extrude_is_nonew(_extrude_gate(m,3,partner)) ||
        _extrude_nonew_shared_planner_error(caller,s)
    entry=_extrude_entity_params(m,3,t,caller)
    p_entry=_extrude_entity_params(m,3,partner,caller)
    entry!==nothing && p_entry!==nothing ||
        _extrude_nonew_shared_planner_error(caller,s)
    params,spec,link=entry
    p_params,p_spec,p_link=p_entry
    link[1]==2 && p_link[1]==2 ||
        _extrude_nonew_shared_planner_error(caller,s)
    spec.type===:translate && p_spec.type===:translate ||
        _extrude_nonew_shared_planner_error(caller,s)
    isequal(params,p_params) && isequal(spec,p_spec) ||
        _extrude_nonew_shared_planner_error(caller,s)
    s_link=get(m.meshing.extrude_sources,(2,s),nothing)
    s_link!==nothing && s_link[1]==1 ||
        _extrude_nonew_shared_planner_error(caller,s)
    isequal(get(m.meshing.extrude,(2,s),nothing),params) &&
        isequal(get(m.meshing.extrude_specs,(2,s),nothing),spec) ||
        _extrude_nonew_shared_planner_error(caller,s)
    source=abs(link[2]);p_source=abs(p_link[2])
    s!=source && s!=p_source ||
        _extrude_nonew_shared_planner_error(caller,s)
    generatrix=abs(s_link[2])
    generatrix in _model_projection_surface_curves(m,source) &&
        generatrix in _model_projection_surface_curves(m,p_source) ||
        _extrude_nonew_shared_planner_error(caller,s)
    sign_t=_extrude_nonew_shell_sign(m,t,s)
    sign_t!=0 && sign_t==-_extrude_nonew_shell_sign(m,partner,s) ||
        _extrude_nonew_shared_planner_error(caller,s)
    return partner
end

# The bounded shared laterals of one Region — surface => partner. The
# current slice admits at most one shared lateral per Region, and only
# inside a literal two-Region pair.
function _extrude_nonew_shared_for(m::GeoModel,t::Int,caller::AbstractString)
    shared=Dict{Int,Int}()
    for shell in m.volumes[t],signed in m.surface_loops[abs(shell)]
        s=abs(signed)
        adjacent=_extrude_quadtri_regions(m,s,caller)
        length(adjacent)==1 && continue
        shared[s]=_extrude_nonew_shared_validate(m,t,s,adjacent,caller)
    end
    length(shared)<=1 || _extrude_nonew_shared_planner_error(
        caller,minimum(keys(shared)))
    return shared
end

# The shared generatrix occupies one side of the completed partner's Source
# cell; locate that side by its level-one endpoint positions.
function _extrude_nonew_shared_side(m::GeoModel,
        completed::_ExtrudeNoNewCompletePlan,generatrix::Int,s::Int,
        caller::AbstractString)
    cell=completed.catalog isa _ExtrudeNoNewCatalog ?
        completed.catalog.source_cell : nothing
    cell!==nothing && length(cell)==4 ||
        _extrude_nonew_shared_planner_error(caller,s)
    chain=_extrude_curve_nodes(m,generatrix,caller)
    length(chain)==2 || _extrude_nonew_shared_planner_error(caller,s)
    pair=_extrude_ekey(chain[1],chain[2])
    side=0
    for k in 1:4
        _extrude_ekey(completed.columns[1,cell[k]],
                      completed.columns[1,cell[mod1(k+1,4)]])==pair || continue
        side==0 || _extrude_nonew_shared_planner_error(caller,s)
        side=k
    end
    side==0 && _extrude_nonew_shared_planner_error(caller,s)
    return side
end

# The per-layer lateral diagonal the completed owner actually emits on the
# shared wall — `nothing` where the wall stays a whole Quad (RecombLaterals
# or an all-boundary problem layer). This is the geometric contract the
# partner's template chain is constrained to reproduce.
function _extrude_nonew_shared_diagonals(completed::_ExtrudeNoNewCompletePlan,
        side::Int,s::Int,caller::AbstractString)
    isempty(completed.faces) && return nothing
    cell=completed.catalog.source_cell
    face=_EXTRUDE_NONEW_HEX_FACES[side]
    intervals=size(completed.faces,2)
    required=Vector{Union{Nothing,NTuple{2,NTuple{3,Float64}}}}(undef,intervals)
    for layer in 1:intervals
        state=completed.faces[side,layer]
        if state==0
            required[layer]=nothing
            continue
        end
        state in (1,2) || _extrude_nonew_shared_planner_error(caller,s)
        v=_extrude_nonew_corners(completed.columns,cell,layer)
        required[layer]=state==1 ? _extrude_ekey(v[face[1]],v[face[3]]) :
                                   _extrude_ekey(v[face[2]],v[face[4]])
    end
    return required
end

# Two emissions of the shared Surface are deterministic over the same
# generatrix chain, transform and diagonal set, so admissible parts are
# bitwise identical — `==` on the mesh structs would compare array
# identity, not content.
function _extrude_nonew_shared_part_equal(a,b)
    if a isa Mesh || b isa Mesh
        return a isa Mesh && b isa Mesh && a.coords==b.coords &&
            a.segs==b.segs && a.tris==b.tris && a.tets==b.tets &&
            a.seg_tag==b.seg_tag && a.tri_tag==b.tri_tag && a.tet_tag==b.tet_tag
    end
    b isa MixedMesh || return false
    a.coords==b.coords && length(a.blocks)==length(b.blocks) || return false
    for (x,y) in zip(a.blocks,b.blocks)
        x===y && continue
        x isa ElementBlock && y isa ElementBlock &&
            x.msh==y.msh && x.nodes==y.nodes && x.tags==y.tags || return false
    end
    return true
end

# Oriented face incidence in one joint position registry — the shared wall
# welds bitwise, so each wall facet must arrive once per Region with
# opposite outward orientation.
function _extrude_nonew_shared_count_block!(counts,nodes,msh,coords,lookup,caller)
    uids=Matrix{Int32}(undef,size(nodes,1),size(nodes,2))
    for index in eachindex(nodes)
        node=nodes[index]
        coordinate=ntuple(k->_model_projection_coordinate_key(coords[k,node]),3)
        uid=get(lookup,coordinate,Int32(0))
        uid>0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared Region volume node is absent from the joint registry"))
        uids[index]=uid
    end
    _extrude_nonew_count_faces!(counts,uids,msh,caller)
    return nothing
end

function _extrude_nonew_shared_count_faces!(counts,volume,lookup,caller)
    if volume isa Mesh
        _extrude_nonew_shared_count_block!(counts,volume.tets,4,volume.coords,lookup,caller)
    else
        for block in volume.blocks
            block isa ElementBlock || continue
            msh_dimension(block.msh)==3 || continue
            _extrude_nonew_shared_count_block!(counts,block.nodes,block.msh,volume.coords,lookup,caller)
        end
    end
    return nothing
end

function _extrude_nonew_shared_positive(volume,caller)
    if volume isa Mesh
        for cell in axes(volume.tets,2)
            _extrude_cell_signed_volume(volume.coords,volume.tets,cell,4)>0 ||
                throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts shared Region volume cell is not positive"))
        end
    else
        for block in volume.blocks
            block isa ElementBlock && msh_dimension(block.msh)==3 || continue
            for cell in axes(block.nodes,2)
                _extrude_cell_signed_volume(volume.coords,block.nodes,cell,block.msh)>0 ||
                    throw(ArgumentError(
                        "$caller: QuadTriNoNewVerts shared Region volume cell is not positive"))
            end
        end
    end
    return nothing
end

# Joint pair certificate in one position registry: each Region counts its
# own oriented volume faces, the welded wall is exactly the facets both
# Regions emitted — each once, with opposite outward orientation — the
# bound shared Surface's faces must equal that welded set, and the
# remaining singly-incident faces must equal the two Regions' exterior
# surfaces. Every emitted volume cell keeps a positive typed volume.
function _extrude_nonew_shared_certify(first::_ExtrudeNoNewCompletePlan,
        second::_ExtrudeNoNewCompletePlan,s::Int,caller::AbstractString)
    lookup=Dict{NTuple{3,Float64},Int32}()
    for volume in (first.volume,second.volume),node in 1:nnodes(volume)
        coordinate=_extrude_nonew_coordinate(volume,node)
        haskey(lookup,coordinate) && continue
        lookup[coordinate]=Int32(length(lookup)+1)
    end
    counts_first=Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence}()
    counts_second=Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence}()
    _extrude_nonew_shared_count_faces!(counts_first,first.volume,lookup,caller)
    _extrude_nonew_shared_count_faces!(counts_second,second.volume,lookup,caller)
    welded=Set{NTuple{4,Int32}}()
    for (key,incidence) in counts_first
        other=get(counts_second,key,nothing)
        other===nothing && continue
        # A wall facet is one boundary face of each Region — a repeated
        # same-position emission would overlap the wall instead.
        incidence.count==0x01 && other.count==0x01 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared wall face $key is emitted " *
            "more than once by one Region"))
        other.orientation==xor(incidence.orientation,0x01) ||
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared wall face $key has " *
                "equal or inconsistent orientation"))
        push!(welded,key)
    end
    bound=Set{NTuple{4,Int32}}()
    _extrude_nonew_surface_faces!(bound,first.surfaces[s],lookup,s,caller)
    bound==welded || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts shared Surface[$s] does not exactly " *
        "cover the internal wall between its two Regions"))
    external=Set{NTuple{4,Int32}}()
    for (key,incidence) in counts_first
        haskey(counts_second,key) && continue
        incidence.count==0x01 && push!(external,key)
    end
    for (key,incidence) in counts_second
        haskey(counts_first,key) && continue
        incidence.count==0x01 && push!(external,key)
    end
    exterior=Set{NTuple{4,Int32}}()
    for (tag,part) in first.surfaces
        tag==s && continue
        _extrude_nonew_surface_faces!(exterior,part,lookup,tag,caller)
    end
    for (tag,part) in second.surfaces
        tag==s && continue
        _extrude_nonew_surface_faces!(exterior,part,lookup,tag,caller)
    end
    exterior==external || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts shared Region exterior surfaces do not " *
        "exactly cover the external boundary"))
    _extrude_nonew_shared_positive(first.volume,caller)
    _extrude_nonew_shared_positive(second.volume,caller)
    return nothing
end

function _extrude_nonew_scope(m::GeoModel,caller::AbstractString;
                               tags=keys(m.volumes),
                               min_angle_deg::Real=25.0,
                               max_periodic_passes=8,
                               size_field::Union{Nothing,AbstractSizeField}=nothing,
                               _working_model::Bool=false,_mesh_input=nothing,
                               _trace_columns::Bool=false)
    selected=Int[]
    for tag in tags
        _extrude_is_nonew(_extrude_gate(m,3,tag)) && push!(selected,tag)
    end
    isempty(selected) && return nothing
    sort!(selected)
    working=_working_model ? m : deepcopy(m)
    # Discover every bounded shared-lateral pair before any planning, so
    # admission failures stay atomic. `shared_of` is Region => (surface,
    # partner); each Region participates in at most one pair.
    selected_set=Set(selected)
    shared_of=Dict{Int,NTuple{2,Int}}()
    for tag in selected
        for (s,partner) in _extrude_nonew_shared_for(working,tag,caller)
            partner in selected_set ||
                _extrude_nonew_shared_planner_error(caller,s)
            shared_of[tag]=(s,partner)
        end
    end
    for (tag,(s,partner)) in shared_of
        get(shared_of,partner,nothing)==(s,tag) ||
            _extrude_nonew_shared_planner_error(caller,s)
    end
    volumes=Dict{Int,_ExtrudeNoNewCompletePlan}()
    surfaces=Dict{Int,Union{Mesh,MixedMesh}}()
    for tag in selected
        info=get(shared_of,tag,nothing)
        diagonals=nothing
        if info!==nothing && info[2]<tag
            # The earlier partner already bound the shared Surface: read
            # its emitted wall diagonals and constrain this plan to them.
            s,partner=info
            owner=volumes[partner]
            generatrix=abs(working.meshing.extrude_sources[(2,s)][2])
            side=_extrude_nonew_shared_side(working,owner,generatrix,s,caller)
            diagonals=_extrude_nonew_shared_diagonals(owner,side,s,caller)
        end
        completed=_extrude_nonew_complete_plan(working,tag,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field,_working_model=true,_mesh_input=_mesh_input,
            _trace_columns=_trace_columns,_shared_diagonals=diagonals)
        info===nothing || completed.catalog isa _ExtrudeNoNewCatalog ||
            _extrude_nonew_shared_planner_error(caller,info[1])
        volumes[tag]=completed
        for (surface,part) in completed.surfaces
            if haskey(surfaces,surface)
                # Only the certified shared lateral may repeat, and only
                # when both Regions emitted bitwise-identical cells.
                (info!==nothing && surface==info[1] &&
                 _extrude_nonew_shared_part_equal(part,surfaces[surface])) ||
                    _extrude_nonew_shared_planner_error(caller,surface)
                continue
            end
            surfaces[surface]=part
        end
    end
    for tag in selected
        info=get(shared_of,tag,nothing)
        info===nothing && continue
        # Certify each pair once, from the earlier Region's plan.
        tag<info[2] || continue
        s,partner=info
        _extrude_nonew_shared_certify(volumes[tag],volumes[partner],s,caller)
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
    owners=Int[]
    for region in _extrude_quadtri_regions(m,t,caller)
        _extrude_is_nonew(_extrude_gate(m,3,region)) || continue
        push!(owners,region)
    end
    if length(owners)==2
        # The bounded literal-lateral pair plans its scope jointly; any
        # other multi-owner incidence keeps the original diagnostic.
        a,b=owners
        get(_extrude_nonew_shared_for(m,a,caller),t,0)==b &&
            get(_extrude_nonew_shared_for(m,b,caller),t,0)==a ||
            _extrude_nonew_shared_planner_error(caller,t)
        scope===nothing || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts Surface[$t] is absent from this mesh scope"))
        built=_extrude_nonew_scope(m,caller;tags=owners,
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field)
        return built.surfaces[t]
    end
    length(owners)>2 &&
        _extrude_nonew_shared_planner_error(caller,t)
    isempty(owners) && return nothing
    owner=owners[1]
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
