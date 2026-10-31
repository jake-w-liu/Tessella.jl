# A completed NoNew plan lives only within the mesh operation that requested
# it. In particular, no model/session cache survives edits to Layers, source
# discretization or geometry.
struct _ExtrudeNoNewCompletePlan
    volume::Union{Mesh,MixedMesh}
    surfaces::Dict{Int,Union{Mesh,MixedMesh}}
    columns::Matrix{NTuple{3,Float64}}
    catalog::Union{_ExtrudeNoNewCatalog,_ExtrudeNoNewTwoTriCatalog,_ExtrudeNoNewQuadPatchCatalog,_ExtrudeNoNewQuadStripCatalog,_ExtrudeNoNewThreeQuadStripCatalog,_ExtrudeNoNewFourQuadStripCatalog,_ExtrudeNoNewRectGridCatalog,_ExtrudeNoNewB4StripCatalog,_ExtrudeNoNewB4Free.Catalog,_ExtrudeNoNewMultiCellCatalog}
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
        # A mesh-level `Reverse Surface` flag only flips the published cell
        # windings — the scope mirrors a reversed Source through its Source
        # mesh (the surface pass already flipped it) and flips each emitted
        # cap/lateral part at collection, so the flag is supported here.
        for loop in m.surfaces[surface],signed in m.loops[loop]
            curve=abs(signed)
            push!(curves,curve)
            a,b=m.curves[curve]
            push!(points,a);push!(points,b)
        end
    end
    # `Reverse Volume`/`Reverse Curve` flip only the published element
    # windings: the volume emitter reverses its completed cells after
    # certification and the curve segments' winding lives on the ordinary
    # curve-meshing path outside this scope.
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

# The generic unmerged-NoNew identity binder certifies every catalog family.
# Each family provides its own face-bucket capacity and whether its typed
# shell carries oriented interior macrofaces; mirror the per-family dispatch
# used by `_extrude_nonew_complete_plan`.
function _extrude_nonew_catalog_face_capacity(catalog,volume)
    catalog isa _ExtrudeNoNewQuadPatchCatalog &&
        return _extrude_nonew_quad_patch_face_capacity(catalog)
    catalog isa _ExtrudeNoNewQuadStripCatalog &&
        return _extrude_nonew_quad_strip_face_capacity(catalog)
    catalog isa _ExtrudeNoNewThreeQuadStripCatalog &&
        return _extrude_nonew_three_quad_strip_face_capacity(catalog)
    catalog isa _ExtrudeNoNewFourQuadStripCatalog &&
        return _extrude_nonew_four_quad_strip_face_capacity(catalog)
    catalog isa _ExtrudeNoNewDynamicGridCatalog &&
        return _extrude_nonew_dynamic_grid_face_capacity(catalog)
    catalog isa _ExtrudeNoNewTwoTriCatalog &&
        return (volume isa Mesh ? 16 : 7)*length(catalog.layer_refs)+2
    catalog isa _ExtrudeNoNewMultiCellCatalog &&
        return 6*length(catalog.source_cells)*length(catalog.layer_refs)+2
    return 0
end

_extrude_nonew_catalog_oriented_internal(catalog)=
    catalog isa Union{_ExtrudeNoNewTwoTriCatalog,_ExtrudeNoNewQuadPatchCatalog,
        _ExtrudeNoNewQuadStripCatalog,_ExtrudeNoNewThreeQuadStripCatalog,
        _ExtrudeNoNewFourQuadStripCatalog,_ExtrudeNoNewDynamicGridCatalog,
        _ExtrudeNoNewMultiCellCatalog}

# The admitted unmerged-NoNew families: every catalog the generic identity
# binder can certify. The generic single-quad catalog carries no face capacity.
const _EXTRUDE_NONEW_UNMERGED_CATALOG_TYPES=Union{
    _ExtrudeNoNewCatalog,_ExtrudeNoNewTwoTriCatalog,
    _ExtrudeNoNewQuadPatchCatalog,_ExtrudeNoNewQuadStripCatalog,
    _ExtrudeNoNewThreeQuadStripCatalog,_ExtrudeNoNewFourQuadStripCatalog,
    _ExtrudeNoNewRectGridCatalog,_ExtrudeNoNewB4StripCatalog,
    _ExtrudeNoNewB4Free.Catalog,_ExtrudeNoNewMultiCellCatalog}

function _extrude_nonew_complete_plan(m::GeoModel,t::Int,
                                      caller::AbstractString;
                                      min_angle_deg::Real=25.0,
                                      max_periodic_passes=8,
                                      size_field::Union{Nothing,AbstractSizeField}=nothing,
                                      _working_model::Bool=false,_mesh_input=nothing,
                                      _trace_columns::Bool=false,
                                      _shared_diagonals=nothing,
                                      _shared_parts=nothing)
    working=_working_model ? m : deepcopy(m)
    _extrude_nonew_scope_contract(working,t,caller)
    plan=_extrude_nonew_plan(working,t,caller;min_angle_deg=min_angle_deg,
        max_periodic_passes=max_periodic_passes,size_field=size_field,_mesh_input=_mesh_input,
        _shared_diagonals=_shared_diagonals,_shared_parts=_shared_parts)
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
    multi_cell=plan.catalog isa _ExtrudeNoNewMultiCellCatalog
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
        (volume isa Mesh ? 16 : 7)*length(plan.catalog.layer_refs)+2 : multi_cell ?
        6*length(plan.catalog.source_cells)*length(plan.catalog.layer_refs)+2 : 0
    if retained_b4
        _extrude_nonew_retained_b4_certificate(working,t,volume,surfaces,plan.catalog,caller;
            face_capacity)
    else
        _extrude_nonew_certify_boundary(volume,surfaces,caller;
            oriented_internal=two_tri || quad_patch || quad_strip || three_quad_strip || four_quad_strip || dynamic_grid || multi_cell,
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

# ─── Shared-surface slice (LAY01) ───────────────────────────────────────────
# Adjacent QuadTriNoNewVerts sweeps may share literal lateral Surfaces when
# every NoNewVerts owner records the identical linear translation operation.
# Upstream (Gmsh 4.15.2 oracle) admits one literal Surface bounded by three
# or more such Regions, a lateral shared with an ordinary layered Region
# (a recombined partner keeps the welded wall in whole quads and sweeps
# around it; an unrecombined partner welds the split wall its prism
# subdivision emits), and a copied cap that is literally the next stacked
# sweep's Source. Each shared
# Surface is emitted once: the earliest NoNewVerts owner's plan binds the
# wall's actual cells on a pure-NoNewVerts wall and transports its per-layer
# diagonals into each later plan as a (side,state) contract, while a wall
# shared with a recombined ordinary partner carries no diagonal — every
# NoNewVerts owner fans the affected layers instead. An unrecombined
# partner keeps the bound wall's split facets: its emitted diagonals seed
# the global `SubdivideExtrudedMesh` edge set (and its swap protection), so
# the partner's own phase-3 tetrahedra arrive on the welded triangles. A
# shared copied cap is simply the top owner's emitted cap, reused bitwise
# as the next sweep's Source mesh. A joint oriented-incidence check
# certifies each wall's per-owner child faces, matching surface cells,
# positive typed volumes and the closed exterior boundary — with three or
# more owners every wall facet must arrive once per Region with mixed
# outward orientation. Everything outside this slice (non-literal sharing,
# partners whose recombination cannot reproduce the welded wall — such as
# mixed recombined/unrecombined owner sets on one wall — unequal sweep
# records on pure walls, unequal/bent cap transforms, non-generic catalogs
# that cannot reproduce the bound emission, Sources welded onto a
# partner's lateral) keeps the original loud diagnostic.

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

# Geometric side of one shared lateral's owners: the wall plane contains
# the generatrix and the sweep translation, so the generatrix × sweep
# normal splits space into the two wall sides and each owner's Source
# boundary centroid carries one signed side. A reversed Source loop (or a
# mesh-level Source reversal) makes its whole Region inside-out — every
# lateral shell sign flips with it — while the Region still sits on the
# opposite wall side upstream welds. Returning the per-owner signed side
# lets the shared validation admit exactly that orientation reversal and
# keep rejecting a literal same-side incidence (overlapping owners).
# `nothing` marks an owner whose Source interior cannot be placed.
function _extrude_nonew_shared_sides(m::GeoModel,adjacent,
                                     generatrix::Int,spec,
                                     caller::AbstractString)
    haskey(m.curves,generatrix) || return nothing
    a,b=m.curves[generatrix]
    (haskey(m.points,a) && haskey(m.points,b)) || return nothing
    ga=m.points[a];gb=m.points[b]
    nx=(gb[2]-ga[2])*spec.T[3]-(gb[3]-ga[3])*spec.T[2]
    ny=(gb[3]-ga[3])*spec.T[1]-(gb[1]-ga[1])*spec.T[3]
    nz=(gb[1]-ga[1])*spec.T[2]-(gb[2]-ga[2])*spec.T[1]
    (nx!=0.0 || ny!=0.0 || nz!=0.0) || return nothing
    sides=Int[]
    for owner in adjacent
        entry=_extrude_entity_params(m,3,owner,caller)
        entry===nothing && return nothing
        _,_,o_link=entry
        o_link[1]==2 || return nothing
        src=abs(o_link[2])
        curves=_model_projection_surface_curves(m,src)
        isempty(curves) && return nothing
        cx=cy=cz=0.0
        for c in curves
            haskey(m.curves,abs(c)) || return nothing
            ca,cb=m.curves[abs(c)]
            corner=c>0 ? ca : cb
            haskey(m.points,corner) || return nothing
            p=m.points[corner]
            cx+=p[1];cy+=p[2];cz+=p[3]
        end
        k=length(curves)
        d=(cx/k-ga[1])*nx+(cy/k-ga[2])*ny+(cz/k-ga[3])*nz
        push!(sides,d>0.0 ? 1 : d<0.0 ? -1 : 0)
    end
    return sides
end

# One boundary Surface of `t` is bounded by more than `t`. Two classes are
# admitted. A shared copied cap (provenance `(2,·)`) is literally one
# owner's copied top and the other's Source — upstream stacks Layered
# sweeps this way and emits the cap once. A shared lateral (provenance
# `(1,·)`) is a generatrix product welded across every owner: each
# QuadTriNoNewVerts owner must record the identical operation, and an
# ordinary layered partner may join the incidence when its layer grid and
# transform match — it then sweeps around the whole-quad wall upstream
# keeps for it. The Surface's own recorded operation must literally be one
# of its owners', and owners must sit on both sides of the wall.
function _extrude_nonew_shared_validate(m::GeoModel,t::Int,s::Int,
                                        adjacent,caller::AbstractString)
    length(adjacent)>=2 && t in adjacent ||
        _extrude_nonew_shared_planner_error(caller,s)
    s_link=get(m.meshing.extrude_sources,(2,s),nothing)
    entry=_extrude_entity_params(m,3,t,caller)
    entry!==nothing || _extrude_nonew_shared_planner_error(caller,s)
    params,spec,link=entry
    link[1]==2 || _extrude_nonew_shared_planner_error(caller,s)
    if s_link!==nothing && s_link[1]==2
        # Shared copied cap: exactly the top owner (whose Source the cap
        # records) and the Source owner (whose literal Source it is).
        # Upstream stacks sweeps of different layer records here but the
        # transforms must agree — a bent stack loses its top face upstream
        # and keeps the planner diagnostic.
        length(adjacent)==2 || _extrude_nonew_shared_planner_error(caller,s)
        partner=adjacent[1]==t ? adjacent[2] : adjacent[1]
        _extrude_is_nonew(_extrude_gate(m,3,partner)) ||
            _extrude_nonew_shared_planner_error(caller,s)
        p_entry=_extrude_entity_params(m,3,partner,caller)
        p_entry!==nothing || _extrude_nonew_shared_planner_error(caller,s)
        p_params,p_spec,p_link=p_entry
        p_link[1]==2 || _extrude_nonew_shared_planner_error(caller,s)
        src_t=abs(link[2])==s;src_p=abs(p_link[2])==s
        xor(src_t,src_p) || _extrude_nonew_shared_planner_error(caller,s)
        top_t=abs(link[2])==abs(s_link[2])
        top_p=abs(p_link[2])==abs(s_link[2])
        xor(top_t,top_p) || _extrude_nonew_shared_planner_error(caller,s)
        (src_t ? top_p : top_t) ||
            _extrude_nonew_shared_planner_error(caller,s)
        top_spec=src_t ? p_spec : spec
        top_params=src_t ? p_params : params
        isequal(spec,p_spec) &&
            isequal(get(m.meshing.extrude,(2,s),nothing),top_params) &&
            isequal(get(m.meshing.extrude_specs,(2,s),nothing),top_spec) ||
            _extrude_nonew_shared_planner_error(caller,s)
        sign_t=_extrude_nonew_shell_sign(m,t,s)
        sign_t!=0 && sign_t==-_extrude_nonew_shell_sign(m,partner,s) ||
            _extrude_nonew_shared_planner_error(caller,s)
        return adjacent
    end
    # Shared lateral: a literal generatrix product welded across its
    # owners. Three or more NoNewVerts owners are admitted (upstream welds
    # a non-manifold wall all of them split identically); ordinary layered
    # partners may join — they recombine the welded quadrangles upstream
    # refuses to split (`Cannot extrude quadrangles without Recombine`)
    # and sweep their own cells around the wall.
    s_link!==nothing && s_link[1]==1 ||
        _extrude_nonew_shared_planner_error(caller,s)
    # The literal lateral slice stays linear-translational — upstream's
    # rotated or bent sweeps lose their shared wall on a different entity
    # and keep the planner diagnostic.
    spec.type===:translate || _extrude_nonew_shared_planner_error(caller,s)
    source=abs(link[2])
    # A literal generatrix product is never also an extrusion Source — a
    # Region sourcing on a partner's lateral keeps the planner diagnostic.
    s==source && _extrude_nonew_shared_planner_error(caller,s)
    generatrix=abs(s_link[2])
    generatrix in _model_projection_surface_curves(m,source) ||
        _extrude_nonew_shared_planner_error(caller,s)
    s_params=get(m.meshing.extrude,(2,s),nothing)
    s_spec=get(m.meshing.extrude_specs,(2,s),nothing)
    owned=isequal(s_params,params) && isequal(s_spec,spec)
    signs=Int[_extrude_nonew_shell_sign(m,t,s)]
    for partner in adjacent
        partner==t && continue
        p_entry=_extrude_entity_params(m,3,partner,caller)
        p_entry===nothing && _extrude_nonew_shared_planner_error(caller,s)
        p_params,p_spec,p_link=p_entry
        p_link[1]==2 || _extrude_nonew_shared_planner_error(caller,s)
        p_source=abs(p_link[2])
        s==p_source && _extrude_nonew_shared_planner_error(caller,s)
        generatrix in _model_projection_surface_curves(m,p_source) ||
            _extrude_nonew_shared_planner_error(caller,s)
        if _extrude_is_nonew(p_params)
            isequal(p_params,params) && isequal(p_spec,spec) ||
                _extrude_nonew_shared_planner_error(caller,s)
        else
            # An ordinary layered partner (`Layers+Recombine`) or an
            # AddVerts partner keeps its whole-quad wall faces upstream
            # refuses to split. An unrecombined partner welds the split
            # wall its own prism subdivision emits — the oracle's mixed
            # NoNewVerts/ordinary shared lateral arrives on the partner's
            # diagonal triangles, which `_extrude_subdivide!` receives as
            # seeded edges. Either way the layer grid and translation must
            # match so the welded wall is geometrically identical.
            ((p_params.recombine && p_params.quad_to_tri in (:none,:add_verts)) ||
             !p_params.recombine) &&
                isequal(p_params.layers,params.layers) &&
                isequal(p_params.heights,params.heights) &&
                isequal(p_params.scale_last,params.scale_last) &&
                isequal(p_spec,spec) ||
                _extrude_nonew_shared_planner_error(caller,s)
        end
        owned=owned || (isequal(s_params,p_params) && isequal(s_spec,p_spec))
        push!(signs,_extrude_nonew_shell_sign(m,partner,s))
    end
    owned || _extrude_nonew_shared_planner_error(caller,s)
    all(!iszero,signs) || _extrude_nonew_shared_planner_error(caller,s)
    if length(unique(signs))==1
        # Reversed-orientation Sources flip every lateral shell sign with
        # their inside-out Region — upstream welds the wall anyway. Admit
        # only when the owners' Source interiors still sit on both
        # geometric sides of the generatrix; a literal same-side
        # incidence keeps the planner diagnostic.
        sides=_extrude_nonew_shared_sides(m,adjacent,generatrix,spec,caller)
        (sides!==nothing && all(!iszero,sides) &&
         length(unique(sides))>1) ||
            _extrude_nonew_shared_planner_error(caller,s)
    end
    return adjacent
end

# The certified shared Surfaces of one Region — surface => its sorted
# complete Region owner list. A Region may bound several shared laterals
# and stacked caps; every wall still requires its own certified owner set.
function _extrude_nonew_shared_for(m::GeoModel,t::Int,caller::AbstractString)
    shared=Dict{Int,Vector{Int}}()
    for shell in m.volumes[t],signed in m.surface_loops[abs(shell)]
        s=abs(signed)
        adjacent=_extrude_quadtri_regions(m,s,caller)
        length(adjacent)==1 && continue
        shared[s]=_extrude_nonew_shared_validate(m,t,s,adjacent,caller)
    end
    return shared
end

# The shared-surface neighbor graph is closed under the owner relation:
# any operation that must plan one member plans the whole connected
# NoNewVerts component, so every bound wall keeps its joint certificate.
# Ordinary layered partners sweep their own plans outside the scope.
function _extrude_nonew_shared_closure(m::GeoModel,seeds,
                                       caller::AbstractString)
    closure=Set{Int}(seeds)
    queue=Int[seeds...]
    while !isempty(queue)
        tag=pop!(queue)
        for (_,owners) in _extrude_nonew_shared_for(m,tag,caller)
            for owner in owners
                owner in closure && continue
                _extrude_is_nonew(_extrude_gate(m,3,owner)) || continue
                push!(closure,owner)
                push!(queue,owner)
            end
        end
    end
    return sort!(collect(closure))
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

# The binding plan's emitted wall decides each layer's diagonal. The
# generic catalog keeps its (side,state) picks in `faces`; other catalogs
# (and the triangle cell, which uses the same struct with three corners)
# do not retain a hex face-state matrix, so their per-layer wall split is
# read back off the emitted Surface part instead — the wall quad on each
# layer is split by exactly the diagonal its two triangles share, or stays
# whole when the part keeps a quadrangle there.
# Per-segment, per-layer wall emission of a completed plan's bound part —
# one vector per generatrix chain segment, each holding `nothing` where
# that segment's wall quadrangle stays whole or its emitted diagonal's
# position pair. This is the multi-segment generalization the per-column
# family needs: a shared generatrix discretized into several chain nodes
# carries one independent quadrangle per segment.
function _extrude_nonew_wall_split_states(m::GeoModel,
        completed::_ExtrudeNoNewCompletePlan,s::Int,
        caller::AbstractString)
    catalog=completed.catalog
    if catalog isa _ExtrudeNoNewCatalog &&
            length(catalog.source_cell)==4 && !isempty(completed.faces)
        generatrix=abs(m.meshing.extrude_sources[(2,s)][2])
        side=_extrude_nonew_shared_side(m,completed,generatrix,s,caller)
        return [_extrude_nonew_shared_diagonals(completed,side,s,caller)]
    end
    part=get(completed.surfaces,s,nothing)
    part===nothing && return nothing
    params=get(m.meshing.extrude,(2,s),nothing)
    spec=get(m.meshing.extrude_specs,(2,s),nothing)
    s_link=get(m.meshing.extrude_sources,(2,s),nothing)
    (params===nothing || spec===nothing ||
     s_link===nothing || s_link[1]!=1) && return nothing
    generatrix=abs(s_link[2])
    chain=_extrude_curve_nodes(m,generatrix,caller)
    a,b=m.curves[generatrix]
    closed=a==b
    nseg=closed ? length(chain) : length(chain)-1
    nseg>=1 || _extrude_nonew_shared_planner_error(caller,s)
    us=_extrude_level_us(params)
    transforms=_extrude_level_transforms(spec,us,caller)
    intervals=length(us)-1
    cells=_extrude_nonew_wall_part_cells(part,s,caller)
    return [begin
        ga=chain[seg]
        gb=chain[closed ? mod1(seg+1,length(chain)) : seg+1]
        required=Vector{Union{Nothing,NTuple{2,NTuple{3,Float64}}}}(
            undef,intervals)
        for layer in 1:intervals
            quad=Set((_extrude_level_at(spec,ga,us[layer],transforms,layer),
                      _extrude_level_at(spec,gb,us[layer],transforms,layer),
                      _extrude_level_at(spec,ga,us[layer+1],transforms,layer+1),
                      _extrude_level_at(spec,gb,us[layer+1],transforms,layer+1)))
            member=Vector{NTuple{3,Float64}}[]
            for cell in cells
                all(p->p in quad,cell) && push!(member,cell)
            end
            required[layer]=_extrude_nonew_wall_quad_state(member,s,caller)
        end
        required
    end for seg in 1:nseg]
end

function _extrude_nonew_wall_diagonals(m::GeoModel,
        completed::_ExtrudeNoNewCompletePlan,s::Int,
        caller::AbstractString)
    states=_extrude_nonew_wall_split_states(m,completed,s,caller)
    states===nothing && return nothing
    # The generic catalog's transported contract encodes one diagonal per
    # layer; a multi-segment wall instead transports through the bound
    # Surface part to a per-column consumer, so this vector stays
    # single-segment.
    length(states)==1 || _extrude_nonew_shared_planner_error(caller,s)
    return states[1]
end

# Two emissions of the shared Surface are deterministic over the same
# generatrix chain, transform and diagonal set, so admissible parts are
# bitwise identical — `==` on the mesh structs would compare array
# identity, not content.
function _extrude_nonew_shared_part_equal(a,b)
    if a isa Mesh || b isa Mesh
        a isa Mesh && b isa Mesh && a.coords==b.coords &&
            a.segs==b.segs && a.tris==b.tris && a.tets==b.tets &&
            a.seg_tag==b.seg_tag && a.tri_tag==b.tri_tag && a.tet_tag==b.tet_tag &&
            return true
    elseif a.coords==b.coords && length(a.blocks)==length(b.blocks) &&
            all(zip(a.blocks,b.blocks)) do (x,y)
                x===y || (x isa ElementBlock && y isa ElementBlock &&
                    x.msh==y.msh && x.nodes==y.nodes && x.tags==y.tags)
            end
        return true
    end
    # A replayed emission can carry a different container or element tags —
    # the authority product of the copied cap versus its consumed Source
    # vector is the same cells and nodes either way. Compare the cell set
    # geometrically as a last resort.
    cells=Dict{Tuple{Int,Vector{NTuple{3,Float64}}},Int}()
    for (part,out) in ((a,1),(b,-1))
        if part isa Mesh
            for cell in axes(part.tris,2)
                key=(2,sort!([ntuple(k->part.coords[k,part.tris[j,cell]],3)
                              for j in 1:3]))
                cells[key]=get(cells,key,0)+(out==1 ? 1 : -1)
            end
            for cell in axes(part.tets,2)
                key=(4,sort!([ntuple(k->part.coords[k,part.tets[j,cell]],3)
                              for j in 1:4]))
                cells[key]=get(cells,key,0)+(out==1 ? 1 : -1)
            end
        else
            for block in part.blocks
                block isa ElementBlock && msh_dimension(block.msh)>=2 || continue
                width=size(block.nodes,1)
                for cell in axes(block.nodes,2)
                    key=(Int(block.msh),sort!([ntuple(
                        k->part.coords[k,block.nodes[j,cell]],3)
                        for j in 1:width]))
                    cells[key]=get(cells,key,0)+(out==1 ? 1 : -1)
                end
            end
        end
    end
    return all(iszero,values(cells))
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

# Joint certificate over the shared Surface's complete NoNewVerts owner
# set in one position registry: each Region counts its own oriented volume
# faces, the welded wall is exactly the facets every owner emitted — each
# once, with outward orientation on both wall sides — the bound shared
# Surface's faces must equal that welded set, and the remaining faces each
# singly incident in exactly one owner must equal the exterior surfaces.
# A facet a strict subset of the owners emits cannot be the certified
# wall and rejects. Every emitted volume cell keeps a positive typed
# volume.
function _extrude_nonew_shared_certify(plans::Vector{_ExtrudeNoNewCompletePlan},
        s::Int,caller::AbstractString)
    length(plans)>=2 || _extrude_nonew_shared_planner_error(caller,s)
    lookup=Dict{NTuple{3,Float64},Int32}()
    for plan in plans,node in 1:nnodes(plan.volume)
        coordinate=_extrude_nonew_coordinate(plan.volume,node)
        haskey(lookup,coordinate) && continue
        lookup[coordinate]=Int32(length(lookup)+1)
    end
    counts=[Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence}()
            for _ in plans]
    for index in eachindex(plans)
        _extrude_nonew_shared_count_faces!(
            counts[index],plans[index].volume,lookup,caller)
    end
    welded=Set{NTuple{4,Int32}}()
    external=Set{NTuple{4,Int32}}()
    seen=Set{NTuple{4,Int32}}()
    for owner in counts,(key,incidence) in owner
        key in seen && continue
        push!(seen,key)
        present=Int[]
        for index in eachindex(counts)
            haskey(counts[index],key) && push!(present,index)
        end
        if length(present)==length(plans)
            # A wall facet is one boundary face of each Region — a repeated
            # same-position emission would overlap the wall instead, and
            # three or more owners must still split it onto both sides.
            for index in present
                counts[index][key].count==0x01 || throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts shared wall face $key is " *
                    "emitted more than once by one Region"))
            end
            orientations=[counts[index][key].orientation for index in present]
            length(unique(orientations))>=2 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared wall face $key has " *
                "equal or inconsistent orientation"))
            push!(welded,key)
        elseif length(present)==1
            incidence.count==0x01 && push!(external,key)
        else
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared Region internal face " *
                "$key is bounded by a strict subset of its owners"))
        end
    end
    bound=Set{NTuple{4,Int32}}()
    _extrude_nonew_surface_faces!(bound,plans[1].surfaces[s],lookup,s,caller)
    bound==welded || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts shared Surface[$s] does not exactly " *
        "cover the internal wall between its Regions"))
    exterior=Set{NTuple{4,Int32}}()
    for plan in plans,(tag,part) in plan.surfaces
        tag==s && continue
        _extrude_nonew_surface_faces!(exterior,part,lookup,tag,caller)
    end
    exterior==external || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts shared Region exterior surfaces do not " *
        "exactly cover the external boundary"))
    for plan in plans
        _extrude_nonew_shared_positive(plan.volume,caller)
    end
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
    # Discover every certified shared Surface before any planning, so
    # admission failures stay atomic. `shared_of` is Region => (surface =>
    # complete sorted owner set) over all its shared walls and stacked
    # caps; every NoNewVerts owner of a certified Surface must join this
    # scope and record the same owner set symmetrically, while ordinary
    # layered partners plan through the ordinary sweep path instead.
    selected_set=Set(selected)
    shared_of=Dict{Int,Dict{Int,Vector{Int}}}()
    for tag in selected
        edges=_extrude_nonew_shared_for(working,tag,caller)
        isempty(edges) && continue
        for (s,owners) in edges
            for owner in owners
                owner==tag && continue
                _extrude_is_nonew(_extrude_gate(working,3,owner)) || continue
                owner in selected_set ||
                    _extrude_nonew_shared_planner_error(caller,s)
            end
        end
        shared_of[tag]=edges
    end
    for (tag,edges) in shared_of,(s,owners) in edges
        for owner in owners
            owner==tag && continue
            _extrude_is_nonew(_extrude_gate(working,3,owner)) || continue
            get(get(shared_of,owner,Dict{Int,Vector{Int}}()),s,Int[])==owners ||
                _extrude_nonew_shared_planner_error(caller,s)
        end
    end
    volumes=Dict{Int,_ExtrudeNoNewCompletePlan}()
    surfaces=Dict{Int,Union{Mesh,MixedMesh}}()
    for tag in selected
        edges=get(shared_of,tag,nothing)
        diagonals=nothing
        parts=nothing
        if edges!==nothing
            for (s,owners) in edges
                s_link=working.meshing.extrude_sources[(2,s)]
                # Copied caps carry no lateral contract — the Source
                # owner's plan consumes the top owner's emitted cap
                # bitwise. Walls shared with a recombined ordinary partner
                # stay in whole quads; every NoNewVerts owner fans them
                # instead.
                s_link[1]==1 || continue
                # The earliest NoNewVerts owner bound the shared Surface —
                # an unrecombined ordinary partner plans outside the scope
                # and holds no contract (its subdivision receives the
                # bound diagonals as seeded edges instead).
                binder_i=findfirst(owner->_extrude_is_nonew(
                    _extrude_gate(working,3,owner)),owners)
                binder_i===nothing && continue
                binder=owners[binder_i]
                if binder!=tag
                    t_link=working.meshing.extrude_sources[(3,tag)]
                    if _extrude_nonew_multi_cell_source(
                            working,abs(t_link[2]),caller)
                        # The per-column planner reproduces the bound wall
                        # facet-for-facet — each segment quadrangle
                        # independently, whole or split — so the binder's
                        # emitted part is the contract itself.
                        parts===nothing &&
                            (parts=Dict{Int,Union{Mesh,MixedMesh}}())
                        parts[s]=volumes[binder].surfaces[s]
                        continue
                    end
                end
                _extrude_nonew_wall_split(working,s,owners,caller) || continue
                binder==tag && continue
                # Only the generic catalog consumes a transported diagonal
                # contract — bounded/strip catalogs emit their wall
                # deterministically and the shared-part bitwise check
                # certifies the weld, so a multi-segment generatrix the
                # contract cannot encode stays weldable here.
                _extrude_nonew_generic_source(working,abs(t_link[2]),caller) ||
                    continue
                # The earliest NoNewVerts owner already bound the shared
                # Surface: read its emitted wall diagonals and constrain
                # this plan to them.
                owner=volumes[binder]
                required=_extrude_nonew_wall_diagonals(working,owner,s,caller)
                required===nothing &&
                    _extrude_nonew_shared_planner_error(caller,s)
                diagonals===nothing && (diagonals=Dict{Int,
                    Vector{Union{Nothing,NTuple{2,NTuple{3,Float64}}}}}())
                diagonals[s]=required
            end
        end
        completed=_extrude_nonew_complete_plan(working,tag,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field,_working_model=true,_mesh_input=_mesh_input,
            _trace_columns=_trace_columns,_shared_diagonals=diagonals,
            _shared_parts=parts)
        volumes[tag]=completed
        for (surface,part) in completed.surfaces
            if haskey(surfaces,surface)
                # Only a certified shared Surface may repeat, and only when
                # every owner emitted bitwise-identical cells.
                (edges!==nothing && haskey(edges,surface) &&
                 _extrude_nonew_shared_part_equal(part,surfaces[surface])) ||
                    _extrude_nonew_shared_planner_error(caller,surface)
                continue
            end
            surfaces[surface]=part
        end
    end
    for tag in selected
        edges=get(shared_of,tag,nothing)
        edges===nothing && continue
        for (s,owners) in edges
            owners=[owner for owner in owners
                    if _extrude_is_nonew(_extrude_gate(working,3,owner))]
            # Certify each shared incidence once, from the earliest
            # NoNewVerts owner — a wall with a single NoNewVerts owner is
            # already certified per-Region by its own boundary check.
            length(owners)>=2 || continue
            tag==owners[1] || continue
            _extrude_nonew_shared_certify(
                [volumes[owner] for owner in owners],s,caller)
        end
    end
    # `Reverse Surface` flips only the published windings — apply it to the
    # stored parts after every certificate, so the projection lookup and
    # the published mesh see the same flipped part.
    for plan in values(volumes),(surface,part) in plan.surfaces
        get(m.meshing.reverse,(2,surface),false) || continue
        plan.surfaces[surface]=_reversed_surface_mesh(part,caller,surface)
    end
    for (surface,part) in surfaces
        get(m.meshing.reverse,(2,surface),false) || continue
        surfaces[surface]=_reversed_surface_mesh(part,caller,surface)
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
                                      size_field::Union{Nothing,AbstractSizeField}=nothing,
                                      _mesh_input=nothing)
    if scope!==nothing
        _extrude_nonew_require_scope(scope,m,caller)
        part=get(scope.surfaces,t,nothing)
        # Scope parts already carry any `Reverse Surface` winding flip —
        # the scope applied it after certification.
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
    if length(owners)>=2
        # The certified shared incidence plans its scope jointly — a
        # welded lateral plans the whole connected component so a middle
        # Region keeps the certificates of its other bound walls; a shared
        # copied cap is emitted once by the top owner's plan. Every owner
        # must record the same certified owner set, which is the complete
        # Region incidence — including any ordinary layered partners.
        adjacent=_extrude_quadtri_regions(m,t,caller)
        for owner in owners
            get(_extrude_nonew_shared_for(m,owner,caller),t,Int[])==adjacent ||
                _extrude_nonew_shared_planner_error(caller,t)
        end
        scope===nothing || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts Surface[$t] is absent from this mesh scope"))
        link=get(m.meshing.extrude_sources,(2,t),nothing)
        if link!==nothing && link[1]==2
            # The stacked Source owner consumes the copied cap the top
            # owner's plan emits — replay that single emission.
            top=0
            for owner in owners
                v_link=get(m.meshing.extrude_sources,(3,owner),nothing)
                v_link!==nothing && abs(v_link[2])==abs(link[2]) &&
                    (top=owner;break)
            end
            top!=0 || _extrude_nonew_shared_planner_error(caller,t)
            part=_extrude_nonew_complete_plan(m,top,caller;
                min_angle_deg=min_angle_deg,
                max_periodic_passes=max_periodic_passes,
                size_field=size_field,_mesh_input=_mesh_input).surfaces[t]
            get(m.meshing.reverse,(2,t),false) &&
                (part=_reversed_surface_mesh(part,caller,t))
            return part
        end
        built=_extrude_nonew_scope(m,caller;
            tags=_extrude_nonew_shared_closure(m,owners,caller),
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field,_mesh_input=_mesh_input)
        # Scope parts already carry any `Reverse Surface` winding flip.
        return built.surfaces[t]
    end
    isempty(owners) && return nothing
    owner=owners[1]
    if scope===nothing
        completed=_extrude_nonew_complete_plan(m,owner,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field,_mesh_input=_mesh_input)
        part=completed.surfaces[t]
        get(m.meshing.reverse,(2,t),false) &&
            (part=_reversed_surface_mesh(part,caller,t))
        return part
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
    # `Reverse Volume` flips only the published cell windings — the stored
    # plan and its certificates keep canonical positive orientation.
    reversed=get(m.meshing.reverse,(3,t),false)
    if scope===nothing
        volume=_extrude_nonew_complete_plan(m,t,caller;
            min_angle_deg=min_angle_deg,max_periodic_passes=max_periodic_passes,
            size_field=size_field).volume
        reversed && (volume=_reversed_volume_mesh(volume,caller,t))
        return volume
    end
    _extrude_nonew_require_scope(scope,m,caller)
    completed=get(scope.volumes,t,nothing)
    completed===nothing && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts Volume[$t] is absent from this mesh scope"))
    volume=completed.volume
    reversed && (volume=_reversed_volume_mesh(volume,caller,t))
    return volume
end
