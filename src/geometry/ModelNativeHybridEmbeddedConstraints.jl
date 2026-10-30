# One deduplicated embedded inventory feeding the shared lower constructor of
# a connected multi-Region component. Upstream `Embedded{}/"X In Volume"` is a
# global, physical-entity constraint set: the same entity declared on several
# Regions contributes its actual lower rows once; each record lands in the
# member Region whose shell strictly contains it; records lying on a shared
# interface Face split its panels and serve both neighbor Regions.
struct _ModelNativeHybridEmbeddedComponentInput
    model::GeoModel
    regions::Vector{Int}
    input::_ModelMeshAuthority
    inventory::Vector{NTuple{2,Int}}
    owners::Vector{Vector{Int}}
    points::Vector{_ModelNativeHybridEmbeddedPointInput}
    curves::Vector{_ModelNativeHybridEmbeddedCurveInput}
    surfaces::Vector{_ModelNativeHybridEmbeddedSurfaceInput}
    records::Vector{NamedTuple}
    face_regions::Dict{Int,Vector{Int}}
    entity_regions::Dict{NTuple{2,Int},Int}
    interface_entities::Vector{NTuple{2,Int}}
    node_regions::Dict{UInt64,Int}
    panel_nodes::Dict{UInt64,Vector{Int}}
    panel_containers::Dict{UInt64,Vector{Int}}
    boundary_records::Vector{UInt64}
    splits::Dict{Int,Vector{NTuple{3,UInt64}}}
    split_faces::Dict{NTuple{3,UInt64},Int}
    segments::Vector{NTuple{2,UInt64}}
    facets::Vector{NTuple{3,UInt64}}
    region_records::Dict{Int,Vector{UInt64}}
end

# A complete actual Region inventory, with one lower witness per declared entity.
# `pool` marks a receipt synthesized by the connected component embedded pool:
# its inventory is the pool-assigned subset rather than the Region's raw
# `Embedded{}` declaration.
struct _ModelNativeHybridEmbeddedConstraintsInput
    model::GeoModel
    region::Int
    input::_ModelMeshAuthority
    inventory::Vector{NTuple{2,Int}}
    points::Vector{_ModelNativeHybridEmbeddedPointInput}
    curves::Vector{_ModelNativeHybridEmbeddedCurveInput}
    surfaces::Vector{_ModelNativeHybridEmbeddedSurfaceInput}
    records::Vector{NamedTuple}
    pool::Union{Nothing,_ModelNativeHybridEmbeddedComponentInput}
end

function _model_native_hybrid_embedded_constraints_input(m,region,input,caller)
    inventory=copy(get(m.embeds,(3,region),NTuple{2,Int}[]))
    haskey(m.volumes,region) && length(inventory)>1 && all(owner->owner[1] in (0,1,2),inventory) ||
        return nothing
    allunique(inventory) || throw(ArgumentError("$caller: hybrid constraint inventory repeats an entity"))
    # Surface witnesses are built first: each embedded Point's sizing receipt
    # must certify every declared sheet entity in the same inventory.
    surfaces=_ModelNativeHybridEmbeddedSurfaceInput[]
    points=_ModelNativeHybridEmbeddedPointInput[]
    curves=_ModelNativeHybridEmbeddedCurveInput[]
    for (dim,tag) in inventory
        if dim==2
            push!(surfaces,_model_native_hybrid_embedded_surface_input(m,region,input,caller;surface=tag))
        end
    end
    for (dim,tag) in inventory
        if dim==0
            push!(points,_model_native_hybrid_embedded_point_input(m,region,input,caller;
                point=tag,surface_witnesses=surfaces))
        elseif dim==1
            push!(curves,_model_native_hybrid_embedded_curve_input(m,region,input,caller;curve=tag))
        end
    end
    records=Dict{UInt64,NamedTuple}()
    function retain(record)
        current=(uid=record.uid,coordinate=record.coordinate,node_owner=record.node_owner)
        prior=get(records,current.uid,nothing)
        prior===nothing || isequal(prior,current) || throw(ArgumentError(
            "$caller: hybrid constraint inventory changes a shared actual node UID/owner/XYZ"))
        records[current.uid]=current
    end
    for point in points;retain(point.record);end
    for curve in curves,record in curve.records;retain(record);end
    for surface in surfaces,record in surface.records;retain(record);end
    return _ModelNativeHybridEmbeddedConstraintsInput(m,region,input,inventory,points,curves,
        surfaces,NamedTuple[records[uid] for uid in sort!(collect(keys(records)))],nothing)
end

function _model_native_hybrid_constraints_source_certificate(receipt,input,uids,coordinates,caller)
    receipt===nothing && return nothing
    if receipt isa _ModelNativeHybridEmbeddedConstraintsInput && receipt.pool!==nothing
        pool=receipt.pool
        pool.model===receipt.model && pool.input===receipt.input &&
            haskey(receipt.model.volumes,receipt.region) &&
            receipt.region in pool.regions || throw(ArgumentError(
            "$caller: hybrid constraints lack their complete component Region inventory"))
        assigned=sort!(NTuple{2,Int}[entity for entity in pool.inventory
            if get(pool.entity_regions,entity,0)==receipt.region])
        sort!(copy(receipt.inventory))==assigned && length(receipt.inventory)>1 ||
            throw(ArgumentError(
            "$caller: hybrid constraints lack their component pool Region inventory"))
        for point in receipt.points
            any(witness->witness===point,pool.points) &&
                get(pool.entity_regions,(0,Int(point.record.owner[2])),0)==receipt.region ||
                throw(ArgumentError("$caller: hybrid Point constraint has a foreign inventory"))
            _model_native_hybrid_point_source_certificate(point,input,uids,coordinates,caller)
        end
        for curve in receipt.curves
            any(witness->witness===curve,pool.curves) &&
                get(pool.entity_regions,(1,curve.curve),0)==receipt.region ||
                throw(ArgumentError("$caller: hybrid Curve constraint has a foreign inventory"))
            _model_native_hybrid_curve_source_certificate(curve,input,uids,coordinates,caller)
        end
        for surface in receipt.surfaces
            any(witness->witness===surface,pool.surfaces) &&
                get(pool.entity_regions,(2,surface.surface),0)==receipt.region ||
                throw(ArgumentError("$caller: hybrid Surface constraint has a foreign inventory"))
            _model_native_hybrid_surface_source_certificate(surface,input,uids,coordinates,caller)
        end
        owners=vcat([(0,Int(point.record.owner[2])) for point in receipt.points],
                    [(1,curve.curve) for curve in receipt.curves],
                    [(2,surface.surface) for surface in receipt.surfaces])
        length(owners)==length(receipt.inventory) && allunique(owners) &&
            Set(owners)==Set(receipt.inventory) || throw(ArgumentError(
            "$caller: hybrid constraints omit or repeat an actual declared entity"))
        fields=_model_native_hybrid_pool_region_fields(pool,receipt.region,caller)
        fresh=fields[4]
        fresh!==nothing && fresh.records==receipt.records || throw(ArgumentError(
            "$caller: hybrid constraint inventory changes its complete actual Source rows"))
        return nothing
    end
    receipt isa _ModelNativeHybridEmbeddedConstraintsInput && receipt.input===input &&
        length(receipt.model.volumes)==1 && haskey(receipt.model.volumes,receipt.region) &&
        get(receipt.model.embeds,(3,receipt.region),NTuple{2,Int}[])==receipt.inventory &&
        length(receipt.inventory)>1 && allunique(receipt.inventory) &&
        all(owner->owner[1] in (0,1,2),receipt.inventory) || throw(ArgumentError(
        "$caller: hybrid constraints lack their complete current Region inventory"))
    owners=vcat([(0,Int(point.record.owner[2])) for point in receipt.points],
                [(1,curve.curve) for curve in receipt.curves],
                [(2,surface.surface) for surface in receipt.surfaces])
    length(owners)==length(receipt.inventory) && allunique(owners) &&
        Set(owners)==Set(receipt.inventory) || throw(ArgumentError(
        "$caller: hybrid constraints omit or repeat an actual declared entity"))
    for point in receipt.points
        point.model===receipt.model && point.region==receipt.region && point.inventory==receipt.inventory ||
            throw(ArgumentError("$caller: hybrid Point constraint has a foreign inventory"))
        _model_native_hybrid_point_source_certificate(point,input,uids,coordinates,caller)
    end
    for curve in receipt.curves
        curve.model===receipt.model && curve.region==receipt.region && curve.inventory==receipt.inventory ||
            throw(ArgumentError("$caller: hybrid Curve constraint has a foreign inventory"))
        _model_native_hybrid_curve_source_certificate(curve,input,uids,coordinates,caller)
    end
    for surface in receipt.surfaces
        surface.model===receipt.model && surface.region==receipt.region &&
            surface.inventory==receipt.inventory || throw(ArgumentError(
            "$caller: hybrid Surface constraint has a foreign inventory"))
        _model_native_hybrid_surface_source_certificate(surface,input,uids,coordinates,caller)
    end
    fresh=_model_native_hybrid_embedded_constraints_input(receipt.model,receipt.region,input,caller)
    fresh!==nothing && fresh.records==receipt.records || throw(ArgumentError(
        "$caller: hybrid constraint inventory changes its complete actual Source rows"))
    return nothing
end

function _model_native_hybrid_constraints_containment_certificate(constructor,caller)
    receipt=constructor.embedded_constraints;receipt===nothing && return nothing
    _model_native_hybrid_constraints_source_certificate(receipt,constructor.input,
        constructor.node_uids,constructor.coordinates,caller)
    for point in receipt.points
        _model_native_hybrid_point_containment_certificate(merge(constructor,(embedded_point=point,)),caller)
    end
    for curve in receipt.curves
        _model_native_hybrid_curve_containment_certificate(merge(constructor,(embedded_curve=curve,)),caller)
    end
    for surface in receipt.surfaces
        _model_native_hybrid_surface_containment_certificate(merge(constructor,(embedded_surface=surface,)),caller)
    end
    return nothing
end

function _model_native_hybrid_constraints_segments(constructor,uids,caller)
    receipt=constructor.embedded_constraints;receipt===nothing && return NTuple{2,Int32}[]
    segments=NTuple{2,Int32}[];seen=Set{NTuple{2,Int32}}()
    for curve in receipt.curves
        for edge in _model_native_hybrid_curve_segments(merge(constructor,(embedded_curve=curve,)),uids,caller)
            key=minmax(edge...)
            # The geometric constraint is one edge even when actual lower
            # cells from different declared Curves reference the same support.
            key in seen && continue
            push!(seen,key);push!(segments,edge)
        end
    end
    return segments
end

# Each declared sheet contributes its complete literal facets and edges like
# the single-Surface route; duplicate cross-entity supports stay protected once.
function _model_native_hybrid_constraints_facets(constructor,uids,caller)
    receipt=constructor.embedded_constraints
    receipt===nothing && return (facets=NTuple{3,Int32}[],edges=NTuple{2,Int32}[])
    facets=NTuple{3,Int32}[];edges=NTuple{2,Int32}[];seen=Set{NTuple{2,Int32}}()
    for surface in receipt.surfaces
        constraints=_model_native_hybrid_surface_constraints(
            merge(constructor,(embedded_surface=surface,)),uids,caller)
        append!(facets,constraints.facets)
        for edge in constraints.edges
            key=minmax(edge...)
            key in seen && continue
            push!(seen,key);push!(edges,edge)
        end
    end
    return (;facets,edges)
end

# Every Surface witness owned by a constructor route, as one flat tuple for
# the global sizing/seed consumers.
function _model_native_hybrid_surface_witnesses(constructor)
    witnesses=_ModelNativeHybridEmbeddedSurfaceInput[]
    constructor.embedded_surface===nothing || push!(witnesses,constructor.embedded_surface)
    receipt=constructor.embedded_constraints
    receipt===nothing || append!(witnesses,receipt.surfaces)
    return Tuple(witnesses)
end

function _model_native_hybrid_embedded_constraints_certificate(mesh,uids,constructor,caller)
    receipt=constructor.embedded_constraints;receipt===nothing && return nothing
    _model_native_hybrid_constraints_source_certificate(receipt,constructor.input,uids,mesh.coords,caller)
    for point in receipt.points
        _model_native_hybrid_embedded_point_certificate(mesh,uids,
            merge(constructor,(embedded_point=point,)),caller)
    end
    for curve in receipt.curves
        _model_native_hybrid_embedded_curve_certificate(mesh,uids,
            merge(constructor,(embedded_curve=curve,)),caller)
    end
    for surface in receipt.surfaces
        _model_native_hybrid_embedded_surface_certificate(mesh,uids,
            merge(constructor,(embedded_surface=surface,)),caller)
    end
    return nothing
end

# ============================================================================
# Connected-component embedded pool
#
# Upstream `Embedded{}/"X In Volume"` declares a global physical constraint
# set: the same entity named on several Regions contributes its actual lower
# rows once; each record lands in the member Region whose shell strictly
# contains it; records lying on a shared interface Face retriangulate that
# Face's panels and serve both neighbor Regions. Tessella keeps that pool
# semantics but certifies a conforming interface: on-panel nodes split their
# containing literal triangles through an exact CDT, on-panel segments must be
# recovered as sub-facet edges, and embedded Surfaces may not overlap or cross
# an interface (upstream emits "overlapping facets" there).
# ============================================================================

_model_native_hybrid_embedded_entity_label(entity)=
    entity[1]==0 ? "embedded Point" : entity[1]==1 ? "embedded Curve" : "embedded Surface"

# Non-throwing strict-interior shell test mirroring `_model_native_hybrid_shell_inside`.
function _model_native_hybrid_pool_shell_inside(point,panels,coordinate)
    bounds=[coordinate(uid) for panel in panels for uid in panel.nodes]
    span=max(1.0,maximum(maximum(b[axis] for b in bounds)-
                         minimum(b[axis] for b in bounds) for axis in 1:3))
    directions=((1.0,0.37139067635410372,0.19518001458970673),
                (0.41421356237309503,1.0,0.27182818284590454),
                (0.23606797749978981,0.57735026918962584,1.0),
                (1.0,1.0,1.0),
                (1.0,-0.70710678118654752,0.61803398874989490),
                (-0.61803398874989490,1.0,-0.41421356237309515))
    for direction in directions
        far=ntuple(axis->point[axis]+direction[axis]*span*4,3)
        crossings=0;degenerate=false
        for panel in panels,triangle in _model_native_hybrid_panel_triangles(panel)
            a,b,c=map(coordinate,triangle)
            relation=_model_native_hybrid_segment_triangle_relation(point,far,a,b,c)
            if relation<0;degenerate=true;break;end
            crossings+=relation
        end
        degenerate && continue
        return isodd(crossings)
    end
    return false
end

# Phase one: deduplicate the member Regions' `Embedded{}` inventory, audit the
# global two-owner bound, and build one lower witness per unique entity. The
# declaring member Region is kept on each witness as provenance; containment is
# classified afterwards against the component constructor.
function _model_native_hybrid_embedded_pool_witnesses(m,regions,input,caller)
    owners=Dict{NTuple{2,Int},Vector{Int}}()
    for tag in regions
        for entity in get(m.embeds,(3,tag),NTuple{2,Int}[])
            entity[1] in (0,1,2) || throw(ArgumentError(
                "$caller: hybrid embedded batch supports Point/Curve/Surface declarations"))
            list=get!(()->Int[],owners,entity);tag in list || push!(list,tag)
        end
    end
    isempty(owners) && return nothing
    # Declarations on Regions outside this connected group still count toward
    # the upstream two-owner bound even though they are not pooled here.
    for ((dim,tag),entities) in m.embeds
        dim==3 || continue
        for entity in entities
            list=get(owners,entity,nothing);list===nothing && continue
            tag in list || push!(list,tag)
        end
    end
    for (entity,list) in owners
        length(list)<=2 || throw(ArgumentError("$caller: " *
            _model_native_hybrid_embedded_entity_label(entity) *
            " spans more than two Region owners"))
    end
    inventory=sort!(collect(keys(owners)))
    surfaces=_ModelNativeHybridEmbeddedSurfaceInput[]
    points=_ModelNativeHybridEmbeddedPointInput[]
    curves=_ModelNativeHybridEmbeddedCurveInput[]
    for entity in inventory
        entity[1]==2 && push!(surfaces,_model_native_hybrid_embedded_surface_input(
            m,first(owners[entity]),input,caller;surface=entity[2]))
    end
    for entity in inventory
        if entity[1]==0
            push!(points,_model_native_hybrid_embedded_point_input(
                m,first(owners[entity]),input,caller;point=entity[2],
                surface_witnesses=surfaces))
        elseif entity[1]==1
            push!(curves,_model_native_hybrid_embedded_curve_input(
                m,first(owners[entity]),input,caller;curve=entity[2]))
        end
    end
    records=Dict{UInt64,NamedTuple}()
    function retain(record)
        prior=get(records,record.uid,nothing)
        prior===nothing || prior.coordinate==record.coordinate || throw(ArgumentError(
            "$caller: hybrid embedded witnesses report inconsistent actual geometry"))
        records[record.uid]=record
    end
    for point in points;retain(point.record);end
    for curve in curves,record in curve.records;retain(record);end
    for surface in surfaces,record in surface.records;retain(record);end
    ordered=[records[uid] for uid in sort!(collect(keys(records)))]
    return (;inventory,owners=[copy(owners[entity]) for entity in inventory],
        points,curves,surfaces,records=ordered,uids=Set{UInt64}(keys(records)))
end

# Retriangulate one whole shared Face around its certified on-panel records —
# the same operation upstream performs by remeshing the surface. The face's
# literal panel vertices are the CDT input points, member embedded nodes join
# as free interior vertices (they are never on the boundary wires: a record on
# a wire edge is rejected earlier because one adjacent panel is not a shared
# Face), certified on-panel segments become internal constraints, and the
# wires are the degree-1 literal edges. The returned sub-facets replace ALL of
# the face's literal panels in both neighbor shells; the certificates below
# prove the outer wire is unchanged so the split stays conforming with every
# adjacent Face's untouched panels and that every constraint survived.
function _model_native_hybrid_pool_split_face(face,tris,members,internals,coordinate,caller)
    !isempty(tris) && all(tri->allunique(collect(tri)),tris) || throw(ArgumentError(
        "$caller: hybrid shared Face[$face] split lacks literal panels"))
    a,b,c=map(coordinate,tris[1])
    u=(b[1]-a[1],b[2]-a[2],b[3]-a[3]);v=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
    normal=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    all(isfinite,normal) && !all(==(0.0),normal) || throw(ArgumentError(
        "$caller: hybrid shared Face[$face] has a degenerate literal triangle"))
    vertices=UInt64[];seen=Set{UInt64}()
    for tri in tris,uid in tri
        uid in seen || (push!(seen,uid);push!(vertices,uid))
    end
    # The 2-D CDT only certifies planar domains: every literal vertex and
    # member record must sit exactly on the face's first-triangle plane.
    all(uid->orient3(a,b,c,coordinate(uid))==0,vertices) || throw(ArgumentError(
        "$caller: embedded interface records on a nonplanar shared Face[$face] are unsupported"))
    member_uids=UInt64[]
    for uid in members
        uid in seen && throw(ArgumentError(
            "$caller: embedded interface node coincides with a shared Face vertex"))
        orient3(a,b,c,coordinate(uid))==0 || throw(ArgumentError(
            "$caller: embedded interface node leaves its shared Face[$face] plane"))
        push!(member_uids,uid)
    end
    drop=argmax(abs.(collect(normal)));axes=filter(!=(drop),(1,2,3))
    project(x)=(x[axes[1]],x[axes[2]])
    uids=vcat(vertices,member_uids)
    points=[project(coordinate(uid)) for uid in uids]
    pos=Dict(uid=>row for (row,uid) in enumerate(uids))
    degree=Dict{NTuple{2,UInt64},Int}()
    for tri in tris,index in 1:3
        key=minmax(tri[index],tri[mod1(index+1,3)])
        degree[key]=get(degree,key,0)+1
    end
    boundary=NTuple{2,Int32}[]
    for (edge,count) in degree
        count==1 || continue
        push!(boundary,(Int32(pos[edge[1]]),Int32(pos[edge[2]])))
    end
    internal=NTuple{2,Int32}[]
    for (x_uid,y_uid) in internals
        haskey(pos,x_uid) && haskey(pos,y_uid) || throw(ArgumentError(
            "$caller: embedded interface segment omits its shared Face endpoint"))
        push!(internal,(Int32(pos[x_uid]),Int32(pos[y_uid])))
    end
    xs=[point[1] for point in points];ys=[point[2] for point in points]
    T=constrained_delaunay(xs,ys,boundary;internal_segments=internal,rng_seed=1)
    mask=classify_interior(T)
    lookup=Dict{Tuple{Float64,Float64},Int}()
    for (k,point) in enumerate(points)
        lookup[(point[1]==0.0 ? 0.0 : point[1],point[2]==0.0 ? 0.0 : point[2])]=k
    end
    subs=NTuple{3,UInt64}[]
    for t in eachindex(T.alive)
        (T.alive[t] && mask[t] && !_is_ghost_tri(T,t)) || continue
        push!(subs,ntuple(slot->begin
            x=T.x[_vert(T,t,slot)];y=T.y[_vert(T,t,slot)]
            uids[lookup[(x==0.0 ? 0.0 : x,y==0.0 ? 0.0 : y)]]
        end,3))
    end
    isempty(subs) && throw(ArgumentError(
        "$caller: hybrid shared Face[$face] retriangulation produced no sub-facets"))
    covered=Set(minmax(sub[index],sub[mod1(index+1,3)]) for sub in subs for index in 1:3)
    for segment in internals
        minmax(segment...) in covered || throw(ArgumentError(
            "$caller: embedded interface segment is not a shared Face sub-facet edge"))
    end
    sub_degree=Dict{NTuple{2,UInt64},Int}()
    for sub in subs,index in 1:3
        key=minmax(sub[index],sub[mod1(index+1,3)])
        sub_degree[key]=get(sub_degree,key,0)+1
    end
    Set(key for (key,count) in sub_degree if count==1)==
        Set(key for (key,count) in degree if count==1) || throw(ArgumentError(
        "$caller: hybrid shared Face[$face] split changes its literal boundary wire"))
    used=Set{UInt64}(uid for sub in subs for uid in sub)
    all(uid->uid in used,uids) || throw(ArgumentError(
        "$caller: hybrid shared Face[$face] split drops a literal or embedded vertex"))
    return subs
end

# Phase two: classify every pool record against the component constructor's
# literal panels, certify each entity's incidence class, and retriangulate the
# shared panels that carry on-panel records. `witnesses` carries the phase-one
# result already audited by the receipt builder.
function _model_native_hybrid_embedded_component_input(m,regions,input,boundaries,
                                                       face_regions,constructor,caller;
                                                       witnesses=nothing)
    witnesses===nothing &&
        (witnesses=_model_native_hybrid_embedded_pool_witnesses(m,regions,input,caller))
    witnesses===nothing && return nothing
    panels=vcat(constructor.triangles,constructor.quads)
    positions=Dict(uid=>row for (row,uid) in enumerate(constructor.node_uids))
    coordinate(uid)=Tuple(constructor.coordinates[:,positions[uid]])
    lit=Tuple{NTuple{3,UInt64},Int}[]
    for panel in panels
        face=Int(panel.owner[2])
        for triangle in _model_native_hybrid_panel_triangles(panel)
            push!(lit,(triangle,face))
        end
    end
    boundary_uids=_model_native_hybrid_shell_boundary(panels)
    bycoord=Dict{NTuple{3,Float64},UInt64}()
    for uid in constructor.node_uids
        c=coordinate(uid);haskey(bycoord,c) || (bycoord[c]=uid)
    end
    region_panels=Dict(tag=>[panel for panel in panels
        if tag in face_regions[Int(panel.owner[2])]] for tag in regions)
    node_regions=Dict{UInt64,Int}();panel_nodes=Dict{UInt64,Vector{Int}}()
    panel_containers=Dict{UInt64,Vector{Int}}();boundary_records=UInt64[]
    boundary_set=Set{UInt64}()
    for record in witnesses.records
        uid=record.uid;point=record.coordinate
        all(isfinite,point) || throw(ArgumentError(
            "$caller: hybrid embedded record carries a non-finite coordinate"))
        twin=get(bycoord,point,UInt64(0))
        twin==0 || twin==uid || throw(ArgumentError(
            "$caller: embedded entity node $uid coincides with a hybrid shell vertex"))
        if uid in boundary_uids
            push!(boundary_records,uid);push!(boundary_set,uid);continue
        end
        # Panel incidence is decided before the strict-interior ray test: a
        # node coplanar with a literal triangle is on the shell even when a
        # degenerate ray direction reports an interior crossing.
        containers=Int[]
        for (index,(triangle,face)) in enumerate(lit)
            a,b,c=map(coordinate,triangle)
            orient3(a,b,c,point)==0 &&
                _model_native_hybrid_coplanar_on_triangle(point,a,b,c) || continue
            push!(containers,index)
        end
        if isempty(containers)
            inside=Int[tag for tag in regions if _model_native_hybrid_pool_shell_inside(
                point,region_panels[tag],coordinate)]
            if length(inside)==1
                node_regions[uid]=only(inside);continue
            end
            length(inside)>1 && throw(ArgumentError(
                "$caller: hybrid embedded node $uid is inside several Region shells"))
            throw(ArgumentError(
                "$caller: embedded entity node $uid lies outside or on the hybrid batch boundary"))
        end
        faces=sort!(unique!(Int[lit[index][2] for index in containers]))
        all(face->length(face_regions[face])==2,faces) || throw(ArgumentError(
            "$caller: embedded entity node $uid lies on the hybrid batch boundary"))
        length(faces)==1 || throw(ArgumentError("$caller: embedded entity node $uid " *
            "lies on a shared Face boundary edge — an ambiguous interface"))
        panel_nodes[uid]=faces;panel_containers[uid]=containers
    end
    entity_records=Dict{NTuple{2,Int},Vector{NamedTuple}}()
    for point in witnesses.points
        entity_records[(0,Int(point.record.owner[2]))]=[point.record]
    end
    for curve in witnesses.curves
        entity_records[(1,curve.curve)]=curve.records
    end
    for surface in witnesses.surfaces
        entity_records[(2,surface.surface)]=surface.records
    end
    entity_regions=Dict{NTuple{2,Int},Int}();interface_entities=NTuple{2,Int}[]
    segments=NTuple{2,UInt64}[];facets=NTuple{3,UInt64}[]
    segment_set=Set{NTuple{2,UInt64}}();facet_set=Set{NTuple{3,UInt64}}()
    region_records=Dict{Int,Vector{UInt64}}()
    register(uid,tag)=(push!(get!(()->UInt64[],region_records,tag),uid))
    face_internals=Dict{Int,Vector{NTuple{2,UInt64}}}()
    # Every Face whose literal closure contains this record: an on-panel
    # member maps to its certified single shared Face, a literal vertex to
    # every incident panel owner. A certified on-panel segment joins the
    # retriangulation of the unique shared Face reaching both endpoints.
    function lit_faces_of(uid)
        haskey(panel_nodes,uid) && return copy(panel_nodes[uid])
        faces=Int[]
        for (triangle,face) in lit
            uid in triangle && push!(faces,face)
        end
        return sort!(unique!(faces))
    end
    for entity in witnesses.inventory
        records=entity_records[entity]
        interior=sort!(unique!(Int[node_regions[record.uid] for record in records
            if haskey(node_regions,record.uid)]))
        on_panel=any(record->haskey(panel_nodes,record.uid),records)
        all_interior=length(interior)==1 &&
            all(record->get(node_regions,record.uid,0)==only(interior),records)
        if all_interior
            region=only(interior);entity_regions[entity]=region
            for record in records;register(record.uid,region);end
        elseif entity[1]==2
            on_panel && throw(ArgumentError("$caller: embedded Surface overlaps or " *
                "touches a shared Region boundary — overlapping facets are unsupported"))
            length(interior)>1 && throw(ArgumentError(
                "$caller: embedded Surface crosses a Region boundary"))
            # All-boundary sheets are the upstream deduplicated no-op only when
            # every literal Tri3 coincides with a shared panel.
            surface=only(Iterators.filter(
                witness->witness.surface==entity[2],witnesses.surfaces))
            for triangle in surface.triangles
                key=Tuple(sort!(collect(triangle.nodes)))
                any(index->Tuple(sort!(collect(lit[index][1])))==key &&
                        length(face_regions[lit[index][2]])==2,eachindex(lit)) || throw(
                    ArgumentError("$caller: embedded Surface overlaps a shared Region " *
                        "boundary panel — overlapping facets are unsupported"))
            end
        else
            # Curve or Point spanning shared interfaces: every record is either
            # strictly inside a member Region, on a shared panel, or an actual
            # shell vertex; every segment must be certifiable in one class.
            push!(interface_entities,entity)
            for record in records
                if haskey(node_regions,record.uid)
                    register(record.uid,node_regions[record.uid])
                elseif haskey(panel_nodes,record.uid)
                    for face in panel_nodes[record.uid],tag in face_regions[face]
                        register(record.uid,tag)
                    end
                end
            end
        end
        if entity[1]==2 && all_interior
            surface=only(Iterators.filter(
                witness->witness.surface==entity[2],witnesses.surfaces))
            for triangle in surface.triangles
                key=Tuple(sort!(collect(triangle.nodes)))
                key in facet_set && continue
                push!(facet_set,key);push!(facets,triangle.nodes)
            end
        end
        if entity[1]==1
            # A sheet's boundary wire rides on its protected facet triples;
            # only Curve cell edges are classified into interior/on-panel
            # segments here.
            cells=NTuple{2,UInt64}[minmax(cell.nodes...) for cell in
                only(Iterators.filter(witness->witness.curve==entity[2],
                    witnesses.curves)).cells]
            for seg in cells
                a_uid,b_uid=seg
                ra=get(node_regions,a_uid,0);rb=get(node_regions,b_uid,0)
                pa=haskey(panel_nodes,a_uid);pb=haskey(panel_nodes,b_uid)
                ba=a_uid in boundary_set;bb=b_uid in boundary_set
                if all_interior
                    # Interior entity: every literal cell edge is an interior
                    # protected segment and must miss every panel.
                    p,q=coordinate(a_uid),coordinate(b_uid)
                    for (triangle,face) in lit
                        ta,tb,tc=map(coordinate,triangle)
                        _model_native_hybrid_segment_triangle_relation(p,q,ta,tb,tc)==0 ||
                            throw(ArgumentError("$caller: " *
                                _model_native_hybrid_embedded_entity_label(entity) *
                                " is outside or on the hybrid shell boundary"))
                    end
                    key=minmax(a_uid,b_uid)
                    key in segment_set ||
                        (push!(segment_set,key);push!(segments,(a_uid,b_uid)))
                    continue
                end
                entity[1]==2 && continue   # sheet edges already handled above
                if ra!=0 && rb!=0
                    ra==rb || throw(ArgumentError("$caller: embedded Curve crosses " *
                        "a Region boundary at a non-shared-boundary point"))
                    # Both endpoints strictly inside the same Region: an
                    # interior protected segment that must miss every panel.
                    p,q=coordinate(a_uid),coordinate(b_uid)
                    for (triangle,face) in lit
                        ta,tb,tc=map(coordinate,triangle)
                        _model_native_hybrid_segment_triangle_relation(p,q,ta,tb,tc)==0 ||
                            throw(ArgumentError("$caller: " *
                                _model_native_hybrid_embedded_entity_label(entity) *
                                " is outside or on the hybrid shell boundary"))
                    end
                    key=minmax(a_uid,b_uid)
                    key in segment_set ||
                        (push!(segment_set,key);push!(segments,(a_uid,b_uid)))
                elseif ra!=0 || rb!=0
                    region=ra!=0 ? ra : rb;other=ra!=0 ? b_uid : a_uid
                    if haskey(panel_nodes,other)
                        any(face->region in face_regions[face],panel_nodes[other]) ||
                            throw(ArgumentError("$caller: embedded Curve ends on a shared " *
                                "Face that does not bound its Region"))
                        allowed=Set{Int}(panel_containers[other])
                    else
                        other in boundary_set || throw(ArgumentError(
                            "$caller: embedded Curve has an unclassified segment endpoint"))
                        any(panel->other in panel.nodes,region_panels[region]) ||
                            throw(ArgumentError("$caller: embedded Curve ends on a shell " *
                                "vertex outside its Region"))
                        allowed=Set(index for (index,(triangle,face)) in enumerate(lit)
                            if other in triangle)
                    end
                    p,q=coordinate(a_uid),coordinate(b_uid)
                    for (index,(triangle,face)) in enumerate(lit)
                        ta,tb,tc=map(coordinate,triangle)
                        relation=_model_native_hybrid_segment_triangle_relation(p,q,ta,tb,tc)
                        (relation==0 || (relation==-1 && index in allowed)) ||
                            throw(ArgumentError("$caller: embedded Curve crosses a Region " *
                                "boundary at a non-shared-boundary point"))
                    end
                    key=minmax(a_uid,b_uid)
                    key in segment_set ||
                        (push!(segment_set,key);push!(segments,(a_uid,b_uid)))
                elseif (pa||ba) && (pb||bb)
                    key=minmax(a_uid,b_uid)
                    common=intersect(lit_faces_of(a_uid),lit_faces_of(b_uid))
                    shared=[face for face in common if length(face_regions[face])==2]
                    length(shared)==1 || throw(ArgumentError("$caller: embedded Curve " *
                        (isempty(shared) ?
                            "crosses a shared Face panel edge at a non-vertex point" :
                            "segment lies on several shared Faces — an ambiguous interface")))
                    push!(get!(()->NTuple{2,UInt64}[],face_internals,only(shared)),key)
                else
                    throw(ArgumentError("$caller: embedded Curve has an unclassified segment"))
                end
            end
        end
    end
    # Retriangulate every shared Face carrying on-panel records or certified
    # on-panel segments; the sub-facets replace that Face's entire literal
    # panel set in both neighbor shells.
    face_members=Dict{Int,Vector{UInt64}}()
    for (uid,faces) in panel_nodes,face in faces
        push!(get!(()->UInt64[],face_members,face),uid)
    end
    face_tris=Dict{Int,Vector{NTuple{3,UInt64}}}()
    for (triangle,face) in lit
        push!(get!(()->NTuple{3,UInt64}[],face_tris,face),triangle)
    end
    splits=Dict{Int,Vector{NTuple{3,UInt64}}}()
    split_faces=Dict{NTuple{3,UInt64},Int}()
    for face in sort!(unique!(vcat(collect(keys(face_members)),
                                   collect(keys(face_internals)))))
        tris=face_tris[face]
        subs=_model_native_hybrid_pool_split_face(face,tris,
            get(face_members,face,UInt64[]),get(face_internals,face,NTuple{2,UInt64}[]),
            coordinate,caller)
        splits[face]=subs
        for triangle in tris
            split_faces[Tuple(sort!(collect(triangle)))]=face
        end
    end
    for tag in keys(region_records);sort!(unique!(region_records[tag]));end
    return _ModelNativeHybridEmbeddedComponentInput(m,copy(regions),input,
        witnesses.inventory,witnesses.owners,witnesses.points,witnesses.curves,
        witnesses.surfaces,witnesses.records,face_regions,entity_regions,
        interface_entities,node_regions,panel_nodes,panel_containers,
        boundary_records,splits,split_faces,segments,facets,region_records)
end

# Semantic freshness comparison for two independently built pools over the
# same authority: mutable model/input members bind by identity while every
# classified inventory, record, split and constraint datum compares by
# content. Witness inputs reduce to their entity tag plus actual record/cell
# rows so a re-derived receipt is not rejected for object identity alone.
function _model_native_hybrid_pool_equal(a,b)
    a===b && return true
    (a===nothing)!=(b===nothing) && return false
    a isa _ModelNativeHybridEmbeddedComponentInput &&
        b isa _ModelNativeHybridEmbeddedComponentInput &&
        a.model===b.model && a.input===b.input || return false
    a.regions==b.regions && a.inventory==b.inventory && a.owners==b.owners &&
        a.records==b.records && a.face_regions==b.face_regions &&
        a.entity_regions==b.entity_regions &&
        a.interface_entities==b.interface_entities &&
        a.node_regions==b.node_regions && a.panel_nodes==b.panel_nodes &&
        a.panel_containers==b.panel_containers &&
        a.boundary_records==b.boundary_records &&
        a.splits==b.splits && a.split_faces==b.split_faces &&
        a.segments==b.segments && a.facets==b.facets &&
        a.region_records==b.region_records || return false
    length(a.points)==length(b.points) && length(a.curves)==length(b.curves) &&
        length(a.surfaces)==length(b.surfaces) || return false
    sort!([(0,Int(point.record.owner[2])) for point in a.points])==
        sort!([(0,Int(point.record.owner[2])) for point in b.points]) &&
        sort!([(1,curve.curve) for curve in a.curves])==
        sort!([(1,curve.curve) for curve in b.curves]) &&
        sort!([(2,surface.surface) for surface in a.surfaces])==
        sort!([(2,surface.surface) for surface in b.surfaces]) || return false
    for (x,y) in zip(sort!(copy(a.points),by=point->Int(point.record.owner[2])),
                     sort!(copy(b.points),by=point->Int(point.record.owner[2])))
        x.record==y.record || return false
    end
    for (x,y) in zip(sort!(copy(a.curves),by=curve->curve.curve),
                     sort!(copy(b.curves),by=curve->curve.curve))
        x.endpoints==y.endpoints && x.records==y.records && x.cells==y.cells ||
            return false
    end
    for (x,y) in zip(sort!(copy(a.surfaces),by=surface->surface.surface),
                     sort!(copy(b.surfaces),by=surface->surface.surface))
        x.loops==y.loops && x.endpoints==y.endpoints && x.triangles==y.triangles &&
            x.records==y.records && x.lower_cells==y.lower_cells || return false
    end
    return true
end

# Per-Region embedded fields of the pool constructor: strictly-contained
# entities keep their full lower witness; interface entities are joint
# interface constraints certified once at the component level.
function _model_native_hybrid_pool_region_fields(pool,tag,caller)
    tag in pool.regions || throw(ArgumentError(
        "$caller: hybrid embedded pool lacks its selected Region"))
    assigned=NTuple{2,Int}[entity for entity in pool.inventory
        if get(pool.entity_regions,entity,0)==tag]
    isempty(assigned) && return (nothing,nothing,nothing,nothing)
    points=[point for point in pool.points if (0,Int(point.record.owner[2])) in assigned]
    curves=[curve for curve in pool.curves if (1,curve.curve) in assigned]
    surfaces=[surface for surface in pool.surfaces if (2,surface.surface) in assigned]
    if length(assigned)==1
        return (isempty(points) ? nothing : only(points),
                isempty(curves) ? nothing : only(curves),
                isempty(surfaces) ? nothing : only(surfaces),nothing)
    end
    records=Dict{UInt64,NamedTuple}()
    function retain(record)
        prior=get(records,record.uid,nothing)
        prior===nothing || prior.coordinate==record.coordinate || throw(ArgumentError(
            "$caller: hybrid embedded pool reports inconsistent actual geometry"))
        records[record.uid]=record
    end
    for point in points;retain(point.record);end
    for curve in curves,record in curve.records;retain(record);end
    for surface in surfaces,record in surface.records;retain(record);end
    receipt=_ModelNativeHybridEmbeddedConstraintsInput(pool.model,tag,pool.input,
        collect(assigned),points,curves,surfaces,
        NamedTuple[records[uid] for uid in sort!(collect(keys(records)))],pool)
    return (nothing,nothing,nothing,receipt)
end

# Kernel-row interior constraints: certified interior segments (including
# interface-touching chains) and interior sheet facets. On-panel segments are
# already recovered as shared sub-facet edges and are not repeated here.
function _model_native_hybrid_pool_kernel_constraints(pool,node_uids,order,caller)
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(node_uids))
    inverse=invperm(order)
    segments=NTuple{2,Int32}[]
    for (a,b) in pool.segments
        haskey(positions,a) && haskey(positions,b) || throw(ArgumentError(
            "$caller: hybrid embedded pool segment omits an actual Source UID"))
        push!(segments,(Int32(inverse[positions[a]]),Int32(inverse[positions[b]])))
    end
    facets=NTuple{3,Int32}[]
    for face in pool.facets
        all(uid->haskey(positions,uid),face) || throw(ArgumentError(
            "$caller: hybrid embedded pool facet omits an actual Source UID"))
        push!(facets,ntuple(slot->Int32(inverse[positions[face[slot]]]),3))
    end
    return segments,facets
end

# Replace every literal facet of a split shared Face by that Face's certified
# sub-facets in the kernel-row recovery facet list, preserving deterministic
# order. `pool.split_faces` maps each literal panel's sorted UID key to its
# owning Face so every literal facet of a split Face is dropped exactly once.
function _model_native_hybrid_pool_split_facets(facets,pool,plan,order,caller)
    isempty(pool.splits) && return facets
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(plan.node_uids))
    inverse=invperm(order)
    literal_faces=Dict{NTuple{3,Int32},Int}()
    for (key,face) in pool.split_faces
        all(uid->haskey(positions,uid),key) || throw(ArgumentError(
            "$caller: hybrid shared Face split omits an actual Source UID"))
        literal_faces[tuple(sort!(Int32[inverse[positions[uid]] for uid in key])...)]=face
    end
    kernel_subs=Dict{Int,Vector{NTuple{3,Int32}}}()
    for (face,subs) in pool.splits
        kernel_subs[face]=[ntuple(slot->begin
            haskey(positions,sub[slot]) || throw(ArgumentError(
                "$caller: hybrid shared Face sub-facet omits an actual Source UID"))
            Int32(inverse[positions[sub[slot]]])
        end,3) for sub in subs]
    end
    result=NTuple{3,Int32}[];emitted=Set{Int}()
    for facet in facets
        face=get(literal_faces,Tuple(sort!(collect(facet))),0)
        if face==0
            push!(result,facet)
        elseif !(face in emitted)
            push!(emitted,face);append!(result,kernel_subs[face])
        end
    end
    length(emitted)==length(pool.splits) || throw(ArgumentError(
        "$caller: hybrid recovery omits a certified shared Face split"))
    return result
end

# Joint component-level certificate over the published per-Region products.
# Every pool record must live in the product of each Region it serves; every
# certified interior segment is a cell edge in its containing product; every
# on-panel node is a shared sub-facet vertex in every neighbor product; a
# certified on-panel segment is a shared sub-facet edge on the interface.
function _model_native_hybrid_embedded_pool_certificate(products,pool,caller)
    uids=Dict(tag=>Set(product.boundary_node_tags) for (tag,product) in products)
    for record in pool.records
        if haskey(pool.node_regions,record.uid)
            tag=pool.node_regions[record.uid]
            record.uid in uids[tag] || throw(ArgumentError(
                "$caller: embedded entity node is missing from its Region product"))
        elseif haskey(pool.panel_nodes,record.uid)
            for face in pool.panel_nodes[record.uid],tag in pool.face_regions[face]
                record.uid in uids[tag] || throw(ArgumentError("$caller: embedded " *
                    "interface node is missing from a neighbor Region product"))
            end
        else
            any(tag->record.uid in uids[tag],pool.regions) || throw(ArgumentError(
                "$caller: embedded boundary node is missing from the Region products"))
        end
    end
    # Per-product edge sets and interface sub-facet edge sets.
    edges=Dict{Int,Set{NTuple{2,Int32}}}();positions=Dict{Int,Dict{UInt64,Int32}}()
    boundary_edges=Dict{Int,Set{NTuple{2,Int32}}}()
    for (tag,product) in products
        mesh=product.mesh
        positions[tag]=Dict(uid=>Int32(row) for (row,uid) in enumerate(product.boundary_node_tags))
        edges[tag]=_volume_edge_set(mesh)
        counts=Dict{Tuple,Int}()
        for block in mesh.blocks,column in axes(block.nodes,2)
            for face in _VOLUME_CELL_FACES[Int(block.msh)]
                key=Tuple(sort!(Int32[block.nodes[slot,column] for slot in face]))
                counts[key]=get(counts,key,0)+1
            end
        end
        boundary_edges[tag]=Set(minmax(face[index],face[mod1(index+1,length(face))])
            for (face,count) in counts if count==1 for index in 1:length(face))
    end
    sub_edge_uids=Set{NTuple{2,UInt64}}()
    for subs in values(pool.splits),sub in subs,index in 1:3
        push!(sub_edge_uids,minmax(sub[index],sub[mod1(index+1,3)]))
    end
    for entity in pool.interface_entities
        label=_model_native_hybrid_embedded_entity_label(entity)
        # Every interior-side segment of an interface entity is a protected
        # edge in the containing Region's product.
        entity_segments=NTuple{2,UInt64}[]
        if entity[1]==1
            curve=only(Iterators.filter(witness->witness.curve==entity[2],pool.curves))
            append!(entity_segments,(minmax(cell.nodes...) for cell in curve.cells))
        end
        for (a_uid,b_uid) in entity_segments
            minmax(a_uid,b_uid) in sub_edge_uids && continue
            ra=get(pool.node_regions,a_uid,0);rb=get(pool.node_regions,b_uid,0)
            region=ra!=0 ? ra : rb
            region==0 && continue
            pa=get(positions[region],a_uid,0);pb=get(positions[region],b_uid,0)
            (pa!=0&&pb!=0) || continue
            minmax(pa,pb) in edges[region] || throw(ArgumentError(
                "$caller: $label segment is not a Region product cell edge"))
        end
        # An on-panel segment is a sub-facet edge on each neighbor's interface.
        for (a_uid,b_uid) in entity_segments
            key=minmax(a_uid,b_uid);key in sub_edge_uids || continue
            pa=get(pool.node_regions,a_uid,0)
            pa!=0 && continue
            faces=unique!(sort!(vcat(
                (get(pool.panel_nodes,a_uid,Int[])),
                (get(pool.panel_nodes,b_uid,Int[])))))
            isempty(faces) && continue
            for face in faces,tag in pool.face_regions[face]
                xa=get(positions[tag],a_uid,0);xb=get(positions[tag],b_uid,0)
                (xa!=0&&xb!=0) || throw(ArgumentError(
                    "$caller: $label interface segment endpoint is missing from a " *
                    "neighbor Region product"))
                minmax(xa,xb) in boundary_edges[tag] || throw(ArgumentError(
                    "$caller: $label interface segment is not a shared sub-facet edge"))
            end
        end
    end
    return nothing
end
