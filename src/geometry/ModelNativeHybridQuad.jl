# Hybrid recovery first creates one actual MVertex per native Quad vector
# entry. Its four temporary facet triangles name that birth explicitly; they
# are not replacement lower Surface elements. Caller-owned Surface receipts
# must attest the Quad UID order and the current ordered corner lists.
struct _ModelNativeHybridQuadBirth
    owner::Tuple{Int,Int32}
    element_uid::UInt64
    base_uids::NTuple{4,UInt64}
    base_rows::NTuple{4,Int32}
    apex_uid::UInt64
    apex_row::Int32
    point::NTuple{3,Float64}
end

struct _ModelNativeHybridQuadPlan
    source_uids::Vector{UInt64}
    source_coordinates::Matrix{Float64}
    node_uids::Vector{UInt64}
    coordinates::Matrix{Float64}
    fan_facets::Matrix{Int32}
    births::Vector{_ModelNativeHybridQuadBirth}
    initial_clock::UInt64
    node_high_water::UInt64
end

struct _ModelNativeHybridPyramidProduction
    nodes::Matrix{Int32}
    element_tags::Vector{UInt64}
    element_high_water::UInt64
    quad_uids::Vector{UInt64}
end

function _model_native_hybrid_quad_plan(source_uids,source_coordinates,quads,
                                      node_clock::UInt64,caller;normal_factor::Float64=0.0)
    isfinite(normal_factor) && normal_factor>=0.0 || throw(ArgumentError(
        "$caller: hybrid Quad normal-apex factor must be finite and nonnegative"))
    size(source_coordinates,1)==3 && size(source_coordinates,2)==length(source_uids) &&
        all(!iszero,source_uids) && allunique(source_uids) && all(isfinite,source_coordinates) &&
        node_clock>=maximum(source_uids;init=UInt64(0)) || throw(ArgumentError(
        "$caller: hybrid Quad input lacks unique current UID/XYZ rows and allocator clock"))
    length(source_uids)<=typemax(Int32)-length(quads) &&
        length(quads)<=fld(typemax(Int32),4) || throw(ArgumentError(
        "$caller: hybrid Quad constructor exceeds native row/facet bounds"))
    last=Base.checked_add(node_clock,UInt64(length(quads)))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_uids))
    cell_uids=Set{UInt64}();supports=Set{NTuple{4,UInt64}}()
    validated=Tuple{Tuple{Int,Int32},UInt64,NTuple{4,UInt64},NTuple{4,Int32}}[]
    for quad in quads
        quad.owner[1]==2 && quad.owner[2]>0 && quad.uid>0 && !(quad.uid in cell_uids) &&
            length(quad.nodes)==4 && allunique(quad.nodes) && all(uid->haskey(positions,uid),quad.nodes) ||
            throw(ArgumentError("$caller: hybrid Quad has foreign, duplicate or missing owned constructor incidence"))
        ids=NTuple{4,UInt64}(Tuple(quad.nodes))
        support=Tuple(sort!(collect(ids)))
        !(support in supports) || throw(ArgumentError(
            "$caller: hybrid Quad constructor aliases another Surface cell"))
        push!(cell_uids,quad.uid);push!(supports,support)
        push!(validated,((2,Int32(quad.owner[2])),UInt64(quad.uid),ids,
            ntuple(slot->positions[ids[slot]],4)))
    end
    coordinates=Matrix{Float64}(undef,3,length(source_uids)+length(quads))
    coordinates[:,1:length(source_uids)].=source_coordinates
    tags=UInt64.(source_uids)
    facets=Matrix{Int32}(undef,3,4length(quads))
    births=_ModelNativeHybridQuadBirth[]
    for (ordinal,(owner,cell,ids,base)) in enumerate(validated)
        apex_uid=node_clock+UInt64(ordinal)
        apex=Int32(length(source_uids)+ordinal)
        # Native SPoint3 uses this exact left-to-right quarter-sum operation.
        point=ntuple(axis->((source_coordinates[axis,base[1]]+source_coordinates[axis,base[2]]+
            source_coordinates[axis,base[3]]+source_coordinates[axis,base[4]])*0.25),3)
        if normal_factor>0.0
            a,b,c,d=map(row->Tuple(source_coordinates[:,row]),base)
            orient3(a,b,c,d)==0 || throw(ArgumentError(
                "$caller: hybrid normal-apex constructor requires planar Quads"))
            u=ntuple(axis->b[axis]-a[axis],3);v=ntuple(axis->c[axis]-a[axis],3)
            # Native normal3points normalizes this first-three-corner cross.
            normal=(u[2]*v[3]-u[3]*v[2],-u[1]*v[3]+u[3]*v[1],u[1]*v[2]-u[2]*v[1])
            magnitude=sqrt(normal[1]*normal[1]+normal[2]*normal[2]+normal[3]*normal[3])
            isfinite(magnitude) && magnitude>0.0 || throw(ArgumentError(
                "$caller: hybrid normal-apex Quad has a degenerate normal"))
            inverse=1.0/magnitude;normal=ntuple(axis->normal[axis]*inverse,3)
            diagonal=ntuple(axis->a[axis]-c[axis],3)
            diameter=sqrt(diagonal[1]*diagonal[1]+diagonal[2]*diagonal[2]+diagonal[3]*diagonal[3])
            offset=normal_factor*diameter
            point=ntuple(axis->point[axis]+offset*normal[axis],3)
        end
        all(isfinite,point) || throw(ArgumentError("$caller: hybrid Quad centroid is not finite"))
        coordinates[:,apex].=point;push!(tags,apex_uid)
        for side in 1:4
            facets[:,4(ordinal-1)+side].=(base[side],base[mod1(side+1,4)],apex)
        end
        push!(births,_ModelNativeHybridQuadBirth(owner,cell,ids,base,apex_uid,apex,point))
    end
    return _ModelNativeHybridQuadPlan(UInt64.(source_uids),copy(source_coordinates),tags,
        coordinates,facets,births,node_clock,last)
end

# A negative option consumes actual Face winding. Restrict admission to a
# convex shell whose literal Quad normals place all recorded births strictly
# inside; do not repair winding or clamp an exterior apex to a different point.
function _model_native_hybrid_normal_apex_certificate(plan,constructor,caller)
    plan.source_uids==constructor.node_uids &&
        isequal(plan.source_coordinates,constructor.coordinates) || throw(ArgumentError(
        "$caller: hybrid normal-apex constructor lost its fixed actual Source"))
    panels=vcat(constructor.triangles,constructor.quads)
    positions=Dict(uid=>row for (row,uid) in enumerate(constructor.node_uids))
    boundary=Set(uid for panel in panels for uid in panel.nodes)
    coordinate(uid)=Tuple(constructor.coordinates[:,positions[uid]])
    allunique(Tuple(plan.coordinates[:,row]) for row in axes(plan.coordinates,2)) ||
        throw(ArgumentError("$caller: hybrid normal-apex birth aliases another constructor vertex"))
    for panel in panels
        triangles=length(panel.nodes)==3 ? (panel.nodes,) :
            ((panel.nodes[1],panel.nodes[2],panel.nodes[3]),
             (panel.nodes[1],panel.nodes[3],panel.nodes[4]))
        for triangle in triangles
            a,b,c=map(coordinate,triangle)
            sides=[orient3(a,b,c,coordinate(uid)) for uid in boundary]
            all(isfinite,sides) || throw(ArgumentError(
                "$caller: hybrid normal-apex shell halfspaces are not finite"))
            signs=Set(sign(side) for side in sides if side!=0)
            length(signs)==1 || throw(ArgumentError(
                "$caller: hybrid normal-apex constructor requires a certified convex shell"))
            for birth in plan.births
                side=orient3(a,b,c,birth.point)
                isfinite(side) && side!=0 && sign(side)==only(signs) || throw(ArgumentError(
                    "$caller: hybrid normal-apex birth is outside or on the shell boundary"))
            end
        end
    end
    return nothing
end

function _model_native_hybrid_normal_apex_fixed_certificate(plan,mesh,uids,caller)
    positions=Dict(uid=>row for (row,uid) in enumerate(uids))
    for birth in plan.births
        row=get(positions,birth.apex_uid,0)
        row>0 && isequal(Tuple(mesh.coords[:,row]),birth.point) || throw(ArgumentError(
            "$caller: hybrid negative-factor apex lost its actual birth UID/XYZ"))
    end
    return nothing
end

# This final constructor consumes already-bound, relocated apex rows. The
# Source coordinate prefix remains immutable. Native MPyramid::reverse swaps
# the first/third corner; the existing factory certifies both base diagonals.
function _model_native_hybrid_pyramid_nodes(plan,coordinates,caller;
                                            birth_order=eachindex(plan.births))
    size(coordinates,1)==3 && size(coordinates,2)>=length(plan.node_uids) &&
        all(isfinite,coordinates) && isequal(coordinates[:,1:length(plan.source_uids)],
            plan.source_coordinates) || throw(ArgumentError(
        "$caller: hybrid Pyramid lost its fixed actual Source UID coordinate prefix"))
    length(birth_order)==length(plan.births) && allunique(birth_order) &&
        all(ordinal->1<=ordinal<=length(plan.births),birth_order) || throw(ArgumentError(
        "$caller: hybrid Pyramid vector omits or aliases an owned Quad birth"))
    nodes=Matrix{Int32}(undef,5,length(plan.births))
    for (column,ordinal) in enumerate(birth_order)
        birth=plan.births[ordinal]
        a,b,c,d=birth.base_rows;apex=birth.apex_row
        sign=orient3(ntuple(axis->coordinates[axis,a],3),
            ntuple(axis->coordinates[axis,b],3),ntuple(axis->coordinates[axis,d],3),
            ntuple(axis->coordinates[axis,apex],3))
        sign!=0 || throw(ArgumentError("$caller: hybrid Pyramid remains flat after native relocation"))
        # Give the existing factory its outward base; its emitted reversal is
        # the precise native current-Volume corner tuple.
        corners=sign<0 ? (c,b,a,d) : (a,b,c,d)
        _write_pyramid!(nodes,column,coordinates,corners[1],corners[2],corners[3],corners[4],apex,caller)
    end
    return nodes
end

# These MElements are born after region insertion and before the later Tet
# optimizer. A subsequent Tet swap starts at this high-water mark; existing
# Tet UIDs and the immutable Pyramid vector keep their own event bindings.
function _model_native_hybrid_pyramid_production(plan,coordinates,element_clock::UInt64,caller;
                                                birth_order=eachindex(plan.births))
    last=Base.checked_add(element_clock,UInt64(length(plan.births)))
    nodes=_model_native_hybrid_pyramid_nodes(plan,coordinates,caller;birth_order)
    labels=UInt64[element_clock+UInt64(column) for column in eachindex(birth_order)]
    quads=UInt64[plan.births[ordinal].element_uid for ordinal in birth_order]
    return _ModelNativeHybridPyramidProduction(nodes,labels,last,quads)
end

# Boundary recovery uses GEntityPtrLessThan's face tag order; Pyramid
# construction later uses the restored Region GFace vector. The constructor
# receipt supplies each face's cell UID order in both traversals;
# winding comes from the current authoritative cell with that exact UID.
# Neither current public label order nor coincident coordinates select cells.
function _model_native_hybrid_constructor_input(authority,boundaries,node_uids,
                                                coordinates,caller;embedded_point=nothing,embedded_curve=nothing)
    authority isa _ModelMeshAuthority && authority.mesh!==nothing &&
        authority.mesh.entity_data!==nothing || throw(ArgumentError(
        "$caller: hybrid recovery requires current actual mesh authority"))
    length(boundaries)>0 && allunique(abs.(boundaries)) &&
        size(coordinates)==(3,length(node_uids)) && allunique(node_uids) &&
        all(!iszero,node_uids) || throw(ArgumentError(
        "$caller: hybrid recovery lacks an owned region/unique input row vector"))
    actual=authority.mesh;data=actual.entity_data
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(node_uids))
    actual_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    length(actual_rows)==nnodes(actual) || throw(ArgumentError(
        "$caller: hybrid authority repeats an actual node UID"))
    for (row,uid) in enumerate(node_uids)
        current=get(actual_rows,uid,Int32(0))
        current>0 && all(axis->isequal(coordinates[axis,row],actual.coords[axis,current]),1:3) ||
            throw(ArgumentError("$caller: hybrid input changes or omits actual node UID $uid"))
    end
    _model_native_hybrid_point_source_certificate(embedded_point,authority,node_uids,coordinates,caller)
    _model_native_hybrid_curve_source_certificate(embedded_curve,authority,node_uids,coordinates,caller)
    quads=NamedTuple[];triangles=NamedTuple[];cell_uids=Set{UInt64}()
    for surface in sort!(abs.(collect(boundaries)))
        signed_surface=surface
        owner=(2,Int32(abs(signed_surface)))
        haskey(authority.products,owner) || throw(ArgumentError(
            "$caller: hybrid region boundary lacks its actual Surface product $owner"))
        product=_model_mesh_product(authority,2,Int(owner[2]);caller)
        constructor=get(authority.surface_constructor_inputs,owner,nothing)
        constructor!==nothing || throw(ArgumentError(
            "$caller: hybrid boundary lacks its owned Surface constructor vector"))
        _model_surface_constructor_binding(product,constructor,caller)
        current=Dict{UInt64,Tuple{Int,NTuple{4,UInt64}}}()
        pdata=product.entity_data
        for (block,labels,owners) in zip(product.blocks,pdata.external_element_tags,pdata.block_entities)
            block.msh in (2,3) && length(labels)==length(owners)==size(block.nodes,2) ||
                throw(ArgumentError("$caller: hybrid boundary requires actual linear Surface cells"))
            width=block.msh==2 ? 3 : 4
            size(block.nodes,1)==width || throw(ArgumentError(
                "$caller: hybrid Surface cell width differs from its declared family"))
            for (column,uid) in enumerate(labels)
                owners[column]==owner[2] && uid>0 && !(uid in cell_uids) &&
                    !haskey(current,uid) || throw(ArgumentError(
                    "$caller: hybrid Surface cells have duplicate or foreign ownership/UID incidence"))
                ids=UInt64[pdata.external_node_tags[row] for row in block.nodes[:,column]]
                allunique(ids) && all(id->haskey(positions,id),ids) || throw(ArgumentError(
                    "$caller: hybrid boundary cell omits a current constructor node UID"))
                current[uid]=(Int(block.msh),ntuple(slot->slot<=width ? ids[slot] : UInt64(0),4))
            end
        end
        cdata=constructor.entity_data;seen=Set{UInt64}()
        # Actual GFace vectors are separate by family even if a carrier's
        # blocks happen to interleave their presentation.
        for family in (2,3)
            for (block,labels) in zip(constructor.blocks,cdata.external_element_tags)
                block.msh==family || continue
                for uid in labels
                    uid in seen && throw(ArgumentError(
                        "$caller: hybrid constructor repeats a Surface cell UID"))
                    push!(seen,uid)
                    entry=get(current,uid,nothing)
                    entry!==nothing && entry[1]==family || throw(ArgumentError(
                        "$caller: hybrid constructor names a foreign actual Surface cell"))
                    ids=entry[2]
                    if family==2
                        push!(triangles,(owner=owner,uid=uid,nodes=(ids[1],ids[2],ids[3]),
                            rows=ntuple(slot->positions[ids[slot]],3)))
                    else
                        push!(quads,(owner=owner,uid=uid,nodes=ids))
                    end
                    push!(cell_uids,uid)
                end
            end
        end
        length(seen)==length(current) || throw(ArgumentError(
            "$caller: hybrid constructor omits actual Surface cells"))
    end
    by_owner=Dict{Int32,Vector{Int}}()
    for (ordinal,quad) in enumerate(quads)
        push!(get!(()->Int[],by_owner,quad.owner[2]),ordinal)
    end
    pyramid_order=Int[ordinal for signed_surface in boundaries for ordinal in
        get(by_owner,Int32(abs(signed_surface)),Int[])]
    length(pyramid_order)==length(quads) && allunique(pyramid_order) || throw(ArgumentError(
        "$caller: restored Region face vector loses an actual Quad birth"))
    return (input=authority,boundaries=copy(boundaries),node_uids=UInt64.(node_uids),
        coordinates=copy(coordinates),triangles=triangles,quads=quads,pyramid_order=pyramid_order,
        embedded_point=embedded_point,embedded_curve=embedded_curve)
end
