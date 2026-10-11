# The native Curve emitter already provides endpoint symbols and an ordered
# parameter vector. Bind those source identities before any dense carrier row
# can be interpreted as a public node label.
function _geo_dim01_producer_trace(m,merged,parts,remaps,caller)
    length(parts)==length(remaps) || throw(ArgumentError(
        "$caller: native Curve products do not align with their emitted rows"))
    owners=Dict{Int,Tuple{Int,Int}}()
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(merged)]
    for ((dimension,tag,part),rows) in zip(parts,remaps)
        dimension in (0,1) && part isa Mesh && length(rows)==nnodes(part) ||
            throw(ArgumentError("$caller: native dim01 producer has invalid entity rows"))
        values=if dimension==0
            haskey(m.points,tag) && length(rows)==1 || throw(ArgumentError(
                "$caller: native Point[$tag] lacks its explicit emitted vertex"))
            Float64[]
        else
            haskey(m.curves,tag) && haskey(m.curve_params,tag) || throw(ArgumentError(
                "$caller: native Curve[$tag] lacks its emitted endpoint/parameter trace"))
            emitted=Float64.(m.curve_params[tag])
            a,b=m.curves[tag]
            if a==b && length(emitted)>1
                _,upper=_model_curve_param_bounds(m,tag,caller)
                emitted[end]>=upper && pop!(emitted)
            end
            length(emitted)==length(rows) || throw(ArgumentError(
                "$caller: native Curve[$tag] parameter births do not align with its rows"))
            emitted
        end
        for (slot,row) in enumerate(rows)
            0<row<=nnodes(merged) || throw(ArgumentError(
                "$caller: native producer references an invalid carrier row"))
            owner=(dimension,tag)
            parameter=dimension==0 ? Float64[] : Float64[values[slot]]
            if dimension==1
                a,b=m.curves[tag]
                if slot==1
                    owner=(0,a);parameter=Float64[]
                elseif a!=b && slot==length(rows)
                    owner=(0,b);parameter=Float64[]
                end
            end
            existing=get(owners,Int(row),owner)
            existing==owner || throw(ArgumentError(
                "$caller: carrier row $row merges distinct native producer identities"))
            # Coordinate equality validates the saved emission row. It never
            # chooses its owner or public identity.
            for axis in 1:3
                isequal(merged.coords[axis,row],part.coords[axis,slot]) || throw(ArgumentError(
                    "$caller: native producer row changed during geometry merge"))
            end
            owners[Int(row)]=owner
            parameters[row]=parameter
        end
    end
    length(owners)==nnodes(merged) || throw(ArgumentError(
        "$caller: native dim01 producer omits a geometry row"))
    order=sort!(unique!([(dimension,tag) for (dimension,tag,_) in parts]))
    return (;owners,parameters,order)
end

# Native Arc vertices and the legacy geometry carrier share the same explicit
# endpoint/Curve-emitter births. The native positive transfinite law samples
# the derivative primitive; the already emitted carrier keeps its old law.
# This native TF primitive evaluates geometry only. It never repeats the
# already completed field or user-callback grading that emitted the carrier.
function _geo_dim01_native_curve_births!(m,merged,parts,remaps,produced,options,caller)
    coordinates=copy(merged.coords)
    parameters=Union{Nothing,Vector{Float64}}[
        value===nothing ? nothing : copy(value) for value in produced.parameters]
    for ((dimension,tag,_),rows) in zip(parts,remaps)
        dimension==1 || continue
        spec=get(m.meshing.transfinite_curves,tag,nothing)
        spec!==nothing && spec.coef>0.0 || continue
        adapter=Model._model_native_arc_mesh_geometry(m,tag,caller)
        adapter===nothing && continue
        values=Model._model_native_arc_transfinite_params(m,tag,spec,options,caller)
        emitted=copy(values)
        a,b=m.curves[tag]
        a==b && length(emitted)>1 && emitted[end]>=adapter.upper && pop!(emitted)
        length(emitted)==length(rows) || throw(ArgumentError(
            "$caller: native and geometry Arc births have different row counts"))
        for (slot,row) in enumerate(rows)
            produced.owners[Int(row)]==(1,tag) || continue
            point=Model._model_high_order_curve_point(adapter,emitted[slot],caller)
            all(isfinite,point) && isfinite(emitted[slot]) || throw(ArgumentError(
                "$caller: native Arc birth has nonfinite coordinates or parameter"))
            for axis in 1:3
                coordinates[axis,row]=point[axis]
            end
            parameters[row]=Float64[emitted[slot]]
        end
        m.curve_params[tag]=values
    end
    return (;coordinates,parameters,geometry_parameters=produced.parameters)
end

# Re-generation of a bound linear dim01 carrier consumes the completed actual
# products. Emitter ranks select row order; public UIDs select row bindings.
# No geometric evaluation or coordinate search occurs in this projection.
function _geo_dim01_generation_carrier!(plan,previous,caller)
    previous!==nothing && previous.carrier isa Mesh &&
        ntris(previous.carrier)==0 && ntets(previous.carrier)==0 || return plan
    authority=plan.model.meshing.mesh_authority
    mesh=authority.mesh;data=mesh.entity_data
    geometry=hasproperty(plan,:geometry_carrier) ? plan.geometry_carrier : mesh
    geometry isa MixedMesh && geometry.entity_data!==nothing ||
        throw(ArgumentError("$caller: completed dim01 geometry lacks explicit vertex identities"))
    geometry_positions=Dict(tag=>row for (row,tag) in
        enumerate(geometry.entity_data.external_node_tags))
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(data.external_node_tags))
    old_slots=Dict(row=>slot for (slot,row) in enumerate(previous.carrier_rows))
    point_order=Tuple{Int,Int32}[]
    curve_order=Tuple{Int,Int32}[]
    for owner in keys(authority.products)
        owner[1]==0 && push!(point_order,owner)
        owner[1]==1 && push!(curve_order,owner)
    end
    sort!(point_order;by=owner->begin
        product=get(previous.products,owner,nothing)
        first_slot=product===nothing ? typemax(Int) :
            minimum((get(old_slots,row,typemax(Int)) for row in product.rows);init=typemax(Int))
        (first_slot,owner[2])
    end)
    sort!(curve_order)
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[]
    node_entities=Vector{Tuple{Int,Int}}[]
    node_identities=Vector{UInt64}[]
    cell_identities=Vector{Vector{UInt64}}[]
    for owner in Iterators.flatten((point_order,curve_order))
        product=authority.products[owner]
        rows=owner[1]==0 ? Int32[row for row in product.rows if data.node_entities[row]==owner] :
            copy(product.rows)
        isempty(rows) && continue
        selected=Int[index for (index,block) in enumerate(product.blocks) if block.msh==1]
        owner[1]==1 && isempty(selected) && continue
        local_rows=Dict(row=>Int32(slot) for (slot,row) in enumerate(rows))
        coordinates=Matrix{Float64}(undef,3,length(rows))
        tags=data.external_node_tags[rows]
        for (slot,tag) in enumerate(tags)
            source=get(geometry_positions,tag,0)
            source>0 || throw(ArgumentError("$caller: geometry omits actual dim01 vertex $tag"))
            coordinates[:,slot]=geometry.coords[:,source]
        end
        count=sum(size(product.blocks[index].nodes,2) for index in selected;init=0)
        segments=Matrix{Int32}(undef,2,count);physical=Int32[];cells=UInt64[];column=0
        for index in selected
            block=product.blocks[index]
            for source in axes(block.nodes,2)
                column+=1
                for slot in 1:2
                    segments[slot,column]=local_rows[product.rows[block.nodes[slot,source]]]
                end
            end
            append!(physical,block.tags);append!(cells,product.element_tags[index])
        end
        push!(parts,(owner[1],Int(owner[2]),Mesh(coordinates;segs=segments,seg_tag=physical)))
        push!(node_entities,Tuple{Int,Int}[(data.node_entities[row][1],Int(data.node_entities[row][2])) for row in rows])
        push!(node_identities,copy(tags));push!(cell_identities,[cells,UInt64[],UInt64[]])
    end
    isempty(parts) && return plan
    carrier,_,remaps,cells=_geo_merge_entity_meshes(parts;node_entities,node_identities,
        cell_identities,return_remaps=true)
    rows=zeros(Int32,nnodes(carrier))
    for (tags,remap) in zip(node_identities,remaps),(slot,row) in enumerate(remap)
        actual=positions[tags[slot]]
        rows[row]==0 || rows[row]==actual || throw(ArgumentError(
            "$caller: geometry merges distinct dim01 producer identities"))
        rows[row]=actual
    end
    all(>(0),rows) && allunique(rows) || throw(ArgumentError(
        "$caller: generated dim01 carrier loses explicit vertex bindings"))
    rebound=Model._model_mesh_authority(mesh,authority.done;caller,
        private_rank=authority.private_rank,carrier,carrier_rows=rows,carrier_cells=cells)
    plan.model.meshing.mesh_authority=Model._ModelMeshAuthority(rebound.mesh,
        rebound.products,rebound.done,rebound.private_rank,authority.native_inputs,
        carrier,rows,authority.producer_rows,rebound.carrier_cells)
    return merge(plan,(geometry_parts=parts,geometry_part_rows=remaps))
end
