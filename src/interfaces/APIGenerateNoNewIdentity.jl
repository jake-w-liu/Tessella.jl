# Unmerged NoNew B4 products retain their literal CAD aliases. Bind only this
# certified operation's emitted position rows to the current lower cell UIDs.
# The selected UID owns the row; other coincident native records stay retained.
function _higher_nonew_identity_classification(assembly,m,mesh,input,scope,dimension,tag,caller)
    Model._extrude_is_nonew(Model._extrude_gate(m,dimension,tag)) || return nothing
    scope isa Model._ExtrudeNoNewScope || return nothing
    Model._extrude_nonew_require_scope(scope,m,caller)
    owners=dimension==3 ? Int[tag] : Model._extrude_quadtri_regions(m,tag,caller)
    selected=Int[volume for volume in owners if
        Model._extrude_nonew_unmerged_b4_incidence(m,volume,caller)]
    isempty(selected) && return nothing
    length(selected)==1 || throw(ArgumentError(
        "$caller: unmerged NoNew product requires one isolated operation"))
    volume=only(selected)
    completed=get(scope.volumes,volume,nothing)
    completed isa Model._ExtrudeNoNewCompletePlan || throw(ArgumentError(
        "$caller: unmerged NoNew product lacks its typed complete scope"))
    completed.catalog isa Union{Model._ExtrudeNoNewB4StripCatalog,Model._ExtrudeNoNewB4Free.Catalog} ||
        throw(ArgumentError("$caller: unmerged NoNew identity requires its indexed B4 constructor"))
    Model._extrude_nonew_scope_contract(m,volume,caller)
    boundary=Model._extrude_nonew_boundary_surfaces(m,volume,caller)
    source=abs(m.meshing.extrude_sources[(3,volume)][2])
    top,laterals=Model._extrude_nonew_region_boundary(m,volume,source,caller)
    Set(keys(completed.surfaces))==Set(abs.(boundary)) || throw(ArgumentError(
        "$caller: unmerged NoNew scope loses its selected Source/cap/laterals"))
    expected=dimension==3 ? completed.volume : get(completed.surfaces,tag,nothing)
    expected!==nothing || throw(ArgumentError(
        "$caller: unmerged NoNew product is absent from its selected scope"))
    emitted=Model._extrude_nonew_projection_lookup(mesh,(volume=expected,),caller;dimension)
    retained_b4=Model._extrude_nonew_retained_b4_catalog(completed.catalog)
    proof=retained_b4 ? Model._extrude_nonew_retained_b4_certificate(
        m,volume,dimension==3 ? mesh : completed.volume,
        completed.surfaces,completed.catalog,caller;
        face_capacity=completed.catalog.face_capacity) : nothing
    retained_b4 || Model._extrude_nonew_certify_boundary(completed.volume,completed.surfaces,caller;
        oriented_internal=true,face_capacity=completed.catalog.face_capacity)
    input isa Model._ModelMeshAuthority && input.mesh!==nothing ||
        throw(ArgumentError("$caller: unmerged NoNew operation lacks current lower authority"))
    actual=input.mesh;data=actual.entity_data
    source_product=Model._model_mesh_product(input,2,source;caller,producer_order=true)
    bound_source=completed.surfaces[source]
    source_product!==nothing && bound_source isa MixedMesh && bound_source.entity_data!==nothing ||
        throw(ArgumentError("$caller: unmerged NoNew Source lacks its current UID constructor"))
    Model._model_surface_constructor_binding(source_product,bound_source,caller)
    primary=(dimension,Int32(tag))
    closure=Set{Tuple{Int,Int32}}((primary,))
    lower=dimension==2 ? [(1,abs(curve)) for loop in m.surfaces[tag] for curve in m.loops[loop]] :
        [(2,source),(2,top),((2,surface) for surface in laterals)...]
    for (dim,entity) in lower
        push!(closure,(dim,Int32(entity)))
        curves=dim==2 ? Model._model_projection_surface_curves(m,entity) : [entity]
        for curve in curves
            push!(closure,(1,Int32(curve)))
            for point in m.curves[curve];push!(closure,(0,Int32(point)));end
        end
    end
    node_tags=zeros(UInt64,nnodes(mesh));retained=falses(nnodes(mesh))
    node_owners=fill(primary,nnodes(mesh))
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    function bind_lower(row)
        uid=data.external_node_tags[row]
        point=Model._extrude_nonew_coordinate(actual,row)
        target=get(emitted,point,Int32(0))
        target>0 || throw(ArgumentError(
            "$caller: unmerged NoNew emitted rows omit lower cell UID $uid"))
        position=get(assembly.node_positions,uid,Int32(0))
        position>0 || throw(ArgumentError(
            "$caller: unmerged NoNew lower UID $uid lacks its native record"))
        current=assembly.nodes[position]
        current.owner==data.node_entities[row] &&
            isequal(Model._extrude_nonew_rect_grid_key(current.point),point) ||
            throw(ArgumentError("$caller: unmerged NoNew lower UID changed ownership or coordinates"))
        current.owner in closure || throw(ArgumentError(
            "$caller: unmerged NoNew lower UID belongs outside the selected operation"))
        if !retained[target]
            node_tags[target]=uid;retained[target]=true;node_owners[target]=current.owner
            parameters[target]=current.parameter===nothing ? nothing : copy(current.parameter)
        end
        return target
    end
    edges=Dict{NTuple{2,Int32},Tuple{Int,Int32}}()
    faces=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    quads=Dict{NTuple{4,Int32},Tuple{Int,Int32}}()
    edge_conflicts=Set{NTuple{2,Int32}}();face_conflicts=Set{NTuple{3,Int32}}()
    quad_conflicts=Set{NTuple{4,Int32}}();cell_owners=Dict{Int,Vector{Int32}}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        isempty(cells) && continue
        msh_dimension(msh)==dimension || throw(ArgumentError(
            "$caller: unmerged NoNew constructor mixes cell dimensions"))
        values=fill(Int32(tag),size(cells,2))
        append!(get!(()->Int32[],cell_owners,msh),values)
        _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,msh,cells,values)
    end
    lower_edges=Set{NTuple{2,Int32}}();lower_faces=Set{Tuple}()
    for (dim,entity) in lower
        owner=(dim,Int32(entity));product=get(input.products,owner,nothing)
        product!==nothing && owner in input.done && product.owner==owner ||
            throw(ArgumentError("$caller: unmerged NoNew lower product $owner is not current and complete"))
        !isempty(product.blocks) && any(block->!isempty(block.nodes),product.blocks) ||
            throw(ArgumentError("$caller: unmerged NoNew lower product has no actual cells"))
        for block in product.blocks
            msh_dimension(block.msh)==dim && block.msh in (1,2,3) ||
                throw(ArgumentError("$caller: unmerged NoNew lower product has unsupported cells"))
            mapped=Matrix{Int32}(undef,size(block.nodes))
            for column in axes(block.nodes,2),slot in axes(block.nodes,1)
                mapped[slot,column]=bind_lower(product.rows[block.nodes[slot,column]])
            end
            for column in axes(mapped,2)
                if dim==1
                    edge=minmax(mapped[1,column],mapped[2,column])
                    edge[1]!=edge[2] && !(edge in lower_edges) || throw(ArgumentError(
                        "$caller: unmerged NoNew lower Curve cells overlap or collapse"))
                    push!(lower_edges,edge)
                else
                    face=Tuple(sort!(collect(@view mapped[:,column])))
                    allunique(face) && !(face in lower_faces) || throw(ArgumentError(
                        "$caller: unmerged NoNew lower Surface cells overlap or collapse"))
                    push!(lower_faces,face)
                end
            end
            _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
                Int(block.msh),mapped,fill(Int32(entity),size(mapped,2)))
        end
    end
    if dimension==3
        # The Region's selected Surface rows have already chosen their UIDs.
        # Resolve lower Curve support in the same Source/cap/lateral order;
        # recorded copied/connecting aliases do not claim that support twice.
        ordered_curves=Int[]
        for (_,surface) in lower,loop in m.surfaces[surface],curve in m.loops[loop]
            abs(curve) in ordered_curves || push!(ordered_curves,abs(curve))
        end
        actual_edges=Model._volume_edge_set(mesh)
        for curve in ordered_curves
            product=get(input.products,(1,Int32(curve)),nothing)
            product!==nothing && (1,Int32(curve)) in input.done || throw(ArgumentError(
                "$caller: unmerged NoNew Curve support lacks its current product"))
            for block in product.blocks
                block.msh==1 || throw(ArgumentError(
                    "$caller: unmerged NoNew Curve support requires actual Lines"))
                for column in axes(block.nodes,2)
                    pair=ntuple(slot->bind_lower(product.rows[block.nodes[slot,column]]),2)
                    edge=minmax(pair...)
                    edge[1]!=edge[2] && edge in actual_edges || throw(ArgumentError(
                        "$caller: unmerged NoNew actual Curve support is absent from the Region"))
                    get(edges,edge,(4,Int32(0)))[1]>1 &&
                        _api_support_insert!(edges,edge_conflicts,edge,(1,Int32(curve)))
                end
            end
        end
    end
    _api_check_support_conflicts(edge_conflicts,face_conflicts,quad_conflicts)
    if dimension==2
        boundary_rows,boundary_edges=Model._surface_boundary_topology(mesh,caller)
        lower_edges==boundary_edges || throw(ArgumentError(
            "$caller: unmerged NoNew actual Curve cells do not cover the typed Face boundary"))
        all(row->!boundary_rows[row] || retained[row],eachindex(boundary_rows)) ||
            throw(ArgumentError("$caller: unmerged NoNew Face boundary loses actual lower UIDs"))
    else
        Model._volume_projection_contract(mesh,caller)
        typed_faces,_=Model._volume_cell_boundary_faces(mesh)
        if retained_b4
            actual_polygons=Set(Model._extrude_nonew_face_key(face) for face in lower_faces)
            actual_polygons==proof.surface_polygons || throw(ArgumentError(
                "$caller: retained NoNew actual Surface cells change the certified macroface coverage"))
            all(retained) || throw(ArgumentError(
                "$caller: retained NoNew Region introduces a node outside its actual lower products"))
            for (facet,surface) in proof.facet_surfaces
                facet[4]==0 || continue
                key=(facet[1],facet[2],facet[3])
                _api_support_insert!(faces,face_conflicts,key,(2,Int32(surface)))
            end
            _api_check_support_conflicts(edge_conflicts,face_conflicts,quad_conflicts)
        else
            lower_faces==typed_faces || throw(ArgumentError(
                "$caller: unmerged NoNew actual Surface cells do not cover the typed Region boundary"))
        end
        all(row->retained[row],(row for face in typed_faces for row in face)) ||
            throw(ArgumentError("$caller: unmerged NoNew Region boundary loses actual lower UIDs"))
    end
    boundaries=Dict(owner=>Int32[abs(child) for child in data.entities[owner].boundaries]
        for owner in closure)
    class=_mixed_classification(mesh,primary,sort!(collect(closure)),node_owners,boundaries,cell_owners;
        edge_entities=edges,face_entities=faces,quad_entities=quads,node_parametric=parameters)
    return (class=class,node_tags=node_tags,retained=retained)
end
