# Projection uses the completed sweep's actual cell complex. Curved laterals
# are not fitted to a plane or replayed through a different curve sampler.
function _extrude_nonew_cell_counts(mesh,remap=nothing;dimension::Int=3)
    counts=Dict{Tuple{Int,Tuple},Int}()
    if mesh isa Mesh
        cells=dimension==3 ? mesh.tets : mesh.tris
        sizehint!(counts,size(cells,2))
        msh=dimension==3 ? 4 : 2
        for column in axes(cells,2)
            cell=ntuple(k->remap===nothing ? cells[k,column] :
                            remap[cells[k,column]],size(cells,1))
            key=(msh,cell)
            counts[key]=get(counts,key,0)+1
        end
    else
        capacity=0
        for block in mesh.blocks
            block isa ElementBlock && msh_dimension(block.msh)==dimension || continue
            capacity=Base.checked_add(capacity,size(block.nodes,2))
        end
        sizehint!(counts,capacity)
        for block in mesh.blocks
            block isa ElementBlock && msh_dimension(block.msh)==dimension || continue
            for column in axes(block.nodes,2)
                cell=ntuple(k->remap===nothing ? block.nodes[k,column] :
                                remap[block.nodes[k,column]],size(block.nodes,1))
                key=(block.msh,cell)
                counts[key]=get(counts,key,0)+1
            end
        end
    end
    return counts
end

function _extrude_nonew_projection_lookup(mesh,completed,
                                           caller::AbstractString;dimension::Int=3)
    nnodes(mesh)==nnodes(completed.volume) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts input node set does not match the completed sweep"))
    lookup=Dict{NTuple{3,Float64},Int32}()
    sizehint!(lookup,nnodes(mesh))
    for node in 1:nnodes(mesh)
        coordinate=_extrude_nonew_coordinate(mesh,node)
        haskey(lookup,coordinate) && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts input has distinct nodes at one coordinate"))
        lookup[coordinate]=Int32(node)
    end
    remap=Vector{Int32}(undef,nnodes(completed.volume))
    for node in eachindex(remap)
        target=get(lookup,_extrude_nonew_coordinate(completed.volume,node),Int32(0))
        target>0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts input coordinates do not match the completed sweep"))
        remap[node]=target
    end
    _extrude_nonew_cell_counts(mesh;dimension=dimension)==
        _extrude_nonew_cell_counts(completed.volume,remap;dimension=dimension) || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts input cell topology or winding " *
            "does not match the completed sweep"))
    return lookup
end

function _extrude_nonew_surface_projection(m::GeoModel,mesh,surface::Int,
                                             caller::AbstractString;scope=nothing)
    _extrude_is_nonew(_extrude_gate(m,2,surface)) || return nothing
    owner=0
    for region in _extrude_quadtri_regions(m,surface,caller)
        _extrude_is_nonew(_extrude_gate(m,3,region)) || continue
        owner==0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared Surface[$surface] " *
            "requires the shared-region planner"))
        owner=region
    end
    owner==0 && return nothing
    working=deepcopy(m)
    completed=if scope===nothing
        _extrude_nonew_complete_plan(working,owner,caller;_working_model=true)
    else
        _extrude_nonew_require_scope(scope,m,caller)
        scope.volumes[owner]
    end
    lookup=_extrude_nonew_projection_lookup(mesh,
        (volume=completed.surfaces[surface],),caller;dimension=2)
    scope=_ExtrudeNoNewScope(working,
        Dict(owner=>completed),completed.surfaces)
    projected=_model_volume_to_mixed(working,completed.volume,owner;
                                      _extrude_scope=scope)
    source_data=projected.entity_data
    closure=Set{Tuple{Int,Int}}(((2,surface),))
    for curve in _model_projection_surface_curves(working,surface)
        push!(closure,(1,curve))
        for point in working.curves[curve]
            push!(closure,(0,point))
        end
    end
    full_to_part=zeros(Int32,nnodes(projected))
    part_to_full=zeros(Int32,nnodes(mesh))
    for node in 1:nnodes(projected)
        target=get(lookup,_extrude_nonew_coordinate(projected,node),Int32(0))
        full_to_part[node]=target
        target==0 || (part_to_full[target]=Int32(node))
    end
    all(>(0),part_to_full) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts surface nodes are absent from its volume"))
    blocks=ElementBlock[]
    owners=Vector{Int32}[]
    for (index,block) in pairs(projected.blocks)
        block isa ElementBlock || continue
        dimension=msh_dimension(block.msh)
        selected=Int[]
        for column in axes(block.nodes,2)
            (dimension,Int(source_data.block_entities[index][column])) in closure &&
                push!(selected,column)
        end
        isempty(selected) && continue
        nodes=Matrix{Int32}(undef,size(block.nodes,1),length(selected))
        for (column,old) in pairs(selected),row in axes(nodes,1)
            node=full_to_part[block.nodes[row,old]]
            node>0 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts surface closure has an absent curve node"))
            nodes[row,column]=node
        end
        push!(blocks,ElementBlock(block.msh,nodes,block.tags[selected]))
        push!(owners,source_data.block_entities[index][selected])
    end
    entities=Dict{Tuple{Int,Int},MixedEntity}(
        key=>source_data.entities[key] for key in closure)
    node_entities=source_data.node_entities[part_to_full]
    all(entity->entity in closure,node_entities) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts surface has a node outside its entity closure"))
    data=MixedEntityData(entities;node_entities=node_entities,
        node_parametric=source_data.node_parametric[part_to_full],
        external_node_tags=UInt64.(1:nnodes(mesh)),block_entities=owners,
        external_element_tags=_model_projection_external_elements(blocks))
    groups=Set{Tuple{Int,Int}}()
    for ((dimension,_),entity) in entities,tag in entity.physical_tags
        push!(groups,(dimension,Int(tag)))
    end
    names=Dict{Tuple{Int,Int},String}(key=>value for (key,value) in
        projected.physical_names if key in groups)
    output=MixedMesh(mesh.coords,blocks;physical_names=names,entity_data=data,
                     elementary_entities=owners)
    diagnostic=validate(output)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid projected QuadTriNoNewVerts surface — " *
        join(diagnostic.messages,"; ")))
    m.curve_params=working.curve_params
    return output
end

function _extrude_nonew_projection_faces(part,lookup,
                                          caller::AbstractString)
    remap=Vector{Int32}(undef,nnodes(part))
    for node in eachindex(remap)
        remap[node]=lookup[_extrude_nonew_coordinate(part,node)]
    end
    faces=(NTuple{N,Int32} where N)[]
    if part isa Mesh
        for column in axes(part.tris,2)
            push!(faces,ntuple(k->remap[part.tris[k,column]],3))
        end
    else
        for block in part.blocks
            block isa ElementBlock && msh_dimension(block.msh)==2 || continue
            for column in axes(block.nodes,2)
                push!(faces,ntuple(k->remap[block.nodes[k,column]],
                                   size(block.nodes,1)))
            end
        end
    end
    isempty(faces) && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts projected surface has no cells"))
    return faces
end

function _extrude_nonew_projection_context(m::GeoModel,mesh,volume::Int,
                                            caller::AbstractString;
                                            scope=nothing)
    working=deepcopy(m)
    completed=if scope===nothing
        _extrude_nonew_complete_plan(working,volume,caller;_working_model=true)
    else
        _extrude_nonew_require_scope(scope,m,caller)
        scope.volumes[volume]
    end
    lookup=_extrude_nonew_projection_lookup(mesh,completed,caller)
    two_tri=completed.catalog isa _ExtrudeNoNewTwoTriCatalog
    quad_patch=completed.catalog isa _ExtrudeNoNewQuadPatchCatalog
    quad_strip=completed.catalog isa _ExtrudeNoNewQuadStripCatalog
    three_quad_strip=completed.catalog isa _ExtrudeNoNewThreeQuadStripCatalog
    intervals=length(completed.catalog.layer_refs)
    face_capacity=quad_patch ? _extrude_nonew_quad_patch_face_capacity(completed.catalog) :
        quad_strip ? _extrude_nonew_quad_strip_face_capacity(completed.catalog) :
        three_quad_strip ? _extrude_nonew_three_quad_strip_face_capacity(completed.catalog) :
        two_tri ? (completed.volume isa Mesh ? 16 : 7)*intervals+2 : 0
    typed_boundary=_extrude_nonew_certify_boundary(mesh,completed.surfaces,caller;
        face_capacity=face_capacity)
    # Retain the audited input-node boundary instead of rebuilding its full
    # volume-face incidence in the classified projection.
    boundary=Set{NTuple{N,Int32} where N}()
    sizehint!(boundary,length(typed_boundary))
    for face in typed_boundary
        push!(boundary,face[4]==0 ? (face[1],face[2],face[3]) : face)
    end
    curves=Set{Int}()
    points=Set{Int}()
    for surface in keys(completed.surfaces)
        for loop in working.surfaces[surface],signed in working.loops[loop]
            curve=abs(signed)
            push!(curves,curve)
            a,b=working.curves[curve]
            push!(points,a);push!(points,b)
        end
    end
    point_nodes=Dict{Int,Int32}()
    for point in points
        coordinate=ntuple(k->_model_projection_coordinate_key(working.points[point][k]),3)
        node=get(lookup,coordinate,Int32(0))
        node>0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts Point[$point] is absent from the completed sweep"))
        point_nodes[point]=node
    end
    # Parameter grading records the native associations; connectivity comes
    # from the retained logical columns, rather than re-evaluated helix chords.
    curve_tags=sort!(collect(curves))
    _model_surface_mesh_curves!(working,curve_tags,caller)
    source=abs(working.meshing.extrude_sources[(3,volume)][2])
    source_curves=Set(_model_projection_surface_curves(working,source))
    source_columns=Dict{NTuple{3,Float64},Int}()
    columns=completed.columns
    for column in axes(columns,2)
        coordinate=ntuple(k->_model_projection_coordinate_key(columns[1,column][k]),3)
        source_columns[coordinate]=column
    end
    edges=_volume_edge_set(mesh)
    entries=Dict{Int,Vector{Tuple{Float64,Int}}}()
    for curve in curve_tags
        a,b=working.curves[curve]
        link=get(working.meshing.extrude_sources,(1,curve),nothing)
        chain=if (quad_patch || quad_strip || three_quad_strip) &&
                (curve in source_curves || (link!==nothing && link[1]==1))
            original=curve in source_curves ? curve : abs(link[2])
            position=findfirst(==(original),completed.catalog.boundary_curves)
            position===nothing && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts Curve[$curve] has no actual source boundary chain"))
            level=curve in source_curves ? 1 : size(columns,1)
            width=(quad_strip || three_quad_strip) ?
                Int(completed.catalog.boundary_widths[position]) : 3
            nodes=Int32[lookup[ntuple(k->_model_projection_coordinate_key(
                columns[level,completed.catalog.boundary_chains[position][index]][k]),3)]
                for index in 1:width]
            nodes[1]==point_nodes[a] || reverse!(nodes)
            nodes
        elseif curve in source_curves || (link!==nothing && link[1]==1)
            Int32[point_nodes[a],point_nodes[b]]
        elseif link!==nothing && link[1]==0
            coordinate=ntuple(k->_model_projection_coordinate_key(
                working.points[abs(link[2])][k]),3)
            column=get(source_columns,coordinate,0)
            column>0 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts Curve[$curve] has no source column"))
            nodes=Int32[lookup[ntuple(k->_model_projection_coordinate_key(p[k]),3)]
                        for p in @view columns[:,column]]
            nodes[1]==point_nodes[a] || reverse!(nodes)
            nodes
        else
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts Curve[$curve] lacks boundary provenance"))
        end
        chain[1]==point_nodes[a] && chain[end]==point_nodes[b] || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts Curve[$curve] endpoints do not match its column"))
        parameters=working.curve_params[curve]
        length(parameters)==length(chain) || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts Curve[$curve] parameters do not match its column"))
        for index in 1:length(chain)-1
            u,v=chain[index],chain[index+1]
            (u<v ? (u,v) : (v,u)) in edges || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts Curve[$curve] column is not a mesh-edge chain"))
        end
        entries[curve]=Tuple{Float64,Int}[(Float64(parameters[index]),Int(chain[index]))
                                        for index in eachindex(chain)]
    end
    surfaces=Dict{Int,Vector{NTuple{N,Int32} where N}}()
    for (surface,part) in completed.surfaces
        surfaces[surface]=_extrude_nonew_projection_faces(part,lookup,caller)
    end
    return (curves=entries,surfaces=surfaces,curve_params=working.curve_params,
            edges=edges,boundary=boundary)
end
