# An operation-local proof relation preserves each actual Surface polygon.
# A retained Quad may own two oriented Region facets on its chosen diagonal;
# this relation creates no Surface triangles, vertices, UIDs or native clocks.

function _extrude_nonew_retained_b4_certificate(m::GeoModel,t::Int,volume,surfaces,
        catalog,caller;face_capacity::Int=0)
    _extrude_nonew_retained_b4_catalog(catalog) &&
        _extrude_nonew_unmerged_b4_incidence(m,t,caller) || throw(ArgumentError(
            "$caller: split Surface coverage lacks its retained B4 operation"))
    _extrude_nonew_scope_contract(m,t,caller)
    params,spec,link=_extrude_entity_params(m,3,t,caller)
    params.recomb_laterals && link[1]==2 && spec.type===:translate ||
        throw(ArgumentError("$caller: split Surface coverage changes its recorded B4 policy"))
    source=abs(link[2])
    top,laterals=_extrude_nonew_region_boundary(m,t,source,caller;sides=4)
    boundary=_extrude_nonew_boundary_surfaces(m,t,caller)
    Set(keys(surfaces))==Set(abs.(boundary)) || throw(ArgumentError(
        "$caller: split Surface coverage loses its Source/cap/lateral constructors"))
    levels,refs=_extrude_nonew_levels(params,caller;
        source_nodes=length(catalog.source.source_coordinates),extra_nodes=0,
        cells_per_interval=6*length(catalog.source.source_cells))
    levels==catalog.levels && refs==catalog.layer_refs || throw(ArgumentError(
        "$caller: split Surface coverage changes its recorded layer constructor"))
    current_source=_extrude_nonew_rect_grid_source(m,source,surfaces[source],spec,caller;
        mode=:b4_strip,bound=catalog.source)
    all(name->isequal(getfield(current_source,name),getfield(catalog.source,name)),
        fieldnames(_ExtrudeNoNewRectGridSource)) || throw(ArgumentError(
            "$caller: split Surface coverage changes its bound Source constructor"))
    nnodes(volume)==length(catalog.source.source_coordinates)*length(catalog.levels) ||
        throw(ArgumentError("$caller: split Surface coverage introduces an unowned body row"))
    columns=_extrude_volume_columns(m,t,source,surfaces[source],params,spec,levels,caller)
    actual_coordinates=Set(_extrude_nonew_coordinate(volume,row) for row in 1:nnodes(volume))
    Set(_extrude_nonew_rect_grid_key(point) for point in columns)==actual_coordinates ||
        throw(ArgumentError("$caller: split Surface coverage changes its current swept coordinates"))
    orientation=Int(catalog.source.source_orientation)*Int(catalog.source.normal_direction)
    orientation in (-1,1) || throw(ArgumentError("$caller: split Surface coverage lacks its source orientation"))
    signs=Dict(abs(surface)=>sign(surface)*orientation for surface in boundary)
    for surface in laterals
        signed_curve=m.meshing.extrude_sources[(2,surface)][2]
        signs[surface]*=sign(signed_curve)
    end
    return _extrude_nonew_split_macroface_certificate(volume,surfaces,signs,caller;face_capacity)
end

function _extrude_nonew_incidence_cycle(key,incidence::_ExtrudeNoNewFaceIncidence)
    incidence.count==1 || throw(ArgumentError("retained boundary: facet is not exterior"))
    orientation=incidence.orientation
    if key[4]==0
        return orientation==1 ? (key[1],key[2],key[3]) : (key[1],key[3],key[2])
    end
    remaining=(key[2],key[3],key[4])
    opposite=remaining[Int(orientation>>1)+1]
    others=filter(!=(opposite),remaining)
    return isodd(orientation) ? (key[1],others[1],opposite,others[2]) :
        (key[1],others[2],opposite,others[1])
end

function _extrude_nonew_oriented_edges!(edges,cycle,sign::Int=1)
    for position in eachindex(cycle)
        first=cycle[position];second=cycle[mod1(position+1,length(cycle))]
        first!=second || throw(ArgumentError("retained boundary: collapsed exterior edge"))
        key=minmax(first,second)
        edges[key]=get(edges,key,0)+(first<second ? sign : -sign)
    end
    return edges
end

function _extrude_nonew_split_macroface_certificate(volume,surfaces,signs,
        caller::AbstractString;face_capacity::Int=0)
    Set(keys(signs))==Set(keys(surfaces)) && all(sign->sign in (-1,1),values(signs)) ||
        throw(ArgumentError("$caller: retained Surface orientation receipts do not align"))
    counts=Dict{NTuple{4,Int32},_ExtrudeNoNewFaceIncidence}()
    face_capacity>0 && sizehint!(counts,face_capacity)
    if volume isa Mesh
        isempty(volume.segs) && isempty(volume.tris) || throw(ArgumentError(
            "$caller: retained Region constructor contains foreign cells"))
        _extrude_nonew_count_faces!(counts,volume.tets,4,caller)
    else
        for block in volume.blocks
            block isa ElementBlock || throw(ArgumentError(
                "$caller: retained Region constructor contains unsupported blocks"))
            isempty(block.nodes) && continue
            msh_dimension(block.msh)==3 || throw(ArgumentError(
                "$caller: retained Region constructor contains foreign cells"))
            expected=get(Dict(4=>4,5=>8,6=>6,7=>5),Int(block.msh),0)
            expected>0 && size(block.nodes,1)==expected || throw(ArgumentError(
                "$caller: retained Region cell rows do not match their native family"))
            _extrude_nonew_count_faces!(counts,block.nodes,Int(block.msh),caller)
        end
    end
    isempty(counts) && throw(ArgumentError("$caller: retained Region constructor has no cells"))
    boundary=Dict(key=>value for (key,value) in counts if value.count==1)
    shell_edges=Dict{NTuple{2,Int32},Int}()
    for (key,incidence) in boundary
        _extrude_nonew_oriented_edges!(shell_edges,_extrude_nonew_incidence_cycle(key,incidence))
    end
    all(iszero,values(shell_edges)) || throw(ArgumentError(
        "$caller: retained Region exterior is not an oriented closed shell"))
    lookup=Dict{NTuple{3,Float64},Int32}()
    for row in 1:nnodes(volume)
        coordinate=_extrude_nonew_coordinate(volume,row)
        haskey(lookup,coordinate) && throw(ArgumentError(
            "$caller: retained Region has distinct rows at one coordinate"))
        lookup[coordinate]=Int32(row)
    end
    actual=Set{NTuple{4,Int32}}();used=Set{NTuple{4,Int32}}()
    relations=Dict{NTuple{4,Int32},NTuple{2,NTuple{4,Int32}}}()
    owners=Dict{NTuple{4,Int32},Int}()
    function admit(nodes,remap,tag)
        width=size(nodes,1)
        width in (3,4) || throw(ArgumentError("$caller: retained Surface has unsupported cell width"))
        for column in axes(nodes,2)
            polygon=ntuple(slot->remap[nodes[slot,column]],width)
            allunique(polygon) || throw(ArgumentError("$caller: retained Surface polygon collapses"))
            key=_extrude_nonew_face_key(polygon)
            key in actual && throw(ArgumentError("$caller: retained Surface polygon is repeated"))
            push!(actual,key)
            parts=NTuple{4,Int32}[]
            if width==4
                for omit in 1:4
                    triangle=Tuple(polygon[slot] for slot in 1:4 if slot!=omit)
                    candidate=_extrude_nonew_face_key(triangle)
                    haskey(boundary,candidate) && push!(parts,candidate)
                end
            end
            if haskey(boundary,key)
                isempty(parts) || throw(ArgumentError(
                    "$caller: retained Surface Quad overlaps split Region facets"))
                push!(parts,key)
            else
                width==4 && length(parts)==2 || throw(ArgumentError(
                    "$caller: retained Surface polygon lacks exact Region facet coverage"))
                left=Set(filter(!iszero,parts[1]));right=Set(filter(!iszero,parts[2]))
                union(left,right)==Set(polygon) && length(intersect(left,right))==2 ||
                    throw(ArgumentError("$caller: retained Quad split loses a corner or its internal diagonal"))
                exterior=Set(minmax(polygon[i],polygon[mod1(i+1,4)]) for i in 1:4)
                diagonal=Tuple(sort!(collect(intersect(left,right))))
                diagonal in exterior && throw(ArgumentError(
                    "$caller: retained Quad split places its shared edge on the exterior"))
                relations[key]=(parts[1],parts[2])
            end
            expected_edges=Dict{NTuple{2,Int32},Int}()
            _extrude_nonew_oriented_edges!(expected_edges,polygon,Int(signs[tag]))
            facet_edges=Dict{NTuple{2,Int32},Int}()
            for part in parts
                part in used && throw(ArgumentError("$caller: retained Region facet is covered twice"))
                push!(used,part);owners[part]=tag
                _extrude_nonew_oriented_edges!(facet_edges,
                    _extrude_nonew_incidence_cycle(part,boundary[part]))
            end
            filter!(entry->!iszero(entry.second),facet_edges)
            expected_edges==facet_edges || throw(ArgumentError(
                "$caller: retained Surface coverage has inconsistent outward orientation"))
        end
    end
    for tag in sort!(collect(keys(surfaces)))
        part=surfaces[tag]
        remap=Int32[get(lookup,_extrude_nonew_coordinate(part,row),Int32(0)) for row in 1:nnodes(part)]
        all(>(0),remap) || throw(ArgumentError(
            "$caller: retained Surface row is absent from its Region constructor"))
        if part isa Mesh
            isempty(part.segs) && isempty(part.tets) || throw(ArgumentError(
                "$caller: retained Surface constructor contains foreign cells"))
            admit(part.tris,remap,tag)
        else
            for block in part.blocks
                block isa ElementBlock || throw(ArgumentError(
                    "$caller: retained Surface constructor contains unsupported blocks"))
                isempty(block.nodes) && continue
                block.msh in (2,3) || throw(ArgumentError(
                    "$caller: retained Surface constructor contains unsupported cells"))
                size(block.nodes,1)==(block.msh==2 ? 3 : 4) || throw(ArgumentError(
                    "$caller: retained Surface cell rows do not match their native family"))
                admit(block.nodes,remap,tag)
            end
        end
    end
    used==Set(keys(boundary)) || throw(ArgumentError(
        "$caller: retained Surface polygons do not completely cover the Region exterior"))
    return (boundary=Set(keys(boundary)),surface_polygons=actual,relations=relations,
        facet_surfaces=owners)
end
