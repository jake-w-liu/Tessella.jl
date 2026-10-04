# Indexed recombined all-boundary strips retain actual source columns and add
# one certified terminal center per original Quad in source-cell order.

function _extrude_nonew_b4_strip_preflight(shape,caller)
    try
        vertices,cells,_,_=_extrude_nonew_rect_grid_sizes(shape[1],shape[2],caller;
            mode=:b4_strip)
        return vertices,Base.checked_mul(7,cells),cells
    catch err
        err isa OverflowError || rethrow()
        _extrude_nonew_b4_strip_error(caller,"source dimensions exceed the indexed size limit")
    end
end

@noinline function _extrude_nonew_b4_strip_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts all-boundary strip $reason"))
end

function _extrude_nonew_b4_strip_emit!(arrays,filled,template,original,
        interval::Int,nlevels::Int,center::Int32,caller)
    for position in 1:Int(template.ncells)
        cell=template.cells[position];family=Int(cell.msh)-3
        1<=family<=4 || _extrude_nonew_b4_strip_error(caller,
            "factory emitted an unsupported cell family")
        nodes=arrays[family];filled[family]+=1;column=filled[family]
        column<=size(nodes,2) || _extrude_nonew_b4_strip_error(caller,
            "factory exceeded its reserved cell capacity")
        for row in axes(nodes,1)
            local_node=Int(cell.nodes[row])
            if local_node==9
                center>0 || _extrude_nonew_b4_strip_error(caller,
                    "factory has an unallocated actual centroid")
                nodes[row,column]=center
            else
                1<=local_node<=8 || _extrude_nonew_b4_strip_error(caller,
                    "factory has an unsupported local vertex")
                source_position=local_node<=4 ? local_node : local_node-4
                level=local_node<=4 ? interval : interval+1
                nodes[row,column]=_extrude_nonew_rect_grid_node(
                    original[source_position],level,nlevels)
            end
        end
    end
    return nothing
end

function _extrude_nonew_b4_strip_finish(m::GeoModel,t::Int,params,source_tag::Int,
        source_mesh,catalog::_ExtrudeNoNewB4StripCatalog,cols,top,laterals,caller)
    nlevels,intervals=_extrude_nonew_rect_grid_dimensions(cols,catalog,caller)
    source=catalog.source;ncells=length(source.source_cells)
    params.recomb_laterals && catalog.recombined ||
        _extrude_nonew_b4_strip_error(caller,
            "free laterals require the constructive all-boundary face phase")
    length(catalog.top_states)==ncells || _extrude_nonew_b4_strip_error(caller,
        "top-state count differs from its actual source")
    column_count=Base.checked_mul(length(source.source_coordinates),nlevels)
    node_count=Base.checked_add(column_count,ncells)
    cell_count=foldl(Base.checked_add,catalog.cell_counts;init=0)
    node_count<=_EXTRUDE_NONEW_MAX_NODES && node_count<=typemax(Int32) &&
        cell_count<=typemax(Int32) || _extrude_nonew_b4_strip_error(caller,
            "indexed output exceeds the checked node or Int32 cell limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in eachindex(source.source_coordinates),level in 1:nlevels
        index=_extrude_nonew_rect_grid_node(Int32(node),level,nlevels)
        point=cols[level,node]
        coordinates[1,index]=point[1]
        coordinates[2,index]=point[2]
        coordinates[3,index]=point[3]
    end
    arrays=ntuple(family->Matrix{Int32}(undef,(4,8,6,5)[family],
        catalog.cell_counts[family]),4)
    filled=zeros(Int,4)
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    sizehint!(edges,catalog.diagonal_count)
    orientation=Int(source.source_orientation)*Int(source.normal_direction)
    whole=_extrude_nonew_templates()[1]
    for source_cell in eachindex(source.source_cells),interval in 1:intervals
        original=source.source_cells[source_cell]
        corners=_extrude_nonew_corners(cols,original,interval)
        if interval==intervals
            state=catalog.top_states[source_cell]
            faces=(0x00,0x00,0x00,0x00,0x00,state)
            _extrude_nonew_add_diagonals!(edges,corners,faces)
            _extrude_nonew_certify(corners,edges,nothing,caller)
            center=_extrude_quadtri_centroid(corners)
            fan=_extrude_nonew_quad_strip_problem_fan(state,caller)
            _extrude_nonew_certify_template((corners...,center),fan,orientation,caller)
            index=Int32(column_count+source_cell)
            coordinates[1,index]=center[1]
            coordinates[2,index]=center[2]
            coordinates[3,index]=center[3]
            _extrude_nonew_b4_strip_emit!(arrays,filled,fan,original,
                interval,nlevels,index,caller)
        else
            _extrude_nonew_certify_template(corners,whole,orientation,caller)
            _extrude_nonew_b4_strip_emit!(arrays,filled,whole,original,
                interval,nlevels,Int32(0),caller)
        end
    end
    Tuple(filled)==catalog.cell_counts && length(edges)==catalog.diagonal_count ||
        _extrude_nonew_b4_strip_error(caller,
            "actual factory output differs from its checked capacities")
    blocks=ElementBlock[];sizehint!(blocks,4)
    for family in 1:4
        catalog.cell_counts[family]>0 || continue
        _extrude_make_positive!(coordinates,arrays[family],family+3)
        push!(blocks,ElementBlock(family+3,arrays[family]))
    end
    volume=MixedMesh(coordinates,blocks)
    sweep=_ExtrudeNoNewIndexedSweep(t,cols,volume)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source_tag,
        top_tag=top,lateral_tags=laterals,catalog=catalog)
end
