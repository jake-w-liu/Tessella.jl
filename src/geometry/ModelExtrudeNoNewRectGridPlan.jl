# Constructive normal-translation NoNew phases for a rectangular B3/B2/B0
# source. Physical source edges are assigned once; per-column cap states carry
# across every layer group. The immutable existing-corner factories are reused.
struct _ExtrudeNoNewRectGridCatalog
    source::_ExtrudeNoNewRectGridSource
    top_states::Vector{UInt8}
    template_indices::Matrix{UInt16}
    cell_counts::NTuple{4,Int}
    face_capacity::Int
    diagonal_count::Int
    recombined::Bool
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@noinline function _extrude_nonew_rect_grid_plan_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts rectangular grid source $reason"))
end

@inline _extrude_nonew_rect_grid_face_capacity(catalog::_ExtrudeNoNewRectGridCatalog)=
    catalog.face_capacity

@inline _extrude_nonew_rect_grid_reverse(state::UInt8)=
    state==0 ? UInt8(0) : UInt8(3-state)

@inline function _extrude_nonew_rect_grid_local_state(source,states,cell,side)
    edge=Int(source.edge_indices[cell][side])
    state=states[edge]
    return source.source_cells[cell][side]==source.edges[edge][1] ? state :
        _extrude_nonew_rect_grid_reverse(state)
end

@inline function _extrude_nonew_rect_grid_laterals(source,states,cell)
    return (_extrude_nonew_rect_grid_local_state(source,states,cell,1),
            _extrude_nonew_rect_grid_local_state(source,states,cell,2),
            _extrude_nonew_rect_grid_local_state(source,states,cell,3),
            _extrude_nonew_rect_grid_local_state(source,states,cell,4))
end

# Gmsh ranks TOP MVertex pointers. A consistent original-source ordinal replaces
# the unobservable pointer order; mesh tags and coordinates never define rank.
# B2/B3 choose only non-boundary corners, B0 chooses among all four corners.
function _extrude_nonew_rect_grid_top_state(cell,mask::UInt8,category::UInt8,caller)
    category in (0x00,0x02,0x03) && mask<=0x0f && count_ones(mask)==category ||
        _extrude_nonew_rect_grid_plan_error(caller,"has an unsupported boundary category")
    category!=2 || mask in (0x03,0x06,0x0c,0x09) ||
        _extrude_nonew_rect_grid_plan_error(caller,"opposite-boundary quads require the boundary-category planner")
    smallest=typemax(Int32)
    position=0
    for corner in 1:4
        if category==0 || (mask & (UInt8(1)<<(corner-1)))==0
            node=cell[corner]
            if node<smallest
                smallest=node
                position=corner
            end
        end
    end
    position!=0 || _extrude_nonew_rect_grid_plan_error(caller,"has no eligible top corner")
    return isodd(position) ? UInt8(1) : UInt8(2)
end

# Cycle (u_lower,v_lower,v_upper,u_upper), with source edges stored minmax.
# One assignment supplies both incident cells, irrespective of their order.
function _extrude_nonew_rect_grid_physical_states!(states,source,recombined,terminal,caller)
    for edge in eachindex(source.edges)
        u,v=source.edges[edge]
        slot=source.edge_curve[edge]
        if source.edge_cells[edge][2]==0
            1<=slot<=4 && source.boundary_vertices[u] && source.boundary_vertices[v] ||
                _extrude_nonew_rect_grid_plan_error(caller,"has an unclassified exterior source edge")
            direction=source.edge_curve_direction[edge]
            direction in (-1,1) ||
                _extrude_nonew_rect_grid_plan_error(caller,"has no native exterior Curve direction")
            states[edge]=recombined ? UInt8(0) : direction>0 ? UInt8(1) : UInt8(2)
        else
            slot==0 || _extrude_nonew_rect_grid_plan_error(caller,"classifies an internal source edge on a Curve")
            ub=source.boundary_vertices[u]
            vb=source.boundary_vertices[v]
            if ub && vb
                _extrude_nonew_rect_grid_plan_error(caller,"an internal edge joins two boundary vertices; the all-boundary phase is required")
            elseif recombined && !terminal
                states[edge]=0
            elseif ub!=vb
                states[edge]=ub ? UInt8(1) : UInt8(2)
            elseif terminal
                # Least original ordinal in the TOP plane joins the OTHER
                # lower endpoint. The opposite direction would break cap rank.
                states[edge]=u<v ? UInt8(2) : UInt8(1)
            else
                states[edge]=0
            end
        end
    end
    return nothing
end

@inline function _extrude_nonew_rect_grid_b2_cap(mask::UInt8,laterals,lower::UInt8,caller)
    lower!=0 && return lower
    side=mask==0x03 ? 1 : mask==0x06 ? 2 : mask==0x0c ? 3 : mask==0x09 ? 4 : 0
    side!=0 || _extrude_nonew_rect_grid_plan_error(caller,"has no adjacent boundary pair")
    state=laterals[side]
    state in (0x01,0x02) ||
        _extrude_nonew_rect_grid_plan_error(caller,"a free boundary face lacks its diagonal")
    anchor=state==1 ? mod1(side+1,4) : side
    return isodd(anchor) ? UInt8(1) : UInt8(2)
end

@inline function _extrude_nonew_rect_grid_mask_code(faces,lower::UInt8,upper::UInt8)
    return 1+Int(faces[1])+3Int(faces[2])+9Int(faces[3])+27Int(faces[4])+
        81Int(lower)+243Int(upper)
end

function _extrude_nonew_rect_grid_plan(source::_ExtrudeNoNewRectGridSource,
        levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},recombined::Bool,caller)
    intervals=length(refs)
    intervals>0 && length(levels)==intervals+1 ||
        _extrude_nonew_rect_grid_plan_error(caller,"requires one more level than positive intervals")
    levels[1]==0.0 && levels[end]==1.0 && all(isfinite,levels) ||
        _extrude_nonew_rect_grid_plan_error(caller,"requires finite normalized layer levels")
    for level in 2:length(levels)
        levels[level]>levels[level-1] ||
            _extrude_nonew_rect_grid_plan_error(caller,"layer levels collapse in Float64")
    end
    a,b=source.grid_shape
    a>=2 && b>=2 || _extrude_nonew_rect_grid_plan_error(caller,"requires at least two cells in each grid direction")
    cells=Base.checked_mul(a,b)
    vertices=Base.checked_mul(Base.checked_add(a,1),Base.checked_add(b,1))
    boundary=Base.checked_mul(2,Base.checked_add(a,b))
    cells==length(source.source_cells) && vertices==length(source.source_coordinates) &&
        length(source.category)==cells && length(source.boundary_masks)==cells &&
        length(source.edge_indices)==cells ||
        _extrude_nonew_rect_grid_plan_error(caller,"source dimensions differ from the certified rectangular graph")
    nodes=Base.checked_mul(vertices,Base.checked_add(intervals,1))
    nodes<=_EXTRUDE_NONEW_MAX_NODES && nodes<=typemax(Int32) ||
        _extrude_nonew_rect_grid_plan_error(caller,"exceeds the node or Int32 limit")
    Base.checked_mul(6,Base.checked_mul(cells,intervals))<=typemax(Int32) ||
        _extrude_nonew_rect_grid_plan_error(caller,"output cell count exceeds Int32")

    top_states=Vector{UInt8}(undef,cells)
    for cell in 1:cells
        top_states[cell]=_extrude_nonew_rect_grid_top_state(source.source_cells[cell],
            source.boundary_masks[cell],source.category[cell],caller)
    end
    canonical=_extrude_nonew_quad_strip_mask_factories()
    templates=_extrude_nonew_templates()
    picks=Matrix{UInt16}(undef,cells,intervals)
    states=Vector{UInt8}(undef,length(source.edges))
    lower_states=zeros(UInt8,cells)
    counts=zeros(Int,4)
    incidences=0
    diagonals=0
    for interval in 1:intervals
        terminal=interval==intervals
        _extrude_nonew_rect_grid_physical_states!(states,source,recombined,terminal,caller)
        for state in states
            diagonals=Base.checked_add(diagonals,Int(state!=0))
        end
        for cell in 1:cells
            laterals=_extrude_nonew_rect_grid_laterals(source,states,cell)
            lower=lower_states[cell]
            category=source.category[cell]
            upper=if terminal
                top_states[cell]
            elseif recombined || category==0
                UInt8(0)
            elseif category==3
                top_states[cell]
            else
                _extrude_nonew_rect_grid_b2_cap(source.boundary_masks[cell],laterals,lower,caller)
            end
            packed=canonical[_extrude_nonew_rect_grid_mask_code(laterals,lower,upper)]
            packed!=0 || _extrude_nonew_rect_grid_plan_error(caller,
                "constructed face states have no existing-corner factory; the boundary-category planner is required")
            template=templates[Int(packed)]
            picks[cell,interval]=packed
            for position in 1:Int(template.ncells)
                family=Int(template.cells[position].msh)-3
                1<=family<=4 || _extrude_nonew_rect_grid_plan_error(caller,"factory has an unsupported cell family")
                counts[family]=Base.checked_add(counts[family],1)
                incidences=Base.checked_add(incidences,family==1 ? 4 : family==2 ? 6 : 5)
            end
            diagonals=Base.checked_add(diagonals,Int(upper!=0))
            lower_states[cell]=upper
        end
    end
    typed_boundary=Base.checked_add(Base.checked_mul(recombined ? 1 : 2,
        Base.checked_mul(boundary,intervals)),Base.checked_mul(3,cells))
    total_faces=Base.checked_add(incidences,typed_boundary)
    iseven(total_faces) || _extrude_nonew_rect_grid_plan_error(caller,"has inconsistent typed face-count parity")
    cell_counts=(counts[1],counts[2],counts[3],counts[4])
    if recombined
        expected=(Base.checked_sub(Base.checked_mul(4,cells),Base.checked_mul(2,boundary)),
            Base.checked_mul(cells,intervals-1),0,Base.checked_add(cells,boundary))
        cell_counts==expected || _extrude_nonew_rect_grid_plan_error(caller,
            "terminal B3/B2/B0 factories differ from the certified family counts")
    end
    return _ExtrudeNoNewRectGridCatalog(source,top_states,picks,cell_counts,
        total_faces÷2,diagonals,recombined,levels,refs)
end

@inline function _extrude_nonew_rect_grid_template(catalog::_ExtrudeNoNewRectGridCatalog,cell,interval)
    return _extrude_nonew_templates()[Int(catalog.template_indices[cell,interval])]
end

@inline _extrude_nonew_rect_grid_faces(catalog::_ExtrudeNoNewRectGridCatalog,cell,interval)=
    _extrude_nonew_rect_grid_template(catalog,cell,interval).faces

@inline function _extrude_nonew_rect_grid_edge_state(catalog::_ExtrudeNoNewRectGridCatalog,edge,interval)
    source=catalog.source
    cell=Int(source.edge_cells[edge][1])
    side=Int(source.edge_sides[edge][1])
    state=_extrude_nonew_rect_grid_faces(catalog,cell,interval)[side]
    return source.source_cells[cell][side]==source.edges[edge][1] ? state :
        _extrude_nonew_rect_grid_reverse(state)
end
