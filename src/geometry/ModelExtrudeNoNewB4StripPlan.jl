# Recombined all-boundary strips use whole earlier Hex8 macros and one
# certified actual terminal center per retained Quad. No layer pick matrix.
struct _ExtrudeNoNewB4StripCatalog <: _ExtrudeNoNewDynamicGridCatalog
    source::_ExtrudeNoNewRectGridSource
    top_states::Vector{UInt8}
    cell_counts::NTuple{4,Int}
    face_capacity::Int
    diagonal_count::Int
    recombined::Bool
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@noinline function _extrude_nonew_b4_strip_plan_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts all-boundary quadrangle strip $reason"))
end

@inline _extrude_nonew_b4_strip_face_capacity(catalog::_ExtrudeNoNewB4StripCatalog)=
    catalog.face_capacity

# Original actual source-column ordinal is one global logical rank. Raw Gmsh
# pointers and public mesh tags are deliberately not used as rank surrogates.
@inline function _extrude_nonew_b4_strip_top_state(cell::NTuple{4,Int32})
    smallest=cell[1]
    position=1
    for corner in 2:4
        if cell[corner]<smallest
            smallest=cell[corner]
            position=corner
        end
    end
    return isodd(position) ? UInt8(1) : UInt8(2)
end

# Count-only preflight precedes the top-state allocation. No emitted centroid,
# factory record or per-layer selection matrix is stored in this planner.
function _extrude_nonew_b4_strip_counts(source::_ExtrudeNoNewRectGridSource,
        intervals::Int,caller)
    a,b=source.grid_shape
    min(a,b)==1 && max(a,b)>=5 ||
        _extrude_nonew_b4_strip_plan_error(caller,"requires one cell across and at least five along")
    intervals>0 || _extrude_nonew_b4_strip_plan_error(caller,"requires positive intervals")
    counts=try
        cells=Base.checked_mul(a,b)
        vertices=Base.checked_mul(Base.checked_add(a,1),Base.checked_add(b,1))
        primary_nodes=Base.checked_mul(vertices,Base.checked_add(intervals,1))
        nodes=Base.checked_add(primary_nodes,cells)
        families=(Base.checked_mul(2,cells),Base.checked_mul(cells,intervals-1),
            0,Base.checked_mul(5,cells))
        total_cells=Base.checked_mul(cells,Base.checked_add(intervals,6))
        faces=Base.checked_add(Base.checked_mul(Base.checked_add(Base.checked_mul(4,cells),1),intervals),
            Base.checked_mul(15,cells))
        (cells,vertices,nodes,families,total_cells,faces)
    catch err
        err isa OverflowError || rethrow()
        _extrude_nonew_b4_strip_plan_error(caller,"size arithmetic overflows Int")
    end
    cells,vertices,nodes,families,total_cells,faces=counts
    nodes<=_EXTRUDE_NONEW_MAX_NODES && nodes<=typemax(Int32) ||
        _extrude_nonew_b4_strip_plan_error(caller,"exceeds the checked node or Int32 limit")
    total_cells<=typemax(Int32) ||
        _extrude_nonew_b4_strip_plan_error(caller,"output cell count exceeds Int32")
    cells==length(source.source_cells) && vertices==length(source.source_coordinates) &&
        length(source.category)==cells && length(source.boundary_masks)==cells ||
        _extrude_nonew_b4_strip_plan_error(caller,"dimensions differ from the certified source graph")
    for cell in 1:cells
        source.category[cell]==0x04 && source.boundary_masks[cell]==0x0f ||
            _extrude_nonew_b4_strip_plan_error(caller,"requires four actual boundary vertices in every source cell")
    end
    return families,faces,cells
end

function _extrude_nonew_b4_strip_plan(source::_ExtrudeNoNewRectGridSource,
        levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},recombined::Bool,caller)
    recombined || _extrude_nonew_b4_strip_plan_error(caller,
        "free lateral faces require the pending all-boundary propagation planner")
    intervals=length(refs)
    intervals>0 && length(levels)==intervals+1 ||
        _extrude_nonew_b4_strip_plan_error(caller,"requires one more level than positive intervals")
    levels[1]==0.0 && levels[end]==1.0 ||
        _extrude_nonew_b4_strip_plan_error(caller,"requires normalized layer levels")
    for level in eachindex(levels)
        isfinite(levels[level]) ||
            _extrude_nonew_b4_strip_plan_error(caller,"requires finite layer levels")
        level==1 || levels[level]>levels[level-1] ||
            _extrude_nonew_b4_strip_plan_error(caller,"layer levels collapse in Float64")
    end
    families,faces,cells=_extrude_nonew_b4_strip_counts(source,intervals,caller)
    top_states=Vector{UInt8}(undef,cells)
    for cell in 1:cells
        top_states[cell]=_extrude_nonew_b4_strip_top_state(source.source_cells[cell])
    end
    return _ExtrudeNoNewB4StripCatalog(source,top_states,families,faces,cells,
        true,levels,refs)
end
