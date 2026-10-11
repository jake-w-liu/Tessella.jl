# Retained copied-Curve incidence keeps the actual lateral Quad products.
# Existing-corner terminal facets are a separate owned Region construction;
# their Quad coverage is certified without publishing Surface subdivisions.

function _extrude_nonew_retained_b4_plan(source::_ExtrudeNoNewRectGridSource,
        levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},caller;top_states=nothing)
    intervals=length(refs)
    intervals>0 && length(levels)==intervals+1 && levels[1]==0.0 && levels[end]==1.0 ||
        throw(ArgumentError("$caller: retained B4 requires normalized positive layers"))
    for index in eachindex(levels)
        isfinite(levels[index]) && (index==1 || levels[index]>levels[index-1]) ||
            throw(ArgumentError("$caller: retained B4 layers collapse or are not finite"))
    end
    _extrude_nonew_b4_strip_counts(source,intervals,caller)
    phase=_ExtrudeNoNewB4Free.SourcePhase.Phase(source,intervals;top_states)
    # Earlier intervals remain whole Hex macros. Only the terminal Region
    # interval may choose lateral facets; all lower cap interfaces stay whole.
    phase.lateral_forbidden[:,1:intervals-1].=true
    phase.cap_forbidden[:,1:intervals].=true
    _ExtrudeNoNewB4Free.SourcePhase.run!(phase)
    any(phase.problems) && throw(ArgumentError(
        "$caller: retained B4 terminal construction requires a new Region vertex"))
    cells=length(source.source_cells)
    indices=Matrix{UInt16}(undef,cells,intervals)
    families=(0,0,0,0);capacity=0
    templates=_extrude_nonew_templates()
    for cell in 1:cells,interval in 1:intervals
        mask=_ExtrudeNoNewB4Free.SourcePhase.mask(phase,cell,interval)
        interval==intervals || all(iszero,mask) || throw(ArgumentError(
            "$caller: retained B4 terminal choice changes an earlier whole macro"))
        index=_ExtrudeNoNewB4Free.FinalFactories.factory_index(mask)
        indices[cell,interval]=index
        template=templates[Int(index)]
        template.faces==mask || throw(ArgumentError(
            "$caller: retained B4 factory differs from its completed physical faces"))
        for position in 1:Int(template.ncells)
            part=template.cells[position];family=Int(part.msh)-3
            1<=family<=4 || throw(ArgumentError("$caller: retained B4 factory has an unsupported family"))
            width=(4,8,6,5)[family]
            all(vertex->1<=vertex<=8,part.nodes[1:width]) || throw(ArgumentError(
                "$caller: retained B4 factory introduces a new vertex"))
            families=Base.setindex(families,Base.checked_add(families[family],1),family)
            capacity=Base.checked_add(capacity,(4,6,5,5)[family])
        end
    end
    diagonal_count=Base.checked_add(count(!iszero,phase.lateral),count(!iszero,phase.cap))
    return _ExtrudeNoNewB4Free.Catalog(source,copy(phase.top_states),indices,
        NTuple{2,Int32}[],NTuple{6,UInt8}[],families,capacity,diagonal_count,true,levels,refs)
end

@inline _extrude_nonew_retained_b4_catalog(catalog)=
    catalog isa _ExtrudeNoNewB4Free.Catalog && catalog.recombined &&
    isempty(catalog.problem_positions) && isempty(catalog.problem_masks)

function _extrude_nonew_retained_b4_finish(m::GeoModel,t::Int,params,source_tag::Int,
        source_mesh,catalog::_ExtrudeNoNewB4Free.Catalog,cols,top,laterals,caller)
    params.recomb_laterals && _extrude_nonew_retained_b4_catalog(catalog) &&
        _extrude_nonew_unmerged_b4_incidence(m,t,caller) || throw(ArgumentError(
            "$caller: retained B4 Region constructor lacks its owned copied-Curve incidence"))
    _extrude_nonew_scope_contract(m,t,caller)
    spec=get(m.meshing.extrude_specs,(3,t),nothing)
    spec!==nothing || throw(ArgumentError("$caller: retained B4 Region lacks its current transform"))
    actual=_ExtrudeNoNewB4Free.emit(catalog,cols,spec,caller)
    size(actual.volume.coords,2)==length(catalog.source.source_coordinates)*length(catalog.levels) ||
        throw(ArgumentError("$caller: retained B4 Region introduces a body row"))
    sweep=_ExtrudeNoNewIndexedSweep(t,cols,actual.volume)
    return (sweep=sweep,edges=actual.edges,source_mesh=source_mesh,source_tag=source_tag,
        top_tag=top,lateral_tags=laterals,catalog=catalog)
end
