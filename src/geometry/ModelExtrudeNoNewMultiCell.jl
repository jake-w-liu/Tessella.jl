# The unstructured multi-cell Source family for QuadTriNoNewVerts — the
# general per-column prism planner upstream runs for every swept triangle.
# Each source cell sweeps its own prism column through the existing-corner
# template relation (`_EXTRUDE_NONEW_PRISM_TEMPLATES`): a retained prism
# when all three lateral quadrangles stay whole, a pyramid+tetrahedron fan
# when exactly one stays whole, and a three-tetrahedron split when all
# three carry diagonals. The two twisted all-split combinations resolve by
# the same diagonal-swap fixpoint `SubdivideExtrudedMesh` phase2 runs
# upstream, so no subdivision ever introduces a vertex.
#
# Whole/split is decided per Source edge, so both incident columns emit
# the identical facets. A boundary edge's quadrangle stays whole exactly
# when the shared wall its generatrix produces stays whole upstream —
# every owner keeps quads, or an unrecombined partner forces the split —
# while a free lateral always takes the chain-preferred split
# (chain-first bottom to chain-second top). An interior edge stays whole
# exactly when it joins two interior cells — Source cells carrying no
# boundary edge — whose two endpoints share one boundary class: the
# oracle keeps the centroid-to-centroid wall of a two-centroid fan
# (interior vertices) and the corner-to-corner chord of an L-shape
# Source (boundary vertices) whole, and splits the corner-to-centroid
# edges between them. A fully surrounded interior cell then retains a
# whole prism, a cell with one such edge fans to a pyramid on the shared
# quad base — the oracle's subregion-pair emission — and the
# two-candidate mixed case never arises on observed oracle meshes. The
# shared wall itself transports through the bound Surface part — each
# segment quadrangle of the earlier owner's emitted wall is reproduced
# facet-for-facet — so multi-segment generatrices weld bitwise across
# owner boundaries, unlike the generic catalog's single
# diagonal-per-layer contract.

struct _ExtrudeNoNewMultiCellCatalog
    source_cells::Vector{NTuple{3,Int32}}
    boundary_curves::Vector{Int}
    boundary_chains::Vector{Vector{Int32}}
    boundary_widths::Vector{Int}
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@noinline _extrude_nonew_multi_cell_error(caller)=
    throw(ArgumentError(
        "$caller: QuadTriNoNewVerts multiple source cells require the " *
        "boundary-category planner"))

# Pre-mesh estimate that a four-sided Source resolves to this family: it
# is not recombined (a recombined Source stays on a bounded quadrangle
# catalog or the generic single-quad route) and its boundary produces more
# than the two-triangle catalog's single strip — either the Source carries
# no transfinite record at all (the native mesher triangulates it freely)
# or a transfinite record carries a generatrix beyond two nodes. Used only
# for shared-wall owner classification; the actual mesh decides admission.
function _extrude_nonew_multi_cell_candidate(m::GeoModel,source::Int,
        sides::Int,caller::AbstractString)
    sides==4 || return false
    _model_surface_recombined(m,source) && return false
    haskey(m.meshing.transfinite_surfaces,source) || return true
    for curve in _model_projection_surface_curves(m,source)
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control===nothing && return false
        _flexible_transfinite_nodes(m,control.num_nodes,curve,caller)>2 &&
            return true
    end
    return false
end

# Whether `source`'s Region resolves to this family for the shared-wall
# split/whole vote. The two-triangle catalog is still its own bounded
# family — its pre-mesh claim excludes the certified 2-node strip — while
# every other bounded catalog is either recombined or transfinite.
function _extrude_nonew_multi_cell_source(m::GeoModel,source::Int,
        caller::AbstractString)
    _extrude_nonew_multi_cell_candidate(m,source,4,caller) || return false
    _extrude_nonew_two_tri_candidate(m,source,4,caller) && return false
    return true
end

# The actual Source mesh decides: two or more linear triangles sweep
# per-cell prism columns. A single triangle stays on the isolated-triangle
# route; any quadrangle or higher-degree cell keeps the bounded catalogs
# or the boundary-category diagnostic.
function _extrude_nonew_multi_cells(source_mesh,caller)
    cells=NTuple{3,Int32}[]
    if source_mesh isa Mesh
        for column in axes(source_mesh.tris,2)
            push!(cells,ntuple(k->source_mesh.tris[k,column],3))
        end
    else
        for block in source_mesh.blocks
            block isa ElementBlock || continue
            msh_dimension(block.msh)==2 || continue
            block.msh==2 || return nothing
            for column in axes(block.nodes,2)
                push!(cells,ntuple(k->block.nodes[k,column],3))
            end
        end
    end
    length(cells)>=2 || return nothing
    total=nnodes(source_mesh)
    all(cell->all(node->Int32(1)<=node<=Int32(total),cell),cells) ||
        return nothing
    return cells
end

# The emitted dim-2 cells of a bound wall part, as position tuples — the
# same read `_extrude_nonew_wall_diagonals` performs for its single-
# segment contract, factored so a multi-segment wall resolves each segment
# quadrangle independently.
function _extrude_nonew_wall_part_cells(part,s::Int,caller::AbstractString)
    cells=Vector{NTuple{3,Float64}}[]
    if part isa Mesh
        for column in axes(part.tris,2)
            push!(cells,[ntuple(k->part.coords[k,part.tris[j,column]],3)
                         for j in 1:3])
        end
    else
        for block in part.blocks
            block isa ElementBlock && msh_dimension(block.msh)==2 ||
                continue
            block.msh in (2,3) ||
                _extrude_nonew_shared_planner_error(caller,s)
            width=size(block.nodes,1)
            for column in axes(block.nodes,2)
                push!(cells,[ntuple(k->part.coords[k,block.nodes[j,column]],3)
                             for j in 1:width])
            end
        end
    end
    return cells
end

# One wall quadrangle's bound emission: `nothing` for a whole quadrangle,
# the shared diagonal position pair for a split pair, else the loud
# diagnostic — a bound wall decomposes each of its quadrangles into
# exactly one quad or two triangles sharing its diagonal.
function _extrude_nonew_wall_quad_state(member,s::Int,caller::AbstractString)
    if length(member)==1 && length(member[1])==4
        return nothing
    end
    if length(member)==2 && all(cell->length(cell)==3,member)
        shared=intersect(Set(member[1]),Set(member[2]))
        length(shared)==2 || _extrude_nonew_shared_planner_error(caller,s)
        a,b=collect(shared)
        return _extrude_ekey(a,b)
    end
    _extrude_nonew_shared_planner_error(caller,s)
end

function _extrude_nonew_multi_cell_finish(m::GeoModel,t::Int,
        params,source::Int,source_mesh,spec,
        cells::Vector{NTuple{3,Int32}},top::Int,laterals,
        shared::Dict{Int,Vector{Int}},caller::AbstractString;
        _shared_parts=nothing)
    levels,refs=_extrude_nonew_levels(params,caller;
        source_nodes=nnodes(source_mesh),extra_nodes=0,
        cells_per_interval=Base.checked_mul(3,length(cells)))
    intervals=length(refs)
    ncells=length(cells)
    cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
    all(p->all(isfinite,p),cols) && allunique(cols) ||
        _extrude_nonew_multi_cell_error(caller)
    position_node=Dict{NTuple{3,Float64},Int32}(
        cols[1,node]=>Int32(node) for node in axes(cols,2))
    # Each generatrix's evaluated chain fixes both the boundary-edge
    # ownership map and the preferred free-lateral diagonal direction
    # (chain-first bottom row to chain-second top row).
    curves=_model_projection_surface_curves(m,source)
    chains=Dict{Int,Vector{NTuple{3,Float64}}}()
    seg_owner=Dict{NTuple{2,NTuple{3,Float64}},Tuple{Int,Int}}()
    boundary_vertex=Set{NTuple{3,Float64}}()
    for curve in curves
        c=abs(curve)
        chain=_extrude_curve_nodes(m,c,caller)
        chains[c]=chain
        a,b=m.curves[c]
        closed=a==b
        nseg=closed ? length(chain) : length(chain)-1
        for seg in 1:nseg
            seg2=closed ? mod1(seg+1,length(chain)) : seg+1
            pa,pb=chain[seg],chain[seg2]
            push!(boundary_vertex,pa)
            push!(boundary_vertex,pb)
            key=_extrude_ekey(pa,pb)
            haskey(seg_owner,key) &&
                _extrude_nonew_multi_cell_error(caller)
            seg_owner[key]=(c,seg)
        end
    end
    # Retain each boundary curve's native chain as source-row indices, the
    # same projection provenance the grid catalogs carry — the volume and
    # surface projections rebuild entity Curve columns through them.
    boundary_curves=Int[]
    boundary_chains=Vector{Int32}[]
    boundary_widths=Int[]
    for curve in curves
        c=abs(curve)
        chain=chains[c]
        indices=Int32[get(position_node,p,Int32(0)) for p in chain]
        all(>(0),indices) || _extrude_nonew_multi_cell_error(caller)
        push!(boundary_curves,c)
        push!(boundary_chains,indices)
        push!(boundary_widths,length(chain))
    end
    lateral_of=Dict{Int,Int}()
    for s in laterals
        s_link=get(m.meshing.extrude_sources,(2,s),nothing)
        (s_link===nothing || s_link[1]!=1) && continue
        generatrix=abs(s_link[2])
        haskey(lateral_of,generatrix) &&
            _extrude_nonew_multi_cell_error(caller)
        lateral_of[generatrix]=s
    end
    # Source-edge incidence: an edge carried by two cells is interior; an
    # edge carried by one must coincide with a generatrix chain segment.
    edge_inc=Dict{NTuple{2,Int32},Int}()
    edge_cells=Dict{NTuple{2,Int32},Vector{Int}}()
    cell_edges=Vector{NTuple{3,NTuple{2,Int32}}}(undef,ncells)
    for (ci,cell) in enumerate(cells)
        es=ntuple(3) do side
            a,b=cell[side],cell[mod1(side+1,3)]
            a<b ? (a,b) : (b,a)
        end
        cell_edges[ci]=es
        for e in es
            edge_inc[e]=get(edge_inc,e,0)+1
            push!(get!(edge_cells,e,Int[]),ci)
        end
    end
    any(count->count>2,values(edge_inc)) &&
        _extrude_nonew_multi_cell_error(caller)
    # A boundary cell is any cell with a boundary (single-incidence) edge —
    # the oracle keeps an interior edge's wall whole exactly when neither
    # incident cell is a boundary cell.
    boundary_cell=falses(ncells)
    for (e,count) in edge_inc
        count==1 || continue
        boundary_cell[only(edge_cells[e])]=true
    end
    edge_meta=Dict{NTuple{2,Int32},Tuple{Int,Int,Int,Int}}()
    n_boundary=0
    for (e,count) in edge_inc
        count==2 && continue
        pa,pb=cols[1,e[1]],cols[1,e[2]]
        hit=get(seg_owner,_extrude_ekey(pa,pb),nothing)
        hit===nothing && _extrude_nonew_multi_cell_error(caller)
        generatrix,seg=hit
        chain=chains[generatrix]
        a,b=m.curves[generatrix]
        seg2=a==b ? mod1(seg+1,length(chain)) : seg+1
        # `first`/`second` are the edge endpoints in chain direction, so
        # the preferred diagonal runs first@layer to second@layer+1.
        first=cols[1,e[1]]==chain[seg] ? e[1] : e[2]
        second=cols[1,e[1]]==chain[seg] ? e[2] : e[1]
        cols[1,first]==chain[seg] && cols[1,second]==chain[seg2] ||
            _extrude_nonew_multi_cell_error(caller)
        edge_meta[e]=(generatrix,seg,first,second)
        n_boundary+=1
    end
    # The triangulation's boundary is exactly the generatrix chain
    # complex — every chain segment is one source boundary edge.
    n_boundary==length(seg_owner) ||
        _extrude_nonew_multi_cell_error(caller)
    # Per-edge, per-layer wall state: `nothing` keeps the quadrangle
    # whole, an edge key fixes the emitted diagonal. `edge_forced` marks
    # the bound-wall contract — those diagonals reproduce an earlier
    # owner's emission and are swap-protected in the twist solve.
    edge_state=Dict{NTuple{2,Int32},
        Vector{Union{Nothing,NTuple{2,NTuple{3,Float64}}}}}()
    edge_forced=Set{NTuple{2,Int32}}()
    part_cells=Dict{Int,Vector{Vector{NTuple{3,Float64}}}}()
    for e in sort!(collect(keys(edge_inc)))
        count=edge_inc[e]
        diag=Vector{Union{Nothing,NTuple{2,NTuple{3,Float64}}}}(
            undef,intervals)
        if count==2
            ci,cj=edge_cells[e][1],edge_cells[e][2]
            # RecombLaterals keeps every lateral quadrangle whole — the
            # oracle's emission is one retained prism per column layer
            # (b_multi_recomblat: 28 prisms, no interior splits) — so
            # interior edges stay whole too, not only between two
            # interior source vertices. Without RecombLaterals the oracle
            # keeps an interior wall whole when both incident cells are
            # interior (no boundary edge) and the edge's endpoints carry
            # the same boundary class — interior-interior and
            # corner-corner edges qualify, the corner-to-centroid edge
            # does not.
            if params.recomb_laterals
                fill!(diag,nothing)
            else
                e1_boundary=cols[1,e[1]] in boundary_vertex
                e2_boundary=cols[1,e[2]] in boundary_vertex
                if !boundary_cell[ci] && !boundary_cell[cj] &&
                        e1_boundary==e2_boundary
                    fill!(diag,nothing)
                else
                    for layer in 1:intervals
                        diag[layer]=_extrude_ekey(cols[layer,e[1]],
                                                  cols[layer+1,e[2]])
                    end
                end
            end
            edge_state[e]=diag
            continue
        end
        generatrix,_,first,second=edge_meta[e]
        s=get(lateral_of,generatrix,0)
        s==0 && _extrude_nonew_multi_cell_error(caller)
        part=haskey(shared,s) && _shared_parts!==nothing ?
            get(_shared_parts,s,nothing) : nothing
        if part!==nothing
            member_cells=get!(part_cells,s,
                _extrude_nonew_wall_part_cells(part,s,caller))
            for layer in 1:intervals
                quad=Set((cols[layer,e[1]],cols[layer,e[2]],
                          cols[layer+1,e[1]],cols[layer+1,e[2]]))
                member=[cell for cell in member_cells
                        if all(p->p in quad,cell)]
                diag[layer]=_extrude_nonew_wall_quad_state(member,s,caller)
            end
            push!(edge_forced,e)
        elseif params.recomb_laterals ||
                (haskey(shared,s) &&
                 !_extrude_nonew_wall_split(m,s,shared[s],caller))
            fill!(diag,nothing)
        else
            for layer in 1:intervals
                diag[layer]=_extrude_ekey(cols[layer,first],
                                          cols[layer+1,second])
            end
        end
        edge_state[e]=diag
    end
    # A column's local template state per lateral face: 0 keeps the face,
    # 1 selects the first-to-top diagonal, 2 the second-to-top, in the
    # cell's own vertex order (`_EXTRUDE_NONEW_PRISM_FACES` lateral k is
    # source edge (k,k+1)).
    function column_state(ci::Int,layer::Int)
        cell=cells[ci]
        es=cell_edges[ci]
        return ntuple(3) do side
            e=es[side]
            d=edge_state[e][layer]
            d===nothing && return UInt8(0)
            p,q=cell[side],cell[mod1(side+1,3)]
            d==_extrude_ekey(cols[layer,p],cols[layer+1,q]) &&
                return UInt8(1)
            d==_extrude_ekey(cols[layer,q],cols[layer+1,p]) &&
                return UInt8(2)
            _extrude_nonew_multi_cell_error(caller)
        end
    end
    # The fixpoint upstream's `SubdivideExtrudedMesh` phase2 generalizes
    # to: any column whose three lateral states have no template is
    # invalid — the two twisted all-split combinations, plus the
    # one-whole combinations whose two split diagonals do not meet at a
    # shared apex corner. Flip the diagonal that lands this column on a
    # valid template without invalidating a currently-valid incident
    # column (a shared interior edge's diagonal toggles in both cells'
    # local frames); when no such flip exists, the first unprotected
    # nonzero face still perturbs the state. A flipped (or
    # contract-forced) face is swap-protected so the fixpoint terminates.
    templates=_extrude_nonew_prism_templates()
    template_faces=Set(template.faces for template in templates)
    locked=Set{Tuple{NTuple{2,Int32},Int}}()
    for e in edge_forced,layer in 1:intervals
        push!(locked,(e,layer))
    end
    flip_edge!(ci::Int,side::Int,layer::Int)=begin
        e=cell_edges[ci][side]
        p,q=cells[ci][side],cells[ci][mod1(side+1,3)]
        first_key=_extrude_ekey(cols[layer,p],cols[layer+1,q])
        edge_state[e][layer]=edge_state[e][layer]==first_key ?
            _extrude_ekey(cols[layer,q],cols[layer+1,p]) : first_key
        push!(locked,(e,layer))
        nothing
    end
    toggle(state::UInt8)=state==0x01 ? UInt8(2) : UInt8(1)
    while true
        invalid=false
        for ci in 1:ncells,layer in 1:intervals
            state=column_state(ci,layer)
            state in template_faces && continue
            invalid=true
            flipped=false
            fallback=0
            for side in 1:3
                e=cell_edges[ci][side]
                edge_state[e][layer]===nothing && continue
                (e,layer) in locked && continue
                fallback==0 && (fallback=side)
                trial=ntuple(k->k==side ? toggle(state[k]) : state[k],3)
                trial in template_faces || continue
                # The shared incident column's state toggles at its own
                # side — skip flips that would turn a valid neighbor
                # invalid when another side already fixes this column.
                hurts=false
                for cj in 1:ncells
                    cj==ci && continue
                    sj=findfirst(isequal(e),cell_edges[cj])
                    sj===nothing && continue
                    sj_state=column_state(cj,layer)
                    sj_state in template_faces || continue
                    trial_j=ntuple(k->k==sj ? toggle(sj_state[k]) :
                                   sj_state[k],3)
                    trial_j in template_faces || (hurts=true)
                end
                hurts && continue
                flip_edge!(ci,side,layer)
                flipped=true
                break
            end
            if !flipped && fallback!=0
                flip_edge!(ci,fallback,layer)
                flipped=true
            end
            flipped || _extrude_nonew_multi_cell_error(caller)
        end
        invalid || break
    end
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    for diag in values(edge_state),d in diag
        d===nothing || push!(edges,d)
    end
    sweep=_ExtrudeVolumeSweep(t,true,cols,
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false])
    for (ci,cell) in enumerate(cells),layer in 1:intervals
        wanted=column_state(ci,layer)
        index=findfirst(template->template.faces==wanted,templates)
        index===nothing && _extrude_nonew_multi_cell_error(caller)
        template=templates[index]
        v=_extrude_nonew_corners(cols,cell,layer)
        _extrude_nonew_certify(v,edges,template,caller)
        _extrude_nonew_emit_template!(sweep,v,template)
    end
    catalog=_ExtrudeNoNewMultiCellCatalog(cells,boundary_curves,
        boundary_chains,boundary_widths,levels,refs)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,
        source_tag=source,top_tag=top,lateral_tags=laterals,catalog=catalog)
end
