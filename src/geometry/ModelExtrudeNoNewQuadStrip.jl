# Two actual native quadrangles share one source edge. All six source vertices
# are boundary vertices; free caps need a joined chain rather than pivot rules.
struct _ExtrudeNoNewStripProblemFan
    ncells::UInt8
    cells::NTuple{7,_ExtrudeNoNewCell}
end

struct _ExtrudeNoNewQuadStripCatalog
    source_cells::NTuple{2,NTuple{4,Int32}}
    source_coordinates::NTuple{6,NTuple{3,Float64}}
    boundary_cycle::NTuple{6,Int32}
    boundary_curves::NTuple{4,Int}
    boundary_chains::NTuple{4,NTuple{3,Int32}}
    boundary_widths::NTuple{4,UInt8}
    edges::NTuple{7,NTuple{2,Int32}}
    edge_indices::NTuple{2,NTuple{4,UInt8}}
    shared_sides::NTuple{2,UInt8}
    preferred::NTuple{2,NTuple{4,UInt8}}
    top_states::NTuple{2,UInt8}
    template_indices::Matrix{UInt16}
    cell_counts::NTuple{4,Int}
    face_capacity::Int
    diagonal_count::Int
    normal_axis::Int8
    normal_direction::Int8
    source_orientation::Int8
    recombined::Bool
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@noinline function _extrude_nonew_quad_strip_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts two-quadrangle strip source $reason"))
end

@inline _extrude_nonew_quad_strip_face_capacity(catalog::_ExtrudeNoNewQuadStripCatalog)=
    catalog.face_capacity

function _extrude_nonew_quad_strip_axis(spec,caller)
    spec.type===:translate || _extrude_nonew_quad_strip_error(caller,
        "requires axis-normal translation; other transforms require the source-grid planner")
    all(isfinite,spec.T) && count(!iszero,spec.T)==1 ||
        _extrude_nonew_quad_strip_error(caller,"requires one finite nonzero translation component")
    axis=spec.T[1]!=0 ? 1 : spec.T[2]!=0 ? 2 : 3
    return axis,spec.T[axis]>0 ? 1 : -1
end

@inline _extrude_nonew_quad_strip_pair(first::Int,second::Int)=1+first+3second
@inline _extrude_nonew_quad_strip_reverse(state::UInt8)=state==0 ? UInt8(0) : UInt8(3-state)
@inline _extrude_nonew_quad_strip_relation_slot(lower,upper,shared)=1+Int(lower)+3Int(upper)+9Int(shared)

struct _ExtrudeNoNewStripTransition
    first::UInt16
    second::UInt16
    changes::UInt8
    alignment::UInt8
end

# Resolve factories only after the complete face mask is known. Distinct
# admissible masks retain adjustable exterior diagonal preferences.
function _extrude_nonew_quad_strip_mask_factories()
    templates=_extrude_nonew_templates()
    length(templates)<=typemax(UInt16) || throw(ErrorException("QuadTriNoNewVerts template index overflow"))
    canonical=zeros(UInt16,729)
    for index in eachindex(templates)
        template=templates[index]
        code=1;power=1
        for state in template.faces
            code+=Int(state)*power;power*=3
        end
        previous=Int(canonical[code])
        if previous==0 || template.kind<templates[previous].kind
            canonical[code]=UInt16(index)
        end
    end
    return canonical
end

function _extrude_nonew_quad_strip_macro_relation(canonical,side::Int,preferred)
    templates=_extrude_nonew_templates()
    relation=zeros(UInt16,27);costs=fill(UInt8(0xff),27)
    for packed in canonical
        packed==0 && continue
        index=Int(packed);template=templates[index]
        admissible=true;changes=0
        for local_side in 1:4
            local_side==side && continue
            state=template.faces[local_side]
            state!=0 || (admissible=false;break)
            changes+=state!=preferred[local_side]
        end
        admissible || continue
        slot=_extrude_nonew_quad_strip_relation_slot(template.faces[5],template.faces[6],template.faces[side])
        previous=Int(relation[slot])
        if changes<costs[slot] || (changes==costs[slot] &&
                (previous==0 || (template.kind,index)<(templates[previous].kind,previous)))
            relation[slot]=packed;costs[slot]=UInt8(changes)
        end
    end
    return relation,costs
end

function _extrude_nonew_quad_strip_chain(intervals::Int,shared_sides,preferred,top_states,caller)
    intervals>0 || _extrude_nonew_quad_strip_error(caller,"requires positive intervals")
    canonical=_extrude_nonew_quad_strip_mask_factories()
    first,first_costs=_extrude_nonew_quad_strip_macro_relation(canonical,Int(shared_sides[1]),preferred[1])
    second,second_costs=_extrude_nonew_quad_strip_macro_relation(canonical,Int(shared_sides[2]),preferred[2])
    templates=_extrude_nonew_templates()
    empty=_ExtrudeNoNewStripTransition(0,0,0xff,0xff)
    relation=fill(empty,81)
    for lower_b in 0:2,lower_a in 0:2,upper_b in 0:2,upper_a in 0:2
        lower=_extrude_nonew_quad_strip_pair(lower_a,lower_b)
        upper=_extrude_nonew_quad_strip_pair(upper_a,upper_b)
        best_rank=(typemax(Int),typemax(Int),typemax(Int),typemax(Int),typemax(Int),typemax(Int))
        for shared in 0:2
            first_slot=_extrude_nonew_quad_strip_relation_slot(lower_a,upper_a,shared)
            second_slot=_extrude_nonew_quad_strip_relation_slot(lower_b,upper_b,
                _extrude_nonew_quad_strip_reverse(UInt8(shared)))
            a=Int(first[first_slot]);b=Int(second[second_slot])
            a!=0 && b!=0 || continue
            changes=Int(first_costs[first_slot])+Int(second_costs[second_slot])
            rank=(changes,Int(templates[a].kind),a,Int(templates[b].kind),b,shared)
            if rank<best_rank
                best_rank=rank
                alignment=Int(upper_a!=top_states[1])+Int(upper_b!=top_states[2])
                relation[lower+9(upper-1)]=_ExtrudeNoNewStripTransition(
                    UInt16(a),UInt16(b),UInt8(changes),UInt8(alignment))
            end
        end
    end
    # Preference changes dominate tentative top-cap alignment. Stable logical
    # state/table order breaks ties; upstream pointer order is not reproduced.
    infinity=(typemax(Int),typemax(Int))
    previous=fill(infinity,9);next=fill(infinity,9);previous[1]=(0,0)
    parents=zeros(UInt8,9,intervals)
    final=_extrude_nonew_quad_strip_pair(Int(top_states[1]),Int(top_states[2]))
    for interval in 1:intervals
        fill!(next,infinity)
        for upper in 1:9
            interval==intervals && upper!=final && continue
            for lower in 1:9
                previous[lower]==infinity && continue
                transition=relation[lower+9(upper-1)]
                transition.first!=0 || continue
                score=(Base.checked_add(previous[lower][1],Int(transition.changes)),
                       Base.checked_add(previous[lower][2],Int(transition.alignment)))
                if score<next[upper]
                    next[upper]=score;parents[upper,interval]=UInt8(lower)
                end
            end
        end
        previous,next=next,previous
    end
    previous[final]!=infinity || _extrude_nonew_quad_strip_error(caller,"has no conforming joined cap chain")
    picks=Matrix{UInt16}(undef,2,intervals)
    upper=final
    for interval in intervals:-1:1
        lower=Int(parents[upper,interval])
        lower>0 || _extrude_nonew_quad_strip_error(caller,"cap-chain predecessor is missing")
        transition=relation[lower+9(upper-1)]
        picks[1,interval]=transition.first;picks[2,interval]=transition.second
        upper=lower
    end
    upper==1 || _extrude_nonew_quad_strip_error(caller,"cap chain does not start with whole source faces")
    return picks
end

function _extrude_nonew_quad_strip_counts(picks,recombined::Bool,intervals::Int,shared_sides,caller)
    if recombined
        counts=(4,Base.checked_mul(2,intervals-1),0,10)
        return counts,Base.checked_add(Base.checked_mul(9,intervals),30),2
    end
    counts=zeros(Int,4);incidences=0
    diagonals=Base.checked_mul(6,intervals)
    templates=_extrude_nonew_templates()
    for interval in 1:intervals,source_cell in 1:2
        template=templates[Int(picks[source_cell,interval])]
        for position in 1:Int(template.ncells)
            family=Int(template.cells[position].msh)-3
            1<=family<=4 || _extrude_nonew_quad_strip_error(caller,"factory has an unsupported cell family")
            counts[family]=Base.checked_add(counts[family],1)
            incidences=Base.checked_add(incidences,family==1 ? 4 : family==2 ? 6 : 5)
        end
        diagonals=Base.checked_add(diagonals,Int(template.faces[6]!=0))
        source_cell==1 && (diagonals=Base.checked_add(diagonals,Int(template.faces[Int(shared_sides[1])]!=0)))
    end
    boundary=Base.checked_add(Base.checked_mul(12,intervals),6)
    all_faces=Base.checked_add(incidences,boundary)
    iseven(all_faces) || _extrude_nonew_quad_strip_error(caller,"has inconsistent face-count parity")
    return Tuple(counts),all_faces÷2,diagonals
end

function _extrude_nonew_quad_strip_catalog(m::GeoModel,source::Int,source_mesh,
        spec,levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},recombined::Bool,caller)
    source_mesh isa MixedMesh && nnodes(source_mesh)==6 && length(source_mesh.blocks)==1 ||
        _extrude_nonew_quad_strip_error(caller,"requires exactly six nodes and two linear quadrangles")
    block=only(source_mesh.blocks)
    block isa ElementBlock && block.msh==3 && size(block.nodes)==(4,2) ||
        _extrude_nonew_quad_strip_error(caller,"requires exactly six nodes and two linear quadrangles")
    haskey(m.meshing.transfinite_surfaces,source) && _model_surface_recombined(m,source) ||
        _extrude_nonew_quad_strip_error(caller,"requires a recombined native transfinite surface")
    _surface_type(m,source) in (:plane,:ruled) && !haskey(m.surface_geometry,source) ||
        _extrude_nonew_quad_strip_error(caller,"requires native planar filling without auxiliary surface geometry")
    loops=m.surfaces[source]
    length(loops)==1 || _extrude_nonew_quad_strip_error(caller,"requires one four-curve boundary loop")
    signed=m.loops[only(loops)]
    length(signed)==4 && allunique(abs.(signed)) ||
        _extrude_nonew_quad_strip_error(caller,"requires four distinct boundary curves")
    axis,direction=_extrude_nonew_quad_strip_axis(spec,caller)
    coordinates=ntuple(node->ntuple(d->source_mesh.coords[d,node],3),6)
    all(p->all(isfinite,p),coordinates) && allunique(coordinates) ||
        _extrude_nonew_quad_strip_error(caller,"has nonfinite or coincident distinct source nodes")
    height=coordinates[1][axis]
    all(p->p[axis]==height,coordinates) || _extrude_nonew_quad_strip_error(caller,
        "is not in an exact axis-aligned plane normal to the translation")
    curves=ntuple(k->abs(signed[k]),4)
    widths=ntuple(4) do number
        curve=curves[number]
        haskey(m.curves,curve) && _curve_type(m,curve)===:line &&
            _occ_geometry(m,curve)===nothing && !(curve in m.meshing.degenerated) ||
            _extrude_nonew_quad_strip_error(caller,"requires four straight native boundary curves")
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control!==nothing || _extrude_nonew_quad_strip_error(caller,"requires controlled boundary curves")
        width=_flexible_transfinite_nodes(m,control.num_nodes,curve,caller)
        width in (2,3) || _extrude_nonew_quad_strip_error(caller,"requires opposite effective curve counts3 and2")
        UInt8(width)
    end
    widths[1]==widths[3] && widths[2]==widths[4] && widths[1]!=widths[2] ||
        _extrude_nonew_quad_strip_error(caller,"requires opposite effective curve counts3 and2")
    chains=ntuple(4) do number
        curve=curves[number];width=Int(widths[number]);a,b=m.curves[curve]
        a!=b || _extrude_nonew_quad_strip_error(caller,"has a closed boundary curve")
        chain=_extrude_curve_nodes(m,curve,caller)
        length(chain)==width && chain[1]==m.points[a] && chain[width]==m.points[b] ||
            _extrude_nonew_quad_strip_error(caller,"boundary chains must retain their actual nodes")
        ntuple(3) do position
            position>width && return Int32(0)
            node=findfirst(==(chain[position]),coordinates)
            node!==nothing || _extrude_nonew_quad_strip_error(caller,
                "actual boundary chain is absent from the retained source mesh")
            Int32(node)
        end
    end
    boundary=Int32[];sizehint!(boundary,6)
    curve_indices=UInt8[];sizehint!(curve_indices,6)
    starts=ntuple(number->chains[number][signed[number]>0 ? 1 : Int(widths[number])],4)
    for number in 1:4
        width=Int(widths[number]);chain=chains[number]
        final=chain[signed[number]>0 ? width : 1]
        final==starts[mod1(number+1,4)] || _extrude_nonew_quad_strip_error(caller,"has an unclosed CAD boundary cycle")
        for position in 1:width-1
            push!(boundary,chain[signed[number]>0 ? position : width-position+1])
            push!(curve_indices,UInt8(number))
        end
    end
    length(boundary)==6 && allunique(boundary) ||
        _extrude_nonew_quad_strip_error(caller,"boundary must use all six distinct source nodes")
    planar=ntuple(node->_extrude_nonew_quad_patch_plane(coordinates[node],axis),6)
    orientation=orient2(planar[starts[1]],planar[starts[2]],planar[starts[3]])
    orientation!=0 || _extrude_nonew_quad_strip_error(caller,"has a collinear CAD boundary turn")
    for corner in 1:4
        orient2(planar[starts[corner]],planar[starts[mod1(corner+1,4)]],planar[starts[mod1(corner+2,4)]])==orientation ||
            _extrude_nonew_quad_strip_error(caller,"CAD corners must be strictly convex")
    end
    for side in 1:6,other in side+1:6
        first_start=boundary[side];first_end=boundary[mod1(side+1,6)]
        second_start=boundary[other];second_end=boundary[mod1(other+1,6)]
        (first_start==second_start || first_start==second_end ||
         first_end==second_start || first_end==second_end) && continue
        !_extrude_nonew_quad_patch_crossing(planar[first_start],planar[first_end],
            planar[second_start],planar[second_end]) ||
            _extrude_nonew_quad_strip_error(caller,"actual sampled boundary is not simple")
    end
    cells=ntuple(column->ntuple(row->block.nodes[row,column],4),2)
    all(cell->all(id->1<=id<=6,cell) && allunique(cell),cells) ||
        _extrude_nonew_quad_strip_error(caller,"contains repeated or foreign quadrangle nodes")
    exterior=ntuple(side->minmax(boundary[side],boundary[mod1(side+1,6)]),6)
    shared=(Int32(0),Int32(0))
    for cell in cells,side in 1:4
        edge_start,edge_end=_extrude_nonew_quad_patch_edge(cell,side)
        edge=minmax(edge_start,edge_end)
        edge in exterior && continue
        if shared==(0,0)
            shared=edge
        else
            edge==shared || _extrude_nonew_quad_strip_error(caller,"contains more than one interior source edge")
        end
    end
    shared!=(0,0) || _extrude_nonew_quad_strip_error(caller,"has no common source edge")
    edges=(exterior...,shared)
    counts=zeros(UInt8,7);first_direction=fill((Int32(0),Int32(0)),7)
    indices=ntuple(2) do number
        cell=cells[number]
        for corner in 1:4
            orient2(planar[cell[corner]],planar[cell[mod1(corner+1,4)]],planar[cell[mod1(corner+2,4)]])==orientation ||
                _extrude_nonew_quad_strip_error(caller,"actual quadrangles must be strictly convex with CAD-relative winding")
        end
        local_indices=ntuple(4) do side
            a,b=_extrude_nonew_quad_patch_edge(cell,side)
            index=findfirst(==(minmax(a,b)),edges)
            index!==nothing || _extrude_nonew_quad_strip_error(caller,"has a foreign source edge")
            counts[index]+=1
            if index<=6
                (a,b)==(boundary[index],boundary[mod1(index+1,6)]) ||
                    _extrude_nonew_quad_strip_error(caller,"exterior quadrangle edges oppose the CAD boundary")
            elseif counts[index]==1
                first_direction[index]=(a,b)
            else
                counts[index]==2 && first_direction[index]==(b,a) ||
                    _extrude_nonew_quad_strip_error(caller,"common source edge requires opposite incidences")
            end
            UInt8(index)
        end
        count(==(UInt8(7)),local_indices)==1 ||
            _extrude_nonew_quad_strip_error(caller,"each quadrangle requires three exterior and one common side")
        local_indices
    end
    all(==(UInt8(1)),@view counts[1:6]) && counts[7]==2 ||
        _extrude_nonew_quad_strip_error(caller,"source-edge incidence does not cover the complete strip")
    shared_sides=ntuple(number->UInt8(findfirst(==(UInt8(7)),indices[number])),2)
    shared_first,shared_second=_extrude_nonew_quad_patch_edge(cells[1],Int(shared_sides[1]))
    for number in 1:2,node in cells[number]
        (node==shared_first || node==shared_second) && continue
        orient2(planar[shared_first],planar[shared_second],planar[node])==(number==1 ? orientation : -orientation) ||
            _extrude_nonew_quad_strip_error(caller,"actual quadrangles must lie on opposite sides of their common edge")
    end
    preferred=ntuple(2) do number
        ntuple(4) do side
            index=Int(indices[number][side])
            index==7 ? UInt8(0) : signed[Int(curve_indices[index])]>0 ? UInt8(1) : UInt8(2)
        end
    end
    top_states=ntuple(number->_extrude_nonew_top_diagonal(cells[number]) ? UInt8(1) : UInt8(2),2)
    intervals=length(refs)
    picks=recombined ? Matrix{UInt16}(undef,2,0) :
        _extrude_nonew_quad_strip_chain(intervals,shared_sides,preferred,top_states,caller)
    cell_counts,face_capacity,diagonals=_extrude_nonew_quad_strip_counts(picks,recombined,intervals,shared_sides,caller)
    return _ExtrudeNoNewQuadStripCatalog(cells,coordinates,Tuple(boundary),curves,chains,widths,
        edges,indices,shared_sides,preferred,top_states,picks,cell_counts,face_capacity,diagonals,
        Int8(axis),Int8(direction),Int8(orientation),recombined,levels,refs)
end

function _extrude_nonew_quad_strip_product_certify(cols,catalog::_ExtrudeNoNewQuadStripCatalog,spec,caller)
    axis,direction=_extrude_nonew_quad_strip_axis(spec,caller)
    axis==catalog.normal_axis && direction==catalog.normal_direction ||
        _extrude_nonew_quad_strip_error(caller,"translation differs from its source catalog")
    size(cols)==(length(catalog.levels),6) && size(cols,1)==length(catalog.layer_refs)+1 &&
        size(cols,1)>=2 || _extrude_nonew_quad_strip_error(caller,"column dimensions differ from its layers")
    previous=catalog.source_coordinates[1][axis]
    for level in axes(cols,1)
        height=cols[level,1][axis]
        isfinite(height) || _extrude_nonew_quad_strip_error(caller,"has a nonfinite normal plane")
        if level==1
            height==previous || _extrude_nonew_quad_strip_error(caller,"first plane differs from its source")
        else
            (direction>0 ? height>previous : height<previous) ||
                _extrude_nonew_quad_strip_error(caller,"actual normal planes must be strictly ordered without collapse")
        end
        for node in 1:6
            point=cols[level,node];original=catalog.source_coordinates[node]
            all(isfinite,point) && point[axis]==height ||
                _extrude_nonew_quad_strip_error(caller,"actual columns do not share one finite normal plane")
            for dimension in 1:3
                dimension==axis && continue
                point[dimension]==original[dimension] ||
                    _extrude_nonew_quad_strip_error(caller,"actual columns change source in-plane coordinates")
            end
            level==1 && point!=original && _extrude_nonew_quad_strip_error(caller,
                "first column row differs from its retained source")
        end
        previous=height
    end
    return nothing
end

@inline function _extrude_nonew_quad_strip_dimensions(cols,catalog,caller)
    nlevels=size(cols,1)
    size(cols,2)==6 && nlevels==length(catalog.levels) &&
        nlevels==length(catalog.layer_refs)+1 && nlevels>=2 ||
        _extrude_nonew_quad_strip_error(caller,"has inconsistent completed columns")
    return nlevels,nlevels-1
end

@inline _extrude_nonew_quad_strip_node(node::Int32,level::Int,nlevels::Int)=
    Int32((Int(node)-1)*nlevels+level)

@inline function _extrude_nonew_quad_strip_faces(catalog,source_cell::Int,interval::Int,intervals::Int)
    if catalog.recombined
        return (0x00,0x00,0x00,0x00,0x00,interval==intervals ? catalog.top_states[source_cell] : 0x00)
    end
    return _extrude_nonew_templates()[Int(catalog.template_indices[source_cell,interval])].faces
end

# Seven cells on nine logical vertices are deliberately separate from the
# unchanged six-slot existing-corner table. Bases are positively oriented on
# the reference Hex8; the actual product orientation is certified at runtime.
@inline function _extrude_nonew_quad_strip_problem_fan(top_state::UInt8,caller)
    top_state in (1,2) || _extrude_nonew_quad_strip_error(caller,"problem top requires one diagonal")
    first=top_state==1 ? _extrude_nonew_cell(4,6,5,7,9) : _extrude_nonew_cell(4,6,5,8,9)
    second=top_state==1 ? _extrude_nonew_cell(4,7,5,8,9) : _extrude_nonew_cell(4,7,6,8,9)
    return _ExtrudeNoNewStripProblemFan(0x07,(
        _extrude_nonew_cell(7,2,1,5,6,9),_extrude_nonew_cell(7,3,2,6,7,9),
        _extrude_nonew_cell(7,4,3,7,8,9),_extrude_nonew_cell(7,1,4,8,5,9),
        _extrude_nonew_cell(7,1,2,3,4,9),first,second))
end

function _extrude_nonew_quad_strip_emit!(arrays,filled,template,original,
                                       interval::Int,nlevels::Int,center::Int32,caller)
    for position in 1:Int(template.ncells)
        cell=template.cells[position];family=Int(cell.msh)-3
        1<=family<=4 || _extrude_nonew_quad_strip_error(caller,"factory emitted an unsupported cell family")
        nodes=arrays[family];filled[family]+=1;column=filled[family]
        column<=size(nodes,2) || _extrude_nonew_quad_strip_error(caller,"factory exceeded its reserved cell capacity")
        for row in axes(nodes,1)
            local_node=Int(cell.nodes[row])
            if local_node==9
                center>0 || _extrude_nonew_quad_strip_error(caller,"factory has an unallocated centroid")
                nodes[row,column]=center
            else
                1<=local_node<=8 || _extrude_nonew_quad_strip_error(caller,"factory has an invalid corner index")
                source_position=local_node<=4 ? local_node : local_node-4
                level=local_node<=4 ? interval : interval+1
                nodes[row,column]=_extrude_nonew_quad_strip_node(original[source_position],level,nlevels)
            end
        end
    end
    return nothing
end

function _extrude_nonew_quad_strip_finish(m::GeoModel,t::Int,params,source::Int,
        source_mesh,catalog::_ExtrudeNoNewQuadStripCatalog,cols,top,laterals,caller)
    nlevels,intervals=_extrude_nonew_quad_strip_dimensions(cols,catalog,caller)
    params.recomb_laterals==catalog.recombined ||
        _extrude_nonew_quad_strip_error(caller,"lateral policy differs from its completed catalog")
    catalog.recombined || size(catalog.template_indices)==(2,intervals) ||
        _extrude_nonew_quad_strip_error(caller,"cap-chain dimensions differ from its layers")
    column_count=Base.checked_mul(6,nlevels)
    node_count=Base.checked_add(column_count,catalog.recombined ? 2 : 0)
    cell_count=foldl(Base.checked_add,catalog.cell_counts;init=0)
    node_count<=_EXTRUDE_NONEW_MAX_NODES && node_count<=typemax(Int32) && cell_count<=typemax(Int32) ||
        _extrude_nonew_quad_strip_error(caller,"indexed output exceeds the node or Int32 cell limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in 1:6,level in 1:nlevels
        index=_extrude_nonew_quad_strip_node(Int32(node),level,nlevels);point=cols[level,node]
        coordinates[1,index]=point[1];coordinates[2,index]=point[2];coordinates[3,index]=point[3]
    end
    arrays=ntuple(family->Matrix{Int32}(undef,(4,8,6,5)[family],catalog.cell_counts[family]),4)
    filled=zeros(Int,4)
    edges=Set{NTuple{2,NTuple{3,Float64}}}();sizehint!(edges,catalog.diagonal_count)
    orientation=Int(catalog.source_orientation)*Int(catalog.normal_direction)
    templates=_extrude_nonew_templates()
    for source_cell in 1:2,interval in 1:intervals
        original=catalog.source_cells[source_cell];corners=_extrude_nonew_corners(cols,original,interval)
        faces=_extrude_nonew_quad_strip_faces(catalog,source_cell,interval,intervals)
        _extrude_nonew_add_diagonals!(edges,corners,faces)
        if catalog.recombined && interval==intervals
            # This path emits a real center. The strict Float64 centroid guard
            # must run on v8, without the unused-logical-centroid fallback.
            _extrude_nonew_certify(corners,edges,nothing,caller)
            center=_extrude_quadtri_centroid(corners)
            fan=_extrude_nonew_quad_strip_problem_fan(catalog.top_states[source_cell],caller)
            vertices=(corners...,center)
            _extrude_nonew_certify_template(vertices,fan,orientation,caller)
            index=Int32(column_count+source_cell)
            coordinates[1,index]=center[1];coordinates[2,index]=center[2];coordinates[3,index]=center[3]
            _extrude_nonew_quad_strip_emit!(arrays,filled,fan,original,interval,nlevels,index,caller)
        else
            template=templates[catalog.recombined ? 1 : Int(catalog.template_indices[source_cell,interval])]
            _extrude_nonew_certify_template(corners,template,orientation,caller)
            _extrude_nonew_quad_strip_emit!(arrays,filled,template,original,interval,nlevels,Int32(0),caller)
        end
    end
    Tuple(filled)==catalog.cell_counts && length(edges)==catalog.diagonal_count ||
        _extrude_nonew_quad_strip_error(caller,"factory output does not match its exact capacities")
    blocks=ElementBlock[];sizehint!(blocks,4)
    for family in 1:4
        catalog.cell_counts[family]>0 || continue
        _extrude_make_positive!(coordinates,arrays[family],family+3)
        push!(blocks,ElementBlock(family+3,arrays[family]))
    end
    # Keep normal constructor validation and copies. No coordinate-cell arrays
    # or weld dictionary participate in direct source-major indexed emission.
    volume=MixedMesh(coordinates,blocks)
    sweep=_ExtrudeNoNewIndexedSweep(t,cols,volume)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source,
            top_tag=top,lateral_tags=laterals,catalog=catalog)
end

function _extrude_nonew_quad_strip_top(m::GeoModel,t::Int,params,
        cols,catalog::_ExtrudeNoNewQuadStripCatalog,caller)
    nlevels,_=_extrude_nonew_quad_strip_dimensions(cols,catalog,caller)
    coordinates=Matrix{Float64}(undef,3,6)
    for node in 1:6
        point=cols[nlevels,node]
        coordinates[1,node]=point[1];coordinates[2,node]=point[2];coordinates[3,node]=point[3]
    end
    nodes=Matrix{Int32}(undef,3,4)
    for source_cell in 1:2
        a,b,c,d=catalog.source_cells[source_cell];first=2source_cell-1;second=2source_cell
        if catalog.top_states[source_cell]==1
            nodes[1,first]=a;nodes[2,first]=b;nodes[3,first]=c
            nodes[1,second]=a;nodes[2,second]=c;nodes[3,second]=d
        else
            nodes[1,first]=a;nodes[2,first]=b;nodes[3,first]=d
            nodes[1,second]=b;nodes[2,second]=c;nodes[3,second]=d
        end
    end
    return Mesh(coordinates;tris=nodes)
end

function _extrude_nonew_quad_strip_lateral(m::GeoModel,t::Int,params::_GeoExtrudeParams,
        gen_signed::Int,cols,catalog::_ExtrudeNoNewQuadStripCatalog,edges,caller)
    gen=abs(gen_signed);number=findfirst(==(gen),catalog.boundary_curves)
    number!==nothing && haskey(m.curves,gen) || _extrude_nonew_quad_strip_error(caller,
        "lateral Surface[$t] has no certified boundary generatrix")
    width=Int(catalog.boundary_widths[number]);chain=catalog.boundary_chains[number]
    a,b=m.curves[gen]
    catalog.source_coordinates[chain[1]]==m.points[a] &&
        catalog.source_coordinates[chain[width]]==m.points[b] ||
        _extrude_nonew_quad_strip_error(caller,"lateral Curve[$gen] endpoints differ from its retained source")
    nlevels,intervals=_extrude_nonew_quad_strip_dimensions(cols,catalog,caller)
    node_count=Base.checked_mul(width,nlevels)
    node_count<=typemax(Int32) || _extrude_nonew_quad_strip_error(caller,"lateral output exceeds the Int32 node limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for position in 1:width,level in 1:nlevels
        point=cols[level,chain[position]];node=(position-1)*nlevels+level
        coordinates[1,node]=point[1];coordinates[2,node]=point[2];coordinates[3,node]=point[3]
    end
    recombine=_extrude_quadtri_lateral_recombine(m,t,params,caller)
    recombine==catalog.recombined || _extrude_nonew_quad_strip_error(caller,
        "lateral Surface[$t] policy differs from its retained catalog")
    nodes=Matrix{Int32}(undef,recombine ? 4 : 3,
        Base.checked_mul(Base.checked_mul(width-1,intervals),recombine ? 1 : 2))
    for segment in 1:width-1,level in 1:intervals
        v0=cols[level,chain[segment]];v1=cols[level,chain[segment+1]]
        v2=cols[level+1,chain[segment]];v3=cols[level+1,chain[segment+1]]
        i0=Int32((segment-1)*nlevels+level);i1=Int32(segment*nlevels+level)
        i2=Int32((segment-1)*nlevels+level+1);i3=Int32(segment*nlevels+level+1)
        column=(segment-1)*intervals+level
        ac=_extrude_ein(edges,v0,v3);bd=_extrude_ein(edges,v1,v2)
        if recombine
            !ac && !bd || _extrude_nonew_quad_strip_error(caller,
                "recombined lateral Surface[$t] carries a subdivision diagonal")
            nodes[1,column]=i0;nodes[2,column]=i1;nodes[3,column]=i3;nodes[4,column]=i2
        else
            xor(ac,bd) || _extrude_nonew_quad_strip_error(caller,
                "free lateral Surface[$t] requires exactly one certified subdivision diagonal")
            first=2column-1;second=2column
            if bd
                nodes[1,first]=i2;nodes[2,first]=i1;nodes[3,first]=i0
                nodes[1,second]=i2;nodes[2,second]=i3;nodes[3,second]=i1
            else
                nodes[1,first]=i2;nodes[2,first]=i3;nodes[3,first]=i0
                nodes[1,second]=i0;nodes[2,second]=i3;nodes[3,second]=i1
            end
        end
    end
    return recombine ? MixedMesh(coordinates,[ElementBlock(3,nodes)]) : Mesh(coordinates;tris=nodes)
end
