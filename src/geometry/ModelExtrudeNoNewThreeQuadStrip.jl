# Three actual boundary quadrangles form a path with two shared source edges.
# Both physical shared-face states participate in one constant-width cap solve.
struct _ExtrudeNoNewThreeQuadStripCatalog
    source_cells::NTuple{3,NTuple{4,Int32}}
    source_coordinates::NTuple{8,NTuple{3,Float64}}
    boundary_cycle::NTuple{8,Int32}
    boundary_curves::NTuple{4,Int}
    boundary_chains::NTuple{4,NTuple{4,Int32}}
    boundary_widths::NTuple{4,UInt8}
    edges::NTuple{10,NTuple{2,Int32}}
    edge_indices::NTuple{3,NTuple{4,UInt8}}
    shared_cells::NTuple{2,NTuple{2,UInt8}}
    shared_sides::NTuple{2,NTuple{2,UInt8}}
    preferred::NTuple{3,NTuple{4,UInt8}}
    top_states::NTuple{3,UInt8}
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

@noinline function _extrude_nonew_three_quad_strip_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts three-quadrangle strip source $reason"))
end

@inline _extrude_nonew_three_quad_strip_face_capacity(catalog::_ExtrudeNoNewThreeQuadStripCatalog)=
    catalog.face_capacity

function _extrude_nonew_three_quad_strip_candidate(m::GeoModel,source::Int,sides::Int,
        caller::AbstractString)
    sides==4 && haskey(m.meshing.transfinite_surfaces,source) &&
        _model_surface_recombined(m,source) || return false
    curves=_model_projection_surface_curves(m,source)
    twos=0;fours=0
    for curve in curves
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control===nothing && return false
        count=_flexible_transfinite_nodes(m,control.num_nodes,curve,caller)
        if count==2
            twos+=1
        elseif count==4
            fours+=1
        else
            return false
        end
    end
    return twos==2 && fours==2
end

function _extrude_nonew_three_quad_strip_axis(spec,caller)
    spec.type===:translate || _extrude_nonew_three_quad_strip_error(caller,
        "requires axis-normal translation; other transforms require the source-grid planner")
    all(isfinite,spec.T) && count(!iszero,spec.T)==1 ||
        _extrude_nonew_three_quad_strip_error(caller,"requires one finite nonzero translation component")
    axis=spec.T[1]!=0 ? 1 : spec.T[2]!=0 ? 2 : 3
    return axis,spec.T[axis]>0 ? 1 : -1
end

@inline _extrude_nonew_three_quad_strip_state(a,b,c)=1+Int(a)+3Int(b)+9Int(c)
@inline _extrude_nonew_three_quad_strip_decode(state::Int)=
    ((state-1)%3,((state-1)÷3)%3,(state-1)÷9)
@inline _extrude_nonew_three_quad_strip_slot(lower,upper,first,second)=
    1+Int(lower)+3Int(upper)+9Int(first)+27Int(second)

struct _ExtrudeNoNewThreeStripTransition
    picks::NTuple{3,UInt16}
    changes::UInt8
    alignment::UInt8
end

function _extrude_nonew_three_quad_strip_macro_sides(cells,sides,number::Int)
    first=UInt8(0);second=UInt8(0)
    for edge in 1:2,incidence in 1:2
        cells[edge][incidence]==number || continue
        if first==0
            first=sides[edge][incidence]
        else
            second=sides[edge][incidence]
        end
    end
    return first,second
end

function _extrude_nonew_three_quad_strip_macro_relation(canonical,sides,preferred)
    templates=_extrude_nonew_templates()
    relation=zeros(UInt16,81);costs=fill(UInt8(0xff),81)
    first,second=Int.(sides)
    for packed in canonical
        packed==0 && continue
        index=Int(packed);template=templates[index]
        admissible=true;changes=0
        for side in 1:4
            (side==first || side==second) && continue
            state=template.faces[side]
            state!=0 || (admissible=false;break)
            changes+=state!=preferred[side]
        end
        admissible || continue
        slot=_extrude_nonew_three_quad_strip_slot(template.faces[5],template.faces[6],
            template.faces[first],second==0 ? UInt8(0) : template.faces[second])
        previous=Int(relation[slot])
        if changes<costs[slot] || (changes==costs[slot] &&
                (previous==0 || (template.kind,index)<(templates[previous].kind,previous)))
            relation[slot]=packed;costs[slot]=UInt8(changes)
        end
    end
    return relation,costs
end

@inline function _extrude_nonew_three_quad_strip_shared_states(cells,physical,number::Int)
    first=UInt8(0);second=UInt8(0);found=false
    for edge in 1:2
        state=physical[edge]
        if cells[edge][1]==number
            local_state=state
        elseif cells[edge][2]==number
            local_state=_extrude_nonew_quad_strip_reverse(state)
        else
            continue
        end
        if found
            second=local_state
        else
            first=local_state;found=true
        end
    end
    return first,second
end

function _extrude_nonew_three_quad_strip_relation(shared_cells,shared_sides,
        preferred,top_states,caller)
    canonical=_extrude_nonew_quad_strip_mask_factories()
    macros=ntuple(number->_extrude_nonew_three_quad_strip_macro_relation(canonical,
        _extrude_nonew_three_quad_strip_macro_sides(shared_cells,shared_sides,number),preferred[number]),3)
    templates=_extrude_nonew_templates()
    empty=_ExtrudeNoNewThreeStripTransition((0,0,0),0xff,0xff)
    relation=fill(empty,729)
    for upper in 1:27,lower in 1:27
        lower_states=_extrude_nonew_three_quad_strip_decode(lower)
        upper_states=_extrude_nonew_three_quad_strip_decode(upper)
        best_rank=ntuple(_->typemax(Int),9)
        for second in UInt8(0):UInt8(2),first in UInt8(0):UInt8(2)
            physical=(first,second)
            slots=ntuple(3) do number
                shared=_extrude_nonew_three_quad_strip_shared_states(shared_cells,physical,number)
                _extrude_nonew_three_quad_strip_slot(lower_states[number],upper_states[number],shared...)
            end
            picks=ntuple(number->macros[number][1][slots[number]],3)
            all(!iszero,picks) || continue
            changes=sum(number->Int(macros[number][2][slots[number]]),1:3)
            a,b,c=Int.(picks)
            rank=(changes,Int(templates[a].kind),a,Int(templates[b].kind),b,
                  Int(templates[c].kind),c,Int(first),Int(second))
            if rank<best_rank
                best_rank=rank
                alignment=sum(number->Int(upper_states[number]!=top_states[number]),1:3)
                relation[lower+27(upper-1)]=_ExtrudeNoNewThreeStripTransition(picks,
                    UInt8(changes),UInt8(alignment))
            end
        end
    end
    return relation
end

function _extrude_nonew_three_quad_strip_chain(intervals::Int,shared_cells,shared_sides,
        preferred,top_states,caller)
    intervals>0 || _extrude_nonew_three_quad_strip_error(caller,"requires positive intervals")
    relation=_extrude_nonew_three_quad_strip_relation(shared_cells,shared_sides,
        preferred,top_states,caller)
    # Changes to preferred exterior directions dominate tentative cap alignment.
    # Original cell, factory and state order give reproducible logical ties.
    infinity=(typemax(Int),typemax(Int))
    previous=fill(infinity,27);next=fill(infinity,27);previous[1]=(0,0)
    parents=zeros(UInt8,27,intervals)
    final=_extrude_nonew_three_quad_strip_state(top_states...)
    for interval in 1:intervals
        fill!(next,infinity)
        for upper in 1:27
            interval==intervals && upper!=final && continue
            for lower in 1:27
                previous[lower]==infinity && continue
                transition=relation[lower+27(upper-1)]
                transition.picks[1]!=0 || continue
                score=(Base.checked_add(previous[lower][1],Int(transition.changes)),
                       Base.checked_add(previous[lower][2],Int(transition.alignment)))
                if score<next[upper]
                    next[upper]=score;parents[upper,interval]=UInt8(lower)
                end
            end
        end
        previous,next=next,previous
    end
    previous[final]!=infinity || _extrude_nonew_three_quad_strip_error(caller,
        "has no conforming joined cap chain")
    picks=Matrix{UInt16}(undef,3,intervals)
    upper=final
    for interval in intervals:-1:1
        lower=Int(parents[upper,interval])
        lower>0 || _extrude_nonew_three_quad_strip_error(caller,"cap-chain predecessor is missing")
        transition=relation[lower+27(upper-1)]
        for number in 1:3
            picks[number,interval]=transition.picks[number]
        end
        upper=lower
    end
    upper==1 || _extrude_nonew_three_quad_strip_error(caller,
        "cap chain does not start with whole source faces")
    return picks
end

function _extrude_nonew_three_quad_strip_counts(picks,recombined::Bool,intervals::Int,
        shared_cells,shared_sides,caller)
    if recombined
        counts=(6,Base.checked_mul(3,intervals-1),0,15)
        return counts,Base.checked_add(Base.checked_mul(13,intervals),45),3
    end
    counts=zeros(Int,4);incidences=0
    diagonals=Base.checked_mul(8,intervals)
    templates=_extrude_nonew_templates()
    for interval in 1:intervals,number in 1:3
        template=templates[Int(picks[number,interval])]
        for position in 1:Int(template.ncells)
            family=Int(template.cells[position].msh)-3
            1<=family<=4 || _extrude_nonew_three_quad_strip_error(caller,
                "factory has an unsupported cell family")
            counts[family]=Base.checked_add(counts[family],1)
            incidences=Base.checked_add(incidences,family==1 ? 4 : family==2 ? 6 : 5)
        end
        diagonals=Base.checked_add(diagonals,Int(template.faces[6]!=0))
        for edge in 1:2
            shared_cells[edge][1]==number || continue
            diagonals=Base.checked_add(diagonals,Int(template.faces[Int(shared_sides[edge][1])]!=0))
        end
    end
    boundary=Base.checked_add(Base.checked_mul(16,intervals),9)
    all_faces=Base.checked_add(incidences,boundary)
    iseven(all_faces) || _extrude_nonew_three_quad_strip_error(caller,
        "has inconsistent face-count parity")
    return Tuple(counts),all_faces÷2,diagonals
end

function _extrude_nonew_three_quad_strip_catalog(m::GeoModel,source::Int,source_mesh,
        spec,levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},recombined::Bool,caller)
    source_mesh isa MixedMesh && nnodes(source_mesh)==8 && length(source_mesh.blocks)==1 ||
        _extrude_nonew_three_quad_strip_error(caller,"requires exactly eight nodes and three linear quadrangles")
    block=only(source_mesh.blocks)
    block isa ElementBlock && block.msh==3 && size(block.nodes)==(4,3) ||
        _extrude_nonew_three_quad_strip_error(caller,"requires exactly eight nodes and three linear quadrangles")
    haskey(m.meshing.transfinite_surfaces,source) && _model_surface_recombined(m,source) ||
        _extrude_nonew_three_quad_strip_error(caller,"requires a recombined native transfinite surface")
    _surface_type(m,source) in (:plane,:ruled) && !haskey(m.surface_geometry,source) ||
        _extrude_nonew_three_quad_strip_error(caller,"requires native planar filling without auxiliary surface geometry")
    loops=m.surfaces[source]
    length(loops)==1 || _extrude_nonew_three_quad_strip_error(caller,"requires one four-curve boundary loop")
    signed=m.loops[only(loops)]
    length(signed)==4 && allunique(abs.(signed)) ||
        _extrude_nonew_three_quad_strip_error(caller,"requires four distinct boundary curves")
    axis,direction=_extrude_nonew_three_quad_strip_axis(spec,caller)
    coordinates=ntuple(node->ntuple(d->source_mesh.coords[d,node],3),8)
    all(p->all(isfinite,p),coordinates) && allunique(coordinates) ||
        _extrude_nonew_three_quad_strip_error(caller,"has nonfinite or coincident distinct source nodes")
    height=coordinates[1][axis]
    all(p->p[axis]==height,coordinates) || _extrude_nonew_three_quad_strip_error(caller,
        "is not in an exact axis-aligned plane normal to the translation")
    curves=ntuple(k->abs(signed[k]),4)
    widths=ntuple(4) do number
        curve=curves[number]
        haskey(m.curves,curve) && _curve_type(m,curve)===:line &&
            _occ_geometry(m,curve)===nothing && !(curve in m.meshing.degenerated) ||
            _extrude_nonew_three_quad_strip_error(caller,"requires four straight native boundary curves")
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control!==nothing || _extrude_nonew_three_quad_strip_error(caller,"requires controlled boundary curves")
        width=_flexible_transfinite_nodes(m,control.num_nodes,curve,caller)
        width in (2,4) || _extrude_nonew_three_quad_strip_error(caller,"requires opposite effective curve counts4 and2")
        UInt8(width)
    end
    widths[1]==widths[3] && widths[2]==widths[4] && widths[1]!=widths[2] ||
        _extrude_nonew_three_quad_strip_error(caller,"requires opposite effective curve counts4 and2")
    chains=ntuple(4) do number
        curve=curves[number];width=Int(widths[number]);a,b=m.curves[curve]
        a!=b || _extrude_nonew_three_quad_strip_error(caller,"has a closed boundary curve")
        chain=_extrude_curve_nodes(m,curve,caller)
        length(chain)==width && chain[1]==m.points[a] && chain[width]==m.points[b] ||
            _extrude_nonew_three_quad_strip_error(caller,"boundary chains must retain their actual nodes")
        ntuple(4) do position
            position>width && return Int32(0)
            node=findfirst(==(chain[position]),coordinates)
            node!==nothing || _extrude_nonew_three_quad_strip_error(caller,
                "actual boundary chain is absent from the retained source mesh")
            Int32(node)
        end
    end
    boundary=Int32[];sizehint!(boundary,8)
    curve_indices=UInt8[];sizehint!(curve_indices,8)
    starts=ntuple(number->chains[number][signed[number]>0 ? 1 : Int(widths[number])],4)
    for number in 1:4
        width=Int(widths[number]);chain=chains[number]
        final=chain[signed[number]>0 ? width : 1]
        final==starts[mod1(number+1,4)] || _extrude_nonew_three_quad_strip_error(caller,"has an unclosed CAD boundary cycle")
        for position in 1:width-1
            push!(boundary,chain[signed[number]>0 ? position : width-position+1])
            push!(curve_indices,UInt8(number))
        end
    end
    length(boundary)==8 && allunique(boundary) ||
        _extrude_nonew_three_quad_strip_error(caller,"boundary must use all eight distinct source nodes")
    planar=ntuple(node->_extrude_nonew_quad_patch_plane(coordinates[node],axis),8)
    orientation=orient2(planar[starts[1]],planar[starts[2]],planar[starts[3]])
    orientation!=0 || _extrude_nonew_three_quad_strip_error(caller,"has a collinear CAD boundary turn")
    for corner in 1:4
        orient2(planar[starts[corner]],planar[starts[mod1(corner+1,4)]],planar[starts[mod1(corner+2,4)]])==orientation ||
            _extrude_nonew_three_quad_strip_error(caller,"CAD corners must be strictly convex")
    end
    for side in 1:8,other in side+1:8
        first_start=boundary[side];first_end=boundary[mod1(side+1,8)]
        second_start=boundary[other];second_end=boundary[mod1(other+1,8)]
        (first_start==second_start || first_start==second_end ||
         first_end==second_start || first_end==second_end) && continue
        !_extrude_nonew_quad_patch_crossing(planar[first_start],planar[first_end],
            planar[second_start],planar[second_end]) ||
            _extrude_nonew_three_quad_strip_error(caller,"actual sampled boundary is not simple")
    end
    cells=ntuple(column->ntuple(row->block.nodes[row,column],4),3)
    all(cell->all(id->1<=id<=8,cell) && allunique(cell),cells) ||
        _extrude_nonew_three_quad_strip_error(caller,"contains repeated or foreign quadrangle nodes")
    exterior=ntuple(side->minmax(boundary[side],boundary[mod1(side+1,8)]),8)
    interior=NTuple{2,Int32}[];sizehint!(interior,2)
    for cell in cells,side in 1:4
        a,b=_extrude_nonew_quad_patch_edge(cell,side);edge=minmax(a,b)
        (edge in exterior || edge in interior) && continue
        length(interior)<2 || _extrude_nonew_three_quad_strip_error(caller,
            "contains more than two interior source edges")
        push!(interior,edge)
    end
    length(interior)==2 || _extrude_nonew_three_quad_strip_error(caller,
        "requires two common source edges")
    edges=(exterior...,interior[1],interior[2])
    counts=zeros(UInt8,10);first_direction=fill((Int32(0),Int32(0)),10)
    indices=ntuple(3) do number
        cell=cells[number]
        for corner in 1:4
            orient2(planar[cell[corner]],planar[cell[mod1(corner+1,4)]],
                planar[cell[mod1(corner+2,4)]])==orientation ||
                _extrude_nonew_three_quad_strip_error(caller,
                    "actual quadrangles must be strictly convex with CAD-relative winding")
        end
        ntuple(4) do side
            a,b=_extrude_nonew_quad_patch_edge(cell,side)
            index=findfirst(==(minmax(a,b)),edges)
            index!==nothing || _extrude_nonew_three_quad_strip_error(caller,"has a foreign source edge")
            counts[index]+=1
            if index<=8
                (a,b)==(boundary[index],boundary[mod1(index+1,8)]) ||
                    _extrude_nonew_three_quad_strip_error(caller,
                        "exterior quadrangle edges oppose the CAD boundary")
            elseif counts[index]==1
                first_direction[index]=(a,b)
            else
                counts[index]==2 && first_direction[index]==(b,a) ||
                    _extrude_nonew_three_quad_strip_error(caller,
                        "common source edges require opposite incidences")
            end
            UInt8(index)
        end
    end
    all(==(UInt8(1)),@view counts[1:8]) && counts[9]==2 && counts[10]==2 ||
        _extrude_nonew_three_quad_strip_error(caller,
            "source-edge incidence does not cover the complete strip")
    shared_cells=ntuple(2) do edge
        numbers=findall(cell->UInt8(edge+8) in cell,indices)
        length(numbers)==2 || _extrude_nonew_three_quad_strip_error(caller,
            "each common edge requires two distinct quadrangles")
        (UInt8(numbers[1]),UInt8(numbers[2]))
    end
    shared_cells[1]!=shared_cells[2] || _extrude_nonew_three_quad_strip_error(caller,
        "quadrangle adjacency must form a path")
    shared_sides=ntuple(edge->ntuple(incidence->UInt8(findfirst(==(UInt8(edge+8)),
        indices[Int(shared_cells[edge][incidence])])),2),2)
    degrees=ntuple(number->count(index->index>8,indices[number]),3)
    count(==(1),degrees)==2 && count(==(2),degrees)==1 ||
        _extrude_nonew_three_quad_strip_error(caller,"quadrangle adjacency must have degrees1/2/1")
    middle=findfirst(==(2),degrees)
    middle_sides=_extrude_nonew_three_quad_strip_macro_sides(shared_cells,shared_sides,middle)
    abs(Int(middle_sides[1])-Int(middle_sides[2]))==2 ||
        _extrude_nonew_three_quad_strip_error(caller,"middle quadrangle requires opposite common sides")
    for edge in 1:2
        first_cell=Int(shared_cells[edge][1]);second_cell=Int(shared_cells[edge][2])
        a,b=_extrude_nonew_quad_patch_edge(cells[first_cell],Int(shared_sides[edge][1]))
        for number in (first_cell,second_cell),node in cells[number]
            (node==a || node==b) && continue
            orient2(planar[a],planar[b],planar[node])==
                (number==first_cell ? orientation : -orientation) ||
                _extrude_nonew_three_quad_strip_error(caller,
                    "adjacent quadrangles must lie on opposite sides of their common edge")
        end
    end
    # Check the complete actual source graph, including both interior edges.
    # Edge crossings plus convex containment checks exclude end-cell overlap.
    for first in 1:10,second in first+1:10
        a,b=edges[first];c,d=edges[second]
        (a==c || a==d || b==c || b==d) && continue
        !_extrude_nonew_quad_patch_crossing(planar[a],planar[b],planar[c],planar[d]) ||
            _extrude_nonew_three_quad_strip_error(caller,"actual source edges cross or touch")
    end
    for first in 1:3,second in first+1:3
        (UInt8(first),UInt8(second)) in shared_cells && continue
        isempty(intersect(cells[first],cells[second])) ||
            _extrude_nonew_three_quad_strip_error(caller,"end quadrangles share an unexpected vertex")
        for (inside,outside) in ((first,second),(second,first)),node in cells[inside]
            contained=true
            for side in 1:4
                a,b=_extrude_nonew_quad_patch_edge(cells[outside],side)
                if orient2(planar[a],planar[b],planar[node])!=orientation
                    contained=false;break
                end
            end
            !contained || _extrude_nonew_three_quad_strip_error(caller,
                "nonadjacent quadrangles overlap by containment")
        end
    end
    preferred=ntuple(3) do number
        ntuple(4) do side
            index=Int(indices[number][side])
            index>8 ? UInt8(0) : signed[Int(curve_indices[index])]>0 ? UInt8(1) : UInt8(2)
        end
    end
    top_states=ntuple(number->_extrude_nonew_top_diagonal(cells[number]) ? UInt8(1) : UInt8(2),3)
    intervals=length(refs)
    picks=recombined ? Matrix{UInt16}(undef,3,0) :
        _extrude_nonew_three_quad_strip_chain(intervals,shared_cells,shared_sides,preferred,top_states,caller)
    cell_counts,face_capacity,diagonals=_extrude_nonew_three_quad_strip_counts(picks,recombined,
        intervals,shared_cells,shared_sides,caller)
    return _ExtrudeNoNewThreeQuadStripCatalog(cells,coordinates,Tuple(boundary),curves,chains,widths,
        edges,indices,shared_cells,shared_sides,preferred,top_states,picks,cell_counts,face_capacity,diagonals,
        Int8(axis),Int8(direction),Int8(orientation),recombined,levels,refs)
end

function _extrude_nonew_three_quad_strip_product_certify(cols,catalog::_ExtrudeNoNewThreeQuadStripCatalog,spec,caller)
    axis,direction=_extrude_nonew_three_quad_strip_axis(spec,caller)
    axis==catalog.normal_axis && direction==catalog.normal_direction ||
        _extrude_nonew_three_quad_strip_error(caller,"translation differs from its source catalog")
    size(cols)==(length(catalog.levels),8) && size(cols,1)==length(catalog.layer_refs)+1 &&
        size(cols,1)>=2 || _extrude_nonew_three_quad_strip_error(caller,"column dimensions differ from its layers")
    previous=catalog.source_coordinates[1][axis]
    for level in axes(cols,1)
        height=cols[level,1][axis]
        isfinite(height) || _extrude_nonew_three_quad_strip_error(caller,"has a nonfinite normal plane")
        if level==1
            height==previous || _extrude_nonew_three_quad_strip_error(caller,"first plane differs from its source")
        else
            (direction>0 ? height>previous : height<previous) ||
                _extrude_nonew_three_quad_strip_error(caller,"actual normal planes must be strictly ordered without collapse")
        end
        for node in 1:8
            point=cols[level,node];original=catalog.source_coordinates[node]
            all(isfinite,point) && point[axis]==height ||
                _extrude_nonew_three_quad_strip_error(caller,"actual columns do not share one finite normal plane")
            for dimension in 1:3
                dimension==axis && continue
                point[dimension]==original[dimension] ||
                    _extrude_nonew_three_quad_strip_error(caller,"actual columns change source in-plane coordinates")
            end
            level==1 && point!=original && _extrude_nonew_three_quad_strip_error(caller,
                "first column row differs from its retained source")
        end
        previous=height
    end
    return nothing
end

@inline function _extrude_nonew_three_quad_strip_dimensions(cols,catalog,caller)
    nlevels=size(cols,1)
    size(cols,2)==8 && nlevels==length(catalog.levels) &&
        nlevels==length(catalog.layer_refs)+1 && nlevels>=2 ||
        _extrude_nonew_three_quad_strip_error(caller,"has inconsistent completed columns")
    return nlevels,nlevels-1
end

@inline _extrude_nonew_three_quad_strip_node(node::Int32,level::Int,nlevels::Int)=
    Int32((Int(node)-1)*nlevels+level)

@inline function _extrude_nonew_three_quad_strip_faces(catalog,source_cell::Int,interval::Int,intervals::Int)
    if catalog.recombined
        return (0x00,0x00,0x00,0x00,0x00,interval==intervals ? catalog.top_states[source_cell] : 0x00)
    end
    return _extrude_nonew_templates()[Int(catalog.template_indices[source_cell,interval])].faces
end

function _extrude_nonew_three_quad_strip_emit!(arrays,filled,template,original,
                                       interval::Int,nlevels::Int,center::Int32,caller)
    for position in 1:Int(template.ncells)
        cell=template.cells[position];family=Int(cell.msh)-3
        1<=family<=4 || _extrude_nonew_three_quad_strip_error(caller,"factory emitted an unsupported cell family")
        nodes=arrays[family];filled[family]+=1;column=filled[family]
        column<=size(nodes,2) || _extrude_nonew_three_quad_strip_error(caller,"factory exceeded its reserved cell capacity")
        for row in axes(nodes,1)
            local_node=Int(cell.nodes[row])
            if local_node==9
                center>0 || _extrude_nonew_three_quad_strip_error(caller,"factory has an unallocated centroid")
                nodes[row,column]=center
            else
                1<=local_node<=8 || _extrude_nonew_three_quad_strip_error(caller,"factory has an invalid corner index")
                source_position=local_node<=4 ? local_node : local_node-4
                level=local_node<=4 ? interval : interval+1
                nodes[row,column]=_extrude_nonew_three_quad_strip_node(original[source_position],level,nlevels)
            end
        end
    end
    return nothing
end

function _extrude_nonew_three_quad_strip_finish(m::GeoModel,t::Int,params,source::Int,
        source_mesh,catalog::_ExtrudeNoNewThreeQuadStripCatalog,cols,top,laterals,caller)
    nlevels,intervals=_extrude_nonew_three_quad_strip_dimensions(cols,catalog,caller)
    params.recomb_laterals==catalog.recombined ||
        _extrude_nonew_three_quad_strip_error(caller,"lateral policy differs from its completed catalog")
    catalog.recombined || size(catalog.template_indices)==(3,intervals) ||
        _extrude_nonew_three_quad_strip_error(caller,"cap-chain dimensions differ from its layers")
    column_count=Base.checked_mul(8,nlevels)
    node_count=Base.checked_add(column_count,catalog.recombined ? 3 : 0)
    cell_count=foldl(Base.checked_add,catalog.cell_counts;init=0)
    node_count<=_EXTRUDE_NONEW_MAX_NODES && node_count<=typemax(Int32) && cell_count<=typemax(Int32) ||
        _extrude_nonew_three_quad_strip_error(caller,"indexed output exceeds the node or Int32 cell limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in 1:8,level in 1:nlevels
        index=_extrude_nonew_three_quad_strip_node(Int32(node),level,nlevels);point=cols[level,node]
        coordinates[1,index]=point[1];coordinates[2,index]=point[2];coordinates[3,index]=point[3]
    end
    arrays=ntuple(family->Matrix{Int32}(undef,(4,8,6,5)[family],catalog.cell_counts[family]),4)
    filled=zeros(Int,4)
    edges=Set{NTuple{2,NTuple{3,Float64}}}();sizehint!(edges,catalog.diagonal_count)
    orientation=Int(catalog.source_orientation)*Int(catalog.normal_direction)
    templates=_extrude_nonew_templates()
    for source_cell in 1:3,interval in 1:intervals
        original=catalog.source_cells[source_cell];corners=_extrude_nonew_corners(cols,original,interval)
        faces=_extrude_nonew_three_quad_strip_faces(catalog,source_cell,interval,intervals)
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
            _extrude_nonew_three_quad_strip_emit!(arrays,filled,fan,original,interval,nlevels,index,caller)
        else
            template=templates[catalog.recombined ? 1 : Int(catalog.template_indices[source_cell,interval])]
            _extrude_nonew_certify_template(corners,template,orientation,caller)
            _extrude_nonew_three_quad_strip_emit!(arrays,filled,template,original,interval,nlevels,Int32(0),caller)
        end
    end
    Tuple(filled)==catalog.cell_counts && length(edges)==catalog.diagonal_count ||
        _extrude_nonew_three_quad_strip_error(caller,"factory output does not match its exact capacities")
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

function _extrude_nonew_three_quad_strip_top(m::GeoModel,t::Int,params,
        cols,catalog::_ExtrudeNoNewThreeQuadStripCatalog,caller)
    nlevels,_=_extrude_nonew_three_quad_strip_dimensions(cols,catalog,caller)
    coordinates=Matrix{Float64}(undef,3,8)
    for node in 1:8
        point=cols[nlevels,node]
        coordinates[1,node]=point[1];coordinates[2,node]=point[2];coordinates[3,node]=point[3]
    end
    nodes=Matrix{Int32}(undef,3,6)
    for source_cell in 1:3
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

function _extrude_nonew_three_quad_strip_lateral(m::GeoModel,t::Int,params::_GeoExtrudeParams,
        gen_signed::Int,cols,catalog::_ExtrudeNoNewThreeQuadStripCatalog,edges,caller)
    gen=abs(gen_signed);number=findfirst(==(gen),catalog.boundary_curves)
    number!==nothing && haskey(m.curves,gen) || _extrude_nonew_three_quad_strip_error(caller,
        "lateral Surface[$t] has no certified boundary generatrix")
    width=Int(catalog.boundary_widths[number]);chain=catalog.boundary_chains[number]
    a,b=m.curves[gen]
    catalog.source_coordinates[chain[1]]==m.points[a] &&
        catalog.source_coordinates[chain[width]]==m.points[b] ||
        _extrude_nonew_three_quad_strip_error(caller,"lateral Curve[$gen] endpoints differ from its retained source")
    nlevels,intervals=_extrude_nonew_three_quad_strip_dimensions(cols,catalog,caller)
    node_count=Base.checked_mul(width,nlevels)
    node_count<=typemax(Int32) || _extrude_nonew_three_quad_strip_error(caller,"lateral output exceeds the Int32 node limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for position in 1:width,level in 1:nlevels
        point=cols[level,chain[position]];node=(position-1)*nlevels+level
        coordinates[1,node]=point[1];coordinates[2,node]=point[2];coordinates[3,node]=point[3]
    end
    recombine=_extrude_quadtri_lateral_recombine(m,t,params,caller)
    recombine==catalog.recombined || _extrude_nonew_three_quad_strip_error(caller,
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
            !ac && !bd || _extrude_nonew_three_quad_strip_error(caller,
                "recombined lateral Surface[$t] carries a subdivision diagonal")
            nodes[1,column]=i0;nodes[2,column]=i1;nodes[3,column]=i3;nodes[4,column]=i2
        else
            xor(ac,bd) || _extrude_nonew_three_quad_strip_error(caller,
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
