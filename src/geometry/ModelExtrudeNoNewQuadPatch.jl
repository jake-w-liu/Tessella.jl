# One native TF3 quadrangle patch: four actual Quad4 cells share one interior
# source pivot. All identities stay in the retained nine-node source array.
struct _ExtrudeNoNewQuadPatchCatalog
    source_cells::NTuple{4,NTuple{4,Int32}}
    source_coordinates::NTuple{9,NTuple{3,Float64}}
    source_pivot::Int32
    boundary_cycle::NTuple{8,Int32}
    boundary_curves::NTuple{4,Int}
    boundary_chains::NTuple{4,NTuple{3,Int32}}
    edges::NTuple{12,NTuple{2,Int32}}
    edge_indices::NTuple{4,NTuple{4,UInt8}}
    preferred::NTuple{4,NTuple{4,UInt8}}
    # First/later free slabs, then the recombined terminal slab. The unchanged
    # whole-Hex template is fixed. Phase choices are independent of layer count.
    phase_indices::NTuple{3,NTuple{4,UInt16}}
    normal_axis::Int8
    normal_direction::Int8
    source_orientation::Int8
    recombined::Bool
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@noinline function _extrude_nonew_quad_patch_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts 2-by-2 quadrangle source $reason"))
end

function _extrude_nonew_quad_patch_axis(spec,caller)
    spec.type===:translate || _extrude_nonew_quad_patch_error(caller,
        "requires axis-normal translation; other transforms require the source-grid planner")
    all(isfinite,spec.T) && count(!iszero,spec.T)==1 ||
        _extrude_nonew_quad_patch_error(caller,
            "requires one finite nonzero translation component")
    axis=spec.T[1]!=0 ? 1 : spec.T[2]!=0 ? 2 : 3
    return axis,spec.T[axis]>0 ? 1 : -1
end

@inline _extrude_nonew_quad_patch_plane(point,axis::Int)=
    axis==1 ? (point[2],point[3]) : axis==2 ? (point[3],point[1]) : (point[1],point[2])

@inline _extrude_nonew_quad_patch_edge(cell,side::Int)=
    (cell[side],cell[side==4 ? 1 : side+1])

@inline function _extrude_nonew_quad_patch_on_segment(a,b,p)
    return orient2(a,b,p)==0 && min(a[1],b[1])<=p[1]<=max(a[1],b[1]) &&
        min(a[2],b[2])<=p[2]<=max(a[2],b[2])
end

function _extrude_nonew_quad_patch_crossing(a,b,c,d)
    ac=orient2(a,b,c);ad=orient2(a,b,d)
    ca=orient2(c,d,a);cb=orient2(c,d,b)
    return (ac*ad<0 && ca*cb<0) ||
        _extrude_nonew_quad_patch_on_segment(a,b,c) ||
        _extrude_nonew_quad_patch_on_segment(a,b,d) ||
        _extrude_nonew_quad_patch_on_segment(c,d,a) ||
        _extrude_nonew_quad_patch_on_segment(c,d,b)
end

@inline function _extrude_nonew_quad_patch_cap(cell,pivot::Int32)
    return cell[1]==pivot || cell[3]==pivot ? UInt8(1) : UInt8(2)
end

# The center pivot fixes physical shared/cap diagonals, not a mandatory fan
# apex. Some deduplicated later-slab templates use an opposite lower corner.
# Preserve createFullHexElems' corner-fan priority before two-prism slicing.
function _extrude_nonew_quad_patch_pick(cell,indices,preferred,pivot,
                                       phase::Int,caller)
    templates=_extrude_nonew_templates()
    cap=_extrude_nonew_quad_patch_cap(cell,pivot)
    best=0;best_kind=typemax(UInt8);best_cost=typemax(Int)
    for index in eachindex(templates)
        template=templates[index]
        template.faces[5]==(phase==2 ? cap : 0) && template.faces[6]==cap || continue
        cost=0;admissible=true
        for side in 1:4
            state=template.faces[side]
            if indices[side]>8
                _,second=_extrude_nonew_quad_patch_edge(cell,side)
                state==(second==pivot ? 1 : 2) || (admissible=false;break)
            elseif phase==3
                state==0 || (admissible=false;break)
            else
                state in (1,2) || (admissible=false;break)
                cost+=state!=preferred[side]
            end
        end
        admissible || continue
        if template.kind<best_kind || (template.kind==best_kind && cost<best_cost)
            best=index;best_kind=template.kind;best_cost=cost
        end
    end
    best>0 && best_kind==1 || _extrude_nonew_quad_patch_error(caller,
        "has no corner-fan factory satisfying its shared and cap faces")
    template=templates[best]
    tets=0;pyramids=0
    for index in 1:Int(template.ncells)
        cell_type=template.cells[index].msh
        cell_type==4 ? (tets+=1) : cell_type==7 ? (pyramids+=1) :
            _extrude_nonew_quad_patch_error(caller,"selected an incompatible factory family")
    end
    expected=phase==1 ? (4,1) : phase==2 ? (6,0) : (0,3)
    (tets,pyramids)==expected || _extrude_nonew_quad_patch_error(caller,
        "factory priority does not preserve the patch transition families")
    best<=typemax(UInt16) || throw(ErrorException("QuadTriNoNewVerts template index overflow"))
    return UInt16(best)
end

function _extrude_nonew_quad_patch_catalog(m::GeoModel,source::Int,source_mesh,
        spec,levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},
        recombined::Bool,caller)
    source_mesh isa MixedMesh && nnodes(source_mesh)==9 &&
        length(source_mesh.blocks)==1 || _extrude_nonew_quad_patch_error(caller,
            "requires exactly nine nodes and four linear quadrangles")
    block=only(source_mesh.blocks)
    block isa ElementBlock && block.msh==3 && size(block.nodes)==(4,4) ||
        _extrude_nonew_quad_patch_error(caller,
            "requires exactly nine nodes and four linear quadrangles")
    haskey(m.meshing.transfinite_surfaces,source) && _model_surface_recombined(m,source) ||
        _extrude_nonew_quad_patch_error(caller,"requires a recombined native transfinite surface")
    _surface_type(m,source) in (:plane,:ruled) && !haskey(m.surface_geometry,source) ||
        _extrude_nonew_quad_patch_error(caller,
            "requires native planar filling without auxiliary surface geometry")
    loops=m.surfaces[source]
    length(loops)==1 || _extrude_nonew_quad_patch_error(caller,
        "requires one four-curve boundary loop")
    signed=m.loops[only(loops)]
    length(signed)==4 && allunique(abs.(signed)) || _extrude_nonew_quad_patch_error(caller,
        "requires four distinct boundary curves")
    axis,direction=_extrude_nonew_quad_patch_axis(spec,caller)
    coordinates=ntuple(node->ntuple(d->source_mesh.coords[d,node],3),9)
    all(p->all(isfinite,p),coordinates) || _extrude_nonew_quad_patch_error(caller,
        "has nonfinite source coordinates")
    for first_node in 1:8,second_node in first_node+1:9
        coordinates[first_node]!=coordinates[second_node] ||
            _extrude_nonew_quad_patch_error(caller,"has coincident distinct source nodes")
    end
    height=coordinates[1][axis]
    all(p->p[axis]==height,coordinates) || _extrude_nonew_quad_patch_error(caller,
        "is not in an exact axis-aligned plane normal to the translation")
    curves=ntuple(k->abs(signed[k]),4)
    chains=ntuple(4) do boundary_number
        curve=curves[boundary_number]
        haskey(m.curves,curve) && _curve_type(m,curve)===:line &&
            _occ_geometry(m,curve)===nothing && !(curve in m.meshing.degenerated) ||
            _extrude_nonew_quad_patch_error(caller,"requires four straight native boundary curves")
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control!==nothing && _flexible_transfinite_nodes(m,control.num_nodes,curve,caller)==3 ||
            _extrude_nonew_quad_patch_error(caller,"requires three-node transfinite boundary curves")
        a,b=m.curves[curve]
        a!=b || _extrude_nonew_quad_patch_error(caller,"has a closed boundary curve")
        chain=_extrude_curve_nodes(m,curve,caller)
        length(chain)==3 && chain[1]==m.points[a] && chain[3]==m.points[b] ||
            _extrude_nonew_quad_patch_error(caller,"boundary chains must retain their three actual nodes")
        ids=ntuple(3) do position
            id=findfirst(==(chain[position]),coordinates)
            id!==nothing || _extrude_nonew_quad_patch_error(caller,
                "actual boundary chain is absent from the retained source mesh")
            Int32(id)
        end
        allunique(ids) || _extrude_nonew_quad_patch_error(caller,
            "has repeated boundary-chain nodes")
        ids
    end
    boundary=ntuple(8) do position
        curve_index=(position+1)÷2
        chain_index=isodd(position) ? (signed[curve_index]>0 ? 1 : 3) : 2
        chains[curve_index][chain_index]
    end
    allunique(boundary) || _extrude_nonew_quad_patch_error(caller,
        "boundary must have eight distinct source nodes")
    for curve_index in 1:4
        final=chains[curve_index][signed[curve_index]>0 ? 3 : 1]
        final==boundary[mod1(2curve_index+1,8)] || _extrude_nonew_quad_patch_error(caller,
            "has an unclosed CAD boundary cycle")
    end
    planar=ntuple(node->_extrude_nonew_quad_patch_plane(coordinates[node],axis),9)
    starts=(boundary[1],boundary[3],boundary[5],boundary[7])
    orientation=orient2(planar[starts[1]],planar[starts[2]],planar[starts[3]])
    orientation!=0 || _extrude_nonew_quad_patch_error(caller,"has a collinear CAD boundary turn")
    for corner in 1:4
        orient2(planar[starts[corner]],planar[starts[mod1(corner+1,4)]],
                planar[starts[mod1(corner+2,4)]])==orientation ||
            _extrude_nonew_quad_patch_error(caller,"CAD corners must be strictly convex")
    end
    pivot_position=findfirst(node->!(Int32(node) in boundary),1:9)
    pivot_position!==nothing || _extrude_nonew_quad_patch_error(caller,
        "requires one interior source pivot")
    pivot=Int32(pivot_position)
    # The actual sampled boundary is simple and star-shaped about the pivot.
    # This certifies its eight-triangle fan even when a rounded Line middle
    # node is not exactly collinear with its two native CAD endpoints.
    for side in 1:8
        a=boundary[side];b=boundary[mod1(side+1,8)]
        orient2(planar[a],planar[b],planar[pivot])==orientation ||
            _extrude_nonew_quad_patch_error(caller,"actual boundary is not strictly star-shaped about its pivot")
        for other in side+1:8
            c=boundary[other];d=boundary[mod1(other+1,8)]
            (a==c || a==d || b==c || b==d) && continue
            !_extrude_nonew_quad_patch_crossing(planar[a],planar[b],planar[c],planar[d]) ||
                _extrude_nonew_quad_patch_error(caller,"actual sampled boundary is not simple")
        end
    end
    cells=ntuple(column->ntuple(row->block.nodes[row,column],4),4)
    all(cell->all(id->1<=id<=9,cell) && allunique(cell),cells) ||
        _extrude_nonew_quad_patch_error(caller,"contains repeated or foreign quadrangle nodes")
    edges=NTuple{2,Int32}[]
    sizehint!(edges,12)
    for side in 1:8
        push!(edges,minmax(boundary[side],boundary[mod1(side+1,8)]))
    end
    for chain in chains
        push!(edges,minmax(pivot,chain[2]))
    end
    allunique(edges) || _extrude_nonew_quad_patch_error(caller,
        "does not have eight exterior and four distinct pivot edges")
    counts=zeros(UInt8,12)
    first_directions=fill((Int32(0),Int32(0)),12)
    indices=ntuple(4) do cell_index
        cell=cells[cell_index]
        pivot in cell && count(id->id in starts,cell)==1 ||
            _extrude_nonew_quad_patch_error(caller,
                "each quadrangle must contain the pivot and one CAD corner")
        for corner in 1:4
            orient2(planar[cell[corner]],planar[cell[mod1(corner+1,4)]],
                    planar[cell[mod1(corner+2,4)]])==orientation ||
                _extrude_nonew_quad_patch_error(caller,
                    "actual quadrangles must be strictly convex with CAD-relative winding")
        end
        local_indices=ntuple(4) do side
            a,b=_extrude_nonew_quad_patch_edge(cell,side)
            index=findfirst(==(minmax(a,b)),edges)
            index!==nothing || _extrude_nonew_quad_patch_error(caller,
                "contains an edge outside the boundary/pivot source complex")
            counts[index]+=1
            if index<=8
                (a,b)==(boundary[index],boundary[mod1(index+1,8)]) ||
                    _extrude_nonew_quad_patch_error(caller,"exterior quadrangle edges oppose the CAD boundary")
            elseif counts[index]==1
                first_directions[index]=(a,b)
            else
                counts[index]==2 && first_directions[index]==(b,a) ||
                    _extrude_nonew_quad_patch_error(caller,"shared source edges require opposite incidences")
            end
            UInt8(index)
        end
        count(i->i<=8,local_indices)==2 || _extrude_nonew_quad_patch_error(caller,
            "each quadrangle requires two exterior and two shared pivot sides")
        local_indices
    end
    all(==(0x01),@view counts[1:8]) && all(==(0x02),@view counts[9:12]) ||
        _extrude_nonew_quad_patch_error(caller,"source-edge incidence does not cover the complete patch")
    preferred=ntuple(4) do cell_index
        ntuple(4) do side
            index=Int(indices[cell_index][side])
            if index>8
                UInt8(0)
            else
                curve_index=(index+1)÷2
                signed[curve_index]>0 ? UInt8(1) : UInt8(2)
            end
        end
    end
    phases=ntuple(3) do phase
        ntuple(cell_index->_extrude_nonew_quad_patch_pick(cells[cell_index],
            indices[cell_index],preferred[cell_index],pivot,phase,caller),4)
    end
    return _ExtrudeNoNewQuadPatchCatalog(cells,coordinates,pivot,boundary,curves,
        chains,Tuple(edges),indices,preferred,phases,Int8(axis),Int8(direction),
        Int8(orientation),recombined,levels,refs)
end

function _extrude_nonew_quad_patch_product_certify(cols,catalog::_ExtrudeNoNewQuadPatchCatalog,
                                                 spec,caller)
    axis,direction=_extrude_nonew_quad_patch_axis(spec,caller)
    axis==catalog.normal_axis && direction==catalog.normal_direction ||
        _extrude_nonew_quad_patch_error(caller,"translation differs from its source catalog")
    size(cols)==(length(catalog.levels),9) && size(cols,1)==length(catalog.layer_refs)+1 &&
        size(cols,1)>=2 || _extrude_nonew_quad_patch_error(caller,"column dimensions differ from its layers")
    previous=catalog.source_coordinates[1][axis]
    for level in axes(cols,1)
        height=cols[level,1][axis]
        isfinite(height) || _extrude_nonew_quad_patch_error(caller,"has a nonfinite normal plane")
        if level==1
            height==previous || _extrude_nonew_quad_patch_error(caller,"first plane differs from its source")
        else
            (direction>0 ? height>previous : height<previous) ||
                _extrude_nonew_quad_patch_error(caller,"actual normal planes must be strictly ordered without collapse")
        end
        for node in 1:9
            point=cols[level,node];original=catalog.source_coordinates[node]
            all(isfinite,point) && point[axis]==height ||
                _extrude_nonew_quad_patch_error(caller,"actual columns do not share one finite normal plane")
            for dimension in 1:3
                dimension==axis && continue
                point[dimension]==original[dimension] ||
                    _extrude_nonew_quad_patch_error(caller,"actual columns change source in-plane coordinates")
            end
            level==1 && point!=original && _extrude_nonew_quad_patch_error(caller,
                "first column row differs from its retained source")
        end
        previous=height
    end
    return nothing
end

@inline function _extrude_nonew_quad_patch_face_capacity(catalog::_ExtrudeNoNewQuadPatchCatalog)
    intervals=length(catalog.layer_refs)
    return catalog.recombined ? Base.checked_add(Base.checked_mul(16,intervals),24) :
                                Base.checked_mul(56,intervals)
end

@inline function _extrude_nonew_quad_patch_dimensions(cols,catalog,caller)
    nlevels=size(cols,1)
    nlevels==length(catalog.levels) && size(cols,2)==9 &&
        nlevels==length(catalog.layer_refs)+1 && nlevels>=2 ||
        _extrude_nonew_quad_patch_error(caller,"has inconsistent completed columns")
    return nlevels,nlevels-1
end

@inline _extrude_nonew_quad_patch_node(node::Int32,level::Int,nlevels::Int)=
    Int32((Int(node)-1)*nlevels+level)

@inline function _extrude_nonew_quad_patch_template(catalog,cell_index::Int,
                                                   interval::Int,intervals::Int)
    index=catalog.recombined ? (interval==intervals ? catalog.phase_indices[3][cell_index] : 1) :
                              catalog.phase_indices[interval==1 ? 1 : 2][cell_index]
    return _extrude_nonew_templates()[Int(index)]
end

function _extrude_nonew_quad_patch_finish(m::GeoModel,t::Int,params,source::Int,
        source_mesh,catalog::_ExtrudeNoNewQuadPatchCatalog,cols,top,laterals,caller)
    nlevels,intervals=_extrude_nonew_quad_patch_dimensions(cols,catalog,caller)
    params.recomb_laterals==catalog.recombined || _extrude_nonew_quad_patch_error(caller,
        "lateral policy differs from its completed catalog")
    node_count=Base.checked_mul(9,nlevels)
    tet_count=catalog.recombined ? 0 : Base.checked_sub(Base.checked_mul(24,intervals),8)
    hex_count=catalog.recombined ? Base.checked_mul(4,intervals-1) : 0
    pyramid_count=catalog.recombined ? 12 : 4
    cell_count=Base.checked_add(Base.checked_add(tet_count,hex_count),pyramid_count)
    node_count<=_EXTRUDE_NONEW_MAX_NODES && node_count<=typemax(Int32) &&
        cell_count<=typemax(Int32) || _extrude_nonew_quad_patch_error(caller,
            "indexed output exceeds the node or Int32 cell limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in 1:9,level in 1:nlevels
        index=_extrude_nonew_quad_patch_node(Int32(node),level,nlevels)
        point=cols[level,node]
        coordinates[1,index]=point[1]
        coordinates[2,index]=point[2]
        coordinates[3,index]=point[3]
    end
    tets=Matrix{Int32}(undef,4,tet_count)
    hexes=Matrix{Int32}(undef,8,hex_count)
    pyramids=Matrix{Int32}(undef,5,pyramid_count)
    filled=zeros(Int,3)
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    diagonal_count=catalog.recombined ? 8 : Base.checked_mul(16,intervals)
    sizehint!(edges,diagonal_count)
    # Exact source partition x actual ordered normal planes supplies the hull
    # proof; no rounded, unused body centroid participates in this route.
    orientation=Int(catalog.source_orientation)*Int(catalog.normal_direction)
    for source_cell in 1:4,interval in 1:intervals
        original=catalog.source_cells[source_cell]
        corners=_extrude_nonew_corners(cols,original,interval)
        template=_extrude_nonew_quad_patch_template(catalog,source_cell,interval,intervals)
        _extrude_nonew_certify_template(corners,template,orientation,caller)
        _extrude_nonew_add_diagonals!(edges,corners,template.faces)
        for local_cell in 1:Int(template.ncells)
            cell=template.cells[local_cell]
            family=cell.msh==4 ? 1 : cell.msh==5 ? 2 : cell.msh==7 ? 3 : 0
            family!=0 || _extrude_nonew_quad_patch_error(caller,"factory emitted an unsupported cell family")
            nodes=family==1 ? tets : family==2 ? hexes : pyramids
            filled[family]+=1
            column=filled[family]
            column<=size(nodes,2) || _extrude_nonew_quad_patch_error(caller,
                "factory exceeded its reserved cell capacity")
            for row in axes(nodes,1)
                local_node=Int(cell.nodes[row])
                1<=local_node<=8 || _extrude_nonew_quad_patch_error(caller,"factory has an invalid corner index")
                source_position=local_node<=4 ? local_node : local_node-4
                node_level=local_node<=4 ? interval : interval+1
                nodes[row,column]=_extrude_nonew_quad_patch_node(original[source_position],node_level,nlevels)
            end
        end
    end
    filled==[tet_count,hex_count,pyramid_count] && length(edges)==diagonal_count ||
        _extrude_nonew_quad_patch_error(caller,"factory output does not match its exact patch capacities")
    blocks=ElementBlock[]
    sizehint!(blocks,3)
    if tet_count>0
        _extrude_make_positive!(coordinates,tets,4)
        push!(blocks,ElementBlock(4,tets))
    end
    if hex_count>0
        _extrude_make_positive!(coordinates,hexes,5)
        push!(blocks,ElementBlock(5,hexes))
    end
    _extrude_make_positive!(coordinates,pyramids,7)
    push!(blocks,ElementBlock(7,pyramids))
    # Public constructors retain their normal structural/geometry validation
    # and copies. The temporary coordinate-cell arrays and weld map are absent.
    volume=MixedMesh(coordinates,blocks)
    sweep=_ExtrudeNoNewIndexedSweep(t,cols,volume)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source,
            top_tag=top,lateral_tags=laterals,catalog=catalog)
end

# The authoritative interior source pivot, rather than minimum native node
# rank, divides all four copied top quadrangles. Keep source node/cell order.
function _extrude_nonew_quad_patch_top(m::GeoModel,t::Int,params,
        cols,catalog::_ExtrudeNoNewQuadPatchCatalog,caller)
    nlevels,_=_extrude_nonew_quad_patch_dimensions(cols,catalog,caller)
    coordinates=Matrix{Float64}(undef,3,9)
    for node in 1:9
        point=cols[nlevels,node]
        coordinates[1,node]=point[1]
        coordinates[2,node]=point[2]
        coordinates[3,node]=point[3]
    end
    nodes=Matrix{Int32}(undef,3,8)
    for cell_index in 1:4
        a,b,c,d=catalog.source_cells[cell_index]
        first=2cell_index-1;second=2cell_index
        if _extrude_nonew_quad_patch_cap(catalog.source_cells[cell_index],catalog.source_pivot)==1
            nodes[1,first]=a;nodes[2,first]=b;nodes[3,first]=c
            nodes[1,second]=a;nodes[2,second]=c;nodes[3,second]=d
        else
            nodes[1,first]=a;nodes[2,first]=b;nodes[3,first]=d
            nodes[1,second]=b;nodes[2,second]=c;nodes[3,second]=d
        end
    end
    return Mesh(coordinates;tris=nodes)
end

function _extrude_nonew_quad_patch_lateral(m::GeoModel,t::Int,params::_GeoExtrudeParams,
        gen_signed::Int,cols,catalog::_ExtrudeNoNewQuadPatchCatalog,edges,caller)
    gen=abs(gen_signed)
    boundary_number=findfirst(==(gen),catalog.boundary_curves)
    boundary_number!==nothing && haskey(m.curves,gen) || _extrude_nonew_quad_patch_error(caller,
        "lateral Surface[$t] has no certified boundary generatrix")
    chain=catalog.boundary_chains[boundary_number]
    a,b=m.curves[gen]
    catalog.source_coordinates[chain[1]]==m.points[a] &&
        catalog.source_coordinates[chain[3]]==m.points[b] || _extrude_nonew_quad_patch_error(caller,
            "lateral Curve[$gen] endpoints differ from its retained source")
    nlevels,intervals=_extrude_nonew_quad_patch_dimensions(cols,catalog,caller)
    node_count=Base.checked_mul(3,nlevels)
    node_count<=typemax(Int32) || _extrude_nonew_quad_patch_error(caller,
        "lateral indexed output exceeds the Int32 node limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for chain_position in 1:3,level in 1:nlevels
        point=cols[level,chain[chain_position]]
        node=(chain_position-1)*nlevels+level
        coordinates[1,node]=point[1]
        coordinates[2,node]=point[2]
        coordinates[3,node]=point[3]
    end
    recombine=_extrude_quadtri_lateral_recombine(m,t,params,caller)
    recombine==catalog.recombined || _extrude_nonew_quad_patch_error(caller,
        "lateral Surface[$t] policy differs from its retained catalog")
    cells_per_interval=recombine ? 2 : 4
    nodes=Matrix{Int32}(undef,recombine ? 4 : 3,Base.checked_mul(cells_per_interval,intervals))
    for segment in 1:2,level in 1:intervals
        first_source=chain[segment];second_source=chain[segment+1]
        v0=cols[level,first_source];v1=cols[level,second_source]
        v2=cols[level+1,first_source];v3=cols[level+1,second_source]
        i0=Int32((segment-1)*nlevels+level);i1=Int32(segment*nlevels+level)
        i2=Int32((segment-1)*nlevels+level+1);i3=Int32(segment*nlevels+level+1)
        column=(segment-1)*intervals+level
        if recombine
            !_extrude_ein(edges,v1,v2) && !_extrude_ein(edges,v0,v3) ||
                _extrude_nonew_quad_patch_error(caller,
                    "recombined lateral Surface[$t] carries a subdivision diagonal")
            nodes[1,column]=i0;nodes[2,column]=i1;nodes[3,column]=i3;nodes[4,column]=i2
        else
            first=2column-1;second=2column
            ac=_extrude_ein(edges,v0,v3);bd=_extrude_ein(edges,v1,v2)
            xor(ac,bd) || _extrude_nonew_quad_patch_error(caller,
                "free lateral Surface[$t] requires exactly one certified subdivision diagonal")
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
