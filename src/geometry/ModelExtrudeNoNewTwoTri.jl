# A bounded two-triangle source complex, with original source mesh identities.
# Its five edges include one internal diagonal; only the four exterior edges
# carry CAD Curve provenance. One joined template pair is shared by all levels.
struct _ExtrudeNoNewTwoTriCatalog
    source_cells::NTuple{2,NTuple{3,Int32}}
    source_coordinates::NTuple{4,NTuple{3,Float64}}
    boundary_cycle::NTuple{4,Int32}
    boundary_curves::NTuple{4,Int}
    edges::NTuple{5,NTuple{2,Int32}}
    edge_indices::NTuple{2,NTuple{3,UInt8}}
    shared_sides::NTuple{2,UInt8}
    preferred::NTuple{2,NTuple{3,UInt8}}
    template_indices::NTuple{2,UInt8}
    normal_axis::Int8
    normal_direction::Int8
    source_orientation::Int8
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

# The exact product proof establishes four distinct source columns at every
# strictly ordered level. Final integer connectivity needs no coordinate weld.
struct _ExtrudeNoNewIndexedSweep{M<:Union{Mesh,MixedMesh}}
    tag::Int
    cols::Matrix{NTuple{3,Float64}}
    volume::M
end

@inline _extrude_volume_part(::GeoModel,sweep::_ExtrudeNoNewIndexedSweep,
                             ::AbstractString)=sweep.volume

@noinline function _extrude_nonew_two_tri_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts two-triangle source $reason"))
end

function _extrude_nonew_two_tri_axis(spec,caller)
    spec.type===:translate || _extrude_nonew_two_tri_error(caller,
        "requires axis-normal translation; other transforms require the source-grid planner")
    all(isfinite,spec.T) || _extrude_nonew_two_tri_error(caller,
        "translation must be finite")
    count(!iszero,spec.T)==1 || _extrude_nonew_two_tri_error(caller,
        "requires exactly one nonzero translation component")
    axis=spec.T[1]!=0 ? 1 : spec.T[2]!=0 ? 2 : 3
    return axis,spec.T[axis]>0 ? 1 : -1
end

@inline _extrude_nonew_two_tri_plane(p,axis::Int)=
    axis==1 ? (p[2],p[3]) : axis==2 ? (p[3],p[1]) : (p[1],p[2])

function _extrude_nonew_two_tri_distinct(points)
    for first in 1:3,second in first+1:4
        points[first]==points[second] && return false
    end
    return true
end

@inline function _extrude_nonew_two_tri_edge(cell,side::Int)
    return (cell[side],cell[side==3 ? 1 : side+1])
end

@inline function _extrude_nonew_two_tri_state(cell,side::Int,state::UInt8)
    a,b=_extrude_nonew_two_tri_edge(cell,side)
    return state==0 || a<b ? state : UInt8(3)-state
end

# A hard shared-face equality is distinct from the four adjustable exterior
# preferences. Rank ties by immutable table order, once for the whole region.
function _extrude_nonew_two_tri_pick(cells,shared,preferred,
                                     recombined::Bool,caller)
    templates=_extrude_nonew_prism_templates()
    best=(UInt8(0),UInt8(0))
    best_cost=typemax(Int)
    for first in eachindex(templates),second in eachindex(templates)
        a,b=templates[first],templates[second]
        if recombined
            all(iszero,a.faces) && all(iszero,b.faces) || continue
        else
            all(!iszero,a.faces) && all(!iszero,b.faces) || continue
        end
        sa=Int(shared[1]);sb=Int(shared[2])
        _extrude_nonew_two_tri_state(cells[1],sa,a.faces[sa])==
            _extrude_nonew_two_tri_state(cells[2],sb,b.faces[sb]) || continue
        cost=0
        for side in 1:3
            side==sa || (cost+=a.faces[side]!=preferred[1][side])
            side==sb || (cost+=b.faces[side]!=preferred[2][side])
        end
        if cost<best_cost
            best=(UInt8(first),UInt8(second))
            best_cost=cost
        end
    end
    best[1]!=0 || _extrude_nonew_two_tri_error(caller,
        "has no conforming joined prism templates")
    return best
end

function _extrude_nonew_two_tri_catalog(m::GeoModel,source::Int,source_mesh,
        spec,levels::Vector{Float64},refs::Vector{NTuple{2,Int32}},
        recombined::Bool,caller)
    source_mesh isa Mesh && nnodes(source_mesh)==4 && ntris(source_mesh)==2 &&
        nsegs(source_mesh)==0 && ntets(source_mesh)==0 ||
        _extrude_nonew_two_tri_error(caller,
            "requires exactly four nodes and two linear triangles")
    haskey(m.meshing.transfinite_surfaces,source) ||
        _extrude_nonew_two_tri_error(caller,"requires a native transfinite surface")
    _model_surface_recombined(m,source) &&
        _extrude_nonew_two_tri_error(caller,
            "must retain its two source triangles")
    _surface_type(m,source) in (:plane,:ruled) &&
        !haskey(m.surface_geometry,source) || _extrude_nonew_two_tri_error(caller,
            "requires native planar filling without auxiliary surface geometry")
    loops=m.surfaces[source]
    length(loops)==1 || _extrude_nonew_two_tri_error(caller,
        "requires one four-curve boundary loop")
    signed=m.loops[only(loops)]
    length(signed)==4 && allunique(abs.(signed)) ||
        _extrude_nonew_two_tri_error(caller,"requires four distinct boundary curves")
    axis,direction=_extrude_nonew_two_tri_axis(spec,caller)
    points=ntuple(node->ntuple(d->source_mesh.coords[d,node],3),4)
    all(p->all(isfinite,p),points) && _extrude_nonew_two_tri_distinct(points) ||
        _extrude_nonew_two_tri_error(caller,
            "has nonfinite or coincident distinct source nodes")
    plane=points[1][axis]
    all(p->p[axis]==plane,points) || _extrude_nonew_two_tri_error(caller,
        "is not in an exact axis-aligned plane normal to the translation")

    # Locate each CAD endpoint in the original source node array. Four distinct
    # coordinates make this constant-sized correspondence one-to-one.
    curves=ntuple(k->abs(signed[k]),4)
    starts=ntuple(4) do k
        curve=curves[k]
        haskey(m.curves,curve) && _curve_type(m,curve)===:line &&
            _occ_geometry(m,curve)===nothing || _extrude_nonew_two_tri_error(caller,
                "requires four straight native boundary curves")
        curve in m.meshing.degenerated && _extrude_nonew_two_tri_error(caller,
            "has a degenerated boundary curve")
        cspec=get(m.meshing.transfinite_curves,curve,nothing)
        cspec!==nothing && _flexible_transfinite_nodes(
            m,cspec.num_nodes,curve,caller)==2 || _extrude_nonew_two_tri_error(caller,
            "requires endpoint-only transfinite boundary curves")
        curve_start,curve_end=m.curves[curve]
        curve_start!=curve_end || _extrude_nonew_two_tri_error(caller,"has a closed boundary curve")
        endpoint=signed[k]>0 ? curve_start : curve_end
        position=findfirst(==(m.points[endpoint]),points)
        position===nothing && _extrude_nonew_two_tri_error(caller,
            "CAD boundary endpoint is absent from the source mesh")
        Int32(position)
    end
    allunique(starts) || _extrude_nonew_two_tri_error(caller,
        "boundary loop does not have four distinct corners")
    preferred_boundary=ntuple(4) do k
        curve=curves[k]
        curve_start,curve_end=m.curves[curve]
        endpoint=signed[k]>0 ? curve_end : curve_start
        points[starts[mod1(k+1,4)]]==m.points[endpoint] ||
            _extrude_nonew_two_tri_error(caller,"has an unclosed CAD boundary cycle")
        chain=_extrude_curve_nodes(m,curve,caller)
        length(chain)==2 && chain[1]==m.points[curve_start] &&
            chain[2]==m.points[curve_end] ||
            _extrude_nonew_two_tri_error(caller,
                "boundary chain is not its two actual CAD endpoints")
        # A local face cycle (a,b,b_top,a_top) prefers diagonal(a,b_top)
        # when the native generatrix runs a->b; reversing it swaps the state.
        signed[k]>0 ? UInt8(1) : UInt8(2)
    end
    orientation=orient2((_extrude_nonew_two_tri_plane(points[starts[k]],axis)
                         for k in 1:3)...)
    orientation!=0 || _extrude_nonew_two_tri_error(caller,
        "has a collinear boundary turn")
    for k in 1:4
        turn_first,turn_second,turn_third=
            starts[k],starts[mod1(k+1,4)],starts[mod1(k+2,4)]
        orient2(_extrude_nonew_two_tri_plane(points[turn_first],axis),
                _extrude_nonew_two_tri_plane(points[turn_second],axis),
                _extrude_nonew_two_tri_plane(points[turn_third],axis))==orientation ||
            _extrude_nonew_two_tri_error(caller,
                "boundary must be a strictly convex four-corner polygon")
    end
    cells=ntuple(cell->ntuple(k->source_mesh.tris[k,cell],3),2)
    all(cell->all(id->1<=id<=4,cell) && allunique(cell),cells) ||
        _extrude_nonew_two_tri_error(caller,"contains repeated or foreign triangle nodes")
    for cell in cells
        orient2(_extrude_nonew_two_tri_plane(points[cell[1]],axis),
                _extrude_nonew_two_tri_plane(points[cell[2]],axis),
                _extrude_nonew_two_tri_plane(points[cell[3]],axis))==orientation ||
            _extrude_nonew_two_tri_error(caller,
                "triangle winding must agree with the signed CAD boundary")
    end
    shared_first=0;shared_second=0;common=0
    for first in 1:3,second in 1:3
        ea=_extrude_nonew_two_tri_edge(cells[1],first)
        eb=_extrude_nonew_two_tri_edge(cells[2],second)
        minmax(ea...)==minmax(eb...) || continue
        common+=1
        ea==(eb[2],eb[1]) || _extrude_nonew_two_tri_error(caller,
            "shared edge incidences must have opposite orientation")
        shared_first=first;shared_second=second
    end
    common==1 || _extrude_nonew_two_tri_error(caller,
        "must have exactly one shared source edge")
    shared=(UInt8(shared_first),UInt8(shared_second))
    diagonal=minmax(_extrude_nonew_two_tri_edge(cells[1],shared_first)...)
    boundary=ntuple(k->minmax(starts[k],starts[mod1(k+1,4)]),4)
    diagonal in boundary && _extrude_nonew_two_tri_error(caller,
        "shared edge must join opposite boundary corners")
    edges=(boundary...,diagonal)
    indices=ntuple(2) do cell
        ntuple(3) do side
            edge=_extrude_nonew_two_tri_edge(cells[cell],side)
            index=findfirst(==(minmax(edge...)),edges)
            index===nothing && _extrude_nonew_two_tri_error(caller,
                "triangle edge is absent from the CAD boundary complex")
            if index<=4
                edge==(starts[index],starts[mod1(index+1,4)]) ||
                    _extrude_nonew_two_tri_error(caller,
                        "exterior triangle edge opposes the CAD boundary")
            end
            UInt8(index)
        end
    end
    for edge in 1:4
        count(==(UInt8(edge)),indices[1])+count(==(UInt8(edge)),indices[2])==1 ||
            _extrude_nonew_two_tri_error(caller,
                "CAD boundary edge must have exactly one source incidence")
    end
    diagonal_first,diagonal_second=diagonal
    first_other=only(id for id in cells[1] if id!=diagonal_first && id!=diagonal_second)
    second_other=only(id for id in cells[2] if id!=diagonal_first && id!=diagonal_second)
    first_other!=second_other &&
        orient2(_extrude_nonew_two_tri_plane(points[diagonal_first],axis),
                _extrude_nonew_two_tri_plane(points[diagonal_second],axis),
                _extrude_nonew_two_tri_plane(points[first_other],axis))==
        -orient2(_extrude_nonew_two_tri_plane(points[diagonal_first],axis),
                 _extrude_nonew_two_tri_plane(points[diagonal_second],axis),
                 _extrude_nonew_two_tri_plane(points[second_other],axis)) ||
        _extrude_nonew_two_tri_error(caller,
            "triangles do not form opposite sides of their shared diagonal")
    preferred=ntuple(2) do cell
        ntuple(3) do side
            index=Int(indices[cell][side])
            index==5 ? UInt8(0) : preferred_boundary[index]
        end
    end
    picks=_extrude_nonew_two_tri_pick(cells,shared,preferred,recombined,caller)
    return _ExtrudeNoNewTwoTriCatalog(cells,points,starts,curves,edges,indices,
        shared,preferred,picks,Int8(axis),Int8(direction),Int8(orientation),levels,refs)
end

# Each logical hull is exactly sourceTriangle x ordered normal interval. The
# certified source complex and this actual-column equality jointly prove all
# cross-source and cross-level contacts, without all-pairs hull comparisons.
function _extrude_nonew_two_tri_product_certify(cols,catalog::_ExtrudeNoNewTwoTriCatalog,
                                               spec,caller)
    axis,direction=_extrude_nonew_two_tri_axis(spec,caller)
    axis==catalog.normal_axis && direction==catalog.normal_direction ||
        _extrude_nonew_two_tri_error(caller,"translation differs from its source catalog")
    size(cols)==(length(catalog.levels),4) &&
        size(cols,1)==length(catalog.layer_refs)+1 && size(cols,1)>=2 ||
        _extrude_nonew_two_tri_error(caller,"column dimensions do not match its layer catalog")
    previous=catalog.source_coordinates[1][axis]
    for level in axes(cols,1)
        height=cols[level,1][axis]
        isfinite(height) || _extrude_nonew_two_tri_error(caller,
            "has a nonfinite normal plane")
        if level==1
            height==previous || _extrude_nonew_two_tri_error(caller,
                "first plane does not match the source")
        else
            (direction>0 ? height>previous : height<previous) ||
                _extrude_nonew_two_tri_error(caller,
                    "actual normal planes collapse or reverse after rounding")
        end
        for node in 1:4
            point=cols[level,node]
            all(isfinite,point) && point[axis]==height ||
                _extrude_nonew_two_tri_error(caller,
                    "columns do not share one finite normal plane at each level")
            original=catalog.source_coordinates[node]
            for dimension in 1:3
                dimension==axis && continue
                point[dimension]==original[dimension] ||
                    _extrude_nonew_two_tri_error(caller,
                        "actual columns do not preserve source in-plane coordinates")
            end
            level==1 && point!=original && _extrude_nonew_two_tri_error(caller,
                "first column row does not match the source mesh")
        end
        previous=height
    end
    return nothing
end

@inline _extrude_nonew_two_tri_node(node::Int32,level::Int,nlevels::Int)=
    Int32((Int(node)-1)*nlevels+level)

function _extrude_nonew_two_tri_finish(m::GeoModel,t::Int,params,source::Int,
        source_mesh,catalog::_ExtrudeNoNewTwoTriCatalog,cols,top,laterals,caller)
    nlevels=size(cols,1);intervals=nlevels-1
    nlevels==length(catalog.levels) && size(cols,2)==4 && intervals>0 ||
        _extrude_nonew_two_tri_error(caller,"has inconsistent finalized columns")
    node_count=Base.checked_mul(4,nlevels)
    cell_count=Base.checked_mul(params.recomb_laterals ? 2 : 6,intervals)
    node_count<=_EXTRUDE_NONEW_MAX_NODES && node_count<=typemax(Int32) &&
        cell_count<=typemax(Int32) || _extrude_nonew_two_tri_error(caller,
            "indexed output exceeds the node or Int32 cell limit")
    templates=_extrude_nonew_prism_templates()
    picks=ntuple(k->templates[Int(catalog.template_indices[k])],2)
    for template in picks
        (params.recomb_laterals ? all(iszero,template.faces) :
                                 all(!iszero,template.faces)) ||
            _extrude_nonew_two_tri_error(caller,"lateral policy differs from its joined templates")
    end
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    params.recomb_laterals || sizehint!(edges,Base.checked_mul(5,intervals))
    for level in 1:intervals,cell in 1:2
        corners=_extrude_nonew_corners(cols,catalog.source_cells[cell],level)
        _extrude_nonew_prism_diagonals!(edges,corners,picks[cell].faces)
    end
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in 1:4,level in 1:nlevels
        index=_extrude_nonew_two_tri_node(Int32(node),level,nlevels)
        point=cols[level,node]
        coordinates[1,index]=point[1]
        coordinates[2,index]=point[2]
        coordinates[3,index]=point[3]
    end
    msh=params.recomb_laterals ? 6 : 4
    width=params.recomb_laterals ? 6 : 4
    nodes=Matrix{Int32}(undef,width,cell_count)
    column=0
    # The exact affine source-complex product already proves logical hull
    # convexity. A rounded centroid need not have a representable position
    # strictly between adjacent stored planes; certify the emitted maps using
    # their exact source orientation instead of an unused centroid surrogate.
    orientation=Int(catalog.source_orientation)*Int(catalog.normal_direction)
    for source_cell in 1:2,level in 1:intervals
        original=catalog.source_cells[source_cell]
        template=picks[source_cell]
        corners=_extrude_nonew_corners(cols,original,level)
        _extrude_nonew_certify_template(corners,template,orientation,caller)
        for local_cell in 1:Int(template.ncells)
            cell=template.cells[local_cell]
            Int(cell.msh)==msh || _extrude_nonew_two_tri_error(caller,
                "joined template violates the finalized cell family")
            column+=1
            for row in 1:width
                local_node=Int(cell.nodes[row])
                1<=local_node<=6 || _extrude_nonew_two_tri_error(caller,
                    "joined template has an invalid local node")
                index=local_node<=3 ? local_node : local_node-3
                node_level=local_node<=3 ? level : level+1
                nodes[row,column]=_extrude_nonew_two_tri_node(
                    original[index],node_level,nlevels)
            end
        end
    end
    column==cell_count || _extrude_nonew_two_tri_error(caller,
        "joined templates emitted an unexpected cell count")
    _extrude_make_positive!(coordinates,nodes,msh)
    volume=params.recomb_laterals ? MixedMesh(coordinates,[ElementBlock(6,nodes)]) :
                                   Mesh(coordinates;tets=nodes)
    sweep=_ExtrudeNoNewIndexedSweep(t,cols,volume)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source,
            top_tag=top,lateral_tags=laterals,catalog=catalog)
end

# The completed product has already certified four distinct source columns
# and strictly ordered planes. Each native endpoint-only generatrix selects
# two columns in its own stored direction. Build their source-major surface
# nodes directly, preserving createQuaTri's cell order and selected diagonal.
# Public Mesh/ElementBlock/MixedMesh constructors retain their normal copies.
function _extrude_nonew_two_tri_lateral(m::GeoModel,t::Int,
        params::_GeoExtrudeParams,gen_signed::Int,
        cols::Matrix{NTuple{3,Float64}},catalog::_ExtrudeNoNewTwoTriCatalog,
        edges,caller::AbstractString)
    gen=abs(gen_signed)
    gen in catalog.boundary_curves && haskey(m.curves,gen) ||
        _extrude_nonew_two_tri_error(caller,
            "lateral Surface[$t] has no certified boundary generatrix")
    a,b=m.curves[gen]
    haskey(m.points,a) && haskey(m.points,b) ||
        _extrude_nonew_two_tri_error(caller,
            "lateral Curve[$gen] has an unknown endpoint")
    first=findfirst(==(m.points[a]),catalog.source_coordinates)
    second=findfirst(==(m.points[b]),catalog.source_coordinates)
    first!==nothing && second!==nothing && first!=second ||
        _extrude_nonew_two_tri_error(caller,
            "lateral Curve[$gen] endpoints differ from the certified source")
    nlevels=size(cols,1)
    intervals=nlevels-1
    nlevels==length(catalog.levels) && size(cols,2)==4 && intervals>0 ||
        _extrude_nonew_two_tri_error(caller,
            "lateral Surface[$t] has inconsistent finalized columns")
    node_count=Base.checked_mul(2,nlevels)
    node_count<=typemax(Int32) || _extrude_nonew_two_tri_error(caller,
        "lateral indexed output exceeds the Int32 node limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for (chain,source) in enumerate((first,second)),level in 1:nlevels
        point=cols[level,source]
        node=(chain-1)*nlevels+level
        coordinates[1,node]=point[1]
        coordinates[2,node]=point[2]
        coordinates[3,node]=point[3]
    end
    recombine=_extrude_quadtri_lateral_recombine(m,t,params,caller)
    if recombine
        nodes=Matrix{Int32}(undef,4,intervals)
        for level in 1:intervals
            v0=cols[level,first];v1=cols[level,second]
            v2=cols[level+1,first];v3=cols[level+1,second]
            !_extrude_ein(edges,v1,v2) && !_extrude_ein(edges,v0,v3) ||
                _extrude_nonew_two_tri_error(caller,
                    "recombined lateral Surface[$t] carries a subdivision diagonal")
            nodes[1,level]=Int32(level)
            nodes[2,level]=Int32(nlevels+level)
            nodes[3,level]=Int32(nlevels+level+1)
            nodes[4,level]=Int32(level+1)
        end
        return MixedMesh(coordinates,[ElementBlock(3,nodes)])
    end
    nodes=Matrix{Int32}(undef,3,Base.checked_mul(2,intervals))
    for level in 1:intervals
        v1=cols[level,second];v2=cols[level+1,first]
        i0=Int32(level);i1=Int32(nlevels+level)
        i2=Int32(level+1);i3=Int32(nlevels+level+1)
        if _extrude_ein(edges,v1,v2)
            nodes[1,2level-1]=i2;nodes[2,2level-1]=i1;nodes[3,2level-1]=i0
            nodes[1,2level]=i2;nodes[2,2level]=i3;nodes[3,2level]=i1
        else
            nodes[1,2level-1]=i2;nodes[2,2level-1]=i3;nodes[3,2level-1]=i0
            nodes[1,2level]=i0;nodes[2,2level]=i3;nodes[3,2level]=i1
        end
    end
    return Mesh(coordinates;tris=nodes)
end
