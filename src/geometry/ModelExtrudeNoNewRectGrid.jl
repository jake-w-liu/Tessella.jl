# A dynamic regular source disk is certified once. Physical face propagation
# selects existing-corner factories in source order, without a global cap search.
include("ModelExtrudeNoNewRectGridSource.jl")
abstract type _ExtrudeNoNewDynamicGridCatalog end
include("ModelExtrudeNoNewRectGridPlan.jl")
include("ModelExtrudeNoNewB4StripPlan.jl")

@inline _extrude_nonew_dynamic_grid_face_capacity(catalog::_ExtrudeNoNewDynamicGridCatalog)=
    catalog.face_capacity

function _extrude_nonew_rect_grid_preflight(shape,caller)
    a,b=shape
    try
        vertices=Base.checked_mul(Base.checked_add(a,1),Base.checked_add(b,1))
        cells=Base.checked_mul(a,b)
        return vertices,Base.checked_mul(6,cells)
    catch err
        err isa OverflowError || rethrow()
        _extrude_nonew_rect_grid_error(caller,"source dimensions exceed the indexed size limit")
    end
end

@inline function _extrude_nonew_rect_grid_dimensions(cols,catalog,caller)
    nlevels=size(cols,1)
    nlevels==length(catalog.levels) &&
        size(cols,2)==length(catalog.source.source_coordinates) &&
        nlevels==length(catalog.layer_refs)+1 && nlevels>=2 ||
        _extrude_nonew_rect_grid_error(caller,"has inconsistent completed columns")
    return nlevels,nlevels-1
end

function _extrude_nonew_rect_grid_product_certify(cols,
        catalog::_ExtrudeNoNewDynamicGridCatalog,spec,caller)
    source=catalog.source
    spec.type===:translate && all(isfinite,spec.T) && count(!iszero,spec.T)==1 ||
        _extrude_nonew_rect_grid_error(caller,"requires one finite axis-normal translation")
    axis=spec.T[1]!=0 ? 1 : spec.T[2]!=0 ? 2 : 3
    direction=spec.T[axis]>0 ? 1 : -1
    axis==source.normal_axis && direction==source.normal_direction ||
        _extrude_nonew_rect_grid_error(caller,"translation differs from its source catalog")
    _extrude_nonew_rect_grid_dimensions(cols,catalog,caller)
    previous=source.source_coordinates[1][axis]
    for level in axes(cols,1)
        height=cols[level,1][axis]
        isfinite(height) || _extrude_nonew_rect_grid_error(caller,"has a nonfinite normal plane")
        if level==1
            height==previous || _extrude_nonew_rect_grid_error(caller,"first plane differs from its source")
        else
            (direction>0 ? height>previous : height<previous) ||
                _extrude_nonew_rect_grid_error(caller,"actual normal planes must be strictly ordered without collapse")
        end
        for node in eachindex(source.source_coordinates)
            point=cols[level,node];original=source.source_coordinates[node]
            all(isfinite,point) && point[axis]==height ||
                _extrude_nonew_rect_grid_error(caller,"actual columns do not share one finite normal plane")
            for dimension in 1:3
                dimension==axis && continue
                point[dimension]==original[dimension] ||
                    _extrude_nonew_rect_grid_error(caller,"actual columns change source in-plane coordinates")
            end
            level==1 && point!=original && _extrude_nonew_rect_grid_error(caller,
                "first column row differs from its retained source")
        end
        previous=height
    end
    return nothing
end

@inline _extrude_nonew_rect_grid_node(node::Int32,level::Int,nlevels::Int)=
    Int32((Int(node)-1)*nlevels+level)

function _extrude_nonew_rect_grid_emit!(arrays,filled,template,original,
        interval::Int,nlevels::Int,caller)
    for position in 1:Int(template.ncells)
        cell=template.cells[position];family=Int(cell.msh)-3
        1<=family<=4 || _extrude_nonew_rect_grid_error(caller,
            "factory emitted an unsupported cell family")
        nodes=arrays[family];filled[family]+=1;column=filled[family]
        column<=size(nodes,2) || _extrude_nonew_rect_grid_error(caller,
            "factory exceeded its reserved cell capacity")
        for row in axes(nodes,1)
            local_node=Int(cell.nodes[row])
            1<=local_node<=8 || _extrude_nonew_rect_grid_error(caller,
                "factory has a non-corner vertex")
            source_position=local_node<=4 ? local_node : local_node-4
            level=local_node<=4 ? interval : interval+1
            nodes[row,column]=_extrude_nonew_rect_grid_node(original[source_position],level,nlevels)
        end
    end
    return nothing
end

function _extrude_nonew_rect_grid_finish(m::GeoModel,t::Int,params,source_tag::Int,
        source_mesh,catalog::_ExtrudeNoNewRectGridCatalog,cols,top,laterals,caller)
    nlevels,intervals=_extrude_nonew_rect_grid_dimensions(cols,catalog,caller)
    source=catalog.source;ncells=length(source.source_cells)
    params.recomb_laterals==catalog.recombined ||
        _extrude_nonew_rect_grid_error(caller,"lateral policy differs from its completed catalog")
    size(catalog.template_indices)==(ncells,intervals) ||
        _extrude_nonew_rect_grid_error(caller,"factory dimensions differ from its source and layers")
    node_count=Base.checked_mul(length(source.source_coordinates),nlevels)
    cell_count=foldl(Base.checked_add,catalog.cell_counts;init=0)
    node_count<=_EXTRUDE_NONEW_MAX_NODES && node_count<=typemax(Int32) &&
        cell_count<=typemax(Int32) || _extrude_nonew_rect_grid_error(caller,
            "indexed output exceeds the node or Int32 cell limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in eachindex(source.source_coordinates),level in 1:nlevels
        index=_extrude_nonew_rect_grid_node(Int32(node),level,nlevels);point=cols[level,node]
        coordinates[1,index]=point[1];coordinates[2,index]=point[2];coordinates[3,index]=point[3]
    end
    arrays=ntuple(family->Matrix{Int32}(undef,(4,8,6,5)[family],catalog.cell_counts[family]),4)
    filled=zeros(Int,4)
    edges=Set{NTuple{2,NTuple{3,Float64}}}();sizehint!(edges,catalog.diagonal_count)
    orientation=Int(source.source_orientation)*Int(source.normal_direction)
    templates=_extrude_nonew_templates()
    for source_cell in 1:ncells,interval in 1:intervals
        original=source.source_cells[source_cell]
        corners=_extrude_nonew_corners(cols,original,interval)
        index=Int(catalog.template_indices[source_cell,interval])
        1<=index<=length(templates) || _extrude_nonew_rect_grid_error(caller,
            "has an invalid selected factory")
        template=templates[index]
        _extrude_nonew_certify_template(corners,template,orientation,caller)
        _extrude_nonew_add_diagonals!(edges,corners,template.faces)
        _extrude_nonew_rect_grid_emit!(arrays,filled,template,original,interval,nlevels,caller)
    end
    Tuple(filled)==catalog.cell_counts && length(edges)==catalog.diagonal_count ||
        _extrude_nonew_rect_grid_error(caller,"factory output does not match its exact capacities")
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

function _extrude_nonew_rect_grid_top(m::GeoModel,t::Int,params,
        cols,catalog::_ExtrudeNoNewDynamicGridCatalog,caller)
    nlevels,_=_extrude_nonew_rect_grid_dimensions(cols,catalog,caller)
    source=catalog.source;ncells=length(source.source_cells)
    coordinates=Matrix{Float64}(undef,3,length(source.source_coordinates))
    for node in eachindex(source.source_coordinates)
        point=cols[nlevels,node]
        coordinates[1,node]=point[1];coordinates[2,node]=point[2];coordinates[3,node]=point[3]
    end
    length(catalog.top_states)==ncells ||
        _extrude_nonew_rect_grid_error(caller,"top-state count differs from its source")
    nodes=Matrix{Int32}(undef,3,Base.checked_mul(2,ncells))
    for source_cell in 1:ncells
        a,b,c,d=source.source_cells[source_cell];first=2source_cell-1;second=2source_cell
        state=catalog.top_states[source_cell]
        state in (0x01,0x02) || _extrude_nonew_rect_grid_error(caller,"top cap is not triangulated")
        if state==1
            nodes[1,first]=a;nodes[2,first]=b;nodes[3,first]=c
            nodes[1,second]=a;nodes[2,second]=c;nodes[3,second]=d
        else
            nodes[1,first]=a;nodes[2,first]=b;nodes[3,first]=d
            nodes[1,second]=b;nodes[2,second]=c;nodes[3,second]=d
        end
    end
    return Mesh(coordinates;tris=nodes)
end

function _extrude_nonew_rect_grid_lateral(m::GeoModel,t::Int,params::_GeoExtrudeParams,
        gen_signed::Int,cols,catalog::_ExtrudeNoNewDynamicGridCatalog,edges,caller)
    source=catalog.source;gen=abs(gen_signed)
    number=findfirst(==(gen),source.boundary_curves)
    number!==nothing && haskey(m.curves,gen) || _extrude_nonew_rect_grid_error(caller,
        "lateral Surface[$t] has no certified boundary generatrix")
    width=source.boundary_widths[number];chain=source.boundary_chains[number]
    a,b=m.curves[gen]
    length(chain)==width && source.source_coordinates[chain[1]]==m.points[a] &&
        source.source_coordinates[chain[width]]==m.points[b] ||
        _extrude_nonew_rect_grid_error(caller,"lateral Curve[$gen] differs from its retained chain")
    nlevels,intervals=_extrude_nonew_rect_grid_dimensions(cols,catalog,caller)
    node_count=Base.checked_mul(width,nlevels)
    node_count<=typemax(Int32) || _extrude_nonew_rect_grid_error(caller,
        "lateral output exceeds the Int32 node limit")
    coordinates=Matrix{Float64}(undef,3,node_count)
    for position in 1:width,level in 1:nlevels
        point=cols[level,chain[position]];node=(position-1)*nlevels+level
        coordinates[1,node]=point[1];coordinates[2,node]=point[2];coordinates[3,node]=point[3]
    end
    recombine=_extrude_quadtri_lateral_recombine(m,t,params,caller)
    recombine==catalog.recombined || _extrude_nonew_rect_grid_error(caller,
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
            !ac && !bd || _extrude_nonew_rect_grid_error(caller,
                "recombined lateral Surface[$t] carries a subdivision diagonal")
            nodes[1,column]=i0;nodes[2,column]=i1;nodes[3,column]=i3;nodes[4,column]=i2
        else
            xor(ac,bd) || _extrude_nonew_rect_grid_error(caller,
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
