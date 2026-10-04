# Isolated NoNewVerts region kernels for one nondegenerate source triangle or
# quadrangle, and bounded two-triangle, two-quadrangle and 2-by-2 quadrangle source grids. Quadrangle face choices
# form a cap chain; triangular caps allow independent prism choices per interval.
# All other source categories retain explicit preflight blockers until their
# region-wide propagation phases are implemented.

include("ModelExtrudeQuadTriNoNewTemplates.jl")
include("ModelExtrudeNoNewPrismTemplates.jl")
include("ModelExtrudeNoNewGlobal.jl")
include("ModelExtrudeNoNewJacobian.jl")
include("ModelExtrudeNoNewNonHexJacobian.jl")

const _EXTRUDE_NONEW_MAX_NODES = 10_000_000

struct _ExtrudeNoNewCatalog
    source_cell::Union{NTuple{3,Int32},NTuple{4,Int32}}
    boundary_mask::UInt8
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@inline function _extrude_nonew_top_diagonal(cell)
    first=minimum(cell)
    return first==cell[1] || first==cell[3]
end

function _extrude_nonew_levels(params,caller;source_nodes::Int=4,
        extra_nodes::Int=1,cells_per_interval::Int=6)
    length(params.layers)==length(params.heights) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts layer groups have inconsistent lengths"))
    intervals=0
    previous=0.0
    for group in eachindex(params.layers)
        count=params.layers[group]
        count>0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts requires positive layer counts"))
        height=params.heights[group]
        isfinite(height) && height>previous || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts requires strictly increasing positive layer heights"))
        intervals=try Base.checked_add(intervals,count) catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("$caller: QuadTriNoNewVerts layer count overflows Int"))
        end
        previous=height
    end
    intervals>0 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts requires at least one layer"))
    last(params.heights)==1.0 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts requires a normalized final layer height of 1.0"))
    # This bound precedes level/column allocation. Quadrangles reserve a possible
    # recorded final-cell centroid; triangles never create one. The node limit
    # matches the structured kernels' default.
    intervals<=(_EXTRUDE_NONEW_MAX_NODES-source_nodes-extra_nodes)÷source_nodes || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts exceeds the $_EXTRUDE_NONEW_MAX_NODES node limit"))
    cells_per_interval*intervals<=typemax(Int32) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts output cell count exceeds Int32"))
    levels=_extrude_level_us(params)
    all(isfinite,levels) && all(i->levels[i]>levels[i-1],2:length(levels)) ||
        throw(ArgumentError(
            "$caller: QuadTriNoNewVerts layer levels collapse in Float64"))
    refs=Vector{NTuple{2,Int32}}(undef,intervals)
    position=0
    for group in eachindex(params.layers),local_layer in 1:params.layers[group]
        position+=1
        refs[position]=(Int32(group-1),Int32(local_layer-1))
    end
    return levels,refs
end

function _extrude_nonew_source_cell(source::Mesh,caller)
    ntris(source)==1 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts multiple source cells require the boundary-category planner"))
    cell=ntuple(k->source.tris[k,1],3)
    nnodes(source)==3 && allunique(cell) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts requires one triangle with three distinct source nodes"))
    return cell
end

function _extrude_nonew_source_cell(source::MixedMesh,caller)
    count=0
    cell=(Int32(0),Int32(0),Int32(0),Int32(0))
    for block in source.blocks
        msh_dimension(block.msh)==2 || continue
        block isa ElementBlock && block.msh in (2,3) || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts currently requires one linear source triangle or quadrangle"))
        for column in axes(block.nodes,2)
            count+=1
            count==1 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts multiple source cells require the boundary-category planner"))
            cell=block.msh==2 ? ntuple(k->block.nodes[k,column],3) :
                ntuple(k->block.nodes[k,column],4)
        end
    end
    count==1 && nnodes(source)==length(cell) && allunique(cell) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts requires one cell with distinct source nodes"))
    return cell
end

function _extrude_nonew_region_boundary(m,t,source,caller;sides::Int=4)
    boundary=abs.(_model_volume_boundary_surfaces(m,t,caller))
    length(boundary)==sides+2 && allunique(boundary) && source in boundary ||
        throw(ArgumentError(
            "$caller: QuadTriNoNewVerts requires an isolated $(sides+2)-face sweep"))
    isempty(get(m.embeds,(3,t),NTuple{2,Int}[])) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts embedded volume constraints require the region planner"))
    top=0
    laterals=Int[]
    for surface in boundary
        isempty(get(m.embeds,(2,surface),NTuple{2,Int}[])) || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts embedded surface constraints require the region planner"))
        adjacent=_extrude_quadtri_regions(m,surface,caller)
        length(adjacent)==1 && adjacent[1]==t || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared surfaces require the neighboring-region planner"))
        surface==source && continue
        link=get(m.meshing.extrude_sources,(2,surface),nothing)
        link===nothing && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts boundary surface lacks extrusion provenance"))
        if link[1]==2 && abs(link[2])==source
            top==0 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts has multiple copied caps"))
            top=surface
        elseif link[1]==1
            push!(laterals,surface)
        else
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts has incompatible boundary extrusion provenance"))
        end
    end
    top!=0 && length(laterals)==sides || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts requires one copied cap and $sides laterals"))
    return top,sort!(laterals)
end

function _extrude_nonew_preferences(m,source,cell,source_mesh,caller)
    curves=_model_projection_surface_curves(m,source)
    sides=length(cell)
    length(curves)==sides || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts source needs $sides boundary generatrices"))
    coordinates=ntuple(k->ntuple(d->source_mesh.coords[d,cell[k]],3),sides)
    choices=zeros(UInt8,sides)
    for curve in curves
        chain=_extrude_curve_nodes(m,curve,caller)
        length(chain)==2 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts source boundary contains more than one segment"))
        found=false
        for side in 1:sides
            next=mod1(side+1,sides)
            if chain[1]==coordinates[side] && chain[2]==coordinates[next]
                choices[side]==0 || throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts has duplicate source boundary incidence"))
                choices[side]=1;found=true;break
            elseif chain[2]==coordinates[side] && chain[1]==coordinates[next]
                choices[side]==0 || throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts has duplicate source boundary incidence"))
                choices[side]=2;found=true;break
            end
        end
        found || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts source corners do not match its boundary chains"))
    end
    all(!iszero,choices) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts has a source corner outside its boundary"))
    return ntuple(k->choices[k],sides)
end

# Build the constant-size cap transition relation once for the actual lateral
# preferences. Lateral triangles are required; their orientation is adjustable.
# A tentative cap aligned to the final cap is lower priority than those
# preferences, matching the upstream retry that removes a trapping cap seed.
function _extrude_nonew_chain(intervals::Int,preferred,top_state::UInt8,caller)
    templates=_extrude_nonew_templates()
    length(templates)<=typemax(UInt16) || throw(ErrorException(
        "QuadTriNoNewVerts template table exceeds its packed index range"))
    relation=zeros(UInt16,3,3)
    costs=fill(typemax(Int),3,3)
    for (index,template) in enumerate(templates)
        all(side->template.faces[side]!=0,1:4) || continue
        lower=Int(template.faces[5])+1
        upper=Int(template.faces[6])+1
        cost=4count(side->template.faces[side]!=preferred[side],1:4)
        cost+=template.faces[6]==top_state ? 0 : 1
        if cost<costs[lower,upper]
            costs[lower,upper]=cost
            relation[lower,upper]=UInt16(index)
        end
    end
    parents=zeros(UInt8,3,intervals)
    selected=zeros(UInt16,3,intervals)
    distance=(0,typemax(Int),typemax(Int))
    for layer in 1:intervals
        next_1=typemax(Int);next_2=typemax(Int);next_3=typemax(Int)
        for upper in 1:3
            layer==intervals && upper!=Int(top_state)+1 && continue
            best=typemax(Int)
            for lower in 1:3
                distance[lower]==typemax(Int) && continue
                index=relation[lower,upper]
                index==0 && continue
                score=distance[lower]+costs[lower,upper]
                if score<best
                    best=score
                    parents[upper,layer]=UInt8(lower)
                    selected[upper,layer]=index
                end
            end
            if upper==1
                next_1=best
            elseif upper==2
                next_2=best
            else
                next_3=best
            end
        end
        distance=(next_1,next_2,next_3)
    end
    upper=Int(top_state)+1
    distance[upper]!=typemax(Int) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts has no conforming source-quad cap chain"))
    result=Vector{UInt16}(undef,intervals)
    for layer in intervals:-1:1
        result[layer]=selected[upper,layer]
        upper=Int(parents[upper,layer])
    end
    return result
end

@inline function _extrude_nonew_corners(cols,cell::NTuple{4,Int32},layer)
    return ntuple(k->k<=4 ? cols[layer,cell[k]] : cols[layer+1,cell[k-4]],8)
end

@inline function _extrude_nonew_corners(cols,cell::NTuple{3,Int32},layer)
    return ntuple(k->k<=3 ? cols[layer,cell[k]] : cols[layer+1,cell[k-3]],6)
end

function _extrude_nonew_add_diagonals!(edges,v,states)
    for side in 1:6
        state=states[side]
        state==0 && continue
        face=_EXTRUDE_NONEW_HEX_FACES[side]
        a,b=state==1 ? (face[1],face[3]) : (face[2],face[4])
        push!(edges,_extrude_ekey(v[a],v[b]))
    end
    return nothing
end

function _extrude_nonew_exact_centroid_face(vertices,center,a,b,c,orientation,caller)
    da=vertices[a].-center;db=vertices[b].-center;dc=vertices[c].-center
    determinant=da[1]*(db[2]*dc[3]-db[3]*dc[2])-
                da[2]*(db[1]*dc[3]-db[3]*dc[1])+
                da[3]*(db[1]*dc[2]-db[2]*dc[1])
    determinant*orientation>0 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts exact logical centroid has a folded or degenerate face"))
    return nothing
end

# A centroid used only as a logical boundary witness need not be a representable
# mesh vertex. Adjacent floating-point layer planes can have no Float64 point
# between them. Retry failed witnesses with exact rational arithmetic, while
# emitted template cells retain all their finite positive-volume/map checks.
function _extrude_nonew_exact_centroid_certify(v::NTuple{N,NTuple{3,Float64}},
                                              edges,caller) where N
    family=N==6 ? 6 : 5
    total=0.0
    for (a,b,c,d) in _EXTRUDE_CELL_TETS[family]
        total+=tet_signed_volume(v[a],v[b],v[c],v[d])
    end
    isfinite(total) && total!=0 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts degenerate or nonfinite logical cell volume"))
    orientation=total>0 ? 1 : -1
    vertices=ntuple(i->ntuple(d->Rational{BigInt}(v[i][d]),3),N)
    center=ntuple(d->sum(vertex[d] for vertex in vertices)/N,3)
    faces=N==6 ? _EXTRUDE_QT_PRISM_FACES : _EXTRUDE_QT_HEX_FACES
    for face in faces
        if length(face)==3
            _extrude_nonew_exact_centroid_face(vertices,center,
                face[1],face[2],face[3],orientation,caller)
            continue
        end
        a,b,c,d=face
        first=_extrude_ein(edges,v[a],v[c])
        second=_extrude_ein(edges,v[b],v[d])
        first && second && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts boundary quadrangle has conflicting diagonals"))
        if first || !second
            _extrude_nonew_exact_centroid_face(vertices,center,a,b,c,orientation,caller)
            _extrude_nonew_exact_centroid_face(vertices,center,a,c,d,orientation,caller)
        end
        if second || !first
            _extrude_nonew_exact_centroid_face(vertices,center,a,b,d,orientation,caller)
            _extrude_nonew_exact_centroid_face(vertices,center,b,c,d,orientation,caller)
        end
    end
    return nothing
end

function _extrude_nonew_certify(v,edges,template,caller)
    center=_extrude_quadtri_centroid(v)
    try
        _extrude_quadtri_certify(v,center,edges)
    catch err
        err isa ArgumentError || rethrow()
        if template===nothing
            # This unsliceable branch emits the Float64 centroid as a real
            # vertex, so its original representability checks remain required.
            throw(ArgumentError("$caller: "*replace(err.msg,
                "QuadTriAddVerts"=>"QuadTriNoNewVerts")))
        end
        _extrude_nonew_exact_centroid_certify(v,edges,caller)
    end
    template===nothing && return nothing
    total=0.0
    family=length(v)==6 ? 6 : 5
    for (a,b,c,d) in _EXTRUDE_CELL_TETS[family]
        total+=tet_signed_volume(v[a],v[b],v[c],v[d])
    end
    orientation=total>0 ? 1 : -1
    return _extrude_nonew_certify_template(v,template,orientation,caller)
end

function _extrude_nonew_certify_template(v,template,orientation::Int,caller)
    for position in 1:Int(template.ncells)
        cell=template.cells[position]
        if cell.msh==5
            corners=ntuple(k->v[cell.nodes[k]],8)
            _extrude_nonew_hex_jacobian_certify(corners,orientation,caller)
        end
        for (a,b,c,d) in _EXTRUDE_CELL_TETS[Int(cell.msh)]
            p,q,r,s=(v[cell.nodes[a]],v[cell.nodes[b]],v[cell.nodes[c]],v[cell.nodes[d]])
            -orient3(p,q,r,s)*orientation>0 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts template has a folded or degenerate tetrahedral partition"))
            volume=tet_signed_volume(p,q,r,s)*orientation
            isfinite(volume) && volume>0 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts template volume is not finite and positive"))
        end
        if cell.msh==7
            corners=ntuple(k->v[cell.nodes[k]],5)
            _extrude_nonew_pyramid_jacobian_certify(corners,caller)
        elseif cell.msh==6
            corners=ntuple(k->v[cell.nodes[k]],6)
            _extrude_nonew_prism_jacobian_certify(corners,caller)
        end
    end
    return nothing
end

# Called on an operation-owned working model by the scope helper. No region
# mesh, cache or allocator state is published until its completed face complex
# and every actual cell have passed certification.
function _extrude_nonew_two_tri_candidate(m::GeoModel,source::Int,sides::Int,
                                          caller::AbstractString)
    sides==4 && haskey(m.meshing.transfinite_surfaces,source) || return false
    curves=_model_projection_surface_curves(m,source)
    # Reject known larger grids before source/column allocation. Missing curve
    # controls remain the ordinary native surface mesher's diagnostic.
    for curve in curves
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control===nothing && return false
        _flexible_transfinite_nodes(m,control.num_nodes,curve,caller)==2 ||
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts multiple source cells require the boundary-category planner"))
    end
    return !_model_surface_recombined(m,source)
end

function _extrude_nonew_quad_patch_candidate(m::GeoModel,source::Int,sides::Int,
                                             caller::AbstractString)
    sides==4 && haskey(m.meshing.transfinite_surfaces,source) &&
        _model_surface_recombined(m,source) || return false
    curves=_model_projection_surface_curves(m,source)
    for curve in curves
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control===nothing && return false
        _flexible_transfinite_nodes(m,control.num_nodes,curve,caller)==3 || return false
    end
    return true
end

function _extrude_nonew_quad_strip_candidate(m::GeoModel,source::Int,sides::Int,
                                             caller::AbstractString)
    sides==4 && haskey(m.meshing.transfinite_surfaces,source) &&
        _model_surface_recombined(m,source) || return false
    curves=_model_projection_surface_curves(m,source)
    twos=0
    threes=0
    for curve in curves
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control===nothing && return false
        count=_flexible_transfinite_nodes(m,control.num_nodes,curve,caller)
        if count==2
            twos+=1
        elseif count==3
            threes+=1
        else
            return false
        end
    end
    # The catalog binds the opposite native chains and both actual Quad4 cells.
    return twos==2 && threes==2
end

function _extrude_nonew_plan(m::GeoModel,t::Int,caller::AbstractString;
        min_angle_deg::Real=25.0,max_periodic_passes=8,
        size_field::Union{Nothing,AbstractSizeField}=nothing)
    params=_extrude_gate(m,3,t)
    _extrude_is_nonew(params) || throw(ArgumentError(
        "$caller: Volume[$t] is not a layered QuadTriNoNewVerts region"))
    spec=get(m.meshing.extrude_specs,(3,t),nothing)
    link=get(m.meshing.extrude_sources,(3,t),nothing)
    spec!==nothing && link!==nothing && link[1]==2 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts volume lacks source provenance"))
    spec.type===:rotate && abs(spec.angle)>=2pi && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts full or multiple revolutions require the cyclic/global self-intersection planner"))
    source=abs(link[2])
    haskey(m.meshing.extrude_sources,(2,source)) && throw(ArgumentError(
        "$caller: QuadTriNoNewVerts copied or chained sources require the dependency planner"))
    sides=length(_model_projection_surface_curves(m,source))
    sides in (3,4) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts other source boundaries require the boundary-category planner"))
    top,laterals=_extrude_nonew_region_boundary(m,t,source,caller;sides)
    quad_patch=_extrude_nonew_quad_patch_candidate(m,source,sides,caller)
    quad_strip=!quad_patch && _extrude_nonew_quad_strip_candidate(m,source,sides,caller)
    three_quad_strip=!quad_patch && !quad_strip &&
        _extrude_nonew_three_quad_strip_candidate(m,source,sides,caller)
    four_quad_strip=!quad_patch && !quad_strip && !three_quad_strip &&
        _extrude_nonew_four_quad_strip_candidate(m,source,sides,caller)
    grid_shape=quad_patch || quad_strip || three_quad_strip || four_quad_strip ?
        nothing : _extrude_nonew_rect_grid_shape(m,source,sides,caller)
    rect_grid=grid_shape!==nothing
    two_tri=!quad_patch && !quad_strip && !three_quad_strip && !four_quad_strip && !rect_grid &&
        _extrude_nonew_two_tri_candidate(m,source,sides,caller)
    grid_nodes,grid_cells_per_interval=rect_grid ?
        _extrude_nonew_rect_grid_preflight(grid_shape,caller) : (0,0)
    levels,refs=_extrude_nonew_levels(params,caller;
        source_nodes=rect_grid ? grid_nodes : quad_patch ? 9 : quad_strip ? 6 : three_quad_strip ? 8 : four_quad_strip ? 10 : sides,
        extra_nodes=rect_grid ? 0 : four_quad_strip ? (params.recomb_laterals ? 4 : 0) :
            three_quad_strip ? (params.recomb_laterals ? 3 : 0) :
            quad_strip ? (params.recomb_laterals ? 2 : 0) :
            sides==4 && !two_tri && !quad_patch ? 1 : 0,
        cells_per_interval=rect_grid ? grid_cells_per_interval : four_quad_strip ? (params.recomb_laterals ? 28 : 24) :
            three_quad_strip ? (params.recomb_laterals ? 21 : 18) :
            quad_patch ? 24 : quad_strip ? (params.recomb_laterals ? 14 : 12) :
            two_tri && params.recomb_laterals ? 2 : sides==4 ? 6 : 3)
    source_mesh=mesh_model_surface(m,source;min_angle_deg=min_angle_deg,
        max_periodic_passes=max_periodic_passes,size_field=size_field)
    if quad_patch
        catalog=_extrude_nonew_quad_patch_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_quad_patch_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_quad_patch_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif quad_strip
        catalog=_extrude_nonew_quad_strip_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_quad_strip_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_quad_strip_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif three_quad_strip
        catalog=_extrude_nonew_three_quad_strip_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_three_quad_strip_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_three_quad_strip_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif four_quad_strip
        catalog=_extrude_nonew_four_quad_strip_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_four_quad_strip_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_four_quad_strip_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif rect_grid
        source_data=_extrude_nonew_rect_grid_source(m,source,source_mesh,spec,caller)
        source_data.grid_shape==grid_shape || _extrude_nonew_rect_grid_error(caller,
            "actual source dimensions differ from their preflight")
        catalog=_extrude_nonew_rect_grid_plan(source_data,levels,refs,
            params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_rect_grid_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_rect_grid_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif two_tri
        catalog=_extrude_nonew_two_tri_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_two_tri_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_two_tri_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    end
    cell=_extrude_nonew_source_cell(source_mesh,caller)
    length(cell)==sides || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts source cell does not match its boundary category"))
    preferred=_extrude_nonew_preferences(m,source,cell,source_mesh,caller)
    cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
    all(p->all(isfinite,p),cols) && allunique(cols) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts collapsed or coincident columns require the degenerate-cell planner"))
    _extrude_nonew_global_certify(cols,cell,spec,caller)
    if length(cell)==3
        return _extrude_nonew_triangle_finish(m,t,params,source,source_mesh,
            cell,preferred,cols,levels,refs,top,laterals,caller)
    end
    intervals=length(refs)
    top_state=UInt8(_extrude_nonew_top_diagonal(cell) ? 1 : 2)
    states=Matrix{UInt8}(undef,6,intervals)
    problems=falses(intervals)
    templates=_extrude_nonew_templates()
    picks=zeros(UInt16,intervals)
    if params.recomb_laterals
        fill!(states,0)
        states[6,intervals]=top_state
        # The five required whole faces plus a divided final cap cannot be
        # sliced by either existing-corner construction. Exhaust the relation
        # before recording the upstream all-boundary problem.
        for layer in 1:intervals
            wanted=ntuple(side->states[side,layer],6)
            index=findfirst(template->template.faces==wanted,templates)
            if index===nothing
                layer==intervals || throw(ErrorException(
                    "QuadTriNoNewVerts missing unchanged-hex template"))
                problems[layer]=true
            else
                picks[layer]=UInt16(index)
            end
        end
    else
        picks=_extrude_nonew_chain(intervals,preferred,top_state,caller)
        for layer in 1:intervals,side in 1:6
            states[side,layer]=templates[picks[layer]].faces[side]
        end
    end
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    for layer in 1:intervals
        v=_extrude_nonew_corners(cols,cell,layer)
        _extrude_nonew_add_diagonals!(edges,v,ntuple(side->states[side,layer],6))
    end
    sweep=_ExtrudeVolumeSweep(t,true,cols,
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false])
    for layer in 1:intervals
        v=_extrude_nonew_corners(cols,cell,layer)
        template=problems[layer] ? nothing : templates[picks[layer]]
        _extrude_nonew_certify(v,edges,template,caller)
        if problems[layer]
            first_pyramid=length(sweep.pyramids)+1
            _extrude_quadtri_fan!(sweep,v,edges)
            for index in first_pyramid:length(sweep.pyramids)
                _extrude_nonew_pyramid_jacobian_certify(sweep.pyramids[index],caller)
            end
        else
            _extrude_nonew_emit_template!(sweep,v,template)
        end
    end
    catalog=_ExtrudeNoNewCatalog(cell,0x0f,levels,refs)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source,
        top_tag=top,lateral_tags=laterals,catalog=catalog,faces=states,
        problem_layers=problems)
end

include("ModelExtrudeNoNewTriangle.jl")
include("ModelExtrudeNoNewTwoTri.jl")
include("ModelExtrudeNoNewQuadPatch.jl")
include("ModelExtrudeNoNewQuadStrip.jl")
include("ModelExtrudeNoNewThreeQuadStrip.jl")
include("ModelExtrudeNoNewFourQuadStrip.jl")
include("ModelExtrudeNoNewRectGrid.jl")
