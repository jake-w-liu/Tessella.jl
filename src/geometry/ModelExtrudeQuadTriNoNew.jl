# NoNewVerts kernels for isolated triangles/quadrangles, bounded joined source
# catalogs, complete rectangular Quad disks and recombined all-boundary strips.
# Physical face propagation retains original source-column identities.
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
        extra_nodes::Int=1,extra_nodes_per_interval::Int=0,cells_per_interval::Int=6)
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
    # Reserve all possible recorded centers before allocating levels or columns.
    source_nodes>0 && extra_nodes>=0 && extra_nodes_per_interval>=0 &&
        cells_per_interval>0 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts requires positive source/cell bounds and nonnegative extra-node bounds"))
    node_limit=_EXTRUDE_NONEW_MAX_NODES
    source_nodes<=node_limit && extra_nodes<=node_limit-source_nodes &&
        extra_nodes_per_interval<=node_limit-source_nodes || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts exceeds the $node_limit node limit"))
    stride=source_nodes+extra_nodes_per_interval
    intervals<=(node_limit-source_nodes-extra_nodes)÷stride || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts exceeds the $node_limit node limit"))
    intervals<=typemax(Int32)÷cells_per_interval || throw(ArgumentError(
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

# An axis-normal B4 operation owns two distinct sets of copied CAD Curves
# when automatic coherence is disabled. Their literal tags need not close a
# CAD shell. Admission here binds the recorded operation without welding any
# entities; the completed plan must still certify its actual typed cell shell.
function _extrude_nonew_unmerged_b4_incidence(m::GeoModel,t::Int,caller)
    params=_extrude_gate(m,3,t)
    link=get(m.meshing.extrude_sources,(3,t),nothing)
    spec=get(m.meshing.extrude_specs,(3,t),nothing)
    source=link!==nothing && link[1]==2 ? abs(link[2]) : 0
    source_loops=get(m.surfaces,source,nothing)
    source_loops!==nothing && length(source_loops)==1 || return false
    original=get(m.loops,only(source_loops),nothing)
    original!==nothing && length(original)==4 && allunique(abs.(original)) &&
        all(curve->haskey(m.curves,abs(curve)),original) || return false
    b4=_extrude_is_nonew(params) && spec!==nothing && spec.type===:translate &&
        haskey(m.surfaces,source) &&
        (_extrude_nonew_rect_grid_shape(m,source,4,caller;mode=:b4_strip)!==nothing ||
         _extrude_nonew_rect_grid_shape(m,source,4,caller;mode=:rectangular)!==nothing)
    b4 || return false
    # Unknown entities, repeated members and other incidence defects retain
    # the original CAD diagnostic. Only actual one-owner copied Curve records
    # supply the alternative operation witness.
    length(m.volumes[t])==1 && only(m.volumes[t])>0 || return false
    shell=get(m.surface_loops,only(m.volumes[t]),nothing)
    shell!==nothing && length(shell)==6 && allunique(abs.(shell)) || return false
    source_curves=Set(abs.(original))
    owners=Dict{Int,Set{Int}}();counts=Dict{Int,Int}()
    for signed_surface in shell
        surface=abs(signed_surface)
        haskey(m.surfaces,surface) && length(m.surfaces[surface])==1 || return false
        curves=get(m.loops,only(m.surfaces[surface]),nothing)
        curves!==nothing && length(curves)==4 && allunique(abs.(curves)) || return false
        for signed_curve in curves
            curve=abs(signed_curve)
            haskey(m.curves,curve) || return false
            counts[curve]=get(counts,curve,0)+1
            push!(get!(Set{Int},owners,curve),surface)
        end
    end
    unmerged=false
    for (curve,count) in counts
        count==2 && length(owners[curve])==2 && continue
        curve_link=get(m.meshing.extrude_sources,(1,curve),nothing)
        count==1 && !(curve in source_curves) && curve_link!==nothing &&
            curve_link[1] in (0,1) || return false
        unmerged=true
    end
    return unmerged
end

function _extrude_nonew_boundary_surfaces(m::GeoModel,t::Int,caller)
    try
        return _model_volume_boundary_surfaces(m,t,caller)
    catch error
        error isa ArgumentError || rethrow()
        _extrude_nonew_unmerged_b4_incidence(m,t,caller) || rethrow()
    end
    params=_extrude_gate(m,3,t)
    link=m.meshing.extrude_sources[(3,t)]
    spec=m.meshing.extrude_specs[(3,t)]
    source=abs(link[2])
    _extrude_nonew_rect_grid_axis(spec,caller)
    boundary=_model_volume_boundary_surfaces(m,t,caller;certify_shells=false)
    length(m.volumes[t])==1 && only(m.volumes[t])>0 &&
        length(boundary)==6 && allunique(abs.(boundary)) && -source in boundary ||
        throw(ArgumentError("$caller: B4 operation requires one recorded six-face shell"))
    original=m.loops[only(m.surfaces[source])]
    top=0;laterals=Set{Int}()
    for signed_surface in boundary
        surface=abs(signed_surface)
        surface==source && continue
        isequal(get(m.meshing.extrude,(2,surface),nothing),params) &&
            isequal(get(m.meshing.extrude_specs,(2,surface),nothing),spec) ||
            throw(ArgumentError("$caller: B4 boundary changes its recorded sweep parameters"))
        surface_link=get(m.meshing.extrude_sources,(2,surface),nothing)
        surface_link===nothing && throw(ArgumentError(
            "$caller: B4 boundary lacks recorded Source provenance"))
        length(m.surfaces[surface])==1 || throw(ArgumentError(
            "$caller: B4 boundary has incompatible recorded Curve loops"))
        loop=only(m.surfaces[surface])
        haskey(m.loops,loop) || throw(ArgumentError(
            "$caller: B4 boundary references an unknown Curve loop"))
        curves=m.loops[loop]
        length(curves)==4 && allunique(abs.(curves)) &&
            all(curve->haskey(m.curves,abs(curve)),curves) ||
            throw(ArgumentError("$caller: B4 boundary requires four recorded Curves"))
        if surface_link==link
            top==0 && signed_surface>0 || throw(ArgumentError(
                "$caller: B4 operation requires one positively oriented copied cap"))
            for (curve,parent) in zip(curves,original)
                get(m.meshing.extrude_sources,(1,abs(curve)),nothing)==(1,parent) &&
                    isequal(get(m.meshing.extrude_specs,(1,abs(curve)),nothing),spec) &&
                    isequal(get(m.meshing.extrude,(1,abs(curve)),nothing),params) ||
                    throw(ArgumentError("$caller: B4 copied cap loses its recorded Source Curve binding"))
            end
            top=surface
        elseif surface_link[1]==1
            parent=surface_link[2]
            parent in original && !(parent in laterals) &&
                sign(signed_surface)==sign(parent) && parent in curves ||
                throw(ArgumentError("$caller: B4 laterals do not bijectively bind the Source Curves"))
            copied=0;connecting=Int[]
            for signed_curve in curves
                signed_curve==parent && continue
                curve=abs(signed_curve)
                curve_link=get(m.meshing.extrude_sources,(1,curve),nothing)
                isequal(get(m.meshing.extrude_specs,(1,curve),nothing),spec) &&
                    isequal(get(m.meshing.extrude,(1,curve),nothing),params) ||
                    throw(ArgumentError("$caller: B4 lateral Curve changes its recorded sweep parameters"))
                if curve_link==(1,parent)
                    copied+=1
                elseif curve_link!==nothing && curve_link[1]==0
                    push!(connecting,abs(curve_link[2]))
                else
                    throw(ArgumentError("$caller: B4 lateral Curve has incompatible Source provenance"))
                end
            end
            copied==1 && sort!(connecting)==sort!(collect(m.curves[abs(parent)])) ||
                throw(ArgumentError("$caller: B4 lateral loses its copied Curve or endpoint bindings"))
            push!(laterals,parent)
        else
            throw(ArgumentError("$caller: B4 boundary has incompatible Source provenance"))
        end
    end
    top!=0 && laterals==Set(original) || throw(ArgumentError(
        "$caller: B4 operation requires its copied cap and all Source laterals"))
    return boundary
end

function _extrude_nonew_region_boundary(m,t,source,caller;sides::Int=4,
                                        shared::Dict{Int,Vector{Int}}=Dict{Int,Vector{Int}}())
    boundary=abs.(_extrude_nonew_boundary_surfaces(m,t,caller))
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
        # `shared` carries the certified shared incidences from the scope
        # (surface => its complete sorted Region owner set); every other
        # multi-Region incidence keeps the original loud diagnostic.
        length(adjacent)==1 && adjacent[1]==t ||
            (haskey(shared,surface) && adjacent==shared[surface]) ||
            throw(ArgumentError(
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

# Which Source-cell side a shared wall's generatrix chain occupies: the
# generatrix's two endpoints match exactly one consecutive corner pair of
# the actual source cell, on either traversal direction. Every certified
# shared lateral's generatrix is a Source boundary curve, so exactly one
# side answers.
function _extrude_nonew_whole_side(m,cell::NTuple{4,Int32},source_mesh,
                                   generatrix::Int,s::Int,caller)
    coordinates=ntuple(k->ntuple(d->source_mesh.coords[d,cell[k]],3),4)
    chain=_extrude_curve_nodes(m,generatrix,caller)
    length(chain)==2 || _extrude_nonew_shared_planner_error(caller,s)
    for side in 1:4
        next=mod1(side+1,4)
        (chain[1]==coordinates[side] && chain[2]==coordinates[next]) ||
            (chain[2]==coordinates[side] && chain[1]==coordinates[next]) ||
            continue
        return side
    end
    _extrude_nonew_shared_planner_error(caller,s)
end

# Build the constant-size cap transition relation once for the actual lateral
# preferences. Lateral triangles are required; their orientation is adjustable.
# A tentative cap aligned to the final cap is lower priority than those
# preferences, matching the upstream retry that removes a trapping cap seed.
function _extrude_nonew_chain(intervals::Int,preferred,top_state::UInt8,caller;
                              forced=nothing)
    templates=_extrude_nonew_templates()
    length(templates)<=typemax(UInt16) || throw(ErrorException(
        "QuadTriNoNewVerts template table exceeds its packed index range"))
    function build_relation(constraints)
        relation=zeros(UInt16,3,3)
        costs=fill(typemax(Int),3,3)
        for (index,template) in enumerate(templates)
            all(k->template.faces[k]!=0,1:4) || continue
            # Each bound shared lateral accepts only the templates whose
            # matching face emits its transported geometric diagonal.
            all(constraint->template.faces[Int(constraint[1])]==Int(constraint[2]),
                constraints) || continue
            lower=Int(template.faces[5])+1
            upper=Int(template.faces[6])+1
            cost=4count(k->template.faces[k]!=preferred[k],1:4)
            cost+=template.faces[6]==top_state ? 0 : 1
            if cost<costs[lower,upper]
                costs[lower,upper]=cost
                relation[lower,upper]=UInt16(index)
            end
        end
        return relation,costs
    end
    relation,costs=build_relation(NTuple{2,UInt8}[])
    layered=nothing
    if forced!==nothing
        length(forced)==intervals || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared Region diagonal contract has inconsistent layers"))
        layered=[build_relation(forced[layer]) for layer in 1:intervals]
    end
    parents=zeros(UInt8,3,intervals)
    selected=zeros(UInt16,3,intervals)
    distance=(0,typemax(Int),typemax(Int))
    for layer in 1:intervals
        layer_relation,layer_costs=
            layered===nothing ? (relation,costs) : layered[layer]
        next_1=typemax(Int);next_2=typemax(Int);next_3=typemax(Int)
        for upper in 1:3
            layer==intervals && upper!=Int(top_state)+1 && continue
            best=typemax(Int)
            for lower in 1:3
                distance[lower]==typemax(Int) && continue
                index=layer_relation[lower,upper]
                index==0 && continue
                score=distance[lower]+layer_costs[lower,upper]
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

# Translate the lateral-wall diagonal edge positions each bound shared
# Surface emits (one unordered position pair per layer, per shared wall)
# into this Region's per-layer (side,state) template constraints: the same
# geometric edge must be the diagonal of each matching hex face on every
# layer. Each pair is a diagonal of exactly one lateral face, so each wall's
# side is unambiguous, and distinct shared walls occupy distinct source
# edges — a Region bounded by two partners may therefore carry two
# constraints on one layer.
function _extrude_nonew_shared_chain_contract(cols,cell::NTuple{4,Int32},
        intervals::Int,required_map,caller)
    forced=[NTuple{2,UInt8}[] for _ in 1:intervals]
    for (_,required) in required_map
        length(required)==intervals || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts shared Region diagonal contract has inconsistent layers"))
        side0=0
        for layer in 1:intervals
            req=required[layer]
            req===nothing && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared Region wall loses its bound diagonal"))
            v=_extrude_nonew_corners(cols,cell,layer)
            side=0;state=0
            for k in 1:4
                face=_EXTRUDE_NONEW_HEX_FACES[k]
                if _extrude_ekey(v[face[1]],v[face[3]])==req
                    side!=0 && throw(ArgumentError(
                        "$caller: QuadTriNoNewVerts shared Region diagonal is not unique to one wall"))
                    side=k;state=1
                elseif _extrude_ekey(v[face[2]],v[face[4]])==req
                    side!=0 && throw(ArgumentError(
                        "$caller: QuadTriNoNewVerts shared Region diagonal is not unique to one wall"))
                    side=k;state=2
                end
            end
            side==0 && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared Region wall does not reach its bound diagonal"))
            layer==1 && (side0=side)
            side==side0 || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared Region wall changes sides between layers"))
            push!(forced[layer],(UInt8(side),UInt8(state)))
        end
    end
    return forced
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
# A certified shared copied cap is a two-triangle Source whenever its
# transfinite root is — the cap's own Surface/Curves are generated copies,
# so the native controls bind on the root and each copied boundary chain
# must be its two endpoints.
function _extrude_nonew_cap_source(m::GeoModel,source::Int,sides::Int,
                                   caller::AbstractString)
    sides==4 || return false
    link=get(m.meshing.extrude_sources,(2,source),nothing)
    (link===nothing || link[1]!=2) && return false
    root=_extrude_quadtri_source_root(m,source)
    haskey(m.meshing.transfinite_surfaces,root) || return false
    for curve in _model_projection_surface_curves(m,source)
        chain=_extrude_curve_nodes(m,curve,caller)
        length(chain)==2 || return false
    end
    return !_model_surface_recombined(m,source)
end

function _extrude_nonew_two_tri_candidate(m::GeoModel,source::Int,sides::Int,
                                          caller::AbstractString)
    sides==4 || return false
    if !haskey(m.meshing.transfinite_surfaces,source)
        return _extrude_nonew_cap_source(m,source,sides,caller)
    end
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

# Whether `source`'s four-sided sweep resolves to the generic single-quad
# catalog — the only family that carries a transported per-side diagonal
# contract. Every bounded/strip/grid/triangle family answers false; an
# evaluation failure of a candidate probe conservatively counts as
# non-generic (the wall then stays whole and the dedup certificate decides).
function _extrude_nonew_generic_source(m::GeoModel,source::Int,
                                       caller::AbstractString)
    length(_model_projection_surface_curves(m,source))==4 || return false
    try
        _extrude_nonew_quad_patch_candidate(m,source,4,caller) && return false
        _extrude_nonew_quad_strip_candidate(m,source,4,caller) && return false
        _extrude_nonew_three_quad_strip_candidate(m,source,4,caller) && return false
        _extrude_nonew_four_quad_strip_candidate(m,source,4,caller) && return false
        _extrude_nonew_rect_grid_shape(m,source,4,caller)===nothing || return false
        _extrude_nonew_rect_grid_shape(m,source,4,caller;mode=:b4_strip)===nothing ||
            return false
        _extrude_nonew_two_tri_candidate(m,source,4,caller) && return false
    catch
        return false
    end
    return true
end

# A welded lateral wall may keep its split emission only when every
# NoNewVerts owner emits the wall split — the generic single-quad catalog
# carries a transported per-side diagonal contract, and every other
# catalog family still splits every boundary lateral quadrangle while it is
# not RecombLaterals (its per-quad subdivision diagonal is certified by the
# retained plan, so its emission either transports the contract or is
# compared bitwise by the shared-part dedup) — and every ordinary partner
# is unrecombined: the partner's prism subdivision emits the same diagonal
# triangles upstream keeps for it, and the bound wall's diagonals seed the
# global `SubdivideExtrudedMesh` edge set so both arrive on one facet set.
# A wall with a recombined ordinary/AddVerts partner or a RecombLaterals
# NoNewVerts owner stays in whole quadrangles upstream (the generic owners
# fan their layers to centroids, so the welded face is the pyramid base) —
# the same end the mixed-partner oracle reaches. The oracle also keeps the
# wall whole when `Recombine Surface{s}` constrains it directly — its
# owners fan — while leaving every other generated lateral's split
# emission alone, so a recombined wall joins the whole-quad side too.
function _extrude_nonew_wall_split(m::GeoModel,s::Int,owners,
                                   caller::AbstractString)
    _model_surface_recombined(m,s) && return false
    for owner in owners
        entry=_extrude_entity_params(m,3,owner,caller)
        entry===nothing && return false
        o_params,_,o_link=entry
        o_link[1]==2 || return false
        if _extrude_is_nonew(o_params)
            (_extrude_nonew_generic_source(m,abs(o_link[2]),caller) ||
             !o_params.recomb_laterals) || return false
        else
            o_params.recombine && return false
        end
    end
    return true
end

# Whether shared lateral `t` is certified between owners that all emit it
# split — the same rule `_extrude_nonew_wall_split` applies when the
# connected component plans the wall's transported contract, read off the
# surface's adjacent Regions instead of a validated owner set. Regions that
# use `t` as their sweep Source are not wall owners and never carry the
# lateral contract, so they are left out of the test.
function _extrude_nonew_split_wall(m::GeoModel,t::Int,caller::AbstractString)
    owners=Int[]
    for r in _extrude_quadtri_regions(m,t,caller)
        source=get(m.meshing.extrude_sources,(3,r),nothing)
        source!==nothing && abs(source[2])==t && continue
        push!(owners,r)
    end
    length(owners)>=2 || return false
    return _extrude_nonew_wall_split(m,t,owners,caller)
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
        size_field::Union{Nothing,AbstractSizeField}=nothing,_mesh_input=nothing,
        _shared_diagonals=nothing)
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
    sides=length(_model_projection_surface_curves(m,source))
    sides in (3,4) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts other source boundaries require the boundary-category planner"))
    # Certified shared incidences of the bounded translation slice (this
    # Region may bound several such walls and one stacked cap); every other
    # incidence still rejects.
    shared=_extrude_nonew_shared_for(m,t,caller)
    source_link=get(m.meshing.extrude_sources,(2,source),nothing)
    # A Source that is itself a generated Surface is admitted only as the
    # certified stacked cap — literally one partner's copied top re-used as
    # this Region's Source. Every other copied or chained source keeps the
    # dependency-planner diagnostic.
    source_link===nothing ||
        (source_link[1]==2 && haskey(shared,source) &&
         shared[source]==sort!(_extrude_quadtri_regions(m,source,caller))) ||
        throw(ArgumentError(
        "$caller: QuadTriNoNewVerts copied or chained sources require the dependency planner"))
    top,laterals=_extrude_nonew_region_boundary(m,t,source,caller;sides,shared)
    quad_patch=_extrude_nonew_quad_patch_candidate(m,source,sides,caller)
    quad_strip=!quad_patch && _extrude_nonew_quad_strip_candidate(m,source,sides,caller)
    three_quad_strip=!quad_patch && !quad_strip &&
        _extrude_nonew_three_quad_strip_candidate(m,source,sides,caller)
    four_quad_strip=!quad_patch && !quad_strip && !three_quad_strip &&
        _extrude_nonew_four_quad_strip_candidate(m,source,sides,caller)
    grid_shape=quad_patch || quad_strip || three_quad_strip || four_quad_strip ?
        nothing : _extrude_nonew_rect_grid_shape(m,source,sides,caller)
    rect_grid=grid_shape!==nothing
    strip_shape=quad_patch || quad_strip || three_quad_strip || four_quad_strip || rect_grid ?
        nothing : _extrude_nonew_rect_grid_shape(m,source,sides,caller;mode=:b4_strip)
    b4_strip=strip_shape!==nothing
    two_tri=!quad_patch && !quad_strip && !three_quad_strip && !four_quad_strip && !rect_grid && !b4_strip &&
        _extrude_nonew_two_tri_candidate(m,source,sides,caller)
    # A copied shared cap is not a lateral and leaves every certified
    # catalog free. A wall that stays in whole quadrangles (an
    # ordinary/AddVerts partner or a RecombLaterals NoNewVerts owner keeps
    # quads upstream) carries no contract either, so a non-generic owner
    # may free-emit it and let the shared-part dedup and joint incidence
    # certificate decide. A split wall is bounded to four-sided sources:
    # the generic catalog still carries the transported contract while a
    # non-generic owner either binds the wall — its own certified split
    # emission becomes the contract its partners reproduce — or free-emits
    # it as a later owner, where the dedup decides whether its boundary
    # diagonals coincide bitwise with the bound part.
    split_walls=[key for key in keys(shared) if
        get(m.meshing.extrude_sources,(2,key),(0,0))[1]==1 &&
        _extrude_nonew_wall_split(m,key,shared[key],caller)]
    isempty(split_walls) || sides==4 ||
        _extrude_nonew_shared_planner_error(caller,minimum(
            isempty(split_walls) ? collect(keys(shared)) : split_walls))
    grid_nodes,grid_cells_per_interval,grid_centers,grid_interval_centers=if rect_grid
        vertices,cells=_extrude_nonew_rect_grid_preflight(grid_shape,caller)
        (vertices,cells,0,0)
    elseif b4_strip && params.recomb_laterals
        vertices,cells,centers=_extrude_nonew_b4_strip_preflight(strip_shape,caller)
        (vertices,cells,centers,0)
    elseif b4_strip
        _ExtrudeNoNewB4Free.preflight(strip_shape,caller)
    else
        (0,0,0,0)
    end
    levels,refs=_extrude_nonew_levels(params,caller;
        extra_nodes_per_interval=grid_interval_centers,
        source_nodes=rect_grid || b4_strip ? grid_nodes : quad_patch ? 9 : quad_strip ? 6 : three_quad_strip ? 8 : four_quad_strip ? 10 : sides,
        extra_nodes=rect_grid || b4_strip ? grid_centers : four_quad_strip ? (params.recomb_laterals ? 4 : 0) :
            three_quad_strip ? (params.recomb_laterals ? 3 : 0) :
            quad_strip ? (params.recomb_laterals ? 2 : 0) :
            sides==4 && !two_tri && !quad_patch ? 1 : 0,
        cells_per_interval=rect_grid || b4_strip ? grid_cells_per_interval : four_quad_strip ? (params.recomb_laterals ? 28 : 24) :
            three_quad_strip ? (params.recomb_laterals ? 21 : 18) :
            quad_patch ? 24 : quad_strip ? (params.recomb_laterals ? 14 : 12) :
            two_tri && params.recomb_laterals ? 2 : sides==4 ? 6 : 3)
    # The copied cap was constructed from the Source emitter vector. Preserve
    # that explicit birth binding when the later Volume phase consumes the
    # completed Source; global Point/Curve allocation rank is a separate order.
    produced_source=_mesh_input isa _ModelMeshAuthority &&
        haskey(_mesh_input.producer_rows,(2,Int32(source))) ?
        _model_mesh_product(_mesh_input,2,source;caller,producer_order=true) : nothing
    source_mesh=produced_source===nothing ?
        mesh_model_surface(m,source;min_angle_deg=min_angle_deg,
            max_periodic_passes=max_periodic_passes,size_field=size_field,
            _mesh_input=_mesh_input,_boundary_regrade=false) : produced_source
    if quad_patch
        catalog=_extrude_nonew_quad_patch_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_quad_patch_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_quad_patch_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif quad_strip
        cols=_mesh_input===nothing ? nothing :
            _extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        actual_cap=_mesh_input===nothing ? nothing :
            _model_mesh_cap_input(_mesh_input,source_mesh,top,cols,caller)
        catalog=_extrude_nonew_quad_strip_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller;
            top_states=actual_cap===nothing ? nothing : actual_cap.top_states)
        cols===nothing && (cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller))
        _extrude_nonew_quad_strip_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_quad_strip_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif three_quad_strip
        cols=_mesh_input===nothing ? nothing :
            _extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        actual_cap=_mesh_input===nothing ? nothing :
            _model_mesh_cap_input(_mesh_input,source_mesh,top,cols,caller)
        catalog=_extrude_nonew_three_quad_strip_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller;
            top_states=actual_cap===nothing ? nothing : actual_cap.top_states)
        cols===nothing && (cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller))
        _extrude_nonew_three_quad_strip_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_three_quad_strip_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif four_quad_strip
        cols=_mesh_input===nothing ? nothing :
            _extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        actual_cap=_mesh_input===nothing ? nothing :
            _model_mesh_cap_input(_mesh_input,source_mesh,top,cols,caller)
        catalog=_extrude_nonew_four_quad_strip_catalog(m,source,source_mesh,spec,
            levels,refs,params.recomb_laterals,caller;
            top_states=actual_cap===nothing ? nothing : actual_cap.top_states)
        cols===nothing && (cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller))
        _extrude_nonew_four_quad_strip_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_four_quad_strip_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif rect_grid
        source_data=_extrude_nonew_rect_grid_source(m,source,source_mesh,spec,caller;_mesh_input)
        source_data.grid_shape==grid_shape || _extrude_nonew_rect_grid_error(caller,
            "actual source dimensions differ from their preflight")
        catalog=_extrude_nonew_rect_grid_plan(source_data,levels,refs,
            params.recomb_laterals,caller)
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        _extrude_nonew_rect_grid_product_certify(cols,catalog,spec,caller)
        return _extrude_nonew_rect_grid_finish(m,t,params,source,source_mesh,
            catalog,cols,top,laterals,caller)
    elseif b4_strip
        source_data=_extrude_nonew_rect_grid_source(m,source,source_mesh,spec,caller;
            mode=:b4_strip,_mesh_input)
        source_data.grid_shape==strip_shape || _extrude_nonew_b4_strip_error(caller,
            "actual source dimensions differ from their preflight")
        cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
        actual_cap=_model_mesh_cap_input(_mesh_input,source_mesh,top,cols,caller)
        top_states=actual_cap===nothing ? nothing : actual_cap.top_states
        retained=params.recomb_laterals && _extrude_nonew_unmerged_b4_incidence(m,t,caller)
        catalog=retained ?
            _extrude_nonew_retained_b4_plan(source_data,levels,refs,caller;top_states) :
            params.recomb_laterals ?
                _extrude_nonew_b4_strip_plan(source_data,levels,refs,true,caller;top_states) :
                _ExtrudeNoNewB4Free.plan(source_data,levels,refs,caller;top_states)
        _extrude_nonew_rect_grid_product_certify(cols,catalog,spec,caller)
        if retained
            return _extrude_nonew_retained_b4_finish(m,t,params,source,source_mesh,
                catalog,cols,top,laterals,caller)
        elseif params.recomb_laterals
            return _extrude_nonew_b4_strip_finish(m,t,params,source,source_mesh,
                catalog,cols,top,laterals,caller)
        end
        return _ExtrudeNoNewB4Free.finish(m,t,params,source,source_mesh,
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
    # A shared lateral whose certified owner set contains a non-NoNewVerts
    # Region stays in whole quadrangles — an ordinary layered or AddVerts
    # partner sweeps quad boundary faces upstream keeps for it — so this
    # Region must not put a diagonal on that wall. No template can (the
    # chain requires every lateral face split); such layers fan to their
    # centroid, which lands the wall face on a pyramid's whole-quad base.
    whole_sides=Set{Int}()
    for (s,owners) in shared
        s_link=get(m.meshing.extrude_sources,(2,s),nothing)
        (s_link===nothing || s_link[1]!=1) && continue
        _extrude_nonew_wall_split(m,s,owners,caller) && continue
        side=_extrude_nonew_whole_side(m,cell,source_mesh,
            abs(s_link[2]),s,caller)
        side in whole_sides && _extrude_nonew_shared_planner_error(caller,s)
        push!(whole_sides,side)
    end
    top_state=UInt8(_extrude_nonew_top_diagonal(cell) ? 1 : 2)
    states=Matrix{UInt8}(undef,6,intervals)
    problems=falses(intervals)
    templates=_extrude_nonew_templates()
    picks=zeros(UInt16,intervals)
    if params.recomb_laterals
        fill!(states,0)
        states[6,intervals]=top_state
        # A transported shared-wall diagonal cannot land on a whole
        # RecombLaterals quad — symmetric parameters keep this defensive.
        _shared_diagonals===nothing ||
            all(required->all(isnothing,required),values(_shared_diagonals)) ||
            throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared Region wall keeps a diagonal its partner does not emit"))
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
        forced=_shared_diagonals===nothing ? nothing :
            _extrude_nonew_shared_chain_contract(
                cols,cell,intervals,_shared_diagonals,caller)
        if isempty(whole_sides)
            picks=_extrude_nonew_chain(intervals,preferred,top_state,caller;
                forced=forced)
            for layer in 1:intervals,side in 1:6
                states[side,layer]=templates[picks[layer]].faces[side]
            end
        else
            # Whole-quad walls cannot come from the template chain (every
            # lateral face must carry a diagonal): every layer fans to its
            # centroid, so the shared side stays the pyramid's whole-quad
            # base. A transported diagonal on another shared wall still
            # binds that side; free laterals take their preferred diagonal;
            # the bottom cap stays whole on its single-quad Source and
            # interior caps stay whole between consecutive fans.
            fill!(problems,true)
            for layer in 1:intervals
                for side in 1:4
                    states[side,layer]=side in whole_sides ? 0 : preferred[side]
                end
                states[5,layer]=0
                states[6,layer]=layer==intervals ? top_state : 0
            end
            forced===nothing || for layer in 1:intervals
                for (side,state) in forced[layer]
                    Int(side) in whole_sides && throw(ArgumentError(
                        "$caller: QuadTriNoNewVerts shared Region wall " *
                        "keeps both a whole-quad and a diagonal contract"))
                    states[Int(side),layer]=state
                end
            end
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
include("ModelExtrudeNoNewB4Strip.jl")
include("ModelExtrudeNoNewB4Free.jl")
