# Operation-owned native surface geometry and producer identities. The actual
# lower-dimensional products determine identity; coordinates only verify an
# existing identity and evaluate its Surface parametrization.
import ..Mesh2D
using ..SizeField: AbstractAnisoField, BoundedSize, GMSH_MAX_SIZE
import ..SizeField: metric_at
struct _GmshNativePlaneFrame
    origin::NTuple{3,Float64}
    first::NTuple{3,Float64}
    second::NTuple{3,Float64}
    normal::NTuple{3,Float64}
    occ::Bool
end
_GmshNativePlaneFrame(origin,first,second,normal)=
    _GmshNativePlaneFrame(origin,first,second,normal,false)
include("ModelPlaneSmallSVD.jl")

function _gmsh_native_plane_frame(samples,tolerance,characteristic_length)
    length(samples)>=3 || return nothing
    anchor=samples[1]
    edge=Mesh2D._gmsh_sub3(samples[2],anchor)
    threshold=sqrt(tolerance)*characteristic_length
    for position in 3:length(samples)
        other=Mesh2D._gmsh_sub3(samples[position],anchor)
        # SVector3::crossprod used by computeMeanPlane has this exact second
        # component expression. prodve below uses its separate literal body.
        normal=(edge[2]*other[3]-edge[3]*other[2],
            -(edge[1]*other[3]-other[1]*edge[3]),
            edge[1]*other[2]-edge[2]*other[1])
        Mesh2D._gmsh_norm3(normal)>threshold || continue
        reference=normal[1]==0.0 ? (1.0,0.0,0.0) :
                  normal[2]==0.0 ? (0.0,1.0,0.0) : (0.0,0.0,1.0)
        first=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(normal,reference))
        second=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(first,normal))
        rhs=Mesh2D._gmsh_dot3(anchor,normal)
        dominant=abs(normal[1])>=abs(normal[2]) && abs(normal[1])>=abs(normal[3]) ? 1 :
            abs(normal[2])>=abs(normal[1]) && abs(normal[2])>=abs(normal[3]) ? 2 : 3
        origin=ntuple(axis->axis==dominant ? rhs/normal[dominant] : 0.0,3)
        all(isfinite,(origin...,first...,second...,normal...)) || throw(ArgumentError(
            "native surface mesher: mean Plane frame is not finite"))
        return _GmshNativePlaneFrame(origin,first,second,normal)
    end
    return nothing
end

@inline function _gmsh_native_surface_parameters(frame::_GmshNativePlaneFrame,point)
    offset=Mesh2D._gmsh_sub3(point,frame.origin)
    return (Mesh2D._gmsh_dot3(offset,frame.first),
            Mesh2D._gmsh_dot3(offset,frame.second))
end

@inline function _surface_metric_lift(frame::_GmshNativePlaneFrame,u,v)
    frame.occ && return ntuple(axis->
        fma(u,frame.first[axis],v*frame.second[axis])+frame.origin[axis],3)
    return ntuple(axis->frame.origin[axis]+frame.first[axis]*u+frame.second[axis]*v,3)
end

function _gmsh_native_surface_curves(m::GeoModel,tag::Int)
    curves=Int[];seen=Set{Int}()
    for loop in m.surfaces[tag],signed in m.loops[loop]
        curve=abs(signed)
        curve in seen && continue
        push!(curves,curve);push!(seen,curve)
    end
    return curves
end

function _gmsh_native_curve_sample(m::GeoModel,curve::Int,parameter,caller)
    if _curve_type(m,curve)===:line && _occ_geometry(m,curve)===nothing
        a,b=m.curves[curve];first=m.points[a];last=m.points[b]
        return ntuple(axis->first[axis]+parameter*(last[axis]-first[axis]),3)
    end
    return _model_curve_point(m,curve,parameter,caller)
end

function _gmsh_native_surface_frame(m::GeoModel,tag::Int,options,caller)
    geometry=get(m.surface_geometry,tag,nothing)
    if geometry!==nothing && hasproperty(geometry,:occ)
        geometry.occ===:plane || throw(ArgumentError(
            "$caller: native Surface[$tag] requires its actual Plane parametrization"))
        return _GmshNativePlaneFrame(geometry.center,geometry.X,geometry.Y,geometry.axis,true)
    end
    receipt=get(m.meshing.native_geo_kernel.plane_frames,tag,nothing)
    receipt===nothing || return _GmshNativePlaneFrame(
        receipt.origin,receipt.first,receipt.second,receipt.normal)
    samples=NTuple{3,Float64}[]
    for curve in _gmsh_native_surface_curves(m,tag)
        low,high=_model_curve_param_bounds(m,curve,caller)
        push!(samples,_gmsh_native_curve_sample(m,curve,low+0.333*(high-low),caller))
        push!(samples,_gmsh_native_curve_sample(m,curve,low+0.666*(high-low),caller))
    end
    frame=_gmsh_native_plane_frame(samples,options.geometry_tolerance,options.ctx_lc)
    frame===nothing && throw(ArgumentError(
        "$caller: native mean Plane SVD fallback for Surface[$tag] is not implemented"))
    return frame
end

function _gmsh_native_surface_canonical_orientation(m,tag,mesh,physical,caller)
    geometry=get(m.surface_geometry,tag,nothing)
    geometry!==nothing && hasproperty(geometry,:occ) && geometry.occ===:plane || return mesh
    cells=if mesh isa Mesh
        mesh.tris
    else
        block=findfirst(block->msh_dimension(block.msh)==2 && !isempty(block.nodes),mesh.blocks)
        block===nothing && return mesh
        mesh.blocks[block].nodes
    end
    isempty(cells) && return mesh
    normal=_occ_surface_normal(geometry,0.0,0.0)
    a=ntuple(axis->physical[axis,cells[1,1]],3)
    b=ntuple(axis->physical[axis,cells[2,1]],3)
    c=ntuple(axis->physical[axis,cells[3,1]],3)
    product=_model_orientation_dot(_model_orientation_cross(
        _model_orientation_sub(b,a),_model_orientation_sub(c,a)),normal)
    return product>0 ? mesh : _reversed_surface_mesh(mesh,caller,tag)
end

struct _ModelSurfaceSource
    x::Vector{Float64}
    y::Vector{Float64}
    tags::Vector{UInt64}
    owners::Vector{Tuple{Int,Int32}}
    parameters::Vector{Union{Nothing,Vector{Float64}}}
    world::Vector{NTuple{3,Float64}}
    positions::Dict{UInt64,Int}
    occurrence::Vector{UInt64}
    segments::Vector{NTuple{3,Int}}
    boundaries::Vector{Tuple{Int,Int}}
    internals::Vector{Tuple{Int,Int}}
    point_positions::Dict{Int,Int}
    # Source rows minted by forced periodic subdivision — actual nodes the
    # surface producer creates for a Curve, not retained lower products.
    forced_rows::Set{Int}
end

_ModelSurfaceSource()=_ModelSurfaceSource(Float64[],Float64[],UInt64[],
    Tuple{Int,Int32}[],Union{Nothing,Vector{Float64}}[],NTuple{3,Float64}[],
    Dict{UInt64,Int}(),UInt64[],NTuple{3,Int}[],Tuple{Int,Int}[],
    Tuple{Int,Int}[],Dict{Int,Int}(),Set{Int}())

function _gmsh_native_surface_add_source!(source,m,mesh,row,frame,caller;
                                         surface_tag,ctx_lc)
    data=mesh.entity_data
    tag=data.external_node_tags[row]
    point=ntuple(axis->mesh.coords[axis,row],3)
    owner=data.node_entities[row]
    previous=get(source.positions,tag,0)
    if previous!=0
        isequal(source.world[previous],point) && source.owners[previous]==owner ||
            throw(ArgumentError("$caller: actual lower node $tag has conflicting identity data"))
        return previous
    end
    all(isfinite,point) || throw(ArgumentError("$caller: actual lower node $tag is not finite"))
    parameter=data.node_parametric[row]
    # Native Curve vertices are reparametrized through their stored Curve U.
    # Their actual mesh coordinates can have been edited independently of U.
    # Native GEO Points on a Plane and discrete Curve vertices instead use
    # their actual coordinates (MVertex.cpp reparamMeshVertexOnFace).
    u,v=if frame.occ && owner[1]==1 && haskey(m.curves,Int(owner[2]))
        value=parameter===nothing || isempty(parameter) ? 0.0 : parameter[1]
        _occ_edge_reparam_on_face(m,m.surface_geometry[surface_tag],
            Int(owner[2]),value,1,ctx_lc,caller)
    elseif frame.occ && owner[1]==0 && haskey(m.points,Int(owner[2]))
        _occ_vertex_reparam_on_face(m,surface_tag,m.surface_geometry[surface_tag],
            Int(owner[2]),1,ctx_lc)
    elseif owner[1]==1 && haskey(m.curves,Int(owner[2]))
        # Original MVertex instances can be associated with a Curve by P2
        # pruning. Their base getParameter writes zero and returns false;
        # native reparamOnFace consumes that scalar before its warning.
        value=parameter===nothing || isempty(parameter) ? 0.0 : parameter[1]
        _gmsh_native_surface_parameters(frame,
            _gmsh_native_curve_sample(m,Int(owner[2]),value,caller))
    else
        _gmsh_native_surface_parameters(frame,point)
    end
    isfinite(u) && isfinite(v) || throw(ArgumentError(
        "$caller: actual lower node $tag has non-finite Surface parameters"))
    push!(source.x,u);push!(source.y,v);push!(source.tags,tag)
    push!(source.owners,owner);push!(source.world,point)
    push!(source.parameters,parameter===nothing ? nothing : copy(parameter))
    position=length(source.tags);source.positions[tag]=position
    owner[1]==0 && (source.point_positions[Int(owner[2])]=position)
    return position
end

function _gmsh_native_surface_curve!(source,m,input,curve,frame,internal,caller;
                                    surface_tag,ctx_lc)
    mesh=_model_mesh_product(input,1,curve;caller)
    mesh!==nothing && mesh.entity_data!==nothing || throw(ArgumentError(
        "$caller: native Surface input lacks actual Curve[$curve] identities"))
    visited=Int32[]
    for block in mesh.blocks
        block.msh==1 || throw(ArgumentError(
            "$caller: native Curve[$curve] source must have linear Line cells"))
        append!(visited,vec(block.nodes))
    end
    isempty(visited) && throw(ArgumentError("$caller: native Curve[$curve] source has no Line cells"))
    unique!(visited)
    sort!(visited;by=row->mesh.entity_data.external_node_tags[row])
    mapping=zeros(Int,nnodes(mesh))
    for row in visited
        mapping[row]=_gmsh_native_surface_add_source!(source,m,mesh,row,frame,caller;
            surface_tag,ctx_lc)
        push!(source.occurrence,mesh.entity_data.external_node_tags[row])
    end
    for block in mesh.blocks,column in axes(block.nodes,2)
        a,b=mapping[block.nodes[:,column]]
        a!=b || throw(ArgumentError("$caller: actual Curve[$curve] has a repeated Line identity"))
        push!(source.segments,(a,b,curve))
        push!(internal ? source.internals : source.boundaries,(a,b))
    end
    return nothing
end

# A periodic relation can synchronize an interior parameter onto a curve
# whose actual lower product has no vertex there — a translated slave copy
# of a master's embedded split is the common producer. `_surface_pslg`
# evaluates such forced parameters into PSLG vertices; the native source
# materializes the same vertices here so the recovered cells carry matching
# slave/master chains. A forced parameter resolving onto a stored Point row
# reuses that vertex — the same coordinate dedup `_add_surface_curve_point!`
# performs against model points.
function _gmsh_native_surface_forced_vertex!(source::_ModelSurfaceSource,m,
        curve::Int,parameter::Float64,frame,next_tag,caller;
        surface_tag,ctx_lc,_periodic_positions=nothing)
    point=_periodic_forced_curve_point(m,curve,parameter,
        _periodic_positions,caller)
    scale=max(1.0,hypot(point...))
    tolerance=128eps(Float64)*scale
    matched=0;matched_tag=0
    for (point_tag,position) in source.point_positions
        haskey(m.points,point_tag) || continue
        candidate=m.points[point_tag]
        hypot(candidate[1]-point[1],candidate[2]-point[2],
              candidate[3]-point[3])<=tolerance || continue
        if matched!=0 && matched!=position
            throw(ArgumentError(
                "$caller: Curve[$curve] forced subdivision coincides with " *
                "Points[$matched_tag] and [$point_tag]"))
        end
        matched=position;matched_tag=point_tag
    end
    matched!=0 && return matched
    u,v=if frame.occ && haskey(m.curves,curve)
        _occ_edge_reparam_on_face(m,m.surface_geometry[surface_tag],
            curve,parameter,1,ctx_lc,caller)
    else
        _gmsh_native_surface_parameters(frame,point)
    end
    isfinite(u) && isfinite(v) || throw(ArgumentError(
        "$caller: Curve[$curve] forced subdivision is not finite on the " *
        "Surface"))
    next_tag===nothing && throw(ArgumentError(
        "$caller: forced Curve subdivision needs a node label allocator"))
    next_tag[]==typemax(UInt64) && throw(ArgumentError(
        "$caller: forced Curve subdivision exhausts UInt64 labels"))
    label=next_tag[]+UInt64(1)
    next_tag[]=label
    push!(source.x,u);push!(source.y,v);push!(source.tags,label)
    push!(source.owners,(1,Int32(curve)));push!(source.world,point)
    push!(source.parameters,Float64[parameter])
    position=length(source.tags);source.positions[label]=position
    push!(source.occurrence,label)
    push!(source.forced_rows,position)
    return position
end

# Expand every constrained source segment whose owning Curve carries a
# forced parameter strictly inside the segment's parameter span. The forced
# lists live in the curve's native frame and already include both bounds, so
# interior values land in exactly one cell; each fragment keeps the
# segment's owner and its boundary/internal classification.
function _gmsh_native_surface_forced_subdivide!(source::_ModelSurfaceSource,m,
        forced,frame,next_tag,caller;surface_tag,ctx_lc,
        _periodic_positions=nothing)
    (forced===nothing || isempty(forced)) && return nothing
    function endpoint_parameter(position,curve,t0,t1,at_start)
        owner=source.owners[position]
        if owner[1]==1 && Int(owner[2])==curve
            parameter=source.parameters[position]
            parameter!==nothing && !isempty(parameter) && return parameter[1]
        end
        a,b=m.curves[curve]
        # A closed curve's wrap vertex is both endpoints: it leaves at t0 as a
        # segment's start and arrives at t1 as its end. Resolving both to t0
        # would re-insert every earlier interior vertex in the closing cell.
        position==get(source.point_positions,a,-1) &&
            return (a==b && !at_start) ? t1 : t0
        position==get(source.point_positions,b,-1) && return t1
        return nothing
    end
    chains=Dict{Tuple{Int,Int},Vector{Int}}()
    for (a,b,owner) in source.segments
        haskey(chains,(a,b)) && continue
        parameters=get(forced,owner,nothing)
        degenerated=owner in m.meshing.degenerated
        chain=Int[a]
        if parameters!==nothing && !degenerated
            t0,t1=_model_curve_param_bounds(m,owner,caller)
            pa=endpoint_parameter(a,owner,t0,t1,true)
            pb=endpoint_parameter(b,owner,t0,t1,false)
            if pa!==nothing && pb!==nothing
                tolerance=128eps(Float64)*max(1.0,abs(t1-t0))
                lo,hi=minmax(pa,pb)
                interior=sort!([u for u in parameters if
                                lo+tolerance<u<hi-tolerance];rev=pa>pb)
                for parameter in interior
                    push!(chain,_gmsh_native_surface_forced_vertex!(
                        source,m,owner,parameter,frame,next_tag,caller;
                        surface_tag,ctx_lc,_periodic_positions))
                end
            end
        end
        push!(chain,b)
        chains[(a,b)]=chain
    end
    segments=NTuple{3,Int}[]
    for (a,b,owner) in source.segments
        chain=chains[(a,b)]
        for index in 1:(length(chain)-1)
            push!(segments,(chain[index],chain[index+1],owner))
        end
    end
    for records in (source.boundaries,source.internals)
        expanded=Tuple{Int,Int}[]
        for (a,b) in records
            chain=get(chains,(a,b),nothing)
            chain===nothing && (chain=Int[a,b])
            for index in 1:(length(chain)-1)
                push!(expanded,(chain[index],chain[index+1]))
            end
        end
        empty!(records);append!(records,expanded)
    end
    copy!(resize!(source.segments,length(segments)),segments)
    return nothing
end

function _gmsh_native_surface_source(m::GeoModel,tag::Int,input,frame,caller;
                                     ctx_lc=0.0,forced=nothing,
                                     next_tag=nothing,
                                     _periodic_positions=nothing)
    source=_ModelSurfaceSource()
    boundaries=_gmsh_native_surface_curves(m,tag)
    for curve in boundaries
        _gmsh_native_surface_curve!(source,m,input,curve,frame,false,caller;
            surface_tag=tag,ctx_lc)
    end
    points,curves=_model_surface_embedding_tags(m,tag,caller)
    isempty(intersect(Set(boundaries),Set(curves))) || throw(ArgumentError(
        "$caller: a Surface boundary Curve cannot also be embedded"))
    for curve in curves
        _gmsh_native_surface_curve!(source,m,input,curve,frame,true,caller;
            surface_tag=tag,ctx_lc)
    end
    for point in points
        mesh=_model_mesh_product(input,0,point;caller)
        mesh!==nothing && mesh.entity_data!==nothing || throw(ArgumentError(
            "$caller: embedded Point[$point] lacks actual source identities"))
        rows=unique(Int32[node for block in mesh.blocks for node in block.nodes])
        length(rows)==1 || throw(ArgumentError(
            "$caller: embedded Point[$point] needs an explicit unique actual producer node"))
        row=only(rows)
        position=_gmsh_native_surface_add_source!(source,m,mesh,row,frame,caller;
            surface_tag=tag,ctx_lc)
        source.point_positions[point]=position
        push!(source.occurrence,mesh.entity_data.external_node_tags[row])
    end
    _gmsh_native_surface_forced_subdivide!(source,m,forced,frame,next_tag,
        caller;surface_tag=tag,ctx_lc,_periodic_positions)
    return source
end

# A source vertex can sit exactly on a supplied segment — an embedded Point
# on an embedded Curve chain is the common producer. Like
# `Mesh2D._insert_segment!`, the constraint then becomes the ordered chain of
# recovered sub-edges; every fragment keeps the segment's owner and its
# boundary/internal classification.
function _gmsh_native_surface_split_constraints!(source::_ModelSurfaceSource)
    chains=Dict{Tuple{Int,Int},Vector{Int}}()
    scratch=Mesh2D._GmshOrientScratch()
    for (a,b,_) in source.segments
        haskey(chains,(a,b)) && continue
        ax,ay=source.x[a],source.y[a]
        bx,by=source.x[b],source.y[b]
        dx=bx-ax;dy=by-ay
        hits=Tuple{Float64,Int}[]
        if dx!=0.0 || dy!=0.0
            lo_x,hi_x=minmax(ax,bx);lo_y,hi_y=minmax(ay,by)
            denom=dx*dx+dy*dy
            for v in eachindex(source.x)
                (v==a || v==b) && continue
                vx,vy=source.x[v],source.y[v]
                (lo_x<=vx<=hi_x && lo_y<=vy<=hi_y) || continue
                Mesh2D._gmsh_orient2_value((ax,ay),(bx,by),(vx,vy),scratch)==
                    0.0 || continue
                push!(hits,(((vx-ax)*dx+(vy-ay)*dy)/denom,v))
            end
            sort!(hits;by=first,alg=Base.Sort.MergeSort)
        end
        chain=Int[a]
        last_x,last_y=ax,ay
        for (_,v) in hits
            # A distinct source row sharing an endpoint's or an earlier hit's
            # coordinates splits nothing — a zero-length fragment is not a
            # legal constraint, and `b` stays the chain's real endpoint.
            (source.x[v]==bx && source.y[v]==by) && continue
            (source.x[v]==last_x && source.y[v]==last_y) && continue
            push!(chain,v)
            last_x,last_y=source.x[v],source.y[v]
        end
        push!(chain,b)
        chains[(a,b)]=chain
    end
    isempty(chains) && return nothing
    segments=NTuple{3,Int}[]
    for (a,b,owner) in source.segments
        chain=chains[(a,b)]
        for index in 1:(length(chain)-1)
            push!(segments,(chain[index],chain[index+1],owner))
        end
    end
    for (records,) in ((source.boundaries,),(source.internals,))
        expanded=Tuple{Int,Int}[]
        for (a,b) in records
            chain=get(chains,(a,b),nothing)
            if chain===nothing
                push!(expanded,(a,b))
            else
                for index in 1:(length(chain)-1)
                    push!(expanded,(chain[index],chain[index+1]))
                end
            end
        end
        empty!(records);append!(records,expanded)
    end
    copy!(resize!(source.segments,length(segments)),segments)
    return nothing
end

# The constrained segments surviving `_gmsh_native_surface_initial` are the
# recovered dim-1 chains: forced periodic subdivisions and on-segment vertex
# splits insert vertices between the consumed cells' endpoints. Emit each
# owning Curve's ordered actual UID chain so the merged dim-1 product can
# receive the subdivision upstream's meshGFace writes back onto each GEdge.
function _gmsh_native_surface_curve_chains(source::_ModelSurfaceSource)
    edges=Dict{Int32,Vector{Tuple{Int,Int}}}()
    order=Int32[]
    for (a,b,owner) in source.segments
        key=Int32(abs(owner))
        if !haskey(edges,key)
            edges[key]=Tuple{Int,Int}[]
            push!(order,key)
        end
        push!(edges[key],(a,b))
    end
    chains=Dict{Int32,Vector{UInt64}}()
    for key in order
        adjacency=Dict{Int,Int}()
        incoming=Set{Int}()
        consistent=true
        for (a,b) in edges[key]
            haskey(adjacency,a) && (consistent=false)
            adjacency[a]=b
            push!(incoming,b)
        end
        consistent || continue
        heads=Int[a for (a,b) in edges[key] if !(a in incoming)]
        length(heads)>1 && continue
        vertex=isempty(heads) ? first(edges[key])[1] : only(heads)
        start=vertex
        chain=UInt64[source.tags[vertex]]
        for _ in eachindex(edges[key])
            vertex=get(adjacency,vertex,0)
            vertex==0 && break
            push!(chain,source.tags[vertex])
            vertex==start && break
        end
        length(chain)==length(edges[key])+1 || continue
        chains[key]=chain
    end
    return chains
end

function _gmsh_native_surface_initial(source::_ModelSurfaceSource,tag::Int;
                                     cleanup_faces=false)
    _gmsh_native_surface_split_constraints!(source)
    rows=Mesh2D._gmsh_initial_uid_order(source.occurrence)
    order=Int[source.positions[tag] for tag in source.occurrence[rows]]
    product=Mesh2D._gmsh_poly_initial(source.x,source.y,source.segments,order;
        # PolyMesh payload is an operation row. Actual UInt64 public labels
        # determine insertion order above and remain separate for MTri ties.
        vertex_labels=1:length(source.tags),surface_tag=tag)
    return Mesh2D._gmsh_poly_triangulation(product,source.boundaries,source.internals;
        cleanup_surface=cleanup_faces ? tag : nothing)
end

struct _ModelSurfaceProduction{M}
    mesh::M
    node_tags::Vector{UInt64}
    node_owners::Vector{Tuple{Int,Int32}}
    node_parameters::Vector{Union{Nothing,Vector{Float64}}}
    retained_nodes::BitVector
    node_high_water::UInt64
    element_tags::Vector{UInt64}
    element_high_water::UInt64
    # Actual constructor traversal can differ from the tensor's row storage.
    birth_rows::Vector{Int32}
    constructor_mesh::Union{Nothing,Mesh,MixedMesh}
    source_node_updates::Dict{UInt64,NTuple{3,Float64}}
    native_recombine_events::Vector{Tuple{Symbol,String}}
    # Ordered actual node UIDs of each consumed Curve chain after forced
    # subdivisions and on-segment vertex splits — the meshGFace write-back
    # payload the merged dim-1 product receives when the recovered chain
    # grew past its consumed elements.
    curve_chains::Dict{Int32,Vector{UInt64}}
end
_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,
    element_tags,element_high_water,birth_rows,constructor_mesh,source_node_updates,
    native_recombine_events)=_ModelSurfaceProduction(mesh,tags,owners,parameters,
        retained,high_water,element_tags,element_high_water,birth_rows,constructor_mesh,
        source_node_updates,native_recombine_events,Dict{Int32,Vector{UInt64}}())
_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,
    element_tags,element_high_water,birth_rows,constructor_mesh,source_node_updates)=
    _ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,
        element_tags,element_high_water,birth_rows,constructor_mesh,source_node_updates,
        Tuple{Symbol,String}[])
_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,
    element_tags,element_high_water,birth_rows,constructor_mesh)=_ModelSurfaceProduction(
        mesh,tags,owners,parameters,retained,high_water,element_tags,element_high_water,
        birth_rows,constructor_mesh,Dict{UInt64,NTuple{3,Float64}}())
_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,
    element_tags,element_high_water,birth_rows)=_ModelSurfaceProduction(mesh,tags,owners,
        parameters,retained,high_water,element_tags,element_high_water,birth_rows,nothing)
_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,
    element_tags,element_high_water)=_ModelSurfaceProduction(mesh,tags,owners,
        parameters,retained,high_water,element_tags,element_high_water,
        Int32.(eachindex(tags)))
_ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water)=
    _ModelSurfaceProduction(mesh,tags,owners,parameters,retained,high_water,UInt64[],UInt64(0))

function _gmsh_native_surface_export(T,inside,numbers,source::_ModelSurfaceSource,
                                      frame,tag::Int,high_water;projected=false,trace=nothing)
    length(numbers)==T.nreal || throw(ArgumentError(
        "native Surface producer labels do not align with actual vertices"))
    used=falses(T.nreal)
    cells=NTuple{3,Int32}[]
    order=Int[triangle for triangle in eachindex(T.alive) if inside[triangle]]
    if trace!==nothing
        sort!(order;lt=(a,b)->begin
            first=trace.radii[Mesh2D._gmsh_triangle_key(T,a)]
            second=trace.radii[Mesh2D._gmsh_triangle_key(T,b)]
            first>second && return true
            first<second && return false
            return Mesh2D._gmsh_triangle_priority_numbers(numbers,T,a)<
                   Mesh2D._gmsh_triangle_priority_numbers(numbers,T,b)
        end)
    end
    for triangle in order
        nodes=Mesh2D._gmsh_surface_triangle_nodes(T,triangle,trace)
        used[collect(nodes)].=true;push!(cells,nodes)
    end
    rows=Int32.(findall(used));mapping=zeros(Int32,T.nreal)
    mapping[rows]=Int32.(1:length(rows))
    coords=Matrix{Float64}(undef,3,length(rows))
    owners=Vector{Tuple{Int,Int32}}(undef,length(rows))
    parameters=Vector{Union{Nothing,Vector{Float64}}}(undef,length(rows))
    retained=falses(length(rows))
    for (position,row) in enumerate(rows)
        if row<=length(source.tags)
            coords[:,position].=source.world[row]
            owners[position]=source.owners[row]
            parameters[position]=source.parameters[row]
            retained[position]=!(row in source.forced_rows)
        else
            coords[:,position].=_surface_metric_lift(frame,T.x[row],T.y[row])
            owners[position]=(2,Int32(tag));parameters[position]=[T.x[row],T.y[row]]
        end
        if projected
            coords[1,position]=T.x[row];coords[2,position]=T.y[row]
            coords[3,position]=0.0
        end
    end
    triangles=Matrix{Int32}(undef,3,length(cells))
    for (column,nodes) in enumerate(cells),slot in 1:3
        triangles[slot,column]=mapping[nodes[slot]]
    end
    labels=trace===nothing ? UInt64[] :
        UInt64[trace.elements[Mesh2D._gmsh_triangle_key(T,triangle)] for triangle in order]
    element_max=trace===nothing ? UInt64(0) : trace.element_high_water
    return _ModelSurfaceProduction(Mesh(coords;tris=triangles),numbers[rows],owners,
        parameters,retained,UInt64(high_water),labels,element_max)
end

function _gmsh_native_surface_lower_input(m,options,input,high_water,caller)
    if input!==nothing
        input.mesh!==nothing || throw(ArgumentError(
            "$caller: native Surface input has no actual lower-dimensional products"))
        maximum_tag=maximum(input.mesh.entity_data.external_node_tags;init=UInt64(0))
        limit=high_water===nothing ? maximum_tag : UInt64(high_water)
        maximum_tag<=limit || throw(ArgumentError(
            "$caller: native Surface high water omits a consumed actual node label"))
        return input,limit
    end
    # Direct Model callers also use the real detached Point/Curve producer.
    # The registered allocator does not read or publish a live API session.
    # Like upstream `meshGFace`, the plan consumes each already-meshed Curve's
    # stored `mesh_vertices` (`m.curve_params`) verbatim — only a Curve still
    # missing its discretization is graded, and then stored once. Regrading
    # stored parameters would re-derive them from an `options.ctx_lc` whose
    # bounding box samples the stored nodes themselves, so every Surface call
    # would ratchet the shared-edge subdivision and crack the exact weld.
    plan=_model_mesh_generate_higher(m,nothing,Tuple{Int,Int,Mesh}[],
        Vector{Int32}[],Tuple{Int,Int32}[],1,options,caller;
        element_order=1,renumber=false,reuse_curve_params=true)
    merge!(m.curve_params,plan.model.curve_params)
    return plan.model.meshing.mesh_authority,plan.max_node_tag
end

struct _GmshNativeSurfaceBackground{B}
    frame::_GmshNativePlaneFrame
    background::B
end
struct _GmshNativeSurfaceLift
    frame::_GmshNativePlaneFrame
end
(lift::_GmshNativeSurfaceLift)(u,v)=_surface_metric_lift(lift.frame,u,v)
function (background::_GmshNativeSurfaceBackground)(u,v)
    return background.background(_surface_metric_lift(background.frame,u,v)...)
end

function _gmsh_native_surface_vertex_sizes(T,inside,source,m,embedded_points)
    sizes=fill(-1.0,T.nreal)
    for triangle in eachindex(T.alive)
        inside[triangle] || continue
        for (first,last) in ((1,2),(1,3),(2,3))
            a=Mesh2D._vert(T,triangle,first);b=Mesh2D._vert(T,triangle,last)
            difference=Mesh2D._gmsh_sub3(source.world[a],source.world[b])
            distance=Mesh2D._gmsh_norm3(difference)
            (sizes[a]<0 || sizes[a]>distance) && (sizes[a]=distance)
            (sizes[b]<0 || sizes[b]>distance) && (sizes[b]=distance)
        end
    end
    # buildMeshGenerationDataStructures caps embedded GVertex sizes before
    # applying the boundary-edge floor. The row comes from its actual Point
    # producer identity, rather than a coordinate lookup into the surface.
    for point in embedded_points
        row=source.point_positions[point]
        sizes[row]=min(sizes[row],_model_vertex_lc(m,point))
    end
    for (a,b) in source.boundaries
        half=0.5*Mesh2D._gmsh_norm3(Mesh2D._gmsh_sub3(source.world[a],source.world[b]))
        sizes[a]<half && (sizes[a]=half)
        sizes[b]<half && (sizes[b]=half)
    end
    return sizes
end

function _gmsh_native_surface_recombine_births(product,mesh,protected,caller)
    isempty(product.element_tags) && return _ModelSurfaceProduction(mesh,
        product.node_tags,product.node_owners,product.node_parameters,
        product.retained_nodes,product.node_high_water)
    source=product.mesh
    length(product.element_tags)==ntris(source) || throw(ArgumentError(
        "$caller: initial Surface triangle births do not align"))
    previous=Dict{NTuple{3,Int32},UInt64}()
    incidence=Dict{NTuple{2,Int32},Int}()
    for column in axes(source.tris,2)
        nodes=source.tris[:,column]
        key=Tuple(sort(nodes))
        haskey(previous,key) && throw(ArgumentError(
            "$caller: initial Surface triangles repeat a primary identity"))
        previous[key]=product.element_tags[column]
        for (a,b) in ((1,2),(2,3),(3,1))
            edge=minmax(nodes[a],nodes[b])
            incidence[edge]=get(incidence,edge,0)+1
        end
    end
    candidates=count(pair->last(pair)==2 &&
        (protected===nothing || !(first(pair) in protected)),incidence)
    quadrangles=sum(size(block.nodes,2) for block in mesh.blocks if block.msh==3;init=0)
    high_water=product.element_high_water
    limit=typemax(UInt64)
    reservation=UInt64(candidates)+UInt64(quadrangles)
    (high_water>limit || reservation>limit-high_water) && throw(ArgumentError(
        "$caller: recombined Surface element births exhaust UInt64 labels"))
    # RecombineTriangle constructs one temporary MQuadrangle per eligible pair
    # before admission. Accepted quads are constructed after that pair pass.
    high_water+=UInt64(candidates)
    births=UInt64[]
    for block in mesh.blocks,column in axes(block.nodes,2)
        if block.msh==2
            key=Tuple(sort(block.nodes[:,column]))
            tag=get(previous,key,UInt64(0))
            tag>0 || throw(ArgumentError(
                "$caller: recombination lost a retained triangle birth"))
            push!(births,tag)
        elseif block.msh==3
            high_water+=UInt64(1)
            push!(births,high_water)
        else
            throw(ArgumentError("$caller: recombination emitted an unsupported cell family"))
        end
    end
    return _ModelSurfaceProduction(mesh,product.node_tags,product.node_owners,
        product.node_parameters,product.retained_nodes,product.node_high_water,
        births,high_water)
end

function _gmsh_native_initial_surface(m::GeoModel,tag::Int,options,input,
                                       high_water,caller;algorithm=3,size_field=nothing,
                                       element_high_water=UInt64(0),
                                       forced=nothing,
                                       _periodic_positions=nothing)
    _surface_type(m,tag)===:plane || throw(ArgumentError(
        "$caller: native initial Surface[$tag] requires its actual Plane parametrization"))
    authority,limit=_gmsh_native_surface_lower_input(m,options,input,high_water,caller)
    fullquad=_model_surface_recombined(m,tag) && m.meshing.recombine_algo in (2,3)
    quad_input=fullquad ? _gmsh_fullquad_input(m,tag,authority,
        element_high_water,options,caller) : nothing
    if quad_input!==nothing
        authority=quad_input.authority
        element_high_water=quad_input.element_high_water
        quad_input.admitted && (options=_gmsh_fullquad_size_options(options))
    end
    frame=_gmsh_native_surface_frame(m,tag,options,caller)
    # Forced periodic subdivisions mint fresh node labels above the consumed
    # lower-product high water; the raised limit feeds later refinements.
    next_tag=Ref{UInt64}(limit)
    source=_gmsh_native_surface_source(m,tag,authority,frame,caller;
        ctx_lc=options.ctx_lc,forced,next_tag,
        _periodic_positions)
    limit=next_tag[]
    T=_gmsh_native_surface_initial(source,tag;cleanup_faces=algorithm==3)
    curve_chains=_gmsh_native_surface_curve_chains(source)
    inside=classify_interior(T)
    numbers=source.tags
    trace=nothing
    if algorithm in (5,6) && length(source.tags)!=3
        embedded_points,_=_model_surface_embedding_tags(m,tag,caller)
        sizes=_gmsh_native_surface_vertex_sizes(T,inside,source,m,embedded_points)
        background=_GmshNativeSurfaceBackground(frame,
            _ModelAlgorithmBackgroundSize(m,2,tag,options,size_field,caller))
        background_sizes=copy(sizes)
        model_vertices=Int32[]
        for point in embedded_points
            row=source.point_positions[point]
            background_sizes[row]=background.background(source.world[row]...)
            push!(model_vertices,Int32(row))
        end
        extend=_resolved_extend_from_boundary(m,2,tag)
        # Tensor background arm: when the resolved size selection is a metric
        # field, hand the trace a per-cell evaluator. Each cavity in-circle
        # test and circumcenter then uses the tensor pulled back at the cell's
        # lifted centroid; the queue's scalar criterion keeps the same
        # directional-size policy as `_metric_surface_algorithm_field`.
        lift=_GmshNativeSurfaceLift(frame)
        trace_metric=let raw=options.inner_field===nothing ? size_field :
                            options.inner_field
            if raw isa AbstractAnisoField
                lc_min=isnan(options.lc_min) ? 0.0 : options.lc_min
                lc_max=isnan(options.lc_max) ? GMSH_MAX_SIZE : options.lc_max
                lo=min(lc_min,lc_max)
                hi=max(lo,min(lc_max,options.ctx_lc))
                wrapped=BoundedSize(raw;size_min=lo,size_max=hi,
                                    factor=options.lc_factor)
                f=frame.first;s=frame.second;entity=(2,tag)
                (u,v)->begin
                    w=lift(u,v)
                    mm=metric_at(wrapped,w[1],w[2],w[3],entity)
                    (mm.m11*f[1]*f[1]+mm.m22*f[2]*f[2]+mm.m33*f[3]*f[3]+
                     2*(mm.m12*f[1]*f[2]+mm.m13*f[1]*f[3]+mm.m23*f[2]*f[3]),
                     mm.m11*f[1]*s[1]+mm.m22*f[2]*s[2]+mm.m33*f[3]*f[3]+
                     mm.m12*(f[1]*s[2]+f[2]*s[1])+mm.m13*(f[1]*s[3]+f[3]*s[1])+
                     mm.m23*(f[2]*s[3]+f[3]*s[2]),
                     mm.m11*s[1]*s[1]+mm.m22*s[2]*s[2]+mm.m33*s[3]*s[3]+
                     2*(mm.m12*s[1]*s[2]+mm.m13*s[1]*s[3]+mm.m23*s[2]*s[3]))
                end
            else
                (Mesh2D._gmsh_dot3(frame.first,frame.first),
                 Mesh2D._gmsh_dot3(frame.first,frame.second),
                 Mesh2D._gmsh_dot3(frame.second,frame.second))
            end
        end
        trace=Mesh2D._gmsh_surface_trace(T,inside,source.world,lift,
            element_high_water,sizes,background_sizes,extend>0 || extend==-2;
            temporary_cells=length(source.tags)==3 ? 0 : count(inside),
            metric=trace_metric)
        insertion=refine!(T;min_angle_deg=0,gmsh_insertion=true,
            gmsh_algorithm=algorithm==5 ? :delaunay : :frontal,
            gmsh_extend_boundary=extend>0 || extend==-2,
            vertex_sizes=sizes,vertex_bgm=background_sizes,bgm_size=background,
            model_vertices=model_vertices,
            gmsh_vertex_numbers=source.tags,gmsh_vertex_high_water=limit,
            return_vertex_numbers=true,_gmsh_surface_trace=trace)
        inside=insertion.interior;numbers=insertion.vertex_numbers
        limit=insertion.vertex_high_water
    end
    product=_gmsh_native_surface_export(T,inside,numbers,source,frame,tag,limit;
        projected=true,trace)
    if algorithm==3 || length(source.tags)==3
        # meshGenerator first constructs the recovered MTriangles, deletes
        # that temporary GFace mesh, then BDS2GMSH constructs the public cells.
        # Its three-vertex shortcut constructs only the public triangle.
        cells=ntris(product.mesh)
        temporary=length(source.tags)==3 ? 0 : cells
        reservation=UInt64(temporary)+UInt64(cells)
        reservation<=typemax(UInt64)-element_high_water || throw(ArgumentError(
            "$caller: initial Surface element births exhaust UInt64 labels"))
        birth_high_water=element_high_water+UInt64(temporary)
        births=UInt64[birth_high_water+UInt64(column) for column in 1:cells]
        product=_ModelSurfaceProduction(product.mesh,product.node_tags,
            product.node_owners,product.node_parameters,product.retained_nodes,
            product.node_high_water,births,birth_high_water+UInt64(cells),product.birth_rows)
    end
    mesh=product.mesh
    if length(source.tags)==3 && ntris(mesh)==1
        nodes=Tuple(mesh.tris[:,1])
        first_block_tags=product.node_tags
        first_slot=argmin(ntuple(slot->first_block_tags[nodes[slot]],3))
        for slot in 1:3
            mesh.tris[slot,1]=nodes[mod1(first_slot+slot-1,3)]
        end
    end
    physical=Matrix{Float64}(undef,3,nnodes(mesh))
    for row in axes(physical,2)
        original=get(source.positions,product.node_tags[row],0)
        physical[:,row].=original>0 ? source.world[original] :
            _surface_metric_lift(frame,mesh.coords[1,row],mesh.coords[2,row])
    end
    if algorithm in (5,6)
        product,physical=_gmsh_native_surface_smooth(m,tag,product,physical,frame,caller)
        mesh=product.mesh
    end
    if _model_surface_recombined(m,tag)
        recombine_algorithm=m.meshing.recombine_algo
        recombine_algorithm in (0,1,2,3,4) || throw(ArgumentError(
            "$caller: native Surface[$tag] experimental recombination is not implemented"))
        protected=_embedded_protected_edges(mesh,source.x,source.y,source.internals)
        if recombine_algorithm==4
            product,physical=_gmsh_native_surface_bipartite(m,tag,product,physical,caller)
        elseif quad_input===nothing
            product=_gmsh_native_surface_recombine(m,tag,product,physical,protected,frame,caller)
        elseif quad_input.admitted
            product,physical=_gmsh_native_surface_fullquad(m,tag,product,physical,protected,
                frame,quad_input,caller;reverse_initial=trace===nothing)
        end
        mesh=product.mesh
    end
    # GFace's per-entity smoothing attribute applies to transfinite meshing.
    # The ordinary native producer uses CTX::mesh.nbSmoothing above.
    # Initial PolyMesh storage is CCW in UV; the native recovered triangles
    # are CW. OCC subsequently orients these cells against its CAD normal in
    # world coordinates, while the built-in mean Plane basis is opposite its
    # declared normal.
    trace===nothing && !(quad_input!==nothing && quad_input.admitted) &&
        (mesh=_reversed_surface_mesh(mesh,caller,tag))
    constructor_mesh=if any(pair->first(pair)[1]==2 &&
            abs(Int(last(pair).master_entity))==tag,m.periodic)
        # copyMesh runs before orientMeshGFace. Keep this operation's actual
        # world rows and cell constructor winding for newly meshed masters.
        mesh isa Mesh ? Mesh(copy(physical);tris=copy(mesh.tris)) :
            MixedMesh(copy(physical),ElementBlock[ElementBlock(block.msh,
                copy(block.nodes),copy(block.tags)) for block in mesh.blocks])
    else
        nothing
    end
    mesh=_gmsh_native_surface_canonical_orientation(m,tag,mesh,physical,caller)
    triangle_count=mesh isa Mesh ? ntris(mesh) :
        sum(size(block.nodes,2) for block in mesh.blocks if block.msh==2;init=0)
    if length(source.tags)==3 && triangle_count==1
        # meshGenerator's three-vertex shortcut iterates its MVertex UID set.
        # Surface orientation can exchange the last two slots, while the
        # smallest born UID remains its constructor's first vertex.
        triangles=mesh isa Mesh ? mesh.tris : only(block.nodes for block in mesh.blocks if block.msh==2)
        wound_nodes=Tuple(triangles[:,1])
        triangle_labels=product.node_tags
        first_slot=argmin(ntuple(slot->triangle_labels[wound_nodes[slot]],3))
        for slot in 1:3
            triangles[slot,1]=wound_nodes[mod1(first_slot+slot-1,3)]
        end
    end
    get(m.meshing.reverse,(2,tag),false) &&
        (mesh=_reversed_surface_mesh(mesh,caller,tag))
    nnodes(mesh)==length(product.node_tags) || throw(ArgumentError(
        "$caller: Surface[$tag] postprocessing changed the actual producer node set"))
    for row in axes(mesh.coords,2)
        original=get(source.positions,product.node_tags[row],0)
        mesh.coords[:,row].=original>0 ? source.world[original] :
            quad_input!==nothing && quad_input.admitted ? physical[:,row] :
            _surface_metric_lift(frame,mesh.coords[1,row],mesh.coords[2,row])
    end
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid native initial Surface[$tag] — "*join(diagnostic.messages,"; ")))
    return _ModelSurfaceProduction(mesh,product.node_tags,product.node_owners,
        product.node_parameters,product.retained_nodes,product.node_high_water,
        product.element_tags,product.element_high_water,product.birth_rows,constructor_mesh,
        product.source_node_updates,product.native_recombine_events,curve_chains)
end

# A pure geometric certificate, independent of the fitted normal's sign.
# Keep candidate selection local so later readonly findfirst closures do not
# capture reassigned anchor/second bindings (Core.Box).
function _gmsh_native_plane_declared_triple(coordinates,projection)
    isempty(coordinates) && return false
    anchor=first(coordinates)
    for second in coordinates
        second==anchor && continue
        for third in coordinates
            orient2((anchor[projection[1]],anchor[projection[2]]),
                    (second[projection[1]],second[projection[2]]),
                    (third[projection[1]],third[projection[2]]))!=0 && return true
        end
        return false
    end
    return false
end

function _gmsh_native_surface_evaluation_frame(m::GeoModel,tag::Int,caller)
    _model_require_plane_surface(m,tag,caller,"Plane parametrization")
    curves=_gmsh_native_surface_curves(m,tag)
    samples=NTuple{3,Float64}[]
    for curve in curves
        low,high=_model_curve_param_bounds(m,curve,caller)
        push!(samples,_gmsh_native_curve_sample(m,curve,low+0.333*(high-low),caller))
        push!(samples,_gmsh_native_curve_sample(m,curve,low+0.666*(high-low),caller))
    end
    options=_model_default_mesh1d_options(m,caller)
    receipt=get(m.meshing.native_geo_kernel.plane_frames,tag,nothing)
    frame=receipt===nothing ? _gmsh_native_plane_frame(samples,options.geometry_tolerance,options.ctx_lc) :
        _GmshNativePlaneFrame(receipt.origin,receipt.first,receipt.second,receipt.normal)
    frame===nothing && return nothing
    point_tags=unique!(_loop_points(m,first(m.surfaces[tag])))
    append!(point_tags,[point for point in _model_points_of(m,[(2,tag)],caller)
        if !(point in point_tags)])
    coordinates=NTuple{3,Float64}[m.points[point] for point in point_tags]
    normal=Mesh2D._gmsh_normalize3(frame.normal)
    dominant=argmax(abs.(frame.normal))
    projection=dominant==3 ? (1,2) : dominant==2 ? (1,3) : (2,3)
    # A curved wire can bound a valid Plane with collinear declared vertices.
    # Keep the already-computed native fit and basis; only its geometric
    # triple certificate needs the same on-curve samples used by the fit.
    if !_gmsh_native_plane_declared_triple(coordinates,projection)
        for (index,coordinate) in enumerate(samples)
            push!(coordinates,coordinate)
            push!(point_tags,m.curves[curves[(index+1)÷2]][1])
        end
    end
    anchor=first(coordinates)
    second_index=findfirst(point->point!=anchor,coordinates)
    second_index===nothing && return nothing
    second=coordinates[second_index]
    third_index=findfirst(coordinates) do candidate
        orient2((anchor[projection[1]],anchor[projection[2]]),
                (second[projection[1]],second[projection[2]]),
                (candidate[projection[1]],candidate[projection[2]]))!=0
    end
    third_index===nothing && return nothing
    result=(point_tags=point_tags,coordinates=coordinates,anchor=anchor,
        second=second,third=coordinates[third_index],projection=projection,
        properties=(normal...,Mesh2D._gmsh_dot3(normal,frame.origin)),
        normal=normal,origin=frame.origin,first_direction=frame.first,
        second_direction=frame.second,native_frame=frame)
    return receipt===nothing ? result : merge(result,(native_bounds=receipt.bounds,))
end

function _gmsh_native_plane_birth_frame(m::GeoModel,tag,tolerance,context_lc)
    caller="native GEO Plane constructor"
    samples=NTuple{3,Float64}[]
    curves=_gmsh_native_surface_curves(m,tag)
    for curve in curves
        low,high=_model_curve_param_bounds(m,curve,caller)
        push!(samples,_gmsh_native_curve_sample(m,curve,low+0.333*(high-low),caller))
        push!(samples,_gmsh_native_curve_sample(m,curve,low+0.666*(high-low),caller))
    end
    frame=_gmsh_native_plane_frame(samples,tolerance,context_lc)
    frame===nothing || return frame
    vertices=sort!(unique!(Int[point for curve in curves for point in m.curves[curve]]))
    points=copy(samples);append!(points,[m.points[point] for point in vertices])
    return _gmsh_native_plane_small_svd_frame(points)
end

function _gmsh_native_plane_boundary_bounds(m::GeoModel,tag,frame)
    caller="native GEO Plane parameter bounds"
    lower=(Inf,Inf);upper=(-Inf,-Inf)
    for curve in _gmsh_native_surface_curves(m,tag)
        low,high=_model_curve_param_bounds(m,curve,caller)
        step=(high-low)/4.0
        for position in 0:4
            point=_gmsh_native_curve_sample(m,curve,low+step*position,caller)
            uv=_gmsh_native_surface_parameters(frame,point)
            lower=(uv[1]<lower[1] ? uv[1] : lower[1],
                   uv[2]<lower[2] ? uv[2] : lower[2])
            upper=(uv[1]>upper[1] ? uv[1] : upper[1],
                   uv[2]>upper[2] ? uv[2] : upper[2])
        end
    end
    all(isfinite,(lower...,upper...)) || throw(ArgumentError(
        "$caller: Surface[$tag] has no finite boundary samples"))
    return lower,upper
end

# Gmsh copyMesh consumes the actual master and its Point/Curve counterpart
# vectors. Coordinates verify a bound identity; they never select a vertex.
# Copying and link publication consume the same actual master row view.
# Tensor producers expose their lattice rows; other producers retain the
# native constructor or storage view already bound to their UID records.
function _model_periodic_surface_source(m,constraint,input,caller)
    master=Int(constraint.master_entity)
    actual=_model_mesh_product(input,2,master;caller,producer_order=true)
    source=get(input.surface_constructor_inputs,(2,Int32(master)),nothing)
    if source===nothing
        source=actual
    else
        actual!==nothing || throw(ArgumentError(
            "$caller: periodic constructor lacks its actual Surface[$master] product"))
        _model_surface_constructor_binding(actual,source,caller)
    end
    return source
end

function _model_periodic_surface_production(m::GeoModel,
        constraint::ModelPeriodicConstraint,input::_ModelMeshAuthority,
        node_high_water,element_high_water,caller::AbstractString)
    slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
    constraint.dim==2 || throw(ArgumentError(
        "$caller: actual periodic Surface producer requires dimension two"))
    source=_model_periodic_surface_source(m,constraint,input,caller)
    completed_empty=source!==nothing && source.entity_data!==nothing &&
        isempty(source.blocks) && (2,Int32(master)) in input.done &&
        haskey(input.products,(2,Int32(master))) &&
        isempty(input.products[(2,Int32(master))].blocks) &&
        all(owner->owner[1]<2,source.entity_data.node_entities)
    source!==nothing && source.entity_data!==nothing &&
        (_surface_cell_count(source)>0 || completed_empty) || throw(ArgumentError(
        "$caller: periodic Surface[$slave] lacks the actual Surface[$master] product"))
    maximum_tag=maximum(input.mesh.entity_data.external_node_tags;init=UInt64(0))
    limit=node_high_water===nothing ? maximum_tag : UInt64(node_high_water)
    limit>=maximum_tag || throw(ArgumentError(
        "$caller: periodic high water omits actual source nodes"))
    cell_limit=UInt64(element_high_water)
    maximum((tag for block in input.mesh.entity_data.external_element_tags
        for tag in block);init=UInt64(0))<=cell_limit || throw(ArgumentError(
        "$caller: periodic high water omits actual source cells"))
    point_map,curve_pairs,reversed=_model_periodic_surface_edge_map(
        m,slave,master,constraint.affine,constraint.atol,caller)
    counterparts=Dict{UInt64,Tuple{MixedMesh,Int}}()
    function bind!(master_product,master_row,slave_product,slave_row)
        master_tag=master_product.entity_data.external_node_tags[master_row]
        previous=get(counterparts,master_tag,nothing)
        if previous!==nothing
            previous_mesh,previous_row=previous
            previous_mesh.entity_data.external_node_tags[previous_row]==
                slave_product.entity_data.external_node_tags[slave_row] ||
                throw(ArgumentError("$caller: periodic source node $master_tag has conflicting counterparts"))
        end
        counterparts[master_tag]=(slave_product,Int(slave_row))
        return nothing
    end
    # Native principal vertices use the first actual Point-owned vertex.
    # A Point cell and an unrelated coincident vertex are separate objects.
    for master_point in sort!(collect(keys(point_map)))
        slave_point=point_map[master_point]
        master_product=_model_mesh_product(input,0,master_point;caller)
        slave_product=_model_mesh_product(input,0,slave_point;caller)
        master_product!==nothing && slave_product!==nothing || throw(ArgumentError(
            "$caller: periodic Point counterpart lacks actual source products"))
        master_rows=findall(==((0,Int32(master_point))),master_product.entity_data.node_entities)
        slave_rows=findall(==((0,Int32(slave_point))),slave_product.entity_data.node_entities)
        (isempty(master_rows) || isempty(slave_rows)) && throw(ArgumentError(
            "$caller: periodic Point counterpart lacks an owned actual vertex"))
        bind!(master_product,first(master_rows),slave_product,first(slave_rows))
    end
    # Private ranks retain native GEdge::mesh_vertices order independently
    # of public labels and Line-cell orientation. Endpoints are mapped above.
    for (index,(slave_curve,master_curve)) in enumerate(curve_pairs)
        master_product=_model_mesh_product(input,1,master_curve;caller)
        slave_product=_model_mesh_product(input,1,slave_curve;caller)
        master_product!==nothing && slave_product!==nothing || throw(ArgumentError(
            "$caller: periodic Curve counterpart lacks actual source products"))
        master_rows=findall(==((1,Int32(master_curve))),master_product.entity_data.node_entities)
        slave_rows=findall(==((1,Int32(slave_curve))),slave_product.entity_data.node_entities)
        length(master_rows)==length(slave_rows) || throw(ArgumentError(
            "$caller: periodic Curve[$slave_curve] owned vertex count differs from Curve[$master_curve]"))
        reversed[index] && reverse!(master_rows)
        for (master_row,slave_row) in zip(master_rows,slave_rows)
            bind!(master_product,master_row,slave_product,slave_row)
        end
    end
    count=nnodes(source);coords=Matrix{Float64}(undef,3,count)
    tags=Vector{UInt64}(undef,count);owners=Vector{Tuple{Int,Int32}}(undef,count)
    parameters=Vector{Union{Nothing,Vector{Float64}}}(undef,count)
    retained=falses(count);data=source.entity_data
    for row in axes(source.coords,2)
        owner=data.node_entities[row]
        if owner==(2,Int32(master))
            # Literal copyMesh accumulation, then target->parFromPoint and
            # target->point. Raw affine XYZ is not the stored slave vertex.
            point=(source.coords[1,row],source.coords[2,row],source.coords[3,row],1.0)
            transformed=ntuple(3) do axis
                value=0.0
                for component in 1:4
                    value+=constraint.affine[4(axis-1)+component]*point[component]
                end
                value
            end
            parameter=model_parametrization(m,2,slave,collect(transformed))
            coordinate=model_value(m,2,slave,parameter)
            length(parameter)==2 && length(coordinate)==3 &&
                all(isfinite,parameter) && all(isfinite,coordinate) || throw(ArgumentError(
                "$caller: periodic Surface[$slave] produced a non-finite parametrized vertex"))
            limit<typemax(UInt64) || throw(ArgumentError(
                "$caller: periodic node allocation exhausts UInt64 public labels"))
            limit+=UInt64(1);tags[row]=limit;coords[:,row].=coordinate
            owners[row]=(2,Int32(slave));parameters[row]=parameter
        else
            counterpart=get(counterparts,data.external_node_tags[row],nothing)
            counterpart!==nothing || throw(ArgumentError(
                "$caller: actual periodic source node $(data.external_node_tags[row]) has no Point/Curve counterpart"))
            product,position=counterpart;target_data=product.entity_data
            tags[row]=target_data.external_node_tags[position]
            coords[:,row].=product.coords[:,position]
            owners[row]=target_data.node_entities[position]
            parameter=target_data.node_parametric[position]
            parameters[row]=parameter===nothing ? nothing : copy(parameter)
            retained[row]=true
        end
    end
    allunique(tags) || throw(ArgumentError(
        "$caller: periodic emitted rows repeat an actual node identity"))
    blocks=ElementBlock[];cell_tags=UInt64[]
    # Native copyMesh creates triangles, then quadrangles, with source winding.
    for (msh,nodes) in sort!(_surface_dim2_blocks(source);by=first)
        msh in (2,3) || throw(ArgumentError(
            "$caller: actual periodic Surface copy requires linear Triangles or Quadrangles"))
        for column in axes(nodes,2)
            cell_limit<typemax(UInt64) || throw(ArgumentError(
                "$caller: periodic cell allocation exhausts UInt64 public labels"))
            cell_limit+=UInt64(1);push!(cell_tags,cell_limit)
        end
        push!(blocks,ElementBlock(msh,copy(nodes)))
    end
    mesh=length(blocks)==1 && only(blocks).msh==2 ?
        Mesh(coords;tris=only(blocks).nodes) : MixedMesh(coords,blocks)
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: actual periodic Surface[$slave] mesh is invalid — " *
        join(diagnostic.messages,"; ")))
    return _ModelSurfaceProduction(mesh,tags,owners,parameters,retained,limit,
        cell_tags,cell_limit)
end
