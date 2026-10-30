# `BoundaryLayer Field =` topological surface pass — upstream's
# `Mesh2DWithBoundaryLayers` (`Generator.cpp`) drives `buildAdditionalPoints2D`
# (`boundaryLayersData.cpp`), which grows vertex columns off each face's
# wall-listed boundary edges before `modifyInitialMeshForBoundaryLayers`
# (`meshGFace.cpp`) emits the strip and re-fills the shrunken face. This file
# transports the resolved field parameters (thickness/ratio/hwall/fan targets)
# into the same bounded construction: `_bl_2d_offsets` emits the certified
# strip and the ordinary constrained-Delaunay fill covers the remainder.
#
# Bounded scope, matching the upstream option set the field record carries:
# - Plane surfaces only (`_model_surface_plane` admission); upstream's curved
#   `parFromPoint` column construction is a documented gap.
# - `Size`/`Ratio`/`Thickness`/`Quads`/`FanPointsList`/`ExcludedSurfacesList`
#   and the `BetaLaw`/`NbLayers`/`Beta` schedule apply; `SizesList`,
#   `FanPointsSizesList`, `IntersectMetrics` and `AnisoMax` stay loud
#   rejections rather than silently dropping parameters.
# - `PointsList`/`NodesList` end-node grading (`addBoundaryLayerPoints` in
#   `meshGEdge.cpp`) is not ported — the field's curve targets drive the
#   strip here; endpoint tags only mark fan membership like upstream.
# - Transfinite, extruded, periodic-slave and parametric (OCC) surfaces keep
#   their own producers; a `BoundaryLayer` field never alters them, matching
#   upstream's dispatch order (transfinite/extruded paths run before the
#   boundary-layer pass and return first).

# Per-selected-field normalized record (`_ModelBoundaryLayerField` is
# declared in `Model.jl` next to `ModelMeshingAttributes`, which stores the
# resolved map). Entity tags are literal model tags; the per-face curve
# intersection upstream computes in `setupFor2d` happens at mesh time in
# `_model_surface_bl_specs` below.
function _ModelBoundaryLayerField(tag::Integer,config,caller::AbstractString)
    isempty(config.fan_sizes) || length(config.fan_sizes)==
        length(config.fan_points) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$tag].FanPointsSizesList must be " *
        "empty or match FanPointsList length"))
    return _ModelBoundaryLayerField(
        Int(tag),Int.(collect(config.curves)),Int.(collect(config.points)),
        Int.(collect(config.fan_points)),Int.(collect(config.fan_sizes)),
        Int.(collect(config.excluded_surfaces)),
        Float64(config.size),Float64.(collect(config.sizes)),
        Float64(config.ratio),Float64(config.size_far),
        Float64(config.thickness),!iszero(config.quads),
        Int(config.fan_elements),!iszero(config.beta_law),
        Float64(config.beta),Int(config.nb_layers),
        !iszero(config.intersect_metrics),Float64(config.aniso_max),
        Set{Int}(Int.(config.curves)),Set{Int}(Int.(config.fan_points)),
        Set{Int}(Int.(config.excluded_surfaces)))
end

# Entity retag for `_model_identity_*` moves — entity-carrying option lists
# remap like every other attribute; the numeric schedule crosses unchanged.
function _model_identity_bl_field(bl::_ModelBoundaryLayerField,
                                  dimension::Int,old_tag::Int,new_tag::Int)
    remap(list,dim)=Int[v==old_tag && dim==dimension ? new_tag : v for v in list]
    curves=remap(bl.curves,1);fans=remap(bl.fan_points,0)
    excluded=remap(bl.excluded_surfaces,2);points=remap(bl.points,0)
    return _ModelBoundaryLayerField(
        bl.tag,curves,points,fans,copy(bl.fan_sizes),excluded,
        bl.hwall,copy(bl.sizes),bl.ratio,bl.hfar,bl.thickness,bl.quads,
        bl.fan_elements,bl.beta_law,bl.beta,bl.nb_layers,
        bl.intersect_metrics,bl.aniso_max,
        Set{Int}(curves),Set{Int}(fans),Set{Int}(excluded))
end

# Geometric row schedule (`L = hWall; while L <= thickness { push L;
# L += hWall·ratio^ith }` — boundaryLayersData.cpp's non-beta branch) and its
# Beta-law alternative (`t[i] = 1+β(1−p)/(1+p)`, `L = hWall·t[i]/t[0]`).
function _bl_field_offsets(bl::_ModelBoundaryLayerField,
                           caller::AbstractString)
    if bl.beta_law
        beta=bl.beta
        (isfinite(beta) && beta>1) || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(bl.tag)].Beta must exceed 1"))
        nl=bl.nb_layers
        nl>=1 || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(bl.tag)].NbLayers must be positive"))
        zlog=log((1+beta)/(beta-1))
        t=Vector{Float64}(undef,nl)
        for i in 1:nl
            p=exp(zlog*(1-i/nl));t[i]=1+beta*(1-p)/(1+p)
        end
        return Float64[bl.hwall*t[i]/t[1] for i in 1:nl]
    end
    hw=bl.hwall
    (isfinite(hw) && hw>0) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)].Size must be positive"))
    (isfinite(bl.ratio) && bl.ratio>1) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)].Ratio must exceed 1"))
    (isfinite(bl.thickness) && bl.thickness>0) || throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)].Thickness must be positive"))
    offsets=Float64[]
    L=hw;k=0
    while L<=bl.thickness
        push!(offsets,L)
        k+=1
        L+=hw*bl.ratio^k
        length(offsets)>100_000 && throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(bl.tag)] schedule does not " *
            "reach Thickness=$(bl.thickness)"))
    end
    isempty(offsets) && throw(ArgumentError(
        "$caller: BoundaryLayer Field[$(bl.tag)] Thickness=" *
        "$(bl.thickness) admits no layer row"))
    return offsets
end

# Fields whose `CurvesList` intersects this face's boundary and whose
# `ExcludedSurfacesList` does not suppress it — `setupFor2d(face)`'s
# `_curveTags ∩ face.edges` test. A field listing only `PointsList` end
# nodes contributes nothing topologically here. Like upstream, a wall
# curve shared by more than one un-excluded surface rejects — only 2-D
# boundary layers are supported.
function _model_surface_bl_specs(m::GeoModel,t::Int)
    fields=m.meshing.boundary_layer_fields
    isempty(fields) && return nothing
    boundary=Set{Int}()
    for loop in m.surfaces[t],signed in m.loops[loop]
        push!(boundary,abs(signed))
    end
    specs=_ModelBoundaryLayerField[]
    for tag in sort!(collect(keys(fields)))
        spec=fields[tag]
        t in spec.excluded_set && continue
        shared=intersect(boundary,spec.curve_set)
        isempty(shared) && continue
        # `setupFor2d` counts the wall curve's adjacent faces that are
        # neither extruded nor excluded; more than one means the layer
        # would grow inside a volume — upstream errors.
        for curve in shared
            numf=0
            for (stag,loops) in m.surfaces
                stag in spec.excluded_set && continue
                any(loop->any(signed->abs(signed)==curve,
                              m.loops[loop]),loops) || continue
                numf+=1
            end
            numf>1 && throw(ArgumentError(
                "Only 2D Boundary Layers are supported (curve $curve " *
                "is adjacent to $numf surfaces)"))
        end
        push!(specs,spec)
    end
    return isempty(specs) ? nothing : specs
end

# `_surface_pslg` twin that additionally records the ordered wire vertex
# lists and each boundary segment's source curve — the wall-run extraction
# below needs the curve ownership of every boundary segment.
function _bl_surface_pslg(m::GeoModel,t::Int,forced,caller::AbstractString;
                          param_sizes::Dict{Tuple{Int,Float64},Float64}=
                              Dict{Tuple{Int,Float64},Float64}(),
                          plane=_model_surface_plane(m,t,caller;
                                                     allow_ruled=true),
                          _periodic_positions=nothing)
    xs=Float64[];ys=Float64[];mesh_sizes=Float64[]
    index=Dict{Int,Int}()
    canonical=Dict{Tuple{Float64,Float64},NTuple{3,Float64}}()
    wires=Vector{Int}[]
    wire_curves=Vector{Int}[]
    vert_curve=Dict{Int,Int}()
    for loop_id in m.surfaces[t]
        _verify_loop_closed(m,loop_id,caller,"Surface[$t]")
        loop_idx=Int[]
        loop_curve=Int[]
        for signed in m.loops[loop_id]
            curve=abs(signed)
            a,b=m.curves[curve]
            t0,t1=_model_curve_param_bounds(m,curve,caller)
            for parameter in _surface_curve_parameters(m,forced,curve,signed,
                                                       caller)
                vertex=if parameter==t0
                    _add_surface_point!(xs,ys,mesh_sizes,index,canonical,
                                        m,a,caller,plane)
                elseif parameter==t1
                    _add_surface_point!(xs,ys,mesh_sizes,index,canonical,
                                        m,b,caller,plane)
                else
                    point=_periodic_forced_curve_point(
                        m,curve,parameter,_periodic_positions,caller)
                    scale=max(1.0,hypot(point...))
                    abs(_plane_offset(plane,point))<=1e-12*scale ||
                        throw(ArgumentError(
                            "$caller: Curve[$curve] subdivision is not " *
                            "coplanar with Surface[$t]"))
                    _add_surface_curve_point!(
                        xs,ys,mesh_sizes,index,canonical,m,point,
                        get(param_sizes,(curve,parameter),
                            _surface_curve_mesh_size(
                                m,curve,parameter,caller)),
                        curve,caller,plane)
                end
                if !isempty(loop_idx) && loop_idx[end]==vertex
                    # Corner shared by two boundary curves: the slot
                    # owns the outgoing segment, so the continuing
                    # curve's tag wins.
                    loop_curve[end]=curve
                else
                    push!(loop_idx,vertex)
                    push!(loop_curve,curve)
                end
                parameter==t0 || parameter==t1 ||
                    haskey(vert_curve,vertex) || (vert_curve[vertex]=curve)
            end
        end
        while length(loop_idx)>1 && loop_idx[end]==loop_idx[1]
            pop!(loop_idx);pop!(loop_curve)
        end
        nloop=length(loop_idx)
        nloop>=3 || throw(ArgumentError(
            "$caller: Loop[$loop_id] needs at least three points"))
        push!(wires,loop_idx)
        push!(wire_curves,loop_curve)
    end
    embedded=get(m.embeds,(2,t),NTuple{2,Int}[])
    internal=Tuple{Int,Int}[]
    for (edim,etag) in embedded
        if edim==0
            _add_surface_point!(xs,ys,mesh_sizes,index,canonical,m,etag,
                                caller,plane)
        elseif edim!=1
            throw(ArgumentError(
                "$caller: unsupported embedding dimension $edim"))
        end
    end
    for (edim,etag) in embedded
        edim==1 || continue
        haskey(m.curves,etag) || throw(ArgumentError(
            "$caller: unknown embedded Curve[$etag]"))
        a,b=m.curves[etag]
        t0,t1=_model_curve_param_bounds(m,etag,caller)
        parameters=etag in m.meshing.degenerated ? nothing :
                   get(forced,etag,nothing)
        curve_nodes=Int[]
        if parameters===nothing
            push!(curve_nodes,_add_surface_point!(
                xs,ys,mesh_sizes,index,canonical,m,a,caller,plane))
            push!(curve_nodes,_add_surface_point!(
                xs,ys,mesh_sizes,index,canonical,m,b,caller,plane))
        else
            for parameter in parameters
                vertex=if parameter==t0
                    _add_surface_point!(xs,ys,mesh_sizes,index,canonical,m,a,
                                        caller,plane)
                elseif parameter==t1
                    _add_surface_point!(xs,ys,mesh_sizes,index,canonical,m,b,
                                        caller,plane)
                else
                    point=_periodic_forced_curve_point(
                        m,etag,parameter,_periodic_positions,caller)
                    scale=max(1.0,hypot(point...))
                    abs(_plane_offset(plane,point))<=1e-12*scale ||
                        throw(ArgumentError(
                            "$caller: embedded Curve[$etag] subdivision is " *
                            "not coplanar with Surface[$t]"))
                    _add_surface_curve_point!(
                        xs,ys,mesh_sizes,index,canonical,m,point,
                        get(param_sizes,(etag,parameter),
                            _surface_curve_mesh_size(
                                m,etag,parameter,caller)),
                        etag,caller,plane)
                end
                push!(curve_nodes,vertex)
            end
        end
        length(curve_nodes)>=2 || throw(ArgumentError(
            "$caller: embedded Curve[$etag] needs two distinct endpoints"))
        for segment_index in 1:(length(curve_nodes)-1)
            first_node=curve_nodes[segment_index]
            second_node=curve_nodes[segment_index+1]
            (first_node==second_node ||
             (xs[first_node]==xs[second_node] &&
              ys[first_node]==ys[second_node])) && throw(ArgumentError(
                "$caller: embedded Curve[$etag] has coincident subdivision nodes"))
            push!(internal,(first_node,second_node))
        end
        for node in curve_nodes
            haskey(vert_curve,node) || (vert_curve[node]=etag)
        end
    end
    return xs,ys,mesh_sizes,wires,wire_curves,embedded,internal,index,
           canonical,vert_curve
end

# Even-odd containment over the face's wires — winding-insensitive, so hole
# wires and reversed outer wires classify identically (`inside` = contained
# by an odd number of wires).
function _bl_ray_cast(px::Float64,py::Float64,
                      xs::Vector{Float64},ys::Vector{Float64},
                      wire::Vector{Int})
    n=length(wire);inside=false
    for k in 1:n
        i=wire[k];j=wire[mod1(k+1,n)]
        y1=ys[i];y2=ys[j]
        (y1>py)!=(y2>py) || continue
        xint=xs[i]+(py-y1)*(xs[j]-xs[i])/(y2-y1)
        xint>px && (inside=!inside)
    end
    return inside
end

function _bl_point_in_domain(px::Float64,py::Float64,
                             xs::Vector{Float64},ys::Vector{Float64},
                             wires::Vector{Vector{Int}})
    count=0
    for wire in wires
        _bl_ray_cast(px,py,xs,ys,wire) && (count+=1)
    end
    return isodd(count)
end

# Walk a run end's covered span along the wire — upstream's
# `addColumnAtTheEndOfTheBL` makes the end vertex's column the adjacent
# edge's own node chain, so the covered span may cross several boundary
# segments (the boundary write-back BL-grades the adjacent curve near the
# end, exactly like upstream's `addBoundaryLayerPoints`). Bounded: the
# span must stay on `adjacent_curve` (one upstream `GEdge`), must not
# reach a wall segment or another end's claimed position ("boundary
# layers must be separate"), and must be straight within tolerance —
# bending end columns are a documented gap. Returns `(top_vertex,
# dir3d)`: the vertex id the end column's top welds onto (an existing
# wire vertex when the consumed depth lands on it, else an interpolated
# `new_vertex!` node) and the unit 3-D direction of the first covered
# segment — the `vertex_dirs` override for that end.
function _bl_walk_covered!(dsign::Int,p_end::Int,adjacent_curve::Int,
                         offend::Float64,spec,caller::AbstractString,
                         t::Int,covered_positions::Set{Int},
                         claimed::Set{Int},
                         seg_cover::Dict{Int,Vector{Tuple{Float64,Float64}}},
                         snapped_coords::Vector{NTuple{3,Float64}},
                         covered_owner::Dict{NTuple{3,Float64},Int32},
                         seg_spec::Vector{Int},loop::Vector{Int},
                         loop_curve::Vector{Int},xs::Vector{Float64},
                         ys::Vector{Float64},
                         coords3d::Vector{NTuple{3,Float64}},
                         start0::Int,nseg::Int,new_vertex!,vertex_id!,ax)
    pmod(p)=mod1(p,nseg)
    rvert(p)=loop[mod1(start0+p-1,nseg)]
    stol=1e-9*max(1.0,offend)
    function claim_seg!(sidx::Int,covlen::Float64,seglen::Float64)
        # Covered interval on this segment measured from its own start
        # vertex (`loop[sidx]`); a forward walk starts at 0, a backward
        # walk ends at seglen.
        lo,hi=dsign>0 ? (0.0,covlen) : (seglen-covlen,seglen)
        spans=get!(seg_cover,sidx,Tuple{Float64,Float64}[])
        for (olo,ohi) in spans
            (lo<ohi-stol && hi>olo+stol) && throw(ArgumentError(
                "$caller: Surface[$t] boundary-layer run ends consume " *
                "overlapping spans of one boundary segment — upstream " *
                "requires separated layers"))
        end
        push!(spans,(lo,hi))
        return nothing
    end
    pos=pmod(p_end)
    traveled=0.0
    dir3=nothing
    d2x=0.0;d2y=0.0
    while true
        npos=pmod(pos+dsign)
        sidx=dsign>0 ? mod1(start0+pos-1,nseg) :
                      mod1(start0+npos-1,nseg)
        seg_spec[sidx]==0 || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)] layer depth " *
            "$offend on Surface[$t] reaches a wall segment — the " *
            "adjacent non-wall run is too short"))
        ctag=loop_curve[sidx]
        ctag==adjacent_curve || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)] layer depth " *
            "$offend on Surface[$t] spans more than the adjacent " *
            "Curve[$adjacent_curve] — upstream attaches the end column " *
            "to one edge"))
        a=rvert(pos);b=rvert(npos)
        s3x=coords3d[b][1]-coords3d[a][1]
        s3y=coords3d[b][2]-coords3d[a][2]
        s3z=coords3d[b][3]-coords3d[a][3]
        sl3=hypot(s3x,s3y,s3z)
        sl3>0 || throw(ArgumentError(
            "$caller: Surface[$t] boundary has a coincident vertex"))
        s2x=xs[b]-xs[a];s2y=ys[b]-ys[a]
        sl2=hypot(s2x,s2y)
        if dir3===nothing
            dir3=(s3x/sl3,s3y/sl3,s3z/sl3)
            d2x=s2x/sl2;d2y=s2y/sl2
        else
            cross=s2x*d2y-s2y*d2x
            abs(cross)<=1e-9*sl2 || throw(ArgumentError(
                "$caller: BoundaryLayer Field[$(spec.tag)] covered span " *
                "on Surface[$t] is not straight — bending end columns " *
                "are not implemented"))
        end
        dnew=traveled+sl2
        if dnew<offend-stol
            npos in claimed && throw(ArgumentError(
                "$caller: Surface[$t] boundary-layer run ends consume " *
                "overlapping spans — upstream requires separated layers"))
            claim_seg!(sidx,sl2,sl2)
            push!(claimed,npos);push!(covered_positions,npos)
            push!(snapped_coords,coords3d[b])
            covered_owner[coords3d[b]]=Int32(ctag)
            traveled=dnew;pos=npos;continue
        end
        npos in claimed && throw(ArgumentError(
            "$caller: Surface[$t] boundary-layer run ends consume " *
            "overlapping spans — upstream requires separated layers"))
        claim_seg!(sidx,offend-traveled,sl2)
        if abs(dnew-offend)<=stol
            # The consumed depth lands on this wire vertex — the end
            # column's top welds onto it and it stays in the wire.
            push!(claimed,npos)
            push!(snapped_coords,coords3d[b])
            return vertex_id!(b),dir3
        end
        frac=(offend-traveled)/sl2
        tc=(coords3d[a][1]+frac*s3x,
            coords3d[a][2]+frac*s3y,
            coords3d[a][3]+frac*s3z)
        for s in snapped_coords
            (abs(tc[1]-s[1])<=stol && abs(tc[2]-s[2])<=stol &&
             abs(tc[3]-s[3])<=stol) && throw(ArgumentError(
                "$caller: Surface[$t] boundary-layer run ends meet at " *
                "one point — upstream requires separated layers"))
        end
        covered_owner[tc]=Int32(ctag)
        push!(snapped_coords,tc)
        return new_vertex!(tc[ax[1]],tc[ax[2]],tc),dir3
    end
end

# Strip + fill construction for one face. Returns Mesh (all triangles) or
# MixedMesh (triangles + kept quadrangles when the field's `Quads` or the
# face's recombine flag retains them, like upstream's
# `splitElementsInBoundaryLayerIfNeeded` gate).
function _model_bl_surface_mesh(m::GeoModel,t::Int,
                                specs::Vector{_ModelBoundaryLayerField},
                                forced,
                                param_sizes::Dict{Tuple{Int,Float64},Float64},
                                plane,caller::AbstractString;
                                _periodic_positions=nothing)
    for spec in specs
        isempty(spec.sizes) || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].SizesList " *
            "(hwall_n_nodes) is not implemented"))
        isempty(spec.fan_sizes) || throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].FanPointsSizesList " *
            "is not implemented"))
        spec.intersect_metrics && throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].IntersectMetrics " *
            "is not implemented"))
        spec.aniso_max<1e9 && throw(ArgumentError(
            "$caller: BoundaryLayer Field[$(spec.tag)].AnisoMax is not " *
            "implemented"))
    end
    xs,ys,mesh_sizes,wires,wire_curves,embedded,internal,index,canonical,
        vert_curve=
        _bl_surface_pslg(m,t,forced,caller;param_sizes=param_sizes,
                         plane=plane,_periodic_positions=_periodic_positions)
    ax=plane.axes
    coords3d=Vector{NTuple{3,Float64}}(undef,length(xs))
    point_of=Dict{Int,Int}()
    for (pid,v) in index
        point_of[v]=pid
    end
    for i in eachindex(xs)
        key=(xs[i]==0 ? 0.0 : xs[i],ys[i]==0 ? 0.0 : ys[i])
        coordinate=get(canonical,key,nothing)
        coordinate===nothing && throw(ErrorException(
            "$caller: Surface[$t] boundary vertex has no canonical " *
            "coordinate"))
        coords3d[i]=coordinate
    end
    # The surface's true unit normal, signed so that a positive rotation
    # around it matches the left-hand side in the (x,y) projection frame.
    nvec=plane.normal
    nlen=plane.norm
    nsurf=(nvec[1]/nlen,nvec[2]/nlen,nvec[3]/nlen)
    nuv=ax==(1,2) ? (0.0,0.0,1.0) : ax==(1,3) ? (0.0,-1.0,0.0) : (1.0,0.0,0.0)
    flip=nsurf[1]*nuv[1]+nsurf[2]*nuv[2]+nsurf[3]*nuv[3]
    plane_normal=flip>0 ? nsurf : (-nsurf[1],-nsurf[2],-nsurf[3])
    blmod=getfield(parentmodule(@__MODULE__),:BoundaryLayer)
    new_xs=Float64[];new_ys=Float64[]
    new_coords=NTuple{3,Float64}[]
    new_owners=Tuple{Int,Int32}[]
    new_params=Vector{Union{Nothing,Vector{Float64}}}()
    vertex_index=Dict{NTuple{2,Float64},Int}()
    # Nodes a run end column emits on the adjacent non-wall segment —
    # upstream's `addColumnAtTheEndOfTheBL` makes that column the adjacent
    # edge's own mesh chain, so those coordinates are owned by the
    # adjacent curve, not the surface.
    covered_owner=Dict{NTuple{3,Float64},Int32}()
    # Exact wire-vertex (or interpolated top) coordinates that strip
    # column nodes must weld onto bitwise — column nodes computed from
    # `vertex_dirs` may differ by a rounding ulp from the covered wire
    # vertices the boundary write-back produced.
    snapped_coords=NTuple{3,Float64}[]
    function vertex_id!(i::Int)
        key=(xs[i]==0 ? 0.0 : xs[i],ys[i]==0 ? 0.0 : ys[i])
        return get!(vertex_index,key) do
            push!(new_xs,xs[i]);push!(new_ys,ys[i])
            push!(new_coords,coords3d[i])
            pid=get(point_of,i,0)
            if pid!=0
                push!(new_owners,(0,Int32(pid)))
                push!(new_params,nothing)
            else
                curve=get(vert_curve,i,0)
                if curve!=0
                    push!(new_owners,(1,Int32(curve)))
                    push!(new_params,nothing)
                else
                    push!(new_owners,(2,Int32(t)))
                    push!(new_params,[xs[i],ys[i]])
                end
            end
            length(new_coords)
        end
    end
    function new_vertex!(x::Float64,y::Float64,coordinate)
        key=(x==0 ? 0.0 : x,y==0 ? 0.0 : y)
        return get!(vertex_index,key) do
            push!(new_xs,x);push!(new_ys,y);push!(new_coords,coordinate)
            ctag=get(covered_owner,coordinate,Int32(0))
            if ctag!=0
                push!(new_owners,(1,ctag))
                push!(new_params,nothing)
            else
                push!(new_owners,(2,Int32(t)))
                push!(new_params,[x,y])
            end
            length(new_coords)
        end
    end
    strip_coords=NTuple{3,Float64}[]
    strip_tris=NTuple{3,Int}[]
    strip_quads=NTuple{4,Int}[]
    reduced_wires=Vector{Int}[]
    # Map each wire's segment index to its field (0 = ordinary boundary
    # segment). Segments carry the tag of the curve that emitted their first
    # vertex — `_bl_surface_pslg`'s per-run bookkeeping.
    for (w,loop) in enumerate(wires)
        loop_curve=wire_curves[w]
        nseg=length(loop)
        seg_spec=Vector{Int}(undef,nseg)
        for k in 1:nseg
            curve=loop_curve[k]
            owner=0
            for (si,spec) in enumerate(specs)
                if curve in spec.curve_set
                    owner=si
                    break
                end
            end
            seg_spec[k]=owner
        end
        all(==(0),seg_spec) && begin
            push!(reduced_wires,[vertex_id!(v) for v in loop])
            continue
        end
        # Seed every wire vertex so wall-chain and run-end nodes keep
        # their point/curve ownership through the strip pool — upstream
        # reuses the same MVertex objects along the covered boundary.
        for v in loop
            vertex_id!(v)
        end
        if all(!=(0),seg_spec) && !all(==(seg_spec[1]),seg_spec)
            throw(ArgumentError(
                "$caller: Surface[$t] wall segments split across " *
                "boundary-layer fields around a closed wire is not " *
                "implemented"))
        end
        # Rotate to a non-wall segment (a wholly-walled wire stays put and
        # produces the closed-ring run).
        start0=1
        if !all(!=(0),seg_spec)
            start0=findfirst(==(0),seg_spec)::Int+1
            start0>nseg && (start0=1)
        end
        # Maximal cyclic runs of same-field wall segments: (spec_index,
        # first_seg, last_seg, closed) in rotated walk order.
        runs=NTuple{4,Int}[]
        k=1
        while k<=nseg
            sk=mod1(start0+k-1,nseg)
            owner=seg_spec[sk]
            if owner==0
                k+=1
                continue
            end
            first=k
            last=k
            while k+1<=nseg && seg_spec[mod1(start0+k,nseg)]==owner
                k+=1;last=k
            end
            push!(runs,(owner,first,last,
                        first==1 && last==nseg ? 1 : 0))
            k+=1
        end
        # A closed run covers the whole wire; otherwise run segs first..last
        # span wall vertices at rotated positions first..last+1.
        top_of=Dict{Int,Int}()
        # Run-end covered spans: wire positions consumed by the strip side
        # (skipped in the reduced wire), positions already claimed by an
        # end walk (two ends covering one vertex = touching layers), exact
        # snap targets for strip column nodes, and per-end-vertex top
        # vertex overrides.
        covered_positions=Set{Int}()
        claimed=Set{Int}()
        seg_cover=Dict{Int,Vector{Tuple{Float64,Float64}}}()
        top_override=Dict{Int,Int}()
        for (si,first,last,closed_flag) in runs
            spec=specs[si]
            closed_run=!iszero(closed_flag)
            chain_positions=closed_run ? collect(1:nseg) :
                collect(first:last+1)
            chain=Int[loop[mod1(start0+p-1,nseg)] for p in chain_positions]
            nchain=length(chain)
            # Domain-side orientation: probe a point off the first wall
            # segment's left normal; flip the chain when the wall's left is
            # not the face interior.
            a=chain[1];b=chain[2]
            dx=xs[b]-xs[a];dy=ys[b]-ys[a]
            seglen=hypot(dx,dy)
            seglen>0 || throw(ArgumentError(
                "$caller: Surface[$t] wall curve has a coincident vertex"))
            offsets=_bl_field_offsets(spec,caller)
            probe=min(0.5*offsets[1],0.05*seglen)
            lx=-dy/seglen;ly=dx/seglen
            mx=0.5*(xs[a]+xs[b])+probe*lx
            my=0.5*(ys[a]+ys[b])+probe*ly
            inside=_bl_point_in_domain(mx,my,xs,ys,wires)
            if !inside
                reverse!(chain)
                a=chain[1];b=chain[2]
            end
            wall_coords=Matrix{Float64}(undef,3,nchain)
            for (i,v) in enumerate(chain)
                wall_coords[1,i],wall_coords[2,i],wall_coords[3,i]=coords3d[v]
            end
            wall_segs=Matrix{Int32}(undef,2,closed_run ? nchain : nchain-1)
            for i in 1:(nchain-1)
                wall_segs[1,i]=i;wall_segs[2,i]=i+1
            end
            if closed_run
                wall_segs[1,nchain]=nchain;wall_segs[2,nchain]=1
            end
            wall_mesh=Mesh(wall_coords;segs=wall_segs)
            # Run-end columns ride the adjacent non-wall segment —
            # upstream's `addColumnAtTheEndOfTheBL` makes the end vertex's
            # column the adjacent edge's own node chain (graded by
            # `addBoundaryLayerPoints` upstream). Bounded: the consumed
            # depth `offsets[end]` must fit inside the single adjacent
            # segment, and two run ends must not consume overlapping
            # spans of the same segment.
            vdirs=nothing
            if !closed_run
                p_start=first
                p_end=last+1
                v_start=loop[mod1(start0+p_start-1,nseg)]
                v_end=loop[mod1(start0+p_end-1,nseg)]
                prev_seg=mod1(start0+first-2,nseg)
                next_seg=mod1(start0+last,nseg)
                seg_spec[prev_seg]==0 || throw(ArgumentError(
                    "$caller: Surface[$t] has adjacent boundary-layer " *
                    "runs — upstream requires separated layers"))
                seg_spec[next_seg]==0 || throw(ArgumentError(
                    "$caller: Surface[$t] has adjacent boundary-layer " *
                    "runs — upstream requires separated layers"))
                # Walk the covered span on the preceding segment (backward
                # from v_start) and on the following segment (forward from
                # v_end). Bounded: the span stays on one adjacent curve,
                # never reaches a wall segment, and stays straight.
                tops_s=_bl_walk_covered!(
                    -1,p_start,loop_curve[prev_seg],offsets[end],
                    spec,caller,t,covered_positions,claimed,seg_cover,
                    snapped_coords,covered_owner,seg_spec,loop,loop_curve,
                    xs,ys,coords3d,start0,nseg,new_vertex!,vertex_id!,ax)
                tops_e=_bl_walk_covered!(
                    +1,p_end,loop_curve[next_seg],offsets[end],
                    spec,caller,t,covered_positions,claimed,seg_cover,
                    snapped_coords,covered_owner,seg_spec,loop,loop_curve,
                    xs,ys,coords3d,start0,nseg,new_vertex!,vertex_id!,ax)
                top_override[v_start]=tops_s[1]
                top_override[v_end]=tops_e[1]
                dir_s=tops_s[2];dir_e=tops_e[2]
                vdirs=Dict{Int,NTuple{3,Float64}}()
                # `chain` may have been reversed for orientation; attach
                # each direction to the slot holding its end vertex.
                for (i,v) in enumerate(chain)
                    v==v_start && (vdirs[i]=dir_s)
                    v==v_end && (vdirs[i]=dir_e)
                end
            end
            fans=Int[]
            for (i,v) in enumerate(chain)
                pid=get(point_of,v,0)
                pid in spec.fan_set || continue
                interior=closed_run || (i>1 && i<nchain)
                interior || continue
                push!(fans,i)
            end
            strip=blmod._bl_2d_offsets(
                wall_mesh,offsets;fans=fans,
                fan_elements=spec.fan_elements,
                vertex_dirs=vdirs,
                plane_normal=plane_normal,caller=caller)
            mesh=strip.mesh
            # tops keys on the wall node's exact 3-D coordinate. Run-end
            # tops are the covered-span top vertices the walk resolved —
            # snap the strip's computed top onto that exact coordinate.
            for (i,v) in enumerate(chain)
                vid=get(top_override,v,0)
                if vid!=0
                    push!(snapped_coords,new_coords[vid])
                    top_of[v]=vid
                else
                    top=strip.tops[coords3d[v]]
                    top_of[v]=new_vertex!(top[ax[1]],top[ax[2]],top)
                end
            end
            base=length(strip_coords)
            for i in axes(mesh.coords,2)
                push!(strip_coords,(mesh.coords[1,i],mesh.coords[2,i],
                                    mesh.coords[3,i]))
            end
            for block in mesh.blocks
                if block.msh==2
                    for c in axes(block.nodes,2)
                        push!(strip_tris,(base+Int(block.nodes[1,c]),
                            base+Int(block.nodes[2,c]),
                            base+Int(block.nodes[3,c])))
                    end
                elseif block.msh==3
                    for c in axes(block.nodes,2)
                        push!(strip_quads,(base+Int(block.nodes[1,c]),
                            base+Int(block.nodes[2,c]),
                            base+Int(block.nodes[3,c]),
                            base+Int(block.nodes[4,c])))
                    end
                else
                    throw(ErrorException(
                        "$caller: boundary-layer strip emitted unexpected " *
                        "cell type $(block.msh)"))
                end
            end
        end
        # Rebuild the wire: walk the rotated order; at each run emit only
        # the top polyline. Run-end vertices are covered by the strip
        # side (their columns ride the adjacent segments), so they drop
        # out of the fill boundary — the top of each end column lands
        # on the adjacent segment and the wire resumes past it.
        new_loop=Int[]
        k=1
        while k<=nseg
            sk=mod1(start0+k-1,nseg)
            run_index=0
            for (ri,run) in enumerate(runs)
                if run[2]<=k<=run[3]
                    run_index=ri
                    break
                end
            end
            if run_index==0
                k in covered_positions ||
                    push!(new_loop,vertex_id!(loop[sk]))
                k+=1
                continue
            end
            si,first,last,closed_flag=runs[run_index]
            if !iszero(closed_flag)
                for p in 1:nseg
                    v=loop[mod1(start0+p-1,nseg)]
                    push!(new_loop,top_of[v])
                end
                break
            end
            for p in first:last+1
                v=loop[mod1(start0+p-1,nseg)]
                push!(new_loop,top_of[v])
            end
            # The run's end vertex sits at rotated position `last+1`,
            # consumed by the strip side — resume on the next segment.
            k=last+2
        end
        push!(reduced_wires,new_loop)
    end
    segs2=Tuple{Int,Int}[]
    for wire in reduced_wires
        nloop=length(wire)
        nloop>=3 || throw(ArgumentError(
            "$caller: Surface[$t] boundary layer leaves a degenerate wire"))
        for k in 1:nloop
            i,j=wire[k],wire[mod1(k+1,nloop)]
            i==j || (new_xs[i]==new_xs[j] && new_ys[i]==new_ys[j]) ||
                push!(segs2,(i,j))
        end
    end
    # Embedded entities ride the same constraints the ordinary fill carries —
    # remap their PSLG ids through `vertex_id!`.
    embedded_vertex_remap=Dict{Int,Int}()
    internal2=Tuple{Int,Int}[]
    for (a,b) in internal
        na=get!(embedded_vertex_remap,a) do
            vertex_id!(a)
        end
        nb=get!(embedded_vertex_remap,b) do
            vertex_id!(b)
        end
        na==nb && continue
        push!(internal2,(na,nb))
    end
    for (edim,etag) in embedded
        edim==0 || continue
        haskey(m.points,etag) || throw(ArgumentError(
            "$caller: unknown embedded Point[$etag]"))
        # `_bl_surface_pslg` already added the point vertex; find it through
        # `index` so it lands in the fill's vertex set.
        vertex_id!(index[etag])
    end
    triangulation=constrained_delaunay(
        new_xs,new_ys,segs2;internal_segments=internal2)
    fill=to_mesh(
        triangulation;interior=classify_interior(triangulation))
    # Merge every vertex class into one coordinate pool: boundary and top
    # vertices first (stable order), then strip interiors.
    pool=Dict{NTuple{3,Float64},Int32}()
    out_coords=NTuple{3,Float64}[]
    out_owners=Tuple{Int,Int32}[]
    out_params=Union{Nothing,Vector{Float64}}[]
    function pool_id!(coordinate,owner,param)
        return get!(pool,coordinate) do
            push!(out_coords,coordinate)
            push!(out_owners,owner);push!(out_params,param)
            Int32(length(out_coords))
        end
    end
    new_to_pool=Vector{Int32}(undef,length(new_coords))
    for i in eachindex(new_coords)
        new_to_pool[i]=pool_id!(new_coords[i],new_owners[i],new_params[i])
    end
    fill_to_pool=Vector{Int32}(undef,nnodes(fill))
    for node in 1:nnodes(fill)
        key=(fill.coords[1,node]==0 ? 0.0 : fill.coords[1,node],
             fill.coords[2,node]==0 ? 0.0 : fill.coords[2,node])
        vid=get(vertex_index,key,nothing)
        vid===nothing && throw(ErrorException(
            "$caller: Surface[$t] fill introduced an unmapped coordinate"))
        fill_to_pool[node]=new_to_pool[vid]
    end
    snap_tol=1e-9*max(1.0,maximum(xs)-minimum(xs),
                      maximum(ys)-minimum(ys))
    strip_to_pool=Vector{Int32}(undef,length(strip_coords))
    for i in eachindex(strip_coords)
        coordinate=strip_coords[i]
        for s in snapped_coords
            (abs(coordinate[1]-s[1])<=snap_tol &&
             abs(coordinate[2]-s[2])<=snap_tol &&
             abs(coordinate[3]-s[3])<=snap_tol) &&
                (coordinate=s;break)
        end
        ctag=get(covered_owner,coordinate,Int32(0))
        strip_to_pool[i]=ctag==0 ?
            pool_id!(coordinate,(2,Int32(t)),
                     [coordinate[ax[1]],coordinate[ax[2]]]) :
            pool_id!(coordinate,(1,ctag),nothing)
    end
    out_tris=NTuple{3,Int32}[]
    out_quads=NTuple{4,Int32}[]
    for tri in strip_tris
        push!(out_tris,(strip_to_pool[tri[1]],strip_to_pool[tri[2]],
                        strip_to_pool[tri[3]]))
    end
    keep_quads=_model_surface_recombined(m,t) ||
        any(spec->spec.quads,specs)
    for quad in strip_quads
        a,b,c,d=quad
        if keep_quads
            push!(out_quads,(strip_to_pool[a],strip_to_pool[b],
                             strip_to_pool[c],strip_to_pool[d]))
        else
            # `quadsToTriangles` on the strip split: the emitted quad's
            # (b1,b2,t2,t1) frame splits along the b1–t2 diagonal.
            push!(out_tris,(strip_to_pool[a],strip_to_pool[b],
                            strip_to_pool[c]))
            push!(out_tris,(strip_to_pool[a],strip_to_pool[c],
                            strip_to_pool[d]))
        end
    end
    for c in 1:ntris(fill)
        push!(out_tris,(fill_to_pool[Int(fill.tris[1,c])],
                        fill_to_pool[Int(fill.tris[2,c])],
                        fill_to_pool[Int(fill.tris[3,c])]))
    end
    coords_out=Matrix{Float64}(undef,3,length(out_coords))
    for i in eachindex(out_coords)
        coords_out[1,i],coords_out[2,i],coords_out[3,i]=out_coords[i]
    end
    isempty(out_tris) && isempty(out_quads) && throw(ErrorException(
        "$caller: Surface[$t] boundary-layer pass produced no cells"))
    result=if isempty(out_quads)
        tris=Matrix{Int32}(undef,3,length(out_tris))
        for (i,tri) in enumerate(out_tris)
            tris[1,i],tris[2,i],tris[3,i]=tri
        end
        Mesh(coords_out;tris=tris)
    else
        tris=Matrix{Int32}(undef,3,length(out_tris))
        for (i,tri) in enumerate(out_tris)
            tris[1,i],tris[2,i],tris[3,i]=tri
        end
        quads=Matrix{Int32}(undef,4,length(out_quads))
        for (i,quad) in enumerate(out_quads)
            quads[1,i],quads[2,i],quads[3,i],quads[4,i]=quad
        end
        blocks=ElementBlock[]
        isempty(out_tris) || push!(blocks,ElementBlock(2,tris))
        push!(blocks,ElementBlock(3,quads))
        MixedMesh(coords_out,blocks)
    end
    diagnostic=validate(result)
    diagnostic.ok || throw(ErrorException(
        "$caller: Surface[$t] boundary-layer mesh is invalid — " *
        join(diagnostic.messages,"; ")))
    return result,out_owners,out_params
end

# `_occ_curved_surface_production` for `Mesh|MixedMesh` strip output —
# boundary vertices reuse the consumed lower product's actual UIDs
# (matched by owning entity and bitwise coordinate, like upstream's
# `meshGFace` vertex reassociation) while strip/fill interiors mint fresh
# labels from the authority's high water.
function _bl_surface_production(m::GeoModel,t::Int,mesh,
                                owners::Vector{Tuple{Int,Int32}},
                                params,input,limit,caller::AbstractString)
    tags=Vector{UInt64}(undef,nnodes(mesh))
    retained=falses(nnodes(mesh))
    lookup=Dict{Tuple{Int,Int32,NTuple{3,Float64}},UInt64}()
    source=input===nothing ? nothing : input.mesh
    if source!==nothing && source.entity_data!==nothing
        data=source.entity_data
        for row in 1:nnodes(source)
            entity=data.node_entities[row]
            key=(entity[1],entity[2],
                 (source.coords[1,row],source.coords[2,row],
                  source.coords[3,row]))
            haskey(lookup,key) || (lookup[key]=data.external_node_tags[row])
        end
    end
    next=UInt64(limit)
    for k in 1:nnodes(mesh)
        owner=owners[k]
        point=(mesh.coords[1,k],mesh.coords[2,k],mesh.coords[3,k])
        tag=get(lookup,(owner[1],owner[2],point),UInt64(0))
        if tag!=0
            tags[k]=tag;retained[k]=true
        else
            next+=UInt64(1)
            iszero(next) && throw(ArgumentError(
                "$caller: Surface[$t] vertex births exhaust UInt64 labels"))
            tags[k]=next
        end
    end
    return _ModelSurfaceProduction(mesh,tags,owners,params,retained,next)
end
