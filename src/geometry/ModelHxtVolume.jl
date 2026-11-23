# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: GPL-2.0-or-later
# Contains Gmsh 4.15.2 adaptations: src/mesh/meshGRegionHxt.cpp (Gmsh2Hxt,
# nodalSizesCallBack, Hxt2Gmsh), src/mesh/Generator.cpp
# (FindConnectedRegions) and src/mesh/BackgroundMeshTools.cpp
# (BGM_MeshSizeWithoutScaling). See LICENSES/README.md and NOTICES.md#gmsh.

# `Mesh.Algorithm3D = 10`: the Regions connected to `t` through shared
# Surfaces (FindConnectedRegions) are meshed by one call of the HXT port, with
# the input meshGRegionHxt builds (Gmsh2Hxt) -- every vertex of the Surfaces'
# triangles, the Curves' lines and the Regions' embedded Points in
# MVertex-number order, the lines Curve by Curve and the triangles Surface by
# Surface in tag order -- and Gmsh's size callback. Each Region then takes its
# own cells (Hxt2Gmsh); the group's run is kept for its other Regions.
function _model_hxt_group(m::GeoModel,t::Int,caller)
    meshed(v)=!haskey(m.meshing.transfinite_volumes,v) && !haskey(m.meshing.extrude,(3,v))
    regions_of=Dict{Int,Vector{Int}}()
    for v in sort!(collect(keys(m.volumes)))
        for s in abs.(_model_volume_boundary_surfaces(m,v,caller))
            push!(get!(()->Int[],regions_of,s),v)
        end
    end
    seen=Set{Int}([t]);stack=[t]
    while !isempty(stack)
        r=pop!(stack)
        for s in abs.(_model_volume_boundary_surfaces(m,r,caller)),other in get(regions_of,s,Int[])
            other in seen || (push!(seen,other);push!(stack,other))
        end
    end
    return sort!([v for v in seen if meshed(v) || v==t])
end

function _model_hxt_volume(m::GeoModel,t::Int,caller;size_field=nothing,
                           _mesh_input=nothing,options,births=nothing)
    _mesh_input isa _ModelMeshAuthority || throw(ArgumentError(
        "$caller: Mesh.Algorithm3D=10 (HXT) needs the generated surface meshes"))
    group=_model_hxt_group(m,t,caller)
    cache=get(task_local_storage(),:tessella_hxt_runs,nothing)
    if !(cache isa Tuple && cache[1]===_mesh_input)
        cache=(_mesh_input,Dict{Vector{Int},Any}())
        task_local_storage(:tessella_hxt_runs,cache)
    end
    run=get(cache[2],group,nothing)
    if run===nothing
        run=(result=_model_hxt_run(m,group,caller;size_field,_mesh_input,options),
             pending=Set(group))
        length(group)>1 && (cache[2][group]=run)
    end
    delete!(run.pending,t)
    isempty(run.pending) && delete!(cache[2],group)
    hm,uids=run.result
    color=UInt32(findfirst(==(t),group)-1)
    # Hxt2Gmsh: the cells of the Region in order; a vertex not on the
    # boundary becomes a new MVertex when a cell first names it. The Region's
    # rows are its own input vertices (UID order), then its new ones.
    nb=length(uids)
    used=falses(nb)
    for c in 0:hm.ntet-1
        hm.color[c+1]==color || continue
        for j in 1:4
            n=hm.node[4c+j]
            n<nb && (used[n+1]=true)
        end
    end
    local_row=zeros(Int32,hm.nv)
    own=Int[]
    for i in 1:nb
        used[i] || continue
        push!(own,i);local_row[i]=Int32(length(own))
    end
    nown=length(own)
    newcoords=NTuple{3,Float64}[]
    cells=NTuple{4,Int32}[]
    sizehint!(cells,hm.ntet)
    for c in 0:hm.ntet-1
        hm.color[c+1]==color || continue
        v=ntuple(4) do j
            n=hm.node[4c+j]
            if local_row[n+1]==0
                push!(newcoords,Hxt.point(hm,n))
                local_row[n+1]=Int32(nown+length(newcoords))
            end
            return local_row[n+1]
        end
        push!(cells,v)
    end
    all_coords=Matrix{Float64}(undef,3,nown+length(newcoords))
    for (r,i) in enumerate(own)
        p=Hxt.point(hm,i-1)
        for k in 1:3;all_coords[k,r]=p[k];end
    end
    for (i,p) in enumerate(newcoords),k in 1:3
        all_coords[k,nown+i]=p[k]
    end
    tets=Matrix{Int32}(undef,4,length(cells))
    for (i,c) in enumerate(cells)
        # HXT's cells have orient3d(a,b,c,d) < 0: positive volume, as Mesh's
        for k in 1:4;tets[k,i]=c[k];end
    end
    mesh=Mesh(all_coords;tets)
    if births!==nothing
        births.source=Mesh(all_coords[:,1:nown])
        high=max(births.node_high_water,maximum(uids;init=UInt64(0)))
        births.node_tags=vcat(uids[own],collect(UInt64,high+1:high+length(newcoords)))
        births.node_high_water=high+UInt64(length(newcoords))
        _model_native_volume_cell_births!(births,mesh)
    end
    return mesh
end

# Gmsh2Hxt and hxtTetMesh for the Regions `group` (tag order): the meshed HXT
# structure and the UIDs of its input vertices (its first rows).
function _model_hxt_run(m::GeoModel,group::Vector{Int},caller;size_field,_mesh_input,options)
    region_faces=Dict{Int,Vector{Int}}()
    for r in group
        faces=abs.(_model_volume_boundary_surfaces(m,r,caller))
        for (d,e) in get(m.embeds,(3,r),NTuple{2,Int}[])
            d==2 && push!(faces,e)
        end
        region_faces[r]=faces
    end
    surfaces=sort!(unique(reduce(vcat,values(region_faces);init=Int[])))
    surface_curves=Dict{Int,Vector{Int}}()
    for s in surfaces
        curves=_model_surface_curves(m,s)
        for (d,e) in get(m.embeds,(2,s),NTuple{2,Int}[])
            d==1 && push!(curves,e)
        end
        surface_curves[s]=curves
    end
    curves=reduce(vcat,values(surface_curves);init=Int[])
    points=Int[]
    for r in group,(d,e) in get(m.embeds,(3,r),NTuple{2,Int}[])
        d==1 && push!(curves,e)
        d==0 && push!(points,e)
    end
    curves=sort!(unique(curves))
    # every vertex, keyed by its born UID
    coordinate=Dict{UInt64,NTuple{3,Float64}}()
    prescribed=Dict{UInt64,Float64}()
    point_nodes=Tuple{UInt64,Int}[]
    lines=Tuple{UInt64,UInt64,Int}[]
    triangles=Tuple{UInt64,UInt64,UInt64,Int}[]
    product_of(dimension,tag)=_model_mesh_product(_mesh_input,dimension,tag;caller,
        producer_order=true,shared=true)
    for p in points
        product=product_of(0,p)
        (product isa MixedMesh && product.entity_data!==nothing) || throw(ArgumentError(
            "$caller: embedded Point $p has no generated product"))
        uid=product.entity_data.external_node_tags
        for block in product.blocks,column in axes(block.nodes,2)
            a=block.nodes[1,column]
            coordinate[uid[a]]=_model_mesh_coordinate(product,a)
            push!(point_nodes,(uid[a],p))
            # Mesh.MeshSizeFromPoints: the Point's prescribed size
            prescribed_size=Float64(get(m.point_size,p,GMSH_MAX_SIZE))
            options.from_points && prescribed_size<GMSH_MAX_SIZE &&
                (prescribed[uid[a]]=prescribed_size)
        end
    end
    for c in curves
        # a degenerate Curve (a sphere's pole) carries no lines
        product=product_of(1,c)
        (product isa MixedMesh && product.entity_data!==nothing) || continue
        uid=product.entity_data.external_node_tags
        for block in product.blocks
            block.msh==1 || continue
            for column in axes(block.nodes,2)
                a=block.nodes[1,column];b=block.nodes[2,column]
                coordinate[uid[a]]=_model_mesh_coordinate(product,a)
                coordinate[uid[b]]=_model_mesh_coordinate(product,b)
                push!(lines,(uid[a],uid[b],c))
            end
        end
    end
    for s in surfaces
        product=product_of(2,s)
        product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
            "$caller: entity (2, $s) has no generated product"))
        uid=product.entity_data.external_node_tags
        for block in product.blocks
            block.msh==2 || throw(ArgumentError(
                "$caller: Surface $s contains elements which are not triangles. " *
                "The HXT 3D meshing algorithm only supports triangles."))
            for column in axes(block.nodes,2)
                v=ntuple(k->block.nodes[k,column],3)
                for k in 1:3
                    coordinate[uid[v[k]]]=_model_mesh_coordinate(product,v[k])
                end
                push!(triangles,(uid[v[1]],uid[v[2]],uid[v[3]],s))
            end
        end
    end
    uids=sort!(collect(keys(coordinate)))
    row=Dict(u=>UInt32(i-1) for (i,u) in enumerate(uids))
    coords=zeros(4length(uids))
    for (i,u) in enumerate(uids)
        for k in 1:3;coords[4(i-1)+k]=coordinate[u][k];end
        coords[4i]=get(prescribed,u,0.0)
    end
    hm=Hxt.Mesh(coords)
    for (u,p) in point_nodes
        push!(hm.point_node,row[u]);push!(hm.point_color,UInt32(p))
    end
    for (a,b,c) in lines
        push!(hm.line_node,row[a],row[b]);push!(hm.line_color,UInt32(c))
    end
    for (a,b,c,s) in triangles
        push!(hm.tri_node,row[a],row[b],row[c]);push!(hm.tri_color,UInt32(s))
    end
    for r in group
        push!(hm.surfaces_per_volume,UInt32.(region_faces[r]))
    end
    for s in surfaces
        push!(hm.curves_per_surface,UInt32.(surface_curves[s]))
    end
    # nodalSizesCallBack: BGM_MeshSizeWithoutScaling in the point's Region
    # (no Min/Max/factor), capped by CTX::lc, and by the interpolated size when
    # the boundary is extended into the volume.
    interpolate=m.meshing.lc_extend_from_boundary>0 || m.meshing.lc_extend_from_boundary==-3
    lc_global=Float64(options.ctx_lc)
    inner_field=options.inner_field;callback=options.callback
    size=if inner_field===nothing && size_field===nothing && callback===nothing
        # no field: every candidate gets CTX::lc (or the interpolated size)
        lc=min(lc_global,Float64(GMSH_MAX_SIZE))
        (x,y,z,interpolated,color)->interpolate && interpolated>0.0 ? min(interpolated,lc) : lc
    else
        function (x,y,z,interpolated,color)
            r=group[Int(color)+1]
            raw=inner_field!==nothing ? Float64(field_value(inner_field,x,y,z,(3,r))) :
                size_field!==nothing ? Float64(size_at(size_field,x,y,z,(3,r))) : GMSH_MAX_SIZE
            callback===nothing || (raw=Float64(_model_algorithm_callback(
                callback,3,r,x,y,z,raw,caller)))
            return (interpolate && interpolated>0.0 ? min(interpolated,min(lc_global,raw)) :
                min(lc_global,raw))::Float64
        end
    end
    recovery=hm->try
        _model_hxt_boundary_recovery(hm)
    catch err
        err isa ModelNativeTetgenArena.Unsupported || rethrow()
        throw(ArgumentError("$caller: HXT boundary recovery of Volumes $group is not supported: $(err.reason)"))
    end
    Hxt.tet_mesh!(hm;min_size=Float64(options.lc_min),max_size=Float64(options.lc_max),
        factor=Float64(options.lc_factor),size_callback=size,
        quality_min=Float64(options.optimize_threshold),optimize=options.optimize>0,
        boundary_recovery=recovery)
    return hm,uids
end

# hxt_boundary_recovery on the HXT mesh: the TetGen/BR port with Gmsh's
# Mesh.ToleranceInitialDelaunay (1e-12) for both of HXT's tolerance switches.
function _model_hxt_boundary_recovery(hm::Hxt.Mesh)
    xyz=NTuple{3,Float64}[Hxt.point(hm,v) for v in 0:hm.nv-1]
    return ModelNativeTetgenArena.BR.hxt_recovery(xyz,hm.node,hm.neigh,hm.ntet,
        hm.tri_node,hm.tri_color,hm.line_node,hm.line_color;epsilon=1e-12,facet_overlap=1e-12)
end
