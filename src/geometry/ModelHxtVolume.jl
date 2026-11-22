# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: GPL-2.0-or-later
# Contains Gmsh 4.15.2 adaptations: src/mesh/meshGRegionHxt.cpp (Gmsh2Hxt,
# nodalSizesCallBack, Hxt2Gmsh) and src/mesh/BackgroundMeshTools.cpp
# (BGM_MeshSizeWithoutScaling). See LICENSES/README.md and NOTICES.md#gmsh.

# `Mesh.Algorithm3D = 10`: one Region meshed by the HXT port, with the input
# meshGRegionHxt builds -- every vertex of the bounding Surfaces' triangles and
# Curves' lines in MVertex-number order, the lines Curve by Curve and the
# triangles Surface by Surface in tag order -- and Gmsh's size callback.
function _model_hxt_volume(m::GeoModel,t::Int,caller;size_field=nothing,
                           _mesh_input=nothing,options,births=nothing)
    _mesh_input isa _ModelMeshAuthority || throw(ArgumentError(
        "$caller: Mesh.Algorithm3D=10 (HXT) needs the generated surface meshes"))
    isempty(get(m.embeds,(3,t),NTuple{2,Int}[])) || throw(ArgumentError(
        "$caller: Mesh.Algorithm3D=10 (HXT) with embedded entities in Volume[$t] is not implemented"))
    boundaries=_model_volume_boundary_surfaces(m,t,caller)
    surfaces=sort!(unique(abs.(boundaries)))
    curves=sort!(unique(reduce(vcat,(_model_surface_curves(m,s) for s in surfaces);init=Int[])))
    # every vertex, keyed by its born UID
    coordinate=Dict{UInt64,NTuple{3,Float64}}()
    lines=Tuple{UInt64,UInt64,Int}[]
    triangles=Tuple{UInt64,UInt64,UInt64,Int}[]
    function blocks(dimension,tag)
        product=_model_mesh_product(_mesh_input,dimension,tag;caller,producer_order=true,shared=true)
        product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
            "$caller: entity ($dimension, $tag) has no generated product"))
        return product
    end
    for c in curves
        # a degenerate Curve (a sphere's pole) carries no lines
        product=_model_mesh_product(_mesh_input,1,c;caller,producer_order=true,shared=true)
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
        product=blocks(2,s)
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
    for (i,u) in enumerate(uids),k in 1:3
        coords[4(i-1)+k]=coordinate[u][k]
    end
    hm=Hxt.Mesh(coords)
    for (a,b,c) in lines
        push!(hm.line_node,row[a],row[b]);push!(hm.line_color,UInt32(c))
    end
    for (a,b,c,s) in triangles
        push!(hm.tri_node,row[a],row[b],row[c]);push!(hm.tri_color,UInt32(s))
    end
    push!(hm.surfaces_per_volume,UInt32.(abs.(boundaries)))
    for s in surfaces
        push!(hm.curves_per_surface,UInt32.(_model_surface_curves(m,s)))
    end
    # nodalSizesCallBack: BGM_MeshSizeWithoutScaling (no Min/Max/factor),
    # capped by CTX::lc, and by the interpolated size when the boundary is
    # extended into the volume.
    interpolate=m.meshing.lc_extend_from_boundary>0 || m.meshing.lc_extend_from_boundary==-3
    lc_global=options.ctx_lc
    function size(x,y,z,interpolated,color)
        raw=options.inner_field!==nothing ? field_value(options.inner_field,x,y,z,(3,t)) :
            size_field!==nothing ? size_at(size_field,x,y,z,(3,t)) : GMSH_MAX_SIZE
        options.callback===nothing || (raw=_model_algorithm_callback(
            options.callback,3,t,x,y,z,raw,caller))
        return interpolate && interpolated>0.0 ? min(interpolated,min(lc_global,raw)) :
            min(lc_global,raw)
    end
    Hxt.tet_mesh!(hm;min_size=options.lc_min,max_size=options.lc_max,factor=options.lc_factor,
        size_callback=size,quality_min=options.optimize_threshold,optimize=options.optimize>0)
    # Hxt2Gmsh: the cells of the Region in order; a vertex not on the
    # boundary becomes a new MVertex when a cell first names it.
    nb=length(uids)
    local_row=Dict{UInt32,Int32}()
    newcoords=NTuple{3,Float64}[]
    cells=NTuple{4,Int32}[]
    for c in 0:hm.ntet-1
        hm.color[c+1]==0 || continue
        v=ntuple(4) do j
            n=hm.node[4c+j]
            n<nb && return Int32(n+1)
            get!(local_row,n) do
                push!(newcoords,Hxt.point(hm,n))
                Int32(nb+length(newcoords))
            end
        end
        push!(cells,v)
    end
    all_coords=Matrix{Float64}(undef,3,nb+length(newcoords))
    for i in 1:nb,k in 1:3
        all_coords[k,i]=coords[4(i-1)+k]
    end
    for (i,p) in enumerate(newcoords),k in 1:3
        all_coords[k,nb+i]=p[k]
    end
    tets=Matrix{Int32}(undef,4,length(cells))
    for (i,c) in enumerate(cells)
        # HXT's cells have orient3d(a,b,c,d) < 0: positive volume, as Mesh's
        for k in 1:4;tets[k,i]=c[k];end
    end
    mesh=Mesh(all_coords;tets)
    if births!==nothing
        births.source=Mesh(all_coords[:,1:nb])
        high=max(births.node_high_water,maximum(uids;init=UInt64(0)))
        births.node_tags=vcat(uids,collect(UInt64,high+1:high+length(newcoords)))
        births.node_high_water=high+UInt64(length(newcoords))
        _model_native_volume_cell_births!(births,mesh)
    end
    return mesh
end
