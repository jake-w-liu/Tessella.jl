# The insertion background and the propagated boundary-node size are
# independent arrays in MTet4. Do not interpolate nearest boundary sizes
# into this background evaluator: the cavity handles their propagation.
function _model_algorithm_size_options(m::GeoModel,options,caller)
    options===nothing && return _model_default_mesh1d_options(m,caller)
    options isa _ModelMesh1DOptions || throw(ArgumentError(
        "$caller: native algorithm sizing requires a mesh option bundle"))
    return options
end

# Existing boundary/embedded vertices retain their actual owner dimension.
# Native public IDs bind these rows to the caller's generation input; new
# volume vertices default to the interior movement role.
function _model_native_volume_optimizer_movable(mesh,births,input)
    fixed=falses(nnodes(mesh))
    for face in first(boundary_faces(mesh.tets)),node in face
        fixed[node]=true
    end
    if births!==nothing && input!==nothing && hasproperty(input,:mesh)
        actual=input.mesh
        data=actual isa MixedMesh ? actual.entity_data : nothing
        if data!==nothing
            rows=Dict(tag=>row for (row,tag) in enumerate(data.external_node_tags))
            for (row,tag) in enumerate(births.node_tags)
                row<=nnodes(mesh) || break
                binding=get(rows,tag,0)
                binding!=0 && data.node_entities[binding][1]<3 && (fixed[row]=true)
            end
        end
    end
    return .!fixed
end

function _model_native_optimize_volume(mesh,options,births,input,faces,edges,caller;
                                       protected_nodes=())
    if births!==nothing
        if births.cell_source===nothing
            # The inherited factory fill creates its cells before optimization.
            # Bind those actual births to the caller's current counter so that
            # removed cells and transient swaps still consume element numbers.
            _model_native_volume_cell_births!(births,mesh)
        elseif births.cell_source!==mesh && births.cell_source.tets==mesh.tets
            # Smoothing copies coordinates without replacing the MElements.
            births.cell_source=mesh
        end
    end
    options.optimize_netgen>0 && ntets(mesh)>0 && throw(ArgumentError(
        "$caller: the native Netgen tetrahedron optimizer is not implemented"))
    options.optimize>0 && options.optimize_threshold>0.0 || return mesh
    receipt=births!==nothing && births.cell_source===mesh ? births : nothing
    movement=_model_native_volume_optimizer_movable(mesh,receipt,input)
    for node in protected_nodes
        movement[node]=false
    end
    for face in faces,node in face
        movement[node]=false
    end
    for edge in edges,node in edge
        movement[node]=false
    end
    for _ in 1:options.optimize
        numbers=receipt===nothing ? collect(UInt64,1:ntets(mesh)) : receipt.element_tags
        high_water=receipt===nothing ? maximum(numbers;init=UInt64(0)) : receipt.element_high_water
        product=_gmsh_optimize_tetrahedra(mesh;movable=movement,
            element_numbers=numbers,element_high_water=high_water,
            threshold=options.optimize_threshold,embedded_faces=faces,embedded_edges=edges)
        mesh=product.mesh
        if receipt!==nothing
            receipt.cell_source=mesh
            receipt.element_tags=product.element_numbers
            receipt.element_high_water=product.element_high_water
        end
    end
    return mesh
end

struct _ModelAlgorithmBackgroundSize
    model::GeoModel
    dimension::Int
    tag::Int
    options::_ModelMesh1DOptions
    field::Union{Nothing,AbstractSizeField}
    caller::String
end

function _ModelAlgorithmBackgroundSize(m::GeoModel,dimension::Int,tag::Int,
                                       options,field,caller)
    return _ModelAlgorithmBackgroundSize(m,dimension,tag,options,field,String(caller))
end

function _model_algorithm_callback(callback,dimension,tag,x,y,z,value,caller)
    return _apply_size_callback(callback,dimension,tag,x,y,z,value,caller)
end

function (background::_ModelAlgorithmBackgroundSize)(x,y,z)
    m=background.model
    options=background.options
    dimension=background.dimension;tag=background.tag
    entity=(dimension,tag)
    # Raw background comes before the callback; CTX::lc, Min/Max and factors
    # are applied after it by BGM_MeshSize. A callback may enlarge the raw
    # minimum as well as tighten it.
    value=if options.inner_field!==nothing
        field_value(options.inner_field,x,y,z,entity)
    elseif background.field!==nothing
        size_at(background.field,x,y,z,entity)
    else
        GMSH_MAX_SIZE
    end
    if dimension==2 && options.from_curvature>0
        parameters=model_parametrization(m,2,tag,[x,y,z])
        curvature=only(model_curvature(m,2,tag,parameters))
        curvature>0 &&
            (value=min(value,2π/curvature/options.from_curvature))
    end
    options.callback===nothing || (value=_model_algorithm_callback(
        options.callback,dimension,tag,x,y,z,value,background.caller))
    value=value<options.ctx_lc ? value : options.ctx_lc
    value=max(value,options.lc_min)
    value=min(value,options.lc_max)
    if value<=0
        _model_mesh_error(options,
            "Wrong mesh element size lc = $value on entity ($dimension,$tag)")
        value=options.ctx_lc
    end
    return value*options.lc_factor
end

struct _ModelSurfaceBackgroundSize
    frame::_SurfaceMetricFrame
    background::_ModelAlgorithmBackgroundSize
end

function (background::_ModelSurfaceBackgroundSize)(x,y)
    point=_surface_metric_lift(background.frame,x,y)
    return background.background(point...)
end

function _model_volume_embedded_sizes(m::GeoModel,tag::Int,mesh::Mesh,
                                      curves,options,caller)
    result=Dict{Int32,Float64}()
    # insertVerticesInRegion seeds prescribed GVertex sizes and the smallest
    # incident actual GEdge link, then setLcs updates them from face/tet edges.
    # These propagated values are separate from BGM_MeshSize.
    for (dimension,entity) in get(m.embeds,(3,tag),NTuple{2,Int}[])
        dimension==0 || continue
        node=_model_projection_embedded_point_node(m,mesh,entity,1e-12,caller)
        result[node]=min(get(result,node,Inf),_model_vertex_lc(m,entity))
    end
    edges=_tet_edge_set(mesh)
    for curve in curves
        entries=_model_projection_tet_curve_nodes(m,mesh,curve,edges,caller)
        closed=m.curves[curve][1]==m.curves[curve][2]
        links=closed ? length(entries) : length(entries)-1
        for index in 1:links
            first=Int32(entries[index][2])
            second=Int32(entries[mod1(index+1,length(entries))][2])
            a=_pt3(mesh,first);b=_pt3(mesh,second)
            value=hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3])
            result[first]=min(get(result,first,Inf),value)
            result[second]=min(get(result,second,Inf),value)
        end
    end
    return result
end

# Native boundary vertices are immutable. A numeric perturbation used by the
# initial tetrahedralizer must not leak into the emitted PLC or trigger a
# centroid fan before constrained boundary recovery has been attempted.
function _model_native_quality_volume_fill(surface,order,kernel,caller;protected_facets=nothing)
    work=ModelNativeTetgenArena.QualityWorkspace()
    ModelNativeTetgenArena.has_bad_quality(kernel,work) || return nothing
    inverse=Vector{Int32}(undef,length(order))
    for (native,row) in enumerate(order);inverse[row]=Int32(native);end
    facets=NTuple{3,Int32}[(inverse[surface.tris[1,column]],
        inverse[surface.tris[2,column]],inverse[surface.tris[3,column]])
        for column in axes(surface.tris,2)]
    if protected_facets!==nothing
        for face in protected_facets
            push!(facets,ntuple(i->inverse[face[i]],3))
        end
    end
    segments=Set{Tuple{Int32,Int32}}()
    for face in facets,index in 1:3
        push!(segments,ModelNativeTetgenArena.edgekey(face[index],face[mod1(index+1,3)]))
    end
    arena=ModelNativeTetgenArena.Arena(kernel;facets,segments)
    ModelNativeTetgenArena.improve_quality32!(arena,work)>0 || return nothing
    ModelNativeTetgenArena.consistency(arena) || throw(ArgumentError(
        "$caller: native initial quality flips lost reciprocal face bonds"))
    cells,_=ModelNativeTetgenArena.real_cells(arena)
    mapped=Matrix{Int32}(undef,4,length(cells))
    for (column,cell) in enumerate(cells)
        vertices=ntuple(slot->Int32(order[cell[slot]]),4)
        sign=ModelNativeInitialTet.orient3(
            (_model_mesh_coordinate(surface,row) for row in vertices)...)
        sign>0 && (vertices=(vertices[2],vertices[1],vertices[3],vertices[4]))
        for slot in 1:4;mapped[slot,column]=vertices[slot];end
    end
    return Mesh(surface.coords;tets=mapped)
end

# A geometric fill can cover a planar Source sheet with the opposite
# diagonal. Actual lower constructor facets also need their exact row
# incidence before the fill can return ahead of constrained recovery.
function _model_native_boundary_facets_match(surface,mesh)
    size(surface.tris,1)==3 && size(mesh.tets,1)==4 || return false
    wanted=Set{NTuple{3,Int32}}()
    for column in axes(surface.tris,2)
        vertices=surface.tris[:,column]
        allunique(vertices) && all(row->0<row<=size(surface.coords,2),vertices) || return false
        face=Tuple(sort!(Int32.(vertices)))
        face in wanted && return false
        push!(wanted,face)
    end
    isempty(wanted) && return false
    counts=Dict{NTuple{3,Int32},Int}()
    for column in axes(mesh.tets,2)
        vertices=mesh.tets[:,column]
        allunique(vertices) && all(row->0<row<=size(mesh.coords,2),vertices) || return false
        for missing in 1:4
            face=Tuple(sort!(Int32[vertices[i] for i in 1:4 if i!=missing]))
            counts[face]=get(counts,face,0)+1
        end
    end
    all(count->count in (1,2),values(counts)) || return false
    return Set(face for (face,count) in counts if count==1)==wanted
end

function _model_native_initial_volume_fill(surface::Mesh,volume::Int,caller;
                                           interior_points=nothing,
                                           node_order=nothing,initial_surface=surface,options=nothing,births=nothing,
                                           protected_facets=nothing)
    extra=interior_points===nothing ? NTuple{3,Float64}[] : interior_points
    failures=String[]
    function certify(mesh)
        passed,why=_certify_surface_fill(surface,mesh)
        passed && node_order!==nothing && !_model_native_boundary_facets_match(surface,mesh) &&
            return false,"actual PLC constructor facet incidence differs"
        return passed,why
    end
    if node_order!==nothing
        native_kernel=Ref{ModelNativeInitialTet.Kernel}()
        mesh=_model_native_initial_delaunay_fill(initial_surface,node_order,options,caller;
            births,kernel_output=native_kernel)
        valid,reason=certify(mesh)
        if valid
            # Tetgen's -Y reconstruction optimizer runs even when the later
            # Gmsh optimizer is disabled. These flips precede final cell births.
            quality=_model_native_quality_volume_fill(initial_surface,node_order,native_kernel[],caller;protected_facets)
            if quality!==nothing
                valid,reason=certify(quality)
                valid || throw(ArgumentError(
                    "$caller: native initial quality fill broke its PLC: $reason"))
                mesh=quality
            end
            finished=_finish_interior(mesh,extra)
            finished.tets==mesh.tets && finished.coords==mesh.coords &&
                _model_native_volume_cell_births!(births,finished)
            return finished
        end
        push!(failures,reason)
        try
            recovered=_model_native_recovered_volume_fill(initial_surface,node_order,native_kernel[],caller;protected_facets)
            valid,reason=certify(recovered)
            if valid
                finished=_finish_interior(recovered,extra)
                finished.tets==recovered.tets && finished.coords==recovered.coords &&
                    _model_native_volume_cell_births!(births,finished)
                return finished
            end
            push!(failures,reason)
        catch error
            error isa ModelNativeTetgenArena.Unsupported || rethrow()
            push!(failures,sprint(showerror,error))
        end
    end
    try
        # Constrained recovery must retain every actual initial vertex,
        # including embedded inputs already present before the cube births.
        T=delaunay3d(vec(copy(initial_surface.coords[1,:])),vec(copy(initial_surface.coords[2,:])),
                     vec(copy(initial_surface.coords[3,:]));perturb=false,
                     preserve_input_rows=node_order!==nothing)
        selected=_classify_by_centroid(T,surface)
        for drop_flat in (false,true)
            keep=drop_flat ? _rb_drop_flats(T,selected) : selected
            mesh=to_mesh3(T;keep,preserve_input_rows=node_order!==nothing)
            valid,reason=certify(mesh)
            valid && return _finish_interior(mesh,extra)
            push!(failures,reason)
        end
    catch error
        error isa InterruptException && rethrow()
        (error isa ArgumentError || error isa ErrorException) || rethrow()
        push!(failures,sprint(showerror,error))
    end
    try
        # Float Steiner recovery is the cheaper general path and must run
        # before the exact-CDT fallback: the Rational{BigInt} kernel can burn
        # through its predicate budget (or the process heap) on benign
        # near-cospherical boundaries the float kernel certifies quickly.
        mesh=recover_boundary(surface;max_seeds=8,steiner=true)
        valid,reason=certify(mesh)
        valid && return _finish_interior(mesh,extra)
        push!(failures,reason)
    catch error
        error isa InterruptException && rethrow()
        (error isa ArgumentError || error isa ErrorException) || rethrow()
        push!(failures,sprint(showerror,error))
    end
    if applicable(_recover_boundary_exact,surface)
        try
            # The exact-CDT backend only returns a mesh after its own
            # `certify_exact` — a Rational{BigInt} audit strictly stronger than
            # the Float64 `_certify_surface_fill` gate above: same crease-chain,
            # pierce, boundary-membership and per-region-area predicates, all in
            # exact arithmetic. Re-running the float gate here can only produce
            # false negatives on its output: exact Steiner vertices split crease
            # and facet chains at rational parameters whose Float64 images land
            # ~1 ulp off the supporting line, so `_rb_onseg`'s exact-zero
            # collinearity probe finds no covering chain and reports "crease
            # segment not recovered" for a segment the exact certificate proved
            # conforming. `_model_native_boundary_facets_match` is likewise
            # unsatisfiable by construction for any Steiner-subdivided boundary —
            # its literal-incidence check guards against `oriented_surface_
            # constructor` diagonal flips, and the exact backend never consumes
            # the constructor mesh. `tetrahedralize` already trusts this same
            # certificate; match it.
            return _finish_interior(_recover_boundary_exact(surface),extra)
        catch error
            error isa InterruptException && rethrow()
            (error isa ArgumentError || error isa ErrorException) || rethrow()
            push!(failures,sprint(showerror,error))
        end
    end
    throw(ArgumentError("$caller: native Volume[$volume] boundary recovery failed — "*
                        join(failures,"; ")))
end
