# ── Built-in ruled `Surface`: the plane producer on the ruled chart ─────────
#
# `meshGenerator` meshes a non-plane built-in face exactly as it meshes a
# plane — PolyMesh initial mesh of the boundary parameters, recovered edges,
# BDS2GMSH, `bowyerWatsonFrontal` (or `bowyerWatson`), `laplaceSmoothing` —
# except that the parameters come from `reparamOnFace` on the ruled patch,
# `delaunayizeBDS` swaps the initial mesh in the first fundamental form, and
# the frontal points are put on the surface (`optimalPointFrontalB`).

# The initial mesh: PolyMesh of the source parameters, then the BDS swap
# pass. Returns the triangulation with the swaps applied, its interior and
# the initial cells in BDS2GMSH order with their MTriangle node order.
function _gmsh_native_ruled_initial(source::_ModelSurfaceSource,tag::Int,ground)
    _gmsh_native_surface_split_constraints!(source)
    rows=Mesh2D._gmsh_initial_uid_order(source.occurrence)
    order=Int[source.positions[label] for label in source.occurrence[rows]]
    product=Mesh2D._gmsh_poly_initial(source.x,source.y,source.segments,order;
        vertex_labels=1:length(source.tags),surface_tag=tag)
    T=Mesh2D._gmsh_poly_triangulation(product,source.boundaries,source.internals)
    inside=classify_interior(T)
    p=product.mesh;n=length(source.tags)
    row_of=zeros(Int32,length(p.x));row_of[product.source_rows]=Int32.(1:n)
    # `_gmsh_poly_triangulation` creates one slot per non-box face, in face order.
    interior_faces=falses(length(p.face_edge));slot=0
    for f in eachindex(p.face_edge)
        h=p.face_edge[f]
        (row_of[p.origin[h]]==0 || row_of[p.origin[p.next[h]]]==0 ||
         row_of[p.origin[p.next[p.next[h]]]]==0) && continue
        slot+=1
        interior_faces[f]=inside[slot]
    end
    # BDS point numbers: the rank of each node in `all_vertices` (ordered by
    # node number); the four box points are -1..-4.
    rank=zeros(Int,n);rank[sortperm(source.tags)]=0:n-1
    ids=Int[v<=4 ? -v : rank[row_of[v]] for v in eachindex(p.x)]
    constrained=Tuple{Int32,Int32}[(product.source_rows[a],product.source_rows[b])
                                   for (a,b,_) in source.segments]
    owners=Tuple{Int,Int}[v<=4 ? (2,0) : (Int(source.owners[row_of[v]][1]),
        Int(source.owners[row_of[v]][2])) for v in eachindex(p.x)]
    bds=Mesh2D._gmsh_bds_delaunayize(p,ids,interior_faces,constrained,ground;owners)
    for (a,b,c,d) in bds.swaps
        ra=row_of[a];rb=row_of[b];rc=row_of[c];rd=row_of[d]
        done=false
        for t in eachindex(T.alive)
            (T.alive[t] && inside[t] && !Mesh2D._is_ghost_tri(T,t)) || continue
            for k in 1:3
                q,r=Mesh2D._edge(T,t,k)
                ((q==ra && r==rb) || (q==rb && r==ra)) || continue
                apex=Mesh2D._vert(T,t,k)
                (apex==rc || apex==rd) || throw(ErrorException(
                    "native initial mesh: BDS swap does not match the triangulation"))
                Mesh2D._flip!(T,t,k)
                done=true
                break
            end
            done && break
        end
        done || throw(ErrorException("native initial mesh: BDS swap edge not found"))
    end
    slots=Dict{NTuple{3,Int32},Int32}()
    for t in eachindex(T.alive)
        (T.alive[t] && inside[t] && !Mesh2D._is_ghost_tri(T,t)) || continue
        slots[Mesh2D._gmsh_triangle_key(T,t)]=Int32(t)
    end
    initial=Tuple{Int32,NTuple{3,Int32}}[]
    for nodes in bds.faces
        mapped=(row_of[nodes[1]],row_of[nodes[2]],row_of[nodes[3]])
        push!(initial,(slots[Mesh2D._gmsh_sorted_key(mapped...)],mapped))
    end
    length(initial)==length(slots) || throw(ErrorException(
        "native initial mesh: BDS faces do not cover the interior"))
    return T,inside,initial
end

# `BGM_MeshSize` at the lifted point of `(u,v)`.
struct _GmshRuledSurfaceBackground{C,B}
    chart::C
    background::B
end
(background::_GmshRuledSurfaceBackground)(u,v)=
    background.background(_surface_metric_lift(background.chart,u,v)...)

# Whether the ruled producer below serves Surface `tag`: a warped built-in
# patch meshed by Delaunay or Frontal-Delaunay, without the recombination
# and periodic paths the plane producer has and this one has not been given.
function _gmsh_native_ruled_admitted(m::GeoModel,tag::Int,algorithm)
    algorithm in (5,6) || return false
    _model_surface_recombined(m,tag) && return false
    any(pair->first(pair)[1]==2 && (first(pair)[2]==tag ||
        abs(Int(last(pair).master_entity))==tag),m.periodic) && return false
    haskey(m.surfaces,tag) && length(m.surfaces[tag])==1 || return false
    length(m.loops[first(m.surfaces[tag])]) in (3,4) || return false
    return true
end

function _gmsh_native_ruled_surface(m::GeoModel,tag::Int,options,input,
                                    high_water,caller;algorithm=6,size_field=nothing,
                                    element_high_water=UInt64(0),forced=nothing,
                                    _periodic_positions=nothing)
    run()=_gmsh_native_ruled_surface_run(m,tag,options,input,high_water,
        caller,algorithm,size_field,element_high_water,forced,_periodic_positions)
    haskey(task_local_storage(),:_ruled_eval_cache) && return run()
    return task_local_storage(run,:_ruled_eval_cache,Dict{Tuple{Symbol,Int},Any}())
end

function _gmsh_native_ruled_surface_run(m::GeoModel,tag::Int,options,input,
                                        high_water,caller,algorithm,size_field,
                                        element_high_water,forced,_periodic_positions)
    authority,limit=_gmsh_native_surface_lower_input(m,options,input,high_water,caller)
    chart=_gmsh_ruled_chart(m,tag,caller)
    next_tag=Ref{UInt64}(limit)
    source=_gmsh_native_surface_source(m,tag,authority,chart,caller;
        ctx_lc=options.ctx_lc,forced,next_tag,_periodic_positions)
    limit=next_tag[]
    ground=chart.form
    T,inside,initial=_gmsh_native_ruled_initial(source,tag,ground)
    curve_chains=_gmsh_native_surface_curve_chains(source)
    numbers=source.tags
    trace=nothing
    if length(source.tags)!=3
        embedded_points,_=_model_surface_embedding_tags(m,tag,caller)
        sizes=_gmsh_native_surface_vertex_sizes(T,inside,source,m,embedded_points)
        background=_GmshRuledSurfaceBackground(chart,
            _ModelAlgorithmBackgroundSize(m,2,tag,options,size_field,caller))
        background_sizes=copy(sizes)
        model_vertices=Int32[]
        for point in embedded_points
            row=source.point_positions[point]
            background_sizes[row]=background.background(source.world[row]...)
            push!(model_vertices,Int32(row))
        end
        extend=_resolved_extend_from_boundary(m,2,tag)
        trace=Mesh2D._gmsh_surface_trace(T,inside,source.world,chart.lift,
            element_high_water,sizes,background_sizes,extend>0 || extend==-2;
            temporary_cells=count(inside),
            metric=Mesh2D._GmshFirstFundamental(ground),initial)
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
    product=_gmsh_native_surface_export(T,inside,numbers,source,chart,tag,limit;
        projected=true,trace)
    if trace===nothing
        # meshGenerator's three-vertex shortcut constructs one triangle.
        cells=ntris(product.mesh)
        births=UInt64[element_high_water+UInt64(column) for column in 1:cells]
        product=_ModelSurfaceProduction(product.mesh,product.node_tags,
            product.node_owners,product.node_parameters,product.retained_nodes,
            product.node_high_water,births,element_high_water+UInt64(cells),
            product.birth_rows)
    end
    mesh=product.mesh
    physical=Matrix{Float64}(undef,3,nnodes(mesh))
    for row in axes(physical,2)
        original=get(source.positions,product.node_tags[row],0)
        physical[:,row].=original>0 ? source.world[original] :
            _surface_metric_lift(chart,mesh.coords[1,row],mesh.coords[2,row])
    end
    if trace!==nothing
        product,physical=_gmsh_native_surface_smooth(m,tag,product,physical,chart,caller)
        mesh=product.mesh
    end
    get(m.meshing.reverse,(2,tag),false) &&
        (mesh=_reversed_surface_mesh(mesh,caller,tag))
    nnodes(mesh)==length(product.node_tags) || throw(ArgumentError(
        "$caller: Surface[$tag] postprocessing changed the actual producer node set"))
    for row in axes(mesh.coords,2)
        mesh.coords[1,row]=physical[1,row];mesh.coords[2,row]=physical[2,row]
        mesh.coords[3,row]=physical[3,row]
    end
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid native ruled Surface[$tag] — "*join(diagnostic.messages,"; ")))
    return _ModelSurfaceProduction(mesh,product.node_tags,product.node_owners,
        product.node_parameters,product.retained_nodes,product.node_high_water,
        product.element_tags,product.element_high_water,product.birth_rows,nothing,
        product.source_node_updates,product.native_recombine_events,curve_chains)
end
