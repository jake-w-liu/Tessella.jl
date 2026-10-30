# ════════════════════════════════════════════════════════════════════════════════
# Native adaptive remeshing — AdaptMesh (`.geo`) and CLI batch mode `-4`.
#
# Bounded Coupez-style loop over the *current* simplex carrier:
#   per iteration  1. collapse short interior metric edges
#                  2. swap tetrahedra whose quality is below the threshold
#                  3. split over-long metric edges via `refine_to_size`
# Entity ownership travels on `seg_tag`/`tri_tag`/`tet_tag` (per-family entity
# tags). Boundary triangles and segments are embedded in the swap pass, so
# boundary topology, interfaces and entity ownership cannot change; refinement
# splits coincident lower-dimensional supports so the boundary stays conforming.
# Vertex relocation is disabled: every surviving node keeps its coordinates, so
# callers can recover provenance (and thereby entity/node bindings) exactly.
# ════════════════════════════════════════════════════════════════════════════════

@inline _adapt_ek(a::Int32,b::Int32)=a<=b ? (a,b) : (b,a)
@inline _adapt_ek(a,b)=_adapt_ek(Int32(a),Int32(b))

@inline function _adapt_gamma(coords,q)
    a=(coords[1,q[1]],coords[2,q[1]],coords[3,q[1]])
    b=(coords[1,q[2]],coords[2,q[2]],coords[3,q[2]])
    c=(coords[1,q[3]],coords[2,q[3]],coords[3,q[3]])
    d=(coords[1,q[4]],coords[2,q[4]],coords[3,q[4]])
    value,_=_gopt_gamma(a,b,c,d)
    return value
end

@inline function _adapt_signed_volume(coords,q)
    a=(coords[1,q[1]],coords[2,q[1]],coords[3,q[1]])
    b=(coords[1,q[2]],coords[2,q[2]],coords[3,q[2]])
    c=(coords[1,q[3]],coords[2,q[3]],coords[3,q[3]])
    d=(coords[1,q[4]],coords[2,q[4]],coords[3,q[4]])
    return tet_signed_volume(a,b,c,d)
end

# Unique edges of the live tetrahedra, in deterministic node order.
function _adapt_tet_edges(tets::Vector{NTuple{4,Int32}},dead)
    seen=Set{Tuple{Int32,Int32}}()
    @inbounds for t in eachindex(tets)
        dead[t] && continue
        q=tets[t]
        for (a,b) in ((q[1],q[2]),(q[1],q[3]),(q[1],q[4]),(q[2],q[3]),(q[2],q[4]),(q[3],q[4]))
            push!(seen,_adapt_ek(a,b))
        end
    end
    return seen
end

"""
Interior edge collapse. Both endpoints must be referenced by tetrahedra only —
never by a boundary/interface segment or triangle — so boundary topology cannot
change. Every live tetrahedron incident to the removed endpoint must share one
region tag (upstream requires both endpoints to sit on volume entities of the
same region); the absorbed vertex keeps its position, cavity tets containing the
collapsed edge die, and all other cavity tets are re-pointed while keeping their
own tag. A collapse commits only when every re-pointed tetrahedron stays
strictly positively oriented and the post-collapse worst cavity quality is at
least `min(previous_worst, qmin)`. Returns the number of committed collapses.
"""
function _adapt_collapse_pass(m::Mesh,field,entities,ratio,qmin)
    tets=NTuple{4,Int32}[ntuple(k->m.tets[k,t],4) for t in axes(m.tets,2)]
    ttag=Int32[m.tet_tag[t] for t in axes(m.tets,2)]
    dead=falses(length(tets))
    boundary=Set{Int32}()
    @inbounds for s in axes(m.segs,2)
        push!(boundary,m.segs[1,s]);push!(boundary,m.segs[2,s])
    end
    @inbounds for f in axes(m.tris,2), k in 1:3
        push!(boundary,m.tris[k,f])
    end
    incidence=Dict{Int32,Vector{Int32}}()
    @inbounds for t in eachindex(tets), k in 1:4
        push!(get!(() -> Int32[],incidence,tets[t][k]),Int32(t))
    end
    # Deterministic candidate order: ascending metric length, then edge key.
    candidates=Tuple{Float64,Int32,Int32}[]
    for (a,b) in _adapt_tet_edges(tets,dead)
        (a in boundary || b in boundary) && continue
        p=(m.coords[1,a],m.coords[2,a],m.coords[3,a])
        q=(m.coords[1,b],m.coords[2,b],m.coords[3,b])
        best=0.0
        if entities===nothing
            best=metric_edge_length(field,p,q)
        else
            for t in get(incidence,b,Int32[])
                dead[t] && continue
                context=get(entities,ttag[t],nothing)
                best=max(best,metric_edge_length(field,p,q;entity=context))
            end
        end
        isfinite(best) || continue
        best<ratio || continue
        # A single region tag must cover the whole removed star.
        ts=[t for t in get(incidence,b,Int32[]) if !dead[t]]
        isempty(ts) && continue
        tag=ttag[ts[1]]
        all(t->ttag[t]==tag,ts) || continue
        push!(candidates,(best,a,b))
    end
    sort!(candidates)
    collapsed=0
    for (_,a,b) in candidates
        # `b` folds into `a`. Survivor choice is deterministic (smaller index
        # absorbs) so repeated runs reproduce.
        a<b || ((a,b)=(b,a))
        haskey(incidence,a) || continue            # `a` was itself collapsed away
        ts=[t for t in get(incidence,b,Int32[]) if !dead[t] && b in tets[t]]
        isempty(ts) && continue
        tag=ttag[ts[1]]
        all(t->ttag[t]==tag,ts) || continue
        worst_old=Inf;survive=NTuple{4,Int32}[];dies=Int32[]
        ok=true
        for t in ts
            q=tets[t]
            worst_old=min(worst_old,_adapt_gamma(m.coords,q))
            if a in q
                push!(dies,t)                     # degenerates to a face — dies
            else
                nq=ntuple(k->q[k]==b ? a : q[k],4)
                _adapt_signed_volume(m.coords,nq)>0 || (ok=false;break)
                worst=_adapt_gamma(m.coords,nq)
                worst+1e-12<min(worst_old,qmin) && (ok=false;break)
                push!(survive,nq)
            end
        end
        ok || continue
        length(dies)<length(ts) || continue       # whole star dying leaves a hole
        si=1
        for t in ts
            if t in dies || (a in tets[t] && b in tets[t])
                dead[t]=true
            else
                tets[t]=survive[si];si+=1
                push!(incidence[a],t)
            end
        end
        delete!(incidence,b)
        collapsed+=1
    end
    collapsed==0 && return m,0
    keep=[t for t in eachindex(tets) if !dead[t]]
    M=Matrix{Int32}(undef,4,length(keep));tags=Vector{Int32}(undef,length(keep))
    @inbounds for (j,t) in enumerate(keep)
        q=tets[t];M[1,j],M[2,j],M[3,j],M[4,j]=q[1],q[2],q[3],q[4];tags[j]=ttag[t]
    end
    out=Mesh(m.coords;segs=m.segs,tris=m.tris,tets=M,
        seg_tag=m.seg_tag,tri_tag=m.tri_tag,tet_tag=tags)
    diag=validate(out)
    diag.ok || return m,0                      # stay safe: never emit an invalid mesh
    return out,collapsed
end

# Remove nodes that no cell references (collapse orphans). Segments/triangles
# cannot reference collapsed interior nodes, but dead cells of any family are
# dropped and node indices compacted for a clean downstream bind.
function _adapt_compact(m::Mesh)
    used=falses(size(m.coords,2))
    @inbounds for t in axes(m.tets,2), k in 1:4;used[m.tets[k,t]]=true;end
    @inbounds for f in axes(m.tris,2), k in 1:3;used[m.tris[k,f]]=true;end
    @inbounds for s in axes(m.segs,2), k in 1:2;used[m.segs[k,s]]=true;end
    all(used) && return m
    nid=zeros(Int32,length(used));nv=0
    @inbounds for i in eachindex(used);used[i] && (nid[i]=(nv+=1));end
    coords=Matrix{Float64}(undef,3,nv)
    @inbounds for i in eachindex(used)
        used[i] || continue
        coords[1,nid[i]]=m.coords[1,i];coords[2,nid[i]]=m.coords[2,i];coords[3,nid[i]]=m.coords[3,i]
    end
    remap(cells)=map(v->nid[v],cells)
    Mesh(coords;segs=remap(m.segs),tris=remap(m.tris),tets=remap(m.tets),
        seg_tag=m.seg_tag,tri_tag=m.tri_tag,tet_tag=m.tet_tag)
end

@inline _adapt_face_key(a,b,c)=Tuple(sort!(Int32[a,b,c]))

"""
`adapt_mesh(mesh, field; ...) -> (mesh, stats)`

Coupez-style adaptive remeshing of the current simplex mesh against a size /
metric field. The mesh's `seg_tag`/`tri_tag`/`tet_tag` carry per-dimension
entity tags: they are propagated through splits, retained by swaps and collapse,
and select per-entity field contexts when `entity_resolver` is supplied (same
resolver contract as [`refine_to_size`](@ref)). The loop is bounded by `niter`
and terminates early once a full pass makes no change.

Guarantees: boundary segments/triangles stay conforming (refinement splits the
coincident supports, collapse never touches a boundary-referenced vertex, and
all boundary cells are embedded in the swap pass); every emitted tetrahedron is
positively oriented (`validate`d); tags are conservative (children inherit the
parent tag). With `strict=false` (default) edges that cannot be split to
Float64 resolution remain unsplit and are counted in `stats.unresolved_edges`
rather than throwing; `strict=true` propagates the `refine_to_size` diagnostic.

`collapse_ratio` collapses interior edges whose metric length is below it;
`quality_threshold` drives the swap pass (upstream `qMin`-style gamma);
`collapse_quality_floor` bounds how much a collapse may degrade the local worst
quality (`sliverLimit`-style).
"""
function adapt_mesh(m::Mesh,field::AbstractSizeField;
        entity_resolver=nothing,vertex_entities=nothing,
        niter::Integer=10,collapse_ratio::Real=0.4,
        quality_threshold::Real=0.3,collapse_quality_floor::Real=0.2,
        strict::Bool=false,max_nodes=PIPELINE_DEFAULT_MAX_NODES,
        max_tets=PIPELINE_DEFAULT_MAX_TETS,
        caller::AbstractString="adapt_mesh")
    m isa Mesh || throw(ArgumentError(
        "$caller: adaptive remeshing requires a simplex carrier mesh"))
    size(m.tets,2)>0 || throw(ArgumentError(
        "$caller: adaptive remeshing currently supports tetrahedral meshes " *
        "(the input carries no tetrahedra)"))
    0<collapse_ratio<1 || throw(ArgumentError(
        "$caller: collapse_ratio must lie in (0,1), got $collapse_ratio"))
    niter>=0 || throw(ArgumentError("$caller: niter must be non-negative"))
    diag=validate(m)
    diag.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(diag.messages,"; ")))
    entities=_mesh3_tag_entities(m.tet_tag,entity_resolver,3,caller;required=true)
    embedded_faces=Set{NTuple{3,Int32}}()
    for f in axes(m.tris,2)
        push!(embedded_faces,_adapt_face_key(m.tris[1,f],m.tris[2,f],m.tris[3,f]))
    end
    embedded_edges=Set{Tuple{Int32,Int32}}()
    for s in axes(m.segs,2)
        push!(embedded_edges,_adapt_ek(m.segs[1,s],m.segs[2,s]))
    end
    mesh=m;splits=0;collapses=0;swaps=0;iterations=0;refined_stable=false
    for iteration in 1:Int(niter)
        changed=false;iterations=iteration
        # 1. collapse short interior metric edges
        mesh,ncollapse=_adapt_collapse_pass(mesh,field,entities,
            Float64(collapse_ratio),Float64(collapse_quality_floor))
        collapses+=ncollapse;changed|=ncollapse>0
        mesh=_adapt_compact(mesh)
        # 2. quality-driven edge swaps (boundary cells embedded/frozen)
        result=_gmsh_optimize_tetrahedra(mesh;
            embedded_faces=embedded_faces,embedded_edges=embedded_edges,
            relocate=false,threshold=Float64(quality_threshold),positive=true)
        nswap=result.edge_swaps
        nswap>0 && (mesh=result.mesh;swaps+=nswap;changed=true)
        # 3. split over-long metric edges. The field-evaluation pass is the
        # expensive step; when the previous refine was already a no-op and this
        # round's collapse/swap changed nothing, the mesh is unchanged since
        # that refine — another one cannot move anything either.
        if refined_stable && !changed
            break
        end
        refined=refine_to_size(mesh,field;
            entity_resolver=entity_resolver,vertex_entities=vertex_entities,
            max_nodes,max_tets,best_effort=!strict)
        moved=(size(refined.coords,2)!=size(mesh.coords,2) ||
               size(refined.tets,2)!=size(mesh.tets,2))
        moved && (splits+=1;changed=true)
        refined_stable=!moved
        mesh=refined
        changed || break
    end
    # Certify the metric postcondition on the emitted mesh. `refine_to_size`
    # guarantees it under `strict`; in best-effort mode residual violators are
    # counted and reported.
    unresolved=0
    if !strict
        dead=falses(size(mesh.tets,2))
        tets=NTuple{4,Int32}[ntuple(k->mesh.tets[k,t],4) for t in axes(mesh.tets,2)]
        for e in _adapt_tet_edges(tets,dead)
            p=(mesh.coords[1,e[1]],mesh.coords[2,e[1]],mesh.coords[3,e[1]])
            q=(mesh.coords[1,e[2]],mesh.coords[2,e[2]],mesh.coords[3,e[2]])
            metric_edge_length(field,p,q)>1+1e-9 && (unresolved+=1)
        end
    end
    diag=validate(mesh)
    diag.ok || throw(ErrorException(
        "$caller: adaptive remeshing produced an invalid mesh — " *
        join(diag.messages,"; ")))
    return mesh,(iterations=iterations,splits=splits,collapses=collapses,
        swaps=swaps,unresolved_edges=unresolved)
end
