# Native insertVerticesInRegion marks actual Surface Tri constructors for
# sizing. The temporary Quad fan is an independent recovery constraint.
module ModelNativeHybridInsertion
const M=parentmodule(@__MODULE__)
function quiescent(mesh,node_uids,element_uids,triangles,caller;minimize=false)
    length(node_uids)==M.nnodes(mesh) && allunique(node_uids) && all(!iszero,node_uids) &&
        length(element_uids)==M.ntets(mesh) && allunique(element_uids) && all(!iszero,element_uids) &&
        size(mesh.tets,1)==4 && all(isfinite,mesh.coords) &&
        all(row->0<row<=M.nnodes(mesh),mesh.tets) || throw(ArgumentError(
        "$caller: hybrid insertion lacks complete current UID/XYZ/cell constructor incidence"))
    sizes=fill(minimize ? Inf : 0.0,M.nnodes(mesh))
    boundary=falses(M.nnodes(mesh));seen=Set{UInt64}()
    function update(first,last)
        a=ntuple(axis->mesh.coords[axis,first],3);b=ntuple(axis->mesh.coords[axis,last],3)
        x=a[1]-b[1];y=a[2]-b[2];z=a[3]-b[3]
        length=sqrt(x*x+y*y+z*z)
        for row in (first,last)
            sizes[row]=minimize ? min(sizes[row],length) : max(sizes[row],length)
        end
    end
    for triangle in triangles
        triangle.owner[1]==2 && triangle.owner[2]>0 && triangle.uid>0 &&
            !(triangle.uid in seen) && length(triangle.rows)==3 && allunique(triangle.rows) &&
            all(row->0<row<=M.nnodes(mesh),triangle.rows) &&
            Tuple(node_uids[collect(triangle.rows)])==triangle.nodes || throw(ArgumentError(
            "$caller: hybrid sizing Tri differs from its actual owned constructor receipt"))
        push!(seen,triangle.uid)
        for row in triangle.rows;boundary[row]=true;end
        for side in 1:3
            update(triangle.rows[side],triangle.rows[mod1(side+1,3)])
        end
    end
    # Native MTet setLcs sees every current Tet edge for vertices absent from
    # the actual Tri marker set. Quad-only inputs therefore retain these sizes.
    for cell in eachcol(mesh.tets),first in 1:4,last in first+1:4
        a=cell[first];b=cell[last]
        p=ntuple(axis->mesh.coords[axis,a],3);q=ntuple(axis->mesh.coords[axis,b],3)
        x=p[1]-q[1];y=p[2]-q[2];z=p[3]-q[3]
        length=sqrt(x*x+y*y+z*z)
        for row in (a,b)
            boundary[row] || (sizes[row]=minimize ? min(sizes[row],length) : max(sizes[row],length))
        end
    end
    all(value->isfinite(value) && value>0.0,sizes) || throw(ArgumentError(
        "$caller: hybrid insertion has no native incident size for a constructor row"))
    tets=copy(mesh.tets);radii=Float64[]
    for column in axes(tets,2)
        points0=ntuple(slot->ntuple(axis->mesh.coords[axis,tets[slot,column]],3),4)
        flip=M.orient3(points0...)>0
        if flip
            tets[1,column],tets[2,column]=tets[2,column],tets[1,column]
        end
        cell=tets[:,column]
        points=flip ? ntuple(slot->ntuple(axis->mesh.coords[axis,cell[slot]],3),4) : points0
        center=M._gmsh3_native_center(points...)
        x=points[1][1]-center[1];y=points[1][2]-center[2];z=points[1][3]-center[3]
        length=0.25*(sizes[cell[1]]+sizes[cell[2]]+sizes[cell[3]]+sizes[cell[4]])
        radius=sqrt(x*x+y*y+z*z)/length
        isfinite(radius) && 0<=radius<1.0 || throw(M.ModelNativeTetgenArena.Unsupported(
            "hybrid normalized-radius constructor requires the native insertion branch"))
        push!(radii,radius)
    end
    permutation=sortperm(eachindex(radii);lt=(a,b)->radii[a]>radii[b] ? true :
        radii[a]<radii[b] ? false : element_uids[a]<element_uids[b])
    return (mesh=M.Mesh(mesh.coords;tets=tets[:,permutation]),
        element_numbers=UInt64.(element_uids[permutation]),node_numbers=UInt64.(node_uids),
        radii=radii[permutation],sizes=sizes)
end

# The scoped driver consumes the actual GFace marker sizes. The generic
# exterior-facet sizing policy and its Arena remain unchanged.
const D=parentmodule(M._refine_gmsh_delaunay3)
const NATIVE_FACES=((1,2,3),(1,3,4),(1,4,2),(2,4,3))

function native_face(T,cell,role,storage)
    rows=ntuple(slot->storage[cell][NATIVE_FACES[role][slot]],3)
    slot=D._nslot_by_face(T,cell,rows...)
    slot>0 || throw(ArgumentError("hybrid cavity constructor face has no current carrier"))
    return (cell=Int32(cell),role=role,rows=rows,neighbor=D._nbr(T,cell,slot))
end

# The MTet constructor's four face slots and breadth-first cavity discovery
# determine new cell births. Triangulation storage has a different local
# permutation; it remains an independent topology carrier.
function native_cavity!(T,seed,point,tags,protected_faces,storage;
                        metric=nothing,entityof=t->nothing)
    previous=length(T.mark);resize!(T.mark,length(T.alive))
    previous<length(T.mark) && fill!(@view(T.mark[previous+1:end]),Int32(0))
    if T.epoch==typemax(Int32);fill!(T.mark,Int32(0));T.epoch=Int32(0);end
    T.epoch+=Int32(1)
    empty!(T.cavity);push!(T.cavity,Int32(seed));T.mark[seed]=T.epoch
    shell=NamedTuple[];cursor=1
    while cursor<=length(T.cavity)
        cell=T.cavity[cursor];cursor+=1
        for role in 1:4
            face=native_face(T,cell,role,storage);neighbor=face.neighbor
            if neighbor==0 || D._sort3t(face.rows...) in protected_faces
                push!(shell,face)
            elseif T.mark[neighbor]!=T.epoch
                if tags[neighbor]==tags[cell] &&
                        native_insphere(T,storage[neighbor],point,metric,
                            entityof(neighbor))
                    T.mark[neighbor]=T.epoch;push!(T.cavity,neighbor)
                else
                    push!(shell,face)
                end
            end
        end
    end
    return shell
end

function native_star!(T,point,tags,protected_faces,storage,shell)
    while true
        wrong=[face for face in shell if M.orient3((D._pt(T,row) for row in face.rows)...,point)>=0]
        isempty(wrong) && break
        for face in wrong
            key=D._sort3t(face.rows...)
            any(entry->D._sort3t(entry.rows...)==key,shell) || continue
            neighbor=face.neighbor
            neighbor>0 && tags[neighbor]==tags[face.cell] &&
                !(key in protected_faces) || return false
            T.mark[neighbor]!=T.epoch || throw(ArgumentError(
                "hybrid cavity shell names an already selected interior carrier"))
            for role in 1:4
                replacement=native_face(T,neighbor,role,storage)
                at=findfirst(entry->D._sort3t(entry.rows...)==D._sort3t(replacement.rows...),shell)
                at===nothing ? push!(shell,replacement) : deleteat!(shell,at)
            end
            T.mark[neighbor]=T.epoch;push!(T.cavity,neighbor)
        end
    end
    empty!(T.boundary)
    for face in shell;push!(T.boundary,(face.rows...,face.neighbor));end
    return true
end

function native_weights(T,rows,point)
    a,b,c,d=ntuple(slot->D._pt(T,rows[slot]),4)
    m00,m10,m20=ntuple(axis->b[axis]-a[axis],3)
    m01,m11,m21=ntuple(axis->c[axis]-a[axis],3)
    m02,m12,m22=ntuple(axis->d[axis]-a[axis],3)
    x,y,z=ntuple(axis->point[axis]-a[axis],3)
    determinant=m00*(m11*m22-m12*m21)-m01*(m10*m22-m12*m20)+m02*(m10*m21-m11*m20)
    determinant!=0 || return (1.0,0.0,0.0,0.0)
    inverse=1.0/determinant
    u=(x*(m11*m22-m12*m21)-m01*(y*m22-m12*z)+m02*(y*m21-m11*z))*inverse
    v=(m00*(y*m22-m12*z)-x*(m10*m22-m12*m20)+m02*(m10*z-y*m20))*inverse
    w=(m00*(m11*z-y*m21)-m01*(m10*z-y*m20)+x*(m10*m21-m11*m20))*inverse
    return (1.0-u-v-w,u,v,w)
end

function native_containing(T,point,storage;reference_tolerance=1e-6)
    for cell in T.cavity
        weights=native_weights(T,storage[cell],point)
        _,u,v,w=weights
        (u < -reference_tolerance || v < -reference_tolerance || w < -reference_tolerance ||
            u > ((1.0+reference_tolerance)-v-w)) && continue
        return cell,weights
    end
    return nothing
end

function positive_constructor(mesh)
    tets=copy(mesh.tets)
    for column in axes(tets,2)
        points=ntuple(slot->ntuple(axis->mesh.coords[axis,tets[slot,column]],3),4)
        if M.orient3(points...)>0
            tets[1,column],tets[2,column]=tets[2,column],tets[1,column]
        end
    end
    return M.Mesh(copy(mesh.coords);tets=tets,tet_tag=copy(mesh.tet_tag))
end

function native_sizes(mesh,node_uids,triangles,caller;minimize=false,seeds=nothing)
    length(node_uids)==M.nnodes(mesh) && allunique(node_uids) &&
        all(!iszero,node_uids) || throw(ArgumentError(
        "$caller: hybrid sizing lacks unique current native node identities"))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(node_uids))
    sizes=fill(minimize ? Inf : 0.0,M.nnodes(mesh));boundary=falses(M.nnodes(mesh))
    seen=Set{UInt64}()
    if seeds!==nothing
        for seed in M._model_native_hybrid_seed_records(seeds)
            row=get(positions,seed.uid,Int32(0));row==0 && continue
            all(axis->isequal(mesh.coords[axis,row],seed.coordinate[axis]),1:3) ||
                throw(ArgumentError("$caller: hybrid embedded sizing changes a selected actual UID coordinate"))
            sizes[row]=sizes[row]==0.0 ? seed.size : min(sizes[row],seed.size)
        end
    end
    for triangle in triangles
        triangle.owner[1]==2 && triangle.owner[2]>0 && triangle.uid>0 &&
            !(triangle.uid in seen) && length(triangle.nodes)==3 && allunique(triangle.nodes) ||
            throw(ArgumentError("$caller: global GFace sizing lacks unique actual Tri incidence"))
        push!(seen,triangle.uid)
        points=hasproperty(triangle,:coordinates) ? triangle.coordinates :
            ntuple(slot->ntuple(axis->mesh.coords[axis,triangle.rows[slot]],3),3)
        length(points)==3 && all(point->length(point)==3 && all(isfinite,point),points) ||
            throw(ArgumentError("$caller: actual GFace sizing coordinates do not align"))
        rows=ntuple(slot->get(positions,triangle.nodes[slot],Int32(0)),3)
        for slot in 1:3
            row=rows[slot];row==0 && continue
            all(axis->isequal(mesh.coords[axis,row],points[slot][axis]),1:3) ||
                throw(ArgumentError("$caller: global GFace sizing changes a selected actual node"))
            boundary[row]=true
        end
        for slot in 1:3
            next=mod1(slot+1,3);a=points[slot];b=points[next]
            x=a[1]-b[1];y=a[2]-b[2];z=a[3]-b[3]
            length=sqrt(x*x+y*y+z*z)
            for row in (rows[slot],rows[next])
                row==0 && continue
                sizes[row]=minimize ? min(sizes[row],length) : max(sizes[row],length)
            end
        end
    end
    for cell in eachcol(mesh.tets),first in 1:4,last in first+1:4
        a=cell[first];b=cell[last]
        x=mesh.coords[1,a]-mesh.coords[1,b]
        y=mesh.coords[2,a]-mesh.coords[2,b]
        z=mesh.coords[3,a]-mesh.coords[3,b]
        length=sqrt(x*x+y*y+z*z)
        for row in (a,b)
            boundary[row] || (sizes[row]=minimize ? min(sizes[row],length) : max(sizes[row],length))
        end
    end
    all(value->isfinite(value) && value>0.0,sizes) || throw(ArgumentError(
        "$caller: native hybrid constructor lacks a current incident size"))
    return sizes
end

# ── Tensor arm: metric-composed cavity/admission ─────────────────────────
# The 3-D counterpart of the 2-D `_GmshSurfaceMetricCarrier` contract: when
# `refine` receives a resolved `AbstractAnisoField`, the queue radius reads
# the metric-normalized circumradius `sqrt(RᵀMR)`, the candidate point is
# the metric circumcentre, cavity growth uses the metric in-sphere test,
# and the `.05` close-point admission reads `sqrt(dᵀMd)` in units of the
# floored local tensor. Scalar vertex `sizes`/`background`, the orient3
# star-shapedness certificate, barycentric containment, weights and all
# element/node bookkeeping keep their Euclidean semantics — the same
# split the 2-D arm keeps between `vsize`/`vbgm` scalars and the carrier
# tensor. Every metric solve that cannot produce a finite positive result
# falls back to the established Euclidean predicate so certificates stay
# total.
const SF=parentmodule(M.metric_at)

@inline _native_metric_quad(m::M.Metric3,x,y,z)=
    m.m11*x*x+m.m22*y*y+m.m33*z*z+2*(m.m12*x*y+m.m13*x*z+m.m23*y*z)

# `entity===nothing` resolves through the four-argument `metric_at`
# overload; a `(3,tag)` entity selects the Region-restricted evaluation
# exactly like the scalar region background's `background_size(x,y,z,tag)`.
# A field sample that cannot yield a finite SPD tensor — a MathEvalAniso
# coordinate singularity or a non-positive-definite evaluation — returns
# `nothing` so every consumer takes its established Euclidean fallback,
# the same totality contract `native_metric_center` gives a singular solve.
@inline function _native_metric_at(metric::M.AbstractAnisoField,x,y,z,entity)
    all(isfinite,(x,y,z)) || return nothing
    value=try
        M.metric_at(metric,x,y,z,entity)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    value isa M.Metric3 || return nothing
    return value
end

# The 3-D `_gmsh_metric_floor` analogue: a scalar local bound `h` folds
# into a resolved SPD metric through the canonical metric intersection —
# every eigendirection whose directional size exceeds `h` takes the
# eigenvalue floor `1/h²` (`λᵢ′ = max(λᵢ, 1/h²)`), the metric-space
# `min(lc_bgm, lc_vsize)` of the scalar arm's `extend_boundary`
# composition. Isotropic inputs therefore keep an isotropic result, and
# stronger eigendirections keep the field's anisotropy exactly like the
# 2-D floor.
function native_metric_floor(metric::M.Metric3,h::Float64)
    (isfinite(h) && h>0) || return metric
    floor_=inv(h*h)
    isfinite(floor_) || return metric
    values,vectors=SF._sym3_eigh(metric.m11,metric.m22,metric.m33,
        metric.m12,metric.m13,metric.m23)
    minimum(values)>=floor_ && return metric
    return SF.metric_from_axes(max(values[1],floor_),max(values[2],floor_),
        max(values[3],floor_),vectors...)
end

# The metric circumcentre of a tet under M: the unique point equidistant
# from all four vertices in the M-norm —
#     2(pᵢ−p₁)ᵀ M c = pᵢᵀ M pᵢ − p₁ᵀ M p₁,  i ∈ {2,3,4}.
# Returns `(center, r2)` with `r2 = (c−p₁)ᵀM(c−p₁)`, or `nothing` when the
# 3×3 system is singular — callers then keep the Euclidean predicate,
# mirroring the scalar arm's exact-arithmetic fallbacks.
function native_metric_center(pa,pb,pc,pd,m::M.Metric3)
    u=(pb[1]-pa[1],pb[2]-pa[2],pb[3]-pa[3])
    v=(pc[1]-pa[1],pc[2]-pa[2],pc[3]-pa[3])
    w=(pd[1]-pa[1],pd[2]-pa[2],pd[3]-pa[3])
    row(d)=(2*(d[1]*m.m11+d[2]*m.m12+d[3]*m.m13),
            2*(d[1]*m.m12+d[2]*m.m22+d[3]*m.m23),
            2*(d[1]*m.m13+d[2]*m.m23+d[3]*m.m33))
    a1=row(u);a2=row(v);a3=row(w)
    rhs=(_native_metric_quad(m,pb[1],pb[2],pb[3])-_native_metric_quad(m,pa[1],pa[2],pa[3]),
         _native_metric_quad(m,pc[1],pc[2],pc[3])-_native_metric_quad(m,pa[1],pa[2],pa[3]),
         _native_metric_quad(m,pd[1],pd[2],pd[3])-_native_metric_quad(m,pa[1],pa[2],pa[3]))
    det=a1[1]*(a2[2]*a3[3]-a2[3]*a3[2])-
        a1[2]*(a2[1]*a3[3]-a2[3]*a3[1])+
        a1[3]*(a2[1]*a3[2]-a2[2]*a3[1])
    bound=(abs(a1[1])+abs(a1[2])+abs(a1[3]))*
          (abs(a2[1])+abs(a2[2])+abs(a2[3]))*
          (abs(a3[1])+abs(a3[2])+abs(a3[3]))
    (isfinite(det) && bound>0 && abs(det)>64eps(Float64)*bound) || return nothing
    center=((rhs[1]*(a2[2]*a3[3]-a2[3]*a3[2])-
             a1[2]*(rhs[2]*a3[3]-a2[3]*rhs[3])+
             a1[3]*(rhs[2]*a3[2]-a2[2]*rhs[3]))/det,
            (a1[1]*(rhs[2]*a3[3]-a2[3]*rhs[3])-
             rhs[1]*(a2[1]*a3[3]-a2[3]*a3[1])+
             a1[3]*(a2[1]*rhs[3]-rhs[2]*a3[1]))/det,
            (a1[1]*(a2[2]*rhs[3]-rhs[2]*a3[2])-
             a1[2]*(a2[1]*rhs[3]-rhs[2]*a3[1])+
             rhs[1]*(a2[1]*a3[2]-a2[2]*a3[1]))/det)
    all(isfinite,center) || return nothing
    return center,_native_metric_quad(m,center[1]-pa[1],center[2]-pa[2],
        center[3]-pa[3])
end

# Metric-normalized circumradius for one cell — the queue-radius analogue
# of the 2-D carrier's `sqrt(RᵀMR)`. The resolved metric is taken at the
# cell centroid (the 2-D arm evaluates its callable at the triangle
# centroid) and is floored by the mean vertex `sizes` when the boundary
# extension applies — the metric-space `min(lc_bgm, mean vsize)` of the
# scalar `lc=min(bgm, sizes)` composition.
function native_metric_radius(T,vertices,sizes,background,extend,metric,entity)
    points=ntuple(slot->D._pt(T,vertices[slot]),4)
    m=_native_metric_at(metric,.25*(points[1][1]+points[2][1]+points[3][1]+points[4][1]),
        .25*(points[1][2]+points[2][2]+points[3][2]+points[4][2]),
        .25*(points[1][3]+points[2][3]+points[3][3]+points[4][3]),entity)
    m===nothing &&
        return D._gmsh3_native_radius(T,vertices,sizes,background,extend)
    # Mirror the scalar `lc=min(bgm,sizes)` composition: the metric-space
    # floor folds the mean vertex background (edge-derived for seeds and
    # the field evaluation for Steiner vertices — the same `vSizesBGM`
    # transport the scalar arm keeps), then the extend propagation
    # tightens it by the mean vertex size.
    h=.25*(background[vertices[1]]+background[vertices[2]]+
           background[vertices[3]]+background[vertices[4]])
    extend && (h=min(h,.25*(sizes[vertices[1]]+sizes[vertices[2]]+
                            sizes[vertices[3]]+sizes[vertices[4]])))
    m=native_metric_floor(m,h)
    solved=native_metric_center(points...,m)
    solved===nothing &&
        return D._gmsh3_native_radius(T,vertices,sizes,background,extend)
    return sqrt(max(solved[2],0.0))
end

function native_queue_tet!(queue,T,t,sizes,background,extend,skipped,numbers,
                           storage,metric,entity)
    T.alive[t] || return nothing
    key=D._gmsh3_key(T,t)
    (key in skipped || key in queue.pending) && return nothing
    normalized=native_metric_radius(T,
        storage===nothing ? ntuple(slot->D._vert(T,t,slot),4) : storage[t],
        sizes,background,extend,metric,entity)
    isfinite(normalized) || throw(ArgumentError(
        "native hybrid metric insertion: normalized circumradius is not finite"))
    normalized>=1 || return nothing
    push!(queue.pending,key)
    D._gmsh3_push!(queue,D._GmshTetEntry(normalized,key,Int32(t),numbers[t]))
    return nothing
end

# Metric in-sphere cavity test: the candidate point is inside the
# neighbour's metric circumball when its squared M-distance to the metric
# circumcentre is strictly below the metric circumradius squared — the
# tensor restatement of the scalar insphere·orient3 sign test (a
# co-spherical point stays outside, matching the strict comparison).
@inline function native_insphere(T,vertices,point,metric,entity)
    points=ntuple(slot->D._pt(T,vertices[slot]),4)
    metric===nothing &&
        return D.insphere(points...,point)*M.orient3(points...)>0
    m=_native_metric_at(metric,.25*(points[1][1]+points[2][1]+points[3][1]+points[4][1]),
        .25*(points[1][2]+points[2][2]+points[3][2]+points[4][2]),
        .25*(points[1][3]+points[2][3]+points[3][3]+points[4][3]),entity)
    m===nothing && return D.insphere(points...,point)*M.orient3(points...)>0
    solved=native_metric_center(points...,m)
    solved===nothing &&
        return D.insphere(points...,point)*M.orient3(points...)>0
    center,r2=solved
    return _native_metric_quad(m,point[1]-center[1],point[2]-center[2],
        point[3]-center[3])<r2
end

# The candidate insertion point under the metric arm — the tet's metric
# circumcentre under the unfloored resolved metric, the 3-D counterpart of
# `_gmsh_surface_delaunay_point` for the metric carrier. A singular solve
# keeps the Euclidean centre.
function native_metric_point(T,vertices,metric,entity,fallback)
    points=ntuple(slot->D._pt(T,vertices[slot]),4)
    m=_native_metric_at(metric,.25*(points[1][1]+points[2][1]+points[3][1]+points[4][1]),
        .25*(points[1][2]+points[2][2]+points[3][2]+points[4][2]),
        .25*(points[1][3]+points[2][3]+points[3][3]+points[4][3]),entity)
    m===nothing && return fallback
    solved=native_metric_center(points...,m)
    solved===nothing && return fallback
    return solved[1]
end

# The `.05·lc` close-point admission in normalized metric units: each
# shell vertex's displacement to the candidate is read under the resolved
# metric at the pair midpoint, floored by the candidate's effective scalar
# bound `lc` — the anisotropic restatement of the scalar arm's
# min-direction distance gate. `lc` already folds the propagated boundary
# extension and the scalarized background exactly like the scalar arm.
function native_metric_admission(T,point,metric,lc,entity)
    for (a,b,c,_) in T.boundary,node in (a,b,c)
        p=D._pt(T,node)
        resolved=_native_metric_at(metric,.5*(p[1]+point[1]),.5*(p[2]+point[2]),
            .5*(p[3]+point[3]),entity)
        if resolved===nothing
            # Unresolvable metric sample — the scalar arm's `.05·lc` gate.
            (isfinite(lc) && lc>0) || continue
            if hypot(point[1]-p[1],point[2]-p[2],point[3]-p[3])<.05*lc
                return false
            end
            continue
        end
        m=native_metric_floor(resolved,lc)
        M._metric_displacement(m,point[1]-p[1],point[2]-p[2],point[3]-p[3],
            "native hybrid metric insertion")<.05 && return false
    end
    return true
end

function refine(mesh,node_uids,element_uids,triangles,background_size,caller;
                minimize=false,extend=true,node_high_water=nothing,element_high_water=nothing,
                region_background=false,seeds=nothing,protected_edges=Set{NTuple{2,Int32}}(),
                protected_faces=Set{NTuple{3,Int32}}(),metric=nothing,metric_entity=nothing)
    owned=positive_constructor(mesh)
    sizes=native_sizes(owned,node_uids,triangles,caller;minimize,seeds)
    return refine_owned(owned,background_size,sizes,caller;
        extend_boundary=extend,minimize_boundary=minimize,
        node_numbers=node_uids,node_high_water,element_numbers=element_uids,element_high_water,
        native_storage_order=true,return_numbers=true,region_background,protected_edges,
        protected_faces,metric,metric_entity)
end

function smooth(product,initial_rows,passes,caller)
    initial_rows<=M.nnodes(product.mesh) && passes>=0 || throw(ArgumentError(
        "$caller: native hybrid smoothing lacks its insertion birth boundary"))
    movement=BitVector(row>initial_rows for row in 1:M.nnodes(product.mesh))
    cells=D._GmshOptimizeTet[]
    for column in axes(product.mesh.tets,2)
        vertices=ntuple(slot->product.mesh.tets[slot,column],4)
        points=ntuple(slot->ntuple(axis->product.mesh.coords[axis,vertices[slot]],3),4)
        quality=product.radii[column]
        push!(cells,D._GmshOptimizeTet(vertices,fill(Int32(0),4),
            product.element_numbers[column],Int32(0),quality,false))
    end
    state=D._GmshOptimizeState(copy(product.mesh.coords),cells,product.element_high_water,
        Set{NTuple{3,Int32}}(),Set{Tuple{Int32,Int32}}(),movement,
        Dict{NTuple{3,Int32},Int}(),Dict{Tuple{Int32,Int32},Int}(),Set{Int32}())
    D._gopt_connect!(state,eachindex(cells))
    for _ in 1:passes,index in eachindex(cells)
        cells[index].quality<.4 || continue
        for role in 1:4
            D._gopt_smooth_vertex!(state,index,role)
        end
    end
    mesh=M.Mesh(state.coords;tets=copy(product.mesh.tets),tet_tag=copy(product.mesh.tet_tag))
    M.validate(mesh).ok && isequal(mesh.coords[:,1:initial_rows],
        product.mesh.coords[:,1:initial_rows]) || throw(ArgumentError(
        "$caller: native insertion smoothing changes retained constructor rows"))
    return merge(product,(mesh=mesh,))
end

function refine_owned(mesh::M.Mesh,background_size,sizing,caller;
        extend_boundary::Bool=true,minimize_boundary::Bool=false,
        embedded_sizes=Dict{Int32,Float64}(),
        protected_faces=Set{NTuple{3,Int32}}(),
        protected_edges=Set{NTuple{2,Int32}}(),
        element_numbers=nothing,element_high_water=nothing,
        node_numbers=nothing,node_high_water=nothing,
        native_storage_order::Bool=false,return_numbers::Bool=false,region_background::Bool=false,
        metric=nothing,metric_entity=nothing,
        max_nodes::Int=D.PIPELINE_DEFAULT_MAX_NODES,
        max_tets::Int=D.PIPELINE_DEFAULT_MAX_TETS,max_iterations::Int=1_000_000)
    native_storage_order && return_numbers || throw(ArgumentError(
        "$caller: scoped native insertion requires explicit constructor and allocation receipts"))
    max_nodes=D._refine_limit3(max_nodes,"max_nodes";minimum=1)
    max_tets=D._refine_limit3(max_tets,"max_tets";minimum=1)
    max_iterations=D._refine_limit3(max_iterations,"max_iterations";minimum=1)
    nn=size(mesh.coords,2);nt=size(mesh.tets,2)
    (nn<=max_nodes && nt<=max_tets) || throw(ArgumentError(
        "$caller: input exceeds the construction limits"))
    diagnostic=M.validate(mesh)
    diagnostic.ok || throw(ArgumentError("$caller: invalid input — "*join(diagnostic.messages,"; ")))
    T,tags,exterior=D._gmsh3_triangulation(mesh)
    region_background && (length(mesh.tet_tag)!=nt || any(<=(0),mesh.tet_tag)) &&
        throw(ArgumentError("$caller: hybrid Region sizing lacks complete actual cell ownership"))
    numbers=D._gmsh3_birth_numbers(element_numbers,nt,caller,"element")
    number=D._gmsh3_birth_maximum(element_high_water,numbers,caller,"element")
    vertex_numbers=D._gmsh3_birth_numbers(node_numbers,nn,caller,"node")
    vertex_number=D._gmsh3_birth_maximum(node_high_water,vertex_numbers,caller,"node")
    storage=native_storage_order ? NTuple{4,Int32}[
        ntuple(slot->mesh.tets[slot,column],4) for column in axes(mesh.tets,2)] : nothing
    length(sizing)==nn && all(value->isfinite(value) && value>0.0,sizing) ||
        throw(ArgumentError("$caller: owned native sizing does not bind every current constructor row"))
    metric===nothing || metric isa M.AbstractAnisoField || throw(ArgumentError(
        "$caller: native hybrid insertion metric must be an AbstractAnisoField"))
    # The per-cell metric entity mirrors the scalar region background: under
    # `region_background` the owning Region tag selects the field exactly
    # like `background_size(x,y,z,tag)`; otherwise the caller's fixed entity
    # (`(3,tag)` or `nothing`) applies to every evaluation.
    entityof=if metric===nothing
        t->nothing
    elseif region_background
        t->(3,Int(tags[t]))
    else
        t->metric_entity
    end
    sizes=copy(sizing);background=copy(sizing)
    queue=D._GmshTetQueue();skipped=Set{NTuple{4,Int32}}()
    for t in eachindex(T.alive)
        if metric===nothing
            D._gmsh3_queue_tet!(queue,T,t,sizes,background,extend_boundary,skipped,numbers,storage)
        else
            native_queue_tet!(queue,T,t,sizes,background,extend_boundary,
                skipped,numbers,storage,metric,entityof(t))
        end
    end
    live=nt;iteration=0;newtets=Int32[]
    while true
        entry=D._gmsh3_pop!(queue)
        entry===nothing && break
        t=entry.tet
        (T.alive[t] && D._gmsh3_key(T,t)==entry.key && !(entry.key in skipped)) || continue
        iteration<max_iterations || throw(ArgumentError("$caller: iteration limit reached"))
        iteration+=1
        point=storage===nothing ? D._gmsh3_center(T,t) :
            D._gmsh3_native_center(T,storage[t])
        if metric!==nothing
            vertices=storage===nothing ?
                ntuple(slot->D._vert(T,t,slot),4) : storage[t]
            point=native_metric_point(T,vertices,metric,entityof(t),point)
        end
        all(isfinite,point) || throw(ArgumentError("$caller: circumcenter overflows Float64"))
        # The stack search MVertex has a birth even when no containing cell is found.
        vertex_number=Base.checked_add(vertex_number,UInt64(1))
        shell=native_cavity!(T,t,point,tags,protected_faces,storage;
            metric,entityof)
        containing=native_containing(T,point,storage)
        if containing===nothing
            push!(skipped,entry.key);continue
        end
        # The candidate MVertex is allocated before star and distance rejection.
        vertex_number=Base.checked_add(vertex_number,UInt64(1))
        if !native_star!(T,point,tags,protected_faces,storage,shell)
            push!(skipped,entry.key);continue
        end
        owner,weights=containing
        propagated=0.0
        for k in 1:4
            propagated+=weights[k]*sizes[storage[owner][k]]
        end
        target=Float64(region_background ? background_size(point...,tags[owner]) : background_size(point...))
        (isfinite(propagated) && propagated>0 && isfinite(target) && target>0) ||
            throw(ArgumentError("$caller: insertion size must be finite and positive"))
        lc=extend_boundary ? min(propagated,target) : target
        # MTetrahedron objects for every shell facet are constructed before
        # insertVertexB's distance rejection. Their operation identities are
        # therefore consumed even when the candidate cavity is discarded.
        candidates=UInt64(length(T.boundary))
        candidates<=typemax(UInt64)-number || throw(ArgumentError(
            "$caller: candidate element number allocation overflows"))
        first_number=number
        number+=candidates
        # Under the metric arm the `.05` close-point gate reads normalized
        # metric displacements; the topology, shell-conservation and
        # protected-incidence certificates inside `_gmsh3_compatible` are
        # metric-independent and run with the scalar distance gate off.
        compatible=metric===nothing ?
            D._gmsh3_compatible(T,point,lc,protected_faces,protected_edges) :
            (D._gmsh3_compatible(T,point,0.0,protected_faces,protected_edges) &&
             native_metric_admission(T,point,metric,lc,entityof(owner)))
        if !compatible
            push!(skipped,entry.key);continue
        end
        T.nreal<max_nodes || throw(ArgumentError("$caller: node limit reached"))
        remaining=live-length(T.cavity)
        length(T.boundary)<=max_tets-remaining || throw(ArgumentError("$caller: tetrahedron limit reached"))
        region=tags[owner]
        push!(T.x,point[1]);push!(T.y,point[2]);push!(T.z,point[3])
        push!(T.vtet,Int32(0));T.nreal+=1
        push!(sizes,propagated);push!(background,target)
        push!(vertex_numbers,vertex_number)
        for old in T.cavity;D._killtet!(T,old);end
        empty!(newtets)
        D._retriangulate3!(T,T.boundary,Int32(T.nreal),newtets)
        resize!(tags,length(T.alive))
        resize!(numbers,length(T.alive))
        storage===nothing || resize!(storage,length(T.alive))
        for (ordinal,new) in enumerate(newtets)
            tags[new]=region
            numbers[new]=first_number+UInt64(ordinal)
            storage===nothing || (storage[new]=(shell[ordinal].rows...,Int32(T.nreal)))
            if metric===nothing
                D._gmsh3_queue_tet!(queue,T,new,sizes,background,extend_boundary,skipped,numbers,storage)
            else
                native_queue_tet!(queue,T,new,sizes,background,extend_boundary,
                    skipped,numbers,storage,metric,entityof(new))
            end
        end
        live=remaining+length(newtets)
    end
    coordinates=Matrix{Float64}(undef,3,T.nreal)
    for node in 1:T.nreal
        coordinates[1,node]=T.x[node];coordinates[2,node]=T.y[node];coordinates[3,node]=T.z[node]
    end
    live_rows=findall(T.alive)
    live=length(live_rows)
    tets=Matrix{Int32}(undef,4,live)
    output_tags=isempty(mesh.tet_tag) ? Int32[] : Vector{Int32}(undef,live)
    output_numbers=return_numbers ? Vector{UInt64}(undef,live) : UInt64[]
    output_rows=live_rows
    radii=storage===nothing ? Float64[] :
        Float64[D._gmsh3_key(T,t) in skipped ? 0. :
            metric===nothing ?
                D._gmsh3_native_radius(T,storage[t],sizes,background,extend_boundary) :
                native_metric_radius(T,storage[t],sizes,background,
                    extend_boundary,metric,entityof(t))
            for t in output_rows]
    if storage!==nothing
        pool_rows=output_rows
        permutation=sortperm(eachindex(pool_rows);lt=(a,b)->radii[a]>radii[b] ? true :
            radii[a]<radii[b] ? false : numbers[pool_rows[a]]<numbers[pool_rows[b]])
        output_rows=output_rows[permutation]
    end
    column=0
    for t in output_rows
        column+=1
        vertices=storage===nothing ? (D._vert(T,t,1),D._vert(T,t,2),D._vert(T,t,4),D._vert(T,t,3)) : storage[t]
        for slot in 1:4;tets[slot,column]=vertices[slot];end
        isempty(output_tags) || (output_tags[column]=tags[t])
        return_numbers && (output_numbers[column]=numbers[t])
    end
    output=M.Mesh(coordinates;segs=copy(mesh.segs),tris=copy(mesh.tris),tets=tets,
        seg_tag=copy(mesh.seg_tag),tri_tag=copy(mesh.tri_tag),tet_tag=output_tags)
    M.validate(output).ok || throw(ErrorException("$caller: invalid replacement mesh"))
    Set(first(M.boundary_faces(output.tets)))==Set(first(M.boundary_faces(mesh.tets))) ||
        throw(ErrorException("$caller: exterior boundary changed"))
    return return_numbers ? (mesh=output,element_numbers=output_numbers,element_high_water=number,
        node_numbers=vertex_numbers,node_high_water=vertex_number,radii=radii[permutation],
        sizes=sizes,background=background) : output
end

end
