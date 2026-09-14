# General affine entity transforms, `Duplicata` copies, and post-transform
# coherence for the native `GeoModel`.
#
# Semantics mirror the Gmsh built-in kernel's `ApplicationOnShapes` +
# `ReplaceAllDuplicates` pipeline (`Geo.cpp`): a transform mutates the shared
# vertices of each listed entity in place, recursively through boundary
# topology, with every vertex moved at most once per statement; a coincident
# entity merge then restores coherence. `Duplicata` deep-copies entities at
# their current positions (fresh tags allocated through Gmsh's shared
# `Geometry.OldNewReg` counter) and is merge-free on its own — the copies are
# normally collapsed or spread by the enclosing transform statement.

const _COHERENCE_RTOL = 1e-8

# A fused affine map plus the exact sequence of 4×4 Gmsh `ApplicationOnShapes`
# matrices that realize it. `linear`/`offset` describe the map semantically
# (encoding representability checks); `steps` are the per-pass matrices Gmsh
# applies to each owned vertex, in order — `Rotate` is three passes
# (translate by -origin, rotate, translate by +origin), so coordinates match
# Gmsh's rounding bit-for-bit.
struct _AffineTransform
    linear::NTuple{9,Float64}   # row-major 3×3
    offset::NTuple{3,Float64}
    steps::Vector{NTuple{16,Float64}}  # row-major 4×4, homogeneous last column
end

# `vecmat4x4` from Geo.cpp: res[i] = Σ_j mat[i][j]·vec[j] with vec=(x,y,z,1),
# accumulated in order.
@inline function _gmsh_matvec4x4(mat::NTuple{16,Float64}, p::NTuple{3,Float64})
    v=(p[1],p[2],p[3],1.0)
    return ntuple(3) do i
        acc=0.0
        for j in 1:4
            acc+=mat[4*(i-1)+j]*v[j]
        end
        acc
    end
end

@inline function _affine_apply_steps(t::_AffineTransform, p::NTuple{3,Float64})
    for step in t.steps
        p=_gmsh_matvec4x4(step,p)
    end
    return p
end

@inline function _affine_apply(t::_AffineTransform, p::NTuple{3,Float64})
    (a,b,c,d,e,f,g,h,i)=t.linear
    return (a*p[1]+b*p[2]+c*p[3]+t.offset[1],
            d*p[1]+e*p[2]+f*p[3]+t.offset[2],
            g*p[1]+h*p[2]+i*p[3]+t.offset[3])
end

@inline function _linear_apply(t::_AffineTransform, v::NTuple{3,Float64})
    (a,b,c,d,e,f,g,h,i)=t.linear
    return (a*v[1]+b*v[2]+c*v[3], d*v[1]+e*v[2]+f*v[3], g*v[1]+h*v[2]+i*v[3])
end

@inline _dot(a::NTuple{3,Float64},b::NTuple{3,Float64}) =
    a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

@inline function _entity_label(dim::Int)
    dim==0 && return "Point"
    return _model_periodic_entity_label(dim)
end

@inline _gmsh_translation_step(d::NTuple{3,Float64}) =
    (1.0,0.0,0.0,d[1], 0.0,1.0,0.0,d[2], 0.0,0.0,1.0,d[3], 0.0,0.0,0.0,1.0)

function _affine_translation(delta, caller)
    d=_finite_vector3(delta,caller,"translation delta")
    return _AffineTransform((1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),d,
                            [_gmsh_translation_step(d)])
end

function _affine_dilation(center, scale, caller)
    c=_finite_vector3(center,caller,"dilation center")
    scales=scale isa Real ?
        let s=_finite_scalar(scale,caller,"dilation scale"); (s,s,s) end :
        _finite_vector3(scale,caller,"dilation scales")
    (sx,sy,sz)=scales
    # SetDilatationMatrix: mat[i][3] = T[i]*(1-scale_i).
    step=(sx,0.0,0.0,c[1]*(1.0-sx),
          0.0,sy,0.0,c[2]*(1.0-sy),
          0.0,0.0,sz,c[3]*(1.0-sz),
          0.0,0.0,0.0,1.0)
    return _AffineTransform((sx,0.0,0.0,0.0,sy,0.0,0.0,0.0,sz),
                            (c[1]*(1-sx),c[2]*(1-sy),c[3]*(1-sz)),[step])
end

# `norme`/`prodve` from Gmsh's Numeric.h: normalization multiplies by the
# reciprocal; the cross product uses Gmsh's component order.
@inline function _gmsh_norme!(v::Vector{Float64})
    mod=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])
    if mod!=0.0
        inv=1.0/mod
        v[1]*=inv;v[2]*=inv;v[3]*=inv
    end
    return v
end

@inline _gmsh_prodve(a::Vector{Float64},b::Vector{Float64}) =
    [a[2]*b[3]-a[3]*b[2], -a[1]*b[3]+a[3]*b[1], a[1]*b[2]-a[2]*b[1]]

# SetRotationMatrix from Geo.cpp: orthonormal basis (axis, t1, t2) built by
# Gmsh's GramSchmidt, an in-basis rotation about the axis, then
# plan'·rot·plan with Gmsh's triple-loop order.
function _gmsh_rotation_matrix(axis::NTuple{3,Float64}, angle::Float64)
    axe=[axis[1],axis[2],axis[3]]
    if axe[1]!=0.0
        t1=[0.0,1.0,0.0];t2=[0.0,0.0,1.0]
    elseif axe[2]!=0.0
        t1=[1.0,0.0,0.0];t2=[0.0,0.0,1.0]
    else
        t1=[1.0,0.0,0.0];t2=[0.0,1.0,0.0]
    end
    # GramSchmidt(axe, t1, t2): v1=axis, v2=t1, v3=t2.
    _gmsh_norme!(axe)
    t1=_gmsh_norme!(_gmsh_prodve(t2,axe))
    t2=_gmsh_norme!(_gmsh_prodve(axe,t1))
    plan=(axe,t1,t2)   # rows
    c,s=cos(angle),sin(angle)
    rot=((1.0,0.0,0.0),(0.0,c,-s),(0.0,s,c))
    interm=ntuple(3) do i
        ntuple(3) do j
            acc=0.0
            for k in 1:3
                acc+=plan[k][i]*rot[k][j]   # invplan[i][k] = plan[k][i]
            end
            acc
        end
    end
    mat=ntuple(3) do i
        ntuple(3) do j
            acc=0.0
            for k in 1:3
                acc+=interm[i][k]*plan[k][j]
            end
            acc
        end
    end
    return (mat[1][1],mat[1][2],mat[1][3],0.0,
            mat[2][1],mat[2][2],mat[2][3],0.0,
            mat[3][1],mat[3][2],mat[3][3],0.0,
            0.0,0.0,0.0,1.0)
end

function _affine_rotation(axis, origin, angle, caller)
    a=_finite_vector3(axis,caller,"rotation axis")
    o=_finite_vector3(origin,caller,"rotation origin")
    θ=_finite_scalar(angle,caller,"rotation angle")
    n=sqrt(a[1]^2+a[2]^2+a[3]^2)
    n>0 || throw(ArgumentError("$caller: rotation axis must be nonzero"))
    rstep=_gmsh_rotation_matrix(a,θ)
    L=(rstep[1],rstep[2],rstep[3],rstep[5],rstep[6],rstep[7],
       rstep[9],rstep[10],rstep[11])
    Ro=_gmsh_matvec4x4(rstep,o)
    return _AffineTransform(L,(o[1]-Ro[1],o[2]-Ro[2],o[3]-Ro[3]),
                            [_gmsh_translation_step((-o[1],-o[2],-o[3])),
                             rstep,
                             _gmsh_translation_step(o)])
end

# Reflection across the plane `A*x + B*y + C*z + D = 0` (Gmsh `Symmetry`). A
# degenerate zero normal takes Gmsh's `p -> 1e-12` floor and acts as identity.
function _affine_symmetry(a, b, c, d, caller)
    A=_finite_scalar(a,caller,"symmetry plane coefficient A")
    B=_finite_scalar(b,caller,"symmetry plane coefficient B")
    C=_finite_scalar(c,caller,"symmetry plane coefficient C")
    D=_finite_scalar(d,caller,"symmetry plane coefficient D")
    p=A*A+B*B+C*C
    p==0.0 && (p=1e-12)
    F=-2.0/p
    step=(1.0+A*A*F, A*B*F, A*C*F, A*D*F,
          A*B*F, 1.0+B*B*F, B*C*F, B*D*F,
          A*C*F, B*C*F, 1.0+C*C*F, C*D*F,
          0.0,0.0,0.0,1.0)
    L=(step[1],step[2],step[3],step[5],step[6],step[7],step[9],step[10],step[11])
    return _AffineTransform(L,(step[4],step[8],step[12]),[step])
end

"""
    transform_entities!(model, transform::_AffineTransform, entities; caller) -> entities

Apply an affine transform in place to the shared vertices of each listed
`(dimension, tag)` entity. Boundary topology is traversed recursively —
surfaces contribute their loops' curves' endpoints and volumes their shells' —
and each point is moved at most once per call, matching Gmsh's per-statement
transformed-point set. Compact volume encodings are updated when the transform
is representable (isotropic for spheres; axis-perpendicular isotropic for
cylinders and cones; recomputed corners for materialized `add_box!` volumes)
and are discarded with the shell topology preserved otherwise. Boolean operand
snapshots are transformed with the result; snapshots recorded for other
Booleans stay untouched, matching Gmsh's operation-time snapshot semantics. A
coincident-entity coherence merge runs afterwards. The call is atomic:
non-representable transforms and non-finite results leave the model unchanged.
"""
function transform_entities!(m::GeoModel, t::_AffineTransform,
                             entities::AbstractVector{<:Tuple{Integer,Integer}};
                             caller::AbstractString="transform_entities!")
    normalized=NTuple{2,Int}[]
    for (dim,tag) in entities
        d=_dimension(dim,caller); tg=_tag(tag,caller,d)
        haskey(m.discrete,(d,tg)) && throw(ArgumentError(
            "$caller: discrete $(_entity_label(d))[$tg] is not transformable"))
        push!(normalized,(d,tg))
    end
    move=Int[]; seen=Set{Int}(); volumes=Int[]
    for (d,tg) in normalized
        for p in _entity_point_tags(m,d,tg,caller)
            p in seen || (push!(seen,p); push!(move,p))
        end
        d==3 && push!(volumes,tg)
    end
    plans=[_plan_volume_transform(m,tg,t,caller) for tg in volumes]
    coords=[(tag,_finite_result(_affine_apply_steps(t,m.points[tag]),caller))
            for tag in move]
    for (tag,p) in coords
        m.points[tag]=p
    end
    for plan in plans
        _apply_volume_plan!(m,plan)
    end
    _reconcile_box_encodings!(m,seen)
    coherence!(m)
    return normalized
end

# Point tags owned by one entity — the transform's atomic unit. Parent entities
# recurse through boundary topology; duplicate endpoints are deduplicated by the
# caller's per-statement set.
function _entity_point_tags(m::GeoModel, d::Int, tag::Int, caller)
    if d==0
        haskey(m.points,tag) || throw(ArgumentError("$caller: unknown Point[$tag]"))
        return (tag,)
    elseif d==1
        haskey(m.curves,tag) || throw(ArgumentError("$caller: unknown Curve[$tag]"))
        a,b=m.curves[tag]
        cps=get(m.curve_control_points,tag,Int[])
        return (a,b,cps...)
    elseif d==2
        haskey(m.surfaces,tag) || throw(ArgumentError(
            "$caller: unknown Surface[$tag]"))
        pts=Int[]
        for l in m.surfaces[tag]
            for c in m.loops[l]
                a,b=m.curves[abs(c)]
                push!(pts,a,b)
                append!(pts,get(m.curve_control_points,abs(c),Int[]))
            end
        end
        return pts
    else
        haskey(m.volumes,tag) || throw(ArgumentError("$caller: unknown Volume[$tag]"))
        pts=collect(_model_volume_owned_points(m,tag))
        for shell in m.volumes[tag], signed_surface in m.surface_loops[shell],
            loop in m.surfaces[abs(signed_surface)], signed_curve in m.loops[loop]
            append!(pts,get(m.curve_control_points,abs(signed_curve),Int[]))
        end
        return pts
    end
end

# A staged volume-encoding update, validated against pre-transform state and
# applied only after every point move commits.
function _plan_volume_transform(m::GeoModel, tag::Int, t::_AffineTransform,
                                caller)
    if haskey(m.cylinders,tag)
        rec=m.cylinders[tag]
        tr=_transform_axis_encoding(t,rec.center,rec.axis,caller,"Cylinder[$tag]")
        return (dict=:cylinders,tag=tag,
                record=(center=tr.center,axis=tr.axis,radius=rec.radius*tr.perp,
                        height=tr.height))
    elseif haskey(m.spheres,tag)
        rec=m.spheres[tag]
        s2=_linear_isotropy(t.linear)
        s2===nothing && throw(ArgumentError(
            "$caller: transform is not representable on implicit Sphere[$tag] — " *
            "it requires an isotropic linear part"))
        return (dict=:spheres,tag=tag,
                record=(center=_finite_result(
                            _affine_apply_steps(t,rec.center),caller),
                        radius=rec.radius*sqrt(s2)))
    elseif haskey(m.cones,tag)
        rec=m.cones[tag]
        tr=_transform_axis_encoding(t,rec.center,rec.axis,caller,"Cone[$tag]")
        return (dict=:cones,tag=tag,
                record=(center=tr.center,axis=tr.axis,r1=rec.r1*tr.perp,
                        r2=rec.r2*tr.perp,height=tr.height))
    elseif haskey(m.booleans,tag)
        A,B=m.boolean_operands[tag]
        return (dict=:boolean_operands,tag=tag,
                record=(_transform_mesh_snapshot(A,t,caller),
                        _transform_mesh_snapshot(B,t,caller)))
    else
        return (dict=nothing,tag=tag,record=nothing)
    end
end

function _apply_volume_plan!(m::GeoModel, plan)
    plan.dict===nothing && return nothing
    getfield(m,plan.dict)[plan.tag]=plan.record
    return nothing
end

function _transform_mesh_snapshot(mesh::Mesh, t::_AffineTransform, caller)
    coords=Matrix{Float64}(undef,3,size(mesh.coords,2))
    for i in axes(mesh.coords,2)
        coords[:,i].=_finite_result(
            _affine_apply_steps(
                t,(mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])),
            caller)
    end
    return Mesh(coords; segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
end

# Squared isotropic scale of the linear part, or `nothing` when its Gram matrix
# is not a scaled identity.
function _linear_isotropy(L::NTuple{9,Float64})
    cols=((L[1],L[4],L[7]),(L[2],L[5],L[8]),(L[3],L[6],L[9]))
    d=(_dot(cols[1],cols[1]),_dot(cols[2],cols[2]),_dot(cols[3],cols[3]))
    s2=sum(d)/3
    s2>0 || return nothing
    tol=1e-12*max(1.0,s2)
    (abs(d[1]-s2)<=tol && abs(d[2]-s2)<=tol && abs(d[3]-s2)<=tol &&
     abs(_dot(cols[1],cols[2]))<=tol && abs(_dot(cols[1],cols[3]))<=tol &&
     abs(_dot(cols[2],cols[3]))<=tol) || return nothing
    return s2
end

# Transform an axis-encoded primitive: the axis is a free vector (linear part
# only, length = height) and the cross-section radii scale by the transform's
# stretch perpendicular to it — representable only when that stretch is
# isotropic.
function _transform_axis_encoding(t::_AffineTransform, center, axis, caller, what)
    a=_linear_apply(t,axis)
    n=sqrt(_dot(a,a))
    n>0 || throw(ArgumentError("$caller: transform collapses the axis of $what"))
    (u,v)=_perp_basis(a ./ n)
    uL=_linear_apply(t,u); vL=_linear_apply(t,v)
    g=(_dot(uL,uL),_dot(uL,vL),_dot(vL,vL))
    s2=(g[1]+g[3])/2
    s2>0 || throw(ArgumentError("$caller: transform collapses $what"))
    tol=1e-12*max(1.0,s2)
    (abs(g[1]-s2)<=tol && abs(g[3]-s2)<=tol && abs(g[2])<=tol) || throw(ArgumentError(
        "$caller: transform is not representable on $what — it would warp the cross-section"))
    return (center=_finite_result(_affine_apply_steps(t,center),caller),
            axis=_finite_result(a,caller),perp=sqrt(s2),height=n)
end

function _perp_basis(n::NTuple{3,Float64})
    ax,ay,az=abs.(n)
    seed=ax<=ay ? (ax<=az ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
                  (ay<=az ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    u=(n[2]*seed[3]-n[3]*seed[2],n[3]*seed[1]-n[1]*seed[3],
       n[1]*seed[2]-n[2]*seed[1])
    u=u ./ sqrt(_dot(u,u))
    v=(n[2]*u[3]-n[3]*u[2],n[3]*u[1]-n[1]*u[3],n[1]*u[2]-n[2]*u[1])
    return (u,v)
end

# Materialized `add_box!` volumes keep their `box_extents` encoding only while
# the shell still owns exactly the encoded corners (`_model_volume_bounds`
# validates this). After a point move, re-derive extents when the moved corners
# still form an axis-aligned box; otherwise drop the encoding so queries and
# meshing fall back to the explicit shell rather than reading stale bounds.
function _reconcile_box_encodings!(m::GeoModel, moved::Set{Int})
    for (tag,ext) in collect(m.box_extents)
        owned=sort!(collect(_model_volume_owned_points(m,tag)))
        any(p->p in moved,owned) || continue
        x0,y0,z0,dx,dy,dz=ext
        scale=max(1.0,abs(x0),abs(y0),abs(z0),abs(dx),abs(dy),abs(dz))
        tol=1e-9*scale
        keep=false
        if length(owned)==8
            mapped=[m.points[p] for p in owned]
            lo=ntuple(k->minimum(getindex.(mapped,k)),3)
            hi=ntuple(k->maximum(getindex.(mapped,k)),3)
            ext2=hi .- lo
            corners=[(lo[1]+ix*ext2[1],lo[2]+iy*ext2[2],lo[3]+iz*ext2[3])
                     for ix in (0,1),iy in (0,1),iz in (0,1)]
            keep=all(e->e>0,ext2) && length(unique(mapped))==8 && all(
                c->any(p->_points_close(p,c,tol),mapped),corners)
            keep && (m.box_extents[tag]=(lo[1],lo[2],lo[3],
                                         ext2[1],ext2[2],ext2[3]))
        end
        keep || delete!(m.box_extents,tag)
    end
    return nothing
end

"""
    coherence!(model; tol) -> Bool

Merge coincident entities the way Gmsh's `ReplaceAllDuplicates` does after a
transform: points within `tol` (relative to the model's coordinate scale) fuse
into the lowest tag, then curves sharing an endpoint pair, then surfaces whose
flattened boundary-curve sequences match. References are rewired through
physical groups, embeddings, periodic relations, discrete boundaries, compound
members, and transfinite corner lists. Automatic tag counters reset to the
largest surviving tag per merged dimension, matching `Geometry.OldNewReg`.
"""
function coherence!(m::GeoModel; tol::Real=_COHERENCE_RTOL)
    scale=isempty(m.points) ? 1.0 :
          max(1.0,maximum(p->maximum(abs.(p)),values(m.points)))
    eps=Float64(tol)*scale
    merged=_merge_points!(m,eps)
    merged |= _merge_curves!(m)
    merged |= _merge_surfaces!(m)
    return merged
end

"""
    merge_vertices!(model, tags; caller)

Gmsh's `Coherence Point{...}` / `mergeVertices`: snap every listed point to
the first tag's coordinates, then run the global coherence merge. Unknown
tags raise an error before any position is changed.
"""
function merge_vertices!(m::GeoModel, tags;
                         caller::AbstractString="merge_vertices!")
    length(tags)<2 && return nothing
    haskey(m.points,tags[1]) || throw(ArgumentError(
        "$caller: unknown Point[$(tags[1])]"))
    for t in tags[2:end]
        haskey(m.points,t) || throw(ArgumentError("$caller: unknown Point[$t]"))
    end
    target=m.points[tags[1]]
    for t in tags[2:end]
        m.points[t]=target
    end
    coherence!(m)
    return nothing
end

function _merge_points!(m::GeoModel,eps)
    isempty(m.points) && return false
    grid=Dict{NTuple{3,Int},Vector{Int}}()
    mapping=Dict{Int,Int}()
    for tag in sort!(collect(keys(m.points)))
        p=m.points[tag]
        cell=ntuple(i->floor(Int,p[i]/eps),3)
        keep=0
        for dx in -1:1, dy in -1:1, dz in -1:1
            for other in get(grid,(cell[1]+dx,cell[2]+dy,cell[3]+dz),Int[])
                if _points_close(p,m.points[other],eps)
                    keep=other; break
                end
            end
            keep!=0 && break
        end
        if keep==0
            push!(get!(grid,cell,Int[]),tag)
        else
            mapping[tag]=keep
        end
    end
    isempty(mapping) && return false
    for (drop,keep) in sort!(collect(mapping))
        _rewire_entity_refs!(m,0,drop,keep)
        for (c,(a,b)) in collect(m.curves)
            (a==drop || b==drop) &&
                (m.curves[c]=(a==drop ? keep : a, b==drop ? keep : b))
        end
        for (c,cps) in collect(m.curve_control_points)
            rewired=[p==drop ? keep : p for p in cps]
            rewired==cps || (m.curve_control_points[c]=rewired)
        end
        delete!(m.points,drop); delete!(m.point_size,drop)
        _drop_entity_state!(m,0,drop)
    end
    m.next_tag[1]=isempty(m.points) ? 0 : maximum(keys(m.points))
    return true
end

function _merge_curves!(m::GeoModel)
    seen=Dict{NTuple{2,Int},Int}()
    mapping=Dict{Int,Int}()
    for tag in sort!(collect(keys(m.curves)))
        key=m.curves[tag]
        haskey(seen,key) ? (mapping[tag]=seen[key]) : (seen[key]=tag)
    end
    isempty(mapping) && return false
    for (drop,keep) in sort!(collect(mapping))
        _rewire_entity_refs!(m,1,drop,keep)
        for l in keys(m.loops)
            m.loops[l]=[v==drop ? keep : (v==-drop ? -keep : v)
                        for v in m.loops[l]]
        end
        delete!(m.curves,drop)
        # The dropped curve's control points stay behind as unreferenced
        # vertices, matching Gmsh's destroyed-copy behavior.
        delete!(m.curve_control_points,drop)
        _drop_entity_state!(m,1,drop)
    end
    m.next_tag[2]=isempty(m.curves) ? 0 : maximum(keys(m.curves))
    return true
end

function _merge_surfaces!(m::GeoModel)
    seen=Dict{Vector{Int},Int}()
    mapping=Dict{Int,Int}()
    for tag in sort!(collect(keys(m.surfaces)))
        flat=Int[]
        for l in m.surfaces[tag]; append!(flat,abs.(m.loops[l])); end
        haskey(seen,flat) ? (mapping[tag]=seen[flat]) : (seen[flat]=tag)
    end
    isempty(mapping) && return false
    for (drop,keep) in sort!(collect(mapping))
        _rewire_entity_refs!(m,2,drop,keep)
        for sl in keys(m.surface_loops)
            m.surface_loops[sl]=[v==drop ? keep : (v==-drop ? -keep : v)
                                 for v in m.surface_loops[sl]]
        end
        orphan_loops=m.surfaces[drop]
        delete!(m.surfaces,drop)
        _drop_entity_state!(m,2,drop)
        referenced=Set{Int}()
        for loops in values(m.surfaces); union!(referenced,loops); end
        for l in orphan_loops
            l in referenced || delete!(m.loops,l)
        end
    end
    m.next_tag[3]=isempty(m.surfaces) ? 0 : maximum(keys(m.surfaces))
    return true
end

# Replace references to the dropped entity in every `(dim, tag)`-indexed or
# member-list store. Signed topology memberships (curve loops, surface loops)
# are rewired by the dedicated merge passes above.
function _rewire_entity_refs!(m::GeoModel, dim::Int, drop::Int, keep::Int)
    for (key,members) in m.physical
        key[1]==dim && drop in members &&
            (m.physical[key]=unique!([v==drop ? keep : v for v in members]))
    end
    embeds=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (key,entities) in m.embeds
        new_key=key==(dim,drop) ? (dim,keep) : key
        merged=unique!([e==(dim,drop) ? (dim,keep) : e for e in entities])
        if haskey(embeds,new_key)
            append!(embeds[new_key],merged); unique!(embeds[new_key])
        else
            embeds[new_key]=merged
        end
    end
    m.embeds=embeds
    periodic=Dict{Tuple{Int,Int},ModelPeriodicConstraint}()
    for (key,c) in m.periodic
        c.dim==dim || (periodic[key]=c; continue)
        Int(c.slave_entity)==drop && continue   # slave entity is gone
        master=Int(c.master_entity)==drop ? keep : Int(c.master_entity)
        periodic[key]=ModelPeriodicConstraint(c.dim,c.slave_entity,
                                              Int32(master),c.affine,
                                              c.reversed,c.atol)
    end
    m.periodic=periodic
    discrete=Dict{Tuple{Int,Int},DiscreteEntity}()
    for (key,ent) in m.discrete
        new_key=key==(dim,drop) ? (dim,keep) : key
        ent.boundary=unique!([b==(dim,drop) ? (dim,keep) : b
                              for b in ent.boundary])
        discrete[new_key]=ent
    end
    m.discrete=discrete
    attached=Dict{Tuple{Int,Int},DiscreteEntity}()
    for (key,ent) in m.meshing.attached
        new_key=key==(dim,drop) ? (dim,keep) : key
        ent.boundary=unique!([b==(dim,drop) ? (dim,keep) : b
                              for b in ent.boundary])
        attached[new_key]=ent
    end
    m.meshing.attached=attached
    m.meshing.compounds=Pair{Int,Vector{Int}}[
        cdim==dim ? cdim=>unique!([t==drop ? keep : t for t in ctags]) :
                    cdim=>copy(ctags)
        for (cdim,ctags) in m.meshing.compounds]
    if dim==0
        for (s,a) in m.meshing.transfinite_surfaces
            drop in a.corners && (m.meshing.transfinite_surfaces[s]=
                (arrangement=a.arrangement,
                 corners=[c==drop ? keep : c for c in a.corners]))
        end
        for v in keys(m.meshing.transfinite_volumes)
            corners=m.meshing.transfinite_volumes[v]
            drop in corners && (m.meshing.transfinite_volumes[v]=
                [c==drop ? keep : c for c in corners])
        end
    end
    for i in eachindex(m.meshing.homology_requests)
        r=m.meshing.homology_requests[i]
        (drop in r.domain || drop in r.subdomain || drop in r.dims) &&
            (m.meshing.homology_requests[i]=
                (kind=r.kind,
                 domain=[t==drop ? keep : t for t in r.domain],
                 subdomain=[t==drop ? keep : t for t in r.subdomain],
                 dims=[t==drop ? keep : t for t in r.dims]))
    end
    return nothing
end

function _drop_entity_state!(m::GeoModel, dim::Int, tag::Int)
    key=(dim,tag)
    delete!(m.entity_names,key)
    delete!(m.entity_visibility,key)
    delete!(m.entity_colors,key)
    delete!(m.meshing.recombine,key)
    delete!(m.meshing.smoothing,key)
    delete!(m.meshing.reverse,key)
    delete!(m.meshing.algorithm,key)
    delete!(m.meshing.size_at_params,key)
    delete!(m.meshing.size_from_boundary,key)
    delete!(m.meshing.attached,key)
    if dim==1
        delete!(m.meshing.transfinite_curves,tag)
        delete!(m.curve_control_points,tag)
    elseif dim==2
        delete!(m.meshing.transfinite_surfaces,tag)
    elseif dim==3
        delete!(m.meshing.transfinite_volumes,tag)
        delete!(m.meshing.outward_orientation,tag)
    end
    return nothing
end

# ── Duplicata ────────────────────────────────────────────────────────────────

# Gmsh's shared `Geometry.OldNewReg` tag source: the maximum of the
# curve/surface/volume counters plus every live curve loop, surface loop,
# physical, and non-point discrete tag.
function _geo_newreg_tag(m::GeoModel)
    mx=max(m.next_tag[2],m.next_tag[3],m.next_tag[4],m.physical_tag_max)
    for t in keys(m.loops); mx=max(mx,t); end
    for (d,t) in keys(m.discrete); d>0 && (mx=max(mx,t)); end
    for t in keys(m.surface_loops); mx=max(mx,t); end
    return mx+1
end

_geo_newreg_alloc!(m::GeoModel,dim::Int,caller) =
    _alloc_tag!(m,dim,_geo_newreg_tag(m),caller)

"""
    duplicate_entities!(model, entities; caller) -> copies

Deep-copy each listed native entity at its current position, mirroring Gmsh's
`Duplicata`: copies carry no physical memberships or names, curve copies
allocate their two coincident endpoint vertices plus Gmsh's transient beg/end
vertex pair, and surface/volume copies recurse through boundary topology with
tags drawn from the `OldNewReg` counter. Coincident results stay unmerged —
an enclosing transform's coherence pass collapses them exactly like Gmsh.
"""
function duplicate_entities!(m::GeoModel,
                             entities::AbstractVector{<:Tuple{Integer,Integer}};
                             caller::AbstractString="duplicate_entities!")
    out=NTuple{2,Int}[]
    for (dim,tag) in entities
        d=_dimension(dim,caller); tg=_tag(tag,caller,d)
        haskey(m.discrete,(d,tg)) && throw(ArgumentError(
            "$caller: discrete $(_entity_label(d))[$tg] cannot be duplicated"))
        push!(out,(d,_duplicate_entity!(m,d,tg,caller)))
    end
    return out
end

function _duplicate_entity!(m::GeoModel, d::Int, tag::Int, caller)
    if d==0
        haskey(m.points,tag) || throw(ArgumentError("$caller: unknown Point[$tag]"))
        return _fresh_point_copy!(m,tag,caller)
    elseif d==1
        return _duplicate_curve!(m,tag,caller)
    elseif d==2
        return _duplicate_surface!(m,tag,caller)
    else
        return _duplicate_volume!(m,tag,caller)
    end
end

function _fresh_point_copy!(m::GeoModel, src::Int, caller)
    t=_alloc_tag!(m,0,0,caller)
    m.points[t]=m.points[src]
    haskey(m.point_size,src) && (m.point_size[t]=m.point_size[src])
    return t
end

# Gmsh's `DuplicateCurve` allocates the curve tag first (shared `NEWREG`
# counter), then one coincident copy per control point — stored on the copy as
# its `curve_control_points` — and finally the wired beg/end vertex copies.
# For a plain line the source control list is its endpoint pair, so each
# curve copy burns four point tags: two orphans, then the wired pair.
function _duplicate_curve!(m::GeoModel, src::Int, caller)
    haskey(m.curves,src) || throw(ArgumentError("$caller: unknown Curve[$src]"))
    a,b=m.curves[src]
    t=_geo_newreg_alloc!(m,1,caller)
    source_cps=get(m.curve_control_points,src,Int[a,b])
    m.curve_control_points[t]=
        [_fresh_point_copy!(m,c,caller) for c in source_cps]
    pa=_fresh_point_copy!(m,a,caller)
    pb=_fresh_point_copy!(m,b,caller)
    m.curves[t]=(pa,pb)
    return t
end

function _duplicate_surface!(m::GeoModel, src::Int, caller)
    haskey(m.surfaces,src) || throw(ArgumentError("$caller: unknown Surface[$src]"))
    t=_geo_newreg_alloc!(m,2,caller)
    loops=Int[]
    for (i,l) in enumerate(m.surfaces[src])
        curves=[sign(c)*_duplicate_curve!(m,abs(c),caller) for c in m.loops[l]]
        lt=(i==1 && !haskey(m.loops,t)) ? t :
           ((isempty(m.loops) ? 0 : maximum(keys(m.loops)))+1)
        m.loops[lt]=curves
        push!(loops,lt)
    end
    m.surfaces[t]=loops
    return t
end

function _duplicate_volume!(m::GeoModel, src::Int, caller)
    haskey(m.volumes,src) || throw(ArgumentError("$caller: unknown Volume[$src]"))
    t=_geo_newreg_alloc!(m,3,caller)
    if !isempty(m.volumes[src])
        shells=Int[]
        for (i,sl) in enumerate(m.volumes[src])
            surfs=[sign(s)*_duplicate_surface!(m,abs(s),caller)
                   for s in m.surface_loops[sl]]
            slt=(i==1 && !haskey(m.surface_loops,t)) ? t :
                ((isempty(m.surface_loops) ? 0 : maximum(keys(m.surface_loops)))+1)
            m.surface_loops[slt]=surfs
            push!(shells,slt)
        end
        m.volumes[t]=shells
    else
        m.volumes[t]=Int[]
    end
    haskey(m.box_extents,src) && (m.box_extents[t]=m.box_extents[src])
    haskey(m.cylinders,src) && (m.cylinders[t]=m.cylinders[src])
    haskey(m.spheres,src) && (m.spheres[t]=m.spheres[src])
    haskey(m.cones,src) && (m.cones[t]=m.cones[src])
    if haskey(m.booleans,src)
        m.booleans[t]=m.booleans[src]
        m.boolean_operands[t]=m.boolean_operands[src]
    end
    return t
end
