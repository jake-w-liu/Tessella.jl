"""
    Tessella

A Julia-native, robust, memory-efficient mesh generator. Goal: *design → mesh
always works*. See `PLAN.md` for the staged architecture, `DEVELOPMENT.md` for the
mandatory CRC (Correctness–Robustness–Completeness) discipline, and `STATUS.md` for
the live stage board.

Pipeline (0-D → 3-D): exact predicates ([`Predicates`](@ref)) underpin a robust
2-D Delaunay/CDT/refinement core ([`Mesh2D`](@ref)), 1-D edge + surface meshing
([`Mesh1D`](@ref), [`MeshSurface`](@ref)) under size fields ([`SizeField`](@ref)),
and a 3-D Delaunay kernel + volume filling ([`Mesh3D`](@ref)), with quality/
optimization ([`Optimize`](@ref)), surface healing ([`Heal`](@ref)), and gmsh
`.msh` / STL I/O ([`IO`](@ref)). The top-level [`mesh_volume`](@ref) ties them
together with the *validated-or-explicit-blocker* contract.
"""
module Tessella

# ── Version / capability banner ────────────────────────────────────────────────
# Highest development stage whose CRC gate (STATUS.md) is green. Stage 6 includes
# high-order elements and solver-consumable I/O; stages 0–5 are prerequisites.
const TESSELLA_STAGE = 6  # see STATUS.md stage board

# ── Submodules (PLAN.md §3) ────────────────────────────────────────────────────
include("core/Predicates.jl")     # Stage 0: adaptive exact orient/incircle/insphere + SoS
include("core/MeshTypes.jl")      # Stage 0: compact SoA mesh, topology, quality, CRC checksum
include("core/Transform.jl")      # P3: validated affine mesh transformations
include("core/Elements.jl")       # P2: general fixed-node Gmsh element/entity model + mixed MSH I/O
include("meshing/Recombine.jl")   # P4: deterministic triangle-to-quad surface recombination
include("meshing/Refine.jl")      # P4: deterministic one-level uniform simplex refinement
include("structured/Transfinite.jl") # P4: validated four-sided planar transfinite patches
include("structured/TransfiniteCurve.jl") # P4: Gmsh straight-curve transfinite laws
include("structured/TransfiniteTriangle.jl") # P4: three-sided structured triangle patches
include("structured/TransfiniteQuad.jl") # P4: recombined four-sided quadrangle patches
include("structured/TransfiniteVolume.jl") # P4: affine six-face transfinite volumes
include("structured/TransfinitePrism.jl") # P4: affine five-face transfinite prisms
include("structured/TransfiniteHex.jl") # P4: affine six-face recombined hexahedra
include("meshing/ExactMesh3D.jl") # Stage 3: exact-coordinate (Rational{BigInt}) 3-D Delaunay
include("interfaces/IO.jl")       # Stage 0: .msh v2/v4 read/write, STL, .geo scan
include("meshing/Mesh2D.jl")      # Stage 1: 2-D Delaunay + CDT + Ruppert refinement
include("fields/SizeField.jl")    # Stage 2/4: size fields
include("meshing/Mesh1D.jl")      # Stage 2: graded 1-D edge meshing
include("meshing/MeshSurface.jl") # Stage 2: planar / cylinder / parametric surface meshing
include("meshing/Mesh3D.jl")      # Stage 3: 3-D Delaunay + volume filling (+ multi-region)
include("meshing/RecoverCDT.jl")  # Stage 3: general conforming-Delaunay boundary recovery (exact)
include("meshing/Optimize.jl")    # Stage 4: tet quality report + Laplacian smoothing
include("geometry/Heal.jl")       # Stage 5: surface-defect detection ("heal, don't fail")
include("geometry/Geometry.jl")   # Stage 5: native constructive primitive surfaces
include("geometry/CAD.jl")        # Stage 5: native analytical geometry (surfaces + exact imprints), no OCC
include("geometry/NURBS.jl")      # P3: native B-spline/NURBS evaluation
include("geometry/BRep.jl")       # P3: ISO-10303-21 STEP / IGES classified-solid and NURBS import
include("meshing/HighOrder.jl")   # Stage 6: quadratic (P2) tet generation + type-11 I/O
include("geometry/Model.jl")      # P2: tagged geometry/entity kernel
include("geometry/GeoExec.jl")    # P3: bounded .geo execution
include("meshing/BoundaryLayer.jl") # P4: prismatic and 2-D quad/fan boundary-layer extrusion
include("meshing/Periodic.jl")    # P4: periodic identification
include("interfaces/Post.jl")     # P5: views and plugins
include("interfaces/API.jl")      # P5: model/mesh/option façade
include("interfaces/CLI.jl")      # P5: command-line entry
include("interfaces/GUI.jl")      # P5: headless GUI state machine

using .MeshTypes: Mesh, validate, mesh_crc, tet_signed_volume, boundary_faces
using .Transform: affine_transform, translate_mesh, rotate_mesh, dilate_mesh, mirror_mesh
using .Elements: ElementSpec, MSH_CATALOG, msh_spec, msh_num_nodes, msh_dimension,
                 msh_order, msh_family, ElementBlock, ElementRef,
                 SpecialElementBlock, MixedEntity, MixedEntityData, MixedMesh,
                 mixed_crc, simplex_to_mixed, mixed_to_simplex, write_mixed_msh,
                 read_mixed_msh, lagrange_nodes, add_block!
using .Recombine: recombine_triangles
using .Refine: refine_uniform
using .Transfinite: mesh_transfinite_patch
using .TransfiniteCurve: transfinite_curve_parameters, transfinite_curve_hwall
using .TransfiniteTriangle: mesh_transfinite_triangle
using .TransfiniteQuad: mesh_transfinite_quad_patch
using .TransfiniteVolume: mesh_transfinite_volume
using .TransfinitePrism: mesh_transfinite_prism
using .TransfiniteHex: mesh_transfinite_hex
using .Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!,
              add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!,
              embed!, translate_volume!, dilate_volume!, rotate_volume!,
              add_physical_group!, mesh_model_surface, mesh_model_volume
using .NURBS: NURBSCurve, NURBSSurface, nurbs_eval, bspline_basis
using .GeoExec: execute_geo
using .BoundaryLayer: mesh_boundary_layer, mesh_boundary_layer_2d, mesh_boundary_layer_filled
using .Periodic: periodic_identify
using .BRep: import_step, import_iges, import_nurbs_step, import_nurbs_iges,
                  export_iges_nurbs
using .SizeField: AbstractField, AbstractSizeField, AbstractAnisoField, ConstantSize, FunctionSize,
                  DistanceField, ThresholdField, BoxField, BallField, CylinderField,
                  FrustumField, MinSize, MaxSize,
                  BoundedSize, field_value, size_at, metric_at, build_geo_size_field,
                  build_geo_boundary_layer_fields,
                  Metric3, isotropic_metric, metric_size, metric_eigenvalues,
                  directional_size, metric_edge_length,
                  MathEvalField, MathEvalAnisoField, GradientField, LaplacianField,
                  MeanField, CurvatureField, MaxEigenHessianField, LonLatField, ParametricField,
                  StructuredField, RestrictField, ConstantField, ExtendField, OctreeField,
                  PostViewField, MinAnisoField, IntersectAnisoField, AttractorAnisoCurveField,
                  BoundaryLayerField, AutomaticMeshSizeField, ExternalProcessField
using .Mesh1D: mesh_segment
using .Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using .Mesh3D: tetrahedralize, tetrahedralize_multi, tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion, recover_boundary, mesh_boolean, mesh_sized_conforming, mesh_cylinder, refine_to_size
using .RecoverCDT: recover_boundary_cdt, recover_partition_cdt, mesh_sized_cdt
using .Optimize: smooth_laplacian, smooth_odt, smooth_optimize, remove_slivers, mesh_quality
using .Heal: is_meshable

# Install Mesh3D's exact-recovery extension only after both modules and their
# public types are available, avoiding a Mesh3D ↔ RecoverCDT include cycle.
Mesh3D._recover_boundary_exact(surface::Mesh) = recover_boundary_cdt(surface)
Mesh3D._recover_partition_exact(surfaces::AbstractVector{Mesh}) = recover_partition_cdt(surfaces)

export mesh_volume, mesh_planar, mesh_sized_extrude, mesh_sized, refine_to_size, stage
# curated re-exports of the public API
export Mesh, validate, mesh_crc, mesh_quality, is_meshable
export affine_transform, translate_mesh, rotate_mesh, dilate_mesh, mirror_mesh
export ElementSpec, MSH_CATALOG, msh_spec, msh_num_nodes, msh_dimension, msh_order,
       msh_family, ElementBlock, ElementRef, SpecialElementBlock, MixedEntity,
       MixedEntityData, MixedMesh, mixed_crc, simplex_to_mixed, mixed_to_simplex,
       write_mixed_msh, read_mixed_msh, lagrange_nodes, add_block!
export recombine_triangles
export refine_uniform
export mesh_transfinite_patch
export transfinite_curve_parameters, transfinite_curve_hwall
export mesh_transfinite_triangle
export mesh_transfinite_quad_patch
export mesh_transfinite_volume
export mesh_transfinite_prism
export mesh_transfinite_hex
export GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!, add_box!
export add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
export embed!, translate_volume!, dilate_volume!, rotate_volume!
export add_physical_group!, mesh_model_surface, mesh_model_volume
export NURBSCurve, NURBSSurface, nurbs_eval, bspline_basis
export execute_geo, mesh_boundary_layer, mesh_boundary_layer_2d,
       mesh_boundary_layer_filled, periodic_identify
export import_step, import_iges, import_nurbs_step, import_nurbs_iges, export_iges_nurbs
export AbstractField, AbstractSizeField, AbstractAnisoField, ConstantSize, FunctionSize, DistanceField,
       ThresholdField, BoxField, BallField, CylinderField, FrustumField,
       MinSize, MaxSize, BoundedSize, field_value, size_at, metric_at, Metric3,
       isotropic_metric, metric_size, metric_eigenvalues, directional_size,
       metric_edge_length
export MathEvalField, MathEvalAnisoField, GradientField, LaplacianField, MeanField,
       CurvatureField, MaxEigenHessianField, LonLatField, ParametricField, StructuredField,
       RestrictField, ConstantField, ExtendField, OctreeField, PostViewField,
       MinAnisoField, IntersectAnisoField, AttractorAnisoCurveField, BoundaryLayerField,
       AutomaticMeshSizeField, ExternalProcessField
export build_geo_size_field, build_geo_boundary_layer_fields
export tetrahedralize, tetrahedralize_multi, tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion, recover_boundary, recover_boundary_cdt, recover_partition_cdt, mesh_sized_cdt, mesh_boolean, mesh_sized_conforming, mesh_cylinder, smooth_laplacian, smooth_odt, smooth_optimize, remove_slivers

"""
    stage() -> Int

Current implemented development stage (see the `STATUS.md` stage board).
"""
stage() = TESSELLA_STAGE

function _planar_entity_context(value,caller::AbstractString)
    value===nothing && return nothing
    value isa Tuple && length(value)==2 && value[1] isa Integer && value[2] isa Integer ||
        throw(ArgumentError("$caller: entity context must be nothing or a (dimension, tag) integer tuple"))
    (value[1] isa Bool || value[2] isa Bool) && throw(ArgumentError(
        "$caller: entity dimension and tag must not be Bool"))
    dim=try Int(value[1]) catch err
        err isa InterruptException && rethrow()
        err isa InexactError || rethrow()
        throw(ArgumentError("$caller: entity dimension is outside the platform Int range"))
    end
    tag=try Int(value[2]) catch err
        err isa InterruptException && rethrow()
        err isa InexactError || rethrow()
        throw(ArgumentError("$caller: entity tag is outside the platform Int range"))
    end
    dim in 0:3 || throw(ArgumentError("$caller: entity dimension must be in 0:3"))
    tag>0 || throw(ArgumentError("$caller: entity tag must be positive"))
    return (dim,tag)
end

function _planar_vertex_context(spec,npoints::Int,index::Int,p)
    spec===nothing && return nothing
    raw=if spec isa Tuple
        spec
    elseif spec isa AbstractVector
        length(spec)==npoints || throw(ArgumentError(
            "mesh_planar: vertex_entities must have one entry per input point"))
        spec[index]
    elseif spec isa AbstractDict
        haskey(spec,index) || throw(ArgumentError(
            "mesh_planar: vertex_entities has no entry for point $index"))
        spec[index]
    elseif applicable(spec,index,p)
        spec(index,p)
    else
        throw(ArgumentError(
            "mesh_planar: vertex_entities must be a point-entity tuple, per-point vector/dictionary, or callable (index, p)"))
    end
    context=_planar_entity_context(raw,"mesh_planar point $index")
    (context===nothing || context[1]==0) || throw(ArgumentError(
        "mesh_planar: vertex_entities point $index must have dimension 0"))
    return context
end

function _mesh_sized_vertex_context(spec,npoints::Int,index::Int,p)
    spec===nothing && return nothing
    raw=if spec isa Tuple
        spec
    elseif spec isa AbstractVector
        length(spec)==npoints || throw(ArgumentError(
            "mesh_sized: vertex_entities must have one entry per surface point"))
        spec[index]
    elseif spec isa AbstractDict
        get(spec,index,nothing)
    elseif applicable(spec,index,p)
        spec(index,p)
    else
        throw(ArgumentError(
            "mesh_sized: vertex_entities must be a point-entity tuple, per-surface-point vector/dictionary, or callable (index, p)"))
    end
    context=_planar_entity_context(raw,"mesh_sized surface point $index")
    (context===nothing || context[1]==0) || throw(ArgumentError(
        "mesh_sized: vertex_entities point $index must have dimension 0"))
    return context
end

function _planar_boundary_context(spec,nsegments::Int,index::Int,p,q,fallback)
    spec===nothing && return fallback
    raw = if spec isa Tuple
        spec
    elseif spec isa AbstractVector
        length(spec)==nsegments || throw(ArgumentError(
            "mesh_planar: boundary_entities must have one entry per segment"))
        spec[index]
    elseif spec isa AbstractDict
        haskey(spec,index) || throw(ArgumentError(
            "mesh_planar: boundary_entities has no entry for segment $index"))
        spec[index]
    elseif applicable(spec,index,p,q)
        spec(index,p,q)
    else
        throw(ArgumentError(
            "mesh_planar: boundary_entities must be an entity tuple, per-segment vector/dictionary, or callable (index, p, q)"))
    end
    return _planar_entity_context(raw,"mesh_planar boundary segment $index")
end

function _grade_planar_constraints(xs::Vector{Float64},ys::Vector{Float64},segments,
                                   field::AbstractSizeField,face_entity,boundary_entities,
                                   vertex_entities)
    outx=copy(xs);outy=copy(ys);outsegments=Tuple{Int,Int}[]
    npoints=length(xs);nsegments=length(segments)
    for (index,segment) in pairs(segments)
        a=Int(segment[1]);b=Int(segment[2])
        (1<=a<=npoints && 1<=b<=npoints && a!=b) || throw(ArgumentError(
            "mesh_planar: segment $index has invalid endpoints ($a,$b) for $npoints points"))
        p=(xs[a],ys[a],0.0);q=(xs[b],ys[b],0.0)
        edge_entity=_planar_boundary_context(boundary_entities,nsegments,Int(index),
                                             p,q,face_entity)
        begin_entity=_planar_vertex_context(vertex_entities,npoints,a,p)
        end_entity=_planar_vertex_context(vertex_entities,npoints,b,q)
        points,_=mesh_segment(p,q,field;entity=edge_entity,
                              endpoint_entities=(begin_entity,end_entity))
        previous=a
        @inbounds for k in 2:length(points)-1
            length(outx)<typemax(Int32) || throw(ArgumentError(
                "mesh_planar: graded boundary node count exceeds Int32 indexing"))
            push!(outx,points[k][1]);push!(outy,points[k][2])
            current=length(outx)
            push!(outsegments,(previous,current));previous=current
        end
        push!(outsegments,(previous,b))
    end
    return outx,outy,outsegments
end

"""
    mesh_planar(xs, ys, segments; min_angle_deg=25.0, max_area=Inf, rng_seed=1,
                field=nothing, entity=nothing, boundary_entities=nothing,
                vertex_entities=nothing) -> Mesh

Quality 2-D triangle mesh of the planar straight-line graph (points `(xs,ys)` +
constraint `segments`, `(i,j)` index pairs): constrained-Delaunay triangulate,
Ruppert-refine to the angle/area bound and optional [`AbstractSizeField`](@ref),
keep the interior of the constrained domain, and return a validated 2-D
[`Mesh`](@ref) (nodes carry `z = 0`). The domain boundary must be closed
constrained loops. `entity` classifies interior field queries. `boundary_entities`
can separately classify each input constraint as one entity tuple, a vector or
dictionary keyed by segment index, or a callable `(index, p, q) -> entity`; this
is required for curve-sensitive `Restrict` and `Constant` fields.
`vertex_entities` similarly accepts a point entity, per-point vector/dictionary,
or `(index, p) -> entity` callable and applies Gmsh's endpoint size rule. This is
the 2-D counterpart of [`mesh_volume`](@ref).
"""
function mesh_planar(xs::Vector{Float64}, ys::Vector{Float64},
                     segments::AbstractVector{<:Tuple{Integer,Integer}};
                     min_angle_deg::Real=25.0, max_area::Real=Inf, rng_seed::Integer=1,
                     field::Union{Nothing,AbstractSizeField}=nothing, entity=nothing,
                     boundary_entities=nothing,vertex_entities=nothing)
    face_entity=_planar_entity_context(entity,"mesh_planar")
    if field===nothing
        (boundary_entities===nothing && vertex_entities===nothing) ||
            throw(ArgumentError(
                "mesh_planar: boundary_entities and vertex_entities require a size field"))
        mesh_xs=xs;mesh_ys=ys;mesh_segments=segments
    else
        mesh_xs,mesh_ys,mesh_segments=_grade_planar_constraints(
            xs,ys,segments,field,face_entity,boundary_entities,vertex_entities)
    end
    T = constrained_delaunay(mesh_xs, mesh_ys, mesh_segments; rng_seed=rng_seed)
    sizefn = (field===nothing || field isa AbstractAnisoField) ? nothing :
             (x,y)->size_at(field,x,y,0.0,face_entity)
    edgefn = field isa AbstractAnisoField ?
             (ax,ay,bx,by)->metric_edge_length(field,(ax,ay,0.0),(bx,by,0.0);
                                                    entity=face_entity) : nothing
    interior = refine!(T; min_angle_deg=min_angle_deg, max_area=max_area,
                       size=sizefn,edge_metric=edgefn)
    m = to_mesh(T; interior=interior)
    size(m.tris, 2) > 0 ||
        throw(ArgumentError("mesh_planar: constrained boundary encloses no triangles"))
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_planar: produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

"""
    mesh_sized_extrude(xs, ys, segments, z0, z1; hmax, min_angle_deg=25.0) -> Mesh

**Uniform size-controlled** tet mesh of the prismatic (extruded) domain whose
cross-section is the polygon `(xs, ys)` with boundary/hole `segments`, extruded from
`z=z0` to `z=z1`, with a **guaranteed maximum edge length `≤ hmax`**. The cross-section
is meshed by the size-controlled 2-D Ruppert mesher (edges `≤ hmax/√2`), then extruded
into `⌈(z1−z0)/(hmax/√2)⌉` layers; each triangular prism is split into three tets by a
**column-global-index diagonal rule**, so every shared quad face picks the same diagonal
in both incident prisms — the result is **conforming** (no face shared by >2 tets),
watertight, provably valid (each tet oriented positive), and of **exact volume**
(area(polygon)·(z1−z0), boundary preserved — no jitter). The `√2` factor bounds the
prism face-diagonal `√(edge²+layer²) ≤ hmax`.

Handles **non-convex** cross-sections and polygons with holes (via the `segments`
constrained-Delaunay + interior classification). This is the extruded/prismatic case of
uniform sizing on an arbitrary domain — complementing [`mesh_box`](@ref) (box) and
[`mesh_cylinder`](@ref) (cylinder). For a general closed faceted surface, use
[`mesh_sized`](@ref).
"""
function mesh_sized_extrude(xs::AbstractVector, ys::AbstractVector,
                            segments::AbstractVector{<:Tuple{Integer,Integer}},
                            z0::Real, z1::Real; hmax::Real, min_angle_deg::Real=25.0)
    hm=try Float64(hmax) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: hmax must be Float64-representable: $(sprint(showerror,err))"))
    end
    za=try Float64(z0) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: z0 must be Float64-representable: $(sprint(showerror,err))"))
    end
    zb=try Float64(z1) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: z1 must be Float64-representable: $(sprint(showerror,err))"))
    end
    (isfinite(hm)&&hm>0) || throw(ArgumentError("mesh_sized_extrude: hmax must be finite and positive (got $hmax)"))
    (isfinite(za)&&isfinite(zb)&&zb>za) ||
        throw(ArgumentError("mesh_sized_extrude: need finite z1 > z0 (got z0=$z0, z1=$z1)"))
    h2 = hm/sqrt(2.0)                                      # prism face-diagonal ≤ hmax
    (isfinite(h2)&&h2>0) || throw(ArgumentError("mesh_sized_extrude: hmax is below Float64 spacing resolution"))
    T2 = constrained_delaunay(Float64.(xs), Float64.(ys), segments)
    interior = refine!(T2; min_angle_deg=min_angle_deg, max_area=0.5*h2*h2, size=(cx,cy)->h2)
    m2 = to_mesh(T2; interior=interior)                    # 2-D mesh: coords (z=0), tris
    nv2 = size(m2.coords, 2)
    nv2 >= 3 && size(m2.tris, 2) >= 1 ||
        throw(ErrorException("mesh_sized_extrude: cross-section produced no interior triangles (check segments orientation/closure)"))
    ratio=(zb-za)/h2
    (isfinite(ratio)&&ratio<=typemax(Int)) ||
        throw(ArgumentError("mesh_sized_extrude: axial layer count exceeds the platform Int limit"))
    nz=max(1,ceil(Int,ratio));dz=(zb-za)/nz
    levels=try Base.checked_add(nz,1) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: axial layer count overflows Int"))
    end
    nout=try Base.checked_mul(levels,nv2) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: node count overflows Int"))
    end
    nout<=typemax(Int32) || throw(ArgumentError("mesh_sized_extrude: $nout nodes exceed Int32"))
    coords = Matrix{Float64}(undef, 3, nout)
    @inbounds for k in 0:nz, i in 1:nv2
        id = k*nv2 + i
        coords[1,id] = m2.coords[1,i]; coords[2,id] = m2.coords[2,i]; coords[3,id] = za + k*dz
    end
    _n(id) = @inbounds (coords[1,id], coords[2,id], coords[3,id])
    ntout=try Base.checked_mul(3,Base.checked_mul(size(m2.tris,2),nz)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: tetrahedron count overflows Int"))
    end
    tets = Matrix{Int32}(undef, 4, ntout); col = 0
    @inbounds for ti in 1:size(m2.tris,2)
        v = sort!(Int[m2.tris[1,ti], m2.tris[2,ti], m2.tris[3,ti]])   # v[1]<v[2]<v[3] (column-index rule)
        for k in 0:nz-1
            b1=k*nv2+v[1]; b2=k*nv2+v[2]; b3=k*nv2+v[3]
            t1=(k+1)*nv2+v[1]; t2=(k+1)*nv2+v[2]; t3=(k+1)*nv2+v[3]
            for te in ((b1,b2,b3,t3), (b1,b2,t3,t2), (b1,t2,t3,t1))
                a,b,c,d = te
                tet_signed_volume(_n(a),_n(b),_n(c),_n(d)) < 0 && ((c,d) = (d,c))   # orient positive
                col += 1; tets[1,col]=a; tets[2,col]=b; tets[3,col]=c; tets[4,col]=d
            end
        end
    end
    m = Mesh(coords; tets=tets)
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_sized_extrude: produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

"""
    mesh_sized(surface::Mesh; hmax) -> Mesh
    mesh_sized(surface::Mesh; field) -> Mesh
    mesh_sized(surface::Mesh, field::AbstractSizeField) -> Mesh

Size-controlled tet mesh of the domain enclosed by an arbitrary closed triangulated
`surface` (curved or polyhedral, convex or non-convex).  Pass exactly one of a uniform
`hmax` or a spatial [`AbstractSizeField`](@ref).  Uniform sizing guarantees maximum edge
length `≤ hmax`; field sizing enforces each edge against the minimum field value sampled
at its endpoints and midpoint. Two stages:

1. **Fill** — [`tetrahedralize`](@ref) (Delaunay + ray-cast interior classification)
   gives a valid conforming tet mesh of the domain; if that is not watertight (a
   non-convex boundary the Delaunay restriction misses), [`recover_boundary`](@ref) is
   used instead.
2. **Size** — [`refine_to_size`](@ref) longest-edge-bisects the fill until every edge is
   `≤ hmax`. Bisection is conforming, validity/volume-preserving, and terminating, so the
   result is a valid, watertight, conforming mesh with `maxedge ≤ hmax` at the domain's
   own (faceted-surface) volume.

Verified on convex boxes, spheres, conical frusta, and a non-convex L-prism (guaranteed
`maxedge ≤ hmax`, valid, watertight, conforming). Returns the mesh, or throws an explicit
blocker if the surface cannot be filled watertight (never a silently bad mesh). For
axis-aligned boxes prefer [`mesh_box`](@ref) and for extrusions [`mesh_sized_extrude`](@ref)
(exact boundary, no faceting); `mesh_sized` is the general fallback for arbitrary surfaces.
`entity` and `entity_resolver` have the classification semantics documented by
[`refine_to_size`](@ref); input surface triangles/segments and their tags are
reattached to the conforming fill before field refinement. `vertex_entities`
classifies the input surface vertices as point entities. It accepts one `(0,tag)`
tuple, a per-surface-point vector, a sparse dictionary keyed by point index, or a
callable `(index, point) -> entity`; point classifications are remapped through
the conforming fill and used to seed the original classified boundary edges;
the point value is not recursively imposed on the resulting child edges.
"""
function mesh_sized(surface::Mesh; hmax::Union{Nothing,Real}=nothing,
                    field::Union{Nothing,AbstractSizeField}=nothing, entity=nothing,
                    entity_resolver=nothing,vertex_entities=nothing)
    (hmax===nothing) == (field===nothing) &&
        throw(ArgumentError("mesh_sized: pass exactly one of hmax or field"))
    target = if hmax===nothing
        field::AbstractSizeField
    else
        hm=try Float64(hmax) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("mesh_sized: hmax must be Float64-representable: $(sprint(showerror,err))"))
        end
        (isfinite(hm)&&hm>0) ||
            throw(ArgumentError("mesh_sized: hmax must be finite and positive (got $hmax)"))
        ConstantSize(hm)
    end
    m0,mapping = _attach_surface_classification(
        _conforming_fill(surface;caller="mesh_sized"),surface,"mesh_sized")
    point_contexts=if vertex_entities===nothing
        nothing
    else
        contexts=Vector{Union{Nothing,Tuple{Int,Int}}}(nothing,size(m0.coords,2))
        @inbounds for index in axes(surface.coords,2)
            point=(surface.coords[1,index],surface.coords[2,index],surface.coords[3,index])
            contexts[mapping[index]]=_mesh_sized_vertex_context(
                vertex_entities,size(surface.coords,2),index,point)
        end
        contexts
    end
    m = refine_to_size(m0,target;entity=entity,entity_resolver=entity_resolver,
                       vertex_entities=point_contexts)
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_sized: refinement produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

mesh_sized(surface::Mesh,field::AbstractSizeField;entity=nothing,
           entity_resolver=nothing,vertex_entities=nothing)=
    mesh_sized(surface;field=field,entity=entity,entity_resolver=entity_resolver,
               vertex_entities=vertex_entities)

@inline _classification_coord_key(x::Float64,y::Float64,z::Float64)=
    (x==0 ? 0.0 : x,y==0 ? 0.0 : y,z==0 ? 0.0 : z)

function _attach_surface_classification(volume::Mesh,surface::Mesh,
                                        caller::AbstractString)
    ids=Dict{NTuple{3,Float64},Int32}()
    @inbounds for i in axes(volume.coords,2)
        key=_classification_coord_key(volume.coords[1,i],volume.coords[2,i],
                                      volume.coords[3,i])
        haskey(ids,key) && throw(ErrorException(
            "$caller: conforming fill contains duplicate represented coordinates"))
        ids[key]=Int32(i)
    end
    mapping=Vector{Int32}(undef,size(surface.coords,2))
    @inbounds for i in eachindex(mapping)
        key=_classification_coord_key(surface.coords[1,i],surface.coords[2,i],
                                      surface.coords[3,i])
        mapped=get(ids,key,Int32(0))
        mapped!=0 || throw(ErrorException(
            "$caller: conforming fill did not preserve surface vertex $i at $key"))
        mapping[i]=mapped
    end
    S=Matrix{Int32}(undef,2,size(surface.segs,2))
    @inbounds for s in axes(surface.segs,2),i in 1:2
        S[i,s]=mapping[surface.segs[i,s]]
    end
    F=Matrix{Int32}(undef,3,size(surface.tris,2))
    @inbounds for f in axes(surface.tris,2),i in 1:3
        F[i,f]=mapping[surface.tris[i,f]]
    end
    out=Mesh(copy(volume.coords);segs=S,tris=F,tets=copy(volume.tets),
             seg_tag=copy(surface.seg_tag),tri_tag=copy(surface.tri_tag),
             tet_tag=copy(volume.tet_tag))
    diagnostic=validate(out)
    diagnostic.ok || throw(ErrorException(
        "$caller: classified conforming fill is invalid — "*join(diagnostic.messages,"; ")))
    return out,mapping
end

# Fill a closed PLC without ever accepting a convex-hull cap as its boundary.
# `tetrahedralize` owns the restriction, star-shaped fan, exact-CDT recovery, and
# bounded Float64 fallback, and returns only after Mesh3D's exact geometric gate.
function _conforming_fill(surface::Mesh; caller::AbstractString,
                          rng_seed::Integer=1, optimize::Bool=false)
    try
        return tetrahedralize(surface; rng_seed=rng_seed, optimize=optimize)
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        throw(ErrorException("$caller: could not construct a conforming volume mesh — " *
                             sprint(showerror, err)))
    end
end

"""
    mesh_volume(surface; smooth=true, smooth_iters=5, rng_seed=1, check=true) -> Mesh

Robust volume mesh of the closed triangulated `surface`, with the "always-valid or
explicit blocker" contract: the surface is first screened by [`is_meshable`](@ref)
(unless `check=false`); a defect raises an `ArgumentError` carrying the precise
[`Heal.SurfaceReport`](@ref) — never a silent bad mesh. Otherwise an initial
[`tetrahedralize`](@ref) fill must pass the exact PLC-boundary conformity gate; a
non-conforming convex-hull cap falls through [`recover_boundary`](@ref) and then
[`recover_boundary_cdt`](@ref). If `smooth`, the conforming fill is Laplacian-smoothed
([`smooth_laplacian`](@ref)). The result is `validate`-checked before return.
"""
function mesh_volume(surface::Mesh; smooth::Bool=true, smooth_iters::Integer=5,
                     optimize::Bool=false, rng_seed::Integer=1, check::Bool=true)
    if check
        ok, report = is_meshable(surface)
        ok || throw(ArgumentError("mesh_volume: input surface is not meshable — " *
                                  join(report.messages, "; ")))
    end
    m = check ? _conforming_fill(surface; caller="mesh_volume", rng_seed=rng_seed,
                                 optimize=optimize) :
                tetrahedralize(surface; rng_seed=rng_seed, optimize=optimize, check=false)
    smooth && (m = smooth_laplacian(m; iters=smooth_iters))
    # optimization-based (min-dihedral) smoothing targets slivers the mean-smoother and
    # topological flips leave behind; only in the optimize path (it is a local search).
    optimize && (m = smooth_optimize(m; iters=smooth_iters))
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_volume: produced an invalid mesh — " *
                                    join(diag.messages, "; ")))
    return m
end

end # module Tessella
