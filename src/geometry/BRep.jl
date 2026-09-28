"""
    BRep

Native ISO-10303-21 STEP and IGES CAD import. Solids that classify as an
axis-aligned block, sphere, right circular cylinder, right circular cone, or
torus are converted to Tessella surfaces and filled. Closed-shell
`MANIFOLD_SOLID_BREP`/`FACETED_BREP`/`BREP_WITH_VOIDS` topology bounded by
planar faces over straight edges imports as the exact polygonal boundary; a
single unbounded `TOROIDAL_SURFACE` face imports as a parametric torus
boundary. IGES solids bounded by trimmed plane surfaces over line edges
(144→142/141→102/110 or the 186/514/510/508/504/502 MSBO chain) and IGES
160/198 tori take the same paths. NURBS curves and surfaces
(STEP `B_SPLINE_*` / IGES 126/128) import as [`NURBSCurve`](@ref) /
[`NURBSSurface`](@ref). Unrecognized topology is an explicit blocker, not a
silent empty mesh.
"""
module BRep

using ..Geometry: box_surface, cylinder_surface, sphere_surface, cone_surface
using ..Mesh3D: tetrahedralize
using ..MeshTypes: Mesh, validate, ntets, tet_volume, node
using ..NURBS: NURBSCurve, NURBSSurface

export import_step, import_iges, parse_step_entities
export import_nurbs_step, import_nurbs_iges, export_iges_nurbs

# A compact CAD record cannot legitimately expand one numeric multiplicity into
# an unbounded allocation. One million knot values already permits models far
# beyond the degree/control counts exercised by the native meshing path.
const _MAX_BREP_KNOT_VALUES=1_000_000

# ── STEP (ISO-10303-21) ───────────────────────────────────────────────────────

struct StepEntity
    id::Int
    kind::String
    args::Vector{Any}
end

@inline function _step_skipws(s, i, last)
    while i<=last && isspace(s[i]); i=nextind(s,i); end
    return i
end

"""
    parse_step_entities(source) -> Dict{Int,StepEntity}

Parse the `DATA` section of an ISO-10303-21 STEP string into entities keyed by
their positive STEP identifier. The bounded native parser accepts the scalar,
reference, enum, string, nested-list, and complex-entity syntax needed by
Tessella's classified primitive and NURBS import paths. Duplicate identifiers,
nonfinite numbers, and malformed syntax raise `ArgumentError`.
"""
function parse_step_entities(source::AbstractString)
    matched=match(r"(?s)\bDATA\s*;(.*)ENDSEC\s*;"i, source)
    matched===nothing && throw(ArgumentError("import_step: missing DATA section"))
    data=matched.captures[1]
    entities=Dict{Int,StepEntity}()
    i=firstindex(data); last=lastindex(data)
    while i<=last
        i=_step_skipws(data,i,last)
        i>last && break
        data[i]=='#' || throw(ArgumentError("import_step: expected entity at $(repr(data[i]))"))
        j=nextind(data,i); idstart=j
        while j<=last && isdigit(data[j]); j=nextind(data,j); end
        j>idstart || throw(ArgumentError("import_step: entity is missing a numeric identifier"))
        id=tryparse(Int,data[idstart:prevind(data,j)])
        (id!==nothing && id>0) || throw(ArgumentError(
            "import_step: entity identifier is outside the positive platform Int range"))
        haskey(entities,id) && throw(ArgumentError("import_step: duplicate entity #$id"))
        j=_step_skipws(data,j,last)
        j<=last && data[j]=='=' || throw(ArgumentError("import_step: entity #$id missing '='"))
        j=_step_skipws(data,nextind(data,j),last)
        if j<=last && data[j]=='('
            args, k=_step_parse_complex(data,j,id)
            kind="COMPLEX"
        else
            k=j
            while k<=last && (isletter(data[k]) || isdigit(data[k]) || data[k]=='_')
                k=nextind(data,k)
            end
            k>j || throw(ArgumentError("import_step: entity #$id is missing a type"))
            kind=uppercase(data[j:prevind(data,k)])
            k=_step_skipws(data,k,last)
            k<=last && data[k]=='(' || throw(ArgumentError("import_step: entity #$id missing '('"))
            args, k=_step_parse_list(data,k)
        end
        k=_step_skipws(data,k,last)
        k<=last && data[k]==';' || throw(ArgumentError("import_step: entity #$id missing ';'"))
        entities[id]=StepEntity(id,kind,args)
        i=nextind(data,k)
    end
    isempty(entities) && throw(ArgumentError("import_step: DATA section is empty"))
    return entities
end

function _step_parse_complex(s, start, id)
    s[start]=='(' || throw(ArgumentError("import_step: expected complex entity"))
    i=nextind(s,start); parts=Any[]; last=lastindex(s)
    while i<=last
        i=_step_skipws(s,i,last)
        i<=last || throw(ArgumentError("import_step: unterminated complex entity #$id"))
        s[i]==')' && return parts, nextind(s,i)
        k=i
        while k<=last && (isletter(s[k]) || isdigit(s[k]) || s[k]=='_')
            k=nextind(s,k)
        end
        k>i || throw(ArgumentError("import_step: complex entity #$id missing type"))
        kind=uppercase(s[i:prevind(s,k)])
        k=_step_skipws(s,k,last)
        args=Any[]
        if k<=last && s[k]=='('
            args, k=_step_parse_list(s,k)
        end
        push!(parts, StepEntity(id,kind,args))
        i=k
    end
    throw(ArgumentError("import_step: unterminated complex entity #$id"))
end

function _step_parse_list(s, start)
    s[start]=='(' || throw(ArgumentError("import_step: expected list"))
    i=nextind(s,start); items=Any[]
    last=lastindex(s)
    while i<=last
        i=_step_skipws(s,i,last)
        i<=last || throw(ArgumentError("import_step: unterminated list"))
        if s[i]==')'
            isempty(items) || throw(ArgumentError("import_step: missing list parameter"))
            return items, nextind(s,i)
        end
        if s[i]=='('
            sub,i=_step_parse_list(s,i); push!(items,sub)
        elseif s[i]=='\''
            j=nextind(s,i)
            while j<=last
                if s[j]=='\''
                    following=nextind(s,j)
                    # STEP represents an embedded apostrophe with two
                    # consecutive apostrophes inside the same string.
                    following<=last && s[following]=='\'' || break
                    j=nextind(s,following)
                else
                    j=nextind(s,j)
                end
            end
            j<=last || throw(ArgumentError("import_step: unterminated string"))
            push!(items,replace(s[nextind(s,i):prevind(s,j)],"''"=>"'"))
            i=nextind(s,j)
        elseif s[i]=='$' || s[i]=='*'
            push!(items,nothing); i=nextind(s,i)
        elseif s[i]=='.'
            j=nextind(s,i)
            while j<=last && s[j]!='.'; j=nextind(s,j); end
            j<=last || throw(ArgumentError("import_step: unterminated enum"))
            push!(items,s[i:j]); i=nextind(s,j)
        elseif s[i]=='#'
            j=nextind(s,i); k=j
            while k<=last && isdigit(s[k]); k=nextind(s,k); end
            k>j || throw(ArgumentError("import_step: reference is missing an identifier"))
            ref=tryparse(Int,s[j:prevind(s,k)])
            (ref!==nothing && ref>0) || throw(ArgumentError(
                "import_step: reference identifier is outside the positive platform Int range"))
            push!(items,(:ref,ref)); i=k
        elseif isletter(s[i]) || s[i]=='_'
            # Typed parameter such as LENGTH_MEASURE(1.E-07): parse it as a
            # nested id-less entity so callers can pattern-match its kind.
            k=i
            while k<=last && (isletter(s[k]) || isdigit(s[k]) || s[k]=='_')
                k=nextind(s,k)
            end
            kind=uppercase(s[i:prevind(s,k)])
            k=_step_skipws(s,k,last)
            k<=last && s[k]=='(' || throw(ArgumentError(
                "import_step: expected '(' after typed parameter $kind"))
            subargs,k=_step_parse_list(s,k)
            push!(items,StepEntity(0,kind,subargs)); i=k
        else
            j=i
            while j<=last && (isdigit(s[j]) || s[j] in ('.','e','E','+','-')); j=nextind(s,j); end
            raw=s[i:prevind(s,j)]
            v=tryparse(Float64,raw)
            v===nothing && throw(ArgumentError("import_step: bad number $raw"))
            isfinite(v) || throw(ArgumentError("import_step: non-finite number $raw"))
            push!(items,v); i=j
        end
        i=_step_skipws(s,i,last)
        i<=last || throw(ArgumentError("import_step: unterminated list"))
        s[i]==')' && return items, nextind(s,i)
        s[i]==',' || throw(ArgumentError("import_step: expected ',' between list parameters"))
        i=nextind(s,i)
    end
    throw(ArgumentError("import_step: unterminated list"))
end

function _step_points(entities)
    points=NTuple{3,Float64}[]
    for ent in values(entities)
        ent.kind=="CARTESIAN_POINT" || continue
        length(ent.args)>=2 || continue
        coords=ent.args[end]
        coords isa Vector && length(coords)>=3 || continue
        x,y,z=Float64(coords[1]),Float64(coords[2]),Float64(coords[3])
        all(isfinite,(x,y,z)) || throw(ArgumentError("import_step: non-finite CARTESIAN_POINT"))
        push!(points,(x,y,z))
    end
    return points
end

function _unique_sorted(values)
    u=sort!(unique(values))
    return u
end

function _as_box(points)
    isempty(points) && return nothing
    uniq=unique(points)
    xs=_unique_sorted([p[1] for p in uniq])
    ys=_unique_sorted([p[2] for p in uniq])
    zs=_unique_sorted([p[3] for p in uniq])
    (length(xs)>=2 && length(ys)>=2 && length(zs)>=2) || return nothing
    x0,x1=xs[1],xs[end]; y0,y1=ys[1],ys[end]; z0,z1=zs[1],zs[end]
    corners=Set((x,y,z) for x in (x0,x1), y in (y0,y1), z in (z0,z1))
    Set(uniq) >= corners || return nothing
    for p in uniq
        (x0<=p[1]<=x1 && y0<=p[2]<=y1 && z0<=p[3]<=z1) || return nothing
        (p[1]==x0 || p[1]==x1 || p[2]==y0 || p[2]==y1 || p[3]==z0 || p[3]==z1) ||
            return nothing
    end
    return (x0,x1,y0,y1,z0,z1)
end

function _first_ref(args)
    for arg in args
        arg isa Tuple && arg[1]===:ref && return arg
    end
    return nothing
end

function _deref_point(entities, arg)
    arg isa Tuple && arg[1]===:ref || return nothing
    ent=get(entities,arg[2],nothing)
    ent===nothing && return nothing
    if ent.kind=="CARTESIAN_POINT"
        c=ent.args[end]
        c isa Vector && length(c)>=3 || return nothing
        return (Float64(c[1]),Float64(c[2]),Float64(c[3]))
    elseif ent.kind=="AXIS2_PLACEMENT_3D"
        for a in ent.args
            p=_deref_point(entities,a)
            p!==nothing && return p
        end
    end
    return nothing
end

function _deref_direction(entities, arg)
    arg isa Tuple && arg[1]===:ref || return nothing
    ent=get(entities,arg[2],nothing)
    ent===nothing && return nothing
    if ent.kind=="DIRECTION"
        d=ent.args[end]
        d isa Vector && length(d)>=3 || return nothing
        return (Float64(d[1]),Float64(d[2]),Float64(d[3]))
    elseif ent.kind=="AXIS2_PLACEMENT_3D" && length(ent.args)>=3
        return _deref_direction(entities, ent.args[3])
    end
    return nothing
end

function _as_sphere(entities)
    for ent in values(entities)
        ent.kind=="SPHERE" || continue
        length(ent.args)>=2 || continue
        r=Float64(ent.args[end])
        r>0 || throw(ArgumentError("import_step: SPHERE radius must be positive"))
        center=_deref_point(entities, _first_ref(ent.args))
        center===nothing && (center=(0.0,0.0,0.0))
        return (center,r)
    end
    return nothing
end

function _as_cylinder(entities)
    for ent in values(entities)
        ent.kind=="RIGHT_CIRCULAR_CYLINDER" || continue
        length(ent.args)>=3 || continue
        height=Float64(ent.args[end-1]); radius=Float64(ent.args[end])
        (height>0 && radius>0) || throw(ArgumentError(
            "import_step: cylinder height and radius must be positive"))
        place=_first_ref(ent.args)
        center=_deref_point(entities, place)
        center===nothing && (center=(0.0,0.0,0.0))
        axis=_deref_direction(entities, place)
        axis===nothing && (axis=(0.0,0.0,1.0))
        return (center,axis,radius,height)
    end
    return nothing
end

function _as_cone(entities)
    for ent in values(entities)
        ent.kind=="RIGHT_CIRCULAR_CONE" || continue
        length(ent.args)>=4 || continue
        height=Float64(ent.args[end-2]); r1=Float64(ent.args[end-1]); r2=Float64(ent.args[end])
        (height>0 && (r1>0 || r2>0) && r1>=0 && r2>=0) || throw(ArgumentError(
            "import_step: cone height must be positive and at least one radius must be positive"))
        place=_first_ref(ent.args)
        center=_deref_point(entities, place)
        center===nothing && (center=(0.0,0.0,0.0))
        axis=_deref_direction(entities, place)
        axis===nothing && (axis=(0.0,0.0,1.0))
        return (center,axis,r1,r2,height)
    end
    return nothing
end

"""
    import_step(path; fill=true) -> Mesh

Import one classified STEP solid: a synthetic axis-aligned eight-corner point
block, `SPHERE`, `RIGHT_CIRCULAR_CYLINDER`, `RIGHT_CIRCULAR_CONE`,
`TOROIDAL_SURFACE`, or a `MANIFOLD_SOLID_BREP`/`FACETED_BREP`/`BREP_WITH_VOIDS`
closed shell bounded by planar faces over straight edges (a single unbounded
`TOROIDAL_SURFACE` face also imports as a full torus). By default the
classified boundary is tetrahedralized; `fill=false` returns its validated
closed triangle surface. Multiple recognized solids, nonplanar or nonlinear
face geometry, open shells, mixed point-cloud topology, and files without a
recognized solid are explicit blockers. Use [`import_nurbs_step`](@ref) for
STEP B-splines.
"""
function import_step(path::AbstractString; fill::Bool=true)
    isfile(path) || throw(ArgumentError("import_step: missing file $path"))
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "import_step: $path is not an ISO-10303-21 STEP file"))
    entities=parse_step_entities(source)
    breps=[ent for ent in sort!(collect(values(entities)); by=e->e.id)
           if ent.kind in _STEP_BREP_SOLID_KINDS]
    if !isempty(breps)
        length(breps)>1 && throw(ArgumentError(
            "import_step: multiple recognized solids are not supported in one import"))
        return _step_brep_solid(only(breps),entities,fill)
    end
    points=_step_points(entities)
    primitive_count=count(ent->ent.kind in ("SPHERE","RIGHT_CIRCULAR_CYLINDER",
                                             "RIGHT_CIRCULAR_CONE","TOROIDAL_SURFACE"),
                          values(entities))
    primitive_count<=1 || throw(ArgumentError(
        "import_step: multiple recognized solids are not supported in one import"))
    sph=_as_sphere(entities)
    if sph!==nothing
        surface=sphere_surface(sph[1],sph[2])
        fill || return surface
        return _filled(surface,"import_step")
    end
    cyl=_as_cylinder(entities)
    if cyl!==nothing
        surface=cylinder_surface(cyl[1],cyl[2],cyl[3],cyl[4])
        fill || return surface
        return _filled(surface,"import_step")
    end
    cone=_as_cone(entities)
    if cone!==nothing
        surface=cone_surface(cone[1],cone[2],cone[3],cone[4],cone[5])
        fill || return surface
        return _filled(surface,"import_step")
    end
    tor=_as_torus(entities)
    if tor!==nothing
        surface=_torus_surface(tor[1],tor[2],tor[3],tor[4],tor[5],"import_step")
        fill || return surface
        return _filled(surface,"import_step")
    end
    point_cloud_only=all(ent->ent.kind=="CARTESIAN_POINT",values(entities))
    box=point_cloud_only ? _as_box(points) : nothing
    if box!==nothing
        surface=box_surface(box...)
        fill || return surface
        return _filled(surface,"import_step")
    end
    kinds=sort!(unique(ent.kind for ent in values(entities)))
    throw(ArgumentError(
        "import_step: no supported solid (axis-aligned 8-corner block, SPHERE, " *
        "RIGHT_CIRCULAR_CYLINDER, RIGHT_CIRCULAR_CONE, TOROIDAL_SURFACE, or " *
        "closed-shell planar polyhedron); saw $(join(kinds, ", ")). " *
        "NURBS curves/surfaces use import_nurbs_step"))
end

function _filled(surface::Mesh, caller)
    mesh=tetrahedralize(surface)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: fill produced an invalid mesh — "*
                                    join(diag.messages,"; ")))
    ntets(mesh)>0 || throw(ErrorException("$caller: fill produced no tetrahedra"))
    return mesh
end

# ── closed-shell BRep solids (planar polyhedra and full tori) ─────────────────

const _STEP_BREP_SOLID_KINDS=("MANIFOLD_SOLID_BREP","FACETED_BREP","BREP_WITH_VOIDS")
const _STEP_BREP_FACE_KINDS=("FACE","FACE_SURFACE","ADVANCED_FACE")
const _STEP_BREP_BOUND_KINDS=("FACE_BOUND","FACE_OUTER_BOUND")
const _STEP_BREP_WRAP_KINDS=("SURFACE_CURVE","SEAM_CURVE","TRIMMED_CURVE")

# A full torus boundary is sampled pole-free on a regular (u,v) grid; the two
# periodic directions are merged by construction, not by tolerance. The grid is
# sized for ~1% fill-volume error while keeping the boundary modest.
const _BREP_TORUS_NU=48
const _BREP_TORUS_NV=24
const _BREP_MAX_SURFACE_NODES=2_000_000
const _BREP_MAX_SURFACE_TRIS=8_000_000

function _deref_entity(entities,arg)
    arg isa Tuple && arg[1]===:ref || return nothing
    return get(entities,arg[2],nothing)
end

@inline _brep_enum(arg)=arg==".T." ? true : arg==".F." ? false : nothing
@inline _brep_vsub(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _brep_vdot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _brep_vcross(a,b)=(a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])

function _brep_vunit(a,caller)
    n=_brep_vdot(a,a)
    (isfinite(n)&&n>0) || throw(ArgumentError(
        "$caller: direction has no positive length"))
    s=sqrt(n)
    return (a[1]/s,a[2]/s,a[3]/s)
end

function _brep_ortho(d,z,caller)
    t=_brep_vdot(d,z)
    w=(d[1]-t*z[1],d[2]-t*z[2],d[3]-t*z[3])
    return _brep_vunit(w,caller)
end

function _brep_default_x(z,caller)
    ax=abs(z[1]); ay=abs(z[2]); az=abs(z[3])
    ref=ax<=ay ? (ax<=az ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
                 (ay<=az ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    return _brep_ortho(ref,z,caller)
end

function _step_vertex_point(entities,arg)
    ent=_deref_entity(entities,arg)
    ent===nothing && return nothing
    if ent.kind=="VERTEX_POINT"
        return _deref_point(entities,_first_ref(ent.args))
    elseif ent.kind=="CARTESIAN_POINT"
        return _deref_point(entities,arg)
    end
    return nothing
end

function _step_straight_edge(entities,ent,depth,curvekinds)
    ent===nothing && (push!(curvekinds,"(missing curve)"); return false)
    ent.kind=="LINE" && return true
    if depth<8 && ent.kind in _STEP_BREP_WRAP_KINDS
        return _step_straight_edge(entities,
            _deref_entity(entities,_first_ref(ent.args)),depth+1,curvekinds)
    end
    push!(curvekinds,ent.kind)
    return false
end

function _step_oriented_edge(entities,oe,curvekinds,caller)
    ecref=nothing; eflag=true
    for a in oe.args
        a isa Tuple && a[1]===:ref && (ecref=a)
        flag=_brep_enum(a)
        flag===nothing || (eflag=flag)
    end
    ec=_deref_entity(entities,ecref)
    (ec!==nothing && ec.kind=="EDGE_CURVE") || throw(ArgumentError(
        "$caller: ORIENTED_EDGE #$(oe.id) does not reference an EDGE_CURVE"))
    length(ec.args)>=4 || throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) is truncated"))
    v1=_step_vertex_point(entities,ec.args[2])
    v2=_step_vertex_point(entities,ec.args[3])
    (v1===nothing || v2===nothing) && throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) has unresolvable vertex geometry"))
    curve=_deref_entity(entities,ec.args[4])
    straight=_step_straight_edge(entities,curve,0,curvekinds)
    return (eflag ? (v1,v2) : (v2,v1)),straight
end

function _chain_ring(edges,caller)
    length(edges)>=3 || throw(ArgumentError(
        "$caller: face boundary loop has fewer than three edges"))
    for i in 1:length(edges)-1
        edges[i][2]==edges[i+1][1] || throw(ArgumentError(
            "$caller: face boundary edges do not form an ordered closed chain"))
    end
    edges[end][2]==edges[1][1] || throw(ArgumentError(
        "$caller: face boundary loop does not close"))
    return NTuple{3,Float64}[e[1] for e in edges]
end

# Returns (ring, straight) or (nothing, false) for an unrecognized loop kind.
function _step_loop_ring(entities,loop,curvekinds,caller)
    if loop.kind=="POLY_LOOP"
        i=_skip_name(loop.args)
        length(loop.args)>=i || throw(ArgumentError(
            "$caller: POLY_LOOP #$(loop.id) is truncated"))
        refs=loop.args[i]
        refs isa Vector || throw(ArgumentError(
            "$caller: POLY_LOOP #$(loop.id) vertex list is malformed"))
        pts=NTuple{3,Float64}[]
        for r in refs
            p=_step_vertex_point(entities,r)
            p===nothing && throw(ArgumentError(
                "$caller: POLY_LOOP #$(loop.id) has unresolvable vertex geometry"))
            push!(pts,p)
        end
        length(pts)>1 && first(pts)==last(pts) && pop!(pts)
        length(pts)>=3 || throw(ArgumentError(
            "$caller: POLY_LOOP #$(loop.id) has fewer than three distinct vertices"))
        return pts,true
    elseif loop.kind=="EDGE_LOOP"
        i=_skip_name(loop.args)
        length(loop.args)>=i || throw(ArgumentError(
            "$caller: EDGE_LOOP #$(loop.id) is truncated"))
        refs=loop.args[i]
        refs isa Vector || throw(ArgumentError(
            "$caller: EDGE_LOOP #$(loop.id) member list is malformed"))
        edges=Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
        straight=true
        for r in refs
            oe=_deref_entity(entities,r)
            (oe!==nothing && oe.kind=="ORIENTED_EDGE") || throw(ArgumentError(
                "$caller: EDGE_LOOP #$(loop.id) member is not an ORIENTED_EDGE"))
            e,st=_step_oriented_edge(entities,oe,curvekinds,caller)
            straight &= st
            push!(edges,e)
        end
        return _chain_ring(edges,caller),straight
    end
    push!(curvekinds,loop.kind)
    return nothing,false
end

# Returns (rings or nothing, surface entity or nothing, sense, straight).
function _step_face_parts(entities,fe,curvekinds,caller)
    i=_skip_name(fe.args)
    length(fe.args)>=i || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has no boundary list"))
    bounds=fe.args[i]
    bounds isa Vector || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) boundary list is malformed"))
    isempty(bounds) && throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has no face boundaries"))
    surf=nothing; sense=true
    if fe.kind!="FACE"
        length(fe.args)>=i+2 || throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) is missing its surface or sense flag"))
        surf=_deref_entity(entities,fe.args[i+1])
        surf===nothing && throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) references a missing surface"))
        flag=_brep_enum(fe.args[i+2])
        flag===nothing || (sense=flag)
    end
    rings=Vector{NTuple{3,Float64}}[]
    straight=true
    for bref in bounds
        be=_deref_entity(entities,bref)
        (be!==nothing && be.kind in _STEP_BREP_BOUND_KINDS) || throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) boundary is not a FACE_BOUND/FACE_OUTER_BOUND"))
        length(be.args)>=3 || throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) is truncated"))
        loop=_deref_entity(entities,be.args[2])
        loop===nothing && throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) references a missing loop"))
        ring,st=_step_loop_ring(entities,loop,curvekinds,caller)
        straight &= st
        ring===nothing && return nothing,surf,sense,straight
        flag=_brep_enum(be.args[3])
        flag===nothing || flag || reverse!(ring)
        push!(rings,ring)
    end
    return rings,surf,sense,straight
end

function _newell_normal(ring::AbstractVector{<:NTuple{3}})
    nx=ny=nz=0.0; n=length(ring)
    for i in 1:n
        p=ring[i]; q=ring[i==n ? 1 : i+1]
        nx+=(p[2]-q[2])*(p[3]+q[3])
        ny+=(p[3]-q[3])*(p[1]+q[1])
        nz+=(p[1]-q[1])*(p[2]+q[2])
    end
    return (nx,ny,nz)
end

function _ring_scale(ring::AbstractVector{<:NTuple{3}})
    p0=ring[1]; s=0.0
    for p in ring
        d=_brep_vsub(p,p0)
        s=max(s,_brep_vdot(d,d))
    end
    return sqrt(s)
end

function _check_coplanar(ring,nrm,p0,tol,fid,caller)
    for p in ring
        d=abs(_brep_vdot(_brep_vsub(p,p0),nrm))
        (isfinite(d) && d<=tol) || throw(ArgumentError(
            "$caller: face #$fid is not planar (vertex off its plane by $d)"))
    end
    return nothing
end

# Coplanarity check for a bare polygonal FACE (no surface record).
function _check_poly_face(ring,fid,caller)
    nv=_newell_normal(ring)
    nn=sqrt(_brep_vdot(nv,nv))
    scale=_ring_scale(ring)
    (isfinite(nn) && nn>1e-12*max(scale*scale,1.0)) || throw(ArgumentError(
        "$caller: face #$fid boundary is degenerate (zero area)"))
    nrm=(nv[1]/nn,nv[2]/nn,nv[3]/nn)
    _check_coplanar(ring,nrm,ring[1],1e-9*max(scale,1.0),fid,caller)
    return nothing
end

function _check_plane_face(entities,surf,ring,sense,fid,caller)
    plcref=_first_ref(surf.args)
    plcref===nothing && throw(ArgumentError(
        "$caller: PLANE #$(surf.id) has no placement"))
    p0=_deref_point(entities,plcref)
    p0===nothing && throw(ArgumentError(
        "$caller: PLANE #$(surf.id) has no location point"))
    z=_deref_direction(entities,plcref)
    z===nothing && throw(ArgumentError(
        "$caller: PLANE #$(surf.id) has no axis direction"))
    nrm=_brep_vunit(z,caller)
    scale=_ring_scale(ring)
    _check_coplanar(ring,nrm,p0,1e-9*max(scale,1.0),fid,caller)
    nv=_newell_normal(ring)
    (isfinite(_brep_vdot(nv,nv)) && _brep_vdot(nv,nv)>0) || throw(ArgumentError(
        "$caller: face #$fid boundary is degenerate (zero area)"))
    _brep_vdot(nv,nrm)*(sense ? 1.0 : -1.0)>0 || throw(ArgumentError(
        "$caller: face #$fid boundary winding contradicts its PLANE normal/sense flag"))
    return nothing
end

_ring_degenerate(r)=all(p->p==first(r),r)
_rings_all_degenerate(rings)=isempty(rings) || all(_ring_degenerate,rings)

function _step_torus_params(entities,surf,caller)
    length(surf.args)>=3 || throw(ArgumentError(
        "$caller: TOROIDAL_SURFACE #$(surf.id) is truncated"))
    R=Float64(surf.args[end-1]); r=Float64(surf.args[end])
    (isfinite(R) && isfinite(r) && R>0 && r>0) || throw(ArgumentError(
        "$caller: TOROIDAL_SURFACE #$(surf.id) radii must be positive"))
    r<R || throw(ArgumentError(
        "$caller: TOROIDAL_SURFACE #$(surf.id) is not a ring torus (minor radius ≥ major radius)"))
    center=(0.0,0.0,0.0); z=(0.0,0.0,1.0); x=nothing
    plcref=_first_ref(surf.args)
    if plcref!==nothing
        plc=_deref_entity(entities,plcref)
        if plc!==nothing && plc.kind=="AXIS2_PLACEMENT_3D"
            p=_deref_point(entities,_first_ref(plc.args))
            p===nothing || (center=p)
            refs=[a for a in plc.args if a isa Tuple && a[1]===:ref]
            if length(refs)>=2
                d=_deref_direction(entities,refs[2]); d===nothing || (z=d)
            end
            if length(refs)>=3
                x=_deref_direction(entities,refs[3])
            end
        end
    end
    all(isfinite,center) || throw(ArgumentError(
        "$caller: TOROIDAL_SURFACE #$(surf.id) has a non-finite center"))
    return center,z,x,R,r
end

function _as_torus(entities)
    for ent in values(entities)
        ent.kind=="TOROIDAL_SURFACE" || continue
        return _step_torus_params(entities,ent,"import_step")
    end
    return nothing
end

function _step_torus_face(entities,surf,fill,caller)
    center,z,x,R,r=_step_torus_params(entities,surf,caller)
    surface=_torus_surface(center,z,x,R,r,caller)
    return fill ? _filled(surface,caller) : surface
end

# Pole-free periodic torus boundary: u is the revolution angle, v the meridian
# angle; the parametric rectangle wraps so seam vertices are shared by index.
function _torus_surface(center,Z,X,R,r,caller)
    z=_brep_vunit(Z,caller)
    x=X===nothing ? _brep_default_x(z,caller) : _brep_ortho(X,z,caller)
    y=_brep_vcross(z,x)
    nu=_BREP_TORUS_NU; nv=_BREP_TORUS_NV
    nn=_checked_brep_mul(nu,nv,caller,"torus node count")
    nn<=_BREP_MAX_SURFACE_NODES || throw(ArgumentError(
        "$caller: torus boundary exceeds the surface node limit"))
    nt=_checked_brep_mul(2,nn,caller,"torus triangle count")
    nt<=_BREP_MAX_SURFACE_TRIS || throw(ArgumentError(
        "$caller: torus boundary exceeds the surface triangle limit"))
    C=Matrix{Float64}(undef,3,nn)
    tp=2π
    for i in 0:nu-1
        u=tp*i/nu; cu=cos(u); su=sin(u)
        for j in 0:nv-1
            v=tp*j/nv; cv=cos(v); sv=sin(v)
            rho=R+r*cv
            idx=i*nv+j+1
            C[1,idx]=center[1]+rho*(cu*x[1]+su*y[1])+r*sv*z[1]
            C[2,idx]=center[2]+rho*(cu*x[2]+su*y[2])+r*sv*z[2]
            C[3,idx]=center[3]+rho*(cu*x[3]+su*y[3])+r*sv*z[3]
        end
    end
    T=Matrix{Int32}(undef,3,nt); k=0
    for i in 0:nu-1, j in 0:nv-1
        a=Int32(i*nv+j+1);        b=Int32(((i+1)%nu)*nv+j+1)
        c=Int32(((i+1)%nu)*nv+((j+1)%nv)+1); d=Int32(i*nv+((j+1)%nv)+1)
        k+=1; T[:,k]=[a;b;c]
        k+=1; T[:,k]=[a;c;d]
    end
    surface=Mesh(C;tris=T)
    diag=validate(surface)
    diag.ok || throw(ErrorException(
        "$caller: torus boundary mesh is invalid — "*join(diag.messages,"; ")))
    return surface
end

# ── polygonal shell assembly ──────────────────────────────────────────────────

function _poly_scale(pts)
    isempty(pts) && return 0.0
    x0=y0=z0=Inf; x1=y1=z1=-Inf
    for p in pts
        x0=min(x0,p[1]); y0=min(y0,p[2]); z0=min(z0,p[3])
        x1=max(x1,p[1]); y1=max(y1,p[2]); z1=max(z1,p[3])
    end
    all(isfinite,(x0,y0,z0,x1,y1,z1)) || return 0.0
    dx=x1-x0; dy=y1-y0; dz=z1-z0
    return sqrt(dx*dx+dy*dy+dz*dz)
end

# Snap-deduplicate ring vertices: sorted clustering merges only points within
# a relative 1e-9 tolerance, so independently-written identical corners merge
# while genuinely distinct vertices stay separate.
function _dedup_ring_nodes(rings,caller)
    pts=NTuple{3,Float64}[]
    for f in eachindex(rings), i in eachindex(rings[f])
        push!(pts,rings[f][i])
    end
    isempty(pts) && throw(ArgumentError("$caller: polyhedral shell has no vertices"))
    tol=1e-9*max(_poly_scale(pts),1.0)
    order=sortperm(pts)
    canon=Vector{Int32}(undef,length(pts))
    nodes=NTuple{3,Float64}[]
    lastpt=nothing
    for oi in order
        p=pts[oi]
        if lastpt===nothing ||
          !(abs(p[1]-lastpt[1])<=tol && abs(p[2]-lastpt[2])<=tol && abs(p[3]-lastpt[3])<=tol)
            push!(nodes,p); lastpt=p
        end
        canon[oi]=Int32(length(nodes))
    end
    length(nodes)<=_BREP_MAX_SURFACE_NODES || throw(ArgumentError(
        "$caller: polyhedral shell exceeds the surface node limit"))
    irings=Vector{Int32}[]
    k=0
    for f in eachindex(rings)
        out=Int32[]; prev=Int32(0)
        for _ in eachindex(rings[f])
            k+=1; v=canon[k]
            v!=prev && push!(out,v); prev=v
        end
        length(out)>1 && first(out)==last(out) && pop!(out)
        push!(irings,out)
    end
    return nodes,irings
end

function _fan_signed_volume(ring,nodes)
    p0=nodes[ring[1]]; s=0.0
    for i in 2:length(ring)-1
        s+=_brep_vdot(_brep_vcross(p0,nodes[ring[i]]),nodes[ring[i+1]])
    end
    return s/6
end

function _signed_area2(xs,ys)
    s=0.0; n=length(xs)
    for i in 1:n
        j=i==n ? 1 : i+1
        s+=xs[i]*ys[j]-xs[j]*ys[i]
    end
    return s/2
end

# Ear-clipping triangulation of a planar polygon. The working list is kept in
# CCW order in the dominant-axis projection; emitted triangles are flipped back
# when the stored ring was clockwise so they match the ring's 3-D orientation.
function _ear_clip!(tris,ring,nodes,caller)
    n=length(ring)
    nv=_newell_normal_idx(ring,nodes)
    (all(isfinite,nv) && _brep_vdot(nv,nv)>0) || throw(ArgumentError(
        "$caller: face polygon is degenerate (zero or non-finite area)"))
    ax=abs(nv[1]); ay=abs(nv[2]); az=abs(nv[3])
    i1,i2=ax>=ay ? (ax>=az ? (2,3) : (1,2)) : (ay>=az ? (1,3) : (1,2))
    xs=[nodes[v][i1] for v in ring]
    ys=[nodes[v][i2] for v in ring]
    W=collect(1:n)
    rev=_signed_area2(xs,ys)<0
    rev && reverse!(W)
    while length(W)>3
        t=_find_ear(W,xs,ys,true)
        if t==0
            t=_find_ear(W,xs,ys,false)
            t!=0 && (deleteat!(W,t); continue)
        end
        t==0 && throw(ArgumentError(
            "$caller: face boundary is not a simple polygon"))
        m=length(W)
        ip=W[t==1 ? m : t-1]; ic=W[t]; inx=W[t==m ? 1 : t+1]
        rev ? push!(tris,(ring[ip],ring[inx],ring[ic])) :
              push!(tris,(ring[ip],ring[ic],ring[inx]))
        deleteat!(W,t)
    end
    a,b,c=ring[W[1]],ring[W[2]],ring[W[3]]
    rev ? push!(tris,(a,c,b)) : push!(tris,(a,b,c))
    return tris
end

function _newell_normal_idx(ring,nodes)
    nx=ny=nz=0.0; n=length(ring)
    for i in 1:n
        p=nodes[ring[i]]; q=nodes[ring[i==n ? 1 : i+1]]
        nx+=(p[2]-q[2])*(p[3]+q[3])
        ny+=(p[3]-q[3])*(p[1]+q[1])
        nz+=(p[1]-q[1])*(p[2]+q[2])
    end
    return (nx,ny,nz)
end

# convex=true: find a convex ear (emits a triangle); convex=false: find a
# straight-line vertex safe to remove without emitting.
function _find_ear(W,xs,ys,convex)
    m=length(W)
    for t in 1:m
        ip=W[t==1 ? m : t-1]; ic=W[t]; inx=W[t==m ? 1 : t+1]
        ax=xs[ip]; ay=ys[ip]; bx=xs[ic]; by=ys[ic]; cx=xs[inx]; cy=ys[inx]
        cross=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
        if convex
            cross>0 || continue
        else
            cross==0 || continue
            dabx=cx-ax; daby=cy-ay
            dot=(bx-ax)*dabx+(by-ay)*daby
            (0<=dot<=dabx*dabx+daby*daby) || continue
        end
        inside=false
        for u in W
            (u==ip||u==ic||u==inx) && continue
            x=xs[u]; y=ys[u]
            c1=(bx-ax)*(y-ay)-(by-ay)*(x-ax)
            c2=(cx-bx)*(y-by)-(cy-by)*(x-bx)
            c3=(ax-cx)*(y-cy)-(ay-cy)*(x-cx)
            if c1>0 && c2>0 && c3>0
                inside=true; break
            end
        end
        inside || return t
    end
    return 0
end

# Assemble a closed, consistently oriented polygonal shell into a surface Mesh.
# rings: one outer boundary polygon per face (points in traversal order);
# roles: +1 outer shell, -1 void shell, one per face.
function _polyhedral_surface(rings,roles,caller)
    nfaces=length(rings)
    nfaces>=4 || throw(ArgumentError(
        "$caller: a closed polyhedral shell needs at least four faces"))
    nodes,irings=_dedup_ring_nodes(rings,caller)
    edge_dir=Dict{Tuple{Int32,Int32},Int}()
    edge_faces=Dict{Tuple{Int32,Int32},Vector{Int}}()
    face_edges=[Tuple{Int32,Int32}[] for _ in 1:nfaces]
    for f in 1:nfaces
        r=irings[f]; n=length(r)
        n>=3 || throw(ArgumentError(
            "$caller: face $f has fewer than three distinct vertices"))
        for i in 1:n
            a=r[i]; b=r[i==n ? 1 : i+1]
            a==b && throw(ArgumentError("$caller: face $f contains a zero-length edge"))
            key=a<b ? (a,b) : (b,a)
            edge_dir[(a,b)]=get(edge_dir,(a,b),0)+1
            fl=get!(edge_faces,key,Int[])
            f in fl || push!(fl,f)
            push!(face_edges[f],(a,b))
        end
    end
    for key in sort!(collect(keys(edge_faces)))
        ab=get(edge_dir,(key[1],key[2]),0)
        ba=get(edge_dir,(key[2],key[1]),0)
        (length(edge_faces[key])==2 && ab+ba==2) || throw(ArgumentError(
            "$caller: polyhedral shell is not closed — edge $(key) belongs to "*
            "$(length(edge_faces[key])) face(s)"))
        (ab==1 && ba==1) || throw(ArgumentError(
            "$caller: polyhedral shell faces are inconsistently oriented across edge $(key)"))
    end
    comp=fill(0,nfaces)
    ncomp=_brep_assign_components!(comp,face_edges,edge_faces,nfaces)
    scale=max(_poly_scale(nodes),1.0)
    vtol=1e-12*scale*scale*scale
    total=0.0
    for c in 1:ncomp
        members=[f for f in 1:nfaces if comp[f]==c]
        role=roles[members[1]]
        all(f->roles[f]==role,members) || throw(ArgumentError(
            "$caller: shell component $c mixes outer and void faces"))
        V=sum(_fan_signed_volume(irings[f],nodes) for f in members)
        abs(V)>vtol || throw(ArgumentError(
            "$caller: polyhedral shell component $c has zero volume"))
        if (role>0)==(V<0)
            for f in members
                reverse!(irings[f])
            end
            V=-V
        end
        total+=role>0 ? V : -abs(V)
    end
    total>vtol || throw(ArgumentError(
        "$caller: polyhedral shells do not bound a positive volume"))
    tris=NTuple{3,Int32}[]
    for f in 1:nfaces
        _ear_clip!(tris,irings[f],nodes,caller)
    end
    length(tris)<=_BREP_MAX_SURFACE_TRIS || throw(ArgumentError(
        "$caller: polyhedral shell exceeds the surface triangle limit"))
    C=Matrix{Float64}(undef,3,length(nodes))
    for (i,p) in enumerate(nodes); C[:,i]=[p...]; end
    T=Matrix{Int32}(undef,3,length(tris))
    for (i,t) in enumerate(tris); T[:,i]=[t...]; end
    surface=Mesh(C;tris=T)
    diag=validate(surface)
    diag.ok || throw(ErrorException(
        "$caller: assembled polyhedral surface is invalid — "*join(diag.messages,"; ")))
    return surface
end

function _brep_assign_components!(comp,face_edges,edge_faces,nfaces)
    nc=0; queue=Int[]
    for seed in 1:nfaces
        comp[seed]!=0 && continue
        nc+=1; comp[seed]=nc; empty!(queue); push!(queue,seed)
        while !isempty(queue)
            f=popfirst!(queue)
            for (a,b) in face_edges[f]
                key=a<b ? (a,b) : (b,a)
                for g in edge_faces[key]
                    g!=f && comp[g]==0 && (comp[g]=nc; push!(queue,g))
                end
            end
        end
    end
    return nc
end

function _step_brep_solid(brep::StepEntity,entities,fill)
    caller="import_step"
    length(brep.args)>=2 || throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shell"))
    shells=Tuple{Int,StepEntity}[]
    for (si,sarg) in enumerate(brep.args[2:end])
        items=sarg isa Vector ? sarg : Any[sarg]
        for item in items
            se=_deref_entity(entities,item)
            se===nothing && throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references a missing shell"))
            se.kind=="OPEN_SHELL" && throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references OPEN_SHELL #$(se.id) — "*
                "the shell is not closed"))
            se.kind=="CLOSED_SHELL" || throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references $(se.kind) #$(se.id) "*
                "where a CLOSED_SHELL was expected"))
            push!(shells,(si==1 ? 1 : -1,se))
        end
    end
    isempty(shells) && throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shells"))
    curvekinds=Set{String}(); seen=Set{String}()
    faces=NamedTuple{(:role,:rings,:surf,:sense,:straight,:id),
                     Tuple{Int,Union{Nothing,Vector{Vector{NTuple{3,Float64}}}},
                           Union{Nothing,StepEntity},Bool,Bool,Int}}[]
    for (role,sh) in shells
        i=_skip_name(sh.args)
        length(sh.args)>=i || throw(ArgumentError(
            "$caller: CLOSED_SHELL #$(sh.id) has no face list"))
        flist=sh.args[i]
        flist isa Vector || throw(ArgumentError(
            "$caller: CLOSED_SHELL #$(sh.id) face list is malformed"))
        isempty(flist) && throw(ArgumentError("$caller: CLOSED_SHELL #$(sh.id) is empty"))
        for fref in flist
            fe=_deref_entity(entities,fref)
            fe===nothing && throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(sh.id) references a missing face"))
            fe.kind in _STEP_BREP_FACE_KINDS || throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(sh.id) member $(fe.kind) #$(fe.id) is not a face"))
            rings,surf,sense,straight=_step_face_parts(entities,fe,curvekinds,caller)
            surf===nothing || push!(seen,surf.kind)
            push!(faces,(role=role,rings=rings,surf=surf,sense=sense,
                         straight=straight,id=fe.id))
        end
    end
    planar=!isempty(faces); ntorus=0; torus_face=nothing
    for face in faces
        face.rings===nothing && (planar=false; continue)
        length(face.rings)>1 && throw(ArgumentError(
            "$caller: face #$(face.id) has $(length(face.rings)) boundary loops — "*
            "inner loops (holes) are not supported"))
        ring=face.rings[1]
        if face.surf===nothing
            _check_poly_face(ring,face.id,caller)
        elseif face.surf.kind=="PLANE"
            _check_plane_face(entities,face.surf,ring,face.sense,face.id,caller)
        elseif face.surf.kind=="TOROIDAL_SURFACE"
            planar=false; ntorus+=1; torus_face=face; continue
        else
            planar=false; push!(curvekinds,"face surface $(face.surf.kind)"); continue
        end
        face.straight || (planar=false)
    end
    if planar
        rings=[only(f.rings) for f in faces]
        surface=_polyhedral_surface(rings,[f.role for f in faces],caller)
        return fill ? _filled(surface,caller) : surface
    end
    if ntorus==1 && length(faces)==1 && torus_face.rings!==nothing &&
       _rings_all_degenerate(torus_face.rings)
        return _step_torus_face(entities,torus_face.surf,fill,caller)
    end
    ntorus==1 && length(faces)==1 && throw(ArgumentError(
        "$caller: trimmed TOROIDAL_SURFACE faces are not supported (full torus only)"))
    kinds=sort!(collect(union(seen,curvekinds)))
    isempty(kinds) && push!(kinds,"unrecognized topology")
    throw(ArgumentError(
        "$caller: unsupported BREP face or edge topology; saw $(join(kinds,", "))"))
end

# ── IGES ──────────────────────────────────────────────────────────────────────

"""
    import_iges(path; fill=true) -> Mesh

Import one classified IGES solid: a type 150 Block, type 158 Sphere, the
supported type 156 Cylinder/Cone layouts, a type 160/198 torus, a full-torus
type 120 Surface of Revolution (full-circle generatrix in an axial plane,
full 2π sweep — the OCC export layout), or a
closed-shell polyhedron bounded by trimmed planar faces over line edges —
either through type 144 trimmed surfaces on type 108/190 planes with
142/141→102/110 boundaries, or through the type 186/514/510/508/504/502
manifold-solid chain. `fill=true` returns a validated tetrahedral mesh and
`fill=false` returns the closed triangle boundary. Multiple recognized solids,
nonplanar/nonlinear face topology, and files without a recognized solid are
explicit blockers. Use [`import_nurbs_iges`](@ref) for IGES 126/128 entities.
"""
function import_iges(path::AbstractString; fill::Bool=true)
    isfile(path) || throw(ArgumentError("import_iges: missing file $path"))
    lines=read(path,String)
    types,params,bad,transforms=_iges_entity_records(lines)
    surface=_iges_structured_solid(types,params,bad,transforms,"import_iges")
    if surface!==nothing
        nprim=count(rec->(!isempty(rec) && isinteger(rec[1]) &&
                          Int(rec[1]) in (150,156,158)),values(params))
        nprim==0 || throw(ArgumentError(
            "import_iges: multiple recognized solids are not supported in one import"))
        return fill ? _filled(surface,"import_iges") : surface
    end
    records=_iges_records(lines)
    isempty(records) && throw(ArgumentError("import_iges: no parameter records"))
    matches=Mesh[]
    for rec in records
        type=_as_int(rec[1],"import_iges","entity type")
        if type==150 && length(rec)>=7
            L,W,H=rec[2],rec[3],rec[4]
            x,y,z=rec[5],rec[6],rec[7]
            (L>0 && W>0 && H>0) || throw(ArgumentError("import_iges: Block extents must be positive"))
            push!(matches,box_surface(x,x+L,y,y+W,z,z+H))
        elseif type==158 && length(rec)>=5
            r=rec[2]; x,y,z=rec[3],rec[4],rec[5]
            r>0 || throw(ArgumentError("import_iges: Sphere radius must be positive"))
            push!(matches,sphere_surface((x,y,z),r))
        elseif type==156 && length(rec)>=10
            h,r1,r2=rec[2],rec[3],rec[4]
            x,y,z=rec[5],rec[6],rec[7]
            zi,zj,zk=rec[8],rec[9],rec[10]
            (h>0 && (r1>0 || r2>0) && r1>=0 && r2>=0) || throw(ArgumentError(
                "import_iges: Cone height must be positive and at least one radius must be positive"))
            push!(matches,cone_surface((x,y,z),(zi,zj,zk),r1,r2,h))
        elseif type==156 && length(rec)>=8
            r,h=rec[2],rec[3]; x,y,z=rec[4],rec[5],rec[6]
            zi,zj,zk=rec[7],rec[8], length(rec)>=9 ? rec[9] : 1.0
            (r>0 && h>0) || throw(ArgumentError("import_iges: Cylinder radius/height must be positive"))
            push!(matches,cylinder_surface((x,y,z),(zi,zj,zk),r,h))
        end
    end
    length(matches)<=1 || throw(ArgumentError(
        "import_iges: multiple recognized solids are not supported in one import"))
    isempty(matches) || return fill ? _filled(only(matches),"import_iges") : only(matches)
    types=sort!(unique(_as_int(rec[1],"import_iges","entity type") for rec in records))
    throw(ArgumentError(
        "import_iges: no supported solid (150 Block, 158 Sphere, 156 Cylinder/Cone); saw types $(types). " *
        "NURBS curves/surfaces use import_nurbs_iges"))
end

function _iges_records(source::AbstractString)
    isascii(source) || throw(ArgumentError("IGES import: input must be ASCII"))
    records=Vector{Vector{Float64}}()
    buf=IOBuffer()
    for raw in split(source, r"\r?\n")
        line=length(raw)>=80 ? raw[1:80] : rpad(raw,80)
        section=line[73]
        section=='P' || continue
        # IGES parameter data occupy columns 1:64. Columns 65:72 are the
        # directory-entry pointer and must never be parsed as parameters.
        body=rstrip(line[1:64])
        print(buf, body)
        if occursin(';', body)
            text=String(take!(buf))
            terminator=findfirst(==(';'),text)
            terminator===nothing && throw(ArgumentError("IGES import: missing record terminator"))
            tail=text[nextind(text,terminator):end]
            isempty(strip(tail)) || throw(ArgumentError(
                "IGES import: unexpected data after parameter-record terminator"))
            payload=text[firstindex(text):prevind(text,terminator)]
            pieces=split(payload, ','; keepempty=true)
            vals=Float64[]
            for piece in pieces
                token=strip(piece)
                isempty(token) && throw(ArgumentError("IGES import: empty numeric parameter"))
                normalized=replace(token,'D'=>'E','d'=>'e')
                v=tryparse(Float64,normalized)
                v===nothing && throw(ArgumentError(
                    "IGES import: invalid numeric parameter $(repr(token))"))
                isfinite(v) || throw(ArgumentError(
                    "IGES import: non-finite numeric parameter $(repr(token))"))
                push!(vals,v)
            end
            isempty(vals) || push!(records,vals)
        end
    end
    position(buf)==0 || throw(ArgumentError("IGES import: unterminated parameter record"))
    return records
end

function _as_int(value, caller, name)
    value isa Real || throw(ArgumentError("$caller: $name must be a real integer"))
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    v=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    (isfinite(v) && isinteger(v)) || throw(ArgumentError(
        "$caller: $name must be a finite integer"))
    (typemin(Int)<=v<=typemax(Int)) || throw(ArgumentError(
        "$caller: $name exceeds the platform Int range"))
    return Int(v)
end

# ── IGES directory-entry-linked entity parse ─────────────────────────────────

# Unlike `_iges_records` (flat, strict) this parser links each parameter record
# to its directory-entry sequence number so entity references can be resolved,
# and captures each entity's DE transform pointer (field 7, columns 49-56) so a
# referenced type 124 matrix can be applied when classifying analytic geometry.
# Records with non-numeric parameters (e.g. Hollerith names) are tolerated and
# marked in `bad`; structural malformation still throws.
function _iges_entity_records(source::AbstractString)
    isascii(source) || throw(ArgumentError("IGES import: input must be ASCII"))
    types=Dict{Int,Int}(); params=Dict{Int,Vector{Float64}}(); bad=Set{Int}()
    transforms=Dict{Int,Int}()
    buf=IOBuffer()
    for raw in split(source, r"\r?\n")
        line=length(raw)>=80 ? raw[1:80] : rpad(raw,80)
        section=line[73]
        if section=='D'
            seq=_iges_int_field(line[74:80]); seq===nothing && continue
            isodd(seq) || continue
            t=_iges_int_field(line[1:8]); t===nothing && continue
            types[seq]=t
            xf=_iges_int_field(line[49:56])
            xf!==nothing && xf!=0 && (transforms[seq]=xf)
        elseif section=='P'
            body=rstrip(line[1:64])
            print(buf,body)
            occursin(';',body) || continue
            text=String(take!(buf))
            terminator=findfirst(==(';'),text)
            terminator===nothing && throw(ArgumentError(
                "IGES import: missing record terminator"))
            tail=text[nextind(text,terminator):end]
            isempty(strip(tail)) || throw(ArgumentError(
                "IGES import: unexpected data after parameter-record terminator"))
            payload=text[firstindex(text):prevind(text,terminator)]
            vals=_iges_param_values(payload)
            de=_iges_int_field(line[65:72])
            if de===nothing
                # record whose owning DE cannot be identified: ignore
            elseif vals===nothing
                push!(bad,de)
            else
                params[de]=vals
            end
        end
    end
    position(buf)==0 || throw(ArgumentError(
        "IGES import: unterminated parameter record"))
    return types,params,bad,transforms
end

_iges_int_field(s)=tryparse(Int,strip(s))

function _iges_param_values(payload)
    vals=Float64[]
    for piece in split(payload,',';keepempty=true)
        token=strip(piece)
        isempty(token) && return nothing
        normalized=replace(token,'D'=>'E','d'=>'e')
        v=tryparse(Float64,normalized)
        (v===nothing || !isfinite(v)) && return nothing
        push!(vals,v)
    end
    return vals
end

function _iges_rec(params,bad,de,caller)
    de in bad && throw(ArgumentError(
        "$caller: IGES entity #$de has unsupported parameter data"))
    rec=get(params,de,nothing)
    rec===nothing && throw(ArgumentError(
        "$caller: IGES entity #$de is missing or unparseable"))
    return rec
end

function _iges_param_int(rec,i,caller,name)
    i<=length(rec) || throw(ArgumentError("$caller: $name — record is truncated"))
    return _as_int(rec[i],caller,name)
end

# ── IGES structured solids ───────────────────────────────────────────────────

function _iges_structured_solid(types,params,bad,transforms,caller)
    msbo=sort!([de for (de,t) in types if t==186])
    trimmed=sort!([de for (de,t) in types if t==144])
    bounded=sort!([de for (de,t) in types if t==143])
    tori=sort!([de for (de,t) in types if t==160 || t==198])
    hasbrep=!isempty(msbo)||!isempty(trimmed)||!isempty(bounded)
    hasbrep || isempty(tori) || return _iges_torus_solid(types,params,bad,tori,caller)
    hasbrep || return nothing
    length(msbo)>1 && throw(ArgumentError(
        "$caller: multiple recognized solids are not supported in one import"))
    isempty(tori) || throw(ArgumentError(
        "$caller: multiple recognized solids are not supported in one import"))
    return _iges_brep_surface(types,params,bad,transforms,msbo,trimmed,bounded,caller)
end

function _iges_torus_solid(types,params,bad,tori,caller)
    length(tori)>1 && throw(ArgumentError(
        "$caller: multiple recognized solids are not supported in one import"))
    return _iges_torus_surface(types,params,bad,only(tori),caller)
end

function _iges502_point(types,params,bad,vde,idx,caller)
    get(types,vde,0)==502 || throw(ArgumentError(
        "$caller: edge references #$vde which is not a Vertex List (type 502)"))
    rec=_iges_rec(params,bad,vde,caller)
    idx>=1 || throw(ArgumentError("$caller: vertex list index must be positive"))
    # rec[1] is the entity type, rec[2] the vertex count — vertex idx starts
    # at 3+3*(idx-1).
    base=3+3*(idx-1)
    length(rec)>=base+2 || throw(ArgumentError(
        "$caller: vertex list #$vde is truncated at index $idx"))
    return (rec[base],rec[base+1],rec[base+2])
end

_pt_close(a,b,tol)=abs(a[1]-b[1])<=tol && abs(a[2]-b[2])<=tol && abs(a[3]-b[3])<=tol

# Directed segments → vertex ring (unclosed; last point != first). IGES
# composites are expected to be ordered; a greedy reconnect tolerates writers
# that list members out of order or with flipped storage.
function _chain_segments(segs,tol,caller)
    length(segs)>=3 || throw(ArgumentError(
        "$caller: boundary loop has fewer than three edges"))
    if _segments_chained(segs,tol)
        _pt_close(segs[end][2],segs[1][1],tol) || throw(ArgumentError(
            "$caller: boundary loop does not close"))
        return NTuple{3,Float64}[s[1] for s in segs]
    end
    remaining=collect(segs)
    s0=popfirst!(remaining)
    chain=NTuple{3,Float64}[s0[1]]
    lastend=s0[2]
    while !isempty(remaining)
        found=false
        for (i,(a,b)) in enumerate(remaining)
            if _pt_close(a,lastend,tol)
                push!(chain,a); lastend=b; deleteat!(remaining,i); found=true; break
            end
        end
        found && continue
        for (i,(a,b)) in enumerate(remaining)
            if _pt_close(b,lastend,tol)
                push!(chain,b); lastend=a; deleteat!(remaining,i); found=true; break
            end
        end
        found || throw(ArgumentError(
            "$caller: boundary edges do not form a closed chain"))
    end
    _pt_close(lastend,chain[1],tol) || throw(ArgumentError(
        "$caller: boundary loop does not close"))
    return chain
end

function _segments_chained(segs,tol)
    for i in 1:length(segs)-1
        _pt_close(segs[i][2],segs[i+1][1],tol) || return false
    end
    return true
end

# Model-space boundary curves: type 110 lines and type 102 composites of them.
# Returns Vector{(p1,p2)} or nothing (recording the curve type in curvekinds).
function _iges_model_segments(types,params,bad,de,curvekinds,caller)
    t=get(types,de,0)
    rec=_iges_rec(params,bad,de,caller)
    if t==110
        length(rec)>=7 || throw(ArgumentError("$caller: IGES line #$de is truncated"))
        return Tuple{NTuple{3,Float64},NTuple{3,Float64}}[
            ((rec[2],rec[3],rec[4]),(rec[5],rec[6],rec[7]))]
    elseif t==102
        n=_iges_param_int(rec,2,caller,"composite curve member count")
        length(rec)>=2+n || throw(ArgumentError(
            "$caller: composite curve #$de is truncated"))
        segs=Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
        for k in 1:n
            part=_iges_model_segments(types,params,bad,
                       _iges_param_int(rec,2+k,caller,"composite curve member"),
                       curvekinds,caller)
            part===nothing && return nothing
            append!(segs,part)
        end
        return segs
    end
    push!(curvekinds,t)
    return nothing
end

# Resolve a 510 Loop (edge-list topology) to a vertex ring.
function _iges_loop508(types,params,bad,de,curvekinds,caller)
    get(types,de,0)==508 || throw(ArgumentError(
        "$caller: face references #$de which is not a Loop (type 508)"))
    rec=_iges_rec(params,bad,de,caller)
    n=_iges_param_int(rec,2,caller,"loop edge count")
    edges=Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
    i=3
    for _ in 1:n
        length(rec)>=i+4 || throw(ArgumentError("$caller: loop #$de is truncated"))
        typ=_iges_param_int(rec,i,caller,"loop entry type")
        ede=_iges_param_int(rec,i+1,caller,"loop edge list")
        ndx=_iges_param_int(rec,i+2,caller,"loop edge index")
        of=_iges_param_int(rec,i+3,caller,"loop edge orientation")
        npc=_iges_param_int(rec,i+4,caller,"loop parameter-curve count")
        npc>=0 || throw(ArgumentError("$caller: loop #$de has a negative parameter count"))
        i+=5+2*npc
        length(rec)>=i-1 || throw(ArgumentError("$caller: loop #$de is truncated"))
        typ==0 || throw(ArgumentError(
            "$caller: loop #$de uses non-edge entries (type $typ) — unsupported"))
        get(types,ede,0)==504 || throw(ArgumentError(
            "$caller: loop #$de edge #$ede is not an Edge List (type 504)"))
        erec=_iges_rec(params,bad,ede,caller)
        ndx>=1 || throw(ArgumentError("$caller: loop #$de edge index must be positive"))
        # erec[1] is the entity type, erec[2] the edge count — edge ndx starts
        # at 3+5*(ndx-1) (CURV,SVP,SVI,TVP,TVI).
        base=3+5*(ndx-1)
        length(erec)>=base+4 || throw(ArgumentError(
            "$caller: edge list #$ede is truncated at edge $ndx"))
        cde=_iges_param_int(erec,base,caller,"edge curve")
        ct=get(types,cde,0)
        ct==110 || (push!(curvekinds,ct); return nothing)
        svp=_iges_param_int(erec,base+1,caller,"edge start vertex list")
        svi=_iges_param_int(erec,base+2,caller,"edge start vertex index")
        tvp=_iges_param_int(erec,base+3,caller,"edge end vertex list")
        tvi=_iges_param_int(erec,base+4,caller,"edge end vertex index")
        sv=_iges502_point(types,params,bad,svp,svi,caller)
        tv=_iges502_point(types,params,bad,tvp,tvi,caller)
        push!(edges,of==0 ? (tv,sv) : (sv,tv))
    end
    isempty(edges) && throw(ArgumentError("$caller: loop #$de is empty"))
    tol=1e-9*max(_poly_scale(vcat([s[1] for s in edges],[s[2] for s in edges])),1.0)
    return _chain_segments(edges,tol,caller)
end

# Resolve a face boundary entity (142 curve-on-surface, 141 boundary, 102/110
# model curves, or a 508 loop) to a vertex ring, or nothing for unsupported
# curve kinds (recorded in curvekinds).
function _iges_boundary(types,params,bad,de,curvekinds,caller)
    t=get(types,de,0)
    rec=_iges_rec(params,bad,de,caller)
    if t==142
        length(rec)>=5 || throw(ArgumentError(
            "$caller: curve-on-surface #$de is truncated"))
        bp=_iges_param_int(rec,4,caller,"curve-on-surface model curve")
        cp=_iges_param_int(rec,5,caller,"curve-on-surface parameter curve")
        tgt=bp!=0 ? bp : cp
        tgt!=0 || throw(ArgumentError(
            "$caller: curve-on-surface #$de has no boundary curve"))
        return _iges_boundary(types,params,bad,tgt,curvekinds,caller)
    elseif t==141
        length(rec)>=5 || throw(ArgumentError("$caller: boundary #$de is truncated"))
        btype=_iges_param_int(rec,2,caller,"boundary type")
        btype in (1,2) || throw(ArgumentError(
            "$caller: boundary #$de uses only parameter-space curves"))
        nb=_iges_param_int(rec,5,caller,"boundary curve count")
        segs=Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
        i=6
        for _ in 1:nb
            length(rec)>=i+2 || throw(ArgumentError("$caller: boundary #$de is truncated"))
            mde=_iges_param_int(rec,i,caller,"boundary model curve")
            mflag=_iges_param_int(rec,i+1,caller,"boundary curve orientation")
            npc=_iges_param_int(rec,i+2,caller,"boundary parameter-curve count")
            npc>=0 || throw(ArgumentError(
                "$caller: boundary #$de has a negative parameter count"))
            i+=3+npc
            length(rec)>=i-1 || throw(ArgumentError("$caller: boundary #$de is truncated"))
            part=_iges_model_segments(types,params,bad,mde,curvekinds,caller)
            part===nothing && return nothing
            if mflag==2
                part=reverse!([(b,a) for (a,b) in part])
            elseif mflag!=1
                throw(ArgumentError(
                    "$caller: boundary #$de uses curve orientation flag $mflag — unsupported"))
            end
            append!(segs,part)
        end
        isempty(segs) && throw(ArgumentError("$caller: boundary #$de is empty"))
        tol=1e-9*max(_poly_scale(vcat([s[1] for s in segs],[s[2] for s in segs])),1.0)
        return _chain_segments(segs,tol,caller)
    elseif t==102 || t==110
        segs=_iges_model_segments(types,params,bad,de,curvekinds,caller)
        segs===nothing && return nothing
        tol=1e-9*max(_poly_scale(vcat([s[1] for s in segs],[s[2] for s in segs])),1.0)
        return _chain_segments(segs,tol,caller)
    elseif t==508
        return _iges_loop508(types,params,bad,de,curvekinds,caller)
    end
    push!(curvekinds,t)
    return nothing
end

# Face over type 144 Trimmed Parametric Surface.
# Returns (rings, surface_type, surface_de); rings is empty for an untrimmed
# natural-boundary face (only meaningful on a torus) or unresolvable curves.
function _iges_face144(types,params,bad,de,curvekinds,caller)
    rec=_iges_rec(params,bad,de,caller)
    length(rec)>=5 || throw(ArgumentError(
        "$caller: trimmed surface #$de is truncated"))
    surf_de=_iges_param_int(rec,2,caller,"trimmed surface base")
    n1=_iges_param_int(rec,3,caller,"outer-boundary flag")
    n2=_iges_param_int(rec,4,caller,"inner-boundary count")
    n2>=0 || throw(ArgumentError("$caller: trimmed surface #$de has a negative inner count"))
    length(rec)>=5+n2 || throw(ArgumentError("$caller: trimmed surface #$de is truncated"))
    n2==0 || throw(ArgumentError(
        "$caller: trimmed surface #$de has $n2 inner boundaries — face holes are not supported"))
    stype=get(types,surf_de,0)
    stype!=0 || throw(ArgumentError(
        "$caller: trimmed surface #$de references missing surface #$surf_de"))
    n1==0 || n1==1 || throw(ArgumentError(
        "$caller: trimmed surface #$de outer-boundary flag $n1 is invalid"))
    if n1==0
        stype==198 || throw(ArgumentError(
            "$caller: trimmed surface #$de on type $stype has no outer boundary"))
        return Vector{NTuple{3,Float64}}[],stype,surf_de
    end
    ptro=_iges_param_int(rec,5,caller,"outer boundary")
    ptro!=0 || throw(ArgumentError(
        "$caller: trimmed surface #$de has no outer boundary"))
    ring=_iges_boundary(types,params,bad,ptro,curvekinds,caller)
    rings=ring===nothing ? Vector{NTuple{3,Float64}}[] : Vector{NTuple{3,Float64}}[ring]
    return rings,stype,surf_de
end

# Face over type 510 Face (manifold-solid boundary).
# Returns (rings, surface_type, surface_de); rings is empty for unresolvable
# curve kinds.
function _iges_face510(types,params,bad,de,curvekinds,caller)
    rec=_iges_rec(params,bad,de,caller)
    length(rec)>=4 || throw(ArgumentError("$caller: face #$de is truncated"))
    surf_de=_iges_param_int(rec,2,caller,"face surface")
    n=_iges_param_int(rec,3,caller,"face loop count")
    _iges_param_int(rec,4,caller,"face outer-loop flag")
    n>=1 || throw(ArgumentError("$caller: face #$de has no boundary loops"))
    n==1 || throw(ArgumentError(
        "$caller: face #$de has $n loops — inner loops (holes) are not supported"))
    length(rec)>=4+n || throw(ArgumentError("$caller: face #$de is truncated"))
    stype=get(types,surf_de,0)
    stype!=0 || throw(ArgumentError(
        "$caller: face #$de references missing surface #$surf_de"))
    ring=_iges_loop508(types,params,bad,
        _iges_param_int(rec,5,caller,"face loop"),curvekinds,caller)
    rings=ring===nothing ? Vector{NTuple{3,Float64}}[] : Vector{NTuple{3,Float64}}[ring]
    return rings,stype,surf_de
end

# Plane parameters for type 108 (implicit plane) and 190 (plane surface).
function _iges_plane_params(types,params,bad,surf_de,caller)
    st=get(types,surf_de,0)
    rec=_iges_rec(params,bad,surf_de,caller)
    if st==108
        length(rec)>=5 || throw(ArgumentError("$caller: plane #$surf_de is truncated"))
        nrm=(rec[2],rec[3],rec[4]); d=rec[5]
        un=_brep_vunit(nrm,caller)
        all(isfinite,(d,)) || throw(ArgumentError(
            "$caller: plane #$surf_de has a non-finite offset"))
        return un,(un[1]*d,un[2]*d,un[3]*d)
    elseif st==190
        length(rec)>=3 || throw(ArgumentError(
            "$caller: plane surface #$surf_de is truncated"))
        pde=_iges_param_int(rec,2,caller,"plane surface point")
        nde=_iges_param_int(rec,3,caller,"plane surface normal")
        get(types,pde,0)==116 || throw(ArgumentError(
            "$caller: plane surface #$surf_de point #$pde is not a Point (type 116)"))
        get(types,nde,0)==123 || throw(ArgumentError(
            "$caller: plane surface #$surf_de normal #$nde is not a Direction (type 123)"))
        prec=_iges_rec(params,bad,pde,caller)
        nrec=_iges_rec(params,bad,nde,caller)
        length(prec)>=4 || throw(ArgumentError("$caller: point #$pde is truncated"))
        length(nrec)>=4 || throw(ArgumentError("$caller: direction #$nde is truncated"))
        p0=(prec[2],prec[3],prec[4]); nrm=(nrec[2],nrec[3],nrec[4])
        all(isfinite,p0) || throw(ArgumentError("$caller: point #$pde is non-finite"))
        return _brep_vunit(nrm,caller),p0
    end
    return nothing
end

# Returns (ring, orientation flip) pairs grouped per face, plus surface types.
function _iges_brep_faces(types,params,bad,msbo,trimmed,bounded,caller)
    facerecs=Tuple{Vector{NTuple{3,Float64}},Int,Int}[]  # (ring, role, surf_de)
    flips=Bool[]
    curvekinds=Set{Int}(); seen=Set{Int}()
    for de in msbo
        rec=_iges_rec(params,bad,de,caller)
        nv=_iges_param_int(rec,4,caller,"MSBO void-shell count")
        nv>=0 || throw(ArgumentError("$caller: MSBO #$de has a negative void count"))
        length(rec)>=4+2*nv || throw(ArgumentError("$caller: MSBO #$de is truncated"))
        shells=Tuple{Int,Int,Int}[]
        push!(shells,(1,_iges_param_int(rec,2,caller,"MSBO shell"),
                      _iges_param_int(rec,3,caller,"MSBO shell orientation")))
        for k in 1:nv
            push!(shells,(-1,_iges_param_int(rec,4+2k-1,caller,"MSBO void shell"),
                          _iges_param_int(rec,4+2k,caller,"MSBO void orientation")))
        end
        for (role,sde,sflag) in shells
            get(types,sde,0)==514 || throw(ArgumentError(
                "$caller: MSBO #$de references #$sde which is not a Shell (type 514)"))
            srec=_iges_rec(params,bad,sde,caller)
            nf=_iges_param_int(srec,2,caller,"shell face count")
            nf>=1 || throw(ArgumentError("$caller: shell #$sde is empty"))
            length(srec)>=2+2*nf || throw(ArgumentError("$caller: shell #$sde is truncated"))
            for k in 1:nf
                fde=_iges_param_int(srec,2+2k-1,caller,"shell face")
                fof=_iges_param_int(srec,2+2k,caller,"shell face orientation")
                get(types,fde,0)==510 || throw(ArgumentError(
                    "$caller: shell #$sde member #$fde is not a Face (type 510)"))
                rings,stype,surf_de=_iges_face510(types,params,bad,fde,curvekinds,caller)
                push!(seen,stype)
                flip=(fof==0)!=(sflag==0)
                for r in rings
                    flip && reverse!(r)
                    push!(facerecs,(r,role,surf_de))
                end
                isempty(rings) && push!(facerecs,(NTuple{3,Float64}[],role,surf_de))
            end
        end
    end
    for de in trimmed
        rings,stype,surf_de=_iges_face144(types,params,bad,de,curvekinds,caller)
        push!(seen,stype)
        if isempty(rings)
            push!(facerecs,(NTuple{3,Float64}[],1,surf_de))
        else
            for r in rings
                push!(facerecs,(r,1,surf_de))
            end
        end
    end
    for de in bounded
        rec=_iges_rec(params,bad,de,caller)
        length(rec)>=4 || throw(ArgumentError("$caller: bounded surface #$de is truncated"))
        _iges_param_int(rec,2,caller,"bounded surface type")
        surf_de=_iges_param_int(rec,3,caller,"bounded surface base")
        nb=_iges_param_int(rec,4,caller,"bounded surface boundary count")
        nb>=1 || throw(ArgumentError("$caller: bounded surface #$de has no boundaries"))
        nb==1 || throw(ArgumentError(
            "$caller: bounded surface #$de has $nb boundaries — face holes are not supported"))
        length(rec)>=4+nb || throw(ArgumentError("$caller: bounded surface #$de is truncated"))
        stype=get(types,surf_de,0)
        stype!=0 || throw(ArgumentError(
            "$caller: bounded surface #$de references missing surface #$surf_de"))
        push!(seen,stype)
        ring=_iges_boundary(types,params,bad,
            _iges_param_int(rec,5,caller,"bounded surface boundary"),curvekinds,caller)
        push!(facerecs,(ring===nothing ? NTuple{3,Float64}[] : ring,1,surf_de))
    end
    return facerecs,curvekinds,seen
end

function _iges_brep_surface(types,params,bad,transforms,msbo,trimmed,bounded,caller)
    facerecs,curvekinds,seen=_iges_brep_faces(types,params,bad,msbo,trimmed,bounded,caller)
    planar=!isempty(facerecs); ntorus=0; torus_de=0
    rings=Vector{NTuple{3,Float64}}[]; roles=Int[]
    for (ring,role,surf_de) in facerecs
        stype=get(types,surf_de,0)
        if stype==198 || stype==120
            planar=false; ntorus+=1; torus_de=surf_de; continue
        end
        if stype==108 || stype==190
            isempty(ring) && (planar=false; continue)
            plane=_iges_plane_params(types,params,bad,surf_de,caller)
            nrm,p0=plane
            scale=_ring_scale(ring)
            _check_coplanar(ring,nrm,p0,1e-9*max(scale,1.0),surf_de,caller)
            push!(rings,ring); push!(roles,role); continue
        end
        planar=false
    end
    if planar
        surface=_polyhedral_surface(rings,roles,caller)
        return surface
    end
    if ntorus==1 && length(facerecs)==1
        ring=facerecs[1][1]
        if get(types,torus_de,0)==120
            return _iges_torus120(
                types,params,bad,transforms,torus_de,ring,caller)
        end
        (isempty(ring) || _ring_degenerate(ring)) || throw(ArgumentError(
            "$caller: trimmed type-198 torus faces are not supported (full torus only)"))
        return _iges_torus_surface(types,params,bad,torus_de,caller)
    end
    kinds=sort!(collect(union(seen,curvekinds)))
    isempty(kinds) && push!(kinds,0)
    throw(ArgumentError(
        "$caller: unsupported IGES BRep face or edge topology; saw entity types $(kinds)"))
end

function _iges_torus_surface(types,params,bad,de,caller)
    t=get(types,de,0)
    rec=_iges_rec(params,bad,de,caller)
    center=(0.0,0.0,0.0); z=(0.0,0.0,1.0); x=nothing
    if t==160
        length(rec)>=9 || throw(ArgumentError("$caller: torus #$de is truncated"))
        R=rec[2]; r=rec[3]
        center=(rec[4],rec[5],rec[6]); z=(rec[7],rec[8],rec[9])
    elseif t==198
        length(rec)>=5 || throw(ArgumentError("$caller: toroidal surface #$de is truncated"))
        pde=_iges_param_int(rec,2,caller,"torus center point")
        ade=_iges_param_int(rec,3,caller,"torus axis direction")
        get(types,pde,0)==116 || throw(ArgumentError(
            "$caller: toroidal surface #$de center #$pde is not a Point (type 116)"))
        get(types,ade,0)==123 || throw(ArgumentError(
            "$caller: toroidal surface #$de axis #$ade is not a Direction (type 123)"))
        prec=_iges_rec(params,bad,pde,caller)
        arec=_iges_rec(params,bad,ade,caller)
        length(prec)>=4 || throw(ArgumentError("$caller: point #$pde is truncated"))
        length(arec)>=4 || throw(ArgumentError("$caller: direction #$ade is truncated"))
        center=(prec[2],prec[3],prec[4]); z=(arec[2],arec[3],arec[4])
        R=rec[4]; r=rec[5]
        if length(rec)>=6
            rde=_iges_param_int(rec,6,caller,"torus reference direction")
            rde!=0 || (x=nothing)
            if rde!=0
                get(types,rde,0)==123 || throw(ArgumentError(
                    "$caller: toroidal surface #$de reference #$rde is not a Direction"))
                rrec=_iges_rec(params,bad,rde,caller)
                length(rrec)>=4 || throw(ArgumentError(
                    "$caller: direction #$rde is truncated"))
                x=(rrec[2],rrec[3],rrec[4])
            end
        end
    else
        throw(ArgumentError("$caller: entity #$de type $t is not a torus"))
    end
    (isfinite(R) && isfinite(r) && R>0 && r>0) || throw(ArgumentError(
        "$caller: torus #$de radii must be positive"))
    r<R || throw(ArgumentError(
        "$caller: torus #$de is not a ring torus (minor radius ≥ major radius)"))
    all(isfinite,center) || throw(ArgumentError("$caller: torus #$de has a non-finite center"))
    all(isfinite,z) || throw(ArgumentError("$caller: torus #$de has a non-finite axis"))
    return _torus_surface(center,z,x,R,r,caller)
end

# Rigid-motion transform (type 124) referenced by an entity's DE transform
# field → ((row1,row2,row3), translation), or nothing for the identity.
# Rotation columns must be orthonormal — anything else is not a rigid motion
# and would distort the classified circle.
function _iges_de_transform(types,params,bad,transforms,de,caller)
    tde=get(transforms,de,0)
    tde==0 && return nothing
    get(types,tde,0)==124 || throw(ArgumentError(
        "$caller: entity #$de references #$tde which is not a Transformation " *
        "Matrix (type 124)"))
    rec=_iges_rec(params,bad,tde,caller)
    length(rec)>=13 || throw(ArgumentError(
        "$caller: transformation matrix #$tde is truncated"))
    # Type 124 stores the translation interleaved: R11,R12,R13,T1, R21,...,T3.
    R=((rec[2],rec[3],rec[4]),(rec[6],rec[7],rec[8]),(rec[10],rec[11],rec[12]))
    t=(rec[5],rec[9],rec[13])
    vals=(R[1]...,R[2]...,R[3]...,t...)
    all(isfinite,vals) || throw(ArgumentError(
        "$caller: transformation matrix #$tde is non-finite"))
    for i in 1:3, j in i:3
        col_i=(R[1][i],R[2][i],R[3][i]); col_j=(R[1][j],R[2][j],R[3][j])
        want=i==j ? 1.0 : 0.0
        abs(_brep_vdot(col_i,col_j)-want)<=1e-6 || throw(ArgumentError(
            "$caller: transformation matrix #$tde is not a rigid motion"))
    end
    return R,t
end

# Apply (R,t) to a point (translate=true) or a direction (translate=false).
function _iges_xf_apply(xf,p,translate::Bool)
    xf===nothing && return p
    R,t=xf
    q=(R[1][1]*p[1]+R[1][2]*p[2]+R[1][3]*p[3],
       R[2][1]*p[1]+R[2][2]*p[2]+R[2][3]*p[3],
       R[3][1]*p[1]+R[3][2]*p[2]+R[3][3]*p[3])
    return translate ? (q[1]+t[1],q[2]+t[2],q[3]+t[3]) : q
end

# Type 120 Surface of Revolution → full ring torus. The generatrix must be a
# complete circular arc (type 100, start == end point) whose plane contains
# the revolution axis (a type 110 line), swept through a full 2π. `ring` is
# the resolved face boundary: OCC-style exports carry the (u,v)-param-box
# seam, whose coordinates sit off the surface and carry no trim — ignored.
# Any ring vertex landing ON the torus is a real 3-D trim → explicit blocker.
function _iges_torus120(types,params,bad,transforms,de,ring,caller)
    rec=_iges_rec(params,bad,de,caller)
    length(rec)>=5 || throw(ArgumentError(
        "$caller: surface of revolution #$de is truncated"))
    lde=_iges_param_int(rec,2,caller,"revolution axis line")
    gde=_iges_param_int(rec,3,caller,"revolution generatrix")
    sa,ea=rec[4],rec[5]
    (isfinite(sa) && isfinite(ea)) || throw(ArgumentError(
        "$caller: surface of revolution #$de has a non-finite sweep"))
    get(types,lde,0)==110 || throw(ArgumentError(
        "$caller: surface of revolution #$de axis #$lde is not a Line " *
        "(type 110)"))
    lrec=_iges_rec(params,bad,lde,caller)
    length(lrec)>=7 || throw(ArgumentError("$caller: line #$lde is truncated"))
    lxf=_iges_de_transform(types,params,bad,transforms,lde,caller)
    p1=_iges_xf_apply(lxf,(lrec[2],lrec[3],lrec[4]),true)
    p2=_iges_xf_apply(lxf,(lrec[5],lrec[6],lrec[7]),true)
    (all(isfinite,p1) && all(isfinite,p2)) || throw(ArgumentError(
        "$caller: line #$lde is non-finite"))
    axis=_brep_vunit(_brep_vsub(p2,p1),caller)
    span=ea-sa
    abs(span-2π)<=1e-6*max(1.0,abs(span)) || throw(ArgumentError(
        "$caller: surface of revolution #$de is not a full 2π sweep"))
    get(types,gde,0)==100 || throw(ArgumentError(
        "$caller: surface of revolution #$de generatrix #$gde is not a " *
        "Circular Arc (type 100)"))
    grec=_iges_rec(params,bad,gde,caller)
    length(grec)>=8 || throw(ArgumentError("$caller: arc #$gde is truncated"))
    zl,xc,yc,x1,y1,x2,y2=grec[2],grec[3],grec[4],
                          grec[5],grec[6],grec[7],grec[8]
    all(isfinite,(zl,xc,yc,x1,y1,x2,y2)) || throw(ArgumentError(
        "$caller: arc #$gde is non-finite"))
    rl=hypot(x1-xc,y1-yc)
    rl>0 || throw(ArgumentError("$caller: arc #$gde has a zero radius"))
    hypot(x1-x2,y1-y2)<=1e-9*max(1.0,rl) || throw(ArgumentError(
        "$caller: revolution #$de generatrix #$gde is not a full circle"))
    gxf=_iges_de_transform(types,params,bad,transforms,gde,caller)
    C=_iges_xf_apply(gxf,(xc,yc,zl),true)
    n=_iges_xf_apply(gxf,(0.0,0.0,1.0),false)
    all(isfinite,C) || throw(ArgumentError("$caller: arc #$gde is non-finite"))
    n=_brep_vunit(n,caller)
    tol=1e-7*max(1.0,rl,sqrt(_brep_vdot(C,C)))
    # A torus generatrix circle lies in a plane through the revolution axis:
    # the axis direction is parallel to the circle plane and the axis itself
    # satisfies the plane equation.
    abs(_brep_vdot(axis,n))<=1e-7 || throw(ArgumentError(
        "$caller: revolution #$de axis is not parallel to the generatrix " *
        "plane — not a torus"))
    abs(_brep_vdot(_brep_vsub(p1,C),n))<=tol || throw(ArgumentError(
        "$caller: revolution #$de axis does not lie in the generatrix " *
        "plane — not a torus"))
    rel=_brep_vsub(C,p1)
    axial=_brep_vdot(rel,axis)
    radial=_brep_vsub(rel,(axis[1]*axial,axis[2]*axial,axis[3]*axial))
    R=sqrt(_brep_vdot(radial,radial))
    (isfinite(R) && R>rl) || throw(ArgumentError(
        "$caller: revolution #$de generatrix reaches the axis — " *
        "not a ring torus"))
    center=(p1[1]+axis[1]*axial,p1[2]+axis[2]*axial,p1[3]+axis[3]*axial)
    x=_brep_vunit(radial,caller)
    r2=rl*rl
    for p in ring
        d=_brep_vsub(p,center)
        av=_brep_vdot(d,axis)
        rv=_brep_vsub(d,(axis[1]*av,axis[2]*av,axis[3]*av))
        rho=sqrt(_brep_vdot(rv,rv))
        abs((rho-R)^2+av*av-r2)<=1e-7*max(1.0,R*R,r2) && throw(ArgumentError(
            "$caller: trimmed surface-of-revolution faces are not " *
            "supported (full torus only)"))
    end
    return _torus_surface(center,axis,x,R,rl,caller)
end

function _checked_brep_add(args...)
    length(args)>=4 || throw(ArgumentError("BRep: internal checked-add contract"))
    caller=args[end-1]; name=args[end]
    total=0
    try
        for value in args[1:end-2]
            value isa Int || throw(ArgumentError("$caller: $name is not an Int count"))
            total=Base.checked_add(total,value)
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("$caller: $name overflows the platform Int range"))
    end
    return total
end

function _checked_brep_mul(a::Int,b::Int,caller,name)
    try
        return Base.checked_mul(a,b)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name overflows the platform Int range"))
    end
end

function _expand_knots(mults, knots, caller)
    mults isa Vector && knots isa Vector || throw(ArgumentError(
        "$caller: knot multiplicities and knots must be lists"))
    length(mults)==length(knots) || throw(ArgumentError(
        "$caller: knot multiplicity count mismatch"))
    expanded=Tuple{Int,Float64}[]
    total=0
    for (m,k) in zip(mults,knots)
        mm=_as_int(m,caller,"knot multiplicity")
        mm>=1 || throw(ArgumentError("$caller: knot multiplicity must be ≥ 1"))
        kk=try Float64(k) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: knot must be Float64-representable"))
        end
        isfinite(kk) || throw(ArgumentError("$caller: non-finite knot"))
        total=try Base.checked_add(total,mm) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: expanded knot count overflows Int"))
        end
        total<=_MAX_BREP_KNOT_VALUES || throw(ArgumentError(
            "$caller: expanded knot count $total exceeds $_MAX_BREP_KNOT_VALUES"))
        push!(expanded,(mm,kk))
    end
    U=Float64[]; sizehint!(U,total)
    for (mm,kk) in expanded, _ in 1:mm
        push!(U,kk)
    end
    return U
end

function _ref_points(entities, arg)
    arg isa Vector || return nothing
    points=NTuple{3,Float64}[]
    for item in arg
        p=_deref_point(entities,item)
        p===nothing && return nothing
        push!(points,p)
    end
    return points
end

function _float_vec(arg, caller)
    arg isa Vector || throw(ArgumentError("$caller: expected a numeric list"))
    vals=Float64[]
    for item in arg
        item isa Real || throw(ArgumentError("$caller: expected numbers"))
        v=try Float64(item) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: number must be Float64-representable"))
        end
        isfinite(v) || throw(ArgumentError("$caller: non-finite number"))
        push!(vals,v)
    end
    return vals
end

function _skip_name(args)
    return !isempty(args) && args[1] isa AbstractString ? 2 : 1
end

function _bspline_curve_from_args(entities, degree_arg, points_arg, mults_arg, knots_arg,
                                  weights_arg, caller)
    degree=_as_int(degree_arg,caller,"degree")
    points=_ref_points(entities, points_arg)
    points===nothing && throw(ArgumentError("$caller: missing control points"))
    knots=_expand_knots(mults_arg, knots_arg, caller)
    weights=weights_arg===nothing ? nothing : _float_vec(weights_arg, caller)
    return NURBSCurve(degree, knots, points, weights)
end

function _complex_parts(ent::StepEntity, caller)
    parts=Dict{String,StepEntity}()
    for part in ent.args
        part isa StepEntity || continue
        haskey(parts,part.kind) && throw(ArgumentError(
            "$caller: duplicate complex component $(part.kind)"))
        parts[part.kind]=part
    end
    return parts
end

function _step_nurbs_curve(ent::StepEntity, entities)
    if ent.kind=="B_SPLINE_CURVE_WITH_KNOTS"
        i=_skip_name(ent.args)
        length(ent.args)>=i+6 || return nothing
        return _bspline_curve_from_args(entities, ent.args[i], ent.args[i+1],
                                        ent.args[i+5], ent.args[i+6], nothing,
                                        "import_nurbs_step")
    elseif ent.kind=="COMPLEX"
        parts=_complex_parts(ent,"import_nurbs_step")
        haskey(parts,"B_SPLINE_CURVE") || return nothing
        haskey(parts,"B_SPLINE_CURVE_WITH_KNOTS") || return nothing
        bc=parts["B_SPLINE_CURVE"]; bk=parts["B_SPLINE_CURVE_WITH_KNOTS"]
        length(bc.args)>=2 && length(bk.args)>=2 || return nothing
        weights=if haskey(parts,"RATIONAL_B_SPLINE_CURVE") && !isempty(parts["RATIONAL_B_SPLINE_CURVE"].args)
            parts["RATIONAL_B_SPLINE_CURVE"].args[1]
        else
            nothing
        end
        return _bspline_curve_from_args(entities, bc.args[1], bc.args[2],
                                        bk.args[1], bk.args[2], weights,
                                        "import_nurbs_step")
    end
    return nothing
end

function _nested_ref_points(entities, arg)
    arg isa Vector || return nothing
    rows=Vector{NTuple{3,Float64}}[]
    for row in arg
        pts=_ref_points(entities,row)
        pts===nothing && return nothing
        push!(rows,pts)
    end
    isempty(rows) && return nothing
    nv=length(rows[1])
    all(length(r)==nv for r in rows) || return nothing
    nu=length(rows)
    C=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    for i in 1:nu, j in 1:nv
        C[i,j]=rows[i][j]
    end
    return C
end

function _nested_float_matrix(arg, caller)
    arg isa Vector || throw(ArgumentError("$caller: weights must be a nested list"))
    isempty(arg) && throw(ArgumentError("$caller: weight matrix must not be empty"))
    rows=Vector{Float64}[]
    for row in arg
        push!(rows,_float_vec(row,caller))
    end
    nv=length(rows[1])
    nv>0 || throw(ArgumentError("$caller: weight rows must not be empty"))
    all(length(row)==nv for row in rows) || throw(ArgumentError(
        "$caller: weight rows must have equal lengths"))
    nu=length(rows)
    W=Matrix{Float64}(undef,nu,nv)
    for i in 1:nu, j in 1:nv
        W[i,j]=rows[i][j]
    end
    return W
end

function _bspline_surface_from_args(entities,du_arg,dv_arg,points_arg,
                                    u_mults_arg,v_mults_arg,u_knots_arg,v_knots_arg,
                                    weights_arg,caller)
    du=_as_int(du_arg,caller,"degree_u")
    dv=_as_int(dv_arg,caller,"degree_v")
    C=_nested_ref_points(entities,points_arg)
    C===nothing && throw(ArgumentError("$caller: missing surface control points"))
    u_knots=_expand_knots(u_mults_arg,u_knots_arg,caller)
    v_knots=_expand_knots(v_mults_arg,v_knots_arg,caller)
    weights=weights_arg===nothing ? nothing : _nested_float_matrix(weights_arg,caller)
    return NURBSSurface(du,dv,u_knots,v_knots,C,weights)
end

function _step_nurbs_surface(ent::StepEntity, entities)
    if ent.kind=="B_SPLINE_SURFACE_WITH_KNOTS"
        i=_skip_name(ent.args)
        length(ent.args)>=i+10 || return nothing
        return _bspline_surface_from_args(entities,ent.args[i],ent.args[i+1],
            ent.args[i+2],ent.args[i+7],ent.args[i+8],ent.args[i+9],ent.args[i+10],
            nothing,"import_nurbs_step")
    elseif ent.kind=="COMPLEX"
        parts=_complex_parts(ent,"import_nurbs_step")
        haskey(parts,"B_SPLINE_SURFACE") || return nothing
        haskey(parts,"B_SPLINE_SURFACE_WITH_KNOTS") || return nothing
        bs=parts["B_SPLINE_SURFACE"]; bk=parts["B_SPLINE_SURFACE_WITH_KNOTS"]
        length(bs.args)>=3 && length(bk.args)>=4 || return nothing
        weights=if haskey(parts,"RATIONAL_B_SPLINE_SURFACE") &&
                   !isempty(parts["RATIONAL_B_SPLINE_SURFACE"].args)
            parts["RATIONAL_B_SPLINE_SURFACE"].args[1]
        else
            nothing
        end
        return _bspline_surface_from_args(entities,bs.args[1],bs.args[2],bs.args[3],
            bk.args[1],bk.args[2],bk.args[3],bk.args[4],weights,
            "import_nurbs_step")
    end
    return nothing
end

"""
    import_nurbs_step(path) -> Vector

Import supported STEP `B_SPLINE_CURVE_WITH_KNOTS` and
`B_SPLINE_SURFACE_WITH_KNOTS` entities as native [`NURBSCurve`](@ref) and
[`NURBSSurface`](@ref) objects, including supported complex rational curve and
surface forms. Returns every recognized object in STEP entity order and blocks
when none are present.
"""
function import_nurbs_step(path::AbstractString)
    isfile(path) || throw(ArgumentError("import_nurbs_step: missing file $path"))
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "import_nurbs_step: $path is not an ISO-10303-21 STEP file"))
    entities=parse_step_entities(source)
    objects=Any[]
    for ent in sort!(collect(values(entities)); by=e->e.id)
        c=_step_nurbs_curve(ent,entities)
        c===nothing || push!(objects,c)
        s=_step_nurbs_surface(ent,entities)
        s===nothing || push!(objects,s)
    end
    isempty(objects) && throw(ArgumentError(
        "import_nurbs_step: no B_SPLINE_CURVE_WITH_KNOTS or B_SPLINE_SURFACE_WITH_KNOTS"))
    return objects
end

function _iges_nurbs_curve(rec)
    _as_int(rec[1],"import_nurbs_iges","entity type")==126 || return nothing
    length(rec)>=8 || throw(ArgumentError("import_nurbs_iges: IGES 126 record is truncated"))
    K=_as_int(rec[2],"import_nurbs_iges","K")
    M=_as_int(rec[3],"import_nurbs_iges","M")
    (K>=1 && M>=1 && K>=M) || throw(ArgumentError(
        "import_nurbs_iges: require K ≥ M ≥ 1 for IGES 126"))
    n=_checked_brep_add(K,1,"import_nurbs_iges","control count")
    nknots=_checked_brep_add(K,M,2,"import_nurbs_iges","knot count")
    (n<=_MAX_BREP_KNOT_VALUES && nknots<=_MAX_BREP_KNOT_VALUES) || throw(ArgumentError(
        "import_nurbs_iges: IGES 126 count exceeds $_MAX_BREP_KNOT_VALUES"))
    i=8
    required=_checked_brep_add(i-1,nknots,n,
                               _checked_brep_mul(3,n,"import_nurbs_iges","control coordinate count"),
                               "import_nurbs_iges","IGES 126 record length")
    length(rec)>=required || throw(ArgumentError(
        "import_nurbs_iges: IGES 126 record is truncated"))
    knots=Float64[rec[i+k-1] for k in 1:nknots]; i+=nknots
    weights=Float64[rec[i+k-1] for k in 1:n]; i+=n
    controls=NTuple{3,Float64}[]
    for _ in 1:n
        push!(controls,(rec[i],rec[i+1],rec[i+2])); i+=3
    end
    return NURBSCurve(M,knots,controls,weights)
end

function _iges_nurbs_surface(rec)
    _as_int(rec[1],"import_nurbs_iges","entity type")==128 || return nothing
    length(rec)>=10 || throw(ArgumentError("import_nurbs_iges: IGES 128 record is truncated"))
    K1=_as_int(rec[2],"import_nurbs_iges","K1")
    K2=_as_int(rec[3],"import_nurbs_iges","K2")
    M1=_as_int(rec[4],"import_nurbs_iges","M1")
    M2=_as_int(rec[5],"import_nurbs_iges","M2")
    (K1>=1 && M1>=1 && K1>=M1 && K2>=1 && M2>=1 && K2>=M2) ||
        throw(ArgumentError("import_nurbs_iges: require K1 ≥ M1 ≥ 1 and K2 ≥ M2 ≥ 1"))
    nu=_checked_brep_add(K1,1,"import_nurbs_iges","u control count")
    nv=_checked_brep_add(K2,1,"import_nurbs_iges","v control count")
    nku=_checked_brep_add(K1,M1,2,"import_nurbs_iges","u knot count")
    nkv=_checked_brep_add(K2,M2,2,"import_nurbs_iges","v knot count")
    maximum((nu,nv,nku,nkv))<=_MAX_BREP_KNOT_VALUES || throw(ArgumentError(
        "import_nurbs_iges: IGES 128 count exceeds $_MAX_BREP_KNOT_VALUES"))
    ncontrols=_checked_brep_mul(nu,nv,"import_nurbs_iges","surface control count")
    ncoords=_checked_brep_mul(3,ncontrols,"import_nurbs_iges","control coordinate count")
    i=11
    required=_checked_brep_add(i-1,nku,nkv,ncontrols,ncoords,
                               "import_nurbs_iges","IGES 128 record length")
    length(rec)>=required || throw(ArgumentError(
        "import_nurbs_iges: IGES 128 record is truncated"))
    knots_u=Float64[rec[i+k-1] for k in 1:nku]; i+=nku
    knots_v=Float64[rec[i+k-1] for k in 1:nkv]; i+=nkv
    W=Matrix{Float64}(undef,nu,nv)
    C=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    for j in 1:nv, iu in 1:nu
        W[iu,j]=rec[i]; i+=1
    end
    for j in 1:nv, iu in 1:nu
        C[iu,j]=(rec[i],rec[i+1],rec[i+2]); i+=3
    end
    return NURBSSurface(M1,M2,knots_u,knots_v,C,W)
end

"""
    import_nurbs_iges(path) -> Vector

Import all supported IGES type 126 B-spline curves and type 128 tensor-product
B-spline surfaces as native NURBS objects. Record counts, finite numeric data,
and allocation bounds are checked before arrays are allocated.
"""
function import_nurbs_iges(path::AbstractString)
    isfile(path) || throw(ArgumentError("import_nurbs_iges: missing file $path"))
    records=_iges_records(read(path,String))
    isempty(records) && throw(ArgumentError("import_nurbs_iges: no parameter records"))
    objects=Any[]
    for rec in records
        isempty(rec) && continue
        c=_iges_nurbs_curve(rec)
        c===nothing || push!(objects,c)
        s=_iges_nurbs_surface(rec)
        s===nothing || push!(objects,s)
    end
    isempty(objects) && throw(ArgumentError(
        "import_nurbs_iges: no IGES 126/128 NURBS records"))
    return objects
end

function _iges126_fields(c::NURBSCurve)
    K=length(c.controls)-1
    polynomial=all(==(1.0), c.weights) ? 1.0 : 0.0
    rec=Float64[126,K,c.degree,0,0,polynomial,0]
    append!(rec,c.knots)
    append!(rec,c.weights)
    for p in c.controls
        append!(rec, (p[1],p[2],p[3]))
    end
    push!(rec,c.knots[c.degree+1],c.knots[end-c.degree],0.0,0.0,0.0)
    return rec
end

function _iges128_fields(s::NURBSSurface)
    nu,nv=size(s.controls)
    K1,K2=nu-1,nv-1
    polynomial=all(==(1.0), s.weights) ? 1.0 : 0.0
    rec=Float64[128,K1,K2,s.degree_u,s.degree_v,0,0,polynomial,0,0]
    append!(rec,s.knots_u)
    append!(rec,s.knots_v)
    for j in 1:nv, i in 1:nu
        push!(rec,s.weights[i,j])
    end
    for j in 1:nv, i in 1:nu
        p=s.controls[i,j]
        append!(rec,(p[1],p[2],p[3]))
    end
    push!(rec, s.knots_u[s.degree_u+1], s.knots_u[end-s.degree_u],
          s.knots_v[s.degree_v+1], s.knots_v[end-s.degree_v])
    return rec
end

function _iges_sequence(section::Char,sequence::Int)
    1<=sequence<=9_999_999 || throw(ArgumentError(
        "export_iges_nurbs: $section sequence $sequence exceeds seven digits"))
    return string(section,lpad(sequence,7,'0'))
end

function _iges_count_field(section::Char,count::Int)
    0<=count<=9_999_999 || throw(ArgumentError(
        "export_iges_nurbs: $section count $count exceeds seven digits"))
    return string(section,lpad(count,7))
end

function _iges_section_line(data::AbstractString,section::Char,sequence::Int)
    isascii(data) || throw(ArgumentError("export_iges_nurbs: IGES sections must be ASCII"))
    ncodeunits(data)<=72 || throw(ArgumentError(
        "export_iges_nurbs: internal $section-section line exceeds 72 columns"))
    return rpad(data,72)*_iges_sequence(section,sequence)
end

_iges_hollerith(value::AbstractString)=string(ncodeunits(value),'H',value)

function _iges_global_lines()
    # A stable compatibility epoch makes otherwise identical exports byte-for-byte
    # reproducible. It is metadata only; geometric parameter records carry no date.
    timestamp="20260824.000000"
    fields=["","",_iges_hollerith("Tessella"),_iges_hollerith("Tessella.iges"),
            _iges_hollerith("Tessella.jl"),_iges_hollerith("Tessella.jl"),
            "32","308","15","308","15","","1.","2",_iges_hollerith("MM"),
            "1","0.01",_iges_hollerith(timestamp),"1E-12","1E308","","",
            "11","0",_iges_hollerith(timestamp),""]
    payload=join(fields,',')*';'
    chunks=String[]
    start=firstindex(payload)
    while start<=lastindex(payload)
        stop=min(start+71,lastindex(payload))
        push!(chunks,payload[start:stop])
        start=stop+1
    end
    return [_iges_section_line(chunk,'G',i) for (i,chunk) in enumerate(chunks)]
end

function _iges_numeric_text(record::Vector{Float64})
    out=IOBuffer()
    for (i,value) in enumerate(record)
        isfinite(value) || throw(ArgumentError(
            "export_iges_nurbs: record contains a non-finite value"))
        i>1 && write(out,',')
        if isinteger(value) && abs(value)<1e12
            print(out,Int(value))
        else
            print(out,value)
        end
    end
    write(out,';')
    return String(take!(out))
end

function _iges_parameter_chunks(record::Vector{Float64})
    encoded=_iges_numeric_text(record)
    chunks=String[]
    start=firstindex(encoded)
    while start<=lastindex(encoded)
        stop=min(start+63,lastindex(encoded))
        push!(chunks,encoded[start:stop])
        start=stop+1
    end
    return chunks
end

function _iges_directory_lines(entity_type::Int,parameter_start::Int,
                               parameter_lines::Int,directory_sequence::Int,
                               status::String)
    length(status)==8 && all(isdigit,status) || throw(ArgumentError(
        "export_iges_nurbs: internal directory status must have eight digits"))
    first_fields=(entity_type,parameter_start,0,0,0,0,0,0)
    all(value->0<=value<=99_999_999,first_fields) || throw(ArgumentError(
        "export_iges_nurbs: directory field exceeds eight digits"))
    first=join(lpad.(string.(first_fields),8))*status*
          _iges_sequence('D',directory_sequence)
    second_fields=Any[entity_type,0,0,parameter_lines,0,nothing,nothing,nothing,0]
    second=join(value===nothing ? " "^8 : lpad(string(value),8)
                for value in second_fields)*_iges_sequence('D',directory_sequence+1)
    return first,second
end

function _iges_output_entities(records)
    entities=NamedTuple{(:record,:status),Tuple{Vector{Float64},String}}[]
    for record in records
        entity_type=_as_int(record[1],"export_iges_nurbs","entity type")
        entity_type in (126,128) || throw(ArgumentError(
            "export_iges_nurbs: internal unsupported entity type $entity_type"))
        if entity_type==128
            # An IGES 128 is a parametric surface definition. A 144 wrapper with
            # no trimming loops exposes its full rectangular parameter domain as
            # a topological surface to independent CAD readers.
            base_index=length(entities)+2
            base_directory_sequence=2base_index-1
            push!(entities,(record=Float64[144,base_directory_sequence,0,0,0],
                            status="00000000"))
            push!(entities,(record=record,status="00010000"))
        else
            push!(entities,(record=record,status="00000000"))
        end
        2length(entities)<=9_999_999 || throw(ArgumentError(
            "export_iges_nurbs: directory entry count exceeds seven digits"))
    end
    return entities
end

function _write_iges(path, records)
    target=abspath(path); parent=dirname(target)
    isdir(parent) || throw(ArgumentError(
        "export_iges_nurbs: parent directory does not exist: $parent"))
    entities=_iges_output_entities(records)
    chunks=[_iges_parameter_chunks(entity.record) for entity in entities]
    parameter_starts=Int[]; next_parameter=1
    for entity_chunks in chunks
        push!(parameter_starts,next_parameter)
        next_parameter=_checked_brep_add(next_parameter,length(entity_chunks),
                                          "export_iges_nurbs","parameter sequence")
        next_parameter-1<=9_999_999 || throw(ArgumentError(
            "export_iges_nurbs: parameter line count exceeds seven digits"))
    end
    global_lines=_iges_global_lines()
    mktemp(parent) do temporary,io
        println(io,_iges_section_line("Tessella.jl NURBS IGES export",'S',1))
        for line in global_lines
            println(io,line)
        end
        for i in eachindex(entities)
            entity_type=_as_int(entities[i].record[1],"export_iges_nurbs","entity type")
            directory_sequence=2i-1
            first,second=_iges_directory_lines(entity_type,parameter_starts[i],
                                                length(chunks[i]),directory_sequence,
                                                entities[i].status)
            println(io,first); println(io,second)
        end
        parameter_sequence=1
        for i in eachindex(entities),chunk in chunks[i]
            directory_sequence=2i-1
            pointer=lpad(string(lpad(directory_sequence,7,'0')),8)
            println(io,rpad(chunk,64)*pointer*_iges_sequence('P',parameter_sequence))
            parameter_sequence+=1
        end
        counts=_iges_count_field('S',1)*_iges_count_field('G',length(global_lines))*
               _iges_count_field('D',2length(entities))*
               _iges_count_field('P',parameter_sequence-1)
        println(io,_iges_section_line(counts,'T',1))
        flush(io); close(io)
        mv(temporary,target;force=true)
    end
    return path
end

"""
    export_iges_nurbs(path, objects) -> path

Validate and atomically write a nonempty collection of native NURBS curves and
surfaces as IGES type 126/128 parameter records. Each surface receives an
untrimmed type-144 wrapper so CAD readers expose it as a topological face.
Unsupported objects are rejected before the destination is replaced.
"""
function export_iges_nurbs(path::AbstractString, objects)
    items=try collect(objects) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("export_iges_nurbs: objects must be an iterable collection"))
    end
    isempty(items) && throw(ArgumentError("export_iges_nurbs: no NURBS objects"))
    records=Vector{Vector{Float64}}()
    for obj in items
        if obj isa NURBSCurve
            checked=NURBSCurve(obj.degree,obj.knots,obj.controls,obj.weights)
            push!(records,_iges126_fields(checked))
        elseif obj isa NURBSSurface
            checked=NURBSSurface(obj.degree_u,obj.degree_v,obj.knots_u,obj.knots_v,
                                 obj.controls,obj.weights)
            push!(records,_iges128_fields(checked))
        else
            throw(ArgumentError("export_iges_nurbs: expected NURBSCurve or NURBSSurface"))
        end
    end
    return _write_iges(path, records)
end

end # module
