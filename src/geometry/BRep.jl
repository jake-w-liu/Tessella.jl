"""
    BRep

Native ISO-10303-21 STEP and IGES CAD import. Solids that classify as an
axis-aligned block, sphere, right circular cylinder, or right circular cone
are converted to Tessella surfaces and filled. NURBS curves and surfaces
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
        while i<=last && (isspace(s[i]) || s[i]==','); i=nextind(s,i); end
        i<=last || throw(ArgumentError("import_step: unterminated list"))
        s[i]==')' && return items, nextind(s,i)
        if s[i]=='('
            sub,i=_step_parse_list(s,i); push!(items,sub)
        elseif s[i]=='\''
            j=nextind(s,i)
            while j<=last && s[j]!='\''; j=nextind(s,j); end
            j<=last || throw(ArgumentError("import_step: unterminated string"))
            push!(items,s[nextind(s,i):prevind(s,j)])
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
        else
            j=i
            while j<=last && (isdigit(s[j]) || s[j] in ('.','e','E','+','-')); j=nextind(s,j); end
            raw=s[i:prevind(s,j)]
            v=tryparse(Float64,raw)
            v===nothing && throw(ArgumentError("import_step: bad number $raw"))
            isfinite(v) || throw(ArgumentError("import_step: non-finite number $raw"))
            push!(items,v); i=j
        end
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
block, `SPHERE`, `RIGHT_CIRCULAR_CYLINDER`, or `RIGHT_CIRCULAR_CONE`. By default
the classified boundary is tetrahedralized; `fill=false` returns its validated
closed triangle surface. Multiple recognized solids, mixed point-cloud topology,
and files without a recognized solid are explicit blockers. Use
[`import_nurbs_step`](@ref) for STEP B-splines.
"""
function import_step(path::AbstractString; fill::Bool=true)
    isfile(path) || throw(ArgumentError("import_step: missing file $path"))
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "import_step: $path is not an ISO-10303-21 STEP file"))
    entities=parse_step_entities(source)
    points=_step_points(entities)
    primitive_count=count(ent->ent.kind in ("SPHERE","RIGHT_CIRCULAR_CYLINDER",
                                             "RIGHT_CIRCULAR_CONE"),
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
        "RIGHT_CIRCULAR_CYLINDER, or RIGHT_CIRCULAR_CONE); saw $(join(kinds, ", ")). " *
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

# ── IGES ──────────────────────────────────────────────────────────────────────

"""
    import_iges(path; fill=true) -> Mesh

Import one classified IGES primitive record: type 150 Block, type 158 Sphere,
or the supported type 156 Cylinder/Cone layouts. `fill=true` returns a validated
tetrahedral mesh and `fill=false` returns the closed triangle boundary. Multiple
recognized solids and files without a recognized solid are explicit blockers. Use
[`import_nurbs_iges`](@ref) for IGES 126/128 entities.
"""
function import_iges(path::AbstractString; fill::Bool=true)
    isfile(path) || throw(ArgumentError("import_iges: missing file $path"))
    lines=read(path,String)
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
