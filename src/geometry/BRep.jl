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
        id=parse(Int,data[idstart:prevind(data,j)])
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
            push!(items,(:ref, parse(Int,s[j:prevind(s,k)]))); i=k
        else
            j=i
            while j<=last && (isdigit(s[j]) || s[j] in ('.','e','E','+','-')); j=nextind(s,j); end
            raw=s[i:prevind(s,j)]
            v=tryparse(Float64,raw)
            v===nothing && throw(ArgumentError("import_step: bad number $raw"))
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

function import_step(path::AbstractString; fill::Bool=true)
    isfile(path) || throw(ArgumentError("import_step: missing file $path"))
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "import_step: $path is not an ISO-10303-21 STEP file"))
    entities=parse_step_entities(source)
    points=_step_points(entities)
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
    box=_as_box(points)
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

function import_iges(path::AbstractString; fill::Bool=true)
    isfile(path) || throw(ArgumentError("import_iges: missing file $path"))
    lines=read(path,String)
    records=_iges_records(lines)
    isempty(records) && throw(ArgumentError("import_iges: no parameter records"))
    for rec in records
        type=Int(rec[1])
        if type==150 && length(rec)>=7
            L,W,H=rec[2],rec[3],rec[4]
            x,y,z=rec[5],rec[6],rec[7]
            (L>0 && W>0 && H>0) || throw(ArgumentError("import_iges: Block extents must be positive"))
            surface=box_surface(x,x+L,y,y+W,z,z+H)
            return fill ? _filled(surface,"import_iges") : surface
        elseif type==158 && length(rec)>=5
            r=rec[2]; x,y,z=rec[3],rec[4],rec[5]
            r>0 || throw(ArgumentError("import_iges: Sphere radius must be positive"))
            surface=sphere_surface((x,y,z),r)
            return fill ? _filled(surface,"import_iges") : surface
        elseif type==156 && length(rec)>=10
            h,r1,r2=rec[2],rec[3],rec[4]
            x,y,z=rec[5],rec[6],rec[7]
            zi,zj,zk=rec[8],rec[9],rec[10]
            (h>0 && (r1>0 || r2>0) && r1>=0 && r2>=0) || throw(ArgumentError(
                "import_iges: Cone height must be positive and at least one radius must be positive"))
            surface=cone_surface((x,y,z),(zi,zj,zk),r1,r2,h)
            return fill ? _filled(surface,"import_iges") : surface
        elseif type==156 && length(rec)>=8
            r,h=rec[2],rec[3]; x,y,z=rec[4],rec[5],rec[6]
            zi,zj,zk=rec[7],rec[8], length(rec)>=9 ? rec[9] : 1.0
            (r>0 && h>0) || throw(ArgumentError("import_iges: Cylinder radius/height must be positive"))
            surface=cylinder_surface((x,y,z),(zi,zj,zk),r,h)
            return fill ? _filled(surface,"import_iges") : surface
        end
    end
    types=sort!(unique(Int(rec[1]) for rec in records))
    throw(ArgumentError(
        "import_iges: no supported solid (150 Block, 158 Sphere, 156 Cylinder/Cone); saw types $(types). " *
        "NURBS curves/surfaces use import_nurbs_iges"))
end

function _iges_records(source::AbstractString)
    records=Vector{Vector{Float64}}()
    buf=IOBuffer()
    current=0
    for raw in split(source, r"\r?\n")
        line=length(raw)>=80 ? raw[1:80] : rpad(raw,80)
        section=line[73]
        section=='P' || continue
        body=rstrip(line[1:72])
        print(buf, body)
        if occursin(';', body)
            text=String(take!(buf))
            text=replace(text, ';'=>"")
            pieces=split(text, ','; keepempty=false)
            vals=Float64[]
            for piece in pieces
                v=tryparse(Float64, strip(piece))
                v===nothing && continue
                push!(vals,v)
            end
            isempty(vals) || push!(records,vals)
        end
    end
    return records
end

function _as_int(value, caller, name)
    x=try Int(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be an integer"))
    end
    return x
end

function _expand_knots(mults, knots, caller)
    mults isa Vector && knots isa Vector || throw(ArgumentError(
        "$caller: knot multiplicities and knots must be lists"))
    length(mults)==length(knots) || throw(ArgumentError(
        "$caller: knot multiplicity count mismatch"))
    U=Float64[]
    for (m,k) in zip(mults,knots)
        mm=_as_int(m,caller,"knot multiplicity")
        mm>=1 || throw(ArgumentError("$caller: knot multiplicity must be ≥ 1"))
        kk=Float64(k)
        isfinite(kk) || throw(ArgumentError("$caller: non-finite knot"))
        for _ in 1:mm
            push!(U,kk)
        end
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
        v=Float64(item)
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

function _step_nurbs_curve(ent::StepEntity, entities)
    if ent.kind=="B_SPLINE_CURVE_WITH_KNOTS"
        i=_skip_name(ent.args)
        length(ent.args)>=i+6 || return nothing
        return _bspline_curve_from_args(entities, ent.args[i], ent.args[i+1],
                                        ent.args[i+5], ent.args[i+6], nothing,
                                        "import_nurbs_step")
    elseif ent.kind=="COMPLEX"
        parts=Dict{String,StepEntity}()
        for part in ent.args
            part isa StepEntity || continue
            parts[part.kind]=part
        end
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

function _step_nurbs_surface(ent::StepEntity, entities)
    ent.kind=="B_SPLINE_SURFACE_WITH_KNOTS" || return nothing
    i=_skip_name(ent.args)
    length(ent.args)>=i+10 || return nothing
    du=_as_int(ent.args[i],"import_nurbs_step","degree_u")
    dv=_as_int(ent.args[i+1],"import_nurbs_step","degree_v")
    C=_nested_ref_points(entities, ent.args[i+2])
    C===nothing && throw(ArgumentError("import_nurbs_step: missing surface control points"))
    u_knots=_expand_knots(ent.args[i+7], ent.args[i+9], "import_nurbs_step")
    v_knots=_expand_knots(ent.args[i+8], ent.args[i+10], "import_nurbs_step")
    return NURBSSurface(du,dv,u_knots,v_knots,C)
end

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
    Int(rec[1])==126 || return nothing
    length(rec)>=8 || throw(ArgumentError("import_nurbs_iges: IGES 126 record is truncated"))
    K=_as_int(rec[2],"import_nurbs_iges","K")
    M=_as_int(rec[3],"import_nurbs_iges","M")
    n=K+1
    nknots=K+M+2
    i=8
    length(rec)>=i+nknots+n+3n-1 || throw(ArgumentError(
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
    Int(rec[1])==128 || return nothing
    length(rec)>=10 || throw(ArgumentError("import_nurbs_iges: IGES 128 record is truncated"))
    K1=_as_int(rec[2],"import_nurbs_iges","K1")
    K2=_as_int(rec[3],"import_nurbs_iges","K2")
    M1=_as_int(rec[4],"import_nurbs_iges","M1")
    M2=_as_int(rec[5],"import_nurbs_iges","M2")
    nu,nv=K1+1,K2+1
    nku,nkv=K1+M1+2,K2+M2+2
    i=11
    length(rec)>=i+nku+nkv+nu*nv+3*nu*nv-1 || throw(ArgumentError(
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
    push!(rec, c.knots[c.degree+1], c.knots[end-c.degree])
    return rec
end

function _iges128_fields(s::NURBSSurface)
    nu,nv=size(s.controls)
    K1,K2=nu-1,nv-1
    polynomial=all(==(1.0), s.weights) ? 1.0 : 0.0
    rec=Float64[128,K1,K2,s.degree_u,s.degree_v,0,0,0,0,polynomial]
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

function _write_iges(path, records)
    open(path,"w") do io
        println(io, rpad("",72)*"S      1")
        println(io, rpad("1H,,1H;,8HTessella,11HTessella.jl,11HTessella.jl,32,8,15,11HTessella.jl,1.",72)*"G      1")
        seq=0
        for rec in records
            buf=IOBuffer()
            for (i,v) in enumerate(rec)
                i>1 && print(buf,',')
                if v==floor(v) && abs(v)<1e12
                    print(buf, Int(v))
                else
                    print(buf, v)
                end
            end
            print(buf,';')
            bytes=take!(buf)
            n=length(bytes)
            pos=1
            while pos<=n
                seq+=1
                stop=min(pos+63,n)
                write(io,view(bytes,pos:stop))
                print(io, " "^(72-(stop-pos+1)))
                println(io, "P", lpad(seq,7))
                pos=stop+1
            end
        end
        println(io, rpad("",72)*"T      1")
    end
    return path
end

function export_iges_nurbs(path::AbstractString, objects)
    isempty(objects) && throw(ArgumentError("export_iges_nurbs: no NURBS objects"))
    records=Vector{Vector{Float64}}()
    for obj in objects
        if obj isa NURBSCurve
            push!(records,_iges126_fields(obj))
        elseif obj isa NURBSSurface
            push!(records,_iges128_fields(obj))
        else
            throw(ArgumentError("export_iges_nurbs: expected NURBSCurve or NURBSSurface"))
        end
    end
    return _write_iges(path, records)
end

end # module
