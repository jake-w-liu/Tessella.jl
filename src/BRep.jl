"""
    BRep

Native ISO-10303-21 STEP and IGES CAD import. Solids that classify as an
axis-aligned block, sphere, or right circular cylinder are converted to Tessella
surfaces and filled. Unrecognized topology is an explicit blocker, not a silent
empty mesh.
"""
module BRep

using ..Geometry: box_surface, cylinder_surface, sphere_surface
using ..Mesh3D: tetrahedralize
using ..MeshTypes: Mesh, validate, ntets, tet_volume, node

export import_step, import_iges, parse_step_entities

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
        k=j
        while k<=last && (isletter(data[k]) || isdigit(data[k]) || data[k]=='_')
            k=nextind(data,k)
        end
        kind=uppercase(data[j:prevind(data,k)])
        k=_step_skipws(data,k,last)
        k<=last && data[k]=='(' || throw(ArgumentError("import_step: entity #$id missing '('"))
        args, k=_step_parse_list(data,k)
        k=_step_skipws(data,k,last)
        k<=last && data[k]==';' || throw(ArgumentError("import_step: entity #$id missing ';'"))
        entities[id]=StepEntity(id,kind,args)
        i=nextind(data,k)
    end
    isempty(entities) && throw(ArgumentError("import_step: DATA section is empty"))
    return entities
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
    box=_as_box(points)
    if box!==nothing
        surface=box_surface(box...)
        fill || return surface
        return _filled(surface,"import_step")
    end
    kinds=sort!(unique(ent.kind for ent in values(entities)))
    throw(ArgumentError(
        "import_step: no supported solid (axis-aligned 8-corner block, SPHERE, " *
        "or RIGHT_CIRCULAR_CYLINDER); saw $(join(kinds, ", "))"))
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
        "import_iges: no supported solid (150 Block, 158 Sphere, 156 Cylinder); saw types $(types)"))
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

end # module
