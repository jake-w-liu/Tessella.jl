# Mesh import retains the source solid partition. These are discrete closed
# boundary entities, not reconstructed analytic CAD Face/Curve parameter maps.
function _brep_assembly(records,fill_solids,caller)
    length(records)>1 || throw(ArgumentError("$caller: assembly requires multiple source solids"))
    length(records)<=typemax(Int32) || throw(ArgumentError("$caller: too many source solids"))
    bounds=[bounding_box(record.surface) for record in records]
    for i in 1:length(records)-1,j in i+1:length(records)
        alo,ahi=bounds[i];blo,bhi=bounds[j]
        any(axis->ahi[axis]<blo[axis] || bhi[axis]<alo[axis],1:3) ||
            throw(ArgumentError("$caller: cannot certify independent solids $(records[i].name) and " *
                "$(records[j].name); touching or overlapping bounding boxes are unsupported"))
    end
    # Certify every source boundary before any fill, and fill each Region
    # independently. No coordinate matching can select a Region or weld UIDs.
    parts=Mesh[];shells=Matrix{Int32}[]
    for record in records
        validate(record.surface).ok || throw(ArgumentError("$caller: invalid source solid $(record.name)"))
        _brep_assembly_shell_certificate(record.surface.coords,record.surface.tris,caller)
        part=fill_solids ? _filled(record.surface,caller) : record.surface
        shell=fill_solids ? _brep_tet_shell(part,caller) : copy(part.tris)
        _brep_assembly_shell_certificate(part.coords,shell,caller)
        push!(parts,part);push!(shells,shell)
    end
    node_count=0;cell_count=0
    for (part,shell) in zip(parts,shells)
        node_count=_checked_brep_add(node_count,nnodes(part),0,caller,"assembly node count")
        cell_count=_checked_brep_add(cell_count,size(shell,2),ntets(part),caller,"assembly cell count")
    end
    node_count<=typemax(Int32) && cell_count<=typemax(Int32) ||
        throw(ArgumentError("$caller: assembly exceeds Int32 topology limits"))
    coords=Matrix{Float64}(undef,3,node_count)
    entities=Dict{Tuple{Int,Int},MixedEntity}()
    names=Dict{Tuple{Int,Int},String}()
    owners=Tuple{Int,Int32}[];blocks=ElementBlock[]
    block_owners=Vector{Int32}[];labels=Vector{UInt64}[]
    node_offset=0;cell_offset=UInt64(0)
    for (region,(record,part,shell)) in enumerate(zip(records,parts,shells))
        tag=Int32(region);lo,hi=bounds[region];bbox=(lo...,hi...)
        entities[(2,region)]=MixedEntity(2,tag,bbox;physical_tags=[tag])
        entities[(3,region)]=MixedEntity(3,tag,bbox;physical_tags=[tag],boundaries=[tag])
        names[(2,region)]=record.name*" boundary"
        names[(3,region)]=record.name
        count=nnodes(part);coords[:,node_offset+1:node_offset+count]=part.coords
        boundary=falses(count)
        for row in shell;boundary[row]=true;end
        append!(owners,[(boundary[row] ? 2 : 3,tag) for row in 1:count])
        for (kind,cells) in ((2,shell),(4,part.tets))
            isempty(cells) && continue
            push!(blocks,ElementBlock(kind,cells .+ Int32(node_offset),fill(tag,size(cells,2))))
            push!(block_owners,fill(tag,size(cells,2)))
            push!(labels,UInt64[cell_offset+UInt64(i) for i in 1:size(cells,2)])
            cell_offset+=UInt64(size(cells,2))
        end
        node_offset+=count
    end
    data=MixedEntityData(entities;node_entities=owners,
        node_parametric=fill(nothing,node_count),external_node_tags=UInt64.(1:node_count),
        block_entities=block_owners,external_element_tags=labels)
    result=MixedMesh(coords,blocks;physical_names=names,entity_data=data,
                     elementary_entities=block_owners)
    validate(result).ok || throw(ErrorException("$caller: invalid classified assembly product"))
    return result
end

function _brep_assembly_shell_certificate(coords,shell,caller)
    isempty(shell) && throw(ArgumentError("$caller: source solid has no closed shell"))
    incidences=Dict{Tuple{Int32,Int32},Tuple{Int,Int}}()
    for face in eachcol(shell)
        a,b,c=face
        for (u,v) in ((a,b),(b,c),(c,a))
            key=minmax(u,v);count,orientation=get(incidences,key,(0,0))
            incidences[key]=(count+1,orientation+(u<v ? 1 : -1))
        end
    end
    all(value->value==(2,0),values(incidences)) || throw(ArgumentError(
        "$caller: source solid boundary is not a consistently oriented closed shell"))
    return nothing
end

# Literal positive-tet face traversal supplies outward boundary cells. A Dict
# counts incidence only; the first-encounter vector controls emitted order.
function _brep_tet_shell(mesh,caller)
    counts=Dict{NTuple{3,Int32},Int}();faces=NTuple{3,Int32}[];keys=NTuple{3,Int32}[]
    for cell in eachcol(mesh.tets)
        a,b,c,d=cell
        for face in ((b,c,d),(a,d,c),(a,b,d),(a,c,b))
            key=Tuple(sort(collect(face)))
            if !haskey(counts,key);push!(keys,key);push!(faces,face);end
            counts[key]=get(counts,key,0)+1
            counts[key]<=2 || throw(ArgumentError("$caller: nonmanifold filled solid boundary"))
        end
    end
    selected=[faces[i] for i in eachindex(faces) if counts[keys[i]]==1]
    isempty(selected) && throw(ArgumentError("$caller: filled solid has no exterior shell"))
    shell=Matrix{Int32}(undef,3,length(selected))
    for (column,face) in enumerate(selected);shell[:,column]=collect(face);end
    return shell
end

function _iges_assembly_transform(surface,types,params,bad,transforms,de,caller)
    xf=_iges_de_transform(types,params,bad,transforms,de,caller)
    xf===nothing && return surface
    coordinates=Matrix{Float64}(undef,size(surface.coords))
    for row in axes(coordinates,2)
        coordinates[:,row]=collect(_iges_xf_apply(xf,node(surface,row),true))
    end
    R,_=xf
    determinant=R[1][1]*(R[2][2]*R[3][3]-R[2][3]*R[3][2])-
        R[1][2]*(R[2][1]*R[3][3]-R[2][3]*R[3][1])+R[1][3]*(R[2][1]*R[3][2]-R[2][2]*R[3][1])
    return Mesh(coordinates;tris=determinant<0 ? surface.tris[[1,3,2],:] : surface.tris)
end

function _step_solid_roots(entities)
    breps=[ent for ent in values(entities) if ent.kind in _STEP_BREP_SOLID_KINDS]
    descendants=Set{Int}();queue=Any[ent.args for ent in breps]
    while !isempty(queue)
        value=pop!(queue)
        if value isa Tuple && length(value)==2 && value[1]===:ref
            id=value[2];id in descendants && continue
            push!(descendants,id);ent=get(entities,id,nothing)
            ent===nothing || push!(queue,ent.args)
        elseif value isa AbstractVector
            append!(queue,value)
        elseif value isa StepEntity
            push!(queue,value.args)
        end
    end
    primitives=[ent for ent in values(entities) if !(ent.id in descendants) &&
        ent.kind in ("SPHERE","RIGHT_CIRCULAR_CYLINDER","RIGHT_CIRCULAR_CONE","TOROIDAL_SURFACE")]
    return sort!(vcat(breps,primitives);by=ent->ent.id)
end

function _step_root_surface(ent,entities)
    ent.kind in _STEP_BREP_SOLID_KINDS && return _step_brep_solid(ent,entities,false)
    # Retain placement references, but expose only this root to the existing
    # single-primitive classifier. Selecting a root never depends on Dict order.
    selected=Dict(id=>other for (id,other) in entities if id==ent.id ||
        !(other.kind in ("SPHERE","RIGHT_CIRCULAR_CYLINDER","RIGHT_CIRCULAR_CONE","TOROIDAL_SURFACE")))
    if ent.kind=="SPHERE"
        params=_as_sphere(selected)
        params===nothing || return sphere_surface(params...)
    elseif ent.kind=="RIGHT_CIRCULAR_CYLINDER"
        params=_as_cylinder(selected)
        params===nothing || return cylinder_surface(params...)
    elseif ent.kind=="RIGHT_CIRCULAR_CONE"
        params=_as_cone(selected)
        params===nothing || return cone_surface(params...)
    elseif ent.kind=="TOROIDAL_SURFACE"
        params=_step_torus_params(entities,ent,"import_step")
        return _torus_surface(params...,"import_step")
    end
    throw(ArgumentError("import_step: source solid #$(ent.id) $(ent.kind) is truncated"))
end

function _step_assembly_records(roots,entities)
    for ent in values(entities)
        endswith(ent.kind,"SOLID") && !(ent.kind in _STEP_BREP_SOLID_KINDS) &&
            throw(ArgumentError("import_step: unsupported solid root #$(ent.id) $(ent.kind)"))
    end
    factors=_step_member_unit_factors(entities)
    return [(name="STEP #$(ent.id) $(ent.kind)"*
                 (!isempty(ent.args) && ent.args[1] isa String && !isempty(ent.args[1]) ? " "*ent.args[1] : ""),
             surface=_step_unit_scaled(_step_root_surface(ent,entities),
                 get(factors,ent.id,1.0))) for ent in roots]
end

function _iges_primitive_surface(rec)
    type=_as_int(rec[1],"import_iges","entity type")
    if type==150 && length(rec)>=7
        L,W,H,x,y,z=rec[2:7]
        (L>0 && W>0 && H>0) || throw(ArgumentError("import_iges: Block extents must be positive"))
        return box_surface(x,x+L,y,y+W,z,z+H)
    elseif type==158 && length(rec)>=5
        r,x,y,z=rec[2:5]
        r>0 || throw(ArgumentError("import_iges: Sphere radius must be positive"))
        return sphere_surface((x,y,z),r)
    elseif type==156 && length(rec)>=10
        h,r1,r2,x,y,z,zi,zj,zk=rec[2:10]
        (h>0 && (r1>0 || r2>0) && r1>=0 && r2>=0) || throw(ArgumentError(
            "import_iges: Cone height must be positive and at least one radius must be positive"))
        return cone_surface((x,y,z),(zi,zj,zk),r1,r2,h)
    elseif type==156 && length(rec)>=8
        r,h,x,y,z,zi,zj=rec[2:8];zk=length(rec)>=9 ? rec[9] : 1.0
        (r>0 && h>0) || throw(ArgumentError("import_iges: Cylinder radius/height must be positive"))
        return cylinder_surface((x,y,z),(zi,zj,zk),r,h)
    end
    type in (150,156,158) && throw(ArgumentError("import_iges: primitive type $type is truncated"))
    return nothing
end

function _iges_linked_primitive_entries(types,params,bad)
    ids=Set(de for (de,type) in types if type in (150,156,158))
    for (de,rec) in params
        !isempty(rec) && rec[1] in (150,156,158) && push!(ids,de)
    end
    return [(id=de,record=_iges_rec(params,bad,de,"import_iges")) for de in sort!(collect(ids))]
end

function _iges_assembly_directory_certificate(source,types)
    # Admit the implemented solid and support layouts, plus property records.
    # An unknown geometric entity must not silently disappear from an assembly.
    supported=(100,102,108,110,116,120,123,124,141,142,143,144,150,156,158,
               160,186,190,198,406,502,504,508,510,514)
    all(type->type in supported,values(types)) || throw(ArgumentError(
        "import_iges: assembly contains an unsupported directory-entry entity type"))
    seen=Set{Int}();root_records=Dict{Int,Int}()
    for raw in split(source,r"\r?\n")
        line=length(raw)>=80 ? raw[1:80] : rpad(raw,80)
        if line[73]=='D'
            de=_iges_int_field(line[74:80])
            de===nothing && throw(ArgumentError("import_iges: invalid directory-entry sequence"))
            de>0 || throw(ArgumentError("import_iges: directory-entry sequences must be positive"))
            de in seen && throw(ArgumentError("import_iges: duplicate directory-entry sequence #$de"))
            push!(seen,de)
        elseif line[73]=='P'
            de=_iges_int_field(line[65:72])
            get(types,de,0) in (150,156,158,160,186,198) || continue
            occursin(';',line[1:64]) || continue
            root_records[de]=get(root_records,de,0)+1
            root_records[de]==1 || throw(ArgumentError(
                "import_iges: multiple parameter records claim solid directory entry #$de"))
        end
    end
    return nothing
end

function _iges_assembly_records(primitive_entries,types,params,bad,transforms;frame_context=nothing)
    caller="import_iges";records=NamedTuple[]
    msbo=sort!([de for (de,type) in types if type==186])
    trimmed=sort!([de for (de,type) in types if type==144])
    bounded=sort!([de for (de,type) in types if type==143])
    referenced_tori=Set{Int}()
    if !isempty(msbo)
        isempty(trimmed) && isempty(bounded) || throw(ArgumentError(
            "$caller: mixed MSBO and ungrouped trimmed Face roots cannot be assigned independent solid identities"))
        for de in msbo
            faces,_,_=_iges_brep_faces(types,params,bad,[de],Int[],Int[],caller;frame_context)
            union!(referenced_tori,(face[3] for face in faces))
            surface=_iges_brep_surface(types,params,bad,transforms,[de],Int[],Int[],caller;frame_context)
            surface=_iges_assembly_transform(surface,types,params,bad,transforms,de,caller)
            push!(records,(id=de,name="IGES #$de type 186",surface=surface))
        end
    elseif !isempty(trimmed) || !isempty(bounded)
        # The legacy face collection declares one closed shell. Without an
        # MSBO root there is no source identity for splitting it into solids.
        surface=_iges_brep_surface(types,params,bad,transforms,Int[],trimmed,bounded,caller;frame_context)
        faces,_,_=_iges_brep_faces(types,params,bad,Int[],trimmed,bounded,caller;frame_context)
        union!(referenced_tori,(face[3] for face in faces))
        id=minimum(vcat(trimmed,bounded))
        push!(records,(id=id,name="IGES #$id closed shell",surface=surface))
    end
    for de in sort!([de for (de,type) in types if type in (160,198) && !(de in referenced_tori)])
        push!(records,(id=de,name="IGES #$de type $(types[de])",
            surface=_iges_assembly_transform(_iges_torus_surface(types,params,bad,de,caller),
                types,params,bad,transforms,de,caller)))
    end
    for entry in primitive_entries
        rec=entry.record;type=_as_int(rec[1],caller,"entity type")
        if !(type in (150,156,158))
            isempty(types) && !(type in (124,406)) && throw(ArgumentError(
                "$caller: assembly contains unsupported parameter entity type $type"))
            continue
        end
        entry.id!==nothing && entry.id>0 && isodd(entry.id) || throw(ArgumentError(
            "$caller: multiple recognized solids require positive odd IGES directory-entry pointers"))
        isempty(types) || get(types,entry.id,0)==type || throw(ArgumentError(
            "$caller: primitive parameter type and directory entry #$(entry.id) disagree"))
        surface=_iges_primitive_surface(rec)
        surface=_iges_assembly_transform(surface,types,params,bad,transforms,entry.id,caller)
        push!(records,(id=entry.id,name="IGES #$(entry.id) type $type",surface=surface))
    end
    allunique(record.id for record in records) || throw(ArgumentError(
        "$caller: multiple recognized solids require unique IGES directory-entry pointers"))
    all(record->record.id>0 && isodd(record.id),records) || throw(ArgumentError(
        "$caller: solid roots require positive odd IGES directory-entry identities"))
    for record in records
        haskey(types,record.id) || continue
        parameters=_iges_rec(params,bad,record.id,caller)
        !isempty(parameters) && _as_int(parameters[1],caller,"solid parameter type")==types[record.id] ||
            throw(ArgumentError("$caller: solid parameter type and directory entry #$(record.id) disagree"))
    end
    return sort!(records;by=record->record.id)
end
