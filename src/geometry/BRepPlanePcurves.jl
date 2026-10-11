# IGES 5.3 section 4.50 supplies this chart for a Form-1 Plane Surface:
# unit normal z, x = unit(refdir - dot(refdir,z)*z), y = unit(cross(z,x)),
# and position C + u*x + v*y. An implicit type-108 plane supplies no chart.
function _iges_uv_dependency(types,params,bad,de,kind,form,caller,context;
                             allow_transform::Bool=false)
    get(types,de,0)==kind || throw(ArgumentError(
        "$caller: UV chart dependency #$de is not entity type $kind"))
    get(context.forms,de,nothing)==form || throw(ArgumentError(
        "$caller: UV chart dependency #$de needs explicit form $form"))
    allow_transform || get(context.transforms,de,0)==0 || throw(ArgumentError(
        "$caller: UV chart dependency #$de has an unsupported directory transform"))
    rec=_iges_rec(params,bad,de,caller)
    _iges_param_int(rec,1,caller,"UV chart entity type")==kind || throw(ArgumentError(
        "$caller: UV chart dependency #$de parameter type disagrees with its directory entry"))
    return rec
end

function _iges_plane_uv_rotation_certificate(rotation,caller,description)
    columns=ntuple(axis->ntuple(row->rotation[row][axis],3),3)
    for i in 1:3,j in i:3
        want=i==j ? 1.0 : 0.0
        abs(_brep_vdot(columns[i],columns[j])-want)<=1e-6 ||
            throw(ArgumentError("$caller: UV chart $description is not rigid"))
    end
    determinant=_brep_vdot(columns[1],_brep_vcross(columns[2],columns[3]))
    determinant>0 && abs(determinant-1)<=1e-6 || throw(ArgumentError(
        "$caller: UV chart $description is not a proper rotation"))
    return nothing
end

function _iges_plane_uv_transform(types,params,bad,surface,caller,context)
    matrix=get(context.transforms,surface,0)
    matrix==0 && return nothing
    # IGES 5.3 section4.21: a matrix's DE field7 points to the next,
    # later transform. Compose parentR*childR and parentT+parentR*childT.
    visited=Set{Int}();owner=surface;combined=nothing
    while matrix!=0
        matrix>0 && isodd(matrix) || throw(ArgumentError(
            "$caller: UV chart defining matrix pointer $matrix is not a positive directory identity"))
        matrix in visited && throw(ArgumentError(
            "$caller: UV chart defining matrices contain a cycle at #$matrix"))
        push!(visited,matrix)
        _iges_uv_dependency(types,params,bad,matrix,124,0,caller,context;
                            allow_transform=true)
        transform=_iges_de_transform(types,params,bad,context.transforms,owner,caller)
        rotation,_=transform
        _iges_plane_uv_rotation_certificate(rotation,caller,"defining matrix #$matrix")
        if combined===nothing
            combined=transform
        else
            child_rotation,child_translation=combined
            composed_rotation=ntuple(row->ntuple(column->sum(
                rotation[row][axis]*child_rotation[axis][column] for axis in 1:3),3),3)
            composed_translation=_iges_xf_apply(transform,child_translation,true)
            all(isfinite,(composed_rotation[1]...,composed_rotation[2]...,
                          composed_rotation[3]...,composed_translation...)) ||
                throw(ArgumentError("$caller: UV chart defining matrix composition is non-finite"))
            _iges_plane_uv_rotation_certificate(composed_rotation,caller,"matrix composition at #$matrix")
            combined=(composed_rotation,composed_translation)
        end
        owner=matrix;matrix=get(context.transforms,matrix,0)
    end
    return combined
end

function _iges_uv_unit(vector,caller)
    magnitude=hypot(vector...)
    isfinite(magnitude) && magnitude>0 || throw(ArgumentError(
        "$caller: UV chart direction has no finite positive length"))
    return ntuple(axis->vector[axis]/magnitude,3)
end

function _iges_plane_uv_frame(types,params,bad,surface,caller,context)
    rec=_iges_uv_dependency(types,params,bad,surface,190,1,caller,context;
                            allow_transform=true)
    length(rec)>=4 || throw(ArgumentError("$caller: parameterized plane #$surface has no reference direction"))
    point=_iges_param_int(rec,2,caller,"UV chart origin")
    normal=_iges_param_int(rec,3,caller,"UV chart normal")
    reference=_iges_param_int(rec,4,caller,"UV chart reference direction")
    p=_iges_uv_dependency(types,params,bad,point,116,0,caller,context)
    n=_iges_uv_dependency(types,params,bad,normal,123,0,caller,context)
    d=_iges_uv_dependency(types,params,bad,reference,123,0,caller,context)
    all(record->length(record)>=4,(p,n,d)) || throw(ArgumentError(
        "$caller: UV chart point or direction record is truncated"))
    origin=(p[2],p[3],p[4])
    z=_iges_uv_unit((n[2],n[3],n[4]),caller)
    ref=_iges_uv_unit((d[2],d[3],d[4]),caller)
    projection=_brep_vdot(ref,z)
    abs(projection)<=1e-9 || throw(ArgumentError(
        "$caller: UV chart reference direction is not approximately perpendicular to its normal"))
    x=_iges_uv_unit(ntuple(axis->ref[axis]-projection*z[axis],3),caller)
    y=_iges_uv_unit(_brep_vcross(z,x),caller)
    transform=_iges_plane_uv_transform(types,params,bad,surface,caller,context)
    return _iges_xf_apply(transform,origin,true),
           _iges_xf_apply(transform,x,false),_iges_xf_apply(transform,y,false)
end

function _iges_plane_uv_segments(types,params,bad,de,caller,context,active)
    de in active && throw(ArgumentError("$caller: UV composite curve references a cycle at #$de"))
    kind=get(types,de,0)
    kind in (110,102) || throw(ArgumentError(
        "$caller: UV-only plane boundaries need literal Line/CompositeCurve records; entity #$de has type $kind"))
    rec=_iges_uv_dependency(types,params,bad,de,kind,0,caller,context)
    if kind==110
        length(rec)>=7 || throw(ArgumentError("$caller: UV line #$de is truncated"))
        rec[4]==0 && rec[7]==0 || throw(ArgumentError("$caller: UV line #$de has nonzero Z coordinates"))
        a=(rec[2],rec[3]);b=(rec[5],rec[6])
        a!=b || throw(ArgumentError("$caller: UV line #$de has zero length"))
        return Tuple{NTuple{2,Float64},NTuple{2,Float64}}[(a,b)]
    end
    count=_iges_param_int(rec,2,caller,"UV composite member count")
    count>=1 && length(rec)>=2+count || throw(ArgumentError(
        "$caller: UV composite #$de is empty or truncated"))
    push!(active,de)
    segments=Tuple{NTuple{2,Float64},NTuple{2,Float64}}[]
    for index in 1:count
        member=_iges_param_int(rec,2+index,caller,"UV composite member")
        append!(segments,_iges_plane_uv_segments(types,params,bad,member,caller,context,active))
    end
    delete!(active,de)
    return segments
end

function _iges_plane_uv_wire(types,params,bad,boundary,parameter,rec,caller,
                             context,expected_surface)
    context!==nothing || throw(ArgumentError(
        "$caller: UV-only boundary #$boundary has no declared plane chart metadata"))
    _iges_uv_dependency(types,params,bad,boundary,142,0,caller,context)
    surface=_iges_param_int(rec,3,caller,"UV boundary surface")
    expected_surface!==nothing && surface==expected_surface || throw(ArgumentError(
        "$caller: UV boundary #$boundary does not bind its declared Face surface"))
    parameter!=0 || throw(ArgumentError("$caller: curve-on-surface #$boundary has no boundary curve"))
    preference=_iges_param_int(rec,6,caller,"UV boundary preference")
    preference in (0,1,2,3) || throw(ArgumentError("$caller: UV boundary #$boundary has invalid preference $preference"))
    origin,x,y=_iges_plane_uv_frame(types,params,bad,surface,caller,context)
    segments=_iges_plane_uv_segments(types,params,bad,parameter,caller,context,Set{Int}())
    length(segments)>=3 || throw(ArgumentError("$caller: UV boundary #$boundary has fewer than three edges"))
    all(index->segments[index][2]==segments[index+1][1],1:length(segments)-1) &&
        last(segments)[2]==first(segments)[1] || throw(ArgumentError(
            "$caller: UV boundary #$boundary is not a literal closed ordered wire"))
    ring=NTuple{3,Float64}[]
    for (uv,_) in segments
        position=ntuple(axis->origin[axis]+uv[1]*x[axis]+uv[2]*y[axis],3)
        all(isfinite,position) || throw(ArgumentError("$caller: mapped UV boundary #$boundary is non-finite"))
        push!(ring,position)
    end
    return ring
end
