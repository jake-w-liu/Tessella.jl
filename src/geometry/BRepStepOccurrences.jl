# Each map expands literal representation items. Source definitions are never
# mistaken for occurrences at their original placement.
const _STEP_OCCURRENCE_REPRESENTATIONS=("SHAPE_REPRESENTATION",
    "ADVANCED_BREP_SHAPE_REPRESENTATION","FACETED_BREP_SHAPE_REPRESENTATION")
const _STEP_OCCURRENCE_IDENTITY=(((1.,0.,0.),(0.,1.,0.),(0.,0.,1.)),(0.,0.,0.))

function _step_occurrence_entity(entities,reference,kind)
    entity=_deref_entity(entities,reference)
    entity!==nothing && entity.kind==kind || throw(ArgumentError(
        "import_step: mapped occurrence requires an actual $kind reference"))
    return entity
end

function _step_occurrence_vector(entity,kind)
    length(entity.args)==2 && entity.args[2] isa AbstractVector &&
        length(entity.args[2])==3 || throw(ArgumentError(
        "import_step: mapped occurrence $kind #$(entity.id) is truncated"))
    all(x->x isa Real && !(x isa Bool) && isfinite(x),entity.args[2]) ||
        throw(ArgumentError("import_step: mapped occurrence $kind has invalid coordinates"))
    return Tuple(Float64.(entity.args[2]))
end

function _step_occurrence_unit(vector)
    length=hypot(vector...)
    isfinite(length) && length>0 || throw(ArgumentError(
        "import_step: mapped occurrence has a zero or nonfinite direction"))
    return ntuple(i->vector[i]/length,3)
end

function _step_occurrence_frame(entities,reference)
    frame=_step_occurrence_entity(entities,reference,"AXIS2_PLACEMENT_3D")
    length(frame.args)==4 || throw(ArgumentError(
        "import_step: mapped occurrence placement #$(frame.id) is truncated"))
    origin=_step_occurrence_vector(
        _step_occurrence_entity(entities,frame.args[2],"CARTESIAN_POINT"),"origin")
    z=frame.args[3]===nothing ? (0.,0.,1.) : _step_occurrence_unit(
        _step_occurrence_vector(_step_occurrence_entity(entities,frame.args[3],"DIRECTION"),"axis"))
    default_x=abs(z[1])==1.0 && z[2]==0.0 && z[3]==0.0 ? (0.,1.,0.) : (1.,0.,0.)
    reference_x=frame.args[4]===nothing ? default_x : _step_occurrence_unit(
        _step_occurrence_vector(_step_occurrence_entity(entities,frame.args[4],"DIRECTION"),"reference"))
    projection=_brep_vdot(reference_x,z)
    x=_step_occurrence_unit(ntuple(i->reference_x[i]-projection*z[i],3))
    y=_step_occurrence_unit(_brep_vcross(z,x))
    return origin,(x,y,z)
end

# A Cartesian target has a direct occurrence frame, unlike the coordinate
# change between two AXIS2 placements. Only a rigid operator is admitted.
function _step_occurrence_cartesian_transform(entities,operator)
    length(operator.args)==8 && all(value->value isa AbstractString,operator.args[1:2]) &&
        (operator.args[3]===nothing || operator.args[3] isa AbstractString) ||
        throw(ArgumentError("import_step: Cartesian mapped target requires its eight actual fields"))
    scale=operator.args[7]
    (scale===nothing || scale isa Real && !(scale isa Bool) &&
        isfinite(scale) && scale==1) || throw(ArgumentError(
        "import_step: Cartesian mapped targets support unit or default scale only"))
    origin=_step_occurrence_vector(_step_occurrence_entity(
        entities,operator.args[6],"CARTESIAN_POINT"),"operator origin")
    function direction(reference,fallback,label)
        reference===nothing && return fallback
        return _step_occurrence_unit(_step_occurrence_vector(
            _step_occurrence_entity(entities,reference,"DIRECTION"),label))
    end
    first_axis=direction(operator.args[4],(1.,0.,0.),"operator axis1")
    second_axis=direction(operator.args[5],(0.,1.,0.),"operator axis2")
    third_axis=operator.args[8]===nothing ?
        _step_occurrence_unit(_brep_vcross(first_axis,second_axis)) :
        direction(operator.args[8],(0.,0.,1.),"operator axis3")
    projection=_brep_vdot(first_axis,third_axis)
    x=_step_occurrence_unit(ntuple(i->first_axis[i]-projection*third_axis[i],3))
    y=_step_occurrence_unit(_brep_vcross(third_axis,x))
    axes=(x,y,third_axis)
    return ntuple(i->ntuple(j->axes[j][i],3),3),origin
end

function _step_occurrence_transform(entities,origin,target)
    operator=_deref_entity(entities,target)
    if operator!==nothing && operator.kind=="CARTESIAN_TRANSFORMATION_OPERATOR_3D"
        # Keep the actual mapping-origin admission even though a Cartesian
        # target does not consume its frame in the pinned native transform.
        _step_occurrence_frame(entities,origin)
        return _step_occurrence_cartesian_transform(entities,operator)
    end
    a,A=_step_occurrence_frame(entities,origin)
    b,B=_step_occurrence_frame(entities,target)
    # The pinned native importer uses a coordinate change from the Target
    # frame to the Origin frame: A' * B, A' * (b-a). A translated Origin
    # followed by a rotated Target distinguishes this from B * A'.
    rotation=ntuple(i->ntuple(j->sum(A[i][k]*B[j][k] for k in 1:3),3),3)
    translation=ntuple(i->sum(A[i][j]*(b[j]-a[j]) for j in 1:3),3)
    return rotation,translation
end

function _step_occurrence_compose(parent,child)
    P,p=parent;C,c=child
    rotation=ntuple(i->ntuple(j->sum(P[i][k]*C[k][j] for k in 1:3),3),3)
    translation=ntuple(i->p[i]+sum(P[i][j]*c[j] for j in 1:3),3)
    all(isfinite,(rotation[1]...,rotation[2]...,rotation[3]...,translation...)) ||
        throw(ArgumentError("import_step: mapped occurrence transform is nonfinite"))
    columns=ntuple(j->ntuple(i->rotation[i][j],3),3)
    all(abs(_brep_vdot(columns[i],columns[j])-(i==j ? 1.0 : 0.0))<=1e-10
        for i in 1:3 for j in i:3) &&
        abs(_brep_vdot(columns[1],_brep_vcross(columns[2],columns[3]))-1.0)<=1e-10 ||
        throw(ArgumentError("import_step: mapped occurrence composition is not a proper rigid transform"))
    return rotation,translation
end

function _step_occurrence_import(entities,fill_solids)
    # These independent assembly relationships need their own producer. Never
    # fall through to the source-solid classifier and ignore their placements.
    unsupported=("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
        "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION","ITEM_DEFINED_TRANSFORMATION")
    for entity in values(entities)
        kinds=entity.kind=="COMPLEX" ? [part.kind for part in entity.args] : [entity.kind]
        any(kind->kind in unsupported,kinds) && throw(ArgumentError(
            "import_step: representation-relationship assembly transforms are unsupported"))
    end
    mapped=[entity for entity in values(entities) if entity.kind=="MAPPED_ITEM"]
    isempty(mapped) && return nothing
    definitions=[entity for entity in values(entities) if entity.kind=="SHAPE_DEFINITION_REPRESENTATION"]
    length(definitions)==1 && length(only(definitions).args)==2 || throw(ArgumentError(
        "import_step: mapped occurrences need one explicit shape-definition representation"))
    root=_deref_entity(entities,only(definitions).args[2])
    root!==nothing && root.kind in _STEP_OCCURRENCE_REPRESENTATIONS || throw(ArgumentError(
        "import_step: mapped occurrences have an unsupported root representation"))
    length(root.args)==3 || throw(ArgumentError("import_step: root representation is truncated"))
    context=root.args[3]
    context_entity=_deref_entity(entities,context)
    context_entity!==nothing || throw(ArgumentError(
        "import_step: mapped occurrence representation context is missing"))
    # A unit-bearing shared representation context declares the actual file
    # units; every source row and transform origin converts by the same
    # length factor. Contexts without unit declarations keep raw rows.
    unit_factor=haskey(_step_relationship_parts(context_entity),
        "GLOBAL_UNIT_ASSIGNED_CONTEXT") ?
        _step_relationship_context(entities,context)[3]["LENGTH_UNIT"].factor : 1.0
    records=NamedTuple[];source_surfaces=Dict{Int,Mesh}();active=Set{Int}();visited=Set{Int}()
    counts=[0,0]
    function expand(representation,transform,path,depth)
        depth<=128 || throw(ArgumentError("import_step: mapped occurrence nesting exceeds 128"))
        representation.id in active && throw(ArgumentError(
            "import_step: mapped representation cycle at #$(representation.id)"))
        representation.kind in _STEP_OCCURRENCE_REPRESENTATIONS &&
            length(representation.args)==3 && representation.args[2] isa AbstractVector &&
            representation.args[3]==context || throw(ArgumentError(
            "import_step: mapped representations need supported items and the same context"))
        push!(active,representation.id)
        emitted=length(records)
        for reference in representation.args[2]
            item=_deref_entity(entities,reference)
            item!==nothing || throw(ArgumentError("import_step: mapped representation item is missing"))
            if item.kind=="AXIS2_PLACEMENT_3D"
                _step_occurrence_frame(entities,reference)
            elseif item.kind=="MAPPED_ITEM"
                length(item.args)==3 || throw(ArgumentError("import_step: mapped item is truncated"))
                push!(visited,item.id)
                map=_step_occurrence_entity(entities,item.args[2],"REPRESENTATION_MAP")
                length(map.args)==2 || throw(ArgumentError("import_step: representation map is truncated"))
                source=_deref_entity(entities,map.args[2])
                source!==nothing || throw(ArgumentError("import_step: mapped source representation is missing"))
                local_transform=_step_occurrence_transform(entities,map.args[1],item.args[3])
                expand(source,_step_occurrence_compose(transform,local_transform),
                    path*"/MAPPED_ITEM #$(item.id)",depth+1)
            elseif item.kind in _STEP_BREP_SOLID_KINDS || item.kind in
                    ("SPHERE","RIGHT_CIRCULAR_CYLINDER","RIGHT_CIRCULAR_CONE","TOROIDAL_SURFACE")
                length(records)<4096 || throw(ArgumentError("import_step: more than 4096 solid occurrences"))
                surface=get!(source_surfaces,item.id) do
                    _step_root_surface(item,entities)
                end
                counts[1]=_checked_brep_add(counts[1],nnodes(surface),0,"import_step","occurrence nodes")
                counts[2]=_checked_brep_add(counts[2],size(surface.tris,2),0,"import_step","occurrence triangles")
                counts[1]<=_BREP_MAX_SURFACE_NODES && counts[2]<=_BREP_MAX_SURFACE_TRIS ||
                    throw(ArgumentError("import_step: mapped occurrences exceed surface storage limits"))
                rotation,translation=transform
                coordinates=Matrix{Float64}(undef,size(surface.coords))
                for column in axes(coordinates,2),i in 1:3
                    coordinates[i,column]=unit_factor*(translation[i]+
                        sum(rotation[i][j]*surface.coords[j,column] for j in 1:3))
                end
                all(isfinite,coordinates) || throw(ArgumentError("import_step: mapped occurrence coordinates are nonfinite"))
                push!(records,(name=path*"/$(item.kind) #$(item.id)",surface=Mesh(coordinates;tris=surface.tris)))
            else
                throw(ArgumentError("import_step: unsupported mapped representation item $(item.kind) #$(item.id)"))
            end
        end
        length(records)>emitted || throw(ArgumentError("import_step: mapped representation contains no solid occurrence"))
        delete!(active,representation.id)
        return nothing
    end
    expand(root,_STEP_OCCURRENCE_IDENTITY,"STEP representation #$(root.id)",0)
    visited==Set(item.id for item in mapped) || throw(ArgumentError(
        "import_step: mapped items outside the selected root are unsupported"))
    length(records)>1 && return _brep_assembly(records,fill_solids,"import_step")
    surface=only(records).surface
    validate(surface).ok || throw(ArgumentError("import_step: invalid mapped source solid"))
    _brep_assembly_shell_certificate(surface.coords,surface.tris,"import_step")
    return fill_solids ? _filled(surface,"import_step") : surface
end
