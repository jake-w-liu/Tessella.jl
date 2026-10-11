# Product-usage placements are a separate producer from representation maps.
# Every selected placement binds literal product and representation references;
# source definitions are never emitted as unplaced independent solids.
const _STEP_RELATIONSHIP_TRIGGERS=("NEXT_ASSEMBLY_USAGE_OCCURRENCE",
    "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION",
    "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION","ITEM_DEFINED_TRANSFORMATION")

function _step_relationship_parts(entity)
    parts=entity.kind=="COMPLEX" ? entity.args : [entity]
    all(part->part isa StepEntity,parts) &&
        allunique(part.kind for part in parts) || throw(ArgumentError(
        "import_step: ambiguous complex assembly entity #$(entity.id)"))
    return Dict(part.kind=>part for part in parts)
end

function _step_relationship_entity(entities,reference,kind)
    reference isa Tuple && length(reference)==2 && reference[1]===:ref ||
        throw(ArgumentError("import_step: assembly $kind requires a literal reference"))
    entity=_deref_entity(entities,reference)
    entity!==nothing && entity.kind==kind || throw(ArgumentError(
        "import_step: assembly requires an actual $kind reference"))
    return entity
end

function _step_relationship_context(entities,reference)
    context=_deref_entity(entities,reference)
    context!==nothing || throw(ArgumentError("import_step: assembly context is missing"))
    parts=_step_relationship_parts(context)
    required=Set(("REPRESENTATION_CONTEXT","GEOMETRIC_REPRESENTATION_CONTEXT",
        "GLOBAL_UNIT_ASSIGNED_CONTEXT"))
    required ⊆ Set(keys(parts)) && Set(keys(parts)) ⊆ union(required,
        Set(("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",))) || throw(ArgumentError(
        "import_step: assembly needs an explicit geometric unit context"))
    labels=parts["REPRESENTATION_CONTEXT"].args
    length(labels)==2 && all(value->value isa AbstractString,labels) ||
        throw(ArgumentError("import_step: assembly context labels are truncated"))
    dimension=parts["GEOMETRIC_REPRESENTATION_CONTEXT"].args
    length(dimension)==1 && dimension[1] isa Real && !(dimension[1] isa Bool) &&
        dimension[1]==3 || throw(ArgumentError("import_step: assembly context must be three-dimensional"))
    assigned=parts["GLOBAL_UNIT_ASSIGNED_CONTEXT"].args
    length(assigned)==1 && assigned[1] isa AbstractVector && length(assigned[1])==3 &&
        allunique(assigned[1]) || throw(ArgumentError("import_step: assembly unit context is truncated"))
    kinds=Set{String}();units=Set{Int}();length_unit=0
    for reference in assigned[1]
        unit=_deref_entity(entities,reference)
        unit!==nothing || throw(ArgumentError("import_step: assembly unit reference is missing"))
        unit_parts=_step_relationship_parts(unit)
        category=intersect(Set(keys(unit_parts)),Set(("LENGTH_UNIT","PLANE_ANGLE_UNIT","SOLID_ANGLE_UNIT")))
        length(category)==1 && Set(keys(unit_parts))==union(category,Set(("NAMED_UNIT","SI_UNIT"))) ||
            throw(ArgumentError("import_step: assembly supports explicit SI unit contexts only"))
        kind=only(category)
        isempty(unit_parts[kind].args) && unit_parts["NAMED_UNIT"].args==Any[nothing] ||
            throw(ArgumentError("import_step: assembly unit declaration is malformed"))
        si=unit_parts["SI_UNIT"].args
        expected=kind=="LENGTH_UNIT" ? ".METRE." : kind=="PLANE_ANGLE_UNIT" ? ".RADIAN." : ".STERADIAN."
        prefix=kind=="LENGTH_UNIT" ? ".MILLI." : nothing
        # The existing solid producer retains literal source coordinates.
        # Equal context references alone do not implement source-unit scaling.
        length(si)==2 && si[2]==expected && si[1]==prefix || throw(ArgumentError(
            "import_step: product usages currently require millimetre/radian/steradian units; conversion is unsupported"))
        kind in kinds && throw(ArgumentError("import_step: duplicate assembly unit category"))
        push!(kinds,kind);push!(units,unit.id)
        kind=="LENGTH_UNIT" && (length_unit=unit.id)
    end
    if haskey(parts,"GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT")
        uncertainty=parts["GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT"].args
        length(uncertainty)==1 && uncertainty[1] isa AbstractVector ||
            throw(ArgumentError("import_step: assembly uncertainty context is truncated"))
        for reference in uncertainty[1]
            measure=_step_relationship_entity(entities,reference,"UNCERTAINTY_MEASURE_WITH_UNIT")
            length(measure.args)==4 && all(value->value isa AbstractString,measure.args[3:4]) ||
                throw(ArgumentError("import_step: assembly uncertainty declaration is truncated"))
            value=measure.args[1];unit=_deref_entity(entities,measure.args[2])
            value isa StepEntity && value.kind=="LENGTH_MEASURE" && length(value.args)==1 &&
                value.args[1] isa Real && !(value.args[1] isa Bool) && isfinite(value.args[1]) &&
                value.args[1]>=0 && unit!==nothing && unit.id==length_unit ||
                throw(ArgumentError("import_step: invalid assembly uncertainty measure"))
        end
    end
    return context.id,units
end

function _step_relationship_product(entities,definition)
    length(definition.args)==4 && all(value->value isa AbstractString,definition.args[1:2]) ||
        throw(ArgumentError("import_step: assembly product definition is truncated"))
    formation=_step_relationship_entity(entities,definition.args[3],"PRODUCT_DEFINITION_FORMATION")
    length(formation.args)==3 && all(value->value isa AbstractString,formation.args[1:2]) ||
        throw(ArgumentError("import_step: assembly product formation is truncated"))
    product=_step_relationship_entity(entities,formation.args[3],"PRODUCT")
    length(product.args)==4 && all(value->value isa AbstractString,product.args[1:3]) &&
        product.args[4] isa AbstractVector && !isempty(product.args[4]) ||
        throw(ArgumentError("import_step: assembly product is truncated"))
    definition_context=_step_relationship_entity(entities,definition.args[4],"PRODUCT_DEFINITION_CONTEXT")
    length(definition_context.args)==3 && all(value->value isa AbstractString,
        definition_context.args[[1,3]]) || throw(ArgumentError("import_step: assembly product context is truncated"))
    application=_step_relationship_entity(entities,definition_context.args[2],"APPLICATION_CONTEXT")
    length(application.args)==1 && application.args[1] isa AbstractString ||
        throw(ArgumentError("import_step: assembly application context is truncated"))
    for reference in product.args[4]
        context=_step_relationship_entity(entities,reference,"PRODUCT_CONTEXT")
        length(context.args)==3 && all(value->value isa AbstractString,context.args[[1,3]]) &&
            context.args[2]==definition_context.args[2] ||
            throw(ArgumentError("import_step: incompatible assembly product contexts"))
    end
    return product,formation
end

function _step_relationship_import(entities,fill_solids)
    kinds(entity)=entity.kind=="COMPLEX" ? [part.kind for part in entity.args] : [entity.kind]
    any(entity->any(kind->kind in _STEP_RELATIONSHIP_TRIGGERS,kinds(entity)),values(entities)) || return nothing
    for entity in values(entities)
        any(kind->kind in ("MAPPED_ITEM","REPRESENTATION_MAP","ASSEMBLY_COMPONENT_USAGE",
            "PRODUCT_DEFINITION_USAGE","PRODUCT_DEFINITION_RELATIONSHIP","SPECIFIED_HIGHER_USAGE_OCCURRENCE"),
            kinds(entity)) && throw(ArgumentError("import_step: unsupported mixed assembly relationship workflow"))
        endswith(entity.kind,"SOLID") && !(entity.kind in _STEP_BREP_SOLID_KINDS) &&
            throw(ArgumentError("import_step: unsupported assembly solid #$(entity.id) $(entity.kind)"))
    end
    selected(kind)=sort!([entity for entity in values(entities) if kind in kinds(entity)];by=entity->entity.id)
    usages=selected("NEXT_ASSEMBLY_USAGE_OCCURRENCE")
    !isempty(usages) && length(usages)<=4096 ||
        throw(ArgumentError("import_step: assembly needs one to 4096 actual direct usages"))
    parents=Set{Int}();children=Set{Int}()
    for usage in usages
        usage.kind=="NEXT_ASSEMBLY_USAGE_OCCURRENCE" && length(usage.args)==6 &&
            all(value->value isa AbstractString,usage.args[1:3]) &&
            (usage.args[6]===nothing || usage.args[6] isa AbstractString) ||
            throw(ArgumentError("import_step: assembly usage is truncated"))
        parent=_step_relationship_entity(entities,usage.args[4],"PRODUCT_DEFINITION")
        child=_step_relationship_entity(entities,usage.args[5],"PRODUCT_DEFINITION")
        push!(parents,parent.id);push!(children,child.id)
    end
    length(parents)==1 && isempty(intersect(parents,children)) ||
        throw(ArgumentError("import_step: only acyclic direct usages from one assembly root are supported"))
    definitions=union(parents,children)
    Set(entity.id for entity in selected("PRODUCT_DEFINITION"))==definitions ||
        throw(ArgumentError("import_step: assembly has unreachable product definitions"))
    products=Dict{Int,StepEntity}();formations=Set{Int}()
    representations=Dict{Int,NamedTuple}();shape_definitions=Set{Int}()
    for definition in selected("PRODUCT_DEFINITION")
        product,formation=_step_relationship_product(entities,definition)
        products[definition.id]=product;push!(formations,formation.id)
    end
    Set(entity.id for entity in selected("PRODUCT"))==Set(product.id for product in values(products)) &&
        Set(entity.id for entity in selected("PRODUCT_DEFINITION_FORMATION"))==formations ||
        throw(ArgumentError("import_step: assembly has unreachable products or formations"))
    for binding in selected("SHAPE_DEFINITION_REPRESENTATION")
        binding.kind=="SHAPE_DEFINITION_REPRESENTATION" && length(binding.args)==2 ||
            throw(ArgumentError("import_step: assembly shape binding is truncated"))
        shape=_step_relationship_entity(entities,binding.args[1],"PRODUCT_DEFINITION_SHAPE")
        length(shape.args)==3 && all(value->value isa AbstractString,shape.args[1:2]) ||
            throw(ArgumentError("import_step: assembly product shape is truncated"))
        definition=_step_relationship_entity(entities,shape.args[3],"PRODUCT_DEFINITION")
        definition.id in definitions && !haskey(representations,definition.id) ||
            throw(ArgumentError("import_step: ambiguous assembly product shape binding"))
        representation=_deref_entity(entities,binding.args[2])
        representation!==nothing && representation.kind in _STEP_OCCURRENCE_REPRESENTATIONS &&
            length(representation.args)==3 && representation.args[1] isa AbstractString &&
            representation.args[2] isa AbstractVector && !isempty(representation.args[2]) &&
            allunique(representation.args[2]) ||
            throw(ArgumentError("import_step: assembly representation is truncated or unsupported"))
        context,units=_step_relationship_context(entities,representation.args[3])
        items=StepEntity[]
        for reference in representation.args[2]
            item=_deref_entity(entities,reference)
            item!==nothing || throw(ArgumentError("import_step: assembly representation item is missing"))
            item.kind=="AXIS2_PLACEMENT_3D" && _step_occurrence_frame(entities,reference)
            push!(items,item)
        end
        representations[definition.id]=(;representation,context,units,items)
        push!(shape_definitions,shape.id)
    end
    Set(keys(representations))==definitions &&
        allunique(info.representation.id for info in values(representations)) ||
        throw(ArgumentError("import_step: assembly product representations are incomplete or shared ambiguously"))
    Set(entity.id for entity in values(entities) if entity.kind in _STEP_OCCURRENCE_REPRESENTATIONS)==
        Set(info.representation.id for info in values(representations)) ||
        throw(ArgumentError("import_step: assembly has unselected shape representations"))
    root=only(parents);root_info=representations[root]
    all(item->item.kind=="AXIS2_PLACEMENT_3D",root_info.items) ||
        throw(ArgumentError("import_step: assembly root must contain placement frames only"))
    solid_roots=Dict(entity.id=>entity for entity in _step_solid_roots(entities))
    solids=Dict{Int,StepEntity}()
    for child in children
        info=representations[child]
        source=[item for item in info.items if haskey(solid_roots,item.id)]
        length(source)==1 && all(item->item.kind=="AXIS2_PLACEMENT_3D" || item.id==only(source).id,info.items) ||
            throw(ArgumentError("import_step: each assembly leaf must bind one actual closed solid"))
        solids[child]=only(source)
    end
    Set(keys(solid_roots))==Set(solid.id for solid in values(solids)) ||
        throw(ArgumentError("import_step: assembly has unplaced closed-solid definitions"))
    placements=Dict{Int,NamedTuple}();relations=Set{Int}();transforms=Set{Int}()
    for placement in selected("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
        placement.kind=="CONTEXT_DEPENDENT_SHAPE_REPRESENTATION" && length(placement.args)==2 ||
            throw(ArgumentError("import_step: assembly placement binding is truncated"))
        shape=_step_relationship_entity(entities,placement.args[2],"PRODUCT_DEFINITION_SHAPE")
        length(shape.args)==3 && all(value->value isa AbstractString,shape.args[1:2]) ||
            throw(ArgumentError("import_step: assembly placement product shape is truncated"))
        usage=_step_relationship_entity(entities,shape.args[3],"NEXT_ASSEMBLY_USAGE_OCCURRENCE")
        usage.id in (item.id for item in usages) && !haskey(placements,usage.id) ||
            throw(ArgumentError("import_step: duplicate or unselected assembly placement usage"))
        relation=_deref_entity(entities,placement.args[1])
        relation!==nothing || throw(ArgumentError("import_step: assembly relationship is missing"))
        parts=_step_relationship_parts(relation)
        Set(keys(parts))==Set(("REPRESENTATION_RELATIONSHIP",
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION","SHAPE_REPRESENTATION_RELATIONSHIP")) ||
            throw(ArgumentError("import_step: assembly needs an actual transformed shape relationship"))
        relationship=parts["REPRESENTATION_RELATIONSHIP"].args
        transformed=parts["REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"].args
        length(relationship)==4 && all(value->value isa AbstractString,relationship[1:2]) &&
            length(transformed)==1 && isempty(parts["SHAPE_REPRESENTATION_RELATIONSHIP"].args) ||
            throw(ArgumentError("import_step: assembly relationship is truncated"))
        child=usage.args[5][2];info=representations[child]
        relationship[3]==(:ref,info.representation.id) &&
            relationship[4]==(:ref,root_info.representation.id) ||
            throw(ArgumentError("import_step: relationship direction disagrees with parent/child usage"))
        info.context!=root_info.context && info.units==root_info.units ||
            throw(ArgumentError("import_step: assembly needs distinct three-dimensional contexts with the same actual units"))
        transform=_step_relationship_entity(entities,transformed[1],"ITEM_DEFINED_TRANSFORMATION")
        length(transform.args)==4 && transform.args[1] isa AbstractString &&
            (transform.args[2]===nothing || transform.args[2] isa AbstractString) ||
            throw(ArgumentError("import_step: assembly item transformation is truncated"))
        source=_step_relationship_entity(entities,transform.args[3],"AXIS2_PLACEMENT_3D")
        target=_step_relationship_entity(entities,transform.args[4],"AXIS2_PLACEMENT_3D")
        source.id in (item.id for item in info.items) && target.id in (item.id for item in root_info.items) ||
            throw(ArgumentError("import_step: assembly transform frames are not literal representation members"))
        local_transform=_step_occurrence_transform(entities,transform.args[3],transform.args[4])
        rigid=_step_occurrence_compose(_STEP_OCCURRENCE_IDENTITY,local_transform)
        placements[usage.id]=(;usage,child,rigid)
        push!(shape_definitions,shape.id);push!(relations,relation.id);push!(transforms,transform.id)
    end
    Set(keys(placements))==Set(usage.id for usage in usages) &&
        Set(entity.id for entity in selected("PRODUCT_DEFINITION_SHAPE"))==shape_definitions &&
        Set(entity.id for entity in selected("ITEM_DEFINED_TRANSFORMATION"))==transforms &&
        Set(entity.id for entity in selected("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"))==relations &&
        Set(entity.id for entity in selected("REPRESENTATION_RELATIONSHIP"))==relations &&
        Set(entity.id for entity in selected("SHAPE_REPRESENTATION_RELATIONSHIP"))==relations ||
        throw(ArgumentError("import_step: assembly placement inventory is incomplete or orphaned"))
    records=NamedTuple[];surfaces=Dict{Int,Mesh}();node_count=0;triangle_count=0
    root_product=products[root]
    for usage in usages
        placement=placements[usage.id];solid=solids[placement.child]
        surface=get!(surfaces,solid.id) do
            _step_root_surface(solid,entities)
        end
        node_count=_checked_brep_add(node_count,nnodes(surface),0,"import_step","assembly occurrence nodes")
        triangle_count=_checked_brep_add(triangle_count,size(surface.tris,2),0,"import_step","assembly occurrence triangles")
        node_count<=_BREP_MAX_SURFACE_NODES && triangle_count<=_BREP_MAX_SURFACE_TRIS ||
            throw(ArgumentError("import_step: product usages exceed assembly surface storage limits"))
        rotation,translation=placement.rigid
        coordinates=Matrix{Float64}(undef,size(surface.coords))
        for column in axes(coordinates,2),i in 1:3
            coordinates[i,column]=translation[i]+sum(rotation[i][j]*surface.coords[j,column] for j in 1:3)
        end
        all(isfinite,coordinates) || throw(ArgumentError("import_step: assembly occurrence coordinates are nonfinite"))
        product=products[placement.child]
        name="STEP PRODUCT #$(root_product.id) $(root_product.args[2])/NAUO #$(usage.id) $(usage.args[2])/" *
            "PRODUCT #$(product.id) $(product.args[2])/$(solid.kind) #$(solid.id)"
        push!(records,(;name,surface=Mesh(coordinates;tris=surface.tris)))
    end
    length(records)>1 && return _brep_assembly(records,fill_solids,"import_step")
    surface=only(records).surface
    validate(surface).ok || throw(ArgumentError("import_step: invalid product-usage solid"))
    _brep_assembly_shell_certificate(surface.coords,surface.tris,"import_step")
    return fill_solids ? _filled(surface,"import_step") : surface
end
