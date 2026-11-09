function _model_import_xao_geometry!(model,path,caller)
    metadata=IO.read_xao(path)
    text=read(path,String)
    header=match(r"<shape\b([^>]*)>",text)
    header===nothing && throw(ArgumentError("$caller: XAO has no shape"))
    attrs=IO._xao_attrs(header.captures[1])
    format=IO._xao_unescape(get(attrs,"format",""))
    format=="BREP" || throw(ArgumentError("$caller: unsupported XAO shape format $(repr(format))"))
    external=IO._xao_unescape(get(attrs,"file",""))
    source=if isempty(external)
        shape=match(r"(?s)<shape\b[^>]*>(.*?)</shape>",text)
        shape===nothing && throw(ArgumentError("$caller: XAO has no shape data"))
        cdata=match(r"(?s)^\s*<!\[CDATA\[(.*)\]\]>\s*$",shape.captures[1])
        cdata===nothing ? IO._xao_unescape(shape.captures[1]) : String(cdata.captures[1])
    else
        shape_path=isabspath(external) ? external : joinpath(dirname(path),external)
        isfile(shape_path) || throw(ArgumentError("$caller: missing XAO shape file $shape_path"))
        read(shape_path,String)
    end
    kernel=model.occ_kernel
    retained=Set(keys(kernel.shapes))
    offsets=ntuple(d->max(kernel.maxima[d+2],model.next_tag[d]),4)
    receipt=Ref{Any}(nothing)
    Base.invokelatest(Model._model_occ_import_brep_source,model,String(source);
        caller=caller,topology_receipt=receipt)
    kernel=model.occ_kernel
    imported=setdiff(Set(keys(kernel.shapes)),retained)
    references=receipt[]
    targets=Dict{UInt64,Tuple{Int,Int}}()
    indices=Dict{Tuple{Int,Int},UInt64}()
    for entity in metadata.entities
        key=(entity.dim,entity.index)
        haskey(indices,key) && throw(ArgumentError("$caller: duplicate XAO topology index $key"))
        bound=get(references,entity.reference,nothing)
        bound!==nothing && bound[1]==entity.dim || throw(ArgumentError(
            "$caller: XAO topology reference $(entity.reference) has the wrong dimension"))
        id=kernel.bindings[bound]
        haskey(targets,id) && throw(ArgumentError("$caller: duplicate XAO geometry reference"))
        entity.index<=typemax(Int32)-offsets[entity.dim+1]-1 || throw(ArgumentError(
            "$caller: XAO entity tags exceed Int32"))
        targets[id]=(entity.dim,offsets[entity.dim+1]+entity.index+1)
        indices[key]=id
    end
    Base.invokelatest(Model._occ_brep_retag_entities!,kernel,targets,imported;caller=caller)
    Base.invokelatest(Model._model_occ_synchronize!,model)
    for entity in metadata.entities
        isempty(entity.name) && continue
        dimension,tag=kernel.shape_tags[indices[(entity.dim,entity.index)]]
        Model.set_entity_name!(model,dimension,tag,entity.name)
    end
    for group in metadata.groups
        tags=Int[]
        for index in group.tags
            id=get(indices,(group.dim,index-1),nothing)
            id===nothing && throw(ArgumentError("$caller: XAO group references an unknown topology index"))
            push!(tags,kernel.shape_tags[id][2])
        end
        unique!(tags)
        isempty(tags) && continue
        key=(group.dim,group.tag)
        if haskey(model.physical,key)
            for tag in tags
                tag in model.physical[key] && continue
                push!(model.physical[key],tag)
                push!(get!(model.entity_physicals,(group.dim,tag),Int[]),group.tag)
            end
        else
            Model.add_physical_group!(model,group.dim,tags;tag=group.tag,name=group.name)
        end
    end
    return nothing
end
