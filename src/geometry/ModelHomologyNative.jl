function _native_homology_entities(m,physical_tags,top)
    entities=Tuple{Int,Int}[]
    if isempty(physical_tags)
        for entity in model_entities(m)
            entity[1]==top && push!(entities,entity)
        end
        return entities
    end
    for tag in physical_tags
        for dim in 0:3
            members=get(m.physical,(dim,tag),nothing)
            members===nothing && continue
            append!(entities,[(dim,entity) for entity in sort(members)])
        end
    end
    return entities
end

function _native_homology_elements(cells,entities)
    elements=Tuple{Int32,Vector{UInt64}}[]
    for entity in entities
        append!(elements,get(cells,entity,Tuple{Int32,Vector{UInt64}}[]))
    end
    return elements
end

function _native_homology_reverse_nodes(cell)
    n=length(cell.nodes)
    slots=cell.dim==0 ? (1,) : cell.dim==1 ? (2,1) :
        cell.dim==2 ? (n==3 ? (1,3,2) : (1,4,3,2)) :
        n==4 ? (2,1,3,4) : n==5 ? (3,2,1,4,5) :
        n==6 ? (2,1,3,5,4,6) : (3,2,1,4,7,6,5,8)
    return UInt64[cell.nodes[i] for i in slots]
end

function _native_homology_publish_chain!(m,complex,chain,dim,name,element_counter,
                                         element_counter_wrapped=Ref(false))
    entity_tag=maximum(last,model_entities(m);init=0)+1
    add_discrete_entity!(m,dim,entity_tag)
    push!(m.meshing.chain_entities,(dim,entity_tag))
    record=m.discrete[(dim,entity_tag)]
    ids=sort!(collect(keys(chain));by=id->_native_homology_key(complex.cells[id].nodes))
    for id in ids
        coefficient=chain[id];coefficient==0 && continue
        cell=complex.cells[id]
        cell.combined && throw(ArgumentError(
            "compute_homology!: representative contains an unexpanded combined cell"))
        nodes=coefficient<0 && dim>0 ? _native_homology_reverse_nodes(cell) : cell.nodes
        for copy_index in 1:abs(coefficient)
            element_counter[]+=UInt64(1)
            element_counter_wrapped[]=element_counter_wrapped[] || iszero(element_counter[])
            _record_append_element!(record,_native_homology_cell_type(cell),
                element_counter[],copy(nodes))
        end
    end
    physical=add_physical_group!(m,dim,[entity_tag];name=name)
    m.physical_names[(dim,physical)]=name
    return dim,physical
end

function _native_homology_error(error_callback,message,legacy_message)
    error_callback===nothing && throw(ArgumentError("compute_homology!: $legacy_message"))
    error_callback(message)
    return nothing
end

function _native_compute_homology!(m,cells,requests;coordinates=nothing,
                                   initial_max_element_tag=UInt64(0),
                                   initial_element_counter_wrapped=false,
                                   element_counter_report=nothing,error_callback=nothing)
    caller="compute_homology!";output=Tuple{Int,Int}[]
    initial_element_counter_wrapped isa Bool || throw(ArgumentError(
        "$caller: initial element counter wrapped flag must be Bool"))
    element_counter_report===nothing ||
        element_counter_report isa Base.RefValue{Tuple{UInt64,Bool}} || throw(ArgumentError(
            "$caller: constructor clock report must be Ref{Tuple{UInt64,Bool}}"))
    initial_max_element_tag isa Integer && !(initial_max_element_tag isa Bool) &&
        0<=initial_max_element_tag<=typemax(UInt64) || throw(ArgumentError(
            "$caller: initial element maximum must be a nonnegative Int32 integer"))
    maximum_element=UInt64(initial_max_element_tag)
    # A wrapped constructor clock can be lower than surviving element UIDs.
    if !initial_element_counter_wrapped
        for (_,record) in _discrete_mesh_records_model(m)
            maximum_element=max(maximum_element,UInt64(maximum(record.element_tags;init=UInt64(0))))
        end
        authority=m.meshing.mesh_authority
        if authority!==nothing && authority.mesh!==nothing
            for tags in authority.mesh.entity_data.external_element_tags
                maximum_element=max(maximum_element,maximum(tags;init=UInt64(0)))
            end
        end
    end
    element_counter=Ref(maximum_element)
    element_counter_wrapped=Ref(initial_element_counter_wrapped)
    pairs=Tuple{Vector{Int},Vector{Int}}[]
    for request in requests
        key=(request.domain,request.subdomain)
        key in pairs || push!(pairs,key)
    end
    sort!(pairs;by=pair->(Tuple(pair[1]),Tuple(pair[2])))
    for key in pairs
        group=[request for request in requests if (request.domain,request.subdomain)==key]
        initial=group[1];top=model_dimension(m)
        domain=_native_homology_entities(m,initial.domain,top)
        subdomain=isempty(initial.domain) || isempty(initial.subdomain) ? Tuple{Int,Int}[] :
            _native_homology_entities(m,initial.subdomain,top)
        domain_elements=_native_homology_elements(cells,domain)
        subdomain_elements=_native_homology_elements(cells,subdomain)
        measure_geometry=any(request->request.kind=="Cohomology",group)
        complex=Ref{Union{Nothing,_NativeHomologyComplex}}(nothing)
        computed_h=Set{Int}();computed_c=Set{Int}()
        for request in group
            kind=request.kind
            if !(kind in ("Homology","Cohomology","Betti"))
                _native_homology_error(error_callback,
                    "Unknown type of homology computation: $kind",
                    "unsupported computation type \"$kind\"")
                continue
            end
            dimensions=isempty(request.dims) ? collect(0:top-1) :
                Int[dim for dim in request.dims if 0<=dim<=top]
            if isempty(dimensions)
                _native_homology_error(error_callback,
                    "Invalid homology computation dimensions given",
                    "invalid homology computation dimensions given")
                continue
            end
            computed=kind=="Cohomology" ? computed_c : computed_h
            kind!="Betti" && all(dim->dim in computed,dimensions) && continue
            if complex[]===nothing
                isempty(domain) && _native_homology_error(error_callback,
                    "Domain is empty","domain is empty")
                bound_coordinates=coordinates===nothing && measure_geometry ?
                    first(_model_discrete_curve_sources(m,caller)) : coordinates
                built=_native_homology_complex(domain_elements,subdomain_elements,
                    bound_coordinates;measure_geometry)
                isempty(built.active[1]) && _native_homology_error(error_callback,
                    "Cell Complex is empty: check the domain and the mesh",
                    "cell complex is empty: check the domain and the mesh")
                complex[]=built
            end
            current=complex[]
            if kind=="Betti"
                _native_homology_restore!(current)
                _native_homology_reduce!(current)
                for dim in 0:2
                    _native_homology_combine!(current,dim;dual=true)
                end
                _native_homology_integer_basis(current,false)
                continue
            end
            union!(computed,dimensions)
            chains=_native_homology_basis_chains!(current,dimensions;dual=kind=="Cohomology")
            domain_name=_homology_domain_string(request.domain,request.subdomain)
            for dim in dimensions
                for (index,chain) in enumerate(chains[dim+1])
                    prefix=kind=="Homology" ? "H_" : "H^"
                    name=prefix*string(dim)*domain_name*string(index)
                    push!(output,_native_homology_publish_chain!(m,current,chain,dim,name,
                                                               element_counter,element_counter_wrapped))
                end
            end
        end
    end
    element_counter_report===nothing ||
        (element_counter_report[]=(element_counter[],element_counter_wrapped[]))
    return output
end

# Completed products already bind every vertex by its public label. Retain
# those labels and cell-slot order directly; no geometry carrier row or
# Cartesian coordinate participates in vertex identity.
function _native_homology_actual_cells(mesh::MixedMesh)
    data=mesh.entity_data
    data===nothing && throw(ArgumentError(
        "compute_homology!: actual mesh requires entity metadata"))
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    sizehint!(coordinates,nnodes(mesh))
    for (row,tag) in enumerate(data.external_node_tags)
        coordinates[tag]=ntuple(axis->mesh.coords[axis,row],3)
    end
    cells=Dict{Tuple{Int,Int},Vector{Tuple{Int32,Vector{UInt64}}}}()
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock || throw(ArgumentError(
            "compute_homology!: unsupported actual cell block"))
        dimension=msh_dimension(block.msh)
        for column in axes(block.nodes,2)
            owner=(dimension,Int(data.block_entities[bi][column]))
            nodes=copy(data.external_node_tags[block.nodes[:,column]])
            push!(get!(()->Tuple{Int32,Vector{UInt64}}[],cells,owner),(Int32(block.msh),nodes))
        end
    end
    return cells,coordinates
end
