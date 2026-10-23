# Model data retains actual UIDs independently of its exported mesh numbering.
# All rows and headers are validated before any file publication; a live-object
# lookup failure has its own literal header/row prefix, separate from validation.
using LinearAlgebra: I
using Printf: @printf

function _native_view_writer_groups(view,path::AbstractString)
    slots=Any[]
    for range in view.model_step_ranges
        any(slot->slot===range.model_slot,slots) || push!(slots,range.model_slot)
    end
    root,extension=splitext(path)
    return [(slot=slot,path=length(slots)==1 ? String(path) : root*"_"*string(i-1)*extension,
        keys=sort!([key for key in keys(view.model_contexts)
            if view.model_contexts[key].model_slot===slot])) for (i,slot) in enumerate(slots)]
end

function _native_view_writer_mesh(slot,caller)
    m=slot.model;cached=slot.mesh;class=slot.class
    cached,class=_dim01_actual_cache(cached,class;physical_names=m.physical_names)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(0),UInt64(0),Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    _dim01_ingest!(assembly,m,cached,class;preserve_source_order=true)
    cells=sort!(collect(eachindex(assembly.cells));by=i->_dim01_cell_order(assembly.cells[i]),
        alg=Base.Sort.MergeSort)
    rows=sort!(collect(eachindex(assembly.nodes));by=i->assembly.nodes[i].owner,
        alg=Base.Sort.MergeSort)
    inverse=zeros(Int32,length(rows))
    for (i,row) in enumerate(rows);inverse[row]=Int32(i);end
    coordinates=Matrix{Float64}(undef,3,length(rows))
    for (i,row) in enumerate(rows),axis in 1:3;coordinates[axis,i]=assembly.nodes[row].point[axis];end
    node_tags=UInt64[assembly.nodes[row].tag for row in rows]
    owners=[assembly.nodes[row].owner for row in rows]
    parameters=[Model._model_mesh_projected_parameter(assembly.nodes[row].owner,
        assembly.nodes[row].parameter) for row in rows]
    families=unique(assembly.cells[i].msh for i in cells)
    blocks=ElementBlock[];block_owners=Vector{Int32}[];element_tags=Vector{UInt64}[]
    for msh in families
        selected=[i for i in cells if assembly.cells[i].msh==msh]
        nodes=Matrix{Int32}(undef,Elements.msh_num_nodes(msh),length(selected))
        for (column,i) in enumerate(selected);nodes[:,column]=inverse[assembly.cells[i].nodes];end
        push!(blocks,ElementBlock(msh,nodes,Int32[assembly.cells[i].physical for i in selected]))
        push!(block_owners,Int32[assembly.cells[i].owner[2] for i in selected])
        push!(element_tags,UInt64[assembly.cells[i].tag for i in selected])
    end
    entities,_=_dim01_entities(m,model_entities(m),assembly,rows)
    partition=nothing
    if cached isa MixedMesh && cached.partition_data!==nothing
        source=cached.partition_data;actual=cached.entity_data
        actual!==nothing || throw(ArgumentError("$caller: partition UID authority is absent"))
        partition_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
        for (key,previous) in source.entities
            current=get(entities,key,nothing)
            current!==nothing || throw(ArgumentError("$caller: partition entity $key is no longer present"))
            referenced=Set{Int32}()
            for (bi,block) in enumerate(blocks),column in axes(block.nodes,2)
                (Elements.msh_dimension(block.msh),Int(block_owners[bi][column]))==key || continue
                union!(referenced,block.nodes[:,column])
            end
            for row in eachindex(owners);owners[row]==(key[1],Int32(key[2])) && push!(referenced,Int32(row));end
            bounds=isempty(referenced) ? current.bbox : ntuple(axis->axis<=3 ?
                minimum(coordinates[axis,row] for row in referenced) :
                maximum(coordinates[axis-3,row] for row in referenced),6)
            partition_entities[key]=Elements.MixedEntity(key[1],key[2],bounds;
                physical_tags=current.physical_tags,boundaries=current.boundaries,
                embedded_curves=current.embedded_curves,
                partitions=previous.partitions,parent=previous.parent)
            delete!(entities,key)
        end
        source_tags=reduce(vcat,actual.external_element_tags;init=UInt64[])
        emitted_tags=reduce(vcat,element_tags;init=UInt64[])
        positions=Dict(uid=>index for (index,uid) in enumerate(emitted_tags))
        ghost_records=Elements.MixedGhostElement[]
        for record in source.ghost_elements
            uid=source_tags[record.element]
            haskey(positions,uid) || throw(ArgumentError("$caller: ghost element UID $uid is no longer present"))
            push!(ghost_records,Elements.MixedGhostElement(positions[uid],record.partition,record.ghost_partitions))
        end
        partition=Elements.MixedPartitionData(source.num_partitions;
            entities=partition_entities,ghost_entities=source.ghost_entities,ghost_elements=ghost_records)
    end
    metadata=Elements._native_mixed_entity_data(entities;node_entities=owners,
        node_parametric=parameters,external_node_tags=node_tags,block_entities=block_owners,
        external_element_tags=element_tags)
    return MixedMesh(coordinates,blocks;entity_data=metadata,elementary_entities=block_owners,
        physical_names=m.physical_names,partition_data=partition)
end

function _native_view_writer_section(view,key;scheme::Bool,name=view.name)
    tags,values,components=view.model_data[key];context=view.model_contexts[key]
    counts=context.counts
    components>0 && length(tags)==length(counts) || throw(ArgumentError(
        "model view writer: incomplete row metadata"))
    all(count->count>0 && count%components==0,counts) || throw(ArgumentError(
        "model view writer: values do not form complete component rows"))
    kind=key[2];nodal=kind=="ElementNodeData"
    row_nodes=nodal ? Int32.(counts.÷components) : fill(Int32(1),length(tags))
    nodal || all(==(components),counts) || throw(ArgumentError(
        "model view writer: non-nodal row widths differ from the component count"))
    strings=String["\""*name*"\""]
    scheme && kind!="NodeData" && push!(strings,"\"INTERPOLATION_SCHEME\"")
    header=_native_view_partition_header(key[1],components,length(tags),
        get(context,:partitions,Set{Int}()))
    return Elements.MshUIDDataSection(kind,strings,[context.time],header,tags,
        row_nodes,UInt64[],nodal,values,4;logical_extent=context.extent)
end

function _native_view_writer_rows(section,node_map,element_map;binary::Bool)
    Elements._msh_uid_validate(section)
    mapping=section.name=="NodeData" ? node_map : element_map
    io=IOBuffer()
    println(io,'$',section.name)
    println(io,length(section.strings));foreach(value->println(io,value),section.strings)
    println(io,length(section.reals));foreach(value->@printf(io,"%.16g\n",value),section.reals)
    println(io,length(section.header));foreach(value->println(io,value),section.header)
    offset=1;components=Int(section.header[2])
    for row in eachindex(section.tags)
        uid=section.tags[row]
        if !haskey(mapping,uid)
            what=section.name=="NodeData" ? "node" : "element"
            return (bytes=take!(io),error="Unknown $what $uid in data")
        end
        emitted=mapping[uid];nodes=Int(section.row_nodes[row]);width=nodes*components
        if binary
            emitted<=typemax(Int32) || throw(ArgumentError(
                "model view writer: binary data UID exceeds Int32"))
            write(io,Int32(emitted))
            section.name=="ElementNodeData" && write(io,Int32(nodes))
            write(io,@view(section.values[offset:offset+width-1]))
        else
            print(io,emitted)
            section.name=="ElementNodeData" && print(io,' ',nodes)
            for value in @view(section.values[offset:offset+width-1]);@printf(io," %.16g",value);end
            println(io)
        end
        offset+=width
    end
    binary && println(io)
    println(io,"\$End",section.name)
    return (bytes=take!(io),error=nothing)
end

# The ordered polynomial lattice defines V*C=I. A fixed LU solve preserves
# floating arithmetic, including the signed zeros in interpolation coefficients.
function _native_view_writer_basis(msh::Int)
    spec=msh_spec(msh);family=spec.family;p=spec.order
    lattice=family===:pnt ? [(0,)] : family===:lin ? Elements._line_monomials(p) :
        family===:tri ? Elements._tri_monomials(p,spec.serendipity) :
        family===:qua ? Elements._qua_monomials(p,spec.serendipity) :
        family===:tet ? Elements._tet_monomials(p,spec.serendipity) :
        family===:hex ? Elements._hex_monomials(p,spec.serendipity) :
        family===:pri ? Elements._pri_monomials(p,spec.serendipity) : nothing
    spec.serendipity && family in (:qua,:hex,:pri) &&
        (lattice=MeshFunctionSpaces._serendipity_exponents(family,p))
    if family===:pyr && p==1
        # A collapsed pyramid uses (1,w,v,u,uv/(1-w)); the family code supplies
        # the collapse denominator in its interpretation of the exponent table.
        return (family=6,coefficients=[.25 -.25 -.25 -.25 .25;
            .25 -.25 -.25 .25 -.25;.25 -.25 .25 .25 .25;
            .25 -.25 .25 -.25 -.25;0. 1. 0. 0. 0.],
            exponents=[0 0 0;0 0 1;0 1 0;1 0 0;1 1 0])
    end
    family===:pyr && p>=2 && return _native_view_writer_pyramid_basis(msh)
    lattice===nothing && throw(ArgumentError(
        "model view writer: interpolation matrices for this family are not implemented"))
    length(lattice)==spec.nnodes || throw(ArgumentError(
        "model view writer: reference lattice does not cover the element"))
    dimension=max(1,Elements.msh_dimension(msh));nodes=Elements.lagrange_nodes(msh)
    exponents=Matrix{Int}(undef,length(lattice),dimension)
    for row in eachindex(lattice),axis in 1:dimension
        exponents[row,axis]=axis<=length(lattice[row]) ? lattice[row][axis] : 0
    end
    matrix=Matrix{Float64}(undef,length(lattice),length(lattice))
    for row in axes(matrix,1),column in axes(matrix,2)
        matrix[row,column]=prod(nodes[axis,column]^exponents[row,axis]
            for axis in 1:dimension)
    end
    coefficients=matrix\Matrix{Float64}(I,length(lattice),length(lattice))
    code=family===:pnt ? 1 : family===:lin ? 2 : family===:tri ? 3 :
        family===:qua ? 4 : family===:tet ? 5 : family===:pri ? 7 : 8
    return (family=code,coefficients,exponents)
end

function _native_view_writer_scheme(mesh;types=nothing,kind=nothing)
    selected=Dict{Int,NamedTuple}()
    current=types===nothing ? (Int(block.msh) for block in mesh.blocks if !isempty(block.nodes)) : types
    for msh in current
        spec=msh_spec(msh)
        basis=if kind=="ElementData" && spec.order>1
            family=spec.family
            code=family===:pnt ? 1 : family===:lin ? 2 : family===:tri ? 3 :
                family===:qua ? 4 : family===:tet ? 5 : family===:pyr ? 6 : family===:pri ? 7 : 8
            (family=code,coefficients=ones(1,1),exponents=zeros(Int,1,3))
        else
            _native_view_writer_basis(Int(msh))
        end
        selected[basis.family]=basis
    end
    isempty(selected) && return UInt8[]
    io=IOBuffer();println(io,"\$InterpolationScheme\n\"INTERPOLATION_SCHEME\"")
    println(io,length(selected))
    for code in sort!(collect(keys(selected)))
        basis=selected[code];println(io,code);println(io,2)
        for matrix in (basis.coefficients,basis.exponents)
            println(io,size(matrix,1),' ',size(matrix,2))
            for row in axes(matrix,1)
                for column in axes(matrix,2);@printf(io,"%.16g ",matrix[row,column]);end
                println(io)
            end
        end
    end
    println(io,"\$EndInterpolationScheme")
    return take!(io)
end

function _native_view_writer_ascii_mesh(io,mesh,version,names)
    mesh.partition_data===nothing || version==4.1 || throw(ArgumentError(
        "model view writer: partition and ghost metadata require MSH 4.1"))
    println(io,"\$MeshFormat");println(io,version," 0 8");println(io,"\$EndMeshFormat")
    Elements._write_mixed_physical_names(io,names,true)
    data=mesh.entity_data
    if version==4.1
        entities=sort!(collect(values(data.entities));by=e->(e.dim,e.tag))
        println(io,"\$Entities")
        println(io,join([count(e->e.dim==dim,entities) for dim in 0:3],' '))
        max_curve=maximum((Int(e.tag) for e in entities if e.dim==1);init=0)
        for entity in entities
            print(io,entity.tag,' ')
            for value in entity.bbox[1:(entity.dim==0 ? 3 : 6)];@printf(io,"%.16g ",value);end
            print(io,length(entity.physical_tags),' ')
            for tag in entity.physical_tags;print(io,tag,' ');end
            if entity.dim>0
                boundaries=Elements._mixed_entity_encoded_boundaries(entity,max_curve)
                print(io,length(boundaries),' ')
                for tag in boundaries;print(io,tag,' ');end
            end
            println(io)
        end
        println(io,"\$EndEntities")
        mesh.partition_data===nothing || Elements._write_mixed_partitioned_entities(
            io,mesh.partition_data,false,true)
    end
    println(io,"\$Nodes")
    if version==2.2
        println(io,size(mesh.coords,2))
        for row in axes(mesh.coords,2)
            @printf(io,"%d %.16g %.16g %.16g\n",row,mesh.coords[1,row],mesh.coords[2,row],mesh.coords[3,row])
        end
    else
        runs=Elements._mixed_metadata_node_runs(data)
        tags=data.external_node_tags
        println(io,length(runs),' ',length(tags),' ',minimum(tags;init=typemax(UInt64)),' ',maximum(tags;init=UInt64(0)))
        for run in runs
            row=first(run);dim,owner=data.node_entities[row]
            parametric=data.node_parametric[row]!==nothing
            println(io,dim,' ',owner,' ',Int(parametric),' ',length(run))
            for row in run;println(io,tags[row]);end
            for row in run
                @printf(io,"%.16g %.16g %.16g",mesh.coords[1,row],mesh.coords[2,row],mesh.coords[3,row])
                if parametric
                    for value in data.node_parametric[row];@printf(io," %.16g",value);end
                end
                println(io)
            end
        end
    end
    println(io,"\$EndNodes\n\$Elements")
    if version==2.2
        cell_count=Elements._assert_mixed_structure(mesh,"model view writer")
        println(io,cell_count);order,emitted=Elements._mixed_v2_order(mesh)
        ids=Elements._v2_entity_ids(mesh)
        for ref in order
            block=mesh.blocks[ref.block];column=ref.cell
            owner=Elements._v2_cell_entity(mesh,ids,ref.block,column)
            print(io,emitted[ref.block][column],' ',block.msh," 2 ",block.tags[column],' ',owner)
            for node in block.nodes[:,column];print(io,' ',node);end
            println(io)
        end
    else
        runs=Elements._mixed_metadata_element_runs(mesh,data)
        tags=UInt64[uid for block in data.external_element_tags for uid in block]
        println(io,length(runs),' ',length(tags),' ',minimum(tags;init=typemax(UInt64)),' ',maximum(tags;init=UInt64(0)))
        for run in runs
            block=mesh.blocks[run.block]
            owner=data.block_entities[run.block][first(run.cells)]
            println(io,Elements.msh_dimension(block.msh),' ',owner,' ',block.msh,' ',length(run.cells))
            for column in run.cells
                print(io,data.external_element_tags[run.block][column],' ')
                for node in block.nodes[:,column];print(io,data.external_node_tags[node],' ');end
                println(io)
            end
        end
    end
    println(io,"\$EndElements")
    (version==4.1 && mesh.partition_data!==nothing) && Elements._write_mixed_ghost_elements(io,
        mesh.partition_data,reduce(vcat,data.external_element_tags;init=UInt64[]),false,true)
    return nothing
end

function _native_view_writer_payload(view,group,mesh;version,binary::Bool,save_mesh::Bool,scheme::Bool,
        lookup=nothing,name=view.name)
    interpolation=scheme ? _native_view_writer_scheme(mesh;
        types=_native_view_interpolation_history(view.model_contexts),kind=view.model_kind) : UInt8[]
    sections=[_native_view_writer_section(view,key;scheme=!isempty(interpolation),name)
              for key in group.keys]
    data=mesh.entity_data
    data!==nothing || throw(ArgumentError("model view writer: actual mesh identity is absent"))
    source_nodes=data.external_node_tags
    source_cells=UInt64[uid for tags in data.external_element_tags for uid in tags]
    node_map=Dict(uid=>uid for uid in source_nodes);element_map=Dict(uid=>uid for uid in source_cells)
    if save_mesh && version==2.2
        node_map=Dict(uid=>UInt64(i) for (i,uid) in enumerate(source_nodes))
        _,emitted=Elements._mixed_v2_order(mesh)
        element_map=Dict(uid=>emitted[bi][column] for (bi,tags) in enumerate(data.external_element_tags)
            for (column,uid) in enumerate(tags))
    end
    if lookup!==nothing
        if view.model_kind=="NodeData"
            filter!(pair->get(lookup.nodes,first(pair),false),node_map)
        else
            filter!(pair->get(lookup.elements,first(pair),nothing)!==nothing,element_map)
        end
    end
    io=IOBuffer()
    if save_mesh
        names=Elements._copy_physical_names(mesh.physical_names,"model view writer")
        Elements._assert_mixed_msh_format(mesh,version,binary)
        if !binary
            _native_view_writer_ascii_mesh(io,mesh,version,names)
        elseif version==2.2
            Elements._write_mixed_v2_binary(io,mesh,names,true,nothing)
        else
            Elements._write_mixed_v4_binary(io,mesh,names,true,nothing,true)
            mesh.partition_data===nothing || Elements._write_mixed_ghost_elements(io,
                mesh.partition_data,source_cells,true,true)
        end
    else
        println(io,"\$MeshFormat");println(io,version,' ',Int(binary)," 8")
        binary && (write(io,Int32(1));println(io))
        println(io,"\$EndMeshFormat")
    end
    write(io,interpolation)
    for section in sections
        result=_native_view_writer_rows(section,node_map,element_map;binary)
        write(io,result.bytes)
        result.error===nothing || return (bytes=take!(io),error=result.error)
    end
    return (bytes=take!(io),error=nothing)
end

function _native_view_writer_lookup!(slot,mesh,kind)
    catalogue=get!(_NativeLookupCatalogue,NATIVE_LOOKUP_CATALOGUES,slot)
    data=mesh.entity_data
    if kind=="NodeData" && !catalogue.nodes_primed
        _native_diagnostic_emit_locked(:debug,"Rebuilding mesh node cache")
        catalogue.nodes=Dict(uid=>true for uid in data.external_node_tags)
        # The constructor high-water mark survives node retirement and renumbering.
        # It determines the native cache representation independently of live labels.
        maximum_tag=max(slot.node_tag_max,maximum(data.external_node_tags;init=UInt64(0)))
        count=length(data.external_node_tags)
        dense=maximum_tag==count || maximum_tag<10*UInt64(count)
        if maximum_tag==count
            _native_diagnostic_emit_locked(:debug,"We have a dense node numbering in the cache")
        elseif dense
            _native_diagnostic_emit_locked(:debug,"We have a fairly dense node numbering - still using cache vector")
        end
        catalogue.node_dense_max=dense ? maximum_tag : nothing
        catalogue.nodes_primed=dense || !isempty(catalogue.nodes)
    elseif kind!="NodeData" && !catalogue.elements_primed
        _native_diagnostic_emit_locked(:debug,"Rebuilding mesh element cache")
        for (bi,block) in enumerate(mesh.blocks),(column,uid) in enumerate(data.external_element_tags[bi])
            catalogue.elements[uid]=_NativeElementLookupBinding(Int32(block.msh),
                Elements.msh_dimension(block.msh),Int(data.block_entities[bi][column]))
        end
        catalogue.elements_primed=!isempty(catalogue.elements)
    end
    return catalogue
end

# Independent primary-corner serialization from owned current mesh identities
# and populated model-view rows. Allocated empty steps do not emit records.
function _native_view_text_table(view,key;name=view.name)
    section=_native_view_writer_section(view,key;scheme=false,name)
    Elements._msh_uid_validate(section)
    rows=Dict{UInt64,Vector{Float64}}();offset=1
    for (uid,count) in zip(section.tags,view.model_contexts[key].counts)
        width=Int(count);rows[uid]=copy(@view section.values[offset:offset+width-1])
        offset+=width
    end
    return (step=key[1],kind=key[2],components=Int(section.header[2]),
        time=section.reals[1],rows)
end

function _native_view_text_cells(mesh;primary_only=true)
    data=mesh.entity_data
    data!==nothing || throw(ArgumentError("model view writer: actual mesh identity is absent"))
    owners=sort!(collect(keys(data.entities)))
    ordinals=Dict(owner=>index-1 for (index,owner) in enumerate(owners))
    per_owner=Dict{Tuple{Int,Int32},Int}();cells=NamedTuple[]
    refs,_=Elements._mixed_v2_order(mesh)
    # Entity order precedes family vectors. Preserve their current column
    # order; numeric element UIDs do not define constructor chronology.
    function order(ref)
        block=mesh.blocks[ref.block];family=msh_spec(block.msh).family
        rank=family===:pnt ? 0 : family===:lin ? 1 : family===:tri ? 2 :
            family===:qua ? 3 : family===:tet ? 4 : family===:hex ? 5 :
            family===:pri ? 6 : family===:pyr ? 7 : 8
        return ((Elements.msh_dimension(block.msh),data.block_entities[ref.block][ref.cell]),rank)
    end
    sort!(refs;by=order,alg=Base.Sort.MergeSort)
    for ref in refs
        block=mesh.blocks[ref.block];column=ref.cell;spec=msh_spec(block.msh)
        owner=(Elements.msh_dimension(block.msh),data.block_entities[ref.block][column])
        haskey(ordinals,owner) || throw(ArgumentError("model view writer: cell owner is absent"))
        ordinal=get(per_owner,owner,0);per_owner[owner]=ordinal+1
        primary=Elements._msh_primary_nodes(spec.family)
        nodes=copy(@view block.nodes[1:(primary_only ? primary : size(block.nodes,1)),column])
        push!(cells,(uid=data.external_element_tags[ref.block][column],family=spec.family,
            msh=Int(block.msh),owner,physical=block.tags[column],
            nodes,node_uids=data.external_node_tags[nodes],entity=ordinals[owner],element=ordinal))
    end
    return cells
end

function _native_view_text_values(table,cell)
    components=table.components
    if table.kind=="NodeData"
        all(uid->haskey(table.rows,uid),cell.node_uids) || return nothing
        return [copy(@view table.rows[uid][1:components]) for uid in cell.node_uids]
    end
    haskey(table.rows,cell.uid) || return nothing
    values=table.rows[cell.uid]
    if table.kind=="ElementData"
        return [copy(@view values[1:components]) for _ in cell.nodes]
    end
    length(values)>=length(cell.nodes)*components || return nothing
    return [copy(@view values[components*(i-1)+1:components*i]) for i in eachindex(cell.nodes)]
end


function _native_view_writer_text_payload(view,plans;extension::AbstractString,
        name=view.name)
    extension in (".pos",".txt",".csv",".vtk") || throw(ArgumentError(
        "model view writer: this text format is not implemented"))
    meshes=IdDict{Any,Any}(plan.group.slot=>plan.mesh for plan in plans)
    tables=NamedTuple[]
    for key in sort!(collect(keys(view.model_contexts)))
        context=view.model_contexts[key]
        mesh=get(meshes,context.model_slot,nothing)
        mesh===nothing && throw(ArgumentError("model view writer: owned model is absent"))
        push!(tables,merge(_native_view_text_table(view,key;name),
            (mesh=mesh,cells=_native_view_text_cells(mesh),slot=context.model_slot)))
    end
    io=IOBuffer()
    if extension==".pos"
        if !isempty(tables) && any(table->table.slot!==first(tables).slot,tables)
            return (bytes=UInt8[],error="Cannot export multi-mesh datasets in .pos format")
        end
        println(io,"View \"",name,"\" {")
        if !isempty(tables)
            mesh=first(tables).mesh
            for cell in first(tables).cells
                values=[_native_view_text_values(table,cell) for table in tables]
                any(isnothing,values) && continue
                family=cell.family
                shape=family===:pnt ? 'P' : family===:lin ? 'L' : family===:tri ? 'T' :
                    family===:qua ? 'Q' : family===:tet ? 'S' : family===:hex ? 'H' :
                    family===:pri ? 'I' : family===:pyr ? 'Y' : nothing
                shape===nothing && throw(ArgumentError("model view writer: unknown POS cell family"))
                components=first(tables).components
                print(io,components==3 ? 'V' : components==9 ? 'T' : 'S',shape,'(')
                for (i,node) in enumerate(cell.nodes),axis in 1:3
                    (i!=1 || axis!=1) && print(io,',')
                    @printf(io,"%.16g",mesh.coords[axis,node])
                end
                print(io,"){");first_value=true
                for step_values in values,row_values in step_values,value in row_values
                    first_value || print(io,',')
                    @printf(io,"%.16g",value);first_value=false
                end
                println(io,"};")
            end
        end
        println(io,"};")
    else
        for table in tables,cell in table.cells
            values=_native_view_text_values(table,cell);values===nothing && continue
            for (node,row_values) in zip(cell.nodes,values)
                print(io,table.step,' ');@printf(io,"%.16g ",table.time)
                print(io,cell.entity,' ',cell.element,' ')
                for axis in 1:3;@printf(io,"%.16g ",table.mesh.coords[axis,node]);end
                for value in row_values;@printf(io,"%.16g ",value);end
            end
            println(io)
        end
    end
    return (bytes=take!(io),error=nothing)
end

function _native_view_writer_text_payload(view,group,mesh;extension::AbstractString,
        name=view.name)
    return _native_view_writer_text_payload(view,[(group=group,mesh=mesh)];
        extension,name)
end

# addListData stores one element at a time: all x coordinates, then all y,
# all z, and component/node values for every time step (Gmsh tutorial x3).
# The bound family scheme — not the stored type name — defines each flat
# record's geometry/value widths, matching PViewDataList's `_stat` decode.
function _native_list_view_records(view,caller)
    view=_view_data_owner(view)
    elements=Post.PosElement[]
    for type_key in view.list_types
        layout=Post._pos_list_record_layout(view.interp,type_key)
        layout===nothing && throw(ArgumentError("$caller: unsupported list data type $type_key"))
        ngeo,nval=layout.ngeo,layout.nval
        count=view.list_counts[type_key];data=view.list_data[type_key]
        if count==0
            isempty(data) || throw(ArgumentError("$caller: zero list element count has data"))
            continue
        end
        length(data)%count==0 || throw(ArgumentError("$caller: list data do not align with numEle"))
        width=length(data)÷count;value_width=width-3ngeo
        value_width>0 && value_width%(nval*layout.ncomp)==0 || throw(ArgumentError(
            "$caller: list data do not contain complete time-step values"))
        steps=value_width÷(layout.nval*layout.ncomp)
        for element in 0:count-1
            offset=element*width
            coords=Matrix{Float64}(undef,3,ngeo)
            for axis in 1:3,node in 1:ngeo
                coords[axis,node]=data[offset+(axis-1)*ngeo+node]
            end
            values=reshape(copy(@view(data[offset+3ngeo+1:offset+width])),layout.ncomp,nval,steps)
            push!(elements,Post.PosElement(layout.kind,coords,values;
                suffix=layout.suffix,scheme=layout.scheme))
        end
    end
    return elements
end

function _native_list_view_pos(view,caller;name=view.name)
    view=_view_data_owner(view)
    elements=_native_list_view_records(view,caller)
    # T2/T3 string blocks (POS text records and Annotate annotations) decode
    # back to PosText records: one list entry per record, '\0'-joined
    # per-step strings; 3 coords per record carry x/y/z with implicit style
    # zero, 4 carry x/y/z/style. Other dataType keys stay unexported.
    texts=Post.PosText[]
    for key in view.list_string_types
        key in ("T2","T3") || continue
        count=view.list_string_counts[key]
        coords=get(view.list_string_coords,key,Float64[])
        strings=get(view.list_strings,key,String[])
        length(strings)==count || throw(ArgumentError(
            "$caller: $key string entries do not align with the record count"))
        count==0 && (isempty(coords) || throw(ArgumentError(
            "$caller: $key record count is zero but coords are nonempty")))
        count==0 && continue
        length(coords)%count==0 || throw(ArgumentError(
            "$caller: $key coords do not align with the record count"))
        width=div(length(coords),count)
        width in (3,4) || throw(ArgumentError(
            "$caller: $key coords must hold 3 or 4 entries per record " *
            "(got $width)"))
        for record in 1:count
            base=width*(record-1)
            xyz=(coords[base+1],coords[base+2],coords[base+3])
            style=width==4 ? coords[base+4] : 0.0
            push!(texts,Post.PosText(key=="T2" ? 2 : 3,xyz,style,
                                     split(strings[record],'\0')))
        end
    end
    return Post.PosView(name;elements,texts,time=view.list_time)
end


# The exported step table: declared physical times, padded by step index when
# shorter (upstream's implicit 0,1,2… numbering), capped at the smallest
# record's step count like `NbTimeStepLIST`. Records must share a common step
# count, as in upstream `PViewDataList` (`different steps` views cannot be
# represented).
function _native_list_view_times(view,elements)
    view=_view_data_owner(view)
    steps=isempty(elements) ? 0 : size(first(elements).values,3)
    all(element->size(element.values,3)==steps,elements) ||
        throw(ArgumentError("API.view.write: list data records do not share " *
            "a common time step count"))
    times=copy(view.list_time)
    length(times)<steps &&
        append!(times,Float64(index) for index in length(times):steps-1)
    length(times)>steps && resize!(times,steps)
    return times
end

# Raw TXT uses the model writer's numeric corner records. CSV separates the
# same fields with commas. A row is one element/time pair, with the step, time,
# list entity (0), element ordinal, XYZ and components repeated per corner.
function _native_list_view_text_payload(elements,times,separator::Char)
    io=IOBuffer()
    for step in eachindex(times),(ordinal,element) in enumerate(elements)
        first_value=true
        # Rows pair each node coordinate with its node value; value-only
        # scheme nodes beyond the stored geometry are not representable.
        for node in 1:min(size(element.coords,2),size(element.values,2))
            for value in (step-1,times[step],0,ordinal-1)
                first_value || print(io,separator)
                @printf(io,"%.17g",value);first_value=false
            end
            for axis in 1:3
                print(io,separator);@printf(io,"%.17g",element.coords[axis,node])
            end
            for component in axes(element.values,1)
                print(io,separator);@printf(io,"%.17g",element.values[component,node,step])
            end
        end
        separator==' ' && print(io,' ')
        println(io)
    end
    return take!(io)
end

function _native_list_view_vtk_values(io,elements,step)
    for element in elements,node in axes(element.values,2)
        for component in axes(element.values,1)
            component==1 || print(io,' ')
            @printf(io,"%.17g",element.values[component,node,step])
        end
        println(io)
    end
end

function _native_list_view_vtk_payload(elements,times,caller)
    components=isempty(elements) ? 1 : size(first(elements).values,1)
    all(element->size(element.values,1)==components,elements) || throw(ArgumentError(
        "$caller: VTK list export requires a common scalar, vector or tensor component count"))
    cell_types=(point=1,line=3,triangle=5,quadrangle=9,tetrahedron=10,
        hexahedron=12,prism=13,pyramid=14)
    all(element->haskey(cell_types,element.kind) && element.scheme===nothing &&
        isempty(element.suffix),elements) || throw(ArgumentError(
        "$caller: VTK list export supports only linear point, line, triangle, quadrangle, "*
        "tetrahedron, hexahedron, prism and pyramid cells"))
    points=foldl((count,element)->Base.checked_add(count,size(element.coords,2)),elements;init=0)
    cell_size=Base.checked_add(points,length(elements))
    io=IOBuffer()
    println(io,"# vtk DataFile Version 2.0\nTessella list view\nASCII\nDATASET UNSTRUCTURED_GRID")
    println(io,"POINTS ",points," double")
    # A list element owns its corner values. Coincident coordinates cannot
    # identify shared points: adjacent elements may carry discontinuous data.
    for element in elements,node in axes(element.coords,2)
        @printf(io,"%.17g %.17g %.17g\n",element.coords[1,node],element.coords[2,node],element.coords[3,node])
    end
    println(io,"CELLS ",length(elements),' ',cell_size)
    offset=0
    for element in elements
        count=size(element.coords,2);print(io,count)
        for node in 0:count-1;print(io,' ',offset+node);end
        println(io);offset+=count
    end
    println(io,"CELL_TYPES ",length(elements))
    for element in elements;println(io,getproperty(cell_types,element.kind));end
    println(io,"FIELD FieldData 1\nTIME 1 ",length(times)," double")
    for (step,time) in enumerate(times)
        step==1 || print(io,' ');@printf(io,"%.17g",time)
    end
    println(io)
    println(io,"POINT_DATA ",points)
    if components==9
        # Legacy TENSORS are symmetric. A FIELD array retains all nine public
        # tensor components, including nonsymmetric input, for every time step.
        println(io,"FIELD ListData ",length(times))
        for step in eachindex(times)
            println(io,"step_",step-1," 9 ",points," double")
            _native_list_view_vtk_values(io,elements,step)
        end
    else
        for step in eachindex(times)
            if components==1
                println(io,"SCALARS step_",step-1," double 1\nLOOKUP_TABLE default")
            else
                println(io,"VECTORS step_",step-1," double")
            end
            _native_list_view_vtk_values(io,elements,step)
        end
    end
    return take!(io)
end

# Record family → MSH catalog family for list-view mesh export.
const _LIST_MSH_FAMILY=(point=:pnt,line=:lin,triangle=:tri,quadrangle=:qua,
    tetrahedron=:tet,hexahedron=:hex,prism=:pri,pyramid=:pyr)

# (family, node count) → Gmsh 4.15.2 type tag. Complete shapes win; the
# incomplete (serendipity) table then fills gaps like the 13-node quadratic
# pyramid (type 19). Higher-order incomplete pyramids have no upstream
# `MElement` class — Gmsh 4.15.2 itself rejects type 125 on import — so they
# are excluded from the table and fail as unsupported rather than emitting
# an unreadable file.
const _LIST_MSH_TYPE_TABLE=let table=Dict{Tuple{Symbol,Int},Int}()
    for serendipity in (false,true)
        for (tag,spec) in Elements.MSH_CATALOG
            spec.serendipity==serendipity || continue
            serendipity && spec.family===:pyr && spec.order>2 && continue
            key=(spec.family,spec.nnodes)
            haskey(table,key) || (table[key]=tag)
        end
    end
    table
end

# Mesh-based MSH 2.2 export of a list view, matching upstream's
# PViewDataList::writeMSH path: first-seen node table sharing, per-record
# element connectivity in stored node order, an $InterpolationScheme section
# when the view binds family matrices, and one $ElementNodeData section per
# step (the section's component count comes from the first populated type in
# canonical order; rows emit the element's flat per-step value prefix).
function _native_list_view_msh_payload(view,elements,times,caller;
        name=view.name,scheme::Bool=true)
    view=_view_data_owner(view)
    interp=scheme ? view.interp : Dict{Int32,Vector{Matrix{Float64}}}()
    io=IOBuffer()
    println(io,"\$MeshFormat");println(io,"2.2 0 8");println(io,"\$EndMeshFormat")
    tags=Dict{NTuple{3,Float64},Int32}();nodes=NTuple{3,Float64}[]
    connectivity=Vector{Int32}[]
    for element in elements
        ids=Int32[]
        for column in axes(element.coords,2)
            key=(element.coords[1,column],element.coords[2,column],
                 element.coords[3,column])
            push!(ids,get!(tags,key) do
                push!(nodes,key);Int32(length(nodes))
            end)
        end
        push!(connectivity,ids)
    end
    println(io,"\$Nodes");println(io,length(nodes))
    for (tag,(x,y,z)) in enumerate(nodes)
        @printf(io,"%d %.16g %.16g %.16g\n",tag,x,y,z)
    end
    println(io,"\$EndNodes\n\$Elements");println(io,length(elements))
    for (tag,(element,ids)) in enumerate(zip(elements,connectivity))
        family=_LIST_MSH_FAMILY[element.kind]
        msh=get(_LIST_MSH_TYPE_TABLE,(family,size(element.coords,2)),0)
        msh==0 && throw(ArgumentError(
            "$caller: no supported MSH element type for a $(element.kind) " *
            "record with $(size(element.coords,2)) geometry nodes"))
        print(io,tag,' ',msh," 2 1 1")
        for node in ids;print(io,' ',node);end
        println(io)
    end
    println(io,"\$EndElements")
    if !isempty(interp)
        println(io,"\$InterpolationScheme\n\"INTERPOLATION_SCHEME\"")
        println(io,length(interp))
        for code in sort!(collect(keys(interp)))
            matrices=interp[code]
            length(matrices)>=2 || throw(ArgumentError(
                "$caller: malformed interpolation matrices for family $code"))
            # Upstream's MSH writer emits the two value matrices.
            println(io,code);println(io,2)
            for matrix in matrices[1:2]
                matrix isa AbstractMatrix && size(matrix,1)>=1 &&
                    size(matrix,2)>=1 || throw(ArgumentError(
                        "$caller: malformed interpolation matrices for family $code"))
                println(io,size(matrix,1),' ',size(matrix,2))
                for row in axes(matrix,1)
                    for column in axes(matrix,2)
                        @printf(io,"%.16g ",matrix[row,column])
                    end
                    println(io)
                end
            end
        end
        println(io,"\$EndInterpolationScheme")
    end
    # The shared section component count is the first populated canonical
    # type's; rows carry that many components per declared value node —
    # extra components of wider records keep upstream's flat-prefix
    # truncation, narrower records are zero-padded.
    components=0
    for type_key in view.list_types
        get(view.list_counts,type_key,0)>0 || continue
        layout=Post._pos_list_record_layout(view.interp,type_key)
        layout===nothing && throw(ArgumentError(
            "$caller: unsupported list data type $type_key"))
        components=layout.ncomp;break
    end
    for (index,time) in enumerate(times)
        println(io,"\$ElementNodeData")
        strings=String["\""*name*"\""]
        isempty(interp) || push!(strings,"\"INTERPOLATION_SCHEME\"")
        println(io,length(strings));foreach(line->println(io,line),strings)
        println(io,1);@printf(io,"%.16g\n",time)
        println(io,3);println(io,index-1);println(io,components)
        println(io,length(elements))
        for (tag,element) in enumerate(elements)
            nval=size(element.values,2)
            print(io,tag,' ',nval)
            width=nval*components
            flat=vec(@view element.values[:,:,index])
            for column in 1:width
                @printf(io," %.16g",column<=length(flat) ? flat[column] : 0.0)
            end
            println(io)
        end
        println(io,"\$EndElementNodeData")
    end
    return take!(io)
end

function _native_list_view_payload(view,extension,caller;name=view.name)
    if extension==".pos"
        parsed=_native_list_view_pos(view,caller;name)
        io=IOBuffer();Post.write_pos(io,parsed)
        return take!(io)
    end
    elements=_native_list_view_records(view,caller)
    times=_native_list_view_times(view,elements)
    extension==".txt" && return _native_list_view_text_payload(elements,times,' ')
    extension==".csv" && return _native_list_view_text_payload(elements,times,',')
    extension==".vtk" && return _native_list_view_vtk_payload(elements,times,caller)
    extension==".msh" || throw(ArgumentError("$caller: unsupported list export format"))
    scheme=!iszero(get(OPTIONS,"PostProcessing.SaveInterpolationMatrices",1.0))
    return _native_list_view_msh_payload(view,elements,times,caller;
        name,scheme)
end

include("APINativeModelViewFormats.jl")

function _view_write(tag,filename,append=false)
    caller="API.view.write"
    filename isa AbstractString || throw(ArgumentError("$caller: fileName must be a string"))
    append_value=_mesh_query_bool(append,caller,"append")
    return lock(STATE_LOCK) do
        _slot_capture_locked!()
        # Aliases export the shared data under their own per-view name.
        name=_session_view(tag,caller).name
        view=_view_data_owner(_session_view(tag,caller))
        extension=lowercase(splitext(filename)[2])
        extension in (".msh",".pos",".txt",".csv",".vtk") || throw(ArgumentError(
            "$caller: this view export format is not implemented"))
        if isempty(view.model_step_ranges)
            extension in (".msh",".pos",".txt",".csv",".vtk") || throw(ArgumentError(
                "$caller: list-based views support .pos, .msh, .txt, .csv and .vtk export"))
            extension==".vtk" && append_value && throw(ArgumentError(
                "$caller: appending list data to a VTK file is not supported"))
            extension==".msh" && append_value && throw(ArgumentError(
                "$caller: appending list data to a MSH file is not supported"))
            if extension==".pos" && !isempty(view.interp)
                message="Cannot export datasets with interpolation matrices in "*
                    "old POS format: consider using the new mesh-based format "*
                    "instead, or select 'Adapt post-processing data' before exporting"
                NATIVE_LOOKUP_QUERY_ERRORS[]+=1
                NATIVE_MSG_LAST_ERROR[]=message
                _native_diagnostic_emit_locked(:error,message)
                throw(ArgumentError(message))
            end
            if extension==".msh"
                msh_version=get(OPTIONS,"Mesh.MshFileVersion",4.1)
                msh_version!=2.2 && _native_diagnostic_emit_locked(:warning,
                    "Mesh-based export of list-based datasets not available "*
                    "with MSH $(msh_version): using MSH 2.2")
            end
            bytes=_native_list_view_payload(view,extension,caller;name)
            Sys.iswindows() && (bytes=Vector{UInt8}(codeunits(replace(String(bytes),"\n"=>"\r\n"))))
            _native_diagnostic_emit_locked(:info,"Writing '$filename'...")
            if !iszero(get(OPTIONS,"PostProcessing.Binary",0.0))
                message=extension==".pos" ? "Only parsed .pos files can be exported for this view type" :
                    "Only ASCII list data can be exported in $extension format"
                _native_diagnostic_emit_locked(:warning,message)
            end
            Base.open(filename,append_value ? "a" : "w") do stream;write(stream,bytes);end
            _native_diagnostic_emit_locked(:info,"Done writing '$filename'")
            return nothing
        end
        isempty(view.list_types) || throw(ArgumentError(
            "$caller: exporting a view containing both list and model data is not implemented"))
        model_msh=extension==".msh"
        version=get(OPTIONS,"Mesh.MshFileVersion",4.1)
        !model_msh || version in (2.2,4.1) || throw(ArgumentError("$caller: unsupported MSH version"))
        binary=!iszero(get(OPTIONS,"PostProcessing.Binary",0.0))
        save_mesh=!iszero(get(OPTIONS,"PostProcessing.SaveMesh",1.0))
        scheme=!iszero(get(OPTIONS,"PostProcessing.SaveInterpolationMatrices",1.0))
        groups=_native_view_writer_groups(view,filename)
        plans=[(group=group,mesh=_native_view_writer_mesh(group.slot,caller)) for group in groups]
        # Prevalidate all held rows, including mutable-vector revalidation,
        # before a target is opened. Live-object errors publish native prefixes.
        for plan in plans,key in plan.group.keys
            _native_view_writer_section(view,key;scheme=model_msh && scheme,name)
        end
        if extension in (".csv",".vtk")
            append_value && throw(ArgumentError("$caller: appending model data in $extension format is not supported"))
            publications=_native_model_view_publications(view,plans,filename,extension,caller)
            _native_diagnostic_emit_locked(:info,"Writing '$filename'...")
            length(publications)>1 && _native_diagnostic_emit_locked(:info,
                "Exporting multi-mesh view in separate files")
            binary && _native_diagnostic_emit_locked(:warning,
                "Only ASCII model data can be exported in $extension format")
            for publication in publications
                bytes=!Sys.iswindows() ? publication.bytes :
                    Vector{UInt8}(codeunits(replace(String(publication.bytes),"\n"=>"\r\n")))
                Base.open(publication.path,"w") do stream;write(stream,bytes);end
            end
            _native_diagnostic_emit_locked(:info,"Done writing '$filename'")
            return nothing
        end
        _native_diagnostic_emit_locked(:info,"Writing '$filename'...")
        if !model_msh
            result=_native_view_writer_text_payload(view,plans;extension,name)
            if result.error!==nothing
                NATIVE_LOOKUP_QUERY_ERRORS[]+=1;NATIVE_MSG_LAST_ERROR[]=result.error
                _native_diagnostic_emit_locked(:error,result.error)
                throw(ArgumentError(result.error))
            end
            if extension==".pos"
                interpolation=_native_view_interpolation_history(view.model_contexts)
                present=interpolation===nothing ?
                    any(plan->any(block->!isempty(block.nodes),plan.mesh.blocks),plans) : !isempty(interpolation)
                present && _native_diagnostic_emit_locked(:warning,
                    "Discarding interpolation matrices when saving in .pos format")
                binary && _native_diagnostic_emit_locked(:warning,
                    "Only parsed .pos files can be exported for this view type")
            end
            bytes=!Sys.iswindows() ? result.bytes :
                Vector{UInt8}(codeunits(replace(String(result.bytes),"\n"=>"\r\n")))
            Base.open(filename,append_value && extension==".pos" ? "a" : "w") do io;write(io,bytes);end
            _native_diagnostic_emit_locked(:info,"Done writing '$filename'")
            return nothing
        end
        length(plans)>1 && _native_diagnostic_emit_locked(:info,
            "Exporting multi-mesh view in separate files")
        for plan in plans
            lookup=_native_view_writer_lookup!(plan.group.slot,plan.mesh,view.model_kind)
            result=_native_view_writer_payload(view,plan.group,plan.mesh;
                version,binary,save_mesh,scheme,lookup,name)
            bytes=binary || !Sys.iswindows() ? result.bytes :
                Vector{UInt8}(codeunits(replace(String(result.bytes),"\n"=>"\r\n")))
            Base.open(plan.group.path,append_value ? "a" : "w") do io;write(io,bytes);end
            if result.error!==nothing
                NATIVE_LOOKUP_QUERY_ERRORS[]+=1;NATIVE_MSG_LAST_ERROR[]=result.error
                _native_diagnostic_emit_locked(:error,result.error)
                throw(ArgumentError(result.error))
            end
        end
        _native_diagnostic_emit_locked(:info,"Done writing '$filename'")
        return nothing
    end
end

# Model-view interpolation belongs to the first successful model-data storage
# creation, and is replaced when a type change recreates that storage.
function _native_view_interpolation_history(contexts)
    history=nothing
    for context in values(contexts)
        value=get(context,:interpolation_msh,nothing)
        value===nothing && continue
        value isa Tuple && all(msh->msh isa Integer && !(msh isa Bool) && msh>0,value) ||
            throw(ArgumentError("model view: malformed interpolation birth history"))
        foreach(msh_spec,value)
        history===nothing ? (history=value) : history==value || throw(ArgumentError(
            "model view: interpolation birth history differs across allocated steps"))
    end
    return history
end

function _native_view_interpolation_types(slot)
    slot===nothing && throw(ArgumentError("model view: interpolation birth has no current model"))
    mesh=_native_view_writer_mesh(slot,"model view interpolation birth")
    return Tuple(Int(block.msh) for block in mesh.blocks if !isempty(block.nodes))
end

# Rational pyramid space: u^i v^j w^k/(1-w)^min(i,j),
# max(i,j)+k <= p. See Bergot, Cohen and Durufle, Higher-order
# finite elements for hybrid meshes using new nodal pyramidal elements.
# The quadratic serendipity space omits the base interior bubble.
function _native_view_writer_pyramid_basis(msh::Int)
    spec=msh_spec(msh)
    spec.family===:pyr && spec.order>=2 || throw(ArgumentError(
        "model view writer: this rational pyramid order is not implemented"))
    # The quadratic serendipity element only omits the (p,p,0) base bubble,
    # which the lattice filter below removes; incomplete pyramids of order >=3
    # omit more and fail the node-count certificate instead.
    spec.serendipity && spec.order!=2 && throw(ArgumentError(
        "model view writer: this rational pyramid order is not implemented"))
    p=spec.order
    lattice=NTuple{3,Int}[(i,j,k) for i in 0:p for j in 0:p
        for k in 0:p-max(i,j) if !(spec.serendipity && i==p && j==p)]
    length(lattice)==spec.nnodes || throw(ArgumentError(
        "model view writer: rational pyramid space does not cover its nodes"))
    nodes=Elements.lagrange_nodes(msh);n=length(lattice)
    matrix=Matrix{Float64}(undef,n,n);exponents=Matrix{Int}(undef,n,3)
    for row in 1:n
        i,j,k=lattice[row];exponents[row,:].=(i,j,k)
        for column in 1:n
            u,v,w=nodes[:,column]
            matrix[row,column]=w==1.0 && max(i,j)>0 ? 0.0 :
                u^i*v^j*w^k/(1.0-w)^min(i,j)
        end
    end
    coefficients=matrix\Matrix{Float64}(I,n,n)
    return (family=6,coefficients,exponents)
end
