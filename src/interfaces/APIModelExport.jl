# Current-model export consumes the actual cache/record union. It never creates
# a post-processing view or changes the live model, cache, tags, or filename.
import ..IO

function _model_export_source_tags(source,class,caller)
    source===nothing && return nothing
    if source isa MixedMesh && source.entity_data!==nothing
        data=source.entity_data
        return (nodes=data.external_node_tags,cells=reduce(vcat,data.external_element_tags;init=UInt64[]))
    end
    class!==nothing && class.public_tags!==nothing || throw(ArgumentError(
        "$caller: retained mesh metadata requires actual public UID authority"))
    return (nodes=class.public_tags.node_tags,cells=class.public_tags.element_tags)
end

function _model_export_data_sections(source,node_map,cell_map,caller)
    sections=Elements.MshDataSection[]
    source isa MixedMesh || return sections
    for section in source.data_sections
        all(index->node_map[index]>0,section.nodes) &&
        all(index->cell_map[index]>0,section.elements) || throw(ArgumentError(
            "$caller: physical selection would omit retained view-data records; use Mesh.SaveAll=1"))
        push!(sections,Elements.MshDataSection(section.name,section.strings,
            section.reals,section.header,cell_map[section.elements],
            node_map[section.nodes],section.row_nodes,section.implicit_nodes,
            section.values,section.anchor))
    end
    return sections
end

function _model_export_periodic(source,node_map,caller)
    links=Elements.MixedPeriodicLink[]
    source isa MixedMesh || return links
    for link in source.periodic_links
        all(index->node_map[index]>0,link.slave_nodes) &&
        all(index->node_map[index]>0,link.master_nodes) || throw(ArgumentError(
            "$caller: physical selection would omit periodic nodes; use Mesh.SaveAll=1"))
        push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
            node_map[link.slave_nodes],node_map[link.master_nodes];affine=link.affine))
    end
    return links
end

function _model_export_mesh(model,source,class,version,save_all,caller)
    source,class=_dim01_actual_cache(source,class;physical_names=model.physical_names)
    if source isa MixedMesh && !isempty(model.meshing.periodic_uid_history)
        source,class=_periodic_uid_history_source(model,source,class,caller)
    end
    complete=_native_view_writer_mesh((model=model,mesh=source,class=class),caller)
    complete.partition_data===nothing || version==4.1 || throw(ArgumentError(
        "$caller: retained partition and ghost metadata requires MSH 4.1"))
    data=complete.entity_data
    declarations=complete.partition_data===nothing ? data.entities :
        merge(data.entities,complete.partition_data.entities)
    physical=Set(key for (key,entity) in declarations if !isempty(entity.physical_tags))
    all_cells=save_all || isempty(physical)
    selections=[findall(column->all_cells ||
        (msh_dimension(block.msh),Int(data.block_entities[bi][column])) in physical,
        axes(block.nodes,2)) for (bi,block) in enumerate(complete.blocks)]
    referenced=falses(nnodes(complete))
    for (block,columns) in zip(complete.blocks,selections)
        referenced[block.nodes[:,columns]].=true
    end
    rows=version==2.2 ? findall(referenced) : findall(row->all_cells || referenced[row] ||
        (data.node_entities[row][1],Int(data.node_entities[row][2])) in physical,
        axes(complete.coords,2))
    inverse=zeros(Int32,nnodes(complete));inverse[rows]=Int32.(eachindex(rows))
    blocks=ElementBlock[];owners=Vector{Int32}[];uids=Vector{UInt64}[]
    complete_cell_uids=reduce(vcat,data.external_element_tags;init=UInt64[])
    complete_cell_map=zeros(Int32,length(complete_cell_uids));offset=0;emitted=0
    for (bi,block) in enumerate(complete.blocks)
        columns=selections[bi]
        for column in columns;emitted+=1;complete_cell_map[offset+column]=Int32(emitted);end
        offset+=size(block.nodes,2)
        isempty(columns) && continue
        weights=version==2.2 && save_all ? zeros(Int32,length(columns)) : block.tags[columns]
        push!(blocks,ElementBlock(block.msh,inverse[block.nodes[:,columns]],weights))
        push!(owners,data.block_entities[bi][columns]);push!(uids,data.external_element_tags[bi][columns])
    end
    # Native MSH2 SaveAll clears element physical associations. The detached
    # entity projection must agree with those cells while names remain present.
    export_entities=version==2.2 && save_all ? Dict(key=>Elements._native_mixed_entity(
        entity.dim,entity.tag,entity.bbox;boundaries=entity.boundaries,
        embedded_curves=entity.embedded_curves,partitions=entity.partitions,parent=entity.parent)
        for (key,entity) in data.entities) : data.entities
    metadata=Elements._native_mixed_entity_data(export_entities;
        node_entities=data.node_entities[rows],node_parametric=data.node_parametric[rows],
        external_node_tags=data.external_node_tags[rows],block_entities=owners,external_element_tags=uids)
    partition=complete.partition_data
    if partition!==nothing
        ghosts=Elements.MixedGhostElement[Elements.MixedGhostElement(
            complete_cell_map[record.element],record.partition,record.ghost_partitions)
            for record in partition.ghost_elements if complete_cell_map[record.element]>0]
        partition=Elements.MixedPartitionData(partition.num_partitions;entities=partition.entities,
            ghost_entities=partition.ghost_entities,ghost_elements=ghosts)
    end
    source_tags=_model_export_source_tags(source,class,caller)
    node_map=Int32[];cell_map=Int32[]
    if source_tags!==nothing
        positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags[rows]))
        cell_positions=Dict(uid=>Int32(index) for (index,uid) in enumerate(reduce(vcat,uids;init=UInt64[])))
        node_map=Int32[get(positions,uid,0) for uid in source_tags.nodes]
        cell_map=Int32[get(cell_positions,uid,0) for uid in source_tags.cells]
    end
    links=_model_export_periodic(source,node_map,caller)
    declared=Set((Int(link.dim),Int(link.slave_entity)) for link in links)
    for ((dim,slave),_) in model.periodic
        dim==3 && continue
        (dim,slave) in declared || throw(ArgumentError(
            "$caller: declared periodic entity lacks an actual retained mesh link"))
    end
    sections=_model_export_data_sections(source,node_map,cell_map,caller)
    ancillary=source isa MixedMesh ? source.ancillary_sections : Elements.MshAncillarySection[]
    return MixedMesh(complete.coords[:,rows],blocks;entity_data=metadata,
        elementary_entities=owners,physical_names=complete.physical_names,
        periodic_links=links,partition_data=partition,ancillary_sections=ancillary,
        data_sections=sections)
end

# Versioned MSH filename aliases. `.msh` stays `Mesh.MshFileVersion`-driven;
# `.msh2`/`.msh22` force MSH 2.2 and `.msh4`/`.msh41` force MSH 4.1. MSH 1.x
# and 3.x never shipped writable support and MSH 4.0 is not a stable public
# dialect, so those aliases reject loudly rather than guessing a version.
const _MSH_EXTENSION_VERSIONS=Dict{String,Union{Nothing,Float64}}(
    ".msh"=>nothing,".msh2"=>2.2,".msh22"=>2.2,".msh4"=>4.1,".msh41"=>4.1)

function _model_export_write(filename)
    caller="API.write"
    filename isa AbstractString && !isempty(filename) && !occursin('\0',filename) ||
        throw(ArgumentError("$caller: filename must be a nonempty NUL-free string"))
    extension=lowercase(splitext(filename)[2])
    if extension in (".msh1",".msh3",".msh40")
        throw(ArgumentError(
            "$caller: unsupported MSH extension `$extension` (supported " *
            "aliases: .msh, .msh2/.msh22 for MSH 2.2, .msh4/.msh41 for MSH 4.1)"))
    end
    (haskey(_MSH_EXTENSION_VERSIONS,extension) ||
        extension in (".unv",".mesh",".stl",".su2",".bdf",".nas",".inp",".neu",
                      ".p3d",".mail",".ir3")) || throw(ArgumentError(
        "$caller: current-model export supports .msh, .msh2/.msh22, " *
        ".msh4/.msh41, .unv, .mesh, .stl, .su2, .bdf, .nas, .inp, .neu, " *
        ".p3d, .mail and .ir3"))
    return lock(STATE_LOCK) do
        model=_model_locked()
        save_all=!iszero(get(OPTIONS,"Mesh.SaveAll",0.0))
        if haskey(_MSH_EXTENSION_VERSIONS,extension)
            version=something(_MSH_EXTENSION_VERSIONS[extension],
                get(OPTIONS,"Mesh.MshFileVersion",4.1))
            version in (2.2,4.1) || throw(ArgumentError(
                "$caller: MSH version must be 2.2 or 4.1"))
            target=_model_export_mesh(model,LAST_MESH[],LAST_MESH_CLASS[],version,
                save_all,caller)
            Elements.write_mixed_msh(filename,target;version,
                binary=!iszero(get(OPTIONS,"Mesh.Binary",0.0)))
        elseif extension==".stl"
            # STL carries only surface facets: no physical selection applies,
            # so the export projection always uses SaveAll semantics.
            target=_model_export_mesh(model,LAST_MESH[],LAST_MESH_CLASS[],4.1,
                true,caller)
            IO.write_stl(filename,target)
        else
            # The foreign writers use the complete MSH4-style projection:
            # actual UIDs, per-cell physical tags, and block entity ownership.
            # UNV/MEDIT carry no group metadata so names are dropped (upstream
            # behaviour); SU2 maps physical groups to MARKER sections with
            # names; BDF carries the physical tag in each card's PID field;
            # INP emits *ELSET/*NSET sections and NEU ELEMENT GROUP/BOUNDARY
            # CONDITIONS sections carrying tags plus physical names. P3D
            # exports only structured transfinite grids upstream, which the
            # projection cannot represent, so it always rejects; MAIL emits
            # the face triangles with trailing ` 0 0 0` edge placeholders;
            # IR3 emits self-describing vertex-count element records whose
            # tag field follows the Mesh.SaveElementTagType semantics.
            target=_model_export_mesh(model,LAST_MESH[],LAST_MESH_CLASS[],4.1,
                save_all,caller)
            # `Mesh.SaveElementTagType` is not a session-tracked option, so
            # the upstream default (1: elementary tag) applies; read it
            # anyway in case the table grows the key. `trunc` reproduces
            # the `(int)` cast for a raw table write.
            tag_type=trunc(Int,get(OPTIONS,"Mesh.SaveElementTagType",1.0))
            extension==".unv" ? IO.write_unv(filename,target) :
                extension==".mesh" ? IO.write_mesh(filename,target) :
                extension==".su2" ? IO.write_su2(filename,target) :
                extension==".inp" ? IO.write_inp(filename,target) :
                extension==".neu" ? IO.write_neu(filename,target) :
                extension==".p3d" ? IO.write_p3d(filename,target) :
                extension==".mail" ? IO.write_mail(filename,target) :
                extension==".ir3" ? IO.write_ir3(filename,target;
                    element_tag_type=tag_type) :
                    IO.write_bdf(filename,target)
        end
        return nothing
    end
end

"""
    write(filename)

Write the current model's actual mesh to a `.msh`, `.msh2`, `.msh22`, `.msh4`,
`.msh41`, `.stl`, `.unv`, `.mesh`, `.su2`, `.bdf`, `.nas`, `.inp`, `.neu`,
`.p3d`, `.mail` or `.ir3` file. `.msh` uses `Mesh.MshFileVersion`
to select 2.2 or 4.1 while the versioned aliases force the corresponding
dialect; `.msh1`, `.msh3` and `.msh40` reject. `Mesh.Binary` selects the
encoding, and `Mesh.SaveAll`
controls physical selection for MSH and the foreign formats. `.stl` writes
ASCII surface facets (quads split on the first-third vertex diagonal) and
rejects volume elements; MSH4 preserves actual UIDs and retained periodic,
partition, ghost, and data metadata; MSH2 uses its native dense numbering.
`.unv` emits I-DEAS datasets 2411/2412 with elementary and physical tags;
`.mesh` emits INRIA/Medit sections with `ref` carrying the physical tag (or
entity tag when no physical tag is set); `.su2` emits NDIME/NPOIN/NELEM with
physical groups as named MARKER sections; `.bdf`/`.nas` emit Nastran
GRID and element cards with the physical tag in the property-id field; `.inp`
emits Abaqus `*NODE`/`*ELEMENT` blocks plus `*ELSET`/`*NSET` physical
sections (entity sets use the `Volume|Surface|Line|Point<tag>` convention);
`.neu` emits a Gambit neutral deck with `ELEMENT GROUP` sections per top-dim
physical tag and `BOUNDARY CONDITIONS` sections resolved from dim-1 physical
cells; `.mail` emits the CEA triangulation deck (` nv nt` header,
`%19.10E` vertex rows, dense index triples and ` 0 0 0` edge placeholders)
keeping only triangle cells; `.ir3` emits the Iridium deck (`33` format
line, `numVertices num2D num3D` header, `%.16g` vertex rows and
self-describing `num tag numVert` records with per-dimension numbering)
whose tag field carries the elementary entity tag (the upstream
`Mesh.SaveElementTagType` default; `IO.write_ir3` also honors values 1–3);
`.p3d` writes only
entities with structured transfinite grids upstream, so every export
rejects with "no structured grids to save". The foreign formats retain
neither partitions, periodic links nor data
sections; `.su2`, `.inp` and `.neu` retain physical names.
Unsupported metadata/format combinations reject before replacing the target.
The model, mesh cache, views, counters, and model filename remain unchanged.
"""
write(filename)=_model_export_write(filename)
