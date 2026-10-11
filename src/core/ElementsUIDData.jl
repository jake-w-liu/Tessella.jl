# Native model views can outlive their source cells. These rows use actual
# UIDs; the existing dense MshDataSection constructor remains separate.
struct MshUIDDataSection
    name::String
    strings::Vector{String}
    reals::Vector{Float64}
    header::Vector{Int64}
    tags::Vector{UInt64}
    row_nodes::Vector{Int32}
    nodes::Vector{UInt64}
    implicit_nodes::Bool
    values::Vector{Float64}
    anchor::Int32
    logical_extent::UInt64
    function MshUIDDataSection(name,strings,reals,header,tags,row_nodes,nodes,
                               implicit_nodes,values,anchor;logical_extent)
        name in ("NodeData","ElementData","ElementNodeData") || throw(ArgumentError(
            "MSH UID data: unsupported section kind"))
        strings isa AbstractVector && all(value->value isa AbstractString,strings) ||
            throw(ArgumentError("MSH UID data: strings must be a vector of strings"))
        all(value->!occursin('\n',value) && !occursin('\r',value),strings) || throw(ArgumentError(
            "MSH UID data: string header must remain on one line"))
        all(value->value isa Integer && !(value isa Bool),header) || throw(ArgumentError(
            "MSH UID data: integer header must contain integers"))
        h=Int64.(header)
        length(h)>=3 && h[1]>=0 && h[2]>0 && h[3]==length(tags) || throw(ArgumentError(
            "MSH UID data: header does not align with its rows"))
        all(value->value isa Integer && !(value isa Bool) && value>=0,tags) && allunique(tags) ||
            throw(ArgumentError("MSH UID data: row UIDs must be nonnegative and distinct"))
        all(value->value isa Integer && !(value isa Bool) && 0<value<=typemax(Int32),row_nodes) ||
            throw(ArgumentError("MSH UID data: invalid per-row node count"))
        count=length(tags)
        length(row_nodes)==count || throw(ArgumentError("MSH UID data: node counts do not align"))
        implicit_nodes isa Bool || throw(ArgumentError("MSH UID data: dialect flag must be Bool"))
        if name=="ElementNodeData"
            expected_nodes=implicit_nodes ? 0 : sum(Int,row_nodes;init=0)
            length(nodes)==expected_nodes || throw(ArgumentError("MSH UID data: explicit node UIDs do not align"))
        else
            !implicit_nodes && isempty(nodes) && all(==(1),row_nodes) || throw(ArgumentError(
                "MSH UID data: non-nodal rows carry unexpected connectivity"))
        end
        all(value->value isa Integer && !(value isa Bool) && value>=0,nodes) || throw(ArgumentError(
            "MSH UID data: explicit node UIDs must be nonnegative"))
        expected=Base.checked_mul(sum(Int,row_nodes;init=0),Int(h[2]))
        length(values)==expected || throw(ArgumentError("MSH UID data: values do not align with complete rows"))
        all(value->value isa Real && !(value isa Bool) && isfinite(value),values) &&
            all(value->value isa Real && !(value isa Bool) && isfinite(value),reals) || throw(ArgumentError(
                "MSH UID data: values and real header must be finite numbers"))
        logical_extent isa Integer && !(logical_extent isa Bool) && logical_extent>=0 ||
            throw(ArgumentError("MSH UID data: invalid logical storage extent"))
        extent=UInt64(logical_extent)
        all(tag->UInt64(tag)<extent,tags) || throw(ArgumentError(
            "MSH UID data: logical storage does not cover all row UIDs"))
        anchor isa Integer && !(anchor isa Bool) && 0<=anchor<=5 || throw(ArgumentError(
            "MSH UID data: invalid section anchor"))
        return new(String(name),String.(strings),Float64.(reals),h,UInt64.(tags),Int32.(row_nodes),
                   UInt64.(nodes),implicit_nodes,Float64.(values),Int32(anchor),extent)
    end
end

function _msh_uid_validate(section::MshUIDDataSection)
    MshUIDDataSection(section.name,section.strings,section.reals,section.header,
        section.tags,section.row_nodes,section.nodes,section.implicit_nodes,section.values,
        section.anchor;logical_extent=section.logical_extent)
    return nothing
end

function _msh_uid_extent(tags,initial::Integer;previous::UInt64=UInt64(0))
    initial>=0 || throw(ArgumentError("MSH UID data: negative initial storage extent"))
    extent=max(UInt64(initial),previous)
    for tag in tags
        UInt64(tag)>=extent && (extent=Base.checked_add(UInt64(tag),UInt64(100)))
    end
    return extent
end

function _msh_uid_from_dense(section::MshDataSection;node_tags,element_tags,
                            initial::Integer=0,logical_extent=nothing)
    source=section.name=="NodeData" ? node_tags : element_tags
    indices=section.name=="NodeData" ? section.nodes : section.elements
    all(index->0<index<=length(source),indices) || throw(ArgumentError(
        "MSH UID data: dense row has no actual source identity"))
    tags=UInt64[source[index] for index in indices]
    counts=section.name=="ElementNodeData" ? copy(section.row_nodes) : fill(Int32(1),length(tags))
    nodes=if section.name=="ElementNodeData" && !section.implicit_nodes
        all(index->0<index<=length(node_tags),section.nodes) || throw(ArgumentError(
            "MSH UID data: explicit node has no actual source identity"))
        UInt64[node_tags[index] for index in section.nodes]
    else;UInt64[];end
    extent=logical_extent===nothing ? _msh_uid_extent(tags,initial) : logical_extent
    return MshUIDDataSection(section.name,section.strings,section.reals,section.header,
        tags,counts,nodes,section.implicit_nodes,section.values,section.anchor;logical_extent=extent)
end

function _msh_uid_numbered(section::MshUIDDataSection,mapping;renumber::Bool=true)
    _msh_uid_validate(section)
    copied()=MshUIDDataSection(section.name,section.strings,section.reals,section.header,
        section.tags,section.row_nodes,section.nodes,section.implicit_nodes,section.values,
        section.anchor;logical_extent=section.logical_extent)
    renumber || return (section=copied(),warning=nothing)
    section.name=="ElementNodeData" && !section.implicit_nodes && throw(ArgumentError(
        "MSH UID data: explicit-node numbering needs a separately bound node map"))
    outside=UInt64[uid for uid in keys(mapping) if UInt64(uid)>=section.logical_extent]
    if !isempty(outside)
        return (section=copied(),warning="Wrong source index $(minimum(outside)) in step data renumbering")
    end
    kept=Int[row for row in eachindex(section.tags) if haskey(mapping,section.tags[row])]
    sort!(kept;by=row->mapping[section.tags[row]],alg=MergeSort)
    offsets=Int[1]
    for count in section.row_nodes
        push!(offsets,Base.checked_add(offsets[end],Base.checked_mul(Int(count),Int(section.header[2]))))
    end
    tags=UInt64[mapping[section.tags[row]] for row in kept]
    values=Float64[value for row in kept for value in @view section.values[offsets[row]:offsets[row+1]-1]]
    header=copy(section.header);header[3]=length(tags)
    extent=isempty(mapping) ? UInt64(1) : Base.checked_add(maximum(UInt64(value) for value in Base.values(mapping)),UInt64(1))
    result=MshUIDDataSection(section.name,section.strings,section.reals,header,tags,
        section.row_nodes[kept],UInt64[],section.implicit_nodes,values,section.anchor;logical_extent=extent)
    return (section=result,warning=nothing)
end

# File insertion updates an existing row's prefix. A one-node nodal write
# keeps its existing multiplicity; wider writes set their own multiplicity.
# Owned arrays retain initialized tails without native out-of-bounds storage.
function _msh_uid_file_rows(name,tags,counts,nodes,implicit,values,components)
    allunique(tags) && return (;tags,counts,nodes,values)
    name=="ElementNodeData" && !implicit && throw(ArgumentError(
        "MSH UID data: duplicate explicit-node rows need a bound insertion dialect"))
    selected=Dict{UInt64,Tuple{Int32,Vector{Float64}}}()
    offset=1
    for row in eachindex(tags)
        count=Int(counts[row]);width=Base.checked_mul(count,components)
        incoming=copy(@view values[offset:offset+width-1])
        previous=get(selected,tags[row],nothing)
        if name=="ElementNodeData" && count==1 && previous!==nothing && previous[1]>1
            retained=copy(previous[2]);retained[1:width]=incoming
            selected[tags[row]]=(previous[1],retained)
        else
            selected[tags[row]]=(counts[row],incoming)
        end
        offset+=width
    end
    ordered=sort!(collect(keys(selected)))
    return (tags=ordered,counts=Int32[selected[tag][1] for tag in ordered],nodes=UInt64[],
            values=Float64[value for tag in ordered for value in selected[tag][2]])
end

# MSH has no storage-capacity field. A held logical extent can accompany an
# in-memory roundtrip; a new import grows fresh storage in file row order.
function _msh_uid_decode(section::MshAncillarySection;initial::Integer=0,
                         logical_extent=nothing,swap::Bool=false)
    section.name in ("NodeData","ElementData","ElementNodeData") || throw(ArgumentError(
        "MSH UID data: ancillary section is not model data"))
    io=IOBuffer(section.payload)
    integer_line()=parse(Int,strip(readline(io)))
    bound=length(section.payload)
    nstrings=integer_line();0<=nstrings<=bound || throw(ArgumentError("MSH UID data: invalid string count"))
    strings=String[readline(io) for _ in 1:nstrings]
    nreals=integer_line();0<=nreals<=bound || throw(ArgumentError("MSH UID data: invalid real count"))
    reals=Float64[parse(Float64,strip(readline(io))) for _ in 1:nreals]
    nheader=integer_line();3<=nheader<=bound || throw(ArgumentError("MSH UID data: invalid integer header count"))
    header=Int64[integer_line() for _ in 1:nheader]
    components=Int(header[2]);entries=Int(header[3])
    components>0 && 0<=entries<=bound || throw(ArgumentError("MSH UID data: invalid row header"))
    tags=UInt64[];counts=Int32[];nodes=UInt64[];values=Float64[]
    implicit=false;dialect=nothing
    for _ in 1:entries
        if section.binary
            raw=read(io,UInt32);swap && (raw=bswap(raw))
            tag=Int(reinterpret(Int32,raw));tag>=0 || throw(ArgumentError("MSH UID data: invalid binary UID"))
            push!(tags,UInt64(tag))
            count=if section.name=="ElementNodeData"
                raw=read(io,UInt32);swap && (raw=bswap(raw));Int(reinterpret(Int32,raw))
            else;1;end
            0<count<=typemax(Int32) || throw(ArgumentError("MSH UID data: invalid binary node count"))
            push!(counts,Int32(count));implicit=section.name=="ElementNodeData"
            for _ in 1:Base.checked_mul(count,components)
                raw=read(io,UInt64);swap && (raw=bswap(raw));push!(values,reinterpret(Float64,raw))
            end
        else
            tokens=split(strip(readline(io)))
            isempty(tokens) && throw(ArgumentError("MSH UID data: missing row"))
            tag=parse(UInt64,tokens[1]);push!(tags,tag)
            count=section.name=="ElementNodeData" ? parse(Int,tokens[2]) : 1
            0<count<=typemax(Int32) || throw(ArgumentError("MSH UID data: invalid row node count"))
            push!(counts,Int32(count))
            if section.name=="ElementNodeData"
                remaining=length(tokens)-2
                if remaining==Base.checked_mul(count,components)
                    row_dialect=true;first_value=3
                elseif remaining==Base.checked_mul(count,components+1)
                    row_dialect=false;first_value=count+3
                    append!(nodes,parse.(UInt64,tokens[3:count+2]))
                else;throw(ArgumentError("MSH UID data: nodal row width does not align"));end
                dialect===nothing || dialect==row_dialect || throw(ArgumentError("MSH UID data: mixed row dialects"))
                dialect=row_dialect;implicit=row_dialect
            else
                length(tokens)==components+1 || throw(ArgumentError("MSH UID data: row width does not align"))
                first_value=2
            end
            append!(values,parse.(Float64,tokens[first_value:end]))
        end
    end
    all(byte->byte in (0x09,0x0a,0x0d,0x20),read(io)) || throw(ArgumentError("MSH UID data: trailing row payload"))
    extent=logical_extent===nothing ? (section.uid_context===nothing ?
        _msh_uid_extent(tags,initial) : section.uid_context.logical_extent) : logical_extent
    rows=_msh_uid_file_rows(section.name,tags,counts,nodes,implicit,values,components)
    header[3]=length(rows.tags)
    return MshUIDDataSection(section.name,strings,reals,header,rows.tags,rows.counts,rows.nodes,
        implicit,rows.values,section.anchor;logical_extent=extent)
end

function _msh_uid_ancillary(section::MshUIDDataSection;binary::Bool=false)
    _msh_uid_validate(section)
    binary && section.name=="ElementNodeData" && !section.implicit_nodes && throw(ArgumentError(
        "MSH UID data: explicit-node dialect cannot be emitted as native binary data"))
    io=IOBuffer()
    println(io,length(section.strings));foreach(value->println(io,value),section.strings)
    println(io,length(section.reals));foreach(value->@printf(io,"%.17g\n",value),section.reals)
    println(io,length(section.header));foreach(value->println(io,value),section.header)
    offset=1;node_offset=1;components=Int(section.header[2])
    for row in eachindex(section.tags)
        tag=section.tags[row];count=Int(section.row_nodes[row]);width=Base.checked_mul(count,components)
        if binary
            tag<=typemax(Int32) || throw(ArgumentError("MSH UID data: native binary UID exceeds Int32"))
            write(io,Int32(tag))
            section.name=="ElementNodeData" && write(io,Int32(count))
            write(io,@view(section.values[offset:offset+width-1]))
        else
            print(io,tag)
            if section.name=="ElementNodeData"
                print(io,' ',count)
                if !section.implicit_nodes
                    for tag in @view section.nodes[node_offset:node_offset+count-1];print(io,' ',tag);end
                    node_offset+=count
                end
            end
            for value in @view section.values[offset:offset+width-1];@printf(io," %.17g",value);end
            println(io)
        end
        offset+=width
    end
    # Ancillary binary payloads exclude the separator before their terminator.
    # The mixed writer supplies that single newline after the encoded rows.
    return MshAncillarySection(section.name,binary,take!(io),section.anchor,
                               _MshUIDContext(section.logical_extent))
end

function _msh_uid_recombine_payload(mesh::MixedMesh;node_tags,element_tags)
    ancillary=MshAncillarySection[]
    for section in mesh.ancillary_sections
        if section.uid_context!==nothing || !(section.name in ("NodeData","ElementData","ElementNodeData"))
            push!(ancillary,section);continue
        end
        initial=section.name=="NodeData" ? length(node_tags) : length(element_tags)
        uid=_msh_uid_decode(section;initial)
        push!(ancillary,_msh_uid_ancillary(uid;binary=section.binary))
    end
    for section in mesh.data_sections
        initial=section.name=="NodeData" ? length(node_tags) : length(element_tags)
        uid=_msh_uid_from_dense(section;node_tags,element_tags,initial)
        push!(ancillary,_msh_uid_ancillary(uid))
    end
    return ancillary
end

function _msh_uid_numbered_payload(mesh::MixedMesh,node_mapping,element_mapping;kind=nothing)
    kind in (nothing,:node,:element) || throw(ArgumentError(
        "MSH UID data: numbering kind must be node or element"))
    ancillary=MshAncillarySection[];events=Tuple{Symbol,String}[]
    grouped=Dict{Tuple{String,String},Vector{Int}}()
    decoded=Dict{Int,MshUIDDataSection}()
    for (index,section) in enumerate(mesh.ancillary_sections)
        section.uid_context===nothing && continue
        kind===nothing || (section.name=="NodeData")==(kind===:node) || continue
        uid=_msh_uid_decode(section);decoded[index]=uid
        name=isempty(uid.strings) ? "" : first(uid.strings)
        push!(get!(()->Int[],grouped,(name,uid.name)),index)
    end
    for key in sort!(collect(keys(grouped)))
        indices=sort!(grouped[key];by=index->decoded[index].header[1],alg=MergeSort)
        steps=1+maximum(decoded[index].header[1] for index in indices)
        description=key[2]=="NodeData" ? "nodal" : "element"
        push!(events,(:info,"Renumbering $description model data ($steps $(steps==1 ? "step" : "steps"))"))
        for index in indices
            section=mesh.ancillary_sections[index];uid=decoded[index]
            mapping=section.name=="NodeData" ? node_mapping : element_mapping
            result=_msh_uid_numbered(uid,mapping)
            result.warning===nothing || push!(events,(:warning,result.warning))
        end
    end
    for (index,section) in enumerate(mesh.ancillary_sections)
        if !haskey(decoded,index)
            push!(ancillary,section);continue
        end
        uid=_msh_uid_decode(section)
        mapping=section.name=="NodeData" ? node_mapping : element_mapping
        result=_msh_uid_numbered(uid,mapping)
        push!(ancillary,_msh_uid_ancillary(result.section;binary=section.binary))
    end
    return (;ancillary,events)
end


# Model-data UIDs refer to current mesh objects. MSH2.2 writes dense node/cell
# labels, so bind the live rows through exactly the geometry writer's maps.
# A consumed cell table has no live map and retains its orphan rows. The mesh
# and its owned logical storage context remain unchanged by serialization.
function _msh_uid_writer_section(section::MshAncillarySection,writer)
    uid=_msh_uid_decode(section)
    data=writer.mesh.entity_data
    if writer.version==2.2 && data!==nothing
        node=uid.name=="NodeData"
        source=node ? data.external_node_tags :
            reduce(vcat,data.external_element_tags;init=UInt64[])
        target=node ? writer.node_tags : writer.element_tags
        length(source)==length(target) || throw(ArgumentError(
            "MSH UID data: writer source identities do not align"))
        if !isempty(source)
            mapping=Dict{UInt64,UInt64}(old=>new for (old,new) in zip(source,target))
            uid=_msh_uid_numbered(uid,mapping).section
        end
    end
    return _msh_uid_ancillary(uid;binary=writer.binary)
end
