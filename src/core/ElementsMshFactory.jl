# Private MSH4.1 constructor replay. Generic readers keep literal positive tags.
struct _MixedMshFactoryContext
    node_max::UInt64
    element_max::UInt64
end
function _mixed_msh_factory_context(node_max,element_max)
    function counter(value,name)
        value isa Integer && !(value isa Bool) && 0<=value<=typemax(UInt64) ||
            throw(ArgumentError("read_mixed_msh: invalid native $name counter"))
        return UInt64(value)
    end
    return _MixedMshFactoryContext(counter(node_max,"node"),counter(element_max,"element"))
end
mutable struct _MixedMshFactoryState
    initial::_MixedMshFactoryContext
    node_max::UInt64
    element_max::UInt64
    final_nodes::Set{UInt64}
    final_cells::Set{UInt64}
    node_events::Vector{Tuple{Int,Int32,UInt64,UInt64}}
    cell_events::Vector{Tuple{Int,Int32,UInt64,UInt64}}
end
_mixed_msh_factory_state(::Nothing)=nothing
_mixed_msh_factory_state(input::_MixedMshFactoryContext)=_MixedMshFactoryState(
    input,input.node_max,input.element_max,Set{UInt64}(),Set{UInt64}(),
    Tuple{Int,Int32,UInt64,UInt64}[],Tuple{Int,Int32,UInt64,UInt64}[])
function _mixed_msh_factory_tag!(state,kind::Symbol,owner,raw::UInt64)
    state===nothing && return raw
    node=kind===:node
    kind in (:node,:cell) || throw(ArgumentError("read_mixed_msh: invalid factory kind"))
    maximum_tag=node ? state.node_max : state.element_max
    final=if iszero(raw)
        maximum_tag<typemax(UInt64) || throw(ArgumentError(
            "read_mixed_msh: native $kind constructor counter overflows UInt64"))
        maximum_tag+UInt64(1)
    else
        raw
    end
    seen=node ? state.final_nodes : state.final_cells
    final in seen && throw(ArgumentError(
        "read_mixed_msh: native $kind constructors collide at tag $final"))
    push!(seen,final)
    if node
        state.node_max=max(maximum_tag,final)
        push!(state.node_events,(Int(owner[1]),Int32(owner[2]),raw,final))
    else
        state.element_max=max(maximum_tag,final)
        push!(state.cell_events,(Int(owner[1]),Int32(owner[2]),raw,final))
    end
    return final
end
function _mixed_msh_factory_bindings(mesh)
    data=mesh.entity_data
    data!==nothing || throw(ArgumentError("read_mixed_msh: factory input lacks entity metadata"))
    nodes=Tuple((data.node_entities[row],data.external_node_tags[row],
        Tuple(mesh.coords[:,row]),data.node_parametric[row]===nothing ? nothing :
        Tuple(data.node_parametric[row])) for row in axes(mesh.coords,2))
    cells=Tuple((Int(block.msh),owners[column],labels[column],block.tags[column],
        Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (block,owners,labels) in zip(mesh.blocks,data.block_entities,
            data.external_element_tags) for column in axes(block.nodes,2))
    return (nodes=nodes,cells=cells)
end
function _read_mixed_msh_factory_with_identities(path,input::_MixedMshFactoryContext;kwargs...)
    source_sha=bytes2hex(sha256(read(path)))
    reference=Ref{Any}(nothing)
    mesh=read_mixed_msh(path;kwargs...,_source_identities=reference,_factory_context=input)
    source_sha==bytes2hex(sha256(read(path))) || throw(ArgumentError(
        "read_mixed_msh: factory input changed while being read"))
    identities=reference[]
    identities!==nothing && identities.factory_events!==nothing || throw(ArgumentError(
        "read_mixed_msh: missing native constructor event receipt"))
    factory=merge(identities.factory_events,(source_sha256=source_sha,input=mesh,
        binding=_mixed_msh_factory_bindings(mesh)))
    return mesh,merge(identities,(factory_events=factory,))
end
function _mixed_msh_factory_binding(mesh,receipt)
    return receipt!==nothing && receipt.input===mesh &&
        isequal(_mixed_msh_factory_bindings(mesh),receipt.binding)
end
