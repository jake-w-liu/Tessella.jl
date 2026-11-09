# ---- ONELAB solver-client protocol over TCP --------------------------------------
#
# Gmsh's ONELAB layer couples a parameter *server* (the gmsh session, started
# with `gmsh -listen host:port` / `Solver.AlwaysListen`) with external solver
# *clients* that connect over a stream socket. This file implements the
# client side of that documented wire protocol — Tessella's counterpart of
# the public standalone `onelab.py`/`onelab` client that ships in the SDK —
# so a Julia session can act as an ONELAB solver client against a real server:
#
#   * framing: every message is `Int32 type, Int32 byteLength, payload` —
#     little-endian on all supported hosts, matching upstream's native `ii`
#     packing on the same platforms;
#   * payload strings are UTF-8; parameter records travel as the version-tagged
#     NUL-separated member stream (`_onelab_wire_encode`/`_onelab_wire_decode`);
#   * the conversation is request/response: the client queries or publishes
#     parameters, sends info/warning/error text, asks the server to launch
#     sub-clients (`GMSH_CONNECT`) and finishes with `GMSH_STOP`.
#
# Async `INFO`/`WARNING`/`ERROR`/`PROGRESS` messages can arrive at any time;
# they are collected in `client.messages` instead of being confused with the
# awaited reply (the reference client prints INFO inline but lets other async
# types fall through — buffering all four keeps the request/response state
# clean on a real server).
#
# `APIOnelabServer.jl` supplies the native server/process role;
# `client_run` below asks its connected server to launch a sub-client.

import Sockets
import Sockets: TCPSocket
using Printf: @sprintf

# Upstream message type codes (`_GMSH_*` in the public client).
const _ONELAB_MSG = (
    START=Int32(1), STOP=Int32(2),
    INFO=Int32(10), WARNING=Int32(11), ERROR=Int32(12), PROGRESS=Int32(13),
    MERGE_FILE=Int32(20), PARSE_STRING=Int32(21),
    PARAMETER=Int32(23), PARAMETER_QUERY=Int32(24),
    PARAMETER_QUERY_ALL=Int32(25), PARAMETER_QUERY_END=Int32(26),
    CONNECT=Int32(27), OLPARSE=Int32(28), PARAMETER_NOT_FOUND=Int32(29),
    PARAMETER_CLEAR=Int32(31), PARAMETER_UPDATE=Int32(32),
    OPEN_PROJECT=Int32(33), CLIENT_CHANGED=Int32(34))
const _ONELAB_WIRE_VERSION="1.3"
const _ONELAB_ASYNC=(_ONELAB_MSG.INFO,_ONELAB_MSG.WARNING,_ONELAB_MSG.ERROR,
                     _ONELAB_MSG.PROGRESS)

# ---- wire parameter record -------------------------------------------------------
#
# One parameter on the wire is a version-tagged NUL-separated member stream.
# Member order is fixed by the protocol (the reference client's `_members`
# table): the shared base members first, then the type-specific tail.
#
#   base:   name, label, help, changedValue, visible, readOnly,
#           attributes (dict string->string), clients (dict string->int)
#   number: values (list float), min, max, step, index,
#           choices (list float), labels (dict float->string)
#   string: values (list string), kind, choices (list string)
#
# Lists encode as a decimal count followed by the items; dicts as a count
# followed by alternating key/value items; floats use "%.16g".

mutable struct _OnelabWireParameter
    kind::Symbol                       # :number or :string
    name::String
    label::String
    help::String
    changed_value::Int
    visible::Bool
    read_only::Bool
    attributes::Dict{String,String}
    attr_order::Vector{String}         # wire/insertion order — upstream emits
    clients::Dict{String,Int}          #   dicts in insertion order, so the
    client_order::Vector{String}       #   parallel key lists keep decode ->
    value_labels::Dict{Float64,String} #   encode byte-exact.
    label_order::Vector{Float64}
    values_f::Vector{Float64}
    values_s::Vector{String}
    min::Float64
    max::Float64
    step::Float64
    index::Int
    choices_f::Vector{Float64}
    choices_s::Vector{String}
    kind_name::String
end

function _OnelabWireParameter(kind::Symbol,name::AbstractString)
    kind in (:number,:string) || throw(ArgumentError(
        "onelab wire parameter kind must be :number or :string"))
    p=_OnelabWireParameter(kind,String(name),"","",31,true,false,
        Dict{String,String}(),String[],
        Dict{String,Int}(),String[],
        Dict{Float64,String}(),Float64[],
        kind===:number ? Float64[0.0] : Float64[],
        kind===:string ? String[""] : String[],
        -1.0e200,1.0e200,0.0,-1,Float64[],String[],
        kind===:number ? "number" : "generic")
    return p
end

# Insertion-ordered dict insert — the wire preserves the sender's order.
@inline function _onelab_dict_put!(dict,order,k,v)
    haskey(dict,k) || push!(order,k)
    dict[k]=v
    return dict
end

# Iterate `order` first, then any keys added by direct dict mutation.
@inline function _onelab_dict_items(dict,order)
    seen=Set{eltype(order)}(order)
    return Iterators.flatten(((Pair(k,dict[k]) for k in order if haskey(dict,k)),
                              (Pair(k,dict[k]) for k in keys(dict) if !(k in seen))))
end

# `%.16g` formatting for the wire (the reference client's float conversion).
function _onelab_wire_f64(v::Real)
    f=Float64(v)
    isfinite(f) || throw(ArgumentError(
        "onelab wire: non-finite parameter value"))
    return @sprintf("%.16g",f)
end

# UTF-8 validity without the generic `isvalid(::AbstractString)` character
# iteration. Under `--compile=min` that fallback allocates per character, so a
# 4 MiB wire value took seconds to validate before the pump could observe its
# deadline. ASCII is checked by Base's fast `isascii(::String)`; other strings
# go through a UTF-aware regex, which is strict (rejects overlongs, surrogates
# and out-of-range code points) and scans 4 MiB in milliseconds.
const _ONELAB_UTF8_CHECK=r"^[\s\S]*$"
_onelab_valid_text(s::AbstractString)=isascii(s) || occursin(_ONELAB_UTF8_CHECK,s)

function _onelab_wire_encode(p::_OnelabWireParameter)
    parts=String[_ONELAB_WIRE_VERSION,p.kind===:number ? "number" : "string"]
    push!(parts,p.name,p.label,p.help,
          string(p.changed_value),string(p.visible ? 1 : 0),
          string(p.read_only ? 1 : 0))
    push!(parts,string(length(p.attributes)))
    for kv in _onelab_dict_items(p.attributes,p.attr_order)
        push!(parts,kv.first,kv.second)
    end
    push!(parts,string(length(p.clients)))
    for kv in _onelab_dict_items(p.clients,p.client_order)
        push!(parts,kv.first,string(kv.second))
    end
    if p.kind===:number
        push!(parts,string(length(p.values_f)))
        append!(parts,(_onelab_wire_f64(v) for v in p.values_f))
        push!(parts,_onelab_wire_f64(p.min),_onelab_wire_f64(p.max),
              _onelab_wire_f64(p.step),string(p.index))
        push!(parts,string(length(p.choices_f)))
        append!(parts,(_onelab_wire_f64(v) for v in p.choices_f))
        push!(parts,string(length(p.value_labels)))
        for kv in _onelab_dict_items(p.value_labels,p.label_order)
            push!(parts,_onelab_wire_f64(kv.first),kv.second)
        end
    else
        push!(parts,string(length(p.values_s)))
        append!(parts,p.values_s)
        push!(parts,p.kind_name)
        push!(parts,string(length(p.choices_s)))
        append!(parts,p.choices_s)
    end
    all(part->_onelab_valid_text(part) && !occursin('\0',part),parts) || throw(ArgumentError(
        "onelab wire: parameter strings must be valid NUL-free UTF-8"))
    return join(parts,'\0')
end

mutable struct _OnelabWireCursor
    items::Vector{SubString{String}}
    position::Int
end

function _onelab_wire_pop!(items,caller::AbstractString,what::AbstractString)
    items.position<=length(items.items) || throw(ArgumentError(
        "$caller: truncated ONELAB parameter payload ($what)"))
    result=String(items.items[items.position]);items.position+=1
    return result
end

function _onelab_wire_count!(items,caller,what;width=1)
    count=parse(Int,_onelab_wire_pop!(items,caller,what))
    0<=count<=(length(items.items)-items.position+1)÷width || throw(ArgumentError(
        "$caller: invalid ONELAB parameter count ($what)"))
    return count
end

function _onelab_wire_pop_list!(items,caller,what,conv)
    n=_onelab_wire_count!(items,caller,what)
    out=conv===nothing ? String[] :
        conv===Float64 ? Float64[] : Any[]
    for _ in 1:n
        raw=_onelab_wire_pop!(items,caller,what)
        conv===nothing && (push!(out,String(raw));continue)
        v=parse(Float64,raw)
        push!(out,v)
    end
    return out
end

function _onelab_wire_decode(msg::AbstractString,caller::AbstractString)
    _onelab_valid_text(msg) || throw(ArgumentError("$caller: ONELAB parameter is not valid UTF-8"))
    items=_OnelabWireCursor(split(String(msg),'\0';keepempty=true),1)
    version=_onelab_wire_pop!(items,caller,"version")
    version==_ONELAB_WIRE_VERSION || throw(ArgumentError(
        "$caller: ONELAB protocol version $(repr(version)) (expected "*
        repr(_ONELAB_WIRE_VERSION)*")"))
    type=_onelab_wire_pop!(items,caller,"type")
    type in ("number","string") || throw(ArgumentError(
        "$caller: unsupported ONELAB parameter type $(repr(type))"))
    p=_OnelabWireParameter(Symbol(type),
        _onelab_wire_pop!(items,caller,"name"))
    p.label=_onelab_wire_pop!(items,caller,"label")
    p.help=_onelab_wire_pop!(items,caller,"help")
    p.changed_value=parse(Int,_onelab_wire_pop!(items,caller,"changedValue"))
    p.visible=parse(Int,_onelab_wire_pop!(items,caller,"visible"))!=0
    p.read_only=parse(Int,_onelab_wire_pop!(items,caller,"readOnly"))!=0
    nattr=_onelab_wire_count!(items,caller,"attributes";width=2)
    for _ in 1:nattr
        k=_onelab_wire_pop!(items,caller,"attributes")
        _onelab_dict_put!(p.attributes,p.attr_order,k,
            _onelab_wire_pop!(items,caller,"attributes"))
    end
    ncli=_onelab_wire_count!(items,caller,"clients";width=2)
    for _ in 1:ncli
        k=_onelab_wire_pop!(items,caller,"clients")
        _onelab_dict_put!(p.clients,p.client_order,k,
            parse(Int,_onelab_wire_pop!(items,caller,"clients")))
    end
    if p.kind===:number
        p.values_f=_onelab_wire_pop_list!(items,caller,"values",Float64)
        p.min=parse(Float64,_onelab_wire_pop!(items,caller,"min"))
        p.max=parse(Float64,_onelab_wire_pop!(items,caller,"max"))
        p.step=parse(Float64,_onelab_wire_pop!(items,caller,"step"))
        p.index=parse(Int,_onelab_wire_pop!(items,caller,"index"))
        p.choices_f=_onelab_wire_pop_list!(items,caller,"choices",Float64)
        nlab=_onelab_wire_count!(items,caller,"labels";width=2)
        for _ in 1:nlab
            k=parse(Float64,_onelab_wire_pop!(items,caller,"labels"))
            _onelab_dict_put!(p.value_labels,p.label_order,k,
                _onelab_wire_pop!(items,caller,"labels"))
        end
    else
        p.values_s=_onelab_wire_pop_list!(items,caller,"values",nothing)
        p.kind_name=_onelab_wire_pop!(items,caller,"kind")
        p.choices_s=_onelab_wire_pop_list!(items,caller,"choices",nothing)
    end
    items.position==length(items.items)+1 || throw(ArgumentError(
        "$caller: trailing ONELAB parameter members"))
    return p
end

# ---- client object -----------------------------------------------------------------

"""An ONELAB solver client connected to a parameter server over TCP — the
client half of the documented upstream socket protocol. Construct with
`API.onelab.client_connect`; close with `API.onelab.client_disconnect`."""
mutable struct OnelabClient
    socket::Union{TCPSocket,Nothing}
    addr::String                    # "host:port" for the reconnect arm
    name::String
    messages::Vector{Tuple{Int32,String}}  # buffered async INFO/WARNING/ERROR/PROGRESS
    subclients::Int                 # outstanding GMSH_CONNECT sub-clients
    OnelabClient()=new(nothing,"","",Tuple{Int32,String}[],0)
end

function _onelab_client_check(c,caller::AbstractString)
    c isa OnelabClient || throw(ArgumentError("$caller: not an ONELAB client"))
    (c.socket===nothing || !Sockets.isopen(c.socket)) && throw(ArgumentError(
        "$caller: ONELAB client is not connected"))
    return nothing
end

function _onelab_client_parse_addr(addr,port,caller::AbstractString)
    if port===nothing
        addr isa AbstractString || throw(ArgumentError(
            "$caller: addr must be a \"host:port\" string"))
        s=String(addr)
        i=findlast(':',s)
        i===nothing && throw(ArgumentError(
            "$caller: addr $(repr(s)) is not a \"host:port\" string"))
        host=s[1:prevind(s,i)]
        p=tryparse(Int,s[nextind(s,i):end])
        (p===nothing || p<0 || p>65535) && throw(ArgumentError(
            "$caller: addr $(repr(s)) has no valid port"))
        return host,p,s
    end
    host=addr isa AbstractString ? String(addr) :
        addr isa Sockets.IPAddr ? string(addr) :
        throw(ArgumentError("$caller: host must be a string or IPAddr"))
    p=Int(port)
    (0<=p<=65535) || throw(ArgumentError("$caller: invalid port $p"))
    return host,p,host*":"*string(p)
end

function _onelab_client_send(c::OnelabClient,type::Int32,
                             payload::AbstractString="")
    c.socket===nothing && return nothing
    data=codeunits(String(payload))
    frame=IOBuffer()
    write(frame,htol(Int32(type)));write(frame,htol(Int32(length(data))))
    write(frame,data)
    bytes=take!(frame)
    for attempt in 1:2
        try
            write(c.socket,bytes)
            flush(c.socket)
            return nothing
        catch err
            attempt==2 && rethrow()
            # Upstream retries the send once on a fresh socket; the frame
            # bytes live in `bytes`, which survives the failed write.
            try
                Sockets.close(c.socket)
                host,port,_=_onelab_client_parse_addr(c.addr,nothing,
                    "onelab client reconnect")
                c.socket=Sockets.connect(host,port)
            catch
                rethrow(err)
            end
        end
    end
    return nothing
end

# Blocking exact-length read: `read(sock,UInt8)` waits for data and throws
# `EOFError` once the peer closes the stream.
function _onelab_client_recv_exact(sock::TCPSocket,n::Int)
    n>=0 || throw(ArgumentError("onelab client: negative byte count"))
    buf=UInt8[];sizehint!(buf,min(n,4096))
    while length(buf)<n
        count=min(n-length(buf),max(1,bytesavailable(sock)),65536)
        data=read(sock,count)
        isempty(data) && throw(EOFError())
        append!(buf,data)
    end
    return buf
end

# Receive one framed message. Async INFO/WARNING/ERROR/PROGRESS payloads are
# buffered on the client and the read continues until a reply-type frame —
# the server may interleave them at any point.
function _onelab_client_receive(c::OnelabClient)
    sock=c.socket
    sock===nothing && throw(ArgumentError("onelab client: not connected"))
    while true
        hdr=_onelab_client_recv_exact(sock,8)
        type=ltoh(reinterpret(Int32,hdr[1:4])[1])
        len=ltoh(reinterpret(Int32,hdr[5:8])[1])
        len<0 && throw(ArgumentError(
            "onelab client: negative frame length $len"))
        payload=len==0 ? "" :
            String(_onelab_client_recv_exact(sock,Int(len)))
        if type in _ONELAB_ASYNC
            push!(c.messages,(type,payload))
            continue
        end
        return type,payload
    end
end

# ---- parameter exchange ---------------------------------------------------------

# Upstream `defineParameter`: QUERY then either PARAMETER (server copy wins —
# we reply PARAMETER_UPDATE and adopt the server's record) or
# PARAMETER_NOT_FOUND (define by sending PARAMETER). Returns the server's
# parameter record.
function _onelab_client_define(c::OnelabClient,p::_OnelabWireParameter,
                               caller::AbstractString)
    _onelab_client_send(c,_ONELAB_MSG.PARAMETER_QUERY,_onelab_wire_encode(p))
    type,msg=_onelab_client_receive(c)
    if type==_ONELAB_MSG.PARAMETER
        server=_onelab_wire_decode(msg,caller)
        _onelab_client_send(c,_ONELAB_MSG.PARAMETER_UPDATE,
            _onelab_wire_encode(server))
        return server
    elseif type==_ONELAB_MSG.PARAMETER_NOT_FOUND
        _onelab_client_send(c,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(p))
        return p
    end
    throw(ArgumentError(
        "$caller: unexpected ONELAB reply type $type"))
end

# QUERY then decode PARAMETER into `p`; returns `false` on NOT_FOUND.
function _onelab_client_fetch(c::OnelabClient,p::_OnelabWireParameter,
                              caller::AbstractString)
    _onelab_client_send(c,_ONELAB_MSG.PARAMETER_QUERY,_onelab_wire_encode(p))
    type,msg=_onelab_client_receive(c)
    if type==_ONELAB_MSG.PARAMETER
        q=_onelab_wire_decode(msg,caller)
        for f in fieldnames(_OnelabWireParameter); setfield!(p,f,getfield(q,f)); end
        return true
    elseif type==_ONELAB_MSG.PARAMETER_NOT_FOUND
        return false
    end
    throw(ArgumentError("$caller: unexpected ONELAB reply type $type"))
end

# ---- public endpoint helpers -----------------------------------------------------

function _onelab_client_connect(addr,port;name="tessella client")
    caller="API.onelab.client_connect"
    host,p,canon=_onelab_client_parse_addr(addr,port,caller)
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    sock=Sockets.connect(host,p)
    c=OnelabClient()
    c.socket=sock;c.addr=canon;c.name=String(name)
    try
        _onelab_client_send(c,_ONELAB_MSG.START,string(getpid()))
    catch err
        close(sock)
        rethrow()
    end
    return c
end

function _onelab_client_disconnect(c)
    caller="API.onelab.client_disconnect"
    c isa OnelabClient || throw(ArgumentError("$caller: not an ONELAB client"))
    if c.socket!==nothing && Sockets.isopen(c.socket)
        try
            _onelab_client_wait(c)
            _onelab_client_send(c,_ONELAB_MSG.STOP,"Goodbye!")
        catch
        end
        Sockets.close(c.socket)
    end
    c.socket=nothing
    return nothing
end

function _onelab_client_get(c,kind::Symbol,name;
                            warn::Bool=true)
    caller=kind===:number ? "API.onelab.client_get_number" :
        "API.onelab.client_get_string"
    _onelab_client_check(c,caller)
    p=_OnelabWireParameter(kind,name)
    found=_onelab_client_fetch(c,p,caller)
    if !found
        warn && push!(c.messages,(_ONELAB_MSG.WARNING,
            "Unknown parameter $name"))
        return nothing
    end
    return p
end

function _onelab_client_get_number(c,name;warn=true)
    p=_onelab_client_get(c,:number,name;warn)
    p===nothing && return nothing
    return copy(p.values_f)
end
function _onelab_client_get_string(c,name;warn=true)
    p=_onelab_client_get(c,:string,name;warn)
    p===nothing && return nothing
    return copy(p.values_s)
end

function _onelab_client_get_parameter(c,kind::Symbol,name)
    caller="API.onelab.client_get_parameter"
    _onelab_client_check(c,caller)
    p=_OnelabWireParameter(kind,name)
    _onelab_client_fetch(c,p,caller) || return nothing
    return p
end

function _onelab_client_set_number(c,name,values)
    caller="API.onelab.client_set_number"
    _onelab_client_check(c,caller)
    vals=_onelab_floats(values,caller)
    p=_OnelabWireParameter(:number,name)
    _onelab_client_fetch(c,p,caller)   # server copy wins when present
    p.values_f=vals
    _onelab_client_send(c,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(p))
    return nothing
end

function _onelab_client_set_string(c,name,values)
    caller="API.onelab.client_set_string"
    _onelab_client_check(c,caller)
    vals=_onelab_strings(values,caller)
    p=_OnelabWireParameter(:string,name)
    _onelab_client_fetch(c,p,caller)
    p.values_s=vals
    _onelab_client_send(c,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(p))
    return nothing
end

# Apply keyword metadata to a wire parameter (define_* kwargs).
function _onelab_client_apply_kwargs(p::_OnelabWireParameter,kwargs,caller)
    for (k,v) in kwargs
        if k===:values
            if p.kind===:number
                p.values_f=_onelab_floats(v,caller)
            else
                p.values_s=_onelab_strings(v,caller)
            end
        elseif k===:value
            if p.kind===:number
                p.values_f=_onelab_floats(v,caller)
            else
                p.values_s=_onelab_strings(v,caller)
            end
        elseif k===:min
            p.min=Float64(v)
        elseif k===:max
            p.max=Float64(v)
        elseif k===:step
            p.step=Float64(v)
        elseif k===:index
            p.index=Int(v)
        elseif k===:label
            p.label=String(v)
        elseif k===:help || k===:short_help || k===:long_help
            p.help=String(v)
        elseif k===:visible
            p.visible=Bool(v)
        elseif k===:read_only || k===:readOnly
            p.read_only=Bool(v)
        elseif k===:attributes
            for (ak,av) in pairs(v)
                _onelab_dict_put!(p.attributes,p.attr_order,
                    String(ak),String(av))
            end
        elseif k===:clients
            for (ak,av) in pairs(v)
                _onelab_dict_put!(p.clients,p.client_order,
                    String(ak),Int(av))
            end
        elseif k===:choices
            if p.kind===:number
                p.choices_f=_onelab_floats(v,caller)
            else
                p.choices_s=_onelab_strings(v,caller)
            end
        elseif k===:labels || k===:value_labels
            for (ak,av) in pairs(v)
                _onelab_dict_put!(p.value_labels,p.label_order,
                    Float64(ak),String(av))
            end
        elseif k===:kind
            p.kind_name=String(v)
        else
            throw(ArgumentError(
                "$caller: unknown parameter attribute $k"))
        end
    end
    return p
end

function _onelab_client_define_number(c,name;kwargs...)
    caller="API.onelab.client_define_number"
    _onelab_client_check(c,caller)
    p=_OnelabWireParameter(:number,name)
    _onelab_client_apply_kwargs(p,Dict(kwargs),caller)
    server=_onelab_client_define(c,p,caller)
    return copy(server.values_f)
end

function _onelab_client_define_string(c,name;kwargs...)
    caller="API.onelab.client_define_string"
    _onelab_client_check(c,caller)
    p=_OnelabWireParameter(:string,name)
    _onelab_client_apply_kwargs(p,Dict(kwargs),caller)
    server=_onelab_client_define(c,p,caller)
    return copy(server.values_s)
end

function _onelab_client_get_all(c)
    caller="API.onelab.client_get_all"
    _onelab_client_check(c,caller)
    out=_OnelabWireParameter[]
    for kind in (:number,:string)
        _onelab_client_send(c,_ONELAB_MSG.PARAMETER_QUERY_ALL,
            _onelab_wire_encode(_OnelabWireParameter(kind,"dummy")))
        type,msg=_onelab_client_receive(c)
        while type==_ONELAB_MSG.PARAMETER_QUERY_ALL
            push!(out,_onelab_wire_decode(msg,caller))
            type,msg=_onelab_client_receive(c)
        end
    end
    return out
end

function _onelab_client_get_names(c;search="")
    caller="API.onelab.client_get_names"
    search isa AbstractString || throw(ArgumentError(
        "$caller: search must be a string"))
    params=_onelab_client_get_all(c)
    return sort!(String[p.name for p in params
                        if _onelab_glob_match(p.name,String(search))])
end

function _onelab_client_clear(c,name)
    caller="API.onelab.client_clear"
    _onelab_client_check(c,caller)
    _onelab_client_send(c,_ONELAB_MSG.PARAMETER_CLEAR,
        _onelab_check_name(name,caller))
    return nothing
end

function _onelab_client_get_changed(c,name)
    caller="API.onelab.client_get_changed"
    _onelab_client_check(c,caller)
    key=_onelab_check_name(name,caller)
    _onelab_client_send(c,_ONELAB_MSG.CLIENT_CHANGED,join(("get",key),'\0'))
    type,msg=_onelab_client_receive(c)
    type==_ONELAB_MSG.CLIENT_CHANGED || throw(ArgumentError(
        "$caller: unexpected ONELAB reply type $type"))
    return msg=="true"
end

function _onelab_client_set_changed(c,name,changed)
    caller="API.onelab.client_set_changed"
    _onelab_client_check(c,caller)
    key=_onelab_check_name(name,caller)
    _onelab_client_send(c,_ONELAB_MSG.CLIENT_CHANGED,
        join(("set",key,changed ? "true" : "false"),'\0'))
    return nothing
end

# ---- commands, files, sub-clients -------------------------------------------------

function _onelab_client_send_command(c,command)
    caller="API.onelab.client_send_command"
    _onelab_client_check(c,caller)
    cmd=_onelab_check_name(command,caller)
    _onelab_client_send(c,_ONELAB_MSG.PARSE_STRING,cmd)
    return nothing
end

function _onelab_client_merge_file(c,path)
    caller="API.onelab.client_merge_file"
    _onelab_client_check(c,caller)
    _onelab_client_send(c,_ONELAB_MSG.MERGE_FILE,String(path))
    return nothing
end

function _onelab_client_open_project(c,path)
    caller="API.onelab.client_open_project"
    _onelab_client_check(c,caller)
    _onelab_client_send(c,_ONELAB_MSG.OPEN_PROJECT,String(path))
    return nothing
end

function _onelab_client_send_message(c,type::Int32,msg,caller)
    _onelab_client_check(c,caller)
    _onelab_client_send(c,type,String(msg))
    return nothing
end
_onelab_client_info(c,msg)=_onelab_client_send_message(c,_ONELAB_MSG.INFO,msg,
    "API.onelab.client_info")
_onelab_client_warning(c,msg)=_onelab_client_send_message(c,_ONELAB_MSG.WARNING,
    msg,"API.onelab.client_warning")
_onelab_client_error(c,msg)=_onelab_client_send_message(c,_ONELAB_MSG.ERROR,
    msg,"API.onelab.client_error")

# Server-side .ol preprocessing request (`OLPARSE`). The reply reports whether
# the template resolved a parameter; plain text can be written with false.
function _onelab_client_preprocess(c,name,filename)
    caller="API.onelab.client_preprocess"
    _onelab_client_check(c,caller)
    _onelab_client_send(c,_ONELAB_MSG.OLPARSE,
        join((String(name),String(filename)),'\0'))
    type,msg=_onelab_client_receive(c)
    return type==_ONELAB_MSG.OLPARSE && msg=="true"
end

# Ask the server to launch a sub-client (upstream `runSubClient`): the request
# travels as `name \0 command` and each outstanding sub-client is reaped by a
# `STOP` frame. `block=false` leaves the reaping to `client_wait`.
function _onelab_client_run(c,name,command;arguments="",block=true)
    caller="API.onelab.client_run"
    _onelab_client_check(c,caller)
    client=_onelab_check_name(name,caller)
    cmd=_onelab_check_name(command,caller)
    isempty(client) && throw(ArgumentError(
        "$caller: sub-client name must be nonempty"))
    isempty(cmd) && throw(ArgumentError(
        "$caller: sub-client command must be nonempty"))
    args=String(arguments)
    payload=isempty(args) ? join((client,cmd),'\0') :
        join((client,cmd*" "*args),'\0')
    _onelab_client_send(c,_ONELAB_MSG.CONNECT,payload)
    c.subclients+=1
    block && _onelab_client_wait(c)
    return nothing
end

# Reap all outstanding sub-clients (one STOP each).
function _onelab_client_wait(c)
    caller="API.onelab.client_wait"
    _onelab_client_check(c,caller)
    while c.subclients>0
        type,msg=_onelab_client_receive(c)
        type==_ONELAB_MSG.STOP && (c.subclients-=1)
    end
    return nothing
end

# Read the client's Action string parameter — the upstream solver-client
# lifecycle (`initialize`/`check`/`compute`/`stop`). `checkParameters` is the
# `check` arm: the client defines its parameter records against the server.
function _onelab_client_action(c)
    caller="API.onelab.client_action"
    _onelab_client_check(c,caller)
    p=_onelab_client_get_parameter(c,:string,c.name*"/Action")
    p===nothing && return ""
    return isempty(p.values_s) ? "" : p.values_s[1]
end

# A convenience `check` step: define a batch of parameter records against the
# server (each `(name, kwargs...)` like `client_define_*`). Returns a Dict of
# adopted server values keyed by name.
function _onelab_client_check_parameters(c,defs)
    caller="API.onelab.client_check_parameters"
    _onelab_client_check(c,caller)
    out=Dict{String,Vector}()
    for d in defs
        if d isa Pair
            name,spec=d
        else
            name=d
            spec=nothing
        end
        kind=:number;kwargs=Dict{Symbol,Any}()
        if spec isa NamedTuple
            kwargs=Dict{Symbol,Any}(pairs(spec))
        elseif spec isa Dict
            kwargs=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in spec)
        elseif spec===nothing
        else
            throw(ArgumentError(
                "$caller: parameter spec must be a NamedTuple or Dict"))
        end
        haskey(kwargs,:kind) && (kind=kwargs[:kind]===:string ? :string : :number;
                                 delete!(kwargs,:kind))
        p=_OnelabWireParameter(kind,String(name))
        _onelab_client_apply_kwargs(p,kwargs,caller)
        server=_onelab_client_define(c,p,caller)
        out[p.name]=kind===:number ? copy(server.values_f) :
            copy(server.values_s)
    end
    return out
end

# Drain currently-buffered async messages (INFO/WARNING/ERROR/PROGRESS) as
# `(type, text)` pairs; upstream prints INFO inline.
function _onelab_client_messages(c)
    caller="API.onelab.client_messages"
    c isa OnelabClient || throw(ArgumentError("$caller: not an ONELAB client"))
    out=copy(c.messages)
    empty!(c.messages)
    return out
end
