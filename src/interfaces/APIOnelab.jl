# ---- session ONELAB parameter store -------------------------------------------
#
# Gmsh 4.15.2 `onelab` module support. Number and string parameters have
# independent namespaces: the same name can hold a record of each type.
# `SESSION_ONELAB_CHANGED` is the compatibility changed-name flag
# set that solver clients poll through `set_changed`/`get_changed`/
# `clear_changed`. Both tables are emptied by `initialize`/`finalize`;
# `clear` preserves the parameter space like upstream's model/view-only reset.
# All state reads/mutations run under `STATE_LOCK`.
#
# `.geo` interop: `execute_geo`'s `_GeoNumericContext.onelab_*` dicts are
# scalar (first-value) mirrors of an operation-owned full parameter store.
# Existing values win `DefineNumber` defaults. The native process server and
# GEO execution share this detached store, whose values, metadata and client
# registry publish only after the file frontend's model preparation succeeds.
#
# The JSON payload is the ONELAB parameter array: objects with "type"
# ("number"/"string"), "name", "values" and optional "min"/"max"/"step",
# "choices", "valueLabels", "attributes", "label", "shortHelp", "longHelp",
# "readOnly", "visible", "closed" and "kind" members. The codec below is a
# minimal self-contained JSON implementation (the package has no JSON
# dependency); unknown object members are ignored on input, matching the
# upstream fromJSON convention.

mutable struct _OnelabParameter
    name::String
    kind::Symbol                     # :number or :string
    values_f::Vector{Float64}        # number values
    values_s::Vector{String}         # string values
    min::Union{Nothing,Float64}
    max::Union{Nothing,Float64}
    step::Union{Nothing,Float64}
    choices_f::Vector{Float64}
    choices_s::Vector{String}
    # "valueLabels" number-key/label members preserved verbatim (keys arrive
    # as JSON object member names, i.e. strings).
    value_labels::Dict{String,String}
    attributes::Dict{String,String}
    label::Union{Nothing,String}
    short_help::Union{Nothing,String}
    long_help::Union{Nothing,String}
    read_only::Bool
    visible::Bool
    closed::Bool
    # Upstream "kind" member — "parameter" is the ordinary default.
    kind_name::String
end

_OnelabParameter(name::String,kind::Symbol)=_OnelabParameter(
    name,kind,Float64[],String[],nothing,nothing,nothing,
    Float64[],String[],Dict{String,String}(),Dict{String,String}(),
    nothing,nothing,nothing,false,true,false,"parameter")

const _OnelabParameterKey=Tuple{Symbol,String}
const SESSION_ONELAB = Dict{_OnelabParameterKey,_OnelabParameter}()
const SESSION_ONELAB_CHANGED = Set{String}()

function _onelab_preferred_key(parameters,name)
    for kind in (:number,:string)
        key=(kind,name)
        haskey(parameters,key) && return key
    end
    return nothing
end

# ---- minimal JSON decoder ------------------------------------------------------

_onelab_json_err(caller::AbstractString,msg::AbstractString)=
    throw(ArgumentError("$caller: $msg"))

function _onelab_json_skipws(s::String,i::Int,n::Int)
    while i<=n
        c=s[i]
        (c==' ' || c=='\t' || c=='\n' || c=='\r') || break
        i=nextind(s,i)
    end
    return i
end

function _onelab_json_literal(s::String,i::Int,n::Int,caller::AbstractString,
                              word::String,val)
    last=i+sizeof(word)-1
    (last<=n && SubString(s,i,last)==word) ||
        _onelab_json_err(caller,"invalid JSON literal near byte $i")
    return val,last+1
end

function _onelab_json_hex4(s::String,i::Int,n::Int,caller::AbstractString)
    value=UInt32(0)
    for _ in 1:4
        i<=n || _onelab_json_err(caller,"truncated \\u escape in JSON string")
        c=s[i]
        d='0'<=c<='9' ? c-'0' : 'a'<=c<='f' ? c-'a'+10 :
           'A'<=c<='F' ? c-'A'+10 :
           _onelab_json_err(caller,"invalid \\u escape in JSON string")
        value=(value<<4)|UInt32(d)
        i=nextind(s,i)
    end
    return value,i
end

function _onelab_json_string(s::String,i::Int,n::Int,caller::AbstractString)
    # s[i] == '"'; returns (String, next_index).
    io=IOBuffer();i=nextind(s,i)
    while i<=n
        c=s[i]
        if c=='"'
            return String(take!(io)),nextind(s,i)
        elseif c=='\\'
            i=nextind(s,i)
            i<=n || _onelab_json_err(caller,"unterminated escape in JSON string")
            e=s[i]
            if e=='"'||e=='\\'||e=='/'
                write(io,e);i=nextind(s,i)
            elseif e=='b';write(io,'\b');i=nextind(s,i)
            elseif e=='f';write(io,'\f');i=nextind(s,i)
            elseif e=='n';write(io,'\n');i=nextind(s,i)
            elseif e=='r';write(io,'\r');i=nextind(s,i)
            elseif e=='t';write(io,'\t');i=nextind(s,i)
            elseif e=='u'
                cp,i=_onelab_json_hex4(s,nextind(s,i),n,caller)
                if 0xD800<=cp<=0xDBFF
                    # High surrogate — the low surrogate must follow as \uXXXX.
                    (i<=n && s[i]=='\\') || _onelab_json_err(caller,
                        "lone high surrogate in JSON string")
                    j=nextind(s,i)
                    (j<=n && s[j]=='u') || _onelab_json_err(caller,
                        "lone high surrogate in JSON string")
                    lo,j=_onelab_json_hex4(s,nextind(s,j),n,caller)
                    0xDC00<=lo<=0xDFFF || _onelab_json_err(caller,
                        "invalid low surrogate in JSON string")
                    cp=0x10000+((cp-0xD800)<<10)+(lo-0xDC00)
                    i=j
                elseif 0xDC00<=cp<=0xDFFF
                    _onelab_json_err(caller,"lone low surrogate in JSON string")
                end
                write(io,Char(cp))
            else
                _onelab_json_err(caller,"invalid escape '\\$e' in JSON string")
            end
        elseif c<' '
            _onelab_json_err(caller,"unescaped control character in JSON string")
        else
            write(io,c);i=nextind(s,i)
        end
    end
    _onelab_json_err(caller,"unterminated JSON string")
end

function _onelab_json_number(s::String,i::Int,n::Int,caller::AbstractString)
    j=i
    (j<=n && s[j]=='-') && (j=nextind(s,j))
    (j<=n && '0'<=s[j]<='9') ||
        _onelab_json_err(caller,"invalid JSON number near byte $i")
    while j<=n && '0'<=s[j]<='9';j=nextind(s,j);end
    if j<=n && s[j]=='.'
        j=nextind(s,j)
        (j<=n && '0'<=s[j]<='9') ||
            _onelab_json_err(caller,"invalid JSON number near byte $i")
        while j<=n && '0'<=s[j]<='9';j=nextind(s,j);end
    end
    if j<=n && (s[j]=='e'||s[j]=='E')
        j=nextind(s,j)
        (j<=n && (s[j]=='+'||s[j]=='-')) && (j=nextind(s,j))
        (j<=n && '0'<=s[j]<='9') ||
            _onelab_json_err(caller,"invalid JSON number near byte $i")
        while j<=n && '0'<=s[j]<='9';j=nextind(s,j);end
    end
    return parse(Float64,SubString(s,i,prevind(s,j))),j
end

function _onelab_json_object(s::String,i::Int,n::Int,caller::AbstractString)
    # s[i] == '{'; returns (Dict{String,Any}, next_index).
    object=Dict{String,Any}()
    i=_onelab_json_skipws(s,nextind(s,i),n)
    (i<=n && s[i]=='}') && return object,nextind(s,i)
    while true
        i=_onelab_json_skipws(s,i,n)
        (i<=n && s[i]=='"') ||
            _onelab_json_err(caller,"JSON object key must be a string")
        key,i=_onelab_json_string(s,i,n,caller)
        i=_onelab_json_skipws(s,i,n)
        (i<=n && s[i]==':') ||
            _onelab_json_err(caller,"expected ':' in JSON object")
        value,i=_onelab_json_value(s,i+1,n,caller)
        object[key]=value
        i=_onelab_json_skipws(s,i,n)
        i<=n || _onelab_json_err(caller,"unterminated JSON object")
        s[i]==',' && (i=nextind(s,i);continue)
        s[i]=='}' && return object,nextind(s,i)
        _onelab_json_err(caller,"expected ',' or '}' in JSON object")
    end
end

function _onelab_json_array(s::String,i::Int,n::Int,caller::AbstractString)
    # s[i] == '['; returns (Vector{Any}, next_index).
    values=Any[]
    i=_onelab_json_skipws(s,nextind(s,i),n)
    (i<=n && s[i]==']') && return values,nextind(s,i)
    while true
        value,i=_onelab_json_value(s,i,n,caller)
        push!(values,value)
        i=_onelab_json_skipws(s,i,n)
        i<=n || _onelab_json_err(caller,"unterminated JSON array")
        s[i]==',' && (i=nextind(s,i);continue)
        s[i]==']' && return values,nextind(s,i)
        _onelab_json_err(caller,"expected ',' or ']' in JSON array")
    end
end

function _onelab_json_value(s::String,i::Int,n::Int,caller::AbstractString)
    i=_onelab_json_skipws(s,i,n)
    i<=n || _onelab_json_err(caller,"unexpected end of JSON input")
    c=s[i]
    c=='{' && return _onelab_json_object(s,i,n,caller)
    c=='[' && return _onelab_json_array(s,i,n,caller)
    c=='"' && return _onelab_json_string(s,i,n,caller)
    c=='t' && return _onelab_json_literal(s,i,n,caller,"true",true)
    c=='f' && return _onelab_json_literal(s,i,n,caller,"false",false)
    c=='n' && return _onelab_json_literal(s,i,n,caller,"null",nothing)
    (c=='-' || '0'<=c<='9') && return _onelab_json_number(s,i,n,caller)
    _onelab_json_err(caller,"unexpected character $(repr(c)) in JSON input")
end

function _onelab_json_parse(src::AbstractString,caller::AbstractString)
    s=String(src);n=lastindex(s)
    value,i=_onelab_json_value(s,firstindex(s),n,caller)
    _onelab_json_skipws(s,i,n)<=n &&
        _onelab_json_err(caller,"trailing data after JSON value")
    return value
end

# ---- JSON encoder ----------------------------------------------------------------

function _onelab_json_quote(io::IO,s::AbstractString)
    write(io,'"')
    for c in s
        if c=='"';write(io,"\\\"")
        elseif c=='\\';write(io,"\\\\")
        elseif c=='\b';write(io,"\\b")
        elseif c=='\f';write(io,"\\f")
        elseif c=='\n';write(io,"\\n")
        elseif c=='\r';write(io,"\\r")
        elseif c=='\t';write(io,"\\t")
        elseif c<'\u20'
            print(io,"\\u",lpad(string(UInt32(c),base=16),4,'0'))
        else
            write(io,c)
        end
    end
    write(io,'"')
    return nothing
end

function _onelab_json_number!(io::IO,v::Float64,caller::AbstractString)
    isfinite(v) || _onelab_json_err(caller,
        "parameter values must be JSON-representable (finite)")
    # Gmsh-style integral compactness: 1 instead of 1.0.
    if isinteger(v) && abs(v)<=9.007199254740992e15
        print(io,Int64(v))
    else
        print(io,v)
    end
    return nothing
end

function _onelab_json_string_map!(io::IO,map::Dict{String,String})
    write(io,"{")
    for (i,key) in enumerate(sort!(collect(keys(map))))
        i>1 && write(io,", ")
        _onelab_json_quote(io,key);write(io,": ")
        _onelab_json_quote(io,map[key])
    end
    write(io,"}")
    return nothing
end

# Serialize one parameter record in the ONELAB parameter-array element shape.
# Optional members are emitted only when set or non-default.
function _onelab_param_json!(io::IO,p::_OnelabParameter,caller::AbstractString)
    write(io,"  {\n    \"type\": ")
    _onelab_json_quote(io,p.kind===:number ? "number" : "string")
    write(io,",\n    \"name\": ");_onelab_json_quote(io,p.name)
    write(io,",\n    \"values\": [")
    if p.kind===:number
        for (i,v) in enumerate(p.values_f)
            i>1 && write(io,", ")
            _onelab_json_number!(io,v,caller)
        end
    else
        for (i,v) in enumerate(p.values_s)
            i>1 && write(io,", ")
            _onelab_json_quote(io,v)
        end
    end
    write(io,"]")
    if p.kind===:number
        p.min===nothing || (write(io,",\n    \"min\": ");
            _onelab_json_number!(io,p.min,caller))
        p.max===nothing || (write(io,",\n    \"max\": ");
            _onelab_json_number!(io,p.max,caller))
        p.step===nothing || (write(io,",\n    \"step\": ");
            _onelab_json_number!(io,p.step,caller))
        if !isempty(p.choices_f)
            write(io,",\n    \"choices\": [")
            for (i,v) in enumerate(p.choices_f)
                i>1 && write(io,", ")
                _onelab_json_number!(io,v,caller)
            end
            write(io,"]")
        end
        if !isempty(p.value_labels)
            write(io,",\n    \"valueLabels\": ")
            _onelab_json_string_map!(io,p.value_labels)
        end
    else
        if !isempty(p.choices_s)
            write(io,",\n    \"choices\": [")
            for (i,v) in enumerate(p.choices_s)
                i>1 && write(io,", ")
                _onelab_json_quote(io,v)
            end
            write(io,"]")
        end
    end
    if !isempty(p.attributes)
        write(io,",\n    \"attributes\": ")
        _onelab_json_string_map!(io,p.attributes)
    end
    p.label===nothing || (write(io,",\n    \"label\": ");
        _onelab_json_quote(io,p.label))
    p.short_help===nothing || (write(io,",\n    \"shortHelp\": ");
        _onelab_json_quote(io,p.short_help))
    p.long_help===nothing || (write(io,",\n    \"longHelp\": ");
        _onelab_json_quote(io,p.long_help))
    p.read_only && write(io,",\n    \"readOnly\": true")
    p.visible || write(io,",\n    \"visible\": false")
    p.closed && write(io,",\n    \"closed\": true")
    p.kind_name=="parameter" ||
        (write(io,",\n    \"kind\": ");_onelab_json_quote(io,p.kind_name))
    write(io,"\n  }")
    return nothing
end

# ---- coercion helpers ------------------------------------------------------------

function _onelab_check_name(name,caller::AbstractString)
    name isa AbstractString ||
        throw(ArgumentError("$caller: name must be a string"))
    return String(name)
end

function _onelab_check_format(format,caller::AbstractString)
    format isa AbstractString ||
        throw(ArgumentError("$caller: format must be a string"))
    f=lowercase(String(format))
    f=="json" || throw(ArgumentError(
        "$caller: unsupported ONELAB data format $(repr(format))"))
    return f
end

function _onelab_floats(values,caller::AbstractString)
    if values isa Real
        values isa Bool &&
            throw(ArgumentError("$caller: values must not be Bool"))
        v=Float64(values)
        isfinite(v) ||
            throw(ArgumentError("$caller: values must be finite"))
        return Float64[v]
    end
    values isa AbstractVector || throw(ArgumentError(
        "$caller: values must be a real number or a vector of real numbers"))
    out=Float64[]
    for v in values
        v isa Bool &&
            throw(ArgumentError("$caller: values must not contain Bool"))
        v isa Real || throw(ArgumentError(
            "$caller: values must be real numbers"))
        f=Float64(v)
        isfinite(f) ||
            throw(ArgumentError("$caller: values must be finite"))
        push!(out,f)
    end
    return out
end

function _onelab_strings(values,caller::AbstractString)
    values isa AbstractString && return String[String(values)]
    values isa AbstractVector || throw(ArgumentError(
        "$caller: values must be a string or a vector of strings"))
    out=String[]
    for v in values
        v isa AbstractString ||
            throw(ArgumentError("$caller: values must be strings"))
        push!(out,String(v))
    end
    return out
end

# JSON-decoded counterparts: values arrive as `Vector{Any}`/`Float64`/String.
function _onelab_json_floats(value,caller::AbstractString,what::String)
    value===nothing && return Float64[]
    value isa Vector ||
        value isa Real || throw(ArgumentError(
            "$caller: $what must be an array of numbers"))
    items=value isa Vector ? value : Any[value]
    out=Float64[]
    for v in items
        (v isa Real && !(v isa Bool)) || throw(ArgumentError(
            "$caller: $what must be an array of numbers"))
        f=Float64(v)
        isfinite(f) || throw(ArgumentError(
            "$caller: $what must be finite"))
        push!(out,f)
    end
    return out
end

function _onelab_json_strings(value,caller::AbstractString,what::String)
    value===nothing && return String[]
    value isa Vector ||
        value isa AbstractString || throw(ArgumentError(
            "$caller: $what must be an array of strings"))
    items=value isa Vector ? value : Any[value]
    out=String[]
    for v in items
        v isa AbstractString || throw(ArgumentError(
            "$caller: $what must be an array of strings"))
        push!(out,String(v))
    end
    return out
end

function _onelab_json_scalar(value,caller::AbstractString,what::String)
    (value isa Real && !(value isa Bool)) || throw(ArgumentError(
        "$caller: $what must be a number"))
    f=Float64(value)
    isfinite(f) || throw(ArgumentError("$caller: $what must be finite"))
    return f
end

_onelab_json_bool(value)=value isa Bool ? value : !(value isa Real) ?
    nothing : !iszero(value)
# `nothing` marks "not a bool" for the caller to ignore/error.

function _onelab_json_optstring(value,caller::AbstractString,what::String)
    value===nothing && return nothing
    value isa AbstractString || throw(ArgumentError(
        "$caller: $what must be a string"))
    return String(value)
end

function _onelab_json_string_map(value,caller::AbstractString,what::String)
    value isa Dict || throw(ArgumentError(
        "$caller: $what must be a JSON object"))
    out=Dict{String,String}()
    for (key,v) in value
        key isa AbstractString || throw(ArgumentError(
            "$caller: $what keys must be strings"))
        out[String(key)]=v isa AbstractString ? String(v) : string(v)
    end
    return out
end

# Build a record from one decoded JSON object. Unknown members are ignored
# (upstream fromJSON convention); "min"/"max"/"step" apply to number records.
function _onelab_param_from_json(obj,caller::AbstractString)
    obj isa Dict{String,Any} || throw(ArgumentError(
        "$caller: ONELAB parameter must be a JSON object"))
    raw_name=get(obj,"name",nothing)
    raw_name isa AbstractString || throw(ArgumentError(
        "$caller: ONELAB parameter requires a \"name\" string"))
    name=String(raw_name)
    isempty(name) && throw(ArgumentError(
        "$caller: ONELAB parameter name must be nonempty"))
    raw_type=get(obj,"type","number")
    raw_type isa AbstractString || throw(ArgumentError(
        "$caller: ONELAB parameter \"type\" must be a string"))
    type=String(raw_type)
    type in ("number","string") || throw(ArgumentError(
        "$caller: unsupported ONELAB parameter type $(repr(type))"))
    p=_OnelabParameter(name,Symbol(type))
    for (key,value) in obj
        (key=="type" || key=="name") && continue
        if key=="values"
            if p.kind===:number
                p.values_f=_onelab_json_floats(value,caller,"\"values\"")
            else
                p.values_s=_onelab_json_strings(value,caller,"\"values\"")
            end
        elseif key=="min" && p.kind===:number
            p.min=_onelab_json_scalar(value,caller,"\"min\"")
        elseif key=="max" && p.kind===:number
            p.max=_onelab_json_scalar(value,caller,"\"max\"")
        elseif key=="step" && p.kind===:number
            p.step=_onelab_json_scalar(value,caller,"\"step\"")
        elseif key=="choices"
            if p.kind===:number
                p.choices_f=_onelab_json_floats(value,caller,"\"choices\"")
            else
                p.choices_s=_onelab_json_strings(value,caller,"\"choices\"")
            end
        elseif key=="valueLabels"
            p.value_labels=_onelab_json_string_map(value,caller,
                "\"valueLabels\"")
        elseif key=="attributes"
            p.attributes=_onelab_json_string_map(value,caller,"\"attributes\"")
        elseif key=="label"
            p.label=_onelab_json_optstring(value,caller,"\"label\"")
        elseif key=="shortHelp"
            p.short_help=_onelab_json_optstring(value,caller,"\"shortHelp\"")
        elseif key=="longHelp"
            p.long_help=_onelab_json_optstring(value,caller,"\"longHelp\"")
        elseif key=="readOnly"
            b=_onelab_json_bool(value)
            b===nothing && throw(ArgumentError(
                "$caller: \"readOnly\" must be a boolean"))
            p.read_only=b
        elseif key=="visible"
            b=_onelab_json_bool(value)
            b===nothing && throw(ArgumentError(
                "$caller: \"visible\" must be a boolean"))
            p.visible=b
        elseif key=="closed"
            b=_onelab_json_bool(value)
            b===nothing && throw(ArgumentError(
                "$caller: \"closed\" must be a boolean"))
            p.closed=b
        elseif key=="kind"
            k=_onelab_json_optstring(value,caller,"\"kind\"")
            k===nothing || (p.kind_name=k)
        end
    end
    return p
end

_onelab_number_record(name::String,values::Vector{Float64})=
    (p=_OnelabParameter(name,:number);p.values_f=values;p)
_onelab_string_record(name::String,values::Vector{String})=
    (p=_OnelabParameter(name,:string);p.values_s=values;p)

# ---- glob matching ---------------------------------------------------------------

# `*` matches any (possibly empty) character run, `?` exactly one character;
# an empty pattern matches everything (upstream `getNames("")` behaviour).
function _onelab_glob_match(name::String,pattern::String)
    isempty(pattern) && return true
    ni=firstindex(name);nn=lastindex(name)
    pi=firstindex(pattern);pn=lastindex(pattern)
    star_p=0;star_n=0
    while ni<=nn
        if pi<=pn && (pattern[pi]=='?' || pattern[pi]==name[ni])
            ni=nextind(name,ni);pi=nextind(pattern,pi)
        elseif pi<=pn && pattern[pi]=='*'
            star_p=pi;star_n=ni;pi=nextind(pattern,pi)
        elseif star_p!=0
            # Backtrack: `*` consumes one more name character and the pattern
            # resumes right after it (`nextind` at `nn` yields the end index,
            # letting `*` absorb the tail and exit the loop).
            pi=nextind(pattern,star_p)
            star_n=nextind(name,star_n)
            ni=star_n
        else
            return false
        end
    end
    while pi<=pn && pattern[pi]=='*';pi=nextind(pattern,pi);end
    return pi>pn
end

# ---- session endpoints --------------------------------------------------------------

function _onelab_set(data,format)
    caller="API.onelab.set"
    data isa AbstractString || throw(ArgumentError(
        "$caller: data must be a JSON string"))
    _onelab_check_format(format,caller)
    lock(STATE_LOCK) do
        _model_locked()
        parsed=_onelab_json_parse(String(data),caller)
        items=parsed isa Vector ? parsed : Any[parsed]
        for item in items
            p=_onelab_param_from_json(item,caller)
            SESSION_ONELAB[(p.kind,p.name)]=p
            _onelab_mark_changed!((p.kind,p.name))
        end
        return nothing
    end
end

function _onelab_get(name,format)
    caller="API.onelab.get"
    key=_onelab_check_name(name,caller)
    _onelab_check_format(format,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        preferred=isempty(key) ? nothing : _onelab_preferred_key(SESSION_ONELAB,key)
        parameters=isempty(key) ? sort!(collect(keys(SESSION_ONELAB));
            by=item->(item[2],item[1]===:number ? 0 : 1)) :
            preferred===nothing ? _OnelabParameterKey[] : [preferred]
        isempty(parameters) && return "[]"
        io=IOBuffer();write(io,"[\n")
        for (i,n) in enumerate(parameters)
            i>1 && write(io,",\n")
            _onelab_param_json!(io,SESSION_ONELAB[n],caller)
        end
        write(io,"\n]")
        return String(take!(io))
    end
end

function _onelab_get_names(search)
    caller="API.onelab.get_names"
    pattern=_onelab_check_name(search,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        return sort!(String[n for (_,n) in keys(SESSION_ONELAB)
                            if _onelab_glob_match(n,pattern)])
    end
end

function _onelab_set_number(name,values)
    caller="API.onelab.set_number"
    key=_onelab_check_name(name,caller)
    isempty(key) && throw(ArgumentError("$caller: name must be nonempty"))
    vals=_onelab_floats(values,caller)
    lock(STATE_LOCK) do
        _model_locked()
        SESSION_ONELAB[(:number,key)]=_onelab_number_record(key,vals)
        _onelab_mark_changed!((:number,key))
        return nothing
    end
end

function _onelab_set_string(name,values)
    caller="API.onelab.set_string"
    key=_onelab_check_name(name,caller)
    isempty(key) && throw(ArgumentError("$caller: name must be nonempty"))
    vals=_onelab_strings(values,caller)
    lock(STATE_LOCK) do
        _model_locked()
        SESSION_ONELAB[(:string,key)]=_onelab_string_record(key,vals)
        _onelab_mark_changed!((:string,key))
        return nothing
    end
end

function _onelab_get_number(name)
    caller="API.onelab.get_number"
    key=_onelab_check_name(name,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        p=get(SESSION_ONELAB,(:number,key),nothing)
        (p===nothing || p.kind!==:number) && throw(ArgumentError(
            "$caller: no ONELAB number parameter $(repr(key))"))
        return copy(p.values_f)
    end
end

function _onelab_get_string(name)
    caller="API.onelab.get_string"
    key=_onelab_check_name(name,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        p=get(SESSION_ONELAB,(:string,key),nothing)
        (p===nothing || p.kind!==:string) && throw(ArgumentError(
            "$caller: no ONELAB string parameter $(repr(key))"))
        return copy(p.values_s)
    end
end

function _onelab_clear(name)
    caller="API.onelab.clear"
    key=_onelab_check_name(name,caller)
    lock(STATE_LOCK) do
        _model_locked()
        if isempty(key)
            for n in collect(keys(SESSION_ONELAB))
                _onelab_store_delete!(_ONELAB_SESSION_STORE,n)
            end
            empty!(SESSION_ONELAB_CHANGED)
        else
            for n in unique!([name for (_,name) in keys(SESSION_ONELAB)])
                _onelab_glob_match(n,key) && _onelab_store_delete!(_ONELAB_SESSION_STORE,n)
            end
            filter!(n->!_onelab_glob_match(n,key),SESSION_ONELAB_CHANGED)
        end
        return nothing
    end
end

_onelab_run(name,command;kwargs...)=_onelab_run_session(name,command;kwargs...)

function _onelab_set_changed(name)
    caller="API.onelab.set_changed"
    key=_onelab_check_name(name,caller)
    lock(STATE_LOCK) do
        _model_locked()
        if haskey(_ONELAB_SESSION_STORE.client_changed,key)
            _onelab_store_set_client_changed!(_ONELAB_SESSION_STORE,key,true)
        else
            _onelab_mark_changed!(key)
        end
        return nothing
    end
end

function _onelab_get_changed(name)
    caller="API.onelab.get_changed"
    key=_onelab_check_name(name,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        return get(_ONELAB_SESSION_STORE.client_changed,key,key in SESSION_ONELAB_CHANGED)
    end
end

function _onelab_clear_changed()
    lock(STATE_LOCK) do
        _model_locked()
        empty!(SESSION_ONELAB_CHANGED)
        empty!(_ONELAB_SESSION_STORE.typed_changed)
        for name in keys(_ONELAB_SESSION_STORE.client_changed)
            _onelab_store_set_client_changed!(_ONELAB_SESSION_STORE,name,false)
        end
        return nothing
    end
end

# Sibling interface modules included here so they land inside `module API`
# like `module onelab` (API.jl's include list is owned elsewhere).
include("APIParser.jl")
include("APIFltk.jl")
include("APIGraphics.jl")
include("APIOnelabClient.jl")
include("APIOnelabServer.jl")
