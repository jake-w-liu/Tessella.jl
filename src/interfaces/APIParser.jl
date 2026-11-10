# `parser` module: the `.geo` parser variable store, matching the gmsh 4.15.2
# `gmshParser*` API surface. Number variables are `Vector{Float64}`; string
# variables are `Vector{String}`; name queries use the upstream `search`
# regular-expression argument. `parse` executes a file against the current
# model and existing parser namespace. Model, mesh and parameter changes use
# the same detached preparation and native publication as ONELAB ParseString.

const _PARSER_NUMBERS = Dict{String,Vector{Float64}}()
const _PARSER_STRINGS = Dict{String,Vector{String}}()

function _parser_reset_session!()
    empty!(_PARSER_NUMBERS);empty!(_PARSER_STRINGS)
    return nothing
end

function _parser_check_name(name,caller::AbstractString)
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    occursin('\0',name) && throw(ArgumentError(
        "$caller: name contains a NUL byte"))
    return String(name)
end

function _parser_set_number(name,value)
    caller="API.parser.set_number"
    key=_parser_check_name(name,caller)
    isempty(key) && throw(ArgumentError("$caller: name must not be empty"))
    raw=value isa Real && !(value isa Bool) ? [value] :
        value isa AbstractVector ? collect(value) : throw(ArgumentError(
            "$caller: value must be a number or a vector of numbers"))
    numbers=Vector{Float64}(undef,length(raw))
    for i in eachindex(raw)
        v=raw[i]
        (v isa Real && !(v isa Bool)) || throw(ArgumentError(
            "$caller: value $i must be a number"))
        f=try Float64(v) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: value $i must be representable as Float64"))
        end
        isfinite(f) || throw(ArgumentError(
            "$caller: value $i is not finite"))
        numbers[i]=f
    end
    lock(STATE_LOCK) do
        _model_locked()
        _PARSER_NUMBERS[key]=numbers
        delete!(_PARSER_STRINGS,key)
        return nothing
    end
end

function _parser_get_number(name)
    caller="API.parser.get_number"
    key=_parser_check_name(name,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        # Upstream returns an empty vector for a missing variable.
        return copy(get(_PARSER_NUMBERS,key,Float64[]))
    end
end

function _parser_set_string(name,value)
    caller="API.parser.set_string"
    key=_parser_check_name(name,caller)
    isempty(key) && throw(ArgumentError("$caller: name must not be empty"))
    raw=value isa AbstractString ? [value] :
        value isa AbstractVector ? collect(value) : throw(ArgumentError(
            "$caller: value must be a string or a vector of strings"))
    strings=Vector{String}(undef,length(raw))
    for i in eachindex(raw)
        s=raw[i]
        s isa AbstractString || throw(ArgumentError(
            "$caller: value $i must be a string"))
        occursin('\0',s) && throw(ArgumentError(
            "$caller: value $i contains a NUL byte"))
        strings[i]=String(s)
    end
    lock(STATE_LOCK) do
        _model_locked()
        _PARSER_STRINGS[key]=strings
        delete!(_PARSER_NUMBERS,key)
        return nothing
    end
end

function _parser_get_string(name)
    caller="API.parser.get_string"
    key=_parser_check_name(name,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        # Upstream returns an empty vector for a missing variable.
        return copy(get(_PARSER_STRINGS,key,String[]))
    end
end

function _parser_get_names(search)
    caller="API.parser.get_names"
    pattern=_parser_check_name(search,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        names=String[]
        if isempty(pattern)
            append!(names,keys(_PARSER_NUMBERS))
            append!(names,keys(_PARSER_STRINGS))
        else
            re=try Regex(pattern) catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError(
                    "$caller: invalid regular expression '$pattern'"))
            end
            for key in keys(_PARSER_NUMBERS)
                occursin(re,key) && push!(names,key)
            end
            for key in keys(_PARSER_STRINGS)
                occursin(re,key) && push!(names,key)
            end
        end
        sort!(names)
        return names
    end
end

function _parser_clear(name)
    caller="API.parser.clear"
    key=_parser_check_name(name,caller)
    lock(STATE_LOCK) do
        _model_locked()
        if isempty(key)
            empty!(_PARSER_NUMBERS)
            empty!(_PARSER_STRINGS)
        else
            delete!(_PARSER_NUMBERS,key)
            delete!(_PARSER_STRINGS,key)
        end
        return nothing
    end
end

# ParseFile consumes the current namespace and model. A successful operation
# publishes the actual geometry/mesh together with its final symbol table.
function _parser_parse(path)
    caller="API.parser.parse"
    source=_parser_check_name(path,caller)
    isempty(source) && throw(ArgumentError("$caller: file name must not be empty"))
    return lock(STATE_LOCK) do
        _model_locked()
        isfile(source) || throw(ArgumentError("$caller: missing file $source"))
        Base.invokelatest(_onelab_live_parse_string!,"";filename=abspath(source))
        return nothing
    end
end

module parser

using ..API: _parser_set_number,_parser_get_number,_parser_set_string,
             _parser_get_string,_parser_get_names,_parser_clear,_parser_parse

"Return the sorted variable names matching the regular expression `search` (all when empty)."
get_names(search="")=_parser_get_names(search)
"Set number variable `name` to `value` (scalar or vector); create or update it."
set_number(name,value)=_parser_set_number(name,value)
"Set string variable `name` to `value` (string or vector); create or update it."
set_string(name,value)=_parser_set_string(name,value)
"Return the `Vector{Float64}` value of number variable `name` (empty when missing)."
get_number(name)=_parser_get_number(name)
"Return the `Vector{String}` value of string variable `name` (empty when missing)."
get_string(name)=_parser_get_string(name)
"Clear all parser variables, or remove the single variable `name`."
clear(name="")=_parser_clear(name)
"Parse `fileName` against the current model and existing parser variables."
parse(fileName)=_parser_parse(fileName)

end # module parser
