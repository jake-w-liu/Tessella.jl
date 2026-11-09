# Native text-template preprocessing. The parameter server does not parse a
# solver's input language: ordinary text passes through and OL.get exchanges
# parameter values; OL.iftrue/ifntrue and OL.include select template text.
mutable struct _OnelabTemplateState
    store::_OnelabServerStore
    client::String
    used::Bool
    bytes::Int
    output_bytes::Int
    active::Set{String}
    diagnostics::Function
end

function _onelab_template_call(source,start,caller)
    matched=match(r"^OL\.([A-Za-z][A-Za-z0-9_]*)\s*\(",SubString(source,start))
    matched===nothing && throw(ArgumentError("$caller: malformed OL command"))
    opening=start+ncodeunits(matched.match)-1;depth=1;index=nextind(source,opening)
    first=index
    while index<=lastindex(source)
        character=source[index]
        character=='(' && (depth+=1)
        character==')' && (depth-=1)
        depth==0 && return String(matched.captures[1]),String(source[first:prevind(source,index)]),nextind(source,index)
        index=nextind(source,index)
    end
    throw(ArgumentError("$caller: unterminated OL command"))
end

function _onelab_template_value(state,raw)
    name=String(strip(raw))
    parameter=get(state.store.parameters,name,nothing)
    if parameter===nothing || parameter.kind===:number && isempty(parameter.values_f) ||
            parameter.kind===:string && isempty(parameter.values_s)
        state.diagnostics("OL.get: unknown or empty parameter $(repr(name))")
        return nothing
    end
    _onelab_store_register_client!(state.store,parameter,state.client)
    state.used=true
    return parameter.kind===:number ? _onelab_wire_f64(first(parameter.values_f)) : first(parameter.values_s)
end

function _onelab_template_substitute(state,source)
    output=IOBuffer();index=firstindex(source)
    while index<=lastindex(source)
        location=findnext("OL.",source,index)
        location===nothing && (write(output,SubString(source,index));break)
        start=first(location)
        write(output,SubString(source,index,prevind(source,start)))
        command,argument,after=_onelab_template_call(source,start,"ONELAB preprocessing")
        command=="get" || throw(ArgumentError("ONELAB preprocessing: unknown inline OL.$command"))
        value=_onelab_template_value(state,argument)
        value===nothing && return "??"
        write(output,value);index=after
    end
    return String(take!(output))
end

function _onelab_template_write!(state,output,text)
    count=ncodeunits(text)
    state.output_bytes<=64*1024*1024-count || throw(ArgumentError("ONELAB preprocessing: output exceeds64MiB"))
    state.output_bytes+=count;write(output,text)
    return nothing
end

function _onelab_template_file!(state,output,path,depth=0)
    depth<64 || throw(ArgumentError("ONELAB preprocessing: include nesting exceeds64"))
    path=abspath(path)
    path in state.active && throw(ArgumentError("ONELAB preprocessing: recursive include $path"))
    isfile(path) || throw(ArgumentError("ONELAB preprocessing: missing template $path"))
    size=filesize(path)
    size<=64*1024*1024-state.bytes || throw(ArgumentError("ONELAB preprocessing: templates exceed64MiB"))
    state.bytes+=size
    push!(state.active,path)
    branches=Tuple{Bool,Bool,Bool}[] # parent active, condition, else seen
    active=true
    try
        Base.open(path,"r") do input
            while !eof(input)
                line=readline(input)
                _onelab_valid_text(line) || throw(ArgumentError("ONELAB preprocessing: template is not valid UTF-8"))
                text=strip(line)
                if startswith(text,"OL.iftrue(") || startswith(text,"OL.ifntrue(")
                    length(branches)<64 || throw(ArgumentError("ONELAB preprocessing: condition nesting exceeds64"))
                    command,argument,after=_onelab_template_call(String(text),firstindex(text),"ONELAB preprocessing")
                    isempty(strip(SubString(text,after))) || throw(ArgumentError("ONELAB preprocessing: text after condition"))
                    value=active ? _onelab_template_value(state,argument) : nothing
                    number=value===nothing ? nothing : tryparse(Float64,value)
                    condition=number===nothing ? value!==nothing && !isempty(value) : !iszero(number)
                    command=="ifntrue" && (condition=!condition)
                    push!(branches,(active,condition,false));active=active && condition
                elseif text in ("OL.else","OL.else()")
                    isempty(branches) && throw(ArgumentError("ONELAB preprocessing: else without condition"))
                    parent,condition,seen=branches[end]
                    seen && throw(ArgumentError("ONELAB preprocessing: duplicate else"))
                    branches[end]=(parent,condition,true);active=parent && !condition
                elseif text in ("OL.endif","OL.endif()")
                    isempty(branches) && throw(ArgumentError("ONELAB preprocessing: endif without condition"))
                    active=pop!(branches)[1]
                elseif !active
                    continue
                elseif startswith(text,"OL.include(")
                    command,argument,after=_onelab_template_call(String(text),firstindex(text),"ONELAB preprocessing")
                    isempty(strip(SubString(text,after))) || throw(ArgumentError("ONELAB preprocessing: text after include"))
                    included=_onelab_template_substitute(state,String(strip(argument)))
                    _onelab_template_file!(state,output,joinpath(dirname(path),included),depth+1)
                else
                    _onelab_template_write!(state,output,_onelab_template_substitute(state,line)*"\n")
                end
            end
        end
        isempty(branches) || throw(ArgumentError("ONELAB preprocessing: unterminated condition"))
        # The reference preprocessor emits the final empty getline record.
        _onelab_template_write!(state,output,"\n")
    finally
        delete!(state.active,path)
    end
    return nothing
end

function _onelab_preprocess!(store,payload,epoch=nothing;directory=pwd(),diagnostics=message->nothing)
    fields=split(payload,'\0';keepempty=true)
    length(fields)==2 || throw(ArgumentError("ONELAB server: malformed preprocessing request"))
    client=String(fields[1]);filename=String(fields[2])
    !isempty(client) && !isempty(filename) || throw(ArgumentError("ONELAB preprocessing: name and filename must be nonempty"))
    path=abspath(joinpath(directory,filename))
    endswith(path,".ol") || throw(ArgumentError("ONELAB preprocessing: template filename must end in .ol"))
    destination=chop(path;tail=3)
    (isempty(basename(destination)) || isdir(destination)) && throw(ArgumentError(
        "ONELAB preprocessing: output must be a file, not a directory"))
    return _onelab_store_locked(store,epoch) do
        state=_OnelabTemplateState(store,client,false,0,0,Set{String}(),diagnostics)
        temporary,output=mktemp(dirname(destination))
        try
            _onelab_template_file!(state,output,path)
            close(output)
            # A file rename replaces a file atomically and cannot recursively
            # remove an existing output directory (unlike mv(force=true)).
            Base.Filesystem.rename(temporary,destination)
            return state.used
        finally
            isopen(output) && close(output)
            isfile(temporary) && rm(temporary)
        end
    end
end
