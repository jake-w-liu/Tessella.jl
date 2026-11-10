# Native ONELAB database records contain a decimal byte count, one space,
# the versioned wire record including its final NUL, and a line ending.
const _ONELAB_DATABASE_BYTES=64*1024*1024

function _onelab_database_read(path;max_bytes::Int=_ONELAB_DATABASE_BYTES)
    caller="ONELAB database"
    max_bytes>0 || throw(ArgumentError("$caller: byte capacity must be positive"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    size=filesize(path)
    size<=max_bytes || throw(ArgumentError("$caller: file exceeds byte capacity"))
    return Base.open(path,"r") do input
        header=readline(input)
        startswith(header,"ONELAB database") || throw(ArgumentError("$caller: invalid header"))
        records=_OnelabWireParameter[]
        while !eof(input)
            count=0;digits=0
            while true
                eof(input) && throw(ArgumentError("$caller: truncated record count"))
                byte=read(input,UInt8)
                byte==UInt8(' ') && break
                UInt8('0')<=byte<=UInt8('9') || throw(ArgumentError("$caller: invalid record count"))
                digit=Int(byte-UInt8('0'))
                count<=div(max_bytes-digit,10) || throw(ArgumentError("$caller: record exceeds byte capacity"))
                count=10count+digit;digits+=1
                digits<=20 || throw(ArgumentError("$caller: invalid record count"))
            end
            digits>0 && count>1 && count<=size-position(input) || throw(ArgumentError(
                "$caller: truncated or empty record"))
            bytes=read(input,count)
            last(bytes)==0 || throw(ArgumentError("$caller: record lacks its NUL terminator"))
            eof(input) && throw(ArgumentError("$caller: record lacks its line ending"))
            ending=read(input,UInt8)
            if ending==UInt8('\r')
                eof(input) && throw(ArgumentError("$caller: truncated line ending"))
                ending=read(input,UInt8)
            end
            ending==UInt8('\n') || throw(ArgumentError("$caller: invalid line ending"))
            pop!(bytes)
            parameter=try
                _onelab_wire_decode(String(bytes),caller)
            catch error
                error isa InterruptException && rethrow()
                error isa ArgumentError && rethrow()
                throw(ArgumentError("$caller: malformed parameter record"))
            end
            isempty(parameter.name) && throw(ArgumentError("$caller: empty parameter name"))
            push!(records,parameter)
        end
        return records
    end
end

function _onelab_database_load!(store,path)
    records=_onelab_database_read(path)
    stage=_OnelabServerStore(copy(store.parameters),copy(store.changed),copy(store.typed_changed),
        deepcopy(store.wire),copy(store.client_changed),copy(store.commands),ReentrantLock())
    for record in records
        _onelab_store_receive!(stage,record,"Gmsh")
    end
    for field in (:parameters,:changed,:typed_changed,:wire,:client_changed)
        destination=getfield(store,field);empty!(destination)
        destination isa AbstractSet ? union!(destination,getfield(stage,field)) :
            merge!(destination,getfield(stage,field))
    end
    return nothing
end

function _onelab_database_check_size(parameters,remaining)
    for parameter in parameters
        for text in (parameter.name,something(parameter.label,""),
                something(parameter.short_help,""))
            count=ncodeunits(text)
            count<=remaining || throw(ArgumentError("ONELAB database: parameter exceeds byte capacity"))
            remaining-=count
        end
        if parameter.kind===:string
            for text in Iterators.flatten(((parameter.kind_name,),parameter.values_s,parameter.choices_s))
                count=ncodeunits(text)
                count<=remaining || throw(ArgumentError("ONELAB database: parameter exceeds byte capacity"))
                remaining-=count
            end
        else
            count=length(parameter.values_f)+length(parameter.choices_f)
            count<=div(remaining,2) || throw(ArgumentError("ONELAB database: parameter exceeds byte capacity"))
            remaining-=2count
        end
        for (key,value) in parameter.attributes
            count=ncodeunits(key)+ncodeunits(value)
            count<=remaining || throw(ArgumentError("ONELAB database: parameter exceeds byte capacity"))
            remaining-=count
        end
        if parameter.kind===:number
            for (key,label) in parameter.value_labels
                value=tryparse(Float64,key)
                value!==nothing && isfinite(value) || continue
                count=ncodeunits(label)
                count<=remaining || throw(ArgumentError("ONELAB database: parameter exceeds byte capacity"))
                remaining-=count
            end
        end
    end
    return nothing
end

function _onelab_database_write(store,path;max_bytes::Int=_ONELAB_DATABASE_BYTES)
    caller="ONELAB database"
    max_bytes>0 || throw(ArgumentError("$caller: byte capacity must be positive"))
    target=abspath(path)
    isdir(target) && throw(ArgumentError("$caller: output must be a file"))
    _onelab_database_check_size(values(store.parameters),max_bytes)
    temporary,output=mktemp(dirname(target))
    try
        header="ONELAB database created by onelab server on "*
            Libc.strftime("%a %b %e %H:%M:%S %Y",time())*"\n"
        ncodeunits(header)<=max_bytes || throw(ArgumentError("$caller: header exceeds byte capacity"))
        write(output,header)
        for key in sort!(collect(keys(store.parameters));by=key->(key[2],key[1]===:number ? 0 : 1))
            remaining=max_bytes-position(output)
            parameter=store.parameters[key]
            wire=_onelab_store_wire(store,parameter,"")
            encoded=_onelab_wire_encode(wire;max_bytes=remaining)
            count=ncodeunits(encoded)+1
            prefix=string(count)*" "
            count+ncodeunits(prefix)+1<=remaining || throw(ArgumentError("$caller: records exceed byte capacity"))
            write(output,prefix,encoded,UInt8(0),UInt8('\n'))
        end
        close(output)
        Base.Filesystem.rename(temporary,target)
    finally
        isopen(output) && close(output)
        isfile(temporary) && rm(temporary)
    end
    return nothing
end
