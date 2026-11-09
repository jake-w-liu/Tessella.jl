# Owned solver processes. Windows assigns the suspended process to a job before
# resuming it, so even a client that exits before its descendants can be reaped.
mutable struct _OnelabOwnedProcess
    process::Any
    handle::Ptr{Cvoid}
    job::Ptr{Cvoid}
    pid::Int
    exit_code::Union{Nothing,Int}
    closed::Bool
    lock::ReentrantLock
end

struct _OnelabWinSecurity
    length::UInt32
    descriptor::Ptr{Cvoid}
    inherit::Int32
end
struct _OnelabWinStartup
    size::UInt32
    reserved::Ptr{UInt16}
    desktop::Ptr{UInt16}
    title::Ptr{UInt16}
    x::UInt32
    y::UInt32
    width::UInt32
    height::UInt32
    xchars::UInt32
    ychars::UInt32
    fill::UInt32
    flags::UInt32
    show::UInt16
    reserved_size::UInt16
    reserved_data::Ptr{UInt8}
    input::Ptr{Cvoid}
    output::Ptr{Cvoid}
    error::Ptr{Cvoid}
end
struct _OnelabWinProcessInfo
    process::Ptr{Cvoid}
    thread::Ptr{Cvoid}
    pid::UInt32
    tid::UInt32
end
struct _OnelabWinBasicLimits
    process_time::Int64
    job_time::Int64
    flags::UInt32
    minimum_working_set::UInt
    maximum_working_set::UInt
    active_processes::UInt32
    affinity::UInt
    priority::UInt32
    scheduling::UInt32
end
struct _OnelabWinIOCounters
    read_operations::UInt64
    write_operations::UInt64
    other_operations::UInt64
    read_bytes::UInt64
    write_bytes::UInt64
    other_bytes::UInt64
end
struct _OnelabWinJobLimits
    basic::_OnelabWinBasicLimits
    io::_OnelabWinIOCounters
    process_memory::UInt
    job_memory::UInt
    peak_process_memory::UInt
    peak_job_memory::UInt
end

_onelab_win_error()=ccall((:GetLastError,"kernel32"),UInt32,())
_onelab_win_close(handle)=handle==C_NULL ? nothing :
    (ccall((:CloseHandle,"kernel32"),Int32,(Ptr{Cvoid},),handle);nothing)

function _onelab_command_text(value,caller,what)
    value isa AbstractString || throw(ArgumentError("$caller: $what must be a string"))
    result=String(value)
    isvalid(result) && !occursin('\0',result) || throw(ArgumentError(
        "$caller: $what must be valid text without NUL bytes"))
    return result
end

function _onelab_windows_quote(argument::String)
    io=IOBuffer();write(io,'"');slashes=0
    for character in argument
        if character=='\\'
            slashes+=1
        elseif character=='"'
            write(io,repeat("\\",2slashes+1));write(io,'"');slashes=0
        else
            write(io,repeat("\\",slashes));write(io,character);slashes=0
        end
    end
    write(io,repeat("\\",2slashes));write(io,'"')
    return String(take!(io))
end

function _onelab_windows_words(command::String)
    count=Ref{Int32}(0)
    array=ccall((:CommandLineToArgvW,"shell32"),Ptr{Ptr{UInt16}},
        (Cwstring,Ref{Int32}),command,count)
    array==C_NULL && throw(ArgumentError(
        "API.onelab.run: cannot parse command line (Windows error $(_onelab_win_error()))"))
    try
        maximum_units=length(transcode(UInt16,command))+1
        0<count[]<=maximum_units || throw(ArgumentError("API.onelab.run: invalid command line"))
        result=String[]
        for index in 1:count[]
            word=unsafe_load(array,index);length_units=0
            while length_units<maximum_units && unsafe_load(word,length_units+1)!=0
                length_units+=1
            end
            length_units<maximum_units || throw(ArgumentError("API.onelab.run: invalid command word"))
            push!(result,transcode(String,unsafe_wrap(Vector{UInt16},word,length_units;own=false)))
        end
        return result
    finally
        ccall((:LocalFree,"kernel32"),Ptr{Cvoid},(Ptr{Cvoid},),array)
    end
end

function _onelab_process_start(command::String,arguments::Vector{String};
        directory=pwd(),quiet=false)
    caller="API.onelab.run"
    command=String(strip(_onelab_command_text(command,caller,"command")))
    isempty(command) && throw(ArgumentError("$caller: command must be nonempty"))
    directory=abspath(_onelab_command_text(directory,caller,"working directory"))
    isdir(directory) || throw(ArgumentError("$caller: missing working directory $directory"))
    foreach(argument->_onelab_command_text(argument,caller,"argument"),arguments)
    if !Sys.iswindows()
        words=Base.shell_split(command)
        isempty(words) && throw(ArgumentError("$caller: command must name an executable"))
        # The Vector{String} constructor does not accept keywords; wrap the
        # resulting command in the keyword-bearing Cmd constructor.
        program=Cmd(Cmd(vcat(words,arguments));dir=directory,detach=true)
        process=try
            run(pipeline(program;stdin=devnull,
                stdout=quiet ? devnull : stdout,
                stderr=quiet ? devnull : stderr);wait=false)
        catch err
            err isa InterruptException && rethrow()
            err isa Base.IOError && throw(ArgumentError(
                "$caller: cannot launch client: $(err.msg)"))
            rethrow()
        end
        return _OnelabOwnedProcess(process,C_NULL,C_NULL,Int(getpid(process)),nothing,false,ReentrantLock())
    end

    words=_onelab_windows_words(command)
    isempty(words) && throw(ArgumentError("$caller: command must name an executable"))
    text=command*join((" "*_onelab_windows_quote(argument) for argument in arguments))
    line=push!(transcode(UInt16,text),0)
    length(line)<=32767 || throw(ArgumentError("$caller: Windows command line exceeds32767 UTF-16 units"))
    location=push!(transcode(UInt16,directory),0)
    # CreateProcess resolves relative executable paths against the parent's
    # directory. Supply the actual file explicitly for GEO-relative programs.
    first_word=words[1]
    relative_program=!isabspath(first_word) && (occursin('/',first_word) || occursin('\\',first_word))
    application=relative_program ? push!(transcode(UInt16,joinpath(directory,first_word)),0) : UInt16[]
    job=ccall((:CreateJobObjectW,"kernel32"),Ptr{Cvoid},(Ptr{Cvoid},Ptr{UInt16}),C_NULL,C_NULL)
    job==C_NULL && throw(ArgumentError("$caller: cannot create owned process job (Windows error $(_onelab_win_error()))"))
    info=Ref(_OnelabWinProcessInfo(C_NULL,C_NULL,0,0));null_handle=Ptr{Cvoid}(0)
    created=false
    try
        limits=Ref(_OnelabWinJobLimits(
            _OnelabWinBasicLimits(0,0,0x00002000,0,0,0,0,0,0),
            _OnelabWinIOCounters(0,0,0,0,0,0),0,0,0,0))
        ccall((:SetInformationJobObject,"kernel32"),Int32,
            (Ptr{Cvoid},Int32,Ref{_OnelabWinJobLimits},UInt32),
            job,9,limits,sizeof(_OnelabWinJobLimits))!=0 || throw(ArgumentError(
                "$caller: cannot configure owned process job (Windows error $(_onelab_win_error()))"))
        security=Ref(_OnelabWinSecurity(sizeof(_OnelabWinSecurity),C_NULL,1))
        null_handle=ccall((:CreateFileW,"kernel32"),Ptr{Cvoid},
            (Cwstring,UInt32,UInt32,Ref{_OnelabWinSecurity},UInt32,UInt32,Ptr{Cvoid}),
            "NUL",0xc0000000,3,security,3,0x80,C_NULL)
        null_handle==Ptr{Cvoid}(typemax(UInt)) && throw(ArgumentError(
            "$caller: cannot open solver standard streams (Windows error $(_onelab_win_error()))"))
        out=quiet ? null_handle : ccall((:GetStdHandle,"kernel32"),Ptr{Cvoid},(UInt32,),reinterpret(UInt32,Int32(-11)))
        err=quiet ? null_handle : ccall((:GetStdHandle,"kernel32"),Ptr{Cvoid},(UInt32,),reinterpret(UInt32,Int32(-12)))
        startup=Ref(_OnelabWinStartup(sizeof(_OnelabWinStartup),C_NULL,C_NULL,C_NULL,
            0,0,0,0,0,0,0,0x00000101,0,0,C_NULL,null_handle,out,err))
        success=GC.@preserve line location application begin
            ccall((:CreateProcessW,"kernel32"),Int32,
                (Ptr{UInt16},Ptr{UInt16},Ptr{Cvoid},Ptr{Cvoid},Int32,UInt32,
                 Ptr{Cvoid},Ptr{UInt16},Ref{_OnelabWinStartup},Ref{_OnelabWinProcessInfo}),
                isempty(application) ? C_NULL : pointer(application),pointer(line),C_NULL,C_NULL,
                1,0x08000204,C_NULL,pointer(location),startup,info)
        end
        success!=0 || throw(ArgumentError("$caller: cannot launch client (Windows error $(_onelab_win_error()))"))
        created=true
        ccall((:AssignProcessToJobObject,"kernel32"),Int32,(Ptr{Cvoid},Ptr{Cvoid}),
            job,info[].process)!=0 || throw(ArgumentError(
                "$caller: cannot own client process tree (Windows error $(_onelab_win_error()))"))
        ccall((:ResumeThread,"kernel32"),UInt32,(Ptr{Cvoid},),info[].thread)!=typemax(UInt32) ||
            throw(ArgumentError("$caller: cannot start client (Windows error $(_onelab_win_error()))"))
        _onelab_win_close(info[].thread)
        return _OnelabOwnedProcess(nothing,info[].process,job,Int(info[].pid),nothing,false,ReentrantLock())
    catch
        created && ccall((:TerminateProcess,"kernel32"),Int32,(Ptr{Cvoid},UInt32),info[].process,1)
        created && _onelab_win_close(info[].thread)
        created && _onelab_win_close(info[].process)
        _onelab_win_close(job)
        rethrow()
    finally
        null_handle!=C_NULL && null_handle!=Ptr{Cvoid}(typemax(UInt)) && _onelab_win_close(null_handle)
    end
end

_onelab_exit_status(code,signal)=iszero(signal) ? Int(code) : 128+Int(signal)

function _onelab_process_exited(process::_OnelabOwnedProcess)
    return lock(process.lock) do
    process.exit_code!==nothing && return true
    if Sys.iswindows()
        ccall((:WaitForSingleObject,"kernel32"),UInt32,(Ptr{Cvoid},UInt32),process.handle,0)==0 || return false
        code=Ref{UInt32}(0)
        ccall((:GetExitCodeProcess,"kernel32"),Int32,(Ptr{Cvoid},Ref{UInt32}),process.handle,code)!=0 ||
            throw(ArgumentError("API.onelab.run: cannot read client exit status"))
        process.exit_code=Int(code[])
    else
        Base.process_exited(process.process) || return false
        process.exit_code=_onelab_exit_status(process.process.exitcode,process.process.termsignal)
    end
    return true
    end
end

function _onelab_process_close!(process::_OnelabOwnedProcess)
    return lock(process.lock) do
    process.closed && return nothing
    process.closed=true
    if Sys.iswindows()
        ccall((:TerminateJobObject,"kernel32"),Int32,(Ptr{Cvoid},UInt32),process.job,1)
        ccall((:WaitForSingleObject,"kernel32"),UInt32,(Ptr{Cvoid},UInt32),process.handle,5000)
        _onelab_win_close(process.handle);_onelab_win_close(process.job)
    else
        # detach=true gives this solver a fresh group, including its workers.
        ccall(:kill,Int32,(Int32,Int32),-process.pid,Base.SIGTERM)
        ccall(:kill,Int32,(Int32,Int32),-process.pid,Base.SIGKILL)
        try wait(process.process) catch end
    end
    process.exit_code===nothing && (process.exit_code=-1)
    return nothing
    end
end
