# Native Msg buffering has process lifetime: initialize/finalize do not stop
# an active logger. Terminal selection only controls the printed Msg channel;
# the experimental producer's direct progress output is a separate channel.
import ..Model
using Printf: @sprintf
const NATIVE_LOGGER_ACTIVE=Ref(false)
const NATIVE_LOGGER_MESSAGES=String[]

function _native_diagnostic_emit_locked(kind::Symbol,message::AbstractString;options=OPTIONS)
    if kind===:trace
        println(stdout,message);flush(stdout)
        return nothing
    end
    threshold=kind===:error ? 1 : kind===:warning ? 2 : kind===:info ? 4 : 99
    get(options,"General.Verbosity",5.0)>=threshold || return nothing
    label=kind===:error ? "Error" : kind===:warning ? "Warning" : kind===:info ? "Info" : "Debug"
    NATIVE_LOGGER_ACTIVE[] && push!(NATIVE_LOGGER_MESSAGES,string(label,": ",message))
    if !iszero(get(options,"General.Terminal",1.0))
        stream=kind in (:error,:warning) ? stderr : stdout
        prefix=kind===:error ? "Error   : " : kind===:warning ? "Warning : " : kind===:info ? "Info    : " : "Debug   : "
        println(stream,prefix,message);flush(stream)
    end
    return nothing
end

_native_detached_diagnostic_emit(kind,message,options)=
    _native_diagnostic_emit_locked(kind,message;options)
_native_detached_diagnostic_errors()=NATIVE_LOOKUP_QUERY_ERRORS[]
function _native_detached_diagnostic_error(message,options)
    NATIVE_LOOKUP_QUERY_ERRORS[]+=1
    NATIVE_MSG_LAST_ERROR[]=String(message)
    _native_diagnostic_emit_locked(:error,message;options)
    abort=get(options,"General.AbortOnError",2.0)
    (abort==2 || abort==3) && throw(ErrorException(message))
    abort==4 && exit(1)
    return nothing
end

function _native_logger_start()
    return lock(STATE_LOCK) do
        _model_locked()
        if NATIVE_LOGGER_ACTIVE[]
            _native_diagnostic_emit_locked(:warning,"Logger already started - ignoring")
        else
            empty!(NATIVE_LOGGER_MESSAGES);NATIVE_LOGGER_ACTIVE[]=true
        end
        nothing
    end
end
function _native_logger_stop()
    return lock(STATE_LOCK) do
        _model_locked()
        if NATIVE_LOGGER_ACTIVE[]
            NATIVE_LOGGER_ACTIVE[]=false;empty!(NATIVE_LOGGER_MESSAGES)
        else
            _native_diagnostic_emit_locked(:warning,"Logger not started - ignoring")
        end
        nothing
    end
end
function _native_logger_get()
    return lock(STATE_LOCK) do
        _model_locked();copy(NATIVE_LOGGER_MESSAGES)
    end
end

function _native_diagnostic_cpu_time()
    if Sys.iswindows()
        creation=Ref{UInt64}(0);exit_time=Ref{UInt64}(0)
        kernel=Ref{UInt64}(0);user_time=Ref{UInt64}(0)
        handle=ccall((:GetCurrentProcess,"kernel32"),Ptr{Cvoid},())
        ok=ccall((:GetProcessTimes,"kernel32"),Int32,
            (Ptr{Cvoid},Ref{UInt64},Ref{UInt64},Ref{UInt64},Ref{UInt64}),
            handle,creation,exit_time,kernel,user_time)
        ok!=0 || throw(ErrorException("Native diagnostic CPU clock failed"))
        return Float64(kernel[]+user_time[])/1.e7
    end
    return Float64(ccall(:clock,Clong,()))/1.e6
end
function _native_diagnostic_recombine_done_locked(wall,cpu;empty=false)
    elapsed=empty ? 0.0 : max(0.0,Float64(time_ns()-wall)/1.e9)
    consumed=empty ? 0.0 : max(0.0,_native_diagnostic_cpu_time()-cpu)
    _native_diagnostic_emit_locked(:info,@sprintf(
        "Done recombining 2D mesh (Wall %gs, CPU %gs)",elapsed,consumed))
end

function _native_detached_diagnostic_done(wall,cpu,options;empty=false)
    elapsed=empty ? 0.0 : max(0.0,Float64(time_ns()-wall)/1.e9)
    consumed=empty ? 0.0 : max(0.0,_native_diagnostic_cpu_time()-cpu)
    _native_diagnostic_emit_locked(:info,@sprintf(
        "Done recombining 2D mesh (Wall %gs, CPU %gs)",elapsed,consumed);options)
end
Model._MODEL_NATIVE_RECOMBINE_DIAGNOSTICS[]=(
    emit=_native_detached_diagnostic_emit,error=_native_detached_diagnostic_error,
    errors=_native_detached_diagnostic_errors,cpu=_native_diagnostic_cpu_time,
    done=_native_detached_diagnostic_done)
