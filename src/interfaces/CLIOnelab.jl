# ONELAB batch execution keeps its parameter server while the metamodel's
# geometry is read afresh for each action. Client commands run through the
# same owned process/socket engine as OnelabRun and API.onelab.run.
function _run_onelab_geo(path,preludes,numeric_options)
    runner=API._onelab_geo_runner()
    temporary=isempty(preludes) ? nothing : _inject_geo(path,preludes)
    source=temporary===nothing ? path : temporary
    solver_options=Dict{String,Float64}()
    recorder=(family,member,value)->begin
        family=="Solver" && (solver_options[member]=Float64(value))
        return nothing
    end
    database=splitext(path)[1]*".db"
    try
        for action in ("","check")
            if action=="check" && !iszero(get(solver_options,"AutoLoadDatabase",0.)) && isfile(database)
                API._onelab_database_load!(runner.store,database)
                API._onelab_geo_refresh!(runner)
            end
            result=Base.invokelatest(execute_geo,source;mesh_dim=0,
                _forced_numeric_options=numeric_options,
                _initial_numeric_option_callback=recorder,
                _initial_onelab_numbers=runner.numbers,
                _initial_onelab_strings=runner.strings,
                _initial_onelab_runner=runner,_initial_onelab_action=action)
            API._onelab_geo_sync!(runner)
            runner.remote===nothing || API._onelab_remote_model_info!(runner.remote;path,source=result.model)
            failure=_cli_geo_batch_error(result)
            failure===nothing || throw(failure)
            result.exit_code===nothing || return path
            API._onelab_geo_refresh!(runner)
        end
        iszero(get(solver_options,"AutoSaveDatabase",1.)) ||
            API._onelab_database_write(runner.store,database)
    finally
        temporary===nothing || rm(temporary;force=true)
    end
    return path
end

const _ONELAB_CLI_BATCH_FLAGS=("-","-parse_and_exit","-1","-2","-3","-0","-run",
    "-save","-refine","-check","-reclassify","-reparam","-barycentric_refine")

function _onelab_cli_value_count(flag,tokens,index)
    flag in ("-setnumber","-setstring") && return 2
    if flag in ("-setlist","-setlistofnumbers")
        index+2<=length(tokens) || return 2
        count=tryparse(Int,tokens[index+2])
        return 2+(count===nothing ? 0 : clamp(count,0,length(tokens)))
    elseif flag in ("-part_weight","-partWeight")
        count=0
        while index+count+1<=length(tokens) && lowercase(tokens[index+count+1]) in
                ("tri","quad","tet","hex","pri","pyr","trih","lin")
            count+=2
        end
        return count
    end
    return haskey(_NUMERIC_FLAGS,flag) || flag in ("-o","-format","-f","-vmsh","-algo",
        "-log","-statreport","-append_statreport","-reclassify","-reparam","-part",
        "-string","-option","-bgm") ? 1 : 0
end

function _onelab_cli_flags(tokens)
    ordinary=String[];selected=nothing;index=1;batch=false;metamodel=false
    while index<=length(tokens)
        flag=tokens[index]
        if flag in ("-onelab","-socket")
            selected===nothing || throw(ArgumentError("tessella: duplicate or conflicting ONELAB client flags"))
            index+=1
            index<=length(tokens) && !startswith(tokens[index],'-') && !isempty(tokens[index]) ||
                throw(ArgumentError("tessella: $flag needs $(flag=="-socket" ? "an address" : "a client name")"))
            name=flag=="-socket" ? "GmshRemote" : tokens[index]
            address=flag=="-socket" ? tokens[index] : nothing
            if flag=="-onelab" && index<length(tokens) && !startswith(tokens[index+1],'-')
                index+=1;address=tokens[index]
            end
            address===nothing || !isempty(address) || throw(ArgumentError("tessella: $flag needs an address"))
            selected=(;name,address,service=flag=="-socket")
        else
            batch|=flag in _ONELAB_CLI_BATCH_FLAGS
            metamodel|=flag=="-run"
            push!(ordinary,flag)
            if flag=="-convert"
                append!(ordinary,@view(tokens[index+1:end]))
                break
            end
            count=min(_onelab_cli_value_count(flag,tokens,index),length(tokens)-index)
            append!(ordinary,@view(tokens[index+1:index+count]))
            index+=count
        end
        index+=1
    end
    return (selected===nothing ? nothing : merge(selected,(;batch,metamodel))),ordinary
end

function _run_onelab_empty(preludes,numeric_options)
    statements=String[]
    for (name,value) in numeric_options
        push!(statements,"$name=$value;")
    end
    for prelude in preludes
        push!(statements,prelude isa String ? prelude : "$(prelude[1])=$(prelude[2]);")
    end
    isempty(statements) || API._onelab_live_parse_string!(join(statements,'\n'))
    return ""
end

function _run_onelab_remote_cli(tokens,selected)
    epoch=lock(API.STATE_LOCK) do
        API.initialize()
        API._ONELAB_SESSION_EPOCH[]
    end
    try
        return API._onelab_remote_scope(selected.name,selected.address,epoch) do session
            result=lock(API.STATE_LOCK) do
                API._onelab_remote_valid(session)
                _run(tokens;_reuse_session=true)
            end
            if !selected.metamodel || !session.model_info_sent
                lock(API.STATE_LOCK) do
                    API._onelab_remote_valid(session)
                    API._onelab_remote_model_info!(session)
                end
            end
            if selected.service && !selected.batch
                code=API._onelab_remote_service!(session)
                code==0 || throw(ArgumentError("tessella: remote GEO Exit $code"))
            end
            return result
        end
    finally
        lock(API.STATE_LOCK) do
            epoch==API._ONELAB_SESSION_EPOCH[] && API.finalize()
        end
    end
end

function _run_onelab_cli(files,preludes,numeric_options)
    if isempty(files)
        isempty(preludes) && isempty(numeric_options) && return ""
        return mktempdir() do directory
            empty_source=joinpath(directory,"onelab.geo");write(empty_source,"")
            _run_onelab_geo(empty_source,preludes,numeric_options)
            return ""
        end
    end
    length(files)==1 || throw(ArgumentError(
        "tessella: -run requires a single GEO metamodel"))
    path=first(only(files))
    lowercase(splitext(path)[2]) in (".geo","") || throw(ArgumentError(
        "tessella: -run requires a GEO metamodel"))
    return _run_onelab_geo(path,preludes,numeric_options)
end
