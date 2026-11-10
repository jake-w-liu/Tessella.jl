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
