module APINativeGeoCPUTimeTestV2
using Tessella,Test
if Sys.iswindows()
    folder=normpath(joinpath(@__DIR__,"..","..","validation","native_geometry"))
    geo_folder=joinpath(folder,"geo_cpu_time_v1")
    include(joinpath(folder,"geo_cpu_time_v1_cpu_time_public_controls_v2.jl"))
    include(joinpath(folder,"geo_cpu_time_v1_defineconstant_basic_public_controls_v1.jl"))
    include(joinpath(folder,"geo_cpu_time_v1_cpu_actual_loaded_helper_controls_v1.jl"))
    # Only input-file locations move; the exact captured GEO bytes and all
    # snapshots, USER brackets, detached lifetime and allocation caps remain.
    for row in CPUTIMEPublicV2.CPU_POINT_ROWS_V1
        row["path"]=joinpath(geo_folder,row["kind"]*".geo")
    end
    for row in DefineConstantBasicPublicV1.DEFINE_CONSTANT_BASIC_ROWS_V1
        row["path"]=joinpath(geo_folder,row["form"]*".geo")
    end
    prefix=get(ENV,"TESSELLA_NATIVE_GEO_CPU_TIME_OUTPUT_PREFIX",joinpath(mktempdir(),"native_geo_cpu_time"))
    failures=String[]
    for (label,body) in (("ClockPublic6",CPUTIMEPublicV2.run),
                         ("DefineConstantBasic18",DefineConstantBasicPublicV1.run),
                         ("CPUActualLoadedHelper32",CPUActualLoadedHelperV1.run))
        println("CPU_PORTABLE_SCOPE_START ",label);flush(stdout)
        try
            body(prefix)
        catch error
            error isa InterruptException && rethrow()
            push!(failures,label)
            showerror(stdout,error,catch_backtrace());println()
        end
        println("CPU_PORTABLE_SCOPE_END ",label);flush(stdout)
    end
    isempty(failures) || error("CPU permanent wrapper failed scopes: "*repr(failures))
else
    @test_skip "Pinned native/public clock controls require Windows GetProcessTimes"
end
end
