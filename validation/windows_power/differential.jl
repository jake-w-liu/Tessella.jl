# Run with --project=<repository> --startup-file=no --check-bounds=yes --threads=4.
# Optional arguments: --python PATH --gmsh-executable PATH. Environment aliases:
# GMSH_PYTHON_EXECUTABLE and GMSH_EXECUTABLE. Python must provide gmsh==4.15.2.
# No reference backend is loaded by production or the ordinary unit tests.
using Tessella
using Test

Sys.iswindows() && Sys.ARCH===:x86_64 || error(
    "Windows power differential requires Windows x86_64; other platforms use libm")

function validation_arguments(args)
    python=get(ENV,"GMSH_PYTHON_EXECUTABLE",something(Sys.which("python"),"python"))
    # Python wheel gmsh.bat is a hosted-DLL launcher and inherits PC53; this
    # oracle intentionally checks native MinGW executable startup at PC64.
    executable=get(ENV,"GMSH_EXECUTABLE",something(Sys.which("gmsh.exe"),""))
    index=1
    while index<=length(args)
        flag=args[index]
        flag in ("--python","--gmsh-executable") || error("unknown argument: $flag")
        index<length(args) || error("missing value for $flag")
        if flag=="--python"
            python=args[index+1]
        else
            executable=args[index+1]
        end
        index+=2
    end
    return python,executable
end

python,executable=validation_arguments(ARGS)
primary=joinpath(@__DIR__,"primary.py")
command=`$python $primary`
if isempty(executable)
    error("Standalone primary gate needs --gmsh-executable PATH or GMSH_EXECUTABLE")
end
run(`$command --gmsh-executable $executable`)

# These tests load production definitions, committed primary bits, and owned
# temporary CLI fixtures; they do not redefine kernels or use audit scratch files.
include(joinpath(@__DIR__,"../../test/core/gmsh_windows_power_test.jl"))
include(joinpath(@__DIR__,"../../test/interfaces/cli_precision_test.jl"))

function resource_report(inputs)
    libm=Tessella.GmshLibm
    saved=libm._win_x87_control_word()
    try
        libm._win_x87_set_control_word!((saved&0xf0ff)|UInt16(0x0200))
        for count in (1000,2000,4000)
            output=Vector{Float64}(undef,count)
            WindowsPowerTest.batch!(output,inputs)
            allocated=Int[];seconds=Float64[]
            for _ in 1:7
                result=@timed WindowsPowerTest.batch!(output,inputs)
                @assert result.value===output
                push!(allocated,result.bytes)
                push!(seconds,result.time)
            end
            sort!(allocated);sort!(seconds)
            @assert allocated[4]==0
            println("power_batch count=",count," median7_alloc_bytes=",allocated[4],
                    " retained_output_bytes=",Base.summarysize(output),
                    " median7_seconds=",seconds[4])
        end
    finally
        libm._win_x87_set_control_word!(saved)
    end
end

rows=[split(row,'\t') for row in readlines(joinpath(@__DIR__,"corpus.tsv"))[2:end]]
inputs=[(reinterpret(Float64,parse(UInt64,row[2];base=16)),
         reinterpret(Float64,parse(UInt64,row[3];base=16))) for row in rows]
resource_report(inputs)
println("windows_power_differential PASS primary_cases=2874 modes=12 cli_cases=132")
