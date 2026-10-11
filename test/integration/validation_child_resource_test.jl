using Test

@testset "Validation descendants keep explicit heap and CLI thread policy" begin
    repository=normpath(joinpath(@__DIR__,"..",".."))
    sources=("validation/run_all.jl","validation/geo_constraints/differential.jl",
        "validation/gmsh_parity/msh_size_t_width.jl")
    commands=Tuple{String,String}[]
    for relative in sources
        for line in eachline(joinpath(repository,relative))
            found=match(r"`\$\(Base\.julia_cmd\(\)\)[^`]*`",line)
            found===nothing && continue
            push!(commands,(relative,found.match))
        end
    end
    @test !isempty(commands)
    @test all(count("--heap-size-hint=2G",command)==1 for (_,command) in commands)
    @test count(occursin("--threads=4",command) for (_,command) in commands)==1
    before=(Base.JLOptions().heap_size_hint,Base.JLOptions().opt_level,
        Base.JLOptions().check_bounds,Threads.nthreads(),Threads.ngcthreads())
    mktempdir() do directory
        child=joinpath(directory,"child.jl")
        write(child,raw"""o=Base.JLOptions(); println("RESOURCE_CHILD heap=",Int(o.heap_size_hint)," threads=",Threads.nthreads()," gc=",Threads.ngcthreads()," opt=",o.opt_level," bounds=",o.check_bounds)""")
        representatives=Tuple{String,String}[]
        push!(representatives,first(filter(row->row[1]==sources[1] && !occursin("--threads=4",row[2]),commands)))
        push!(representatives,only(filter(row->occursin("--threads=4",row[2]),commands)))
        for relative in sources[2:3]
            push!(representatives,first(filter(row->row[1]==relative,commands)))
        end
        for (relative,source_command) in representatives
            # Transport the actual command AST, replacing fixture arguments
            # with a harmless stdlib child. Explicit O2/bounds makes this test
            # independent of the user's parent compiler/resource settings.
            scope=Module(gensym(:ValidationChildResource),true,true)
            for token in eachmatch(r"\$([A-Za-z_]\w*)",source_command)
                variable=token.captures[1]
                value=variable=="project" || endswith(variable,"_project") ? repository :
                    variable=="paths" ? [child] : child
                Core.eval(scope,Expr(:(=),Symbol(variable),QuoteNode(value)))
            end
            transport=replace(source_command,raw"$(Base.julia_cmd())"=>
                raw"$(Base.julia_cmd()) -O2 --check-bounds=yes";count=1)
            command=Core.eval(scope,Meta.parse(transport))
            @test command isa Cmd
            @test count(==("--heap-size-hint=2G"),command.exec)==1
            expected_threads=occursin("--threads=4",source_command) ? 4 : 1
            configured=addenv(command,"JULIA_NUM_THREADS"=>"1","JULIA_NUM_GC_THREADS"=>"1")
            actual=strip(read(configured,String))
            @test actual=="RESOURCE_CHILD heap=2147483648 threads=$expected_threads gc=1 opt=2 bounds=1"
        end
    end
    after=(Base.JLOptions().heap_size_hint,Base.JLOptions().opt_level,
        Base.JLOptions().check_bounds,Threads.nthreads(),Threads.ngcthreads())
    @test after==before
end
